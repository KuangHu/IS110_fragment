#!/usr/bin/env python3
"""Prepare annotated GenBank files for visual review of strict rearrangement calls.

For each selected case (inversion, translocation, duplication) emits paired
ref.gbk + tgt.gbk files in
    review_cases/{type}/case_{N}_{ref_short}_vs_{tgt}/
with feature annotations:
  - IS110 transposase positions (CDS features from is_hits.tsv)
  - rearrangement boundaries (misc_feature with /note describing the call)
"""
import argparse, csv, gzip, json, os, re, shutil
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--selection", required=True)
    p.add_argument("--out-dir",   required=True)
    p.add_argument("--runs-dir",  required=True,
                   help="cross_ref_is_runs/results — for is_hits + records per species")
    p.add_argument("--src-dirs",  required=True,
                   help="JSON: species → source FASTA dir")
    p.add_argument("--focal-window", type=int, default=50000,
                   help="bp of context around anchors / inversion segment")
    return p.parse_args()


def load_is_hits(path):
    """assembly → contig (local) → [(s, e, strand)]"""
    idx = defaultdict(lambda: defaultdict(list))
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            asm = r["assembly"]
            contig = r["contig"]
            local = contig.split("|", 1)[1] if "|" in contig else contig
            idx[asm][local].append((int(r["tnp_start"]), int(r["tnp_end"]),
                                    r.get("tnp_strand", "+")))
    return idx


def load_locus(records_path, ref_id):
    with open(records_path) as f:
        recs = json.load(f)
    for r in recs:
        rid = r.get("ref_id") or r.get("is110_id")
        if rid != ref_id: continue
        src = r.get("source", {})
        ie = r.get("is_element") or src.get("is_element") or {}
        contig = src.get("contig")
        local = contig.split("|", 1)[1] if "|" in contig else contig
        return (local,
                int(ie.get("source_start") or src.get("transposase_start")),
                int(ie.get("source_end")   or src.get("transposase_end")))
    return None


def load_assembly_seqs(src_dir, asm):
    for ext in (".fna.gz", ".fa.gz", ".fna", ".fa"):
        path = os.path.join(src_dir, asm + ext)
        if os.path.exists(path):
            opener = gzip.open if path.endswith(".gz") else open
            mode = "rt"
            recs = {}
            with opener(path, mode) as fh:
                for rec in SeqIO.parse(fh, "fasta"):
                    local = rec.id.split("|", 1)[1] if "|" in rec.id else rec.id
                    recs[local] = rec
                    recs[rec.id] = rec  # also key by full name
            return recs
    return {}


ANCHOR_RE = re.compile(r"([\w.|]+):(\d+)-(\d+)\(([+-])\)")
def parse_anchors(details):
    # the contig field in details is already in the local-after-assembly form
    # (e.g., "ENA|CM132337|CM132337.1") — this is the same form is_hits uses,
    # so do NOT strip further pipes.
    return [(m.group(1), int(m.group(2)), int(m.group(3)), m.group(4))
            for m in ANCHOR_RE.finditer(details)]


def make_region_record(seq_record, sub_start, sub_end, name_suffix):
    """Slice a SeqRecord to [sub_start, sub_end] (0-based)."""
    sub = seq_record[sub_start:sub_end]
    sub.id = f"{seq_record.id[:15]}_{name_suffix}"[:16]
    sub.name = sub.id
    sub.description = f"{seq_record.description} | region {sub_start+1}-{sub_end}"
    sub.annotations["molecule_type"] = "DNA"
    return sub, sub_start  # also return offset for adjusting feature coords


def add_is_features(sub_rec, offset, is_hits_for_contig):
    """Add IS110 transposase CDS features in the region."""
    sub_len = len(sub_rec.seq)
    for (s, e, strand) in is_hits_for_contig:
        # only if overlaps region
        if e <= offset or s >= offset + sub_len:
            continue
        s_rel = max(0, s - offset)
        e_rel = min(sub_len, e - offset)
        strand_int = 1 if strand == "+" else -1
        feat = SeqFeature(FeatureLocation(s_rel, e_rel, strand=strand_int),
                          type="CDS",
                          qualifiers={"product": ["IS110 transposase"],
                                      "gene": ["IS110"],
                                      "locus_tag": [f"IS110_{s}"]})
        sub_rec.features.append(feat)


def add_misc_feature(sub_rec, offset, start, end, note, ftype="misc_feature"):
    sub_len = len(sub_rec.seq)
    if end <= offset or start >= offset + sub_len:
        return
    s_rel = max(0, start - offset)
    e_rel = min(sub_len, end - offset)
    feat = SeqFeature(FeatureLocation(s_rel, e_rel),
                      type=ftype, qualifiers={"note": [note]})
    sub_rec.features.append(feat)


def process_inversion(case_n, case, paths, out_dir, window):
    """Inversion: extract ref ±window around IS, tgt anchor union ±window."""
    species = case['species']
    ref_id  = case['ref_id']
    tgt_asm = case['tgt_asm']
    rev_bp  = int(case['rev_block_bp'])
    rev_id  = float(case['rev_block_ident'])
    side    = case['best_side']

    case_dir = os.path.join(out_dir, "inversions",
                            f"case_{case_n}_{species[:5]}_{ref_id.split('|')[0]}_vs_{tgt_asm}")
    os.makedirs(case_dir, exist_ok=True)

    ref_asm = ref_id.split("|")[0]
    runs = paths['runs']
    records_path = os.path.join(runs, species, "records", "records.json")
    is_hits_path = os.path.join(runs, species, "is_hits.tsv")
    src_dir = paths['src'][species]
    locus = load_locus(records_path, ref_id)
    if not locus:
        print(f"  case_{case_n}: locus not found")
        return
    ref_contig, is_s, is_e = locus
    is_hits = load_is_hits(is_hits_path)
    ref_seqs = load_assembly_seqs(src_dir, ref_asm)
    tgt_seqs = load_assembly_seqs(src_dir, tgt_asm)
    if ref_contig not in ref_seqs:
        print(f"  case_{case_n}: ref contig {ref_contig} not found")
        return

    # REF: focal window around IS
    ref_rec = ref_seqs[ref_contig]
    rsx = max(0, is_s - window)
    rex = min(len(ref_rec.seq), is_e + window)
    ref_sub, ref_off = make_region_record(ref_rec, rsx, rex, "REF")
    add_is_features(ref_sub, ref_off, is_hits.get(ref_asm, {}).get(ref_contig, []))
    add_misc_feature(ref_sub, ref_off, is_s, is_e,
                     f"Reference IS110 (id={ref_id})", ftype="repeat_region")
    add_misc_feature(ref_sub, ref_off,
                     max(0, is_s - window if side == "up" else is_e),
                     min(len(ref_rec.seq), is_s if side == "up" else is_e + window),
                     f"Inversion MIDDLE-proxy ({side}-side, 30 kb)", ftype="misc_feature")

    # TGT: extract from target contig union of anchors
    anchors = parse_anchors(case.get('details', ''))
    if not anchors:
        # synthesize from details column of the all-candidates table
        cand_path = os.path.join(runs, species, "inversion_candidates_dedup.tsv")
        with open(cand_path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("ref_id") == ref_id:
                    anchors = parse_anchors(r.get("details", ""))
                    break
    if anchors:
        # Pick the contig with the most anchors
        bycontig = defaultdict(list)
        for (c, s, e, st) in anchors:
            bycontig[c].append((s, e, st))
        tgt_contig = max(bycontig.keys(), key=lambda k: len(bycontig[k]))
        tgt_rec = tgt_seqs.get(tgt_contig)
        if tgt_rec is None:
            print(f"  case_{case_n}: tgt contig {tgt_contig} not found")
            return
        s_min = min(s for s, _, _ in bycontig[tgt_contig])
        e_max = max(e for _, e, _ in bycontig[tgt_contig])
        tsx = max(0, s_min - window)
        tex = min(len(tgt_rec.seq), e_max + window)
        tgt_sub, tgt_off = make_region_record(tgt_rec, tsx, tex, "TGT")
        add_is_features(tgt_sub, tgt_off, is_hits.get(tgt_asm, {}).get(tgt_contig, []))
        for (s, e, st) in bycontig[tgt_contig]:
            add_misc_feature(tgt_sub, tgt_off, s, e,
                             f"anchor hit ({st} strand) — inversion endpoint",
                             ftype="misc_feature")
    else:
        tgt_sub = None

    with open(os.path.join(case_dir, "ref.gbk"), "w") as fh:
        SeqIO.write([ref_sub], fh, "genbank")
    if tgt_sub:
        with open(os.path.join(case_dir, "tgt.gbk"), "w") as fh:
            SeqIO.write([tgt_sub], fh, "genbank")
    with open(os.path.join(case_dir, "README.md"), "w") as fh:
        fh.write(f"# Inversion case {case_n}\n\n")
        fh.write(f"- species: {species}\n")
        fh.write(f"- ref_id: `{ref_id}`\n")
        fh.write(f"- ref_asm: `{ref_asm}` contig `{ref_contig}` IS at {is_s}-{is_e}\n")
        fh.write(f"- tgt_asm: `{tgt_asm}`\n")
        fh.write(f"- best_side (MIDDLE proxy): {side}\n")
        fh.write(f"- rev_block_bp: {rev_bp} ({rev_bp/300:.0f}% of 30 kb proxy)\n")
        fh.write(f"- rev_block_ident: {rev_id:.1f}%\n")
        if anchors:
            fh.write(f"- anchor positions in tgt: {len(anchors)} hits (see misc_features in tgt.gbk)\n")
    print(f"  case_{case_n}: wrote {case_dir}")


def process_translocation(case_n, case, paths, out_dir, window):
    species = case['species']
    ref_id  = case['ref_id']
    tgt_asm = case['tgt_asm']
    depth   = case.get('min_anchor_depth', '?')

    case_dir = os.path.join(out_dir, "translocations",
                            f"case_{case_n}_{species[:5]}_{ref_id.split('|')[0]}_vs_{tgt_asm}")
    os.makedirs(case_dir, exist_ok=True)

    ref_asm = ref_id.split("|")[0]
    runs = paths['runs']
    records_path = os.path.join(runs, species, "records", "records.json")
    is_hits_path = os.path.join(runs, species, "is_hits.tsv")
    src_dir = paths['src'][species]
    locus = load_locus(records_path, ref_id)
    if not locus:
        print(f"  case_{case_n}: locus not found"); return
    ref_contig, is_s, is_e = locus
    is_hits = load_is_hits(is_hits_path)
    ref_seqs = load_assembly_seqs(src_dir, ref_asm)
    tgt_seqs = load_assembly_seqs(src_dir, tgt_asm)
    if ref_contig not in ref_seqs:
        print(f"  case_{case_n}: ref contig not found"); return

    ref_rec = ref_seqs[ref_contig]
    rsx = max(0, is_s - window)
    rex = min(len(ref_rec.seq), is_e + window)
    ref_sub, ref_off = make_region_record(ref_rec, rsx, rex, "REF")
    add_is_features(ref_sub, ref_off, is_hits.get(ref_asm, {}).get(ref_contig, []))
    add_misc_feature(ref_sub, ref_off, is_s, is_e,
                     f"Reference IS110 (translocation candidate)", ftype="repeat_region")

    anchors = parse_anchors(case.get('details', ''))
    tgt_records = []
    if anchors:
        bycontig = defaultdict(list)
        for (c, s, e, st) in anchors:
            bycontig[c].append((s, e, st))
        for tgt_contig, hits in bycontig.items():
            tgt_rec = tgt_seqs.get(tgt_contig)
            if tgt_rec is None: continue
            s_min = min(s for s, _, _ in hits)
            e_max = max(e for _, e, _ in hits)
            tsx = max(0, s_min - window)
            tex = min(len(tgt_rec.seq), e_max + window)
            tgt_sub, tgt_off = make_region_record(tgt_rec, tsx, tex, "TGT")
            add_is_features(tgt_sub, tgt_off,
                            is_hits.get(tgt_asm, {}).get(tgt_contig, []))
            for (s, e, st) in hits:
                add_misc_feature(tgt_sub, tgt_off, s, e,
                                 f"anchor hit ({st}) — translocation endpoint",
                                 ftype="misc_feature")
            tgt_records.append(tgt_sub)

    with open(os.path.join(case_dir, "ref.gbk"), "w") as fh:
        SeqIO.write([ref_sub], fh, "genbank")
    if tgt_records:
        with open(os.path.join(case_dir, "tgt.gbk"), "w") as fh:
            SeqIO.write(tgt_records, fh, "genbank")
    with open(os.path.join(case_dir, "README.md"), "w") as fh:
        fh.write(f"# Translocation case {case_n}\n\n")
        fh.write(f"- species: {species}\n")
        fh.write(f"- ref_id: `{ref_id}`\n")
        fh.write(f"- ref_asm: `{ref_asm}` contig `{ref_contig}` IS at {is_s}-{is_e}\n")
        fh.write(f"- tgt_asm: `{tgt_asm}` (chromosome-level / complete)\n")
        fh.write(f"- min_anchor_depth: {depth} bp\n")
        if anchors:
            tgt_contigs = sorted({c for c,_,_,_ in anchors})
            fh.write(f"- tgt contigs hit: {tgt_contigs} (translocation = anchor pairs on different contigs)\n")
    print(f"  case_{case_n}: wrote {case_dir}")


def process_duplication(case_n, case, paths, out_dir, window):
    species = case['species']
    ref_id  = case['ref_id']
    tgt_asm = case['tgt_asm']
    side    = case['best_side']
    n_flank = case['n_is_flanked']
    n_pos   = case['n_anchor_positions']
    n_blk   = case['n_fwd_blocks']
    best_bp = case['best_block_bp']

    case_dir = os.path.join(out_dir, "duplications",
                            f"case_{case_n}_{species[:5]}_{ref_id.split('|')[0]}_vs_{tgt_asm}")
    os.makedirs(case_dir, exist_ok=True)

    ref_asm = ref_id.split("|")[0]
    runs = paths['runs']
    records_path = os.path.join(runs, species, "records", "records.json")
    is_hits_path = os.path.join(runs, species, "is_hits.tsv")
    src_dir = paths['src'][species]
    locus = load_locus(records_path, ref_id)
    if not locus:
        print(f"  case_{case_n}: locus not found"); return
    ref_contig, is_s, is_e = locus
    is_hits = load_is_hits(is_hits_path)
    ref_seqs = load_assembly_seqs(src_dir, ref_asm)
    tgt_seqs = load_assembly_seqs(src_dir, tgt_asm)
    if ref_contig not in ref_seqs:
        print(f"  case_{case_n}: ref contig not found"); return

    ref_rec = ref_seqs[ref_contig]
    rsx = max(0, is_s - window)
    rex = min(len(ref_rec.seq), is_e + window)
    ref_sub, ref_off = make_region_record(ref_rec, rsx, rex, "REF")
    add_is_features(ref_sub, ref_off, is_hits.get(ref_asm, {}).get(ref_contig, []))
    add_misc_feature(ref_sub, ref_off, is_s, is_e,
                     f"Reference IS110 (duplication candidate)", ftype="repeat_region")

    # Look up the duplication's anchor positions from rearrangements_duplication_all.tsv
    dup_path = os.path.join(runs, species, "rearrangements_duplication_all.tsv")
    anchors = []
    with open(dup_path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("ref_id") == ref_id and r.get("assembly") == tgt_asm:
                anchors = parse_anchors(r.get("details", ""))
                break

    tgt_records = []
    if anchors:
        bycontig = defaultdict(list)
        for (c, s, e, st) in anchors:
            bycontig[c].append((s, e, st))
        for tgt_contig, hits in bycontig.items():
            tgt_rec = tgt_seqs.get(tgt_contig)
            if tgt_rec is None: continue
            # extend to cover ALL hits + window
            s_min = min(s for s, _, _ in hits)
            e_max = max(e for _, e, _ in hits)
            tsx = max(0, s_min - window)
            tex = min(len(tgt_rec.seq), e_max + window)
            # If the span is HUGE (> 2 Mb), just extract a window around each cluster
            if tex - tsx > 2_000_000:
                # write multiple regions, one per anchor hit
                hits_sorted = sorted(hits, key=lambda x: x[0])
                for j, (s, e, st) in enumerate(hits_sorted):
                    tsx_j = max(0, s - window)
                    tex_j = min(len(tgt_rec.seq), e + window)
                    sub, off = make_region_record(tgt_rec, tsx_j, tex_j, f"TGT{j+1}")
                    add_is_features(sub, off,
                                    is_hits.get(tgt_asm, {}).get(tgt_contig, []))
                    add_misc_feature(sub, off, s, e,
                                     f"anchor hit ({st}) — duplicated copy #{j+1}",
                                     ftype="misc_feature")
                    tgt_records.append(sub)
            else:
                sub, off = make_region_record(tgt_rec, tsx, tex, "TGT")
                add_is_features(sub, off,
                                is_hits.get(tgt_asm, {}).get(tgt_contig, []))
                for (s, e, st) in hits:
                    add_misc_feature(sub, off, s, e,
                                     f"anchor hit ({st}) — duplicated copy",
                                     ftype="misc_feature")
                tgt_records.append(sub)

    with open(os.path.join(case_dir, "ref.gbk"), "w") as fh:
        SeqIO.write([ref_sub], fh, "genbank")
    if tgt_records:
        with open(os.path.join(case_dir, "tgt.gbk"), "w") as fh:
            SeqIO.write(tgt_records, fh, "genbank")
    with open(os.path.join(case_dir, "README.md"), "w") as fh:
        fh.write(f"# Duplication case {case_n}\n\n")
        fh.write(f"- species: {species}\n")
        fh.write(f"- ref_id: `{ref_id}`\n")
        fh.write(f"- ref_asm: `{ref_asm}` contig `{ref_contig}` IS at {is_s}-{is_e}\n")
        fh.write(f"- tgt_asm: `{tgt_asm}`\n")
        fh.write(f"- best_side: {side}\n")
        fh.write(f"- IS-flanked anchor positions: {n_flank}/{n_pos}\n")
        fh.write(f"- forward alignment blocks (MIDDLE→target): {n_blk}\n")
        fh.write(f"- best block bp: {best_bp}\n")
    print(f"  case_{case_n}: wrote {case_dir}")


def main():
    args = parse_args()
    sel = json.load(open(args.selection))
    src_dirs = json.load(open(args.src_dirs))
    paths = {"runs": args.runs_dir, "src": src_dirs}
    os.makedirs(args.out_dir, exist_ok=True)
    print("=== INVERSIONS ===")
    for i, case in enumerate(sel.get("inversions", []), start=1):
        process_inversion(i, case, paths, args.out_dir, args.focal_window)
    print("=== TRANSLOCATIONS ===")
    for i, case in enumerate(sel.get("translocations", []), start=1):
        process_translocation(i, case, paths, args.out_dir, args.focal_window)
    print("=== DUPLICATIONS ===")
    for i, case in enumerate(sel.get("duplications", []), start=1):
        process_duplication(i, case, paths, args.out_dir, args.focal_window)


if __name__ == "__main__":
    main()
