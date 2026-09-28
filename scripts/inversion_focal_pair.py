#!/usr/bin/env python3
"""Generate paired GenBank files for one IS-mediated inversion candidate,
with EXACTLY 5 features per record:

  GENOME 1 (ref):  [region_A] [IS_a] [inversion_DNA →] [IS_b] [region_B]
  GENOME 2 (tgt):  [region_A] [IS_a] [inversion_DNA ←] [IS_b] [region_B]

Key points
- Detect tgt slice orientation vs ref. If tgt is the revcomp of ref at this
  locus (= the two assemblies stored the same chromosome in opposite strands),
  reverse-complement tgt before slicing so the conserved flanks read in the
  same direction.
- Define region_A / region_B from the ACTUAL forward-strand alignment between
  ref and tgt — not from arbitrary flank windows. They're trimmed to the part
  that truly matches.
- inversion_DNA = the region between the two ISes (= the segment that aligns
  reverse-strand between the two genomes).
"""
import argparse, csv, gzip, json, os, re, shutil, subprocess, tempfile
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--species",  required=True)
    p.add_argument("--ref-id",   required=True)
    p.add_argument("--tgt-asm",  required=True)
    p.add_argument("--runs-dir", required=True)
    p.add_argument("--src-dir",  required=True)
    p.add_argument("--out-dir",  required=True)
    p.add_argument("--flank",    type=int, default=5000,
                   help="bp of flank to initially extract on each side of the IS pair")
    return p.parse_args()


ANCHOR_RE = re.compile(r"([\w.|]+):(\d+)-(\d+)\(([+-])\)")


def load_is_hits(path):
    idx = defaultdict(lambda: defaultdict(list))
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            asm = r["assembly"]; contig = r["contig"]
            local = contig.split("|", 1)[1] if "|" in contig else contig
            hit = (int(r["tnp_start"]), int(r["tnp_end"]),
                   r.get("tnp_strand", "+"))
            idx[asm][local].append(hit)
            tail = local.split("|")[-1]
            if tail != local:
                idx[asm][tail].append(hit)
    return idx


def load_locus(records_path, ref_id):
    with open(records_path) as f:
        for r in json.load(f):
            if (r.get("ref_id") or r.get("is110_id")) != ref_id: continue
            src = r.get("source", {})
            ie  = r.get("is_element") or src.get("is_element") or {}
            contig = src.get("contig")
            local = contig.split("|", 1)[1] if "|" in contig else contig
            return (local,
                    int(ie.get("source_start") or src.get("transposase_start")),
                    int(ie.get("source_end")   or src.get("transposase_end")))
    return None


def load_seqs(src_dir, asm):
    for ext in (".fna.gz", ".fa.gz", ".fna", ".fa"):
        path = os.path.join(src_dir, asm + ext)
        if os.path.exists(path):
            opener = gzip.open if path.endswith(".gz") else open
            recs = {}
            with opener(path, "rt") as fh:
                for rec in SeqIO.parse(fh, "fasta"):
                    local = rec.id.split("|", 1)[1] if "|" in rec.id else rec.id
                    recs[local] = rec
                    recs[rec.id] = rec
            return recs
    return {}


def merge_is(hits, gap=500):
    if not hits: return []
    sh = sorted(hits, key=lambda x: x[0])
    merged = [sh[0]]
    for h in sh[1:]:
        if h[0] <= merged[-1][1] + gap:
            merged[-1] = (min(merged[-1][0], h[0]),
                          max(merged[-1][1], h[1]),
                          merged[-1][2])
        else:
            merged.append(h)
    return merged


def find_nearest_is(is_list, pos, max_dist=5000):
    best = None
    for (s, e, st) in is_list:
        d = 0 if s <= pos <= e else min(abs(pos - s), abs(pos - e))
        if d <= max_dist and (best is None or d < best[3]):
            best = (s, e, st, d)
    return best


def write_fasta(label, seq, path):
    with open(path, "w") as fh:
        fh.write(f">{label}\n{seq}\n")


def align_paf(ref_path, query_path, work, threads=4):
    paf = os.path.join(work, "aln.paf")
    subprocess.run(["minimap2", "-x", "asm10", "-c", "--eqx",
                    "-t", str(threads), ref_path, query_path, "-o", paf],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    blocks = []
    with open(paf) as fh:
        for line in fh:
            c = line.split("\t")
            if len(c) < 12: continue
            blocks.append({
                "qname": c[0], "qlen": int(c[1]), "qs": int(c[2]), "qe": int(c[3]),
                "strand": c[4],
                "tname": c[5], "tlen": int(c[6]), "ts": int(c[7]), "te": int(c[8]),
                "matches": int(c[9]), "block": int(c[10]),
            })
    return blocks


def is_flipped(blocks, ref_len):
    """Decide if tgt is reverse-complement of ref overall, by comparing total
    forward-strand alignment vs reverse-strand alignment of the EDGES (flanks).
    If the edge-flanks align REVERSE, tgt is flipped."""
    edge_pad = ref_len // 4
    fwd_edge = sum(b["block"] for b in blocks if b["strand"] == "+"
                   and (b["ts"] < edge_pad or b["te"] > ref_len - edge_pad))
    rev_edge = sum(b["block"] for b in blocks if b["strand"] == "-"
                   and (b["ts"] < edge_pad or b["te"] > ref_len - edge_pad))
    return rev_edge > fwd_edge


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    runs = args.runs_dir
    records = os.path.join(runs, args.species, "records", "records.json")
    is_hits = load_is_hits(os.path.join(runs, args.species, "is_hits.tsv"))
    locus = load_locus(records, args.ref_id)
    if not locus:
        raise SystemExit(f"locus not found for {args.ref_id}")
    ref_contig, ref_is_s, ref_is_e = locus

    cand_path = os.path.join(runs, args.species, "inversion_candidates_dedup.tsv")
    details = None
    for r in csv.DictReader(open(cand_path), delimiter="\t"):
        if r.get("ref_id") == args.ref_id and r.get("assembly") == args.tgt_asm:
            details = r.get("details", ""); break
    if not details: raise SystemExit("anchor details not found")
    anchors = [(m.group(1), int(m.group(2)), int(m.group(3)), m.group(4))
               for m in ANCHOR_RE.finditer(details)]
    bycontig = defaultdict(list)
    for (c, s, e, st) in anchors:
        bycontig[c].append((s, e, st))
    tgt_contig = max(bycontig.keys(), key=lambda k: len(bycontig[k]))
    tgt_hits = sorted(bycontig[tgt_contig], key=lambda x: x[0])
    h_left, h_right = tgt_hits[0], tgt_hits[-1]
    span = h_right[1] - h_left[0]

    # ---- Load ref and tgt sequences ----
    ref_asm = args.ref_id.split("|")[0]
    ref_seqs = load_seqs(args.src_dir, ref_asm)
    tgt_seqs = load_seqs(args.src_dir, args.tgt_asm)
    ref_rec = ref_seqs.get(ref_contig)
    tgt_rec = tgt_seqs.get(tgt_contig)
    if ref_rec is None or tgt_rec is None:
        raise SystemExit(f"sequence load failed: ref={ref_contig in ref_seqs} tgt={tgt_contig in tgt_seqs}")

    # ---- Find IS pair in target ----
    tgt_is_merged = merge_is(is_hits.get(args.tgt_asm, {}).get(tgt_contig, []))
    is_at_left  = find_nearest_is(tgt_is_merged, h_left[0],  max_dist=5000)
    is_at_right = find_nearest_is(tgt_is_merged, h_right[1], max_dist=5000)
    if is_at_left is None or is_at_right is None:
        raise SystemExit("could not locate target IS pair")
    tgt_is_a = (is_at_left[0],  is_at_left[1],  is_at_left[2])
    tgt_is_b = (is_at_right[0], is_at_right[1], is_at_right[2])
    if tgt_is_a[0] > tgt_is_b[0]:
        tgt_is_a, tgt_is_b = tgt_is_b, tgt_is_a

    # ---- Find IS pair in ref ----
    ref_is_merged = merge_is(is_hits.get(ref_asm, {}).get(ref_contig, []))
    is_a_ref = next(((s,e,st) for (s,e,st) in ref_is_merged
                     if s <= ref_is_s <= e or abs(s - ref_is_s) < 200), None)
    if is_a_ref is None:
        is_a_ref = (ref_is_s, ref_is_e, "+")
    partner = None
    for (s, e, st) in ref_is_merged:
        if (s, e) == (is_a_ref[0], is_a_ref[1]): continue
        d = abs(s - is_a_ref[0])
        if abs(d - span) < 5000:
            partner = (s, e, st); break
    if partner is None:
        nearest = sorted(((s,e,st,abs(s-is_a_ref[0])) for (s,e,st) in ref_is_merged
                          if (s,e) != (is_a_ref[0], is_a_ref[1])), key=lambda x: x[3])
        if nearest: partner = (nearest[0][0], nearest[0][1], nearest[0][2])
    if partner is None: raise SystemExit("no partner IS in ref")
    is_b_ref = partner
    if is_a_ref[0] > is_b_ref[0]:
        is_a_ref, is_b_ref = is_b_ref, is_a_ref

    # ---- Extract initial slices ----
    ref_lo = max(0, is_a_ref[0] - args.flank)
    ref_hi = min(len(ref_rec.seq), is_b_ref[1] + args.flank)
    tgt_lo = max(0, tgt_is_a[0] - args.flank)
    tgt_hi = min(len(tgt_rec.seq), tgt_is_b[1] + args.flank)
    ref_seq = str(ref_rec.seq[ref_lo:ref_hi])
    tgt_seq = str(tgt_rec.seq[tgt_lo:tgt_hi])

    # ---- Detect tgt orientation vs ref; revcomp if needed ----
    work = tempfile.mkdtemp(prefix="invpair_")
    write_fasta("REF", ref_seq, os.path.join(work, "ref.fa"))
    write_fasta("TGT", tgt_seq, os.path.join(work, "tgt.fa"))
    blocks = align_paf(os.path.join(work, "ref.fa"),
                       os.path.join(work, "tgt.fa"), work)
    tgt_flipped = is_flipped(blocks, len(ref_seq))
    print(f"tgt orientation vs ref: {'FLIPPED (revcomp needed)' if tgt_flipped else 'same'}")
    if tgt_flipped:
        # revcomp the tgt slice and remap the IS positions to revcomp slice coords
        from Bio.Seq import Seq
        tgt_seq = str(Seq(tgt_seq).reverse_complement())
        L = tgt_hi - tgt_lo
        # original slice-relative positions
        a_s_rel = tgt_is_a[0] - tgt_lo; a_e_rel = tgt_is_a[1] - tgt_lo
        b_s_rel = tgt_is_b[0] - tgt_lo; b_e_rel = tgt_is_b[1] - tgt_lo
        # On revcomp, [s, e) on slice becomes [L - e, L - s). Original IS_b (rightmost)
        # becomes the new leftmost = IS_a on revcomp; original IS_a becomes new IS_b.
        tgt_is_a_rel = (L - b_e_rel, L - b_s_rel,
                        "+" if tgt_is_b[2] == "-" else "-")
        tgt_is_b_rel = (L - a_e_rel, L - a_s_rel,
                        "+" if tgt_is_a[2] == "-" else "-")
        write_fasta("TGT", tgt_seq, os.path.join(work, "tgt.fa"))
        blocks = align_paf(os.path.join(work, "ref.fa"),
                           os.path.join(work, "tgt.fa"), work)
    else:
        tgt_is_a_rel = (tgt_is_a[0] - tgt_lo, tgt_is_a[1] - tgt_lo, tgt_is_a[2])
        tgt_is_b_rel = (tgt_is_b[0] - tgt_lo, tgt_is_b[1] - tgt_lo, tgt_is_b[2])

    # IS positions in ref (relative to slice)
    ref_is_a_rel = (is_a_ref[0] - ref_lo, is_a_ref[1] - ref_lo, is_a_ref[2])
    ref_is_b_rel = (is_b_ref[0] - ref_lo, is_b_ref[1] - ref_lo, is_b_ref[2])

    # ---- Re-derive region_A / region_B from REAL alignment ----
    # After (possibly) revcomping, find:
    #   - largest FORWARD block at the 5' edge (region_A boundary)
    #   - largest FORWARD block at the 3' edge (region_B boundary)
    #   - largest REVERSE block in the middle (= inversion_DNA boundary in target)
    # Then trim region_A end / region_B start to the actual alignment boundaries.
    fwd = [b for b in blocks if b["strand"] == "+"]
    rev = [b for b in blocks if b["strand"] == "-"]
    if not fwd:
        raise SystemExit("no forward-strand alignment between ref and tgt — wrong case")
    # Sort fwd blocks by ts position
    fwd.sort(key=lambda b: b["ts"])
    # region_A boundary: the forward block whose ts < ref_is_a position
    fwd_left  = max((b for b in fwd if b["ts"] < ref_is_a_rel[0]),
                    key=lambda b: b["block"], default=None)
    fwd_right = max((b for b in fwd if b["te"] > ref_is_b_rel[1]),
                    key=lambda b: b["block"], default=None)
    rev_mid   = max((b for b in rev if b["ts"] >= ref_is_a_rel[0] - 500
                                  and b["te"] <= ref_is_b_rel[1] + 500),
                    key=lambda b: b["block"], default=None)

    # ---- Build 5-feature records ----
    def make_5feature(seq, label_prefix, is_a, is_b, fwd_l, fwd_r, rev_m, src_id, src_start):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        sub = SeqRecord(Seq(seq), id=f"{label_prefix}_{src_id[:10]}"[:16],
                        name=f"{label_prefix}_{src_id[:10]}"[:16],
                        description=f"{src_id} slice {src_start+1}-{src_start+len(seq)} ({label_prefix})")
        sub.annotations["molecule_type"] = "DNA"
        L = len(seq)
        # region_A bounds
        rA_end = fwd_l["te"] if fwd_l else is_a[0]
        if rA_end <= 0: rA_end = max(1, is_a[0])
        rA_end = min(rA_end, is_a[0])  # never overlap IS
        sub.features.append(SeqFeature(
            FeatureLocation(0, rA_end), type="misc_feature",
            qualifiers={"label": ["region_A"],
                        "note": [f"region_A — conserved 5' flank ({rA_end} bp aligned)"]}))
        # IS_a
        sub.features.append(SeqFeature(
            FeatureLocation(is_a[0], is_a[1],
                            strand=1 if is_a[2] == "+" else -1),
            type="repeat_region",
            qualifiers={"label": ["IS_a"], "rpt_family": ["IS110"],
                        "note": ["IS110 transposase (IS_a)"]}))
        # inversion_DNA
        invDNA_start = is_a[1]
        invDNA_end   = is_b[0]
        arrow = "forward" if label_prefix == "REF" else "reverse-complement (vs REF)"
        rev_note = ""
        if rev_m:
            rev_note = f" — minimap2 reverse block: t[{rev_m['ts']}-{rev_m['te']}] vs q[{rev_m['qs']}-{rev_m['qe']}] {rev_m['block']} bp"
        sub.features.append(SeqFeature(
            FeatureLocation(invDNA_start, invDNA_end), type="misc_feature",
            qualifiers={"label": ["inversion_DNA"],
                        "note": [f"inversion_DNA ({arrow}){rev_note}"]}))
        # IS_b
        sub.features.append(SeqFeature(
            FeatureLocation(is_b[0], is_b[1],
                            strand=1 if is_b[2] == "+" else -1),
            type="repeat_region",
            qualifiers={"label": ["IS_b"], "rpt_family": ["IS110"],
                        "note": ["IS110 transposase (IS_b)"]}))
        # region_B
        rB_start = fwd_r["ts"] if fwd_r else is_b[1]
        rB_start = max(rB_start, is_b[1])
        sub.features.append(SeqFeature(
            FeatureLocation(rB_start, L), type="misc_feature",
            qualifiers={"label": ["region_B"],
                        "note": [f"region_B — conserved 3' flank ({L - rB_start} bp aligned)"]}))
        return sub

    # In the ref slice, the forward/reverse alignments are with respect to tgt.
    # We use the same fwd_l, fwd_r, rev_m boundaries on the ref slice (ref coords)
    # AND map them to tgt coords via the alignment.
    ref_sub = make_5feature(ref_seq, "REF",
                            ref_is_a_rel, ref_is_b_rel,
                            fwd_l=fwd_left, fwd_r=fwd_right, rev_m=rev_mid,
                            src_id=ref_contig, src_start=ref_lo)

    # For tgt: map the alignment boundaries to tgt coordinates
    # fwd_left.qs/qe → tgt coords of the left forward block (region_A end in tgt)
    # fwd_right.qs/qe → tgt coords of the right forward block (region_B start in tgt)
    class T: pass
    tgt_fwd_left  = T(); tgt_fwd_right = T(); tgt_rev_mid   = T()
    tgt_fwd_left.te  = fwd_left["qe"]  if fwd_left  else tgt_is_a_rel[0]
    tgt_fwd_right.ts = fwd_right["qs"] if fwd_right else tgt_is_b_rel[1]
    tgt_rev_mid_dict = None
    if rev_mid:
        tgt_rev_mid_dict = {"ts": rev_mid["qs"], "te": rev_mid["qe"],
                            "qs": rev_mid["ts"], "qe": rev_mid["te"],
                            "block": rev_mid["block"]}
    # use simple proxies for label_writer
    class FL:
        def __init__(self, te): self.te = te
    class FR:
        def __init__(self, ts): self.ts = ts
    fwd_left_t  = {"te": tgt_fwd_left.te,  "block": fwd_left["block"]  if fwd_left  else 0,
                   "ts": fwd_left["qs"]   if fwd_left  else 0,
                   "qs": fwd_left["ts"]   if fwd_left  else 0,
                   "qe": fwd_left["te"]   if fwd_left  else 0}
    fwd_right_t = {"ts": tgt_fwd_right.ts, "block": fwd_right["block"] if fwd_right else 0,
                   "te": fwd_right["qe"]  if fwd_right else 0,
                   "qs": fwd_right["ts"]  if fwd_right else 0,
                   "qe": fwd_right["te"]  if fwd_right else 0}
    tgt_sub = make_5feature(tgt_seq, "TGT",
                            tgt_is_a_rel, tgt_is_b_rel,
                            fwd_l=fwd_left_t, fwd_r=fwd_right_t, rev_m=tgt_rev_mid_dict,
                            src_id=tgt_contig, src_start=tgt_lo)

    with open(os.path.join(args.out_dir, "ref.gbk"), "w") as fh:
        SeqIO.write([ref_sub], fh, "genbank")
    with open(os.path.join(args.out_dir, "tgt.gbk"), "w") as fh:
        SeqIO.write([tgt_sub], fh, "genbank")
    print(f"\nWrote {args.out_dir}/ref.gbk ({len(ref_seq):,} bp)")
    print(f"Wrote {args.out_dir}/tgt.gbk ({len(tgt_seq):,} bp)  {'(reverse-complemented from source)' if tgt_flipped else ''}")
    shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
