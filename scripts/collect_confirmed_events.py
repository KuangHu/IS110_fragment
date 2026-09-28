#!/usr/bin/env python3
"""Collect CONFIRMED inversion / duplication / translocation events into a
publish-ready folder structure.

Writes VCF v4.3 (the de-facto standard, accepted by dbVar / DGVa) plus a flat
TSV with IS-pipeline metadata plus per-event evidence packs containing the
reference locus FASTA, target region FASTA, dnadiff .report, and progressiveMauve
XMFA — everything needed to reproduce the verdict.

CONFIRMED set (per category):
  inversion_only      : verdict == CONFIRMED_BOTH   (dnadiff + Mauve agree)
  duplication         : verdict == CONFIRMED_DNADIFF (Mauve n/a for dups)
  translocation_only  : verdict == CONFIRMED_BOTH AND deep-anchor (not artifact)

Inputs (per species):
  <results>/<sp>/rearr_validate/validation.tsv      — validator output
  <results>/<sp>/records/records.json               — reference IS loci
  <results>/<sp>/rearrangements_examples.tsv        — anchor positions / details
  <by_species>/<sp>/<assembly>.fna.gz               — assembly source

Usage:
  collect_confirmed_events.py --runs <results_dir> --src <by_species_dir> \\
      --out <output_dir>  [--species ...]
"""
import argparse, csv, gzip, json, os, re, shutil, subprocess, sys, tempfile
from collections import defaultdict

import pysam


CONFIRMED = {
    "inversion_only":     {"CONFIRMED_BOTH"},
    "duplication":        {"CONFIRMED_DNADIFF", "CONFIRMED_BOTH"},
    "translocation_only": {"CONFIRMED_BOTH"},   # already excludes contig-break artifacts
}
SVTYPE = {"inversion_only": "INV", "duplication": "DUP",
          "translocation_only": "BND"}
EVENT_PREFIX = {"inversion_only": "INV", "duplication": "DUP",
                "translocation_only": "TRANSLOC"}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--runs", required=True, help="pipeline results dir (parent of <sp>/)")
    p.add_argument("--src", required=True, help="by_species dir with <asm>.fna.gz")
    p.add_argument("--out", required=True)
    p.add_argument("--species", nargs="*",
                   help="restrict to these species (default: all that have validation.tsv)")
    p.add_argument("--locus-window", type=int, default=30000,
                   help="bp of context around the IS locus written into evidence")
    p.add_argument("--mauve",
                   default="/global/home/users/kh36969/.conda/envs/mauve-env/bin/progressiveMauve",
                   help="progressiveMauve binary for regenerating evidence")
    p.add_argument("--regen-evidence", action="store_true",
                   help="re-run dnadiff + Mauve on each event and save outputs")
    p.add_argument("--threads", type=int, default=4)
    return p.parse_args()


def load_ref_loci(records_json):
    loci = {}
    with open(records_json) as f:
        recs = json.load(f)
    for r in recs:
        rid = r.get("ref_id") or r.get("is110_id")
        src = r.get("source", {})
        ie = r.get("is_element", {})
        contig = src.get("contig")
        s = ie.get("source_start") or src.get("transposase_start")
        e = ie.get("source_end") or src.get("transposase_end")
        strand = src.get("transposase_strand", "+")
        if rid and contig and s and e:
            local = contig.split("|", 1)[1] if "|" in contig else contig
            loci[rid] = {"contig": local, "start": int(s), "end": int(e), "strand": strand}
    return loci


def load_examples_details(tsv):
    """ref_id+tgt_asm+category -> details_string."""
    out = {}
    with open(tsv) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            key = (row["category"], row["ref_id"], row["assembly"])
            out[key] = row["details"]
    return out


def anchor_targets(details):
    out = []
    for m in re.finditer(r"([\w.|]+):(\d+)-(\d+)", details):
        out.append((m.group(1), int(m.group(2)), int(m.group(3))))
    return out


def decompress_assembly(src_dir, asm, work):
    for ext in (".fna.gz", ".fa.gz", ".fna", ".fa"):
        cand = os.path.join(src_dir, asm + ext)
        if os.path.exists(cand):
            out = os.path.join(work, asm + ".fa")
            if cand.endswith(".gz"):
                with gzip.open(cand, "rt") as fi, open(out, "w") as fo:
                    shutil.copyfileobj(fi, fo)
            else:
                shutil.copy(cand, out)
            return out
    return None


def extract_region(fa, contig, start, end, window, label, out_fh):
    try:
        clen = fa.get_reference_length(contig)
    except (KeyError, ValueError):
        return None
    s = max(0, start - window)
    e = min(clen, end + window)
    seq = fa.fetch(contig, s, e)
    if not seq:
        return None
    out_fh.write(f">{label}\n")
    for i in range(0, len(seq), 80):
        out_fh.write(seq[i:i+80] + "\n")
    return (s, e, clen)


def write_vcf_header(fh, contigs, species):
    fh.write("##fileformat=VCFv4.3\n")
    fh.write("##source=Cross_reference_IS pipeline + dnadiff(MUMmer4) + progressiveMauve\n")
    fh.write(f"##species={species}\n")
    fh.write("##ALT=<ID=INV,Description=\"Inversion\">\n")
    fh.write("##ALT=<ID=DUP,Description=\"Duplication\">\n")
    fh.write("##ALT=<ID=BND,Description=\"Translocation breakend\">\n")
    fh.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
    fh.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the SV\">\n")
    fh.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the SV\">\n")
    fh.write("##INFO=<ID=IS_MEDIATED,Number=0,Type=Flag,Description=\"Event is IS-element-mediated\">\n")
    fh.write("##INFO=<ID=IS_FAMILY,Number=1,Type=String,Description=\"IS element family\">\n")
    fh.write("##INFO=<ID=REF_ID,Number=1,Type=String,Description=\"Cross_reference_IS pipeline ref_id\">\n")
    fh.write("##INFO=<ID=TGT_ASM,Number=1,Type=String,Description=\"Target assembly accession\">\n")
    fh.write("##INFO=<ID=VAL_DNADIFF,Number=1,Type=String,Description=\"MUMmer dnadiff verdict\">\n")
    fh.write("##INFO=<ID=VAL_MAUVE,Number=1,Type=String,Description=\"progressiveMauve verdict\">\n")
    fh.write("##INFO=<ID=MAUVE_METRIC,Number=1,Type=String,Description=\"Mauve LCB metric\">\n")
    fh.write("##INFO=<ID=MIN_ANCHOR_END_DIST,Number=1,Type=Integer,Description=\"Min bp from anchor to target contig end\">\n")
    fh.write("##INFO=<ID=CONFIRMED_BY,Number=.,Type=String,Description=\"Gold-standard tool(s) that independently confirmed this event\">\n")
    fh.write("##INFO=<ID=MATEID,Number=1,Type=String,Description=\"Mate breakend ID (BND only)\">\n")
    fh.write("##FILTER=<ID=PASS,Description=\"Confirmed by both gold-standard tools (or dnadiff-only for DUP)\">\n")
    for name, length in contigs:
        fh.write(f"##contig=<ID={name},length={length}>\n")
    fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")


def collect_for_species(sp, args, out_root):
    runs = args.runs
    sp_dir = os.path.join(runs, sp)
    val_tsv = os.path.join(sp_dir, "rearr_validate", "validation.tsv")
    rec_json = os.path.join(sp_dir, "records", "records.json")
    ex_tsv = os.path.join(sp_dir, "rearrangements_examples.tsv")
    if not (os.path.exists(val_tsv) and os.path.exists(rec_json) and os.path.exists(ex_tsv)):
        print(f"[skip] {sp}: missing inputs", file=sys.stderr)
        return None

    ref_loci = load_ref_loci(rec_json)
    details_by = load_examples_details(ex_tsv)

    sp_out = os.path.join(out_root, sp)
    os.makedirs(sp_out, exist_ok=True)
    ev_dir = os.path.join(sp_out, "evidence")
    os.makedirs(ev_dir, exist_ok=True)

    rows = list(csv.DictReader(open(val_tsv), delimiter="\t"))
    confirmed = [r for r in rows if r["verdict"] in CONFIRMED.get(r["category"], set())]
    print(f"  {sp}: {len(confirmed)} confirmed of {len(rows)} validated", file=sys.stderr)
    if not confirmed:
        return {"species": sp, "n_confirmed": 0, "n_validated": len(rows)}

    src_dir = os.path.join(args.src, sp)
    work = tempfile.mkdtemp(prefix=f"collect_{sp}_", dir=out_root)

    flat = []
    vcf_records = []
    cache_fa = {}

    def open_asm(asm):
        if asm in cache_fa:
            return cache_fa[asm]
        fa_path = decompress_assembly(src_dir, asm, work)
        if not fa_path:
            cache_fa[asm] = None; return None
        cache_fa[asm] = pysam.FastaFile(fa_path)
        cache_fa[f"{asm}__path"] = fa_path
        return cache_fa[asm]

    for idx, r in enumerate(confirmed, 1):
        cat = r["category"]
        ref_id = r["ref_id"]; ref_asm = r["ref_asm"]; tgt_asm = r["tgt_asm"]
        locus = ref_loci.get(ref_id)
        if not locus:
            continue
        details = details_by.get((cat, ref_id, tgt_asm), "")
        anchors = anchor_targets(details)
        event_id = f"{EVENT_PREFIX[cat]}_{idx:04d}"

        # --- per-event evidence pack ---
        ev = os.path.join(ev_dir, event_id)
        os.makedirs(ev, exist_ok=True)
        ref_fa_obj = open_asm(ref_asm)
        qry_fa_obj = open_asm(tgt_asm)
        if not ref_fa_obj or not qry_fa_obj:
            continue
        ref_region_fa = os.path.join(ev, "ref_locus.fa")
        with open(ref_region_fa, "w") as fh:
            ref_span = extract_region(ref_fa_obj, locus["contig"], locus["start"],
                                      locus["end"], args.locus_window, "REF_LOCUS", fh)
        if not ref_span:
            continue
        qry_region_fa = os.path.join(ev, "target_region.fa")
        anchor_records = []
        with open(qry_region_fa, "w") as fh:
            seen_c = {}
            for (c, s, e) in anchors:
                seen_c.setdefault(c, [s, e])
                seen_c[c][0] = min(seen_c[c][0], s); seen_c[c][1] = max(seen_c[c][1], e)
            for c, (s, e) in seen_c.items():
                sp_span = extract_region(qry_fa_obj, c, s, e,
                                         args.locus_window, f"TGT_{c}", fh)
                if sp_span:
                    anchor_records.append({"contig": c, "anchor_start": s,
                                           "anchor_end": e, "extracted_span": sp_span})

        # --- regenerate dnadiff + Mauve evidence (optional) ---
        dnadiff_report = ""
        mauve_xmfa = ""
        if args.regen_evidence:
            dd_pre = os.path.join(ev, "dnadiff")
            try:
                subprocess.run(["dnadiff", "-p", dd_pre, ref_region_fa, qry_region_fa],
                               check=True, stdout=subprocess.DEVNULL,
                               stderr=subprocess.DEVNULL, timeout=300)
                if os.path.exists(dd_pre + ".report"):
                    dnadiff_report = dd_pre + ".report"
            except Exception:
                pass
            try:
                xmfa = os.path.join(ev, "mauve.xmfa")
                subprocess.run([args.mauve, f"--output={xmfa}",
                                ref_region_fa, qry_region_fa],
                               check=True, stdout=subprocess.DEVNULL,
                               stderr=subprocess.DEVNULL, cwd=ev, timeout=900)
                if os.path.exists(xmfa):
                    mauve_xmfa = xmfa
            except Exception:
                pass

        # --- event.json (structured metadata) ---
        ev_meta = {
            "event_id": event_id,
            "species": sp,
            "category": cat,
            "svtype": SVTYPE[cat],
            "ref_assembly": ref_asm,
            "target_assembly": tgt_asm,
            "ref_id": ref_id,
            "is_family": "IS110",
            "is_locus": {"contig": locus["contig"],
                         "start": locus["start"], "end": locus["end"],
                         "strand": locus["strand"]},
            "target_anchors": anchor_records,
            "validation": {
                "dnadiff": r.get("dnadiff"),
                "mauve": r.get("mauve"),
                "mauve_metric": r.get("mauve_metric"),
                "min_anchor_end_dist_bp": int(r.get("min_anchor_end_dist", -1) or -1),
                "verdict": r["verdict"],
                "confirmed_by_software": [
                    n for n, ok in [
                        ("MUMmer4_dnadiff_v4.0.1",  r.get("dnadiff") == "YES"),
                        ("progressiveMauve_v2.4.0", r.get("mauve")   == "YES"),
                    ] if ok
                ],
                "tools": [
                    {"name": "MUMmer4 dnadiff", "version": "4.0.1",
                     "ref": "Marçais et al. 2018, PLOS Comput Biol"},
                    {"name": "progressiveMauve", "version": "2.4.0",
                     "ref": "Darling et al. 2010, PLOS ONE"},
                ],
                "params": {"locus_window_bp": args.locus_window,
                           "contig_break_window_bp": 5000,
                           "min_identity_pct": 95},
            },
        }
        with open(os.path.join(ev, "event.json"), "w") as f:
            json.dump(ev_meta, f, indent=2)

        # --- VCF record ---
        sv_chrom = locus["contig"]
        sv_pos = locus["start"]
        sv_end = locus["end"]
        sv_len = sv_end - sv_pos
        # tools that backed this event (one per "YES" verdict)
        _cby = []
        if r.get("dnadiff") == "YES": _cby.append("MUMmer4_dnadiff_v4.0.1")
        if r.get("mauve")   == "YES": _cby.append("progressiveMauve_v2.4.0")
        cby_vcf = ",".join(_cby) if _cby else "."

        info_common = (f"IS_MEDIATED;IS_FAMILY=IS110;REF_ID={ref_id};"
                       f"TGT_ASM={tgt_asm};VAL_DNADIFF={r.get('dnadiff','')};"
                       f"VAL_MAUVE={r.get('mauve','')};"
                       f"MAUVE_METRIC={r.get('mauve_metric','')};"
                       f"MIN_ANCHOR_END_DIST={r.get('min_anchor_end_dist','')};"
                       f"CONFIRMED_BY={cby_vcf}")

        if cat == "translocation_only":
            # paired BND records, MATEID links them
            for bi, (c, s, e) in enumerate(anchors[:2]):
                m_id = f"{event_id}_BND_{bi+1}"
                mate = f"{event_id}_BND_{2-bi}"
                alt = f"N[{c}:{s}["  # generic breakend; direction simplified
                info = (f"SVTYPE=BND;{info_common};MATEID={mate}")
                vcf_records.append((sv_chrom, sv_pos + bi, m_id, "N", alt, ".", "PASS", info))
        else:
            vid = event_id
            alt = f"<{SVTYPE[cat]}>"
            info = (f"SVTYPE={SVTYPE[cat]};END={sv_end};SVLEN={sv_len};{info_common}")
            vcf_records.append((sv_chrom, sv_pos, vid, "N", alt, ".", "PASS", info))

        # which gold-standard tools backed this event
        confirmed_by = []
        if r.get("dnadiff") == "YES":
            confirmed_by.append("MUMmer4_dnadiff_v4.0.1")
        if r.get("mauve") == "YES":
            confirmed_by.append("progressiveMauve_v2.4.0")
        confirmed_by_str = ";".join(confirmed_by) if confirmed_by else ""

        # --- flat TSV row ---
        flat.append({
            "event_id": event_id, "species": sp, "category": cat,
            "svtype": SVTYPE[cat], "ref_id": ref_id, "ref_asm": ref_asm,
            "tgt_asm": tgt_asm, "is_contig": locus["contig"],
            "is_start": locus["start"], "is_end": locus["end"],
            "is_strand": locus["strand"],
            "target_anchor_contigs": ";".join(a["contig"] for a in anchor_records),
            "n_target_contigs": len(anchor_records),
            "val_dnadiff": r.get("dnadiff", ""),
            "val_mauve": r.get("mauve", ""),
            "mauve_metric": r.get("mauve_metric", ""),
            "min_anchor_end_dist_bp": r.get("min_anchor_end_dist", ""),
            "verdict": r["verdict"],
            "confirmed_by_software": confirmed_by_str,
            "evidence_dir": os.path.relpath(ev, sp_out),
        })

    # --- write outputs for this species ---
    # contigs touched (for VCF header)
    contigs = sorted({(r[0], 0) for r in vcf_records})
    vcf_path = os.path.join(sp_out, "confirmed.vcf")
    with open(vcf_path, "w") as fh:
        write_vcf_header(fh, contigs, sp)
        for rec in sorted(vcf_records, key=lambda x: (x[0], x[1])):
            fh.write("\t".join(str(x) for x in rec) + "\n")

    tsv_path = os.path.join(sp_out, "confirmed.tsv")
    if flat:
        with open(tsv_path, "w") as fh:
            cols = list(flat[0].keys())
            fh.write("\t".join(cols) + "\n")
            for row in flat:
                fh.write("\t".join(str(row[c]) for c in cols) + "\n")

    # close pysam handles + cleanup tmp
    for k, v in cache_fa.items():
        if hasattr(v, "close"):
            v.close()
    shutil.rmtree(work, ignore_errors=True)

    return {"species": sp, "n_confirmed": len(flat), "n_validated": len(rows),
            "by_category": {cat: sum(1 for r in flat if r["category"] == cat)
                            for cat in CONFIRMED}}


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    all_species = sorted(d for d in os.listdir(args.runs)
                         if os.path.isdir(os.path.join(args.runs, d, "rearr_validate")))
    target = args.species or all_species
    print(f"Collecting from {len(target)} species ...", file=sys.stderr)

    summary = []
    for sp in target:
        s = collect_for_species(sp, args, args.out)
        if s:
            summary.append(s)

    # all-species summary TSV
    summ_path = os.path.join(args.out, "all_species_summary.tsv")
    with open(summ_path, "w") as f:
        f.write("species\tn_validated\tn_confirmed\tn_inversion\tn_duplication\tn_translocation\n")
        for s in summary:
            bc = s.get("by_category", {})
            f.write(f"{s['species']}\t{s['n_validated']}\t{s['n_confirmed']}\t"
                    f"{bc.get('inversion_only',0)}\t{bc.get('duplication',0)}\t"
                    f"{bc.get('translocation_only',0)}\n")

    # combined VCF + TSV
    all_vcf = os.path.join(args.out, "all_confirmed.vcf")
    all_tsv = os.path.join(args.out, "all_confirmed.tsv")
    with open(all_vcf, "w") as fout:
        fout.write("##fileformat=VCFv4.3\n")
        fout.write("##source=Cross_reference_IS pipeline (dnadiff+progressiveMauve)\n")
        fout.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for s in summary:
            sp = s["species"]
            spv = os.path.join(args.out, sp, "confirmed.vcf")
            if not os.path.exists(spv): continue
            for line in open(spv):
                if line.startswith("#"): continue
                fout.write(line)
    with open(all_tsv, "w") as fout:
        header_written = False
        for s in summary:
            sp = s["species"]
            spt = os.path.join(args.out, sp, "confirmed.tsv")
            if not os.path.exists(spt): continue
            with open(spt) as fin:
                hdr = fin.readline()
                if not header_written:
                    fout.write(hdr); header_written = True
                for line in fin:
                    fout.write(line)

    print(f"\nWrote: {summ_path}\n       {all_vcf}\n       {all_tsv}", file=sys.stderr)
    for s in summary:
        bc = s.get("by_category", {})
        print(f"  {s['species']:<28} confirmed={s['n_confirmed']:<4} "
              f"INV={bc.get('inversion_only',0)} DUP={bc.get('duplication',0)} "
              f"TRANSLOC={bc.get('translocation_only',0)}", file=sys.stderr)


if __name__ == "__main__":
    main()
