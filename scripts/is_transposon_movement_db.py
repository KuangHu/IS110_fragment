#!/usr/bin/env python3
"""Scalable composite-transposon mobilization detector using the species DB.

For each ref candidate composite transposon B:
  1. minimap2 B against the concatenated species DB (= all assemblies in one FASTA).
  2. Each PAF block's tname is "ASSEMBLY|CONTIG" → assembly identifies the target.
  3. Cluster nearby blocks per (assembly, contig) → one cluster per location.
  4. For each cluster passing thresholds:
        a. Extract tgt's left + right flanks at the cluster position (from the
           indexed species DB).
        b. Compare to ref's A_flank and B_flank → SYNTENIC / MOVED / PARTIAL.

This is ~22,000× faster than aligning B to each target FASTA individually,
because minimap2 indexes the DB once and scans all hits in one pass.
"""
import argparse, csv, os, shutil, sys, tempfile
from collections import defaultdict
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib_movement as lm  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--ref-asm", required=True)
    p.add_argument("--ref-fa",  required=True,
                   help="reference FASTA file (uncompressed or .gz, must have .fai)")
    p.add_argument("--species-db", required=True,
                   help="concatenated species FASTA (uncompressed, indexed with .fai)")
    p.add_argument("--species-mmi", default="",
                   help="optional pre-built minimap2 index (.mmi) for the species DB; "
                        "if not given, uses --species-db directly (slower)")
    p.add_argument("--is-hits",  required=True)
    p.add_argument("--out",      required=True)
    p.add_argument("--min-segment-bp", type=int, default=10000)
    p.add_argument("--max-segment-bp", type=int, default=30000)
    p.add_argument("--flank",          type=int, default=5000)
    p.add_argument("--min-coverage",   type=float, default=0.50)
    p.add_argument("--min-identity",   type=float, default=80.0)
    p.add_argument("--flank-match-cov", type=float, default=0.50)
    p.add_argument("--flank-match-id",  type=float, default=90.0)
    p.add_argument("--cluster-gap",    type=int, default=10000)
    p.add_argument("--require-same-strand", action="store_true")
    p.add_argument("--strict-coverage", type=float, default=0.95)
    p.add_argument("--strict-identity", type=float, default=99.0)
    p.add_argument("--all-pairs", action="store_true",
                   help="Also test non-adjacent IS pairs (A..C over B). Default: "
                        "adjacent pairs only -- with all pairs one relocated region "
                        "gave up to three overlapping MOVED rows")
    p.add_argument("--max-secondary", type=int, default=20000,
                   help="minimap2 -N (copies reported per query; default 20000)")
    p.add_argument("--threads", type=int, default=32)
    return p.parse_args()


def load_is_hits(path):
    idx = defaultdict(lambda: defaultdict(list))
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            asm = r["assembly"]; contig = r["contig"]
            local = contig.split("|", 1)[1] if "|" in contig else contig
            hit = (int(r["tnp_start"]), int(r["tnp_end"]),
                   r.get("tnp_strand", "+"))
            idx[asm][local].append(hit)
    return idx


def merge_is(hits, gap=500):
    if not hits: return []
    sh = sorted(hits, key=lambda x: x[0])
    merged = [sh[0]]
    for h in sh[1:]:
        if h[0] <= merged[-1][1] + gap:
            # disagreeing strands are not "the first one's strand"
            merged[-1] = (min(merged[-1][0], h[0]),
                          max(merged[-1][1], h[1]),
                          merged[-1][2] if merged[-1][2] == h[2] else "?")
        else:
            merged.append(h)
    return merged


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)
    # a rerun that finds nothing must not leave the previous table in place
    _old = os.path.join(args.out, f"{args.ref_asm}_movement.tsv")
    if os.path.exists(_old):
        os.remove(_old)
    is_hits = load_is_hits(args.is_hits)

    work_root = tempfile.mkdtemp(prefix="movedb_")

    ref_fa = pysam.FastaFile(args.ref_fa)
    db_fa  = pysam.FastaFile(args.species_db)
    print(f"[{args.ref_asm}] loaded ref + species DB", flush=True)

    # Build composite-transposon candidates from ref.
    # The query ("body") is the WHOLE unit IS_a..IS_b, with the 5 kb outside
    # it as flanks. It used to be only the segment between the two IS, while
    # the reference flanks were taken outside them: the target "flank" at the
    # segment edge was then mostly the IS copy itself (coverage 0.76 at best
    # for an identical site).
    candidates = []
    skipped = defaultdict(int)
    ref_names = set(ref_fa.references)
    for contig, hits in sorted(is_hits.get(args.ref_asm, {}).items()):
        merged = merge_is(hits)
        # records / is_hits name contigs "assembly|contig"; a per-assembly
        # FASTA has the bare name. Try as given, bare, then prefixed (the
        # retry used to build "asm|asm|contig" and skip every candidate).
        for name in (contig, contig.split("|", 1)[-1],
                     f"{args.ref_asm}|{contig.split('|', 1)[-1]}"):
            if name in ref_names:
                contig = name
                break
        else:
            skipped["contig not in --ref-fa"] += 1
            continue
        clen = ref_fa.get_reference_length(contig)
        for i in range(len(merged)):
            for j in (range(i + 1, len(merged)) if args.all_pairs
                      else range(i + 1, min(i + 2, len(merged)))):
                a, b = merged[i], merged[j]
                # 1-based inclusive ORF coordinates: bases strictly between
                seg_len = b[0] - a[1] - 1
                if seg_len < args.min_segment_bp or seg_len > args.max_segment_bp:
                    continue
                if args.require_same_strand and a[2] != b[2]:
                    continue
                u_s, u_e = a[0] - 1, b[1]            # 0-based half-open unit
                body = ref_fa.fetch(contig, u_s, u_e).upper()
                if not body:
                    skipped["empty sequence"] += 1
                    continue
                candidates.append({
                    "id": f"c{len(candidates)}", "contig": contig,
                    "is_a": a, "is_b": b, "seg_len": seg_len, "is_seq": body,
                    "lf": ref_fa.fetch(contig, max(0, u_s - args.flank), u_s).upper(),
                    "rf": ref_fa.fetch(contig, u_e, min(clen, u_e + args.flank)).upper(),
                })
    print(f"[{args.ref_asm}] {len(candidates)} composite-transposon candidates"
          f"{'; SKIPPED ' + str(dict(skipped)) if skipped else ''}", flush=True)
    out_path = os.path.join(args.out, f"{args.ref_asm}_movement.tsv")
    try:
        if candidates:
            lm.search_and_judge(
                args, candidates, db_fa, work_root, out_path,
                "ref_contig\tis_a_start\tis_b_start\tsegment_bp",
                lambda c: f"{c['contig']}\t{c['is_a'][0]}\t{c['is_b'][0]}\t{c['seg_len']}",
                "composite unit")
    finally:
        ref_fa.close(); db_fa.close()
        shutil.rmtree(work_root, ignore_errors=True)


if __name__ == "__main__":
    main()
