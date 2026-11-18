#!/usr/bin/env python3
"""LARGE IS110 element movement detector — uses whole-element boundaries from
records.json (not transposase ORF).

Per ref assembly:
  1. Read records-style TSV (assembly, contig, is_start, is_end, length, ref_id)
     filtered to rows where assembly == --ref-asm.
  2. Extract whole IS body + 5 kb A_flank + 5 kb B_flank.
  3. minimap2 all IS bodies vs species DB (.mmi) in one call.
  4. Cluster nearby PAF blocks, filter cov >= 0.50 / id >= 95%, drop self-hits.
  5. Flank-check each passing cluster (A_flank vs tgt-upstream, B_flank vs
     tgt-downstream — match = >= 50% of 5 kb at >= 90% id).
  6. Verdict: MOVED / SYNTENIC / PARTIAL_MATCH.

Strict layer: cov_IS >= 0.95 AND id_IS >= 99% (essentially identical large element).
"""
import argparse, csv, os, shutil, subprocess, sys, tempfile
from collections import defaultdict
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib_movement as lm  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--ref-asm",  required=True)
    p.add_argument("--ref-fa",   required=True)
    p.add_argument("--species-db", required=True)
    p.add_argument("--species-mmi", default="")
    p.add_argument("--records-tsv", required=True,
                   help="TSV w/ columns ref_id, assembly, contig, is_start, is_end, length")
    p.add_argument("--out", required=True)
    p.add_argument("--flank",           type=int,   default=5000)
    p.add_argument("--min-coverage",    type=float, default=0.50)
    p.add_argument("--min-identity",    type=float, default=95.0)
    p.add_argument("--strict-coverage", type=float, default=0.95)
    p.add_argument("--strict-identity", type=float, default=99.0)
    p.add_argument("--flank-match-cov", type=float, default=0.50)
    p.add_argument("--flank-match-id",  type=float, default=90.0)
    p.add_argument("--cluster-gap",     type=int,   default=2000)
    p.add_argument("--max-secondary", type=int, default=20000,
                   help="minimap2 -N: copies reported per query. Must exceed the "
                        "number of copies in the DB; queries that reach it are "
                        "listed as truncated (default 20000; minimap2's own "
                        "default of 50 saw only the 50 best copies)")
    p.add_argument("--threads", type=int, default=8)
    return p.parse_args()


def load_refs(path, ref_asm):
    out = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["assembly"] != ref_asm: continue
            out.append({
                "ref_id":   r["ref_id"],
                "contig":   r["contig"],
                "is_start": int(r["is_start"]),
                "is_end":   int(r["is_end"]),
                "length":   int(r["length"]),
            })
    return out


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)
    work_root = tempfile.mkdtemp(prefix="largemov_")

    rows = load_refs(args.records_tsv, args.ref_asm)
    if not rows:
        print(f"[{args.ref_asm}] no large IS in records; nothing to do.", flush=True)
        shutil.rmtree(work_root, ignore_errors=True); return

    ref_fa = pysam.FastaFile(args.ref_fa)
    db_fa  = pysam.FastaFile(args.species_db)
    print(f"[{args.ref_asm}] {len(rows)} large IS in ref", flush=True)

    # Build candidates with whole-element body + 5 kb flanks
    candidates = []
    skipped = defaultdict(int)
    for r in rows:
        contig = r["contig"]
        try:
            clen = ref_fa.get_reference_length(contig)
        except (KeyError, ValueError):
            alt = f"{args.ref_asm}|{contig}"
            try:
                clen = ref_fa.get_reference_length(alt); contig = alt
            except (KeyError, ValueError):
                skipped["contig not in --ref-fa"] += 1
                continue
        # is_start/is_end are records.json source_start/source_end: 1-based
        # inclusive. Used as 0-based they dropped the element's first base
        # and put it in the left flank.
        is_s = max(0, r["is_start"] - 1); is_e = min(clen, r["is_end"])
        is_seq = ref_fa.fetch(contig, is_s, is_e).upper()
        lf_seq = ref_fa.fetch(contig, max(0, is_s - args.flank), is_s).upper()
        rf_seq = ref_fa.fetch(contig, is_e, min(clen, is_e + args.flank)).upper()
        if not is_seq:
            skipped["empty element sequence"] += 1
            continue
        cid = f"el{len(candidates)}"
        candidates.append({
            "id": cid, "ref_id": r["ref_id"], "contig": contig,
            "is_start": is_s + 1, "is_end": is_e, "length": r["length"],
            "is_seq": is_seq, "lf": lf_seq, "rf": rf_seq,
        })
    print(f"[{args.ref_asm}] {len(candidates)} large-IS candidates extracted"
          f"{'; SKIPPED ' + str(dict(skipped)) if skipped else ''}", flush=True)
    if not candidates:
        ref_fa.close(); db_fa.close()
        shutil.rmtree(work_root, ignore_errors=True); return

    out_path = os.path.join(args.out, f"{args.ref_asm}_large_movement.tsv")
    try:
        lm.search_and_judge(
            args, candidates, db_fa, work_root, out_path,
            "ref_id\tref_contig\tis_start\tis_end\telement_len",
            lambda c: f"{c['ref_id']}\t{c['contig']}\t{c['is_start']}\t{c['is_end']}\t{c['length']}", "large IS")
    finally:
        ref_fa.close(); db_fa.close()
        shutil.rmtree(work_root, ignore_errors=True)


if __name__ == "__main__":
    main()
