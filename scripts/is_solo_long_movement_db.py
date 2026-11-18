#!/usr/bin/env python3
"""Solo long-IS movement detector across a species DB.

For each long full-length IS in ref (tnp_len >= min, both Pfam domains present):
  1. Extract the IS body + 5 kb upstream (A_flank) + 5 kb downstream (B_flank).
  2. Bulk-align all IS bodies vs species DB index in one minimap2 call.
  3. Cluster nearby forward blocks per (qname, tname), filter by cov / identity.
  4. Drop self-hits (cluster tgt assembly == ref assembly).
  5. For each passing cluster, do the flank check (A_flank vs tgt-upstream,
     B_flank vs tgt-downstream).
  6. Verdict:
        SYNTENIC      both flanks match  (= same context, not mobilized)
        MOVED         neither matches    (= IS body in a different neighborhood)
        PARTIAL_MATCH exactly one matches

Strict layer: cov_IS >= 0.95 AND identity_IS >= 99% (essentially identical IS).
              Loose layer: cov_IS >= 0.50 AND identity_IS >= 95% (per IS110 thresh).
"""
import argparse, csv, os, shutil, subprocess, sys, tempfile
from collections import defaultdict
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib_movement as lm  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--ref-asm", required=True)
    p.add_argument("--ref-fa",  required=True,
                   help="reference FASTA (indexed; same DB the ref lives in is fine)")
    p.add_argument("--species-db", required=True,
                   help="species concatenated FASTA, indexed with .fai")
    p.add_argument("--species-mmi", default="",
                   help="optional pre-built .mmi for the species DB")
    p.add_argument("--is-hits",  required=True)
    p.add_argument("--out",      required=True)
    p.add_argument("--min-tnp-len", type=int, default=1000)
    p.add_argument("--require-both-domains", action=argparse.BooleanOptionalAction,
                   default=True, help="(--no-require-both-domains to turn off; "
                                      "the old store_true flag could not be)")
    p.add_argument("--flank",            type=int,   default=5000)
    p.add_argument("--is-pad",           type=int,   default=50,
                   help="extra bp on each side of tnp ORF when extracting IS body")
    p.add_argument("--min-coverage",     type=float, default=0.50)
    p.add_argument("--min-identity",     type=float, default=95.0)
    p.add_argument("--strict-coverage",  type=float, default=0.95)
    p.add_argument("--strict-identity",  type=float, default=99.0)
    p.add_argument("--flank-match-cov",  type=float, default=0.50)
    p.add_argument("--flank-match-id",   type=float, default=90.0)
    p.add_argument("--cluster-gap",      type=int,   default=2000)
    p.add_argument("--max-secondary", type=int, default=20000,
                   help="minimap2 -N: copies reported per query. Must exceed the "
                        "number of copies in the DB; queries that reach it are "
                        "listed as truncated (default 20000; minimap2's own "
                        "default of 50 saw only the 50 best copies)")
    p.add_argument("--threads", type=int, default=8)
    return p.parse_args()


def load_long_is(path, ref_asm, min_tnp_len, require_both):
    """Return list of (contig, tnp_start, tnp_end, tnp_strand, tnp_len, is_id)."""
    out = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["assembly"] != ref_asm: continue
            try: tlen = int(r["tnp_len"])
            except ValueError: continue
            if tlen < min_tnp_len: continue
            dom = r.get("domains_hit", "")
            if require_both and not ("PF01548" in dom and "PF02371" in dom):
                continue
            out.append((r["contig"], int(r["tnp_start"]), int(r["tnp_end"]),
                        r.get("tnp_strand", "+"), tlen, r["is_id"]))
    return out


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)
    work_root = tempfile.mkdtemp(prefix="solomov_")

    long_is = load_long_is(args.is_hits, args.ref_asm,
                           args.min_tnp_len, args.require_both_domains)
    if not long_is:
        print(f"[{args.ref_asm}] no long IS found; nothing to do.", flush=True)
        shutil.rmtree(work_root, ignore_errors=True); return

    ref_fa = pysam.FastaFile(args.ref_fa)
    db_fa  = pysam.FastaFile(args.species_db)
    print(f"[{args.ref_asm}] {len(long_is)} long IS in ref", flush=True)

    # Build candidates: IS body + A_flank + B_flank
    candidates = []
    skipped = defaultdict(int)
    for (contig, ts, te, strand, tlen, is_id) in long_is:
        try:
            clen = ref_fa.get_reference_length(contig)
        except (KeyError, ValueError):
            alt = f"{args.ref_asm}|{contig}"
            try:
                clen = ref_fa.get_reference_length(alt); contig = alt
            except (KeyError, ValueError):
                skipped["contig not in --ref-fa"] += 1
                continue
        # tnp_start/tnp_end are 1-based inclusive (is_hits.tsv)
        is_s = max(0, ts - 1 - args.is_pad); is_e = min(clen, te + args.is_pad)
        is_seq = ref_fa.fetch(contig, is_s, is_e).upper()
        # flanks start where the padded body ends, so body + flanks are
        # contiguous and the projection in lib_movement lands on the flank
        lf_seq = ref_fa.fetch(contig, max(0, is_s - args.flank), is_s).upper()
        rf_seq = ref_fa.fetch(contig, is_e, min(clen, is_e + args.flank)).upper()
        if not is_seq:
            skipped["empty element sequence"] += 1
            continue
        cid = f"is{len(candidates)}"
        candidates.append({
            "id": cid, "is_id": is_id, "contig": contig,
            "is_start": ts, "is_end": te, "is_strand": strand, "tnp_len": tlen,
            "is_seq": is_seq, "lf": lf_seq, "rf": rf_seq,
        })
    print(f"[{args.ref_asm}] {len(candidates)} solo-IS candidates extracted"
          f"{'; SKIPPED ' + str(dict(skipped)) if skipped else ''}", flush=True)
    if not candidates:
        ref_fa.close(); db_fa.close()
        shutil.rmtree(work_root, ignore_errors=True); return

    out_path = os.path.join(args.out, f"{args.ref_asm}_solo_movement.tsv")
    try:
        lm.search_and_judge(
            args, candidates, db_fa, work_root, out_path,
            "ref_is_id\tref_contig\tis_start\tis_end\ttnp_len",
            lambda c: f"{c['is_id']}\t{c['contig']}\t{c['is_start']}\t{c['is_end']}\t{c['tnp_len']}", "solo IS")
    finally:
        ref_fa.close(); db_fa.close()
        shutil.rmtree(work_root, ignore_errors=True)


if __name__ == "__main__":
    main()
