#!/usr/bin/env python3
"""Split IS110 transposase CDS: two parts of one CDS on the same contig,
same strand, in CDS order, meeting on the CDS, separated on the contig.

This is the one implementation behind hmm_search_is110.sh (step 7),
run_cds_minimap2.sh and search_split_vs_gtdb.sh. Each of them used to carry
its own copy of the pairing loop; those copies accepted an inverted second
half as "split", and kept one row per query (or per query x contig), so a
second locus was lost. The pairing, the junction test and the per-locus
de-duplication are those of detect_is110_variants.py.

Usage:
    detect_split_all.py <in.paf> <out.tsv> [--min-gap 50] [--max-gap 20000]
                        [--min-identity 80] [--min-cov 50]

Output: one row per physical locus (PAF coordinates, 0-based half-open).
"""
import argparse
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import detect_is110_variants as dv  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("paf")
    p.add_argument("out_tsv")
    p.add_argument("--min-gap", type=int, default=50, help="min inserted length on the contig")
    p.add_argument("--max-gap", type=int, default=20000)
    p.add_argument("--min-identity", type=float, default=80)
    p.add_argument("--min-cov", type=float, default=50, help="min %% of the CDS in the two parts")
    p.add_argument("--max-query-gap", type=int, default=50,
                   help="the two parts must meet on the CDS within this many bp")
    p.add_argument("--top", type=int, default=50, help="rows to print")
    a = p.parse_args()

    args = argparse.Namespace(
        min_query_cov=a.min_cov, ins_min_target_gap=a.min_gap, ins_max_target_gap=a.max_gap,
        ins_max_query_gap=a.max_query_gap, del_min_query_gap=10 ** 9, del_max_target_gap=0,
        max_junction_overlap=100, locus_tol=20)
    blocks = dv.parse_paf(a.paf, a.min_identity, a.min_gap, None)
    events, _ = dv.detect_events(blocks, args)
    events = [r for r in dv.deduplicate(events, args.locus_tol) if r["type"] == "insertion"]
    events.sort(key=lambda r: (-r["cov_pct"], r["qname"], r["tname"], r["hi"]["ts"]))

    print(f"\nSplit CDS loci (gap {a.min_gap}-{a.max_gap} bp): {len(events)}")
    print(f"  queries: {len(set(r['qname'] for r in events))}, "
          f"contigs: {len(set(r['tname'] for r in events))}\n")
    for r in events[:a.top]:
        hi, hj = r["hi"], r["hj"]
        print(f"  {r['qname']} (qlen={r['qlen']}) on {r['tname']} ({r['tlen']:,} bp)")
        print(f"    Part1: Q[{hi['qs']}-{hi['qe']}] -> T[{hi['ts']:,}-{hi['te']:,}] {hi['strand']} {hi['ident']:.0f}%")
        print(f"    Part2: Q[{hj['qs']}-{hj['qe']}] -> T[{hj['ts']:,}-{hj['te']:,}] {hj['strand']} {hj['ident']:.0f}%")
        print(f"    Gap: {r['target_gap']:,} bp | Coverage: {r['cov_pct']:.0f}%")

    with open(a.out_tsv, "w") as f:
        f.write("query\tqlen\tsplit_type\tcontig1\ttlen1\tpart1_qs\tpart1_qe\tpart1_ts\tpart1_te\t"
                "part1_strand\tpart1_ident\tcontig2\ttlen2\tpart2_qs\tpart2_qe\tpart2_ts\tpart2_te\t"
                "part2_strand\tpart2_ident\tgap\tcds_coverage\n")
        for r in events:
            hi, hj = r["hi"], r["hj"]
            f.write(f"{r['qname']}\t{r['qlen']}\tsame_contig\t"
                    f"{hi['tname']}\t{hi['tlen']}\t{hi['qs']}\t{hi['qe']}\t{hi['ts']}\t{hi['te']}\t{hi['strand']}\t{hi['ident']:.1f}\t"
                    f"{hj['tname']}\t{hj['tlen']}\t{hj['qs']}\t{hj['qe']}\t{hj['ts']}\t{hj['te']}\t{hj['strand']}\t{hj['ident']:.1f}\t"
                    f"{r['target_gap']}\t{r['cov_pct']:.1f}\n")
    print(f"\nResults written to {a.out_tsv}")
    print("\n=== Split loci per query (top 20) ===")
    for q, c in Counter(r["qname"] for r in events).most_common(20):
        print(f"  {q}: {c}")


if __name__ == "__main__":
    main()
