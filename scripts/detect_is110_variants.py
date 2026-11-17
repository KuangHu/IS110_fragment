#!/usr/bin/env python3
"""Detect IS110 structural variants from minimap2 PAF alignments.

Detects two scenarios:

Scenario 1 — INSERTION in genome (genome longer than reference):
    Reference:  [========= IS110 CDS =========]
    Genome:     [IS110 part1]---gap---[IS110 part2]
    Signal: two hits cover different query portions, large GAP on TARGET,
            small or no gap on query.

Scenario 2 — DELETION in genome (genome shorter than reference):
    Reference:  [IS110 part1]---DNA---[IS110 part2]
    Genome:     [IS110 part1][IS110 part2]
    Signal: two hits cover different query portions, large GAP on QUERY,
            small or no gap on target (parts are adjacent on genome).

Usage:
    python detect_is110_variants.py <input.paf> <output_prefix> [options]
"""

import argparse
import os
import sys
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("paf", help="Input PAF file from minimap2")
    p.add_argument("output_prefix", help="Output prefix for result files")
    p.add_argument("--min-identity", type=float, default=95,
                   help="Min alignment identity %% (default: 95)")
    p.add_argument("--min-query-cov", type=float, default=50,
                   help="Min combined query coverage %% (default: 50)")
    # Scenario 1: insertion
    p.add_argument("--ins-min-target-gap", type=int, default=500,
                   help="Scenario 1: min gap on target (default: 500)")
    p.add_argument("--ins-max-target-gap", type=int, default=50000,
                   help="Scenario 1: max gap on target (default: 50000)")
    # Scenario 2: deletion
    p.add_argument("--del-min-query-gap", type=int, default=200,
                   help="Scenario 2: min gap on query (default: 200)")
    p.add_argument("--del-max-target-gap", type=int, default=200,
                   help="Scenario 2: max distance between parts on target (default: 200)")
    return p.parse_args()


def parse_paf(path, min_identity):
    """Parse PAF, filter by identity, group by query."""
    hits_by_query = defaultdict(list)
    total = 0
    kept = 0
    with open(path) as fh:
        for line in fh:
            c = line.strip().split("\t")
            if len(c) < 12:
                continue
            total += 1
            matches, block = int(c[9]), int(c[10])
            ident = matches / block * 100 if block > 0 else 0
            if ident < min_identity:
                continue
            kept += 1
            hits_by_query[c[0]].append({
                "qname": c[0],
                "qlen": int(c[1]),
                "qs": int(c[2]),
                "qe": int(c[3]),
                "strand": c[4],
                "tname": c[5],
                "tlen": int(c[6]),
                "ts": int(c[7]),
                "te": int(c[8]),
                "ident": ident,
            })
    print(f"Parsed {total} alignments, {kept} pass identity >= {min_identity}%")
    print(f"Unique queries with hits: {len(hits_by_query)}")
    return hits_by_query


def query_overlap(hi, hj):
    """Fraction of overlap on query between two hits."""
    ov = max(0, min(hi["qe"], hj["qe"]) - max(hi["qs"], hj["qs"]))
    min_span = min(hi["qe"] - hi["qs"], hj["qe"] - hj["qs"])
    if min_span == 0:
        return 1.0
    return ov / min_span


def query_coverage(hi, hj):
    """Combined query coverage percentage."""
    covered = set(range(hi["qs"], hi["qe"])) | set(range(hj["qs"], hj["qe"]))
    return len(covered) / hi["qlen"] * 100


def query_gap(hi, hj):
    """Gap on query between two hits (sorted by query start)."""
    a, b = sorted([hi, hj], key=lambda x: x["qs"])
    return b["qs"] - a["qe"]


def detect_insertions(hits_by_query, args):
    """Scenario 1: genome has insertion splitting the IS110.

    Two query portions map to distant spots on the SAME target contig.
    Large gap on target, query portions are complementary.
    """
    results = []
    for qname, hlist in hits_by_query.items():
        # Group by target contig
        by_contig = defaultdict(list)
        for h in hlist:
            by_contig[h["tname"]].append(h)

        for tname, ths in by_contig.items():
            if len(ths) < 2:
                continue
            ths.sort(key=lambda x: x["ts"])
            for i in range(len(ths)):
                for j in range(i + 1, len(ths)):
                    hi, hj = ths[i], ths[j]
                    tgap = hj["ts"] - hi["te"]
                    if tgap < args.ins_min_target_gap or tgap > args.ins_max_target_gap:
                        continue
                    if query_overlap(hi, hj) > 0.5:
                        continue
                    cov = query_coverage(hi, hj)
                    if cov < args.min_query_cov:
                        continue
                    qgap = query_gap(hi, hj)
                    results.append({
                        "type": "insertion",
                        "qname": qname,
                        "qlen": hi["qlen"],
                        "tname": tname,
                        "tlen": hi["tlen"],
                        "hi": hi,
                        "hj": hj,
                        "target_gap": tgap,
                        "query_gap": qgap,
                        "cov_pct": cov,
                    })
    return results


def detect_deletions(hits_by_query, args):
    """Scenario 2: genome has deletion — IS110 is shorter than reference.

    Two query portions (far apart on query) map to ADJACENT spots on target.
    Large gap on query, small gap on target.
    """
    results = []
    for qname, hlist in hits_by_query.items():
        # Group by target contig
        by_contig = defaultdict(list)
        for h in hlist:
            by_contig[h["tname"]].append(h)

        for tname, ths in by_contig.items():
            if len(ths) < 2:
                continue
            # Sort by query start to find query gaps
            ths_qsort = sorted(ths, key=lambda x: x["qs"])
            for i in range(len(ths_qsort)):
                for j in range(i + 1, len(ths_qsort)):
                    hi, hj = ths_qsort[i], ths_qsort[j]
                    # Must cover different query portions
                    if query_overlap(hi, hj) > 0.5:
                        continue
                    # Large gap on query
                    qgap = hj["qs"] - hi["qe"]
                    if qgap < args.del_min_query_gap:
                        continue
                    # Small gap on target (parts are adjacent)
                    tgap = abs(min(hj["ts"], hj["te"]) - max(hi["ts"], hi["te"]))
                    # Also check the other ordering on target
                    tgap_alt = abs(min(hi["ts"], hi["te"]) - max(hj["ts"], hj["te"]))
                    tgap = min(tgap, tgap_alt)
                    if tgap > args.del_max_target_gap:
                        continue
                    cov = query_coverage(hi, hj)
                    if cov < args.min_query_cov:
                        continue
                    results.append({
                        "type": "deletion",
                        "qname": qname,
                        "qlen": hi["qlen"],
                        "tname": tname,
                        "tlen": hi["tlen"],
                        "hi": hi,
                        "hj": hj,
                        "target_gap": tgap,
                        "query_gap": qgap,
                        "cov_pct": cov,
                    })
    return results


def deduplicate(results):
    """Keep best hit per (query, target_contig, type)."""
    best = {}
    for r in results:
        key = (r["type"], r["qname"], r["tname"])
        if key not in best or r["cov_pct"] > best[key]["cov_pct"]:
            best[key] = r
    return sorted(best.values(), key=lambda x: (-x["cov_pct"], x["qname"]))


def write_tsv(results, path):
    """Write results to TSV."""
    with open(path, "w") as f:
        f.write(
            "type\tquery\tqlen\ttarget_contig\ttlen\t"
            "part1_qs\tpart1_qe\tpart1_ts\tpart1_te\tpart1_strand\tpart1_ident\t"
            "part2_qs\tpart2_qe\tpart2_ts\tpart2_te\tpart2_strand\tpart2_ident\t"
            "target_gap\tquery_gap\tcds_coverage\n"
        )
        for r in results:
            hi, hj = r["hi"], r["hj"]
            f.write(
                f"{r['type']}\t{r['qname']}\t{r['qlen']}\t{r['tname']}\t{r['tlen']}\t"
                f"{hi['qs']}\t{hi['qe']}\t{hi['ts']}\t{hi['te']}\t{hi['strand']}\t{hi['ident']:.1f}\t"
                f"{hj['qs']}\t{hj['qe']}\t{hj['ts']}\t{hj['te']}\t{hj['strand']}\t{hj['ident']:.1f}\t"
                f"{r['target_gap']}\t{r['query_gap']}\t{r['cov_pct']:.1f}\n"
            )


def print_summary(results, label, n=20):
    """Print top results."""
    subset = [r for r in results if r["type"] == label]
    if not subset:
        print(f"  No {label} cases found.")
        return

    # Unique queries and targets
    uq = len(set(r["qname"] for r in subset))
    ut = len(set(r["tname"] for r in subset))
    print(f"  {len(subset)} cases ({uq} unique queries, {ut} unique targets)")

    for r in subset[:n]:
        hi, hj = r["hi"], r["hj"]
        print(f"\n  {r['qname']} (qlen={r['qlen']})")
        print(f"    Target: {r['tname']} ({r['tlen']:,} bp)")
        print(f"    Part1: Q[{hi['qs']}-{hi['qe']}] -> T[{hi['ts']:,}-{hi['te']:,}] {hi['strand']} {hi['ident']:.0f}%")
        print(f"    Part2: Q[{hj['qs']}-{hj['qe']}] -> T[{hj['ts']:,}-{hj['te']:,}] {hj['strand']} {hj['ident']:.0f}%")
        if label == "insertion":
            print(f"    Target gap: {r['target_gap']:,} bp | Query coverage: {r['cov_pct']:.0f}%")
        else:
            print(f"    Query gap: {r['query_gap']:,} bp (deleted) | Target gap: {r['target_gap']} bp | Query coverage: {r['cov_pct']:.0f}%")


def main():
    args = parse_args()

    print("=" * 60)
    print("IS110 Structural Variant Detection")
    print("=" * 60)
    print(f"Input:  {args.paf}")
    print(f"Output: {args.output_prefix}_*.tsv")
    print()

    hits_by_query = parse_paf(args.paf, args.min_identity)
    print()

    # Scenario 1: insertions
    print("--- Scenario 1: INSERTION in genome (target gap) ---")
    print(f"    Target gap range: {args.ins_min_target_gap}-{args.ins_max_target_gap} bp")
    insertions = detect_insertions(hits_by_query, args)
    insertions = deduplicate(insertions)
    print_summary(insertions, "insertion")
    print()

    # Scenario 2: deletions
    print("--- Scenario 2: DELETION in genome (query gap) ---")
    print(f"    Query gap min: {args.del_min_query_gap} bp")
    print(f"    Target gap max: {args.del_max_target_gap} bp")
    deletions = detect_deletions(hits_by_query, args)
    deletions = deduplicate(deletions)
    print_summary(deletions, "deletion")
    print()

    # Write outputs
    all_results = insertions + deletions
    if all_results:
        combined_path = f"{args.output_prefix}_all.tsv"
        write_tsv(all_results, combined_path)
        print(f"Combined results: {combined_path} ({len(all_results)} rows)")

    if insertions:
        ins_path = f"{args.output_prefix}_insertions.tsv"
        write_tsv(insertions, ins_path)
        print(f"Insertions:       {ins_path} ({len(insertions)} rows)")

    if deletions:
        del_path = f"{args.output_prefix}_deletions.tsv"
        write_tsv(deletions, del_path)
        print(f"Deletions:        {del_path} ({len(deletions)} rows)")

    print()
    print(f"Total: {len(insertions)} insertions + {len(deletions)} deletions = {len(all_results)} variants")


if __name__ == "__main__":
    main()
