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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paf_utils import aligned_identity, large_indels  # noqa: E402


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
    p.add_argument("--ins-max-query-gap", type=int, default=50,
                   help="Scenario 1: the two halves must meet on the IS within "
                        "this many bp (default: 50)")
    p.add_argument("--ins-max-target-gap", type=int, default=50000,
                   help="Scenario 1: max gap on target (default: 50000)")
    # Scenario 2: deletion
    p.add_argument("--del-min-query-gap", type=int, default=200,
                   help="Scenario 2: min gap on query (default: 200)")
    p.add_argument("--del-max-target-gap", type=int, default=200,
                   help="Scenario 2: max distance between parts on target (default: 200)")
    return p.parse_args()


def parse_paf(path, min_identity, ins_min=500, del_min=200, single_row=None):
    """Parse PAF, filter by identity, group by query."""
    if single_row is None:
        single_row = []
    hits_by_query = defaultdict(list)
    total = 0
    kept = 0
    with open(path) as fh:
        for line in fh:
            c = line.strip().split("\t")
            if len(c) < 12:
                continue
            total += 1
            ident = aligned_identity(c)        # gaps excluded (see paf_utils)
            if ident < min_identity:
                continue
            kept += 1
            # an indel bridged INSIDE this alignment row is an event too
            for kind, n, qp, tp in large_indels(c, min(ins_min, del_min)):
                single_row.append((kind, n, qp, tp, c, ident))
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
                    # ONE interrupted copy: both halves on the same strand, in
                    # the same order on the IS as on the contig, and meeting
                    # on the IS (nothing of the IS missing or doubled at the
                    # junction). Without these tests an inverted second half,
                    # two halves in swapped order, or two partial copies 30 kb
                    # apart were all "1 insertion".
                    if hi["strand"] != hj["strand"]:
                        continue
                    first, second = (hi, hj) if hi["strand"] == "+" else (hj, hi)
                    if not first["qs"] < second["qs"]:
                        continue
                    if abs(second["qs"] - first["qe"]) > args.ins_max_query_gap:
                        continue
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
                    # same strand, and in the same order on the contig as on
                    # the IS (hi precedes hj on the query)
                    if hi["strand"] != hj["strand"]:
                        continue
                    if (hi["ts"] < hj["ts"]) != (hi["strand"] == "+"):
                        continue
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


def row_event(ev):
    """An indel inside ONE alignment row, as a two-part result like those of
    detect_insertions / detect_deletions (hi = target-left part)."""
    kind, n, qp, tp, c, ident = ev
    qlen, qs, qe, strand = int(c[1]), int(c[2]), int(c[3]), c[4]
    tname, tlen, ts, te = c[5], int(c[6]), int(c[7]), int(c[8])
    base = {"qname": c[0], "qlen": qlen, "strand": strand, "tname": tname,
            "tlen": tlen, "ident": ident}
    if kind == "D":                       # extra bases in the genome
        t_left, t_right = (ts, tp), (tp + n, te)
        q_first, q_second = (qs, qp), (qp, qe)
    else:                                 # bases of the IS missing in the genome
        t_left, t_right = (ts, tp), (tp, te)
        q_first, q_second = (qs, qp), (qp + n, qe)
    q_left, q_right = (q_first, q_second) if strand == "+" else (q_second, q_first)
    hi = dict(base, qs=q_left[0], qe=q_left[1], ts=t_left[0], te=t_left[1])
    hj = dict(base, qs=q_right[0], qe=q_right[1], ts=t_right[0], te=t_right[1])
    cov = ((q_first[1] - q_first[0]) + (q_second[1] - q_second[0])) / qlen * 100
    return {"type": "insertion" if kind == "D" else "deletion", "qname": c[0], "qlen": qlen,
            "tname": tname, "tlen": tlen, "hi": hi, "hj": hj,
            "target_gap": n if kind == "D" else 0,
            "query_gap": 0 if kind == "D" else n, "cov_pct": cov}


def deduplicate(results, tol=20):
    """One row per LOCUS: (type, contig, junction positions within `tol` bp),
    keeping the best-covering reference.

    The key used to be (type, query, contig): the same interrupted copy was
    reported once for every near-identical reference IS (and its insert then
    written once per reference to the fragment FASTA, inflating cluster sizes),
    while a SECOND event on the same contig was dropped."""
    def junction(r):
        a, b = sorted((r["hi"], r["hj"]), key=lambda h: h["ts"])
        return a["te"], b["ts"]
    kept = []
    for r in sorted(results, key=lambda x: (-x["cov_pct"], x["qname"], x["tname"])):
        j = junction(r)
        if any(k["type"] == r["type"] and k["tname"] == r["tname"]
               and abs(junction(k)[0] - j[0]) <= tol and abs(junction(k)[1] - j[1]) <= tol
               for k in kept):
            continue
        kept.append(r)
    return kept


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

    single_row = []
    hits_by_query = parse_paf(args.paf, args.min_identity,
                              args.ins_min_target_gap, args.del_min_query_gap, single_row)
    print()

    # Scenario 1: insertions
    print("--- Scenario 1: INSERTION in genome (target gap) ---")
    print(f"    Target gap range: {args.ins_min_target_gap}-{args.ins_max_target_gap} bp")
    insertions = detect_insertions(hits_by_query, args)
    insertions += [r for r in map(row_event, single_row)
                   if r["type"] == "insertion"
                   and args.ins_min_target_gap <= r["target_gap"] <= args.ins_max_target_gap
                   and r["cov_pct"] >= args.min_query_cov]
    insertions = deduplicate(insertions)
    print_summary(insertions, "insertion")
    print()

    # Scenario 2: deletions
    print("--- Scenario 2: DELETION in genome (query gap) ---")
    print(f"    Query gap min: {args.del_min_query_gap} bp")
    print(f"    Target gap max: {args.del_max_target_gap} bp")
    deletions = detect_deletions(hits_by_query, args)
    deletions += [r for r in map(row_event, single_row)
                  if r["type"] == "deletion" and r["query_gap"] >= args.del_min_query_gap
                  and r["cov_pct"] >= args.min_query_cov]
    deletions = deduplicate(deletions)
    print_summary(deletions, "deletion")
    print()

    # Write outputs
    # Always written, header-only when empty: a file left by an earlier run
    # was otherwise read by the next step as this run's result.
    all_results = insertions + deletions
    for label, rows, suffix in (("Combined results:", all_results, "all"),
                                ("Insertions:      ", insertions, "insertions"),
                                ("Deletions:       ", deletions, "deletions")):
        out_path = f"{args.output_prefix}_{suffix}.tsv"
        write_tsv(rows, out_path)
        print(f"{label} {out_path} ({len(rows)} rows)")

    print()
    print(f"Total: {len(insertions)} insertions + {len(deletions)} deletions = {len(all_results)} variants")


if __name__ == "__main__":
    main()
