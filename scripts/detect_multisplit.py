#!/usr/bin/env python3
"""Detect multi-fragment IS110 patterns from minimap2 PAF alignments.

For each (query, target_contig) pair, find all high-identity alignments and
classify how many times the reference IS110 is split.

A single split event can be either:
  - INSERTION: non-overlapping query portions separated by a target gap
  - DELETION:  non-overlapping query portions adjacent on target (query gap only)

We count the number of split events:
  - 0 splits: intact copy (single alignment covering >= min_query_cov)
  - 1 split:  one disruption  (2 segments)
  - 2 splits: two disruptions (3 segments)
  - N splits: (N+1) segments

Usage:
    python detect_multisplit.py <input.paf> <output_prefix> [options]
"""

import argparse
import os
import sys
import csv
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paf_utils import aligned_identity  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("paf", help="Input PAF file from minimap2")
    p.add_argument("output_prefix", help="Output prefix for result files")
    p.add_argument("--min-identity", type=float, default=95,
                   help="Min alignment identity %% (default: 95)")
    p.add_argument("--min-aln-len", type=int, default=200,
                   help="Min alignment length on query (default: 200)")
    p.add_argument("--min-query-cov", type=float, default=50,
                   help="Min combined query coverage %% (default: 50)")
    p.add_argument("--min-qgap-deletion", type=int, default=200,
                   help="Min query gap to call a deletion (default: 200)")
    p.add_argument("--max-tgap-deletion", type=int, default=200,
                   help="Max target gap for deletion (adjacency) (default: 200)")
    p.add_argument("--min-tgap-insertion", type=int, default=500,
                   help="Min target gap to call an insertion (default: 500)")
    p.add_argument("--max-tgap-insertion", type=int, default=50000,
                   help="Max target gap for insertion (default: 50000)")
    p.add_argument("--max-query-overlap", type=float, default=0.3,
                   help="Max fraction of query overlap between consecutive segments (default: 0.3)")
    return p.parse_args()


def parse_paf(path, min_identity, min_aln_len):
    """Parse PAF, filter, group by (query, target)."""
    by_qt = defaultdict(list)
    total = 0
    kept = 0
    with open(path) as fh:
        for line in fh:
            c = line.strip().split("\t")
            if len(c) < 12:
                continue
            total += 1
            ident = aligned_identity(c)   # gap columns excluded (paf_utils)
            if ident < min_identity:
                continue
            qs, qe = int(c[2]), int(c[3])
            if qe - qs < min_aln_len:
                continue
            kept += 1
            key = (c[0], c[5])
            by_qt[key].append({
                "qname": c[0],
                "qlen": int(c[1]),
                "qs": qs,
                "qe": qe,
                "strand": c[4],
                "tname": c[5],
                "tlen": int(c[6]),
                "ts": int(c[7]),
                "te": int(c[8]),
                "ident": ident,
            })
    print(f"Parsed {total:,} alignments, {kept:,} pass identity >= {min_identity}% and len >= {min_aln_len}")
    print(f"Unique (query, target) pairs: {len(by_qt):,}")
    return by_qt


def select_colinear_chain(hits, args):
    """Pick the best co-linear chain of non-overlapping query segments.

    Uses greedy approach grouped by strand: keep all hits on the dominant strand,
    sort by query start, remove overlapping hits (keep larger one).
    """
    if not hits:
        return []

    # Group by strand, pick dominant
    fwd = [h for h in hits if h["strand"] == "+"]
    rev = [h for h in hits if h["strand"] == "-"]
    fwd_qcov = sum(h["qe"] - h["qs"] for h in fwd)
    rev_qcov = sum(h["qe"] - h["qs"] for h in rev)
    chosen = fwd if fwd_qcov >= rev_qcov else rev

    # Sort by query start
    chosen = sorted(chosen, key=lambda x: x["qs"])

    # Greedy: remove overlapping hits (keep larger span)
    result = []
    for h in chosen:
        if not result:
            result.append(h)
            continue
        prev = result[-1]
        ov = max(0, min(h["qe"], prev["qe"]) - max(h["qs"], prev["qs"]))
        min_span = min(h["qe"] - h["qs"], prev["qe"] - prev["qs"])
        if min_span > 0 and ov / min_span > args.max_query_overlap:
            # Overlap on query; keep the one with larger span
            if (h["qe"] - h["qs"]) > (prev["qe"] - prev["qs"]):
                result[-1] = h
        else:
            result.append(h)
    return result


def classify_splits(chain, args):
    """From an ordered chain of non-overlapping segments, classify each gap."""
    if len(chain) < 2:
        return []

    # On target, determine if direction is increasing (same strand + co-linear)
    events = []
    for i in range(len(chain) - 1):
        a, b = chain[i], chain[i + 1]
        qgap = b["qs"] - a["qe"]
        # Target gap depends on strand direction
        if a["strand"] == "+" == b["strand"]:
            tgap = b["ts"] - a["te"]
        elif a["strand"] == "-" == b["strand"]:
            tgap = a["ts"] - b["te"]
        else:
            tgap = None  # strand flip; skip classification
            continue

        # Classify
        if qgap >= args.min_qgap_deletion and abs(tgap) <= args.max_tgap_deletion:
            events.append(("deletion", qgap, tgap))
        elif qgap < args.min_qgap_deletion and args.min_tgap_insertion <= tgap <= args.max_tgap_insertion:
            events.append(("insertion", qgap, tgap))
        elif qgap >= args.min_qgap_deletion and tgap >= args.min_tgap_insertion:
            events.append(("complex", qgap, tgap))  # both query and target gap
        else:
            events.append(("other", qgap, tgap))  # e.g. too-close fragments
    return events


def query_coverage(chain):
    """Combined query coverage percentage using interval merge (memory efficient)."""
    if not chain:
        return 0
    intervals = sorted([(h["qs"], h["qe"]) for h in chain])
    merged_len = 0
    cur_s, cur_e = intervals[0]
    for s, e in intervals[1:]:
        if s <= cur_e:
            cur_e = max(cur_e, e)
        else:
            merged_len += cur_e - cur_s
            cur_s, cur_e = s, e
    merged_len += cur_e - cur_s
    qlen = chain[0]["qlen"]
    return merged_len / qlen * 100 if qlen > 0 else 0


def main():
    args = parse_args()

    print("=" * 70)
    print("IS110 Multi-fragment Detection")
    print("=" * 70)
    print(f"Input:  {args.paf}")
    print(f"Identity threshold: {args.min_identity}%")
    print()

    by_qt = parse_paf(args.paf, args.min_identity, args.min_aln_len)
    print()

    # Process each (query, target) pair
    all_results = []
    for (qname, tname), hits in by_qt.items():
        chain = select_colinear_chain(hits, args)
        if not chain:
            continue
        cov = query_coverage(chain)
        if cov < args.min_query_cov:
            continue

        events = classify_splits(chain, args)
        n_real = sum(1 for t, _, _ in events if t in ("insertion", "deletion", "complex"))
        event_types = [t for t, _, _ in events]

        all_results.append({
            "qname": qname,
            "qlen": chain[0]["qlen"],
            "tname": tname,
            "tlen": chain[0]["tlen"],
            "n_segments": len(chain),
            "n_splits": n_real,
            "event_types": ",".join(event_types) if event_types else "intact",
            "query_coverage": cov,
            "strand": chain[0]["strand"],
            "chain": chain,
            "events": events,
        })

    # ── Summary by split count ──
    from collections import Counter
    split_counts = Counter(r["n_splits"] for r in all_results)
    print("=== Distribution of split events per (query, target_contig) ===")
    for k in sorted(split_counts.keys()):
        n = split_counts[k]
        label = "intact" if k == 0 else f"{k} split{'s' if k>1 else ''}"
        print(f"  {k} splits ({label}): {n:,} ({n/len(all_results)*100:.1f}%)")
    print(f"  Total: {len(all_results):,}")
    print()

    # Breakdown of multi-split cases by event-type pattern
    multi = [r for r in all_results if r["n_splits"] >= 2]
    if multi:
        pattern_counts = Counter(r["event_types"] for r in multi)
        print(f"=== Top event patterns for multi-split cases (n={len(multi):,}) ===")
        for pat, n in pattern_counts.most_common(15):
            print(f"  {pat}: {n:,}")
        print()

    # ── Write outputs ──
    # Full results (one row per chain, events comma-separated)
    out_path = f"{args.output_prefix}_chains.tsv"
    with open(out_path, "w") as f:
        f.write("query\tqlen\ttarget\ttlen\tstrand\tn_segments\tn_splits\t"
                "event_types\tquery_coverage\tsegments\n")
        for r in all_results:
            segs = ";".join(
                f"Q[{h['qs']}-{h['qe']}]->T[{h['ts']}-{h['te']}]({h['ident']:.0f}%)"
                for h in r["chain"]
            )
            f.write(
                f"{r['qname']}\t{r['qlen']}\t{r['tname']}\t{r['tlen']}\t"
                f"{r['strand']}\t{r['n_segments']}\t{r['n_splits']}\t"
                f"{r['event_types']}\t{r['query_coverage']:.1f}\t{segs}\n"
            )
    print(f"Full chain results:   {out_path} ({len(all_results):,} rows)")

    # Multi-split only
    if multi:
        multi_path = f"{args.output_prefix}_multisplit.tsv"
        with open(multi_path, "w") as f:
            f.write("query\tqlen\ttarget\ttlen\tstrand\tn_segments\tn_splits\t"
                    "event_types\tquery_coverage\tsegments\n")
            for r in sorted(multi, key=lambda x: (-x["n_splits"], -x["query_coverage"])):
                segs = ";".join(
                    f"Q[{h['qs']}-{h['qe']}]->T[{h['ts']}-{h['te']}]({h['ident']:.0f}%)"
                    for h in r["chain"]
                )
                f.write(
                    f"{r['qname']}\t{r['qlen']}\t{r['tname']}\t{r['tlen']}\t"
                    f"{r['strand']}\t{r['n_segments']}\t{r['n_splits']}\t"
                    f"{r['event_types']}\t{r['query_coverage']:.1f}\t{segs}\n"
                )
        print(f"Multi-split (>=2):    {multi_path} ({len(multi):,} rows)")

    # Show top examples of triple+ splits
    triple = [r for r in all_results if r["n_splits"] >= 2]
    triple.sort(key=lambda x: (-x["n_splits"], -x["query_coverage"]))
    if triple:
        print(f"\n=== Top multi-split examples ===")
        for r in triple[:10]:
            print(f"\n  {r['qname'][:20]}... on {r['tname']}")
            print(f"    qlen={r['qlen']} tlen={r['tlen']} strand={r['strand']} "
                  f"n_splits={r['n_splits']} cov={r['query_coverage']:.0f}%")
            for h in r["chain"]:
                print(f"      Q[{h['qs']:>6}-{h['qe']:<6}] -> T[{h['ts']:>9,}-{h['te']:<9,}] "
                      f"{h['strand']} {h['ident']:.0f}%")
            for t, qg, tg in r["events"]:
                print(f"      → {t}: qgap={qg} tgap={tg}")


if __name__ == "__main__":
    main()
