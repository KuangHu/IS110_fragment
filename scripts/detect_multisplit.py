#!/usr/bin/env python3
"""Detect multi-fragment IS110 patterns from minimap2 PAF alignments.

For each locus chain (aligned blocks of one reference IS on one contig and
strand, linked by insertion / deletion junctions; see
detect_is110_variants.build_chains) count how many times the IS is split.
Rows are per (query, locus): the same physical copy appears once for every
reference IS that matches it.

A single split event can be either:
  - INSERTION: non-overlapping query portions separated by a target gap
  - DELETION:  non-overlapping query portions adjacent on target (query gap only)

We count the number of split events:
  - 0 splits: intact copy (aligned coverage >= --min-intact-cov), otherwise
              status 'partial' (not counted)
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
import detect_is110_variants as dv  # noqa: E402


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
    p.add_argument("--min-intact-cov", type=float, default=95,
                   help="Min aligned query coverage %% to call a chain without "
                        "events 'intact' (default: 95)")
    p.add_argument("--max-query-overlap", type=float, default=0.3,
                   help="Max fraction of query overlap between consecutive segments (default: 0.3)")
    return p.parse_args()


def locus_rows(paf, args):
    """One row per LOCUS CHAIN (see detect_is110_variants.build_chains).

    The old code built one chain per (query, contig): it kept the hits of the
    "dominant" strand only and one hit per query position, so two copies on a
    contig were one row, a copy on the other strand vanished, halves from two
    different loci became one "100% coverage" row, and a single half-length
    hit was "0 splits = intact". Indels inside one alignment row were not
    seen at all."""
    dargs = argparse.Namespace(
        min_query_cov=0, ins_min_target_gap=args.min_tgap_insertion,
        ins_max_target_gap=args.max_tgap_insertion, ins_max_query_gap=args.min_qgap_deletion - 1,
        del_min_query_gap=args.min_qgap_deletion, del_max_target_gap=args.max_tgap_deletion,
        max_junction_overlap=100)
    blocks = dv.parse_paf(paf, args.min_identity, args.min_tgap_insertion, args.min_qgap_deletion)
    blocks = {k: [b for b in v if b["qe"] - b["qs"] >= args.min_aln_len]
              for k, v in blocks.items()}
    chains, _ = dv.build_chains({k: v for k, v in blocks.items() if v}, dargs)
    rows = []
    for ch in chains:
        cov = ch["cov_pct"]
        events = [(t, qg, tg) for _, _, (t, tg, qg) in ch["links"]]
        if events:
            status = "variant"
        elif cov >= args.min_intact_cov:
            status = "intact"
        else:
            status = "partial"          # a piece of the IS, nothing known about the rest
        if status != "partial" and cov < args.min_query_cov:
            status = "partial"
        b0 = ch["blocks"][0]
        rows.append({
            "qname": b0["qname"], "qlen": b0["qlen"], "tname": b0["tname"], "tlen": b0["tlen"],
            "n_segments": len(ch["blocks"]), "n_splits": len(events),
            "event_types": ",".join(t for t, _, _ in events) if events else status,
            "query_coverage": cov, "strand": b0["strand"], "status": status,
            "chain": sorted(ch["blocks"], key=lambda h: h["qs"]), "events": events,
        })
    return rows


def main():
    args = parse_args()

    print("=" * 70)
    print("IS110 Multi-fragment Detection")
    print("=" * 70)
    print(f"Input:  {args.paf}")
    print(f"Identity threshold: {args.min_identity}%")
    print()

    every = locus_rows(args.paf, args)
    n_partial = sum(1 for r in every if r["status"] == "partial")
    all_results = [r for r in every if r["status"] != "partial"]
    print(f"Locus chains: {len(every):,}  (intact or variant: {len(all_results):,}; "
          f"partial, not counted below: {n_partial:,})")
    print()
    if not all_results:
        print("No intact or variant locus.")

    # ── Summary by split count ──
    from collections import Counter
    split_counts = Counter(r["n_splits"] for r in all_results)
    print("=== Distribution of split events per (query, locus) ===")
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
    # Full results (one row per chain, events comma-separated); partial
    # chains are written too, with status=partial
    out_path = f"{args.output_prefix}_chains.tsv"
    all_written = all_results + [r for r in every if r["status"] == "partial"]
    with open(out_path, "w") as f:
        f.write("query\tqlen\ttarget\ttlen\tstrand\tn_segments\tn_splits\t"
                "event_types\tquery_coverage\tsegments\tstatus\n")
        for r in all_written:
            segs = ";".join(
                f"Q[{h['qs']}-{h['qe']}]->T[{h['ts']}-{h['te']}]({h['ident']:.0f}%)"
                for h in r["chain"]
            )
            f.write(
                f"{r['qname']}\t{r['qlen']}\t{r['tname']}\t{r['tlen']}\t"
                f"{r['strand']}\t{r['n_segments']}\t{r['n_splits']}\t"
                f"{r['event_types']}\t{r['query_coverage']:.1f}\t{segs}\t{r['status']}\n"
            )
    print(f"Full chain results:   {out_path} ({len(all_results):,} rows)")

    # Multi-split only
    # always written (header only when empty): an old file must not pass for this run's
    multi_path = f"{args.output_prefix}_multisplit.tsv"
    with open(multi_path, "w") as f:
        f.write("query\tqlen\ttarget\ttlen\tstrand\tn_segments\tn_splits\t"
                "event_types\tquery_coverage\tsegments\tstatus\n")
        for r in sorted(multi, key=lambda x: (-x["n_splits"], -x["query_coverage"])):
            segs = ";".join(
                f"Q[{h['qs']}-{h['qe']}]->T[{h['ts']}-{h['te']}]({h['ident']:.0f}%)"
                for h in r["chain"]
            )
            f.write(
                f"{r['qname']}\t{r['qlen']}\t{r['tname']}\t{r['tlen']}\t"
                f"{r['strand']}\t{r['n_segments']}\t{r['n_splits']}\t"
                f"{r['event_types']}\t{r['query_coverage']:.1f}\t{segs}\t{r['status']}\n"
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
