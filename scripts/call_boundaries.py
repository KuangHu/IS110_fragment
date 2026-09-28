#!/usr/bin/env python3
"""Analyze multi-distance anchor pairs to discover IS110 element boundaries.

For each (ref_id, distance D), tally up empty/filled site counts.
The smallest D where empty sites appear gives the IS110 boundary in that direction.

We also support asymmetric pairing: up-anchor at one distance, down-anchor at another.
But initially focus on symmetric: same D on both sides.

Usage:
    call_boundaries.py --paf anchors_vs_db.paf --ref-table ref_table.tsv --out OUT
"""

import argparse
import os
import statistics
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_alleles import (load_anchor_hits, place_per_assembly,  # noqa: E402
                         gap_tolerance, classify_gap)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--paf", required=True, help="PAF: anchors vs DB")
    p.add_argument("--ref-table", required=True,
                   help="ref_table.tsv from extract_anchors.py")
    p.add_argument("--out", required=True,
                   help="Output prefix (writes <out>_pairs.tsv, <out>_summary.tsv)")
    p.add_argument("--min-identity", type=float, default=95)
    p.add_argument("--min-coverage", type=float, default=80)
    p.add_argument("--max-pair-distance", type=int, default=200000,
                   help="Max plausible distance between anchors on target (default: 200kb)")
    p.add_argument("--histogram-bin", type=int, default=50,
                   help="Distance histogram bin size in bp (default: 50)")
    p.add_argument("--min-is-size", type=int, default=800,
                   help="Min expected IS element size (bp) for peak filter")
    p.add_argument("--max-is-size", type=int, default=200000,
                   help="Max expected IS element size (bp); larger removals are "
                        "'very_short'")
    p.add_argument("--max-extra-bp", type=int, default=50000,
                   help="Max extra DNA beyond the source state for a valid pair")
    p.add_argument("--tol-bp", type=int, default=100,
                   help="Absolute distance tolerance (default: 100). The old "
                        "tolerance was 10%% of 2D -- +/-8 kb at D=40 kb")
    p.add_argument("--tol-frac", type=float, default=0.002,
                   help="Extra tolerance as a fraction of the expected distance")
    p.add_argument("--margin", type=float, default=0.05,
                   help="Min relative score margin over the runner-up pair")
    p.add_argument("--max-secondary", type=int, default=0,
                   help="The -N used for the anchor search. Anchors with -N+1 "
                        "PAF rows hit the cap and saw only a SAMPLE of the DB; "
                        "they are flagged in the summary (0 = don't flag)")
    return p.parse_args()


def main():
    args = parse_args()

    ref_info = {}
    with open(args.ref_table) as f:
        next(f)  # header
        for line in f:
            parts = line.rstrip().split("\t")
            ref_id, ref_len, anchor_len, tnp_len, flank_size, distances = parts
            ref_info[ref_id] = {
                "ref_len": int(ref_len),
                "anchor_len": int(anchor_len),
                "tnp_len": int(tnp_len),
                "flank_size": int(flank_size),
                "distances": [int(d) for d in distances.split(",")],
            }

    row_counts = {}
    hits, n_total, n_kept = load_anchor_hits(args.paf, args.min_identity,
                                             args.min_coverage, row_counts)
    cap = args.max_secondary + 1 if args.max_secondary > 0 else None
    if cap:
        n_capped = sum(1 for n in row_counts.values() if n >= cap)
        print(f"Anchors at the minimap2 -N cap ({cap} rows): {n_capped:,} / "
              f"{len(row_counts):,} -- their counts are a sample, not a census")
    print(f"PAF rows: {n_total:,}; passed filter: {n_kept:,}")

    out_pairs = open(f"{args.out}_pairs.tsv", "w")
    out_pairs.write("ref_id\tdistance_D\ttnp_len\tflank_size\texpected_v0\texpected_empty\t"
                    "target_contig\ttlen\tup_strand\tup_ts\tup_te\t"
                    "down_strand\tdown_ts\tdown_te\tobserved_distance\tcategory\t"
                    "removed_bp\tassembly\tplacement_status\tn_candidate_pairs\tmargin\n")

    summary = defaultdict(lambda: defaultdict(int))  # (ref_id, D) -> {category: count}
    removed = defaultdict(list)                      # (ref_id, D) -> removed bp at empties

    for ref_id, info in ref_info.items():
        for D in info["distances"]:
            # up anchor ends D bp before the transposase, down anchor starts D
            # bp after it, so the source (filled) gap is 2D + tnp_len. At an
            # empty site the gap is shorter by exactly the element length.
            expected_v0 = 2 * D + info["tnp_len"]
            expected_empty = 2 * D
            tol = gap_tolerance(expected_v0, args.tol_bp, args.tol_frac)
            placements = place_per_assembly(
                hits.get((ref_id, "up", D), {}), hits.get((ref_id, "down", D), {}),
                min_gap=-50, max_gap=min(args.max_pair_distance,
                                         expected_v0 + args.max_extra_bp),
                margin_ratio=args.margin)
            for pl in placements.values():
                if pl["status"] == "ambiguous":
                    category = "ambiguous"
                else:
                    category = classify_gap(pl["gap"], expected_v0, info["tnp_len"],
                                            tol, max_removed=args.max_is_size)
                rem = expected_v0 - pl["gap"]
                if category == "empty":
                    removed[(ref_id, D)].append(rem)
                u, d = pl["up"], pl["down"]
                out_pairs.write(
                    f"{ref_id}\t{D}\t{info['tnp_len']}\t{info['flank_size']}\t"
                    f"{expected_v0}\t{expected_empty}\t"
                    f"{pl['tname']}\t{u.tlen}\t"
                    f"{u.strand}\t{u.ts}\t{u.te}\t"
                    f"{d.strand}\t{d.ts}\t{d.te}\t"
                    f"{pl['gap']}\t{category}\t{rem}\t"
                    f"{pl['assembly']}\t{pl['status']}\t{pl['n_pairs']}\t"
                    f"{pl['margin']:.3f}\n")
                summary[(ref_id, D)][category] += 1
                summary[(ref_id, D)]["_total"] += 1

    out_pairs.close()
    print(f"Wrote pairs: {args.out}_pairs.tsv")

    # One row per (ref_id, D). Counts are DISTINCT ASSEMBLIES, not contigs.
    with open(f"{args.out}_summary.tsv", "w") as f:
        f.write("ref_id\tdistance_D\ttnp_len\tn_pairs\tn_empty\tn_v0\tn_v1plus\t"
                "n_intermediate\tn_very_short\tn_ambiguous\thas_empty\thas_v0\tbimodal\t"
                "is_length_mode\tis_length_mode_support\tanchor_search_capped\n")
        for (ref_id, D), cats in sorted(summary.items()):
            info = ref_info[ref_id]
            n_empty = cats.get("empty", 0)
            n_v0 = cats.get("v0_filled", 0)
            n_v1plus = cats.get("v1plus_filled", 0)
            has_empty = n_empty >= 1
            has_v0 = n_v0 >= 1
            bimodal = has_empty and (has_v0 or n_v1plus >= 1)
            mode, support = length_mode(removed[(ref_id, D)], args.histogram_bin)
            capped = int(bool(cap) and any(
                row_counts.get(f"{ref_id}__{side}{D}", 0) >= cap
                for side in ("up", "down")))
            f.write(f"{ref_id}\t{D}\t{info['tnp_len']}\t{cats.get('_total', 0)}\t"
                    f"{n_empty}\t{n_v0}\t{n_v1plus}\t{cats.get('intermediate', 0)}\t"
                    f"{cats.get('very_short', 0)}\t{cats.get('ambiguous', 0)}\t"
                    f"{int(has_empty)}\t{int(has_v0)}\t{int(bimodal)}\t"
                    f"{'' if mode is None else mode}\t{support}\t{capped}\n")

    print(f"Wrote summary: {args.out}_summary.tsv")

    # Per-reference: smallest D with empty sites bounds the element extent
    print(f"\n=== IS boundary estimates per reference ===")
    print(f"(Smallest D with empty sites suggests the element doesn't extend past D)")
    print(f"")
    by_ref = defaultdict(list)
    for (ref_id, D), cats in summary.items():
        if cats.get("empty", 0) > 0:
            by_ref[ref_id].append(D)

    for ref_id in sorted(ref_info.keys()):
        if ref_id in by_ref:
            min_d = min(by_ref[ref_id])
            mode, _ = length_mode(removed[(ref_id, min_d)], args.histogram_bin)
            print(f"  {ref_id[:40]:40s}  boundary <= {min_d:>6d} bp from tnp; "
                  f"element ~{mode} bp")
        else:
            print(f"  {ref_id[:40]:40s}  boundary > max distance (or no detectable empty sites)")


def length_mode(values, bin_size):
    """Median of the most populated bin -> (length, n_in_bin)."""
    if not values:
        return None, 0
    bins = defaultdict(list)
    for v in values:
        bins[v // bin_size].append(v)
    best = max(bins.values(), key=len)
    return int(statistics.median(best)), len(best)


if __name__ == "__main__":
    main()
