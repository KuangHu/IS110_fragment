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
import re
from collections import defaultdict


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
                   help="Max expected IS element size (bp) for peak filter")
    return p.parse_args()


# Anchor name regex: <ref_id>__up<D>  or  <ref_id>__down<D>
ANCHOR_RE = re.compile(r"^(.+)__(up|down)(\d+)$")


def main():
    args = parse_args()

    # Load reference table
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

    # Parse PAF: organize hits by (ref_id, side, distance)
    hits = defaultdict(lambda: defaultdict(list))  # (ref_id, side, dist) -> {tname: [hits]}

    n_total = 0
    n_kept = 0
    with open(args.paf) as f:
        for line in f:
            c = line.rstrip().split("\t")
            if len(c) < 12:
                continue
            n_total += 1
            qname = c[0]
            qlen = int(c[1])
            qs, qe = int(c[2]), int(c[3])
            strand = c[4]
            tname = c[5]
            tlen = int(c[6])
            ts, te = int(c[7]), int(c[8])
            matches, block = int(c[9]), int(c[10])
            ident = matches / block * 100 if block > 0 else 0
            cov = (qe - qs) / qlen * 100
            if ident < args.min_identity or cov < args.min_coverage:
                continue
            n_kept += 1

            m = ANCHOR_RE.match(qname)
            if not m:
                continue
            ref_id = m.group(1)
            side = m.group(2)
            dist = int(m.group(3))

            hits[(ref_id, side, dist)][tname].append({
                "strand": strand,
                "ts": ts,
                "te": te,
                "ident": ident,
                "tlen": tlen,
            })

    print(f"PAF rows: {n_total:,}; passed filter: {n_kept:,}")

    # Pair up_D with down_D per ref, per target contig
    # For each (ref, D, target_contig), find all up-down pairs and their distances
    out_pairs = open(f"{args.out}_pairs.tsv", "w")
    out_pairs.write("ref_id\tdistance_D\ttnp_len\tflank_size\texpected_v0\texpected_empty\t"
                    "target_contig\ttlen\tup_strand\tup_ts\tup_te\t"
                    "down_strand\tdown_ts\tdown_te\tobserved_distance\tcategory\n")

    summary = defaultdict(lambda: defaultdict(int))  # (ref_id, D) -> {category: count}

    for ref_id, info in ref_info.items():
        for D in info["distances"]:
            up_hits = hits.get((ref_id, "up", D), {})
            down_hits = hits.get((ref_id, "down", D), {})

            # Expected distances:
            # In reference, up_anchor ends at flank_size - D, down_anchor starts at flank_size + tnp_len + D
            # V0 (filled) = (flank_size + tnp_len + D) - (flank_size - D) = 2D + tnp_len
            # Empty (IS110 excised, only transposase removed) = 2D
            # But the IS110 might extend beyond transposase; this is what we're trying to discover.
            expected_v0 = 2 * D + info["tnp_len"]
            expected_empty = 2 * D
            tol_v0 = max(200, expected_v0 * 0.10)
            tol_empty = max(200, expected_empty * 0.10)

            shared_targets = set(up_hits.keys()) & set(down_hits.keys())
            for tname in shared_targets:
                tlen = up_hits[tname][0]["tlen"]
                for up in up_hits[tname]:
                    for down in down_hits[tname]:
                        if up["strand"] != down["strand"]:
                            continue
                        if up["strand"] == "+":
                            distance = down["ts"] - up["te"]
                        else:
                            distance = up["ts"] - down["te"]
                        if abs(distance) > args.max_pair_distance:
                            continue

                        # Categorize
                        if abs(distance - expected_empty) <= tol_empty:
                            category = "empty"
                        elif abs(distance - expected_v0) <= tol_v0:
                            category = "v0_filled"
                        elif distance > expected_v0 + tol_v0:
                            category = "v1plus_filled"
                        elif expected_empty + tol_empty < distance < expected_v0 - tol_v0:
                            # IS110 element bigger than transposase but smaller than full V0?
                            # Or partial deletion? Mark as intermediate (could be IS110 with no cargo)
                            category = "intermediate"
                        elif distance < expected_empty - tol_empty:
                            category = "very_short"
                        else:
                            category = "ambiguous"

                        out_pairs.write(
                            f"{ref_id}\t{D}\t{info['tnp_len']}\t{info['flank_size']}\t"
                            f"{expected_v0}\t{expected_empty}\t"
                            f"{tname}\t{tlen}\t"
                            f"{up['strand']}\t{up['ts']}\t{up['te']}\t"
                            f"{down['strand']}\t{down['ts']}\t{down['te']}\t"
                            f"{distance}\t{category}\n"
                        )
                        summary[(ref_id, D)][category] += 1
                        summary[(ref_id, D)]["_total"] += 1

    out_pairs.close()
    print(f"Wrote pairs: {args.out}_pairs.tsv")

    # Write summary: one row per (ref_id, D)
    with open(f"{args.out}_summary.tsv", "w") as f:
        f.write("ref_id\tdistance_D\ttnp_len\tn_pairs\tn_empty\tn_v0\tn_v1plus\t"
                "n_intermediate\tn_very_short\tn_ambiguous\thas_empty\thas_v0\tbimodal\n")
        for (ref_id, D), cats in sorted(summary.items()):
            info = ref_info[ref_id]
            n_total = cats.get("_total", 0)
            n_empty = cats.get("empty", 0)
            n_v0 = cats.get("v0_filled", 0)
            n_v1plus = cats.get("v1plus_filled", 0)
            n_inter = cats.get("intermediate", 0)
            n_short = cats.get("very_short", 0)
            n_amb = cats.get("ambiguous", 0)
            has_empty = n_empty >= 1
            has_v0 = n_v0 >= 1
            bimodal = has_empty and (has_v0 or n_v1plus >= 1)
            f.write(f"{ref_id}\t{D}\t{info['tnp_len']}\t{n_total}\t"
                    f"{n_empty}\t{n_v0}\t{n_v1plus}\t{n_inter}\t{n_short}\t{n_amb}\t"
                    f"{1 if has_empty else 0}\t{1 if has_v0 else 0}\t"
                    f"{1 if bimodal else 0}\n")

    print(f"Wrote summary: {args.out}_summary.tsv")

    # Per-reference: find smallest D with empty sites
    # That's our estimate of IS110 boundary
    print(f"\n=== IS110 boundary estimates per reference ===")
    print(f"(Smallest D with empty sites suggests IS110 doesn't extend past D)")
    print(f"")
    by_ref = defaultdict(list)
    for (ref_id, D), cats in summary.items():
        if cats.get("empty", 0) > 0:
            by_ref[ref_id].append(D)

    for ref_id in sorted(ref_info.keys()):
        if ref_id in by_ref:
            min_d = min(by_ref[ref_id])
            print(f"  {ref_id[:40]:40s}  IS110 boundary <= {min_d:>6d} bp from tnp")
        else:
            print(f"  {ref_id[:40]:40s}  IS110 boundary > max distance (or no detectable empty sites)")


if __name__ == "__main__":
    main()
