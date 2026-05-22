#!/usr/bin/env python3
"""Pair up_anchor / down_anchor hits per (reference, target_genome) and
compute the inter-anchor distance distribution.

The distance distribution per reference tells us:
  - peak near 0:   EMPTY site (no IS110 between flanks)
  - peak near ref's IS110 length: FILLED V0 (intact IS110)
  - > IS110 length: FILLED V1+ (extra DNA inserted)
  - between 0 and IS110 length: deletion within IS110

Input PAF: anchors aligned against a genome database (minimap2).
Anchor names: <ref_id>__up  and  <ref_id>__down

Usage:
    find_anchor_pairs.py --paf anchors_vs_db.paf --ref-table ref_table.tsv --out OUT
"""

import argparse
import csv
import os
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--paf", required=True, help="PAF: anchors vs genome database")
    p.add_argument("--ref-table", required=True,
                   help="Ref table from extract_anchors.py "
                        "(is_id, ref_len, anchor_len, tnp_len, flank_size, distances)")
    p.add_argument("--out", required=True,
                   help="Output prefix. Writes <out>_pairs.tsv, <out>_summary.tsv")
    p.add_argument("--min-identity", type=float, default=95,
                   help="Min %% identity (default: 95)")
    p.add_argument("--min-coverage", type=float, default=80,
                   help="Min %% query coverage (default: 80)")
    p.add_argument("--max-pair-distance", type=int, default=200000,
                   help="Max distance to consider a valid anchor pair (default: 200000)")
    return p.parse_args()


def main():
    args = parse_args()

    # Load reference info from ref_table.tsv
    # Columns: is_id, ref_len, anchor_len, tnp_len, flank_size, distances_csv
    # For a single fixed D (legacy mode), expected_v0 and expected_empty are
    # computed from ref_len + anchor_len; but with multi-D anchors, we let
    # call_boundaries.py do the proper per-D peak detection.
    ref_info = {}
    with open(args.ref_table) as f:
        header = f.readline().rstrip().split("\t")
        col = {n: i for i, n in enumerate(header)}
        for line in f:
            parts = line.rstrip().split("\t")
            if len(parts) < 5: continue
            is_id = parts[col["is_id"]] if "is_id" in col else parts[0]
            ref_len = int(parts[col.get("ref_len", 1)])
            anchor_len = int(parts[col.get("anchor_len", 2)])
            tnp_len = int(parts[col.get("tnp_len", 3)])
            flank = int(parts[col.get("flank_size", 4)])
            # IS element size unknown yet; use a placeholder
            is_len_estimate = tnp_len
            expected_v0 = ref_len - 2 * anchor_len
            expected_empty = expected_v0 - is_len_estimate
            ref_info[is_id] = {
                "ref_len": ref_len,
                "anchor_len": anchor_len,
                "tnp_len": tnp_len,
                "flank": flank,
                "is_len": is_len_estimate,
                "expected_v0": expected_v0,
                "expected_empty": expected_empty,
            }

    # Parse PAF, group hits by (ref_id, anchor_side, target_contig)
    # hits = {(ref_id, side): {target_contig: [(strand, ts, te, ident)...]}}
    hits = defaultdict(lambda: defaultdict(list))

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

            # Anchor name format: <ref_id>__up or <ref_id>__down
            if "__up" in qname:
                ref_id = qname.rsplit("__up", 1)[0]
                side = "up"
            elif "__down" in qname:
                ref_id = qname.rsplit("__down", 1)[0]
                side = "down"
            else:
                continue

            hits[(ref_id, side)][tname].append({
                "strand": strand,
                "ts": ts,
                "te": te,
                "ident": ident,
                "tlen": tlen,
            })

    print(f"PAF rows: {n_total:,}; passed filter: {n_kept:,}")

    # Pair up anchors per (ref_id, target_contig)
    # An "anchor pair" = one up_hit + one down_hit on the same target_contig,
    # consistent strand, within max_pair_distance.
    out_pairs = open(f"{args.out}_pairs.tsv", "w")
    out_pairs.write("ref_id\tref_len\texpected_v0\texpected_empty\ttarget_contig\ttlen\t"
                    "up_strand\tup_ts\tup_te\t"
                    "down_strand\tdown_ts\tdown_te\t"
                    "distance\tcategory\n")

    # Aggregate per ref
    ref_summary = defaultdict(lambda: {
        "n_pairs": 0,
        "n_empty": 0,
        "n_v0": 0,
        "n_v1plus": 0,
        "n_deletion": 0,
        "n_ambiguous": 0,
        "distances": [],
    })

    all_refs = set(rid for rid, _ in hits.keys())
    for ref_id in all_refs:
        if ref_id not in ref_info:
            continue
        info = ref_info[ref_id]
        expected_v0 = info["expected_v0"]
        expected_empty = info["expected_empty"]
        tol_v0 = max(200, expected_v0 * 0.10)
        tol_empty = max(200, expected_empty * 0.15)

        up_hits = hits.get((ref_id, "up"), {})
        down_hits = hits.get((ref_id, "down"), {})

        for tname in set(up_hits.keys()) & set(down_hits.keys()):
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

                    # Categorize relative to expected_v0 and expected_empty
                    if abs(distance - expected_empty) <= tol_empty:
                        category = "empty"
                    elif abs(distance - expected_v0) <= tol_v0:
                        category = "v0_filled"
                    elif distance > expected_v0 + tol_v0:
                        category = "v1plus_filled"
                    elif expected_empty + tol_empty < distance < expected_v0 - tol_v0:
                        category = "intermediate"
                    elif distance < expected_empty - tol_empty:
                        category = "very_short"
                    else:
                        category = "ambiguous"

                    out_pairs.write(
                        f"{ref_id}\t{info['ref_len']}\t{expected_v0}\t{expected_empty}\t"
                        f"{tname}\t{tlen}\t"
                        f"{up['strand']}\t{up['ts']}\t{up['te']}\t"
                        f"{down['strand']}\t{down['ts']}\t{down['te']}\t"
                        f"{distance}\t{category}\n"
                    )

                    ref_summary[ref_id]["n_pairs"] += 1
                    ref_summary[ref_id]["distances"].append(distance)
                    key = f"n_{category}"
                    ref_summary[ref_id][key] = ref_summary[ref_id].get(key, 0) + 1

    out_pairs.close()
    print(f"Wrote pairs: {args.out}_pairs.tsv")

    # Write summary
    with open(f"{args.out}_summary.tsv", "w") as f:
        f.write("ref_id\tref_len\texpected_v0\texpected_empty\tn_pairs\t"
                "n_empty\tn_v0\tn_v1plus\tn_intermediate\tn_very_short\tn_ambiguous\t"
                "has_empty\thas_v0\tbimodal\n")
        n_with_empty = 0
        n_with_v0 = 0
        n_bimodal = 0
        for ref_id, summary in sorted(ref_summary.items()):
            info = ref_info.get(ref_id, {})
            n_empty = summary.get("n_empty", 0)
            n_v0 = summary.get("n_v0_filled", 0)
            n_v1plus = summary.get("n_v1plus_filled", 0)
            n_inter = summary.get("n_intermediate", 0)
            n_short = summary.get("n_very_short", 0)
            n_amb = summary.get("n_ambiguous", 0)
            has_empty = n_empty >= 1
            has_v0 = n_v0 >= 1
            bimodal = has_empty and (has_v0 or n_v1plus >= 1)
            if has_empty: n_with_empty += 1
            if has_v0: n_with_v0 += 1
            if bimodal: n_bimodal += 1
            f.write(f"{ref_id}\t{info.get('ref_len', '')}\t"
                    f"{info.get('expected_v0', '')}\t{info.get('expected_empty', '')}\t"
                    f"{summary['n_pairs']}\t{n_empty}\t{n_v0}\t{n_v1plus}\t"
                    f"{n_inter}\t{n_short}\t{n_amb}\t"
                    f"{1 if has_empty else 0}\t{1 if has_v0 else 0}\t"
                    f"{1 if bimodal else 0}\n")

    print(f"Wrote summary: {args.out}_summary.tsv")
    print(f"\n=== Cross-reference summary ===")
    print(f"Refs with any pairs:        {len(ref_summary):,}")
    print(f"Refs with empty sites:      {n_with_empty:,}")
    print(f"Refs with V0 sites:         {n_with_v0:,}")
    print(f"Refs with both (bimodal):   {n_bimodal:,}  <-- can refine boundaries")


if __name__ == "__main__":
    main()
