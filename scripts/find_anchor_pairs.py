#!/usr/bin/env python3
"""Pair up_anchor / down_anchor hits per (reference, target assembly) and
compute the inter-anchor distance distribution.

The distance relative to the source (filled) state tells us:
  - same as source:            FILLED V0 (intact element)
  - longer than source:        FILLED V1+ (extra DNA inserted)
  - shorter by >= tnp length:  EMPTY site (element absent)
  - shorter by less:           intermediate (partial deletion)

One placement per distinct assembly (GCA/GCF twins collapsed), chosen by
lib_alleles.place_pair; ambiguous placements are written but not counted as
empty/filled.

Input PAF: anchors aligned against a genome database (minimap2).
Anchor names: <ref_id>__up<D> / <ref_id>__down<D>  (or legacy <ref_id>__up)

Usage:
    find_anchor_pairs.py --paf anchors_vs_db.paf --ref-table ref_table.tsv --out OUT
"""

import argparse
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_alleles import (load_anchor_hits, place_per_assembly,  # noqa: E402
                         gap_tolerance, classify_gap)


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
    p.add_argument("--max-extra-bp", type=int, default=50000,
                   help="Max extra DNA beyond the source state for a valid pair "
                        "(default: 50000)")
    p.add_argument("--tol-bp", type=int, default=100,
                   help="Absolute distance tolerance (default: 100)")
    p.add_argument("--tol-frac", type=float, default=0.002,
                   help="Extra tolerance as a fraction of expected distance "
                        "(default: 0.002)")
    p.add_argument("--margin", type=float, default=0.05,
                   help="Min relative score margin over the runner-up pair, else "
                        "the placement is ambiguous (default: 0.05)")
    return p.parse_args()


def main():
    args = parse_args()

    # Columns: is_id, ref_len, anchor_len, tnp_len, flank_size, distances_csv
    ref_info = {}
    with open(args.ref_table) as f:
        header = f.readline().rstrip().split("\t")
        col = {n: i for i, n in enumerate(header)}
        for line in f:
            parts = line.rstrip().split("\t")
            if len(parts) < 5:
                continue
            is_id = parts[col["is_id"]] if "is_id" in col else parts[0]
            ref_info[is_id] = {
                "ref_len": int(parts[col.get("ref_len", 1)]),
                "anchor_len": int(parts[col.get("anchor_len", 2)]),
                "tnp_len": int(parts[col.get("tnp_len", 3)]),
            }

    hits, n_total, n_kept = load_anchor_hits(args.paf, args.min_identity,
                                             args.min_coverage)
    print(f"PAF rows: {n_total:,}; passed filter: {n_kept:,}")

    out_pairs = open(f"{args.out}_pairs.tsv", "w")
    out_pairs.write("ref_id\tref_len\texpected_v0\texpected_empty\ttarget_contig\ttlen\t"
                    "up_strand\tup_ts\tup_te\t"
                    "down_strand\tdown_ts\tdown_te\t"
                    "distance\tcategory\t"
                    "distance_D\tassembly\tplacement_status\tn_candidate_pairs\tmargin\n")

    ref_summary = defaultdict(lambda: defaultdict(int))
    keys = sorted({(r, D) for (r, _, D) in hits}, key=lambda k: (k[0], k[1] or 0))
    for ref_id, D in keys:
        if ref_id not in ref_info:
            continue
        info = ref_info[ref_id]
        if D is None:   # legacy single-D reference
            expected_v0 = info["ref_len"] - 2 * info["anchor_len"]
        else:
            expected_v0 = 2 * D + info["tnp_len"]
        expected_empty = expected_v0 - info["tnp_len"]
        tol = gap_tolerance(expected_v0, args.tol_bp, args.tol_frac)

        placements = place_per_assembly(
            hits.get((ref_id, "up", D), {}), hits.get((ref_id, "down", D), {}),
            min_gap=-50, max_gap=expected_v0 + args.max_extra_bp,
            margin_ratio=args.margin)
        for pl in placements.values():
            if pl["status"] == "ambiguous":
                category = "ambiguous"
            else:
                category = classify_gap(pl["gap"], expected_v0, info["tnp_len"], tol)
            u, d = pl["up"], pl["down"]
            out_pairs.write(
                f"{ref_id}\t{info['ref_len']}\t{expected_v0}\t{expected_empty}\t"
                f"{pl['tname']}\t{u.tlen}\t"
                f"{u.strand}\t{u.ts}\t{u.te}\t"
                f"{d.strand}\t{d.ts}\t{d.te}\t"
                f"{pl['gap']}\t{category}\t"
                f"{'' if D is None else D}\t{pl['assembly']}\t{pl['status']}\t"
                f"{pl['n_pairs']}\t{pl['margin']:.3f}\n")
            s = ref_summary[ref_id]
            s["n_pairs"] += 1
            s[category] += 1

    out_pairs.close()
    print(f"Wrote pairs: {args.out}_pairs.tsv")

    # Per-ref totals pooled over D (call_boundaries.py reports per D)
    with open(f"{args.out}_summary.tsv", "w") as f:
        f.write("ref_id\tref_len\tn_pairs\t"
                "n_empty\tn_v0\tn_v1plus\tn_intermediate\tn_ambiguous\t"
                "has_empty\thas_v0\tbimodal\n")
        n_with_empty = n_with_v0 = n_bimodal = 0
        for ref_id, s in sorted(ref_summary.items()):
            n_empty, n_v0, n_v1plus = s["empty"], s["v0_filled"], s["v1plus_filled"]
            has_empty, has_v0 = n_empty >= 1, n_v0 >= 1
            bimodal = has_empty and (has_v0 or n_v1plus >= 1)
            n_with_empty += has_empty
            n_with_v0 += has_v0
            n_bimodal += bimodal
            f.write(f"{ref_id}\t{ref_info[ref_id]['ref_len']}\t{s['n_pairs']}\t"
                    f"{n_empty}\t{n_v0}\t{n_v1plus}\t"
                    f"{s['intermediate']}\t{s['ambiguous']}\t"
                    f"{int(has_empty)}\t{int(has_v0)}\t{int(bimodal)}\n")

    print(f"Wrote summary: {args.out}_summary.tsv")
    print(f"\n=== Cross-reference summary (distinct assemblies) ===")
    print(f"Refs with any pairs:        {len(ref_summary):,}")
    print(f"Refs with empty sites:      {n_with_empty:,}")
    print(f"Refs with V0 sites:         {n_with_v0:,}")
    print(f"Refs with both (bimodal):   {n_bimodal:,}  <-- can refine boundaries")


if __name__ == "__main__":
    main()
