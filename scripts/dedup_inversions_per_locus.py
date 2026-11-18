#!/usr/bin/env python3
"""Deduplicate inversion candidates to one event per IS locus (ref_id).

Reads `<species>_inversion_all.tsv` (output of find_rearrangements.py with the
all-inversion catalogue), filters to estimated_inverted_bp >= --min-bp, and
collapses rows with the same ref_id into a single representative — the one
with the largest estimated_inverted_bp.  Also records n_supporting_targets =
number of distinct target assemblies that independently showed the inversion
signature at this locus.

Output schema (compatible with validate_rearrangements_gold.py's --examples):
  category  ref_id  assembly  up_contigs  down_contigs  details
plus extras:
  estimated_inverted_bp  n_supporting_targets

Usage:
  dedup_inversions_per_locus.py --in <species>_inversion_all.tsv \\
      --min-bp 10000 --out inversion_candidates_dedup.tsv
"""
import argparse, csv
from collections import defaultdict
import os, sys


sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib_rearr as lr  # noqa: E402
from lib_alleles import assembly_core  # noqa: E402


def event_points(row):
    """(assembly core, sorted [(contig, position)] of the first up/down anchor)."""
    d = lr.parse_details(row.get("details", ""))
    return assembly_core(row["assembly"]), sorted((c, s) for c, s, e, st in d["up"][:1] + d["down"][:1])


def write_events(rows_by_ref, path, tol=6000):
    """One row per EVENT in a target genome.

    The same inversion / translocation is reached from IS_a and from IS_b
    (anchor positions a few kb apart: one IS plus the anchor offsets), from
    every reference strain carrying the element, and from GCA/GCF twins. Rows
    whose endpoints are ALL within `tol` bp of each other on the same contigs
    of the same assembly core are one event. (Rounding positions into bins
    split one inversion seen from its two IS copies and could merge events
    on either side of a bin edge; two real events closer than `tol` at both
    ends are still merged -- reported as n_merged_rows.)"""
    items = []
    for ref_id, rows in rows_by_ref.items():
        for r in rows:
            core, pts = event_points(r)
            if pts:
                items.append((core, tuple(c for c, _ in pts), tuple(p for _, p in pts), ref_id))
    items.sort()
    events = []          # [core, contigs, [positions...], {ref_ids}, n_rows]
    for core, contigs, pos, ref_id in items:
        for ev in events[::-1]:
            if ev[0] != core or ev[1] != contigs:
                continue
            if all(abs(a - b) <= tol for a, b in zip(ev[2], pos)):
                ev[3].add(ref_id); ev[4] += 1
                break
        else:
            events.append([core, contigs, pos, {ref_id}, 1])
    with open(path, "w") as fh:
        fh.write("target_assembly_core\tanchor_positions\tn_ref_ids\tn_merged_rows\tref_ids\n")
        for core, contigs, pos, refs, n in events:
            where = ";".join(f"{c}:~{p}" for c, p in zip(contigs, pos))
            fh.write(f"{core}\t{where}\t{len(refs)}\t{n}\t{','.join(sorted(refs))}\n")
    return events


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--in", dest="inp", required=True,
                   help="rearrangements_inversion_all.tsv from find_rearrangements.py")
    p.add_argument("--min-bp", type=int, default=10000,
                   help="estimated_inverted_bp floor (default 10000)")
    p.add_argument("--out", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    # Group rows by ref_id, keep best per ref
    by_ref = defaultdict(list)
    with open(args.inp) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            try:
                est = int(row.get("estimated_inverted_bp", "0") or 0)
            except ValueError:
                est = 0
            if est < args.min_bp:
                continue
            row["_est"] = est
            by_ref[row["ref_id"]].append(row)

    # one representative per ref_id (largest est inverted span)
    # ALSO record every supporting target assembly and its estimate, so no
    # genome that witnessed the same inversion is lost.
    out_cols = ["category", "ref_id", "assembly", "up_contigs", "down_contigs",
                "details", "estimated_inverted_bp", "n_supporting_targets",
                "supporting_targets"]
    n_loci = 0
    n_pairs_total = 0

    # Separate JSON file with the FULL per-target list (assembly, est_bp, contig,
    # anchor coords from `details`) so the supporting set is fully recoverable.
    full_json_path = os.path.splitext(args.out)[0] + "_supporting.json"
    full_dict = {}

    with open(args.out, "w") as fh:
        fh.write("\t".join(out_cols) + "\n")
        for ref_id, rows in sorted(by_ref.items()):
            rows.sort(key=lambda r: (-r["_est"], r["assembly"]))
            rep = rows[0]
            n_loci += 1
            n_pairs_total += len(rows)
            # compact TSV summary: "asm:est_bp;asm:est_bp;..."
            support_compact = ";".join(f"{r['assembly']}:{r['_est']}" for r in rows)
            fh.write("\t".join([
                rep["category"], rep["ref_id"], rep["assembly"],
                rep["up_contigs"], rep["down_contigs"], rep["details"],
                str(rep["_est"]), str(len({assembly_core(r["assembly"]) for r in rows})), support_compact,
            ]) + "\n")
            # detailed JSON entry
            full_dict[ref_id] = {
                "representative_assembly": rep["assembly"],
                "representative_estimated_inverted_bp": rep["_est"],
                "n_supporting_targets": len(rows),
                "supporting_targets": [
                    {
                        "assembly": r["assembly"],
                        "estimated_inverted_bp": r["_est"],
                        "up_contigs": r["up_contigs"],
                        "down_contigs": r["down_contigs"],
                        "details": r["details"],
                    } for r in rows
                ],
            }

    import json
    with open(full_json_path, "w") as f:
        json.dump(full_dict, f, indent=2)

    ev_path = os.path.splitext(args.out)[0] + "_events.tsv"
    ev = write_events(by_ref, ev_path)
    print(f"input pairs (>= {args.min_bp} bp): {n_pairs_total:,}")
    print(f"reference IS loci with a candidate:  {n_loci:,}   (NOT a count of inversions)")
    print(f"distinct target events:              {len(ev):,}   -> {ev_path}")
    print(f"compression ratio:                  {n_pairs_total / max(n_loci,1):.1f}x")
    print(f"output TSV:                         {args.out}")
    print(f"output supporting JSON (full):      {full_json_path}")


if __name__ == "__main__":
    main()
