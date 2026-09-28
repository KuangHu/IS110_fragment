#!/usr/bin/env python3
"""Filter translocation candidates to chromosome-level targets only, then
deduplicate per IS locus (one row per ref_id) keeping the full supporting list.

Reads `<species>_translocation_all.tsv` (output of find_rearrangements.py with
the new --fai annotation), and a `complete_targets.txt` file with one assembly
accession per line (passing the completeness criterion).

Workflow (matches the Burkholderia paper 2018 design):
  1. Drop every row whose target assembly is NOT in the complete set.
  2. Group surviving rows by ref_id.
  3. Within each group, pick the representative with the deepest anchor pair
     (max `min_anchor_depth`).
  4. Emit a TSV row (validator-compatible columns) + the full supporting
     target list (compact column + detailed JSON).

Output:
  <out>.tsv                 dedup'd TSV ready for validate_rearrangements_gold.py
  <out>_supporting.json     full per-locus supporting-target list (recoverable)

Usage:
  dedup_translocations_per_locus.py --in <sp>_translocation_all.tsv \\
      --complete-targets complete_targets.txt --out translocation_candidates_dedup.tsv
"""
import argparse, csv, json, os
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--in", dest="inp", required=True)
    p.add_argument("--complete-targets", required=True,
                   help="text file: one assembly accession per line")
    p.add_argument("--min-anchor-depth", type=int, default=0,
                   help="extra filter: require min_anchor_depth >= this (default 0; "
                        "set e.g. 5000 to also keep only deep-anchor candidates)")
    p.add_argument("--out", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    complete = set()
    with open(args.complete_targets) as f:
        for line in f:
            t = line.strip()
            if t:
                complete.add(t)
    print(f"complete-target set: {len(complete):,} assemblies")

    by_ref = defaultdict(list)
    n_in = 0
    n_after_complete = 0
    n_after_depth = 0
    with open(args.inp) as f:
        rd = csv.DictReader(f, delimiter="\t")
        for row in rd:
            n_in += 1
            if row["assembly"] not in complete:
                continue
            n_after_complete += 1
            try:
                depth = int(row.get("min_anchor_depth", "0") or 0)
            except ValueError:
                depth = 0
            if depth < args.min_anchor_depth:
                continue
            n_after_depth += 1
            row["_depth"] = depth
            by_ref[row["ref_id"]].append(row)

    print(f"input rows:                              {n_in:,}")
    print(f"after complete-target filter:            {n_after_complete:,}")
    print(f"after depth filter (>= {args.min_anchor_depth}):  {n_after_depth:,}")
    print(f"unique IS loci (groups):                 {len(by_ref):,}")

    out_cols = ["category", "ref_id", "assembly", "up_contigs", "down_contigs",
                "details", "min_anchor_depth", "up_contig_len", "down_contig_len",
                "n_supporting_targets", "supporting_targets"]
    n_pairs = 0
    full = {}
    with open(args.out, "w") as fh:
        fh.write("\t".join(out_cols) + "\n")
        for ref_id, rows in by_ref.items():
            rows.sort(key=lambda r: -r["_depth"])
            rep = rows[0]
            n_pairs += len(rows)
            support_compact = ";".join(f"{r['assembly']}:{r['_depth']}" for r in rows)
            fh.write("\t".join([
                rep["category"], rep["ref_id"], rep["assembly"],
                rep["up_contigs"], rep["down_contigs"], rep["details"],
                str(rep["_depth"]), str(rep.get("up_contig_len","")),
                str(rep.get("down_contig_len","")), str(len(rows)),
                support_compact,
            ]) + "\n")
            full[ref_id] = {
                "representative_assembly": rep["assembly"],
                "representative_min_anchor_depth": rep["_depth"],
                "n_supporting_targets": len(rows),
                "supporting_targets": [
                    {
                        "assembly": r["assembly"],
                        "min_anchor_depth": r["_depth"],
                        "up_contigs": r["up_contigs"],
                        "down_contigs": r["down_contigs"],
                        "details": r["details"],
                        "up_contig_len": r.get("up_contig_len",""),
                        "down_contig_len": r.get("down_contig_len",""),
                    } for r in rows
                ],
            }
    side = args.out.replace(".tsv", "_supporting.json")
    with open(side, "w") as f:
        json.dump(full, f, indent=2)

    print(f"\noutput TSV: {args.out}")
    print(f"output supporting JSON: {side}")


if __name__ == "__main__":
    main()
