#!/usr/bin/env python3
"""Was the IS of each reference found? (used by search_refs_in_logan.sh)

A reference is transposase + 2 kb of flank on each side. A hit that lies in a
flank only says that the neighbouring genes occur in the target, not the IS:
the old summary counted any >=95% hit, so a 2-kb hit entirely in the left
flank was "IS detected". Here the IS counts as detected only if one alignment
row (one locus) of >= --min-identity aligns >= --min-core-cov % of the CORE.

The core interval comes from `core=<start>-<end>` in the reference header
(0-based half-open, written by extract_is110_refs.sh). Without it the whole
reference is taken as core and a warning is printed.

Usage: ref_core_coverage.py <refs_vs_target.paf> <refs.fa> <out.tsv>
"""
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paf_utils import aligned_identity, aligned_blocks  # noqa: E402


def read_refs(path):
    refs, name = {}, None
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].split()[0]
                m = re.search(r"\bcore=(\d+)-(\d+)", line)
                refs[name] = {"len": 0, "core": (int(m.group(1)), int(m.group(2))) if m else None}
            elif name:
                refs[name]["len"] += len(line.strip())
    return refs


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("paf")
    p.add_argument("refs")
    p.add_argument("out_tsv")
    p.add_argument("--min-identity", type=float, default=95)
    p.add_argument("--min-core-cov", type=float, default=90)
    a = p.parse_args()

    refs = read_refs(a.refs)
    no_core = sum(1 for r in refs.values() if r["core"] is None)
    if no_core:
        print(f"WARNING: {no_core} references have no core= in the header; the whole "
              "sequence is used as core for them")
    for r in refs.values():
        if r["core"] is None:
            r["core"] = (0, r["len"])
        r.update(n_hits=0, n_hits_id=0, best_core=0.0, best_any=0.0)

    with open(a.paf) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 12 or c[0] not in refs:
                continue
            r = refs[c[0]]
            r["n_hits"] += 1
            if aligned_identity(c) < a.min_identity:
                continue
            r["n_hits_id"] += 1
            cs, ce = r["core"]
            blocks = aligned_blocks(c, 1, 1)          # cut at every indel: aligned bases only
            core_bp = sum(max(0, min(qe, ce) - max(qs, cs)) for qs, qe, _, _ in blocks)
            all_bp = sum(qe - qs for qs, qe, _, _ in blocks)
            r["best_core"] = max(r["best_core"], 100.0 * core_bp / max(1, ce - cs))
            r["best_any"] = max(r["best_any"], 100.0 * all_bp / max(1, r["len"]))

    n = len(refs)
    detected = sum(1 for r in refs.values() if r["best_core"] >= a.min_core_cov)
    partial = sum(1 for r in refs.values() if 0 < r["best_core"] < a.min_core_cov)
    flank_only = sum(1 for r in refs.values() if r["n_hits_id"] and r["best_core"] == 0)
    low_id = sum(1 for r in refs.values() if r["n_hits"] and not r["n_hits_id"])
    none = sum(1 for r in refs.values() if not r["n_hits"])
    with open(a.out_tsv, "w") as out:
        out.write("reference\tref_len\tcore_start\tcore_end\tn_hits\tn_hits_identity\t"
                  "best_core_cov_pct\tbest_ref_cov_pct\tstatus\n")
        for name in sorted(refs):
            r = refs[name]
            status = ("core_detected" if r["best_core"] >= a.min_core_cov else
                      "core_partial" if r["best_core"] > 0 else
                      "flank_only" if r["n_hits_id"] else
                      "low_identity_only" if r["n_hits"] else "no_hit")
            out.write(f"{name}\t{r['len']}\t{r['core'][0]}\t{r['core'][1]}\t{r['n_hits']}\t"
                      f"{r['n_hits_id']}\t{r['best_core']:.1f}\t{r['best_any']:.1f}\t{status}\n")

    def pct(x):
        return f"{x:,} ({100.0 * x / n:.1f}%)" if n else "0"
    print(f"\n=== {n:,} references ===")
    print(f"  IS core detected (>= {a.min_core_cov:g}% of the core in one row, "
          f">= {a.min_identity:g}% identity): {pct(detected)}")
    print(f"  core partly aligned:                 {pct(partial)}")
    print(f"  hit in the flanks only (NOT the IS): {pct(flank_only)}")
    print(f"  hits below the identity threshold:   {pct(low_id)}")
    print(f"  no hit:                              {pct(none)}")
    print(f"  per reference: {a.out_tsv}")


if __name__ == "__main__":
    main()
