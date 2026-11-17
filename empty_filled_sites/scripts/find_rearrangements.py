#!/usr/bin/env python3
"""Find genome rearrangements from existing anchor PAF.

For each (ref_id, target_genome), look at anchor pair patterns:
  - Same contig, same strand   → normal layout (already analyzed)
  - Same contig, OPPOSITE strands → INVERSION at that locus
  - Different contigs           → TRANSLOCATION or contig break
  - Anchor matches at multiple loci → DUPLICATION

Output examples per category for review.

Usage:
    find_rearrangements.py <paf> <out_prefix>
"""

import argparse
import re
import sys
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("paf")
    p.add_argument("out_prefix")
    p.add_argument("--min-identity", type=float, default=95)
    p.add_argument("--min-cov", type=float, default=90)
    p.add_argument("--anchor-D", type=int, default=1000,
                   help="Use this D for boundary tightness (default 1000)")
    return p.parse_args()


ANCHOR_RE = re.compile(r"^(.+)__(up|down)(\d+)$")


def main():
    args = parse_args()

    # Stream PAF, collect per (ref, side) hits at the chosen D
    # Group by target genome assembly to find pairs
    # Anchor name format: ref_id__up<D>  or  ref_id__down<D>

    # Index: (ref_id, target_assembly) -> {up: [hits], down: [hits]}
    obs = defaultdict(lambda: {"up": defaultdict(list), "down": defaultdict(list)})
    # second-level key is contig

    n_total = 0
    n_kept = 0
    n_skip_d = 0
    with open(args.paf) as f:
        for line in f:
            n_total += 1
            if n_total % 10_000_000 == 0:
                print(f"  read {n_total:,} lines", flush=True, file=sys.stderr)
            c = line.split("\t")
            if len(c) < 12:
                continue
            qname = c[0]
            qlen = int(c[1])
            qs, qe = int(c[2]), int(c[3])
            strand = c[4]
            tname = c[5]
            ts, te = int(c[7]), int(c[8])
            matches, block = int(c[9]), int(c[10])
            ident = matches / block * 100 if block > 0 else 0
            cov = (qe - qs) / qlen * 100
            if ident < args.min_identity or cov < args.min_cov:
                continue
            m = ANCHOR_RE.match(qname)
            if not m:
                continue
            ref_id, side, D = m.group(1), m.group(2), int(m.group(3))
            if D != args.anchor_D:
                n_skip_d += 1
                continue
            n_kept += 1

            # Parse target: ASSEMBLY|CONTIG
            if "|" in tname:
                assembly, contig = tname.split("|", 1)
            else:
                assembly, contig = "unknown", tname

            obs[(ref_id, assembly)][side][contig].append({
                "strand": strand, "ts": ts, "te": te, "ident": ident, "contig": contig,
            })

    print(f"PAF: {n_total:,} lines, {n_kept:,} kept (D={args.anchor_D}), {n_skip_d:,} other D",
          file=sys.stderr)
    print(f"Unique (ref, target_assembly) pairs: {len(obs):,}", file=sys.stderr)

    # Classify each (ref, assembly) interaction
    counters = defaultdict(int)
    examples = defaultdict(list)
    max_examples = 50

    for (ref_id, assembly), sides in obs.items():
        up_hits = sides["up"]
        down_hits = sides["down"]
        if not up_hits or not down_hits:
            continue

        # Did we find anchor pairs on the same contig with same strand?
        has_normal = False
        has_inversion = False
        has_translocation = False

        # Check same-contig pairings
        for contig in set(up_hits.keys()) & set(down_hits.keys()):
            for u in up_hits[contig]:
                for d in down_hits[contig]:
                    if u["strand"] == d["strand"]:
                        has_normal = True
                    else:
                        has_inversion = True

        # Different-contig pairings = potential translocation
        all_up_contigs = set(up_hits.keys())
        all_down_contigs = set(down_hits.keys())
        if all_up_contigs and all_down_contigs:
            if all_up_contigs.isdisjoint(all_down_contigs):
                has_translocation = True

        if has_normal:
            counters["normal"] += 1
        if has_inversion and not has_normal:
            counters["inversion_only"] += 1
            if len(examples["inversion_only"]) < max_examples:
                examples["inversion_only"].append((ref_id, assembly, sides))
        elif has_inversion:
            counters["inversion_plus_normal"] += 1
        if has_translocation and not has_normal and not has_inversion:
            counters["translocation_only"] += 1
            if len(examples["translocation_only"]) < max_examples:
                examples["translocation_only"].append((ref_id, assembly, sides))
        elif has_translocation:
            counters["translocation_plus"] += 1

        # Duplications: any contig with multiple matches of the same anchor side
        n_dup = sum(1 for c, hs in up_hits.items() if len(hs) > 1)
        n_dup += sum(1 for c, hs in down_hits.items() if len(hs) > 1)
        if n_dup > 0:
            counters["has_duplication"] += 1
            if len(examples["duplication"]) < max_examples:
                examples["duplication"].append((ref_id, assembly, sides))

    print(f"\n=== Rearrangement signature counts (per ref-assembly pair) ===")
    for k, v in sorted(counters.items(), key=lambda x: -x[1]):
        print(f"  {k}: {v:,}")

    # Write examples
    out_path = f"{args.out_prefix}_examples.tsv"
    with open(out_path, "w") as f:
        f.write("category\tref_id\tassembly\tup_contigs\tdown_contigs\tdetails\n")
        for cat, exs in examples.items():
            for ref_id, assembly, sides in exs:
                up_contigs = ",".join(sides["up"].keys())
                down_contigs = ",".join(sides["down"].keys())
                # Detail: top up + down hit
                up_detail = []
                for contig, hs in sides["up"].items():
                    for h in hs[:2]:
                        up_detail.append(f"{contig}:{h['ts']}-{h['te']}({h['strand']})")
                down_detail = []
                for contig, hs in sides["down"].items():
                    for h in hs[:2]:
                        down_detail.append(f"{contig}:{h['ts']}-{h['te']}({h['strand']})")
                detail = f"up=[{';'.join(up_detail)}] down=[{';'.join(down_detail)}]"
                f.write(f"{cat}\t{ref_id}\t{assembly}\t{up_contigs}\t{down_contigs}\t{detail}\n")
    print(f"\nWrote examples: {out_path}")


if __name__ == "__main__":
    main()
