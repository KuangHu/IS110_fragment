#!/usr/bin/env python3
"""Where else does each inserted fragment occur? (find_fragment_origin_local.sh)

One row per query fragment, also when it has no hit. A hit counts only if it
has >= --min-identity AND aligns >= --min-cov % of the fragment (a 300-bp
shared segment of a 2-kb fragment used to count like a full copy), and if it
is not the fragment's own locus (query names are <contig>__<start>_<end>).

What this table is NOT: the searches run with `minimap2 -N 5`, so a count
stops at 6 per organism (shown as "6+"); and "also present in organism X" is
presence, not origin. Organisms that tie for the best identity are all named.

Usage: fragment_origin_summary.py <query.fasta> <out.tsv> <org>=<paf> [<org>=<paf> ...]
"""
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paf_utils import aligned_identity, aligned_query_bases  # noqa: E402

CAP = 6     # primary + 5 secondaries


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("query_fasta")
    p.add_argument("out_tsv")
    p.add_argument("pafs", nargs="+", help="org=path")
    p.add_argument("--min-identity", type=float, default=90)
    p.add_argument("--min-cov", type=float, default=80)
    a = p.parse_args()

    queries = [l[1:].split()[0] for l in open(a.query_fasta) if l.startswith(">")]
    orgs, counts, best = [], {q: {} for q in queries}, {q: {} for q in queries}
    for spec in a.pafs:
        org, path = spec.split("=", 1)
        orgs.append(org)
        if not os.path.exists(path):
            print(f"  {org}: no PAF ({path}); reported as NA")
            for q in queries:
                counts[q][org] = None
            continue
        for q in queries:
            counts[q][org] = 0
        with open(path) as f:
            for line in f:
                c = line.rstrip("\n").split("\t")
                if len(c) < 12 or c[0] not in counts:
                    continue
                ident = aligned_identity(c)
                if ident < a.min_identity:
                    continue
                if 100.0 * aligned_query_bases(c) / max(1, int(c[1])) < a.min_cov:
                    continue
                m = re.match(r"^(.*)__(\d+)_(\d+)$", c[0])
                if m and c[5] == m.group(1):
                    s, e = int(m.group(2)), int(m.group(3))
                    if min(int(c[8]), e) - max(int(c[7]), s) > 0.5 * (e - s):
                        continue                      # the fragment itself
                counts[c[0]][org] += 1
                best[c[0]][org] = max(best[c[0]].get(org, 0.0), ident)

    def show(n):
        return "NA" if n is None else (f"{CAP}+" if n >= CAP else str(n))
    with open(a.out_tsv, "w") as out:
        out.write("query\t" + "\t".join(f"{o}_hits" for o in orgs)
                  + "\tbest_match_organism\tbest_identity\ttotal_hits\n")
        for q in queries:
            b = best[q]
            top = max(b.values()) if b else 0.0
            # ties are all named (the first organism read used to win)
            best_org = ",".join(o for o in orgs if o in b and abs(b[o] - top) < 0.05) or "none"
            total = sum(n for n in counts[q].values() if n)
            capped = any(n is not None and n >= CAP for n in counts[q].values())
            out.write(q + "\t" + "\t".join(show(counts[q][o]) for o in orgs)
                      + f"\t{best_org}\t{top:.1f}\t{total}{'+' if capped else ''}\n")
    n_hit = sum(1 for q in queries if best[q])
    print(f"Summary written: {a.out_tsv}")
    print(f"  {len(queries)} fragments; {n_hit} with a hit outside their own locus "
          f"(identity >= {a.min_identity:g}%, >= {a.min_cov:g}% of the fragment aligned); "
          f"{len(queries) - n_hit} without")
    print(f"  counts stop at {CAP} per organism (minimap2 -N 5)")


if __name__ == "__main__":
    main()
