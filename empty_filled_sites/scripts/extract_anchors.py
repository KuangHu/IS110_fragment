#!/usr/bin/env python3
"""Extract anchor sequences from IS110 references for empty/filled site analysis.

For each IS110 reference (5.4 kb = [2kb flank][IS110 ~1.4kb][2kb flank]):
  - up_anchor   = first ANCHOR_LEN bp (outermost upstream flank)
  - down_anchor = last ANCHOR_LEN bp (outermost downstream flank)

Output one FASTA with both anchors per reference, named:
  <ref_id>__up
  <ref_id>__down

Usage:
    extract_anchors.py <refs.fa> <output_anchors.fa> [--anchor-len 500] [--top-n N]
"""

import argparse
import sys


def parse_fasta(path):
    """Yield (name, header, seq) tuples."""
    name = None
    header = None
    parts = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    yield name, header, "".join(parts)
                header = line[1:].rstrip()
                name = header.split()[0]
                parts = []
            else:
                parts.append(line.strip())
    if name is not None:
        yield name, header, "".join(parts)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("refs_fa", help="Input IS110 reference FASTA (with flanks)")
    p.add_argument("out_fa", help="Output anchor FASTA")
    p.add_argument("--anchor-len", type=int, default=500,
                   help="Length of each anchor (default: 500)")
    p.add_argument("--min-ref-len", type=int, default=4000,
                   help="Skip refs shorter than this (default: 4000)")
    p.add_argument("--top-n", type=int, default=0,
                   help="Only process first N refs (0 = all; default: 0)")
    args = p.parse_args()

    n_written = 0
    n_skipped = 0
    with open(args.out_fa, "w") as out:
        for i, (name, header, seq) in enumerate(parse_fasta(args.refs_fa)):
            if args.top_n > 0 and i >= args.top_n:
                break
            if len(seq) < args.min_ref_len:
                n_skipped += 1
                continue
            up = seq[:args.anchor_len]
            down = seq[-args.anchor_len:]
            out.write(f">{name}__up reflen={len(seq)} anchorlen={args.anchor_len}\n{up}\n")
            out.write(f">{name}__down reflen={len(seq)} anchorlen={args.anchor_len}\n{down}\n")
            n_written += 2

    print(f"Wrote {n_written} anchor sequences to {args.out_fa}", file=sys.stderr)
    print(f"Skipped {n_skipped} short references", file=sys.stderr)


if __name__ == "__main__":
    main()
