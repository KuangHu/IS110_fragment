#!/usr/bin/env python3
"""Stage 3: align the anchors against the genome database.

    search_anchors.py --db genomes.fa --anchors anchors.fa --out anchors_vs_db.paf \
                      --max-secondary 10000 --threads 32
    search_anchors.py --from-paf raw.paf --out anchors_vs_db.paf      # convert an old PAF

Writes two files:

  <out>            the alignments, first 12 PAF columns, SORTED so that all
                   rows of one element are together. Stages 4-6 read it one
                   element at a time (iter_anchor_groups); the raw PAF of a
                   multi-part index is in index-part order and could only be
                   used by loading it whole (98 GB of memory for the 21.6-GB
                   E. coli PAF at -N 50).
  <out>.rows.tsv   per anchor: rows, index parts it has rows in, and the
                   largest number of rows in ONE part. minimap2 applies -N per
                   index part, so only the last number can show that the cap
                   was reached; after sorting it is no longer recoverable from
                   the PAF, which is why it is recorded here.

The sort is stable, so the order of rows within an anchor is minimap2's.
"""
import argparse
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_alleles import ROWS_SUFFIX, _RunCounter  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--db", help="Target genome FASTA (or .mmi)")
    p.add_argument("--anchors", help="Anchor FASTA from extract_anchors.py")
    p.add_argument("--from-paf", help="Sort an existing raw PAF instead of searching")
    p.add_argument("--out", required=True, help="Output PAF (sorted, 12 columns)")
    p.add_argument("--max-secondary", type=int, default=50,
                   help="minimap2 -N: secondary alignments per anchor and index part")
    p.add_argument("--secondary-ratio", type=float, default=0.5, help="minimap2 -p")
    p.add_argument("--index-batch", default=None,
                   help="minimap2 -I (e.g. 16G); default: minimap2's own (8G)")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--sort-mem", default="8G", help="sort -S")
    p.add_argument("--tmp-dir", default=None,
                   help="Directory for sort temporaries (default: next to --out)")
    return p.parse_args()


def main():
    args = parse_args()
    if bool(args.from_paf) == bool(args.db and args.anchors):
        sys.exit("ERROR: give either --db and --anchors, or --from-paf")
    out_dir = os.path.dirname(os.path.abspath(args.out))
    tmp_dir = args.tmp_dir or out_dir
    tmp_out = args.out + ".sorting"

    env = dict(os.environ, LC_ALL="C")
    sorter = subprocess.Popen(
        ["sort", "-s", "-t", "\t", "-k1,1", "-S", args.sort_mem,
         "--parallel", str(max(1, min(args.threads, 16))), "-T", tmp_dir,
         "-o", tmp_out],
        stdin=subprocess.PIPE, env=env)

    mm2 = None
    if args.from_paf:
        src = open(args.from_paf, "rb")
    else:
        cmd = ["minimap2", "-x", "sr", "-c", "--eqx", "--secondary=yes",
               "-N", str(args.max_secondary), "-p", str(args.secondary_ratio),
               "-t", str(args.threads)]
        if args.index_batch:
            cmd += ["-I", args.index_batch]
        cmd += [args.db, args.anchors]
        print(" ".join(cmd), file=sys.stderr, flush=True)
        mm2 = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        src = mm2.stdout

    counter = _RunCounter()
    n = 0
    w = sorter.stdin
    for line in src:
        c = line.split(b"\t", 12)
        if len(c) < 12:
            continue
        counter.feed(c[0])
        if len(c) > 12:
            w.write(b"\t".join(c[:12]) + b"\n")
        else:
            w.write(line if line.endswith(b"\n") else line + b"\n")
        n += 1
    counter.flush()
    src.close()
    w.close()
    if mm2 is not None and mm2.wait() != 0:
        sorter.kill()
        sys.exit(f"ERROR: minimap2 exited {mm2.returncode}")
    if sorter.wait() != 0:
        sys.exit(f"ERROR: sort exited {sorter.returncode}")

    rows_tmp = args.out + ROWS_SUFFIX + ".tmp"
    with open(rows_tmp, "w") as f:
        f.write("anchor\trows\tindex_parts\tmax_rows_one_part\n")
        for q in sorted(counter.total):
            f.write(f"{q.decode()}\t{counter.total[q]}\t{counter.n_runs[q]}\t"
                    f"{counter.max_run[q]}\n")
    # sidecar first: a PAF without it would be read as "runs unknown"
    os.replace(rows_tmp, args.out + ROWS_SUFFIX)
    os.replace(tmp_out, args.out)

    cap = args.max_secondary + 1
    if not args.from_paf:
        at_cap = sum(1 for v in counter.max_run.values() if v >= cap)
        print(f"anchors with alignments: {len(counter.total):,}; at the -N cap "
              f"({cap} rows in one index part): {at_cap:,} "
              f"({100 * at_cap / max(1, len(counter.total)):.1f}%)", file=sys.stderr)
    print(f"wrote {n:,} rows -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
