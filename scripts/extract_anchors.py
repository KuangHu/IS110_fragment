#!/usr/bin/env python3
"""Extract IS-family reference regions and multi-distance anchors.

For each IS transposase position, extract:
  reference = [LARGE_FLANK upstream][transposase][LARGE_FLANK downstream]

Then sample anchors at multiple offsets from the transposase boundaries:
  D in {1000, 5000, 20000, 40000, 80000}  (default)
  up_anchor at distance D   = ref[flank - D - anchor_len : flank - D]
  down_anchor at distance D = ref[flank + tnp_len + D : flank + tnp_len + D + anchor_len]

Anchor naming: <is_id>__up<D>  and  <is_id>__down<D>

Outputs (one prefix --out):
  - <out>_refs.fa     : full large-flank references
  - <out>_anchors.fa  : all anchors at all distances
  - <out>_table.tsv   : is_id, ref_len, anchor_len, tnp_len, flank_size, distances_csv

Input: a single indexed FASTA (with .fai). Sequence extraction is via samtools.

Usage:
    extract_anchors.py --hits is_hits.tsv \\
                       --genome-db genome_db.fa \\
                       --out my_run/anchors \\
                       --flank 80000 \\
                       --anchor-length 500 \\
                       --distances 1000,5000,20000,40000,80000
"""
import argparse
import os
import subprocess
import sys
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--hits", required=True,
                   help="IS hits TSV (from is_detect.py): is_id, assembly, contig, "
                        "tnp_start, tnp_end, tnp_strand, tnp_len, ...")
    p.add_argument("--genome-db", required=True,
                   help="Indexed FASTA (.fai present) containing all contigs referenced "
                        "in --hits")
    p.add_argument("--out", required=True,
                   help="Output prefix. Produces <out>_refs.fa, <out>_anchors.fa, "
                        "<out>_table.tsv")
    p.add_argument("--flank", type=int, default=80000,
                   help="Bases of genomic flank to extract on each side (default 80000)")
    p.add_argument("--anchor-length", type=int, default=500,
                   help="Anchor length in bp (default 500)")
    p.add_argument("--distances", default="1000,5000,20000,40000,80000",
                   help="Comma-separated distances from transposase edge")
    p.add_argument("--keep-ids", default="",
                   help="Optional FASTA whose headers list is_id values to keep")
    return p.parse_args()


def revcomp(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def extract_region(db_fa, contig, start_1, end_1):
    """samtools faidx region (1-based inclusive). Returns uppercase string or None."""
    region = f"{contig}:{max(1, start_1)}-{end_1}"
    res = subprocess.run(["samtools", "faidx", db_fa, region],
                         capture_output=True, text=True)
    if res.returncode != 0:
        return None
    lines = res.stdout.split("\n")
    return "".join(lines[1:]).strip().upper()


def get_contig_length(db_fa, contig):
    """Read .fai for contig length."""
    fai = db_fa + ".fai"
    if not os.path.exists(fai):
        raise FileNotFoundError(f"FASTA index missing: {fai}. Run samtools faidx first.")
    with open(fai) as f:
        for line in f:
            parts = line.split("\t")
            if parts and parts[0] == contig:
                return int(parts[1])
    return None


def main():
    args = parse_args()
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    distances = [int(d) for d in args.distances.split(",")]

    keep_ids = None
    if args.keep_ids:
        keep_ids = set()
        with open(args.keep_ids) as f:
            for line in f:
                if line.startswith(">"):
                    keep_ids.add(line[1:].split()[0])
        print(f"Restricting to {len(keep_ids)} is_id values from {args.keep_ids}",
              file=sys.stderr)

    # Cache contig lengths
    contig_lens = {}

    refs_out = open(args.out + "_refs.fa", "w")
    anchors_out = open(args.out + "_anchors.fa", "w")
    table_out = open(args.out + "_table.tsv", "w")
    table_out.write("is_id\tref_len\tanchor_len\ttnp_len\tflank_size\tdistances\n")

    n_refs = 0
    n_anchors = 0
    n_skipped = 0
    n_input = 0

    with open(args.hits) as fhits:
        header = fhits.readline().rstrip().split("\t")
        # Build column index: tolerate is_id|assembly|contig|tnp_start|tnp_end|tnp_strand|tnp_len
        col = {name: i for i, name in enumerate(header)}
        for required in ("is_id", "contig", "tnp_start", "tnp_end", "tnp_strand"):
            if required not in col:
                sys.exit(f"--hits TSV missing column: {required}")

        for line in fhits:
            parts = line.rstrip().split("\t")
            if len(parts) < len(header):
                continue
            n_input += 1
            is_id = parts[col["is_id"]]
            if keep_ids and is_id not in keep_ids:
                continue
            contig = parts[col["contig"]]
            tnp_start = int(parts[col["tnp_start"]])  # 1-based inclusive
            tnp_end = int(parts[col["tnp_end"]])      # 1-based inclusive
            strand = parts[col["tnp_strand"]]
            tnp_len = tnp_end - tnp_start + 1

            if contig not in contig_lens:
                cl = get_contig_length(args.genome_db, contig)
                if cl is None:
                    n_skipped += 1
                    continue
                contig_lens[contig] = cl
            clen = contig_lens[contig]

            ref_start_1 = tnp_start - args.flank
            ref_end_1 = tnp_end + args.flank
            if ref_start_1 < 1 or ref_end_1 > clen:
                n_skipped += 1
                continue

            ref_seq = extract_region(args.genome_db, contig, ref_start_1, ref_end_1)
            if not ref_seq or len(ref_seq) != args.flank + tnp_len + args.flank:
                n_skipped += 1
                continue
            if strand == "-":
                ref_seq = revcomp(ref_seq)

            flank_size = args.flank
            tnp_offset = flank_size
            tnp_end_in_ref = flank_size + tnp_len
            ref_len = len(ref_seq)

            refs_out.write(f">{is_id} ctg={contig} flank={flank_size} "
                           f"tnp_len={tnp_len} strand={strand}\n")
            for i in range(0, ref_len, 80):
                refs_out.write(ref_seq[i:i+80] + "\n")

            for d in distances:
                up_start = tnp_offset - d - args.anchor_length
                up_end = tnp_offset - d
                down_start = tnp_end_in_ref + d
                down_end = tnp_end_in_ref + d + args.anchor_length
                if up_start < 0 or down_end > ref_len:
                    continue
                up_seq = ref_seq[up_start:up_end]
                down_seq = ref_seq[down_start:down_end]
                anchors_out.write(f">{is_id}__up{d} d={d} "
                                  f"anchor_len={args.anchor_length}\n{up_seq}\n")
                anchors_out.write(f">{is_id}__down{d} d={d} "
                                  f"anchor_len={args.anchor_length}\n{down_seq}\n")
                n_anchors += 2

            table_out.write(f"{is_id}\t{ref_len}\t{args.anchor_length}\t{tnp_len}\t"
                            f"{flank_size}\t{','.join(str(d) for d in distances)}\n")
            n_refs += 1

    refs_out.close()
    anchors_out.close()
    table_out.close()

    print(f"Input hits:        {n_input:,}", file=sys.stderr)
    print(f"References built:  {n_refs:,}", file=sys.stderr)
    print(f"Anchors generated: {n_anchors:,}", file=sys.stderr)
    print(f"Skipped (insufficient flank or missing contig): {n_skipped:,}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
