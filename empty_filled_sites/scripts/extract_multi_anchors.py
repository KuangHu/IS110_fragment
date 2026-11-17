#!/usr/bin/env python3
"""Extract IS110 references with large flanks AND generate multi-distance anchors.

For each IS110 transposase position, extract:
  reference = [LARGE_FLANK upstream][transposase][LARGE_FLANK downstream]

Then sample anchors at multiple offsets from the transposase boundaries:
  D in [1000, 5000, 10000, 20000, 30000]
  up_anchor at distance D   = ref[flank_size - D - anchor_len : flank_size - D]
  down_anchor at distance D = ref[flank_size + tnp_len + D : flank_size + tnp_len + D + anchor_len]

Anchor naming: <ref_id>__up<D>  and  <ref_id>__down<D>

Outputs:
  - large_refs.fa   : full large-flank references (one record per IS110 position)
  - anchors.fa      : all anchors at all distances
  - ref_table.tsv   : ref_id, ref_len, anchor_len, tnp_len, flank_size, distances_csv

Usage:
    extract_multi_anchors.py --positions <is110_positions.tsv>
                             --genome-dir <NCBI_genome_dir>
                             --refs-fa <pilot_refs.fa>  (optional: filter to these IDs)
                             --out-dir <out_dir>
                             --flank 30000
                             --anchor-len 500
                             --distances 1000,5000,10000,20000,30000
"""

import argparse
import os
import sys
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--positions", required=True,
                   help="is110_positions.tsv (from HMM pipeline)")
    p.add_argument("--genome-dir", required=True,
                   help="Directory containing NCBI assembly subdirs with *_genomic.fna")
    p.add_argument("--refs-fa", default="",
                   help="If given, restrict to ref IDs (protein IDs) present in this FASTA")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--flank", type=int, default=30000)
    p.add_argument("--anchor-len", type=int, default=500)
    p.add_argument("--distances", default="1000,5000,10000,20000,30000",
                   help="Comma-separated distances from transposase edge")
    return p.parse_args()


def revcomp(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def load_fasta_index(fna):
    """Load a single FASTA into {contig: sequence}."""
    contigs = {}
    name = None
    parts = []
    with open(fna) as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    contigs[name] = "".join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.strip())
        if name:
            contigs[name] = "".join(parts)
    return contigs


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    distances = [int(d) for d in args.distances.split(",")]

    # Optional: restrict to a subset of ref IDs
    keep_ids = None
    if args.refs_fa:
        keep_ids = set()
        with open(args.refs_fa) as f:
            for line in f:
                if line.startswith(">"):
                    keep_ids.add(line[1:].split()[0])
        print(f"Restricting to {len(keep_ids)} ref IDs from {args.refs_fa}", file=sys.stderr)

    # Group positions by assembly
    by_asm = defaultdict(list)
    with open(args.positions) as f:
        next(f)
        for line in f:
            parts = line.rstrip().split("\t")
            if len(parts) < 7:
                continue
            protein_id, assembly, contig, start, end, strand, _ = parts
            if keep_ids and protein_id not in keep_ids:
                continue
            by_asm[assembly].append({
                "protein_id": protein_id,
                "contig": contig,
                "start": int(start),
                "end": int(end),
                "strand": strand,
            })

    print(f"Assemblies to process: {len(by_asm)}", file=sys.stderr)

    refs_out = open(os.path.join(args.out_dir, "large_refs.fa"), "w")
    anchors_out = open(os.path.join(args.out_dir, "anchors.fa"), "w")
    table_out = open(os.path.join(args.out_dir, "ref_table.tsv"), "w")
    table_out.write("ref_id\tref_len\tanchor_len\ttnp_len\tflank_size\tdistances\n")

    n_refs = 0
    n_anchors = 0
    n_skipped = 0

    for assembly, positions in by_asm.items():
        asm_dir = os.path.join(args.genome_dir, assembly)
        if not os.path.isdir(asm_dir):
            n_skipped += len(positions)
            continue
        fna = None
        for f in os.listdir(asm_dir):
            if f.endswith("_genomic.fna"):
                fna = os.path.join(asm_dir, f)
                break
        if not fna:
            n_skipped += len(positions)
            continue

        contigs = load_fasta_index(fna)

        for pos in positions:
            seq = contigs.get(pos["contig"])
            if not seq:
                n_skipped += 1
                continue
            tnp_start_0 = pos["start"] - 1  # GFF is 1-based inclusive
            tnp_end_0 = pos["end"]
            tnp_len = tnp_end_0 - tnp_start_0
            ref_start = tnp_start_0 - args.flank
            ref_end = tnp_end_0 + args.flank

            # Need full flank on both sides; otherwise skip
            if ref_start < 0 or ref_end > len(seq):
                n_skipped += 1
                continue

            ref_seq = seq[ref_start:ref_end]
            if pos["strand"] == "-":
                ref_seq = revcomp(ref_seq)

            ref_id = pos["protein_id"]
            ref_len = len(ref_seq)
            flank_size = args.flank
            tnp_offset = flank_size  # transposase starts at flank_size in ref_seq
            tnp_end_in_ref = flank_size + tnp_len

            # Write full reference
            refs_out.write(f">{ref_id} asm={assembly} ctg={pos['contig']} "
                           f"flank={flank_size} tnp_len={tnp_len}\n")
            for i in range(0, ref_len, 80):
                refs_out.write(ref_seq[i:i+80] + "\n")

            # Generate multi-distance anchors
            for d in distances:
                up_start = tnp_offset - d - args.anchor_len
                up_end = tnp_offset - d
                down_start = tnp_end_in_ref + d
                down_end = tnp_end_in_ref + d + args.anchor_len
                if up_start < 0 or down_end > ref_len:
                    continue
                up_seq = ref_seq[up_start:up_end]
                down_seq = ref_seq[down_start:down_end]
                anchors_out.write(f">{ref_id}__up{d} d={d} anchor_len={args.anchor_len}\n{up_seq}\n")
                anchors_out.write(f">{ref_id}__down{d} d={d} anchor_len={args.anchor_len}\n{down_seq}\n")
                n_anchors += 2

            table_out.write(f"{ref_id}\t{ref_len}\t{args.anchor_len}\t{tnp_len}\t{flank_size}\t"
                            f"{','.join(str(d) for d in distances)}\n")
            n_refs += 1

    refs_out.close()
    anchors_out.close()
    table_out.close()

    print(f"References extracted:  {n_refs}", file=sys.stderr)
    print(f"Anchors generated:     {n_anchors}", file=sys.stderr)
    print(f"Skipped (no flank room): {n_skipped}", file=sys.stderr)


if __name__ == "__main__":
    main()
