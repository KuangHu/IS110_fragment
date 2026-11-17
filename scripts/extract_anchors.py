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
                        (distances_csv = the D values generated for that ref;
                        flanks are clamped at contig ends, so a ref near an
                        end keeps its small-D anchors)

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
import sys
from collections import defaultdict
from multiprocessing import Pool

import pysam


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
    p.add_argument("--threads", type=int, default=32,
                   help="Worker processes for parallel extraction (default 32)")
    return p.parse_args()


sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_alleles import revcomp  # noqa: E402  (IUPAC-aware)


# ---- multiprocessing worker (pysam, in-process FASTA reads) ----------------
_FA = None        # per-worker pysam.FastaFile
_PARAMS = None    # per-worker (flank, anchor_length, distances)


def _worker_init(genome_db, flank, anchor_length, distances):
    global _FA, _PARAMS
    _FA = pysam.FastaFile(genome_db)
    _PARAMS = (flank, anchor_length, distances)


def _extract_one(hit):
    """hit = (is_id, contig, tnp_start_1, tnp_end_1, strand).
    Returns (ref_record_str, anchors_str, table_row_str) or None (skipped)."""
    is_id, contig, tnp_start, tnp_end, strand = hit
    flank, anchor_length, distances = _PARAMS
    tnp_len = tnp_end - tnp_start + 1

    try:
        clen = _FA.get_reference_length(contig)
    except (KeyError, ValueError):
        return None

    # Clamp the flanks to the contig instead of dropping the hit: requiring
    # the full --flank on both sides discarded every IS within `flank` bp of a
    # contig end (4.5% of E. coli hits at flank=20 kb, far more on plasmids
    # and draft assemblies), even when the small-D anchors fit.
    left_avail = min(flank, tnp_start - 1)
    right_avail = min(flank, clen - tnp_end)
    if left_avail < 0 or right_avail < 0 or tnp_len <= 0:
        return None

    # pysam.fetch is 0-based half-open: [start, end)
    try:
        ref_seq = _FA.fetch(contig, tnp_start - 1 - left_avail, tnp_end + right_avail)
    except (KeyError, ValueError):
        return None
    if not ref_seq:
        return None
    ref_seq = ref_seq.upper()
    if len(ref_seq) != left_avail + tnp_len + right_avail:
        return None
    if strand == "-":
        ref_seq = revcomp(ref_seq)
        up_flank, down_flank = right_avail, left_avail
    else:
        up_flank, down_flank = left_avail, right_avail

    tnp_offset = up_flank
    tnp_end_in_ref = up_flank + tnp_len
    ref_len = len(ref_seq)

    anchor_parts = []
    used = []
    for d in distances:
        up_start = tnp_offset - d - anchor_length
        up_end = tnp_offset - d
        down_start = tnp_end_in_ref + d
        down_end = tnp_end_in_ref + d + anchor_length
        if up_start < 0 or down_end > ref_len:
            continue
        up_seq = ref_seq[up_start:up_end]
        down_seq = ref_seq[down_start:down_end]
        anchor_parts.append(f">{is_id}__up{d} d={d} anchor_len={anchor_length}\n{up_seq}")
        anchor_parts.append(f">{is_id}__down{d} d={d} anchor_len={anchor_length}\n{down_seq}")
        used.append(d)
    if not used:
        return None            # not even the smallest D fits on this contig
    anchors_str = "\n".join(anchor_parts) + "\n"

    ref_lines = [f">{is_id} ctg={contig} flank={up_flank} down_flank={down_flank} "
                 f"tnp_len={tnp_len} strand={strand}"]
    for i in range(0, ref_len, 80):
        ref_lines.append(ref_seq[i:i+80])
    ref_record = "\n".join(ref_lines) + "\n"

    # `distances` lists only the D values that were actually generated for
    # this reference, so downstream stages never expect anchors that do not
    # exist. flank_size is the requested flank (the ref may be shorter).
    table_row = (f"{is_id}\t{ref_len}\t{anchor_length}\t{tnp_len}\t"
                 f"{flank}\t{','.join(str(d) for d in used)}\n")
    return (ref_record, anchors_str, table_row, 2 * len(used))


def main():
    args = parse_args()
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    distances = [int(d) for d in args.distances.split(",")]
    unreachable = [d for d in distances if d + args.anchor_length > args.flank]
    if unreachable:
        # An anchor at distance D needs D + anchor_length bp of flank. These
        # used to be skipped without a word (default flank 80000 with
        # D=80000: the largest distance never existed).
        print(f"WARNING: --flank {args.flank} is too short for D="
              f"{','.join(map(str, unreachable))} (needs D + anchor length "
              f"{args.anchor_length}); no anchors will be made at these distances. "
              f"Raise --flank to >= {max(distances) + args.anchor_length}.",
              file=sys.stderr, flush=True)

    keep_ids = None
    if args.keep_ids:
        keep_ids = set()
        with open(args.keep_ids) as f:
            for line in f:
                if line.startswith(">"):
                    keep_ids.add(line[1:].split()[0])
        print(f"Restricting to {len(keep_ids)} is_id values from {args.keep_ids}",
              file=sys.stderr)

    # Parse all hits up front (cheap; just the TSV)
    hits = []
    n_input = 0
    with open(args.hits) as fhits:
        header = fhits.readline().rstrip().split("\t")
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
            if keep_ids is not None and is_id not in keep_ids:
                continue
            hits.append((is_id, parts[col["contig"]],
                         int(parts[col["tnp_start"]]),
                         int(parts[col["tnp_end"]]),
                         parts[col["tnp_strand"]]))

    print(f"Input hits: {n_input:,}; to process: {len(hits):,} "
          f"(threads={args.threads})", file=sys.stderr, flush=True)

    refs_out = open(args.out + "_refs.fa", "w")
    anchors_out = open(args.out + "_anchors.fa", "w")
    table_out = open(args.out + "_table.tsv", "w")
    table_out.write("is_id\tref_len\tanchor_len\ttnp_len\tflank_size\tdistances\n")

    n_refs = n_anchors = n_skipped = 0
    with Pool(args.threads, initializer=_worker_init,
              initargs=(args.genome_db, args.flank, args.anchor_length, distances)) as pool:
        for i, res in enumerate(pool.imap_unordered(_extract_one, hits, chunksize=16)):
            if i % 5000 == 0:
                print(f"  processed {i:,}/{len(hits):,}", file=sys.stderr, flush=True)
            if res is None:
                n_skipped += 1
                continue
            ref_record, anchors_str, table_row, na = res
            refs_out.write(ref_record)
            if anchors_str:
                anchors_out.write(anchors_str)
            table_out.write(table_row)
            n_refs += 1
            n_anchors += na

    refs_out.close()
    anchors_out.close()
    table_out.close()

    print(f"Input hits:        {n_input:,}", file=sys.stderr)
    print(f"References built:  {n_refs:,}", file=sys.stderr)
    print(f"Anchors generated: {n_anchors:,}", file=sys.stderr)
    print(f"Skipped (smallest D does not fit on the contig, or missing contig): {n_skipped:,}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
