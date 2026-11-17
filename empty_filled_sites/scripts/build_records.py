#!/usr/bin/env python3
"""Build per-IS110 records: define boundaries, extract sequences, find best empty matches.

Pipeline:
  1. Parse PAF (anchors vs DB) to get all anchor pair observations + identities
  2. For each ref:
     - Find V0 peak (cluster of distances near 2*D + tnp_len)
     - Find empty peak (smaller cluster of distances)
     - IS110 length = V0_peak - empty_peak
  3. Extract from source NCBI genome:
     - transposase CDS
     - IS110 element (transposase + extension on both sides based on inferred boundary)
     - upstream flank (genomic, 5 kb default)
     - downstream flank (genomic, 5 kb default)
  4. For empty observations, compute up+down anchor mean identity
     Keep top-5 by similarity, extract their empty junction sequences
  5. Write single JSON file + supplementary FASTAs

Usage:
    build_records.py --paf <anchors_vs_db.paf>
                     --ref-table <ref_table.tsv>
                     --positions <is110_positions.tsv>
                     --genome-dir <NCBI_genome_dir>
                     --out-dir <out_dir>
                     [--anchor-distance 1000]
                     [--top-empty 5]
                     [--flank-out-len 5000]
"""

import argparse
import json
import os
import re
import sys
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--paf", required=True)
    p.add_argument("--ref-table", required=True)
    p.add_argument("--positions", required=True)
    p.add_argument("--genome-dir", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--anchor-distance", type=int, default=1000,
                   help="Which D to use for boundary calling (default: 1000)")
    p.add_argument("--top-empty", type=int, default=5,
                   help="Keep top-N empty observations by anchor identity (default: 5)")
    p.add_argument("--flank-out-len", type=int, default=5000,
                   help="Length of flanks to embed in JSON (default: 5000)")
    p.add_argument("--min-identity", type=float, default=95)
    p.add_argument("--min-cov", type=float, default=90)
    return p.parse_args()


ANCHOR_RE = re.compile(r"^(.+)__(up|down)(\d+)$")


def revcomp(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def load_fasta_index(fna):
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


def find_genome_fna(genome_dir, assembly):
    asm_dir = os.path.join(genome_dir, assembly)
    if not os.path.isdir(asm_dir):
        return None
    for f in os.listdir(asm_dir):
        if f.endswith("_genomic.fna"):
            return os.path.join(asm_dir, f)
    return None


def histogram_peaks(distances, bin_size=50):
    """Return list of (peak_center, peak_count) sorted by count desc."""
    if not distances:
        return []
    bins = defaultdict(int)
    for d in distances:
        bins[d // bin_size] += 1
    sorted_bins = sorted(bins.items(), key=lambda x: -x[1])
    peaks = []
    for k, count in sorted_bins:
        center = k * bin_size + bin_size // 2
        peaks.append((center, count))
    return peaks


def median(values):
    if not values:
        return None
    s = sorted(values)
    return s[len(s) // 2]


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # Load reference table
    ref_info = {}
    with open(args.ref_table) as f:
        next(f)
        for line in f:
            parts = line.rstrip().split("\t")
            ref_id, ref_len, anchor_len, tnp_len, flank_size, distances = parts
            ref_info[ref_id] = {
                "ref_len": int(ref_len),
                "anchor_len": int(anchor_len),
                "tnp_len": int(tnp_len),
                "flank_size": int(flank_size),
            }
    print(f"Loaded {len(ref_info)} refs from ref_table", file=sys.stderr)

    # Load IS110 positions (source genome coords)
    source_positions = {}
    with open(args.positions) as f:
        next(f)
        for line in f:
            parts = line.rstrip().split("\t")
            if len(parts) < 7:
                continue
            protein_id, assembly, contig, start, end, strand, plen = parts
            source_positions[protein_id] = {
                "assembly": assembly,
                "contig": contig,
                "start": int(start),
                "end": int(end),
                "strand": strand,
            }

    # Parse PAF: collect anchor hits per (ref, side, D, target_contig)
    # hits[(ref, side, D)][target_contig] = list of (strand, ts, te, ident)
    hits = defaultdict(lambda: defaultdict(list))
    D_target = args.anchor_distance

    n_total = 0
    n_kept = 0
    with open(args.paf) as f:
        for line in f:
            c = line.rstrip().split("\t")
            if len(c) < 12:
                continue
            n_total += 1
            qname = c[0]
            qlen = int(c[1])
            qs, qe = int(c[2]), int(c[3])
            strand = c[4]
            tname = c[5]
            tlen = int(c[6])
            ts, te = int(c[7]), int(c[8])
            matches, block = int(c[9]), int(c[10])
            ident = matches / block * 100 if block > 0 else 0
            cov = (qe - qs) / qlen * 100
            if ident < args.min_identity or cov < args.min_cov:
                continue
            n_kept += 1

            m = ANCHOR_RE.match(qname)
            if not m:
                continue
            ref_id, side, D = m.group(1), m.group(2), int(m.group(3))
            if D != D_target:
                continue
            hits[(ref_id, side)][tname].append({
                "strand": strand, "ts": ts, "te": te, "ident": ident, "tlen": tlen,
            })

    print(f"PAF: {n_total:,} rows, kept {n_kept:,} after filters", file=sys.stderr)

    # For each ref: compute pairs, identify peaks, build record
    all_records = []
    is_elem_seqs = []
    tnp_seqs = []
    up_flank_seqs = []
    down_flank_seqs = []
    empty_junction_seqs = []

    # Cache loaded genomes to avoid re-reading
    genome_cache = {}

    def get_genome(assembly):
        if assembly in genome_cache:
            return genome_cache[assembly]
        fna = find_genome_fna(args.genome_dir, assembly)
        if not fna:
            genome_cache[assembly] = None
            return None
        contigs = load_fasta_index(fna)
        genome_cache[assembly] = contigs
        return contigs

    for ref_id, info in ref_info.items():
        D = D_target
        up_hits = hits.get((ref_id, "up"), {})
        down_hits = hits.get((ref_id, "down"), {})
        if not up_hits or not down_hits:
            continue

        # Build list of all anchor pairs with metadata
        pairs = []
        for tname in set(up_hits.keys()) & set(down_hits.keys()):
            tlen = up_hits[tname][0]["tlen"]
            for up in up_hits[tname]:
                for down in down_hits[tname]:
                    if up["strand"] != down["strand"]:
                        continue
                    if up["strand"] == "+":
                        distance = down["ts"] - up["te"]
                    else:
                        distance = up["ts"] - down["te"]
                    if distance < -100 or distance > 200000:
                        continue
                    mean_ident = (up["ident"] + down["ident"]) / 2
                    # Parse assembly from contig name (format: ASSEMBLY|CHROM)
                    if "|" in tname:
                        assembly, contig = tname.split("|", 1)
                    else:
                        assembly, contig = "unknown", tname
                    pairs.append({
                        "assembly": assembly, "contig": contig, "tname": tname,
                        "strand": up["strand"], "distance": distance,
                        "up_ts": up["ts"], "up_te": up["te"],
                        "down_ts": down["ts"], "down_te": down["te"],
                        "up_ident": up["ident"], "down_ident": down["ident"],
                        "mean_ident": mean_ident,
                    })

        if not pairs:
            continue

        # Find peaks via histogram
        all_distances = [p["distance"] for p in pairs]
        peaks = histogram_peaks(all_distances, bin_size=50)
        if not peaks:
            continue

        # Top peak should be V0 (matches reference) ~= 2*D + tnp_len
        # Empty peak should be smaller, ~= 2*D
        expected_v0 = 2 * D + info["tnp_len"]
        expected_empty_min = 2 * D - info["tnp_len"]  # IS110 element could be tnp_len bigger than tnp ORF
        expected_empty_max = 2 * D + 200             # or up to a small overshoot

        # Find best V0 peak: closest to expected_v0
        v0_peak = None
        for center, count in peaks:
            if abs(center - expected_v0) < expected_v0 * 0.10 and count >= 2:
                v0_peak = (center, count)
                break

        # Find best empty peak: in the empty range, with at least 2 supporting pairs
        empty_peak = None
        for center, count in peaks:
            if expected_empty_min - 500 < center < expected_empty_max + 200 and count >= 2:
                if v0_peak and abs(center - v0_peak[0]) < 200:
                    continue  # this is the V0 peak, not empty
                empty_peak = (center, count)
                break

        # If we have both peaks, IS110 size = V0 - empty
        if v0_peak and empty_peak:
            is110_length = v0_peak[0] - empty_peak[0]
            # Distribute the extension symmetrically around transposase
            ext_per_side = (is110_length - info["tnp_len"]) // 2
            boundary_5p_offset = -ext_per_side
            boundary_3p_offset = is110_length - info["tnp_len"] - ext_per_side
        else:
            is110_length = info["tnp_len"]
            boundary_5p_offset = 0
            boundary_3p_offset = 0

        # Categorize each pair
        v0_obs = []
        empty_obs = []
        v1plus_obs = []
        for p in pairs:
            d = p["distance"]
            if v0_peak and abs(d - v0_peak[0]) <= 250:
                category = "v0"
                v0_obs.append(p)
            elif empty_peak and abs(d - empty_peak[0]) <= 250:
                category = "empty"
                empty_obs.append(p)
            elif v0_peak and d > v0_peak[0] + 500:
                category = "v1plus"
                v1plus_obs.append(p)

        # Sort empty observations by mean_ident desc, keep top-N
        empty_obs.sort(key=lambda x: -x["mean_ident"])
        top_empty = empty_obs[:args.top_empty]

        # Extract source sequences
        src = source_positions.get(ref_id)
        if not src:
            print(f"  WARNING: no source position for {ref_id}", file=sys.stderr)
            continue
        contigs = get_genome(src["assembly"])
        if not contigs:
            continue
        contig_seq = contigs.get(src["contig"])
        if not contig_seq:
            continue

        # Position handling: prodigal/HMM gives 1-based inclusive coords
        tnp_start_0 = src["start"] - 1
        tnp_end_0 = src["end"]

        # IS110 element bounds (in source genome coords)
        is_start_0 = tnp_start_0 + boundary_5p_offset
        is_end_0 = tnp_end_0 + boundary_3p_offset

        # Flanks (5 kb each)
        up_flank_start = max(0, is_start_0 - args.flank_out_len)
        up_flank_end = is_start_0
        down_flank_start = is_end_0
        down_flank_end = min(len(contig_seq), is_end_0 + args.flank_out_len)

        # Extract sequences (forward strand of contig)
        tnp_seq = contig_seq[tnp_start_0:tnp_end_0]
        is_seq = contig_seq[is_start_0:is_end_0]
        up_flank_seq = contig_seq[up_flank_start:up_flank_end]
        down_flank_seq = contig_seq[down_flank_start:down_flank_end]

        if src["strand"] == "-":
            tnp_seq = revcomp(tnp_seq)
            is_seq = revcomp(is_seq)
            # On - strand, "upstream" and "downstream" semantics flip
            up_flank_seq, down_flank_seq = revcomp(down_flank_seq), revcomp(up_flank_seq)

        # Extract empty junction sequences from top-5 empty obs
        for eo in top_empty:
            eo_contigs = get_genome(eo["assembly"])
            if not eo_contigs:
                continue
            eo_contig_seq = eo_contigs.get(eo["contig"]) or eo_contigs.get(eo["tname"])
            if not eo_contig_seq:
                continue
            if eo["strand"] == "+":
                junction_start = eo["up_te"]
                junction_end = eo["down_ts"]
                junction = eo_contig_seq[junction_start:junction_end]
            else:
                junction_start = eo["down_te"]
                junction_end = eo["up_ts"]
                junction = revcomp(eo_contig_seq[junction_start:junction_end])
            eo["junction_sequence"] = junction
            eo["junction_length"] = len(junction)

        # Build the record
        record = {
            "ref_id": ref_id,
            "source": {
                "assembly": src["assembly"],
                "contig": src["contig"],
                "transposase_start": src["start"],
                "transposase_end": src["end"],
                "transposase_strand": src["strand"],
            },
            "is_element": {
                "length": is_end_0 - is_start_0,
                "start_offset_5p": boundary_5p_offset,
                "end_offset_3p": boundary_3p_offset,
                "source_start": is_start_0 + 1,  # back to 1-based
                "source_end": is_end_0,
                "sequence": is_seq,
            },
            "transposase_cds": {
                "length": len(tnp_seq),
                "sequence": tnp_seq,
            },
            "upstream_flank": {
                "length": len(up_flank_seq),
                "sequence": up_flank_seq,
            },
            "downstream_flank": {
                "length": len(down_flank_seq),
                "sequence": down_flank_seq,
            },
            "boundary_evidence": {
                "method": "empty_vs_filled",
                "anchor_distance_D": D,
                "v0_peak_distance": v0_peak[0] if v0_peak else None,
                "v0_peak_count": v0_peak[1] if v0_peak else 0,
                "empty_peak_distance": empty_peak[0] if empty_peak else None,
                "empty_peak_count": empty_peak[1] if empty_peak else 0,
                "is_element_length_inferred": is110_length,
                "n_v0_observations": len(v0_obs),
                "n_empty_observations": len(empty_obs),
                "n_v1plus_observations": len(v1plus_obs),
                "n_top_empty_kept": len(top_empty),
            },
            "filled_observations": [
                {k: v for k, v in p.items() if k != "tname"} for p in v0_obs[:50]
            ],
            "empty_observations": [
                {k: v for k, v in p.items() if k != "tname"} for p in top_empty
            ],
            "v1plus_observations": [
                {k: v for k, v in p.items() if k != "tname"} for p in v1plus_obs[:50]
            ],
        }

        all_records.append(record)
        is_elem_seqs.append((ref_id, is_seq))
        tnp_seqs.append((ref_id, tnp_seq))
        up_flank_seqs.append((ref_id, up_flank_seq))
        down_flank_seqs.append((ref_id, down_flank_seq))
        for i, eo in enumerate(top_empty):
            if "junction_sequence" in eo:
                empty_junction_seqs.append(
                    (f"{ref_id}__empty{i}_{eo['assembly']}", eo["junction_sequence"])
                )

    # Write outputs
    json_path = os.path.join(args.out_dir, "records.json")
    with open(json_path, "w") as f:
        json.dump(all_records, f, indent=2)
    print(f"Wrote {len(all_records)} records to {json_path}", file=sys.stderr)

    def write_fa(path, seqs):
        with open(path, "w") as f:
            for name, seq in seqs:
                f.write(f">{name}\n")
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + "\n")

    write_fa(os.path.join(args.out_dir, "is_elements.fa"), is_elem_seqs)
    write_fa(os.path.join(args.out_dir, "transposase_cds.fa"), tnp_seqs)
    write_fa(os.path.join(args.out_dir, "upstream_flanks.fa"), up_flank_seqs)
    write_fa(os.path.join(args.out_dir, "downstream_flanks.fa"), down_flank_seqs)
    write_fa(os.path.join(args.out_dir, "empty_junctions.fa"), empty_junction_seqs)

    # Summary TSV
    with open(os.path.join(args.out_dir, "summary.tsv"), "w") as f:
        f.write("ref_id\tassembly\tcontig\ttnp_start\ttnp_end\ttnp_strand\t"
                "is_element_len\ttranposase_len\tv0_peak\tempty_peak\t"
                "n_v0\tn_empty\tn_v1plus\tn_top_empty\n")
        for r in all_records:
            be = r["boundary_evidence"]
            f.write(f"{r['ref_id']}\t{r['source']['assembly']}\t{r['source']['contig']}\t"
                    f"{r['source']['transposase_start']}\t{r['source']['transposase_end']}\t"
                    f"{r['source']['transposase_strand']}\t"
                    f"{r['is_element']['length']}\t{r['transposase_cds']['length']}\t"
                    f"{be['v0_peak_distance']}\t{be['empty_peak_distance']}\t"
                    f"{be['n_v0_observations']}\t{be['n_empty_observations']}\t"
                    f"{be['n_v1plus_observations']}\t{be['n_top_empty_kept']}\n")

    print(f"All outputs written to {args.out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
