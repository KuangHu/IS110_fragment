#!/usr/bin/env python3
"""Build per-IS110 records using multi-distance anchor data.

Algorithm per ref:
  1. Parse PAF(s): collect anchor pair observations at all D values
  2. For each D, find pairs and classify into v0_filled / empty / etc.
  3. Pick BEST D = smallest D where we found at least 1 empty AND at least 1 V0
     (this gives most direct evidence; smaller D = more genome-specific anchors)
  4. If no D gives both, fall back to D with the most empties (use expected V0)
  5. Compute IS110 size = v0_distance - empty_distance
  6. Extract source sequences using the inferred boundaries
  7. Keep top-N empty observations by mean anchor identity
  8. If NCBI gave no empties, also include Logan empties in the record

Usage:
    build_records_multi_d.py --paf <ncbi.paf> --logan-paf <logan.paf>
                             --ref-table <ref_table.tsv>
                             --positions <is110_positions.tsv>
                             --genome-dir <NCBI_genome_dir>
                             --logan-db <logan.fa>
                             --out-dir <out_dir>
"""

import argparse
import json
import os
import re
import subprocess
import sys
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--paf", required=True, help="NCBI PAF")
    p.add_argument("--logan-paf", default="", help="Logan PAF (optional)")
    p.add_argument("--ref-table", required=True)
    p.add_argument("--positions", required=True)
    p.add_argument("--genome-dir", required=True)
    p.add_argument("--logan-db", default="")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--top-empty", type=int, default=5)
    p.add_argument("--flank-out-len", type=int, default=5000)
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


def parse_paf_into_hits(paf_path, min_identity, min_cov):
    """Return dict: (ref_id, side, D) -> {tname: [hit_dict, ...]}"""
    hits = defaultdict(lambda: defaultdict(list))
    n_total = n_kept = 0
    with open(paf_path) as f:
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
            if ident < min_identity or cov < min_cov:
                continue
            n_kept += 1
            m = ANCHOR_RE.match(qname)
            if not m:
                continue
            ref_id, side, D = m.group(1), m.group(2), int(m.group(3))
            hits[(ref_id, side, D)][tname].append({
                "strand": strand, "ts": ts, "te": te, "ident": ident, "tlen": tlen,
            })
    return hits, n_total, n_kept


def histogram_peaks(distances, bin_size=50):
    if not distances:
        return []
    bins = defaultdict(int)
    for d in distances:
        bins[d // bin_size] += 1
    sorted_bins = sorted(bins.items(), key=lambda x: -x[1])
    return [(k * bin_size + bin_size // 2, count) for k, count in sorted_bins]


def pair_anchors(hits, ref_id, D, info, max_pair_distance=200000):
    """Build anchor pairs and classify by distance."""
    up_hits = hits.get((ref_id, "up", D), {})
    down_hits = hits.get((ref_id, "down", D), {})
    pairs = []
    expected_v0 = 2 * D + info["tnp_len"]

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
                if distance < -100 or distance > max_pair_distance:
                    continue
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
                    "mean_ident": (up["ident"] + down["ident"]) / 2,
                    "D": D,
                })
    return pairs, expected_v0


def categorize_pairs(pairs, expected_v0, D):
    """Identify V0 peak and empty peak from pair distances."""
    distances = [p["distance"] for p in pairs]
    peaks = histogram_peaks(distances, bin_size=50)
    if not peaks:
        return None, None, [], [], []

    # V0 peak: closest to expected_v0 (within 10%)
    v0_peak = None
    for center, count in peaks:
        if abs(center - expected_v0) <= max(200, expected_v0 * 0.10) and count >= 1:
            v0_peak = (center, count)
            break

    # Empty peak: distance < expected_v0 - tnp_len/2 (at least tnp_len/2 shorter)
    # Practical: empty range is roughly [2*D - some_extension, expected_v0 - tnp_len + extension]
    # Use any peak with center < expected_v0 - 500 and count >= 1
    empty_peak = None
    for center, count in peaks:
        if v0_peak and abs(center - v0_peak[0]) < 100:
            continue
        if center < expected_v0 - 500 and count >= 1:
            empty_peak = (center, count)
            break

    # Categorize each pair
    v0_obs = []
    empty_obs = []
    v1plus_obs = []
    for p in pairs:
        d = p["distance"]
        if v0_peak and abs(d - v0_peak[0]) <= 250:
            v0_obs.append(p)
        elif empty_peak and abs(d - empty_peak[0]) <= 250:
            empty_obs.append(p)
        elif v0_peak and d > v0_peak[0] + 500:
            v1plus_obs.append(p)

    return v0_peak, empty_peak, v0_obs, empty_obs, v1plus_obs


def pick_best_d(hits, ref_id, info, distances):
    """Find the best D for this ref: prefer smallest D with both V0 and empty."""
    candidates = []
    for D in distances:
        pairs, expected_v0 = pair_anchors(hits, ref_id, D, info)
        if not pairs:
            continue
        v0_peak, empty_peak, v0_obs, empty_obs, v1plus_obs = categorize_pairs(pairs, expected_v0, D)
        candidates.append({
            "D": D, "pairs": pairs, "expected_v0": expected_v0,
            "v0_peak": v0_peak, "empty_peak": empty_peak,
            "v0_obs": v0_obs, "empty_obs": empty_obs, "v1plus_obs": v1plus_obs,
        })

    if not candidates:
        return None

    # Prefer smallest D with both V0 and empty (highest-quality boundary call)
    for c in sorted(candidates, key=lambda x: x["D"]):
        if c["v0_peak"] and c["empty_peak"]:
            return c

    # Else: any D with empty
    for c in sorted(candidates, key=lambda x: x["D"]):
        if c["empty_peak"]:
            return c

    # Else: just pick the first one (only V0)
    return candidates[0]


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # Reference table
    ref_info = {}
    distances_set = set()
    with open(args.ref_table) as f:
        next(f)
        for line in f:
            parts = line.rstrip().split("\t")
            ref_id, ref_len, anchor_len, tnp_len, flank_size, distances_str = parts
            ds = [int(d) for d in distances_str.split(",")]
            distances_set.update(ds)
            ref_info[ref_id] = {
                "ref_len": int(ref_len), "anchor_len": int(anchor_len),
                "tnp_len": int(tnp_len), "flank_size": int(flank_size),
                "distances": ds,
            }
    distances = sorted(distances_set)

    # Source positions
    source_positions = {}
    with open(args.positions) as f:
        next(f)
        for line in f:
            parts = line.rstrip().split("\t")
            if len(parts) < 7:
                continue
            protein_id, assembly, contig, start, end, strand, _ = parts
            source_positions[protein_id] = {
                "assembly": assembly, "contig": contig,
                "start": int(start), "end": int(end), "strand": strand,
            }

    # Parse NCBI PAF
    print(f"Parsing NCBI PAF...", file=sys.stderr)
    ncbi_hits, n_total, n_kept = parse_paf_into_hits(args.paf, args.min_identity, args.min_cov)
    print(f"  NCBI: {n_total:,} rows, {n_kept:,} kept", file=sys.stderr)

    # Parse Logan PAF if provided
    logan_hits = {}
    if args.logan_paf and os.path.exists(args.logan_paf):
        print(f"Parsing Logan PAF...", file=sys.stderr)
        logan_hits, n_l_total, n_l_kept = parse_paf_into_hits(args.logan_paf, args.min_identity, args.min_cov)
        print(f"  Logan: {n_l_total:,} rows, {n_l_kept:,} kept", file=sys.stderr)

    # Cache loaded NCBI genomes
    genome_cache = {}
    def get_ncbi_genome(assembly):
        if assembly in genome_cache:
            return genome_cache[assembly]
        fna = find_genome_fna(args.genome_dir, assembly)
        if not fna:
            genome_cache[assembly] = None
            return None
        contigs = load_fasta_index(fna)
        genome_cache[assembly] = contigs
        return contigs

    # Logan DB lookup via samtools faidx
    def get_logan_seq(tname, ts, te):
        if not args.logan_db or not os.path.exists(args.logan_db):
            return None
        cmd = ["samtools", "faidx", args.logan_db, f"{tname}:{max(1, ts+1)}-{te}"]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            return None
        seq = "".join(line.strip() for line in result.stdout.split("\n")
                      if not line.startswith(">") and line)
        return seq

    # Build records
    all_records = []
    is_elem_seqs = []
    tnp_seqs = []
    up_flank_seqs = []
    down_flank_seqs = []
    empty_junction_seqs = []

    n_ncbi = n_logan = n_neither = 0

    for ref_id, info in ref_info.items():
        # Try NCBI first
        best = pick_best_d(ncbi_hits, ref_id, info, distances)
        source_db = "ncbi"
        if not best or not best["empty_obs"]:
            # Fall back to Logan
            if logan_hits:
                best_logan = pick_best_d(logan_hits, ref_id, info, distances)
                if best_logan and best_logan["empty_obs"]:
                    best = best_logan
                    source_db = "logan"

        if not best:
            n_neither += 1
            continue

        if source_db == "ncbi":
            n_ncbi += 1
        else:
            n_logan += 1

        # Compute IS110 size
        if best["v0_peak"] and best["empty_peak"]:
            is110_length = best["v0_peak"][0] - best["empty_peak"][0]
        elif best["empty_peak"]:
            is110_length = best["expected_v0"] - best["empty_peak"][0]
        else:
            is110_length = info["tnp_len"]

        ext_per_side = (is110_length - info["tnp_len"]) // 2
        boundary_5p_offset = -ext_per_side
        boundary_3p_offset = is110_length - info["tnp_len"] - ext_per_side

        # Extract source sequences
        src = source_positions.get(ref_id)
        if not src:
            continue
        contigs = get_ncbi_genome(src["assembly"])
        if not contigs:
            continue
        contig_seq = contigs.get(src["contig"])
        if not contig_seq:
            continue

        tnp_start_0 = src["start"] - 1
        tnp_end_0 = src["end"]
        is_start_0 = max(0, tnp_start_0 + boundary_5p_offset)
        is_end_0 = min(len(contig_seq), tnp_end_0 + boundary_3p_offset)

        up_flank_start = max(0, is_start_0 - args.flank_out_len)
        up_flank_end = is_start_0
        down_flank_start = is_end_0
        down_flank_end = min(len(contig_seq), is_end_0 + args.flank_out_len)

        tnp_seq = contig_seq[tnp_start_0:tnp_end_0]
        is_seq = contig_seq[is_start_0:is_end_0]
        up_flank_seq = contig_seq[up_flank_start:up_flank_end]
        down_flank_seq = contig_seq[down_flank_start:down_flank_end]

        if src["strand"] == "-":
            tnp_seq = revcomp(tnp_seq)
            is_seq = revcomp(is_seq)
            up_flank_seq, down_flank_seq = revcomp(down_flank_seq), revcomp(up_flank_seq)

        # Top-N empty observations by mean identity
        empty_obs_sorted = sorted(best["empty_obs"], key=lambda x: -x["mean_ident"])
        top_empty = empty_obs_sorted[:args.top_empty]

        # Extract empty junction sequences
        for eo in top_empty:
            junction = None
            if source_db == "ncbi":
                eo_contigs = get_ncbi_genome(eo["assembly"])
                if eo_contigs:
                    eo_seq = eo_contigs.get(eo["contig"]) or eo_contigs.get(eo["tname"])
                    if eo_seq:
                        if eo["strand"] == "+":
                            junction = eo_seq[eo["up_te"]:eo["down_ts"]]
                        else:
                            junction = revcomp(eo_seq[eo["down_te"]:eo["up_ts"]])
            else:
                if eo["strand"] == "+":
                    junction = get_logan_seq(eo["tname"], eo["up_te"], eo["down_ts"])
                else:
                    raw = get_logan_seq(eo["tname"], eo["down_te"], eo["up_ts"])
                    junction = revcomp(raw) if raw else None
            if junction:
                eo["junction_sequence"] = junction
                eo["junction_length"] = len(junction)

        record = {
            "ref_id": ref_id,
            "source": {
                "assembly": src["assembly"], "contig": src["contig"],
                "transposase_start": src["start"], "transposase_end": src["end"],
                "transposase_strand": src["strand"],
            },
            "is_element": {
                "length": is_end_0 - is_start_0,
                "start_offset_5p": boundary_5p_offset,
                "end_offset_3p": boundary_3p_offset,
                "source_start": is_start_0 + 1, "source_end": is_end_0,
                "sequence": is_seq,
            },
            "transposase_cds": {"length": len(tnp_seq), "sequence": tnp_seq},
            "upstream_flank": {"length": len(up_flank_seq), "sequence": up_flank_seq},
            "downstream_flank": {"length": len(down_flank_seq), "sequence": down_flank_seq},
            "boundary_evidence": {
                "method": "empty_vs_filled",
                "source_db": source_db,
                "anchor_distance_D": best["D"],
                "v0_peak_distance": best["v0_peak"][0] if best["v0_peak"] else None,
                "v0_peak_count": best["v0_peak"][1] if best["v0_peak"] else 0,
                "empty_peak_distance": best["empty_peak"][0] if best["empty_peak"] else None,
                "empty_peak_count": best["empty_peak"][1] if best["empty_peak"] else 0,
                "is_element_length_inferred": is110_length,
                "n_v0_observations": len(best["v0_obs"]),
                "n_empty_observations": len(best["empty_obs"]),
                "n_v1plus_observations": len(best["v1plus_obs"]),
                "n_top_empty_kept": len(top_empty),
            },
            "filled_observations": [
                {k: v for k, v in p.items() if k != "tname"} for p in best["v0_obs"][:50]
            ],
            "empty_observations": [
                {k: v for k, v in p.items() if k != "tname"} for p in top_empty
            ],
            "v1plus_observations": [
                {k: v for k, v in p.items() if k != "tname"} for p in best["v1plus_obs"][:50]
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
                "is_element_len\ttnp_len\tD_used\tsource_db\tv0_peak\tempty_peak\t"
                "n_v0\tn_empty\tn_v1plus\tn_top_empty\n")
        for r in all_records:
            be = r["boundary_evidence"]
            f.write(f"{r['ref_id']}\t{r['source']['assembly']}\t{r['source']['contig']}\t"
                    f"{r['source']['transposase_start']}\t{r['source']['transposase_end']}\t"
                    f"{r['source']['transposase_strand']}\t"
                    f"{r['is_element']['length']}\t{r['transposase_cds']['length']}\t"
                    f"{be['anchor_distance_D']}\t{be['source_db']}\t"
                    f"{be['v0_peak_distance']}\t{be['empty_peak_distance']}\t"
                    f"{be['n_v0_observations']}\t{be['n_empty_observations']}\t"
                    f"{be['n_v1plus_observations']}\t{be['n_top_empty_kept']}\n")

    print(f"\n=== Build summary ===", file=sys.stderr)
    print(f"  Used NCBI evidence:  {n_ncbi}", file=sys.stderr)
    print(f"  Used Logan evidence: {n_logan}", file=sys.stderr)
    print(f"  No empty found:      {n_neither}", file=sys.stderr)
    print(f"  All outputs in: {args.out_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
