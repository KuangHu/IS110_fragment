#!/usr/bin/env python3
"""Validate inversion/translocation/duplication calls using nucmer whole-genome alignment.

For each top example:
  1. Find source genome FASTA (where the IS110 reference came from)
  2. Find target genome FASTA (where rearrangement was detected)
  3. Extract relevant chromosomal regions (around the anchor positions)
  4. Run nucmer on those regions
  5. Parse delta output and confirm rearrangement pattern

For each example, output a verdict:
  - CONFIRMED: nucmer alignment shows the rearrangement
  - REJECTED: nucmer shows simple co-linear alignment (false positive)
  - AMBIGUOUS: alignment is complex / inconclusive

Usage:
    validate_rearrangements.py <examples_tsv> <genome_dir> <out_dir> [--n-per-category N]
"""

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("examples_tsv")
    p.add_argument("genome_dir")
    p.add_argument("out_dir")
    p.add_argument("--n-per-category", type=int, default=5)
    p.add_argument("--region-flank", type=int, default=50000,
                   help="bp around anchor positions to extract for nucmer")
    return p.parse_args()


def find_genome_fna(genome_dir, assembly):
    asm_dir = os.path.join(genome_dir, assembly)
    if not os.path.isdir(asm_dir):
        return None
    for f in os.listdir(asm_dir):
        if f.endswith("_genomic.fna"):
            return os.path.join(asm_dir, f)
    return None


def parse_detail(detail_str):
    """Parse 'up=[contig:start-end(strand);...] down=[...]' into structured form."""
    result = {"up": [], "down": []}
    for side in ("up", "down"):
        m = re.search(rf"{side}=\[([^\]]*)\]", detail_str)
        if not m:
            continue
        for hit in m.group(1).split(";"):
            if not hit:
                continue
            mm = re.match(r"(.+):(\d+)-(\d+)\(([+-])\)", hit)
            if mm:
                result[side].append({
                    "contig": mm.group(1),
                    "start": int(mm.group(2)),
                    "end": int(mm.group(3)),
                    "strand": mm.group(4),
                })
    return result


def extract_region(fna, contig, start, end, out_path):
    """Extract a single region using samtools faidx."""
    if not os.path.exists(fna + ".fai"):
        subprocess.run(["samtools", "faidx", fna], check=True)
    region = f"{contig}:{max(1, start)}-{end}"
    result = subprocess.run(
        ["samtools", "faidx", fna, region],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        return False
    with open(out_path, "w") as f:
        f.write(result.stdout)
    return True


def run_nucmer(ref_fa, query_fa, out_prefix):
    """Run nucmer + show-coords. Returns parsed alignment blocks."""
    subprocess.run(
        ["nucmer", "--maxmatch", "-p", out_prefix, ref_fa, query_fa],
        check=False, capture_output=True
    )
    delta = out_prefix + ".delta"
    if not os.path.exists(delta):
        return []
    result = subprocess.run(
        ["show-coords", "-r", "-c", "-l", "-T", delta],
        capture_output=True, text=True
    )
    blocks = []
    for line in result.stdout.split("\n"):
        if not line or line.startswith("=") or line.startswith("/") or line.startswith("NUCMER") or line.startswith("[S1]"):
            continue
        f = line.split("\t")
        if len(f) < 9:
            continue
        try:
            blocks.append({
                "rs": int(f[0]), "re": int(f[1]),
                "qs": int(f[2]), "qe": int(f[3]),
                "ralen": int(f[4]), "qalen": int(f[5]),
                "ident": float(f[6]),
                "rlen": int(f[7]), "qlen": int(f[8]),
                "ref_name": f[-2] if len(f) >= 11 else "",
                "qry_name": f[-1] if len(f) >= 11 else "",
            })
        except (ValueError, IndexError):
            continue
    return blocks


def classify_alignment(blocks, expected):
    """Classify nucmer alignment: confirms inversion/translocation/duplication?"""
    if not blocks:
        return "no_alignment", []

    # Inversion signature: alignment block where (qe < qs) -- minus strand on query
    # i.e. ref forward maps to query reverse
    inverted = [b for b in blocks if b["qs"] > b["qe"]]
    forward = [b for b in blocks if b["qs"] < b["qe"]]

    if expected == "inversion_only":
        if inverted and len(inverted) >= 1:
            inv_total = sum(b["ralen"] for b in inverted)
            fwd_total = sum(b["ralen"] for b in forward)
            if inv_total > fwd_total * 0.5:
                return "CONFIRMED_inversion", blocks
        return "REJECTED_or_simple", blocks

    elif expected == "translocation_only":
        # Translocation: blocks come from different ref or query ranges,
        # but here we extract regions specifically — use number of distinct
        # orientation blocks as proxy
        return "INSPECT_manually", blocks

    elif expected == "duplication":
        # Duplication: same query region maps to multiple ref regions
        if len(forward) + len(inverted) >= 2:
            return "CONFIRMED_duplication", blocks
        return "REJECTED_single_alignment", blocks

    return "UNKNOWN", blocks


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    examples = defaultdict(list)
    with open(args.examples_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            examples[row["category"]].append(row)

    report_path = os.path.join(args.out_dir, "validation_report.tsv")
    out = open(report_path, "w")
    out.write("category\tref_id\tassembly\tup_loc\tdown_loc\tverdict\tnucmer_summary\tdotplot_path\n")

    categories = ["inversion_only", "translocation_only", "duplication"]

    for category in categories:
        eg_list = examples.get(category, [])
        if not eg_list:
            continue
        print(f"\n=== Validating {category} ({len(eg_list)} candidates, processing top {args.n_per_category}) ===",
              file=sys.stderr)

        for i, eg in enumerate(eg_list[:args.n_per_category]):
            ref_id = eg["ref_id"]
            target_assembly = eg["assembly"]

            # Source: first part before first | of ref_id
            source_assembly = ref_id.split("|")[0]

            # Find genomes
            source_fna = find_genome_fna(args.genome_dir, source_assembly)
            target_fna = find_genome_fna(args.genome_dir, target_assembly)

            if not source_fna or not target_fna:
                out.write(f"{category}\t{ref_id}\t{target_assembly}\tMISSING_GENOME\t-\tno_validation\t-\t-\n")
                continue

            # Parse anchor positions
            details = parse_detail(eg["details"])
            up_hits = details["up"]
            down_hits = details["down"]
            if not up_hits or not down_hits:
                continue

            # Pick the first up and down hit for region extraction
            up_hit = up_hits[0]
            down_hit = down_hits[0]

            # Extract surrounding region from target genome
            tgt_region_dir = os.path.join(args.out_dir, f"{category}_{i}")
            os.makedirs(tgt_region_dir, exist_ok=True)

            # For inversion/duplication: same contig, extract whole region between anchors with flank
            # For translocation: different contigs, extract from each contig

            up_loc = f"{up_hit['contig']}:{up_hit['start']}-{up_hit['end']}({up_hit['strand']})"
            dn_loc = f"{down_hit['contig']}:{down_hit['start']}-{down_hit['end']}({down_hit['strand']})"

            target_region_fa = os.path.join(tgt_region_dir, "target_region.fa")

            if up_hit["contig"] == down_hit["contig"]:
                # Same contig: extract spanning region with flanks
                positions = sorted([up_hit["start"], up_hit["end"], down_hit["start"], down_hit["end"]])
                rstart = max(1, positions[0] - args.region_flank)
                rend = positions[-1] + args.region_flank
                extract_region(target_fna, up_hit["contig"], rstart, rend, target_region_fa)
            else:
                # Different contigs: concatenate two regions
                # Just use up contig for simplicity
                rstart = max(1, up_hit["start"] - args.region_flank)
                rend = up_hit["end"] + args.region_flank
                extract_region(target_fna, up_hit["contig"], rstart, rend, target_region_fa)

            # Extract source region (around the IS110 transposase position)
            # For source we just use the whole genome - small enough

            # Run nucmer: source as ref, target_region as query
            out_prefix = os.path.join(tgt_region_dir, "alignment")
            blocks = run_nucmer(source_fna, target_region_fa, out_prefix)

            # Classify
            verdict, _ = classify_alignment(blocks, category)

            # Summary
            n_blocks = len(blocks)
            n_inv = sum(1 for b in blocks if b.get("qs", 0) > b.get("qe", 0))
            n_fwd = n_blocks - n_inv
            summary = f"n_blocks={n_blocks} fwd={n_fwd} inv={n_inv}"

            out.write(f"{category}\t{ref_id}\t{target_assembly}\t{up_loc}\t{dn_loc}\t{verdict}\t{summary}\t{tgt_region_dir}\n")
            print(f"  [{i+1}] {ref_id[:30]}... → {target_assembly}  →  {verdict}  ({summary})",
                  file=sys.stderr)

    out.close()
    print(f"\nReport: {report_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
