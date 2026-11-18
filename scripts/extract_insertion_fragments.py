#!/usr/bin/env python3
"""Extract inserted DNA fragments from split IS110 cases.

Reads the insertions TSV produced by detect_is110_variants.py, locates the
corresponding target contig in the Logan dataset, and writes each inserted
fragment (the target gap region between IS110 part1 and part2) as FASTA.
"""

import argparse
import csv
import os
import subprocess
import sys
from collections import defaultdict


LOCUS_TOL = 20   # as detect_is110_variants.py --locus-tol


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("insertions_tsv",
                   help="Path to *_insertions.tsv from detect_is110_variants.py")
    p.add_argument("logan_dir", help="Directory with Logan *.contigs.fa.zst files")
    p.add_argument("output_fasta", help="Output FASTA file for inserted fragments")
    p.add_argument("--pad", type=int, default=0,
                   help="Extra flanking bases to include (default: 0)")
    p.add_argument("--min-len", type=int, default=300,
                   help="Minimum fragment length (default: 300)")
    p.add_argument("--max-len", type=int, default=50000,
                   help="Maximum fragment length (default: 50000)")
    return p.parse_args()


def accession_from_contig(contig_name):
    """SRR15695699_126 -> SRR15695699"""
    return contig_name.rsplit("_", 1)[0]


def load_contig_sequence(zst_path, contig_name):
    """Decompress the .zst file and return the specific contig sequence."""
    try:
        result = subprocess.run(
            ["zstd", "-dcq", zst_path],
            capture_output=True, text=True, check=True
        )
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None

    current_name = None
    current_seq = []
    for line in result.stdout.split("\n"):
        if line.startswith(">"):
            if current_name == contig_name and current_seq:
                return "".join(current_seq)
            current_name = line[1:].split()[0]
            current_seq = []
        else:
            current_seq.append(line.strip())
    if current_name == contig_name and current_seq:
        return "".join(current_seq)
    return None


def main():
    args = parse_args()

    # Group insertions by accession to minimize file reads
    by_accession = defaultdict(list)
    with open(args.insertions_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            target = r["target_contig"]
            acc = accession_from_contig(target)

            # Figure out the inserted region bounds (target gap)
            p1_ts, p1_te = int(r["part1_ts"]), int(r["part1_te"])
            p2_ts, p2_te = int(r["part2_ts"]), int(r["part2_te"])
            if p1_ts < p2_ts:
                ins_start, ins_end = p1_te, p2_ts
            else:
                ins_start, ins_end = p2_te, p1_ts

            ins_len = ins_end - ins_start
            if ins_len < args.min_len or ins_len > args.max_len:
                continue

            # Apply padding
            ins_start = max(0, ins_start - args.pad)
            ins_end = ins_end + args.pad

            by_accession[acc].append({
                "target": target,
                "ins_start": ins_start,
                "ins_end": ins_end,
                "query": r["query"],
                "ident": r.get("part1_ident", "?"),
            })

    print(f"Accessions to process: {len(by_accession):,}")

    written = 0
    seen_loci = {}
    n_dup = 0
    failed_files = 0
    failed_contigs = 0

    with open(args.output_fasta, "w") as out:
        for i, (acc, frags) in enumerate(by_accession.items()):
            if i % 500 == 0:
                print(f"  Progress: {i}/{len(by_accession)} accessions, {written} fragments written",
                      flush=True)

            zst_file = os.path.join(args.logan_dir, f"{acc}.contigs.fa.zst")
            if not os.path.exists(zst_file):
                failed_files += len(frags)
                continue

            # Decompress once, extract all needed contigs
            try:
                result = subprocess.run(
                    ["zstd", "-dcq", zst_file],
                    capture_output=True, text=True, check=True
                )
            except subprocess.CalledProcessError:
                failed_files += len(frags)
                continue

            # Build dict of contigs from this file
            wanted = set(f["target"] for f in frags)
            contigs = {}
            current_name = None
            current_seq = []
            for line in result.stdout.split("\n"):
                if line.startswith(">"):
                    if current_name in wanted and current_seq:
                        contigs[current_name] = "".join(current_seq)
                    current_name = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line.strip())
            if current_name in wanted and current_seq:
                contigs[current_name] = "".join(current_seq)

            for f in frags:
                seq = contigs.get(f["target"])
                if seq is None:
                    failed_contigs += 1
                    continue
                # clamp (a padded interval can run past the contig end) and
                # write each locus ONCE: the same insert was written once per
                # reference IS that matched it, which inflated cluster sizes
                # the stored coordinates are the clamped ones (the header said
                # __1000_3600 for 2000 bases of a 3000-bp contig)
                f["ins_start"] = max(0, f["ins_start"])
                f["ins_end"] = min(len(seq), f["ins_end"])
                frag = seq[f["ins_start"]:f["ins_end"]]
                if len(frag) < args.min_len:
                    continue
                # same locus = both ends within LOCUS_TOL bp of a locus already
                # written (floor bins split 1019 / 1020 into two "loci")
                near = seen_loci.setdefault(f["target"], [])
                if any(abs(a - f["ins_start"]) <= LOCUS_TOL and abs(b - f["ins_end"]) <= LOCUS_TOL
                       for a, b in near):
                    n_dup += 1
                    continue
                near.append((f["ins_start"], f["ins_end"]))

                hdr = (f">{f['target']}__{f['ins_start']}_{f['ins_end']}"
                       f" query={f['query']}"
                       f" len={len(frag)}")
                out.write(hdr + "\n")
                for j in range(0, len(frag), 80):
                    out.write(frag[j:j+80] + "\n")
                written += 1

    print(f"\nTotal fragments written: {written:,}  (duplicate loci skipped: {n_dup:,})")
    print(f"Failed file lookups: {failed_files:,}")
    print(f"Failed contig lookups: {failed_contigs:,}")


if __name__ == "__main__":
    main()
