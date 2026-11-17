#!/usr/bin/env python3
"""Find IS110 element boundaries using the transposase protein/CDS as anchor.

Approach:
  1. For each IS110 ref, take the transposase CDS as the anchor.
  2. minimap2 the transposase CDS against the genome database to find all hits.
  3. For each hit, extract surrounding context from the target genome.
  4. Align all these contexts against the original full reference (which has
     30 kb flanks).
  5. Compute coverage depth at each position of the reference.
  6. The IS110 boundary = positions where coverage drops sharply on both sides
     of the transposase.

Conservation = many genomes share this position (IS110 element + maybe some
flanking that's conserved).
Variable = few genomes share this position (genomic flank, host-specific).

Usage:
    protein_anchor_boundary.py --large-refs <large_refs.fa>
                               --ref-table <ref_table.tsv>
                               --db <ncbi_db.fa>
                               --out-dir <out_dir>
"""

import argparse
import os
import subprocess
import sys
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--large-refs", required=True,
                   help="FASTA of large-flank references (from extract_multi_anchors.py)")
    p.add_argument("--ref-table", required=True)
    p.add_argument("--db", required=True,
                   help="NCBI E. coli concat FASTA")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--context-flank", type=int, default=15000,
                   help="Length of context to extract around each transposase hit")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--min-tnp-identity", type=float, default=90,
                   help="Min %% identity for transposase hits in DB")
    p.add_argument("--min-tnp-cov", type=float, default=80,
                   help="Min %% query coverage for transposase hits")
    return p.parse_args()


def parse_fasta_index(path):
    """Yield (name, seq) tuples."""
    name = None
    parts = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    yield name, "".join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.strip())
        if name:
            yield name, "".join(parts)


def revcomp(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # Load reference table
    ref_info = {}
    with open(args.ref_table) as f:
        next(f)
        for line in f:
            parts = line.rstrip().split("\t")
            ref_id, ref_len, anchor_len, tnp_len, flank_size, _ = parts
            ref_info[ref_id] = {
                "ref_len": int(ref_len), "anchor_len": int(anchor_len),
                "tnp_len": int(tnp_len), "flank_size": int(flank_size),
            }

    # Load reference sequences
    ref_seqs = {}
    for name, seq in parse_fasta_index(args.large_refs):
        if name in ref_info:
            ref_seqs[name] = seq

    print(f"Loaded {len(ref_seqs)} references", file=sys.stderr)

    # Extract transposase CDS from each ref (it sits at flank_size : flank_size + tnp_len)
    tnp_fa = os.path.join(args.out_dir, "transposases.fa")
    with open(tnp_fa, "w") as f:
        for ref_id, seq in ref_seqs.items():
            info = ref_info[ref_id]
            tnp_seq = seq[info["flank_size"] : info["flank_size"] + info["tnp_len"]]
            f.write(f">{ref_id}\n{tnp_seq}\n")

    n_tnp = sum(1 for _ in parse_fasta_index(tnp_fa))
    print(f"Wrote {n_tnp} transposase CDS sequences to {tnp_fa}", file=sys.stderr)

    # Step 1: search transposases against NCBI db
    tnp_paf = os.path.join(args.out_dir, "tnp_vs_db.paf")
    if os.path.exists(tnp_paf) and os.path.getsize(tnp_paf) > 0:
        print(f"Reusing existing PAF: {tnp_paf}", file=sys.stderr)
    else:
        print(f"minimap2 transposase CDS vs DB ({args.threads} threads)...", file=sys.stderr)
        subprocess.run([
            "minimap2", "-x", "asm5", "-c", "--eqx",
            "--secondary=yes", "-N", "200", "-p", "0.5",
            "-t", str(args.threads),
            args.db, tnp_fa,
            "-o", tnp_paf,
        ], check=True)

    # Step 2: parse hits, organize by ref_id
    hits_by_ref = defaultdict(list)
    n_total = 0
    n_kept = 0
    with open(tnp_paf) as f:
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
            if ident < args.min_tnp_identity or cov < args.min_tnp_cov:
                continue
            n_kept += 1
            hits_by_ref[qname].append({
                "tname": tname, "tlen": tlen, "ts": ts, "te": te,
                "strand": strand, "ident": ident, "cov": cov,
            })

    print(f"Transposase hits: {n_total:,} total, {n_kept:,} pass filter", file=sys.stderr)
    print(f"Refs with hits: {len(hits_by_ref)}", file=sys.stderr)

    # Step 3: for each ref, extract context from each hit and align back to ref
    # Build a contexts FASTA per ref, then align with minimap2
    summary_path = os.path.join(args.out_dir, "boundary_summary.tsv")
    summary_out = open(summary_path, "w")
    summary_out.write("ref_id\tn_hits\ttnp_len\tflank_size\tref_len\t"
                      "left_boundary\tright_boundary\tis_element_length\t"
                      "left_drop_off\tright_drop_off\n")

    # Cache target genome sequences
    db_index = None  # lazy load on first hit

    def load_db_seq(tname, ts, te):
        """Load specific region from DB. We need an indexed DB - use samtools faidx on demand."""
        # Use samtools faidx (assumes db.fa.fai exists or is created)
        cmd = ["samtools", "faidx", args.db, f"{tname}:{max(1, ts+1)}-{te}"]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            return None
        seq_parts = []
        for line in result.stdout.split("\n"):
            if not line.startswith(">") and line:
                seq_parts.append(line.strip())
        return "".join(seq_parts)

    # Make sure DB is indexed
    if not os.path.exists(args.db + ".fai"):
        print(f"Indexing DB with samtools faidx...", file=sys.stderr)
        subprocess.run(["samtools", "faidx", args.db], check=True)

    for ref_id, hits in hits_by_ref.items():
        info = ref_info[ref_id]
        ref_len = info["ref_len"]
        ref_seq = ref_seqs[ref_id]

        # Extract contexts: per hit, get +/- context_flank around the transposase region
        ctx_fa = os.path.join(args.out_dir, f"context_{ref_id.replace('|', '_')}.fa")
        with open(ctx_fa, "w") as f:
            for i, h in enumerate(hits):
                ctx_start = max(0, h["ts"] - args.context_flank)
                ctx_end = min(h["tlen"], h["te"] + args.context_flank)
                ctx_seq = load_db_seq(h["tname"], ctx_start, ctx_end)
                if not ctx_seq:
                    continue
                if h["strand"] == "-":
                    ctx_seq = revcomp(ctx_seq)
                f.write(f">hit_{i}\n{ctx_seq}\n")

        # Align contexts back to the original reference
        ctx_paf = os.path.join(args.out_dir, f"context_{ref_id.replace('|', '_')}_vs_ref.paf")
        # Make a per-ref FASTA target
        ref_fa = os.path.join(args.out_dir, f"ref_{ref_id.replace('|', '_')}.fa")
        with open(ref_fa, "w") as f:
            f.write(f">{ref_id}\n{ref_seq}\n")

        subprocess.run([
            "minimap2", "-x", "asm5", "-c", "--eqx",
            "--secondary=yes", "-N", "5", "-p", "0.5",
            "-t", str(args.threads),
            ref_fa, ctx_fa,
            "-o", ctx_paf,
        ], check=False)

        # Compute coverage depth at each position of the reference
        depth = [0] * ref_len
        n_aligned = 0
        with open(ctx_paf) as f:
            for line in f:
                c = line.rstrip().split("\t")
                if len(c) < 12:
                    continue
                ts, te = int(c[7]), int(c[8])
                matches, block = int(c[9]), int(c[10])
                ident = matches / block * 100 if block > 0 else 0
                if ident < 90:
                    continue
                n_aligned += 1
                for p in range(ts, min(te, ref_len)):
                    depth[p] += 1

        # Find boundaries: scan outward from transposase, find where depth drops
        # to less than half its peak in the transposase region
        tnp_start = info["flank_size"]
        tnp_end = info["flank_size"] + info["tnp_len"]
        # Peak = max depth around transposase
        peak = max(depth[tnp_start:tnp_end]) if tnp_start < tnp_end else 0
        threshold = max(2, peak // 2)

        # Walk left from tnp_start
        left_boundary = tnp_start
        for p in range(tnp_start - 1, -1, -1):
            if depth[p] < threshold:
                left_boundary = p + 1
                break
        else:
            left_boundary = 0

        # Walk right from tnp_end
        right_boundary = tnp_end
        for p in range(tnp_end, ref_len):
            if depth[p] < threshold:
                right_boundary = p
                break
        else:
            right_boundary = ref_len

        is_element_length = right_boundary - left_boundary
        left_drop_off = left_boundary - tnp_start
        right_drop_off = right_boundary - tnp_end

        summary_out.write(f"{ref_id}\t{n_aligned}\t{info['tnp_len']}\t{info['flank_size']}\t"
                          f"{ref_len}\t{left_boundary}\t{right_boundary}\t{is_element_length}\t"
                          f"{left_drop_off}\t{right_drop_off}\n")

        # Save coverage profile (downsampled for visualization)
        prof_path = os.path.join(args.out_dir, f"depth_{ref_id.replace('|', '_')}.tsv")
        bin_size = 100
        with open(prof_path, "w") as pf:
            pf.write("position\tdepth\n")
            for p in range(0, ref_len, bin_size):
                bin_end = min(p + bin_size, ref_len)
                avg = sum(depth[p:bin_end]) / (bin_end - p)
                pf.write(f"{p}\t{avg:.1f}\n")

        print(f"  {ref_id}: {n_aligned} contexts aligned, peak depth={peak}, "
              f"boundaries=[{left_boundary},{right_boundary}], "
              f"IS110 length = {is_element_length} bp "
              f"(extends {-left_drop_off} bp left, {right_drop_off} bp right of transposase)",
              file=sys.stderr)

    summary_out.close()
    print(f"Summary: {summary_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
