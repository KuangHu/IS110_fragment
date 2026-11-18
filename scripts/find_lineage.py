#!/usr/bin/env python3
"""Find IS110 size/fragmentation lineages across genomes.

For each IS110 query that has multi-split observations, check whether we see
a progressive series V1 -> V2 -> V3 -> V4 where each version adds an insertion
to the previous one (nested pattern).

Steps:
  1. Parse the multisplit chain TSV (from detect_multisplit.py)
  2. Group observations per IS110 query
  3. Within each query, bucket target contigs by split count (0, 1, 2, ...)
  4. For each non-zero split group, extract inserted sequences from Logan contigs
  5. Check nesting: do N-split insertions include (N-1)-split insertions?

Usage:
    python find_lineage.py <chains_tsv> <logan_dir> <output_prefix>
"""

import argparse
import csv
import os
import sys
import re
import subprocess
import tempfile
from collections import defaultdict, Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paf_utils import aligned_identity  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("chains_tsv", help="TSV of multi-split chains (from detect_multisplit.py)")
    p.add_argument("logan_dir", help="Directory of Logan .contigs.fa.zst files")
    p.add_argument("output_prefix", help="Output prefix for result files")
    p.add_argument("--min-ladder-depth", type=int, default=3,
                   help="Min number of distinct split counts (e.g. 3 means V1/V2/V3) (default: 3)")
    p.add_argument("--cluster-identity", type=float, default=95,
                   help="Min identity %% to consider two inserted seqs the same (default: 95)")
    p.add_argument("--max-ins-examples", type=int, default=5,
                   help="Max example target contigs per split-count group (default: 5)")
    p.add_argument("--threads", type=int, default=4)
    return p.parse_args()


SEG_RE = re.compile(r"Q\[(\d+)-(\d+)\]->T\[(\d+)-(\d+)\]\((\d+)%\)")


def parse_segments(segments_str):
    """Parse the 'segments' column: 'Q[a-b]->T[c-d](99%);Q[...]->T[...](99%)'"""
    segs = []
    for m in SEG_RE.finditer(segments_str):
        segs.append({
            "qs": int(m.group(1)), "qe": int(m.group(2)),
            "ts": int(m.group(3)), "te": int(m.group(4)),
            "ident": int(m.group(5)),
        })
    return segs


def extract_insertions(segments, strand):
    """Given a sorted chain of segments, return list of insertion spans on target.

    For a co-linear chain on a single strand, insertion at position i is:
      - target region between segment[i].te and segment[i+1].ts (forward)
      - target region between segment[i+1].te and segment[i].ts (reverse)
    """
    insertions = []
    # Sort by query position (should already be sorted)
    segs = sorted(segments, key=lambda x: x["qs"])
    for i in range(len(segs) - 1):
        a, b = segs[i], segs[i + 1]
        if strand == "+":
            ts, te = a["te"], b["ts"]
        else:
            ts, te = b["te"], a["ts"]
        if te > ts:
            insertions.append({"ts": ts, "te": te, "len": te - ts,
                               "qs_before": a["qe"], "qs_after": b["qs"]})
    return insertions


def accession_from_contig(contig_name):
    return contig_name.rsplit("_", 1)[0]


def load_contig_sequence(zst_path, contig_name):
    """Extract a single contig sequence from a zst fasta."""
    result = subprocess.run(
        ["zstd", "-dcq", zst_path],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        return None
    current = None
    found_parts = []
    for line in result.stdout.split("\n"):
        if line.startswith(">"):
            if current == contig_name and found_parts:
                return "".join(found_parts)
            current = line[1:].split()[0]
            found_parts = []
        elif line and current == contig_name:
            found_parts.append(line.strip())
    if current == contig_name and found_parts:
        return "".join(found_parts)
    return None


def cluster_sequences(seqs, identity_thresh, threads):
    """Cluster sequences using minimap2 pairwise, return cluster assignments.

    Returns: list of cluster IDs (same length as seqs), and dict of representatives.
    """
    if not seqs:
        return [], {}
    if len(seqs) == 1:
        return [0], {0: seqs[0]}

    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as fa:
        for i, s in enumerate(seqs):
            fa.write(f">s{i}\n{s}\n")
        fa_path = fa.name

    try:
        # All-vs-all minimap2
        result = subprocess.run(
            ["minimap2", "-x", "ava-ont", "-c", "--eqx",
             "-t", str(threads), fa_path, fa_path],
            capture_output=True, text=True
        )
        # Union-find clustering
        parent = list(range(len(seqs)))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(a, b):
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb

        for line in result.stdout.split("\n"):
            c = line.strip().split("\t")
            if len(c) < 12:
                continue
            qname, tname = c[0], c[5]
            if qname == tname:
                continue
            ident = aligned_identity(c)   # gap columns excluded (paf_utils)
            if ident < identity_thresh:
                continue
            qi = int(qname[1:])
            ti = int(tname[1:])
            # require substantial overlap of the smaller sequence
            qs, qe = int(c[2]), int(c[3])
            qlen = int(c[1])
            if (qe - qs) / qlen < 0.7:
                continue
            union(qi, ti)

        # Assign cluster IDs
        roots = {}
        cluster_ids = []
        for i in range(len(seqs)):
            r = find(i)
            if r not in roots:
                roots[r] = len(roots)
            cluster_ids.append(roots[r])

        # Get representatives (longest sequence per cluster)
        reps = {}
        for i, cid in enumerate(cluster_ids):
            if cid not in reps or len(seqs[i]) > len(reps[cid]):
                reps[cid] = seqs[i]

        return cluster_ids, reps
    finally:
        os.remove(fa_path)


def main():
    args = parse_args()

    print(f"Loading chains from {args.chains_tsv}")
    # Parse chains TSV
    rows_by_query = defaultdict(list)
    total_rows = 0
    with open(args.chains_tsv) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            total_rows += 1
            n_splits = int(r["n_splits"])
            # Skip zero-split for now (load them separately if needed)
            rows_by_query[r["query"]].append({
                "tname": r["target"],
                "tlen": int(r["tlen"]),
                "strand": r["strand"],
                "n_splits": n_splits,
                "segments": parse_segments(r["segments"]),
            })
    print(f"  Loaded {total_rows:,} rows, {len(rows_by_query):,} unique queries")

    # Filter to queries with a "ladder": observations at >= min_ladder_depth different split counts
    ladders = {}
    for qname, obs in rows_by_query.items():
        counts = Counter(o["n_splits"] for o in obs)
        if len(counts) >= args.min_ladder_depth:
            ladders[qname] = {k: [o for o in obs if o["n_splits"] == k]
                              for k in sorted(counts.keys())}

    print(f"  Queries with ladder depth >= {args.min_ladder_depth}: {len(ladders)}")
    print()

    # Write a summary of ladders
    ladder_path = f"{args.output_prefix}_ladders_summary.tsv"
    with open(ladder_path, "w") as f:
        f.write("query\tmax_splits\tn_split_levels\tobs_per_level\ttotal_observations\n")
        for qname, levels in sorted(ladders.items(), key=lambda x: -max(x[1].keys())):
            obs_counts = ",".join(f"V{k}:{len(v)}" for k, v in levels.items())
            total = sum(len(v) for v in levels.values())
            f.write(f"{qname}\t{max(levels.keys())}\t{len(levels)}\t{obs_counts}\t{total}\n")
    print(f"Ladder summary: {ladder_path}")

    # For each ladder, extract insertions and check for nesting
    # Process the deepest ladders first
    nesting_path = f"{args.output_prefix}_nesting.tsv"
    with open(nesting_path, "w") as out:
        out.write("query\tmax_splits\tlevel\tn_example\texample_target\tstrand\t"
                  "insertion_lengths\tinsertion_clusters\tinsertion_sequences\n")

        for qname, levels in sorted(ladders.items(),
                                    key=lambda x: (-max(x[1].keys()), -sum(len(v) for v in x[1].values())))[:30]:
            print(f"\n=== {qname[:20]}... (max splits: {max(levels.keys())}, levels: {len(levels)}) ===")
            for n_split in sorted(levels.keys()):
                if n_split == 0:
                    # Intact: no insertions; skip but record
                    obs = levels[n_split][:args.max_ins_examples]
                    for i, o in enumerate(obs):
                        out.write(f"{qname}\t{max(levels.keys())}\t{n_split}\t{i}\t{o['tname']}\t"
                                  f"{o['strand']}\tintact\tintact\t\n")
                    continue

                obs = levels[n_split][:args.max_ins_examples]
                # Extract inserted sequences from each observation
                all_ins_seqs = []
                all_ins_meta = []
                for i, o in enumerate(obs):
                    sample = accession_from_contig(o["tname"])
                    zst = os.path.join(args.logan_dir, f"{sample}.contigs.fa.zst")
                    if not os.path.exists(zst):
                        continue
                    contig_seq = load_contig_sequence(zst, o["tname"])
                    if not contig_seq:
                        continue
                    inss = extract_insertions(o["segments"], o["strand"])
                    ins_seqs = []
                    ins_lens = []
                    for ins in inss:
                        ins_lens.append(ins["len"])
                        s = contig_seq[ins["ts"]:ins["te"]]
                        if o["strand"] == "-":
                            comp = str.maketrans("ACGTacgt", "TGCAtgca")
                            s = s.translate(comp)[::-1]
                        ins_seqs.append(s)
                    all_ins_seqs.append((i, o, ins_seqs, ins_lens))

                if not all_ins_seqs:
                    continue

                # Flatten all insertion sequences across this level
                flat_seqs = []
                index_map = []  # (obs_idx, within_obs_idx)
                for obs_idx, o, seqs, _ in all_ins_seqs:
                    for j, s in enumerate(seqs):
                        if len(s) >= 200:
                            flat_seqs.append(s)
                            index_map.append((obs_idx, j))

                # Cluster them
                cluster_ids, reps = cluster_sequences(flat_seqs, args.cluster_identity, args.threads)

                # Report
                print(f"\n  V{n_split}: {len(levels[n_split])} observations (examples: {len(obs)})")
                cluster_summary = Counter(cluster_ids) if cluster_ids else Counter()
                print(f"    Inserted sequence clusters: {len(cluster_summary)}")
                for cid, cnt in cluster_summary.most_common():
                    print(f"      Cluster {cid}: {cnt} occurrences (rep len={len(reps[cid])})")

                # Write per-observation row
                # Build cluster membership per observation
                per_obs_clusters = defaultdict(list)
                for idx, (obs_idx, j) in enumerate(index_map):
                    per_obs_clusters[obs_idx].append(cluster_ids[idx])

                for obs_idx, o, seqs, ins_lens in all_ins_seqs:
                    clusters_str = ",".join(str(c) for c in per_obs_clusters.get(obs_idx, []))
                    lens_str = ",".join(str(l) for l in ins_lens)
                    seqs_preview = ";".join(s[:40] + "..." if len(s) > 40 else s for s in seqs)
                    out.write(f"{qname}\t{max(levels.keys())}\t{n_split}\t{obs_idx}\t"
                              f"{o['tname']}\t{o['strand']}\t{lens_str}\t{clusters_str}\t{seqs_preview}\n")

    print(f"\nNesting details: {nesting_path}")


if __name__ == "__main__":
    main()
