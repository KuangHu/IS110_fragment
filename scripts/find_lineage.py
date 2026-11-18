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
  5. Check nesting: an observation with N insertions is "supported" by one
     with N-1 when it carries all of that one's insertions (same sequence
     cluster, same place in the IS). Insertions of all levels are clustered
     together. A different number of insertions alone is NOT a lineage.

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
from paf_utils import aligned_identity, aligned_query_bases  # noqa: E402


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


MIN_BOTH_COV = 0.9     # of the LONGER of two inserts must be aligned


def matching_pairs(seqs, identity_thresh, threads, min_both_cov=MIN_BOTH_COV):
    """Pairs (i, j), i < j, of sequences that are the same insert: one
    alignment row with identity >= threshold whose aligned bases are at least
    `min_both_cov` of the LONGER sequence (so of both).

    The test used to be 70% of the SHORTER sequence: two 1500-bp inserts that
    share an 1100-bp segment and nothing else were "the same insert", and so
    was any insert with a fragment of it. The measure is symmetric, so the
    result does not depend on the order of the input (ava-ont prints each
    pair in one direction only)."""
    if len(seqs) < 2:
        return set()
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as fa:
        for i, s in enumerate(seqs):
            fa.write(f">s{i}\n{s}\n")
        fa_path = fa.name
    try:
        result = subprocess.run(
            ["minimap2", "-x", "ava-ont", "-c", "--eqx",
             "-t", str(threads), fa_path, fa_path],
            capture_output=True, text=True
        )
        if result.returncode != 0:
            # a failed run gave "every sequence its own cluster"
            raise RuntimeError(f"minimap2 failed ({result.returncode}): {result.stderr[-500:]}")
    finally:
        os.remove(fa_path)
    pairs = set()
    for line in result.stdout.split("\n"):
        c = line.strip().split("\t")
        if len(c) < 12 or c[0] == c[5]:
            continue
        if aligned_identity(c) < identity_thresh:      # gap columns excluded (paf_utils)
            continue
        if aligned_query_bases(c) / max(int(c[1]), int(c[6])) < min_both_cov:
            continue
        i, j = sorted((int(c[0][1:]), int(c[5][1:])))
        pairs.add((i, j))
    return pairs


def cluster_sequences(seqs, identity_thresh, threads, pairs=None):
    """Single-linkage clusters of `matching_pairs`, for REPORTING.

    Returns: list of cluster IDs (same length as seqs), and dict of
    representatives. Membership of one cluster is transitive (A~B and B~C put
    A and C together); the nesting test therefore uses the pairs themselves,
    not the cluster ids."""
    if not seqs:
        return [], {}
    if pairs is None:
        pairs = matching_pairs(seqs, identity_thresh, threads)
    parent = list(range(len(seqs)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, j in sorted(pairs):
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[max(ri, rj)] = min(ri, rj)
    roots = {}
    cluster_ids = []
    for i in range(len(seqs)):
        r = find(i)
        if r not in roots:
            roots[r] = len(roots)
        cluster_ids.append(roots[r])
    reps = {}
    for i, cid in enumerate(cluster_ids):
        if cid not in reps or len(seqs[i]) > len(reps[cid]):
            reps[cid] = seqs[i]
    return cluster_ids, reps


JUNCTION_TOL = 30      # bp on the IS: same insertion point


def nesting_parents(child, parents, same=None):
    """Indices of ALL observations in `parents` whose insertions are all
    present in `child` (the same insert at the same place in the IS), child
    having exactly one more. An insertion is (key, position on the IS);
    `same(key_a, key_b)` says whether two inserts are the same sequence
    (default: equal keys)."""
    if same is None:
        same = lambda x, y: x == y
    found = []
    for k, par in enumerate(parents):
        if len(child) != len(par) + 1:
            continue
        left = list(child)
        ok = True
        for cid, pos in par:
            hit = next((x for x in left if same(x[0], cid) and abs(x[1] - pos) <= JUNCTION_TOL), None)
            if hit is None:
                ok = False
                break
            left.remove(hit)
        if ok:
            found.append(k)
    return found


def nesting_parent(child, parents, same=None):
    """The first of `nesting_parents`, or None."""
    found = nesting_parents(child, parents, same)
    return found[0] if found else None


def nesting_graph(levels, same=None):
    """levels: {number of insertions: [observation, ...]}, an observation
    being its list of (key, position) or None when it could not be examined.

    Returns (parents, reach): parents[(level, i)] = indices of ALL supporting
    observations of level - 1, the one with the longest supported chain
    first; reach[(level, i)] = length of the longest supported chain ending
    in that observation. Taking the first compatible parent made the depth
    depend on the order of the observations: ABC was hung under BC (no
    parent, depth 2) although A -> AB -> ABC was there (depth 3)."""
    parents, reach = {}, {}
    for n in sorted(levels):
        prev = levels.get(n - 1) or []
        cand = [(k, p) for k, p in enumerate(prev) if p is not None]
        for i, obs in enumerate(levels[n]):
            parents[(n, i)] = []
            reach[(n, i)] = 1
            if obs is None or not cand:
                continue
            ks = [cand[k][0] for k in nesting_parents(obs, [p for _, p in cand], same)]
            ks.sort(key=lambda k: (-reach[(n - 1, k)], k))
            parents[(n, i)] = ks
            if ks:
                reach[(n, i)] = reach[(n - 1, ks[0])] + 1
    return parents, reach


def main():
    args = parse_args()

    print(f"Loading chains from {args.chains_tsv}")
    rows_by_query = defaultdict(list)
    total_rows = skipped = 0
    with open(args.chains_tsv) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            total_rows += 1
            status = r.get("status") or ("intact" if int(r["n_splits"]) == 0 else "variant")
            types = [t for t in r["event_types"].split(",") if t]
            # A ladder is a series of INSERTIONS. Rows with a deletion or an
            # unresolved junction, and partial copies, are not rungs of it
            # (their split count used to be put in the same buckets).
            if status == "partial" or any(t not in ("insertion", "intact") for t in types):
                skipped += 1
                continue
            rows_by_query[r["query"]].append({
                "tname": r["target"],
                "tlen": int(r["tlen"]),
                "strand": r["strand"],
                "n_splits": sum(t == "insertion" for t in types),
                "segments": parse_segments(r["segments"]),
            })
    print(f"  Loaded {total_rows:,} rows ({skipped:,} partial / deletion rows not used), "
          f"{len(rows_by_query):,} unique queries")

    # Candidate ladders: observations at >= min_ladder_depth different insertion counts.
    # This is only a candidate: different counts say nothing about inheritance.
    ladders = {}
    for qname, obs in rows_by_query.items():
        counts = Counter(o["n_splits"] for o in obs)
        if len(counts) >= args.min_ladder_depth:
            ladders[qname] = {k: [o for o in obs if o["n_splits"] == k]
                              for k in sorted(counts.keys())}

    print(f"  Queries with >= {args.min_ladder_depth} different insertion counts "
          f"(candidate ladders): {len(ladders)}")
    print()

    order = sorted(ladders.items(),
                   key=lambda x: (-max(x[1].keys()), -sum(len(v) for v in x[1].values()), x[0]))
    support = {}        # qname -> (supported depth, n supported edges) for examined ladders

    nesting_path = f"{args.output_prefix}_nesting.tsv"
    with open(nesting_path, "w") as out:
        out.write("query\tmax_splits\tlevel\tn_example\texample_target\tstrand\t"
                  "insertion_lengths\tinsertion_clusters\tinsertion_positions\t"
                  "parent_example\tinsertion_sequences\n")

        for qname, levels in order[:30]:
            print(f"\n=== {qname[:20]}... (max insertions: {max(levels.keys())}, levels: {len(levels)}) ===")
            # 1. collect the insertions of the examined observations of ALL levels
            examined = {}                      # level -> [obs dict with 'ins']
            flat_seqs, flat_ref = [], []       # sequences to cluster, (level, obs idx, ins idx)
            for n_split in sorted(levels):
                examined[n_split] = []
                for o in levels[n_split][:args.max_ins_examples]:
                    rec = {"o": o, "ins": [], "ok": n_split == 0}
                    if n_split > 0:
                        sample = accession_from_contig(o["tname"])
                        zst = os.path.join(args.logan_dir, f"{sample}.contigs.fa.zst")
                        contig_seq = load_contig_sequence(zst, o["tname"]) if os.path.exists(zst) else None
                        if contig_seq:
                            rec["ok"] = True
                            for ins in extract_insertions(o["segments"], o["strand"]):
                                seq = contig_seq[ins["ts"]:ins["te"]]
                                if o["strand"] == "-":
                                    seq = seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]
                                rec["ins"].append({"len": ins["len"], "pos": ins["qs_before"],
                                                   "seq": seq, "cluster": None})
                    examined[n_split].append(rec)
                    for j, ins in enumerate(rec["ins"]):
                        if len(ins["seq"]) >= 200:
                            flat_seqs.append(ins["seq"])
                            flat_ref.append(ins)
            # 2. ONE comparison over all levels, so that inserts can be
            #    compared between levels (each level used to be clustered on
            #    its own: "cluster 0" of V1 and of V2 were unrelated)
            pairs = matching_pairs(flat_seqs, args.cluster_identity, args.threads)
            cluster_ids, _ = cluster_sequences(flat_seqs, args.cluster_identity, args.threads, pairs)
            for k, (ins, cid) in enumerate(zip(flat_ref, cluster_ids)):
                ins["cluster"] = cid
                ins["key"] = k
            n_extra = 0
            for recs in examined.values():
                for rec in recs:
                    for ins in rec["ins"]:
                        if ins["cluster"] is None:        # too short to compare: its own group
                            n_extra += 1
                            ins["cluster"] = ins["key"] = -n_extra
            # 3. nesting: an observation with N insertions is supported by one
            #    with N-1 if it contains all of that one's insertions. "The
            #    same insert" is a DIRECT match of the two sequences (not
            #    membership of one single-linkage cluster).
            same = lambda x, y: x == y or (min(x, y), max(x, y)) in pairs
            parents, reach = nesting_graph(
                {n: [[(x["key"], x["pos"]) for x in rec["ins"]] if rec["ok"] else None
                     for rec in recs] for n, recs in examined.items()}, same)
            edges = sum(len(v) for v in parents.values())
            for (n_split, i), ks in parents.items():
                examined[n_split][i]["parent"] = ks
            depth = max(reach.values()) if reach else 0
            support[qname] = (depth, edges)
            print(f"  supported nesting edges: {edges}; longest supported chain: {depth} levels")
            for n_split in sorted(examined):
                for i, rec in enumerate(examined[n_split]):
                    o = rec["o"]
                    if n_split == 0:
                        out.write(f"{qname}\t{max(levels)}\t0\t{i}\t{o['tname']}\t{o['strand']}\t"
                                  f"intact\tintact\t\t\t\n")
                        continue
                    if not rec["ok"]:
                        out.write(f"{qname}\t{max(levels)}\t{n_split}\t{i}\t{o['tname']}\t{o['strand']}\t"
                                  f"NA\tNA\tNA\tcontig_not_found\t\n")
                        continue
                    par = ",".join(map(str, rec["parent"]))   # all parents, deepest chain first
                    out.write(
                        f"{qname}\t{max(levels)}\t{n_split}\t{i}\t{o['tname']}\t{o['strand']}\t"
                        + ",".join(str(x["len"]) for x in rec["ins"]) + "\t"
                        + ",".join(str(x["cluster"]) for x in rec["ins"]) + "\t"
                        + ",".join(str(x["pos"]) for x in rec["ins"]) + "\t" + par + "\t"
                        + ";".join(x["seq"][:40] + ("..." if len(x["seq"]) > 40 else "")
                                   for x in rec["ins"]) + "\n")

    # Ladder summary. `supported_depth` = longest chain of observations in
    # which each contains all insertions of the one before plus one, over ALL
    # compatible parents; `supported_edges` = all such parent-child pairs (only for
    # the 30 deepest candidates, and only among the first --max-ins-examples
    # observations per level; NA = not examined).
    ladder_path = f"{args.output_prefix}_ladders_summary.tsv"
    with open(ladder_path, "w") as f:
        f.write("query\tmax_splits\tn_split_levels\tobs_per_level\ttotal_observations\t"
                "supported_depth\tsupported_edges\n")
        for qname, levels in order:
            obs_counts = ",".join(f"V{k}:{len(v)}" for k, v in levels.items())
            total = sum(len(v) for v in levels.values())
            d, e = support.get(qname, ("NA", "NA"))
            f.write(f"{qname}\t{max(levels.keys())}\t{len(levels)}\t{obs_counts}\t{total}\t{d}\t{e}\n")
    print(f"\nLadder summary: {ladder_path}")
    n_sup = sum(1 for d, _ in support.values() if d >= args.min_ladder_depth)
    print(f"Candidate ladders examined: {len(support)}; with a supported chain of "
          f">= {args.min_ladder_depth} levels: {n_sup}")
    print(f"Nesting details: {nesting_path}")


if __name__ == "__main__":
    main()
