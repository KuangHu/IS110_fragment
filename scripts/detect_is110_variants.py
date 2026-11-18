#!/usr/bin/env python3
"""Detect IS110 structural variants from minimap2 PAF alignments.

Detects two scenarios:

Scenario 1 — INSERTION in genome (genome longer than reference):
    Reference:  [========= IS110 CDS =========]
    Genome:     [IS110 part1]---gap---[IS110 part2]
    Signal: two hits cover different query portions, large GAP on TARGET,
            small or no gap on query.

Scenario 2 — DELETION in genome (genome shorter than reference):
    Reference:  [IS110 part1]---DNA---[IS110 part2]
    Genome:     [IS110 part1][IS110 part2]
    Signal: two hits cover different query portions, large GAP on QUERY,
            small or no gap on target (parts are adjacent on genome).

Usage:
    python detect_is110_variants.py <input.paf> <output_prefix> [options]
"""

import argparse
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paf_utils import aligned_identity, aligned_blocks_iv, union_length  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("paf", help="Input PAF file from minimap2")
    p.add_argument("output_prefix", help="Output prefix for result files")
    p.add_argument("--min-identity", type=float, default=95,
                   help="Min alignment identity %% (default: 95)")
    p.add_argument("--min-query-cov", type=float, default=50,
                   help="Min combined query coverage %% (default: 50)")
    # Scenario 1: insertion
    p.add_argument("--ins-min-target-gap", type=int, default=500,
                   help="Scenario 1: min NET inserted length = target gap minus "
                        "query gap (default: 500)")
    p.add_argument("--ins-max-query-gap", type=int, default=50,
                   help="Scenario 1: the two halves must meet on the IS within "
                        "this many bp (default: 50)")
    p.add_argument("--ins-max-target-gap", type=int, default=50000,
                   help="Scenario 1: max gap on target (default: 50000)")
    # Scenario 2: deletion
    p.add_argument("--del-min-query-gap", type=int, default=200,
                   help="Scenario 2: min NET deleted length = query gap minus "
                        "target gap (default: 200)")
    p.add_argument("--del-max-target-gap", type=int, default=200,
                   help="Scenario 2: max distance between parts on target (default: 200)")
    p.add_argument("--max-junction-overlap", type=int, default=100,
                   help="Both scenarios: the two parts may overlap by up to this "
                        "many bp on the IS or on the contig (default: 100)")
    p.add_argument("--max-arm-overlap-frac", type=float, default=0.5,
                   help="Both scenarios: the overlap of the two parts on the IS may "
                        "be at most this fraction of the shorter part (default: 0.5)")
    p.add_argument("--locus-tol", type=int, default=20,
                   help="Records of different reference IS whose junctions lie within "
                        "this many bp of the group's first junction are one locus "
                        "(default: 20)")
    return p.parse_args()


def parse_paf(path, min_identity, ins_min=500, del_min=200, counts=None):
    """Parse PAF, filter by identity, cut every row at its large indels.

    Returns {(query, contig, strand): [block, ...]}. A block is a stretch of
    the IS that is really aligned: a row with a 600-bp insertion gives two
    blocks, exactly like two separate rows would. Query bases inside a large
    'I' run are in no block (they are absent from the genome)."""
    blocks = defaultdict(list)
    total = kept = 0
    row_id = 0
    with open(path) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 12:
                continue
            total += 1
            ident = aligned_identity(c)        # gaps excluded (see paf_utils)
            if ident < min_identity:
                continue
            kept += 1
            row_id += 1
            for qs, qe, ts, te, ivs in aligned_blocks_iv(c, ins_min, del_min):
                blocks[(c[0], c[5], c[4])].append({
                    "qname": c[0], "qlen": int(c[1]), "qs": qs, "qe": qe,
                    "strand": c[4], "tname": c[5], "tlen": int(c[6]),
                    "ts": ts, "te": te, "ident": ident, "row": row_id,
                    "aln": ivs, "n_aln": sum(e - s for s, e in ivs),
                })
    print(f"Parsed {total} alignments, {kept} pass identity >= {min_identity}%")
    print(f"Unique queries with hits: {len(set(k[0] for k in blocks))}")
    if counts is not None:
        counts.update(total=total, kept=kept)
    return blocks


def classify_pair(a, b, args):
    """a = target-left block, b = target-right block, same query/contig/strand.

    Returns (type, target_gap, query_gap) or None. Gaps are signed: negative =
    the two alignments overlap (sequence repeated at the junction, e.g. a
    target-site duplication or junction homology; minimap2 also gives -1..-4
    at clean junctions). The event size is the NET length change, so an
    overlap does not make a 500-bp insertion "499 bp, rejected"."""
    first, second = (a, b) if a["strand"] == "+" else (b, a)   # order on the IS
    if not (first["qs"] < second["qs"] and first["qe"] < second["qe"]):
        return None
    qgap = second["qs"] - first["qe"]
    tgap = b["ts"] - a["te"]
    ov = args.max_junction_overlap
    if qgap < -ov or tgap < -ov:
        return None
    # the overlap must stay a small part of either arm
    if -qgap > getattr(args, "max_arm_overlap_frac", 0.5) * min(first["qe"] - first["qs"], second["qe"] - second["qs"]):
        return None
    net_ins = tgap - qgap            # extra bases in the genome
    net_del = qgap - tgap            # IS bases missing from the genome
    if (qgap <= args.ins_max_query_gap
            and args.ins_min_target_gap <= net_ins <= args.ins_max_target_gap):
        return ("insertion", tgap, qgap)
    if tgap <= args.del_max_target_gap and net_del >= args.del_min_query_gap:
        return ("deletion", tgap, qgap)
    return None


def drop_contained(blocks, slop=5):
    """Remove alternative alignments that lie INSIDE another alignment.

    A reference IS with an internal repeat aligns to an intact copy once in
    full and again, shifted, with each repeat unit (minimap2 reports these
    as secondary rows). The shifted pieces of two intact copies paired up
    as an "interrupted copy" spanning the DNA between them. A block whose
    target interval lies within the target interval of a block from another
    row of the same reference on the same contig (either strand) carries no
    information of its own about that stretch of the contig and is dropped.

    Returns (blocks, n_dropped)."""
    by_qc = defaultdict(list)
    for (q, contig, strand), bl in blocks.items():
        by_qc[(q, contig)] += bl
    out = defaultdict(list)
    dropped = 0
    for (q, contig) in sorted(by_qc):
        bl = sorted(by_qc[(q, contig)], key=lambda x: (x["ts"], -x["te"], x["row"], x["qs"]))
        active = []                          # earlier blocks still open on the target
        for b in bl:
            active = [a for a in active if a["te"] + slop >= b["te"] or a["te"] > b["ts"]]
            inside = False
            for a in active:
                if a["row"] == b["row"] or a["ts"] - slop > b["ts"] or a["te"] + slop < b["te"]:
                    continue
                la, lb = a["te"] - a["ts"], b["te"] - b["ts"]
                if la > lb or (la == lb and (a["n_aln"], -a["row"]) > (b["n_aln"], -b["row"])):
                    inside = True
                    break
            if inside:
                dropped += 1
                continue
            active.append(b)
            out[(q, contig, b["strand"])].append(b)
    return dict(out), dropped


def _best_paths(nodes, edges, bl):
    """Split one connected group of blocks into chains.

    Repeatedly takes the path (blocks in target order, joined by candidate
    junctions) that covers the most aligned IS bases, then continues with the
    blocks that are left. Ties: more junctions inside one alignment row, then
    the smaller total target gap, then the leftmost blocks.

    The earlier rule gave every block its NEAREST partner first: a 200-bp
    match inside an insert took the left arm, the pair failed the coverage
    filter, and the real right arm was left without a partner."""
    left = set(nodes)
    paths = []
    while left:
        order = sorted(left)
        best = {}                            # node -> (score tuple, previous node)
        for j in order:
            cur = ((bl[j]["n_aln"], 0, 0), None)
            for i, ev in edges.get(j, ()):
                if i not in left:
                    continue
                (sc, same, gap), _ = best[i]
                cand = (sc + bl[j]["n_aln"] - max(0, -ev[2]),
                        same + (bl[i]["row"] == bl[j]["row"]),
                        gap - max(0, ev[1]))
                if cand > cur[0]:
                    cur = (cand, i)
            best[j] = cur
        end = max(order, key=lambda j: (best[j][0], -j))
        path = [end]
        while best[path[-1]][1] is not None:
            path.append(best[path[-1]][1])
        path.reverse()
        paths.append(path)
        left -= set(path)
    return paths


def build_chains(blocks, args):
    """Chain blocks into loci.

    Per (query, contig, strand), blocks are taken in target order; every pair
    that forms an insertion or a deletion is a candidate junction. The blocks
    are then split into chains by `_best_paths`: each block is in exactly one
    chain, so two interrupted copies on one contig do not give a third
    "insertion" (left arm of copy 1 + right arm of copy 2).

    Returns (chains, stats). A chain = {"blocks": [...target order...],
    "links": [(a, b, (type, target_gap, query_gap)), ...], "cov_pct": float,
    "id": int}; a block with no partner is a chain of one. Coverage is the
    union of the ALIGNED query intervals of the whole chain: a copy with two
    insertions is three blocks, two events, one coverage, and small 'I' runs
    inside a block are not counted as covered."""
    blocks, n_contained = drop_contained(blocks)
    chains = []
    stats = {"competing": 0, "contained": n_contained, "replacement_like": 0}
    # a junction may leave up to ins_max_query_gap unaligned IS bases, so the
    # target gap of a net insertion of exactly the maximum is larger than it
    window = args.ins_max_target_gap + max(0, args.ins_max_query_gap)
    for key in sorted(blocks):
        bl = sorted(blocks[key], key=lambda x: (x["ts"], x["te"], x["qs"], x["row"]))
        edges = defaultdict(list)            # right block -> [(left block, event)]
        nbr = defaultdict(set)
        n_cand = 0
        for i, a in enumerate(bl):
            for j in range(i + 1, len(bl)):
                b = bl[j]
                if b["ts"] - a["te"] > window:
                    break
                ev = classify_pair(a, b, args)
                if ev:
                    edges[j].append((i, ev))
                    nbr[i].add(j)
                    nbr[j].add(i)
                    n_cand += 1
                elif is_replacement_like(a, b, args):
                    stats["replacement_like"] += 1
        seen = set()
        n_used = 0
        for start in range(len(bl)):
            if start in seen:
                continue
            comp, stack = [], [start]
            seen.add(start)
            while stack:
                k = stack.pop()
                comp.append(k)
                for m in sorted(nbr[k]):
                    if m not in seen:
                        seen.add(m)
                        stack.append(m)
            for path in _best_paths(comp, edges, bl):
                links = []
                for i, j in zip(path, path[1:]):
                    ev = next(e for k, e in edges[j] if k == i)
                    links.append((bl[i], bl[j], ev))
                n_used += len(links)
                covered = union_length([iv for m in path for iv in bl[m]["aln"]])
                chains.append({"blocks": [bl[m] for m in path], "links": links,
                               "cov_pct": covered / bl[path[0]]["qlen"] * 100})
        stats["competing"] += n_cand - n_used
    chains.sort(key=lambda ch: (ch["blocks"][0]["qname"], ch["blocks"][0]["tname"],
                                ch["blocks"][0]["ts"], ch["blocks"][0]["strand"],
                                ch["blocks"][0]["row"]))
    for n, ch in enumerate(chains):
        ch["id"] = n
    return chains, stats


def is_replacement_like(a, b, args):
    """Two colinear parts with BOTH a large gap on the IS and a large gap on
    the contig (IS sequence replaced by other DNA). Neither scenario of this
    script; such pairs are counted in the log and not reported."""
    first, second = (a, b) if a["strand"] == "+" else (b, a)
    if not (first["qs"] < second["qs"] and first["qe"] < second["qe"]):
        return False
    qgap = second["qs"] - first["qe"]
    tgap = b["ts"] - a["te"]
    return (qgap >= args.del_min_query_gap and args.ins_min_target_gap <= tgap
            and tgap > args.del_max_target_gap)


def detect_events(blocks, args):
    """The junctions of all chains that cover enough of the IS."""
    chains, stats = build_chains(blocks, args)
    results = []
    for ch in chains:
        if ch["cov_pct"] < args.min_query_cov:
            continue
        for a, b, (etype, tgap, qgap) in ch["links"]:
            results.append({
                "type": etype, "qname": a["qname"], "qlen": a["qlen"],
                "tname": a["tname"], "tlen": a["tlen"], "hi": a, "hj": b,
                "target_gap": tgap, "query_gap": qgap, "cov_pct": ch["cov_pct"],
                "chain": ch["id"],
            })
    return results, stats


def junction(r):
    a, b = sorted((r["hi"], r["hj"]), key=lambda h: h["ts"])
    return a["te"], b["ts"]


def deduplicate(results, tol=20):
    """One row per PHYSICAL JUNCTION: the same junction is reported once for
    every near-identical reference IS, and those records are merged; the one
    with the best coverage is kept.

    Records are taken in order of their junction position on the contig. A
    record joins the first group whose FIRST junction lies within `tol` bp on
    both sides, unless that group already holds a junction of the same chain:
    two junctions of one copy (`500=200I10=200I500=`: deletions at 500 and
    510) are two events however close they are. Measuring from the first
    member bounds a group to `tol` bp; chaining neighbours without a bound
    merged junctions 40 bp apart through a third record between them.
    The result does not depend on the order of the input."""
    by_contig = defaultdict(list)
    for r in results:
        by_contig[r["tname"]].append(r)
    kept = []
    for tname in sorted(by_contig):
        rs = sorted(by_contig[tname],
                    key=lambda r: (junction(r), -r["cov_pct"], r["qname"], r["type"], r["chain"]))
        groups = []
        for r in rs:
            j = junction(r)
            for g in groups:
                j0 = junction(g[0])
                if (abs(j0[0] - j[0]) <= tol and abs(j0[1] - j[1]) <= tol
                        and all(k["chain"] != r["chain"] for k in g)):
                    g.append(r)
                    break
            else:
                groups.append([r])
        for g in groups:
            kept.append(min(g, key=lambda x: (-x["cov_pct"], x["qname"], x["type"],
                                              junction(x))))
    return kept


def write_tsv(results, path):
    """Write results to TSV."""
    with open(path, "w") as f:
        f.write(
            "type\tquery\tqlen\ttarget_contig\ttlen\t"
            "part1_qs\tpart1_qe\tpart1_ts\tpart1_te\tpart1_strand\tpart1_ident\t"
            "part2_qs\tpart2_qe\tpart2_ts\tpart2_te\tpart2_strand\tpart2_ident\t"
            "target_gap\tquery_gap\tcds_coverage\n"
        )
        for r in results:
            hi, hj = r["hi"], r["hj"]
            f.write(
                f"{r['type']}\t{r['qname']}\t{r['qlen']}\t{r['tname']}\t{r['tlen']}\t"
                f"{hi['qs']}\t{hi['qe']}\t{hi['ts']}\t{hi['te']}\t{hi['strand']}\t{hi['ident']:.1f}\t"
                f"{hj['qs']}\t{hj['qe']}\t{hj['ts']}\t{hj['te']}\t{hj['strand']}\t{hj['ident']:.1f}\t"
                f"{r['target_gap']}\t{r['query_gap']}\t{r['cov_pct']:.1f}\n"
            )


def print_summary(results, label, n=20):
    """Print top results."""
    subset = [r for r in results if r["type"] == label]
    if not subset:
        print(f"  No {label} cases found.")
        return

    # Unique queries and targets
    uq = len(set(r["qname"] for r in subset))
    ut = len(set(r["tname"] for r in subset))
    print(f"  {len(subset)} cases ({uq} unique queries, {ut} unique targets)")

    for r in subset[:n]:
        hi, hj = r["hi"], r["hj"]
        print(f"\n  {r['qname']} (qlen={r['qlen']})")
        print(f"    Target: {r['tname']} ({r['tlen']:,} bp)")
        print(f"    Part1: Q[{hi['qs']}-{hi['qe']}] -> T[{hi['ts']:,}-{hi['te']:,}] {hi['strand']} {hi['ident']:.0f}%")
        print(f"    Part2: Q[{hj['qs']}-{hj['qe']}] -> T[{hj['ts']:,}-{hj['te']:,}] {hj['strand']} {hj['ident']:.0f}%")
        if label == "insertion":
            print(f"    Target gap: {r['target_gap']:,} bp | Query coverage: {r['cov_pct']:.0f}%")
        else:
            print(f"    Query gap: {r['query_gap']:,} bp (deleted) | Target gap: {r['target_gap']} bp | Query coverage: {r['cov_pct']:.0f}%")


def main():
    args = parse_args()

    print("=" * 60)
    print("IS110 Structural Variant Detection")
    print("=" * 60)
    print(f"Input:  {args.paf}")
    print(f"Output: {args.output_prefix}_*.tsv")
    print()

    blocks = parse_paf(args.paf, args.min_identity,
                       args.ins_min_target_gap, args.del_min_query_gap)
    print()

    events, stats = detect_events(blocks, args)
    n_before = len(events)
    events = deduplicate(events, args.locus_tol)
    print(f"Events: {n_before} reference x locus records -> {len(events)} physical loci "
          f"(within {args.locus_tol} bp)")
    if stats["contained"]:
        print(f"  {stats['contained']} alignment pieces lay inside another alignment of the "
              "same reference (alternative placements of a repeat) and were not used")
    if stats["competing"]:
        print(f"  {stats['competing']} further pairings were possible but used a block "
              "that is part of a better-covered chain (not reported)")
    if stats["replacement_like"]:
        print(f"  {stats['replacement_like']} pairs have a large gap on the IS AND on the "
              "contig (replacement); neither scenario, not reported")
    print()

    print("--- Scenario 1: INSERTION in genome (target gap) ---")
    print(f"    Net inserted length: {args.ins_min_target_gap}-{args.ins_max_target_gap} bp")
    insertions = [r for r in events if r["type"] == "insertion"]
    print_summary(insertions, "insertion")
    print()

    print("--- Scenario 2: DELETION in genome (query gap) ---")
    print(f"    Net deleted length min: {args.del_min_query_gap} bp")
    print(f"    Target gap max: {args.del_max_target_gap} bp")
    deletions = [r for r in events if r["type"] == "deletion"]
    print_summary(deletions, "deletion")
    print()

    # Write outputs
    # Always written, header-only when empty: a file left by an earlier run
    # was otherwise read by the next step as this run's result.
    all_results = insertions + deletions
    for label, rows, suffix in (("Combined results:", all_results, "all"),
                                ("Insertions:      ", insertions, "insertions"),
                                ("Deletions:       ", deletions, "deletions")):
        out_path = f"{args.output_prefix}_{suffix}.tsv"
        write_tsv(rows, out_path)
        print(f"{label} {out_path} ({len(rows)} rows)")

    print()
    print(f"Total: {len(insertions)} insertions + {len(deletions)} deletions = {len(all_results)} variants")


if __name__ == "__main__":
    main()
