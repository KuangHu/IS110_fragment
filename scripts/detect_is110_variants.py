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
from paf_utils import aligned_identity, aligned_blocks  # noqa: E402


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
    p.add_argument("--locus-tol", type=int, default=20,
                   help="Records whose junctions lie within this many bp are one "
                        "locus (default: 20)")
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
            for qs, qe, ts, te in aligned_blocks(c, ins_min, del_min):
                blocks[(c[0], c[5], c[4])].append({
                    "qname": c[0], "qlen": int(c[1]), "qs": qs, "qe": qe,
                    "strand": c[4], "tname": c[5], "tlen": int(c[6]),
                    "ts": ts, "te": te, "ident": ident, "row": row_id,
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
    if -qgap > 0.5 * min(first["qe"] - first["qs"], second["qe"] - second["qs"]):
        return None
    net_ins = tgap - qgap            # extra bases in the genome
    net_del = qgap - tgap            # IS bases missing from the genome
    if (qgap <= args.ins_max_query_gap
            and args.ins_min_target_gap <= net_ins <= args.ins_max_target_gap):
        return ("insertion", tgap, qgap)
    if tgap <= args.del_max_target_gap and net_del >= args.del_min_query_gap:
        return ("deletion", tgap, qgap)
    return None


def build_chains(blocks, args):
    """Chain blocks into loci.

    Per (query, contig, strand), blocks are taken in target order. Each block
    gets at most ONE right-hand partner and one left-hand partner: the nearest
    block on the target with which it forms an insertion or a deletion (two
    blocks cut from the same row first). Every permissible pair used to be
    reported, so two interrupted copies on one contig gave a third
    "insertion" (left arm of copy 1 + right arm of copy 2).

    Returns (chains, stats). A chain = {"blocks": [...target order...],
    "links": [(a, b, (type, target_gap, query_gap)), ...], "cov_pct": float};
    a block with no partner is a chain of one. Coverage is the union of the
    aligned query intervals of the whole chain (a copy with two insertions is
    three blocks, two events, one coverage)."""
    chains = []
    stats = {"competing": 0}
    for key in sorted(blocks):
        bl = sorted(blocks[key], key=lambda x: (x["ts"], x["te"]))
        cand = []
        for i, a in enumerate(bl):
            for j in range(i + 1, len(bl)):
                b = bl[j]
                if b["ts"] - a["te"] > args.ins_max_target_gap:
                    break
                ev = classify_pair(a, b, args)
                if ev:
                    cand.append((a["row"] != b["row"], max(0, ev[1]), i, j, ev))
        cand.sort(key=lambda x: x[:4])
        right_of, left_of, link_ev = {}, {}, {}
        for _, _, i, j, ev in cand:
            if i in right_of or j in left_of:
                stats["competing"] += 1
                continue
            right_of[i] = j
            left_of[j] = i
            link_ev[i] = ev
        for i in range(len(bl)):
            if i in left_of:
                continue
            k, members, links = i, [i], []
            while k in right_of:
                links.append((bl[k], bl[right_of[k]], link_ev[k]))
                k = right_of[k]
                members.append(k)
            iv = sorted((bl[m]["qs"], bl[m]["qe"]) for m in members)
            covered, end = 0, 0
            for qs, qe in iv:
                covered += max(0, qe - max(qs, end))
                end = max(end, qe)
            chains.append({"blocks": [bl[m] for m in members], "links": links,
                           "cov_pct": covered / bl[i]["qlen"] * 100})
    return chains, stats


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
            })
    return results, stats


def deduplicate(results, tol=20):
    """One row per PHYSICAL LOCUS: records on the same contig whose two
    junction positions lie within `tol` bp of another member of the group
    (single linkage, so offsets 0 / 20 / 21 are one locus), whatever the
    reference IS and whatever the event type it suggests. The record with the
    best coverage is kept.

    The key used to be (type, query, contig): the same interrupted copy was
    reported once for every near-identical reference IS, while a SECOND event
    on the same contig was dropped."""
    def junction(r):
        a, b = sorted((r["hi"], r["hj"]), key=lambda h: h["ts"])
        return a["te"], b["ts"]
    by_contig = defaultdict(list)
    for r in results:
        by_contig[r["tname"]].append(r)
    kept = []
    for tname in sorted(by_contig):
        rs = sorted(by_contig[tname], key=junction)
        groups = []
        for r in rs:
            j = junction(r)
            for g in groups:
                if any(abs(junction(k)[0] - j[0]) <= tol and abs(junction(k)[1] - j[1]) <= tol
                       for k in g):
                    g.append(r)
                    break
            else:
                groups.append([r])
        # single linkage: join groups that a later record connects
        merged = True
        while merged:
            merged = False
            for x in range(len(groups)):
                for y in range(x + 1, len(groups)):
                    if any(abs(junction(k)[0] - junction(m)[0]) <= tol
                           and abs(junction(k)[1] - junction(m)[1]) <= tol
                           for k in groups[x] for m in groups[y]):
                        groups[x] += groups.pop(y)
                        merged = True
                        break
                if merged:
                    break
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
    if stats["competing"]:
        print(f"  {stats['competing']} further pairings were possible but used a block "
              "that already had a nearer partner (not reported)")
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
