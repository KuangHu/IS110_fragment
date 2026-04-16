#!/usr/bin/env python3
"""Detect split IS110 transposases across ALL contigs (not restricted to same contig)."""

import os
from collections import defaultdict

paf = "/groups/rubin/projects/kuang/out/IS110_fragment/hmm_search/is110_vs_gtdb.paf"
outdir = "/groups/rubin/projects/kuang/out/IS110_fragment/hmm_search"
MAX_GAP = 20000
MIN_IDENT = 80
MIN_COV = 50  # minimum combined query coverage %

# Parse all hits, group by query only
hits_by_query = defaultdict(list)
with open(paf) as fh:
    for line in fh:
        c = line.strip().split("\t")
        if len(c) < 12:
            continue
        matches, block = int(c[9]), int(c[10])
        ident = matches / block * 100 if block > 0 else 0
        if ident < MIN_IDENT:
            continue
        hits_by_query[c[0]].append({
            "qlen": int(c[1]),
            "qs": int(c[2]),
            "qe": int(c[3]),
            "tname": c[5],
            "tlen": int(c[6]),
            "ts": int(c[7]),
            "te": int(c[8]),
            "strand": c[4],
            "ident": ident,
        })

print(f"Queries with hits: {len(hits_by_query)}")

# For each query, find partial hits that together cover the query but are on different genomic locations
results = []
for qname, hlist in hits_by_query.items():
    qlen = hlist[0]["qlen"]

    # Look at all pairs of hits
    for i in range(len(hlist)):
        for j in range(i + 1, len(hlist)):
            hi, hj = hlist[i], hlist[j]

            # Check query overlap — the two parts should cover different portions
            qov = max(0, min(hi["qe"], hj["qe"]) - max(hi["qs"], hj["qs"]))
            min_span = min(hi["qe"] - hi["qs"], hj["qe"] - hj["qs"])
            if min_span == 0:
                continue
            if qov > min_span * 0.5:
                continue

            # Combined query coverage
            covered = set(range(hi["qs"], hi["qe"])) | set(range(hj["qs"], hj["qe"]))
            cov_pct = len(covered) / qlen * 100
            if cov_pct < MIN_COV:
                continue

            same_contig = hi["tname"] == hj["tname"]

            if same_contig:
                # Same contig: gap must be <= 20kb
                gap = abs(hj["ts"] - hi["te"]) if hj["ts"] > hi["te"] else abs(hi["ts"] - hj["te"])
                if gap > MAX_GAP:
                    continue
                if gap < 50:  # too close = not a real split
                    continue
                split_type = "same_contig"
            else:
                # Different contigs: can't measure gap, skip
                continue

            results.append({
                "qname": qname,
                "qlen": qlen,
                "tname1": hi["tname"],
                "tname2": hj["tname"],
                "hi": hi,
                "hj": hj,
                "gap": gap,
                "cov_pct": cov_pct,
                "split_type": split_type,
            })

# Deduplicate: for each query, keep the best hit (highest coverage)
best_per_query = {}
for r in results:
    key = r["qname"]
    if key not in best_per_query or r["cov_pct"] > best_per_query[key]["cov_pct"]:
        best_per_query[key] = r

results_dedup = sorted(best_per_query.values(), key=lambda x: -x["cov_pct"])

print(f"\nTotal split IS110 cases (gap 50bp-{MAX_GAP}bp): {len(results_dedup)}")
print()

# Print top results
for r in results_dedup[:50]:
    hi, hj = r["hi"], r["hj"]
    print(f"  {r['qname']} (qlen={r['qlen']})")
    if r["split_type"] == "same_contig":
        print(f"    SAME CONTIG: {hi['tname']} ({hi['tlen']:,}bp)")
        print(f"    Part1: Q[{hi['qs']}-{hi['qe']}] -> T[{hi['ts']:,}-{hi['te']:,}] {hi['strand']} {hi['ident']:.0f}%")
        print(f"    Part2: Q[{hj['qs']}-{hj['qe']}] -> T[{hj['ts']:,}-{hj['te']:,}] {hj['strand']} {hj['ident']:.0f}%")
        print(f"    Gap: {r['gap']:,}bp | Coverage: {r['cov_pct']:.0f}%")
    else:
        print(f"    DIFF CONTIG: {hi['tname']} + {hj['tname']}")
        print(f"    Part1: Q[{hi['qs']}-{hi['qe']}] -> {hi['tname']}:[{hi['ts']:,}-{hi['te']:,}] {hi['strand']} {hi['ident']:.0f}%")
        print(f"    Part2: Q[{hj['qs']}-{hj['qe']}] -> {hj['tname']}:[{hj['ts']:,}-{hj['te']:,}] {hj['strand']} {hj['ident']:.0f}%")
        print(f"    Coverage: {r['cov_pct']:.0f}%")
    print()

# Write full results
out_tsv = os.path.join(outdir, "split_is110_all.tsv")
with open(out_tsv, "w") as f:
    f.write("query\tqlen\tsplit_type\tcontig1\ttlen1\tpart1_qs\tpart1_qe\tpart1_ts\tpart1_te\tpart1_strand\tpart1_ident\tcontig2\ttlen2\tpart2_qs\tpart2_qe\tpart2_ts\tpart2_te\tpart2_strand\tpart2_ident\tgap\tcds_coverage\n")
    for r in results_dedup:
        hi, hj = r["hi"], r["hj"]
        f.write(f"{r['qname']}\t{r['qlen']}\t{r['split_type']}\t"
                f"{hi['tname']}\t{hi['tlen']}\t{hi['qs']}\t{hi['qe']}\t{hi['ts']}\t{hi['te']}\t{hi['strand']}\t{hi['ident']:.1f}\t"
                f"{hj['tname']}\t{hj['tlen']}\t{hj['qs']}\t{hj['qe']}\t{hj['ts']}\t{hj['te']}\t{hj['strand']}\t{hj['ident']:.1f}\t"
                f"{r['gap']}\t{r['cov_pct']:.1f}\n")

print(f"Results written to {out_tsv}")
