#!/usr/bin/env python3
"""Find whether IS110 insertions have a source copy elsewhere in the same genome.

For each insertion event [IS110_part1]--[inserted_DNA]--[IS110_part2]:
  1. Extract the inserted sequence from the target contig
  2. Search it against all contigs of the SAME Logan sample/accession
  3. Exclude:
     - Self-hit (same contig, same position)
     - Hits that land inside another known insertion region in the same sample
       (these are co-occurring [IS110][X][IS110] events, not sources)
  4. Remaining hits = candidate autologous source copies

Usage:
    python find_insertion_source.py <insertion_tsv> <logan_dir> <output_tsv>
"""

import argparse
import csv
import os
import sys
import subprocess
import tempfile
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paf_utils import aligned_identity, aligned_query_bases  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("insertion_tsv", help="TSV of detected insertions (from detect_is110_variants.py)")
    p.add_argument("logan_dir", help="Directory of Logan .contigs.fa.zst files")
    p.add_argument("output_tsv", help="Output TSV with source annotations")
    p.add_argument("--min-identity", type=float, default=95,
                   help="Min identity %% for source hits (default: 95)")
    p.add_argument("--min-coverage", type=float, default=80,
                   help="Min query coverage %% for source hits (default: 80)")
    p.add_argument("--min-ins-len", type=int, default=200,
                   help="Min inserted length to bother searching (default: 200)")
    p.add_argument("--threads", type=int, default=4,
                   help="Threads for minimap2 (default: 4)")
    p.add_argument("--limit", type=int, default=0,
                   help="Only process first N insertions (for testing; 0 = all)")
    return p.parse_args()


def accession_from_contig(contig_name):
    """Extract accession from Logan contig name: SRR12345_126 -> SRR12345."""
    return contig_name.rsplit("_", 1)[0]


def load_insertions(tsv_path):
    """Load insertion TSV and group by sample accession."""
    insertions = []
    with open(tsv_path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r.get("type") != "insertion":
                continue
            # Compute inserted region bounds on target
            p1_ts, p1_te = int(r["part1_ts"]), int(r["part1_te"])
            p2_ts, p2_te = int(r["part2_ts"]), int(r["part2_te"])
            if p1_ts < p2_ts:
                ins_start, ins_end = p1_te, p2_ts
            else:
                ins_start, ins_end = p2_te, p1_ts
            if ins_end <= ins_start:
                continue
            ins_len = ins_end - ins_start
            insertions.append({
                "query": r["query"],
                "tname": r["target_contig"],
                "ins_start": ins_start,
                "ins_end": ins_end,
                "ins_len": ins_len,
                "sample": accession_from_contig(r["target_contig"]),
            })
    return insertions


def load_sample_contigs(zst_path):
    """Decompress zst and return dict of {contig_name: sequence}."""
    result = subprocess.run(
        ["zstd", "-dcq", zst_path],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        return {}
    contigs = {}
    name = None
    seq_parts = []
    for line in result.stdout.split("\n"):
        if line.startswith(">"):
            if name and seq_parts:
                contigs[name] = "".join(seq_parts)
            name = line[1:].split()[0]
            seq_parts = []
        elif line:
            seq_parts.append(line.strip())
    if name and seq_parts:
        contigs[name] = "".join(seq_parts)
    return contigs


def write_fasta(path, seqs):
    """Write dict of {name: seq} to fasta."""
    with open(path, "w") as f:
        for name, seq in seqs.items():
            f.write(f">{name}\n{seq}\n")


def run_minimap2(query_fa, target_fa, threads):
    """Run minimap2, return list of (qname, qlen, qs, qe, tname, tlen, ts, te, strand, ident)."""
    result = subprocess.run(
        ["minimap2", "-x", "asm10", "-c", "--eqx", "--secondary=yes",
         "-N", "100000", "-p", "0.05", "-t", str(threads),
         target_fa, query_fa],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        # a killed/failed minimap2 used to read as "No source found" for all
        raise RuntimeError(f"minimap2 failed ({result.returncode}): {result.stderr[-500:]}")
    hits = []
    for line in result.stdout.split("\n"):
        c = line.strip().split("\t")
        if len(c) < 12:
            continue
        ident = aligned_identity(c)   # gap columns excluded (paf_utils)
        hits.append({
            "qname": c[0],
            "qlen": int(c[1]),
            "qs": int(c[2]),
            "qe": int(c[3]),
            "strand": c[4],
            "tname": c[5],
            "tlen": int(c[6]),
            "ts": int(c[7]),
            "te": int(c[8]),
            "ident": ident,
            # query bases really aligned: the span qe-qs also contains the
            # parts of the insert that are MISSING from the candidate source
            "qaln": aligned_query_bases(c),
        })
    return hits


def main():
    args = parse_args()

    print(f"Loading insertions from {args.insertion_tsv}")
    insertions = load_insertions(args.insertion_tsv)
    print(f"  Loaded {len(insertions)} insertions")

    # Filter by min length
    insertions = [i for i in insertions if i["ins_len"] >= args.min_ins_len]
    print(f"  After filtering (>= {args.min_ins_len} bp): {len(insertions)}")

    if args.limit > 0:
        insertions = insertions[:args.limit]
        print(f"  Limited to first {len(insertions)} for testing")

    # Group by sample
    by_sample = defaultdict(list)
    for ins in insertions:
        by_sample[ins["sample"]].append(ins)

    print(f"  Spread across {len(by_sample)} samples")
    print()

    # Build per-sample forbidden regions (all known insertion regions in that sample)
    # Format: {sample: {contig: [(ins_start, ins_end), ...]}}
    forbidden = defaultdict(lambda: defaultdict(list))
    for ins in insertions:
        forbidden[ins["sample"]][ins["tname"]].append((ins["ins_start"], ins["ins_end"]))

    # Process each sample
    results = []
    tmpdir = tempfile.mkdtemp(prefix="is110_src_")
    query_fa = os.path.join(tmpdir, "query.fa")
    target_fa = os.path.join(tmpdir, "target.fa")

    for i, (sample, sample_inss) in enumerate(by_sample.items(), 1):
        zst = os.path.join(args.logan_dir, f"{sample}.contigs.fa.zst")
        if not os.path.exists(zst):
            for ins in sample_inss:
                results.append({**ins, "n_src_hits": 0, "src_hits": "no_sample_file"})
            continue

        contigs = load_sample_contigs(zst)
        if not contigs:
            for ins in sample_inss:
                results.append({**ins, "n_src_hits": 0, "src_hits": "empty_sample"})
            continue

        # Write sample contigs as target db
        write_fasta(target_fa, contigs)

        # Build queries: all inserted sequences from this sample
        queries = {}
        ins_id_map = {}
        for j, ins in enumerate(sample_inss):
            # every input locus gets an output row (these two cases used to
            # vanish from the table and from the percentages)
            src_contig = contigs.get(ins["tname"])
            if not src_contig:
                results.append({**ins, "n_src_hits": 0, "src_hits": "no_contig"})
                continue
            ins_seq = src_contig[ins["ins_start"]:ins["ins_end"]]
            if len(ins_seq) < args.min_ins_len:
                results.append({**ins, "n_src_hits": 0, "src_hits": "no_seq"})
                continue
            qid = f"ins_{j}"
            queries[qid] = ins_seq
            ins_id_map[qid] = ins

        if not queries:
            continue

        write_fasta(query_fa, queries)

        # Run minimap2
        hits = run_minimap2(query_fa, target_fa, args.threads)

        # Process hits per insertion
        by_query = defaultdict(list)
        for h in hits:
            if h["ident"] < args.min_identity:
                continue
            # coverage from aligned bases (a source with only 70% of the
            # insert, bridged in one gapped row, was "100% coverage")
            cov = h["qaln"] / h["qlen"] * 100 if h["qlen"] > 0 else 0
            if cov < args.min_coverage:
                continue
            by_query[h["qname"]].append(h)

        for qid, ins in ins_id_map.items():
            qhits = by_query.get(qid, [])
            kept = []
            for h in qhits:
                # Exclude self-hit: same contig, overlapping the insertion region
                if h["tname"] == ins["tname"]:
                    # Check overlap with the insertion region on the same contig
                    ov = max(0, min(h["te"], ins["ins_end"]) - max(h["ts"], ins["ins_start"]))
                    span = h["te"] - h["ts"]
                    if span > 0 and ov / span > 0.5:
                        continue  # self

                # Exclude hits that land inside ANY known insertion region in this sample
                skip = False
                for fs, fe in forbidden[ins["sample"]].get(h["tname"], []):
                    ov = max(0, min(h["te"], fe) - max(h["ts"], fs))
                    span = h["te"] - h["ts"]
                    if span > 0 and ov / span > 0.5:
                        skip = True
                        break
                if skip:
                    continue

                kept.append(h)

            src_str = ";".join(
                f"{h['tname']}:{h['ts']}-{h['te']}({h['strand']},{h['ident']:.0f}%)"
                for h in kept[:5]
            )
            if len(kept) > 5:
                src_str += f";...(+{len(kept)-5})"
            results.append({
                **ins,
                "n_src_hits": len(kept),
                "src_hits": src_str if kept else "none",
            })

        if i % 50 == 0 or i == len(by_sample):
            print(f"  Processed {i}/{len(by_sample)} samples, {len(results)} insertions done")

    # Cleanup
    for f in [query_fa, target_fa]:
        if os.path.exists(f):
            os.remove(f)
    os.rmdir(tmpdir)

    # Write output
    with open(args.output_tsv, "w") as f:
        f.write("query\tsample\ttarget_contig\tins_start\tins_end\tins_len\tn_src_hits\tsrc_hits\n")
        for r in results:
            f.write(f"{r['query']}\t{r['sample']}\t{r['tname']}\t"
                    f"{r['ins_start']}\t{r['ins_end']}\t{r['ins_len']}\t"
                    f"{r['n_src_hits']}\t{r['src_hits']}\n")

    # Summary
    n_total = len(results)
    n_with_src = sum(1 for r in results if r["n_src_hits"] > 0)
    n_none = sum(1 for r in results if r["src_hits"] == "none")
    n_errors = sum(1 for r in results
                   if r["src_hits"] in ("no_sample_file", "empty_sample", "no_seq", "no_contig"))
    assert n_total == len(insertions), (n_total, len(insertions))

    print()
    print("=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"Total insertions processed: {n_total:,}")
    denom = max(1, n_total - n_errors)      # percentages over assessable loci
    print(f"With autologous source:     {n_with_src:,} ({n_with_src/denom*100:.1f}% of assessable)")
    print(f"No source found:            {n_none:,} ({n_none/denom*100:.1f}% of assessable)")
    print(f"Errors/missing:             {n_errors:,}")
    print(f"\nOutput: {args.output_tsv}")


if __name__ == "__main__":
    main()
