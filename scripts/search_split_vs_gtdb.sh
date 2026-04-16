#!/bin/bash
#SBATCH --job-name=split_vs_gtdb
#SBATCH --partition=standard
#SBATCH --cpus-per-task=48
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/groups/rubin/projects/kuang/out/IS110_fragment/hmm_search/split_vs_gtdb_%j.out
#SBATCH --error=/groups/rubin/projects/kuang/out/IS110_fragment/hmm_search/split_vs_gtdb_%j.err

set -euo pipefail

OUTDIR="/groups/rubin/projects/kuang/out/IS110_fragment/hmm_search"
GTDB_FNA="/groups/rubin/databases/GTDB/gtdb_genomes_v2024latest/gtdb_genomes_reps_r220/allfilestogether"
CPUS=48

echo "=== Search 123 split IS110 CDS vs all GTDB ==="
echo "Start: $(date)"

# Step 1: Extract the 123 query sequences
echo "Step 1: Extracting split IS110 query sequences..."
tail -n +2 "$OUTDIR/split_is110_all.tsv" | cut -f1 | sort -u > "$OUTDIR/split_query_ids.txt"

python3 -c "
from Bio import SeqIO
ids = set(open('$OUTDIR/split_query_ids.txt').read().split())
count = 0
with open('$OUTDIR/split_123_queries.fna', 'w') as out:
    for rec in SeqIO.parse('$OUTDIR/is110_cds_nr95.fna', 'fasta'):
        if rec.id in ids:
            SeqIO.write(rec, out, 'fasta')
            count += 1
# some queries may be in the non-dereplicated file
if count < len(ids):
    found = set(r.id for r in SeqIO.parse('$OUTDIR/split_123_queries.fna', 'fasta'))
    missing = ids - found
    for rec in SeqIO.parse('$OUTDIR/is110_cds_nuc.fna', 'fasta'):
        if rec.id in missing:
            with open('$OUTDIR/split_123_queries.fna', 'a') as out:
                SeqIO.write(rec, out, 'fasta')
            count += 1
print(f'  Extracted {count} / {len(ids)} query sequences')
"

# Step 2: minimap2 search against all GTDB
echo "Step 2: Searching against all GTDB genomes..."

find "$GTDB_FNA" -name "*.fna" -print0 | xargs -0 cat \
    | minimap2 -x asm10 -c --eqx --secondary=yes -N 1000 -p 0.3 \
        -t "$CPUS" \
        /dev/stdin "$OUTDIR/split_123_queries.fna" \
        -o "$OUTDIR/split_123_vs_gtdb.paf"

echo "  Total hits: $(wc -l < "$OUTDIR/split_123_vs_gtdb.paf")"
echo "Step 2 done: $(date)"

# Step 3: Detect split cases across all contigs
echo "Step 3: Detecting split cases..."

python3 - << 'PYEOF'
import os
from collections import defaultdict

paf = "/groups/rubin/projects/kuang/out/IS110_fragment/hmm_search/split_123_vs_gtdb.paf"
outdir = "/groups/rubin/projects/kuang/out/IS110_fragment/hmm_search"
MAX_GAP = 20000
MIN_IDENT = 80
MIN_COV = 50

hits_by_query = defaultdict(list)
with open(paf) as fh:
    for line in fh:
        c = line.strip().split("\t")
        if len(c) < 12: continue
        matches, block = int(c[9]), int(c[10])
        ident = matches / block * 100 if block > 0 else 0
        if ident < MIN_IDENT: continue
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

results = []
for qname, hlist in hits_by_query.items():
    qlen = hlist[0]["qlen"]
    # Group by target contig
    by_contig = defaultdict(list)
    for h in hlist:
        by_contig[h["tname"]].append(h)

    for tname, ths in by_contig.items():
        if len(ths) < 2: continue
        ths.sort(key=lambda x: x["ts"])
        for i in range(len(ths)):
            for j in range(i+1, len(ths)):
                hi, hj = ths[i], ths[j]
                gap = hj["ts"] - hi["te"]
                if gap < 50 or gap > MAX_GAP: continue
                qov = max(0, min(hi["qe"], hj["qe"]) - max(hi["qs"], hj["qs"]))
                min_span = min(hi["qe"]-hi["qs"], hj["qe"]-hj["qs"])
                if min_span == 0: continue
                if qov > min_span * 0.5: continue
                covered = set(range(hi["qs"], hi["qe"])) | set(range(hj["qs"], hj["qe"]))
                cov_pct = len(covered) / qlen * 100
                if cov_pct < MIN_COV: continue
                results.append({
                    "qname": qname, "qlen": qlen,
                    "tname": tname, "tlen": hi["tlen"],
                    "hi": hi, "hj": hj,
                    "gap": gap, "cov_pct": cov_pct,
                })

print(f"Total split hits: {len(results)}")

# Deduplicate: unique (query, target_contig)
seen = set()
dedup = []
for r in sorted(results, key=lambda x: -x["cov_pct"]):
    key = (r["qname"], r["tname"])
    if key in seen: continue
    seen.add(key)
    dedup.append(r)

# Count unique queries and targets
uq = len(set(r["qname"] for r in dedup))
ut = len(set(r["tname"] for r in dedup))
print(f"Unique query-target pairs: {len(dedup)}")
print(f"Unique queries with splits: {uq}")
print(f"Unique target contigs: {ut}")

# Write results
out_tsv = os.path.join(outdir, "split_123_all_hits.tsv")
with open(out_tsv, "w") as f:
    f.write("query\tqlen\ttarget_contig\ttlen\tpart1_qs\tpart1_qe\tpart1_ts\tpart1_te\tpart1_strand\tpart1_ident\tpart2_qs\tpart2_qe\tpart2_ts\tpart2_te\tpart2_strand\tpart2_ident\tgap\tcds_coverage\n")
    for r in dedup:
        hi, hj = r["hi"], r["hj"]
        f.write(f"{r['qname']}\t{r['qlen']}\t{r['tname']}\t{r['tlen']}\t"
                f"{hi['qs']}\t{hi['qe']}\t{hi['ts']}\t{hi['te']}\t{hi['strand']}\t{hi['ident']:.1f}\t"
                f"{hj['qs']}\t{hj['qe']}\t{hj['ts']}\t{hj['te']}\t{hj['strand']}\t{hj['ident']:.1f}\t"
                f"{r['gap']}\t{r['cov_pct']:.1f}\n")

print(f"\nResults written to {out_tsv}")

# Summary per query
print(f"\n=== Hits per query (top 20) ===")
from collections import Counter
qcounts = Counter(r["qname"] for r in dedup)
for q, c in qcounts.most_common(20):
    print(f"  {q}: {c} split locations")
PYEOF

echo ""
echo "=== Done ==="
echo "End: $(date)"
