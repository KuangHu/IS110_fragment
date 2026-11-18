#!/bin/bash
#SBATCH --job-name=cds_split_search
#SBATCH --partition=standard
#SBATCH --cpus-per-task=48
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=cds_minimap2_%j.out
#SBATCH --error=cds_minimap2_%j.err

# Search individual IS110 CDS sequences against GTDB to find split CDS cases.
# Uses asm10 preset (more sensitive for shorter/divergent sequences).

set -euo pipefail

QUERY="/home/kuangh/tools/IS110_fragment/results/is110_cds_queries.fna"
GTDB_DIR="/groups/rubin/databases/GTDB/gtdb_genomes_v2024latest/gtdb_genomes_reps_r220/allfilestogether"
OUTDIR="/home/kuangh/tools/IS110_fragment/results/cds_search"
CPUS=48
MIN_IDENT=0.80

mkdir -p "$OUTDIR"

echo "=== CDS-level minimap2 search ==="
echo "Query:    $QUERY ($(grep -c '^>' "$QUERY") CDS sequences)"
echo "Database: $GTDB_DIR"
echo "Output:   $OUTDIR"
echo "CPUs:     $CPUS"
echo "Preset:   asm10 (more sensitive for shorter sequences)"
echo "Start:    $(date)"
echo ""

# Run minimap2 with asm10 preset, report all secondary/supplementary alignments
# -N/-p: every copy, not the 100 best scoring >= half of the best one -- each
# half of an interrupted copy is a low-scoring secondary of any intact copy and
# was never reported (0 of 18 test configurations with an intact copy present).
# -g/-r 200: do not chain the two halves across the insert into one alignment.
find "$GTDB_DIR" -name "*.fna" -exec cat {} + \
    | minimap2 -x asm10 -c --eqx --secondary=yes -N 1000000 -p 0.05 -g 200 -r 200,200 \
        -t "$CPUS" \
        /dev/stdin "$QUERY" \
        -o "$OUTDIR/all_cds_hits.paf"

echo ""
echo "End: $(date)"
echo "Total raw hits: $(wc -l < "$OUTDIR/all_cds_hits.paf")"
echo ""

# Filter by identity and split per CDS query
echo "Splitting results per CDS..."
while IFS= read -r header; do
    CDS_NAME=$(echo "$header" | sed 's/>//' | awk '{print $1}')
    HITS_FILE="$OUTDIR/${CDS_NAME}_hits.paf"
    
    awk -F'\t' -v qid="$CDS_NAME" -v min_id="$MIN_IDENT" \
        '$1 == qid && $10/$11 >= min_id' "$OUTDIR/all_cds_hits.paf" > "$HITS_FILE"
    
    N=$(wc -l < "$HITS_FILE")
    echo "  ${CDS_NAME}: ${N} hits"
done < <(grep '^>' "$QUERY")

echo ""
echo "Done. Per-CDS results in: $OUTDIR/<cds_name>_hits.paf"

# Quick analysis: find split CDS cases
echo ""
echo "=== Scanning for split CDS (same CDS -> distant spots on same contig) ==="
python3 - << 'PYEND'
import os
from collections import defaultdict

outdir = "/home/kuangh/tools/IS110_fragment/results/cds_search"
MIN_GAP = 6000

paf = os.path.join(outdir, "all_cds_hits.paf")
if not os.path.exists(paf):
    print("No PAF file found")
    exit()

# Group by (query_cds, target_contig)
hits = defaultdict(list)
with open(paf) as fh:
    for line in fh:
        c = line.strip().split("\t")
        if len(c) < 12: continue
        qname = c[0]
        qlen = int(c[1])
        qs, qe = int(c[2]), int(c[3])
        tname = c[5]
        tlen = int(c[6])
        ts, te = int(c[7]), int(c[8])
        matches, block = int(c[9]), int(c[10])
        # identity over aligned columns: PAF column 11 counts gap columns, so a
        # copy with an indel fell below the threshold and was dropped
        _cg = next((t[5:] for t in c[12:] if t.startswith('cg:Z:')), '')
        _gaps = sum(int(n) for n, op in __import__('re').findall(r'(\d+)([IDN])', _cg))
        ident = matches / (block - _gaps) * 100 if block - _gaps > 0 else 0
        if ident < 80:
            continue
        hits[(qname, tname)].append({
            "qs": qs, "qe": qe, "qlen": qlen,
            "ts": ts, "te": te, "tlen": tlen,
            "strand": c[4], "ident": ident,
        })

count = 0
for (qname, tname), hlist in hits.items():
    if len(hlist) < 2: continue
    hlist.sort(key=lambda x: x["ts"])
    for i in range(len(hlist)):
        for j in range(i+1, len(hlist)):
            hi, hj = hlist[i], hlist[j]
            tgap = hj["ts"] - hi["te"]
            if tgap < MIN_GAP: continue
            # Different query portions?
            qov = max(0, min(hi["qe"], hj["qe"]) - max(hi["qs"], hj["qs"]))
            min_span = min(hi["qe"]-hi["qs"], hj["qe"]-hj["qs"])
            if qov > min_span * 0.5: continue
            
            total_cov = len(set(range(hi["qs"], hi["qe"])) | set(range(hj["qs"], hj["qe"])))
            cov_pct = total_cov / hi["qlen"] * 100
            
            count += 1
            if count <= 50:
                print(f"  {qname} on {tname} ({hj['tlen']:,}bp)")
                print(f"    Part1: Q[{hi['qs']}-{hi['qe']}] -> T[{hi['ts']:,}-{hi['te']:,}] {hi['strand']} {hi['ident']:.0f}%")
                print(f"    Part2: Q[{hj['qs']}-{hj['qe']}] -> T[{hj['ts']:,}-{hj['te']:,}] {hj['strand']} {hj['ident']:.0f}%")
                print(f"    Gap: {tgap:,}bp | CDS covered: {cov_pct:.0f}%")
            break
        else: continue
        break

print(f"\nTotal split CDS cases found: {count}")
PYEND
