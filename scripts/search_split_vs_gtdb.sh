#!/bin/bash
#SBATCH --job-name=split_vs_gtdb
#SBATCH --partition=standard
#SBATCH --cpus-per-task=48
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/groups/rubin/projects/kuang/out/IS110_fragment/hmm_search/split_vs_gtdb_%j.out
#SBATCH --error=/groups/rubin/projects/kuang/out/IS110_fragment/hmm_search/split_vs_gtdb_%j.err

set -euo pipefail

OUTDIR="${OUTDIR:-/groups/rubin/projects/kuang/out/IS110_fragment/hmm_search}"
GTDB_FNA="${GTDB_FNA:-/groups/rubin/databases/GTDB/gtdb_genomes_v2024latest/gtdb_genomes_reps_r220/allfilestogether}"
CPUS=48
SCRIPTS="${IS110_FRAGMENT_SCRIPTS:-/global/home/users/kh36969/IS110_ecology/IS110_fragment/repo/scripts}"

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
    | minimap2 -x asm10 -c --eqx --secondary=yes -N 1000000 -p 0.05 -g 200 -r 200,200 \
        -t "$CPUS" \
        /dev/stdin "$OUTDIR/split_123_queries.fna" \
        -o "$OUTDIR/split_123_vs_gtdb.paf"

echo "  Total hits: $(wc -l < "$OUTDIR/split_123_vs_gtdb.paf")"
echo "Step 2 done: $(date)"

# Step 3: Detect split cases across all contigs
echo "Step 3: Detecting split cases..."

python3 "$SCRIPTS/detect_split_all.py" "$OUTDIR/split_123_vs_gtdb.paf" "$OUTDIR/split_123_all_hits.tsv" \
    --min-gap 50 --max-gap 20000 --min-identity 80 --min-cov 50

echo ""
echo "=== Done ==="
echo "End: $(date)"
