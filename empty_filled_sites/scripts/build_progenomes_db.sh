#!/bin/bash
#SBATCH --job-name=pg_db
#SBATCH --account=pc_jbei005
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --output=/global/scratch/users/kh36969/progenomes_ecoli/build_db_%j.out
#SBATCH --error=/global/scratch/users/kh36969/progenomes_ecoli/build_db_%j.err

# Concatenate proGenomes E. coli chunks into a single search DB and index.
# SERIAL concat (no race condition).

set -euo pipefail
export PATH="$HOME/.conda/envs/claude-env/bin:$PATH"

CHUNK_DIR="/global/scratch/users/kh36969/progenomes_ecoli/chunks"
DB="/global/scratch/users/kh36969/progenomes_ecoli/progenome_ecoli_db.fa"

echo "=== Build proGenomes search DB ==="
echo "Chunks: $(ls $CHUNK_DIR/*.fna | wc -l)"
echo "Start: $(date)"

# Phase A: serial concat (no race)
> "$DB"
for f in "$CHUNK_DIR"/*.fna; do
    cat "$f" >> "$DB"
done

# Verify
N_HEADERS=$(grep -c '^>' "$DB" || true)
N_CORRUPT=$(grep -c '[ACGT]>' "$DB" || true)
echo "  Headers: $N_HEADERS"
echo "  Corrupt records: $N_CORRUPT"
echo "  DB size: $(du -h "$DB" | cut -f1)"

# Index
echo "  Indexing..."
samtools faidx "$DB"

echo "End: $(date)"
