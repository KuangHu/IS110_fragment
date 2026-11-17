#!/bin/bash
#SBATCH --job-name=rebuild_db
#SBATCH --account=pc_jbei005
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=/global/scratch/users/kh36969/empty_filled_pilot/rebuild_%j.out
#SBATCH --error=/global/scratch/users/kh36969/empty_filled_pilot/rebuild_%j.err

# Rebuild the NCBI E. coli concat database WITHOUT race conditions.
# Previous build used `xargs -P >> file` which interleaved output.
# This build uses per-file output then a serial cat.

set -euo pipefail
export PATH="$HOME/.conda/envs/claude-env/bin:$PATH"

GENOME_DIR="/global/scratch/users/kh36969/ncbi_refs/ecoli_complete/ecoli_complete_data/ncbi_dataset/data"
WORK="/global/scratch/users/kh36969/empty_filled_pilot"
DB="$WORK/ncbi_ecoli_db.fa"
TMP_DIR="$WORK/db_per_file_tmp"
CPUS=${SLURM_CPUS_PER_TASK:-32}

mkdir -p "$TMP_DIR"

echo "Rebuilding NCBI db: $DB"
echo "Start: $(date)"

# Phase A: each genome file processed to its own temp output (parallel-safe)
process_one() {
    local fna="$1"
    local acc=$(basename "$(dirname "$fna")")
    local out="${TMP_DIR}/${acc}.fa"
    awk -v acc="$acc" '/^>/ { sub(/^>/, ">" acc "|"); print; next } { print }' "$fna" > "$out"
}
export -f process_one
export TMP_DIR

find "$GENOME_DIR" -name "*_genomic.fna" -print0 \
    | xargs -0 -P "$CPUS" -I{} bash -c 'process_one "$@"' _ {}

N_FILES=$(ls "$TMP_DIR"/*.fa 2>/dev/null | wc -l)
echo "Per-file outputs: $N_FILES"

# Phase B: serial cat (no race condition)
echo "Concatenating serially..."
> "$DB"
for f in "$TMP_DIR"/*.fa; do
    cat "$f" >> "$DB"
done

# Verify integrity
N_HEADERS=$(grep -c '^>' "$DB")
N_CORRUPT=$(grep -c '[ACGT]>' "$DB" || true)
echo "  Headers: $N_HEADERS"
echo "  Corrupt (seq>header on same line): $N_CORRUPT"
echo "  DB size: $(du -h "$DB" | cut -f1)"

# Index for samtools faidx use
echo "Indexing..."
samtools faidx "$DB"

# Cleanup tmp
echo "Cleaning up tmp dir..."
rm -rf "$TMP_DIR"

# Invalidate downstream PAFs that used the old DB (they need regeneration)
rm -f "$WORK/anchors_vs_ncbi.paf"
rm -f /global/scratch/users/kh36969/empty_filled_multi_pilot/anchors_vs_ncbi.paf
echo "Cleared downstream PAFs (will be regenerated when pipelines re-run)"

echo "End: $(date)"
