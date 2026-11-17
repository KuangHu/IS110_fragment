#!/bin/bash
#SBATCH --job-name=rebuild_logan
#SBATCH --account=pc_jbei005
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=8:00:00
#SBATCH --output=/global/scratch/users/kh36969/logan_results/rebuild_logan_%j.out
#SBATCH --error=/global/scratch/users/kh36969/logan_results/rebuild_logan_%j.err

# Rebuild the Logan E. coli concatenated FASTA without race conditions.
# Same fix as for NCBI: per-file decompression to temp, then serial concat.

set -euo pipefail
export PATH="$HOME/.conda/envs/claude-env/bin:$PATH"

LOGAN_DIR="/global/scratch/users/kh36969/logan_ecoli"
OUT="/global/scratch/users/kh36969/logan_results/ecoli/all_contigs.fa"
TMP_DIR="/global/scratch/users/kh36969/logan_results/ecoli/per_file_tmp"
CPUS=${SLURM_CPUS_PER_TASK:-32}

mkdir -p "$TMP_DIR"

echo "Rebuilding Logan E. coli concat (no race condition)..."
echo "Start: $(date)"

# Phase A: each .zst -> own .fa (parallel-safe)
ZSTD_BIN=$(which zstd)
process_one() {
    local zst="$1"
    local base=$(basename "$zst" .contigs.fa.zst)
    local out="${TMP_DIR}/${base}.fa"
    "$ZSTD_BIN" -dcq "$zst" > "$out"
}
export -f process_one
export TMP_DIR ZSTD_BIN

find "$LOGAN_DIR" -name "*.contigs.fa.zst" -print0 \
    | xargs -0 -P "$CPUS" -I{} bash -c 'process_one "$@"' _ {}

N_FILES=$(ls "$TMP_DIR"/*.fa 2>/dev/null | wc -l)
echo "  Decompressed files: $N_FILES"

# Phase B: serial cat
echo "  Concatenating serially..."
> "$OUT"
for f in "$TMP_DIR"/*.fa; do
    cat "$f" >> "$OUT"
done

# Verify integrity
N_HEADERS=$(grep -c '^>' "$OUT" || true)
N_CORRUPT=$(grep -c '[ACGT]>' "$OUT" || true)
echo "  Headers: $N_HEADERS"
echo "  Corrupt: $N_CORRUPT"
echo "  Size: $(du -h "$OUT" | cut -f1)"

# Index for samtools
echo "  Indexing..."
samtools faidx "$OUT"

# Cleanup tmp
echo "  Cleaning up tmp dir..."
rm -rf "$TMP_DIR"

echo "End: $(date)"
