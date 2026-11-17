#!/bin/bash
#SBATCH --job-name=build_recs
#SBATCH --account=pc_jbei005
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=/global/home/users/kh36969/IS110_ecology/IS110_fragment/empty_filled_sites/results/build_records_%j.out
#SBATCH --error=/global/home/users/kh36969/IS110_ecology/IS110_fragment/empty_filled_sites/results/build_records_%j.err

# Build records JSON + FASTAs from the multi-anchor pilot results.
# Usage: sbatch scripts/run_build_records.sh [pilot|all]

set -euo pipefail
export PATH="$HOME/.conda/envs/claude-env/bin:$PATH"

MODE=${1:-pilot}

PROJ_DIR="/global/home/users/kh36969/IS110_ecology/IS110_fragment/empty_filled_sites"
SCRIPTS="$PROJ_DIR/scripts"
WORK="/global/scratch/users/kh36969/empty_filled_multi_${MODE}"
OUT="/global/scratch/users/kh36969/is110_records_${MODE}"

POSITIONS="/global/scratch/users/kh36969/logan_hmm/ncbi_ecoli/is110_positions.tsv"
GENOME_DIR="/global/scratch/users/kh36969/ncbi_refs/ecoli_complete/ecoli_complete_data/ncbi_dataset/data"

mkdir -p "$OUT"

echo "Build records: mode=$MODE"
echo "Start: $(date)"

python3 "$SCRIPTS/build_records.py" \
    --paf "$WORK/anchors_vs_ncbi.paf" \
    --ref-table "$WORK/ref_table.tsv" \
    --positions "$POSITIONS" \
    --genome-dir "$GENOME_DIR" \
    --out-dir "$OUT" \
    --anchor-distance 1000 \
    --top-empty 5 \
    --flank-out-len 5000 \
    --min-identity 95 --min-cov 90

echo "Output: $OUT"
ls -lh "$OUT"
echo "End: $(date)"
