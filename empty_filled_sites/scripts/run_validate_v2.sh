#!/bin/bash
#SBATCH --job-name=validate_v2
#SBATCH --account=pc_jbei005
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/global/home/users/kh36969/IS110_ecology/IS110_fragment/empty_filled_sites/results/validate_v2_%j.out
#SBATCH --error=/global/home/users/kh36969/IS110_ecology/IS110_fragment/empty_filled_sites/results/validate_v2_%j.err

# Standard-method validation of rearrangements: nucmer + show-diff + IS-mediated check.
# Usage: sbatch scripts/run_validate_v2.sh [events_tsv] [n_per_category]

set -euo pipefail
export PATH="$HOME/.conda/envs/claude-env/bin:$PATH"

PROJ_DIR="/global/home/users/kh36969/IS110_ecology/IS110_fragment/empty_filled_sites"
EVENTS=${1:-$PROJ_DIR/results/rearrangements_22204011_examples.tsv}
N_PER_CAT=${2:-15}

GENOME_DIR="/global/scratch/users/kh36969/ncbi_refs/ecoli_complete/ecoli_complete_data/ncbi_dataset/data"
POSITIONS="/global/scratch/users/kh36969/logan_hmm/ncbi_ecoli/is110_positions.tsv"
OUT_DIR="$PROJ_DIR/results/validate_v2_${SLURM_JOB_ID}"

echo "Events: $EVENTS"
echo "N per category: $N_PER_CAT"
echo "Out: $OUT_DIR"
echo "Start: $(date)"

python3 "$PROJ_DIR/scripts/rearrangement_validator.py" \
    "$EVENTS" "$GENOME_DIR" "$POSITIONS" "$OUT_DIR" \
    --n-per-cat "$N_PER_CAT" \
    --region-flank 50000

echo ""
echo "==== VERDICT TABLE ===="
column -t -s $'\t' "$OUT_DIR/validation_report.tsv" | head -60

echo ""
echo "End: $(date)"
