#!/bin/bash
#SBATCH --job-name=validate_rearr
#SBATCH --account=pc_jbei005
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=/global/home/users/kh36969/IS110_ecology/IS110_fragment/empty_filled_sites/results/validate_%j.out
#SBATCH --error=/global/home/users/kh36969/IS110_ecology/IS110_fragment/empty_filled_sites/results/validate_%j.err

set -euo pipefail
export PATH="$HOME/.conda/envs/claude-env/bin:$PATH"

PROJ_DIR="/global/home/users/kh36969/IS110_ecology/IS110_fragment/empty_filled_sites"
EXAMPLES="$PROJ_DIR/results/rearrangements_22204011_examples.tsv"
GENOME_DIR="/global/scratch/users/kh36969/ncbi_refs/ecoli_complete/ecoli_complete_data/ncbi_dataset/data"
OUT_DIR="$PROJ_DIR/results/validate_${SLURM_JOB_ID}"

python3 "$PROJ_DIR/scripts/validate_rearrangements.py" \
    "$EXAMPLES" "$GENOME_DIR" "$OUT_DIR" \
    --n-per-category 5 \
    --region-flank 50000

echo ""
echo "Report:"
cat "$OUT_DIR/validation_report.tsv"
