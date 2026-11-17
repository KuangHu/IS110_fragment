#!/bin/bash
#SBATCH --job-name=rearrangements
#SBATCH --account=pc_jbei005
#SBATCH --partition=lr_bigmem
#SBATCH --qos=lr_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=400G
#SBATCH --time=4:00:00
#SBATCH --output=/global/home/users/kh36969/IS110_ecology/IS110_fragment/empty_filled_sites/results/rearrangements_%j.out
#SBATCH --error=/global/home/users/kh36969/IS110_ecology/IS110_fragment/empty_filled_sites/results/rearrangements_%j.err

set -euo pipefail
export PATH="$HOME/.conda/envs/claude-env/bin:$PATH"

PROJ_DIR="/global/home/users/kh36969/IS110_ecology/IS110_fragment/empty_filled_sites"
PAF="/global/scratch/users/kh36969/empty_filled_all/anchors_vs_ncbi.paf"
OUT_PREFIX="$PROJ_DIR/results/rearrangements_${SLURM_JOB_ID}"

python3 "$PROJ_DIR/scripts/find_rearrangements.py" \
    "$PAF" "$OUT_PREFIX" \
    --min-identity 95 --min-cov 90 \
    --anchor-D 1000

echo "End: $(date)"
