#!/bin/bash
#SBATCH --job-name=detect95
#SBATCH --account=pc_jbei005
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=8:00:00
#SBATCH --output=/global/scratch/users/kh36969/logan_results/detect95_%j.out
#SBATCH --error=/global/scratch/users/kh36969/logan_results/detect95_%j.err

# Re-run variant detection at 95% identity.
set -euo pipefail
export PATH="$HOME/.conda/envs/claude-env/bin:$PATH"

ORG=${1:?Usage: sbatch detect_95.sh <organism>}
RESULTS_DIR="/global/scratch/users/kh36969/logan_results/${ORG}"
SCRIPTS="/global/home/users/kh36969/IS110_ecology/IS110_fragment/repo/scripts"
PAF="${RESULTS_DIR}/is110_vs_${ORG}.paf"

echo "=== 95% identity detection: ${ORG} ==="
echo "Start: $(date)"
python3 "${SCRIPTS}/detect_is110_variants.py" \
    "$PAF" \
    "${RESULTS_DIR}/is110_${ORG}_95" \
    --min-identity 95 --min-query-cov 50
echo "End: $(date)"
