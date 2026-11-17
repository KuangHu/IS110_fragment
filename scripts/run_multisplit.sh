#!/bin/bash
#SBATCH --job-name=multisplit
#SBATCH --account=pc_jbei005
#SBATCH --partition=lr_bigmem
#SBATCH --qos=lr_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=1000G
#SBATCH --time=12:00:00
#SBATCH --output=/global/scratch/users/kh36969/logan_results/multisplit_%j.out
#SBATCH --error=/global/scratch/users/kh36969/logan_results/multisplit_%j.err

# Run multi-split detection on one organism's PAF.
# Usage: sbatch scripts/run_multisplit.sh <organism>

set -euo pipefail

export PATH="$HOME/.conda/envs/claude-env/bin:$PATH"

ORG=${1:?Usage: sbatch run_multisplit.sh <organism>}
RESULTS_DIR="/global/scratch/users/kh36969/logan_results/${ORG}"
SCRIPTS="/global/home/users/kh36969/IS110_ecology/IS110_fragment/repo/scripts"

PAF="${RESULTS_DIR}/is110_vs_${ORG}.paf"

if [ ! -f "$PAF" ]; then
    echo "ERROR: PAF not found: $PAF"
    exit 1
fi

echo "=== Multi-split detection: ${ORG} ==="
echo "PAF: $PAF ($(du -h "$PAF" | cut -f1))"
echo "Start: $(date)"
echo ""

python3 "${SCRIPTS}/detect_multisplit.py" \
    "$PAF" \
    "${RESULTS_DIR}/is110_${ORG}_multisplit" \
    --min-identity 95 \
    --min-query-cov 50

echo ""
echo "End: $(date)"
