#!/bin/bash
#SBATCH --job-name=lineage
#SBATCH --account=pc_jbei005
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/global/scratch/users/kh36969/logan_results/lineage_%j.out
#SBATCH --error=/global/scratch/users/kh36969/logan_results/lineage_%j.err

# Find IS110 fragmentation lineages (V1 -> V2 -> V3 -> ...)
# Usage: sbatch scripts/run_find_lineage.sh <organism>

set -euo pipefail
export PATH="$HOME/.conda/envs/claude-env/bin:$PATH"

ORG=${1:?Usage: sbatch run_find_lineage.sh <organism>}
SCRATCH="/global/scratch/users/kh36969"
RESULTS_DIR="${SCRATCH}/logan_results/${ORG}"
LOGAN_DIR="${SCRATCH}/logan_${ORG}"
SCRIPT="/global/home/users/kh36969/IS110_ecology/IS110_fragment/repo/scripts/find_lineage.py"

CHAINS_TSV="${RESULTS_DIR}/is110_${ORG}_multisplit_chains.tsv"
OUT_PREFIX="${RESULTS_DIR}/is110_${ORG}_lineage"

if [ ! -f "$CHAINS_TSV" ]; then
    echo "ERROR: chains TSV not found: $CHAINS_TSV"
    exit 1
fi

echo "=== Find lineage: ${ORG} ==="
echo "Chains: $CHAINS_TSV"
echo "Output: ${OUT_PREFIX}_*.tsv"
echo "Start:  $(date)"
echo ""

python3 "$SCRIPT" \
    "$CHAINS_TSV" \
    "$LOGAN_DIR" \
    "$OUT_PREFIX" \
    --min-ladder-depth 3 \
    --cluster-identity 95 \
    --max-ins-examples 5 \
    --threads "$SLURM_CPUS_PER_TASK"

echo ""
echo "End: $(date)"
