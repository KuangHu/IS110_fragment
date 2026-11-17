#!/bin/bash
#SBATCH --job-name=find_src
#SBATCH --account=pc_jbei005
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/global/scratch/users/kh36969/logan_results/find_src_%j.out
#SBATCH --error=/global/scratch/users/kh36969/logan_results/find_src_%j.err

# Find autologous source of IS110 insertions in the same genome.
# Usage: sbatch scripts/run_find_source.sh <organism>

set -euo pipefail
export PATH="$HOME/.conda/envs/claude-env/bin:$PATH"

ORG=${1:?Usage: sbatch run_find_source.sh <organism>}
SCRATCH="/global/scratch/users/kh36969"
RESULTS_DIR="${SCRATCH}/logan_results/${ORG}"
LOGAN_DIR="${SCRATCH}/logan_${ORG}"
SCRIPT="/global/home/users/kh36969/IS110_ecology/IS110_fragment/repo/scripts/find_insertion_source.py"

INS_TSV="${RESULTS_DIR}/is110_${ORG}_95_insertions.tsv"
OUT_TSV="${RESULTS_DIR}/is110_${ORG}_source.tsv"

if [ ! -f "$INS_TSV" ]; then
    echo "ERROR: insertion TSV not found: $INS_TSV"
    exit 1
fi

echo "=== Find autologous source: ${ORG} ==="
echo "Insertions: $INS_TSV ($(wc -l < "$INS_TSV") rows)"
echo "Logan dir:  $LOGAN_DIR"
echo "Output:     $OUT_TSV"
echo "Start:      $(date)"
echo ""

python3 "$SCRIPT" \
    "$INS_TSV" \
    "$LOGAN_DIR" \
    "$OUT_TSV" \
    --threads "$SLURM_CPUS_PER_TASK" \
    --min-identity 95 \
    --min-coverage 80 \
    --min-ins-len 200

echo ""
echo "End: $(date)"
