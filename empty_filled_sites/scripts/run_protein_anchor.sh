#!/bin/bash
#SBATCH --job-name=tnp_anchor
#SBATCH --account=pc_jbei005
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/global/home/users/kh36969/IS110_ecology/IS110_fragment/empty_filled_sites/results/tnp_anchor_%j.out
#SBATCH --error=/global/home/users/kh36969/IS110_ecology/IS110_fragment/empty_filled_sites/results/tnp_anchor_%j.err

# Define IS110 element boundaries using the transposase CDS as anchor.
# Approach: search transposase against NCBI db, extract context per hit,
# align contexts back to the source reference, find the conservation cliff.
#
# Usage: sbatch scripts/run_protein_anchor.sh [pilot|all]

set -euo pipefail
export PATH="$HOME/.conda/envs/claude-env/bin:$PATH"

MODE=${1:-pilot}

PROJ_DIR="/global/home/users/kh36969/IS110_ecology/IS110_fragment/empty_filled_sites"
SCRIPTS="$PROJ_DIR/scripts"
WORK="/global/scratch/users/kh36969/empty_filled_multi_${MODE}"
NCBI_DB="/global/scratch/users/kh36969/empty_filled_pilot/ncbi_ecoli_db.fa"
OUT="/global/scratch/users/kh36969/protein_anchor_${MODE}"

mkdir -p "$OUT"

if [ ! -s "$WORK/large_refs.fa" ] || [ ! -s "$WORK/ref_table.tsv" ]; then
    echo "ERROR: missing $WORK/large_refs.fa or ref_table.tsv (run run_multi_pilot.sh first)"
    exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
    echo "Installing samtools..."
    conda install -y -c bioconda samtools 2>&1 | tail -3
fi

echo "Start: $(date)"

python3 "$SCRIPTS/protein_anchor_boundary.py" \
    --large-refs "$WORK/large_refs.fa" \
    --ref-table "$WORK/ref_table.tsv" \
    --db "$NCBI_DB" \
    --out-dir "$OUT" \
    --context-flank 15000 \
    --threads "$SLURM_CPUS_PER_TASK" \
    --min-tnp-identity 90 \
    --min-tnp-cov 80

echo "Output: $OUT"
ls -lh "$OUT" | head -20
echo ""
echo "Summary:"
cat "$OUT/boundary_summary.tsv"
echo ""
echo "End: $(date)"
