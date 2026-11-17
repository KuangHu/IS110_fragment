#!/bin/bash
#SBATCH --job-name=ef_pilot_logan
#SBATCH --account=pc_jbei005
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=/global/home/users/kh36969/IS110_ecology/IS110_fragment/empty_filled_sites/results/pilot_logan_%j.out
#SBATCH --error=/global/home/users/kh36969/IS110_ecology/IS110_fragment/empty_filled_sites/results/pilot_logan_%j.err

# Phase 0 pilot — Logan E. coli edition.
# Uses the same 10 pilot IS110 references and anchors as run_pilot.sh,
# but searches against Logan E. coli contigs (2.9 TB) instead of NCBI.

set -euo pipefail
export PATH="$HOME/.conda/envs/claude-env/bin:$PATH"

PROJ_DIR="/global/home/users/kh36969/IS110_ecology/IS110_fragment/empty_filled_sites"
SCRIPTS="$PROJ_DIR/scripts"
RESULTS="$PROJ_DIR/results"
WORK="/global/scratch/users/kh36969/empty_filled_pilot"  # share with NCBI pilot
LOGAN_FA="/global/scratch/users/kh36969/logan_results/ecoli/all_contigs.fa"

CPUS=${SLURM_CPUS_PER_TASK:-32}

# Reuse pilot refs / anchors / ref_lengths from NCBI pilot
ANCHORS="$WORK/anchors.fa"
REF_LENGTHS="$WORK/ref_lengths.tsv"

if [ ! -s "$ANCHORS" ] || [ ! -s "$REF_LENGTHS" ]; then
    echo "ERROR: anchors or ref_lengths missing. Run run_pilot.sh first."
    exit 1
fi

mkdir -p "$RESULTS"

echo "============================================"
echo "Empty-vs-Filled Pilot vs LOGAN E. coli"
echo "============================================"
echo "Anchors:    $ANCHORS ($(grep -c '^>' "$ANCHORS"))"
echo "Logan DB:   $LOGAN_FA ($(du -h "$LOGAN_FA" | cut -f1))"
echo "Start:      $(date)"
echo ""

# minimap2 against Logan
PAF="$WORK/anchors_vs_logan.paf"
if [ -s "$PAF" ]; then
    echo "PAF exists, skipping minimap2."
else
    echo "Running minimap2 (asm5, sensitive for short anchors)..."
    minimap2 -x asm5 -c --eqx \
        --secondary=yes -N 50 -p 0.5 \
        -t "$CPUS" \
        "$LOGAN_FA" "$ANCHORS" \
        -o "$PAF"
fi
echo "  PAF rows: $(wc -l < "$PAF")"
echo "minimap2 done: $(date)"
echo ""

# Pair anchors and analyze
OUT_PREFIX="$RESULTS/pilot_logan_${SLURM_JOB_ID}"
echo "Pairing anchors and computing distances..."
python3 "$SCRIPTS/find_anchor_pairs.py" \
    "$PAF" "$REF_LENGTHS" "$OUT_PREFIX" \
    --min-identity 95 --min-cov 90 \
    --max-pair-distance 20000

echo ""
echo "Results:"
echo "  ${OUT_PREFIX}_pairs.tsv"
echo "  ${OUT_PREFIX}_summary.tsv"
echo ""
echo "End: $(date)"
