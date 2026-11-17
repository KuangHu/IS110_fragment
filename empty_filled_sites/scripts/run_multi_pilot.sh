#!/bin/bash
#SBATCH --job-name=ef_multi
#SBATCH --account=pc_jbei005
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=/global/home/users/kh36969/IS110_ecology/IS110_fragment/empty_filled_sites/results/multi_%j.out
#SBATCH --error=/global/home/users/kh36969/IS110_ecology/IS110_fragment/empty_filled_sites/results/multi_%j.err

# Phase 0 multi-distance pilot: re-extract IS110 references with 30 kb flanks,
# generate anchors at distances [1k, 5k, 10k, 20k, 30k] from transposase,
# search vs NCBI E. coli, find smallest distance with empty sites = IS110 boundary.
#
# Usage: sbatch scripts/run_multi_pilot.sh [pilot|all]
#        pilot = top-10 from prior pilot
#        all   = all 2,945 95%-cluster representatives

set -euo pipefail
export PATH="$HOME/.conda/envs/claude-env/bin:$PATH"

MODE=${1:-pilot}

PROJ_DIR="/global/home/users/kh36969/IS110_ecology/IS110_fragment/empty_filled_sites"
SCRIPTS="$PROJ_DIR/scripts"
RESULTS="$PROJ_DIR/results"
WORK="/global/scratch/users/kh36969/empty_filled_multi_${MODE}"

POSITIONS="/global/scratch/users/kh36969/logan_hmm/ncbi_ecoli/is110_positions.tsv"
GENOME_DIR="/global/scratch/users/kh36969/ncbi_refs/ecoli_complete/ecoli_complete_data/ncbi_dataset/data"
NCBI_DB="/global/scratch/users/kh36969/empty_filled_pilot/ncbi_ecoli_db.fa"

# For "pilot": use the same 10 IDs as before
PILOT_IDS_PRIOR="/global/scratch/users/kh36969/empty_filled_pilot/pilot_ref_ids.txt"
PILOT_REFS_PRIOR="/global/scratch/users/kh36969/empty_filled_pilot/pilot_refs.fa"
ALL_REFS="/global/scratch/users/kh36969/is110_refs/is110_reps_95.fa"

CPUS=${SLURM_CPUS_PER_TASK:-32}

mkdir -p "$WORK" "$RESULTS"

echo "============================================"
echo "Multi-distance anchor pilot: mode=$MODE"
echo "============================================"
echo "Start: $(date)"

# Choose ref subset
if [ "$MODE" = "pilot" ]; then
    REFS_FILTER="$PILOT_REFS_PRIOR"
elif [ "$MODE" = "all" ]; then
    REFS_FILTER="$ALL_REFS"
else
    echo "ERROR: unknown mode '$MODE'"
    exit 1
fi

echo "Restricting to refs in: $REFS_FILTER ($(grep -c '^>' "$REFS_FILTER") sequences)"
echo ""

# ────────────────────────────────────────────────
# STEP 1: Extract large-flank refs + multi-distance anchors
# ────────────────────────────────────────────────
ANCHORS="$WORK/anchors.fa"
REF_TABLE="$WORK/ref_table.tsv"

if [ -s "$ANCHORS" ] && [ -s "$REF_TABLE" ]; then
    echo "Step 1: extraction outputs exist, skipping."
else
    echo "Step 1: Extracting refs with 30kb flanks + multi-distance anchors..."
    python3 "$SCRIPTS/extract_multi_anchors.py" \
        --positions "$POSITIONS" \
        --genome-dir "$GENOME_DIR" \
        --refs-fa "$REFS_FILTER" \
        --out-dir "$WORK" \
        --flank 30000 --anchor-len 500 \
        --distances 1000,5000,10000,20000,30000
fi

echo "  Anchors: $(grep -c '^>' "$ANCHORS")"
echo "  Refs:    $(wc -l < "$REF_TABLE")"
echo "Step 1 done: $(date)"
echo ""

# ────────────────────────────────────────────────
# STEP 2: minimap2 anchors against NCBI db (reuse pre-built db)
# ────────────────────────────────────────────────
PAF="$WORK/anchors_vs_ncbi.paf"

if [ ! -s "$NCBI_DB" ]; then
    echo "ERROR: NCBI DB not found at $NCBI_DB. Run run_pilot.sh first to build it."
    exit 1
fi

if [ -s "$PAF" ]; then
    echo "Step 2: PAF exists, skipping minimap2."
else
    echo "Step 2: minimap2 (asm5)..."
    minimap2 -x asm5 -c --eqx \
        --secondary=yes -N 50 -p 0.5 \
        -t "$CPUS" \
        "$NCBI_DB" "$ANCHORS" \
        -o "$PAF"
fi
echo "  PAF rows: $(wc -l < "$PAF")"
echo "Step 2 done: $(date)"
echo ""

# ────────────────────────────────────────────────
# STEP 3: Pair anchors per (ref, distance) and analyze
# ────────────────────────────────────────────────
OUT_PREFIX="$RESULTS/multi_${MODE}_${SLURM_JOB_ID}"
echo "Step 3: Analyzing multi-distance anchor pairs..."
python3 "$SCRIPTS/analyze_multi_anchors.py" \
    "$PAF" "$REF_TABLE" "$OUT_PREFIX" \
    --min-identity 95 --min-cov 90

echo ""
echo "Results:"
echo "  ${OUT_PREFIX}_pairs.tsv"
echo "  ${OUT_PREFIX}_summary.tsv"
echo "End: $(date)"
