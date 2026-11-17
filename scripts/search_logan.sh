#!/bin/bash
#SBATCH --job-name=is110_logan
#SBATCH --account=pc_jbei005
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --output=/global/scratch/users/kh36969/logan_results/logan_%j.out
#SBATCH --error=/global/scratch/users/kh36969/logan_results/logan_%j.err

# Search IS110 elements against Logan contigs for a single organism.
#
# Usage:
#   sbatch scripts/search_logan.sh ecoli
#   sbatch scripts/search_logan.sh efaecium
#   sbatch scripts/search_logan.sh kpneumoniae
#   sbatch scripts/search_logan.sh mavium
#   sbatch scripts/search_logan.sh belkanii
#
# Or run all 5:
#   for org in ecoli efaecium kpneumoniae mavium belkanii; do
#     sbatch scripts/search_logan.sh $org
#   done

set -euo pipefail

# ── Ensure conda tools (minimap2, zstd) are in PATH ──
export PATH="$HOME/.conda/envs/claude-env/bin:$PATH"

# ── Arguments ──
ORG=${1:?Usage: sbatch search_logan.sh <organism>}

SCRATCH="/global/scratch/users/kh36969"
LOGAN_DIR="${SCRATCH}/logan_${ORG}"
RESULTS_DIR="${SCRATCH}/logan_results/${ORG}"
QUERY="/global/scratch/users/kh36969/logan_results/is110_queries.fasta"
SCRIPTS="/global/home/users/kh36969/IS110_ecology/IS110_fragment/repo/scripts"
CPUS=${SLURM_CPUS_PER_TASK:-32}

mkdir -p "$RESULTS_DIR"

echo "============================================"
echo "IS110 vs Logan: ${ORG}"
echo "============================================"
echo "Logan dir:  ${LOGAN_DIR}"
echo "Query:      ${QUERY}"
echo "Output:     ${RESULTS_DIR}"
echo "CPUs:       ${CPUS}"
echo "Start:      $(date)"
echo ""

# ── Verify inputs ──
if [ ! -d "$LOGAN_DIR" ]; then
    echo "ERROR: Logan directory not found: $LOGAN_DIR"
    exit 1
fi
if [ ! -f "$QUERY" ]; then
    echo "ERROR: Query FASTA not found: $QUERY"
    exit 1
fi

N_FILES=$(find "$LOGAN_DIR" -name "*.contigs.fa.zst" | wc -l)
echo "Contig files: ${N_FILES}"
echo ""

# ── Step 1: Decompress and concatenate contigs ──
CONCAT="${RESULTS_DIR}/all_contigs.fa"

if [ -f "$CONCAT" ] && [ "$(stat -c%s "$CONCAT" 2>/dev/null || echo 0)" -gt 1000 ]; then
    echo "Step 1: Concatenated file exists, skipping decompression."
    echo "  Size: $(du -h "$CONCAT" | cut -f1)"
else
    echo "Step 1: Decompressing and concatenating ${N_FILES} contig files..."
    rm -f "$CONCAT"

    # Process in parallel batches, append to single file
    ZSTD_BIN=$(which zstd)
    find "$LOGAN_DIR" -name "*.contigs.fa.zst" -print0 \
        | xargs -0 -P "$CPUS" -I{} "$ZSTD_BIN" -dcq {} \
        >> "$CONCAT"

    echo "  Size: $(du -h "$CONCAT" | cut -f1)"
fi
echo "Step 1 done: $(date)"
echo ""

# ── Step 2: minimap2 search ──
PAF="${RESULTS_DIR}/is110_vs_${ORG}.paf"

if [ -f "$PAF" ] && [ "$(stat -c%s "$PAF" 2>/dev/null || echo 0)" -gt 0 ]; then
    echo "Step 2: PAF file exists, skipping minimap2."
else
    echo "Step 2: Running minimap2 (asm10, sensitive)..."
    minimap2 -x asm10 -c --eqx \
        --secondary=yes -N 100 -p 0.5 \
        -t "$CPUS" \
        "$CONCAT" "$QUERY" \
        -o "$PAF"
fi

echo "  Total hits: $(wc -l < "$PAF")"
echo "Step 2 done: $(date)"
echo ""

# ── Step 3: Detect variants (both insertion and deletion) ──
echo "Step 3: Detecting IS110 structural variants..."
python3 "${SCRIPTS}/detect_is110_variants.py" \
    "$PAF" \
    "${RESULTS_DIR}/is110_${ORG}" \
    --min-identity 95 \
    --min-query-cov 50 \
    --ins-min-target-gap 500 \
    --ins-max-target-gap 50000 \
    --del-min-query-gap 200 \
    --del-max-target-gap 200

echo ""
echo "Step 3 done: $(date)"
echo ""

# ── Step 4: Clean up large concat file (optional) ──
# Uncomment to save disk space after search:
# rm -f "$CONCAT"
# echo "Cleaned up concatenated file."

echo "============================================"
echo "Pipeline complete for ${ORG}"
echo "End: $(date)"
echo "============================================"
