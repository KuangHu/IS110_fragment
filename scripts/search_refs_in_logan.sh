#!/bin/bash
#SBATCH --job-name=refs_in_logan
#SBATCH --account=pc_jbei005
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=/global/scratch/users/kh36969/logan_results/refs_in_logan_%j.out
#SBATCH --error=/global/scratch/users/kh36969/logan_results/refs_in_logan_%j.err

# Search the 2,945 NCBI-derived IS110 reference variants (with flanks) against
# Logan E. coli contigs to validate detection coverage.

set -euo pipefail
export PATH="$HOME/.conda/envs/claude-env/bin:$PATH"

REFS="/global/scratch/users/kh36969/is110_refs/is110_reps_95.fa"
LOGAN_FA="/global/scratch/users/kh36969/logan_results/ecoli/all_contigs.fa"
OUTDIR="/global/scratch/users/kh36969/is110_refs/ecoli_validation"
CPUS=${SLURM_CPUS_PER_TASK:-32}

SCRIPTS="${IS110_FRAGMENT_SCRIPTS:-/global/home/users/kh36969/IS110_ecology/IS110_fragment/repo/scripts}"
source "$SCRIPTS/stage_lib.sh"

mkdir -p "$OUTDIR"

echo "============================================"
echo "Search NCBI IS110 refs against Logan E. coli"
echo "============================================"
echo "Refs:       $REFS ($(grep -c '^>' "$REFS") sequences)"
echo "Logan:      $LOGAN_FA ($(du -h "$LOGAN_FA" | cut -f1))"
echo "Output:     $OUTDIR"
echo "Start:      $(date)"
echo ""

PAF="${OUTDIR}/refs_vs_logan.paf"

MM2_OPTS="-x asm10 -c --eqx --secondary=yes -N 1000000 -p 0.05 -g 200 -r 200,200"
SIG=$(stage_sig "minimap2" "$MM2_OPTS" "$(md5sum < "$REFS")" "$(file_fp "$LOGAN_FA")")
if stage_done "$PAF.done" "$SIG"; then
    echo "PAF exists, skipping."
else
    echo "Running minimap2 (asm10, sensitive)..."
# -N/-p: every copy, not the 100 best scoring >= half of the best one -- each
# half of an interrupted copy is a low-scoring secondary of any intact copy and
# was never reported (0 of 18 test configurations with an intact copy present).
# -g/-r 200: do not chain the two halves across the insert into one alignment.
    # shellcheck disable=SC2086
    minimap2 $MM2_OPTS \
        -t "$CPUS" \
        "$LOGAN_FA" "$REFS" \
        -o "$PAF.tmp"
    mv "$PAF.tmp" "$PAF"
    stage_mark "$PAF.done" "$SIG"
fi

echo "  Total hits: $(wc -l < "$PAF")"
echo "Search done: $(date)"
echo ""

# Was the IS (the core of each reference) found? A hit in the flanks alone is
# not a detection (see ref_core_coverage.py).
echo "Coverage analysis..."
python3 "$SCRIPTS/ref_core_coverage.py" "$PAF" "$REFS" "${OUTDIR}/ref_core_coverage.tsv" \
    --min-identity 95 --min-core-cov 90

echo ""
echo "End: $(date)"
