#!/bin/bash
#SBATCH --job-name=frag_local
#SBATCH --account=pc_jbei005
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/global/scratch/users/kh36969/logan_results/fraglocal_%j.out
#SBATCH --error=/global/scratch/users/kh36969/logan_results/fraglocal_%j.err

# Search fragment cluster reps against all Logan contigs (all organisms)
# to find where each fragment originates.
#
# Usage: sbatch scripts/find_fragment_origin_local.sh <organism>

set -euo pipefail

export PATH="$HOME/.conda/envs/claude-env/bin:$PATH"

ORG=${1:?Usage: sbatch find_fragment_origin_local.sh <organism>}
SCRATCH="/global/scratch/users/kh36969"
RESULTS_DIR="${SCRATCH}/logan_results/${ORG}"

QUERY="${RESULTS_DIR}/fragments_${ORG}_cluster_top50.fasta"
CPUS=${SLURM_CPUS_PER_TASK:-32}
SCRIPTS="${IS110_FRAGMENT_SCRIPTS:-/global/home/users/kh36969/IS110_ecology/IS110_fragment/repo/scripts}"

if [ ! -f "$QUERY" ]; then
    echo "ERROR: $QUERY not found"
    exit 1
fi

echo "=== Local fragment origin search: ${ORG} ==="
echo "Query: $QUERY ($(grep -c '^>' "$QUERY") cluster reps)"
echo "Start: $(date)"
echo ""

# Search against all 5 organism contig databases
for target_org in ecoli efaecium kpneumoniae mavium belkanii; do
    TARGET_DB="${SCRATCH}/logan_results/${target_org}/all_contigs.fa"
    OUT_PAF="${RESULTS_DIR}/fragments_${ORG}_vs_${target_org}.paf"

    if [ ! -f "$TARGET_DB" ]; then
        echo "  ${target_org}: database not found, skipping"
        rm -f "$OUT_PAF"        # no PAF of an earlier run in the summary
        continue
    fi

    echo "--- Searching against ${target_org} ($(du -h "$TARGET_DB" | cut -f1)) ---"
    minimap2 -x asm5 -c --eqx \
        --secondary=yes -N 5 -p 0.8 \
        -t "$CPUS" \
        "$TARGET_DB" "$QUERY" \
        -o "$OUT_PAF.tmp"
    mv "$OUT_PAF.tmp" "$OUT_PAF"
    echo "  Hits: $(wc -l < "$OUT_PAF")"
    echo "  $(date)"
done

echo ""
echo "=== Summary per cluster representative ==="
# presence in other datasets, not origin; see fragment_origin_summary.py
PAF_ARGS=()
for target_org in ecoli efaecium kpneumoniae mavium belkanii; do
    PAF_ARGS+=("${target_org}=${RESULTS_DIR}/fragments_${ORG}_vs_${target_org}.paf")
done
python3 "$SCRIPTS/fragment_origin_summary.py" "$QUERY" \
    "${RESULTS_DIR}/fragments_${ORG}_origin_summary.tsv" "${PAF_ARGS[@]}" \
    --min-identity 90 --min-cov 80

echo ""
echo "End: $(date)"
