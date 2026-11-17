#!/bin/bash
#SBATCH --job-name=build_ncbi
#SBATCH --account=pc_jbei005
#SBATCH --partition=lr_bigmem
#SBATCH --qos=lr_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=900G
#SBATCH --time=12:00:00
#SBATCH --output=/global/home/users/kh36969/IS110_ecology/IS110_fragment/empty_filled_sites/results/build_ncbi_%j.out
#SBATCH --error=/global/home/users/kh36969/IS110_ecology/IS110_fragment/empty_filled_sites/results/build_ncbi_%j.err

# Build records.json from the all-refs NCBI analysis (skipping Logan).
set -euo pipefail
export PATH="$HOME/.conda/envs/claude-env/bin:$PATH"

PROJ_DIR="/global/home/users/kh36969/IS110_ecology/IS110_fragment/empty_filled_sites"
SCRIPTS="$PROJ_DIR/scripts"
WORK="/global/scratch/users/kh36969/empty_filled_all"
GENOME_DIR="/global/scratch/users/kh36969/ncbi_refs/ecoli_complete/ecoli_complete_data/ncbi_dataset/data"
POSITIONS="/global/scratch/users/kh36969/logan_hmm/ncbi_ecoli/is110_positions.tsv"
RECORDS_DIR="/global/scratch/users/kh36969/is110_records_all"

mkdir -p "$RECORDS_DIR"

echo "Building records (NCBI only)..."
echo "Start: $(date)"

python3 "$SCRIPTS/build_records_multi_d.py" \
    --paf "$WORK/anchors_vs_ncbi.paf" \
    --ref-table "$WORK/ref_table.tsv" \
    --positions "$POSITIONS" \
    --genome-dir "$GENOME_DIR" \
    --out-dir "$RECORDS_DIR" \
    --top-empty 5 \
    --flank-out-len 5000 \
    --min-identity 95 --min-cov 90

echo "Output: $RECORDS_DIR"
ls -lh "$RECORDS_DIR"
echo "End: $(date)"
