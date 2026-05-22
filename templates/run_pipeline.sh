#!/bin/bash
#SBATCH --job-name=is_boundary
#SBATCH --account=pc_jbei005
#SBATCH --partition=lr_bigmem
#SBATCH --qos=lr_normal
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# Cross-Reference IS Boundary Pipeline runner template.
#
# Edit the variables below or pass them on the command line, e.g.:
#   sbatch run_pipeline.sh \
#       /path/to/IS_transposase.hmm \
#       /path/to/genome_db.fa \
#       /path/to/output_dir/

set -euo pipefail
export PATH="$HOME/.conda/envs/claude-env/bin:$PATH"

# === Required arguments ===
HMM_PROFILES="${1:?Usage: $0 HMM_PROFILES GENOME_DB OUTPUT_DIR}"
GENOME_DB="${2:?Need GENOME_DB}"
OUT_DIR="${3:?Need OUTPUT_DIR}"

# === Optional tuning ===
DISTANCES="${DISTANCES:-1000,5000,20000,40000,80000}"
ANCHOR_LEN="${ANCHOR_LEN:-500}"
MIN_ID="${MIN_ID:-95}"
MIN_COV="${MIN_COV:-80}"
HIST_BIN="${HIST_BIN:-50}"
MIN_IS_SIZE="${MIN_IS_SIZE:-800}"
MAX_IS_SIZE="${MAX_IS_SIZE:-200000}"
THREADS="${SLURM_CPUS_PER_TASK:-32}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/scripts"

python3 "$SCRIPT_DIR/pipeline.py" \
    --hmm "$HMM_PROFILES" \
    --db "$GENOME_DB" \
    --out "$OUT_DIR" \
    --distances "$DISTANCES" \
    --anchor-length "$ANCHOR_LEN" \
    --min-identity "$MIN_ID" \
    --min-coverage "$MIN_COV" \
    --histogram-bin "$HIST_BIN" \
    --min-is-size "$MIN_IS_SIZE" \
    --max-is-size "$MAX_IS_SIZE" \
    --threads "$THREADS"
