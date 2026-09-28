#!/bin/bash
#SBATCH --job-name=xref_tests
#SBATCH --account=pc_jbei005
#SBATCH --partition=lr_bigmem
#SBATCH --qos=lr_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --output=xref_tests_%j.out

# Unit tests + synthetic end-to-end run with known IS boundaries.
#   sbatch tests/run_tests.sh [WORK_DIR]        (from the repo root)
set -euo pipefail
export PATH="$HOME/.conda/envs/claude-env/bin:$PATH"

REPO="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
WORK="${1:-${TMPDIR:-/tmp}/xref_tests_$$}"
mkdir -p "$WORK"

echo "=== unit tests ==="
python3 "$REPO/tests/test_lib_alleles.py"
python3 "$REPO/tests/test_self_align.py"

echo "=== synthetic data ==="
python3 "$REPO/tests/make_synthetic.py" "$WORK/syn"

echo "=== pipeline on synthetic data ==="
RUN="$WORK/run"
rm -rf "$RUN"; mkdir -p "$RUN"
cp "$WORK/syn/is_hits.tsv" "$RUN/"          # Stage 1 (HMM) is skipped: known hits
python3 "$REPO/scripts/pipeline.py" \
    --hmm unused.hmm --db "$WORK/syn/db.fa" --out "$RUN" \
    --distances 1000,5000 --flank 6000 --threads 4 \
    --detect-rearrangements --detect-growing-tandem --no-render-figures

echo "=== checks ==="
python3 "$REPO/tests/check_synthetic.py" "$WORK/syn" "$RUN"
