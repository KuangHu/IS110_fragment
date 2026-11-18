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
python3 "$REPO/tests/test_stage_units.py"
python3 "$REPO/tests/test_downstream.py"
python3 "$REPO/tests/test_movement.py"
python3 "$REPO/tests/test_rearr.py"

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

echo "=== same result under another hash seed (stages 2-6) ==="
# Python randomises str hashes per process; any output that depends on set
# iteration order differs between runs. Run the pipeline again with other
# seeds and compare every stage 2 and stage 6 file byte for byte.
for seed in 1 2; do
    R2="$WORK/run_seed$seed"
    rm -rf "$R2"; mkdir -p "$R2"
    cp "$WORK/syn/is_hits.tsv" "$R2/"
    PYTHONHASHSEED=$seed python3 "$REPO/scripts/pipeline.py" \
        --hmm unused.hmm --db "$WORK/syn/db.fa" --out "$R2" \
        --distances 1000,5000 --flank 6000 --threads 4 --no-render-figures > "$R2/log.txt" 2>&1
done
for f in anchors_anchors.fa anchors_refs.fa anchors_table.tsv boundaries_summary.tsv pairs_summary.tsv \
         records/records.json records/summary.tsv records/is_elements.fa records/empty_junctions.fa; do
    cmp "$WORK/run_seed1/$f" "$WORK/run_seed2/$f" || { echo "FAIL hash-seed dependence: $f"; exit 1; }
done
echo "PASS outputs identical under PYTHONHASHSEED=1 and 2 (9 files)"
