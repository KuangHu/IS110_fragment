#!/bin/bash
# Download progenome E. coli from Hugging Face dataset hukuang/progenome_ecoli.
# Run on lrc-xfer (no SU cost) inside tmux.
#
# Usage:
#   ssh lrc-xfer.lbl.gov
#   tmux new -s progenomes_dl
#   bash /global/home/users/kh36969/IS110_ecology/IS110_fragment/empty_filled_sites/scripts/download_progenomes_hf.sh

set -euo pipefail
export PATH="$HOME/.conda/envs/claude-env/bin:$PATH"

OUTDIR="/global/scratch/users/kh36969/progenomes_ecoli"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "=== HuggingFace download: hukuang/progenome_ecoli ==="
echo "Target:  $OUTDIR"
echo "Start:   $(date)"
echo ""

# Step 1: install huggingface_hub if missing
if ! command -v hf >/dev/null 2>&1; then
    echo "Installing huggingface_hub..."
    pip install -q -U huggingface_hub
fi

# Step 2: download (resumable, parallel)
hf download hukuang/progenome_ecoli \
    --repo-type dataset \
    --local-dir "$OUTDIR" \
    --max-workers 8

echo ""
echo "=== Download summary ==="
du -sh "$OUTDIR"
ls -lh "$OUTDIR"
ls "$OUTDIR/chunks" 2>/dev/null | head -10
echo ""
echo "End: $(date)"
