#!/bin/bash
# Download all complete E. coli genomes from NCBI Datasets.
# Run on lrc-xfer (no SU cost) inside tmux.
#
# Usage:
#   ssh lrc-xfer.lbl.gov
#   tmux new -s ecoli_dl
#   bash scripts/download_ncbi_ecoli.sh

set -euo pipefail
export PATH="$HOME/.conda/envs/claude-env/bin:$PATH"

OUTDIR="/global/scratch/users/kh36969/ncbi_refs/ecoli_complete"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "Downloading all complete E. coli genomes from NCBI..."
echo "Target: $OUTDIR"
echo "Start: $(date)"

# Use NCBI datasets to download all complete E. coli genomes
# --include genome: just the FASTA files (no proteins/GFF/etc., we'll predict ourselves)
# --dehydrated: get list first, then rehydrate (more reliable for large downloads)

datasets download genome taxon "Escherichia coli" \
    --assembly-level complete \
    --include genome \
    --filename ecoli_complete.zip \
    --dehydrated

unzip -o ecoli_complete.zip -d ecoli_complete_data
datasets rehydrate --directory ecoli_complete_data --max-workers 8

echo ""
echo "Download complete: $(date)"
echo "Genome files:"
find ecoli_complete_data/ncbi_dataset/data -name "*.fna" | wc -l
echo "Total size:"
du -sh ecoli_complete_data
