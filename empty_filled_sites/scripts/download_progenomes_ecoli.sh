#!/bin/bash
# Download proGenomes E. coli genomes.
# Run on lrc-xfer (no SU cost) inside tmux.
#
# Usage:
#   ssh lrc-xfer.lbl.gov
#   tmux new -s progenomes_dl
#   bash /global/home/users/kh36969/IS110_ecology/IS110_fragment/empty_filled_sites/scripts/download_progenomes_ecoli.sh

set -euo pipefail

OUTDIR="/global/scratch/users/kh36969/progenomes_ecoli"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "=== proGenomes v4 download ==="
echo "Target: $OUTDIR"
echo "Start: $(date)"
echo ""

# Step 1: Download proGenomes v4 representatives (all species, ~1.9M genomes representatives)
# We'll then filter to E. coli using taxonomy
echo "Step 1: Downloading proGenomes v4 representative genomes (all species)..."
echo "  WARNING: Full file is large (~50+ GB compressed)."
echo ""

# Try common known URLs (proGenomes 4 / 3 / 2)
URLS=(
    "https://progenomes.embl.de/data/repGenomes/pg4_genomes_representatives.fna.gz"
    "https://progenomes.embl.de/data/proGenomes3.contigs-representatives.fasta.bz2"
    "http://progenomes2.embl.de/data/proGenomes2.1.specI_clustering.tab"
)

for url in "${URLS[@]}"; do
    echo "  Trying: $url"
    code=$(curl -s -o /dev/null -w "%{http_code}" --max-time 30 "$url" || echo "000")
    echo "    HTTP $code"
done

echo ""
echo "Step 2: Try downloading representative genomes (this takes a while)..."
# This will fetch the all-species representative set, then filter
wget -c "https://progenomes.embl.de/data/repGenomes/pg4_genomes_representatives.fna.gz" -O pg4_genomes_representatives.fna.gz || \
    echo "  Direct URL failed; please go to https://progenomes.embl.de/download.cgi to find the right URL"

# Step 3: Get taxonomy mapping so we can filter
wget -c "https://progenomes.embl.de/data/repGenomes/pg4_specI_lineageNCBI.tab" -O pg4_specI_lineageNCBI.tab 2>/dev/null || true
wget -c "https://progenomes.embl.de/data/repGenomes/pg4_specI_clustering.tab" -O pg4_specI_clustering.tab 2>/dev/null || true

echo ""
echo "Step 4: Filter representatives for E. coli (taxid 562)"
if [ -f "pg4_genomes_representatives.fna.gz" ] && [ -f "pg4_specI_lineageNCBI.tab" ]; then
    # Find specI clusters that include taxid 562 (E. coli)
    grep -E "^[^\t]+\t.*\t562\t|^[^\t]+\t.*Escherichia coli" pg4_specI_lineageNCBI.tab > ecoli_specI.tab || true
    echo "  E. coli specI clusters: $(wc -l < ecoli_specI.tab)"
fi

echo ""
echo "Files downloaded:"
ls -lh

echo ""
echo "End: $(date)"
