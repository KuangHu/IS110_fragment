#!/bin/bash
#SBATCH --job-name=cds_split_search
#SBATCH --partition=standard
#SBATCH --cpus-per-task=48
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=cds_minimap2_%j.out
#SBATCH --error=cds_minimap2_%j.err

# Search individual IS110 CDS sequences against GTDB to find split CDS cases.
# Uses asm10 preset (more sensitive for shorter/divergent sequences).

set -euo pipefail

QUERY="${QUERY:-/home/kuangh/tools/IS110_fragment/results/is110_cds_queries.fna}"
GTDB_DIR="${GTDB_DIR:-/groups/rubin/databases/GTDB/gtdb_genomes_v2024latest/gtdb_genomes_reps_r220/allfilestogether}"
OUTDIR="${OUTDIR:-/home/kuangh/tools/IS110_fragment/results/cds_search}"
CPUS="${CPUS:-48}"
SCRIPTS="${IS110_FRAGMENT_SCRIPTS:-/global/home/users/kh36969/IS110_ecology/IS110_fragment/repo/scripts}"
MIN_IDENT=0.80

mkdir -p "$OUTDIR"

echo "=== CDS-level minimap2 search ==="
echo "Query:    $QUERY ($(grep -c '^>' "$QUERY") CDS sequences)"
echo "Database: $GTDB_DIR"
echo "Output:   $OUTDIR"
echo "CPUs:     $CPUS"
echo "Preset:   asm10 (more sensitive for shorter sequences)"
echo "Start:    $(date)"
echo ""

# Run minimap2 with asm10 preset, report all secondary/supplementary alignments
# -N/-p: every copy, not the 100 best scoring >= half of the best one -- each
# half of an interrupted copy is a low-scoring secondary of any intact copy and
# was never reported (0 of 18 test configurations with an intact copy present).
# -g/-r 200: do not chain the two halves across the insert into one alignment.
# sequence names must be unique across the concatenated files (stage_lib.sh)
source "$SCRIPTS/stage_lib.sh"
require_unique_seq_names "$GTDB_DIR"
find "$GTDB_DIR" -name "*.fna" -exec cat {} + \
    | minimap2 -x asm10 -c --eqx --secondary=yes -N 1000000 -p 0.05 -g 200 -r 200,200 \
        -t "$CPUS" \
        /dev/stdin "$QUERY" \
        -o "$OUTDIR/all_cds_hits.paf"

echo ""
echo "End: $(date)"
echo "Total raw hits: $(wc -l < "$OUTDIR/all_cds_hits.paf")"
echo ""

# Filter by identity and split per CDS query
echo "Splitting results per CDS..."
while IFS= read -r header; do
    CDS_NAME=$(echo "$header" | sed 's/>//' | awk '{print $1}')
    HITS_FILE="$OUTDIR/${CDS_NAME}_hits.paf"
    
    # aligned identity (gap columns excluded), as every other reader
    python3 "$SCRIPTS/paf_utils.py" filter --min-identity "$MIN_IDENT" --query "$CDS_NAME" \
        < "$OUTDIR/all_cds_hits.paf" > "$HITS_FILE"
    
    N=$(wc -l < "$HITS_FILE")
    echo "  ${CDS_NAME}: ${N} hits"
done < <(grep '^>' "$QUERY")

echo ""
echo "Done. Per-CDS results in: $OUTDIR/<cds_name>_hits.paf"

# Quick analysis: find split CDS cases
echo ""
echo "=== Scanning for split CDS (same CDS -> distant spots on same contig) ==="
python3 "$SCRIPTS/detect_split_all.py" "$OUTDIR/all_cds_hits.paf" "$OUTDIR/split_cds.tsv" \
    --min-gap 6000 --max-gap 1000000000 --min-identity 80 --min-cov 0
