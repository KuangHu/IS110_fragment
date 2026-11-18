#!/bin/bash
#SBATCH --job-name=minimap2_search
#SBATCH --partition=standard
#SBATCH --cpus-per-task=48
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=minimap2_%j.out
#SBATCH --error=minimap2_%j.err

# Usage:
#   Single query:    sbatch scripts/run_minimap2.sh <query.gbk> <gtdb_dir> <output_dir> [cpus] [min_identity]
#   Multiple queries: sbatch scripts/run_minimap2.sh <gbk_dir>   <gtdb_dir> <output_dir> [cpus] [min_identity]
#
# If the first argument is a directory, all *.gbk files in it are searched in one minimap2 run.
# If the first argument is a file, only that query is searched.
#
# Arguments:
#   query         Path to a .gbk file OR a directory containing .gbk files
#   gtdb_dir      Path to directory containing *.fna genome files
#   output_dir    Output directory for results
#   cpus          CPUs for minimap2 (default: 48)
#   min_identity  Minimum identity threshold 0-1 (default: 0.80)

set -euo pipefail

# ── Parse arguments ──
if [ $# -lt 3 ]; then
    echo "Usage: sbatch scripts/run_minimap2.sh <query.gbk|gbk_dir> <gtdb_dir> <output_dir> [cpus=48] [min_identity=0.80]"
    exit 1
fi

INPUT=$(realpath "$1")
GTDB_DIR=$(realpath "$2")
OUTDIR=$(realpath "$3")
CPUS=${4:-48}
MIN_IDENT=${5:-0.80}

# ── Collect GBK files ──
if [ -d "$INPUT" ]; then
    GBK_FILES=("$INPUT"/*.gbk)
    if [ ${#GBK_FILES[@]} -eq 0 ]; then
        echo "ERROR: No .gbk files found in $INPUT" >&2
        exit 1
    fi
elif [ -f "$INPUT" ]; then
    GBK_FILES=("$INPUT")
else
    echo "ERROR: Input not found: $INPUT" >&2
    exit 1
fi

if [ ! -d "$GTDB_DIR" ]; then
    echo "ERROR: GTDB directory not found: $GTDB_DIR" >&2
    exit 1
fi

# ── Setup ──
mkdir -p "$OUTDIR"
SCRIPTS="${IS110_FRAGMENT_SCRIPTS:-/global/home/users/kh36969/IS110_ecology/IS110_fragment/repo/scripts}"
COMBINED_QUERY="${OUTDIR}/combined_query.fna"
RAW_PAF="${OUTDIR}/all_raw.paf"
> "$COMBINED_QUERY"

echo "=== minimap2 search ==="
echo "Input:        $INPUT (${#GBK_FILES[@]} GBK files)"
echo "Database:     $GTDB_DIR"
echo "Output:       $OUTDIR"
echo "CPUs:         $CPUS"
echo "Min identity: $MIN_IDENT"
echo ""

# ── Extract FASTA from each GBK ──
echo "Extracting FASTA from GBK files..."
for gbk in "${GBK_FILES[@]}"; do
    python -c "
from Bio import SeqIO
rec = SeqIO.read('${gbk}', 'genbank')
SeqIO.write(rec, open('${COMBINED_QUERY}', 'a'), 'fasta')
print(f'  {rec.id}: {len(rec)} bp')
"
done
echo ""
echo "Combined query: $(grep -c '^>' "$COMBINED_QUERY") sequences"
echo "Start: $(date)"
echo ""

# ── Run minimap2 (one pass for all queries) ──
# -c --eqx: CIGAR with =/X for precise identity calculation
# --secondary=yes -N 1000000 -p 0.05 -g 200 -r 200,200: report multiple/split alignments
# Sequence names must be unique across the files that are concatenated:
# two genomes that both have a contig "ctg" become one target, and halves of
# an IS in two different genomes were joined into one "insertion".
source "$SCRIPTS/stage_lib.sh"
require_unique_seq_names "$GTDB_DIR"
find "$GTDB_DIR" -name "*.fna" -exec cat {} + \
    | minimap2 -x asm5 -c --eqx --secondary=yes -N 1000000 -p 0.05 -g 200 -r 200,200 \
        -t "$CPUS" \
        /dev/stdin "$COMBINED_QUERY" \
        -o "$RAW_PAF"

echo ""
echo "End: $(date)"
echo "Total raw hits: $(wc -l < "$RAW_PAF")"
echo ""

# ── Split and filter per query ──
echo "Splitting results per query..."
for gbk in "${GBK_FILES[@]}"; do
    NAME=$(basename "$gbk" .gbk)
    QUERY_ID=$(python -c "from Bio import SeqIO; print(SeqIO.read('${gbk}', 'genbank').id)")
    HITS_PAF="${OUTDIR}/${NAME}_hits.paf"

    # aligned identity (gap columns excluded): `$10/$11` dropped a perfect
    # copy with a 600-bp insert (77%)
    python3 "$SCRIPTS/paf_utils.py" filter --min-identity "$MIN_IDENT" --query "$QUERY_ID" \
        < "$RAW_PAF" > "$HITS_PAF"

    N=$(wc -l < "$HITS_PAF")
    echo "  ${NAME}: ${N} hits (>= ${MIN_IDENT} identity)"
done

echo ""
echo "Done. Per-query results in: $OUTDIR/<name>_hits.paf"
