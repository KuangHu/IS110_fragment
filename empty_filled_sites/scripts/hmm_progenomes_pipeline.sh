#!/bin/bash
#SBATCH --job-name=hmm_pg
#SBATCH --account=pc_jbei005
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/global/scratch/users/kh36969/progenomes_ecoli/hmm_pg_%j.out
#SBATCH --error=/global/scratch/users/kh36969/progenomes_ecoli/hmm_pg_%j.err

# Find IS110 transposases in proGenomes E. coli.
# Pipeline: prodigal -p meta on each chunk -> HMM (DEDD + Tnp20) -> map back.
#
# proGenomes data is non-redundant per-species, so NO dereplication step.
# Run prodigal in parallel per chunk (93 chunks, each ~13 GB).

set -euo pipefail
export PATH="$HOME/.conda/envs/claude-env/bin:$PATH"

CHUNK_DIR="/global/scratch/users/kh36969/progenomes_ecoli/chunks"
OUTDIR="/global/scratch/users/kh36969/progenomes_ecoli/hmm_search"
HMM_DIR="/global/home/users/kh36969/IS110_ecology/IS110_fragment/repo/hmm"
CPUS=${SLURM_CPUS_PER_TASK:-32}

mkdir -p "$OUTDIR/per_chunk_faa"

echo "============================================"
echo "HMM IS110 search on proGenomes E. coli"
echo "============================================"
echo "Chunks: $(ls $CHUNK_DIR/*.fna | wc -l)"
echo "Start: $(date)"
echo ""

# ────────────────────────────────────────────────
# STEP 1: Length-filter to >= 1000 bp (skip short contigs that can't have full ORFs)
# ────────────────────────────────────────────────
FILTERED_DIR="$OUTDIR/filtered_chunks"
mkdir -p "$FILTERED_DIR"

filter_chunk() {
    local fna="$1"
    local base=$(basename "$fna" .fna)
    local out="${FILTERED_DIR}/${base}.fa"
    if [ -s "$out" ]; then return; fi
    awk -v min=1000 '
        /^>/ {
            if (seq && length(seq) >= min) { print name; print seq }
            name = $0; seq = ""; next
        }
        { seq = seq $0 }
        END { if (seq && length(seq) >= min) { print name; print seq } }
    ' "$fna" > "$out"
}
export -f filter_chunk
export FILTERED_DIR

echo "Step 1: Length-filtering chunks to >= 1000 bp..."
find "$CHUNK_DIR" -name "*.fna" -print0 | xargs -0 -P "$CPUS" -I{} bash -c 'filter_chunk "$@"' _ {}
echo "  Filtered chunks: $(ls $FILTERED_DIR/*.fa | wc -l)"
echo "  Size: $(du -sh $FILTERED_DIR | cut -f1)"
echo "Step 1 done: $(date)"
echo ""

# ────────────────────────────────────────────────
# STEP 2: Prodigal -p meta in parallel per chunk
# ────────────────────────────────────────────────
PROT_DIR="$OUTDIR/per_chunk_faa"

run_prodigal() {
    local fa="$1"
    local base=$(basename "$fa" .fa)
    local out="${PROT_DIR}/${base}.faa"
    local gff="${PROT_DIR}/${base}.gff"
    if [ -s "$out" ] && [ -s "$gff" ]; then return; fi
    prodigal -i "$fa" -a "$out" -o "$gff" -f gff -p meta -q > /dev/null 2>&1 || true
}
export -f run_prodigal
export PROT_DIR

echo "Step 2: Running prodigal -p meta in parallel..."
find "$FILTERED_DIR" -name "*.fa" -print0 | xargs -0 -P "$CPUS" -I{} bash -c 'run_prodigal "$@"' _ {}
N_FAA=$(ls "$PROT_DIR"/*.faa 2>/dev/null | wc -l)
echo "  Per-chunk protein FASTAs: $N_FAA"
echo "Step 2 done: $(date)"
echo ""

# ────────────────────────────────────────────────
# STEP 3: Concatenate proteins (serial, no race)
# ────────────────────────────────────────────────
PROTEINS="$OUTDIR/proteins.faa"
GFF="$OUTDIR/proteins.gff"

echo "Step 3: Concatenating proteins..."
> "$PROTEINS"
> "$GFF"
for faa in "$PROT_DIR"/*.faa; do
    cat "$faa" >> "$PROTEINS"
done
for gff in "$PROT_DIR"/*.gff; do
    cat "$gff" >> "$GFF"
done

N_PROT=$(grep -c '^>' "$PROTEINS")
echo "  Total proteins: $N_PROT ($(du -h "$PROTEINS" | cut -f1))"
echo "Step 3 done: $(date)"
echo ""

# ────────────────────────────────────────────────
# STEP 4: HMM search for DEDD + Tnp20
# ────────────────────────────────────────────────
DEDD_HITS="$OUTDIR/dedd_hits.domtbl"
TNP20_HITS="$OUTDIR/tnp20_hits.domtbl"

if [ ! -s "$DEDD_HITS" ]; then
    echo "Step 4a: hmmsearch DEDD (PF01548)..."
    hmmsearch --domtblout "$DEDD_HITS" --cpu "$CPUS" \
        -E 1e-5 --domE 1e-5 \
        "${HMM_DIR}/PF01548.hmm" "$PROTEINS" > /dev/null
fi
echo "  DEDD hits: $(grep -cv '^#' "$DEDD_HITS")"

if [ ! -s "$TNP20_HITS" ]; then
    echo "Step 4b: hmmsearch Tnp20 (PF02371)..."
    hmmsearch --domtblout "$TNP20_HITS" --cpu "$CPUS" \
        -E 1e-5 --domE 1e-5 \
        "${HMM_DIR}/PF02371.hmm" "$PROTEINS" > /dev/null
fi
echo "  Tnp20 hits: $(grep -cv '^#' "$TNP20_HITS")"
echo ""

# ────────────────────────────────────────────────
# STEP 5: Map proteins with both domains back to contigs
# ────────────────────────────────────────────────
echo "Step 5: Mapping IS110 proteins to contigs..."
python3 << 'PYEOF'
import os, re
from collections import defaultdict, Counter

outdir = "/global/scratch/users/kh36969/progenomes_ecoli/hmm_search"

def parse_domtbl(path):
    hits = set()
    with open(path) as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.split()
            if parts:
                hits.add(parts[0])
    return hits

dedd = parse_domtbl(f"{outdir}/dedd_hits.domtbl")
tnp20 = parse_domtbl(f"{outdir}/tnp20_hits.domtbl")
both = dedd & tnp20

print(f"  DEDD proteins:  {len(dedd):,}")
print(f"  Tnp20 proteins: {len(tnp20):,}")
print(f"  Both (IS110):   {len(both):,}")

# Parse protein fasta header for genome + position
prot_info = {}
with open(f"{outdir}/proteins.faa") as f:
    for line in f:
        if not line.startswith(">"):
            continue
        header = line[1:].strip()
        fields = header.split(" # ")
        if len(fields) < 4:
            continue
        pid = fields[0]
        # Prodigal format: <contig>_<idx>
        contig = pid.rsplit("_", 1)[0]
        try:
            start = int(fields[1])
            end = int(fields[2])
            strand = "+" if fields[3] == "1" else "-"
        except ValueError:
            continue
        prot_info[pid] = (contig, start, end, strand)

# Write IS110 positions
out_path = f"{outdir}/is110_positions.tsv"
contigs = set()
n = 0
with open(out_path, "w") as f:
    f.write("protein_id\tcontig\tstart\tend\tstrand\tprotein_len\n")
    for pid in sorted(both):
        if pid not in prot_info:
            continue
        contig, start, end, strand = prot_info[pid]
        plen = (end - start + 1) // 3
        contigs.add(contig)
        f.write(f"{pid}\t{contig}\t{start}\t{end}\t{strand}\t{plen}\n")
        n += 1

print(f"  IS110 positions: {n:,}")
print(f"  Unique contigs with IS110: {len(contigs):,}")
print(f"  Output: {out_path}")
PYEOF

echo ""
echo "Step 5 done: $(date)"
echo ""
echo "End: $(date)"
