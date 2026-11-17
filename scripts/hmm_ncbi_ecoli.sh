#!/bin/bash
#SBATCH --job-name=hmm_ncbi_ecoli
#SBATCH --account=pc_jbei005
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/global/scratch/users/kh36969/logan_results/hmm_ncbi_ecoli_%j.out
#SBATCH --error=/global/scratch/users/kh36969/logan_results/hmm_ncbi_ecoli_%j.err

# Find IS110 transposases in NCBI complete E. coli genomes (validation set).
# Pipeline: prodigal -p single per genome (parallel) -> HMM search.

set -euo pipefail
export PATH="$HOME/.conda/envs/claude-env/bin:$PATH"

GENOME_DIR="/global/scratch/users/kh36969/ncbi_refs/ecoli_complete/ecoli_complete_data/ncbi_dataset/data"
OUTDIR="/global/scratch/users/kh36969/logan_hmm/ncbi_ecoli"
PROT_DIR="${OUTDIR}/per_genome_faa"
HMM_DIR="/global/home/users/kh36969/IS110_ecology/IS110_fragment/repo/hmm"
CPUS=${SLURM_CPUS_PER_TASK:-32}

mkdir -p "$OUTDIR" "$PROT_DIR"
cd "$OUTDIR"

echo "============================================"
echo "HMM IS110 search: NCBI complete E. coli"
echo "============================================"
echo "Start: $(date)"
echo ""

# ────────────────────────────────────────────────
# STEP 1: Run prodigal -p single per genome (fast parallel)
# ────────────────────────────────────────────────
PROTEINS="${OUTDIR}/proteins.faa"
GFF="${OUTDIR}/proteins.gff"

if [ -s "$PROTEINS" ]; then
    echo "Step 1: proteins exist, skipping."
else
    echo "Step 1: Running prodigal -p single per genome (parallel)..."

    # Helper function: takes one genome path, runs prodigal
    run_prodigal() {
        local fna="$1"
        local acc=$(basename "$(dirname "$fna")")
        local out="${PROT_DIR}/${acc}.faa"
        local gff="${PROT_DIR}/${acc}.gff"
        prodigal -i "$fna" -a "$out" -o "$gff" -f gff -p single -q > /dev/null 2>&1 || true
    }
    export -f run_prodigal
    export PROT_DIR

    # Find all genome files and process in parallel
    N_GENOMES=$(find "$GENOME_DIR" -name "*_genomic.fna" | wc -l)
    echo "  Genomes to process: $N_GENOMES"

    find "$GENOME_DIR" -name "*_genomic.fna" -print0 \
        | xargs -0 -P "$CPUS" -I{} bash -c 'run_prodigal "$@"' _ {}

    N_DONE=$(ls "$PROT_DIR"/*.faa 2>/dev/null | wc -l)
    echo "  Prodigal completed: $N_DONE / $N_GENOMES genomes"

    # Concatenate all per-genome proteins, prefixing protein IDs with assembly accession
    echo "  Concatenating per-genome proteins..."
    > "$PROTEINS"
    > "$GFF"
    for faa in "$PROT_DIR"/*.faa; do
        [ -f "$faa" ] || continue
        acc=$(basename "$faa" .faa)
        # Prefix protein IDs in headers with assembly accession
        awk -v acc="$acc" '
            /^>/ { sub(/^>/, ">" acc "|"); print; next }
            { print }
        ' "$faa" >> "$PROTEINS"
        # Same for GFF (prefix seqid with accession)
        gff="${faa%.faa}.gff"
        if [ -f "$gff" ]; then
            awk -v acc="$acc" 'BEGIN{OFS="\t"}
                /^#/ { print; next }
                { $1 = acc "|" $1; print }
            ' "$gff" >> "$GFF"
        fi
    done
fi

N_PROT=$(grep -c '^>' "$PROTEINS")
echo "  Total proteins: $N_PROT ($(du -h "$PROTEINS" | cut -f1))"
echo "Step 1 done: $(date)"
echo ""

# ────────────────────────────────────────────────
# STEP 2: HMM search for DEDD + Tnp20
# ────────────────────────────────────────────────
DEDD_HITS="${OUTDIR}/dedd_hits.domtbl"
TNP20_HITS="${OUTDIR}/tnp20_hits.domtbl"

if [ ! -s "$DEDD_HITS" ]; then
    echo "Step 2a: hmmsearch DEDD (PF01548)..."
    hmmsearch --domtblout "$DEDD_HITS" \
        --cpu "$CPUS" -E 1e-5 --domE 1e-5 \
        "${HMM_DIR}/PF01548.hmm" "$PROTEINS" > /dev/null
fi
echo "  DEDD hits: $(grep -cv '^#' "$DEDD_HITS")"

if [ ! -s "$TNP20_HITS" ]; then
    echo "Step 2b: hmmsearch Tnp20 (PF02371)..."
    hmmsearch --domtblout "$TNP20_HITS" \
        --cpu "$CPUS" -E 1e-5 --domE 1e-5 \
        "${HMM_DIR}/PF02371.hmm" "$PROTEINS" > /dev/null
fi
echo "  Tnp20 hits: $(grep -cv '^#' "$TNP20_HITS")"
echo ""

# ────────────────────────────────────────────────
# STEP 3: Find IS110 (both domains) and map back
# ────────────────────────────────────────────────
echo "Step 3: Mapping IS110 proteins to genomes..."
python3 << 'PYEOF'
import os
from collections import defaultdict, Counter

outdir = "/global/scratch/users/kh36969/logan_hmm/ncbi_ecoli"

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

# Parse protein fasta header to get genome + contig + position
prot_info = {}
with open(f"{outdir}/proteins.faa") as f:
    for line in f:
        if not line.startswith(">"):
            continue
        header = line[1:].strip()
        fields = header.split(" # ")
        if len(fields) < 4: continue
        pid = fields[0]
        # Format: ASSEMBLY|contig_idx
        if "|" in pid:
            assembly, contig_part = pid.split("|", 1)
        else:
            assembly, contig_part = "unknown", pid
        try:
            start = int(fields[1])
            end = int(fields[2])
            strand = "+" if fields[3] == "1" else "-"
        except ValueError:
            continue
        # contig name is everything before the trailing _idx that prodigal adds
        contig = contig_part.rsplit("_", 1)[0]
        prot_info[pid] = (assembly, contig, start, end, strand)

# Write IS110 positions
out_path = f"{outdir}/is110_positions.tsv"
assemblies = set()
n_written = 0
with open(out_path, "w") as f:
    f.write("protein_id\tassembly\tcontig\tstart\tend\tstrand\tprotein_len\n")
    for pid in sorted(both):
        if pid not in prot_info: continue
        asm, contig, start, end, strand = prot_info[pid]
        plen = (end - start + 1) // 3
        assemblies.add(asm)
        f.write(f"{pid}\t{asm}\t{contig}\t{start}\t{end}\t{strand}\t{plen}\n")
        n_written += 1

print(f"  IS110 positions written: {n_written:,}")
print(f"  Unique assemblies with IS110: {len(assemblies):,}")
print(f"  Output: {out_path}")

# Distribution of IS110 count per assembly
counts = Counter()
for pid in both:
    if pid in prot_info:
        counts[prot_info[pid][0]] += 1

print(f"\n  IS110 copies per assembly distribution:")
copy_dist = Counter(counts.values())
for n_copies in sorted(copy_dist.keys()):
    print(f"    {n_copies} copies: {copy_dist[n_copies]:,} assemblies")
PYEOF

echo ""
echo "End: $(date)"
