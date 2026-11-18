#!/bin/bash
# hmmsearch uses --cut_ga (the Pfam gathering thresholds in PF01548.hmm / PF02371.hmm):
# an E-value cutoff depends on the number of proteins searched, so the same
# protein could pass in one dataset and fail in another.
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

GENOME_DIR="${GENOME_DIR:-/global/scratch/users/kh36969/ncbi_refs/ecoli_complete/ecoli_complete_data/ncbi_dataset/data}"
OUTDIR="${OUTDIR:-/global/scratch/users/kh36969/logan_hmm/ncbi_ecoli}"
PROT_DIR="${OUTDIR}/per_genome_faa"
HMM_DIR="/global/home/users/kh36969/IS110_ecology/IS110_fragment/repo/hmm"
CPUS=${SLURM_CPUS_PER_TASK:-32}

mkdir -p "$OUTDIR" "$PROT_DIR"
source "${IS110_FRAGMENT_SCRIPTS:-/global/home/users/kh36969/IS110_ecology/IS110_fragment/repo/scripts}/stage_lib.sh"
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

# One copy of each genome, from an explicit manifest (accession <TAB> path).
# The download holds most assemblies twice, as GCA_x (GenBank) and GCF_x
# (RefSeq): 5,304 of the 13,027 directories were the twin of another one, so
# "assemblies with IS110", copies per assembly and cluster sizes all counted
# those genomes double. Rule: GCA and GCF with the same number are one
# assembly WHATEVER the version suffix; the GCF copy is kept, and of several
# versions the highest.
# The accession is read from the path (any directory or the file name), not
# assumed to be the parent directory: in a flat layout every genome got the
# same "accession" and all but one were lost.
MANIFEST="${OUTDIR}/genomes_one_per_assembly.tsv"
find "$GENOME_DIR" -name "*_genomic.fna" | sort | awk '
    { if (!match($0, /GC[AF]_[0-9]+\.[0-9]+/)) { bad++; next }
      acc = substr($0, RSTART, RLENGTH)
      # last accession in the path = the most specific one
      rest = $0; pos = 0
      while (match(rest, /GC[AF]_[0-9]+\.[0-9]+/)) { acc = substr(rest, RSTART, RLENGTH); rest = substr(rest, RSTART + RLENGTH) }
      if (acc in path) { dup++; next }
      path[acc] = $0
      split(substr(acc, 5), nv, "."); num = nv[1]; ver = nv[2] + 0
      rank = (acc ~ /^GCF_/ ? 1000000 : 0) + ver
      if (!(num in best) || rank > bestrank[num]) { best[num] = acc; bestrank[num] = rank }
    }
    END { for (num in best) print best[num] "\t" path[best[num]]
          if (bad) print bad " genome files without an accession in their path" > "/dev/stderr"
          if (dup) print dup " genome files repeat an accession" > "/dev/stderr" }' \
    | sort > "$MANIFEST.tmp"
mv "$MANIFEST.tmp" "$MANIFEST"
N_ALL=$(find "$GENOME_DIR" -name "*_genomic.fna" | wc -l)
N_GENOMES=$(wc -l < "$MANIFEST")
echo "  Genomes to process: $N_GENOMES  (of $N_ALL files; GCA/GCF twins and older versions skipped)"
if [ "$N_GENOMES" -eq 0 ]; then
    echo "ERROR: no genome with an accession found under $GENOME_DIR" >&2
    exit 1
fi

# the signature covers the genome FILES (size and modification time of each),
# not only the list of their paths: a genome replaced at the same path left
# the old proteins and everything downstream in place
SIG1=$(stage_sig "prodigal single" "$(md5sum < "$MANIFEST")" "$(cut -f2 "$MANIFEST" | files_fp)")
if stage_done "$PROTEINS.done" "$SIG1" "$PROTEINS" "$GFF"; then
    echo "Step 1: proteins exist for this manifest, skipping."
else
    echo "Step 1: Running prodigal -p single per genome (parallel)..."

    # Helper function: "accession<TAB>path" -> prodigal
    run_prodigal() {
        local acc="${1%%$'\t'*}"
        local fna="${1#*$'\t'}"
        local out="${PROT_DIR}/${acc}.faa"
        local gff="${PROT_DIR}/${acc}.gff"
        # a failure is recorded, not hidden ("|| true" left a missing or
        # partial .faa that was counted as a completed genome)
        prodigal -i "$fna" -a "$out.tmp" -o "$gff" -f gff -p single -q > /dev/null 2>&1 \
            && mv "$out.tmp" "$out" || { rm -f "$out.tmp"; echo "$acc" >> "${PROT_DIR}/FAILED.txt"; }
    }
    export -f run_prodigal
    export PROT_DIR

    rm -f "${PROT_DIR}/FAILED.txt"
    tr '\n' '\0' < "$MANIFEST" \
        | xargs -0 -P "$CPUS" -I{} bash -c 'run_prodigal "$@"' _ {}

    if [ -s "${PROT_DIR}/FAILED.txt" ]; then
        echo "ERROR: prodigal failed for $(wc -l < "${PROT_DIR}/FAILED.txt") genomes (${PROT_DIR}/FAILED.txt)" >&2
        exit 1
    fi

    # Concatenate the proteins of the MANIFEST genomes only (a *.faa glob also
    # took the files an earlier run left for genomes no longer selected, e.g.
    # the GCA twins), prefixing protein IDs with the assembly accession
    echo "  Concatenating per-genome proteins..."
    : > "$PROTEINS.tmp"
    : > "$GFF.tmp"
    while IFS=$'\t' read -r acc _path; do
        faa="${PROT_DIR}/${acc}.faa"
        if [ ! -f "$faa" ]; then
            echo "ERROR: no protein file for $acc" >&2
            exit 1
        fi
        awk -v acc="$acc" '
            /^>/ { sub(/^>/, ">" acc "|"); print; next }
            { print }
        ' "$faa" >> "$PROTEINS.tmp"
        gff="${faa%.faa}.gff"
        if [ -f "$gff" ]; then
            awk -v acc="$acc" 'BEGIN{OFS="\t"}
                /^#/ { print; next }
                { $1 = acc "|" $1; print }
            ' "$gff" >> "$GFF.tmp"
        fi
    done < "$MANIFEST"
    mv "$PROTEINS.tmp" "$PROTEINS"
    mv "$GFF.tmp" "$GFF"
    stage_mark "$PROTEINS.done" "$SIG1"
fi

N_PROT=$(grep -c '^>' "$PROTEINS" || true)
echo "  Total proteins: $N_PROT ($(du -h "$PROTEINS" | cut -f1))"
echo "Step 1 done: $(date)"
echo ""

# ────────────────────────────────────────────────
# STEP 2: HMM search for DEDD + Tnp20
# ────────────────────────────────────────────────
DEDD_HITS="${OUTDIR}/dedd_hits.domtbl"
TNP20_HITS="${OUTDIR}/tnp20_hits.domtbl"

SIG2A=$(stage_sig "hmmsearch cut_ga" "$(md5sum < "${HMM_DIR}/PF01548.hmm")" "$SIG1")
if ! stage_done "$DEDD_HITS.done" "$SIG2A"; then
    echo "Step 2a: hmmsearch DEDD (PF01548)..."
    hmmsearch --domtblout "$DEDD_HITS.tmp" \
        --cpu "$CPUS" --cut_ga \
        "${HMM_DIR}/PF01548.hmm" "$PROTEINS" > /dev/null
    mv "$DEDD_HITS.tmp" "$DEDD_HITS"
    stage_mark "$DEDD_HITS.done" "$SIG2A"
fi
echo "  DEDD hits: $(count_rows "$DEDD_HITS")"

SIG2B=$(stage_sig "hmmsearch cut_ga" "$(md5sum < "${HMM_DIR}/PF02371.hmm")" "$SIG1")
if ! stage_done "$TNP20_HITS.done" "$SIG2B"; then
    echo "Step 2b: hmmsearch Tnp20 (PF02371)..."
    hmmsearch --domtblout "$TNP20_HITS.tmp" \
        --cpu "$CPUS" --cut_ga \
        "${HMM_DIR}/PF02371.hmm" "$PROTEINS" > /dev/null
    mv "$TNP20_HITS.tmp" "$TNP20_HITS"
    stage_mark "$TNP20_HITS.done" "$SIG2B"
fi
echo "  Tnp20 hits: $(count_rows "$TNP20_HITS")"
echo ""

# ────────────────────────────────────────────────
# STEP 3: Find IS110 (both domains) and map back
# ────────────────────────────────────────────────
echo "Step 3: Mapping IS110 proteins to genomes..."
export OUTDIR
python3 << 'PYEOF'
import os
from collections import defaultdict, Counter

import os
outdir = os.environ["OUTDIR"]

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
