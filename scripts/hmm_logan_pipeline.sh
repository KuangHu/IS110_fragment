#!/bin/bash
# hmmsearch uses --cut_ga (the Pfam gathering thresholds in PF01548.hmm / PF02371.hmm):
# an E-value cutoff depends on the number of proteins searched, so the same
# protein could pass in one dataset and fail in another.
#SBATCH --job-name=hmm_logan
#SBATCH --account=pc_jbei005
#SBATCH --partition=lr_bigmem
#SBATCH --qos=lr_normal
#SBATCH --cpus-per-task=32
#SBATCH --mem=1000G
#SBATCH --time=3-00:00:00
#SBATCH --output=/global/scratch/users/kh36969/logan_results/hmm_logan_%j.out
#SBATCH --error=/global/scratch/users/kh36969/logan_results/hmm_logan_%j.err

# Find IS110 transposase positions in Logan contigs via:
#   1. Contig dereplication (mmseqs2 linclust at 99% id)
#   2. ORF prediction on reps (prodigal)
#   3. Protein dereplication (mmseqs2 linclust at 100% aa id)
#   4. HMM search with DEDD (PF01548) + Tnp20 (PF02371)
#   5. Map IS110-transposase proteins back to original contigs/samples
#
# Usage: sbatch scripts/hmm_logan_pipeline.sh <organism>

set -euo pipefail
export PATH="$HOME/.conda/envs/claude-env/bin:$PATH"

ORG=${1:?Usage: sbatch hmm_logan_pipeline.sh <organism>}
SCRATCH="/global/scratch/users/kh36969"
LOGAN_DIR="${SCRATCH}/logan_${ORG}"
CONCAT="${SCRATCH}/logan_results/${ORG}/all_contigs.fa"
OUTDIR="${SCRATCH}/logan_hmm/${ORG}"
HMM_DIR="/global/home/users/kh36969/IS110_ecology/IS110_fragment/repo/hmm"
CPUS=${SLURM_CPUS_PER_TASK:-32}

source "${IS110_FRAGMENT_SCRIPTS:-/global/home/users/kh36969/IS110_ecology/IS110_fragment/repo/scripts}/stage_lib.sh"
# 1 = collapse contigs at 99% identity before gene prediction (the old default).
# Members of a contig cluster are then NOT counted: two samples with the same
# contig gave one position row for one sample. Default 0: every contig is
# searched; identical proteins are still searched once (step 3) and expanded.
DEREP_CONTIGS=${DEREP_CONTIGS:-0}

mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "============================================"
echo "HMM IS110 search: ${ORG}"
echo "============================================"
echo "Logan dir:  ${LOGAN_DIR}"
echo "Concat:     ${CONCAT}"
echo "Output:     ${OUTDIR}"
echo "CPUs:       ${CPUS}"
echo "Start:      $(date)"
echo ""

# Verify concat exists (created by search_logan.sh Step 1)
if [ ! -f "$CONCAT" ]; then
    echo "ERROR: concatenated fasta not found: $CONCAT"
    echo "Run search_logan.sh first to create it."
    exit 1
fi

SIZE=$(du -b "$CONCAT" | cut -f1)
echo "Concat size: $(du -h "$CONCAT" | cut -f1) ($SIZE bytes)"
echo ""

# ────────────────────────────────────────────────
# STEP 0: Length-filter contigs (Logan has many short unitigs)
# ────────────────────────────────────────────────
MIN_LEN=1000
FILTERED="${OUTDIR}/contigs_gte${MIN_LEN}.fa"

SIG0=$(stage_sig "filter" "$MIN_LEN" "$(file_fp "$CONCAT")")
if stage_done "$FILTERED.done" "$SIG0"; then
    echo "Step 0: length-filtered contigs already exist, skipping."
else
    echo "Step 0: Filtering contigs to >= ${MIN_LEN} bp..."
    awk -v min="$MIN_LEN" '
        /^>/ {
            if (seq && length(seq) >= min) {
                print name
                print seq
            }
            name = $0
            seq = ""
            next
        }
        { seq = seq $0 }
        END {
            if (seq && length(seq) >= min) {
                print name
                print seq
            }
        }
    ' "$CONCAT" > "$FILTERED.tmp"
    mv "$FILTERED.tmp" "$FILTERED"
    stage_mark "$FILTERED.done" "$SIG0"
fi

N_FILT=$(count_fasta "$FILTERED")
if [ "$N_FILT" -eq 0 ]; then
    echo "No contig >= ${MIN_LEN} bp: nothing to search."
    printf 'protein_id\tcontig\tsample\tstart\tend\tstrand\n' > "${OUTDIR}/is110_positions.tsv"
    exit 0
fi
echo "  Contigs >= ${MIN_LEN}bp: $N_FILT ($(du -h "$FILTERED" | cut -f1))"
echo "Step 0 done: $(date)"
echo ""

# ────────────────────────────────────────────────
# STEP 1: Dereplicate contigs at 99% nt identity
# ────────────────────────────────────────────────
REPS="${OUTDIR}/contig_reps.fa"
CLUSTER_TSV="${OUTDIR}/contig_cluster.tsv"

SIG1=$(stage_sig "derep" "$DEREP_CONTIGS" "0.99 0.95 cov1" "$SIG0")
if stage_done "$REPS.done" "$SIG1"; then
    echo "Step 1: contig set already exists, skipping."
elif [ "$DEREP_CONTIGS" != "1" ]; then
    echo "Step 1: no contig dereplication (every contig is searched)."
    ln -sfn "$FILTERED" "$REPS"
    rm -f "$CLUSTER_TSV"
    stage_mark "$REPS.done" "$SIG1"
else
    echo "WARNING: DEREP_CONTIGS=1 -- samples whose contig is a cluster member are not counted."
    echo "Step 1: Dereplicating contigs at 99% identity (mmseqs2 easy-linclust)..."
    mmseqs easy-linclust \
        "$FILTERED" \
        "${OUTDIR}/contig_cluster" \
        "${OUTDIR}/mmseqs_tmp_contigs" \
        --min-seq-id 0.99 -c 0.95 --cov-mode 1 \
        --threads "$CPUS" \
        --split-memory-limit 800G \
        > "${OUTDIR}/mmseqs_contig.log" 2>&1
    mv "${OUTDIR}/contig_cluster_rep_seq.fasta" "$REPS"
    mv "${OUTDIR}/contig_cluster_cluster.tsv" "$CLUSTER_TSV"
    rm -rf "${OUTDIR}/mmseqs_tmp_contigs"
    rm -f "${OUTDIR}/contig_cluster_all_seqs.fasta"
    stage_mark "$REPS.done" "$SIG1"
fi

N_REPS=$(count_fasta "$REPS")
echo "  Contig reps: $N_REPS ($(du -h "$REPS" | cut -f1))"
echo "Step 1 done: $(date)"
echo ""

# ────────────────────────────────────────────────
# STEP 2: Predict ORFs with prodigal
# ────────────────────────────────────────────────
PROTEINS="${OUTDIR}/proteins.faa"
GFF="${OUTDIR}/proteins.gff"

SIG2=$(stage_sig "prodigal meta" "$SIG1")
if stage_done "$PROTEINS.done" "$SIG2"; then
    echo "Step 2: proteins already exist, skipping."
else
    echo "Step 2: Predicting ORFs (prodigal --meta)..."
    # Split reps into chunks for parallel prodigal. A FRESH directory: chunks
    # left by a failed attempt (other chunk count) were concatenated as well.
    SPLIT_DIR="${OUTDIR}/splits"
    rm -rf "$SPLIT_DIR"
    mkdir -p "$SPLIT_DIR"
    # Split into $CPUS chunks
    awk -v cpus="$CPUS" -v dir="$SPLIT_DIR" '
        BEGIN { total_seqs = 0 }
        /^>/ { total_seqs++ }
        END { print total_seqs > "/dev/stderr" }
    ' "$REPS" 2> "${OUTDIR}/n_reps.txt"

    N_SEQS=$(cat "${OUTDIR}/n_reps.txt")
    PER_CHUNK=$(( (N_SEQS + CPUS - 1) / CPUS ))
    echo "  Splitting $N_SEQS sequences into chunks of $PER_CHUNK..."

    awk -v per="$PER_CHUNK" -v dir="$SPLIT_DIR" '
        BEGIN { chunk = 0; count = 0; out = sprintf("%s/chunk_0.fa", dir) }
        /^>/ {
            if (count >= per) {
                close(out)
                chunk++
                count = 0
                out = sprintf("%s/chunk_%d.fa", dir, chunk)
            }
            count++
        }
        { print > out }
    ' "$REPS"

    # Run prodigal in parallel
    echo "  Running prodigal on each chunk..."
    # each chunk is checked: a bare `wait` returns 0 even if prodigal crashed
    # on a chunk, whose short .faa was then concatenated as if complete
    PIDS=()
    for f in "$SPLIT_DIR"/chunk_*.fa; do
        (prodigal -i "$f" -a "${f%.fa}.faa" -o "${f%.fa}.gff" -f gff -p meta -q > "${f%.fa}.log" 2>&1) &
        PIDS+=($!)
    done
    for pid in "${PIDS[@]}"; do
        wait "$pid" || { echo "ERROR: prodigal failed on a chunk (see $SPLIT_DIR/*.log)" >&2; exit 1; }
    done

    # Concatenate
    cat "$SPLIT_DIR"/chunk_*.faa > "$PROTEINS.tmp"
    cat "$SPLIT_DIR"/chunk_*.gff > "$GFF.tmp"
    mv "$PROTEINS.tmp" "$PROTEINS"
    mv "$GFF.tmp" "$GFF"
    rm -rf "$SPLIT_DIR"
    stage_mark "$PROTEINS.done" "$SIG2"
fi

N_PROT=$(count_fasta "$PROTEINS")
if [ "$N_PROT" -eq 0 ]; then
    echo "No protein predicted: no IS110."
    printf 'protein_id\tcontig\tsample\tstart\tend\tstrand\n' > "${OUTDIR}/is110_positions.tsv"
    exit 0
fi
echo "  Predicted proteins: $N_PROT"
echo "Step 2 done: $(date)"
echo ""

# ────────────────────────────────────────────────
# STEP 3: Dereplicate proteins at 100% identity
# ────────────────────────────────────────────────
PROT_REPS="${OUTDIR}/prot_reps.faa"
PROT_CLUSTER="${OUTDIR}/prot_cluster.tsv"

SIG3=$(stage_sig "prot derep 1.0 1.0 cov0" "$SIG2")
if stage_done "$PROT_REPS.done" "$SIG3"; then
    echo "Step 3: protein reps already exist, skipping."
else
    echo "Step 3: Dereplicating proteins at 100% identity..."
    mmseqs easy-linclust \
        "$PROTEINS" \
        "${OUTDIR}/prot_cluster" \
        "${OUTDIR}/mmseqs_tmp_prot" \
        --min-seq-id 1.0 -c 1.0 --cov-mode 0 \
        --threads "$CPUS" \
        > "${OUTDIR}/mmseqs_prot.log" 2>&1
    mv "${OUTDIR}/prot_cluster_rep_seq.fasta" "$PROT_REPS"
    mv "${OUTDIR}/prot_cluster_cluster.tsv" "$PROT_CLUSTER"
    rm -rf "${OUTDIR}/mmseqs_tmp_prot"
    rm -f "${OUTDIR}/prot_cluster_all_seqs.fasta"
    stage_mark "$PROT_REPS.done" "$SIG3"
fi

N_PROT_REPS=$(count_fasta "$PROT_REPS")
echo "  Unique proteins: $N_PROT_REPS"
echo "Step 3 done: $(date)"
echo ""

# ────────────────────────────────────────────────
# STEP 4: HMM search for DEDD + Tnp20
# ────────────────────────────────────────────────
DEDD_HITS="${OUTDIR}/dedd_hits.domtbl"
TNP20_HITS="${OUTDIR}/tnp20_hits.domtbl"

# (an interrupted hmmsearch leaves a non-empty, header-only table: "-s" skipped it)
SIG4A=$(stage_sig "hmmsearch cut_ga" "$(md5sum < "${HMM_DIR}/PF01548.hmm")" "$SIG3")
if stage_done "$DEDD_HITS.done" "$SIG4A"; then
    echo "Step 4a: DEDD hits already exist."
else
    echo "Step 4a: hmmsearch DEDD (PF01548)..."
    hmmsearch --domtblout "$DEDD_HITS.tmp" \
        --cpu "$CPUS" --cut_ga \
        "${HMM_DIR}/PF01548.hmm" "$PROT_REPS" > /dev/null
    mv "$DEDD_HITS.tmp" "$DEDD_HITS"
    stage_mark "$DEDD_HITS.done" "$SIG4A"
fi
echo "  DEDD hits: $(count_rows "$DEDD_HITS")"

SIG4B=$(stage_sig "hmmsearch cut_ga" "$(md5sum < "${HMM_DIR}/PF02371.hmm")" "$SIG3")
if stage_done "$TNP20_HITS.done" "$SIG4B"; then
    echo "Step 4b: Tnp20 hits already exist."
else
    echo "Step 4b: hmmsearch Tnp20 (PF02371)..."
    hmmsearch --domtblout "$TNP20_HITS.tmp" \
        --cpu "$CPUS" --cut_ga \
        "${HMM_DIR}/PF02371.hmm" "$PROT_REPS" > /dev/null
    mv "$TNP20_HITS.tmp" "$TNP20_HITS"
    stage_mark "$TNP20_HITS.done" "$SIG4B"
fi
echo "  Tnp20 hits: $(count_rows "$TNP20_HITS")"

# ────────────────────────────────────────────────
# STEP 5: Find proteins with BOTH domains, map back to contigs
# ────────────────────────────────────────────────
echo ""
echo "Step 5: Mapping IS110 proteins back to contigs/samples..."

python3 << PYEOF
import os, re
from collections import defaultdict

outdir = "$OUTDIR"

def parse_domtbl(path):
    hits = set()
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 1:
                hits.add(parts[0])
    return hits

dedd = parse_domtbl(os.path.join(outdir, "dedd_hits.domtbl"))
tnp20 = parse_domtbl(os.path.join(outdir, "tnp20_hits.domtbl"))
is110_prot_ids = dedd & tnp20

print(f"  DEDD proteins:  {len(dedd):,}")
print(f"  Tnp20 proteins: {len(tnp20):,}")
print(f"  Both (IS110):   {len(is110_prot_ids):,}")

with open(os.path.join(outdir, "is110_rep_protein_ids.txt"), "w") as f:
    for pid in sorted(is110_prot_ids):
        f.write(pid + "\n")

# Expand via protein cluster membership: which full proteins map to reps?
prot_members = defaultdict(set)  # rep -> {member}
with open(os.path.join(outdir, "prot_cluster.tsv")) as f:
    for line in f:
        rep, member = line.rstrip().split("\t")
        prot_members[rep].add(member)

all_is110_prots = set()
for pid in is110_prot_ids:
    all_is110_prots.update(prot_members.get(pid, {pid}))

print(f"  All IS110 protein instances (via cluster expansion): {len(all_is110_prots):,}")

# Parse GFF to get protein -> contig + position
prot_info = {}  # protein_id -> (contig, start, end, strand)
gff_path = os.path.join(outdir, "proteins.gff")
with open(gff_path) as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.rstrip().split("\t")
        if len(parts) < 9 or parts[2] != "CDS":
            continue
        contig = parts[0]
        start = int(parts[3])
        end = int(parts[4])
        strand = parts[6]
        attrs = parts[8]
        m = re.search(r"ID=([^;]+)", attrs)
        if not m:
            continue
        # Prodigal naming: <contig>_<n>
        # but the protein fasta uses <contig>_<n> as seqid
        # match ID_<contig>_<n> format
        pid_raw = m.group(1)
        # Actual protein name in faa: last field after underscore
        # Use contig + number from GFF
        # Prodigal IDs in faa look like: contig_idx (idx matches a counter)
        # We'll match via contig and start position later
        prot_name = f"{contig}_{pid_raw.split('_')[-1]}"
        prot_info[prot_name] = (contig, start, end, strand)

# Also parse protein fasta headers to get reliable naming
# Prodigal format: >contig_idx # start # end # strand # ...
faa_path = os.path.join(outdir, "proteins.faa")
prot_info2 = {}
with open(faa_path) as f:
    for line in f:
        if not line.startswith(">"):
            continue
        header = line[1:].strip()
        fields = header.split(" # ")
        if len(fields) < 4:
            continue
        pid = fields[0]
        # contig = everything before last _N
        parts = pid.rsplit("_", 1)
        if len(parts) != 2:
            continue
        contig = parts[0]
        try:
            start = int(fields[1])
            end = int(fields[2])
            strand = "+" if fields[3] == "1" else "-"
        except ValueError:
            continue
        prot_info2[pid] = (contig, start, end, strand)

print(f"  Proteins with position info: {len(prot_info2):,}")

# Write IS110 positions
# Sample = accession = contig prefix before last _
out_path = os.path.join(outdir, "is110_positions.tsv")
samples = set()
n_written = 0
with open(out_path, "w") as f:
    f.write("protein_id\tcontig\tsample\tstart\tend\tstrand\n")
    for pid in sorted(all_is110_prots):
        if pid not in prot_info2:
            continue
        contig, start, end, strand = prot_info2[pid]
        sample = contig.rsplit("_", 1)[0]
        samples.add(sample)
        f.write(f"{pid}\t{contig}\t{sample}\t{start}\t{end}\t{strand}\n")
        n_written += 1

print(f"  IS110 positions written: {n_written:,}")
print(f"  Unique samples with IS110: {len(samples):,}")
print(f"  Output: {out_path}")
PYEOF

echo ""
echo "Step 5 done: $(date)"
echo ""
echo "============================================"
echo "Pipeline complete for ${ORG}"
echo "End: $(date)"
echo "============================================"
