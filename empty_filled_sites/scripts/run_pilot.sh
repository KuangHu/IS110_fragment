#!/bin/bash
#SBATCH --job-name=ef_pilot
#SBATCH --account=pc_jbei005
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=/global/home/users/kh36969/IS110_ecology/IS110_fragment/empty_filled_sites/results/pilot_%j.out
#SBATCH --error=/global/home/users/kh36969/IS110_ecology/IS110_fragment/empty_filled_sites/results/pilot_%j.err

# Phase 0 pilot: empty-vs-filled site detection on a subset of IS110 references.
#
# Steps:
#   1. Take top-N most common IS110 reference variants (largest 95% clusters)
#   2. Extract anchors (outer 500 bp on each side)
#   3. Concatenate all NCBI complete E. coli genomes as search target
#   4. minimap2 anchors against the database
#   5. Pair anchors per genome and compute distance distribution
#   6. Report which refs have empty + filled (boundary-discoverable)

set -euo pipefail
export PATH="$HOME/.conda/envs/claude-env/bin:$PATH"

PROJ_DIR="/global/home/users/kh36969/IS110_ecology/IS110_fragment/empty_filled_sites"
SCRIPTS="$PROJ_DIR/scripts"
RESULTS="$PROJ_DIR/results"
REFS_FA="/global/scratch/users/kh36969/is110_refs/is110_reps_95.fa"
CLUSTER_TSV="/global/scratch/users/kh36969/is110_refs/is110_cluster_95.tsv"
GENOME_DIR="/global/scratch/users/kh36969/ncbi_refs/ecoli_complete/ecoli_complete_data/ncbi_dataset/data"
WORK="/global/scratch/users/kh36969/empty_filled_pilot"

CPUS=${SLURM_CPUS_PER_TASK:-32}
TOP_N=10
ANCHOR_LEN=500

mkdir -p "$WORK" "$RESULTS"

echo "============================================"
echo "Empty-vs-Filled Pilot (top $TOP_N refs)"
echo "============================================"
echo "Start: $(date)"
echo ""

# ────────────────────────────────────────────────
# STEP 1: Pick top-N largest clusters as pilot refs
# ────────────────────────────────────────────────
PILOT_IDS="$WORK/pilot_ref_ids.txt"
echo "Step 1: Selecting top-$TOP_N largest clusters..."
# Use temp file to avoid SIGPIPE from head closing the pipe under pipefail
TMPCNT="$WORK/cluster_counts.tsv"
awk '{ counts[$1]++ } END { for (k in counts) print counts[k]"\t"k }' "$CLUSTER_TSV" \
    | sort -k1,1 -rn > "$TMPCNT"
head -n "$TOP_N" "$TMPCNT" | cut -f2 > "$PILOT_IDS"

echo "  Selected ${TOP_N} pilot refs:"
cat "$PILOT_IDS" | head -10
echo ""

# Extract their full sequences from the refs FASTA
PILOT_REFS="$WORK/pilot_refs.fa"
python3 << PYEOF
ids = set(line.strip() for line in open("$PILOT_IDS"))
keep = False
out = open("$PILOT_REFS", "w")
with open("$REFS_FA") as f:
    for line in f:
        if line.startswith(">"):
            name = line[1:].split()[0]
            keep = name in ids
        if keep:
            out.write(line)
out.close()
PYEOF

echo "  Pilot refs FASTA: $(grep -c '^>' "$PILOT_REFS") sequences"
echo "Step 1 done: $(date)"
echo ""

# ────────────────────────────────────────────────
# STEP 2: Extract anchors (up + down) from pilot refs
# ────────────────────────────────────────────────
ANCHORS="$WORK/anchors.fa"
REF_LENGTHS="$WORK/ref_lengths.tsv"

echo "Step 2: Extracting anchors..."
python3 "$SCRIPTS/extract_anchors.py" "$PILOT_REFS" "$ANCHORS" --anchor-len "$ANCHOR_LEN"

# Build ref-length table: ref_id, total_len, anchor_len, is110_len
# IS110 length = total_len - 2*flank (flank=2000 from extract_is110_refs.sh)
FLANK_BP=2000
python3 << PYEOF
out = open("$REF_LENGTHS", "w")
name = None
seq_len = 0
with open("$PILOT_REFS") as f:
    for line in f:
        if line.startswith(">"):
            if name:
                is110_len = seq_len - 2 * $FLANK_BP
                out.write(f"{name}\t{seq_len}\t$ANCHOR_LEN\t{is110_len}\n")
            name = line[1:].split()[0]
            seq_len = 0
        else:
            seq_len += len(line.strip())
    if name:
        is110_len = seq_len - 2 * $FLANK_BP
        out.write(f"{name}\t{seq_len}\t$ANCHOR_LEN\t{is110_len}\n")
out.close()
PYEOF
echo "  Anchors: $(grep -c '^>' "$ANCHORS")"
echo "  Ref lengths: $(wc -l < "$REF_LENGTHS")"
echo "Step 2 done: $(date)"
echo ""

# ────────────────────────────────────────────────
# STEP 3: Build NCBI E. coli search database (concat all genomes)
# ────────────────────────────────────────────────
DB="$WORK/ncbi_ecoli_db.fa"

if [ -s "$DB" ]; then
    echo "Step 3: NCBI DB exists, skipping."
else
    echo "Step 3: Concatenating NCBI E. coli genomes (parallel)..."
    > "$DB"
    # Use parallel cat with file accession prefixed in headers
    find "$GENOME_DIR" -name "*_genomic.fna" -print0 \
        | xargs -0 -P "$CPUS" -I{} bash -c '
            fna="$1"
            acc=$(basename "$(dirname "$fna")")
            awk -v acc="$acc" "/^>/ { sub(/^>/, \">\" acc \"|\"); print; next } { print }" "$fna"
        ' _ {} >> "$DB"
fi

echo "  DB sequences: $(grep -c '^>' "$DB")"
echo "  DB size: $(du -h "$DB" | cut -f1)"
echo "Step 3 done: $(date)"
echo ""

# ────────────────────────────────────────────────
# STEP 4: minimap2 anchors against NCBI DB
# ────────────────────────────────────────────────
PAF="$WORK/anchors_vs_ncbi.paf"
if [ -s "$PAF" ]; then
    echo "Step 4: PAF exists, skipping."
else
    echo "Step 4: minimap2 (asm5, sensitive for short anchors)..."
    minimap2 -x asm5 -c --eqx \
        --secondary=yes -N 50 -p 0.5 \
        -t "$CPUS" \
        "$DB" "$ANCHORS" \
        -o "$PAF"
fi
echo "  PAF rows: $(wc -l < "$PAF")"
echo "Step 4 done: $(date)"
echo ""

# ────────────────────────────────────────────────
# STEP 5: Pair anchors and analyze
# ────────────────────────────────────────────────
OUT_PREFIX="$RESULTS/pilot_${SLURM_JOB_ID}"
echo "Step 5: Pairing anchors and computing distances..."
python3 "$SCRIPTS/find_anchor_pairs.py" \
    "$PAF" "$REF_LENGTHS" "$OUT_PREFIX" \
    --min-identity 95 --min-cov 90 \
    --max-pair-distance 20000 --empty-distance 100

echo ""
echo "Step 5 done: $(date)"
echo ""
echo "Results:"
echo "  $OUT_PREFIX""_pairs.tsv"
echo "  $OUT_PREFIX""_summary.tsv"
echo ""
echo "End: $(date)"
