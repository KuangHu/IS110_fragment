#!/bin/bash
#SBATCH --job-name=ef_all
#SBATCH --account=pc_jbei005
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/global/home/users/kh36969/IS110_ecology/IS110_fragment/empty_filled_sites/results/all_refs_%j.out
#SBATCH --error=/global/home/users/kh36969/IS110_ecology/IS110_fragment/empty_filled_sites/results/all_refs_%j.err

# Full pipeline: all 2,945 IS110 references with multi-distance anchors
# supporting IS110 elements up to 60 kb in size.
#
# 1. Re-extract refs with 80kb flanks (so we can put anchors out to 60kb)
# 2. Generate anchors at D = 1k, 5k, 10k, 20k, 40k, 60k, 80k from each tnp edge
# 3. minimap2 all anchors vs NCBI E. coli DB
# 4. Pair anchors per ref/D, compute IS110 sizes from empty/v0 differences
# 5. For refs with no empty in NCBI, run against Logan E. coli as fallback
# 6. Build final records.json

set -euo pipefail
export PATH="$HOME/.conda/envs/claude-env/bin:$PATH"

PROJ_DIR="/global/home/users/kh36969/IS110_ecology/IS110_fragment/empty_filled_sites"
SCRIPTS="$PROJ_DIR/scripts"
RESULTS="$PROJ_DIR/results"
WORK="/global/scratch/users/kh36969/empty_filled_all"
NCBI_DB="/global/scratch/users/kh36969/empty_filled_pilot/ncbi_ecoli_db.fa"
LOGAN_DB="/global/scratch/users/kh36969/logan_results/ecoli/all_contigs.fa"
GENOME_DIR="/global/scratch/users/kh36969/ncbi_refs/ecoli_complete/ecoli_complete_data/ncbi_dataset/data"
POSITIONS="/global/scratch/users/kh36969/logan_hmm/ncbi_ecoli/is110_positions.tsv"
ALL_REFS=""  # empty = use ALL 42,473 transposase positions (no deduplication)

CPUS=${SLURM_CPUS_PER_TASK:-32}
FLANK=80000  # 80 kb each side; total ref size = 160 kb + transposase
DISTANCES="1000,5000,10000,20000,40000,60000,80000"

mkdir -p "$WORK" "$RESULTS"

echo "============================================"
echo "Empty/Filled pipeline for ALL 2,945 IS110 refs"
echo "============================================"
echo "Flank size: $FLANK bp each side"
echo "Distances:  $DISTANCES"
echo "NCBI DB:    $NCBI_DB"
echo "Logan DB:   $LOGAN_DB"
echo "Start:      $(date)"
echo ""

# ────────────────────────────────────────────────
# STEP 1: Extract refs with large flanks + multi-D anchors
# ────────────────────────────────────────────────
ANCHORS="$WORK/anchors.fa"
REF_TABLE="$WORK/ref_table.tsv"
LARGE_REFS="$WORK/large_refs.fa"

if [ -s "$ANCHORS" ] && [ -s "$REF_TABLE" ]; then
    echo "Step 1: extraction outputs exist, skipping."
else
    echo "Step 1: Extracting refs with ${FLANK}bp flanks + multi-D anchors (all 42,473 instances)..."
    if [ -z "$ALL_REFS" ]; then
        python3 "$SCRIPTS/extract_multi_anchors.py" \
            --positions "$POSITIONS" \
            --genome-dir "$GENOME_DIR" \
            --out-dir "$WORK" \
            --flank "$FLANK" --anchor-len 500 \
            --distances "$DISTANCES"
    else
        python3 "$SCRIPTS/extract_multi_anchors.py" \
            --positions "$POSITIONS" \
            --genome-dir "$GENOME_DIR" \
            --refs-fa "$ALL_REFS" \
            --out-dir "$WORK" \
            --flank "$FLANK" --anchor-len 500 \
            --distances "$DISTANCES"
    fi
fi

N_REFS=$(wc -l < "$REF_TABLE")
N_ANCHORS=$(grep -c '^>' "$ANCHORS")
echo "  Refs:    $((N_REFS - 1))"
echo "  Anchors: $N_ANCHORS"
echo "Step 1 done: $(date)"
echo ""

# ────────────────────────────────────────────────
# STEP 2: Search anchors vs NCBI E. coli DB
# ────────────────────────────────────────────────
NCBI_PAF="$WORK/anchors_vs_ncbi.paf"

if [ -s "$NCBI_PAF" ]; then
    echo "Step 2: NCBI PAF exists, skipping."
else
    echo "Step 2: minimap2 anchors vs NCBI DB..."
    minimap2 -x asm5 -c --eqx \
        --secondary=yes -N 100 -p 0.5 \
        -t "$CPUS" \
        "$NCBI_DB" "$ANCHORS" \
        -o "$NCBI_PAF"
fi
echo "  PAF rows: $(wc -l < "$NCBI_PAF")"
echo "Step 2 done: $(date)"
echo ""

# ────────────────────────────────────────────────
# STEP 3: Analyze NCBI hits and identify boundaries
# ────────────────────────────────────────────────
NCBI_PREFIX="$RESULTS/all_ncbi_${SLURM_JOB_ID}"
echo "Step 3: Analyzing NCBI anchor pairs..."
python3 "$SCRIPTS/analyze_multi_anchors.py" \
    "$NCBI_PAF" "$REF_TABLE" "$NCBI_PREFIX" \
    --min-identity 95 --min-cov 90 \
    --max-pair-distance 200000

echo "Step 3 done: $(date)"
echo ""

# ────────────────────────────────────────────────
# STEP 4: Identify refs needing Logan fallback
# ────────────────────────────────────────────────
NEEDS_LOGAN_IDS="$WORK/needs_logan_refs.txt"

echo "Step 4: Identifying refs with no empties in NCBI..."
python3 << PYEOF
import csv
from collections import defaultdict

# Aggregate per-ref: any empty hit at any D?
has_empty = defaultdict(bool)
all_refs = set()
with open("${NCBI_PREFIX}_summary.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        rid = row["ref_id"]
        all_refs.add(rid)
        if int(row["n_empty"]) >= 1:
            has_empty[rid] = True

needs = sorted(r for r in all_refs if not has_empty[r])
with open("$NEEDS_LOGAN_IDS", "w") as f:
    for r in needs:
        f.write(r + "\n")

print(f"  Total refs: {len(all_refs)}")
print(f"  With empty in NCBI: {sum(has_empty.values())}")
print(f"  Need Logan fallback: {len(needs)}")
PYEOF
echo "Step 4 done: $(date)"
echo ""

# ────────────────────────────────────────────────
# STEP 5 (conditional): Search Logan for refs without NCBI empties
# ────────────────────────────────────────────────
N_NEEDS=$(wc -l < "$NEEDS_LOGAN_IDS")
LOGAN_ANCHORS="$WORK/anchors_for_logan.fa"
LOGAN_PAF="$WORK/anchors_vs_logan.paf"
LOGAN_PREFIX="$RESULTS/all_logan_${SLURM_JOB_ID}"

if [ "$N_NEEDS" -eq 0 ]; then
    echo "Step 5: All refs found empties in NCBI; skipping Logan."
elif [ ! -s "$LOGAN_DB" ]; then
    echo "Step 5: Logan DB not found at $LOGAN_DB; skipping fallback."
else
    echo "Step 5: Searching Logan for $N_NEEDS refs without NCBI empties..."

    if [ -s "$LOGAN_ANCHORS" ]; then
        echo "  Logan anchors subset exists, skipping extraction."
    else
        echo "  Extracting anchors for $N_NEEDS refs..."
        python3 << PYEOF
keep = set(line.strip() for line in open("$NEEDS_LOGAN_IDS"))
out = open("$LOGAN_ANCHORS", "w")
keep_now = False
with open("$ANCHORS") as f:
    for line in f:
        if line.startswith(">"):
            name = line[1:].split()[0]
            ref_id = name.rsplit("__", 1)[0]
            keep_now = ref_id in keep
        if keep_now:
            out.write(line)
out.close()
PYEOF
    fi
    echo "  Logan anchors: $(grep -c '^>' "$LOGAN_ANCHORS")"

    if [ -s "$LOGAN_PAF" ]; then
        echo "  Logan PAF exists, skipping minimap2."
    else
        echo "  minimap2 anchors vs Logan DB..."
        minimap2 -x asm5 -c --eqx \
            --secondary=yes -N 100 -p 0.5 \
            -t "$CPUS" \
            "$LOGAN_DB" "$LOGAN_ANCHORS" \
            -o "$LOGAN_PAF"
    fi
    echo "  Logan PAF rows: $(wc -l < "$LOGAN_PAF")"

    echo "  Analyzing Logan hits..."
    python3 "$SCRIPTS/analyze_multi_anchors.py" \
        "$LOGAN_PAF" "$REF_TABLE" "$LOGAN_PREFIX" \
        --min-identity 95 --min-cov 90 \
        --max-pair-distance 200000
fi
echo "Step 5 done: $(date)"
echo ""

# ────────────────────────────────────────────────
# STEP 6: Build final records JSON
# ────────────────────────────────────────────────
RECORDS_DIR="/global/scratch/users/kh36969/is110_records_all"
mkdir -p "$RECORDS_DIR"

# We use the smallest D where empty was observed for boundary calling.
# build_records.py currently requires --anchor-distance; we need to extend it
# to pick best-D per ref. For now just run at D=1000 first for the small IS110s,
# then extend in a follow-up.

echo "Step 6: Building records (best-D per ref)..."
python3 "$SCRIPTS/build_records_multi_d.py" \
    --paf "$NCBI_PAF" \
    --logan-paf "${LOGAN_PAF:-}" \
    --ref-table "$REF_TABLE" \
    --positions "$POSITIONS" \
    --genome-dir "$GENOME_DIR" \
    --logan-db "$LOGAN_DB" \
    --out-dir "$RECORDS_DIR" \
    --top-empty 5 \
    --flank-out-len 5000 \
    --min-identity 95 --min-cov 90

echo "Step 6 done: $(date)"
echo ""
echo "Records output: $RECORDS_DIR"
echo "End: $(date)"
