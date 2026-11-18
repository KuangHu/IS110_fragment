#!/bin/bash
#SBATCH --job-name=frag_origin
#SBATCH --account=pc_jbei005
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/global/scratch/users/kh36969/logan_results/fragorig_%j.out
#SBATCH --error=/global/scratch/users/kh36969/logan_results/fragorig_%j.err

# Extract insertion fragments, cluster them, and write representatives.
# Usage: sbatch scripts/find_fragment_origin.sh <organism>

set -euo pipefail

export PATH="$HOME/.conda/envs/claude-env/bin:$PATH"
module load bio/blast-plus/2.14.1-gcc-11.4.0 2>/dev/null || true

ORG=${1:?Usage: sbatch find_fragment_origin.sh <organism>}
SCRATCH="${SCRATCH:-/global/scratch/users/kh36969}"
RESULTS_DIR="${RESULTS_DIR:-${SCRATCH}/logan_results/${ORG}}"
LOGAN_DIR="${LOGAN_DIR:-${SCRATCH}/logan_${ORG}}"
SCRIPTS="${IS110_FRAGMENT_SCRIPTS:-/global/home/users/kh36969/IS110_ecology/IS110_fragment/repo/scripts}"
source "$SCRIPTS/stage_lib.sh"

INS_TSV="${RESULTS_DIR}/is110_${ORG}_95_insertions.tsv"
FRAG_FASTA="${RESULTS_DIR}/fragments_${ORG}.fasta"
CLUSTER_PREFIX="${RESULTS_DIR}/fragments_${ORG}_cluster"

echo "=== Fragment origin analysis: ${ORG} ==="
echo "Start: $(date)"
echo ""

# ── Step 1: Extract inserted fragments from Logan contigs ──
# A step is done for a signature of its inputs and options (stage_lib.sh):
# a FASTA from an earlier detection run was otherwise kept, and so were the
# representatives clustered from an earlier FASTA.
SIG1=$(stage_sig "extract fragments" "--min-len 300 --max-len 50000" "$(file_fp "$INS_TSV")" \
       "$(md5sum < "${SCRIPTS}/extract_insertion_fragments.py")")
if stage_done "$FRAG_FASTA.done" "$SIG1"; then
    echo "Step 1: Fragment FASTA exists for this insertion table, skipping."
else
    echo "Step 1: Extracting inserted fragment sequences..."
    rm -f "$FRAG_FASTA.done"
    python3 "${SCRIPTS}/extract_insertion_fragments.py" \
        "$INS_TSV" \
        "$LOGAN_DIR" \
        "$FRAG_FASTA.tmp" \
        --min-len 300 \
        --max-len 50000
    mv "$FRAG_FASTA.tmp" "$FRAG_FASTA"
    stage_mark "$FRAG_FASTA.done" "$SIG1"
fi
echo "  Fragments: $(grep -c '^>' "$FRAG_FASTA" 2>/dev/null || true)"
echo "  Size: $(du -h "$FRAG_FASTA" | cut -f1)"
echo "Step 1 done: $(date)"
echo ""

# ── Step 2: Cluster fragments at 95% identity with cd-hit-est ──
REP_FASTA="${CLUSTER_PREFIX}_rep.fasta"
CLSTR="${CLUSTER_PREFIX}_rep.fasta.clstr"

CDHIT_OPTS="-c 0.95 -n 10 -aS 0.9 -aL 0.9 -M 0 -d 0 -g 1"
SIG2=$(stage_sig "cd-hit-est" "$CDHIT_OPTS" "$SIG1" "$(file_fp "$FRAG_FASTA")")
if stage_done "$REP_FASTA.done" "$SIG2" "$REP_FASTA" "$CLSTR"; then
    echo "Step 2: Representatives exist for these fragments, skipping clustering."
else
    echo "Step 2: Clustering fragments at 95% identity..."
    rm -f "$REP_FASTA.done"
    if ! command -v cd-hit-est &>/dev/null; then
        echo "ERROR: cd-hit-est not found (conda install -c bioconda cd-hit)" >&2
        exit 1
    fi
    # 95% identity (the project threshold) over >= 90% of BOTH sequences:
    # without -aL/-aS a 300-bp fragment clustered into any 5-kb one containing it
    cd-hit-est -i "$FRAG_FASTA" -o "$REP_FASTA.tmp" $CDHIT_OPTS -T "${SLURM_CPUS_PER_TASK:-16}" \
        > "${CLUSTER_PREFIX}.log" 2>&1
    mv "$REP_FASTA.tmp.clstr" "$CLSTR"
    mv "$REP_FASTA.tmp" "$REP_FASTA"
    stage_mark "$REP_FASTA.done" "$SIG2"
fi

N_REPS=$(grep -c '^>' "$REP_FASTA" 2>/dev/null || true)
echo "  Representative clusters: $N_REPS"
echo "Step 2 done: $(date)"
echo ""

# ── Step 3: Summarize cluster sizes ──
echo "Step 3: Top cluster sizes (potential shared origins)..."
python3 - << PYEOF
import re
from collections import defaultdict

clusters = defaultdict(list)
current = None
with open("${CLSTR}") as f:
    for line in f:
        line = line.strip()
        if line.startswith(">Cluster"):
            current = line.split()[1]
        else:
            m = re.search(r'>(\S+)\.\.\.\s+(\*|at\s+)', line)
            if m:
                clusters[current].append(m.group(1))

sizes = [(c, len(m)) for c, m in clusters.items()]
sizes.sort(key=lambda x: -x[1])

print(f"Total clusters: {len(sizes):,}")
print(f"Singletons: {sum(1 for _, n in sizes if n == 1):,}")
print(f"Multi-member: {sum(1 for _, n in sizes if n > 1):,}")
print()
print("Top 20 clusters by size:")
for c, n in sizes[:20]:
    print(f"  Cluster {c}: {n} fragments")

with open("${CLUSTER_PREFIX}_sizes.tsv", "w") as out:
    out.write("cluster\tsize\n")
    for c, n in sizes:
        out.write(f"{c}\t{n}\n")
PYEOF
echo "Step 3 done: $(date)"
echo ""

# ── Step 4: Keep top N largest clusters for BLAST ──
echo "Step 4: Preparing top 50 cluster representatives for BLAST..."
python3 - << PYEOF
import re
from collections import defaultdict

# Re-parse cluster file to get rep->size mapping
clusters = defaultdict(list)
cluster_rep = {}
current = None
with open("${CLSTR}") as f:
    for line in f:
        line = line.strip()
        if line.startswith(">Cluster"):
            current = line.split()[1]
        else:
            m = re.search(r'>(\S+)\.\.\.\s+(\*|at\s+)', line)
            if m:
                name = m.group(1)
                clusters[current].append(name)
                if m.group(2) == '*':
                    cluster_rep[current] = name

# Sort clusters by size, keep top 50
sizes = [(c, len(m)) for c, m in clusters.items()]
sizes.sort(key=lambda x: -x[1])
top = {cluster_rep[c] for c, _ in sizes[:50] if c in cluster_rep}

# Filter the rep fasta down to top 50
written = 0
with open("${REP_FASTA}") as f, open("${CLUSTER_PREFIX}_top50.fasta", "w") as out:
    keep = False
    for line in f:
        if line.startswith(">"):
            name = line[1:].split()[0]
            keep = name in top
            if keep:
                # Add cluster size to header for context
                for c, n in sizes[:50]:
                    if cluster_rep.get(c) == name:
                        line = f">{name} cluster_size={n}\n"
                        break
                written += 1
        if keep:
            out.write(line)

print(f"Wrote {written} top cluster representatives to ${CLUSTER_PREFIX}_top50.fasta")
PYEOF
echo "Step 4 done: $(date)"
echo ""

echo "=== Fragment origin prep complete ==="
echo "Top 50 representatives ready for BLAST:"
echo "  ${CLUSTER_PREFIX}_top50.fasta"
echo "End: $(date)"
