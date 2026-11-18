#!/bin/bash
#SBATCH --job-name=frag_local
#SBATCH --account=pc_jbei005
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/global/scratch/users/kh36969/logan_results/fraglocal_%j.out
#SBATCH --error=/global/scratch/users/kh36969/logan_results/fraglocal_%j.err

# Search fragment cluster reps against all Logan contigs (all organisms)
# to find where each fragment originates.
#
# Usage: sbatch scripts/find_fragment_origin_local.sh <organism>

set -euo pipefail

export PATH="$HOME/.conda/envs/claude-env/bin:$PATH"

ORG=${1:?Usage: sbatch find_fragment_origin_local.sh <organism>}
SCRATCH="/global/scratch/users/kh36969"
RESULTS_DIR="${SCRATCH}/logan_results/${ORG}"

QUERY="${RESULTS_DIR}/fragments_${ORG}_cluster_top50.fasta"
CPUS=${SLURM_CPUS_PER_TASK:-32}

if [ ! -f "$QUERY" ]; then
    echo "ERROR: $QUERY not found"
    exit 1
fi

echo "=== Local fragment origin search: ${ORG} ==="
echo "Query: $QUERY ($(grep -c '^>' "$QUERY") cluster reps)"
echo "Start: $(date)"
echo ""

# Search against all 5 organism contig databases
for target_org in ecoli efaecium kpneumoniae mavium belkanii; do
    TARGET_DB="${SCRATCH}/logan_results/${target_org}/all_contigs.fa"
    OUT_PAF="${RESULTS_DIR}/fragments_${ORG}_vs_${target_org}.paf"

    if [ ! -f "$TARGET_DB" ]; then
        echo "  ${target_org}: database not found, skipping"
        continue
    fi

    echo "--- Searching against ${target_org} ($(du -h "$TARGET_DB" | cut -f1)) ---"
    minimap2 -x asm5 -c --eqx \
        --secondary=yes -N 5 -p 0.8 \
        -t "$CPUS" \
        "$TARGET_DB" "$QUERY" \
        -o "$OUT_PAF" 2>/dev/null
    echo "  Hits: $(wc -l < "$OUT_PAF")"
    echo "  $(date)"
done

echo ""
echo "=== Summary: top hits per cluster representative ==="
python3 - << PYEOF
from collections import defaultdict

organisms = ['ecoli', 'efaecium', 'kpneumoniae', 'mavium', 'belkanii']
hits_per_query_per_org = defaultdict(lambda: defaultdict(int))
best_identity = defaultdict(dict)

for org in organisms:
    paf = "${RESULTS_DIR}/fragments_${ORG}_vs_" + org + ".paf"
    try:
        with open(paf) as f:
            for line in f:
                c = line.strip().split('\t')
                if len(c) < 12: continue
                qname = c[0]
                matches, block = int(c[9]), int(c[10])
                # identity over aligned columns: PAF column 11 counts gap columns, so a
                # copy with an indel fell below the threshold and was dropped
                _cg = next((t[5:] for t in c[12:] if t.startswith('cg:Z:')), '')
                _gaps = sum(int(n) for n, op in __import__('re').findall(r'(\d+)([IDN])', _cg))
                ident = matches / (block - _gaps) * 100 if block - _gaps > 0 else 0
                if ident < 90: continue  # high identity only
                hits_per_query_per_org[qname][org] += 1
                cur_best = best_identity[qname].get(org, 0)
                if ident > cur_best:
                    best_identity[qname][org] = ident
    except FileNotFoundError:
        continue

# Print summary
with open("${RESULTS_DIR}/fragments_${ORG}_origin_summary.tsv", "w") as out:
    out.write("query\tecoli_hits\tefaecium_hits\tkpneumoniae_hits\tmavium_hits\tbelkanii_hits\tbest_match_organism\tbest_identity\ttotal_hits\n")
    for qname in sorted(hits_per_query_per_org.keys()):
        counts = hits_per_query_per_org[qname]
        bests = best_identity[qname]
        best_org = max(bests.keys(), key=lambda o: bests[o]) if bests else 'none'
        best_id = bests.get(best_org, 0)
        total = sum(counts.values())
        out.write(f"{qname}\t{counts.get('ecoli',0)}\t{counts.get('efaecium',0)}\t{counts.get('kpneumoniae',0)}\t{counts.get('mavium',0)}\t{counts.get('belkanii',0)}\t{best_org}\t{best_id:.1f}\t{total}\n")

print(f"Summary written: ${RESULTS_DIR}/fragments_${ORG}_origin_summary.tsv")

# Print top findings sorted by total cross-organism hits
print()
print("Top 20 fragments by hit diversity (how widespread):")
print(f"{'query':<50} {'ec':>6} {'ef':>6} {'kp':>6} {'ma':>6} {'be':>6} {'best':<15} {'id%':>6}")
items = [(q, sum(hits_per_query_per_org[q].values()), hits_per_query_per_org[q], best_identity[q])
         for q in hits_per_query_per_org]
items.sort(key=lambda x: -x[1])
for q, total, counts, bests in items[:20]:
    best_org = max(bests.keys(), key=lambda o: bests[o]) if bests else 'none'
    best_id = bests.get(best_org, 0)
    print(f"{q[:50]:<50} {counts.get('ecoli',0):>6} {counts.get('efaecium',0):>6} "
          f"{counts.get('kpneumoniae',0):>6} {counts.get('mavium',0):>6} {counts.get('belkanii',0):>6} "
          f"{best_org:<15} {best_id:>6.1f}")
PYEOF

echo ""
echo "End: $(date)"
