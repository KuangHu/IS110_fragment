#!/bin/bash
#SBATCH --job-name=refs_in_logan
#SBATCH --account=pc_jbei005
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=/global/scratch/users/kh36969/logan_results/refs_in_logan_%j.out
#SBATCH --error=/global/scratch/users/kh36969/logan_results/refs_in_logan_%j.err

# Search the 2,945 NCBI-derived IS110 reference variants (with flanks) against
# Logan E. coli contigs to validate detection coverage.

set -euo pipefail
export PATH="$HOME/.conda/envs/claude-env/bin:$PATH"

REFS="/global/scratch/users/kh36969/is110_refs/is110_reps_95.fa"
LOGAN_FA="/global/scratch/users/kh36969/logan_results/ecoli/all_contigs.fa"
OUTDIR="/global/scratch/users/kh36969/is110_refs/ecoli_validation"
CPUS=${SLURM_CPUS_PER_TASK:-32}

mkdir -p "$OUTDIR"

echo "============================================"
echo "Search NCBI IS110 refs against Logan E. coli"
echo "============================================"
echo "Refs:       $REFS ($(grep -c '^>' "$REFS") sequences)"
echo "Logan:      $LOGAN_FA ($(du -h "$LOGAN_FA" | cut -f1))"
echo "Output:     $OUTDIR"
echo "Start:      $(date)"
echo ""

PAF="${OUTDIR}/refs_vs_logan.paf"

if [ -s "$PAF" ]; then
    echo "PAF exists, skipping."
else
    echo "Running minimap2 (asm10, sensitive)..."
# -N/-p: every copy, not the 100 best scoring >= half of the best one -- each
# half of an interrupted copy is a low-scoring secondary of any intact copy and
# was never reported (0 of 18 test configurations with an intact copy present).
# -g/-r 200: do not chain the two halves across the insert into one alignment.
    minimap2 -x asm10 -c --eqx \
        --secondary=yes -N 1000000 -p 0.05 -g 200 -r 200,200 \
        -t "$CPUS" \
        "$LOGAN_FA" "$REFS" \
        -o "$PAF"
fi

echo "  Total hits: $(wc -l < "$PAF")"
echo "Search done: $(date)"
echo ""

# Quick coverage analysis: how many of the 2,945 references have a hit in Logan?
echo "Coverage analysis..."
python3 << PYEOF
from collections import defaultdict

paf = "$PAF"
refs = "$REFS"

# Get all reference names
ref_names = set()
ref_lengths = {}
with open(refs) as f:
    for line in f:
        if line.startswith(">"):
            name = line[1:].split()[0]
            ref_names.add(name)
        else:
            if name in ref_lengths:
                ref_lengths[name] += len(line.strip())
            else:
                ref_lengths[name] = len(line.strip())

# Parse PAF: identity and best query coverage per query
best_coverage = {}  # qname -> best coverage %
hit_count = defaultdict(int)
high_id_count = defaultdict(int)
with open(paf) as f:
    for line in f:
        c = line.strip().split("\t")
        if len(c) < 12: continue
        qname = c[0]
        qlen = int(c[1])
        qs, qe = int(c[2]), int(c[3])
        matches, block = int(c[9]), int(c[10])
        # identity over aligned columns: PAF column 11 counts gap columns, so a
        # copy with an indel fell below the threshold and was dropped
        _cg = next((t[5:] for t in c[12:] if t.startswith('cg:Z:')), '')
        _gaps = sum(int(n) for n, op in __import__('re').findall(r'(\d+)([IDN])', _cg))
        ident = matches / (block - _gaps) * 100 if block - _gaps > 0 else 0
        cov = (qe - qs) / qlen * 100
        hit_count[qname] += 1
        if ident >= 95:
            high_id_count[qname] += 1
        if qname not in best_coverage or cov > best_coverage[qname]:
            best_coverage[qname] = cov

# Summary
total_refs = len(ref_names)
refs_with_any_hit = len(hit_count)
refs_with_high_id_hit = len(high_id_count)
refs_no_hit = ref_names - set(hit_count.keys())

# Coverage buckets
cov_buckets = [0, 0, 0, 0, 0]  # 0-20, 20-50, 50-80, 80-95, 95-100
for q, cov in best_coverage.items():
    if cov < 20: cov_buckets[0] += 1
    elif cov < 50: cov_buckets[1] += 1
    elif cov < 80: cov_buckets[2] += 1
    elif cov < 95: cov_buckets[3] += 1
    else: cov_buckets[4] += 1

print(f"\n=== Coverage of {total_refs:,} NCBI IS110 references in Logan E. coli ===")
print(f"  References with any hit:        {refs_with_any_hit:,} ({refs_with_any_hit/total_refs*100:.1f}%)")
print(f"  References with >=95% id hit:   {refs_with_high_id_hit:,} ({refs_with_high_id_hit/total_refs*100:.1f}%)")
print(f"  References with NO hit at all:  {len(refs_no_hit):,} ({len(refs_no_hit)/total_refs*100:.1f}%)")
print(f"\n  Best query-coverage distribution (any identity):")
print(f"    < 20%:    {cov_buckets[0]:,}")
print(f"    20-50%:   {cov_buckets[1]:,}")
print(f"    50-80%:   {cov_buckets[2]:,}")
print(f"    80-95%:   {cov_buckets[3]:,}")
print(f"    >= 95%:   {cov_buckets[4]:,}")

# Show some example references with NO hits (could be cargo artifacts in references, or truly missing in Logan)
if refs_no_hit:
    print(f"\n  Example refs with NO hit in Logan (top 5):")
    for r in list(refs_no_hit)[:5]:
        print(f"    {r}")
PYEOF

echo ""
echo "End: $(date)"
