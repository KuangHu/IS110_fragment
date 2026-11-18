#!/bin/bash
#SBATCH --job-name=extract_is110
#SBATCH --account=pc_jbei005
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/global/scratch/users/kh36969/logan_results/extract_is110_%j.out
#SBATCH --error=/global/scratch/users/kh36969/logan_results/extract_is110_%j.err

# Extract IS110 transposase + flanks from NCBI E. coli genomes,
# then cluster to define IS110 family groups.

set -euo pipefail
export PATH="$HOME/.conda/envs/claude-env/bin:$PATH"

POSITIONS_TSV="${POSITIONS_TSV:-/global/scratch/users/kh36969/logan_hmm/ncbi_ecoli/is110_positions.tsv}"
GENOME_DIR="${GENOME_DIR:-/global/scratch/users/kh36969/ncbi_refs/ecoli_complete/ecoli_complete_data/ncbi_dataset/data}"
OUTDIR="${OUTDIR:-/global/scratch/users/kh36969/is110_refs}"
# accession <TAB> genome file, written by hmm_ncbi_ecoli.sh next to the
# positions table: the genomes the positions were called on
MANIFEST="${MANIFEST:-$(dirname "$POSITIONS_TSV")/genomes_one_per_assembly.tsv}"
source "${IS110_FRAGMENT_SCRIPTS:-/global/home/users/kh36969/IS110_ecology/IS110_fragment/repo/scripts}/stage_lib.sh"
CPUS=${SLURM_CPUS_PER_TASK:-16}
FLANK_BP=2000

mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "============================================"
echo "Extract IS110 + ${FLANK_BP}bp flanks from NCBI E. coli"
echo "============================================"
echo "Start: $(date)"
echo ""

# ────────────────────────────────────────────────
# STEP 1: Extract IS110 + flanks per genome
# ────────────────────────────────────────────────
EXTRACTED="${OUTDIR}/is110_with_flanks.fa"
# the transposase gene alone (the "core"), same ids: this is what is clustered
CORE="${OUTDIR}/is110_core.fa"

# A step is done for a signature of its inputs and options (stage_lib.sh),
# not because its output exists: new positions left the old references, and
# new references left the old representatives.
SIG1=$(stage_sig "extract refs" "$FLANK_BP" "$(file_fp "$POSITIONS_TSV")" \
       "$( [ -f "$MANIFEST" ] && md5sum < "$MANIFEST" || echo no-manifest )" \
       "$( [ -f "$MANIFEST" ] && cut -f2 "$MANIFEST" | files_fp || echo "$GENOME_DIR" )")
if stage_done "$EXTRACTED.done" "$SIG1" "$EXTRACTED" "$CORE"; then
    echo "Step 1: extracted sequences exist for these inputs, skipping."
else
    echo "Step 1: Extracting IS110 + ${FLANK_BP}bp flanks..."
    python3 << PYEOF
import os
from collections import defaultdict

positions_tsv = "$POSITIONS_TSV"
genome_dir = "$GENOME_DIR"
manifest = "$MANIFEST"
allow_missing = "${ALLOW_MISSING_GENOMES:-0}" == "1"
out_path = "$EXTRACTED"
core_path = "$CORE"
flank = $FLANK_BP

# Group positions by assembly accession
by_asm = defaultdict(list)
with open(positions_tsv) as f:
    next(f)  # header
    for line in f:
        parts = line.rstrip().split("\t")
        if len(parts) < 7:
            continue
        protein_id, assembly, contig, start, end, strand, plen = parts
        by_asm[assembly].append({
            "protein_id": protein_id,
            "contig": contig,
            "start": int(start),
            "end": int(end),
            "strand": strand,
        })

print(f"  Assemblies with IS110: {len(by_asm):,}")

def revcomp(seq):
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]

# assembly -> genome file. From the manifest of the IS110 census when there
# is one (any directory layout); otherwise <genome_dir>/<assembly>/*_genomic.fna.
# The path used to be rebuilt from the accession only: with a flat layout
# every assembly was skipped and the script wrote an empty result.
genome_of = {}
if os.path.exists(manifest):
    with open(manifest) as f:
        for line in f:
            acc, path = line.rstrip("\n").split("\t")[:2]
            genome_of[acc] = path
def find_genome(assembly):
    if assembly in genome_of:
        return genome_of[assembly] if os.path.exists(genome_of[assembly]) else None
    asm_dir = os.path.join(genome_dir, assembly)
    if os.path.isdir(asm_dir):
        for f in sorted(os.listdir(asm_dir)):
            if f.endswith("_genomic.fna"):
                return os.path.join(asm_dir, f)
    return None

n_extracted = 0
missing_asm, missing_contig = [], 0
with open(out_path + ".tmp", "w") as out, open(core_path + ".tmp", "w") as core_out:
    for assembly, positions in sorted(by_asm.items()):
        fna = find_genome(assembly)
        if not fna:
            missing_asm.append(assembly)
            continue

        # Load contigs from this genome
        contigs = {}
        cur_name = None
        cur_parts = []
        with open(fna) as fh:
            for line in fh:
                if line.startswith(">"):
                    if cur_name:
                        contigs[cur_name] = "".join(cur_parts)
                    cur_name = line[1:].split()[0]
                    cur_parts = []
                else:
                    cur_parts.append(line.strip())
            if cur_name:
                contigs[cur_name] = "".join(cur_parts)

        # Extract each IS110 with flanks
        for p in positions:
            seq = contigs.get(p["contig"])
            if not seq:
                missing_contig += 1
                continue
            ext_start = max(0, p["start"] - 1 - flank)
            ext_end = min(len(seq), p["end"] + flank)
            extracted = seq[ext_start:ext_end]

            # Always orient transposase on + strand
            if p["strand"] == "-":
                extracted = revcomp(extracted)
                core_s, core_e = ext_end - p["end"], ext_end - (p["start"] - 1)
            else:
                core_s, core_e = (p["start"] - 1) - ext_start, p["end"] - ext_start

            # core=: the transposase gene inside this record (0-based half-open),
            # read by ref_core_coverage.py
            header = (f">{p['protein_id']} asm={assembly} contig={p['contig']} "
                      f"tnp_pos={p['start']}-{p['end']} strand={p['strand']} "
                      f"len={len(extracted)} flank={flank} core={core_s}-{core_e}")
            out.write(header + "\n")
            for i in range(0, len(extracted), 80):
                out.write(extracted[i:i+80] + "\n")
            core_out.write(f">{p['protein_id']}\n{extracted[core_s:core_e]}\n")
            n_extracted += 1

print(f"  Sequences extracted: {n_extracted:,}")
if missing_asm or missing_contig:
    print(f"  assemblies without a genome file: {len(missing_asm):,} "
          f"(e.g. {', '.join(missing_asm[:3])}); positions on a contig not in the file: {missing_contig:,}")
    if not allow_missing:
        os.remove(out_path + ".tmp"); os.remove(core_path + ".tmp")
        raise SystemExit("ERROR: positions without a sequence; nothing written "
                         "(ALLOW_MISSING_GENOMES=1 to write the rest)")
os.replace(out_path + ".tmp", out_path)
os.replace(core_path + ".tmp", core_path)
PYEOF
    stage_mark "$EXTRACTED.done" "$SIG1"
fi

N_EXT=$(grep -c '^>' "$EXTRACTED" || true)
echo "  Total: $N_EXT ($(du -h "$EXTRACTED" | cut -f1))"
echo "Step 1 done: $(date)"
echo ""

# ────────────────────────────────────────────────
# STEP 2: Cluster IS110 family groups
# ────────────────────────────────────────────────
# Cluster at 95% identity (strain-level variants, our standard threshold)
for IDENT in 0.95; do
    PCT=$(echo "$IDENT * 100" | bc | cut -d. -f1)
    REPS="${OUTDIR}/is110_reps_${PCT}.fa"
    CLUSTER="${OUTDIR}/is110_cluster_${PCT}.tsv"

    SIG2=$(stage_sig "cluster cores" "$IDENT" "-c 0.80 --cov-mode 0" "$SIG1" "$(file_fp "$CORE")")
    if stage_done "$REPS.done" "$SIG2" "$REPS" "$CLUSTER"; then
        echo "Step 2 (${PCT}% id): reps exist for these cores, skipping."
        continue
    fi
    rm -f "$REPS.done"
    rm -rf "${OUTDIR}/mmseqs_tmp_${PCT}"

    # Clustered on the transposase gene, not on gene + 4 kb of flank: with the
    # flanks, two identical IS in different places were two "families" and two
    # IS 10% apart in the same place were one.
    echo "Step 2 (${PCT}% id): Clustering the IS core with mmseqs2..."
    mmseqs easy-linclust \
        "$CORE" \
        "${OUTDIR}/cluster_${PCT}" \
        "${OUTDIR}/mmseqs_tmp_${PCT}" \
        --min-seq-id "$IDENT" -c 0.80 --cov-mode 0 \
        --threads "$CPUS" \
        > "${OUTDIR}/mmseqs_${PCT}.log" 2>&1
    mv "${OUTDIR}/cluster_${PCT}_cluster.tsv" "$CLUSTER"
    # representatives = the flanked records of the representative cores
    cut -f1 "$CLUSTER" | sort -u > "${OUTDIR}/cluster_${PCT}_rep_ids.txt"
    awk 'NR==FNR { keep[$1]=1; next }
         /^>/ { on = (substr($1,2) in keep) } on' \
        "${OUTDIR}/cluster_${PCT}_rep_ids.txt" "$EXTRACTED" > "$REPS.tmp"
    mv "$REPS.tmp" "$REPS"
    stage_mark "$REPS.done" "$SIG2"
    rm -f "${OUTDIR}/cluster_${PCT}_rep_seq.fasta"
    rm -rf "${OUTDIR}/mmseqs_tmp_${PCT}"
    rm -f "${OUTDIR}/cluster_${PCT}_all_seqs.fasta"

    N_REPS=$(grep -c '^>' "$REPS" || true)
    echo "  ${PCT}% id reps: $N_REPS"
done

echo ""
echo "Step 2 done: $(date)"
echo ""

# ────────────────────────────────────────────────
# STEP 3: Summarize cluster sizes
# ────────────────────────────────────────────────
echo "Step 3: Cluster size distributions..."
export OUTDIR
python3 << 'PYEOF'
import os
from collections import Counter
outdir = os.environ["OUTDIR"]

for pct in [95]:
    cluster_file = f"{outdir}/is110_cluster_{pct}.tsv"
    sizes = Counter()
    with open(cluster_file) as f:
        for line in f:
            rep, _ = line.rstrip().split("\t")
            sizes[rep] += 1

    print(f"\n  === {pct}% identity clusters ===")
    print(f"  Total clusters: {len(sizes):,}")
    size_dist = Counter(sizes.values())
    print(f"  Cluster size distribution (top 15):")
    for size, n_clusters in sorted(size_dist.items(), key=lambda x: -x[0])[:15]:
        print(f"    Size {size}: {n_clusters} clusters")

    # Top 10 largest clusters
    print(f"  Top 10 largest clusters:")
    for rep, size in sorted(sizes.items(), key=lambda x: -x[1])[:10]:
        print(f"    {rep[:30]}... -> {size:,} members")
PYEOF

echo ""
echo "End: $(date)"
