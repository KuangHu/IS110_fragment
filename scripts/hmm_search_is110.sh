#!/bin/bash
#SBATCH --job-name=hmm_IS110
#SBATCH --partition=standard
#SBATCH --cpus-per-task=48
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/groups/rubin/projects/kuang/out/IS110_fragment/hmm_search/hmm_IS110_%j.out
#SBATCH --error=/groups/rubin/projects/kuang/out/IS110_fragment/hmm_search/hmm_IS110_%j.err

# Pipeline: Find all IS110 transposases in GTDB using Pfam HMMs,
# then search them back against GTDB to find split CDS cases.

set -euo pipefail

GTDB_FNA="/groups/rubin/databases/GTDB/gtdb_genomes_v2024latest/gtdb_genomes_reps_r220/allfilestogether"
GTDB_FAA="/groups/rubin/projects/kuang/out/embedding_project/protein_collection/GTDB/GTDB_faa"
HMM_DIR="/home/kuangh/tools/IS110_fragment/hmm"
OUTDIR="/groups/rubin/projects/kuang/out/IS110_fragment/hmm_search"
CPUS=48

mkdir -p "$OUTDIR"

echo "=== IS110 HMM search pipeline ==="
echo "Start: $(date)"
echo ""

# ──────────────────────────────────────────────
# Step 1: Concatenate pre-computed GTDB proteins
# ──────────────────────────────────────────────
PROTEIN_FILE="$OUTDIR/all_proteins.faa"

if [ -f "$PROTEIN_FILE" ] && [ "$(stat -c%s "$PROTEIN_FILE" 2>/dev/null || echo 0)" -gt 80000000000 ]; then
    echo "Step 1: Protein file already exists, skipping."
else
    echo "Step 1: Concatenating pre-computed GTDB proteins..."
    rm -f "$PROTEIN_FILE"
    find "$GTDB_FAA" -name "*.faa" -print0 | xargs -0 cat > "$PROTEIN_FILE"
    echo "  File size: $(du -h "$PROTEIN_FILE" | cut -f1)"
fi
echo "Step 1 done: $(date)"
echo ""

# ──────────────────────────────────────────────
# Step 2: hmmsearch with both Pfam HMMs
# ──────────────────────────────────────────────
echo "Step 2: Running hmmsearch..."

DEDD_HITS="$OUTDIR/dedd_hits.domtbl"
TNP20_HITS="$OUTDIR/tnp20_hits.domtbl"

if [ ! -f "$DEDD_HITS" ] || [ ! -s "$DEDD_HITS" ]; then
    hmmsearch --domtblout "$DEDD_HITS" \
        --cpu "$CPUS" -E 1e-5 --domE 1e-5 \
        "$HMM_DIR/PF01548.hmm" "$PROTEIN_FILE" > /dev/null
fi
echo "  DEDD hits: $(grep -cv '^#' "$DEDD_HITS") domain hits"

if [ ! -f "$TNP20_HITS" ] || [ ! -s "$TNP20_HITS" ]; then
    hmmsearch --domtblout "$TNP20_HITS" \
        --cpu "$CPUS" -E 1e-5 --domE 1e-5 \
        "$HMM_DIR/PF02371.hmm" "$PROTEIN_FILE" > /dev/null
fi
echo "  Tnp20 hits: $(grep -cv '^#' "$TNP20_HITS") domain hits"

echo "Step 2 done: $(date)"
echo ""

# ──────────────────────────────────────────────
# Step 3: Find proteins with BOTH domains
# ──────────────────────────────────────────────
echo "Step 3: Finding proteins with both DEDD + Tnp20 domains..."

python3 - << 'PYEOF'
import os
outdir = "/groups/rubin/projects/kuang/out/IS110_fragment/hmm_search"

def parse_domtbl(path):
    hits = set()
    with open(path) as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.split()
            if len(parts) >= 1:
                hits.add(parts[0])
    return hits

dedd = parse_domtbl(os.path.join(outdir, "dedd_hits.domtbl"))
tnp20 = parse_domtbl(os.path.join(outdir, "tnp20_hits.domtbl"))
both = dedd & tnp20

print(f"  DEDD domain proteins:  {len(dedd)}")
print(f"  Tnp20 domain proteins: {len(tnp20)}")
print(f"  Both domains (IS110):  {len(both)}")

with open(os.path.join(outdir, "is110_protein_ids.txt"), "w") as f:
    for pid in sorted(both):
        f.write(pid + "\n")
PYEOF

echo "Step 3 done: $(date)"
echo ""

# ──────────────────────────────────────────────
# Step 4: Extract CDS nucleotide sequences
# ──────────────────────────────────────────────
echo "Step 4: Extracting IS110 CDS nucleotide sequences..."

python3 - << 'PYEOF'
import os
from collections import defaultdict
from Bio import SeqIO

outdir = "/groups/rubin/projects/kuang/out/IS110_fragment/hmm_search"
gtdb_faa = "/groups/rubin/projects/kuang/out/embedding_project/protein_collection/GTDB/GTDB_faa"
gtdb_fna = "/groups/rubin/databases/GTDB/gtdb_genomes_v2024latest/gtdb_genomes_reps_r220/allfilestogether"

with open(os.path.join(outdir, "is110_protein_ids.txt")) as f:
    target_ids = set(line.strip() for line in f if line.strip())

if not target_ids:
    print("  No IS110 proteins found!")
    exit()

print(f"  IS110 proteins to extract: {len(target_ids)}")

# Scan .faa files to find coordinates for each target protein
protein_info = {}
found_in_genome = defaultdict(list)

for faa_file in os.listdir(gtdb_faa):
    if not faa_file.endswith(".faa"):
        continue
    genome_name = faa_file.replace(".faa", "")
    faa_path = os.path.join(gtdb_faa, faa_file)

    with open(faa_path) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            protein_id = header.split()[0]
            if protein_id not in target_ids:
                continue

            parts = header.split(" # ")
            if len(parts) < 4:
                continue

            contig = parts[0].rsplit("_", 1)[0]
            start = int(parts[1])
            end = int(parts[2])
            strand = int(parts[3])

            protein_info[protein_id] = {
                "genome": genome_name,
                "contig": contig,
                "start": start,
                "end": end,
                "strand": strand,
            }
            found_in_genome[genome_name].append(protein_id)

    if len(protein_info) == len(target_ids):
        break

print(f"  Found coordinates for {len(protein_info)} / {len(target_ids)} proteins")
print(f"  Across {len(found_in_genome)} genomes")

# Extract nucleotide sequences
out_file = os.path.join(outdir, "is110_cds_nuc.fna")
count = 0
skipped_short = 0
with open(out_file, "w") as fout:
    for genome_name, pids in found_in_genome.items():
        parts = genome_name.split("_")
        accession = f"{parts[0]}_{parts[1]}"
        fna_path = os.path.join(gtdb_fna, f"{accession}_genomic.fna")

        if not os.path.exists(fna_path):
            print(f"  WARNING: genome not found: {fna_path}")
            continue

        contigs = {}
        for rec in SeqIO.parse(fna_path, "fasta"):
            contigs[rec.id] = rec.seq

        for pid in pids:
            info = protein_info[pid]
            contig_seq = contigs.get(info["contig"])
            if contig_seq is None:
                continue

            cds_seq = contig_seq[info["start"]-1:info["end"]]
            if info["strand"] == -1:
                cds_seq = cds_seq.reverse_complement()

            if len(cds_seq) < 800:
                skipped_short += 1
                continue

            fout.write(f">{pid}__{info['genome']}__{info['contig']}_{info['start']}_{info['end']}_{'+' if info['strand']==1 else '-'}\n")
            fout.write(str(cds_seq) + "\n")
            count += 1

print(f"  Extracted {count} IS110 CDS (>= 800bp)")
print(f"  Skipped {skipped_short} short (<800bp) sequences")
PYEOF

echo "Step 4 done: $(date)"
echo ""

# ──────────────────────────────────────────────
# Step 5: Dereplicate with cd-hit at 95% identity
# ──────────────────────────────────────────────
echo "Step 5: Dereplicating with cd-hit-est (95% identity)..."

cd-hit-est -i "$OUTDIR/is110_cds_nuc.fna" \
    -o "$OUTDIR/is110_cds_nr95.fna" \
    -c 0.95 -n 10 -T "$CPUS" -M 0 -d 0 \
    > "$OUTDIR/cdhit.log" 2>&1

NR=$(grep -c '^>' "$OUTDIR/is110_cds_nr95.fna")
echo "  Non-redundant IS110 CDS (95%): $NR sequences"

echo "Step 5 done: $(date)"
echo ""

# ──────────────────────────────────────────────
# Step 6: minimap2 search back against GTDB
# ──────────────────────────────────────────────
echo "Step 6: Searching dereplicated IS110 CDS against GTDB..."

SEARCH_OUT="$OUTDIR/is110_vs_gtdb.paf"

find "$GTDB_FNA" -name "*.fna" -print0 | xargs -0 cat \
    | minimap2 -x asm10 -c --eqx --secondary=yes -N 100 -p 0.5 \
        -t "$CPUS" \
        /dev/stdin "$OUTDIR/is110_cds_nr95.fna" \
        -o "$SEARCH_OUT"

echo "  Total hits: $(wc -l < "$SEARCH_OUT")"
echo "Step 6 done: $(date)"
echo ""

# ──────────────────────────────────────────────
# Step 7: Detect split CDS
# ──────────────────────────────────────────────
echo "Step 7: Detecting split IS110 transposases..."

python3 - << 'PYEOF'
import os
from collections import defaultdict

paf = "/groups/rubin/projects/kuang/out/IS110_fragment/hmm_search/is110_vs_gtdb.paf"
MIN_GAP = 6000
MIN_IDENT = 80

hits = defaultdict(list)
with open(paf) as fh:
    for line in fh:
        c = line.strip().split("\t")
        if len(c) < 12: continue
        matches, block = int(c[9]), int(c[10])
        ident = matches / block * 100 if block > 0 else 0
        if ident < MIN_IDENT: continue
        hits[(c[0], c[5])].append({
            "qs": int(c[2]), "qe": int(c[3]), "qlen": int(c[1]),
            "ts": int(c[7]), "te": int(c[8]), "tlen": int(c[6]),
            "strand": c[4], "ident": ident,
        })

count = 0
results = []
for (qname, tname), hlist in hits.items():
    if len(hlist) < 2: continue
    hlist.sort(key=lambda x: x["ts"])
    for i in range(len(hlist)):
        for j in range(i+1, len(hlist)):
            hi, hj = hlist[i], hlist[j]
            tgap = hj["ts"] - hi["te"]
            if tgap < MIN_GAP: continue
            qov = max(0, min(hi["qe"], hj["qe"]) - max(hi["qs"], hj["qs"]))
            min_span = min(hi["qe"]-hi["qs"], hj["qe"]-hj["qs"])
            if qov > min_span * 0.5: continue
            total_cov = len(set(range(hi["qs"], hi["qe"])) | set(range(hj["qs"], hj["qe"])))
            cov_pct = total_cov / hi["qlen"] * 100
            count += 1
            results.append((qname, tname, hi, hj, tgap, cov_pct))
            break
        else: continue
        break

print(f"\nTotal split IS110 transposase cases: {count}")
print("")
for qname, tname, hi, hj, tgap, cov_pct in sorted(results, key=lambda x: -x[5])[:100]:
    print(f"  {qname} on {tname} ({hj['tlen']:,}bp)")
    print(f"    Part1: Q[{hi['qs']}-{hi['qe']}]/{hi['qlen']} -> T[{hi['ts']:,}-{hi['te']:,}] {hi['strand']} {hi['ident']:.0f}%")
    print(f"    Part2: Q[{hj['qs']}-{hj['qe']}]/{hj['qlen']} -> T[{hj['ts']:,}-{hj['te']:,}] {hj['strand']} {hj['ident']:.0f}%")
    print(f"    Gap: {tgap:,}bp | CDS covered: {cov_pct:.0f}%")

outdir = "/groups/rubin/projects/kuang/out/IS110_fragment/hmm_search"
with open(os.path.join(outdir, "split_is110_results.tsv"), "w") as f:
    f.write("query\ttarget\ttlen\tpart1_qs\tpart1_qe\tpart1_ts\tpart1_te\tpart1_strand\tpart1_ident\tpart2_qs\tpart2_qe\tpart2_ts\tpart2_te\tpart2_strand\tpart2_ident\tgap\tcds_coverage\n")
    for qname, tname, hi, hj, tgap, cov_pct in results:
        f.write(f"{qname}\t{tname}\t{hj['tlen']}\t{hi['qs']}\t{hi['qe']}\t{hi['ts']}\t{hi['te']}\t{hi['strand']}\t{hi['ident']:.1f}\t{hj['qs']}\t{hj['qe']}\t{hj['ts']}\t{hj['te']}\t{hj['strand']}\t{hj['ident']:.1f}\t{tgap}\t{cov_pct:.1f}\n")

PYEOF

echo ""
echo "=== Pipeline complete ==="
echo "End: $(date)"
