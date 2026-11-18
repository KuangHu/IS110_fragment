#!/bin/bash
#SBATCH --job-name=blast_frag
#SBATCH --account=pc_jbei005
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/global/scratch/users/kh36969/logan_results/blast_%j.out
#SBATCH --error=/global/scratch/users/kh36969/logan_results/blast_%j.err

# Remote-BLAST cluster representatives against NCBI nt to identify origin.
# Usage: sbatch scripts/blast_fragments.sh <organism>

set -euo pipefail

module load bio/blast-plus/2.14.1-gcc-11.4.0

ORG=${1:?Usage: sbatch blast_fragments.sh <organism>}
RESULTS_DIR="/global/scratch/users/kh36969/logan_results/${ORG}"
QUERY="${RESULTS_DIR}/fragments_${ORG}_cluster_top50.fasta"
OUT="${RESULTS_DIR}/fragments_${ORG}_blast.tsv"

if [ ! -f "$QUERY" ]; then
    echo "ERROR: $QUERY not found"
    exit 1
fi

echo "=== BLAST fragment cluster reps: ${ORG} ==="
echo "Query: $QUERY ($(grep -c '^>' "$QUERY") sequences)"
echo "Start: $(date)"
echo ""

# BLAST each sequence separately so one hang doesn't kill the batch
python3 - << PYEOF
import subprocess, os, time

queries = []
current = None
seq = []
with open("$QUERY") as f:
    for line in f:
        if line.startswith(">"):
            if current:
                queries.append((current, "".join(seq)))
            current = line[1:].strip()
            seq = []
        else:
            seq.append(line.strip())
    if current:
        queries.append((current, "".join(seq)))

print(f"Running {len(queries)} BLAST searches (max 3 hits each)...")

with open("$OUT", "w") as out:
    out.write("query\tsseqid\tstitle\tpident\tlength\tqcovs\tevalue\tbitscore\n")
    for i, (name, s) in enumerate(queries):
        query_file = f"/tmp/q_{os.getpid()}_{i}.fa"
        with open(query_file, "w") as qf:
            qf.write(f">{name}\n{s}\n")
        print(f"  [{i+1}/{len(queries)}] {name[:60]}... ({len(s)}bp)", flush=True)
        try:
            r = subprocess.run(
                ["blastn", "-query", query_file, "-db", "nt", "-remote",
                 "-evalue", "1e-20", "-max_target_seqs", "3",
                 "-outfmt", "6 sseqid stitle pident length qcovs evalue bitscore"],
                capture_output=True, text=True, timeout=300
            )
            if r.returncode != 0:
                # a network / remote failure is not "no hits"
                out.write(f"{name.split()[0]}\tERROR\t{r.stderr.strip()[:200]}\t-\t-\t-\t-\t-\n")
            elif r.stdout.strip():
                for ln in r.stdout.strip().split("\n"):
                    out.write(f"{name.split()[0]}\t{ln}\n")
            else:
                out.write(f"{name.split()[0]}\tNO_HITS\t-\t-\t-\t-\t-\t-\n")
            out.flush()
        except subprocess.TimeoutExpired:
            print(f"    TIMEOUT")
            out.write(f"{name.split()[0]}\tTIMEOUT\t-\t-\t-\t-\t-\t-\n")
        except Exception as e:
            print(f"    ERROR: {e}")
            out.write(f"{name.split()[0]}\tERROR\t{e}\t-\t-\t-\t-\t-\n")
        finally:
            try:
                os.remove(query_file)
            except OSError:
                pass
        time.sleep(2)  # be polite to NCBI

print("Done.")
PYEOF

echo ""
echo "Results: $OUT"
wc -l "$OUT"
echo "End: $(date)"
