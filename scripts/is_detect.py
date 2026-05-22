#!/usr/bin/env python3
"""IS detection by HMM search.

Given a genome database FASTA and one or more HMM profiles, find candidate IS
positions (transposase ORF coordinates + strand + assembly + contig).

Workflow:
  1. ORF prediction (pyrodigal or prodigal -p meta)
  2. hmmsearch each HMM profile against predicted proteins
  3. Keep proteins passing ALL required HMM domains (default)
     OR ANY (set --any-domain)
  4. Map protein IDs back to nucleotide coordinates (contig + start/end/strand)

Output: TSV
  is_id  assembly  contig  tnp_start  tnp_end  tnp_strand  tnp_len  domains_hit

Usage:
  is_detect.py --db genome_db.fa --hmm PF01548.hmm,PF02371.hmm --out is_hits.tsv \
               [--any-domain] [--e-value 1e-5] [--threads 32]
"""
import argparse, csv, os, subprocess, sys, re
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--db", required=True, help="Target genome DB FASTA")
    p.add_argument("--hmm", required=True,
                   help="Comma-separated HMM profile paths (.hmm files)")
    p.add_argument("--out", required=True, help="Output TSV path")
    p.add_argument("--work-dir", default=None,
                   help="Working directory (proteins, GFF, domtbl). Defaults to <out>_work/")
    p.add_argument("--any-domain", action="store_true",
                   help="Keep proteins matching ANY of the HMMs (default: ALL)")
    p.add_argument("--e-value", type=float, default=1e-5)
    p.add_argument("--threads", type=int, default=32)
    p.add_argument("--orf-caller", default="pyrodigal",
                   choices=["pyrodigal", "prodigal"])
    p.add_argument("--skip-orf", action="store_true",
                   help="Skip ORF prediction (expects <work_dir>/proteins.faa to exist)")
    return p.parse_args()


def run_orfs_pyrodigal(db_fa, out_faa, out_gff, threads):
    """Run pyrodigal -meta in parallel via Python multiprocessing."""
    import pyrodigal
    from multiprocessing import Pool

    def process_one(record_tuple):
        name, seq = record_tuple
        finder = pyrodigal.GeneFinder(meta=True)
        try:
            genes = finder.find_genes(seq.encode())
            faa_lines = []
            gff_lines = []
            import io
            faa_io = io.StringIO()
            gff_io = io.StringIO()
            genes.write_translations(faa_io, sequence_id=name)
            genes.write_gff(gff_io, sequence_id=name)
            return faa_io.getvalue(), gff_io.getvalue()
        except Exception as e:
            return "", ""

    # Read FASTA into list
    records = []
    cur_name, cur_seq = None, []
    with open(db_fa) as f:
        for line in f:
            if line.startswith(">"):
                if cur_name:
                    records.append((cur_name, "".join(cur_seq)))
                cur_name = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line.strip())
    if cur_name:
        records.append((cur_name, "".join(cur_seq)))

    print(f"  ORF calling: {len(records):,} contigs with pyrodigal (n={threads})",
          file=sys.stderr, flush=True)
    with open(out_faa, "w") as ffaa, open(out_gff, "w") as fgff:
        with Pool(threads) as pool:
            for i, (faa_part, gff_part) in enumerate(pool.imap_unordered(process_one, records)):
                ffaa.write(faa_part)
                fgff.write(gff_part)
                if i % 1000 == 0:
                    print(f"    {i}/{len(records)}", file=sys.stderr, flush=True)


def run_hmmsearch(hmm, faa, domtbl, e_value, threads):
    """Run hmmsearch and emit domtbl."""
    cmd = ["hmmsearch", "--domtblout", domtbl, "--cpu", str(threads),
           "-E", str(e_value), "--domE", str(e_value), hmm, faa]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)


def parse_domtbl(domtbl):
    """Yield protein IDs that have a hit in this HMM."""
    hits = set()
    with open(domtbl) as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.split()
            if not parts: continue
            hits.add(parts[0])
    return hits


def parse_gff(gff_path):
    """Parse Prodigal/Pyrodigal GFF: protein_id -> (contig, start, end, strand)."""
    info = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            parts = line.rstrip().split("\t")
            if len(parts) < 9: continue
            if parts[2] != "CDS": continue
            contig = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            # Extract ID=...
            m = re.search(r"ID=([^;]+)", parts[8])
            if not m: continue
            pid = m.group(1)
            # pyrodigal uses sequence_id_N for protein ID in faa
            # check if pid format matches the faa header convention
            info[pid] = (contig, start, end, strand)
    return info


def main():
    args = parse_args()
    work_dir = args.work_dir or args.out + "_work"
    os.makedirs(work_dir, exist_ok=True)

    proteins_faa = os.path.join(work_dir, "proteins.faa")
    proteins_gff = os.path.join(work_dir, "proteins.gff")

    # Step 1: ORFs
    if not args.skip_orf:
        print(f"Step 1: ORF prediction ({args.orf_caller})...", file=sys.stderr)
        if args.orf_caller == "pyrodigal":
            run_orfs_pyrodigal(args.db, proteins_faa, proteins_gff, args.threads)
        else:
            cmd = ["prodigal", "-p", "meta", "-i", args.db,
                   "-a", proteins_faa, "-o", proteins_gff, "-f", "gff"]
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)

    # Step 2: HMM searches
    hmm_files = args.hmm.split(",")
    hmm_hits = []  # list of (hmm_name, set_of_protein_ids)
    for hmm in hmm_files:
        hmm_name = os.path.basename(hmm).replace(".hmm", "")
        domtbl = os.path.join(work_dir, f"{hmm_name}.domtbl")
        if not os.path.exists(domtbl):
            print(f"Step 2: hmmsearch {hmm_name} ...", file=sys.stderr)
            run_hmmsearch(hmm, proteins_faa, domtbl, args.e_value, args.threads)
        hits = parse_domtbl(domtbl)
        hmm_hits.append((hmm_name, hits))
        print(f"  {hmm_name}: {len(hits):,} hits", file=sys.stderr)

    # Step 3: Combine
    if args.any_domain:
        combined = set().union(*(h for _, h in hmm_hits))
    else:
        combined = set.intersection(*(h for _, h in hmm_hits)) if hmm_hits else set()
    print(f"  Combined ({'ANY' if args.any_domain else 'ALL'} domains): "
          f"{len(combined):,}", file=sys.stderr)

    # Step 4: Map back to contig coordinates via GFF
    print(f"Step 3: parsing GFF for coordinates ...", file=sys.stderr)
    gff_info = parse_gff(proteins_gff)
    print(f"  GFF entries: {len(gff_info):,}", file=sys.stderr)

    # Write output TSV
    domain_names = ",".join(n for n, _ in hmm_hits)
    with open(args.out, "w") as fout:
        fout.write("is_id\tassembly\tcontig\ttnp_start\ttnp_end\ttnp_strand"
                   "\ttnp_len\tdomains_hit\n")
        n_out = 0
        for pid in sorted(combined):
            if pid not in gff_info: continue
            contig, start, end, strand = gff_info[pid]
            tnp_len = end - start + 1
            # Build is_id from contig + protein index (last underscore part)
            # pyrodigal protein IDs look like CONTIG_NNN
            # Try to extract assembly from contig if pipe-delimited
            if "|" in contig:
                assembly = contig.split("|")[0]
            else:
                assembly = contig
            is_id = f"{contig}_{pid.rsplit('_', 1)[-1]}"
            fout.write(f"{is_id}\t{assembly}\t{contig}\t{start}\t{end}"
                       f"\t{strand}\t{tnp_len}\t{domain_names}\n")
            n_out += 1
        print(f"\nWrote {n_out:,} IS hits to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
