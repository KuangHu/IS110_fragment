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
    p.add_argument("--max-protein-len", type=int, default=10000,
                   help="Drop proteins longer than N aa before hmmsearch (hmmer's "
                        "pipeline rejects sequences > 100K). Default 10000.")
    return p.parse_args()


def _pyrodigal_process_one(record_tuple):
    """Module-level worker so it can be pickled for multiprocessing.Pool."""
    import io
    import pyrodigal
    name, seq = record_tuple
    finder = pyrodigal.GeneFinder(meta=True)
    try:
        genes = finder.find_genes(seq.encode())
        faa_io = io.StringIO()
        gff_io = io.StringIO()
        genes.write_translations(faa_io, sequence_id=name)
        genes.write_gff(gff_io, sequence_id=name)
        return faa_io.getvalue(), gff_io.getvalue()
    except Exception:
        return "", ""


def _iter_fasta_records(db_fa):
    """Stream (name, seq) tuples from a FASTA without holding all in memory."""
    cur_name, cur_seq = None, []
    with open(db_fa) as f:
        for line in f:
            if line.startswith(">"):
                if cur_name is not None:
                    yield (cur_name, "".join(cur_seq))
                cur_name = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line.strip())
    if cur_name is not None:
        yield (cur_name, "".join(cur_seq))


def run_orfs_pyrodigal(db_fa, out_faa, out_gff, threads):
    """Run pyrodigal -meta in parallel via Python multiprocessing (streaming)."""
    from multiprocessing import Pool

    print(f"  ORF calling: streaming contigs with pyrodigal (n={threads})",
          file=sys.stderr, flush=True)
    with open(out_faa, "w") as ffaa, open(out_gff, "w") as fgff:
        with Pool(threads) as pool:
            for i, (faa_part, gff_part) in enumerate(
                    pool.imap_unordered(_pyrodigal_process_one,
                                        _iter_fasta_records(db_fa),
                                        chunksize=4)):
                ffaa.write(faa_part)
                fgff.write(gff_part)
                if i % 5000 == 0:
                    print(f"    contigs processed: {i:,}", file=sys.stderr, flush=True)


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


def filter_proteins_by_length(in_faa, out_faa, max_len):
    """Stream-filter a protein FASTA, dropping records longer than max_len aa."""
    kept = dropped = 0
    with open(in_faa) as fin, open(out_faa, "w") as fout:
        header = None
        seq_parts = []

        def flush():
            nonlocal kept, dropped
            if header is None:
                return
            seq = "".join(seq_parts)
            if len(seq) > max_len:
                dropped += 1
            else:
                fout.write(header)
                # Write seq in 60-char lines (already is; just emit as-is)
                fout.write(seq + "\n")
                kept += 1

        for line in fin:
            if line.startswith(">"):
                flush()
                header = line
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        flush()
    return kept, dropped


def main():
    args = parse_args()
    work_dir = args.work_dir or args.out + "_work"
    os.makedirs(work_dir, exist_ok=True)

    proteins_faa = os.path.join(work_dir, "proteins.faa")
    proteins_gff = os.path.join(work_dir, "proteins.gff")
    orf_done_marker = os.path.join(work_dir, ".orf_complete")

    # Step 1: ORFs. Cache only if a completion MARKER exists — a partial/killed
    # pyrodigal leaves a non-empty proteins.faa that would otherwise be mistaken
    # for a finished run (this silently truncated a 2.3M-contig DB to 15K once).
    if not args.skip_orf:
        if os.path.exists(orf_done_marker) and os.path.getsize(proteins_faa) > 0:
            print(f"Step 1: completed ORF outputs found (marker present), "
                  f"skipping pyrodigal.", file=sys.stderr, flush=True)
        else:
            if os.path.exists(proteins_faa):
                print(f"Step 1: existing proteins.faa has no completion marker "
                      f"(partial/stale) — rerunning pyrodigal.", file=sys.stderr,
                      flush=True)
            print(f"Step 1: ORF prediction ({args.orf_caller})...", file=sys.stderr)
            if args.orf_caller == "pyrodigal":
                run_orfs_pyrodigal(args.db, proteins_faa, proteins_gff, args.threads)
            else:
                cmd = ["prodigal", "-p", "meta", "-i", args.db,
                       "-a", proteins_faa, "-o", proteins_gff, "-f", "gff"]
                subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
            # mark completion only after the ORF caller returns successfully
            open(orf_done_marker, "w").close()

    # Step 1b: filter proteins for hmmsearch (drop sequences > --max-protein-len)
    # hmmsearch hard-rejects sequences > 100K aa; we use a tighter threshold by default.
    filtered_faa = os.path.join(work_dir, "proteins_filtered.faa")
    if not (os.path.exists(filtered_faa) and os.path.getsize(filtered_faa) > 0):
        print(f"Step 1b: filtering proteins (max length {args.max_protein_len:,} aa) ...",
              file=sys.stderr, flush=True)
        kept, dropped = filter_proteins_by_length(proteins_faa, filtered_faa,
                                                  args.max_protein_len)
        print(f"  kept {kept:,}, dropped {dropped:,} oversized proteins",
              file=sys.stderr, flush=True)

    # Step 2: HMM searches (use the filtered file)
    hmm_files = args.hmm.split(",")
    hmm_hits = []  # list of (hmm_name, set_of_protein_ids)
    for hmm in hmm_files:
        hmm_name = os.path.basename(hmm).replace(".hmm", "")
        domtbl = os.path.join(work_dir, f"{hmm_name}.domtbl")
        if not os.path.exists(domtbl):
            print(f"Step 2: hmmsearch {hmm_name} ...", file=sys.stderr)
            run_hmmsearch(hmm, filtered_faa, domtbl, args.e_value, args.threads)
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
