#!/usr/bin/env python3
"""Filter a FASTA of IS candidates (any family) to IS110-only using two-domain
HMM search: PF01548 (DEDD) + PF02371 (Tnp20). One protein must hit BOTH.

Input:  multi-record FASTA of IS element sequences.
Output:
  is110_filtered.fa       — subset FASTA of IS110 elements
  is110_hits.tsv          — per-IS: is_id, seq_len, tnp orf coords, domains_hit, n_orfs

Two-domain rule matches Cross_reference_IS is_detect.py (set intersection).
"""
import argparse, os, subprocess, tempfile, shutil, sys
from collections import defaultdict
import pyrodigal


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--in-fa",   required=True, help="Input multi-record IS FASTA")
    p.add_argument("--hmm-dedd", default="/global/home/users/kh36969/Cross_reference_IS/hmm/PF01548.hmm")
    p.add_argument("--hmm-tnp20", default="/global/home/users/kh36969/Cross_reference_IS/hmm/PF02371.hmm")
    p.add_argument("--out-fa",  required=True, help="Output IS110-only FASTA")
    p.add_argument("--out-tsv", required=True, help="Per-IS hit summary TSV")
    p.add_argument("--threads", type=int, default=16)
    p.add_argument("--e-value", type=float, default=1e-5)
    p.add_argument("--any-domain", action="store_true",
                   help="Accept IS with ANY domain hit (default: require BOTH)")
    return p.parse_args()


def read_fasta(path):
    """Yield (header, seq) tuples."""
    hdr = None; seq_parts = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if hdr is not None:
                    yield hdr, "".join(seq_parts)
                hdr = line[1:].split()[0]  # first token
                seq_parts = []
            else:
                seq_parts.append(line)
        if hdr is not None:
            yield hdr, "".join(seq_parts)


def predict_orfs(sequences, threads):
    """Run pyrodigal in meta mode across all sequences, write proteins.
    Returns dict {protein_id: (parent_id, orf_num, start, end, strand)}."""
    orf_finder = pyrodigal.GeneFinder(meta=True)
    proteins = []
    orf_info = {}
    for parent_id, seq in sequences:
        genes = orf_finder.find_genes(seq.encode())
        for i, g in enumerate(genes):
            pid = f"{parent_id}_orf{i}"
            aa = g.translate()
            proteins.append((pid, aa))
            orf_info[pid] = (parent_id, i, g.begin, g.end, g.strand)
    return proteins, orf_info


def run_hmmscan(proteins_faa, hmm_db_path, out_domtbl, threads, e_value):
    """Run hmmscan; produce domain table (--domtblout)."""
    cmd = ["hmmscan",
           "--domtblout", out_domtbl,
           "--noali",
           "-E", str(e_value),
           "--cpu", str(threads),
           hmm_db_path,
           proteins_faa]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)


def parse_domtbl(path):
    """Parse hmmscan --domtblout. Return dict {protein_id: set(domain_names)}."""
    hits = defaultdict(set)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"): continue
            fields = line.split()
            if len(fields) < 22: continue
            target = fields[0]   # domain name (from HMM)
            query  = fields[3]   # protein id
            e_value = float(fields[6])
            hits[query].add(target)
    return hits


def main():
    args = parse_args()
    tmp = tempfile.mkdtemp(prefix="isfilt_")
    try:
        # 1. Read input FASTA
        print(f"Reading {args.in_fa} ...", file=sys.stderr, flush=True)
        seqs = list(read_fasta(args.in_fa))
        print(f"  {len(seqs):,} sequences", file=sys.stderr)

        # 2. Predict ORFs
        print(f"Predicting ORFs with pyrodigal (meta) ...", file=sys.stderr, flush=True)
        proteins, orf_info = predict_orfs(seqs, args.threads)
        print(f"  {len(proteins):,} predicted proteins", file=sys.stderr)

        proteins_faa = f"{tmp}/proteins.faa"
        with open(proteins_faa, "w") as fh:
            for pid, aa in proteins:
                fh.write(f">{pid}\n{aa}\n")

        # 3. Build combined HMM db (concatenate + hmmpress)
        print(f"Preparing combined HMM db ...", file=sys.stderr, flush=True)
        hmm_db = f"{tmp}/domains.hmm"
        with open(hmm_db, "w") as fout:
            for h in [args.hmm_dedd, args.hmm_tnp20]:
                with open(h) as fin: fout.write(fin.read())
        subprocess.run(["hmmpress", "-f", hmm_db], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # 4. Run hmmscan
        print(f"Running hmmscan ({args.threads} threads, E<{args.e_value}) ...",
              file=sys.stderr, flush=True)
        domtbl = f"{tmp}/hits.domtbl"
        run_hmmscan(proteins_faa, hmm_db, domtbl, args.threads, args.e_value)

        # 5. Parse hits
        hits_by_protein = parse_domtbl(domtbl)
        print(f"  {len(hits_by_protein):,} proteins with >=1 HMM hit",
              file=sys.stderr)

        # 6. Determine which IS elements pass
        # Domain names in the .hmm files (short names, not accessions)
        target_domains = set()
        for prots in hits_by_protein.values():
            target_domains |= prots
        print(f"  HMM target names present: {target_domains}",
              file=sys.stderr)

        # Group protein->IS
        is_domains = defaultdict(dict)  # is_id -> {orf_num: set of domains}
        for pid, doms in hits_by_protein.items():
            parent_id, orf_num, start, end, strand = orf_info[pid]
            is_domains[parent_id][pid] = {"doms": doms, "start": start, "end": end,
                                          "strand": strand, "orf_num": orf_num}

        # An IS passes if any ORF has hit in BOTH categories (any_domain=False).
        # Since we have exactly 2 HMMs, per-protein test is len(doms) >= 2.
        passing_is = {}
        for is_id, orfs in is_domains.items():
            best_orf = None
            for pid, info in orfs.items():
                doms = info["doms"]
                ok = len(doms) >= (1 if args.any_domain else len(target_domains))
                if not ok: continue
                if best_orf is None or (info["end"] - info["start"]) > (best_orf["end"] - best_orf["start"]):
                    best_orf = {**info, "protein_id": pid}
            if best_orf is not None:
                passing_is[is_id] = best_orf

        print(f"  IS110-passing (>=1 ORF with all {len(target_domains)} domains): "
              f"{len(passing_is):,} / {len(seqs):,} "
              f"({len(passing_is)*100/len(seqs):.1f}%)",
              file=sys.stderr)

        # 7. Write output FASTA + TSV
        seq_by_id = dict(seqs)
        with open(args.out_fa, "w") as fout:
            for is_id in sorted(passing_is):
                fout.write(f">{is_id} length={len(seq_by_id[is_id])}\n")
                seq = seq_by_id[is_id]
                for i in range(0, len(seq), 80):
                    fout.write(seq[i:i+80] + "\n")

        with open(args.out_tsv, "w") as fout:
            fout.write("is_id\tis_length\ttnp_orf_num\ttnp_start\ttnp_end\ttnp_strand\ttnp_len\tdomains_hit\tn_orfs\n")
            for is_id in sorted(passing_is):
                info = passing_is[is_id]
                n_orfs = len(is_domains[is_id])
                doms = ",".join(sorted(info["doms"]))
                fout.write(f"{is_id}\t{len(seq_by_id[is_id])}\t"
                           f"{info['orf_num']}\t{info['start']}\t{info['end']}\t"
                           f"{'+' if info['strand'] == 1 else '-'}\t"
                           f"{info['end']-info['start']}\t{doms}\t{n_orfs}\n")

        print(f"Wrote {args.out_fa}", file=sys.stderr)
        print(f"Wrote {args.out_tsv}", file=sys.stderr)

    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
