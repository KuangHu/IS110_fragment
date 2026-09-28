#!/usr/bin/env python3
"""Cluster long IS110 bodies in a species at <identity>% nucleotide identity,
output cluster representatives + their is_hits.tsv subset.

Inputs:
  --is-hits   per-species is_hits.tsv (assembly, contig, tnp_start, tnp_end, tnp_strand, tnp_len, domains_hit)
  --db        species concatenated FASTA, indexed with .fai
  --min-tnp-len 1000 bp (default — "long" IS110)
  --identity    0.99 (default)

Outputs (under --out):
  long_is_bodies.fa            all long IS body sequences (one record per is_id)
  reps.clstr                   cd-hit-est cluster file
  reps.fa                      one rep per cluster
  reps_is_hits.tsv             is_hits.tsv subset containing only the cluster reps
  rep_ref_assemblies.txt       distinct ref assemblies that contribute a rep
  cluster_summary.tsv          per-cluster: rep_id, n_members, member_asms
"""
import argparse, csv, os, subprocess, sys
from collections import defaultdict
import pysam


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--is-hits", required=True)
    p.add_argument("--db",      required=True, help="species FASTA (indexed)")
    p.add_argument("--out",     required=True)
    p.add_argument("--min-tnp-len",     type=int,   default=1000)
    p.add_argument("--require-both-domains", action="store_true", default=True)
    p.add_argument("--is-pad",          type=int,   default=50)
    p.add_argument("--identity",        type=float, default=0.99)
    p.add_argument("--coverage",        type=float, default=0.95,
                   help="cd-hit -aS (alignment coverage of shorter seq)")
    p.add_argument("--threads",         type=int,   default=8)
    p.add_argument("--word-size",       type=int,   default=10,
                   help="cd-hit -n (10 ok for ≥0.95 identity)")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)
    fa = pysam.FastaFile(args.db)

    # 1. Load long IS records
    rows = []
    with open(args.is_hits) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try: tlen = int(r["tnp_len"])
            except ValueError: continue
            if tlen < args.min_tnp_len: continue
            dom = r.get("domains_hit", "")
            if args.require_both_domains and not ("PF01548" in dom and "PF02371" in dom):
                continue
            rows.append(r)
    print(f"[cluster] {len(rows)} long IS records (≥{args.min_tnp_len} bp + both Pfam domains)",
          flush=True)

    # 2. Extract IS body FASTA (with is_pad on each side)
    fa_path = os.path.join(args.out, "long_is_bodies.fa")
    is_id_to_record = {}
    n_written = 0
    with open(fa_path, "w") as out:
        for r in rows:
            asm = r["assembly"]; contig = r["contig"]
            ts, te = int(r["tnp_start"]), int(r["tnp_end"])
            try: clen = fa.get_reference_length(contig)
            except (KeyError, ValueError):
                alt = f"{asm}|{contig}"
                try: clen = fa.get_reference_length(alt); contig = alt
                except (KeyError, ValueError): continue
            is_s = max(0, ts - args.is_pad); is_e = min(clen, te + args.is_pad)
            seq = fa.fetch(contig, is_s, is_e)
            if not seq: continue
            is_id = r["is_id"]
            out.write(f">{is_id}\n")
            for i in range(0, len(seq), 80): out.write(seq[i:i+80] + "\n")
            is_id_to_record[is_id] = r
            n_written += 1
    fa.close()
    print(f"[cluster] wrote {n_written} IS bodies to {fa_path}", flush=True)
    if n_written == 0:
        sys.exit("No long IS to cluster.")

    # 3. cd-hit-est
    reps_fa  = os.path.join(args.out, "reps.fa")
    cmd = ["cd-hit-est", "-i", fa_path, "-o", reps_fa,
           "-c", str(args.identity), "-aS", str(args.coverage),
           "-n", str(args.word_size), "-T", str(args.threads),
           "-M", "16000", "-d", "0"]
    print(f"[cluster] cd-hit-est: {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    clstr_path = reps_fa + ".clstr"
    print(f"[cluster] cd-hit done. Reps: {reps_fa}", flush=True)

    # 4. Parse cluster file → cluster_id → (rep_id, members[])
    clusters = []
    cur_cluster = None
    with open(clstr_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">Cluster"):
                if cur_cluster: clusters.append(cur_cluster)
                cur_cluster = {"id": line.split()[1], "rep": None, "members": []}
            else:
                # Format: 0    1234nt, >some_id... *   (or)   ... at +/99.10%
                parts = line.split(", ")
                if len(parts) < 2: continue
                meta = parts[1]   # ">some_id... *" or ">some_id... at..."
                seq_id = meta.split(">", 1)[1].split("...")[0]
                cur_cluster["members"].append(seq_id)
                if "*" in line:
                    cur_cluster["rep"] = seq_id
        if cur_cluster: clusters.append(cur_cluster)
    print(f"[cluster] {len(clusters)} clusters at {args.identity*100:.0f}% identity", flush=True)

    # 5. Write reps-only is_hits.tsv
    reps_tsv = os.path.join(args.out, "reps_is_hits.tsv")
    fieldnames = list(rows[0].keys())
    with open(reps_tsv, "w") as fh:
        fh.write("\t".join(fieldnames) + "\n")
        for cl in clusters:
            rep_id = cl["rep"]
            if rep_id not in is_id_to_record: continue
            r = is_id_to_record[rep_id]
            fh.write("\t".join(r[k] for k in fieldnames) + "\n")
    print(f"[cluster] wrote reps_is_hits.tsv with {len(clusters)} reps", flush=True)

    # 6. Distinct ref assemblies
    ref_asms = set()
    with open(reps_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ref_asms.add(r["assembly"])
    with open(os.path.join(args.out, "rep_ref_assemblies.txt"), "w") as fh:
        for a in sorted(ref_asms): fh.write(a + "\n")
    print(f"[cluster] {len(ref_asms)} distinct ref assemblies host the reps", flush=True)

    # 7. Cluster summary
    summary = os.path.join(args.out, "cluster_summary.tsv")
    with open(summary, "w") as fh:
        fh.write("cluster_id\trep_id\tn_members\tn_distinct_asms\tmember_ids\n")
        for cl in clusters:
            asms = set()
            for m in cl["members"]:
                rec = is_id_to_record.get(m)
                if rec: asms.add(rec["assembly"])
            fh.write(f"{cl['id']}\t{cl['rep']}\t{len(cl['members'])}\t"
                     f"{len(asms)}\t{','.join(cl['members'][:10])}"
                     f"{'...' if len(cl['members'])>10 else ''}\n")
    print(f"[cluster] cluster_summary.tsv written.", flush=True)


if __name__ == "__main__":
    main()
