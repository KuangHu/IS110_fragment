#!/usr/bin/env python3
"""Self-detect tandem composite-transposon duplications: [IS][DNA1][IS][DNA2][IS]
where DNA1 ≈ DNA2 at high identity.

For each assembly:
  1. For each contig with >= 3 same-strand IS hits, examine consecutive triplets.
  2. Each triplet (IS_1, IS_2, IS_3) defines:
       DNA1 = contig[IS_1.end : IS_2.start]
       DNA2 = contig[IS_2.end : IS_3.start]
  3. Filter:
       - both sizes in [--min-dna-bp, --max-dna-bp]  (default 10-200 kb)
       - |DNA1_bp - DNA2_bp| <= --size-tol  (default 5 kb)
  4. Align DNA1 vs DNA2 with minimap2 (forward strand only).
  5. PASS = forward block covers >= --min-cov of min(|DNA1|, |DNA2|)
            AND identity >= --min-identity

Reads sequences directly from the species DB (one indexed FASTA), so no need
to load 22,886 per-assembly FASTAs.
"""
import argparse, csv, os, re, shutil, subprocess, tempfile
from collections import defaultdict
import pysam


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--is-hits",   required=True)
    p.add_argument("--species-db", required=True,
                   help="concatenated indexed FASTA of all assemblies in the species "
                        "(contig names should look like 'ASSEMBLY|CONTIG')")
    p.add_argument("--out",       required=True)
    p.add_argument("--min-dna-bp", type=int, default=10000)
    p.add_argument("--max-dna-bp", type=int, default=200000)
    p.add_argument("--size-tol",   type=int, default=5000)
    p.add_argument("--min-cov",    type=float, default=0.80)
    p.add_argument("--min-identity", type=float, default=95.0)
    p.add_argument("--require-same-strand", action="store_true", default=True)
    p.add_argument("--allow-any-strand", action="store_true",
                   help="if set, do not require all 3 ISes to be on the same strand")
    p.add_argument("--threads", type=int, default=8)
    return p.parse_args()


def load_is_hits(path):
    """assembly → contig (local) → sorted [(s, e, strand)]"""
    idx = defaultdict(lambda: defaultdict(list))
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            asm = r["assembly"]; contig = r["contig"]
            local = contig.split("|", 1)[1] if "|" in contig else contig
            idx[asm][local].append((int(r["tnp_start"]), int(r["tnp_end"]),
                                    r.get("tnp_strand", "+")))
    for asm in idx:
        for c in idx[asm]:
            idx[asm][c].sort()
    return idx


def merge_overlapping(hits, gap=500):
    """Merge IS hits that overlap or are adjacent (HMM domain artifacts)."""
    if not hits: return []
    sh = sorted(hits, key=lambda x: x[0])
    merged = [list(sh[0])]
    for h in sh[1:]:
        if h[0] <= merged[-1][1] + gap:
            merged[-1][1] = max(merged[-1][1], h[1])
        else:
            merged.append(list(h))
    return [tuple(m) for m in merged]


def write_fasta(label, seq, path):
    with open(path, "w") as fh:
        fh.write(f">{label}\n")
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i+80] + "\n")


def align_dna1_vs_dna2(d1_path, d2_path, work, threads):
    """minimap2 d1 vs d2, return best FORWARD block (block_bp, identity)."""
    paf = os.path.join(work, "aln.paf")
    subprocess.run(["minimap2", "-x", "asm10", "-c", "--eqx",
                    "-t", str(threads), d1_path, d2_path, "-o", paf],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    best_block = 0
    best_id = 0.0
    with open(paf) as fh:
        for line in fh:
            c = line.split("\t")
            if len(c) < 12: continue
            if c[4] != "+": continue
            matches, block_len = int(c[9]), int(c[10])
            if block_len > best_block:
                best_block = block_len
                best_id = matches / block_len * 100 if block_len > 0 else 0
    return best_block, best_id


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)
    is_hits = load_is_hits(args.is_hits)
    db = pysam.FastaFile(args.species_db)
    print(f"loaded is_hits: {len(is_hits)} assemblies")
    print(f"opened species DB: {args.species_db}")

    work_root = tempfile.mkdtemp(prefix="tandem_self_")
    out_path = os.path.join(args.out, "tandem_self_dup.tsv")
    n_triplets = 0
    n_pass = 0
    n_size_fail = 0
    n_align_fail = 0
    require_same_strand = args.require_same_strand and not args.allow_any_strand

    with open(out_path, "w") as out:
        out.write("assembly\tcontig\tIS1_start\tIS1_end\tIS2_start\tIS2_end\t"
                  "IS3_start\tIS3_end\tstrand\tDNA1_bp\tDNA2_bp\t"
                  "fwd_block_bp\tidentity\tcov\tpasses\n")
        for ai, (asm, contigs) in enumerate(is_hits.items()):
            if ai % 1000 == 0:
                print(f"  asm {ai}/{len(is_hits)}  triplets={n_triplets}  pass={n_pass}",
                      flush=True)
            for contig_local, raw_hits in contigs.items():
                hits = merge_overlapping(raw_hits)
                if len(hits) < 3:
                    continue
                # Consider every consecutive triplet
                for i in range(len(hits) - 2):
                    is1 = hits[i]; is2 = hits[i+1]; is3 = hits[i+2]
                    if require_same_strand:
                        if not (is1[2] == is2[2] == is3[2]):
                            continue
                    dna1_s, dna1_e = is1[1], is2[0]
                    dna2_s, dna2_e = is2[1], is3[0]
                    dna1_len = dna1_e - dna1_s
                    dna2_len = dna2_e - dna2_s
                    if dna1_len < args.min_dna_bp or dna1_len > args.max_dna_bp:
                        continue
                    if dna2_len < args.min_dna_bp or dna2_len > args.max_dna_bp:
                        continue
                    if abs(dna1_len - dna2_len) > args.size_tol:
                        n_size_fail += 1
                        continue
                    n_triplets += 1
                    # Build the species-DB contig key
                    # ENA-style: "asm|ENA|XXX|XXX.1"; ecoli-style: "asm|XXX.1"
                    db_key_full = f"{asm}|{contig_local}"
                    try:
                        dna1 = db.fetch(db_key_full, dna1_s, dna1_e)
                        dna2 = db.fetch(db_key_full, dna2_s, dna2_e)
                    except (KeyError, ValueError):
                        continue
                    if not dna1 or not dna2:
                        continue
                    work = tempfile.mkdtemp(prefix=f"t{n_triplets}_", dir=work_root)
                    d1p = os.path.join(work, "d1.fa")
                    d2p = os.path.join(work, "d2.fa")
                    write_fasta("DNA1", dna1, d1p)
                    write_fasta("DNA2", dna2, d2p)
                    try:
                        bp, id_ = align_dna1_vs_dna2(d1p, d2p, work, args.threads)
                    except Exception:
                        bp, id_ = 0, 0.0
                    shutil.rmtree(work, ignore_errors=True)
                    cov = bp / min(dna1_len, dna2_len) if min(dna1_len, dna2_len) else 0
                    passes = cov >= args.min_cov and id_ >= args.min_identity
                    if passes:
                        n_pass += 1
                    else:
                        n_align_fail += 1
                    out.write(f"{asm}\t{contig_local}\t{is1[0]}\t{is1[1]}\t"
                              f"{is2[0]}\t{is2[1]}\t{is3[0]}\t{is3[1]}\t"
                              f"{is1[2]}\t{dna1_len}\t{dna2_len}\t"
                              f"{bp}\t{id_:.1f}\t{cov:.3f}\t{passes}\n")
    db.close()
    shutil.rmtree(work_root, ignore_errors=True)
    print(f"\nDONE.")
    print(f"  total triplets evaluated:  {n_triplets}")
    print(f"  size-tol failures:         {n_size_fail}")
    print(f"  alignment failures:        {n_align_fail}")
    print(f"  PASS (tandem duplication): {n_pass}")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
