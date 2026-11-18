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
import argparse, csv, os, re, shutil, subprocess, tempfile, sys
from collections import defaultdict
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_alleles import assembly_core, _asm_rank  # noqa: E402


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
    n_nofetch = 0
    arrays = []
    require_same_strand = args.require_same_strand and not args.allow_any_strand

    with open(out_path, "w") as out:
        out.write("assembly\tcontig\tIS1_start\tIS1_end\tIS2_start\tIS2_end\t"
                  "IS3_start\tIS3_end\tstrand\tDNA1_bp\tDNA2_bp\t"
                  "fwd_block_bp\tidentity\tcov\tpasses\n")
        # one of each GCA/GCF pair
        best_twin = {}
        for asm in is_hits:
            c = assembly_core(asm)
            if c not in best_twin or _asm_rank(asm) > _asm_rank(best_twin[c]):
                best_twin[c] = asm
        keep_asm = set(best_twin.values())
        for ai, (asm, contigs) in enumerate(sorted(is_hits.items())):
            if asm not in keep_asm:
                continue
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
                    # tnp coordinates are 1-based inclusive: the bases strictly
                    # between two IS are [end, next_start - 1) 0-based
                    dna1_s, dna1_e = is1[1], is2[0] - 1
                    dna2_s, dna2_e = is2[1], is3[0] - 1
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
                        n_nofetch += 1
                        continue
                    if not dna1 or not dna2:
                        n_nofetch += 1
                        continue
                    work = tempfile.mkdtemp(prefix=f"t{n_triplets}_", dir=work_root)
                    d1p = os.path.join(work, "d1.fa")
                    d2p = os.path.join(work, "d2.fa")
                    write_fasta("DNA1", dna1, d1p)
                    write_fasta("DNA2", dna2, d2p)
                    # a minimap2 failure raises (it was scored as "no alignment")
                    bp, id_ = align_dna1_vs_dna2(d1p, d2p, work, args.threads)
                    shutil.rmtree(work, ignore_errors=True)
                    cov = bp / min(dna1_len, dna2_len) if min(dna1_len, dna2_len) else 0
                    passes = cov >= args.min_cov and id_ >= args.min_identity
                    if passes:
                        n_pass += 1
                        # consecutive passing triplets on one contig are ONE
                        # array: IS-S-IS-S-IS-S-IS gave two "duplications"
                        k = (asm, contig_local)
                        if arrays and arrays[-1]["key"] == k and arrays[-1]["last_i"] == i - 1:
                            arrays[-1]["last_i"] = i
                            arrays[-1]["n_copies"] += 1
                            arrays[-1]["end"] = is3[1]
                        else:
                            arrays.append({"key": k, "last_i": i, "n_copies": 2,
                                           "start": is1[0], "end": is3[1],
                                           "unit_bp": dna1_len})
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
    print(f"  not similar enough:        (counted above as 'alignment failures')")
    print(f"  sequence not in DB:        {n_nofetch}")
    print(f"  PASS triplets:             {n_pass}")
    arr_path = os.path.join(args.out, "tandem_self_arrays.tsv")
    with open(arr_path, "w") as fh:
        fh.write("assembly\tcontig\tarray_start\tarray_end\tn_copies\tunit_dna_bp\n")
        for a in arrays:
            fh.write(f"{a['key'][0]}\t{a['key'][1]}\t{a['start']}\t{a['end']}\t"
                     f"{a['n_copies']}\t{a['unit_bp']}\n")
    print(f"  tandem duplication EVENTS: {len(arrays)}  (one per array; a k-copy array "
          f"is k-1 passing triplets)  -> {arr_path}")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
