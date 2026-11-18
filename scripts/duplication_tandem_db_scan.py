#!/usr/bin/env python3
"""DB-driven tandem composite-transposon duplication detector.

For each contig in the species DB (independent of HMM IS detection):
  1. Self-align the contig with minimap2 -X (no self-self trivial alignment).
  2. Find forward-strand alignment blocks that look like tandem repeats:
       - qs < ts                                 (first copy precedes second copy)
       - ts - qe in [0, --adjacency]            (truly adjacent, ~IS-size gap)
       - (qe-qs) in [--min-dna-bp, --max-dna-bp] (size constraint, default 10-200 kb)
       - identity >= --min-identity              (default 95%)
  3. Cross-reference is_hits.tsv: count how many of the three expected IS
     positions (left of block, in gap between copies, right of block) have an
     HMM-detected IS within --is-tol bp.
  4. Verdict:
       - VERIFIED_CANONICAL  : 3/3 ISes detected
       - LIKELY_CANONICAL    : 2/3 (one IS likely degraded / truncated / missed)
       - CANDIDATE           : 1/3
       - REPEAT_NO_IS        : 0/3 (tandem repeat with no detected IS boundary)
"""
import argparse, csv, os, re, shutil, subprocess, tempfile, sys
from collections import defaultdict
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from detect_tandem import find_tandem_unit  # noqa: E402
from lib_alleles import representative_asm  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--species-db", required=True,
                   help="indexed FASTA of the species DB (contigs named ASSEMBLY|CONTIG)")
    p.add_argument("--is-hits",    required=True,
                   help="is_hits.tsv (for cross-reference; detection itself doesn't depend on it)")
    p.add_argument("--out",        required=True)
    p.add_argument("--assembly-chunk-start", type=int, default=0,
                   help="for SLURM array: starting assembly index (0-based)")
    p.add_argument("--assembly-chunk-end",   type=int, default=999999999,
                   help="for SLURM array: end assembly index (exclusive)")
    p.add_argument("--min-contig-bp", type=int, default=30000,
                   help="skip contigs shorter than this (need to fit 2x DNA + 3 IS)")
    p.add_argument("--min-dna-bp", type=int, default=10000)
    p.add_argument("--max-dna-bp", type=int, default=200000)
    p.add_argument("--adjacency",  type=int, default=5000,
                   help="max gap in tgt between qe (end of copy 1) and ts (start of copy 2)")
    p.add_argument("--max-overlap", type=int, default=3000,
                   help="max overlap (= |gap| when ts < qe) allowed. The central IS bridging "
                        "the two DNA copies often gets included in BOTH copy alignments, "
                        "producing a negative gap of ~IS-size. Default 3 kb covers IS110 + slack.")
    p.add_argument("--min-identity", type=float, default=95.0)
    p.add_argument("--is-tol",     type=int, default=3000)
    p.add_argument("--threads",    type=int, default=8)
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


def is_within(hits, pos, tol):
    for (s, e, _st) in hits:
        d = 0 if s <= pos <= e else min(abs(pos - s), abs(pos - e))
        if d <= tol: return True
    return False


def write_fasta(label, seq, path):
    with open(path, "w") as fh:
        fh.write(f">{label}\n")
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i+80] + "\n")


def self_align(contig_fa, work, threads):
    paf = os.path.join(work, "self.paf")
    subprocess.run(["minimap2", "-x", "asm10", "-c", "--eqx", "-X",
                    "-f", "0",  # keep ALL minimizers: the default drops the most frequent
                    # ones, which in a self-alignment are exactly the repeats sought
                    # (fna project: 0 of 60 IS1 copies found without it, 57 with)
                    "-t", str(threads), contig_fa, contig_fa, "-o", paf],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    blocks = []
    with open(paf) as fh:
        for line in fh:
            c = line.split("\t")
            if len(c) < 12: continue
            blocks.append({
                "strand": c[4],
                "qs": int(c[2]), "qe": int(c[3]),
                "ts": int(c[7]), "te": int(c[8]),
                "matches": int(c[9]), "block": int(c[10]),
            })
    return blocks


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)
    is_hits = load_is_hits(args.is_hits)
    db = pysam.FastaFile(args.species_db)
    print(f"loaded is_hits: {len(is_hits)} assemblies")
    print(f"species DB: {len(db.references)} contigs")

    # Group contigs by assembly
    contigs_by_asm = defaultdict(list)
    for name in db.references:
        if "|" not in name: continue
        asm, _ = name.split("|", 1)
        contigs_by_asm[asm].append(name)
    # one of each GCA/GCF pair: the same contig under both accessions gave
    # two rows. NOTE this changes the assembly indices used for chunking.
    asm_list = sorted(set(representative_asm(contigs_by_asm.keys()).values()))
    from lib_alleles import assembly_core
    _by_core = defaultdict(list)
    for a_ in set(contigs_by_asm) | set(is_hits):
        _by_core[assembly_core(a_)].append(a_)
    twins_of = {a_: [b_ for b_ in v if b_ != a_] for v in _by_core.values() for a_ in v}
    start = args.assembly_chunk_start
    end = min(args.assembly_chunk_end, len(asm_list))
    print(f"processing assemblies [{start}:{end}] of {len(asm_list)}")

    work_root = tempfile.mkdtemp(prefix="tandem_db_")
    out_path = os.path.join(args.out, f"tandem_db_{start}_{end}.tsv")
    n_evaluated = 0
    n_repeat_found = 0
    n_verified = 0  # 3/3 IS
    n_likely   = 0  # 2/3
    n_candidate = 0  # 1/3
    n_no_is = 0
    n_align_failed = 0

    with open(out_path, "w") as out:
        out.write("assembly\tcontig\tqs\tqe\tts\tte\tcopy_bp\tgap_bp\t"
                  "identity\tIS_left\tIS_middle\tIS_right\tn_is\tverdict\t"
                  "n_copies\tarray_start\tarray_end\n")
        for ai in range(start, end):
            asm = asm_list[ai]
            if ai % 200 == 0 and ai > start:
                print(f"  asm {ai-start}/{end-start} (abs={ai}/{len(asm_list)})  "
                      f"evaluated_contigs={n_evaluated}  repeats={n_repeat_found}",
                      flush=True)
            for full_contig in contigs_by_asm[asm]:
                try:
                    clen = db.get_reference_length(full_contig)
                except (KeyError, ValueError):
                    continue
                if clen < args.min_contig_bp:
                    continue
                n_evaluated += 1
                local_contig = full_contig.split("|", 1)[1] if "|" in full_contig else full_contig
                seq = db.fetch(full_contig)
                if not seq: continue
                work = tempfile.mkdtemp(prefix=f"a{ai}_", dir=work_root)
                cfp = os.path.join(work, "c.fa")
                write_fasta(local_contig, seq, cfp)
                try:
                    blocks = self_align(cfp, work, args.threads)
                except subprocess.CalledProcessError:
                    # counted, and the run fails at the end: a table from a
                    # run where minimap2 kept failing read as "no tandems"
                    n_align_failed += 1
                    shutil.rmtree(work, ignore_errors=True)
                    continue
                shutil.rmtree(work, ignore_errors=True)
                tgt_hits = is_hits.get(asm, {}).get(local_contig, [])
                if not tgt_hits:
                    # is_hits.tsv keeps ONE of a GCA/GCF pair, not necessarily
                    # the one selected here: look the contig up under the twin
                    bare = local_contig[3:] if local_contig.startswith("NZ_") else local_contig
                    for twin in twins_of.get(asm, ()):
                        th = is_hits.get(twin, {})
                        tgt_hits = (th.get(local_contig) or th.get(bare)
                                    or th.get("NZ_" + bare) or [])
                        if tgt_hits:
                            break
                # Arrays by OFFSET (detect_tandem.find_tandem_unit): unit =
                # offset, copies = array span / unit. The old filter looked at
                # one self-hit at a time and required a small gap between its
                # two halves: for 3 copies the hit overlaps itself by a whole
                # unit (dropped), for 4 the surviving hit covered two units
                # (reported as one copy twice as long), and two copies with
                # 4 kb of unrelated sequence between them were "canonical".
                hits = [(b["qs"], b["qe"], b["ts"], b["te"], b["strand"],
                         b["matches"] / b["block"] * 100 if b["block"] else 0.0) for b in blocks]
                for arr in find_tandem_unit(hits, args.min_dna_bp, args.min_identity):
                    if arr["repeat_class"] != "tandem":
                        continue
                    unit = arr["repeat_unit_bp"]
                    if unit < args.min_dna_bp or unit > args.max_dna_bp + 5000:
                        continue
                    n_repeat_found += 1
                    a0, a1, n = arr["span_start"], arr["span_end"], arr["n_copies"]
                    is_left = is_within(tgt_hits, a0, args.is_tol)
                    is_right = is_within(tgt_hits, a1, args.is_tol)
                    # an IS at EVERY internal copy boundary
                    is_middle = all(is_within(tgt_hits, a0 + k * unit, args.is_tol)
                                    for k in range(1, n))
                    n_is = int(is_left) + int(is_middle) + int(is_right)
                    if   n_is == 3: verdict = "VERIFIED_CANONICAL"; n_verified += 1
                    elif n_is == 2: verdict = "LIKELY_CANONICAL"; n_likely += 1
                    elif n_is == 1: verdict = "CANDIDATE"; n_candidate += 1
                    else:           verdict = "REPEAT_NO_IS"; n_no_is += 1
                    out.write(f"{asm}\t{local_contig}\t{a0}\t{a0 + unit}\t"
                              f"{a0 + unit}\t{min(a1, a0 + 2 * unit)}\t{unit}\t0\t"
                              f"{arr.get('identity', float('nan')):.1f}\t"
                              f"{is_left}\t{is_middle}\t{is_right}\t{n_is}\t{verdict}\t"
                              f"{n}\t{a0}\t{a1}\n")
    db.close()
    shutil.rmtree(work_root, ignore_errors=True)
    print(f"\nDONE.")
    print(f"  contigs evaluated:        {n_evaluated}")
    print(f"  tandem repeat candidates: {n_repeat_found}")
    print(f"  VERIFIED_CANONICAL (3/3): {n_verified}")
    print(f"  LIKELY_CANONICAL  (2/3):  {n_likely}")
    print(f"  CANDIDATE         (1/3):  {n_candidate}")
    print(f"  REPEAT_NO_IS      (0/3):  {n_no_is}")
    print(f"Wrote {out_path}")
    if n_align_failed:
        sys.exit(f"ERROR: self-alignment failed for {n_align_failed} contigs; "
                 f"the table is incomplete")


if __name__ == "__main__":
    main()
