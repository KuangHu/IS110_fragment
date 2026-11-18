#!/usr/bin/env python3
"""Strict re-test for confirmed inversions: do the two flanks of the element
follow their anchors onto OPPOSITE strands of the target?

    REF:    --U--[IS_a]-- M --[IS_b]--D--
    TARGET: --U--[IS_a]-- rev(M) --[IS_b]--D--

Only M flips. Whichever copy our element is, one of its flanks lies in M and
the other outside, so in the target the two flanks are on opposite strands of
one contig -- which is what the Stage 7 candidate says about the two 500-bp
anchors. This script tests it on the flanks themselves:

  1. up / down reference flanks (in transposase orientation, like the anchors),
     --middle-window bp, capped at the estimated inverted span when the
     candidates file gives one (a flank longer than M cannot be 80% inverted);
  2. each is aligned to the target window around ITS anchor hit;
  3. PASS if each flank is covered >= --min-middle-cov (union) at
     >= --min-identity on the SAME strand as its anchor, and the two anchor
     strands differ.

The previous version took the larger reverse-strand block of either flank and
passed on that alone: a target contig that is simply stored reverse-
complemented (no inversion) passed exactly like a real inversion, and anchors
were looked up by ref_id only (the last target's anchors were used for all).
"""
import argparse, csv, gzip, json, os, re, shutil, subprocess, tempfile, sys
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib_rearr as lr  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--validation", required=True)
    p.add_argument("--candidates", required=True,
                   help="input candidates TSV (has 'details' column with anchor coords)")
    p.add_argument("--records",    required=True)
    p.add_argument("--src-dir",    required=True)
    p.add_argument("--out",        required=True)
    p.add_argument("--middle-window", type=int, default=30000,
                   help="size of the ref MIDDLE proxy (30 kb default)")
    p.add_argument("--tgt-context", type=int, default=30000,
                   help="bp of context around target anchors")
    p.add_argument("--min-middle-cov", type=float, default=0.80,
                   help="fraction of ref MIDDLE proxy that must be covered "
                        "by the longest reverse-strand block")
    p.add_argument("--min-identity",   type=float, default=95.0,
                   help="identity (%) of the longest reverse-strand block")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--sample",  type=int, default=0)
    return p.parse_args()


def find_assembly_fa(src_dir, asm, work):
    for ext in (".fna.gz", ".fa.gz", ".fna", ".fa"):
        cand = os.path.join(src_dir, asm + ext)
        if os.path.exists(cand):
            out = os.path.join(work, asm + ".fa")
            if cand.endswith(".gz"):
                with gzip.open(cand, "rt") as fi, open(out, "w") as fo:
                    shutil.copyfileobj(fi, fo)
            else:
                shutil.copy(cand, out)
            return out
    return None



def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)
    loci = lr.load_loci(args.records)

    # Join validation rows (no 'details' col) with the input candidates file,
    # per (ref_id, target assembly)
    cand = {}
    for r in csv.DictReader(open(args.candidates), delimiter="\t"):
        if r.get("details") and r.get("ref_id"):
            cand[(r["ref_id"], r.get("assembly") or r.get("tgt_asm"))] = r

    rows, n_nodetails = [], 0
    for r in csv.DictReader(open(args.validation), delimiter="\t"):
        if r.get("category") != "inversion_only" or r.get("verdict") != "CONFIRMED_BOTH":
            continue
        c = cand.get((r.get("ref_id"), r.get("tgt_asm") or r.get("assembly")))
        if c is None:
            n_nodetails += 1
            continue
        r["details"] = c["details"]
        r["_span"] = int(c.get("estimated_inverted_bp") or 0)
        rows.append(r)
    if args.sample > 0 and len(rows) > args.sample:
        import random
        random.seed(13)
        rows = random.sample(rows, args.sample)
    print(f"Checking {len(rows)} CONFIRMED_BOTH inversion rows "
          f"({n_nodetails} more have no candidate row for that ref+target: SKIPPED)",
          flush=True)

    out_path = os.path.join(args.out, "per_candidate_revcomp.tsv")
    work_root = tempfile.mkdtemp(prefix="revcheck_")
    fa_cache = {}
    def open_asm(asm):
        if asm in fa_cache: return fa_cache[asm]
        path = find_assembly_fa(args.src_dir, asm, work_root)
        fa_cache[asm] = pysam.FastaFile(path) if path else None
        return fa_cache[asm]

    n_pass = n_assessed = 0
    with open(out_path, "w") as out:
        out.write("ref_id\ttgt_asm\tproxy_len\tup_anchor_strand\tdown_anchor_strand\t"
                  "up_bp\tup_frac\tdown_bp\tdown_frac\tpasses\tstatus\n")

        def skip(ref_id, tgt_asm, why):
            out.write(f"{ref_id}\t{tgt_asm}\t0\t\t\t0\t0.000\t0\t0.000\tNA\t{why}\n")

        for i, r in enumerate(rows):
            if i % 50 == 0:
                print(f"  {i}/{len(rows)}", flush=True)
            ref_id = r["ref_id"]
            tgt_asm = r.get("tgt_asm") or r.get("assembly")
            locus = loci.get(ref_id)
            if not locus:
                skip(ref_id, tgt_asm, "ref_id not in records"); continue
            ref_contig, is_s0, is_e, strand = locus
            ref_fa = open_asm(ref_id.split("|")[0])
            tgt_fa = open_asm(tgt_asm)
            if not ref_fa or not tgt_fa:
                skip(ref_id, tgt_asm, "assembly FASTA not found"); continue
            d = lr.parse_details(r["details"])
            pair = next(((u, w) for u in d["up"] for w in d["down"]
                         if u[0] == w[0] and u[3] != w[3]), None)
            if pair is None:
                skip(ref_id, tgt_asm, "no opposite-strand anchor pair on one contig"); continue
            (c, us, ue, ust), (_, ds, de, dst) = pair
            window = args.middle_window
            if r["_span"] > 4000:
                window = min(window, r["_span"] - 3000)
            try:
                up_seq, down_seq = lr.ref_proxies(ref_fa, ref_contig, is_s0, is_e, strand, window)
            except (KeyError, ValueError):
                skip(ref_id, tgt_asm, "ref contig not in FASTA"); continue
            if min(len(up_seq), len(down_seq)) < 1000:
                skip(ref_id, tgt_asm, "reference flank < 1 kb (contig end)"); continue
            t_up, _ = lr.fetch_window(tgt_fa, c, us, ue, args.tgt_context)
            t_dn, _ = lr.fetch_window(tgt_fa, c, ds, de, args.tgt_context)
            if not t_up or not t_dn:
                skip(ref_id, tgt_asm, "target contig not in FASTA"); continue
            work = tempfile.mkdtemp(prefix="rc_", dir=work_root)
            up_bp = lr.strand_cover(lr.align(up_seq, t_up, work, args.threads), ust,
                                    args.min_identity)
            dn_bp = lr.strand_cover(lr.align(down_seq, t_dn, work, args.threads), dst,
                                    args.min_identity)
            shutil.rmtree(work, ignore_errors=True)
            uf, df = up_bp / len(up_seq), dn_bp / len(down_seq)
            passes = uf >= args.min_middle_cov and df >= args.min_middle_cov
            n_assessed += 1
            if passes: n_pass += 1
            out.write(f"{ref_id}\t{tgt_asm}\t{window}\t{ust}\t{dst}\t{up_bp}\t{uf:.3f}\t"
                      f"{dn_bp}\t{df:.3f}\t{passes}\tassessed\n")

    for fa in fa_cache.values():
        if fa: fa.close()
    shutil.rmtree(work_root, ignore_errors=True)
    print(f"\nDONE: {n_pass}/{n_assessed} assessed rows pass (each flank >= "
          f"{args.min_middle_cov} covered at >= {args.min_identity}% on its anchor's "
          f"strand, strands opposite); {len(rows) - n_assessed} not assessable")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
