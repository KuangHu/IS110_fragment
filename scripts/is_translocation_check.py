#!/usr/bin/env python3
"""Strict IS-mediated translocation check (bilateral forward-block test).

For a candidate (ref_id, tgt_asm) flagged as `translocation_only` by
find_rearrangements.py:

    REF:    --upstream(A)-- [IS_ref] --downstream(B)--
    TGT (translocation):
      contig 1: --upstream(A)-- [IS_tgt_1] --[C, different]--
      contig 2: --[X, different]-- [IS_tgt_2] --downstream(B)--

The 1 kb anchors confirm:
  - up_anchor (last 1 kb of A) hits tgt contig 1 near IS_tgt_1
  - down_anchor (first 1 kb of B) hits tgt contig 2 near IS_tgt_2

To distinguish IS-mediated HR translocation (where the segments A and B are
truly conserved between ref and tgt, just shuffled) from simple IS transposition
(where the two tgt IS copies are independent events with no segment conservation),
we run a BILATERAL FORWARD-BLOCK TEST:

  1. Extract ref A-proxy = ref[IS_start - middle_window : IS_start]
     (30 kb of conserved upstream sequence in ref)
  2. Extract ref B-proxy = ref[IS_end : IS_end + middle_window]
     (30 kb of conserved downstream sequence in ref)
  3. Extract tgt contig 1 region around up_anchor +/- tgt_context
  4. Extract tgt contig 2 region around down_anchor +/- tgt_context
  5. minimap2 -x asm10 -c ref_A_proxy vs tgt_contig1_region
       → find best FORWARD-strand block
  6. minimap2 -x asm10 -c ref_B_proxy vs tgt_contig2_region
       → find best FORWARD-strand block
  7. PASS if BOTH forward blocks cover >= --min-block-cov of the 30 kb proxy
     at >= --min-identity (default 0.80 cov, 95.0 ident)

Pass = bilateral conservation = HR-mediated translocation (Mech B).
Fail = anchor-only / asymmetric conservation = could be simple transposition
       or partial homology rather than true segment exchange.
"""
import argparse, csv, gzip, json, os, re, shutil, subprocess, tempfile, sys
from collections import defaultdict
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib_rearr as lr  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--candidates", required=True,
                   help="strict_translocation_catalogue.tsv (or per-species pre-filtered TSV)")
    p.add_argument("--records",    required=True)
    p.add_argument("--src-dir",    required=True)
    p.add_argument("--out",        required=True)
    p.add_argument("--middle-window", type=int, default=30000,
                   help="size of the ref A-proxy and B-proxy (30 kb default)")
    p.add_argument("--tgt-context",   type=int, default=30000,
                   help="bp of context around each target anchor")
    p.add_argument("--min-block-cov", type=float, default=0.80,
                   help="each forward block must cover >= this fraction of ref proxy")
    p.add_argument("--min-identity",  type=float, default=95.0)
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

    rows = []
    n_in = n_same_contig = 0
    for r in csv.DictReader(open(args.candidates), delimiter="\t"):
        if r.get("ref_id") == "ref_id":
            continue                      # repeated header line in a concatenated input
        # accept either pre-filtered CONFIRMED_BOTH catalogue or raw with verdict col
        verdict = r.get("verdict")
        if verdict and verdict != "CONFIRMED_BOTH":
            continue
        n_in += 1
        d = lr.parse_details(r.get("details", ""))
        # one up and one down anchor on DIFFERENT contigs (the translocation
        # signature). ups[0]/downs[0] could be on the same contig.
        pair = next(((u, w) for u in d["up"] for w in d["down"] if u[0] != w[0]), None)
        if pair is None:
            n_same_contig += 1
            continue
        r["_up"], r["_down"] = pair
        rows.append(r)
    if args.sample > 0 and len(rows) > args.sample:
        import random
        random.seed(13)
        rows = random.sample(rows, args.sample)
    print(f"Checking {len(rows)} translocation candidates "
          f"({n_same_contig} of {n_in} skipped: no up/down pair on different contigs)",
          flush=True)

    out_path = os.path.join(args.out, "per_candidate_translocation.tsv")
    work_root = tempfile.mkdtemp(prefix="trcheck_")
    fa_cache = {}
    def open_asm(asm):
        if asm in fa_cache: return fa_cache[asm]
        path = find_assembly_fa(args.src_dir, asm, work_root)
        fa_cache[asm] = pysam.FastaFile(path) if path else None
        return fa_cache[asm]

    n_pass = n_assessed = 0
    with open(out_path, "w") as out:
        out.write("ref_id\ttgt_asm\tA_proxy_len\tB_proxy_len\t"
                  "A_fwd_block_bp\tA_fwd_block_ident\tA_frac_of_proxy\t"
                  "B_fwd_block_bp\tB_fwd_block_ident\tB_frac_of_proxy\t"
                  "A_passes\tB_passes\tpasses\tstatus\n")

        def skip(ref_id, tgt_asm, why):
            # every candidate gets a row; these used to disappear
            out.write(f"{ref_id}\t{tgt_asm}\t0\t0\t0\t0.0\t0.000\t0\t0.0\t0.000\t"
                      f"NA\tNA\tNA\t{why}\n")

        for i, r in enumerate(rows):
            if i % 50 == 0:
                print(f"  {i}/{len(rows)}", flush=True)
            ref_id  = r.get("ref_id")
            tgt_asm = r.get("tgt_asm") or r.get("assembly")
            locus = loci.get(ref_id)
            if not locus:
                skip(ref_id, tgt_asm, "ref_id not in records"); continue
            ref_contig, is_s0, is_e, strand = locus
            ref_fa = open_asm(ref_id.split("|")[0])
            tgt_fa = open_asm(tgt_asm)
            if not ref_fa or not tgt_fa:
                skip(ref_id, tgt_asm, "assembly FASTA not found"); continue
            try:
                # A = UPSTREAM, B = DOWNSTREAM of the element in transposase
                # orientation -- the frame the up/down anchors are in. They were
                # "genomic left/right", so for a minus-strand element A was
                # aligned to the contig holding the DOWN anchor.
                A_seq, B_seq = lr.ref_proxies(ref_fa, ref_contig, is_s0, is_e, strand,
                                              args.middle_window)
            except (KeyError, ValueError):
                skip(ref_id, tgt_asm, "ref contig not in FASTA"); continue
            (up_c, up_s, up_e, up_st), (dn_c, dn_s, dn_e, dn_st) = r["_up"], r["_down"]
            t1, _ = lr.fetch_window(tgt_fa, up_c, up_s, up_e, args.tgt_context)
            t2, _ = lr.fetch_window(tgt_fa, dn_c, dn_s, dn_e, args.tgt_context)
            if not t1 or not t2:
                skip(ref_id, tgt_asm, "target contig not in FASTA"); continue
            if min(len(A_seq), len(B_seq)) < 1000:
                skip(ref_id, tgt_asm, "reference flank < 1 kb (contig end)"); continue
            # each target contig must continue past the element on the
            # element's side of its anchor; otherwise this is a contig break
            # inside or next to the element, not a translocation
            need = 1000 + (is_e - is_s0) + 5000
            room_u = (tgt_fa.get_reference_length(up_c) - up_e) if up_st == "+" else up_s
            room_d = dn_s if dn_st == "+" else (tgt_fa.get_reference_length(dn_c) - dn_e)
            if min(room_u, room_d) < need:
                skip(ref_id, tgt_asm, "contig ends within the element region (assembly break)")
                continue

            work = tempfile.mkdtemp(prefix="tr_", dir=work_root)
            # the flank must lie on the SAME target strand as its anchor hit;
            # coverage is a union over the proxy (was: the single largest '+'
            # block against a fixed 24 kb)
            ba = [x for x in lr.align(A_seq, t1, work, args.threads) if x["strand"] == up_st]
            bb = [x for x in lr.align(B_seq, t2, work, args.threads) if x["strand"] == dn_st]
            shutil.rmtree(work, ignore_errors=True)
            A_block = lr.strand_cover(ba, up_st, args.min_identity)
            B_block = lr.strand_cover(bb, dn_st, args.min_identity)
            A_ident = max((x["ident"] for x in ba), default=0.0)
            B_ident = max((x["ident"] for x in bb), default=0.0)
            A_len, B_len = len(A_seq), len(B_seq)
            A_frac, B_frac = A_block / A_len, B_block / B_len
            A_pass = A_frac >= args.min_block_cov
            B_pass = B_frac >= args.min_block_cov
            passes = A_pass and B_pass
            n_assessed += 1
            if passes: n_pass += 1
            out.write(f"{ref_id}\t{tgt_asm}\t{A_len}\t{B_len}\t"
                      f"{A_block}\t{A_ident:.1f}\t{A_frac:.3f}\t"
                      f"{B_block}\t{B_ident:.1f}\t{B_frac:.3f}\t"
                      f"{A_pass}\t{B_pass}\t{passes}\tassessed\n")

    for fa in fa_cache.values():
        if fa: fa.close()
    shutil.rmtree(work_root, ignore_errors=True)
    print(f"\nDONE: {n_pass}/{n_assessed} assessed candidates pass the BILATERAL check "
          f"(both flanks >= {args.min_block_cov} of the proxy at >= {args.min_identity}% ident, "
          f"on the anchor's strand); {len(rows) - n_assessed} not assessable")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
