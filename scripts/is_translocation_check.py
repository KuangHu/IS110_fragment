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
import argparse, csv, gzip, json, os, re, shutil, subprocess, tempfile
from collections import defaultdict
import pysam


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


def load_loci(records_path):
    loci = {}
    with open(records_path) as f:
        recs = json.load(f)
    for r in recs:
        rid = r.get("ref_id") or r.get("is110_id")
        src = r.get("source", {})
        ie  = r.get("is_element") or src.get("is_element") or {}
        contig = src.get("contig")
        s = ie.get("source_start") or src.get("transposase_start")
        e = ie.get("source_end")   or src.get("transposase_end")
        if rid and contig and s and e:
            local = contig.split("|", 1)[1] if "|" in contig else contig
            loci[rid] = (local, int(s), int(e))
    return loci


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


def write_fasta(label, seq, path):
    with open(path, "w") as fh:
        fh.write(f">{label}\n")
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i+80] + "\n")


ANCHOR_RE = re.compile(r"([\w.|]+):(\d+)-(\d+)\(([+-])\)")
def parse_anchors_block(details, side_tag):
    """Parse all anchor hits in the side_tag (=up or down) block of details.
    Returns list of (contig, s, e, strand)."""
    # details ~ "up=[...] down=[...]"
    m = re.search(rf"{side_tag}=\[([^\]]*)\]", details)
    if not m:
        return []
    return [(mm.group(1), int(mm.group(2)), int(mm.group(3)), mm.group(4))
            for mm in ANCHOR_RE.finditer(m.group(1))]


def write_target_region(tgt_fa, contig, s_min, e_max, window, out_path, label):
    try:
        clen = tgt_fa.get_reference_length(contig)
    except (KeyError, ValueError):
        return 0
    sx = max(0, s_min - window)
    ex = min(clen, e_max + window)
    seq = tgt_fa.fetch(contig, sx, ex)
    if not seq:
        return 0
    with open(out_path, "w") as fh:
        fh.write(f">{label}_{contig}_{sx}_{ex}\n")
        for k in range(0, len(seq), 80):
            fh.write(seq[k:k+80] + "\n")
    return len(seq)


def best_forward_block(ref_fa_path, tgt_fa_path, work, threads):
    """Return (best_fwd_block_bp, best_fwd_block_ident) — biggest + strand block."""
    paf = os.path.join(work, "aln.paf")
    subprocess.run(["minimap2", "-x", "asm10", "-c", "--eqx",
                    "-t", str(threads), ref_fa_path, tgt_fa_path, "-o", paf],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    best_block = 0
    best_ident = 0.0
    with open(paf) as f:
        for line in f:
            c = line.split("\t")
            if len(c) < 12: continue
            if c[4] != "+": continue
            matches, block_len = int(c[9]), int(c[10])
            if block_len > best_block:
                best_block = block_len
                best_ident = matches / block_len * 100 if block_len > 0 else 0
    return best_block, best_ident


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)
    loci = load_loci(args.records)

    rows = []
    for r in csv.DictReader(open(args.candidates), delimiter="\t"):
        # accept either pre-filtered CONFIRMED_BOTH catalogue or raw with verdict col
        verdict = r.get("verdict")
        if verdict and verdict != "CONFIRMED_BOTH":
            continue
        details = r.get("details", "")
        if not details: continue
        ups   = parse_anchors_block(details, "up")
        downs = parse_anchors_block(details, "down")
        if not ups or not downs: continue
        # require they be on DIFFERENT contigs (translocation signature)
        up_contigs = {a[0] for a in ups}
        down_contigs = {a[0] for a in downs}
        if up_contigs == down_contigs: continue
        r["_up"] = ups[0]      # pick the first up anchor
        r["_down"] = downs[0]  # pick the first down anchor
        rows.append(r)
    if args.sample > 0 and len(rows) > args.sample:
        import random
        random.seed(13)
        rows = random.sample(rows, args.sample)
    print(f"Checking {len(rows)} translocation candidates ...", flush=True)

    out_path = os.path.join(args.out, "per_candidate_translocation.tsv")
    work_root = tempfile.mkdtemp(prefix="trcheck_")
    fa_cache = {}
    def open_asm(asm):
        if asm in fa_cache: return fa_cache[asm]
        path = find_assembly_fa(args.src_dir, asm, work_root)
        fa_cache[asm] = pysam.FastaFile(path) if path else None
        return fa_cache[asm]

    min_block_bp = int(args.middle_window * args.min_block_cov)
    n_pass = 0
    with open(out_path, "w") as out:
        out.write("ref_id\ttgt_asm\tA_proxy_len\tB_proxy_len\t"
                  "A_fwd_block_bp\tA_fwd_block_ident\tA_frac_of_proxy\t"
                  "B_fwd_block_bp\tB_fwd_block_ident\tB_frac_of_proxy\t"
                  "A_passes\tB_passes\tpasses\n")
        for i, r in enumerate(rows):
            if i % 50 == 0:
                print(f"  {i}/{len(rows)}", flush=True)
            ref_id  = r.get("ref_id")
            tgt_asm = r.get("tgt_asm") or r.get("assembly")
            locus = loci.get(ref_id)
            if not locus: continue
            ref_contig, is_s, is_e = locus
            ref_asm = ref_id.split("|")[0]
            ref_fa = open_asm(ref_asm)
            tgt_fa = open_asm(tgt_asm)
            if not ref_fa or not tgt_fa: continue
            try:
                clen = ref_fa.get_reference_length(ref_contig)
            except (KeyError, ValueError):
                continue

            work = tempfile.mkdtemp(prefix="tr_", dir=work_root)
            # ref A-proxy (upstream of IS) and B-proxy (downstream of IS)
            A_s = max(0, is_s - args.middle_window); A_e = is_s
            B_s = is_e; B_e = min(clen, is_e + args.middle_window)
            A_seq = ref_fa.fetch(ref_contig, A_s, A_e) if A_e > A_s else ""
            B_seq = ref_fa.fetch(ref_contig, B_s, B_e) if B_e > B_s else ""
            A_ref_path = os.path.join(work, "ref_A.fa")
            B_ref_path = os.path.join(work, "ref_B.fa")
            if A_seq: write_fasta("REF_A", A_seq, A_ref_path)
            if B_seq: write_fasta("REF_B", B_seq, B_ref_path)

            # target regions around the two anchors
            up_c, up_s, up_e, _ = r["_up"]
            down_c, down_s, down_e, _ = r["_down"]
            tgt1_path = os.path.join(work, "tgt1.fa")
            tgt2_path = os.path.join(work, "tgt2.fa")
            tgt1_bp = write_target_region(tgt_fa, up_c, up_s, up_e,
                                          args.tgt_context, tgt1_path, "TGT1")
            tgt2_bp = write_target_region(tgt_fa, down_c, down_s, down_e,
                                          args.tgt_context, tgt2_path, "TGT2")

            A_block = A_ident = 0; B_block = B_ident = 0
            try:
                if A_seq and tgt1_bp > 0:
                    A_block, A_ident = best_forward_block(A_ref_path, tgt1_path,
                                                          work, args.threads)
                if B_seq and tgt2_bp > 0:
                    B_block, B_ident = best_forward_block(B_ref_path, tgt2_path,
                                                          work, args.threads)
            except Exception:
                pass

            A_len = len(A_seq); B_len = len(B_seq)
            A_frac = A_block / A_len if A_len else 0
            B_frac = B_block / B_len if B_len else 0
            A_pass = (A_block >= min_block_bp and A_ident >= args.min_identity)
            B_pass = (B_block >= min_block_bp and B_ident >= args.min_identity)
            passes = A_pass and B_pass
            if passes: n_pass += 1
            out.write(f"{ref_id}\t{tgt_asm}\t{A_len}\t{B_len}\t"
                      f"{A_block}\t{A_ident:.1f}\t{A_frac:.3f}\t"
                      f"{B_block}\t{B_ident:.1f}\t{B_frac:.3f}\t"
                      f"{A_pass}\t{B_pass}\t{passes}\n")
            shutil.rmtree(work, ignore_errors=True)

    for fa in fa_cache.values():
        if fa: fa.close()
    shutil.rmtree(work_root, ignore_errors=True)
    print(f"\nDONE: {n_pass}/{len(rows)} pass BILATERAL forward-block check "
          f"(both A and B proxies >= {args.min_block_cov} cov at >= {args.min_identity}% ident)")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
