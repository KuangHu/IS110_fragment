#!/usr/bin/env python3
"""Strict re-test for confirmed inversions using the MIDDLE-segment geometry.

Biological model for an IS-mediated inversion:

    REF:    --upstream--[IS_a]-- MIDDLE --[IS_b]--downstream--
    TARGET: --upstream--[IS_a]-- rev(MIDDLE) --[IS_b]--downstream--

Only the MIDDLE segment flips. The ref IS we track could be either IS_a or IS_b,
so we test BOTH 30 kb flanks (downstream-of-IS as MIDDLE-proxy if our IS is IS_a,
upstream-of-IS as MIDDLE-proxy if our IS is IS_b) and take the side with the
better reverse-strand alignment.

For each CONFIRMED_BOTH inversion in a validation.tsv:
  1. Extract REF DOWN PROXY = ref_contig[IS_end : IS_end + window]
     Extract REF UP PROXY   = ref_contig[IS_start - window : IS_start]
  2. Extract TGT REGION = union of up_anchor + down_anchor contig ranges +/- window
  3. minimap2 -x asm10 -c each ref proxy vs target region
  4. Find largest reverse-strand block in each alignment
  5. Pick the side (down or up) with the larger reverse block
  6. PASS if best block coverage >= --min-middle-cov (default 0.80 of ref proxy)
          AND identity                  >= --min-identity (default 95.0)
"""
import argparse, csv, gzip, json, os, re, shutil, subprocess, tempfile
import pysam


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


def load_loci(records_path):
    loci = {}
    with open(records_path) as f:
        recs = json.load(f)
    for r in recs:
        rid = r.get("ref_id") or r.get("is110_id")
        src = r.get("source", {})
        ie = r.get("is_element") or src.get("is_element") or {}
        contig = src.get("contig")
        s = ie.get("source_start") or src.get("transposase_start")
        e = ie.get("source_end") or src.get("transposase_end")
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


ANCHOR_RE = re.compile(r"([\w.|]+):(\d+)-(\d+)")
def anchor_targets(details):
    return [(m.group(1), int(m.group(2)), int(m.group(3))) for m in ANCHOR_RE.finditer(details)]


def best_rev_block(ref_fa_path, tgt_fa_path, work, threads):
    """Run minimap2 ref vs tgt, return (best_rev_block_bp, best_rev_block_ident)."""
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
            if c[4] != "-": continue
            matches, block = int(c[9]), int(c[10])
            if block > best_block:
                best_block = block
                best_ident = matches / block * 100 if block > 0 else 0
    return best_block, best_ident


def write_target_region(tgt_fa, details, window, out_path):
    """Union of anchor contig ranges +/- window. Returns total bp written."""
    anchors = anchor_targets(details)
    if not anchors:
        return 0
    # group anchors by contig and take min/max
    spans = {}
    for (c, s, e) in anchors:
        if c not in spans:
            spans[c] = [s, e]
        else:
            spans[c][0] = min(spans[c][0], s)
            spans[c][1] = max(spans[c][1], e)
    total = 0
    with open(out_path, "w") as fh:
        for c, (s, e) in spans.items():
            try:
                clen = tgt_fa.get_reference_length(c)
            except (KeyError, ValueError):
                continue
            sx = max(0, s - window)
            ex = min(clen, e + window)
            seq = tgt_fa.fetch(c, sx, ex)
            if seq:
                fh.write(f">TGT_{c}\n")
                for k in range(0, len(seq), 80):
                    fh.write(seq[k:k+80] + "\n")
                total += len(seq)
    return total


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)
    loci = load_loci(args.records)

    # Join validation rows (no 'details' col) with the input candidates file
    details_by_ref = {}
    for r in csv.DictReader(open(args.candidates), delimiter="\t"):
        if r.get("details") and r.get("ref_id"):
            details_by_ref[r["ref_id"]] = r["details"]

    rows = []
    for r in csv.DictReader(open(args.validation), delimiter="\t"):
        if r.get("category") != "inversion_only" or r.get("verdict") != "CONFIRMED_BOTH":
            continue
        r["details"] = details_by_ref.get(r.get("ref_id"), "")
        if r["details"]:
            rows.append(r)
    if args.sample > 0 and len(rows) > args.sample:
        import random
        random.seed(13)
        rows = random.sample(rows, args.sample)
    print(f"Checking {len(rows)} CONFIRMED_BOTH inversion rows ...", flush=True)

    out_path = os.path.join(args.out, "per_candidate_revcomp.tsv")
    work_root = tempfile.mkdtemp(prefix="revcheck_")
    fa_cache = {}
    def open_asm(asm):
        if asm in fa_cache: return fa_cache[asm]
        path = find_assembly_fa(args.src_dir, asm, work_root)
        fa_cache[asm] = pysam.FastaFile(path) if path else None
        return fa_cache[asm]

    n_pass = 0
    with open(out_path, "w") as out:
        out.write("ref_id\ttgt_asm\tbest_side\tref_middle_len\trev_block_bp\t"
                  "rev_block_ident\trev_block_frac_of_middle\tdown_bp\tdown_ident\t"
                  "up_bp\tup_ident\tpasses\n")
        for i, r in enumerate(rows):
            if i % 50 == 0:
                print(f"  {i}/{len(rows)}", flush=True)
            ref_id = r["ref_id"]
            tgt_asm = r.get("tgt_asm") or r.get("assembly")
            locus = loci.get(ref_id)
            if not locus: continue
            ref_contig, is_start, is_end = locus
            ref_asm = ref_id.split("|")[0]
            ref_fa = open_asm(ref_asm)
            tgt_fa = open_asm(tgt_asm)
            if not ref_fa or not tgt_fa: continue
            try:
                clen = ref_fa.get_reference_length(ref_contig)
            except (KeyError, ValueError):
                continue

            work = tempfile.mkdtemp(prefix="rc_", dir=work_root)
            # ref MIDDLE proxies: 30kb downstream of IS, and 30kb upstream of IS
            down_s = is_end
            down_e = min(clen, is_end + args.middle_window)
            up_s = max(0, is_start - args.middle_window)
            up_e = is_start
            down_seq = ref_fa.fetch(ref_contig, down_s, down_e) if down_e > down_s else ""
            up_seq   = ref_fa.fetch(ref_contig, up_s, up_e)     if up_e   > up_s   else ""
            ref_down = os.path.join(work, "ref_down.fa")
            ref_up   = os.path.join(work, "ref_up.fa")
            if down_seq: write_fasta("REF_DOWN", down_seq, ref_down)
            if up_seq:   write_fasta("REF_UP",   up_seq,   ref_up)
            # target region: union of anchor contig spans +/- context
            tgt_region = os.path.join(work, "tgt.fa")
            tgt_bp = write_target_region(tgt_fa, r["details"], args.tgt_context, tgt_region)
            if tgt_bp == 0:
                shutil.rmtree(work, ignore_errors=True); continue

            down_block = down_ident = 0
            up_block = up_ident = 0
            try:
                if down_seq:
                    down_block, down_ident = best_rev_block(ref_down, tgt_region, work, args.threads)
                if up_seq:
                    up_block, up_ident = best_rev_block(ref_up, tgt_region, work, args.threads)
            except Exception:
                shutil.rmtree(work, ignore_errors=True); continue

            # pick the side with the larger rev block as the MIDDLE candidate
            if down_block >= up_block:
                side = "down"; block = down_block; ident = down_ident
                ref_len = len(down_seq)
            else:
                side = "up";   block = up_block;   ident = up_ident
                ref_len = len(up_seq)

            frac = block / ref_len if ref_len else 0
            passes = (frac >= args.min_middle_cov) and (ident >= args.min_identity)
            if passes: n_pass += 1
            out.write(f"{ref_id}\t{tgt_asm}\t{side}\t{ref_len}\t{block}\t"
                      f"{ident:.1f}\t{frac:.3f}\t{down_block}\t{down_ident:.1f}\t"
                      f"{up_block}\t{up_ident:.1f}\t{passes}\n")
            shutil.rmtree(work, ignore_errors=True)

    for fa in fa_cache.values():
        if fa: fa.close()
    shutil.rmtree(work_root, ignore_errors=True)
    print(f"\nDONE: {n_pass}/{len(rows)} pass MIDDLE-segment revcomp check "
          f"(cov >= {args.min_middle_cov} of {args.middle_window} bp at ident >= {args.min_identity}%)")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
