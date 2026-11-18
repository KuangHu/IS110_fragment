#!/usr/bin/env python3
"""Final 3-way validator for IS-mediated inversion candidates (alignment-driven).

Tests the canonical signature:
    REF:  [region_A] [IS_a] [DNA-forward] [IS_b] [region_B]
    TGT:  [region_A] [IS_a] [DNA-reverse-complement] [IS_b] [region_B]

IS strand direction does NOT need to match between ref and tgt — only the
structural pattern matters (IS presence on both sides of the inverted DNA).

For each Tier-2 strict inversion candidate:
  1. Slice ref around candidate IS with a wide window covering the inversion span.
  2. Slice tgt over the union of anchor positions.
  3. Align ref vs tgt to detect orientation; revcomp tgt slice if needed.
  4. Re-align; identify boundaries from the alignment itself:
       - largest forward block at the 5' edge → region_A in both
       - largest reverse block in the middle → DNA (inverted segment)
       - largest forward block at the 3' edge → region_B in both
  5. Confirm an IS is present in tgt within --is-tol bp of EACH alignment
     transition boundary (= the IS_a and IS_b in tgt).
  6. Confirm an IS is present in ref at the candidate position (already known).
  7. Measure coverage + identity of each of A, middle, B.

PASS criteria (configurable):
  - region_A forward block covers >= --min-flank-cov of min(A_ref, A_tgt)
    at >= --min-identity
  - region_B forward block covers >= --min-flank-cov at >= --min-identity
  - DNA reverse block covers >= --min-middle-cov at >= --min-identity
  - At least 2 IS hits found in tgt within --is-tol of the alignment transitions
"""
import argparse, csv, gzip, json, os, re, shutil, subprocess, tempfile, sys
from collections import defaultdict
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib_rearr as lr  # noqa: E402
from lib_alleles import revcomp as _rc  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--candidates", required=True)
    p.add_argument("--anchor-tsv", required=True)
    p.add_argument("--records",    required=True)
    p.add_argument("--is-hits",    required=True)
    p.add_argument("--src-dir",    required=True)
    p.add_argument("--out",        required=True)
    p.add_argument("--flank",      type=int, default=10000)
    p.add_argument("--max-span",   type=int, default=2_000_000)
    p.add_argument("--is-tol",     type=int, default=3000,
                   help="max distance from alignment transition to nearest IS in tgt")
    p.add_argument("--min-flank-cov",  type=float, default=0.80)
    p.add_argument("--min-middle-cov", type=float, default=0.80)
    p.add_argument("--min-identity",   type=float, default=95.0)
    p.add_argument("--min-block-bp",   type=int,   default=500,
                   help="ignore alignment blocks smaller than this when finding boundaries")
    p.add_argument("--threads",  type=int, default=8)
    p.add_argument("--sample",   type=int, default=0)
    return p.parse_args()


ANCHOR_RE = re.compile(r"([\w.|]+):(\d+)-(\d+)\(([+-])\)")


def load_loci(records_path):
    loci = {}
    with open(records_path) as f:
        for r in json.load(f):
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


def load_is_hits(path):
    idx = defaultdict(lambda: defaultdict(list))
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            asm = r["assembly"]; contig = r["contig"]
            local = contig.split("|", 1)[1] if "|" in contig else contig
            hit = (int(r["tnp_start"]), int(r["tnp_end"]),
                   r.get("tnp_strand", "+"))
            idx[asm][local].append(hit)
            tail = local.split("|")[-1]
            if tail != local:
                idx[asm][tail].append(hit)
    return idx


def merge_is(hits, gap=500):
    if not hits: return []
    sh = sorted(hits, key=lambda x: x[0])
    merged = [sh[0]]
    for h in sh[1:]:
        if h[0] <= merged[-1][1] + gap:
            merged[-1] = (min(merged[-1][0], h[0]),
                          max(merged[-1][1], h[1]),
                          merged[-1][2])
        else:
            merged.append(h)
    return merged


def is_within(is_list, pos, tol):
    """True if any IS hit is within `tol` bp of pos."""
    for (s, e, _st) in is_list:
        d = 0 if s <= pos <= e else min(abs(pos - s), abs(pos - e))
        if d <= tol:
            return True
    return False


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


def parse_paf(paf_path):
    blocks = []
    with open(paf_path) as fh:
        for line in fh:
            c = line.split("\t")
            if len(c) < 12: continue
            blocks.append({
                "qs": int(c[2]), "qe": int(c[3]), "strand": c[4],
                "ts": int(c[7]), "te": int(c[8]),
                "matches": int(c[9]), "block": int(c[10]),
                "ident": int(c[9]) / int(c[10]) * 100 if int(c[10]) > 0 else 0,
            })
    return blocks


def align(ref_fa, tgt_fa, work, threads):
    paf = os.path.join(work, "aln.paf")
    subprocess.run(["minimap2", "-x", "asm10", "-c", "--eqx",
                    "-t", str(threads), ref_fa, tgt_fa, "-o", paf],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return parse_paf(paf)


def revcomp(seq):
    tab = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(tab)[::-1]


def parse_anchors(details):
    return [(m.group(1), int(m.group(2)), int(m.group(3)), m.group(4))
            for m in ANCHOR_RE.finditer(details)]


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)
    loci = lr.load_loci(args.records)
    is_hits = load_is_hits(args.is_hits)

    anchor_by_pair = {}
    with open(args.anchor_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            anchor_by_pair[(r['ref_id'], r['assembly'])] = r.get('details', '')

    rows = list(csv.DictReader(open(args.candidates), delimiter="\t"))
    if args.sample > 0 and len(rows) > args.sample:
        import random
        random.seed(13)
        rows = random.sample(rows, args.sample)
    print(f"Checking {len(rows)} Tier-2 strict inversion candidates...", flush=True)

    out_path = os.path.join(args.out, "per_candidate_inversion_full.tsv")
    work_root = tempfile.mkdtemp(prefix="invfull_")
    fa_cache = {}
    def open_asm(asm):
        if asm in fa_cache: return fa_cache[asm]
        path = find_assembly_fa(args.src_dir, asm, work_root)
        fa_cache[asm] = (path, pysam.FastaFile(path) if path else None)
        return fa_cache[asm]

    n_pass = 0
    n_skipped = 0
    n_assessed = 0
    n_too_small = 0
    with open(out_path, "w") as out:
        out.write("ref_id\ttgt_asm\ttgt_flipped\tspan\t"
                  "A_ref_len\tA_tgt_len\tA_fwd_bp\tA_fwd_ident\tA_cov\tA_pass\t"
                  "B_ref_len\tB_tgt_len\tB_fwd_bp\tB_fwd_ident\tB_cov\tB_pass\t"
                  "M_ref_len\tM_tgt_len\tM_rev_bp\tM_rev_ident\tM_cov\tM_pass\t"
                  "tgt_IS_at_both_transitions\tall_pass\n")
        for i, r in enumerate(rows):
            if i % 50 == 0:
                print(f"  {i}/{len(rows)}", flush=True)
            ref_id  = r['ref_id']
            tgt_asm = r['tgt_asm']
            details = anchor_by_pair.get((ref_id, tgt_asm), '')
            if not details:
                n_skipped += 1; continue
            d_ = lr.parse_details(details)
            pair = next(((u, w) for u in d_["up"] for w in d_["down"]
                         if u[0] == w[0] and u[3] != w[3]), None)
            if pair is None:
                n_skipped += 1; continue
            up_hit, down_hit = pair
            tgt_contig = up_hit[0]
            h_left, h_right = sorted([up_hit[1:], down_hit[1:]], key=lambda x: x[0])
            span = h_right[1] - h_left[0]
            if span > args.max_span or span < 1000:
                n_skipped += 1; continue

            locus = loci.get(ref_id)
            if not locus:
                n_skipped += 1; continue
            ref_contig, ref_is_s, ref_is_e, tnp_strand = locus
            ref_asm = ref_id.split("|")[0]
            ref_path, ref_fa = open_asm(ref_asm)
            tgt_path, tgt_fa = open_asm(tgt_asm)
            if not ref_fa or not tgt_fa:
                n_skipped += 1; continue
            try:
                ref_clen = ref_fa.get_reference_length(ref_contig)
                tgt_clen = tgt_fa.get_reference_length(tgt_contig)
            except (KeyError, ValueError):
                n_skipped += 1; continue

            # The test is made in the TARGET's frame, where both ends of the
            # event are known (the two anchor hits):
            #
            #   out_L | anchor_L ........ inner ........ anchor_R | out_R
            #
            #   M  the inner interval must align to the reference on the
            #      OPPOSITE strand to the flanks (it is the inverted segment)
            #   A/B the two outer flanks, taken PAD bp beyond the anchors so
            #      they are clear of the IS copies, must align on the flank
            #      strand
            # Orientation comes from the up anchor (anchors are in transposase
            # orientation), not from which strand has more alignment near the
            # edges of the slice -- that heuristic flipped the target whenever
            # the inverted middle outweighed the 10-kb flanks. Coverage is a
            # union over the region being tested: the old M_cov divided a
            # block by its own extent (always >= 1), so 11.5 kb inverted in a
            # 1.17-Mb span passed.
            is_len = ref_is_e - ref_is_s
            PAD = 1500 + is_len + 1000
            flipped = up_hit[3] != tnp_strand        # target stored reverse to the reference
            want_flank = "-" if flipped else "+"
            want_mid = "+" if flipped else "-"
            t_lo = max(0, h_left[0] - PAD - args.flank)
            t_hi = min(tgt_clen, h_right[1] + PAD + args.flank)
            tgt_seq = tgt_fa.fetch(tgt_contig, t_lo, t_hi).upper()
            ref_lo = max(0, ref_is_s - span - PAD - args.flank)
            ref_hi = min(ref_clen, ref_is_e + span + PAD + args.flank)
            ref_seq = ref_fa.fetch(ref_contig, ref_lo, ref_hi).upper()

            work = tempfile.mkdtemp(prefix="inv_", dir=work_root)
            # a failed alignment raises: it used to be counted as "skipped"
            blocks = [x for x in lr.align(tgt_seq, ref_seq, work, args.threads)
                      if x["ident"] >= args.min_identity
                      and x["qe"] - x["qs"] >= args.min_block_bp]
            shutil.rmtree(work, ignore_errors=True)

            def cover(strand, lo, hi):
                """query bp in [lo, hi) covered by blocks on `strand`."""
                return lr.union_bp([(max(x["qs"], lo), min(x["qe"], hi)) for x in blocks
                                    if x["strand"] == strand and min(x["qe"], hi) > max(x["qs"], lo)])

            aL, aR = h_left[0] - t_lo, h_right[1] - t_lo         # anchor span in the slice
            regA = (max(0, aL - PAD - args.flank), max(0, aL - PAD))
            regB = (min(len(tgt_seq), aR + PAD), min(len(tgt_seq), aR + PAD + args.flank))
            PAD_IN = is_len + 1500     # clear the IS copy next to the unflipped anchor
            regM = (h_left[1] - t_lo + PAD_IN, h_right[0] - t_lo - PAD_IN)
            A_ref_len = A_tgt_len = regA[1] - regA[0]
            B_ref_len = B_tgt_len = regB[1] - regB[0]
            M_ref_len = M_tgt_len = max(0, regM[1] - regM[0])
            if M_ref_len < 1000:
                # inverted interval too short to test once the IS copies are
                # set aside: not assessable (it used to be written as a FAIL)
                n_skipped += 1; n_too_small += 1
                continue
            A_bp, B_bp = cover(want_flank, *regA), cover(want_flank, *regB)
            M_bp = cover(want_mid, *regM) if M_ref_len > 0 else 0
            A_id = B_id = M_id = max((x["ident"] for x in blocks), default=0.0)
            A_cov = A_bp / A_ref_len if A_ref_len >= 1000 else 0
            B_cov = B_bp / B_ref_len if B_ref_len >= 1000 else 0
            M_cov = M_bp / M_ref_len if M_ref_len >= 1000 else 0
            A_pass = A_cov >= args.min_flank_cov
            B_pass = B_cov >= args.min_flank_cov
            M_pass = M_cov >= args.min_middle_cov

            # an IS copy at each end of the inverted interval in the target
            tgt_is = merge_is(is_hits.get(tgt_asm, {}).get(tgt_contig, []))
            left_mid  = (h_left[0]  + h_left[1])  // 2
            right_mid = (h_right[0] + h_right[1]) // 2
            tgt_is_at_both = (is_within(tgt_is, left_mid,  args.is_tol)
                              and is_within(tgt_is, right_mid, args.is_tol))

            all_pass = A_pass and B_pass and M_pass and tgt_is_at_both
            if all_pass: n_pass += 1
            n_assessed += 1

            out.write(f"{ref_id}\t{tgt_asm}\t{flipped}\t{span}\t"
                      f"{A_ref_len}\t{A_tgt_len}\t{A_bp}\t{A_id:.1f}\t{A_cov:.3f}\t{A_pass}\t"
                      f"{B_ref_len}\t{B_tgt_len}\t{B_bp}\t{B_id:.1f}\t{B_cov:.3f}\t{B_pass}\t"
                      f"{M_ref_len}\t{M_tgt_len}\t{M_bp}\t{M_id:.1f}\t{M_cov:.3f}\t{M_pass}\t"
                      f"{tgt_is_at_both}\t{all_pass}\n")
            shutil.rmtree(work, ignore_errors=True)

    for _, fa in fa_cache.values():
        if fa: fa.close()
    shutil.rmtree(work_root, ignore_errors=True)
    print(f"\nDONE: {n_pass}/{n_assessed} assessed candidates pass the FULL 3-way check + "
          f"IS-at-transitions. {n_skipped} of {len(rows)} could not be assessed (no row; "
          f"{n_too_small} of them: inverted interval too short to test).")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
