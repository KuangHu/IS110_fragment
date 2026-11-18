#!/usr/bin/env python3
"""Classify MOVED composite-transposon hits as CUT_PASTE or COPY_PASTE.

Given v3 movement output (MOVED rows), for each (ref, tgt) hit:
  1. Extract ref's A_flank (5 kb upstream of IS_a) and B_flank (5 kb downstream of IS_b)
  2. minimap2 both flanks against tgt → find best forward block of each
  3. If both flanks land on the SAME tgt contig:
        dist_in_tgt = | start of B's landing - end of A's landing |
        expected_with_DNA = |IS_a| + |DNA| + |IS_b|
     a) dist_in_tgt ≈ expected_with_DNA  → COPY_PASTE (original site still has element)
     b) dist_in_tgt < 5 kb                → CUT_PASTE (flanks joined; element excised)
     c) anything else                      → AMBIGUOUS (other rearrangement at original site)
  4. If flanks don't both land on same tgt contig → ORIGIN_DIVERGED (unrelated genomic context)
"""
import argparse, csv, gzip, json, os, re, shutil, subprocess, tempfile
from collections import defaultdict
import pysam


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--movement-tsv", required=True,
                   help="v3 per_candidate_movement.tsv output")
    p.add_argument("--is-hits",      required=True)
    p.add_argument("--src-dir",      required=True)
    p.add_argument("--ref-asm",      required=True)
    p.add_argument("--out",          required=True)
    p.add_argument("--flank",        type=int, default=5000)
    p.add_argument("--threads",      type=int, default=8)
    p.add_argument("--copy-paste-margin", type=int, default=2000,
                   help="dist_in_tgt within ±this of expected → COPY_PASTE")
    p.add_argument("--cut-paste-max", type=int, default=5000,
                   help="dist_in_tgt < this → CUT_PASTE (flanks joined directly)")
    return p.parse_args()


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


def find_assembly_fa(src_dir, asm, work):
    for ext in (".fna.gz", ".fa.gz", ".fna", ".fa"):
        path = os.path.join(src_dir, asm + ext)
        if os.path.exists(path):
            out = os.path.join(work, asm + ".fa")
            if path.endswith(".gz"):
                with gzip.open(path, "rt") as fi, open(out, "w") as fo:
                    shutil.copyfileobj(fi, fo)
            else:
                shutil.copy(path, out)
            return out
    return None


def write_fasta(label, seq, path):
    with open(path, "w") as fh:
        fh.write(f">{label}\n")
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i+80] + "\n")


def flank_landing(tgt_path, query_path, work, threads, prox_end, min_ident=90.0,
                  min_frac=0.5):
    """Where does a reference flank land in the target assembly?

    Returns (contig, strand, pos) or None. `pos` is where the flank's
    ELEMENT-PROXIMAL end (prox_end = 'end' for the left flank A, 'start' for
    the right flank B) projects onto the target. Only a block at >= min_ident
    covering >= min_frac of the flank counts, on either strand. (It used to be
    the longest '+' block of any quality: a partial or repeat hit gave a
    meaningless distance, and a target contig stored in the other orientation
    was 'ORIGIN_DIVERGED'.)"""
    paf = os.path.join(work, "_aln.paf")
    subprocess.run(["minimap2", "-x", "asm10", "-c", "--eqx",
                    "-t", str(threads), tgt_path, query_path, "-o", paf],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    best = None
    with open(paf) as fh:
        for line in fh:
            c = line.split("\t")
            if len(c) < 12: continue
            qlen, qs, qe = int(c[1]), int(c[2]), int(c[3])
            ident = int(c[9]) / int(c[10]) * 100 if int(c[10]) else 0
            if ident < min_ident or (qe - qs) < min_frac * qlen:
                continue
            ts, te, strand = int(c[7]), int(c[8]), c[4]
            if prox_end == "end":
                pos = te + (qlen - qe) if strand == "+" else ts - (qlen - qe)
            else:
                pos = ts - qs if strand == "+" else te + qs
            if best is None or (qe - qs) > best[3]:
                best = (c[5], strand, pos, qe - qs)
    return best[:3] if best else None


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)
    is_hits = load_is_hits(args.is_hits)

    moved_rows = [r for r in csv.DictReader(open(args.movement_tsv), delimiter="\t")
                  if r['verdict'] == 'MOVED']
    print(f"Classifying {len(moved_rows)} MOVED hits...", flush=True)

    work_root = tempfile.mkdtemp(prefix="mech_")
    ref_fa_path = find_assembly_fa(args.src_dir, args.ref_asm, work_root)
    if not ref_fa_path:
        raise SystemExit(f"ref FASTA not found")
    ref_fa = pysam.FastaFile(ref_fa_path)

    # Pre-extract ref A_flank and B_flank per unique ref candidate
    flank_cache = {}
    for r in moved_rows:
        key = (r['ref_contig'], int(r['ref_is_a_start']), int(r['ref_is_b_start']))
        if key in flank_cache: continue
        contig, is_a_s, is_b_s = key
        merged = merge_is(is_hits.get(args.ref_asm, {}).get(contig, []))
        is_a = next((h for h in merged if abs(h[0] - is_a_s) < 200), None)
        is_b = next((h for h in merged if abs(h[0] - is_b_s) < 200), None)
        if is_a is None or is_b is None:
            flank_cache[key] = None; continue
        try:
            clen = ref_fa.get_reference_length(contig)
        except (KeyError, ValueError):
            flank_cache[key] = None; continue
        lf_s = max(0, is_a[0] - 1 - args.flank); lf_e = is_a[0] - 1   # 1-based start
        rf_s = is_b[1]; rf_e = min(clen, is_b[1] + args.flank)
        A = ref_fa.fetch(contig, lf_s, lf_e) if lf_e > lf_s else ""
        B = ref_fa.fetch(contig, rf_s, rf_e) if rf_e > rf_s else ""
        expected_with_DNA = (is_b[0] - is_a[1])  # = |IS_a|+|DNA|+|IS_b|, but is_a[1] is IS_a end, is_b[0] is IS_b start
        # actually we want full IS_a-DNA-IS_b length = is_b[1] - is_a[0]
        expected_with_DNA = is_b[1] - is_a[0] + 1
        flank_cache[key] = {"A_seq": A, "B_seq": B,
                            "is_a": is_a, "is_b": is_b,
                            "expected_with_DNA": expected_with_DNA}

    out_path = os.path.join(args.out, "mechanism_classified.tsv")
    n = {"CUT_PASTE":0, "COPY_PASTE":0, "AMBIGUOUS":0, "ORIGIN_DIVERGED":0, "FAILED":0}
    with open(out_path, "w") as out:
        # Write all original columns + extra mechanism columns
        out.write("ref_asm\tref_contig\tref_is_a_start\tref_is_b_start\tsegment_bp\t"
                  "tgt_asm\tB_cov\tB_ident\ttgt_contig\ttgt_pos\t"
                  "Aflank_tgt_pos\tBflank_tgt_pos\tdist_in_tgt\texpected_with_DNA\t"
                  "mechanism\n")
        for i, r in enumerate(moved_rows):
            if i % 10 == 0:
                print(f"  {i}/{len(moved_rows)}", flush=True)
            key = (r['ref_contig'], int(r['ref_is_a_start']), int(r['ref_is_b_start']))
            fc = flank_cache.get(key)
            Apos = Bpos = "NA"; dist = "NA"
            tgt_path = find_assembly_fa(args.src_dir, r['tgt_asm'], work_root) if fc else None
            # rows that cannot be classified are WRITTEN (they used to be
            # counted and dropped from the table)
            if fc is None or min(len(fc["A_seq"]), len(fc["B_seq"])) < 1000:
                mechanism = "FAILED_NO_REF_FLANK"
            elif not tgt_path:
                mechanism = "FAILED_NO_TARGET_FASTA"
            else:
                work = tempfile.mkdtemp(prefix=f"m{i}_", dir=work_root)
                A_path = os.path.join(work, "A.fa")
                B_path = os.path.join(work, "B.fa")
                write_fasta("A", fc["A_seq"], A_path)
                write_fasta("B", fc["B_seq"], B_path)
                # (no try/except: a tool failure is not a class of result)
                A_l = flank_landing(tgt_path, A_path, work, args.threads, "end")
                B_l = flank_landing(tgt_path, B_path, work, args.threads, "start")
                shutil.rmtree(work, ignore_errors=True)
                if A_l: Apos = f"{A_l[0]}:{A_l[2]}({A_l[1]})"
                if B_l: Bpos = f"{B_l[0]}:{B_l[2]}({B_l[1]})"
                if A_l is None or B_l is None:
                    mechanism = "ORIGIN_NOT_FOUND"          # a flank is absent/diverged
                elif A_l[0] != B_l[0]:
                    mechanism = "ORIGIN_ON_TWO_CONTIGS"     # assembly break or rearrangement
                elif A_l[1] != B_l[1]:
                    mechanism = "ORIGIN_REARRANGED"
                else:
                    # SIGNED distance from the end of A to the start of B along
                    # the flank orientation; abs() hid B landing upstream of A
                    dist = (B_l[2] - A_l[2]) if A_l[1] == "+" else (A_l[2] - B_l[2])
                    if dist < -500:
                        mechanism = "ORIGIN_REARRANGED"
                    elif abs(dist - fc["expected_with_DNA"]) <= args.copy_paste_margin:
                        mechanism = "COPY_PASTE"            # origin still occupied
                    elif dist <= 500:
                        mechanism = "CUT_PASTE"             # clean empty site
                    elif dist <= args.cut_paste_max:
                        mechanism = "CUT_PASTE_IS_LEFT"     # one IS copy remains at the origin
                    else:
                        mechanism = "AMBIGUOUS"
            n[mechanism] = n.get(mechanism, 0) + 1
            out.write(f"{r['ref_asm']}\t{r['ref_contig']}\t{r['ref_is_a_start']}\t"
                      f"{r['ref_is_b_start']}\t{r['segment_bp']}\t"
                      f"{r['tgt_asm']}\t{r['B_cov']}\t{r['B_ident']}\t"
                      f"{r['tgt_contig']}\t{r['tgt_pos']}\t"
                      f"{Apos}\t{Bpos}\t{dist}\t{fc['expected_with_DNA'] if fc else 'NA'}\t"
                      f"{mechanism}\n")

    ref_fa.close()
    shutil.rmtree(work_root, ignore_errors=True)
    print(f"\nDONE. {n}")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
