#!/usr/bin/env python3
"""Reciprocal X-C test for IS-mediated HR translocations.

For a Tier-3b case where we've already confirmed:
  - REF: [A]-[IS_ref]-[B]
  - TGT contig 1: [A]-[IS_tgt_1]-[C]
  - TGT contig 2: [X]-[IS_tgt_2]-[B]
  (= bilateral A and B conservation)

A true HR translocation predicts the RECIPROCAL exchange: the [C] and [X]
segments displaced A and B in target should be findable in ref, adjacent to
ANOTHER IS copy in ref (the partner IS that participated in the HR).

This script:
  1. Identifies IS_tgt_1 and IS_tgt_2 positions in target (nearest IS to each anchor).
  2. Extracts tgt's C = 30 kb downstream of IS_tgt_1 (on the contig opposite to A).
  3. Extracts tgt's X = 30 kb upstream of IS_tgt_2 (on the contig opposite to B).
  4. Aligns C and X against the WHOLE ref assembly.
  5. For each best-forward-block hit in ref, checks distance to nearest IS in ref.

PASS_strict: BOTH C and X align to ref at ≥ 80 % cov, ≥ 95 % identity, AND
             both alignment positions are within 10 kb of an IS in ref.
PASS_lenient: same but only ONE of C/X needs to be IS-adjacent in ref.
"""
import argparse, csv, gzip, json, os, re, shutil, subprocess, tempfile
from collections import defaultdict
import pysam


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--candidates", required=True,
                   help="strict_hr_translocation_catalogue.tsv (Tier-3b passes)")
    p.add_argument("--records-runs", required=True,
                   help="cross_ref_is_runs/results — for per-species records.json + is_hits.tsv")
    p.add_argument("--src-dirs",  required=True,
                   help="JSON file: species → source FASTA dir")
    p.add_argument("--candidates-tier3a", required=True,
                   help="strict_translocation_catalogue.tsv — for original anchor details")
    p.add_argument("--out", required=True)
    p.add_argument("--proxy-window", type=int, default=30000)
    p.add_argument("--min-block-cov", type=float, default=0.80)
    p.add_argument("--min-identity",  type=float, default=95.0)
    p.add_argument("--is-tol", type=int, default=10000,
                   help="max distance from C/X alignment block to nearest IS in ref")
    p.add_argument("--threads", type=int, default=8)
    return p.parse_args()


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
    """assembly -> contig (stripped: 'ENA|XXX|XXX.1' or just 'XXX.1') -> [(s, e, strand)].
    Anchor details and PAF parsing store contig in stripped form (assembly prefix removed),
    so we strip here too. Also keep an alias under the bare last token for ecoli-style names."""
    idx = defaultdict(lambda: defaultdict(list))
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            asm = r["assembly"]; contig = r["contig"]
            local = contig.split("|", 1)[1] if "|" in contig else contig
            hit = (int(r["tnp_start"]), int(r["tnp_end"]),
                   r.get("tnp_strand", "+"))
            idx[asm][local].append(hit)
            # also store under the very-last token (works for "CP044314.1"-style ecoli)
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


def find_nearest_is(is_list, pos, max_dist=5000):
    best = None
    for (s, e, st) in is_list:
        d = 0 if s <= pos <= e else min(abs(pos - s), abs(pos - e))
        if d <= max_dist and (best is None or d < best[3]):
            best = (s, e, st, d)
    return best


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


def best_forward_block(ref_index_fa, query_fa, work, threads):
    """Align query against the whole ref assembly. Return (qname, tname, ts, te,
    block_bp, ident) for the best forward-strand block, or None."""
    paf = os.path.join(work, "aln.paf")
    subprocess.run(["minimap2", "-x", "asm10", "-c", "--eqx",
                    "-t", str(threads), ref_index_fa, query_fa, "-o", paf],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    best = None
    with open(paf) as fh:
        for line in fh:
            c = line.split("\t")
            if len(c) < 12: continue
            if c[4] != "+": continue
            matches, block_len = int(c[9]), int(c[10])
            ident = matches / block_len * 100 if block_len > 0 else 0
            if best is None or block_len > best[4]:
                best = (c[0], c[5], int(c[7]), int(c[8]), block_len, ident)
    return best


ANCHOR_RE = re.compile(r"([\w.|]+):(\d+)-(\d+)\(([+-])\)")
def parse_anchors_block(details, side):
    m = re.search(rf"{side}=\[([^\]]*)\]", details)
    if not m: return []
    return [(mm.group(1), int(mm.group(2)), int(mm.group(3)), mm.group(4))
            for mm in ANCHOR_RE.finditer(m.group(1))]


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)
    src_dirs = json.load(open(args.src_dirs))

    # Index Tier-3a candidates by (ref_id, tgt_asm) for the details column
    tier3a_details = {}
    with open(args.candidates_tier3a) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tier3a_details[(r['ref_id'], r['tgt_asm'])] = r.get('details', '')

    # Per-species lazy-load caches
    loci_by_sp = {}
    is_by_sp = {}
    def loci_for(sp):
        if sp not in loci_by_sp:
            loci_by_sp[sp] = load_loci(os.path.join(args.records_runs, sp, "records", "records.json"))
        return loci_by_sp[sp]
    def is_for(sp):
        if sp not in is_by_sp:
            is_by_sp[sp] = load_is_hits(os.path.join(args.records_runs, sp, "is_hits.tsv"))
        return is_by_sp[sp]

    rows = list(csv.DictReader(open(args.candidates), delimiter="\t"))
    print(f"Reciprocal-testing {len(rows)} Tier-3b candidates...", flush=True)
    out_path = os.path.join(args.out, "reciprocal_results.tsv")
    work_root = tempfile.mkdtemp(prefix="recip_")
    fa_cache = {}
    def open_asm(species, asm):
        key = (species, asm)
        if key in fa_cache: return fa_cache[key]
        path = find_assembly_fa(src_dirs[species], asm, work_root)
        fa_cache[key] = (path, pysam.FastaFile(path) if path else None)
        return fa_cache[key]

    min_block_bp = int(args.proxy_window * args.min_block_cov)

    n_both_is_adjacent = 0  # PASS_strict
    n_one_is_adjacent  = 0
    n_both_align       = 0  # both C and X align to ref at ≥thresh, regardless of IS adjacency
    n_neither_align    = 0

    with open(out_path, "w") as out:
        out.write("species\tref_id\ttgt_asm\tC_block_bp\tC_ident\tC_align_pos\tC_dist_to_is\t"
                  "X_block_bp\tX_ident\tX_align_pos\tX_dist_to_is\t"
                  "C_align_passes\tX_align_passes\tC_is_adjacent\tX_is_adjacent\tboth_is_adjacent\n")
        for i, r in enumerate(rows):
            if i % 10 == 0:
                print(f"  {i}/{len(rows)}", flush=True)
            species = r['species']
            ref_id  = r['ref_id']
            tgt_asm = r['tgt_asm']
            details = tier3a_details.get((ref_id, tgt_asm), '')
            if not details: continue
            ups   = parse_anchors_block(details, "up")
            downs = parse_anchors_block(details, "down")
            if not ups or not downs: continue
            up_c, up_s, up_e, _ = ups[0]
            down_c, down_s, down_e, _ = downs[0]

            ih = is_for(species)
            ref_asm = ref_id.split("|")[0]
            ref_path, ref_fa = open_asm(species, ref_asm)
            _, tgt_fa = open_asm(species, tgt_asm)
            if not ref_fa or not tgt_fa: continue

            # Find IS_tgt_1 near up_anchor, IS_tgt_2 near down_anchor
            tgt_is_at_up   = find_nearest_is(merge_is(ih[tgt_asm].get(up_c, [])),
                                              (up_s + up_e) // 2, max_dist=5000)
            tgt_is_at_down = find_nearest_is(merge_is(ih[tgt_asm].get(down_c, [])),
                                              (down_s + down_e) // 2, max_dist=5000)
            if tgt_is_at_up is None or tgt_is_at_down is None:
                continue

            # C = 30 kb DOWNSTREAM of IS_tgt_1 (= away from the up_anchor side of IS_tgt_1)
            # X = 30 kb UPSTREAM of IS_tgt_2 (= away from the down_anchor side of IS_tgt_2)
            # Direction: up_anchor matches ref A which is upstream of ref IS; on target contig 1,
            # the up_anchor sits NEAR IS_tgt_1; we need the OPPOSITE side of IS_tgt_1 — i.e., the side
            # away from up_anchor — which carries C.
            up_anchor_mid = (up_s + up_e) // 2
            if up_anchor_mid < tgt_is_at_up[0]:
                # up_anchor is upstream of IS_tgt_1 → C = downstream of IS_tgt_1
                C_s = tgt_is_at_up[1]
                C_e = C_s + args.proxy_window
            else:
                # up_anchor is downstream → C = upstream of IS_tgt_1
                C_e = tgt_is_at_up[0]
                C_s = max(0, C_e - args.proxy_window)
            down_anchor_mid = (down_s + down_e) // 2
            if down_anchor_mid > tgt_is_at_down[1]:
                # down_anchor is downstream of IS_tgt_2 → X = upstream of IS_tgt_2
                X_e = tgt_is_at_down[0]
                X_s = max(0, X_e - args.proxy_window)
            else:
                X_s = tgt_is_at_down[1]
                X_e = X_s + args.proxy_window

            try:
                clen_up = tgt_fa.get_reference_length(up_c)
                clen_dn = tgt_fa.get_reference_length(down_c)
            except (KeyError, ValueError):
                continue
            C_e = min(C_e, clen_up); X_e = min(X_e, clen_dn)
            C_seq = tgt_fa.fetch(up_c,   C_s, C_e) if C_e > C_s else ""
            X_seq = tgt_fa.fetch(down_c, X_s, X_e) if X_e > X_s else ""

            work = tempfile.mkdtemp(prefix="rc_", dir=work_root)
            C_path = os.path.join(work, "tgt_C.fa")
            X_path = os.path.join(work, "tgt_X.fa")
            if C_seq: write_fasta("TGT_C", C_seq, C_path)
            if X_seq: write_fasta("TGT_X", X_seq, X_path)

            C_block = X_block = None
            try:
                if C_seq:
                    C_block = best_forward_block(ref_path, C_path, work, args.threads)
                if X_seq:
                    X_block = best_forward_block(ref_path, X_path, work, args.threads)
            except Exception:
                pass

            ref_is_full = ih.get(ref_asm, {})

            def block_summary(blk, ref_len):
                if blk is None: return (0, 0.0, "", 999999, False, False)
                qname, tname, ts, te, bp, ident = blk
                # local contig name for IS lookup
                contig_full = tname
                hits = merge_is(ref_is_full.get(contig_full, []))
                if not hits and "|" in tname:
                    # try stripping
                    local = tname.split("|", 1)[1]
                    hits = merge_is(ref_is_full.get(local, []))
                if not hits and not tname.startswith(ref_asm + "|"):
                    hits = merge_is(ref_is_full.get(f"{ref_asm}|{tname}", []))
                aligns = (bp >= min_block_bp and ident >= args.min_identity)
                nis = find_nearest_is(hits, (ts + te) // 2, max_dist=args.is_tol)
                is_adj = nis is not None
                d_to_is = nis[3] if nis else 999999
                return (bp, ident, f"{tname}:{ts}-{te}", d_to_is, aligns, is_adj)

            C_bp, C_id, C_pos, C_d, C_ok, C_isadj = block_summary(C_block, args.proxy_window)
            X_bp, X_id, X_pos, X_d, X_ok, X_isadj = block_summary(X_block, args.proxy_window)
            both_isadj = (C_ok and C_isadj and X_ok and X_isadj)

            if C_ok and X_ok: n_both_align += 1
            elif not C_ok and not X_ok: n_neither_align += 1
            if both_isadj: n_both_is_adjacent += 1
            elif (C_ok and C_isadj) or (X_ok and X_isadj): n_one_is_adjacent += 1

            out.write(f"{species}\t{ref_id}\t{tgt_asm}\t"
                      f"{C_bp}\t{C_id:.1f}\t{C_pos}\t{C_d}\t"
                      f"{X_bp}\t{X_id:.1f}\t{X_pos}\t{X_d}\t"
                      f"{C_ok}\t{X_ok}\t{C_isadj}\t{X_isadj}\t{both_isadj}\n")
            shutil.rmtree(work, ignore_errors=True)

    print(f"\nDONE: reciprocal X-C scoring on {len(rows)} cases")
    print(f"  both C+X align to ref at >= {args.min_block_cov} cov / >= {args.min_identity}% ident: {n_both_align}")
    print(f"  both align AND both adjacent to a ref IS within {args.is_tol} bp: {n_both_is_adjacent}")
    print(f"  only one C/X is IS-adjacent: {n_one_is_adjacent}")
    print(f"  neither C nor X aligns at threshold: {n_neither_align}")
    print(f"Wrote {out_path}")

    # cleanup
    for _, fa in fa_cache.values():
        if fa: fa.close()
    shutil.rmtree(work_root, ignore_errors=True)


if __name__ == "__main__":
    main()
