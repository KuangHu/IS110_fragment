#!/usr/bin/env python3
"""Strict canonical [A][IS][DNA][IS][B] → [A][IS][DNA][IS][DNA][IS][B] check
for IS-mediated duplication candidates.

For each Tier-4 strict_is_mediated_dup candidate (ref_id, tgt_asm):
  1. Locate ref's IS pair from the candidate (anchor positions).
  2. Extract ref slice: A_flank (5 kb upstream of IS_a) + [IS][DNA][IS] + B_flank.
  3. Align A_flank against tgt + B_flank against tgt.
  4. PASS canonical if:
        - A_flank forward-aligns to tgt at >= --min-flank-cov of A, >= --min-identity
        - B_flank forward-aligns to tgt at >= --min-flank-cov of B, >= --min-identity
        - DNA (between ref's IS pair) appears >= 2 times in tgt at >= --min-dna-cov + identity
        - IS hits exist near each DNA-copy boundary in tgt (uses is_hits.tsv)
"""
import argparse, csv, gzip, json, os, re, shutil, subprocess, tempfile
from collections import defaultdict
import pysam


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--candidates", required=True,
                   help="strict_duplication_catalogue.tsv (Tier-4)")
    p.add_argument("--dup-all-input", required=True,
                   help="rearrangements_duplication_all.tsv (has 'details' anchor positions)")
    p.add_argument("--records",      required=True)
    p.add_argument("--is-hits",      required=True)
    p.add_argument("--src-dir",      required=True)
    p.add_argument("--out",          required=True)
    p.add_argument("--flank",        type=int, default=5000)
    p.add_argument("--min-flank-cov",  type=float, default=0.80)
    p.add_argument("--min-dna-cov",    type=float, default=0.50)
    p.add_argument("--min-identity",   type=float, default=95.0)
    p.add_argument("--is-tol",         type=int,   default=3000)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--sample",  type=int, default=0)
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
    idx = defaultdict(lambda: defaultdict(list))
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            asm = r["assembly"]; contig = r["contig"]
            local = contig.split("|", 1)[1] if "|" in contig else contig
            hit = (int(r["tnp_start"]), int(r["tnp_end"]))
            idx[asm][local].append(hit)
            tail = local.split("|")[-1]
            if tail != local:
                idx[asm][tail].append(hit)
    return idx


def is_within(hits, pos, tol):
    for (s, e) in hits:
        d = 0 if s <= pos <= e else min(abs(pos - s), abs(pos - e))
        if d <= tol: return True
    return False


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


def best_forward_block(ref_path, query_path, work, threads):
    paf = os.path.join(work, "_aln.paf")
    subprocess.run(["minimap2", "-x", "asm10", "-c", "--eqx",
                    "-t", str(threads), ref_path, query_path, "-o", paf],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    best = None
    with open(paf) as fh:
        for line in fh:
            c = line.split("\t")
            if len(c) < 12: continue
            if c[4] != "+": continue
            matches, block_len = int(c[9]), int(c[10])
            if best is None or block_len > best[3]:
                ident = matches / block_len * 100 if block_len > 0 else 0
                best = (c[5], int(c[7]), int(c[8]), block_len, ident)
    return best


def all_forward_clusters(ref_path, query_path, work, threads, cluster_gap=10000):
    paf = os.path.join(work, "_aln_all.paf")
    subprocess.run(["minimap2", "-x", "asm10", "-c", "--eqx",
                    "-t", str(threads), ref_path, query_path, "-o", paf],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    by_t = defaultdict(list)
    qlen = None
    with open(paf) as fh:
        for line in fh:
            c = line.split("\t")
            if len(c) < 12: continue
            if c[4] != "+": continue
            qlen = int(c[1])
            by_t[c[5]].append({"qs": int(c[2]), "qe": int(c[3]),
                               "ts": int(c[7]), "te": int(c[8]),
                               "matches": int(c[9]), "block": int(c[10])})
    clusters = []
    for tname, blks in by_t.items():
        blks.sort(key=lambda x: x["ts"])
        cur = [blks[0]]
        groups = []
        for b in blks[1:]:
            if b["ts"] - cur[-1]["te"] <= cluster_gap:
                cur.append(b)
            else:
                groups.append(cur); cur = [b]
        groups.append(cur)
        for g in groups:
            ts = min(b["ts"] for b in g); te = max(b["te"] for b in g)
            q_iv = sorted([(b["qs"], b["qe"]) for b in g])
            merged_q = [list(q_iv[0])]
            for s, e in q_iv[1:]:
                if s <= merged_q[-1][1]:
                    merged_q[-1][1] = max(merged_q[-1][1], e)
                else:
                    merged_q.append([s, e])
            total_q = sum(e-s for s, e in merged_q)
            matches = sum(b["matches"] for b in g)
            block = sum(b["block"] for b in g)
            wid = matches/block*100 if block > 0 else 0
            clusters.append({"tname": tname, "ts": ts, "te": te,
                             "cov": total_q/qlen if qlen else 0,
                             "wid": wid})
    return clusters


ANCHOR_RE = re.compile(r"([\w.|]+):(\d+)-(\d+)\(([+-])\)")
def parse_anchors(details):
    return [(m.group(1), int(m.group(2)), int(m.group(3)), m.group(4))
            for m in ANCHOR_RE.finditer(details)]


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)
    loci = load_loci(args.records)
    is_hits = load_is_hits(args.is_hits)

    # index dup_all by (ref_id, tgt_asm) for the anchor details
    details_by_pair = {}
    with open(args.dup_all_input) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            details_by_pair[(r['ref_id'], r['assembly'])] = r.get('details', '')

    rows = list(csv.DictReader(open(args.candidates), delimiter="\t"))
    if args.sample > 0 and len(rows) > args.sample:
        import random
        random.seed(13)
        rows = random.sample(rows, args.sample)
    print(f"Checking {len(rows)} duplication candidates...", flush=True)

    out_path = os.path.join(args.out, "duplication_canonical_check.tsv")
    work_root = tempfile.mkdtemp(prefix="dupcc_")
    fa_cache = {}
    def open_asm(asm):
        if asm in fa_cache: return fa_cache[asm]
        path = find_assembly_fa(args.src_dir, asm, work_root)
        fa_cache[asm] = (path, pysam.FastaFile(path) if path else None)
        return fa_cache[asm]

    n_pass = 0
    with open(out_path, "w") as out:
        out.write("species\tref_id\ttgt_asm\tA_cov\tA_id\tA_pass\t"
                  "B_cov\tB_id\tB_pass\tn_dna_clusters\tDNA_best_cov\tDNA_best_id\t"
                  "n_is_flanked_clusters\tcanonical_pass\n")
        for i, r in enumerate(rows):
            if i % 100 == 0:
                print(f"  {i}/{len(rows)}", flush=True)
            ref_id = r['ref_id']
            tgt_asm = r['tgt_asm']
            species = r.get('species', 'klebsiella_pneumoniae')
            locus = loci.get(ref_id)
            if not locus: continue
            ref_contig, is_start, is_end = locus
            ref_asm = ref_id.split("|")[0]
            ref_path, ref_fa = open_asm(ref_asm)
            tgt_path, tgt_fa = open_asm(tgt_asm)
            if not ref_fa or not tgt_fa: continue
            try:
                clen = ref_fa.get_reference_length(ref_contig)
            except (KeyError, ValueError):
                continue

            # Get anchor positions
            details = details_by_pair.get((ref_id, tgt_asm), '')
            if not details: continue
            anchors = parse_anchors(details)
            if not anchors: continue
            bycontig = defaultdict(list)
            for (c, s, e, st) in anchors:
                bycontig[c].append((s, e, st))
            tgt_contig = max(bycontig.keys(), key=lambda k: len(bycontig[k]))

            # Extract ref's A flank, DNA, B flank
            # DNA is the IS-locus region in ref (we use the candidate IS ± window for DNA)
            A_s = max(0, is_start - args.flank); A_e = is_start
            B_s = is_end; B_e = min(clen, is_end + args.flank)
            DNA_s = max(0, is_start - args.flank); DNA_e = min(clen, is_end + args.flank)
            A_seq = ref_fa.fetch(ref_contig, A_s, A_e) if A_e > A_s else ""
            B_seq = ref_fa.fetch(ref_contig, B_s, B_e) if B_e > B_s else ""
            DNA_seq = ref_fa.fetch(ref_contig, is_start, is_end) if is_end > is_start else ""

            work = tempfile.mkdtemp(prefix=f"d{i}_", dir=work_root)
            A_path = os.path.join(work, "A.fa")
            B_path = os.path.join(work, "B.fa")
            DNA_path = os.path.join(work, "DNA.fa")
            if A_seq: write_fasta("A", A_seq, A_path)
            if B_seq: write_fasta("B", B_seq, B_path)
            if DNA_seq: write_fasta("DNA", DNA_seq, DNA_path)
            if not (A_seq and B_seq and DNA_seq):
                shutil.rmtree(work, ignore_errors=True); continue

            # Test A and B forward-alignment against tgt
            try:
                A_blk = best_forward_block(tgt_path, A_path, work, args.threads)
                B_blk = best_forward_block(tgt_path, B_path, work, args.threads)
                DNA_clusters = all_forward_clusters(tgt_path, DNA_path, work, args.threads)
            except Exception:
                shutil.rmtree(work, ignore_errors=True); continue

            A_cov = A_blk[3]/len(A_seq) if A_blk else 0
            A_id  = A_blk[4] if A_blk else 0
            B_cov = B_blk[3]/len(B_seq) if B_blk else 0
            B_id  = B_blk[4] if B_blk else 0
            A_pass = A_cov >= args.min_flank_cov and A_id >= args.min_identity
            B_pass = B_cov >= args.min_flank_cov and B_id >= args.min_identity

            passing_dna = [c for c in DNA_clusters
                           if c["cov"] >= args.min_dna_cov and c["wid"] >= args.min_identity]
            best_dna = max(passing_dna, key=lambda c: c["cov"]) if passing_dna else None
            best_dna_cov = best_dna["cov"] if best_dna else 0
            best_dna_id  = best_dna["wid"] if best_dna else 0

            # Check IS-flanking of each DNA copy
            tgt_is_hits = is_hits.get(tgt_asm, {}).get(tgt_contig, [])
            n_is_flanked = 0
            for c in passing_dna:
                # need IS at both edges of the DNA cluster
                # extract local contig from tname
                local = c["tname"].split("|", 1)[1] if "|" in c["tname"] else c["tname"]
                hits = is_hits.get(tgt_asm, {}).get(local, [])
                if is_within(hits, c["ts"], args.is_tol) and is_within(hits, c["te"], args.is_tol):
                    n_is_flanked += 1

            canonical_pass = (A_pass and B_pass and len(passing_dna) >= 2 and n_is_flanked >= 2)
            if canonical_pass: n_pass += 1

            out.write(f"{species}\t{ref_id}\t{tgt_asm}\t"
                      f"{A_cov:.3f}\t{A_id:.1f}\t{A_pass}\t"
                      f"{B_cov:.3f}\t{B_id:.1f}\t{B_pass}\t"
                      f"{len(passing_dna)}\t{best_dna_cov:.3f}\t{best_dna_id:.1f}\t"
                      f"{n_is_flanked}\t{canonical_pass}\n")
            shutil.rmtree(work, ignore_errors=True)

    for _, fa in fa_cache.values():
        if fa: fa.close()
    shutil.rmtree(work_root, ignore_errors=True)
    print(f"\nDONE: {n_pass}/{len(rows)} pass canonical duplication check")
    print(f"  (A and B flanks forward-aligned at >= {args.min_flank_cov} cov, "
          f">= 2 DNA copies forward, both IS-flanked)")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
