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
import argparse, csv, gzip, json, os, re, shutil, subprocess, tempfile, sys
from collections import defaultdict
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib_rearr as lr  # noqa: E402


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

    # CANONICAL = on ONE target contig and strand:
    #     O-flank  [IS]  D  [IS]  D  [IS]   (D = the duplicated flank, O = the other)
    #   * the duplicated flank of the element is present at >= 2 separate places
    #   * an IS copy lies between consecutive copies
    #   * the other flank is present once
    # The previous test searched for the IS ELEMENT ITSELF across the whole
    # target and asked whether each copy was near an IS hit (always true): a
    # target with no duplication but two IS copies was 'canonical'.
    loci = lr.load_loci(args.records)
    n_pass = n_assessed = 0
    n_skip = defaultdict(int)
    with open(out_path, "w") as out:
        out.write("species\tref_id\ttgt_asm\tA_cov\tA_id\tA_pass\t"
                  "B_cov\tB_id\tB_pass\tn_dna_clusters\tDNA_best_cov\tDNA_best_id\t"
                  "n_is_flanked_clusters\tcanonical_pass\tdup_side\ttgt_contig\n")
        for i, r in enumerate(rows):
            if i % 100 == 0:
                print(f"  {i}/{len(rows)}", flush=True)
            ref_id = r['ref_id']
            tgt_asm = r['tgt_asm']
            if 'species' not in r:
                raise SystemExit("--candidates has no 'species' column (it used to be "
                                 "silently labelled klebsiella_pneumoniae)")
            species = r['species']
            locus = loci.get(ref_id)
            if not locus:
                n_skip["ref_id not in records"] += 1; continue
            ref_contig, is_s0, is_e, tnp_strand = locus
            _, ref_fa = open_asm(ref_id.split("|")[0])
            _, tgt_fa = open_asm(tgt_asm)
            if not ref_fa or not tgt_fa:
                n_skip["assembly FASTA not found"] += 1; continue
            d = lr.parse_details(details_by_pair.get((ref_id, tgt_asm), ''))
            anchors = d["up"] + d["down"]
            groups = defaultdict(list)
            for (c, s_, e_, st) in anchors:
                groups[(c, st)].append((s_, e_))
            groups = {k: v for k, v in groups.items() if len(set(v)) >= 2}
            if not groups:
                n_skip["no contig with two same-strand anchor hits"] += 1; continue
            (tgt_contig, st), pos = max(sorted(groups.items()), key=lambda kv: len(set(kv[1])))
            try:
                A_seq, B_seq = lr.ref_proxies(ref_fa, ref_contig, is_s0, is_e, tnp_strand,
                                              args.flank)
            except (KeyError, ValueError):
                n_skip["ref contig not in FASTA"] += 1; continue
            if min(len(A_seq), len(B_seq)) < 1000:
                n_skip["reference flank < 1 kb (contig end)"] += 1; continue
            tseq, toff = lr.fetch_window(tgt_fa, tgt_contig, min(p[0] for p in pos),
                                         max(p[1] for p in pos), args.flank + 30000)
            if not tseq:
                n_skip["target contig not in FASTA"] += 1; continue
            hits = (is_hits.get(tgt_asm, {}).get(tgt_contig)
                    or is_hits.get(tgt_asm, {}).get(tgt_contig.split("|", 1)[-1]) or [])

            def non_is(blocks):
                keep = []
                for b in blocks:
                    ts_, te_ = b["ts"] + toff, b["te"] + toff
                    is_bp = lr.union_bp([(max(ts_, h[0] - 300), min(te_, h[1] + 300))
                                         for h in hits if min(te_, h[1] + 300) > max(ts_, h[0] - 300)])
                    if is_bp < 0.8 * (te_ - ts_):
                        keep.append(b)
                return keep

            work = tempfile.mkdtemp(prefix=f"d{i}_", dir=work_root)
            res = {}
            for name, seq in (("A", A_seq), ("B", B_seq)):
                bl = non_is(lr.align(seq, tseq, work, args.threads))
                dup_bp, copies, places = lr.copy_depth(bl, 1000, args.min_identity, strand=st)
                res[name] = {"cov": lr.strand_cover(bl, st, args.min_identity) / len(seq),
                             "id": max((x["ident"] for x in bl if x["strand"] == st), default=0.0),
                             "dup_bp": dup_bp, "copies": copies, "places": places,
                             "len": len(seq)}
            shutil.rmtree(work, ignore_errors=True)
            dup = max("AB", key=lambda k: (res[k]["copies"], res[k]["dup_bp"]))
            other = "B" if dup == "A" else "A"
            pl = sorted(res[dup]["places"], key=lambda x: x["ts"])
            n_is_between = 0
            for x, y in zip(pl, pl[1:]):
                lo, hi = x["te"] + toff - args.is_tol, y["ts"] + toff + args.is_tol
                if any(h[0] >= lo and h[1] <= hi for h in hits):
                    n_is_between += 1
            A_pass = res["A"]["cov"] >= args.min_flank_cov
            B_pass = res["B"]["cov"] >= args.min_flank_cov
            canonical_pass = (res[dup]["copies"] >= 2
                              and res[dup]["dup_bp"] >= args.min_dna_cov * min(res[dup]["len"], 5000)
                              and n_is_between >= 1
                              and res[other]["cov"] >= args.min_flank_cov)
            n_assessed += 1
            if canonical_pass: n_pass += 1
            out.write(f"{species}\t{ref_id}\t{tgt_asm}\t"
                      f"{res['A']['cov']:.3f}\t{res['A']['id']:.1f}\t{A_pass}\t"
                      f"{res['B']['cov']:.3f}\t{res['B']['id']:.1f}\t{B_pass}\t"
                      f"{res[dup]['copies']}\t{res[dup]['dup_bp'] / res[dup]['len']:.3f}\t"
                      f"{res[dup]['id']:.1f}\t{n_is_between}\t{canonical_pass}\t"
                      f"{'up' if dup == 'A' else 'down'}\t{tgt_contig}\n")

    for _, fa in fa_cache.values():
        if fa: fa.close()
    shutil.rmtree(work_root, ignore_errors=True)
    print(f"\nDONE: {n_pass}/{n_assessed} assessed candidates pass the canonical duplication "
          f"check; NOT assessed (no row): {dict(n_skip)}")
    print("  (one flank present at >= 2 places on one contig and strand with an IS "
          "copy between them, the other flank present once)")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
