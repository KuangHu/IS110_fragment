#!/usr/bin/env python3
"""Scalable composite-transposon mobilization detector using the species DB.

For each ref candidate composite transposon B:
  1. minimap2 B against the concatenated species DB (= all assemblies in one FASTA).
  2. Each PAF block's tname is "ASSEMBLY|CONTIG" → assembly identifies the target.
  3. Cluster nearby blocks per (assembly, contig) → one cluster per location.
  4. For each cluster passing thresholds:
        a. Extract tgt's left + right flanks at the cluster position (from the
           indexed species DB).
        b. Compare to ref's A_flank and B_flank → SYNTENIC / MOVED / PARTIAL.

This is ~22,000× faster than aligning B to each target FASTA individually,
because minimap2 indexes the DB once and scans all hits in one pass.
"""
import argparse, csv, gzip, json, os, re, shutil, subprocess, tempfile
from collections import defaultdict
import pysam


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--ref-asm", required=True)
    p.add_argument("--ref-fa",  required=True,
                   help="reference FASTA file (uncompressed or .gz, must have .fai)")
    p.add_argument("--species-db", required=True,
                   help="concatenated species FASTA (uncompressed, indexed with .fai)")
    p.add_argument("--species-mmi", default="",
                   help="optional pre-built minimap2 index (.mmi) for the species DB; "
                        "if not given, uses --species-db directly (slower)")
    p.add_argument("--is-hits",  required=True)
    p.add_argument("--out",      required=True)
    p.add_argument("--min-segment-bp", type=int, default=10000)
    p.add_argument("--max-segment-bp", type=int, default=30000)
    p.add_argument("--flank",          type=int, default=5000)
    p.add_argument("--min-coverage",   type=float, default=0.50)
    p.add_argument("--min-identity",   type=float, default=80.0)
    p.add_argument("--flank-match-cov", type=float, default=0.50)
    p.add_argument("--flank-match-id",  type=float, default=90.0)
    p.add_argument("--cluster-gap",    type=int, default=10000)
    p.add_argument("--require-same-strand", action="store_true")
    p.add_argument("--threads", type=int, default=32)
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


def write_fasta(label, seq, path):
    with open(path, "w") as fh:
        fh.write(f">{label}\n")
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i+80] + "\n")


def parse_clusters(paf_path, qlen_lookup, cluster_gap):
    """Parse PAF and cluster nearby forward-strand blocks per (qname, tname).
    qlen_lookup: dict qname → query length."""
    by_qt = defaultdict(list)
    with open(paf_path) as fh:
        for line in fh:
            c = line.split("\t")
            if len(c) < 12: continue
            if c[4] != "+": continue
            qname = c[0]
            tname = c[5]
            by_qt[(qname, tname)].append({
                "qs": int(c[2]), "qe": int(c[3]),
                "ts": int(c[7]), "te": int(c[8]),
                "matches": int(c[9]), "block": int(c[10]),
            })
    clusters = []
    for (qname, tname), blks in by_qt.items():
        blks.sort(key=lambda x: x["ts"])
        cur = [blks[0]]
        groups = []
        for b in blks[1:]:
            if b["ts"] - cur[-1]["te"] <= cluster_gap:
                cur.append(b)
            else:
                groups.append(cur); cur = [b]
        groups.append(cur)
        qlen = qlen_lookup.get(qname, 0)
        for g in groups:
            ts = min(b["ts"] for b in g)
            te = max(b["te"] for b in g)
            q_iv = sorted([(b["qs"], b["qe"]) for b in g])
            merged_q = [list(q_iv[0])]
            for s, e in q_iv[1:]:
                if s <= merged_q[-1][1]:
                    merged_q[-1][1] = max(merged_q[-1][1], e)
                else:
                    merged_q.append([s, e])
            total_q = sum(e - s for s, e in merged_q)
            matches = sum(b["matches"] for b in g)
            block   = sum(b["block"]   for b in g)
            wid = matches / block * 100 if block > 0 else 0
            clusters.append({
                "qname": qname, "tname": tname, "ts": ts, "te": te,
                "cov": total_q / qlen if qlen else 0,
                "weighted_ident": wid, "n_blocks": len(g),
            })
    return clusters


def best_forward_block_simple(ref_path, query_path, work, threads):
    """One-vs-one alignment for flank checks (used after we have a cluster location)."""
    paf = os.path.join(work, "_flank.paf")
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


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)
    is_hits = load_is_hits(args.is_hits)

    work_root = tempfile.mkdtemp(prefix="movedb_")

    ref_fa = pysam.FastaFile(args.ref_fa)
    db_fa  = pysam.FastaFile(args.species_db)
    print(f"[{args.ref_asm}] loaded ref + species DB", flush=True)

    # Build composite-transposon candidates from ref
    candidates = []
    for contig, hits in is_hits.get(args.ref_asm, {}).items():
        merged = merge_is(hits)
        try:
            clen = ref_fa.get_reference_length(contig)
        except (KeyError, ValueError):
            try:
                clen = ref_fa.get_reference_length(f"{args.ref_asm}|{contig}")
                contig = f"{args.ref_asm}|{contig}"
            except (KeyError, ValueError):
                continue
        for i in range(len(merged)):
            for j in range(i+1, len(merged)):
                a, b = merged[i], merged[j]
                seg_len = b[0] - a[1]
                if seg_len < args.min_segment_bp or seg_len > args.max_segment_bp:
                    continue
                if args.require_same_strand and a[2] != b[2]:
                    continue
                lf_s = max(0, a[0] - args.flank); lf_e = a[0]
                rf_s = b[1]; rf_e = min(clen, b[1] + args.flank)
                seg_seq = ref_fa.fetch(contig, a[1], b[0])
                lf_seq  = ref_fa.fetch(contig, lf_s, lf_e) if lf_e > lf_s else ""
                rf_seq  = ref_fa.fetch(contig, rf_s, rf_e) if rf_e > rf_s else ""
                if not seg_seq or not lf_seq or not rf_seq:
                    continue
                cid = f"c{len(candidates)}"
                candidates.append({
                    "id": cid, "contig": contig, "is_a": a, "is_b": b,
                    "seg": seg_seq, "lf": lf_seq, "rf": rf_seq,
                    "expected_with_DNA": b[1] - a[0],
                })
    print(f"[{args.ref_asm}] {len(candidates)} composite-transposon candidates", flush=True)
    if not candidates:
        ref_fa.close(); db_fa.close()
        shutil.rmtree(work_root, ignore_errors=True)
        return

    # Write all candidate segments to a single multi-record FASTA
    seg_fa = os.path.join(work_root, "candidates.fa")
    qlen_lookup = {}
    with open(seg_fa, "w") as fh:
        for c in candidates:
            fh.write(f">{c['id']}\n")
            seq = c["seg"]
            qlen_lookup[c["id"]] = len(seq)
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + "\n")

    # One minimap2 call: candidates × species DB
    paf_path = os.path.join(work_root, "all.paf")
    db_input = args.species_mmi if args.species_mmi else args.species_db
    print(f"[{args.ref_asm}] minimap2 candidates vs species DB ({db_input})...", flush=True)
    subprocess.run(["minimap2", "-x", "asm10", "-c", "--eqx",
                    "-t", str(args.threads), db_input, seg_fa,
                    "-o", paf_path],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    print(f"[{args.ref_asm}] PAF generated: {os.path.getsize(paf_path)/1e6:.1f} MB", flush=True)

    clusters = parse_clusters(paf_path, qlen_lookup, args.cluster_gap)
    print(f"[{args.ref_asm}] {len(clusters)} total forward clusters across all tgts", flush=True)

    # Filter clusters by cov + identity
    passing = [c for c in clusters
               if c["cov"] >= args.min_coverage and c["weighted_ident"] >= args.min_identity]
    # Also exclude self-hits (ref hitting its own contig in the DB)
    passing = [c for c in passing
               if not c["tname"].startswith(args.ref_asm + "|")]
    print(f"[{args.ref_asm}] {len(passing)} clusters pass cov/id thresholds (excl self)", flush=True)

    # For each passing cluster, do the flank check
    cand_by_id = {c["id"]: c for c in candidates}
    out_path = os.path.join(args.out, f"{args.ref_asm}_movement.tsv")
    n = {"MOVED":0, "SYNTENIC":0, "PARTIAL_MATCH":0}

    # Write per-candidate flank FASTAs once for reuse
    flank_dir = os.path.join(work_root, "flanks")
    os.makedirs(flank_dir, exist_ok=True)
    for c in candidates:
        write_fasta(f"{c['id']}_lf", c["lf"], os.path.join(flank_dir, f"{c['id']}_lf.fa"))
        write_fasta(f"{c['id']}_rf", c["rf"], os.path.join(flank_dir, f"{c['id']}_rf.fa"))

    with open(out_path, "w") as out:
        out.write("ref_asm\tref_contig\tis_a_start\tis_b_start\tsegment_bp\t"
                  "tgt_asm\ttgt_contig\tts\tte\tn_blocks\tB_cov\tB_ident\t"
                  "left_flank_match\tright_flank_match\tverdict\n")
        for ki, cl in enumerate(passing):
            if ki % 100 == 0:
                print(f"  flank-checking {ki}/{len(passing)}", flush=True)
            cand = cand_by_id[cl["qname"]]
            tname = cl["tname"]
            if "|" not in tname: continue
            tgt_asm, tgt_contig_local = tname.split("|", 1)
            try:
                tlen = db_fa.get_reference_length(tname)
            except (KeyError, ValueError):
                continue
            ts, te = cl["ts"], cl["te"]
            tgt_lf_s = max(0, ts - args.flank); tgt_lf_e = ts
            tgt_rf_s = te; tgt_rf_e = min(tlen, te + args.flank)
            tgt_lf = db_fa.fetch(tname, tgt_lf_s, tgt_lf_e) if tgt_lf_e > tgt_lf_s else ""
            tgt_rf = db_fa.fetch(tname, tgt_rf_s, tgt_rf_e) if tgt_rf_e > tgt_rf_s else ""

            work = tempfile.mkdtemp(prefix=f"fc{ki}_", dir=work_root)
            tlf_path = os.path.join(work, "tlf.fa")
            trf_path = os.path.join(work, "trf.fa")
            if tgt_lf: write_fasta("TLF", tgt_lf, tlf_path)
            if tgt_rf: write_fasta("TRF", tgt_rf, trf_path)

            def flanks_match(refp, tgtp):
                if not os.path.exists(tgtp): return (False, 0, 0)
                try:
                    b = best_forward_block_simple(refp, tgtp, work, args.threads)
                except Exception:
                    return (False, 0, 0)
                if b is None: return (False, 0, 0)
                _, _, _, bp2, id2 = b
                cov2 = bp2 / args.flank
                return (cov2 >= args.flank_match_cov and id2 >= args.flank_match_id,
                        cov2, id2)
            lf_ok, lf_cov, lf_id = flanks_match(
                os.path.join(flank_dir, f"{cand['id']}_lf.fa"), tlf_path)
            rf_ok, rf_cov, rf_id = flanks_match(
                os.path.join(flank_dir, f"{cand['id']}_rf.fa"), trf_path)
            if lf_ok and rf_ok:
                verdict = "SYNTENIC"
            elif not lf_ok and not rf_ok:
                verdict = "MOVED"
            else:
                verdict = "PARTIAL_MATCH"
            n[verdict] += 1

            out.write(f"{args.ref_asm}\t{cand['contig']}\t{cand['is_a'][0]}\t"
                      f"{cand['is_b'][0]}\t{len(cand['seg'])}\t"
                      f"{tgt_asm}\t{tgt_contig_local}\t{ts}\t{te}\t{cl['n_blocks']}\t"
                      f"{cl['cov']:.3f}\t{cl['weighted_ident']:.1f}\t"
                      f"{lf_ok}({lf_cov:.2f}/{lf_id:.0f})\t{rf_ok}({rf_cov:.2f}/{rf_id:.0f})\t"
                      f"{verdict}\n")
            shutil.rmtree(work, ignore_errors=True)

    ref_fa.close(); db_fa.close()
    shutil.rmtree(work_root, ignore_errors=True)
    print(f"[{args.ref_asm}] DONE. {n}")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
