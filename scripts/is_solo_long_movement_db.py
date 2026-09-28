#!/usr/bin/env python3
"""Solo long-IS movement detector across a species DB.

For each long full-length IS in ref (tnp_len >= min, both Pfam domains present):
  1. Extract the IS body + 5 kb upstream (A_flank) + 5 kb downstream (B_flank).
  2. Bulk-align all IS bodies vs species DB index in one minimap2 call.
  3. Cluster nearby forward blocks per (qname, tname), filter by cov / identity.
  4. Drop self-hits (cluster tgt assembly == ref assembly).
  5. For each passing cluster, do the flank check (A_flank vs tgt-upstream,
     B_flank vs tgt-downstream).
  6. Verdict:
        SYNTENIC      both flanks match  (= same context, not mobilized)
        MOVED         neither matches    (= IS body in a different neighborhood)
        PARTIAL_MATCH exactly one matches

Strict layer: cov_IS >= 0.95 AND identity_IS >= 99% (essentially identical IS).
              Loose layer: cov_IS >= 0.50 AND identity_IS >= 95% (per IS110 thresh).
"""
import argparse, csv, os, shutil, subprocess, tempfile
from collections import defaultdict
import pysam


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--ref-asm", required=True)
    p.add_argument("--ref-fa",  required=True,
                   help="reference FASTA (indexed; same DB the ref lives in is fine)")
    p.add_argument("--species-db", required=True,
                   help="species concatenated FASTA, indexed with .fai")
    p.add_argument("--species-mmi", default="",
                   help="optional pre-built .mmi for the species DB")
    p.add_argument("--is-hits",  required=True)
    p.add_argument("--out",      required=True)
    p.add_argument("--min-tnp-len", type=int, default=1000)
    p.add_argument("--require-both-domains", action="store_true", default=True)
    p.add_argument("--flank",            type=int,   default=5000)
    p.add_argument("--is-pad",           type=int,   default=50,
                   help="extra bp on each side of tnp ORF when extracting IS body")
    p.add_argument("--min-coverage",     type=float, default=0.50)
    p.add_argument("--min-identity",     type=float, default=95.0)
    p.add_argument("--strict-coverage",  type=float, default=0.95)
    p.add_argument("--strict-identity",  type=float, default=99.0)
    p.add_argument("--flank-match-cov",  type=float, default=0.50)
    p.add_argument("--flank-match-id",   type=float, default=90.0)
    p.add_argument("--cluster-gap",      type=int,   default=2000)
    p.add_argument("--threads", type=int, default=8)
    return p.parse_args()


def load_long_is(path, ref_asm, min_tnp_len, require_both):
    """Return list of (contig, tnp_start, tnp_end, tnp_strand, tnp_len, is_id)."""
    out = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["assembly"] != ref_asm: continue
            try: tlen = int(r["tnp_len"])
            except ValueError: continue
            if tlen < min_tnp_len: continue
            dom = r.get("domains_hit", "")
            if require_both and not ("PF01548" in dom and "PF02371" in dom):
                continue
            out.append((r["contig"], int(r["tnp_start"]), int(r["tnp_end"]),
                        r.get("tnp_strand", "+"), tlen, r["is_id"]))
    return out


def write_fasta(label, seq, path):
    with open(path, "w") as fh:
        fh.write(f">{label}\n")
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i+80] + "\n")


def parse_clusters(paf_path, qlen_lookup, cluster_gap):
    """Parse PAF and cluster nearby forward blocks per (qname, tname)."""
    by_qt = defaultdict(list)
    with open(paf_path) as fh:
        for line in fh:
            c = line.split("\t")
            if len(c) < 12: continue
            if c[4] != "+": continue
            qname = c[0]; tname = c[5]
            by_qt[(qname, tname)].append({
                "qs": int(c[2]), "qe": int(c[3]),
                "ts": int(c[7]), "te": int(c[8]),
                "matches": int(c[9]), "block": int(c[10]),
            })
    clusters = []
    for (qname, tname), blks in by_qt.items():
        blks.sort(key=lambda x: x["ts"])
        cur = [blks[0]]; groups = []
        for b in blks[1:]:
            if b["ts"] - cur[-1]["te"] <= cluster_gap: cur.append(b)
            else: groups.append(cur); cur = [b]
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
    work_root = tempfile.mkdtemp(prefix="solomov_")

    long_is = load_long_is(args.is_hits, args.ref_asm,
                           args.min_tnp_len, args.require_both_domains)
    if not long_is:
        print(f"[{args.ref_asm}] no long IS found; nothing to do.", flush=True)
        shutil.rmtree(work_root, ignore_errors=True); return

    ref_fa = pysam.FastaFile(args.ref_fa)
    db_fa  = pysam.FastaFile(args.species_db)
    print(f"[{args.ref_asm}] {len(long_is)} long IS in ref", flush=True)

    # Build candidates: IS body + A_flank + B_flank
    candidates = []
    for (contig, ts, te, strand, tlen, is_id) in long_is:
        try:
            clen = ref_fa.get_reference_length(contig)
        except (KeyError, ValueError):
            alt = f"{args.ref_asm}|{contig}"
            try:
                clen = ref_fa.get_reference_length(alt); contig = alt
            except (KeyError, ValueError):
                continue
        is_s = max(0, ts - args.is_pad); is_e = min(clen, te + args.is_pad)
        lf_s = max(0, ts - args.flank); lf_e = ts
        rf_s = te; rf_e = min(clen, te + args.flank)
        is_seq = ref_fa.fetch(contig, is_s, is_e)
        lf_seq = ref_fa.fetch(contig, lf_s, lf_e) if lf_e > lf_s else ""
        rf_seq = ref_fa.fetch(contig, rf_s, rf_e) if rf_e > rf_s else ""
        if not is_seq or not lf_seq or not rf_seq: continue
        cid = f"is{len(candidates)}"
        candidates.append({
            "id": cid, "is_id": is_id, "contig": contig,
            "is_start": ts, "is_end": te, "is_strand": strand, "tnp_len": tlen,
            "is_seq": is_seq, "lf": lf_seq, "rf": rf_seq,
        })
    print(f"[{args.ref_asm}] {len(candidates)} solo-IS candidates extracted", flush=True)
    if not candidates:
        ref_fa.close(); db_fa.close()
        shutil.rmtree(work_root, ignore_errors=True); return

    # One minimap2 call: all IS bodies vs species DB
    is_fa = os.path.join(work_root, "is_bodies.fa")
    qlen_lookup = {}
    with open(is_fa, "w") as fh:
        for c in candidates:
            fh.write(f">{c['id']}\n")
            seq = c["is_seq"]; qlen_lookup[c["id"]] = len(seq)
            for i in range(0, len(seq), 80): fh.write(seq[i:i+80] + "\n")

    paf_path = os.path.join(work_root, "all.paf")
    db_input = args.species_mmi if args.species_mmi else args.species_db
    print(f"[{args.ref_asm}] minimap2 IS bodies vs species DB ({db_input})...", flush=True)
    subprocess.run(["minimap2", "-x", "asm10", "-c", "--eqx",
                    "-t", str(args.threads), db_input, is_fa, "-o", paf_path],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    print(f"[{args.ref_asm}] PAF: {os.path.getsize(paf_path)/1e6:.1f} MB", flush=True)

    clusters = parse_clusters(paf_path, qlen_lookup, args.cluster_gap)
    print(f"[{args.ref_asm}] {len(clusters)} total forward clusters", flush=True)

    passing = [c for c in clusters
               if c["cov"] >= args.min_coverage
               and c["weighted_ident"] >= args.min_identity
               and not c["tname"].startswith(args.ref_asm + "|")]
    print(f"[{args.ref_asm}] {len(passing)} clusters pass cov>={args.min_coverage} "
          f"id>={args.min_identity}% (excl self)", flush=True)

    cand_by_id = {c["id"]: c for c in candidates}

    # Per-candidate flank FASTAs written once
    flank_dir = os.path.join(work_root, "flanks"); os.makedirs(flank_dir, exist_ok=True)
    for c in candidates:
        write_fasta(f"{c['id']}_lf", c["lf"], os.path.join(flank_dir, f"{c['id']}_lf.fa"))
        write_fasta(f"{c['id']}_rf", c["rf"], os.path.join(flank_dir, f"{c['id']}_rf.fa"))

    out_path = os.path.join(args.out, f"{args.ref_asm}_solo_movement.tsv")
    n = {"MOVED":0, "SYNTENIC":0, "PARTIAL_MATCH":0}
    n_strict = {"MOVED":0, "SYNTENIC":0, "PARTIAL_MATCH":0}

    with open(out_path, "w") as out:
        out.write("ref_asm\tref_is_id\tref_contig\tis_start\tis_end\ttnp_len\t"
                  "tgt_asm\ttgt_contig\tts\tte\tn_blocks\tIS_cov\tIS_ident\t"
                  "left_flank_match\tright_flank_match\tverdict\tstrict_IS\n")
        for ki, cl in enumerate(passing):
            if ki % 500 == 0:
                print(f"  flank-checking {ki}/{len(passing)}", flush=True)
            cand = cand_by_id[cl["qname"]]
            tname = cl["tname"]
            if "|" not in tname: continue
            tgt_asm, tgt_contig_local = tname.split("|", 1)
            try: tlen = db_fa.get_reference_length(tname)
            except (KeyError, ValueError): continue
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
                try: b = best_forward_block_simple(refp, tgtp, work, args.threads)
                except Exception: return (False, 0, 0)
                if b is None: return (False, 0, 0)
                _, _, _, bp2, id2 = b
                cov2 = bp2 / args.flank
                return (cov2 >= args.flank_match_cov and id2 >= args.flank_match_id,
                        cov2, id2)

            lf_ok, lf_cov, lf_id = flanks_match(
                os.path.join(flank_dir, f"{cand['id']}_lf.fa"), tlf_path)
            rf_ok, rf_cov, rf_id = flanks_match(
                os.path.join(flank_dir, f"{cand['id']}_rf.fa"), trf_path)
            if lf_ok and rf_ok: verdict = "SYNTENIC"
            elif not lf_ok and not rf_ok: verdict = "MOVED"
            else: verdict = "PARTIAL_MATCH"
            n[verdict] += 1

            strict = (cl["cov"] >= args.strict_coverage and
                      cl["weighted_ident"] >= args.strict_identity)
            if strict: n_strict[verdict] += 1

            out.write(f"{args.ref_asm}\t{cand['is_id']}\t{cand['contig']}\t"
                      f"{cand['is_start']}\t{cand['is_end']}\t{cand['tnp_len']}\t"
                      f"{tgt_asm}\t{tgt_contig_local}\t{ts}\t{te}\t{cl['n_blocks']}\t"
                      f"{cl['cov']:.3f}\t{cl['weighted_ident']:.1f}\t"
                      f"{lf_ok}({lf_cov:.2f}/{lf_id:.0f})\t{rf_ok}({rf_cov:.2f}/{rf_id:.0f})\t"
                      f"{verdict}\t{strict}\n")
            shutil.rmtree(work, ignore_errors=True)

    ref_fa.close(); db_fa.close()
    shutil.rmtree(work_root, ignore_errors=True)
    print(f"[{args.ref_asm}] DONE. loose: {n}    strict: {n_strict}")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
