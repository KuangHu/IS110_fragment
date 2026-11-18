#!/usr/bin/env python3
"""Detect IS-mediated composite-transposon mobilization in bacterial genomes.

Biological model (composite transposon movement):

    REF:    --[A]--[IS_a →]--[ segment B ]--[IS_b →]--[C]--[D]--

    After HR between two same-orientation ISes (or transposase-mediated cut-paste),
    the [IS_a-B-IS_b] composite leaves its original position and reinserts elsewhere:

    TGT:    --[A]--[IS]--[C]--[IS]--[ B ]--[IS]--[D]--
                          ^^^^^^^                ^^^
                  B is GONE from its original location              B appears at a NEW location

The signature in target is: the segment B sequence appears at a DIFFERENT
chromosomal location with DIFFERENT immediate flanks than in ref.

Algorithm per ref assembly:
  1. Find pairs of same-orientation IS hits on the same ref contig, with
     inter-IS distance in [--min-segment-bp, --max-segment-bp].
  2. For each pair, extract the segment B (between IS_a end and IS_b start)
     plus ref left-flank (5 kb before IS_a) and right-flank (5 kb after IS_b).

For each (ref candidate B, tgt assembly) pair:
  3. minimap2 B against tgt assembly → find best forward block.
  4. PASS if block covers >= --min-coverage of |B| at >= --min-identity %.
  5. If pass, extract tgt's left and right flanks of B's landing position.
  6. Compare tgt's flanks to ref's flanks:
        - if MATCH (>= --flank-match-identity over >= 1 kb) → SYNTENIC (B in same place)
        - if DIFFER → MOVED (composite transposon mobilization confirmed)
        - if alignment of B fails → ABSENT (B not found in tgt)
"""
import argparse, csv, gzip, json, os, re, shutil, subprocess, sys, tempfile
from collections import defaultdict
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib_movement as lm  # noqa: E402
from lib_alleles import revcomp  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--ref-asm", required=True,
                   help="reference assembly accession (e.g., GCF_023066685.1)")
    p.add_argument("--tgt-list", required=True,
                   help="file with one tgt assembly per line")
    p.add_argument("--is-hits",  required=True,
                   help="is_hits.tsv (per-species IS positions, all assemblies)")
    p.add_argument("--src-dir",  required=True,
                   help="directory with assembly FASTAs (.fna / .fna.gz / .fa)")
    p.add_argument("--out",      required=True)
    p.add_argument("--min-segment-bp", type=int, default=10000)
    p.add_argument("--max-segment-bp", type=int, default=200000)
    p.add_argument("--flank",          type=int, default=5000)
    p.add_argument("--min-coverage",   type=float, default=0.90)
    p.add_argument("--min-identity",   type=float, default=90.0)
    p.add_argument("--flank-match-cov",  type=float, default=0.50,
                   help="min fraction of flank that must align to count as 'matching'")
    p.add_argument("--flank-match-id",   type=float, default=90.0,
                   help="min identity for flank match")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--require-same-strand", action="store_true",
                   help="require the two flanking ISes to be in the SAME orientation in ref "
                        "(canonical composite-transposon substrate for HR mobilization)")
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


def best_forward_block(ref_path, query_path, work, threads):
    """Return best forward-strand block: (tname, ts, te, block_bp, ident) or None.
    Used for flank matching where we only care about the single best alignment."""
    paf = os.path.join(work, "aln.paf")
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
    """Run minimap2 and return CLUSTERS of nearby forward-strand blocks.

    Each cluster represents one putative location where the query is present
    in the target. Within a cluster, blocks are merged (overlap-corrected) to
    compute total query coverage and weighted identity.

    Returns: list of dicts with keys:
        tname, ts, te (cluster span in tgt),
        q_blocks (union of query intervals covered),
        total_q_cov_bp (unique query bp covered),
        cov (fraction of query covered),
        weighted_ident
    """
    paf = os.path.join(work, "aln.paf")
    subprocess.run(["minimap2", "-x", "asm10", "-c", "--eqx",
                    "-t", str(threads), ref_path, query_path, "-o", paf],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    blocks = []
    qlen = None
    with open(paf) as fh:
        for line in fh:
            c = line.split("\t")
            if len(c) < 12: continue
            if c[4] != "+": continue
            qs, qe = int(c[2]), int(c[3])
            qlen = int(c[1])
            tname = c[5]; ts, te = int(c[7]), int(c[8])
            matches, block_len = int(c[9]), int(c[10])
            blocks.append({"tname": tname, "ts": ts, "te": te,
                           "qs": qs, "qe": qe,
                           "matches": matches, "block": block_len})
    if not blocks or qlen is None:
        return []
    # group by tname then sort by ts and cluster within cluster_gap
    by_tname = {}
    for b in blocks:
        by_tname.setdefault(b["tname"], []).append(b)
    clusters = []
    for tname, blks in by_tname.items():
        blks.sort(key=lambda x: x["ts"])
        cluster_blks = [blks[0]]
        for b in blks[1:]:
            if b["ts"] - cluster_blks[-1]["te"] <= cluster_gap:
                cluster_blks.append(b)
            else:
                clusters.append((tname, cluster_blks))
                cluster_blks = [b]
        clusters.append((tname, cluster_blks))
    # compute per-cluster metrics
    result = []
    for tname, blks in clusters:
        ts = min(b["ts"] for b in blks)
        te = max(b["te"] for b in blks)
        # union of query intervals (overlap-corrected coverage)
        q_iv = sorted([(b["qs"], b["qe"]) for b in blks])
        merged_q = [list(q_iv[0])]
        for s, e in q_iv[1:]:
            if s <= merged_q[-1][1]:
                merged_q[-1][1] = max(merged_q[-1][1], e)
            else:
                merged_q.append([s, e])
        total_q = sum(e - s for s, e in merged_q)
        total_matches = sum(b["matches"] for b in blks)
        total_blk     = sum(b["block"]   for b in blks)
        wid = total_matches / total_blk * 100 if total_blk > 0 else 0
        result.append({
            "tname": tname, "ts": ts, "te": te,
            "total_q_cov_bp": total_q,
            "cov": total_q / qlen if qlen else 0,
            "weighted_ident": wid,
            "n_blocks": len(blks),
        })
    return result


def find_composite_candidates(ref_is_index, ref_fa, args):
    """For each ref contig, yield candidate composite transposons:
       (contig, is_a, is_b, segment_seq, left_flank_seq, right_flank_seq)
    """
    candidates = []
    for contig, hits in ref_is_index.items():
        merged = merge_is(hits)
        try:
            clen = ref_fa.get_reference_length(contig)
        except (KeyError, ValueError):
            continue
        for i in range(len(merged)):
            # adjacent pairs only: with every pair, IS at A, B, C gave the
            # overlapping candidates (A,B), (B,C) and (A,C) and one relocated
            # region up to three MOVED rows
            for j in range(i + 1, min(i + 2, len(merged))):
                is_a = merged[i]; is_b = merged[j]
                # ORF coordinates are 1-based inclusive; segment = the bases
                # strictly between the two IS (0-based half-open)
                seg_start = is_a[1]; seg_end = is_b[0] - 1
                seg_len = seg_end - seg_start
                if seg_len < args.min_segment_bp or seg_len > args.max_segment_bp:
                    continue
                if args.require_same_strand and is_a[2] != is_b[2]:
                    continue
                lf_s = max(0, is_a[0] - 1 - args.flank); lf_e = is_a[0] - 1
                rf_s = is_b[1]; rf_e = min(clen, is_b[1] + args.flank)
                seg = ref_fa.fetch(contig, seg_start, seg_end)
                lf = ref_fa.fetch(contig, lf_s, lf_e) if lf_e > lf_s else ""
                rf = ref_fa.fetch(contig, rf_s, rf_e) if rf_e > rf_s else ""
                if not seg or not lf or not rf:
                    continue
                candidates.append({
                    "contig": contig, "is_a": is_a, "is_b": is_b,
                    "segment_pos": (seg_start, seg_end),
                    "segment_seq": seg, "left_flank": lf, "right_flank": rf,
                })
    return candidates


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)
    is_hits = load_is_hits(args.is_hits)

    work_root = tempfile.mkdtemp(prefix="movecheck_")
    ref_fa_path = find_assembly_fa(args.src_dir, args.ref_asm, work_root)
    if not ref_fa_path:
        raise SystemExit(f"ref {args.ref_asm} FASTA not found in {args.src_dir}")
    ref_fa = pysam.FastaFile(ref_fa_path)

    candidates = find_composite_candidates(is_hits.get(args.ref_asm, {}), ref_fa, args)
    print(f"ref {args.ref_asm}: {len(candidates)} composite-transposon candidates "
          f"(segment {args.min_segment_bp}-{args.max_segment_bp} bp"
          f"{', same-orient' if args.require_same_strand else ''})", flush=True)
    if not candidates:
        print("no candidates — exiting"); return

    tgts = [t.strip() for t in open(args.tgt_list) if t.strip()]
    # skip the reference AND its GCA/GCF twin (the twin counted as another genome)
    tgts = [t for t in tgts if not lm.same_genome(t, args.ref_asm)]
    print(f"testing {len(candidates)} candidates against {len(tgts)} target assemblies", flush=True)

    out_path = os.path.join(args.out, "per_candidate_movement.tsv")
    n = defaultdict(int)
    n_multi = 0   # candidates with >= 2 clusters in one tgt

    def tgt_flanks(fa, cl, la, lb):
        """Target sequence outside the segment, long enough to reach past the
        bracketing IS copies: the reference flanks lie OUTSIDE IS_a / IS_b,
        so the target window must be flank + IS long (query orientation)."""
        tlen = fa.get_reference_length(cl["tname"])
        a, b = lm.projected_ends(cl)
        a, b = max(0, min(tlen, a)), max(0, min(tlen, b))
        if cl["strand"] == "+":
            return (fa.fetch(cl["tname"], max(0, a - args.flank - la), a).upper(),
                    fa.fetch(cl["tname"], b, min(tlen, b + args.flank + lb)).upper())
        return (revcomp(fa.fetch(cl["tname"], b, min(tlen, b + args.flank + la)).upper()),
                revcomp(fa.fetch(cl["tname"], max(0, a - args.flank - lb), a).upper()))

    with open(out_path, "w") as out:
        out.write("ref_asm\tref_contig\tref_is_a_start\tref_is_b_start\tsegment_bp\t"
                  "tgt_asm\tcluster_idx\tn_clusters_in_tgt\tB_cov\tB_ident\tn_aln_blocks\t"
                  "tgt_contig\ttgt_pos\t"
                  "left_flank_match\tright_flank_match\tverdict\ttgt_strand\n")
        for ci, cand in enumerate(candidates):
            work = tempfile.mkdtemp(prefix=f"c{ci}_", dir=work_root)
            seg_path = os.path.join(work, "segment.fa")
            write_fasta("SEG", cand["segment_seq"], seg_path)
            seg_len = len(cand["segment_seq"])
            la = cand["is_a"][1] - cand["is_a"][0] + 1
            lb = cand["is_b"][1] - cand["is_b"][0] + 1
            head = (f"{args.ref_asm}\t{cand['contig']}\t"
                    f"{cand['is_a'][0]}\t{cand['is_b'][0]}\t{seg_len}\t")
            print(f"\n  cand {ci+1}/{len(candidates)}: "
                  f"{cand['contig']}:{cand['segment_pos'][0]}-{cand['segment_pos'][1]} "
                  f"({seg_len} bp)", flush=True)

            for ti, tgt in enumerate(tgts):
                tgt_path = find_assembly_fa(args.src_dir, tgt, work_root)
                if not tgt_path:
                    # an explicit row: these used to vanish (only a counter)
                    out.write(f"{head}{tgt}\t0\t0\t0.000\t0.0\t0\t\t\tNA\tNA\tNO_FASTA\t\n")
                    n["NO_FASTA"] += 1; continue
                paf = os.path.join(work, f"t{ti}.paf")
                # a minimap2 failure raises: caught as "failed" it produced no row
                subprocess.run(["minimap2", "-x", "asm10", "-c", "--eqx",
                                "-t", str(args.threads), tgt_path, seg_path, "-o", paf],
                               check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                # both strands: a target contig assembled in the other
                # orientation used to be reported ABSENT
                clusters = lm.parse_clusters(paf, {"SEG": seg_len}, 10000)
                os.remove(paf)
                passing = [c for c in clusters
                           if c["cov"] >= args.min_coverage
                           and c["weighted_ident"] >= args.min_identity]
                if not passing:
                    out.write(f"{head}{tgt}\t0\t0\t0.000\t0.0\t0\t\t\tNA\tNA\tABSENT\t\n")
                    n["ABSENT"] += 1; continue
                if len(passing) > 1:
                    n_multi += 1
                tgt_fa = pysam.FastaFile(tgt_path)
                for ki, cl in enumerate(passing):
                    # flanks are compared as SEQUENCES per cluster. They were
                    # files named by cluster index and reused on existence: a
                    # target with an empty flank was judged on the previous
                    # target's flank file.
                    tlf, trf = tgt_flanks(tgt_fa, cl, la, lb)
                    L = lm.flank_match(cand["left_flank"].upper(), tlf, args.flank_match_cov,
                                       args.flank_match_id, work, args.threads)
                    R = lm.flank_match(cand["right_flank"].upper(), trf, args.flank_match_cov,
                                       args.flank_match_id, work, args.threads)
                    verdict = lm.site_verdict(L[0], R[0])
                    n[verdict] += 1
                    out.write(f"{head}{tgt}\t{ki}\t{len(passing)}\t"
                              f"{cl['cov']:.3f}\t{cl['weighted_ident']:.1f}\t{cl['n_blocks']}\t"
                              f"{cl['tname']}\t{cl['ts']}\t"
                              f"{lm.fmt(*L)}\t{lm.fmt(*R)}\t{verdict}\t{cl['strand']}\n")
                tgt_fa.close()
            shutil.rmtree(work, ignore_errors=True)

    ref_fa.close()
    shutil.rmtree(work_root, ignore_errors=True)
    print(f"\nDONE. {dict(n)}  multi-cluster tgts={n_multi}")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
