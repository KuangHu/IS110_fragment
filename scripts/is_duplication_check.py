#!/usr/bin/env python3
"""Strict IS-mediated duplication check (v2: IS-flanking + forward-block evidence).

Biological model — IS-mediated duplication via HR between two IS copies:

    REF:     --upstream--[IS_a]-- segment --[IS_b]--downstream--
    TARGET:  --upstream--[IS_a]-- segment --[IS_?]-- segment --[IS_b]--downstream--

The signature in target is: the segment between IS_a and IS_b appears TWICE,
separated by an IS copy, all in direct (same-strand) orientation.

Detection signal:
  (a) the same ref-side anchor (1 kb) hits target at 2+ DISTINCT positions on
      the same target contig (= find_rearrangements `has_duplication`)
  (b) each anchor hit position in target is FLANKED by a confirmed IS110 copy
      (using `is_hits.tsv` for the target assembly) within --is-tol bp
  (c) bonus: a forward-strand alignment block of ref MIDDLE (the side facing
      the duplication) covers each anchor hit location at >= --min-identity

PASS criteria (configurable):
  - >= 2 distinct anchor hit positions (sep >= --min-copy-sep)
  - >= --min-flanked-positions of those positions have an IS within --is-tol bp
  - (optional) at least --min-fwd-blocks forward-strand blocks of ref MIDDLE
    appearing in target at >= --min-identity identity and >= --min-block-bp bp
"""
import argparse, csv, gzip, json, os, re, shutil, subprocess, tempfile, sys
from collections import defaultdict
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib_rearr as lr  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--candidates", required=True)
    p.add_argument("--records",    required=True)
    p.add_argument("--src-dir",    required=True)
    p.add_argument("--is-hits",    required=True,
                   help="is_hits.tsv (per-species IS detections in the DB)")
    p.add_argument("--out",        required=True)
    p.add_argument("--middle-window", type=int, default=30000)
    p.add_argument("--tgt-context",   type=int, default=30000)
    p.add_argument("--is-tol",         type=int, default=2000,
                   help="anchor position must be within this many bp of an IS hit in target")
    p.add_argument("--min-copy-sep",   type=int, default=5000,
                   help="anchor hits must be at least this far apart on target contig")
    p.add_argument("--min-flanked",    type=int, default=2,
                   help="minimum number of anchor hits that must be IS-flanked")
    p.add_argument("--min-fwd-blocks", type=int, default=2,
                   help="minimum forward-strand alignment blocks of ref MIDDLE in target")
    p.add_argument("--min-block-bp",   type=int, default=1000,
                   help="forward block min size (bp)")
    p.add_argument("--min-identity",   type=float, default=95.0)
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


def load_is_hits(is_hits_path):
    """Index: assembly -> contig(local-name) -> [(start, end, strand), ...]."""
    idx = defaultdict(lambda: defaultdict(list))
    with open(is_hits_path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            asm = r["assembly"]
            contig = r["contig"]
            local = contig.split("|", 1)[1] if "|" in contig else contig
            idx[asm][local].append((int(r["tnp_start"]), int(r["tnp_end"]),
                                    r.get("tnp_strand", "+")))
    return idx


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
def parse_anchors(details):
    return [(m.group(1), int(m.group(2)), int(m.group(3)), m.group(4))
            for m in ANCHOR_RE.finditer(details)]


def write_target_region(tgt_fa, anchors, window, out_path):
    if not anchors:
        return 0
    spans = {}
    for (c, s, e, _strand) in anchors:
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
            sx = max(0, s - window); ex = min(clen, e + window)
            seq = tgt_fa.fetch(c, sx, ex)
            if seq:
                fh.write(f">TGT_{c}_{sx}_{ex}\n")
                for k in range(0, len(seq), 80):
                    fh.write(seq[k:k+80] + "\n")
                total += len(seq)
    return total


def find_forward_blocks(ref_fa_path, tgt_fa_path, work, threads,
                        min_block_bp, min_identity):
    paf = os.path.join(work, "aln.paf")
    subprocess.run(["minimap2", "-x", "asm10", "-c", "--eqx",
                    "-t", str(threads), ref_fa_path, tgt_fa_path, "-o", paf],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    blocks = []
    with open(paf) as f:
        for line in f:
            c = line.split("\t")
            if len(c) < 12: continue
            if c[4] != "+": continue
            qs, qe = int(c[2]), int(c[3])
            matches, block_len = int(c[9]), int(c[10])
            if block_len < min_block_bp: continue
            ident = matches / block_len * 100 if block_len > 0 else 0
            if ident < min_identity: continue
            # parse q_name = "TGT_<contig>_<sx>_<ex>" so we can map q positions
            # back to absolute target-contig coordinates
            qname = c[0]
            blocks.append({
                "qname": qname, "q_start": qs, "q_end": qe,
                "block_bp": block_len, "identity": ident,
            })
    return blocks


def parse_tgt_qname(qname):
    """TGT_<contig-with-_/.>_<sx>_<ex>  → (contig, sx, ex)"""
    if not qname.startswith("TGT_"):
        return None, None, None
    rest = qname[4:]
    # the last two _-separated tokens are integers
    parts = rest.rsplit("_", 2)
    if len(parts) < 3:
        return None, None, None
    contig, sx, ex = parts
    try:
        return contig, int(sx), int(ex)
    except ValueError:
        return None, None, None


def nearest_is(is_hits_for_asm, contig, pos, tol):
    """Return (is_start, is_end, distance) of the nearest IS hit within `tol` of pos
    on target contig. None if none found within tolerance."""
    hits = is_hits_for_asm.get(contig, [])
    if not hits:
        return None
    best = None
    for (s, e, strand) in hits:
        # distance from pos to interval [s, e]
        if s <= pos <= e:
            d = 0
        else:
            d = min(abs(pos - s), abs(pos - e))
        if d <= tol and (best is None or d < best[2]):
            best = (s, e, d)
    return best


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)
    loci = lr.load_loci(args.records)
    is_hits = load_is_hits(args.is_hits)
    print(f"loaded is_hits: {sum(sum(len(v) for v in d.values()) for d in is_hits.values())} hits "
          f"across {len(is_hits)} assemblies", flush=True)

    rows = []
    for r in csv.DictReader(open(args.candidates), delimiter="\t"):
        if r.get("category") != "duplication":
            continue
        side = (r.get("best_side") or "").strip()
        if side not in ("up", "down"):
            continue
        d_ = lr.parse_details(r.get("details", ""))
        anchors = d_["up"] + d_["down"]
        if len(anchors) < 2:
            continue
        # Group by contig and keep only contigs with >= 2 distinct positions
        bycontig = defaultdict(list)
        for a in anchors:
            bycontig[a[0]].append(a)
        anchors_dup = []
        for c, hs in sorted(bycontig.items()):
            # thin positions PER STRAND and with a total order. A set sorted
            # by start alone kept an arbitrary strand when an anchor hit both
            # strands at one position, so the verdict changed with
            # PYTHONHASHSEED (True / False / NA for the same input).
            for st in "+-":
                sorted_pos = sorted({(h[1], h[2]) for h in hs if h[3] == st})
                if not sorted_pos:
                    continue
                distinct = [sorted_pos[0]]
                for p in sorted_pos[1:]:
                    if abs(p[0] - distinct[-1][0]) >= args.min_copy_sep:
                        distinct.append(p)
                if len(distinct) >= 2:
                    for s, e in distinct:
                        anchors_dup.append((c, s, e, st))
        # no two same-strand hits far enough apart: kept, and written below
        # as a not-assessable row (it used to vanish from the output)
        r["_anchors"] = anchors_dup if len(anchors_dup) >= 2 else []
        r["_best_side"] = side
        rows.append(r)
    if args.sample > 0 and len(rows) > args.sample:
        import random
        random.seed(13)
        rows = random.sample(rows, args.sample)
    print(f"Checking {len(rows)} duplication candidates "
          f"(min {args.min_copy_sep} bp separation on same contig) ...", flush=True)

    out_path = os.path.join(args.out, "per_candidate_duplication.tsv")
    work_root = tempfile.mkdtemp(prefix="dupcheck_")
    fa_cache = {}
    def open_asm(asm):
        if asm in fa_cache: return fa_cache[asm]
        path = find_assembly_fa(args.src_dir, asm, work_root)
        fa_cache[asm] = pysam.FastaFile(path) if path else None
        return fa_cache[asm]

    n_pass = n_assessed = 0
    with open(out_path, "w") as out:
        out.write("ref_id\ttgt_asm\tbest_side\tn_anchor_positions\t"
                  "n_is_flanked\tcopy_separation\tn_fwd_blocks\t"
                  "best_block_bp\tbest_block_ident\tflanked_positions\tpasses\t"
                  "dup_bp\ttgt_contig\tstatus\n")

        def skip(ref_id, tgt_asm, side, why):
            out.write(f"{ref_id}\t{tgt_asm}\t{side}\t0\t0\t0\t0\t0\t0.0\t\tNA\t0\t\t{why}\n")

        for i, r in enumerate(rows):
            if i % 50 == 0:
                print(f"  {i}/{len(rows)}", flush=True)
            ref_id, tgt_asm, side = r["ref_id"], r["assembly"], r["_best_side"]
            locus = loci.get(ref_id)
            if not locus:
                skip(ref_id, tgt_asm, side, "ref_id not in records"); continue
            ref_contig, is_s0, is_e, tnp_strand = locus
            ref_fa = open_asm(ref_id.split("|")[0])
            tgt_fa = open_asm(tgt_asm)
            if not ref_fa or not tgt_fa:
                skip(ref_id, tgt_asm, side, "assembly FASTA not found"); continue
            try:
                # MIDDLE = the flank on the anchor's side, in TRANSPOSASE
                # orientation. best_side is up/down relative to the transposase;
                # it was applied as genomic left/right, so every minus-strand
                # element was tested with the flank on the wrong side.
                up_seq, down_seq = lr.ref_proxies(ref_fa, ref_contig, is_s0, is_e,
                                                  tnp_strand, args.middle_window)
            except (KeyError, ValueError):
                skip(ref_id, tgt_asm, side, "ref contig not in FASTA"); continue
            mid_seq = down_seq if side == "down" else up_seq
            if len(mid_seq) < args.min_block_bp:
                skip(ref_id, tgt_asm, side, "reference flank too short (contig end)"); continue

            # How many copies of this flank does the REFERENCE itself have
            # near the element? The target must have MORE: a genome that
            # already carries the segment twice (A-E-D-X-E-D-B, unchanged in
            # the target) passed as a new duplication.
            work0 = tempfile.mkdtemp(prefix="dupref_", dir=work_root)
            clen_r = ref_fa.get_reference_length(ref_contig)
            r_lo = max(0, is_s0 - args.middle_window - args.tgt_context)
            r_hi = min(clen_r, is_e + args.middle_window + args.tgt_context)
            ref_win = ref_fa.fetch(ref_contig, r_lo, r_hi).upper()
            ref_copies = max(
                lr.copy_depth([b for b in lr.align(mid_seq, ref_win, work0, args.threads)],
                              args.min_block_bp, args.min_identity, strand=st_)[1]
                for st_ in "+-")
            shutil.rmtree(work0, ignore_errors=True)
            target_is_idx = is_hits.get(tgt_asm, {})
            # Evaluate ONE contig and ONE strand at a time. Positions on
            # different contigs were pooled: the 'separation' was a difference
            # of coordinates on two contigs, one flanked hit on each of two
            # contigs satisfied --min-flanked 2, and a second copy in INVERTED
            # orientation passed a test that requires direct repeats.
            groups = defaultdict(list)
            for (c, s_, e_, st) in r["_anchors"]:
                groups[(c, st)].append((s_, e_))
            best = None
            work = tempfile.mkdtemp(prefix="dup_", dir=work_root)
            for (c, st), pos in sorted(groups.items()):
                pos = sorted(set(pos))
                if len(pos) < 2:
                    continue
                sep = pos[-1][0] - pos[0][0]
                flanked = []
                for (s_, e_) in pos:
                    nis = nearest_is(target_is_idx, c, (s_ + e_) // 2, args.is_tol)
                    if nis is not None:
                        flanked.append((c, s_, e_, nis[0], nis[1], nis[2]))
                tseq, toff = lr.fetch_window(tgt_fa, c, pos[0][0], pos[-1][1], args.tgt_context)
                if not tseq:
                    continue
                # the segment must be present at >= 2 separate places, on the
                # anchors' strand; blocks that are just an IS copy do not count
                blocks = []
                for b in lr.align(mid_seq, tseq, work, args.threads):
                    ts_, te_ = b["ts"] + toff, b["te"] + toff
                    is_bp = lr.union_bp([(max(ts_, h[0] - 300), min(te_, h[1] + 300))
                                         for h in target_is_idx.get(c, [])
                                         if min(te_, h[1] + 300) > max(ts_, h[0] - 300)])
                    if is_bp < 0.8 * (te_ - ts_):
                        blocks.append(b)
                dup_bp, copies, places = lr.copy_depth(blocks, args.min_block_bp,
                                                       args.min_identity, strand=st)
                cand = {"contig": c, "n_pos": len(pos), "sep": sep, "flanked": flanked,
                        "dup_bp": dup_bp, "copies": copies,
                        "best_bp": max((p["qe"] - p["qs"] for p in places), default=0),
                        "best_id": max((p["ident"] for p in places), default=0.0)}
                cand["ref_copies"] = ref_copies
                cand["passes"] = (len(flanked) >= args.min_flanked
                                  and sep >= args.min_copy_sep
                                  and copies >= args.min_fwd_blocks
                                  and copies > ref_copies          # a GAIN
                                  and dup_bp >= args.min_block_bp)
                if best is None or (cand["passes"], cand["dup_bp"]) > (best["passes"], best["dup_bp"]):
                    best = cand
            shutil.rmtree(work, ignore_errors=True)
            if best is None:
                skip(ref_id, tgt_asm, side, "no contig with two same-strand anchor hits"); continue
            n_assessed += 1
            if best["passes"]: n_pass += 1
            flank_str = ";".join(f"{c}:{s_}-{e_}~IS:{iss}-{ise}@{d}bp"
                                 for (c, s_, e_, iss, ise, d) in best["flanked"][:3])
            out.write(f"{ref_id}\t{tgt_asm}\t{side}\t{best['n_pos']}\t"
                      f"{len(best['flanked'])}\t{best['sep']}\t{best['copies']}\t"
                      f"{best['best_bp']}\t{best['best_id']:.1f}\t{flank_str}\t"
                      f"{best['passes']}\t{best['dup_bp']}\t{best['contig']}\tassessed\n")

    for fa in fa_cache.values():
        if fa: fa.close()
    shutil.rmtree(work_root, ignore_errors=True)
    print(f"\nDONE: {n_pass}/{n_assessed} assessed candidates pass the IS-mediated "
          f"duplication test (>= {args.min_flanked} anchor hits IS-flanked within "
          f"{args.is_tol} bp, sep >= {args.min_copy_sep} bp, >= {args.min_block_bp} bp of the "
          f"flank present at >= {args.min_fwd_blocks} places on one contig and strand); "
          f"{len(rows) - n_assessed} not assessable")
    print("  NOTE: rows are per (reference element, target genome). The same target "
          "duplication is reported from every reference that carries the element.")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
