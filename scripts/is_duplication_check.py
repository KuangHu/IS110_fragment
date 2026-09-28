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
import argparse, csv, gzip, json, os, re, shutil, subprocess, tempfile
from collections import defaultdict
import pysam


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
    loci = load_loci(args.records)
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
        anchors = parse_anchors(r.get("details", ""))
        if len(anchors) < 2:
            continue
        # Group by contig and keep only contigs with >= 2 distinct positions
        bycontig = defaultdict(list)
        for a in anchors:
            bycontig[a[0]].append(a)
        anchors_dup = []
        for c, hs in bycontig.items():
            positions = list({(h[1], h[2], h[3]) for h in hs})
            # need 2+ distinct same-contig positions separated by >= min-copy-sep
            sorted_pos = sorted(positions, key=lambda x: x[0])
            distinct = [sorted_pos[0]]
            for p in sorted_pos[1:]:
                if abs(p[0] - distinct[-1][0]) >= args.min_copy_sep:
                    distinct.append(p)
            if len(distinct) >= 2:
                for s, e, strand in distinct:
                    anchors_dup.append((c, s, e, strand))
        if len(anchors_dup) < 2:
            continue
        r["_anchors"] = anchors_dup
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

    n_pass = 0
    with open(out_path, "w") as out:
        out.write("ref_id\ttgt_asm\tbest_side\tn_anchor_positions\t"
                  "n_is_flanked\tcopy_separation\tn_fwd_blocks\t"
                  "best_block_bp\tbest_block_ident\tflanked_positions\tpasses\n")
        for i, r in enumerate(rows):
            if i % 50 == 0:
                print(f"  {i}/{len(rows)}", flush=True)
            ref_id  = r["ref_id"]
            tgt_asm = r["assembly"]
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

            anchors = r["_anchors"]
            target_is_idx = is_hits.get(tgt_asm, {})

            # Test (b): IS-flanking at each anchor hit position
            flanked = []
            for (c, s, e, strand) in anchors:
                mid = (s + e) // 2
                nis = nearest_is(target_is_idx, c, mid, args.is_tol)
                if nis is not None:
                    flanked.append((c, s, e, nis[0], nis[1], nis[2]))
            n_flanked = len(flanked)

            # Position separation
            positions = sorted({(a[0], a[1]) for a in anchors})
            sep = positions[-1][1] - positions[0][1] if len(positions) >= 2 else 0

            # Test (c): forward alignment blocks of MIDDLE in target region
            work = tempfile.mkdtemp(prefix="dup_", dir=work_root)
            if r["_best_side"] == "down":
                ms = is_end; me = min(clen, is_end + args.middle_window)
            else:
                ms = max(0, is_start - args.middle_window); me = is_start
            mid_seq = ref_fa.fetch(ref_contig, ms, me) if me > ms else ""
            ref_mid_path = os.path.join(work, "ref_middle.fa")
            if mid_seq:
                write_fasta("REF_MID", mid_seq, ref_mid_path)
            tgt_region_path = os.path.join(work, "tgt.fa")
            tgt_bp = write_target_region(tgt_fa, anchors,
                                         args.tgt_context, tgt_region_path)
            blocks = []
            if mid_seq and tgt_bp > 0:
                try:
                    blocks = find_forward_blocks(
                        ref_mid_path, tgt_region_path, work, args.threads,
                        min_block_bp=args.min_block_bp, min_identity=args.min_identity,
                    )
                except Exception:
                    blocks = []
            n_fwd = len(blocks)
            best_block_bp = max((b["block_bp"] for b in blocks), default=0)
            best_ident    = max((b["identity"]    for b in blocks), default=0.0)

            # PASS criterion: IS-flanked anchor positions AND (if requested)
            # enough forward-strand alignment blocks of MIDDLE in target
            passes = (n_flanked >= args.min_flanked
                      and len(positions) >= 2
                      and sep >= args.min_copy_sep
                      and n_fwd >= args.min_fwd_blocks)
            if passes: n_pass += 1
            flank_str = ";".join(f"{c}:{s}-{e}~IS:{iss}-{ise}@{d}bp"
                                 for (c, s, e, iss, ise, d) in flanked[:3])
            out.write(f"{ref_id}\t{tgt_asm}\t{r['_best_side']}\t{len(positions)}\t"
                      f"{n_flanked}\t{sep}\t{n_fwd}\t{best_block_bp}\t{best_ident:.1f}\t"
                      f"{flank_str}\t{passes}\n")
            shutil.rmtree(work, ignore_errors=True)

    for fa in fa_cache.values():
        if fa: fa.close()
    shutil.rmtree(work_root, ignore_errors=True)
    print(f"\nDONE: {n_pass}/{len(rows)} pass IS-mediated duplication test "
          f"(>= {args.min_flanked} anchor hits IS-flanked within {args.is_tol} bp, "
          f"sep >= {args.min_copy_sep} bp)")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
