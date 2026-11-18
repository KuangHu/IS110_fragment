#!/usr/bin/env python3
"""Validate IS-mediated rearrangement calls with TWO gold-standard tools.

For each candidate event from find_rearrangements.py (rearrangements_examples.tsv),
align the reference assembly against the target assembly with two independent
whole-genome comparison tools and ask whether each one independently calls a
rearrangement of the matching class:

  Tool 1 — MUMmer4 dnadiff   (nucmer-based; .report breakpoint counts + .qdiff)
  Tool 2 — syri              (minimap2-based; syri.out structural-variant table)

A candidate is CONFIRMED only if BOTH tools report the matching event class
between the two genomes (genome-level), and IS-mediated if a tool breakpoint
falls within --bp-window of the candidate's anchor contigs (locus-level).

Input assemblies live as <src-dir>/<assembly>.fna(.gz). ref assembly is parsed
from ref_id (the part before '|'); target assembly is the `assembly` column.

Usage:
  validate_rearrangements_gold.py --examples rearrangements_examples.tsv \\
      --src-dir /path/to/by_species/<species>/ --out validate/ \\
      --n-per-cat 10 --threads 8
"""
# VERDICTS (changed 2026-11): CONFIRMED_BOTH needs dnadiff's class count AND
#   inversion      Mauve: some-but-not-all LCBs inverted, >= --min-inverted-bp
#                  of them adjacent to the element, target contigs first put
#                  in reference orientation
#   duplication    >= --min-dup-bp of the reference window (outside the IS)
#                  aligned at two places in the target window
#   translocation  both reference flanks follow their anchors (>= 80%) onto
#                  the two contigs; anchors not at a contig end
# CONFIRMED_DNADIFF is no longer produced.
import argparse, gzip, json, os, re, shutil, subprocess, sys, tempfile
from collections import defaultdict

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib_rearr as lr  # noqa: E402
from lib_alleles import revcomp  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--examples", required=True, help="rearrangements_examples.tsv")
    p.add_argument("--src-dir", required=True,
                   help="dir with <assembly>.fna.gz for this species")
    p.add_argument("--records", required=True,
                   help="records.json (gives the reference IS locus per ref_id)")
    p.add_argument("--out", required=True)
    p.add_argument("--n-per-cat", type=int, default=10,
                   help="candidates to validate per category (default 10)")
    p.add_argument("--locus-window", type=int, default=30000,
                   help="bp of context to extract on each side of the locus "
                        "(default 30000). Restricting alignment to the locus makes "
                        "BOTH tools' calls inherently locus-level.")
    p.add_argument("--bp-window", type=int, default=5000,
                   help="a translocation anchor within this many bp of its target "
                        "contig end is treated as a contig-break artifact (default 5000)")
    p.add_argument("--min-inverted-bp", type=int, default=5000,
                   help="inversion: inverted LCBs adjacent to the element must "
                        "total at least this many bp (default 5000)")
    p.add_argument("--min-dup-bp", type=int, default=1000,
                   help="duplication: reference bp (outside the element) present at "
                        ">= 2 places in the target window (default 1000)")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--mauve", default="progressiveMauve",
                   help="progressiveMauve binary")
    p.add_argument("--keep-tmp", action="store_true")
    return p.parse_args()


def load_ref_loci(records_path):
    """ref_id -> (local_contig, start0, end, transposase strand)."""
    return lr.load_loci(records_path)


def extract_region(fa_obj, contig, start, end, window, label, out_handle, rc=False):
    """Write [start-window, end+window) of `contig` from an open pysam.FastaFile
    (reverse-complemented when rc)."""
    try:
        clen = fa_obj.get_reference_length(contig)
    except (KeyError, ValueError):
        return False
    s = max(0, start - window)
    e = min(clen, end + window)
    seq = fa_obj.fetch(contig, s, e)
    if not seq:
        return False
    if rc:
        seq = revcomp(seq.upper())
    out_handle.write(f">{label}\n")
    for i in range(0, len(seq), 80):
        out_handle.write(seq[i:i+80] + "\n")
    return True


def find_assembly_fa(src_dir, asm, work):
    """Locate <asm>.fna(.gz) in src_dir; return path to a plain .fa (decompressed)."""
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


# ----------------------------- tool 1: dnadiff ------------------------------
def run_dnadiff(ref_fa, qry_fa, prefix):
    """Run dnadiff; return dict of breakpoint-class counts from the .report."""
    subprocess.run(["dnadiff", "-p", prefix, ref_fa, qry_fa],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    report = prefix + ".report"
    counts = {"Inversions": 0, "Relocations": 0, "Translocations": 0,
              "Insertions": 0, "TandemIns": 0}
    qry = dict.fromkeys(counts, 0)
    if os.path.exists(report):
        with open(report) as f:
            for line in f:
                parts = line.split()
                if not parts:
                    continue
                key = parts[0]
                if key in counts and len(parts) >= 2:
                    # take the REF-column value (first number)
                    nums = [p for p in parts[1:] if re.match(r"^\d+", p)]
                    if nums:
                        counts[key] = int(re.match(r"^\d+", nums[0]).group())
                    if len(nums) > 1:
                        qry[key] = int(re.match(r"^\d+", nums[1]).group())
    # A duplication in the TARGET is extra sequence in the query: dnadiff
    # lists it in the QRY column (Insertions / TandemIns). Only the REF
    # column was read, so the class that "confirmed" a duplication was
    # extra sequence in the REFERENCE -- i.e. a deletion in the target.
    counts["QryInsertions"] = qry["Insertions"] + qry["TandemIns"]
    return counts


def parse_dnadiff_breakpoints(prefix):
    """Parse .qdiff: list of (qry_contig, qstart, qend, gap_type)."""
    bps = []
    qdiff = prefix + ".qdiff"
    if os.path.exists(qdiff):
        with open(qdiff) as f:
            for line in f:
                c = line.split("\t")
                if len(c) < 5:
                    continue
                qry = c[0]; typ = c[1]
                try:
                    s, e = int(c[2]), int(c[3])
                except ValueError:
                    continue
                bps.append((qry, min(s, e), max(s, e), typ))
    return bps


# ----------------------------- tool 2: progressiveMauve ---------------------
def parse_xmfa_lcbs(xmfa):
    """Parse progressiveMauve XMFA -> list of LCBs; each = {seqnum: (start,end,strand)}."""
    lcbs = []
    cur = {}
    with open(xmfa) as f:
        for line in f:
            if line.startswith(">"):
                # > seqnum:start-end strand /path
                m = re.match(r">\s*(\d+):(\d+)-(\d+)\s+([+-])", line)
                if m:
                    seqn = int(m.group(1)); s = int(m.group(2))
                    e = int(m.group(3)); strand = m.group(4)
                    if s != 0 or e != 0:
                        cur[seqn] = (s, e, strand)
            elif line.startswith("="):
                if cur:
                    lcbs.append(cur)
                cur = {}
        if cur:
            lcbs.append(cur)
    return lcbs


def run_mauve(mauve_bin, ref_fa, qry_fa, work):
    """Run progressiveMauve; analyse LCB structure for inversions / reordering.
    Returns dict {INV: n_inverted_lcbs, REARR: n_order_breaks, n_lcbs, ok} or {ok:False}."""
    # absolute paths: Mauve runs with cwd=work, so a relative --out made
    # every call fail -- and a failure was scored as "no signal"
    work, ref_fa, qry_fa = map(os.path.abspath, (work, ref_fa, qry_fa))
    xmfa = os.path.join(work, "mauve.xmfa")
    try:
        subprocess.run([mauve_bin, f"--output={xmfa}", ref_fa, qry_fa],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                       cwd=work, timeout=1800)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return {"ok": False}
    if not os.path.exists(xmfa):
        return {"ok": False}
    lcbs = parse_xmfa_lcbs(xmfa)
    # keep LCBs present in BOTH genomes (1 = ref, 2 = qry)
    shared = [l for l in lcbs if 1 in l and 2 in l]
    if not shared:
        return {"ok": True, "INV": 0, "REARR": 0, "n_lcbs": 0}
    # inversion: opposite strand between the two genomes
    inv = [l for l in shared if l[1][2] != l[2][2]]
    n_inv = len(inv)
    # reordering: sort by genome-1 start; count discordances in genome-2 order
    order = sorted(shared, key=lambda l: l[1][0])
    g2_starts = [l[2][0] for l in order]
    n_breaks = sum(1 for i in range(1, len(g2_starts))
                   if g2_starts[i] < g2_starts[i - 1])
    return {"ok": True, "INV": n_inv, "REARR": n_breaks, "n_lcbs": len(shared),
            # reference intervals (1-based, in the region FASTA) of inverted LCBs
            "inv_ref": [(l[1][0], l[1][1]) for l in inv]}


# class mapping per candidate category
DNADIFF_CLASS = {"inversion_only": "Inversions",
                 "translocation_only": "Translocations",
                 "duplication": "QryInsertions"}  # extra sequence in the target
# progressiveMauve signal per category: which LCB-derived metric must be > 0
# Translocation: the LCB ORDER metric (REARR) is not evidence. With only the
# two anchor contigs in the query, a true A-IS-C / X-IS-B exchange is collinear
# with the reference A-IS-B, and the count that did come out depended on which
# contig was written first (i.e. on the strand of the reference element): the
# same genome pair gave CONFIRMED_ONE or CONFIRMED_BOTH by listing order.
# The second line of evidence is now the bilateral flank test.
MAUVE_SIGNAL = {"inversion_only": "INV",
                "translocation_only": None,
                "duplication": None}  # Mauve LCBs don't call duplication


def duplicated_ref_bp(ref_seq, tgt_seq, mask, work, threads, min_block=1000, min_ident=95.0):
    """bp of the reference region present at >= 2 non-overlapping places in the
    target region, outside `mask` (the IS element itself -- IS copies are
    expected to be repeated). This is the duplication; 'dnadiff reports an
    insertion somewhere in the window' confirmed a 300-bp deletion 13 kb away
    as a duplication."""
    blocks = [b for b in lr.align(ref_seq, tgt_seq, work, threads)
              if b["ident"] >= min_ident and b["qe"] - b["qs"] >= min_block]
    events = []
    for i, a in enumerate(blocks):
        for b in blocks[i + 1:]:
            if min(a["te"], b["te"]) - max(a["ts"], b["ts"]) > 0:
                continue                      # same place in the target
            s0, e0 = max(a["qs"], b["qs"]), min(a["qe"], b["qe"])
            if e0 - s0 > 0:
                events.append((s0, e0))
    covered = lr.union_bp(events)
    in_mask = lr.union_bp([(max(s0, mask[0]), min(e0, mask[1])) for s0, e0 in events
                           if min(e0, mask[1]) > max(s0, mask[0])])
    return covered - in_mask


def anchor_targets(details):
    """Extract [(contig, start, end)] from the details string.
    Contig names may contain '|' (e.g. ENA|CM003196|CM003196.1)."""
    out = []
    d = lr.parse_details(details)
    for c, s0, e0, st in d["up"] + d["down"]:
        out.append((c, s0, e0))
    return out


def locus_hit(bps, anchors, window):
    """True if any breakpoint falls on an anchor contig within window of its span."""
    for (contig, s, e) in anchors:
        for (bc, bs, be, _t) in bps:
            if bc == contig and bs - window <= e and be + window >= s:
                return True
    return False


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    # pick N per category
    by_cat = defaultdict(list)
    with open(args.examples) as f:
        header = f.readline()
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 6:
                continue
            by_cat[c[0]].append(c)
    candidates = []
    for cat, rows in by_cat.items():
        candidates.extend(rows[:args.n_per_cat])
    print(f"Validating {len(candidates)} candidates "
          f"({', '.join(f'{k}:{min(len(v),args.n_per_cat)}' for k,v in by_cat.items())})",
          file=sys.stderr, flush=True)

    ref_loci = load_ref_loci(args.records)
    print(f"Reference loci loaded: {len(ref_loci):,}", file=sys.stderr)
    W = args.locus_window

    results = []
    for i, c in enumerate(candidates):
        cat, ref_id, tgt_asm, up_c, down_c, details = c[:6]
        ref_asm = ref_id.split("|")[0]
        work = tempfile.mkdtemp(prefix="rgval_", dir=args.out)
        rec = {"category": cat, "ref_id": ref_id, "ref_asm": ref_asm,
               "tgt_asm": tgt_asm, "dnadiff": "NA", "mauve": "NA",
               "n_target_contigs": 0, "min_anchor_end_dist": -1,
               "mauve_metric": "", "verdict": "ERROR"}
        try:
            locus = ref_loci.get(ref_id)
            if not locus:
                rec["verdict"] = "NO_REF_LOCUS"
                results.append(rec); continue
            ref_asm_fa = find_assembly_fa(args.src_dir, ref_asm, work)
            qry_asm_fa = find_assembly_fa(args.src_dir, tgt_asm, work)
            if not ref_asm_fa or not qry_asm_fa:
                rec["verdict"] = "MISSING_ASM"
                results.append(rec); continue

            # ---- build locus-restricted region FASTAs ----
            ref_region = os.path.join(work, "ref_locus.fa")
            qry_region = os.path.join(work, "tgt_locus.fa")
            ref_fo = pysam.FastaFile(ref_asm_fa)
            qry_fo = pysam.FastaFile(qry_asm_fa)
            ref_contig, rs, re_, tnp_strand = locus
            with open(ref_region, "w") as fh:
                if not extract_region(ref_fo, ref_contig, rs, re_, W, "REF_LOCUS", fh):
                    rec["verdict"] = "REF_REGION_FAIL"; results.append(rec); continue
            ref_off = max(0, rs - W)                 # region offset in the contig
            det = lr.parse_details(details)
            # Stage 7 writes only the duplicated side for a duplication
            # (up=[] down=[a;b]); requiring both sides sent every real
            # duplication row to NO_ANCHOR_PAIR
            if (cat == "duplication" and not (det["up"] or det["down"])) or \
               (cat != "duplication" and not (det["up"] and det["down"])):
                rec["verdict"] = "NO_ANCHOR_PAIR"; results.append(rec); continue
            # Is each target contig stored in the same orientation as the
            # reference? Anchors are in transposase orientation, so a collinear
            # same-orientation contig has anchor strand == transposase strand.
            # Each contig is written in REFERENCE orientation (by its first
            # anchor), so "opposite strand" in the tools below means inverted
            # relative to the flank, not "this contig happens to be stored
            # reverse-complemented" (which scored as an inversion in Mauve).
            orient = {}
            for c_, s_, e_, st_ in det["up"] + det["down"]:
                orient.setdefault(c_, st_ != tnp_strand)
            # target: the anchor contigs (dedup), each ±W around the anchor span
            anchors = anchor_targets(details)
            seen_c = {}
            for (contig, s, e) in anchors:
                seen_c.setdefault(contig, [s, e])
                seen_c[contig][0] = min(seen_c[contig][0], s)
                seen_c[contig][1] = max(seen_c[contig][1], e)
            with open(qry_region, "w") as fh:
                ntc = 0
                for contig, (s, e) in seen_c.items():
                    if extract_region(qry_fo, contig, s, e, W, f"TGT_{contig}", fh,
                                      rc=orient.get(contig, False)):
                        ntc += 1
            rec["n_target_contigs"] = ntc
            if ntc == 0:
                rec["verdict"] = "TGT_REGION_FAIL"; results.append(rec); continue

            # anchor distance-to-contig-end: small => the anchor sits at a contig
            # break (the translocation 'signal' is then an assembly artifact, not
            # biology). Track the minimum across this candidate's anchors.
            min_end_dist = None
            for contig, (s, e) in seen_c.items():
                try:
                    clen = qry_fo.get_reference_length(contig)
                except (KeyError, ValueError):
                    continue
                d = min(s, clen - e)
                min_end_dist = d if min_end_dist is None else min(min_end_dist, d)
            rec["min_anchor_end_dist"] = min_end_dist if min_end_dist is not None else -1

            # Tool 1: dnadiff (MUMmer) on the locus regions
            dd_prefix = os.path.join(work, "dd")
            dd_counts = run_dnadiff(ref_region, qry_region, dd_prefix)
            dd_class = DNADIFF_CLASS.get(cat)
            dd_present = dd_counts.get(dd_class, 0) > 0 if dd_class else False
            rec["dnadiff"] = "YES" if dd_present else "no"

            # Tool 2: progressiveMauve (LCB structure) on the locus regions
            signal = MAUVE_SIGNAL.get(cat)
            second_ok = None                    # second line of evidence
            if signal == "INV":
                mv = run_mauve(args.mauve, ref_region, qry_region, work)
                if not mv.get("ok"):
                    rec["mauve"] = "FAILED"
                else:
                    # An inversion: SOME but not all LCBs inverted (targets are
                    # already in reference orientation), the inverted LCBs
                    # abut the element, and they are not a small unrelated
                    # inversion elsewhere in the window (any inversion of any
                    # size anywhere in +/-30 kb used to confirm).
                    is_a, is_b = rs - ref_off, re_ - ref_off
                    adj = [(a, b) for a, b in mv["inv_ref"]
                           if a <= is_b + 2000 + (re_ - rs) and b >= is_a - 2000 - (re_ - rs)]
                    inv_bp = lr.union_bp(adj)
                    partial = 0 < mv["INV"] < mv["n_lcbs"]
                    rec["mauve_metric"] = (f"INV:{mv['INV']}/{mv['n_lcbs']}LCB;"
                                           f"inv_bp_at_IS:{inv_bp}")
                    second_ok = partial and inv_bp >= args.min_inverted_bp
                    rec["mauve"] = "YES" if second_ok else "no"
            elif cat == "duplication":
                ref_seq = ref_fo.fetch(ref_contig, ref_off,
                                       min(ref_fo.get_reference_length(ref_contig), re_ + W)).upper()
                tgt_seq = "".join(
                    qry_fo.fetch(c_, max(0, s_ - W),
                                 min(qry_fo.get_reference_length(c_), e_ + W)).upper() + "N" * 200
                    for c_, (s_, e_) in seen_c.items())
                dup_bp = duplicated_ref_bp(ref_seq, tgt_seq, (rs - ref_off, re_ - ref_off),
                                           work, args.threads)
                # the same measure of the reference against ITSELF: sequence
                # already repeated in the reference is not a new duplication
                self_bp = duplicated_ref_bp(ref_seq, ref_seq, (rs - ref_off, re_ - ref_off),
                                            work, args.threads)
                gain = dup_bp - self_bp
                rec["mauve_metric"] = f"dup_bp_outside_IS:{dup_bp}|already_in_ref:{self_bp}"
                rec["mauve"] = "n/a"
                second_ok = gain >= args.min_dup_bp
            elif cat == "translocation_only":
                # bilateral flank test (see is_translocation_check.py)
                pair = next(((u, w) for u in det["up"] for w in det["down"] if u[0] != w[0]), None)
                rec["mauve"] = "n/a"
                if pair is not None:
                    A_seq, B_seq = lr.ref_proxies(ref_fo, ref_contig, rs, re_, tnp_strand, W)
                    fr = []
                    for prox, (c_, s_, e_, st_) in ((A_seq, pair[0]), (B_seq, pair[1])):
                        t_, _ = lr.fetch_window(qry_fo, c_, s_, e_, W)
                        fr.append(lr.strand_cover(lr.align(prox, t_, work, args.threads),
                                                  st_, 95.0) / max(1, len(prox))
                                  if len(prox) >= 1000 and t_ else 0.0)
                    rec["mauve_metric"] = f"flank_cov:{fr[0]:.2f}/{fr[1]:.2f}"
                    second_ok = min(fr) >= 0.8

            # verdict (duplication uses dnadiff only, since Mauve is n/a there)
            dd_ok = rec["dnadiff"] == "YES"
            mv_ok = bool(second_ok)
            # translocations across 2 contigs are unresolvable when an anchor sits
            # within --bp-window of a contig end (= contig-break artifact, not a
            # confirmable rearrangement in a draft target)
            # ... and the contig must also run PAST THE ELEMENT on the
            # element's side of the anchor (anchor 1 kb + element + window).
            # One genome split in the middle of a 16-kb element (A-E[:8k] /
            # E[8k:]-B) had deep anchors and was a confirmed translocation.
            through = True
            if cat == "translocation_only":
                need = 1000 + (re_ - rs) + args.bp_window
                for side, toward_plus in (("up", True), ("down", False)):
                    for c_, s_, e_, st_ in det[side][:1]:
                        try:
                            clen_ = qry_fo.get_reference_length(c_)
                        except (KeyError, ValueError):
                            continue
                        # the element lies downstream of the up anchor and
                        # upstream of the down anchor, in ANCHOR orientation
                        room = (clen_ - e_) if (st_ == "+") == toward_plus else s_
                        through = through and room >= need
            artifact = (cat == "translocation_only" and rec["n_target_contigs"] >= 2
                        and (0 <= rec.get("min_anchor_end_dist", -1) < args.bp_window
                             or not through))
            if artifact:
                rec["verdict"] = "UNRESOLVABLE_CONTIG_BREAK"
            elif rec["mauve"] == "FAILED":
                # a tool failure is an error, not half a confirmation
                rec["verdict"] = "ERROR:MAUVE"
            elif dd_ok and mv_ok:
                rec["verdict"] = "CONFIRMED_BOTH"
            elif dd_ok or mv_ok:
                rec["verdict"] = "CONFIRMED_ONE"
            else:
                rec["verdict"] = "REJECTED"
        except Exception as e:
            rec["verdict"] = f"ERROR:{type(e).__name__}"
        finally:
            if not args.keep_tmp:
                shutil.rmtree(work, ignore_errors=True)
        results.append(rec)
        if (i + 1) % 5 == 0:
            print(f"  {i+1}/{len(candidates)}", file=sys.stderr, flush=True)

    # write per-candidate + summary
    cols = ["category", "ref_id", "ref_asm", "tgt_asm", "n_target_contigs",
            "min_anchor_end_dist", "dnadiff", "mauve", "mauve_metric", "verdict"]
    with open(os.path.join(args.out, "validation.tsv"), "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in results:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    print("\n=== Validation summary (per category) ===", file=sys.stderr)
    for cat in by_cat:
        rs = [r for r in results if r["category"] == cat]
        n = len(rs)
        if not n:
            continue
        both = sum(1 for r in rs if r["verdict"] == "CONFIRMED_BOTH")
        one = sum(1 for r in rs if r["verdict"] == "CONFIRMED_ONE")
        dd_only = sum(1 for r in rs if r["verdict"] == "CONFIRMED_DNADIFF")
        rej = sum(1 for r in rs if r["verdict"] == "REJECTED")
        cb = sum(1 for r in rs if r["verdict"] == "UNRESOLVABLE_CONTIG_BREAK")
        err = sum(1 for r in rs if r["verdict"].startswith(("ERROR", "NO_", "MISSING", "REF_", "TGT_")))
        print(f"  {cat:<20} n={n:<3} both={both} one={one} dnadiff_only={dd_only} "
              f"rejected={rej} contig-break={cb} err/skip={err}", file=sys.stderr)
    print(f"\nWrote {args.out}/validation.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()
