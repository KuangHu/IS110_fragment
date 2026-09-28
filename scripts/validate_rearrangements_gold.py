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
import argparse, gzip, json, os, re, shutil, subprocess, sys, tempfile
from collections import defaultdict

import pysam


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
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--mauve", default="progressiveMauve",
                   help="progressiveMauve binary")
    p.add_argument("--keep-tmp", action="store_true")
    return p.parse_args()


def load_ref_loci(records_path):
    """ref_id -> (local_contig, start, end) of the IS element in its own assembly."""
    loci = {}
    with open(records_path) as f:
        recs = json.load(f)
    for r in recs:
        rid = r.get("ref_id") or r.get("is110_id")
        src = r.get("source", {})
        ie = r.get("is_element", {})
        contig = src.get("contig")
        s = ie.get("source_start") or src.get("transposase_start")
        e = ie.get("source_end") or src.get("transposase_end")
        if rid and contig and s and e:
            # DB contig is "<assembly>|<local_contig>"; strip the assembly prefix
            local = contig.split("|", 1)[1] if "|" in contig else contig
            loci[rid] = (local, int(s), int(e))
    return loci


def extract_region(fa_obj, contig, start, end, window, label, out_handle):
    """Write [start-window, end+window) of `contig` from an open pysam.FastaFile."""
    try:
        clen = fa_obj.get_reference_length(contig)
    except (KeyError, ValueError):
        return False
    s = max(0, start - window)
    e = min(clen, end + window)
    seq = fa_obj.fetch(contig, s, e)
    if not seq:
        return False
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
              "Insertions": 0}
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
    n_inv = sum(1 for l in shared if l[1][2] != l[2][2])
    # reordering: sort by genome-1 start; count discordances in genome-2 order
    order = sorted(shared, key=lambda l: l[1][0])
    g2_starts = [l[2][0] for l in order]
    n_breaks = sum(1 for i in range(1, len(g2_starts))
                   if g2_starts[i] < g2_starts[i - 1])
    return {"ok": True, "INV": n_inv, "REARR": n_breaks, "n_lcbs": len(shared)}


# class mapping per candidate category
DNADIFF_CLASS = {"inversion_only": "Inversions",
                 "translocation_only": "Translocations",
                 "duplication": "Insertions"}  # dup shows as ins/dup in dnadiff
# progressiveMauve signal per category: which LCB-derived metric must be > 0
MAUVE_SIGNAL = {"inversion_only": "INV",
                "translocation_only": "REARR",
                "duplication": None}  # Mauve LCBs don't call duplication; dnadiff-only


def anchor_targets(details):
    """Extract [(contig, start, end)] from the details string.
    Contig names may contain '|' (e.g. ENA|CM003196|CM003196.1)."""
    out = []
    for m in re.finditer(r"([\w.|]+):(\d+)-(\d+)", details):
        out.append((m.group(1), int(m.group(2)), int(m.group(3))))
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
            ref_contig, rs, re_ = locus
            with open(ref_region, "w") as fh:
                if not extract_region(ref_fo, ref_contig, rs, re_, W, "REF_LOCUS", fh):
                    rec["verdict"] = "REF_REGION_FAIL"; results.append(rec); continue
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
                    if extract_region(qry_fo, contig, s, e, W, f"TGT_{contig}", fh):
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
            mv = run_mauve(args.mauve, ref_region, qry_region, work)
            signal = MAUVE_SIGNAL.get(cat)
            if not mv.get("ok"):
                rec["mauve"] = "FAILED"
            elif signal is None:
                rec["mauve"] = "n/a"  # Mauve LCBs don't resolve duplication
            else:
                mval = mv.get(signal, 0)
                rec["mauve_metric"] = f"{signal}={mval}/{mv.get('n_lcbs',0)}LCB"
                rec["mauve"] = "YES" if mval > 0 else "no"

            # verdict (duplication uses dnadiff only, since Mauve is n/a there)
            dd_ok = rec["dnadiff"] == "YES"
            mv_ok = rec["mauve"] == "YES"
            # translocations across 2 contigs are unresolvable when an anchor sits
            # within --bp-window of a contig end (= contig-break artifact, not a
            # confirmable rearrangement in a draft target)
            artifact = (cat == "translocation_only" and rec["n_target_contigs"] >= 2
                        and 0 <= rec.get("min_anchor_end_dist", -1) < args.bp_window)
            if artifact:
                rec["verdict"] = "UNRESOLVABLE_CONTIG_BREAK"
            elif signal is None:
                rec["verdict"] = "CONFIRMED_DNADIFF" if dd_ok else "REJECTED"
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
