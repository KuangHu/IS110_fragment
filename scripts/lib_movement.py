#!/usr/bin/env python3
"""Shared logic for the "did this element move?" scripts
(is_large_element_movement_db, is_solo_long_movement_db,
is_transposon_movement_db, is_transposon_movement_check).

Each of them carried its own copy of the same four steps -- search the body,
cluster the hits, compare the flanks, give a verdict -- and the same defects:

  * minus-strand hits were discarded, so a target contig assembled in the
    other orientation gave no row (or 'ABSENT');
  * the search used minimap2's default secondary cap, so only the ~50 best
    copies in the whole database were ever seen;
  * the reference's own GCA/GCF twin counted as another genome;
  * a flank cut short by a contig end was divided by the nominal flank length
    and read as "flank differs" -> MOVED;
  * target flanks were taken at the alignment ends, not at the element ends,
    so a partly aligned copy at the SAME site was MOVED;
  * a copy at another site of the reference (a paralog) that sits at the same
    site in the target was MOVED, once per ordered pair of copies.
"""
import os
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_alleles import assembly_core, revcomp, split_target  # noqa: E402

MIN_FLANK = 1000        # shorter than this on either side: not assessable


def write_fasta(label, seq, path):
    with open(path, "w") as fh:
        fh.write(f">{label}\n")
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i+80] + "\n")


def census_args(max_secondary, min_coverage):
    """minimap2 options for 'every copy in the database'.

    asm10 keeps at most 50 secondary chains scoring >= 0.8 of the best one.
    The reference's own full-length copy is always in the database, so
    anything below 80% of it -- the whole 'loose' tier -- was never reported,
    and of the full copies only 50."""
    return ["--secondary=yes", "-N", str(max_secondary),
            "-p", f"{min(0.5, 0.8 * min_coverage):.2f}"]


def same_genome(asm_a, asm_b):
    """True for one assembly and for its GCA/GCF twin."""
    return assembly_core(asm_a) == assembly_core(asm_b)


def parse_clusters(paf_path, qlen_lookup, cluster_gap, rows_per_query=None):
    """PAF -> clusters of nearby blocks per (query, target, STRAND).

    Blocks are chained on the running maximum target end (chaining on the
    previous block's end split a locus when a short block was contained in a
    longer one). Coverage is the union of query intervals; identity is
    matches / alignment columns over the non-redundant blocks."""
    by_qts = defaultdict(list)
    with open(paf_path) as fh:
        for line in fh:
            c = line.split("\t")
            if len(c) < 12:
                continue
            if rows_per_query is not None:
                rows_per_query[c[0]] = rows_per_query.get(c[0], 0) + 1
            by_qts[(c[0], c[5], c[4])].append({
                "qs": int(c[2]), "qe": int(c[3]), "ts": int(c[7]), "te": int(c[8]),
                "matches": int(c[9]), "block": int(c[10])})
    clusters = []
    for (qname, tname, strand) in sorted(by_qts):
        blks = sorted(by_qts[(qname, tname, strand)], key=lambda x: (x["ts"], x["te"]))
        groups, cur, cur_te = [], [blks[0]], blks[0]["te"]
        for b in blks[1:]:
            # one COPY = blocks that also advance along the query in step with
            # the target ('-' strand: backwards). Proximity alone merged two
            # adjacent copies (A-E-E-B) into one cluster and gave a copy with
            # swapped halves coverage 1.0 / "strict".
            prev = cur[-1]
            in_order = (b["qs"] >= prev["qs"] + 0.5 * (prev["qe"] - prev["qs"])) if strand == "+" \
                else (b["qe"] <= prev["qe"] - 0.5 * (prev["qe"] - prev["qs"]))
            contained = prev["ts"] <= b["ts"] and b["te"] <= cur_te
            if b["ts"] - cur_te <= cluster_gap and (in_order or contained):
                cur.append(b)
                cur_te = max(cur_te, b["te"])
            else:
                groups.append(cur)
                cur, cur_te = [b], b["te"]
        groups.append(cur)
        qlen = qlen_lookup.get(qname, 0)
        for g in groups:
            covered, c_s, c_e = 0, None, None
            for s, e in sorted((b["qs"], b["qe"]) for b in g):
                if c_e is None or s > c_e:
                    if c_e is not None:
                        covered += c_e - c_s
                    c_s, c_e = s, e
                else:
                    c_e = max(c_e, e)
            covered += c_e - c_s
            # identity without double counting: drop blocks contained in
            # another block's query interval
            keep = [b for b in g if not any(
                o is not b and o["qs"] <= b["qs"] and b["qe"] <= o["qe"]
                and (o["qe"] - o["qs"]) > (b["qe"] - b["qs"]) for o in g)]
            matches = sum(b["matches"] for b in keep)
            block = sum(b["block"] for b in keep)
            clusters.append({
                "qname": qname, "tname": tname, "strand": strand,
                "ts": min(b["ts"] for b in g), "te": max(b["te"] for b in g),
                "qs": min(b["qs"] for b in g), "qe": max(b["qe"] for b in g),
                "qlen": qlen,
                "cov": covered / qlen if qlen else 0,
                "weighted_ident": matches / block * 100 if block else 0,
                "n_blocks": len(g)})
    return clusters


def projected_ends(cl):
    """Target coordinates of the two ELEMENT ends, projected through the
    unaligned query overhang (0-based half-open). With flanks taken at the
    alignment ends, a copy whose outer 2.8 kb is diverged had 2.8 kb of
    element interior as its 'flank' and was MOVED at the same site."""
    head, tail = cl["qs"], cl["qlen"] - cl["qe"]
    if cl["strand"] == "+":
        return cl["ts"] - head, cl["te"] + tail
    return cl["ts"] - tail, cl["te"] + head


def target_flanks(fa, cl, flank):
    """(left, right) flanks of the cluster in QUERY orientation: `left` is
    what lies before the query's first base. Returns '' for a missing side."""
    tlen = fa.get_reference_length(cl["tname"])
    a, b = projected_ends(cl)
    a, b = max(0, min(tlen, a)), max(0, min(tlen, b))
    lo = fa.fetch(cl["tname"], max(0, a - flank), a) if a > 0 else ""
    hi = fa.fetch(cl["tname"], b, min(tlen, b + flank)) if b < tlen else ""
    if cl["strand"] == "+":
        return lo.upper(), hi.upper()
    return revcomp(hi.upper()), revcomp(lo.upper())


def _best_block(ref_seq, qry_seq, work, threads):
    """Longest same-strand alignment of qry on ref: (query bp aligned, identity)."""
    d = tempfile.mkdtemp(prefix="fl_", dir=work)
    try:
        r, q = os.path.join(d, "r.fa"), os.path.join(d, "q.fa")
        write_fasta("r", ref_seq, r)
        write_fasta("q", qry_seq, q)
        res = subprocess.run(["minimap2", "-x", "asm10", "-c", "--eqx",
                              "-t", str(threads), r, q],
                             check=True, capture_output=True, text=True)
        best = (0, 0.0)
        for line in res.stdout.splitlines():
            c = line.split("\t")
            if len(c) < 12 or c[4] != "+":
                continue
            span, block = int(c[3]) - int(c[2]), int(c[10])
            if span > best[0]:
                best = (span, int(c[9]) / block * 100 if block else 0.0)
        return best
    finally:
        shutil.rmtree(d, ignore_errors=True)


def flank_match(ref_flank, tgt_flank, min_cov, min_id, work, threads):
    """Compare one flank. Returns (state, coverage, identity) with state
    True / False / None. None = NOT ASSESSABLE: one of the two flanks is
    shorter than MIN_FLANK (element at a contig end). Coverage is over the
    shorter of the two flanks, not the nominal flank length; a minimap2
    failure raises instead of reading as 'differs'."""
    denom = min(len(ref_flank), len(tgt_flank))
    if denom < MIN_FLANK:
        return None, 0.0, 0.0
    span, ident = _best_block(ref_flank, tgt_flank, work, threads)
    cov = span / denom
    return (cov >= min_cov and ident >= min_id), cov, ident


def site_verdict(left, right):
    """left/right: states from flank_match()."""
    if left is None or right is None:
        return "FLANK_UNASSESSABLE"
    if left and right:
        return "SYNTENIC"
    if not left and not right:
        return "MOVED"
    return "PARTIAL_MATCH"


def fmt(state, cov, ident):
    return f"{'NA' if state is None else state}({cov:.2f}/{ident:.0f})"


def ref_site_states(ref_fa_path, moved, work, threads, max_slack=2000):
    """What does the REFERENCE genome have at each 'MOVED' target site?

    moved: list of (key, left_flank, right_flank, element_len[, body_seq]) in
    query orientation. The two target flanks are aligned to the whole
    reference assembly:

      empty      both flanks land next to each other -> the reference has the
                 site without the element: the evidence of movement
      occupied   both flanks land ~element_len apart AND the reference
                 sequence between them is the element (when body_seq is
                 given): the reference carries ANOTHER COPY at this very site,
                 the target copy is syntenic with it. (A reference with copies
                 at X and Y against a target with the same two gave 'X moved
                 to Y' and 'Y moved to X'.)
      occupied_other   ~element_len apart but the sequence there is something
                 else (a different element of similar size was replaced)
      absent     the flanks do not both map: the target site is not in the
                 reference -- no statement about the reference is possible
    The gap tolerance scales with the element: with a fixed 2 kb, a 1.5-kb
    IS still sitting at the site read as 'empty'.
    Returns {key: state}."""
    if not moved:
        return {}
    import pysam
    d = tempfile.mkdtemp(prefix="refsite_", dir=work)
    try:
        q = os.path.join(d, "flanks.fa")
        with open(q, "w") as fh:
            for i, m in enumerate(moved):
                fh.write(f">{i}_L\n{m[1]}\n>{i}_R\n{m[2]}\n")
        res = subprocess.run(["minimap2", "-x", "asm10", "-c", "--eqx", "-N", "20",
                              "-t", str(threads), ref_fa_path, q],
                             check=True, capture_output=True, text=True)
        hits = defaultdict(list)
        for line in res.stdout.splitlines():
            c = line.split("\t")
            if len(c) < 12:
                continue
            qlen, qs, qe = int(c[1]), int(c[2]), int(c[3])
            if (qe - qs) < 0.5 * qlen or int(c[9]) / max(1, int(c[10])) < 0.9:
                continue
            ts, te = int(c[7]), int(c[8])
            # project to the flank's element-proximal end
            idx, side = c[0].rsplit("_", 1)
            if c[4] == "+":
                edge = te + (qlen - qe) if side == "L" else ts - qs
            else:
                edge = ts - (qlen - qe) if side == "L" else te + qs
            hits[(int(idx), side)].append((c[5], c[4], edge))
        ref_fa = pysam.FastaFile(ref_fa_path)
        out = {}
        for i, m in enumerate(moved):
            key, elen = m[0], m[3]
            body = m[4] if len(m) > 4 else None
            empty_max = min(max_slack, max(200, 0.3 * elen))
            occ_tol = max(200, 0.2 * elen)
            state = "absent"
            for (ct, st, e_l) in sorted(hits.get((i, "L"), ())):
                for (ct2, st2, e_r) in sorted(hits.get((i, "R"), ())):
                    if ct != ct2 or st != st2:
                        continue
                    gap = (e_r - e_l) if st == "+" else (e_l - e_r)
                    if abs(gap - elen) <= occ_tol:
                        lo, hi = sorted((e_l, e_r))
                        same = True
                        if body:
                            there = ref_fa.fetch(ct, max(0, lo), hi).upper()
                            span, ident = _best_any_strand(there, body, d, threads)
                            same = span >= 0.5 * len(body) and ident >= 90
                        if same:
                            state = "occupied"
                        elif state not in ("occupied", "empty"):
                            state = "occupied_other"
                    elif -200 <= gap <= empty_max and state != "occupied":
                        state = "empty"
            out[key] = state
        ref_fa.close()
        return out
    finally:
        shutil.rmtree(d, ignore_errors=True)


def _best_any_strand(ref_seq, qry_seq, work, threads):
    """(query bp aligned, identity) of the longest alignment, either strand."""
    dd = tempfile.mkdtemp(prefix="bd_", dir=work)
    try:
        r, q = os.path.join(dd, "r.fa"), os.path.join(dd, "q.fa")
        write_fasta("r", ref_seq, r)
        write_fasta("q", qry_seq, q)
        res = subprocess.run(["minimap2", "-x", "asm10", "-c", "--eqx", "-t", str(threads), r, q],
                             check=True, capture_output=True, text=True)
        cov, m_, b_ = [], 0, 0
        for line in res.stdout.splitlines():
            c = line.split("\t")
            if len(c) < 12:
                continue
            cov.append((int(c[2]), int(c[3]))); m_ += int(c[9]); b_ += int(c[10])
        tot, cs, ce = 0, None, None
        for s0, e0 in sorted(cov):
            if ce is None or s0 > ce:
                if ce is not None: tot += ce - cs
                cs, ce = s0, e0
            else:
                ce = max(ce, e0)
        tot += (ce - cs) if ce is not None else 0
        return tot, (m_ / b_ * 100 if b_ else 0.0)
    finally:
        shutil.rmtree(dd, ignore_errors=True)


def refine(verdict, ref_site):
    """A 'MOVED' copy that sits where the reference has another copy is not
    a movement."""
    if verdict != "MOVED":
        return verdict
    if ref_site == "occupied":
        return "SYNTENIC_OTHER_REF_COPY"
    if ref_site in ("empty", "occupied_other"):
        return "MOVED"            # the reference has this site WITHOUT the element
    # Flanks differ, but the target site cannot be found in the reference.
    # That is also what unrelated insertions in both flanks, a deleted element
    # end, or a strain-specific region look like -- not evidence of movement.
    return "DIFFERENT_CONTEXT"


def split_tname(tname):
    return split_target(tname)


def search_and_judge(args, candidates, db_fa, work_root, out_path,
                     id_header, id_values, tag):
    """Steps 3-6 for the *_movement_db scripts.

    candidates: dicts with id, is_seq (body), lf, rf (reference flanks, same
    strand as the body). id_header / id_values(cand) give the leading columns.
    Needs args: ref_asm, ref_fa, species_db, species_mmi, flank, min_coverage,
    min_identity, strict_coverage, strict_identity, flank_match_cov,
    flank_match_id, cluster_gap, max_secondary, threads.
    """
    is_fa = os.path.join(work_root, "bodies.fa")
    qlen = {}
    with open(is_fa, "w") as fh:
        for c in candidates:
            fh.write(f">{c['id']}\n{c['is_seq']}\n")
            qlen[c["id"]] = len(c["is_seq"])
    paf = os.path.join(work_root, "all.paf")
    db_input = args.species_mmi or args.species_db
    print(f"[{args.ref_asm}] minimap2 {tag} bodies vs species DB ({db_input})...",
          flush=True)
    subprocess.run(["minimap2", "-x", "asm10", "-c", "--eqx"]
                   + census_args(args.max_secondary, args.min_coverage)
                   + ["-t", str(args.threads), db_input, is_fa, "-o", paf],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    rows_per_q = {}
    clusters = parse_clusters(paf, qlen, args.cluster_gap, rows_per_q)
    truncated = sorted(q for q, n in rows_per_q.items() if n > args.max_secondary)
    if truncated:
        print(f"[{args.ref_asm}] WARNING: {len(truncated)} queries reached "
              f"--max-secondary {args.max_secondary}; their copy census is "
              f"TRUNCATED (raise it): {truncated[:5]}", flush=True)
    passing = [c for c in clusters
               if c["cov"] >= args.min_coverage
               and c["weighted_ident"] >= args.min_identity
               and not same_genome(split_target(c["tname"])[0], args.ref_asm)]
    print(f"[{args.ref_asm}] {len(clusters)} clusters (both strands), "
          f"{len(passing)} pass cov>={args.min_coverage} id>={args.min_identity}% "
          f"(reference assembly and its GCA/GCF twin excluded)", flush=True)

    cand_by_id = {c["id"]: c for c in candidates}
    rows, moved = [], []
    for ki, cl in enumerate(passing):
        if ki % 500 == 0:
            print(f"  flank-checking {ki}/{len(passing)}", flush=True)
        cand = cand_by_id[cl["qname"]]
        tlf, trf = target_flanks(db_fa, cl, args.flank)
        L = flank_match(cand["lf"], tlf, args.flank_match_cov, args.flank_match_id,
                        work_root, args.threads)
        R = flank_match(cand["rf"], trf, args.flank_match_cov, args.flank_match_id,
                        work_root, args.threads)
        verdict = site_verdict(L[0], R[0])
        if verdict == "MOVED":
            moved.append((ki, tlf, trf, qlen[cl["qname"]], cand["is_seq"]))
        rows.append((ki, cand, cl, L, R, verdict))
    ref_site = ref_site_states(args.ref_fa, moved, work_root, args.threads)

    verdicts = ("MOVED", "DIFFERENT_CONTEXT", "SYNTENIC", "SYNTENIC_OTHER_REF_COPY",
                "PARTIAL_MATCH", "FLANK_UNASSESSABLE")
    n = dict.fromkeys(verdicts, 0)
    n_strict = dict.fromkeys(verdicts, 0)
    with open(out_path + ".tmp", "w") as out:
        out.write("ref_asm\t" + id_header + "\ttgt_asm\ttgt_contig\tts\tte\tstrand\t"
                  "n_blocks\tIS_cov\tIS_ident\tleft_flank_match\tright_flank_match\t"
                  "verdict\tstrict_IS\tref_site\n")
        for ki, cand, cl, L, R, verdict in rows:
            rs = ref_site.get(ki, "na")
            verdict = refine(verdict, rs)
            strict = (cl["cov"] >= args.strict_coverage
                      and cl["weighted_ident"] >= args.strict_identity)
            n[verdict] += 1
            n_strict[verdict] += strict
            tgt_asm, tgt_contig = split_target(cl["tname"])
            out.write(f"{args.ref_asm}\t{id_values(cand)}\t{tgt_asm}\t{tgt_contig}\t"
                      f"{cl['ts']}\t{cl['te']}\t{cl['strand']}\t{cl['n_blocks']}\t"
                      f"{cl['cov']:.3f}\t{cl['weighted_ident']:.1f}\t"
                      f"{fmt(*L)}\t{fmt(*R)}\t{verdict}\t{strict}\t{rs}\n")
    os.replace(out_path + ".tmp", out_path)
    print(f"[{args.ref_asm}] DONE. loose: {n}    strict: {n_strict}")
    print("  MOVED rows are per (reference copy, target copy): one relocation "
          "seen from several references or twins is several rows -- collapse "
          "by target site before counting events.")
    print(f"Wrote {out_path}")
    return n, n_strict
