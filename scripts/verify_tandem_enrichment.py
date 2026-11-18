#!/usr/bin/env python3
"""Verify tandem-repeat enrichment in IS elements vs length-matched genomic background.

Claim tested: tandem repeats occur more often inside IS elements than in normal
genomic DNA. We use TWO independent gold-standard tandem callers (TRF + ULTRA)
and a length-matched random-window control sampled from the same genome DB.

Design
------
  Test set      : IS-element sequences from records.json (one per record).
  Control set 1 : "genomic" — for each IS of length L, sample --control-ratio random
                  windows of length L from the same DB, rejecting windows overlapping
                  a known IS locus. Tests enrichment over the genome at large
                  (confounds base composition + mechanism).
  Control set 2 : "shuffle" — a dinucleotide-preserving shuffle of each IS (same 1-mer
                  AND 2-mer frequencies, same first/last base). Tandem ORDER is
                  destroyed but composition preserved, so enrichment over this control
                  cannot be explained by base composition — it isolates a mechanistic
                  signal (active duplication).
  Callers       : TRF and ULTRA, identical parameters on every set.
  Metric      : per-sequence tandem density = (bp covered by a tandem array) / len.
                Computed for ALL periods and for "large" units (period >= --min-large-period).
  Stats       : Fisher exact (presence/absence), Mann-Whitney U (density),
                label-permutation test (mean-density difference), with fold-enrichment.
  Agreement   : a sequence is "tandem+" by intersection (called by BOTH tools).

Outputs (in --out dir):
  test.fa, control.fa
  trf_test.dat, trf_control.dat, ultra_test.json, ultra_control.json
  per_sequence.tsv   (seq_id, set, length, trf_bp, ultra_bp, both_bp, ... )
  enrichment_report.txt / enrichment_report.json

Usage:
  verify_tandem_enrichment.py --records records.json --genome-db genomes.fa \\
      --out tandem_verify/ --control-ratio 3 --threads 16
"""
import argparse, json, os, random, subprocess, sys
from collections import defaultdict

import pysam
import numpy as np
from scipy import stats


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--records", required=True, help="records.json (IS element seqs)")
    p.add_argument("--genome-db", required=True, help="Indexed FASTA (.fai present)")
    p.add_argument("--out", required=True)
    p.add_argument("--control-ratio", type=int, default=3,
                   help="Control windows per IS (default 3)")
    p.add_argument("--min-large-period", type=int, default=50,
                   help="Period >= this counts as a 'large' (structural) tandem unit")
    p.add_argument("--min-seq-len", type=int, default=300,
                   help="Skip IS sequences shorter than this")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--trf", default="trf")
    p.add_argument("--ultra", default="ultra")
    p.add_argument("--threads", type=int, default=16,
                   help="Threads for ULTRA (TRF is single-threaded)")
    p.add_argument("--n-permute", type=int, default=10000)
    p.add_argument("--cluster-identity", type=float, default=0.95,
                   help="one test sequence per cd-hit-est cluster at this identity "
                        "(default 0.95, the project's IS110 threshold; 1.0 = exact "
                        "sequences only)")
    p.add_argument("--keep-duplicates", action="store_true",
                   help="Test every record, even identical element sequences "
                        "(pseudo-replicated; default: one per distinct sequence)")
    p.add_argument("--all-records", action="store_true",
                   help="Also test records without a called boundary "
                        "(transposase-only / symmetric fallback)")
    return p.parse_args()


# ----------------------------- sequence sets --------------------------------
_RC = str.maketrans("ACGTN", "TGCAN")


def _revcomp(s):
    return s.translate(_RC)[::-1]


def _contig_key(name):
    c = name.split("|")[-1]
    return c[3:] if c.startswith("NZ_") else c


def cluster_representatives(seqs, identity, work_dir, threads):
    """One sequence per cd-hit-est cluster (both strands, >= `identity`,
    length within the same fraction). Exact-string de-duplication left an
    element and its one-SNP relatives as separate "independent" observations:
    Fisher's p moved from 0.5 to 1e-59 by cloning one sequence."""
    if len(seqs) < 2:
        return seqs
    fa, out = os.path.join(work_dir, "dedup_in.fa"), os.path.join(work_dir, "dedup_out.fa")
    ids = list(seqs)
    with open(fa, "w") as f:
        for i, k in enumerate(ids):
            f.write(f">s{i}\n{seqs[k]}\n")
    word = 10 if identity >= 0.95 else 8
    subprocess.run(["cd-hit-est", "-i", fa, "-o", out, "-c", str(identity),
                    "-s", str(identity), "-n", str(word), "-d", "0",
                    "-T", str(threads), "-M", "0"],
                   check=True, stdout=subprocess.DEVNULL)
    keep = {}
    with open(out) as f:
        for line in f:
            if line.startswith(">"):
                k = ids[int(line[2:].split()[0])]
                keep[k] = seqs[k]
    return keep



def load_is_sequences(records_path, min_len, called_only=True, dedup=True):
    """Return (seqs dict {id: seq}, is_loci {contig: [(s,e)]}, stats).

    called_only  keep only records whose boundaries were CALLED
                 (method empty_vs_filled_decomposition). For the others the
                 "element" is the bare transposase ORF or a symmetric guess,
                 so its ends -- where IS tandem repeats sit -- are not in the
                 sequence at all. Ignored for records without that label
                 (older record files).
    dedup        one sequence per distinct element. A catalogue holds the
                 same IS copy from thousands of genomes; treated as
                 independent observations they made every p-value
                 meaningless (pseudo-replication), while each control window
                 is a genuinely independent draw."""
    with open(records_path) as f:
        records = json.load(f)
    seqs = {}
    is_loci = defaultdict(list)
    stats_ = defaultdict(int)
    seen_seq = set()
    # "labelled" = the file uses the current method vocabulary at all. Testing
    # for a successful call meant that a file with ONLY fallback records was
    # taken for an old unlabelled one and every record was kept.
    # (any method label counts: a list of known labels went stale as soon as
    # a new one was added, and a file holding only the new label was read as
    # an unlabelled old file -- every record kept)
    labelled = any((r.get("boundary_evidence") or {}).get("method") for r in records)
    for r in records:
        rid = r.get("ref_id") or r.get("is110_id")
        ie = r.get("is_element") or (r.get("source", {}) or {}).get("is_element") or {}
        seq = ie.get("sequence", "")
        # genomic locus to exclude from control sampling
        src = r.get("source", {})
        contig = src.get("contig")
        s = ie.get("source_start") or ie.get("start")      # either record schema
        e = ie.get("source_end") or ie.get("end")
        if contig and s and e:
            # keyed by the bare contig accession, so the GCA and GCF copies of
            # one genome (CP012378.1 / NZ_CP012378.1) share the exclusion --
            # the twin's copy of the very element under test was a "control"
            is_loci[_contig_key(contig)].append((min(s, e), max(s, e)))
        if not (rid and seq and len(seq) >= min_len):
            stats_["too_short_or_missing"] += 1
            continue
        if (called_only and labelled and (r.get("boundary_evidence") or {}).get("method")
                != "empty_vs_filled_decomposition"):
            stats_["boundary_not_called"] += 1
            continue
        seq = seq.upper()
        if dedup:
            # one orientation: an element and its reverse complement are the
            # same sequence
            seq = min(seq, _revcomp(seq))
            if seq in seen_seq:
                stats_["duplicate_sequence"] += 1
                continue
            seen_seq.add(seq)
        seqs[rid] = seq
    stats_["kept"] = len(seqs)
    stats_["records"] = len(records)
    return seqs, is_loci, dict(stats_)


def sample_controls(fa, is_seqs, is_loci, ratio, min_len, rng):
    """Length-matched random windows avoiding known IS loci.
    Returns ({id: seq}, n_short) -- n_short = controls that could not be drawn.

    A window is drawn uniformly over all ELIGIBLE START POSITIONS for that
    length (contig weight = clen - L + 1). Weighting by contig length gave a
    1-kb contig one third of the 1-kb windows next to a 2-kb contig (it has
    1 of 1,002 possible starts). No window is used twice."""
    contigs = list(fa.references)
    lengths = np.array([fa.get_reference_length(c) for c in contigs], dtype=np.int64)

    def overlaps_is(contig, s, e):
        for (a, b) in is_loci.get(_contig_key(contig), []):
            if s <= b and a <= e:
                return True
        return False

    controls, used, n_short = {}, set(), 0
    weights_by_len = {}
    for rid, seq in is_seqs.items():
        L = len(seq)
        if L not in weights_by_len:
            w = np.clip(lengths - L + 1, 0, None).astype(float)
            weights_by_len[L] = (w / w.sum()) if w.sum() > 0 else None
        w = weights_by_len[L]
        made = attempts = 0
        while w is not None and made < ratio and attempts < ratio * 50:
            attempts += 1
            ci = int(rng.choice(len(contigs), p=w))
            contig, clen = contigs[ci], int(lengths[ci])
            start = int(rng.integers(0, clen - L + 1))
            end = start + L
            if (contig, start) in used or overlaps_is(contig, start + 1, end):
                continue
            cseq = fa.fetch(contig, start, end).upper()
            if not cseq or len(cseq) != L:
                continue
            # skip windows that are mostly N
            if cseq.count("N") > 0.1 * L:
                continue
            used.add((contig, start))
            controls[f"{rid}__ctrl{made}"] = cseq
            made += 1
        n_short += ratio - made
    return controls, n_short


def write_fasta(path, seqs):
    with open(path, "w") as f:
        for name, seq in seqs.items():
            f.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")


def dinuc_shuffle(seq, pyrng, max_tries=25):
    """Altschul-Erikson dinucleotide-preserving shuffle.

    Returns a permutation of `seq` with identical dinucleotide composition and
    identical first/last base. This is the composition-matched null: it destroys
    any tandem ORDER while preserving 1-mer AND 2-mer frequencies, so any tandem
    enrichment over this control cannot be explained by base composition alone.
    Sequences containing non-ACGT or shorter than 4 bp are returned unchanged.
    """
    s = seq.upper()
    n = len(s)
    if n < 4 or set(s) - set("ACGT"):
        return s
    first, last = s[0], s[-1]
    # out-edge multiset per vertex
    base_edges = defaultdict(list)
    for i in range(n - 1):
        base_edges[s[i]].append(s[i + 1])

    for _ in range(max(max_tries, 2000)):
        # 1) pick a random "last edge" for every vertex except `last`
        last_edge = {}
        for v, es in base_edges.items():
            if v == last:
                continue
            last_edge[v] = es[pyrng.randrange(len(es))]
        # 2) verify those last-edges form a tree rooted at `last` (no cycles)
        ok = True
        for v in base_edges:
            if v == last:
                continue
            seen = set()
            cur = v
            while cur != last:
                if cur in seen or cur not in last_edge:
                    ok = False
                    break
                seen.add(cur)
                cur = last_edge[cur]
            if not ok:
                break
        if not ok:
            continue
        # 3) per vertex: shuffle remaining edges, append the chosen last-edge
        out = {}
        for v, es in base_edges.items():
            edges = list(es)
            if v != last:
                edges.remove(last_edge[v])
                pyrng.shuffle(edges)
                edges.append(last_edge[v])
            else:
                pyrng.shuffle(edges)
            out[v] = edges
        # 4) traverse the Eulerian path from `first`
        res = [first]
        idx = defaultdict(int)
        cur = first
        for _ in range(n - 1):
            nxt = out[cur][idx[cur]]
            idx[cur] += 1
            res.append(nxt)
            cur = nxt
        return "".join(res)
    # No valid shuffle found. Returning the ORIGINAL here put the unshuffled
    # element into its own null set (86 of 100 tries for a low-complexity
    # sequence at 25 attempts); the caller drops the pair instead.
    return None


# ----------------------------- tandem callers -------------------------------
def run_trf(trf_bin, fa_path, out_dat):
    """Run TRF in -ngs mode. Returns {seq_id: [(start,end,period), ...]} (1-based)."""
    # match mismatch delta PM PI minscore maxperiod
    cmd = [trf_bin, fa_path, "2", "7", "7", "80", "10", "50", "500", "-h", "-ngs"]
    with open(out_dat, "w") as fout:
        subprocess.run(cmd, check=True, stdout=fout, stderr=subprocess.DEVNULL)
    hits = defaultdict(list)
    cur = None
    with open(out_dat) as f:
        for line in f:
            if line.startswith("@"):
                cur = line[1:].strip().split()[0]
            elif cur:
                parts = line.split()
                if len(parts) < 4:
                    continue
                try:
                    s = int(parts[0]); e = int(parts[1]); period = int(parts[2])
                except ValueError:
                    continue
                hits[cur].append((s, e, period))
    return hits


def run_ultra(ultra_bin, fa_path, out_tsv, threads=1):
    """Run ULTRA (TSV output: SeqID Start End Period Score ...).
    Returns {seq_id: [(start,end,period), ...]} as 1-based inclusive intervals."""
    # -p 500: same largest period as TRF above (ULTRA's default is 100, so
    # every "large" unit of 100-500 bp could be called by TRF only and never
    # counted as agreed)
    cmd = [ultra_bin, "--tsv", "-p", "500", "-t", str(threads), "-o", out_tsv, fa_path]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    hits = defaultdict(list)
    with open(out_tsv) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4 or parts[0] == "SeqID":
                continue
            sid = parts[0].split()[0]
            try:
                start0 = int(parts[1]); end = int(parts[2]); period = int(parts[3])
            except ValueError:
                continue
            # ULTRA Start is 0-based half-open [start,end) -> 1-based inclusive
            hits[sid].append((start0 + 1, end, period))
    return hits


def covered_bp(intervals, min_period=0):
    """Union length of intervals with period >= min_period."""
    iv = sorted((s, e) for (s, e, p) in intervals if p >= min_period)
    if not iv:
        return 0
    total = 0
    cs, ce = iv[0]
    for s, e in iv[1:]:
        if s <= ce:
            ce = max(ce, e)
        else:
            total += ce - cs + 1
            cs, ce = s, e
    total += ce - cs + 1
    return total


# ----------------------------- statistics -----------------------------------
def fold_and_test(test_vals, ctrl_vals, n_permute, rng):
    test_vals = np.asarray(test_vals, float)
    ctrl_vals = np.asarray(ctrl_vals, float)
    mt, mc = test_vals.mean(), ctrl_vals.mean()
    # 0/0 is "no estimate", not infinite enrichment; inf is not JSON either
    fold = (mt / mc) if mc > 0 else None
    # Mann-Whitney
    try:
        u, p_mw = stats.mannwhitneyu(test_vals, ctrl_vals, alternative="greater")
    except ValueError:
        p_mw = float("nan")
    # permutation on mean difference
    obs = mt - mc
    allv = np.concatenate([test_vals, ctrl_vals])
    n_t = len(test_vals)
    ge = 0
    for _ in range(n_permute):
        rng.shuffle(allv)
        if (allv[:n_t].mean() - allv[n_t:].mean()) >= obs:
            ge += 1
    p_perm = (ge + 1) / (n_permute + 1)
    return mt, mc, fold, p_mw, p_perm


def paired_sign_flip(test_vals, ctrl_vals, n_permute, rng):
    """One-sided p for mean(test - control) > 0 when each control belongs to
    one test sequence (the dinucleotide shuffle OF that element): labels are
    exchanged within pairs. The unpaired permutation ignored the pairing
    (p = 0.47 where the exact paired test gives 0.002)."""
    d = np.asarray(test_vals, float) - np.asarray(ctrl_vals, float)
    obs = d.mean()
    signs = rng.choice([-1.0, 1.0], size=(n_permute, len(d)))
    ge = int(((signs * d).mean(axis=1) >= obs).sum())
    return (ge + 1) / (n_permute + 1)


def fisher_presence(test_pos, test_n, ctrl_pos, ctrl_n):
    table = [[test_pos, test_n - test_pos], [ctrl_pos, ctrl_n - ctrl_pos]]
    odds, p = stats.fisher_exact(table, alternative="greater")
    return odds, p, table


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    print("Loading IS sequences ...", file=sys.stderr, flush=True)
    is_seqs, is_loci, load_stats = load_is_sequences(
        args.records, args.min_seq_len, called_only=not args.all_records,
        dedup=not args.keep_duplicates)
    if not args.keep_duplicates and args.cluster_identity < 1.0 and is_seqs:
        n0 = len(is_seqs)
        is_seqs = cluster_representatives(is_seqs, args.cluster_identity, args.out,
                                          args.threads)
        load_stats["near_identical_collapsed"] = n0 - len(is_seqs)
        load_stats["kept"] = len(is_seqs)
    print(f"  IS sequences (>= {args.min_seq_len} bp): {len(is_seqs):,}  {load_stats}",
          file=sys.stderr)
    if not is_seqs:
        sys.exit("No IS sequences to test.")

    import random as _random
    pyrng = _random.Random(args.seed)

    print("Sampling length-matched control windows ...", file=sys.stderr, flush=True)
    fa = pysam.FastaFile(args.genome_db)
    genomic, n_short = sample_controls(fa, is_seqs, is_loci, args.control_ratio,
                                       args.min_seq_len, rng)
    print(f"  genomic control windows: {len(genomic):,} "
          f"(target {len(is_seqs)*args.control_ratio:,}; {n_short:,} could not be drawn)",
          file=sys.stderr)

    print("Generating dinucleotide-shuffle control ...", file=sys.stderr, flush=True)
    shuffle = {f"{sid}__dishuf": dinuc_shuffle(seq, pyrng)
               for sid, seq in is_seqs.items()}
    n_noshuf = sum(1 for v in shuffle.values() if v is None)
    shuffle = {k: v for k, v in shuffle.items() if v is not None}
    print(f"  dinuc-shuffle controls: {len(shuffle):,} "
          f"({n_noshuf} elements could not be shuffled and have no shuffle control)",
          file=sys.stderr)

    # Write FASTAs with tool-safe integer IDs. ULTRA truncates sequence names at
    # '|' (and our ref_ids contain '|'), which silently zeroed all ULTRA lookups;
    # using safe IDs and mapping back fixes it for any caller's name quirks.
    orig_sets = {"IS": is_seqs, "genomic": genomic, "shuffle": shuffle}
    prefixes = {"IS": "i", "genomic": "g", "shuffle": "s"}
    pathnames = {"IS": "test.fa", "genomic": "control_genomic.fa",
                 "shuffle": "control_shuffle.fa"}
    safe_seqs_by_set = {}   # setname -> {safe_id: seq}
    safe2orig_by_set = {}   # setname -> {safe_id: orig_id}
    fastas = {}             # setname -> fasta path
    for setname, seqs in orig_sets.items():
        safe, s2o = {}, {}
        for i, (orig, seq) in enumerate(seqs.items()):
            sid = f"{prefixes[setname]}{i}"
            safe[sid] = seq
            s2o[sid] = orig
        safe_seqs_by_set[setname] = safe
        safe2orig_by_set[setname] = s2o
        path = os.path.join(args.out, pathnames[setname])
        fastas[setname] = path
        write_fasta(path, safe)

    mlp = args.min_large_period

    def build_rows(seqs, trf_hits, ultra_hits, setname, s2o):
        rows = []
        for sid, seq in seqs.items():
            L = len(seq)
            trf_all = covered_bp(trf_hits.get(sid, []), 0)
            trf_lg = covered_bp(trf_hits.get(sid, []), mlp)
            ul_all = covered_bp(ultra_hits.get(sid, []), 0)
            ul_lg = covered_bp(ultra_hits.get(sid, []), mlp)
            rows.append({
                "seq_id": s2o.get(sid, sid), "set": setname, "length": L,
                "trf_density_all": trf_all / L, "trf_density_large": trf_lg / L,
                "ultra_density_all": ul_all / L, "ultra_density_large": ul_lg / L,
                "both_present_all": (trf_all > 0 and ul_all > 0),
                "both_present_large": (trf_lg > 0 and ul_lg > 0),
            })
        return rows

    # Run both callers on each set (keyed by safe IDs)
    rows_by_set = {}
    for setname in ("IS", "genomic", "shuffle"):
        seqs = safe_seqs_by_set[setname]
        path = fastas[setname]
        if not seqs:
            rows_by_set[setname] = []
            continue
        print(f"Running TRF + ULTRA on '{setname}' ({len(seqs):,} seqs) ...",
              file=sys.stderr, flush=True)
        trf_hits = run_trf(args.trf, path, os.path.join(args.out, f"trf_{setname}.dat"))
        ultra_hits = run_ultra(args.ultra, path,
                               os.path.join(args.out, f"ultra_{setname}.tsv"),
                               threads=args.threads)
        rows_by_set[setname] = build_rows(seqs, trf_hits, ultra_hits, setname,
                                          safe2orig_by_set[setname])

    test_rows = rows_by_set["IS"]

    # per-sequence table (all sets)
    with open(os.path.join(args.out, "per_sequence.tsv"), "w") as f:
        cols = ["seq_id", "set", "length", "trf_density_all", "trf_density_large",
                "ultra_density_all", "ultra_density_large",
                "both_present_all", "both_present_large"]
        f.write("\t".join(cols) + "\n")
        for setname in ("IS", "genomic", "shuffle"):
            for r in rows_by_set[setname]:
                f.write("\t".join(str(r[c]) for c in cols) + "\n")

    metrics = ["trf_density_all", "ultra_density_all",
               "trf_density_large", "ultra_density_large"]
    pres_flags = [("both_present_all", "any period"),
                  ("both_present_large", f"period>={mlp}bp")]

    report = {"n_IS": len(test_rows), "input_filter": load_stats,
              "deduplicated": not args.keep_duplicates,
              "called_boundaries_only": not args.all_records,
              "note": "p-values are per test, not corrected for the "
                      "12 comparisons reported",
              "min_large_period": mlp,
              "control_ratio": args.control_ratio, "comparisons": {}}

    # IS vs each control set
    for ctrl_set in ("genomic", "shuffle"):
        ctrl_rows = rows_by_set[ctrl_set]
        comp = {"n_control": len(ctrl_rows)}
        if ctrl_rows:
            for metric in metrics:
                mt, mc, fold, p_mw, p_perm = fold_and_test(
                    [r[metric] for r in test_rows], [r[metric] for r in ctrl_rows],
                    args.n_permute, rng)
                comp[metric] = {"mean_IS": mt, "mean_control": mc,
                                "fold_enrichment": fold,
                                "p_mannwhitney": p_mw, "p_permutation": p_perm}
                if ctrl_set == "shuffle":
                    # paired: each shuffle belongs to one element
                    by_id = {r["seq_id"]: r[metric] for r in test_rows}
                    pairs = [(by_id[r["seq_id"][:-len("__dishuf")]], r[metric])
                             for r in ctrl_rows
                             if r["seq_id"][:-len("__dishuf")] in by_id]
                    if pairs:
                        comp[metric]["p_paired_sign_flip"] = paired_sign_flip(
                            [a for a, _ in pairs], [b for _, b in pairs],
                            args.n_permute, rng)
                        comp[metric]["n_pairs"] = len(pairs)
            for flag, _ in pres_flags:
                tp = sum(1 for r in test_rows if r[flag])
                cp = sum(1 for r in ctrl_rows if r[flag])
                odds, p, _ = fisher_presence(tp, len(test_rows), cp, len(ctrl_rows))
                comp[flag] = {"IS_frac": tp / len(test_rows),
                              "control_frac": cp / len(ctrl_rows),
                              "IS_positive": tp, "control_positive": cp,
                              "odds_ratio": odds if np.isfinite(odds) else None,
                              "p_fisher": p}
        report["comparisons"][ctrl_set] = comp

    with open(os.path.join(args.out, "enrichment_report.json"), "w") as f:
        json.dump(report, f, indent=2)

    # human-readable
    lines = [f"Tandem enrichment  |  IS n={report['n_IS']:,}  "
             f"(large period >= {mlp} bp)"]
    ctrl_labels = {"genomic": "length-matched random genomic windows "
                              "(composition + mechanism)",
                   "shuffle": "dinucleotide-shuffle of each IS "
                              "(isolates mechanism beyond composition)"}
    for ctrl_set in ("genomic", "shuffle"):
        comp = report["comparisons"][ctrl_set]
        lines.append("")
        lines.append(f"=== IS vs {ctrl_set} control — {ctrl_labels[ctrl_set]} "
                     f"(n={comp.get('n_control', 0):,}) ===")
        if not comp.get("n_control"):
            lines.append("  (no control sequences)")
            continue
        lines.append("  Presence (tandem called by BOTH TRF and ULTRA):")
        for flag, lab in pres_flags:
            d = comp[flag]
            lines.append(f"    [{lab:<14}] IS {d['IS_frac']*100:5.1f}%  vs "
                         f"{d['control_frac']*100:5.1f}%   OR={'n/a' if d['odds_ratio'] is None else format(d['odds_ratio'], '.2f')}  "
                         f"p_fisher={d['p_fisher']:.1e}")
        lines.append("  Density (tandem bp / length):")
        for metric in metrics:
            d = comp[metric]
            lines.append(f"    {metric:<20} IS {d['mean_IS']*100:5.2f}%  vs "
                         f"{d['mean_control']*100:5.2f}%   fold="
                         f"{'n/a' if d['fold_enrichment'] is None else format(d['fold_enrichment'], '.2f') + 'x'}  "
                         + (f"p_paired={d['p_paired_sign_flip']:.1e}" if "p_paired_sign_flip" in d
                            else f"p_perm={d['p_permutation']:.1e}"))
    text = "\n".join(lines)
    with open(os.path.join(args.out, "enrichment_report.txt"), "w") as f:
        f.write(text + "\n")
    print("\n" + text, file=sys.stderr)


if __name__ == "__main__":
    main()
