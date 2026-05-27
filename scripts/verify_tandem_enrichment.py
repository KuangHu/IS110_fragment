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
    return p.parse_args()


# ----------------------------- sequence sets --------------------------------
def load_is_sequences(records_path, min_len):
    """Return (seqs dict {id: seq}, is_loci {contig: [(s,e)]})."""
    with open(records_path) as f:
        records = json.load(f)
    seqs = {}
    is_loci = defaultdict(list)
    for r in records:
        rid = r.get("ref_id") or r.get("is110_id")
        ie = r.get("is_element") or (r.get("source", {}) or {}).get("is_element") or {}
        seq = ie.get("sequence", "")
        # genomic locus to exclude from control sampling
        src = r.get("source", {})
        contig = src.get("contig")
        s = ie.get("source_start")
        e = ie.get("source_end")
        if contig and s and e:
            is_loci[contig].append((min(s, e), max(s, e)))
        if rid and seq and len(seq) >= min_len:
            seqs[rid] = seq.upper()
    return seqs, is_loci


def sample_controls(fa, is_seqs, is_loci, ratio, min_len, rng):
    """Sample length-matched random windows avoiding IS loci. Returns {id: seq}."""
    contigs = list(fa.references)
    lengths = np.array([fa.get_reference_length(c) for c in contigs], dtype=np.int64)
    usable = lengths >= min_len
    contigs = [c for c, u in zip(contigs, usable) if u]
    clens = lengths[usable].astype(float)
    if not contigs:
        return {}
    weights = clens / clens.sum()
    cidx = np.arange(len(contigs))

    # interval lookup for exclusion
    def overlaps_is(contig, s, e):
        for (a, b) in is_loci.get(contig, []):
            if s < b and a < e:
                return True
        return False

    controls = {}
    for rid, seq in is_seqs.items():
        L = len(seq)
        made = 0
        attempts = 0
        while made < ratio and attempts < ratio * 50:
            attempts += 1
            ci = int(rng.choice(cidx, p=weights))
            contig = contigs[ci]
            clen = int(fa.get_reference_length(contig))
            if clen < L:
                continue
            start = rng.integers(0, clen - L + 1)
            end = start + L
            if overlaps_is(contig, start + 1, end):
                continue
            cseq = fa.fetch(contig, int(start), int(end)).upper()
            if not cseq or len(cseq) != L:
                continue
            # skip windows that are mostly N
            if cseq.count("N") > 0.1 * L:
                continue
            controls[f"{rid}__ctrl{made}"] = cseq
            made += 1
    return controls


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

    for _ in range(max_tries):
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
    return s  # could not build a valid shuffle; return original


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
    cmd = [ultra_bin, "--tsv", "-t", str(threads), "-o", out_tsv, fa_path]
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
    fold = (mt / mc) if mc > 0 else float("inf")
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


def fisher_presence(test_pos, test_n, ctrl_pos, ctrl_n):
    table = [[test_pos, test_n - test_pos], [ctrl_pos, ctrl_n - ctrl_pos]]
    odds, p = stats.fisher_exact(table, alternative="greater")
    return odds, p, table


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    print("Loading IS sequences ...", file=sys.stderr, flush=True)
    is_seqs, is_loci = load_is_sequences(args.records, args.min_seq_len)
    print(f"  IS sequences (>= {args.min_seq_len} bp): {len(is_seqs):,}", file=sys.stderr)
    if not is_seqs:
        sys.exit("No IS sequences to test.")

    import random as _random
    pyrng = _random.Random(args.seed)

    print("Sampling length-matched control windows ...", file=sys.stderr, flush=True)
    fa = pysam.FastaFile(args.genome_db)
    genomic = sample_controls(fa, is_seqs, is_loci, args.control_ratio,
                              args.min_seq_len, rng)
    print(f"  genomic control windows: {len(genomic):,} "
          f"(target {len(is_seqs)*args.control_ratio:,})", file=sys.stderr)

    print("Generating dinucleotide-shuffle control ...", file=sys.stderr, flush=True)
    shuffle = {f"{sid}__dishuf": dinuc_shuffle(seq, pyrng)
               for sid, seq in is_seqs.items()}
    print(f"  dinuc-shuffle controls: {len(shuffle):,}", file=sys.stderr)

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

    report = {"n_IS": len(test_rows), "min_large_period": mlp,
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
            for flag, _ in pres_flags:
                tp = sum(1 for r in test_rows if r[flag])
                cp = sum(1 for r in ctrl_rows if r[flag])
                odds, p, _ = fisher_presence(tp, len(test_rows), cp, len(ctrl_rows))
                comp[flag] = {"IS_frac": tp / len(test_rows),
                              "control_frac": cp / len(ctrl_rows),
                              "IS_positive": tp, "control_positive": cp,
                              "odds_ratio": odds, "p_fisher": p}
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
                         f"{d['control_frac']*100:5.1f}%   OR={d['odds_ratio']:.2f}  "
                         f"p_fisher={d['p_fisher']:.1e}")
        lines.append("  Density (tandem bp / length):")
        for metric in metrics:
            d = comp[metric]
            lines.append(f"    {metric:<20} IS {d['mean_IS']*100:5.2f}%  vs "
                         f"{d['mean_control']*100:5.2f}%   fold={d['fold_enrichment']:.2f}x  "
                         f"p_perm={d['p_permutation']:.1e}")
    text = "\n".join(lines)
    with open(os.path.join(args.out, "enrichment_report.txt"), "w") as f:
        f.write(text + "\n")
    print("\n" + text, file=sys.stderr)


if __name__ == "__main__":
    main()
