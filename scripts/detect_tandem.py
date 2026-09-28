#!/usr/bin/env python3
"""Detect tandem repeats inside a single sequence via self-alignment.

For each sequence in the input set (V_ref + variants from observations), run
minimap2 self-vs-self. Filter out the trivial diagonal hits (qs==ts) and
inspect the remaining alignments — they reveal regions that occur multiple
times tandemly within the sequence.

For each non-diagonal hit, compute the offset = target_start - query_start.
Same-offset hits with consistent query positions ⇒ tandem repeats.

Outputs:
  tandem_findings.json — per-sequence summary of tandem regions
  tandem_findings.tsv  — flat table: seq_id, repeat_unit_bp, n_copies, ...
"""
import argparse, json, os, subprocess, sys, tempfile
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--records", required=True,
                   help="records_final/records.json (for V_ref sequences)")
    p.add_argument("--obs", required=True,
                   help="observations.json (for variant sequences)")
    p.add_argument("--out", required=True, help="output directory")
    p.add_argument("--min-seq-len", type=int, default=1500,
                   help="Skip sequences shorter than this (default 1500)")
    p.add_argument("--min-repeat-len", type=int, default=50,
                   help="Min length of a tandem unit (default 50)")
    p.add_argument("--min-identity", type=float, default=95,
                   help="Min alignment identity (default 95)")
    p.add_argument("--min-copies", type=int, default=2,
                   help="Min number of tandem copies to report (default 2 — 1 extra beyond original)")
    p.add_argument("--threads", type=int, default=8)
    return p.parse_args()


def self_minimap2(seq, work_dir, threads):
    """Self-align a single sequence. Return PAF lines as parsed tuples."""
    fa = os.path.join(work_dir, "s.fa")
    paf = os.path.join(work_dir, "s.paf")
    with open(fa, "w") as f:
        f.write(">s\n" + seq + "\n")
    subprocess.run(["minimap2", "-x", "asm10", "-c", "--eqx",
                    "-f", "0",  # keep ALL minimizers: the default drops the most frequent
                    # ones, which in a self-alignment are exactly the repeats sought
                    # (fna project: 0 of 60 IS1 copies found without it, 57 with)
                    "-X",  # skip self-self diagonal
                    "-t", str(threads), fa, fa, "-o", paf],
                   check=True, capture_output=True)
    hits = []
    with open(paf) as f:
        for line in f:
            c = line.split("\t")
            if len(c) < 12: continue
            qs, qe = int(c[2]), int(c[3])
            strand = c[4]
            ts, te = int(c[7]), int(c[8])
            matches, block = int(c[9]), int(c[10])
            ident = matches / block * 100 if block > 0 else 0
            hits.append((qs, qe, ts, te, strand, ident))
    return hits


def find_tandem_unit(hits, min_repeat_len, min_identity):
    """Group non-diagonal hits by offset to find tandem repeats.

    Returns list of {repeat_unit_bp, n_copies, span_start, span_end, offset}
    """
    # Only forward-strand hits with offset > 0
    fwd_hits = [(qs, qe, ts, te, ident) for (qs, qe, ts, te, strand, ident) in hits
                if strand == "+" and ident >= min_identity
                and ts > qs and qe - qs >= min_repeat_len]
    if not fwd_hits:
        return []

    # Group by offset (rounded to nearest 50)
    by_offset = defaultdict(list)
    for h in fwd_hits:
        qs, qe, ts, te, ident = h
        offset = ts - qs
        bucket = round(offset / 50) * 50
        by_offset[bucket].append(h)

    results = []
    for offset, group in by_offset.items():
        if offset < min_repeat_len: continue
        # Combined query span
        q_starts = [h[0] for h in group]
        q_ends = [h[1] for h in group]
        span_s = min(q_starts)
        span_e = max(q_ends)
        # Number of tandem copies = n_hits + 1 (each hit shows pair of copies)
        # But hits can chain — better: compute total hit span / offset (= number of repeats)
        # Simpler heuristic: # tandem copies ≈ (max_te - min_qs) / offset + 1
        max_te = max(h[3] for h in group)
        min_qs = min(h[0] for h in group)
        n_copies_est = round((max_te - min_qs) / offset) + 1
        if n_copies_est < 2: continue
        results.append({
            "repeat_unit_bp": offset,
            "n_copies": n_copies_est,
            "span_start": min_qs,
            "span_end": max_te,
            "total_span_bp": max_te - min_qs,
            "n_hits": len(group),
        })
    # Sort: most repeats first, larger unit first
    results.sort(key=lambda r: (-r["n_copies"], -r["repeat_unit_bp"]))

    # Deduplicate: if a larger offset finding's span is contained in a smaller
    # offset finding's span AND the larger is a near-multiple of the smaller,
    # drop the larger (it's just a derivative within the same tandem).
    # Smallest unit = true repeat.
    by_unit = sorted(results, key=lambda r: r["repeat_unit_bp"])
    kept = []
    for r in by_unit:
        is_derivative = False
        for k in kept:
            # is r's offset a multiple of k's offset (within tolerance)?
            ratio = r["repeat_unit_bp"] / k["repeat_unit_bp"]
            if abs(ratio - round(ratio)) > 0.1: continue
            if round(ratio) < 2: continue
            # And r's span overlaps k's span substantially
            ov_s = max(r["span_start"], k["span_start"])
            ov_e = min(r["span_end"], k["span_end"])
            if ov_e - ov_s > 0.5 * (r["span_end"] - r["span_start"]):
                is_derivative = True
                break
        if not is_derivative:
            kept.append(r)
    kept.sort(key=lambda r: (-r["n_copies"], -r["repeat_unit_bp"]))
    return kept


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    print("Collecting sequences...", file=sys.stderr, flush=True)
    seq_set = {}  # name -> {seq, source, kind, cds_in_seq}

    # V_ref sequences + their CDS position (in V_ref-internal coords)
    v_ref_cds = {}  # v1_id -> (cds_s, cds_e) in V_ref-internal coords
    with open(args.records) as f:
        records = json.load(f)
    for r in records:
        v1_id = r.get("is110_id") or r.get("ref_id")
        if not v1_id: continue
        src = r.get("source", {})
        ie = src.get("is_element") or r.get("is_element") or {}
        tn = src.get("transposase_cds") or r.get("transposase_cds") or {}
        seq = ie.get("sequence", "")
        is_len = ie.get("length", 0)
        if "start" in ie and "start" in tn:
            strand = ie.get("strand", "+")
            if strand == "+":
                cs = max(0, tn["start"] - ie["start"])
                ce = min(is_len, tn["end"] - ie["start"] + 1)
            else:
                cs = max(0, ie["end"] - tn["end"])
                ce = min(is_len, ie["end"] - tn["start"] + 1)
        else:
            off5 = ie.get("start_offset_5p", 0)
            tnp_len = tn.get("length", 0)
            cs = max(0, -off5)
            ce = min(is_len, cs + tnp_len)
        v_ref_cds[v1_id] = (cs, ce)
        if not seq or len(seq) < args.min_seq_len: continue
        seq_set[f"vref_{v1_id}"] = {
            "seq": seq, "len": len(seq),
            "kind": "v_ref", "v1_parent": v1_id,
            "source_target": v1_id,
            "cds_intervals": [(cs, ce)],
        }

    # Variant sequences (insertion + deletion). CDS mapped from V_ref via blocks.
    if not os.path.exists(args.obs):
        # No variants stage was run; only V_ref sequences available
        obs = []
    else:
        with open(args.obs) as f:
            obs = json.load(f)
    seen_vids = set()
    for o in obs:
        if o["category"] not in ("insertion", "deletion"): continue
        vid = o.get("variant_id")
        seq = o.get("variant_seq", "")
        if not vid or not seq or len(seq) < args.min_seq_len: continue
        if (o["v1_parent_id"], vid) in seen_vids: continue
        seen_vids.add((o["v1_parent_id"], vid))
        # Map V_ref CDS to variant coords using blocks
        cds_intervals = []
        cds_range = v_ref_cds.get(o["v1_parent_id"])
        if cds_range:
            cs, ce = cds_range
            for b in o["comparison_to_v1"]["blocks"]:
                v1_s, v1_e = b["v1_pos"]
                b_s, b_e = b["between_pos"]
                ov_s = max(v1_s, cs)
                ov_e = min(v1_e, ce)
                if ov_e <= ov_s: continue
                cds_intervals.append((b_s + (ov_s - v1_s), b_s + (ov_e - v1_s)))
        seq_set[f"var_{o['v1_parent_id']}__{vid}"] = {
            "seq": seq, "len": len(seq),
            "kind": o["category"], "v1_parent": o["v1_parent_id"],
            "source_target": o["target"],
            "cds_intervals": cds_intervals,
        }

    print(f"  Total sequences (>= {args.min_seq_len} bp): {len(seq_set):,}",
          file=sys.stderr)

    # Self-align each sequence
    findings = []
    work_root = tempfile.mkdtemp(prefix="self_align_")
    for i, (name, info) in enumerate(seq_set.items()):
        if i % 2000 == 0:
            print(f"  {i:,}/{len(seq_set):,}", file=sys.stderr, flush=True)
        wd = os.path.join(work_root, f"s_{i}")
        os.makedirs(wd, exist_ok=True)
        try:
            hits = self_minimap2(info["seq"], wd, args.threads)
        except Exception as e:
            continue
        repeats = find_tandem_unit(hits, args.min_repeat_len, args.min_identity)
        # Filter to require >= min_copies
        repeats = [r for r in repeats if r["n_copies"] >= args.min_copies]
        if repeats:
            findings.append({
                "seq_id": name,
                "kind": info["kind"],
                "v1_parent": info["v1_parent"],
                "source_target": info["source_target"],
                "seq_len": info["len"],
                "cds_intervals": info.get("cds_intervals", []),
                "tandem_repeats": repeats[:5],  # top 5
            })

    subprocess.run(["rm", "-rf", work_root])

    print(f"\nSequences with tandem repeats: {len(findings):,}", file=sys.stderr)

    with open(os.path.join(args.out, "tandem_findings.json"), "w") as f:
        json.dump(findings, f, indent=2)
    with open(os.path.join(args.out, "tandem_findings.tsv"), "w") as f:
        f.write("seq_id\tkind\tv1_parent\tsource_target\tseq_len\t"
                "repeat_unit_bp\tn_copies\ttotal_span_bp\tspan_start\tspan_end\n")
        for x in findings:
            for r in x["tandem_repeats"]:
                f.write(f"{x['seq_id']}\t{x['kind']}\t{x['v1_parent']}\t"
                        f"{x['source_target']}\t{x['seq_len']}\t"
                        f"{r['repeat_unit_bp']}\t{r['n_copies']}\t"
                        f"{r['total_span_bp']}\t{r['span_start']}\t{r['span_end']}\n")

    # Summary stats
    n_v_ref = sum(1 for x in findings if x["kind"] == "v_ref")
    n_var = sum(1 for x in findings if x["kind"] != "v_ref")
    print(f"  V_refs with tandem: {n_v_ref:,}", file=sys.stderr)
    print(f"  Variants with tandem: {n_var:,}", file=sys.stderr)


if __name__ == "__main__":
    main()
