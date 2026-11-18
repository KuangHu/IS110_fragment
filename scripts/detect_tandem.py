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
import argparse, json, os, shutil, subprocess, sys, tempfile
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
                   help="Min length of a tandem unit (default 50). NOTE: the "
                        "self-alignment (minimap2 asm10, k=19) does not reliably "
                        "see arrays of units shorter than ~150 bp; use TRF/ULTRA "
                        "(verify_tandem_enrichment.py) for short units")
    p.add_argument("--min-identity", type=float, default=95,
                   help="Min alignment identity (default 95)")
    p.add_argument("--min-copies", type=int, default=2,
                   help="Min number of tandem copies to report (default 2 — 1 extra beyond original)")
    p.add_argument("--include-dispersed", action="store_true",
                   help="Also report dispersed direct repeats (two copies with "
                        "other sequence between them). Default: tandem arrays only")
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


def find_tandem_unit(hits, min_repeat_len, min_identity, adjacency=0.8):
    """Group off-diagonal self-hits by offset; classify each repeat.

    A self-hit q[qs,qe) ~ t[ts,te) with offset = ts - qs says the sequence
    repeats itself `offset` bp further on. Whether that is a TANDEM array
    depends on how long the repeated stretch is relative to the offset:

      tandem     the copies are adjacent, so the hit runs for (n-1) units:
                 aligned span >= adjacency * offset. unit = offset,
                 n_copies = round(array_span / offset) with
                 array_span = max_te - min_qs  (n units long, so NO "+1":
                 the old formula reported 3 copies for 2 and so on).
      dispersed  a short repeat (aligned span < adjacency * offset) with
                 unrelated sequence in between: two copies `offset` apart.
                 These were reported as 2-copy tandems with unit = offset
                 (E. coli: 217 of 350 findings).

    Offsets are grouped by closeness (10% or 25 bp), not by rounding to 50:
    the reported unit is the median true offset, not a multiple of 50.

    Returns list of {repeat_class, repeat_unit_bp, n_copies, span_start,
    span_end, total_span_bp, repeat_len_bp, n_hits}.
    """
    # A self-hit and its mirror image (query and target swapped) are the same
    # repeat; minimap2 reports either. Requiring ts > qs threw away arrays
    # whose only report was the mirror. Normalise, then de-duplicate.
    norm = set()
    for (qs, qe, ts, te, strand, ident) in hits:
        if strand != "+" or ident < min_identity or qe - qs < min_repeat_len or ts == qs:
            continue
        if ts < qs:
            qs, qe, ts, te = ts, te, qs, qe
        norm.add((ts - qs, qs, qe, ts, te, round(ident, 2)))
    fwd = sorted(norm)
    groups = []
    for h in fwd:
        if groups and h[0] - groups[-1][-1][0] <= max(25, 0.1 * groups[-1][0][0]):
            groups[-1].append(h)
        else:
            groups.append([h])

    # One offset group can hold several unconnected repeats (two dispersed
    # pairs with the same spacing; two separate arrays of the same unit).
    # Treated as one, their outer span was divided by the offset: two 300-bp
    # pairs 3 kb apart became a 4-copy 1-kb "tandem". Split a group wherever
    # the next hit starts beyond the reach (target end) of the hits so far.
    connected = []
    for g in groups:
        cur, reach = [], None
        for h in sorted(g, key=lambda h: h[1]):
            if reach is not None and h[1] > reach:
                connected.append(cur); cur = []
                reach = None
            cur.append(h)
            reach = max(reach or 0, h[4])
        if cur:
            connected.append(cur)

    results = []
    for g in connected:
        offs = sorted(h[0] for h in g)
        offset = offs[len(offs) // 2]
        if offset < min_repeat_len:
            continue
        min_qs, max_qe = min(h[1] for h in g), max(h[2] for h in g)
        max_te = max(h[4] for h in g)
        aligned = max_qe - min_qs              # length of the repeated stretch
        if aligned >= adjacency * offset:
            cls, n = "tandem", max(2, round((max_te - min_qs) / offset))
        else:
            cls, n = "dispersed", 2
        results.append({
            "repeat_class": cls,
            "repeat_unit_bp": offset if cls == "tandem" else aligned,
            "n_copies": n,
            "span_start": min_qs, "span_end": max_te,
            "total_span_bp": max_te - min_qs,
            "repeat_len_bp": aligned,
            "copy_spacing_bp": offset,
            "identity": sum(h[5] * (h[2] - h[1]) for h in g) / max(1, sum(h[2] - h[1] for h in g)),
            "n_hits": len(g),
        })

    # A tandem array of unit U also self-aligns at offsets 2U, 3U, ...: drop
    # those derivatives (smallest unit = true repeat).
    kept = []
    for r in sorted(results, key=lambda r: r["copy_spacing_bp"]):
        derivative = False
        for k in kept:
            if k["repeat_class"] != "tandem":
                continue
            ratio = r["copy_spacing_bp"] / k["copy_spacing_bp"]
            if round(ratio) < 2 or abs(ratio - round(ratio)) > 0.1:
                continue
            ov = min(r["span_end"], k["span_end"]) - max(r["span_start"], k["span_start"])
            if ov > 0.5 * (r["span_end"] - r["span_start"]):
                derivative = True
                break
        if not derivative:
            kept.append(r)
    kept.sort(key=lambda r: (r["repeat_class"] != "tandem", -r["n_copies"],
                             -r["repeat_unit_bp"], r["span_start"]))
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
                if b.get("strand", "+") == "+":
                    cds_intervals.append((b_s + (ov_s - v1_s), b_s + (ov_e - v1_s)))
                else:      # reverse block: V1 start maps to the block's END
                    cds_intervals.append((b_e - (ov_e - v1_s), b_e - (ov_s - v1_s)))
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
    n_failed = 0
    failed_ids = []
    work_root = tempfile.mkdtemp(prefix="self_align_")
    for i, (name, info) in enumerate(seq_set.items()):
        if i % 2000 == 0:
            print(f"  {i:,}/{len(seq_set):,}", file=sys.stderr, flush=True)
        wd = os.path.join(work_root, f"s_{i}")
        os.makedirs(wd, exist_ok=True)
        try:
            hits = self_minimap2(info["seq"], wd, args.threads)
        except subprocess.CalledProcessError as e:
            # counted and reported; these used to vanish without a trace
            n_failed += 1
            failed_ids.append(name)
            if n_failed <= 5:
                print(f"  WARNING: self-alignment failed for {name}: "
                      f"{(e.stderr or b'')[-200:]!r}", file=sys.stderr)
            continue
        finally:
            shutil.rmtree(wd, ignore_errors=True)
        repeats = find_tandem_unit(hits, args.min_repeat_len, args.min_identity)
        # Filter to require >= min_copies
        repeats = [r for r in repeats if r["n_copies"] >= args.min_copies]
        if not args.include_dispersed:
            repeats = [r for r in repeats if r["repeat_class"] == "tandem"]
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
    failed_path = os.path.join(args.out, "tandem_failed.tsv")

    if n_failed:
        print(f"  WARNING: {n_failed:,} sequences could not be self-aligned",
              file=sys.stderr)
    print(f"\nSequences with tandem repeats: {len(findings):,}", file=sys.stderr)

    with open(os.path.join(args.out, "tandem_findings.json"), "w") as f:
        json.dump(findings, f, indent=2)
    with open(os.path.join(args.out, "tandem_findings.tsv"), "w") as f:
        f.write("seq_id\tkind\tv1_parent\tsource_target\tseq_len\t"
                "repeat_unit_bp\tn_copies\ttotal_span_bp\tspan_start\tspan_end\t"
                "repeat_class\tcopy_spacing_bp\n")
        for x in findings:
            for r in x["tandem_repeats"]:
                f.write(f"{x['seq_id']}\t{x['kind']}\t{x['v1_parent']}\t"
                        f"{x['source_target']}\t{x['seq_len']}\t"
                        f"{r['repeat_unit_bp']}\t{r['n_copies']}\t"
                        f"{r['total_span_bp']}\t{r['span_start']}\t{r['span_end']}\t"
                        f"{r['repeat_class']}\t{r['copy_spacing_bp']}\n")

    # sequences that could not be examined are NOT negatives: listed, and the
    # stage exits non-zero so the pipeline does not mark it complete
    with open(failed_path, "w") as f:
        f.write("seq_id\n" + "".join(x + "\n" for x in failed_ids))
    # Summary stats
    n_v_ref = sum(1 for x in findings if x["kind"] == "v_ref")
    n_var = sum(1 for x in findings if x["kind"] != "v_ref")
    print(f"  V_refs with tandem: {n_v_ref:,}", file=sys.stderr)
    print(f"  Variants with tandem: {n_var:,}", file=sys.stderr)
    if n_failed:
        sys.exit(f"ERROR: {n_failed} sequences could not be self-aligned "
                 f"(see {failed_path}); results are incomplete")


if __name__ == "__main__":
    main()
