#!/usr/bin/env python3
"""Validate V_0 (shortest variant) by checking V_ref's transposase CDS is covered.

Uses already-computed data:
  - records_final/records.json: gives V_ref's transposase position
  - observations.json: gives V_0's alignment blocks to V_ref

For each lineage, check:
  - V_ref's transposase region [tnp_s, tnp_e] (in V_ref coords)
  - V_0's alignment blocks (mapping V_0 -> V_ref positions)
  - PASS if the transposase region is covered by V_0's alignment blocks at
    >= --min-cds-coverage% (default 95%)

V_0 = V_REF (clonal observation) is automatically PASS (HMM already validated upstream).

Usage:
    validate_v0_cds_coverage.py --lineages lineages.json --obs observations.json \\
        --records records.json --out out_dir
"""
import argparse, json, os, sys


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--lineages", required=True)
    p.add_argument("--obs", required=True)
    p.add_argument("--records", required=True,
                   help="records_final/records.json with V_ref CDS positions")
    p.add_argument("--out", required=True)
    p.add_argument("--min-cds-coverage", type=float, default=95,
                   help="Min %% of V_ref's transposase that must be covered in V_0 (default 95)")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    print(f"Loading records_final/records.json...", file=sys.stderr, flush=True)
    with open(args.records) as f:
        records = json.load(f)

    # Build dict: v1_id -> {tnp_start_in_v_ref, tnp_end_in_v_ref, v_ref_len}
    # Two record schemas supported:
    #   (A) records_final/records.json: source.{is_element, transposase_cds} with genomic coords
    #   (B) build_records.py (Cross_reference_IS Stage 6): top-level is_element + transposase_cds
    tnp_pos = {}
    for r in records:
        v1_id = r.get("is110_id") or r.get("ref_id")
        if not v1_id:
            continue
        src = r.get("source", {})
        # is_element + transposase may be at top level OR under source
        is_elem = src.get("is_element") or r.get("is_element") or {}
        tnp = src.get("transposase_cds") or r.get("transposase_cds") or {}
        if not is_elem or not tnp:
            continue
        is_len = is_elem.get("length", 0)
        if is_len <= 0:
            continue
        # Schema A: is_element has start/end/strand (genomic), transposase has start/end (genomic)
        if "start" in is_elem and "start" in tnp:
            is_start = is_elem["start"]
            is_end = is_elem.get("end", is_start + is_len)
            strand = is_elem.get("strand", "+")
            tnp_start = tnp["start"]
            tnp_end = tnp["end"]
            if strand == "+":
                t_s_in_v = max(0, tnp_start - is_start)
                t_e_in_v = min(is_len, tnp_end - is_start + 1)
            else:
                t_s_in_v = max(0, is_end - tnp_end)
                t_e_in_v = min(is_len, is_end - tnp_start + 1)
        else:
            # Schema B: derive from start_offset_5p + transposase length
            #   transposase sits at position [-start_offset_5p, -start_offset_5p + tnp_len]
            #   in V_ref-internal coords (start_offset_5p is negative offset from tnp start)
            off5 = is_elem.get("start_offset_5p", 0)
            tnp_len = tnp.get("length", 0)
            t_s_in_v = max(0, -off5)
            t_e_in_v = min(is_len, t_s_in_v + tnp_len)
        if t_e_in_v <= t_s_in_v:
            continue
        tnp_pos[v1_id] = {
            "tnp_start_in_v_ref": t_s_in_v,
            "tnp_end_in_v_ref": t_e_in_v,
            "v_ref_len": is_len,
            "tnp_len": t_e_in_v - t_s_in_v,
        }
    print(f"  V_ref entries with CDS info: {len(tnp_pos):,}", file=sys.stderr)

    print(f"Loading observations...", file=sys.stderr, flush=True)
    with open(args.obs) as f:
        obs = json.load(f)
    # Index blocks by (v1_id, variant_id)
    blocks_by_vid = {}
    for o in obs:
        vid = o.get("variant_id")
        if not vid: continue
        key = (o["v1_parent_id"], vid)
        if key not in blocks_by_vid:
            blocks_by_vid[key] = o["comparison_to_v1"]["blocks"]

    print(f"Loading lineages...", file=sys.stderr, flush=True)
    with open(args.lineages) as f:
        lineages = json.load(f)

    # For each lineage chain, check V_0
    n_pass = 0
    n_fail = 0
    n_no_info = 0
    report = []
    filtered = {}
    for v1_id, info in lineages.items():
        tnp_info = tnp_pos.get(v1_id)
        if not tnp_info:
            n_no_info += len(info["lineages"])
            continue
        t_s = tnp_info["tnp_start_in_v_ref"]
        t_e = tnp_info["tnp_end_in_v_ref"]
        tnp_total = t_e - t_s
        if tnp_total <= 0:
            n_no_info += len(info["lineages"])
            continue

        for li, chain in enumerate(info["lineages"]):
            v0_step = next((s for s in chain if s["type"] != "empty"), None)
            if not v0_step: continue
            vid = v0_step.get("variant_id")
            if not vid: continue

            # V_REF clonal case → automatic pass
            if v0_step["type"] == "v_ref" or vid == "V_REF":
                cds_covered_pct = 100.0
                passed = True
            else:
                # Look up V_0's blocks
                blocks = blocks_by_vid.get((v1_id, vid), [])
                if not blocks:
                    n_no_info += 1
                    continue
                # Compute total V_ref-coverage within tnp region
                # Each block has 'v1_pos': [start, end]
                covered_bp = 0
                for b in blocks:
                    bs, be = b["v1_pos"]
                    # Intersect with [t_s, t_e]
                    ov_s = max(bs, t_s)
                    ov_e = min(be, t_e)
                    if ov_e > ov_s:
                        covered_bp += ov_e - ov_s
                cds_covered_pct = covered_bp / tnp_total * 100
                passed = cds_covered_pct >= args.min_cds_coverage

            report.append({
                "v1_id": v1_id,
                "lineage_idx": li,
                "v0_variant_id": vid,
                "v0_length": v0_step["length"],
                "tnp_total_bp": tnp_total,
                "cds_coverage_pct": round(cds_covered_pct, 1),
                "passes": passed,
            })
            if passed:
                n_pass += 1
                if v1_id not in filtered:
                    filtered[v1_id] = {
                        "v1_id": v1_id,
                        "v_ref_length": info["v_ref_length"],
                        "lineages": [],
                    }
                filtered[v1_id]["lineages"].append(chain)
            else:
                n_fail += 1

    print(f"\nResults:", file=sys.stderr)
    print(f"  Lineages with V_0 PASSING CDS coverage check: {n_pass:,}",
          file=sys.stderr)
    print(f"  Lineages with V_0 FAILING (CDS too truncated): {n_fail:,}",
          file=sys.stderr)
    print(f"  Lineages with no CDS/block info available:    {n_no_info:,}",
          file=sys.stderr)
    print(f"  V_1 refs with at least one passing lineage:   {len(filtered):,}",
          file=sys.stderr)

    with open(os.path.join(args.out, "lineages_validated.json"), "w") as f:
        json.dump(filtered, f, indent=2)
    with open(os.path.join(args.out, "v0_cds_report.tsv"), "w") as f:
        f.write("v1_id\tlineage_idx\tv0_variant_id\tv0_length\t"
                "tnp_total_bp\tcds_coverage_pct\tpasses\n")
        for r in report:
            f.write(f"{r['v1_id']}\t{r['lineage_idx']}\t{r['v0_variant_id']}\t"
                    f"{r['v0_length']}\t{r['tnp_total_bp']}\t"
                    f"{r['cds_coverage_pct']}\t{r['passes']}\n")
    print(f"\nSaved: {args.out}/lineages_validated.json", file=sys.stderr)
    print(f"Saved: {args.out}/v0_cds_report.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()
