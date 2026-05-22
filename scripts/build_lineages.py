#!/usr/bin/env python3
"""Lineage builder v3 — corrected framing: every step is an insertion.

For each V1 anchor site, sort all observed IS variants by length (shortest to
longest), with V0=empty at the start. Each step must contain the previous step
with high identity AND coverage (this defines a true nested insertion).

Output per V1:
  [V0 empty,  V1 shortest,  V2 = V1 + insertion,  V3 = V2 + insertion, ...]

Each adjacent step is validated:
  - identity_to_prev_step >= 95%
  - prev_step is >= 80% covered by current step

If validation fails between V_n and V_(n+1), the chain splits.

Inputs:
  observations.json  (from downstream_finder.py)

Outputs:
  lineages.json  — per V1 ref, list of validated lineage chains
  lineages.tsv   — flat table with step-by-step coverage/identity numbers
"""
import argparse, hashlib, json, os, subprocess, sys, tempfile
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--obs", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--anchor-identity", type=float, default=95)
    p.add_argument("--cluster-id", type=float, default=0.99,
                   help="cd-hit-est identity for collapsing variants (default 0.99)")
    p.add_argument("--step-identity", type=float, default=95,
                   help="Min %% identity between adjacent steps (default 95)")
    p.add_argument("--step-coverage", type=float, default=80,
                   help="Min %% of shorter step covered by longer step (default 80)")
    p.add_argument("--min-step-bp", type=int, default=200,
                   help="Min length difference between adjacent steps (smaller diffs "
                        "are treated as the same step, just SNP variants; default 200)")
    p.add_argument("--threads", type=int, default=8)
    return p.parse_args()


def short_hash(s):
    return hashlib.md5(s.encode()).hexdigest()[:10]


def write_fa(path, seqs):
    with open(path, "w") as f:
        for sid, seq in seqs:
            f.write(f">{sid}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")


def cluster_variants(seqs, identity, work_dir):
    """cd-hit-est; return mapping seq_id -> cluster_rep_id."""
    if len(seqs) <= 1:
        return {sid: sid for sid, _ in seqs}
    fa = os.path.join(work_dir, "in.fa")
    out_fa = os.path.join(work_dir, "out.fa")
    write_fa(fa, seqs)
    word = 10 if identity >= 0.95 else 8
    r = subprocess.run(["cd-hit-est", "-i", fa, "-o", out_fa,
                        "-c", str(identity), "-n", str(word),
                        "-T", "0", "-M", "0", "-d", "0"],
                       capture_output=True, text=True)
    if r.returncode != 0:
        return {sid: sid for sid, _ in seqs}
    mapping = {}
    rep = None
    with open(out_fa + ".clstr") as f:
        for line in f:
            if line.startswith(">Cluster"):
                rep = None
                continue
            parts = line.split(">")
            if len(parts) < 2: continue
            sid = parts[1].split("...")[0].strip()
            if "*" in line:
                rep = sid
                mapping[sid] = sid
            elif rep:
                mapping[sid] = rep
    return mapping


def pairwise_minimap2(seqs, work_dir, threads):
    """All-vs-all minimap2. Return dict (q,t) -> {cov_of_q, ident}."""
    if len(seqs) < 2:
        return {}
    fa = os.path.join(work_dir, "all.fa")
    write_fa(fa, seqs)
    paf = os.path.join(work_dir, "av.paf")
    subprocess.run(["minimap2", "-x", "asm10", "-c", "--eqx",
                    "-t", str(threads), fa, fa, "-o", paf],
                   check=True, capture_output=True)
    seq_lens = {sid: len(s) for sid, s in seqs}
    match_bp = defaultdict(int)
    ident_sum = defaultdict(float)
    ident_count = defaultdict(int)
    with open(paf) as f:
        for line in f:
            c = line.split("\t")
            if len(c) < 12: continue
            q, t = c[0], c[5]
            if q == t: continue
            qs, qe = int(c[2]), int(c[3])
            matches, block = int(c[9]), int(c[10])
            ident = matches / block * 100 if block > 0 else 0
            match_bp[(q, t)] += qe - qs
            ident_sum[(q, t)] += ident * (qe - qs)
            ident_count[(q, t)] += qe - qs
    out = {}
    for (q, t), m in match_bp.items():
        qlen = seq_lens.get(q, 1)
        cov_q_in_t = m / qlen * 100
        avg_id = ident_sum[(q, t)] / ident_count[(q, t)] if ident_count[(q, t)] > 0 else 0
        out[(q, t)] = {"cov_of_q": cov_q_in_t, "ident": avg_id}
    return out


def build_chains(variant_list, pair_metrics, step_identity, step_coverage):
    """variant_list is sorted shortest-to-longest. Build chains where each
    step n+1 contains step n with sufficient identity + coverage.

    Returns list of chains (each chain is a list of variant indices into variant_list).
    """
    n = len(variant_list)
    if n == 0: return []
    chains = [[0]]  # start one chain with smallest
    for i in range(1, n):
        cur = variant_list[i]
        # Try to extend any existing chain whose last element is contained in cur
        extended_any = False
        for chain in chains:
            prev = variant_list[chain[-1]]
            metrics = pair_metrics.get((prev["vid"], cur["vid"]))
            if metrics and metrics["cov_of_q"] >= step_coverage \
                       and metrics["ident"] >= step_identity:
                chain.append(i)
                extended_any = True
        if not extended_any:
            # Start a new chain (this variant isn't a child of any existing chain)
            chains.append([i])
    return chains


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)
    work_root = tempfile.mkdtemp(prefix="lineage_v3_")
    print(f"Loading {args.obs}...", file=sys.stderr, flush=True)
    with open(args.obs) as f:
        obs = json.load(f)
    print(f"  {len(obs):,} observations", file=sys.stderr)

    per_v1 = defaultdict(list)
    for o in obs:
        per_v1[o["v1_parent_id"]].append(o)
    print(f"  {len(per_v1):,} V1 references", file=sys.stderr)

    out_lineages = {}
    n_v1 = 0
    n_multistep = 0

    for v1_id, observations in per_v1.items():
        n_v1 += 1
        if n_v1 % 1000 == 0:
            print(f"    {n_v1:,}/{len(per_v1):,}", file=sys.stderr, flush=True)

        v1_len = observations[0]["comparison_to_v1"]["v1_len"]
        observations = [o for o in observations
                        if o["anchor_site"]["up_anchor_identity"] >= args.anchor_identity
                        and o["anchor_site"]["down_anchor_identity"] >= args.anchor_identity]
        if not observations: continue

        n_empty = sum(1 for o in observations if o["category"] == "empty")
        n_clonal = sum(1 for o in observations if o["category"] == "clonal_V1")

        # Collect all variant sequences (no category distinction now)
        # Each unique variant_id -> {seq, len, n_obs, example_target}
        variants_raw = {}
        for o in observations:
            cat = o["category"]
            if cat not in ("deletion", "insertion"): continue
            vid = o.get("variant_id")
            seq = o.get("variant_seq", "")
            if not vid or not seq: continue
            if vid not in variants_raw:
                variants_raw[vid] = {
                    "vid": vid, "seq": seq, "len": o.get("variant_len", len(seq)),
                    "n_observations": 0, "example_target": o["target"],
                    "v_ref_cov_pct": round(
                        o["comparison_to_v1"]["total_matched_bp"] / v1_len * 100, 1
                    ),
                }
            variants_raw[vid]["n_observations"] += 1

        # Add V_ref itself as a "variant" entry (if observed clonally)
        v_ref_seq = None
        # Find V_ref seq from records.json? It's not directly in obs, but the
        # V_ref blocks correspond to V1's content. We can extract V_ref seq from
        # the clonal_V1 observation if available, otherwise from the matched blocks
        # of a sufficiently-covered insertion variant.
        # Easier: just record V_ref as a placeholder with length=v1_len and no seq.
        if n_clonal > 0:
            variants_raw["V_REF"] = {
                "vid": "V_REF", "seq": "", "len": v1_len,
                "n_observations": n_clonal, "example_target": v1_id,
                "v_ref_cov_pct": 100.0,
            }

        if not variants_raw and n_empty == 0:
            continue

        # Cluster at 99% identity
        if len(variants_raw) > 1:
            work_dir = os.path.join(work_root, v1_id.replace("|", "_").replace(".", "_"))
            os.makedirs(work_dir, exist_ok=True)
            seqs = [(v["vid"], v["seq"]) for v in variants_raw.values() if v["seq"]]
            # Skip V_REF in clustering since it has no sequence
            if seqs:
                mapping = cluster_variants(seqs, args.cluster_id, work_dir)
                aggregated = {}
                for vid, v in variants_raw.items():
                    if not v["seq"]:
                        # V_REF kept as-is
                        aggregated[vid] = v.copy()
                        aggregated[vid]["cluster_members"] = [vid]
                        continue
                    rep = mapping.get(vid, vid)
                    if rep not in aggregated:
                        aggregated[rep] = v.copy()
                        aggregated[rep]["cluster_members"] = []
                        aggregated[rep]["n_observations"] = 0
                    aggregated[rep]["cluster_members"].append(vid)
                    aggregated[rep]["n_observations"] += v["n_observations"]
                variants_raw = aggregated

        # Sort variants shortest to longest
        variants_sorted = sorted(variants_raw.values(), key=lambda v: v["len"])

        # Length-bucket: collapse variants within --min-step-bp of each other into
        # a single representative (the one with the most observations)
        if args.min_step_bp > 0 and len(variants_sorted) > 1:
            buckets = []
            for v in variants_sorted:
                if buckets and v["len"] - buckets[-1][-1]["len"] < args.min_step_bp:
                    buckets[-1].append(v)
                else:
                    buckets.append([v])
            collapsed = []
            for bucket in buckets:
                # Pick the representative with most observations, or longest
                bucket.sort(key=lambda v: (-v["n_observations"], -v["len"]))
                rep = bucket[0].copy()
                rep["n_observations"] = sum(v["n_observations"] for v in bucket)
                rep["bucket_size"] = len(bucket)
                rep["bucket_lengths"] = [v["len"] for v in bucket]
                collapsed.append(rep)
            variants_sorted = collapsed

        # Pairwise minimap2 to validate containment
        seqs_for_align = [(v["vid"], v["seq"]) for v in variants_sorted if v["seq"]]
        if len(seqs_for_align) >= 2:
            work_dir = os.path.join(work_root, v1_id.replace("|", "_").replace(".", "_"))
            os.makedirs(work_dir, exist_ok=True)
            pair_metrics = pairwise_minimap2(seqs_for_align, work_dir, args.threads)
        else:
            pair_metrics = {}

        # Build chains
        chains = build_chains(variants_sorted, pair_metrics,
                              args.step_identity, args.step_coverage)

        # Format each chain with the empty prefix
        lineages = []
        for chain in chains:
            steps = []
            if n_empty > 0:
                steps.append({
                    "step_index": 0, "label": "V0",
                    "type": "empty", "length": 0,
                    "n_observations": n_empty,
                    "identity_to_prev": None,
                    "coverage_to_prev": None,
                })
            for i, idx in enumerate(chain):
                v = variants_sorted[idx]
                # identity / coverage relative to previous step
                if i == 0:
                    id_to_prev = None
                    cov_to_prev = None
                else:
                    prev_v = variants_sorted[chain[i - 1]]
                    metrics = pair_metrics.get((prev_v["vid"], v["vid"]), {})
                    id_to_prev = round(metrics.get("ident", 0), 1) if metrics else None
                    cov_to_prev = round(metrics.get("cov_of_q", 0), 1) if metrics else None
                step_idx = len(steps)
                steps.append({
                    "step_index": step_idx,
                    "label": f"V{step_idx}",
                    "type": "v_ref" if v["vid"] == "V_REF" else "variant",
                    "length": v["len"],
                    "variant_id": v["vid"],
                    "n_observations": v["n_observations"],
                    "v_ref_cov_pct": v["v_ref_cov_pct"],
                    "example_target": v["example_target"],
                    "identity_to_prev": id_to_prev,
                    "coverage_to_prev": cov_to_prev,
                })
            if len(steps) >= 2:
                lineages.append(steps)

        if not lineages: continue
        if any(len(l) >= 3 for l in lineages): n_multistep += 1

        out_lineages[v1_id] = {
            "v1_id": v1_id,
            "v_ref_length": v1_len,
            "n_lineages": len(lineages),
            "lineages": lineages,
        }

    print(f"\n  V1 refs with at least one lineage: {len(out_lineages):,}",
          file=sys.stderr)
    print(f"  V1 refs with a lineage of >= 3 steps: {n_multistep:,}",
          file=sys.stderr)

    with open(os.path.join(args.out, "lineages.json"), "w") as f:
        json.dump(out_lineages, f, indent=2)

    with open(os.path.join(args.out, "lineages.tsv"), "w") as f:
        f.write("v1_id\tv_ref_length\tlineage_idx\tstep\tlabel\ttype\tlength\t"
                "variant_id\tn_observations\tv_ref_cov_pct\t"
                "identity_to_prev\tcoverage_to_prev\texample_target\n")
        for v1_id, info in out_lineages.items():
            for li, lineage in enumerate(info["lineages"]):
                for s in lineage:
                    f.write(f"{v1_id}\t{info['v_ref_length']}\t{li}\t"
                            f"{s['step_index']}\t{s['label']}\t{s['type']}\t"
                            f"{s['length']}\t{s.get('variant_id', '')}\t"
                            f"{s['n_observations']}\t{s.get('v_ref_cov_pct', '')}\t"
                            f"{s.get('identity_to_prev', '')}\t"
                            f"{s.get('coverage_to_prev', '')}\t"
                            f"{s.get('example_target', '')}\n")

    print(f"\nSaved: {args.out}/lineages.json", file=sys.stderr)
    print(f"Saved: {args.out}/lineages.tsv", file=sys.stderr)
    subprocess.run(["rm", "-rf", work_root])


if __name__ == "__main__":
    main()
