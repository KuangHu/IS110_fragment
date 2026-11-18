#!/usr/bin/env python3
"""Lineage builder v3 — nested containment series per anchor site.

For each V1 anchor site, sort all observed IS variants by length (shortest to
longest), with V0=empty at the start. Each step must contain the previous step
with high identity AND coverage (this defines a nested series).

DIRECTION IS NOT INFERRED. A containment series V0 < V1 < V2 is equally
consistent with successive insertions (V0 -> V1 -> V2) and with successive
deletions (V2 -> V1 -> V0): `L-X-R` vs `L-R` cannot tell insertion from
deletion. In the fna project, deciding direction by majority/frequency gave a
different ancestral allele on 55.7% of loci when the outgroup set changed;
only a local-flank genealogy with outgroups made the call stable. So every
lineage is emitted with "direction": "unpolarized", and nothing downstream
may read step order as time order until a polarity test has been run.

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
    p.add_argument("--records", default="",
                   help="Optional records.json — used to attach the real V_ref "
                        "DNA sequence to v_ref steps (otherwise sequence stays empty).")
    p.add_argument("--anchor-identity", type=float, default=95)
    p.add_argument("--cluster-id", type=float, default=0.99,
                   help="cd-hit-est identity for collapsing variants (default 0.99)")
    p.add_argument("--step-identity", type=float, default=95,
                   help="Min %% identity between adjacent steps (default 95)")
    p.add_argument("--step-coverage", type=float, default=95,
                   help="Min %% of shorter step covered by longer step (default 95)")
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


def cluster_variants(seqs, identity, work_dir, min_step_bp=0):
    """cd-hit-est; return mapping seq_id -> cluster_rep_id."""
    if len(seqs) <= 1:
        return {sid: sid for sid, _ in seqs}
    fa = os.path.join(work_dir, "in.fa")
    out_fa = os.path.join(work_dir, "out.fa")
    write_fa(fa, seqs)
    word = 10 if identity >= 0.95 else 8
    # -s/-aS: cd-hit's identity is over the SHORTER sequence, so without a
    # length cutoff a variant fully contained in a longer one (exactly the
    # nested series this script is meant to find) clusters with it and the
    # step disappears. Only near-equal-length sequences may merge.
    r = subprocess.run(["cd-hit-est", "-i", fa, "-o", out_fa,
                        "-c", str(identity), "-n", str(word),
                        "-s", str(identity), "-aS", str(identity),
                        "-T", "0", "-M", "0", "-d", "0"],
                       capture_output=True, text=True)
    if r.returncode != 0:
        # was: silently "no clustering", which changes every count downstream
        raise RuntimeError(f"cd-hit-est failed: {r.stderr.strip()[:500]}")
    # Read whole clusters first: the starred representative is not always
    # the first line, and members listed before it used to be left unmapped.
    clusters, cur = [], None
    with open(out_fa + ".clstr") as f:
        for line in f:
            if line.startswith(">Cluster"):
                cur = {"rep": None, "members": []}
                clusters.append(cur)
                continue
            parts = line.split(">")
            if len(parts) < 2 or cur is None: continue
            sid = parts[1].split("...")[0].strip()
            cur["members"].append(sid)
            if line.rstrip().endswith("*"):
                cur["rep"] = sid
    lens = {sid: len(sq) for sid, sq in seqs}
    mapping = {}
    for c in clusters:
        rep_id = c["rep"] or c["members"][0]
        for sid in c["members"]:
            # a member that differs from the representative by a whole step
            # (>= min_step_bp) is its own variant even at 99% identity:
            # 30,000 bp vs 30,250 bp is a 250-bp step, not a SNP variant
            if min_step_bp and abs(lens[sid] - lens[rep_id]) >= min_step_bp:
                mapping[sid] = sid
            else:
                mapping[sid] = rep_id
    return mapping


def pairwise_minimap2(seqs, work_dir, threads):
    """All-vs-all minimap2. Return dict (q,t) -> {cov_of_q, ident}.

    cov_of_q is the UNION of q's aligned intervals (summing alignment lengths
    let overlapping/repeated hits exceed 100% and fake a containment).
    -N is raised to the number of sequences: with the default (5 secondary
    targets per query) a site with more than six variants lost pairs at
    random, and a missing pair breaks the chain."""
    if len(seqs) < 2:
        return {}
    fa = os.path.join(work_dir, "all.fa")
    write_fa(fa, seqs)
    paf = os.path.join(work_dir, "av.paf")
    subprocess.run(["minimap2", "-x", "asm10", "-c", "--eqx",
                    "-f", "0",  # keep ALL minimizers: the default drops the most frequent
                    # ones, which in a self-alignment are exactly the repeats sought
                    # (fna project: 0 of 60 IS1 copies found without it, 57 with)
                    "-N", str(max(50, len(seqs) + 5)), "-p", "0.2",
                    "-t", str(threads), fa, fa, "-o", paf],
                   check=True, capture_output=True)
    seq_lens = {sid: len(s) for sid, s in seqs}
    ivs = defaultdict(list)
    with open(paf) as f:
        for line in f:
            c = line.split("\t")
            if len(c) < 12: continue
            q, t = c[0], c[5]
            if q == t: continue
            if c[4] != "+":
                continue           # containment is same-strand
            ivs[(q, t)].append((int(c[2]), int(c[3]), int(c[7]), int(c[8]),
                                int(c[9]), int(c[10])))
    out = {}
    for (q, t), iv in ivs.items():
        # Containment needs ONE consistent path: blocks in the same order in
        # both sequences, no stretch of t used twice. Coverage alone said
        # X+Y is "contained" in Y+X, and X+X in a sequence with one X.
        kept = []
        for b in sorted(iv, key=lambda b: -(b[1] - b[0])):
            if all(((b[0] < k[0]) == (b[2] < k[2]))
                   and min(b[3], k[3]) - max(b[2], k[2]) <= 50
                   and min(b[1], k[1]) - max(b[0], k[0]) <= 50 for k in kept):
                kept.append(b)
        covered, cur_s, cur_e = 0, None, None
        for s0, e0 in sorted((b[0], b[1]) for b in kept):
            if cur_e is None or s0 > cur_e:
                if cur_e is not None:
                    covered += cur_e - cur_s
                cur_s, cur_e = s0, e0
            else:
                cur_e = max(cur_e, e0)
        covered += cur_e - cur_s
        out[(q, t)] = {"cov_of_q": covered / max(1, seq_lens.get(q, 1)) * 100,
                       "ident": sum(b[4] for b in kept) / max(1, sum(b[5] for b in kept)) * 100}
    return out


MAX_CHAINS = 500


def build_chains(variant_list, pair_metrics, step_identity, step_coverage):
    """variant_list is sorted shortest-to-longest. Build chains where each
    step n+1 contains step n with sufficient identity + coverage.

    Returns list of chains (each chain is a list of variant indices into variant_list).
    """
    n = len(variant_list)
    if n == 0: return []

    def contains(i_small, i_big):
        m = pair_metrics.get((variant_list[i_small]["vid"], variant_list[i_big]["vid"]))
        return bool(m and m["cov_of_q"] >= step_coverage and m["ident"] >= step_identity)

    chains = [[0]]  # start one chain with smallest
    for i in range(1, n):
        new, seen = [], set()
        for chain in chains:
            # deepest element of this chain that `i` contains. Only the TIP
            # was tested, so with A < B and a variant C that contains A but
            # not B, the relation A < C was lost and C started from nothing.
            j = next((k for k in range(len(chain) - 1, -1, -1)
                      if contains(chain[k], i)), None)
            if j is None:
                continue
            if j == len(chain) - 1:
                chain.append(i)
                seen.add(tuple(chain))
            else:
                branch = tuple(chain[:j + 1] + [i])
                if branch not in seen:
                    seen.add(branch)
                    new.append(list(branch))
        if not seen:
            new.append([i])   # not a child of anything seen so far
        # Branching can double the number of chains at every layer (4,096
        # chains for 24 variants; memory error at 40). Past MAX_CHAINS no new
        # branches are opened: existing chains still grow.
        chains.extend(new[:max(0, MAX_CHAINS - len(chains))] if seen else new)
    # drop chains that are a prefix of another (left behind by branching)
    uniq = []
    for c in sorted(set(map(tuple, chains)), key=lambda c: (-len(c), c)):
        if not any(u[:len(c)] == c for u in uniq):
            uniq.append(c)
    return [list(c) for c in sorted(uniq)]


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)
    work_root = tempfile.mkdtemp(prefix="lineage_v3_")
    import atexit, shutil
    atexit.register(shutil.rmtree, work_root, ignore_errors=True)
    print(f"Loading {args.obs}...", file=sys.stderr, flush=True)
    with open(args.obs) as f:
        obs = json.load(f)
    print(f"  {len(obs):,} observations", file=sys.stderr)

    v_ref_seqs = {}      # ref_id -> IS-element nt sequence (for v_ref steps)
    v_upstreams = {}     # ref_id -> upstream flank sequence (shared by every step in the chain)
    v_downstreams = {}   # ref_id -> downstream flank sequence
    if args.records and os.path.exists(args.records):
        print(f"Loading V_ref + flank sequences from {args.records}...", file=sys.stderr, flush=True)
        with open(args.records) as f:
            recs = json.load(f)
        for r in recs:
            rid = r.get("ref_id") or r.get("is110_id")
            src = r.get("source", {})
            ie = r.get("is_element") or src.get("is_element") or {}
            up = r.get("upstream_flank") or src.get("upstream_flank") or {}
            dn = r.get("downstream_flank") or src.get("downstream_flank") or {}
            if not rid:
                continue
            if ie.get("sequence"): v_ref_seqs[rid] = ie["sequence"]
            if up.get("sequence"): v_upstreams[rid] = up["sequence"]
            if dn.get("sequence"): v_downstreams[rid] = dn["sequence"]
        print(f"  V_ref/flanks loaded: ref={len(v_ref_seqs):,}, "
              f"up={len(v_upstreams):,}, down={len(v_downstreams):,}",
              file=sys.stderr)

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
                    "nested_event_keys": set(),
                    "v_ref_cov_pct": round(
                        o["comparison_to_v1"]["total_matched_bp"] / v1_len * 100, 1
                    ),
                }
            variants_raw[vid]["n_observations"] += 1
            ni = o.get("nested_insert")
            if ni:
                variants_raw[vid]["nested_event_keys"].add(ni["event_key"])

        # Add V_ref itself as a "variant" entry (if observed clonally)
        v_ref_seq = None
        # Find V_ref seq from records.json? It's not directly in obs, but the
        # V_ref blocks correspond to V1's content. We can extract V_ref seq from
        # the clonal_V1 observation if available, otherwise from the matched blocks
        # of a sufficiently-covered insertion variant.
        # Easier: just record V_ref as a placeholder with length=v1_len and no seq.
        if n_clonal > 0:
            variants_raw["V_REF"] = {
                "vid": "V_REF",
                "seq": v_ref_seqs.get(v1_id, ""),
                "len": v1_len,
                "n_observations": n_clonal, "example_target": v1_id,
                "v_ref_cov_pct": 100.0,
            }

        if not variants_raw and n_empty == 0:
            continue

        # Cluster at 99% identity
        if len(variants_raw) > 1:
            work_dir = os.path.join(work_root, v1_id.replace("|", "_").replace(".", "_"))
            os.makedirs(work_dir, exist_ok=True)
            # Exclude V_REF from clustering — it represents the canonical IS and
            # must stay as its own entry, even though it now carries a sequence.
            seqs = [(v["vid"], v["seq"]) for v in variants_raw.values()
                    if v["seq"] and v["vid"] != "V_REF"]
            if seqs:
                mapping = cluster_variants(seqs, args.cluster_id, work_dir,
                                           min_step_bp=args.min_step_bp)
                aggregated = {}
                for vid, v in variants_raw.items():
                    if vid == "V_REF" or not v["seq"]:
                        # V_REF kept as-is (even with sequence populated)
                        aggregated[vid] = v.copy()
                        aggregated[vid]["cluster_members"] = [vid]
                        continue
                    rep = mapping.get(vid, vid)
                    if rep not in aggregated:
                        # the REPRESENTATIVE's own record (was: whichever
                        # member came first, filed under the rep's key)
                        aggregated[rep] = variants_raw.get(rep, v).copy()
                        aggregated[rep]["cluster_members"] = []
                        aggregated[rep]["n_observations"] = 0
                        aggregated[rep]["nested_event_keys"] = set()
                    aggregated[rep]["cluster_members"].append(vid)
                    aggregated[rep]["nested_event_keys"] = (
                        set(aggregated[rep].get("nested_event_keys", ()))
                        | v.get("nested_event_keys", set()))
                    aggregated[rep]["n_observations"] += v["n_observations"]
                variants_raw = aggregated

        # Sort variants shortest to longest
        variants_sorted = sorted(variants_raw.values(), key=lambda v: (v["len"], v["vid"]))

        # Pairwise minimap2 (containment + similarity), before any merging
        seqs_for_align = [(v["vid"], v["seq"]) for v in variants_sorted if v["seq"]]
        if len(seqs_for_align) >= 2:
            work_dir = os.path.join(work_root, v1_id.replace("|", "_").replace(".", "_"))
            os.makedirs(work_dir, exist_ok=True)
            pair_metrics = pairwise_minimap2(seqs_for_align, work_dir, args.threads)
        else:
            pair_metrics = {}

        # Length-bucket: collapse variants within --min-step-bp of each other
        # into one representative (most observations) -- but only if they are
        # also the same sequence. Length alone merged unrelated variants that
        # happen to be equally long (two different cargos of the same size)
        # and added their observation counts together.
        def same_seq(a, b):
            m1 = pair_metrics.get((a["vid"], b["vid"]))
            m2 = pair_metrics.get((b["vid"], a["vid"]))
            return bool(m1 and m2 and min(m1["cov_of_q"], m2["cov_of_q"]) >= args.step_coverage
                        and min(m1["ident"], m2["ident"]) >= args.step_identity)

        if args.min_step_bp > 0 and len(variants_sorted) > 1:
            buckets = []
            for v in variants_sorted:
                # similar to EVERY member (complete linkage): similarity to
                # the first member alone chained B and C together through A
                home = next((bk for bk in reversed(buckets)
                             if v["len"] - bk[0]["len"] < args.min_step_bp
                             and all(same_seq(m, v) for m in bk)), None)
                if home is not None:
                    home.append(v)
                else:
                    buckets.append([v])
            collapsed = []
            for bucket in buckets:
                bucket.sort(key=lambda v: (-v["n_observations"], -v["len"], v["vid"]))
                rep = bucket[0].copy()
                rep["n_observations"] = sum(v["n_observations"] for v in bucket)
                rep["nested_event_keys"] = set().union(
                    *(v.get("nested_event_keys", set()) for v in bucket))
                rep["bucket_size"] = len(bucket)
                rep["bucket_lengths"] = [v["len"] for v in bucket]
                collapsed.append(rep)
            variants_sorted = sorted(collapsed, key=lambda v: (v["len"], v["vid"]))

        # Build chains
        chains = build_chains(variants_sorted, pair_metrics,
                              args.step_identity, args.step_coverage)

        # flanks for this anchor site (shared across every step of every chain)
        upstream_seq   = v_upstreams.get(v1_id, "")
        downstream_seq = v_downstreams.get(v1_id, "")

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
                    "upstream":   upstream_seq,
                    "insertion":  "",
                    "downstream": downstream_seq,
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
                    "nested_event_keys": sorted(v.get("nested_event_keys", ())),
                    "identity_to_prev": id_to_prev,
                    "coverage_to_prev": cov_to_prev,
                    "upstream":   upstream_seq,
                    "insertion":  v.get("seq", ""),
                    "downstream": downstream_seq,
                })
            if len(steps) >= 2:
                lineages.append(steps)

        if not lineages: continue
        if any(len(l) >= 3 for l in lineages): n_multistep += 1

        out_lineages[v1_id] = {
            "v1_id": v1_id,
            "v_ref_length": v1_len,
            "n_lineages": len(lineages),
            "chains_truncated": len(chains) >= MAX_CHAINS,
            "direction": "unpolarized",
            "direction_note": "steps are ordered by length and nested by "
                              "containment; this is not a time order",
            "counts_are": "distinct assemblies (GCA/GCF twins collapsed)",
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
                "identity_to_prev\tcoverage_to_prev\texample_target\tdirection\n")
        for v1_id, info in out_lineages.items():
            for li, lineage in enumerate(info["lineages"]):
                for s in lineage:
                    f.write(f"{v1_id}\t{info['v_ref_length']}\t{li}\t"
                            f"{s['step_index']}\t{s['label']}\t{s['type']}\t"
                            f"{s['length']}\t{s.get('variant_id', '')}\t"
                            f"{s['n_observations']}\t{s.get('v_ref_cov_pct', '')}\t"
                            f"{s.get('identity_to_prev', '')}\t"
                            f"{s.get('coverage_to_prev', '')}\t"
                            f"{s.get('example_target', '')}\t"
                            f"{info['direction']}\n")

    print(f"\nSaved: {args.out}/lineages.json", file=sys.stderr)
    print(f"Saved: {args.out}/lineages.tsv", file=sys.stderr)
    subprocess.run(["rm", "-rf", work_root])


if __name__ == "__main__":
    main()
