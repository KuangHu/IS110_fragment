#!/usr/bin/env python3
"""Lineage visualizer v9: tracks sorted strictly by length.

Difference from v8:
  - V_REF is inserted into the chain at the position determined by its length
    (instead of always being the first non-empty track).
  - Tracks are guaranteed to go shortest → longest.

So a lineage might look like:
   V0 EMPTY
   V1 = 1,391 bp variant
   V2 = 1,477 bp variant
   V3 = V_REF (1,800 bp)
   V4 = 6,968 bp variant

Coloring still uses inheritance via pairwise minimap2 between adjacent tracks.
CDS overlay (red arrow) shown on every track via its V_ref-block mapping.
"""
import argparse, hashlib, json, os, subprocess, sys, tempfile
os.environ["MPLBACKEND"] = "Agg"
from pygenomeviz import GenomeViz
import matplotlib.cm as cm


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--lineages", required=True)
    p.add_argument("--obs", required=True)
    p.add_argument("--records", required=True)
    p.add_argument("--v1-id", required=True)
    p.add_argument("--lineage-idx", type=int, default=0)
    p.add_argument("--out", required=True)
    p.add_argument("--flank-bp", type=int, default=500)
    p.add_argument("--min-cargo-bp", type=int, default=20)
    p.add_argument("--threads", type=int, default=4)
    return p.parse_args()


def load_obs(obs_path, v1_id):
    with open(obs_path) as f:
        obs = json.load(f)
    info = {}
    empties = []
    for o in obs:
        if o["v1_parent_id"] != v1_id: continue
        if o["category"] == "empty":
            empties.append(o)
        vid = o.get("variant_id")
        if vid and vid not in info:
            info[vid] = {
                "target": o["target"],
                "between_start": o["anchor_site"]["between_start"],
                "between_end": o["anchor_site"]["between_end"],
                "strand": o["target_strand"],
                "blocks": o["comparison_to_v1"]["blocks"],
                "variant_seq": o.get("variant_seq", ""),
            }
    return info, empties


def load_v_ref(records_path, v1_id):
    with open(records_path) as f:
        records = json.load(f)
    for r in records:
        rid = r.get("is110_id") or r.get("ref_id")
        if rid != v1_id: continue
        src = r.get("source", {})
        ie = src.get("is_element") or r.get("is_element") or {}
        tn = src.get("transposase_cds") or r.get("transposase_cds") or {}
        is_len = ie.get("length", 0)
        if "start" in ie and "start" in tn:
            strand = ie.get("strand", "+")
            if strand == "+":
                t_s = max(0, tn["start"] - ie["start"])
                t_e = min(is_len, tn["end"] - ie["start"] + 1)
            else:
                t_s = max(0, ie["end"] - tn["end"])
                t_e = min(is_len, ie["end"] - tn["start"] + 1)
        else:
            off5 = ie.get("start_offset_5p", 0)
            tnp_len = tn.get("length", 0)
            t_s = max(0, -off5)
            t_e = min(is_len, t_s + tnp_len)
        return {
            "seq": ie.get("sequence", ""),
            "len": is_len,
            "cds_start": t_s,
            "cds_end": t_e,
            "assembly": src.get("assembly", ""),
            "contig": src.get("contig", ""),
        }
    return None


def palette(n):
    if n <= 10: cmap = cm.get_cmap("tab10")(range(max(2, n)))
    elif n <= 20: cmap = cm.get_cmap("tab20")(range(n))
    else: cmap = cm.get_cmap("hsv")([i / n for i in range(n)])
    return ["#{:02x}{:02x}{:02x}".format(int(c[0]*255), int(c[1]*255), int(c[2]*255))
            for c in cmap]


def pairwise_align(seq_a, seq_b, work_dir, threads):
    fa_q = os.path.join(work_dir, "q.fa")
    fa_t = os.path.join(work_dir, "t.fa")
    paf = os.path.join(work_dir, "qt.paf")
    with open(fa_q, "w") as f: f.write(">q\n" + seq_a + "\n")
    with open(fa_t, "w") as f: f.write(">t\n" + seq_b + "\n")
    subprocess.run(["minimap2", "-x", "asm10", "-c", "--eqx",
                    "-t", str(threads), fa_t, fa_q, "-o", paf],
                   check=True, capture_output=True)
    blocks = []
    with open(paf) as f:
        for line in f:
            c = line.split("\t")
            if len(c) < 12: continue
            qs, qe = int(c[2]), int(c[3])
            strand = c[4]
            ts, te = int(c[7]), int(c[8])
            blocks.append((qs, qe, ts, te, strand))
    return blocks


def lookup_color_at(intervals, pos):
    for s, e, tag in intervals:
        if s <= pos < e: return tag
    return None


def propagate_coloring(prev_intervals, prev_seq, new_seq, work_dir, threads,
                       new_color_tag):
    new_len = len(new_seq)
    blocks = pairwise_align(prev_seq, new_seq, work_dir, threads)
    color_at = [new_color_tag] * new_len
    for (qs, qe, ts, te, strand) in blocks:
        if strand == "+":
            for i in range(te - ts):
                v_n_pos = ts + i
                if v_n_pos >= new_len: break
                v_prev_pos = qs + i
                if v_prev_pos >= len(prev_seq): break
                col = lookup_color_at(prev_intervals, v_prev_pos)
                if col is not None:
                    color_at[v_n_pos] = col
        else:
            for i in range(te - ts):
                v_n_pos = ts + i
                if v_n_pos >= new_len: break
                v_prev_pos = qe - 1 - i
                if v_prev_pos < 0 or v_prev_pos >= len(prev_seq): break
                col = lookup_color_at(prev_intervals, v_prev_pos)
                if col is not None:
                    color_at[v_n_pos] = col
    intervals = []
    if not color_at: return [(0, new_len, new_color_tag)]
    cur_start = 0
    cur_tag = color_at[0]
    for i in range(1, new_len):
        if color_at[i] != cur_tag:
            intervals.append((cur_start, i, cur_tag))
            cur_start = i
            cur_tag = color_at[i]
    intervals.append((cur_start, new_len, cur_tag))
    return intervals


def map_v_ref_to_variant(v_ref_s, v_ref_e, variant_blocks):
    """Map V_ref positions [v_ref_s, v_ref_e] into variant coords via blocks."""
    out = []
    for b in variant_blocks:
        v1_s, v1_e = b["v1_pos"]
        b_s, b_e = b["between_pos"]
        ov_s = max(v1_s, v_ref_s)
        ov_e = min(v1_e, v_ref_e)
        if ov_e <= ov_s: continue
        offset_s = ov_s - v1_s
        offset_e = ov_e - v1_s
        out.append((b_s + offset_s, b_s + offset_e))
    return out


def main():
    args = parse_args()
    with open(args.lineages) as f:
        lineages = json.load(f)
    info = lineages[args.v1_id]
    chain = info["lineages"][args.lineage_idx]

    v_ref = load_v_ref(args.records, args.v1_id)
    if not v_ref or not v_ref["seq"]:
        sys.exit(f"No V_ref sequence for {args.v1_id}")

    obs_info, empty_obs = load_obs(args.obs, args.v1_id)
    empty_present = any(s["type"] == "empty" for s in chain)
    real_steps = [s for s in chain if s["type"] != "empty"]

    # Gather all variants (real_steps) + V_REF itself, then sort by length
    all_variants = []
    for s in real_steps:
        vid = s.get("variant_id")
        if vid == "V_REF": continue  # we'll add V_REF separately
        t = obs_info.get(vid, {})
        all_variants.append({
            "kind": "observed",
            "vid": vid,
            "seq": t.get("variant_seq", ""),
            "len": s["length"],
            "blocks": t.get("blocks", []),
            "target": t.get("target", "?"),
            "pos_start": t.get("between_start", "?"),
            "pos_end": t.get("between_end", "?"),
            "strand": t.get("strand", "?"),
            "n_obs": s["n_observations"],
        })
    # Add V_REF
    all_variants.append({
        "kind": "v_ref",
        "vid": "V_REF",
        "seq": v_ref["seq"],
        "len": v_ref["len"],
        "blocks": [{"v1_pos": [0, v_ref["len"]],
                    "between_pos": [0, v_ref["len"]], "ident": 100.0}],
        "target": f"{v_ref['assembly']}|{v_ref['contig']}",
        "pos_start": "—",
        "pos_end": "—",
        "strand": "+",
        "n_obs": 1,
    })

    # Sort by length
    all_variants.sort(key=lambda v: v["len"])

    # Inheritance coloring: start with shortest, propagate forward
    work_root = tempfile.mkdtemp(prefix="viz_v9_")
    intervals_per_variant = []
    # First variant: simple V_ref-block-based coloring
    first = all_variants[0]
    first_blocks = sorted(first["blocks"], key=lambda b: b["between_pos"][0])
    cur = 0
    cargo_idx = 0
    initial = []
    for b in first_blocks:
        bs, be = b["between_pos"]
        if bs > cur:
            if bs - cur >= args.min_cargo_bp:
                initial.append((cur, bs, f"v1_cargo_{cargo_idx}"))
                cargo_idx += 1
            else:
                initial.append((cur, bs, "v_ref"))
        initial.append((bs, be, "v_ref"))
        cur = max(cur, be)
    if cur < first["len"]:
        if first["len"] - cur >= args.min_cargo_bp:
            initial.append((cur, first["len"], f"v1_cargo_{cargo_idx}"))
        else:
            initial.append((cur, first["len"], "v_ref"))
    if not initial: initial = [(0, first["len"], "v_ref")]
    intervals_per_variant.append(initial)

    # For each subsequent variant: propagate from previous
    for i in range(1, len(all_variants)):
        prev = all_variants[i-1]
        cur_var = all_variants[i]
        wd = os.path.join(work_root, f"v{i}")
        os.makedirs(wd, exist_ok=True)
        new_tag = f"new_at_V{i + (1 if empty_present else 0) + 1}"  # account for V0 empty
        new_intervals = propagate_coloring(
            intervals_per_variant[i-1], prev["seq"], cur_var["seq"],
            wd, args.threads, new_tag)
        intervals_per_variant.append(new_intervals)

    # Color map
    all_tags = []
    seen = set()
    for ivs in intervals_per_variant:
        for _, _, t in ivs:
            if t not in seen:
                seen.add(t)
                all_tags.append(t)
    color_map = {"v_ref": "#aedef5"}
    other = [t for t in all_tags if t != "v_ref"]
    cmap_colors = palette(max(2, len(other)))
    for i, t in enumerate(other):
        color_map[t] = cmap_colors[i % len(cmap_colors)]

    # Build tracks
    labels = []
    track_data = []
    track_idx = 0
    if empty_present:
        ex = empty_obs[0] if empty_obs else None
        tgt = ex["target"] if ex else "?"
        pos = (f"{ex['anchor_site']['between_start']:,}-"
               f"{ex['anchor_site']['between_end']:,}") if ex else ""
        n_e = sum(s["n_observations"] for s in chain if s["type"] == "empty")
        labels.append(f"V{track_idx}  EMPTY ({n_e} obs)\n{tgt} {pos}")
        track_data.append({"kind": "empty", "len": 200, "intervals": None,
                           "cds_intervals": []})
        track_idx += 1

    for v in all_variants:
        if v["kind"] == "v_ref":
            cds_intervals = [(v_ref["cds_start"], v_ref["cds_end"])]
            lbl = (f"V{track_idx}  V_REF ({v['len']:,} bp)  CDS={v_ref['cds_start']}-{v_ref['cds_end']}\n"
                   f"{v['target']}")
        else:
            cds_intervals = map_v_ref_to_variant(v_ref["cds_start"], v_ref["cds_end"],
                                                 v["blocks"])
            lbl = (f"V{track_idx}  ({v['len']:,} bp)  {v['n_obs']} obs\n"
                   f"{v['target']}:{v['pos_start']}-{v['pos_end']} ({v['strand']})")
        labels.append(lbl)
        track_data.append({
            "kind": v["kind"],
            "len": v["len"],
            "intervals": intervals_per_variant[all_variants.index(v)],
            "cds_intervals": cds_intervals,
        })
        track_idx += 1

    max_len = max(t["len"] for t in track_data if t["kind"] != "empty")
    track_len = args.flank_bp + max_len + args.flank_bp

    gv = GenomeViz(fig_width=16, fig_track_height=1.05,
                   track_align_type="left",
                   feature_track_ratio=0.45, link_track_ratio=0.55)

    for i, td in enumerate(track_data):
        track = gv.add_feature_track(labels[i], track_len, labelsize=9)
        flank = args.flank_bp
        track.add_feature(0, flank, plotstyle="bigbox", fc="#dddddd", lw=0.3)
        if td["kind"] == "empty":
            track.add_feature(flank, flank + 200, plotstyle="bigbox",
                              fc="#ff4444", lw=0.5, hatch="//",
                              label="× IS absent",
                              text_kws={"size": 9, "color": "darkred"})
            track.add_feature(flank + 200, flank + 200 + flank,
                              plotstyle="bigbox", fc="#dddddd", lw=0.3)
            continue
        for (s, e, tag) in td["intervals"]:
            color = color_map.get(tag, "#888888")
            label_txt = tag if (tag != "v_ref" and e - s >= 100) else ""
            track.add_feature(flank + s, flank + e, plotstyle="bigbox",
                              fc=color, lw=0.3, label=label_txt,
                              text_kws={"size": 6})
        track.add_feature(flank + td["len"], flank + td["len"] + flank,
                          plotstyle="bigbox", fc="#dddddd", lw=0.3)
        for (cs, ce) in td["cds_intervals"]:
            if ce > cs:
                track.add_feature(flank + cs, flank + ce, plotstyle="bigarrow",
                                  fc="#cc1500", lw=0.8,
                                  label="IS110 CDS" if ce-cs >= 300 else "",
                                  text_kws={"size": 7})

    for i in range(len(track_data) - 1):
        gv.add_link((labels[i], 0, args.flank_bp),
                    (labels[i+1], 0, args.flank_bp), color="lightgrey")
        gv.add_link((labels[i], track_len - args.flank_bp, track_len),
                    (labels[i+1], track_len - args.flank_bp, track_len),
                    color="lightgrey")

    fig = gv.plotfig()
    fig.suptitle(
        f"IS110 lineage at {args.v1_id} (sorted shortest → longest)\n"
        f"V_REF = {v_ref['len']:,} bp (HMM-validated CDS as red arrow), "
        f"now placed by length in the chain",
        fontsize=10)
    fig.savefig(args.out + ".png", dpi=180, bbox_inches="tight")
    fig.savefig(args.out + ".svg", bbox_inches="tight")
    print(f"Saved {args.out}.png ({len(track_data)} tracks)")
    subprocess.run(["rm", "-rf", work_root])


if __name__ == "__main__":
    main()
