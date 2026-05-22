#!/usr/bin/env python3
"""Visualize all sequences with tandem repeats.

Each sequence becomes one figure:
  - A horizontal bar showing the sequence
  - Non-repeat regions in light blue (or gray)
  - Each tandem repeat region in a distinct color (per repeat unit)
  - Individual repeat units within the tandem shown as separate adjacent blocks
  - Labels show accession + sequence length + repeat unit size + copy count

For batch output: one PNG per sequence in --out directory.
"""
import argparse, json, os, sys
os.environ["MPLBACKEND"] = "Agg"
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.cm as cm


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--findings", required=True,
                   help="tandem_findings.json from self_align_tandem.py")
    p.add_argument("--out", required=True, help="Output directory")
    p.add_argument("--max-figures", type=int, default=400,
                   help="Limit number of figures (default 400 = enough for ~371 findings)")
    return p.parse_args()


def palette(n):
    if n <= 10: cmap = cm.get_cmap("tab10")(range(max(2, n)))
    elif n <= 20: cmap = cm.get_cmap("tab20")(range(n))
    else: cmap = cm.get_cmap("hsv")([i / n for i in range(n)])
    return [(c[0], c[1], c[2]) for c in cmap]


def draw_one(rec, out_png):
    """Render a single sequence with its tandem repeat highlighted."""
    seq_len = rec["seq_len"]
    repeats = rec["tandem_repeats"]
    fig_w = 14
    fig_h = 2.2
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # Pick colors for each repeat finding (top 3 most prominent)
    rep_colors = palette(max(2, len(repeats)))

    # Backbone first
    ax.add_patch(patches.Rectangle((0, 0.45), seq_len, 0.1,
                                    facecolor="#aedef5",
                                    edgecolor="#333", lw=0.5,
                                    label="sequence body"))

    # Overlay each repeat region's tandem units
    for ri, r in enumerate(repeats[:3]):
        span_s = r["span_start"]
        span_e = r["span_end"]
        unit_bp = r["repeat_unit_bp"]
        n_copies = r["n_copies"]
        color = rep_colors[ri]
        # Draw individual tandem units (best-effort)
        # We don't know exact unit boundaries, so estimate from span/n_copies
        unit_len = max(1, (span_e - span_s) / n_copies)
        for k in range(n_copies):
            us = span_s + k * unit_len
            ue = min(span_e, us + unit_len)
            ax.add_patch(patches.Rectangle((us, 0.42 - ri * 0.05), ue - us, 0.16,
                                            facecolor=color, edgecolor="#222",
                                            lw=0.3, alpha=0.85))

    # CDS overlay (red arrows) — drawn ABOVE the backbone
    cds_intervals = rec.get("cds_intervals", [])
    for (cs, ce) in cds_intervals:
        if ce > cs:
            # Draw an arrow shape pointing right
            ax.add_patch(patches.FancyArrow(
                cs, 0.75, ce - cs, 0,
                width=0.18,
                head_width=0.28,
                head_length=min(seq_len * 0.015, (ce - cs) * 0.3),
                length_includes_head=True,
                facecolor="#cc1500", edgecolor="#330000", lw=0.5))
            # Label
            mid = (cs + ce) / 2
            ax.annotate("IS110 CDS", (mid, 1.05), ha="center",
                        fontsize=7, color="darkred")

    # Axis settings
    ax.set_xlim(-seq_len * 0.02, seq_len * 1.02)
    ax.set_ylim(0, 1.3)
    ax.set_yticks([])
    ax.set_xlabel("position (bp)")

    title_lines = [
        f"{rec['seq_id']}  (len {seq_len:,} bp, kind={rec['kind']})",
        f"V1 parent: {rec['v1_parent']}     source: {rec['source_target']}",
    ]
    for ri, r in enumerate(repeats[:3]):
        title_lines.append(
            f"Tandem #{ri+1}: unit={r['repeat_unit_bp']:,} bp × {r['n_copies']} copies "
            f"in [{r['span_start']:,}-{r['span_end']:,}] (span {r['total_span_bp']:,} bp)"
        )
    ax.set_title("\n".join(title_lines), fontsize=9, loc="left")

    plt.tight_layout()
    plt.savefig(out_png, dpi=120, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)
    with open(args.findings) as f:
        findings = json.load(f)
    findings.sort(key=lambda x: -max((r["n_copies"] for r in x["tandem_repeats"]), default=0))
    n = min(len(findings), args.max_figures)
    print(f"Rendering {n} tandem-repeat figures...", file=sys.stderr, flush=True)
    for i, rec in enumerate(findings[:n]):
        if i % 50 == 0:
            print(f"  {i}/{n}", file=sys.stderr, flush=True)
        safe = rec["seq_id"].replace("|", "_").replace(".", "_").replace("/", "_")
        # Prefix with copy count for browsability
        top = rec["tandem_repeats"][0]
        prefix = f"copies{top['n_copies']:03d}_unit{top['repeat_unit_bp']:05d}"
        out_png = os.path.join(args.out, f"{prefix}_{safe}.png")
        try:
            draw_one(rec, out_png)
        except Exception as e:
            print(f"  Error on {safe}: {e}", file=sys.stderr)
    print(f"\nDone. {n} figures saved to {args.out}/", file=sys.stderr)


if __name__ == "__main__":
    main()
