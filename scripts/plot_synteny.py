"""Plot synteny ribbon diagrams between an IS element and a target genome.

Reads alignment hits from a PAF file, optionally annotates the query track
with CDS features from a GenBank file, and auto-flips the target region
if the dominant alignment strand is reverse.

Usage:
    python plot_synteny.py hits.paf --target NZ_CP084874.1 \
        --gbk data/query.gbk --output-prefix results/synteny
"""

import argparse
import os
import sys

# Allow importing parse_paf from same directory
sys.path.insert(0, os.path.dirname(__file__))
from parse_paf import filter_hits, parse_paf


def is_transposase(note: str) -> bool:
    """Check if CDS is IS110 transposase by domain annotation."""
    return "DEDD" in note or "Tnp20" in note


def cds_color(strand: int, note: str) -> str:
    """Color CDS by type."""
    if is_transposase(note):
        return "#EE6677"  # red - transposase
    if strand == 1:
        return "#4477AA"  # blue - forward
    return "#66CCEE"      # cyan - reverse


def parse_gbk(gbk_path: str) -> tuple[int, list[tuple]]:
    """Parse CDS features from a GenBank file.

    Returns (sequence_length, list of (start, end, strand, label, note)).
    """
    from Bio import SeqIO

    rec = SeqIO.read(gbk_path, "genbank")
    features = []
    for f in rec.features:
        if f.type == "CDS":
            start = int(f.location.start)
            end = int(f.location.end)
            strand = f.location.strand
            label = f.qualifiers.get("label", [""])[0]
            note = f.qualifiers.get("note", [""])[0]
            features.append((start, end, strand, label, note))
    return len(rec), features


def plot(hits, query_len, target_name, cds_features=None,
         pad=2000, fig_width=20, auto_flip=True, output_prefix="synteny"):
    """Build and save the synteny ribbon plot."""
    from pygenomeviz import GenomeViz

    if not hits:
        print("No hits to plot.", file=sys.stderr)
        return

    # ── Auto-flip detection ──
    fwd_bases = sum(h.query_end - h.query_start for h in hits if h.strand == "+")
    rev_bases = sum(h.query_end - h.query_start for h in hits if h.strand == "-")
    flip_target = auto_flip and (rev_bases > fwd_bases)
    if flip_target:
        print(f"Auto-flip: {rev_bases} bp (-) vs {fwd_bases} bp (+) → flipping target")

    # ── Target region ──
    tmin = min(h.target_start for h in hits)
    tmax = max(h.target_end for h in hits)
    tstart = tmin - pad
    tend = tmax + pad
    tlen = tend - tstart

    if flip_target:
        track2_name = f"{target_name} ({tend:,}→{tstart:,}) [rev comp]"
    else:
        track2_name = f"{target_name} ({tstart:,}→{tend:,})"

    # ── Query track name ──
    query_track_name = hits[0].query_name

    # ── Build plot ──
    gv = GenomeViz(
        fig_width=fig_width,
        fig_track_height=1.0,
        track_align_type="left",
        feature_track_ratio=0.5 if cds_features else 0.25,
        link_track_ratio=1.0,
    )

    # Track 1: Query
    track1 = gv.add_feature_track(query_track_name, query_len)

    if cds_features:
        for start, end, strand, label, note in cds_features:
            show_label = ""
            if is_transposase(note):
                show_label = "IS110 (DEDD/Tnp20)"
            elif (end - start) > 1500:
                show_label = label

            track1.add_feature(
                start, end, strand=strand,
                plotstyle="arrow",
                fc=cds_color(strand, note),
                ec="black", lw=0.3,
                label=show_label,
                text_kws=dict(size=7, color="black"),
            )

    # Track 2: Target genome (zoomed region)
    track2 = gv.add_feature_track(track2_name, tlen)

    for h in hits:
        t1 = h.target_start - tstart
        t2 = h.target_end - tstart
        strand = h.strand
        if flip_target:
            t1_flip = tlen - t1
            t2_flip = tlen - t2
            t1, t2 = t2_flip, t1_flip
            strand = "+" if strand == "-" else "-"
        st = -1 if strand == "-" else 1
        track2.add_feature(min(t1, t2), max(t1, t2), strand=st,
                           plotstyle="arrow", fc="#CCCCCC", ec="black", lw=0.3)

    # ── Add links ──
    min_ident = min(h.identity for h in hits)
    min_ident_pct = int(min_ident * 100)

    for h in hits:
        identity_pct = round(h.identity * 100, 1)
        t1 = h.target_start - tstart
        t2 = h.target_end - tstart
        strand = h.strand

        if flip_target:
            t1_flip = tlen - t1
            t2_flip = tlen - t2
            t1, t2 = t2_flip, t1_flip
            strand = "+" if strand == "-" else "-"

        link1 = (query_track_name, h.query_start, h.query_end)
        if strand == "-":
            link2 = (track2_name, t2, t1)
        else:
            link2 = (track2_name, t1, t2)

        gv.add_link(
            link1, link2,
            color="blue",
            inverted_color="red",
            v=identity_pct,
            vmin=min_ident_pct,
            vmax=100,
            curve=True,
            alpha=0.6,
        )

    # ── Colorbar ──
    colors = ["blue"]
    has_inv = any(
        (h.strand == "-") != flip_target for h in hits
    )
    if has_inv:
        colors.append("red")

    gv.set_colorbar(
        colors=colors,
        vmin=min_ident_pct,
        vmax=100,
        bar_label="Identity (%)",
        bar_height=0.15,
        bar_width=0.01,
        bar_labelsize=10,
        tick_labelsize=8,
    )

    # ── Save ──
    gv.savefig(f"{output_prefix}.png", dpi=200)
    gv.savefig(f"{output_prefix}.pdf")
    print(f"Saved: {output_prefix}.png / .pdf")


def main():
    parser = argparse.ArgumentParser(
        description="Plot synteny ribbon between IS element and target genome")
    parser.add_argument("paf_file", help="PAF alignment file")
    parser.add_argument("gbk_file", help="GenBank file for query sequence and CDS annotations")
    parser.add_argument("--target", required=True,
                        help="Target genome name to plot against")
    parser.add_argument("--min-identity", type=float, default=0.0,
                        help="Minimum identity threshold (0-1, default: 0.0)")
    parser.add_argument("--pad", type=int, default=2000,
                        help="Padding (bp) around target region (default: 2000)")
    parser.add_argument("--output-prefix", default="synteny",
                        help="Output prefix for PNG/PDF (default: synteny)")
    parser.add_argument("--fig-width", type=float, default=20,
                        help="Figure width in inches (default: 20)")
    parser.add_argument("--no-flip", action="store_true",
                        help="Disable auto-flip of target track")

    args = parser.parse_args()

    # Load and filter hits
    hits = parse_paf(args.paf_file)
    hits = filter_hits(hits, min_identity=args.min_identity, target=args.target)

    if not hits:
        print(f"No hits found for target '{args.target}'", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(hits)} hits against {args.target}")

    # Parse GBK for query length and CDS annotations
    query_len, cds_features = parse_gbk(args.gbk_file)
    print(f"Loaded {len(cds_features)} CDS features from {args.gbk_file}")

    plot(
        hits=hits,
        query_len=query_len,
        target_name=args.target,
        cds_features=cds_features,
        pad=args.pad,
        fig_width=args.fig_width,
        auto_flip=not args.no_flip,
        output_prefix=args.output_prefix,
    )


if __name__ == "__main__":
    main()
