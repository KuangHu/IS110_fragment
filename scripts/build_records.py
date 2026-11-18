#!/usr/bin/env python3
"""Build per-IS110 records using multi-distance anchor data.

Algorithm per ref:
  1. Parse PAF(s): collect anchor pair observations at all D values
  2. For each D, find pairs and classify into v0_filled / empty / etc.
  3. Pick BEST D = smallest D where we found at least 1 empty AND at least 1 V0
     (this gives most direct evidence; smaller D = more genome-specific anchors)
  4. If no D gives both, fall back to D with the most empties (use expected V0)
  5. Compute IS110 size = v0_distance - empty_distance
  5b. BOUNDARY BY DECOMPOSITION: at the smallest D <= --decompose-max-d with
     empty sites, extract the filled interval (source) and each empty interval
     (target) between the anchors and explain filled = empty + insert
     (lib_alleles.decompose_alleles). The insert coordinates ARE the element
     boundaries, per side, with junction microhomology/TSD reported. The
     consensus over empties is used; if nothing decomposes, fall back to the
     peak difference split symmetrically around the transposase and say so
     in boundary_evidence.method (peak_symmetric_fallback). A ref with no
     empty site at all gets method no_empty_tnp_only: its "element" is the
     transposase ORF, not a boundary call.
  6. Extract source sequences using the inferred boundaries
  7. Keep top-N empty observations by mean anchor identity

Usage:
    build_records.py --paf anchors_vs_db.paf \\
                     --ref-table ref_table.tsv \\
                     --hits is_hits.tsv \\
                     --genome-db genome_db.fa \\
                     --out out_dir/
"""

import argparse
import json
import os
import sys
from collections import Counter, defaultdict
from multiprocessing import Pool

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_alleles import (load_anchor_hits, place_per_assembly, source_core,  # noqa: E402
                         gap_tolerance, classify_gap, decompose_alleles,
                         find_tsd, canon_key, revcomp)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--paf", required=True, help="Anchors-vs-DB PAF")
    p.add_argument("--ref-table", required=True,
                   help="ref_table.tsv from extract_anchors.py")
    p.add_argument("--hits", required=True,
                   help="IS hits TSV from is_detect.py "
                        "(is_id, assembly, contig, tnp_start, tnp_end, tnp_strand, ...)")
    p.add_argument("--genome-db", required=True,
                   help="Indexed FASTA (.fai present) containing all source contigs")
    p.add_argument("--source-db", default="",
                   help="Indexed FASTA holding the SOURCE contigs of --hits when "
                        "IS were detected in a different DB than --genome-db "
                        "(default: --genome-db)")
    p.add_argument("--min-is-size", type=int, default=0,
                   help="Smallest removed length accepted as an empty site")
    p.add_argument("--max-is-size", type=int, default=0,
                   help="Largest removed length accepted as an empty site (0 = no "
                        "limit). Same meaning as in call_boundaries.py")
    p.add_argument("--boundaries", default="",
                   help="Optional boundaries TSV from call_boundaries.py")
    p.add_argument("--out", required=True,
                   help="Output directory (writes records.json + sequence FASTAs)")
    p.add_argument("--top-empty", type=int, default=5)
    p.add_argument("--flank-out-len", type=int, default=5000)
    p.add_argument("--min-identity", type=float, default=95)
    p.add_argument("--min-coverage", type=float, default=80)
    p.add_argument("--decompose-max-d", type=int, default=5000,
                   help="Largest anchor distance D used for boundary decomposition "
                        "(default 5000)")
    p.add_argument("--max-align-cells", type=float, default=2e7,
                   help="Skip the SNP-tolerant alignment when len(empty)*len(filled) "
                        "exceeds this (exact decomposition, then the linear-time "
                        "hamming split, are used instead)")
    p.add_argument("--tol-bp", type=int, default=100)
    p.add_argument("--tol-frac", type=float, default=0.002)
    p.add_argument("--margin", type=float, default=0.05)
    p.add_argument("--drop-unsupported", action="store_true",
                   help="write no record for an element with neither a filled nor "
                        "an empty observation in another genome (default: keep it, "
                        "method no_observation_tnp_only)")
    p.add_argument("--max-extra-bp", type=int, default=50000)
    p.add_argument("--threads", type=int, default=32,
                   help="Number of worker processes for sequence extraction (default 32)")
    return p.parse_args()


# ---- multiprocessing worker for parallel record build ----------------------
_FA = None   # per-worker pysam.FastaFile handle: TARGET genomes (empties)
_SRC = None  # per-worker handle for SOURCE contigs; None -> same as _FA


def _worker_init(genome_db, source_db=None):
    """Open the indexed FASTA(s) once per worker process (not fork-safe to share).

    Sources and targets are separate handles: with --source-db every source
    contig used to be looked up in the target DB, failed, and the record was
    dropped without a word (records.json = [] on a clean exit)."""
    global _FA, _SRC
    _FA = pysam.FastaFile(genome_db)
    _SRC = pysam.FastaFile(source_db) if source_db and source_db != genome_db else None


def _fetch(contig, start_0, end_0, source=False):
    """Fetch [start_0, end_0); return uppercase or None. source=True reads the
    source DB."""
    if end_0 <= start_0:
        return ""
    fa = _SRC if (source and _SRC is not None) else _FA
    try:
        seq = fa.fetch(contig, max(0, start_0), end_0)
    except (KeyError, ValueError):
        return None
    return seq.upper() if seq else None


def _empty_interval(eo):
    """Target sequence between the two anchors' inner edges, in transposase
    orientation. gap_start/gap_end are the PROJECTED edges (lib_alleles.
    place_pair), already in genomic order for either strand."""
    if eo["distance"] <= 0:
        return ""
    raw = _fetch(eo["tname"], eo["gap_start"], eo["gap_end"])
    if raw is None:
        return None
    return raw if eo["strand"] == "+" else revcomp(raw)


def _build_one_record(packed):
    """Build a single IS record. Returns tuple or None on failure."""
    ((ref_id, info, best, src), genome_db, flank_out_len, top_empty_n, db_label,
     max_align_cells) = packed

    # Compute IS size from the peaks (used only by the fallback)
    has_empty = bool(best["empty_peak"])
    if best["v0_peak"] and best["empty_peak"]:
        is_length = best["v0_peak"][0] - best["empty_peak"][0]
    elif best["empty_peak"]:
        is_length = best["expected_v0"] - best["empty_peak"][0]
    else:
        is_length = info["tnp_len"]
    # classify_gap accepts an empty when removed >= tnp_len - tol, so the peak
    # difference can come out up to `tol` SHORTER than the transposase; the
    # symmetric split then put both "boundaries" inside the ORF. An element
    # cannot be shorter than the ORF that defines it.
    is_length = max(is_length, info["tnp_len"])

    # Need full contig length to clamp; pysam provides .get_reference_length
    try:
        contig_len = (_SRC or _FA).get_reference_length(src["contig"])
    except (KeyError, ValueError):
        return None

    tnp_start_0 = src["start"] - 1
    tnp_end_0 = src["end"]

    dec = None
    if best.get("decompose"):
        dec = decompose_boundary(src, info, best["decompose"], top_empty_n,
                                 max_align_cells)
    # An "empty" is first called on length alone (removed >= tnp_len). When
    # decomposition shows that what the target lacks does NOT contain the
    # transposase, it is a flank deletion next to an intact element, not an
    # empty site: drop it from the evidence (it used to yield a 1,500-bp
    # "element" from a 1,500-bp upstream deletion).
    n_rejected = 0
    if dec and dec.get("rejected_assemblies"):
        rej = set(dec["rejected_assemblies"])
        # Only the top-N empties are decomposed. The others were left in the
        # evidence untested: six genomes with the same 1,500-bp flank
        # deletion gave 5 rejected + 1 "empty" and a 1,500-bp element. An
        # unexamined empty with the same gap as a rejected one (and unlike
        # every accepted one) is the same deletion.
        acc = set(dec.get("accepted_assemblies", ()))
        rej_gaps = [e["distance"] for e in best["empty_obs"] if e.get("assembly") in rej]
        acc_gaps = [e["distance"] for e in best["empty_obs"] if e.get("assembly") in acc]
        near = lambda g, gs: any(abs(g - x) <= 2 for x in gs)   # noqa: E731
        kept = [e for e in best["empty_obs"]
                if e.get("assembly") not in rej
                and (e.get("assembly") in acc or near(e["distance"], acc_gaps)
                     or not near(e["distance"], rej_gaps))]
        n_rejected = len(best["empty_obs"]) - len(kept)
        best = dict(best, empty_obs=kept,
                    empty_peak=length_mode([e["distance"] for e in kept]))
        has_empty = bool(best["empty_peak"])
        if best["v0_peak"] and best["empty_peak"]:
            is_length = best["v0_peak"][0] - best["empty_peak"][0]
        elif best["empty_peak"]:
            is_length = best["expected_v0"] - best["empty_peak"][0]
        else:
            is_length = info["tnp_len"]
        is_length = max(is_length, info["tnp_len"])

    if dec and dec["ok"]:
        method = "empty_vs_filled_decomposition"
        is_start_0, is_end_0 = dec["is_start_0"], dec["is_end_0"]
        boundary_5p_offset, boundary_3p_offset = dec["off5"], dec["off3"]
    else:
        # FALLBACK: the peak difference only gives a length, so the extension
        # is split evenly between the two sides. It is not a boundary call.
        # With no empty site at all there is no length evidence either: the
        # "element" is just the transposase ORF. Label it so it is never
        # mistaken for a boundary call (it used to share the fallback label).
        # no_empty_tnp_only    : seen filled in >= 1 OTHER genome, never empty
        # no_observation_tnp_only: neither a filled nor an empty observation in
        #   any other genome (no anchor pair at all, or only pairs that are
        #   neither -- larger than the filled state, or overlapping anchors).
        #   Such records used to be labelled no_empty_tnp_only, as if a filled
        #   observation stood behind them (or, with no pair at all, were
        #   dropped without a record).
        if has_empty:
            method = "peak_symmetric_fallback"
        elif best["v0_obs"]:
            method = "no_empty_tnp_only"
        else:
            method = "no_observation_tnp_only"
        ext5 = (is_length - info["tnp_len"]) // 2
        ext3 = is_length - info["tnp_len"] - ext5
        # offsets are in TRANSPOSASE orientation: on the minus strand the 5'
        # side is the genomic right. They used to be applied left/right
        # regardless of strand, so coordinates and reported offsets disagreed
        # by 1 bp for odd extensions.
        left, right = (ext5, ext3) if src["strand"] == "+" else (ext3, ext5)
        is_start_0 = max(0, tnp_start_0 - left)
        is_end_0 = min(contig_len, tnp_end_0 + right)
        left, right = tnp_start_0 - is_start_0, is_end_0 - tnp_end_0   # after clamping
        ext5, ext3 = (left, right) if src["strand"] == "+" else (right, left)
        boundary_5p_offset, boundary_3p_offset = -ext5, ext3

    up_flank_start = max(0, is_start_0 - flank_out_len)
    up_flank_end = is_start_0
    down_flank_start = is_end_0
    down_flank_end = min(contig_len, is_end_0 + flank_out_len)

    tnp_seq = _fetch(src["contig"], tnp_start_0, tnp_end_0, source=True)
    is_seq = _fetch(src["contig"], is_start_0, is_end_0, source=True)
    up_seq = _fetch(src["contig"], up_flank_start, up_flank_end, source=True)
    down_seq = _fetch(src["contig"], down_flank_start, down_flank_end, source=True)
    if tnp_seq is None or is_seq is None or up_seq is None or down_seq is None:
        return None

    if src["strand"] == "-":
        tnp_seq = revcomp(tnp_seq)
        is_seq = revcomp(is_seq)
        up_seq, down_seq = revcomp(down_seq), revcomp(up_seq)

    empty_obs_sorted = sorted(best["empty_obs"],
                              key=lambda x: (-x["mean_ident"], x["assembly"]))
    top_empty = empty_obs_sorted[:top_empty_n]

    # Empty junction sequences
    empties_out = []
    for i, eo in enumerate(top_empty):
        target = eo.get("tname") or eo.get("contig")
        if not target:
            continue
        junction = _empty_interval(eo)
        if junction:
            eo["junction_sequence"] = junction
            eo["junction_length"] = len(junction)
            tag = eo.get("assembly") or target or "x"
            empties_out.append((f"{ref_id}__empty{i}_{tag}", junction))

    record = {
        "ref_id": ref_id,
        # orientation-invariant key of the element + 100 bp either side:
        # the same insertion seen in several genomes shares it
        "event_key": canon_key((up_seq[-100:] if up_seq else "") + is_seq
                               + (down_seq[:100] if down_seq else "")),
        "source": {
            "assembly": src["assembly"], "contig": src["contig"],
            "transposase_start": src["start"], "transposase_end": src["end"],
            "transposase_strand": src["strand"],
        },
        "is_element": {
            "length": is_end_0 - is_start_0,
            "start_offset_5p": boundary_5p_offset,
            "end_offset_3p": boundary_3p_offset,
            "source_start": is_start_0 + 1, "source_end": is_end_0,
            "sequence": is_seq,
        },
        "transposase_cds": {"length": len(tnp_seq), "sequence": tnp_seq},
        "upstream_flank": {"length": len(up_seq), "sequence": up_seq},
        "downstream_flank": {"length": len(down_seq), "sequence": down_seq},
        "boundary_evidence": {
            "method": method,
            "decomposition": dec,
            "source_db": db_label,
            "anchor_distance_D": best["D"],
            "v0_peak_distance": best["v0_peak"][0] if best["v0_peak"] else None,
            "v0_peak_count": best["v0_peak"][1] if best["v0_peak"] else 0,
            "empty_peak_distance": best["empty_peak"][0] if best["empty_peak"] else None,
            "empty_peak_count": best["empty_peak"][1] if best["empty_peak"] else 0,
            # how many 50-bp bins share the top count (1 = unique mode)
            "v0_peak_n_tied": mode_ties([p["distance"] for p in best["v0_obs"]]),
            "empty_peak_n_tied": mode_ties([p["distance"] for p in best["empty_obs"]]),
            # support of the call that decided the element, and how many
            # candidates had that same support (1 = no tie):
            #   decomposition : empties voting for the chosen (off5, off3) /
            #                   distinct calls with as many votes
            #   fallback      : observations in the chosen empty-length bin /
            #                   bins with as many
            #   no empty site : 0 / 0 (nothing was chosen)
            "call_support": (dec["n_consensus"] if dec and dec["ok"] else
                             (best["empty_peak"][1] if has_empty else 0)),
            "call_n_tied": (dec["n_tied_calls"] if dec and dec["ok"] else
                            (mode_ties([p["distance"] for p in best["empty_obs"]])
                             if has_empty else 0)),
            "is_element_length_inferred": is_length,
            "call_confidence": None,       # filled in below from the two fields above
            "n_v0_observations": len(best["v0_obs"]),
            "n_empty_observations": len(best["empty_obs"]),
            "n_empty_rejected_flank_deletion": n_rejected,
            "n_v1plus_observations": len(best["v1plus_obs"]),
            "n_ambiguous_placements": best.get("n_ambiguous", 0),
            "n_empty_anchor_overlap": best.get("n_anchor_overlap_all_d", 0),
            "counts_are": "distinct assemblies (GCA/GCF twins collapsed)",
            "n_top_empty_kept": len(top_empty),
        },
        "filled_observations": [
            {k: v for k, v in p.items() if k != "tname"} for p in best["v0_obs"][:50]
        ],
        "empty_observations": [
            {k: v for k, v in p.items() if k != "tname"} for p in top_empty
        ],
        "v1plus_observations": [
            {k: v for k, v in p.items() if k != "tname"} for p in best["v1plus_obs"][:50]
        ],
    }
    _be = record["boundary_evidence"]
    _be["call_confidence"] = call_confidence(method, _be["call_support"], _be["call_n_tied"])
    return (record, is_seq, tnp_seq, up_seq, down_seq, empties_out)


def length_mode(values, bin_size=50):
    """Median of the most populated bin -> (length, n_in_bin), or None."""
    if not values:
        return None
    bins = defaultdict(list)
    for v in values:
        bins[v // bin_size].append(v)
    best = max(bins.values(), key=len)
    best.sort()
    return (best[len(best) // 2], len(best))


def call_confidence(method, support, n_tied):
    """One flag from call_support / call_n_tied (reporting only):
         tied      another candidate had as much support as the chosen one
         support1  the chosen call rests on a single observation
         ok        support >= 2 and no tie
         no_call   no boundary or length was called (transposase-only records)
    'tied' wins when both apply."""
    if method not in ("empty_vs_filled_decomposition", "peak_symmetric_fallback"):
        return "no_call"
    if n_tied > 1:
        return "tied"
    return "support1" if support <= 1 else "ok"


def mode_ties(values, bin_size=50):
    """Number of bins that share the largest size (1 = the mode is unique).
    length_mode() takes the FIRST such bin in the order of `values`, i.e. the
    bin of the first observation in sorted assembly order -- deterministic,
    but not a judgement. Reporting only: this function chooses nothing."""
    if not values:
        return 0
    bins = Counter(v // bin_size for v in values)
    top = max(bins.values())
    return sum(1 for n in bins.values() if n == top)


def pair_anchors(hits, ref_id, D, info, opts):
    """One scored placement per distinct assembly (lib_alleles.place_pair).

    Replaces an every-up-x-every-down loop that counted a genome once per
    anchor-copy combination and could pair anchors on different repeat
    copies. Ambiguous placements are counted but never used as evidence."""
    expected_v0 = 2 * D + info["tnp_len"]
    placements = place_per_assembly(
        hits.get((ref_id, "up", D), {}), hits.get((ref_id, "down", D), {}),
        exclude_core=source_core(ref_id),
        min_gap=-50, max_gap=expected_v0 + opts["max_extra_bp"],
        margin_ratio=opts["margin"])
    pairs, n_ambiguous = [], 0
    for pl in placements.values():
        if pl["status"] == "ambiguous":
            n_ambiguous += 1
            continue
        u, d = pl["up"], pl["down"]
        pairs.append({
            "assembly": pl["assembly"], "contig": pl["tname"].split("|", 1)[-1],
            "tname": pl["tname"], "strand": pl["strand"], "distance": pl["gap"],
            "up_ts": u.ts, "up_te": u.te, "down_ts": d.ts, "down_te": d.te,
            "gap_start": pl["start"], "gap_end": pl["end"],
            "up_ident": u.ident, "down_ident": d.ident,
            "mean_ident": (u.ident + d.ident) / 2, "D": D,
            "placement_status": pl["status"], "placement_margin": round(pl["margin"], 3),
        })
    return pairs, expected_v0, n_ambiguous


def categorize_pairs(pairs, expected_v0, tnp_len, opts):
    """Split placements into V0 / empty / V1+ against the source state.

    The old version called an 'empty peak' any histogram bin > 500 bp short
    of the source and then took everything within 250 bp of it; with a
    tolerance tied to D it also let filled sites read as empty at large D.
    Classification is now lib_alleles.classify_gap with an absolute
    tolerance; peaks are summaries of the classified sets."""
    tol = gap_tolerance(expected_v0, opts["tol_bp"], opts["tol_frac"])
    v0_obs, empty_obs, v1plus_obs = [], [], []
    n_overlap = 0
    for p in pairs:
        cat = classify_gap(p["distance"], expected_v0, tnp_len, tol,
                           max_removed=opts.get("max_is_size") or None,
                           min_removed=opts.get("min_is_size") or None)
        if cat == "v0_filled":
            v0_obs.append(p)
        elif cat == "empty":
            # A NEGATIVE projected gap means the two anchors overlap on the
            # target: the element reaches into the anchors at this D. The
            # interval between the anchors then holds only part of the
            # element, and decomposing it against an empty junction of ""
            # returned that part as an exact call (3,200 bp for a 3,240-bp
            # element, consensus_support 1.0). Not usable at this D -- a
            # larger D, whose anchors lie outside the element, is.
            if p["distance"] < 0:
                n_overlap += 1
            else:
                empty_obs.append(p)
        elif cat == "v1plus_filled":
            v1plus_obs.append(p)
    v0_peak = length_mode([p["distance"] for p in v0_obs])
    empty_peak = length_mode([p["distance"] for p in empty_obs])
    return v0_peak, empty_peak, v0_obs, empty_obs, v1plus_obs, n_overlap


def pick_best_d(hits, ref_id, info, distances, opts):
    """Find the best D for this ref: prefer smallest D with both V0 and empty.

    Also returns, as best['decompose'], the smallest D <= decompose_max_d that
    has empty sites -- the interval used for base-level decomposition."""
    candidates = []
    for D in distances:
        pairs, expected_v0, n_amb = pair_anchors(hits, ref_id, D, info, opts)
        if not pairs:
            continue
        v0_peak, empty_peak, v0_obs, empty_obs, v1plus_obs, n_overlap = categorize_pairs(
            pairs, expected_v0, info["tnp_len"], opts)
        candidates.append({
            "D": D, "pairs": pairs, "expected_v0": expected_v0,
            "n_ambiguous": n_amb, "n_anchor_overlap": n_overlap,
            "v0_peak": v0_peak, "empty_peak": empty_peak,
            "v0_obs": v0_obs, "empty_obs": empty_obs, "v1plus_obs": v1plus_obs,
        })

    if not candidates:
        # No anchor pair in any other genome. The transposase hit is real, so
        # the element is still catalogued -- as an explicitly unsupported
        # record (method no_observation_tnp_only), unless --drop-unsupported.
        D0 = distances[0]
        return {"D": D0, "pairs": [], "expected_v0": 2 * D0 + info["tnp_len"],
                "n_ambiguous": 0, "n_anchor_overlap": 0, "n_anchor_overlap_all_d": 0,
                "v0_peak": None, "empty_peak": None, "v0_obs": [], "empty_obs": [],
                "v1plus_obs": [], "decompose": None}

    best = None
    # Prefer smallest D with both V0 and empty (highest-quality boundary call)
    for c in sorted(candidates, key=lambda x: x["D"]):
        if c["v0_peak"] and c["empty_peak"]:
            best = c
            break
    if best is None:
        # Else: any D with empty
        for c in sorted(candidates, key=lambda x: x["D"]):
            if c["empty_peak"]:
                best = c
                break
    if best is None:
        # Else: just pick the first one (only V0)
        best = candidates[0]

    best = dict(best)
    # empties seen only with overlapping anchors, at any D (see categorize_pairs)
    best["n_anchor_overlap_all_d"] = sum(c["n_anchor_overlap"] for c in candidates)
    best["decompose"] = None
    for c in sorted(candidates, key=lambda x: x["D"]):
        if c["D"] <= opts["decompose_max_d"] and c["empty_obs"]:
            best["decompose"] = {"D": c["D"], "empty_obs": c["empty_obs"]}
            break
    return best


ORF_SLACK = 10


def decompose_boundary(src, info, dec, top_n, max_align_cells):
    """Base-level element boundaries from empty-vs-filled decomposition.

    filled = source interval between the D-anchors, oriented with the
    transposase on +; empty = each target's interval between the same anchors.
    decompose_alleles() explains filled as empty + insert; the insert's
    coordinates in the filled interval are the element boundaries. Each empty
    votes (off5, off3) relative to the transposase; the most common call wins.
    Returns a dict (or None when no empty decomposes to an insert that
    contains the whole transposase)."""
    D, tnp_len = dec["D"], info["tnp_len"]
    g0 = src["start"] - 1 - D                 # genomic 0-based start of filled
    g1 = src["end"] + D
    filled = _fetch(src["contig"], g0, g1, source=True)
    if not filled or len(filled) != g1 - g0:
        return None
    if src["strand"] == "-":
        filled = revcomp(filled)
    cells_ok = lambda e: len(e) * len(filled) <= max_align_cells  # noqa: E731

    votes = Counter()
    calls = []
    methods = Counter()
    n_try = 0
    rejected = []
    # ties broken by assembly name: the top-N cut must not depend on dict/set
    # order (identical inputs gave 1,400 bp or 1,600 bp by PYTHONHASHSEED)
    for eo in sorted(dec["empty_obs"],
                     key=lambda x: (-x["mean_ident"], x["assembly"]))[:top_n]:
        n_try += 1
        empty = _empty_interval(eo)
        if empty is None:
            continue
        r = decompose_alleles(empty, filled, min_insert=max(50, tnp_len - 20),
                              tolerant=cells_ok(empty))
        methods[r["method"]] += 1
        if r["method"] == "none":
            continue
        st, ln = r["insert_start"], r["insert_len"]
        # Canonical placement before voting. With junction microhomology the
        # same insert can sit anywhere in the repeat; the exact method returns
        # the rightmost placement, the alignment-based ones whatever the
        # aligner chose. Unnormalised, identical calls split their votes
        # (lower consensus_support, more n_distinct_calls) and the reported
        # boundary depended on which empty happened to be first. Slide every
        # call to the rightmost placement (= the exact method's convention).
        while st + ln < len(filled) and filled[st] == filled[st + ln]:
            st += 1
        lo = st
        while lo > 0 and filled[lo - 1] == filled[lo - 1 + ln]:
            lo -= 1
        # The element must contain the transposase ORF. Among the equivalent
        # placements [lo, st] take the rightmost one that still starts at or
        # before the ORF, so sliding never pushes the boundary into the ORF.
        st = max(lo, min(st, D))
        # ORF_SLACK: the ORF caller may extend the start codon a few bp past
        # the true element end, so a call that clips <= ORF_SLACK bp of the
        # predicted ORF is kept -- and REPORTED (orf_clipped_5p/3p), never
        # silently presented as containing the whole ORF.
        if st > D + ORF_SLACK or st + ln < D + tnp_len - ORF_SLACK:
            methods["insert_misses_tnp"] += 1
            rejected.append(eo.get("assembly", ""))
            continue
        off5, off3 = st - D, (st + ln) - (D + tnp_len)
        votes[(off5, off3)] += 1
        calls.append((off5, off3, st, ln, r, eo.get("assembly", "")))
    if not votes:
        return {"ok": False, "D": D, "n_attempted": n_try,
                "method_counts": dict(methods), "rejected_assemblies": rejected,
                "accepted_assemblies": []}
    # most votes; ties -> the call that keeps the whole ORF, then the shorter
    # element, then the offsets themselves (fully deterministic)
    (off5, off3), n_cons = min(
        votes.items(),
        key=lambda kv: (-kv[1], max(0, kv[0][0]) + max(0, -kv[0][1]),
                        kv[0][1] - kv[0][0], kv[0]))
    n_tied = sum(1 for v in votes.values() if v == n_cons)
    off5, off3, st, ln, r, asm = next(c for c in calls if (c[0], c[1]) == (off5, off3))
    ins = filled[st:st + ln]
    tsd_len, tsd_seq, tsd_side, tsd_conf = find_tsd(ins, filled[:st], filled[st + ln:])
    # genomic coordinates of the element (0-based half-open)
    if src["strand"] == "+":
        is0, is1 = g0 + st, g0 + st + ln
    else:
        is0, is1 = g1 - (st + ln), g1 - st
    return {"ok": True, "D": D, "is_start_0": is0, "is_end_0": is1,
            "off5": off5, "off3": off3, "length": ln,
            "n_attempted": n_try, "n_decomposed": len(calls),
            "n_consensus": n_cons,
            "consensus_support": round(n_cons / len(calls), 3),
            "n_distinct_calls": len(votes),
            "consensus_tied": n_tied > 1,
            "n_tied_calls": n_tied,
            "orf_clipped_5p": max(0, off5), "orf_clipped_3p": max(0, -off3),
            "rejected_assemblies": rejected,
            "accepted_assemblies": [c[5] for c in calls],
            "method_counts": dict(methods),
            "consensus_method": r["method"], "event_class": r["event_class"],
            "empty_offset": r["offset"] + (st - r["insert_start"]),
            "placement": "rightmost_within_microhomology",
            "junction_microhomology_bp": r["junction_microhomology_bp"],
            "target_lost_bp": r["target_lost_bp"],
            "tsd": {"length": tsd_len, "sequence": tsd_seq, "side": tsd_side,
                    "confidence": tsd_conf},
            "example_empty_assembly": asm}


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    # Reference table
    ref_info = {}
    distances_set = set()
    with open(args.ref_table) as f:
        next(f)
        for line in f:
            parts = line.rstrip().split("\t")
            ref_id, ref_len, anchor_len, tnp_len, flank_size, distances_str = parts
            ds = [int(d) for d in distances_str.split(",")]
            distances_set.update(ds)
            ref_info[ref_id] = {
                "ref_len": int(ref_len), "anchor_len": int(anchor_len),
                "tnp_len": int(tnp_len), "flank_size": int(flank_size),
                "distances": ds,
            }
    distances = sorted(distances_set)

    # Source positions: load is_hits.tsv from is_detect.py
    source_positions = {}
    with open(args.hits) as f:
        header = f.readline().rstrip().split("\t")
        col = {n: i for i, n in enumerate(header)}
        for line in f:
            parts = line.rstrip().split("\t")
            if len(parts) < len(header): continue
            is_id = parts[col["is_id"]]
            source_positions[is_id] = {
                "assembly": parts[col["assembly"]] if "assembly" in col else parts[col["contig"]],
                "contig": parts[col["contig"]],
                "start": int(parts[col["tnp_start"]]),
                "end": int(parts[col["tnp_end"]]),
                "strand": parts[col["tnp_strand"]],
            }

    # Parse anchor PAF
    print(f"Parsing anchor PAF...", file=sys.stderr)
    ncbi_hits, n_total, n_kept = load_anchor_hits(args.paf, args.min_identity,
                                                  args.min_coverage)
    print(f"  PAF: {n_total:,} rows, {n_kept:,} kept", file=sys.stderr)

    # Build records — parallelize per-ref using pysam (in-process FASTA reads)
    n_done = n_neither = 0
    db_label = os.path.basename(args.genome_db)

    # Pre-filter refs to those with a valid pick (cheap; no FASTA I/O)
    work_items = []
    opts = {"tol_bp": args.tol_bp, "tol_frac": args.tol_frac,
            "margin": args.margin, "max_extra_bp": args.max_extra_bp,
            "decompose_max_d": args.decompose_max_d,
            "min_is_size": args.min_is_size, "max_is_size": args.max_is_size}
    for ref_id, info in ref_info.items():
        best = pick_best_d(ncbi_hits, ref_id, info, distances, opts)
        if not best["v0_obs"] and not best["empty_obs"]:
            n_neither += 1
            if args.drop_unsupported:
                continue
        src = source_positions.get(ref_id)
        if not src:
            continue
        n_done += 1
        work_items.append((ref_id, info, best, src))

    print(f"  Refs to process: {len(work_items):,} ({n_neither:,} without any filled or "
          f"empty observation in another genome: "
          f"{'dropped' if args.drop_unsupported else 'kept as no_observation_tnp_only'})",
          file=sys.stderr, flush=True)

    # Worker pool: each worker opens its own pysam.FastaFile (not fork-safe to share)
    all_records = []
    is_elem_seqs = []
    tnp_seqs = []
    up_flank_seqs = []
    down_flank_seqs = []
    empty_junction_seqs = []

    pool_args = [(item, args.genome_db, args.flank_out_len, args.top_empty, db_label,
                  args.max_align_cells)
                 for item in work_items]
    n_failed = 0
    with Pool(args.threads, initializer=_worker_init,
              initargs=(args.genome_db, args.source_db)) as pool:
        for i, packed in enumerate(pool.imap_unordered(_build_one_record, pool_args,
                                                       chunksize=8)):
            if i % 1000 == 0:
                print(f"    built {i:,}/{len(pool_args):,}", file=sys.stderr, flush=True)
            if packed is None:
                n_failed += 1
                continue
            record, is_seq, tnp_seq, up_seq, down_seq, empties = packed
            all_records.append(record)
            ref_id = record["ref_id"]
            is_elem_seqs.append((ref_id, is_seq))
            tnp_seqs.append((ref_id, tnp_seq))
            up_flank_seqs.append((ref_id, up_seq))
            down_flank_seqs.append((ref_id, down_seq))
            for tag, seq in empties:
                empty_junction_seqs.append((tag, seq))

    # Write outputs
    # imap_unordered returns in completion order; sort so the files are
    # reproducible run to run
    order = sorted(range(len(all_records)), key=lambda i: all_records[i]["ref_id"])
    all_records = [all_records[i] for i in order]
    is_elem_seqs = [is_elem_seqs[i] for i in order]
    tnp_seqs = [tnp_seqs[i] for i in order]
    up_flank_seqs = [up_flank_seqs[i] for i in order]
    down_flank_seqs = [down_flank_seqs[i] for i in order]
    empty_junction_seqs.sort()
    if n_failed:
        print(f"WARNING: {n_failed:,} of {len(pool_args):,} refs produced no record "
              f"(source contig or sequence not found in the FASTA -- wrong "
              f"--source-db / --genome-db?)", file=sys.stderr, flush=True)
    if pool_args and not all_records:
        sys.exit("ERROR: no record could be built; refusing to write an empty "
                 "records.json")
    json_path = os.path.join(args.out, "records.json")
    with open(json_path + ".tmp", "w") as f:
        json.dump(all_records, f, indent=2)
    os.replace(json_path + ".tmp", json_path)
    print(f"Wrote {len(all_records)} records to {json_path}", file=sys.stderr)

    def write_fa(path, seqs):
        with open(path, "w") as f:
            for name, seq in seqs:
                f.write(f">{name}\n")
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + "\n")

    write_fa(os.path.join(args.out, "is_elements.fa"), is_elem_seqs)
    write_fa(os.path.join(args.out, "transposase_cds.fa"), tnp_seqs)
    write_fa(os.path.join(args.out, "upstream_flanks.fa"), up_flank_seqs)
    write_fa(os.path.join(args.out, "downstream_flanks.fa"), down_flank_seqs)
    write_fa(os.path.join(args.out, "empty_junctions.fa"), empty_junction_seqs)

    # Summary TSV
    with open(os.path.join(args.out, "summary.tsv"), "w") as f:
        f.write("ref_id\tassembly\tcontig\ttnp_start\ttnp_end\ttnp_strand\t"
                "is_element_len\ttnp_len\tD_used\tsource_db\tv0_peak\tempty_peak\t"
                "n_v0\tn_empty\tn_v1plus\tn_top_empty\t"
                "boundary_method\toffset_5p\toffset_3p\tdecomp_D\tdecomp_support\t"
                "junction_microhomology_bp\ttsd_len\tevent_key\t"
                "call_support\tcall_n_tied\tv0_peak_count\tv0_peak_n_tied\t"
                "empty_peak_count\tempty_peak_n_tied\tcall_confidence\n")
        for r in all_records:
            be = r["boundary_evidence"]
            dc = be.get("decomposition") or {}
            f.write(f"{r['ref_id']}\t{r['source']['assembly']}\t{r['source']['contig']}\t"
                    f"{r['source']['transposase_start']}\t{r['source']['transposase_end']}\t"
                    f"{r['source']['transposase_strand']}\t"
                    f"{r['is_element']['length']}\t{r['transposase_cds']['length']}\t"
                    f"{be['anchor_distance_D']}\t{be['source_db']}\t"
                    f"{be['v0_peak_distance']}\t{be['empty_peak_distance']}\t"
                    f"{be['n_v0_observations']}\t{be['n_empty_observations']}\t"
                    f"{be['n_v1plus_observations']}\t{be['n_top_empty_kept']}\t"
                    f"{be['method']}\t{r['is_element']['start_offset_5p']}\t"
                    f"{r['is_element']['end_offset_3p']}\t"
                    f"{dc.get('D', '')}\t{dc.get('consensus_support', '')}\t"
                    f"{dc.get('junction_microhomology_bp', '')}\t"
                    f"{(dc.get('tsd') or {}).get('length', '')}\t{r['event_key']}\t"
                    f"{be['call_support']}\t{be['call_n_tied']}\t"
                    f"{be['v0_peak_count']}\t{be['v0_peak_n_tied']}\t"
                    f"{be['empty_peak_count']}\t{be['empty_peak_n_tied']}\t"
                    f"{be['call_confidence']}\n")

    print(f"\n=== Build summary ===", file=sys.stderr)
    print(f"  Records built:       {len(all_records):,} of {n_done:,} refs with "
          f"anchor placements", file=sys.stderr)
    n_dec = sum(1 for r in all_records
                if r["boundary_evidence"]["method"] == "empty_vs_filled_decomposition")
    n_fb = sum(1 for r in all_records
               if r["boundary_evidence"]["method"] == "peak_symmetric_fallback")
    print(f"  Boundary by decomposition: {n_dec:,} / {len(all_records):,}; "
          f"symmetric fallback (length only, NOT a boundary): {n_fb:,}; "
          f"no empty site (transposase ORF only): "
          f"{len(all_records) - n_dec - n_fb:,}", file=sys.stderr)
    # event_key hashes the exact element + 100 bp of context: one SNP in the
    # flank gives a new key, so this is an UPPER bound on insertion events.
    print(f"  Distinct exact element+context keys (event_key; upper bound on "
          f"insertion events): {len({r['event_key'] for r in all_records}):,}",
          file=sys.stderr)
    print(f"  No observation in another genome: {n_neither}", file=sys.stderr)
    print(f"  All outputs in: {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
