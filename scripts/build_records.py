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
     in boundary_evidence.method.
  6. Extract source sequences using the inferred boundaries
  7. Keep top-N empty observations by mean anchor identity
  8. If NCBI gave no empties, also include Logan empties in the record

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
import re
import sys
from collections import Counter, defaultdict
from multiprocessing import Pool

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_alleles import (load_anchor_hits, place_per_assembly,  # noqa: E402
                         gap_tolerance, classify_gap, decompose_alleles,
                         find_tsd, canon_key)


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
                        "exceeds this (exact decomposition is still tried)")
    p.add_argument("--tol-bp", type=int, default=100)
    p.add_argument("--tol-frac", type=float, default=0.002)
    p.add_argument("--margin", type=float, default=0.05)
    p.add_argument("--max-extra-bp", type=int, default=50000)
    p.add_argument("--threads", type=int, default=32,
                   help="Number of worker processes for sequence extraction (default 32)")
    return p.parse_args()


ANCHOR_RE = re.compile(r"^(.+)__(up|down)(\d+)$")


def revcomp(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


# ---- multiprocessing worker for parallel record build ----------------------
_FA = None  # per-worker pysam.FastaFile handle


def _worker_init(genome_db):
    """Open the indexed FASTA once per worker process (not fork-safe to share)."""
    global _FA
    _FA = pysam.FastaFile(genome_db)


def _fetch(contig, start_0, end_0):
    """Fetch [start_0, end_0) from pysam handle; return uppercase or None."""
    if end_0 <= start_0:
        return ""
    try:
        seq = _FA.fetch(contig, max(0, start_0), end_0)
    except (KeyError, ValueError):
        return None
    return seq.upper() if seq else None


def _build_one_record(packed):
    """Build a single IS record. Returns tuple or None on failure."""
    ((ref_id, info, best, src), genome_db, flank_out_len, top_empty_n, db_label,
     max_align_cells) = packed

    # Compute IS size from the peaks (used only by the fallback)
    if best["v0_peak"] and best["empty_peak"]:
        is_length = best["v0_peak"][0] - best["empty_peak"][0]
    elif best["empty_peak"]:
        is_length = best["expected_v0"] - best["empty_peak"][0]
    else:
        is_length = info["tnp_len"]

    # Need full contig length to clamp; pysam provides .get_reference_length
    try:
        contig_len = _FA.get_reference_length(src["contig"])
    except (KeyError, ValueError):
        return None

    tnp_start_0 = src["start"] - 1
    tnp_end_0 = src["end"]

    dec = None
    if best.get("decompose"):
        dec = decompose_boundary(src, info, best["decompose"], top_empty_n,
                                 max_align_cells)
    if dec and dec["ok"]:
        method = "empty_vs_filled_decomposition"
        is_start_0, is_end_0 = dec["is_start_0"], dec["is_end_0"]
        boundary_5p_offset, boundary_3p_offset = dec["off5"], dec["off3"]
    else:
        # FALLBACK: the peak difference only gives a length, so the extension
        # is split evenly between the two sides. It is not a boundary call.
        method = "peak_symmetric_fallback"
        ext_per_side = (is_length - info["tnp_len"]) // 2
        boundary_5p_offset = -ext_per_side
        boundary_3p_offset = is_length - info["tnp_len"] - ext_per_side
        is_start_0 = max(0, tnp_start_0 + boundary_5p_offset)
        is_end_0 = min(contig_len, tnp_end_0 + boundary_3p_offset)

    up_flank_start = max(0, is_start_0 - flank_out_len)
    up_flank_end = is_start_0
    down_flank_start = is_end_0
    down_flank_end = min(contig_len, is_end_0 + flank_out_len)

    tnp_seq = _fetch(src["contig"], tnp_start_0, tnp_end_0)
    is_seq = _fetch(src["contig"], is_start_0, is_end_0)
    up_seq = _fetch(src["contig"], up_flank_start, up_flank_end)
    down_seq = _fetch(src["contig"], down_flank_start, down_flank_end)
    if tnp_seq is None or is_seq is None or up_seq is None or down_seq is None:
        return None

    if src["strand"] == "-":
        tnp_seq = revcomp(tnp_seq)
        is_seq = revcomp(is_seq)
        up_seq, down_seq = revcomp(down_seq), revcomp(up_seq)

    empty_obs_sorted = sorted(best["empty_obs"], key=lambda x: -x["mean_ident"])
    top_empty = empty_obs_sorted[:top_empty_n]

    # Empty junction sequences
    empties_out = []
    for i, eo in enumerate(top_empty):
        target = eo.get("tname") or eo.get("contig")
        if not target:
            continue
        if eo["strand"] == "+":
            junction = _fetch(target, eo["up_te"], eo["down_ts"])
        else:
            raw = _fetch(target, eo["down_te"], eo["up_ts"])
            junction = revcomp(raw) if raw else None
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
            "is_element_length_inferred": is_length,
            "n_v0_observations": len(best["v0_obs"]),
            "n_empty_observations": len(best["empty_obs"]),
            "n_v1plus_observations": len(best["v1plus_obs"]),
            "n_ambiguous_placements": best.get("n_ambiguous", 0),
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
    return (record, is_seq, tnp_seq, up_seq, down_seq, empties_out)


def load_fasta_index(fna):
    contigs = {}
    name = None
    parts = []
    with open(fna) as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    contigs[name] = "".join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.strip())
        if name:
            contigs[name] = "".join(parts)
    return contigs


def find_genome_fna(genome_dir, assembly):
    asm_dir = os.path.join(genome_dir, assembly)
    if not os.path.isdir(asm_dir):
        return None
    for f in os.listdir(asm_dir):
        if f.endswith("_genomic.fna"):
            return os.path.join(asm_dir, f)
    return None


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


def pair_anchors(hits, ref_id, D, info, opts):
    """One scored placement per distinct assembly (lib_alleles.place_pair).

    Replaces an every-up-x-every-down loop that counted a genome once per
    anchor-copy combination and could pair anchors on different repeat
    copies. Ambiguous placements are counted but never used as evidence."""
    expected_v0 = 2 * D + info["tnp_len"]
    placements = place_per_assembly(
        hits.get((ref_id, "up", D), {}), hits.get((ref_id, "down", D), {}),
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
    for p in pairs:
        cat = classify_gap(p["distance"], expected_v0, tnp_len, tol)
        if cat == "v0_filled":
            v0_obs.append(p)
        elif cat == "empty":
            empty_obs.append(p)
        elif cat == "v1plus_filled":
            v1plus_obs.append(p)
    v0_peak = length_mode([p["distance"] for p in v0_obs])
    empty_peak = length_mode([p["distance"] for p in empty_obs])
    return v0_peak, empty_peak, v0_obs, empty_obs, v1plus_obs


def pick_best_d(hits, ref_id, info, distances, opts):
    """Find the best D for this ref: prefer smallest D with both V0 and empty.

    Also returns, as best['decompose'], the smallest D <= decompose_max_d that
    has empty sites -- the interval used for base-level decomposition."""
    candidates = []
    for D in distances:
        pairs, expected_v0, n_amb = pair_anchors(hits, ref_id, D, info, opts)
        if not pairs:
            continue
        v0_peak, empty_peak, v0_obs, empty_obs, v1plus_obs = categorize_pairs(
            pairs, expected_v0, info["tnp_len"], opts)
        candidates.append({
            "D": D, "pairs": pairs, "expected_v0": expected_v0,
            "n_ambiguous": n_amb,
            "v0_peak": v0_peak, "empty_peak": empty_peak,
            "v0_obs": v0_obs, "empty_obs": empty_obs, "v1plus_obs": v1plus_obs,
        })

    if not candidates:
        return None

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
    best["decompose"] = None
    for c in sorted(candidates, key=lambda x: x["D"]):
        if c["D"] <= opts["decompose_max_d"] and c["empty_obs"]:
            best["decompose"] = {"D": c["D"], "empty_obs": c["empty_obs"]}
            break
    return best


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
    filled = _fetch(src["contig"], g0, g1)
    if not filled or len(filled) != g1 - g0:
        return None
    if src["strand"] == "-":
        filled = revcomp(filled)
    cells_ok = lambda e: len(e) * len(filled) <= max_align_cells  # noqa: E731

    votes = Counter()
    calls = []
    methods = Counter()
    n_try = 0
    for eo in sorted(dec["empty_obs"], key=lambda x: -x["mean_ident"])[:top_n]:
        n_try += 1
        if eo["distance"] <= 0:
            empty = ""
        elif eo["strand"] == "+":
            empty = _fetch(eo["tname"], eo["up_te"], eo["down_ts"])
        else:
            raw = _fetch(eo["tname"], eo["down_te"], eo["up_ts"])
            empty = revcomp(raw) if raw else None
        if empty is None:
            continue
        r = decompose_alleles(empty, filled, min_insert=max(50, tnp_len - 20),
                              tolerant=cells_ok(empty))
        methods[r["method"]] += 1
        if r["method"] == "none":
            continue
        st, ln = r["insert_start"], r["insert_len"]
        # the element must contain the whole transposase ORF
        if st > D + 10 or st + ln < D + tnp_len - 10:
            methods["insert_misses_tnp"] += 1
            continue
        off5, off3 = st - D, (st + ln) - (D + tnp_len)
        votes[(off5, off3)] += 1
        calls.append((off5, off3, st, ln, r, eo.get("assembly", "")))
    if not votes:
        return {"ok": False, "D": D, "n_attempted": n_try,
                "method_counts": dict(methods)}
    (off5, off3), n_cons = votes.most_common(1)[0]
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
            "method_counts": dict(methods),
            "consensus_method": r["method"], "event_class": r["event_class"],
            "empty_offset": r["offset"],
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
    logan_hits = {}  # legacy compatibility; unused with single-DB design

    # Build records — parallelize per-ref using pysam (in-process FASTA reads)
    n_done = n_neither = 0
    db_label = os.path.basename(args.genome_db)

    # Pre-filter refs to those with a valid pick (cheap; no FASTA I/O)
    work_items = []
    opts = {"tol_bp": args.tol_bp, "tol_frac": args.tol_frac,
            "margin": args.margin, "max_extra_bp": args.max_extra_bp,
            "decompose_max_d": args.decompose_max_d}
    for ref_id, info in ref_info.items():
        best = pick_best_d(ncbi_hits, ref_id, info, distances, opts)
        if not best:
            n_neither += 1
            continue
        src = source_positions.get(ref_id)
        if not src:
            continue
        n_done += 1
        work_items.append((ref_id, info, best, src))

    print(f"  Refs to process: {len(work_items):,} (skipped {n_neither:,} with no usable D)",
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
    with Pool(args.threads, initializer=_worker_init,
              initargs=(args.genome_db,)) as pool:
        for i, packed in enumerate(pool.imap_unordered(_build_one_record, pool_args,
                                                       chunksize=8)):
            if i % 1000 == 0:
                print(f"    built {i:,}/{len(pool_args):,}", file=sys.stderr, flush=True)
            if packed is None:
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
    json_path = os.path.join(args.out, "records.json")
    with open(json_path, "w") as f:
        json.dump(all_records, f, indent=2)
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
                "junction_microhomology_bp\ttsd_len\tevent_key\n")
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
                    f"{(dc.get('tsd') or {}).get('length', '')}\t{r['event_key']}\n")

    print(f"\n=== Build summary ===", file=sys.stderr)
    print(f"  Records built:       {n_done}", file=sys.stderr)
    n_dec = sum(1 for r in all_records
                if r["boundary_evidence"]["method"] == "empty_vs_filled_decomposition")
    print(f"  Boundary by decomposition: {n_dec:,} / {len(all_records):,} "
          f"(rest: peak_symmetric_fallback)", file=sys.stderr)
    print(f"  Distinct insertion events (event_key): "
          f"{len({r['event_key'] for r in all_records}):,}", file=sys.stderr)
    print(f"  No empty found:      {n_neither}", file=sys.stderr)
    print(f"  All outputs in: {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
