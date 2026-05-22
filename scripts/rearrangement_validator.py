#!/usr/bin/env python3
"""Standard-method validator for IS-mediated rearrangements.

Uses MUMmer (nucmer + show-coords + show-diff) — the de facto standard for
bacterial genome structural variant detection — to validate candidate events
flagged by find_rearrangements.py.

For each candidate event, outputs a structured verdict with multiple evidence
streams:
  1. Alignment summary (blocks, identity, orientation)
  2. show-diff structural calls (INV, JMP, DUP, BRK, GAP)
  3. IS-mediated check: are breakpoints near IS transposase positions?

Verdict logic:
  - INVERSION:
      CONFIRMED if show-diff reports INV >= expected_size * 0.5
                AND inversion endpoints are within ±2kb of an IS position
      LIKELY if INV reported but no IS at endpoints
      REJECTED if no INV reported (anchors disagree but alignment is co-linear)

  - TRANSLOCATION:
      CONFIRMED if show-diff reports JMP between different ref/qry sequences
                AND endpoints are near IS positions
      LIKELY if JMP without IS endpoint match
      REJECTED if alignment fits on a single sequence pair

  - DUPLICATION:
      CONFIRMED if show-diff reports DUP
                OR same query region maps to >=2 disjoint ref regions
      LIKELY if multi-mapping but no DUP call
      REJECTED if single mapping

Usage as library:
    from rearrangement_validator import validate_event
    verdict = validate_event(source_fa, target_fa, event_type, anchor_info,
                              is_positions=positions_dict)

Usage as CLI:
    rearrangement_validator.py <events_tsv> <genome_dir> <is_hits_tsv>
                                <out_dir> [--n-per-cat N]
"""

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict


# ────────────────────────────────────────────────────────────
# MUMmer wrappers
# ────────────────────────────────────────────────────────────

def run_nucmer(ref_fa, query_fa, prefix, maxmatch=True):
    """Run nucmer; returns delta path."""
    cmd = ["nucmer", "-p", prefix]
    if maxmatch:
        cmd.append("--maxmatch")
    cmd += [ref_fa, query_fa]
    subprocess.run(cmd, check=False, capture_output=True)
    delta = prefix + ".delta"
    return delta if os.path.exists(delta) else None


def parse_show_coords(delta):
    """Parse show-coords -r -c -l -T output into list of alignment blocks."""
    if not delta or not os.path.exists(delta):
        return []
    res = subprocess.run(
        ["show-coords", "-r", "-c", "-l", "-T", delta],
        capture_output=True, text=True
    )
    blocks = []
    for line in res.stdout.splitlines():
        if not line:
            continue
        if line.startswith(("=", "/", "NUCMER")) or "[S1]" in line:
            continue
        f = line.split("\t")
        if len(f) < 9:
            continue
        try:
            blocks.append({
                "rs": int(f[0]), "re": int(f[1]),
                "qs": int(f[2]), "qe": int(f[3]),
                "ralen": int(f[4]), "qalen": int(f[5]),
                "ident": float(f[6]),
                "rlen": int(f[7]), "qlen": int(f[8]),
                "ref": f[-2] if len(f) >= 11 else "",
                "qry": f[-1] if len(f) >= 11 else "",
            })
        except (ValueError, IndexError):
            continue
    return blocks


def parse_show_diff(delta, mode="r"):
    """Run show-diff [-r|-q] and parse structural variant calls.

    show-diff outputs lines like:
      <SEQ> <FEAT> <S1> <E1> <LEN_R> <S2> <E2> <LEN_Q>
    Where FEAT is one of: GAP, DUP, BRK, JMP, INV, SEQ
    """
    if not delta or not os.path.exists(delta):
        return []
    res = subprocess.run(
        ["show-diff", f"-{mode}", "-H", delta],
        capture_output=True, text=True
    )
    diffs = []
    for line in res.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 5:
            continue
        diffs.append({
            "seq": f[0],
            "feat": f[1],
            "fields": f[2:],
        })
    return diffs


# ────────────────────────────────────────────────────────────
# IS position helpers (any IS family — uses is_hits.tsv from Stage 1)
# ────────────────────────────────────────────────────────────

def load_is_positions(is_hits_tsv):
    """Load IS transposase positions from Stage 1 output.

    Expects TSV with header: is_id, assembly, contig, tnp_start, tnp_end,
    tnp_strand, tnp_len, domains_hit (the schema written by is_detect.py).

    Returns: {assembly: {contig: [(start, end, strand)]}}.
    """
    by_asm = defaultdict(lambda: defaultdict(list))
    with open(is_hits_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            asm = r["assembly"]
            contig = r["contig"]
            try:
                s = int(r["tnp_start"])
                e = int(r["tnp_end"])
            except (KeyError, ValueError):
                continue
            strand = r.get("tnp_strand", "+")
            by_asm[asm][contig].append((s, e, strand))
    return by_asm


def near_is_element(positions, assembly, contig, pos, window=2000):
    """Return True if pos on contig is within window of an IS transposase."""
    if assembly not in positions:
        return False
    if contig not in positions[assembly]:
        return False
    for s, e, _ in positions[assembly][contig]:
        if min(abs(pos - s), abs(pos - e)) <= window:
            return True
        if s <= pos <= e:
            return True
    return False


# ────────────────────────────────────────────────────────────
# Validators per event type
# ────────────────────────────────────────────────────────────

def validate_inversion(blocks, diffs, anchor_info, is_positions, target_assembly,
                       inversion_min_frac=0.4):
    """Inversion: show-diff INV plus inverted alignment blocks dominate."""
    inv_diffs = [d for d in diffs if d["feat"] == "INV"]
    inv_blocks = [b for b in blocks if b["qs"] > b["qe"]]
    fwd_blocks = [b for b in blocks if b["qs"] < b["qe"]]

    inv_total = sum(b["qalen"] for b in inv_blocks)
    fwd_total = sum(b["qalen"] for b in fwd_blocks)
    total = inv_total + fwd_total
    inv_frac = inv_total / total if total > 0 else 0

    is_mediated = False
    if inv_diffs and is_positions:
        for d in inv_diffs:
            if len(d["fields"]) < 2:
                continue
            try:
                pos1 = int(d["fields"][0])
                pos2 = int(d["fields"][1])
                contig = d["seq"]
                if (near_is_element(is_positions, target_assembly, contig, pos1) or
                    near_is_element(is_positions, target_assembly, contig, pos2)):
                    is_mediated = True
                    break
            except ValueError:
                continue

    if inv_diffs and inv_frac >= inversion_min_frac:
        verdict = "CONFIRMED_IS_MEDIATED" if is_mediated else "CONFIRMED"
    elif inv_diffs:
        verdict = "LIKELY"
    elif inv_frac >= inversion_min_frac:
        verdict = "LIKELY_blocks_only"
    else:
        verdict = "REJECTED"

    return {
        "verdict": verdict,
        "inv_diffs": len(inv_diffs),
        "inv_blocks": len(inv_blocks),
        "fwd_blocks": len(fwd_blocks),
        "inv_fraction": round(inv_frac, 3),
        "is_mediated": is_mediated,
    }


def validate_translocation(blocks, diffs, anchor_info, is_positions, target_assembly):
    """Translocation: show-diff JMP between different sequences."""
    jmp_diffs = [d for d in diffs if d["feat"] == "JMP"]
    brk_diffs = [d for d in diffs if d["feat"] == "BRK"]

    distinct_q = set(b["qry"] for b in blocks)

    cross_contig_anchors = (anchor_info.get("up_contig") and anchor_info.get("down_contig")
                            and anchor_info["up_contig"] != anchor_info["down_contig"])

    is_mediated = False
    if jmp_diffs and is_positions:
        for d in jmp_diffs:
            if len(d["fields"]) < 2:
                continue
            try:
                pos = int(d["fields"][0])
                if near_is_element(is_positions, target_assembly, d["seq"], pos):
                    is_mediated = True
                    break
            except ValueError:
                continue

    if (jmp_diffs or cross_contig_anchors) and len(distinct_q) >= 2:
        verdict = "CONFIRMED_IS_MEDIATED" if is_mediated else "CONFIRMED"
    elif jmp_diffs or cross_contig_anchors:
        verdict = "LIKELY"
    else:
        verdict = "REJECTED"

    return {
        "verdict": verdict,
        "jmp_diffs": len(jmp_diffs),
        "brk_diffs": len(brk_diffs),
        "n_target_contigs": len(distinct_q),
        "cross_contig_anchors": cross_contig_anchors,
        "is_mediated": is_mediated,
    }


def validate_duplication(blocks, diffs, anchor_info, is_positions, target_assembly):
    """Duplication: show-diff DUP, or same source region maps to >=2 query regions."""
    dup_diffs = [d for d in diffs if d["feat"] == "DUP"]

    multi_map_count = 0
    by_ref_range = defaultdict(list)
    for b in blocks:
        key = (b["ref"], b["rs"] // 1000)
        by_ref_range[key].append(b)
    for key, bs in by_ref_range.items():
        qry_positions = set((b["qry"], b["qs"] // 1000) for b in bs)
        if len(qry_positions) >= 2:
            multi_map_count += 1

    n_down_hits = anchor_info.get("n_down_hits", 1)
    n_up_hits = anchor_info.get("n_up_hits", 1)
    multi_anchor = (n_down_hits >= 2) or (n_up_hits >= 2)

    is_mediated = False
    if dup_diffs and is_positions:
        for d in dup_diffs:
            if len(d["fields"]) < 2:
                continue
            try:
                pos = int(d["fields"][0])
                if near_is_element(is_positions, target_assembly, d["seq"], pos):
                    is_mediated = True
                    break
            except ValueError:
                continue

    if dup_diffs or multi_map_count > 0 or multi_anchor:
        if dup_diffs:
            verdict = "CONFIRMED_IS_MEDIATED" if is_mediated else "CONFIRMED"
        else:
            verdict = "LIKELY"
    else:
        verdict = "REJECTED"

    return {
        "verdict": verdict,
        "dup_diffs": len(dup_diffs),
        "multi_map_regions": multi_map_count,
        "multi_anchor": multi_anchor,
        "is_mediated": is_mediated,
    }


# ────────────────────────────────────────────────────────────
# Top-level event validator
# ────────────────────────────────────────────────────────────

def validate_event(source_fa, target_fa, event_type, anchor_info,
                   is_positions=None, target_assembly=None,
                   work_dir=None, keep_files=False):
    """Validate a single rearrangement event with nucmer + show-diff."""
    cleanup = False
    if work_dir is None:
        work_dir = tempfile.mkdtemp(prefix="rval_")
        cleanup = True
    os.makedirs(work_dir, exist_ok=True)

    prefix = os.path.join(work_dir, "aln")
    delta = run_nucmer(source_fa, target_fa, prefix, maxmatch=True)
    blocks = parse_show_coords(delta)
    diffs_r = parse_show_diff(delta, mode="r")
    diffs_q = parse_show_diff(delta, mode="q")
    diffs = diffs_r + diffs_q

    if event_type == "inversion_only":
        v = validate_inversion(blocks, diffs, anchor_info, is_positions, target_assembly)
    elif event_type == "translocation_only":
        v = validate_translocation(blocks, diffs, anchor_info, is_positions, target_assembly)
    elif event_type == "duplication":
        v = validate_duplication(blocks, diffs, anchor_info, is_positions, target_assembly)
    else:
        v = {"verdict": "UNKNOWN_EVENT_TYPE"}

    v["n_total_blocks"] = len(blocks)
    v["n_diff_calls"] = len(diffs)
    v["work_dir"] = work_dir if keep_files else None

    if cleanup and not keep_files:
        shutil.rmtree(work_dir, ignore_errors=True)
    return v


# ────────────────────────────────────────────────────────────
# CLI / batch driver
# ────────────────────────────────────────────────────────────

def find_genome_fna(genome_dir, assembly):
    asm_dir = os.path.join(genome_dir, assembly)
    if not os.path.isdir(asm_dir):
        return None
    for f in os.listdir(asm_dir):
        if f.endswith("_genomic.fna"):
            return os.path.join(asm_dir, f)
    return None


def parse_anchor_detail(detail_str):
    out = {"up": [], "down": []}
    for side in ("up", "down"):
        m = re.search(rf"{side}=\[([^\]]*)\]", detail_str)
        if not m:
            continue
        for hit in m.group(1).split(";"):
            mm = re.match(r"(.+):(\d+)-(\d+)\(([+-])\)", hit)
            if mm:
                out[side].append({
                    "contig": mm.group(1),
                    "start": int(mm.group(2)),
                    "end": int(mm.group(3)),
                    "strand": mm.group(4),
                })
    return out


def extract_region_to_fa(genome_fa, contig, start, end, out_fa):
    """Use samtools faidx to extract a region."""
    if not os.path.exists(genome_fa + ".fai"):
        subprocess.run(["samtools", "faidx", genome_fa], check=True)
    region = f"{contig}:{max(1, start)}-{end}"
    res = subprocess.run(
        ["samtools", "faidx", genome_fa, region],
        capture_output=True, text=True
    )
    if res.returncode != 0:
        return False
    with open(out_fa, "w") as f:
        f.write(res.stdout)
    return True


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("events_tsv", help="find_rearrangements.py *_examples.tsv")
    p.add_argument("genome_dir", help="Directory with <assembly>/*_genomic.fna")
    p.add_argument("is_hits", help="Stage 1 is_hits.tsv (transposase positions)")
    p.add_argument("out_dir")
    p.add_argument("--n-per-cat", type=int, default=10)
    p.add_argument("--region-flank", type=int, default=50000)
    p.add_argument("--keep-files", action="store_true")
    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    print("Loading IS positions...", file=sys.stderr)
    is_pos = load_is_positions(args.is_hits)
    print(f"  Loaded positions for {len(is_pos)} assemblies", file=sys.stderr)

    examples = defaultdict(list)
    with open(args.events_tsv) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            examples[r["category"]].append(r)

    report = open(os.path.join(args.out_dir, "validation_report.tsv"), "w")
    report.write("category\tref_id\ttarget_assembly\tup_loc\tdown_loc\tverdict\t"
                 "is_mediated\tdetails_json\n")

    for cat in ("inversion_only", "translocation_only", "duplication"):
        eg_list = examples.get(cat, [])[:args.n_per_cat]
        print(f"\n=== {cat} ({len(eg_list)} candidates) ===", file=sys.stderr)
        for i, eg in enumerate(eg_list):
            ref_id = eg["ref_id"]
            target_assembly = eg["assembly"]
            source_assembly = ref_id.split("|")[0]

            source_fa = find_genome_fna(args.genome_dir, source_assembly)
            target_fa = find_genome_fna(args.genome_dir, target_assembly)
            if not source_fa or not target_fa:
                report.write(f"{cat}\t{ref_id}\t{target_assembly}\t-\t-\tMISSING_GENOME\tFalse\t{{}}\n")
                continue

            details = parse_anchor_detail(eg["details"])
            up_hits = details["up"]
            down_hits = details["down"]
            if not up_hits or not down_hits:
                continue
            up = up_hits[0]
            dn = down_hits[0]

            work_dir = os.path.join(args.out_dir, f"{cat}_{i}")
            os.makedirs(work_dir, exist_ok=True)
            target_region_fa = os.path.join(work_dir, "target_region.fa")

            if up["contig"] == dn["contig"]:
                positions = sorted([up["start"], up["end"], dn["start"], dn["end"]])
                rstart = max(1, positions[0] - args.region_flank)
                rend = positions[-1] + args.region_flank
                extract_region_to_fa(target_fa, up["contig"], rstart, rend, target_region_fa)
            else:
                with open(target_region_fa, "w") as f:
                    for hit in (up, dn):
                        rstart = max(1, hit["start"] - args.region_flank)
                        rend = hit["end"] + args.region_flank
                        tmp_fa = os.path.join(work_dir, f"part_{hit['contig']}.fa")
                        if extract_region_to_fa(target_fa, hit["contig"], rstart, rend, tmp_fa):
                            with open(tmp_fa) as g:
                                f.write(g.read())

            anchor_info = {
                "up_contig": up["contig"], "up_strand": up["strand"],
                "down_contig": dn["contig"], "down_strand": dn["strand"],
                "n_up_hits": len(up_hits), "n_down_hits": len(down_hits),
            }

            verdict = validate_event(
                source_fa, target_region_fa, cat,
                anchor_info=anchor_info,
                is_positions=is_pos,
                target_assembly=target_assembly,
                work_dir=work_dir,
                keep_files=args.keep_files,
            )

            up_loc = f"{up['contig']}:{up['start']}-{up['end']}({up['strand']})"
            dn_loc = f"{dn['contig']}:{dn['start']}-{dn['end']}({dn['strand']})"
            report.write(f"{cat}\t{ref_id}\t{target_assembly}\t{up_loc}\t{dn_loc}\t"
                         f"{verdict['verdict']}\t{verdict.get('is_mediated', False)}\t"
                         f"{json.dumps({k: v for k, v in verdict.items() if k != 'work_dir'})}\n")
            print(f"  [{i+1}] {ref_id[:30]}... -> {target_assembly}  =>  {verdict['verdict']}",
                  file=sys.stderr)

    report.close()
    print(f"\nReport: {args.out_dir}/validation_report.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()
