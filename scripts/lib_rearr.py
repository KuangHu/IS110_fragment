#!/usr/bin/env python3
"""Shared helpers for the scripts that confirm Stage 7 rearrangement candidates.

One convention drives everything here and was missed by most of the callers:
ANCHORS ARE IN TRANSPOSASE ORIENTATION. extract_anchors.py reverse-complements
the reference for a minus-strand transposase, so

  * the "up" anchor is the GENOMIC RIGHT flank of a minus-strand element;
  * the strand written in `details` (contig:ts-te(strand)) is the target
    strand relative to the transposase-oriented anchor.

So a reference flank ("proxy") must be taken on the transposase-relative side
and oriented with the transposase; aligned to the target it must then lie on
the SAME strand as the corresponding anchor hit. Taking "left of the IS" for
"up" and accepting only '+' alignments failed every minus-strand element and
every target contig stored in the other orientation.
"""
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_alleles import revcomp  # noqa: E402

_SIDE_RE = {s: re.compile(rf"{s}=\[([^\]]*)\]") for s in ("up", "down")}
_HIT_RE = re.compile(r"^(.+):(\d+)-(\d+)\(([+-])\)$")


def parse_details(details):
    """'up=[c:1-2(+);c:5-9(-)] down=[...]' -> {'up': [(contig, s, e, strand)], 'down': [...]}.

    Parsed structurally. The old pattern [\\w.|]+ cut contig names at '-' or
    ':' ('contig-12:100-600(+)' became contig '12')."""
    out = {"up": [], "down": []}
    for side, rx in _SIDE_RE.items():
        m = rx.search(details or "")
        if not m:
            continue
        for tok in m.group(1).split(";"):
            h = _HIT_RE.match(tok.strip())
            if h:
                out[side].append((h.group(1), int(h.group(2)), int(h.group(3)), h.group(4)))
    return out


def load_loci(records_path):
    """ref_id -> (local contig, start0, end, strand): the element as a 0-based
    half-open interval plus the transposase strand. source_start is 1-based
    inclusive in records.json."""
    loci = {}
    with open(records_path) as f:
        recs = json.load(f)
    for r in recs:
        rid = r.get("ref_id") or r.get("is110_id")
        src = r.get("source", {})
        ie = r.get("is_element") or src.get("is_element") or {}
        contig = src.get("contig")
        s = ie.get("source_start") or src.get("transposase_start")
        e = ie.get("source_end") or src.get("transposase_end")
        if rid and contig and s is not None and e is not None:
            local = contig.split("|", 1)[1] if "|" in contig else contig
            loci[rid] = (local, int(s) - 1, int(e),
                         src.get("transposase_strand") or ie.get("strand") or "+")
    return loci


def ref_proxies(ref_fa, contig, s0, e, strand, window):
    """(up, down) reference flanks of `window` bp in TRANSPOSASE orientation:
    `up` ends at the element's 5' end, `down` starts at its 3' end."""
    clen = ref_fa.get_reference_length(contig)
    left = ref_fa.fetch(contig, max(0, s0 - window), s0).upper()
    right = ref_fa.fetch(contig, e, min(clen, e + window)).upper()
    if strand == "-":
        return revcomp(right), revcomp(left)
    return left, right


def flip(strand):
    return "-" if strand == "+" else "+"


def align(query_seq, target_seq, work, threads, extra=()):
    """minimap2 asm10 of query on target. Returns a list of dicts
    (qs, qe, ts, te, strand, ident) with identity = matches / columns.
    Raises on failure -- a tool error must never read as 'no alignment'."""
    d = tempfile.mkdtemp(prefix="al_", dir=work)
    try:
        q, t = os.path.join(d, "q.fa"), os.path.join(d, "t.fa")
        with open(q, "w") as fh:
            fh.write(f">q\n{query_seq}\n")
        with open(t, "w") as fh:
            fh.write(f">t\n{target_seq}\n")
        # --no-long-join -r 500: do not bridge large indels. asm10 otherwise
        # joins across a 13 kb insertion into ONE block whose column identity
        # (82%) fails every threshold, and a duplicated segment never shows up
        # as two blocks. Split blocks + union coverage handle both.
        res = subprocess.run(["minimap2", "-x", "asm10", "-c", "--eqx", "-p", "0.1",
                              "-N", "50", "--no-long-join", "-r", "500",
                              "-t", str(threads), *extra, t, q],
                             check=True, capture_output=True, text=True)
        out = []
        for line in res.stdout.splitlines():
            c = line.split("\t")
            if len(c) < 12:
                continue
            out.append({"qs": int(c[2]), "qe": int(c[3]), "ts": int(c[7]), "te": int(c[8]),
                        "strand": c[4],
                        "ident": int(c[9]) / int(c[10]) * 100 if int(c[10]) else 0.0})
        return out
    finally:
        shutil.rmtree(d, ignore_errors=True)


def union_bp(intervals):
    total, cs, ce = 0, None, None
    for s, e in sorted(intervals):
        if ce is None or s > ce:
            if ce is not None:
                total += ce - cs
            cs, ce = s, e
        else:
            ce = max(ce, e)
    return total + (ce - cs if ce is not None else 0)


def strand_cover(blocks, strand, min_ident, key="q"):
    """bp of the query (key='q') or target ('t') covered by blocks on `strand`
    at >= min_ident: a UNION, so one indel or IS inside a 30 kb flank does not
    reduce 'the single largest block' below the threshold."""
    a, b = (key + "s", key + "e")
    return union_bp([(x[a], x[b]) for x in blocks
                     if x["strand"] == strand and x["ident"] >= min_ident])


def fetch_window(fa, contig, s, e, window):
    """Target window around [s, e); returns (sequence, offset) or ('', 0)."""
    try:
        clen = fa.get_reference_length(contig)
    except (KeyError, ValueError):
        return "", 0
    a, b = max(0, s - window), min(clen, e + window)
    return fa.fetch(contig, a, b).upper(), a


def copy_depth(blocks, min_block, min_ident, strand=None):
    """How many separate places in the target carry the same stretch of the
    query? blocks: output of align(). Blocks that overlap ON THE TARGET are one
    place (the longer is kept). Returns (dup_bp, max_copies, places):
      dup_bp      query bp present at >= 2 places
      max_copies  largest number of places covering one stretch >= min_block
      places      the blocks used
    This is the test for 'the segment appears twice'. Counting forward blocks
    >= 1 kb anywhere in the window (what the duplication check did) passed a
    target with no duplicated segment at all."""
    use = sorted((b for b in blocks if b["ident"] >= min_ident
                  and b["qe"] - b["qs"] >= min_block
                  and (strand is None or b["strand"] == strand)),
                 key=lambda b: -(b["qe"] - b["qs"]))
    places = []
    for b in use:
        if not any(min(b["te"], p["te"]) - max(b["ts"], p["ts"]) > 0 for p in places):
            places.append(b)
    pts = sorted([(p["qs"], 1) for p in places] + [(p["qe"], -1) for p in places])
    dup_bp, max_copies, depth, last = 0, 0, 0, None
    for pos, d in pts:
        if last is not None and pos > last:
            if depth >= 2:
                dup_bp += pos - last
            if pos - last >= min_block:
                max_copies = max(max_copies, depth)
        depth += d
        last = pos
    return dup_bp, max_copies, places
