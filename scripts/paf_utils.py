#!/usr/bin/env python3
"""Shared PAF helpers for the IS110_fragment detectors.

Identity: PAF column 11 (alignment block length) INCLUDES gap columns. minimap2
bridges moderate indels inside one alignment, so `matches / block` for an IS
copy with a 600-bp insertion is 71% and the row was thrown away by the 95%
filter -- the event was neither called nor counted as an intact copy. Identity
is taken over aligned (non-gap) columns here.
"""
import re

_CG = re.compile(r"(\d+)([MIDNSHP=X])")


def cigar_of(cols):
    for t in cols[12:]:
        if t.startswith("cg:Z:"):
            return t[5:]
    return ""


def aligned_identity(cols):
    """Percent identity over aligned columns (gaps excluded)."""
    matches, block = int(cols[9]), int(cols[10])
    cg = cigar_of(cols)
    if cg:
        gaps = sum(int(n) for n, op in _CG.findall(cg) if op in "IDN")
        aligned = block - gaps
        return 100.0 * matches / aligned if aligned > 0 else 0.0
    for t in cols[12:]:
        if t.startswith("de:f:"):          # gap-compressed divergence
            return 100.0 * (1.0 - float(t[5:]))
    return 100.0 * matches / block if block > 0 else 0.0


def large_indels(cols, min_len):
    """Indels >= min_len inside one alignment row.

    Yields (kind, length, q_pos, t_pos): kind 'D' = bases only in the TARGET
    (an insertion in the genome when the IS is the query), 'I' = bases only in
    the QUERY (a deletion in the genome). q_pos / t_pos are 0-based positions
    of the gap on the query's forward strand and on the target."""
    cg = cigar_of(cols)
    if not cg:
        return
    qlen, qs, qe, strand, ts = int(cols[1]), int(cols[2]), int(cols[3]), cols[4], int(cols[7])
    q, t = (qs, ts) if strand == "+" else (qe, ts)
    for n, op in _CG.findall(cg):
        n = int(n)
        if op in "ID" and n >= min_len:
            if op == "D":
                yield ("D", n, q, t)
            else:
                yield ("I", n, q if strand == "+" else q - n, t)
        if op in "M=XD":
            t += n
        if op in "M=XI":
            q += n if strand == "+" else -n
