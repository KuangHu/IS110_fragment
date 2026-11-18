#!/usr/bin/env python3
"""Shared PAF helpers for the IS110_fragment detectors.

Identity
    PAF column 11 (alignment block length) INCLUDES gap columns. minimap2
    bridges moderate indels inside one alignment, so `matches / block` for an
    IS copy with a 600-bp insertion is 71% and the row was thrown away by the
    95% filter. Identity here is matches / aligned columns, where the aligned
    columns are the M/=/X operations of the CIGAR. (Not `block - gaps`:
    minimap2 leaves ambiguous bases out of the block length, also those inside
    a gap, so that difference can be smaller than the number of matches.)

    Every reader REQUIRES a CIGAR with =/X operations (minimap2 -c --eqx):
    a row without one raises an error. The earlier fallbacks (1 - de, or
    matches / block) were approximations and are gone -- every minimap2 call
    in this repository uses --eqx.

Coordinates
    PAF is 0-based half-open. Query coordinates are always on the query's
    forward strand; on a '-' row the target-left end of the alignment is the
    query's right end.

Command line (replaces the `awk '$10/$11 >= min'` filters in the shell scripts):
    paf_utils.py filter --min-identity 95 [--query NAME] < in.paf > out.paf
"""
import re
import sys

_CG = re.compile(r"(\d+)([MIDNSHP=X])")


def cigar_of(cols):
    for t in cols[12:]:
        if t.startswith("cg:Z:"):
            return t[5:]
    return ""


def require_eqx(cols):
    """The row's CIGAR; raises unless it is an =/X CIGAR."""
    cg = cigar_of(cols)
    if not cg or "M" in cg:
        raise ValueError(f"PAF row without an =/X CIGAR (query {cols[0]}): "
                         "run minimap2 with -c --eqx")
    return cg


def aligned_identity(cols):
    """Percent identity over aligned columns: '=' / ('=' + 'X'). Exact."""
    eq = x = 0
    for n, op in _CG.findall(require_eqx(cols)):
        if op == "=":
            eq += int(n)
        elif op == "X":
            x += int(n)
    return 100.0 * eq / (eq + x) if eq + x else 0.0


def aligned_blocks(cols, min_target_only=None, min_query_only=None):
    """Cut one PAF row at its large indels.

    Returns a list of (qs, qe, ts, te) in TARGET order; qs/qe on the query's
    forward strand. A 'D' run (bases only in the target) of at least
    `min_target_only`, or an 'I' run (bases only in the query) of at least
    `min_query_only`, ends a block; smaller indels stay inside it. With both
    limits None, or without a CIGAR, the row is one block.
    The bases of a large 'I' belong to NO block: they are not aligned."""
    qs, qe, strand, ts, te = int(cols[2]), int(cols[3]), cols[4], int(cols[7]), int(cols[8])
    cg = require_eqx(cols)
    if min_target_only is None and min_query_only is None:
        return [(qs, qe, ts, te)]
    plus = strand == "+"
    q = qs if plus else qe          # query position at the target-left end
    t = ts
    blocks = []
    bq, bt = q, t                   # start of the open block
    for n, op in _CG.findall(cg):
        n = int(n)
        cut = ((op == "D" and min_target_only is not None and n >= min_target_only) or
               (op == "I" and min_query_only is not None and n >= min_query_only))
        if cut:
            if t > bt:
                blocks.append((min(bq, q), max(bq, q), bt, t))
            if op == "D":
                t += n
            else:
                q += n if plus else -n
            bq, bt = q, t
            continue
        if op in "M=XDN":
            t += n
        if op in "M=XI":
            q += n if plus else -n
    if t > bt:
        blocks.append((min(bq, q), max(bq, q), bt, t))
    return blocks


def aligned_query_bases(cols):
    """Query bases that are aligned to the target (CIGAR M/=/X). A query span
    qe - qs also contains the 'I' runs, i.e. query sequence that is ABSENT
    from the target. Without a CIGAR the span is returned."""
    return sum(int(n) for n, op in _CG.findall(require_eqx(cols)) if op in "=X")


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


def _main(argv):
    import argparse
    p = argparse.ArgumentParser(description="Filter a PAF by aligned identity")
    p.add_argument("cmd", choices=["filter"])
    p.add_argument("--min-identity", type=float, required=True,
                   help="percent; a value <= 1 is read as a fraction")
    p.add_argument("--query", default=None, help="keep only this query name")
    a = p.parse_args(argv)
    thr = a.min_identity * 100 if a.min_identity <= 1 else a.min_identity
    for line in sys.stdin:
        c = line.rstrip("\n").split("\t")
        if len(c) < 12 or (a.query is not None and c[0] != a.query):
            continue
        if aligned_identity(c) >= thr:
            sys.stdout.write(line)


if __name__ == "__main__":
    _main(sys.argv[1:])
