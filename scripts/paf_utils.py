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
    """The row's CIGAR; raises unless it is an =/X CIGAR.

    'N' operations (reference skips of the splice presets) are refused as
    well: the walkers below would otherwise step over them as if they were
    aligned target bases. The assembly presets used here do not emit them."""
    cg = cigar_of(cols)
    if not cg or "M" in cg:
        raise ValueError(f"PAF row without an =/X CIGAR (query {cols[0]}): "
                         "run minimap2 with -c --eqx")
    if "N" in cg:
        raise ValueError(f"PAF row with an 'N' operation in its CIGAR (query {cols[0]}): "
                         "not supported, use an assembly preset")
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


def eqx_runs(cols):
    """Every '=' / 'X' run of the row as (q_lo, q_hi, t_lo, t_hi, op), in
    target order; q on the query's forward strand. These are the ALIGNED
    columns: no indel, of whatever size, is inside a run."""
    qs, qe, strand, ts = int(cols[2]), int(cols[3]), cols[4], int(cols[7])
    plus = strand == "+"
    q = qs if plus else qe
    t = ts
    for n, op in _CG.findall(require_eqx(cols)):
        n = int(n)
        if op in "=X":
            yield ((q, q + n) if plus else (q - n, q)) + (t, t + n, op)
            t += n
            q += n if plus else -n
        elif op == "D":
            t += n
        elif op == "I":
            q += n if plus else -n


def _add_iv(ivs, lo, hi):
    """Append a query interval, joining it to the last one when they touch."""
    if ivs and (ivs[-1][1] == lo or ivs[-1][0] == hi):
        ivs[-1] = (min(ivs[-1][0], lo), max(ivs[-1][1], hi))
    else:
        ivs.append((lo, hi))


def aligned_blocks_iv(cols, min_target_only=None, min_query_only=None):
    """Cut one PAF row at its large indels.

    Returns a list of (qs, qe, ts, te, intervals) in TARGET order; qs/qe on
    the query's forward strand. A 'D' run (bases only in the target) of at
    least `min_target_only`, or an 'I' run (bases only in the query) of at
    least `min_query_only`, ends a block; smaller indels stay inside it.
    The bases of a large 'I' belong to NO block.

    `intervals` are the query intervals of the block that are really aligned
    ('=' / 'X'). An 'I' run BELOW the cutting limit stays inside qs..qe but is
    in no interval: coverage must be counted from the intervals, not from
    qe - qs (three 100-bp 'I' runs made a 90%-aligned copy "100% covered")."""
    cg = require_eqx(cols)
    strand = cols[4]
    plus = strand == "+"
    q = int(cols[2]) if plus else int(cols[3])   # query position at the target-left end
    t = int(cols[7])
    blocks = []
    bq, bt, ivs = q, t, []                       # the open block

    def close():
        if t > bt:
            blocks.append((min(bq, q), max(bq, q), bt, t, sorted(ivs)))

    for n, op in _CG.findall(cg):
        n = int(n)
        cut = ((op == "D" and min_target_only is not None and n >= min_target_only) or
               (op == "I" and min_query_only is not None and n >= min_query_only))
        if cut:
            close()
            if op == "D":
                t += n
            else:
                q += n if plus else -n
            bq, bt, ivs = q, t, []
            continue
        if op in "=X":
            _add_iv(ivs, *((q, q + n) if plus else (q - n, q)))
        if op in "=XD":
            t += n
        if op in "=XI":
            q += n if plus else -n
    close()
    return blocks


def aligned_blocks(cols, min_target_only=None, min_query_only=None):
    """aligned_blocks_iv without the interval lists: (qs, qe, ts, te)."""
    return [b[:4] for b in aligned_blocks_iv(cols, min_target_only, min_query_only)]


def union_length(intervals):
    """Number of positions covered by a list of half-open intervals."""
    covered, end = 0, None
    for s, e in sorted(intervals):
        if end is None or s > end:
            covered += e - s
            end = e
        elif e > end:
            covered += e - end
            end = e
    return covered


def aligned_query_bases(cols):
    """Query bases that are aligned to the target (CIGAR =/X). A query span
    qe - qs also contains the 'I' runs, i.e. query sequence that is ABSENT
    from the target."""
    return sum(int(n) for n, op in _CG.findall(require_eqx(cols)) if op in "=X")


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
