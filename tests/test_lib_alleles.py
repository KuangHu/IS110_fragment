#!/usr/bin/env python3
"""Unit tests for scripts/lib_alleles.py. Plain asserts; run directly:

    python3 tests/test_lib_alleles.py
"""
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
from lib_alleles import (Hit, assembly_core, canonical_insert_key,  # noqa: E402
                         classify_gap, decompose_alleles, find_tsd, gap_tolerance,
                         place_pair, representative_asm, revcomp)

rng = random.Random(7)


def rseq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(s, rate):
    out = list(s)
    for i in range(len(out)):
        if rng.random() < rate:
            out[i] = rng.choice([b for b in "ACGT" if b != out[i]])
    return "".join(out)


def hit(tname, strand, ts, te, qlen=500, nmatch=500):
    return Hit("q", qlen, 0, qlen, strand, tname, 10 ** 6, ts, te, nmatch, qlen,
               nmatch / qlen * 100, 100.0)


def test_assemblies():
    assert assembly_core("GCA_000005845.2") == assembly_core("GCF_000005845.1")
    assert assembly_core("my_genome") == "my_genome"
    reps = representative_asm(["GCA_000005845.2", "GCF_000005845.1", "GCA_000009999.1"])
    assert reps == {"000005845": "GCF_000005845.1", "000009999": "GCA_000009999.1"}
    # order-independent
    assert reps == representative_asm(["GCA_000009999.1", "GCF_000005845.1",
                                       "GCA_000005845.2"])


def test_place_pair():
    # clean empty site: anchors abut (gap 0) -- used to be discarded
    pl = place_pair([hit("a|c", "+", 0, 500)], [hit("a|c", "+", 500, 1000)])
    assert pl and pl["gap"] == 0 and pl["status"] == "unique"
    # TSD overlap: gap -5 is allowed
    pl = place_pair([hit("a|c", "+", 0, 500)], [hit("a|c", "+", 495, 995)])
    assert pl and pl["gap"] == -5
    # minus strand swaps sides
    pl = place_pair([hit("a|c", "-", 2000, 2500)], [hit("a|c", "-", 0, 500)])
    assert pl and pl["gap"] == 1500 and pl["start"] == 500 and pl["end"] == 2000
    # two equally good down copies -> ambiguous, not forced
    pl = place_pair([hit("a|c", "+", 0, 500)],
                    [hit("a|c", "+", 2000, 2500), hit("a|c", "+", 9000, 9500)])
    assert pl["status"] == "ambiguous" and pl["n_pairs"] == 2
    # a clearly better pair wins
    pl = place_pair([hit("a|c", "+", 0, 500)],
                    [hit("a|c", "+", 2000, 2500), hit("a|c", "+", 9000, 9500, nmatch=400)])
    assert pl["status"] == "resolved" and pl["end"] == 2000
    # beyond max_gap -> nothing
    assert place_pair([hit("a|c", "+", 0, 500)], [hit("a|c", "+", 90000, 90500)],
                      max_gap=50000) is None


def qhit(tname, strand, ts, te, qs, qe, qlen=500):
    n = qe - qs
    return Hit("q", qlen, qs, qe, strand, tname, 10 ** 6, ts, te, n, n, 100.0, n / qlen * 100)


def test_place_pair_projects_clipped_anchors():
    """Anchors pass at 80% coverage; the gap is measured between the anchors'
    true inner edges, not between the alignment ends."""
    # + strand: up anchor at 1000-1500 aligned only on query [0,400) -> te=1400;
    # down anchor at 4700-5200 aligned on query [100,500) -> ts=4800. True gap 3200.
    pl = place_pair([qhit("a|c", "+", 1000, 1400, 0, 400)],
                    [qhit("a|c", "+", 4800, 5200, 100, 500)])
    assert pl["gap"] == 3200 and (pl["start"], pl["end"]) == (1500, 4700), pl
    # - strand: down anchor occupies 1000-1500 (reversed), up anchor 4700-5200.
    # down aligned on query [100,500) -> target [1000,1400); up aligned on
    # query [0,400) -> target [4800,5200).
    pl = place_pair([qhit("a|c", "-", 4800, 5200, 0, 400)],
                    [qhit("a|c", "-", 1000, 1400, 100, 500)])
    assert pl["gap"] == 3200 and (pl["start"], pl["end"]) == (1500, 4700), pl


def test_place_pair_no_cross_copy():
    """Two loci; at each, one anchor is full length and the other clipped. The
    best-scoring combination is the CROSS pair (full up of locus 1 + full down
    of locus 2) -- it must not be returned."""
    ups = [qhit("a|c", "+", 1000, 1500, 0, 500),          # locus 1, full
           qhit("a|c", "+", 10000, 10400, 0, 400)]        # locus 2, clipped outer... inner
    downs = [qhit("a|c", "+", 4800, 5200, 100, 500),      # locus 1, clipped
             qhit("a|c", "+", 13700, 14200, 0, 500)]      # locus 2, full
    pl = place_pair(ups, downs, max_gap=50000)
    assert pl["gap"] == 3200, pl                          # a same-locus pair
    assert pl["status"] == "ambiguous"                    # two loci, same score
    assert pl["n_pairs"] == 2


def test_place_per_assembly_twin_fallback():
    from lib_alleles import place_per_assembly
    up = {"GCF_000000009.1|x1": [hit("GCF_000000009.1|x1", "+", 0, 500)],
          "GCA_000000009.1|y": [hit("GCA_000000009.1|y", "+", 0, 500)]}
    down = {"GCF_000000009.1|x2": [hit("GCF_000000009.1|x2", "+", 900, 1400)],
            "GCA_000000009.1|y": [hit("GCA_000000009.1|y", "+", 900, 1400)]}
    out = place_per_assembly(up, down)
    assert list(out) == ["000000009"] and out["000000009"]["assembly"] == "GCA_000000009.1"


def test_tolerant_terminal_gap_and_json():
    import json
    from lib_alleles import tolerant_decompose
    L, R, ins, extra = rseq(900), rseq(900), rseq(1500), rseq(1200)
    # a second 1,200-bp extra block at the START of the long allele: not one
    # dominant insertion -> must be refused
    td = tolerant_decompose(L + R, extra + L + ins + R, 0.85, 0.95, 50, 3.0)
    assert td is None or not td["ok"], td
    # SNP-divergent empty -> alignment path; result must be JSON-serialisable
    r = decompose_alleles(mutate(L, 0.01) + mutate(R, 0.01), L + ins + R)
    assert r["method"] == "tolerant_alignment" and r["insert_len"] == 1500
    json.dumps(r)


def test_classify_gap_min_max():
    tol = gap_tolerance(3200)
    assert classify_gap(2000, 3200, 1200, tol) == "empty"
    assert classify_gap(2000, 3200, 1200, tol, min_removed=1500) == "intermediate"
    assert classify_gap(2000, 3200, 1200, tol, max_removed=1000) == "very_short"


def test_classify_gap_large_D():
    # D = 40 kb, tnp 1200, filled source gap = 81,200. A target with an
    # intact 1.45 kb IS110 must NOT read as empty (old tolerance: +/-8 kb).
    D, tnp = 40000, 1200
    exp_v0 = 2 * D + tnp
    tol = gap_tolerance(exp_v0)
    assert classify_gap(exp_v0 + 30, exp_v0, tnp, tol) == "v0_filled"
    assert classify_gap(exp_v0 - 1450, exp_v0, tnp, tol) == "empty"
    assert classify_gap(exp_v0 - 600, exp_v0, tnp, tol) == "intermediate"
    assert classify_gap(exp_v0 + 5000, exp_v0, tnp, tol) == "v1plus_filled"
    # the old rule: |gap - 2D| <= max(200, 0.1*2D) called this filled site empty
    assert abs(exp_v0 - 2 * D) <= max(200, 0.1 * 2 * D)


def test_decompose_exact_and_offset():
    left, right, ins = rseq(800), rseq(800), rseq(1450)
    short, long_ = left + right, left + ins + right
    r = decompose_alleles(short, long_)
    mh = r["junction_microhomology_bp"]
    assert r["method"] == "exact"
    assert r["insert_len"] == 1450
    assert 800 <= r["insert_start"] <= 800 + mh
    # reconstruction: short[:offset] + insert + short[offset:] == long
    o = r["offset"]
    assert short[:o] + r["insert_seq"] + short[o:] == long_


def test_decompose_tolerant_with_snps():
    left, right, ins = rseq(900), rseq(900), rseq(1300)
    long_ = left + ins + right
    short = mutate(left, 0.01) + mutate(right, 0.01)
    r = decompose_alleles(short, long_)
    assert r["method"] == "tolerant_alignment", r
    assert abs(r["insert_len"] - 1300) <= 5
    assert abs(r["insert_start"] - 900) <= 5 + r["junction_microhomology_bp"]
    # empty interval of zero length -> pure insertion
    r0 = decompose_alleles("", ins)
    assert r0["event_class"] == "pure_insertion" and r0["insert_len"] == 1300


def test_decompose_refuses_substitution():
    a = rseq(1000)
    assert decompose_alleles(a, mutate(a, 0.02))["event_class"] == "substitution_only"


def test_canonical_insert_key_invariance():
    # insert flanked by a 6 bp direct repeat: several equally valid boundaries
    left, right, core = rseq(300), rseq(300), rseq(1000)
    dr = "GATTAC"
    long_ = left + dr + core + dr + right
    ilen = len(dr) + len(core)
    k_left = canonical_insert_key(long_, 300, ilen)
    k_right = canonical_insert_key(long_, 306, ilen)
    rc = revcomp(long_)
    k_rc = canonical_insert_key(rc, len(rc) - (300 + ilen), ilen)
    assert k_left == k_right == k_rc


def test_find_tsd():
    left, right, core = rseq(50), rseq(50), rseq(500)
    tsd = "ACGTTG"
    n, seq, side, _ = find_tsd(core + tsd, left + tsd, right)
    assert n >= 6 and seq.endswith(tsd)


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"PASS {t.__name__}")
    print(f"{len(tests)} tests passed")
