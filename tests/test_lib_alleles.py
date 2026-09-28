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
