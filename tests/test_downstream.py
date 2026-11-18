#!/usr/bin/env python3
"""Unit tests for stages 7-11 (variants, lineages, V0 check, tandem).

Needs minimap2, samtools, cd-hit-est on PATH:  python3 tests/test_downstream.py
"""
import os
import random
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
import build_lineages  # noqa: E402
import detect_tandem  # noqa: E402
import find_variants  # noqa: E402

rng = random.Random(11)
r = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731


def tandem(seq):
    with tempfile.TemporaryDirectory() as wd:
        hits = detect_tandem.self_minimap2(seq, wd, 2)
    return detect_tandem.find_tandem_unit(hits, 50, 95)


def test_tandem_copy_number_is_exact():
    u = r(300)
    for n in (2, 3, 7):
        res = tandem(r(1500) + u * n + r(1500))
        assert res and res[0]["repeat_class"] == "tandem", res
        assert res[0]["n_copies"] == n, (n, res)            # was n + 1
        assert abs(res[0]["repeat_unit_bp"] - 300) <= 3, res  # was rounded to 50
    u = r(730)
    res = tandem(r(800) + u * 4 + r(800))
    assert res[0]["n_copies"] == 4 and abs(res[0]["repeat_unit_bp"] - 730) <= 3, res


def test_dispersed_repeat_is_not_tandem():
    w = r(600)
    res = tandem(r(500) + w + r(800) + w + r(500))
    assert res and all(x["repeat_class"] == "dispersed" for x in res), res
    assert abs(res[0]["copy_spacing_bp"] - 1400) <= 3
    assert tandem(r(4000)) == []


def test_union_len():
    assert find_variants.union_len([(0, 100), (50, 150), (300, 400)]) == 250
    assert find_variants.union_len([]) == 0


def test_own_v1_alignment_is_not_outcompeted():
    # 30 near-identical elements; each between-seq must get blocks to ITS V1
    base = r(1500)
    with tempfile.TemporaryDirectory() as wd:
        for i in range(3):
            v1 = base[:100 + i] + "A" + base[101 + i:]
            rid, res = find_variants.align_to_own_v1((f"v{i}", v1, [("q", v1)], wd))
            assert res and find_variants.union_len(
                [(b["v1_start"], b["v1_end"]) for b in res["q"]]) >= 1450, res


def test_fetch_batch_duplicates_and_failure():
    with tempfile.TemporaryDirectory() as wd:
        fa = os.path.join(wd, "g.fa")
        seq = r(5000)
        with open(fa, "w") as f:
            f.write(">asm|c1\n" + seq + "\n")
        subprocess.run(["samtools", "faidx", fa], check=True)
        reg = ("asm|c1", 101, 200)
        got = find_variants.extract_target_batch(fa, [reg, reg, ("asm|c1", 1, 50)])
        assert got[reg] == seq[100:200] and len(got) == 2
        try:
            find_variants.extract_target_batch(fa, [("nope", 1, 10)])
        except RuntimeError:
            pass
        else:
            raise AssertionError("missing contig must raise, not be dropped")


def test_build_chains_branches_from_inner_step():
    # A < B, and C contains A but not B: A < C must be kept
    vl = [{"vid": "A"}, {"vid": "B"}, {"vid": "C"}]
    ok = {"cov_of_q": 100, "ident": 100}
    chains = build_lineages.build_chains(vl, {("A", "B"): ok, ("A", "C"): ok}, 95, 95)
    assert sorted(chains) == [[0, 1], [0, 2]], chains
    chains = build_lineages.build_chains(vl, {("A", "B"): ok, ("B", "C"): ok,
                                              ("A", "C"): ok}, 95, 95)
    assert chains == [[0, 1, 2]], chains


def test_nested_variants_are_not_clustered_together():
    core = r(1200)
    longer = core[:600] + r(1500) + core[600:]      # contains `core` entirely
    twin = core[:300] + "T" + core[301:]            # same length, 1 SNP
    with tempfile.TemporaryDirectory() as wd:
        m = build_lineages.cluster_variants(
            [("core", core), ("longer", longer), ("twin", twin)], 0.99, wd)
    assert m["core"] != m["longer"], m               # used to merge
    assert m["core"] == m["twin"], m


def test_pairwise_coverage_is_a_union():
    a = r(1000)
    b = a + r(500) + a                               # a occurs twice in b
    with tempfile.TemporaryDirectory() as wd:
        pm = build_lineages.pairwise_minimap2([("a", a), ("b", b)], wd, 2)
    assert 95 <= pm[("a", "b")]["cov_of_q"] <= 100.01, pm


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"PASS {t.__name__}")
    print(f"{len(tests)} tests passed")
