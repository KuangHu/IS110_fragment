#!/usr/bin/env python3
"""Tests for the findings of the final review (docs/reviews/
2026-11-18_codex_final.md) and for the per-element PAF reader.

    python3 tests/test_review_final.py
"""
import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(HERE, "..", "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)
import build_records  # noqa: E402
import collect_confirmed_events as cce  # noqa: E402
import pipeline  # noqa: E402
import lib_alleles as la  # noqa: E402
from lib_alleles import Hit  # noqa: E402
from test_stage_units import _anchor_hits, _FakeFasta, rseq  # noqa: E402

OPTS = {"tol_bp": 100, "tol_frac": 0.002, "margin": 0.05, "max_extra_bp": 50000,
        "decompose_max_d": 5000, "min_is_size": 0, "max_is_size": 0}


def _paf_row(q, t, ts, te, strand="+", qlen=500):
    return "\t".join(map(str, [q, qlen, 0, qlen, strand, t, 900000, ts, te,
                               qlen, qlen, 60, "tp:A:P", "cg:Z:500="])) + "\n"


def _pair(q_ref, D, t, up_ts, gap, strand="+"):
    return [_paf_row(f"{q_ref}__up{D}", t, up_ts, up_ts + 500, strand),
            _paf_row(f"{q_ref}__down{D}", t, up_ts + 500 + gap, up_ts + 1000 + gap, strand)]


# ---------------------------------------------------------------- finding 4
def test_source_core_without_assembly_prefix():
    """'<contig>_<n>' comes from target '<contig>': same genome."""
    assert la.source_core("isolateA_7") == "isolateA"
    assert la.source_core("GCA_000001.1|CP1.1_7") == "000001"
    assert la.source_core("lab|scaffold_3_12") == "lab"
    hits = _anchor_hits("isolateA_7", 1000, ["isolateA", "isolateB"])
    pairs, _, _ = build_records.pair_anchors(hits, "isolateA_7", 1000, {"tnp_len": 1000}, OPTS)
    assert [p["assembly"] for p in pairs] == ["isolateB"], pairs


# ---------------------------------------------------------------- finding 6
def test_filled_observation_at_a_larger_d_is_used():
    """D=1000 only has a pair larger than the filled state, D=5000 has the
    exact filled one: the record is 'seen filled', not 'no observation'."""
    rid, t = "GCA_000001.1|CP1.1_7", "GCA_000002.1|X"
    hits = _anchor_hits(rid, 5000, [t], tnp_len=1200)
    hits[(rid, "up", 1000)] = {t: [Hit("u", 500, 0, 500, "+", t, 900000, 10000, 10500,
                                       500, 500, 100.0, 100.0)]}
    hits[(rid, "down", 1000)] = {t: [Hit("d", 500, 0, 500, "+", t, 900000, 20500, 21000,
                                         500, 500, 100.0, 100.0)]}       # gap 10,000
    best = build_records.pick_best_d(hits, rid, {"tnp_len": 1200}, [1000, 5000], OPTS)
    assert best["D"] == 5000 and len(best["v0_obs"]) == 1, (best["D"], best["v0_obs"])


# --------------------------------------------------------------- finding 11
def test_stub_keeps_ambiguous_placements():
    rid, t = "GCA_000001.1|CP1.1_7", "GCA_000002.1|X"
    h = lambda q, ts: Hit(q, 500, 0, 500, "+", t, 900000, ts, ts + 500, 500, 500, 100.0, 100.0)  # noqa: E731
    hits = {(rid, "up", 1000): {t: [h("u", 10000), h("u", 400000)]},
            (rid, "down", 1000): {t: [h("d", 13500), h("d", 403500)]}}
    _, _, n_amb = build_records.pair_anchors(hits, rid, 1000, {"tnp_len": 1000}, OPTS)
    assert n_amb == 1
    best = build_records.pick_best_d(hits, rid, {"tnp_len": 1000}, [1000], OPTS)
    assert best["pairs"] == [] and best["n_ambiguous"] == 1, best


# --------------------------------------------------------------- finding 12
def test_legacy_anchor_names_reach_stage_6():
    """'<ref>__up' / '__down' (no distance) are the anchors of a single-D
    reference: stages 5/6 must see their observation."""
    rid = "GCA_000001.1|CP1.1_7"
    named = _anchor_hits(rid, 1000, ["GCA_000002.1|X"], tnp_len=1200)
    legacy = {(rid, "up", None): named[(rid, "up", 1000)],
              (rid, "down", None): named[(rid, "down", 1000)]}
    info = {"tnp_len": 1200, "distances": [1000]}
    best = build_records.pick_best_d(legacy, rid, info, [1000], OPTS)
    assert len(best["v0_obs"]) == 1, best
    try:
        la.anchor_side(legacy, rid, "up", 1000, 2)
    except ValueError:
        pass
    else:
        raise AssertionError("legacy anchors of a multi-D reference must be refused")


# ----------------------------------------------------------- findings 5, 8
def _run_stage6(td, lines, ref_rows, hit_rows, seqs, extra=()):
    """Run build_records.py on a tiny synthetic input; return its records."""
    fa = os.path.join(td, "db.fa")
    with open(fa, "w") as f:
        for name, s in seqs.items():
            f.write(f">{name}\n{s}\n")
    subprocess.run(["samtools", "faidx", fa], check=True)
    paf = os.path.join(td, "a.paf")
    open(paf, "w").writelines(lines)
    tab = os.path.join(td, "table.tsv")
    with open(tab, "w") as f:
        f.write("is_id\tref_len\tanchor_len\ttnp_len\tflank_size\tdistances\n")
        f.writelines(ref_rows)
    hits = os.path.join(td, "is_hits.tsv")
    with open(hits, "w") as f:
        f.write("is_id\tassembly\tcontig\ttnp_start\ttnp_end\ttnp_strand\n")
        f.writelines(hit_rows)
    out = os.path.join(td, "rec")
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "build_records.py"),
                        "--paf", paf, "--ref-table", tab, "--hits", hits,
                        "--genome-db", fa, "--out", out, "--threads", "2",
                        "--max-secondary", "50", *extra],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    with open(os.path.join(out, "records.json")) as f:
        return {x["ref_id"]: x for x in json.load(f)}, r.stderr


def test_hit_without_anchors_gets_a_record():
    """A transposase on a contig too short for any anchor pair is in is_hits
    but not in the reference table. It used to vanish between stages."""
    with tempfile.TemporaryDirectory() as td:
        seqs = {"A|c1": rseq(20000), "A|short": rseq(5000), "B|c1": rseq(20000)}
        lines = _pair("A|c1_1", 1000, "B|c1", 5000, 3200)          # filled elsewhere
        recs, _ = _run_stage6(
            td, lines,
            ["A|c1_1\t4200\t500\t1200\t1500\t1000\n"],
            ["A|c1_1\tA\tA|c1\t8001\t9200\t+\n", "A|short_1\tA\tA|short\t101\t1300\t+\n"],
            seqs)
        assert recs["A|c1_1"]["boundary_evidence"]["method"] == "no_empty_tnp_only"
        short = recs["A|short_1"]
        assert short["boundary_evidence"]["method"] == "no_anchors_tnp_only", short["boundary_evidence"]
        assert short["transposase_cds"]["sequence"] == seqs["A|short"][100:1300]
        assert short["boundary_evidence"]["call_confidence"] == "no_call"
        recs, _ = _run_stage6(
            td, lines, ["A|c1_1\t4200\t500\t1200\t1500\t1000\n"],
            ["A|c1_1\tA\tA|c1\t8001\t9200\t+\n", "A|short_1\tA\tA|short\t101\t1300\t+\n"],
            seqs, extra=["--drop-unsupported"])
        assert set(recs) == {"A|c1_1"}, set(recs)


def test_drop_unsupported_applies_after_rejection():
    """The only 'empty' site is a 1,500-bp flank deletion next to an intact
    element. Decomposition rejects it; with --drop-unsupported the record
    must then go (it used to pass the filter before the rejection)."""
    with tempfile.TemporaryDirectory() as td:
        up, tnp, down = rseq(8000), rseq(1200), rseq(8000)
        src = up + tnp + down
        # target: same locus with 1,500 bp of the upstream flank deleted
        # (src 6300..7800), transposase kept. Anchors at D=5000: up =
        # src[2500:3000], down = src[14200:14700]; on the target the down
        # anchor sits 1,500 bp further left.
        tgt = up[:6300] + up[7800:] + tnp + down
        seqs = {"A|c1": src, "B|c1": tgt}
        lines = [_paf_row("A|c1_1__up5000", "B|c1", 2500, 3000),
                 _paf_row("A|c1_1__down5000", "B|c1", 14200 - 1500, 14700 - 1500)]
        args = (td, lines, ["A|c1_1\t12200\t500\t1200\t5500\t5000\n"],
                ["A|c1_1\tA\tA|c1\t8001\t9200\t+\n"], seqs)
        recs, _ = _run_stage6(*args)
        be = recs["A|c1_1"]["boundary_evidence"]
        assert be["method"] == "no_observation_tnp_only", be["method"]
        assert be["n_empty_rejected_flank_deletion"] == 1 and be["n_empty_observations"] == 0
        recs, _ = _run_stage6(*args, extra=["--drop-unsupported"])
        assert recs == {}, recs


# ---------------------------------------------------------------- finding 9
def test_cap_is_judged_per_index_part():
    """4 rows from two index parts of 2 rows each is not a capped search at
    -N 2 (cap = 3 rows in ONE part); 3 rows in one part is."""
    with tempfile.TemporaryDirectory() as td:
        paf = os.path.join(td, "raw.paf")
        part = lambda ts: [_paf_row("r__up1000", f"G{i}|c", t, t + 500) for i, t in ts]  # noqa: E731
        # part 1: queries up, down; part 2: up, down again
        rows = (part([(1, 100), (2, 100)]) + [_paf_row("r__down1000", "G1|c", 4000, 4500)]
                + part([(3, 100), (4, 100)]) + [_paf_row("r__down1000", "G3|c", 4000, 4500)])
        open(paf, "w").writelines(rows)
        # raw multi-part PAF: elements are not contiguous? (one element here,
        # so it is) -- run counter from the stream
        rc = {}
        la.load_anchor_hits(paf, 0, 0, rc)
        assert rc["r__up1000"] == 2 and rc["r__down1000"] == 1, rc
        # sorted by search_anchors.py: the runs are gone from the PAF and
        # come from the sidecar instead
        out = os.path.join(td, "sorted.paf")
        subprocess.run([sys.executable, os.path.join(SCRIPTS, "search_anchors.py"),
                        "--from-paf", paf, "--out", out], check=True, capture_output=True)
        side = open(out + ".rows.tsv").read().split("\n")
        assert side[0].split("\t") == ["anchor", "rows", "index_parts", "max_rows_one_part"]
        assert "r__up1000\t4\t2\t2" in side, side
        rc = {}
        groups = list(la.iter_anchor_groups(out, 0, 0, row_counts=rc))
        assert rc["r__up1000"] == 2 and len(groups) == 1
        assert len(groups[0][1][("r", "up", 1000)]) == 4          # all four targets kept
        assert all(len(l.rstrip("\n").split("\t")) == 12 for l in open(out))


def test_multipart_raw_paf_is_rejected_then_sorted():
    """Two elements, two index parts: element rows are interleaved in the raw
    PAF. The per-element reader must refuse it, and read the sorted one."""
    with tempfile.TemporaryDirectory() as td:
        raw = os.path.join(td, "raw.paf")
        rows = (_pair("A|c_1", 1000, "G1|c", 100, 3000) + _pair("B|c_2", 1000, "G1|c", 9000, 3000)
                + _pair("A|c_1", 1000, "G2|c", 100, 3000) + _pair("B|c_2", 1000, "G2|c", 9000, 3000))
        open(raw, "w").writelines(rows)
        try:
            list(la.iter_anchor_groups(raw, 0, 0))
        except ValueError as e:
            assert "not contiguous" in str(e)
        else:
            raise AssertionError("interleaved elements must be refused")
        out = os.path.join(td, "s.paf")
        subprocess.run([sys.executable, os.path.join(SCRIPTS, "search_anchors.py"),
                        "--from-paf", raw, "--out", out], check=True, capture_output=True)
        got = {r: {k: sorted(v) for k, v in h.items()}
               for r, h in la.iter_anchor_groups(out, 0, 0)}
        whole, _, _ = la.load_anchor_hits(raw, 0, 0)
        assert set(got) == {"A|c_1", "B|c_2"}
        for key, by_t in whole.items():
            assert sorted(by_t) == got[key[0]][key], key


# ---------------------------------------------------------------- finding 3
def test_stage_signature_follows_the_code():
    with tempfile.TemporaryDirectory() as td:
        script = os.path.join(td, "worker.py")
        out = os.path.join(td, "o.txt")
        open(script, "w").write(f"open({out!r}, 'a').write('1')\n")
        cmd = [sys.executable, script]
        run = lambda: pipeline._run_stage(td, "4_pairs", "t", cmd, [out])  # noqa: E731
        run(); run()
        assert open(out).read() == "1"                       # cached
        open(script, "w").write(f"open({out!r}, 'a').write('2')\n")
        run()
        assert open(out).read() == "12", "changed code must re-run the stage"


# ----------------------------------------------------------- findings 2, 10
def test_target_events_compare_anchor_placements():
    row = lambda eid, asm, hits, cat="inversion_only": {  # noqa: E731
        "event_id": eid, "category": cat, "tgt_asm": asm, "_hits": hits}
    # nested inversions: envelopes overlap, no anchor coincides -> 2 events
    flat = [row("e1", "GCA_1.1", [("c", 10000, 10500, "+"), ("c", 100000, 100500, "-")]),
            row("e2", "GCA_1.1", [("c", 40000, 40500, "+"), ("c", 70000, 70500, "-")])]
    cce.group_target_events(flat, "sp")
    assert flat[0]["target_event"] != flat[1]["target_event"], flat
    # the same placements seen from another reference copy -> 1 event,
    # also when that copy lies in the opposite orientation
    flat = [row("e1", "GCA_1.1", [("c", 10000, 10500, "+"), ("c", 100000, 100500, "-")]),
            row("e2", "GCA_1.1", [("c", 10020, 10520, "+"), ("c", 100010, 100510, "-")]),
            row("e3", "GCA_1.1", [("c", 100000, 100500, "+"), ("c", 10000, 10500, "-")])]
    cce.group_target_events(flat, "sp")
    assert len({r["target_event"] for r in flat}) == 1, flat
    assert flat[0]["n_records_same_target_event"] == 3
    # GCA/GCF twins of the target are one genome
    flat = [row("e1", "GCA_2.1", [("CP1.1", 100, 600, "+"), ("CP1.1", 9000, 9500, "-")]),
            row("e2", "GCF_2.1", [("NZ_CP1.1", 100, 600, "+"), ("NZ_CP1.1", 9000, 9500, "-")]),
            row("e3", "GCF_3.1", [("NZ_CP1.1", 100, 600, "+"), ("NZ_CP1.1", 9000, 9500, "-")])]
    cce.group_target_events(flat, "sp")
    assert flat[0]["target_event"] == flat[1]["target_event"] != flat[2]["target_event"]
    # a chain of near-misses does not merge its ends
    flat = [row(f"e{i}", "GCA_1.1", [("c", 10000 + 900 * i, 10500 + 900 * i, "+")])
            for i in range(4)]
    cce.group_target_events(flat, "sp")
    assert flat[0]["target_event"] != flat[3]["target_event"], flat


# ---------------------------------------------------------------- finding 7
def test_tandem_filter_knows_every_method_label():
    import verify_tandem_enrichment as vte
    with tempfile.TemporaryDirectory() as td:
        for method, want in (("no_observation_tnp_only", 0), ("no_anchors_tnp_only", 0),
                             ("some_future_label", 0), ("empty_vs_filled_decomposition", 1)):
            p = os.path.join(td, "r.json")
            json.dump([{"ref_id": "r", "is_element": {"sequence": rseq(2000), "length": 2000},
                        "source": {"contig": "c"},
                        "boundary_evidence": {"method": method}}], open(p, "w"))
            seqs, _, stats = vte.load_is_sequences(p, 1500, called_only=True, dedup=False)
            assert len(seqs) == want, (method, len(seqs), dict(stats))


# --------------------------------------------------------------- finding 13
def test_review_case_locus_is_zero_based():
    import prepare_review_cases as prc
    with tempfile.TemporaryDirectory() as td:
        p = os.path.join(td, "r.json")
        json.dump([{"ref_id": "r", "source": {"contig": "A|c", "transposase_start": 101,
                                               "transposase_end": 1300},
                    "is_element": {"source_start": 101, "source_end": 1300}}], open(p, "w"))
        contig, s0, e0 = prc.load_locus(p, "r")
        assert (contig, s0, e0) == ("c", 100, 1300) and e0 - s0 == 1200


# ---------------------------------------------------------------- finding 1
def test_reference_step_survives_without_clonal_observation():
    """Source element + one other genome carrying it with a 250-bp insertion:
    the lineage is reference -> variant, although no OTHER genome carries the
    reference state (the source genome itself is not an observation)."""
    src = open(os.path.join(SCRIPTS, "build_lineages.py")).read()
    a = src.index("        # The reference element is a KNOWN state of this locus")
    b = src.index("        # Cluster at 99% identity")
    block = src[a:b]
    import textwrap
    ns = {"variants_raw": {"v1": {"vid": "v1"}}, "n_empty": 0, "n_clonal": 0,
          "v_ref_seqs": {"R": "ACGT"}, "v1_id": "R", "v1_len": 6000}
    code = "for _ in (0,):\n" + textwrap.indent(textwrap.dedent(block), "    ")
    exec(code, ns)
    v = ns["variants_raw"]["V_REF"]
    assert v["n_observations"] == 0 and v["observed_in_other_genome"] is False
    assert v["len"] == 6000 and v["seq"] == "ACGT"
    # nothing at all for this element -> still skipped
    ns = {"variants_raw": {}, "n_empty": 0, "n_clonal": 0,
          "v_ref_seqs": {}, "v1_id": "R", "v1_len": 6000}
    exec(code, ns)
    assert "V_REF" not in ns["variants_raw"]


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"PASS {name}")
    print("ALL PASS")
