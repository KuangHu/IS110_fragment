#!/usr/bin/env python3
"""Checks for the 2026-11 review fixes (synthetic PAF rows / random DNA only).
Run:  python3 tests/test_fragment.py"""
import os, random, subprocess, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
SCR = os.path.join(HERE, "..", "scripts")
sys.path.insert(0, SCR)
import paf_utils as pu

fails = []
def check(name, cond, info=""):
    print(("PASS " if cond else "FAIL ") + name + (f"  [{info}]" if info and not cond else ""))
    if not cond: fails.append(name)

def row(q, qlen, qs, qe, strand, t, tlen, ts, te, cg=None, matches=None, extra=()):
    block = max(qe - qs, te - ts)
    m = matches if matches is not None else min(qe - qs, te - ts)
    c = [q, qlen, qs, qe, strand, t, tlen, ts, te, m, block, 60]
    if cg: c.append("cg:Z:" + cg)
    return "\t".join(map(str, list(c) + list(extra)))

def detect(rows, *opts):
    with tempfile.TemporaryDirectory() as d:
        paf = os.path.join(d, "x.paf")
        open(paf, "w").write("\n".join(rows) + "\n")
        r = subprocess.run([sys.executable, os.path.join(SCR, "detect_is110_variants.py"),
                            paf, os.path.join(d, "o"), *opts], capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        out = [l.split("\t") for l in open(os.path.join(d, "o_all.tsv")).read().splitlines()[1:]]
        return out

def test_identity():
    c = row("q", 2000, 0, 2000, "+", "t", 9000, 0, 2600, "1000=600D1000=", 2000).split("\t")
    check("identity ignores the gap", abs(pu.aligned_identity(c) - 100) < 1e-9)
    c = "q 2000 0 2000 + t 9000 0 2600 2000 2000 60 cg:Z:1000=600D1000= nn:i:600".split()
    check("identity never above 100", pu.aligned_identity(c) == 100.0)
    c = row("q", 2400, 0, 2400, "+", "t", 9000, 0, 1200, "400=600I400=600I400=", 1200).split("\t")
    check("aligned query bases exclude I", pu.aligned_query_bases(c) == 1200)
    check("blocks of a multi-indel row",
          pu.aligned_blocks(c, 500, 200) == [(0, 400, 0, 400), (1000, 1400, 400, 800), (2000, 2400, 800, 1200)])
    c = row("q", 2000, 0, 2000, "-", "t", 9000, 100, 2700, "1200=600D800=", 2000).split("\t")
    check("blocks on the minus strand",
          pu.aligned_blocks(c, 500, 200) == [(800, 2000, 100, 1300), (0, 800, 1900, 2700)],
          str(pu.aligned_blocks(c, 500, 200)))

def test_eqx_is_required():
    for cg in ("1000M", None):
        c = row("q", 1000, 0, 1000, "+", "t", 9000, 0, 1000, cg).split("\t")
        for fn in (pu.aligned_identity, pu.aligned_query_bases, lambda x: pu.aligned_blocks(x, 500, 200)):
            try:
                fn(c); ok = False
            except ValueError:
                ok = True
            check(f"row with CIGAR {cg!r} is refused (no approximate identity)", ok)
    c = row("q", 1000, 0, 1000, "+", "t", 9000, 0, 1000, "950=50X").split("\t")
    check("identity is '=' / ('=' + 'X')", abs(pu.aligned_identity(c) - 95.0) < 1e-9)

def test_output_independent_of_hash_seed():
    rows = []
    for k in range(5):
        for q in ("qa", "qb", "qc"):
            rows += [row(q, 2000, 0, 1000 + k, "+", f"t{k}", 9000, 0, 1000 + k, f"{1000 + k}="),
                     row(q, 2000, 1000 + k, 2000, "+", f"t{k}", 9000, 1600 + k, 2600, f"{1000 - k}=")]
    outs = []
    with tempfile.TemporaryDirectory() as d:
        open(os.path.join(d, "x.paf"), "w").write("\n".join(rows) + "\n")
        for seed in ("1", "2"):
            for script, suffix in (("detect_is110_variants.py", "_all.tsv"), ("detect_multisplit.py", "_chains.tsv")):
                subprocess.run([sys.executable, os.path.join(SCR, script), os.path.join(d, "x.paf"),
                                os.path.join(d, "o" + seed)], env=dict(os.environ, PYTHONHASHSEED=seed),
                               capture_output=True, check=True)
            outs.append(open(os.path.join(d, f"o{seed}_all.tsv")).read() + open(os.path.join(d, f"o{seed}_chains.tsv")).read())
    check("detector and multisplit tables identical under PYTHONHASHSEED=1 and 2",
          outs[0] == outs[1] and outs[0].count("insertion") >= 5, str(len(outs[0])))

def test_detector():
    two = [row("q", 2000, 0, 1000, "+", "t", 9000, 0, 1000, "1000="),
           row("q", 2000, 1000, 2000, "+", "t", 9000, 1600, 2600, "1000="),
           row("q", 2000, 0, 1000, "+", "t", 9000, 4000, 5000, "1000="),
           row("q", 2000, 1000, 2000, "+", "t", 9000, 5600, 6600, "1000=")]
    out = detect(two)
    check("two interrupted copies -> two insertions", len(out) == 2 and
          sorted((int(r[8]), int(r[13])) for r in out) == [(1000, 1600), (5000, 5600)], str(out))
    out = detect([row("q", 2400, 0, 2400, "+", "t", 9000, 0, 1200, "400=600I400=600I400=", 1200)],
                 "--min-query-cov", "60")
    check("multi-deletion row: 50% aligned does not pass 60%", out == [], str(out))
    out = detect([row("q", 2400, 0, 2400, "+", "t", 9000, 0, 1200, "400=600I400=600I400=", 1200)],
                 "--min-query-cov", "40")
    check("multi-deletion row: two deletions, 50% coverage each",
          len(out) == 2 and all(r[0] == "deletion" and r[19] == "50.0" for r in out), str(out))
    # overlapping junctions (as minimap2 reports them)
    out = detect([row("q", 2000, 0, 1000, "+", "t", 9000, 0, 1000, "1000="),
                  row("q", 2000, 999, 2000, "+", "t", 9000, 1499, 2500, "1001=")])
    check("500-bp insertion with a 1-bp junction overlap", len(out) == 1 and out[0][0] == "insertion", str(out))
    out = detect([row("q", 2000, 0, 900, "+", "t", 9000, 0, 900, "900="),
                  row("q", 2000, 1097, 2000, "+", "t", 9000, 897, 1800, "903=")])
    check("200-bp deletion with a 3-bp target overlap", len(out) == 1 and out[0][0] == "deletion", str(out))
    out = detect([row("q", 2000, 0, 1060, "+", "t", 9000, 0, 1060, "1060="),
                  row("q", 2000, 1000, 2000, "+", "t", 9000, 1600, 2600, "1000=")])
    check("600-bp insertion with 60-bp junction homology", len(out) == 1 and out[0][0] == "insertion", str(out))
    # non-events
    out = detect([row("q", 2000, 0, 1000, "+", "t", 9000, 0, 1000, "1000="),
                  row("q", 2000, 1000, 2000, "-", "t", 9000, 1600, 2600, "1000=")])
    check("inverted second half is not an insertion", out == [])
    out = detect([row("q", 2000, 1000, 2000, "+", "t", 9000, 0, 1000, "1000="),
                  row("q", 2000, 0, 1000, "+", "t", 9000, 1600, 2600, "1000=")])
    check("swapped halves are not an insertion", out == [])
    out = detect([row("q", 2000, 1000, 2000, "-", "t", 9000, 0, 1000, "1000="),
                  row("q", 2000, 0, 1000, "-", "t", 9000, 1600, 2600, "1000=")])
    check("minus-strand insertion", len(out) == 1 and out[0][0] == "insertion", str(out))
    # one locus seen by three references with junction offsets 0 / 20 / 21
    rows = []
    for k, off in enumerate((0, 20, 21)):
        rows += [row(f"q{k}", 2000, 0, 1000 + off, "+", "t", 9000, 0, 1000 + off, f"{1000 + off}="),
                 row(f"q{k}", 2000, 1000 + off, 2000, "+", "t", 9000, 1600 + off, 2600, f"{1000 - off}=")]
    out = detect(rows, "--max-junction-overlap", "100")
    import itertools
    outs = {tuple(map(tuple, detect([x for k in perm for x in rows[2 * k:2 * k + 2]])))
            for perm in itertools.permutations(range(3))}
    check("junction offsets 0/20/21: 0 and 20 are one locus, 21 is outside the 20-bp bound, "
          "whatever the input order", len(out) == 2 and len(outs) == 1, f"{len(out)} {len(outs)}")

def mm2_detect(query, contigs, script="detect_is110_variants.py", suffix="_all.tsv", opts=()):
    """real minimap2 with the pipeline options, both presets -> {preset: rows}"""
    res = {}
    with tempfile.TemporaryDirectory() as d:
        q, t, paf = (os.path.join(d, n) for n in ("q.fa", "t.fa", "x.paf"))
        open(q, "w").write(f">IS\n{query}\n")
        open(t, "w").write("".join(f">{k}\n{v}\n" for k, v in contigs.items()))
        for preset in ("asm10", "asm5"):
            r = subprocess.run(["minimap2", "-x", preset, "-c", "--eqx", "--secondary=yes", "-N", "1000000",
                                "-p", "0.05", "-g", "200", "-r", "200,200", t, q, "-o", paf],
                               capture_output=True, text=True)
            assert r.returncode == 0, r.stderr
            subprocess.run([sys.executable, os.path.join(SCR, script), paf, os.path.join(d, "o"), *opts],
                           capture_output=True, text=True, check=True)
            res[preset] = [l.split("\t") for l in open(os.path.join(d, "o" + suffix)).read().splitlines()[1:]]
    return res

def test_final_review_detector():
    rnd = random.Random(23)
    dna = lambda n: "".join(rnd.choice("ACGT") for _ in range(n))
    rc = lambda x: x.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    # 1. reference with an internal repeat, two INTACT copies on the contig
    R = dna(500); IS = R + R + dna(800)
    for name, tgt in (("plus", dna(300) + IS + dna(1000) + IS + dna(300)),
                      ("minus", dna(300) + rc(IS) + dna(1000) + rc(IS) + dna(300))):
        got = mm2_detect(IS, {"t": tgt})
        check(f"repeat pieces inside two intact copies ({name}) are no event",
              got["asm10"] == [] and got["asm5"] == [], str(got))
    # 2. a short competing match must not take the arm of the real junction
    A = row("q", 2000, 0, 500, "+", "t", 9000, 0, 500, "500=")
    B = row("q", 2000, 500, 700, "+", "t", 9000, 1100, 1300, "200=")
    C = row("q", 2000, 500, 2000, "+", "t", 9000, 1800, 3300, "1500=")
    for order in ((A, B, C), (C, B, A), (B, A, C)):
        out = detect(list(order))
        check("competing 200-bp match: the full-coverage insertion is kept",
              len(out) == 1 and out[0][0] == "insertion" and (out[0][8], out[0][13], out[0][19]) == ("500", "1800", "100.0"),
              str(out))
    IS2 = dna(2000)
    got = mm2_detect(IS2, {"p": dna(500) + IS2[:500] + dna(600) + IS2[500:700] + dna(500) + IS2[500:] + dna(500),
                           "m": dna(500) + rc(IS2[:500] + dna(600) + IS2[500:700] + dna(500) + IS2[500:]) + dna(500)})
    for preset, rows in got.items():
        check(f"minimap2 {preset}: insert carrying a 200-bp piece of the IS, one insertion per strand",
              sorted((r[3], r[0], r[19]) for r in rows) == [("m", "insertion", "100.0"), ("p", "insertion", "100.0")], str(rows))
    # 3. small I runs are not aligned sequence
    r1 = row("q", 3000, 0, 1400, "+", "t", 9000, 0, 1100, "400=100I300=100I300=100I100=", 1100)
    r2 = row("q", 3000, 1400, 1600, "+", "t", 9000, 1700, 1900, "200=")
    out = detect([r1, r2])
    check("three 100-bp I runs: 43% aligned does not pass the 50% filter", out == [], str(out))
    out = detect([r1, r2], "--min-query-cov", "40")
    check("... and is reported as 43.3% when the filter allows it", len(out) == 1 and out[0][19] == "43.3", str(out))
    out = multisplit([row("q", 3000, 0, 3000, "+", "t", 9000, 0, 2700, "698=100I700=100I700=100I602=", 2700)])
    check("full-span row with 300 unaligned bases is partial (90%), not intact",
          [(r["status"], r["query_coverage"]) for r in out] == [("partial", "90.0")], str(out))
    import parse_paf as pp
    with tempfile.TemporaryDirectory() as d:
        open(os.path.join(d, "x.paf"), "w").write(row("q", 2400, 0, 2400, "+", "t", 9000, 0, 1200, "400=600I400=600I400=", 1200) + "\n")
        g = pp.summarize_groups(pp.group_by_genome(pp.parse_paf(os.path.join(d, "x.paf"))))
    check("parse_paf summary counts aligned bases", (g[0]["total_aligned_bp"], g[0]["query_coverage_pct"]) == (1200, 50.0), str(g))
    # 4. two junctions of one chain stay two events; no drift across references
    out = detect([row("q", 1410, 0, 1410, "+", "t", 9000, 0, 1010, "500=200I10=200I500=", 1010)])
    check("deletions at 500 and 510 of one copy are two events",
          sorted(r[8] for r in out) == ["500", "510"] and all(r[0] == "deletion" for r in out), str(out))
    rows = []
    for k, off in enumerate((0, 20, 40)):
        rows += [row(f"q{k}", 2000, 0, 1000, "+", "t", 90000, 0, 1000 + off, f"{1000 + off}="),
                 row(f"q{k}", 2000, 1000, 2000, "+", "t", 90000, 1600 + off, 2600 + off, "1000=")]
    out = detect(rows)
    check("junctions 0 and 40 bp are not merged through one at 20", len(out) == 2, str(len(out)))
    # 9. net insertion of exactly the maximum with unaligned IS bases at the junction
    out = detect([row("q", 2000, 0, 1002, "+", "t", 90000, 0, 1002, "1002="),
                  row("q", 2000, 1030, 2000, "+", "t", 90000, 51030, 52000, "970=")])
    check("net insertion of 50,000 bp (target gap 50,028) is accepted",
          len(out) == 1 and out[0][0] == "insertion" and out[0][17] == "50028", str(out))
    out = detect([row("q", 2000, 0, 1002, "+", "t", 90000, 0, 1002, "1002="),
                  row("q", 2000, 1030, 2000, "+", "t", 90000, 51031, 52001, "970=")])
    check("net insertion of 50,001 bp is not", out == [], str(out))
    # 13. --max-query-overlap of detect_multisplit
    arms = [row("q", 1100, 0, 1000, "+", "t", 9000, 0, 1000, "1000="),
            row("q", 1100, 900, 1100, "+", "t", 9000, 1600, 1800, "200=")]
    st = {o: sorted(r["status"] for r in multisplit(arms, "--max-query-overlap", o)) for o in ("0", "0.3", "0.5")}
    check("multisplit enforces --max-query-overlap",
          st == {"0": ["partial", "partial"], "0.3": ["partial", "partial"], "0.5": ["variant"]}, str(st))
    # 17. N operations
    c = row("q", 2000, 0, 2000, "+", "t", 9000, 0, 2600, "1000=600N1000=", 2000).split("\t")
    try:
        pu.aligned_blocks(c, 500, 200); ok = False
    except ValueError:
        ok = True
    check("a CIGAR with an N operation is refused", ok)

def multisplit(rows, *opts):
    with tempfile.TemporaryDirectory() as d:
        paf = os.path.join(d, "x.paf")
        open(paf, "w").write("\n".join(rows) + "\n")
        r = subprocess.run([sys.executable, os.path.join(SCR, "detect_multisplit.py"),
                            paf, os.path.join(d, "o"), *opts], capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        import csv
        return list(csv.DictReader(open(os.path.join(d, "o_chains.tsv")), delimiter="\t"))

def test_multisplit():
    out = multisplit([row("q", 2000, 0, 2000, "+", "t", 90000, 0, 2000, "2000="),
                      row("q", 2000, 0, 2000, "+", "t", 90000, 30000, 32000, "2000="),
                      row("q", 2000, 0, 2000, "-", "t", 90000, 60000, 62000, "2000=")])
    check("three intact copies (two strands) are three rows",
          len(out) == 3 and all(r["status"] == "intact" for r in out), str(len(out)))
    out = multisplit([row("q", 2000, 0, 1000, "+", "t", 90000, 0, 1000, "1000=")])
    check("a half-length hit is partial, not intact", [r["status"] for r in out] == ["partial"])
    out = multisplit([row("q", 2000, 1000, 2000, "+", "t", 90000, 0, 1000, "1000="),
                      row("q", 2000, 0, 1000, "+", "t", 90000, 40000, 41000, "1000=")])
    check("halves of two loci are not one copy", sorted(r["status"] for r in out) == ["partial", "partial"])
    out = multisplit([row("q", 2400, 0, 2400, "+", "t", 90000, 0, 1200, "400=600I400=600I400=", 1200)],
                     "--min-query-cov", "40")
    check("two deletions inside one row are two splits",
          len(out) == 1 and out[0]["n_splits"] == "2" and out[0]["query_coverage"] == "50.0", str(out))

def test_lineage():
    import find_lineage as fl
    A, B, C = (0, 500), (1, 900), (2, 1300)
    check("nested observation finds its parent", fl.nesting_parent([A, (1, 910)], [[C], [A]]) == 1)
    check("unrelated insertions give no edge", fl.nesting_parent([B, C], [[A]]) is None)
    check("same cluster elsewhere in the IS gives no edge", fl.nesting_parent([(0, 1500), B], [[A]]) is None)
    rnd = random.Random(5)
    big = "".join(rnd.choice("ACGT") for _ in range(2500))
    small = big[700:1700]
    other = "".join(rnd.choice("ACGT") for _ in range(1500))
    c1, _ = fl.cluster_sequences([small, big, other], 95, 1)
    c2, _ = fl.cluster_sequences([big, small, other], 95, 1)
    check("clustering does not depend on input order",
          c1[0] == c1[1] != c1[2] and c2[0] == c2[1] != c2[2], f"{c1} {c2}")

def run(script, *args):
    r = subprocess.run([sys.executable, os.path.join(SCR, script), *args], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return r.stdout

def test_ref_core_and_origin():
    import csv
    with tempfile.TemporaryDirectory() as d:
        refs = os.path.join(d, "r.fa")
        open(refs, "w").write(">flankhit core=2000-3200\n" + "A" * 5200 + "\n>corehit core=2000-3200\n" + "A" * 5200 + "\n>nohit core=2000-3200\n" + "A" * 5200 + "\n")
        paf = os.path.join(d, "r.paf")
        open(paf, "w").write("\n".join([
            row("flankhit", 5200, 0, 2000, "+", "t", 90000, 0, 2000, "2000="),
            row("corehit", 5200, 1500, 3500, "+", "t", 90000, 5000, 7000, "2000=")]) + "\n")
        out = os.path.join(d, "o.tsv")
        run("ref_core_coverage.py", paf, refs, out)
        st = {r["reference"]: r["status"] for r in csv.DictReader(open(out), delimiter="\t")}
        check("flank-only hit is not an IS detection",
              st == {"flankhit": "flank_only", "corehit": "core_detected", "nohit": "no_hit"}, str(st))
        q = os.path.join(d, "q.fa")
        open(q, "w").write(">c1__100_2100\nA\n>lonely__0_500\nA\n")
        full = [row("c1__100_2100", 2000, 0, 2000, "+", f"x{i}", 9000, 0, 2000, "2000=") for i in range(20)]
        p1, p2 = os.path.join(d, "a.paf"), os.path.join(d, "b.paf")
        open(p1, "w").write("\n".join(full + [row("c1__100_2100", 2000, 0, 2000, "+", "c1", 9000, 100, 2100, "2000=")]) + "\n")
        open(p2, "w").write(row("c1__100_2100", 2000, 0, 300, "+", "y", 9000, 0, 300, "300=") + "\n")
        o2 = os.path.join(d, "s.tsv")
        run("fragment_origin_summary.py", q, o2, f"ecoli={p1}", f"kp={p2}", f"gone={d}/none.paf")
        rows = {r["query"]: r for r in csv.DictReader(open(o2), delimiter="\t")}
        r = rows["c1__100_2100"]
        check("origin summary: capped count, short hit and own locus not counted, no-hit row kept",
              r["ecoli_hits"] == "20+" and r["kp_hits"] == "0" and r["gone_hits"] == "NA"
              and r["best_match_organism"] == "ecoli" and rows["lonely__0_500"]["best_match_organism"] == "none",
              str(rows))

def test_final_review_tables():
    import csv
    with tempfile.TemporaryDirectory() as d:
        # 11. identity over the core, not over the row
        refs = os.path.join(d, "r.fa")
        open(refs, "w").write(">lowcore core=2000-3200\n" + "A" * 5200 + "\n>goodcore core=2000-3200\n" + "A" * 5200 + "\n")
        paf = os.path.join(d, "r.paf")
        open(paf, "w").write("\n".join([
            row("lowcore", 5200, 0, 5200, "+", "t", 90000, 0, 5200, "2000=180X3020="),
            row("goodcore", 5200, 0, 5200, "+", "t", 90000, 9000, 14200, "1800=180X3220=")]) + "\n")
        out = os.path.join(d, "o.tsv")
        run("ref_core_coverage.py", paf, refs, out)
        st = {r["reference"]: (r["status"], r["best_core_identity"]) for r in csv.DictReader(open(out), delimiter="\t")}
        check("a core at 85% identity inside a 96.5% row is not core_detected",
              st == {"lowcore": ("core_low_identity", "NA"), "goodcore": ("core_detected", "100.0")}, str(st))
        # 15. cap decided on raw rows
        q = os.path.join(d, "q.fa")
        open(q, "w").write(">c1__100_2100\nA\n")
        rows = [row("c1__100_2100", 2000, 0, 2000, "+", "c1", 9000, 100, 2100, "2000=")]
        rows += [row("c1__100_2100", 2000, 0, 2000, "+", f"x{i}", 9000, 0, 2000, "1980=20X") for i in range(5)]
        p1, p2 = os.path.join(d, "a.paf"), os.path.join(d, "b.paf")
        open(p1, "w").write("\n".join(rows) + "\n")
        open(p2, "w").write("\n".join(rows[:4]) + "\n")
        o2 = os.path.join(d, "s.tsv")
        run("fragment_origin_summary.py", q, o2, f"capped={p1}", f"open={p2}")
        r = next(csv.DictReader(open(o2), delimiter="\t"))
        check("self + 5 external rows at -N 5 is '5+', 4 rows is an exact 3",
              (r["capped_hits"], r["open_hits"], r["total_hits"]) == ("5+", "3", "8+"), str(r))
        # 14. extractor dedup does not depend on the order of the table
        rnd = random.Random(3)
        logan = os.path.join(d, "logan"); os.mkdir(logan)
        fa = os.path.join(d, "S1.contigs.fa")
        open(fa, "w").write(">S1_1\n" + "".join(rnd.choice("ACGT") for _ in range(5000)) + "\n")
        subprocess.run(["zstd", "-q", fa, "-o", os.path.join(logan, "S1.contigs.fa.zst")], check=True)
        hdr = ("type query qlen target_contig tlen part1_qs part1_qe part1_ts part1_te part1_strand part1_ident "
               "part2_qs part2_qe part2_ts part2_te part2_strand part2_ident target_gap query_gap cds_coverage").split()
        def n_frag(offsets):
            t = os.path.join(d, "ins.tsv")
            with open(t, "w") as f:
                f.write("\t".join(hdr) + "\n")
                for k, off in enumerate(offsets):
                    f.write("\t".join(map(str, ["insertion", f"q{k}", 2000, "S1_1", 5000, 0, 1000, 0, 1000 + off, "+", 100,
                                                1000, 2000, 2000 + off, 3000, "+", 100, 1000, 0, 100])) + "\n")
            o = os.path.join(d, "f.fa")
            run("extract_insertion_fragments.py", t, logan, o)
            return sorted(l.split()[0] for l in open(o) if l.startswith(">"))
        res = [n_frag(o) for o in ((0, 20, 21), (20, 0, 21), (21, 20, 0))]
        check("extractor: offsets 0/20/21 give the same two loci in any order",
              res[0] == res[1] == res[2] == [">S1_1__1000_2000", ">S1_1__1021_2021"], str(res))

def test_stage_lib():
    with tempfile.TemporaryDirectory() as d:
        sh = f"""set -euo pipefail
source "{SCR}/stage_lib.sh"
: > "{d}/empty.fa"
N=$(count_fasta "{d}/empty.fa"); echo "n=$N"
S1=$(stage_sig a b); S2=$(stage_sig a c)
stage_done "{d}/m" "$S1" && echo before=done || echo before=todo
stage_mark "{d}/m" "$S1"
stage_done "{d}/m" "$S1" && echo same=done || echo same=todo
stage_done "{d}/m" "$S2" && echo changed=done || echo changed=todo
"""
        r = subprocess.run(["bash", "-c", sh], capture_output=True, text=True)
        check("stage markers and zero counts",
              r.returncode == 0 and r.stdout.split() == ["n=0", "before=todo", "same=done", "changed=todo"],
              r.stdout + r.stderr)

def sh(script, env, *args):
    e = dict(os.environ, IS110_FRAGMENT_SCRIPTS=SCR, SLURM_CPUS_PER_TASK="2", **env)
    return subprocess.run(["bash", os.path.join(SCR, script), *args], capture_output=True, text=True, env=e)

def test_final_review_shell():
    rnd = random.Random(31)
    dna = lambda n: "".join(rnd.choice("ACGT") for _ in range(n))
    with tempfile.TemporaryDirectory() as d:
        # 16. a marker without its product is not "done"; 10. files_fp sees a replaced file
        g = os.path.join(d, "g.fna"); open(g, "w").write(">c\n" + "A" * 2000 + "\n")
        code = f"""set -euo pipefail
source "{SCR}/stage_lib.sh"
stage_mark "{d}/out.done" S
stage_done "{d}/out.done" S && echo missing=done || echo missing=todo
: > "{d}/out"
stage_done "{d}/out.done" S && echo empty=done || echo empty=todo
stage_done "{d}/out.done" S "{d}/out" "{d}/companion" && echo comp=done || echo comp=todo
A=$(echo "{g}" | files_fp)
printf '>c\\n%s\\n' "$(head -c 3000 /dev/zero | tr '\\0' A)" > "{g}"
B=$(echo "{g}" | files_fp)
[ "$A" != "$B" ] && echo fp=changed || echo fp=same
"""
        r = subprocess.run(["bash", "-c", code], capture_output=True, text=True)
        check("stage_done needs its products; files_fp changes with the file",
              r.stdout.split() == ["missing=todo", "empty=done", "comp=todo", "fp=changed"], r.stdout + r.stderr)
        check("hmm_ncbi_ecoli signs the genome files, not only their paths",
              'cut -f2 "$MANIFEST" | files_fp' in open(os.path.join(SCR, "hmm_ncbi_ecoli.sh")).read())
        # 8. the same contig name in two genome files stops both searches
        IS = dna(2000)
        gt = os.path.join(d, "gtdb"); os.mkdir(gt)
        open(os.path.join(gt, "a.fna"), "w").write(">ctg\n" + dna(500) + IS[:1000] + dna(6500) + "\n")
        open(os.path.join(gt, "b.fna"), "w").write(">ctg\n" + dna(7999) + IS[1000:] + dna(500) + "\n")
        q = os.path.join(d, "q.fna"); open(q, "w").write(">IS\n" + IS + "\n")
        o1 = os.path.join(d, "cds"); os.mkdir(o1)
        r = sh("run_cds_minimap2.sh", dict(QUERY=q, GTDB_DIR=gt, OUTDIR=o1, CPUS="2"))
        check("run_cds_minimap2.sh refuses duplicate contig names",
              r.returncode != 0 and "more than one genome file" in r.stderr and not os.path.exists(os.path.join(o1, "split_cds.tsv")),
              r.stderr[-300:])
        o2 = os.path.join(d, "gt2"); os.mkdir(o2)
        open(os.path.join(o2, "split_is110_all.tsv"), "w").write("query\nIS\n")
        for n in ("is110_cds_nr95.fna", "is110_cds_nuc.fna"):
            open(os.path.join(o2, n), "w").write(">IS\n" + IS + "\n")
        r = sh("search_split_vs_gtdb.sh", dict(OUTDIR=o2, GTDB_FNA=gt, CPUS="2"))
        check("search_split_vs_gtdb.sh refuses duplicate contig names",
              r.returncode != 0 and "more than one genome file" in r.stderr and not os.path.exists(os.path.join(o2, "split_123_all_hits.tsv")),
              r.stderr[-300:])
        # 12 + 7. reference extraction: flat genome layout through the manifest; new positions -> new representatives
        gd = os.path.join(d, "genomes"); os.mkdir(gd)
        contig = dna(20000)
        fna = os.path.join(gd, "GCF_000000001.1_genomic.fna"); open(fna, "w").write(">c1\n" + contig + "\n")
        hd = os.path.join(d, "hmm"); os.mkdir(hd)
        pos = os.path.join(hd, "is110_positions.tsv")
        def positions(rows):
            open(pos, "w").write("protein_id\tassembly\tcontig\tstart\tend\tstrand\tlen\n" +
                                 "".join(f"{pid}\tGCF_000000001.1\tc1\t{a}\t{b}\t+\t400\n" for pid, a, b in rows))
        positions([("p1", 3001, 4200), ("p2", 12001, 13200)])
        od = os.path.join(d, "refs")
        env = dict(POSITIONS_TSV=pos, GENOME_DIR=gd, OUTDIR=od)
        r = sh("extract_is110_refs.sh", env)
        check("flat genome layout without a manifest: an error, not an empty result",
              r.returncode != 0 and "positions without a sequence" in r.stderr and not os.path.exists(os.path.join(od, "is110_with_flanks.fa")),
              r.stdout[-300:] + r.stderr[-300:])
        open(os.path.join(hd, "genomes_one_per_assembly.tsv"), "w").write(f"GCF_000000001.1\t{fna}\n")
        r = sh("extract_is110_refs.sh", env)
        ids = lambda f: sorted(l[1:].split()[0] for l in open(os.path.join(od, f)) if l.startswith(">"))
        check("flat genome layout with the census manifest: both records",
              r.returncode == 0 and ids("is110_with_flanks.fa") == ["p1", "p2"] and ids("is110_reps_95.fa") == ["p1", "p2"],
              r.stdout[-300:] + r.stderr[-300:])
        r = sh("extract_is110_refs.sh", env)
        check("unchanged inputs: both steps skipped", r.returncode == 0 and r.stdout.count("skipping") == 2, r.stdout[-400:])
        positions([("p3", 7001, 8200)])
        r = sh("extract_is110_refs.sh", env)
        check("changed positions: references AND representatives are rebuilt",
              r.returncode == 0 and ids("is110_with_flanks.fa") == ["p3"] and ids("is110_reps_95.fa") == ["p3"],
              r.stdout[-300:] + r.stderr[-300:])
        # 7. fragment origin: a new insertion table gives new fragments and new representatives
        res = os.path.join(d, "res"); os.mkdir(res)
        logan = os.path.join(d, "logan"); os.mkdir(logan)
        fa = os.path.join(d, "S1.contigs.fa"); open(fa, "w").write(">S1_1\n" + dna(9000) + "\n")
        subprocess.run(["zstd", "-q", fa, "-o", os.path.join(logan, "S1.contigs.fa.zst")], check=True)
        hdr = ("type query qlen target_contig tlen part1_qs part1_qe part1_ts part1_te part1_strand part1_ident "
               "part2_qs part2_qe part2_ts part2_te part2_strand part2_ident target_gap query_gap cds_coverage").split()
        def table(te, ts):
            with open(os.path.join(res, "is110_x_95_insertions.tsv"), "w") as f:
                f.write("\t".join(hdr) + "\n" + "\t".join(map(str, ["insertion", "q", 2000, "S1_1", 9000, 0, 1000, te - 1000, te, "+", 100,
                                                                     1000, 2000, ts, ts + 1000, "+", 100, ts - te, 0, 100])) + "\n")
        env = dict(RESULTS_DIR=res, LOGAN_DIR=logan)
        reps = lambda: sorted(l[1:].split()[0] for l in open(os.path.join(res, "fragments_x_cluster_rep.fasta")) if l.startswith(">"))
        table(1000, 2000)
        r1 = sh("find_fragment_origin.sh", env, "x"); a = reps() if r1.returncode == 0 else r1.stderr[-300:]
        r2 = sh("find_fragment_origin.sh", env, "x")
        table(4000, 5200)
        r3 = sh("find_fragment_origin.sh", env, "x"); b = reps() if r3.returncode == 0 else r3.stderr[-300:]
        check("fragment origin: same table skips, new table rebuilds fragments and representatives",
              a == ["S1_1__1000_2000"] and r2.stdout.count("skipping") == 2 and b == ["S1_1__4000_5200"], f"{a} {b} {r2.stdout[-300:]}")

def test_minimap2_round_trip():
    """Real minimap2 (the options of the pipeline scripts) -> detector."""
    rnd = random.Random(11)
    dna = lambda n: "".join(rnd.choice("ACGT") for _ in range(n))
    rc = lambda x: x.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    IS = dna(1500)
    contigs = {
        "ins_plus": dna(3000) + IS[:700] + dna(600) + IS[700:] + dna(3000),
        "ins_minus": dna(3000) + rc(IS[:900] + dna(2500) + IS[900:]) + dna(3000),
        "del250": dna(3000) + IS[:600] + IS[850:] + dna(3000),
        "intact": dna(3000) + IS + dna(3000),
        "two": dna(2000) + IS[:700] + dna(800) + IS[700:] + dna(5000) + IS[:400] + dna(1200) + IS[400:] + dna(2000),
        "inverted_half": dna(3000) + IS[:700] + dna(600) + rc(IS[700:]) + dna(3000),
    }
    with tempfile.TemporaryDirectory() as d:
        q, t, paf = (os.path.join(d, n) for n in ("q.fa", "t.fa", "x.paf"))
        open(q, "w").write(f">IS\n{IS}\n")
        open(t, "w").write("".join(f">{k}\n{v}\n" for k, v in contigs.items()))
        for preset in ("asm10", "asm5"):
            r = subprocess.run(["minimap2", "-x", preset, "-c", "--eqx", "--secondary=yes", "-N", "1000000",
                                "-p", "0.05", "-g", "200", "-r", "200,200", t, q, "-o", paf],
                               capture_output=True, text=True)
            assert r.returncode == 0, r.stderr
            subprocess.run([sys.executable, os.path.join(SCR, "detect_is110_variants.py"), paf,
                            os.path.join(d, "o")], capture_output=True, text=True, check=True)
            got = sorted((l.split("\t")[3], l.split("\t")[0])
                         for l in open(os.path.join(d, "o_all.tsv")).read().splitlines()[1:])
            want = [("del250", "deletion"), ("ins_minus", "insertion"), ("ins_plus", "insertion"),
                    ("two", "insertion"), ("two", "insertion")]
            check(f"minimap2 {preset} round trip: events on the right contigs only", got == want, str(got))

for t in (test_identity, test_eqx_is_required, test_output_independent_of_hash_seed, test_detector, test_final_review_detector, test_multisplit, test_lineage, test_ref_core_and_origin, test_final_review_tables, test_stage_lib, test_final_review_shell,
          test_minimap2_round_trip):


    t()
print("\nALL PASS" if not fails else f"\n{len(fails)} FAILED: {fails}")
sys.exit(1 if fails else 0)
