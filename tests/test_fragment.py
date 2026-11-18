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
            rows += [row(q, 2000, 0, 1000, "+", f"t{k}", 9000, 0, 1000 + k, f"{1000 + k}="),
                     row(q, 2000, 1000, 2000, "+", f"t{k}", 9000, 1600 + k, 2600, f"{1000 - k}=")]
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
        rows += [row(f"q{k}", 2000, 0, 1000, "+", "t", 9000, 0, 1000 + off, f"{1000 + off}="),
                 row(f"q{k}", 2000, 1000, 2000, "+", "t", 9000, 1600 + off, 2600, f"{1000 - off}=")]
    out = detect(rows, "--max-junction-overlap", "100")
    check("junction offsets 0/20/21 are one locus", len(out) == 1, str(len(out)))

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
              r["ecoli_hits"] == "6+" and r["kp_hits"] == "0" and r["gone_hits"] == "NA"
              and r["best_match_organism"] == "ecoli" and rows["lonely__0_500"]["best_match_organism"] == "none",
              str(rows))

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

for t in (test_identity, test_eqx_is_required, test_output_independent_of_hash_seed, test_detector, test_multisplit, test_lineage, test_ref_core_and_origin, test_stage_lib,
          test_minimap2_round_trip):


    t()
print("\nALL PASS" if not fails else f"\n{len(fails)} FAILED: {fails}")
sys.exit(1 if fails else 0)
