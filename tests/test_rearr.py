#!/usr/bin/env python3
"""Confirmation scripts for Stage 7 candidates, on synthetic genomes with
known answers. Needs minimap2 + samtools:  python3 tests/test_rearr.py
"""
import csv
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(HERE, "..", "scripts")
sys.path.insert(0, SCRIPTS)
from lib_alleles import revcomp  # noqa: E402
import lib_rearr as lr  # noqa: E402

rng = random.Random(21)
r = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
A, B, C, X, IS = r(40000), r(40000), r(40000), r(40000), r(1500)


def find(hay, needle):
    """(start, end, strand) of needle in hay, either strand."""
    i = hay.find(needle)
    if i >= 0:
        return i, i + len(needle), "+"
    i = hay.find(revcomp(needle))
    assert i >= 0
    return i, i + len(needle), "-"


def details(genome, up_anchor, down_anchor):
    """What Stage 7 would write: where the two (transposase-oriented) anchors land."""
    out = {}
    for side, anc in (("up", up_anchor), ("down", down_anchor)):
        for c, s in genome.items():
            if anc in s or revcomp(anc) in s:
                a, b, st = find(s, anc)
                out[side] = f"{c}:{a}-{b}({st})"
    return f"up=[{out['up']}] down=[{out['down']}]"


def setup(wd, ref_seq, is_strand, targets):
    """Reference contig c1 = A + IS + B (IS possibly on the minus strand)."""
    src = os.path.join(wd, "src"); os.makedirs(src, exist_ok=True)
    with open(os.path.join(src, "GCA_000000001.1.fa"), "w") as f:
        f.write(f">c1\n{ref_seq}\n")
    for asm, ctgs in targets.items():
        with open(os.path.join(src, asm + ".fa"), "w") as f:
            for c, s in ctgs.items():
                f.write(f">{c}\n{s}\n")
    s0 = len(A)
    rec = [{"ref_id": "GCA_000000001.1|c1_1",
            "source": {"assembly": "GCA_000000001.1", "contig": "GCA_000000001.1|c1",
                       "transposase_start": s0 + 1, "transposase_end": s0 + len(IS),
                       "transposase_strand": is_strand},
            "is_element": {"source_start": s0 + 1, "source_end": s0 + len(IS),
                           "length": len(IS)}}]
    rp = os.path.join(wd, "records.json")
    json.dump(rec, open(rp, "w"))
    # anchors as extract_anchors makes them: 500 bp, 1 kb from the element,
    # in transposase orientation
    left, right = ref_seq[s0 - 1500:s0 - 1000], ref_seq[s0 + len(IS) + 1000:s0 + len(IS) + 1500]
    up, down = (left, right) if is_strand == "+" else (revcomp(right), revcomp(left))
    return src, rp, up, down


def run_translocation(is_strand, flip_c2):
    ref = A + (IS if is_strand == "+" else revcomp(IS)) + B
    c2 = X + IS + B
    tgt = {"t1": A + IS + C, "t2": revcomp(c2) if flip_c2 else c2}
    with tempfile.TemporaryDirectory() as wd:
        src, rp, up, down = setup(wd, ref, is_strand, {"GCA_000000002.1": tgt,
                                                       "GCA_000000003.1": {"t1": ref}})
        cand = os.path.join(wd, "cand.tsv")
        with open(cand, "w") as f:
            f.write("category\tref_id\tassembly\tdetails\n")
            f.write(f"translocation_only\tGCA_000000001.1|c1_1\tGCA_000000002.1\t"
                    f"{details(tgt, up, down)}\n")
        out = os.path.join(wd, "out")
        subprocess.run([sys.executable, os.path.join(SCRIPTS, "is_translocation_check.py"),
                        "--candidates", cand, "--records", rp, "--src-dir", src,
                        "--out", out, "--threads", "2"], check=True, stdout=subprocess.DEVNULL)
        return list(csv.DictReader(open(os.path.join(out, "per_candidate_translocation.tsv")),
                                   delimiter="\t"))[0]


def test_translocation_bilateral_all_orientations():
    for is_strand in "+-":
        for flip_c2 in (False, True):
            row = run_translocation(is_strand, flip_c2)
            assert row["passes"] == "True", (is_strand, flip_c2, row)
            print(f"PASS translocation: IS strand {is_strand}, contig flipped={flip_c2}")


def test_parse_details_and_loci():
    d = lr.parse_details("up=[contig-12:100-600(+);lcl|c1:plasmid:1-5(-)] down=[]")
    assert d["up"] == [("contig-12", 100, 600, "+"), ("lcl|c1:plasmid", 1, 5, "-")], d
    assert d["down"] == []
    print("PASS parse_details keeps contig names with '-' and ':'")


def run_inversion(target_seq, is_strand="+", which="a"):
    """ref = A IS M rc(IS) B. Reference element = IS_a (which='a') or IS_b."""
    M = C
    ref = A + IS + M + revcomp(IS) + B
    with tempfile.TemporaryDirectory() as wd:
        src = os.path.join(wd, "src"); os.makedirs(src)
        with open(os.path.join(src, "GCA_000000001.1.fa"), "w") as f:
            f.write(f">c1\n{ref}\n")
        with open(os.path.join(src, "GCA_000000002.1.fa"), "w") as f:
            f.write(f">t1\n{target_seq}\n")
        if which == "a":
            s0, strand = len(A), "+"
        else:
            s0, strand = len(A) + len(IS) + len(M), "-"
        rec = [{"ref_id": "GCA_000000001.1|c1_1",
                "source": {"assembly": "GCA_000000001.1", "contig": "GCA_000000001.1|c1",
                           "transposase_strand": strand},
                "is_element": {"source_start": s0 + 1, "source_end": s0 + len(IS)}}]
        rp = os.path.join(wd, "records.json"); json.dump(rec, open(rp, "w"))
        left, right = ref[s0 - 1500:s0 - 1000], ref[s0 + len(IS) + 1000:s0 + len(IS) + 1500]
        up, down = (left, right) if strand == "+" else (revcomp(right), revcomp(left))
        det = details({"t1": target_seq}, up, down)
        cand = os.path.join(wd, "cand.tsv"); val = os.path.join(wd, "val.tsv")
        with open(cand, "w") as f:
            f.write("category\tref_id\tassembly\tdetails\testimated_inverted_bp\n")
            # a second target for the same ref: must not replace this one's anchors
            f.write(f"inversion_only\tGCA_000000001.1|c1_1\tGCA_000000002.1\t{det}\t{len(M)}\n")
            f.write("inversion_only\tGCA_000000001.1|c1_1\tGCA_000000009.1\t"
                    "up=[zz:1-500(+)] down=[zz:900-1400(-)]\t5000\n")
        with open(val, "w") as f:
            f.write("category\tref_id\ttgt_asm\tverdict\n")
            f.write("inversion_only\tGCA_000000001.1|c1_1\tGCA_000000002.1\tCONFIRMED_BOTH\n")
        out = os.path.join(wd, "out")
        subprocess.run([sys.executable, os.path.join(SCRIPTS, "inversion_revcomp_check.py"),
                        "--validation", val, "--candidates", cand, "--records", rp,
                        "--src-dir", src, "--out", out, "--threads", "2"],
                       check=True, stdout=subprocess.DEVNULL)
        return list(csv.DictReader(open(os.path.join(out, "per_candidate_revcomp.tsv")),
                                   delimiter="\t"))[0]


def test_inversion_flank_check():
    M = C
    ref = A + IS + M + revcomp(IS) + B
    inv = A + IS + revcomp(M) + revcomp(IS) + B
    for which in "ab":
        for name, tgt, want in (("true inversion", inv, "True"),
                                ("true inversion, contig reversed", revcomp(inv), "True"),
                                ("no inversion", ref, "NA"),
                                ("no inversion, contig reversed", revcomp(ref), "NA")):
            row = run_inversion(tgt, which=which)
            assert row["passes"] == want, (which, name, row)
            print(f"PASS inversion check, ref element IS_{which}: {name} -> {row['passes']}")


def is_positions(seq):
    """1-based (start, end, strand) of every IS copy in seq."""
    out = []
    for q, st in ((IS, "+"), (revcomp(IS), "-")):
        i = seq.find(q)
        while i >= 0:
            out.append((i + 1, i + len(IS), st)); i = seq.find(q, i + 1)
    return sorted(out)


def run_dup(ref, s0, strand, tgt_seq, canonical=False):
    with tempfile.TemporaryDirectory() as wd:
        src = os.path.join(wd, "src"); os.makedirs(src)
        with open(os.path.join(src, "GCA_000000001.1.fa"), "w") as f:
            f.write(f">c1\n{ref}\n")
        with open(os.path.join(src, "GCA_000000002.1.fa"), "w") as f:
            f.write(f">t1\n{tgt_seq}\n")
        rec = [{"ref_id": "GCA_000000001.1|c1_1",
                "source": {"assembly": "GCA_000000001.1", "contig": "GCA_000000001.1|c1",
                           "transposase_strand": strand},
                "is_element": {"source_start": s0 + 1, "source_end": s0 + len(IS)}}]
        rp = os.path.join(wd, "records.json"); json.dump(rec, open(rp, "w"))
        ih = os.path.join(wd, "is_hits.tsv")
        with open(ih, "w") as f:
            f.write("is_id\tassembly\tcontig\ttnp_start\ttnp_end\ttnp_strand\n")
            for k, (a_, b_, st) in enumerate(is_positions(tgt_seq)):
                f.write(f"x{k}\tGCA_000000002.1\tGCA_000000002.1|t1\t{a_}\t{b_}\t{st}\n")
        left, right = ref[s0 - 1500:s0 - 1000], ref[s0 + len(IS) + 1000:s0 + len(IS) + 1500]
        up, down = (left, right) if strand == "+" else (revcomp(right), revcomp(left))
        best_side, hits = None, []
        for side, anc in (("up", up), ("down", down)):
            h = []
            for q, st in ((anc, "+"), (revcomp(anc), "-")):
                i = tgt_seq.find(q)
                while i >= 0:
                    h.append(f"t1:{i}-{i + 500}({st})"); i = tgt_seq.find(q, i + 1)
            if len(h) >= 2:
                best_side, hits = side, h
        assert best_side, "fixture has no anchor multiplicity"
        det = f"up=[{';'.join(hits)}] down=[]" if best_side == "up" else f"up=[] down=[{';'.join(hits)}]"
        cand = os.path.join(wd, "cand.tsv")
        with open(cand, "w") as f:
            f.write("category\tref_id\tassembly\tdetails\tbest_side\n")
            f.write(f"duplication\tGCA_000000001.1|c1_1\tGCA_000000002.1\t{det}\t{best_side}\n")
        out = os.path.join(wd, "out")
        if canonical:
            c2 = os.path.join(wd, "strict.tsv")
            with open(c2, "w") as f:
                f.write("species\tref_id\ttgt_asm\n")
                f.write("test_species\tGCA_000000001.1|c1_1\tGCA_000000002.1\n")
            subprocess.run([sys.executable,
                            os.path.join(SCRIPTS, "duplication_strict_canonical_check.py"),
                            "--candidates", c2, "--dup-all-input", cand, "--records", rp,
                            "--is-hits", ih, "--src-dir", src, "--out", out, "--threads", "2"],
                           check=True, stdout=subprocess.DEVNULL)
            rows = list(csv.DictReader(open(os.path.join(out, "duplication_canonical_check.tsv")),
                                       delimiter="\t"))
            return rows[0] if rows else {"canonical_pass": "NA"}
        subprocess.run([sys.executable, os.path.join(SCRIPTS, "is_duplication_check.py"),
                        "--candidates", cand, "--records", rp, "--src-dir", src,
                        "--is-hits", ih, "--out", out, "--threads", "2"],
                       check=True, stdout=subprocess.DEVNULL)
        return list(csv.DictReader(open(os.path.join(out, "per_candidate_duplication.tsv")),
                                   delimiter="\t"))[0]


def test_duplication_check():
    S = B[:12000]
    rest = B[12000:]
    for strand in "+-":
        el = IS if strand == "+" else revcomp(IS)
        # reference: A [IS] S rest      (for '-': the mirror image, so the
        # duplicated segment is still on the element's DOWN side)
        ref = A + el + S + rest if strand == "+" else revcomp(A + IS + S + rest)
        s0 = ref.find(el)
        dup = A + IS + S + IS + S + IS + rest
        cases = (("true duplication", dup, "True"),
                 ("true duplication, contig reversed", revcomp(dup), "True"),
                 # no duplicated segment: just a second copy of the anchor
                 # region next to another IS far away
                 ("anchor region copied only", A + IS + S + rest + C + IS + S[500:2500] + X, "False"),
                 # second copy inverted: not a direct repeat
                 ("second copy inverted", A + IS + S + IS + revcomp(S) + IS + rest, "NA"))
        for name, tgt, want in cases:
            row = run_dup(ref, s0, strand, tgt)
            assert row["passes"] == want, (strand, name, row)
            print(f"PASS duplication check, element strand {strand}: {name} -> {row['passes']} "
                  f"(dup_bp={row['dup_bp']})")
            crow = run_dup(ref, s0, strand, tgt, canonical=True)
            assert crow["canonical_pass"] == want, (strand, name, crow)
            print(f"PASS canonical check,   element strand {strand}: {name} -> "
                  f"{crow['canonical_pass']}")


def run_inv_full(ref, s0, strand, tgt):
    with tempfile.TemporaryDirectory() as wd:
        src = os.path.join(wd, "src"); os.makedirs(src)
        with open(os.path.join(src, "GCA_000000001.1.fa"), "w") as f:
            f.write(f">c1\n{ref}\n")
        with open(os.path.join(src, "GCA_000000002.1.fa"), "w") as f:
            f.write(f">t1\n{tgt}\n")
        rec = [{"ref_id": "GCA_000000001.1|c1_1",
                "source": {"assembly": "GCA_000000001.1", "contig": "GCA_000000001.1|c1",
                           "transposase_strand": strand},
                "is_element": {"source_start": s0 + 1, "source_end": s0 + len(IS)}}]
        rp = os.path.join(wd, "records.json"); json.dump(rec, open(rp, "w"))
        ih = os.path.join(wd, "is_hits.tsv")
        with open(ih, "w") as f:
            f.write("is_id\tassembly\tcontig\ttnp_start\ttnp_end\ttnp_strand\n")
            for k, (a_, b_, st) in enumerate(is_positions(tgt)):
                f.write(f"x{k}\tGCA_000000002.1\tGCA_000000002.1|t1\t{a_}\t{b_}\t{st}\n")
        left, right = ref[s0 - 1500:s0 - 1000], ref[s0 + len(IS) + 1000:s0 + len(IS) + 1500]
        up, down = (left, right) if strand == "+" else (revcomp(right), revcomp(left))
        anc = os.path.join(wd, "anc.tsv"); cand = os.path.join(wd, "cand.tsv")
        with open(anc, "w") as f:
            f.write("ref_id\tassembly\tdetails\n")
            f.write(f"GCA_000000001.1|c1_1\tGCA_000000002.1\t{details({'t1': tgt}, up, down)}\n")
        with open(cand, "w") as f:
            f.write("ref_id\ttgt_asm\nGCA_000000001.1|c1_1\tGCA_000000002.1\n")
        out = os.path.join(wd, "out")
        subprocess.run([sys.executable, os.path.join(SCRIPTS, "inversion_full_check.py"),
                        "--candidates", cand, "--anchor-tsv", anc, "--records", rp,
                        "--is-hits", ih, "--src-dir", src, "--out", out, "--threads", "2"],
                       check=True, stdout=subprocess.DEVNULL)
        rows = list(csv.DictReader(open(os.path.join(out, "per_candidate_inversion_full.tsv")),
                                   delimiter="\t"))
        return rows[0] if rows else {"all_pass": "NA"}


def test_inversion_full_check():
    n = 0
    for mlen in (8000, 12000, 20000, 40000):
        M = C[:mlen]
        ref = A + IS + M + revcomp(IS) + B
        inv = A + IS + revcomp(M) + revcomp(IS) + B
        for which, (s0, strand) in (("a", (len(A), "+")), ("b", (len(A) + len(IS) + mlen, "-"))):
            for tgt in (inv, revcomp(inv)):
                row = run_inv_full(ref, s0, strand, tgt)
                assert row["all_pass"] == "True", (mlen, which, row)
                n += 1
    print(f"PASS full inversion check: {n}/16 true inversions pass (was 0/16)")
    # only 3 kb of a 40 kb interval inverted: must not pass
    M = C
    ref = A + IS + M + revcomp(IS) + B
    part = A + IS + M[:1000] + revcomp(M[1000:4000]) + M[4000:] + revcomp(IS) + B
    try:
        row = run_inv_full(ref, len(A), "+", part)
    except AssertionError:
        row = {"all_pass": "NA"}
    assert row["all_pass"] != "True", row
    print(f"PASS full inversion check: partial 3 kb inversion -> {row['all_pass']}")


def test_tandem_db_scan():
    S = B[:12000]
    genomes = {"GCA_000000021.1": A + (IS + S) * 2 + IS + C,       # 2 copies
               "GCA_000000022.1": A + (IS + S) * 3 + IS + C,       # 3 (was: no row)
               "GCA_000000023.1": A + (IS + S) * 4 + IS + C,       # 4 (was: unit doubled)
               "GCF_000000023.1": A + (IS + S) * 4 + IS + C,       # twin of the above
               "GCA_000000024.1": A + IS + S + IS + X[:4000] + S + IS + C,  # not adjacent
               "GCA_000000025.1": A + IS + S + IS + C}             # no duplication
    with tempfile.TemporaryDirectory() as wd:
        db = os.path.join(wd, "db.fa"); ih = os.path.join(wd, "is_hits.tsv")
        with open(db, "w") as f, open(ih, "w") as g:
            g.write("is_id\tassembly\tcontig\ttnp_start\ttnp_end\ttnp_strand\n")
            for asm, sq in genomes.items():
                f.write(f">{asm}|c1\n{sq}\n")
                for k, (a_, b_, st) in enumerate(is_positions(sq)):
                    g.write(f"{asm}_{k}\t{asm}\t{asm}|c1\t{a_}\t{b_}\t{st}\n")
        subprocess.run(["samtools", "faidx", db], check=True)
        out = os.path.join(wd, "out")
        subprocess.run([sys.executable, os.path.join(SCRIPTS, "duplication_tandem_db_scan.py"),
                        "--species-db", db, "--is-hits", ih, "--out", out, "--threads", "2"],
                       check=True, stdout=subprocess.DEVNULL)
        f = [x for x in os.listdir(out) if x.endswith(".tsv")][0]
        rows = list(csv.DictReader(open(os.path.join(out, f)), delimiter="\t"))
    got = {(r_["assembly"][4:], int(r_["n_copies"]), r_["verdict"]) for r_ in rows}
    cores = sorted(r_["assembly"][4:] for r_ in rows)
    assert cores == ["000000021.1", "000000022.1", "000000023.1"], rows
    for core, n in (("000000021.1", 2), ("000000022.1", 3), ("000000023.1", 4)):
        assert (core, n, "VERIFIED_CANONICAL") in got, (core, n, got)
    assert all(abs(int(r_["copy_bp"]) - 13500) <= 20 for r_ in rows), rows
    print("PASS tandem DB scan: 2/3/4-copy arrays exact, twin collapsed, "
          "separated copies and single copy not reported")


MAUVE = os.path.expanduser("~/.conda/envs/mauve-env/bin/progressiveMauve")


def run_gold(cat, ref, tgt, s0, strand="+", order=("up", "down")):
    """validate_rearrangements_gold.py on one candidate; returns its row."""
    with tempfile.TemporaryDirectory() as wd:
        src = os.path.join(wd, "src"); os.makedirs(src)
        with open(os.path.join(src, "GCA_000000001.1.fa"), "w") as f:
            f.write(f">c1\n{ref}\n")
        with open(os.path.join(src, "GCA_000000002.1.fa"), "w") as f:
            for c, sq in tgt.items():
                f.write(f">{c}\n{sq}\n")
        rec = [{"ref_id": "GCA_000000001.1|c1_1",
                "source": {"assembly": "GCA_000000001.1", "contig": "GCA_000000001.1|c1",
                           "transposase_strand": strand},
                "is_element": {"source_start": s0 + 1, "source_end": s0 + len(IS)}}]
        rp = os.path.join(wd, "records.json"); json.dump(rec, open(rp, "w"))
        left, right = ref[s0 - 1500:s0 - 1000], ref[s0 + len(IS) + 1000:s0 + len(IS) + 1500]
        up, down = (left, right) if strand == "+" else (revcomp(right), revcomp(left))
        if cat == "duplication":       # the down anchor hits twice
            hits = []
            for c, sq in tgt.items():
                i = sq.find(down)
                while i >= 0:
                    hits.append(f"{c}:{i}-{i + 500}(+)"); i = sq.find(down, i + 1)
            a, b, st = find(next(iter(tgt.values())), up)
            det = f"up=[t1:{a}-{b}({st})] down=[{';'.join(hits)}]"
        else:
            det = details(tgt, up, down)
        ex = os.path.join(wd, "ex.tsv")
        with open(ex, "w") as f:
            f.write("category\tref_id\tassembly\tup_contigs\tdown_contigs\tdetails\n")
            f.write(f"{cat}\tGCA_000000001.1|c1_1\tGCA_000000002.1\t-\t-\t{det}\n")
        out = os.path.join(wd, "out")          # absolute; a relative one is tested below
        subprocess.run([sys.executable, os.path.join(SCRIPTS, "validate_rearrangements_gold.py"),
                        "--examples", ex, "--src-dir", src, "--records", rp, "--out", out,
                        "--mauve", MAUVE], check=True, capture_output=True)
        return list(csv.DictReader(open(os.path.join(out, "validation.tsv")), delimiter="\t"))[0]


def test_gold_validator():
    if not os.path.exists(MAUVE) or not shutil.which("dnadiff"):
        print("SKIP gold validator (progressiveMauve / dnadiff not available)")
        return
    M = C
    ref = A + IS + M + revcomp(IS) + B
    s0 = len(A)
    cases = [
        ("inversion_only", "true 40 kb inversion", {"t1": A + IS + revcomp(M) + revcomp(IS) + B}, True),
        ("inversion_only", "true inversion, contig stored reversed",
         {"t1": revcomp(A + IS + revcomp(M) + revcomp(IS) + B)}, True),
        ("inversion_only", "only 3 kb inverted, 18 kb from the IS",
         {"t1": A + IS + M[:18000] + revcomp(M[18000:21000]) + M[21000:] + revcomp(IS) + B}, False),
    ]
    for cat, name, tgt, want in cases:
        # the 3 kb case needs anchors as a candidate would have them: fake by
        # using the true-inversion anchors' positions when the down anchor is not inverted
        try:
            row = run_gold(cat, ref, tgt, s0)
        except AssertionError:
            continue
        got = row["verdict"] == "CONFIRMED_BOTH"
        assert got == want, (name, row)
        print(f"PASS gold {cat}: {name} -> {row['verdict']} [{row['mauve_metric']}]")
    # duplication: real one vs a window with just a 300-bp deletion
    ref2 = A + IS + B
    S = B[:12000]
    dup = A + IS + S + IS + S + B[12000:]
    row = run_gold("duplication", ref2, {"t1": dup}, len(A))
    assert row["verdict"] == "CONFIRMED_BOTH", row
    print(f"PASS gold duplication: true 12 kb duplication -> {row['verdict']} [{row['mauve_metric']}]")
    nodup = A + IS + B[:13500] + B[13800:]
    row = run_gold("duplication", ref2, {"t1": nodup}, len(A))
    assert row["verdict"] != "CONFIRMED_BOTH" and not row["verdict"].startswith("CONFIRMED_DNADIFF"), row
    print(f"PASS gold duplication: 300-bp deletion only -> {row['verdict']} [{row['mauve_metric']}]")
    # translocation: same verdict whatever the element strand (= contig listing order)
    for strand in "+-":
        refx = A + (IS if strand == "+" else revcomp(IS)) + B
        row = run_gold("translocation_only", refx, {"t1": A + IS + C, "t2": X + IS + B},
                       len(A), strand)
        assert row["verdict"] == "CONFIRMED_BOTH", (strand, row)
        print(f"PASS gold translocation, element strand {strand} -> {row['verdict']} "
              f"[{row['mauve_metric']}]")


if __name__ == "__main__":
    test_parse_details_and_loci()
    test_translocation_bilateral_all_orientations()
    test_inversion_flank_check()
    test_inversion_full_check()
    test_duplication_check()
    test_tandem_db_scan()
    test_gold_validator()
