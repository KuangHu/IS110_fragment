#!/usr/bin/env python3
"""Confirmation scripts for Stage 7 candidates, on synthetic genomes with
known answers. Needs minimap2 + samtools:  python3 tests/test_rearr.py
"""
import csv
import json
import os
import random
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


if __name__ == "__main__":
    test_parse_details_and_loci()
    test_translocation_bilateral_all_orientations()
    test_inversion_flank_check()
