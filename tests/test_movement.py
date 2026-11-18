#!/usr/bin/env python3
"""is_large_element_movement_db.py on a synthetic species DB with known answers.

Needs minimap2 + samtools on PATH:  python3 tests/test_movement.py
"""
import csv
import os
import random
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(HERE, "..", "scripts")
sys.path.insert(0, SCRIPTS)
from lib_alleles import revcomp  # noqa: E402

rng = random.Random(7)
r = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731

A, B, C, D, E_ = r(20000), r(20000), r(20000), r(20000), r(20000)   # backbone blocks
EL = r(8000)                                                    # the element
OTHER = r(30000)

genomes = {
    # reference: element between A and B; site C|D is empty
    "GCA_000000001.1": {"c1": A + EL + B + C + D + E_},
    # its RefSeq twin: identical -> must be ignored entirely
    "GCF_000000001.1": {"c1": A + EL + B + C + D + E_},
    # same site, same orientation
    "GCA_000000002.1": {"c1": A + EL + B + C + D + E_},
    # same site, contig assembled in the other orientation (used to vanish)
    "GCA_000000003.1": {"c1": revcomp(A + EL + B + C + D + E_)},
    # element at C|D instead; A|B empty  -> MOVED, reference site empty
    "GCA_000000004.1": {"c1": A + B + C + EL + D + E_},
    # same site but the contig ends 300 bp after the element -> not assessable
    "GCA_000000005.1": {"c1": A + EL + B[:300]},
    # moved AND reverse orientation of the whole contig
    "GCA_000000006.1": {"c1": revcomp(A + B + C + EL + D + E_)},
    # element in unrelated sequence the reference does not have
    "GCA_000000007.1": {"c1": OTHER[:15000] + EL + OTHER[15000:]},
}


def main():
    with tempfile.TemporaryDirectory() as wd:
        db = os.path.join(wd, "db.fa")
        with open(db, "w") as f:
            for asm, ctgs in genomes.items():
                for c, s in ctgs.items():
                    f.write(f">{asm}|{c}\n{s}\n")
        subprocess.run(["samtools", "faidx", db], check=True)
        ref = os.path.join(wd, "ref.fa")
        with open(ref, "w") as f:
            f.write(f">c1\n{genomes['GCA_000000001.1']['c1']}\n")
        subprocess.run(["samtools", "faidx", ref], check=True)
        rec = os.path.join(wd, "recs.tsv")
        with open(rec, "w") as f:
            f.write("ref_id\tassembly\tcontig\tis_start\tis_end\tlength\n")
            f.write(f"R1\tGCA_000000001.1\tc1\t{len(A) + 1}\t{len(A) + len(EL)}\t{len(EL)}\n")
        out = os.path.join(wd, "out")
        subprocess.run([sys.executable, os.path.join(SCRIPTS, "is_large_element_movement_db.py"),
                        "--ref-asm", "GCA_000000001.1", "--ref-fa", ref, "--species-db", db,
                        "--records-tsv", rec, "--out", out, "--threads", "2"],
                       check=True, stdout=subprocess.DEVNULL)
        rows = list(csv.DictReader(open(os.path.join(out, "GCA_000000001.1_large_movement.tsv")),
                                   delimiter="\t"))
    got = {r_["tgt_asm"]: r_ for r_ in rows}
    assert len(rows) == len(got), "one row per target copy"
    want = {"GCA_000000002.1": "SYNTENIC", "GCA_000000003.1": "SYNTENIC",
            "GCA_000000004.1": "MOVED", "GCA_000000005.1": "FLANK_UNASSESSABLE",
            "GCA_000000006.1": "MOVED", "GCA_000000007.1": "MOVED"}
    assert "GCF_000000001.1" not in got and "GCA_000000001.1" not in got, \
        "reference and its GCF twin must be excluded"
    for asm, v in want.items():
        assert asm in got, f"{asm}: no row (minus-strand copies used to be dropped)"
        assert got[asm]["verdict"] == v, (asm, got[asm])
        print(f"PASS {asm}: {v} strand={got[asm]['strand']} ref_site={got[asm]['ref_site']}")
    assert got["GCA_000000004.1"]["ref_site"] == "empty"
    assert got["GCA_000000006.1"]["ref_site"] == "empty"
    assert got["GCA_000000007.1"]["ref_site"] == "absent"
    assert got["GCA_000000003.1"]["strand"] == "-"
    assert all(float(r_["IS_cov"]) > 0.99 for r_ in rows), "1-based start: whole element queried"
    print("PASS ref_site states, twin exclusion, coverage")


def paralog():
    """Reference has copies at X and Y; target has the same two. Neither moved."""
    g = {"GCA_000000011.1": {"c1": A + EL + B + C + EL + D},
         "GCA_000000012.1": {"c1": A + EL + B + C + EL + D}}
    with tempfile.TemporaryDirectory() as wd:
        db = os.path.join(wd, "db.fa")
        with open(db, "w") as f:
            for asm, ctgs in g.items():
                f.write(f">{asm}|c1\n{ctgs['c1']}\n")
        subprocess.run(["samtools", "faidx", db], check=True)
        ref = os.path.join(wd, "ref.fa")
        with open(ref, "w") as f:
            f.write(f">c1\n{g['GCA_000000011.1']['c1']}\n")
        subprocess.run(["samtools", "faidx", ref], check=True)
        rec = os.path.join(wd, "recs.tsv")
        y = len(A) + len(EL) + len(B) + len(C)
        with open(rec, "w") as f:
            f.write("ref_id\tassembly\tcontig\tis_start\tis_end\tlength\n")
            f.write(f"X\tGCA_000000011.1\tc1\t{len(A) + 1}\t{len(A) + len(EL)}\t{len(EL)}\n")
            f.write(f"Y\tGCA_000000011.1\tc1\t{y + 1}\t{y + len(EL)}\t{len(EL)}\n")
        out = os.path.join(wd, "out")
        subprocess.run([sys.executable, os.path.join(SCRIPTS, "is_large_element_movement_db.py"),
                        "--ref-asm", "GCA_000000011.1", "--ref-fa", ref, "--species-db", db,
                        "--records-tsv", rec, "--out", out, "--threads", "2"],
                       check=True, stdout=subprocess.DEVNULL)
        rows = list(csv.DictReader(open(os.path.join(out, "GCA_000000011.1_large_movement.tsv")),
                                   delimiter="\t"))
    v = sorted(r_["verdict"] for r_ in rows)
    assert v == ["SYNTENIC", "SYNTENIC", "SYNTENIC_OTHER_REF_COPY", "SYNTENIC_OTHER_REF_COPY"], v
    print("PASS two copies in both genomes: no MOVED (was 2 MOVED)")


if __name__ == "__main__":
    main()
    paralog()
