#!/usr/bin/env python3
"""Unit tests for stage-level logic that the synthetic end-to-end run does not
reach (contig ends, prodigal IDs, cache markers, vote normalisation).

    python3 tests/test_stage_units.py
"""
import os
import random
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
import build_records  # noqa: E402
import extract_anchors  # noqa: E402
import is_detect  # noqa: E402
import pipeline  # noqa: E402
from lib_alleles import decompose_alleles, hamming_split, revcomp  # noqa: E402

rng = random.Random(3)


def rseq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(s, rate):
    out = list(s)
    for i in range(len(out)):
        if rng.random() < rate:
            out[i] = rng.choice([b for b in "ACGT" if b != out[i]])
    return "".join(out)


def test_revcomp_iupac():
    assert revcomp("ACGTRYKMSWBDHVN") == "NBDHVWSKMRYACGT"
    assert revcomp(revcomp("ARYKMNNacgt")) == "ARYKMNNacgt"


def test_hamming_split():
    empty = rseq(10000)
    ins = rseq(1500)
    off = 4200
    filled = empty[:off] + ins + empty[off:]
    snp = mutate(empty, 0.004)
    hs = hamming_split(snp, filled, min_ident=0.95, min_insert=50)
    assert hs and hs["insert_len"] == 1500
    assert abs(hs["offset"] - off) <= hs["ambiguity"] + 1, hs
    # through decompose_alleles, with the quadratic alignment switched off
    r = decompose_alleles(snp, filled, min_insert=50, tolerant=False)
    assert r["method"] == "hamming_split" and r["insert_len"] == 1500
    assert filled[r["insert_start"]:r["insert_start"] + r["insert_len"]] == r["insert_seq"]
    # and it stays off when asked
    assert decompose_alleles(snp, filled, tolerant=False,
                             hamming_fallback=False)["method"] == "none"
    # unrelated sequence: refused
    assert hamming_split(rseq(10000), filled) is None
    # a 3-bp deletion elsewhere shifts the frame -> dense mismatch run -> refused
    # (accepting it would make the insert 3 bp too long)
    shifted = empty[:2000] + empty[2003:]
    assert hamming_split(shifted, filled, min_insert=50) is None
    # too short an insert
    assert hamming_split(empty, empty[:500] + "ACGTACGTAC" + empty[500:]) is None


def test_extract_anchors_contig_end():
    """A hit near a contig end keeps the D values that fit (it used to be dropped),
    on both strands, and the anchors sit exactly D bp from the transposase."""
    with tempfile.TemporaryDirectory() as td:
        g = rseq(30000)
        fa = os.path.join(td, "g.fa")
        with open(fa, "w") as f:
            f.write(">A|c1\n" + g + "\n")
        subprocess.run(["samtools", "faidx", fa], check=True)
        extract_anchors._worker_init(fa, 20000, 500, [1000, 5000, 10000])
        # tnp at 3001-4200 (1-based): 3000 bp on the left, 25800 on the right
        for strand in "+-":
            res = extract_anchors._extract_one(("x", "A|c1", 3001, 4200, strand))
            assert res is not None, "hit near contig end must not be dropped"
            ref, anchors, row, n = res
            assert row.rstrip().split("\t")[5] == "1000", row      # only D=1000 fits
            assert n == 2
            seqs = dict(zip(*[iter(anchors.strip().split("\n"))] * 2))
            up = [v for k, v in seqs.items() if "__up1000" in k][0]
            down = [v for k, v in seqs.items() if "__down1000" in k][0]
            left, right = g[1500:2000], g[5200:5700]     # genomic anchors at D=1000
            if strand == "+":
                assert (up, down) == (left, right)
            else:
                assert (up, down) == (revcomp(right), revcomp(left))
        # nothing fits -> skipped
        assert extract_anchors._extract_one(("y", "A|c1", 801, 2000, "+")) is None
        # plenty of room -> all distances
        res = extract_anchors._extract_one(("z", "A|c1", 14001, 15200, "+"))
        assert res[2].rstrip().split("\t")[5] == "1000,5000,10000"


def test_parse_gff_both_callers():
    with tempfile.TemporaryDirectory() as td:
        p = os.path.join(td, "x.gff")
        with open(p, "w") as f:
            # prodigal CLI: ID is <contig ordinal>_<gene index>
            f.write("GCA_1.1|ctgA\tProdigal_v2.6.3\tCDS\t10\t99\t1.0\t+\t0\tID=1_1;partial=00;\n")
            f.write("GCA_1.1|ctgB\tProdigal_v2.6.3\tCDS\t5\t64\t1.0\t-\t0\tID=2_1;partial=00;\n")
            # pyrodigal: ID is <contig>_<gene index>
            f.write("GCA_2.1|c_x\tpyrodigal_v3\tCDS\t7\t66\t1.0\t+\t0\tID=GCA_2.1|c_x_3;partial=00;\n")
        info = is_detect.parse_gff(p)
        # keys must equal the protein FASTA headers of both callers
        assert info["GCA_1.1|ctgA_1"] == ("GCA_1.1|ctgA", 10, 99, "+")
        assert info["GCA_1.1|ctgB_1"] == ("GCA_1.1|ctgB", 5, 64, "-")
        assert info["GCA_2.1|c_x_3"] == ("GCA_2.1|c_x", 7, 66, "+")


def test_stage_markers():
    with tempfile.TemporaryDirectory() as td:
        out = os.path.join(td, "o.txt")
        later = os.path.join(td, ".stage_4_pairs.done")
        cmd = [sys.executable, "-c", f"open({out!r}, 'a').write('x')"]
        run = lambda **kw: pipeline._run_stage(td, "3_minimap2", "t", cmd, [out], **kw)  # noqa: E731
        # a partial output with no marker is NOT trusted: the stage runs
        open(out, "w").write("partial")
        open(later, "w").close()
        run()
        assert open(out).read() == "partialx"
        assert os.path.exists(os.path.join(td, ".stage_3_minimap2.done"))
        assert not os.path.exists(later), "downstream marker must be invalidated"
        run()                                   # marker present -> skipped
        assert open(out).read() == "partialx"
        os.remove(os.path.join(td, ".stage_3_minimap2.done"))
        run(trust_existing=True)                # adopted, not re-run
        assert open(out).read() == "partialx"
        # a failing stage leaves no marker
        os.remove(os.path.join(td, ".stage_3_minimap2.done"))
        bad = [sys.executable, "-c", "import sys; sys.exit(3)"]
        try:
            pipeline._run_stage(td, "3_minimap2", "t", bad, [out])
            raise AssertionError("expected failure")
        except subprocess.CalledProcessError:
            pass
        assert not os.path.exists(os.path.join(td, ".stage_3_minimap2.done"))
        # atomic -o: the tool writes <path>.tmp, renamed only on success
        o2 = os.path.join(td, "paf")
        c2 = [sys.executable, "-c", "import sys; open(sys.argv[1],'w').write('ok')", o2]
        pipeline._run_stage(td, "3_minimap2", "t", c2, [o2], atomic_o=3)
        assert open(o2).read() == "ok" and not os.path.exists(o2 + ".tmp")


class _FakeFasta:
    def __init__(self, seqs):
        self.seqs = seqs

    def fetch(self, contig, start, end):
        return self.seqs[contig][start:end]

    def get_reference_length(self, contig):
        return len(self.seqs[contig])


def test_votes_are_placement_invariant():
    """With junction microhomology the exact and the alignment-based methods
    place the same insert differently; both must vote for ONE boundary."""
    D, mh = 1000, "GATTACAGG"
    left, right = rseq(1500), rseq(1500)
    body = rseq(1300)
    tnp = body[100:1200]                       # element = 100 + 1100 tnp + 100
    # element followed by a copy of its own first 9 bp => 9 bp slide
    element = mh + body
    filled_g = left + element + mh + right
    empty_exact = left + mh + right
    empty_snp = mutate(left, 0.004) + mh + mutate(right, 0.004)
    tnp_start0 = len(left) + len(mh) + 100
    src = {"contig": "S", "start": tnp_start0 + 1, "end": tnp_start0 + len(tnp), "strand": "+"}
    build_records._FA = _FakeFasta({"S": filled_g, "E1": empty_exact, "E2": empty_snp})

    g0 = tnp_start0 - D
    obs = []
    for name, seq, ident in (("E1", empty_exact, 100.0), ("E2", empty_snp, 99.6)):
        # the empty interval between the same anchors: [g0, len - (flank right of anchor))
        e = {"tname": name, "strand": "+", "mean_ident": ident, "assembly": name,
             "up_te": g0, "down_ts": len(seq) - (len(filled_g) - (src["end"] + D))}
        e["distance"] = e["down_ts"] - e["up_te"]
        obs.append(e)
    dec = build_records.decompose_boundary(src, {"tnp_len": len(tnp)},
                                           {"D": D, "empty_obs": obs}, 5, 2e7)
    assert dec and dec["ok"], dec
    assert dec["n_decomposed"] == 2 and dec["n_distinct_calls"] == 1, dec
    assert dec["consensus_support"] == 1.0
    assert dec["length"] == len(element)
    assert dec["junction_microhomology_bp"] >= len(mh)


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"PASS {name}")
    print("all stage unit tests passed")
