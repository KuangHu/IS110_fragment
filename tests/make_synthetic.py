#!/usr/bin/env python3
"""Build a small synthetic genome DB with IS elements at KNOWN boundaries.

    background G (60 kb)
    IS1 at 30,000 on +  : 150 bp 5' end + 1200 bp "tnp" + 100 bp 3' end
    IS2 at 45,000 on -  : 80 bp 5' end + 1100 bp "tnp" + 220 bp 3' end

    GCA_000000001.1  filled (both IS)            } GenBank/RefSeq twins:
    GCF_000000001.1  identical copy of the above } must count ONCE
    GCA_000000002.1  empty, 0.3% SNPs            -> tolerant decomposition
    GCA_000000003.1  empty, exact
    GCA_000000005.1  IS1 carries an 800 bp nested cargo insertion; IS2 filled

Writes <out>/db.fa (+ .fai via samtools), <out>/is_hits.tsv, <out>/truth.json.
"""
import json
import os
import random
import subprocess
import sys

out = sys.argv[1]
os.makedirs(out, exist_ok=True)
rng = random.Random(11)


def rseq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def mutate(s, rate):
    out_ = list(s)
    for i in range(len(out_)):
        if rng.random() < rate:
            out_[i] = rng.choice([b for b in "ACGT" if b != out_[i]])
    return "".join(out_)


G = rseq(60000)
is1 = {"five": rseq(150), "tnp": rseq(1200), "three": rseq(100), "site": 30000, "strand": "+"}
is2 = {"five": rseq(80), "tnp": rseq(1100), "three": rseq(220), "site": 45000, "strand": "-"}
cargo = rseq(800)
CARGO_AT = 700                       # offset inside IS1


def element(e, with_cargo=False):
    s = e["five"] + e["tnp"] + e["three"]
    if with_cargo:
        s = s[:CARGO_AT] + cargo + s[CARGO_AT:]
    return s


def build(bg, cargo_is1=False, empty=False):
    """Insert IS2 first (downstream) so IS1's coordinates are unaffected."""
    if empty:
        return bg
    s1 = element(is1, cargo_is1)
    s2 = revcomp(element(is2))
    return (bg[:is1["site"]] + s1 + bg[is1["site"]:is2["site"]] + s2 + bg[is2["site"]:])


genomes = {
    "GCA_000000001.1|ctgA": build(G),
    "GCF_000000001.1|NZ_ctgA": build(G),
    "GCA_000000002.1|ctgB": mutate(G, 0.003),
    "GCA_000000003.1|ctgC": G,
    "GCA_000000005.1|ctgE": build(G, cargo_is1=True),
}
with open(os.path.join(out, "db.fa"), "w") as f:
    for name, seq in genomes.items():
        f.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i + 80] + "\n")
subprocess.run(["samtools", "faidx", os.path.join(out, "db.fa")], check=True)

# transposase coordinates in the filled genome A (1-based, inclusive)
len1 = len(element(is1))
t1s = is1["site"] + len(is1["five"]) + 1
t1e = t1s + len(is1["tnp"]) - 1
# IS2 is reverse-complemented; it starts at site + len1 (after IS1 insert)
is2_start0 = is2["site"] + len1
len2 = len(element(is2))
# on the + strand the element reads revcomp(three) + revcomp(tnp) + revcomp(five)
t2s = is2_start0 + len(is2["three"]) + 1
t2e = t2s + len(is2["tnp"]) - 1

with open(os.path.join(out, "is_hits.tsv"), "w") as f:
    f.write("is_id\tassembly\tcontig\ttnp_start\ttnp_end\ttnp_strand\ttnp_len\tdomains_hit\n")
    f.write(f"IS1\tGCA_000000001.1\tGCA_000000001.1|ctgA\t{t1s}\t{t1e}\t+\t{t1e - t1s + 1}\tsynthetic\n")
    f.write(f"IS2\tGCA_000000001.1\tGCA_000000001.1|ctgA\t{t2s}\t{t2e}\t-\t{t2e - t2s + 1}\tsynthetic\n")

truth = {
    "IS1": {"off5": -len(is1["five"]), "off3": len(is1["three"]), "length": len1,
            "source_start_1": is1["site"] + 1, "source_end": is1["site"] + len1},
    "IS2": {"off5": -len(is2["five"]), "off3": len(is2["three"]), "length": len2,
            "source_start_1": is2_start0 + 1, "source_end": is2_start0 + len2},
    "cargo_len": len(cargo),
    "n_distinct_assemblies": 4,
}
json.dump(truth, open(os.path.join(out, "truth.json"), "w"), indent=2)
print(json.dumps(truth, indent=2))
