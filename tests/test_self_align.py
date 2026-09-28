#!/usr/bin/env python3
"""detect_tandem.self_minimap2 must see a high-copy tandem array (-f 0).

Needs minimap2 on PATH:  python3 tests/test_self_align.py
"""
import os
import random
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))
import detect_tandem  # noqa: E402

rng = random.Random(3)
r = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
unit = r(300)
seq = r(3000) + unit * 40 + r(3000)          # 40-copy tandem array
with tempfile.TemporaryDirectory() as wd:
    hits = detect_tandem.self_minimap2(seq, wd, 2)
arr0, arr1 = 3000, 3000 + 40 * 300
inside = [h for h in hits if arr0 <= h[0] and h[1] <= arr1 + 50]
assert inside, "no self-alignments inside the 40-copy tandem array"
print(f"PASS self-alignment finds the tandem array ({len(inside)} hits)")
