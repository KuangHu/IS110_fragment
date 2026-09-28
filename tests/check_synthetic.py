#!/usr/bin/env python3
"""Check a pipeline run on tests/make_synthetic.py data against the truth.

    check_synthetic.py <synthetic_dir> <run_dir>
"""
import csv
import json
import os
import sys

syn, run = sys.argv[1], sys.argv[2]
truth = json.load(open(os.path.join(syn, "truth.json")))
fails = []


def check(cond, msg):
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        fails.append(msg)


# --- Stage 5: counts are distinct assemblies; tolerance does not scale with D
rows = list(csv.DictReader(open(os.path.join(run, "boundaries_summary.tsv")), delimiter="\t"))
by = {(r["ref_id"], int(r["distance_D"])): r for r in rows}
for D in (1000, 5000):
    r1, r2 = by.get(("IS1", D)), by.get(("IS2", D))
    check(r1 is not None and r2 is not None, f"stage5 rows exist at D={D}")
    if r1 and r2:
        check(int(r1["n_empty"]) == 2, f"IS1 D={D}: 2 empty assemblies (got {r1['n_empty']})")
        check(int(r1["n_v0"]) == 1, f"IS1 D={D}: 1 filled assembly, GCA/GCF twin counted once (got {r1['n_v0']})")
        check(int(r1["n_v1plus"]) == 1, f"IS1 D={D}: 1 cargo-carrying assembly (got {r1['n_v1plus']})")
        check(int(r2["n_empty"]) == 2 and int(r2["n_v0"]) == 2,
              f"IS2 D={D}: 2 empty + 2 filled (got {r2['n_empty']}/{r2['n_v0']})")
        check(int(r1["n_pairs"]) == truth["n_distinct_assemblies"],
              f"IS1 D={D}: one placement per distinct assembly (got {r1['n_pairs']})")
        mode = int(r1["is_length_mode"] or 0)
        check(abs(mode - truth["IS1"]["length"]) <= 5,
              f"IS1 D={D}: element length from empties ~{truth['IS1']['length']} (got {mode})")

# --- Stage 6: base-level boundaries by decomposition, both strands
recs = {r["ref_id"]: r for r in json.load(open(os.path.join(run, "records", "records.json")))}
for rid in ("IS1", "IS2"):
    r = recs.get(rid)
    check(r is not None, f"{rid}: record built")
    if not r:
        continue
    be, ie, t = r["boundary_evidence"], r["is_element"], truth[rid]
    check(be["method"] == "empty_vs_filled_decomposition",
          f"{rid}: boundary by decomposition (got {be['method']})")
    mh = (be.get("decomposition") or {}).get("junction_microhomology_bp", 0) or 0
    check(ie["length"] == t["length"], f"{rid}: element length {t['length']} (got {ie['length']})")
    d5 = ie["start_offset_5p"] - t["off5"]
    check(0 <= d5 <= mh, f"{rid}: 5' offset {t['off5']} within microhomology {mh} "
                         f"(got {ie['start_offset_5p']})")
    d3 = ie["end_offset_3p"] - t["off3"]
    check(d3 == d5, f"{rid}: 3' offset shifts with 5' (slide) (got {ie['end_offset_3p']})")
    check(0 <= ie["source_start"] - t["source_start_1"] <= mh
          or 0 <= t["source_end"] - ie["source_end"] <= mh,
          f"{rid}: genomic coords match truth ({ie['source_start']}-{ie['source_end']} vs "
          f"{t['source_start_1']}-{t['source_end']})")
    check(be["n_empty_observations"] == 2, f"{rid}: 2 empty assemblies in record")

# --- Stage 8: abutting-anchor empties are kept; twins collapsed; nested insert keyed
obs = json.load(open(os.path.join(run, "variants", "observations.json")))
for rid in ("IS1", "IS2"):
    o = [x for x in obs if x["v1_parent_id"] == rid]
    cats = [x["category"] for x in o]
    check(cats.count("empty") == 2, f"{rid}: 2 empty observations incl. 0-bp gaps (got {cats.count('empty')})")
    check(len({x["assembly"] for x in o}) == len(o), f"{rid}: one observation per assembly")
ins = [x for x in obs if x["v1_parent_id"] == "IS1" and x["category"] == "insertion"]
check(len(ins) == 1, f"IS1: 1 nested-insertion variant (got {len(ins)})")
if ins:
    ni = ins[0].get("nested_insert") or {}
    check(abs(ni.get("insert_len", 0) - truth["cargo_len"]) <= ni.get("junction_microhomology_bp", 0),
          f"IS1: nested insert {truth['cargo_len']} bp (got {ni.get('insert_len')})")
    check(bool(ni.get("event_key")), "IS1: nested insert has an event_key")

# --- Stage 9: lineages carry no direction claim
lin = json.load(open(os.path.join(run, "lineages", "lineages.json")))
check(all(v.get("direction") == "unpolarized" for v in lin.values()) and lin,
      "lineages are labelled unpolarized")

print()
print(f"{'ALL CHECKS PASSED' if not fails else f'{len(fails)} CHECK(S) FAILED'}")
sys.exit(1 if fails else 0)
