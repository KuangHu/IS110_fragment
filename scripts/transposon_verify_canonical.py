#!/usr/bin/env python3
"""Verify the canonical [IS][DNA][IS] structure of MOVED events in tgt.

For each MOVED verdict, check whether tgt has IS hits within --is-tol bp of
BOTH the start and end of the moved DNA segment.

  PASS_CANONICAL: IS within --is-tol of BOTH segment boundaries in tgt
                  → confirmed composite transposon structure
  PASS_PARTIAL:   IS within --is-tol of only ONE boundary
                  → unilateral IS-flanking; ambiguous mechanism
  FAIL:           no IS within --is-tol of either boundary
                  → DNA at new location but NOT bracketed by IS → not a
                    composite transposon mobilization (could be HGT cargo etc.)
"""
import argparse, csv, os, sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_alleles import assembly_core  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--moved-tsv", required=True,
                   help="klebs_transposon_MOVED_catalogue.tsv")
    p.add_argument("--is-hits",   required=True)
    p.add_argument("--out",       required=True)
    p.add_argument("--is-tol", type=int, default=3000)
    return p.parse_args()


def load_is_hits(path):
    """assembly → contig (local) → sorted list of (s, e)"""
    idx = defaultdict(lambda: defaultdict(list))
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            asm = r["assembly"]; contig = r["contig"]
            local = contig.split("|", 1)[1] if "|" in contig else contig
            idx[asm][local].append((int(r["tnp_start"]), int(r["tnp_end"])))
    for asm in idx:
        for c in idx[asm]:
            idx[asm][c].sort()
    return idx


def bracketing_is(hits, ts, te, tol):
    """Is the interval [ts, te] bracketed by TWO DISTINCT IS copies?

    start: an IS hit overlapping [ts - tol, ts + tol]; end: a DIFFERENT hit
    overlapping [te - tol, te + tol] that lies after the first. The old test
    looked for 'any hit within tol of ts' and 'any hit within tol of te'
    independently: one IS in the middle of a 5.5 kb interval satisfied both
    and the event was CANONICAL."""
    near_s = [h for h in hits if h[0] <= ts + tol and h[1] >= ts - tol]
    near_e = [h for h in hits if h[0] <= te + tol and h[1] >= te - tol]
    both = any(a != b and a[1] < b[0] for a in near_s for b in near_e)
    return bool(near_s), bool(near_e), both


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)
    is_hits = load_is_hits(args.is_hits)
    print(f"loaded is_hits for {len(is_hits)} assemblies", flush=True)

    out_path = os.path.join(args.out, "moved_canonical_verified.tsv")
    n = {"CANONICAL":0, "PARTIAL_IS":0, "NO_IS_AT_BOUNDARY":0, "NO_HITS_FOR_TGT":0}
    with open(args.moved_tsv) as fin, open(out_path, "w") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        hdr = reader.fieldnames + ["tgt_IS_at_start", "tgt_IS_at_end", "canonical_verdict"]
        fout.write("\t".join(hdr) + "\n")
        for r in reader:
            tgt_asm = r["tgt_asm"]
            tgt_contig = r["tgt_contig"]
            ts = int(r["ts"]); te = int(r["te"])
            tgt_hits = is_hits.get(tgt_asm, {}).get(tgt_contig, [])
            if not tgt_hits:
                # is_hits.tsv keeps one of a GCA/GCF pair: look up the twin
                twin = next((a for a in is_hits if a != tgt_asm
                             and assembly_core(a) == assembly_core(tgt_asm)), None)
                if twin:
                    tgt_hits = (is_hits[twin].get(tgt_contig)
                                or is_hits[twin].get("NZ_" + tgt_contig)
                                or is_hits[twin].get(tgt_contig[3:] if tgt_contig.startswith("NZ_") else "")
                                or [])
            if not tgt_hits:
                v_start = v_end = False
                verdict = "NO_HITS_FOR_TGT"
            else:
                start_ok, end_ok, both = bracketing_is(tgt_hits, ts, te, args.is_tol)
                v_start, v_end = start_ok, end_ok
                if both:
                    verdict = "CANONICAL"
                elif start_ok or end_ok:
                    verdict = "PARTIAL_IS"
                else:
                    verdict = "NO_IS_AT_BOUNDARY"
            n[verdict] += 1
            fout.write("\t".join(r[c] for c in reader.fieldnames)
                       + f"\t{v_start}\t{v_end}\t{verdict}\n")

    print(f"\nDONE.")
    total = sum(n.values())
    for v, c in sorted(n.items(), key=lambda x: -x[1]):
        pct = 100 * c / total if total else 0.0
        print(f"  {v:<20} {c:>7} ({pct:.1f}%)")
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
