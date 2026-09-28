#!/usr/bin/env python3
"""Generate paired GenBank files for visual review of one STRICT CANONICAL
composite-transposon movement event:

  REF:  [region_A] [IS_a] [DNA →] [IS_b] [region_B]
  TGT:  [region_C] [IS_a'] [DNA →] [IS_b'] [region_D]

The DNA is identical between ref and tgt (forward match). The flanking regions
A/B differ from C/D — that's the proof of movement to a new context.

Exactly 5 features per record: region_C/A, IS_a, DNA, IS_b, region_D/B.
"""
import argparse, csv, gzip, json, os
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--ref-asm", required=True)
    p.add_argument("--ref-contig", required=True)
    p.add_argument("--is-a-start", type=int, required=True)
    p.add_argument("--is-b-start", type=int, required=True)
    p.add_argument("--tgt-asm", required=True)
    p.add_argument("--tgt-contig", required=True)
    p.add_argument("--ts", type=int, required=True, help="DNA start in tgt (from MOVED catalogue)")
    p.add_argument("--te", type=int, required=True, help="DNA end in tgt")
    p.add_argument("--is-hits", required=True)
    p.add_argument("--src-dir", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--flank", type=int, default=5000)
    return p.parse_args()


def load_is_hits(path):
    idx = defaultdict(lambda: defaultdict(list))
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            asm = r["assembly"]; contig = r["contig"]
            local = contig.split("|", 1)[1] if "|" in contig else contig
            idx[asm][local].append((int(r["tnp_start"]), int(r["tnp_end"]),
                                    r.get("tnp_strand", "+")))
    return idx


def merge_is(hits, gap=500):
    if not hits: return []
    sh = sorted(hits, key=lambda x: x[0])
    merged = [sh[0]]
    for h in sh[1:]:
        if h[0] <= merged[-1][1] + gap:
            merged[-1] = (min(merged[-1][0], h[0]),
                          max(merged[-1][1], h[1]),
                          merged[-1][2])
        else:
            merged.append(h)
    return merged


def load_seqs(src_dir, asm):
    for ext in (".fna.gz", ".fa.gz", ".fna", ".fa"):
        path = os.path.join(src_dir, asm + ext)
        if os.path.exists(path):
            opener = gzip.open if path.endswith(".gz") else open
            recs = {}
            with opener(path, "rt") as fh:
                for rec in SeqIO.parse(fh, "fasta"):
                    local = rec.id.split("|", 1)[1] if "|" in rec.id else rec.id
                    recs[local] = rec
                    recs[rec.id] = rec
            return recs
    return {}


def slice_5feature(rec, sub_s, sub_e, is_a_rel, is_b_rel, label_prefix,
                   region_a_name, region_b_name):
    sub = rec[sub_s:sub_e]
    sub.id = (label_prefix + "_" + rec.id[:10])[:16]
    sub.name = sub.id
    sub.description = f"{rec.id} slice {sub_s+1}-{sub_e}"
    sub.annotations["molecule_type"] = "DNA"
    L = len(sub.seq)

    # region_A or region_C
    if is_a_rel[0] > 0:
        sub.features.append(SeqFeature(
            FeatureLocation(0, is_a_rel[0]), type="misc_feature",
            qualifiers={"label": [region_a_name],
                        "note": [f"{region_a_name} — conserved 5' flank"]}))
    # IS_a
    sub.features.append(SeqFeature(
        FeatureLocation(is_a_rel[0], is_a_rel[1],
                        strand=1 if is_a_rel[2] == "+" else -1),
        type="repeat_region",
        qualifiers={"label": ["IS_a"], "rpt_family": ["IS110"],
                    "note": ["IS110 transposase (IS_a)"]}))
    # DNA
    sub.features.append(SeqFeature(
        FeatureLocation(is_a_rel[1], is_b_rel[0]),
        type="misc_feature",
        qualifiers={"label": ["DNA"],
                    "note": ["mobile DNA cargo (= the composite transposon's middle)"]}))
    # IS_b
    sub.features.append(SeqFeature(
        FeatureLocation(is_b_rel[0], is_b_rel[1],
                        strand=1 if is_b_rel[2] == "+" else -1),
        type="repeat_region",
        qualifiers={"label": ["IS_b"], "rpt_family": ["IS110"],
                    "note": ["IS110 transposase (IS_b)"]}))
    # region_B or region_D
    if is_b_rel[1] < L:
        sub.features.append(SeqFeature(
            FeatureLocation(is_b_rel[1], L), type="misc_feature",
            qualifiers={"label": [region_b_name],
                        "note": [f"{region_b_name} — conserved 3' flank"]}))
    return sub


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    is_hits = load_is_hits(args.is_hits)

    # --- REF side ---
    ref_seqs = load_seqs(args.src_dir, args.ref_asm)
    ref_rec = ref_seqs.get(args.ref_contig)
    if ref_rec is None:
        raise SystemExit(f"ref contig {args.ref_contig} not found")
    ref_is_merged = merge_is(is_hits.get(args.ref_asm, {}).get(args.ref_contig, []))
    is_a_ref = next((h for h in ref_is_merged if abs(h[0] - args.is_a_start) < 500), None)
    is_b_ref = next((h for h in ref_is_merged if abs(h[0] - args.is_b_start) < 500), None)
    if is_a_ref is None or is_b_ref is None:
        raise SystemExit(f"ref IS not found: a={is_a_ref} b={is_b_ref}")
    if is_a_ref[0] > is_b_ref[0]:
        is_a_ref, is_b_ref = is_b_ref, is_a_ref
    ref_lo = max(0, is_a_ref[0] - args.flank)
    ref_hi = min(len(ref_rec.seq), is_b_ref[1] + args.flank)
    ref_is_a_rel = (is_a_ref[0]-ref_lo, is_a_ref[1]-ref_lo, is_a_ref[2])
    ref_is_b_rel = (is_b_ref[0]-ref_lo, is_b_ref[1]-ref_lo, is_b_ref[2])
    ref_sub = slice_5feature(ref_rec, ref_lo, ref_hi,
                             ref_is_a_rel, ref_is_b_rel,
                             "REF", "region_A", "region_B")
    print(f"REF slice: {ref_lo}-{ref_hi} ({ref_hi-ref_lo} bp)")
    print(f"  IS_a: {is_a_ref}  IS_b: {is_b_ref}")

    # --- TGT side ---
    tgt_seqs = load_seqs(args.src_dir, args.tgt_asm)
    tgt_rec = tgt_seqs.get(args.tgt_contig)
    if tgt_rec is None:
        raise SystemExit(f"tgt contig {args.tgt_contig} not found")
    tgt_is_merged = merge_is(is_hits.get(args.tgt_asm, {}).get(args.tgt_contig, []))
    # Find IS_a near ts and IS_b near te
    def near(hits, pos, tol=5000):
        best = None
        for h in hits:
            d = 0 if h[0] <= pos <= h[1] else min(abs(pos-h[0]), abs(pos-h[1]))
            if d <= tol and (best is None or d < best[3]):
                best = (h[0], h[1], h[2], d)
        return best
    is_a_tgt = near(tgt_is_merged, args.ts)
    is_b_tgt = near(tgt_is_merged, args.te)
    if is_a_tgt is None or is_b_tgt is None:
        raise SystemExit(f"tgt IS not found near boundaries: a={is_a_tgt} b={is_b_tgt}")
    is_a_tgt = (is_a_tgt[0], is_a_tgt[1], is_a_tgt[2])
    is_b_tgt = (is_b_tgt[0], is_b_tgt[1], is_b_tgt[2])
    if is_a_tgt[0] > is_b_tgt[0]:
        is_a_tgt, is_b_tgt = is_b_tgt, is_a_tgt
    tgt_lo = max(0, is_a_tgt[0] - args.flank)
    tgt_hi = min(len(tgt_rec.seq), is_b_tgt[1] + args.flank)
    tgt_is_a_rel = (is_a_tgt[0]-tgt_lo, is_a_tgt[1]-tgt_lo, is_a_tgt[2])
    tgt_is_b_rel = (is_b_tgt[0]-tgt_lo, is_b_tgt[1]-tgt_lo, is_b_tgt[2])
    tgt_sub = slice_5feature(tgt_rec, tgt_lo, tgt_hi,
                             tgt_is_a_rel, tgt_is_b_rel,
                             "TGT", "region_C", "region_D")
    print(f"TGT slice: {tgt_lo}-{tgt_hi} ({tgt_hi-tgt_lo} bp)")
    print(f"  IS_a: {is_a_tgt}  IS_b: {is_b_tgt}")

    # Write
    with open(os.path.join(args.out_dir, "ref.gbk"), "w") as fh:
        SeqIO.write([ref_sub], fh, "genbank")
    with open(os.path.join(args.out_dir, "tgt.gbk"), "w") as fh:
        SeqIO.write([tgt_sub], fh, "genbank")
    print(f"\nWrote {args.out_dir}/ref.gbk ({ref_hi-ref_lo} bp)")
    print(f"Wrote {args.out_dir}/tgt.gbk ({tgt_hi-tgt_lo} bp)")


if __name__ == "__main__":
    main()
