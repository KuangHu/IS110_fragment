#!/usr/bin/env python3
"""Cluster LARGE IS110 elements (whole element, not just transposase) by
nucleotide identity using cd-hit-est.

Input: records.json from the empty-vs-filled IS boundary pipeline.
       Each record has is_element.source_start / source_end / length / sequence.

Filter: is_element.length >= --min-element-bp (default 5000 bp).
Output: cluster reps + reps_records.tsv subset + rep_ref_assemblies.txt.
"""
import argparse, json, os, subprocess, sys
import pysam


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--records", required=True,
                   help="records.json from empty-vs-filled pipeline")
    p.add_argument("--db", required=True,
                   help="species concatenated FASTA, indexed with .fai")
    p.add_argument("--out", required=True)
    p.add_argument("--min-element-bp", type=int, default=5000)
    p.add_argument("--identity",  type=float, default=0.99)
    p.add_argument("--coverage",  type=float, default=0.95,
                   help="cd-hit-est -aS (cov of shorter seq)")
    p.add_argument("--threads",   type=int,   default=8)
    p.add_argument("--word-size", type=int,   default=10)
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    with open(args.records) as fh:
        recs = json.load(fh)

    # Filter to whole-element length >= threshold
    kept = []
    for r in recs:
        ie = r.get("is_element", {})
        L = ie.get("length", 0)
        if L < args.min_element_bp: continue
        src = r.get("source", {})
        if not all(k in src for k in ("assembly", "contig")):
            continue
        if not all(k in ie for k in ("source_start", "source_end")):
            continue
        kept.append({
            "ref_id": r["ref_id"],
            "assembly": src["assembly"],
            "contig":   src["contig"],
            "is_start": int(ie["source_start"]),
            "is_end":   int(ie["source_end"]),
            "length":   int(L),
            "sequence": ie.get("sequence", ""),
        })
    print(f"[cluster] {len(kept)} IS elements with length >= {args.min_element_bp} bp",
          flush=True)
    if not kept:
        sys.exit("No large IS elements found.")

    # Extract sequence from DB if missing (sanity-check)
    fa = pysam.FastaFile(args.db)
    fa_path = os.path.join(args.out, "large_is_bodies.fa")
    n_written = 0
    with open(fa_path, "w") as out:
        for r in kept:
            seq = r["sequence"]
            if not seq:
                try:
                    seq = fa.fetch(r["contig"], r["is_start"], r["is_end"])
                except (KeyError, ValueError):
                    continue
            if not seq: continue
            out.write(f">{r['ref_id']}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + "\n")
            n_written += 1
    fa.close()
    print(f"[cluster] wrote {n_written} body sequences to {fa_path}", flush=True)

    # cd-hit-est
    reps_fa = os.path.join(args.out, "reps.fa")
    cmd = ["cd-hit-est", "-i", fa_path, "-o", reps_fa,
           "-c", str(args.identity), "-aS", str(args.coverage),
           "-n", str(args.word_size), "-T", str(args.threads),
           "-M", "16000", "-d", "0"]
    print(f"[cluster] cd-hit-est: {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    clstr_path = reps_fa + ".clstr"

    # Parse cluster file
    clusters = []; cur = None
    with open(clstr_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">Cluster"):
                if cur: clusters.append(cur)
                cur = {"id": line.split()[1], "rep": None, "members": []}
            else:
                parts = line.split(", ")
                if len(parts) < 2: continue
                meta = parts[1]
                seq_id = meta.split(">", 1)[1].split("...")[0]
                cur["members"].append(seq_id)
                if "*" in line: cur["rep"] = seq_id
        if cur: clusters.append(cur)
    print(f"[cluster] {len(clusters)} clusters at {args.identity*100:.0f}% identity",
          flush=True)

    # Write reps-only records TSV (assembly, contig, is_start, is_end, length, ref_id)
    by_refid = {r["ref_id"]: r for r in kept}
    reps_tsv = os.path.join(args.out, "reps_records.tsv")
    with open(reps_tsv, "w") as fh:
        fh.write("ref_id\tassembly\tcontig\tis_start\tis_end\tlength\n")
        for cl in clusters:
            r = by_refid.get(cl["rep"])
            if not r: continue
            fh.write(f"{r['ref_id']}\t{r['assembly']}\t{r['contig']}\t"
                     f"{r['is_start']}\t{r['is_end']}\t{r['length']}\n")
    print(f"[cluster] wrote reps_records.tsv with {len(clusters)} reps", flush=True)

    # Distinct ref assemblies
    refs = set()
    with open(reps_tsv) as fh:
        next(fh)
        for line in fh:
            refs.add(line.split("\t")[1])
    with open(os.path.join(args.out, "rep_ref_assemblies.txt"), "w") as fh:
        for a in sorted(refs): fh.write(a + "\n")
    print(f"[cluster] {len(refs)} distinct ref assemblies host the reps", flush=True)

    # Cluster summary
    summary = os.path.join(args.out, "cluster_summary.tsv")
    with open(summary, "w") as fh:
        fh.write("cluster_id\trep_id\tn_members\tn_distinct_asms\tmember_lengths\n")
        for cl in clusters:
            asms, lens = set(), []
            for m in cl["members"]:
                if m in by_refid:
                    asms.add(by_refid[m]["assembly"])
                    lens.append(by_refid[m]["length"])
            lens_str = ",".join(str(x) for x in lens[:5]) + ("..." if len(lens)>5 else "")
            fh.write(f"{cl['id']}\t{cl['rep']}\t{len(cl['members'])}\t"
                     f"{len(asms)}\t{lens_str}\n")
    print(f"[cluster] cluster_summary.tsv written.", flush=True)


if __name__ == "__main__":
    main()
