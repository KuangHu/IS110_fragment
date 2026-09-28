#!/usr/bin/env python3
"""Parse per-target PAF files from cross-organism large-IS alignment and
produce a filtered HGT catalogue.

Each PAF row has query = "SOURCE_SP|ref_id" (added prefix so we can filter
by source), target = "ASM|CONTIG" (from the species DB).

For each hit:
  - source species = query prefix
  - target species = the species whose mmi this PAF was searched against
  - drop hits where source == target (self-species)
  - filter by identity + coverage
  - annotate with taxonomic-distance tier
"""
import argparse, csv, glob, os, sys
from collections import defaultdict


# Taxonomic groups (family-level)
TAXO = {
    "klebsiella_pneumoniae":    ("Enterobacteriaceae", "Gammaproteobacteria", "Proteobacteria"),
    "escherichia_coli":         ("Enterobacteriaceae", "Gammaproteobacteria", "Proteobacteria"),
    "salmonella_enterica":      ("Enterobacteriaceae", "Gammaproteobacteria", "Proteobacteria"),
    "enterobacter_hormaechei":  ("Enterobacteriaceae", "Gammaproteobacteria", "Proteobacteria"),
    "yersinia_enterocolitica":  ("Yersiniaceae",       "Gammaproteobacteria", "Proteobacteria"),
    "burkholderia_pseudomallei":("Burkholderiaceae",   "Betaproteobacteria",  "Proteobacteria"),
    "enterococcus_faecium":     ("Enterococcaceae",    "Bacilli",             "Firmicutes"),
    "streptococcus_suis":       ("Streptococcaceae",   "Bacilli",             "Firmicutes"),
    "mycobacterium_tuberculosis":("Mycobacteriaceae",  "Actinomycetia",       "Actinomycetota"),
    "leptospira_interrogans":   ("Leptospiraceae",     "Spirochaetia",        "Spirochaetota"),
}


def distance_tier(sp1, sp2):
    """Return taxonomic-distance category between two species."""
    if sp1 == sp2: return "same_species"
    t1 = TAXO.get(sp1); t2 = TAXO.get(sp2)
    if not t1 or not t2: return "unknown"
    if t1[0] == t2[0]: return "same_family"
    if t1[1] == t2[1]: return "same_class"
    if t1[2] == t2[2]: return "same_phylum"
    return "different_phylum"


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--paf-dir", required=True,
                   help="dir containing hits_vs_<target>.paf")
    p.add_argument("--out", required=True)
    p.add_argument("--min-identity", type=float, default=95.0)
    p.add_argument("--min-cov",      type=float, default=0.50,
                   help="min fraction of query covered by the hit")
    p.add_argument("--strict-identity", type=float, default=99.0)
    p.add_argument("--strict-cov",   type=float, default=0.90)
    return p.parse_args()


def main():
    args = parse_args()

    hits = []
    for paf in sorted(glob.glob(f"{args.paf_dir}/hits_vs_*.paf")):
        target_sp = os.path.basename(paf).replace("hits_vs_", "").replace(".paf", "")
        with open(paf) as fh:
            for line in fh:
                c = line.rstrip().split("\t")
                if len(c) < 12: continue
                qname = c[0]                # "SOURCE_SP|ref_id"
                if "|" not in qname: continue
                source_sp, ref_id = qname.split("|", 1)
                # Skip self-species hits
                if source_sp == target_sp: continue
                qlen = int(c[1]); qs = int(c[2]); qe = int(c[3])
                tname = c[5]                # "ASM|CONTIG"
                if "|" not in tname: continue
                tgt_asm, tgt_contig = tname.split("|", 1)
                # Skip cases where query source matches target assembly (shouldn't happen but safety)
                strand = c[4]
                ts = int(c[7]); te = int(c[8])
                matches = int(c[9]); block = int(c[10])
                ident = matches / block * 100 if block > 0 else 0
                cov = (qe - qs) / qlen if qlen > 0 else 0
                if ident < args.min_identity: continue
                if cov < args.min_cov: continue
                is_strict = (ident >= args.strict_identity and cov >= args.strict_cov)
                hits.append({
                    "source_species": source_sp,
                    "ref_id": ref_id,
                    "query_bp": qlen,
                    "target_species": target_sp,
                    "target_assembly": tgt_asm,
                    "target_contig": tgt_contig,
                    "strand": strand,
                    "qs": qs, "qe": qe, "ts": ts, "te": te,
                    "cov": cov, "identity": ident,
                    "strict": is_strict,
                    "distance": distance_tier(source_sp, target_sp),
                })

    # Write catalogue
    with open(args.out, "w") as fo:
        fo.write("source_species\tref_id\tquery_bp\ttarget_species\ttarget_assembly\t"
                 "target_contig\tstrand\tqs\tqe\tts\tte\tcov\tidentity\tstrict\tdistance\n")
        for h in hits:
            fo.write("\t".join([
                h["source_species"], h["ref_id"], str(h["query_bp"]),
                h["target_species"], h["target_assembly"], h["target_contig"],
                h["strand"], str(h["qs"]), str(h["qe"]),
                str(h["ts"]), str(h["te"]),
                f"{h['cov']:.3f}", f"{h['identity']:.1f}",
                str(h["strict"]), h["distance"],
            ]) + "\n")
    print(f"Total cross-species hits (id>={args.min_identity}%, cov>={args.min_cov}): {len(hits)}")
    print(f"Wrote {args.out}")

    # Summarise by distance tier
    from collections import Counter
    dist_all = Counter(h["distance"] for h in hits)
    dist_strict = Counter(h["distance"] for h in hits if h["strict"])
    print(f"\nBy distance tier:")
    for tier in ["same_family", "same_class", "same_phylum", "different_phylum"]:
        print(f"  {tier:<20}: all={dist_all.get(tier, 0):>6}  strict={dist_strict.get(tier, 0):>5}")

    # Cross-species pair matrix
    pairs = Counter()
    strict_pairs = Counter()
    for h in hits:
        pairs[(h["source_species"], h["target_species"])] += 1
        if h["strict"]:
            strict_pairs[(h["source_species"], h["target_species"])] += 1

    print(f"\nTop species pairs by strict-HGT hit count:")
    for (sp1, sp2), n in sorted(strict_pairs.items(), key=lambda x: -x[1])[:20]:
        t = distance_tier(sp1, sp2)
        print(f"  {sp1:<24} -> {sp2:<24} {n:>4} ({t})")


if __name__ == "__main__":
    main()
