#!/usr/bin/env python3
"""Find genome rearrangements from the anchor PAF.

For each (ref_id, target_assembly), classify anchor pair patterns:
  - Same contig, same strand    → normal layout (already analyzed by stages 4–6)
  - Same contig, OPPOSITE strands → INVERSION at that locus
  - Different contigs             → TRANSLOCATION or contig break
  - Anchor matches at multiple loci → DUPLICATION

Reads the PAF produced by Stage 3 (anchors vs DB) and the ref-table from
Stage 2 — no IS-family-specific logic.

Usage:
    find_rearrangements.py --paf anchors_vs_db.paf --out rearrangements
"""
import argparse
import os
import re
import sys
from collections import defaultdict


ANCHOR_RE = re.compile(r"^(.+)__(up|down)(\d+)$")


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--paf", required=True, help="Anchors vs DB PAF (Stage 3 output)")
    p.add_argument("--out", required=True, help="Output prefix")
    p.add_argument("--min-identity", type=float, default=95)
    p.add_argument("--min-coverage", type=float, default=90)
    p.add_argument("--anchor-D", type=int, default=1000,
                   help="Use this anchor distance for classification (default 1000)")
    p.add_argument("--max-examples", type=int, default=50)
    p.add_argument("--fai", default="",
                   help="Optional .fai of the species DB — enables anchor-depth + "
                        "target-contig-length annotation for translocation candidates "
                        "(needed for the chromosome-level target filter).")
    return p.parse_args()


def main():
    args = parse_args()

    obs = defaultdict(lambda: {"up": defaultdict(list), "down": defaultdict(list)})

    n_total = 0
    n_kept = 0
    with open(args.paf) as f:
        for line in f:
            n_total += 1
            if n_total % 10_000_000 == 0:
                print(f"  read {n_total:,} lines", flush=True, file=sys.stderr)
            c = line.split("\t")
            if len(c) < 12:
                continue
            qname = c[0]
            qlen = int(c[1])
            qs, qe = int(c[2]), int(c[3])
            strand = c[4]
            tname = c[5]
            ts, te = int(c[7]), int(c[8])
            matches, block = int(c[9]), int(c[10])
            ident = matches / block * 100 if block > 0 else 0
            cov = (qe - qs) / qlen * 100
            if ident < args.min_identity or cov < args.min_coverage:
                continue
            m = ANCHOR_RE.match(qname)
            if not m:
                continue
            ref_id, side, D = m.group(1), m.group(2), int(m.group(3))
            if D != args.anchor_D:
                continue
            n_kept += 1

            if "|" in tname:
                assembly, contig = tname.split("|", 1)
            else:
                assembly, contig = "unknown", tname

            obs[(ref_id, assembly)][side][contig].append({
                "strand": strand, "ts": ts, "te": te, "ident": ident, "contig": contig,
            })

    print(f"PAF: {n_total:,} lines, {n_kept:,} kept (D={args.anchor_D})", file=sys.stderr)
    print(f"Unique (ref, target_assembly) pairs: {len(obs):,}", file=sys.stderr)

    counters = defaultdict(int)
    examples = defaultdict(list)
    inversion_all = []   # every inversion_only candidate (no cap) + estimated span
    translocation_all = []   # every translocation_only candidate + anchor depths + contig lens
    duplication_all = []   # every has_duplication candidate + (up_mult, down_mult, copy_span_estimate)

    # Load target-contig lengths once (so we can score anchor depth for translocations)
    contig_lens = {}
    if getattr(args, "fai", None) and os.path.exists(args.fai):
        print(f"Loading contig lengths from {args.fai}...", file=sys.stderr)
        with open(args.fai) as fh:
            for line in fh:
                p = line.split("\t")
                if len(p) >= 2:
                    contig_lens[p[0]] = int(p[1])
        print(f"  {len(contig_lens):,} contigs loaded", file=sys.stderr)

    for (ref_id, assembly), sides in obs.items():
        # Dedup hits per (contig, ts, te, strand): same anchor hitting the same
        # position multiple times = minimap2 secondary alignments OR Logan-derived
        # contig redundancy in the DB, not a real biological duplication.
        for side in ("up", "down"):
            for contig, hs in sides[side].items():
                seen = set()
                uniq = []
                for h in hs:
                    key = (h["ts"], h["te"], h["strand"])
                    if key in seen:
                        continue
                    seen.add(key)
                    uniq.append(h)
                sides[side][contig] = uniq
        up_hits = sides["up"]
        down_hits = sides["down"]
        if not up_hits or not down_hits:
            continue

        has_normal = False
        has_inversion = False
        has_translocation = False

        for contig in set(up_hits.keys()) & set(down_hits.keys()):
            for u in up_hits[contig]:
                for d in down_hits[contig]:
                    if u["strand"] == d["strand"]:
                        has_normal = True
                    else:
                        has_inversion = True

        all_up_contigs = set(up_hits.keys())
        all_down_contigs = set(down_hits.keys())
        if all_up_contigs and all_down_contigs:
            if all_up_contigs.isdisjoint(all_down_contigs):
                has_translocation = True

        if has_normal:
            counters["normal"] += 1
        if has_inversion and not has_normal:
            counters["inversion_only"] += 1
            if len(examples["inversion_only"]) < args.max_examples:
                examples["inversion_only"].append((ref_id, assembly, sides))
            # ALSO record EVERY inversion candidate (no cap) with its estimated
            # inverted-segment span. For each opposite-strand u/d pair on the same
            # target contig, the inverted region is at least the outer envelope
            # between the two anchor positions. Keep the LARGEST estimate per event.
            best_span = 0
            best_pair = None
            for contig in set(up_hits.keys()) & set(down_hits.keys()):
                for u in up_hits[contig]:
                    for d in down_hits[contig]:
                        if u["strand"] == d["strand"]:
                            continue
                        span = max(u["te"], d["te"]) - min(u["ts"], d["ts"])
                        if span > best_span:
                            best_span = span
                            best_pair = (contig, u, d)
            inversion_all.append((ref_id, assembly, sides, best_span, best_pair))
        elif has_inversion:
            counters["inversion_plus_normal"] += 1
        if has_translocation and not has_normal and not has_inversion:
            counters["translocation_only"] += 1
            if len(examples["translocation_only"]) < args.max_examples:
                examples["translocation_only"].append((ref_id, assembly, sides))
            # ALSO record every translocation candidate (no cap) with the BEST
            # opposite-contig anchor pair: pick the up/down hit whose minimum
            # anchor-to-contig-end depth is largest (= deepest into contigs,
            # most likely a real translocation rather than a contig-break artifact).
            best_depth = -1
            best_pair = None
            # .fai is keyed by full "<assembly>|<contig>"; the up/down dicts are
            # keyed by stripped contig name. Look up the full key, fall back to bare.
            def _len(contig):
                full = f"{assembly}|{contig}"
                return contig_lens.get(full) or contig_lens.get(contig) or 0
            for u_contig, u_hits in up_hits.items():
                for d_contig, d_hits in down_hits.items():
                    if u_contig == d_contig:
                        continue   # same-contig: not a translocation signature
                    u_len = _len(u_contig)
                    d_len = _len(d_contig)
                    for u in u_hits:
                        u_depth = min(u["ts"], max(0, u_len - u["te"])) if u_len else 0
                        for d in d_hits:
                            d_depth = min(d["ts"], max(0, d_len - d["te"])) if d_len else 0
                            min_depth = min(u_depth, d_depth)
                            if min_depth > best_depth:
                                best_depth = min_depth
                                best_pair = (u_contig, u, u_len, d_contig, d, d_len)
            translocation_all.append((ref_id, assembly, sides, best_depth, best_pair))
        elif has_translocation:
            counters["translocation_plus"] += 1

        n_dup = sum(1 for c, hs in up_hits.items() if len(hs) > 1)
        n_dup += sum(1 for c, hs in down_hits.items() if len(hs) > 1)
        if n_dup > 0:
            counters["has_duplication"] += 1
            if len(examples["duplication"]) < args.max_examples:
                examples["duplication"].append((ref_id, assembly, sides))
            # ALSO record EVERY duplication candidate. For each contig with >=2 up
            # OR >=2 down anchor hits, pick the LARGEST same-strand pair distance
            # as the candidate inter-copy spacing (proxy for duplicated segment size).
            up_mult = max((len(hs) for hs in up_hits.values()), default=0)
            down_mult = max((len(hs) for hs in down_hits.values()), default=0)
            best_span = 0
            best_pair = None
            best_side = None
            best_contig = None
            for side, hits in (("up", up_hits), ("down", down_hits)):
                for contig, hs in hits.items():
                    if len(hs) < 2:
                        continue
                    # pairwise distances between same-strand hits
                    same_strand = {"+": [], "-": []}
                    for h in hs:
                        same_strand[h["strand"]].append(h)
                    for st in ("+", "-"):
                        L = same_strand[st]
                        for i in range(len(L)):
                            for j in range(i + 1, len(L)):
                                span = abs(L[j]["ts"] - L[i]["ts"])
                                if span > best_span:
                                    best_span = span
                                    best_pair = (L[i], L[j])
                                    best_side = side
                                    best_contig = contig
            duplication_all.append((ref_id, assembly, sides, up_mult, down_mult,
                                    best_span, best_pair, best_side, best_contig))

    print(f"\n=== Rearrangement signature counts (per ref-assembly pair) ===")
    for k, v in sorted(counters.items(), key=lambda x: -x[1]):
        print(f"  {k}: {v:,}")

    out_path = f"{args.out}_examples.tsv"
    with open(out_path, "w") as f:
        f.write("category\tref_id\tassembly\tup_contigs\tdown_contigs\tdetails\n")
        for cat, exs in examples.items():
            for ref_id, assembly, sides in exs:
                up_contigs = ",".join(sides["up"].keys())
                down_contigs = ",".join(sides["down"].keys())
                up_detail = []
                for contig, hs in sides["up"].items():
                    for h in hs[:2]:
                        up_detail.append(f"{contig}:{h['ts']}-{h['te']}({h['strand']})")
                down_detail = []
                for contig, hs in sides["down"].items():
                    for h in hs[:2]:
                        down_detail.append(f"{contig}:{h['ts']}-{h['te']}({h['strand']})")
                detail = f"up=[{';'.join(up_detail)}] down=[{';'.join(down_detail)}]"
                f.write(f"{cat}\t{ref_id}\t{assembly}\t{up_contigs}\t{down_contigs}\t{detail}\n")
    print(f"\nWrote examples: {out_path}")

    # ALL inversion candidates with estimated span — uncapped, in the same format
    # as examples.tsv so the validator can consume it directly.
    all_inv_path = f"{args.out}_inversion_all.tsv"
    with open(all_inv_path, "w") as f:
        f.write("category\tref_id\tassembly\tup_contigs\tdown_contigs\tdetails\t"
                "estimated_inverted_bp\n")
        for ref_id, assembly, sides, span, best_pair in inversion_all:
            up_contigs = ",".join(sides["up"].keys())
            down_contigs = ",".join(sides["down"].keys())
            # write the BEST opposite-strand u/d pair as the details (so the
            # validator extracts the right target region)
            if best_pair is not None:
                contig, u, d = best_pair
                up_detail = f"{contig}:{u['ts']}-{u['te']}({u['strand']})"
                down_detail = f"{contig}:{d['ts']}-{d['te']}({d['strand']})"
            else:
                up_detail = down_detail = ""
            detail = f"up=[{up_detail}] down=[{down_detail}]"
            f.write(f"inversion_only\t{ref_id}\t{assembly}\t{up_contigs}\t"
                    f"{down_contigs}\t{detail}\t{span}\n")
    print(f"Wrote all-inversion catalogue: {all_inv_path}  "
          f"({len(inversion_all):,} events)")

    # ALL translocation candidates with anchor-depth + target-contig-length annotation
    all_tr_path = f"{args.out}_translocation_all.tsv"
    with open(all_tr_path, "w") as f:
        f.write("category\tref_id\tassembly\tup_contigs\tdown_contigs\tdetails\t"
                "min_anchor_depth\tup_contig_len\tdown_contig_len\n")
        for ref_id, assembly, sides, depth, best_pair in translocation_all:
            up_contigs = ",".join(sides["up"].keys())
            down_contigs = ",".join(sides["down"].keys())
            if best_pair is not None:
                u_contig, u, u_len, d_contig, d, d_len = best_pair
                up_detail = f"{u_contig}:{u['ts']}-{u['te']}({u['strand']})"
                down_detail = f"{d_contig}:{d['ts']}-{d['te']}({d['strand']})"
            else:
                u_len = d_len = 0
                up_detail = down_detail = ""
            detail = f"up=[{up_detail}] down=[{down_detail}]"
            f.write(f"translocation_only\t{ref_id}\t{assembly}\t{up_contigs}\t"
                    f"{down_contigs}\t{detail}\t{depth}\t{u_len}\t{d_len}\n")
    print(f"Wrote all-translocation catalogue: {all_tr_path}  "
          f"({len(translocation_all):,} events)")

    # ALL duplication candidates — anchor multiplicity tells us the same anchor
    # hit two+ positions on the same target contig, the classic IS-HR signature.
    # Write ALL anchor hit positions (not just one pair) so the strict-check
    # downstream can extract the full target region.
    all_dup_path = f"{args.out}_duplication_all.tsv"
    with open(all_dup_path, "w") as f:
        f.write("category\tref_id\tassembly\tup_contigs\tdown_contigs\tdetails\t"
                "up_anchor_mult\tdown_anchor_mult\tcopy_span_estimate\tbest_side\n")
        for (ref_id, assembly, sides, up_mult, down_mult,
             span, best_pair, best_side, best_contig) in duplication_all:
            up_contigs = ",".join(sides["up"].keys())
            down_contigs = ",".join(sides["down"].keys())
            # Determine which side carries the multiplicity (best_side) — if
            # tied or ambiguous, prefer the side with more total hits.
            side_for_detail = best_side
            if side_for_detail is None:
                if up_mult >= 2:   side_for_detail = "up"
                elif down_mult >= 2: side_for_detail = "down"
            if side_for_detail == "up":
                hits_dict = sides["up"]
            elif side_for_detail == "down":
                hits_dict = sides["down"]
            else:
                hits_dict = {}
            # Only contigs with multiplicity (>=2 distinct hits on the same contig)
            # carry the real duplication signal. Cross-contig hits at the same
            # position are assembly-redundancy artifacts.
            dup_contigs = {c: hs for c, hs in hits_dict.items() if len(hs) >= 2}
            all_hits = []
            for c, hs in dup_contigs.items():
                for h in hs:
                    all_hits.append(f"{c}:{h['ts']}-{h['te']}({h['strand']})")
            uniq = list(dict.fromkeys(all_hits))
            pair_detail = ";".join(uniq)
            if side_for_detail == "up":
                detail = f"up=[{pair_detail}] down=[]"
            elif side_for_detail == "down":
                detail = f"up=[] down=[{pair_detail}]"
            else:
                detail = "up=[] down=[]"
            f.write(f"duplication\t{ref_id}\t{assembly}\t{up_contigs}\t"
                    f"{down_contigs}\t{detail}\t{up_mult}\t{down_mult}\t"
                    f"{span}\t{side_for_detail or ''}\n")
    print(f"Wrote all-duplication catalogue: {all_dup_path}  "
          f"({len(duplication_all):,} events)")


if __name__ == "__main__":
    main()
