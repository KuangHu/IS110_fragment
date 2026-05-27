#!/usr/bin/env python3
"""Find "downstream" IS110 variants at the same anchor site as V1.

Requirements:
  1. Same anchor: V1's 5 kb upstream + 5 kb downstream flanks both align in target
     at close together positions (so target has matching chromosomal site)
  2. Many 95% identity regions: the between-anchor sequence in target shares
     many blocks (>=95% identity) with V1's IS110 sequence

For each V1, find:
  - empty: between-anchor in target is ~0 bp
  - V1 (clonal): between-anchor matches V1 ~fully
  - downstream: between-anchor has V1's content + extra (>=95% blocks match V1)
  - unrelated: between-anchor doesn't match V1 well

Iterative: a downstream becomes the next-round "V1" query.

Usage:
  downstream_finder.py --records V1_records.json --db target.fa --out OUT_DIR
"""
import argparse, json, os, re, subprocess, sys, hashlib
from collections import defaultdict

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--records", required=True,
                   help="JSON file with V1 records (records_final/records.json format)")
    p.add_argument("--db", required=True,
                   help="Target database FASTA (indexed, has .fai)")
    p.add_argument("--out", required=True, help="Output directory")
    p.add_argument("--threads", type=int, default=32)
    p.add_argument("--anchor-len", type=int, default=5000,
                   help="Length of upstream/downstream anchor (default 5000)")
    p.add_argument("--anchor-identity", type=float, default=95,
                   help="Min identity for anchor hits")
    p.add_argument("--anchor-cov", type=float, default=80,
                   help="Min query coverage for anchor hits")
    p.add_argument("--max-pair-dist", type=int, default=500000,
                   help="Max distance between paired anchor hits (bp)")
    p.add_argument("--min-similarity-blocks", type=int, default=2,
                   help="Min # of 95%% identity blocks needed for 'downstream' call")
    p.add_argument("--min-extra-bp", type=int, default=200,
                   help="Min extra bp in target (vs V1) to call 'downstream' (vs clonal V1)")
    p.add_argument("--min-v1-coverage", type=float, default=80,
                   help="Min %% of V1 sequence that must be covered in insertion variant (default 80%%)")
    p.add_argument("--min-deletion-match", type=int, default=200,
                   help="Min bp of V1 that must be present in a 'deletion' variant (default 200)")
    p.add_argument("--min-deletion-shortening", type=int, default=200,
                   help="Min bp shorter than V1 to classify as deletion (default 200)")
    p.add_argument("--max-records", type=int, default=0,
                   help="Limit input V1 records for testing (0 = all)")
    return p.parse_args()


def revcomp(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def extract_target(db, contig, start, end):
    """samtools faidx region."""
    region = f"{contig}:{max(1, start)}-{end}"
    res = subprocess.run(["samtools", "faidx", db, region],
                         capture_output=True, text=True)
    if res.returncode != 0: return None
    return "".join(res.stdout.split("\n")[1:]).strip().upper()


def extract_target_batch(db, regions, batch_size=1000):
    """Batch samtools faidx - takes many regions at once.

    regions: list of (contig, start, end) tuples.
    Returns dict: (contig, start, end) -> sequence string.
    """
    results = {}
    for i in range(0, len(regions), batch_size):
        batch = regions[i:i+batch_size]
        region_strs = [f"{c}:{max(1,s)}-{e}" for c, s, e in batch]
        res = subprocess.run(["samtools", "faidx", db] + region_strs,
                             capture_output=True, text=True)
        if res.returncode != 0:
            continue
        # Parse multi-FASTA output
        cur_header, cur_seq = None, []
        for line in res.stdout.split("\n"):
            if line.startswith(">"):
                if cur_header:
                    # Parse "CONTIG:S-E" back to tuple
                    m = re.match(r"(.+):(\d+)-(\d+)", cur_header)
                    if m:
                        key = (m.group(1), int(m.group(2)), int(m.group(3)))
                        # Match to original key (start may be max(1,start))
                        for c, s, e in batch:
                            if c == key[0] and max(1, s) == key[1] and e == key[2]:
                                results[(c, s, e)] = "".join(cur_seq).upper()
                                break
                cur_header = line[1:].strip()
                cur_seq = []
            else:
                cur_seq.append(line.strip())
        if cur_header:
            m = re.match(r"(.+):(\d+)-(\d+)", cur_header)
            if m:
                key = (m.group(1), int(m.group(2)), int(m.group(3)))
                for c, s, e in batch:
                    if c == key[0] and max(1, s) == key[1] and e == key[2]:
                        results[(c, s, e)] = "".join(cur_seq).upper()
                        break
    return results


def write_fasta(items, path):
    with open(path, "w") as f:
        for header, seq in items:
            f.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")


def short_hash(s):
    return hashlib.md5(s.encode()).hexdigest()[:10]


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)

    print(f"Loading V1 records from {args.records} ...", file=sys.stderr, flush=True)
    with open(args.records) as f:
        records = json.load(f)
    if args.max_records > 0:
        records = records[:args.max_records]
    print(f"  {len(records):,} V1 records", file=sys.stderr)

    # === Step 1: Build anchor FASTA ===
    # For each V1, write 2 anchor seqs (up + down)
    anchor_fa = f"{args.out}/anchors.fa"
    is110_fa = f"{args.out}/v1_is110.fa"

    anchors_list = []
    is110_list = []
    ref_meta = {}  # v1_id -> metadata
    for r in records:
        # Field name compatibility: records_final uses "is110_id"; build_records.py
        # (Cross_reference_IS Stage 6) uses "ref_id". Accept either.
        rid = r.get("is110_id") or r.get("ref_id")
        if not rid:
            continue
        # Layout compatibility:
        #   records_final: source.{upstream_flank, downstream_flank, is_element}
        #   build_records.py: top-level {upstream_flank, downstream_flank, is_element}
        src = r.get("source", {})
        up_obj = (src.get("upstream_flank") or r.get("upstream_flank") or {})
        down_obj = (src.get("downstream_flank") or r.get("downstream_flank") or {})
        is_obj = (src.get("is_element") or r.get("is_element") or {})
        up = up_obj.get("sequence", "")
        down = down_obj.get("sequence", "")
        is110_seq = is_obj.get("sequence", "")
        if not (up and down and is110_seq):
            continue
        if len(up) < args.anchor_len or len(down) < args.anchor_len:
            continue
        # Anchor = innermost 5 kb (closest to IS)
        up_anchor = up[-args.anchor_len:]      # 5 kb immediately upstream of IS
        down_anchor = down[:args.anchor_len]   # 5 kb immediately downstream
        anchors_list.append((f"{rid}__up", up_anchor))
        anchors_list.append((f"{rid}__down", down_anchor))
        is110_list.append((rid, is110_seq))
        ref_meta[rid] = {
            "is110_len": is_obj.get("length") or len(is110_seq),
            "up_anchor_seq": up_anchor,
            "down_anchor_seq": down_anchor,
            "is110_seq": is110_seq,
            "source": {k: v for k, v in src.items()
                       if k in ("assembly", "contig")},
        }
    print(f"  Built {len(anchors_list)//2:,} anchor pairs", file=sys.stderr)
    write_fasta(anchors_list, anchor_fa)
    write_fasta(is110_list, is110_fa)

    # === Step 2: minimap2 anchors against DB ===
    anchor_paf = f"{args.out}/anchors_vs_db.paf"
    if not os.path.exists(anchor_paf) or os.path.getsize(anchor_paf) == 0:
        print(f"Running minimap2 anchors vs DB ...", file=sys.stderr, flush=True)
        subprocess.run([
            "minimap2", "-x", "asm10", "-c", "--eqx",
            "--secondary=yes", "-N", "20", "-p", "0.5",
            "-t", str(args.threads), args.db, anchor_fa, "-o", anchor_paf
        ], check=True)
    print(f"  anchor PAF: {anchor_paf} ({os.path.getsize(anchor_paf)/1e6:.1f} MB)",
          file=sys.stderr)

    # === Step 3: Parse anchor hits and find pairs ===
    print(f"Parsing anchor hits ...", file=sys.stderr, flush=True)
    # anchor_id -> [(target, target_start, target_end, strand, ident, cov)]
    anchor_hits = defaultdict(list)
    with open(anchor_paf) as f:
        for line in f:
            c = line.split("\t")
            if len(c) < 12: continue
            qname, qlen = c[0], int(c[1])
            qs, qe = int(c[2]), int(c[3])
            strand = c[4]
            tname = c[5]
            ts, te = int(c[7]), int(c[8])
            matches, block = int(c[9]), int(c[10])
            ident = matches / block * 100 if block > 0 else 0
            cov = (qe - qs) / qlen * 100
            if ident < args.anchor_identity or cov < args.anchor_cov: continue
            anchor_hits[qname].append({
                "target": tname, "ts": ts, "te": te, "strand": strand,
                "ident": ident, "cov": cov,
            })

    # === Step 4: For each V1, find matching anchor pairs (up + down on same contig) ===
    print(f"Pairing anchors per V1 ...", file=sys.stderr, flush=True)
    v1_pairs = defaultdict(list)  # v1_id -> list of {target, strand, between_start, between_end}
    for rid in ref_meta:
        ups = anchor_hits.get(f"{rid}__up", [])
        downs = anchor_hits.get(f"{rid}__down", [])
        if not ups or not downs: continue
        # Group by target
        ups_by_t = defaultdict(list)
        downs_by_t = defaultdict(list)
        for h in ups: ups_by_t[h["target"]].append(h)
        for h in downs: downs_by_t[h["target"]].append(h)
        for tname in set(ups_by_t.keys()) & set(downs_by_t.keys()):
            for u in ups_by_t[tname]:
                for d in downs_by_t[tname]:
                    if u["strand"] != d["strand"]: continue  # skip inversions
                    # Determine between-anchor region
                    if u["strand"] == "+":
                        # up anchor ends at u["te"], down anchor starts at d["ts"]
                        between_start = u["te"]
                        between_end = d["ts"]
                    else:
                        between_start = d["te"]
                        between_end = u["ts"]
                    if between_end <= between_start: continue
                    dist = between_end - between_start
                    if dist > args.max_pair_dist: continue
                    v1_pairs[rid].append({
                        "target": tname, "strand": u["strand"],
                        "between_start": between_start, "between_end": between_end,
                        "between_len": dist,
                        "up_ident": u["ident"], "down_ident": d["ident"],
                    })

    n_pairs = sum(len(v) for v in v1_pairs.values())
    print(f"  Found {n_pairs:,} anchor pairs across {len(v1_pairs):,} V1s",
          file=sys.stderr)

    # === Step 5: Extract between-anchor sequences (batched) ===
    between_fa = f"{args.out}/between_anchors.fa"
    print(f"Extracting between-anchor sequences (batched) ...",
          file=sys.stderr, flush=True)

    # Build region list
    all_regions = []  # (contig, start, end)
    region_to_pair = {}  # (contig, start, end) -> (rid, i, p)
    for rid, pairs in v1_pairs.items():
        for i, p in enumerate(pairs):
            region = (p["target"], p["between_start"] + 1, p["between_end"])
            all_regions.append(region)
            region_to_pair[region] = (rid, i, p)
    print(f"  Total regions to extract: {len(all_regions):,}",
          file=sys.stderr, flush=True)

    # Batch fetch
    seqs = extract_target_batch(args.db, all_regions, batch_size=500)
    print(f"  Fetched {len(seqs):,} sequences", file=sys.stderr, flush=True)

    extracted = []
    for region, seq in seqs.items():
        if region not in region_to_pair: continue
        rid, i, p = region_to_pair[region]
        if p["strand"] == "-": seq = revcomp(seq)
        p["between_seq"] = seq
        p["pair_id"] = f"{rid}__pair{i}__{p['target']}_{p['between_start']}_{p['between_end']}"
        extracted.append((p["pair_id"], seq))
    write_fasta(extracted, between_fa)
    print(f"  Wrote {len(extracted):,} between-anchor sequences",
          file=sys.stderr)

    # === Step 6: Compare each between-seq to its V1 IS110 using minimap2 ===
    print(f"Aligning between-seqs to V1 IS110 sequences ...", file=sys.stderr, flush=True)
    compare_paf = f"{args.out}/between_vs_v1.paf"
    # Build query (each between seq) vs target (each V1 IS110)
    # Use minimap2 with high sensitivity
    subprocess.run([
        "minimap2", "-x", "asm10", "-c", "--eqx",
        "-t", str(args.threads), is110_fa, between_fa, "-o", compare_paf
    ], check=True, capture_output=True)
    print(f"  comparison PAF: {os.path.getsize(compare_paf)/1e6:.1f} MB",
          file=sys.stderr)

    # Parse: for each between seq, list of alignments to V1
    between_aligns = defaultdict(list)
    with open(compare_paf) as f:
        for line in f:
            c = line.split("\t")
            if len(c) < 12: continue
            qname, qlen = c[0], int(c[1])
            qs, qe = int(c[2]), int(c[3])
            strand = c[4]
            tname = c[5]
            ts, te = int(c[7]), int(c[8])
            matches, block = int(c[9]), int(c[10])
            ident = matches / block * 100 if block > 0 else 0
            if ident < 95: continue
            between_aligns[qname].append({
                "v1_id": tname, "v1_start": ts, "v1_end": te,
                "between_start": qs, "between_end": qe,
                "strand": strand, "ident": ident,
                "block_len": qe - qs,
            })

    # === Step 7: Classify each observation ===
    print(f"Classifying observations ...", file=sys.stderr, flush=True)
    obs_records = []
    counts = defaultdict(int)
    for rid, pairs in v1_pairs.items():
        v1_len = ref_meta[rid]["is110_len"]
        for p in pairs:
            if "between_seq" not in p: continue
            blocks = between_aligns.get(p["pair_id"], [])
            # Restrict to alignments to OUR V1 (rid)
            blocks = [b for b in blocks if b["v1_id"] == rid]

            between_len = p["between_len"]
            total_matched = sum(b["block_len"] for b in blocks)

            v1_cov_frac = total_matched / v1_len if v1_len > 0 else 0
            min_cov_frac = args.min_v1_coverage / 100.0

            # Classification:
            #  empty     : nothing between anchors (< 100 bp)
            #  deletion  : shorter than V1 but contains V1 fragments (>= min_deletion_match bp)
            #  clonal_V1 : matches V1 fully (~same length, >=95% coverage)
            #  insertion : V1's content + extra DNA (>=80% V1 coverage, >=min_extra_bp larger)
            #  unrelated : virtually no V1 content (<10% coverage)
            #  partial   : everything else (ambiguous)
            if between_len < 100:
                category = "empty"
            elif (between_len < v1_len - args.min_deletion_shortening
                  and total_matched >= args.min_deletion_match):
                category = "deletion"
            elif v1_cov_frac >= 0.95 and between_len < v1_len + args.min_extra_bp:
                category = "clonal_V1"
            elif (len(blocks) >= args.min_similarity_blocks
                  and v1_cov_frac >= min_cov_frac
                  and between_len >= v1_len + args.min_extra_bp):
                category = "insertion"
            elif v1_cov_frac < 0.1:
                category = "unrelated"
            else:
                category = "partial"

            counts[category] += 1

            obs = {
                "v1_parent_id": rid,
                "category": category,
                "target": p["target"],
                "target_strand": p["strand"],
                "anchor_site": {
                    "between_start": p["between_start"],
                    "between_end": p["between_end"],
                    "between_len": between_len,
                    "up_anchor_identity": round(p["up_ident"], 1),
                    "down_anchor_identity": round(p["down_ident"], 1),
                },
                "comparison_to_v1": {
                    "v1_len": v1_len,
                    "n_match_blocks": len(blocks),
                    "total_matched_bp": total_matched,
                    "blocks": [{"v1_pos": [b["v1_start"], b["v1_end"]],
                                "between_pos": [b["between_start"], b["between_end"]],
                                "ident": round(b["ident"], 1)} for b in blocks],
                },
            }
            if category == "insertion":
                # Save the full insertion variant sequence + ID
                obs["variant_id"] = f"ins_{short_hash(p['between_seq'])}"
                obs["variant_seq"] = p["between_seq"]
                obs["variant_len"] = between_len
            elif category == "deletion":
                # Save the deletion variant (shorter than V1) + ID
                obs["variant_id"] = f"del_{short_hash(p['between_seq'])}"
                obs["variant_seq"] = p["between_seq"]
                obs["variant_len"] = between_len
            obs_records.append(obs)

    print(f"\n=== Classification counts ===", file=sys.stderr)
    for k in ("empty", "deletion", "clonal_V1", "insertion", "partial", "unrelated"):
        print(f"  {k:>12}: {counts[k]:>8,}", file=sys.stderr)

    # === Step 8: Save all observations + extract variant cluster reps ===
    with open(f"{args.out}/observations.json", "w") as f:
        json.dump(obs_records, f, indent=2)

    # Collect unique insertion + deletion variants (separate FASTAs)
    seen_ins = {}
    seen_del = {}
    for obs in obs_records:
        if obs["category"] == "insertion":
            d = seen_ins
        elif obs["category"] == "deletion":
            d = seen_del
        else:
            continue
        vid = obs["variant_id"]
        if vid not in d:
            d[vid] = obs["variant_seq"]

    if seen_ins:
        next_query = f"{args.out}/insertion_variants.fa"
        write_fasta(list(seen_ins.items()), next_query)
        print(f"\n  Unique insertion variants: {len(seen_ins):,}", file=sys.stderr)
        print(f"  Saved: {next_query}", file=sys.stderr)
    if seen_del:
        del_query = f"{args.out}/deletion_variants.fa"
        write_fasta(list(seen_del.items()), del_query)
        print(f"  Unique deletion variants:  {len(seen_del):,}", file=sys.stderr)
        print(f"  Saved: {del_query}", file=sys.stderr)

    # Summary
    summary = {
        "n_v1_records": len(records),
        "n_anchor_pairs": n_pairs,
        "n_observations": len(obs_records),
        "counts": dict(counts),
        "n_unique_insertion_variants": len(seen_ins),
        "n_unique_deletion_variants": len(seen_del),
    }
    with open(f"{args.out}/summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nDone. Summary: {args.out}/summary.json", file=sys.stderr)


if __name__ == "__main__":
    main()
