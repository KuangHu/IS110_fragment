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
import argparse, json, os, re, shutil, subprocess, sys, hashlib, tempfile
from multiprocessing.pool import ThreadPool
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_alleles import (parse_paf_line, place_per_assembly,  # noqa: E402
                         decompose_alleles, canonical_insert_key, canon_key,
                         revcomp)

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
    p.add_argument("--max-pair-dist", type=int, default=100000,
                   help="Max distance between paired anchor hits (bp, default "
                        "100000; was 500000, which let anchors pair across "
                        "distant repeat copies)")
    p.add_argument("--max-secondary", type=int, default=20,
                   help="minimap2 -N for the anchor search (targets per anchor "
                        "<= N+1; raise it for a census of large DBs)")
    p.add_argument("--max-anchor-overlap", type=int, default=50,
                   help="Anchors may overlap by this much: at an empty site the "
                        "anchors abut, or overlap by the TSD (default 50)")
    p.add_argument("--margin", type=float, default=0.05,
                   help="Min relative score margin of the best anchor pair over "
                        "the runner-up, else the site is ambiguous (default 0.05)")
    p.add_argument("--max-align-cells", type=float, default=2e7,
                   help="Skip SNP-tolerant decomposition above this "
                        "len(V1)*len(variant) (default 2e7)")
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


def union_len(intervals):
    """Total length covered by a set of [start, end) intervals."""
    total, cur_s, cur_e = 0, None, None
    for s, e in sorted(intervals):
        if cur_e is None or s > cur_e:
            if cur_e is not None:
                total += cur_e - cur_s
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        total += cur_e - cur_s
    return total


_CG = re.compile(r"(\d+)([=XID])")


def split_at_indels(qlen, qs, qe, strand, ts, te, cigar, min_gap=10):
    """One PAF alignment -> gap-free(ish) blocks, split at every insertion or
    deletion of >= min_gap bp.

    A block's [v1_start, v1_end) then really is V1 sequence present in the
    variant. Using the alignment's outer span counted a 300-bp deletion inside
    the transposase as covered (CDS 'coverage' 100%, V0 check passed), and a
    variant = V1 + a 300-bp insert came back as ONE block, so the insertion
    was never seen.
    Coordinates: v1 = target, between = query, both 0-based half-open on
    their own forward strands; `strand` is the query's strand vs V1."""
    out = []
    if not cigar:
        return [{"v1_start": ts, "v1_end": te, "between_start": qs, "between_end": qe,
                 "strand": strand, "ident": 0.0, "block_len": qe - qs}]
    t = ts
    q = qs if strand == "+" else qe          # '-' : query consumed from its end
    cur = None

    def close():
        if cur and cur["n"] > 0:
            b0, b1 = sorted((cur["q0"], cur["q1"]))
            out.append({"v1_start": cur["t0"], "v1_end": cur["t1"],
                        "between_start": b0, "between_end": b1, "strand": strand,
                        "ident": cur["m"] / cur["n"] * 100, "block_len": b1 - b0})

    for n, op in ((int(n), op) for n, op in _CG.findall(cigar)):
        if op in "=X" or n < min_gap:
            if cur is None:
                cur = {"t0": t, "t1": t, "q0": q, "q1": q, "m": 0, "n": 0}
            cur["n"] += n
            if op == "=":
                cur["m"] += n
        else:
            close(); cur = None
        if op in "=XD":
            t += n
        if op in "=XI":
            q += n if strand == "+" else -n
        if cur is not None:
            cur["t1"], cur["q1"] = t, q
    close()
    return out


def collinear_plus_cover(blocks):
    """bp of V1 covered by '+' blocks that are in the same ORDER in the
    variant as in V1 (largest first, each kept only if consistent with those
    already kept). 'Clonal' needs this: an inverted or rotated copy has the
    same sequence content and used to be called clonal_V1."""
    kept = []
    for b in sorted((b for b in blocks if b["strand"] == "+"),
                    key=lambda b: -(b["v1_end"] - b["v1_start"])):
        ok = True
        for k in kept:
            before_v1 = b["v1_start"] < k["v1_start"]
            before_bt = b["between_start"] < k["between_start"]
            if before_v1 != before_bt:
                ok = False
                break
        if ok:
            kept.append(b)
    return union_len([(b["v1_start"], b["v1_end"]) for b in kept]), kept


def align_to_own_v1(job):
    """minimap2 of one V1's between-anchor sequences against THAT V1 only.

    All between-sequences used to be mapped in one run against a FASTA of
    every V1 element. A database holds thousands of near-identical copies of
    the same IS, minimap2 reports the best few targets per query, and only
    hits to the query's own V1 were kept -- so whether a perfectly clonal
    site got any block was a lottery (E. coli: 41,204 between-sequences that
    cover a V1 element at >=95% identity had no hit to their own V1 and
    were classed 'unrelated'/'partial'). One target per run removes the
    competition."""
    rid, v1_seq, queries, work_dir = job
    d = tempfile.mkdtemp(prefix="v1_", dir=work_dir)
    try:
        ref, qry = os.path.join(d, "v1.fa"), os.path.join(d, "q.fa")
        write_fasta([("v1", v1_seq)], ref)
        write_fasta(queries, qry)
        res = subprocess.run(["minimap2", "-x", "asm10", "-c", "--eqx", "-t", "1",
                              ref, qry], capture_output=True, text=True)
        if res.returncode != 0:
            return rid, None
        out = defaultdict(list)
        for line in res.stdout.splitlines():
            c = line.split("\t")
            if len(c) < 12:
                continue
            cg = next((x[5:] for x in c[12:] if x.startswith("cg:Z:")), "")
            for blk in split_at_indels(int(c[1]), int(c[2]), int(c[3]), c[4],
                                       int(c[7]), int(c[8]), cg):
                if blk["ident"] >= 95:
                    blk["v1_id"] = rid
                    out[c[0]].append(blk)
        return rid, dict(out)
    finally:
        shutil.rmtree(d, ignore_errors=True)


def extract_target(db, contig, start, end):
    """samtools faidx region."""
    region = f"{contig}:{max(1, start)}-{end}"
    res = subprocess.run(["samtools", "faidx", db, region],
                         capture_output=True, text=True)
    if res.returncode != 0: return None
    return "".join(res.stdout.split("\n")[1:]).strip().upper()


def extract_target_batch(db, regions, batch_size=1000):
    """Batch samtools faidx. regions: iterable of (contig, start, end), 1-based
    inclusive. Returns dict (contig, start, end) -> sequence.

    A failed samtools call raises: it used to `continue`, silently dropping
    up to `batch_size` observations per failure."""
    regions = sorted(set(regions))
    results = {}
    for i in range(0, len(regions), batch_size):
        batch = regions[i:i+batch_size]
        by_name = {f"{c}:{max(1, s)}-{e}": (c, s, e) for c, s, e in batch}
        res = subprocess.run(["samtools", "faidx", db] + list(by_name),
                             capture_output=True, text=True)
        if res.returncode != 0:
            raise RuntimeError(f"samtools faidx failed on a batch of {len(batch)} "
                               f"regions: {res.stderr.strip()[:500]}")
        cur, buf = None, []
        for line in res.stdout.split("\n") + [">"]:
            if line.startswith(">"):
                if cur in by_name:
                    results[by_name[cur]] = "".join(buf).upper()
                cur, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
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
    n_skip = defaultdict(int)
    for r in records:
        # Field name compatibility: records_final uses "is110_id"; build_records.py
        # (Cross_reference_IS Stage 6) uses "ref_id". Accept either.
        rid = r.get("is110_id") or r.get("ref_id")
        if not rid:
            n_skip["record without an id"] += 1
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
            n_skip["missing flank or element sequence"] += 1
            continue
        if len(up) < args.anchor_len or len(down) < args.anchor_len:
            # element within anchor_len of a contig end: no variant search
            n_skip[f"flank shorter than --anchor-len {args.anchor_len}"] += 1
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
    for why, n in sorted(n_skip.items()):
        print(f"  NOT searched: {n:,} records ({why})", file=sys.stderr)
    write_fasta(anchors_list, anchor_fa)
    write_fasta(is110_list, is110_fa)

    # === Step 2: minimap2 anchors against DB ===
    anchor_paf = f"{args.out}/anchors_vs_db.paf"
    # The cached PAF is only valid for THESE anchors, database and -N. It was
    # reused on mere existence: after records.json changed (a Stage 6 rerun)
    # the new anchors were paired from the old alignments.
    with open(anchor_fa, "rb") as fh:
        paf_key = (f"{hashlib.md5(fh.read()).hexdigest()} {os.path.abspath(args.db)} "
                   f"{os.path.getsize(args.db)} {os.path.getmtime(args.db)} "
                   f"N={args.max_secondary}")
    key_file = anchor_paf + ".key"
    cached = (os.path.exists(anchor_paf) and os.path.getsize(anchor_paf) > 0
              and os.path.exists(key_file) and open(key_file).read().strip() == paf_key)
    if not cached:
        print(f"Running minimap2 anchors vs DB ...", file=sys.stderr, flush=True)
        subprocess.run([
            "minimap2", "-x", "asm10", "-c", "--eqx",
            "--secondary=yes", "-N", str(args.max_secondary), "-p", "0.5",
            "-t", str(args.threads), args.db, anchor_fa, "-o", anchor_paf + ".tmp"
        ], check=True)
        os.replace(anchor_paf + ".tmp", anchor_paf)
        with open(key_file, "w") as fh:
            fh.write(paf_key + "\n")
    print(f"  anchor PAF: {anchor_paf} ({os.path.getsize(anchor_paf)/1e6:.1f} MB)",
          file=sys.stderr)

    # === Step 3: Parse anchor hits and find pairs ===
    print(f"Parsing anchor hits ...", file=sys.stderr, flush=True)
    # anchor_id -> [(target, target_start, target_end, strand, ident, cov)]
    anchor_hits = defaultdict(lambda: defaultdict(list))  # qname -> tname -> [Hit]
    paf_rows = defaultdict(int)
    with open(anchor_paf) as f:
        for line in f:
            h = parse_paf_line(line)
            if h is None: continue
            paf_rows[h.qname] += 1
            if h.ident < args.anchor_identity or h.cov < args.anchor_cov: continue
            anchor_hits[h.qname][h.tname].append(h)

    n_capped = sum(1 for n in paf_rows.values() if n > args.max_secondary)
    if n_capped:
        print(f"  WARNING: {n_capped:,} of {len(paf_rows):,} anchors reached "
              f"--max-secondary {args.max_secondary}: the genomes seen are the best "
              f"{args.max_secondary + 1} hits, NOT a census. Counts per category are "
              f"a sample for those sites.", file=sys.stderr)

    # === Step 4: one scored anchor placement per distinct assembly ===
    # Every up x down combination used to become an observation, so a genome
    # with two anchor copies was counted several times and anchors could pair
    # across repeat copies; and `between_end <= between_start` was skipped,
    # which discarded exactly the clean empty sites (anchors abut, gap 0, or
    # overlap by the TSD). See lib_alleles.place_pair.
    print(f"Pairing anchors per V1 ...", file=sys.stderr, flush=True)
    v1_pairs = defaultdict(list)
    n_ambiguous = 0
    for rid in ref_meta:
        ups = anchor_hits.get(f"{rid}__up", {})
        downs = anchor_hits.get(f"{rid}__down", {})
        if not ups or not downs: continue
        placements = place_per_assembly(
            ups, downs, min_gap=-args.max_anchor_overlap,
            max_gap=args.max_pair_dist, margin_ratio=args.margin)
        for pl in placements.values():
            if pl["status"] == "ambiguous":
                n_ambiguous += 1
                continue
            v1_pairs[rid].append({
                "target": pl["tname"], "assembly": pl["assembly"],
                "strand": pl["strand"],
                "between_start": pl["start"], "between_end": pl["end"],
                "between_len": max(0, pl["gap"]),
                "anchor_overlap": max(0, -pl["gap"]),
                "placement_status": pl["status"],
                "placement_margin": round(pl["margin"], 3),
                "up_ident": pl["up"].ident, "down_ident": pl["down"].ident,
            })
    print(f"  Ambiguous placements (not used): {n_ambiguous:,}", file=sys.stderr)

    n_pairs = sum(len(v) for v in v1_pairs.values())
    print(f"  Found {n_pairs:,} anchor pairs across {len(v1_pairs):,} V1s",
          file=sys.stderr)

    # === Step 5: Extract between-anchor sequences (batched) ===
    between_fa = f"{args.out}/between_anchors.fa"
    print(f"Extracting between-anchor sequences (batched) ...",
          file=sys.stderr, flush=True)

    # Build region list
    all_regions = []  # (contig, start, end)
    # Different V1s (the same insertion seen from different source genomes)
    # routinely place on the SAME target interval. A dict keyed by region kept
    # only the last of them; every other V1 silently lost that observation.
    region_to_pair = defaultdict(list)  # (contig, start, end) -> [(rid, i, p)]
    for rid, pairs in v1_pairs.items():
        for i, p in enumerate(pairs):
            p["pair_id"] = f"{rid}__pair{i}__{p['target']}_{p['between_start']}_{p['between_end']}"
            if p["between_len"] == 0:
                p["between_seq"] = ""        # abutting / TSD-overlapping anchors
                continue
            region = (p["target"], p["between_start"] + 1, p["between_end"])
            all_regions.append(region)
            region_to_pair[region].append((rid, i, p))
    print(f"  Total regions to extract: {len(all_regions):,}",
          file=sys.stderr, flush=True)

    # Batch fetch
    seqs = extract_target_batch(args.db, all_regions, batch_size=500)
    print(f"  Fetched {len(seqs):,} sequences", file=sys.stderr, flush=True)

    extracted = []
    n_wrong_len = 0
    for region, seq in seqs.items():
        for rid, i, p in region_to_pair.get(region, ()):
            if len(seq) != p["between_len"]:
                n_wrong_len += 1          # interval runs off the contig
                continue
            sq = revcomp(seq) if p["strand"] == "-" else seq
            p["between_seq"] = sq
            extracted.append((p["pair_id"], sq))
    write_fasta(extracted, between_fa)
    n_lost = sum(1 for ps in v1_pairs.values() for p in ps if "between_seq" not in p)
    print(f"  Wrote {len(extracted):,} between-anchor sequences; "
          f"{n_lost:,} placements without a sequence "
          f"({n_wrong_len:,} of them: fetched length != interval)", file=sys.stderr)

    # === Step 6: Compare each between-seq to ITS OWN V1 (see align_to_own_v1) ===
    print(f"Aligning between-seqs to their own V1 ...", file=sys.stderr, flush=True)
    jobs = []
    for rid, pairs in v1_pairs.items():
        q = [(p["pair_id"], p["between_seq"]) for p in pairs if p.get("between_seq")]
        if q:
            jobs.append((rid, ref_meta[rid]["is110_seq"], q, args.out))
    between_aligns = {}
    n_failed = 0
    with ThreadPool(max(1, args.threads)) as pool:
        for rid, res in pool.imap_unordered(align_to_own_v1, jobs, chunksize=4):
            if res is None:
                n_failed += 1
                continue
            between_aligns.update(res)
    if n_failed:
        raise RuntimeError(f"minimap2 failed for {n_failed} of {len(jobs)} V1s")
    print(f"  {len(between_aligns):,} between-seqs with >=95% blocks to their V1",
          file=sys.stderr)

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
            # bp of V1 covered, as a UNION on V1 coordinates. Summing block
            # lengths counted overlapping alignments twice and a variant with
            # two copies of the element as 200% of V1.
            any_cov = union_len([(b["v1_start"], b["v1_end"]) for b in blocks])
            # ... and only blocks on the SAME strand and in the same order as V1
            total_matched, chain = collinear_plus_cover(blocks)
            rearranged = (any_cov / v1_len >= 0.95 > total_matched / v1_len) if v1_len else False

            v1_cov_frac = total_matched / v1_len if v1_len > 0 else 0
            min_cov_frac = args.min_v1_coverage / 100.0

            # Classification:
            #  empty     : nothing between anchors (< 100 bp)
            #  deletion  : shorter than V1 but contains V1 fragments (>= min_deletion_match bp)
            #  clonal_V1 : matches V1 fully (~same length, >=95% coverage)
            #  insertion : V1's content + extra DNA (>=80% V1 coverage, >=min_extra_bp larger)
            #              OR the variant decomposes cleanly as V1 + one insert
            #              (minimap2 often aligns only one side of V1 across a
            #              large nested insert, leaving too few blocks)
            #  unrelated : virtually no V1 content (<10% coverage)
            #  partial   : everything else (ambiguous)
            decomp = None
            if between_len >= v1_len + args.min_extra_bp:
                v1_seq, seq = ref_meta[rid]["is110_seq"], p["between_seq"]
                decomp = decompose_alleles(
                    v1_seq, seq, min_insert=args.min_extra_bp,
                    tolerant=len(v1_seq) * len(seq) <= args.max_align_cells)
                if decomp["method"] == "none" or decomp["event_class"] not in (
                        "insertion_target_retained", "replacement"):
                    decomp = None
            if between_len < 100:
                category = "empty"
            elif (between_len <= v1_len - args.min_deletion_shortening
                  and total_matched >= args.min_deletion_match):
                category = "deletion"
            elif v1_cov_frac >= 0.95 and between_len < v1_len + args.min_extra_bp:
                category = "clonal_V1"
            elif rearranged and between_len < v1_len + args.min_extra_bp:
                # all of V1 is there, but inverted or in another order
                category = "rearranged"
            elif (len(blocks) >= args.min_similarity_blocks
                  and v1_cov_frac >= min_cov_frac
                  and between_len >= v1_len + args.min_extra_bp):
                category = "insertion"
            elif decomp is not None:
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
                "assembly": p["assembly"],
                "target_strand": p["strand"],
                "anchor_site": {
                    "between_start": p["between_start"],
                    "between_end": p["between_end"],
                    "between_len": between_len,
                    "anchor_overlap": p["anchor_overlap"],
                    "placement_status": p["placement_status"],
                    "placement_margin": p["placement_margin"],
                    "up_anchor_identity": round(p["up_ident"], 1),
                    "down_anchor_identity": round(p["down_ident"], 1),
                },
                "comparison_to_v1": {
                    "v1_len": v1_len,
                    "n_match_blocks": len(blocks),
                    "total_matched_bp": total_matched,
                    "any_strand_matched_bp": any_cov,
                    "blocks": [{"v1_pos": [b["v1_start"], b["v1_end"]],
                                "between_pos": [b["between_start"], b["between_end"]],
                                "strand": b["strand"],
                                "ident": round(b["ident"], 1)} for b in blocks],
                },
            }
            if category == "insertion":
                # Save the full insertion variant sequence + ID
                obs["variant_id"] = f"ins_{short_hash(p['between_seq'])}"
                obs["variant_seq"] = p["between_seq"]
                obs["variant_len"] = between_len
                # Explain the variant as V1 + an insert. The insert key is
                # orientation- and junction-invariant, and the event key adds
                # the V1 context at the insertion point, so the same nested
                # insertion keys identically whichever genome/strand it is in.
                v1_seq, seq = ref_meta[rid]["is110_seq"], p["between_seq"]
                d = decomp
                if d is not None:
                    ctx = v1_seq[max(0, d["offset"] - 100):d["offset"] + 100]
                    ikey = canonical_insert_key(seq, d["insert_start"], d["insert_len"])
                    obs["nested_insert"] = {
                        "method": d["method"], "event_class": d["event_class"],
                        "offset_in_v1": d["offset"],
                        "insert_start": d["insert_start"],
                        "insert_len": d["insert_len"],
                        "junction_microhomology_bp": d["junction_microhomology_bp"],
                        "target_lost_bp": d["target_lost_bp"],
                        "insert_key": ikey,
                        "event_key": f"{ikey}:{canon_key(ctx)}",
                    }
            elif category == "deletion":
                # Save the deletion variant (shorter than V1) + ID
                obs["variant_id"] = f"del_{short_hash(p['between_seq'])}"
                obs["variant_seq"] = p["between_seq"]
                obs["variant_len"] = between_len
            obs_records.append(obs)

    print(f"\n=== Classification counts ===", file=sys.stderr)
    for k in ("empty", "deletion", "clonal_V1", "insertion", "rearranged", "partial", "unrelated"):
        print(f"  {k:>12}: {counts[k]:>8,}", file=sys.stderr)
    n_ev = len({o["nested_insert"]["event_key"] for o in obs_records
                if "nested_insert" in o})
    print(f"  distinct nested-insertion events (event_key): {n_ev:,}", file=sys.stderr)

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
        "n_ambiguous_placements": n_ambiguous,
        "n_anchors_at_secondary_cap": n_capped,
        "max_secondary": args.max_secondary,
        "n_records_not_searched": dict(n_skip),
        "n_placements_without_sequence": n_lost,
        "counts_are": "distinct assemblies (GCA/GCF twins collapsed)",
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
