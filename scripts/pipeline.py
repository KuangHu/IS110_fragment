#!/usr/bin/env python3
"""Orchestrates the full Cross-Reference IS boundary detection pipeline.

Stages:
  1. is_detect.py             — HMM-based IS discovery (extension point: protein/blastp)
  2. extract_anchors.py       — multi-D anchor FASTA from IS hits + indexed FASTA
  3. minimap2                 — anchors vs target DB
  4. find_anchor_pairs.py     — pair up/down anchors per target
  5. call_boundaries.py       — peak-detect V0 + empty distances
  6. build_records.py         — final per-IS JSON records (boundaries called)
  7. find_rearrangements.py   — (optional --detect-rearrangements)
                                  inversion/translocation/duplication signatures
  8. find_variants.py         — (optional --detect-growing-tandem)
                                  per-IS anchor-site variant discovery (empty / clonal
                                  / deletion / insertion)
  9. build_lineages.py        — size-sorted nested lineages per anchor site
 10a. validate_v0_cds.py      — V_0 must retain ≥95% of V_ref's IS110 CDS
 10b. detect_tandem.py        — self-alignment tandem-repeat detection
 11. viz_lineage.py/viz_tandem.py — figures (optional, on by default with --detect-growing-tandem)

Each stage is cacheable. A stage is skipped only when its completion marker
(<out>/.stage_<name>.done, written after the stage exits successfully) AND its
outputs exist. Output existence alone is not enough: a stage killed mid-write
(SLURM time limit, OOM) leaves a partial file -- a truncated PAF silently
removes anchors from every later count. Outputs without a marker are re-made;
pass --trust-existing to adopt the outputs of a run made before markers
existed. Stage 1 is the exception: is_hits.tsv is written atomically by
is_detect.py and may also be supplied by hand, so its existence is trusted.

CLI:
    # Boundary call only:
    python3 pipeline.py --hmm A.hmm,B.hmm --db genomes.fa --out my_run/

    # Boundary + rearrangements + growing/tandem analysis + figures:
    python3 pipeline.py --hmm A.hmm,B.hmm --db genomes.fa --out my_run/ \\
        --detect-rearrangements --detect-growing-tandem

Python API:
    from pipeline import find_is_boundaries, detect_growing_and_tandem
"""
import argparse, os, subprocess, sys


_STAGE_ORDER = ["2_anchors", "3_minimap2", "4_pairs", "5_boundaries", "6_records",
                "7_rearrangements", "8_variants", "9_lineages", "10a_v0_cds",
                "10b_tandem"]


_RAN = set()   # output dirs in which some stage has executed in this process


def _run_stage(output_dir, name, title, cmd, outputs, trust_existing=False,
               atomic_o=None, inputs=()):
    """Run one cached stage.

    outputs   files that must all exist for a skip
    inputs    files whose modification invalidates the stage (newer than marker)
    atomic_o  index in `cmd` of an output path the tool writes itself; it is
              written to <path>.tmp and renamed on success

    The marker stores the exact command line. A stage is skipped only if the
    marker exists, records the SAME command (thresholds, DB paths, ...), all
    outputs exist and no input is newer than the marker. A marker that only
    said "done" let a rerun with --min-identity 100 return the 95% results.
    """
    import json
    marker = os.path.join(output_dir, f".stage_{name}.done")
    # Written when an EARLIER stage re-runs; survives a crash/restart, so a
    # later `--trust-existing` cannot adopt outputs that the previous
    # process had already invalidated (that state used to live only in
    # memory).
    invalid = os.path.join(output_dir, f".stage_{name}.invalid")
    sig = json.dumps([os.path.basename(cmd[0])] + list(cmd[1:]))
    have = all(os.path.exists(o) for o in outputs)
    if have and os.path.exists(marker):
        with open(marker) as f:
            old_sig = f.read().strip()
        stale_in = [i for i in inputs
                    if os.path.exists(i) and os.path.getmtime(i) > os.path.getmtime(marker)]
        if old_sig == sig and not stale_in:
            return
        why = (f"input changed: {os.path.basename(stale_in[0])}" if stale_in
               else "parameters changed since the cached run")
        print(f"=== {title}: {why} -- re-running ===", flush=True)
    elif (have and trust_existing and output_dir not in _RAN
          and not os.path.exists(invalid)
          and all(os.path.getsize(o) > 0 for o in outputs)):
        # Adoption is only for an untouched legacy directory. Once any stage
        # has run in this invocation its downstream outputs are stale and
        # must never be adopted.
        print(f"=== {title}: adopting existing outputs (--trust-existing) ===",
              flush=True)
        with open(marker, "w") as f:
            f.write(sig + "\n")
        return
    elif any(os.path.exists(o) for o in outputs):
        print(f"=== {title}: outputs present but no completion marker "
              f"(interrupted or pre-marker run) -- re-running ===", flush=True)
    # This stage is about to change its outputs, so every later stage's
    # result is stale: drop their markers too (they re-run in order).
    for later in _STAGE_ORDER[_STAGE_ORDER.index(name):]:
        m = os.path.join(output_dir, f".stage_{later}.done")
        if os.path.exists(m):
            os.remove(m)
        open(os.path.join(output_dir, f".stage_{later}.invalid"), "w").close()
    _RAN.add(output_dir)
    print(f"=== {title} ===", flush=True)
    if atomic_o is not None:
        final = cmd[atomic_o]
        run_cmd = list(cmd)
        run_cmd[atomic_o] = final + ".tmp"
        subprocess.run(run_cmd, check=True)
        os.replace(final + ".tmp", final)
    else:
        subprocess.run(cmd, check=True)
    missing = [o for o in outputs if not os.path.exists(o)]
    if missing:
        raise RuntimeError(f"{title}: finished without writing {missing}")
    with open(marker, "w") as f:
        f.write(sig + "\n")
    if os.path.exists(invalid):
        os.remove(invalid)


def find_is_boundaries(
    hmm_profiles,           # list[str] of HMM profile paths
    genome_db,              # indexed FASTA (.fai present)
    output_dir,
    distances=(1000, 5000, 20000, 40000, 80000),
    anchor_length=500,
    flank=80000,
    min_identity=95,
    min_coverage=80,
    histogram_bin=50,
    expected_is_size_range=(800, 200000),
    require_all_domains=True,
    threads=32,
    source_db=None,         # optional — defaults to genome_db
    detect_rearrangements=False,
    anchor_d_rearrangements=1000,
    max_secondary=50,
    trust_existing=False,
):
    os.makedirs(output_dir, exist_ok=True)
    # An anchor at distance D needs D + anchor_length bp of flank; with the
    # old defaults (flank 80000, D up to 80000) the largest distance was
    # silently never generated.
    need = max(distances) + anchor_length
    if flank < need:
        print(f"NOTE: --flank {flank} < max distance + anchor length; using {need}",
              flush=True)
        flank = need
    script_dir = os.path.dirname(os.path.abspath(__file__))

    is_hits_tsv = os.path.join(output_dir, "is_hits.tsv")
    anchors_prefix = os.path.join(output_dir, "anchors")
    refs_fa = anchors_prefix + "_refs.fa"
    anchors_fa = anchors_prefix + "_anchors.fa"
    ref_table_tsv = anchors_prefix + "_table.tsv"
    anchors_paf = os.path.join(output_dir, "anchors_vs_db.paf")
    pairs_prefix = os.path.join(output_dir, "pairs")
    boundaries_prefix = os.path.join(output_dir, "boundaries")
    records_dir = os.path.join(output_dir, "records")

    src_db = source_db or genome_db

    # Stage 1
    if not os.path.exists(is_hits_tsv):
        cmd = ["python3", os.path.join(script_dir, "is_detect.py"),
               "--db", src_db,
               "--hmm", ",".join(hmm_profiles),
               "--out", is_hits_tsv,
               "--threads", str(threads)]
        if not require_all_domains:
            cmd.append("--any-domain")
        print("=== Stage 1: IS detection ===", flush=True)
        subprocess.run(cmd, check=True)

    def stage(name, title, cmd, outputs, **kw):
        _run_stage(output_dir, name, title, cmd, outputs,
                   trust_existing=trust_existing, **kw)

    # Stage 2
    stage("2_anchors", "Stage 2: anchor extraction",
          ["python3", os.path.join(script_dir, "extract_anchors.py"),
           "--hits", is_hits_tsv,
           "--genome-db", src_db,
           "--out", anchors_prefix,
           "--flank", str(flank),
           "--anchor-length", str(anchor_length),
           "--distances", ",".join(str(d) for d in distances),
           "--threads", str(threads)],
          [anchors_fa, refs_fa, ref_table_tsv],
          # Stage 1 is cached on existence; a regenerated or hand-replaced
          # is_hits.tsv must still invalidate everything built from it
          # ... and so must a rebuilt database (its .fai is rewritten)
          inputs=[is_hits_tsv, src_db + ".fai"])

    # Stage 3
    # Use 'sr' preset for short anchor queries (500 bp).
    # -N caps secondary alignments PER ANCHOR: each anchor reports at most
    # N+1 targets. On the 13,027-assembly E. coli DB with -N 50, 91% of
    # anchors (1,986,882 / 2,179,807) hit exactly 51 rows -- every count
    # downstream is a sample of the DB, biased to the closest matches.
    # Raise --max-secondary to at least the number of distinct assemblies
    # for a census (PAF size grows ~linearly with it).
    cmd = ["minimap2", "-x", "sr", "-c", "--eqx",
           "--secondary=yes", "-N", str(max_secondary), "-p", "0.5",
           "-t", str(threads),
           genome_db, anchors_fa, "-o", anchors_paf]
    stage("3_minimap2", "Stage 3: anchors vs DB (minimap2)", cmd, [anchors_paf],
          atomic_o=len(cmd) - 1, inputs=[anchors_fa, genome_db + ".fai"])

    # Stage 4
    stage("4_pairs", "Stage 4: anchor pair finder",
          ["python3", os.path.join(script_dir, "find_anchor_pairs.py"),
           "--paf", anchors_paf,
           "--ref-table", ref_table_tsv,
           "--out", pairs_prefix,
           "--min-identity", str(min_identity),
           "--min-coverage", str(min_coverage)],
          [pairs_prefix + "_pairs.tsv", pairs_prefix + "_summary.tsv"],
          inputs=[anchors_paf, ref_table_tsv])

    # Stage 5
    stage("5_boundaries", "Stage 5: boundary peak caller",
          ["python3", os.path.join(script_dir, "call_boundaries.py"),
           "--paf", anchors_paf,
           "--ref-table", ref_table_tsv,
           "--out", boundaries_prefix,
           "--min-identity", str(min_identity),
           "--min-coverage", str(min_coverage),
           "--histogram-bin", str(histogram_bin),
           "--min-is-size", str(expected_is_size_range[0]),
           "--max-is-size", str(expected_is_size_range[1]),
           "--max-secondary", str(max_secondary)],
          [boundaries_prefix + "_pairs.tsv", boundaries_prefix + "_summary.tsv"],
          inputs=[anchors_paf, ref_table_tsv])

    # Stage 6
    records_json = os.path.join(records_dir, "records.json")
    os.makedirs(records_dir, exist_ok=True)
    stage("6_records", "Stage 6: build records",
          ["python3", os.path.join(script_dir, "build_records.py"),
           "--paf", anchors_paf,
           "--ref-table", ref_table_tsv,
           "--hits", is_hits_tsv,
           "--genome-db", genome_db,
           "--out", records_dir,
           "--source-db", src_db,
           "--min-identity", str(min_identity),
           "--min-coverage", str(min_coverage),
           # same size limits as Stage 5, so the two stages cannot disagree
           # about what an empty site is
           "--min-is-size", str(expected_is_size_range[0]),
           "--max-is-size", str(expected_is_size_range[1]),
           "--threads", str(threads)],
          [records_json] + [os.path.join(records_dir, f) for f in (
              "summary.tsv", "is_elements.fa", "transposase_cds.fa",
              "upstream_flanks.fa", "downstream_flanks.fa", "empty_junctions.fa")],
          inputs=[anchors_paf, ref_table_tsv, is_hits_tsv,
                  genome_db + ".fai", src_db + ".fai"])

    # Stage 7 (optional): rearrangement detection from the same anchor PAF
    if detect_rearrangements:
        rearr_prefix = os.path.join(output_dir, "rearrangements")
        stage("7_rearrangements", "Stage 7: rearrangement detection",
              ["python3", os.path.join(script_dir, "find_rearrangements.py"),
               "--paf", anchors_paf,
               "--out", rearr_prefix,
               "--min-identity", str(min_identity),
               "--min-coverage", str(min_coverage),
               "--anchor-D", str(anchor_d_rearrangements),
               # contig lengths: without them every translocation candidate
               # got anchor depth 0 and could not be told from a contig break
               "--fai", genome_db + ".fai"],
              [rearr_prefix + "_examples.tsv"], inputs=[anchors_paf])

    return records_json


def detect_growing_and_tandem(
    records_json,             # records.json from Stage 6
    genome_db,                # indexed FASTA (same as Stage 1 target)
    output_dir,               # same out dir as Stage 1-7 (sub-dirs added)
    threads=32,
    min_v1_coverage=80,       # require >=80% V_ref covered for "downstream" variant
    min_step_bp=10,           # in lineage, min length diff between adjacent steps
    cluster_id=0.99,          # cluster near-identical variants
    min_cds_coverage=95,      # V_0 must retain >=95% of V_ref's CDS to be valid
    min_tandem_repeat_len=50,
    min_tandem_copies=2,
    render_figures=True,
    trust_existing=False,
):
    """Stages 8-11: variant discovery, lineage building, CDS validation, tandem
    detection, optional figure rendering. Outputs land in `output_dir/`."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.makedirs(output_dir, exist_ok=True)

    variants_dir = os.path.join(output_dir, "variants")
    lineage_dir = os.path.join(output_dir, "lineages")
    validated_dir = os.path.join(output_dir, "v0_validated")
    tandem_dir = os.path.join(output_dir, "tandem")
    figs_dir = os.path.join(output_dir, "figures")

    def stage(name, title, cmd, outputs):
        _run_stage(output_dir, name, title, cmd, outputs,
                   trust_existing=trust_existing)

    # Stage 8: find anchor-site variants (V_ref vs target → empty/clonal/deletion/insertion)
    obs_json = os.path.join(variants_dir, "observations.json")
    os.makedirs(variants_dir, exist_ok=True)
    stage("8_variants", "Stage 8: anchor-site variant discovery",
          ["python3", os.path.join(script_dir, "find_variants.py"),
           "--records", records_json,
           "--db", genome_db,
           "--out", variants_dir,
           "--threads", str(threads),
           "--min-v1-coverage", str(min_v1_coverage)],
          [obs_json])

    # Stage 9: build size-sorted nested lineages
    lineages_json = os.path.join(lineage_dir, "lineages.json")
    os.makedirs(lineage_dir, exist_ok=True)
    stage("9_lineages", "Stage 9: build lineages",
          ["python3", os.path.join(script_dir, "build_lineages.py"),
           "--obs", obs_json,
           "--out", lineage_dir,
           "--records", records_json,
           "--cluster-id", str(cluster_id),
           "--min-step-bp", str(min_step_bp),
           "--threads", str(min(threads, 8))],
          [lineages_json])

    # Stage 10a: validate V_0 retains the IS110 CDS
    validated_json = os.path.join(validated_dir, "lineages_validated.json")
    os.makedirs(validated_dir, exist_ok=True)
    stage("10a_v0_cds", "Stage 10a: CDS validation of V_0",
          ["python3", os.path.join(script_dir, "validate_v0_cds.py"),
           "--lineages", lineages_json,
           "--obs", obs_json,
           "--records", records_json,
           "--out", validated_dir,
           "--min-cds-coverage", str(min_cds_coverage)],
          [validated_json])

    # Stage 10b: tandem detection (self-alignment, catches even single-observation tandems)
    tandem_json = os.path.join(tandem_dir, "tandem_findings.json")
    os.makedirs(tandem_dir, exist_ok=True)
    stage("10b_tandem", "Stage 10b: tandem-repeat detection",
          ["python3", os.path.join(script_dir, "detect_tandem.py"),
           "--records", records_json,
           "--obs", obs_json,
           "--out", tandem_dir,
           "--min-repeat-len", str(min_tandem_repeat_len),
           "--min-copies", str(min_tandem_copies),
           "--threads", str(min(threads, 8))],
          [tandem_json])

    # Stage 11 (optional): render figures
    if render_figures:
        os.makedirs(figs_dir, exist_ok=True)
        # Tandem figures (all)
        tandem_fig_dir = os.path.join(figs_dir, "tandem")
        if not os.path.isdir(tandem_fig_dir) or not os.listdir(tandem_fig_dir):
            os.makedirs(tandem_fig_dir, exist_ok=True)
            cmd = ["python3", os.path.join(script_dir, "viz_tandem.py"),
                   "--findings", tandem_json,
                   "--out", tandem_fig_dir]
            print("=== Stage 11a: render tandem figures ===", flush=True)
            subprocess.run(cmd, check=True)
        # Lineage figures (multi-step only, validated)
        import json as _json
        with open(validated_json) as f:
            validated = _json.load(f)
        lineage_fig_dir = os.path.join(figs_dir, "lineages")
        os.makedirs(lineage_fig_dir, exist_ok=True)
        n_rendered = 0
        for v1_id, info in validated.items():
            for li, chain in enumerate(info["lineages"]):
                if len([s for s in chain if s.get("type") != "empty"]) < 2:
                    continue
                safe = v1_id.replace("|", "_").replace(".", "_")
                out_prefix = os.path.join(lineage_fig_dir, f"{safe}_lin{li}")
                if os.path.exists(out_prefix + ".png"):
                    continue
                cmd = ["python3", os.path.join(script_dir, "viz_lineage.py"),
                       "--lineages", validated_json,
                       "--obs", obs_json,
                       "--records", records_json,
                       "--v1-id", v1_id,
                       "--lineage-idx", str(li),
                       "--out", out_prefix,
                       "--threads", "4"]
                try:
                    subprocess.run(cmd, check=True, capture_output=True)
                    n_rendered += 1
                except Exception:
                    pass
        print(f"=== Stage 11b: rendered {n_rendered} lineage figures ===", flush=True)

    return {
        "observations": obs_json,
        "lineages": lineages_json,
        "lineages_validated": validated_json,
        "tandem_findings": tandem_json,
        "figures_dir": figs_dir if render_figures else None,
    }


def main():
    p = argparse.ArgumentParser(description="Cross-Reference IS boundary pipeline")
    p.add_argument("--hmm", required=True, help="Comma-separated HMM profile paths")
    p.add_argument("--db", required=True, help="Target genome DB FASTA (indexed)")
    p.add_argument("--source-db", help="Where to detect IS (defaults to --db)")
    p.add_argument("--out", required=True, help="Output directory")
    p.add_argument("--distances", default="1000,5000,20000,40000,80000")
    p.add_argument("--anchor-length", type=int, default=500)
    p.add_argument("--flank", type=int, default=80000)
    p.add_argument("--min-identity", type=float, default=95)
    p.add_argument("--min-coverage", type=float, default=80)
    p.add_argument("--histogram-bin", type=int, default=50)
    p.add_argument("--min-is-size", type=int, default=800)
    p.add_argument("--max-is-size", type=int, default=200000)
    p.add_argument("--any-domain", action="store_true")
    p.add_argument("--threads", type=int, default=32)
    p.add_argument("--max-secondary", type=int, default=50,
                   help="minimap2 -N for the anchor search: max target HITS per "
                        "anchor is N+1. Default 50 samples the DB. For a census "
                        "it must exceed the number of DB sequences an anchor can "
                        "hit -- GCA/GCF twins and extra locus copies each use a "
                        "slot, so the number of distinct assemblies is NOT "
                        "enough; raise it until Stage 5 reports no capped anchors")
    p.add_argument("--trust-existing", action="store_true",
                   help="Adopt stage outputs that exist but have no completion "
                        "marker (runs made before markers existed). Without it "
                        "such outputs are treated as interrupted and re-made")
    p.add_argument("--detect-rearrangements", action="store_true",
                   help="Run Stage 7: classify inversions/translocations/duplications "
                        "from the anchor PAF (off by default)")
    p.add_argument("--anchor-D-rearrangements", type=int, default=1000,
                   help="Anchor distance to use for rearrangement classification "
                        "(default 1000; must be one of --distances)")
    # Stages 8-11 (variant/lineage/tandem analysis)
    p.add_argument("--detect-growing-tandem", action="store_true",
                   help="Run Stages 8-11: per-anchor-site variant discovery, "
                        "lineage building, CDS-validation, tandem detection, figures.")
    p.add_argument("--min-v1-coverage", type=float, default=80)
    p.add_argument("--min-step-bp", type=int, default=10)
    p.add_argument("--cluster-id", type=float, default=0.99)
    p.add_argument("--min-cds-coverage", type=float, default=95)
    p.add_argument("--min-tandem-repeat-len", type=int, default=50)
    p.add_argument("--min-tandem-copies", type=int, default=2)
    p.add_argument("--no-render-figures", action="store_true",
                   help="Skip the figure-rendering stage (Stage 11).")
    args = p.parse_args()

    result = find_is_boundaries(
        hmm_profiles=args.hmm.split(","),
        genome_db=args.db,
        source_db=args.source_db,
        output_dir=args.out,
        distances=[int(d) for d in args.distances.split(",")],
        anchor_length=args.anchor_length,
        flank=args.flank,
        min_identity=args.min_identity,
        min_coverage=args.min_coverage,
        histogram_bin=args.histogram_bin,
        expected_is_size_range=(args.min_is_size, args.max_is_size),
        require_all_domains=not args.any_domain,
        threads=args.threads,
        detect_rearrangements=args.detect_rearrangements,
        anchor_d_rearrangements=args.anchor_D_rearrangements,
        max_secondary=args.max_secondary,
        trust_existing=args.trust_existing,
    )
    print(f"\nFinal records: {result}", flush=True)

    if args.detect_growing_tandem:
        ext = detect_growing_and_tandem(
            records_json=result,
            genome_db=args.db,
            output_dir=args.out,
            threads=args.threads,
            min_v1_coverage=args.min_v1_coverage,
            min_step_bp=args.min_step_bp,
            cluster_id=args.cluster_id,
            min_cds_coverage=args.min_cds_coverage,
            min_tandem_repeat_len=args.min_tandem_repeat_len,
            min_tandem_copies=args.min_tandem_copies,
            render_figures=not args.no_render_figures,
            trust_existing=args.trust_existing,
        )
        print(f"\nStages 8-11 outputs:", flush=True)
        for k, v in ext.items():
            print(f"  {k}: {v}", flush=True)


if __name__ == "__main__":
    main()
