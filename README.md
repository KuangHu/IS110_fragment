# Cross-Reference IS Boundary Detection

A reusable pipeline that defines the **boundary of any IS element family**
using empty-vs-filled comparative genomics — works for IS110, IS3, IS5, IS200,
IS630, IS66, or any IS family you can identify by an HMM profile.

## Core idea

The transposase HMM tells you *where* an IS element lives. But the precise
**boundaries** of the element (where the IS starts and ends, including any
cargo) are not in the HMM. To find them, we compare two genome states:

- **Filled (V0)**: a genome where the IS is at this locus
- **Empty**: a genome where the IS has been excised, leaving only flanking DNA

Place a pair of "anchors" at distance D from the transposase (one upstream, one
downstream). Search those anchors across a database of genomes:

- Anchors paired close together → empty (no IS between them)
- Anchors paired far apart → V0 (IS still there, plus any cargo)

The difference between the two peak distances **is the IS element length**.
The element's **boundaries** come from decomposing the filled interval against
each empty interval (`filled = empty[:offset] + insert + empty[offset:]`): the
insert's coordinates are the 5' and 3' ends, per side, to the base (up to the
junction microhomology, which is reported).

## Pipeline stages

```
┌──────────────────────────────────────────────────────────────────┐
│ 1. IS detection  (HMM search of target HMM profile across DB)    │
│    → list of (genome, contig, position, strand) for each IS hit  │
│    extension point: protein/blastp or other identifier backends  │
│                                                                  │
│ 2. Anchor extraction  (multi-D: 1, 5, 20, 40, 80 kb)             │
│    → FASTA: ref_id__up<D>, ref_id__down<D>                       │
│                                                                  │
│ 3. Anchor search  (minimap2 anchors vs target genome DB)         │
│    → PAF                                                         │
│                                                                  │
│ 4. Anchor pair finder  (ONE scored placement per assembly)       │
│    → TSV: ref_id, D, target, anchor pair distance, category     │
│                                                                  │
│ 5. Boundary caller  (per ref/D: V0 / empty / V1+ counts,         │
│    element length from empties)                                  │
│    → TSV: ref_id, D, counts, IS_length, anchor_search_capped     │
│                                                                  │
│ 6. Records builder  (empty-vs-filled DECOMPOSITION → per-side    │
│    boundaries, microhomology, TSD, event_key)                    │
│    → JSON: { is_id, boundary, source, empties, fills }           │
│                                                                  │
│ 7. Rearrangement detection  (optional, --detect-rearrangements)  │
│    → TSV: per-(ref, target) inversion/translocation/duplication  │
│                                                                  │
│ 8. Variant discovery       (--detect-growing-tandem)             │
│    For each IS in records.json, search its anchors against DB    │
│    and classify whatever lies between paired anchors:            │
│       empty / clonal / deletion / insertion variant              │
│    → observations.json                                           │
│                                                                  │
│ 9. Lineage building                                              │
│    Per anchor site, sort variants shortest → longest and verify  │
│    each step contains the previous (≥95% id, ≥80% cov).          │
│    → lineages.json (V0 EMPTY ⊂ V1 ⊂ V2 … containment series,     │
│      "direction": "unpolarized" — NOT a time order)              │
│                                                                  │
│10a. V_0 CDS validation                                           │
│    Require the smallest variant of each lineage to still         │
│    contain ≥95% of V_ref's IS110 transposase CDS.                │
│    → lineages_validated.json                                     │
│                                                                  │
│10b. Tandem-repeat detection (self-alignment)                     │
│    For each IS / variant sequence, self-align via minimap2 -X    │
│    to find internal tandem repeats (works even on single seqs).  │
│    → tandem_findings.json/.tsv                                   │
│                                                                  │
│11. Visualisation (optional, on by default)                       │
│    PNG + SVG: lineages (sorted shortest→longest, CDS arrow,      │
│    color-inherited cargo) + tandem-repeat catalog.               │
│    → figures/lineages/, figures/tandem/                          │
└──────────────────────────────────────────────────────────────────┘
```

## Full end-to-end usage (input: genomes → output: lineages, tandems, figures)

```bash
sbatch templates/run_pipeline.sh \
  hmm/PF01548.hmm,hmm/PF02371.hmm \
  /path/to/genome_db.fa \
  /path/to/output_dir/ \
  --detect-rearrangements \
  --detect-growing-tandem
```

Output tree:
```
output_dir/
├── is_hits.tsv                       (stage 1)
├── anchors_refs.fa, anchors_anchors.fa, anchors_table.tsv  (stage 2)
├── anchors_vs_db.paf                 (stage 3)
├── pairs_pairs.tsv, pairs_summary.tsv (stage 4)
├── boundaries_pairs.tsv, boundaries_summary.tsv (stage 5)
├── records/records.json              (stage 6)
├── rearrangements_examples.tsv       (stage 7, if --detect-rearrangements)
├── variants/observations.json        (stage 8)
├── lineages/lineages.json            (stage 9)
├── v0_validated/lineages_validated.json (stage 10a)
├── tandem/tandem_findings.json       (stage 10b)
└── figures/
    ├── lineages/                     (per-lineage PNG/SVG, validated)
    └── tandem/                       (per-tandem PNG, sorted by copy count)
```

**On Stage 1 generality:** Stage 1's job is to define *what counts as an IS hit*.
The current implementation uses HMM profiles (one or more), but this is the
pluggable boundary of the pipeline. Future backends — protein query +
`blastp` against predicted ORFs, nucleotide search via `tblastn`, or any other
identifier — slot in by producing the same Stage 1 output TSV
(`is_id, assembly, contig, tnp_start, tnp_end, tnp_strand, tnp_len, domains_hit`).
Stages 2–7 are IS-family-agnostic and consume only that TSV.

## Correctness notes (2026-09 upgrade, ported from fna_based_mgefinder_project)

Each item below was a measured bug in the fna project and was present here.

| # | was | now |
|---|---|---|
| 1 | `find_variants.py` skipped `between_end <= between_start`, so clean empty sites (anchors abut, gap 0, or overlap by the TSD) were **discarded** | gaps down to `-max_anchor_overlap` are kept; 0-bp empties are counted |
| 2 | every up × down hit pair on a contig was emitted and counted; anchors could pair across repeat copies (up to 200–500 kb apart) | `lib_alleles.place_pair`: one length-normalised, scored placement per assembly, interval cap, runner-up margin → `ambiguous` placements are reported and never used as evidence |
| 3 | empty tolerance = max(200, 10% of 2D): ±8 kb at D=40 kb, so filled IS110 sites read as empty; boundaries were a 50-bp-bin length split **symmetrically** around the transposase | absolute tolerance (`--tol-bp 100` + `--tol-frac 0.002`); boundaries by decomposition (`boundary_evidence.method = empty_vs_filled_decomposition`, fallback labelled `peak_symmetric_fallback`) |
| 4 | self-alignments (`detect_tandem`, `duplication_tandem_db_scan`, lineage all-vs-all) ran minimap2 without `-f 0`, which drops the most frequent minimizers — i.e. the repeats being sought | `-f 0` |
| 5 | GCA_x / GCF_x twins counted twice (E. coli DB: 13,027 assemblies = 7,723 distinct) and counts were per contig | Stage 1 keeps one representative per assembly; every count is per distinct assembly |
| 6 | lineages were described as successive insertions | labelled `"direction": "unpolarized"`; length order is not time order (see `build_lineages.py`) |
| 7 | variants keyed by an exact, orientation-dependent hash | `event_key` (records) and `nested_insert.event_key` (variants): orientation- and junction-invariant insert key + canonical context |

**Anchor-search cap.** minimap2 `-N` limits each anchor to N+1 targets. On the
E. coli DB with the old `-N 50`, 91% of anchors (1,986,882 / 2,179,807) hit
exactly 51 rows, so every count was a *sample* of the DB biased to the closest
matches. `--max-secondary` (env `MAX_SECONDARY` in the template) now sets it;
Stage 5 flags capped anchors in `anchor_search_capped`. For a census set it to
at least the number of distinct assemblies — PAF size grows with it.

Tests: `sbatch tests/run_tests.sh [WORK_DIR]` runs unit tests and a synthetic
end-to-end run with known boundaries (both strands, a GCA/GCF twin, a SNP-divergent
empty and a nested cargo insertion).

## Usage

### Quick start (use defaults — IS110)

```bash
cd /global/home/users/kh36969/Cross_reference_IS
sbatch templates/run_pipeline.sh \
  --hmm hmm/PF01548.hmm,hmm/PF02371.hmm \
  --db /path/to/genome_db.fa \
  --out my_run/
```

### Custom IS family (e.g., IS3)

```bash
sbatch templates/run_pipeline.sh \
  --hmm /path/to/IS3_transposase.hmm \
  --db /path/to/genome_db.fa \
  --out is3_run/
```

### As a Python library (programmatic)

```python
from scripts.pipeline import find_is_boundaries

results = find_is_boundaries(
    hmm_profiles=["PF01548.hmm", "PF02371.hmm"],
    genome_db="genomes.fa",
    output_dir="my_results/",
    distances=[1000, 5000, 20000, 40000, 80000],
    min_identity=95,
    threads=32,
)
# results: dict of {is_id: {boundary, v0_peak, empty_peak, ...}}
```

## Module layout

```
Cross_reference_IS/
├── README.md                         — this file
├── scripts/
│   ├── lib_alleles.py                — shared placement / decomposition / dedup helpers
│   ├── is_detect.py                  — HMM-based IS detection (Stage 1)
│   ├── extract_anchors.py            — multi-D anchor FASTA builder
│   ├── find_anchor_pairs.py          — pair up/down anchors per target
│   ├── call_boundaries.py            — histogram → V0/empty peaks → IS size
│   ├── build_records.py              — final JSON records
│   ├── find_rearrangements.py        — inversion/translocation/duplication (Stage 7)
│   ├── rearrangement_validator.py    — MUMmer-based verdict on Stage 7 events
│   ├── find_variants.py              — per-anchor-site variant discovery (Stage 8)
│   ├── build_lineages.py             — nested lineages from variants (Stage 9)
│   ├── validate_v0_cds.py            — confirm V_0 retains the IS110 CDS (Stage 10a)
│   ├── detect_tandem.py              — self-alignment tandem-repeat scan (Stage 10b)
│   ├── viz_lineage.py                — per-lineage PNG/SVG (Stage 11)
│   ├── viz_tandem.py                 — per-tandem PNG (Stage 11)
│   ├── verify_tandem_enrichment.py   — verification: TRF+ULTRA tandem enrichment vs controls
│   └── pipeline.py                   — full orchestration (Python API + CLI)
├── templates/
│   ├── run_pipeline.sh               — SLURM submitter template
│   └── run_pipeline_array.sh         — array-job template for big runs
├── examples/
│   ├── is110_ecoli.config.json       — config used for the original IS110 run
│   └── is3_example.config.json       — example for IS3
└── hmm/                              — bundled HMM profiles (symlinked from repo)
```

## Rearrangement detection (optional)

The same anchor PAF that drives boundary calling also encodes structural
rearrangement signals. Two scripts handle this:

**`find_rearrangements.py`** — fast PAF-only classifier. For each
`(ref_id, target_assembly)` pair:

- Same contig + same strand → normal (already handled by Stages 4–6)
- Same contig + opposite strands → **inversion**
- Anchors on different contigs → **translocation** or contig break
- Multiple hits per anchor side → **duplication**

Run via the pipeline:

```bash
python3 scripts/pipeline.py --hmm hmm/PF01548.hmm,hmm/PF02371.hmm \
    --db genomes.fa --out my_run/ --detect-rearrangements
```

Or standalone on an existing PAF:

```bash
python3 scripts/find_rearrangements.py --paf my_run/anchors_vs_db.paf \
    --out my_run/rearrangements --anchor-D 1000
```

Output: `<out>_examples.tsv` — top examples per category for review.

**`rearrangement_validator.py`** — MUMmer-based confirmation. Takes the
`_examples.tsv` plus per-assembly genome FASTAs and runs `nucmer + show-diff`
on each candidate event, producing a `CONFIRMED / LIKELY / REJECTED` verdict
plus an `is_mediated` flag (true if the breakpoint sits within ±2 kb of an IS
transposase from Stage 1's `is_hits.tsv`).

```bash
python3 scripts/rearrangement_validator.py \
    my_run/rearrangements_examples.tsv \
    /path/to/genome_dir/ \
    my_run/is_hits.tsv \
    my_run/validation/ --n-per-cat 10
```

The validator is intentionally **not** wired into the pipeline — it needs
MUMmer installed and a per-assembly genome directory layout
(`<genome_dir>/<assembly>/*_genomic.fna`) which is too project-specific to
assume.

## Verification modules

The pipeline *discovers* events (rearrangements, tandem repeats); the
verification modules *prove* them with independent gold-standard tools. Each
discovery type has a matching verifier.

### Rearrangements → `rearrangement_validator.py`
Confirms Stage 7 inversion / translocation / duplication calls with MUMmer
(`nucmer` + `show-diff`) and flags whether the breakpoint is IS-mediated. See the
"Rearrangement detection" section above. (Roadmap: add a second caller — `syri`
— so an event is confirmed only on agreement of two independent tools.)

### Tandem repeats → `verify_tandem_enrichment.py`
Tests the claim **"tandem repeats are enriched in IS elements vs normal DNA."**
A single tandem call proves nothing — tandems occur in all DNA — so this module
measures the *background rate* and asks whether IS elements exceed it.

- **Two independent callers**: TRF (alignment-score model) + ULTRA (HMM model).
  A sequence is counted tandem-positive only if **both** agree.
- **Two null controls** per IS element:
  - `genomic` — length-matched random windows from the same DB (excluding IS
    loci). Tests enrichment over the genome at large (composition + mechanism).
  - `shuffle` — dinucleotide-preserving shuffle of each IS (Altschul–Erikson;
    identical 1-mer + 2-mer frequencies, scrambled order). Isolates a
    **mechanistic** signal beyond base composition.
- **Metrics**: per-sequence tandem density + presence, split into *any period*
  and *large* (period ≥ 50 bp, the structural-duplication track).
- **Statistics**: Fisher exact (presence → odds ratio), Mann-Whitney U (density),
  label-permutation (10⁴×), reported as **fold-enrichment + p-values**.

```bash
python3 scripts/verify_tandem_enrichment.py \
    --records my_run/records/records.json \
    --genome-db genomes.fa \
    --out my_run/tandem_verify/ \
    --control-ratio 3 --min-large-period 50 --threads 16
```

Output: `tandem_verify/enrichment_report.txt` (+ `.json`, `per_sequence.tsv`).
Interpretation: enriched vs **both** controls (especially `shuffle`) ⇒ the IS
process itself generates tandem arrays beyond what sequence composition predicts.

## Required tools

Core pipeline (Stages 1–11):
- minimap2 ≥ 2.30
- hmmer (hmmsearch)
- samtools (with .fai index of DB)
- prodigal or pyrodigal (ORF calling for HMM on raw genomes)
- **pysam** (in-process FASTA extraction in Stages 2 & 6 — required, not optional)
- Python 3.8+; matplotlib + pyGenomeViz for visualisation (Stage 11)

Verification modules (optional, run separately):
- **TRF** + **ULTRA** — tandem callers for `verify_tandem_enrichment.py`
- **scipy** + **numpy** — statistics for `verify_tandem_enrichment.py`
- **MUMmer** (nucmer, show-diff) — for `rearrangement_validator.py`

> **Performance note:** Stages 2 (extract_anchors) and 6 (build_records) use
> `pysam` for in-process FASTA reads, parallelised across `--threads`. This is
> essential at scale — the earlier per-record `samtools faidx` subprocess design
> took >24 h on a 150 GB / 3 M-contig DB; pysam brings it to minutes. Stage 1
> (`is_detect.py`) streams contigs to pyrodigal, caches `proteins.faa` on rerun,
> and filters proteins > 10 kb aa before hmmsearch (hmmer rejects > 100 kb aa).

## Configuration

See `examples/is110_ecoli.config.json` for the full schema. Key fields:

```json
{
  "hmm_profiles": ["hmm/PF01548.hmm", "hmm/PF02371.hmm"],
  "require_all_domains": true,
  "anchor_distances_bp": [1000, 5000, 20000, 40000, 80000],
  "anchor_length_bp": 500,
  "min_identity_pct": 95,
  "min_coverage_pct": 80,
  "histogram_bin_bp": 50,
  "expected_is_size_range": [800, 200000]
}
```

## Origin

Extracted from the IS110_fragment project (April 2026). The original pipeline
analyzed 34,406 confident IS110 boundaries across 13K NCBI E. coli genomes.
This is a generalized version for any IS family.
