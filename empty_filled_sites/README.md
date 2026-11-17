# Empty-vs-Filled IS110 Boundary Discovery

Side project: discover precise IS110 element boundaries by comparing genomes
that have an IS110 inserted at a given locus ("filled") with genomes that
have the same locus but without the IS110 ("empty").

## Motivation

Current pipeline uses HMM-derived transposase positions + arbitrary 2 kb flanks
as IS110 references. The 2 kb is a guess — the true IS110 boundary is unknown.
We need an unbiased, data-driven way to define exactly where the IS110 starts
and ends, including its cargo region (which we currently treat as part of the
flank).

The empty-vs-filled approach uses comparative genomics to bracket the IS110:
if we can find genomes where the same flanking sequences are directly adjacent
(no IS110 between them), the gap between those flanks in IS110-containing
genomes IS the IS110 element.

## Existing tools surveyed

| Tool | What it does | Why it doesn't fit |
|---|---|---|
| panISa | Finds IS insertions in BAM (short reads) via clipped reads | Short-read; needs raw read alignments; we have assemblies |
| ISMapper | Finds known IS positions from short reads | Short-read; database-driven; not boundary discovery |
| Flanker | Clusters gene flanking regions | Doesn't pair empty/filled sites |
| ISsaga, ISEScan | IS annotation | Don't do empty/filled comparison |

Conclusion: build it ourselves. The algorithm is straightforward.

## Algorithm

For each IS110 reference (currently a 5.4 kb sequence: `[2 kb flank][IS110 ~1.4 kb][2 kb flank]`):

1. **Extract anchor regions** — short, unique sequences at the outer edges of the flanks:
    - `up_anchor` = positions [0, 500] (outer 500 bp of upstream flank; clearly NOT IS110)
    - `down_anchor` = positions [4900, 5400] (outer 500 bp of downstream flank)

2. **Search both anchors** against all available complete genomes (NCBI E. coli first):
    - minimap2 each anchor against the genome database
    - Filter hits at >=95% identity, >=90% query coverage

3. **Pair anchors per genome**:
    - For each genome, find all up_anchor hits and all down_anchor hits
    - Pair them up if:
        * Same contig
        * Correct relative orientation (consistent strand)
        * Distance between them is plausible (e.g., 0–20 kb)
    - If multiple valid pairs, record all of them (paralogous loci)

4. **Compute the inter-anchor distance**:
    - `distance` = (down_anchor target start) - (up_anchor target end)
    - Distribution interpretation:
        * `~0 bp`: **EMPTY site** — no IS110 ever inserted, or excised cleanly
        * `~4.4 kb`: **FILLED V0 site** — matches the original reference layout
        * `> 4.4 kb`: **FILLED V1+ site** — extra DNA inserted within the IS110
        * `< 4.4 kb` and > 0: deletion within the IS110 (V-deletion case)
        * `< 0` (anchors swapped): inversion or rearrangement; flag separately

5. **Define the IS110 boundary**:
    - From an empty site, the boundary = the position where flanks meet directly.
    - Compare empty vs filled site sequences character-by-character:
        * The aligned-bases-on-each-side define the IS110 ends
        * Any small flanking duplication = target site duplication (TSD)
    - Result: refined `[upstream_flank]|[IS110_BODY]|[downstream_flank]` boundaries.

6. **Export refined IS110 references**:
    - Use empty-site-derived boundaries to trim each reference to its true IS110 body
    - These become much higher-quality queries for downstream multi-split analysis.

## Inputs

- IS110 reference set: `/global/scratch/users/kh36969/is110_refs/is110_reps_95.fa`
  (2,945 references at 95% nt identity, with 2 kb flanks)
- Genome database to search:
    - **Phase A** (development/validation): NCBI complete E. coli — 13,027 genomes (61 GB)
    - **Phase B** (scale-up): Logan E. coli contigs — 480K samples (2.9 TB)

## Outputs

For each IS110 reference:
- `<ref_id>_anchor_pairs.tsv`: per-genome distance distribution
- `<ref_id>_boundaries.json`: refined IS110 start/end coordinates and TSDs
- Aggregate report:
    - Fraction of references with informative empty/filled distributions
    - Refined reference FASTA: `is110_refined.fa` (replaces `is110_reps_95.fa` for downstream)

## Validation criteria

- **Reference must show bimodal distance distribution** (empty + filled peaks).
  If only filled, no empty sites exist → can't refine boundary, keep as-is.
- **Empty-site count >= 5** to call a confident boundary.
- **TSD detection**: small repeat (3–10 bp) flanking the IS110 body in filled
  sites that's seamless in empty sites.

## Phases

### Phase 0 — Pilot (quick, validate concept)
- Pick **10 IS110 references** with the largest cluster sizes (most common variants)
- Run anchor search against the 13K NCBI E. coli genomes
- Generate distance distribution plots and confirm bimodality
- **Decision point**: does the approach produce empty + filled peaks?

### Phase 1 — Scale to all 2,945 references on NCBI
- Run for all references against NCBI complete genomes
- Build refined `is110_refined.fa`

### Phase 2 — Scale to Logan
- Run on Logan E. coli contigs to enrich for rare empty sites
- Logan has 35× more samples; should reveal more excision events

### Phase 3 — Hand-off to main pipeline
- Replace `is110_reps_95.fa` with `is110_refined.fa` in the main pipeline
- Re-run multi-split detection with the refined references

## Edge cases

- **Multi-copy IS110**: same anchor pair might match multiple genomic loci. Filter to within-plausible-distance pairs and report all.
- **Repetitive flanks**: anchors may be non-unique. Quality flag = `n_hits_per_anchor`. Discard refs where anchors are too repetitive.
- **Anchor falls inside another IS110**: anchor itself might be repetitive (some IS elements share flanking motifs). Cross-check against IS110 transposase HMM hits.
- **Inverted insertions**: anchor pairs in opposite orientation indicate IS110 in reverse on the target. Record but treat separately.

## Project structure

```
empty_filled_sites/
├── README.md           # this file
├── scripts/            # implementation
│   ├── extract_anchors.py
│   ├── find_anchor_pairs.py
│   ├── refine_boundary.py
│   └── run_pipeline.sh
├── refs/               # input refs and intermediate FASTAs
└── results/            # per-reference distance distributions, boundaries
```
