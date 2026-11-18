The current outputs are not reliable enough for the paper’s copy counts or lineage claims. I verified the failures below against HEAD `696e5e6` on `review/2026-11`. The repository is unchanged; synthetic fixtures and reproducers are under `<tmp>/`.

1. **High — `scripts/detect_is110_variants.py:131`: pairs can join different loci into a fictitious insertion.**  
   Input: two genuine copies, each split into Q`[0,1000)` and Q`[1000,2000)`, at T`[0,1000)/[1600,2600)` and T`[4000,5000)/[5600,6600)`. Actual: **three insertions**, including the first copy’s left arm joined to the second copy’s right arm. Expected: two. Verified in `test_core.py`: `two_events` returned three records. Strand, order and query-junction checks cannot establish locus membership. **Fix:** select compatible, non-reusing chains per physical locus; reject or flag competing assignments instead of emitting every permissible pair.

2. **High — `scripts/detect_is110_variants.py:238`: multiple CIGAR deletions inflate coverage and fabricate the reported “parts.”**  
   Input: `400=600I400=600I400=`, query length 2400, `--min-query-cov 60`. Only 1200 query bases align: 50% coverage. Actual: two deletions, each reported at **75% coverage**, on both strands. Expected: neither passes 60%. Verified through the real CLI in `test_followups.py`. Each `row_event()` subtracts only its own `I`; other deleted regions remain inside its invented query intervals. **Fix:** walk the CIGAR once into actual aligned blocks, exclude every `I` from aligned query coverage, and derive each event from its adjacent blocks.

3. **High — `scripts/detect_multisplit.py:100`, `:192`: one chain per contig loses copies and manufactures “intact” observations.**  
   Verified in `test_core.py`:
   - Two intact copies → one retained.
   - Opposite-strand intact copies → the entire nondominant strand disappears.
   - One 50%-length hit → zero splits, `event_types=intact`.
   - Non-collinear halves from separate loci → a 100%-coverage chain with `other`, still counted in the zero-split distribution.
   - A single row containing two 600-bp deletions → 100% coverage and zero splits.

   **Fix:** retain multiple locus chains on both strands, expand CIGAR events, distinguish `partial` and `unresolved` from `intact`, and require near-full aligned coverage for intact calls.

4. **High — `scripts/find_lineage.py:208`, `:239`, `:285`: nesting is never tested, and cluster IDs cannot be compared across levels.**  
   Input: intact IS, one insertion A, then two unrelated insertions B/C. Actual: a three-level ladder; A gets cluster `0`, while B/C independently get `0,1`. Expected: no supported inheritance edge. Verified in `test_lineage.py`; clustering A/B/C together returned `[0,1,2]`. The script also buckets deletion/complex split counts as possible insertion ladders. **Fix:** collect insertions across all levels, assign global sequence clusters, match inherited fragments at corresponding query junctions, and require a connected chain of supported observations before reporting a ladder.

5. **High — `scripts/hmm_logan_pipeline.sh:274`, `:349`: contig cluster members never re-enter the census.**  
   Input: identical contigs from SRR1 and SRR2 clustered under SRR1; its predicted protein has both domains. Actual: one position/sample, SRR1. Expected: two samples. Verified by executing the embedded mapping block in `test_pipeline.py`: `contig_members=2 position_rows=1`. **Fix:** predict/search all original contigs, or expand contig membership with coordinate-aware alignment and verify transferred ORFs/domains. Protein-cluster expansion alone cannot restore discarded contigs.

6. **High — `scripts/extract_is110_refs.sh:150`: reference clustering measures flanking context rather than IS identity.**  
   Using the script’s actual MMseqs options:
   - Identical 1200-bp IS sequences with unrelated 2000-bp flanks → **two clusters**.
   - IS sequences differing at 10% of bases, with identical flanks → **one 95% cluster**.

   Verified with MMseqs in `test_pipeline.py`. Expected family membership should follow the IS sequence. **Fix:** cluster defined element/core sequences; retain flanks as annotations and searchable context.

7. **High — `scripts/search_refs_in_logan.sh:95`: flank-only hits count as successful IS validation.**  
   Input: a 5200-bp reference comprising flank/core/flank; its only hit is Q`[0,2000)`, entirely within the left flank. Actual: `References with >=95% id hit: 1 (100.0%)`. Expected: zero validated IS detections. Verified by executing the embedded analysis in `test_pipeline.py`. **Fix:** record the oriented core’s query coordinates and require sufficient core coverage at a coherent locus. Full-reference coverage buckets do not establish core detection.

8. **High — `scripts/find_insertion_source.py:224`: gap-excluded identity is combined with gap-inclusive query coverage, creating false sources.**  
   Input: a 2000-bp inserted fragment and another contig retaining only 1400 bp, with an internal 600-bp deletion. Actual: minimap2 reports a gapped full-span row; the script reports a source at **100% identity and 100% coverage**, despite only 70% being present. Expected: rejection under the default 80% coverage threshold. Verified end to end: `n_src_hits=1`, source `SRR300_2:200-1600(+,100%)`. **Fix:** compute source coverage from aligned query bases, excluding CIGAR `I` operations.

9. **High — restart handling reuses incomplete or obsolete results, defeating the fixes.**  
   Relevant locations: `scripts/search_logan.sh:70`, `:95`; `scripts/hmm_ncbi_ecoli.sh:42`, `:90`, `:120`; `scripts/hmm_logan_pipeline.sh:147`. Verified:
   - Header-only interrupted `.domtbl` → `SKIP`.
   - Nonempty partial PAF → `SKIP`.
   - New one-assembly manifest plus stale GCA/GCF protein files → two proteins concatenated for one selected assembly.
   - Failed-run chunk directory reused with fewer chunks → **five stored records for two input records**.

   Expected: only complete products from the current inputs/options. **Fix:** use input/options manifests, completion markers and atomic publication; concatenate only manifest-listed outputs; create a fresh chunk directory per attempt.

10. **Medium — `scripts/run_minimap2.sh:89`, `scripts/hmm_search_is110.sh:242`: concatenation can erase genome identity.**  
    Input: two different genome records both named `ctg`, one containing the left IS half and the other the right half. Actual: the detector joins their PAF rows into an insertion; my minimap2→detector experiment reported a 599-bp gap. Expected: no same-contig event. **Fix:** prefix target sequence IDs with assembly/sample identity before concatenation and preserve a reverse mapping. This matters whenever input headers are not globally unique.

11. **Medium — `scripts/find_fragment_origin_local.sh:49`, `:81`, `:95`: the origin summary is neither a full-length-match census nor evidence of origin.**  
    Input: 20 exact copies in one dataset and a 300-bp shared segment of a 2000-bp fragment in another. Actual: counts **6 and 1**; the 15%-coverage match contributes to organism distribution. Verified in `test_pipeline.py`. Self/known insertion hits remain eligible, equal identities select whichever organism is encountered first, and queries without hits disappear. **Fix:** require aligned query coverage, exclude known insertion/self loci where claiming sources, report tied best matches and explicit no-hit rows, and label capped counts as capped.

12. **Medium — `scripts/find_lineage.py:163`: clustering changes when input order changes.**  
    Input: a 1000-bp sequence contained in a 2500-bp sequence. Actual:
    ```
    lengths [1000,2500] → clusters [0,0]
    lengths [2500,1000] → clusters [0,1]
    ```
    Verified using actual minimap2 in `test_lineage.py`. `ava-ont` suppresses one query–target direction, while the coverage check examines only the query. **Fix:** use a symmetric coverage criterion; use both-sequence coverage for fragment equivalence. Also set `check=True`: a simulated minimap2 failure currently returns ordinary singleton assignments.

13. **Medium — `scripts/paf_utils.py:27`: subtracting full CIGAR gaps from column 11 can produce identity above 100%.**  
    Input: valid minimap2 counter conventions for `1000=600D1000=`, with 600 ambiguous bases in the gap: matches=2000, block=2000, `nn:i:600`. Actual: **142.857%**; expected: 100%. Verified in `test_misc.py`. Minimap2 excludes ambiguous bases from `blen`, including ambiguous gap bases, before printing column 11. [Minimap2 2.30 implementation](https://github.com/lh3/minimap2/blob/v2.30/align.c#L242-L290). **Fix:** derive the aligned-column denominator from `M/=/X` operations, with an explicit ambiguity policy; do not infer it by subtracting full gaps from `blen`.

14. **Medium — identity thresholds remain inconsistent between representations and readers.**  
    Locations: `scripts/paf_utils.py:30`, `scripts/run_minimap2.sh:107`, `scripts/run_cds_minimap2.sh:54`, `scripts/parse_paf.py:59`. Verified:
    - A 95%-aligned-identity row → 95% with CIGAR, **94.95%** with its gap-compressed `de` value.
    - A perfect 600-bp-gapped alignment → 100% with CIGAR, **76.923%** without tags.
    - The shell `$10/$11` filter retained one of two actual minimap2 rows, dropping the gapped perfect copy.

    `de` includes gap-opening events and is not the same metric. **Fix:** share one documented identity calculation across readers; mark approximations explicitly or require CIGAR for exact filtering. Exact gap-excluded identity cannot generally be reconstructed from tagless PAF.

15. **Medium — `scripts/detect_is110_variants.py:145`, `:148`, `:203`: raw-gap thresholds lose boundary events and repeat-associated insertions.**  
    Actual minimap2 experiments:
    - True 500-bp insertion → `(target_gap, query_gap)=(499,-1)` under `asm10`, `(496,-4)` under `asm5`; both rejected.
    - True 200-bp deletion → query gap 197 and target overlap 3; rejected.
    - True 600-bp insertion in an IS with terminal repeats and 60-bp junction homology → rejected on both strands; corresponding 5-bp controls pass.

    **Fix:** normalize overlapping junction alignments, use the net length change with signed gaps, and represent breakpoint uncertainty. Separate permitted overlap from permitted missing query sequence instead of applying `abs(query_gap)` to both.

16. **Medium — `scripts/detect_is110_variants.py:254`, `scripts/extract_insertion_fragments.py:153`: deduplication does not provide stable physical-locus counts.**  
    Verified:
    - Same-event reference junction offsets 0/20/21 → two detector records.
    - Extractor intervals `[1019,1619)` and `[1020,1620)` → two fragments, despite a 1-bp difference.
    - Different event types at identical target junctions remain separate records.

    The extractor’s floor bins do not implement distance tolerance. **Fix:** normalize breakpoint uncertainty into stable locus IDs shared by detector and readers; retain reference/type interpretations as annotations on a physical locus. Do not silently equate record counts with locus counts.

17. **Medium — `scripts/hmm_ncbi_ecoli.sh:67`: assembly selection depends on a particular directory layout and guessed pairing.**  
    Verified in `test_shell.py`:
    - Flat layout with two genomes → the last path printed twice; first genome lost.
    - Nested common `genome/` directories → the same failure.
    - GCA `.1` and corresponding-number GCF `.2` → both retained.

    The worker also derives its output accession from the immediate parent, creating collisions in those layouts. **Fix:** use an explicit assembly/path manifest and NCBI pairing metadata, with a declared version policy. Exact same-number/version GCA/GCF pairs and GCF-only assemblies worked.

18. **Medium — legacy split callers still accept inversions and discard loci.**  
    Locations: `scripts/hmm_search_is110.sh:290`, `scripts/run_cds_minimap2.sh:111`, `scripts/search_split_vs_gtdb.sh:108`, `scripts/detect_split_all.py:100`. An inverted second half produced **one split** in each embedded shell detector. Separately, two legitimate loci on different contigs produced **one row** in `detect_split_all.py`, because deduplication keys only on query. Verified in `test_legacy.py`. **Fix:** route callers through one tested strand/order/junction implementation and deduplicate by physical locus.

19. **Medium — `scripts/find_insertion_source.py:199`, `:294`: rows disappear, and empty input crashes.**  
    Input: two insertion rows in one sample, one referring to a missing contig. Actual: exit 0 with one output row; expected: two rows, including an explicit error. Header-only insertion input produces `ZeroDivisionError`. Verified in `test_lineage.py` and `test_core.py`. **Fix:** emit one result/status per input locus and guard zero denominators.

20. **Medium — `scripts/hmm_logan_pipeline.sh:89`, `:180`: legitimate zero-result stages abort under `set -e`.**  
    `N=$(grep -c '^>' empty_output)` prints zero but returns status 1; the assignment aborts the script. Verified in `test_shell.py`. Expected: a defined empty-result outcome with current header-only tables. **Fix:** count with a command that succeeds on zero matches, then explicitly handle the empty case. The `grep` substitutions inside `echo` did **not** abort in my test; they should not be confused with assignments.

21. **Medium — `scripts/search_logan.sh:118` versus `scripts/find_fragment_origin.sh:26` and `scripts/run_find_source.sh:24`: default producer and consumers use different result names.**  
    The search writes `is110_ecoli_insertions.tsv`; consumers read `is110_ecoli_95_insertions.tsv`. Verified that the producer file exists while the consumer path does not. Running the separate `detect_95.sh` supplies the latter, but direct search→origin/source execution fails or reads an old result. **Fix:** unify prefixes or pass the exact produced TSV through the workflow.

22. **Medium — `scripts/parse_paf.py:136`, `scripts/plot_synteny.py:210`: multi-query results are combined as one query.**  
    Two full copies yielded **200% query coverage** in the parser summary. Plotting target hits from queries of lengths 2000 and 3000 against the 2000-bp query track raised:
    ```
    LinkRangeError ... upper_end=3000 ... seg_size=2000
    ```
    Verified with the installed plotting library. **Fix:** group/filter by query and target; calculate coverage from merged intervals per query/locus, and keep copy burden separate.

23. **Medium — hard-coded foreign-cluster paths make several entry points unusable here.**  
    Examples: `scripts/detect_split_all.py:11`, `scripts/hmm_search_is110.sh:18`, `scripts/run_cds_minimap2.sh:15`, `scripts/search_split_vs_gtdb.sh:12`. Running the standalone detector produced `FileNotFoundError` for `/groups/rubin/...`. **Fix:** accept paths as arguments/configuration and pass them into embedded Python; changing only the shell variable is insufficient where Python repeats the constant.

24. **Low — `scripts/extract_insertion_fragments.py:150`, `:159`: clamped sequence and recorded coordinates disagree.**  
    Input: a 3000-bp contig with padded requested interval `[1000,3600)`. Actual: 2000 bases written under header `__1000_3600`; expected end coordinate 3000. Verified in `test_followups.py`. **Fix:** clamp the stored coordinates before constructing the deduplication key and header.

25. **Low — `scripts/blast_fragments.sh:71`: multiline error messages corrupt the TSV.**  
    A simulated failed BLAST with newline/tab-bearing stderr produced field counts **`[8,3,7]`** instead of a header and one eight-field error row. Verified in `test_misc.py`, without making a remote request. **Fix:** normalize stderr whitespace before writing, or use a consistently quoted TSV writer and matching reader.

The minimal corrected locus-chain and breakpoint-aware nesting functions are in [proposed.py](<tmp>/proposed.py). They were tested on two intact copies, two interrupted copies, opposite strands, minus-strand events and reversed target order. They must receive CIGAR-normalized blocks; ambiguous competing chains still need explicit handling.

For the remaining classification and clustering changes, the essential corrections are:

```python
# Intact classification: use aligned-block coverage, not PAF query span.
status = (
    "variant" if n_real
    else "intact" if cov >= 95 and not events
    else "partial_or_unresolved"
)

# Preserve the existing "70% of the smaller sequence" clustering policy,
# but make it independent of the emitted query/target direction.
paired = sum(int(n) for n, op in _CG.findall(cigar_of(c))
             if op in "M=X")
if paired / min(int(c[1]), int(c[6])) < 0.7:
    continue
```

For nesting, preserve query breakpoints from `extract_insertions()`, cluster across **all** levels once, and apply `nesting_edge()` from the linked file. Use stronger both-sequence coverage for inherited-fragment equivalence, and require a connected path of supported edges.

**(a) Verdict on the HEAD commit message**

| Commit item | Verdict |
|---|---|
| Aligned identity and within-row indels | **Incomplete** — ordinary ACGT CIGAR arithmetic works; ambiguity handling and multi-indel coverage fail. |
| Minimap2 census/splitting options | **Incomplete** — recover split arms, but introduce sensitivity gaps and do not guarantee every copy. |
| Sorted, single-writer Logan concatenation | **Correct** — verified byte-for-byte, including paths containing spaces. Existing cached files remain a problem. |
| Variant strand/order/junction tests; per-locus output; always write outputs | **Incomplete** — ordinary controls pass; cross-locus calls and deduplication failures remain. Empty outputs are correctly overwritten. |
| `detect_split_all.py` strand/order/distance | **Incomplete** — those changes work; per-query deduplication still loses loci. |
| Fragment deduplication and end clamping | **Incomplete** — slicing clamps correctly; bin deduplication and coordinate metadata fail. |
| Source minimap2 failures and census options | **Incomplete** — failures now raise; coverage creates new false positives. |
| HMM gathering thresholds, chunk failure handling, deterministic ordering | **Incomplete** — active profiles support GA and waits detect failure; restart products remain unsafe. |
| One assembly copy; recorded/fatal failures | **Incomplete** — normal layout works; layout, pairing and stale-output cases do not. |
| Fragment clustering at 95%, both coverages ≥90% | **Correct** — a contained 300-bp/5000-bp pair changed from one representative to two. |
| Failed remote BLAST becomes `ERROR` | **Incomplete** — classification is correct; error serialization can corrupt the table. |

The most consequential “Not changed” problems are **per-contig multisplit chaining, unsupported nesting, and missing contig-cluster expansion**. They directly undermine copy counts and lineage evidence. Flank-dominated reference clustering and flank-only validation independently undermine reference coverage. All were reproduced above, with concrete fixes. The cross-organism summary cannot establish origin in its present form.

**(b) Checked and found correct**

- `D` means extra target/genome sequence; `I` means query sequence missing from the genome. `large_indels()` positions and `row_event()` left/right assignment are correct for single events on both strands. Walking successive indels produces correct gap positions.
- Ordinary plus/minus insertions, minus-strand deletions, 5-bp query overlaps and events whose visible right arm ends exactly at the contig end pass. Opposite-strand/swapped-order controls are rejected by the main detector.
- Both active Pfam files are valid HMMs with GA thresholds; `hmmsearch --cut_ga` succeeded. The two named alias files fail `hmmstat`, but current scripts reference the valid `PF01548.hmm`/`PF02371.hmm` files.
- Quoted `xargs -I` handling preserved paths containing spaces. Two hundred parallel failure appends produced 200 unique complete records. I did not reproduce a single-run append race.
- Chunk PID waits detect a failed child. Bash array/PID arithmetic tested successfully.
- Reference extraction’s `start-1:end` conversion and strand orientation are correct. Variant TSV headers and rows have the expected 20 columns, matching the Python readers.

For minimap2 specifically, both presets split a 199-bp indel under the new options. `-N` limits secondaries and `-p` sets their score ratio; neither overrides seeding/alignment thresholds. Forced multi-part stdin indexing with `-N 1` emitted **67 rows and 34 primaries**; `--split-prefix` reduced this to two rows and one primary. These are per-part reporting effects, not 34 independent globally selected primaries. [Minimap2 manual](https://lh3.github.io/minimap2/minimap2.html).

A 3000-copy synthetic target emitted 3000 rows for one query. That demonstrates output scaling; I did **not** establish an actual out-of-memory failure at production scale.

**(c) Could not verify**

- Production event frequencies, reported assembly totals, actual pairing metadata, and whether existing cached outputs are corrupted; no real sequence files or large scratch data were opened.
- Production-scale runtime/memory, SLURM behavior, remote BLAST availability, or concurrent independent jobs sharing an output directory.
- Biological ancestry or insertion origin from these alignments alone.

The three highest-value tests to add are:

1. **Structural-event round trip:** synthetic sequences → minimap2 → detectors → extraction/source readers, covering both strands, repeats, overlaps, threshold boundaries, several indels, two loci and non-events.
2. **Census and failed-resume invariance:** assembly layouts/versions, clustered sample expansion, stale GCA files, obsolete chunks, interrupted outputs and zero results.
3. **Lineage/validation negative controls:** unrelated 0/1/2-split observations must not form an inherited ladder; sequence permutation must not change membership; flank-only hits must not validate an IS.
---

## What was done with these findings (2026-11-17)

Tests: `python3 tests/test_fragment.py` (synthetic PAF rows, random DNA, and a
real minimap2 → detector round trip with both presets).

| # | Status | Change |
|---|---|---|
| 1 | fixed | `detect_is110_variants.build_chains`: each aligned block gets at most one partner on each side, the nearest; the number of pairings refused for that reason is printed |
| 2 | fixed | rows are cut into aligned blocks (`paf_utils.aligned_blocks`); coverage = union of the chain's aligned query intervals |
| 3 | fixed | `detect_multisplit.py` reports locus chains on both strands; `status` = intact / variant / partial; partial chains are not counted as intact |
| 4 | fixed | `find_lineage.py`: insertion-only ladders, one clustering over all levels, nesting test (`supported_depth`, `supported_edges`, `parent_example`) |
| 5 | changed default | `hmm_logan_pipeline.sh`: contigs are no longer dereplicated (`DEREP_CONTIGS=0`); `DEREP_CONTIGS=1` restores the old behaviour with a warning |
| 6 | fixed | `extract_is110_refs.sh` clusters the transposase gene (`is110_core.fa`); representatives are the flanked records of the representative genes |
| 7 | fixed | `ref_core_coverage.py`: detection needs ≥ 90% of the core in one row at ≥ 95% identity; flank-only hits are reported as such |
| 8 | fixed | source coverage from aligned query bases |
| 9 | fixed | `stage_lib.sh` markers (inputs + options) and tmp-then-rename in `search_logan.sh`, `hmm_logan_pipeline.sh`, `hmm_ncbi_ecoli.sh`, `search_refs_in_logan.sh`; fresh chunk directory; manifest-only concatenation |
| 10 | check added | `run_minimap2.sh`, `hmm_search_is110.sh` stop when a sequence name occurs in two genome files (names are not rewritten) |
| 11 | fixed | `fragment_origin_summary.py`: coverage filter, own locus excluded, ties named, no-hit rows, capped counts marked |
| 12 | fixed | symmetric coverage; minimap2 failure raises |
| 13 | fixed | identity = matches / (M, =, X columns), capped at 100 |
| 14 | partly | one calculation shared (`paf_utils.py filter` replaces the awk filters; `parse_paf.py` uses the CIGAR). Without a CIGAR the value stays an approximation; the detector and the filter print how many rows that concerns |
| 15 | fixed | net length change with signed gaps; `--max-junction-overlap` |
| 16 | fixed | single-linkage locus groups across references and event types (`--locus-tol`); extractor uses the same tolerance |
| 17 | fixed | accession taken from anywhere in the path; one assembly per number, GCF preferred, highest version; explicit manifest |
| 18 | fixed | the three embedded detectors call `detect_split_all.py`, which uses the tested chain code and reports one row per locus |
| 19 | fixed | one output row per input locus (`no_contig`, `no_seq`); empty input no longer divides by zero |
| 20 | fixed | `count_fasta` / `count_rows`; zero-result stages write a header-only table |
| 21 | fixed | `search_logan.sh` writes `is110_<org>_95_*`, the name the next steps read |
| 22 | fixed | summary per (genome, query) with merged-interval coverage; `plot_synteny.py --query` |
| 23 | partly | paths can be set by environment (`GTDB_FNA`, `OUTDIR`, …, `IS110_FRAGMENT_SCRIPTS`); `detect_split_all.py` takes arguments. The `#SBATCH --output` lines still name the first cluster |
| 24 | fixed | clamped coordinates in key and header |
| 25 | fixed | error text flattened to one field |

Behaviour that changes results on a rerun: event sizes are net lengths, so
insertions/deletions at overlapping junctions are now called; one record per
physical locus; multisplit rows are per locus and partial copies are no longer
"intact"; the Logan HMM census searches every contig; reference clusters
follow the transposase gene.

Not done: competing pairings are counted but not written per record; breakpoint
uncertainty is not stored (raw `target_gap` / `query_gap` stay signed, the net
size is their difference); `find_lineage.py` still examines only the first
`--max-ins-examples` observations per level of the 30 deepest candidates;
`empty_filled_sites/` was not reviewed.
