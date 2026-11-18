# Final adversarial review (codex), 2026-11-18, and response

Review of HEAD `45f24df` against `696e5e6`. The reviewer used synthetic data
only; names such as `probe_core.py` refer to its temporary reproducers, which
are not kept. Every finding below was reproduced as a test in
`tests/test_fragment.py` before it was fixed.

## Response

| Finding | Fixed in | Test (tests/test_fragment.py) |
|---|---|---|
| 1 repeat pieces inside intact copies | `82e731a` | test_final_review_detector: "repeat pieces inside two intact copies" (real minimap2, both strands) |
| 2 short competing match takes an arm | `82e731a` | "competing 200-bp match", "insert carrying a 200-bp piece of the IS" (real minimap2) |
| 3 small I runs counted as coverage | `82e731a` | "three 100-bp I runs", "full-span row with 300 unaligned bases", "parse_paf summary counts aligned bases" |
| 4 dedup merges distinct junctions | `82e731a` | "deletions at 500 and 510", "junctions 0 and 40 bp are not merged through one at 20" |
| 5 shared segment = same insert | `1e1ed22` | test_lineage: "two 1500-bp inserts sharing 1100 bp", "A~B and B~C" |
| 6 first parent understates depth | `1e1ed22` | test_lineage "lineage depth 3 ... in either order"; test_lineage_cli |
| 7 obsolete representatives kept | `a9b2128` | test_final_review_shell: "changed positions", "fragment origin: ... new table rebuilds" |
| 8 concatenation erases genome identity | `a9b2128` | "run_cds_minimap2.sh refuses ...", "search_split_vs_gtdb.sh refuses ..." |
| 9 raw target gap cutoff | `82e731a` | "net insertion of 50,000 bp", "net insertion of 50,001 bp is not" |
| 10 restart signature without genome contents | `a9b2128` | "files_fp changes with the file", "hmm_ncbi_ecoli signs the genome files" |
| 11 flanks carry a low-identity core | `949ac7d` | test_final_review_tables: "a core at 85% identity" |
| 12 flat genome layout lost | `a9b2128` | "flat genome layout without / with the census manifest" |
| 13 --max-query-overlap ignored | `82e731a` | "multisplit enforces --max-query-overlap" |
| 14 extractor dedup order dependent | `949ac7d` | "extractor: offsets 0/20/21 give the same two loci in any order" |
| 15 filtering hides the search cap | `949ac7d` | "self + 5 external rows at -N 5 is '5+'" |
| 16 marker without product is done | `a9b2128` | "stage_done needs its products" |
| 17 N operations | `82e731a` | "a CIGAR with an N operation is refused" |

Choices made where a finding needed a definition (to be confirmed):

- **Locus grouping (4, 14):** records of different references are one locus
  when both junction ends lie within `--locus-tol` (20 bp) of the group's
  first junction in coordinate order. Offsets 0 / 20 / 21 are therefore two
  loci (the previous round made them one by unbounded single linkage).
- **Same insert (5):** one alignment covering at least 90% of the longer
  insert at the cluster identity. A fragment contained in a longer insert is
  no longer "the same insert".
- **Contained pieces (1):** an alignment piece lying within another alignment
  of the same reference on the same contig (5-bp slop, either strand) is not
  used for chaining.
- **Replacement junctions:** pairs with a large gap on the IS and on the
  contig are neither scenario; they are counted in the detector log, not
  reported.
- **N operations (17):** refused with an error.
- **Finding 10:** the signature uses size and modification time of each
  genome file, not a checksum of 7,700 genomes.
- `detect_multisplit.py --max-query-overlap` (default 0.3) now takes effect;
  the detector keeps 0.5 (`--max-arm-overlap-frac`).

Not measured: how often any of these cases occurs in the real data.

## Review text

**No — I would not start the paper rerun with this code.** I verified false insertion calls, missed clean insertions, inflated coverage, incorrect lineage support, and reuse of obsolete representatives. The existing tests pass despite these failures.

Reviewed HEAD `45f24df` against `696e5e6`. The repository is unchanged. All experiments used synthetic data under `<reviewer tmp>/`; reproducer names below refer to that directory.

1. **High — `scripts/detect_is110_variants.py:147`: alternative matches inside intact copies become fictitious interrupted copies.**

   Input: query `R500 + R500 + U800`, where the two 500-bp regions are identical; target contains **two intact copies**, separated by 1000 random bases. Minimap2 reports both full alignments plus shorter repeat alignments.

   **Actual:** one additional insertion per orientation, including junction T`[1300,3099)`, at `55.6%` coverage. **Expected:** zero structural events.

   **Verified:** `probe_non_events.py`, real minimap2 with both pipeline presets, on both strands:
   ```text
   asm10: insertion 1300→3099, coverage 55.6
   asm5:  insertion 1300→3099, coverage 55.6
   ```
   The greedy assignment leaves full-copy alignments and independently chains their contained alternative matches.

   **Suggested fix:** resolve physical-copy membership and contained alternative alignments before assigning event links. A nearest permissible partner alone does not establish a locus.

2. **High — `scripts/detect_is110_variants.py:157`: a short competing match consumes an arm before coverage is evaluated.**

   Input, query length 2000:

   - A: Q`[0,500)` → T`[0,500)`
   - B: Q`[500,700)` → T`[1100,1300)`, a matching segment inside the insert
   - C: Q`[500,2000)` → T`[1800,3300)`

   **Actual:** A–B wins, covers only 35%, and is discarded. C remains unpaired; **zero events** survive. **Expected:** retain the full-coverage A–C insertion, or explicitly report competing interpretations.

   **Verified:** `probe_roundtrip.py` and `probe_non_events.py`. Both presets and strands give zero events; removing only the 200-bp competing alignment restores:
   ```text
   insertion, coverage 100.0
   ```

   **Suggested fix:** evaluate complete chains and their coverage before committing exclusive partners; preserve alternatives when assignment remains ambiguous.

3. **High — `scripts/paf_utils.py:79`, `scripts/detect_is110_variants.py:172`, `scripts/detect_multisplit.py:86`: smaller `I` runs still inflate “aligned” coverage.**

   Runs below the splitting threshold remain inside query intervals, so unioning those intervals counts absent query bases.

   **Verified with real minimap2 under both presets:**

   - A 3000-bp query has 1300 aligned bases, **43.3%** coverage. The first row contains three `100I` runs; the second supplies another 200 aligned bases. **Actual:** an insertion passes the default 50% filter as **53.3%**. **Expected:** reject it.
   - A full-span 3000-bp row with CIGAR `698=100I700=100I700=100I602=` has **90%** aligned coverage. **Actual:** `status=intact`, coverage **100.0**. **Expected:** partial under the default 95% intact threshold.

   Saved inputs/results are in `coverage_real2/` and `small_indel_real/`.

   The related parser problem remains at `scripts/parse_paf.py:137`: `400=600I400=600I400=` reports `total_aligned_bp=2400`, coverage `100.0`, although only 1200 bases align (`probe_core.py`).

   **Suggested fix:** calculate coverage from all `=/X` intervals independently of which indels define structural events.

4. **High — `scripts/detect_is110_variants.py:221`: deduplication merges distinct junctions, including events from the same chain.**

   Input: `500=200I10=200I500=`. There are two deletion junctions, at target positions 500 and 510.

   **Actual:** two detected deletions become **one** output deletion. **Expected:** preserve both supported junctions.

   Single-linkage drift also merges junctions 500 and 540 when another reference contributes an interpretation at 520:
   ```text
   Events: 3 reference x locus records -> 1 physical loci (within 20 bp)
   ```

   **Verified:** `probe_core.py` and `probe_more.py`.

   **Suggested fix:** preserve distinct junctions within an observed chain; match equivalent events across references using bounded breakpoint uncertainty rather than unrestricted single linkage.

5. **High — `scripts/find_lineage.py:170`: fragment clusters do not establish inherited-fragment equivalence.**

   Input: two 1500-bp inserts with an identical 1100-bp internal segment and unrelated 200-bp ends on both sides.

   **Actual:** clusters `[0,0]`; placing these inserts at the same query junction creates an inheritance edge. **Expected:** sharing a local segment should not establish inheritance of the complete insert.

   **Verified:** `probe_lineage.py`, using actual minimap2:
   ```text
   cargo_cluster_lengths [1500,1500,1000]
   shared_bp=1100 clusters=[0,0,1] nesting=0
   ```

   **Suggested fix:** use stronger both-sequence coverage for inheritance, and verify inherited fragments directly. A shared local region—or transitive cluster membership—is insufficient.

6. **High — `scripts/find_lineage.py:319`: choosing the first compatible parent understates the longest supported lineage.**

   Input observations:

   - V1: A
   - V2, in this order: BC, AB
   - V3: ABC

   ABC matches both V2 observations. BC has no V1 parent; AB does.

   **Actual:** `supported_depth=2`. Reordering the two V2 rows changes it to **3**. **Expected:** depth 3 in either order, through A → AB → ABC.

   **Verified:** `probe_lineage.py`, through the real CLI with synthetic compressed contigs and actual clustering:
   ```text
   bad_order:  supported_depth=2
   good_order: supported_depth=3
   ```

   **Suggested fix:** consider every compatible parent and select the maximum reachable depth; retain all supported edges if reporting graph support.

7. **High — `scripts/extract_is110_refs.sh:156`, `scripts/find_fragment_origin.sh:57`: obsolete representatives survive changed inputs.**

   Reference extraction can produce new core sequences, but existing nonempty `is110_reps_95.fa` causes clustering to be skipped. Fragment extraction now refreshes when the insertion table changes, while its clustering step still skips whenever the representative file exists.

   **Actual, verified:**
   ```text
   Step 2 (95% id): reps exist, skipping.
   Step 2: Representatives exist, skipping clustering.
   old_rep_header_retained=True
   ```
   **Expected:** representatives and memberships regenerated from the current cores/fragments.

   **Verified:** `probe_shell.py` and `probe_resumes.py`, executing the actual clustering-control blocks.

   **Suggested fix:** sign inputs, options and dependent outputs throughout these workflows. Invalidate all downstream products when extraction changes. Existing output directories require particular care during this migration.

8. **High — `scripts/run_cds_minimap2.sh:38`, `scripts/search_split_vs_gtdb.sh:50`: two remaining concatenation entry points erase genome identity.**

   Input: two genome files both name their contig `ctg`; one contains only the left query half, the other only the right half.

   **Actual:** both shell scripts exit successfully and report:
   ```text
   ctg: part1_te=1000, part2_ts=7999, gap=6999
   ```
   **Expected:** zero same-contig events.

   **Verified:** actual executions of both scripts with synthetic files and temporary output directories. `probe_shell.py` reproduces the CDS case.

   **Suggested fix:** prefix target IDs with assembly identity, or apply the duplicate-name rejection consistently across every concatenation entry point.

9. **Medium — `scripts/detect_is110_variants.py:150`: the search cutoff still uses raw target gap despite the net-size definition.**

   Input: a 50,030-bp insertion replacing 30 query bases, net change **50,000 bp**.

   Actual minimap2 junctions:
   ```text
   left:  Q[0,1002)    T[0,1002)
   right: Q[1030,2000) T[51030,52000)
   ```
   Thus `tgap=50028`, `qgap=28`, net insertion `50000`.

   **Actual:** zero events because candidate enumeration stops at raw gap >50,000. **Expected:** accept the configured inclusive net upper boundary.

   **Verified:** `probe_more.py`, real minimap2.

   **Suggested fix:** make candidate pruning consistent with the signed-gap/net-size predicate, allowing permitted positive query gaps.

10. **Medium — `scripts/hmm_ncbi_ecoli.sh:79`: the restart signature omits genome contents.**

    Input: replace a synthetic genome at the same manifest path, changing its length from 2000 to 3000 bases.

    **Actual:** both runs produce signature:
    ```text
    13e582d8f6fd6a07181d0b9684a16a85
    ```
    A previous marker therefore skips gene prediction and its downstream searches. **Expected:** changed genome input invalidates the census.

    **Verified:** `probe_shell.py`, executing the actual manifest and signature block.

    **Suggested fix:** incorporate fingerprints of every manifest-listed genome, not just accession/path text.

11. **Medium — `scripts/ref_core_coverage.py:65`: matching flanks can qualify a below-threshold core.**

    Input: 5200-bp reference, core Q`[2000,3200)`, CIGAR `2000=180X3020=`. Core identity is **85%**; row identity is approximately **96.5%**.

    **Actual:** `best_core_cov_pct=100.0`, `status=core_detected` under `--min-identity 95`. **Expected:** this should not count as a ≥95%-identity IS-core validation.

    **Verified:** `probe_refs_origin.py`.

    **Suggested fix:** compute identity over the core’s aligned columns, alongside core coverage.

12. **Medium — `scripts/extract_is110_refs.sh:81`: layouts accepted by the new manifest are silently lost downstream.**

    Input: a valid flat-layout genome named `GCF_000000001.1_genomic.fna`, with two valid position records.

    **Actual:** the manifest accepts the file, but reference extraction requires `GENOME_DIR/<assembly>/` and writes **zero records**, successfully. The accession-directory control writes both records correctly.

    **Verified:** `probe_refs_origin.py`:
    ```text
    extract normal records=2
    extract flat   records=0
    ```

    **Suggested fix:** carry the assembly/path manifest into extraction rather than reconstructing paths from accessions.

13. **Medium — `scripts/detect_multisplit.py:57`, `:75`: `--max-query-overlap` is ignored.**

    Input: arms Q`[0,1000)` and Q`[900,1100)`, with a 600-bp target gap. Query overlap is 50% of the smaller arm.

    **Actual:** both `--max-query-overlap 0` and the default `0.3` report one variant. **Expected:** both reject this overlap under their stated limits.

    **Verified:** `check_controls.py`:
    ```text
    max_overlap 0:   variant, 1 split
    max_overlap 0.3: variant, 1 split
    actual_overlap_fraction=0.5
    ```

    **Suggested fix:** pass and enforce the requested fractional overlap threshold in chain classification.

14. **Medium — `scripts/extract_insertion_fragments.py:163`: extractor deduplication remains input-order dependent.**

    Input: otherwise identical intervals with endpoint offsets 0, 20 and 21.

    **Actual:**
    ```text
    order 0,20,21: 2 fragments
    order 20,0,21: 1 fragment
    ```
    **Expected:** invariant grouping, consistent with the detector’s declared locus policy. Discarded duplicates never enter the extractor’s proximity set, breaking transitivity.

    **Verified:** `check_final.py`, real extraction from a synthetic compressed sample.

    **Suggested fix:** resolve interval groups before writing fragments, using the same corrected locus identities as the detector.

15. **Medium — `scripts/fragment_origin_summary.py:66`: filtering conceals capped searches.**

    Input: one exact self-copy and 20 eligible external copies at 99% identity. Actual `minimap2 -N 5` returns self plus five external hits.

    **Actual:** after excluding self, the summary prints `o_hits=5`, `total_hits=5`, with no lower-bound indicator. **Expected:** capped/unknown completeness, such as `5+`; there are 20 external copies.

    **Verified:** `probe_more.py`, actual minimap2:
    ```text
    raw_rows=6 self_present=True
    o_hits=5 total_hits=5
    ```

    **Suggested fix:** track truncation before identity, coverage and self-hit filtering; propagate uncertainty after filtering.

16. **Medium — `scripts/stage_lib.sh:19`: a matching marker can skip a missing product.**

    Input: valid marker holding the current signature, with its output file absent.

    **Actual:** `stage_done` returns success:
    ```text
    missing_product SKIP product_exists=False
    ```
    **Expected:** rebuild the missing product. Current callers can instead fail later when opening the supposedly completed output.

    **Verified:** `probe_shell.py`, using the shared helper under `set -euo pipefail`.

    **Suggested fix:** completion checks must validate every required output, including companion membership/GFF files; allow legitimately empty products.

17. **Low — `scripts/paf_utils.py:79`: `N` operations are advanced through but never split.**

    Input: `1000=600N1000=`.

    **Actual:** one block spanning T`[0,2600)`, zero events. **Expected:** two aligned blocks separated by the 600 target-only bases, or an explicit rejection of unsupported `N` operations.

    **Verified:** `probe_core.py`; `check_final.py` also checks both strands with mixed `N` and `D`. The unused `large_indels()` helper additionally places a subsequent `D` 600 bp too early because it does not advance through `N`.

    **Suggested fix:** define and enforce `N` handling consistently. This does not affect the prescribed assembly presets, which did not emit `N` in my experiments.

**(a) Commit verdicts**

| Commit | Verdict | Reason |
|---|---|---|
| `3a59f42` | **Incomplete** | Substantial fixes work, but chaining, coverage, nesting and restart correctness remain insufficient; several new failure modes above follow from the replacement implementations. |
| `63ca14c` | **Correct for current callers** | Exact `=/X` identity and rejection of missing/`M` CIGARs work. All repository minimap2 invocations examined supply `-c --eqx`. |
| `39e8cdf` | **Correct within its scope** | The changed log paths and environment reads are consistent. Required variables are exported before both `hmm_search_is110.sh` heredocs. |
| `45f24df` | **Correct as a narrow regression** | The hash-seed test passes. It does not establish invariance to observation/input order or biological correctness. |

**(b) Checked and found correct**

- `tests/test_fragment.py` completed with **ALL PASS**.
- Mixed `I`, `D`, `=` and `X` coordinate walking works on both strands when every indel is cut: 1750 aligned query bases and target endpoint 3050 in both controls.
- The tested 500-bp insertion and 200-bp deletion lower boundaries with junction overlaps pass.
- Ordinary multiple insertion/deletion chains, two interrupted copies, and intact copies on opposite strands pass.
- Same-row/two-row representations and two references for one ordinary junction deduplicate to one result.
- Global lineage clustering now shares IDs across examined levels; ordinary unrelated-fragment controls reject inheritance.
- Reference core extraction matched expected sequences on both strands near contig ends. Representative-selection awk selected the correct flanked record.
- Manifest awk selected the preferred GCF and highest version in flat and nested layouts.
- `DEREP_CONTIGS=0` creates a usable representative symlink.
- A failed HMM stage exited 9, left only its temporary output, and published neither output nor completion marker.
- Source searches reject a candidate retaining only 70% of the insert, retain explicit missing-contig rows, and handle empty input.
- Multi-query plotting selected the requested query and produced PNG/PDF successfully.
- All 18 shell scripts passed `bash -n`.

The test suite misses the competing-block, contained-repeat, short-`I` coverage, adjacent-junction, branching-lineage and downstream-cache cases above. Its hash-seed fixture also contains CIGAR/query-span inconsistencies at lines 62–63, limiting its value as an alignment-validity control.

**(c) What I could not verify**

I did not inspect real sequence files, run the full datasets, touch `/global/scratch`, submit SLURM jobs, or measure how frequently these failures occur in the paper’s inputs. Shell verification used actual temporary executions or extracted script blocks, not a complete production HMM run.

Separate insertion and deletion junctions in one copy worked. A replacement junction with both a large query gap and large target gap remains unclassified by the current predicates; its treatment needs to match the paper’s stated event scope.