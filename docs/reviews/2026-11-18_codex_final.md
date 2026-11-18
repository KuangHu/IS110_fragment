**No. I would not start the paper rerun with this HEAD.** Verified failures remove valid lineages, merge distinct rearrangements, retain self-observations for supported non-NCBI names, and allow cached stages to mix old and new counting rules.

The repository remains unchanged. Reproductions used synthetic sequences only. Commands and outputs are collected in evidence.txt; the reproduction runner identifies the individual probes.

1. **High — `scripts/build_lineages.py:312`, triggered by `scripts/find_variants.py:366`: self-exclusion removes the reference step and can erase the entire lineage.**

   `V_REF` is added only when another observation is `clonal_V1`. After excluding the source genome, a known reference element can disappear from its own containment series.

   **Failing input → actual → expected:** source element 6,000 bp; another genome contains that element plus a 250-bp insertion. Stage 8 correctly reports one insertion. Stage 9 produces **zero lineages**, despite receiving the reference sequence through `--records`. Expected: the known 6,000-bp reference and its 6,250-bp variant form a two-step series, without counting the source as an external observation.

   **Verified:** `probe_lineage.py`, using real minimap2:
   ```
   LINEAGES_SELF_EXCLUDED n_refs=0
   LINEAGES_CONTROL n_refs=1 steps=[('v_ref', 6000), ('variant', 6250)]
   ```
   The control restores the former source clonal observation.

   **Fix:** retain the known reference as a comparison step independently of its external observation count.

2. **High — `scripts/collect_confirmed_events.py:181`: overlapping envelopes merge distinct target events.**

   The collector reduces all anchors on one contig to their outer envelope, then groups overlapping envelopes. Actual endpoint correspondence is lost.

   **Failing input → actual → expected:** one inversion has anchors at `[10000,10500)` and `[100000,100500)`; another has anchors at `[40000,40500)` and `[70000,70500)` on the same chromosome. None of their actual anchor intervals coincide. Nevertheless, their envelopes overlap, producing **one target event instead of two**.

   **Verified:** `probe_core.py`:
   ```
   TARGET_NESTED events=1 group_sizes=[2, 2]
   ```

   **Fix:** preserve and compare corresponding anchor endpoints, including orientation. Envelope overlap does not establish event identity; transitive overlap grouping can compound the conflation.

3. **High — `scripts/pipeline.py:73`: upgrading this code can retain old stage-4/5 counts while rebuilding stage 6 with self-exclusion.**

   Cache signatures contain command arguments, but no implementation version. The stage-4/5 commands did not change with self-exclusion; the stage-6 command subsequently changed through `--max-secondary`.

   **Failing input → actual → expected:** reuse a completed output directory after these code updates, with unchanged input files and parameters. Stages 4/5 can remain cached with self-observations, while stage 6 rebuilds without them. Expected: all affected stages recompute under the same counting policy.

   **Verified:** `probe_cache.py` exercises the real `_run_stage` with a temporary worker whose implementation changes:
   ```
   CODE_CHANGE_CACHE stage4_old_count=2 stage6_new_count=1
   ```

   **Fix:** include implementation fingerprints in stage dependencies/signatures. For this regeneration, use fresh output directories or explicitly invalidate every affected stage.

4. **High — `scripts/lib_alleles.py:306`: bare non-NCBI reference IDs do not exclude their source genome.**

   For an unprefixed contig, stage 1 generates `<contig>_<protein-number>`. `source_core()` leaves that suffix attached because it only strips the assembly prefix before a pipe.

   **Failing input → actual → expected:** reference `isolateA_7`, originating from target contig/genome `isolateA`. The computed source core is `isolateA_7`, while the target core is `isolateA`; the self-placement survives. Expected: zero external observations.

   **Verified:** `probe_core.py`:
   ```
   SOURCE_CORE isolateA_7 core=isolateA_7 own_placements=1
   ```
   Standard NCBI IDs and prefixed custom names with additional contig pipes correctly returned zero.

   **Fix:** use explicit source assembly metadata where available; otherwise recover the source contig by removing the generated final protein-number suffix before parsing its assembly.

5. **Medium — `scripts/extract_anchors.py:136`, `scripts/build_records.py:692`: unsupported retention still loses hits before stage 6.**

   This is a **remaining loss path**, rather than a new regression. Stage 6 iterates the reference table, so its new stub cannot rescue hits that stage 2 omitted because no anchor pair fits.

   **Failing input → actual → expected:** a valid 1,200-bp ORF at positions 101–1300 on a 5,000-bp contig, with 500-bp anchors and distances 1000/5000. Extraction returns `None`; no reference-table row or final record follows. Expected under the new retention policy: a `no_observation_tnp_only` record.

   **Verified:** `probe_core.py`:
   ```
   EDGE_LOSS ORF_101_1300=None ORF_2001_3200_built=True
   ```

   Stage 2 reports an aggregate skipped count, but the element remains absent from the catalogue.

   **Fix:** seed record retention from `is_hits.tsv`, carrying an explicit “anchors unavailable” status.

6. **Medium — `scripts/build_records.py:500`: the fallback D can falsely label an externally observed element unsupported.**

   When no distance has empties, the function chooses the first candidate—even if it contains only V1+ placements and another distance contains valid V0 observations.

   **Failing input → actual → expected:** transposase length 1,200; D=1000 has a 10,000-bp gap, while D=5000 has the exact filled gap of 11,200 bp. Actual: `no_observation_tnp_only`, `n_v0=0`; with `--drop-unsupported`, the record disappears. Expected: retain the filled observation at D=5000.

   **Verified:** `probe_core.py` and `probe_main.py`:
   ```
   D_SELECTION_DROP records=0
   D_SELECTION_stage4_D filled_observations=1
   D_SELECTION_stage5_D filled_observations=1
   ```

   **Fix:** prefer a candidate with V0 evidence before falling back to a candidate with neither V0 nor empty evidence.

7. **Medium — `scripts/verify_tandem_enrichment.py:139`: the new method bypasses the default called-boundary filter.**

   Its vocabulary omits `no_observation_tnp_only`. A file containing only that label is therefore treated as an unlabelled legacy file.

   **Failing input → actual → expected:** one sufficiently long unsupported record, with default `called_only=True`. Actual: retained for tandem enrichment. Expected: excluded because no boundary was called.

   **Verified:** `probe_core.py`:
   ```
   TANDEM_CALLED_ONLY retained=1 boundary_not_called=0
   TANDEM_OLD_LABEL retained=0 boundary_not_called=1
   ```
   Only the method label changed between those cases.

   **Fix:** recognize the new method, preferably by detecting the presence of the method field rather than maintaining separate incomplete vocabularies.

8. **Medium — `scripts/build_records.py:695`: `--drop-unsupported` runs before evidence rejection.**

   A record may pass this filter, lose every supposed empty observation during flank-deletion rejection, and still be written as unsupported.

   **Failing input → actual → expected:** the only other genome retains the transposase but loses 1,500 bp of upstream flank. Decomposition correctly rejects that observation. With `--drop-unsupported`, the final unsupported record nevertheless remains.

   **Verified:** `probe_main.py` executes `main()` and the real record builder with a sequential worker adapter and synthetic in-memory sequences:
   ```
   DROP_UNSUPPORTED records=1 methods=['no_observation_tnp_only'] n_empty=[0] n_rejected=[1]
   ```
   Expected: zero records.

   **Fix:** apply the final unsupported filter after decomposition and rejection.

9. **Medium — `scripts/build_records.py:673`: multipart indexes produce false cap flags.**

   Counts are summed across index parts, although `-N` applies within each part.

   **Failing input → actual → expected:** four exact target placements, split across two index parts containing two targets each, with `-N 2`. No part reaches three rows, and all four placements are retained. The aggregate count nevertheless marks the anchor capped.

   **Verified:** real minimap2 in `probe_cap.py`:
   ```
   MULTIPART N=2 index_parts=2 raw_rows=4 flag=yes distinct_targets=4
   MULTIPART N=100 index_parts=2 raw_rows=4 flag=no distinct_targets=4
   MULTIPART_IDENTICAL_TARGETS True
   ```

   **Fix:** determine saturation per index part, or report uncertainty for multipart searches. Aggregate PAF row counts cannot establish that the cap was reached.

10. **Medium — `scripts/collect_confirmed_events.py:176`: target-event identity does not collapse GCA/GCF twins.**

    Grouping uses the literal target assembly name. Representative selection upstream occurs per reference element, so different references can retain different usable twins of the same target genome.

    **Failing input → actual → expected:** identical target anchor coordinates on identical contig names, once under `GCA_2.1` and once under `GCF_2.1`. Actual: two target events. Expected: one under the stated twin convention.

    **Verified:** `probe_core.py`:
    ```
    TARGET_TWINS events=2
    ```

    **Fix:** group using canonical assembly identity and validated contig correspondence.

11. **Medium — `scripts/build_records.py:482`: the stub erases known ambiguous placements.**

    `pair_anchors()` counts ambiguous assemblies, but `pick_best_d()` discards that count when no unambiguous pairs remain and hard-codes zero in the stub.

    **Failing input → actual → expected:** another genome has two equally scoring anchor-pair loci. Actual record: `n_ambiguous_placements=0`. Expected: one ambiguous assembly, distinguishable from no alignment evidence.

    **Verified:** `probe_core.py`:
    ```
    AMBIGUITY pair_anchors=1 stub_n_ambiguous=0
    ```

    **Fix:** preserve ambiguity information from the attempted distances when constructing the stub.

12. **Medium — `scripts/build_records.py:399`, `scripts/lib_alleles.py:152`: legacy no-D anchors are accepted but silently lose observations in stages 5/6.**

    This is a **remaining compatibility gap**. The shared parser stores legacy hits under `D=None`; stage 4 handles that key, while stages 5/6 request numeric distances only.

    **Failing input → actual → expected:** valid legacy `__up`/`__down` alignments at the filled distance, with a single-D reference table. Stage 4 reports one filled observation; stage 5 reports none and stage 6 constructs an unsupported stub. Expected: normalize the unambiguous single distance, or reject the format explicitly.

    **Verified:** `probe_legacy.py`:
    ```
    LEGACY_stage4_legacy filled=1
    LEGACY_stage5_legacy filled=0
    LEGACY_STAGE6 v0_observations=0 stub=True
    ```

    **Fix:** explicitly translate single-D legacy hits or fail on unsupported legacy input. The cap regex itself attributes these names correctly.

13. **Low — `scripts/prepare_review_cases.py:56`: focal element annotations still omit their first base.**

    The hit loader now converts coordinates correctly, but `load_locus()` still returns a 1-based start to callers constructing zero-based features.

    **Failing input → actual → expected:** source element and transposase both at 101–1300. The CDS feature is correctly `[100,1300)`, but the focal repeat feature is `[101,1300)`, one base shorter.

    **Verified:** `probe_aid.py`:
    ```
    AID_FEATURES [('CDS', 100, 1300, 1200), ('repeat_region', 101, 1300, 1199)]
    ```

    **Fix:** convert record-derived starts before slicing or constructing features.

**(a) Commit verdicts**

| Commit | Verdict | Reason |
|---|---|---|
| `76b5de7` | **Incomplete** | Changed hit conversions, structural parsing and distinct IDs work; record-derived focal annotations remain wrong (#13). |
| `da572a1` | **New error** | Event grouping conflates distinct endpoint pairs and misses twin equivalence (#2, #10). |
| `73b2699` | **Correct in checked cases** | Tie/support reporting matches normalized votes and retained evidence. |
| `7ba1aae` | **New error** | Breaks reference-step retention downstream (#1); source exclusion is also incomplete (#4). |
| `d79f8c8` | **New error** | Reader vocabulary regression and filtering/counting failures (#6–8, #11); retention remains incomplete (#5). |
| `d56ae43` | **Correct in checked cases** | Confidence precedence matches its documented rules. |
| `b31eff7` | **New error** | Multipart row aggregation falsely reports cap saturation (#9). |
| `e79bfa1` | **Correct in checked cases** | No remaining hash/completion-order failure demonstrated in the exercised fixtures. |

**(b) Checked and found correct**

- Standard NCBI reference IDs, additional pipes inside contig names, and prefixed non-NCBI names exclude the source core correctly.
- Excluding a second locus—or another element—in the **same source genome is intentional under this commit’s assembly-level policy**. Those placements are lost as intragenomic-copy evidence; these observation lists consequently cannot serve as an intragenomic-copy census.
- The movement tests passed, including twin exclusion and avoiding false movement between existing copies.
- Decomposition normalization combined two equivalent calls into one vote category. A genuine 3-versus-3 tie produced `call_support=3`, `call_n_tied=2`, `call_confidence=tied`.
- Mode tie counting uses the same binning as `length_mode`; flank rejection recomputes the retained empty peak and tie inputs.
- Cap row counting occurs before identity/coverage filtering, and attribution accepts legacy names. Stage 6’s aggregation across **all D values** is explicit; it does not specifically describe saturation at `D_used`.
- Existing tests passed: 17 allele tests, 14 stage-unit tests, eight downstream tests, self-alignment, and movement.
- Hash seed 1 with one worker matched seed 2 with four workers for 40 synthetic stage-1/2/6 worker cases and real stage-8/9 outputs.

**(c) Could not verify**

I did not inspect real sequences or `records.json`, estimate affected paper-data counts, rerun external rearrangement confirmations, or execute the complete end-to-end harness. Determinism checks were targeted, not proof across every stage-1–11 input or rendered figure. The inherited weaknesses documented in earlier reviews remain outside these commit verdicts unless explicitly reproduced above.
---

## Response (2026-11-18)

| Finding | Fix | Test (`tests/test_review_final.py`) |
|---|---|---|
| 1 reference step lost | `build_lineages.py`: V_REF is always a step; `n_observations` counts other genomes only, `observed_in_other_genome` says whether any | `test_reference_step_survives_without_clonal_observation` |
| 2 envelopes merge events | `collect_confirmed_events.py`: anchor placements compared one by one, with orientation; a row joins a group only if it matches the group's first row | `test_target_events_compare_anchor_placements` |
| 3 cache ignores code | `pipeline.py`: stage signature includes an md5 of the stage script and `lib_*.py` | `test_stage_signature_follows_the_code` |
| 4 bare IDs keep own genome | `lib_alleles.source_core` strips the protein number when there is no assembly prefix | `test_source_core_without_assembly_prefix` |
| 5 hits without anchors lost | `build_records.py` seeds from `is_hits.tsv`; method `no_anchors_tnp_only` | `test_hit_without_anchors_gets_a_record` |
| 6 wrong fallback D | smallest D with a filled observation is preferred over the first candidate | `test_filled_observation_at_a_larger_d_is_used` |
| 7 tandem filter vocabulary | any method label marks a labelled file (`verify_tandem_enrichment.py`, `cluster_large_is_elements.py`) | `test_tandem_filter_knows_every_method_label` |
| 8 drop before rejection | `--drop-unsupported` is applied after flank-deletion rejection | `test_drop_unsupported_applies_after_rejection` |
| 9 cap summed over parts | rows are counted per index part (`search_anchors.py` sidecar, `_RunCounter`) | `test_cap_is_judged_per_index_part` |
| 10 twins are two events | grouping key is the assembly core; `NZ_` contig prefix ignored (NC_ chromosomes are not matched) | `test_target_events_compare_anchor_placements` |
| 11 stub drops ambiguity | the stub keeps the largest ambiguous-assembly count over the distances | `test_stub_keeps_ambiguous_placements` |
| 12 legacy anchor names | `lib_alleles.anchor_side`: a name without D answers for a single-D reference, is refused otherwise | `test_legacy_anchor_names_reach_stage_6` |
| 13 focal feature off by one | `prepare_review_cases.load_locus` returns a 0-based start | `test_review_case_locus_is_zero_based` |

Also in this round: stages 4–7 read the anchor PAF one element at a time
(`test_multipart_raw_paf_is_rejected_then_sorted`), and the default
`--max-secondary` is 10,000.

Not changed: finding 1's test exercises the reference-step block only, not a
full stage 8–9 run with minimap2; NC_-named RefSeq chromosomes are still two
target genomes for event grouping.
