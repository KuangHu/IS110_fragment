# Adversarial review of the core pipeline, round 2 (commit 217aea1)

Reviewer: OpenAI codex CLI, 2026-11-17, on `git diff 3a2e9b2 217aea1`, focus: pairing by locus.
Fixes: the commit that adds this file. Not fixed, documented in README "Known limits": #9 (projection through deleted anchor bases), residual part of #3 and #5.

I found **12 reproducible issues**. The supplied suite passes. I reviewed commit `217aea1` from a temporary snapshot and left the repository unchanged.

All probe scripts, generated FASTAs/PAFs, and logs are under `<tmp>/`. Coordinates below are 0-based, half-open unless stated otherwise.

1. **High — `scripts/lib_alleles.py:240`: mutual-nearest pairing still confidently pairs neighbouring copies.**

   Choosing the smallest signed gap does not establish that two anchors belong to the same locus.

   **Failing input:** Actual minimap2 alignments contained upstream hits at `[1100,1500)` and `[3000,3500)`, and downstream hits at `[4700,5200)` and `[6700,7100)`. The outer partners were clipped to 400/500 query bases. The two intended same-copy gaps were both 3,200 bp.

   **Actual:** The second upstream anchor paired with the first downstream anchor:

   ```text
   start=3500 end=4700 gap=1200
   score=2.0 second_score=1.8 margin=0.10 status=resolved
   ```

   That manufactures a 2,000 bp removal relative to the filled state. **Expected:** retain uncertainty between repeat copies; do not resolve this cross-copy combination as an empty allele.

   **Verified:** `align_probes.py`, fixture `nearby`; actual PAF is `align/nearby/hits.paf`.

   **Suggested fix:** Require additional locus compatibility or synteny evidence for interleaving anchors. Mutual-nearest status alone cannot certify orthology.

2. **High — `scripts/build_records.py:122`: negative gaps become falsely exact, truncated element calls.**

   `_empty_interval` converts every nonpositive distance to `""`. Decomposition then treats the entire source interval between anchors as the insertion, even when that interval excludes element ends.

   **Failing input:** A 3,240 bp element contains a 1,200 bp transposase and 1,020 bp on either side. At the default `D=1000`, each anchor extends 20 bp into the element. Actual minimap2 alignments against the empty genome project to:

   ```text
   gap_start=2020 gap_end=1980 gap=-40
   ```

   **Actual full pipeline:**

   ```text
   is_element.length=3200
   is_element_length_inferred=3240
   method=empty_vs_filled_decomposition
   consensus_support=1.0
   junction_microhomology_bp=0
   ```

   Forty element bases are omitted despite a confident decomposition label. **Expected:** recover all 3,240 bp using additional flanking sequence, or refuse this distance and try larger anchors.

   **Verified:** `end_to_end_probes.py`, fixture `inside_IS_D1000`. Stages 4 and 5 both classify the negative gap as empty; Stage 6 introduces the truncation.

   **Suggested fix:** Do not use an empty string as exact boundary evidence for overlapping projected anchors. Extend the comparison interval or retry a larger D.

3. **High — `scripts/lib_alleles.py:445`: the clean-junction rule still accepts flank indels as insertion sequence.**

   A shifted flank need not produce mismatches next to the junction: periodic or homopolymeric sequence can conceal the shift.

   **Failing input:** A 1,500 bp insertion, a 3 bp empty-allele deletion 45 bp upstream, and a 40 bp homopolymer immediately before the junction.

   **Actual full pipeline:**

   ```text
   length=1503
   source_start=5348 source_end=6850
   consensus_method=hamming_split
   consensus_support=1.0
   junction_microhomology_bp=0
   target_lost_bp=0
   ```

   **Expected:** the 1,500 bp insertion at source `[5350,6850)`, with the separate flank difference accounted for, or refusal.

   An indel-aware alignment of these same alleles returns:

   ```text
   method=tolerant_alignment insert_start=4850 insert_len=1500 identity=1.0
   ```

   **Verified:** `hamming_pipeline_probes.py`, fixture `concealed_indel`.

   **Suggested fix:** Validate proposed splits and suspicious flank differences with an indel-aware alignment. Mismatch-free junction windows do not prove indel-free flanks.

4. **High — `scripts/lib_alleles.py:253`: the stray rule discards an equally strong competing placement.**

   The XOR condition excludes pairs whose **two** anchors are outside the selected mutual-nearest loci.

   **Failing PAF geometry:** All hits belong to one assembly. Full hits score 1 per anchor; clipped hits score 0.8.

   ```text
   contig x: full up inner edge 1500; full down inner edge 4700
   contig y: full up edge 1500; clipped up edge 3500
             clipped down edge 4700; full down edge 6700
   ```

   Contig y’s internal clipped pair becomes its mutual-nearest locus. Its two full outer anchors form another valid, score-2.0 placement, but that pair is excluded because neither endpoint belongs to a retained locus.

   **Actual:**

   ```text
   winner=x score=2.0 second_score=1.8 margin=0.10 status=resolved
   ```

   **Expected:** ambiguous: x and y both provide compatible score-2.0 complete placements.

   **Verified:** `probes.py`, output `hidden_stray_pair`.

   **Suggested fix:** Preserve competing complete pairings during ambiguity assessment, including pairs with two unmatched endpoints.

5. **High — `scripts/pipeline.py:201`: changing the genome database still silently reuses stale results.**

   Command signatures contain paths, not input identities. Stage 3 supplies no `inputs`; Stage 2 supplies only `is_hits.tsv`. Database and code changes therefore escape validation. Stage 1 also remains cached solely on output existence.

   **Failing sequence:** Complete the synthetic pipeline, remove one empty genome from the FASTA at the same path, rebuild its index, and rerun the identical command.

   ```text
   initial:           IS1=2 empties, IS2=2 empties
   changed_db_cached: IS1=2 empties, IS2=2 empties; stages executed=0
   changed_db_fresh:  IS1=1 empty,   IS2=1 empty
   ```

   **Expected:** invalidate the alignment and dependent results when the database changes.

   **Verified:** `cache_probes.py`, using actual pipeline runs.

   **Suggested fix:** Record dependencies for every stage, including source/target FASTAs, indexes, profiles, scripts and relevant parameters. Give automatically generated Stage 1 results a manifest too.

6. **High — `scripts/is_detect.py:248`: the new HMM disambiguation introduces a fresh cache-name collision.**

   Appending the profile’s position does not create a collision-free namespace.

   **Failing input, fresh working directory:**

   ```text
   A/model.hmm       matches protein P
   B/model.hmm       matches protein P
   C/model.1.hmm     independently matches nothing
   ```

   The first profile is renamed `model.1`, colliding with the third profile’s original basename.

   **Actual ALL-domains output:**

   ```text
   P retained; domains_hit=model.1,model.2,model.1
   ```

   **Expected:** no retained protein, because C does not match.

   Changing a single profile from `A/model.hmm` to a nonmatching `D/model.hmm` in an existing working directory also reuses the old search.

   **Verified:** `hmm_probes.py`, using real `hmmbuild` and `hmmsearch`; see `collision_three` and `third_independent`.

   **Suggested fix:** Use cache identities derived from profile contents, searched proteins and search parameters. Keep display names separate.

7. **High — `scripts/build_records.py:438`: flank-deletion rejection examines only five observations but treats the remainder as empty evidence.**

   **Failing input:** Six distinct target assemblies each delete the same 1,500 bp upstream flank segment and retain the entire transposase.

   **Actual full pipeline:**

   ```text
   n_empty_rejected_flank_deletion=5
   n_empty_observations=1
   method=peak_symmetric_fallback
   is_element.length=1500
   ```

   **Expected:** six rejected flank deletions, zero empty observations, and no inferred 1,500 bp element.

   **Verified:** `end_to_end_probes.py`, fixture `six_flank_deletions_pipeline`.

   **Suggested fix:** Separate boundary-vote sampling from empty-site validation. Unexamined distance candidates must not remain confirmed empty observations.

8. **High — `scripts/pipeline.py:80`: restarting with `--trust-existing` can adopt outputs invalidated by the preceding process.**

   `_RAN` preserves invalidation only in memory.

   **Failing sequence:**

   1. Old upstream and downstream outputs exist.
   2. Upstream successfully regenerates `NEWup`, invalidating downstream markers.
   3. Downstream fails while its old output remains.
   4. Restart in a new process with `--trust-existing`.

   **Actual:**

   ```text
   NEWup OLDdown
   === down: adopting existing outputs (--trust-existing) ===
   ```

   **Expected:** regenerate the invalidated downstream output.

   **Verified:** `cache_probes.py`, output `restart_trust_existing`; the restart uses a separate Python process.

   **Suggested fix:** Persist invalidation state, or quarantine invalidated outputs so legacy adoption cannot certify them after restart.

9. **Medium — `scripts/lib_alleles.py:215`: projecting clipped query bases assumes sequence that may actually be deleted.**

   The strand formulas are arithmetically correct for conserved, unaligned sequence. They are unjustified when clipping reflects an indel.

   **Failing input:** A true 1,500 bp element removal at `D=5000`, plus deletion of 80 bp from each anchor’s inward end. Actual minimap2 hits pass the identity and coverage filters.

   ```text
   physical target inter-anchor sequence=9700 bp
   projected gap=9540 bp
   ```

   **Actual full pipeline:** a 1,660 bp `peak_symmetric_fallback` element instead of the known 1,500 bp element.

   **Expected:** establish the inward-edge correspondence or mark the gap uncertain.

   **Verified:** `align_probes.py` and `end_to_end_probes.py`, fixture `anchor_deleted_D5000`.

   **Suggested fix:** Require aligned inward edges or validate an extension through the clipped region before assigning exact projected coordinates.

10. **Medium — `scripts/lib_alleles.py:240`: equal-distance nearest choices depend on PAF row order and can change record inclusion.**

    Nearest-partner ties retain the first encountered hit before deterministic candidate sorting runs.

    **Failing input:** One full upstream hit; two downstream representations with the same projected inner edge, one full-length and one clipped to 400/500 bases.

    ```text
    clipped first -> score=1.8 second=2.0 margin=-0.111 status=ambiguous
    full first    -> score=2.0 second=1.8 margin= 0.100 status=resolved
    ```

    Stage 6 excludes the first ordering and accepts the second.

    **Expected:** identical results under input permutation.

    **Verified:** `probes.py`, output `tie_order_changes_acceptance`. Equal-score placements on different contigs also change the reported winning contig under reversal.

    **Suggested fix:** Define deterministic nearest-partner tie-breaking before filtering loci, using score and a complete stable hit key.

11. **Medium — `scripts/lib_alleles.py:153`: duplicate PAF hits become fictitious competing loci.**

    Hits are retained as separate list entries; locus membership uses their indices.

    **Failing input:** One upstream/downstream pair, then the identical upstream PAF row repeated.

    ```text
    original:      n_pairs=1 status=unique
    duplicate up: n_pairs=2 status=ambiguous
    ```

    Repeating both rows gives three candidate pairs and remains ambiguous.

    **Expected:** one biological placement regardless of repeated identical alignment records.

    **Verified:** `probes.py`, outputs `single`, `dupe_up`, and `dupe_both`.

    **Suggested fix:** Deduplicate alignment identities before pairing and collapse redundant representations of the same physical placement.

12. **Low — `tests/check_synthetic.py:64`: the strengthened checker still accepts completely corrupted element sequence.**

    Its sequence validation checks length only.

    **Failing input:** Replace IS1’s JSON element sequence with an equally long all-A string.

    **Actual:**

    ```text
    exit=0
    ALL CHECKS PASSED
    ```

    **Expected:** failure because the sequence disagrees with its source interval and emitted FASTA.

    **Verified:** `coverage_probes.py`; output is preserved in `coverage/checker.log`.

    **Suggested fix:** Compare sequences with the strand-oriented source interval and corresponding FASTA record, and validate coordinate/offset relationships directly.

**(a) Verdicts on the 21 round-1 findings**

Here, **a = correct**, **b = incomplete**, **c = introduced a new error**.

| Round 1 | Finding | Verdict | Assessment |
|---:|---|:---:|---|
| 1 | Clipped-anchor gap measurement | **b** | Conserved-clipping formulas work; deleted clipped sequence remains unsafe. |
| 2 | Cross-copy pairing | **c** | Distant-copy fixture improves; filtering introduces hidden competing placements and still resolves nearby cross-copy pairs. |
| 3 | Hamming acceptance of nearby indels | **b** | Random-flank example is refused; a homopolymer conceals the same error. |
| 4 | Flank deletions counted as empties | **b** | Rejection works only for examined observations. |
| 5 | Decomposition excludes ORF bases | **b** | Up to 10 bp per end remains permitted; clipping is now reported and documented. |
| 6 | Hash-dependent boundary consensus | **a** | Ordering and consensus ties are deterministic in the tested cases. Pairing has a separate input-order defect. |
| 7 | Markers omit input/parameter identities | **b** | Command changes are detected; database, profile and code identities remain incomplete. |
| 8 | `--trust-existing` bypasses invalidation | **b** | Same-process case fixed; restart case remains unsafe. |
| 9 | Separate source DB loses records | **a** | Actual separate-source/target pipeline built both records correctly. |
| 10 | Same-basename HMM cache collision | **c** | Original two-profile case fixed; generated aliases introduce another fresh collision. |
| 11 | Stage 5 size filters ignored by Stage 6 | **a** | Positive minimum and maximum limits propagate and work in actual pipeline runs. |
| 12 | NumPy coordinates break JSON | **a** | Alignment-derived results serialize successfully. |
| 13 | Terminal indels escape dominance test | **a** | Terminal competing blocks are counted and the supplied adversarial shape is refused. |
| 14 | Preferred twin has no usable pair | **a** | Split-contig preferred twin correctly falls through to the other assembly. |
| 15 | Minus-strand fallback offsets | **a** | Odd-length extensions agree with coordinates and sequence orientation. |
| 16 | Flank SNP splits `event_key` | **b** | Documented rather than corrected; exact-context keys still split. |
| 17 | Missing FASTAs do not invalidate cache | **a** | Deleted element FASTA and reference FASTA were regenerated. |
| 18 | Census guidance ignores twins/copies | **a** | README/CLI guidance now requires accounting for target hits and checking the cap flag. Automatic census completion remains manual. |
| 19 | Empty keep-list retains everything | **a** | Actual extraction retained zero references. |
| 20 | Checker accepts corrupt coordinates/duplicates | **b** | Both original corruptions now fail; sequence corruption still passes. |
| 21 | `domains_hit` includes nonmatching profiles | **a** | Independent ANY-domain run lists only the matching profile. |

**(b) Checked and found correct**

- Actual minimap2 tests covered inward, outward and two-ended clipping on both target strands, plus outward clipping at contig ends. Conserved-sequence projections returned the expected edges without off-by-one errors.
- Positive-gap `gap_start`/`gap_end` are in genomic order on both strands; `_empty_interval` correctly reverses minus-strand sequence. The negative-gap failure is described above.
- Distant duplicated loci remain ambiguous in the original cross-copy fixture. Opposite-strand copies do not cross-pair; anchors on different contigs do not form a pair.
- Four `PYTHONHASHSEED` values—0, 1, 7 and 99—produced identical direct-probe outputs with input order held fixed.
- A 100 bp microhomology fixture crossing the ORF start correctly clamps the selected placement to the ORF start on both strands. An explicit 1:1 boundary vote chooses the documented shorter call and reports `consensus_tied=True`, support `0.5`.
- Atomic-output failure preserves the previous complete file and leaves no completion marker; successful replacement removes the temporary file.

The tests have meaningful coverage gaps:

- The new cross-copy test exercises separated loci, missing interleaving anchors.
- The flank-deletion test has one target, missing the five-observation cutoff.
- Marker tests pass `inputs` directly, missing the pipeline callers’ absent dependencies.
- Replacing `consensus_tied` with constant `False` in a temporary mutant still passed every stage unit test.
- A genuine single SNP two bases from an otherwise indel-free junction is now refused. In the tested asymmetric element, that forces a symmetric fallback shifted 100 bp from the true boundaries. This is a documented conservative refusal with a measurable loss of boundary accuracy.

**(c) Could not verify**

- Production prevalence, affected paper records, or error rates: no production database or independently curated boundary truth was supplied.
- Whether minimap2 itself emits identical duplicate rows in the deployment configuration; duplicate-row robustness was tested with explicit PAF inputs.
- Large-scale SLURM interruptions, concurrent writers, filesystem crash durability, and biological HMM sensitivity/specificity.