# Adversarial review of the core pipeline (stages 1-6), round 1

Reviewer: OpenAI codex CLI, 2026-11-17, against commit 3a2e9b2.
Fixes for these findings: commit 217aea1 (see README "Known limits" for what was documented rather than changed).
Line numbers refer to 3a2e9b2. Absolute paths were rewritten to repo-relative.

Reviewed all scoped files, tests, README, and HEAD (`3a2e9b2`). No repository files changed. The supplied test suite passed; the failures below came from separate adversarial probes. Fixtures and logs are under `<tmp>`.

1. **High — [scripts/lib_alleles.py:212](../../scripts/lib_alleles.py#L212): clipped anchor alignments shift measured gaps.**

   Placement uses alignment endpoints without accounting for unaligned query bases, although 80% coverage is accepted.

   **Failure/verification:** Actual minimap2 alignments of 500-bp anchors covering query intervals `[0,400)` and `[100,500)` made a true 3,200-bp filled gap measure 3,400 bp. Both target orientations were incorrectly classified `v1plus_filled`.

   **Fix:** Require alignment through each anchor’s inward-facing edge, or explicitly map that edge using query coordinates and alignment operations.

2. **High — [scripts/lib_alleles.py:218](../../scripts/lib_alleles.py#L218): scoring can confidently pair different repeat copies.**

   Selecting one highest-scoring combination does not prevent the cross-copy pairing described as fixed in the comments.

   **Failure/verification:** An actual minimap2 fixture contained two loci, each with a 3,200-bp gap. The first locus’s upstream anchor and second locus’s downstream anchor aligned completely; their partners aligned at 80% coverage. `place_pair` selected the cross-copy 12,200-bp gap as `resolved`, margin `0.10`, manufacturing a V1+ allele.

   **Fix:** Add locus compatibility checks, including intervening anchor placements; reject unresolved cross-copy combinations.

3. **High — [scripts/lib_alleles.py:363](../../scripts/lib_alleles.py#L363): `hamming_split` accepts nearby indels and reports incorrect exact boundaries.**

   The mismatch-window heuristic does not establish that the flanks are indel-free.

   **Failure/verification:** With `Random(42)`, a 10,000-bp background, a 1,500-bp insertion at offset 5,000, and an empty allele deleting background bases `[4990,4993)`, Stage 6 returned a **1,503-bp** element at `[4997,6500)`, support `1.0`, microhomology `0`, and target loss `0`.

   **Fix:** Validate the junction and flanks with an indel-aware alignment; separate the small indel from the dominant insertion or refuse the call.

4. **High — [scripts/build_records.py:291](../../scripts/build_records.py#L291): unrelated flank deletions become empty-site evidence.**

   Length classification is retained even when decomposition establishes that the missing sequence does not contain the transposase.

   **Failure/verification:** At D=5,000, deleting 1,500 bp upstream while retaining the entire 1,200-bp transposase produced `n_empty_observations=1`. Decomposition reported `insert_misses_tnp`, but the record still emitted a 1,500-bp `peak_symmetric_fallback` element.

   **Fix:** Reclassify these observations as flank deletions; distinguish unverified distance candidates from confirmed empty sites.

5. **High — [scripts/build_records.py:401](../../scripts/build_records.py#L401): accepted decomposition boundaries can exclude transposase bases.**

   The “whole transposase” check permits ten bases missing at either end. HEAD’s length clamp only protects the fallback.

   **Failure/verification:** For transposase interval `[3000,4200)`, an empty allele removing `[3005,4195)` produced an accepted 1,190-bp element with offsets `+5,-5` and support `1.0`.

   **Fix:** Require the selected placement to contain the complete ORF. Represent partial deletions separately.

6. **High — [scripts/build_records.py:374](../../scripts/build_records.py#L374): boundary consensus depends on Python’s hash seed.**

   Assembly set iteration supplies an unstable order; sorting equal identities preserves it, and truncating to five observations changes the winning vote.

   **Failure/verification:** Six equally scoring empty assemblies supported three 1,400-bp and three 1,600-bp calls. Identical inputs returned **1,600 bp with seed 0** and **1,400 bp with seed 1**, each reporting support `0.6`.

   **Fix:** Use deterministic ordering and define explicit behavior for tied observations and competing boundary calls.

7. **High — [scripts/pipeline.py:56](../../scripts/pipeline.py#L56): completion markers do not identify the inputs or parameters they certify.**

   Changes to thresholds, databases, hits, HMMs, or code can silently reuse previous results. Stage 1 reruns also do not invalidate downstream markers.

   **Failure/verification:** Rerunning the supplied fixture with identity changed from 95% to 100% executed zero stages and retained two empty observations per record. A fresh 100% run found one per record.

   **Fix:** Store and validate stage manifests containing input identities, relevant arguments, and code version; propagate invalidation from Stage 1.

8. **High — [scripts/pipeline.py:58](../../scripts/pipeline.py#L58): `--trust-existing` defeats downstream invalidation.**

   An upstream rerun removes downstream markers, after which the adoption branch trusts the now-stale downstream files.

   **Failure/verification:** A two-stage `_run_stage` probe regenerated `NEW anchors`, invalidated the PAF marker, then adopted `OLD PAF` instead of executing its replacement command.

   **Fix:** Never adopt outputs invalidated during the current invocation. Restrict adoption to an initial, consistent legacy snapshot.

9. **High — [scripts/pipeline.py:204](../../scripts/pipeline.py#L204): a separate source database silently loses every source-only record.**

   Stage 2 reads sources from `src_db`, but Stage 6 receives only the target database and tries to fetch source contigs there.

   **Failure/verification:** Running the supplied source against a target containing only its two empty genomes exited successfully with `records.json=[]`, despite Stage 5 finding both elements. The summary misleadingly printed `Records built: 2`.

   **Fix:** Give Stage 6 separate source and target FASTA handles; fail or explicitly report missing source sequences.

10. **High — [scripts/is_detect.py:243](../../scripts/is_detect.py#L243): different HMMs with the same basename share one search result.**

    The domtbl cache key is only the basename, so an ALL-domains requirement can become repeated checks of the first profile.

    **Failure/verification:** Using real `hmmbuild` and `hmmsearch`, `/first/model.hmm` matched one protein while `/second/model.hmm` independently matched none. Stage 1 nevertheless emitted that protein as passing ALL domains, with `model,model`.

    **Fix:** Key searches by profile content/path and search parameters, and verify each required profile independently.

11. **High — [scripts/pipeline.py:199](../../scripts/pipeline.py#L199): Stage 5 filters do not constrain final records.**

    Stage 6 recomputes calls without consuming Stage 5 output or receiving its size limits. `--min-is-size` is also unused inside Stage 5.

    **Failure/verification:** A fresh pipeline run with `--max-is-size 1000` classified the fixture’s 1,400/1,450-bp removals as `very_short` in Stage 5, yet Stage 6 emitted both elements with two empty observations each.

    **Fix:** Share filtering configuration and accepted observations across stages; enforce both size bounds.

12. **High — [scripts/build_records.py:535](../../scripts/build_records.py#L535): SNP-tolerant calls can crash JSON serialization.**

    Alignment coordinates remain NumPy integers and propagate into records.

    **Failure/verification:** Removing only the exact empty genome from the supplied synthetic database left SNP-divergent empties. The actual pipeline crashed in Stage 6 with `TypeError: Object of type int64 is not JSON serializable`, leaving partial JSON.

    **Fix:** Convert alignment-derived coordinates and metadata to native Python scalar types before constructing records.

13. **Medium — [scripts/lib_alleles.py:309](../../scripts/lib_alleles.py#L309): tolerant decomposition ignores terminal indels in its dominance check.**

    Only gaps between aligned blocks contribute to insertion/deletion accounting.

    **Failure/verification:** For `short=L900+R900` and `long=extra1200+L900+insert1500+R900`, the function accepted the 1,500-bp insertion with identity/coverage `1.0`, despite another 1,200-bp extra block violating dominance `3`.

    **Fix:** Account for leading and trailing gaps when calculating competing indels, coverage, and event class.

14. **Medium — [scripts/lib_alleles.py:248](../../scripts/lib_alleles.py#L248): representative selection can discard the only valid assembly placement.**

    The preferred twin is chosen before checking whether its anchors form a usable pair.

    **Failure/verification:** A GCA assembly had a valid same-contig pair; its preferred GCF twin had upstream/downstream hits on different contigs. `place_per_assembly` returned `{}`.

    **Fix:** Evaluate placements first, then select an eligible representative with an explicit fallback policy.

15. **Medium — [scripts/build_records.py:149](../../scripts/build_records.py#L149): minus-strand fallback coordinates disagree with reported offsets.**

    Transposase-oriented offsets are applied directly to genomic left/right coordinates.

    **Failure/verification:** A minus-strand 1,200-bp ORF with inferred element length 1,401 reported offsets `-100,+101`, but emitted coordinates implied `-101,+100`.

    **Fix:** Convert offsets according to strand before constructing genomic coordinates; recompute offsets after any clamping.

16. **Medium — [scripts/build_records.py:193](../../scripts/build_records.py#L193): `event_key` splits a known shared insertion after a neutral flank SNP.**

    The key hashes the exact element plus 100-bp flanks, rather than identifying an insertion at an orthologous locus.

    **Failure/verification:** Two cloned source sequences differing by one upstream SNP produced identical element sequences and coordinates, but keys `3a6a2026f500bc64` and `212caba8c8e0cea5`. The distinct-event summary therefore counted two.

    **Fix:** Establish locus identity through homologous context and separate it from exact allele identity; do not interpret unique exact-context hashes as insertion-event counts.

17. **Medium — [scripts/pipeline.py:209](../../scripts/pipeline.py#L209): Stage 6 caching does not require its sequence FASTAs.**

    Only JSON and summary TSV are checked, although `is_elements.fa` is a primary output.

    **Failure/verification:** Deleting `records/is_elements.fa` from a completed temporary run and rerunning exited successfully without recreating it. Stage 2 similarly omits `anchors_refs.fa` from its required outputs.

    **Fix:** Include every promised output in completion validation.

18. **Medium — [scripts/pipeline.py:373](../../scripts/pipeline.py#L373): the documented census setting can still omit most assemblies.**

    Setting `-N` to the number of distinct assemblies does not allow for GCA/GCF twins or multiple locus copies, which consume separate alignment slots.

    **Failure/verification:** An actual minimap2 database containing 40 assembly cores and 80 twin genomes yielded **17** paired cores at `-N 40`; `-N 80` yielded all 40.

    **Fix:** Deduplicate targets before searching or increase the cap adaptively until searches are uncapped; correct the census guidance.

19. **Medium — [scripts/extract_anchors.py:194](../../scripts/extract_anchors.py#L194): an empty keep-list retains every hit.**

    `if keep_ids` conflates an explicitly empty restriction with no restriction.

    **Failure/verification:** Running extraction with an empty `--keep-ids` FASTA retained both supplied synthetic references; it should retain zero.

    **Fix:** Test `keep_ids is not None`.

20. **Medium — [tests/check_synthetic.py:57](../../tests/check_synthetic.py#L57): the end-to-end checker accepts coordinate corruption and duplicate records.**

    The coordinate assertion joins endpoint checks with `or`; the dictionary construction at line 41 silently collapses duplicate IDs.

    **Failure/verification:** I increased one record’s `source_end` by **100,000 bp** and appended a duplicate record. The checker still exited zero with `ALL CHECKS PASSED`.

    **Fix:** Validate both coordinates and their relationship to strand, length, offsets, and sequence; assert record cardinality and ID uniqueness before constructing dictionaries. Add SNP-only serialization, clipping, cache-change, and competing-consensus cases.

21. **Low — [scripts/is_detect.py:279](../../scripts/is_detect.py#L279): `domains_hit` lists profiles that did not match.**

    Every output row receives all requested profile names.

    **Failure/verification:** A real `--any-domain` run logged `first: 1 hits`, `second: 0 hits`, but emitted `domains_hit=first,second`.

    **Fix:** Construct the field from each protein’s actual profile memberships.

**(a) Checked and found correct**

- Full-length anchor extraction and decomposition coordinate conversion on both strands, including clamped flanks.
- Actual Prodigal and Pyrodigal protein/GFF mapping: 50/50 proteins mapped for each caller.
- Exact decomposition, rightmost placement, `_slide` extent, and `canonical_insert_key` reverse-complement invariance across 1,000 randomized cases.
- The supplied microhomology vote-normalization case, basic TSD detection, and ordinary marker failure/invalidation tests.

The complete supplied suite passed; see test.log (temporary file, not kept).

**(b) Could not verify**

No production database or independently curated boundary truth was supplied, so I could not quantify affected published records or biological error rates. Large-scale SLURM interruption, concurrent writers, and HMM sensitivity/specificity were not assessed.