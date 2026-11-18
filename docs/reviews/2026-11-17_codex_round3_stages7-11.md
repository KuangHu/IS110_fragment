# Adversarial review of stages 7-11 (commit 6ab9c9e)

Reviewer: OpenAI codex CLI, 2026-11-17. Fixes: the commit that adds this file.
Not fixed (see README "Known limits"): #17 cross-contig duplication candidates, #19/#20 viz_lineage colouring and added V_REF track, #28 event-key context not junction-invariant, #31 short tandem units.

Reviewed commit **`6ab9c9e`**, which was HEAD when the review started. HEAD subsequently advanced to an unrelated movement-script commit; downstream implementation files remained unchanged. The test runner gained one additional test invocation.

**The existing test suite passes, but the counterexamples below reproduce wrong classifications, false CDS validation, missed observations, and invalid statistical assumptions.** I made no repository edits.

Reproduction scripts, inputs, outputs, and logs are preserved under `<tmp>/`. Each finding identifies its reproduction script. Unless explicitly described as a supplied-PAF or isolated-transform test, reproductions used the installed tools.

1. **High — `scripts/find_variants.py:391`, `scripts/validate_v0_cds.py:148`: alignment envelopes count deleted CDS bases as present.**  
   A 10,000 bp reference with a 900 bp CDS loses 300 bp inside that CDS. Minimap2 produces one 97%-identity alignment spanning reference `[0,10000)`. Actual output: `total_matched_bp=10000`, followed by CDS coverage `100.0`, `passes=True`. Expected CDS coverage: **66.7%, FAIL**.  
   Verified through stages 8–10a in `probe_stages.py`; raw alignment is in `core.log`. Unioning envelopes fixes overlap double counting but still includes CIGAR deletions. **Fix:** derive covered reference intervals from aligned CIGAR operations, excluding reference deletions.

2. **High — `scripts/find_variants.py:420`: inversions and rearranged reference content become “clonal.”**  
   With unchanged flanking anchors, both `revcomp(V1)` and a rotation `V1[5000:]+V1[:5000]` are classified `clonal_V1`. Expected: neither should contribute a clonal observation; use `partial` or an explicit rearranged category.  
   `probe_stages.py` produced `('inverted','clonal_V1',10000,1)` and `('rotated','clonal_V1',10000,2)`. **Fix:** require compatible strand, order, and coordinate correspondence for clonal calls.

3. **High — `scripts/build_lineages.py:138`: query coverage does not establish sequence containment.**  
   For independent random 1,500 bp sequences `X,Y`, the metrics report 100% coverage and identity for both:
   - query `X+Y` against target `Y+X+random300`;
   - query `X+X` against target `X+random1800`.

   The latter reuses one target copy to explain two query copies. Neither establishes nested retention of the shorter sequence. Verified with real minimap2 in `probe_core.py`, outputs `PAIR_SCRAMBLED` and `PAIR_DUPLICATE`. **Fix:** validate a strand-consistent, ordered, non-reusing alignment path.

4. **High — `scripts/detect_tandem.py:108`: disconnected hits become one continuous array.**  
   Supplied valid forward hits `[0,300)→[1000,1300)` and `[3000,3300)→[4000,4300)` describe two dispersed repeat pairs. Actual: a **four-copy, 1,000 bp tandem** spanning `[0,4300)`. Two separate two-copy 300 bp arrays similarly become one **12-copy array**. Expected: separate findings without counting intervening sequence.  
   Verified in `probe_more.py`. **Fix:** partition offset groups into connected arrays and calculate support from actual aligned intervals.

5. **High — `scripts/detect_tandem.py:94`: valid repeats are discarded when minimap2 reports the reciprocal coordinate direction.**  
   A real sequence containing two independent, exact, two-copy 300 bp arrays returns no findings. Its valid repeat hits have `ts < qs`; the only other alignment fails identity. Expected: two arrays.  
   Verified in `probe_core.py`, output `TANDEM_SEPARATE`, including the actual PAF-derived hits. **Fix:** normalize reciprocal forward hits before grouping; do not require minimap2 to choose `ts > qs`.

6. **High — `scripts/pipeline.py:307`: downstream cache checks have no input dependencies.**  
   Run `detect_growing_and_tandem()` with empty records, then replace that same file with a valid record. Actual: all downstream stages remain cached, stage-8 marker is unchanged, observations remain `[]`. Expected: recomputation.  
   Verified using the real pipeline API in `probe_stages.py`, output `PIPELINE_CACHE`. **Fix:** supply stage-specific `inputs`, including records, database/index, observations, and lineages.

7. **High — `scripts/find_variants.py:247`: a same-size database replacement reuses stale anchor placements.**  
   After a valid run, replace both genomic anchors with unrelated sequences while preserving FASTA path and byte size. Actual cached run still reports one observation; a fresh output directory reports zero.  
   Verified in `probe_stages.py`: `CACHE_SAME_SIZE 1 cached observations vs 0 fresh`. **Fix:** fingerprint database content or use a reliable immutable database identifier; size alone is insufficient.

8. **High — `scripts/build_lineages.py:94`: `.clstr` parsing assumes the representative precedes its members.**  
   Real cd-hit output can list the shorter member before the starred representative. Input order `[small,large]` returns only `{'large':'large'}`; reversing the same sequences returns mappings for both. Main then treats an omitted member as its own cluster.  
   Verified in `probe_stages.py`, including both actual `.clstr` files. **Fix:** collect the entire cluster, find its starred representative, then map every member.

9. **High — `scripts/build_lineages.py:84`: the proportional length cutoff still erases qualifying absolute-length steps.**  
   An exact 30,000 bp sequence and that sequence plus 250 bp cluster together at default identity. The 250 bp difference exceeds both the standalone 200 bp and pipeline 10 bp minimum step. Expected: retain the length step.  
   Verified with real cd-hit in `probe_stages.py`, representative-first order. **Fix:** enforce the absolute step-size constraint during clustering, rather than relying on `-s 0.99`.

10. **High — `scripts/pipeline.py:314`, `scripts/find_variants.py:48`: stage 8 cannot honor the pipeline’s census setting.**  
    Sixty genomes with identical valid clonal loci produce **9 observations** using stage 8’s default `-N 20`, versus **60** with `--max-secondary 100`. The pipeline’s `--max-secondary` is not forwarded to this stage. Expected: a requested census must reach stage 8 and disclose capped searches.  
    Verified in `probe_more.py`, outputs `CENSUS`. **Fix:** forward the setting and report saturation.

11. **High — `scripts/find_variants.py:422`: counting PAF records rejects genuine insertions.**  
    A 10 kb reference gains a 300 bp insertion and loses one unrelated base. Minimap2 returns one 97.08%-identity block covering the reference; quadratic decomposition is skipped and the gap-free fallback refuses the extra indel. Actual: `partial`. Expected: insertion with 299 bp net growth and essentially complete reference retention.  
    Verified through stage 8 in `probe_stages.py`. **Fix:** inspect insertion evidence within CIGAR/alignment structure instead of requiring two PAF records.

12. **High — `scripts/find_rearrangements.py:103`: selecting a preferred twin before checking usable evidence loses candidates.**  
    Supply a complete opposite-strand anchor pair in GCA and only the upstream anchor in its preferred GCF twin. Actual: GCA is deleted, GCF is subsequently skipped, and all catalogues are empty. Expected: retain the usable GCA candidate.  
    Verified with supplied valid PAF in `probe_rearr.py`, case `twin_incomplete`. **Fix:** evaluate usability before twin selection, with fallback as in `place_per_assembly`.

13. **High — `scripts/verify_tandem_enrichment.py:114`: exact-string deduplication does not remove the stated pseudoreplication problem.**  
    An element, its one-SNP relative, and its reverse complement all remain separate tests; only the identical duplicate disappears. Genomic controls can also repeat the same physical window. Expected: independence must be established or dependence incorporated into inference.  
    Verified in `probe_stats.py`: retained IDs `['a','near','rc']`; three controls contained one unique sequence. `probe_stats_extra.py` demonstrates Fisher p changing from `0.5` to `1.10e-59` under clone replication. **Fix:** canonicalize orientation, cluster near-identical elements, track control coordinates, and use appropriate grouped inference.

14. **Medium — `scripts/build_lineages.py:356`: length buckets require similarity only to their first member.**  
    Three equal-length sequences have A–B and A–C identity approximately 96%, while B–C has no qualifying relationship. Actual: all six observations merge into representative B. Expected: sequences failing the bucket’s mutual-similarity requirement must remain separate.  
    Verified through real clustering/alignment and lineage construction in `probe_more.py`, `BUCKET_RESULT`. **Fix:** validate every member against the eventual representative or use complete-linkage membership.

15. **Medium — `scripts/find_rearrangements.py:138`: projected-start collapsing still creates false duplications across indels.**  
    A full alignment `50=1D450=` and its clipped `450=` suffix describe one placement. Their approximate projected starts differ by one base. Actual: duplication with multiplicity two and estimated spacing 51 bp. Expected: one placement.  
    Verified with supplied CIGAR-bearing PAF in `probe_rearr.py`, `indel_clipped`. **Fix:** collapse using agreement of overlapping query-to-target mappings, not an ungapped projection alone.

16. **Medium — `scripts/find_rearrangements.py:154`: “normal” ignores order and distance, and can suppress an inversion.**  
    Upstream at 10 kb and downstream at 1 kb, both `+`, produces `normal: 1`. Adding a distant same-strand downstream hit to a valid local opposite-strand pair removes that inversion from the all-inversion catalogue. Expected: validate a compatible normal locus before using it to suppress other evidence.  
    Verified in `probe_rearr.py`, `normal_wrong_order` and `inversion_suppressed`. **Fix:** classify compatible placements/loci rather than arbitrary cross-products.

17. **Medium — `scripts/find_rearrangements.py:227`: copies on separate contigs cannot produce a duplication candidate.**  
    Two contigs in one assembly each carry a complete anchor pair at distinct positions. Actual: `normal: 1`, no duplication. Expected: retain a cross-contig duplication candidate or explicit ambiguity; contig separation does not establish assembly redundancy.  
    Verified in `probe_rearr.py`, `dup_two_contigs`. **Fix:** consider assembly-wide placements and assess redundancy separately.

18. **Medium — `scripts/find_variants.py:455`, `scripts/detect_tandem.py:215`, `scripts/viz_lineage.py:178`: CDS mapping loses strand and assumes gap-free offsets.**  
    For reverse alignment reference `[0,1000)` to variant `[100,1100)`, CDS `[100,300)` should map to `[800,1000)`. Actual mapping: `[200,400)`. Stage 8 discards the alignment strand needed to fix this.  
    Verified directly and through tandem metadata construction in `probe_final.py`; the latter isolates the transform with supplied repeat hits. **Fix:** retain strand/CIGAR and map aligned coordinates explicitly.

19. **Medium — `scripts/viz_lineage.py:137`: inherited coloring moves cargo to the wrong position across gaps.**  
    Insert 100 bp at position 1,500 in a 3 kb sequence. Actual coloring marks `[0,3000)` inherited and `[3000,3100)` new. Expected new cargo: `[1500,1600)`.  
    Verified with real minimap2 in `probe_final.py`, `GAPPED_COLOR`. **Fix:** propagate colors through CIGAR-aligned blocks.

20. **Medium — `scripts/viz_lineage.py:216`: the figure adds a reference step absent from the selected chain.**  
    Supply a valid two-step chain of nested 1,500/2,000 bp variants unrelated to the 3,000 bp reference. Actual figure has three tracks, adding V_REF as the largest chain step. Expected: preserve the selected chain; any comparison reference must be shown separately.  
    Verified by actual rendering and captured track labels in `probe_final.py`, `VIZ_ADDED_VREF`.

21. **Medium — `scripts/verify_tandem_enrichment.py:93`: “called boundaries only” stops filtering when no record was successfully called.**  
    Input records explicitly labelled transposase-only and symmetric fallback, with no decomposition call. Actual: both are kept under default options. Expected: zero called elements.  
    Verified in `probe_stats.py`, `CALLED_ONLY_NO_CALLS`. **Fix:** distinguish an unlabelled legacy schema from a labelled dataset containing zero successful calls.

22. **Medium — `scripts/verify_tandem_enrichment.py:102`: controls can be the IS sequences being tested.**  
    In supported schema A, loci use `start/end`; this loader only reads `source_start/source_end`, leaving exclusion empty. Separately, a valid schema-B GCA locus leaves its identical GCF copy eligible. Both fixtures produced three controls identical to the test IS.  
    Verified in `probe_stats.py` and `probe_edges.py`. Expected: exclude those loci/copies. **Fix:** normalize schemas and screen known/homologous IS copies, including twins.

23. **Medium — `scripts/verify_tandem_enrichment.py:134`: control sampling biases eligible starts and can lose length matching.**  
    For a 1 kb and 2 kb contig, sampling 1 kb windows assigns about one-third of controls to the short contig, although it has only one of 1,002 eligible starts. Actual: **331/1000**, versus approximately **1/1000** for uniform eligible starts. A fragmented database also yielded three short controls but only one long control.  
    Verified in `probe_stats_extra.py`. **Fix:** weight eligible starts for each requested length and explicitly handle insufficient matched controls.

24. **Medium — `scripts/verify_tandem_enrichment.py:274`: ULTRA’s default period ceiling censors the large-repeat comparison.**  
    A four-copy 200 bp array is detected by TRF and missed by the current ULTRA invocation. Explicit ULTRA `-p 500` detects it. Expected: align or explicitly restrict the detectable period ranges before interpreting caller agreement.  
    Verified with both real callers in `probe_stats.py`, `MAX_PERIOD_DEFAULT` and `MAX_PERIOD_500`. Installed ULTRA defaults to maximum period 100; TRF is configured for 500.

25. **Medium — `scripts/verify_tandem_enrichment.py:243`: shuffle retry exhaustion silently substitutes the original sequence.**  
    For `A200+C200+A200+C200`, default shuffling returns the original in **86/100** trials. With 10,000 retries, **0/100** remain unchanged and all outputs are distinct. Expected: a valid randomized null or an explicit shuffle failure.  
    Verified in `probe_core.py` and `probe_stats_extra.py`. **Fix:** use a robust spanning-tree sampler or explicit failure handling instead of silently returning input.

26. **Medium — `scripts/verify_tandem_enrichment.py:325`: unrestricted permutation ignores the matched shuffle design.**  
    For nine matched pairs `test=control+0.01`, with heterogeneous baseline densities, current permutation gives `p=0.4685`; exact within-pair label exchange gives `p=0.001953125`. Expected: preserve the pairing/dependence of each IS and its composition-matched shuffle.  
    Verified in `probe_stats_extra.py`, `PAIRED_NULL`. **Fix:** use paired randomization and corresponding matched analyses for shuffle controls.

27. **Medium — `scripts/detect_tandem.py:237`: failed alignments still produce a successful negative catalogue.**  
    Force minimap2 to exit 1 for a valid element. Actual: process exit 0 and `tandem_findings.json=[]`, with only a stderr warning. The pipeline can mark this result complete. Expected: failed/unknown observations must remain distinguishable from successful negatives.  
    Verified using a temporary failing executable in `probe_final.py`, `FAILED_ALIGN_SUCCESS`. **Fix:** fail the stage or emit machine-readable failures that prevent complete-result caching.

28. **Medium — `scripts/find_variants.py:472`: the complete nested event key is not junction/orientation invariant.**  
    An insertion with two-base junction ambiguity and its reverse-complement representation get identical insert keys but different context hashes: `…:2d4f63e053a8366e` versus `…:31b1ecea5945186d`. Expected: one event key.  
    Verified using the exact stage-8 key construction in `probe_edges.py`, `EVENT_CONTEXT_KEYS`. **Fix:** canonicalize the context over equivalent junction placements as well as the inserted sequence.

29. **Medium — `scripts/build_lineages.py:173`: the branching fix introduces exponential materialization.**  
    A layered containment graph with two alternatives per layer produces 256, 1,024, and 4,096 chains for 16, 20, and 24 variants. Forty variants raise `MemoryError` under a 128 MiB limit in about 1.4 seconds. Expected: bounded or explicitly controlled processing of plausible multi-variant sites.  
    Verified in `probe_more.py` and `probe_chain_memory.py`. **Fix:** retain a DAG and enumerate chains lazily with explicit limits/status.

30. **Medium — `scripts/find_rearrangements.py:177`: hash seed changes the selected inversion locus.**  
    Two equal-span candidate pairs on different contigs select `c2` under `PYTHONHASHSEED=0` and `c1` under seed 1. These are the coordinates passed to downstream validation. Expected: deterministic selection or retention of both.  
    Verified in `probe_rearr.py`, `HASH_SEEDS`. **Fix:** sort candidates and apply explicit coordinate tie-breakers.

31. **Medium — `scripts/detect_tandem.py:48`: the advertised minimum repeat size is not supported by the alignment settings.**  
    A 2,180 bp sequence containing three exact copies of a random 60 bp unit yields no self-alignment hits and no finding. It exceeds the default sequence-length cutoff and its unit exceeds `--min-repeat-len 50`. Expected: detection.  
    Verified with real minimap2 in `probe_more.py`, `SHORT_UNIT`. **Fix:** use settings/callers sensitive to the requested unit and array lengths.

32. **Low — `scripts/verify_tandem_enrichment.py:314`: zero versus zero becomes infinite enrichment.**  
    Actual complete report: `mean_IS=0`, `mean_control=0`, `fold_enrichment=Infinity`. Expected: undefined/no-estimable enrichment, not an infinite effect. The emitted `Infinity` is also nonstandard JSON.  
    Verified end to end in `probe_stats.py`, `ENRICHMENT_ZERO`. **Fix:** explicitly handle 0/0 and serialize undefined/infinite values consistently.

33. **Low — `scripts/viz_tandem.py:63`: dispersed findings are drawn and labelled as contiguous tandems.**  
    A dispersed finding with two 500 bp copies spaced 1,500 bp apart is drawn as two adjacent 1,000 bp blocks spanning the gap, labelled “Tandem.” Expected: separate copies with an intervening gap.  
    Verified by actual rendering and captured rectangle geometry in `probe_final.py`, `VIZ_DISPERSED_RECTS`. **Fix:** branch on `repeat_class` and use spacing/copy intervals.

34. **Low — `scripts/find_variants.py:417`: the deletion threshold excludes equality.**  
    With default minimum shortening 200 bp, a 200 bp deletion is `clonal_V1`; a 201 bp deletion is `deletion`. Expected: the stated minimum includes 200 bp.  
    Verified through stage 8 in `probe_edges.py`, `DELETION_THRESHOLD`. **Fix:** use an inclusive threshold.

For the individual claims in the starting HEAD commit message:

| Commit-message item | Verdict |
|---|---|
| Align each between-sequence only to its own V1 | **Correct** targeting fix; classification limitations remain |
| Preserve observations from two V1s sharing an interval | **Correct**; reproduced six observations for each reference |
| Union V1 coverage on reference coordinates | **Incomplete**: gapped envelopes still count absent bases |
| Raise on failed samtools batches | **Correct**; existing failure test passes |
| Key anchor cache on anchors, DB, and `-N` | **Incomplete**: DB identity is effectively path plus size |
| Count skipped records | **Incomplete**: records lacking either reference-ID field still disappear without a skip counter |
| cd-hit length cutoff preserves nested variants | **Incomplete**: absolute qualifying steps can still collapse |
| Raise on cd-hit failure | **Correct** |
| Use representative’s own record | **Correct** when the cluster mapping is complete |
| Union pairwise coverage | **Incomplete**: gapped, reordered, and reused-target matches remain unsafe |
| Scale pairwise `-N` with variant count | **Correct improvement**; not proof of exhaustive placement reporting |
| Branch chains from deepest contained step | **Correct for the tested simple branches; new error:** unbounded exponential expansion |
| Merge length buckets only when mutually similar | **Incomplete**: only similarity to the initial member is checked |
| Union CDS coverage | **Incomplete**: overlap arithmetic is fixed, missing CIGAR bases are not |
| Carry direction labels into validated output | **Correct** |
| Require adjacent tandem copies; default-off dispersed findings | **Incomplete**: disconnected groups and tolerance-based false calls remain |
| Correct `n_copies` | **Correct for clean single arrays; incomplete** for disconnected groups |
| Report true offset rather than rounded units | **Correct for the tested clean arrays** |
| Count and report self-alignment failures | **Correct as a stderr diagnostic; incomplete** as successful-result handling |
| Collapse rearrangement hits per physical placement | **Incomplete** across indels |
| Parse assembly as in stages 4–6 | **Correct** consistency change |
| Pass `.fai` for translocation depth | **Correct**; reproduced depth 2,000 bp with the expected contig lengths |
| One enrichment test per distinct element | **Incomplete**: exact strings are not independent elements |
| Called boundaries only by default | **Incomplete** when no successful boundary call exists |
| Add eight downstream tests | **Correct**, but substantial failure modes remain untested |

I checked and found the following working: the full supplied suite and synthetic integration run; clean tandem copy counts and offsets; exclusion of tested inverted repeats; tested mildly diverged tandem copies; union arithmetic for overlapping intervals; duplicate samtools region extraction; the shared-region observation fix; existing strand-oriented record handling; and preservation of `direction="unpolarized"`. An isolated duplicated-half-CDS test correctly reports **50%, FAIL**, confirming that the new union arithmetic itself works.

The important test gap is visible in `tests/check_synthetic.py:100`: downstream checking ends with the lineage direction label. It does not validate CDS outcomes, tandem findings, rearrangement catalogues, figures, or enrichment statistics. The pairwise test in `tests/test_downstream.py:103` checks that coverage stays below 100%, but not whether it establishes containment.

I could not verify the quoted E. coli before/after counts, quantify effects on the paper’s tables, or assess contamination rates in its production controls: those production inputs/results were not supplied. I also did not establish exhaustive chain correctness, shuffle uniformity across all compositions, or behavior across other external-tool versions. The reproduced shuffle composition checks did preserve first/last bases and dinucleotide counts on 100 random sequences; that does not cure the retry fallback bias.