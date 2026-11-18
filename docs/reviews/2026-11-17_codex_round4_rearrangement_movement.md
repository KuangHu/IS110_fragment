# Adversarial review of the rearrangement / movement scripts (commit a7842d2)

Reviewer: OpenAI codex CLI, 2026-11-17. Fixes: commits 7835a00, 3b1b3fb and the
commit that adds this file.

NOT fixed -- these verdicts remain weaker than their names (see README "Known limits"):
- #3  inversion checks accept coverage assembled from inverted PARALOGS: two identical
      genomes that both carry M and rc(M) can pass. A coherent single alignment path
      through the changed junctions is not required.
- #6/#7 movement: flank insertions or a deleted element end no longer give MOVED (they
      are DIFFERENT_CONTEXT), but they are not recognised as "same site" either.
- #15/#16 rearrangement_validator.py (MUMmer show-diff): misses a clean tandem
      duplication; is_mediated mixes reference and query coordinates. Use
      validate_rearrangements_gold.py; treat is_mediated as unreliable.
- #21 inversions shorter than ~6 kb (full check) are reported as not assessable.
- #23 duplication_tandem_self_detect.py misses arrays whose unit contains an IS.
- #33 one-base IS start offset in the three review/figure scripts
      (prepare_review_cases.py, inversion_focal_pair.py, composite_transposon_review_pair.py),
      which also still use the old contig-name parser (#24, two of them).
- #31 (part) event identity in the collector is (category, ref_id, target): the same
      target event reached from several reference elements is several records.
- #1 is handled by requiring each contig to run past the element; a translocation whose
      contig ends just beyond the element is therefore "not assessable", not rejected.

Reviewed HEAD `a7842d2`. **The revised scripts still produce unsupported confirmed events, miss clean events, and corrupt some deposited metadata.** Both supplied test suites pass.

The repository is unchanged. All fixtures, commands, outputs, and reproducers are in the evidence archive (temporary file, not kept). The runner (temporary file, not kept) reproduces the review from this checkout.

Below, `E` denotes a synthetic element, `D` a synthetic segment, and other letters independent random backbone sequences. Lengths and complete FASTAs are preserved with each case.

1. **Critical — `scripts/validate_rearrangements_gold.py:437`; `scripts/is_translocation_check.py:177`: an assembly break becomes a confirmed translocation.**  
   Input: reference `A–E–B`, with a 16-kb E; target contains exactly the same sequence split into `A–E[:8000]` and `E[8000:]–B`. Actual: bilateral `passes=True`; gold `CONFIRMED_BOTH`, `min_anchor_end_dist=9000`. Expected: no biological translocation; assembly-break uncertainty. Anchor depth alone does not establish junction continuity through a long element. **Fix:** assess continuity and complete element/junction evidence, accounting for element length. Verified with real minimap2 and dnadiff: case output (temporary file, not kept).

2. **Critical — `scripts/validate_rearrangements_gold.py:237,415`; `scripts/is_duplication_check.py:336`; `scripts/duplication_strict_canonical_check.py:284`: duplication confirmation does not require copy gain relative to the reference.**  
   An unchanged genome `A–E–D–X–E–D–B` passes both duplication checks (`dup_bp=12000`). Separately, reference and target each contain two D copies, but target has an unrelated 6-kb insertion between them: gold reports `CONFIRMED_BOTH`, `dup_bp_outside_IS:27000`. Expected: no new duplication. **Fix:** compare orthologous reference/target copy structure and establish a gain, rather than detecting repeated sequence already present in both. Unchanged-genome run (temporary file, not kept), gold run (temporary file, not kept).

3. **High — `scripts/inversion_revcomp_check.py:157`; `scripts/inversion_full_check.py:310`: inverted paralogs can confirm an inversion between identical genomes.**  
   Identical reference/target genomes containing M and reverse-complement M pass the flank recheck. An unchanged `A–E–M–rc(E)–rc(M)–E–B` also passes the full check: `M_cov=0.927`, both transition checks true, `all_pass=True`. Expected: reject an inversion between identical genomes. Coverage can be assembled from different reference copies. **Fix:** require one coherent orthologous alignment path supporting the changed junctions. Flank run (temporary file, not kept), full run (temporary file, not kept).

4. **High — `scripts/lib_movement.py:246,254`; `scripts/transposon_mechanism_classify.py:207`: similar gap length is mistaken for the same element.**  
   Reference has E at X and unrelated, same-length F at Y; target moves E to replace F. Actual: `SYNTENIC_OTHER_REF_COPY`, `ref_site=occupied`. A corresponding origin replacement is classified `COPY_PASTE` merely because its gap remains 15,000 bp. Expected: movement/replacement, with the original element absent from its origin. **Fix:** verify the intervening body’s sequence identity and structure before declaring another reference copy or retained origin. Movement run (temporary file, not kept), mechanism run (temporary file, not kept).

5. **High — `scripts/lib_movement.py:244`: short occupied sites are classified empty.**  
   The empty-site interval accepts gaps through 2,000 bp before checking occupancy. Two unchanged 1,500-bp IS copies with shared 50-bp padding therefore produce two **strict** `MOVED` rows with `ref_site=empty`. Expected: both cross-copy comparisons recognized as existing reference copies. **Fix:** make empty/occupied hypotheses distinguishable for short elements and verify the body. Run (temporary file, not kept).

6. **High — `scripts/lib_movement.py:163,177`: unrelated flank insertions become element movement.**  
   E remains between the same backbone sequences; insert 6 kb of unrelated DNA 2.2 kb away on each side. Actual: both flank coverages `0.44`, `MOVED`, `strict_IS=True`, `ref_site=absent`. Expected: same locus with flank insertions, or unresolved context—not confirmed element relocation. **Fix:** establish origin/destination correspondence using longer backbone alignments and junctions; flank mismatch alone is insufficient. Run (temporary file, not kept).

7. **High — `scripts/lib_movement.py:116`: projection assumes missing query ends remain physically present.**  
   Reference `A–E12kb–B`; target `A–E[3000:9000]–B`. Actual: projected flanks extend into backbone, both match fractions become `0.40`, verdict `MOVED`. Expected: an element deletion at the same site. **Fix:** distinguish divergence from deleted overhangs; require boundary evidence or mark projected boundaries uncertain. Run (temporary file, not kept).

8. **High — `scripts/lib_movement.py:80`: proximity clustering conflates copies and accepts scrambled bodies as intact.**  
   Target `A–E–E–B` produces one cluster spanning both copies, `n_blocks=2`, `SYNTENIC`. Target containing swapped E halves also gets `IS_cov=1.000`, `IS_ident=100.0`, `strict_IS=True`. Expected: two copy placements in the first case; a rearranged body in the second. **Fix:** chain blocks by query order, orientation and consistent placement; split query-coordinate resets. Adjacent-copy run (temporary file, not kept), scrambled-body run (temporary file, not kept).

9. **High — `scripts/is_transposon_movement_check.py:347`: the pairwise checker omits reference-paralog refinement.**  
   Identical genomes carrying the same two composite units yield two `MOVED` rows, one in each cross-copy comparison. The DB checker correctly changes these to `SYNTENIC_OTHER_REF_COPY`. Expected: zero movements in either implementation. **Fix:** apply reference-site comparison consistently to the pairwise checker. Run (temporary file, not kept).

10. **High — `scripts/duplication_strict_canonical_check.py:280,284`: “canonical” does not establish `IS–D–IS–D–IS`.**  
    Target `A–D–E–D–B`, containing only the middle IS, passes `canonical_pass=True`. Two widely dispersed D copies with one IS somewhere between them also pass. Expected: fail canonical topology. **Fix:** require three distinct boundary IS copies, appropriate order/adjacency, and coherent outer-flank placement around the duplicated unit. Missing-outer-IS run (temporary file, not kept), dispersed run (temporary file, not kept).

11. **High — `scripts/validate_rearrangements_gold.py:333`: actual Stage 7 duplication rows cannot reach validation.**  
    Stage 7’s duplication catalogue intentionally writes only the duplicated anchor side, e.g. `up=[] down=[…;…]`. A clean duplication in this format returns `NO_ANCHOR_PAIR`; both standalone duplication checks pass it. Expected: duplication evidence assessed. **Fix:** make anchor requirements category-specific and consume the actual producer schema. Run (temporary file, not kept).

12. **High — `scripts/lib_rearr.py:159,173`: `copy_depth` returns incorrect depth and depends on block order.**  
    Real minimap2 alignments with query intervals `[0,1500]`, `[0,1500]`, `[749,2250]` at separate target places return `(dup_bp=1500, max_copies=0)`. Expected: at least two copies over 1,500 bp. Separately, reordering equal-length overlapping target blocks changes `(0,1)` to `(1500,2)`. **Fix:** accumulate contiguous coverage at depth ≥k across intervening endpoints, and replace order-dependent greedy overlap removal. Depth run (temporary file, not kept), order run (temporary file, not kept).

13. **High — `scripts/is_duplication_check.py:230`: verdicts depend on `PYTHONHASHSEED`.**  
    A clean duplication whose 500-bp anchor is palindromic has valid `+` and `-` hits at each position. Actual: seed 0 gives `True`, seed 2 `False`, seed 6 `NA`. The set followed by sorting only on start arbitrarily discards strand alternatives. Expected: deterministic assessment retaining the supported placement. **Fix:** preserve alternatives and evaluate strand groups before thinning positions. Twelve-seed run (temporary file, not kept).

14. **High — `scripts/dedup_inversions_per_locus.py:37`; `scripts/dedup_translocations_per_locus.py:42`: anchor bins are not event identities.**  
    A real 40-kb inversion represented from IS_a and IS_b produces **two** `_events.tsv` rows. Conversely, distinct translocation pairs at `(10000,20000)` and `(10900,20900)` collapse into one. Expected: one inversion; preservation of distinguishable breakpoint pairs. **Fix:** deduplicate resolved junction/event representations with explicit uncertainty, rather than rounded raw anchor starts. Real inversion run (temporary file, not kept), collision run (temporary file, not kept).

15. **High — `scripts/rearrangement_validator.py:291,316`: a clean canonical tandem duplication is rejected.**  
    Reference `A–E–D–E–B`; target `A–E–D–E–D–E–B`. Actual: `REJECTED`, `dup_diffs=0`, `multi_map_regions=0`; show-diff describes the change as GAPs, and overlapping reference intervals start in different bins. Expected: detect the duplication. **Fix:** compare overlapping reference intervals at independent target placements; do not require their starts to share a 1-kb bin or the feature label to be DUP. Run (temporary file, not kept).

16. **High — `scripts/rearrangement_validator.py:352,203`: reference breakpoints are tested in target coordinates.**  
    A genuine 80-kb inversion has target endpoints near 100,000 and 180,000; its only target IS is at 40,001–41,500. Actual: `CONFIRMED_IS_MEDIATED`, because a **reference** INV breakpoint near 40,000 is matched against that target IS. Every query INV endpoint is far from it. Expected: inversion confirmed without target IS-adjacency support. **Fix:** retain the reference/query coordinate frame and corresponding assembly for every diff. Run (temporary file, not kept).

17. **High — `scripts/is_translocation_reciprocal.py:319`: independently IS-adjacent C and X are treated as reciprocal evidence.**  
    C aligns beside an IS on reference c3; X beside a different IS on c2. No reference `X–IS–C` partner exists. Actual: `both_is_adjacent=True`, the documented strict-pass condition. Expected: reciprocal linkage unproven. **Fix:** require C and X to flank the same partner IS with compatible orientation and geometry. Run (temporary file, not kept).

18. **High — `scripts/cluster_large_is_elements.py:46,53`: explicitly uncalled boundaries pass when no called record is present.**  
    Input containing only a 7,500-bp `peak_symmetric_fallback` record produces a representative in `reps_records.tsv`. Expected: exclude this explicitly labelled boundary guess. Its eligibility currently depends on other records in the file. **Fix:** inspect each record’s method; distinguish missing legacy labels from explicitly uncalled labels. Run (temporary file, not kept).

19. **High — `scripts/parse_cross_organism_hgt.py:106`: remote fragments are merged into one strict placement.**  
    Two 45%-coverage fragments land approximately 890 kb apart on one contig. Actual: one placement spanning 10,000–904,500, `cov=0.900`, `identity=100.0`, `strict=True`. Expected: two partial placements, neither meeting the coverage threshold. **Fix:** cluster by physical placement and compatible query/target order before calculating union coverage. Run (temporary file, not kept).

20. **High — `scripts/hmm_filter_is110.py:61,178`: duplicate FASTA IDs can substitute a nonmatching sequence into filtered output.**  
    Two sequences share an ID: the first matches both configured artificial profiles; the second matches neither. Actual: the second sequence is emitted as passing, with ORF coordinates overwritten by the second record. Tested alone, that sequence produces no passing row. Expected: reject duplicate IDs or retain record-specific identity throughout. **Fix:** use unique internal IDs and reject ambiguous output identifiers. Run (temporary file, not kept).

21. **Medium — `scripts/inversion_full_check.py:289,298`; `scripts/inversion_revcomp_check.py:138`: padding erases clean small inversions.**  
    A clean 5-kb inversion gets `M_ref_len=0`, `M_pass=False` in the full checker. A clean 2.5-kb inversion is tested against 30-kb proxies by the flank checker and fails. Expected: assessment of the actual inverted segment, or explicit insufficient-resolution status. **Fix:** derive boundaries from alignments/IS positions and cap windows consistently for small spans. Run (temporary file, not kept).

22. **Medium — `scripts/duplication_tandem_db_scan.py:123,171`: twin selection loses IS annotations.**  
    The DB contains identical GCA/GCF twins, while `is_hits.tsv` contains the GCA detections. The scanner selects GCF and reports a clean array `REPEAT_NO_IS`, with all three IS flags false. Expected: `VERIFIED_CANONICAL`. **Fix:** reconcile annotations with the selected assembly, including contig correspondence. Run (temporary file, not kept).

23. **Medium — `scripts/duplication_tandem_self_detect.py:142`: an internal IS hides a genuine tandem array.**  
    Two identical 13.5-kb DNA units, each containing an internal IS, are bracketed by three outer IS copies. Actual: no output rows, because every consecutive inter-IS gap is below 10 kb. Expected: detection of the complete repeated units. **Fix:** consider compatible nonconsecutive boundaries and deduplicate resulting arrays. Run (temporary file, not kept).

24. **Medium — `scripts/is_translocation_reciprocal.py:152`; `scripts/prepare_review_cases.py:76`; `scripts/inversion_focal_pair.py:38`: old contig parsers remain.**  
    `contig-1:10-510(+)` parses as contig `1`; `contig:part` becomes `part`. Renaming reciprocal-test target contigs to hyphen/colon names changes an assessed case into zero output rows. Expected: preserved names and an assessed/status row. **Fix:** use the structural shared parser everywhere. Reciprocal run (temporary file, not kept), review-parser run (temporary file, not kept).

25. **Medium — `scripts/is_large_element_movement_db.py:88`; `scripts/is_solo_long_movement_db.py:100`: prefixed records do not work with ordinary assembly FASTAs.**  
    Input contig `GCA_…|c1`, as produced by records/is_hits, and reference FASTA contig `c1`. Actual: lookup retries `GCA_…|GCA_…|c1`, skips all candidates and produces no table. Expected: resolve `c1`. **Fix:** explicitly support qualified and local names without blindly adding a prefix. Run (temporary file, not kept).

26. **Medium — `scripts/transposon_mechanism_classify.py:220`: the new failure-row path crashes.**  
    A MOVED row whose reference IS pair cannot be found sets `fc=None` and selects `FAILED_NO_REF_FLANK`; output then accesses `fc['expected_with_DNA']`. Actual: `TypeError`, header-only output. Expected: the promised failure row. **Fix:** write nullable fields without dereferencing missing flank metadata. Run (temporary file, not kept).

27. **Medium — `scripts/rearrangement_validator.py:76,108`: show-tool failures become absent evidence.**  
    A malformed delta makes real show-coords and show-diff exit 255. Both parser wrappers return `[]`. Expected: an explicit tool error, preventing a negative verdict. **Fix:** check subprocess return codes and propagate stderr. Parser run (temporary file, not kept), exit-code verification (temporary file, not kept).

28. **Medium — `scripts/collect_confirmed_events.py:42,142,338`: generated records are not uniformly valid VCF 4.3.**  
    Observed outputs include `<TRA>` as a nonpermitted first-level symbolic SV type; `VAL_METRIC` declared `Number=1` but parsed as two comma-separated values; and literal `REF=R`. Expected: valid symbolic types, matching cardinalities and permitted REF bases. **Fix:** use a valid uncertain breakend representation, percent-encode metadata delimiters, and convert ambiguous REF bases according to the specification. These requirements follow [VCF 4.3 §§1.4.5 and 1.6.1](https://samtools.github.io/hts-specs/VCFv4.3.pdf). Generated records and pysam parsing (temporary file, not kept).

29. **Medium — `scripts/collect_confirmed_events.py:223,227`: an evidence directory can contain another event’s report.**  
    First collect target GCA_2 with regenerated evidence; then collect GCA_3 into the same output without regeneration. Actual: `event.json` names GCA_3 while the dnadiff report retains exactly its previous SHA-256. Expected: evidence corresponding to the current event, or absent evidence. **Fix:** use content-derived stable event IDs and replace/check evidence directories against input provenance. Run (temporary file, not kept).

30. **Medium — `scripts/collect_confirmed_events.py:252`: missing target regions do not prevent publication.**  
    Retain validation/details but supply a target FASTA lacking both named contigs. Actual: a confirmed TSV/VCF record is still written with `n_target_contigs=0` and empty target evidence. Expected: explicit extraction failure and no publish-ready evidence claim. **Fix:** require all necessary target regions before emitting the event. Run (temporary file, not kept).

31. **Medium — `scripts/dedup_translocations_per_locus.py:130`; `scripts/cluster_long_is_reps.py:158`; `scripts/collect_confirmed_events.py:204`: counts still confuse rows/accessions with independent support or events.**  
    GCA/GCF twins yield `n_supporting_targets=2`; the long-IS cluster reports `n_distinct_asms=2` for the same twin pair. Repeating one validation row produces two deposited event IDs. Expected: one independent genome/supporting event in these examples. **Fix:** distinguish raw rows, normalized genome identities and event identities, and deduplicate exact repeated validation evidence. Twin counts (temporary file, not kept), collector duplication (temporary file, not kept).

32. **Medium — `scripts/duplication_tandem_db_scan.py:201`: reported identity is the configured threshold.**  
    Real self-alignments of identical repeat units have 100% identity; the table writes `identity=95.0`. Expected: measured identity. **Fix:** retain and report the supporting alignment identity rather than `args.min_identity`. Run (temporary file, not kept).

33. **Medium — `scripts/prepare_review_cases.py:102`; `scripts/inversion_focal_pair.py:251`; `scripts/composite_transposon_review_pair.py:141`: review annotations retain a one-base coordinate error.**  
    All three treat a 1-based inclusive IS start as a zero-based feature start. A 1,500-bp IS is annotated as 1,499 bp in both pair-generation scripts. The small helper example extracts `CGT` instead of `ACGT`. Expected: the complete IS interval. **Fix:** convert starts to `start-1` before slicing/feature construction, including reverse-complement remapping. Executed helper and CLI outputs (temporary file, not kept).

34. **Medium — `scripts/is_large_element_movement_db.py:71`: a successful empty rerun preserves stale results.**  
    Run once with one element, then rerun with an empty records table and the same output directory. Actual: exit 0, “nothing to do,” and the old movement table remains unchanged with one row. Expected: an empty current-run table or explicit invalidation. **Fix:** replace outputs even when no candidates remain. Other movement entry points have analogous early returns. Verified large-element run (temporary file, not kept).

35. **Low — `scripts/hmm_filter_is110.py:190`: `n_orfs` is not the predicted ORF count.**  
    The synthetic sequence has three predicted ORFs; output reports `n_orfs=1`, counting only ORFs with domain hits. Expected: three predicted ORFs, with hit-bearing ORFs reported separately if useful. **Fix:** count from ORF prediction metadata. Run (temporary file, not kept).

(a) Commit verdicts follow. “Correct” here applies to the stated fix, not certification of the entire script.

| Commit | Verdict | Assessment |
|---|---|---|
| `0ffb3a4` | **New error / incomplete** | Strand handling, census options and basic clipping fixes work. New reference-site refinement suppresses real replacements and misclassifies short occupied sites; called-boundary filtering and twin handling remain incomplete. |
| `cc1182e` | **Incomplete** | Shared parsing and transposase-oriented flanks work. Coverage-based checks still accept inverted paralogs; small-span handling is incomplete. |
| `29ca058` | **New error / incomplete** | Orientation normalization and query-side dnadiff handling work. The new two-anchor gate excludes actual duplication catalogue rows; event deduplication and reciprocal linkage remain unsound. |
| `23a520f` | **New error / incomplete** | The same-segment requirement improves the old block-count test, but the new depth helper can return zero for two copies. Hash-seed dependence and missing reference-copy comparison remain. |
| `9f9d911` | **Incomplete** | It queries flank DNA rather than the IS body, but still does not establish the claimed canonical topology or a copy-number gain. |
| `2bfd2d0` | **New error / incomplete** | Full headers, species-qualified IDs and real reference-base fetching work. `<TRA>` is invalid; metadata cardinality, evidence reuse and event identity remain problematic. |
| `9270f80` | **Correct for checked scoped claims** | Deterministic Stage 7 pair selection/schema and the exercised tandem normalization changes behave as described. Its unrelated stage 8–11 claims were not comprehensively audited here. |
| `308a6a0` | **Incomplete** | The advertised 16 positive orientation cases pass. Fixed-region coverage still permits paralog-based false confirmation and erases small inversions. |
| `94997f2` | **New error / incomplete** | Clean 2/3/4-copy arrays are counted correctly. Twin selection can discard all IS annotation; identity reporting is incorrect. |
| `a628f41` | **New error / incomplete** | Strand-aware signed flank projection works in exercised cases. The new unclassifiable-row path crashes; HGT aggregation merges remote fragments; origin occupancy remains inferred from length alone. |
| `a7842d2` | **New error / incomplete** | Nucmer failures raise and vacuous anchor terms were removed. The restored IS lookup now exposes reference/query coordinate mixing and false `IS_MEDIATED` labels. |

(b) Checks that held:

- Both `tests/test_rearr.py` and `tests/test_movement.py` completed successfully.
- Clean inversion orientation controls passed from **both IS_a and IS_b**, including reversed target contigs; unchanged controls were rejected by the gold validator.
- Shared `parse_details` preserves pipes, colons, hyphens and dots. Minus-strand proxies, clipped-window offsets and the simple minus-strand projection control were correct.
- The census options recovered all 76 tested placements on both strands, including a 56.25%-coverage copy.
- Both clustering scripts kept nested 6-kb and 12-kb bodies separate with the new coverage options.
- The canonical movement verifier accepts the DB output schema and rejects reuse of a single IS hit for both ends.
- Artificial-profile HMM controls required both profiles on the same ORF; profiles split across ORFs did not pass.
- Ordinary emitted VCF fields have header declarations. Colons in contig names are legal under VCF 4.3 and were not counted as an error.

The tests leave the demonstrated defects undetected: they use only an 8-kb movement body, do not test copy-gain relative to an already repetitive reference, do not exercise palindromic-anchor hash seeds, use a duplication fixture schema different from Stage 7’s uncapped catalogue, and contain no collector/VCF, reciprocal-linkage, clustering-filter or HMM duplicate-ID regression coverage.

(c) Limits of verification:

Every scoped file was read and exercised through its CLI or relevant helpers. I did **not** establish how many deposited real-data rows are affected, test a 20,000-copy census, exhaust every parameter combination, or validate biological IS110 specificity. Artificial profiles tested HMM bookkeeping and acceptance logic.

Direct DB-output → mechanism-classifier and pairwise-output → canonical-verifier runs also raise schema `KeyError`s; their intended adapters were not supplied. DB-output → canonical-verifier works. Claims in `9270f80` concerning unrelated stages were outside this standalone-script review.