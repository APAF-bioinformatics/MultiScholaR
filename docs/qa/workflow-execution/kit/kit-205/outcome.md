# KIT-205 outcome: the enrichment background is the proteins tested in each comparison

Ticket r1.2 (`tickets/kit/KIT-205.md`; r1.2 applies the re-check's minors K205-Q01 to Q07 and what merged since: KIT-203
at 241914cb, KIT-208's clause E-7). Branch `kit-205-enrichment-background`, cut from `v0.9-staging` at 241914cb (KIT-203
merged, VAL-157). Every run is listed in `runs.tsv`.

## Commits

| Commit | What |
|---|---|
| e43dca86 | tests only: `test-kit205-enrichment-background.R` (B1 to B9) and 10.2's moved expectations; nothing under `R/` |
| 04c466b4 | tests only: B2's helper runs both comparisons in one `processEnrichments` call (so X (m1), a mask pooled over comparisons, can show); no expected value changed |
| 325f9522 | the fix: `R/`, `man/processEnrichments.Rd`, the four report templates, NEWS, ledger row `VAL-168` |
| cd5d590f | R1 records and R2's scripts (`control-c.R`, `w113-cause.R`, `replay-k205.R`, `mutate-k205.R`), `wording.md`; nothing under `R/` |
| b7bd3af8 | R2 records |
| e1b88862 | merge of `origin/v0.9-staging` at 91e3cdcb (KIT-207 merged, VAL-158): only NEWS and the ledger shared; the ledger row placed after VAL-158 |

## Step 0 and the gate

- The ticket was revised to r1.2 in the main checkout (history row, closure table section 18): Q01 to Q07 with their
  texts; every cited line of KIT-203's files re-read at 241914cb; KIT-203's own moved expectations named in 10.2 (S6,
  S9, S10, R3 A107-03); the legend's file corrected (`R/mod_omic_enrich.R`, not `mod_omic_enrich_panels.R`); the
  label's two further readers named (the completed status's "Significance:" line and the download summary); KIT-208's
  clause in 3 item 5; R2's scripts job at 150 minutes (Q02) and KIT-203's G1 in it (stop rule 4 had no run).
- **The owner has not given his word for the UI re-freeze (Q01).** The default of 6.3 applies: `omicEnrichW1FrozenUi`
  is not re-frozen; its test (a) fails on head at `html`, `help` and `choices` (SD-1's words), and the plan section 11
  row below is written. The W113 cases are not re-frozen either (ED-2's default).
- Q05's check (read at 241914cb, before R1): the label, note and notice never reach the TSV, `@result` or the plot
  inputs: `setEnrichmentPlotInputs` stores the token `significance_rule` (`R/func_general_enrichment_engine_helpers.R:
  227-234`) and the TSV carries the column `significance_rule` with the token (`R/func_general_enrichment_orchestration.R:
  674-677` at 241914cb). No stop.

## What was observed

- **R1 (run 45369191776, e43dca86; 30.6 job-minutes).** The new file: 103 passed, 125 failed, 0 errors (`r1/`). Every
  test fails at the value the ticket states:
  - B1: `@universe` holds the 12 untested proteins (length 120 against 108 tested annotated, k = 12); `family_size` on
    base is 9 against the independent count of 6 over the tested ids (3 terms fall below 10 tested members, the
    10-member term among them).
  - B2 (a), (b): the custom background holds every row's id, for up and down, in both comparisons; (c) passes (the
    shared name stays; its second assertion, an untested row's own name, fails on base as expected).
  - B4: first part passes; second part fails (the row with a q-value and no p-value is in the up list and the
    background, on both routes). B5: first part passes; second part fails (universe 120, not 119).
  - B6 (b): the recorder sees the option unset (`"<unset>"`) for both records; (c) `fixed` is NULL.
  - B7: 97 failures (the marker not read; every paragraph prints the old sentences, "measured proteins" in all 24).
  - B9: "measured" in the label, note, notice and legend.
  - Guards passing on the tests-only commit: B3, B4 first part, B5 first part, B6 (a), B8.
  - 10.2's moved expectations fail at their old values: `test-workflow-analysis-record.R` (names, sources: 2),
    `test-omics-enrich-core.R` (D6 test 1: 3; D6 test 4, `fixed` NULL: 1), `test-go-propagation.R` (T7(c): 1; T11(a): 4;
    T11(d): 12), `test-omic-enrich-module.R` (R3 P-28 names: 1; R3 A107-03 notice: 1), `test-enrichment-significance.R`
    (S6: 3; S9: 4; S10: 1), `test-workflow-omics-binding-steps.R` (the strip list: 1). The other failures in these files
    fail on base already: the clusterProfiler captures since KIT-173, and a lipidomics input digest recorded on the
    laptop in the KIT-094 characterisation.
  - Stop rule 1 does not fire.
  - After R1, B2's helper was changed to run both comparisons in one call (04c466b4): mechanics only, no expected value
    changed. R1's B2 values are from the per-comparison helper; on base both helpers send every row, so B2 fails there
    either way. X (m1) shows the new helper sees a pooled mask.
- **R2 (run 45374523563, cd5d590f, R/ = 325f9522, against the branch point 241914cb; 348.1 job-minutes).**
  - Tests (`r2/compare.tsv`): the new file only on head, 228 of 228. Every other file `same`, except
    `test-omic-enrich-module.R` (`new-failure`: test (a), the UI capture, at `html`, `help`, `choices`; R3 P-28's
    `expect_identical(r3$d6, omicEnrichW113Frozen$gprofiler)`, the W113 marker; `message-changed`: the W113 test, every
    case with a completed state now differing by the marker (provenance `rest` digest; record entry `parameters$
    enrichment_background "fixed ; fixed ; tested"`), and `(b) clusterprofiler`, `(b) kernel_clusterprofiler`, failing
    since KIT-173, now with SD-1's words in the notice, the "Significance:" line and the download summary's digest) and
    `test-omics-enrich-core.R` (`message-changed`: the two W2 clusterProfiler captures, failing since KIT-173, with
    SD-1's notice). All of these are 6.3's rows; 10.2's edited expectations pass. Stop rule 3 does not fire.
  - Check: 0 errors, 3 warnings, 3 notes on both sides, identical (`r2/check-compare.txt`).
  - Control C (`control-c.tsv`, `control-c-verdict.txt`): 26 of 26 rows equal. (a) four cases without untested rows
    (clusterProfiler with direct annotations, two comparisons; with the shipped parent terms; KIT-189's form without
    `raw_pvalue`; g:Profiler, two comparisons): head = base for the TSV SHA-256, `identical()` `@result`, plot inputs
    and `@universe`, and for every argument `perform_enrichment` received. (b) five cases with untested rows (12 rows: 6
    in the 10-member term and 6 spread, or 12 spread; 2 rows with a q-value and no p-value; 2 rows without a q-value
    in KIT-189's form; g:Profiler, two comparisons): head on the table = base on the table without its untested rows,
    byte for byte. The largest |log10 p| change against base on the whole table is 0.50 (direct annotations, up) and
    0.52 (parent terms, up); 4 terms tested by base were not tested by head (parent terms, down: INFERRED, they fell
    below 10 tested members and left the size-filtered family). The 2.5 rows recomputed with `stats::phyper` agree with the
    ticket's (`control-c-phyper.tsv`: 0.0028992, 0.0027978, 0.0055023, 0.00015488, 0.0032837). Stop rule 2 does not fire.
  - W113 cause control (`w113-cause.tsv`, `w113-detail.tsv`, `w113-verdict.txt`): 11 cases on both sides; in the 7 with
    a completed state the provenance differs only by `fixed/enrichment_background` (7 leaves) and the record step only
    by the one entry `enrichment_background` (source `fixed`, provenance `fixed`, value `"tested"`; 42 leaves); 106
    provenance and 468 record leaves equal; 0 other. Stop rule 5 does not fire.
  - A3 (ii) (`r2/scripts/replay-k205.tsv`, `replay-k205-base.tsv`): fixture = KIT-173's T7 project with six proteins
    (P03959, P0A6P7, P0A766, P0A7S3 of `fill`, P00547, P00561 of `p_term`) without values in group A; on base the
    published DA table has 6 untested rows, all 6 annotated (so the shipped DA writes such rows with no p-value end to
    end). The record made on base has no marker and replays on head byte-identical to base (result tables exact, plot
    inputs, every TSV's SHA-256); the replayed provenance omits the marker. A fresh head run on the same input (equal
    DA input digests) carries the marker, differs from base, and equals `processEnrichments()` under the earlier rule
    on the recorded tables without their 6 untested rows. The change: up 73 terms to 65 (largest |log10 p| change 1.37;
    BgRatio denominator 69 to 63), down 53 to 53 (0.43). 14 of 14 checks.
  - UI dumps (`dumps/`, `od5-k203.R` unchanged): the base dump equals the base tree's stored block (test (a) passes in
    the base shard); the two head dumps, separate processes, are byte-identical (SHA-256 1b2cca1a...); the base-to-head
    diff (`dumps/base-head.diff`) moves only the HTML digest, the default's note, the legend and the selector's default
    label, each by "measured" to "tested". Not re-frozen (above).
  - X (`r2/scripts/x-verdict.tsv`): B2 fails under (m1) (5), B6 (b, c) under (m2) (2), B5 under (m3) (2; B4 also, 4);
    the unmutated control passes every test.
  - KIT-203's G1 on head (`r2/scripts/g1-summary.tsv`): match on every list (400 of 400; default 4 lists, 13 calls;
    q-value 92 lists, 12,151 calls; m = 1,249), environment equal to KIT-189's provenance. Stop rule 4 does not fire.
- **R3 (run 45383752326, the merge e1b88862 against staging 91e3cdcb, 20 shards per side, no check; 1,215.4
  job-minutes).** Attempt 1: two base shards (10 and 16) failed while installing the library (package downloads
  failed, then Glimma's install), so 45 files had no base result: an instrument defect with no result (VAL-54);
  `gh run rerun --failed` re-ran those two shards and the results job (instrument rerun 1 of 3; nothing in the request
  changed). Attempt 2 (`r3/compare.tsv`): 462 files `same`; the new file only on head, 228 of 228; `test-omic-enrich-
  module.R` `new-failure` and `test-omics-enrich-core.R` `message-changed`, with failure texts identical to R2's
  (`compare-failures.txt` of both runs differ by 0 lines): the UI capture (a), R3 P-28's W113 copy, the W113 marker,
  and SD-1's words in the clusterProfiler captures failing since KIT-173. The four GUI-entry files are `same` (the
  proteomics and transcriptomics lanes fail on both sides as before). No new failure elsewhere, no timeout, no crash.

- **R2b (run 45404531700; audit K205-I01; head db09ddfa = the merge with staging 8cb536cc plus the I02/I03 text
  repairs; base 8cb536cc; 326.1 job-minutes, 21.74 core-hours of 36.00).** R2b 45404531700 on the merge with staging
  8cb536cc (KIT-196, KIT-159 merged; only the ledger hand-resolved): 6.3's verdicts only; check identical. R3 stands
  (section 8). In detail (`r2b/`): 23 files same; the new file only on head, 228 of 228 (with I02's clause);
  `test-omic-enrich-module.R` `new-failure` and `test-omics-enrich-core.R` `message-changed` with failure texts byte-equal
  to R2's (`compare-failures.txt` differs by 0 lines); check 0/3/3 on both sides, identical. Dump control
  (`dumps/r2b/comparison.md`; nothing re-frozen, with the owner): for the UI block and the five W113 cases `gprofiler`,
  `kernel_gprofiler`, `restore`, `organism_filter/no_mapping`, `organism_filter/mapped`, the base dump equals the stored
  block, the two head dumps are byte-identical, and the only leaves that differ are SD-1's words (UI: `html`, `help`,
  the selector's default label) and the marker (W113: `provenance/rest`, `record/parameters$enrichment_background`).

## Acceptance criteria

| # | Met? | Evidence |
|---|---|---|
| A1 | met | R1: every new test fails at its stated value; guards B3, B4 first part, B5 first part, B6 (a), B8 pass; 10.2's expectations fail at their old values. R2: all pass on the fix |
| A2 | met | control C 26 of 26 (`control-c.tsv`) |
| A3 | met | (i) B6 (b, c) passes on head, fails under (m2); (ii) 14 of 14 (`r2/scripts/replay-k205.tsv`) |
| A4 | met | B7: 159 expectations pass (24 clusterProfiler paragraphs, and the g:Profiler bullet in each marker case) |
| A5 | met | R2: only 6.3's moves and 10.2's edits; check identical |
| A6 | met | R3 attempt 2: only 6.3's moves; 462 files same |
| A7 | met (default) | W113 cause control: only the marker; the UI: no re-freeze (no owner's word), the dump diff shows only SD-1's words |
| A8 | met | NEWS, ledger, `man/processEnrichments.Rd`; no version, `NAMESPACE` or lock change; this file |

## Plan section 11 rows (default; the coordinator writes them)

- "| KIT-205 | `omicEnrichW113Frozen` provenance/record of every case: the `enrichment_background` fixed entry (cause
  control `kit-205/w113-cause.tsv`) | on merge | from KIT-205 R3 |". Also `test-omic-enrich-module.R` R3 P-28's
  comparison with `omicEnrichW113Frozen$gprofiler` (`:4321` at the fix). The non-clusterProfiler W113 cases pass on base
  in R2's base shard (ED-2's precondition holds); the owner's word was not given, so they are not re-frozen.
- The KIT-179 row for `omicEnrichW1FrozenUi` gains "and KIT-205's "tested" words (SD-1 (b))".

## Recorded divergences and notes

- The 2.5 rows are in `control-c-phyper.tsv`, not in `control-c.tsv` (one schema per file).
- "A record that has an enrichment step" is read by the templates as `record_present` without `enrichment_background`
  (ticket 3 item 5, r1.2).
- KIT-208's clause (E-7) is printed on every g:Profiler paragraph. It describes the app's tables, which always carry
  `gene_name`; a script call of `processEnrichments` with accessions as the id column on g:Profiler is not described
  by it (KIT-208's matter, not changed here).
- The replay re-executes and never compares results with the record (SD-3): the old-record guarantee rests on B6, X (m2)
  and A3 (ii).
- Line numbers at the fix (325f9522): `enrichmentTestedRows` `R/func_general_enrichment_orchestration.R:150`; the
  reductions `:279` (g:Profiler), `:574` (clusterProfiler); the background `:331`, `:591`; the marker
  `R/utils_workflow_analysis_record.R:654`; the replay wrapper `:922`; the report source `:1169`; the declaration's
  `fixed` `R/utils_workflow_dia_da.R:672`; the templates' marker reading `DIANN_report.rmd:787`.

## Hand-offs (the coordinator sends them)

- **KIT-193** (6.2 item 5 (b), if this ticket merges first): `workflowReplayRuleMarkers` has three rows after KIT-190
  (`enrichment_background`, option `MultiScholaR.enrichment_background`, current `"tested"`, earlier `"all_rows"`);
  the record replay has a second wrapper site, the enrichment step at `R/utils_workflow_analysis_record.R:922` (the
  loop's at `:788`); `processEnrichments` (row T4) changed (the reduction; KIT-203's flag line moved from `:669` to
  `:706`); the docstrings by KIT-193's 5.3 item 4. (Audit K205-I05:) Row T4 is filled on `kit-193-lean-proofs`
  (cd26bb37, digest of `processEnrichments` at KIT-203's merge). KIT-205 changes `processEnrichments`
  (`enrichmentTestedRows` at :279, :574), so the T4 digest is re-computed with `bindings-t4.R` on the merged code and
  pair (ii) is re-run there. The tripwire `test-kit193-bindings.R` fails until then. The statement `QvalueBh.flag` does
  not change: the reduction precedes `enricher`, and q and BH are formed from enricher's own p-values. The citations
  `:669` become `:706`.
- **KIT-208** (E-6): whichever of KIT-205 and KIT-208 merges second runs R2b on the merge. The clause of E-7 is in the
  four templates (KIT-208 need not edit them); after audit K205-I02 it reads "(the first gene name of each protein,
  from the protein annotation or, where the annotation gives none, from the identifier table; proteins without a gene
  name were left out of the lists and the background)". (Audit K205-I04:) If KIT-205 merges first, KIT-208's R2b aligns
  two counts with the tested background: `workflowEnrichmentUnnamedProteins` and `diaEnrichmentNoIdentifiers`' refusal
  count over `enrichmentTestedRows(table)$rows` (which also honours the replay option), and the refusal says 'tested
  proteins' (wording for the owner). Add one case with named untested and unnamed tested rows (refused; n = tested rows)
  and one with an unnamed untested row (not counted). T5, T7 ('25 of 70', (10, 45)) and T8 are unchanged: every protein
  of its fixture is tested. KIT-205's `test-kit205-enrichment-background.R` is in KIT-208's R2 grep scope; B2's recorder
  takes `...`, so KIT-208's two new `perform_enrichment` arguments pass.
- **KIT-200** (8; audit K205-I06): if this ticket merges first, KIT-200's actual interim texts on its branch are the
  ones to remove or replace (the ticket quoted older ones): its NEWS sentence (`NEWS.md:41-42` on the branch) "Until the
  enrichment background is changed in its own update, such a protein is still counted in that comparison's enrichment
  background, as proteins without a p-value are today." is deleted, and its ledger clause "Untested proteins remain in
  the enrichment background until the D-B ticket (VAL-129)" becomes "Untested proteins are not in the enrichment
  background (KIT-205, VAL-<KIT-205's row>)".
- **Plan** (audit K205-I06): section 2 row (ticket 12; 106.28 used of 469.33), order step 4 text (ticket 12), section
  13 rows SD-2 (pool X129, open), SD-3 (X129), SD-4, F2 (ticket 12's table); and the section 11 rows above.
- **KIT-168, KIT-179**: the revision-row text of ticket 12.

## For the owner's wording pass (VAL-41; drafts, nothing final)

Every changed user-facing text, old and new, verbatim (also in `wording.md`):

1. Selector, default choice (also the completed status's "Significance:" line and the download summary): old
   "Benjamini-Hochberg over every term of 10 to 500 measured proteins (default)"; new "Benjamini-Hochberg over every
   term of 10 to 500 tested proteins (default)".
2. Default rule's note, first sentence: old "Each term's p-value is adjusted together with every Gene Ontology term
   (biological process, molecular function and cellular component together) that has 10 to 500 of your measured
   proteins, including terms with none of the listed proteins (they count as p = 1)."; new the same with "tested".
3. Run notice (clusterProfiler, default rule): old "For this organism, terms are called significant by the
   Benjamini-Hochberg adjusted p-value over every term of 10 to 500 measured proteins."; new the same with "tested".
4. Results table legend, second sentence: old "p_adjust_family: Benjamini-Hochberg over every term of 10 to 500
   measured proteins."; new "p_adjust_family: Benjamini-Hochberg over every term of 10 to 500 tested proteins."
5. Reports, clusterProfiler family (default rule): with the rule "terms with 10 to 500 annotated tested proteins,
   whether or not a listed protein falls in them (terms without one enter with p = 1)."; earlier records "terms with 10
   to 500 annotated proteins of that comparison's differential-abundance table, whether or not a listed protein falls in
   them (terms without one enter with p = 1)." (old: "... annotated measured proteins, ...").
6. Reports, clusterProfiler family (q-value rule): with the rule "terms with 10 to 500 annotated tested proteins and at
   least one protein in the list."; earlier records "terms with 10 to 500 annotated proteins of that comparison's
   differential-abundance table and at least one protein in the list." (old: "... annotated measured proteins ...").
7. Reports, clusterProfiler background: with the rule "The background was the proteins tested in that comparison
   (those with a p-value) with at least one Gene Ontology annotation{background_clause}."; earlier records "The
   background was every protein listed in that comparison's differential-abundance table, including any without a
   p-value, with at least one Gene Ontology annotation{background_clause}." (old: "The background was the measured
   proteins with at least one Gene Ontology annotation{background_clause}.").
8. Reports, g:Profiler bullet: with the rule "• The background was the identifiers of the proteins tested in each
   comparison (those with a p-value in its differential-abundance table; the first gene name of each protein, from the protein annotation or, where the annotation gives none, from the identifier table; proteins without a gene name were left out of the lists and the background) that g:Profiler
   recognises (a custom background)."; earlier records "• The background was the identifiers of each comparison's
   differential-abundance table (the first gene name of each protein, from the protein annotation or, where the annotation gives none, from the identifier table; proteins without a gene
   name were left out of the lists and the background) that g:Profiler recognises (a custom background)." (old: "•
   The background was the identifiers of each comparison's differential-abundance table (gene names where the table
   has them, otherwise protein accessions) that g:Profiler recognises (a custom background).").
9. `processEnrichments` help: `significance_rule` "measured" to "tested"; the new Details paragraph (`wording.md`).
10. NEWS: the new entry (ticket 12 with Q03's worked example) and KIT-203's entry "measured" to "tested".
11. Log line (both routes, no notice): "Using %d unique proteins as background for enrichment analysis (%d rows without
    a p-value left out)".

## Not checked

- R CMD check on the merge with 91e3cdcb itself (R2b checked the later merge with 8cb536cc, which contains it).
- Real data and the incidence of untested rows on it (KIT-160 to KIT-162); Windows.
- The g:Profiler service with a custom background (no service call; the route is checked through the arguments sent).
- A non-numeric p or q column (ticket 13).
- That KIT-193's row T2/T4 digests follow (KIT-193 not merged; handed off).
- (Audit K205-I08) Per-comparison masking on the clusterProfiler route is shown once, by control C (b)
  (`b_cluster_none` runs two comparisons with different untested sets in one call); no suite test or mutation covers it
  (X (m1) mutates the g:Profiler site only).
- (Audit K205-I09) B2's failure on base was observed with the per-comparison helper only (R1); with the helper of
  04c466b4 it is INFERRED to fail too (base sends every row either way); X (m1) shows the new helper detects pooling.
- (Audit K205-I10) Two of R2's 26 scope files run 0 tests (Q-163 quarantine: `test-prot-11-enrichment-module-contracts.R`,
  `test-prot-11a-enrichment-module-characterization.R`), so 24 files are informative.
- (Audit K205-I11) Ticket 3 item 1's "every background row could enter a list" holds only where `log2FC` is finite on
  tested rows: true of package tables (limma gives a finite logFC wherever p is finite), not of every scripted table; a
  tested row with an NA `log2FC` stays in the background and in no list.

## For the change report

**What was wrong, who was affected.** Proteomics enrichment (both the g:Profiler route for its 13 organisms and the
clusterProfiler route for every other organism) took the background (universe) of each comparison from every row of
that comparison's differential-abundance table, including proteins that had no p-value there (no value in a compared
group; after KIT-200 also none in the analysed samples). Such a protein could never be in an up or down list, yet it was
counted in the universe, which changes the hypergeometric test's N and K with no chance of selection: the test is biased
in either direction and can exceed its nominal level. In the ticket's worked example (a 50-protein term none of whose
proteins was untested, 300 of 3,000 proteins untested, lists of 150 with no real change) the old background called the
term on 5.6 % of lists at the 5 % level (exact); with 20 of the term's 50 untested it lost almost all its power (0.52 %).
The reports said the background was "the measured proteins". Affected: every proteomics enrichment comparison with at
least one annotated untested protein (on the default chain, designs of three or more groups, or samples removed after
protein QC; DIA with limpa has none).

**Code.** Base 241914cb: `custom_bg` from every row, `R/func_general_enrichment_orchestration.R:301-302`;
`background_IDs` from every row, `:556`, passed as `universe` at `:563` and `:580`. Fix 325f9522:
`enrichmentTestedRows` (`:150`) applied per comparison before the lists and the background on both routes (`:279`,
`:574`); the replay rule marker `enrichment_background` (`R/utils_workflow_analysis_record.R:654`), fixed by the
enrichment declaration (`R/utils_workflow_dia_da.R:672`) and applied around the record replay's enrichment step (`:922`);
the report source (`:1169`) and the four templates' sentences by marker (`DIANN_report.rmd:787` and the paragraph
below it, likewise in DIANN_limpa, LFQ and TMT); KIT-203's texts "tested".

**What changed in results.** A comparison in which every protein has a p-value gives byte-identical results (control
C (a)). A comparison with untested proteins gives exactly what the earlier code gives on its table without them
(control C (b)): every p-value, both Benjamini-Hochberg columns, the q-value, KIT-203's whole-family value and its family
size, and so the flags can change. On the fixtures: largest |log10 p| change 0.50 (12 of 120 untested, direct
annotations) and 0.52 (with parent terms), some terms leaving the tested set; on KIT-173's project with 6 of 70
untested, up list 73 to 65 terms, largest |log10 p| change 1.37. Analyses recorded before this change replay with every
row, byte for byte, and their reports describe that background; new records carry the marker.

**Evidence and limits.** R1 45369191776 (failing first), R2 45374523563 (control C, W113 cause control, A3 (ii), X, G1,
regression scope with check), R3 45383752326 (full suite on the merge). Fixtures only: KIT-173's committed E. coli
UniProt annotation (2026_04) and T7 project; no real data, no g:Profiler call (its route is checked through the
arguments sent to the service). The size figures are exact hypergeometric computations, not simulations.

**Exact sizes (audit K205-I07).** Exact sizes of the old test (lists of 150, 10 % untested, none in the term, nominal
5 %): 5.6 % for a term of 50, 8.4 % for 200, 10.8 % for 400; at 5 % untested 4.5 %, 5.7 %, 6.4 % (exact hypergeometric,
stdlib; audit K205-I07). At 2 % untested there is no inflation (3.96 %, 4.53 %, 4.56 %). These are exact computations for
those settings (lists drawn at random from the tested proteins, the old test with N = 3,000), not simulations; under the
tested universe each is at most 5 %.

**Remaining limit.** On g:Profiler the domain is `custom`, so tested proteins that g:Profiler recognises but that carry
no annotation in a source still count in that source's background; the GO route counts annotated proteins only (SD-2,
pool X129, open).

## Budget

Core-hours (job minutes x 4 / 60): R1 2.04 (30.6 min), R2 23.21 (348.1 min), R2b 21.74 (326.1 min; outside the 469.33 ceiling, as ticket 11 allows; its own ceiling 36.00), R3 81.03 (1,215.4 min, of which the
instrument rerun of two base shards and the results job 60.1 min); total 106.28 within the 469.33 ceiling (128.02 with R2b) against 88.43 expected
(R3's row ceiling 186.67; R2's 46.00). Instrument reruns (VAL-54): R3's two base shards, 1 of 3. No correction
pass was used.


## Merge record (coordinator)

Merged into `v0.9-staging` on 2026-11-25 from 85816589 after the implementation audit and its re-check 1 (register `AUDIT-KIT-205-impl-20261125`: mergeable now). The ledger row is VAL-168. R2b (run 45404531700, with R CMD check, on the merge with staging 8cb536cc) and the full suite R3 (45383752326) stand. Since 8cb536cc staging gained the validation harness's records and KIT-206 (proteomics session export and Summary archive, merge 9e59b90a): KIT-206 and this ticket share only NEWS and the ledger; no run was made on this merge. **Known failing tests on staging from this merge, by design, until the owner allows the re-freeze (decisions pool X163):** 9 expectations in 3 tests of `tests/testthat/test-omic-enrich-module.R`: test (a) of the tab's stored screen at `:209` (html, help, choices), the five stored-record cases at `:4201` (`gprofiler`, `kernel_gprofiler`, `restore`, `organism_filter/mapped`, `organism_filter/no_mapping`) and the comparison at `:4321`; the file's counts go from 1051 passing / 19 failing to 1042 / 28. The control for the re-freeze is in `dumps/r2b/`: the base dump equals the stored block, the head dumps are byte-identical, and only the declared words and the one rule marker move; the stored screen block is lines 138 to 206. Correction to `dumps/r2b/comparison.md` (re-check finding K205-R1-01): outside the screen block the two head files differ in three lines (2009, 2043, 2113), all the `gprofiler_plot` digest of the refusal sessions (`pre_da/run_before_da`, `pre_da/da`, `pre_da/open`), the known placeholder-plot flake; not in `kernel_gprofiler` and not one line. KIT-208 merges after this ticket and aligns its two counts to the tested rows. Ticket total 128.02 core-hours including R2b.
