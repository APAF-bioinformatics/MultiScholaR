# KIT-214 outcome

Ticket `tickets/kit/KIT-214.md` r1.5 (owner decision VAL-173 (1); coordinator E-1 to E-8; r1.4 K202-I07; r1.5 the
marker's snapshot prediction, written before R2). Branch `kit-214` from
`v0.9-staging` 74b1caee.

## Step 0 (before R1)

- Every line number and statement of the ticket's sections 2, 3, 6, 8 and 10 re-read at 74b1caee. 74b1caee differs
  from 9f3b7c49 (the ticket's reading) only in `docs/methods-decision-ledger.md`; every citation holds. The register's
  RC1-01 to RC1-03 are in r1.2 (verified, not re-applied).
- One contradiction with the code, written into r1.3: G7's base column said "same", but RC1-02 made its head
  expectation "`Samples` counts it" (the text design sample). On base `Samples` counts the numeric columns
  (`R/mod_omic_qc_tabs.R:662`), so the text sample is not counted. G7 is split into (a) `Samples` (fails on base by
  assertion; added to R1's failing-first list and A2) and (b) the rest (passes on both).
- R2 scope grep (section 10) at 74b1caee: 33 files with `workflowReplayRule`, 32 without (r1.0 recorded 32 at 91e3cdcb
  without it).

### Step-0 survey (ticket section 6)

Searched `tests/testthat` at 74b1caee for (1) mocks or injected functions with a fixed formals list of a function whose
call this ticket changes (the two duplicate kernels, `resolveDuplicateFeaturesByIntensityFn`, the ITSD builders, the
`norm$itsd_table` and `qc$itsd$metrics` lambdas), and (2) tests that run C1's declared step,
`resolveMetabDuplicateAssayData` or a QC display on an object with a numeric non-sample column against a stored value.

- `test-lipid-01g-qc-duplicates-module-characterization.R:153-157`, `:355-358` and
  `test-metab-01u-qc-duplicates-module-characterization.R:177-181` mock the kernel with `(assay_tibble, id_col,
  sample_cols)` and assert `sample_cols == c("Sample1", "Sample2")`: their designs' samples, so D, and the kernel's call
  shape is unchanged under E-7 (B). Their sample columns are double and their ids text, so the callers' conversion of
  every numeric column is a no-op on what the mock receives.
- `test-lipid-01c-qc-duplicate-module-seams.R:577` mocks the lipid kernel likewise: the file is quarantined (skip on
  line 1).
- `test-metab-01p-qc-duplicates-module-characterization.R:141-144`, `test-metab-01u-...:310`,
  `test-omics-audit-carry-overs.R:187-190`, `test-module-ci-metab-qc.R:80`, `test-omics-qc-core.R:134` call
  `resolveMetabDuplicateAssayData` without `designSamples` (NULL keeps today's rule).
- `test-omic-gui-config.R:692-696` mocks the metric kernels with `sample_columns = NULL` and calls the lambdas with four
  arguments (`:697-700`); `:933-941` mocks the two builders with `(assay_data, <id>, annotation_cols = NULL,
  is_patterns = NULL)` and calls the lambdas with three arguments. The `...` form passes nothing on those calls.
- The frozen QC displays (`test-omic-qc-module.R` `omicQcFrozen`, `omicQcItsdFrozen`; `test-omics-qc-core.R`
  `qcCoreFrozenLipidItsd`, `qcCoreFrozenMetabItsd`, `qcCoreFrozenFigures`) and `test-workflow-omics-application-reads.R`
  use the module-CI or declarations objects: text ids, every numeric column a design sample (D = N).

Result: no hit. Nothing is put to the owner (VAL-130 (4)).

## Commits

- cd88f95f: tests only (nothing under `R/`, `man/`, `NAMESPACE`): `tests/testthat/test-kit214-sample-columns.R` (E1 to
  E10, T6 (i) to (v), G1, G2, G3, G7 (a) and (b)), `kit-214/g0-base-values.R` (G0), `kit-214/dup-samples-control.R`
  (S1).
- fece3cbe: G0 and R1 records.
- 49011d1f: merge of staging 185ebf75 (KIT-202 in); no conflict.
- 8e0465e7: the fix (13 files: 11 under `R/`, `NEWS.md`, the ledger row `VAL-<train5-2>`). `man/`, `NAMESPACE`,
  `DESCRIPTION` unchanged; the two kernels (`R/func_omics_qc_shared.R:116-160`, `R/func_lipid_qc_support_helpers.R`)
  and the batch adaptors unchanged.
- the R1b/R2 records commit, then bb36097e: merge of staging 67569959 (ledger rows only since 185ebf75; the ledger
  conflict resolved with staging's rows first and this ticket's `VAL-<train5-2>` row after them).

## The fix (8e0465e7)

- **Own helpers (K202-I07, ticket 3 item 8).** `R/utils_workflow_omics_declarations.R`: `omicsDesignSampleColumnsOf`
  (`:202`, rule-free D), `omicsDesignNumericColumnsOf` (D that hold numbers), `omicsDuplicateResolutionSampleColumns`
  (`:227`, its own marker and refusal class `omics_duplicate_resolution_samples_invalid`, ticket 7's text),
  `omicsDuplicateResolutionInput` (`:251`, the whole-matrix callers' conversion of every numeric column and D, only with
  a duplicated id). KIT-202's `omicsIntensityFilterSampleColumns` is neither called nor edited.
- **C1.** Lipid adaptor (`R/utils_workflow_lipid_declarations.R:65`); `resolveMetabDuplicateAssayData` gains
  `designSamples = NULL` (`R/func_metab_s4_qc_methods.R:898`, `:921`), passed by the metabolite adaptor; partitioned
  writer passes D to the batch kernel (`R/utils_workflow_omics_declarations.R:339`) and keeps `as_double` over every
  numeric column; the declaration fixes `duplicate_resolution_samples = "design"` (`:315`); the marker row is last in
  `workflowReplayRuleMarkers` (`R/utils_workflow_analysis_record.R:663`).
- **C2, C3.** `R/mod_omic_qc_tabs.R:621` (count of D), `:666` (`Samples` = D; missingness over D that hold numbers);
  `omicsAssayMissingPercent(scope, columns = NULL)` (`R/utils_workflow_omics_declarations.R:1361`), given D that hold
  numbers by `readWorkflowOmicsShell`.
- **C4.** `analyzeWorkflowOmicsItsd` passes D that hold numbers to `itsd$metrics` and `itsd$quant`
  (`R/func_omics_qc_shared.R:256`); the two metric lambdas gain `...`.
- **C5.** The two builders gain `sample_columns = NULL` (`@noRd`); the two `itsd_table` lambdas gain `...`;
  `omicNormItsdTables` passes D that hold numbers per page (`R/mod_omic_norm_panels.R:1258`).
- **C6.** Unchanged (E-3): no caller.
- Departures from the ticket's letter, no result changed on any tested input: without a duplicated id and with no design sample in the table, the whole-matrix caller hands the kernel every numeric column (no new warning; the kernel returns the table unchanged); when no design sample of an assay holds numbers, C4 and C5 keep today's columns (the kernels' and the builder's own rule). A third `@noRd` helper, `omicsDesignNumericColumnsOf`, gives D ∩ numeric. Two new `@noRd` helpers lacked some `@param` lines (R.md:277); added in the audit-repair records commit (roxygen comments only, `@noRd`, no `man/` page, no built file changed).
- **Divergence recorded (N6).** The new code reads S4 slots with `methods::slot()` and `@`, as the surrounding package
  code does, against `R-bioinformatics.md`'s accessor rule (existing-code scope rule).

## G0 and R1: run 45472096981 on cd88f95f (one request: a scripts job and one test job; base none)

`git diff 74b1caee cd88f95f -- R/ man/ NAMESPACE` is empty, so the scripts job ran the base's package code (G0 binds
the base by content). Job minutes: tests 11.6, scripts 9.1 (20.7 in all, 1.38 core-hours; ceilings 45 and 60).

Observed, G0 (`g0/g0-base-values.tsv`, `g0/console.txt`; exit 0, all 23 fixtures built):
- C1 base kept rows (`.row`), `whole` and `partitioned` identical on every fixture: E1 tables row 1 (A), E2 tables row 1
  (A), the tie row 1 (the first in input order: the oracle's tie rule holds), the code-built module-CI
  `duplicate_feature = TRUE` object row 2, through the import also row 2, the seed with a duplicated `Alignment ID`
  row 1, A2 lipid_canonical row 12 (the x 2 copy), A2 metab_lc row 5 (the x 2 copy).
- Vendor columns: the module-CI `duplicate_feature = TRUE` object through the lipid MS-DIAL import keeps only `Peak ID`
  (character), `Name` and the four samples (so S1 classes it unaffected at run time); built in code it keeps `RT (min)`,
  `Precursor m/z` and `Total Score` (numeric).
- C2 / C3: metabolite seed 5 / Samples 5, Missingness 0 (frame and shell); seed with ids 0 to 7: 5 / 5, 2.5; lipid
  MS-DIAL with `Peak ID` 1:3: 5 / 5, 13.3.
- C4 on the seed, pattern "ITSD": `ITSD_Caffeine_d9` mean 38,721.4, CV 55.893 %; `ITSD_Leucine_d10` mean 25,641.6, CV
  55.886 %; "IS with CV > 30%: 2"; long-table samples include `Alignment ID`. Objects without standards record the
  refusal "No internal standards found" (expected; no test reads them).
- C5: E9's table counts its numeric annotation (`Average Mz`) on base.
- Every hand-computed base value of section 10 equals G0's (the leucine mean, not stated in the ticket, is 25,641.6).
  No correction to the ticket's values.

Observed, R1 (`r1/tests.tsv`, `r1/failures.txt`, `r1/test-kit214-sample-columns.log`): 122 passed, 44 failed, 0 errors.
Every failure is an assertion; none is an error:
- E1 (4 cases) kept S1 = 50 (head 60); E2 (4) 10 (head 12); E3 and E9 the argument is refused (result NULL); E4 "5";
  E5 Samples 5 (both paths; Missingness 0 passes); E6 Samples 5, Missingness 2.5; E7 Samples 5, Missingness 13.3; E8
  mean 38,721.4 and CV 55.89 % (and the CV line and long table); E10 no classed error (both omics).
- T6 (i) marker row absent; (iii) and (v) no marker in the record; (iv) non-vacuity fails (both rules keep row A on
  base) and record (b) has no marker; the replays themselves report no mismatch.
- G7 (a): `Samples` 3 (head 4), both paths, both omics.
- Passing on base as the ticket expects: T6 (ii), G1 (six objects, C1 to C5 against the oracles), G2, G3, G7 (b).

Inferred: stop rule 1 does not apply (no E row passes on base); the ticket goes on to the fix commit after KIT-202's
merge, with R1b on the merged state.

## R1b: run 45518385430 on 49011d1f (tests-only state merged with staging 185ebf75)

Observed (`r1b/`): 122 passed, 44 failed, 0 errors; `failures.txt` names the same tests with the same actual values as
R1. Line numbers of section 8's files re-read on the merge (KIT-202 moved `resolveMetabDuplicateAssayData` to `:891`,
the duplicate declaration to `:243`, `omicsAssayMissingPercent` to `:1284`; 2.3's definition is KIT-202's 7.1, so stop
rule 4 does not apply).

## R2: run 45520213660, head 8e0465e7 against base 49011d1f (49 files, 3 shards per side, check, S1)

Observed (`r2/compare.tsv`, `r2/check-compare.txt`, `s1/`); 280.1 job-minutes (ceiling 675):
- `test-kit214-sample-columns.R`: base 122 passed, 44 failed; head 176 passed, 0 failed (`fixed`): every E row, T6 (i)
  to (v), G1 to G3 and G7 pass on head.
- `test-workflow-omics-binding-steps.R` (58/1 both sides) and `test-workflow-omics-declarations.R` (45/23 both sides):
  `message-changed`, counts unchanged, as predicted in r1.5 before the run: the KIT-094 snapshot's duplicate step
  (steps[[4]]) records 1 parameter instead of 0 in both omics (`r2/compare-failures-binding-steps.txt`), the KIT-073
  registry the declaration's `fixed`. Both snapshots were recorded on the laptop and already fail on base; replaced at
  the end-of-project re-record (VAL-181, VAL-189).
- The other 46 files `same`, including `test-kit202-intensity-filter-samples.R` (G5), `test-metab-02f-...` (G6), the
  three mocking files of the step-0 survey, the frozen QC families and `test-omics-norm-itsd-selection.R`.
- Check: 0 errors, 3 warnings, 3 notes on both sides; findings identical (6 on both), none only on head.
- S1 (`s1/compare.tsv`, `s1/dup-identity.tsv`, `s1/changed-rows.tsv`; both sides exit 0):
  (i) base equal to head on 32 of 32 (A)/(A2) rows (both variants); non-vacuous for lipid and metabolite (module-CI
  objects with 2 duplicated ids per assay, the A2 tables with 1); (ii) head equal to the oracle over D on 12 of 12 (B)
  rows, base equal to the oracle over every numeric column on 12 of 12; (ii)'s rows list each input once per variant (whole, then partitioned); the variant is not printed (`setting` blank).
  (iii) head `whole` equal to `partitioned` on 22
  of 22; (iv) head under the earlier option equal to base on 12 of 12. The run-time-classed input (module-CI
  `duplicate_feature = TRUE` through the import) was classed unaffected (its vendor columns are dropped, as G0 showed).
  Changed kept rows (`whole`, base against head): 5 ids: 2.2's two tables in both omics (base row 1, head row 2) and the
  code-built module-CI MS-DIAL table (`P001`: base row 2, head row 1). The seed with a duplicated `Alignment ID` keeps
  the same row under both rules (equal sample counts per candidate), as 2.2 predicts.

Inferred: A2 (head side), A3 and A4 are met; nothing unpredicted, so stop rules 2 and 3 do not apply.

## R3

Run 45573131036: full suite, head bb36097e against staging 67569959, 20 shards per side, 1,270.2 job-minutes
(ceiling 2,800). Observed (`r3/compare.tsv`): 481 files; 478 `same`; `test-kit214-sample-columns.R` `only-on-head`
(176 passed); `test-workflow-omics-binding-steps.R` and `test-workflow-omics-declarations.R` `message-changed` with
counts unchanged, as in R2. No `new-failure`, no other `message-changed`, no timeout or crash. Inferred: A5 met.

## Budget

Runs 45472096981, 45518385430, 45520213660, 45573131036: 1,582.8 job-minutes = 105.5 core-hours (4 cores),
against 103.70 expected and 488.33 ceiling. No correction pass and no instrument rerun.

## Acceptance criteria

- A1: met (G0).
- A2: met (R1, R1b base side; R2 head side).
- A3: met (S1 in R2).
- A4: met (R2), with the two predicted `message-changed` files of r1.5.
- A5: met (R3: no difference other than A4's).
- A6: NEWS and ledger row (`VAL-<train5-2>`, S1's run and count filled) added; `git diff 67569959 bb36097e -- man/ NAMESPACE DESCRIPTION` (the staging merged last) and `git diff 49011d1f 8e0465e7 -- man/ NAMESPACE DESCRIPTION` (the fix) are empty; the kernels' bodies unchanged; this file ends with "For the change report".

## Not checked

- K242-P01 (the declarations and KIT-172 fixtures never apply ITSD normalisation) does not touch this ticket's evidence: S1 and the new tests run the duplicate step and the QC displays, never the ITSD normalisation, and E8, G1 and G7 assert that the QC tab found its standards. That a changed kept row of a standard changes its assay's ITSD normalisation stays inferred; the repository's declaration fixtures could not show it.
- Real data; Windows (no platform code).
- The downstream reach of a changed kept row (ITSD normalisation, RUV, limma): S1 measures the step only.
- Whether any KIT-159 harness record holds QC display values on the metabolomics corpora (C2 to C4 move them: an id of
  0 is no longer a missing cell, the sample count is n, not n + 1).
- The application driven in a browser: the displays are checked through their functions (`omicQcDataSummaryUi`,
  `omicQcAssayStatsFrame`, `analyzeWorkflowOmicsItsd`, the `itsd_table` lambda), not through `omicNormItsdTables`'
  output with an affected object.
- Import Existing Design with numeric annotation columns: not executed; it reaches the same declared step (inferred).

## For the change report

**What was wrong, and who was affected.** In lipidomics and metabolomics, when a feature identifier appears on more than
one row, the duplicate-resolution step keeps the row with the highest average intensity. It averaged every numeric
column of the table, not only the samples. A numeric identifier, a numeric annotation column (retention time, m/z,
score), or the numeric columns kept by Import Existing Design therefore entered the average and could decide which row
was kept. This reached analyses with a duplicated identifier and a numeric non-sample column: custom imports with a
repeated numeric identifier or a numeric annotation column, Import Existing Design, and objects built in code. MS-DIAL
imports at their default identifier (unique) and the metabolomics workbooks (which call the S4 method that already
averaged the design's samples) were not affected. The kept row then feeds every later step (ITSD normalisation if it is
a standard, between-sample normalisation, RUV control selection, limma's shared prior and BH family). Separately, the
Quality Control summary's "Number of Samples", the per-assay `Samples` and `Missingness` (where an identifier of 0 was a
missing cell), the ITSD means and CVs in QC, and the normalisation's ITSD selection table counted every numeric column
(display only): an MS-DIAL import showed one sample too many.

**Code.** In the legacy pipeline (`main` 436cdbe8) the same choices sit in the QC modules: C1 `R/mod_lipid_qc_duplicates.R:261`, `R/mod_metab_qc_duplicates.R:261`; C2 `R/mod_lipid_qc_s4.R:175`, `R/mod_metab_qc_s4.R:175`; C3 `:251`, `:258-263` of both; C4 `R/mod_lipid_qc_itsd.R:156-160`, `:179`, `R/mod_metab_qc_itsd.R:156-160`, `:179`; C5 `R/func_lipid_norm.R:81-82`, `R/func_metab_norm.R:81-82`. There too the design builder kept only the id, annotation and assigned samples (`R/mod_lipid_design_builder.R:1070-1083`), and the S4 `resolveDuplicateFeatures` the workbooks used averaged the design's samples (`R/func_metab_s4_objects.R:4738-4740`). Base 185ebf75 (after KIT-202): `R/utils_workflow_lipid_declarations.R:61`, `R/func_metab_s4_qc_methods.R:906`,
`R/utils_workflow_omics_declarations.R:259-266` (C1); `R/mod_omic_qc_tabs.R:619`, `:662-669`,
`R/utils_workflow_omics_declarations.R:1246`, `:1284` (C2, C3); `R/func_omics_qc_shared.R:243-267`,
`R/utils_omic_gui_config.R:2502-2526`, `:2728-2751` (C4, C5 lambdas); `R/func_lipid_norm_support_helpers.R:38-39`,
`R/func_metab_norm_support_helpers.R:130-131` (C5). Fix 8e0465e7: see "The fix" above.

**What changed in results.** Only the kept row of an identifier listed more than once, and only where a numeric
non-sample column exists: on the synthetic fixtures 5 identifiers changed their kept row (S1, `s1/changed-rows.tsv`);
on these fixtures the non-sample values (an id of 1000, an m/z of 100 or 875) are as large as the intensities or larger. On usual intensity scales the earlier average, which added each row's non-sample values once, lowered the average of a sparsely measured row more than that of a fully measured one. So on such tables, away from near-ties, the new rule changes the kept row only to a row measured in fewer samples whose measured values average higher. For example, with retention time and m/z as the extra columns, 2 of 6 samples at a mean of 210,000 now win against 6 of 6 at 150,000, where before the fully measured row was kept. The rule itself ('highest mean of the measured values') is the owner's open item N1 (plan section 13);
all other inputs, including the GUI-entry lanes, the module-CI objects and real-shaped lane tables with a duplicated
identifier, are byte-identical base against head. Display values: the seed MS-DIAL import shows 4 samples, not 5; its
ITSD caffeine standard mean 48,400 and CV 0.38 %, not 38,721.4 and 55.89 %; the CV > 30 % count 0, not 2; with an
identifier of 0 the missingness 0.0 %, not 2.5 %. The kept row brings its own annotation columns with it: on the code-built MS-DIAL table the reported lipid of `P001`
changes from TG 52:3 (base, row 2) to PC 34:1 (head, row 1) (AUDIT-KIT-214-impl-20261125 I01). Analyses recorded before
this change replay with the earlier rule (the replay marker `duplicate_resolution_samples`). **Open (owner):** whether
the rule itself (N1) should change, (a) keep it, (b) the row measured in most samples first, then the highest mean,
(c) the highest mean with missing values counted in a stated way; any change is a new ticket with its own marker.

**Evidence.** G0 and R1 run 45472096981; R1b 45518385430; R2 and S1 45520213660 (synthetic and repository fixtures
only: 2.2's two tables, the module-CI and GUI-entry objects, the MS-DIAL seed, lane tables with copied rows); R3
45573131036. Limits: no real data; the step measured, not the downstream analysis.
