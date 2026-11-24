# KIT-202 outcome (R2 judged; R3 and the check on the merge running)

Ticket `tickets/kit/KIT-202.md` r1.3 (owner decisions VAL-130). Branch `kit-202-intensity-filter-samples` from
`v0.9-staging` 9453e79b.

## Commits

- f9d7d910: tests only (nothing under `R/` or `man/`): `tests/testthat/test-kit202-intensity-filter-samples.R` (T1 to
  T7), `tests/testthat/test-metab-qc-shared.R` (T8, D7), `kit-202/filter-samples-control.R` (S1, runs in R2).

## R1: run 45326839837 on f9d7d910 (tests only, base none, 1 job, 10.8 job-minutes, ceiling 45)

Observed (`r1/tests.tsv`, `r1/failures.txt`, the two logs):
- `test-kit202-intensity-filter-samples.R`: 123 passed, 32 failed, 0 errors. `test-metab-qc-shared.R`: 18 passed, 2
  failed, 0 errors. Every failure is an assertion; none is an error.
- Values at which the tests fail, each as the ticket states for the base:
  - T1 (a) metabolite seed: threshold 4.9 (head 16410), kept 1 to 8 (head 1, 2, 3, 4, 7, 8). G0's value.
  - T1 (b) lipid: threshold 2.2 (head 100.9), kept {3} (head {2, 3}). The kept set was hand-computed in the ticket.
  - T2, both variants, both omics: kept as T1's base sets.
  - T3 (a): the helper refuses `sample_columns` (result NULL), both omics.
  - T3 (b): no `omics_intensity_filter_samples_guessed` warning; the metabolite case keeps 5 to 8 (the diff prints
    `actual[1:3]` 5 6 7 against `expected[1:7]`), the lipid case keeps {3}. Hand values.
  - T5 (a): no `omics_intensity_filter_samples_invalid` error from either S4 method (absent and text sample) or from the
    partitioned step (text sample).
  - T6 (i): the marker row is absent. T6 (iii) and (v): no marker in the record (in process; managed, lipid and metab).
  - T6 (iv): non-vacuity fails (both rules give the same kept set on base), and record (b) and its replay carry no
    marker. The replay itself has no mismatch on either owner, and record (a) carries no marker, on base.
  - T7: kept 3 to 8 (head 1, 2, 3, 4, 7, 8). Hand value.
  - T8: the mocked helper receives `sample_columns` NULL (head c("S1", "S2")), and the text-sample assay is not refused.
- Passing on base, as the ticket expects: T1 (b)'s guard at 0.5, T4 (six unaffected objects against the oracle at
  percentiles 1 and 10 and proportions 0.2 and 0.5), T5 (b), T6 (ii).

Inferred: R1's stated base values hold, so the ticket went on to the fix commit (r1.3 decoupled it from KIT-190).

## Fix and R2

Commits: 2622fcff (true merge with `v0.9-staging` 9f3b7c49; no 7.3 file changed except the marker list, so no R1b);
6e02de15 (the fix); 5d322d9c (S1 instrument repair, docs only).

**R2, run 45469464049 on 6e02de15**, base e305fab3 (staging head when sent; it differs from the merged 9f3b7c49 only in
`tests/testthat/test-omic-enrich-module.R`, KIT-205's re-freeze, and the ledger; neither in R2's scope, `R/` identical).
363.1 job-minutes (3 + 3 shards, check on both sides, scripts), ceiling 675. Observed (`r2/compare.tsv`):
- 39 files `same`; `test-kit202-intensity-filter-samples.R` `only-on-head`, 163 passed, 0 failed; `test-metab-qc-shared.R`
  `count-changed` (11 -> 12 passed in the edited test), no failure: as 6.2 predicts.
- `test-workflow-omics-declarations.R` `message-changed` (base 46 / 22 failed, head 45 / 23): the head adds one failing
  expectation, `expect_identical(current$registry, frozen$registry)`, whose differences are the intensity filter
  declarations' `fixed` (a list on head, NULL frozen) and `global_statistics` text, both omics (`r2/moved-declarations-registry.txt`).
  Decided by the owner (VAL-181): no exception by file; `fixtures/kit073-omics-declarations-frozen.rds` is deferred to
  the end-of-project re-record; it does not stop the ticket.
- `test-workflow-omics-binding-steps.R` `message-changed` (58 / 1 failed on both sides), **not predicted**: the one
  failing expectation (KIT-094's characterisation against `fixtures/kit094-characterisation.rds`) fails on base for input
  sha256 values only; on head it also shows `steps[[3]]$parameters` of length 3 against 2, both omics: the intensity
  filter step's record carries the marker entry (`r2/moved-binding-steps.txt`). 6.2 predicted no record snapshot moves;
  this one does. By 6.2's rule the ticket stopped before merge. Base message: the four `inputs[[i]]$sha256` differences
  (both omics), then `records`. Head message: the same four, plus `actual$<omic>$steps[[3]]$parameters` length 3 against
  2 and `parameters[[3]]` a list against absent, both omics. Predicted by the ticket: no. Decided: VAL-189 (the coordinator, applying the standing re-record rule and the owner's VAL-181 to the same kind of
  stored reference; told to the owner on 2026-11-25): the fixture is not edited and is deferred to the end-of-project re-record, as VAL-181; it does not stop
  the ticket.
- Check: findings identical on both sides (0 errors, 3 warnings, 3 notes each); none only on head (`r2/check-compare.txt`).
- S1 died in R2 (instrument defect: the child's loop overwrote its output-directory variable; no rows written, no
  result). Repaired in 5d322d9c (VAL-54 instrument run 1 of 3; no expected value, comparison or stop rule changed).

**S1, run 45482129676 on 5d322d9c** (base argument e305fab3), 16.3 job-minutes, ceiling 75. Every criterion equal
(`s1/compare.tsv`, `s1/filter-identity.tsv`): (i) unaffected base = head 836/836 (11 assays x 76 settings); (ii) head =
design-sample oracle 228/228 and base = numeric-column oracle 228/228 (3 affected objects x 76); (iii) whole =
partitioned 42/42 on head and on base, unaffected base = head 33/33 per variant; (iv) head under the earlier option =
base 228/228.

Stored references deferred to the end-of-project re-record: `tests/testthat/fixtures/kit073-omics-declarations-frozen.rds`
(`test-workflow-omics-declarations.R`, VAL-181) and `tests/testthat/fixtures/kit094-characterisation.rds`
(`test-workflow-omics-binding-steps.R`, VAL-189).

## Check on the merge (R2b, audit I02) and staging after 393ebb85

R2b's check-only run, b20af066 against 393ebb85, declared in `runs.tsv` before it was sent: result below when back.
S1 is not repeated: KIT-167 (merged at d593e351) changed no file of 7.3 and no lipid or metabolite filter, declaration
or design-save code; its `R/mod_omic_import.R` hunk runs only behind `referenceMappingImportOffered` (the gate on, and
proteomics or transcriptomics).

Staging moved on to 0cf0e440 (KIT-208) while R3 ran against 393ebb85. `git diff --stat 393ebb85..0cf0e440 -- R/ tests/`:
`R/func_general_enrichment_engine_helpers.R`, `R/func_general_enrichment_orchestration.R`,
`R/func_prot_da_results_methods.R`, `R/utils_workflow_dia_da.R`, `R/utils_workflow_enrichment.R` and three test files
(`test-general-enrichment-shared.R`, `test-kit208-gprofiler-identifiers.R`, `test-omic-da-module.R`). None is a file of
7.3, a lipid or metabolite filter, import, design-save or declaration file, or a file of this ticket's tests
(inferred from the diff: R3 stands for the merge with 0cf0e440 unless R3's failing set lists one of these test files).

## Hand-offs (audit I07, I08)

- **KIT-214:** do not call `omicsIntensityFilterSampleColumns` as-is for the duplicate step. If it did, a replay under
  the intensity filter's earlier option would also switch the duplicate step's columns, and its refusal would name the
  wrong step. Give the definition a marker argument and a step label, or write the duplicate step's own call. Then
  KIT-214's T6 (ii) discriminates on its own option.
- **KIT-217:** it edits the same hunks and adds `intensity_filter_missing` to the same declaration. The machinery sets
  one option per fixed marker (`workflowReplayRuleOptions`), so the two rules replay independently. The KIT-073 registry
  and the KIT-094 snapshot will move again; cite VAL-181 and VAL-189, or add a row.
- **The marker's position** in `workflowReplayRuleMarkers` does not matter: the list is read only by name, so a merge
  conflicts only as text. Whoever merges second keeps all rows; each ticket's T6 (i) catches a lost row.

Divergence from `conventions/R.md` recorded (audit I11): 43 added R lines exceed 100 characters, the files' own idiom.

## Acceptance criteria

- A1: met. Base side R1; head side R2: every new test passes (163), T8's file passes.
- A2: met (S1, run 45482129676).
- A3: met with two declared, decided differences: `test-workflow-omics-declarations.R` (stated in the ticket's
  6.2 from r1.3; VAL-181) and `test-workflow-omics-binding-steps.R` (predicted by the ticket: no; VAL-189, the coordinator's ruling told to the owner), both
  `message-changed` by the marker and declaration text only. Every other verdict as 6.2. No `compare-error`, no
  `new-failure`.
- A4: met (check findings identical).
- A5: pending (R3 on the merge with staging 393ebb85, d593e351; R1b not needed: the merge changes no 7.3 file and no
  lipid or metabolite filter, import or design-save file).
- A6: met (`git diff 9f3b7c49..6e02de15 -- R/` touches only 7.3's files; the comparison lines unchanged).
- A7: met (NEWS, ledger row `VAL-<merge>`; `DESCRIPTION`, `NAMESPACE`, lock files unchanged; `man/` as 7.3).

## Not checked

- R3 (full suite) and the merge with the current staging.
- Whether other stored record snapshots outside R2's scope carry the intensity filter step's parameters (R3 would show).
- The managed worker route on an affected table (T6 (v) runs the worker on the text-id fixtures; it shows the marker only).
- Import Existing Design with numeric non-sample columns, or with a run read back as not numeric (no test or script runs
  that route).
- Which MS-DIAL releases write integer identifiers; real data.

## For the change report

**What was wrong, who was affected.** The lipidomics and metabolomics intensity filter removes a feature that is below a threshold in at least the chosen proportion of samples. The threshold is the chosen percentile of the assay's intensities. Both the percentile and the count took every numeric column of the assay table as a sample. A numeric column that is not a sample therefore entered both: an integer identifier, a numeric annotation, retention time, m/z, a score, or a run that is in the table but not in the design. With m such columns beside n samples, a feature's share became (k + a)/(n + m), where a ≤ m of its non-sample values were below the threshold. The threshold was a percentile of a mixture of intensities and other numbers. Identifiers and retention times are small next to intensities, so they pulled the threshold down, and intensities seldom fell below it. A feature could also be removed because its own identifier was below the threshold, so removal could depend on its row number.

Affected:
- MS-DIAL imports, lipid and metabolite, with their default identifier column ('Peak ID' or 'Alignment ID', whole numbers). Shown on the application's own import path by KIT-198's G0, run 45314460795.
- Custom-format imports whose chosen identifier or annotation column is numeric.
- Designs loaded with Import Existing Design whose tables hold numeric non-sample columns (inferred from the code, not run).
- The metabolomics workbooks, which keep every vendor column (shown, T7).
- Objects built in code.
- Scripts calling `lipidIntensityFilteringHelper()` or `metaboliteIntensityFilteringHelper()`.
- The KIT harness's five metabolomics corpora (`Alignment ID` 0 to p − 1).

Not affected: tables whose numeric columns are exactly the design's samples. These include LipidSearch imports, custom imports with text identifier and annotation columns, and the four GUI-entry lanes (shown, S1); KIT-165's and KIT-166's inputs are inferred to be unaffected.

**Code.** Base 393ebb85 (the same lines at 9f3b7c49):
- the S4 methods' sample columns: `R/func_lipid_s4_normalization_methods.R:929`, `R/func_metab_s4_qc_methods.R:84`;
- the exported helpers: `R/func_lipid_qc.R:47`, `R/func_metab_qc_filtering_helpers.R:20`;
- the partitioned writer, the declared step's default variant: `R/utils_workflow_omics_declarations.R:154-159`.

Fix 6e02de15:
- `omicsIntensityFilterSampleColumns` (`R/utils_workflow_omics_declarations.R:160`): the design's samples less the id column, in table order; every numeric column under a record's earlier rule;
- `omicsIntensityFilterRefuse` (`:133`) and `omicsIntensityFilterHelperColumns` (`:175`);
- the S4 methods: `R/func_lipid_s4_normalization_methods.R:912`, `:931`, `:964`; `R/func_metab_s4_qc_methods.R:67`, `:86`, `:119`;
- the partitioned writer (`R/utils_workflow_omics_declarations.R:223-231`) and the page filters (`R/utils_workflow_lipid_declarations.R:97`, `R/utils_workflow_metab_declarations.R:85`);
- the helpers' new argument `sample_columns` (`R/func_lipid_qc.R:50`, `R/func_metab_qc_filtering_helpers.R:23`);
- the replay marker `intensity_filter_samples` (`R/utils_workflow_analysis_record.R:659`), fixed by the step's declaration (`R/utils_workflow_omics_declarations.R:204`).

**What changed in results.** Tables without such a column give byte-identical results (S1 (i): 836 of 836 comparisons, 11 assays at 76 settings; S1 (iii): 33 of 33 per variant; the four GUI-entry test files `same` in R2).

On tables with such a column, the threshold is now the percentile of the design's sample intensities and the share counts the design's samples only. The features kept change, and so does everything downstream (normalisation, RUV, differential abundance, integration). On the fixtures:
- MS-DIAL metabolite seed (8 features, 4 samples), 10th percentile: threshold 4.9 becomes 16,410; features removed at 0.5 go from 0 to 2.
- MS-DIAL lipid table (3 features, 4 samples): threshold 2.2 becomes 100.9; at 0.2 the kept set goes from {3} to {2, 3}.
- Workbook shape (id, retention time and m/z kept): threshold 2.725 becomes 16,410; at 0.25 the kept set goes from features 3 to 8 to {1, 2, 3, 4, 7, 8}.

Over S1's 76 settings per fixture, the kept set changed at 46, 21 and 42 settings. The new rule kept fewer features at 32, 12 and 26 settings, because the threshold rises to the intensity scale. It kept more at 11, 9 and 12 settings, because a feature no longer counts its own identifier or retention time as a value below the threshold. The largest change was 5, 2 and 4 features (computed from `s1/filter-identity.tsv`). Sizes on real data were not measured.

Analyses recorded before this change replay with the earlier rule: their step carries no marker (T6 (iv); S1 (iv), 228 of 228). New records carry the marker. A design sample that is missing from the table or not numeric is now refused. In the application this can happen only through Import Existing Design, when the folder's column mapping lists a run whose column holds no value at all. Without `sample_columns`, the exported helpers count every numeric column except the identifier, and warn.

**Evidence and limits.**
- R1 45326839837: 32 + 2 assertions failing at their stated base values on the tests-only commit f9d7d910.
- R2 45469464049: 43 files, base e305fab3 against head 6e02de15; check identical.
- S1 45482129676: every criterion equal.
- R3 `<R3 id>`: full suite on the merge with 393ebb85.

Fixtures only: the repository's MS-DIAL metabolite seed, the module-CI lipid table with integer ids, the GUI-entry lanes and the module-CI QC objects. No real data. That current MS-DIAL releases write integer identifiers was not checked against a named export. The Import Existing Design route and the managed worker were not run on an affected table; they use the same code (inferred).
