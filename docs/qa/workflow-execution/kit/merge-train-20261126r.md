# Merge train 25 of 2026-11-26: KIT-255 alone

Branched from `origin/v0.9-staging` at `7f17547f99f5440f1330d7e760722be994460611` (fetched twice, before building and
again before sending; it had not moved; it holds trains 15 to 24, the latest being train 23 and train 24's landing). One
member. The branch `origin/kit-255b` carries the correction pass and, after the closure check, the four record fixes
C-R1 to C-R4.

## Train

Base: `7f17547f99f5440f1330d7e760722be994460611` (staging).

| Order | Ticket | Branch head | Merge commit | Ledger row | Audit |
|---|---|---|---|---|---|
| 1 | KIT-255 | `origin/kit-255b` 0ab6d30d (record fixes on 9a34269d, which holds staging 0ebce762) | 1567f188 | `VAL-<train25-1>` (replaces `VAL-<merge>`) | AUDIT-KIT-255-impl-20261126 (closure check: ready once the four record fixes land; they are in 0ab6d30d) |

Head (bound by the run): `1567f188d706e8dedb3d2760e08fff82e88a83c6`. This file is in the declaration commit after it, under
`docs/` only.

### Conflicts met and their resolutions

- **`docs/methods-decision-ledger.md`** only. Kept every staging row (through VAL-343) in order; KIT-255's row placed after
  VAL-343 as `VAL-<train25-1>` (the placeholder, to be numbered at landing), its status cell now reading
  "merged (<commit>), train 25" (the `<commit>` placeholder is also left for landing). Nothing else in the ledger differs
  from staging.
- No conflict under `R/`, `tests/`, `inst/`, `man/`, `DESCRIPTION` or `NAMESPACE`.
- **`NEWS.md`** merged by itself: KIT-255's one bullet (25 lines) added, nothing removed. It sits where the branch put it
  (after KIT-265's bullet at the end of the proteomics/omics entries before the KIT-259 pair), so every staging bullet is
  kept and staging's own order is unchanged; the coordinator places it in landing order.
- Auto-merged code: `R/func_general_filemgmt_study_parameter_writer.R`, `R/func_lipid_s4_correlation_methods.R`,
  `R/func_metab_s4_qc_methods.R`. `git diff --shortstat origin/v0.9-staging HEAD -- R inst man tests`: 22 files,
  1370 insertions, 25 deletions (KIT-255's change only).
- **The replay marker list** (`workflowReplayRuleMarkers`, `R/utils_workflow_analysis_record.R`, read in the merged tree):
  `missing_value_comparison, enrichment_background, intensity_filter_samples, duplicate_resolution_samples,
  ruv_feature_levels, treat_moderation, ruv_control_grouping, ruv_k_bound, intensity_filter_missing,
  duplicate_resolution_rank, correlation_exclusion`. It ends as required; the doc comment (line 696) names KIT-217,
  KIT-258 and KIT-255.

## Declaration (written before the run is sent)

**Purpose.** KIT-255's merge check on the merged tree, head against staging: the correction pass's files on the tree
that holds KIT-190, KIT-217, KIT-258, KIT-259, KIT-264, KIT-265 and VAL-325's edit. The full suite is not repeated.

**Scope: 54 files plus R CMD check.**
- R2's 50 files (request `kit255-r2-c25e1335`), all present on the merged tree.
- `test-kit248-correlation-skip-record.R`, `test-kit217-intensity-filter-missing.R`, `test-ruv-control-anova-grouping.R`,
  `test-omics-module-globals.R` (the audit's closure check, section 6).

**Scope recount on the merged tree.** The ticket's grep
(`git grep -lE "filterSamplesByLipidCorrelationThreshold|filterSamplesByMetaboliteCorrelationThreshold|filter_correlation|_pearson_selection|pearson-selection|validateLipidDaAnalysisPreflight|validateMetabDesignDaPreflight|validateOmicsDesignAlignment|omicsIntegrationInputSource|writeOmicsIntegrationInput|runWorkflowOmicsReferenceCorrelation|runWorkflowOmicsManagedCorrelation" -- tests/testthat`,
helpers excluded) gives **31 test files** on the merged tree (the closure check's 31 at 176f9d3f; **30** on the staging
base, because `test-kit255-omics-correlation-assays.R` is the 31st). One helper still matches
(`helper-module-ci-metab-design.R`, not run as a file). All 31 are in the scope, `test-kit248-correlation-skip-record.R`
included. What changed since R2's recount: that file (it joined the grep on staging). The coordinator's four added files
are in addition to the grep.
**Not in scope, noted:** the heatmap grep (`generateLipidDAHeatmap|generateMetabDAHeatmap`) gives 17 test files on the merged
tree; seven are outside R2's list and outside this scope as instructed: `test-kit220-da-threshold-consumers.R`,
`test-metab-02aa-da-display-filter-characterization.R`, `test-metab-02ab-da-module-characterization.R`,
`test-metab-02c-da-results-output-characterization.R`, `test-metab-da-results-output-shared.R`,
`test-workflow-omics-application-reads.R`, `test-workflow-omics-figures.R`. (Some joined after R2; R2's own four heatmap files are in scope.) Reported
for the coordinator, not added.

**Run settings.** Copied from KIT-255's R3 request (`kit255-r3-23fe97d7`): `file_timeout_minutes` 45, no extra packages,
`blas` reference, `os` ubuntu-22.04, `check` true, `fetch_refs` empty, no scripts, `env`
`MULTISCHOLAR_TRX_MANIFEST_NOTE=1` and `MULTISCHOLAR_KIT105_W3_DUMP=1`. Changed from R3: **4 shards per side** (R3: 2)
for 54 files. Taken from R2, because R3 did not run these files: `file_timeouts` 120 minutes for
`test-prot-norm-core.R`, `test-gui-entry-lipidomics.R`, `test-gui-entry-metabolomics.R` (R2's values, unchanged).
Base = `7f17547f` (staging as branched), head = `1567f188`; one run.

**Ceiling.** 600 job-minutes = 40.0 core-hours (job-minutes x 4 / 60). Expected about 480 job-minutes (32 core-hours): R2's
50 files with 4 shards per side took 439.4 job-minutes; the four added files are small (the three non-omics-module files
about 10 to 20 s each in earlier trains; `test-kit217` and `test-ruv-control-anova-grouping` a few minutes). The scope
fits; if the first shard times show it will not, the run is still left to finish (an observation is not a failure).

### Predictions

Counts are passed/failed/errors/skipped read from the newest committed table that holds the file (head side of the newest
train table: `merge-train-20261126q-compare.tsv`, `merge-train-20261126p/compare.tsv`, `...m`, `...k` and earlier, or
`kit-255/r2` and `r3b`). The base count is the staging base's; where a file's newest table is older than later staging
changes, the prediction is the verdict (`same`) and a base count that differs is read from the run.

| File | Base | Head | Verdict |
|---|---|---|---|
| test-kit255-omics-correlation-assays.R | absent | 360/0/0/0 (R3b) | only-on-head |
| test-omic-norm-module.R | 639/0 (train 23) | 631/8: new failures exactly the eight digests, lipid `:671`, `:695`, `:712`, `:743`; metab `:896`, `:920`, `:937`, `:968` (VAL-303) | new-failure |
| test-workflow-omics-declarations.R | 43/25 (table ...m) | 43/25; only the selection step's `fixed` lines differ | message-changed |
| test-workflow-omics-binding-steps.R | 57/2 (train 23) | 57/2; only the selection step's parameters (3 against 2) differ | message-changed |
| test-metab-s4-qc-shared.R | 32/0 (train 23) | 35/0 | count-changed |
| test-kit236-report-sentences.R | 289/0 (train 24) | 293/0, all passed | count-changed |
| test-kit254-omics-correlation-skip.R | 379/0 (train 24) | 379/0 | same |
| test-omics-module-globals.R | **14/0** (newest committed table: `kit-213/r3`, `checkpoint-20261126`) | same as base | same |
| test-gui-entry-lipidomics.R | 19/16 (train 23) | same as base | same |
| test-gui-entry-metabolomics.R | 57/48 (train 23) | same as base | same |
| the other 43 files | as the newest table | same as base | same |
| R CMD check | 0 ERROR, 3 WARNING, 3 NOTE | identical findings, nothing only on head | identical |

Where this differs from the audit's: `test-omics-module-globals.R` is predicted at the base's committed 14/0, not 13/1
(no committed table holds 13/1; the verdict stays `same`). Other base counts agree with the audit's.
Failing on both sides and not KIT-255's, with identical lists expected: `test-prot-norm-core.R` (1547/307 at train 23),
`test-workflow-saved-da.R` (755/196), `test-omic-da-module.R` (1288/81), `test-workflow-operation-registry.R` (1107/2),
`test-workflow-kit-harness.R` (734/1), `test-kit127-omics-ruv-factors.R` (425/1), and the two GUI-entry files.
`test-metab-03c-norm-module-seams.R` runs no test on either side (0/0/0/0).

### Decision and attempts

- Every prediction holds: TRAIN GREEN for KIT-255.
- Any other outcome is a result: recorded with its files and counts, not repaired, not re-run, handed back.
- One run; a re-send only for an instrument defect with no result judged, with the coordinator's word. An HTTP 500/502 at
  dispatch with no run created: started once by hand with the README's command, recorded.

## What the coordinator resolves at landing

- The ledger: the row `VAL-<train25-1>` takes its number (next free after the rows landed before it; this train's base
  ends at VAL-343); its status cell "merged (<commit>), train 25" takes the merge commit (the placeholder
  `<commit>` is in the same row, `docs/methods-decision-ledger.md` line 2516 on this branch).
- VAL-283's status: "implemented by KIT-255 (VAL-<n>), merged by train 25 (run <id>)". VAL-299's dated note (the replay
  statement under the current rule also names samples set missing in one assay), as the audit's I05 worded it.
  VAL-330's note (C-R4) is already in the row, from `kit-255b`.
- NEWS: KIT-255's bullet in landing order.
- No other placeholder in `docs/qa/workflow-execution/kit/kit-255/` (searched for `<merge>`, `<commit>`, `<date>`, `<N>`).
- Candidates and records owed (closure check C.2 (7)): KIT-255.md, DISPATCH-STATE, C1 (the DA tab and tables mark not-tested
  comparisons), the KIT-242 hand-off, I10, RC2-06, I11, I12.

## Change after the run was sent (VAL-344)

The owner's decision VAL-344 (staging 31a4d684) fixed the RUV sentence of KIT-255's NEWS and `wording.md` item 8 as the
text quoted in its (OQ-A), and turned VAL-330's dated note into a pointer. Run 45726836882 was already out, so it was not
disturbed: one documents-only commit follows the declaration on this branch (and the same on `kit-255b`, aa184024).
**The run's head `1567f188` differs from this branch's tip by documents only**: `NEWS.md`,
`docs/methods-decision-ledger.md` (VAL-330's note; the KIT-255 row's "Information for the owner" clause),
`docs/qa/workflow-execution/kit/kit-255/outcome.md` and `wording.md`, plus this file. No `R/`, `tests/`, `inst/` or `man/`
file changed. The base stays `7f17547f` (staging has only gained VAL-344's ledger row since; the merge into staging at
landing will meet it in the ledger).

## Result

Run 45726836882 (request `kit255-train25-1567f188`, request commit 8f30e413 on `run-tests`), created by the dispatcher
without error (no hand dispatch); green (infrastructure). Head `1567f188` against base `7f17547f`; fetched once. My
`gh run watch` ended on a TLS handshake timeout from the API (not a run failure); the run was then read by single status
calls and the one fetch. Tables: `merge-train-20261126r-compare.tsv`, `-check-compare.txt`, `-tests.tsv`,
`-tests-base.tsv` (this directory).

**Job-minutes 464.6 = 30.98 core-hours** (job-minutes x 4 / 60) against the ceiling of 600 (40.0 core-hours) and the
expected about 480: tests head 59.3/63.1/55.6/40.4, tests base 53.6/61.0/53.8/52.0, check head 9.2, check base 15.8,
plan 0.3, results 0.5 (job start to end from the API).

| File | Predicted | Observed base | Observed head | Verdict |
|---|---|---|---|---|
| test-kit255-omics-correlation-assays.R | only-on-head, 360/0 | absent | 360/0/0/0 | only-on-head: HELD |
| test-omic-norm-module.R | new-failure, 631/8, the eight digests | 639/0 | 631/8 | new-failure: HELD. The eight failures are the (b) lipid and metab "skip" scenarios' correlation, final_plot, export and reset provenance digests (actual 15e689a4 x3, c0ed83bb; 9b3649ec x3, 62d62f5a: the predicted values); the stored-file line numbers :671 to :968 are not shown by the tables (the traces point at `expectOmicNormScenario`, test lines 1560 and 1565) |
| test-workflow-omics-declarations.R | message-changed, 43/25, selection step's `fixed` lines only | 43/25 | 43/25 | message-changed: HELD (the head side differs from the base side only by `fixed` is a list against NULL, for `lipid_pearson_selection` and `metab_pearson_selection`) |
| test-workflow-omics-binding-steps.R | message-changed, 57/2, parameters 3 against 2 | 57/2 | 57/2 | message-changed: HELD (only `steps[[10]]$parameters` length 3 against 2, lipidomics and metabolomics) |
| test-metab-s4-qc-shared.R | count-changed 32 to 35 | 32/0 | 35/0 | count-changed: HELD |
| test-kit236-report-sentences.R | count-changed 289 to 293, all passed | 289/0 | 293/0 | count-changed: HELD |
| test-kit254-omics-correlation-skip.R | same, 379/0 | 379/0 | 379/0 | same: HELD |
| test-omics-module-globals.R | same, **14/0** | 13/1 | 13/1 | same: verdict HELD; **the declared count was wrong** (see below) |
| test-gui-entry-lipidomics.R / -metabolomics.R | same at the base's counts | 19/16, 57/48 | 19/16, 57/48 | same: HELD |
| the other 43 files | same | same | same | same: HELD |
| R CMD check | identical | 0 E / 3 W / 3 N | 0 E / 3 W / 3 N | identical (6 findings on both; only on head: none; only on base: none): HELD |

Totals: 54 files, 48 `same`, 1 `only-on-head`, 1 `new-failure`, 2 `message-changed`, 2 `count-changed`; no verdict outside the
predicted set.

**Prediction error of the declaration (verdict unaffected).** `test-omics-module-globals.R` was declared at 14/0 (newest
committed table); it is 13/1 on both sides. Coordinator's note, quoted: "tests/testthat/test-omics-module-globals.R is 13
passed / 1 failed on BOTH sides of your comparison (staging holds KIT-261's `shown_idle <<- idle`; KIT-278's fix is not on
your base), verdict same." The verdict `same` held; KIT-255 adds no match to that file's patterns.

**Decision.** Every verdict prediction held. TRAIN GREEN for KIT-255. The count prediction for
`test-omics-module-globals.R` missed (13/1, not 14/0): an error of the declaration, not of the change.
Not covered: the seven heatmap-consumer test files named in the scope section; the full suite (the next checkpoint).
