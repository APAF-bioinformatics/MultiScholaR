# KIT-219 predictions (committed before each run; none is changed afterwards)

Budget (ticket 11): expected 111.7 core-hours, ceiling 485.33; core-hours = job-minutes x 4 / 60; running total in
`runs.tsv`. Rows G, R1, R2, R2b/C (conditional), R3, Correction (one pass), Instrument (up to three per item).

## Run G (G0 on the base): `g0.R` on the tests-only commit (R/ identical to base 0bdf72131)

Purpose: for the real kernel on both routes (g:Profiler with `gost` answered by a fixed result and no service call, GO
offline), list the files a run writes, the names the rule derives from the stored result, and the two contrast-name
sources. Ceiling 60 min (4.0 core-h), expected 2.0. Decides stop rule 2. Binding: the commit named in `runs.tsv`.

Prediction (by reading the kernel: `R/func_general_enrichment_orchestration.R`, `R/func_general_enrichment_engine_helpers.R`):
- run1, both routes: no error; `written_not_derived` and `derived_not_written` are both empty (every file written is a name
  the rule gives, and every name the rule gives is written, since both directions have a result with terms).
- g:Profiler route: `source_da_data_names` equals `source_contrast_names_argument` (both `groupB-groupA`, the contrasts
  table's raw name) and `source_manifest_results_names`; so the plots, named from the `contrast_names` argument, carry the
  sanitised name of the stored result's contrast.
- run2 (down list empty): the second run's derived names exclude the down files of the first run, which stay on disk
  (`first_run_files_not_in_second_run_names` = the three down files = `of_which_still_on_disk`).
- `package_current_files`: `<function absent>` (base).
A violation of the first or second line is stop rule 2 (report; the coordinator decides). An error in a route is an
instrument question first (mechanics of the script, no result), decided from the console log.

## Run R1: the new test file on the tests-only commit (no base comparison)

Purpose: each failing-first row fails at its stated reason; each G row passes. Ceiling 60 min (4.0), expected 2.0.
Decides stop rule 1 (a failing-first row that passes on base is dropped with its route recorded; if E1 to E4 all pass,
stop before the fix). One file, one shard, no check.

Base values (the first failing expectation of each row; the rest of the row also fails or is not reached):
| Row | Prediction on the tests-only commit |
|---|---|
| E1 | fails: the enrichment workbook has 3 result sheets (the planted table is in it), not 2; the planted files are in `Enrichment_Plots`; no record |
| E2 | fails: the workbook has 2 result sheets and `Enrichment_Plots` holds all six files; no record |
| E3 | fails: the workbook has 2 sheets; both tables are copied; no record |
| E4 | fails: `renderFn` is called once (render = 1), no refusal |
| E5 | fails: the workbook has 2 sheets; both tables are copied; no record |
| E6 | fails: `renderFn` is called once |
| E7a | fails: no line `enrichment_ran <- `, the two chunk headers lack `eval=enrichment_ran`, no "No enrichment analysis was run" line, the figure chunk holds `cat(glue::glue("Warning:` (all five expectations, both templates) |
| E7b | fails: the chunk `enrichment_presence` does not exist (`length(code) > 0` FALSE) in both templates |
| E8 | fails: with a file in `pathway_dir` and none in the publication folder `enrichment_ran` is TRUE in all four templates |
| G1, G2, G3 | pass (the unchanged path; G3 for lipidomics and metabolomics) |
| G5 | fails on base (the limpa template's STRING chunks read `current_project_paths$pathway_dir`): the template change's own row |
Instrument risks named in advance (VAL-54; not results): a helper that does not run on the runner (the owner fixture, the
module-CI helpers); decided from `failures.txt`.

## Run R1h (head check, declared 2026-11-27 before the request): the new test file and G0 on the fix commit

Purpose: all rows of the new file pass on the fix (A2's "all pass on the fix") before the large R2 is sent, and G0 is
run on head for G4's check with the package function (`package_current_files`) and for the g:Profiler route's (ii),
whose mock was repaired (the kernel retries a NULL answer; instrument defect, no result changed). One request: one tests
job (the new file) and one scripts job (`g0.R`). Ceiling 8.0 core-h (two jobs of 60 min), expected 3.3, drawn from the
total's slack (ticket total ceiling 485.33; used so far 6.1).
Predictions: every test of the file passes (E1 to E8, G1 to G3, G5). G0 on head: run1 both routes `package_equals_derived`
TRUE and `written_not_derived` empty; g:Profiler run2 derived excludes the three down files, which stay on disk
(`first_run_files_not_in_second_run_names` = `of_which_still_on_disk` = 3 names); GO run2 as on the base.
A failing row on head is a result of this change (the ticket's one correction pass); a row whose mechanics fail without a
result is an instrument defect (VAL-54).

## R2 scope and price (not yet requested)

The ticket's regression scope formula gives 49 files at the branch point (39 by the grep, helpers dropped, with the new
file, +4 GUI-entry lane files, +1 `test-prot-summary-report-templates.R`, +5 template-naming tests): above the ticket's 40,
so the coordinator re-prices R2 before it is sent (ticket 10). Not sent by this implementer.

## Run R2 (declared before the request; coordinator's rulings 1 to 3 of 2026-11-27)

Merge: `origin/v0.9-staging` a67c4161 merged into kit-219 once (conflict only in the ledger: both sides kept, staging's rows first).
Base commit `a70b352c5` = staging a67c4161 plus the tests-only commit (branch `kit-219-base`); head = the tip of kit-219 named in
the request. Scope by the ticket's formula on the merged tree: 50 files (`r2-predictions.tsv`: the grep's files, the four GUI-entry
lanes, the templates test, the template-naming tests, the e2e and go-propagation files; `test-kit266-limpa-combined-qc-figure.R`
enters by the grep since the merge, `test-kit219` and `test-kit216` run on both sides). `check: true` on both sides.
Plan: 4 shards per side from per-file seconds (train 35's full-suite tables, `run-tests-45762830948`, base side; kit216 585 s and
kit219 121 s from head records; the planned load per shard 45.3, 38.5, 38.2, 38.3 minutes, the list ordered so the heavy files
fall in different shards). Only `test-kit266-limpa-combined-qc-figure.R` is above 20 minutes (2,638 s = 44 min): explicit timeout
66 minutes (1.5 x). Planned job-minutes: 8 test jobs x (shard + 10 set-up) = 401, check 2 x 90 = 180, the G0 (iii)/(iv) scripts job
20: **601**. Ceiling **1.6 x 601 = 962 job-minutes = 64.1 core-hours** (not above 1,400). Env `MULTISCHOLAR_TRX_MANIFEST_NOTE=1`
(a transcriptomics file is in scope); no file of the diff matches `^(mod|func)_(lipid|metab|omic|omics)...R$`, so
`test-omics-module-globals.R` is not added.
Predictions per file: `r2-predictions.tsv` (verdict `same` for every file except `test-kit219-enrichment-files-current.R`: `fixed`).
Check: no finding only on head. Failures identical on both sides are "unchanged".
G0 (iii) and (iv) (`g0b.R`, scripts job, fix commit): prediction: no W1 fixture is a listed leaf (the lipid and trx fixtures hold no
enrichment table; the metab fixtures with DA files are refused by KIT-212 first; the in-memory proteomics fakes get NULL and the
proteomics owner is refused by KIT-212's unmade DA files; the `none` metab case has no files). (iv): every line is reviewed by
hand; a line asserting the contents of `Enrichment_Plots` or the workbook's sheets for a case where this change omits files is a
leaf (OD-3, stop). The GUI-entry lanes are decided by R2.

## Run R2b (the ticket's one correction pass; declared before the request)

R2 (run 45785931266, 50 files) showed two files moved by this change: `test-omics-summary-core.R` (25 to 34 failures) and
`test-omic-summary-module.R` (28 to 54): the W1 harness's recording `copyToResultsSummary` has fixed formals, and
`copyWorkflowSummaryResults` passed `enrichment_files = character()` to it for the reference-state lipid and metab owners even
when the enrichment folder held nothing to leave out ("unused argument", the Copy notice changed). Correction: the argument is
passed only when the folder holds a file outside the current list (otherwise the copy is the unchanged one: same files, same
workbook). G0b's mechanics (the working directory under `test_file`) repaired as instrument repair 1 of 3 for that item.
Rerun, base and head, of the files that can reach the Copy call (26: the 22 files that name the copy or the W1/summary harness,
`r2b-predictions.tsv`, plus the four GUI-entry lanes), `check: true`, and `g0b.R`. Plan: 3 shards per side (45.1, 40.0, 39.7
min), 6 x (shard + 10) = 310, check 2 x 90 = 180, g0b 20: 510 job-minutes; ceiling 1.6 x 510 = **816 job-minutes (54.4 core-h)**.
Predictions: every file `same` (the two moved files return to their base failure lists: 25 and 28 failures on both sides) except
`test-kit219-enrichment-files-current.R` `fixed`; check findings identical; G0b: no listed leaf, as declared for R2.
Files not rerun keep the R2 evidence (they do not reach the changed line): the other 24.

## Run R2c (infrastructure and instrument rerun; declared before the request)

In R2b the base side's shard 3 lost its job after 7 minutes (no result for 8 files: `test-kit212-da-current.R`,
`test-kit243-summary-current-outputs.R`, `test-kit257-norm-own-files.R`, `test-kit271-limpa-table-name.R`,
`test-kit259b-second-correction.R`, `test-kit207-summary-head-guard.R`, `test-kit219-enrichment-files-current.R`,
`test-gui-entry-lipidomics.R`), and `g0b.R` stopped on a call expression whose target is not a name (instrument repair 2 of 3 for
G0b; no expectation changed). Rerun of those 8 files, base and head, 2 shards per side (22.2 and 17.5 min), no check (R2b's
check ran on both sides and is identical), and `g0b.R`. Plan: 2 x (32.2 + 27.5) = 119 + scripts 15 = 134 job-minutes; ceiling
1.6 x 134 = **214 job-minutes (14.3 core-h)**. Order: kit212, kit243, kit257, kit271, kit259b, kit207, kit219, gui-entry-lipidomics.
Predictions: all `same` (the head counts of R2b: kit212 309, kit243 246, kit257 133, kit271 78, kit259b 90, kit207 292, gui-entry-
lipidomics 19/16) except kit219 `fixed`. G0b: as declared for R2.

## Run R2d (coordinator's ruling I03 of 2026-11-26; declared before the request)

One new test row, tests only, on unchanged code (`test-kit219-enrichment-files-current.R`): E9, the decided loop through the real
entry points: a stale file is in the publication folder beside the current files; `assertWorkflowSummaryEnrichmentCurrent` raises class
`workflow_summary_enrichment_not_current` with T2's text (by value) and the Report is refused with the notice "[error] Report
generation failed: <T2>", nothing rendered, the bundle recorder not called; Copy again leaves the stale file out (still on disk), the
folder holds the current files only, and the notice is T1 naming "1 enrichment file(s)" and the file; the Report then renders. E4 also
asserts the bundle recorder is not called (`report$bundle`, a counter added to `k219Report`; an added assertion, none removed).
Prediction on HEAD: E9 PASSES, E4 passes, every other row as in R1h. On the BASE (staging plus this test file): E9 fails: the
function `assertWorkflowSummaryEnrichmentCurrent` is absent (its expectations are skipped), the first Report renders (render 1, not 0),
the notice is absent, Copy again copies the stale file (the folder holds 7 files) and no T1; E4's new assertion passes on base (the
bundle recorder is not reached only if the Report renders: on base E4 renders, so `report$bundle` is 1 and the assertion fails). If
E9 fails on head it is a RESULT: stop and hand back, no repair. Run: the new file, base against head, 1 shard per side, no check.
Plan 2 x (3 + 10) = 26 job-minutes; ceiling **6.0 core-hours (90 job-minutes)**. Used before it: 66.32.
