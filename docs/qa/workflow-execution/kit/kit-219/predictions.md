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
