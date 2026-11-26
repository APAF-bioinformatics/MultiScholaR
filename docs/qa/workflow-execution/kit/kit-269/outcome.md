# KIT-269 outcome

Ticket r1.1, base `v0.9-staging` 4d588ba3ed5c6fee76a5a1450c028b409ba419ba (fetched 2026-11-26), branch `kit-269`.

## Bounded work item

Outcome: a vector-valued setting is listed in `study_parameters.txt` on the path without an analysis record, every
value comma-joined and in full, and a listing that still fails says so in the file. Files: `formatConfigList`
(`R/func_general_filemgmt_config_io_helpers.R`), the writer's handler
(`R/func_general_filemgmt_study_parameter_writer.R`), the new test file, NEWS, ledger, these records. Exit: ticket
section 8.

## Budget declared before the first run

Core-hours = job-minutes x 4 / 60. Expected 28.57 core-hours in total (R1 1.33, R2 27.23), ceiling 99.0
(ticket section 10). R1: one shard, ceiling 45 minutes (3.00 core-hours), expected 20 minutes. R2: 3 shards per side
plus check on both, ceiling 630 minutes (42.00 core-hours), expected 408.5 minutes (27.23). One correction pass for a
result; instrument defects under the brief's four conditions, at most three.

## Predictions, committed before R1 (test file at the tests-only commit; no R/ change)

R1 binds the tests-only commit, no base (`base_commit` none), the one new file. It decides stop rule 1 (the base
block empty with no `Analysis Record:` line; T1 at the `map_chr` error).

| Test | Base prediction |
|---|---|
| T1 | error "Result must be length 1, not 2" (1 error) |
| T2 | error, same message (1 error) |
| T3 | 4 expectations pass (scalar, pair, empty list, skipped name); the seven-element row fails by value (`...` after five) |
| T4 | error, same message (1 error) |
| T5 | `Analysis Record:` absent (pass), one header, `Yaxis Limit: -3, 3` absent from the block (1 fail); the printed first six lines after the dashes are the block that follows the Additional block, not settings; no `Error formatting config` line (pass) |
| T6 | console line "WORKFLOW ARGS: formatConfigList failed: boom" present (pass), file returned (pass); for each of the three faults the line after the dashes is not the draft failure line (3 fails) |
| T7 | error, same message (1 error) |

File totals predicted: 15 passed, 5 failed, 4 errors, 0 skipped.

R2 (head against the tests-only commit as base, 46 files: the 44 grep files, `test-gui-entry-proteomics.R`, the new
file; check on both): every file `same` except the new file (`only-on-head`, all pass: 20 expectations); in particular
`test-omic-summary-module.R` same with identical failures, `test-general-filemgmt-config-roundtrip-contracts.R`,
`test-kit235-lipid-metab-report-parser.R` and `test-general-filemgmt-config-direct-shared.R` same; R CMD check
identical on base and head. Note: with the tests-only commit as base, the new file exists on base, so it is expected
to be `fixed` (not `only-on-head`) with the base failures above.

## Not checked (so far)

Everything: no R was run; no run has been sent.

## R1 (run 45707147607, tests-only commit 79680e76, no base; 11.3 minutes wall, one job)

File result 15 passed, 5 failed, 4 errors: exactly the prediction. T1, T2, T4: purrr error, "Result must be length 1, not
2" (T1, T4 at index 1 and 17 as predicted; T2 "not 3"). T3: only the seven-element row failed ("Alpha: v1, v2, v3, v4,
v5, ..."). T5: no `Analysis Record:` line; the file held 116 lines ending at the Additional dashes (the printed
"first six lines" showed NA, see the helper note); `Yaxis Limit: -3, 3` absent (1 failure). T6: console line present
(pass), file returned (pass), three failures for the line after the dashes. T7: error. Stop rule 1 not triggered.

Instrument defect, recorded (brief condition 3): the helper `k269After` indexed `lines[(header + 2):length(lines)]`, which
on an empty block returns `NA` and the dashes (a descending sequence). It changed no expected value; it was replaced by
`lines[-seq_len(header + 1L)]` (commit 50077f3a). The base values of T5 and T6 are unchanged in kind (they fail by
value); R2's base run re-shows them with the final file. No R1 re-send.

Plan re-check 1 (coordinator's message): R03 changed the failure reason to the condition's message with ANSI stripped and
white space collapsed (`conditionMessage(e)` already holds an rlang condition's parents), so T6's parent case asserts
the cause appears exactly once (commit 8f5c30e8, tests only). R1 ran before this edit; its expectation there was the
parent-only line, which failed on base by value, as the changed assertion does. Ticket text not edited; for the
coordinator: R01 (the dependency note's premise about KIT-263's test file: it runs the record path) and R02
(section 2 `:45`, stale after F5) stand as stale ticket sentences.

## R2 (run 45709024838; head 220b1e2e, base 50077f3a = the final tests-only commit; 46 files; check on both)

Observed (`compare.tsv`, `check-compare.txt`): 45 files `same`, 1 `fixed` (the new file), 0 other verdicts. The new file:
base 17 passed / 6 failed / 4 errors (T1, T2, T4, T7 error; T3 seven-element row, T5 `Yaxis Limit`, T6 three line
values fail; the same failures as R1 with the final file), head 33 passed / 0 / 0. Pre-existing failures identical on
both sides: `test-omic-summary-module.R` 1675/28 (the predicted 28), `test-gui-entry-proteomics.R` 168/25,
`test-omics-summary-core.R` 676/25, `test-integration-input-transcriptomics.R` 67/3, `test-module-ci-trx-summary.R`
83/2, `test-kit224-ruv-feature-levels.R` 65/1. Named predictions: `test-kit235-lipid-metab-report-parser.R` 54/0 same;
`test-general-filemgmt-config-direct-shared.R` 33/0 same; `test-general-filemgmt-config-roundtrip-contracts.R` 24/0
same; no stored digest moved. R CMD check: base and head each 0 ERRORs, 3 WARNINGs, 3 NOTEs, identical findings
(6 on both, none only on either side). Prediction check: the committed head count was 20 (90a1f289), then 24
(56872624); counted per call it is 31; observed 33 because `expect_match` evaluates its object twice, so T6's parent
case runs `line()` (two expectations) twice. The base was predicted 15/5/4; observed 17/6/4 (8f5c30e8's count row adds
one failure; the double evaluation adds two passes). Every verdict prediction held.

## Acceptance criteria

| # | Criterion | Status | Evidence |
|---|---|---|---|
| A1 | tests-only commit fails as stated; base block empty; no `Analysis Record:` line | met | R1 45707147607 (15/5/4 as predicted) and R2 base (17/6/4 with the final file); T5's base file ends at the dashes (116 lines) |
| A2 | all rows pass on the fix | met | R2 head 33/0/0 |
| A3 | every file `same` except the new file; check identical | met | R2: 45 same, 1 fixed; check identical (0/3/3 both sides) |
| A4 | printed form (a+); failure wording draft | met, wording draft | T1-T4 on head; `wording.md` |
| A5 | NEWS, ledger row, outcome with change report | met | `NEWS.md`, `docs/methods-decision-ledger.md` (`VAL-<merge>`), this file |

## Core-hours

R1 about 12 job-minutes (one tests job plus plan and results); R2 301.8 job-minutes summed from the job start and end times
(3+3 test jobs, 2 check jobs, plan, results). Total about 314 minutes = about 20.9 core-hours against expected 28.57
and ceiling 99.0. No correction pass, no further instrument run.

## Deviations

- Instrument defect in R1 (the `k269After` helper; section "R1"): fixed in commit 50077f3a without changing an expected value.
- R03 of the plan re-check changed the failure reason after R1 was sent (T6's parent case; commit 8f5c30e8). R1 was
  not re-sent; the base failure is by value either way and R2's base shows the final file.
- The R2 base is the final tests-only commit (50077f3a), so the new file's verdict is `fixed` (the ticket predicted
  "only on head" for a base that lacks the file).
- Stale ticket sentences for the coordinator (not edited): R01 (dependency note on KIT-263's test file: it runs the
  record path) and R02 (section 2 `:45`, replaced by F5).
- `cli::` is used in the handler as in `R/utils_workflow_summary.R:231`; `cli` is not listed in `DESCRIPTION` Imports
  or Suggests (pre-existing for about 1,000 uses; check identical). Unrelated; left alone.

## Not checked

- No R on the laptop; everything above is from the runner.
- Whether the lipidomics and metabolomics report templates read the newly printed block (only `test-kit235` cases M and L2, same).
- Windows. Real data.
- Merge with the partner branches by `git merge-tree` only: kit-269 against `origin/merge-train-20261126k` merges
  cleanly in the two R files; NEWS.md and docs/methods-decision-ledger.md conflict textually (adjacent insertions at the
  same anchor), for the coordinator at merge. KIT-213 and KIT-217 were not tried.

## For the change report

What was wrong and who was affected. The "Additional Configuration Parameters" block of `study_parameters.txt` was empty,
with no message, whenever the file was written without an analysis record, because a setting with several numbers could
not be listed: the shipped `[plotRle] yaxis_limit = -3,3` in proteomics, the RUV control vector `ctrl` in lipidomics
and metabolomics. The exported `createStudyParametersFile()` stopped on the shipped configuration ("Result must be
length 1, not 2"). In `main` (the lipidomics and metabolomics files after RUV too, by reading: `main:R/mod_lipid_summary.R:232`, `main:R/mod_metab_summary.R:232`, `ctrl` stored at `main:R/func_lipid_s4_objects.R:2552`) and the published release 0.3.6.3 (`v0.3.6.3:R/func_general_filemgmt.R:1155`, `:4047`; tags before it were not read), every proteomics project that uses the shipped configuration
wrote the file with an empty block (main: `main:R/func_general_filemgmt.R:1755` and `main:R/func_general_helpers.R:1072`
for `formatConfigList`, the call at `:5388` inside a `tryCatch` whose handler (`:5399`) loses its assignment; the
conversion of `yaxis_limit` at `main:R/func_general_filemgmt.R:1250-1252`). In 0.9 the application takes the record path
and is affected only when the record cannot be built; scripted callers are affected. Text settings of six or more values
were cut after five with "..." and are now printed in full.

Code references. Base 4d588ba3: `R/func_general_filemgmt_config_io_helpers.R:483-488` (five-value cut, character vectors
only; plain `as.character` for others), `R/func_general_filemgmt_study_parameter_writer.R:973-976` (handler assigned with
`<-`). Fix 220b1e2e (`R/func_general_filemgmt_config_io_helpers.R:483-485`, `R/func_general_filemgmt_study_parameter_writer.R:973-979`): the plain-value branch joins every value of length above one with ", " in full; the handler
assigns with `<<-` and writes "The configuration settings could not be listed here: <message>" (draft wording).

Results: unchanged (display only). Evidence: runs 45707147607 (R1) and 45709024838 (R2) above; base block observed
empty in R1 and R2's base (T5, 116 lines ending at the dashes).

## Correction pass (coordinator's ruling 2026-11-26 on audit I01, option (a)): named vectors keep their names

On the path without an analysis record a named vector is printed `name=value`, as the record path prints it (the RUV
control vector `ctrl` is stored named per feature; without names the controls cannot be identified from
`study_parameters.txt`). One change in `formatConfigList`'s vector branch; unnamed vectors print as before. Budget: at
most 4.0 core-hours (60 job-minutes) for both runs together; KIT-269 had used 20.86 of its 99.0 ceiling.

### Predictions, committed and pushed before run 1

Test added (T2): `formatConfigList(list(ctrl = list(assay_pos = c(Feat_A = TRUE, Feat_B = FALSE, Feat_C = TRUE))))`
must equal `c("Ctrl:", "  Assay Pos: Feat_A=TRUE, Feat_B=FALSE, Feat_C=TRUE")`.

- Run 1: the test file alone, no base, on the commit that adds the row, with `R/` as at 9e6c3fc3 (the current head's
  code). Predicted 33 passed, 1 failed, 0 errors, 0 skipped: the new row fails by value (the names are dropped:
  `Assay Pos: TRUE, FALSE, TRUE`), every other expectation as in R2's head (33 passed). Ceiling 45 job-minutes
  (expected about 11).
- Run 2 (after the fix): base 50077f3a (the base R2 used) against the fix commit, no check, four files: the KIT-269
  file (base 17 passed / 6 failed / 4 errors, head 34 / 0 / 0, `fixed`), `test-kit235-lipid-metab-report-parser.R`
  (54/0, `same`), `test-general-filemgmt-config-direct-shared.R` (33/0, `same`),
  `test-general-filemgmt-config-roundtrip-contracts.R` (24/0, `same`). Ceiling 45 job-minutes (expected about 30).
  `test-kit213-study-parameters-labels.R` is not on this branch and is left to the train.
- Decision: a missed prediction is a result: recorded, no repair, stop and hand back.

Change report, correction pass: the fix commit (below) adds, in the same vector branch, `name=value` for a named vector
(`R/func_general_filemgmt_config_io_helpers.R:483-487`), so the RUV control vector in proteomics, lipidomics and
metabolomics prints each feature id with its flag, as the record path does. Display only.

### Run 1 (45726785537, tests-only commit edf5d314, no base): prediction held

33 passed, 1 failed, 0 errors, 0 skipped (`r-corr1-tests.tsv`): the new T2 row failed by value, actual `Assay Pos: TRUE,
FALSE, TRUE`, expected `Assay Pos: Feat_A=TRUE, Feat_B=FALSE, Feat_C=TRUE`. 7.5 job-minutes (plan 0.28, tests 6.72,
results 0.53) = 0.50 core-hours.
