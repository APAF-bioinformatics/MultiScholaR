# KIT-269 outcome (in progress)

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

## R2 predictions (unchanged from above; base = 50077f3a, head = 220b1e2e)

Base is the final tests-only commit, so the new file is on both sides: verdict `fixed` for it (base 15/5/4 or equal
by value; head 24 passed expected, 0 failed). Every other file `same`. Check identical.
