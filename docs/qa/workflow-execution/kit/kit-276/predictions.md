# KIT-276 predictions (committed before R1 and before R2)

R1 = the tests-only commit against itself on the runner (one job, `base_commit` none): each new row at the stated value.
"By value" = an expectation fails with a value, no error before it. Ruling of 2026-11-26: T1 asserts the record's steps,
T2 the report-level text; T18 and T17 assert both mismatch readings.

| Row | Base (tests-only commit) | Head (fix) |
|---|---|---|
| T1 (in process, worker) | history ends `cf, ruv_covariates`: FAILS by value (history, record steps); the DA guard (`w_terms` empty, `ruv_mode` none) passes | passes |
| T2 | covariates Methods paragraph and Results text printed: FAILS by value | passes |
| T3 | history ends `ruv_covariates, ruv_corrected`; two RUV steps in the record: FAILS by value | passes |
| T4 | head stays `ruv_covariates`; old text: FAILS by value | passes |
| T5 | head stays `ruv_corrected`, DA mode "correction": FAILS by value | passes |
| T6 | precondition passes (DA current, enrichment product present); after the Skip status "current": FAILS by value | passes |
| T7 | no refusal (old Skip text): FAILS by value | passes |
| T7b | nothing dispatched (`list()`): FAILS by value | passes |
| T8 | summary assigned never, holders keep the run's values: FAILS by value; the first-correlation guard and the builds guard pass | passes |
| T8m | `finish` rejects the new argument: the captured error string FAILS the `NA` expectation by value | passes |
| T19 | recorded call `list(final = "FINAL")` for the RUV head: FAILS by value; the no-RUV head call passes | passes |
| T8b | `ruv_mode` "covariates", applied TRUE, `ruv_k` set: FAILS by value | passes |
| T8c (three call sites) | file stays, study parameters print the section: FAILS by value (each of three) | passes |
| T9 | limpa summary stays: FAILS by value; the other three routes pass (untouched) | passes |
| T12, T13 | pass | pass |
| T14 | first part passes; the DA and the correlation leave `ruv_covariates`: FAILS by value | passes |
| T15 | `ruv_covariates` still listed, selection allowed: FAILS | passes |
| T16 | passes | passes |
| T17 | replay refused (`workflow_replay_invalid` captured as the value): FAILS by value | replay states equal the history; mismatches over the origin map exactly `ruv_covariates/missing`; over the history's states none |
| T18 | history `cf, rc, rk`; replay refused (class captured): FAILS by value | replay states equal the history; (i) listed none; (ii) exactly `ruv_corrected/missing` |
| note 1.5 | `ruv_covariates_cancor.png` stays: FAILS by value | passes |
| `test-workflow-integration-input.R:703` | `ruv_corrected` is still in the history: FAILS by value | passes |

Stop rules: T1, T4 or T18 passing on the base; a guard (T12, T13, T16) failing on the head; any other stored value moving.
R2 (base = the tests-only commit, head = the fix; the R2 grep list of 97 files at the request; R CMD check both sides):
the new file only-on-head... passes completely on head; `test-prot-norm-core.R`, `test-kit244-limpa-pass-through.R`,
`test-kit206-...`, `test-kit248-...`, `test-kit257-...`, `test-kit129-run-path.R`, `test-dia-limpa-ruv-rules.R`,
`test-workflow-analysis-record.R`, `test-module-ci-prot-dia-limpa.R`, every other file: `same`;
`test-workflow-integration-input.R`: `fixed` (the declared edit); R CMD check identical.
Instrument risks recorded in advance (repaired only as VAL-54 mechanics): the mocked non-function binding in T19, the
enrichment recipe in T6, the replay harness paths in T17 and T18.

## Addendum before R3 (instrument repair 2 of 3, declared before sending)

R2 (45784213557) showed T2's last expectation failing on the head because the fixture's DPC state has no step record, so the
template did not take the limpa route (its R1 base failure count of 1 already showed two of T2's three expectations passing
on the base for the same reason). The repair gives the report chunks the DPC step as test-kit227's E7 does; no expectation
or text changes. Predicted with it: T2 on the base (tests-only e6f83d43) fails all three expectations by value (covariates
Methods paragraph, covariates Results caption, the two-stage text absent); on the head (fix 79e77621) it passes. R3 runs the
new file only, base against head.

## Addendum: the further bounded correction (VAL-389), committed before the small run

Base = 7649808b (the branch with row T8d and without the one line); head = the branch with
`metadata["ruv_optimization_result"] <- list(NULL)` in the `run$ruv_mode`-is-NULL branch of `refreshDIANativeSummaryArchive`.
- T8d on the base: the sweep is kept: `metadata$ruv_optimization_result` is not NULL, `ruv_optimization_results.RDS` is back and
  the study parameters print the section: FAILS by value (the fields `ruv_mode` "none", `ruv_applied` FALSE and `ruv_k` NULL pass).
  On the head: passes.
- The precondition rows of T8d (the export wrote the file; the Skip removed it) pass on both; `protSummaryRefreshArchive`
  returns NULL (refreshed) on both. If it returns the stale text instead (the KIT-206 guard comparing the head with the
  export's RUV reference), T8d fails on both sides at that expectation: that is a miss, and I stop and hand back.
- Every other row of the new file: as in R3 (195/0 on the head). `test-kit216-prot-session-after-enrichment.R`: `same`
  (its refresh rows take the `run$ruv_mode` branch).
