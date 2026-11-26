# KIT-271 step 0 (before any test or code)

Base: origin/v0.9-staging 894b9cc821a70b03baab1e1fa47b3d47440c7ed5 (fetched 2026-11-26), branch `kit-271`, worktree
`MultiScholaR-kit-271`. Authority: ticket r1.1 and AUDIT-KIT-271-plan-20261126 "Re-check 1", which governs where it
differs. Owner's decision: ledger VAL-317 (Q1 = N1; Copy's other deliveries are KIT-273; the stored values wait).

## Bounded work item

Outcome: on DIA with limpa the correlation filter's table is stored, exported and delivered as
`normalised_results_cln_with_replicates.tsv` / `normalised_results.tsv`; an earlier build's project (only the `ruv_`
file on disk) keeps working. Files: `R/utils_workflow_dia_norm.R` (a new internal function; the pinned resolver is not
edited), `R/utils_workflow_dia_products.R`, `R/mod_prot_norm_observer_helpers.R`, `R/utils_workflow_dia_sessions.R`
(`collectDIANativeSessionMetadata` call and `collectDIANativeGuiFiles` fallback), `R/func_general_filemgmt_results_summary_helpers.R`
(one Copy branch), `tests/testthat/test-kit271-limpa-table-name.R`, NEWS, ledger, these records. No stored value edited.
Exit: the ticket's section 8 on R1 and R2 (A0 to A6).

## Cited lines re-read on the base (all as the ticket reads them)

- Resolver `R/utils_workflow_dia_norm.R:742-768`: NULL results fall through to `ruv_normalised_...` (the "applied" branch);
  definition digest pinned at `tests/testthat/test-prot-norm-core.R:31` (read, unchanged by this ticket).
- Limpa pass-through `utils_workflow_dia_norm.R:702-726`: `norm_data$ruv_optimization_result <- NULL` at `:716`.
- Callers: `utils_workflow_dia_products.R:1976-1979` (writer), `mod_prot_norm_observer_helpers.R:226-228`,
  `utils_workflow_dia_sessions.R:1917-1918` (`results_file`); the fourth (`saveProtNormCorrelationResults`) is untouched.
- `collectDIANativeGuiFiles` `utils_workflow_dia_sessions.R:1675-1681` (the only reader of `results_file`; `git grep results_file -- R`
  confirms); Copy `func_general_filemgmt_results_summary_helpers.R:373-405` (one `getDIAActiveManager()` use, `:403`).
- W1 captures: `test-prot-norm-core.R:862, 897, 933, 992, 1027, 1063` (`files`), `:875, :1005` (resolver `calls = 1L`),
  counters at `:356-359, :447-459`; `test-workflow-saved-da.R:1136`; `tests/testdata/gui-entry/prot_dia_limpa.json`
  13 lines (`:2760 ... :5667`, ten steps).
- Scope grep at the branch point: 29 files (list as the ticket), plus `test-gui-entry-proteomics.R`,
  `test-prot-norm-worker-figures.R` and the new file = 32 files.

## Re-check 1 repairs applied to the plan

- R01 (major). The new function `resolveProtNormCorrelationTableNames(manager, normRuvOptimizationResult,
  workflowRuvOptimizationResult, messagePrefix, messageFn)`: on `workflow_type == "DIA_LIMPA"` it calls the pinned resolver
  once as `resolveProtNormCorrelationResultFilenames(list(ruv_skipped = TRUE), NULL, messagePrefix, messageFn = function(...) invisible(NULL))`,
  returns that list and writes its own console line; elsewhere it returns the resolver's own value (one resolver call on every
  route). So the W1 `calls` (`:875, :1005`) do not move and the resolver stays the only holder of the names. The resolver body
  is not edited. T6 gets DIA_LIMPA rows (NULL, applied, skipped), run only where the function exists.
- R02. Copy's fallback and the export fallback take the new name first, then the old. Test T8b (both present, head only).
- R03. Step 0's "already fails on runners" comes from ONE existing run's artifact fetched with `fetch.sh` (run 45650328815,
  train 20261126c, the newest run holding all three files; no new run). The W1 dump variable is
  `MULTISCHOLAR_KIT111_W1_DUMP=/dev/stdout` (a file path), `MULTISCHOLAR_KIT112_W1_DUMP=1`.
- R04. T6's new-function rows are in a `test_that` of their own with `skip_if_not(<function exists>)`: skipped on base, run on head.
- R06. The export fallback sits in `collectDIANativeGuiFiles` (covers the new export, the Summary's archive refresh and metadata
  from builds before KIT-206): on a DIA-with-limpa owner it lists the new table if it is on disk, else the old one if it is,
  else the name the metadata gives. No edit beside `:1917` other than the new function's call.
- R05 (stale ticket text, for the coordinator; the ticket is not edited here): the description's "the resolver learns the
  route, four callers pass it" (superseded by r1.1); the dependency note's "Merge order: any" (now: before KIT-266, any with
  KIT-259); section 2.4's four rows (Copy "none to code" -> fallback 2; metadata "no change needed" -> fallback 1; stored
  references "five file listings in one JSON" -> eight places; tests "the new file; none else changes" -> two test files'
  stored values).

## A0: the eight stored places and the two `calls`, from run 45650328815 (train 20261126c, head 4baea49f, base 828dd105)

| Place | Already fails on runners at base? | Evidence |
|---|---|---|
| `test-prot-norm-core.R:862` (`limpa_covariates_correlation files`) | yes, in process and worker | failure text: actual `926838a9...` vs expected `34e64dc1...` |
| `:897`, `:933` (`limpa_covariates_ruv`, `_ruv_repeat` files) | yes, both paths | same list |
| `:992`, `:1027`, `:1063` (`limpa_correction_*` files) | yes, both paths | same list |
| `test-workflow-saved-da.R:1136` (`limpa_one_contrast files`) | yes, in process and worker | actual `8fc2b731...` vs expected `2ce30818...` |
| `prot_dia_limpa.json` (13 lines, ten steps) | NO: the lane passes. The 25 failures of `test-gui-entry-proteomics.R` (168 passed) are in the other four lanes (prot_dia 7, prot_lfq_fragpipe 6, prot_lfq 6, prot_tmt 6); no failure names prot_dia_limpa | failures.txt of the run; the same 168/25 on every train 20261126a to e |
| W1 `calls` of the resolver (`:875`, `:1005`) | no failure on `calls` at base | they stay 1 by R01 |

The ticket's section 6 expected the gui-entry file to be `message-changed` or `count-changed`; step 0 shows the limpa lane
passes today, so its ten changed steps will be NEW failures.

## Predictions

R1 (run of the tests-only commit, new file only, base code): per row, on the base.
- T1 (in process and worker): FAIL at value (present `ruv_normalised_...tsv`, expected `normalised_...tsv`); the setup errors row passes.
- T2 (covariates, correction): the name row FAILS (after = `ruv_` name); `after == before`, the content digest and the
  adjusted matrix rows pass; state rows pass.
- T3 (in process, worker): `after == before` passes, name row FAILS.
- T4 (both paths): `normalised_results.tsv %in% tables` FAILS, `RUV_normalised_results.tsv` absent FAILS, the digest row
  passes (the delivered RUV file holds the same values: A5's base digest is printed).
- T5: `results_file` row, destination-new row and destination-old-absent row FAIL at value; the error row passes.
- T5b (both modes): `results_file`, destination rows FAIL at value; error and archive rows pass (export succeeds on base).
- T6 a passes. T6 b is SKIPPED (function absent). T7 passes (both routes).
- T8 (export): passes (guard). T8 (Copy): both rows FAIL (`normalised_results.tsv` absent, `RUV_normalised_results.tsv` present).
- T8b: SKIPPED.
- Stop rule 1 would fire if T1's base shows `normalised_` or the reverse listing is not the base.
If a row errors for its mechanics (the export fixture's native flag, Copy under a scope), that is an instrument defect
(brief rule 3, at most three fixes) unless it is a failure at a value.

R2 (base = tests-only commit, head = the fix), by file, 32 files, check on both sides.
- `test-kit271-limpa-table-name.R`: `fixed` (base fails the rows above, head passes every row; T6 b and T8b run on head).
- `test-prot-norm-core.R`: `message-changed`; same counts (the six W1 (a) limpa `files` and their W2 worker repeats already
  fail on base; only their actual digest changes); the `calls` fields, `:31`'s definition-digest block and the pass-through
  block `same`.
- `test-workflow-saved-da.R`: `message-changed`, limited to `limpa_one_contrast files` (W1 (a) and W2), counts the same.
- `test-gui-entry-proteomics.R`: `new-failure`: lane prot_dia_limpa, the ten steps with the changed keys (168/25 becomes
  158/35); no other lane changes.
- `test-kit244-limpa-pass-through.R`, `test-kit257-norm-own-files.R`, `test-prot-norm-worker-figures.R`: `same`.
- Every other file (the remaining 26 of the scope list): `same`.
- R CMD check: identical (the new function is `@noRd`; no Rd or NAMESPACE edit).
- Declared stored values, none edited: the six W1 limpa `files` digests (message-changed each, W1 and W2), `limpa_one_contrast`
  `files` (message-changed), the 13 gui-entry lines (new failure on ten steps). Dumps: `MULTISCHOLAR_KIT111_W1_DUMP=/dev/stdout`,
  `MULTISCHOLAR_KIT112_W1_DUMP=1` so the head values are in the job logs.
Any other verdict stops the work (stop rule 2): a missed prediction is reported, not repaired.

## Budget (core-hours = job-minutes x 4 / 60)

Expected 30.93, ceiling 103.0 (ticket section 10): R1 20 min / 45 ceiling (1 job); R2 444 expected / 660 ceiling
(3 shards per side + 2 checks); one correction pass 705 ceiling; instrument allowance 135 ceiling. Core-hours used are
reported in outcome.md against these.
