# KIT-276 outcome (branch kit-276; implementation, for the independent audit)

State: done and ready for its implementation audit. Branch tip is the last commit of `kit-276` (fix 79e77621; tests-only
base for the comparisons: branch `kit-276-tests` e6f83d43 = e77f096c + two mechanics repairs of the new file, no change
under `R/`). Staging branch point 38443af8. Nothing merged. Running total 83.3 of 172.00 core-hours (restated expected 84.2).

## What changed (by file)

- `R/utils_workflow_dia_state.R`: `workflowHistoryAfterCommit()` (one table, `DIA_LIMPA` only) and its call in
  `commitWorkflowReferenceState` (the committed state moves to the end; a `correlation_filtered` commit removes the RUV
  states; a commit of one RUV mode removes the other; the step provenance is cut through `truncateDIAStepProvenance`).
- `R/mod_prot_norm_workflow_helpers.R`: the Skip of a RUV head (`runProtNormRuvSkipAtRuvHead`: busy refusal, select
  `correlation_filtered`, holders, leftover file, trigger, texts, Figure 1 in process or dispatched through
  `dispatchProtNormLimpaComposite`), `clearProtNormRuvHolders`, `protNormRuvSkipTexts`, `completeProtNormCorrelationAndCleanUp`
  (flag passed only when TRUE, then the leftover file), the three call sites, the dispatcher's hand-over to `finish`.
- `R/mod_prot_norm_correlation_helpers.R`: `ruvRemoved` argument of `completeProtNormCorrelationWorkflow`.
- `R/mod_prot_norm_session_helpers.R`: `runProtNormResetWorkflow` clears `ruv_covariates_summary` on the limpa route.
- `R/utils_workflow_dia_sessions.R`: `diaLimpaNoRuvMode`; archive `ruv_mode` "none", `ruv_applied` FALSE, `ruv_k` from the
  cleared holder (export), same rule in the refresh.
- `R/utils_workflow_dia_norm.R`: `removeDiaLimpaRuvMarkers` also removes `ruv_covariates_cancor.png` (note 1.5);
  `removeProtNormLeftoverRuvResults`.
- `R/utils_workflow_analysis_record.R`: the replay's two `protein_ruv_covariates` branches read `correlation_filtered`.
- Tests: `tests/testthat/test-kit276-ruv-superseded.R` (new, 21 rows); `test-workflow-integration-input.R:703` declared edit.
- `NEWS.md` (one entry; two sentences removed), `docs/methods-decision-ledger.md` (row `VAL-<merge>`), this directory.

## Runs (runs.tsv; tables committed under r1, r1b, r2, r3)

| Run | What | Result against the prediction |
|---|---|---|
| R1 45776354953 | tests-only e77f096c | 66 predicted failures at the stated values, all by value; T8b and T14 errored on test mechanics (instrument defect 1); `test-workflow-integration-input.R:703` failed as predicted; guards T12, T13, T16 passed |
| R1b 45780547662 | instrument repair 1 (T8b marks the owner a native session; T14's fixtures live in the test's frame) | T8b and T14 failed by value as predicted; 0 errors |
| R2 45784213557 | f847c433 vs fix 2f76cca9, 97 files, R CMD check both sides | 95 same, 2 fixed, 0 new-failure, 0 message-changed; check identical (0 ERROR, 3 WARNING, 3 NOTE both sides); `test-prot-norm-core.R` same 1555/307 (the ticket's 1547 was the older record; base and head agree); one expectation of the new file failed on the head: T2's last (below) |
| R3 45799233856 | instrument repair 2 (T2 gives the report chunks the DPC step as test-kit227's E7 does), new file, e6f83d43 vs 72ba4edc | base 123 passed / 72 failed (T2 fails all three expectations); head 195 passed / 0 failed |

Instrument repairs used: 2 of 3 (both mechanics only, declared before sending, no expectation or text changed; the R2
finding for T2 is recorded here as it was: a T2 expectation failed on the head because the fixture's DPC state has no step
record, so the report took the non-limpa route; the repaired row then fails on the base and passes on the head). Correction
pass: not used. Cost: R1 1.89, R1b 1.49, R2 77.0, R3 2.8 = 83.2 core-hours.

## Acceptance criteria

| # | Evidence |
|---|---|
| A1 | R1/R1b: every failing-first row failed by value at its stated value; guards passed. Met |
| A2 | R2/R3 head: all 195 rows pass in process and on the worker lane; the DA after each limpa sequence carries the state the record, tab and archive state. Met |
| A3 | R2: every file same except the new file and the declared edit; `test-prot-norm-core.R` same; check identical. Met |
| A4 | Skip is S1; texts are drafts (wording.md below). Met |
| A5 | T14, T16 pass: an earlier build's history opens and records as before and heals at the next DA that commits and the next correlation; classic route unchanged (T16; the omics routes are outside the table by construction, not run). Met, omics guard by reading |
| A6 | T18 head: replay from `correlation_filtered`, history-state comparison empty. Met |
| A7 | NEWS, ledger row, this file. Met |

## Departures from the ticket's text (rulings of 2026-11-27; step0.md)

T18 and T17 assert both mismatch readings (the history's states: empty; the origin's whole map: exactly
`ruv_corrected/missing` and `ruv_covariates/missing`); the completion's flag is passed only when TRUE (fixed-formals stubs at
`test-prot-05b-norm-module-contracts.R:2302`, `:2369`); `ruv_covariates_summary` is cleared in `runProtNormResetWorkflow`
(limpa only), not `resetProtNormOutputs` (pinned by KIT-122 F10; T9 shows the other three routes untouched); the leftover
file is removed at the call sites through one helper (T8c: apply, correlation Skip, RUV Skip; the managed finish in T8m);
T1 asserts the record's steps and T2 the report; T6 produces an enrichment first (the test-kit213 recipe) and asserts the
product present before and absent after; a focused recording-runner row (T19) and the managed `finish` row (T8m). Own
departures: T8b's refresh half is not tested (the KIT-206 head guard refuses a refresh from a head that is not the exported
one, so the rewrite cannot be reached in the re-apply state); the refresh code is written per the ticket and covered by
reading. Dispatch of the managed Skip is tested with a recording runner (T7b: one composite dispatch, none built in the
application).

## Removed NEWS sentences (A7)

KIT-273: "The file follows what the analysis used: pressing "Skip RUV" after a RUV run does not yet undo that run (the next
analysis still uses the factors or the corrected matrix, whatever the tab says), so the file is then still delivered."
KIT-266: "One limit: if you re-apply the sample-correlation filter after the RUV step, Figure 1 is rebuilt with two columns, but
the report's Methods and Results text still describe the earlier RUV step as part of the analysis (this was so before this
change and is not changed by it)." The ledger's VAL-379 "Known limit" text is left as history; the new row supersedes it.

## Stored values that move

None frozen. Declared hand-written edit: `tests/testthat/test-workflow-integration-input.R:703`. R2 showed no stored
capture, no `gui-entry` lane and no frozen `dia_normalization_retry` field moving (all `same`). The end-of-project re-record
list is unchanged by this ticket.

## Drafts for the wording pass (VAL-41)

The three user texts in `protNormRuvSkipTexts()` (Skip, second sentence "Run the differential analysis again.", correlation
re-applied) and the busy refusal, and the NEWS entry. The tests assert them as drafts.

## For the change report

What was wrong: on DIA with limpa, after a RUV run, applying the correlation filter again left the earlier RUV step in the
history, so the analysis record, the study parameters, the report (Methods paragraph, Results text) and the session file
described covariates the analysis did not carry; the RUV tab still showed the run's sweep and text; the archive said covariates
mode and k. "Skip RUV" after a covariates or correction run committed nothing, so the next differential analysis still used
the W terms or the corrected matrix while the tab said it carried none. Affected: users of the RUV step on this route (new in
0.9, tag v0.9.0 and later; none on `main`). Code (base 38443af8 -> fix): `commitWorkflowReferenceState`
(`R/utils_workflow_dia_state.R:853`, history written at `:891`), the Skip branch (`R/mod_prot_norm_workflow_helpers.R:80-91`),
`completeProtNormCorrelationWorkflow` (`R/mod_prot_norm_correlation_helpers.R:382`), the archive
(`R/utils_workflow_dia_sessions.R:1972`, refresh `:1800`), the replay (`R/utils_workflow_analysis_record.R:1058-1066`).
What changed: history and record list only the steps that produced the current table; Skip takes the RUV result out; the
tab, holders, archive fields and texts follow; no number any step computes changed (after Skip the next differential analysis
reads the table without RUV, as the tab said). Evidence: R1, R1b, R2, R3 above (fixture: the small recipe, in process and
through the worker); limits: no real data, no managed chain end to end.

## Not checked

The managed chain end to end (module test server, KIT-266's file) for the Apply hand-over and the managed Skip beyond the
recording-runner rows (lines read: `R/mod_prot_norm_workflow_helpers.R` dispatcher and observers); the archive refresh rewrite;
lipidomics, metabolomics and transcriptomics histories (outside the table; the R2 same verdicts of their files only);
Windows; real data; a Skip with the RUV design panel on in a real session; the wording.

## Unrelated discoveries

`test-module-ci-trx-qc.R` errors identically on base and head in R2 (`filtered[, design_order]`, subscript out of bounds);
left untouched.
