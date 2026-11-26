# KIT-276 outcome (branch kit-276; implementation, for the independent audit)

State: done and ready for its implementation audit. Branch tip is the last commit of `kit-276` (fix 79e77621; tests-only
base for the comparisons: branch `kit-276-tests` e6f83d43 = e77f096c + two mechanics repairs of the new file, no change
under `R/`). Staging branch point 38443af8. Nothing merged. Running total 87.64 of 172.00 core-hours (83.15 before R4) (restated expected 84.2).

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
| R1 45776354953 | tests-only e77f096c | 64 failures and 2 errors; 63 rows at the stated values; T8b and T14 errored on test mechanics (instrument defect 1, repaired in R1b); T2 failed on its third expectation only, its stated values not reached because the fixture's record had no DPC step (repaired after R2: the correction pass, below); `test-workflow-integration-input.R:703` failed as predicted; guards T12, T13, T16 passed |
| R1b 45780547662 | instrument repair 1 (T8b marks the owner a native session; T14's fixtures live in the test's frame) | T8b and T14 failed by value as predicted; 0 errors |
| R2 45784213557 | f847c433 vs fix 2f76cca9, 97 files, R CMD check both sides | 95 same, 2 fixed, 0 new-failure, 0 message-changed; check identical (0 ERROR, 3 WARNING, 3 NOTE both sides); `test-prot-norm-core.R` same 1555/307 (the ticket's 1547 was the older record; base and head agree); one expectation of the new file failed on the head: T2's last (below) |
| R3 45799233856 | THE TICKET'S ONE CORRECTION PASS (spent; coordinator's ruling on AUDIT-KIT-276-impl I02): T2 gives the report chunks the DPC step as test-kit227's E7 does, expectations unchanged, new file, e6f83d43 vs 72ba4edc | base 123 passed / 72 failed (T2 fails all three expectations); head 195 passed / 0 failed |

Instrument repairs used: 1 of 3 (R1b: T8b and T14 errored before any expectation; mechanics only, declared before sending).
Correction pass: SPENT at R3. T2 had run to completion and been judged before it was changed (R1/R1b: it failed on its third
expectation only, its first two passing on the base against the prediction; R2: it failed on the fixed code), so under VAL-54
the T2 change is a correction, not an instrument repair; I first described it as an instrument repair, and the coordinator
ruled otherwise. Its content is accepted (base fails all three expectations, head passes). Stop rule 3 (an unpredicted R2
verdict: the new file ended 194/1) was passed over: I repaired and sent R3 before reporting; it was disclosed at hand-back.
Cost: R1 1.89, R1b 1.49, R2 77.0, R3 2.8 = 83.15 core-hours (1,247.2 job-minutes).

## Acceptance criteria

| # | Evidence |
|---|---|
| A1 | Met at R1b except T2; T2 at R3 (R1: 63 rows at the stated values, T8b and T14 errored on mechanics and are shown failing by value in R1b) |
| A2 | R2/R3 head: all 195 rows pass in process and on the worker lane; the DA after each limpa sequence carries the state the record, tab and archive state. Met |
| A3 | R2: every file same except the new file and the declared edit; `test-prot-norm-core.R` same; check identical. Met |
| A4 | Skip is S1; texts are drafts (wording.md). Met |
| A5 | T14, T16 pass: an earlier build's history opens and records as before and heals at the next DA that commits and the next correlation; classic route unchanged (T16; the omics routes are outside the table by construction, not run). Met, omics guard by reading |
| A6 | T18 head: replay from `correlation_filtered`, history-state comparison empty. Met |
| A7 | NEWS, ledger row, this file. Met |

## Departures from the ticket's text (rulings of 2026-11-26; step0.md)

T18 and T17 assert both mismatch readings (the history's states: empty; the origin's whole map: exactly
`ruv_corrected/missing` and `ruv_covariates/missing`); the completion's flag is passed only when TRUE (fixed-formals stubs at
`test-prot-05b-norm-module-contracts.R:2302`, `:2369`); `ruv_covariates_summary` is cleared in `runProtNormResetWorkflow`
(limpa only), not `resetProtNormOutputs` (pinned by KIT-122 F10; T9 shows the other three routes untouched); the leftover
file is removed at the call sites through one helper (T8c: apply, correlation Skip, RUV Skip; the managed finish in T8m);
T1 asserts the record's steps and T2 the report; T6 produces an enrichment first (the test-kit213 recipe) and asserts the
product present before and absent after; a focused recording-runner row (T19) and the managed `finish` row (T8m). Own
departures: T8b's refresh half is not tested. KNOWN DEFECT awaiting the owner's word (AUDIT-KIT-276-impl I01): after
"Skip RUV" the Summary's refresh of the session file (`refreshDIANativeSummaryArchive`) keeps the export's
`ruv_optimization_result`, so `saveProtNormExportMetadataFiles` writes `ruv_optimization_results.RDS` back, the archive lists
it and a reopen restores the sweep; the refresh branch is reachable after a Skip (my earlier statement that it could not be
reached was wrong). The fix is prepared below, not applied. Dispatch of the managed Skip is tested with a recording runner (T7b: one composite dispatch, none built in the
application).

## Removed NEWS sentences (A7)

KIT-273: "The file follows what the analysis used: pressing "Skip RUV" after a RUV run does not yet undo that run (the next
analysis still uses the factors or the corrected matrix, whatever the tab says), so the file is then still delivered."
KIT-266: "One limit: if you re-apply the sample-correlation filter after the RUV step, Figure 1 is rebuilt with two columns, but
the report's Methods and Results text still describe the earlier RUV step as part of the analysis (this was so before this
change and is not changed by it)." The ledger's VAL-379 "Known limit" text is left as history; the new row supersedes it.

## Stored values that move

None frozen. Declared hand-written edit: `tests/testthat/test-workflow-integration-input.R:703` (against staging's unedited file the train will see `count-changed`, 240/0 to 241/0, not `fixed`: R2's base was the tests-only commit). R2 showed no stored
capture and no frozen `dia_normalization_retry` field moving (all `same`); `test-gui-entry-proteomics.R` was not in R2: by reading no
leaf of `prot_dia_limpa.json` moves. The end-of-project re-record
list is unchanged by this ticket.

## Drafts for the wording pass (VAL-41)

See `wording.md`. The three user texts in `protNormRuvSkipTexts()` (Skip, second sentence "Run the differential analysis again.", correlation
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
the non-native (legacy) export path, which still writes the tab's `input$ruv_mode` (`R/mod_prot_norm_session_helpers.R:144-145`; native is the production default); `test-gui-entry-proteomics.R`; Windows; real data; a Skip with the RUV design panel on in a real session; the wording.

## Unrelated discoveries

`test-module-ci-trx-qc.R` errors identically on base and head in R2 (`filtered[, design_order]`, subscript out of bounds);
left untouched.

## Further bounded correction (VAL-389, owner-allowed): R4 45863987179 is a MISS; stopped, no repair

Done as prepared: tests-only commit 7649808b (row T8d), predictions committed before the run, then the one line
`metadata["ruv_optimization_result"] <- list(NULL)` in the `run$ruv_mode`-is-NULL branch of `refreshDIANativeSummaryArchive`
(head 3e4cb94d), then one run (new file and `test-kit216-prot-session-after-enrichment.R`, base against head). Result: the
new file is 202 passed / 5 failed on BOTH sides (verdict `same`), `test-kit216-...` same 152/0/1. T8d fails identically on the
base and the head at its first non-precondition expectation: `protSummaryRefreshArchive` returns the text "The session file
was not updated: the current table is not the one the last session export saved ..." instead of NULL. That is the KIT-206
head guard (`assertProtNormExportHeadCurrent`) refusing the refresh: after the export at the RUV head the export's stored
correlation reference is the RUV state (the RUV publication sets `correlation_filtered_obj` to it,
`R/utils_workflow_dia_products.R:1824`), and after the Skip the head is `correlation_filtered`, whose table components (the
design with its W columns, or the corrected matrix) differ; the guard's exemption for RUV heads does not apply (the head is
no longer one). So in the sequence AUDIT-KIT-276-impl I01 describes, the refresh is refused (the user is told the session file
was not updated) and does not write the sweep back; the audit's step that the guard passes did not hold in this run. The
remaining expectations (archive fields, absent file) then ran on the unrefreshed export and failed on both sides for that
reason; no result for the line itself was produced. Per the instruction a miss is not repaired: the I01 line stays in the code
as written (harmless: it only clears a value in the branch that is reached when a refresh happens with no RUV description),
T8d stays as written and failing on both sides, and the case returns to the owner. A re-apply at the same threshold was not
examined. Cost R4: 67.3 job-minutes = 4.49 core-hours; total 87.64 of 172.00.
Note: KIT-277 (branch kit-277, not merged) also edits `refreshDIANativeSummaryArchive` a few lines away (it takes the
annotation and FASTA tables from the current bus); whichever lands second resolves by hand. I did not merge it.
