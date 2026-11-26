# KIT-274 outcome

R1 run 45723980071 on tests-only 050a0402 (tables: r1-tests.tsv, r1-failures.txt).
Fail by value as predicted: T1 (3 failures), T2 (3), T3 (4), T4 (4), T5 (1), T6 (4), T7 (3), T12 (3). Pass as predicted: T8
(digest of the written file ba839e69b0d2c199d37a7e68ed0207679b0138f7751c5a6897f4e3c233da6108, product trx.tmm_logcpm,
normalised, not corrected), T9 (k 2, corrected, trx.ruviii), T11. X1 base line printed: `ruv_applied=FALSE export=ok entry=listed k=0 text=RUVIII (ruv::RUVIII) was applied ... k = 0 ...`.
MISSED PREDICTION: T10 (guard, predicted pass on base) failed: route A and route B both listed TRUE on the base.
Likely cause (inferred, not shown): k274Run builds a fresh omicNormReactiveState() per run where test-kit264 E1/E2 share one `nd`, so the
second run may not have run; the test has no precondition that the second run made a new trx_ruvs entity. No error text in the log.
Stop rule: a missed prediction stops the work; reported, not corrected.
Not checked: everything after R1 (the fix is drafted in the worktree, uncommitted; R2; R CMD check).

## R1b (VAL-54 instrument re-send 1 of 3), after the coordinator's ruling (R1 above is kept as it happened)

Ruling quoted: "treated as an instrument defect, re-send 1 of 3 (VAL-54) ... the behaviour T10 guards (a later run supersedes the
listed export) is what KIT-264's E1/E2 assert on the same base, and that file passed 85/0 in your own R1 with one shared `nd`; your T10
differs only in building a fresh omicNormReactiveState() per run, and it has no precondition showing the second run happened. So the run
gave no result on the guarded behaviour." Terms: tests-only; (1) one `nd` per project in k274Run; (2) a precondition in T10 that the
second run committed a new trx_ruvs entity; (3) no assertion changed. Sharing `nd` changes nothing T1 to T9, T11, T12 exercise (only T10
calls k274Run). Predictions for R1b: T10 passes on base with the precondition met; all others as in R1 (T1-T7, T12 fail by value; T8, T9,
T11 pass; X1 base line as in R1). If T10 still fails with the precondition met: a result, stop. If the precondition fails: hand back.
T4 deviation (accepted by the coordinator): the press is made in the summary hook, not through `ruviii = "group"`, so the product set is
read before and after it; it asserts the same observables as the plan register's P04.

R1b result: run 45726224100 on tests-only 60565dd4 (r1b-tests.tsv, r1b-failures.txt): 53 passed, 25 failed, 0 errors (185.5 s). As predicted, except that T11 was not in the file (see the correction below):
T1 (3), T2 (3), T3 (4), T4 (4), T5 (1), T6 (4), T7 (3), T12 (3) fail by value; T8, T9, T10 pass. T10: both preconditions met
(the second run committed a new trx_ruvs entity), `route A listed: FALSE route B listed: FALSE` on the base: R1's T10 failure was the
instrument (no second run).

Correction (implementation audit I01), 2026-11-26: commit 60565dd4, ruled tests-only with no assertion changed, also deleted
T11 (050a0402 :480-484, two `saved_entity_scope` expectations); its message, the file header and the R1b and R2
records did not say so, and the records stating that T11 passed in R1b and R2 were wrong. T11 passed in R1 (51 = 46 +
T10's 3 + T11's 2). R1b (53/25) and R2 (82/0) ran without it. Why it was removed, from my own commits: by mistake, while
restructuring T10. My edit replaced the text from the start of T10 to the next section header, and T11 sat between them; it
never failed anywhere and sharing one `nd` did not touch it (T11 uses no k274Run). The behaviour it guarded is asserted by
test-workflow-trx-declarations.R:410-420 (refusals 4 and 5), which passed on both sides of R2. No other expectation changed in
60565dd4: T10's assertions are as in R1, with two preconditions added as ruled. The deletion is not accepted (coordinator's
ruling): T11 was restored verbatim from 050a0402 (nothing differs, not one character) in 27a42adb, and run in R2b (below).

Cost, measured from the jobs API (implementation audit I02): R1 15.98, R1b 14.68, R2 352.32 job-minutes; total 382.98 =
25.53 core-hours (expected 28.20, ceiling 103.0), before R2b.

## R2: run 45728402078, head 87106bc3 (fix, merged with staging 7f17547f) against base edb15841 (the same tree, fix reverted)

The base is a revert commit on kit-274 itself (only that branch may be pushed); the head is the revert of the revert (tree equal to the
merge ba2ad567). Scope recounted on the head: the ticket's grep gives 36 files (34 at 75acd421 + test-kit213-study-parameters-labels.R from
train 23 + the new file) + test-omics-module-globals.R by instruction = 37. Tables: r2-compare.tsv, r2-tests-head.tsv, r2-tests-base.tsv,
r2-compare-failures.txt, r2-check-compare.txt. Result: same 36, fixed 1, new-failure 0, message-changed 0, count-changed 0. "Clean under the rule".
As predicted. The new file went 53/25 to 82/0 (the predicted 78 was the base's total; the head runs 4 more expectations inside the `if (length(notes))` blocks that the refusal makes true). test-kit264-trx-export-current.R 85/0 both sides. The pre-existing failures (test-omic-summary-module 28, test-omics-summary-core 25, test-workflow-trx-declarations 31, test-module-ci-trx-summary 1, test-workflow-kit-harness 1, test-workflow-omics-binding-steps 2, test-workflow-operation-registry 2, test-omics-module-globals 1) are identical on both sides with identical messages (stored laptop references, as the brief says): unchanged. Lipidomics and metabolomics record and replay tests: same.
R CMD check: base and head each 0 ERRORs, 3 WARNINGs, 3 NOTEs; findings identical (6), none only on head.
X1 line before (R1): `export=ok entry=listed k=0 text=RUVIII (ruv::RUVIII) was applied ... k = 0 ...`; after (R2 head): `export=trx_ruviii_no_ruv/trx_error/rlang_error/error/condition: no RUVIII export is made: this analysis applied no correction for unwanted variation (k = 0; the design none: no correction for unwanted variation), so RUVIII would return the TMM logCPM unchanged; the integration uses the TMM logCPM. entry=absent ...`.
Stored values: none moved; none edited (test-kit264 file unedited, no stored reference touched).

## Acceptance

| # | Criterion | Result |
|---|---|---|
| A1 | tests-only commit: T1-T7, T12 fail by value; T8, T9 pass; T10 passes in R1b; T11 passed in R1 and was then deleted (I01), restored and run in R2b; base X1 line printed | met (R1, R1b, R2b; T10's first-run failure was an instrument defect, ruled by the coordinator) |
| A2 | every row passes, in process, managed and on replay | met (R2 head 82/0; T4 managed on the chain mocks, T7 replay) |
| A3 | R2 base against head: every file same except the new file; failures identical; check identical; X1 line before and after | met |
| A4 | refusal and replay statement are the decided ones (drafts for the wording pass) | met (wording.md) |
| A5 | integration input and file unchanged in the no-RUV state | met: T8 passes on both sides; the written file's digest ba839e69b0d2c199d37a7e68ed0207679b0138f7751c5a6897f4e3c233da6108 on R1, R1b and R2 head |
| A6 | NEWS, ledger row, outcome with "For the change report" | met (ledger row VAL-<merge>) |

## Not checked
- The real managed worker (T4 uses the chain mocks, the journey helper's runner); the real application and tab by eye; Windows; real data.
- That `study_parameters.txt` and the workbook take the export only through `buildWorkflowAnalysisRecord` (by reading at the writer's call only).
- A project whose report was generated before keeps its old template copy (press Save again); a leftover ruviii_before_and_after.png is still delivered by Copy (F1).
- I07 (a refused Save leaves an earlier Save's integration input) and I06/N5 are candidates, not changed. I07 candidate as the ticket records it: a falsifier first (Save, run with another k, DA, Save refused, then read that the earlier integration/transcriptomics folder is still there); the owner's ruling that the refusal removes the earlier input like the other three classes; one more class in the `inherits(error, ...)` of summaryFailed (R/utils_workflow_summary.R, tests in test-kit207-summary-head-guard.R); a read of the integration reader for a source-against-head check; the same for proteomics and lipidomics.

## For the change report
- New to 0.9 (all of transcriptomics is); none on `main`. Who: a user of the RUV design section (the installer's test features; off in the shipped default) who chose the design none, or a design whose k rule found 0, and pressed "Run RUVIII Export".
- What was wrong: the export was accepted with k = 0 and 0 controls; the analysis record, workbook and study parameters listed it, and the report said RUVIII "was applied ... with the RUVs inferred negative controls and k = 0, for the integration export only" beside its own "No correction for unwanted variation was applied". With k = 0 RUVIII returns its input; the integration file in that state is the TMM logCPM. No number was wrong; the description of the analysis was.
- Code (base 75acd4219 and 7f17547f; fix commit 971b54bf): R/utils_workflow_trx_declarations.R acceptTrxRuviiiExport (:661) and runTrxReferenceRuviiiExport (:1332) now call the new trxAssertRuviiiExportApplied (the integration's own trxIntegrationRuvState condition); R/utils_workflow_summary.R trxRuviiiExportParameters (:714, the managed request) and runWorkflowSummaryExportStep (:618, managed setup inside the failure handler); R/utils_workflow_analysis_record.R export loop (:546, by the record's own last trx_ruvs step) and the replay (skip with statement and warning of class workflow_replay_export_skipped); inst/reports/transcriptomics/transcriptomics_report.rmd:49 (guard).
- What changed in results: no value. The refusal, no listing, and the replay skip.
- In tag v0.9.0 the same acceptance exists (`v0.9.0:R/utils_workflow_trx_declarations.R:661-673`), with the design none
  (`:491`, `ruv_applied = d$applied`) and the report line `export <- record$analyses$trx_ruviii_export` (`:47`; its figure shown
  whenever the file exists, `:126`): a v0.9.0 build that used the RUV design section's none could record and report the k = 0
  export the same way. v0.9.0 is not a release; `main` (0.4.1.2) has no transcriptomics.
- Evidence: runs 45723980071 (R1), 45726224100 (R1b), 45728402078 (R2); the X1 line before and after above.

## R2b: run 45739818528 (T11 restored; coordinator's ruling on implementation audit I01)

Head e4335b58 (staging ec27f69e merged, fix present, T11 restored in 27a42adb) against base 7ca71057 (the same tree with the fix reverted, on kit-274),
the new file only, one shard per side, no check. Table: r2b-compare.tsv. Result as predicted (predictions.md, pushed before the request):
base 55 passed / 25 failed, head 84 passed / 0 failed, verdict fixed; T11 passes on both sides (no T11 failure in either failure list);
T8's written-file digest ba839e69b0d2c199d37a7e68ed0207679b0138f7751c5a6897f4e3c233da6108 again. Merge of origin/v0.9-staging ec27f69e into kit-274: no conflict
(NEWS and the ledger merged automatically; my NEWS bullet stays first under Fixes, my ledger row stays last as VAL-<merge>; nothing under R/, tests/ or inst/ conflicted).
Job-minutes of R2b were not read from the API (one shard per side; ceiling 40).
