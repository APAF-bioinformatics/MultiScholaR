# KIT-274 step 0: the ticket's citations re-read on the base

Base: origin/v0.9-staging 75acd4219 (merge train 21 of 2026-11-26; contains train 19 with KIT-264, trains 20 and 22).
Branch kit-274 in the worktree /home/doktersmol/Projects/APAF/MultiScholaR-kit-274. Read by reading only; no R.
Bounded item: the ticket KIT-274 r1.2 (owner decisions VAL-318, VAL-332); the plan register
AUDIT-KIT-274-plan-20261126 governs where the ticket is unclear (P01, P02, P04 applied).

## Lines re-read (ticket cites train l b6b4ed6b; base is 75acd4219)

| Ticket cites | On the base | Moved? |
|---|---|---|
| utils_workflow_trx_declarations.R:661-673 `acceptTrxRuviiiExport` (owner check :662) | :661-673, same text | no |
| :1332-1342 `runTrxReferenceRuviiiExport` | :1332-1342 | no |
| :990-997 `trxIntegrationRuvState` (applied = ruv_applied and k above 0) | :990-997 | no |
| :644-655 `trxRuviiiExportOperation`; :442-474 design none and k rule reasons; :491 `ruv_not_applied_reason` | same | no |
| utils_workflow_summary.R:618-643 `runWorkflowSummaryExportStep` (managed branch :621-634, `step$parameters` :624 outside any handler; in-process tryCatch :635-642) | :618-643, same | no |
| utils_workflow_summary.R:714-717 `trxRuviiiExportParameters`; :703-708 `trxSummaryIntegrationRequires` | same | no |
| utils_workflow_analysis_record.R:346-351 `workflowRecordExportIsCurrent`; export loop :546-549 | predicate :346-351; loop at :546-549 (the `next` guard on :547) | no |
| record replay loop for exports :1077-1092; statements pattern :1093-1100 | loop `for (analysis in ...)` at :1083-1096; `replayed` list :1098; design_statement :1099-1102; correlation_statement :1105-1109 | yes, shifted by about 6 lines (KIT-258's marker, added by train 20, is above in `workflowReplayRuleMarkers`) |
| transcriptomics_report.rmd:44-46 `no_ruv`, :49 `export <-` | :44-46, :49 | no |
| :115-124 the "ruviii-text" chunk | chunk opens at :117 (text :118-125) | yes, 2 lines |
| :170-175 the "ruviii-figure" chunk | :172-175 | no |
| test-kit264-trx-export-current.R X1 :496-517, helpers :47-152 | X1 at :496-517 region, helpers :47-152 (read; not edited) | no |
| helper-module-ci-trx-tabs.R:200-254 (journey, managed runner) | same | no |

Between origin/merge-train-20261126l (a2f69738) and the base the only change in the four target files is the
KIT-258 replay-marker text in utils_workflow_analysis_record.R (:663, :701-706); not this ticket's hunks.

## Confirmed

- The three call sites exist as the ticket says: in process and on replay `runTrxReferenceRuviiiExport` (replay entry
  `ruviii_export` in `trxRouteReplay`, :747-748); managed `trxRuviiiExportParameters` called at
  utils_workflow_summary.R:624 outside any handler; `acceptTrxRuviiiExport` the backstop every route ends in.
- The record export loop :546-549 keeps an export when `workflowRecordExportIsCurrent`; steps are built before it (the
  `steps` list, :416-482), each step's `outcomes` carry `ruv_applied` and `k` (trxRuvOutcomes :404; the declared design's
  outcomes :491-492).
- The template guard's place: `no_ruv` is defined at :46 before `export <-` at :49.
- Nothing in the design no longer fits; no stop.

## Scope recount (the ticket's grep, helper files excluded, data/benchmark excluded)

34 files on the base (as the ticket predicted for train l; the new file makes 35). The caller asked to include
`tests/testthat/test-omics-module-globals.R` (it is not in the grep: it scans `R/(mod|func)_(lipid|metab|omic|omics)*`
for globals and `<<-`; this ticket's edits are in `utils_*.R` files and add no `<<-`, no global): R2 scope = 34 + new
file + test-omics-module-globals.R = 36 files.

## Stored values, digests, pinned texts or oracles the change could move (whole tests/ tree, data/benchmark excluded)

None moves. Each is a legacy-grid export (k at least 1), which is unchanged:
- tests/testdata/gui-entry/trx_screen.json:1668, :1834, :2013, :2208, :2412: `ruviii_status` "RUVIII export complete (integration only; DA results unchanged)." after a legacy run.
- tests/testdata/mci/manifest.json (MCI-032.1 ruviii-export): legacy grid.
- tests/testthat/test-omic-summary-module.R:2548 (a technical-replicate refusal, class trx_ruviii_replicates_invalid, with the
  tab's own prefix; unchanged by the new class), :2570, :2601, :2632, :2678, :2716, :2840, :2996, :3036 (frozen status values "k = 2, 18 inferred...": legacy runs).
- tests/testthat/test-kit264-trx-export-current.R: X1 prints a different line, never fails (its own tryCatch; no edit, P01); header texts :10, :23-24 and rows E8 :471-489 name the applied sentence after a legacy run.
- tests/testthat/test-omic-gui-config.R:1728, :1742-1747 and test-omic-summary-module.R:262: the step's help text and its entries, unchanged (the help is not edited).
- tests/testthat/test-omics-summary-core.R:94, :105-109, :169: the entries stay in utils_workflow_summary.R; the file must contain no `shiny::` outside comments; the managed-branch edit uses none. The `ruviii` steps there drive legacy runs.
- tests/testthat/test-module-ci-trx-norm.R:84, test-trx-norm-module.R:641, :682: normalisation texts, not touched.
No test file other than the new one is edited.
