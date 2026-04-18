# Metabolomics Summary Module Seam Map

## Goal

Introduce low-risk top-level seams in R/mod_metab_summary.R while keeping the public mod_metab_summary_ui() and mod_metab_summary_server() entry points behaviorally stable.

## Current Position In The Flow

- target: `R/mod_metab_summary.R`
- classification: `complete`
- checkpoint reached: `wave 2 summary observer-shell helper live apply checkpoint`
- next step: `Treat R/mod_metab_summary.R as the completed public wrapper identity for this manual target and archive this handover; keep tests/testthat/test-metab-05-summary-module-characterization.R as the regression surface and do not reopen the wrapper unless a real regression appears.`

## Existing Safety Net

- `tests/testthat/test-metab-05-summary-module-characterization.R`

## Notes

Manual bucket 0 metabolomics summary-module stabilization target.

- This target previously had no active handover; this file now records the
  first metabolomics summary stabilization stop point.
- Classification refresh on April 18, 2026 kept
  `R/mod_metab_summary.R`
  at `635` lines with `2` top-level functions, a `505` line largest
  top-level function, and labels `high-risk-wrapper` /
  `needs-seam-introduction` before this checkpoint.
- The focused characterization gate now lives in
  `tests/testthat/test-metab-05-summary-module-characterization.R`.
- The first bounded seam now lives in-place in
  `R/mod_metab_summary.R`
  as `setupMetabSummaryTemplateStatusOutput()`.
- `output$template_status`
  now routes through `setupMetabSummaryTemplateStatusOutput()`, so the
  report-template status render registration stays in one top-level stop point
  instead of remaining inline in `mod_metab_summary_server()`.
- The focused characterization gate now freezes the current template-status
  text variants for missing omic directories, missing template files, and
  present template files, plus the wrapper handoff from
  `mod_metab_summary_server()` into
  `setupMetabSummaryTemplateStatusOutput()`.
- The focused gate reran green after this checkpoint via direct
  `testthat::test_file()` invocation because this worktree does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- Post-checkpoint classification now keeps
  `R/mod_metab_summary.R`
  at `643` lines with `3` top-level functions, a `486` line largest
  top-level function, and labels `high-risk-wrapper` /
  `needs-seam-introduction`.
- The second bounded seam now lives in-place in
  `R/mod_metab_summary.R`
  as `setupMetabSummaryBootstrapOutputs()`.
- `output$session_summary`
  and the initial `output$report_ready` bootstrap registration now route
  through `setupMetabSummaryBootstrapOutputs()`, so the wrapper tail stops
  owning that default-output cluster inline.
- The focused characterization gate now freezes the initial summary text, the
  initial `report_ready` hidden-state bootstrap, and the wrapper handoff from
  `mod_metab_summary_server()` into
  `setupMetabSummaryBootstrapOutputs()`.
- The focused gate reran green after this checkpoint via direct
  `testthat::test_file()` invocation because this worktree does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- Post-checkpoint classification now keeps
  `R/mod_metab_summary.R`
  at `650` lines with `4` top-level functions, a `479` line largest
  top-level function, and labels `high-risk-wrapper` /
  `needs-seam-introduction`.
- The third bounded seam now lives in-place in
  `R/mod_metab_summary.R`
  as `runMetabSummaryExportSessionObserverShell()`.
- The `export_session_state` observer now routes through
  `runMetabSummaryExportSessionObserverShell()`, so the export path assembly,
  payload snapshot, `saveRDS()` call, and success/error notification tail are
  no longer owned inline by `mod_metab_summary_server()`.
- The focused characterization gate now also freezes the export-session
  payload, success/error notification and logging contract, and the wrapper
  handoff from `mod_metab_summary_server()` into
  `runMetabSummaryExportSessionObserverShell()`.
- The focused gate reran green after this checkpoint via direct
  `testthat::test_file()` invocation because this worktree does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- Post-checkpoint classification now keeps
  `R/mod_metab_summary.R`
  at `696` lines with `5` top-level functions, a `457` line largest
  top-level function, and labels `high-risk-wrapper` /
  `needs-seam-introduction`.
- The fourth bounded seam now lives in-place in
  `R/mod_metab_summary.R`
  as `runMetabSummaryGithubPushObserverShell()`.
- The `push_to_github` observer now routes through
  `runMetabSummaryGithubPushObserverShell()`, so the GitHub `req()` gate,
  progress wrapper, option registration, push dispatch, success/error
  notification and logging contract, and session-summary update are no longer
  owned inline by `mod_metab_summary_server()`.
- The focused characterization gate now also freezes the GitHub observer
  shell's option payload, push dispatch arguments, success/error reporting
  contract, and the wrapper handoff from `mod_metab_summary_server()` into
  `runMetabSummaryGithubPushObserverShell()`.
- The focused gate reran green after this checkpoint via direct
  `testthat::test_file()` invocation because this worktree does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- Post-checkpoint classification now keeps
  `R/mod_metab_summary.R`
  at `748` lines with `6` top-level functions, a `437` line largest
  top-level function, and labels `high-risk-wrapper` /
  `needs-seam-introduction`.
- The fifth bounded seam now lives in-place in
  `R/mod_metab_summary.R`
  as `runMetabSummaryGenerateReportObserverShell()`.
- The `generate_report` observer now routes through
  `runMetabSummaryGenerateReportObserverShell()`, so the template
  presence/download branch, `RenderReport()` dispatch, `report_ready`
  bootstrap, download-handler registration, session-summary update, and
  success/error notification tail are no longer owned inline by
  `mod_metab_summary_server()`.
- The focused characterization gate now also freezes the package-template
  render success contract, the invalid-project-dir guard, the render failure
  reporting contract, and the wrapper handoff from
  `mod_metab_summary_server()` into
  `runMetabSummaryGenerateReportObserverShell()`.
- The focused gate reran green after this checkpoint via direct
  `testthat::test_file()` invocation because this worktree does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- Post-checkpoint classification now keeps
  `R/mod_metab_summary.R`
  at `874` lines with `8` top-level functions, a `311` line largest
  top-level function, and labels `high-risk-wrapper` /
  `needs-seam-introduction`.
- The sixth bounded seam now lives in-place in
  `R/mod_metab_summary.R`
  as `runMetabSummaryCopyToPublicationObserverShell()`.
- The `copy_to_publication` observer now routes through
  `runMetabSummaryCopyToPublicationObserverShell()`, so the fallback
  `study_parameters.txt` creation, workflow-data or file-backed
  `design_matrix`/`contrasts_tbl` resolution, global `project_dirs`
  assignment, `copyToResultsSummary()` dispatch, and success/error
  notification tail are no longer owned inline by
  `mod_metab_summary_server()`.
- The focused characterization gate now also freezes the fallback
  `study_parameters.txt` creation contract, the file-backed copy success
  contract including global `project_dirs` assignment, the copy error
  reporting contract, and the wrapper handoff from
  `mod_metab_summary_server()` into
  `runMetabSummaryCopyToPublicationObserverShell()`.
- The focused gate reran green after this checkpoint via direct
  `testthat::test_file()` invocation because this worktree does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- Post-checkpoint classification now keeps
  `R/mod_metab_summary.R`
  at `978` lines with `9` top-level functions, a `215` line largest
  top-level function, and labels `high-risk-wrapper` /
  `needs-seam-introduction`.
- The seventh bounded seam now lives in-place in
  `R/mod_metab_summary.R`
  as `runMetabSummarySaveWorkflowArgsObserverShell()`.
- The `save_workflow_args` observer now routes through
  `runMetabSummarySaveWorkflowArgsObserverShell()`, so the state-manager S4
  selection, `config_list` global assignment, `createWorkflowArgsFromConfig()`
  dispatch, integration-object persistence, fallback `study_parameters.txt`
  warning path, and success/error notification tail are no longer owned inline
  by `mod_metab_summary_server()`.
- The focused characterization gate now also freezes the save-workflow
  observer shell's S4 export contract, fallback warning contract, and the
  wrapper handoff from `mod_metab_summary_server()` into
  `runMetabSummarySaveWorkflowArgsObserverShell()`.
- The focused gate reran green after this checkpoint via direct
  `testthat::test_file()` invocation because this worktree does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- Post-checkpoint classification now keeps
  `R/mod_metab_summary.R`
  at `1074` lines with `10` top-level functions, a `101` line largest
  top-level function, and labels `high-risk-wrapper` /
  `needs-seam-introduction`.
- The eighth bounded seam now lives in-place in
  `R/mod_metab_summary.R`
  as `setupMetabSummaryServerBootstrapState()`.
- The initial `experiment_label` autopopulate and `reactiveValues()`
  bootstrap now route through `setupMetabSummaryServerBootstrapState()`, so
  the wrapper head no longer owns that initialization block inline.
- The focused characterization gate now also freezes the bootstrap-state
  helper's experiment-label update contract, default reactive-state payload,
  and the wrapper handoff from `mod_metab_summary_server()` into
  `setupMetabSummaryServerBootstrapState()`.
- The focused gate reran green after this checkpoint via direct
  `testthat::test_file()` invocation because this worktree does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- Post-checkpoint classification now keeps
  `R/mod_metab_summary.R`
  at `1084` lines with `11` top-level functions, a `93` line largest
  top-level function, and labels `high-risk-wrapper` /
  `needs-seam-introduction`.
- The ninth bounded seam now lives in-place in
  `R/mod_metab_summary.R`
  at line `985`
  as `registerMetabSummaryServerObservers()`.
- The remaining inline `observeEvent()` wiring for
  `save_workflow_args`,
  `copy_to_publication`,
  `generate_report`,
  `push_to_github`,
  and `export_session_state`
  now routes through `registerMetabSummaryServerObservers()`, so
  `mod_metab_summary_server()` no longer owns that observer-registration block
  inline.
- The focused characterization gate now also freezes the observer-registration
  helper's event ordering, the `ignoreInit = TRUE` copy-to-publication
  contract, and the wrapper handoff from `mod_metab_summary_server()` into
  `registerMetabSummaryServerObservers()`.
- The focused gate reran green after this checkpoint via direct
  `testthat::test_file()` invocation because this worktree does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- Post-checkpoint classification now reports
  `R/mod_metab_summary.R`
  at `1104` lines with `12` top-level functions, an `87` line largest
  top-level function, `0` observers, `2` renderers, and label `review`.
- All summary-module observer handlers and their registration wiring now
  delegate through top-level seams; the next safe target is a staged
  exact-source wave for the stabilized summary helpers, and this iteration
  stops before staging that wave.
- The first exact-source metabolomics summary server-setup wave now lives in
  `tools/refactor/manifest-metab-summary-wave1.yml`
  and stages the wrapper setup cluster from
  `R/mod_metab_summary.R`
  into
  `tools/refactor/staging/wave1_metabolomics_summary_server_setup/R/mod_metab_summary_server_helpers.R`
  without rewriting live sources before review.
- The staged wave covers
  `setupMetabSummaryTemplateStatusOutput()`,
  `setupMetabSummaryBootstrapOutputs()`,
  `setupMetabSummaryServerBootstrapState()`,
  and
  `registerMetabSummaryServerObservers()`.
- `tools/refactor/verify_refactor.R` passed for
  `tools/refactor/manifest-metab-summary-wave1.yml`
  before staging the reviewed server-setup helper target.
- The focused summary-module loader in
  `tests/testthat/test-metab-05-summary-module-characterization.R`
  now widens the filename-coupled seam surface by loading
  `R/mod_metab_summary_server_helpers.R`
  before
  `R/mod_metab_summary.R`,
  so the focused gate survives the live server-helper apply boundary.
- The focused summary-module gate reran green after the staged wave via direct
  `testthat::test_file()` invocation because this worktree does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- Live classification remains
  `review`
  at `1104` lines with `12` top-level functions, an `87` line largest
  top-level function, `0` observers, and `2` renderers because the reviewed
  wave is staged only and
  `R/mod_metab_summary.R`
  remains the live wrapper.
- The next safe stop point is a reviewed live apply of wave 1 into
  `R/mod_metab_summary_server_helpers.R`,
  with `DESCRIPTION` `Collate:` loading the new helper ahead of
  `R/mod_metab_summary.R`,
  then `tools/refactor/check_wave_apply.R`
  and the focused summary-module gate rerun before any further wrapper
  extraction.
- That reviewed live apply is now complete:
  `R/mod_metab_summary_server_helpers.R`
  now owns
  `setupMetabSummaryTemplateStatusOutput()`,
  `setupMetabSummaryBootstrapOutputs()`,
  `setupMetabSummaryServerBootstrapState()`,
  and
  `registerMetabSummaryServerObservers()`
  in live `R/`, so
  `R/mod_metab_summary.R`
  no longer carries those definitions inline.
- `DESCRIPTION` `Collate:` now loads
  `R/mod_metab_summary_server_helpers.R`
  ahead of
  `R/mod_metab_summary.R`,
  matching the new live helper boundary.
- `tools/refactor/check_wave_apply.R` passed for
  `tools/refactor/manifest-metab-summary-wave1.yml`
  after the live apply.
- The focused summary-module gate reran green after the live apply via direct
  `testthat::test_file()` invocation because this worktree does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- Post-apply classification now keeps
  `R/mod_metab_summary.R`
  at `962` lines with `8` top-level functions, an `87` line largest
  top-level function, `0` observers, `0` renderers, and label `review`.
- The next safe stop point is a reviewed staged wave for the remaining
  observer-shell cluster in
  `R/mod_metab_summary.R`,
  with the focused summary-module gate rerun after that next bounded
  checkpoint.
- The next reviewed exact-source observer-shell manifest now exists at
  `tools/refactor/manifest-metab-summary-wave2.yml`,
  targeting the remaining helper block from
  `runMetabSummaryExportSessionObserverShell()` through
  `runMetabSummaryGithubPushObserverShell()` into the new helper target
  `R/mod_metab_summary_observer_helpers.R`.
- `Rscript tools/refactor/verify_refactor.R --manifest
  tools/refactor/manifest-metab-summary-wave2.yml` passed on April 18, 2026,
  so the reviewed selectors for the remaining observer-shell cluster are
  verified against the current `R/mod_metab_summary.R` source tree before any
  live apply.
- The reviewed observer-shell manifest has now been staged into
  `tools/refactor/staging/wave2_metabolomics_summary_observer_shell_helpers/`,
  producing `R/mod_metab_summary_observer_helpers.R` plus
  `collate-metab-summary-wave2.txt` without rewriting live sources.
- Staged-output review confirmed that `R/mod_metab_summary_observer_helpers.R`
  holds the exact-source observer-shell helper block from
  `runMetabSummaryExportSessionObserverShell()` through
  `runMetabSummaryGithubPushObserverShell()`, and the collate fragment orders
  the new observer helper before the existing
  `R/mod_metab_summary_server_helpers.R`.
- The focused summary-module loader now widens the filename-coupled seam
  surface by loading `R/mod_metab_summary_observer_helpers.R` before
  `R/mod_metab_summary_server_helpers.R` and `R/mod_metab_summary.R`, so the
  focused gate survives the future live observer-helper apply boundary.
- The focused summary-module gate reran green after the staged wave via direct
  `testthat::test_file()` invocation because this worktree does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- Classification refreshed on April 18, 2026 keeps
  `R/mod_metab_summary.R`
  at `962` lines with `8` top-level functions, an `87` line largest
  top-level function, `0` observers, `0` renderers, and label `review`
  because the reviewed wave is staged only and
  `R/mod_metab_summary.R`
  remains the live wrapper.
- The next safe stop point now advances from staging the observer-shell helper
  block to applying the reviewed wave-2 manifest live, updating `DESCRIPTION`
  so `R/mod_metab_summary_observer_helpers.R` collates before
  `R/mod_metab_summary_server_helpers.R` and `R/mod_metab_summary.R`, then
  rerunning `tools/refactor/check_wave_apply.R` and the focused summary-module
  gate.
- The reviewed exact-source observer-shell wave is now live via
  `tools/refactor/manifest-metab-summary-wave2.yml`
  into
  `R/mod_metab_summary_observer_helpers.R`,
  moving
  `runMetabSummaryExportSessionObserverShell()`,
  `runMetabSummarySaveWorkflowArgsObserverShell()`,
  `runMetabSummaryCopyToPublicationObserverShell()`,
  `runMetabSummaryGenerateReportObserverShell()`,
  and
  `runMetabSummaryGithubPushObserverShell()`
  out of
  `R/mod_metab_summary.R`
  without hand-rewriting helper bodies.
- `DESCRIPTION` `Collate:` now loads
  `mod_metab_summary_observer_helpers.R`
  before
  `mod_metab_summary_server_helpers.R`
  and
  `mod_metab_summary.R`,
  with the live load-order artifact recorded in
  `tools/refactor/collate-metab-summary-wave2.txt`.
- `Rscript tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-summary-wave2.yml`
  passed after the live apply.
- The focused summary-module gate reran green again after the live apply via
  direct `testthat::test_file()` invocation because this worktree does not
  include `renv/activate.R` for `tools/test_with_renv.R`.
- Post-apply classification now records
  `R/mod_metab_summary.R`
  at `165` lines with `2` top-level functions, an `87` line largest
  top-level function, `0` observers, `0` renderers, and label `review`;
  the live file now only carries
  `mod_metab_summary_ui()`
  plus the review-frozen
  `mod_metab_summary_server()`
  entry shell.
- Treat the manual metabolomics summary-module target as complete:
  helper bodies now live in
  `R/mod_metab_summary_server_helpers.R`
  and
  `R/mod_metab_summary_observer_helpers.R`,
  while
  `R/mod_metab_summary.R`
  remains the public wrapper identity for this module and bucket 0 can move on
  to the next manual target.
