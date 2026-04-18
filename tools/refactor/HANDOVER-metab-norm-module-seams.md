# Metabolomics Normalization Module Seam Map

## Goal

Introduce bounded top-level seams in `R/mod_metab_norm.R` while keeping the
live metabolomics normalization module behavior frozen behind the existing
public entry points.

## Current Position In The Flow

- target: `R/mod_metab_norm.R`
- classification: `direct-extraction-ready`
- checkpoint reached: `wrapper entrypoint apply wave`
- next step: `Treat R/mod_metab_norm.R as the completed breadcrumb stub for this manual target and archive this handover.`

## Existing Safety Net

- `tests/testthat/test-metab-03b-norm-module-helper-characterization.R`

## Notes

Manual bucket 0 metabolomics normalization module stabilization target.

- This target previously had no active handover; this file now records the
  current metabolomics normalization module seam stop point.
- Classification refresh on April 17, 2026 now reports
  `R/mod_metab_norm.R`
  at `69` lines; the current classifier output reports `0` top-level
  functions, `0` module servers, and label
  `direct-extraction-ready`.
- Historical seam-introduction bullets below preserve the original in-file
  introduction points; the current live helper ownership after the wave-2
  through wave-15 applies is summarized near the end of this handover.
- The focused gate lives in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R`
  and now freezes the current contracts for
  `buildMetabNormCorrelationFilterSummary()`
  and
  `resolveMetabNormFinalQcRenderState()`
  and
  `buildMetabNormFinalQcPcaPlot()`
  and
  `renderMetabNormFinalQcPlot()`
  and
  `appendMetabNormNormalizationLog()`
  and
  `renderMetabNormNormalizationLog()`
  and
  `resolveMetabNormExportSourceDir()`
  and
  `collectMetabNormFeatureCountsPerAssay()`
  and
  `buildMetabNormExportSessionData()`
  and
  `saveMetabNormExportSessionRdsFiles()`
  and
  `saveMetabNormExportMetadataFiles()`
  and
  `saveMetabNormExportSummaryFile()`
  and
  `runMetabNormExportSessionWorkflow()`
  and
  `checkMetabNormExportSessionReady()`
  and
  `dispatchMetabNormExportSession()`
  and
  `handleMetabNormExportSessionOutcome()`
  and
  `runMetabNormExportSessionObserverShell()`
  and
  `runMetabNormNormalizationObserverShell()`
  and
  `renderMetabNormItsdSelectionUi()`
  and
  `renderMetabNormItsdSelectionTable()`
  and
  `registerMetabNormItsdSelectionTracking()`
  and
  `runMetabNormItsdSelectionTableObserverShell()`
  and
  `renderMetabNormRuvQcUi()`
  and
  `renderMetabNormRuvCancorPlot()`
  and
  `renderMetabNormRuvOptimizationSummary()`
  and
  `renderMetabNormRuvResultsTable()`
  and
  `resolveMetabNormManualItsdFeatureIds()`
  and
  `runMetabNormPreNormalizationQcStep()`
  and
  `generateMetabNormPreNormalizationQc()`
  and
  `runMetabNormItsdProgressApplyShell()`
  and
  `runMetabNormItsdNormalizationStep()`
  and
  `runMetabNormLog2ProgressApplyShell()`
  and
  `runMetabNormLog2TransformationStep()`
  and
  `runMetabNormBetweenSampleProgressApplyShell()`
  and
  `runMetabNormBetweenSampleNormalizationStep()`
  and
  `runMetabNormRuvOptimizationStep()`
  and
  `runMetabNormRuvCorrectionStep()`
  and
  `runMetabNormRuvQcStep()`
  and
  `runMetabNormCompositeQcFigureStep()`
  and
  `runMetabNormCompositeQcRefreshShell()`
  and
  `runMetabNormNormalizationPipelineShell()`
  and
  `runMetabNormResetNormalizationObserverShell()`
  and
  `resolveMetabNormSkipCorrelationInputObject()`
  and
  `completeMetabNormSkipCorrelationState()`
  and
  `handleMetabNormSkipCorrelationOutcome()`
  and
  `dispatchMetabNormSkipCorrelation()`
  and
  `runMetabNormSkipCorrelationObserverEntry()`
  and
  `runMetabNormSkipCorrelationObserverShell()`
  and
  `runMetabNormApplyCorrelationWorkflow()`
  and
  `dispatchMetabNormApplyCorrelation()`
  and
  `handleMetabNormApplyCorrelationOutcome()`
  and
  `runMetabNormApplyCorrelationObserverEntry()`
  and
  `updateMetabNormDesignDrivenChoices()`
  and
  `initializeMetabNormAssayNames()`
  and
  `runMetabNormApplyCorrelationObserverShell()`
  and
  `runMetabNormAssayLabelBindingShell()`
  and
  `runMetabNormQcImageBindingShell()`
  and
  `runMetabNormItsdSelectionTrackingObserverShell()`.
- The newest bounded live seam now sits in
  `R/mod_metab_norm.R`
  at line `991`
  as `renderMetabNormRuvResultsTable()`, which owns the per-assay RUV results
  table render path's successful optimization-results lookup, datatable
  handoff, and null exit before control returns to `mod_metab_norm_server()`.
- The per-assay RUV results table binding in
  `R/mod_metab_norm.R`
  at line `1360`
  now delegates through that helper instead of keeping the
  `DT::renderDataTable` body inline inside the RUV render observer.
- The focused metabolomics normalization module helper gate now also freezes
  the per-assay RUV results table render seam contract in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R`
  at line `1454`
  and reran green again via direct `testthat::test_file()` because this
  worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `1432` lines with `14` top-level functions after the per-assay RUV
  results table render seam in `R/mod_metab_norm.R`.
- The newest bounded live seam now sits in
  `R/mod_metab_norm.R`
  at line `875`
  as `renderMetabNormItsdSelectionTable()`, which owns the per-assay ITSD
  selection table render path's assay lookup, `buildItsdSelectionTable()`
  handoff, candidate preselection, datatable formatting, and null exit before
  control returns to `mod_metab_norm_server()`.
- The per-assay ITSD selection table binding in
  `R/mod_metab_norm.R`
  at line `1274`
  now delegates through that helper instead of keeping the
  `DT::renderDataTable` body inline inside the ITSD render observer.
- The focused metabolomics normalization module helper gate now also freezes
  the per-assay ITSD selection table render seam contract in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R`
  at line `1247`
  and reran green again via direct `testthat::test_file()` because this
  worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `1456` lines with `15` top-level functions after the per-assay ITSD
  selection table render seam in `R/mod_metab_norm.R`.
- The newest bounded live seam now sits in
  `R/mod_metab_norm.R`
  at line `928`
  as `registerMetabNormItsdSelectionTracking()`, which owns the per-assay
  ITSD selection tracking path's observer registration, normalized input-id
  construction, selection-state updates, and logging before control returns to
  `mod_metab_norm_server()`.
- The ITSD selection tracking observe block in
  `R/mod_metab_norm.R`
  at line `1312`
  now delegates through that helper instead of keeping the
  `purrr::walk(norm_data$assay_names, ...)`
  block and nested `shiny::observeEvent()` registration inline inside the DT
  selection observer.
- The focused metabolomics normalization module helper gate now also freezes
  the per-assay ITSD selection tracking seam contract in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R`
  at lines `1390` and `1453`
  and reran green again via direct `testthat::test_file()` because this
  worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `1476` lines with `17` top-level functions after the per-assay ITSD
  selection tracking seam in `R/mod_metab_norm.R`.
- The newest bounded live seam now sits in
  `R/mod_metab_norm.R`
  at line `1133`
  as `runMetabNormRuvBindingObserverShell()`, which owns the per-assay RUV
  plot/summary/table binding path's assay-name `req()` guard, non-empty
  optimization-results `req()` guard, safe output-id derivation, and
  delegation into `renderMetabNormRuvCancorPlot()`,
  `renderMetabNormRuvOptimizationSummary()`, and
  `renderMetabNormRuvResultsTable()` before control returns to
  `mod_metab_norm_server()`.
- The per-assay RUV binding observe block in
  `R/mod_metab_norm.R`
  at line `1443`
  now delegates through that helper instead of keeping the
  `shiny::req(length(norm_data$ruv_optimization_results) > 0)`
  guard and
  `purrr::walk(norm_data$assay_names, ...)`
  plot/summary/table output binding block inline inside
  `mod_metab_norm_server()`.
- The focused metabolomics normalization module helper gate now also freezes
  the per-assay RUV binding observer-shell seam contract in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R`
  at lines `1913` and `1981`
  and reran green again via direct `testthat::test_file()` because this
  worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `1515` lines with `22` top-level functions after the per-assay RUV
  binding observer-shell seam in `R/mod_metab_norm.R`.
- Next safe stop point:
  introduce one bounded top-level seam for the static assay-label output
  binding cluster, beginning at `R/mod_metab_norm.R:1338` around the eight
  `renderMetabNormAssayLabel()` registrations, then rerun the focused gate.
- The newest bounded live seam now sits in
  `R/mod_metab_norm.R`
  at line `963`
  as `renderMetabNormRuvOptimizationSummary()`, which owns the per-assay RUV
  optimization summary render path's success formatting, failed-run message
  fallback, and pending-message fallback before control returns to
  `mod_metab_norm_server()`.
- The per-assay RUV summary binding in
  `R/mod_metab_norm.R`
  at line `1330`
  now delegates through that helper instead of keeping the `renderText` body
  inline inside the RUV render observer.
- The focused metabolomics normalization module helper gate now also freezes
  the per-assay RUV optimization summary render seam contract in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R`
  at line `1401`
  and reran green again via direct `testthat::test_file()` because this
  worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `1416` lines with `13` top-level functions after the per-assay RUV
  optimization summary render seam in `R/mod_metab_norm.R`.
- The first bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `buildMetabNormCorrelationFilterSummary()`, which now owns the
  `correlation_filter_summary` render path's empty-result fallback, per-assay
  Pearson summary formatting, and sample-count delta summary before the wrapper
  falls through to the remaining final-QC and export paths.
- The `correlation_filter_summary` render binding in
  `R/mod_metab_norm.R`
  now delegates through that helper instead of keeping the formatting block
  inline inside `mod_metab_norm_server()` at line `2175`.
- The second bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `resolveMetabNormFinalQcRenderState()`, which now owns the
  `final_qc_plot` render path's source-object priority order and reusable
  empty-state fallback plot before the wrapper falls through to PCA rendering.
- The `final_qc_plot` render binding in
  `R/mod_metab_norm.R`
  now delegates through that helper for source resolution and no-data fallback
  instead of keeping that decision block inline inside
  `mod_metab_norm_server()` at line `2188`.
- The third bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `buildMetabNormFinalQcPcaPlot()`, which now owns the
  `final_qc_plot` render path's `plotPca()` argument forwarding, multi-plot
  return-shape normalization, and error fallback plot before the wrapper falls
  through to the export observer.
- The `final_qc_plot` render binding in
  `R/mod_metab_norm.R`
  now delegates through that helper for PCA rendering and error handling
  instead of keeping that normalization block inline inside
  `mod_metab_norm_server()` at line `2201`.
- The fourth bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `resolveMetabNormExportSourceDir()`, which now owns the export-session
  observer's source-directory priority order, `export_dir` fallback, and
  missing fallback-directory creation before the wrapper falls through to the
  larger session payload assembly.
- The `input$export_session` observer binding in
  `R/mod_metab_norm.R`
  now delegates through that helper for export-directory resolution instead of
  keeping that fallback block inline inside `mod_metab_norm_server()` at line
  `2226`.
- The fifth bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `collectMetabNormFeatureCountsPerAssay()`, which now owns the
  export-session observer's per-assay feature and sample counting from the
  current `MetaboliteAssayData` object before the wrapper falls through to the
  larger session payload assembly.
- The export-session session-data helper in
  `R/mod_metab_norm.R`
  at line `815`
  now delegates through that helper for per-assay feature-count collection
  instead of keeping that loop inline inside `mod_metab_norm_server()`.
- The sixth bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `buildMetabNormExportSessionData()`, which now owns the export-session
  observer's current-state lookup, payload field assembly, and reuse of
  `collectMetabNormFeatureCountsPerAssay()` before the wrapper falls through to
  the downstream file writes.
- The export-session orchestration helper in
  `R/mod_metab_norm.R`
  at line `1019`
  now delegates through that helper for session-data payload assembly instead
  of keeping that list construction inline inside `mod_metab_norm_server()`.
- The seventh bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `saveMetabNormExportSessionRdsFiles()`, which now owns the
  export-session observer's timestamped main-session filename generation,
  main/latest session RDS writes, progress-step forwarding, and save logging
  before the wrapper falls through to the downstream metadata and summary
  writers.
- The export-session orchestration helper in
  `R/mod_metab_norm.R`
  at line `1033`
  now delegates through that helper for main/latest session RDS writing
  instead of keeping that block inline inside `mod_metab_norm_server()`.
- The eighth bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `saveMetabNormExportMetadataFiles()`, which now owns the export-session
  observer's metadata-file redundancy writes for per-assay RUV optimization
  results, ITSD selections, QC parameters, and warning-only fallback handling
  before the wrapper falls through to the downstream summary writer.
- The export-session orchestration helper in
  `R/mod_metab_norm.R`
  at line `1040`
  now delegates through that helper for metadata-file redundancy writes
  instead of keeping that try-catch block inline inside
  `mod_metab_norm_server()`.
- The ninth bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `saveMetabNormExportSummaryFile()`, which now owns the export-session
  observer's human-readable summary assembly, successful per-assay RUV summary
  filtering, feature-count summary formatting, summary-file write, and summary
  logging before the wrapper falls through to the remaining final
  notifications.
- The export-session orchestration helper in
  `R/mod_metab_norm.R`
  at line `1046`
  now delegates through that helper for human-readable summary-file creation
  instead of keeping that block inline inside `mod_metab_norm_server()`.
- The tenth bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `runMetabNormExportSessionWorkflow()`, which now owns the export-session
  observer's `withProgress()` shell, progress-step orchestration, session-data
  export logging, and delegation to the existing RDS, metadata, and summary
  writer seams before the wrapper falls through to the remaining completion
  notifications and error handling.
- The `input$export_session` observer binding in
  `R/mod_metab_norm.R`
  now delegates through that helper for export-session orchestration instead
  of keeping that `withProgress()` block inline inside
  `mod_metab_norm_server()` at line `2228`.
- The eleventh bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `handleMetabNormExportSessionOutcome()`, which now owns the
  export-session observer's success log, success notification, completion log,
  and mirrored error log/notification tail after the export workflow helper
  returns or throws.
- The `input$export_session` observer binding in
  `R/mod_metab_norm.R`
  now delegates through that helper for final export outcome reporting instead
  of keeping that success/error notification tail inline inside
  `mod_metab_norm_server()` at line `2283`.
- The twelfth bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `checkMetabNormExportSessionReady()`, which now owns the
  export-session observer's normalization-complete warning notification gate
  before the wrapper falls through to the export-session dispatch helper seam.
- The `input$export_session` observer binding in
  `R/mod_metab_norm.R`
  now delegates through that helper for normalization-complete gating instead
  of keeping that warning notification inline inside
  `mod_metab_norm_server()` at line `2303`.
- The thirteenth bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `dispatchMetabNormExportSession()`, which now owns the
  export-session observer's remaining `tryCatch` dispatch shell, source-dir
  resolution handoff, workflow helper call, and outcome helper delegation
  before the wrapper falls through to the next skip-correlation observer seam.
- The `input$export_session` observer binding in
  `R/mod_metab_norm.R`
  now delegates through that helper for export-session dispatch instead of
  keeping that `tryCatch` shell inline inside
  `mod_metab_norm_server()` at line `2314`.
- The fourteenth bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `completeMetabNormSkipCorrelationState()`, which now owns the
  skip-correlation-filter observer's normalization-complete `saveState()`
  call and tab-status completion update before the wrapper falls through to
  the remaining success log and notification tail.
- The `input$skip_correlation_filter` observer binding in
  `R/mod_metab_norm.R`
  now delegates through that helper for the normalization-complete
  state-save/status update shell instead of keeping that block inline inside
  `mod_metab_norm_server()` at line `2275`.
- The fifteenth bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `handleMetabNormSkipCorrelationOutcome()`, which now owns the
  skip-correlation-filter observer's success log entry and ready-for-DE
  notification tail after the state-save/status helper returns.
- The `input$skip_correlation_filter` observer binding in
  `R/mod_metab_norm.R`
  now delegates through that helper for the success feedback tail instead of
  keeping that add-log/show-notification shell inline inside
  `mod_metab_norm_server()` at line `2280`.
- The sixteenth bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `resolveMetabNormSkipCorrelationInputObject()`, which now owns the
  skip-correlation-filter observer's current normalized-object priority order
  and `NULL` fallback before the wrapper falls through to the state-save/status
  and success-feedback helpers.
- The `input$skip_correlation_filter` observer binding in
  `R/mod_metab_norm.R`
  now delegates through that helper for current normalized-object resolution
  instead of keeping that selection block inline inside
  `mod_metab_norm_server()` at line `2266`.
- The seventeenth bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `dispatchMetabNormSkipCorrelation()`, which now owns the
  skip-correlation-filter observer's remaining `NULL` short-circuit,
  state-save/status helper delegation, and success-feedback helper delegation
  before the wrapper falls through to the observer-entry seam.
- The observer-entry helper in
  `R/mod_metab_norm.R`
  now delegates through that helper for skip-correlation dispatch instead of
  keeping that remaining shell inline before the wrapper returns to
  `mod_metab_norm_server()`.
- The eighteenth bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `runMetabNormSkipCorrelationObserverEntry()`, which now owns the
  skip-correlation-filter observer's remaining current normalized-object
  resolution and dispatch handoff before the wrapper falls through to the
  remaining apply-correlation-filter observer extraction seam.
- The `input$skip_correlation_filter` observer binding in
  `R/mod_metab_norm.R`
  now delegates through that helper instead of keeping that entry shell inline
  inside `mod_metab_norm_server()` at line `2319`.
- The nineteenth bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `runMetabNormApplyCorrelationObserverEntry()`, which now owns the
  apply-correlation-filter observer's current normalized-object resolution,
  threshold/grouping setup capture, working-notification bootstrap, and
  observer-state construction before the wrapper falls through to the dispatch
  handoff seam.
- The `input$apply_correlation_filter` observer binding in
  `R/mod_metab_norm.R`
  now delegates through that helper instead of keeping that entry shell inline
  inside `mod_metab_norm_server()` at line `2258`.
- The twentieth bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `dispatchMetabNormApplyCorrelation()`, which now owns the
  apply-correlation-filter observer's remaining `tryCatch` dispatch shell,
  apply-correlation workflow-helper delegation, and outcome-helper delegation
  for the wrapper before the observer-entry helper returns.
- The observer-entry helper in
  `R/mod_metab_norm.R`
  now delegates through that helper instead of keeping that `tryCatch` shell
  inline before control returns to `mod_metab_norm_server()`.
- The twenty-first bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `runMetabNormApplyCorrelationObserverEntry()`, which now also owns the
  workflow-data/norm-data forwarding, remove-notification forwarding, and the
  final handoff into `dispatchMetabNormApplyCorrelation()` so the server
  observer is reduced to a single top-level helper call.
- The `input$apply_correlation_filter` observer binding in
  `R/mod_metab_norm.R`
  now delegates fully through that helper instead of keeping the last
  dispatch-handoff shell inline inside `mod_metab_norm_server()` at line
  `2258`.
- The twenty-second bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `handleMetabNormApplyCorrelationOutcome()`, which now owns the
  apply-correlation-filter observer's correlation-filtered state persistence,
  tab-status completion update, success log/notification cleanup, and mirrored
  error cleanup after the apply-correlation workflow helper returns or throws.
- The dispatch helper in
  `R/mod_metab_norm.R`
  at line `2386`
  now delegates through that helper for the completion/outcome tail instead of
  keeping that success/error cleanup inline before control returns to the
  observer-entry helper.
- The twenty-third bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `runMetabNormApplyCorrelationWorkflow()`, which now owns the
  apply-correlation-filter observer's remaining `req()` gate, Pearson
  correlation calculation, and correlation-threshold filter shell before the
  dispatch helper falls through to the outcome helper.
- The dispatch helper in
  `R/mod_metab_norm.R`
  at line `2386`
  now delegates through that helper for the remaining correlation
  calculation/filter shell instead of keeping those calls inline inside
  `dispatchMetabNormApplyCorrelation()`.
- Because this worktree does not include `renv/activate.R`,
  `tools/test_with_renv.R`
  is not available for this target; the focused gate reran green via direct
  `testthat::test_file()`.
- Classification refresh on April 17, 2026 still reports
  `R/mod_metab_norm.R`
  as `high-risk-wrapper` / `needs-seam-introduction` after the staged
  apply-correlation wave and live apply.
- The first exact-source metabolomics normalization module wave now lives in
  `tools/refactor/manifest-metab-norm-module-wave1.yml`
  and its reviewed staging artifacts remain in
  `tools/refactor/staging/wave1_metabolomics_norm_module_apply_correlation_helpers/R/mod_metab_norm_server_helpers.R`.
- The applied wave now lives in
  `R/mod_metab_norm_server_helpers.R`
  and owns
  `runMetabNormApplyCorrelationWorkflow()`,
  `dispatchMetabNormApplyCorrelation()`,
  `handleMetabNormApplyCorrelationOutcome()`,
  and
  `runMetabNormApplyCorrelationObserverEntry()` after exact-source removal from
  `R/mod_metab_norm.R`.
- `tools/refactor/verify_refactor.R` passed for
  `tools/refactor/manifest-metab-norm-module-wave1.yml`
  before the live apply, and `tools/refactor/check_wave_apply.R` passed after
  the live apply.
- `DESCRIPTION` `Collate:` now loads
  `mod_metab_norm_server_helpers.R`
  before
  `mod_metab_norm.R`.
- The focused gate loader in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R`
  now widens the future apply boundary by loading
  `R/mod_metab_norm_server_helpers.R`
  before
  `R/mod_metab_norm.R`.
- The twenty-fourth bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `runMetabNormApplyCorrelationObserverShell()`, which now owns the
  remaining `input$apply_correlation_filter` observer's
  `workflow_data$state_manager` gate, normalization-or-RUV readiness gate, and
  handoff into `runMetabNormApplyCorrelationObserverEntry()` before the wrapper
  falls through to the mirrored skip-correlation observer shell.
- The `input$apply_correlation_filter` observer binding in
  `R/mod_metab_norm.R`
  now delegates through that helper instead of keeping that remaining shell
  inline inside `mod_metab_norm_server()` at line `2332`.
- The focused gate in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R`
  now also freezes the apply-correlation observer shell contract for
  `runMetabNormApplyCorrelationObserverShell()`.
- The twenty-fifth bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `runMetabNormSkipCorrelationObserverShell()`, which now owns the
  remaining `input$skip_correlation_filter` observer's
  `workflow_data$state_manager` gate, normalization-or-RUV readiness gate, and
  handoff into `runMetabNormSkipCorrelationObserverEntry()` before the wrapper
  falls through to the remaining export-session observer shell.
- The `input$skip_correlation_filter` observer binding in
  `R/mod_metab_norm.R`
  now delegates through that helper instead of keeping that remaining shell
  inline inside `mod_metab_norm_server()` at line `2349`.
- The focused gate in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R`
  now also freezes the skip-correlation observer shell contract for
  `runMetabNormSkipCorrelationObserverShell()`.
- The twenty-sixth bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `runMetabNormExportSessionObserverShell()`, which now owns the remaining
  `input$export_session` observer's export-button log line,
  `workflow_data$state_manager` gate, normalization-complete readiness gate,
  and handoff into `dispatchMetabNormExportSession()` before the wrapper falls
  through to the next reset observer shell target.
- The `input$export_session` observer binding in
  `R/mod_metab_norm.R`
  now delegates through that helper instead of keeping that remaining shell
  inline inside `mod_metab_norm_server()` at line `2412`.
- The focused gate in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R`
  now also freezes the export-session observer shell contract for
  `runMetabNormExportSessionObserverShell()`.
- The focused metabolomics normalization module helper gate reran green after
  the live apply via direct `testthat::test_file()` because this worktree
  does not include `renv/activate.R` for `tools/test_with_renv.R`.
- The twenty-seventh bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `runMetabNormResetNormalizationObserverShell()`, which now owns the
  remaining `input$reset_normalization` observer's
  `workflow_data$state_manager` gate, filtered-state save shell, local reset
  state cleanup, and mirrored error-notification tail before the wrapper falls
  through to the remaining render/observe blocks.
- The `input$reset_normalization` observer binding in
  `R/mod_metab_norm.R`
  now delegates through that helper instead of keeping that remaining reset
  shell inline inside `mod_metab_norm_server()` at line `2291`.
- The focused gate in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R`
  now also freezes the reset observer shell's success and error contracts for
  `runMetabNormResetNormalizationObserverShell()`.
- The focused metabolomics normalization module helper gate reran green again
  after the reset observer shell seam via direct `testthat::test_file()`
  because this worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- The second exact-source metabolomics normalization module wave now lives in
  `tools/refactor/manifest-metab-norm-module-wave2.yml`.
- The live apply now extends
  `R/mod_metab_norm_server_helpers.R`
  and owns
  `checkMetabNormExportSessionReady()`,
  `dispatchMetabNormExportSession()`,
  `handleMetabNormExportSessionOutcome()`,
  `runMetabNormExportSessionObserverShell()`,
  `runMetabNormResetNormalizationObserverShell()`,
  `resolveMetabNormSkipCorrelationInputObject()`,
  `completeMetabNormSkipCorrelationState()`,
  `handleMetabNormSkipCorrelationOutcome()`,
  `dispatchMetabNormSkipCorrelation()`,
  `runMetabNormSkipCorrelationObserverEntry()`,
  and
  `runMetabNormSkipCorrelationObserverShell()`
  after exact-source removal from
  `R/mod_metab_norm.R`.
- `tools/refactor/verify_refactor.R` passed for
  `tools/refactor/manifest-metab-norm-module-wave2.yml`
  before the live apply, and `tools/refactor/check_wave_apply.R` passed after
  the live apply.
- `DESCRIPTION` `Collate:` already loaded
  `mod_metab_norm_server_helpers.R`
  before
  `mod_metab_norm.R`,
  so no additional load-order update was required for the wave-2 apply.
- The focused metabolomics normalization module helper gate reran green after
  the wave-2 apply via direct `testthat::test_file()` because this worktree
  still does not include `renv/activate.R` for `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `2175` lines with `14` top-level functions after the wave-2
  export/reset/skip observer-shell apply in
  `R/mod_metab_norm.R`.
- The third exact-source metabolomics normalization module wave now lives in
  `tools/refactor/manifest-metab-norm-module-wave3.yml`,
  with the reviewed staging artifacts in
  `tools/refactor/staging/wave3_metabolomics_norm_module_export_workflow_helpers/R/mod_metab_norm_server_helpers.R`
  and
  `tools/refactor/staging/wave3_metabolomics_norm_module_export_workflow_helpers/collate-metab-norm-module-wave3.txt`.
- The live apply now further extends
  `R/mod_metab_norm_server_helpers.R`
  and owns
  `resolveMetabNormExportSourceDir()`,
  `collectMetabNormFeatureCountsPerAssay()`,
  `buildMetabNormExportSessionData()`,
  `saveMetabNormExportSessionRdsFiles()`,
  `saveMetabNormExportMetadataFiles()`,
  `saveMetabNormExportSummaryFile()`,
  and
  `runMetabNormExportSessionWorkflow()`
  after exact-source removal from
  `R/mod_metab_norm.R`.
- `tools/refactor/verify_refactor.R` passed for
  `tools/refactor/manifest-metab-norm-module-wave3.yml`
  before the live apply, and `tools/refactor/check_wave_apply.R` passed after
  the live apply.
- `DESCRIPTION` `Collate:` already loaded
  `mod_metab_norm_server_helpers.R`
  before
  `mod_metab_norm.R`,
  so no additional load-order update was required for the wave-3 apply.
- The focused metabolomics normalization module helper gate reran green after
  the wave-3 apply via direct `testthat::test_file()` because this worktree
  still does not include `renv/activate.R` for `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `1901` lines with `7` top-level functions after the wave-3
  export-workflow helper apply in
  `R/mod_metab_norm.R`.
- Next safe stop point:
  stage the next exact-source metabolomics normalization module helper wave
  for the remaining final-QC render helper trio into
  `R/mod_metab_norm_server_helpers.R`,
  then rerun the focused gate.
- The fourth exact-source metabolomics normalization module wave now lives in
  `tools/refactor/manifest-metab-norm-module-wave4.yml`,
  with the reviewed staging artifacts in
  `tools/refactor/staging/wave4_metabolomics_norm_module_final_qc_render_helpers/R/mod_metab_norm_server_helpers.R`
  and
  `tools/refactor/staging/wave4_metabolomics_norm_module_final_qc_render_helpers/collate-metab-norm-module-wave4.txt`.
- The live apply now further extends
  `R/mod_metab_norm_server_helpers.R`
  and owns
  `buildMetabNormCorrelationFilterSummary()`,
  `resolveMetabNormFinalQcRenderState()`,
  and
  `buildMetabNormFinalQcPcaPlot()`
  after exact-source removal from
  `R/mod_metab_norm.R`.
- `tools/refactor/verify_refactor.R` passed for
  `tools/refactor/manifest-metab-norm-module-wave4.yml`
  before the live apply, and `tools/refactor/check_wave_apply.R` passed after
  the live apply.
- `DESCRIPTION` `Collate:` already loaded
  `mod_metab_norm_server_helpers.R`
  before
  `mod_metab_norm.R`,
  so no additional load-order update was required for the wave-4 apply.
- The focused metabolomics normalization module helper gate reran green after
  the wave-4 apply via direct `testthat::test_file()` because this worktree
  still does not include `renv/activate.R` for `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `1809` lines with `4` top-level functions after the wave-4
  final-QC render helper apply in
  `R/mod_metab_norm.R`.
- The fifth exact-source metabolomics normalization module wave now lives in
  `tools/refactor/manifest-metab-norm-module-wave5.yml`,
  with the reviewed staging artifacts in
  `tools/refactor/staging/wave5_metabolomics_norm_module_apply_correlation_observer_shell/R/mod_metab_norm_server_helpers.R`
  and
  `tools/refactor/staging/wave5_metabolomics_norm_module_apply_correlation_observer_shell/collate-metab-norm-module-wave5.txt`.
- `tools/refactor/verify_refactor.R` passed for
  `tools/refactor/manifest-metab-norm-module-wave5.yml`
  before the staging run.
- The staged helper file cleanly extracts
  `runMetabNormApplyCorrelationObserverShell()`
  from
  `R/mod_metab_norm.R`
  into
  `R/mod_metab_norm_server_helpers.R`
  without mutating live package sources.
- The focused metabolomics normalization module helper gate reran green again
  after the wave-5 staging run via direct `testthat::test_file()` because this
  worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `1809` lines with `4` top-level functions after the wave-5 staging-only
  checkpoint in
  `R/mod_metab_norm.R`.
- The fifth exact-source metabolomics normalization module wave is now live via
  `tools/refactor/manifest-metab-norm-module-wave5.yml`.
- The live apply now further extends
  `R/mod_metab_norm_server_helpers.R`
  and owns
  `runMetabNormApplyCorrelationObserverShell()`
  after exact-source removal from
  `R/mod_metab_norm.R`.
- `tools/refactor/verify_refactor.R` passed again for
  `tools/refactor/manifest-metab-norm-module-wave5.yml`
  before the live apply, and `tools/refactor/check_wave_apply.R` passed after
  the live apply.
- `DESCRIPTION` `Collate:` already loaded
  `mod_metab_norm_server_helpers.R`
  before
  `mod_metab_norm.R`,
  so no additional load-order update was required for the wave-5 apply.
- The focused metabolomics normalization module helper gate reran green again
  after the wave-5 apply via direct `testthat::test_file()` because this
  worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `1783` lines with `2` top-level functions after the wave-5
  apply-correlation observer-shell apply in
  `R/mod_metab_norm.R`.
- The focused gate now also freezes the run-normalization observer-shell
  contract in
  `R/mod_metab_norm.R`
  as `runMetabNormNormalizationObserverShell()`.
- The twenty-eighth bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `runMetabNormNormalizationObserverShell()`, which now owns the remaining
  `input$run_normalization` observer's `workflow_data$state_manager` gate,
  top-level `withProgress()` shell, and mirrored completion/error
  notification tail while the inline callback still owns the detailed
  normalization pipeline body.
- The `input$run_normalization` observer binding in
  `R/mod_metab_norm.R`
  now delegates through that helper instead of keeping that remaining shell
  inline inside `mod_metab_norm_server()` at line `1364`.
- The focused metabolomics normalization module helper gate reran green again
  after the live run-normalization observer-shell seam via direct
  `testthat::test_file()` because this worktree still does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- The sixth exact-source metabolomics normalization module wave now lives in
  `tools/refactor/manifest-metab-norm-module-wave6.yml`,
  with the reviewed staging artifacts in
  `tools/refactor/staging/wave6_metabolomics_norm_module_run_normalization_observer_shell/R/mod_metab_norm_server_helpers.R`
  and
  `tools/refactor/staging/wave6_metabolomics_norm_module_run_normalization_observer_shell/collate-metab-norm-module-wave6.txt`.
- `tools/refactor/check_wave_apply.R` passed for the live apply of
  `tools/refactor/manifest-metab-norm-module-wave6.yml`.
- The reviewed wave-6 apply moved
  `runMetabNormNormalizationObserverShell()`
  from
  `R/mod_metab_norm.R`
  into
  `R/mod_metab_norm_server_helpers.R`
  while keeping the inline normalization pipeline callback inside
  `mod_metab_norm_server()`.
- The focused gate now also freezes the run-normalization observer-shell
  contract in
  `R/mod_metab_norm_server_helpers.R`
  as `runMetabNormNormalizationObserverShell()`.
- The `input$run_normalization` observer binding in
  `R/mod_metab_norm.R`
  still delegates through that helper while the inline callback owns the
  detailed normalization pipeline body.
- The focused metabolomics normalization module helper gate reran green again
  after the wave-6 live apply via direct `testthat::test_file()` because this
  worktree still does not include `renv/activate.R` for `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `1778` lines with `2` top-level functions after the wave-6 live apply in
  `R/mod_metab_norm.R`.
- The twenty-ninth bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `resolveMetabNormManualItsdFeatureIds()`, which now owns the manual ITSD
  feature-ID resolution block while the inline callback still owns the
  surrounding ITSD normalization step.
- The inline `runPipelineFn` callback passed to
  `runMetabNormNormalizationObserverShell()`
  now delegates through that helper instead of rebuilding manual feature IDs
  inline inside `mod_metab_norm_server()`.
- The focused gate now also freezes the manual ITSD feature-ID seam contract
  in
  `R/mod_metab_norm.R`
  as `resolveMetabNormManualItsdFeatureIds()`.
- The focused metabolomics normalization module helper gate reran green again
  after the manual ITSD feature-ID seam via direct `testthat::test_file()`
  because this worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Classification refresh on April 17, 2026 still reports
  `R/mod_metab_norm.R`
  as `high-risk-wrapper` / `needs-seam-introduction` at `1778` lines with `4`
  top-level functions after the manual ITSD feature-ID seam checkpoint.
- The thirtieth bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `runMetabNormItsdNormalizationStep()`, which now owns the ITSD
  normalization apply/saveState block while the inline callback still owns the
  surrounding normalization pipeline orchestration.
- The inline `runPipelineFn` callback passed to
  `runMetabNormNormalizationObserverShell()`
  now delegates through that helper instead of keeping the ITSD apply/saveState
  block inline inside `mod_metab_norm_server()` at line `1420`.
- The focused gate now also freezes the ITSD normalization apply/saveState seam
  contract in
  `R/mod_metab_norm.R`
  as `runMetabNormItsdNormalizationStep()`.
- The focused metabolomics normalization module helper gate reran green again
  after the ITSD normalization apply/saveState seam via direct
  `testthat::test_file()` because this worktree still does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- Classification refresh on April 17, 2026 now reports
  `R/mod_metab_norm.R`
  as `high-risk-wrapper` / `needs-seam-introduction` at `1811` lines with `6`
  top-level functions after the ITSD normalization apply/saveState seam
  checkpoint.
- The thirty-first bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `runMetabNormLog2TransformationStep()`, which now owns the log2
  transformation apply/saveState block while the inline callback still owns
  the surrounding normalization pipeline orchestration.
- The inline `runPipelineFn` callback passed to
  `runMetabNormNormalizationObserverShell()`
  now delegates through that helper instead of keeping the log2
  transformation apply/saveState block inline inside `mod_metab_norm_server()`
  at line `1473`.
- The focused gate now also freezes the log2 transformation apply/saveState
  seam contract in
  `R/mod_metab_norm.R`
  as `runMetabNormLog2TransformationStep()`.
- The focused metabolomics normalization module helper gate reran green again
  after the log2 transformation apply/saveState seam via direct
  `testthat::test_file()` because this worktree still does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- Classification refresh on April 17, 2026 now reports
  `R/mod_metab_norm.R`
  as `high-risk-wrapper` / `needs-seam-introduction` at `1842` lines with `8`
  top-level functions after the log2 transformation apply/saveState seam
  checkpoint.
- The thirty-second bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `runMetabNormBetweenSampleNormalizationStep()`, which now owns the
  between-sample normalization apply/saveState block while the inline callback
  still owns the surrounding normalization pipeline orchestration.
- The inline `runPipelineFn` callback passed to
  `runMetabNormNormalizationObserverShell()`
  now delegates through that helper instead of keeping the between-sample
  normalization apply/saveState block inline inside `mod_metab_norm_server()`
  at line `1528`.
- The focused gate now also freezes the between-sample normalization
  apply/saveState seam contract in
  `R/mod_metab_norm.R`
  as `runMetabNormBetweenSampleNormalizationStep()`.
- The focused metabolomics normalization module helper gate reran green again
  after the between-sample normalization apply/saveState seam via direct
  `testthat::test_file()` because this worktree still does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- Classification refresh on April 17, 2026 now reports
  `R/mod_metab_norm.R`
  as `high-risk-wrapper` / `needs-seam-introduction` at `1876` lines with `10`
  top-level functions after the between-sample normalization apply/saveState
  seam checkpoint.
- The thirty-third bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `runMetabNormPostNormalizationQcStep()`, which now owns the
  post-normalization QC generation block while the inline callback still owns
  the surrounding normalization pipeline orchestration.
- The inline `runPipelineFn` callback passed to
  `runMetabNormNormalizationObserverShell()`
  now delegates through that helper instead of keeping the post-normalization
  QC generation block inline inside `mod_metab_norm_server()` at line `1566`.
- The focused gate now also freezes the post-normalization QC seam contract in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R`
  as `runMetabNormPostNormalizationQcStep()`.
- The focused metabolomics normalization module helper gate reran green again
  after the post-normalization QC seam via direct `testthat::test_file()`
  because this worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Classification refresh on April 17, 2026 now reports
  `R/mod_metab_norm.R`
  as `high-risk-wrapper` / `needs-seam-introduction` at `1904` lines with `12`
  top-level functions after the post-normalization QC seam checkpoint.
- The thirty-fourth bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `runMetabNormRuvOptimizationStep()`, which now owns the opening RUV-III
  optimization block while the inline callback still owns the surrounding
  normalization pipeline orchestration.
- The inline `runPipelineFn` callback passed to
  `runMetabNormNormalizationObserverShell()`
  now delegates through that helper instead of keeping the `ruv_params` /
  `runPerAssayRuvOptimization()` block inline inside `mod_metab_norm_server()`
  at line `1638`.
- The focused gate now also freezes the opening RUV-III optimization seam
  contract in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R`
  as `runMetabNormRuvOptimizationStep()`.
- The focused metabolomics normalization module helper gate reran green again
  after the opening RUV-III optimization seam via direct
  `testthat::test_file()` because this worktree still does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- Classification refresh on April 17, 2026 now reports
  `R/mod_metab_norm.R`
  as `high-risk-wrapper` / `needs-seam-introduction` at `1953` lines with `14`
  top-level functions after the opening RUV-III optimization seam checkpoint.
- The thirty-fifth bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `runMetabNormRuvCorrectionStep()`, which now owns the RUV-III
  apply/saveState block while the inline callback still owns the downstream
  RUV QC plot generation block.
- The inline `runPipelineFn` callback passed to
  `runMetabNormNormalizationObserverShell()`
  now delegates through that helper instead of keeping the
  `ruvIII_C_Varying()` / `saveState()` block inline inside
  `mod_metab_norm_server()` at line `1697`.
- The focused gate now also freezes the RUV-III apply/saveState seam contract
  in `tests/testthat/test-metab-03b-norm-module-helper-characterization.R`
  as `runMetabNormRuvCorrectionStep()`.
- The focused metabolomics normalization module helper gate reran green again
  after the RUV-III apply/saveState seam via direct `testthat::test_file()`
  because this worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Classification refresh on April 17, 2026 now reports
  `R/mod_metab_norm.R`
  as `high-risk-wrapper` / `needs-seam-introduction` at `1987` lines with `16`
  top-level functions after the RUV-III apply/saveState seam checkpoint.
- The thirty-sixth bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `runMetabNormRuvQcStep()`, which now owns the RUV QC progress update,
  saved-plot generation call, and completion log entry while the inline
  callback still owns the downstream composite QC figure block.
- The inline `runPipelineFn` callback passed to
  `runMetabNormNormalizationObserverShell()`
  now delegates through that helper instead of keeping the
  `shiny::incProgress(..., detail = "Generating RUV QC plots...")` /
  `generateMetabQcPlots()` block inline inside `mod_metab_norm_server()` at
  line `1743`.
- The focused gate now also freezes the RUV QC generation seam contract in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R`
  as `runMetabNormRuvQcStep()`.
- The focused metabolomics normalization module helper gate reran green again
  after the RUV QC generation seam via direct `testthat::test_file()`
  because this worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Classification refresh on April 17, 2026 now reports
  `R/mod_metab_norm.R`
  as `high-risk-wrapper` / `needs-seam-introduction` at `2021` lines with `18`
  top-level functions after the RUV QC generation seam checkpoint.
- The thirty-seventh bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `runMetabNormCompositeQcFigureStep()`, which now owns the composite-QC
  opening log entry, per-assay saved-image file-list and row-label assembly,
  `generateCompositeFromFiles()` delegation, savePlot shell, and warning/log
  fallback while the inline callback keeps only the downstream plot-refresh
  tail.
- The inline `runPipelineFn` callback passed to
  `runMetabNormNormalizationObserverShell()`
  now delegates through that helper instead of keeping the
  `add_log("Generating composite QC figure...")` block inline inside
  `mod_metab_norm_server()` at line `1887`.
- The focused gate now also freezes the composite QC figure seam contract in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R`
  as `runMetabNormCompositeQcFigureStep()`.
- The focused metabolomics normalization module helper gate reran green again
  after the composite QC figure seam via direct `testthat::test_file()`
  because this worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Classification refresh on April 17, 2026 now reports
  `R/mod_metab_norm.R`
  as `high-risk-wrapper` / `needs-seam-introduction` at `2075` lines with `20`
  top-level functions after the composite QC figure seam checkpoint.
- The thirty-eighth bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `runMetabNormPreNormalizationQcStep()`, which now owns the opening
  pre-normalization progress shell, post-filter state capture, and
  pre-normalization QC plot generation while the inline callback keeps the
  downstream ITSD branch and later normalization steps.
- The inline `runPipelineFn` callback passed to
  `runMetabNormNormalizationObserverShell()`
  now delegates through that helper instead of keeping the
  `shiny::incProgress(..., detail = "Capturing pre-normalization state...")` /
  `generateMetabQcPlots()` block inline inside `mod_metab_norm_server()` at
  line `1799`.
- The focused gate now also freezes the pre-normalization capture/pre-QC seam
  contract in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R`
  as `runMetabNormPreNormalizationQcStep()`.
- The focused metabolomics normalization module helper gate reran green again
  after the pre-normalization capture/pre-QC seam via direct
  `testthat::test_file()`
  because this worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Classification refresh on April 17, 2026 now reports
  `R/mod_metab_norm.R`
  as `high-risk-wrapper` / `needs-seam-introduction` at `2112` lines with `22`
  top-level functions after the pre-normalization capture/pre-QC seam
  checkpoint.
- The thirty-ninth bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `runMetabNormItsdProgressApplyShell()`, which now owns the ITSD-step
  progress update, apply-branch log entry, manual ITSD feature-ID resolution
  handoff, skipped-branch log entry, and delegation to
  `runMetabNormItsdNormalizationStep()` while the inline callback keeps the
  downstream log2 branch and later normalization steps.
- The inline `runPipelineFn` callback passed to
  `runMetabNormNormalizationObserverShell()`
  now delegates through that helper instead of keeping the
  `shiny::incProgress(..., detail = "Applying ITSD normalization...")` /
  `if (isTRUE(input$apply_itsd))` block inline inside `mod_metab_norm_server()`
  at line `1873`.
- The focused gate now also freezes the ITSD progress/apply shell contract in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R`
  as `runMetabNormItsdProgressApplyShell()`.
- The focused metabolomics normalization module helper gate reran green again
  after the ITSD progress/apply shell seam via direct `testthat::test_file()`
  because this worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Classification refresh on April 17, 2026 now reports
  `R/mod_metab_norm.R`
  as `high-risk-wrapper` / `needs-seam-introduction` at `2163` lines with `24`
  top-level functions after the ITSD progress/apply shell checkpoint.
- The fortieth bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `runMetabNormLog2ProgressApplyShell()`, which now owns the log2-step
  progress update, offset log entry, and delegation to
  `runMetabNormLog2TransformationStep()` while the inline callback keeps the
  downstream between-sample normalization branch and later normalization
  steps.
- The inline `runPipelineFn` callback passed to
  `runMetabNormNormalizationObserverShell()`
  now delegates through that helper instead of keeping the
  `shiny::incProgress(..., detail = "Applying log2 transformation...")` /
  `add_log(paste("Applying log2 transformation..."))` block inline inside
  `mod_metab_norm_server()` at line `1920`.
- The focused gate now also freezes the log2 progress/apply shell contract in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R`
  as `runMetabNormLog2ProgressApplyShell()`.
- The focused metabolomics normalization module helper gate reran green again
  after the log2 progress/apply shell seam via direct `testthat::test_file()`
  because this worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Classification refresh on April 17, 2026 now reports
  `R/mod_metab_norm.R`
  as `high-risk-wrapper` / `needs-seam-introduction` at `2195` lines with `26`
  top-level functions after the log2 progress/apply shell checkpoint.
- The forty-first bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `runMetabNormBetweenSampleProgressApplyShell()`, which now owns the
  between-sample step progress update, conditional apply log entry, and
  delegation to `runMetabNormBetweenSampleNormalizationStep()` while the
  inline callback keeps the downstream post-normalization QC and RUV branches.
- The inline `runPipelineFn` callback passed to
  `runMetabNormNormalizationObserverShell()`
  now delegates through that helper instead of keeping the
  `shiny::incProgress(..., detail = "Applying between-sample normalization...")`
  / conditional
  `add_log(paste("Applying between-sample normalization..."))` block inline
  inside `mod_metab_norm_server()` at line `1973`.
- The focused gate now also freezes the between-sample normalization
  progress/apply shell contract in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R`
  as `runMetabNormBetweenSampleProgressApplyShell()`.
- The focused metabolomics normalization module helper gate reran green again
  after the between-sample normalization progress/apply shell seam via direct
  `testthat::test_file()`
  because this worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Classification refresh on April 17, 2026 now reports
  `R/mod_metab_norm.R`
  as `high-risk-wrapper` / `needs-seam-introduction` at `2233` lines with `28`
  top-level functions after the between-sample normalization progress/apply
  shell checkpoint.
- Next safe stop point:
  introduce one low-risk top-level seam for the RUV-III progress/apply shell
  inside the inline `runPipelineFn` callback, beginning with the
  `shiny::incProgress(..., detail = "Running RUV-III batch correction...")` /
  `if (input$ruv_mode != "skip")` block at line `1994`, then rerun the
  focused gate.
- The forty-second bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `runMetabNormRuvProgressApplyShell()`, which now owns the RUV-III step
  progress update, conditional apply log entry, delegation to
  `runMetabNormRuvOptimizationStep()`,
  `runMetabNormRuvCorrectionStep()`, and
  `runMetabNormRuvQcStep()`, plus the skipped-branch
  `ruv_corrected_obj` / `ruv_complete` state capture while the inline callback
  keeps only the downstream composite-QC helper call and plot-refresh tail.
- The inline `runPipelineFn` callback passed to
  `runMetabNormNormalizationObserverShell()`
  now delegates through that helper instead of keeping the
  `shiny::incProgress(..., detail = "Running RUV-III batch correction...")` /
  `if (input$ruv_mode != "skip")` block inline inside
  `mod_metab_norm_server()` at line `2086`.
- The focused gate now also freezes the RUV-III progress/apply shell contract
  in `tests/testthat/test-metab-03b-norm-module-helper-characterization.R`
  as `runMetabNormRuvProgressApplyShell()` across both apply and skip
  branches.
- The focused metabolomics normalization module helper gate reran green again
  after the RUV-III progress/apply shell seam via direct `testthat::test_file()`
  because this worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Classification refresh on April 17, 2026 now reports
  `R/mod_metab_norm.R`
  as `high-risk-wrapper` / `needs-seam-introduction` at `2296` lines with `30`
  top-level functions after the RUV-III progress/apply shell checkpoint.
- Next safe stop point:
  introduce one low-risk top-level seam for the final composite-QC /
  plot-refresh tail inside the inline `runPipelineFn` callback, beginning with
  `runMetabNormCompositeQcFigureStep()` at line `2108` and the
  `norm_data$plot_refresh_trigger <- norm_data$plot_refresh_trigger + 1` tail
  at line `2120`, then rerun the focused gate.
- The forty-third bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `runMetabNormCompositeQcRefreshShell()`, which now owns the final
  composite-QC helper delegation plus the closing `plot_refresh_trigger`
  increment while the inline callback keeps only the remaining orchestration
  sequence across the already-extracted pipeline helpers.
- The inline `runPipelineFn` callback passed to
  `runMetabNormNormalizationObserverShell()`
  now delegates through that helper instead of keeping the final
  `runMetabNormCompositeQcFigureStep()` call and
  `norm_data$plot_refresh_trigger <- norm_data$plot_refresh_trigger + 1`
  tail inline inside `mod_metab_norm_server()` at line `2143`.
- The focused gate now also freezes the final composite-QC / plot-refresh tail
  contract in `tests/testthat/test-metab-03b-norm-module-helper-characterization.R`
  as `runMetabNormCompositeQcRefreshShell()`.
- The focused metabolomics normalization module helper gate reran green again
  after the final composite-QC / plot-refresh tail seam via direct
  `testthat::test_file()`
  because this worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- The forty-fourth bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `runMetabNormNormalizationPipelineShell()`, which now owns the remaining
  `runPipelineFn` orchestration body: current-state lookup, plot-aesthetics
  capture, sequential pre-QC / ITSD / log2 / between-sample / post-QC / RUV
  handoffs, and the final `runMetabNormCompositeQcRefreshShell()` delegation
  while `runMetabNormNormalizationObserverShell()` keeps only the outer
  progress / completion / error wrapper.
- The inline `runPipelineFn` callback passed to
  `runMetabNormNormalizationObserverShell()`
  now delegates through that helper instead of keeping the remaining
  normalization orchestration shell inline inside `mod_metab_norm_server()`
  at line `2181`.
- The focused gate now also freezes the run-normalization pipeline shell
  contract in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R`
  as `runMetabNormNormalizationPipelineShell()`.
- The focused metabolomics normalization module helper gate reran green again
  after the run-normalization pipeline shell seam via direct
  `testthat::test_file()`
  because this worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Classification refresh on April 17, 2026 now reports
  `R/mod_metab_norm.R`
  as `high-risk-wrapper` / `needs-seam-introduction` at `2378` lines with `34`
  top-level functions after the runPipelineFn orchestration shell
  checkpoint.
- The seventh exact-source metabolomics normalization module wave now lives in
  `tools/refactor/manifest-metab-norm-module-wave7.yml`,
  with the reviewed staging artifacts in
  `tools/refactor/staging/wave7_metabolomics_norm_module_run_normalization_pipeline_shell/R/mod_metab_norm_server_helpers.R`
  and
  `tools/refactor/staging/wave7_metabolomics_norm_module_run_normalization_pipeline_shell/collate-metab-norm-module-wave7.txt`.
- `tools/refactor/check_wave_apply.R` passed for the live apply of
  `tools/refactor/manifest-metab-norm-module-wave7.yml`.
- The reviewed wave-7 apply moved
  `runMetabNormNormalizationPipelineShell()`
  from
  `R/mod_metab_norm.R`
  into
  `R/mod_metab_norm_server_helpers.R`
  while leaving the individual normalization step helpers live in
  `R/mod_metab_norm.R`.
- The focused gate now also freezes the run-normalization pipeline shell
  contract in
  `R/mod_metab_norm_server_helpers.R`
  as `runMetabNormNormalizationPipelineShell()`.
- The focused metabolomics normalization module helper gate reran green again
  after the wave-7 live apply via direct `testthat::test_file()` because this
  worktree still does not include `renv/activate.R` for `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `2257` lines with `32` top-level functions after the wave-7 live apply
  in `R/mod_metab_norm.R`.
- Next safe stop point:
  stage one reviewed exact-source wave for the opening run-normalization step
  helper cluster:
  `resolveMetabNormManualItsdFeatureIds()`,
  `runMetabNormPreNormalizationQcStep()`,
  `runMetabNormItsdProgressApplyShell()`, and
  `runMetabNormItsdNormalizationStep()`
  into `R/mod_metab_norm_server_helpers.R`, then rerun the focused gate.
- The eighth exact-source metabolomics normalization module wave now lives in
  `tools/refactor/manifest-metab-norm-module-wave8.yml`,
  with the reviewed staging artifacts in
  `tools/refactor/staging/wave8_metabolomics_norm_module_run_normalization_opening_step_helpers/R/mod_metab_norm_server_helpers.R`
  and
  `tools/refactor/staging/wave8_metabolomics_norm_module_run_normalization_opening_step_helpers/collate-metab-norm-module-wave8.txt`.
- `tools/refactor/check_wave_apply.R` passed for the live apply of
  `tools/refactor/manifest-metab-norm-module-wave8.yml`.
- The reviewed wave-8 apply moved
  `resolveMetabNormManualItsdFeatureIds()`,
  `runMetabNormPreNormalizationQcStep()`,
  `runMetabNormItsdProgressApplyShell()`, and
  `runMetabNormItsdNormalizationStep()`
  from `R/mod_metab_norm.R`
  into `R/mod_metab_norm_server_helpers.R`
  while leaving the remaining pre-RUV / RUV / composite run-normalization step
  helpers live in `R/mod_metab_norm.R`.
- `DESCRIPTION` `Collate:` already loaded
  `mod_metab_norm_server_helpers.R`
  before `mod_metab_norm.R`, so no additional load-order update was required
  for the wave-8 apply.
- The focused metabolomics normalization module helper gate reran green again
  after the wave-8 live apply via direct `testthat::test_file()` because this
  worktree still does not include `renv/activate.R` for `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `2070` lines with `24` top-level functions after the wave-8 live apply
  in `R/mod_metab_norm.R`.
- Next safe stop point:
  stage one reviewed exact-source wave for the remaining pre-RUV
  run-normalization step helper cluster:
  `runMetabNormLog2ProgressApplyShell()`,
  `runMetabNormLog2TransformationStep()`,
  `runMetabNormBetweenSampleProgressApplyShell()`,
  `runMetabNormBetweenSampleNormalizationStep()`, and
  `runMetabNormPostNormalizationQcStep()`
  into `R/mod_metab_norm_server_helpers.R`, then rerun the focused gate.
- The ninth exact-source metabolomics normalization module wave now lives in
  `tools/refactor/manifest-metab-norm-module-wave9.yml`,
  with the reviewed staging artifacts in
  `tools/refactor/staging/wave9_metabolomics_norm_module_run_normalization_pre_ruv_step_helpers/R/mod_metab_norm_server_helpers.R`
  and
  `tools/refactor/staging/wave9_metabolomics_norm_module_run_normalization_pre_ruv_step_helpers/collate-metab-norm-module-wave9.txt`.
- `tools/refactor/check_wave_apply.R` passed for the live apply of
  `tools/refactor/manifest-metab-norm-module-wave9.yml`.
- The reviewed wave-9 apply moved
  `runMetabNormLog2ProgressApplyShell()`,
  `runMetabNormLog2TransformationStep()`,
  `runMetabNormBetweenSampleProgressApplyShell()`,
  `runMetabNormBetweenSampleNormalizationStep()`, and
  `runMetabNormPostNormalizationQcStep()`
  from `R/mod_metab_norm.R`
  into `R/mod_metab_norm_server_helpers.R`
  while leaving the remaining RUV / composite run-normalization step helpers
  live in `R/mod_metab_norm.R`.
- `DESCRIPTION` `Collate:` already loaded
  `mod_metab_norm_server_helpers.R`
  before `mod_metab_norm.R`, so no additional load-order update was required
  for the wave-9 apply.
- The focused metabolomics normalization module helper gate reran green again
  after the wave-9 live apply via direct `testthat::test_file()` because this
  worktree still does not include `renv/activate.R` for `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `1895` lines with `14` top-level functions after the wave-9 live apply
  in `R/mod_metab_norm.R`.
- Next safe stop point:
  stage one reviewed exact-source wave for the remaining RUV / composite
  run-normalization step helper cluster:
  `runMetabNormRuvProgressApplyShell()`,
  `runMetabNormRuvOptimizationStep()`,
  `runMetabNormRuvCorrectionStep()`,
  `runMetabNormRuvQcStep()`,
  `runMetabNormCompositeQcFigureStep()`, and
  `runMetabNormCompositeQcRefreshShell()`
  into `R/mod_metab_norm_server_helpers.R`, then rerun the focused gate.
- The tenth exact-source metabolomics normalization module wave now lives in
  `tools/refactor/manifest-metab-norm-module-wave10.yml`,
  with the reviewed staging artifacts in
  `tools/refactor/staging/wave10_metabolomics_norm_module_run_normalization_ruv_composite_step_helpers/R/mod_metab_norm_server_helpers.R`.
- `tools/refactor/check_wave_apply.R` passed for the live apply of
  `tools/refactor/manifest-metab-norm-module-wave10.yml`.
- The reviewed wave-10 apply moved
  `runMetabNormRuvProgressApplyShell()`,
  `runMetabNormRuvOptimizationStep()`,
  `runMetabNormRuvCorrectionStep()`,
  `runMetabNormRuvQcStep()`,
  `runMetabNormCompositeQcFigureStep()`, and
  `runMetabNormCompositeQcRefreshShell()`
  from `R/mod_metab_norm.R`
  into `R/mod_metab_norm_server_helpers.R`,
  leaving only the public `mod_metab_norm_ui()` and
  `mod_metab_norm_server()` top-level functions live in
  `R/mod_metab_norm.R`.
- `DESCRIPTION` `Collate:` already loaded
  `mod_metab_norm_server_helpers.R`
  before `mod_metab_norm.R`, so no additional load-order update was required
  for the wave-10 apply.
- The focused metabolomics normalization module helper gate reran green again
  after the wave-10 live apply via direct `testthat::test_file()` because this
  worktree still does not include `renv/activate.R` for `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `1513` lines with `2` top-level functions after the wave-10 live apply
  in `R/mod_metab_norm.R`.
- Next safe stop point:
  introduce one bounded top-level seam in `R/mod_metab_norm.R` for the nested
  composite image-assembly helper cluster inside `mod_metab_norm_server()`:
  `generateCompositeFromFiles()`,
  `createLabelPlot()`,
  `createTitlePlot()`, and
  `loadImageAsPlot()`,
  then rerun the focused gate.
- The next bounded seam then went live in `R/mod_metab_norm.R` as
  `buildMetabNormLabelPlot()`,
  `buildMetabNormTitlePlot()`,
  `loadMetabNormImageAsPlot()`, and
  `generateMetabNormCompositeFromFiles()`, and
  `mod_metab_norm_server()` now delegated the composite image-assembly cluster
  through `generateMetabNormCompositeFromFiles()` instead of keeping the
  helper bodies nested in the wrapper.
- The focused gate already covered that composite image-assembly seam contract,
  and reran green again after the seam via direct `testthat::test_file()`
  because this worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Live classification then remained
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `1514` lines with `7` top-level functions after the composite
  image-assembly seam checkpoint in `R/mod_metab_norm.R`.
- Next safe stop point:
  stage one reviewed exact-source wave for the new composite image-assembly
  helper cluster
  `buildMetabNormLabelPlot()`,
  `buildMetabNormTitlePlot()`,
  `loadMetabNormImageAsPlot()`, and
  `generateMetabNormCompositeFromFiles()`
  into `R/mod_metab_norm_server_helpers.R`, then rerun the focused gate.
- The eleventh exact-source metabolomics normalization module wave now lives in
  `tools/refactor/manifest-metab-norm-module-wave11.yml`,
  with the reviewed staging artifacts in
  `tools/refactor/staging/wave11_metabolomics_norm_module_composite_image_assembly_helpers/R/mod_metab_norm_server_helpers.R`.
- `tools/refactor/check_wave_apply.R` passed for the live apply of
  `tools/refactor/manifest-metab-norm-module-wave11.yml`.
- The reviewed wave-11 apply moved
  `buildMetabNormLabelPlot()`,
  `buildMetabNormTitlePlot()`,
  `loadMetabNormImageAsPlot()`, and
  `generateMetabNormCompositeFromFiles()`
  from `R/mod_metab_norm.R`
  into `R/mod_metab_norm_server_helpers.R`, and
  `mod_metab_norm_server()` at line `1186` now delegates the composite
  image-assembly path through the helper file instead of keeping the helper
  cluster in the wrapper source file.
- `DESCRIPTION` `Collate:` already loaded
  `mod_metab_norm_server_helpers.R`
  before `mod_metab_norm.R`, so no additional load-order update was required
  for the wave-11 apply.
- The focused metabolomics normalization module helper gate reran green again
  after the wave-11 live apply via direct `testthat::test_file()` because this
  worktree still does not include `renv/activate.R` for `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `1360` lines with `2` top-level functions after the wave-11 live apply
  in `R/mod_metab_norm.R`.
- The next bounded seam now lives in
  `R/mod_metab_norm.R`
  as
  `generateMetabNormPreNormalizationQc()`, and the selected-tab auto-trigger
  block in `mod_metab_norm_server()` now delegates through that helper instead
  of keeping the nested `generatePreNormalizationQc()` body inline in the
  wrapper.
- The focused gate now also freezes the auto pre-normalization QC seam
  contracts in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R`
  for `generateMetabNormPreNormalizationQc()`, and reran green again after the
  seam via direct `testthat::test_file()` because this worktree still does not
  include `renv/activate.R` for `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `1380` lines with `4` top-level functions after the pre-normalization QC
  helper seam checkpoint in `R/mod_metab_norm.R`.
- The twelfth exact-source metabolomics normalization module wave now lives in
  `tools/refactor/manifest-metab-norm-module-wave12.yml`,
  with the reviewed staging artifacts in
  `tools/refactor/staging/wave12_metabolomics_norm_module_pre_normalization_qc_helper/R/mod_metab_norm_server_helpers.R`
  and
  `tools/refactor/staging/wave12_metabolomics_norm_module_pre_normalization_qc_helper/collate-metab-norm-module-wave12.txt`.
- `tools/refactor/check_wave_apply.R` passed for the live apply of
  `tools/refactor/manifest-metab-norm-module-wave12.yml`.
- The reviewed wave-12 apply moved
  `generateMetabNormPreNormalizationQc()`
  from `R/mod_metab_norm.R`
  into `R/mod_metab_norm_server_helpers.R`, and
  `mod_metab_norm_server()` at line `855` now delegates the auto
  pre-normalization QC path through the helper file instead of keeping the
  helper body live in the wrapper source file.
- `DESCRIPTION` `Collate:` already loaded
  `mod_metab_norm_server_helpers.R`
  before `mod_metab_norm.R`, so no additional load-order update was required
  for the wave-12 apply.
- The focused metabolomics normalization module helper gate reran green again
  after the wave-12 live apply via direct `testthat::test_file()` because this
  worktree still does not include `renv/activate.R` for `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `1328` lines with `2` top-level functions after the wave-12 live apply
  in `R/mod_metab_norm.R`.
- The next bounded live seam now lives in
  `R/mod_metab_norm.R`
  at line `680`
  as
  `getPlotAesthetics()`, which now owns the module's plot color/shape default
  fallback and explicit input passthrough for the wrapper.
- The auto pre-normalization QC observer at line `853`, the normalization
  pipeline shell handoff at line `1134`, and the final QC PCA render at line
  `1295` now delegate through that top-level helper instead of keeping the
  plot-aesthetics resolution inline inside `mod_metab_norm_server()`.
- The focused metabolomics normalization module helper gate now also freezes
  the plot-aesthetics seam contract in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R`,
  and reran green again via direct `testthat::test_file()` because this
  worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `1339` lines with `3` top-level functions after the plot-aesthetics seam
  checkpoint in `R/mod_metab_norm.R`.
- Next safe stop point:
  introduce one bounded top-level seam in `R/mod_metab_norm.R` for the nested
  assay-label renderer `render_assay_label()` at line `1022`, then rerun the
  focused gate.
- The next bounded live seam now lives in
  `R/mod_metab_norm.R`
  at line `700`
  as
  `renderMetabNormAssayLabel()`, which now owns the static assay-label render
  binding for the module's PCA, density, RLE, and correlation tabs.
- The assay-label output bindings in
  `R/mod_metab_norm.R`
  at lines `1035` through `1045` now delegate through that top-level helper
  instead of keeping the nested `render_assay_label()` helper inside
  `mod_metab_norm_server()`.
- The focused metabolomics normalization module helper gate now also freezes
  the assay-label seam contract in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R`
  and reran green again via direct `testthat::test_file()` because this
  worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `1346` lines with `4` top-level functions after the assay-label seam
  checkpoint in `R/mod_metab_norm.R`.
- Next safe stop point:
  introduce one bounded top-level seam in `R/mod_metab_norm.R` for the nested
  QC image renderer `render_qc_image_for_assay()` at line `788`, then rerun
  the focused gate.
- The next bounded live seam now lives in
  `R/mod_metab_norm.R`
  at line `718`
  as
  `renderMetabNormQcImageForAssay()`, which now owns the static QC image
  render binding for the module's PCA, density, RLE, and correlation assay
  outputs.
- The QC image output bindings in
  `R/mod_metab_norm.R`
  at line `1010`
  now delegate through that top-level helper instead of keeping the nested
  `render_qc_image_for_assay()` helper inside `mod_metab_norm_server()`.
- The focused metabolomics normalization module helper gate now also freezes
  the QC image seam contract in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R`
  at line `152`
  and reran green again via direct `testthat::test_file()` because this
  worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- The next bounded live seam now lives in
  `R/mod_metab_norm.R`
  at line `765`
  as `appendMetabNormNormalizationLog()`, which now owns the normalization log
  timestamp formatting, entry assembly, and reactive log-vector append for the
  wrapper.
- The nested `add_log()` helper in
  `R/mod_metab_norm.R`
  at line `841`
  now delegates through that top-level helper instead of keeping the
  timestamp/append block inline inside `mod_metab_norm_server()`.
- The focused metabolomics normalization module helper gate now also freezes
  the normalization log append seam contract in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R`
  at line `120`
  and reran green again via direct `testthat::test_file()` because this
  worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `1370` lines with `8` top-level functions after the normalization log
  append seam checkpoint in `R/mod_metab_norm.R`.
- Next safe stop point:
  introduce one bounded top-level seam in `R/mod_metab_norm.R` for the nested
  normalization log `renderText` binding at line `960`, then rerun the focused
  gate.
- The next bounded live seam now lives in
  `R/mod_metab_norm.R`
  at line `784`
  as `renderMetabNormNormalizationLog()`, which now owns the normalization log
  empty-state fallback and newline-collapsed `renderText` output for the
  wrapper.
- The `output$norm_log` binding in
  `R/mod_metab_norm.R`
  at line `978`
  now delegates through that top-level helper instead of keeping the
  normalization log `renderText` body inline inside `mod_metab_norm_server()`.
- The focused metabolomics normalization module helper gate now also freezes
  the normalization log render seam contract in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R`
  at line `141`
  and reran green again via direct `testthat::test_file()` because this
  worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `1385` lines with `9` top-level functions after the normalization log
  render seam checkpoint in `R/mod_metab_norm.R`.
- The next bounded live seam now lives in
  `R/mod_metab_norm.R`
  at line `802`
  as `renderMetabNormCorrelationFilterSummary()`, which now owns the
  correlation-filter summary incomplete-state fallback, filtered/original
  source-object selection, and handoff into
  `buildMetabNormCorrelationFilterSummary()`.
- The `output$correlation_filter_summary` binding in
  `R/mod_metab_norm.R`
  at line `1347`
  now delegates through that top-level helper instead of keeping the
  correlation filter summary `renderText` body inline inside
  `mod_metab_norm_server()`.
- The focused metabolomics normalization module helper gate now also freezes
  the correlation filter summary render seam contract in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R`
  at line `168`
  and reran green again via direct `testthat::test_file()` because this
  worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `1398` lines with `10` top-level functions after the correlation filter
  summary render seam checkpoint in `R/mod_metab_norm.R`.
- The next bounded live seam now lives in
  `R/mod_metab_norm.R`
  at line `831`
  as `renderMetabNormFinalQcPlot()`, which now owns the final-QC readiness
  gate, render-state resolution handoff, fallback-plot return, plot-aesthetics
  lookup, and PCA plot builder delegation for the wrapper.
- The `output$final_qc_plot` binding in
  `R/mod_metab_norm.R`
  at line `1392`
  now delegates through that top-level helper instead of keeping the final QC
  plot `renderPlot` body inline inside `mod_metab_norm_server()`.
- The focused metabolomics normalization module helper gate now also freezes
  the final-QC render seam contract in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R`
  at line `241`
  and reran green again via direct `testthat::test_file()` because this
  worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `1416` lines with `12` top-level functions after the final QC plot render
  seam checkpoint in `R/mod_metab_norm.R`.
- The thirteenth exact-source metabolomics normalization module wave now lives
  in `tools/refactor/manifest-metab-norm-module-wave13.yml`, with the reviewed
  staging artifacts in
  `tools/refactor/staging/wave13_metabolomics_norm_module_render_helper_cluster/R/mod_metab_norm_server_helpers.R`
  and
  `tools/refactor/staging/wave13_metabolomics_norm_module_render_helper_cluster/collate-metab-norm-module-wave13.txt`.
- `tools/refactor/check_wave_apply.R` passed for the live apply of
  `tools/refactor/manifest-metab-norm-module-wave13.yml`.
- The reviewed wave-13 apply moved
  `getPlotAesthetics()`,
  `renderMetabNormAssayLabel()`,
  `renderMetabNormQcImageForAssay()`,
  `appendMetabNormNormalizationLog()`,
  `renderMetabNormNormalizationLog()`,
  `renderMetabNormCorrelationFilterSummary()`, and
  `renderMetabNormFinalQcPlot()`
  from `R/mod_metab_norm.R`
  into `R/mod_metab_norm_server_helpers.R`,
  leaving only the public `mod_metab_norm_ui()` and
  `mod_metab_norm_server()` top-level functions live in
  `R/mod_metab_norm.R`.
- `DESCRIPTION` `Collate:` already loaded
  `mod_metab_norm_server_helpers.R`
  before `mod_metab_norm.R`, so no additional load-order update was required
  for the wave-13 apply.
- The focused metabolomics normalization module helper gate reran green again
  after the wave-13 live apply via direct `testthat::test_file()` because this
  worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- The next bounded live seam now sits in
  `R/mod_metab_norm.R`
  as `updateMetabNormDesignDrivenChoices()`, which now owns the
  design-matrix-driven plot-aesthetic and RUV grouping input-choice updates
  before the wrapper falls through to the remaining render/UI and observer
  registrations.
- The design-matrix observer binding in
  `R/mod_metab_norm.R`
  now delegates through that helper instead of keeping the input-choice update
  block inline inside `mod_metab_norm_server()` at line `883`.
- The focused metabolomics normalization module helper gate reran green again
  after the design-driven choice observer seam via direct
  `testthat::test_file()` because this worktree still does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `1264` lines with `3` top-level functions after the design-driven choice
  observer seam in `R/mod_metab_norm.R`.
- The next bounded live seam now sits in
  `R/mod_metab_norm.R`
  at line `737`
  as `initializeMetabNormAssayNames()`, which now owns the assay-name
  detection, `MetaboliteAssayData` guard, reactive assay-name assignment, and
  warning-logged fallback for the wrapper before control returns to the
  remaining auto-trigger and observer registrations.
- The assay-name initialization observer in
  `R/mod_metab_norm.R`
  at line `846`
  now delegates through that helper instead of keeping the state-manager
  lookup, assay detection, and warning handling inline inside
  `mod_metab_norm_server()`.
- The focused metabolomics normalization module helper gate now also freezes
  the assay-name initialization seam contract in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R`
  and reran green again via direct `testthat::test_file()` because this
  worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `1284` lines with `5` top-level functions after the assay-name
  initialization observer seam in `R/mod_metab_norm.R`.
- The next bounded live seam now sits in
  `R/mod_metab_norm.R`
  at line `773`
  as `runMetabNormAutoPreNormalizationQcObserverShell()`, which now owns the
  selected-tab `"norm"` gate, pre-generated short-circuit, state-manager
  lookup, `MetaboliteAssayData` guard, and `withProgress()` handoff into the
  existing pre-normalization QC generator before control returns to the
  remaining render/UI and observer registrations.
- The selected-tab auto-trigger pre-normalization QC observer in
  `R/mod_metab_norm.R`
  at line `929`
  now delegates through that helper instead of keeping the selected-tab gate,
  state-manager `req()`, S4 guard, and progress-wrapped helper handoff inline
  inside `mod_metab_norm_server()`.
- The focused metabolomics normalization module helper gate now also freezes
  the selected-tab auto pre-normalization QC observer seam contract in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R`
  and reran green again via direct `testthat::test_file()` because this
  worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `1326` lines with `8` top-level functions after the auto-trigger
  pre-normalization QC observer seam in `R/mod_metab_norm.R`.
- Next safe stop point:
  introduce one bounded top-level seam in `R/mod_metab_norm.R` for the dynamic
  ITSD selection tables `renderUI` at line `962`, then rerun the focused gate.
- Update on April 17, 2026:
  the next bounded live seam now sits in `R/mod_metab_norm.R` at line `1170`
  as `runMetabNormAssayLabelBindingShell()`, which now owns the eight static
  assay-label output bindings across the PCA, density, RLE, and correlation
  tabs before control returns to `mod_metab_norm_server()`.
- The static assay-label output binding cluster in `R/mod_metab_norm.R` at
  line `1362` now delegates through that helper instead of keeping the eight
  `renderMetabNormAssayLabel()` registrations inline inside
  `mod_metab_norm_server()`.
- The focused metabolomics normalization module helper gate now also freezes
  the static assay-label binding shell contract in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R` at
  line `2020` and reran green again via direct `testthat::test_file()`
  because this worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `1525` lines with `23` top-level functions after the static assay-label
  binding seam checkpoint in `R/mod_metab_norm.R`.
- Next safe stop point:
  introduce one bounded top-level seam in `R/mod_metab_norm.R` for the static
  QC plot image output binding cluster at line `1323`, then rerun the focused
  gate.
- Update on April 17, 2026:
  the next bounded live seam now sits in `R/mod_metab_norm.R` at line `1170`
  as `runMetabNormQcImageBindingShell()`, which now owns the twenty-four
  static QC plot image output bindings across the PCA, density, RLE, and
  correlation tabs before control returns to `mod_metab_norm_server()`.
- The static QC plot image output binding cluster in `R/mod_metab_norm.R` at
  line `1365` now delegates through that helper instead of keeping the
  twenty-four `renderMetabNormQcImageForAssay()` registrations inline inside
  `mod_metab_norm_server()`.
- The focused metabolomics normalization module helper gate now also freezes
  the static QC image binding shell contract in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R` at
  line `2062` and reran green again via direct `testthat::test_file()`
  because this worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `1537` lines with `24` top-level functions after the static QC image
  binding seam checkpoint in `R/mod_metab_norm.R`.
- Next safe stop point:
  introduce one bounded top-level seam in `R/mod_metab_norm.R` for the ITSD
  selection tracking observer wrapper at line `1393`, then rerun the focused
  gate.
- Update on April 17, 2026:
  the next bounded live seam now sits in `R/mod_metab_norm.R` at line `956`
  as `runMetabNormItsdSelectionTrackingObserverShell()`, which now owns the
  ITSD selection tracking observer wrapper's assay-name `req()` gate and
  delegation into `registerMetabNormItsdSelectionTracking()` before control
  returns to `mod_metab_norm_server()`.
- The ITSD selection tracking observe block in `R/mod_metab_norm.R` at line
  `1411` now delegates through that helper instead of keeping the inline
  `shiny::req(norm_data$assay_names)` gate and
  `registerMetabNormItsdSelectionTracking()` handoff inside
  `mod_metab_norm_server()`.
- The focused metabolomics normalization module helper gate now also freezes
  the ITSD selection tracking observer-shell seam contract in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R` at
  lines `1508` and `1540` and reran green again via direct
  `testthat::test_file()` because this worktree still does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `1553` lines with `25` top-level functions after the ITSD selection
  tracking observer-shell seam checkpoint in `R/mod_metab_norm.R`.
- Next safe stop point:
  introduce one bounded top-level seam in `R/mod_metab_norm.R` for the main
  normalization pipeline observer wrapper at line `1421`, then rerun the
  focused gate.
- Update on April 17, 2026:
  the next bounded live seam now sits in `R/mod_metab_norm.R` at line `1249`
  as `runMetabNormNormalizationObserverWrapper()`, which now owns the main
  run-normalization observer wrapper's
  `runMetabNormNormalizationObserverShell()` handoff, `inputValues` assembly,
  plot-aesthetics callback construction, and delegation into
  `runMetabNormNormalizationPipelineShell()` before control returns to
  `mod_metab_norm_server()`.
- The `input$run_normalization` observe block in `R/mod_metab_norm.R` at line
  `1483` now delegates through that helper instead of keeping the inline
  `runMetabNormNormalizationObserverShell()` handoff and
  `runMetabNormNormalizationPipelineShell()` input assembly inside
  `mod_metab_norm_server()`.
- The focused metabolomics normalization module helper gate now also freezes
  the run-normalization observer-wrapper seam contract in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R` at
  line `5148` and reran green again via direct `testthat::test_file()`
  because this worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `1585` lines with `27` top-level functions after the main
  normalization pipeline observer-wrapper seam checkpoint in
  `R/mod_metab_norm.R`.
- Next safe stop point:
  introduce one bounded top-level seam in `R/mod_metab_norm.R` for the reset
  normalization observer wrapper at line `1500`, then rerun the focused gate.
- Update on April 17, 2026:
  the next bounded live seam now sits in `R/mod_metab_norm.R` at line `1311`
  as `runMetabNormResetNormalizationObserverWrapper()`, which now owns the
  reset-normalization observer wrapper's
  `runMetabNormResetNormalizationObserverShell()` handoff and dependency
  forwarding before control returns to `mod_metab_norm_server()`.
- The `input$reset_normalization` observe block in `R/mod_metab_norm.R` at
  line `1520` now delegates through that helper instead of keeping the inline
  `runMetabNormResetNormalizationObserverShell()` handoff inside
  `mod_metab_norm_server()`.
- The focused metabolomics normalization module helper gate now also freezes
  the reset-normalization observer-wrapper seam contract in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R` at
  line `5289` and reran green again via direct `testthat::test_file()`
  because this worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `1604` lines with `29` top-level functions after the reset
  normalization observer-wrapper seam checkpoint in `R/mod_metab_norm.R`.
- Next safe stop point:
  introduce one bounded top-level seam in `R/mod_metab_norm.R` for the
  apply-correlation observer wrapper at line `1542`, then rerun the focused
  gate.
- Update on April 17, 2026:
  the next bounded live seam now sits in `R/mod_metab_norm.R` at line `1330`
  as `runMetabNormApplyCorrelationObserverWrapper()`, which now owns the
  apply-correlation observer wrapper's
  `runMetabNormApplyCorrelationObserverShell()` handoff, threshold/grouping
  input forwarding, and observer-entry dependency forwarding before control
  returns to `mod_metab_norm_server()`.
- The `input$apply_correlation_filter` observe block in
  `R/mod_metab_norm.R` at line `1568` now delegates through that helper
  instead of keeping the inline
  `runMetabNormApplyCorrelationObserverShell()` handoff inside
  `mod_metab_norm_server()`.
- The focused metabolomics normalization module helper gate now also freezes
  the apply-correlation observer-wrapper seam contract in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R` at
  line `5337` and reran green again via direct `testthat::test_file()`
  because this worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `1628` lines with `31` top-level functions after the
  apply-correlation observer-wrapper seam checkpoint in
  `R/mod_metab_norm.R`.
- Next safe stop point:
  introduce one bounded top-level seam in `R/mod_metab_norm.R` for the
  skip-correlation observer wrapper at line `1583`, then rerun the focused
  gate.
- Update on April 17, 2026:
  the next bounded live seam now sits in `R/mod_metab_norm.R` at line `1356`
  as `runMetabNormSkipCorrelationObserverWrapper()`, which now owns the
  skip-correlation observer wrapper's
  `runMetabNormSkipCorrelationObserverShell()` handoff and observer-entry
  dependency forwarding before control returns to `mod_metab_norm_server()`.
- The `input$skip_correlation_filter` observe block in
  `R/mod_metab_norm.R` at line `1604` now delegates through that helper
  instead of keeping the inline
  `runMetabNormSkipCorrelationObserverShell()` handoff inside
  `mod_metab_norm_server()`.
- The focused metabolomics normalization module helper gate now also freezes
  the skip-correlation observer-wrapper seam contract in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R` at
  line `5412` and reran green again via direct `testthat::test_file()`
  because this worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `1648` lines with `33` top-level functions after the
  skip-correlation observer-wrapper seam checkpoint in
  `R/mod_metab_norm.R`.
- Next safe stop point:
  introduce one bounded top-level seam in `R/mod_metab_norm.R` for the
  export-session observer wrapper at line `1633`, then rerun the focused
  gate.
- Update on April 17, 2026:
  the next bounded live seam now sits in `R/mod_metab_norm.R` at line `1377`
  as `runMetabNormExportSessionObserverWrapper()`, which now owns the
  export-session observer wrapper's
  `runMetabNormExportSessionObserverShell()` handoff, export input snapshot
  assembly, and dependency forwarding before control returns to
  `mod_metab_norm_server()`.
- The `input$export_session` observe block in `R/mod_metab_norm.R` at line
  `1673` now delegates through that helper instead of keeping the inline
  `runMetabNormExportSessionObserverShell()` handoff inside
  `mod_metab_norm_server()`.
- The focused metabolomics normalization module helper gate now also freezes
  the export-session observer-wrapper seam contract in
  `tests/testthat/test-metab-03b-norm-module-helper-characterization.R` at
  line `5471` and reran green again via direct `testthat::test_file()`
  because this worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `1686` lines with `35` top-level functions after the
  export-session observer-wrapper seam checkpoint in
  `R/mod_metab_norm.R`.
- The fifteenth exact-source metabolomics normalization module wave now lives
  in `tools/refactor/manifest-metab-norm-module-wave15.yml`,
  with reviewed staging artifacts in
  `tools/refactor/staging/wave15_metabolomics_norm_module_ui_binding_helper_cluster/R/mod_metab_norm_server_helpers.R`
  and
  `tools/refactor/staging/wave15_metabolomics_norm_module_ui_binding_helper_cluster/collate-metab-norm-module-wave15.txt`.
- The wave-15 live apply now further extends
  `R/mod_metab_norm_server_helpers.R`
  and owns the remaining top-level UI/binding helper cluster from
  `updateMetabNormDesignDrivenChoices()` through
  `runMetabNormAssayLabelBindingShell()`
  after exact-source removal from
  `R/mod_metab_norm.R`.
- The module server bindings in
  `R/mod_metab_norm.R`
  at lines `836`, `847`, `863`, `879`, `887`, `896`, `905`, `914`, `925`,
  and `965` continue delegating through those helpers after the live apply.
- `tools/refactor/verify_refactor.R` passed again for
  `tools/refactor/manifest-metab-norm-module-wave15.yml`
  before the live apply, and `tools/refactor/check_wave_apply.R` passed after
  the live apply.
- `DESCRIPTION` `Collate:` already loads
  `mod_metab_norm_server_helpers.R`
  before
  `mod_metab_norm.R`,
  so no additional load-order update was required for the wave-15 apply.
- The focused metabolomics normalization module helper gate reran green again
  after the wave-15 apply via direct `testthat::test_file()`
  because this worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction` at `1031` lines with `4` top-level functions
  after the wave-15 UI/binding helper apply in `R/mod_metab_norm.R`.
- Next safe stop point:
  stage one reviewed exact-source helper wave for the remaining top-level
  normalization observer wrapper exposed in `R/mod_metab_norm.R` at line
  `696`, into
  `R/mod_metab_norm_server_helpers.R`, then rerun the focused gate.
- The sixteenth exact-source metabolomics normalization module wave now lives
  in `tools/refactor/manifest-metab-norm-module-wave16.yml`,
  with reviewed staging artifacts in
  `tools/refactor/staging/wave16_metabolomics_norm_module_normalization_observer_wrapper/R/mod_metab_norm_server_helpers.R`
  and
  `tools/refactor/staging/wave16_metabolomics_norm_module_normalization_observer_wrapper/collate-metab-norm-module-wave16.txt`.
- The wave-16 live apply now further extends
  `R/mod_metab_norm_server_helpers.R`
  and owns the remaining top-level normalization observer wrapper
  `runMetabNormNormalizationObserverWrapper()` after exact-source removal from
  `R/mod_metab_norm.R`.
- The `input$run_normalization` observe block in `R/mod_metab_norm.R` at line
  `874` continues delegating through
  `runMetabNormNormalizationObserverWrapper()` after the live apply, with the
  helper now defined in `R/mod_metab_norm_server_helpers.R` at line `2876`.
- `tools/refactor/verify_refactor.R` passed again for
  `tools/refactor/manifest-metab-norm-module-wave16.yml`
  before the live apply, and `tools/refactor/check_wave_apply.R` passed after
  the live apply.
- The focused metabolomics normalization module helper gate reran green again
  after the wave-16 apply via direct `testthat::test_file()`
  because this worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction` at `970` lines with `2` top-level functions after
  the wave-16 normalization observer wrapper apply in `R/mod_metab_norm.R`.
- The seventeenth exact-source metabolomics normalization module wave now
  lives in `tools/refactor/manifest-metab-norm-module-wave17.yml`, with
  reviewed staging artifacts in
  `tools/refactor/staging/wave17_metabolomics_norm_module_wrapper_entrypoints/R/mod_metab_norm_ui.R`,
  `tools/refactor/staging/wave17_metabolomics_norm_module_wrapper_entrypoints/R/mod_metab_norm_server.R`,
  and
  `tools/refactor/staging/wave17_metabolomics_norm_module_wrapper_entrypoints/collate-metab-norm-module-wave17.txt`.
- The wave-17 live apply created `R/mod_metab_norm_ui.R` and
  `R/mod_metab_norm_server.R`, moved `mod_metab_norm_ui()` and
  `mod_metab_norm_server()` out of `R/mod_metab_norm.R`, and left
  `R/mod_metab_norm.R` as the breadcrumb stub for the module identity.
- `DESCRIPTION` `Collate:` now loads `mod_metab_norm_server_helpers.R`,
  `mod_metab_norm_ui.R`, `mod_metab_norm_server.R`, and then
  `mod_metab_norm.R` so package load order matches the final wrapper split.
- `tools/refactor/verify_refactor.R` passed again for
  `tools/refactor/manifest-metab-norm-module-wave17.yml` before the live
  apply, and `tools/refactor/check_wave_apply.R` passed after the live apply.
- The focused metabolomics normalization module helper gate reran green again
  after the wave-17 apply via direct `testthat::test_file()`
  because this worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Live classification now records `R/mod_metab_norm.R` at
  `direct-extraction-ready` with `69` lines and `0` top-level functions after
  the wave-17 wrapper entrypoint apply.
- Treat the manual metabolomics normalization wrapper target as complete:
  `R/mod_metab_norm.R` is now a breadcrumb stub and bucket 0 can move on to
  the next manual target.
