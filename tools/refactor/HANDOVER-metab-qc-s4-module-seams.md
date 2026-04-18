# Metabolomics QC S4 Module Seam Map

## Goal

Introduce bounded top-level seams in `R/mod_metab_qc_s4.R` while keeping the
live metabolomics QC finalization module behavior frozen behind the existing
public entry points.

## Current Position In The Flow

- target: `/home/doktersmol/Documents/MultiScholaR-metab-safe/R/mod_metab_qc_s4.R`
- classification: `review` (post-checkpoint structural classifier)
- checkpoint reached: `reviewed wave-2 exact-source apply moved runMetabQcS4ServerBody() into R/mod_metab_qc_s4_server_helpers.R while the public wrapper stayed stable`
- next step: `Treat the metabolomics QC S4 module wrapper target as complete and archived; keep the focused QC S4-module gate as the regression surface and do not reopen R/mod_metab_qc_s4.R unless a real regression appears.`

## Existing Safety Net

- `tests/testthat/test-metab-01s-qc-s4-module-characterization.R`

## Notes

Manual bucket 0 metabolomics QC finalization-module stabilization target.

- This file records the current metabolomics QC S4 module seam stop point.
- Classification refreshed on April 18, 2026 now reports
  `R/mod_metab_qc_s4.R`
  at `171` lines with `2` top-level functions, a `54` line largest
  top-level function, and labels
  `high-risk-wrapper`
  plus
  `needs-seam-introduction`.
- The seam-introduction bullets below preserve the original in-file landing
  points for the frozen helper surface; after the reviewed wave-1 live apply,
  those `21` helper definitions now live in
  `R/mod_metab_qc_s4_server_helpers.R`
  while
  `R/mod_metab_qc_s4.R`
  retains only
  `mod_metab_qc_s4_ui()`
  and
  `mod_metab_qc_s4_server()`.
- The focused gate now lives in
  `tests/testthat/test-metab-01s-qc-s4-module-characterization.R`
  and now loads
  `R/mod_metab_qc_s4_server_helpers.R`
  plus
  `R/mod_metab_qc_s4.R`
  directly; it freezes the live helper-surface contracts for
  `buildMetabQcS4DataSummaryUi()`
  plus
  `getMetabQcS4DataSummaryState()`
  plus
  `buildMetabQcS4DataSummaryRenderOutput()`
  plus
  `buildMetabQcS4StateHistoryUi()`
  plus
  `getMetabQcS4StateHistory()`
  plus
  `buildMetabQcS4StateHistoryRenderOutput()`
  plus
  `buildMetabQcS4AssayStatsDatatable()`
  plus
  `getMetabQcS4AssayStatsState()`
  plus
  `buildMetabQcS4AssayStatsRenderOutput()`
  plus
  `buildMetabQcS4FilterPlotRenderOutput()`
  plus
  `renderMetabQcS4FilterPlot()`
  plus
  `buildMetabQcS4FinalizeResultsText()`
  plus
  `getMetabQcS4FinalizeState()`
  plus
  `validateMetabQcS4FinalizeState()`
  plus
  `saveMetabQcS4CompletedState()`
  plus
  `completeMetabQcS4TabStatus()`
  plus
  `getMetabQcS4FinalizeHistory()`
  plus
  `updateMetabQcS4TrackingPlot()`
  plus
  `reportMetabQcS4FinalizeSuccess()`
  plus
  `reportMetabQcS4FinalizeError()`
  plus
  `runMetabQcS4FinalizeWorkflow()`
  plus
  `runMetabQcS4ServerBody()`
  and the public-wrapper handoff from
  `mod_metab_qc_s4_server()`
  into the server-body seam plus the downstream helper delegation from
  `runMetabQcS4ServerBody()`
  into the `state_history` render path, the `state_history` render-output
  seam, the `data_summary` render path, the `data_summary` render-path
  current-state fetch seam, the `data_summary` render-output seam, the
  `assay_stats_table` render path, the `assay_stats_table` render-output seam,
  the `filter_plot` render-output seam, and the finalize observer workflow
  seam.
- The first bounded live seam now sits in
  `R/mod_metab_qc_s4.R`
  at line `94`
  as `buildMetabQcS4DataSummaryUi()`, which now owns the invalid-state
  placeholder, per-assay metabolite counting, numeric-sample counting,
  design-group counting, and summary-table assembly before control returns to
  the metabolomics QC finalization module wrapper.
- The `data_summary` render path in
  `R/mod_metab_qc_s4.R`
  at line `649`
  now delegates through that helper instead of building the current-data
  summary UI inline inside `mod_metab_qc_s4_server()`.
- The second bounded live seam now sits in
  `R/mod_metab_qc_s4.R`
  at line `176`
  as `buildMetabQcS4StateHistoryUi()`, which now owns the empty-history
  placeholder, ordered state-history list assembly, per-step icon selection,
  current-state bolding, and current-marker span before control returns to the
  metabolomics QC finalization module wrapper.
- The `state_history` render path in
  `R/mod_metab_qc_s4.R`
  at line `640`
  now delegates through that helper instead of building the processing-history
  list inline inside `mod_metab_qc_s4_server()`.
- The eleventh bounded live seam now sits in
  `R/mod_metab_qc_s4.R`
  at line `228`
  as `getMetabQcS4StateHistory()`, which now owns the
  `state_manager$getHistory()` handoff plus `character(0)` fallback that feeds
  `buildMetabQcS4StateHistoryUi()` before control returns to the metabolomics
  QC finalization module wrapper.
- The `state_history` render path in
  `R/mod_metab_qc_s4.R`
  at line `640`
  now delegates through that helper instead of wiring the history-fetch
  `tryCatch()` inline inside `mod_metab_qc_s4_server()`.
- The twelfth bounded live seam now sits in
  `R/mod_metab_qc_s4.R`
  at line `490`
  as `getMetabQcS4DataSummaryState()`, which now owns the
  `state_manager$getState()` handoff plus `NULL` fallback that feeds
  `buildMetabQcS4DataSummaryUi()` before control returns to the metabolomics
  QC finalization module wrapper.
- The `data_summary` render path in
  `R/mod_metab_qc_s4.R`
  at line `649`
  now delegates through that helper instead of wiring the current-state fetch
  `tryCatch()` inline inside `mod_metab_qc_s4_server()`.
- The thirteenth bounded live seam now sits in
  `R/mod_metab_qc_s4.R`
  at line `504`
  as `getMetabQcS4AssayStatsState()`, which now owns the
  `state_manager$getState()` handoff plus `NULL` fallback that feeds
  `buildMetabQcS4AssayStatsDatatable()` before control returns to the
  metabolomics QC finalization module wrapper.
- The `assay_stats_table` render path in
  `R/mod_metab_qc_s4.R`
  at line `660`
  now delegates through that helper instead of wiring the current-state fetch
  `tryCatch()` inline inside `mod_metab_qc_s4_server()`.
- The third bounded live seam now sits in
  `R/mod_metab_qc_s4.R`
  at line `255`
  as `buildMetabQcS4AssayStatsDatatable()`, which now owns the invalid-state
  null exit, per-assay metabolite counting, numeric-sample counting,
  assay-level missingness calculation, and datatable assembly before control
  returns to the metabolomics QC finalization module wrapper.
- The `assay_stats_table` render path in
  `R/mod_metab_qc_s4.R`
  at line `660`
  now delegates through that helper instead of building the per-assay
  statistics datatable inline inside `mod_metab_qc_s4_server()`.
- The fourth bounded live seam now sits in
  `R/mod_metab_qc_s4.R`
  at line `322`
  as `renderMetabQcS4FilterPlot()`, which now owns the reactive plot
  `req()`, grob/gtable dispatch through `grid::grid.draw()`, and ggplot
  dispatch through `print()` before control returns to the metabolomics QC
  finalization module wrapper.
- The `filter_plot` render path in
  `R/mod_metab_qc_s4.R`
  at line `683`
  now delegates through that helper instead of inlining the QC-progress plot
  dispatch inside `mod_metab_qc_s4_server()`.
- The fifth bounded live seam now sits in
  `R/mod_metab_qc_s4.R`
  at line `344`
  as `buildMetabQcS4FinalizeResultsText()`, which now owns per-assay retained
  metabolite counting, processing-history text assembly, and finalization
  summary string creation before control returns to the metabolomics QC
  finalization module wrapper.
- The sixth bounded live seam now sits in
  `R/mod_metab_qc_s4.R`
  at line `460`
  as `reportMetabQcS4FinalizeSuccess()`, which now owns the
  `buildMetabQcS4FinalizeResultsText()` handoff, `finalize_results`
  render binding, completion log, and success toast before control returns to
  the metabolomics QC finalization module wrapper.
- The seventh bounded live seam now sits in
  `R/mod_metab_qc_s4.R`
  at line `398`
  as `saveMetabQcS4CompletedState()`, which now owns the
  `state_manager$saveState()` handoff, finalized-state name reuse, config
  wiring, and completion-description persistence before control returns to the
  metabolomics QC finalization module wrapper.
- The eighth bounded live seam now sits in
  `R/mod_metab_qc_s4.R`
  at line `443`
  as `updateMetabQcS4TrackingPlot()`, which now owns the
  `updateMetaboliteFiltering()` handoff, finalized-step name reuse,
  `filter_plot` reactive update, and null fallback before control returns to
  the metabolomics QC finalization module wrapper.
- The ninth bounded live seam now sits in
  `R/mod_metab_qc_s4.R`
  at line `418`
  as `completeMetabQcS4TabStatus()`, which now owns the
  `workflow_data$tab_status` list replacement, quality-control completion
  status update, and reactive assignment before control returns to the
  metabolomics QC finalization module wrapper.
- The tenth bounded live seam now sits in
  `R/mod_metab_qc_s4.R`
  at line `433`
  as `getMetabQcS4FinalizeHistory()`, which now owns the
  `state_manager$getHistory()` handoff that feeds
  `reportMetabQcS4FinalizeSuccess()` before control returns to the
  metabolomics QC finalization module wrapper.
- The fourteenth bounded live seam now sits in
  `R/mod_metab_qc_s4.R`
  at line `389`
  as `getMetabQcS4FinalizeState()`, which now owns the
  `state_manager$getState()` handoff that feeds
  `validateMetabQcS4FinalizeState()` before downstream finalize helper
  handoffs return control to the metabolomics QC finalization module wrapper.
- The fifteenth bounded live seam now sits in
  `R/mod_metab_qc_s4.R`
  at line `399`
  as `validateMetabQcS4FinalizeState()`, which now owns the
  `shiny::req(current_s4)` handoff plus the `MetaboliteAssayData` guard before
  downstream finalize helper handoffs return control to the metabolomics QC
  finalization module wrapper.
- The sixteenth bounded live seam now sits in
  `R/mod_metab_qc_s4.R`
  at line `519`
  as `reportMetabQcS4FinalizeError()`, which now owns the
  `paste("Error finalizing QC:", e$message)` message assembly,
  `logger::log_error()` handoff, and `shiny::showNotification(..., type =
  "error")` failure reporting before control returns to the metabolomics QC
  finalization module wrapper.
- The seventeenth bounded live seam now sits in
  `R/mod_metab_qc_s4.R`
  at line `538`
  as `runMetabQcS4FinalizeWorkflow()`, which now owns the remaining
  `finalize_qc` observer `tryCatch()` envelope and the sequencing across
  `getMetabQcS4FinalizeState()`, `validateMetabQcS4FinalizeState()`,
  `saveMetabQcS4CompletedState()`, `updateMetabQcS4TrackingPlot()`,
  `completeMetabQcS4TabStatus()`, `getMetabQcS4FinalizeHistory()`,
  `reportMetabQcS4FinalizeSuccess()`, and `reportMetabQcS4FinalizeError()`
  before control returns to the metabolomics QC finalization module wrapper.
- The finalize observer in
  `R/mod_metab_qc_s4.R`
  at line `671`
  now delegates through that helper instead of wiring the `tryCatch()`
  envelope and downstream finalize-helper chain inline inside
  `mod_metab_qc_s4_server()`.
- The eighteenth bounded live seam now sits in
  `R/mod_metab_qc_s4.R`
  at line `617`
  as `buildMetabQcS4AssayStatsRenderOutput()`, which now owns the remaining
  `assay_stats_table` render-path sequencing across
  `getMetabQcS4AssayStatsState()` and
  `buildMetabQcS4AssayStatsDatatable()` before control returns to the
  metabolomics QC finalization module wrapper.
- The `assay_stats_table` render path in
  `R/mod_metab_qc_s4.R`
  at line `660`
  now delegates through that helper instead of wiring the current-state fetch
  plus datatable-builder chain inline inside `mod_metab_qc_s4_server()`.
- The nineteenth bounded live seam now sits in
  `R/mod_metab_qc_s4.R`
  at line `242`
  as `buildMetabQcS4StateHistoryRenderOutput()`, which now owns the remaining
  `state_history` render-path sequencing across
  `getMetabQcS4StateHistory()` and
  `buildMetabQcS4StateHistoryUi()` before control returns to the metabolomics
  QC finalization module wrapper.
- The `state_history` render path in
  `R/mod_metab_qc_s4.R`
  at line `640`
  now delegates through that helper instead of wiring the history-fetch plus
  state-history UI-builder chain inline inside `mod_metab_qc_s4_server()`.
- The twentieth bounded live seam now sits in
  `R/mod_metab_qc_s4.R`
  at line `613`
  as `buildMetabQcS4DataSummaryRenderOutput()`, which now owns the remaining
  `data_summary` render-path sequencing across
  `getMetabQcS4DataSummaryState()` and
  `buildMetabQcS4DataSummaryUi()` before control returns to the metabolomics
  QC finalization module wrapper.
- The `data_summary` render path in
  `R/mod_metab_qc_s4.R`
  at line `675`
  now delegates through that helper instead of wiring the current-state fetch
  plus data-summary UI-builder chain inline inside `mod_metab_qc_s4_server()`.
- The twenty-first bounded live seam now sits in
  `R/mod_metab_qc_s4.R`
  at line `344`
  as `buildMetabQcS4FilterPlotRenderOutput()`, which now owns the remaining
  `filter_plot` render-path sequencing across
  `renderMetabQcS4FilterPlot()` before control returns to the metabolomics QC
  finalization module wrapper.
- The `filter_plot` render path in
  `R/mod_metab_qc_s4.R`
  at line `702`
  now delegates through that helper instead of wiring the render-helper
  handoff inline inside `mod_metab_qc_s4_server()`.
- Because this worktree does not include `renv/activate.R`,
  `tools/test_with_renv.R` is not available for this target and the focused
  gate reran green via direct `testthat::test_file()`.
- The first reviewed exact-source helper-surface manifest now exists at
  `tools/refactor/manifest-metab-qc-s4-module-wave1.yml`,
  targeting the frozen top-level QC S4 helper surface from
  `R/mod_metab_qc_s4.R`
  into
  `R/mod_metab_qc_s4_server_helpers.R`
  while the public wrapper stays live and unchanged.
- `Rscript tools/refactor/verify_refactor.R --manifest tools/refactor/manifest-metab-qc-s4-module-wave1.yml`
  passed against the current source tree before staging.
- The reviewed helper-surface wave has now been staged into
  `tools/refactor/staging/wave1_metabolomics_qc_s4_module_helper_surface/`,
  producing
  `R/mod_metab_qc_s4_server_helpers.R`
  plus
  `tools/refactor/collate-metab-qc-s4-module-wave1.txt`
  without rewriting live sources.
- Staged-output review confirmed that
  `R/mod_metab_qc_s4_server_helpers.R`
  contains all `21` frozen top-level helper definitions from
  `buildMetabQcS4DataSummaryUi()`
  through
  `buildMetabQcS4AssayStatsRenderOutput()`,
  and the staged helper file parses cleanly.
- The focused gate reran green after staging via direct `testthat::test_file()`
  because `tools/test_with_renv.R` still cannot run in this worktree without
  `renv/activate.R`.
- The reviewed exact-source helper-surface wave is now applied live via
  `tools/refactor/manifest-metab-qc-s4-module-wave1.yml`
  into
  `R/mod_metab_qc_s4_server_helpers.R`
  after moving all `21` frozen top-level helper definitions out of
  `R/mod_metab_qc_s4.R`
  while keeping
  `mod_metab_qc_s4_server()`
  as the public wrapper identity.
- `DESCRIPTION` `Collate:` now loads
  `mod_metab_qc_s4_server_helpers.R`
  before
  `mod_metab_qc_s4.R`
  so package load order matches the live helper surface.
- `Rscript tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-qc-s4-module-wave1.yml`
  passed after the live apply.
- The focused gate reran green again after the live apply via direct
  `testthat::test_file()` because `tools/test_with_renv.R` still cannot run in
  this worktree without `renv/activate.R`.
- The twenty-second bounded live seam now sits in
  `R/mod_metab_qc_s4.R`
  at line `115`
  as `runMetabQcS4ServerBody()`, which now owns the module callback's
  `filterPlot` reactive setup, the `state_history`, `data_summary`,
  `assay_stats_table`, and `filter_plot` render bindings, plus the
  `finalize_qc` observer registration before control returns to the public QC
  S4 server wrapper.
- `mod_metab_qc_s4_server()` in
  `R/mod_metab_qc_s4.R`
  at line `191`
  now delegates through `runMetabQcS4ServerBody()` instead of wiring the full
  `moduleServer()` callback inline inside the public wrapper.
- The focused gate in
  `tests/testthat/test-metab-01s-qc-s4-module-characterization.R`
  now additionally freezes the direct seam contract for
  `runMetabQcS4ServerBody()` plus the public-wrapper handoff into the
  server-body seam, and reran green again via direct `testthat::test_file()`
  because this worktree still does not include `renv/activate.R`.
- Post-checkpoint classification now records
  `R/mod_metab_qc_s4.R`
  at `202` lines with `3` top-level functions, a `54` line largest top-level
  function, and label `review`.
- The second reviewed exact-source QC S4 helper wave is now live on April 18,
  2026 via `tools/refactor/manifest-metab-qc-s4-module-wave2.yml`, with live
  collate artifact `tools/refactor/collate-metab-qc-s4-module-wave2.txt` and
  staged review artifacts retained under
  `tools/refactor/staging/wave2_metabolomics_qc_s4_module_server_body/`.
- The wave-2 live apply moved `runMetabQcS4ServerBody()` out of
  `R/mod_metab_qc_s4.R` and into `R/mod_metab_qc_s4_server_helpers.R` without
  changing the public `mod_metab_qc_s4_server()` wrapper signature or handoff.
- `Rscript tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-qc-s4-module-wave2.yml`
  passed after the live apply, and the focused gate reran green again via
  direct `testthat::test_file()` because this worktree still does not include
  `renv/activate.R`.
- Post-wave-2 classification now records `R/mod_metab_qc_s4.R` at `131` lines
  with `2` top-level functions, a `54` line largest top-level function, and
  label `review`.
- Treat `tools/refactor/HANDOVER-metab-qc-s4-module-seams.md` as complete and
  archived; future metabolomics QC S4 maintenance should run against the
  focused gate and touch `R/mod_metab_qc_s4_server_helpers.R` or the public
  wrapper only for real regressions, not for more structural reopening.
