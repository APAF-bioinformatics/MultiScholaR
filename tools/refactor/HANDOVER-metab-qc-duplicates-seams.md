# Metabolomics QC Duplicates Module Seam Map

## Goal

Introduce bounded top-level seams in R/mod_metab_qc_duplicates.R while keeping the live metabolomics duplicate-resolution module behavior frozen behind the existing public entry points.

## Current Position In The Flow

- target: `/home/doktersmol/Documents/MultiScholaR-metab-safe/R/mod_metab_qc_duplicates.R`
- classification: `direct-extraction-ready`
- checkpoint reached: `metab_qc_duplicates_wrapper_entrypoints_wave7_applied moved mod_metab_qc_duplicates_ui() and mod_metab_qc_duplicates_server() into dedicated live entrypoint files and reran the focused duplicate-module gate`
- next step: `Treat R/mod_metab_qc_duplicates.R as the completed breadcrumb stub for this manual target and archive this handover.`

## Existing Safety Net

- `tests/testthat/test-metab-01p-qc-duplicates-module-characterization.R`
- `tests/testthat/test-metab-01o-qc-module-characterization.R`

## Notes

Manual bucket 0 metabolomics QC duplicate-module stabilization target.

- This target previously had no active handover; this file now records the
  current metabolomics QC duplicate-module seam stop point.
- Classification refresh on April 18, 2026 now reports
  `R/mod_metab_qc_duplicates.R`
  at `49` lines; the current classifier output reports `0` top-level
  functions, `0` module servers, and label
  `direct-extraction-ready`.
- Historical seam-introduction bullets below preserve the earlier wrapper stop
  points; the current live helper and entrypoint ownership after the wave-1
  through wave-7 applies is summarized near the end of this handover.
- Classification refreshed on April 17, 2026 keeps
  `R/mod_metab_qc_duplicates.R`
  at `377` lines with `6` top-level functions, a `135` line largest
  top-level function, and labels
  `high-risk-wrapper`
  plus
  `needs-seam-introduction`.
- The focused gate now lives in
  `tests/testthat/test-metab-01p-qc-duplicates-module-characterization.R`
  and now freezes the direct seam contracts for
  `resolveMetabDuplicateAssayData()`
  plus
  `buildMetabDuplicateResolutionSummary()`
  plus
  `buildMetabDuplicateSummaryUi()`
  plus
  `buildMetabDuplicateTablesUi()`
  plus
  `registerMetabDuplicateTableRenderers()`
  plus
  `detectMetabDuplicateFeatures()`
  plus
  `revertMetabDuplicateResolution()`
  plus
  `renderMetabDuplicateFilterPlot()`
  plus the wrapper handoff from
  `mod_metab_qc_duplicates_server()`
  into the `detect_duplicates` observer,
  the duplicate-table renderer seam,
  into the summary render seam,
  the duplicate-tables render seam,
  the `resolve_duplicates` observer seam,
  the `revert_duplicates` observer seam,
  and the `filter_plot` render seam.
- The first bounded live seam now sits in
  `R/mod_metab_qc_duplicates.R`
  at line `96`
  as `resolveMetabDuplicateAssayData()`, which now owns the per-assay
  numeric-column detection, duplicate-resolution delegation, no-numeric-data
  warning path, and stats assembly before control returns to
  `mod_metab_qc_duplicates_server()`.
- The `resolve_duplicates` observer in
  `R/mod_metab_qc_duplicates.R`
  at line `441`
  now delegates through that helper instead of keeping the full per-assay loop
  inline inside the observer body.
- The second bounded live seam now sits in
  `R/mod_metab_qc_duplicates.R`
  at line `148`
  as `buildMetabDuplicateResolutionSummary()`, which now owns the
  per-assay summary-line assembly, total-removed counting, and resolution
  result-text construction before control returns to
  `mod_metab_qc_duplicates_server()`.
- The `resolve_duplicates` observer in
  `R/mod_metab_qc_duplicates.R`
  at line `441`
  now delegates through that helper instead of building the resolution
  summary text inline inside the observer body.
- The third bounded live seam now sits in
  `R/mod_metab_qc_duplicates.R`
  at line `191`
  as `buildMetabDuplicateSummaryUi()`, which now owns the duplicate-summary
  placeholder message, per-assay duplicate counts, and status-icon selection
  before control returns to `mod_metab_qc_duplicates_server()`.
- The `duplicate_summary` render path in
  `R/mod_metab_qc_duplicates.R`
  at line `420`
  now delegates through that helper instead of building the duplicate summary
  UI inline inside the wrapper.
- The fourth bounded live seam now sits in
  `R/mod_metab_qc_duplicates.R`
  at line `223`
  as `buildMetabDuplicateTablesUi()`, which now owns the null passthrough,
  zero-duplicate empty-state well, per-assay duplicate filtering, duplicate
  tab assembly, and table-output namespacing before control returns to
  `mod_metab_qc_duplicates_server()`.
- The `duplicate_tables` render path in
  `R/mod_metab_qc_duplicates.R`
  at line `425`
  now delegates through that helper instead of building the duplicate tables
  UI inline inside the wrapper.
- The fifth bounded live seam now sits in
  `R/mod_metab_qc_duplicates.R`
  at line `264`
  as `registerMetabDuplicateTableRenderers()`, which now owns the
  populated-assay filtering, sanitized duplicate-table output id construction,
  and DT renderer registration before control returns to
  `mod_metab_qc_duplicates_server()`.
- The duplicate-table registration observer in
  `R/mod_metab_qc_duplicates.R`
  at line `433`
  now delegates through that helper instead of registering the per-assay
  duplicate table renderers inline inside the wrapper.
- The sixth bounded live seam now sits in
  `R/mod_metab_qc_duplicates.R`
  at line `304`
  as `detectMetabDuplicateFeatures()`, which now owns the state-manager
  presence check, metabolomics S4 retrieval and class validation,
  duplicate-detection delegation, and total-duplicate counting before control
  returns to `mod_metab_qc_duplicates_server()`.
- The `detect_duplicates` observer in
  `R/mod_metab_qc_duplicates.R`
  at line `391`
  now delegates through that helper instead of retrieving the current state,
  validating the S4 class, invoking `findMetabDuplicateFeatureIDs()`, and
  counting duplicates inline inside the observer body.
- The seventh bounded live seam now sits in
  `R/mod_metab_qc_duplicates.R`
  at line `335`
  as `revertMetabDuplicateResolution()`, which now owns the state-manager
  presence check, history lookup, previous-state selection, revert
  delegation, and revert result-text assembly before control returns to
  `mod_metab_qc_duplicates_server()`.
- The `revert_duplicates` observer in
  `R/mod_metab_qc_duplicates.R`
  at line `520`
  now delegates through that helper instead of retrieving history,
  selecting the previous state, reverting inline, and building the revert
  result text inside the observer body.
- The eighth bounded live seam now sits in
  `R/mod_metab_qc_duplicates.R`
  at line `357`
  as `renderMetabDuplicateFilterPlot()`, which now owns the filter-plot
  reactive lookup, req gating, grob/gtable draw dispatch, and ggplot print
  dispatch before control returns to `mod_metab_qc_duplicates_server()`.
- The `filter_plot` render path in
  `R/mod_metab_qc_duplicates.R`
  at line `542`
  now delegates through that helper instead of reading the reactive value and
  branching on plot classes inline inside the wrapper.
- The existing parent QC wrapper gate in
  `tests/testthat/test-metab-01o-qc-module-characterization.R`
  remained covered from the prior checkpoint, preserving the upstream wrapper
  contract that still registers the duplicate-resolution submodule.
- The first exact-source duplicate-module wave now lives in
  `tools/refactor/manifest-metab-qc-duplicates-wave1.yml`,
  with reviewed staging artifacts in
  `tools/refactor/staging/wave1_metabolomics_qc_duplicates_ui_render_helpers/R/mod_metab_qc_duplicates_ui_helpers.R`,
  `tools/refactor/staging/wave1_metabolomics_qc_duplicates_ui_render_helpers/R/mod_metab_qc_duplicates_render_helpers.R`,
  and
  `tools/refactor/staging/wave1_metabolomics_qc_duplicates_ui_render_helpers/collate-metab-qc-duplicates-wave1.txt`.
- The staged wave covers
  `buildMetabDuplicateSummaryUi()`,
  `buildMetabDuplicateTablesUi()`,
  `registerMetabDuplicateTableRenderers()`, and
  `renderMetabDuplicateFilterPlot()`
  from
  `R/mod_metab_qc_duplicates.R`
  without rewriting live sources.
- `tools/refactor/verify_refactor.R` passed for
  `tools/refactor/manifest-metab-qc-duplicates-wave1.yml`
  before staging the reviewed helper targets.
- The focused duplicate-module gate loader in
  `tests/testthat/test-metab-01p-qc-duplicates-module-characterization.R`
  now probes
  `R/mod_metab_qc_duplicates_ui_helpers.R`
  and
  `R/mod_metab_qc_duplicates_render_helpers.R`
  ahead of
  `R/mod_metab_qc_duplicates.R`,
  so the characterization surface survives the live apply boundary.
- The focused gate reran green after the staged wave via direct
  `testthat::test_file()` invocation because this worktree does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- Live classification remains
  `high-risk-wrapper`
  /
  `needs-seam-introduction`
  at `548` lines with `10` top-level functions and a `167` line largest
  top-level function because the reviewed wave is staged only and
  `R/mod_metab_qc_duplicates.R`
  remains the live wrapper.
- The next safe stop point is a reviewed live apply of wave 1 into
  `R/mod_metab_qc_duplicates_ui_helpers.R`
  and
  `R/mod_metab_qc_duplicates_render_helpers.R`,
  with `DESCRIPTION` `Collate:` loading both helper files ahead of
  `R/mod_metab_qc_duplicates.R`,
  then `tools/refactor/check_wave_apply.R`
  and the focused gate rerun before any further wrapper extraction.
- Wave 1 is now live in
  `R/mod_metab_qc_duplicates_ui_helpers.R`
  and
  `R/mod_metab_qc_duplicates_render_helpers.R`
  after applying
  `tools/refactor/manifest-metab-qc-duplicates-wave1.yml`
  through the repo refactor engine.
- `tools/refactor/collate-metab-qc-duplicates-wave1.txt` now records the
  applied collate output, and `DESCRIPTION` `Collate:` now loads
  `mod_metab_qc_duplicates_ui_helpers.R`
  plus
  `mod_metab_qc_duplicates_render_helpers.R`
  ahead of
  `mod_metab_qc_duplicates.R`.
- `tools/refactor/check_wave_apply.R` passed during the live apply checkpoint
  for
  `tools/refactor/manifest-metab-qc-duplicates-wave1.yml`.
- The focused duplicate-module gate reran green again after the live apply via
  direct `testthat::test_file()` invocation because this worktree still does
  not include `renv/activate.R` for `tools/test_with_renv.R`.
- Post-apply classification now records
  `R/mod_metab_qc_duplicates.R`
  at `420` lines with `6` top-level functions, a `167` line largest
  top-level function, and labels
  `high-risk-wrapper`
  plus
  `needs-seam-introduction`.
- The second exact-source duplicate-module wave now lives in
  `tools/refactor/manifest-metab-qc-duplicates-wave2.yml`,
  with reviewed staging artifacts in
  `tools/refactor/staging/wave2_metabolomics_qc_duplicates_server_helpers/R/mod_metab_qc_duplicates_server_helpers.R`
  and
  `tools/refactor/staging/wave2_metabolomics_qc_duplicates_server_helpers/collate-metab-qc-duplicates-wave2.txt`.
- The staged wave covers
  `resolveMetabDuplicateAssayData()`,
  `buildMetabDuplicateResolutionSummary()`,
  `detectMetabDuplicateFeatures()`, and
  `revertMetabDuplicateResolution()`
  from
  `R/mod_metab_qc_duplicates.R`
  into
  `R/mod_metab_qc_duplicates_server_helpers.R`
  without hand-rewriting helper bodies.
- `tools/refactor/verify_refactor.R` passed for
  `tools/refactor/manifest-metab-qc-duplicates-wave2.yml`
  before staging and before the live apply.
- The focused duplicate-module gate loader in
  `tests/testthat/test-metab-01p-qc-duplicates-module-characterization.R`
  now probes
  `R/mod_metab_qc_duplicates_ui_helpers.R`,
  `R/mod_metab_qc_duplicates_render_helpers.R`,
  and
  `R/mod_metab_qc_duplicates_server_helpers.R`
  ahead of
  `R/mod_metab_qc_duplicates.R`,
  so the characterization surface survives the second live apply boundary.
- Wave 2 is now live in
  `R/mod_metab_qc_duplicates_server_helpers.R`
  after applying
  `tools/refactor/manifest-metab-qc-duplicates-wave2.yml`
  through the repo refactor engine.
- `DESCRIPTION` `Collate:` now loads
  `mod_metab_qc_duplicates_ui_helpers.R`,
  `mod_metab_qc_duplicates_render_helpers.R`,
  and
  `mod_metab_qc_duplicates_server_helpers.R`
  ahead of
  `mod_metab_qc_duplicates.R`,
  with the applied load-order artifact recorded in
  `tools/refactor/collate-metab-qc-duplicates-wave2.txt`.
- `tools/refactor/check_wave_apply.R` passed during the live apply checkpoint
  for
  `tools/refactor/manifest-metab-qc-duplicates-wave2.yml`.
- The focused duplicate-module gate reran green again after the live apply via
  direct `testthat::test_file()` invocation because this worktree still does
  not include `renv/activate.R` for `tools/test_with_renv.R`.
- Post-apply classification now records
  `R/mod_metab_qc_duplicates.R`
  at `276` lines with `2` top-level functions, a `167` line largest
  top-level function, and labels
  `high-risk-wrapper`
  plus
  `needs-seam-introduction`.
- The exact-source top-level extraction surface is now exhausted in
  `R/mod_metab_qc_duplicates.R`;
  the next safe stop point is one bounded in-place observer seam inside
  `mod_metab_qc_duplicates_server()`,
  starting with the `detect_duplicates` observer notification/logging path,
  then rerunning the focused duplicate-module gate before any further
  observer extraction.
- The ninth bounded live seam now sits in
  `R/mod_metab_qc_duplicates.R`
  at line `96`
  as `reportMetabDuplicateDetection()`, which now owns the
  `detect_duplicates` observer success-path log text, notification text, and
  notification-type selection before control returns to
  `mod_metab_qc_duplicates_server()`.
- The `detect_duplicates` observer in
  `R/mod_metab_qc_duplicates.R`
  at line `146`
  now delegates through that helper instead of logging the duplicate count
  and choosing the success notification type inline inside the observer body.
- The focused duplicate-module gate in
  `tests/testthat/test-metab-01p-qc-duplicates-module-characterization.R`
  now additionally freezes the direct seam contract for
  `reportMetabDuplicateDetection()`
  plus the wrapper handoff from
  `mod_metab_qc_duplicates_server()`
  into the `detect_duplicates` observer reporting seam.
- The focused duplicate-module gate reran green again after the seam via
  direct `testthat::test_file()` invocation because this worktree still does
  not include `renv/activate.R` for `tools/test_with_renv.R`.
- Post-checkpoint classification now records
  `R/mod_metab_qc_duplicates.R`
  at `292` lines with `3` top-level functions, a `156` line largest
  top-level function, and labels
  `high-risk-wrapper`
  plus
  `needs-seam-introduction`.
- The next safe stop point is one bounded in-place observer seam inside
  `mod_metab_qc_duplicates_server()`,
  continuing within the `detect_duplicates` observer by extracting the error
  notification/logging path, then rerunning the focused duplicate-module gate
  before any further observer extraction.
- The tenth bounded live seam now sits in
  `R/mod_metab_qc_duplicates.R`
  at line `123`
  as `reportMetabDuplicateDetectionError()`, which now owns the
  `detect_duplicates` observer error-path message assembly, error logging, and
  error notification dispatch before control returns to
  `mod_metab_qc_duplicates_server()`.
- The `detect_duplicates` observer in
  `R/mod_metab_qc_duplicates.R`
  at line `165`
  now delegates through that helper instead of constructing the error message
  inline and calling `logger::log_error()` plus
  `shiny::showNotification()` directly inside the observer body.
- The focused duplicate-module gate in
  `tests/testthat/test-metab-01p-qc-duplicates-module-characterization.R`
  now additionally freezes the direct seam contract for
  `reportMetabDuplicateDetectionError()`
  plus the wrapper handoff from
  `mod_metab_qc_duplicates_server()`
  into the `detect_duplicates` observer error seam.
- The focused duplicate-module gate reran green again after the seam via
  direct `testthat::test_file()` invocation because this worktree still does
  not include `renv/activate.R` for `tools/test_with_renv.R`.
- Post-checkpoint classification now records
  `R/mod_metab_qc_duplicates.R`
  at `311` lines with `4` top-level functions, a `156` line largest
  top-level function, and labels
  `high-risk-wrapper`
  plus
  `needs-seam-introduction`.
- The next safe stop point is one bounded in-place observer shell seam inside
  `mod_metab_qc_duplicates_server()`,
  extracting the remaining `detect_duplicates` observer orchestration around
  duplicate-info state update plus success/error helper dispatch, then
  rerunning the focused duplicate-module gate before any further observer
  extraction.
- The eleventh bounded live seam now sits in
  `R/mod_metab_qc_duplicates.R`
  at line `142`
  as `runMetabDuplicateDetectionObserverShell()`, which now owns the
  remaining `detect_duplicates` observer shell orchestration around
  duplicate-info state update plus success/error helper dispatch before
  control returns to `mod_metab_qc_duplicates_server()`.
- The `detect_duplicates` observer in
  `R/mod_metab_qc_duplicates.R`
  at line `198`
  now delegates through that helper instead of coordinating the duplicate-info
  reactive write plus the success and error reporting helpers inline inside
  the observer body.
- The focused duplicate-module gate in
  `tests/testthat/test-metab-01p-qc-duplicates-module-characterization.R`
  now additionally freezes the direct seam contract for
  `runMetabDuplicateDetectionObserverShell()`
  plus the wrapper handoff from
  `mod_metab_qc_duplicates_server()`
  into the `detect_duplicates` observer shell seam.
- The focused duplicate-module gate reran green again after the seam via
  direct `testthat::test_file()` invocation because this worktree still does
  not include `renv/activate.R` for `tools/test_with_renv.R`.
- Post-checkpoint classification now records
  `R/mod_metab_qc_duplicates.R`
  at `335` lines with `5` top-level functions, a `146` line largest
  top-level function, and labels
  `high-risk-wrapper`
  plus
  `needs-seam-introduction`.
- The next safe stop point is one bounded in-place observer shell seam inside
  `mod_metab_qc_duplicates_server()`,
  shifting to the `resolve_duplicates` observer by extracting the
  success/error dispatch around resolution-results rendering, duplicate-info
  cleanup, and working-notification teardown, then rerunning the focused
  duplicate-module gate before any further observer extraction.
- The twelfth bounded live seam now sits in
  `R/mod_metab_qc_duplicates.R`
  at line `176`
  as `runMetabDuplicateResolutionObserverShell()`, which now owns the
  `resolve_duplicates` observer success/error dispatch around
  `resolution_results` rendering, duplicate-info cleanup, success/error
  notification dispatch, and working-notification teardown before control
  returns to `mod_metab_qc_duplicates_server()`.
- The `resolve_duplicates` observer in
  `R/mod_metab_qc_duplicates.R`
  at line `280`
  now delegates through that helper instead of rendering the resolution text,
  clearing duplicate-info state, and managing the resolve success/error
  notifications inline inside the observer body.
- The focused duplicate-module gate in
  `tests/testthat/test-metab-01p-qc-duplicates-module-characterization.R`
  now additionally freezes the direct seam contract for
  `runMetabDuplicateResolutionObserverShell()`
  plus the wrapper handoff from
  `mod_metab_qc_duplicates_server()`
  into the `resolve_duplicates` observer shell seam.
- The focused duplicate-module gate reran green again after the seam via
  direct `testthat::test_file()` invocation because this worktree still does
  not include `renv/activate.R` for `tools/test_with_renv.R`.
- Post-checkpoint classification now records
  `R/mod_metab_qc_duplicates.R`
  at `377` lines with `6` top-level functions, a `135` line largest
  top-level function, and labels
  `high-risk-wrapper`
  plus
  `needs-seam-introduction`.
- The next safe stop point is one bounded in-place observer shell seam inside
  `mod_metab_qc_duplicates_server()`,
  shifting to the `revert_duplicates` observer by extracting the result-text
  rendering, reactive reset, and success/error notification dispatch, then
  rerunning the focused duplicate-module gate before any further observer
  extraction.
- The thirteenth bounded live seam now sits in
  `R/mod_metab_qc_duplicates.R`
  at line `229`
  as `runMetabDuplicateRevertObserverShell()`, which now owns the
  `revert_duplicates` observer success/error dispatch around
  `resolution_results` rendering, reactive reset, and success/error
  notification dispatch before control returns to
  `mod_metab_qc_duplicates_server()`.
- The `revert_duplicates` observer in
  `R/mod_metab_qc_duplicates.R`
  at line `396`
  now delegates through that helper instead of rendering the revert result,
  clearing the duplicate-module reactive state, and coordinating the revert
  success/error notifications inline inside the observer body.
- The focused duplicate-module gate in
  `tests/testthat/test-metab-01p-qc-duplicates-module-characterization.R`
  now additionally freezes the direct seam contract for
  `runMetabDuplicateRevertObserverShell()`
  plus the wrapper handoff from
  `mod_metab_qc_duplicates_server()`
  into the `revert_duplicates` observer shell seam.
- The focused duplicate-module gate reran green again after the seam via
  direct `testthat::test_file()` invocation because this worktree still does
  not include `renv/activate.R` for `tools/test_with_renv.R`.
- Post-checkpoint classification now records
  `R/mod_metab_qc_duplicates.R`
  at `416` lines with `7` top-level functions, a `128` line largest
  top-level function, and labels
  `high-risk-wrapper`
  plus
  `needs-seam-introduction`.
- The next safe stop point is one reviewed exact-source helper wave for the
  duplicate observer shell cluster in
  `R/mod_metab_qc_duplicates.R`,
  staging
  `reportMetabDuplicateDetection()`,
  `reportMetabDuplicateDetectionError()`,
  `runMetabDuplicateDetectionObserverShell()`,
  `runMetabDuplicateResolutionObserverShell()`,
  and
  `runMetabDuplicateRevertObserverShell()`
  into
  `R/mod_metab_qc_duplicates_server_helpers.R`,
  then verifying/applying that wave and rerunning the focused duplicate-module
  gate before any further wrapper extraction.
- The third exact-source duplicate-module wave now lives in
  `tools/refactor/manifest-metab-qc-duplicates-wave3.yml`,
  with reviewed staging artifacts in
  `tools/refactor/staging/wave3_metabolomics_qc_duplicates_observer_shell_helpers/R/mod_metab_qc_duplicates_server_helpers.R`
  and
  `tools/refactor/staging/wave3_metabolomics_qc_duplicates_observer_shell_helpers/collate-metab-qc-duplicates-wave3.txt`.
- The wave covers
  `reportMetabDuplicateDetection()`,
  `reportMetabDuplicateDetectionError()`,
  `runMetabDuplicateDetectionObserverShell()`,
  `runMetabDuplicateResolutionObserverShell()`,
  and
  `runMetabDuplicateRevertObserverShell()`
  from
  `R/mod_metab_qc_duplicates.R`
  into
  `R/mod_metab_qc_duplicates_server_helpers.R`
  without hand-rewriting helper bodies.
- `tools/refactor/verify_refactor.R` and
  `tools/refactor/check_wave_apply.R` both passed for
  `tools/refactor/manifest-metab-qc-duplicates-wave3.yml`.
- Wave 3 is now live in
  `R/mod_metab_qc_duplicates_server_helpers.R`
  after applying
  `tools/refactor/manifest-metab-qc-duplicates-wave3.yml`
  through the repo refactor engine.
- `DESCRIPTION` `Collate:` continues to load
  `mod_metab_qc_duplicates_ui_helpers.R`,
  `mod_metab_qc_duplicates_render_helpers.R`,
  and
  `mod_metab_qc_duplicates_server_helpers.R`
  ahead of
  `mod_metab_qc_duplicates.R`,
  with the matching load-order artifact recorded in
  `tools/refactor/collate-metab-qc-duplicates-wave3.txt`.
- The focused duplicate-module gate reran green again after the live apply via
  direct `testthat::test_file()` invocation because this worktree still does
  not include `renv/activate.R` for `tools/test_with_renv.R`.
- Post-apply classification now records
  `R/mod_metab_qc_duplicates.R`
  at `242` lines with `2` top-level functions, a `128` line largest
  top-level function, and labels
  `high-risk-wrapper`
  plus
  `needs-seam-introduction`.
- The exact-source top-level extraction surface is now exhausted again in
  `R/mod_metab_qc_duplicates.R`;
  the next safe stop point is one bounded in-place seam inside
  `mod_metab_qc_duplicates_server()`,
  starting with the `resolve_duplicates` observer shell's current-state
  preflight plus duplicate-resolution dispatch segment before the
  QC-update/save-state tail, then rerunning the focused duplicate-module gate
  before any further wrapper extraction.
- The next bounded live seam now sits in
  `R/mod_metab_qc_duplicates.R`
  at line `96`
  as `prepareMetabDuplicateResolutionState()`, which now owns the
  `resolve_duplicates` observer shell's current-state lookup, `req` gating,
  metabolomics S4 class validation, duplicate-resolution dispatch, and
  resolved-assay writeback before control returns to
  `mod_metab_qc_duplicates_server()`.
- The `resolve_duplicates` observer in
  `R/mod_metab_qc_duplicates.R`
  at line `190`
  now delegates through that helper instead of retrieving the current state,
  validating the S4 class, invoking `resolveMetabDuplicateAssayData()`, and
  writing the resolved assay list inline inside the observer shell.
- The focused duplicate-module gate in
  `tests/testthat/test-metab-01p-qc-duplicates-module-characterization.R`
  now additionally freezes the direct seam contract for
  `prepareMetabDuplicateResolutionState()`
  plus the wrapper handoff through the `resolve_duplicates` observer
  preflight seam, and reran green again via direct `testthat::test_file()`
  because this worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Post-checkpoint classification now records
  `R/mod_metab_qc_duplicates.R`
  at `259` lines with `3` top-level functions, a `118` line largest
  top-level function, and labels
  `high-risk-wrapper`
  plus
  `needs-seam-introduction`.
- The next safe stop point is one bounded in-place seam inside
  `mod_metab_qc_duplicates_server()`, shifting to the remaining
  `resolve_duplicates` observer shell tail around the `resolution_stats()`
  reactive write, `saveState()` call, `updateMetaboliteFiltering()` refresh,
  and `filter_plot()` update before the final summary return, then rerunning
  the focused duplicate-module gate before any further wrapper extraction.
- The next bounded live seam now sits in
  `R/mod_metab_qc_duplicates.R`
  at line `123`
  as `applyMetabDuplicateResolutionState()`, which now owns the
  `resolution_stats()` reactive write, workflow `saveState()` call,
  `updateMetaboliteFiltering()` refresh, warning-backed QC plot fallback, and
  `filter_plot()` reactive update before control returns to
  `mod_metab_qc_duplicates_server()`.
- The `resolve_duplicates` observer in
  `R/mod_metab_qc_duplicates.R`
  at line `223`
  now delegates through both
  `prepareMetabDuplicateResolutionState()`
  and
  `applyMetabDuplicateResolutionState()`
  instead of coordinating the resolved-state write, QC refresh, and
  workflow-state persistence inline inside the observer shell.
- The focused duplicate-module gate in
  `tests/testthat/test-metab-01p-qc-duplicates-module-characterization.R`
  now additionally freezes the direct seam contract for
  `applyMetabDuplicateResolutionState()`
  plus the wrapper handoff through the `resolve_duplicates` observer
  apply-state seam, and reran green again via direct
  `testthat::test_file()` because this worktree still does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- Post-checkpoint classification now records
  `R/mod_metab_qc_duplicates.R`
  at `286` lines with `4` top-level functions, a `102` line largest
  top-level function, and labels
  `high-risk-wrapper`
  plus
  `needs-seam-introduction`.
- The next safe stop point is one reviewed exact-source helper wave for
  `prepareMetabDuplicateResolutionState()`
  plus
  `applyMetabDuplicateResolutionState()`
  from
  `R/mod_metab_qc_duplicates.R`
  into
  `R/mod_metab_qc_duplicates_server_helpers.R`,
  then manifest verification/live apply and the focused duplicate-module gate
  rerun before any further wrapper extraction.
- Wave 4 manifest now lives in
  `tools/refactor/manifest-metab-qc-duplicates-wave4.yml`
  and stages the exact-source duplicate-resolution state helper pair from
  `R/mod_metab_qc_duplicates.R`
  into
  `R/mod_metab_qc_duplicates_server_helpers.R`
  without rewriting helper bodies during review.
- The reviewed staging artifacts for wave 4 now live in
  `tools/refactor/staging/wave4_metabolomics_qc_duplicates_resolution_state_helpers/R/mod_metab_qc_duplicates_server_helpers.R`
  and
  `tools/refactor/staging/wave4_metabolomics_qc_duplicates_resolution_state_helpers/collate-metab-qc-duplicates-wave4.txt`.
- Wave 4 is now live after applying
  `tools/refactor/manifest-metab-qc-duplicates-wave4.yml`
  through the repo refactor engine.
- `tools/refactor/check_wave_apply.R` passed during the wave 4 live apply,
  and
  `tools/refactor/collate-metab-qc-duplicates-wave4.txt`
  now records the applied explicit collate order for the duplicate-module
  helper set.
- The focused duplicate-module gate in
  `tests/testthat/test-metab-01p-qc-duplicates-module-characterization.R`
  continued to load
  `R/mod_metab_qc_duplicates_server_helpers.R`
  before
  `R/mod_metab_qc_duplicates.R`,
  so the characterization surface stayed frozen across the live apply
  boundary, and the gate reran green again via direct `testthat::test_file()`
  because this worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Post-apply classification now records
  `R/mod_metab_qc_duplicates.R`
  at `218` lines with `2` top-level functions, a `102` line largest
  top-level function, and labels
  `high-risk-wrapper`
  plus
  `needs-seam-introduction`.
- The next safe stop point is one bounded in-place seam inside
  `mod_metab_qc_duplicates_server()`
  for the remaining `resolve_duplicates` dispatch closure that composes
  `prepareMetabDuplicateResolutionState()`,
  `applyMetabDuplicateResolutionState()`,
  and
  `buildMetabDuplicateResolutionSummary()`
  before control returns to
  `runMetabDuplicateResolutionObserverShell()`,
  then the focused duplicate-module gate rerun before any further wrapper
  extraction.
- The next bounded live seam now sits in
  `R/mod_metab_qc_duplicates.R`
  as `runMetabDuplicateResolutionWorkflow()`, which now owns the remaining
  `resolve_duplicates` dispatch closure that composes
  `prepareMetabDuplicateResolutionState()`,
  `applyMetabDuplicateResolutionState()`,
  and
  `buildMetabDuplicateResolutionSummary()`
  before control returns to
  `runMetabDuplicateResolutionObserverShell()`.
- The `resolve_duplicates` observer in
  `R/mod_metab_qc_duplicates.R`
  now delegates through
  `runMetabDuplicateResolutionWorkflow()`
  instead of composing the resolution preflight, apply-state, and final
  summary flow inline inside
  `mod_metab_qc_duplicates_server()`.
- The focused duplicate-module gate in
  `tests/testthat/test-metab-01p-qc-duplicates-module-characterization.R`
  now additionally freezes the direct seam contract for
  `runMetabDuplicateResolutionWorkflow()`
  plus the wrapper handoff through the `resolve_duplicates` observer workflow
  seam, and reran green again via direct `testthat::test_file()` because this
  worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Post-checkpoint classification now records
  `R/mod_metab_qc_duplicates.R`
  at `235` lines with `3` top-level functions, a `84` line largest top-level
  function, and labels
  `high-risk-wrapper`
  plus
  `needs-seam-introduction`.
- The next safe stop point is one reviewed exact-source helper wave for
  `runMetabDuplicateResolutionWorkflow()`
  from
  `R/mod_metab_qc_duplicates.R`
  into
  `R/mod_metab_qc_duplicates_server_helpers.R`,
  then manifest verification/live apply and the focused duplicate-module gate
  rerun before any further wrapper extraction.
- April 18, 2026 wave 5 manifest now lives in
  `tools/refactor/manifest-metab-qc-duplicates-wave5.yml`
  and stages the exact-source duplicate-resolution workflow helper from
  `R/mod_metab_qc_duplicates.R`
  into
  `R/mod_metab_qc_duplicates_server_helpers.R`
  without hand-rewriting the helper body during review.
- The reviewed staging artifacts for wave 5 now live in
  `tools/refactor/staging/wave5_metabolomics_qc_duplicates_resolution_workflow_helper/R/mod_metab_qc_duplicates_server_helpers.R`
  and
  `tools/refactor/staging/wave5_metabolomics_qc_duplicates_resolution_workflow_helper/collate-metab-qc-duplicates-wave5.txt`.
- `tools/refactor/verify_refactor.R` passed for
  `tools/refactor/manifest-metab-qc-duplicates-wave5.yml`
  before staging the reviewed workflow helper target.
- The focused duplicate-module gate in
  `tests/testthat/test-metab-01p-qc-duplicates-module-characterization.R`
  reran green again via direct `testthat::test_file()` because this worktree
  still does not include `renv/activate.R` for `tools/test_with_renv.R`.
- Live classification remains
  `R/mod_metab_qc_duplicates.R`
  at `235` lines with `3` top-level functions, a `84` line largest top-level
  function, and labels
  `high-risk-wrapper`
  plus
  `needs-seam-introduction`
  because the reviewed wave is staged only and the live wrapper is unchanged.
- The next safe stop point is a reviewed live apply of wave 5 from
  `tools/refactor/manifest-metab-qc-duplicates-wave5.yml`,
  materializing
  `runMetabDuplicateResolutionWorkflow()`
  in
  `R/mod_metab_qc_duplicates_server_helpers.R`,
  recording the applied explicit collate order, then running
  `tools/refactor/check_wave_apply.R`
  and the focused duplicate-module gate before any further wrapper
  extraction.
- April 18, 2026 wave 5 now applies live via
  `tools/refactor/manifest-metab-qc-duplicates-wave5.yml`,
  materializing
  `runMetabDuplicateResolutionWorkflow()`
  in
  `R/mod_metab_qc_duplicates_server_helpers.R`
  while removing that exact helper body from the live wrapper without
  hand-rewriting it during the apply step.
- `tools/refactor/check_wave_apply.R` passed for
  `tools/refactor/manifest-metab-qc-duplicates-wave5.yml`,
  and
  `tools/refactor/collate-metab-qc-duplicates-wave5.txt`
  now records the applied explicit collate order for the duplicate-module
  helper surface.
- The focused duplicate-module gate in
  `tests/testthat/test-metab-01p-qc-duplicates-module-characterization.R`
  reran green again after the live apply via direct `testthat::test_file()`
  because this worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Post-apply classification now records
  `R/mod_metab_qc_duplicates.R`
  at `201` lines with `2` top-level functions, a `84` line largest
  top-level function, and labels
  `high-risk-wrapper`
  plus
  `needs-seam-introduction`.
- The next safe stop point is one bounded in-place seam inside
  `mod_metab_qc_duplicates_server()`
  for the remaining `resolve_duplicates` observer wrapper that still owns
  the `shiny::req()` guard, working-notification setup, and
  `runMetabDuplicateResolutionObserverShell()` handoff around
  `runMetabDuplicateResolutionWorkflow()`,
  then the focused duplicate-module gate rerun before any further wrapper
  extraction.
- The next bounded live seam now sits in
  `R/mod_metab_qc_duplicates.R`
  at line `96`
  as `runMetabDuplicateResolutionObserver()`, which now owns the
  `shiny::req()` guard, working-notification setup, and
  `runMetabDuplicateResolutionObserverShell()` handoff around
  `runMetabDuplicateResolutionWorkflow()`
  before control returns to `mod_metab_qc_duplicates_server()`.
- The `resolve_duplicates` observer in
  `R/mod_metab_qc_duplicates.R`
  at line `190`
  now delegates through
  `runMetabDuplicateResolutionObserver()`
  instead of keeping that notification-shell and workflow-handoff logic
  inline inside `mod_metab_qc_duplicates_server()`.
- The focused duplicate-module gate in
  `tests/testthat/test-metab-01p-qc-duplicates-module-characterization.R`
  now additionally freezes the direct seam contract for
  `runMetabDuplicateResolutionObserver()`
  plus the wrapper handoff through the `resolve_duplicates` observer seam,
  and reran green again via direct `testthat::test_file()` because this
  worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Post-checkpoint classification now records
  `R/mod_metab_qc_duplicates.R`
  at `223` lines with `4` top-level functions, a `72` line largest
  top-level function, and labels
  `high-risk-wrapper`
  plus
  `needs-seam-introduction`.
- The next safe stop point is one reviewed exact-source helper wave for
  `runMetabDuplicateResolutionObserver()`
  from
  `R/mod_metab_qc_duplicates.R`
  into
  `R/mod_metab_qc_duplicates_server_helpers.R`,
  then manifest verification/live apply and the focused duplicate-module gate
  rerun before any further wrapper extraction.
- April 18, 2026 wave 6 manifest now lives in
  `tools/refactor/manifest-metab-qc-duplicates-wave6.yml`
  and stages the exact-source duplicate-resolution observer wrapper helper
  from
  `R/mod_metab_qc_duplicates.R`
  into
  `R/mod_metab_qc_duplicates_server_helpers.R`
  without hand-rewriting the helper body during review.
- Reviewed wave 6 staging artifacts now live in
  `tools/refactor/staging/wave6_metabolomics_qc_duplicates_resolution_observer_wrapper_helper/R/mod_metab_qc_duplicates_server_helpers.R`
  and
  `tools/refactor/staging/wave6_metabolomics_qc_duplicates_resolution_observer_wrapper_helper/collate-metab-qc-duplicates-wave6.txt`.
- April 18, 2026 wave 6 now applies live via
  `tools/refactor/manifest-metab-qc-duplicates-wave6.yml`,
  materializing
  `runMetabDuplicateResolutionObserver()`
  in
  `R/mod_metab_qc_duplicates_server_helpers.R`
  while removing that exact helper body from the live wrapper without
  hand-rewriting it during the apply step.
- `tools/refactor/verify_refactor.R` and
  `tools/refactor/check_wave_apply.R` both passed for
  `tools/refactor/manifest-metab-qc-duplicates-wave6.yml`,
  and
  `tools/refactor/collate-metab-qc-duplicates-wave6.txt`
  now records the applied explicit collate order for the duplicate-module
  helper surface.
- The focused duplicate-module gate in
  `tests/testthat/test-metab-01p-qc-duplicates-module-characterization.R`
  reran green again after the live apply via direct `testthat::test_file()`
  because this worktree still does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Post-apply classification now records
  `R/mod_metab_qc_duplicates.R`
  at `190` lines with `2` top-level functions, a `72` line largest
  top-level function, and labels
  `high-risk-wrapper`
  plus
  `needs-seam-introduction`.
- Despite the conservative classifier labels, the manual metabolomics QC
  duplicate-module target is now within the stabilization budget and no
  longer blocks bucket 0; any later extraction from the public wrapper
  entrypoints is optional cleanup rather than required stabilization work.
- April 18, 2026 wave 7 manifest now lives in
  `tools/refactor/manifest-metab-qc-duplicates-wave7.yml`
  and stages the remaining public duplicate-module wrapper entrypoints from
  `R/mod_metab_qc_duplicates.R`
  into
  `R/mod_metab_qc_duplicates_ui.R`
  and
  `R/mod_metab_qc_duplicates_server.R`
  without hand-rewriting either entrypoint during review.
- Reviewed wave 7 staging artifacts now live in
  `tools/refactor/staging/wave7_metabolomics_qc_duplicates_wrapper_entrypoints/R/mod_metab_qc_duplicates_ui.R`,
  `tools/refactor/staging/wave7_metabolomics_qc_duplicates_wrapper_entrypoints/R/mod_metab_qc_duplicates_server.R`,
  and
  `tools/refactor/staging/wave7_metabolomics_qc_duplicates_wrapper_entrypoints/collate-metab-qc-duplicates-wave7.txt`.
- April 18, 2026 wave 7 now applies live via
  `tools/refactor/manifest-metab-qc-duplicates-wave7.yml`,
  materializing
  `mod_metab_qc_duplicates_ui()`
  in
  `R/mod_metab_qc_duplicates_ui.R`
  and
  `mod_metab_qc_duplicates_server()`
  in
  `R/mod_metab_qc_duplicates_server.R`
  while leaving
  `R/mod_metab_qc_duplicates.R`
  behind as the breadcrumb stub for the module identity.
- `DESCRIPTION` `Collate:` now loads
  `mod_metab_qc_duplicates_ui_helpers.R`,
  `mod_metab_qc_duplicates_render_helpers.R`,
  `mod_metab_qc_duplicates_server_helpers.R`,
  `mod_metab_qc_duplicates_ui.R`,
  `mod_metab_qc_duplicates_server.R`,
  and then
  `mod_metab_qc_duplicates.R`,
  with the applied load-order artifact recorded in
  `tools/refactor/collate-metab-qc-duplicates-wave7.txt`.
- The focused duplicate-module gate loader in
  `tests/testthat/test-metab-01p-qc-duplicates-module-characterization.R`
  now probes
  `R/mod_metab_qc_duplicates_ui.R`
  and
  `R/mod_metab_qc_duplicates_server.R`
  ahead of
  `R/mod_metab_qc_duplicates.R`
  so the wrapper characterization surface survives the final live apply
  boundary.
- `tools/refactor/verify_refactor.R` passed again for
  `tools/refactor/manifest-metab-qc-duplicates-wave7.yml` before the live
  apply, and `tools/refactor/check_wave_apply.R` passed after the live apply.
- The focused duplicate-module gate reran green again after the wave-7 apply
  via direct `testthat::test_file()` because this worktree still does not
  include `renv/activate.R` for `tools/test_with_renv.R`.
- Live classification now records
  `R/mod_metab_qc_duplicates.R`
  at `49` lines with `0` top-level functions after the wave-7 wrapper
  entrypoint apply.
- Treat the manual metabolomics QC duplicate-module target as complete:
  `R/mod_metab_qc_duplicates.R`
  is now a breadcrumb stub and bucket 0 can move on to the next manual
  target.
