# Metabolomics QC Intensity Module Seam Map

## Goal

Introduce bounded top-level seams in R/mod_metab_qc_intensity.R while keeping the live metabolomics intensity filter behavior frozen behind the existing public entry points.

## Current Position In The Flow

- target: `R/mod_metab_qc_intensity.R`
- classification: `review`
- checkpoint reached: `eighth bounded revert-reporting seam is now live in R/mod_metab_qc_intensity.R`
- next step: `Treat R/mod_metab_qc_intensity.R as the completed metabolomics QC intensity wrapper identity for this bucket; any later extraction from the remaining revert/apply orchestration shell is optional cleanup rather than required stabilization work.`

## Existing Safety Net

- `tests/testthat/test-metab-01q-qc-intensity-module-characterization.R`

## Notes

Manual bucket 0 metabolomics QC intensity-module stabilization target.

- This file now records the current metabolomics QC intensity-module seam stop
  point after two bounded render-path extractions plus one summary-text
  extraction plus one per-assay stats/count extraction plus one QC plot/update
  extraction plus one workflow state-save extraction plus one grouped
  success-reporting extraction plus one revert reset/reporting extraction.
- Classification refreshed on April 18, 2026 keeps
  `R/mod_metab_qc_intensity.R`
  at `483` lines with `10` top-level functions, a `114` line largest top-level
  function, and label `review`.
- The focused gate now lives in
  `tests/testthat/test-metab-01q-qc-intensity-module-characterization.R`
  and freezes the direct seam contracts for
  `buildMetabIntensityAssayTabsUi()`
  and
  `renderMetabIntensityFilterPlot()`
  and
  `buildMetabIntensityFilterStats()`
  and
  `updateMetabIntensityFilterQcPlot()`
  and
  `saveMetabIntensityFilterState()`
  and
  `buildMetabIntensityFilterSummary()`
  and
  `reportMetabIntensityFilterSuccess()`
  and
  `reportMetabIntensityFilterRevertSuccess()`
  plus the wrapper handoff from
  `mod_metab_qc_intensity_server()`
  into those seams across the `assay_results_tabs` and `filter_plot` render
  paths plus the `apply_filter` observer's stats/count, QC plot/update,
  workflow state-save, and success-reporting handoffs plus the
  `revert_filter` observer's reset/reporting handoff.
- The gate loader probes
  `R/mod_metab_qc_intensity_ui_helpers.R`
  ahead of
  `R/mod_metab_qc_intensity.R`,
  so the characterization surface survives a later exact-source helper apply.
- The first bounded live seam now sits in
  `R/mod_metab_qc_intensity.R`
  at line `105`
  as `buildMetabIntensityAssayTabsUi()`, which now owns the null passthrough,
  per-assay tab assembly, stat-row labeling, and tabset namespacing before
  control returns to the intensity-filter module wrapper.
- The `assay_results_tabs` render path in
  `R/mod_metab_qc_intensity.R`
  at line `471`
  now delegates through that helper instead of building the full per-assay
  results UI inline inside `mod_metab_qc_intensity_server()`.
- The second bounded live seam now sits in
  `R/mod_metab_qc_intensity.R`
  at line `153`
  as `renderMetabIntensityFilterPlot()`, which now owns the plot `req()`
  guard plus the `grob`/`gtable` and `ggplot` dispatch before control returns
  to the intensity-filter module wrapper.
- The `filter_plot` render path in
  `R/mod_metab_qc_intensity.R`
  at line `479`
  now delegates through that helper instead of branching on plot classes
  inline inside `mod_metab_qc_intensity_server()`.
- The third bounded live seam now sits in
  `R/mod_metab_qc_intensity.R`
  at line `250`
  as `buildMetabIntensityFilterSummary()`, which now owns the total-metabolite
  aggregation plus the result-text assembly before control returns to the
  intensity-filter module wrapper.
- The new success-reporting seam in
  `R/mod_metab_qc_intensity.R`
  now delegates through that helper instead of rebuilding the filter result
  text inline in the post-save completion tail.
- The fourth bounded live seam now sits in
  `R/mod_metab_qc_intensity.R`
  at line `172`
  as `buildMetabIntensityFilterStats()`, which now owns the metabolite-id
  unique-count fallback handling, per-assay stats data-frame assembly, and
  percent-retained calculation before control returns to the intensity-filter
  module wrapper.
- The `apply_filter` observer in
  `R/mod_metab_qc_intensity.R`
  at line `411`
  now delegates through that helper instead of counting original and filtered
  per-assay metabolites inline inside `mod_metab_qc_intensity_server()`.
- The fifth bounded live seam now sits in
  `R/mod_metab_qc_intensity.R`
  at line `201`
  as `updateMetabIntensityFilterQcPlot()`, which now owns the
  `updateMetaboliteFiltering()` handoff, warning fallback, and reactive
  filter-plot update before control returns to the intensity-filter module
  wrapper.
- The `apply_filter` observer in
  `R/mod_metab_qc_intensity.R`
  at line `417`
  now delegates through that helper instead of refreshing the QC tracking plot
  inline inside `mod_metab_qc_intensity_server()`.
- The sixth bounded live seam now sits in
  `R/mod_metab_qc_intensity.R`
  at line `227`
  as `saveMetabIntensityFilterState()`, which now owns the
  `state_manager$saveState()` handoff, state-name reuse, config wiring, and
  description formatting before control returns to the intensity-filter module
  wrapper.
- The `apply_filter` observer in
  `R/mod_metab_qc_intensity.R`
  at line `423`
  now delegates through that helper instead of persisting the filtered
  workflow state inline inside `mod_metab_qc_intensity_server()`.
- The seventh bounded live seam now sits in
  `R/mod_metab_qc_intensity.R`
  at line `297`
  as `reportMetabIntensityFilterSuccess()`, which now owns the summary-helper
  handoff, `filter_results` render binding, completion log, working-notification
  dismissal, and success toast before control returns to the intensity-filter
  module wrapper.
- The `apply_filter` observer in
  `R/mod_metab_qc_intensity.R`
  at line `431`
  now delegates through that helper instead of wiring the result render,
  completion log, and success notifications inline inside
  `mod_metab_qc_intensity_server()`.
- The eighth bounded live seam now sits in
  `R/mod_metab_qc_intensity.R`
  at line `339`
  as `reportMetabIntensityFilterRevertSuccess()`, which now owns the revert
  result render binding, reactive reset, completion log, and success toast
  before control returns to the metabolomics intensity-filter module wrapper.
- The `revert_filter` observer in
  `R/mod_metab_qc_intensity.R`
  at line `447`
  now delegates through that helper instead of wiring the result render,
  reactive reset, completion log, and success notification inline inside
  `mod_metab_qc_intensity_server()`.
- Because this worktree does not include `renv/activate.R`,
  `tools/test_with_renv.R`
  is not available for this target and the focused gate reran green via direct
  `testthat::test_file()`.
- The wrapper now sits inside the playbook's ideal size band as a public
  orchestrator shell, so this manual target no longer blocks bucket 0
  stabilization work.
