# Metabolomics QC ITSD Module Seam Map

## Goal

Introduce bounded top-level seams in `R/mod_metab_qc_itsd.R` while keeping the
live metabolomics internal-standard QC module behavior frozen behind the
existing public entry points.

## Current Position In The Flow

- target: `R/mod_metab_qc_itsd.R`
- classification: `direct-extraction-ready` (structural classifier)
- checkpoint reached: `wrapper-entrypoint wave 1 is applied live for R/mod_metab_qc_itsd.R`
- next step: `Optional cleanup only: if later decomposition is wanted, draft and stage an exact-source helper wave from R/mod_metab_qc_itsd.R while keeping the applied UI/server entrypoints stable.`

## Existing Safety Net

- `tests/testthat/test-metab-01r-qc-itsd-module-characterization.R`
- `tests/testthat/test-metab-01o-qc-module-characterization.R`

## Notes

Manual bucket 0 metabolomics QC ITSD-module stabilization target.

- Classification refreshed on April 18, 2026 records
  `R/mod_metab_qc_itsd.R`
  at `444` lines with `6` top-level functions, a `58` line largest
  top-level function, and label `direct-extraction-ready`.
- The focused gate now lives in
  `tests/testthat/test-metab-01r-qc-itsd-module-characterization.R`
  and now loads
  `R/mod_metab_qc_itsd.R`,
  `R/mod_metab_qc_itsd_ui.R`,
  and
  `R/mod_metab_qc_itsd_server.R`
  so it continues to freeze the direct seam contract for
  `mod_metab_qc_itsd_ui()`,
  `analyzeMetabQcItsdData()`,
  `buildMetabQcItsdSummaryUi()`,
  `buildMetabQcItsdVizTabsUi()`,
  `buildMetabQcItsdCvPlot()`,
  `buildMetabQcItsdIntensityPlot()`,
  and
  `runMetabQcItsdServerBody()`
  plus the public-wrapper handoff from
  `mod_metab_qc_itsd_server()`
  into the server-body seam and the downstream helper delegation for the
  analyze observer and summary, visualization-tab, CV-plot, and
  intensity-plot bindings after the live wrapper-entrypoint apply.
- The first bounded live seam now sits in
  `R/mod_metab_qc_itsd.R`
  at line `33`
  as `analyzeMetabQcItsdData()`, which now owns the
  `MetaboliteAssayData` guard, internal-standard pattern resolution,
  per-assay ID-column search order, metrics and long-data assembly, and
  result-text construction before control returns to the ITSD module wrapper.
- The second bounded live seam now sits in
  `R/mod_metab_qc_itsd.R`
  at line `176`
  as `buildMetabQcItsdSummaryUi()`, which now owns the null-state placeholder,
  per-assay median-CV summary-row assembly, status-icon and color selection,
  and explanatory threshold footer before control returns to the ITSD module
  wrapper.
- The third bounded live seam now sits in
  `R/mod_metab_qc_itsd.R`
  at line `226`
  as `buildMetabQcItsdVizTabsUi()`, which now owns the null-state exit,
  visualization tab-list assembly, namespaced plot-output registration, and
  tabset handoff before control returns to the ITSD module wrapper.
- The fourth bounded live seam now sits in
  `R/mod_metab_qc_itsd.R`
  at line `261`
  as `buildMetabQcItsdCvPlot()`, which now owns the CV-based ID reorder,
  threshold-status bucketing, and lollipop ggplot assembly before control
  returns to the ITSD module wrapper.
- The fifth bounded live seam now sits in
  `R/mod_metab_qc_itsd.R`
  at line `322`
  as `buildMetabQcItsdIntensityPlot()`, which now owns the long-data ggplot
  assembly, assay facet wiring, axis/legend labels, and rotated sample-axis
  theme settings before control returns to the ITSD module wrapper.
- The sixth bounded live seam now sits in
  `R/mod_metab_qc_itsd.R`
  at line `350`
  as `runMetabQcItsdServerBody()`, which now owns the module callback
  reactive-state setup, `input$analyze_is` observer registration, and the
  summary, visualization-tab, CV-plot, and intensity-plot render bindings
  before control returns to the public ITSD server wrapper.
- The reviewed exact-source wrapper-entrypoint wave from
  `tools/refactor/manifest-metab-qc-itsd-wave1.yml`
  is now applied live, where
  `R/mod_metab_qc_itsd_ui.R`
  at line `5`
  and
  `R/mod_metab_qc_itsd_server.R`
  at line `6`
  preserve the current public entrypoints without hand-rewriting their bodies.
- `DESCRIPTION`
  `Collate:`
  now orders
  `R/mod_metab_qc_itsd.R`,
  `R/mod_metab_qc_itsd_ui.R`,
  and
  `R/mod_metab_qc_itsd_server.R`
  so the helper-bearing source stays loaded before the dedicated server
  entrypoint file.
- Post-apply verification passed via
  `Rscript tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-qc-itsd-wave1.yml`.
- Because this worktree does not include `renv/activate.R`,
  `tools/test_with_renv.R`
  is not available for this target and the focused gate reran green via direct
  `testthat::test_file()`.
- The metabolomics QC ITSD wrapper identity now sits inside the playbook ideal
  size band and no longer blocks bucket 0 stabilization work; any later
  helper-only extraction from `R/mod_metab_qc_itsd.R` is optional cleanup
  rather than required stabilization.
