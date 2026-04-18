# Lipid QC ITSD Module Seam Map

## Goal

Document the wrapper-level stabilization stop point for mod_lipid_qc_itsd.R while keeping the public internal-standard QC module contract stable.

## Current Position In The Flow

- target: `/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_itsd.R`
- classification: `review`
- active stop point: the third live wrapper seam is now in
  [R/mod_lipid_qc_itsd.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_itsd.R:154)
  as `registerLipidItsdCvPlotOutput()`, and
  [mod_lipid_qc_itsd_server()](</home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_itsd.R:209>)
  now delegates the `output$cv_plot` registration through that seam at
  [R/mod_lipid_qc_itsd.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_itsd.R:399).
- next step: `No further live seam is currently required for this stabilization lane; if future structural work resumes, start from the existing focused characterization gate and review whether the remaining intensity plot should move.`

## Existing Safety Net

- focused wrapper gate:
  - [tests/testthat/test-lipid-01e-qc-itsd-module-seams.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tests/testthat/test-lipid-01e-qc-itsd-module-seams.R:1)
- replay command:
  - `Rscript -e "testthat::test_file('tests/testthat/test-lipid-01e-qc-itsd-module-seams.R', stop_on_failure = TRUE)"`

## Notes

- No prior target handover existed for this wrapper lane.
- April 16, 2026 introduced one bounded live seam in
  [R/mod_lipid_qc_itsd.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_itsd.R:81)
  by extracting the `output$is_summary` render registration into
  `registerLipidItsdSummaryOutput()`.
- The live wrapper now delegates the empty-state prompt, per-assay median-CV
  summary list, and legend shell through that helper without changing the
  public `mod_lipid_qc_itsd_server()` entrypoint.
- Focused seam characterization in
  [tests/testthat/test-lipid-01e-qc-itsd-module-seams.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tests/testthat/test-lipid-01e-qc-itsd-module-seams.R:1)
  now freezes:
  - the empty-state summary prompt and muted help styling
  - the per-assay summary rendering for good, acceptable, and review CV bands
  - the live module-server handoff into
    `registerLipidItsdSummaryOutput()`
- The repo-local `tools/test_with_renv.R` wrapper is not usable in this
  checkout because `renv/activate.R` is absent, so the focused gate currently
  reruns via direct `testthat::test_file(...)`.
- April 16, 2026 introduced a second bounded live seam in
  [R/mod_lipid_qc_itsd.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_itsd.R:123)
  by extracting the `output$is_viz_tabs` render registration into
  `registerLipidItsdVisualizationTabsOutput()`.
- The live wrapper now delegates the visualization-tab shell, tabset id, and
  plot output bindings through that helper without changing the public
  `mod_lipid_qc_itsd_server()` entrypoint.
- Focused seam characterization in
  [tests/testthat/test-lipid-01e-qc-itsd-module-seams.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tests/testthat/test-lipid-01e-qc-itsd-module-seams.R:1)
  now also freezes:
  - the empty-state visualization-tab null return
  - the stable visualization-tab labels, namespaced tabset id, and plot output
    ids
  - the live module-server handoff into
    `registerLipidItsdVisualizationTabsOutput()`
- The focused wrapper gate reran green after the live seam with `14`
  source-based passes.
- Post-seam classification for
  [R/mod_lipid_qc_itsd.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_itsd.R:1)
  is `418` lines, `4` top-level functions, max top-level function length
  `260`, `1` observer, `3` renderers, and labels
  `high-risk-wrapper` plus `needs-seam-introduction`.
- The focused wrapper gate reran green after the second live seam with `25`
  source-based passes.
- This target remains `in_progress`: the summary and visualization-tab
  registrations are now top-level, but the wrapper still carries both plot
  renders and the analysis observer, so the next clean stop point is the CV
  plot registration rather than staged extraction.
- April 16, 2026 introduced a third bounded live seam in
  [R/mod_lipid_qc_itsd.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_itsd.R:154)
  by extracting the `output$cv_plot` render registration into
  `registerLipidItsdCvPlotOutput()`.
- The live wrapper now delegates the CV lollipop-plot registration, CV-band
  status mapping, factor reordering, and plot shell through that helper
  without changing the public `mod_lipid_qc_itsd_server()` entrypoint.
- Focused seam characterization in
  [tests/testthat/test-lipid-01e-qc-itsd-module-seams.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tests/testthat/test-lipid-01e-qc-itsd-module-seams.R:1)
  now also freezes:
  - the reordered `is_id` factor levels used by the CV plot
  - the stable CV status bands, labels, facet shell, and flipped coordinates
  - the live module-server handoff into
    `registerLipidItsdCvPlotOutput()`
- The focused wrapper gate reran green after the third live seam with `38`
  source-based passes via direct `testthat::test_file(...)`.
- Post-third-seam classification for
  [R/mod_lipid_qc_itsd.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_itsd.R:1)
  is `426` lines, `5` top-level functions, max top-level function length
  `218`, `1` observer, `2` renderers, and label `review`.
- The classifier still emits the generic next step
  `Add focused characterization before structural edits.` even though the
  focused wrapper gate now covers the three live seams already introduced.
- Lipid-QC ITSD wrapper stabilization is now complete for this backlog
  target; no further live seam is currently required in
  [R/mod_lipid_qc_itsd.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_itsd.R:1).
