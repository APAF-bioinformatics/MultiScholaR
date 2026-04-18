# Lipid Design Module Seam Map

## Goal

Document the wrapper-level stabilization stop point for
[R/mod_lipid_design.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_design.R:1)
while keeping the public lipid-design module contract stable.

## Current Position In The Flow

- target:
  `/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_design.R`
- classification:
  `direct-extraction-ready`
- active stop point:
  The live helper-wave apply and the reviewed entrypoint-wave apply are now in
  place. The preview, builder, and import helper set from
  [tools/refactor/manifest-lipid-design-module-wave1.yml](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/manifest-lipid-design-module-wave1.yml:1)
  now lives in
  [R/mod_lipid_design_preview_helpers.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_design_preview_helpers.R:1),
  [R/mod_lipid_design_builder_helpers.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_design_builder_helpers.R:1),
  and
  [R/mod_lipid_design_import_helpers.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_design_import_helpers.R:1),
  while the remaining public UI and server entrypoints from
  [tools/refactor/manifest-lipid-design-module-wave2.yml](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/manifest-lipid-design-module-wave2.yml:1)
  now live in
  [R/mod_lipid_design_ui.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_design_ui.R:1)
  and
  [R/mod_lipid_design_server.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_design_server.R:1),
  leaving
  [R/mod_lipid_design.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_design.R:1)
  as a breadcrumb wrapper.
- next step:
  No further stabilization work is required for this target unless later
  documentation regeneration or dependency-complete package verification is
  requested. `scripts/classify_target.py --json R/mod_lipid_design.R` now
  reports `direct-extraction-ready` at `52` lines and `0` top-level
  functions.

## Existing Safety Net

- focused wrapper gate:
  - `tests/testthat/test-lipid-14-design-module-seams.R`
- replay command:
  - `Rscript -e "source('R/mod_lipid_design_preview_helpers.R'); source('R/mod_lipid_design_builder_helpers.R'); source('R/mod_lipid_design_import_helpers.R'); source('R/mod_lipid_design_ui.R'); source('R/mod_lipid_design_server.R'); source('R/mod_lipid_design.R'); testthat::test_file('tests/testthat/test-lipid-14-design-module-seams.R', stop_on_failure = TRUE)"`

## Notes

- April 15, 2026 stabilize-mode checkpoint applied the second bounded
  entrypoint wave live via
  [tools/refactor/manifest-lipid-design-module-wave2.yml](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/manifest-lipid-design-module-wave2.yml:1)
  into
  [R/mod_lipid_design_ui.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_design_ui.R:1)
  and
  [R/mod_lipid_design_server.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_design_server.R:1).
- The live collate artifact now exists at
  [tools/refactor/collate-lipid-design-module-wave2.txt](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/collate-lipid-design-module-wave2.txt:1),
  `DESCRIPTION` now collates the helper plus UI/server entrypoint files
  immediately before `R/mod_lipid_design.R`, and the focused source-based gate
  now loads the extracted UI/server entrypoints before the wrapper
  breadcrumb.
- Post-apply `check_wave_apply.R --manifest
  tools/refactor/manifest-lipid-design-module-wave2.yml` passed, and
  `scripts/classify_target.py --json R/mod_lipid_design.R` now reports
  `direct-extraction-ready` at `52` lines and `0` top-level functions while
  `R/mod_lipid_design_ui.R` and `R/mod_lipid_design_server.R` land at `75`
  and `79` lines.
- The public lipid design module surface now lives in dedicated helper/UI/server
  files while `R/mod_lipid_design.R` remains a breadcrumb stub, so this target
  no longer needs additional stabilization work unless later doc regeneration
  is requested in a dependency-complete environment.
- April 15, 2026 completed one bounded staged-wave checkpoint by verifying and
  staging
  [tools/refactor/manifest-lipid-design-module-wave1.yml](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/manifest-lipid-design-module-wave1.yml:1),
  extracting the seam-ready preview, builder, and import helper set from live
  [R/mod_lipid_design.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_design.R:1)
  into staged files
  [R/mod_lipid_design_preview_helpers.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/staging/wave1_lipidomics_design_module_host_helpers/R/mod_lipid_design_preview_helpers.R:1),
  [R/mod_lipid_design_builder_helpers.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/staging/wave1_lipidomics_design_module_host_helpers/R/mod_lipid_design_builder_helpers.R:1),
  and
  [R/mod_lipid_design_import_helpers.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/staging/wave1_lipidomics_design_module_host_helpers/R/mod_lipid_design_import_helpers.R:1).
- The staged collate artifact now exists at
  [collate-lipid-design-module-wave1.txt](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/staging/wave1_lipidomics_design_module_host_helpers/collate-lipid-design-module-wave1.txt:1),
  ordering the staged helper files before `mod_lipid_design.R` for a later
  live apply.
- Live post-staging classification for
  [R/mod_lipid_design.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_design.R:1)
  remains `860` lines, `12` top-level functions, and max top-level function
  length `70`; the staged helper files land at `31`, `200`, and `435` lines.
- The focused lipid-design wrapper gate reran green after the staged-wave
  checkpoint, and the target remains `in_progress` because the next clean stop
  point is staged-wave review for live apply readiness rather than another
  live seam.
- April 15, 2026 completed one bounded live-apply checkpoint by applying
  [tools/refactor/manifest-lipid-design-module-wave1.yml](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/manifest-lipid-design-module-wave1.yml:1)
  into
  [R/mod_lipid_design_preview_helpers.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_design_preview_helpers.R:1),
  [R/mod_lipid_design_builder_helpers.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_design_builder_helpers.R:1),
  and
  [R/mod_lipid_design_import_helpers.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_design_import_helpers.R:1).
- The live collate artifact now exists at
  [tools/refactor/collate-lipid-design-module-wave1.txt](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/collate-lipid-design-module-wave1.txt:1),
  `DESCRIPTION` now collates the three helper files immediately before
  `R/mod_lipid_design.R`, and the focused source-based gate now loads the
  helper files before the wrapper.
- Post-apply `check_wave_apply.R --manifest
  tools/refactor/manifest-lipid-design-module-wave1.yml` passed, and
  `scripts/classify_target.py --json R/mod_lipid_design.R` now reports
  `review` at `204` lines and `2` top-level functions.
- April 15, 2026 completed one bounded staged-wave checkpoint by verifying and
  staging
  [tools/refactor/manifest-lipid-design-module-wave2.yml](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/manifest-lipid-design-module-wave2.yml:1)
  for the remaining public entrypoint surface in
  [R/mod_lipid_design.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_design.R:1),
  materializing staged files
  [R/mod_lipid_design_ui.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/staging/wave2_lipidomics_design_module_entrypoints/R/mod_lipid_design_ui.R:1)
  and
  [R/mod_lipid_design_server.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/staging/wave2_lipidomics_design_module_entrypoints/R/mod_lipid_design_server.R:1).
- The staged collate artifact now exists at
  [collate-lipid-design-module-wave2.txt](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/staging/wave2_lipidomics_design_module_entrypoints/collate-lipid-design-module-wave2.txt:1),
  ordering the live helper files, staged UI/server entrypoints, and
  `R/mod_lipid_design.R` for later live apply review.
- The focused lipid-design wrapper gate reran green after the staged-wave
  checkpoint via direct `testthat::test_file(...)` with `58` passes because
  this worktree still lacks `renv/activate.R`.
- Live post-staging classification for
  [R/mod_lipid_design.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_design.R:1)
  remains `204` lines, `2` top-level functions, and label `review`; the staged
  UI and server entrypoints land at `75` and `79` lines, so the target remains
  `in_progress` with the next clean stop point at staged-wave review for live
  apply readiness rather than another live seam.
- Keep future work review-first around the remaining public wrapper/UI/server
  surface; do not create ad hoc live helper files outside exact-source waves.
- April 15, 2026 introduced one bounded live seam in
  [R/mod_lipid_design.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_design.R:117)
  by extracting the saved-preview render path into
  `formatLipidDesignAssaysPreview()` and
  `registerLipidDesignPreviewOutputs()`.
- The compact wrapper now delegates `design_matrix_preview`,
  `contrasts_preview`, and `assays_preview` through that seam at
  [R/mod_lipid_design.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_design.R:763),
  keeping the current `workflow_data` handoff stable through one top-level
  preview-registration seam.
- April 15, 2026 introduced one bounded live seam in
  [R/mod_lipid_design.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_design.R:148)
  by extracting the builder-module registration fallback into
  `registerLipidDesignBuilderModule()`.
- The compact wrapper now delegates the `mod_lipid_design_builder_server()`
  registration and missing-builder `reactiveVal(NULL)` fallback through that
  seam at
  [R/mod_lipid_design.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_design.R:748),
  keeping the current builder handoff stable through one top-level
  registration seam.
- April 15, 2026 introduced one bounded live seam in
  [R/mod_lipid_design.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_design.R:174)
  by extracting the builder-results save observer into
  `runLipidDesignBuilderObserverShell()` and
  `registerLipidDesignBuilderResultsObserver()`.
- The compact wrapper now delegates the `builder_results_rv()` observer
  registration and its save-to-workflow/save-to-disk shell through that seam
  at
  [R/mod_lipid_design.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_design.R:813),
  keeping the current builder result handoff stable through one top-level
  observer seam.
- April 15, 2026 introduced one bounded live seam in
  [R/mod_lipid_design.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_design.R:348)
  by extracting the import bootstrap and modal shell into
  `initializeLipidDesignImportBootstrap()`,
  `buildLipidDesignImportModal()`, and
  `registerLipidDesignImportModalShell()`.
- The compact wrapper now delegates the `resolved_volumes` bootstrap,
  `shinyDirChoose()` registration, import modal UI construction, and selected
  directory preview handoff through that seam at
  [R/mod_lipid_design.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_design.R:450),
  keeping the current import observer body untouched for the next checkpoint.
- The focused lipid-design wrapper gate now exists in
  [tests/testthat/test-lipid-14-design-module-seams.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tests/testthat/test-lipid-14-design-module-seams.R:1)
  and now also characterizes the builder registration handoff, fallback,
  builder-results observer handoff, import bootstrap/modal shell, and import
  confirmation observer handoff;
  it reran green via direct `testthat::test_file(...)` because this worktree
  still lacks `renv/activate.R`.
- Post-seam classification for
  [R/mod_lipid_design.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_design.R:1)
  is `860` lines, `12` top-level functions, and max top-level function length
  `70`; `scripts/classify_target.py` still reports `review`, so the backlog
  target remains `in_progress` even though the wrapper has no remaining live
  observer bodies.
- April 15, 2026 introduced one bounded live seam in
  [R/mod_lipid_design.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_design.R:434)
  by extracting the import-confirmation flow into
  `runLipidDesignImportConfirmationShell()` and
  `registerLipidDesignImportConfirmationObserver()`.
- The compact wrapper now delegates the `input$confirm_import` observer
  registration, import-path resolution, and imported-design hydration through
  that seam at
  [R/mod_lipid_design.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_design.R:818),
  leaving the file at a clean top-level seam boundary for staged extraction
  planning.
