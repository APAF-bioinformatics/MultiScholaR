# Metabolomics Design Module Seam Map

## Goal

Introduce low-risk top-level seams in `R/mod_metab_design.R` while keeping the
public `mod_metab_design_ui()` and `mod_metab_design_server()` entry points
behaviorally stable.

## Current Position In The Flow

- target: `R/mod_metab_design.R`
- classification: `applied import-wrapper helper wave; manual wrapper target complete`
- checkpoint reached: `live import-wrapper helper wave for mod_metab_design_server()`
- next step: `No further stabilization work remains in R/mod_metab_design.R; if bucket 0 design work continues, shift to the oversized builder twin target in R/mod_metab_design_builder.R.`

## Existing Safety Net

- `tests/testthat/test-metab-04-design-module-characterization.R`

## Notes

Manual bucket 0 metabolomics design module stabilization target.

- This target previously had no active handover; this file now records the
  first metabolomics design stabilization stop point.
- Classification refresh on April 16, 2026 kept
  `R/mod_metab_design.R`
  at `717` lines with `2` top-level functions, a `593` line largest top-level
  function, and labels `high-risk-wrapper` / `needs-seam-introduction`
  before this checkpoint.
- The focused characterization gate now lives in
  `tests/testthat/test-metab-04-design-module-characterization.R`
  and freezes the current preview, state-output, and builder-module
  registration contracts for
  `registerMetabDesignPreviewOutputs()` and
  `registerMetabDesignStateOutputs()`, plus
  `registerMetabDesignBuilderModule()`.
- The first bounded seam now also lives in-place in
  `R/mod_metab_design.R`
  as `registerMetabDesignPreviewOutputs()`.
- `output$design_matrix_preview`,
  `output$contrasts_preview`, and
  `output$assays_preview`
  now route through `registerMetabDesignPreviewOutputs()`, so the saved-preview
  render registration stays in one top-level stop point instead of remaining
  inline in the wrapper body.
- The second bounded seam now also lives in-place in
  `R/mod_metab_design.R`
  as `registerMetabDesignStateOutputs()`.
- `output$data_available` and
  `output$design_matrix_exists`
  now route through `registerMetabDesignStateOutputs()`, so the reactive
  availability registration and `outputOptions()` wiring stay in one top-level
  stop point instead of remaining inline in the wrapper body.
- The third bounded seam now also lives in-place in
  `R/mod_metab_design.R`
  as `registerMetabDesignBuilderModule()`.
- `builder_results_rv`
  now routes through `registerMetabDesignBuilderModule()`, so the optional
  builder-server registration path and missing-builder fallback stay in one
  top-level stop point instead of remaining inline in the wrapper body.
- The fourth bounded seam now also lives in-place in
  `R/mod_metab_design.R`
  as `registerMetabDesignBuilderResultsObserver()`.
- the wrapper now routes the builder-results observer/save path through
  `registerMetabDesignBuilderResultsObserver()` at:
  `R/mod_metab_design.R:760`
- the focused characterization gate now also freezes the observer-registration
  seam and delegated builder-results handoff for
  `registerMetabDesignBuilderResultsObserver()`.
- The focused characterization gate now also freezes the remaining import
  wrapper shell in `mod_metab_design_server()` by asserting the project-base
  volume injection, `shinyDirChoose()` registration, namespaced import modal
  controls, `import_dir_path` resolution, and top-level observer registration
  contract before the next structural seam.
- The fifth bounded seam now also lives in-place in
  `R/mod_metab_design.R`
  as `initializeMetabDesignImportBootstrap()`.
- the wrapper now routes `resolved_volumes` setup and `shinyDirChoose()`
  registration through `initializeMetabDesignImportBootstrap()` at:
  `R/mod_metab_design.R:367`
- the focused characterization gate now also freezes the import bootstrap seam
  directly, covering project-base volume injection, chooser registration, and
  the wrapper handoff from `mod_metab_design_server()` into
  `initializeMetabDesignImportBootstrap()`.
- The focused gate reran green after this checkpoint via a direct
  `testthat::test_file()` invocation because this worktree does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- Post-checkpoint classification now moves
  `R/mod_metab_design.R`
  to `796` lines with `8` top-level functions, a `394` line largest top-level
  function, and label `review`.
- The sixth bounded seam now also lives in-place in
  `R/mod_metab_design.R`
  as `registerMetabDesignImportModalShell()`.
- the wrapper now routes the `show_import_modal` observer and
  `output$import_dir_path` render registration through
  `registerMetabDesignImportModalShell()` at:
  `R/mod_metab_design.R:403`
- the focused characterization gate now also freezes the import modal-shell
  seam directly, covering modal registration, namespaced control ids, import
  path rendering, and the wrapper handoff from `mod_metab_design_server()` into
  `registerMetabDesignImportModalShell()`.
- Post-checkpoint classification now moves
  `R/mod_metab_design.R`
  to `823` lines with `9` top-level functions, a `379` line largest top-level
  function, and label `review`.
- The seventh bounded seam now also lives in-place in
  `R/mod_metab_design.R`
  as `resolveMetabDesignImportPreflight()`.
- the wrapper now routes `confirm_import` path resolution plus required
  `design_matrix.tab` / `assay_manifest.txt` validation through
  `resolveMetabDesignImportPreflight()` at:
  `R/mod_metab_design.R:445`
- the focused characterization gate now also freezes the import preflight seam
  directly, covering the resolved import-path contract, returned artifact paths,
  missing-file failure messages, and the wrapper handoff from
  `mod_metab_design_server()` into `resolveMetabDesignImportPreflight()`.
- The focused gate reran green after this checkpoint via a direct
  `testthat::test_file()` invocation because this worktree does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- Post-checkpoint classification now keeps
  `R/mod_metab_design.R`
  at `review` with `858` lines, `10` top-level functions, and a `371` line
  largest top-level function.
- The eighth bounded seam now also lives in-place in
  `R/mod_metab_design.R`
  as `hydrateMetabDesignImportArtifacts()`.
- the wrapper now routes `confirm_import` config resolution plus
  `design_matrix.tab` / `assay_manifest.txt` / `data_cln_*.tab` hydration
  through `hydrateMetabDesignImportArtifacts()` at:
  `R/mod_metab_design.R:488`
- the focused characterization gate now also freezes the import hydration seam
  directly, covering config loading, design/manifest/assay artifact hydration,
  missing-assay failure messaging, and the wrapper handoff from
  `mod_metab_design_server()` into `hydrateMetabDesignImportArtifacts()`.
- The focused gate reran green after this checkpoint via a direct
  `testthat::test_file()` invocation because this worktree does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- The ninth bounded seam now also lives in-place in
  `R/mod_metab_design.R`
  as `hydrateMetabDesignImportMetadata()`.
- the wrapper now routes `confirm_import` metadata hydration for
  `column_mapping.json` loading or inference plus `contrast_strings.tab`
  parsing through `hydrateMetabDesignImportMetadata()` at:
  `R/mod_metab_design.R:579`
- the focused characterization gate now also freezes the import metadata seam
  directly for both inferred and JSON-backed column-mapping paths, and freezes
  the wrapper handoff from `mod_metab_design_server()` into
  `hydrateMetabDesignImportMetadata()`.
- The focused gate reran green after this checkpoint via a direct
  `testthat::test_file()` invocation because this worktree does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- Post-checkpoint classification now keeps
  `R/mod_metab_design.R`
  at `review` with `911` lines, `13` top-level functions, and a `220` line
  largest top-level function.
- The tenth bounded seam now also lives in-place in
  `R/mod_metab_design.R`
  as `hydrateMetabDesignImportWorkflowState()`.
- the wrapper now routes `confirm_import` post-import workflow-state hydration
  for design-matrix / contrasts / assay assignment plus `tech_rep_group`
  mutation through `hydrateMetabDesignImportWorkflowState()` at:
  `R/mod_metab_design.R:692`
- the focused characterization gate now also freezes the import workflow-state
  seam directly, including `contrasts_tbl` global-assignment handoff and
  `tech_rep_group` mutation, and freezes the wrapper handoff from
  `mod_metab_design_server()` into `hydrateMetabDesignImportWorkflowState()`.
- The focused gate reran green after this checkpoint via a direct
  `testthat::test_file()` invocation because this worktree does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- Post-checkpoint classification now keeps
  `R/mod_metab_design.R`
  at `review` with `934` lines, `14` top-level functions, and a `212` line
  largest top-level function.
- The eleventh bounded seam now also lives in-place in
  `R/mod_metab_design.R`
  as `createMetabDesignImportedS4Object()`.
- the wrapper now routes `confirm_import` imported S4 argument assembly,
  debug logging, and `createMetaboliteAssayData()` invocation through
  `createMetabDesignImportedS4Object()` at:
  `R/mod_metab_design.R:723`
- the focused characterization gate now also freezes the imported S4 seam
  directly, including blank-annotation normalization to `NA_character_`, and
  freezes the wrapper handoff from `mod_metab_design_server()` into
  `createMetabDesignImportedS4Object()`.
- The focused gate reran green after this checkpoint via a direct
  `testthat::test_file()` invocation because this worktree does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- Post-checkpoint classification now keeps
  `R/mod_metab_design.R`
  at `review` with `953` lines, `15` top-level functions, and a `176` line
  largest top-level function.
- The twelfth bounded seam now also lives in-place in
  `R/mod_metab_design.R`
  as `saveMetabDesignImportedS4State()`.
- the wrapper now routes `confirm_import` WorkflowState initialization plus
  `saveState()` persistence through `saveMetabDesignImportedS4State()` at:
  `R/mod_metab_design.R:778`
- the focused characterization gate now also freezes the imported
  state-manager seam directly, covering create-on-demand WorkflowState
  initialization, `metab_raw_data_s4` save arguments, and the wrapper handoff
  from `mod_metab_design_server()` into `saveMetabDesignImportedS4State()`.
- The focused gate reran green after this checkpoint via a direct
  `testthat::test_file()` invocation because this worktree does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- Post-checkpoint classification now keeps
  `R/mod_metab_design.R`
  at `review` with `966` lines, `17` top-level functions, and a `163` line
  largest top-level function.
- The thirteenth bounded seam now also lives in-place in
  `R/mod_metab_design.R`
  as `initializeMetabDesignImportedQcBaseline()`.
- the wrapper now routes `confirm_import` raw-data
  `updateMetaboliteFiltering()` baseline setup through
  `initializeMetabDesignImportedQcBaseline()` at:
  `R/mod_metab_design.R:804`
- the focused characterization gate now also freezes the imported QC-baseline
  seam directly, covering the `1_Raw_Data` metabolomics update contract, the
  warning-on-failure behavior, and the wrapper handoff from
  `mod_metab_design_server()` into
  `initializeMetabDesignImportedQcBaseline()`.
- The focused gate reran green after this checkpoint via a direct
  `testthat::test_file()` invocation because this worktree does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- Post-checkpoint classification now keeps
  `R/mod_metab_design.R`
  at `review` with `978` lines, `18` top-level functions, and a `153` line
  largest top-level function.
- The fourteenth bounded seam now also lives in-place in
  `R/mod_metab_design.R`
  as `completeMetabDesignImportedPostCheckpoint()`.
- the wrapper now routes the remaining imported completion handoff for
  `qc_trigger(TRUE)`, `tab_status$design_matrix <- "complete"`, and
  success-notification cleanup through
  `completeMetabDesignImportedPostCheckpoint()` at:
  `R/mod_metab_design.R:826`
- the focused characterization gate now also freezes the imported completion
  seam directly, covering the success log messages, conditional `qcTrigger`
  handoff, returned completion payload, and the wrapper handoff from
  `mod_metab_design_server()` into
  `completeMetabDesignImportedPostCheckpoint()`.
- The focused gate reran green after this checkpoint via a direct
  `testthat::test_file()` invocation because this worktree does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- Post-checkpoint classification now keeps
  `R/mod_metab_design.R`
  at `review` with `1003` lines, `19` top-level functions, and a `146` line
  largest top-level function.
- The fifteenth bounded seam now also lives in-place in
  `R/mod_metab_design.R`
  as `registerMetabDesignImportObserverShell()`.
- the wrapper now routes the remaining `confirm_import` observer shell through
  `registerMetabDesignImportObserverShell()` at:
  `R/mod_metab_design.R:864`
- the focused characterization gate now also freezes the direct observer-shell
  seam contract for success-path orchestration plus shared import-error
  notification handoff, and freezes the wrapper handoff from
  `mod_metab_design_server()` into `registerMetabDesignImportObserverShell()`.
- The focused gate reran green after this checkpoint via a direct
  `testthat::test_file()` invocation because this worktree does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- Post-checkpoint classification now keeps
  `R/mod_metab_design.R`
  at `review` with `1022` lines, `20` top-level functions, and a `146` line
  largest top-level function.
- The sixteenth bounded checkpoint now stages one exact-source import-wrapper
  helper wave via
  `tools/refactor/manifest-metab-design-wave1.yml`.
- the staged review target now lives in
  `tools/refactor/staging/wave1_metabolomics_design_import_helpers/R/mod_metab_design_import_helpers.R`
  with the emitted collate preview in
  `tools/refactor/staging/wave1_metabolomics_design_import_helpers/collate-metab-design-wave1.txt`.
- `verify_refactor.R` passed for the staged import-helper wave, freezing the
  exact selector set before any live apply step.
- The focused gate reran green after the staged-wave checkpoint via a direct
  `testthat::test_file()` invocation because this worktree does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- The seventeenth bounded checkpoint now applies
  `tools/refactor/manifest-metab-design-wave1.yml` into the live helper target
  `R/mod_metab_design_import_helpers.R`.
- the reviewed import-wrapper helper cluster now lives in
  `R/mod_metab_design_import_helpers.R`, and the extracted symbols no longer
  live inline in `R/mod_metab_design.R`.
- `DESCRIPTION` `Collate:` now loads `mod_metab_design_import_helpers.R`
  between `mod_metab_design_builder.R` and `mod_metab_design.R`, preserving the
  live helper-before-wrapper load order for the applied wave.
- the focused characterization gate now loads both
  `R/mod_metab_design_import_helpers.R` and `R/mod_metab_design.R`, and reran
  green after apply via a direct `testthat::test_file()` invocation because
  this worktree does not include `renv/activate.R` for `tools/test_with_renv.R`.
- Post-apply classification now records
  `R/mod_metab_design.R`
  at `439` lines with `7` top-level functions, a `146` line largest top-level
  function, and label `review`; the manual wrapper target is now below the
  stabilization budget and no longer needs further wrapper-seam work.
- the next structural follow-up, if design stabilization continues in bucket 0,
  belongs to `R/mod_metab_design_builder.R` rather than this wrapper target.
