# Metabolomics Import Module Seam Map

## Goal

Introduce bounded top-level seams in `R/mod_metab_import.R` while keeping the
live metabolomics import module behavior frozen behind the existing public
entry points.

## Current Position In The Flow

- target: `R/mod_metab_import.R`
- classification: `review`
- checkpoint reached: `display/status apply wave`
- next step: `Stage the final exact-source metabolomics import wrapper reduction rooted at resolveMetabImportColumnName(), resolveMetabImportSampleColumns(), mod_metab_import_ui(), and mod_metab_import_server() so R/mod_metab_import.R can collapse into a breadcrumb wrapper.`

## Existing Safety Net

- `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
- `tests/testthat/test-metab-00-import-detection-characterization.R`

## Notes

Manual bucket 0 metabolomics import module stabilization target.

- This target previously had no active handover; this file now records the
  current metabolomics import module seam stop point.
- `R/mod_metab_import.R`
  now sits at `517` lines with `4` top-level functions and remains in
  `review` after the reviewed display/status apply wave.
- The first focused gate now lives in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  and freezes the current contracts for:
  `resolveMetabImportColumnName()`
  and
  `resolveMetabImportSampleColumns()`,
  plus
  `buildMetabImportWorkflowPayload()`
  and
  `applyMetabImportWorkflowPayload()`,
  plus
  `prepareMetabImportAssaySelectionState()`,
  plus
  `finalizeMetabImportProcessingFeedback()`.
- The bounded live seams now sit in
  `R/mod_metab_import.R`
  and route the custom status outputs plus
  `get_metabolite_id_col()`,
  `get_annotation_col()`,
  and
  `get_sample_columns()`
  through the first helper pair, while the `process_import` observer now pushes
  assay assembly, optional sample sanitization, and workflow payload creation
  through `buildMetabImportWorkflowPayload()`.
- The third bounded live seam now sits in
  `R/mod_metab_import.R`
  as `applyMetabImportWorkflowPayload()`, which now owns workflow-data
  population, workflow-type synchronization, processing-log setup import
  storage, tab-completion state updates, and completion-summary logging for the
  `process_import` observer.
- The fourth bounded live seam now sits in
  `R/mod_metab_import.R`
  as `finalizeMetabImportProcessingFeedback()`, which now owns the
  `process_import` working-notification teardown plus the success toast and
  error logging/error-toast finalization branch.
- The fifth bounded live seam now sits in
  `R/mod_metab_import.R`
  as `prepareMetabImportAssaySelectionState()`, which now owns the nested
  `import_data()` header-read, format detection, importer selection, and
  select-input default preparation path before the observer applies the Shiny
  updates.
- The sixth bounded live seam now sits in
  `R/mod_metab_import.R`
  as `applyMetabImportAssaySelectionState()`, which now owns the nested
  `import_data()` local reactive-value commit, metabolite/annotation select
  input updates, optional internal-standard pattern text update, and import
  summary logging path before the observer falls through to its remaining error
  branch.
- The seventh bounded live seam now sits in
  `R/mod_metab_import.R`
  as `finalizeMetabImportAssaySelectionError()`, which now owns the nested
  `import_data()` error-message normalization plus the import-error logging and
  Shiny notification path for the wrapper.
- The eighth bounded live seam now sits in
  `R/mod_metab_import.R`
  as `handleMetabImportAssayFileSelection()`, which now owns the nested
  `assay1_file` observeEvent parseFilePaths lookup, local-path/output update,
  and post-selection import trigger handoff for the wrapper.
- The same file-selection seam now also owns the nested `assay2_file`
  observeEvent parseFilePaths lookup plus the assay-2 local-path/output update
  path for the wrapper, so both shinyFiles observers delegate through one
  top-level helper instead of keeping separate inline parse/update branches.
- The ninth bounded live seam now sits in
  `R/mod_metab_import.R`
  as `runMetabImportAssaySelection()`, which now owns the local
  `import_data()` wrapper's `req()` gate plus the prepare/apply/error-finalize
  helper handoff before the observer tree continues into the output renderers.
- The tenth bounded live seam now sits in
  `R/mod_metab_import.R`
  as `buildMetabImportFormatDetectionStatus()`, which now owns the
  `format_detection_status` render path's `req()` gate, confidence threshold
  coloring, vendor-label mapping, and alert-body assembly before the wrapper's
  render block falls through to the remaining validation outputs.
- The eleventh bounded live seam now sits in
  `R/mod_metab_import.R`
  as `buildMetabImportMetaboliteIdStatus()`, which now owns the
  `metabolite_id_status` render path's `req()` gate, unique-id counting, and
  missing-column status assembly before the wrapper falls through to the
  remaining annotation and custom validation outputs.
- The twelfth bounded live seam now sits in
  `R/mod_metab_import.R`
  as `buildMetabImportAnnotationStatus()`, which now owns the
  `annotation_status` render path's `req()` gate, present-column success state,
  missing-column fallback, and optional-column fallback before the wrapper
  falls through to the sample-column display path.
- The thirteenth bounded live seam now sits in
  `R/mod_metab_import.R`
  as `buildMetabImportSampleColumnsDisplay()`, which now owns the
  `sample_columns_display` render path's `req()` gate, over-ten truncation, and
  collapsed full-list fallback before the wrapper falls through to the
  available-columns display path.
- The fourteenth bounded live seam now sits in
  `R/mod_metab_import.R`
  as `buildMetabImportAvailableColumnsDisplay()`, which now owns the
  `available_columns_display` render path's `req()` gate and collapsed
  header-list formatting before the wrapper falls through to the remaining
  custom validation outputs.
- The fifteenth bounded live seam now sits in
  `R/mod_metab_import.R`
  as `buildMetabImportCustomMetaboliteIdStatus()`, which now owns the
  `metabolite_id_status_custom` render path's `req()` gate, empty-input
  prompt, case-insensitive column resolution, and unique-id status assembly
  before the wrapper falls through to the remaining custom annotation output.
- The sixteenth bounded live seam now sits in
  `R/mod_metab_import.R`
  as `buildMetabImportCustomAnnotationStatus()`, which now owns the
  `annotation_status_custom` render path's `req()` gate, optional-column
  fallback, case-insensitive column resolution, and found/missing status
  assembly before the wrapper falls through to the validation summary output.
- The seventeenth bounded live seam now sits in
  `R/mod_metab_import.R`
  as `buildMetabImportValidationSummary()`, which now owns the
  `validation_summary` render path's `req()` gates, metabolite/sample-column
  accessor delegation, validation helper handoff, success-summary assembly,
  warning-list fallback, and failure-list rendering before the wrapper falls
  through to the remaining import-status output.
- The eighteenth bounded live seam now sits in
  `R/mod_metab_import.R`
  as `buildMetabImportStatus()`, which now owns the `import_status` render
  path's completion-state gate, processing-log summary formatting, uppercase
  format display, and success-alert assembly before the wrapper falls through
  to the last inline observer shell.
- The nineteenth bounded live seam now sits in
  `R/mod_metab_import.R`
  as `runMetabImportProcessing()`, which now owns the `process_import`
  observer shell's assay-data/metabolite `req()` gates, accessor delegation,
  working-notification setup, sanitize-name logging/message branch, and
  workflow build/apply/finalize helper handoff.
- The twentieth bounded live seam now sits in
  `R/mod_metab_import.R`
  as `setupMetabImportShinyFiles()`, which now owns the remaining shinyFiles
  setup shell's volume fallback, chooser registration, assay-1 import trigger
  observer wiring, assay-2 path observer wiring, and parse-error logging path
  before the wrapper falls through to the local `import_data()` shell and
  remaining effective column-accessor reactives.
- The twenty-first bounded live seam now sits in
  `R/mod_metab_import.R`
  as `setupMetabImportColumnAccessors()`, which now owns the dropdown/custom
  metabolite-column resolution reactive, the dropdown/custom annotation-column
  resolution reactive, and the sample-column reactive's `req()` and resolver
  delegation before the wrapper falls through to the local `import_data()`
  shell and remaining output/observer registration.
- The focused gate reran green after the new seam via direct
  `testthat::test_file()` invocation because this worktree does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- The focused gate now also freezes the state-application seam contract in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  for `applyMetabImportAssaySelectionState()`.
- The focused gate now also freezes the assay-file selection seam contract in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  for `handleMetabImportAssayFileSelection()`.
- The focused gate now also freezes the assay-2 routing contract in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  by asserting the same helper updates `assay2_file` and `assay2_path` without
  firing the assay-1 import callback path.
- The focused gate now also freezes the local `import_data()` orchestration seam
  contract in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  by asserting `runMetabImportAssaySelection()` preserves the `req()` check,
  the prepare/apply delegation path, and the error-finalization handoff.
- The focused gate now also freezes the assay-selection error-finalization seam
  contract in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  for `finalizeMetabImportAssaySelectionError()`.
- The focused gate now also freezes the format-detection status render seam
  contract in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  for `buildMetabImportFormatDetectionStatus()` by asserting the warning/danger
  threshold mapping, vendor-label formatting, and `req()` gating contract.
- The focused gate now also freezes the metabolite-id validation status render
  seam contract in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  for `buildMetabImportMetaboliteIdStatus()` by asserting the success-path
  unique-id count, missing-column fallback, and `req()` gating contract.
- The focused gate now also freezes the annotation validation status render
  seam contract in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  for `buildMetabImportAnnotationStatus()` by asserting the found, missing,
  optional, and `req()` gating contracts.
- The focused gate now also freezes the sample-columns display render seam
  contract in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  for `buildMetabImportSampleColumnsDisplay()` by asserting the truncation,
  full-list fallback, and `req()` gating contracts.
- The focused gate now also freezes the available-columns display render seam
  contract in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  for `buildMetabImportAvailableColumnsDisplay()` by asserting the collapsed
  header-list formatting and `req()` gating contracts.
- The focused gate now also freezes the custom metabolite-id status render seam
  contract in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  for `buildMetabImportCustomMetaboliteIdStatus()` by asserting the empty-input
  prompt, case-insensitive resolution, missing-column fallback, and `req()`
  gating contracts.
- The focused gate now also freezes the custom annotation status render seam
  contract in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  for `buildMetabImportCustomAnnotationStatus()` by asserting the optional
  fallback, case-insensitive resolution, missing-column fallback, and `req()`
  gating contracts.
- The focused gate now also freezes the validation summary render seam
  contract in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  for `buildMetabImportValidationSummary()` by asserting the metabolite/sample
  accessor delegation, success/warning/failure report assembly, and `req()`
  gating contracts.
- The focused gate now also freezes the import-status render seam contract in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  for `buildMetabImportStatus()` by asserting the completion-summary alert
  assembly, uppercase format delegation, and incomplete-state `NULL`
  fallback contract.
- The focused gate now also freezes the `process_import` observer seam
  contract in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  for `runMetabImportProcessing()` by asserting accessor delegation,
  working-notification setup, sanitize-name messaging, payload/apply/finalize
  handoff, error-finalization fallback, and `req()` gating contracts.
- The focused gate now also freezes the shinyFiles setup seam contract in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  for `setupMetabImportShinyFiles()` by asserting the volume fallback,
  chooser registration, assay-1 import trigger delegation, assay-2 parse-error
  logging, and observer wiring contracts.
- The focused gate now also freezes the column-accessor setup seam contract in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  for `setupMetabImportColumnAccessors()` by asserting the custom
  metabolite/annotation resolution delegation, the dropdown bypass contract,
  and the sample-column reactive `req()` and resolver handoff.
- The twenty-second bounded live seam now sits in
  `R/mod_metab_import.R`
  as `setupMetabImportAssaySelectionCallback()`, which now owns the remaining
  local `import_data()` callback shell's current `local_data$assay1_file`
  lookup plus its delegation into `runMetabImportAssaySelection()` before the
  wrapper falls through to the remaining output and observer registration.
- The focused gate now also freezes the `import_data()` callback seam contract
  in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  for `setupMetabImportAssaySelectionCallback()` by asserting the callback
  reads the current assay-1 file at invocation time and preserves the
  `local_data`/`session` delegation handoff.
- The twenty-third bounded live seam now sits in
  `R/mod_metab_import.R`,
  as `setupMetabImportProcessingObserver()`, which now owns the remaining
  `process_import` observer registration shell's current `input`,
  `local_data`, `columnAccessors`, and `workflow_data` lookup before the
  wrapper falls through to the remaining output registration.
- The focused gate now also freezes the `process_import` observer registration
  seam contract in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  for `setupMetabImportProcessingObserver()` by asserting the observer handler
  reads the current import inputs and local assay state at trigger time and
  preserves the `runMetabImportProcessing()` delegation handoff.
- The twenty-fourth bounded live seam now sits in
  `R/mod_metab_import.R`,
  as `setupMetabImportStatusOutput()`, which now owns the remaining
  `import_status` render registration shell's current `workflow_data`
  completion-state and processing-log lookup before the wrapper falls through
  to the remaining output registration.
- The focused gate now also freezes the `import_status` render registration
  seam contract in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  for `setupMetabImportStatusOutput()` by asserting the render callback reads
  the current workflow state at render time and preserves the
  `buildMetabImportStatus()` delegation handoff.
- The twenty-fifth bounded live seam now sits in
  `R/mod_metab_import.R`,
  as `setupMetabImportValidationSummaryOutput()`, which now owns the remaining
  `validation_summary` render registration shell's current
  `local_data$assay1_data` and `columnAccessors` lookup before the wrapper
  falls through to the remaining custom output and observer registration.
- The focused gate now also freezes the `validation_summary` render
  registration seam contract in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  for `setupMetabImportValidationSummaryOutput()` by asserting the render
  callback reads the current assay data and column accessors at render time and
  preserves the `buildMetabImportValidationSummary()` delegation handoff.
- The twenty-sixth bounded live seam now sits in
  `R/mod_metab_import.R`,
  as `setupMetabImportCustomAnnotationStatusOutput()`, which now owns the
  remaining `annotation_status_custom` render registration shell's current
  `local_data$assay1_data` and `input$annotation_col_custom` lookup before the
  wrapper falls through to the remaining validation-summary and observer
  registration.
- The focused gate now also freezes the `annotation_status_custom` render
  registration seam contract in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  for `setupMetabImportCustomAnnotationStatusOutput()` by asserting the render
  callback reads the current assay data and custom annotation column input at
  render time and preserves the
  `buildMetabImportCustomAnnotationStatus()` delegation handoff.
- The twenty-seventh bounded live seam now sits in
  `R/mod_metab_import.R`,
  as `setupMetabImportCustomMetaboliteIdStatusOutput()`, which now owns the
  remaining `metabolite_id_status_custom` render registration shell's current
  `local_data$assay1_data` and `input$metabolite_id_col_custom` lookup before
  the wrapper falls through to the remaining output registration.
- The focused gate now also freezes the `metabolite_id_status_custom` render
  registration seam contract in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  for `setupMetabImportCustomMetaboliteIdStatusOutput()` by asserting the
  render callback reads the current assay data and custom metabolite column
  input at render time and preserves the
  `buildMetabImportCustomMetaboliteIdStatus()` delegation handoff.
- The twenty-eighth bounded live seam now sits in
  `R/mod_metab_import.R`,
  as `setupMetabImportAvailableColumnsDisplayOutput()`, which now owns the
  remaining `available_columns_display` render registration shell's current
  `local_data$all_headers` lookup before the wrapper falls through to the
  remaining output registration.
- The focused gate now also freezes the `available_columns_display` render
  registration seam contract in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  for `setupMetabImportAvailableColumnsDisplayOutput()` by asserting the
  render callback reads the current header vector at render time and preserves
  the `buildMetabImportAvailableColumnsDisplay()` delegation handoff.
- The twenty-ninth bounded live seam now sits in
  `R/mod_metab_import.R`,
  as `setupMetabImportSampleColumnsDisplayOutput()`, which now owns the
  remaining `sample_columns_display` render registration shell's current
  `local_data$assay1_import_result` lookup before the wrapper falls through to
  the remaining inline render registrations.
- The focused gate now also freezes the `sample_columns_display` render
  registration seam contract in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  for `setupMetabImportSampleColumnsDisplayOutput()` by asserting the render
  callback reads the current import result at render time and preserves the
  `buildMetabImportSampleColumnsDisplay()` delegation handoff.
- The thirtieth bounded live seam now sits in
  `R/mod_metab_import.R`,
  as `setupMetabImportAnnotationStatusOutput()`, which now owns the former
  inline `annotation_status` render registration shell's current
  `local_data$assay1_data` and `input$annotation_col` lookup while preserving
  the delegation into `buildMetabImportAnnotationStatus()`.
- The focused gate now also freezes the `annotation_status` render
  registration seam contract in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  for `setupMetabImportAnnotationStatusOutput()` by asserting the render
  callback reads the current assay data and annotation column input at render
  time and preserves the `buildMetabImportAnnotationStatus()` delegation
  handoff.
- The thirty-first bounded live seam now sits in
  `R/mod_metab_import.R`,
  as `setupMetabImportMetaboliteIdStatusOutput()`, which now owns the former
  inline `metabolite_id_status` render registration shell's current
  `local_data$assay1_data` and `input$metabolite_id_col` lookup while
  preserving the delegation into `buildMetabImportMetaboliteIdStatus()`.
- The focused gate now also freezes the `metabolite_id_status` render
  registration seam contract in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  for `setupMetabImportMetaboliteIdStatusOutput()` by asserting the render
  callback reads the current assay data and metabolite column input at render
  time and preserves the `buildMetabImportMetaboliteIdStatus()` delegation
  handoff.
- The thirty-second bounded live seam now sits in
  `R/mod_metab_import.R`,
  as `setupMetabImportFormatDetectionStatusOutput()`, which now owns the
  former inline `format_detection_status` render registration shell's current
  `local_data$detected_format` and `local_data$format_confidence` lookup while
  preserving the delegation into `buildMetabImportFormatDetectionStatus()`.
- The focused gate now also freezes the `format_detection_status` render
  registration seam contract in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  for `setupMetabImportFormatDetectionStatusOutput()` by asserting the render
  callback reads the current detected format and confidence at render time and
  preserves the `buildMetabImportFormatDetectionStatus()` delegation handoff.
- The thirty-third bounded live seam now sits in
  `R/mod_metab_import.R`,
  as `setupMetabImportFileLoadedOutput()`, which now owns the former inline
  `file_loaded` reactive registration shell's current
  `local_data$assay1_data` lookup plus the unsuspended
  `shiny::outputOptions()` handoff.
- The focused gate now also freezes the `file_loaded` reactive registration
  seam contract in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  for `setupMetabImportFileLoadedOutput()` by asserting the reactive callback
  reads the current assay data at evaluation time and preserves the
  unsuspended `file_loaded` output registration via `shiny::outputOptions()`.
- The focused metabolomics import module helper gate reran green after this
  reactive-registration seam via direct `testthat::test_file()` invocation in
  this worktree.
- Post-checkpoint classification now records
  `R/mod_metab_import.R`
  at `review` with `1561` lines, `43` top-level functions, and a `298` line
  largest top-level function.
- The first exact-source metabolomics import output-registration wave now
  lives in
  `tools/refactor/manifest-metab-import-wave2.yml`
  and stages the late-file registration cluster rooted at
  `R/mod_metab_import.R:1286`
  into
  `tools/refactor/staging/wave2_metabolomics_import_output_registrations/R/mod_metab_import_registration_helpers.R`
  without rewriting live sources.
- The staged wave covers
  `setupMetabImportFileLoadedOutput()`,
  `setupMetabImportStatusOutput()`,
  `setupMetabImportValidationSummaryOutput()`,
  `setupMetabImportFormatDetectionStatusOutput()`,
  `setupMetabImportMetaboliteIdStatusOutput()`,
  `setupMetabImportAnnotationStatusOutput()`,
  `setupMetabImportSampleColumnsDisplayOutput()`,
  `setupMetabImportAvailableColumnsDisplayOutput()`,
  `setupMetabImportCustomMetaboliteIdStatusOutput()`,
  and
  `setupMetabImportCustomAnnotationStatusOutput()`.
- `tools/refactor/verify_refactor.R` passed for
  `tools/refactor/manifest-metab-import-wave2.yml`
  before staging the reviewed output-registration helper target.
- The focused metabolomics import module helper gate reran green after the
  staged wave via direct
  `Rscript -e "testthat::test_file(...)"` invocation because this worktree
  does not currently have `renv/activate.R` available for
  `tools/test_with_renv.R`.
- The helper-characterization loader in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  now widens the filename-coupled seam surface by loading
  `R/mod_metab_import_registration_helpers.R`
  before
  `R/mod_metab_import.R`,
  so the focused gate survives the live apply boundary.
- The first exact-source metabolomics import output-registration wave is now
  live in
  `R/mod_metab_import_registration_helpers.R`
  after applying
  `tools/refactor/manifest-metab-import-wave2.yml`;
  the extracted registration helpers are now called from
  `R/mod_metab_import.R`
  instead of remaining duplicated inline.
- `tools/refactor/check_wave_apply.R` passed for
  `tools/refactor/manifest-metab-import-wave2.yml`
  after the live apply.
- `DESCRIPTION` `Collate:` now loads
  `mod_metab_import_registration_helpers.R`
  before
  `mod_metab_import.R`
  so package load order matches the new helper dependency.
- The focused metabolomics import module helper gate reran green after the
  live apply via direct
  `Rscript -e "testthat::test_file(...)"` invocation because this worktree
  still does not have `renv/activate.R` available for
  `tools/test_with_renv.R`.
- The next exact-source metabolomics import server-setup wave now lives in
  `tools/refactor/manifest-metab-import-wave3.yml`
  and stages the late-wrapper setup cluster from
  `R/mod_metab_import.R`
  into
  `tools/refactor/staging/wave3_metabolomics_import_server_setup/R/mod_metab_import_server_helpers.R`
  without rewriting live sources before review.
- The helper-characterization loader in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  now widens the filename-coupled seam surface by loading
  `R/mod_metab_import_server_helpers.R`
  between
  `R/mod_metab_import_registration_helpers.R`
  and
  `R/mod_metab_import.R`,
  so the focused gate survives the live server-helper apply boundary.
- The reviewed server-setup helper wave is now live in
  `R/mod_metab_import_server_helpers.R`
  after applying
  `tools/refactor/manifest-metab-import-wave3.yml`;
  `setupMetabImportShinyFiles()`,
  `setupMetabImportColumnAccessors()`,
  `setupMetabImportAssaySelectionCallback()`,
  and
  `setupMetabImportProcessingObserver()`
  are now called from
  `R/mod_metab_import.R`
  instead of remaining duplicated inline.
- `tools/refactor/check_wave_apply.R` passed for
  `tools/refactor/manifest-metab-import-wave3.yml`
  after the live apply.
- `DESCRIPTION` `Collate:` now loads
  `mod_metab_import_server_helpers.R`
  before
  `mod_metab_import.R`
  so package load order matches the new helper dependency.
- The next exact-source metabolomics import processing/state wave now lives in
  `tools/refactor/manifest-metab-import-wave4.yml`
  and stages the helper cluster into
  `tools/refactor/staging/wave4_metabolomics_import_processing_state/R/mod_metab_import_processing_helpers.R`
  before the reviewed live apply.
- The helper-characterization loader in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  now widens the filename-coupled seam surface by loading
  `R/mod_metab_import_processing_helpers.R`
  between
  `R/mod_metab_import_registration_helpers.R`
  and
  `R/mod_metab_import_server_helpers.R`,
  so the focused gate survives the live processing-helper apply boundary.
- The reviewed processing/state helper wave is now live in
  `R/mod_metab_import_processing_helpers.R`
  after applying
  `tools/refactor/manifest-metab-import-wave4.yml`;
  the extracted processing/state helpers are now called from
  `R/mod_metab_import.R`
  instead of remaining duplicated inline.
- `tools/refactor/check_wave_apply.R` passed for
  `tools/refactor/manifest-metab-import-wave4.yml`
  after the live apply.
- `DESCRIPTION` `Collate:` now loads
  `mod_metab_import_processing_helpers.R`
  before
  `mod_metab_import_server_helpers.R`
  and
  `mod_metab_import.R`
  so package load order matches the new processing-helper dependency.
- The focused metabolomics import module helper gate reran green after the
  live apply via direct
  `Rscript -e "testthat::test_file(...)"` invocation because this worktree
  still does not have `renv/activate.R` available for
  `tools/test_with_renv.R`,
  and the adjacent metabolomics import detection guardrail also reran green.
- The reviewed display/status helper wave is now live in
  `R/mod_metab_import_display_helpers.R`
  after applying
  `tools/refactor/manifest-metab-import-wave5.yml`;
  the extracted file-selection and display/status helpers are now called from
  `R/mod_metab_import.R`
  instead of remaining duplicated inline.
- `tools/refactor/check_wave_apply.R` passed for
  `tools/refactor/manifest-metab-import-wave5.yml`
  after the live apply.
- `DESCRIPTION` `Collate:` now loads
  `mod_metab_import_display_helpers.R`
  before
  `mod_metab_import_registration_helpers.R`,
  `mod_metab_import_processing_helpers.R`,
  `mod_metab_import_server_helpers.R`,
  and
  `mod_metab_import.R`
  so package load order matches the new display/status helper dependency.
- The helper-characterization loader in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  now widens the filename-coupled seam surface by loading
  `R/mod_metab_import_display_helpers.R`
  before the existing registration, processing, server, and wrapper files, so
  the focused gate survives the live display/status helper apply boundary.
- The focused metabolomics import module helper gate reran green after the
  live apply via direct
  `Rscript -e "testthat::test_file(...)"` invocation because this worktree
  still does not have `renv/activate.R` available for
  `tools/test_with_renv.R`,
  and the adjacent metabolomics import detection guardrail also reran green.
- Post-checkpoint classification now records
  `R/mod_metab_import.R`
  at `review` with `517` lines, `4` top-level functions, and a `298` line
  largest top-level function.
- The final exact-source metabolomics import wrapper-reduction wave is now
  live after applying
  `tools/refactor/manifest-metab-import-wave6.yml`;
  `resolveMetabImportColumnName()` and
  `resolveMetabImportSampleColumns()` now live in
  `R/mod_metab_import_column_helpers.R`,
  and the two public entry points now live in
  `R/mod_metab_import_ui.R`
  and
  `R/mod_metab_import_server.R`.
- `tools/refactor/check_wave_apply.R` passed for
  `tools/refactor/manifest-metab-import-wave6.yml`
  after the live apply.
- `DESCRIPTION` `Collate:` now loads
  `mod_metab_import_column_helpers.R`,
  `mod_metab_import_ui.R`,
  `mod_metab_import_server_helpers.R`,
  `mod_metab_import_server.R`,
  and then
  `mod_metab_import.R`
  so package load order matches the final wrapper split.
- The helper-characterization loader in
  `tests/testthat/test-metab-00b-import-module-helper-characterization.R`
  now widens the filename-coupled seam surface by loading
  `R/mod_metab_import_column_helpers.R`
  before the existing display, registration, processing, server-helper, and
  wrapper files, so the focused gate survives the final wrapper apply
  boundary.
- The focused metabolomics import module helper gate reran green after the
  live apply via direct
  `Rscript -e "testthat::test_file(...)"` invocation because this worktree
  still does not have `renv/activate.R` available for
  `tools/test_with_renv.R`,
  and the adjacent metabolomics import detection guardrail also reran green.
- Post-checkpoint classification now records
  `R/mod_metab_import.R`
  at `direct-extraction-ready` with `66` lines, `0` top-level functions, and
  a `0` line largest top-level function.
- Treat the metabolomics import wrapper target as complete:
  `R/mod_metab_import.R`
  is now the breadcrumb stub for the public module identity, and this
  handover can be considered archival for the completed manual target.
