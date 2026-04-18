# Lipid Design Builder Seams

## Goal

Stabilize `R/mod_lipid_design_builder.R` with bounded exact-source helper waves while preserving the public module entrypoints.

## Current Position In The Flow

- target: `R/mod_lipid_design_builder.R`
- classification: `direct-extraction-ready`
- next step: archive this handover; the live helper/UI/server waves are applied and the wrapper is now a breadcrumb stub

## Existing Safety Net

- `Rscript -e "source('R/mod_lipid_design_builder_display_helpers.R'); source('R/mod_lipid_design_builder_action_helpers.R'); source('R/mod_lipid_design_builder_state_helpers.R'); source('R/mod_lipid_design_builder_ui.R'); source('R/mod_lipid_design_builder_server.R'); source('R/mod_lipid_design_builder.R'); testthat::test_file('tests/testthat/test-lipid-13-design-builder-seams.R', stop_on_failure = TRUE)"`

## Notes

- April 15, 2026 stabilize-mode checkpoint introduced the first bounded summary-output seam in `R/mod_lipid_design_builder.R` with:
  - `formatLipidDesignTechRepSummary()`
  - `formatLipidDesignRemovedSamplesDisplay()`
  - `formatLipidDesignContrastFactorsInfo()`
  - `registerLipidDesignSummaryOutputShells()`
- The live wrapper now routes `tech_rep_summary`, `removed_samples_display`, and `contrast_factors_info` through that helper seam while keeping the public module entrypoints unchanged.
- April 15, 2026 stabilize-mode checkpoint introduced the second bounded adjacent output seam in `R/mod_lipid_design_builder.R` with:
  - `formatLipidDesignAvailableFactorsDisplay()`
  - `formatLipidDesignDefinedContrastLines()`
  - `formatLipidDesignRangePreview()`
  - `buildLipidDesignReplicateInputLabel()`
  - `registerLipidDesignAdjacentOutputShells()`
- The live wrapper now routes `available_factors_display`, `defined_contrasts_display`, `range_preview`, and `replicate_inputs` through that helper seam while keeping the public module entrypoints unchanged.
- April 15, 2026 stabilize-mode checkpoint introduced the third bounded sample-edit bookkeeping seam in `R/mod_lipid_design_builder.R` with:
  - `renameLipidDesignDesignMatrixRuns()`
  - `renameLipidDesignAssayColumns()`
  - `renameLipidDesignTrackedSamples()`
  - `applyLipidDesignSampleRenameMap()`
  - `buildLipidDesignBulkRenameMap()`
  - `registerLipidDesignSampleRenameShells()`
- The live wrapper now routes `rename_sample` and `bulk_rename` through that helper seam while keeping the public module entrypoints unchanged.
- April 15, 2026 stabilize-mode checkpoint introduced the fourth bounded metadata-assignment seam in `R/mod_lipid_design_builder.R` with:
  - `appendLipidDesignFactorName()`
  - `buildLipidDesignMetadataReplicateNumbers()`
  - `buildLipidDesignMetadataGroupName()`
  - `listLipidDesignAssignedGroups()`
  - `applyLipidDesignMetadataAssignment()`
  - `registerLipidDesignMetadataAssignmentShells()`
- The live wrapper now routes `add_factor` and `assign_metadata` through that helper seam while keeping the public module entrypoints unchanged.
- April 15, 2026 stabilize-mode checkpoint introduced the fifth bounded technical-replicate seam in `R/mod_lipid_design_builder.R` with:
  - `applyLipidDesignTechRepAssignment()`
  - `registerLipidDesignTechRepAssignmentShells()`
- The live wrapper now routes `assign_tech_reps` through that helper seam while keeping the public module entrypoints unchanged.
- April 15, 2026 stabilize-mode checkpoint introduced the sixth bounded contrast-management seam in `R/mod_lipid_design_builder.R` with:
  - `appendLipidDesignContrast()`
  - `registerLipidDesignContrastManagementShells()`
- The live wrapper now routes `add_contrast` through that helper seam while keeping the public module entrypoints unchanged.
- April 15, 2026 stabilize-mode checkpoint introduced the seventh bounded sample-removal seam in `R/mod_lipid_design_builder.R` with:
  - `appendLipidDesignRemovedSamples()`
  - `registerLipidDesignSampleRemovalShells()`
- The live wrapper now routes `remove_samples` through that helper seam while keeping the public module entrypoints unchanged.
- April 15, 2026 stabilize-mode checkpoint introduced the eighth bounded reset-request seam in `R/mod_lipid_design_builder.R` with:
  - `buildLipidDesignResetConfirmationBody()`
  - `buildLipidDesignResetConfirmationModal()`
  - `registerLipidDesignResetRequestShells()`
- The live wrapper now routes `reset_changes` through that helper seam while keeping the public module entrypoints unchanged.
- April 15, 2026 stabilize-mode checkpoint introduced the ninth bounded reset-confirmation seam in `R/mod_lipid_design_builder.R` with:
  - `applyLipidDesignResetState()`
  - `runLipidDesignResetConfirmationShell()`
  - `registerLipidDesignResetConfirmationShells()`
- The live wrapper now routes `confirm_reset` through that helper seam while keeping the public module entrypoints unchanged.
- April 15, 2026 stabilize-mode checkpoint introduced the tenth bounded save-results seam in `R/mod_lipid_design_builder.R` with:
  - `buildLipidDesignSaveResultsDesignMatrix()`
  - `buildLipidDesignSaveResultsDataList()`
  - `buildLipidDesignSaveResultsContrastsTable()`
  - `buildLipidDesignSaveResultsPayload()`
  - `runLipidDesignSaveResultsShell()`
  - `registerLipidDesignSaveResultsShells()`
- The live wrapper now routes `save_results` through that helper seam while keeping the public module entrypoints unchanged.
- The focused lipid design seam gate now also characterizes the sample-edit bookkeeping helpers, metadata-assignment helpers, technical-replicate helper/shell, contrast-management helper/shell, sample-removal helper/shell, and reset-request helper/shell so the current stop point is frozen before any staged wave planning.
- The focused lipid design seam gate now also characterizes the reset-confirmation state-apply helper and observer shell so the adjacent reset cluster is frozen before moving on.
- The focused lipid design seam gate now also characterizes the save-results design-matrix filter, assay-column selection, payload assembly, and observer shell so the final live save path is frozen before staged-wave review.
- April 15, 2026 stabilize-mode checkpoint introduced the eleventh bounded factor/group dropdown observer-registration seam in `R/mod_lipid_design_builder.R` with:
  - `runLipidDesignFactorGroupDropdownShell()`
  - `registerLipidDesignFactorGroupDropdownShells()`
- The live wrapper now routes `factor1_select`, `factor2_select`, `factor3_select`, `contrast_group1`, and `contrast_group2` refreshes through that helper seam while keeping the public module entrypoints unchanged.
- The focused lipid design seam gate now also characterizes the factor/group dropdown shell and observer handoff so the first remaining inline selection-refresh observer is frozen before moving on.
- April 15, 2026 stabilize-mode checkpoint introduced the twelfth bounded sample-selection input observer-registration seam in `R/mod_lipid_design_builder.R` with:
  - `runLipidDesignSampleSelectionInputShell()`
  - `registerLipidDesignSampleSelectionInputShells()`
- The live wrapper now routes `sample_to_rename`, `selected_runs`, `samples_to_transform`, `tech_rep_samples`, and `samples_to_remove` refreshes through that helper seam while keeping the public module entrypoints unchanged.
- The focused lipid design seam gate now also characterizes the sample-selection input shell and observer handoff so the remaining inline selection-refresh cluster is frozen before any staged-wave review.
- April 15, 2026 stabilize-mode checkpoint introduced the thirteenth bounded initial-state reset observer seam in `R/mod_lipid_design_builder.R` with:
  - `runLipidDesignInitialStateShell()`
  - `registerLipidDesignInitialStateShells()`
- The live wrapper now routes initial-state reactive resets plus the `sample_to_rename`, `selected_runs`, `samples_to_transform`, `tech_rep_samples`, `samples_to_remove`, and `formula_string` rehydration through that helper seam while keeping the public module entrypoints unchanged.
- The focused lipid design seam gate now also characterizes the initial-state reset shell and bind-event handoff so the reset-triggered state rehydration path is frozen before moving on.
- April 15, 2026 stabilize-mode checkpoint introduced the fourteenth bounded data-table render/refresh seam in `R/mod_lipid_design_builder.R` with:
  - `buildLipidDesignActiveDataTable()`
  - `registerLipidDesignDataTableShells()`
- The live wrapper now routes `data_table` rendering plus proxy refresh replacement through that helper seam while keeping the public module entrypoints unchanged.
- The focused lipid design seam gate now also characterizes the active-table filter helper plus the render/refresh shell handoff so the remaining inline data-table path is frozen before reopening staged-wave review.
- April 15, 2026 stabilize-mode checkpoint introduced the fifteenth bounded initial-state builder seam in `R/mod_lipid_design_builder.R` with:
  - `getLipidDesignSampleColumns()`
  - `buildLipidDesignContrastState()`
  - `buildLipidDesignInitialState()`
- The live wrapper now routes sample-column detection plus imported-versus-fresh initial-state assembly through that helper seam while keeping the public module entrypoints unchanged.
- The focused lipid design seam gate now also characterizes sample-column detection and imported/fresh initial-state assembly so the last nested initialization cluster is frozen before staged-wave review resumes.
- April 15, 2026 stabilize-mode checkpoint introduced the first bounded staged-wave review artifact set for `R/mod_lipid_design_builder.R` with:
  - `tools/refactor/manifest-lipid-design-builder-wave1.yml`
  - `tools/refactor/staging/wave1_lipidomics_design_builder_pre_server_helpers/R/mod_lipid_design_builder_display_helpers.R`
  - `tools/refactor/staging/wave1_lipidomics_design_builder_pre_server_helpers/R/mod_lipid_design_builder_action_helpers.R`
  - `tools/refactor/staging/wave1_lipidomics_design_builder_pre_server_helpers/R/mod_lipid_design_builder_state_helpers.R`
- The staged wave verifies and materializes the accumulated pre-server helper cluster without changing live `R/mod_lipid_design_builder.R`, leaving the wrapper at `1909` lines while the staged helper files land at `366`, `428`, and `658` lines.
- The staged wave collate artifact now exists at `tools/refactor/staging/wave1_lipidomics_design_builder_pre_server_helpers/tools/refactor/collate-lipid-design-builder-wave1.txt`.
- The focused lipid design seam gate reran green after the staged-wave checkpoint via direct `testthat::test_file(...)` with `202` passes.
- Post-checkpoint `classify_target.py --json` still reports `review` at `1909` lines and `56` top-level functions, so any staged extraction should remain review-first in this worktree.
- April 15, 2026 stabilize-mode checkpoint applied the first bounded pre-server helper wave live via `tools/refactor/manifest-lipid-design-builder-wave1.yml` into:
  - `R/mod_lipid_design_builder_display_helpers.R`
  - `R/mod_lipid_design_builder_action_helpers.R`
  - `R/mod_lipid_design_builder_state_helpers.R`
- The live collate artifact now exists at `tools/refactor/collate-lipid-design-builder-wave1.txt`, and `DESCRIPTION` now collates the three new helper files immediately before `R/mod_lipid_design_builder.R`.
- Live `R/mod_lipid_design_builder.R` now measures `507` lines with `2` top-level functions while the helper files hold the extracted pre-server seams at `366`, `428`, and `658` lines respectively.
- The focused lipid design seam gate reran green after the live apply once the source-based test harness loaded the new helper files first, and `check_wave_apply.R --manifest tools/refactor/manifest-lipid-design-builder-wave1.yml` passed on the live rewrite.
- Post-apply `classify_target.py --json` still reports `review`, now at `507` lines and `2` top-level functions, so the remaining work is the public wrapper/UI/server surface rather than more inline observer seam carving.
- `tools/test_with_renv.R` is not usable in this worktree right now because `renv/activate.R` is absent, so the focused gate currently reruns via direct `testthat::test_file(...)`.
- Keep future work on this lane review-first around the remaining wrapper/UI/server surface; do not create ad hoc live helper files outside exact-source waves.
- April 15, 2026 stabilize-mode checkpoint introduced the second bounded staged-wave review artifact set for `R/mod_lipid_design_builder.R` with:
  - `tools/refactor/manifest-lipid-design-builder-wave2.yml`
  - `tools/refactor/staging/wave2_lipidomics_design_builder_entrypoints/R/mod_lipid_design_builder_ui.R`
  - `tools/refactor/staging/wave2_lipidomics_design_builder_entrypoints/R/mod_lipid_design_builder_server.R`
- The staged wave verifies and materializes the remaining public entrypoint surface into dedicated UI and server files without changing live `R/mod_lipid_design_builder.R`, leaving the live wrapper at `507` lines while the staged entrypoint files land at `240` and `172` lines.
- The staged wave collate artifact now exists at `tools/refactor/staging/wave2_lipidomics_design_builder_entrypoints/collate-lipid-design-builder-wave2.txt`.
- The focused lipid design seam gate reran green after the staged-wave checkpoint via direct `testthat::test_file(...)` with `202` passes.
- April 15, 2026 stabilize-mode checkpoint applied the second bounded entrypoint wave live via `tools/refactor/manifest-lipid-design-builder-wave2.yml` into:
  - `R/mod_lipid_design_builder_ui.R`
  - `R/mod_lipid_design_builder_server.R`
- The live collate artifact now exists at `tools/refactor/collate-lipid-design-builder-wave2.txt`, `DESCRIPTION` now collates the UI and server entrypoints immediately before `R/mod_lipid_design_builder.R`, and the focused source-based gate now loads the helper files plus the extracted UI/server entrypoints before the wrapper breadcrumb.
- Post-apply `check_wave_apply.R --manifest tools/refactor/manifest-lipid-design-builder-wave2.yml` passed, and `scripts/classify_target.py --json R/mod_lipid_design_builder.R` now reports `direct-extraction-ready` at `97` lines and `0` top-level functions.
- The public lipid design builder surface now lives in dedicated helper/UI/server files while `R/mod_lipid_design_builder.R` remains a breadcrumb stub, so this target no longer needs additional stabilization work unless later doc regeneration is requested in a dependency-complete environment.
