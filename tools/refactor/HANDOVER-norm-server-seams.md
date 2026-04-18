# Normalization Server Seam Map

## Goal

Start stabilization of the proteomics normalization wrapper while keeping
[mod_prot_norm_server()](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm_server.R:9)
behavior frozen behind the existing public entry points.

## Current Position In The Flow

- Proteomics import stabilization is complete and archived in
  [HANDOVER-import-server-seams.md](/home/doktersmol/Documents/MultiScholaR/tools/refactor/HANDOVER-import-server-seams.md:1).
- Survey and classification were rerun on April 11, 2026 for the normalization
  pair:
  - [mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:1)
    is `high-risk-wrapper` and `needs-seam-introduction`
  - [func_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_norm.R:1)
    is `direct-extraction-ready`
- Classification was rerun again on April 11, 2026 after the wave 3 apply:
  - [mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:1)
    remains `high-risk-wrapper` and `needs-seam-introduction`
- Classification was rerun again on April 11, 2026 after the QC image seam:
  - [mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:1)
    still remains `high-risk-wrapper` and `needs-seam-introduction`
- Classification was rerun again on April 11, 2026 after the render-shell seam:
  - [mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:1)
    still remains `high-risk-wrapper` and `needs-seam-introduction`
- Classification was rerun again on April 12, 2026 after the tab-entry and
  normalization-shell seams:
  - [mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:1)
    remains `high-risk-wrapper` and `needs-seam-introduction`
- Classification was rerun again on April 12, 2026 after the late-observer
  seam:
  - [mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:1)
    remains `high-risk-wrapper` and `needs-seam-introduction`
- Classification was rerun again on April 12, 2026 after the wave 4 entrypoint
  apply:
  - [R/mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:1)
    is now `direct-extraction-ready`
  - [R/mod_prot_norm_server.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm_server.R:1)
    is now a compact `review` wrapper with no inline observers
- The first normalization seam is now introduced in-place in
  [mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:488)
  with:
  - `buildProtNormQcImagePayload()`
  - `getProtNormFilteringSummaryText()`
  - `buildProtNormRuvOptimizationSummary()`
  - `prepareProtNormOptimizationResultsTable()`
- The second bottom-render/default seam is now also introduced in-place in
  [mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:577)
  with:
  - `resolveProtNormFinalQcRenderState()`
  - `getProtNormRuvCanonicalCorrelationPlot()`
  - `getProtNormDefaultCorrelationFilterSummaryText()`
- The third plot-aesthetic/grouping seam is now also introduced in-place in
  [mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:615)
  with:
  - `updateProtNormDesignDrivenChoices()`
  - `getProtNormPlotAesthetics()`
  - `getProtNormRuvGroupingVariable()`
- The fourth composite-builder seam is now also introduced in-place in
  [mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:684)
  with:
  - `buildProtNormLabelPlot()`
  - `buildProtNormTitlePlot()`
  - `loadProtNormImageAsPlot()`
  - `generateProtNormCompositeFromFiles()`
- The fifth pre-observer decision seam is now also introduced in-place in
  [mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:822)
  with:
  - `shouldProtNormAutoGeneratePreQc()`
  - `regenerateProtNormQcForAestheticChange()`
- The sixth QC generation support seam is now also introduced in-place in
  [mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:868)
  with:
  - `initializeProtNormQcPlotPaths()`
  - `saveProtNormQcPlotArtifact()`
  - `recordProtNormQcPlotPath()`
  - `buildProtNormDensityPlot()`
  - `buildProtNormCorrelationPlot()`
- The seventh QC generation wrapper seam is now also introduced in-place in
  [mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:962)
  with:
  - `resolveProtNormQcStateObject()`
  - `generateProtNormPreNormalizationQcArtifacts()`
  - `generateProtNormPostNormalizationQcArtifacts()`
  - `generateProtNormRuvCorrectedQcArtifacts()`
- The first normalization helper apply wave is now live via
  [manifest-norm-server-wave1.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-norm-server-wave1.yml:1):
  - [R/mod_prot_norm_qc_support_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm_qc_support_helpers.R:1)
  - [R/mod_prot_norm_qc_generation_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm_qc_generation_helpers.R:1)
- The next bounded main-observer seam is now introduced in-place in
  [mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:491)
  with:
  - `prepareProtNormNormalizationRun()`
  - `runProtNormBetweenSamplesStep()`
  - `runProtNormPostNormalizationQcStep()`
  - `buildProtNormSkippedRuvResult()`
  - `applyProtNormSkippedRuvState()`
- The next RUV-core observer seam is now introduced in-place in
  [mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:638)
  with:
  - `buildProtNormManualRuvResult()`
  - `updateProtNormRuvAuditTrail()`
  - `persistProtNormRuvResult()`
  - `resolveProtNormRuvParameters()`
  - `applyProtNormRuvCorrectionStep()`
  - `finalizeProtNormRuvCleanupStep()`
- The remaining step-6/completion shell is now also introduced in-place in
  [mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:1019)
  with:
  - `resolveProtNormStep6QcObject()`
  - `runProtNormStep6RuvQc()`
  - `resolveProtNormCompositeFigureInputs()`
  - `generateProtNormCompositeQcFigure()`
  - `finalizeProtNormWorkflowState()`
  - `buildProtNormCompletionNotification()`
- The second normalization helper apply wave is now live via
  [manifest-norm-server-wave2.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-norm-server-wave2.yml:1):
  - [R/mod_prot_norm_workflow_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm_workflow_helpers.R:1)
  - [R/mod_prot_norm_ruv_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm_ruv_helpers.R:1)
- The first correlation-completion seam is now introduced in-place in
  [mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:508)
  with:
  - `resolveProtNormCorrelationResultFilenames()`
  - `saveProtNormCorrelationResults()`
  - `updateProtNormFinalFilteringPlot()`
  - `updateProtNormFinalQcPlot()`
  - `finalizeProtNormCorrelationWorkflowState()`
  - `buildProtNormCorrelationSummaryText()`
- The next core correlation seam is now also introduced in-place in
  [mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:754)
  with:
  - `resolveProtNormCorrelationInputObject()`
  - `runProtNormCorrelationVectorStep()`
  - `runProtNormCorrelationFilterStep()`
  - `prepareProtNormSkippedCorrelationState()`
- The outer correlation observer shell seam is now also introduced in-place in
  [mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:859)
  with:
  - `runProtNormApplyCorrelationWorkflow()`
  - `runProtNormSkipCorrelationWorkflow()`
  - `completeProtNormCorrelationWorkflow()`
  - `handleProtNormCorrelationError()`
- The export session seam is now also introduced in-place in
  [mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:1068)
  with:
  - `canProtNormExportFilteredSession()`
  - `resolveProtNormExportSourceDir()`
  - `collectProtNormExportSessionData()`
  - `buildProtNormExportSummaryContent()`
  - `saveProtNormExportMetadataFiles()`
  - `saveProtNormExportArtifacts()`
  - `runProtNormExportSessionWorkflow()`
  - `handleProtNormExportError()`
- The reset normalization seam is now also introduced in-place in
  [mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:1334)
  with:
  - `resolveProtNormPreNormalizationState()`
  - `revertProtNormStateManagerToPreNormalization()`
  - `resetProtNormReactiveState()`
  - `resetProtNormOutputs()`
  - `buildProtNormResetNotificationMessage()`
  - `runProtNormResetWorkflow()`
  - `handleProtNormResetError()`
- The third normalization helper apply wave is now live via
  [manifest-norm-server-wave3.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-norm-server-wave3.yml:1):
  - [R/mod_prot_norm_correlation_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm_correlation_helpers.R:1)
  - [R/mod_prot_norm_session_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm_session_helpers.R:1)
- The QC image registration seam is now also introduced in-place in
  [mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:537)
  with:
  - `renderProtNormQcImage()`
  - `registerProtNormQcImageOutputs()`
- The render output registration seam is now also introduced in-place in
  [mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:634)
  with:
  - `renderProtNormPostNormFilteringSummary()`
  - `renderProtNormFilteringSummaryText()`
  - `renderProtNormFinalQcPlot()`
  - `renderProtNormRuvCanonicalCorrelationPlot()`
  - `renderProtNormRuvOptimizationSummary()`
  - `renderProtNormRuvOptimizationTable()`
  - `registerProtNormRenderOutputs()`
- The tab-entry observer seam is now also introduced in-place in
  [mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:826)
  with:
  - `notifyProtNormNormalizationPrereqWarning()`
  - `handleProtNormPreQcGenerationError()`
  - `runProtNormTabEntryWorkflow()`
- The remaining normalization observer shell is now also introduced in-place in
  [mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:913)
  with:
  - `handleProtNormNormalizationError()`
  - `runProtNormNormalizationWorkflow()`
- The late correlation/export/reset observer seam is now also introduced
  in-place in
  [mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:1093)
  with:
  - `runProtNormApplyCorrelationObserver()`
  - `runProtNormSkipCorrelationObserver()`
  - `notifyProtNormExportSessionPrereqWarning()`
  - `runProtNormExportObserver()`
  - `runProtNormResetObserver()`
- The fourth normalization helper apply wave is now live via
  [manifest-norm-server-wave4.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-norm-server-wave4.yml:1):
  - [R/mod_prot_norm_ui.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm_ui.R:1)
  - [R/mod_prot_norm_server_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm_server_helpers.R:1)
  - [R/mod_prot_norm_server.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm_server.R:1)
- [R/mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:1)
  is now a `94` line breadcrumb stub with no live functions
- The helperized normalization workflow now lives behind:
  - [R/mod_prot_norm_qc_support_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm_qc_support_helpers.R:1)
  - [R/mod_prot_norm_qc_generation_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm_qc_generation_helpers.R:1)
  - [R/mod_prot_norm_workflow_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm_workflow_helpers.R:1)
  - [R/mod_prot_norm_ruv_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm_ruv_helpers.R:1)
  - [R/mod_prot_norm_correlation_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm_correlation_helpers.R:1)
  - [R/mod_prot_norm_session_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm_session_helpers.R:1)
  - [R/mod_prot_norm_server_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm_server_helpers.R:1)
  - [R/mod_prot_norm_server.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm_server.R:1)
  - [R/mod_prot_norm_ui.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm_ui.R:1)
- Direct helper contracts now live in
  [test-prot-05b-norm-module-contracts.R](/home/doktersmol/Documents/MultiScholaR/tests/testthat/test-prot-05b-norm-module-contracts.R:1).

## Existing Safety Net

Current normalization-focused coverage:

- [test-prot-05-normalisation.R](/home/doktersmol/Documents/MultiScholaR/tests/testthat/test-prot-05-normalisation.R:1)
- [test-prot-05b-norm-module-contracts.R](/home/doktersmol/Documents/MultiScholaR/tests/testthat/test-prot-05b-norm-module-contracts.R:1)
- [test-prot-06-ruv.R](/home/doktersmol/Documents/MultiScholaR/tests/testthat/test-prot-06-ruv.R:1)

These currently cover:

- `scaleCenterAndFillMissing()`
- `findBestK()`
- `getRuvIIIReplicateMatrixHelper()`
- `getNegCtrlProtAnovaHelper()`
- QC image payload helper contract
- filtering summary text helper contract
- RUV optimization summary helper contract
- optimization table preparation helper contract
- final QC render-state selector helper contract
- canonical correlation plot selector helper contract
- correlation-summary default helper contract
- design-driven choice update helper contract
- plot-aesthetic fallback helper contract
- RUV grouping fallback helper contract
- label-plot helper contract
- title-plot helper contract
- image-loader helper contract
- composite-builder helper contract
- pre-QC auto-trigger decision helper contract
- aesthetic-change regeneration helper contract
- QC path initialization helper contract
- QC plot artifact saver helper contract
- QC plot path recorder helper contract
- QC density-plot helper contract
- QC correlation-plot helper contract
- QC state resolver helper contract
- pre-normalization QC wrapper helper contract
- post-normalization QC wrapper helper contract
- RUV-corrected QC wrapper helper contract
- normalization run bootstrap helper contract
- between-samples normalization step helper contract
- post-normalization QC step helper contract
- skipped-RUV result helper contract
- skipped-RUV state application helper contract
- manual-RUV result builder helper contract
- RUV audit-trail helper contract
- RUV result persistence helper contract
- RUV parameter resolution helper contract
- RUV correction apply helper contract
- post-RUV cleanup/state-save helper contract
- step-6 QC object resolution helper contract
- step-6 RUV QC/cancor helper contract
- composite-figure input resolver helper contract
- composite-figure generation helper contract
- workflow-finalization helper contract
- completion-notification helper contract
- correlation-result filename helper contract
- correlation-result persistence helper contract
- final filtering-plot update helper contract
- final QC-plot update helper contract
- correlation workflow-finalization helper contract
- correlation summary-text helper contract
- correlation input-object helper contract
- correlation vector step helper contract
- correlation filter step helper contract
- skipped-correlation state helper contract
- apply-correlation workflow shell helper contract
- skip-correlation workflow shell helper contract
- correlation workflow completion helper contract
- correlation error handler helper contract
- export-readiness helper contract
- export source-dir resolver helper contract
- export session-data collector helper contract
- export summary builder helper contract
- export artifact saver helper contract
- export workflow runner helper contract
- export error handler helper contract
- pre-normalization reset-state resolver helper contract
- state-manager reset/revert helper contract
- normalization reactive reset helper contract
- normalization reset output helper contract
- normalization reset workflow shell helper contract
- normalization reset error handler helper contract
- QC image render helper contract
- QC image output-registration helper contract
- post-normalization filtering-summary render helper contract
- filtering-summary text render helper contract
- final-QC render helper contract
- canonical-correlation render helper contract
- RUV summary render helper contract
- RUV optimization-table render helper contract
- remaining render-output registration helper contract
- normalization-entry prereq warning helper contract
- pre-normalization QC error handler helper contract
- tab-entry auto-pre-QC workflow helper contract
- normalization workflow error handler helper contract
- normalization workflow shell helper contract
- apply-correlation observer shell helper contract
- skip-correlation observer shell helper contract
- export prereq warning helper contract
- export observer shell helper contract
- reset observer shell helper contract

## Gate Status

The focused normalization gate was rerun on April 12, 2026 after the
late-observer seam and stayed green:

- `test-prot-05-normalisation.R` passed with one skip because
  `cp05_normalised.rds` is a Git LFS pointer in this checkout
- `test-prot-05b-norm-module-contracts.R` passed
- `test-prot-06-ruv.R` passed with one skip because
  `cp06_ruv_corrected.rds` is a Git LFS pointer in this checkout

## Why The Next Split Is Different

- [mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:1)
  is still a large live wrapper at `1550` lines with one dominant server body
  and multiple nested observers, but its largest top-level function is now
  down to `458` lines.
- [func_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/func_prot_norm.R:1)
  is already direct-extraction-ready and should not drive the first wrapper
  seams.
- The first three helper waves are now live, so the next work shifts from
  helper extraction to deciding whether the remaining observer shells
  should be helperized further or whether normalization is now "small enough"
  relative to the rest of the backlog.

## Candidate Internal Seams

### 1. Bottom QC and RUV render helpers

Current state:

- complete and green

Responsibility:

- QC image payload resolution
- filtering summary fallback text
- RUV optimization summary string formatting
- RUV optimization table preparation
- final QC render-state selection
- canonical correlation plot selection
- reset-time correlation-summary default text

### 2. Remaining bottom renderers and reset defaults

Responsibility:

- final QC plot rendering shell
- cancor plot placeholder/render shell
- reset-time default render text shells

Current state:

- complete and green

### 3. Plot-aesthetic and grouping-choice helpers

Responsibility:

- design-driven choice updates
- plot aesthetic fallback resolution
- RUV grouping variable fallback resolution

Current state:

- complete and green

### 4. Main normalization workflow observer

Responsibility:

- normalization dispatch
- post-normalization QC generation
- RUV optimization and apply flow
- state-manager updates
- composite QC generation

Current state:

- the composite builder, the pre-observer auto-trigger/regeneration decisions,
  the shared QC generation support/wrapper helpers, and the full
  `run_normalization` helper cluster are now extracted behind
  [R/mod_prot_norm_qc_support_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm_qc_support_helpers.R:1)
  ,
  [R/mod_prot_norm_qc_generation_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm_qc_generation_helpers.R:1)
  ,
  [R/mod_prot_norm_workflow_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm_workflow_helpers.R:1)
  and
  [R/mod_prot_norm_ruv_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm_ruv_helpers.R:1)
- [runProtNormNormalizationWorkflow()](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:928)
  now owns the main normalization observer shell, and
  [observeEvent(input$run_normalization, ...)](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:1262)
  is now just a wrapper call
- [runProtNormTabEntryWorkflow()](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:859)
  now owns the normalization-tab auto-pre-QC shell, and the corresponding
  tab-selection observer is no longer the next live seam target
- [runProtNormApplyCorrelationObserver()](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:1093),
  [runProtNormSkipCorrelationObserver()](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:1154),
  [runProtNormExportObserver()](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:1225),
  and
  [runProtNormResetObserver()](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:1274)
  now own the late observer shells, and the corresponding button observers are
  now wrapper calls

Next safe target:

- introduce one more top-level helper seam in
  [mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:1),
  now most likely around the remaining observer/render registration shell and
  local QC generator closures in
  [mod_prot_norm_server()](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:1308)

### 5. Correlation, reset, and export observers

Current target:

- completed and applied live via
  [manifest-norm-server-wave3.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-norm-server-wave3.yml:1)

Responsibility:

- correlation filtering apply/skip branches
- post-filter state updates
- normalization reset and state-manager rollback flow
- export filtered session flow

Current state:

- normalization wrapper stabilization is complete in live `R/`
- the public entrypoints now live in dedicated files:
  - [R/mod_prot_norm_ui.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm_ui.R:1)
  - [R/mod_prot_norm_server.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm_server.R:1)
- the remaining observer shell now delegates through
  `registerProtNormServerObservers()` in
  [R/mod_prot_norm_server_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm_server_helpers.R:1)
- the target file [R/mod_prot_norm.R](/home/doktersmol/Documents/MultiScholaR/R/mod_prot_norm.R:1)
  is archived as a breadcrumb stub and no longer blocks proteomics work
- normalization gate reruns after the wave 4 apply stayed green:
  - `test-prot-05-normalisation`
  - `test-prot-05b-norm-module-contracts`
  - `test-prot-06-ruv`

## Safe Next Step

Move to the next proteomics bucket: peptide/QC-rollup.

## Resume Commands

Normalization archive checkpoint:

```bash
Rscript tools/test_with_renv.R tests/testthat/test-prot-05-normalisation.R
Rscript tools/test_with_renv.R tests/testthat/test-prot-05b-norm-module-contracts.R
Rscript tools/test_with_renv.R tests/testthat/test-prot-06-ruv.R
```

Then continue with the next backlog target rather than resuming normalization.

```text
Use $god-module-stabilization to move on from archived proteomics normalization.
Treat tools/refactor/HANDOVER-norm-server-seams.md as complete and archived, keep the normalization gate as the regression surface for future changes, and pick up the next proteomics backlog item instead of reopening mod_prot_norm.R.
```

## Guardrails

- keep helper names in camelCase
- treat unresolved Git LFS snapshot fixtures as skips, not behavior regressions
- future normalization edits should touch the dedicated helper/server files, not
  re-expand `mod_prot_norm.R`
- this handover is archival after the green wave 4 apply
