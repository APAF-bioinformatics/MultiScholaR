# Lipid QC Duplicates Module Seam Map

## Goal

Document the wrapper-level stabilization stop point for
[mod_lipid_qc_duplicates.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_duplicates.R:1)
while keeping the public duplicate-resolution module contract stable.

## Current Position In The Flow

- target: `/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_duplicates.R`
- classification: `review`
- active stop point: wave 1 is now applied live via
  [tools/refactor/manifest-lipid-qc-duplicates-module-wave1.yml](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/manifest-lipid-qc-duplicates-module-wave1.yml:1)
  into
  [R/mod_lipid_qc_duplicates_server_helpers.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_duplicates_server_helpers.R:1)
  while the public `mod_lipid_qc_duplicates_ui()` /
  `mod_lipid_qc_duplicates_server()` identity stays in
  [R/mod_lipid_qc_duplicates.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_duplicates.R:1).
- next step: none for this wrapper target; treat future helper-level cleanup as
  optional follow-up outside the god-module backlog lane.

## Existing Safety Net

- focused wrapper gate:
  - [tests/testthat/test-lipid-01c-qc-duplicate-module-seams.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tests/testthat/test-lipid-01c-qc-duplicate-module-seams.R:1)
- replay command:
  - `Rscript -e "testthat::test_file('tests/testthat/test-lipid-01c-qc-duplicate-module-seams.R', stop_on_failure = TRUE)"`
- gate coverage now characterizes:
  - `buildLipidDuplicateSummaryUi()` preserving the pre-detection prompt
  - `buildLipidDuplicateSummaryUi()` preserving per-assay duplicate counts and
    status-icon rendering
  - `registerLipidDuplicateTablesOutput()` preserving the duplicate-table
    render registration shell handoff contract
  - `buildLipidDuplicateTablesUi()` returning `NULL` before duplicate detection
  - the no-duplicates fallback panel text
  - the duplicate-tab tabset id, assay label, and sanitized output id contract
  - `registerLipidDuplicateSummaryOutput()` preserving the summary render
    registration shell handoff contract
  - `registerLipidDuplicateFilterPlotOutput()` preserving the QC-progress grob
    render registration shell handoff contract
  - `registerLipidDuplicateFilterPlotOutput()` preserving the QC-progress
    ggplot render registration shell handoff contract
  - `initializeLipidDuplicateServerState()` preserving the three reactive-value
    default shell contract for duplicate info, resolution stats, and filter
    plot state
  - `registerLipidDuplicateServerBindings()` preserving the live server
    registration shell handoff contract
  - `runLipidDuplicateModuleServerShell()` preserving the moduleServer body
    handoff contract for `session$ns`, duplicate state initialization, and live
    binding registration
  - `registerLipidDuplicateTableOutputs()` registering only non-empty duplicate
    assays with sanitized output ids
  - `registerLipidDuplicateTableObserver()` preserving the table-registration
    observer shell handoff contract
  - `handleLipidDuplicateDetection()` preserving the detect-summary workflow
    contract
  - `handleLipidDuplicateDetection()` preserving the zero-duplicate
    notification contract
  - `applyLipidDuplicateDetectionResult()` preserving the post-detection apply
    contract
  - `registerLipidDuplicateDetectObserver()` preserving the detect observer
    shell handoff contract
  - `registerLipidDuplicateDetectObserver()` preserving the detect failure
    notification contract
  - `handleLipidDuplicateResolution()` preserving the resolve-save-summary
    workflow contract
  - `handleLipidDuplicateResolution()` preserving the no-numeric assay fallback
    and QC-plot warning contract
  - `applyLipidDuplicateResolutionResult()` preserving the post-resolution
    reset and notification contract
  - `registerLipidDuplicateResolveObserver()` preserving the resolve observer
    shell handoff contract
  - `registerLipidDuplicateResolveObserver()` preserving the failure cleanup
    contract
  - `applyLipidDuplicateRevertResult()` preserving the post-revert reset and
    notification contract
  - `registerLipidDuplicateRevertObserver()` preserving the revert observer
    shell handoff contract
  - `registerLipidDuplicateRevertObserver()` preserving the revert failure
    notification contract
  - `handleLipidDuplicateRevert()` preserving the revert workflow contract
  - `handleLipidDuplicateRevert()` preserving the no-history failure contract

## Notes

- No prior target handover existed for this wrapper lane.
- April 16, 2026 introduced one bounded live seam in
  [mod_lipid_qc_duplicates.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_duplicates.R:96)
  by extracting the duplicate-table tabset builder into
  `buildLipidDuplicateTablesUi()`.
- The live wrapper now delegates the duplicate-table rendering branch through
  that helper without changing the outer module API.
- April 16, 2026 introduced a second bounded live seam in
  [mod_lipid_qc_duplicates.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_duplicates.R:175)
  by extracting the duplicate-summary UI builder into
  `buildLipidDuplicateSummaryUi()`.
- The live wrapper now delegates the duplicate-summary rendering branch through
  that helper without changing the outer module API.
- April 16, 2026 introduced a third bounded live seam in
  [mod_lipid_qc_duplicates.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_duplicates.R:143)
  by extracting the per-assay duplicate-table registration observer into
  `registerLipidDuplicateTableOutputs()` with shared
  `sanitizeLipidDuplicateTableOutputId()` and `buildLipidDuplicateDatatable()`
  contracts.
- The live wrapper now delegates the duplicate-table `renderDT()` registration
  path through that seam without changing the outer module API.
- April 16, 2026 introduced a fourth bounded live seam in
  [mod_lipid_qc_duplicates.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_duplicates.R:241)
  by extracting the duplicate-resolution workflow into
  `handleLipidDuplicateResolution()`.
- The live wrapper now delegates per-assay duplicate resolution, state save,
  QC-plot refresh, and result-summary text generation through that helper
  without changing the outer module API.
- April 16, 2026 introduced a fifth bounded live seam in
  [mod_lipid_qc_duplicates.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_duplicates.R:200)
  by extracting the duplicate-detection workflow into
  `handleLipidDuplicateDetection()`.
- The live wrapper now delegates duplicate detection, duplicate counting,
  notification payload selection, and detection logging through that helper
  without changing the outer module API.
- April 16, 2026 introduced a sixth bounded live seam in
  [mod_lipid_qc_duplicates.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_duplicates.R:351)
  by extracting the duplicate-revert workflow into
  `handleLipidDuplicateRevert()`.
- The live wrapper now delegates state-history lookup, revert dispatch, revert
  result text generation, and success notification payload selection through
  that helper without changing the outer module API.
- April 16, 2026 introduced a seventh bounded live seam in
  [mod_lipid_qc_duplicates.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_duplicates.R:351)
  by extracting the resolve observer's post-resolution reset and notification
  workflow into `applyLipidDuplicateResolutionResult()`.
- The live wrapper now delegates post-resolution reactive state updates, result
  text rendering, duplicate-table reset, success logging, and working/success
  notification dispatch through that helper without changing the outer module
  API.
- April 16, 2026 introduced an eighth bounded live seam in
  [mod_lipid_qc_duplicates.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_duplicates.R:388)
  by extracting the revert observer's post-revert reset and notification
  workflow into `applyLipidDuplicateRevertResult()`.
- The live wrapper now delegates post-revert result rendering, duplicate-state
  reset, QC-plot reset, and success-notification dispatch through that helper
  without changing the outer module API.
- April 16, 2026 introduced a ninth bounded live seam in
  [mod_lipid_qc_duplicates.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_duplicates.R:351)
  by extracting the detect observer's post-detection reactive-state update and
  notification workflow into `applyLipidDuplicateDetectionResult()`.
- The live wrapper now delegates duplicate-info state updates and detection
  notification dispatch through that helper without changing the outer module
  API.
- April 16, 2026 introduced a tenth bounded live seam in
  [mod_lipid_qc_duplicates.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_duplicates.R:406)
  by extracting the resolve observer shell into
  `registerLipidDuplicateResolveObserver()`.
- The live wrapper now delegates working-notification setup, resolution/apply
  handoff, and failure cleanup through that helper without changing the outer
  module API.
- April 16, 2026 introduced an eleventh bounded live seam in
  [mod_lipid_qc_duplicates.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_duplicates.R:369)
  by extracting the detect observer shell into
  `registerLipidDuplicateDetectObserver()`.
- The live wrapper now delegates `detect_duplicates` req gating, detect/apply
  handoff, and error-notification dispatch through that helper without
  changing the outer module API.
- April 16, 2026 introduced a twelfth bounded live seam in
  [mod_lipid_qc_duplicates.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_duplicates.R:517)
  by extracting the revert observer shell into
  `registerLipidDuplicateRevertObserver()`.
- The live wrapper now delegates `revert_duplicates` revert/apply handoff and
  error-notification dispatch through that helper without changing the outer
  module API.
- April 16, 2026 introduced a thirteenth bounded live seam in
  [mod_lipid_qc_duplicates.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_duplicates.R:175)
  by extracting the duplicate-table registration observer shell into
  `registerLipidDuplicateTableObserver()`.
- The live wrapper now delegates req-gated `duplicate_info()` handoff and
  duplicate-table registration through that helper without changing the outer
  module API.
- April 16, 2026 introduced a fourteenth bounded live seam in
  [mod_lipid_qc_duplicates.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_duplicates.R:220)
  by extracting the duplicate-summary render registration shell into
  `registerLipidDuplicateSummaryOutput()`.
- The live wrapper now delegates `duplicate_info()` handoff and
  `output$duplicate_summary` registration through that helper without changing
  the outer module API.
- April 16, 2026 introduced a fifteenth bounded live seam in
  [mod_lipid_qc_duplicates.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_duplicates.R:233)
  by extracting the duplicate-tables render registration shell into
  `registerLipidDuplicateTablesOutput()`.
- The live wrapper now delegates `duplicate_info()` and `ns` handoff plus
  `output$duplicate_tables` registration through that helper without changing
  the outer module API.
- April 16, 2026 introduced a sixteenth bounded live seam in
  [mod_lipid_qc_duplicates.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_duplicates.R:250)
  by extracting the QC-progress render registration shell into
  `registerLipidDuplicateFilterPlotOutput()`.
- The live wrapper now delegates `filter_plot()` handoff plus grob/ggplot
  dispatch through that helper without changing the outer module API.
- April 16, 2026 introduced a seventeenth bounded live seam in
  [mod_lipid_qc_duplicates.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_duplicates.R:661)
  by extracting the live server registration branch into
  `registerLipidDuplicateServerBindings()`.
- The live wrapper now delegates detect, summary, tables, per-assay table
  observer, resolve, revert, and QC-progress output binding through that
  helper without changing the outer module API.
- April 16, 2026 introduced an eighteenth bounded live seam in
  [mod_lipid_qc_duplicates.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_duplicates.R:96)
  by extracting the `moduleServer()` reactive-value shell setup into
  `initializeLipidDuplicateServerState()`.
- The live wrapper now delegates duplicate-info, resolution-stats, and
  QC-progress reactive-value initialization through that helper without
  changing the outer module API.
- April 16, 2026 introduced a nineteenth bounded live seam in
  [mod_lipid_qc_duplicates.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_duplicates.R:727)
  by extracting the `moduleServer()` body shell into
  `runLipidDuplicateModuleServerShell()`.
- The live wrapper now delegates `session$ns` capture, duplicate-state
  initialization, and live server-binding registration through that helper
  without changing the outer module API.
- The focused duplicate-wrapper gate reran green after the live seam.
- Post-seam classification metrics are `769` lines, `23` top-level functions,
  and max top-level function length `58`.
- The classifier now reports `review`, so this target remained
  `in_progress` at a review/staging stop point rather than a further live-seam
  stop point once the moduleServer body was centralized.
- April 16, 2026 completed one bounded staged-wave checkpoint by verifying and
  staging
  [tools/refactor/manifest-lipid-qc-duplicates-module-wave1.yml](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/manifest-lipid-qc-duplicates-module-wave1.yml:1)
  for the seam-ready duplicate-module helper surface in
  [R/mod_lipid_qc_duplicates.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_duplicates.R:1),
  materializing
  [R/mod_lipid_qc_duplicates_server_helpers.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/staging/wave1_lipidomics_qc_duplicates_module_server_helpers/R/mod_lipid_qc_duplicates_server_helpers.R:1)
  and
  [tools/refactor/collate-lipid-qc-duplicates-module-wave1.txt](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/collate-lipid-qc-duplicates-module-wave1.txt:1).
- The staged helper wave lifts `21` seam-ready helper and observer functions
  into a dedicated staged helper file, while the live wrapper remains
  unchanged and the public `mod_lipid_qc_duplicates_ui()` /
  `mod_lipid_qc_duplicates_server()` identity stays in
  [R/mod_lipid_qc_duplicates.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_duplicates.R:1).
- The focused duplicate-wrapper gate reran green after the staged-wave
  checkpoint with `203` source-based passes.
- Live post-staging classification for
  [R/mod_lipid_qc_duplicates.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_duplicates.R:1)
  remains `769` lines, `23` top-level functions, max top-level function length
  `58`, and label `review`, while the staged helper target lands at `657`
  lines, `21` top-level functions, max top-level function length `31`, and
  label `direct-extraction-ready`.
- April 16, 2026 completed one bounded live-apply checkpoint by applying
  [tools/refactor/manifest-lipid-qc-duplicates-module-wave1.yml](/home/doktersmol/Documents/MultiScholaR-lipid-lane/tools/refactor/manifest-lipid-qc-duplicates-module-wave1.yml:1)
  into
  [R/mod_lipid_qc_duplicates_server_helpers.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_duplicates_server_helpers.R:1)
  and the trimmed public wrapper
  [R/mod_lipid_qc_duplicates.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_duplicates.R:1).
- The live apply moved `21` seam-ready helper and observer functions out of the
  public wrapper, refreshed the live collate order in `DESCRIPTION`, and kept
  the public module entrypoints stable in `mod_lipid_qc_duplicates.R`.
- The focused duplicate-wrapper gate reran green after the live apply once the
  characterization test sourced both the helper file and wrapper file, again
  landing at `203` source-based passes.
- Live post-apply classification for
  [R/mod_lipid_qc_duplicates.R](/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc_duplicates.R:1)
  is `133` lines, `2` top-level functions, max top-level function length `58`,
  and label `review`.
- Lipid-QC duplicates wrapper stabilization is now complete for this backlog
  target; no further live seam or staged wave is required in
  `mod_lipid_qc_duplicates.R`.
