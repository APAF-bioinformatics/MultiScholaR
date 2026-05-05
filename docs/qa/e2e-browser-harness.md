# GUI E2E Browser Harness

This document defines the tracked E2E-003 browser automation contract for
MultiScholaR's Shiny GUI. The goal is to make every workflow test read as a
scenario, not as a pile of raw browser-driver calls.

## Scope

The shared harness in `tests/testthat/helper-e2e-browser.R` owns the reusable
browser mechanics for all current production omics:

- proteomics
- metabolomics
- lipidomics

Transcriptomics and integration tiles are intentionally excluded from the active
workflow matrix until backing production modules exist.

## Launch Contract

The harness supports both launch surfaces required by the E2E plan:

- `launch_surface = "run_app"` starts `MultiScholaR::run_app(test_mode = TRUE)`.
- `launch_surface = "app_dir"` starts `inst/app/app.R` with
  `multischolar.test_mode` and `shiny.testmode` enabled in the child process.

Every browser launch uses deterministic defaults:

- fixed random seed
- explicit load timeout and interaction timeout
- fixed browser viewport
- `MULTISCHOLAR_E2E_ARTIFACT_DIR` override for CI artifacts
- stable snapshot variant based on OS and R version

## Selector Contract

Scenario tests should use semantic helpers instead of raw CSS wherever possible.

- `e2e_selector()` maps a `data-testid` to a CSS selector.
- `e2e_click_testid()` clicks semantic controls and waits for idle.
- `e2e_input_id()` maps omic, step, and input names to namespaced Shiny input IDs.
- `e2e_switch_workflow_tab()` switches workflow tabs by stable Shiny input IDs.

The current namespace mapping is:

- proteomics import: `proteomics_workflow-setup_import-*`
- proteomics design: `proteomics_workflow-design_matrix-*`
- proteomics normalization: `proteomics_workflow-normalization-*`
- proteomics DA: `proteomics_workflow-differential_expression-*`
- proteomics summary: `proteomics_workflow-session_summary-*`
- metabolomics import: `metabolomics_workflow-import-*`
- metabolomics design: `metabolomics_workflow-design-*`
- metabolomics normalization: `metabolomics_workflow-norm-*`
- metabolomics DA: `metabolomics_workflow-de-*`
- metabolomics summary: `metabolomics_workflow-summary-*`
- lipidomics import: `lipidomics_workflow-import-*`
- lipidomics design: `lipidomics_workflow-design-*`
- lipidomics normalization: `lipidomics_workflow-norm-*`
- lipidomics DA: `lipidomics_workflow-de-*`
- lipidomics summary: `lipidomics_workflow-summary-*`

## Import Contract

Test mode must not depend on native OS file pickers. Raw import modules therefore
need standard browser file inputs and matching server observers.

- Proteomics uses `search_results_standard` and `fasta_file_standard`.
- Metabolomics uses `assay1_file_std` and `assay2_file_std`.
- Lipidomics uses `assay1_file_std` and `assay2_file_std`.

The harness uses `e2e_upload_lane_inputs()` to map fixture-manifest lanes onto
the correct current input IDs. Full workflow tickets should not duplicate this
mapping.

## Workflow Scenario Shape

One browser test should map to one readable workflow scenario:

1. Launch the app in E2E mode.
2. Select one or more omics with `e2e_select_omics()`.
3. Complete the setup modal with `e2e_complete_project_setup()`.
4. Upload fixture inputs with `e2e_upload_lane_inputs()`.
5. Drive import, design, normalization, DA, enrichment, summary, and report
   actions with the semantic action helpers.
6. Assert the state digest with `e2e_assert_state_digest()`.
7. Assert browser cleanliness with `e2e_assert_no_console_errors()` and
   `e2e_assert_no_error_notifications()`.
8. Validate downloads and generated artifacts with `e2e_get_download()` and
   `e2e_assert_file_nonempty()`.

Scenario tests should prefer explicit parameter assertions over snapshots for
data integrity. Screenshots and DOM snapshots are for forensics, not the primary
oracle for scientific correctness.

## Failure Artifacts

`e2e_capture_failure_artifacts()` writes a case-specific artifact directory with:

- `failure-reason.txt`
- `dom.html`
- `browser-logs.json`
- `screenshot.png`

CI should upload `MULTISCHOLAR_E2E_ARTIFACT_DIR` on failure for every smoke,
nightly, and release-gate browser job.

## State Digest Assertions

The app-level `test_state_digest` is the cross-workflow truth source for GUI
state. Browser tests should assert at least:

- selected omics
- initialized omics
- required project directory keys
- workflow type per omic
- active tab or step status after major transitions
- export paths after normalization-session export
- report fingerprints after summary/report generation

This avoids fragile DOM scraping and directly detects wiring regressions from
the de-monolithing work.

## E2E-003 Verification

The tracked unit tests in `tests/testthat/test-e2e-browser-harness.R` verify:

- selector construction
- omic namespace mapping
- launch-surface mapping
- driver action dispatch
- upload helper dispatch
- state digest parsing and invariant assertions
- console error detection
- failure artifact capture
- an optional real `shinytest2` smoke path when browser dependencies are present

The tracked integration tests in `tests/testthat/test-e2e-test-mode-integration.R`
also verify that metabolomics and lipidomics expose browser-upload controls in
test mode and that their server-side standard upload handlers preserve selected
datapaths.
