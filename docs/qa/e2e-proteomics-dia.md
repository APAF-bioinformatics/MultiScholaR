# E2E-004 Proteomics DIA Canonical Workflow

This document records the executable browser contract for the canonical
proteomics DIA lane. The tracked scenario lives in
`tests/testthat/test-e2e-proteomics-dia.R`.

## Scenario

The test drives the real Shiny GUI through the production module wiring:

- launch `run_app(test_mode = TRUE)`
- select proteomics
- create a project in a temporary output directory
- upload the DIA-NN `prot_dia/report.tsv` fixture
- force the import format to DIA-NN
- complete the design builder from `prot_dia/design.tsv`
- run every DIA peptide QC stage
- run every DIA protein QC stage
- run normalization with RUV skipped for deterministic CI runtime
- skip final correlation filtering and export the filtered session
- load the exported filtered session in DA
- run DA
- save workflow arguments
- copy publication artifacts
- seed and render the DIA report template
- download the generated report
- export the summary session state

## Required Invariants

The browser test asserts the following production handoffs:

- the state digest resolves workflow type as `DIA`
- proteomics project directories exist under the generated project root
- import, design, QC, normalization, and DA step statuses advance through the live
  `tab_status` keys
- `filtered_session_data_latest.rds` and `filtered_session_summary.txt` are
  created in `scripts/proteomics`
- the exported filtered session contains the full R6 state payload required by
  the DA loader: `r6_current_state_name`, `r6_complete_states`, and
  `r6_state_history`
- the DA completion status updates `differential_expression`, the key observed by
  the proteomics workflow shell
- `DIANN_report.rmd` is present in `scripts/proteomics`
- `study_parameters.txt`, the rendered report, and the downloaded report are
  non-empty
- browser console errors and error notifications are absent at closeout

## Determinism

Live report-template downloads are avoided by copying the tracked E2E stub from
`tests/testdata/e2e/report_templates/DIANN_report.rmd` into the generated project
before report rendering. Enrichment is intentionally not executed in this ticket;
the backend matrix is covered separately by E2E-007.

The scenario skips only when local `shinytest2`, `chromote`, or Chrome are not
available. Harness-level unit tests still validate the action dispatch and artifact
assertions without requiring a browser.

## Closeout Notes

E2E-004 now guards the browser-visible regressions found while wiring the DIA lane:

- import completion must unlock design without skipping the design builder
- DIA protein replicate filtering must execute once, complete QC status, and avoid
  parallel-cluster setup in test mode
- extracted protein S4 QC and normalization methods must remain registered after
  package collation
- workflow status observers must be monotonic so completed QC/normalization states
  are not downgraded by older observers
- normalization export must preserve the R6 current state, complete state map, and
  state history for DA reload
- limma normalization functions must resolve in the app process while still allowing
  method-level test doubles
- PCA QC fallback must tolerate degenerate small fixtures with constant features
- browser console assertions must distinguish real JavaScript failures from benign
  Shiny/DataTables lifecycle logs

Verified on 2026-05-05 with:

```r
Sys.setenv(NOT_CRAN = "true")
pkgload::load_all(export_all = TRUE, helpers = FALSE, attach_testthat = FALSE)
testthat::test_dir(
  "tests/testthat",
  filter = "^(app-shell-shared|e2e-browser-harness|e2e-manifest|e2e-test-mode-integration|e2e-proteomics-dia|general-pca-rle-direct-shared|general-plotting-pca-helper-characterization|prot-01c-import-module-contracts|prot-01d-import-module-characterization|prot-02c-qc-protein-module-contracts|prot-qc-coordinator-direct-shared|prot-05b-norm-module-contracts|prot-05d-norm-module-characterization|prot-s4-norm-methods-direct-shared|prot-07b-da-handlers-compat|prot-s4-qc-direct-shared)$",
  reporter = "summary"
)
```
