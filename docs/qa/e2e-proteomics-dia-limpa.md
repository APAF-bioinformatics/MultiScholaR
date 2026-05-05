# E2E-006 Proteomics DIA limpa Imputation Lane

E2E-006 covers the DIA branch that performs protein rollup through limpa
DPC-Quant instead of the canonical IQ MaxLFQ route. The purpose is not to
numerically revalidate limpa itself; it is to prove that the GUI wiring can
select the limpa path and preserve that decision through state, export, DA
reload, report-parameter serialization, and branch-specific report routing.

## Browser Contract

The tracked scenario lives in
`tests/testthat/test-e2e-proteomics-dia-limpa.R` and drives the production app
through these handoffs:

- launches `run_app(test_mode = TRUE)` through the shared browser harness;
- selects proteomics and creates an isolated temporary project;
- uploads `prot_dia_limpa/report.tsv` and the lane FASTA through the import UI;
- completes the design builder from `prot_dia_limpa/design.tsv`;
- runs the full DIA peptide QC sequence;
- selects `limpa DPC-Quant` in the protein rollup module before applying rollup;
- continues through accession cleanup, protein filtering, deduplication, and
  protein replicate filtering;
- runs deterministic normalization with RUV skipped;
- skips final correlation filtering and exports the filtered session;
- reloads that filtered session in DA and runs DA;
- saves workflow arguments, copies publication artifacts, renders the report,
  downloads the report, and exports summary session state.

## Required Evidence

The test fails unless the limpa branch is observable in all of the following
places:

- `filtered_session_data_latest.rds` has `workflow_type == "DIA"`,
  `report_template == "DIANN_limpa_report.rmd"`, `limpa_applied == TRUE`, and
  `use_limpa == TRUE`;
- the current exported S4 object has `@args$globalParameters$use_limpa == TRUE`
  and `@args$globalParameters$report_template == "DIANN_limpa_report.rmd"`;
- the current exported S4 object contains
  `@args$proteinMissingValueImputationLimpa`;
- the current exported S4 object contains
  `@args$limpa_dpc_quant_results$dpc_method == "limpa_dpc_quant"`;
- the DPC quantification method starts with `limpa_dpc_quant`;
- `limpa_parameters.RDS` and `limpa_dpc_quant_results.RDS` are written beside
  the filtered session export;
- `filtered_session_summary.txt` records `Report Template:
  DIANN_limpa_report.rmd` and `limpa DPC-Quant: yes`;
- `study_parameters.txt` records `DIANN_limpa_report.rmd`,
  `proteinMissingValueImputationLimpa`, and `limpa_dpc_quant`;
- the source directory contains `DIANN_limpa_report.rmd` and does not create
  `DIANN_report.rmd` for this lane;
- rendered report outputs match
  `DIANN_limpa_report_proteomics_.*\\.(docx|html|pdf)`;
- `session_state_*.RDS` is exported from the summary module;
- browser console errors and Shiny error notifications are absent.

## Test-Mode Behavior

The CI/browser lane does not require the optional `limpa` package to be installed.
When `limpa` is available, the rollup module calls
`proteinMissingValueImputationLimpa()`. In test mode without the package, it uses
a deterministic fallback that builds a protein-level S4 object from peptide
quantities while preserving the same branch metadata:

- `globalParameters$use_limpa`;
- `globalParameters$report_template`;
- `proteinMissingValueImputationLimpa`;
- `limpa_dpc_quant_results`;
- `quantification_method`;
- DPC peptide/protein count summaries.

This keeps the GUI wiring test deterministic while still requiring the same
downstream DA/report contracts as the real limpa path.

## Regression Surface

E2E-006 specifically guards regressions in these touch points:

- the protein QC tab must expose a user-selectable limpa rollup option;
- the selected rollup branch must be passed through the server observer seam;
- the state manager must save `protein_s4_created` from the limpa branch;
- exported normalization sessions must include explicit limpa fields rather than
  relying on a human to inspect nested S4 args;
- DA reload must tolerate a filtered session containing DPC-Quant metadata;
- DA wrapper DPC-Quant detection must not silently miss the branch because of
  stale variable names;
- summary report resolution must choose `DIANN_limpa_report.rmd` when the branch
  is selected;
- fixture manifest metadata must point the `prot_dia_limpa` lane at the limpa
  report template.

## Focused Runner

Run the ticket suite with:

```r
Sys.setenv(NOT_CRAN = "true")
pkgload::load_all(export_all = TRUE, helpers = FALSE, attach_testthat = FALSE)
testthat::test_dir(
  "tests/testthat",
  filter = "^(e2e-manifest|e2e-fixture-integrity|e2e-browser-harness|prot-02c-qc-protein-module-contracts|prot-05b-norm-module-contracts|prot-12-summary-module-contracts|e2e-proteomics-dia-limpa)$",
  reporter = "summary"
)
```
