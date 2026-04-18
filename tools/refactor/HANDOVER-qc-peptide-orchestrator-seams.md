# Peptide QC Orchestrator Seam Map

## Goal

Freeze and stabilize the public wrapper in R/mod_prot_qc_peptide.R after the earlier helper/file extractions landed.

## Current Position In The Flow

- target: `/home/doktersmol/Documents/MultiScholaR/R/mod_prot_qc_peptide.R`
- classification: `review`
- checkpoint reached: `wrapper characterization plus one live module-spec seam landed in R/mod_prot_qc_peptide.R`
- next step: `Archive this wrapper target unless a future behavior change reopens the orchestrator contract.`

## Existing Safety Net

- `tests/testthat/test-prot-02-qc-filtering.R`
- `tests/testthat/test-prot-02-qc-peptide-groupaware.R`
- `tests/testthat/test-prot-02b-qc-peptide-module-contracts.R`
- `tests/testthat/test-prot-03-rollup.R`

## Notes

- Manual target bucket 0 follow-up for the peptide QC public wrapper after the
  earlier helper/module extractions moved the heavy logic into dedicated files.
- April 16, 2026 stabilize iteration added direct wrapper characterization in
  `tests/testthat/test-prot-02b-qc-peptide-module-contracts.R` to freeze:
  `mod_prot_qc_peptide_ui()` tab labels and namespaced ids, plus
  `mod_prot_qc_peptide_server()` DIA-only submodule dispatch.
- The same checkpoint landed one bounded live seam in
  `R/mod_prot_qc_peptide.R` via:
  `getProtQcPeptideModuleSpecs()`,
  `buildProtQcPeptideTab()`,
  and
  `runProtQcPeptideSubmodule()`,
  collapsing the repeated wrapper fan-out while preserving the existing public
  tab ids, fallback labels, and submodule call order.
- Focused gate reran green for:
  - `test-prot-02-qc-filtering`
  - `test-prot-02-qc-peptide-groupaware`
  - `test-prot-02b-qc-peptide-module-contracts`
  - `test-prot-03-rollup`
- `R/mod_prot_qc_peptide.R` now measures `128` lines and is a stabilized
  orchestrator wrapper; no additional extraction/staging work is needed for
  this manual target.
