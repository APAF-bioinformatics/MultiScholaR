# Protein QC Orchestrator Seam Map

## Goal

Freeze and stabilize the public wrapper in `R/mod_prot_qc_protein.R` after the
earlier helper/file extractions landed elsewhere in the QC and rollup bucket.

## Current Position In The Flow

- target: `/home/doktersmol/Documents/MultiScholaR/R/mod_prot_qc_protein.R`
- classification: `review`
- checkpoint reached: `wrapper characterization plus one live module-spec seam landed in R/mod_prot_qc_protein.R`
- next step: `Archive this wrapper target unless a future behavior change reopens the orchestrator contract.`

## Existing Safety Net

- `tests/testthat/test-prot-02-qc-filtering.R`
- `tests/testthat/test-prot-02-qc-peptide-groupaware.R`
- `tests/testthat/test-prot-02c-qc-protein-module-contracts.R`
- `tests/testthat/test-prot-03-rollup.R`

## Notes

- Manual target bucket 0 follow-up for the protein QC public wrapper after the
  earlier helper/module extractions moved the heavy logic into dedicated files.
- April 16, 2026 stabilize iteration added direct wrapper characterization in
  `tests/testthat/test-prot-02c-qc-protein-module-contracts.R` to freeze:
  `mod_prot_qc_protein_ui()` workflow-dependent tab labels and namespaced ids,
  plus `mod_prot_qc_protein_server()` DIA rollup dispatch and common-module
  fan-out.
- The same checkpoint landed one bounded live module-spec seam in
  `R/mod_prot_qc_protein.R` via:
  `getProtQcProteinRollupModuleSpec()`,
  `getProtQcProteinCommonModuleSpecs()`,
  `getProtQcProteinModuleSpecs()`,
  `buildProtQcProteinTab()`,
  and
  `runProtQcProteinSubmodule()`,
  collapsing the repeated inline existence checks while preserving the public
  tab ids, fallback labels, and submodule call order.
- Focused gate reran green for:
  - `test-prot-02-qc-filtering`
  - `test-prot-02-qc-peptide-groupaware`
  - `test-prot-02c-qc-protein-module-contracts`
  - `test-prot-03-rollup`
- `R/mod_prot_qc_protein.R` now measures `143` lines and is a stabilized
  orchestrator wrapper; no additional extraction/staging work is needed for
  this manual target.
