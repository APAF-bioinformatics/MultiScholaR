# Peptide Rollup Module Seam Map

## Goal

Freeze and stabilize the public wrapper in `R/mod_prot_qc_peptide_rollup.R`
with bounded in-file seams before any larger wrapper decomposition.

## Current Position In The Flow

- target: `/home/doktersmol/Documents/MultiScholaR/R/mod_prot_qc_peptide_rollup.R`
- classification: `review`
- checkpoint reached: `wrapper revert-path seam landed in R/mod_prot_qc_peptide_rollup.R`
- next step: `Wrapper stabilization stop point is now complete; keep the focused peptide-QC gate as the regression surface and do not reopen this file unless a real regression appears.`

## Existing Safety Net

- `tests/testthat/test-prot-02-qc-filtering.R`
- `tests/testthat/test-prot-02-qc-peptide-groupaware.R`
- `tests/testthat/test-prot-02b-qc-peptide-module-contracts.R`
- `tests/testthat/test-prot-03-rollup.R`

## Notes

- Manual target bucket 0 follow-up for the peptide rollup wrapper after the
  adjacent peptide QC wrappers were stabilized.
- April 16, 2026 stabilize iteration landed one bounded live apply-path seam
  in `R/mod_prot_qc_peptide_rollup.R` via:
  - `runPeptideRollupApplyStep()`
  - `updatePeptideRollupOutputs()`
  - `runPeptideRollupApplyObserver()`
- The new helper cluster now owns the rollup S4 transformation call, QC
  parameter tracking, state-save metadata, result text assembly, plot refresh,
  and apply success/error notification handling while the wrapper keeps the
  public module id, inline revert observer wiring, and plot render binding
  unchanged.
- The same checkpoint rewired the inline `apply_rollup` observer through the
  new top-level helper seam without changing the public apply button id,
  result output id, or plot output id.
- The same checkpoint added direct characterization coverage in
  `tests/testthat/test-prot-02b-qc-peptide-module-contracts.R` for the new
  apply step, output-refresh helper, apply observer happy/error paths, and
  wrapper-level apply delegation.
- April 16, 2026 stabilize iteration then landed one bounded live revert-path
  seam in `R/mod_prot_qc_peptide_rollup.R` via:
  - `runPeptideRollupRevertStep()`
  - `runPeptideRollupRevertObserver()`
- The same checkpoint moved the previous inline `revert_rollup` history
  handling, result text rendering, success notification, and error reporting
  behind a top-level helper seam without changing the public revert button id
  or result output id.
- The same checkpoint extended
  `tests/testthat/test-prot-02b-qc-peptide-module-contracts.R` with direct
  revert-step characterization, revert observer happy/error path coverage, and
  wrapper-level revert delegation coverage for
  `mod_prot_qc_peptide_rollup_server()`.
- Focused gate reran green again for:
  - `test-prot-02-qc-filtering`
  - `test-prot-02-qc-peptide-groupaware`
  - `test-prot-02b-qc-peptide-module-contracts`
  - `test-prot-03-rollup`
- Focused gate reran green for:
  - `test-prot-02-qc-filtering`
  - `test-prot-02-qc-peptide-groupaware`
  - `test-prot-02b-qc-peptide-module-contracts`
  - `test-prot-03-rollup`
- `R/mod_prot_qc_peptide_rollup.R` now measures `260` live lines with `7`
  top-level functions and still heuristically classifies as `review`, but the
  wrapper is now a stabilized thin wrapper/manual target complete.
