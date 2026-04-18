# Peptide Replicate QC Seam Map

## Goal

Introduce one safe in-file seam at a time in R/mod_prot_qc_peptide_replicate.R
while preserving the module's public contract.

## Current Position In The Flow

- target: `/home/doktersmol/Documents/MultiScholaR/R/mod_prot_qc_peptide_replicate.R`
- classification: `review`
- checkpoint reached: `revert-observer helper seam landed live in R/mod_prot_qc_peptide_replicate.R`
- next step: `Manual target complete; keep the focused peptide-QC gate as the regression surface and do not reopen this wrapper unless a real regression appears.`

## Existing Safety Net

- `tests/testthat/test-prot-02-qc-filtering.R`
- `tests/testthat/test-prot-02-qc-peptide-groupaware.R`
- `tests/testthat/test-prot-02b-qc-peptide-module-contracts.R`
- `tests/testthat/test-prot-03-rollup.R`

## Notes

- This manual target previously had no dedicated handover; this file now marks
  the first explicit stop point for the peptide replicate QC wrapper.
- April 16, 2026 stabilize iteration landed one bounded live revert-observer
  seam in `R/mod_prot_qc_peptide_replicate.R` via
  `runPeptideReplicateRevertStep()` and
  `runPeptideReplicateRevertObserver()`.
- The new helper pair now owns history lookup, state-manager revert
  delegation, revert result text assembly, success notification, and revert
  error logging/notification while the wrapper keeps the public module id,
  apply observer wiring, and plot render binding unchanged.
- The same checkpoint rewired the inline `revert_replicate` observer through
  the new helper seam without changing the public revert button id, result
  output id, or plot output id.
- The same checkpoint added direct characterization coverage in
  `tests/testthat/test-prot-02b-qc-peptide-module-contracts.R` for the new
  revert step and observer helper happy/error paths plus wrapper-level revert
  delegation.
- Focused gate reran green after the live seam for:
  - `tests/testthat/test-prot-02-qc-filtering.R`
  - `tests/testthat/test-prot-02-qc-peptide-groupaware.R`
  - `tests/testthat/test-prot-02b-qc-peptide-module-contracts.R`
  - `tests/testthat/test-prot-03-rollup.R`
- April 16, 2026 stabilize iteration landed one bounded live apply-path seam
  in `R/mod_prot_qc_peptide_replicate.R` via:
  `runPeptideReplicateApplyStep()`,
  `updatePeptideReplicateOutputs()`,
  and
  `runPeptideReplicateApplyObserver()`.
- The new helper cluster now owns the S4 filter call, QC parameter tracking,
  state-save metadata, result text assembly, plot refresh, and apply
  success/error notification handling while the wrapper keeps the public module
  id, revert observer wiring, and plot render binding unchanged.
- The same checkpoint rewired the inline `apply_replicate_filter` observer
  through the new helper seam without changing the public apply button id,
  result output id, or plot output id.
- The same checkpoint added direct characterization coverage in
  `tests/testthat/test-prot-02b-qc-peptide-module-contracts.R`
  for the new apply step, output-refresh helper, apply observer happy/error
  paths, and wrapper-level apply delegation.
- Focused gate reran green after the live seam for:
  - `tests/testthat/test-prot-02-qc-filtering.R`
  - `tests/testthat/test-prot-02-qc-peptide-groupaware.R`
  - `tests/testthat/test-prot-02b-qc-peptide-module-contracts.R`
  - `tests/testthat/test-prot-03-rollup.R`
- `R/mod_prot_qc_peptide_replicate.R` now measures `282` live lines with `7`
  top-level functions and currently classifies as `review`.
- `R/mod_prot_qc_peptide_replicate.R` is now a stabilized thin
  wrapper/manual target complete; keep the focused peptide-QC gate as the
  regression surface and do not reopen this file for more stabilization work
  unless a real regression appears.
