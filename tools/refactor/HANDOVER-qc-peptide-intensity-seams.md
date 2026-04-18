# Peptide Intensity Module Seam Map

## Goal

Freeze and stabilize the public wrapper in R/mod_prot_qc_peptide_intensity.R with bounded in-file seams before any larger observer extraction.

## Current Position In The Flow

- target: `/home/doktersmol/Documents/MultiScholaR/R/mod_prot_qc_peptide_intensity.R`
- classification: `review`
- checkpoint reached: `wrapper bootstrap seam landed in R/mod_prot_qc_peptide_intensity.R`
- next step: `Manual target complete; keep the focused peptide-QC gate as the regression surface and do not reopen this file unless a real regression appears.`

## Existing Safety Net

- `tests/testthat/test-prot-02-qc-filtering.R`
- `tests/testthat/test-prot-02-qc-peptide-groupaware.R`
- `tests/testthat/test-prot-02b-qc-peptide-module-contracts.R`
- `tests/testthat/test-prot-03-rollup.R`

## Notes

- Manual target bucket 0 follow-up for the peptide intensity wrapper after the
  adjacent peptide QC orchestrator and imputation wrappers were stabilized.
- April 16, 2026 stabilize iteration landed one bounded live preview seam in
  `R/mod_prot_qc_peptide_intensity.R` via
  `buildPeptideIntensityThresholdPreview()`.
- The new helper now owns the state-manager lookup, temporary
  `updateMissingValueParameters()` call, and both formatted cutoff preview
  strings while the wrapper keeps the public module id, apply/revert observer
  wiring, and plot binding unchanged.
- The same checkpoint rewired `calculated_groupwise_percent` and
  `calculated_max_groups_percent` through one shared reactive preview so the
  temporary threshold calculation is no longer duplicated across two inline
  renderers.
- The same checkpoint added direct characterization coverage in
  `tests/testthat/test-prot-02b-qc-peptide-module-contracts.R` for the helper's
  state-manager lookup, delegated temporary update call, and formatted cutoff
  strings.
- Focused gate reran green for:
  - `test-prot-02b-qc-peptide-module-contracts.R`
- April 16, 2026 stabilize iteration landed one bounded revert-path seam in
  `R/mod_prot_qc_peptide_intensity.R` via
  `runPeptideIntensityRevertStep()` and
  `runPeptideIntensityRevertObserver()`.
- The new helper pair now owns history lookup, state-manager revert delegation,
  revert result text assembly, success notification, and revert error
  logging/notification while the wrapper keeps the public module id, apply
  observer wiring, and plot binding unchanged.
- The same checkpoint rewired the inline `revert_intensity` observer to the new
  top-level helper shell without changing the public revert button contract.
- The same checkpoint added direct characterization coverage in
  `tests/testthat/test-prot-02b-qc-peptide-module-contracts.R` for the new
  revert step and observer helper happy/error paths.
- Focused gate reran green for:
  - `test-prot-02b-qc-peptide-module-contracts.R`
- April 16, 2026 stabilize iteration landed one bounded apply-path seam in
  `R/mod_prot_qc_peptide_intensity.R` via
  `runPeptideIntensityApplyStep()`,
  `updatePeptideIntensityOutputs()`, and
  `runPeptideIntensityApplyObserver()`.
- The new helper cluster now owns strict/flexible parameter updates, QC
  parameter tracking, state save metadata, result text assembly, plot refresh,
  and apply success/error notification handling while the wrapper keeps the
  public module id, preview reactive wiring, revert observer wiring, and plot
  binding unchanged.
- The same checkpoint rewired the inline `apply_intensity_filter` observer to
  the new top-level helper shell without changing the public apply button
  contract.
- The same checkpoint added direct characterization coverage in
  `tests/testthat/test-prot-02b-qc-peptide-module-contracts.R` for the new
  apply step, output refresh helper, and apply observer happy/error paths.
- Focused gate reran green for:
  - `test-prot-02b-qc-peptide-module-contracts.R`
- April 16, 2026 stabilize iteration landed one bounded wrapper-bootstrap seam
  in `R/mod_prot_qc_peptide_intensity.R` via
  `setupPeptideIntensityServerBootstrap()`.
- The new helper now owns the shared preview reactive, preview output render
  bindings, apply/revert observer registration, and plot render binding while
  the public module wrapper keeps only the `reactiveVal()` bootstrap and module
  entrypoint identity.
- The same checkpoint rewired the remaining inline server wiring through one
  top-level helper shell without changing the public module id, preview text
  contract, apply/revert button contract, or plot output contract.
- The same checkpoint added direct characterization coverage in
  `tests/testthat/test-prot-02b-qc-peptide-module-contracts.R` for the new
  bootstrap helper's preview/output/observer/plot delegation plus the thin
  wrapper's forwarded bootstrap arguments.
- Focused gate reran green for:
  - `test-prot-02-qc-filtering.R`
  - `test-prot-02-qc-peptide-groupaware.R`
  - `test-prot-02b-qc-peptide-module-contracts.R`
  - `test-prot-03-rollup.R`
- `R/mod_prot_qc_peptide_intensity.R` now measures `452` lines with `9`
  top-level functions and is now a stabilized thin wrapper/manual target
  complete; keep the focused peptide-QC gate as the regression surface and do
  not reopen this file for more stabilization work unless a real regression
  appears.
