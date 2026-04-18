# Protein Accession Cleanup Seam Map

## Goal

Introduce one safe in-file seam at a time in `R/mod_prot_qc_protein_cleanup.R`
while preserving the module's public contract.

## Current Position In The Flow

- target: `/home/doktersmol/Documents/MultiScholaR/R/mod_prot_qc_protein_cleanup.R`
- classification: `review`
- checkpoint reached: `protein_cleanup_wrapper_completion`
- next step: `Treat the protein accession cleanup wrapper target as complete and archived. Keep test-prot-02c-qc-protein-module-contracts.R as the regression surface and do not reopen mod_prot_qc_protein_cleanup.R unless a real behavior regression appears.`

## Existing Safety Net

- `tests/testthat/test-prot-02c-qc-protein-module-contracts.R`

## Notes

- This manual target previously had no dedicated handover; this file now marks
  the first explicit stop point for the protein accession cleanup wrapper.
- April 16, 2026 stabilize iteration landed one bounded live apply-path seam in
  `R/mod_prot_qc_protein_cleanup.R` via:
  `runProteinAccessionCleanupStep()`,
  `updateProteinAccessionCleanupOutputs()`,
  and
  `runProteinAccessionCleanupApplyObserver()`.
- The new helper cluster now owns the FASTA-backed accession cleanup step,
  workflow/QC result tracking, state-save metadata, result text assembly, plot
  refresh, and apply success/error notification handling while the wrapper
  keeps the public module id, revert observer wiring, and plot render binding
  unchanged.
- The same checkpoint rewired the inline `apply_accession_cleanup` observer
  through the new top-level helper shell without changing the public apply
  button contract.
- The same checkpoint added direct characterization coverage in
  `tests/testthat/test-prot-02c-qc-protein-module-contracts.R`
  for the new apply step, output-refresh helper, apply observer happy/error
  paths, and wrapper-level apply delegation.
- April 16, 2026 stabilize iteration landed one bounded live revert-path seam
  in `R/mod_prot_qc_protein_cleanup.R` via:
  `runProteinAccessionCleanupRevertStep()`
  and
  `runProteinAccessionCleanupRevertObserver()`.
- The new helper pair now owns the revert-path history lookup, previous-state
  selection, state-manager rollback, summary text construction, success log,
  and success/error notification handling while the wrapper keeps the public
  module id, apply observer wiring, and plot render binding unchanged.
- The same checkpoint rewired the inline `revert_accession_cleanup` observer
  through the new top-level helper shell without changing the public revert
  button contract.
- The same checkpoint added direct characterization coverage in
  `tests/testthat/test-prot-02c-qc-protein-module-contracts.R`
  for the revert helper's success and no-history error paths, revert observer
  happy/error paths, and wrapper-level revert delegation.
- Focused gate reran green after the live seam for:
  - `tests/testthat/test-prot-02c-qc-protein-module-contracts.R`
- Live `R/mod_prot_qc_protein_cleanup.R` now measures `397` lines with `7`
  top-level functions and currently classifies as `review`.
- April 16, 2026 bounded stabilize follow-up extracted
  `bindProteinAccessionCleanupPlot()` from the remaining `renderPlot()`
  binding in `R/mod_prot_qc_protein_cleanup.R`.
- The new helper now owns the plot-output render binding, stored-grid
  requirement, and `grid.draw()` delegation while leaving the wrapper as a
  thin module shell that only creates the shared reactive plot state and wires
  the public apply/revert observers.
- The same checkpoint added direct characterization coverage in
  `tests/testthat/test-prot-02c-qc-protein-module-contracts.R` for the new
  render-binding helper's `renderPlot()` assignment and stored-grid draw path,
  plus wrapper-level delegation of the bind helper.
- Focused gate reran green after the render-path seam for:
  - `tests/testthat/test-prot-02c-qc-protein-module-contracts.R`
- Live `R/mod_prot_qc_protein_cleanup.R` now measures `403` lines with `8`
  top-level functions and still classifies as `review`, but this manual target
  is complete: the wrapper is now a stabilized thin shell with helper-level
  seams and wrapper-level contract coverage, so future work should move to the
  next backlog item instead of reopening `mod_prot_qc_protein_cleanup.R` for
  more stabilization.
