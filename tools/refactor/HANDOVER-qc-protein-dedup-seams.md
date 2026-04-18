# Protein Duplicate Removal Seam Map

## Goal

Introduce one safe in-file seam at a time in `R/mod_prot_qc_protein_dedup.R`
while preserving the module's public contract.

## Current Position In The Flow

- target: `/home/doktersmol/Documents/MultiScholaR/R/mod_prot_qc_protein_dedup.R`
- classification: `review`
- checkpoint reached: `protein_duplicate_removal_render_path_seam`
- next step: `Manual target complete; keep this wrapper archived unless a future staged extraction is explicitly restarted.`

## Existing Safety Net

- `tests/testthat/test-prot-02c-qc-protein-module-contracts.R`

## Notes

- This manual target previously had no dedicated handover; this file now marks
  the first explicit stop point for the duplicate-removal wrapper.
- April 16, 2026 stabilize iteration landed one bounded live apply-path seam in
  `R/mod_prot_qc_protein_dedup.R` via:
  `runProteinDuplicateRemovalStep()`,
  `updateProteinDuplicateRemovalOutputs()`,
  and
  `runProteinDuplicateRemovalApplyObserver()`.
- The new helper cluster now owns duplicate detection, aggregation-method
  resolution, state-save metadata, summary-text assembly, plot refresh, and
  apply success/error notification handling while the wrapper keeps the public
  module id, revert observer wiring, and plot render binding unchanged.
- The same checkpoint rewired the inline `apply_duplicate_removal` observer
  through the new top-level helper shell without changing the public apply
  button contract.
- The same checkpoint added direct characterization coverage in
  `tests/testthat/test-prot-02c-qc-protein-module-contracts.R`
  for the new apply step, output-refresh helper, apply observer happy/error
  paths, and wrapper-level apply delegation.
- April 16, 2026 stabilize follow-up landed one bounded live revert-path seam
  in `R/mod_prot_qc_protein_dedup.R` via
  `runProteinDuplicateRemovalRevertStep()` and
  `runProteinDuplicateRemovalRevertObserver()`.
- The new helper pair now owns revert-path history lookup, prior-state
  selection, state-manager rollback, result-text assembly, success logging, and
  revert success/error notification handling while the wrapper keeps the public
  module id, apply observer wiring, and plot render binding unchanged.
- The same checkpoint rewired the inline `revert_duplicate_removal` observer
  through the new top-level helper shell without changing the public revert
  button contract.
- The same checkpoint added direct characterization coverage in
  `tests/testthat/test-prot-02c-qc-protein-module-contracts.R`
  for the revert helper's success and no-history error paths, revert observer
  happy/error paths, and wrapper-level revert delegation.
- Focused gate reran green after the live seam for:
  - `test-prot-02c-qc-protein-module-contracts`
- April 16, 2026 stabilize follow-up landed one bounded live render-path seam
  in `R/mod_prot_qc_protein_dedup.R` via
  `bindProteinDuplicateRemovalPlot()`.
- The new helper now owns the `renderPlot()` output binding, stored-grid
  requirement, and `grid.draw()` delegation while the wrapper keeps the public
  module id plus the apply and revert observer wiring unchanged.
- The same checkpoint rewired the wrapper's final `duplicate_removal_plot`
  binding through the new top-level helper shell without changing the public
  plot output contract.
- The same checkpoint added direct characterization coverage in
  `tests/testthat/test-prot-02c-qc-protein-module-contracts.R`
  for the bind helper's render assignment and stored-grid draw path, plus
  wrapper-level delegation of the bind helper.
- Focused gate reran green after the live seam for:
  - `test-prot-02c-qc-protein-module-contracts`
- `R/mod_prot_qc_protein_dedup.R` now measures `297` live lines with `8`
  top-level functions and remains classified as `review`, but this manual
  target is complete and should stay archived behind this handover.
