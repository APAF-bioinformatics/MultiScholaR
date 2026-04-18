# Metabolomics QC Module Seam Map

## Goal

Introduce bounded top-level seams in R/mod_metab_qc.R while keeping the live metabolomics QC module behavior frozen behind the existing public entry points.

## Current Position In The Flow

- target: `/home/doktersmol/Documents/MultiScholaR-metab-safe/R/mod_metab_qc.R`
- classification: `review`
- checkpoint reached: `first bounded sub-module initialization seam is now live in R/mod_metab_qc.R`
- next step: `Treat R/mod_metab_qc.R as the completed metabolomics QC wrapper identity for this bucket; any later extraction from the dynamic-tab render shell is optional cleanup rather than required stabilization work.`

## Existing Safety Net

- `tests/testthat/test-metab-01o-qc-module-characterization.R`

## Notes

Manual bucket 0 metabolomics QC module stabilization target.

- This target previously had no active handover; this file now records the
  first metabolomics QC module seam stop point.
- Classification refreshed on April 17, 2026 keeps
  `R/mod_metab_qc.R`
  at `167` lines with `3` top-level functions, a `81` line largest top-level
  function, and label `review`.
- The first focused gate now lives in
  `tests/testthat/test-metab-01o-qc-module-characterization.R`
  and freezes the current contracts for
  `initializeMetabQcSubmodules()`
  plus the wrapper handoff from
  `mod_metab_qc_server()`
  into that seam across both the `qc_trigger()` path and the state-detected
  auto-initialization path.
- The first bounded live seam now sits in
  `R/mod_metab_qc.R`
  at line `51`
  as `initializeMetabQcSubmodules()`, which now owns the duplicated
  intensity/duplicates/internal-standard/finalization server registration path
  before control returns to the two initialization branches inside
  `mod_metab_qc_server()`.
- The `qc_trigger()` observer branch in
  `R/mod_metab_qc.R`
  at line `123`
  now delegates through that helper instead of keeping the full four-module
  registration block inline inside the observer body.
- The state-detected auto-initialization branch in
  `R/mod_metab_qc.R`
  at line `141`
  now delegates through the same helper instead of repeating the same
  registration block inline inside the fallback observer.
- The focused gate reran green after the new seam via direct
  `testthat::test_file()` invocation because this worktree does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- The metabolomics QC wrapper now sits comfortably inside the playbook's ideal
  size band as a public orchestrator shell, so this manual target no longer
  blocks bucket 0 stabilization work.
