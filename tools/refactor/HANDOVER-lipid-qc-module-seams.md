# Lipid QC Module Seam Map

## Goal

Document the wrapper-level stabilization stop point for mod_lipid_qc.R while keeping the public lipid-QC module contract stable.

## Current Position In The Flow

- target: `/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc.R`
- classification: `review`
- active stop point: `initializeLipidQcSubmodules()` is now live in the wrapper. Both the startup-trigger and auto-init paths route through that seam.
- next step: `No further live seam work is required for this backlog target.`

## Existing Safety Net

- focused wrapper gate:
  - `tests/testthat/test-lipid-01b-qc-module-contracts.R`
- replay command:
  - `Rscript -e "testthat::test_file('tests/testthat/test-lipid-01b-qc-module-contracts.R', stop_on_failure = TRUE)"`

## Notes

- No prior handover existed for this wrapper lane.
- The new source-driven gate freezes:
  - the UI heading and `dynamic_qc_tabs` placeholder contract
  - the empty-state alert when no lipid workflow state is available
  - the auto-init path when the state manager already holds a
    `LipidomicsAssayData` object
  - the trigger-init path when `qc_trigger` starts `TRUE`
- The harness intentionally evaluates the live source with a local
  `moduleServer()` capture so the wrapper body can run under
  `MockShinySession` without mutating live `R/`.
- April 15, 2026 introduced one bounded live seam in
  `/home/doktersmol/Documents/MultiScholaR-lipid-lane/R/mod_lipid_qc.R`
  by extracting the duplicated four-submodule initialization block into
  `initializeLipidQcSubmodules()`.
- The trigger-driven and state-detected initialization branches now share that
  helper without changing the public wrapper contract or the existing log flow.
- The focused wrapper gate reran green after the live seam.
- Post-seam classification metrics are `153` lines, `3` top-level functions,
  and max top-level function length `69`.
- The classifier still reports `review`, but this wrapper is now inside the
  ideal size band and no further live seam is required for the current backlog
  target.
