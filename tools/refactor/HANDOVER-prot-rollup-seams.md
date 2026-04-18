# Proteomics Rollup Seam Map

## Goal

Stage exact-source extraction checkpoints for R/func_prot_rollup.R while keeping the live rollup entry points behaviorally frozen.

## Current Position In The Flow

- target: `/home/doktersmol/Documents/MultiScholaR/R/func_prot_rollup.R`
- classification: `direct-extraction-ready`
- checkpoint reached: `wave2_proteomics_rollup_methods_live_apply`
- next step: `Manual rollup target is complete; keep this handover as the archival seam record.`

## Existing Safety Net

- `tests/testthat/test-prot-03-rollup.R`

## Notes

- Wave 1 manifest:
  `tools/refactor/manifest-prot-rollup-wave1.yml`
  now applies live into
  `R/func_prot_rollup_count_helpers.R`,
  covering:
  `calcPeptidesPerProtein()`,
  `calcTotalPeptides()`,
  `countPeptidesPerRun()`,
  `count_num_peptides()`,
  `countProteinsPerRun()`,
  `countUniqueProteins()`,
  `count_num_proteins()`,
  and
  `count_num_samples()`.
- The staged review artifact remains at
  `tools/refactor/staging/wave1_proteomics_rollup_count_helpers/R/func_prot_rollup_count_helpers.R`,
  and the apply-time collate artifact now exists at
  `tools/refactor/collate-prot-rollup-wave1.txt`.
- `DESCRIPTION` `Collate:` now includes `R/func_prot_rollup_count_helpers.R`
  immediately before `R/func_prot_rollup.R`.
- The extracted count/reporting helper file is `252` lines and the trimmed live
  `R/func_prot_rollup.R` tail is now `229` lines.
- Focused gate reran green after the live apply checkpoint:
  - `tests/testthat/test-prot-03-rollup.R`
- The remaining live tail now holds only `rollUpPrecursorToPeptideHelper()` and
  the `rollUpPrecursorToPeptide` S4 method shell, so the next safe stop point
  is a reviewed second extraction wave for that pair rather than another count
  helper move.
- Wave 2 manifest:
  `tools/refactor/manifest-prot-rollup-wave2.yml`
  now applies live into
  `R/func_prot_rollup_methods.R`,
  covering:
  `rollUpPrecursorToPeptideHelper()`
  and the
  `rollUpPrecursorToPeptide` S4 method shell.
- The staged review artifact now exists at
  `tools/refactor/staging/wave2_proteomics_rollup_methods/R/func_prot_rollup_methods.R`,
  and the apply-time collate artifact now exists at
  `tools/refactor/collate-prot-rollup-wave2.txt`.
- `DESCRIPTION` `Collate:` now includes `R/func_prot_rollup_methods.R`
  immediately before `R/func_prot_rollup.R`.
- The extracted rollup methods file is `95` lines and
  `R/func_prot_rollup.R` is now a `28` line breadcrumb stub.
- Focused gate reran green after the live apply checkpoint:
  - `tests/testthat/test-prot-03-rollup.R`
- The manual rollup target is now complete in live `R/`.
