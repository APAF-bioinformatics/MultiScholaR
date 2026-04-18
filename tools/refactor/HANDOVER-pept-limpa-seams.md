# Peptide Limpa Seam Map

## Goal

Reconcile R/func_pept_limpa.R as a breadcrumb stub while keeping live peptide limpa behavior frozen behind the current public APIs.

## Current Position In The Flow

- target: `/home/doktersmol/Documents/MultiScholaR/R/func_pept_limpa.R`
- classification: `direct-extraction-ready`
- checkpoint reached: `breadcrumb stub reconciled after the peptide missingness extraction so R/func_pept_limpa.R now points only at the live peptide limpa owners`
- next step: `No further seams remain in R/func_pept_limpa.R; if later cleanup is needed, classify R/func_peptide_qc_imputation.R separately for the remaining legacy imputation inventory.`

## Existing Safety Net

- `tests/testthat/test-prot-04-design.R`

## Notes

Manual target bucket 0 for peptide limpa stabilization.

- This manual target previously had no active handover; this file now records
  the first peptide-limpa stabilization stop point.
- The exact `peptideMissingValueImputationLimpa()` `PeptideQuantitativeData`
  method is already live in `R/func_pept_s4_missingness.R`.
- `DESCRIPTION` still collates `func_pept_limpa.R` as the breadcrumb identity
  for this manual target, after `func_pept_s4_missingness.R` and before the
  broader imputation helper file.
- The focused design gate reran green with `1572` passes and the same expected
  Git LFS snapshot skip:
  `Rscript tools/test_with_renv.R tests/testthat/test-prot-04-design.R`
- Live `R/func_pept_limpa.R` is now a `38`-line breadcrumb stub with `0`
  top-level functions and no stale future-extraction note; it points only at
  `R/func_pept_s4_missingness.R` for the active peptide limpa method and
  `R/func_peptide_qc_imputation.R` for the shared imputation helper surface.
- `R/func_pept_limpa.R` no longer needs additional stabilization work.
  Any later retirement of stale duplicate peptide-limpa declarations should
  move to `R/func_peptide_qc_imputation.R` under a fresh classification and
  handover pass.
