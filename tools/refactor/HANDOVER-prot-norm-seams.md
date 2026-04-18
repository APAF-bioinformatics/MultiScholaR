# Proteomics Normalization Function Seam Map

## Goal

Stage exact-source helper extraction checkpoints for `R/func_prot_norm.R` while keeping live normalization and RUV behavior frozen behind the current exported function surface.

## Current Position In The Flow

- target: `/home/doktersmol/Documents/MultiScholaR/R/func_prot_norm.R`
- classification: `direct-extraction-ready`
- checkpoint reached: `prot_norm_wave3_live_apply`
- next step: `No further func_prot_norm.R extraction work is pending; future normalization stabilization should resume from the remaining mod_prot_norm.R observer/render shell seams.`

## Existing Safety Net

- `tests/testthat/test-prot-05-normalisation.R`
- `tests/testthat/test-prot-05b-norm-module-contracts.R`
- `tests/testthat/test-prot-06-ruv.R`

## Notes

Manual target bucket 0 for proteomics normalization function stabilization.

- This target previously had no dedicated handover; this file now records the first `func_prot_norm.R` stop point.
- April 16, 2026 wave 1 was staged, reviewed, and applied live from
  `tools/refactor/manifest-prot-norm-wave1.yml` into
  `tools/refactor/staging/prot-norm-wave1/R/func_prot_norm_optimization_helpers.R`
  with emitted collate artifact
  `tools/refactor/collate-prot-norm-wave1.txt`.
- The live wave extracts the RUV optimization helper cluster from
  `R/func_prot_norm.R`:
  - `findBestK()`
  - `findBestKForAssayList()`
  - `calculateSeparationScore()`
  - `calculateCompositeScore()`
  - `calculateAdaptiveMaxK()`
- `R/func_prot_norm_optimization_helpers.R` is now live at `283` lines and
  `R/func_prot_norm.R` is trimmed to `691` lines while preserving the existing
  exported normalization surface.
- `DESCRIPTION` `Collate:` now includes
  `R/func_prot_norm_optimization_helpers.R` immediately before
  `R/func_prot_norm.R`.
- The focused normalization gate stayed green before apply and after the live
  apply checkpoint, with the same expected snapshot skips for `cp05` and
  `cp06`.
- April 16, 2026 wave 2 was staged, reviewed, and applied live from
  `tools/refactor/manifest-prot-norm-wave2.yml` with staged review artifact
  `tools/refactor/staging/prot-norm-wave2/R/func_prot_norm_optimization_helpers.R`
  and emitted collate artifact
  `tools/refactor/collate-prot-norm-wave2.txt`.
- The live wave extracts the remaining helper tail from `R/func_prot_norm.R`:
  - `updateRuvParameters()`
  - `getRuvIIIReplicateMatrixHelper()`
  - `getNegCtrlProtAnovaHelper()`
  - `extractRuvResults()`
  - `scaleCenterAndFillMissing()`
- `R/func_prot_norm_optimization_helpers.R` is now live at `523` lines and
  `R/func_prot_norm.R` is trimmed to `456` lines, with the latter now carrying
  only the remaining `findBestNegCtrlPercentage()` orchestration body plus the
  historical file-header notes.
- The focused normalization gate stayed green after the wave 2 live apply, with
  the same expected snapshot skips for `cp05` and `cp06`.
- April 16, 2026 wave 3 was staged, reviewed, and applied live from
  `tools/refactor/manifest-prot-norm-wave3.yml` with staged review artifact
  `tools/refactor/staging/prot-norm-wave3/R/func_prot_norm_optimization_helpers.R`
  and emitted collate artifact
  `tools/refactor/collate-prot-norm-wave3.txt`.
- The live wave extracts the remaining orchestration body from
  `R/func_prot_norm.R`:
  - `findBestNegCtrlPercentage()`
- `R/func_prot_norm_optimization_helpers.R` is now live at `810` lines and
  `R/func_prot_norm.R` is trimmed to `170` lines, with the live source now
  reduced to historical header and TODO breadcrumb notes only.
- The focused normalization gate reran green after the wave 3 live apply, with
  the same expected snapshot skips for `cp05` and `cp06`.
- `func_prot_norm.R` is now complete as a stabilization target; the remaining
  proteomics normalization backlog work stays in `R/mod_prot_norm.R`.
