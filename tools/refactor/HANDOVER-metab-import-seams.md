# Metabolomics Import Seam Map

## Goal

Stage exact-source extraction checkpoints for `R/func_metab_import.R` while
keeping live metabolomics import behavior frozen behind the current public
helpers and module entry points.

## Current Position In The Flow

- target: `R/func_metab_import.R`
- classification: `direct-extraction-ready`
- checkpoint reached: `wave 1 import detection and column-mapping helper checkpoint is now applied live into R/func_metab_import_detection.R`
- next step: `Manual target complete; any further extraction from R/func_metab_import.R is optional follow-up rather than required god-module stabilization work.`

## Existing Safety Net

- `tests/testthat/test-metab-00-import-detection-characterization.R`
- `tests/testthat/test-metab-01-qc-metrics-characterization.R`
- `tests/testthat/test-metab-01e-qc-variability-helper-characterization.R`
- `tests/testthat/test-metab-02f-da-quant-data-characterization.R`

## Notes

Manual bucket 0 metabolomics import stabilization target.

- This target previously had no active handover; this file now records the
  first metabolomics import stabilization stop point.
- Classification refresh on April 16, 2026 keeps
  `R/func_metab_import.R`
  at `683` lines with `9` top-level functions, a `143` line largest top-level
  function, and label `direct-extraction-ready`.
- The first focused characterization gate now lives in
  `tests/testthat/test-metab-00-import-detection-characterization.R`
  and freezes the current contracts for:
  `detectMetabolomicsFormat()`,
  `findMetabMatchingColumn()`,
  `getMetabolomicsColumnDefaults()`,
  `validateColumnMapping()`,
  and
  `validateMetabColumnMapping()`.
- Wave 1 manifest:
  `tools/refactor/manifest-metab-import-wave1.yml`
  now applies the low-risk metabolomics import detection and column-mapping
  helper cluster live into
  `R/func_metab_import_detection.R`
  and rewrites `R/func_metab_import.R` to remove the extracted definitions.
- The staged wave-1 review artifacts remain available at
  `tools/refactor/staging/wave1_metabolomics_import_detection/R/func_metab_import_detection.R`
  and
  `tools/refactor/staging/wave1_metabolomics_import_detection/collate-metab-import-wave1.txt`.
- The staged wave-1 collate artifact now exists at
  `tools/refactor/staging/wave1_metabolomics_import_detection/collate-metab-import-wave1.txt`.
- Manifest verification passed through
  `tools/refactor/verify_refactor.R --manifest tools/refactor/manifest-metab-import-wave1.yml`.
- `DESCRIPTION` `Collate:` now loads
  `R/func_metab_import_detection.R`
  ahead of
  `R/func_metab_import.R`
  so the live helper layout matches the applied extraction shape.
- Post-apply checks passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-import-wave1.yml`.
- After the wave-1 apply checkpoint,
  `R/func_metab_import.R`
  is down to `356` lines with `4` top-level functions and is now within the
  default size budget.
- The new live helper file
  `R/func_metab_import_detection.R`
  is `332` lines with `5` top-level functions and stays within the default
  size budget.
- The focused gate reran green after the live wave-1 checkpoint via a direct
  `testthat::test_file()` invocation because this worktree does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- This manual metabolomics import target is now below the stabilization soft
  cap and is no longer a blocker in the god-module backlog.
