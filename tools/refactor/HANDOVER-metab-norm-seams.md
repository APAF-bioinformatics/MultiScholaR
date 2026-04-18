# Metabolomics Normalization Seam Map

## Goal

Stage exact-source extraction checkpoints for `R/func_metab_norm.R` while
keeping live metabolomics normalization helpers frozen behind the current
public APIs.

## Current Position In The Flow

- target: `R/func_metab_norm.R`
- classification: `direct-extraction-ready`
- checkpoint reached: `wave 2 ITSD selection helper checkpoint is now applied live into R/func_metab_norm_support_helpers.R`
- next step: `Metabolomics normalization helper stabilization is complete for the current backlog target; any later extraction of generateMetabQcPlots() is optional cleanup rather than a required stop point.`

## Existing Safety Net

- `tests/testthat/test-metab-03-norm-support-characterization.R`

## Notes

Manual bucket 0 metabolomics normalization stabilization target.

- This target previously had no active handover; this file now records the
  first metabolomics normalization stabilization stop point.
- Classification refresh on April 16, 2026 kept
  `R/func_metab_norm.R`
  at `685` lines with `7` top-level functions and label
  `direct-extraction-ready`
  before the first extraction checkpoint.
- The first focused characterization gate now lives in
  `tests/testthat/test-metab-03-norm-support-characterization.R`
  and freezes the current contracts for:
  `extractBestKPerAssay()`,
  `extractCtrlPerAssay()`,
  `buildCombinedRuvTable()`,
  and
  `buildNormConfig()`.
- Wave 1 manifest:
  `tools/refactor/manifest-metab-norm-wave1.yml`
  now applies the low-risk metabolomics normalization support helper cluster
  live into
  `R/func_metab_norm_support_helpers.R`
  and rewrites `R/func_metab_norm.R` to remove the extracted definitions.
- The staged wave-1 review artifacts remain available at
  `tools/refactor/staging/wave1_metabolomics_norm_support_helpers/R/func_metab_norm_support_helpers.R`
  and
  `tools/refactor/staging/wave1_metabolomics_norm_support_helpers/collate-metab-norm-wave1.txt`.
- Manifest verification passed through
  `tools/refactor/verify_refactor.R --manifest tools/refactor/manifest-metab-norm-wave1.yml`.
- `DESCRIPTION` `Collate:` now loads
  `R/func_metab_norm_support_helpers.R`
  ahead of
  `R/func_metab_norm.R`
  so the live helper layout matches the applied extraction shape.
- Post-apply checks passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-norm-wave1.yml`.
- After the wave-1 apply checkpoint,
  `R/func_metab_norm.R`
  is down to `597` lines with `3` top-level functions and a `209` line largest
  top-level function while remaining `direct-extraction-ready`.
- The new live helper file
  `R/func_metab_norm_support_helpers.R`
  is `92` lines with `4` top-level functions and stays within the default size
  budget.
- The focused gate reran green after the live wave-1 checkpoint via a direct
  `testthat::test_file()` invocation because this worktree does not include
  `renv/activate.R` for `tools/test_with_renv.R`.
- The focused characterization gate now also freezes the current contracts for
  `buildItsdSelectionTable()`,
  including the multi-annotation ordering path and the empty-sample warning
  fallback.
- Wave 2 manifest:
  `tools/refactor/manifest-metab-norm-wave2.yml`
  now applies the bounded ITSD selection helper checkpoint live into the
  existing helper file
  `R/func_metab_norm_support_helpers.R`
  and rewrites `R/func_metab_norm.R` to remove the extracted definition.
- The staged wave-2 review artifacts remain available at
  `tools/refactor/staging/wave2_metabolomics_norm_itsd_selection_table/R/func_metab_norm_support_helpers.R`
  and
  `tools/refactor/staging/wave2_metabolomics_norm_itsd_selection_table/collate-metab-norm-wave2.txt`.
- Manifest verification passed through
  `tools/refactor/verify_refactor.R --manifest tools/refactor/manifest-metab-norm-wave2.yml`.
- Post-apply checks passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-norm-wave2.yml`.
- After the wave-2 apply checkpoint,
  `R/func_metab_norm.R`
  is down to `473` lines with `2` top-level functions while
  `R/func_metab_norm_support_helpers.R`
  is up to `217` lines with `5` top-level functions.
- The focused gate reran green again after the live wave-2 checkpoint via a
  direct `testthat::test_file()` invocation because this worktree does not
  include `renv/activate.R` for `tools/test_with_renv.R`.
- With the live wrapper file now inside the default `150-500` line budget,
  this manual bucket 0 metabolomics normalization helper target is no longer a
  stabilization blocker; any remaining direct extraction, such as
  `generateMetabQcPlots()`, is optional follow-up rather than required work.
