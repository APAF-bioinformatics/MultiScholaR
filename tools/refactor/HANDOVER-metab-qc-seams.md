# Metabolomics QC Seam Map

## Goal

Stage exact-source extraction checkpoints for `R/func_metab_qc.R` while keeping
live metabolomics QC behavior frozen behind the current public APIs.

## Current Position In The Flow

- target: `R/func_metab_qc.R`
- classification: `direct-extraction-ready`
- checkpoint reached: `wave 6 plotting helper checkpoint is now applied live into R/func_metab_qc_plotting_helpers.R`
- next step: `Manual target complete; any remaining extraction from R/func_metab_qc.R is now optional cleanup rather than required god-module stabilization work.`

## Existing Safety Net

- `tests/testthat/test-metab-01-qc-metrics-characterization.R`
- `tests/testthat/test-metab-01b-qc-filtering-helpers-characterization.R`
- `tests/testthat/test-metab-01c-qc-correlation-helper-characterization.R`
- `tests/testthat/test-metab-01d-qc-progress-helper-characterization.R`
- `tests/testthat/test-glimma-plot.R`

## Notes

Manual bucket 0 metabolomics QC stabilization target.

- This target previously had no active handover; this file now records the
  first metabolomics QC stabilization stop point.
- Classification refresh on April 16, 2026 keeps
  `R/func_metab_qc.R` at `1970` lines with `27` top-level functions,
  `391` lines for the largest top-level function, and labels
  `review` plus `direct-extraction-ready`.
- The first characterization gate now lives in
  `tests/testthat/test-metab-01-qc-metrics-characterization.R`
  and freezes the current contracts for:
  `countUniqueMetabolites()`,
  `countMetabolitesPerSample()`,
  `calculateMissingness()`,
  `calculateSumIntensityPerSample()`,
  and
  `calculateTotalUniqueMetabolitesAcrossAssays()`.
- Wave 1 manifest:
  `tools/refactor/manifest-metab-qc-wave1.yml`
  now applies the low-risk metrics helper cluster live into
  `R/func_metab_qc_metrics_helpers.R`
  and rewrites `R/func_metab_qc.R` to remove the extracted definitions.
- The staged wave-1 review artifacts remain available at
  `tools/refactor/staging/wave1_metabolomics_qc_metrics_helpers/R/func_metab_qc_metrics_helpers.R`
  and
  `tools/refactor/staging/wave1_metabolomics_qc_metrics_helpers/collate-metab-qc-wave1.txt`.
- The focused gate for this checkpoint is:
  - `tests/testthat/test-metab-01-qc-metrics-characterization.R`
- The focused gate loader now reads the extracted helper file first and then
  `R/func_metab_qc.R`, so the characterization test follows the live layout
  after apply.
- Post-apply checks passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-qc-wave1.yml`.
- The focused gate reran green after apply via a direct `testthat::test_file()`
  invocation because this worktree does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Wave 2 manifest:
  `tools/refactor/manifest-metab-qc-wave2.yml`
  now applies the low-risk filtering helper cluster live into
  `R/func_metab_qc_filtering_helpers.R`
  and rewrites `R/func_metab_qc.R` to remove:
  `metaboliteIntensityFilteringHelper()`
  and
  `resolveDuplicateFeaturesByIntensity()`.
- The staged wave-2 review artifacts remain available at
  `tools/refactor/staging/wave2_metabolomics_qc_filtering_helpers/R/func_metab_qc_filtering_helpers.R`
  and
  `tools/refactor/staging/wave2_metabolomics_qc_filtering_helpers/collate-metab-qc-wave2.txt`.
- The second focused gate now lives in
  `tests/testthat/test-metab-01b-qc-filtering-helpers-characterization.R`
  and freezes the current contracts for
  `metaboliteIntensityFilteringHelper()`
  and
  `resolveDuplicateFeaturesByIntensity()`.
- `DESCRIPTION` `Collate:` now loads
  `R/func_metab_qc_metrics_helpers.R`
  and
  `R/func_metab_qc_filtering_helpers.R`
  ahead of
  `R/func_metab_qc.R`
  so the live helper layout matches the staged extraction shape.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-qc-wave2.yml`.
- The third focused gate now lives in
  `tests/testthat/test-metab-01c-qc-correlation-helper-characterization.R`
  and freezes the current contract for
  `calculateMetabolitePairCorrelation()`.
- Wave 3 manifest:
  `tools/refactor/manifest-metab-qc-wave3.yml`
  now applies the bounded correlation helper checkpoint live into
  `R/func_metab_qc_correlation_helpers.R`
  and rewrites `R/func_metab_qc.R` to remove
  `calculateMetabolitePairCorrelation()`.
- The staged wave-3 review artifacts remain available at
  `tools/refactor/staging/wave3_metabolomics_qc_correlation_helpers/R/func_metab_qc_correlation_helpers.R`
  and
  `tools/refactor/staging/wave3_metabolomics_qc_correlation_helpers/collate-metab-qc-wave3.txt`.
- `DESCRIPTION` `Collate:` now loads
  `R/func_metab_qc_metrics_helpers.R`,
  `R/func_metab_qc_filtering_helpers.R`,
  and
  `R/func_metab_qc_correlation_helpers.R`
  ahead of
  `R/func_metab_qc.R`
  so the live helper layout matches the staged extraction shape.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-qc-wave3.yml`.
- After the wave-3 apply checkpoint, the classifier reports
  `R/func_metab_qc.R`
  at `1556` lines with `19` top-level functions while keeping
  `review` plus `direct-extraction-ready`.
- The fourth focused gate now lives in
  `tests/testthat/test-metab-01d-qc-progress-helper-characterization.R`
  and freezes the current contracts for
  `getFilteringProgressMetabolomics()`
  and
  `updateFilteringProgressMetabolomics()`.
- Wave 4 manifest:
  `tools/refactor/manifest-metab-qc-wave4.yml`
  now applies the bounded progress helper checkpoint live into
  `R/func_metab_qc_progress_helpers.R`
  and rewrites `R/func_metab_qc.R` to remove:
  `getFilteringProgressMetabolomics()`
  and
  `updateFilteringProgressMetabolomics()`.
- The staged wave-4 review artifacts remain available at
  `tools/refactor/staging/wave4_metabolomics_qc_progress_helpers/R/func_metab_qc_progress_helpers.R`
  and
  `tools/refactor/staging/wave4_metabolomics_qc_progress_helpers/collate-metab-qc-wave4.txt`.
- `DESCRIPTION` `Collate:` now loads
  `R/func_metab_qc_metrics_helpers.R`,
  `R/func_metab_qc_filtering_helpers.R`,
  `R/func_metab_qc_correlation_helpers.R`,
  and
  `R/func_metab_qc_progress_helpers.R`
  ahead of
  `R/func_metab_qc.R`
  so the live helper layout continues to match the staged extraction shape.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-qc-wave4.yml`.
- The focused gate reran green after apply via a direct `testthat::test_file()`
  invocation because this worktree does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- After the wave-4 apply checkpoint, the classifier reports
  `R/func_metab_qc.R`
  at `1478` lines with `17` top-level functions while keeping
  `review` plus `direct-extraction-ready`.
- The fifth focused gate now lives in
  `tests/testthat/test-metab-01e-qc-variability-helper-characterization.R`
  and freezes the current contracts for
  `calculateMetaboliteCVs()`
  and
  `getInternalStandardMetrics()`.
- Wave 5 manifest:
  `tools/refactor/manifest-metab-qc-wave5.yml`
  now applies the bounded variability helper checkpoint live into
  `R/func_metab_qc_variability_helpers.R`
  and rewrites `R/func_metab_qc.R` to remove:
  `calculateMetaboliteCVs()`
  and
  `getInternalStandardMetrics()`.
- The staged wave-5 review artifacts remain available at
  `tools/refactor/staging/wave5_metabolomics_qc_variability_helpers/R/func_metab_qc_variability_helpers.R`
  and
  `tools/refactor/staging/wave5_metabolomics_qc_variability_helpers/collate-metab-qc-wave5.txt`.
- `DESCRIPTION` `Collate:` now loads
  `R/func_metab_qc_metrics_helpers.R`,
  `R/func_metab_qc_filtering_helpers.R`,
  `R/func_metab_qc_correlation_helpers.R`,
  `R/func_metab_qc_progress_helpers.R`,
  and
  `R/func_metab_qc_variability_helpers.R`
  ahead of
  `R/func_metab_qc.R`
  so the live helper layout continues to match the staged extraction shape.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-qc-wave5.yml`.
- The focused gate reran green after apply via a direct `testthat::test_file()`
  invocation because this worktree does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- After the wave-5 apply checkpoint, the classifier reports
  `R/func_metab_qc.R`
  at `1173` lines with `13` top-level functions while keeping
  `review` plus `direct-extraction-ready`.
- The sixth focused gate now lives in
  `tests/testthat/test-metab-01f-qc-plotting-helper-characterization.R`
  and freezes the current plot-surface contracts for
  `generateMetaboliteFilteringPlots()`.
- Wave 6 manifest:
  `tools/refactor/manifest-metab-qc-wave6.yml`
  now applies the bounded plotting helper checkpoint live into
  `R/func_metab_qc_plotting_helpers.R`
  and rewrites `R/func_metab_qc.R` to remove
  `generateMetaboliteFilteringPlots()`.
- The staged wave-6 review artifacts remain available at
  `tools/refactor/staging/wave6_metabolomics_qc_plotting_helpers/R/func_metab_qc_plotting_helpers.R`
  and
  `tools/refactor/staging/wave6_metabolomics_qc_plotting_helpers/collate-metab-qc-wave6.txt`.
- `DESCRIPTION` `Collate:` now loads
  `R/func_metab_qc_metrics_helpers.R`,
  `R/func_metab_qc_filtering_helpers.R`,
  `R/func_metab_qc_correlation_helpers.R`,
  `R/func_metab_qc_progress_helpers.R`,
  `R/func_metab_qc_variability_helpers.R`,
  and
  `R/func_metab_qc_plotting_helpers.R`
  ahead of
  `R/func_metab_qc.R`
  so the live helper layout continues to match the staged extraction shape.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-qc-wave6.yml`.
- The focused gate reran green after apply via a direct `testthat::test_file()`
  invocation because this worktree does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- After the wave-6 apply checkpoint, the classifier reports
  `R/func_metab_qc.R`
  at `750` lines with `12` top-level functions and a `110` line largest
  top-level function while keeping `direct-extraction-ready`.
- This manual metabolomics QC target is now below the stabilization soft cap
  and is no longer a blocker in the god-module backlog.
