# Metabolomics DA Seam Map

## Goal

Stage exact-source extraction checkpoints for R/func_metab_da.R while keeping live metabolomics DA behavior frozen behind the current public APIs.

## Current Position In The Flow

- target: `/home/doktersmol/Documents/MultiScholaR-metab-safe/R/func_metab_da.R`
- classification: `direct-extraction-ready`
- checkpoint reached: `wrapper breadcrumb cleanup checkpoint reduces R/func_metab_da.R to a minimal metabolomics DA identity shell`
- next step: `No further metabolomics DA stabilization work is required for this bucket; treat the handover as archival and move later loop work to other targets.`

## Existing Safety Net

- `tests/testthat/test-metab-02-da-results-long-format-characterization.R`
- `tests/testthat/test-metab-02b-da-volcano-static-characterization.R`
- `tests/testthat/test-metab-02c-da-results-output-characterization.R`
- `tests/testthat/test-metab-02d-da-heatmap-characterization.R`
- `tests/testthat/test-metab-02e-da-volcano-glimma-characterization.R`
- `tests/testthat/test-metab-02f-da-quant-data-characterization.R`
- `tests/testthat/test-metab-02g-da-core-engine-characterization.R`
- `tests/testthat/test-metab-02h-da-orchestration-characterization.R`
- `tests/testthat/test-metab-02i-da-counts-table-characterization.R`

## Notes

Manual bucket 0 metabolomics DA stabilization target.

- This target previously had no active handover; this file now records the
  first metabolomics DA stabilization stop point.
- Classification refresh on April 16, 2026 now kept
  `R/func_metab_da.R`
  at `111` lines with `1` top-level function and label
  `direct-extraction-ready`
  after the eighth bounded helper extraction.
- A cleanup-only characterization gate now lives in
  `tests/testthat/test-metab-02i-da-counts-table-characterization.R`
  and freezes the current `getCountsTable()` contract by loading the
  metabolite S4 helper definition first when present and then falling back to
  `R/func_metab_da.R`.
- The duplicate `getCountsTable()` definition has now been removed from
  `R/func_metab_da.R`
  so the metabolomics DA wrapper shell relies on the shared helper that already
  collates from
  `R/func_metab_s4_da_results.R`.
- Classification refresh on April 16, 2026 now reduces
  `R/func_metab_da.R`
  to a `50` line wrapper breadcrumb with `0` top-level functions after the
  cleanup-only counts-table checkpoint.
- A final breadcrumb cleanup on April 16, 2026 trims
  `R/func_metab_da.R`
  to `27` lines while preserving only the public wrapper-identity comment
  surface and keeping `0` top-level functions.
- The focused characterization gate now lives in
  `tests/testthat/test-metab-02-da-results-long-format-characterization.R`
  and freezes the current contract for
  `createMetabDaResultsLongFormat()`
  by loading the extracted helper file first when present and then falling back
  to `R/func_metab_da.R`.
- The next focused characterization gate now lives in
  `tests/testthat/test-metab-02b-da-volcano-static-characterization.R`
  and freezes the current contract for
  `generateMetabDAVolcanoStatic()`
  by loading the extracted helper file first when present and then falling back
  to `R/func_metab_da.R`.
- The next focused characterization gate now lives in
  `tests/testthat/test-metab-02c-da-results-output-characterization.R`
  and freezes the current contract for
  `outputMetabDaResultsAllContrasts()`
  by loading the extracted helper file first when present, loading the live
  static-volcano helper when present, and stubbing the heatmap branch so the
  output/export surface stays characterized in this worktree without optional
  `ComplexHeatmap` and `circlize`.
- The next focused characterization gate now lives in
  `tests/testthat/test-metab-02d-da-heatmap-characterization.R`
  and freezes the current early-exit and assay-filtering contract for
  `generateMetabDAHeatmap()`
  by loading the extracted helper file first when present and exercising the
  live pre-`ComplexHeatmap`/`circlize` branches that remain available in this
  worktree without those optional packages.
- The next focused characterization gate now lives in
  `tests/testthat/test-metab-02e-da-volcano-glimma-characterization.R`
  and freezes the current early-exit, combined-view, missing-data, and
  missing-coefficient guards for
  `generateMetabDAVolcanoPlotGlimma()`
  by loading the extracted helper file first when present and exercising the
  live pre-widget branches that remain available in this worktree without
  optional Glimma stack packages.
- Wave 1 manifest:
  `tools/refactor/manifest-metab-da-wave1.yml`
  now applies the bounded metabolomics DA long-format helper checkpoint live into
  `R/func_metab_da_results_long_format.R`
  and rewrites `R/func_metab_da.R` to remove
  `createMetabDaResultsLongFormat()`.
- The staged wave-1 review artifacts remain available at
  `tools/refactor/staging/wave1_metabolomics_da_results_long_format/R/func_metab_da_results_long_format.R`
  and
  `tools/refactor/staging/wave1_metabolomics_da_results_long_format/collate-metab-da-wave1.txt`.
- Wave 2 manifest:
  `tools/refactor/manifest-metab-da-wave2.yml`
  now applies the bounded metabolomics DA static-volcano helper checkpoint live
  into
  `R/func_metab_da_volcano_static.R`
  and rewrites `R/func_metab_da.R` to remove
  `generateMetabDAVolcanoStatic()`.
- The staged wave-2 review artifacts remain available at
  `tools/refactor/staging/wave2_metabolomics_da_volcano_static/R/func_metab_da_volcano_static.R`
  and
  `tools/refactor/staging/wave2_metabolomics_da_volcano_static/collate-metab-da-wave2.txt`.
- Wave 3 manifest:
  `tools/refactor/manifest-metab-da-wave3.yml`
  now applies the bounded metabolomics DA results-output helper checkpoint live
  into
  `R/func_metab_da_results_output.R`
  and rewrites `R/func_metab_da.R` to remove
  `outputMetabDaResultsAllContrasts()`.
- The staged wave-3 review artifacts remain available at
  `tools/refactor/staging/wave3_metabolomics_da_results_output/R/func_metab_da_results_output.R`
  and
  `tools/refactor/staging/wave3_metabolomics_da_results_output/collate-metab-da-wave3.txt`.
- Wave 4 manifest:
  `tools/refactor/manifest-metab-da-wave4.yml`
  now applies the bounded metabolomics DA heatmap helper checkpoint live into
  `R/func_metab_da_heatmap.R`
  and rewrites `R/func_metab_da.R` to remove
  `generateMetabDAHeatmap()`.
- The staged wave-4 review artifacts remain available at
  `tools/refactor/staging/wave4_metabolomics_da_heatmap/R/func_metab_da_heatmap.R`
  and
  `tools/refactor/staging/wave4_metabolomics_da_heatmap/collate-metab-da-wave4.txt`.
- Wave 5 manifest:
  `tools/refactor/manifest-metab-da-wave5.yml`
  now applies the bounded metabolomics DA Glimma-volcano helper checkpoint live into
  `R/func_metab_da_volcano_glimma.R`
  and rewrites `R/func_metab_da.R` to remove
  `generateMetabDAVolcanoPlotGlimma()`.
- The staged wave-5 review artifacts remain available at
  `tools/refactor/staging/wave5_metabolomics_da_volcano_glimma/R/func_metab_da_volcano_glimma.R`
  and
  `tools/refactor/staging/wave5_metabolomics_da_volcano_glimma/collate-metab-da-wave5.txt`.
- The next focused characterization gate now lives in
  `tests/testthat/test-metab-02f-da-quant-data-characterization.R`
  and freezes the current metadata-filtering and numeric-column contract for
  `getMetaboliteQuantData()`
  by loading the extracted helper file first when present and then falling back
  to `R/func_metab_da.R`.
- Wave 6 manifest:
  `tools/refactor/manifest-metab-da-wave6.yml`
  now applies the bounded metabolomics DA quant-data helper checkpoint live into
  `R/func_metab_da_quant_data.R`
  and rewrites `R/func_metab_da.R` to remove
  `getMetaboliteQuantData()`.
- The staged wave-6 review artifacts remain available at
  `tools/refactor/staging/wave6_metabolomics_da_quant_data/R/func_metab_da_quant_data.R`
  and
  `tools/refactor/staging/wave6_metabolomics_da_quant_data/collate-metab-da-wave6.txt`.
- The next focused characterization gate now lives in
  `tests/testthat/test-metab-02g-da-core-engine-characterization.R`
  and freezes the current pre-limma no-common-samples and
  contrast-level-validation guards for
  `runTestsContrastsMetabDA()`
  by loading the extracted helper file first when present and then falling back
  to `R/func_metab_da.R`.
- The next focused characterization gate now lives in
  `tests/testthat/test-metab-02h-da-orchestration-characterization.R`
  and freezes the current input guard, assay aggregation, unmatched-assay skip,
  and qvalue-warning capture contract for
  `runMetabolitesDA()`
  by loading the extracted helper file first when present and then falling back
  to `R/func_metab_da.R`.
- Wave 7 manifest:
  `tools/refactor/manifest-metab-da-wave7.yml`
  now applies the bounded metabolomics DA core-engine helper checkpoint live into
  `R/func_metab_da_core_engine.R`
  and rewrites `R/func_metab_da.R` to remove
  `runTestsContrastsMetabDA()`.
- The staged wave-7 review artifacts remain available at
  `tools/refactor/staging/wave7_metabolomics_da_core_engine/R/func_metab_da_core_engine.R`
  and
  `tools/refactor/staging/wave7_metabolomics_da_core_engine/collate-metab-da-wave7.txt`.
- Wave 8 manifest:
  `tools/refactor/manifest-metab-da-wave8.yml`
  now applies the bounded metabolomics DA orchestration helper checkpoint live into
  `R/func_metab_da_orchestration.R`
  and rewrites `R/func_metab_da.R` to remove
  `runMetabolitesDA()`.
- The staged wave-8 review artifacts remain available at
  `tools/refactor/staging/wave8_metabolomics_da_orchestration/R/func_metab_da_orchestration.R`
  and
  `tools/refactor/staging/wave8_metabolomics_da_orchestration/collate-metab-da-wave8.txt`.
- `DESCRIPTION` `Collate:` now loads
  `R/func_metab_da_results_long_format.R`
  and
  `R/func_metab_da_volcano_static.R`
  and
  `R/func_metab_da_volcano_glimma.R`
  and
  `R/func_metab_da_heatmap.R`
  and
  `R/func_metab_da_results_output.R`
  and
  `R/func_metab_da_quant_data.R`
  and
  `R/func_metab_da_core_engine.R`
  ahead of
  `R/func_metab_da.R`
  so the live helper layout matches the staged extraction shape.
- Post-apply checks passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-da-wave1.yml`.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-da-wave2.yml`.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-da-wave3.yml`.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-da-wave4.yml`.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-da-wave5.yml`.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-da-wave6.yml`.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-da-wave7.yml`.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-da-wave8.yml`.
- The focused gate reran green after apply via a direct `testthat::test_file()`
  invocation because the existing full-package Glimma gate cannot load in this
  worktree without optional packages such as `limma` and `Glimma`.
- The focused cleanup gate reran green after the breadcrumb trim via a direct
  `testthat::test_file()` invocation of
  `tests/testthat/test-metab-02i-da-counts-table-characterization.R`.
- `DESCRIPTION` `Collate:` now also loads
  `R/func_metab_da_orchestration.R`
  ahead of
  `R/func_metab_da.R`
  so the live helper layout matches the staged extraction shape for the
  orchestration seam.
