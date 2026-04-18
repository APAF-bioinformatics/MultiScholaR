# Metabolomics S4 Seam Map

## Goal

Stage exact-source extraction checkpoints for R/func_metab_s4_objects.R while keeping live metabolomics S4 behavior frozen behind the current public APIs.

## Current Position In The Flow

- target: `R/func_metab_s4_objects.R`
- classification: `direct-extraction-ready`
- checkpoint reached: `the bounded metabolomics S4 ruvIII_C_Varying() normalization-method wave is now live from R/func_metab_s4_objects.R into R/func_metab_s4_norm_methods.R`
- next step: `No additional live method blocks remain in R/func_metab_s4_objects.R; treat the file as a stabilized breadcrumb/generics shell and shift later bucket work to a different target if needed.`

## Existing Safety Net

- `tests/testthat/test-metab-02i-da-counts-table-characterization.R`
- `tests/testthat/test-metab-02j-da-results-class-characterization.R`
- `tests/testthat/test-metab-02k-da-numsig-barplot-characterization.R`
- `tests/testthat/test-metab-02l-da-volcano-plot-characterization.R`
- `tests/testthat/test-metab-02m-da-results-wide-format-characterization.R`
- `tests/testthat/test-metab-02n-da-results-long-format-characterization.R`
- `tests/testthat/test-metab-02o-da-interactive-volcano-characterization.R`
- `tests/testthat/test-metab-02p-da-analysis-characterization.R`
- `tests/testthat/test-metab-02q-da-analysis-helper-characterization.R`
- `tests/testthat/test-metab-02r-log-transform-assays-characterization.R`
- `tests/testthat/test-metab-02s-normalise-untransformed-data-characterization.R`
- `tests/testthat/test-metab-01d-qc-progress-helper-characterization.R`
- `tests/testthat/test-metab-01f-qc-plotting-helper-characterization.R`
- `tests/testthat/test-metab-01g-s4-plot-pca-characterization.R`
- `tests/testthat/test-metab-01j-s4-pearson-correlation-characterization.R`
- `tests/testthat/test-metab-01k-s4-plot-pearson-characterization.R`
- `tests/testthat/test-metab-01l-s4-calculate-pair-correlation-characterization.R`
- `tests/testthat/test-metab-01m-s4-resolve-duplicate-intensity-characterization.R`
- `tests/testthat/test-metab-01n-s4-find-duplicate-feature-ids-characterization.R`
- `tests/testthat/test-metab-02z-ruv-cancor-characterization.R`
- `tests/testthat/test-metab-03a-ruviii-c-varying-characterization.R`

## Notes

Manual bucket 0 metabolomics S4 stabilization target.

- This target previously had no active handover; this file now records the
  first metabolomics S4 stabilization stop point.
- Classification refresh on April 16, 2026 keeps
  `R/func_metab_s4_objects.R`
  at `4894` lines with `38` top-level functions while retaining
  `review` plus `direct-extraction-ready`.
- The first focused characterization gate now lives in
  `tests/testthat/test-metab-02i-da-counts-table-characterization.R`
  and freezes the current `getCountsTable()` contract by loading the extracted
  helper first when present, then falling back to
  `R/func_metab_s4_objects.R`
  and
  `R/func_metab_da.R`.
- Wave 1 manifest:
  `tools/refactor/manifest-metab-s4-wave1.yml`
  now applies the bounded metabolomics S4 counts-table helper checkpoint live
  into
  `R/func_metab_s4_da_results.R`
  and rewrites
  `R/func_metab_s4_objects.R`
  to remove
  `getCountsTable()`.
- The staged wave-1 review artifacts remain available at
  `tools/refactor/staging/R/func_metab_s4_da_results.R`
  and
  `tools/refactor/staging/wave1_metabolomics_s4_da_results/collate-metab-s4-wave1.txt`.
- The second focused characterization gate now lives in
  `tests/testthat/test-metab-02j-da-results-class-characterization.R`
  and freezes the
  `MetabolomicsDifferentialAbundanceResults`
  S4 class layout by loading the extracted class first when present, then
  falling back to
  `R/func_metab_s4_objects.R`.
- Wave 2 manifest:
  `tools/refactor/manifest-metab-s4-wave2.yml`
  now applies the bounded metabolomics S4 DA-results class checkpoint live
  into
  `R/func_metab_s4_da_results.R`
  and rewrites
  `R/func_metab_s4_objects.R`
  to remove
  `MetabolomicsDifferentialAbundanceResults`.
- The staged wave-2 review artifact now lives at
  `tools/refactor/staging/wave2_metabolomics_s4_da_results_class/R/func_metab_s4_da_results.R`.
- `DESCRIPTION` `Collate:` already loaded
  `R/func_metab_s4_da_results.R`
  ahead of
  `R/func_metab_s4_objects.R`,
  so wave 2 required no additional collate change.
- `DESCRIPTION` `Collate:` now loads
  `R/func_metab_s4_da_results.R`
  ahead of
  `R/func_metab_s4_objects.R`
  so the live helper layout matches the extracted shape.
- Post-apply checks passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-s4-wave1.yml`.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-s4-wave2.yml`.
- The third focused characterization gate now lives in
  `tests/testthat/test-metab-02k-da-numsig-barplot-characterization.R`
  and freezes the
  `plotNumSigDiffExpBarPlot()`
  list-method by loading the extracted helper file first when present, then
  falling back to
  `R/func_metab_s4_objects.R`.
- Wave 3 manifest:
  `tools/refactor/manifest-metab-s4-wave3.yml`
  now applies the bounded metabolomics S4 DA-results summary barplot checkpoint
  live into
  `R/func_metab_s4_da_results.R`
  and rewrites
  `R/func_metab_s4_objects.R`
  to remove
  `plotNumSigDiffExpBarPlot()`.
- The staged wave-3 review artifacts now live at
  `tools/refactor/staging/wave3_metabolomics_s4_da_results_numsig_barplot/R/func_metab_s4_da_results.R`
  and
  `tools/refactor/staging/wave3_metabolomics_s4_da_results_numsig_barplot/collate-metab-s4-wave3.txt`.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-s4-wave3.yml`.
- The focused gate reran green after apply via a direct `testthat::test_file()`
  invocation because this worktree does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- The fourth focused characterization gate now lives in
  `tests/testthat/test-metab-02l-da-volcano-plot-characterization.R`
  and freezes the
  `plotVolcanoS4()`
  list-method by loading the extracted helper file first when present, then
  falling back to
  `R/func_metab_s4_objects.R`.
- Wave 4 manifest:
  `tools/refactor/manifest-metab-s4-wave4.yml`
  now applies the bounded metabolomics S4 DA-results volcano plot checkpoint
  live into
  `R/func_metab_s4_da_results.R`
  and rewrites
  `R/func_metab_s4_objects.R`
  to remove
  `plotVolcanoS4()`.
- The staged wave-4 review artifacts now live at
  `tools/refactor/staging/wave4_metabolomics_s4_da_results_volcano_plot/R/func_metab_s4_da_results.R`
  and
  `tools/refactor/staging/wave4_metabolomics_s4_da_results_volcano_plot/collate-metab-s4-wave4.txt`.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-s4-wave4.yml`.
- Classification refresh on April 16, 2026 now measures
  `R/func_metab_s4_objects.R`
  at `4845` lines with `37` top-level functions while retaining
  `review` plus `direct-extraction-ready`.
- The fifth focused characterization gate now lives in
  `tests/testthat/test-metab-02m-da-results-wide-format-characterization.R`
  and freezes the
  `getDaResultsWideFormat()`
  list-method by loading the extracted helper file first when present, then
  falling back to
  `R/func_metab_s4_objects.R`.
- Wave 5 manifest:
  `tools/refactor/manifest-metab-s4-wave5.yml`
  now applies the bounded metabolomics S4 DA-results wide-format checkpoint
  live into
  `R/func_metab_s4_da_results.R`
  and rewrites
  `R/func_metab_s4_objects.R`
  to remove
  `getDaResultsWideFormat()`.
- The staged wave-5 review artifacts now live at
  `tools/refactor/staging/wave5_metabolomics_s4_da_results_wide_format/R/func_metab_s4_da_results.R`
  and
  `tools/refactor/staging/wave5_metabolomics_s4_da_results_wide_format/collate-metab-s4-wave5.txt`.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-s4-wave5.yml`.
- The sixth focused characterization gate now lives in
  `tests/testthat/test-metab-02n-da-results-long-format-characterization.R`
  and freezes the
  `getDaResultsLongFormat()`
  list-method by loading the extracted helper file first when present, then
  falling back to
  `R/func_metab_s4_objects.R`.
- Wave 6 manifest:
  `tools/refactor/manifest-metab-s4-wave6.yml`
  now applies the bounded metabolomics S4 DA-results long-format checkpoint
  live into
  `R/func_metab_s4_da_results.R`
  and rewrites
  `R/func_metab_s4_objects.R`
  to remove
  `getDaResultsLongFormat()`.
- The staged wave-6 review artifacts now live at
  `tools/refactor/staging/wave6_metabolomics_s4_da_results_long_format/R/func_metab_s4_da_results.R`
  and
  `tools/refactor/staging/wave6_metabolomics_s4_da_results_long_format/collate-metab-s4-wave6.txt`.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-s4-wave6.yml`.
- Classification refresh on April 16, 2026 now measures
  `R/func_metab_s4_objects.R`
  at `4805` lines with `36` top-level functions while retaining
  `review` plus `direct-extraction-ready`.
- The focused metabolomics S4 long-format gate reran green after the live
  apply via a direct `testthat::test_file()` invocation because this worktree
  does not include `renv/activate.R` for `tools/test_with_renv.R`.
- The seventh focused characterization gate now lives in
  `tests/testthat/test-metab-02o-da-interactive-volcano-characterization.R`
  and freezes the
  `plotInteractiveVolcano()`
  list-method with a zero-coefficient fallback run plus source-shape checks for
  the qvalue and Glimma hooks by loading the extracted helper file first when
  present, then falling back to
  `R/func_metab_s4_objects.R`.
- Wave 7 manifest:
  `tools/refactor/manifest-metab-s4-wave7.yml`
  now applies the bounded metabolomics S4 interactive-volcano checkpoint live
  into
  `R/func_metab_s4_da_results.R`
  and rewrites
  `R/func_metab_s4_objects.R`
  to remove
  `plotInteractiveVolcano()`.
- The staged wave-7 review artifacts now live at
  `tools/refactor/staging/wave7_metabolomics_s4_interactive_volcano/R/func_metab_s4_da_results.R`
  and
  `tools/refactor/staging/wave7_metabolomics_s4_interactive_volcano/collate-metab-s4-wave7.txt`.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-s4-wave7.yml`.
- Classification refresh on April 16, 2026 now measures
  `R/func_metab_s4_objects.R`
  at `4674` lines with `35` top-level functions while retaining
  `review` plus `direct-extraction-ready`.
- The focused metabolomics S4 interactive-volcano gate reran green after the
  live apply via a direct `testthat::test_file()` invocation because this
  worktree does not include `renv/activate.R` for `tools/test_with_renv.R`.
- The eighth focused characterization gate now lives in
  `tests/testthat/test-metab-02p-da-analysis-characterization.R`
  and freezes the
  `differentialAbundanceAnalysis()`
  list-method by loading the extracted helper file first when present, then
  falling back to
  `R/func_metab_s4_objects.R`.
- Wave 8 manifest:
  `tools/refactor/manifest-metab-s4-wave8.yml`
  now applies the bounded metabolomics S4 DA list-method checkpoint live
  into the dedicated helper file
  `R/func_metab_s4_da_methods.R`
  and rewrites
  `R/func_metab_s4_objects.R`
  to remove
  `differentialAbundanceAnalysis()`.
- The staged wave-8 review artifacts now live at
  `tools/refactor/staging/wave8_metabolomics_s4_da_method/R/func_metab_s4_da_methods.R`
  and
  `tools/refactor/staging/wave8_metabolomics_s4_da_method/collate-metab-s4-wave8.txt`.
- `DESCRIPTION` `Collate:` now loads
  `R/func_metab_s4_da_methods.R`
  between
  `R/func_metab_s4_da_results.R`
  and
  `R/func_metab_s4_objects.R`
  so the live package load order matches the extracted DA-method layout.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-s4-wave8.yml`.
- Classification refresh on April 16, 2026 now measures
  `R/func_metab_s4_objects.R`
  at `4625` lines with `34` top-level functions while retaining
  `review` plus `direct-extraction-ready`.
- The focused metabolomics S4 DA-analysis gate reran green after the live
  apply via a direct `testthat::test_file()` invocation because this
  worktree does not include `renv/activate.R` for `tools/test_with_renv.R`.
- The ninth focused characterization gate now lives in
  `tests/testthat/test-metab-02q-da-analysis-helper-characterization.R`
  and freezes the
  `differentialAbundanceAnalysisHelper()`
  S4 method with a numeric-leading group sanitization run plus source-shape
  checks by loading the extracted helper file first when present, then falling
  back to
  `R/func_metab_s4_objects.R`.
- Wave 9 manifest:
  `tools/refactor/manifest-metab-s4-wave9.yml`
  now applies the bounded metabolomics S4 DA helper-method checkpoint live
  into the existing helper file
  `R/func_metab_s4_da_methods.R`
  and rewrites
  `R/func_metab_s4_objects.R`
  to remove
  `differentialAbundanceAnalysisHelper()`.
- The staged wave-9 review artifacts now live at
  `tools/refactor/staging/wave9_metabolomics_s4_da_helper_method/R/func_metab_s4_da_methods.R`
  and
  `tools/refactor/staging/wave9_metabolomics_s4_da_helper_method/collate-metab-s4-wave9.txt`.
- `R/func_metab_s4_da_methods.R`
  already sat in the live `DESCRIPTION` collate order ahead of
  `R/func_metab_s4_objects.R`,
  so wave 9 required no additional collate change.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-s4-wave9.yml`.
- Classification refresh on April 16, 2026 now measures
  `R/func_metab_s4_objects.R`
  at `4459` lines with `33` top-level functions while retaining
  `review` plus `direct-extraction-ready`.
- The focused metabolomics S4 DA-helper gate reran green after the live
  apply via a direct `testthat::test_file()` invocation because this
  worktree does not include `renv/activate.R` for `tools/test_with_renv.R`.
- The metabolomics S4 DA-method family now lives in
  `R/func_metab_s4_da_methods.R`
  for both
  `differentialAbundanceAnalysis()`
  and
  `differentialAbundanceAnalysisHelper()`.
- The tenth focused characterization gate now lives in
  `tests/testthat/test-metab-02r-log-transform-assays-characterization.R`
  and freezes the
  `logTransformAssays()`
  S4 method with direct assay-transformation coverage plus source-shape checks
  by loading the extracted helper file first when present, then falling back to
  `R/func_metab_s4_objects.R`.
- Wave 10 manifest:
  `tools/refactor/manifest-metab-s4-wave10.yml`
  now applies the bounded metabolomics S4 log-transform normalization-method
  checkpoint live into the new helper file
  `R/func_metab_s4_norm_methods.R`
  and rewrites
  `R/func_metab_s4_objects.R`
  to remove
  `logTransformAssays()`.
- The staged wave-10 review artifacts now live at
  `tools/refactor/staging/wave10_metabolomics_s4_log_transform_method/R/func_metab_s4_norm_methods.R`
  and
  `tools/refactor/staging/wave10_metabolomics_s4_log_transform_method/collate-metab-s4-wave10.txt`.
- `DESCRIPTION` `Collate:` now loads
  `R/func_metab_s4_norm_methods.R`
  between
  `R/func_metab_s4_da_methods.R`
  and
  `R/func_metab_s4_objects.R`
  so the live package load order matches the extracted normalization-method
  layout.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-s4-wave10.yml`.
- Classification refresh on April 16, 2026 now measures
  `R/func_metab_s4_objects.R`
  at `4305` lines with `32` top-level functions while retaining
  `review` plus `direct-extraction-ready`.
- The focused metabolomics S4 log-transform gate reran green after the live
  apply via a direct `testthat::test_file()` invocation because this
  worktree does not include `renv/activate.R` for `tools/test_with_renv.R`.
- The metabolomics S4 normalization-method family now starts in
  `R/func_metab_s4_norm_methods.R`
  for
  `logTransformAssays()`.
- The eleventh focused characterization gate now lives in
  `tests/testthat/test-metab-02s-normalise-untransformed-data-characterization.R`
  and freezes the
  `normaliseUntransformedData()`
  S4 method with direct average-centered ITSD normalization coverage plus
  source-shape checks by loading the extracted helper file first when present,
  then falling back to
  `R/func_metab_s4_objects.R`.
- Wave 11 manifest:
  `tools/refactor/manifest-metab-s4-wave11.yml`
  now applies the bounded metabolomics S4 ITSD normalization-method checkpoint
  live into the existing helper file
  `R/func_metab_s4_norm_methods.R`
  and rewrites
  `R/func_metab_s4_objects.R`
  to remove
  `normaliseUntransformedData()`.
- The staged wave-11 review artifacts now live at
  `tools/refactor/staging/wave11_metabolomics_s4_normalise_untransformed_data_method/R/func_metab_s4_norm_methods.R`
  and
  `tools/refactor/staging/wave11_metabolomics_s4_normalise_untransformed_data_method/collate-metab-s4-wave11.txt`.
- `DESCRIPTION` `Collate:` already loaded
  `R/func_metab_s4_norm_methods.R`
  ahead of
  `R/func_metab_s4_objects.R`,
  so wave 11 required no additional collate change.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-s4-wave11.yml`.
- Classification refresh on April 16, 2026 now measures
  `R/func_metab_s4_objects.R`
  at `3922` lines with `29` top-level functions while retaining
  `review` plus `direct-extraction-ready`.
- The focused metabolomics S4 ITSD normalization gate reran green after the
  live apply via a direct `testthat::test_file()` invocation because this
  worktree does not include `renv/activate.R` for `tools/test_with_renv.R`.
- The previously extracted metabolomics S4 log-transform gate also reran green
  after the helper update for
  `R/func_metab_s4_norm_methods.R`.
- The metabolomics S4 normalization-method family now lives in
  `R/func_metab_s4_norm_methods.R`
  for both
  `logTransformAssays()`
  and
  `normaliseUntransformedData()`.
- The twelfth focused characterization gate now lives in
  `tests/testthat/test-metab-02t-metabolite-intensity-filtering-characterization.R`
  and freezes the
  `metaboliteIntensityFiltering()`
  S4 method with config-comment cleanup coverage plus source-shape checks by
  loading the extracted helper file first when present, then falling back to
  `R/func_metab_s4_objects.R`.
- Wave 12 manifest:
  `tools/refactor/manifest-metab-s4-wave12.yml`
  now applies the bounded metabolomics S4 QC-method checkpoint live into the
  new helper file
  `R/func_metab_s4_qc_methods.R`
  and rewrites
  `R/func_metab_s4_objects.R`
  to remove
  `metaboliteIntensityFiltering()`.
- The staged wave-12 review artifacts now live at
  `tools/refactor/staging/wave12_metabolomics_s4_metabolite_intensity_filtering_method/R/func_metab_s4_qc_methods.R`
  and
  `tools/refactor/staging/wave12_metabolomics_s4_metabolite_intensity_filtering_method/collate-metab-s4-wave12.txt`.
- `DESCRIPTION` `Collate:` now loads
  `R/func_metab_s4_qc_methods.R`
  ahead of
  `R/func_metab_s4_objects.R`
  so the live package load order matches the extracted metabolomics S4 QC-method
  layout.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-s4-wave12.yml`.
- Classification refresh on April 16, 2026 now measures
  `R/func_metab_s4_objects.R`
  at `3788` lines with `28` top-level functions while retaining
  `review` plus `direct-extraction-ready`.
- The focused metabolomics S4 intensity-filtering gate reran green after the
  live apply via a direct `testthat::test_file()` invocation because this
  worktree does not include `renv/activate.R` for `tools/test_with_renv.R`.
- The metabolomics S4 QC-method family now starts in
  `R/func_metab_s4_qc_methods.R`
  for
  `metaboliteIntensityFiltering()`.
- A thirteenth direct characterization gate now lives in
  `tests/testthat/test-metab-02u-resolve-duplicate-features-characterization.R`
  for the second bounded metabolomics S4 QC-method checkpoint:
  `resolveDuplicateFeatures()`.
- Wave 13 manifest:
  `tools/refactor/manifest-metab-s4-wave13.yml`
  now applies the bounded metabolomics S4 QC-method checkpoint live into the
  existing helper file
  `R/func_metab_s4_qc_methods.R`
  and rewrites
  `R/func_metab_s4_objects.R`
  to remove
  `resolveDuplicateFeatures()`.
- The staged wave-13 review artifacts now live at
  `tools/refactor/staging/wave13_metabolomics_s4_resolve_duplicate_features_method/R/func_metab_s4_qc_methods.R`
  and
  `tools/refactor/staging/wave13_metabolomics_s4_resolve_duplicate_features_method/collate-metab-s4-wave13.txt`.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-s4-wave13.yml`.
- Classification refresh on April 16, 2026 now measures
  `R/func_metab_s4_objects.R`
  at `3540` lines with `27` top-level functions while retaining
  `review` plus `direct-extraction-ready`.
- The focused metabolomics S4 duplicate-resolution gate reran green after the
  live apply via a direct `testthat::test_file()` invocation because this
  worktree does not include `renv/activate.R` for `tools/test_with_renv.R`.
- The metabolomics S4 QC-method family now lives in
  `R/func_metab_s4_qc_methods.R`
  for both
  `metaboliteIntensityFiltering()`
  and
  `resolveDuplicateFeatures()`.
- A fourteenth direct characterization gate now lives in
  `tests/testthat/test-metab-02v-sample-correlation-filter-characterization.R`
  for the third bounded metabolomics S4 QC-method checkpoint:
  `filterSamplesByMetaboliteCorrelationThreshold()`.
- Wave 14 manifest:
  `tools/refactor/manifest-metab-s4-wave14.yml`
  now applies the bounded metabolomics S4 QC-method checkpoint live into the
  existing helper file
  `R/func_metab_s4_qc_methods.R`
  and rewrites
  `R/func_metab_s4_objects.R`
  to remove
  `filterSamplesByMetaboliteCorrelationThreshold()`.
- The staged wave-14 review artifacts now live at
  `tools/refactor/staging/wave14_metabolomics_s4_sample_correlation_filter_method/R/func_metab_s4_qc_methods.R`
  and
  `tools/refactor/staging/wave14_metabolomics_s4_sample_correlation_filter_method/collate-metab-s4-wave14.txt`.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-s4-wave14.yml`.
- Classification refresh on April 16, 2026 now measures
  `R/func_metab_s4_objects.R`
  at `3439` lines with `26` top-level functions while retaining
  `review` plus `direct-extraction-ready`.
- The focused metabolomics S4 sample-correlation gate reran green after the
  live apply via a direct `testthat::test_file()` invocation because this
  worktree does not include `renv/activate.R` for `tools/test_with_renv.R`.
- The metabolomics S4 QC-method family now lives in
  `R/func_metab_s4_qc_methods.R`
  for
  `metaboliteIntensityFiltering()`,
  `resolveDuplicateFeatures()`,
  and
  `filterSamplesByMetaboliteCorrelationThreshold()`.
- The fifteenth bounded checkpoint now lives in
  `tools/refactor/manifest-metab-s4-wave15.yml`
  and stages the remaining `FilteringProgressMetabolomics` class from
  `R/func_metab_s4_objects.R` into the existing helper file
  `R/func_metab_qc_progress_helpers.R`.
- The staged wave-15 review artifacts now live at
  `tools/refactor/staging/wave15_metabolomics_s4_filtering_progress_class/R/func_metab_qc_progress_helpers.R`
  and
  `tools/refactor/staging/wave15_metabolomics_s4_filtering_progress_class/collate-metab-s4-wave15.txt`.
- Wave 15 applies live into `R/func_metab_qc_progress_helpers.R`, keeps the
  existing `getFilteringProgressMetabolomics()` and
  `updateFilteringProgressMetabolomics()` helpers in place, and removes the
  duplicate `FilteringProgressMetabolomics` class/accessor tail from
  `R/func_metab_s4_objects.R`.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-s4-wave15.yml`.
- Classification refresh on April 16, 2026 now measures
  `R/func_metab_s4_objects.R`
  at `3389` lines with `25` top-level functions while retaining
  `review` plus `direct-extraction-ready`.
- The focused metabolomics QC progress and plotting gates reran green after
  the live apply via direct `testthat::test_file()` invocations:
  `tests/testthat/test-metab-01d-qc-progress-helper-characterization.R`
  and
  `tests/testthat/test-metab-01f-qc-plotting-helper-characterization.R`.
- The `FilteringProgressMetabolomics` class plus its global access/update
  helpers now live together in
  `R/func_metab_qc_progress_helpers.R`, so the metabolomics S4 object wrapper
  no longer carries the duplicate QC progress tail.
- The sixteenth bounded checkpoint now lives in
  `tools/refactor/manifest-metab-s4-wave16.yml`
  and stages the `plotPca()` S4 method from
  `R/func_metab_s4_objects.R`
  into the new helper file
  `R/func_metab_s4_plotting_methods.R`.
- The staged wave-16 review artifacts now live at
  `tools/refactor/staging/wave16_metabolomics_s4_plot_pca_method/R/func_metab_s4_plotting_methods.R`
  and
  `tools/refactor/staging/wave16_metabolomics_s4_plot_pca_method/collate-metab-s4-wave16.txt`.
- Wave 16 applies live into `R/func_metab_s4_plotting_methods.R`, updates
  `DESCRIPTION` `Collate:` so the new helper loads before
  `R/func_metab_s4_objects.R`, and removes the extracted `plotPca()` method
  from the wrapper file.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-s4-wave16.yml`.
- Classification refresh on April 16, 2026 now measures
  `R/func_metab_s4_objects.R`
  at `3239` lines with `24` top-level functions while retaining
  `review` plus `direct-extraction-ready`.
- The focused metabolomics S4 plotting gates reran green after the live apply
  via direct `testthat::test_file()` invocations:
  `tests/testthat/test-metab-01g-s4-plot-pca-characterization.R`
  and
  `tests/testthat/test-metab-01f-qc-plotting-helper-characterization.R`.
- The metabolomics S4 plotting-method family now starts in
  `R/func_metab_s4_plotting_methods.R`
  for
  `plotPca()`.
- The seventeenth focused characterization gate now lives in
  `tests/testthat/test-metab-01h-s4-plot-rle-characterization.R`
  and freezes the
  `plotRle()`
  list-method by loading the extracted plotting helper file first when
  present, then falling back to
  `R/func_metab_s4_objects.R`.
- The seventeenth bounded checkpoint now lives in
  `tools/refactor/manifest-metab-s4-wave17.yml`
  and stages the `plotRle()` S4 method from
  `R/func_metab_s4_objects.R`
  into the existing helper file
  `R/func_metab_s4_plotting_methods.R`.
- The staged wave-17 review artifacts now live at
  `tools/refactor/staging/wave17_metabolomics_s4_plot_rle_method/R/func_metab_s4_plotting_methods.R`
  and
  `tools/refactor/staging/wave17_metabolomics_s4_plot_rle_method/collate-metab-s4-wave17.txt`.
- Wave 17 applies live into `R/func_metab_s4_plotting_methods.R`, keeps the
  existing `DESCRIPTION` `Collate:` order intact, and removes the extracted
  `plotRle()` method from the wrapper file.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-s4-wave17.yml`.
- Classification refresh on April 16, 2026 now measures
  `R/func_metab_s4_objects.R`
  at `3080` lines with `23` top-level functions while retaining
  `review` plus `direct-extraction-ready`.
- The focused metabolomics S4 plotting gates reran green after the live apply
  via direct `testthat::test_file()` invocations:
  `tests/testthat/test-metab-01h-s4-plot-rle-characterization.R`
  and
  `tests/testthat/test-metab-01f-qc-plotting-helper-characterization.R`.
- The metabolomics S4 plotting-method family now lives in
  `R/func_metab_s4_plotting_methods.R`
  for
  `plotPca()`
  and
  `plotRle()`.
- The eighteenth focused characterization gate now lives in
  `tests/testthat/test-metab-01i-s4-plot-density-characterization.R`
  and freezes both `plotDensity()` S4 methods by loading the extracted
  plotting helper file first when present, then falling back to
  `R/func_metab_s4_objects.R`.
- The eighteenth bounded checkpoint now lives in
  `tools/refactor/manifest-metab-s4-wave18.yml`
  and stages both `plotDensity()` S4 methods from
  `R/func_metab_s4_objects.R`
  into the existing helper file
  `R/func_metab_s4_plotting_methods.R`.
- The staged wave-18 review artifacts now live at
  `tools/refactor/staging/wave18_metabolomics_s4_plot_density_methods/R/func_metab_s4_plotting_methods.R`
  and
  `tools/refactor/collate-metab-s4-wave18.txt`.
- Wave 18 applies live into `R/func_metab_s4_plotting_methods.R`, keeps the
  existing `DESCRIPTION` `Collate:` order intact, and removes both extracted
  `plotDensity()` methods from the wrapper file.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-s4-wave18.yml`.
- Classification refresh on April 16, 2026 now measures
  `R/func_metab_s4_objects.R`
  at `2753` lines with `21` top-level functions while retaining
  `review` plus `direct-extraction-ready`.
- The focused metabolomics S4 plotting gates reran green after the live apply
  via direct `testthat::test_file()` invocations:
  `tests/testthat/test-metab-01f-qc-plotting-helper-characterization.R`,
  `tests/testthat/test-metab-01g-s4-plot-pca-characterization.R`,
  `tests/testthat/test-metab-01h-s4-plot-rle-characterization.R`,
  and
  `tests/testthat/test-metab-01i-s4-plot-density-characterization.R`.
- The metabolomics S4 plotting-method family now lives in
  `R/func_metab_s4_plotting_methods.R`
  for
  `plotPca()`,
  `plotRle()`,
  and
  `plotDensity()`.
- A nineteenth direct characterization gate now lives in
  `tests/testthat/test-metab-01j-s4-pearson-correlation-characterization.R`
  and freezes the current
  `pearsonCorForSamplePairs()`
  S4 method by loading the extracted QC-method helper first when present, then
  falling back to
  `R/func_metab_s4_objects.R`.
- Wave 19 manifest:
  `tools/refactor/manifest-metab-s4-wave19.yml`
  now applies the bounded metabolomics S4 Pearson-correlation QC-method
  checkpoint live into
  `R/func_metab_s4_qc_methods.R`
  and rewrites
  `R/func_metab_s4_objects.R`
  to remove
  `pearsonCorForSamplePairs()`.
- The staged wave-19 review artifacts now live at
  `tools/refactor/staging/wave19_metabolomics_s4_pearson_cor_for_sample_pairs/R/func_metab_s4_qc_methods.R`
  and
  `tools/refactor/collate-metab-s4-wave19.txt`.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-s4-wave19.yml`.
- Classification refresh on April 16, 2026 now measures
  `R/func_metab_s4_objects.R`
  at `2488` lines with `20` top-level functions while retaining
  `review` plus `direct-extraction-ready`.
- The focused metabolomics S4 Pearson-correlation gate reran green after the
  live apply via a direct `testthat::test_file()` invocation because this
  worktree does not include `renv/activate.R` for `tools/test_with_renv.R`.
- The metabolomics S4 QC-method helper file now also lives in
  `R/func_metab_s4_qc_methods.R`
  for
  `pearsonCorForSamplePairs()`.
- A twentieth direct characterization gate now lives in
  `tests/testthat/test-metab-01k-s4-plot-pearson-characterization.R`
  and freezes the current
  `plotPearson()`
  S4 method by loading the extracted QC-method helper first when present, then
  falling back to
  `R/func_metab_s4_objects.R`.
- Wave 20 manifest:
  `tools/refactor/manifest-metab-s4-wave20.yml`
  now applies the bounded metabolomics S4 Pearson-plot QC-method checkpoint
  live into
  `R/func_metab_s4_qc_methods.R`
  and rewrites
  `R/func_metab_s4_objects.R`
  to remove
  `plotPearson()`.
- The staged wave-20 review artifacts now live at
  `tools/refactor/staging/wave20_metabolomics_s4_plot_pearson/R/func_metab_s4_qc_methods.R`
  and
  `tools/refactor/collate-metab-s4-wave20.txt`.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-s4-wave20.yml`.
- Classification refresh on April 16, 2026 now measures
  `R/func_metab_s4_objects.R`
  at `2398` lines with `19` top-level functions while retaining
  `review` plus `direct-extraction-ready`.
- The focused metabolomics S4 Pearson-plot gate reran green after the
  live apply via a direct `testthat::test_file()` invocation because this
  worktree does not include `renv/activate.R` for `tools/test_with_renv.R`.
- `DESCRIPTION` `Collate:` already loads
  `R/func_metab_s4_qc_methods.R`
  ahead of
  `R/func_metab_s4_objects.R`
  so wave 20 required no additional collate change.
- The metabolomics S4 QC-method helper file now also lives in
  `R/func_metab_s4_qc_methods.R`
  for
  `pearsonCorForSamplePairs()`
  and
  `plotPearson()`.
- A twenty-first direct characterization gate now lives in
  `tests/testthat/test-metab-01l-s4-calculate-pair-correlation-characterization.R`
  and freezes the current
  `calculateMetabolitePairCorrelation()`
  helper by loading the extracted QC-method helper file first when present,
  then falling back to
  `R/func_metab_s4_objects.R`.
- Wave 21 manifest:
  `tools/refactor/manifest-metab-s4-wave21.yml`
  now applies the bounded metabolomics S4 pair-correlation helper checkpoint
  live into
  `R/func_metab_s4_qc_methods.R`
  and rewrites
  `R/func_metab_s4_objects.R`
  to remove
  `calculateMetabolitePairCorrelation()`.
- The staged wave-21 review artifacts now live at
  `tools/refactor/staging/wave21_metabolomics_s4_calculate_metabolite_pair_correlation/R/func_metab_s4_qc_methods.R`
  and
  `tools/refactor/collate-metab-s4-wave21.txt`.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-s4-wave21.yml`.
- Classification refresh on April 16, 2026 now measures
  `R/func_metab_s4_objects.R`
  at `2315` lines with `18` top-level functions while retaining
  `review` plus `direct-extraction-ready`.
- The focused metabolomics S4 pair-correlation helper gate reran green after
  the live apply via a direct `testthat::test_file()` invocation because this
  worktree does not include `renv/activate.R` for `tools/test_with_renv.R`.
- `DESCRIPTION` `Collate:` already loads
  `R/func_metab_s4_qc_methods.R`
  ahead of
  `R/func_metab_s4_objects.R`
  so wave 21 required no additional collate change.
- The metabolomics S4 QC-method family now also lives in
  `R/func_metab_s4_qc_methods.R`
  for
  `pearsonCorForSamplePairs()`,
  `plotPearson()`,
  and
  `calculateMetabolitePairCorrelation()`.
- A twenty-second direct characterization gate now lives in
  `tests/testthat/test-metab-01m-s4-resolve-duplicate-intensity-characterization.R`
  and freezes the current
  `resolveDuplicateFeaturesByIntensity()`
  helper by loading the extracted QC-method helper file first when present,
  then falling back to
  `R/func_metab_s4_objects.R`.
- The focused metabolomics S4 duplicate-intensity helper gate reran green via a
  direct `testthat::test_file()` invocation because this worktree does not
  include `renv/activate.R` for `tools/test_with_renv.R`.
- Wave 22 manifest:
  `tools/refactor/manifest-metab-s4-wave22.yml`
  now applies the bounded metabolomics S4 duplicate-intensity helper checkpoint
  live into
  `R/func_metab_s4_qc_methods.R`
  and rewrites
  `R/func_metab_s4_objects.R`
  to remove
  `resolveDuplicateFeaturesByIntensity()`.
- The staged wave-22 review artifacts now live at
  `tools/refactor/staging/wave22_metabolomics_s4_resolve_duplicate_features_by_intensity/R/func_metab_s4_qc_methods.R`
  and
  `tools/refactor/collate-metab-s4-wave22.txt`.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-s4-wave22.yml`.
- Classification refresh on April 16, 2026 now measures
  `R/func_metab_s4_objects.R`
  at `2247` lines with `17` top-level functions while retaining
  `review` plus `direct-extraction-ready`.
- The focused metabolomics S4 duplicate-intensity helper gate reran green after
  the live apply via a direct `testthat::test_file()` invocation because this
  worktree does not include `renv/activate.R` for `tools/test_with_renv.R`.
- `DESCRIPTION` `Collate:` already loads
  `R/func_metab_s4_qc_methods.R`
  ahead of
  `R/func_metab_s4_objects.R`
  so wave 22 required no additional collate change.
- The metabolomics S4 QC-method family now also lives in
  `R/func_metab_s4_qc_methods.R`
  for
  `pearsonCorForSamplePairs()`,
  `plotPearson()`,
  `calculateMetabolitePairCorrelation()`,
  and
  `resolveDuplicateFeaturesByIntensity()`.
- The twenty-third direct characterization gate now lives in
  `tests/testthat/test-metab-01n-s4-find-duplicate-feature-ids-characterization.R`
  and freezes
  `findMetabDuplicateFeatureIDs()`
  by loading the extracted QC-helper file first when present, then falling back
  to
  `R/func_metab_s4_objects.R`
  while preserving assay-specific duplicate counts, unnamed-assay fallback
  names, missing-column warnings, and the class guard.
- Wave 23 manifest:
  `tools/refactor/manifest-metab-s4-wave23.yml`
  now applies the bounded metabolomics S4 duplicate-ID helper checkpoint live
  from
  `R/func_metab_s4_objects.R`
  into
  `R/func_metab_s4_qc_methods.R`
  while rewriting
  `R/func_metab_s4_objects.R`
  to remove
  `findMetabDuplicateFeatureIDs()`.
- The staged wave-23 review artifacts now live at
  `tools/refactor/staging/wave23_metabolomics_s4_find_duplicate_feature_ids/R/func_metab_s4_qc_methods.R`
  and
  `tools/refactor/collate-metab-s4-wave23.txt`.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-s4-wave23.yml`.
- Classification refresh on April 16, 2026 now measures
  `R/func_metab_s4_objects.R`
  at `2165` lines with `16` top-level functions while retaining
  `review` plus `direct-extraction-ready`.
- The focused metabolomics S4 duplicate-ID helper gate reran green after
  the live apply via a direct `testthat::test_file()` invocation because this
  worktree does not include `renv/activate.R` for `tools/test_with_renv.R`.
- `DESCRIPTION` `Collate:` already loads
  `R/func_metab_s4_qc_methods.R`
  ahead of
  `R/func_metab_s4_objects.R`
  so wave 23 required no additional collate change.
- The metabolomics S4 QC-method family now also lives in
  `R/func_metab_s4_qc_methods.R`
  for
  `findMetabDuplicateFeatureIDs()`
  alongside
  `pearsonCorForSamplePairs()`,
  `plotPearson()`,
  `calculateMetabolitePairCorrelation()`,
  `resolveDuplicateFeaturesByIntensity()`,
  `metaboliteIntensityFiltering()`,
  and
  `filterSamplesByMetaboliteCorrelationThreshold()`.
- Next stop point for this handover is staging the bounded metabolomics S4
  metabolite-intensity filtering helper checkpoint for
  `metaboliteIntensityFilteringHelper()`
  from
  `R/func_metab_s4_objects.R`
  into
  `R/func_metab_s4_qc_methods.R`
  using
  `tests/testthat/test-metab-02t-metabolite-intensity-filtering-characterization.R`
  as the focused gate.
- The focused metabolomics S4 intensity-filtering characterization gate reran
  green before staging wave 24 via a direct `testthat::test_file()` invocation
  because this worktree does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Wave 24 manifest:
  `tools/refactor/manifest-metab-s4-wave24.yml`
  now stages the bounded metabolomics S4
  `metaboliteIntensityFilteringHelper()`
  checkpoint from
  `R/func_metab_s4_objects.R`
  into
  `R/func_metab_s4_qc_methods.R`
  without rewriting live sources.
- The staged wave-24 review artifacts now live at
  `tools/refactor/staging/wave24_metabolomics_s4_metabolite_intensity_filtering_helper/R/func_metab_s4_qc_methods.R`
  and
  `tools/refactor/collate-metab-s4-wave24.txt`.
- The focused metabolomics S4 intensity-filtering characterization gate reran
  green after staging wave 24 via a direct `testthat::test_file()` invocation
  because this worktree does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Wave 24 manifest:
  `tools/refactor/manifest-metab-s4-wave24.yml`
  now applies the bounded metabolomics S4
  `metaboliteIntensityFilteringHelper()`
  checkpoint live from
  `R/func_metab_s4_objects.R`
  into
  `R/func_metab_s4_qc_methods.R`
  while rewriting
  `R/func_metab_s4_objects.R`
  to remove
  `metaboliteIntensityFilteringHelper()`.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-s4-wave24.yml`.
- Classification refresh on April 16, 2026 now measures
  `R/func_metab_s4_objects.R`
  at `2114` lines with `15` top-level functions while retaining
  `review` plus `direct-extraction-ready`.
- The focused metabolomics S4 intensity-filtering characterization gate reran
  green after the live apply via a direct `testthat::test_file()` invocation
  because this worktree does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- `DESCRIPTION` `Collate:` already loads
  `R/func_metab_s4_qc_methods.R`
  ahead of
  `R/func_metab_s4_objects.R`
  so wave 24 required no additional collate change.
- The metabolomics S4 QC-method family now also lives in
  `R/func_metab_s4_qc_methods.R`
  for
  `metaboliteIntensityFilteringHelper()`
  alongside
  `pearsonCorForSamplePairs()`,
  `plotPearson()`,
  `calculateMetabolitePairCorrelation()`,
  `resolveDuplicateFeaturesByIntensity()`,
  `metaboliteIntensityFiltering()`,
  `filterSamplesByMetaboliteCorrelationThreshold()`,
  and
  `findMetabDuplicateFeatureIDs()`.
- A twenty-fifth direct characterization gate now lives in
  `tests/testthat/test-metab-02w-normalise-between-samples-characterization.R`
  and freezes the current metabolomics S4 between-sample normalization method
  by loading
  `R/func_metab_s4_norm_methods.R`
  first when present and otherwise falling back to
  `R/func_metab_s4_objects.R`.
- The focused metabolomics S4 between-sample normalization characterization
  gate reran green before staging via a direct `testthat::test_file()`
  invocation because this worktree does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- Wave 25 manifest:
  `tools/refactor/manifest-metab-s4-wave25.yml`
  now stages the bounded metabolomics S4
  `normaliseBetweenSamples()`
  checkpoint from
  `R/func_metab_s4_objects.R`
  into
  `R/func_metab_s4_norm_methods.R`.
- The staged wave-25 review artifacts now live at
  `tools/refactor/staging/wave25_metabolomics_s4_normalise_between_samples_method/R/func_metab_s4_norm_methods.R`
  and
  `tools/refactor/collate-metab-s4-wave25.txt`.
- Wave 25 manifest:
  `tools/refactor/manifest-metab-s4-wave25.yml`
  now applies the bounded metabolomics S4
  `normaliseBetweenSamples()`
  checkpoint live from
  `R/func_metab_s4_objects.R`
  into
  `R/func_metab_s4_norm_methods.R`
  while rewriting
  `R/func_metab_s4_objects.R`
  to remove
  `normaliseBetweenSamples()`.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-s4-wave25.yml`.
- Classification refresh on April 16, 2026 now measures
  `R/func_metab_s4_objects.R`
  at `1893` lines with `14` top-level functions while retaining
  `review` plus `direct-extraction-ready`.
- The focused metabolomics S4 between-sample normalization characterization
  gate reran green after the live apply via a direct `testthat::test_file()`
  invocation because this worktree does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- `DESCRIPTION` `Collate:` already loads
  `R/func_metab_s4_norm_methods.R`
  ahead of
  `R/func_metab_s4_objects.R`
  so wave 25 required no additional collate change.
- The metabolomics S4 normalization-method family now also lives in
  `R/func_metab_s4_norm_methods.R`
  for
  `normaliseBetweenSamples()`
  alongside
  `logTransformAssays()`
  and
  `normaliseUntransformedData()`.
- A twenty-sixth direct characterization gate now lives in
  `tests/testthat/test-metab-02x-clean-design-matrix-characterization.R`
  and freezes the current metabolomics S4 design-cleanup method by loading
  `R/func_metab_s4_norm_methods.R`
  first when present and otherwise falling back to
  `R/func_metab_s4_objects.R`.
- The focused metabolomics S4 design-cleanup characterization gate reran green
  before staging via a direct `testthat::test_file()` invocation because this
  worktree does not include `renv/activate.R` for `tools/test_with_renv.R`.
- Wave 26 manifest:
  `tools/refactor/manifest-metab-s4-wave26.yml`
  now stages the bounded metabolomics S4
  `cleanDesignMatrix()`
  checkpoint from
  `R/func_metab_s4_objects.R`
  into
  `R/func_metab_s4_norm_methods.R`.
- The staged wave-26 review artifacts now live at
  `tools/refactor/staging/wave26_metabolomics_s4_clean_design_matrix_method/R/func_metab_s4_norm_methods.R`
  and
  `tools/refactor/collate-metab-s4-wave26.txt`.
- Wave 26 manifest:
  `tools/refactor/manifest-metab-s4-wave26.yml`
  now applies the bounded metabolomics S4
  `cleanDesignMatrix()`
  checkpoint live from
  `R/func_metab_s4_objects.R`
  into
  `R/func_metab_s4_norm_methods.R`
  while rewriting
  `R/func_metab_s4_objects.R`
  to remove
  `cleanDesignMatrix()`.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-s4-wave26.yml`.
- Classification refresh on April 16, 2026 now measures
  `R/func_metab_s4_objects.R`
  at `1812` lines with `13` top-level functions while retaining
  `review` plus `direct-extraction-ready`.
- The focused metabolomics S4 design-cleanup characterization gate reran green
  after the live apply via a direct `testthat::test_file()` invocation because
  this worktree does not include `renv/activate.R` for `tools/test_with_renv.R`.
- The focused metabolomics S4 between-sample normalization characterization gate
  also reran green after the live apply via a direct `testthat::test_file()`
  invocation, confirming the adjacent normalization method still calls the
  extracted design-cleanup method correctly.
- `DESCRIPTION` `Collate:` already loads
  `R/func_metab_s4_norm_methods.R`
  ahead of
  `R/func_metab_s4_objects.R`
  so wave 26 required no additional collate change.
- The metabolomics S4 normalization-method family now also lives in
  `R/func_metab_s4_norm_methods.R`
  for
  `cleanDesignMatrix()`
  and
  `normaliseBetweenSamples()`
  alongside
  `logTransformAssays()`
  and
  `normaliseUntransformedData()`.
- A twenty-seventh direct characterization gate now lives in
  `tests/testthat/test-metab-02y-get-neg-ctrl-metab-anova-characterization.R`
  and freezes the current metabolomics S4 negative-control selector by loading
  `R/func_metab_s4_norm_methods.R`
  first when present and otherwise falling back to
  `R/func_metab_s4_objects.R`.
- The focused metabolomics S4 negative-control characterization gate reran
  green before staging via a direct `testthat::test_file()` invocation because
  this worktree does not include `renv/activate.R` for `tools/test_with_renv.R`.
- Wave 27 manifest:
  `tools/refactor/manifest-metab-s4-wave27.yml`
  now stages the bounded metabolomics S4
  `getNegCtrlMetabAnova()`
  checkpoint from
  `R/func_metab_s4_objects.R`
  into
  `R/func_metab_s4_norm_methods.R`.
- The staged wave-27 review artifacts now live at
  `tools/refactor/staging/wave27_metabolomics_s4_get_neg_ctrl_metab_anova/R/func_metab_s4_norm_methods.R`
  and
  `tools/refactor/collate-metab-s4-wave27.txt`.
- Wave 27 manifest:
  `tools/refactor/manifest-metab-s4-wave27.yml`
  now applies the bounded metabolomics S4
  `getNegCtrlMetabAnova()`
  checkpoint live from
  `R/func_metab_s4_objects.R`
  into
  `R/func_metab_s4_norm_methods.R`
  while rewriting
  `R/func_metab_s4_objects.R`
  to remove
  `getNegCtrlMetabAnova()`.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-s4-wave27.yml`.
- Classification refresh on April 16, 2026 now measures
  `R/func_metab_s4_objects.R`
  at `1441` lines with `9` top-level functions while retaining
  `review` plus `direct-extraction-ready`.
- The focused metabolomics S4 negative-control characterization gate reran
  green after the live apply via a direct `testthat::test_file()` invocation
  because this worktree does not include `renv/activate.R` for
  `tools/test_with_renv.R`.
- `DESCRIPTION` `Collate:` already loads
  `R/func_metab_s4_norm_methods.R`
  ahead of
  `R/func_metab_s4_objects.R`
  so wave 27 required no additional collate change.
- The metabolomics S4 normalization-method family now also lives in
  `R/func_metab_s4_norm_methods.R`
  for
  `getNegCtrlMetabAnova()`,
  `cleanDesignMatrix()`,
  and
  `normaliseBetweenSamples()`
  alongside
  `logTransformAssays()`
  and
  `normaliseUntransformedData()`.
- A twenty-eighth direct characterization gate now exists in
  `tests/testthat/test-metab-02z-ruv-cancor-characterization.R`
  and freezes the current metabolomics S4 RUV canonical-correlation method by
  loading
  `R/func_metab_s4_norm_methods.R`
  first when present and otherwise falling back to
  `R/func_metab_s4_objects.R`.
- The focused metabolomics S4 RUV canonical-correlation gate reran green
  before staging via a direct `testthat::test_file()` invocation because this
  worktree does not include `renv/activate.R` for `tools/test_with_renv.R`.
- Wave 28 manifest:
  `tools/refactor/manifest-metab-s4-wave28.yml`
  now stages the bounded metabolomics S4
  `ruvCancor()`
  checkpoint from
  `R/func_metab_s4_objects.R`
  into
  `R/func_metab_s4_norm_methods.R`.
- The staged wave-28 review artifacts now live at
  `tools/refactor/staging/wave28_metabolomics_s4_ruv_cancor/R/func_metab_s4_norm_methods.R`
  and
  `tools/refactor/collate-metab-s4-wave28.txt`.
- Wave 28 manifest:
  `tools/refactor/manifest-metab-s4-wave28.yml`
  now applies the bounded metabolomics S4
  `ruvCancor()`
  checkpoint live from
  `R/func_metab_s4_objects.R`
  into
  `R/func_metab_s4_norm_methods.R`
  while rewriting
  `R/func_metab_s4_objects.R`
  to remove
  `ruvCancor()`.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-s4-wave28.yml`.
- Classification refresh on April 16, 2026 now measures
  `R/func_metab_s4_objects.R`
  at `1088` lines with `8` top-level functions while retaining
  `direct-extraction-ready`.
- The focused metabolomics S4 RUV canonical-correlation gate reran green after
  the live apply via a direct `testthat::test_file()` invocation because this
  worktree does not include `renv/activate.R` for `tools/test_with_renv.R`.
- `DESCRIPTION` `Collate:` already loads
  `R/func_metab_s4_norm_methods.R`
  ahead of
  `R/func_metab_s4_objects.R`
  so wave 28 required no additional collate change.
- The metabolomics S4 normalization-method family now also lives in
  `R/func_metab_s4_norm_methods.R`
  for
  `ruvCancor()`,
  `getNegCtrlMetabAnova()`,
  `cleanDesignMatrix()`,
  and
  `normaliseBetweenSamples()`
  alongside
  `logTransformAssays()`
  and
  `normaliseUntransformedData()`.
- A twenty-ninth direct characterization gate now exists in
  `tests/testthat/test-metab-03a-ruviii-c-varying-characterization.R`
  and freezes the current metabolomics S4 RUV-III varying-k method by loading
  `R/func_metab_s4_norm_methods.R`
  first when present and otherwise falling back to
  `R/func_metab_s4_objects.R`.
- The focused metabolomics S4 RUV-III varying-k gate reran green before
  staging via a direct `testthat::test_file()` invocation because this
  worktree does not include `renv/activate.R` for `tools/test_with_renv.R`.
- Wave 29 manifest:
  `tools/refactor/manifest-metab-s4-wave29.yml`
  now stages the bounded metabolomics S4
  `ruvIII_C_Varying()`
  checkpoint from
  `R/func_metab_s4_objects.R`
  into
  `R/func_metab_s4_norm_methods.R`.
- The staged wave-29 review artifacts now live at
  `tools/refactor/staging/wave29_metabolomics_s4_ruviii_c_varying/R/func_metab_s4_norm_methods.R`
  and
  `tools/refactor/collate-metab-s4-wave29.txt`.
- Wave 29 manifest:
  `tools/refactor/manifest-metab-s4-wave29.yml`
  now applies the bounded metabolomics S4
  `ruvIII_C_Varying()`
  checkpoint live from
  `R/func_metab_s4_objects.R`
  into
  `R/func_metab_s4_norm_methods.R`
  while rewriting
  `R/func_metab_s4_objects.R`
  to remove
  `ruvIII_C_Varying()`.
- Post-apply checks also passed through
  `tools/refactor/check_wave_apply.R --manifest tools/refactor/manifest-metab-s4-wave29.yml`.
- Classification refresh on April 16, 2026 now measures
  `R/func_metab_s4_objects.R`
  at `511` lines with `2` top-level functions while retaining
  `direct-extraction-ready`.
- The focused metabolomics S4 RUV-III varying-k gate reran green after the
  live apply via a direct `testthat::test_file()` invocation because this
  worktree does not include `renv/activate.R` for `tools/test_with_renv.R`.
- `DESCRIPTION` `Collate:` already loads
  `R/func_metab_s4_norm_methods.R`
  ahead of
  `R/func_metab_s4_objects.R`
  so wave 29 required no additional collate change.
- The metabolomics S4 normalization-method family now also lives in
  `R/func_metab_s4_norm_methods.R`
  for
  `ruvIII_C_Varying()`,
  `ruvCancor()`,
  `getNegCtrlMetabAnova()`,
  `cleanDesignMatrix()`,
  and
  `normaliseBetweenSamples()`
  alongside
  `logTransformAssays()`
  and
  `normaliseUntransformedData()`.
- No additional live method blocks remain in
  `R/func_metab_s4_objects.R`;
  treat the file as a stabilized breadcrumb/generics shell for this target and
  move later bucket work to a different backlog target if further bucket work
  is needed.
