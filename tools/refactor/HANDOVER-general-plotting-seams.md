# General Plotting Seam Map

## Goal

Freeze bounded extraction slices for `func_general_plotting.R` while staged
exact-source waves drain the file without changing the public plotting surface.

## Current Position In The Flow

- target: `/home/doktersmol/Documents/MultiScholaR/R/func_general_plotting.R`
- classification: `complete (archival)`
- checkpoint reached: `general_plotting_reporting_dispatch_apply_wave7`
- next step: `Manual target is complete; keep this handover as the archival seam record for bucket 13 and do not reopen R/func_general_plotting.R unless a regression specifically lands there.`

## Existing Safety Net

- `tests/testthat/test-prot-04-design.R`
- `tests/testthat/test-prot-08-volcano.R`
- `tests/testthat/test-prot-09-heatmap.R`

## Notes

Manual target bucket 13 for shared plotting stabilization.

- This target previously had no dedicated handover; this file now records the
  first safe stop point for `R/func_general_plotting.R`.
- April 18, 2026 classification stayed `direct-extraction-ready` with `41`
  top-level functions and a `94`-line maximum helper span.
- April 18, 2026
  [tools/refactor/manifest-general-plotting-pca-rle-wave1.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-plotting-pca-rle-wave1.yml:1)
  was authored and verified cleanly against the live source tree.
- The first bounded design/QC helper wave staged cleanly into
  [tools/refactor/staging/general-plotting-pca-rle-wave1/R/func_general_plotting_pca_rle_helpers.R](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-plotting-pca-rle-wave1/R/func_general_plotting_pca_rle_helpers.R:1)
  for:
  - `plotPcaHelper()`
  - `plotPcaListHelper()`
  - `plotPcaGgpairs()`
  - `plotRleHelper()`
  - `getMaxMinBoxplot()`
  - `rlePcaPlotList()`
- The staged helper file currently measures `548` lines, which is above the
  ideal band but within the acceptable `501-800` range for a first bounded
  general-plotting slice.
- The staged collate preview is recorded in
  [tools/refactor/staging/general-plotting-pca-rle-wave1/collate-general-plotting-pca-rle-wave1.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-plotting-pca-rle-wave1/collate-general-plotting-pca-rle-wave1.txt:1).
- The focused gate reran green after staging through
  `Rscript tools/test_with_renv.R tests/testthat/test-prot-04-design.R`
  with `1572` passing expectations and `1` skipped Git LFS-backed snapshot
  check.
- April 18, 2026 the reviewed PCA/RLE wave was applied live from
  [tools/refactor/manifest-general-plotting-pca-rle-wave1.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-plotting-pca-rle-wave1.yml:1)
  into:
  - [R/func_general_plotting_pca_rle_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_plotting_pca_rle_helpers.R:1)
  - [R/func_general_plotting.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_plotting.R:1)
- The live apply materialized the first bounded PCA/RLE helper slice out of
  [R/func_general_plotting.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_plotting.R:1)
  for:
  - `plotPcaHelper()`
  - `plotPcaListHelper()`
  - `plotPcaGgpairs()`
  - `plotRleHelper()`
  - `getMaxMinBoxplot()`
  - `rlePcaPlotList()`
- The emitted collate artifact now exists at
  [tools/refactor/collate-general-plotting-pca-rle-wave1.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/collate-general-plotting-pca-rle-wave1.txt:1),
  and `DESCRIPTION` now collates
  `R/func_general_plotting_pca_rle_helpers.R` immediately before
  `R/func_general_plotting.R`.
- The focused gate stayed green before apply and reran green after the live
  apply through `Rscript tools/test_with_renv.R tests/testthat/test-prot-04-design.R`
  with the same `1572` passing expectations and `1` skipped Git LFS-backed
  snapshot check.
- Post-apply classification remains `direct-extraction-ready`, and
  [R/func_general_plotting.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_plotting.R:1)
  now reports `2296` lines, `32` top-level functions, and a `94`-line maximum
  helper span.
- April 18, 2026
  [tools/refactor/manifest-general-plotting-qc-support-wave2.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-plotting-qc-support-wave2.yml:1)
  was authored and verified cleanly against the current live source tree.
- The second bounded general-plotting follow-up wave staged cleanly into
  [tools/refactor/staging/general-plotting-qc-support-wave2/R/func_general_plotting_qc_support_helpers.R](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-plotting-qc-support-wave2/R/func_general_plotting_qc_support_helpers.R:1)
  for:
  - `plotDensityOfProteinIntensityPerSample()`
  - `plotPercentSamplesVsProteinQuantified()`
  - `plotNumMissingValues()`
  - `plotNumOfValues()`
  - `plotNumOfValuesNoLog()`
- The staged QC support helper file measures `167` lines, which stays inside
  the ideal sub-`500` line band for a bounded follow-up slice.
- The staged collate preview is recorded in
  [tools/refactor/staging/general-plotting-qc-support-wave2/collate-general-plotting-qc-support-wave2.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-plotting-qc-support-wave2/collate-general-plotting-qc-support-wave2.txt:1).
- The focused gate reran green after staging through
  `Rscript tools/test_with_renv.R tests/testthat/test-prot-04-design.R`
  with `1572` passing expectations and `1` skipped Git LFS-backed snapshot
  check.
- This checkpoint stops at a staged wave only, so live
  [R/func_general_plotting.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_plotting.R:1)
  remains at `2296` lines, `32` top-level functions, and a `94`-line maximum
  helper span.
- This stop point now keeps correlation heatmaps, protein heatmaps, volcano,
  enrichment, and file-writer helpers live in `R/func_general_plotting.R`
  pending review and any later live apply of the staged QC support wave.
- April 18, 2026 the reviewed QC support wave was applied live from
  [tools/refactor/manifest-general-plotting-qc-support-wave2.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-plotting-qc-support-wave2.yml:1)
  into:
  - [R/func_general_plotting_qc_support_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_plotting_qc_support_helpers.R:1)
  - [R/func_general_plotting.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_plotting.R:1)
- The live apply materialized the bounded QC support helper slice for:
  - `plotDensityOfProteinIntensityPerSample()`
  - `plotPercentSamplesVsProteinQuantified()`
  - `plotNumMissingValues()`
  - `plotNumOfValues()`
  - `plotNumOfValuesNoLog()`
- The emitted collate artifact now exists at
  [tools/refactor/collate-general-plotting-qc-support-wave2.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/collate-general-plotting-qc-support-wave2.txt:1),
  and `DESCRIPTION` now collates
  `R/func_general_plotting_qc_support_helpers.R` immediately before
  `R/func_general_plotting.R`.
- The focused gate stayed green before apply and reran green after the live
  apply through `Rscript tools/test_with_renv.R tests/testthat/test-prot-04-design.R`
  with the same `1572` passing expectations and `1` skipped Git LFS-backed
  snapshot check.
- Post-apply classification remains `direct-extraction-ready`, and
  [R/func_general_plotting.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_plotting.R:1)
  now reports `2134` lines, `26` top-level functions, and a `94`-line maximum
  helper span.
- This stop point now keeps heatmap, volcano, palette/theme, and save helpers
  live in `R/func_general_plotting.R`; the next bounded wave should stage the
  heatmap/save-support cluster without rewriting live sources yet.
- April 18, 2026
  [tools/refactor/manifest-general-plotting-heatmap-support-wave3.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-plotting-heatmap-support-wave3.yml:1)
  was authored and verified cleanly against the current live source tree.
- The third bounded general-plotting follow-up wave staged cleanly into
  [tools/refactor/staging/general-plotting-heatmap-support-wave3/R/func_general_plotting_heatmap_support_helpers.R](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-plotting-heatmap-support-wave3/R/func_general_plotting_heatmap_support_helpers.R:1)
  for:
  - `getSamplesCorrelationHeatMap()`
  - `getProteinsHeatMap()`
  - `save_heatmap_products()`
- The staged heatmap/save-support helper file measures `320` lines, which stays
  inside the ideal sub-`500` line band for a bounded follow-up slice.
- The staged collate preview is recorded in
  [tools/refactor/staging/general-plotting-heatmap-support-wave3/collate-general-plotting-heatmap-support-wave3.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-plotting-heatmap-support-wave3/collate-general-plotting-heatmap-support-wave3.txt:1).
- The focused gate reran green after staging through
  `Rscript tools/test_with_renv.R tests/testthat/test-prot-04-design.R`
  with `1572` passing expectations and `1` skipped Git LFS-backed snapshot
  check.
- This checkpoint stops at a staged wave only, so live
  [R/func_general_plotting.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_plotting.R:1)
  remains at `2134` lines, `26` top-level functions, and a `94`-line maximum
  helper span.
- This stop point now keeps volcano/Glimma, enrichment, and palette/theme
  helpers live in `R/func_general_plotting.R` pending review and any later
  live apply of the staged heatmap/save-support wave.
- April 18, 2026 the reviewed heatmap/save-support wave was applied live from
  [tools/refactor/manifest-general-plotting-heatmap-support-wave3.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-plotting-heatmap-support-wave3.yml:1)
  into:
  - [R/func_general_plotting_heatmap_support_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_plotting_heatmap_support_helpers.R:1)
  - [R/func_general_plotting.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_plotting.R:1)
- The live apply materialized the bounded heatmap/save-support helper slice
  for:
  - `getSamplesCorrelationHeatMap()`
  - `getProteinsHeatMap()`
  - `save_heatmap_products()`
- The emitted collate artifact now exists at
  [tools/refactor/collate-general-plotting-heatmap-support-wave3.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/collate-general-plotting-heatmap-support-wave3.txt:1),
  and `DESCRIPTION` now collates
  `R/func_general_plotting_qc_support_helpers.R` and
  `R/func_general_plotting_heatmap_support_helpers.R` immediately before
  `R/func_general_plotting.R`.
- The focused gate stayed green before apply and reran green after the live
  apply through `Rscript tools/test_with_renv.R tests/testthat/test-prot-04-design.R`
  with the same `1572` passing expectations and `1` skipped Git LFS-backed
  snapshot check.
- Post-apply classification remains `direct-extraction-ready`, and
  [R/func_general_plotting.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_plotting.R:1)
  now reports `1817` lines, `23` top-level functions, and a `52`-line maximum
  helper span.
- This stop point now keeps the volcano/Glimma helpers, remaining QC utility
  plots, and palette/theme helpers live in `R/func_general_plotting.R`; the
  next bounded checkpoint should stage a dedicated volcano/Glimma follow-up
  wave without rewriting live sources yet.
- April 18, 2026
  [tools/refactor/manifest-general-plotting-volcano-glimma-wave4.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-plotting-volcano-glimma-wave4.yml:1)
  was authored and verified cleanly against the current live source tree.
- The fourth bounded general-plotting follow-up wave staged cleanly into
  [tools/refactor/staging/general-plotting-volcano-glimma-wave4/R/func_general_plotting_volcano_glimma_helpers.R](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-plotting-volcano-glimma-wave4/R/func_general_plotting_volcano_glimma_helpers.R:1)
  for:
  - `plotOneVolcano()`
  - `plotOneVolcanoNoVerticalLines()`
  - `printOneVolcanoPlotWithProteinLabel()`
  - `getGlimmaVolcanoProteomics()`
  - `getGlimmaVolcanoProteomicsWidget()`
  - `getGlimmaVolcanoPhosphoproteomics()`
  - `plotVolcano()`
- The staged volcano/Glimma helper file measures `658` lines, which lands in
  the acceptable `501-800` band for a dense but bounded shared plotting slice.
- The staged collate preview is recorded in
  [tools/refactor/staging/general-plotting-volcano-glimma-wave4/collate-general-plotting-volcano-glimma-wave4.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-plotting-volcano-glimma-wave4/collate-general-plotting-volcano-glimma-wave4.txt:1).
- The focused gate reran green after staging through
  `Rscript tools/test_with_renv.R tests/testthat/test-prot-04-design.R`
  with `1572` passing expectations and `1` skipped Git LFS-backed snapshot
  check.
- Supplemental volcano snapshot coverage was rerun through
  `Rscript tools/test_with_renv.R tests/testthat/test-prot-08-volcano.R`, but
  the `cp08_volcano_input.rds` fixture is unreadable in this workspace
  (`readRDS()` unknown input format), so that safety-net test remains
  environment-blocked rather than code-regressed.
- This checkpoint stops at a staged wave only, so live
  [R/func_general_plotting.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_plotting.R:1)
  remains at `1817` lines, `23` top-level functions, and a `52`-line maximum
  helper span.
- This stop point now keeps the staged volcano/Glimma helper file pending
  review while the remaining QC utility plots and palette/theme helpers stay
  live in `R/func_general_plotting.R`; the next bounded checkpoint should
  review and apply the staged volcano/Glimma wave without opening a second
  seam.
- April 18, 2026 the reviewed volcano/Glimma wave was applied live from
  [tools/refactor/manifest-general-plotting-volcano-glimma-wave4.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-plotting-volcano-glimma-wave4.yml:1)
  into:
  - [R/func_general_plotting_volcano_glimma_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_plotting_volcano_glimma_helpers.R:1)
  - [R/func_general_plotting.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_plotting.R:1)
- The live apply materialized the bounded volcano/Glimma helper slice for:
  - `plotOneVolcano()`
  - `plotOneVolcanoNoVerticalLines()`
  - `printOneVolcanoPlotWithProteinLabel()`
  - `getGlimmaVolcanoProteomics()`
  - `getGlimmaVolcanoProteomicsWidget()`
  - `getGlimmaVolcanoPhosphoproteomics()`
  - `plotVolcano()`
- The emitted collate artifact now exists at
  [tools/refactor/collate-general-plotting-volcano-glimma-wave4.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/collate-general-plotting-volcano-glimma-wave4.txt:1),
  and `DESCRIPTION` now collates
  `R/func_general_plotting_qc_support_helpers.R`,
  `R/func_general_plotting_heatmap_support_helpers.R`, and
  `R/func_general_plotting_volcano_glimma_helpers.R` immediately before
  `R/func_general_plotting.R`.
- The focused gate stayed green after the live apply through
  `Rscript tools/test_with_renv.R tests/testthat/test-prot-04-design.R`
  with `1572` passing expectations and `1` skipped Git LFS-backed snapshot
  check.
- Post-apply classification remains `direct-extraction-ready`, and
  [R/func_general_plotting.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_plotting.R:1)
  now reports `1166` lines, `16` top-level functions, and a `52`-line maximum
  helper span.
- This stop point now keeps the remaining QC utility helpers, palette/theme
  helpers, reporting helpers, and `plotPcaDispatch()` live in
  `R/func_general_plotting.R`; the next bounded checkpoint should stage a
  dedicated palette/theme helper wave for:
  - `getCategoricalColourPalette()`
  - `getOneContinousPalette()`
  - `getContinousColourRules()`
  - `getCategoricalAndContinuousColourRules()`
  - `apafTheme()`
  - `get_color_palette()`
  without opening a second live-apply step in the same iteration.
- April 18, 2026
  [tools/refactor/manifest-general-plotting-palette-theme-wave5.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-plotting-palette-theme-wave5.yml:1)
  was authored and verified cleanly against the current live source tree.
- The fifth bounded general-plotting follow-up wave staged cleanly into
  [tools/refactor/staging/general-plotting-palette-theme-wave5/R/func_general_plotting_palette_theme_helpers.R](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-plotting-palette-theme-wave5/R/func_general_plotting_palette_theme_helpers.R:1)
  for:
  - `getCategoricalColourPalette()`
  - `getOneContinousPalette()`
  - `getContinousColourRules()`
  - `getCategoricalAndContinuousColourRules()`
  - `apafTheme()`
  - `get_color_palette()`
- The staged palette/theme helper file measures `210` lines, which stays
  inside the ideal sub-`500` line band for a bounded follow-up slice.
- The staged collate preview is recorded in
  [tools/refactor/staging/general-plotting-palette-theme-wave5/collate-general-plotting-palette-theme-wave5.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-plotting-palette-theme-wave5/collate-general-plotting-palette-theme-wave5.txt:1).
- The focused gate reran green after staging through
  `Rscript tools/test_with_renv.R tests/testthat/test-prot-04-design.R`
  with `1572` passing expectations and `1` skipped Git LFS-backed snapshot
  check.
- This checkpoint stops at a staged wave only, so live
  [R/func_general_plotting.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_plotting.R:1)
  remains at `1166` lines, `16` top-level functions, and a `52`-line maximum
  helper span.
- This stop point now keeps the staged palette/theme helper file pending
  review while the remaining QC utility, reporting, and dispatch helpers stay
  live in `R/func_general_plotting.R`; the next bounded checkpoint should
  review and apply the staged palette/theme wave without opening a second seam
  in the same iteration.
- April 18, 2026 the reviewed palette/theme wave was applied live from
  [tools/refactor/manifest-general-plotting-palette-theme-wave5.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-plotting-palette-theme-wave5.yml:1)
  into:
  - [R/func_general_plotting_palette_theme_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_plotting_palette_theme_helpers.R:1)
  - [R/func_general_plotting.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_plotting.R:1)
- The live apply materialized the bounded palette/theme helper slice for:
  - `getCategoricalColourPalette()`
  - `getOneContinousPalette()`
  - `getContinousColourRules()`
  - `getCategoricalAndContinuousColourRules()`
  - `apafTheme()`
  - `get_color_palette()`
- The emitted collate artifact now exists at
  [tools/refactor/collate-general-plotting-palette-theme-wave5.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/collate-general-plotting-palette-theme-wave5.txt:1),
  and `DESCRIPTION` now collates
  `R/func_general_plotting_palette_theme_helpers.R` immediately before
  `R/func_general_plotting.R`.
- The focused gate reran green after the live apply through
  `Rscript tools/test_with_renv.R tests/testthat/test-prot-04-design.R`
  with `1572` passing expectations and `1` skipped Git LFS-backed snapshot
  check.
- Post-apply classification remains `direct-extraction-ready`, and
  [R/func_general_plotting.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_plotting.R:1)
  now reports `962` lines, `10` top-level functions, and a `52`-line maximum
  helper span.
- This stop point now keeps the remaining QC utility helpers, reporting
  helpers, and `plotPcaDispatch()` live in `R/func_general_plotting.R`; the
  next bounded checkpoint should stage a dedicated QC utility helper wave for:
  - `plotPeptidesProteinsCountsPerSampleHelper()`
  - `plotHistogramOfPercentMissingPerIndvidual()`
  - `getOneRlePlotData()`
  - `plotRleQc()`
  - `compareUmapComponentsPairs()`
  - `umap_factor_plot()`
  while keeping `printPValuesDistribution()`, `gg_save_logging()`,
  `summarizeQCPlot()`, and `plotPcaDispatch()` live until that staged review
  lands.
- April 18, 2026
  [tools/refactor/manifest-general-plotting-qc-utility-wave6.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-plotting-qc-utility-wave6.yml:1)
  was authored and verified cleanly against the current live source tree.
- The sixth bounded general-plotting follow-up wave staged cleanly into
  [tools/refactor/staging/general-plotting-qc-utility-wave6/R/func_general_plotting_qc_utility_helpers.R](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-plotting-qc-utility-wave6/R/func_general_plotting_qc_utility_helpers.R:1)
  for:
  - `plotPeptidesProteinsCountsPerSampleHelper()`
  - `plotHistogramOfPercentMissingPerIndvidual()`
  - `getOneRlePlotData()`
  - `plotRleQc()`
  - `compareUmapComponentsPairs()`
  - `umap_factor_plot()`
- The staged QC utility helper file measures `186` lines, which stays inside
  the ideal sub-`500` line band for a bounded follow-up slice.
- The staged collate preview is recorded in
  [tools/refactor/staging/general-plotting-qc-utility-wave6/collate-general-plotting-qc-utility-wave6.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-plotting-qc-utility-wave6/collate-general-plotting-qc-utility-wave6.txt:1).
- The focused gate reran green after staging through
  `Rscript tools/test_with_renv.R tests/testthat/test-prot-04-design.R`
  with `1572` passing expectations and `1` skipped Git LFS-backed snapshot
  check.
- This checkpoint stops at a staged wave only, so live
  [R/func_general_plotting.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_plotting.R:1)
  remains `962` lines with `10` top-level functions, a `52`-line maximum
  helper span, and `direct-extraction-ready` classification.
- This stop point now keeps the staged QC utility helper file pending review
  while `printPValuesDistribution()`, `gg_save_logging()`,
  `summarizeQCPlot()`, and `plotPcaDispatch()` stay live in
  `R/func_general_plotting.R`; the next bounded checkpoint should review and
  apply the staged QC utility wave without opening a second seam in the same
  iteration.
- April 18, 2026 the reviewed QC utility wave was applied live from
  [tools/refactor/manifest-general-plotting-qc-utility-wave6.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-plotting-qc-utility-wave6.yml:1)
  into:
  - [R/func_general_plotting_qc_utility_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_plotting_qc_utility_helpers.R:1)
  - [R/func_general_plotting.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_plotting.R:1)
- The live apply materialized the bounded QC utility helper slice for:
  - `plotPeptidesProteinsCountsPerSampleHelper()`
  - `plotHistogramOfPercentMissingPerIndvidual()`
  - `getOneRlePlotData()`
  - `plotRleQc()`
  - `compareUmapComponentsPairs()`
  - `umap_factor_plot()`
- The emitted collate artifact now exists at
  [tools/refactor/collate-general-plotting-qc-utility-wave6.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/collate-general-plotting-qc-utility-wave6.txt:1),
  and `DESCRIPTION` now collates
  `R/func_general_plotting_qc_utility_helpers.R` immediately before
  `R/func_general_plotting.R`.
- The focused gate stayed green before apply and reran green after the live
  apply through `Rscript tools/test_with_renv.R tests/testthat/test-prot-04-design.R`
  with `1572` passing expectations and `1` skipped Git LFS-backed snapshot
  check.
- Post-apply classification remains `direct-extraction-ready`, and
  [R/func_general_plotting.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_plotting.R:1)
  now reports `782` lines, `4` top-level functions, and a `52`-line maximum
  helper span.
- This stop point now keeps only the reporting helpers and
  `plotPcaDispatch()` live in `R/func_general_plotting.R`; the next bounded
  checkpoint should stage a final reporting/dispatch helper wave for:
  - `printPValuesDistribution()`
  - `gg_save_logging()`
  - `summarizeQCPlot()`
  - `plotPcaDispatch()`
- April 18, 2026
  [tools/refactor/manifest-general-plotting-reporting-dispatch-wave7.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-plotting-reporting-dispatch-wave7.yml:1)
  was authored and verified cleanly against the current live source tree.
- The seventh bounded general-plotting follow-up wave staged cleanly into
  [tools/refactor/staging/general-plotting-reporting-dispatch-wave7/R/func_general_plotting_reporting_dispatch_helpers.R](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-plotting-reporting-dispatch-wave7/R/func_general_plotting_reporting_dispatch_helpers.R:1)
  for:
  - `printPValuesDistribution()`
  - `gg_save_logging()`
  - `summarizeQCPlot()`
  - `plotPcaDispatch()`
- The staged reporting/dispatch helper file measures `149` lines, which stays
  inside the ideal sub-`500` line band for a bounded final slice.
- The staged collate preview is recorded in
  [tools/refactor/staging/general-plotting-reporting-dispatch-wave7/collate-general-plotting-reporting-dispatch-wave7.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/staging/general-plotting-reporting-dispatch-wave7/collate-general-plotting-reporting-dispatch-wave7.txt:1).
- The focused gate reran green after staging through
  `Rscript tools/test_with_renv.R tests/testthat/test-prot-04-design.R`
  with `1572` passing expectations and `1` skipped Git LFS-backed snapshot
  check.
- This checkpoint stops at a staged wave only, so live
  [R/func_general_plotting.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_plotting.R:1)
  remains `782` lines with `4` top-level functions, a `52`-line maximum
  helper span, and `direct-extraction-ready` classification.
- This stop point now keeps the staged final reporting/dispatch helper file
  pending review; the next bounded checkpoint should review and apply the
  staged final wave so no other live general-plotting extraction candidates
  remain in `R/func_general_plotting.R`.
- April 18, 2026 the reviewed reporting/dispatch wave was applied live from
  [tools/refactor/manifest-general-plotting-reporting-dispatch-wave7.yml](/home/doktersmol/Documents/MultiScholaR/tools/refactor/manifest-general-plotting-reporting-dispatch-wave7.yml:1)
  into:
  - [R/func_general_plotting_reporting_dispatch_helpers.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_plotting_reporting_dispatch_helpers.R:1)
  - [R/func_general_plotting.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_plotting.R:1)
- The live apply materialized the final bounded reporting/dispatch helper
  slice for:
  - `printPValuesDistribution()`
  - `gg_save_logging()`
  - `summarizeQCPlot()`
  - `plotPcaDispatch()`
- The emitted collate artifact now exists at
  [tools/refactor/collate-general-plotting-reporting-dispatch-wave7.txt](/home/doktersmol/Documents/MultiScholaR/tools/refactor/collate-general-plotting-reporting-dispatch-wave7.txt:1),
  and `DESCRIPTION` now collates
  `R/func_general_plotting_reporting_dispatch_helpers.R` immediately before
  `R/func_general_plotting.R`.
- The focused gate stayed green before apply and reran green after the live
  apply through `Rscript tools/test_with_renv.R tests/testthat/test-prot-04-design.R`
  with the same `1572` passing expectations and `1` skipped Git LFS-backed
  snapshot check.
- [R/func_general_plotting.R](/home/doktersmol/Documents/MultiScholaR/R/func_general_plotting.R:1)
  now contains `0` top-level functions and no remaining live
  general-plotting extraction candidates; the file is retained as a
  `637`-line breadcrumb scaffold for the shared plotting surface.
- Bucket 13 manual target is complete; keep this handover as archival and move
  general-plotting follow-up work only if a later regression requires it.
