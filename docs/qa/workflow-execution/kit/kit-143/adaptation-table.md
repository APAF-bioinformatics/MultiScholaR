# KIT-143: the reference workbook adapted into the integration module

Reference: `APAF-bioinformatics/e_shivani_mofa`, `scripts/mofa_ruvIII.Rmd` and the script it sources,
`scripts/rank_enrichment_stringdb.R` (read as code only, 2026-11-21; no data file opened). Chunk numbers are the
reference workbook's chunk labels. The package's workbook is `inst/workbooks/integration/integration_workflow_starter.rmd`
(sections 1 to 10); the app's Integration tab calls the same functions.

**Claim (audit SCI-14, user 2026-11-21).** The module reproduces the reference workbook's MOFA2 configuration exactly
on synthetic data (`tests/testthat/test-integration-workbook-mofa.R`). This is implementation conformance: it says
nothing about the module's defaults, the collaborator's data or the collaborator's results.

**Input matrix (audit SCI-04).** The reference workbook's inputs were RUV-corrected with the groups as replicate sets
and then compared between the same groups. Settled by INT-23: the module takes each omic's normalised matrix before
RUV correction by default; the corrected matrix stays selectable (KIT-153, added to the workbook when it lands), and
a factor-against-group test is withheld when an input was corrected with the groups (the claim limit, KIT-060 D10).

## Steps

| # | Reference workbook step (chunk) | Package equivalent | Adopted? |
|---|---|---|---|
| 1 | `installMultiScholaR()`, `pacman` loads, runtime installs (first chunks) | Dependencies declared by the package and its installer | Not adopted: no installs at run time (conventions) |
| 2 | `setupDirectories(omic_types = "integration")` (3) | Workbook section 2: the same folder layout (`results/integration/mofa/...`, `results_summary/integration/Study_report`); the app copies the workbook into `scripts/integration/` | Adopted (layout) |
| 3 | `vroom` reads of four files named by hand (4) | Each omic's Summary writes `integration/<omic>/` (contract v2, Parquet; INT-10, KIT-147); `integrationScanProject`, `importIntegrationInputs` (section 3) | Replaced |
| 4 | `column_to_rownames("Protein.Ids")` (6) | The manifest names the feature id column; `readIntegrationInput` | Replaced |
| 5 | `limpa::dpcImpute` (Seer), `missForest` (lipidomics, metabolomics) | none | Not adopted, because the imputed values were never given to the model (chunk 7 uses the unimputed tables); MOFA2 handles missing values; imputation is out of scope (plan). Stated as a limitation |
| 6 | Long table of imputed values written to `imputed_data.tsv` | none | Not adopted, because unused |
| 7 | A sample lost by one omic removed from all omics by hand (7) | `missing_samples`: `keep_missing` by default (KIT-060 D4, INT-010), `drop_incomplete` selectable (the reference behaviour) | Adopted as an option |
| 8 | Samples matched by sort order: `x[, sort(colnames(x))]` per omic and `sample_names_lookup` (7) | Matching by the biological sample each omic saved (INT-011), with an editable match table (`integrationSampleMatchTemplate`, `integrationReadSampleMatch`, `integrationAlignmentPreview`; section 4) | Not adopted, because it pairs samples by position: on the corpus it joins different biological samples when run names differ (`tests/testthat/test-integration-workbook.R`) |
| 9 | Technical replicates: none in the reference data | Runs of one biological sample averaged over observed values (INT-2) | Package addition; the upward bias under missing-not-at-random values is stated (SCI-19) |
| 10 | Internal standards (`ITSD`) removed in the metabolomics preparation (prose) | `internal_standards = drop` (INT-3) | Adopted |
| 11 | Protein display names from RefSeq/UniProt lookups (8, not evaluated) | `display_name` from each omic's own feature table (KIT-060 D11) | Replaced |
| 12 | Every feature z-scored, `t(scale(t(x)))` (9) | `scaling = feature_zscore` (default on, INT-24), after the feature selection; features without variance and features seen in too few samples left out | Adopted, with the selection first |
| 13 | No feature selection | `feature_selection = variance_above_trend`, about 1,000 per data set (INT-12); `none` selectable | Package addition (balances data sets of different size) |
| 14 | Design and group vector `Y` by `arrange(Run)` (10) | The group from each omic's saved samples table, matched by biological sample | Replaced |
| 15 | `create_mofa(list of matrices)` (11) | `integrationTrainMofa`: the same call | Adopted |
| 16 | Data options: `center_groups = FALSE`, `scale_views` left FALSE (12) | `center_groups` TRUE by default, FALSE allowed only with the z-score (plan finding 3, KIT-151); `scale_views` TRUE (INT-24) | Adopted as options; defaults changed |
| 17 | Model options: `num_factors = 5` (12) | `num_factors` default 5 | Adopted |
| 18 | Training options left at MOFA2's defaults: seed 42, `fast`, 1,000 rounds, single precision (12) | A declared seed (default 20261119), `slow` to 10,000 rounds (INT-15), five training seeds with each factor's agreement (INT-14), double precision (INT-6, INT-17), the maths-library pin (INT-17, INT-25); every reference value can be set | Adopted as options; defaults changed |
| 19 | `prepare_mofa` (13) | `integrationTrainMofa`: the same call | Adopted |
| 20 | `reticulate::use_virtualenv`, macOS HDF5 paths (14) | The bundled Python (`RETICULATE_PYTHON`), probed (`integrationMofaCapability`); no basilisk fallback (INT-012, INT-036) | Not adopted, because machine-specific |
| 21 | `run_mofa` skipped when a model file exists; model cached with `saveRDS` (14) | A saved model is looked up by its identity (inputs, every setting, MOFA2/mofapy2/numpy versions) and retrained when any changes (KIT-148, KIT-151); the model file is kept in the store and the record bundle (`model.hdf5`) | Not adopted, because a changed input or setting reused the old model |
| 22 | `samples_metadata` with the condition (15) | The aligned sample table with its group | Replaced |
| 23 | "Understanding MOFA factor derivation" demonstration, a correlation of data and reconstruction printed as "R-squared" | none | Not adopted, because the statistic is mislabelled (a correlation, not R²) and the chunk is teaching material |
| 24 | Variance explained tables and plots (16 to 18) | `integrationPlotVariance`, `integrationPlotVarianceTotal`, `integrationResultsVarianceTable` (section 6), stated as variance of the selected features (SCI-12) | Adopted |
| 25 | `plot_factor` by condition with violins (19) | `integrationPlotFactorsByGroup` and the association table (ANOVA or Kruskal-Wallis, adjusted over every factor, INT-20) beside the claim limit | Adopted, with the claim limit |
| 26 | `ggpairs` pair plot recoloured by `help_fix_plot` | `integrationPlotFactorPairs` (ggplot2 only; KIT-152) | Adopted |
| 27 | `plotMofaFeatureWeightsVsRank`: `plot_weights` for Factor 1 per data set (20) | `integrationPlotWeights(results, factor, view)` for every factor and data set | Adopted |
| 28 | `plotMofaWeights`: top 20 by absolute weight, Factor 1 (25) | `integrationPlotWeights`, `integrationResultsWeightsTable` (`top_features`, default 10, ties by order) | Adopted; the old package function's Factor 1 defect not carried (INT-033) |
| 29 | `processOmicForPlot`: boxplots of features matching a search term (27, its use commented out) | none | Not adopted, because study-specific; the top-feature scatter and the saved values cover it |
| 30 | `plot_data_scatter`, top 20 features, Factor 1 (28) | `integrationPlotTopFeatureScatter`, `integrationPlotTopFeatureHeatmap` on the values given to the model (KIT-152) | Adopted |
| 31 | `run_umap(n_neighbors = 4)`, `run_tsne(perplexity = 1)` after `set.seed(42)` (29) | `integrationFactorMap`: off by default, explicit seed, refused at settings invalid for the sample count, a caution beside it (INT-9, KIT-152) | Adopted as an option |
| 32 | Factors, weights and data written as tables (30) | `integrationSaveResults` (section 8) | Adopted (the data handed to the model stay in the store) |
| 33 | `create_pca_plot`: a PCA per omic | none | Not adopted: per-data-set PCA is out of scope (plan); each omic's own QC shows it |
| 34 | Factor 1 against Factor 2 labelled "% variance" from the two factors' variances | none | Not adopted, because the statistic is mislabelled (a share of two factors' variance, not variance explained) |
| 35 | `patchwork` panel and `pca_results.RDS` | none | Not adopted |
| 36 | STRING ranked enrichment on Factor 1's raw weights for Seer and SomaScan (`runOneStringDbRankEnrichmentMofa`, species 9606, FDR 0.05), a cached result reused | KIT-141's enrichment (`runIntegrationEnrichment`): every measured feature ranked by its correlation with the factor (INT-18), KEGG and Reactome for metabolites and lipids (INT-19), key read from `STRING_API_KEY` and never recorded (INT-21) | Replaced (section added when KIT-141 lands) |
| 37 | SomaScan isoforms collapsed to the largest absolute weight | KIT-141's declared identifier mapping per omic | Replaced |
| 38 | `printStringDbFunctionalEnrichmentBarGraph`, top-5 tables, `plot_and_save` | KIT-141's enrichment tables and figures | Replaced (pending KIT-141) |
| 39 | Two inputs of one omic (Seer and SomaScan, both proteomics) | One input per omic | Not adopted: out of scope (plan); stated as a limitation |
| 40 | Inputs RUV-corrected with the groups (file names `*_ruv_normalised_*`) | INT-23 (above) | Not adopted as the default (SCI-04) |
| 41 | Differential-correlation network (`mofa_ruvIII_correlation_network.Rmd`) | none | Not adopted: out of scope (plan) |
| 42 | none | Analysis record, report and methods paragraph (`integrationWriteRecordBundle`, section 9); repeat from the record (`replayIntegrationProject`, section 10) | Package addition |

## Custom functions

| Function (file) | Package equivalent or reason |
|---|---|
| `installMultiScholaR` (mofa_ruvIII.Rmd) | Not adopted: no runtime installs |
| `help_fix_plot` (mofa_ruvIII.Rmd) | `integrationPlotFactorPairs` |
| `plotMofaFeatureWeightsVsRank` (mofa_ruvIII.Rmd) | `integrationPlotWeights` |
| `plotMofaWeights` (mofa_ruvIII.Rmd) | `integrationPlotWeights`, `integrationResultsWeightsTable` |
| `processOmicForPlot` (mofa_ruvIII.Rmd) | Not adopted: study-specific search |
| `create_pca_plot` (mofa_ruvIII.Rmd) | Not adopted: per-data-set PCA out of scope |
| `plot_and_save` (mofa_ruvIII.Rmd) | KIT-141's enrichment figures (pending) |
| `submitStringDBEnrichment`, `retrieveStringDBEnrichmentResults`, `downloadStringDBResultsFile`, `downloadStringDBGraph` (rank_enrichment_stringdb.R) | KIT-141's STRING client on `httr2`: the key sent by POST only, never in a URL, record or log; a network failure refused with no partial result |
| `runOneStringDbRankEnrichment`, `runOneStringDbRankEnrichmentMofa` (rank_enrichment_stringdb.R) | `runIntegrationEnrichment` (KIT-141; pending) |
| `printStringDbFunctionalEnrichmentBarGraph` (rank_enrichment_stringdb.R) | KIT-141's enrichment figure (pending) |
