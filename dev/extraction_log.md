# Function Extraction Log

**Generated:** 2025-11-27 16:10:21

## Summary

- **ERROR**: 13
- **EXTRACTED**: 25
- **SKIPPED**: 315

## Details

### R/func_prot_import.R

| Function | Source | Status | Details |
|----------|--------|--------|---------|
| `detectProteomicsFormat` | mod_prot_import.R | SKIPPED | Already exists in target |
| `importDIANNData` | mod_prot_import.R | SKIPPED | Already exists in target |
| `importSpectronautData` | mod_prot_import.R | SKIPPED | Already exists in target |
| `importFragPipeData` | mod_prot_import.R | SKIPPED | Already exists in target |
| `importMaxQuantData` | mod_prot_import.R | SKIPPED | Already exists in target |
| `getDefaultProteomicsConfig` | mod_prot_import.R | SKIPPED | Already exists in target |
| `importProteomeDiscovererTMTData` | file_management.R | SKIPPED | Already exists in target |
| `formatDIANN` | file_management.R | SKIPPED | Already exists in target |
| `formatDIANNParquet` | file_management.R | SKIPPED | Already exists in target |
| `PeptideQuantitativeDataDiann` | peptideVsSamplesS4Objects.R | SKIPPED | Already exists in target |

### R/func_prot_qc_peptide.R

| Function | Source | Status | Details |
|----------|--------|--------|---------|
| `peptideIntensityFiltering` | peptideVsSamplesS4Objects.R | EXTRACTED | Lines 301-370 (70 lines, roxygen docs) |
| `removePeptidesWithMissingValuesPercent` | peptideVsSamplesS4Objects.R | EXTRACTED | Lines 374-434 (61 lines, roxygen docs) |
| `removePeptidesWithOnlyOneReplicate` | peptideVsSamplesS4Objects.R | EXTRACTED | Lines 528-558 (31 lines, roxygen docs) |
| `filterMinNumPeptidesPerProtein` | peptideVsSamplesS4Objects.R | EXTRACTED | Lines 441-484 (44 lines, roxygen docs) |
| `filterMinNumPeptidesPerSample` | peptideVsSamplesS4Objects.R | EXTRACTED | Lines 489-523 (35 lines, roxygen docs) |
| `srlQvalueProteotypicPeptideClean` | peptideVsSamplesS4Objects.R | EXTRACTED | Lines 194-254 (61 lines, roxygen docs) |
| `peptideIntensityFilteringHelper` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `removePeptidesWithMissingValuesPercentHelper` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `removePeptidesWithOnlyOneReplicateHelper` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `filterMinNumPeptidesPerProteinHelper` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `filterMinNumPeptidesPerSampleHelper` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `srlQvalueProteotypicPeptideCleanHelper` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `checkPeptideNAPercentages` | helper_functions.R | SKIPPED | Already exists in target |
| `removePeptidesOnlyInHek293` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `compareTwoPeptideDataObjects` | peptideVsSamplesS4Objects.R | SKIPPED | Already exists in target |
| `summarisePeptideObject` | peptideVsSamplesS4Objects.R | SKIPPED | Already exists in target |
| `calculatePeptidePearsonCorrelation` | peptideVsSamplesS4Objects.R | SKIPPED | Already exists in target |

### R/func_prot_qc.R

| Function | Source | Status | Details |
|----------|--------|--------|---------|
| `proteinIntensityFiltering` | proteinVsSamplesS4Objects.R | ERROR | Function not found: proteinIntensityFiltering in R/proteinVsSamplesS4Objects.R |
| `removeProteinsWithOnlyOneReplicate` | proteinVsSamplesS4Objects.R | ERROR | Function not found: removeProteinsWithOnlyOneReplicate in R/proteinVsSamplesS4Objects.R |
| `removeRowsWithMissingValuesPercent` | proteinVsSamplesS4Objects.R | ERROR | Function not found: removeRowsWithMissingValuesPercent in R/proteinVsSamplesS4Objects.R |
| `removeProteinsWithOnlyOneReplicateHelper` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `removeEmptyRows` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `removeRowsWithMissingValues` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `removeRowsWithMissingValuesPercentHelper` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `checkProteinNAPercentages` | helper_functions.R | SKIPPED | Already exists in target |
| `getProteinNARecommendations` | helper_functions.R | SKIPPED | Already exists in target |
| `validatePostImputationProteinData` | helper_functions.R | SKIPPED | Already exists in target |
| `getPairsOfSamplesTable` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `calulatePearsonCorrelation` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `calculatePearsonCorrelationOptimized` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `calulatePearsonCorrelationForSamplePairsHelper` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `filterSamplesByPeptideCorrelationThreshold` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `findSamplesPairBelowPeptideCorrelationThreshold` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `filterSamplesByProteinCorrelationThresholdHelper` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `removeProteinWithOnlyOneReplicate` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `calculatePercentMissingPeptidePerReplicate` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `calculatePercentMissingProteinPerReplicate` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `getSamplesCorrelationMatrix` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `avgReplicateProteinIntensity` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `calculatePercentMissingPerProtein` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `calculateMissingValuesPerProteinFishersTest` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `getRowsToKeepList` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `averageValuesFromReplicates` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `proteinTechRepCorrelationHelper` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `updateProteinFiltering` | QC_visualisation.R | SKIPPED | Already exists in target |

### R/func_prot_rollup.R

| Function | Source | Status | Details |
|----------|--------|--------|---------|
| `rollUpPrecursorToPeptide` | peptideVsSamplesS4Objects.R | EXTRACTED | Lines 259-297 (39 lines, roxygen docs) |
| `rollUpPrecursorToPeptideHelper` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `calcPeptidesPerProtein` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `calcTotalPeptides` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `countPeptidesPerRun` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `count_num_peptides` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `countProteinsPerRun` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `countUniqueProteins` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `count_num_proteins` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `count_num_samples` | qc_and_rollup.R | SKIPPED | Already exists in target |

### R/func_prot_norm.R

| Function | Source | Status | Details |
|----------|--------|--------|---------|
| `normaliseBetweenSamples` | proteinVsSamplesS4Objects.R | ERROR | Function not found: normaliseBetweenSamples in R/proteinVsSamplesS4Objects.R |
| `normaliseUntransformedData` | proteinVsSamplesS4Objects.R | ERROR | Function not found: normaliseUntransformedData in R/proteinVsSamplesS4Objects.R |
| `ruvIII_C_Varying` | proteinVsSamplesS4Objects.R | ERROR | Function not found: ruvIII_C_Varying in R/proteinVsSamplesS4Objects.R |
| `ruvCancor` | proteinVsSamplesS4Objects.R | ERROR | Function not found: ruvCancor in R/proteinVsSamplesS4Objects.R |
| `ruvCancorFast` | proteinVsSamplesS4Objects.R | ERROR | Function not found: ruvCancorFast in R/proteinVsSamplesS4Objects.R |
| `getNegCtrlProtAnova` | proteinVsSamplesS4Objects.R | ERROR | Function not found: getNegCtrlProtAnova in R/proteinVsSamplesS4Objects.R |
| `getRuvIIIReplicateMatrixHelper` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `getNegCtrlProtAnovaHelper` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `findBestK` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `findBestKForAssayList` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `findBestNegCtrlPercentage` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `extractRuvResults` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `scaleCenterAndFillMissing` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `updateRuvParameters` | helper_functions.R | SKIPPED | Already exists in target |
| `calculateSeparationScore` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `calculateCompositeScore` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `calculateAdaptiveMaxK` | de_proteins_functions.R | SKIPPED | Already exists in target |

### R/func_pept_norm.R

| Function | Source | Status | Details |
|----------|--------|--------|---------|
| `findBestNegCtrlPercentagePeptides` | peptideVsSamplesS4Objects.R | EXTRACTED | Lines 1112-1347 (236 lines, roxygen docs) |
| `getNegCtrlProtAnovaPeptides` | peptideVsSamplesS4Objects.R | EXTRACTED | Lines 767-813 (47 lines, roxygen docs) |
| `log2TransformPeptideMatrix` | peptideVsSamplesS4Objects.R | EXTRACTED | Lines 2155-2218 (64 lines, roxygen docs) |
| `log2Transformation` | qc_and_rollup.R | SKIPPED | Already exists in target |

### R/func_peptide_qc_imputation.R

| Function | Source | Status | Details |
|----------|--------|--------|---------|
| `peptideMissingValueImputation` | peptideVsSamplesS4Objects.R | EXTRACTED | Lines 584-626 (43 lines, roxygen docs) |
| `peptideMissingValueImputationLimpa` | peptideVsSamplesS4Objects.R | EXTRACTED | Lines 1551-1759 (209 lines, roxygen docs) |
| `proteinMissingValueImputationLimpa` | limpa_functions.R | EXTRACTED | Lines 3-216 (214 lines, roxygen docs) |
| `imputePerCol` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `validatePostImputationData` | helper_functions.R | SKIPPED | Already exists in target |
| `peptideMissingValueImputationHelper` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `proteinMissingValueImputation` | qc_and_rollup.R | SKIPPED | Already exists in target |

### R/func_prot_de.R

| Function | Source | Status | Details |
|----------|--------|--------|---------|
| `differentialExpressionAnalysis` | protein_de_analysis_wrapper.R | EXTRACTED | Lines 16-85 (70 lines, roxygen docs) |
| `differentialExpressionAnalysisHelper` | protein_de_analysis_wrapper.R | EXTRACTED | Lines 89-556 (468 lines, roxygen docs) |
| `outputDeResultsAllContrasts` | protein_de_analysis_wrapper.R | EXTRACTED | Lines 896-1149 (254 lines, roxygen docs) |
| `generateVolcanoPlotGlimma` | protein_de_analysis_wrapper.R | SKIPPED | Already exists in target |
| `generateDEHeatmap` | protein_de_analysis_wrapper.R | SKIPPED | Already exists in target |
| `deAnalysisWrapperFunction` | de_analysis_function_wrapper.R | SKIPPED | Already exists in target |
| `outputDeAnalysisResults` | de_analysis_function_wrapper.R | SKIPPED | Already exists in target |
| `getDeResultsLongFormat` | metaboliteVsSamplesS4Objects.R | EXTRACTED | Lines 2716-2756 (41 lines, roxygen docs) |
| `getDeResultsWideFormat` | metaboliteVsSamplesS4Objects.R | EXTRACTED | Lines 2666-2711 (46 lines, roxygen docs) |
| `prepareDataForVolcanoPlot` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `ebFit` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `runTest` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `runTests` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `runTestsContrasts` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `saveDeProteinList` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `createDeResultsLongFormat` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `countStatDeGenes` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `countStatDeGenesHelper` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `printCountDeGenesTable` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `getSignificantData` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `getTypeOfGrouping` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `extractResults` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `writeInteractiveVolcanoPlotProteomics` | de_analysis_function_wrapper.R | SKIPPED | Already exists in target |
| `writeInteractiveVolcanoPlotProteomicsWidget` | de_analysis_function_wrapper.R | SKIPPED | Already exists in target |
| `writeInteractiveVolcanoPlotProteomicsMain` | de_analysis_function_wrapper.R | SKIPPED | Already exists in target |
| `getDataMatrix` | de_analysis_function_wrapper.R | SKIPPED | Already exists in target |

### R/func_metab_import.R

| Function | Source | Status | Details |
|----------|--------|--------|---------|
| `createMetaboliteAssayData` | metaboliteVsSamplesS4Objects.R | SKIPPED | Already exists in target |
| `getMetaboliteQuantData` | QC_visualisation.R | SKIPPED | Already exists in target |

### R/func_metab_qc.R

| Function | Source | Status | Details |
|----------|--------|--------|---------|
| `metaboliteIntensityFiltering` | metabolite_qc.R | EXTRACTED | Lines 46-172 (127 lines, roxygen docs) |
| `metaboliteIntensityFilteringHelper` | metabolite_qc.R | SKIPPED | Already exists in target |
| `findDuplicateFeatureIDs` | metabolite_qc.R | SKIPPED | Already exists in target |
| `resolveDuplicateFeatures` | metabolite_qc.R | ERROR | Function not found: resolveDuplicateFeatures in R/metabolite_qc.R |
| `resolveDuplicateFeaturesByIntensity` | metabolite_qc.R | SKIPPED | Already exists in target |
| `updateMetaboliteFiltering` | QC_visualisation.R | SKIPPED | Already exists in target |
| `getFilteringProgressMetabolomics` | QC_visualisation.R | SKIPPED | Already exists in target |
| `updateFilteringProgressMetabolomics` | QC_visualisation.R | SKIPPED | Already exists in target |
| `countUniqueMetabolites` | QC_visualisation.R | SKIPPED | Already exists in target |
| `countMetabolitesPerSample` | QC_visualisation.R | SKIPPED | Already exists in target |
| `calculateMissingness` | QC_visualisation.R | SKIPPED | Already exists in target |
| `calculateSumIntensityPerSample` | QC_visualisation.R | SKIPPED | Already exists in target |
| `calculateMetaboliteCVs` | QC_visualisation.R | SKIPPED | Already exists in target |
| `getInternalStandardMetrics` | QC_visualisation.R | SKIPPED | Already exists in target |
| `calculateTotalUniqueMetabolitesAcrossAssays` | QC_visualisation.R | SKIPPED | Already exists in target |
| `calculateMetabolitePairCorrelation` | metaboliteVsSamplesS4Objects.R | SKIPPED | Already exists in target |
| `generateMetaboliteFilteringPlots` | QC_visualisation.R | SKIPPED | Already exists in target |

### R/func_metab_norm.R

| Function | Source | Status | Details |
|----------|--------|--------|---------|
| `logTransformAssays` | metabolite_normalization.R | EXTRACTED | Lines 3-149 (147 lines, roxygen docs) |
| `getNegCtrlMetabAnova` | metaboliteVsSamplesS4Objects.R | EXTRACTED | Lines 1546-1838 (293 lines, roxygen docs) |

### R/func_metab_de.R

| Function | Source | Status | Details |
|----------|--------|--------|---------|
| `differentialAbundanceAnalysis` | metaboliteVsSamplesS4Objects.R | ERROR | Function not found: differentialAbundanceAnalysis in R/metaboliteVsSamplesS4Objects.R |
| `differentialAbundanceAnalysisHelper` | metaboliteVsSamplesS4Objects.R | ERROR | Function not found: differentialAbundanceAnalysisHelper in R/metaboliteVsSamplesS4Objects.R |
| `getCountsTable` | metabolite_de_analysis_wrapper.R | SKIPPED | Already exists in target |

### R/func_multiomics_mofa.R

| Function | Source | Status | Details |
|----------|--------|--------|---------|
| `plotMofaWeights` | multiomics_functions_MOFA.R | SKIPPED | Already exists in target |

### R/func_multiomics_enrich.R

| Function | Source | Status | Details |
|----------|--------|--------|---------|
| `submitStringDBEnrichment` | multiomics_enrichment_functions.R | SKIPPED | Already exists in target |
| `downloadStringDBGraph` | multiomics_enrichment_functions.R | SKIPPED | Already exists in target |
| `downloadStringDBResultsFile` | multiomics_enrichment_functions.R | SKIPPED | Already exists in target |
| `retrieveStringDBEnrichmentResults` | multiomics_enrichment_functions.R | SKIPPED | Already exists in target |
| `runOneStringDbRankEnrichment` | multiomics_enrichment_functions.R | SKIPPED | Already exists in target |
| `runOneStringDbRankEnrichmentMofa` | multiomics_enrichment_functions.R | SKIPPED | Already exists in target |
| `runStringDbEnrichmentFromDEResults` | multiomics_enrichment_functions.R | SKIPPED | Already exists in target |
| `runStringDbEnrichmentFromDEResultsMultiple` | multiomics_enrichment_functions.R | SKIPPED | Already exists in target |
| `runMetabolomicsEnrichmentAnalysis` | multiomics_enrichment_functions.R | SKIPPED | Already exists in target |
| `runMetabolomicsPathwayEnrichment` | multiomics_enrichment_functions.R | SKIPPED | Already exists in target |
| `printStringDbFunctionalEnrichmentBarGraph` | multiomics_enrichment_functions.R | SKIPPED | Already exists in target |
| `getStringDbSpecies` | multiomics_enrichment_functions.R | SKIPPED | Already exists in target |
| `searchStringDbSpecies` | multiomics_enrichment_functions.R | SKIPPED | Already exists in target |
| `runStringDbEnrichmentAllContrasts` | string_enrichment_functions_refactored.R | SKIPPED | Already exists in target |
| `plotStringDbEnrichmentResults` | string_enrichment_functions_refactored.R | SKIPPED | Already exists in target |
| `runKeggEnrichment` | multiomics_enrichment_functions.R | SKIPPED | Already exists in target |
| `runReactomeEnrichment` | multiomics_enrichment_functions.R | SKIPPED | Already exists in target |

### R/func_general_design.R

| Function | Source | Status | Details |
|----------|--------|--------|---------|
| `cleanDesignMatrix` | proteinVsSamplesS4Objects.R | ERROR | Function not found: cleanDesignMatrix in R/proteinVsSamplesS4Objects.R |
| `cleanDesignMatrixPeptide` | peptideVsSamplesS4Objects.R | EXTRACTED | Lines 174-189 (16 lines, roxygen docs) |
| `cleanDesignMatrixCleanCategories` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `cleanDesignMatrixCleanCategoriesMap` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `cleanDesignMatrixCreateEachVersusAllColumns` | qc_and_rollup.R | SKIPPED | Already exists in target |

### R/func_general_filemgmt.R

| Function | Source | Status | Details |
|----------|--------|--------|---------|
| `setupDirectories` | file_management.R | SKIPPED | Already exists in target |
| `saveListOfPdfs` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `getProjectPaths` | helper_functions.R | SKIPPED | Already exists in target |
| `createDirectoryIfNotExists` | helper_functions.R | SKIPPED | Already exists in target |
| `createDirIfNotExists` | helper_functions.R | SKIPPED | Already exists in target |
| `sourceRmdFileSimple` | helper_functions.R | SKIPPED | Already exists in target |
| `sourceRmdFile` | helper_functions.R | SKIPPED | Already exists in target |
| `createOutputDir` | helper_functions.R | SKIPPED | Already exists in target |
| `testRequiredFiles` | helper_functions.R | SKIPPED | Already exists in target |
| `testRequiredFilesWarning` | helper_functions.R | SKIPPED | Already exists in target |
| `testRequiredArguments` | helper_functions.R | SKIPPED | Already exists in target |
| `parseType` | helper_functions.R | SKIPPED | Already exists in target |
| `parseString` | helper_functions.R | SKIPPED | Already exists in target |
| `parseList` | helper_functions.R | SKIPPED | Already exists in target |
| `isArgumentDefined` | helper_functions.R | SKIPPED | Already exists in target |
| `savePlot` | helper_functions.R | SKIPPED | Already exists in target |
| `save_plot` | helper_functions.R | SKIPPED | Already exists in target |
| `write_results` | helper_functions.R | SKIPPED | Already exists in target |
| `readConfigFile` | helper_functions.R | SKIPPED | Already exists in target |
| `readConfigFileSection` | helper_functions.R | SKIPPED | Already exists in target |
| `loadDependencies` | helper_functions.R | SKIPPED | Already exists in target |
| `setupAndShowDirectories` | helper_functions.R | SKIPPED | Already exists in target |
| `copyToResultsSummary` | helper_functions.R | SKIPPED | Already exists in target |
| `downloadReportTemplate` | helper_functions.R | SKIPPED | Already exists in target |
| `RenderReport` | helper_functions.R | SKIPPED | Already exists in target |
| `createStudyParametersFile` | helper_functions.R | SKIPPED | Already exists in target |
| `createWorkflowArgsFromConfig` | helper_functions.R | SKIPPED | Already exists in target |

### R/func_general_helpers.R

| Function | Source | Status | Details |
|----------|--------|--------|---------|
| `saveTimeRecord` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `changeToCategorical` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `peptidesIntensityMatrixPivotLonger` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `proteinIntensityMatrixPivotLonger` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `createIdToAttributeHash` | helper_functions.R | SKIPPED | Already exists in target |
| `convertKeyToAttribute` | helper_functions.R | SKIPPED | Already exists in target |
| `setArgsDefault` | helper_functions.R | SKIPPED | Already exists in target |
| `getFunctionName` | helper_functions.R | SKIPPED | Already exists in target |
| `getFunctionNameSecondLevel` | helper_functions.R | SKIPPED | Already exists in target |
| `checkParamsObjectFunctionSimplify` | helper_functions.R | SKIPPED | Already exists in target |
| `checkParamsObjectFunctionSimplifyAcceptNull` | helper_functions.R | SKIPPED | Already exists in target |
| `updateParamInObject` | helper_functions.R | SKIPPED | Already exists in target |
| `updateConfigParameter` | helper_functions.R | SKIPPED | Already exists in target |
| `extract_experiment` | helper_functions.R | SKIPPED | Already exists in target |
| `formatConfigList` | helper_functions.R | SKIPPED | Already exists in target |
| `updateMissingValueParameters` | helper_functions.R | SKIPPED | Already exists in target |
| `chooseBestProteinAccession_s3` | helper_functions.R | SKIPPED | Already exists in target |
| `calcHtSize` | qc_and_rollup.R | SKIPPED | Already exists in target |

### R/func_general_plotting.R

| Function | Source | Status | Details |
|----------|--------|--------|---------|
| `plotPeptidesProteinsCountsPerSampleHelper` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `plotHistogramOfPercentMissingPerIndvidual` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `getOneRlePlotData` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `plotRleQc` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `compareUmapComponentsPairs` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `umap_factor_plot` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `getCategoricalColourPalette` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `getOneContinousPalette` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `getContinousColourRules` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `getCategoricalAndContinuousColourRules` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `getSamplesCorrelationHeatMap` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `plotDensityOfProteinIntensityPerSample` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `plotPercentSamplesVsProteinQuantified` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `getProteinsHeatMap` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `apafTheme` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `get_color_palette` | qc_and_rollup.R | SKIPPED | Already exists in target |
| `plotNumMissingValues` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `plotNumOfValues` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `plotNumOfValuesNoLog` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `plotPcaHelper` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `plotPcaListHelper` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `plotPcaGgpairs` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `plotRleHelper` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `getMaxMinBoxplot` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `rlePcaPlotList` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `plotVolcano` | de_proteins_functions.R | EXTRACTED | Lines 1147-1187 (41 lines, roxygen docs) |
| `plotOneVolcano` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `plotOneVolcanoNoVerticalLines` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `printOneVolcanoPlotWithProteinLabel` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `getGlimmaVolcanoProteomics` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `getGlimmaVolcanoProteomicsWidget` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `getGlimmaVolcanoPhosphoproteomics` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `printPValuesDistribution` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `gg_save_logging` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `summarizeQCPlot` | helper_functions.R | SKIPPED | Already exists in target |
| `plotPca` | peptideVsSamplesS4Objects.R | EXTRACTED | Lines 1798-1849 (52 lines, none docs) |

### R/func_general_s4_objects.R

| Function | Source | Status | Details |
|----------|--------|--------|---------|
| `FilteringProgress` | qc_and_rollup.R | EXTRACTED | Lines 2256-2322 (67 lines, roxygen docs) |

### R/func_prot_annotation.R

| Function | Source | Status | Details |
|----------|--------|--------|---------|
| `cleanIsoformNumber` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `subsetQuery` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `batchQueryEvidenceHelper` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `batchQueryEvidence` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `batchQueryEvidenceHelperGeneId` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `batchQueryEvidenceGeneId` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `goIdToTerm` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `uniprotGoIdToTerm` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `getUniprotAnnotations` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `directUniprotDownload` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `standardizeUniprotColumns` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `createEmptyUniprotTable` | de_proteins_functions.R | SKIPPED | Already exists in target |
| `getUniProtAnnotation` | annotation.R | SKIPPED | Already exists in target |
| `matchAnnotations` | annotation.R | SKIPPED | Already exists in target |
| `detectEnsemblIds` | annotation.R | SKIPPED | Already exists in target |
| `taxonIdToGprofilerOrganism` | annotation.R | SKIPPED | Already exists in target |
| `convertEnsemblToUniprot` | annotation.R | SKIPPED | Already exists in target |
| `getUniprotAnnotationsFull` | annotation.R | SKIPPED | Already exists in target |
| `getFastaFields` | get_best_accession_helper.R | SKIPPED | Already exists in target |
| `parseFastaObject` | get_best_accession_helper.R | SKIPPED | Already exists in target |
| `parseFastaFile` | get_best_accession_helper.R | SKIPPED | Already exists in target |
| `chooseBestPhosphositeAccession` | get_best_accession_helper.R | SKIPPED | Already exists in target |
| `chooseBestProteinAccessionHelper` | get_best_accession_helper.R | SKIPPED | Already exists in target |
| `rankProteinAccessionHelper` | get_best_accession_helper.R | SKIPPED | Already exists in target |
| `processFastaFile` | get_best_accession_helper.R | SKIPPED | Already exists in target |
| `updateProteinIDs` | get_best_accession_helper.R | SKIPPED | Already exists in target |
| `cleanMaxQuantProteins` | get_best_accession_helper.R | SKIPPED | Already exists in target |
| `processAndFilterData` | get_best_accession_helper.R | SKIPPED | Already exists in target |
| `saveResults` | get_best_accession_helper.R | SKIPPED | Already exists in target |

### R/func_general_enrichment.R

| Function | Source | Status | Details |
|----------|--------|--------|---------|
| `parseNumList` | enrichment_functions.R | SKIPPED | Already exists in target |
| `convertIdToAnnotation` | enrichment_functions.R | SKIPPED | Already exists in target |
| `oneGoEnrichment` | enrichment_functions.R | SKIPPED | Already exists in target |
| `runOneGoEnrichmentInOutFunction` | enrichment_functions.R | SKIPPED | Already exists in target |
| `convertProteinAccToGeneSymbol` | enrichment_functions.R | SKIPPED | Already exists in target |
| `buildAnnotationIdToAnnotationNameDictionary` | enrichment_functions.R | SKIPPED | Already exists in target |
| `buildOneProteinToAnnotationList` | enrichment_functions.R | SKIPPED | Already exists in target |
| `listifyTableByColumn` | enrichment_functions.R | SKIPPED | Already exists in target |
| `runGsea` | enrichment_functions.R | SKIPPED | Already exists in target |
| `runEnricher` | enrichment_functions.R | SKIPPED | Already exists in target |
| `getUniprotAccToGeneSymbolDictionary` | enrichment_functions.R | SKIPPED | Already exists in target |
| `queryRevigo` | enrichment_functions.R | SKIPPED | Already exists in target |
| `clusterPathways` | enrichment_functions.R | SKIPPED | Already exists in target |
| `getEnrichmentHeatmap` | enrichment_functions.R | SKIPPED | Already exists in target |
| `readEnrichmentResultFiles` | enrichment_functions.R | SKIPPED | Already exists in target |
| `filterResultsWithRevigo` | enrichment_functions.R | SKIPPED | Already exists in target |
| `filterResultsWithRevigoScholar` | enrichment_functions.R | SKIPPED | Already exists in target |
| `saveFilteredFunctionalEnrichmentTable` | enrichment_functions.R | SKIPPED | Already exists in target |
| `evaluateBestMinMaxGeneSetSize` | enrichment_functions.R | SKIPPED | Already exists in target |
| `drawListOfFunctionalEnrichmentHeatmaps` | enrichment_functions.R | SKIPPED | Already exists in target |
| `drawListOfFunctionalEnrichmentHeatmapsScholar` | enrichment_functions.R | SKIPPED | Already exists in target |
| `saveListOfFunctionalEnrichmentHeatmaps` | enrichment_functions.R | SKIPPED | Already exists in target |
| `enrichedPathwayBarPlot` | enrichment_functions.R | SKIPPED | Already exists in target |
| `enrichedGoTermBarPlot` | enrichment_functions.R | SKIPPED | Already exists in target |
| `createWordCloudDataFrame` | enrichment_functions.R | SKIPPED | Already exists in target |
| `cleanDuplicatesEnrichment` | enrichment_functions.R | SKIPPED | Already exists in target |
| `plotEnrichmentBarplot` | enrichment_functions.R | SKIPPED | Already exists in target |
| `list2df` | enrichment_functions.R | SKIPPED | Already exists in target |
| `list2graph` | enrichment_functions.R | SKIPPED | Already exists in target |
| `get_param_change_message` | enrichment_functions.R | SKIPPED | Already exists in target |
| `node_add_alpha` | enrichment_functions.R | SKIPPED | Already exists in target |
| `get_enrichplot_color` | enrichment_functions.R | SKIPPED | Already exists in target |
| `set_enrichplot_color` | enrichment_functions.R | SKIPPED | Already exists in target |
| `add_node_label` | enrichment_functions.R | SKIPPED | Already exists in target |
| `get_ggrepel_segsize` | enrichment_functions.R | SKIPPED | Already exists in target |
| `cnetplotEdited` | enrichment_functions.R | SKIPPED | Already exists in target |
| `fc_readable` | enrichment_functions.R | SKIPPED | Already exists in target |
| `update_n` | enrichment_functions.R | SKIPPED | Already exists in target |
| `extract_geneSets` | enrichment_functions.R | SKIPPED | Already exists in target |
| `enrichProteinsPathwaysHelper` | enrichment_functions.R | SKIPPED | Already exists in target |
| `enrichProteinsPathways` | enrichment_functions.R | SKIPPED | Already exists in target |
| `download_uniprot_data` | enrichment_functions.R | SKIPPED | Already exists in target |
| `uniprotGoIdToTermSimple` | enrichment_functions.R | SKIPPED | Already exists in target |
| `createDEResultsForEnrichment` | functional_enrichment.R | SKIPPED | Already exists in target |
| `createEnrichmentResults` | functional_enrichment.R | SKIPPED | Already exists in target |
| `perform_enrichment` | functional_enrichment.R | SKIPPED | Already exists in target |
| `generate_enrichment_plots` | functional_enrichment.R | SKIPPED | Already exists in target |
| `summarize_enrichment` | functional_enrichment.R | SKIPPED | Already exists in target |
| `processEnrichments` | functional_enrichment.R | SKIPPED | Already exists in target |
| `getEnrichmentResult` | functional_enrichment.R | SKIPPED | Already exists in target |
| `getEnrichmentPlotly` | functional_enrichment.R | SKIPPED | Already exists in target |
| `getEnrichmentSummary` | functional_enrichment.R | SKIPPED | Already exists in target |

### R/func_phospho_annotation.R

| Function | Source | Status | Details |
|----------|--------|--------|---------|
| `addColumnsToEvidenceTbl` | phosphoproteomics_helper.R | SKIPPED | Already exists in target |
| `getMaxProb` | phosphoproteomics_helper.R | SKIPPED | Already exists in target |
| `getMaxProbFutureMap` | phosphoproteomics_helper.R | SKIPPED | Already exists in target |
| `getBestPosition` | phosphoproteomics_helper.R | SKIPPED | Already exists in target |
| `getBestPositionFutureMap` | phosphoproteomics_helper.R | SKIPPED | Already exists in target |
| `getPosString` | phosphoproteomics_helper.R | SKIPPED | Already exists in target |
| `getXMerString` | phosphoproteomics_helper.R | SKIPPED | Already exists in target |
| `getXMersList` | phosphoproteomics_helper.R | SKIPPED | Already exists in target |
| `formatPhosphositePosition` | phosphoproteomics_helper.R | SKIPPED | Already exists in target |
| `removePeptidesWithoutAbundances` | phosphoproteomics_helper.R | SKIPPED | Already exists in target |
| `filterPeptideAndExtractProbabilities` | phosphoproteomics_helper.R | SKIPPED | Already exists in target |
| `addPeptideStartAndEnd` | phosphoproteomics_helper.R | SKIPPED | Already exists in target |
| `addPhosphositesPositionsString` | phosphoproteomics_helper.R | SKIPPED | Already exists in target |
| `addXMerStrings` | phosphoproteomics_helper.R | SKIPPED | Already exists in target |
| `filterByScoreAndGetSimilarPeptides` | phosphoproteomics_helper.R | SKIPPED | Already exists in target |
| `allPhosphositesPivotLonger` | phosphoproteomics_helper.R | SKIPPED | Already exists in target |
| `groupParalogPeptides` | phosphoproteomics_helper.R | SKIPPED | Already exists in target |
| `allPhosphositesPivotWider` | phosphoproteomics_helper.R | SKIPPED | Already exists in target |
| `uniquePhosphositesSummariseLongList` | phosphoproteomics_helper.R | SKIPPED | Already exists in target |
| `uniquePhosphositesSummariseWideList` | phosphoproteomics_helper.R | SKIPPED | Already exists in target |
| `processMultisiteEvidence` | phosphoproteomics_helper.R | SKIPPED | Already exists in target |
| `getUniprotAccRankFromSitesId` | phosphoproteomics_helper.R | SKIPPED | Already exists in target |

### R/func_prot_limpa.R

| Function | Source | Status | Details |
|----------|--------|--------|---------|
| `generateLimpaQCPlots` | limpa_functions.R | SKIPPED | Already exists in target |
| `convertDpcDEToStandardFormat` | limpa_functions.R | SKIPPED | Already exists in target |

