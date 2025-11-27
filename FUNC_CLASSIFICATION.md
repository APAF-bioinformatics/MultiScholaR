# Function Classification Matrix

This document maps all functions in MultiScholaR to their target `func_*.R` files for golem-ification.

## Classification Rules

- **Omics Type**: prot (proteomics), pept (peptides), metab (metabolomics), lipid (lipidomics), multiomics, general
- **Analysis Stage**: import, qc, norm (normalization), de (differential expression), enrich (enrichment)
- **Function Type**: S4 method (exportMethods), regular function (export), internal function

## Function Mappings

### Proteomics Functions

#### func_prot_import.R
- `importDIANNData()` - R/mod_prot_import.R
- `importFragPipeData()` - R/mod_prot_import.R
- `importProteomeDiscovererTMTData()` - R/file_management.R
- `importMaxQuantData()` - R/mod_prot_import.R
- `importSpectronautData()` - R/mod_prot_import.R
- `formatDIANN()` - R/file_management.R
- `formatDIANNParquet()` - R/file_management.R
- `detectProteomicsFormat()` - R/mod_prot_import.R

#### func_prot_qc_peptide.R
- `peptideIntensityFiltering()` - S4 method, R/peptideVsSamplesS4Objects.R
- `peptideIntensityFilteringHelper()` - R/peptideVsSamplesS4Objects.R
- `removePeptidesWithMissingValuesPercent()` - S4 method, R/peptideVsSamplesS4Objects.R
- `removePeptidesWithMissingValuesPercentHelper()` - R/peptideVsSamplesS4Objects.R
- `removePeptidesWithOnlyOneReplicate()` - S4 method, R/peptideVsSamplesS4Objects.R
- `removePeptidesWithOnlyOneReplicateHelper()` - R/peptideVsSamplesS4Objects.R
- `filterMinNumPeptidesPerProtein()` - S4 method, R/peptideVsSamplesS4Objects.R
- `filterMinNumPeptidesPerProteinHelper()` - R/peptideVsSamplesS4Objects.R
- `filterMinNumPeptidesPerSample()` - S4 method, R/peptideVsSamplesS4Objects.R
- `filterMinNumPeptidesPerSampleHelper()` - R/peptideVsSamplesS4Objects.R
- `srlQvalueProteotypicPeptideClean()` - S4 method, R/peptideVsSamplesS4Objects.R
- `srlQvalueProteotypicPeptideCleanHelper()` - R/peptideVsSamplesS4Objects.R
- `checkPeptideNAPercentages()` - R/de_proteins_functions.R
- `removePeptidesOnlyInHek293()` - R/de_proteins_functions.R
- `removePeptidesWithoutAbundances()` - R/de_proteins_functions.R
- `filterByScoreAndGetSimilarPeptides()` - R/enrichment_functions.R
- `filterPeptideAndExtractProbabilities()` - R/enrichment_functions.R
- `groupParalogPeptides()` - R/enrichment_functions.R

#### func_prot_qc.R
- `proteinIntensityFiltering()` - S4 method, R/proteinVsSamplesS4Objects.R
- `proteinIntensityFilteringHelper()` - R/proteinVsSamplesS4Objects.R
- `removeProteinsWithOnlyOneReplicate()` - S4 method, R/proteinVsSamplesS4Objects.R
- `removeProteinWithOnlyOneReplicate()` - R/de_proteins_functions.R
- `removeProteinsWithOnlyOneReplicateHelper()` - R/proteinVsSamplesS4Objects.R
- `removeRowsWithMissingValuesPercent()` - S4 method, R/proteinVsSamplesS4Objects.R
- `removeRowsWithMissingValuesPercentHelper()` - R/proteinVsSamplesS4Objects.R
- `removeEmptyRows()` - R/de_proteins_functions.R
- `removeRowsWithMissingValues()` - R/helper_functions.R
- `calculatePercentMissingPerProtein()` - R/de_proteins_functions.R
- `calculatePercentMissingProteinPerReplicate()` - R/de_proteins_functions.R
- `calculatePercentMissingPeptidePerReplicate()` - R/de_proteins_functions.R
- `calculateMissingValuesPerProteinFishersTest()` - R/de_proteins_functions.R
- `checkProteinNAPercentages()` - R/de_proteins_functions.R
- `getProteinNARecommendations()` - R/de_proteins_functions.R

#### func_prot_rollup.R
- `rollUpPrecursorToPeptide()` - S4 method, R/peptideVsSamplesS4Objects.R
- `rollUpPrecursorToPeptideHelper()` - R/peptideVsSamplesS4Objects.R
- `calcPeptidesPerProtein()` - R/peptideVsSamplesS4Objects.R
- `calcTotalPeptides()` - R/peptideVsSamplesS4Objects.R
- `countPeptidesPerRun()` - R/peptideVsSamplesS4Objects.R
- `countProteinsPerRun()` - R/proteinVsSamplesS4Objects.R
- `countUniqueProteins()` - R/proteinVsSamplesS4Objects.R
- `count_num_peptides()` - R/peptideVsSamplesS4Objects.R
- `count_num_proteins()` - R/proteinVsSamplesS4Objects.R
- `count_num_samples()` - R/proteinVsSamplesS4Objects.R

#### func_prot_norm.R
- `normaliseBetweenSamples()` - S4 method, R/proteinVsSamplesS4Objects.R
- `normaliseUntransformedData()` - S4 method, R/proteinVsSamplesS4Objects.R
- `ruvIII_C_Varying()` (protein method) - S4 method, R/proteinVsSamplesS4Objects.R
- `ruvCancor()` - S4 method, R/proteinVsSamplesS4Objects.R
- `ruvCancorFast()` - S4 method, R/proteinVsSamplesS4Objects.R
- `findBestNegCtrlPercentage()` - R/proteinVsSamplesS4Objects.R
- `findBestK()` - R/proteinVsSamplesS4Objects.R
- `findBestKForAssayList()` - R/proteinVsSamplesS4Objects.R
- `getNegCtrlProtAnova()` - S4 method, R/proteinVsSamplesS4Objects.R
- `getNegCtrlProtAnovaHelper()` - R/proteinVsSamplesS4Objects.R
- `getRuvIIIReplicateMatrixHelper()` - R/proteinVsSamplesS4Objects.R
- `extractRuvResults()` - R/helper_functions.R
- `updateRuvParameters()` - R/helper_functions.R
- `scaleCenterAndFillMissing()` - R/helper_functions.R

#### func_pept_norm.R
- `ruvIII_C_Varying()` (peptide method) - S4 method, R/peptideVsSamplesS4Objects.R
- `findBestNegCtrlPercentagePeptides()` - S4 method, R/peptideVsSamplesS4Objects.R
- `getNegCtrlProtAnovaPeptides()` - S4 method, R/peptideVsSamplesS4Objects.R
- `log2TransformPeptideMatrix()` - S4 method, R/peptideVsSamplesS4Objects.R
- `log2Transformation()` - R/helper_functions.R

#### func_peptide_qc_imputation.R
- `peptideMissingValueImputation()` - S4 method, R/peptideVsSamplesS4Objects.R
- `peptideMissingValueImputationLimpa()` - S4 method, R/peptideVsSamplesS4Objects.R
- `peptideMissingValueImputationHelper()` - R/peptideVsSamplesS4Objects.R
- `imputePerCol()` - R/helper_functions.R
- `validatePostImputationData()` - R/peptideVsSamplesS4Objects.R
- `validatePostImputationProteinData()` - R/proteinVsSamplesS4Objects.R
- `proteinMissingValueImputation()` - R/proteinVsSamplesS4Objects.R
- `proteinMissingValueImputationLimpa()` - S4 method, R/limpa_functions.R

#### func_prot_de.R
- `differentialExpressionAnalysis()` (protein method) - S4 method, R/protein_de_analysis_wrapper.R
- `differentialExpressionAnalysisHelper()` - S4 method, R/protein_de_analysis_wrapper.R
- `deAnalysisWrapperFunction()` - R/de_analysis_function_wrapper.R
- `outputDeResultsAllContrasts()` - S4 method, R/protein_de_analysis_wrapper.R
- `outputDeAnalysisResults()` - R/de_analysis_function_wrapper.R
- `generateVolcanoPlotGlimma()` - R/protein_de_analysis_wrapper.R
- `generateDEHeatmap()` - R/protein_de_analysis_wrapper.R
- `createDeResultsLongFormat()` - R/de_analysis_function_wrapper.R
- `getDeResultsLongFormat()` - S4 method, R/metaboliteVsSamplesS4Objects.R
- `getDeResultsWideFormat()` - S4 method, R/metaboliteVsSamplesS4Objects.R
- `prepareDataForVolcanoPlot()` - R/de_proteins_functions.R
- `ebFit()` - R/de_proteins_functions.R
- `runTest()` - R/de_proteins_functions.R
- `runTests()` - R/de_proteins_functions.R
- `runTestsContrasts()` - R/de_proteins_functions.R
- `saveDeProteinList()` - R/de_proteins_functions.R

### Metabolomics Functions

#### func_metab_import.R
- `createMetaboliteAssayData()` - R/metaboliteVsSamplesS4Objects.R
- `getMetaboliteQuantData()` - R/QC_visualisation.R

#### func_metab_qc.R
- `metaboliteIntensityFiltering()` - S4 method, R/metabolite_qc.R
- `metaboliteIntensityFilteringHelper()` - R/metabolite_qc.R
- `updateMetaboliteFiltering()` - R/QC_visualisation.R
- `getFilteringProgressMetabolomics()` - R/QC_visualisation.R
- `updateFilteringProgressMetabolomics()` - R/QC_visualisation.R
- `countUniqueMetabolites()` - R/QC_visualisation.R
- `countMetabolitesPerSample()` - R/QC_visualisation.R
- `calculateMissingness()` - R/QC_visualisation.R
- `calculateSumIntensityPerSample()` - R/QC_visualisation.R
- `calculateMetaboliteCVs()` - R/QC_visualisation.R
- `calculateMetabolitePairCorrelation()` - R/metaboliteVsSamplesS4Objects.R
- `getInternalStandardMetrics()` - R/QC_visualisation.R
- `calculateTotalUniqueMetabolitesAcrossAssays()` - R/QC_visualisation.R
- `findDuplicateFeatureIDs()` - R/metabolite_qc.R
- `resolveDuplicateFeatures()` - S4 method, R/metabolite_qc.R
- `resolveDuplicateFeaturesByIntensity()` - R/metabolite_qc.R

#### func_metab_norm.R
- `logTransformAssays()` (metabolite method) - S4 method, R/metabolite_normalization.R
- `normaliseUntransformedData()` (metabolite method) - S4 method, R/metabolite_normalization.R
- `normaliseBetweenSamples()` (metabolite method) - S4 method, R/proteinVsSamplesS4Objects.R
- `getNegCtrlMetabAnova()` - S4 method, R/metaboliteVsSamplesS4Objects.R

#### func_metab_de.R
- `differentialAbundanceAnalysis()` (metabolite method) - S4 method, R/metaboliteVsSamplesS4Objects.R
- `differentialAbundanceAnalysisHelper()` - S4 method, R/metaboliteVsSamplesS4Objects.R
- `getCountsTable()` - R/metabolite_de_analysis_wrapper.R, R/metabolites_de_analysis_wrapper.R

### Lipidomics Functions (Placeholders)

#### func_lipid_import.R
- Placeholder for future lipidomics import functions

#### func_lipid_qc.R
- Placeholder for future lipidomics QC functions

#### func_lipid_norm.R
- Placeholder for future lipidomics normalization functions

#### func_lipid_de.R
- Placeholder for future lipidomics DE functions

### Multiomics Functions

#### func_multiomics_mofa.R
- `plotMofaWeights()` - R/multiomics_functions_MOFA.R
- Functions from `multiomics_functions_MOFA.R`

#### func_multiomics_enrich.R
- `submitStringDBEnrichment()` - R/multiomics_enrichment_functions.R
- `downloadStringDBGraph()` - R/multiomics_enrichment_functions.R
- `downloadStringDBResultsFile()` - R/multiomics_enrichment_functions.R
- `retrieveStringDBEnrichmentResults()` - R/multiomics_enrichment_functions.R
- `runOneStringDbRankEnrichment()` - R/multiomics_enrichment_functions.R
- `runOneStringDbRankEnrichmentMofa()` - R/multiomics_enrichment_functions.R
- `runStringDbEnrichmentAllContrasts()` - R/string_enrichment_functions_refactored.R
- `runStringDbEnrichmentFromDEResults()` - R/multiomics_enrichment_functions.R
- `runStringDbEnrichmentFromDEResultsMultiple()` - R/multiomics_enrichment_functions.R
- `runMetabolomicsEnrichmentAnalysis()` - R/multiomics_enrichment_functions.R
- `runMetabolomicsPathwayEnrichment()` - R/multiomics_enrichment_functions.R
- `printStringDbFunctionalEnrichmentBarGraph()` - R/multiomics_enrichment_functions.R
- `getStringDbSpecies()` - R/multiomics_enrichment_functions.R
- `searchStringDbSpecies()` - R/multiomics_enrichment_functions.R

### General Functions

#### func_general_plotting.R
See plan for complete list - includes all plotting utilities

#### func_general_filemgmt.R
See plan for complete list - includes all file management functions

#### func_general_helpers.R
See plan for complete list - includes all helper/utility functions

#### func_general_annotation.R
See plan for complete list - includes all annotation functions

#### func_general_enrichment.R
See plan for complete list - includes all enrichment analysis functions

#### func_general_design.R
- `cleanDesignMatrix()` - S4 method, R/proteinVsSamplesS4Objects.R
- `cleanDesignMatrixPeptide()` - S4 method, R/peptideVsSamplesS4Objects.R
- `cleanDesignMatrixCleanCategories()` - R/helper_functions.R
- `cleanDesignMatrixCleanCategoriesMap()` - R/helper_functions.R
- `cleanDesignMatrixCreateEachVersusAllColumns()` - R/helper_functions.R

#### func_general_s4_generics.R
- All setGeneric() definitions from R/allGenerics.R
- Must be loaded early before setMethod() calls

#### func_prot_s4_objects.R
- `ProteinQuantitativeData` class definition - R/proteinVsSamplesS4Objects.R
- `PeptideQuantitativeData` class definition - R/peptideVsSamplesS4Objects.R
- `PeptideQuantitativeDataDiann()` constructor - R/peptideVsSamplesS4Objects.R
- All proteomics-specific S4 methods from R/proteinVsSamplesS4Objects.R and R/peptideVsSamplesS4Objects.R

#### func_metab_s4_objects.R
- `MetaboliteAssayData` class definition - R/metaboliteVsSamplesS4Objects.R
- `MetabolomicsDifferentialAbundanceResults` class definition - R/metaboliteVsSamplesS4Objects.R
- `createMetaboliteAssayData()` constructor - R/metaboliteVsSamplesS4Objects.R
- All metabolomics-specific S4 methods from R/metaboliteVsSamplesS4Objects.R, R/metabolite_qc.R, R/metabolite_normalization.R

#### func_lipid_s4_objects.R
- Placeholder for future lipidomics S4 classes

#### func_general_s4_objects.R
- `FilteringProgress` class definition - R/qc_and_rollup.R
- `FilteringProgressMetabolomics` class definition - R/qc_and_rollup.R
- Shared S4 utilities

