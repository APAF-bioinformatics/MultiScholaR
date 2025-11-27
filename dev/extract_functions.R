# ============================================================================
# extract_functions.R
# ============================================================================
# Purpose: Safely extract functions from monolithic R files to organized
#          func_*.R files following golem conventions.
#
# Usage:
#   source("dev/extract_functions.R")
#   extract_all(dry_run = TRUE)   # Preview what will be extracted
#   extract_all(dry_run = FALSE)  # Execute extraction
#
# Author: MultiScholaR Development Team
# ============================================================================

# --- Configuration -----------------------------------------------------------

#' Extraction configuration list
#' Each entry specifies: target file, source file, and functions to extract
EXTRACTIONS <- list(
  
  # ==========================================================================
  # PROTEOMICS FUNCTIONS
  # ==========================================================================
  
  # --- func_prot_import.R - Proteomics data import functions ---
  list(
    target = "R/func_prot_import.R"
    , source = "R/mod_prot_import.R"
    , functions = c(
      "detectProteomicsFormat"
      , "importDIANNData"
      , "importSpectronautData"
      , "importFragPipeData"
      , "importMaxQuantData"
      , "getDefaultProteomicsConfig"
    )
  )
  , list(
    target = "R/func_prot_import.R"
    , source = "R/file_management.R"
    , functions = c(
      "importProteomeDiscovererTMTData"
      , "formatDIANN"
      , "formatDIANNParquet"
    )
  )
  
  # --- func_prot_qc_peptide.R - Peptide-level QC filtering ---
  # S4 methods from peptideVsSamplesS4Objects.R
  , list(
    target = "R/func_prot_qc_peptide.R"
    , source = "R/peptideVsSamplesS4Objects.R"
    , functions = c(
      "peptideIntensityFiltering"
      , "removePeptidesWithMissingValuesPercent"
      , "removePeptidesWithOnlyOneReplicate"
      , "filterMinNumPeptidesPerProtein"
      , "filterMinNumPeptidesPerSample"
      , "srlQvalueProteotypicPeptideClean"
    )
  )
  # Helper functions from qc_and_rollup.R (CORRECTED SOURCE)
  , list(
    target = "R/func_prot_qc_peptide.R"
    , source = "R/qc_and_rollup.R"
    , functions = c(
      "peptideIntensityFilteringHelper"
      , "removePeptidesWithMissingValuesPercentHelper"
      , "removePeptidesWithOnlyOneReplicateHelper"
      , "filterMinNumPeptidesPerProteinHelper"
      , "filterMinNumPeptidesPerSampleHelper"
      , "srlQvalueProteotypicPeptideCleanHelper"
    )
  )
  # NA checking functions from helper_functions.R (CORRECTED SOURCE)
  , list(
    target = "R/func_prot_qc_peptide.R"
    , source = "R/helper_functions.R"
    , functions = c(
      "checkPeptideNAPercentages"
    )
  )
  
  # --- func_prot_qc.R - Protein-level QC filtering ---
  # S4 methods from proteinVsSamplesS4Objects.R
  , list(
    target = "R/func_prot_qc.R"
    , source = "R/proteinVsSamplesS4Objects.R"
    , functions = c(
      "proteinIntensityFiltering"
      , "removeProteinsWithOnlyOneReplicate"
      , "removeRowsWithMissingValuesPercent"
    )
  )
  # Helper functions from qc_and_rollup.R (CORRECTED SOURCE)
  # Note: proteinIntensityFilteringHelper does not exist as standalone function
  , list(
    target = "R/func_prot_qc.R"
    , source = "R/qc_and_rollup.R"
    , functions = c(
      "removeProteinsWithOnlyOneReplicateHelper"
    )
  )
  # Functions from de_proteins_functions.R
  , list(
    target = "R/func_prot_qc.R"
    , source = "R/de_proteins_functions.R"
    , functions = c(
      "removeEmptyRows"
      , "removeRowsWithMissingValues"
      , "removeRowsWithMissingValuesPercentHelper"
    )
  )
  # NA checking functions from helper_functions.R (CORRECTED SOURCE)
  , list(
    target = "R/func_prot_qc.R"
    , source = "R/helper_functions.R"
    , functions = c(
      "checkProteinNAPercentages"
      , "getProteinNARecommendations"
      , "validatePostImputationProteinData"
    )
  )
  
  # --- func_prot_rollup.R - Peptide-to-protein rollup ---
  # S4 methods from peptideVsSamplesS4Objects.R
  , list(
    target = "R/func_prot_rollup.R"
    , source = "R/peptideVsSamplesS4Objects.R"
    , functions = c(
      "rollUpPrecursorToPeptide"
    )
  )
  # Helper functions from qc_and_rollup.R (CORRECTED SOURCE)
  , list(
    target = "R/func_prot_rollup.R"
    , source = "R/qc_and_rollup.R"
    , functions = c(
      "rollUpPrecursorToPeptideHelper"
      , "calcPeptidesPerProtein"
      , "calcTotalPeptides"
      , "countPeptidesPerRun"
      , "count_num_peptides"
      , "countProteinsPerRun"
      , "countUniqueProteins"
      , "count_num_proteins"
      , "count_num_samples"
    )
  )
  
  # --- func_prot_norm.R - Protein normalization ---
  # S4 methods from proteinVsSamplesS4Objects.R
  , list(
    target = "R/func_prot_norm.R"
    , source = "R/proteinVsSamplesS4Objects.R"
    , functions = c(
      "normaliseBetweenSamples"
      , "normaliseUntransformedData"
      , "ruvIII_C_Varying"
      , "ruvCancor"
      , "ruvCancorFast"
      , "getNegCtrlProtAnova"
    )
  )
  # Helper functions from de_proteins_functions.R (CORRECTED SOURCE)
  , list(
    target = "R/func_prot_norm.R"
    , source = "R/de_proteins_functions.R"
    , functions = c(
      "getRuvIIIReplicateMatrixHelper"
      , "getNegCtrlProtAnovaHelper"
      , "findBestK"
      , "findBestKForAssayList"
      , "findBestNegCtrlPercentage"
    )
  )
  # Helper functions from de_proteins_functions.R (CORRECTED SOURCE)
  , list(
    target = "R/func_prot_norm.R"
    , source = "R/de_proteins_functions.R"
    , functions = c(
      "extractRuvResults"
    )
  )
  # Helper functions from qc_and_rollup.R (CORRECTED SOURCE)
  , list(
    target = "R/func_prot_norm.R"
    , source = "R/qc_and_rollup.R"
    , functions = c(
      "scaleCenterAndFillMissing"
    )
  )
  # Helper functions from helper_functions.R
  , list(
    target = "R/func_prot_norm.R"
    , source = "R/helper_functions.R"
    , functions = c(
      "updateRuvParameters"
    )
  )
  
  # --- func_pept_norm.R - Peptide normalization ---
  , list(
    target = "R/func_pept_norm.R"
    , source = "R/peptideVsSamplesS4Objects.R"
    , functions = c(
      "findBestNegCtrlPercentagePeptides"
      , "getNegCtrlProtAnovaPeptides"
      , "log2TransformPeptideMatrix"
    )
  )
  # From qc_and_rollup.R (CORRECTED SOURCE)
  , list(
    target = "R/func_pept_norm.R"
    , source = "R/qc_and_rollup.R"
    , functions = c(
      "log2Transformation"
    )
  )
  
  # --- func_peptide_qc_imputation.R - Missing value imputation ---
  # S4 methods from peptideVsSamplesS4Objects.R
  , list(
    target = "R/func_peptide_qc_imputation.R"
    , source = "R/peptideVsSamplesS4Objects.R"
    , functions = c(
      "peptideMissingValueImputation"
      , "peptideMissingValueImputationLimpa"
    )
  )
  # From limpa_functions.R
  , list(
    target = "R/func_peptide_qc_imputation.R"
    , source = "R/limpa_functions.R"
    , functions = c(
      "proteinMissingValueImputationLimpa"
    )
  )
  # From de_proteins_functions.R (CORRECTED SOURCE)
  , list(
    target = "R/func_peptide_qc_imputation.R"
    , source = "R/de_proteins_functions.R"
    , functions = c(
      "imputePerCol"
    )
  )
  # Validation functions from helper_functions.R (CORRECTED SOURCE)
  , list(
    target = "R/func_peptide_qc_imputation.R"
    , source = "R/helper_functions.R"
    , functions = c(
      "validatePostImputationData"
    )
  )
  
  # --- func_prot_de.R - Protein differential expression ---
  , list(
    target = "R/func_prot_de.R"
    , source = "R/protein_de_analysis_wrapper.R"
    , functions = c(
      "differentialExpressionAnalysis"
      , "differentialExpressionAnalysisHelper"
      , "outputDeResultsAllContrasts"
      , "generateVolcanoPlotGlimma"
      , "generateDEHeatmap"
    )
  )
  , list(
    target = "R/func_prot_de.R"
    , source = "R/de_analysis_function_wrapper.R"
    , functions = c(
      "deAnalysisWrapperFunction"
      , "outputDeAnalysisResults"
    )
  )
  , list(
    target = "R/func_prot_de.R"
    , source = "R/metaboliteVsSamplesS4Objects.R"
    , functions = c(
      "getDeResultsLongFormat"
      , "getDeResultsWideFormat"
    )
  )
  # DE functions from de_proteins_functions.R (including CORRECTED createDeResultsLongFormat)
  , list(
    target = "R/func_prot_de.R"
    , source = "R/de_proteins_functions.R"
    , functions = c(
      "prepareDataForVolcanoPlot"
      , "ebFit"
      , "runTest"
      , "runTests"
      , "runTestsContrasts"
      , "saveDeProteinList"
      , "createDeResultsLongFormat"
    )
  )
  
  # ==========================================================================
  # METABOLOMICS FUNCTIONS
  # ==========================================================================
  
  # --- func_metab_import.R - Metabolomics data import ---
  , list(
    target = "R/func_metab_import.R"
    , source = "R/metaboliteVsSamplesS4Objects.R"
    , functions = c(
      "createMetaboliteAssayData"
    )
  )
  , list(
    target = "R/func_metab_import.R"
    , source = "R/QC_visualisation.R"
    , functions = c(
      "getMetaboliteQuantData"
    )
  )
  
  # --- func_metab_qc.R - Metabolomics QC filtering ---
  , list(
    target = "R/func_metab_qc.R"
    , source = "R/metabolite_qc.R"
    , functions = c(
      "metaboliteIntensityFiltering"
      , "metaboliteIntensityFilteringHelper"
      , "findDuplicateFeatureIDs"
      , "resolveDuplicateFeatures"
      , "resolveDuplicateFeaturesByIntensity"
    )
  )
  , list(
    target = "R/func_metab_qc.R"
    , source = "R/QC_visualisation.R"
    , functions = c(
      "updateMetaboliteFiltering"
      , "getFilteringProgressMetabolomics"
      , "updateFilteringProgressMetabolomics"
      , "countUniqueMetabolites"
      , "countMetabolitesPerSample"
      , "calculateMissingness"
      , "calculateSumIntensityPerSample"
      , "calculateMetaboliteCVs"
      , "getInternalStandardMetrics"
      , "calculateTotalUniqueMetabolitesAcrossAssays"
    )
  )
  , list(
    target = "R/func_metab_qc.R"
    , source = "R/metaboliteVsSamplesS4Objects.R"
    , functions = c(
      "calculateMetabolitePairCorrelation"
    )
  )
  
  # --- func_metab_norm.R - Metabolomics normalization ---
  , list(
    target = "R/func_metab_norm.R"
    , source = "R/metabolite_normalization.R"
    , functions = c(
      "logTransformAssays"
    )
  )
  , list(
    target = "R/func_metab_norm.R"
    , source = "R/metaboliteVsSamplesS4Objects.R"
    , functions = c(
      "getNegCtrlMetabAnova"
    )
  )
  
  # --- func_metab_de.R - Metabolomics differential abundance ---
  , list(
    target = "R/func_metab_de.R"
    , source = "R/metaboliteVsSamplesS4Objects.R"
    , functions = c(
      "differentialAbundanceAnalysis"
      , "differentialAbundanceAnalysisHelper"
    )
  )
  , list(
    target = "R/func_metab_de.R"
    , source = "R/metabolite_de_analysis_wrapper.R"
    , functions = c(
      "getCountsTable"
    )
  )
  
  # ==========================================================================
  # MULTIOMICS FUNCTIONS
  # ==========================================================================
  
  # --- func_multiomics_mofa.R - MOFA integration ---
  , list(
    target = "R/func_multiomics_mofa.R"
    , source = "R/multiomics_functions_MOFA.R"
    , functions = c(
      "plotMofaWeights"
    )
  )
  
  # --- func_multiomics_enrich.R - Multiomics enrichment ---
  , list(
    target = "R/func_multiomics_enrich.R"
    , source = "R/multiomics_enrichment_functions.R"
    , functions = c(
      "submitStringDBEnrichment"
      , "downloadStringDBGraph"
      , "downloadStringDBResultsFile"
      , "retrieveStringDBEnrichmentResults"
      , "runOneStringDbRankEnrichment"
      , "runOneStringDbRankEnrichmentMofa"
      , "runStringDbEnrichmentFromDEResults"
      , "runStringDbEnrichmentFromDEResultsMultiple"
      , "runMetabolomicsEnrichmentAnalysis"
      , "runMetabolomicsPathwayEnrichment"
      , "printStringDbFunctionalEnrichmentBarGraph"
      , "getStringDbSpecies"
      , "searchStringDbSpecies"
    )
  )
  , list(
    target = "R/func_multiomics_enrich.R"
    , source = "R/string_enrichment_functions_refactored.R"
    , functions = c(
      "runStringDbEnrichmentAllContrasts"
    )
  )
  
  # ==========================================================================
  # GENERAL FUNCTIONS
  # ==========================================================================
  
  # --- func_general_design.R - Design matrix functions ---
  , list(
    target = "R/func_general_design.R"
    , source = "R/proteinVsSamplesS4Objects.R"
    , functions = c(
      "cleanDesignMatrix"
    )
  )
  , list(
    target = "R/func_general_design.R"
    , source = "R/peptideVsSamplesS4Objects.R"
    , functions = c(
      "cleanDesignMatrixPeptide"
    )
  )
  # Helper functions from qc_and_rollup.R (CORRECTED SOURCE)
  , list(
    target = "R/func_general_design.R"
    , source = "R/qc_and_rollup.R"
    , functions = c(
      "cleanDesignMatrixCleanCategories"
      , "cleanDesignMatrixCleanCategoriesMap"
      , "cleanDesignMatrixCreateEachVersusAllColumns"
    )
  )
  
  # --- func_general_filemgmt.R - File management ---
  , list(
    target = "R/func_general_filemgmt.R"
    , source = "R/file_management.R"
    , functions = c(
      "setupDirectories"
    )
  )
  , list(
    target = "R/func_general_filemgmt.R"
    , source = "R/qc_and_rollup.R"
    , functions = c(
      "saveListOfPdfs"
    )
  )
  
  # --- func_general_helpers.R - General utility functions ---
  , list(
    target = "R/func_general_helpers.R"
    , source = "R/qc_and_rollup.R"
    , functions = c(
      "saveTimeRecord"
      , "changeToCategorical"
      , "peptidesIntensityMatrixPivotLonger"
      , "proteinIntensityMatrixPivotLonger"
    )
  )
  
  # --- func_general_plotting.R - Plotting functions ---
  , list(
    target = "R/func_general_plotting.R"
    , source = "R/qc_and_rollup.R"
    , functions = c(
      "plotPeptidesProteinsCountsPerSampleHelper"
      , "plotHistogramOfPercentMissingPerIndvidual"
      , "getOneRlePlotData"
      , "plotRleQc"
      , "compareUmapComponentsPairs"
      , "umap_factor_plot"
      , "getCategoricalColourPalette"
      , "getOneContinousPalette"
      , "getContinousColourRules"
      , "getCategoricalAndContinuousColourRules"
      , "getSamplesCorrelationHeatMap"
      , "plotDensityOfProteinIntensityPerSample"
      , "plotPercentSamplesVsProteinQuantified"
      , "getProteinsHeatMap"
      , "apafTheme"
      , "get_color_palette"
    )
  )
  
  # --- func_general_s4_objects.R - S4 class definitions ---
  , list(
    target = "R/func_general_s4_objects.R"
    , source = "R/qc_and_rollup.R"
    , functions = c(
      "FilteringProgress"
    )
  )
  
  # --- Additional func_prot_qc.R functions from qc_and_rollup.R ---
  , list(
    target = "R/func_prot_qc.R"
    , source = "R/qc_and_rollup.R"
    , functions = c(
      "getPairsOfSamplesTable"
      , "calulatePearsonCorrelation"
      , "calculatePearsonCorrelationOptimized"
      , "calulatePearsonCorrelationForSamplePairsHelper"
      , "filterSamplesByPeptideCorrelationThreshold"
      , "findSamplesPairBelowPeptideCorrelationThreshold"
      , "filterSamplesByProteinCorrelationThresholdHelper"
      , "removeProteinWithOnlyOneReplicate"
      , "calculatePercentMissingPeptidePerReplicate"
      , "calculatePercentMissingProteinPerReplicate"
      , "getSamplesCorrelationMatrix"
      , "avgReplicateProteinIntensity"
      , "calculatePercentMissingPerProtein"
      , "calculateMissingValuesPerProteinFishersTest"
    )
  )
  
  # --- Additional func_prot_qc_peptide.R functions from qc_and_rollup.R ---
  , list(
    target = "R/func_prot_qc_peptide.R"
    , source = "R/qc_and_rollup.R"
    , functions = c(
      "removePeptidesOnlyInHek293"
    )
  )
  
  # --- Additional func_peptide_qc_imputation.R functions from qc_and_rollup.R ---
  , list(
    target = "R/func_peptide_qc_imputation.R"
    , source = "R/qc_and_rollup.R"
    , functions = c(
      "peptideMissingValueImputationHelper"
      , "proteinMissingValueImputation"
    )
  )
  
  # ==========================================================================
  # PHASE 2: Additional extractions from de_proteins_functions.R
  # ==========================================================================
  
  # --- func_general_plotting.R - Plotting functions from de_proteins_functions.R ---
  , list(
    target = "R/func_general_plotting.R"
    , source = "R/de_proteins_functions.R"
    , functions = c(
      "plotNumMissingValues"
      , "plotNumOfValues"
      , "plotNumOfValuesNoLog"
      , "plotPcaHelper"
      , "plotPcaListHelper"
      , "plotPcaGgpairs"
      , "plotRleHelper"
      , "getMaxMinBoxplot"
      , "rlePcaPlotList"
      , "plotVolcano"
      , "plotOneVolcano"
      , "plotOneVolcanoNoVerticalLines"
      , "printOneVolcanoPlotWithProteinLabel"
      , "getGlimmaVolcanoProteomics"
      , "getGlimmaVolcanoProteomicsWidget"
      , "getGlimmaVolcanoPhosphoproteomics"
      , "printPValuesDistribution"
      , "gg_save_logging"
    )
  )
  
  # --- func_prot_de.R - Additional DE functions from de_proteins_functions.R ---
  , list(
    target = "R/func_prot_de.R"
    , source = "R/de_proteins_functions.R"
    , functions = c(
      "countStatDeGenes"
      , "countStatDeGenesHelper"
      , "printCountDeGenesTable"
      , "getSignificantData"
      , "getTypeOfGrouping"
      , "extractResults"
    )
  )
  
  # --- func_prot_qc.R - Additional QC functions from de_proteins_functions.R ---
  , list(
    target = "R/func_prot_qc.R"
    , source = "R/de_proteins_functions.R"
    , functions = c(
      "getRowsToKeepList"
      , "averageValuesFromReplicates"
      , "proteinTechRepCorrelationHelper"
    )
  )
  
  # --- func_prot_norm.R - Normalization helpers from de_proteins_functions.R ---
  , list(
    target = "R/func_prot_norm.R"
    , source = "R/de_proteins_functions.R"
    , functions = c(
      "calculateSeparationScore"
      , "calculateCompositeScore"
      , "calculateAdaptiveMaxK"
    )
  )
  
  # --- func_prot_annotation.R - Annotation functions from de_proteins_functions.R ---
  , list(
    target = "R/func_prot_annotation.R"
    , source = "R/de_proteins_functions.R"
    , functions = c(
      "cleanIsoformNumber"
      , "subsetQuery"
      , "batchQueryEvidenceHelper"
      , "batchQueryEvidence"
      , "batchQueryEvidenceHelperGeneId"
      , "batchQueryEvidenceGeneId"
      , "goIdToTerm"
      , "uniprotGoIdToTerm"
      , "getUniprotAnnotations"
      , "directUniprotDownload"
      , "standardizeUniprotColumns"
      , "createEmptyUniprotTable"
    )
  )
  
  # ==========================================================================
  # PHASE 3: Additional extractions from helper_functions.R
  # ==========================================================================
  
  # --- func_general_helpers.R - General utility functions from helper_functions.R ---
  , list(
    target = "R/func_general_helpers.R"
    , source = "R/helper_functions.R"
    , functions = c(
      "createIdToAttributeHash"
      , "convertKeyToAttribute"
      , "setArgsDefault"
      , "getFunctionName"
      , "getFunctionNameSecondLevel"
      , "checkParamsObjectFunctionSimplify"
      , "checkParamsObjectFunctionSimplifyAcceptNull"
      , "updateParamInObject"
      , "updateConfigParameter"
      , "extract_experiment"
      , "formatConfigList"
      , "updateMissingValueParameters"
      , "chooseBestProteinAccession_s3"
    )
  )
  
  # --- func_general_filemgmt.R - File management functions from helper_functions.R ---
  , list(
    target = "R/func_general_filemgmt.R"
    , source = "R/helper_functions.R"
    , functions = c(
      "getProjectPaths"
      , "createDirectoryIfNotExists"
      , "createDirIfNotExists"
      , "sourceRmdFileSimple"
      , "sourceRmdFile"
      , "createOutputDir"
      , "testRequiredFiles"
      , "testRequiredFilesWarning"
      , "testRequiredArguments"
      , "parseType"
      , "parseString"
      , "parseList"
      , "isArgumentDefined"
      , "savePlot"
      , "save_plot"
      , "write_results"
      , "readConfigFile"
      , "readConfigFileSection"
      , "loadDependencies"
      , "setupAndShowDirectories"
      , "copyToResultsSummary"
      , "downloadReportTemplate"
      , "RenderReport"
      , "createStudyParametersFile"
      , "createWorkflowArgsFromConfig"
    )
  )
  
  # --- func_general_plotting.R - Plotting functions from helper_functions.R ---
  , list(
    target = "R/func_general_plotting.R"
    , source = "R/helper_functions.R"
    , functions = c(
      "summarizeQCPlot"
    )
  )
  
  # ==========================================================================
  # PHASE 4: Extractions from enrichment_functions.R
  # ==========================================================================
  
  # --- func_general_enrichment.R - Enrichment analysis functions ---
  , list(
    target = "R/func_general_enrichment.R"
    , source = "R/enrichment_functions.R"
    , functions = c(
      "parseNumList"
      , "convertIdToAnnotation"
      , "oneGoEnrichment"
      , "runOneGoEnrichmentInOutFunction"
      , "convertProteinAccToGeneSymbol"
      , "buildAnnotationIdToAnnotationNameDictionary"
      , "buildOneProteinToAnnotationList"
      , "listifyTableByColumn"
      , "runGsea"
      , "runEnricher"
      , "getUniprotAccToGeneSymbolDictionary"
      , "queryRevigo"
      , "clusterPathways"
      , "getEnrichmentHeatmap"
      , "readEnrichmentResultFiles"
      , "filterResultsWithRevigo"
      , "filterResultsWithRevigoScholar"
      , "saveFilteredFunctionalEnrichmentTable"
      , "evaluateBestMinMaxGeneSetSize"
      , "drawListOfFunctionalEnrichmentHeatmaps"
      , "drawListOfFunctionalEnrichmentHeatmapsScholar"
      , "saveListOfFunctionalEnrichmentHeatmaps"
      , "enrichedPathwayBarPlot"
      , "enrichedGoTermBarPlot"
      , "createWordCloudDataFrame"
      , "cleanDuplicatesEnrichment"
      , "plotEnrichmentBarplot"
      , "list2df"
      , "list2graph"
      , "get_param_change_message"
      , "node_add_alpha"
      , "get_enrichplot_color"
      , "set_enrichplot_color"
      , "add_node_label"
      , "get_ggrepel_segsize"
      , "cnetplotEdited"
      , "fc_readable"
      , "update_n"
      , "extract_geneSets"
      , "enrichProteinsPathwaysHelper"
      , "enrichProteinsPathways"
      , "download_uniprot_data"
      , "uniprotGoIdToTermSimple"
    )
  )
  
  # ==========================================================================
  # PHASE 5: Extractions from annotation.R
  # ==========================================================================
  
  # --- func_prot_annotation.R - Annotation functions from annotation.R ---
  , list(
    target = "R/func_prot_annotation.R"
    , source = "R/annotation.R"
    , functions = c(
      "getUniProtAnnotation"
      , "matchAnnotations"
      , "detectEnsemblIds"
      , "taxonIdToGprofilerOrganism"
      , "convertEnsemblToUniprot"
      , "getUniprotAnnotationsFull"
    )
  )
  
  # ==========================================================================
  # PHASE 6: Extractions from get_best_accession_helper.R
  # ==========================================================================
  
  # --- func_prot_annotation.R - FASTA and accession functions ---
  , list(
    target = "R/func_prot_annotation.R"
    , source = "R/get_best_accession_helper.R"
    , functions = c(
      "getFastaFields"
      , "parseFastaObject"
      , "parseFastaFile"
      , "chooseBestPhosphositeAccession"
      , "chooseBestProteinAccessionHelper"
      , "rankProteinAccessionHelper"
      , "processFastaFile"
      , "updateProteinIDs"
      , "cleanMaxQuantProteins"
      , "processAndFilterData"
      , "saveResults"
    )
  )
  
  # ==========================================================================
  # PHASE 7: Extractions from phosphoproteomics_helper.R
  # ==========================================================================
  
  # --- func_phospho_annotation.R - Phosphoproteomics functions ---
  , list(
    target = "R/func_phospho_annotation.R"
    , source = "R/phosphoproteomics_helper.R"
    , functions = c(
      "addColumnsToEvidenceTbl"
      , "getMaxProb"
      , "getMaxProbFutureMap"
      , "getBestPosition"
      , "getBestPositionFutureMap"
      , "getPosString"
      , "getXMerString"
      , "getXMersList"
      , "formatPhosphositePosition"
      , "removePeptidesWithoutAbundances"
      , "filterPeptideAndExtractProbabilities"
      , "addPeptideStartAndEnd"
      , "addPhosphositesPositionsString"
      , "addXMerStrings"
      , "filterByScoreAndGetSimilarPeptides"
      , "allPhosphositesPivotLonger"
      , "groupParalogPeptides"
      , "allPhosphositesPivotWider"
      , "uniquePhosphositesSummariseLongList"
      , "uniquePhosphositesSummariseWideList"
      , "processMultisiteEvidence"
      , "getUniprotAccRankFromSitesId"
    )
  )
  
  # ==========================================================================
  # PHASE 8: Extractions from functional_enrichment.R
  # ==========================================================================
  
  # --- func_general_enrichment.R - Enrichment functions ---
  , list(
    target = "R/func_general_enrichment.R"
    , source = "R/functional_enrichment.R"
    , functions = c(
      "createDEResultsForEnrichment"
      , "createEnrichmentResults"
      , "perform_enrichment"
      , "generate_enrichment_plots"
      , "summarize_enrichment"
      , "processEnrichments"
      , "getEnrichmentResult"
      , "getEnrichmentPlotly"
      , "getEnrichmentSummary"
    )
  )
  
  # ==========================================================================
  # PHASE 9: Extractions from limpa_functions.R
  # ==========================================================================
  
  # --- func_prot_limpa.R - Limpa helper functions ---
  , list(
    target = "R/func_prot_limpa.R"
    , source = "R/limpa_functions.R"
    , functions = c(
      "generateLimpaQCPlots"
      , "convertDpcDEToStandardFormat"
    )
  )
  
  # ==========================================================================
  # PHASE 10: Additional missed functions (Phase 3 verification)
  # ==========================================================================
  
  # NOTE: plotVolcano is already in Phase 2 extraction list (line ~579)
  
  # --- func_general_helpers.R - Helper functions from qc_and_rollup.R ---
  , list(
    target = "R/func_general_helpers.R"
    , source = "R/qc_and_rollup.R"
    , functions = c(
      "calcHtSize"
    )
  )
  
  # --- func_prot_qc.R - Protein filtering from QC_visualisation.R ---
  , list(
    target = "R/func_prot_qc.R"
    , source = "R/QC_visualisation.R"
    , functions = c(
      "updateProteinFiltering"
    )
  )
  
  # --- func_metab_qc.R - Metabolite filtering plots from QC_visualisation.R ---
  , list(
    target = "R/func_metab_qc.R"
    , source = "R/QC_visualisation.R"
    , functions = c(
      "generateMetaboliteFilteringPlots"
    )
  )
  
  # --- func_multiomics_enrich.R - StringDB and enrichment functions ---
  , list(
    target = "R/func_multiomics_enrich.R"
    , source = "R/string_enrichment_functions_refactored.R"
    , functions = c(
      "plotStringDbEnrichmentResults"
    )
  )
  , list(
    target = "R/func_multiomics_enrich.R"
    , source = "R/multiomics_enrichment_functions.R"
    , functions = c(
      "runKeggEnrichment"
      , "runReactomeEnrichment"
    )
  )
  
  # --- func_prot_de.R - Interactive volcano and data matrix functions ---
  , list(
    target = "R/func_prot_de.R"
    , source = "R/de_analysis_function_wrapper.R"
    , functions = c(
      "writeInteractiveVolcanoPlotProteomics"
      , "writeInteractiveVolcanoPlotProteomicsWidget"
      , "writeInteractiveVolcanoPlotProteomicsMain"
      , "getDataMatrix"
    )
  )
  
  # --- func_prot_import.R - Peptide data constructor ---
  , list(
    target = "R/func_prot_import.R"
    , source = "R/peptideVsSamplesS4Objects.R"
    , functions = c(
      "PeptideQuantitativeDataDiann"
    )
  )
  
  # --- func_prot_qc_peptide.R - Peptide QC helper functions ---
  , list(
    target = "R/func_prot_qc_peptide.R"
    , source = "R/peptideVsSamplesS4Objects.R"
    , functions = c(
      "compareTwoPeptideDataObjects"
      , "summarisePeptideObject"
      , "calculatePeptidePearsonCorrelation"
    )
  )
  
  # --- func_general_plotting.R - PCA plotting ---
  , list(
    target = "R/func_general_plotting.R"
    , source = "R/peptideVsSamplesS4Objects.R"
    , functions = c(
      "plotPca"
    )
  )
  
  # --- RunApplet stays in shiny_applets.R - Shiny utility for applets ---
  # Note: RunApplet is a Shiny utility, keeping in shiny_applets.R is appropriate
)

# --- Core Extraction Functions -----------------------------------------------

#' Find the start line of a function definition
#'
#' Searches for function definitions matching common R patterns:
#' - functionName <- function(
#' - setMethod(f = "functionName"
#' - setGeneric(name = "functionName"
#' - setClass("ClassName"
#'
#' @param lines Character vector of file lines
#' @param func_name Name of function to find
#' @return Integer line number where function definition starts, or NA if not found
#' @noRd
find_function_start <- function(lines, func_name) {
  # Pattern 1: Standard function assignment (functionName <- function)
  pattern_standard <- paste0("^\\s*", func_name, "\\s*(<-|=)\\s*function\\s*\\(")
  

  # Pattern 2: setMethod with f = "funcName" or f="funcName"
  pattern_setmethod <- paste0('setMethod\\s*\\(\\s*f\\s*=\\s*["\']', func_name, '["\']')
  
  # Pattern 3: setGeneric with name = "funcName"
  pattern_setgeneric <- paste0('setGeneric\\s*\\(\\s*name\\s*=\\s*["\']', func_name, '["\']')
  
  # Pattern 4: setClass for class definitions

  pattern_setclass <- paste0('setClass\\s*\\(\\s*["\']', func_name, '["\']')
  
  for (i in seq_along(lines)) {
    line <- lines[i]
    # Skip commented lines - they're not actual function definitions
    if (grepl("^\\s*#", line, perl = TRUE)) {
      next
    }
    if (grepl(pattern_standard, line, perl = TRUE) ||
        grepl(pattern_setmethod, line, perl = TRUE) ||
        grepl(pattern_setgeneric, line, perl = TRUE) ||
        grepl(pattern_setclass, line, perl = TRUE)) {
      return(i)
    }
  }
  
  return(NA_integer_)
}

#' Find the start of comment block preceding a function
#'
#' Looks backwards from function definition to find roxygen (#') or
#' regular (#) comment blocks. Stops at blank line or non-comment code.
#'
#' @param lines Character vector of file lines
#' @param func_start_line Line number where function definition starts
#' @return Integer line number where comment block starts
#' @noRd
find_comment_start <- function(lines, func_start_line) {
  if (func_start_line <= 1) {
    return(1L)
  }
  
  comment_start <- func_start_line
  
  # Look backwards from the line before function definition
  for (i in (func_start_line - 1):1) {
    line <- lines[i]
    trimmed <- trimws(line)
    
    # Empty line signals end of comment block
    if (trimmed == "") {
      break
    }
    
    # Check if it's a comment line (roxygen #' or regular #)
    if (grepl("^#", trimmed)) {
      comment_start <- i
    } else {
      # Non-comment, non-empty line - stop here
      break
    }
  }
  
  return(comment_start)
}

#' Find the end of a function definition
#'
#' Tracks brace depth to find where function body ends.
#' Handles nested braces correctly.
#'
#' @param lines Character vector of file lines
#' @param func_start_line Line number where function definition starts
#' @return Integer line number where function ends
#' @noRd
find_function_end <- function(lines, func_start_line) {
  brace_depth <- 0
  in_function <- FALSE
  in_string <- FALSE
  string_char <- NULL
  
  for (i in func_start_line:length(lines)) {
    line <- lines[i]
    
    # Process character by character to handle strings and braces
    j <- 1
    while (j <= nchar(line)) {
      char <- substr(line, j, j)
      prev_char <- if (j > 1) substr(line, j - 1, j - 1) else ""
      
      # Handle string boundaries (skip escaped quotes)
      if (!in_string && (char == '"' || char == "'")) {
        in_string <- TRUE
        string_char <- char
      } else if (in_string && char == string_char && prev_char != "\\") {
        in_string <- FALSE
        string_char <- NULL
      }
      
      # Only count braces outside strings
      if (!in_string) {
        # Check for comments - ignore rest of line
        if (char == "#") {
          break
        }
        
        if (char == "{") {
          brace_depth <- brace_depth + 1
          in_function <- TRUE
        } else if (char == "}") {
          brace_depth <- brace_depth - 1
          
          # When we return to depth 0 after being inside, function is complete
          if (in_function && brace_depth == 0) {
            return(i)
          }
        }
      }
      
      j <- j + 1
    }
  }
  
  # If we get here, we didn't find a closing brace - return last line
  warning("Could not find closing brace for function starting at line ", func_start_line)
  return(length(lines))
}

#' Extract a single function with its comments from source file
#'
#' @param source_file Path to source R file
#' @param func_name Name of function to extract
#' @return List with success status, extracted code, and line info
#' @noRd
extract_function <- function(source_file, func_name) {
  if (!file.exists(source_file)) {
    return(list(
      success = FALSE
      , error = paste("Source file not found:", source_file)
      , code = NULL
      , start_line = NA
      , end_line = NA
    ))
  }
  
  lines <- readLines(source_file, warn = FALSE)
  
  # Find function definition
  func_start <- find_function_start(lines, func_name)
  
  if (is.na(func_start)) {
    return(list(
      success = FALSE
      , error = paste("Function not found:", func_name, "in", source_file)
      , code = NULL
      , start_line = NA
      , end_line = NA
    ))
  }
  
  # Find comment block start (before function)
  comment_start <- find_comment_start(lines, func_start)
  
  # Find function end
  func_end <- find_function_end(lines, func_start)
  
  # Extract the complete block
  extracted_lines <- lines[comment_start:func_end]
  
  return(list(
    success = TRUE
    , error = NULL
    , code = extracted_lines
    , start_line = comment_start
    , end_line = func_end
    , func_def_line = func_start
  ))
}

#' Check if function already exists in target file
#'
#' @param target_file Path to target R file
#' @param func_name Name of function to check
#' @return Logical TRUE if function exists
#' @noRd
function_exists_in_target <- function(target_file, func_name) {
  if (!file.exists(target_file)) {
    return(FALSE)
  }
  
  lines <- readLines(target_file, warn = FALSE)
  func_start <- find_function_start(lines, func_name)
  
  return(!is.na(func_start))
}

#' Append extracted function to target file
#'
#' @param target_file Path to target R file
#' @param code Character vector of code lines to append
#' @param func_name Name of function (for logging)
#' @return Logical TRUE if successful
#' @noRd
append_to_target <- function(target_file, code, func_name) {
  # Read existing content
  existing <- if (file.exists(target_file)) {
    readLines(target_file, warn = FALSE)
  } else {
    character(0)
  }
  
  # Add separator and function
  separator <- c(
    ""
    , paste0("# ", paste(rep("-", 76), collapse = ""))
    , paste0("# ", func_name)
    , paste0("# ", paste(rep("-", 76), collapse = ""))
  )
  
  new_content <- c(existing, separator, code, "")
  
  writeLines(new_content, target_file)
  
  return(TRUE)
}

# --- Main Extraction Functions -----------------------------------------------

#' Extract all functions according to configuration
#'
#' @param dry_run If TRUE, only preview what will be done (default: TRUE)
#' @param extractions List of extraction configurations (default: EXTRACTIONS)
#' @param log_file Path to log file (default: "dev/extraction_log.md")
#' @return Invisible data frame with extraction results
#' @export
extract_all <- function(dry_run = TRUE
                        , extractions = EXTRACTIONS
                        , log_file = "dev/extraction_log.md") {
  
  results <- data.frame(
    target = character()
    , source = character()
    , func_name = character()
    , status = character()
    , message = character()
    , start_line = integer()
    , end_line = integer()
    , stringsAsFactors = FALSE
  )
  
  cat("\n")
  cat("==========================================================\n")
  cat(if (dry_run) "DRY RUN: " else "EXECUTING: ", "Function Extraction\n")
  cat("==========================================================\n")
  cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  
  for (config in extractions) {
    target <- config$target
    source <- config$source
    functions <- config$functions
    
    cat("Target:", target, "\n")
    cat("Source:", source, "\n")
    cat(paste(rep("-", 50), collapse = ""), "\n")
    
    for (func_name in functions) {
      # Check if already exists in target
      if (function_exists_in_target(target, func_name)) {
        status <- "SKIPPED"
        message <- "Already exists in target"
        cat(sprintf("  [%s] %s - %s\n", status, func_name, message))
        
        results <- rbind(results, data.frame(
          target = target
          , source = source
          , func_name = func_name
          , status = status
          , message = message
          , start_line = NA_integer_
          , end_line = NA_integer_
          , stringsAsFactors = FALSE
        ))
        next
      }
      
      # Extract function
      extraction <- extract_function(source, func_name)
      
      if (!extraction$success) {
        status <- "ERROR"
        message <- extraction$error
        cat(sprintf("  [%s] %s - %s\n", status, func_name, message))
        
        results <- rbind(results, data.frame(
          target = target
          , source = source
          , func_name = func_name
          , status = status
          , message = message
          , start_line = NA_integer_
          , end_line = NA_integer_
          , stringsAsFactors = FALSE
        ))
        next
      }
      
      # Report extraction details
      num_lines <- length(extraction$code)
      has_roxygen <- any(grepl("^#'", extraction$code))
      doc_type <- if (has_roxygen) "roxygen" else if (any(grepl("^#", extraction$code))) "comments" else "none"
      
      if (dry_run) {
        status <- "WOULD_EXTRACT"
        message <- sprintf("Lines %d-%d (%d lines, %s docs)"
                           , extraction$start_line
                           , extraction$end_line
                           , num_lines
                           , doc_type)
        cat(sprintf("  [%s] %s - %s\n", status, func_name, message))
      } else {
        # Actually append to target
        tryCatch({
          append_to_target(target, extraction$code, func_name)
          status <- "EXTRACTED"
          message <- sprintf("Lines %d-%d (%d lines, %s docs)"
                             , extraction$start_line
                             , extraction$end_line
                             , num_lines
                             , doc_type)
          cat(sprintf("  [%s] %s - %s\n", status, func_name, message))
        }, error = function(e) {
          status <<- "ERROR"
          message <<- paste("Write failed:", e$message)
          cat(sprintf("  [%s] %s - %s\n", status, func_name, message))
        })
      }
      
      results <- rbind(results, data.frame(
        target = target
        , source = source
        , func_name = func_name
        , status = status
        , message = message
        , start_line = extraction$start_line
        , end_line = extraction$end_line
        , stringsAsFactors = FALSE
      ))
    }
    cat("\n")
  }
  
  # Summary
  cat("==========================================================\n")
  cat("SUMMARY\n")
  cat("==========================================================\n")
  
  status_counts <- table(results$status)
  for (s in names(status_counts)) {
    cat(sprintf("  %s: %d\n", s, status_counts[s]))
  }
  cat("\n")
  
  # Write logs if not dry run
  if (!dry_run) {
    write_extraction_log(results, log_file)
    cat("Log written to:", log_file, "\n")
    
    # Write error log for functions not found
    error_log_file <- sub("\\.md$", "_errors.md", log_file)
    write_error_log(results, error_log_file)
    cat("Error log written to:", error_log_file, "\n")
  }
  
  invisible(results)
}

#' Write extraction log to markdown file
#'
#' @param results Data frame of extraction results
#' @param log_file Path to log file
#' @noRd
write_extraction_log <- function(results, log_file) {
  log_lines <- c(
    "# Function Extraction Log"
    , ""
    , paste("**Generated:**", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
    , ""
    , "## Summary"
    , ""
  )
  
  status_counts <- table(results$status)
  for (s in names(status_counts)) {
    log_lines <- c(log_lines, sprintf("- **%s**: %d", s, status_counts[s]))
  }
  
  log_lines <- c(log_lines, "", "## Details", "")
  
  # Group by target
  targets <- unique(results$target)
  for (target in targets) {
    log_lines <- c(log_lines, paste("###", target), "")
    
    target_results <- results[results$target == target, ]
    
    log_lines <- c(log_lines
                   , "| Function | Source | Status | Details |"
                   , "|----------|--------|--------|---------|")
    
    for (i in seq_len(nrow(target_results))) {
      row <- target_results[i, ]
      log_lines <- c(log_lines
                     , sprintf("| `%s` | %s | %s | %s |"
                               , row$func_name
                               , basename(row$source)
                               , row$status
                               , row$message))
    }
    
    log_lines <- c(log_lines, "")
  }
  
  writeLines(log_lines, log_file)
}

#' Write error log for functions not found
#'
#' Creates a markdown file listing only the functions that couldn't be found,
#' organized for easy investigation and reclassification.
#'
#' @param results Data frame of extraction results
#' @param error_log_file Path to error log file
#' @noRd
write_error_log <- function(results, error_log_file) {
  # Filter to only errors

error_results <- results[results$status == "ERROR", ]
  
  if (nrow(error_results) == 0) {
    cat("No errors to log.\n")
    return(invisible(NULL))
  }
  
  log_lines <- c(
    "# Functions Not Found - Investigation Required"
    , ""
    , paste("**Generated:**", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
    , ""
    , sprintf("**Total functions not found:** %d", nrow(error_results))
    , ""
    , "These functions were listed in FUNC_CLASSIFICATION.md but could not be"
    , "located in their expected source files. Possible reasons:"
    , ""
    , "1. Function name differs from classification (typo, case sensitivity)"
    , "2. Function is in a different file than documented"
    , "3. Function was renamed or removed"
    , "4. Function is defined inline within another function"
    , "5. Function uses a different definition pattern (e.g., = instead of <-)"
    , ""
    , "---"
    , ""
    , "## Functions by Target File"
    , ""
  )
  
  # Group by target
  targets <- unique(error_results$target)
  for (target in targets) {
    target_errors <- error_results[error_results$target == target, ]
    
    log_lines <- c(log_lines
                   , paste("###", basename(target))
                   , ""
                   , "| Function | Expected Source | Action Needed |"
                   , "|----------|-----------------|---------------|")
    
    for (i in seq_len(nrow(target_errors))) {
      row <- target_errors[i, ]
      log_lines <- c(log_lines
                     , sprintf("| `%s` | %s | [ ] Locate/Reclassify |"
                               , row$func_name
                               , basename(row$source)))
    }
    
    log_lines <- c(log_lines, "")
  }
  
  # Add a section grouped by source file for easier investigation
  log_lines <- c(log_lines
                 , "---"
                 , ""
                 , "## Functions by Source File (for grep investigation)"
                 , ""
                 , "Use these commands to search for functions in the codebase:"
                 , ""
                 , "```bash"
                 , "# Search for a function definition")
  
  sources <- unique(error_results$source)
  for (source in sources) {
    source_errors <- error_results[error_results$source == source, ]
    func_list <- paste(source_errors$func_name, collapse = "|")
    log_lines <- c(log_lines
                   , sprintf("grep -n -E '(%s)' R/*.R", func_list))
  }
  
  log_lines <- c(log_lines
                 , "```"
                 , ""
                 , "## Quick Reference: All Missing Functions"
                 , ""
                 , "```")
  
  for (i in seq_len(nrow(error_results))) {
    row <- error_results[i, ]
    log_lines <- c(log_lines
                   , sprintf("%s -> %s (expected in %s)"
                             , row$func_name
                             , basename(row$target)
                             , basename(row$source)))
  }
  
  log_lines <- c(log_lines, "```", "")
  
  writeLines(log_lines, error_log_file)
}

#' Preview a single function extraction
#'
#' Useful for debugging or inspecting what will be extracted
#'
#' @param source_file Path to source file
#' @param func_name Name of function
#' @param show_code If TRUE, print the extracted code (default: TRUE)
#' @return Invisible extraction result list
#' @export
preview_function <- function(source_file, func_name, show_code = TRUE) {
  result <- extract_function(source_file, func_name)
  
  if (!result$success) {
    cat("ERROR:", result$error, "\n")
    return(invisible(result))
  }
  
  cat("Function:", func_name, "\n")
  cat("Source:", source_file, "\n")
  cat("Lines:", result$start_line, "-", result$end_line, "\n")
  cat("Function def at line:", result$func_def_line, "\n")
  cat("Total lines:", length(result$code), "\n")
  
  has_roxygen <- any(grepl("^#'", result$code))
  cat("Has roxygen:", has_roxygen, "\n")
  
  if (show_code) {
    cat("\n--- Extracted Code ---\n")
    cat(result$code, sep = "\n")
    cat("--- End ---\n")
  }
  
  invisible(result)
}

# --- Entry Point Message -----------------------------------------------------
cat("Function extraction script loaded.\n")
cat("Usage:\n")
cat("  extract_all(dry_run = TRUE)   # Preview extractions\n")
cat("  extract_all(dry_run = FALSE)  # Execute extractions\n")
cat("  preview_function('R/file.R', 'funcName')  # Preview single function\n")
cat("\n")

