# R/allGenerics.R
# This file contains all setGeneric definitions for the ProteomeScholaR package.
# Loading this file early ensures that generic functions are defined before
# their corresponding setMethod definitions are encountered during package loading.

# --- From proteinVsSamplesS4Objects.R ---

setGeneric( name ="setProteinData"
            , def=function( theObject, protein_quant_table, protein_id_column) {
                standardGeneric("setProteinData")
            })

setGeneric(name ="cleanDesignMatrix"
           , def=function( theObject) {
             standardGeneric("cleanDesignMatrix")
           })

setGeneric(name="proteinIntensityFiltering"
           , def=function( theObject
                           , proteins_intensity_cutoff_percentile = NULL
                           , proteins_proportion_of_samples_below_cutoff = NULL
                           , core_utilisation = NULL) {
             standardGeneric("proteinIntensityFiltering")
           }
           , signature=c("theObject", "proteins_intensity_cutoff_percentile", "proteins_proportion_of_samples_below_cutoff", "core_utilisation"))

setGeneric(name="removeProteinsWithOnlyOneReplicate"
           , def=function( theObject, core_utilisation = NULL, grouping_variable = NULL) {
             standardGeneric("removeProteinsWithOnlyOneReplicate")
           }
           , signature=c("theObject", "core_utilisation", "grouping_variable"))

setGeneric(name="proteinMissingValueImputationLimpa"
           , def=function(theObject, 
                          dpc_results = NULL,
                          dpc_slope = 0.8,
                          quantified_protein_column = NULL,
                          verbose = TRUE,
                          chunk = 1000) {
             standardGeneric("proteinMissingValueImputationLimpa")
           }
           , signature=c("theObject"))

setGeneric(name="plotRle"
           , def=function( theObject, grouping_variable, yaxis_limit = c(), sample_label = NULL) {
             standardGeneric("plotRle")
           }
           , signature=c("theObject", "grouping_variable", "yaxis_limit", "sample_label"))

setGeneric(name="plotRleList"
           , def=function( theObject, list_of_columns, yaxis_limit = c()) {
             standardGeneric("plotRleList")
           }
           , signature=c("theObject", "list_of_columns", "yaxis_limit"))

setGeneric(name="plotPca"
           , def=function( theObject, grouping_variable, shape_variable = NULL, label_column, title, font_size, cv_percentile = 0.90 ) {
             standardGeneric("plotPca")
           }
           , signature=c("theObject", "grouping_variable", "shape_variable", "label_column", "title", "font_size", "cv_percentile"))

setGeneric(name="plotPcaList"
           , def=function( theObject, grouping_variables_list, label_column, title, font_size, cv_percentile = 0.90 ) {
             standardGeneric("plotPcaList")
           }
           , signature=c("theObject", "grouping_variables_list", "label_column", "title", "font_size", "cv_percentile"))

setGeneric(name="getPcaMatrix"
           , def=function( theObject) {
             standardGeneric("getPcaMatrix")
           }
           , signature=c("theObject"))

setGeneric(name="proteinTechRepCorrelation"
           , def=function( theObject,  tech_rep_num_column = NULL, tech_rep_remove_regex = NULL) {
             standardGeneric("proteinTechRepCorrelation")
           }
           , signature=c("theObject", "tech_rep_num_column", "tech_rep_remove_regex"))

setGeneric(name="plotPearson",
           def=function(theObject, tech_rep_remove_regex = NULL, correlation_group = NA, exclude_pool_samples = TRUE) {
             standardGeneric("plotPearson")
           },
           signature=c("theObject", "tech_rep_remove_regex", "correlation_group", "exclude_pool_samples"))

setGeneric("InitialiseGrid", function(dummy = NULL) {
  standardGeneric("InitialiseGrid")
})

setGeneric(name = "createGridQC",
           def = function(theObject, pca_titles, density_titles, rle_titles, pearson_titles, save_path = NULL, file_name = "pca_density_rle_pearson_corr_plots_merged") {
             standardGeneric("createGridQC")
           },
           signature = c("theObject", "pca_titles", "density_titles", "rle_titles", "pearson_titles", "save_path", "file_name"))

setGeneric(name="normaliseBetweenSamples"
           , def=function( theObject, normalisation_method = NULL) {
             standardGeneric("normaliseBetweenSamples")
           }
           , signature=c("theObject", "normalisation_method"))

setGeneric(name="pearsonCorForSamplePairs"
           , def=function( theObject, tech_rep_remove_regex = NULL, correlation_group = NA, exclude_pool_samples = TRUE) {
             standardGeneric("pearsonCorForSamplePairs")
           }
           , signature=c("theObject", "tech_rep_remove_regex", "correlation_group", "exclude_pool_samples"))

setGeneric(name="getNegCtrlProtAnova"
           , def=function( theObject
                           , ruv_grouping_variable  = NULL
                           , percentage_as_neg_ctrl  = NULL
                           , num_neg_ctrl  = NULL
                           , ruv_qval_cutoff = NULL
                           , ruv_fdr_method = NULL ) {
             standardGeneric("getNegCtrlProtAnova")
           }
           , signature=c("theObject", "ruv_grouping_variable", "num_neg_ctrl", "ruv_qval_cutoff", "ruv_fdr_method"))

setGeneric(name="getLowCoefficientOfVariationProteins"
           , def=function( theObject
                           , percentage_as_neg_ctrl = NULL
                           , num_neg_ctrl = NULL ) {
             standardGeneric("getLowCoefficientOfVariationProteins")
           }
           , signature=c("theObject", "percentage_as_neg_ctrl", "num_neg_ctrl"))

setGeneric(name="ruvCancor"
           , def=function( theObject, ctrl= NULL, num_components_to_impute=NULL, ruv_grouping_variable = NULL ) {
             standardGeneric("ruvCancor")
           }
           , signature=c("theObject", "ctrl", "num_components_to_impute", "ruv_grouping_variable"))

setGeneric(name="ruvCancorFast"
           , def=function( theObject, ctrl= NULL, num_components_to_impute=NULL, ruv_grouping_variable = NULL, simple_imputation_method = "none" ) {
             standardGeneric("ruvCancorFast")
           }
           , signature=c("theObject", "ctrl", "num_components_to_impute", "ruv_grouping_variable", "simple_imputation_method"))

setGeneric(name="getRuvIIIReplicateMatrix"
           , def=function( theObject,  ruv_grouping_variable = NULL) {
             standardGeneric("getRuvIIIReplicateMatrix")
           }
           , signature=c("theObject", "ruv_grouping_variable"))

setGeneric(name="ruvIII_C_Varying"
           , def=function( theObject, ruv_grouping_variable = NULL, ruv_number_k = NULL, ctrl = NULL)  {
             standardGeneric("ruvIII_C_Varying")
           }
           , signature=c("theObject", "ruv_grouping_variable", "ruv_number_k", "ctrl"))

setGeneric(name="removeRowsWithMissingValuesPercent"
           , def=function( theObject
                           , ruv_grouping_variable = NULL
                           , groupwise_percentage_cutoff = NULL
                           , max_groups_percentage_cutoff = NULL
                           , proteins_intensity_cutoff_percentile = NULL ) {
             standardGeneric("removeRowsWithMissingValuesPercent")
           }
           , signature=c("theObject"
                         , "ruv_grouping_variable"
                         , "groupwise_percentage_cutoff"
                         , "max_groups_percentage_cutoff"
                         , "proteins_intensity_cutoff_percentile" ))

setGeneric(name="averageTechReps"
           , def=function( theObject, design_matrix_columns, biological_replicate_column = NULL ) {
             standardGeneric("averageTechReps")
           }
           , signature=c("theObject", "design_matrix_columns", "biological_replicate_column" ))

setGeneric(name="preservePeptideNaValues"
           , def=function( peptide_obj, protein_obj)  {
             standardGeneric("preservePeptideNaValues")
           }
           , signature=c("peptide_obj", "protein_obj" ))

setGeneric(name="chooseBestProteinAccession"
           , def=function(theObject, delim=NULL, seqinr_obj=NULL, seqinr_accession_column=NULL, replace_zero_with_na = NULL, aggregation_method = NULL) {
             standardGeneric("chooseBestProteinAccession")
           }
           , signature=c("theObject", "delim", "seqinr_obj", "seqinr_accession_column"))

setGeneric(name="chooseBestProteinAccessionSumDuplicates" # Note: This might be redundant or specific, verify usage
           , def=function(theObject, delim=NULL, seqinr_obj=NULL, seqinr_accession_column=NULL, replace_zero_with_na = NULL, aggregation_method = NULL) {
             standardGeneric("chooseBestProteinAccessionSumDuplicates")
           }
           , signature=c("theObject", "delim", "seqinr_obj", "seqinr_accession_column"))

setGeneric(name="filterSamplesByProteinCorrelationThreshold"
           , def=function( theObject, threshold = NULL, correlation_group = NULL, tech_rep_remove_regex = NULL) {
             standardGeneric("filterSamplesByProteinCorrelationThreshold")
           }
           , signature=c("theObject", "threshold", "correlation_group", "tech_rep_remove_regex"))

setGeneric(name="plotDensity"
           , def=function( theObject, grouping_variable, title = "", font_size = 8) {
             standardGeneric("plotDensity")
           }
           , signature=c("theObject", "grouping_variable", "title", "font_size")) # Base signature might need refinement based on methods

setGeneric(name="plotPcaBox"
           , def=function( theObject, grouping_variable, title = "", font_size = 8, show_legend = FALSE) {
             standardGeneric("plotPcaBox")
           }
           , signature=c("theObject")) # Dispatch only on theObject to avoid S4 coercion issues

setGeneric(name="plotDensityList"
           , def=function( theObject, grouping_variables_list, title = "", font_size = 8) {
             standardGeneric("plotDensityList")
           }
           , signature=c("theObject", "grouping_variables_list", "title", "font_size"))

# --- From protein_de_analysis_wrapper.R ---

setGeneric(name="differentialExpressionAnalysis"
           , def=function( theObject
                           , contrasts_tbl = NULL
                           , formula_string = NULL
                           , group_id = NULL
                           , de_q_val_thresh = NULL
                           , treat_lfc_cutoff = NULL
                           , eBayes_trend = NULL
                           , eBayes_robust = NULL
                           , args_group_pattern = NULL
                           , args_row_id = NULL
                           , qvalue_column = "fdr_qvalue"
                           , raw_pvalue_column = "raw_pvalue" ) {
             standardGeneric("differentialExpressionAnalysis")
           }
           , signature=c("theObject"))

setGeneric(name="differentialExpressionAnalysisHelper"
           , def=function( theObject
                           , contrasts_tbl = NULL
                           , formula_string = NULL
                           , group_id = NULL
                           , de_q_val_thresh = NULL
                           , treat_lfc_cutoff = NULL
                           , eBayes_trend = NULL
                           , eBayes_robust = NULL
                           , args_group_pattern = NULL
                           , args_row_id = NULL
                           , qvalue_column = "fdr_qvalue"
                           , raw_pvalue_column = "raw_pvalue" ) {
             standardGeneric("differentialExpressionAnalysisHelper")
           }
           , signature=c("theObject"))

setGeneric(name="outputDeResultsAllContrasts"
           , def=function(theObject,
                          de_results_list_all_contrasts = NULL,
                          uniprot_tbl = NULL,
                          de_output_dir = NULL,
                          publication_graphs_dir = NULL,
                          file_prefix = "de_proteins",
                          args_row_id = NULL,
                          gene_names_column = "gene_names",
                          uniprot_id_column = "Entry") {
             standardGeneric("outputDeResultsAllContrasts")
           }
           , signature=c("theObject"))


# --- From peptideVsSamplesS4Objects.R ---

setGeneric(name ="cleanDesignMatrixPeptide"
           , def=function( theObject) {
             standardGeneric("cleanDesignMatrixPeptide")
           })

setGeneric(name="calcPeptideMatrix",
           def=function(theObject) {
             standardGeneric("calcPeptideMatrix")
           },
           signature=c("theObject"))


setGeneric(name="srlQvalueProteotypicPeptideClean"
           , def=function( theObject, qvalue_threshold = NULL, global_qvalue_threshold = NULL, choose_only_proteotypic_peptide = NULL, input_matrix_column_ids = NULL) {
             standardGeneric("srlQvalueProteotypicPeptideClean")
           }
           , signature = c("theObject", "qvalue_threshold", "global_qvalue_threshold", "choose_only_proteotypic_peptide", "input_matrix_column_ids") )

setGeneric(name="rollUpPrecursorToPeptide"
           , def=function( theObject, core_utilisation = NULL) {
             standardGeneric("rollUpPrecursorToPeptide")
           }
           , signature=c("theObject", "core_utilisation"))

setGeneric(name="peptideIntensityFiltering"
           , def=function( theObject, peptides_intensity_cutoff_percentile = NULL, peptides_proportion_of_samples_below_cutoff = NULL, core_utilisation = NULL) {
             standardGeneric("peptideIntensityFiltering")
           }
           , signature=c("theObject", "peptides_intensity_cutoff_percentile", "peptides_proportion_of_samples_below_cutoff", "core_utilisation"))

setGeneric(name="removePeptidesWithMissingValuesPercent"
           , def=function( theObject
                           , grouping_variable = NULL
                           , groupwise_percentage_cutoff = NULL
                           , max_groups_percentage_cutoff = NULL
                           , peptides_intensity_cutoff_percentile = NULL) {
             standardGeneric("removePeptidesWithMissingValuesPercent")
           }
           , signature=c("theObject"
                         , "grouping_variable"
                         , "groupwise_percentage_cutoff"
                         , "max_groups_percentage_cutoff"
                         , "peptides_intensity_cutoff_percentile" ))

setGeneric(name="filterMinNumPeptidesPerProtein"
           , def=function( theObject, ...) {
             standardGeneric("filterMinNumPeptidesPerProtein")
           }
           , signature=c("theObject"))

setGeneric( name="filterMinNumPeptidesPerSample"
            , def=function( theObject, peptides_per_sample_cutoff = NULL, core_utilisation = NULL, inclusion_list = NULL) {
              standardGeneric("filterMinNumPeptidesPerSample")
           }
           , signature=c("theObject", "peptides_per_sample_cutoff", "core_utilisation", "inclusion_list" ))

setGeneric( name="removePeptidesWithOnlyOneReplicate"
            , def=function( theObject, replicate_group_column = NULL, core_utilisation = NULL) {
              standardGeneric("removePeptidesWithOnlyOneReplicate")
            }
            , signature=c("theObject", "replicate_group_column", "core_utilisation" ))

setGeneric( name="plotPeptidesProteinsCountsPerSample"
            , def=function( theObject ) {
              standardGeneric("plotPeptidesProteinsCountsPerSample")
            }
            , signature=c("theObject" ))

setGeneric( name="peptideMissingValueImputation"
            , def=function( theObject,  imputed_value_column = NULL, proportion_missing_values = NULL, core_utilisation = NULL) {
              standardGeneric("peptideMissingValueImputation")
            }
            , signature=c("theObject", "imputed_value_column", "proportion_missing_values", "core_utilisation" ))

setGeneric(name="getNegCtrlProtAnovaPeptides"
           , def=function( theObject
                           , ruv_grouping_variable = NULL
                           , percentage_as_neg_ctrl = NULL
                           , num_neg_ctrl = NULL
                           , ruv_qval_cutoff = NULL
                           , ruv_fdr_method = NULL ) {
             standardGeneric("getNegCtrlProtAnovaPeptides")
           }
           , signature=c("theObject", "ruv_grouping_variable", "percentage_as_neg_ctrl", "num_neg_ctrl", "ruv_qval_cutoff", "ruv_fdr_method"))

setGeneric(name="peptideMissingValueImputationLimpa"
           , def=function(theObject, 
                          imputed_value_column = NULL, 
                          use_log2_transform = TRUE,
                          verbose = TRUE,
                          ensure_matrix = TRUE) {
             standardGeneric("peptideMissingValueImputationLimpa")
           }
           , signature=c("theObject"))

setGeneric(name="findBestNegCtrlPercentagePeptides"
           , def=function(theObject,
                          percentage_range = NULL,
                          num_components_to_impute = 5,
                          ruv_grouping_variable = "group",
                          ruv_qval_cutoff = 0.05,
                          ruv_fdr_method = "qvalue",
                          separation_metric = "max_difference",
                          k_penalty_weight = 0.5,
                          max_acceptable_k = 3,
                          adaptive_k_penalty = TRUE,
                          verbose = TRUE,
                          ensure_matrix = TRUE) {
             standardGeneric("findBestNegCtrlPercentagePeptides")
           }
           , signature=c("theObject"))

setGeneric(name="log2TransformPeptideMatrix"
           , def=function(theObject) {
             standardGeneric("log2TransformPeptideMatrix")
           }
           , signature=c("theObject"))


# --- From metabolite_qc.R ---

setGeneric(name="metaboliteIntensityFiltering"
           , def=function( theObject, metabolites_intensity_cutoff_percentile = NULL, metabolites_proportion_of_samples_below_cutoff = NULL) {
             standardGeneric("metaboliteIntensityFiltering")
           }
           , signature=c("theObject")) # Dispatch only on the object

setGeneric(name = "resolveDuplicateFeatures",
           def = function(theObject, itsd_pattern_columns = NULL) {
             standardGeneric("resolveDuplicateFeatures")
           }
)

# --- From metabolite_normalization.R ---

setGeneric(name = "logTransformAssays",
           def = function(theObject, offset = 1, ...) {
             standardGeneric("logTransformAssays")
           }
)

setGeneric(name = "normaliseUntransformedData",
           def = function(theObject, method = "ITSD", ...) {
             standardGeneric("normaliseUntransformedData")
           }
)

# --- From metaboliteVsSamplesS4Objects.R ---
# (Only generics defined uniquely in this file)

setGeneric(name = "getNegCtrlMetabAnova",
           def = function(theObject,
                          ruv_grouping_variable = NULL,
                          percentage_as_neg_ctrl = NULL, # Allow list/vector/single value
                          num_neg_ctrl = NULL,
                          ruv_qval_cutoff = NULL,
                          ruv_fdr_method = NULL) {
             standardGeneric("getNegCtrlMetabAnova")
           },
           signature = c("theObject")) # Primary dispatch on object 

###### Metabolomics


setGeneric(name="plotVolcano",
           def=function(objectsList,
                        de_q_val_thresh = 0.05,
                        qvalue_column = "q_value",
                        log2fc_column = "logFC") {
             standardGeneric("plotVolcano")
           },
           signature=c("objectsList"))


setGeneric(name="plotVolcanoS4",
           def=function(objectsList,
                        de_q_val_thresh = 0.05,
                        qvalue_column = "fdr_qvalue",
                        log2fc_column = "logFC") {
             standardGeneric("plotVolcanoS4")
           },
           signature=c("objectsList"))


setGeneric(name="getDeResultsWideFormat"
           , def=function(objectsList
                          , qvalue_column = "fdr_qvalue"
                          , raw_pvalue_column = "raw_pvalue"
                          , log2fc_column = "logFC") {
             standardGeneric("getDeResultsWideFormat")
           },
           signature=c("objectsList"))



setGeneric(name="getDeResultsLongFormat"
           , def=function(objectsList) {
             standardGeneric("getDeResultsLongFormat")
           },
           signature=c("objectsList"))



setGeneric(name="plotInteractiveVolcano"
           , def=function(objectsList, anno_list = NULL) {
             standardGeneric("plotInteractiveVolcano")
           },
           signature=c("objectsList"))



setGeneric(name = "createGridQCMetabolomics",
           def = function(theObject, pca_titles, density_titles, rle_titles, pearson_titles, save_path = NULL, file_name = "pca_density_rle_pearson_corr_plots_merged") {
             standardGeneric("createGridQCMetabolomics")
           },
           signature = c("theObject"))


setGeneric(name="plotNumSigDiffExpBarPlot",
           def=function(objectsList ) {
             standardGeneric("plotNumSigDiffExpBarPlot")
           },
           signature=c("objectsList"  ))


setGeneric( name ="differentialAbundanceAnalysisHelper"
            , def=function(theObject
                           , contrasts_tbl = NULL
                           , formula_string = NULL
                           , group_id = NULL
                           , de_q_val_thresh = NULL
                           , treat_lfc_cutoff = NULL
                           , eBayes_trend = NULL
                           , eBayes_robust = NULL
                           , args_group_pattern = NULL) {
              standardGeneric("differentialAbundanceAnalysisHelper")
            }
            , signature=c("theObject"))


setGeneric( name ="differentialAbundanceAnalysis"
            , def=function(objectsList
                           , contrasts_tbl = NULL
                           , formula_string = NULL
                           , group_id = NULL
                           , de_q_val_thresh = NULL
                           , treat_lfc_cutoff = NULL
                           , eBayes_trend = NULL
                           , eBayes_robust = NULL
                           , args_group_pattern = NULL) {
              standardGeneric("differentialAbundanceAnalysis")
            }
            , signature=c("objectsList"))
