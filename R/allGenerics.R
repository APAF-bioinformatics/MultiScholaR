# R/allGenerics.R
# This file contains all setGeneric definitions for the ProteomeScholaR package.
# Loading this file early ensures that generic functions are defined before
# their corresponding setMethod definitions are encountered during package loading.

# --- From proteinVsSamplesS4Objects.R ---

#' @name setProteinData
#' @rdname setProteinData
#' @title Set Protein Quantification Data
#' @description A generic function to set or update the protein quantification data
#' within a proteomics data object.
#' @param theObject The S4 object into which the data will be loaded.
#' @param protein_quant_table A data frame or matrix of protein quantification data.
#' @param protein_id_column The name of the column containing protein identifiers.
#' @return The modified S4 object.
#' @export
setGeneric( name ="setProteinData"
            , def=function( theObject, protein_quant_table, protein_id_column) {
                standardGeneric("setProteinData")
            })

#' @name cleanDesignMatrix
#' @rdname cleanDesignMatrix
#' @title Clean the Design Matrix
#' @description A generic function to perform cleaning operations on the design matrix
#' of a data object, such as removing unused factor levels or standardizing columns.
#' @param theObject The S4 object containing the design matrix.
#' @return The modified S4 object with a cleaned design matrix.
#' @export
setGeneric(name ="cleanDesignMatrix"
           , def=function( theObject) {
             standardGeneric("cleanDesignMatrix")
           })

#' @name proteinIntensityFiltering
#' @rdname proteinIntensityFiltering
#' @title Filter Proteins Based on Intensity
#' @description A generic function to filter low-abundance proteins based on intensity
#' across samples.
#' @param theObject The S4 object with protein data.
#' @param proteins_intensity_cutoff_percentile Numeric (0-100). The percentile of
#' non-zero expression values to use as an intensity cutoff for each protein.
#' @param proteins_proportion_of_samples_below_cutoff Numeric (0-1). The maximum
#' proportion of samples that can be below the intensity cutoff for a protein to be kept.
#' @param core_utilisation Number of CPU cores for parallel processing.
#' @return The filtered S4 object.
#' @export
setGeneric(name="proteinIntensityFiltering"
           , def=function( theObject
                           , proteins_intensity_cutoff_percentile = NULL
                           , proteins_proportion_of_samples_below_cutoff = NULL
                           , core_utilisation = NULL) {
             standardGeneric("proteinIntensityFiltering")
           }
           , signature=c("theObject", "proteins_intensity_cutoff_percentile", "proteins_proportion_of_samples_below_cutoff", "core_utilisation"))

#' @name removeProteinsWithOnlyOneReplicate
#' @rdname removeProteinsWithOnlyOneReplicate
#' @title Remove Proteins with Only a Single Replicate
#' @description A generic function to filter out proteins that are identified in only
#' a single replicate within their experimental group.
#' @param theObject The S4 object.
#' @param core_utilisation Number of CPU cores for parallel processing.
#' @param grouping_variable The column name in the design matrix that defines experimental groups.
#' @return The filtered S4 object.
#' @export
setGeneric(name="removeProteinsWithOnlyOneReplicate"
           , def=function( theObject, core_utilisation = NULL, grouping_variable = NULL) {
             standardGeneric("removeProteinsWithOnlyOneReplicate")
           }
           , signature=c("theObject", "core_utilisation", "grouping_variable"))

#' @name plotRle
#' @rdname plotRle
#' @title Plot Relative Log Expression (RLE)
#' @description A generic function to generate a Relative Log Expression (RLE) plot
#' for visualizing sample-level normalization effects.
#' @param theObject The S4 object.
#' @param grouping_variable The column in the design matrix to use for grouping/coloring.
#' @param yaxis_limit A numeric vector of length 2 defining the y-axis limits.
#' @param sample_label The column in the design matrix to use for labeling samples.
#' @return A ggplot object representing the RLE plot.
#' @export
setGeneric(name="plotRle"
           , def=function( theObject, grouping_variable, yaxis_limit = c(), sample_label = NULL) {
             standardGeneric("plotRle")
           }
           , signature=c("theObject", "grouping_variable", "yaxis_limit", "sample_label"))

#' @name plotRleList
#' @rdname plotRleList
#' @title Plot RLE for a List of Grouping Variables
#' @description A generic function to generate a list of RLE plots, one for each
#' specified grouping variable.
#' @param theObject The S4 object.
#' @param list_of_columns A character vector of column names in the design matrix
#' to use for generating separate RLE plots.
#' @param yaxis_limit A numeric vector of length 2 defining the y-axis limits for all plots.
#' @return A list of ggplot objects.
#' @export
setGeneric(name="plotRleList"
           , def=function( theObject, list_of_columns, yaxis_limit = c()) {
             standardGeneric("plotRleList")
           }
           , signature=c("theObject", "list_of_columns", "yaxis_limit"))

#' @name plotPca
#' @rdname plotPca
#' @title Plot Principal Component Analysis (PCA)
#' @description A generic function to perform PCA and generate a scatter plot of the
#' principal components.
#' @param theObject The S4 object.
#' @param grouping_variable The column in the design matrix for coloring points.
#' @param shape_variable Optional. The column in the design matrix for shaping points.
#' @param label_column The column in the design matrix for labeling points.
#' @param title The plot title.
#' @param font_size The font size for plot text.
#' @return A ggplot object representing the PCA plot.
#' @export
setGeneric(name="plotPca"
           , def=function( theObject, grouping_variable, shape_variable = NULL, label_column, title, font_size ) {
             standardGeneric("plotPca")
           }
           , signature=c("theObject", "grouping_variable", "shape_variable", "label_column", "title", "font_size"))

#' @name plotPcaList
#' @rdname plotPcaList
#' @title Plot PCA for a List of Grouping Variables
#' @description A generic function to generate a list of PCA plots, each colored by a
#' different grouping variable.
#' @param theObject The S4 object.
#' @param grouping_variables_list A character vector of column names in the design
#' matrix to use for coloring separate PCA plots.
#' @param label_column The column in the design matrix for labeling points.
#' @param title The base plot title.
#' @param font_size The font size for plot text.
#' @return A list of ggplot objects.
#' @export
setGeneric(name="plotPcaList"
           , def=function( theObject, grouping_variables_list, label_column, title, font_size ) {
             standardGeneric("plotPcaList")
           }
           , signature=c("theObject", "grouping_variables_list", "label_column", "title", "font_size"))

#' @name getPcaMatrix
#' @rdname getPcaMatrix
#' @title Get PCA Matrix
#' @description A generic function to compute and retrieve the PCA rotation matrix.
#' @param theObject The S4 object.
#' @return A matrix containing the principal components.
#' @export
setGeneric(name="getPcaMatrix"
           , def=function( theObject) {
             standardGeneric("getPcaMatrix")
           }
           , signature=c("theObject"))

#' @name proteinTechRepCorrelation
#' @rdname proteinTechRepCorrelation
#' @title Calculate Technical Replicate Correlation for Proteins
#' @description A generic function to calculate the correlation between technical replicates
#' at the protein level.
#' @param theObject The S4 object.
#' @param tech_rep_num_column The column in the design matrix identifying technical replicates.
#' @param tech_rep_remove_regex A regex to identify and process technical replicate names.
#' @return A data frame containing correlation results.
#' @export
setGeneric(name="proteinTechRepCorrelation"
           , def=function( theObject,  tech_rep_num_column = NULL, tech_rep_remove_regex = NULL) {
             standardGeneric("proteinTechRepCorrelation")
           }
           , signature=c("theObject", "tech_rep_num_column", "tech_rep_remove_regex"))

#' @name plotPearson
#' @rdname plotPearson
#' @title Plot Pearson Correlation Heatmap
#' @description A generic function to generate a heatmap of Pearson correlations
#' between samples.
#' @param theObject The S4 object.
#' @param tech_rep_remove_regex A regex to process sample names for display.
#' @param correlation_group Optional. A column in the design matrix to group samples for correlation.
#' @return A heatmap plot object.
#' @export
setGeneric(name="plotPearson",
           def=function(theObject, tech_rep_remove_regex, correlation_group = NA  ) {
             standardGeneric("plotPearson")
           },
           signature=c("theObject", "tech_rep_remove_regex", "correlation_group" ))

#' @name InitialiseGrid
#' @rdname InitialiseGrid
#' @title Initialize a QC Plot Grid
#' @description A generic function to initialize an empty grid for arranging QC plots.
#' @param dummy A placeholder argument.
#' @return An initialized grid object.
#' @export
setGeneric("InitialiseGrid", function(dummy = NULL) {
  standardGeneric("InitialiseGrid")
})

#' @name createGridQC
#' @rdname createGridQC
#' @title Create a Grid of QC Plots
#' @description A generic function to combine multiple QC plots (PCA, density, RLE, Pearson)
#' into a single grid layout.
#' @param theObject The S4 object.
#' @param pca_titles Titles for the PCA plots.
#' @param density_titles Titles for the density plots.
#' @param rle_titles Titles for the RLE plots.
#' @param pearson_titles Titles for the Pearson correlation plots.
#' @param save_path The directory path to save the output file.
#' @param file_name The name of the output file.
#' @return A grid plot object and saves it to a file.
#' @export
setGeneric(name = "createGridQC",
           def = function(theObject, pca_titles, density_titles, rle_titles, pearson_titles, save_path = NULL, file_name = "pca_density_rle_pearson_corr_plots_merged") {
             standardGeneric("createGridQC")
           },
           signature = c("theObject", "pca_titles", "density_titles", "rle_titles", "pearson_titles", "save_path", "file_name"))

#' @name normaliseBetweenSamples
#' @rdname normaliseBetweenSamples
#' @title Normalize Data Between Samples
#' @description A generic function to apply between-sample normalization methods
#' like quantile normalization.
#' @param theObject The S4 object.
#' @param normalisation_method The name of the normalization method to apply (e.g., 'quantile').
#' @return The S4 object with normalized data.
#' @export
setGeneric(name="normaliseBetweenSamples"
           , def=function( theObject, normalisation_method = NULL) {
             standardGeneric("normaliseBetweenSamples")
           }
           , signature=c("theObject", "normalisation_method"))

#' @name pearsonCorForSamplePairs
#' @rdname pearsonCorForSamplePairs
#' @title Calculate Pearson Correlation for Sample Pairs
#' @description A generic function to calculate pairwise Pearson correlations between all samples.
#' @param theObject The S4 object.
#' @param tech_rep_remove_regex A regex to process sample names.
#' @param correlation_group Optional. A column in the design matrix to group samples.
#' @return A data frame of pairwise correlation values.
#' @export
setGeneric(name="pearsonCorForSamplePairs"
           , def=function( theObject,   tech_rep_remove_regex = NULL, correlation_group = NA ) {
             standardGeneric("pearsonCorForSamplePairs")
           }
           , signature=c("theObject", "tech_rep_remove_regex", "correlation_group"))

#' @name getNegCtrlProtAnova
#' @rdname getNegCtrlProtAnova
#' @title Identify Negative Control Proteins using ANOVA
#' @description A generic function to identify a set of negative control proteins
#' (i.e., proteins that do not change across conditions) using an ANOVA-based approach.
#' @param theObject The S4 object.
#' @param ruv_grouping_variable The column in the design matrix defining experimental groups for ANOVA.
#' @param percentage_as_neg_ctrl The percentage of proteins to select as negative controls.
#' @param num_neg_ctrl The absolute number of proteins to select as negative controls.
#' @param ruv_qval_cutoff The q-value threshold for determining non-changing proteins.
#' @param ruv_fdr_method The FDR method for p-value adjustment.
#' @return A vector of protein identifiers to be used as negative controls.
#' @export
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

#' @name getLowCoefficientOfVariationProteins
#' @rdname getLowCoefficientOfVariationProteins
#' @title Identify Negative Control Proteins by Low CV
#' @description A generic function to identify negative control proteins based on their
#' low coefficient of variation (CV) across samples.
#' @param theObject The S4 object.
#' @param percentage_as_neg_ctrl The percentage of proteins with the lowest CV to select.
#' @param num_neg_ctrl The absolute number of proteins with the lowest CV to select.
#' @return A vector of protein identifiers to be used as negative controls.
#' @export
setGeneric(name="getLowCoefficientOfVariationProteins"
           , def=function( theObject
                           , percentage_as_neg_ctrl = NULL
                           , num_neg_ctrl = NULL ) {
             standardGeneric("getLowCoefficientOfVariationProteins")
           }
           , signature=c("theObject", "percentage_as_neg_ctrl", "num_neg_ctrl"))

#' @name ruvCancor
#' @rdname ruvCancor
#' @title Apply RUV-Cancor Normalization
#' @description A generic function to apply the RUV-Cancor normalization method to
#' remove unwanted variation.
#' @param theObject The S4 object.
#' @param ctrl A vector of negative control protein/feature identifiers.
#' @param num_components_to_impute The number of unwanted variation factors to remove.
#' @param ruv_grouping_variable The column in the design matrix defining experimental groups.
#' @return The S4 object with RUV-Cancor normalized data.
#' @export
setGeneric(name="ruvCancor"
           , def=function( theObject, ctrl= NULL, num_components_to_impute=NULL, ruv_grouping_variable = NULL ) {
             standardGeneric("ruvCancor")
           }
           , signature=c("theObject", "ctrl", "num_components_to_impute", "ruv_grouping_variable"))

#' @name getRuvIIIReplicateMatrix
#' @rdname getRuvIIIReplicateMatrix
#' @title Get RUV-III Replicate Matrix
#' @description A generic function to create the replicate matrix required for
#' RUV-III normalization.
#' @param theObject The S4 object.
#' @param ruv_grouping_variable The column in the design matrix defining replicate groups.
#' @return A replicate matrix for use with RUV-III.
#' @export
setGeneric(name="getRuvIIIReplicateMatrix"
           , def=function( theObject,  ruv_grouping_variable = NULL) {
             standardGeneric("getRuvIIIReplicateMatrix")
           }
           , signature=c("theObject", "ruv_grouping_variable"))

#' @name ruvIII_C_Varying
#' @rdname ruvIII_C_Varying
#' @title Apply RUV-III Normalization
#' @description A generic function to apply the RUV-III normalization method using
#' a varying number of unwanted variation factors (k).
#' @param theObject The S4 object.
#' @param ruv_grouping_variable The column in the design matrix defining replicate groups.
#' @param ruv_number_k The number of unwanted variation factors (k) to remove.
#' @param ctrl A vector of negative control protein/feature identifiers.
#' @return The S4 object with RUV-III normalized data.
#' @export
setGeneric(name="ruvIII_C_Varying"
           , def=function( theObject, ruv_grouping_variable = NULL, ruv_number_k = NULL, ctrl = NULL)  {
             standardGeneric("ruvIII_C_Varying")
           }
           , signature=c("theObject", "ruv_grouping_variable", "ruv_number_k", "ctrl"))

#' @name removeRowsWithMissingValuesPercent
#' @rdname removeRowsWithMissingValuesPercent
#' @title Filter Rows by Missing Value Percentage
#' @description A generic function to filter rows (proteins/peptides) based on the
#' percentage of missing values within experimental groups.
#' @param theObject The S4 object.
#' @param ruv_grouping_variable The column in the design matrix for group-wise filtering.
#' @param groupwise_percentage_cutoff The maximum allowed percentage of missing values within any single group.
#' @param max_groups_percentage_cutoff The maximum number of groups that can exceed the `groupwise_percentage_cutoff`.
#' @param proteins_intensity_cutoff_percentile An intensity percentile to define what is considered a missing value.
#' @return The filtered S4 object.
#' @export
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

#' @name averageTechReps
#' @rdname averageTechReps
#' @title Average Technical Replicates
#' @description A generic function to average the quantification values of technical replicates.
#' @param theObject The S4 object.
#' @param design_matrix_columns Columns in the design matrix that identify the technical replicates.
#' @return The S4 object with technical replicates averaged.
#' @export
setGeneric(name="averageTechReps"
           , def=function( theObject, design_matrix_columns ) {
             standardGeneric("averageTechReps")
           }
           , signature=c("theObject", "design_matrix_columns" ))

#' @name preservePeptideNaValues
#' @rdname preservePeptideNaValues
#' @title Preserve Peptide NA Values in Protein Data
#' @description After protein-level operations (like normalization or imputation), this
#' generic function re-introduces NA values from the original peptide data where appropriate.
#' @param peptide_obj The S4 object containing the original peptide-level data.
#' @param protein_obj The S4 object containing the processed protein-level data.
#' @return The protein-level S4 object with NA values from peptide data restored.
#' @export
setGeneric(name="preservePeptideNaValues"
           , def=function( peptide_obj, protein_obj)  {
             standardGeneric("preservePeptideNaValues")
           }
           , signature=c("peptide_obj", "protein_obj" ))

#' @name chooseBestProteinAccession
#' @rdname chooseBestProteinAccession
#' @title Choose the Best Protein Accession
#' @description A generic function to resolve protein ambiguity by selecting a single
#' representative protein accession from a group of accessions.
#' @param theObject The S4 object.
#' @param delim The delimiter used to separate protein accessions.
#' @param seqinr_obj A sequence database object used for selection criteria (e.g., length).
#' @param seqinr_accession_column The column in `seqinr_obj` with accession numbers.
#' @param replace_zero_with_na Whether to replace zero values with NA before aggregation.
#' @param aggregation_method The method to aggregate data for duplicate accessions (e.g., 'sum').
#' @return The S4 object with resolved protein accessions.
#' @export
setGeneric(name="chooseBestProteinAccession"
           , def=function(theObject, delim=NULL, seqinr_obj=NULL, seqinr_accession_column=NULL, replace_zero_with_na = NULL, aggregation_method = NULL) {
             standardGeneric("chooseBestProteinAccession")
           }
           , signature=c("theObject", "delim", "seqinr_obj", "seqinr_accession_column"))

#' @name chooseBestProteinAccessionSumDuplicates
#' @rdname chooseBestProteinAccessionSumDuplicates
#' @title Choose Best Protein Accession and Sum Duplicates
#' @description A generic function similar to `chooseBestProteinAccession` but with a specific
#' emphasis on summing duplicate entries after selection.
#' @param theObject The S4 object.
#' @param delim The delimiter used to separate protein accessions.
#' @param seqinr_obj A sequence database object.
#' @param seqinr_accession_column The column in `seqinr_obj` with accession numbers.
#' @param replace_zero_with_na Whether to replace zero values with NA.
#' @param aggregation_method The aggregation method, typically 'sum'.
#' @return The S4 object with resolved and aggregated protein accessions.
#' @export
setGeneric(name="chooseBestProteinAccessionSumDuplicates" # Note: This might be redundant or specific, verify usage
           , def=function(theObject, delim=NULL, seqinr_obj=NULL, seqinr_accession_column=NULL, replace_zero_with_na = NULL, aggregation_method = NULL) {
             standardGeneric("chooseBestProteinAccessionSumDuplicates")
           }
           , signature=c("theObject", "delim", "seqinr_obj", "seqinr_accession_column"))

#' @name filterSamplesByProteinCorrelationThreshold
#' @rdname filterSamplesByProteinCorrelationThreshold
#' @title Filter Samples by Protein Correlation Threshold
#' @description A generic function to remove samples that have a low correlation
#' with other samples in their group.
#' @param theObject The S4 object.
#' @param threshold The minimum Pearson correlation threshold.
#' @param correlation_group The column in the design matrix defining sample groups.
#' @param tech_rep_remove_regex A regex to process sample names.
#' @return The filtered S4 object.
#' @export
setGeneric(name="filterSamplesByProteinCorrelationThreshold"
           , def=function( theObject, threshold = NULL, correlation_group = NULL, tech_rep_remove_regex = NULL) {
             standardGeneric("filterSamplesByProteinCorrelationThreshold")
           }
           , signature=c("theObject", "threshold", "correlation_group", "tech_rep_remove_regex"))

#' @name plotDensity
#' @rdname plotDensity
#' @title Plot Sample Density
#' @description A generic function to generate density plots for each sample, often used
#' to visualize data distributions before and after normalization.
#' @param theObject The S4 object.
#' @param grouping_variable The column in the design matrix for coloring density lines.
#' @param title The plot title.
#' @param font_size The font size for plot text.
#' @return A ggplot object representing the density plot.
#' @export
setGeneric(name="plotDensity"
           , def=function( theObject, grouping_variable, title = "", font_size = 8) {
             standardGeneric("plotDensity")
           }
           , signature=c("theObject", "grouping_variable", "title", "font_size"))

#' @name plotDensityList
#' @rdname plotDensityList
#' @title Plot Density for a List of Grouping Variables
#' @description A generic function to generate a list of density plots, each colored
#' by a different grouping variable.
#' @param theObject The S4 object.
#' @param grouping_variables_list A character vector of column names for coloring plots.
#' @param title The base plot title.
#' @param font_size The font size for plot text.
#' @return A list of ggplot objects.
#' @export
setGeneric(name="plotDensityList"
           , def=function( theObject, grouping_variables_list, title = "", font_size = 8) {
             standardGeneric("plotDensityList")
           }
           , signature=c("theObject", "grouping_variables_list", "title", "font_size"))


# --- From peptideVsSamplesS4Objects.R ---

#' @name cleanDesignMatrixPeptide
#' @rdname cleanDesignMatrixPeptide
#' @title Clean Peptide Design Matrix
#' @description A generic function to clean the design matrix, specifically for
#' peptide-level data.
#' @param theObject The S4 object containing peptide data.
#' @return The modified S4 object with a cleaned design matrix.
#' @export
setGeneric(name ="cleanDesignMatrixPeptide"
           , def=function( theObject) {
             standardGeneric("cleanDesignMatrixPeptide")
           })

#' @name srlQvalueProteotypicPeptideClean
#' @rdname srlQvalueProteotypicPeptideClean
#' @title Filter Peptides by Q-value and Proteotypic Status
#' @description A generic function to filter peptides based on Spectronaut's
#' q-value (SRL) and optionally select only proteotypic peptides.
#' @param theObject The S4 object.
#' @param qvalue_threshold The q-value threshold for keeping peptides.
#' @param global_qvalue_threshold The global q-value threshold.
#' @param choose_only_proteotypic_peptide Logical, whether to keep only proteotypic peptides.
#' @param input_matrix_column_ids Column identifiers in the input matrix.
#' @return The filtered S4 object.
#' @export
setGeneric(name="srlQvalueProteotypicPeptideClean"
           , def=function( theObject, qvalue_threshold = NULL, global_qvalue_threshold = NULL, choose_only_proteotypic_peptide = NULL, input_matrix_column_ids = NULL) {
             standardGeneric("srlQvalueProteotypicPeptideClean")
           }
           , signature = c("theObject", "qvalue_threshold", "global_qvalue_threshold", "choose_only_proteotypic_peptide", "input_matrix_column_ids") )

#' @name rollUpPrecursorToPeptide
#' @rdname rollUpPrecursorToPeptide
#' @title Roll Up Precursor Data to Peptide Level
#' @description A generic function to aggregate precursor-level quantification data
#' to the peptide level.
#' @param theObject The S4 object with precursor data.
#' @param core_utilisation Number of CPU cores for parallel processing.
#' @return The S4 object with data aggregated to the peptide level.
#' @export
setGeneric(name="rollUpPrecursorToPeptide"
           , def=function( theObject, core_utilisation = NULL) {
             standardGeneric("rollUpPrecursorToPeptide")
           }
           , signature=c("theObject", "core_utilisation"))

#' @name peptideIntensityFiltering
#' @rdname peptideIntensityFiltering
#' @title Filter Peptides Based on Intensity
#' @description A generic function to filter low-abundance peptides based on intensity
#' across samples.
#' @param theObject The S4 object with peptide data.
#' @param peptides_intensity_cutoff_percentile The percentile for the intensity cutoff.
#' @param peptides_proportion_of_samples_below_cutoff The proportion of samples threshold.
#' @param core_utilisation Number of CPU cores for parallel processing.
#' @return The filtered S4 object.
#' @export
setGeneric(name="peptideIntensityFiltering"
           , def=function( theObject, peptides_intensity_cutoff_percentile = NULL, peptides_proportion_of_samples_below_cutoff = NULL, core_utilisation = NULL) {
             standardGeneric("peptideIntensityFiltering")
           }
           , signature=c("theObject", "peptides_intensity_cutoff_percentile", "peptides_proportion_of_samples_below_cutoff", "core_utilisation"))

#' @name removePeptidesWithMissingValuesPercent
#' @rdname removePeptidesWithMissingValuesPercent
#' @title Filter Peptides by Missing Value Percentage
#' @description A generic function to filter peptides based on the percentage of missing
#' values within experimental groups.
#' @param theObject The S4 object.
#' @param grouping_variable The column for group-wise filtering.
#' @param groupwise_percentage_cutoff The allowed missing value percentage within a group.
#' @param max_groups_percentage_cutoff The max number of groups exceeding the cutoff.
#' @param peptides_intensity_cutoff_percentile An intensity percentile to define missingness.
#' @return The filtered S4 object.
#' @export
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

#' @name filterMinNumPeptidesPerProtein
#' @rdname filterMinNumPeptidesPerProtein
#' @title Filter by Minimum Number of Peptides per Protein
#' @description A generic function to ensure that each protein is supported by a minimum
#' number of peptides or peptidoforms.
#' @param theObject The S4 object.
#' @param num_peptides_per_protein_thresh The minimum number of peptides required per protein.
#' @param num_peptidoforms_per_protein_thresh The minimum number of peptidoforms required per protein.
#' @param core_utilisation Number of CPU cores for parallel processing.
#' @return The filtered S4 object.
#' @export
setGeneric(name="filterMinNumPeptidesPerProtein"
           , def=function( theObject
                           , num_peptides_per_protein_thresh = NULL
                           , num_peptidoforms_per_protein_thresh = NULL
                           , core_utilisation = NULL ) {
             standardGeneric("filterMinNumPeptidesPerProtein")
           }
           , signature=c("theObject"
                         , "num_peptides_per_protein_thresh"
                         , "num_peptidoforms_per_protein_thresh"
                         , "core_utilisation"))

#' @name filterMinNumPeptidesPerSample
#' @rdname filterMinNumPeptidesPerSample
#' @title Filter by Minimum Number of Peptides per Sample
#' @description A generic function to remove samples that have too few identified peptides.
#' @param theObject The S4 object.
#' @param peptides_per_sample_cutoff The minimum number of peptides required per sample.
#' @param core_utilisation Number of CPU cores for parallel processing.
#' @param inclusion_list A list of samples to always keep, regardless of peptide count.
#' @return The filtered S4 object.
#' @export
setGeneric( name="filterMinNumPeptidesPerSample"
            , def=function( theObject, peptides_per_sample_cutoff = NULL, core_utilisation = NULL, inclusion_list = NULL) {
              standardGeneric("filterMinNumPeptidesPerSample")
           }
           , signature=c("theObject", "peptides_per_sample_cutoff", "core_utilisation", "inclusion_list" ))

#' @name removePeptidesWithOnlyOneReplicate
#' @rdname removePeptidesWithOnlyOneReplicate
#' @title Remove Peptides with Only a Single Replicate
#' @description A generic function to filter out peptides that are identified in only
#' a single replicate within their experimental group.
#' @param theObject The S4 object.
#' @param replicate_group_column The column identifying replicate groups.
#' @param core_utilisation Number of CPU cores for parallel processing.
#' @return The filtered S4 object.
#' @export
setGeneric( name="removePeptidesWithOnlyOneReplicate"
            , def=function( theObject, replicate_group_column = NULL, core_utilisation = NULL) {
              standardGeneric("removePeptidesWithOnlyOneReplicate")
            }
            , signature=c("theObject", "replicate_group_column", "core_utilisation" ))

#' @name plotPeptidesProteinsCountsPerSample
#' @rdname plotPeptidesProteinsCountsPerSample
#' @title Plot Peptide and Protein Counts per Sample
#' @description A generic function to generate a bar plot showing the number of
#' identified peptides and proteins for each sample.
#' @param theObject The S4 object.
#' @return A ggplot object.
#' @export
setGeneric( name="plotPeptidesProteinsCountsPerSample"
            , def=function( theObject ) {
              standardGeneric("plotPeptidesProteinsCountsPerSample")
            }
            , signature=c("theObject" ))

#' @name peptideMissingValueImputation
#' @rdname peptideMissingValueImputation
#' @title Impute Missing Peptide Values
#' @description A generic function to perform missing value imputation at the peptide level.
#' @param theObject The S4 object.
#' @param imputed_value_column The column to store imputed values.
#' @param proportion_missing_values The proportion of missing values to guide imputation.
#' @param core_utilisation Number of CPU cores for parallel processing.
#' @return The S4 object with imputed peptide data.
#' @export
setGeneric( name="peptideMissingValueImputation"
            , def=function( theObject,  imputed_value_column = NULL, proportion_missing_values = NULL, core_utilisation = NULL) {
              standardGeneric("peptideMissingValueImputation")
            }
            , signature=c("theObject", "imputed_value_column", "proportion_missing_values", "core_utilisation" ))


# --- From metabolite_qc.R ---

#' @name metaboliteIntensityFiltering
#' @rdname metaboliteIntensityFiltering
#' @title Filter Metabolites Based on Intensity
#' @description A generic function to filter low-abundance metabolites based on
#' intensity across samples.
#' @param theObject The S4 object with metabolite data.
#' @param metabolites_intensity_cutoff_percentile The percentile for the intensity cutoff.
#' @param metabolites_proportion_of_samples_below_cutoff The proportion of samples threshold.
#' @return The filtered S4 object.
#' @export
setGeneric(name="metaboliteIntensityFiltering"
           , def=function( theObject, metabolites_intensity_cutoff_percentile = NULL, metabolites_proportion_of_samples_below_cutoff = NULL) {
             standardGeneric("metaboliteIntensityFiltering")
           }
           , signature=c("theObject")) # Dispatch only on the object

#' @name resolveDuplicateFeatures
#' @rdname resolveDuplicateFeatures
#' @title Resolve Duplicate Metabolite Features
#' @description A generic function to handle duplicate metabolite feature identifiers,
#' for instance by aggregating them or selecting the most intense one.
#' @param theObject The S4 object.
#' @param itsd_pattern_columns Columns related to internal standards that might be
#' used in the resolution logic.
#' @return The S4 object with duplicate features resolved.
#' @export
setGeneric(name = "resolveDuplicateFeatures",
           def = function(theObject, itsd_pattern_columns = NULL) {
             standardGeneric("resolveDuplicateFeatures")
           }
)

# --- From metabolite_normalization.R ---

#' @name logTransformAssays
#' @rdname logTransformAssays
#' @title Log-Transform Assay Data
#' @description A generic function to apply a log transformation (e.g., log2) to the
#' quantitative assay data in an object.
#' @param theObject The S4 object.
#' @param offset A small value to add to the data before transformation to avoid log(0).
#' @param ... Additional arguments passed to methods.
#' @return The S4 object with log-transformed data.
#' @export
setGeneric(name = "logTransformAssays",
           def = function(theObject, offset = 1, ...) {
             standardGeneric("logTransformAssays")
           }
)

#' @name normaliseUntransformedData
#' @rdname normaliseUntransformedData
#' @title Normalize Untransformed Data
#' @description A generic function to apply normalization methods, such as ITSD
#' (Internal Standard) normalization, to the untransformed data.
#' @param theObject The S4 object.
#' @param method The normalization method to apply (e.g., "ITSD").
#' @param ... Additional arguments passed to methods.
#' @return The S4 object with normalized data.
#' @export
setGeneric(name = "normaliseUntransformedData",
           def = function(theObject, method = "ITSD", ...) {
             standardGeneric("normaliseUntransformedData")
           }
)

# --- From metaboliteVsSamplesS4Objects.R ---
# (Only generics defined uniquely in this file)

#' @name getNegCtrlMetabAnova
#' @rdname getNegCtrlMetabAnova
#' @title Identify Negative Control Metabolites using ANOVA
#' @description A generic function to identify a set of negative control metabolites
#' (i.e., metabolites that do not change across conditions) using an ANOVA-based approach.
#' @param theObject The S4 object with metabolite data.
#' @param ruv_grouping_variable The column in the design matrix defining experimental groups for ANOVA.
#' @param percentage_as_neg_ctrl The percentage of metabolites to select as negative controls.
#' @param num_neg_ctrl The absolute number of metabolites to select as negative controls.
#' @param ruv_qval_cutoff The q-value threshold for determining non-changing metabolites.
#' @param ruv_fdr_method The FDR method for p-value adjustment.
#' @return A vector of metabolite identifiers to be used as negative controls.
#' @export
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