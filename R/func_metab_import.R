# ============================================================================
# func_metab_import.R
# ============================================================================
# Purpose: Metabolomics data import functions
# 
# This file contains functions for importing metabolomics data from various
# platforms and formats. Functions in this file are used by metabolomics
# import modules and related workflows.
#
# Functions to extract here:
# - createMetaboliteAssayData(): Creates MetaboliteAssayData S4 object
# - getMetaboliteQuantData(): Extracts quantitative data from assay tibbles
# - Additional metabolomics import helper functions
#
# Dependencies:
# - dplyr, tidyr
# - func_general_s4_objects.R (for S4 class definitions)
# ============================================================================

# TODO: Extract the following functions from their current locations:

# Function 1: createMetaboliteAssayData()
# Current location: R/metaboliteVsSamplesS4Objects.R
# Description: Creates a MetaboliteAssayData S4 object from assay data
# createMetaboliteAssayData <- function(...) {
#   # Extract from R/metaboliteVsSamplesS4Objects.R
# }

# Function 2: getMetaboliteQuantData()
# Current location: R/QC_visualisation.R
# Description: Extracts quantitative data and sample names from assay tibble
# getMetaboliteQuantData <- function(assay_data) {
#   # Extract from R/QC_visualisation.R
# }


# ----------------------------------------------------------------------------
# createMetaboliteAssayData
# ----------------------------------------------------------------------------
#' Create MetaboliteAssayData Object
#'
#' Constructor function for the MetaboliteAssayData class.
#'
#' @param metabolite_data Named list of data frames (assays).
#' @param design_matrix Experimental design data frame.
#' @param metabolite_id_column Name of the **primary feature ID** column within assays (e.g., `"database_identifier"`).
#' @param annotation_id_column Name of the **annotation ID** column (e.g., `"metabolite_identification"`).
#' @param sample_id Name of the sample ID column in design_matrix and assays.
#' @param group_id Name of the group column in design_matrix.
#' @param technical_replicate_id Name of the technical replicate column in design_matrix (use NA_character_ if none).
#' @param database_identifier_type Type of identifier in the `annotation_id_column` (e.g., `"Mixed_CHEBI_Unknown"`).
#' @param internal_standard_regex Regex to identify internal standards. Use `NA_character_` or `""` if none.
#' @param args List of arguments (e.g., from config).
#'
#' @return A MetaboliteAssayData object.
#' @export
#' @examples
#' \dontrun{
#' # Assuming lcms_pos_df, lcms_neg_df, gcms_df are data frames
#' # with 'Metabolite' as ID column and samples as other columns
#' # Assuming design_df has 'SampleID', 'Group', 'Replicate' columns
#' assays_list <- list(
#'     LCMS_Pos = lcms_pos_df,
#'     LCMS_Neg = lcms_neg_df,
#'     GCMS = gcms_df
#' )
#' config <- list(...) # Your config list
#'
#' met_assay_obj <- createMetaboliteAssayData(
#'     metabolite_data = assays_list,
#'     design_matrix = design_df,
#'     metabolite_id_column = "Metabolite",
#'     sample_id = "SampleID",
#'     group_id = "Group",
#'     technical_replicate_id = "Replicate",
#'     database_identifier_type = "InternalName",
#'     internal_standard_regex = "^IS_",
#'     args = config
#' )
#' }
createMetaboliteAssayData <- function(
    metabolite_data,
    design_matrix,
    metabolite_id_column = "database_identifier",
    annotation_id_column = "metabolite_identification",
    sample_id = "Sample_ID",
    group_id = "group",
    technical_replicate_id = NA_character_,
    database_identifier_type = "Unknown",
    internal_standard_regex = NA_character_,
    args = list()) {
    # Perform basic checks before creating the object
    stopifnot(is.list(metabolite_data))
    stopifnot(all(sapply(metabolite_data, is.data.frame)))
    stopifnot(is.data.frame(design_matrix))
    # Add more checks as needed...

    obj <- new("MetaboliteAssayData",
        metabolite_data = metabolite_data,
        metabolite_id_column = metabolite_id_column,
        annotation_id_column = annotation_id_column,
        database_identifier_type = database_identifier_type,
        internal_standard_regex = internal_standard_regex,
        design_matrix = design_matrix,
        sample_id = sample_id,
        group_id = group_id,
        technical_replicate_id = technical_replicate_id,
        args = args
    )
    # Validity check is automatically called by 'new'
    return(obj)
}


# ----------------------------------------------------------------------------
# getMetaboliteQuantData
# ----------------------------------------------------------------------------
#' @title Extract Quantitative Data and Sample Names from Assay Tibble
#' @description Separates annotation columns from quantitative data columns
#'              in a metabolomics assay tibble.
#'
#' @param assay_data A tibble/data.frame representing one metabolomics assay,
#'                   with metabolite annotations and sample intensity columns.
#' @param sample_id_col The name of the column in the design matrix that contains
#'                      the sample identifiers (which match column names in assay_data).
#'                      ***Note: This argument is currently unused but kept for consistency.
#'                      The function infers sample columns based on numeric type.***
#'
#' @return A list containing:
#'         - `quant_data`: A data frame with only the quantitative (numeric) columns.
#'         - `sample_names`: A character vector of the sample column names.
#'         - `annotation_data`: A data frame with the non-numeric annotation columns.
#'
#' @importFrom dplyr select where
#' @keywords internal
#' @noRd
#' @export 
getMetaboliteQuantData <- function(assay_data) {
    # Identify quantitative (numeric) columns - assumes sample columns are numeric
    quant_cols <- sapply(assay_data, is.numeric)
    quant_data <- assay_data[, quant_cols, drop = FALSE]
    sample_names <- colnames(quant_data)

    # Identify annotation (non-numeric) columns
    annotation_data <- assay_data[, !quant_cols, drop = FALSE]

    return(list(
        quant_data = as.data.frame(quant_data), # Convert tibble to df if needed downstream
        sample_names = sample_names,
        annotation_data = as.data.frame(annotation_data)
    ))
}

