# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

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
#' @param sample_columns Optional character vector of explicit sample column names.
#'                       When provided, these columns are used directly instead of
#'                       guessing based on numeric type. This is the preferred method
#'                       as MS-DIAL data contains many numeric annotation columns
#'                       (scores, m/z, RT) that are NOT sample data.
#'
#' @return A list containing:
#'         - `quant_data`: A data frame with only the quantitative (sample) columns.
#'         - `sample_names`: A character vector of the sample column names.
#'         - `annotation_data`: A data frame with the non-sample annotation columns.
#'
#' @importFrom dplyr select where
#' @keywords internal
#' @noRd
#' @export 
getMetaboliteQuantData <- function(assay_data, sample_columns = NULL) {
    if (!is.null(sample_columns) && length(sample_columns) > 0) {
        # Use explicit sample columns (preferred - avoids including numeric annotation cols)
        valid_cols <- intersect(sample_columns, colnames(assay_data))
        if (length(valid_cols) == 0) {
            warning("None of the provided sample_columns exist in assay_data. Falling back to numeric detection.")
            quant_cols <- sapply(assay_data, is.numeric)
            quant_data <- assay_data[, quant_cols, drop = FALSE]
            sample_names <- colnames(quant_data)
        } else {
            quant_data <- assay_data[, valid_cols, drop = FALSE]
            sample_names <- valid_cols
        }
    } else {
        # Fallback: guess by numeric type (unreliable for MS-DIAL data with numeric scores)
        quant_cols <- sapply(assay_data, is.numeric)
        quant_data <- assay_data[, quant_cols, drop = FALSE]
        sample_names <- colnames(quant_data)
    }

    # Identify annotation columns (everything NOT in sample_names)
    annotation_data <- assay_data[, setdiff(colnames(assay_data), sample_names), drop = FALSE]

    return(list(
        quant_data = as.data.frame(quant_data)
        , sample_names = sample_names
        , annotation_data = as.data.frame(annotation_data)
    ))
}








# ----------------------------------------------------------------------------
# findMatchingColumn
# ----------------------------------------------------------------------------
#' @title Find Matching Column in Headers (Case-Insensitive)
#' @description Searches for a column name in headers using case-insensitive matching.
#'              Used for auto-populating column mapping inputs.
#'
#' @param headers Character vector of column headers from the data file.
#' @param candidates Character vector of candidate column names to search for.
#'
#' @return The first matching header name (preserving original case), or NULL if no match.
#'
#' @export
# findMetabMatchingColumn <- function(headers, candidates) {
#     # Deprecated in favor of shared findMatchingColumn
#     findMatchingColumn(headers, candidates)
# }


# ----------------------------------------------------------------------------
# importMSDIALData
# ----------------------------------------------------------------------------
#' @title Import MS-DIAL Data (Alias)
#' @description Alias for importMetabMSDIALData for backward compatibility.
#' @export
importMSDIALData <- function(...) {
    importMetabMSDIALData(...)
}


# ----------------------------------------------------------------------------
# importMSDIALData (Implementation)
# ----------------------------------------------------------------------------
#' @title Import MS-DIAL Data
#' @description Parses MS-DIAL exported data files (height or area matrix).
#'              Handles the wide format with metabolite annotations and sample columns.
#'
#' @param filepath Path to the MS-DIAL export file (TSV or CSV).
#' @param metabolite_id_column Column to use as metabolite ID (default: auto-detect).
#' @param annotation_column Column containing metabolite annotations (default: auto-detect).
#' @param skip_rows Number of header rows to skip (default: 0).
#'
#' @return A list containing:
#'   - `data`: The imported data frame
#'   - `headers`: Column headers
#'   - `detected_columns`: List of auto-detected column assignments
#'   - `sample_columns`: Character vector of sample column names
#'   - `annotation_columns`: Character vector of annotation column names
#'
#' @importFrom vroom vroom
#' @importFrom logger log_info log_error
#' @export
importMetabMSDIALData <- function(
    filepath
    , metabolite_id_column = NULL
    , annotation_column = NULL
    , skip_rows = 0
) {
    logger::log_info(paste("Starting MS-DIAL import from:", filepath))
    
    # Check file exists
    if (!file.exists(filepath)) {
        stop("File not found: ", filepath)
    }
    
    # Detect delimiter from file extension
    file_ext <- tolower(tools::file_ext(filepath))
    delim <- if (file_ext == "csv") "," else "\t"
    
    # Read data
    data <- tryCatch({
        vroom::vroom(
            filepath
            , delim = delim
            , skip = skip_rows
            , show_col_types = FALSE
        )
    }, error = function(e) {
        stop("Failed to read file: ", e$message)
    })
    
    headers <- names(data)
    logger::log_info(sprintf("Read %d rows and %d columns", nrow(data), ncol(data)))
    
    # Auto-detect column mappings using defaults
    defaults <- getMetabolomicsColumnDefaults("msdial")
    
    detected_columns <- list(
        metabolite_id = if (!is.null(metabolite_id_column)) {
            metabolite_id_column
        } else {
            findMetabMatchingColumn(headers, defaults$metabolite_id)
        }
        , annotation = if (!is.null(annotation_column)) {
            annotation_column
        } else {
            findMetabMatchingColumn(headers, defaults$annotation)
        }
        , rt = findMetabMatchingColumn(headers, defaults$rt)
        , mz = findMetabMatchingColumn(headers, defaults$mz)
        , adduct = findMetabMatchingColumn(headers, defaults$adduct)
    )
    
    # Identify sample columns (numeric columns not in known annotation columns)
    # These patterns match MS-DIAL annotation columns that should NOT be treated as samples
    known_annotation_patterns <- c(
        # ID and name columns
        "peak id", "alignment id", "^name$", "metabolite name"
        # Scan and RT columns
        , "^scan$", "rt left", "^rt \\(min\\)$", "rt right", "average rt", "reference rt"
        # Mass columns
        , "precursor m/z", "average mz", "reference m/z", "model masses"
        # Intensity metrics (these are per-feature, not per-sample)
        , "^height$", "^area$"
        # Adduct and isotope
        , "^adduct$", "adduct type", "^isotope$"
        # Identification columns
        , "formula", "ontology", "inchikey", "smiles", "^comment$"
        # Annotation tags and matching
        , "annotation tag", "rt matched", "m/z matched", "ms/ms matched"
        , "rt similarity", "m/z similarity"
        # Scoring columns
        , "simple dot product", "weighted dot product", "reverse dot product"
        , "matched peaks", "total score", "s/n", "snratio"
        # Spectral data
        , "ms1 isotopes", "msms spectrum", "ms/ms"
        # Legacy patterns
        , "post curation", "fill", "reference", "peak shape", "spectral"
        , "^row$", "^class$"
    )
    
    sample_columns <- character(0)
    annotation_columns <- character(0)
    
    for (col in headers) {
        col_lower <- tolower(col)
        is_annotation <- any(sapply(known_annotation_patterns, function(p) grepl(p, col_lower)))
        
        if (!is_annotation && is.numeric(data[[col]])) {
            sample_columns <- c(sample_columns, col)
        } else {
            annotation_columns <- c(annotation_columns, col)
        }
    }
    
    logger::log_info(sprintf(
        "Detected %d sample columns and %d annotation columns"
        , length(sample_columns)
        , length(annotation_columns)
    ))
    
    return(list(
        data = as.data.frame(data)
        , headers = headers
        , detected_columns = detected_columns
        , sample_columns = sample_columns
        , annotation_columns = annotation_columns
        , format = "msdial"
        , is_pattern = defaults$is_pattern
    ))
}




