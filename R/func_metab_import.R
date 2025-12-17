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


# ----------------------------------------------------------------------------
# detectMetabolomicsFormat
# ----------------------------------------------------------------------------
#' @title Detect Metabolomics Data Format
#' @description Auto-detect the vendor format of metabolomics data based on
#'              column headers and filename patterns.
#'
#' @param headers Character vector of column headers from the data file.
#' @param filename Name of the data file (optional, used for additional hints).
#'
#' @return A list containing:
#'   - `format`: Character string identifying the format ("msdial", "progenesis", "xcms", "compound_discoverer", "unknown")
#'   - `confidence`: Numeric score (0-1) indicating detection confidence
#'   - `all_scores`: Named numeric vector of scores for all formats
#'
#' @export
detectMetabolomicsFormat <- function(headers, filename = NULL) {
    # Convert headers to lowercase for case-insensitive matching
    headers_lower <- tolower(headers)
    filename_lower <- if (!is.null(filename)) tolower(filename) else ""
    
    # --- MS-DIAL detection ---
    # Supports both older (Alignment ID, Average Rt) and newer (Peak ID, RT (min)) formats
    msdial_markers <- c(
        # Core identifiers (old and new formats)
        "peak id"
        , "alignment id"
        , "name"
        , "metabolite name"
        # Retention time variants
        , "rt (min)"
        , "rt left(min)"
        , "rt right (min)"
        , "average rt(min)"
        , "reference rt"
        # Mass variants
        , "precursor m/z"
        , "average mz"
        , "reference m/z"
        # Intensity
        , "height"
        , "area"
        # Adduct
        , "adduct"
        , "adduct type"
        # Annotation/identification
        , "formula"
        , "ontology"
        , "inchikey"
        , "smiles"
        , "annotation tag (vs1.0)"
        , "annotation tag (vs. ms-dial)"
        # Matching scores
        , "rt matched"
        , "m/z matched"
        , "ms/ms matched"
        , "total score"
        , "s/n"
        , "msms spectrum"
    )
    msdial_found <- sum(msdial_markers %in% headers_lower)
    msdial_score <- msdial_found / min(length(msdial_markers), 10)  # Normalize to top 10 markers
    
    # Filename bonus
    if (grepl("msdial|ms-dial|height|area", filename_lower)) {
        msdial_score <- min(1.0, msdial_score + 0.15)
    }
    
    # --- Progenesis QI detection ---
    progenesis_markers <- c(
        "compound"
        , "neutral mass (da)"
        , "m/z"
        , "charge"
        , "retention time (min)"
        , "chromatographic peak width"
        , "max abundance"
        , "raw abundance"
        , "normalised abundance"
        , "highest mean condition"
        , "score"
        , "fragmentation score"
        , "mass error (ppm)"
    )
    progenesis_found <- sum(progenesis_markers %in% headers_lower)
    progenesis_score <- progenesis_found / min(length(progenesis_markers), 6)
    
    if (grepl("progenesis|qi", filename_lower)) {
        progenesis_score <- min(1.0, progenesis_score + 0.15)
    }
    
    # --- XCMS detection ---
    xcms_markers <- c(
        "featureid"
        , "mz"
        , "mzmin"
        , "mzmax"
        , "rt"
        , "rtmin"
        , "rtmax"
        , "npeaks"
        , "peakidx"
        , "into"
        , "intb"
        , "maxo"
        , "sn"
    )
    xcms_found <- sum(xcms_markers %in% headers_lower)
    xcms_score <- xcms_found / min(length(xcms_markers), 6)
    
    if (grepl("xcms|peak_?list", filename_lower)) {
        xcms_score <- min(1.0, xcms_score + 0.15)
    }
    
    # --- Compound Discoverer detection ---
    cd_markers <- c(
        "compound name"
        , "molecular formula"
        , "molecular weight"
        , "rt [min]"
        , "area"
        , "height"
        , "mz"
        , "delta mass [ppm]"
        , "msi annotation level"
        , "# compounds"
        , "# features"
        , "ms2 spectrum"
        , "best match"
    )
    cd_found <- sum(cd_markers %in% headers_lower)
    cd_score <- cd_found / min(length(cd_markers), 6)
    
    if (grepl("compound.?discoverer|cd_", filename_lower)) {
        cd_score <- min(1.0, cd_score + 0.15)
    }
    
    # Determine best match
    scores <- c(
        msdial = msdial_score
        , progenesis = progenesis_score
        , xcms = xcms_score
        , compound_discoverer = cd_score
    )
    
    best_format <- names(which.max(scores))
    best_score <- max(scores)
    
    # Minimum confidence threshold
    if (best_score < 0.2) {
        best_format <- "unknown"
    }
    
    return(list(
        format = best_format
        , confidence = best_score
        , all_scores = scores
    ))
}


# ----------------------------------------------------------------------------
# getMetabolomicsColumnDefaults
# ----------------------------------------------------------------------------
#' @title Get Default Column Mappings for Metabolomics Formats
#' @description Returns the expected column names for each supported vendor format.
#'              Used for auto-populating column mapping dropdowns.
#'
#' @param format Character string specifying the format ("msdial", "progenesis", "xcms", "compound_discoverer").
#'
#' @return A list containing expected column names:
#'   - `metabolite_id`: Primary identifier column
#'   - `annotation`: Metabolite name/annotation column
#'   - `rt`: Retention time column
#'   - `mz`: Mass-to-charge ratio column
#'   - `adduct`: Adduct type column (if applicable)
#'   - `is_pattern`: Suggested regex for internal standards
#'
#' @export
getMetabolomicsColumnDefaults <- function(format) {
    defaults <- switch(tolower(format)
        , "msdial" = list(
            # Support both old format (Alignment ID) and new format (Peak ID)
            metabolite_id = c("Peak ID", "peak id", "Alignment ID", "alignment id")
            , annotation = c("Name", "name", "Metabolite name", "metabolite name")
            , rt = c("RT (min)", "rt (min)", "Average Rt(min)", "average rt(min)")
            , mz = c("Precursor m/z", "precursor m/z", "Average Mz", "average mz")
            , adduct = c("Adduct", "adduct", "Adduct type", "adduct type")
            , height = c("Height", "height")
            , area = c("Area", "area")
            , formula = c("Formula", "formula")
            , ontology = c("Ontology", "ontology")
            , is_pattern = "^IS_|_d[0-9]+$|ISTD|IS-|w/o MS2:"
        )
        , "progenesis" = list(
            metabolite_id = c("Compound", "compound")
            , annotation = c("Compound", "compound")
            , rt = c("Retention time (min)", "retention time (min)")
            , mz = c("m/z", "M/z")
            , adduct = NULL
            , is_pattern = "^IS_|ISTD"
        )
        , "xcms" = list(
            metabolite_id = c("featureid", "featureId", "FeatureID")
            , annotation = NULL
            , rt = c("rt", "RT")
            , mz = c("mz", "MZ", "m/z")
            , adduct = NULL
            , is_pattern = "^IS_|ISTD"
        )
        , "compound_discoverer" = list(
            metabolite_id = c("Compound Name", "compound name")
            , annotation = c("Compound Name", "compound name")
            , rt = c("RT [min]", "rt [min]")
            , mz = c("MZ", "mz", "m/z")
            , adduct = NULL
            , is_pattern = "^IS_|ISTD|Internal"
        )
        , "custom" = list(
            metabolite_id = NULL
            , annotation = NULL
            , rt = NULL
            , mz = NULL
            , adduct = NULL
            , is_pattern = NA_character_
        )
        , list(
            metabolite_id = NULL
            , annotation = NULL
            , rt = NULL
            , mz = NULL
            , adduct = NULL
            , is_pattern = NA_character_
        )
    )
    
    return(defaults)
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
findMatchingColumn <- function(headers, candidates) {
    if (is.null(candidates) || length(candidates) == 0) {
        return(NULL)
    }
    
    headers_lower <- tolower(headers)
    candidates_lower <- tolower(candidates)
    
    for (cand in candidates_lower) {
        match_idx <- which(headers_lower == cand)
        if (length(match_idx) > 0) {
            return(headers[match_idx[1]])
        }
    }
    
    return(NULL)
}


# ----------------------------------------------------------------------------
# importMSDIALData
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
importMSDIALData <- function(
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
            findMatchingColumn(headers, defaults$metabolite_id)
        }
        , annotation = if (!is.null(annotation_column)) {
            annotation_column
        } else {
            findMatchingColumn(headers, defaults$annotation)
        }
        , rt = findMatchingColumn(headers, defaults$rt)
        , mz = findMatchingColumn(headers, defaults$mz)
        , adduct = findMatchingColumn(headers, defaults$adduct)
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


# ----------------------------------------------------------------------------
# validateColumnMapping
# ----------------------------------------------------------------------------
#' @title Validate Column Mapping for Metabolomics Data
#' @description Checks that required columns exist in the data and returns validation status.
#'
#' @param data Data frame to validate.
#' @param metabolite_id_column Name of the metabolite ID column.
#' @param sample_columns Character vector of sample column names.
#'
#' @return A list containing:
#'   - `valid`: Logical indicating overall validity
#'   - `errors`: Character vector of error messages
#'   - `warnings`: Character vector of warning messages
#'   - `summary`: Summary statistics
#'
#' @export
validateColumnMapping <- function(data, metabolite_id_column, sample_columns) {
    errors <- character(0)
    warnings <- character(0)
    
    # Check metabolite ID column
    if (is.null(metabolite_id_column) || !nzchar(metabolite_id_column)) {
        errors <- c(errors, "Metabolite ID column is not specified")
    } else if (!metabolite_id_column %in% names(data)) {
        errors <- c(errors, paste("Metabolite ID column not found:", metabolite_id_column))
    }
    
    # Check sample columns
    if (length(sample_columns) == 0) {
        errors <- c(errors, "No sample columns specified")
    } else {
        missing_samples <- sample_columns[!sample_columns %in% names(data)]
        if (length(missing_samples) > 0) {
            errors <- c(errors, paste("Sample columns not found:", paste(missing_samples, collapse = ", ")))
        }
    }
    
    # Calculate summary if valid so far
    summary_stats <- list()
    if (length(errors) == 0) {
        summary_stats$n_metabolites <- length(unique(data[[metabolite_id_column]]))
        summary_stats$n_samples <- length(sample_columns)
        
        # Check for duplicates
        n_duplicates <- sum(duplicated(data[[metabolite_id_column]]))
        if (n_duplicates > 0) {
            warnings <- c(warnings, sprintf("%d duplicate metabolite IDs detected", n_duplicates))
        }
        
        # Check for missing values
        sample_data <- data[, sample_columns, drop = FALSE]
        n_missing <- sum(is.na(sample_data) | sample_data == 0)
        total_cells <- nrow(sample_data) * ncol(sample_data)
        pct_missing <- round((n_missing / total_cells) * 100, 1)
        summary_stats$pct_missing <- pct_missing
        
        if (pct_missing > 50) {
            warnings <- c(warnings, sprintf("High proportion of missing values: %.1f%%", pct_missing))
        }
    }
    
    return(list(
        valid = length(errors) == 0
        , errors = errors
        , warnings = warnings
        , summary = summary_stats
    ))
}