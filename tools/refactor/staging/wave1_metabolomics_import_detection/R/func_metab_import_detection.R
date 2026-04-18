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
        format = best_format,
        confidence = best_score,
        all_scores = scores
    ))
}

# ----------------------------------------------------------------------------
# findMatchingColumn (Alias)
# ----------------------------------------------------------------------------
#' @title Find Matching Column (Alias)
#' @description Alias for findMatchingColumn for backward compatibility.
#' @export
findMetabMatchingColumn <- function(headers, candidates) {
    if (!requireNamespace("MultiScholaR", quietly = TRUE)) {
        # Fallback if package not fully loaded
        return(findMatchingColumn(headers, candidates))
    }
    MultiScholaR::findMatchingColumn(headers, candidates)
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
# validateColumnMapping (Alias)
# ----------------------------------------------------------------------------
#' @title Validate Column Mapping (Alias)
#' @description Alias for validateMetabColumnMapping for backward compatibility.
#' @export
validateColumnMapping <- function(data, metabolite_id_column, sample_columns) {
    validateMetabColumnMapping(data, metabolite_id_column, sample_columns)
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
validateMetabColumnMapping <- function(data, metabolite_id_column, sample_columns) {
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

