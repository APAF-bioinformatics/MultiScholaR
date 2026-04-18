# ----------------------------------------------------------------------------
# importMSDIALData
# ----------------------------------------------------------------------------
#' @title Import MS-DIAL Data (Alias)
#' @description Alias for importLipidMSDIALData for backward compatibility.
#' @export
importMSDIALData <- function(...) {
    importLipidMSDIALData(...)
}

# ----------------------------------------------------------------------------
# importMSDIALData (Implementation)
# ----------------------------------------------------------------------------
#' @title Import MS-DIAL Data
#' @description Parses MS-DIAL exported data files (height or area matrix).
#'              Handles the wide format with lipid annotations and sample columns.
#'
#' @param filepath Path to the MS-DIAL export file (TSV or CSV).
#' @param lipid_id_column Column to use as lipid ID (default: auto-detect).
#' @param annotation_column Column containing lipid annotations (default: auto-detect).
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
importLipidMSDIALData <- function(
  filepath,
  lipid_id_column = NULL,
  annotation_column = NULL,
  skip_rows = 0
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
    data <- tryCatch(
        {
            vroom::vroom(
                filepath,
                delim = delim,
                skip = skip_rows,
                show_col_types = FALSE
            )
        },
        error = function(e) {
            stop("Failed to read file: ", e$message)
        }
    )

    headers <- names(data)
    logger::log_info(sprintf("Read %d rows and %d columns", nrow(data), ncol(data)))

    # Auto-detect column mappings using defaults
    defaults <- getLipidomicsColumnDefaults("msdial")

    detected_columns <- list(
        lipid_id = if (!is.null(lipid_id_column)) {
            lipid_id_column
        } else {
            findLipidMatchingColumn(headers, defaults$lipid_id)
        },
        annotation = if (!is.null(annotation_column)) {
            annotation_column
        } else {
            findLipidMatchingColumn(headers, defaults$annotation)
        },
        rt = findLipidMatchingColumn(headers, defaults$rt),
        mz = findLipidMatchingColumn(headers, defaults$mz),
        adduct = findLipidMatchingColumn(headers, defaults$adduct)
    )

    # Identify sample columns (numeric columns not in known annotation columns)
    # These patterns match MS-DIAL annotation columns that should NOT be treated as samples
    known_annotation_patterns <- c(
        # ID and name columns
        "peak id", "alignment id", "^name$", "lipid name"
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
        , "annotation tag", "rt matched", "m/z matched", "ms/ms matched",
        "rt similarity", "m/z similarity"
        # Scoring columns
        , "simple dot product", "weighted dot product", "reverse dot product",
        "matched peaks", "total score", "s/n", "snratio"
        # Spectral data
        , "ms1 isotopes", "msms spectrum", "ms/ms"
        # Legacy patterns
        , "post curation", "fill", "reference", "peak shape", "spectral",
        "^row$", "^class$"
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
        "Detected %d sample columns and %d annotation columns",
        length(sample_columns),
        length(annotation_columns)
    ))

    return(list(
        data = as.data.frame(data),
        headers = headers,
        detected_columns = detected_columns,
        sample_columns = sample_columns,
        annotation_columns = annotation_columns,
        format = "msdial",
        is_pattern = defaults$is_pattern
    ))
}

# ----------------------------------------------------------------------------
# importLipidSearchData
# ----------------------------------------------------------------------------
#' @title Import LipidSearch Data
#' @description Parses LipidSearch exported data files (Alignment format).
#'
#' @param filepath Path to the LipidSearch export file (CSV or TSV).
#' @param lipid_id_column Column to use as lipid ID (default: auto-detect).
#' @param annotation_column Column containing lipid annotations (default: auto-detect).
#' @param skip_rows Number of header rows to skip (default: 0).
#'
#' @return A list containing the imported data and metadata.
#'
#' @importFrom vroom vroom
#' @importFrom logger log_info log_error
#' @export
importLipidSearchData <- function(
  filepath,
  lipid_id_column = NULL,
  annotation_column = NULL,
  skip_rows = 0
) {
    logger::log_info(paste("Starting LipidSearch import from:", filepath))

    if (!file.exists(filepath)) {
        stop("File not found: ", filepath)
    }

    # LipidSearch usually exports CSV
    file_ext <- tolower(tools::file_ext(filepath))
    delim <- if (file_ext == "tsv" || file_ext == "txt") "\t" else ","

    data <- tryCatch(
        {
            vroom::vroom(
                filepath,
                delim = delim,
                skip = skip_rows,
                show_col_types = FALSE
            )
        },
        error = function(e) {
            stop("Failed to read file: ", e$message)
        }
    )

    headers <- names(data)
    logger::log_info(sprintf("Read %d rows and %d columns", nrow(data), ncol(data)))

    # Auto-detect column mappings
    defaults <- getLipidomicsColumnDefaults("lipidsearch")

    detected_columns <- list(
        lipid_id = if (!is.null(lipid_id_column)) {
            lipid_id_column
        } else {
            findMatchingColumn(headers, defaults$lipid_id)
        },
        annotation = if (!is.null(annotation_column)) {
            annotation_column
        } else {
            findMatchingColumn(headers, defaults$annotation)
        },
        rt = findMatchingColumn(headers, defaults$rt),
        mz = findMatchingColumn(headers, defaults$mz),
        adduct = findMatchingColumn(headers, defaults$adduct)
    )

    # Identify sample columns
    # Exclude known metadata columns
    known_metadata <- c(
        "lipidclass", "lipidname", "lipidion", "lipidgroup",
        "class", "fattyacid", "iontype", "adduct",
        "calcmz", "basert", "mz", "rt",
        "formula", "grade", "rejection",
        "prediction", "si", "idx", "id",
        "peakarea", "peakheight" # Generic headers, but LipidSearch usually has specific sample columns
        , "f-score", "c-score", "meanarea", "meanheight"
    )

    sample_columns <- character(0)
    annotation_columns <- character(0)

    for (col in headers) {
        col_lower <- tolower(col)
        # Check if it matches any metadata pattern
        is_metadata <- any(sapply(known_metadata, function(m) col_lower == m))

        # Also check for regex patterns (e.g. "Norm" columns if we only want raw)
        # Note: Keeping Norm columns as potential samples for now, user can exclude in design

        if (!is_metadata && is.numeric(data[[col]])) {
            sample_columns <- c(sample_columns, col)
        } else {
            annotation_columns <- c(annotation_columns, col)
        }
    }

    logger::log_info(sprintf(
        "Detected %d sample columns and %d annotation columns",
        length(sample_columns),
        length(annotation_columns)
    ))

    return(list(
        data = as.data.frame(data),
        headers = headers,
        detected_columns = detected_columns,
        sample_columns = sample_columns,
        annotation_columns = annotation_columns,
        format = "lipidsearch",
        is_pattern = defaults$is_pattern
    ))
}

