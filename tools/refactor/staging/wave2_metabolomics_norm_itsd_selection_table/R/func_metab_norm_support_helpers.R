#' Build ITSD Selection Table for DT Display
#'
#' @title Build ITSD selection table for DT display
#' @description Creates a data frame suitable for DT::datatable display with
#'              pre-selection of candidate internal standards based on pattern matching.
#'
#' @param assay_data Data frame for a single assay
#' @param metabolite_id_col Column name for metabolite IDs
#' @param annotation_cols Character vector of columns to search for IS patterns
#' @param is_patterns Named list of regex patterns for IS detection
#'
#' @return Data frame with columns: feature_id, annotation, mean_intensity, cv_percent, is_candidate
#'
#' @importFrom dplyr select mutate across all_of any_of filter arrange desc
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_detect str_to_lower
#' @importFrom stats sd
#' @noRd
buildItsdSelectionTable <- function(
    assay_data
    , metabolite_id_col
    , annotation_cols = NULL
    , is_patterns = list(
        itsd = "(?i)itsd|istd|internal.?standard"
        , rt = "(?i)^rt_|_rt$|retention.?time"
        , isd = "(?i)^is_|_is$|^isd_|_isd$"
        , labeled = "(?i)_d[0-9]+$|_13c|_15n|deuterated|labeled"
    )
) {
    # --- Input validation ---
    stopifnot(
        is.data.frame(assay_data)
        , is.character(metabolite_id_col)
        , metabolite_id_col %in% colnames(assay_data)
    )

    # --- Identify sample columns (numeric) vs metadata columns ---
    numeric_cols <- names(assay_data)[sapply(assay_data, is.numeric)]
    sample_cols <- setdiff(numeric_cols, metabolite_id_col)

    if (length(sample_cols) == 0) {
        warning("No sample columns found in assay_data")
        return(data.frame(
            feature_id = character(0)
            , annotation = character(0)
            , mean_intensity = numeric(0)
            , cv_percent = numeric(0)
            , is_candidate = logical(0)
        ))
    }

    # --- Determine annotation columns to search ---
    if (is.null(annotation_cols)) {
        # Default: search metabolite ID column and common annotation columns
        potential_anno_cols <- c(
            metabolite_id_col
            , "metabolite"
            , "metabolite_identification"
            , "annotation"
            , "compound_name"
            , "name"
        )
        annotation_cols <- intersect(potential_anno_cols, colnames(assay_data))
    }

    if (length(annotation_cols) == 0) {
        annotation_cols <- metabolite_id_col
    }

    # --- Calculate summary statistics per feature ---
    feature_stats <- assay_data |>
        dplyr::mutate(
            feature_id = as.character(.data[[metabolite_id_col]])
            , mean_intensity = rowMeans(
                dplyr::across(dplyr::all_of(sample_cols))
                , na.rm = TRUE
            )
            , sd_intensity = apply(
                dplyr::across(dplyr::all_of(sample_cols))
                , 1
                , sd
                , na.rm = TRUE
            )
        ) |>
        dplyr::mutate(
            cv_percent = ifelse(
                mean_intensity > 0
                , (sd_intensity / mean_intensity) * 100
                , NA_real_
            )
        )

    # --- Create annotation string for pattern matching ---
    if (length(annotation_cols) > 1) {
        feature_stats$annotation <- apply(
            feature_stats[, annotation_cols, drop = FALSE]
            , 1
            , \(x) paste(na.omit(x), collapse = " | ")
        )
    } else {
        feature_stats$annotation <- as.character(feature_stats[[annotation_cols[1]]])
    }

    # --- Detect IS candidates using patterns ---
    combined_pattern <- paste(unlist(is_patterns), collapse = "|")

    feature_stats$is_candidate <- stringr::str_detect(
        stringr::str_to_lower(feature_stats$annotation)
        , combined_pattern
    )

    # --- Select and arrange output columns ---
    result <- feature_stats |>
        dplyr::select(
            feature_id
            , annotation
            , mean_intensity
            , cv_percent
            , is_candidate
        ) |>
        dplyr::arrange(dplyr::desc(is_candidate), cv_percent)

    return(result)
}

