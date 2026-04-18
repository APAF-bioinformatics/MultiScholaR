#' Extract Best K Values Per Assay from RUV Results
#'
#' @title Extract best k values per assay from RUV results
#' @param ruv_results List of per-assay RUV optimization results
#' @return Named list of k values (or NA for failed assays)
#' @noRd
extractBestKPerAssay <- function(ruv_results) {
    purrr::map(ruv_results, \(x) {
        if (isTRUE(x$success)) x$best_k else NA_integer_
    })
}

#' Extract Control Feature Indices Per Assay from RUV Results
#'
#' @title Extract control feature indices per assay from RUV results
#' @param ruv_results List of per-assay RUV optimization results
#' @return Named list of control indices (or NULL for failed assays)
#' @noRd
extractCtrlPerAssay <- function(ruv_results) {
    purrr::map(ruv_results, \(x) {
        if (isTRUE(x$success)) x$control_genes_index else NULL
    })
}

#' Build Combined RUV Optimization Results Table
#'
#' @title Build combined RUV optimization results table
#' @param ruv_results List of per-assay RUV optimization results
#' @return Data frame with Assay column and optimization metrics
#' @importFrom dplyr bind_rows mutate
#' @noRd
buildCombinedRuvTable <- function(ruv_results) {
    rows <- purrr::imap(ruv_results, \(result, assay_name) {
        if (isTRUE(result$success)) {
            data.frame(
                Assay = assay_name
                , Best_K = result$best_k
                , Best_Percentage = result$best_percentage
                , Separation_Score = round(result$separation_score, 4)
                , Num_Controls = sum(result$control_genes_index, na.rm = TRUE)
                , Status = "Success"
                , stringsAsFactors = FALSE
            )
        } else {
            data.frame(
                Assay = assay_name
                , Best_K = NA_integer_
                , Best_Percentage = NA_real_
                , Separation_Score = NA_real_
                , Num_Controls = NA_integer_
                , Status = paste("Failed:", result$error)
                , stringsAsFactors = FALSE
            )
        }
    })

    dplyr::bind_rows(rows)
}

#' Build Normalization Configuration Object
#'
#' @title Build normalization configuration object
#' @param input Shiny input object
#' @return List of configuration parameters for state saving
#' @noRd
buildNormConfig <- function(input) {
    list(
        itsd = list(
            applied = input$apply_itsd %||% FALSE
            , method = input$itsd_method %||% "median"
        )
        , log2 = list(
            offset = input$log_offset %||% 1
        )
        , normalization = list(
            method = input$norm_method %||% "cyclicloess"
        )
        , ruv = list(
            mode = input$ruv_mode %||% "skip"
            , grouping_variable = input$ruv_grouping_variable
            , auto_percentage_min = input$auto_percentage_min
            , auto_percentage_max = input$auto_percentage_max
            , separation_metric = input$separation_metric
            , k_penalty_weight = input$k_penalty_weight
            , adaptive_k_penalty = input$adaptive_k_penalty
            , manual_k = input$ruv_k
            , manual_percentage = input$ruv_percentage
        )
        , timestamp = Sys.time()
    )
}

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

