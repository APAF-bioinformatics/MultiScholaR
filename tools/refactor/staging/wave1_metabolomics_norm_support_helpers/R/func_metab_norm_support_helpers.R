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

