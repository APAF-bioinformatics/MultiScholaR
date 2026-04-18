# ----------------------------------------------------------------------------
# generateMetabDAVolcanoPlotGlimma
# ----------------------------------------------------------------------------
#' Generate interactive Glimma volcano plot for metabolomics DA results
#'
#' @description Creates an interactive volcano plot using the Glimma package
#'   for metabolomics differential abundance results. Supports per-assay
#'   or combined viewing.
#'
#' @param da_results_list Results list from `runMetabolitesDA()`.
#' @param selected_contrast Contrast to display (from friendly_name or comparison).
#' @param selected_assay Optional assay to filter (NULL for combined view).
#' @param da_q_val_thresh Q-value threshold for significance marking.
#' @param metabolite_id_column Column name for metabolite IDs.
#' @param annotation_column Column name for metabolite annotations (for labels).
#'
#' @return A Glimma HTML widget or NULL if generation fails.
#'
#' @importFrom Glimma glimmaVolcano
#' @importFrom dplyr filter mutate case_when select
#' @importFrom rlang sym
#' @importFrom stringr str_extract
#' @importFrom logger log_info log_error log_warn
#' @export
generateMetabDAVolcanoPlotGlimma <- function(
  da_results_list,
  selected_contrast = NULL,
  selected_assay = NULL,
  da_q_val_thresh = 0.05,
  metabolite_id_column = "metabolite_id",
  annotation_column = "metabolite_name"
) {
    # [D66:START] -------------------------
    d66_log <- function(...) message(sprintf("[D66] %s", paste0(...)))
    d66_log("=== ENTER generateMetabDAVolcanoPlotGlimma ===")
    d66_log("  selected_contrast = ", if (is.null(selected_contrast)) "NULL" else selected_contrast)
    d66_log("  selected_assay = ", if (is.null(selected_assay)) "NULL" else selected_assay)
    # [D66:END] ---------------------------

    logger::log_info("--- Entering generateMetabDAVolcanoPlotGlimma ---")
    logger::log_info(sprintf("   selected_contrast = %s", selected_contrast))
    logger::log_info(sprintf("   selected_assay = %s", ifelse(is.null(selected_assay), "NULL (combined)", selected_assay)))

    if (is.null(da_results_list) || is.null(da_results_list$da_metabolites_long)) {
        # [D66:START]
        d66_log("  ERROR: da_results_list or da_metabolites_long is NULL")
        d66_log("    da_results_list is NULL = ", is.null(da_results_list))
        if (!is.null(da_results_list)) {
            d66_log("    da_results_list names = ", paste(names(da_results_list), collapse = ", "))
            d66_log("    da_metabolites_long is NULL = ", is.null(da_results_list$da_metabolites_long))
        }
        # [D66:END]
        logger::log_warn("   No DE results available")
        return(NULL)
    }

    if (is.null(selected_contrast)) {
        logger::log_warn("   No contrast selected")
        return(NULL)
    }

    # Get data
    da_metabolites_long <- da_results_list$da_metabolites_long
    contrasts_results <- da_results_list$contrasts_results

    # [D66:START] -------------------------
    d66_log("  da_metabolites_long dims = ", nrow(da_metabolites_long), " x ", ncol(da_metabolites_long))
    d66_log("  da_metabolites_long columns = ", paste(colnames(da_metabolites_long), collapse = ", "))
    if (nrow(da_metabolites_long) > 0 && "comparison" %in% colnames(da_metabolites_long)) {
        d66_log("  unique comparisons = ", paste(unique(da_metabolites_long$comparison), collapse = ", "))
    }
    if (nrow(da_metabolites_long) > 0 && "friendly_name" %in% colnames(da_metabolites_long)) {
        d66_log("  unique friendly_names = ", paste(unique(da_metabolites_long$friendly_name), collapse = ", "))
    }
    # [D66:END] ---------------------------

    # Extract comparison name (handle "=" format from full_format)
    comparison_to_search <- stringr::str_extract(selected_contrast, "^[^=]+")
    if (is.na(comparison_to_search)) {
        comparison_to_search <- selected_contrast
    }

    # [D66:START]
    d66_log("  comparison_to_search = ", comparison_to_search)
    # [D66:END]

    # Glimma requires matching dimensions between fit object and data.
    # For "Combined" view across multiple assays, this is not possible since
    # each assay has its own fit object with different genes.
    # Return NULL with a message - user should use static plot for combined view.
    if (is.null(selected_assay) || selected_assay == "Combined") {
        d66_log("  NOTE: Combined view not supported for Glimma - use static plot")
        logger::log_info("   Glimma not supported for Combined view (dimension mismatch). Use static plot.")
        return(NULL)
    }

    # Filter for selected contrast
    contrast_data <- da_metabolites_long |>
        dplyr::filter(comparison == comparison_to_search | friendly_name == comparison_to_search)

    # Filter by selected assay (required for Glimma)
    contrast_data <- contrast_data |>
        dplyr::filter(assay == selected_assay)

    if (nrow(contrast_data) == 0) {
        logger::log_warn(sprintf("   No data found for contrast: %s", comparison_to_search))
        return(NULL)
    }

    logger::log_info(sprintf("   Found %d rows for contrast", nrow(contrast_data)))

    # Determine which assay to use for the fit object
    if (!is.null(selected_assay) && selected_assay != "Combined" && selected_assay %in% names(contrasts_results)) {
        fit_assay <- selected_assay
    } else {
        # Use first available assay for combined view
        fit_assay <- names(contrasts_results)[1]
    }

    if (is.null(contrasts_results[[fit_assay]])) {
        logger::log_warn(sprintf("   No fit object for assay: %s", fit_assay))
        return(NULL)
    }

    fit_obj <- contrasts_results[[fit_assay]]$fit.eb

    # [D66:START] -------------------------
    d66_log("  fit_assay = ", fit_assay)
    d66_log("  fit_obj class = ", class(fit_obj)[1])
    # [D66:END] ---------------------------

    # Find coefficient index
    coef_names <- colnames(fit_obj$coefficients)

    # [D66:START] -------------------------
    d66_log("  Looking for coefficient:")
    d66_log("    comparison_to_search (friendly) = ", comparison_to_search)
    d66_log("    coef_names (actual) = ", paste(coef_names, collapse = ", "))
    # [D66:END] ---------------------------

    # Extract comparison strings from data for matching
    available_contrasts_in_data <- unique(contrast_data$comparison)
    d66_log("    Available contrasts in data: ", paste(available_contrasts_in_data, collapse = ", "))

    coef_index <- integer(0)
    # 1. Try exact match on friendly name
    coef_index <- which(coef_names == comparison_to_search)

    # 2. Try exact match on actual comparison strings found in data
    if (length(coef_index) == 0 && length(available_contrasts_in_data) > 0) {
        for (actual_comp in available_contrasts_in_data) {
            match_idx <- which(coef_names == actual_comp)
            if (length(match_idx) > 0) {
                coef_index <- match_idx
                d66_log("    Match by actual contrast from data: ", actual_comp)
                break
            }
        }
    }

    # 3. Try fuzzy match (remove spaces and non-alphanumeric)
    if (length(coef_index) == 0) {
        clean_target <- gsub("[^A-Za-z0-9]", "", comparison_to_search)
        clean_coefs <- gsub("[^A-Za-z0-9]", "", coef_names)
        coef_index <- which(clean_coefs == clean_target)
        if (length(coef_index) > 0) d66_log("    Match by fuzzy name: ", coef_names[coef_index[1]])
    }

    # 4. Try fuzzy match on actual comparison strings
    if (length(coef_index) == 0 && length(available_contrasts_in_data) > 0) {
        clean_comps <- gsub("[^A-Za-z0-9]", "", available_contrasts_in_data)
        clean_coefs <- gsub("[^A-Za-z0-9]", "", coef_names)
        for (i in seq_along(clean_comps)) {
            match_idx <- which(clean_coefs == clean_comps[i])
            if (length(match_idx) > 0) {
                coef_index <- match_idx
                d66_log("    Match by fuzzy actual contrast: ", available_contrasts_in_data[i])
                break
            }
        }
    }

    if (length(coef_index) == 0) {
        d66_log("  ERROR: No coefficient found!")
        logger::log_warn(sprintf("   No coefficient found for: %s", comparison_to_search))
        logger::log_info(sprintf("   Available coefficients: %s", paste(coef_names, collapse = ", ")))
        return(NULL)
    }

    coef_index <- coef_index[1]
    d66_log("  FINAL coef_index = ", coef_index, " (", coef_names[coef_index], ")")
    logger::log_info(sprintf("   Using coefficient index %d: %s", coef_index, coef_names[coef_index]))

    # Prepare volcano plot annotation table with Inf/NaN protection
    volcano_tab <- contrast_data |>
        dplyr::mutate(
            # Protect against zero FDR causing Inf in -log10
            fdr_safe = ifelse(fdr_qvalue == 0, min(fdr_qvalue[fdr_qvalue > 0], na.rm = TRUE) * 0.1, fdr_qvalue),
            # Handle cases where all FDR might be zero (unlikely but safe)
            fdr_safe = ifelse(is.na(fdr_safe), 1e-10, fdr_safe),
            lqm = -log10(fdr_safe),
            # Ensure logFC is numeric and non-NA for status calculation
            logFC_safe = as.numeric(ifelse(is.na(logFC) | is.nan(logFC), 0, logFC)),
            label = dplyr::case_when(
                abs(logFC_safe) >= 1 & fdr_safe >= as.double(da_q_val_thresh) ~ "Not sig., |logFC| >= 1",
                abs(logFC_safe) >= 1 & fdr_safe < as.double(da_q_val_thresh) ~ "Sig., |logFC| >= 1",
                abs(logFC_safe) < 1 & fdr_safe < as.double(da_q_val_thresh) ~ "Sig., |logFC| < 1",
                TRUE ~ "Not sig."
            ),
            display_name = ifelse(
                !is.na(metabolite_name) & metabolite_name != "",
                metabolite_name,
                metabolite_id
            )
        ) |>
        # FINAL Inf/NaN scrub for plotting columns
        dplyr::filter(!is.na(logFC_safe), !is.nan(logFC_safe), !is.infinite(logFC_safe)) |>
        dplyr::filter(!is.na(lqm), !is.nan(lqm), !is.infinite(lqm)) |>
        dplyr::select(
            metabolite_id,
            metabolite_name,
            display_name,
            assay,
            logFC = logFC_safe,
            raw_pvalue,
            fdr_qvalue,
            lqm,
            label,
            significant
        )

    # Get counts matrix for the selected assay
    theObject <- da_results_list$theObject
    assay_data <- theObject@metabolite_data[[fit_assay]]
    sample_cols <- intersect(colnames(assay_data), theObject@design_matrix[[theObject@sample_id]])

    counts_mat <- as.matrix(assay_data[, sample_cols, drop = FALSE])
    rownames(counts_mat) <- assay_data[[theObject@metabolite_id_column]]

    # Protect Glimma Javascript from crashing due to NA/NaN/Inf in counts matrix
    if (any(is.na(counts_mat) | is.nan(counts_mat) | is.infinite(counts_mat))) {
        logger::log_info("   Sanitizing counts_mat for Glimma by replacing NA/NaN/Inf with 0")
        counts_mat[is.na(counts_mat) | is.nan(counts_mat) | is.infinite(counts_mat)] <- 0
    }

    # Get groups
    dm <- theObject@design_matrix
    rownames(dm) <- dm[[theObject@sample_id]]
    groups <- dm[sample_cols, theObject@group_id]

    # [D66:START] -------------------------
    d66_log("  Preparing Glimma inputs:")
    d66_log("    counts_mat dims = ", nrow(counts_mat), " x ", ncol(counts_mat))
    d66_log("    groups = ", paste(groups, collapse = ", "))
    d66_log("    fit_obj nrow = ", nrow(fit_obj$coefficients))
    d66_log("    coef_index = ", coef_index)
    # [D66:END] ---------------------------

    # =====================================================================
    # CRITICAL: Reorder volcano_tab to match fit_obj row order
    # Glimma maps 'status' and 'anno' vectors by index. If volcano_tab is
    # sorted differently than the fit object, colors and tooltips will be
    # scrambled.
    # =====================================================================
    fit_feature_ids <- rownames(fit_obj$coefficients)
    volcano_tab <- volcano_tab[match(fit_feature_ids, volcano_tab$metabolite_id), , drop = FALSE]

    # Handle any potential mismatches (though match IDs should exist)
    if (any(is.na(volcano_tab$metabolite_id))) {
        logger::log_warn("   Some features in fit object not found in volcano_tab. Removing NAs.")
        volcano_tab <- volcano_tab[!is.na(volcano_tab$metabolite_id), ]
        # Also subset other inputs if necessary, but usually fit_obj is the reference
    }

    # [FIX]: Use glimmaXY for better stability and explicit control over data linking
    # 1. Sync counts_mat to volcano_tab exactly
    if (!is.null(counts_mat)) {
        counts_mat <- counts_mat[as.character(volcano_tab$metabolite_id), , drop = FALSE]
    }

    # 2. Sanitize the annotation dataframe to prevent D3.js / DataTables crashing
    # Convert to character, replace NAs with empty strings, and rename columns to avoid dots
    clean_anno <- volcano_tab |>
        dplyr::select(-any_of(c("logFC", "raw_pvalue", "fdr_qvalue", "lqm", "label", "significant", "log2FC", "negLog10FDR"))) |>
        dplyr::rename_with(~ gsub("[\\. ]", "_", .)) |>
        dplyr::mutate(dplyr::across(everything(), as.character)) |>
        dplyr::mutate(dplyr::across(everything(), ~ ifelse(is.na(.), "", .))) |>
        as.data.frame()
    
    # Ensure ID is first and rownames are set for linking
    if ("metabolite_id" %in% names(clean_anno)) {
        clean_anno <- clean_anno |> dplyr::relocate(metabolite_id)
    }
    rownames(clean_anno) <- as.character(volcano_tab$metabolite_id)

    # 3. Prepare status vector (1 = Up, -1 = Down, 0 = Not Sig)
    status_vec <- ifelse(volcano_tab$significant == "Up", 1,
                    ifelse(volcano_tab$significant == "Down", -1, 0))

    # [D66:START] -------------------------
    d66_log("  !! PRE-GLIMMA SANITY CHECK !!")
    d66_log("    x has NA: ", any(is.na(volcano_tab$logFC)), " / Inf: ", any(is.infinite(volcano_tab$logFC)))
    d66_log("    y has NA: ", any(is.na(volcano_tab$lqm)), " / Inf: ", any(is.infinite(volcano_tab$lqm)))
    d66_log("    status has NA: ", any(is.na(status_vec)))
    if (!is.null(counts_mat)) d66_log("    counts_mat has NA: ", any(is.na(counts_mat)), " / Inf: ", any(is.infinite(counts_mat)))
    d66_log("    groups has NA: ", any(is.na(groups)))
    d66_log("    clean_anno has NA in log2FC: ", any(is.na(clean_anno$log2FC)))
    d66_log("    clean_anno has NA in negLog10FDR: ", any(is.na(clean_anno$negLog10FDR)))
    d66_log("    checking clean_anno char cols for literal 'NA': ", any(clean_anno == "NA", na.rm=TRUE))
    # [D66:END] ---------------------------

    # Generate Glimma widget
    tryCatch(
        {
            logger::log_info("   Generating Glimma widget using glimmaXY...")

            # Use glimmaXY directly as it's more stable for these custom objects in Shiny
            glimma_widget <- Glimma::glimmaXY(
                x = volcano_tab$logFC,
                y = volcano_tab$lqm,
                xlab = "log2FC",
                ylab = "negLog10FDR",
                status = status_vec,
                counts = counts_mat,
                groups = groups,
                transform.counts = "none",
                anno = clean_anno,
                status.cols = c("#1052bd", "silver", "#cc212f"),
                main = paste("Interactive Volcano Plot:", selected_contrast)
            )

            # CRITICAL: Inject CSS fix to re-unify the split DataTables layout
            css_fix <- "<style> .dataTables_wrapper { overflow-x: auto !important; } .dataTables_scroll { display: table !important; width: auto !important; min-width: 100% !important; } .dataTables_scrollHead, .dataTables_scrollBody, .dataTables_scrollHeadInner, .dataTables_scrollHeadInner > table, .dataTables_scrollBody > table { display: contents !important; } .dataTables_scrollBody thead { display: none !important; } .dataTable th, .dataTable td { white-space: nowrap !important; } </style>"
            glimma_widget <- htmlwidgets::prependContent(glimma_widget, htmltools::HTML(css_fix))

            logger::log_info("--- Exiting generateMetabDAVolcanoPlotGlimma (success) ---")
            return(glimma_widget)
        },
        error = function(e) {
            logger::log_error(sprintf("   Glimma widget generation failed: %s", e$message))
            return(NULL)
        }
    )
}

