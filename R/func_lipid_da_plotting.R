# ----------------------------------------------------------------------------
# generateLipidDAVolcanoPlotGlimma
# ----------------------------------------------------------------------------
#' Generate interactive Glimma volcano plot for lipidomics DA results
#'
#' @description Creates an interactive volcano plot using the Glimma package
#'   for lipidomics differential abundance results. Supports per-assay
#'   or combined viewing.
#'
#' @param da_results_list Results list from `runLipidsDA()`.
#' @param selected_contrast Contrast to display (from friendly_name or comparison).
#' @param selected_assay Optional assay to filter (NULL for combined view).
#' @param da_q_val_thresh Q-value threshold for significance marking.
#' @param lipid_id_column Column name for lipid IDs.
#' @param annotation_column Column name for lipid annotations (for labels).
#'
#' @return A Glimma HTML widget or NULL if generation fails.
#'
#' @importFrom Glimma glimmaVolcano
#' @importFrom dplyr filter mutate case_when select
#' @importFrom rlang sym
#' @importFrom stringr str_extract
#' @importFrom logger log_info log_error log_warn
#' @export
generateLipidDAVolcanoPlotGlimma <- function(
  da_results_list,
  selected_contrast = NULL,
  selected_assay = NULL,
  da_q_val_thresh = 0.05,
  lipid_id_column = "lipid_id",
  annotation_column = "lipid_name"
) {
    # [D66:START] -------------------------
    d66_log <- function(...) message(sprintf("[D66] %s", paste0(...)))
    d66_log("=== ENTER generateLipidDAVolcanoPlotGlimma ===")
    d66_log("  selected_contrast = ", if (is.null(selected_contrast)) "NULL" else selected_contrast)
    d66_log("  selected_assay = ", if (is.null(selected_assay)) "NULL" else selected_assay)
    # [D66:END] ---------------------------

    logger::log_info("--- Entering generateLipidDAVolcanoPlotGlimma ---")
    logger::log_info(sprintf("   selected_contrast = %s", selected_contrast))
    logger::log_info(sprintf("   selected_assay = %s", ifelse(is.null(selected_assay), "NULL (combined)", selected_assay)))

    if (is.null(da_results_list) || is.null(da_results_list$da_lipids_long)) {
        # [D66:START]
        d66_log("  ERROR: da_results_list or da_lipids_long is NULL")
        d66_log("    da_results_list is NULL = ", is.null(da_results_list))
        if (!is.null(da_results_list)) {
            d66_log("    da_results_list names = ", paste(names(da_results_list), collapse = ", "))
            d66_log("    da_lipids_long is NULL = ", is.null(da_results_list$da_lipids_long))
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
    da_lipids_long <- da_results_list$da_lipids_long
    contrasts_results <- da_results_list$contrasts_results

    # [D66:START] -------------------------
    d66_log("  da_lipids_long dims = ", nrow(da_lipids_long), " x ", ncol(da_lipids_long))
    d66_log("  da_lipids_long columns = ", paste(colnames(da_lipids_long), collapse = ", "))
    if (nrow(da_lipids_long) > 0 && "comparison" %in% colnames(da_lipids_long)) {
        d66_log("  unique comparisons = ", paste(unique(da_lipids_long$comparison), collapse = ", "))
    }
    if (nrow(da_lipids_long) > 0 && "friendly_name" %in% colnames(da_lipids_long)) {
        d66_log("  unique friendly_names = ", paste(unique(da_lipids_long$friendly_name), collapse = ", "))
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
    contrast_data <- da_lipids_long |>
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
                !is.na(lipid_name) & lipid_name != "",
                lipid_name,
                lipid_id
            )
        ) |>
        # FINAL Inf/NaN scrub for plotting columns
        dplyr::filter(!is.na(logFC_safe), !is.nan(logFC_safe), !is.infinite(logFC_safe)) |>
        dplyr::filter(!is.na(lqm), !is.nan(lqm), !is.infinite(lqm)) |>
        dplyr::select(
            lipid_id,
            lipid_name,
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
    assay_data <- theObject@lipid_data[[fit_assay]]
    sample_cols <- intersect(colnames(assay_data), theObject@design_matrix[[theObject@sample_id]])

    counts_mat <- as.matrix(assay_data[, sample_cols, drop = FALSE])
    rownames(counts_mat) <- assay_data[[theObject@lipid_id_column]]

    # Protect Glimma Javascript from crashing due to NA/NaN/Inf in counts matrix
    if (any(is.na(counts_mat) | is.nan(counts_mat) | is.infinite(counts_mat))) {
        logger::log_info("   Sanitizing counts_mat for Glimma by replacing NA/NaN/Inf with 0")
        counts_mat[is.na(counts_mat) | is.nan(counts_mat) | is.infinite(counts_mat)] <- 0
    }

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
    volcano_tab <- volcano_tab[match(fit_feature_ids, volcano_tab$lipid_id), , drop = FALSE]

    # Handle any potential mismatches (though match IDs should exist)
    if (any(is.na(volcano_tab$lipid_id))) {
        logger::log_warn("   Some features in fit object not found in volcano_tab. Removing NAs.")
        volcano_tab <- volcano_tab[!is.na(volcano_tab$lipid_id), ]
        # Also subset other inputs if necessary, but usually fit_obj is the reference
    }

    # [FIX]: Use glimmaXY for better stability and explicit control over data linking
    # 1. Sync counts_mat to volcano_tab exactly
    if (!is.null(counts_mat)) {
        counts_mat <- counts_mat[as.character(volcano_tab$lipid_id), , drop = FALSE]
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
    if ("lipid_id" %in% names(clean_anno)) {
        clean_anno <- clean_anno |> dplyr::relocate(lipid_id)
    }
    rownames(clean_anno) <- as.character(volcano_tab$lipid_id)

    # 3. Prepare status vector (1 = Up, -1 = Down, 0 = Not Sig)
    status_vec <- ifelse(volcano_tab$significant == "Up", 1,
                    ifelse(volcano_tab$significant == "Down", -1, 0))

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

            logger::log_info("--- Exiting generateLipidDAVolcanoPlotGlimma (success) ---")
            return(glimma_widget)
        },
        error = function(e) {
            logger::log_error(sprintf("   Glimma widget generation failed: %s", e$message))
            return(NULL)
        }
    )
}

# ----------------------------------------------------------------------------
# generateLipidDAHeatmap
# ----------------------------------------------------------------------------
#' Generate heatmap for lipidomics DA results
#'
#' @description Creates a heatmap of top differentially abundant lipids
#'   with customizable clustering and scaling options.
#'
#' @param da_results_list Results list from `runLipidsDA()`.
#' @param selected_contrast Contrast to display.
#' @param selected_assay Optional assay filter (NULL for combined).
#' @param top_n Number of top lipids to include (by |logFC|).
#' @param clustering_method Hierarchical clustering method.
#' @param distance_method Distance metric for clustering.
#' @param cluster_rows Logical, cluster rows.
#' @param cluster_cols Logical, cluster columns.
#' @param scale_data Scaling option: "row", "column", "both", or "none".
#' @param color_scheme Color palette name.
#' @param show_lipid_names Logical, show lipid name labels.
#' @param da_q_val_thresh Q-value threshold for significance.
#'
#' @param tree_cut_method Method for cutting the row tree: "none", "k_clusters", "height_cutoff", "dynamic".
#' @param n_clusters Number of clusters (for k_clusters).
#' @param cut_height Height cutoff (for height_cutoff).
#' @param min_cluster_size Minimum cluster size (for dynamic).
#'
#' @return A list containing the Heatmap object and cluster assignments, or NULL.
#'
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom circlize colorRamp2
#' @importFrom stats hclust dist cor cutree
#' @importFrom dplyr filter arrange desc slice_head
#' @importFrom logger log_info log_error log_warn
#' @export
generateLipidDAHeatmap <- function(
  da_results_list,
  selected_contrast = NULL,
  selected_assay = NULL,
  top_n = 50,
  clustering_method = "ward.D2",
  distance_method = "euclidean",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale_data = "row",
  color_scheme = "RdBu",
  show_lipid_names = FALSE,
  da_q_val_thresh = 0.05,
  tree_cut_method = "none",
  n_clusters = 4,
  cut_height = 10,
  min_cluster_size = 5
) {
    logger::log_info("--- Entering generateLipidDAHeatmap ---")

    # Handle potential NA/NULL parameters from Shiny
    cluster_rows <- if (is.null(cluster_rows) || is.na(cluster_rows)) TRUE else cluster_rows
    cluster_cols <- if (is.null(cluster_cols) || is.na(cluster_cols)) TRUE else cluster_cols
    tree_cut_method <- if (is.null(tree_cut_method) || is.na(tree_cut_method)) "none" else tree_cut_method
    top_n <- if (is.null(top_n) || is.na(top_n)) 50 else top_n
    logger::log_info(sprintf("   selected_contrast = %s, top_n = %d", selected_contrast, top_n))

    if (is.null(da_results_list) || is.null(da_results_list$da_lipids_long)) {
        logger::log_warn("   No DE results available")
        return(NULL)
    }

    if (is.null(selected_contrast)) {
        logger::log_warn("   No contrast selected")
        return(NULL)
    }

    # Get data
    da_lipids_long <- da_results_list$da_lipids_long

    # Extract comparison name
    comparison_to_search <- stringr::str_extract(selected_contrast, "^[^=]+")
    if (is.na(comparison_to_search)) {
        comparison_to_search <- selected_contrast
    }

    # Filter for significant results in selected contrast
    contrast_data <- da_lipids_long |>
        dplyr::filter(comparison == comparison_to_search | friendly_name == comparison_to_search) |>
        dplyr::filter(fdr_qvalue < da_q_val_thresh)

    # Optionally filter by assay
    if (!is.null(selected_assay) && !is.na(selected_assay) && selected_assay != "Combined") {
        contrast_data <- contrast_data |>
            dplyr::filter(assay == selected_assay)
    }

    if (nrow(contrast_data) == 0) {
        logger::log_warn("   No significant lipids found")
        return(NULL)
    }

    # Get top N by absolute logFC
    top_lipids <- contrast_data |>
        dplyr::arrange(dplyr::desc(abs(logFC))) |>
        dplyr::slice_head(n = top_n)

    logger::log_info(sprintf("   Selected %d top lipids", nrow(top_lipids)))

    # Get expression matrix
    theObject <- da_results_list$theObject

    # Determine which assay(s) to use
    if (!is.null(selected_assay) && !is.na(selected_assay) && selected_assay != "Combined") {
        assays_to_use <- selected_assay
    } else {
        assays_to_use <- unique(top_lipids$assay)
    }

    # Build expression matrix from relevant assays
    expr_list <- lapply(assays_to_use, function(assay_name) {
        assay_data <- theObject@lipid_data[[assay_name]]
        if (is.null(assay_data)) {
            return(NULL)
        }

        # Get lipid IDs for this assay
        assay_lipid_ids <- top_lipids |>
            dplyr::filter(assay == assay_name) |>
            dplyr::pull(lipid_id)

        if (length(assay_lipid_ids) == 0) {
            return(NULL)
        }

        # Filter to selected lipids
        rows_to_keep <- assay_data[[theObject@lipid_id_column]] %in% assay_lipid_ids
        assay_subset <- assay_data[rows_to_keep, , drop = FALSE]

        # Get sample columns
        sample_cols <- intersect(colnames(assay_subset), theObject@design_matrix[[theObject@sample_id]])

        # Build matrix
        mat <- as.matrix(assay_subset[, sample_cols, drop = FALSE])
        rownames(mat) <- assay_subset[[theObject@lipid_id_column]]

        return(mat)
    })

    # Combine matrices (if multiple assays)
    expr_matrix <- do.call(rbind, expr_list[!sapply(expr_list, is.null)])

    if (is.null(expr_matrix) || nrow(expr_matrix) == 0) {
        logger::log_warn("   Could not build expression matrix")
        return(NULL)
    }

    logger::log_info(sprintf("   Expression matrix: %d x %d", nrow(expr_matrix), ncol(expr_matrix)))

    # Apply scaling
    if (scale_data == "row") {
        expr_matrix <- t(scale(t(expr_matrix)))
    } else if (scale_data == "column") {
        expr_matrix <- scale(expr_matrix)
    } else if (scale_data == "both") {
        expr_matrix <- t(scale(t(expr_matrix)))
        expr_matrix <- scale(expr_matrix)
    }

    # Handle NA/Inf from scaling
    expr_matrix[is.na(expr_matrix)] <- 0
    expr_matrix[is.infinite(expr_matrix)] <- 0

    # Build row labels (lipid names if requested)
    if (show_lipid_names) {
        # Map IDs to names
        id_to_name <- stats::setNames(top_lipids$lipid_name, top_lipids$lipid_id)
        row_labels <- id_to_name[rownames(expr_matrix)]
        row_labels[is.na(row_labels)] <- rownames(expr_matrix)[is.na(row_labels)]
    } else {
        row_labels <- rownames(expr_matrix)
    }

    # Calculate clustering
    row_clust <- NULL
    col_clust <- NULL

    if (cluster_rows && nrow(expr_matrix) > 1) {
        if (distance_method %in% c("pearson", "spearman")) {
            row_dist <- stats::as.dist(1 - stats::cor(t(expr_matrix), method = distance_method, use = "pairwise.complete.obs"))
        } else {
            row_dist <- stats::dist(expr_matrix, method = distance_method)
        }
        row_clust <- stats::hclust(row_dist, method = clustering_method)
    }

    if (cluster_cols && ncol(expr_matrix) > 1) {
        if (distance_method %in% c("pearson", "spearman")) {
            col_dist <- stats::as.dist(1 - stats::cor(expr_matrix, method = distance_method, use = "pairwise.complete.obs"))
        } else {
            col_dist <- stats::dist(t(expr_matrix), method = distance_method)
        }
        col_clust <- stats::hclust(col_dist, method = clustering_method)
    }

    # Tree Cutting Logic
    row_clusters <- NULL
    if (isTRUE(cluster_rows) && !is.null(row_clust) && !is.na(tree_cut_method) && tree_cut_method != "none") {
        logger::log_info(sprintf("   Applying tree cutting method: %s", tree_cut_method))

        if (tree_cut_method == "k_clusters") {
            if (!is.na(n_clusters) && n_clusters > 1 && n_clusters <= nrow(expr_matrix)) {
                row_clusters <- stats::cutree(row_clust, k = n_clusters)
            } else {
                logger::log_warn("   Invalid n_clusters for k_clusters method")
            }
        } else if (tree_cut_method == "height_cutoff") {
            if (!is.na(cut_height)) {
                row_clusters <- stats::cutree(row_clust, h = cut_height)
            }
        } else if (tree_cut_method == "dynamic") {
            if (requireNamespace("dynamicTreeCut", quietly = TRUE)) {
                # dynamicTreeCut requires a dissimilarity matrix, not just the tree
                # We need to recalculate dist if it wasn't preserved, but we have row_dist
                if (exists("row_dist")) {
                    row_clusters <- dynamicTreeCut::cutreeDynamic(
                        dendro = row_clust,
                        distM = as.matrix(row_dist),
                        deepSplit = 2,
                        pamRespectsDendro = FALSE,
                        minClusterSize = min_cluster_size,
                        verbose = 0
                    )
                    names(row_clusters) <- row_clust$labels
                }
            } else {
                logger::log_warn("   dynamicTreeCut package not installed, skipping dynamic cutting")
            }
        }
    }

    # Color palette
    color_fn <- switch(color_scheme,
        "RdBu" = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        "RdYlBu" = circlize::colorRamp2(c(-2, 0, 2), c("blue", "yellow", "red")),
        "coolwarm" = circlize::colorRamp2(c(-2, 0, 2), c("#3B4CC0", "white", "#B40426")),
        "viridis" = circlize::colorRamp2(seq(-2, 2, length.out = 256), viridisLite::viridis(256)),
        "plasma" = circlize::colorRamp2(seq(-2, 2, length.out = 256), viridisLite::plasma(256)),
        "inferno" = circlize::colorRamp2(seq(-2, 2, length.out = 256), viridisLite::inferno(256)),
        circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
    )

    # Get column annotations (groups)
    dm <- theObject@design_matrix
    rownames(dm) <- dm[[theObject@sample_id]]
    col_groups <- dm[colnames(expr_matrix), theObject@group_id]

    # Dynamic group color mapping
    unique_groups <- sort(unique(col_groups))
    group_colors <- if (length(unique_groups) <= 9) {
        RColorBrewer::brewer.pal(max(3, length(unique_groups)), "Set1")[1:length(unique_groups)]
    } else {
        grDevices::rainbow(length(unique_groups))
    }
    names(group_colors) <- unique_groups

    # Add cluster annotation if available
    left_annotation <- NULL
    if (!is.null(row_clusters)) {
        # Create a safe color palette for clusters
        n_clust <- length(unique(row_clusters))
        if (n_clust <= 8) {
            clust_colors <- RColorBrewer::brewer.pal(max(3, n_clust), "Dark2")[1:n_clust]
        } else {
            clust_colors <- grDevices::rainbow(n_clust)
        }
        names(clust_colors) <- unique(row_clusters)

        left_annotation <- ComplexHeatmap::HeatmapAnnotation(
            Cluster = as.character(row_clusters),
            col = list(Cluster = clust_colors),
            which = "row"
        )
    }

    # Create heatmap
    tryCatch(
        {
            hm <- ComplexHeatmap::Heatmap(
                expr_matrix,
                name = "Z-score",
                col = color_fn,
                cluster_rows = if (is.null(row_clust)) cluster_rows else row_clust,
                cluster_columns = if (is.null(col_clust)) cluster_cols else col_clust,
                show_row_names = show_lipid_names,
                row_labels = row_labels,
                show_column_names = TRUE,
                column_title = paste("Top", nrow(expr_matrix), "DE Lipids:", selected_contrast),
                row_title = "Lipids",
                top_annotation = ComplexHeatmap::HeatmapAnnotation(
                    Group = col_groups,
                    col = list(Group = group_colors)
                ),
                left_annotation = left_annotation
            )

            logger::log_info("--- Exiting generateLipidDAHeatmap (success) ---")
            return(list(
                plot = hm,
                row_clusters = row_clusters,
                col_clusters = col_clust
            ))
        },
        error = function(e) {
            logger::log_error(sprintf("   Heatmap generation failed: %s", e$message))
            return(NULL)
        }
    )
}

# ----------------------------------------------------------------------------
# generateLipidDAVolcanoStatic
# ----------------------------------------------------------------------------
#' Generate static ggplot2 volcano plot for lipidomics DA results
#'
#' @description Creates a static volcano plot using ggplot2 as a fallback
#'   when Glimma is not available or for export purposes.
#'
#' @param da_results_list Results list from `runLipidsDA()`.
#' @param selected_contrast Contrast to display.
#' @param selected_assay Optional assay filter (NULL for combined/faceted).
#' @param da_q_val_thresh Q-value threshold for significance.
#' @param lfc_threshold Log fold-change threshold for significance lines.
#' @param show_labels Logical, label top significant lipids.
#' @param n_labels Number of top lipids to label.
#'
#' @return A ggplot2 object.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline labs theme_minimal scale_color_manual facet_wrap
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr filter arrange slice_head
#' @export
generateLipidDAVolcanoStatic <- function(
  da_results_list,
  selected_contrast = NULL,
  selected_assay = NULL,
  da_q_val_thresh = 0.05,
  lfc_threshold = 1,
  show_labels = TRUE,
  n_labels = 10
) {
    if (is.null(da_results_list) || is.null(da_results_list$da_lipids_long)) {
        return(NULL)
    }

    if (is.null(selected_contrast)) {
        return(NULL)
    }

    # Get data
    da_lipids_long <- da_results_list$da_lipids_long

    # Extract comparison name
    comparison_to_search <- stringr::str_extract(selected_contrast, "^[^=]+")
    if (is.na(comparison_to_search)) {
        comparison_to_search <- selected_contrast
    }

    # Filter
    plot_data <- da_lipids_long |>
        dplyr::filter(comparison == comparison_to_search | friendly_name == comparison_to_search)

    if (!is.null(selected_assay) && selected_assay != "Combined") {
        plot_data <- plot_data |>
            dplyr::filter(assay == selected_assay)
    }

    if (nrow(plot_data) == 0) {
        return(NULL)
    }

    # Add plot columns
    plot_data <- plot_data |>
        dplyr::mutate(
            neg_log10_q = -log10(fdr_qvalue),
            display_name = ifelse(!is.na(lipid_name) & lipid_name != "",
                lipid_name, lipid_id
            )
        )

    # Get top lipids for labeling
    if (show_labels) {
        top_to_label <- plot_data |>
            dplyr::filter(significant != "NS") |>
            dplyr::arrange(fdr_qvalue) |>
            dplyr::slice_head(n = n_labels)
    }

    # Build plot
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(
        x = logFC,
        y = neg_log10_q,
        color = significant
    )) +
        ggplot2::geom_point(alpha = 0.6, size = 2) +
        ggplot2::geom_hline(
            yintercept = -log10(da_q_val_thresh),
            linetype = "dashed",
            color = "gray50"
        ) +
        ggplot2::geom_vline(
            xintercept = c(-lfc_threshold, lfc_threshold),
            linetype = "dashed",
            color = "gray50"
        ) +
        ggplot2::scale_color_manual(
            values = c("Up" = "#E64B35", "Down" = "#4DBBD5", "NS" = "gray70"),
            name = "Significance"
        ) +
        ggplot2::labs(
            title = paste("Volcano Plot:", selected_contrast),
            x = "Log2 Fold Change",
            y = "-log10(Q-value)"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "bottom")

    # Add labels
    if (show_labels && nrow(top_to_label) > 0) {
        p <- p + ggrepel::geom_text_repel(
            data = top_to_label,
            ggplot2::aes(label = display_name),
            size = 3,
            max.overlaps = 20
        )
    }

    # Facet if multiple assays
    if (is.null(selected_assay) || selected_assay == "Combined") {
        if (length(unique(plot_data$assay)) > 1) {
            p <- p + ggplot2::facet_wrap(~assay, scales = "free")
        }
    }

    return(p)
}

