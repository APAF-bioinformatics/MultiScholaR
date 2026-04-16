# ----------------------------------------------------------------------------
# generateMetabDAHeatmap
# ----------------------------------------------------------------------------
#' Generate heatmap for metabolomics DA results
#'
#' @description Creates a heatmap of top differentially expressed metabolites
#'   with customizable clustering and scaling options.
#'
#' @param da_results_list Results list from `runMetabolitesDA()`.
#' @param selected_contrast Contrast to display.
#' @param selected_assay Optional assay filter (NULL for combined).
#' @param top_n Number of top metabolites to include (by |logFC|).
#' @param clustering_method Hierarchical clustering method.
#' @param distance_method Distance metric for clustering.
#' @param cluster_rows Logical, cluster rows.
#' @param cluster_cols Logical, cluster columns.
#' @param scale_data Scaling option: "row", "column", "both", or "none".
#' @param color_scheme Color palette name.
#' @param show_metabolite_names Logical, show metabolite name labels.
#' @param da_q_val_thresh Q-value threshold for significance.
#'
#' @return A ggplot2 or ComplexHeatmap object, or NULL.
#'
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom circlize colorRamp2
#' @importFrom stats hclust dist cor
#' @importFrom dplyr filter arrange desc slice_head
#' @importFrom logger log_info log_error log_warn
#' @export
generateMetabDAHeatmap <- function(
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
  show_metabolite_names = FALSE,
  da_q_val_thresh = 0.05,
  tree_cut_method = "none",
  n_clusters = 4,
  cut_height = 10,
  min_cluster_size = 5
) {
    logger::log_info("--- Entering generateMetabDAHeatmap ---")

    # Handle potential NA/NULL parameters from Shiny
    cluster_rows <- if (is.null(cluster_rows) || is.na(cluster_rows)) TRUE else cluster_rows
    cluster_cols <- if (is.null(cluster_cols) || is.na(cluster_cols)) TRUE else cluster_cols
    tree_cut_method <- if (is.null(tree_cut_method) || is.na(tree_cut_method)) "none" else tree_cut_method
    top_n <- if (is.null(top_n) || is.na(top_n)) 50 else top_n
    logger::log_info(sprintf("   selected_contrast = %s, top_n = %d", selected_contrast, top_n))

    if (is.null(da_results_list) || is.null(da_results_list$da_metabolites_long)) {
        logger::log_warn("   No DE results available")
        return(NULL)
    }

    if (is.null(selected_contrast)) {
        logger::log_warn("   No contrast selected")
        return(NULL)
    }

    # Get data
    da_metabolites_long <- da_results_list$da_metabolites_long

    # Extract comparison name
    comparison_to_search <- stringr::str_extract(selected_contrast, "^[^=]+")
    if (is.na(comparison_to_search)) {
        comparison_to_search <- selected_contrast
    }

    # Filter for significant results in selected contrast
    contrast_data <- da_metabolites_long |>
        dplyr::filter(comparison == comparison_to_search | friendly_name == comparison_to_search) |>
        dplyr::filter(fdr_qvalue < da_q_val_thresh)

    # Optionally filter by assay
    if (!is.null(selected_assay) && !is.na(selected_assay) && selected_assay != "Combined") {
        contrast_data <- contrast_data |>
            dplyr::filter(assay == selected_assay)
    }

    if (nrow(contrast_data) == 0) {
        logger::log_warn("   No significant metabolites found")
        return(NULL)
    }

    # Get top N by absolute logFC
    top_metabolites <- contrast_data |>
        dplyr::arrange(dplyr::desc(abs(logFC))) |>
        dplyr::slice_head(n = top_n)

    logger::log_info(sprintf("   Selected %d top metabolites", nrow(top_metabolites)))

    # Get expression matrix
    theObject <- da_results_list$theObject

    # Determine which assay(s) to use
    if (!is.null(selected_assay) && !is.na(selected_assay) && selected_assay != "Combined") {
        assays_to_use <- selected_assay
    } else {
        assays_to_use <- unique(top_metabolites$assay)
    }

    # Build expression matrix from relevant assays
    expr_list <- lapply(assays_to_use, function(assay_name) {
        assay_data <- theObject@metabolite_data[[assay_name]]
        if (is.null(assay_data)) {
            return(NULL)
        }

        # Get metabolite IDs for this assay
        assay_metab_ids <- top_metabolites |>
            dplyr::filter(assay == assay_name) |>
            dplyr::pull(metabolite_id)

        if (length(assay_metab_ids) == 0) {
            return(NULL)
        }

        # Filter to selected metabolites
        rows_to_keep <- assay_data[[theObject@metabolite_id_column]] %in% assay_metab_ids
        assay_subset <- assay_data[rows_to_keep, , drop = FALSE]

        # Get sample columns
        sample_cols <- intersect(colnames(assay_subset), theObject@design_matrix[[theObject@sample_id]])

        # Build matrix
        mat <- as.matrix(assay_subset[, sample_cols, drop = FALSE])
        rownames(mat) <- assay_subset[[theObject@metabolite_id_column]]

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

    # Build row labels (metabolite names if requested)
    if (show_metabolite_names) {
        # Map IDs to names
        id_to_name <- stats::setNames(top_metabolites$metabolite_name, top_metabolites$metabolite_id)
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
                show_row_names = show_metabolite_names,
                row_labels = row_labels,
                show_column_names = TRUE,
                column_title = paste("Top", nrow(expr_matrix), "DE Metabolites:", selected_contrast),
                row_title = "Metabolites",
                top_annotation = ComplexHeatmap::HeatmapAnnotation(
                    Group = col_groups,
                    col = list(Group = group_colors)
                ),
                left_annotation = left_annotation
            )

            logger::log_info("--- Exiting generateMetabDAHeatmap (success) ---")
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

