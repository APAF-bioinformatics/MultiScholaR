# ----------------------------------------------------------------------------
# generateProtDAHeatmap
# ----------------------------------------------------------------------------
#' Generate DA Results Heatmap with Advanced Clustering
#'
#' @description Creates a heatmap of top differentially abundant proteins
#'   with customizable clustering and scaling options. Uses ComplexHeatmap
#'   for professional-quality visualization with group annotations.
#'
#' @param da_results_list Results list from `differentialAbundanceAnalysis()`.
#' @param selected_contrast Contrast to display.
#' @param top_n_genes Number of top proteins to include (by |log2FC|).
#' @param clustering_method Hierarchical clustering method.
#' @param distance_method Distance metric for clustering.
#' @param cluster_rows Logical, cluster rows.
#' @param cluster_cols Logical, cluster columns.
#' @param scale_data Scaling option: "row", "column", "both", or "none".
#' @param color_scheme Color palette name.
#' @param tree_cut_method Method for tree cutting: "k_clusters", "height_cutoff", "dynamic", or "none".
#' @param n_clusters Number of clusters for k_clusters method.
#' @param cut_height Height for heigh_cutoff method.
#' @param min_cluster_size Minimum cluster size for dynamic method.
#'
#' @return A list containing:
#'   \item{plot}{A ComplexHeatmap object, or NULL if no data.}
#'   \item{row_clusters}{Named vector of row cluster assignments.}
#'   \item{col_clusters}{Named vector of column cluster assignments.}
#'
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @importFrom circlize colorRamp2
#' @importFrom stats hclust dist cor setNames cutree
#' @importFrom dplyr filter arrange desc slice_head pull
#' @importFrom logger log_info log_error log_warn
#' @export
generateProtDAHeatmap <- function(
  da_results_list,
  selected_contrast = NULL,
  top_n_genes = 50,
  clustering_method = "ward.D2",
  distance_method = "euclidean",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale_data = "row",
  color_scheme = "RdBu",
  show_gene_names = FALSE,
  da_q_val_thresh = 0.05,
  qvalue_column = "fdr_qvalue",
  log2fc_column = "log2FC",
  tree_cut_method = "none",
  n_clusters = 4,
  cut_height = 0.5,
  min_cluster_size = 3
) {
  logger::log_info("--- Entering generateProtDAHeatmap ---")
  logger::log_info(sprintf("   selected_contrast = %s, top_n_genes = %d", selected_contrast, top_n_genes))

  if (is.null(da_results_list) || is.null(da_results_list$da_proteins_long)) {
    logger::log_warn("   No DA results available")
    return(NULL)
  }

  if (is.null(selected_contrast)) {
    logger::log_warn("   No contrast selected")
    return(NULL)
  }

  # Get the contrast-specific results
  da_proteins_long <- da_results_list$da_proteins_long

  # Extract comparison name (handle both full_format and friendly name)
  comparison_to_search <- stringr::str_extract(selected_contrast, "^[^=]+")
  if (is.na(comparison_to_search)) {
    comparison_to_search <- selected_contrast
  }

  # Filter for significant results in selected contrast
  contrast_data <- da_proteins_long |>
    dplyr::filter(comparison == comparison_to_search | comparison == selected_contrast) |>
    dplyr::filter(!!rlang::sym(qvalue_column) < as.double(da_q_val_thresh))

  if (nrow(contrast_data) == 0) {
    logger::log_warn(sprintf("   No significant proteins found for contrast: %s", selected_contrast))
    return(NULL)
  }

  logger::log_info(sprintf("   Found %d significant proteins", nrow(contrast_data)))

  # Get top N by absolute log2FC
  top_proteins <- contrast_data |>
    dplyr::arrange(dplyr::desc(abs(!!rlang::sym(log2fc_column)))) |>
    dplyr::slice_head(n = top_n_genes)

  logger::log_info(sprintf("   Selected top %d proteins for heatmap", nrow(top_proteins)))

  # Get the S4 object
  theObject <- da_results_list$theObject
  protein_id_col <- theObject@protein_id_column

  # Extract protein IDs
  protein_ids <- top_proteins |> dplyr::pull(!!rlang::sym(protein_id_col))

  # Get expression data from the S4 object
  quant_table <- theObject@protein_quant_table

  # Filter to selected proteins
  rows_to_keep <- quant_table[[protein_id_col]] %in% protein_ids
  quant_subset <- quant_table[rows_to_keep, , drop = FALSE]

  if (nrow(quant_subset) == 0) {
    logger::log_warn("   No expression data found for selected proteins")
    return(NULL)
  }

  # Get sample columns (intersection with design matrix samples)
  sample_cols <- intersect(colnames(quant_subset), theObject@design_matrix[[theObject@sample_id]])

  # Build expression matrix
  expr_matrix <- as.matrix(quant_subset[, sample_cols, drop = FALSE])
  rownames(expr_matrix) <- quant_subset[[protein_id_col]]

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

  # Build row labels (gene names if requested)
  if (show_gene_names) {
    # Try to map protein IDs to gene names from DA results
    if ("gene_names" %in% names(top_proteins)) {
      id_to_name <- stats::setNames(top_proteins$gene_names, top_proteins[[protein_id_col]])
      row_labels <- id_to_name[rownames(expr_matrix)]
      row_labels[is.na(row_labels)] <- rownames(expr_matrix)[is.na(row_labels)]
    } else {
      row_labels <- rownames(expr_matrix)
    }
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

  # Perform tree cutting if requested
  row_clusters <- NULL
  if (!is.null(row_clust) && tree_cut_method != "none") {
    logger::log_info(sprintf("   Performing tree cutting: method=%s", tree_cut_method))

    tryCatch(
      {
        if (tree_cut_method == "k_clusters") {
          row_clusters <- stats::cutree(row_clust, k = min(n_clusters, nrow(expr_matrix)))
        } else if (tree_cut_method == "height_cutoff") {
          row_clusters <- stats::cutree(row_clust, h = cut_height)
        } else if (tree_cut_method == "dynamic") {
          if (requireNamespace("dynamicTreeCut", quietly = TRUE)) {
            row_clusters <- dynamicTreeCut::cutreeDynamic(
              dendro = row_clust,
              distM = as.matrix(row_dist),
              deepSplit = 2,
              pamRespectsDendro = FALSE,
              minClusterSize = min_cluster_size
            )
            names(row_clusters) <- row_clust$labels
          } else {
            logger::log_warn("   dynamicTreeCut package not installed, falling back to k=4")
            row_clusters <- stats::cutree(row_clust, k = min(4, nrow(expr_matrix)))
          }
        }
        logger::log_info(sprintf("   Tree cutting complete: %d clusters formed", length(unique(row_clusters))))
      },
      error = function(e) {
        logger::log_error(sprintf("   Tree cutting failed: %s", e$message))
      }
    )
  }

  # Color palette using circlize for proper scaling
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

  # Generate dynamic group colors
  unique_groups <- unique(col_groups)
  group_colors <- stats::setNames(
    grDevices::hcl.colors(length(unique_groups), "Set2"),
    unique_groups
  )

  # Add cluster annotation if available
  left_annotation <- NULL
  if (!is.null(row_clusters)) {
    cluster_colors <- stats::setNames(
      grDevices::rainbow(length(unique(row_clusters))),
      unique(row_clusters)
    )
    left_annotation <- ComplexHeatmap::HeatmapAnnotation(
      Cluster = as.character(row_clusters),
      col = list(Cluster = cluster_colors),
      which = "row"
    )
  }

  # Create heatmap with tryCatch for graceful error handling
  tryCatch(
    {
      hm <- ComplexHeatmap::Heatmap(
        expr_matrix,
        name = "Z-score",
        col = color_fn,
        cluster_rows = if (is.null(row_clust)) cluster_rows else row_clust,
        cluster_columns = if (is.null(col_clust)) cluster_cols else col_clust,
        show_row_names = show_gene_names,
        row_labels = row_labels,
        show_column_names = TRUE,
        column_title = paste("Top", nrow(expr_matrix), "DA Proteins:", selected_contrast),
        row_title = "Proteins",
        top_annotation = ComplexHeatmap::HeatmapAnnotation(
          Group = col_groups,
          col = list(Group = group_colors)
        ),
        left_annotation = left_annotation
      )

      logger::log_info("--- Exiting generateProtDAHeatmap (success) ---")

      # Return list with plot and clusters
      return(list(
        plot = hm,
        row_clusters = row_clusters,
        col_clusters = NULL
      ))
    },
    error = function(e) {
      logger::log_error(sprintf("   Heatmap generation failed: %s", e$message))
      return(NULL)
    }
  )
}

# ----------------------------------------------------------------------------
# getDataMatrix
# ----------------------------------------------------------------------------
# Helper function to get data matrix
getDataMatrix <- function(obj) {
  if (inherits(obj, "MetaboliteAssayData")) {
    message(sprintf("   Getting data matrix for object of class: %s", class(obj)[1]))
    message(sprintf("   Processing MetaboliteAssayData"))
    message(sprintf(
      "   Metabolite data dimensions: %d rows, %d cols",
      nrow(obj@metabolite_data), ncol(obj@metabolite_data)
    ))
    matrix_data <- as.matrix(obj@metabolite_data[, -1]) # Exclude Name column
    colnames(matrix_data) <- colnames(obj@metabolite_data)[-1]
    rownames(matrix_data) <- obj@metabolite_data$Name
    message(sprintf(
      "   Created matrix with dimensions: %d rows, %d cols",
      nrow(matrix_data), ncol(matrix_data)
    ))
    matrix_data
  } else if (inherits(obj, "ProteinQuantitativeData")) {
    message(sprintf("   Processing ProteinQuantitativeData"))
    message(sprintf(
      "   Protein quant table dimensions: %d rows, %d cols",
      nrow(obj@protein_quant_table), ncol(obj@protein_quant_table)
    ))
    result <- as.matrix(column_to_rownames(obj@protein_quant_table, obj@protein_id_column))
    message(sprintf(
      "   Created matrix with dimensions: %d rows, %d cols",
      nrow(result), ncol(result)
    ))
    result
  } else {
    message(sprintf("   ERROR: Unsupported object type: %s", class(obj)[1]))
    stop("Unsupported object type")
  }
}

