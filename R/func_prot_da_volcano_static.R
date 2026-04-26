# ----------------------------------------------------------------------------
# generateProtDAVolcanoStatic
# ----------------------------------------------------------------------------
#' Generate Static Volcano Plot for Proteomics DA Results
#'
#' @description Creates a static volcano plot using ggplot2 and ggrepel.
#'
#' @param da_results_list A list containing DA results (da_proteins_long, etc.)
#' @param selected_contrast Character string of the contrast to plot.
#' @param da_q_val_thresh Numeric threshold for q-value significance (default: 0.05).
#' @param lfc_threshold Numeric threshold for log2 fold-change (default: 1).
#' @param show_labels Logical, whether to label top significant proteins.
#' @param n_labels Number of top proteins to label (default: 10).
#'
#' @return A ggplot2 object.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline labs theme_minimal scale_color_manual
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr filter arrange slice_head
#' @export
generateProtDAVolcanoStatic <- function(
  da_results_list,
  selected_contrast = NULL,
  da_q_val_thresh = 0.05,
  lfc_threshold = 1,
  show_labels = TRUE,
  n_labels = 10
) {
  if (is.null(da_results_list) || is.null(da_results_list$da_proteins_long)) {
    return(NULL)
  }

  if (is.null(selected_contrast)) {
    return(NULL)
  }

  # Extract comparison name
  comparison_to_search <- stringr::str_extract(selected_contrast, "^[^=]+")
  if (is.na(comparison_to_search)) {
    comparison_to_search <- selected_contrast
  }

  # Filter data for selected contrast
  plot_data <- da_results_list$da_proteins_long |>
    dplyr::filter(comparison == comparison_to_search)

  if (nrow(plot_data) == 0) {
    # Try fuzzy match if exact match fails
    available_comparisons <- unique(da_results_list$da_proteins_long$comparison)
    matching_key <- which(stringr::str_detect(available_comparisons, stringr::fixed(comparison_to_search)))
    if (length(matching_key) > 0) {
      comparison_to_search <- available_comparisons[matching_key[1]]
      plot_data <- da_results_list$da_proteins_long |>
        dplyr::filter(comparison == comparison_to_search)
    }
  }

  if (nrow(plot_data) == 0) {
    return(NULL)
  }

  # Add plotting columns
  # CRITICAL FIX: Look up the ID column dynamically instead of hardcoding uniprot_acc
  id_col <- intersect(c("Protein.Ids", "Protein.ID", "Entry", "uniprot_acc", "sites_id"), names(plot_data))
  if (length(id_col) > 0) {
    id_col <- id_col[1]
    message(paste("   generateProtDAVolcanoStatic: Using ID column =", id_col))
  } else {
    # Fallback to first column if none of the above match
    id_col <- names(plot_data)[1]
    message(paste("   generateProtDAVolcanoStatic: Fallback to first column as ID =", id_col))
  }

  plot_data <- plot_data |>
    dplyr::mutate(
      neg_log10_q = -log10(fdr_qvalue),
      display_name = if ("gene_name" %in% names(plot_data)) {
        ifelse(!is.na(gene_name) & gene_name != "", gene_name, !!sym(id_col))
      } else {
        !!sym(id_col)
      },
      significant_label = case_when(
        fdr_qvalue < as.double(da_q_val_thresh) & log2FC >= lfc_threshold ~ "Up",
        fdr_qvalue < as.double(da_q_val_thresh) & log2FC <= -lfc_threshold ~ "Down",
        TRUE ~ "NS"
      )
    )

  # Define colors
  volcano_colors <- c("Up" = "#d9534f", "Down" = "#5bc0de", "NS" = "grey70")

  # Create base plot
  p <- ggplot2::ggplot(plot_data, aes(x = log2FC, y = neg_log10_q, color = significant_label)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.5) +
    ggplot2::geom_hline(yintercept = -log10(da_q_val_thresh), linetype = "dashed", color = "grey50") +
    ggplot2::geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed", color = "grey50") +
    ggplot2::scale_color_manual(values = volcano_colors) +
    ggplot2::labs(
      title = paste("Volcano Plot:", comparison_to_search),
      x = "Log2 Fold Change",
      y = "-Log10 Q-value",
      color = "Status"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  # Add labels for top proteins
  if (show_labels) {
    top_to_label <- plot_data |>
      dplyr::filter(significant_label != "NS") |>
      dplyr::arrange(fdr_qvalue) |>
      dplyr::slice_head(n = n_labels)

    if (nrow(top_to_label) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = top_to_label,
        aes(label = display_name),
        size = 3,
        box.padding = 0.5,
        max.overlaps = 15,
        show.legend = FALSE
      )
    }
  }

  return(p)
}
