# ----------------------------------------------------------------------------
# generateMetabDAVolcanoStatic
# ----------------------------------------------------------------------------
#' Generate static ggplot2 volcano plot for metabolomics DE results
#'
#' @description Creates a static volcano plot using ggplot2 as a fallback
#'   when Glimma is not available or for export purposes.
#'
#' @param da_results_list Results list from `runMetabolitesDA()`.
#' @param selected_contrast Contrast to display.
#' @param selected_assay Optional assay filter (NULL for combined/faceted).
#' @param da_q_val_thresh Q-value threshold for significance.
#' @param lfc_threshold Log fold-change threshold for significance lines.
#' @param show_labels Logical, label top significant metabolites.
#' @param n_labels Number of top metabolites to label.
#'
#' @return A ggplot2 object.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline labs theme_minimal scale_color_manual facet_wrap
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr filter arrange slice_head
#' @export
generateMetabDAVolcanoStatic <- function(
  da_results_list,
  selected_contrast = NULL,
  selected_assay = NULL,
  da_q_val_thresh = 0.05,
  lfc_threshold = 1,
  show_labels = TRUE,
  n_labels = 10
) {
    if (is.null(da_results_list) || is.null(da_results_list$da_metabolites_long)) {
        return(NULL)
    }

    if (is.null(selected_contrast)) {
        return(NULL)
    }

    # Get data
    da_metabolites_long <- da_results_list$da_metabolites_long

    # Extract comparison name
    comparison_to_search <- stringr::str_extract(selected_contrast, "^[^=]+")
    if (is.na(comparison_to_search)) {
        comparison_to_search <- selected_contrast
    }

    # Filter
    plot_data <- da_metabolites_long |>
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
            display_name = ifelse(!is.na(metabolite_name) & metabolite_name != "",
                metabolite_name, metabolite_id
            )
        )

    # Get top metabolites for labeling
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

