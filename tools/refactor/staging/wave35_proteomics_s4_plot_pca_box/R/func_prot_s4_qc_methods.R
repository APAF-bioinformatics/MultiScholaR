#' Create PC1/PC2 Boxplots from PCA ggplot Object
#'
#' @description Extracts PCA data from a ggplot object and creates boxplots
#' for PC1 and PC2 grouped by a specified variable. Works with both classic
#' S3 ggplot objects and new S7-based ggplot2 (v3.5+) objects.
#'
#' @param theObject A ggplot object containing PCA data with PC1 and PC2 columns
#' @param grouping_variable Character string specifying the grouping variable
#' @param title Plot title (default: "")
#' @param font_size Font size for plot text (default: 8)
#'
#' @return A patchwork combined plot with PC1 and PC2 boxplots
#' @export
setMethod(
  f = "plotPcaBox",
  signature = "ANY",
  definition = function(theObject, grouping_variable, title = "", font_size = 8, show_legend = FALSE) {
    # Validate input is a ggplot-like object (works with both S3 and S7 ggplot)
    if (!inherits(theObject, c("gg", "ggplot"))) {
      stop("theObject must be a ggplot object. Got class: ", paste(class(theObject), collapse = ", "))
    }

    # Extract data directly from the ggplot object
    if (!is.null(theObject$data) && is.data.frame(theObject$data)) {
      pca_data <- as_tibble(theObject$data)
    } else {
      # Fall back to other extraction methods
      pca_data <- as_tibble(ggplot_build(theObject)$data[[1]])

      # If the data doesn't have PC1/PC2, try to extract from the plot's environment
      if (!("PC1" %in% colnames(pca_data) && "PC2" %in% colnames(pca_data))) {
        # Try to get the data from the plot's environment
        if (exists("data", envir = environment(theObject$mapping$x))) {
          pca_data <- as_tibble(get("data", envir = environment(theObject$mapping$x)))
        } else {
          stop("Could not extract PCA data from the ggplot object")
        }
      }
    }

    # Check if grouping variable exists in the data
    if (!grouping_variable %in% colnames(pca_data)) {
      stop(sprintf("grouping_variable '%s' not found in the data", grouping_variable))
    }

    # Determine legend position
    legend_pos <- if (show_legend) "right" else "none"

    # Create PC1 boxplot
    pc1_box <- ggplot(pca_data, aes(x = !!sym(grouping_variable), y = PC1, fill = !!sym(grouping_variable))) +
      geom_boxplot(notch = TRUE) +
      theme_bw() +
      labs(
        title = title,
        x = "",
        y = "PC1"
      ) +
      theme(
        legend.position = legend_pos,
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = font_size),
        plot.margin = margin(b = 0, t = 5, l = 5, r = 5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
      )

    # Add explicit fill scale to support >6 discrete levels
    categorical_colors <- getCategoricalColourPalette()
    pc1_box <- pc1_box + scale_fill_manual(values = categorical_colors)

    # Create PC2 boxplot
    pc2_box <- ggplot(pca_data, aes(x = !!sym(grouping_variable), y = PC2, fill = !!sym(grouping_variable))) +
      geom_boxplot(notch = TRUE) +
      theme_bw() +
      labs(
        x = "",
        y = "PC2"
      ) +
      theme(
        legend.position = legend_pos,
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = font_size),
        plot.margin = margin(t = 0, b = 5, l = 5, r = 5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
      )

    # Add explicit fill scale to support >6 discrete levels
    pc2_box <- pc2_box + scale_fill_manual(values = categorical_colors)

    # Combine plots with minimal spacing
    # If legend is enabled, we might want to collect guides to avoid duplication
    # but for now let's keep it simple as patchwork/cowplot handles it
    combined_plot <- pc1_box / pc2_box +
      plot_layout(heights = c(1, 1), guides = if (show_legend) "collect" else NULL) +
      plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))

    return(combined_plot)
  }
)

