#' @export
setMethod(
  f = "plotRle",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, grouping_variable, yaxis_limit = c(), sample_label = NULL) {
    protein_quant_table <- theObject@protein_quant_table
    protein_id_column <- theObject@protein_id_column
    design_matrix <- theObject@design_matrix
    sample_id <- theObject@sample_id

    frozen_protein_matrix <- protein_quant_table |>
      column_to_rownames(protein_id_column) |>
      as.matrix()

    design_matrix <- as.data.frame(design_matrix)

    if (!is.null(sample_label)) {
      if (sample_label %in% colnames(design_matrix)) {
        rownames(design_matrix) <- design_matrix[, sample_label]
        colnames(frozen_protein_matrix) <- design_matrix[, sample_label]
      }
    } else {
      rownames(design_matrix) <- design_matrix[, sample_id]
    }

    # print( design_matrix)

    rowinfo_vector <- NA
    if (!is.na(grouping_variable)) {
      rowinfo_vector <- design_matrix[colnames(frozen_protein_matrix), grouping_variable]
    }

    print(rownames(design_matrix))
    print(colnames(frozen_protein_matrix))
    print(rowinfo_vector)
    # Handle missing/non-finite values
    working_matrix <- frozen_protein_matrix
    working_matrix[!is.finite(working_matrix)] <- NA

    rle_plot_before_cyclic_loess <- plotRleHelper(t(working_matrix),
      rowinfo = rowinfo_vector,
      yaxis_limit = yaxis_limit
    )

    return(rle_plot_before_cyclic_loess)
  }
)

#' @export
setMethod(
  f = "plotRleList",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, list_of_columns, yaxis_limit = c()) {
    protein_quant_table <- theObject@protein_quant_table
    protein_id_column <- theObject@protein_id_column
    design_matrix <- theObject@design_matrix
    sample_id <- theObject@sample_id

    frozen_protein_matrix <- protein_quant_table |>
      column_to_rownames(protein_id_column) |>
      as.matrix()

    design_matrix <- as.data.frame(design_matrix)
    rownames(design_matrix) <- design_matrix[, sample_id]

    # print( design_matrix)

    runOneRle <- function(column_name) {
      rowinfo_vector <- NA

      if (column_name %in% colnames(design_matrix)) {
        rowinfo_vector <- design_matrix[colnames(frozen_protein_matrix), column_name]
      }

      rle_plot_before_cyclic_loess <- plotRleHelper(t(frozen_protein_matrix),
        rowinfo = rowinfo_vector,
        yaxis_limit = yaxis_limit
      )

      return(rle_plot_before_cyclic_loess)
    }

    list_of_rle_plots <- purrr::map(list_of_columns, runOneRle)

    names(list_of_rle_plots) <- list_of_columns

    return(list_of_rle_plots)
  }
)

#' @export
savePlotRleList <- function(input_list, prefix = "RLE", suffix = c("png", "pdf"), output_dir) {
  list_of_filenames <- expand_grid(column = names(input_list), suffix = suffix) |>
    mutate(filename = paste0("RLE", "_", column, ".", suffix)) |>
    left_join(
      tibble(
        column = names(input_list),
        plots = input_list
      ),
      by = join_by(column)
    )


  purrr::walk2(
    list_of_filenames$plots,
    list_of_filenames$filename,
    \(.x, .y){
      ggsave(plot = .x, filename = file.path(output_dir, .y))
    }
  )

  list_of_filenames
}

#' @export
setMethod(
  f = "plotPca",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, grouping_variable, shape_variable = NULL, label_column, title, font_size = 8, cv_percentile = 0.90) {
    # Defensive checks
    if (!is.character(grouping_variable) || length(grouping_variable) != 1) {
      stop("grouping_variable must be a single character string")
    }

    if (!is.null(shape_variable) && (!is.character(shape_variable) || length(shape_variable) != 1)) {
      stop("shape_variable must be NULL or a single character string")
    }

    if (!grouping_variable %in% colnames(theObject@design_matrix)) {
      stop(sprintf("grouping_variable '%s' not found in design matrix", grouping_variable))
    }

    if (!is.null(shape_variable) && !shape_variable %in% colnames(theObject@design_matrix)) {
      stop(sprintf("shape_variable '%s' not found in design matrix", shape_variable))
    }

    protein_quant_table <- theObject@protein_quant_table
    protein_id_column <- theObject@protein_id_column
    design_matrix <- theObject@design_matrix
    sample_id <- theObject@sample_id

    frozen_protein_matrix <- protein_quant_table |>
      column_to_rownames(protein_id_column) |>
      as.matrix()

    frozen_protein_matrix_pca <- frozen_protein_matrix
    frozen_protein_matrix_pca[!is.finite(frozen_protein_matrix_pca)] <- NA

    if (is.null(label_column) ||
        length(label_column) == 0L ||
        is.na(label_column) ||
        label_column == "") {
      label_column <- ""
    }

    required_cols <- c(sample_id, grouping_variable)
    if (!is.null(shape_variable)) {
      required_cols <- c(required_cols, shape_variable)
    }
    missing_cols <- setdiff(required_cols, colnames(design_matrix))
    if (length(missing_cols) > 0) {
      stop(sprintf("Missing columns in design matrix: %s", paste(missing_cols, collapse = ", ")))
    }

    tryCatch(
      {
        pca_plot <- plotPcaHelper(frozen_protein_matrix_pca,
          design_matrix,
          sample_id_column = sample_id,
          grouping_variable = grouping_variable,
          shape_variable = shape_variable,
          label_column = label_column,
          title = title,
          geom.text.size = font_size
        )
        return(pca_plot)
      },
      error = function(e) {
        stop(sprintf("Error in plotPcaHelper: %s", e$message))
      }
    )
  }
)

#' @export
setMethod(
  f = "plotPcaList",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, grouping_variables_list, label_column, title, font_size = 8, cv_percentile = 0.90) {
    protein_quant_table <- theObject@protein_quant_table
    protein_id_column <- theObject@protein_id_column
    design_matrix <- theObject@design_matrix
    sample_id <- theObject@sample_id

    frozen_protein_matrix <- protein_quant_table |>
      column_to_rownames(protein_id_column) |>
      as.matrix()

    frozen_protein_matrix_pca <- frozen_protein_matrix
    frozen_protein_matrix_pca[!is.finite(frozen_protein_matrix_pca)] <- NA

    if (is.null(label_column) ||
        length(label_column) == 0L ||
        is.na(label_column) ||
        label_column == "") {
      label_column <- ""
    }

    pca_plots_list <- plotPcaListHelper(frozen_protein_matrix_pca,
      design_matrix,
      sample_id_column = sample_id,
      grouping_variables_list = grouping_variables_list,
      label_column = label_column,
      title = title,
      geom.text.size = font_size
    )

    return(pca_plots_list)
  }
)

#' @export
setMethod(
  f = "plotDensity",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, grouping_variable, title = "", font_size = 8) {
    pca_plot <- plotPca(
      theObject,
      grouping_variable = grouping_variable,
      label_column = NULL,
      title = title,
      font_size = font_size
    )
    plotDensity(
      pca_plot,
      grouping_variable = grouping_variable,
      title = title,
      font_size = font_size
    )
  }
)

#' @export
savePlotPcaList <- function(input_list, prefix = "PCA", suffix = c("png", "pdf"), output_dir) {
  list_of_filenames <- expand_grid(column = names(input_list), suffix = suffix) |>
    mutate(filename = paste0("RLE", "_", column, ".", suffix)) |>
    left_join(
      tibble(
        column = names(input_list),
        plots = input_list
      ),
      by = join_by(column)
    )


  purrr::walk2(
    list_of_filenames$plots,
    list_of_filenames$filename,
    \(.x, .y){
      ggsave(plot = .x, filename = file.path(output_dir, .y))
    }
  )

  list_of_filenames
}

#' @export
setMethod(
  f = "getPcaMatrix",
  signature = "ProteinQuantitativeData",
  definition = function(theObject) {
    protein_quant_table <- theObject@protein_quant_table
    protein_id_column <- theObject@protein_id_column
    design_matrix <- theObject@design_matrix
    sample_id <- theObject@sample_id


    frozen_protein_matrix <- protein_quant_table |>
      column_to_rownames(protein_id_column) |>
      as.matrix()

    frozen_protein_matrix_pca <- frozen_protein_matrix
    frozen_protein_matrix_pca[!is.finite(frozen_protein_matrix_pca)] <- NA


    pca_mixomics_before_cyclic_loess <- mixOmics::pca(t(as.matrix(frozen_protein_matrix_pca)))$variates$X |>
      as.data.frame() |>
      rownames_to_column(sample_id) |>
      left_join(design_matrix, by = sample_id)


    return(pca_mixomics_before_cyclic_loess)
  }
)

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

#' @export
setMethod(
  f = "plotDensityList",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, grouping_variables_list, title = "", font_size = 8) {
    # Create a list of density plots for each grouping variable
    density_plots_list <- purrr::map(grouping_variables_list, function(group_var) {
      tryCatch(
        {
          plotDensity(theObject,
            grouping_variable = group_var,
            title = title,
            font_size = font_size
          )
        },
        error = function(e) {
          warning(sprintf("Error creating density plot for %s: %s", group_var, e$message))
          return(NULL)
        }
      )
    })

    # Name the list elements with the grouping variables
    names(density_plots_list) <- grouping_variables_list

    # Remove any NULL elements (failed plots)
    density_plots_list <- density_plots_list[!sapply(density_plots_list, is.null)]

    return(density_plots_list)
  }
)

#' @export
savePlotDensityList <- function(input_list, prefix = "Density", suffix = c("png", "pdf"), output_dir) {
  list_of_filenames <- expand_grid(column = names(input_list), suffix = suffix) |>
    mutate(filename = paste0(prefix, "_", column, ".", suffix)) |>
    left_join(
      tibble(
        column = names(input_list),
        plots = input_list
      ),
      by = join_by(column)
    )

  purrr::walk2(
    list_of_filenames$plots,
    list_of_filenames$filename,
    \(.x, .y) {
      ggsave(plot = .x, filename = file.path(output_dir, .y))
    }
  )

  list_of_filenames
}

