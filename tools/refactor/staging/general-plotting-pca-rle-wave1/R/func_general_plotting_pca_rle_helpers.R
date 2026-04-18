# ----------------------------------------------------------------------------
# plotPcaHelper
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Plot PCA Helper
#'
#' @description
#' This function performs a Principal Component Analysis (PCA) on the given data and generates a ggplot object for visualization.
#' It includes filtering steps to select the most variable features before running PCA.
#'
#' @param data A matrix or data frame of quantitative data, with features in rows and samples in columns.
#' @param design_matrix A data frame containing sample metadata.
#' @param sample_id_column The name of the column in `design_matrix` that identifies samples. Defaults to "Sample_ID".
#' @param grouping_variable The column in `design_matrix` used for coloring points in the PCA plot. Defaults to "group".
#' @param shape_variable An optional column in `design_matrix` to be used for the shape aesthetic of points.
#' @param label_column An optional column in `design_matrix` to be used for labeling points.
#' @param title The title of the plot.
#' @param geom.text.size The size of the text labels if `label_column` is provided. Defaults to 11.
#' @param ncomp The number of principal components to compute. Defaults to 2.
#' @param cv_percentile The percentile threshold for selecting features based on coefficient of variation. Defaults to 0.90 (top 10% most variable features).
#' @param ... Additional arguments passed to other methods (not currently used).
#'
#' @return A `ggplot` object representing the PCA plot.
#'
#' @importFrom ggplot2 ggplot aes geom_point xlab ylab labs theme element_blank scale_x_continuous scale_y_continuous coord_cartesian
#' @importFrom ggrepel geom_text_repel
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join
#' @importFrom rlang sym
#'
#' @examples
#' # Create dummy data
#' data <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' colnames(data) <- paste0("Sample", 1:10)
#' rownames(data) <- paste0("Protein", 1:100)
#'
#' # Create dummy design matrix
#' design_matrix <- data.frame(
#'   Sample_ID = paste0("Sample", 1:10),
#'   group = rep(c("A", "B"), each = 5),
#'   batch = rep(c("X", "Y"), times = 5)
#' )
#'
#' # Generate PCA plot
#' plotPcaHelper(
#'   data = data,
#'   design_matrix = design_matrix,
#'   title = "PCA of Dummy Data",
#'   grouping_variable = "group",
#'   shape_variable = "batch"
#' )
#'
#' @export
plotPcaHelper <- function(data,
                          design_matrix,
                          sample_id_column = "Sample_ID",
                          grouping_variable = "group",
                          shape_variable = NULL,
                          label_column = NULL,
                          title, geom.text.size = 11, ncomp = 2,
                          cv_percentile = 0.90,
                          ...) {
  # DEBUG66: Memory helper using BOTH R heap AND process memory (Task Manager view)
  checkMem <- function(step) {
    checkMemoryBoth(step, context = "plotPcaHelper")
  }

  message("+===========================================================================+")
  message("|  DEBUG66: Entering plotPcaHelper                                          |")
  message("+===========================================================================+")
  entry_mem <- checkMem("Entry")
  message(sprintf("   DEBUG66 [plotPcaHelper] Arg: nrow(data) = %d, ncol(data) = %d", nrow(data), ncol(data)))
  message(sprintf("   DEBUG66 [plotPcaHelper] Arg: nrow(design_matrix) = %d", nrow(design_matrix)))
  message(sprintf("   DEBUG66 [plotPcaHelper] Arg: grouping_variable = '%s'", grouping_variable))
  message(sprintf("   DEBUG66 [plotPcaHelper] Arg: shape_variable = '%s'", ifelse(is.null(shape_variable), "NULL", shape_variable)))
  message(sprintf("   DEBUG66 [plotPcaHelper] Arg: cv_percentile = %.2f", cv_percentile))

  # Ensure design_matrix is a data frame
  design_matrix <- as.data.frame(design_matrix)

  # --- START of modification ---
  if (nrow(data) > 1) {
    # 1. Filter by presence in samples
    # Keep features that have values in at least 5% of the samples.
    num_samples <- ncol(data)
    min_samples_present <- floor(0.05 * num_samples)
    data_abundant <- data[rowSums(!is.na(data)) >= min_samples_present, ]

    # 2. Filter by Coefficient of Variation (CV) on the already filtered data
    if (nrow(data_abundant) > 1) {
      cvs <- apply(data_abundant, 1, function(x) {
        if (all(is.na(x))) {
          return(NA)
        }
        s <- sd(x, na.rm = TRUE)
        m <- mean(x, na.rm = TRUE)
        if (is.na(s) || is.na(m) || abs(m) < 1e-6) {
          return(0)
        }
        s / m
      })

      # Find the CV threshold based on the specified percentile
      cv_threshold <- quantile(cvs, cv_percentile, na.rm = TRUE)

      if (!is.na(cv_threshold) && cv_threshold > 0) {
        data_filtered <- data_abundant[which(cvs >= cv_threshold), ]
        message(sprintf("   [plotPcaHelper] Filtering by CV. Rows before: %d, Rows after: %d", nrow(data), nrow(data_filtered)))
      } else {
        data_filtered <- data_abundant
      }
    } else {
      data_filtered <- data_abundant
    }
  } else {
    data_filtered <- data
  }
  # --- END of modification ---

  checkMem("After CV filtering")
  message(sprintf("   DEBUG66 [plotPcaHelper] data_filtered dims: %d x %d", nrow(data_filtered), ncol(data_filtered)))

  checkMem("Before mixOmics::pca")
  pca_input <- t(as.matrix(data_filtered))
  if (requireNamespace("mixOmics", quietly = TRUE)) {
    message("   DEBUG66 [plotPcaHelper] Step: Calling mixOmics::pca()...")
    pca.res <- mixOmics::pca(pca_input, ncomp = ncomp)
  } else {
    message("   DEBUG66 [plotPcaHelper] mixOmics not available; falling back to stats::prcomp()")
    pca_fallback <- stats::prcomp(pca_input, center = TRUE, scale. = TRUE)
    component_names <- paste0("PC", seq_len(min(ncomp, ncol(pca_fallback$x))))
    pca_scores <- as.data.frame(pca_fallback$x[, component_names, drop = FALSE])
    explained_variance <- (pca_fallback$sdev^2) / sum(pca_fallback$sdev^2)
    pca.res <- list(
      variates = list(X = pca_scores),
      prop_expl_var = list(
        X = stats::setNames(
          explained_variance[seq_along(component_names)],
          component_names
        )
      )
    )
  }
  checkMem("After mixOmics::pca")
  message(sprintf("   DEBUG66 [plotPcaHelper] pca.res class: %s", class(pca.res)[1]))
  proportion_explained <- pca.res$prop_expl_var

  checkMem("Before temp_tbl creation")
  message("   DEBUG66 [plotPcaHelper] Step: Creating temp_tbl with left_join...")
  temp_tbl <- pca.res$variates$X |>
    as.data.frame() |>
    rownames_to_column(var = sample_id_column) |>
    left_join(design_matrix, by = sample_id_column)
  checkMem("After temp_tbl creation")
  message(sprintf("   DEBUG66 [plotPcaHelper] temp_tbl dims: %d x %d", nrow(temp_tbl), ncol(temp_tbl)))

  # --- START of fix for shape aesthetic ---
  # Ensure the shape variable is treated as a discrete factor for plotting
  if (!is.null(shape_variable) && shape_variable %in% colnames(temp_tbl)) {
    temp_tbl[[shape_variable]] <- as.factor(temp_tbl[[shape_variable]])
  }
  # --- END of fix for shape aesthetic ---

  # More defensive check for grouping variables
  if (!grouping_variable %in% colnames(temp_tbl)) {
    stop(sprintf("Grouping variable '%s' not found in the data", grouping_variable))
  }

  if (!is.null(shape_variable) && !shape_variable %in% colnames(temp_tbl)) {
    stop(sprintf("Shape variable '%s' not found in the data", shape_variable))
  }

  checkMem("Before ggplot creation")
  message("   DEBUG66 [plotPcaHelper] Step: Building ggplot object...")

  # Create base plot with appropriate aesthetics based on whether shape_variable is NULL
  if (is.null(label_column) || label_column == "") {
    if (is.null(shape_variable)) {
      # No shape variation, only color
      base_plot <- temp_tbl |>
        ggplot(aes(PC1, PC2, color = !!sym(grouping_variable)))
    } else {
      # Both color and shape
      base_plot <- temp_tbl |>
        ggplot(aes(PC1, PC2, color = !!sym(grouping_variable), shape = !!sym(shape_variable)))
    }
  } else {
    if (!label_column %in% colnames(temp_tbl)) {
      stop(sprintf("Label column '%s' not found in the data", label_column))
    }

    if (is.null(shape_variable)) {
      # No shape variation, only color, with labels
      base_plot <- temp_tbl |>
        ggplot(aes(PC1, PC2, color = !!sym(grouping_variable), label = !!sym(label_column)))
    } else {
      # Both color and shape, with labels
      base_plot <- temp_tbl |>
        ggplot(aes(PC1, PC2,
          color = !!sym(grouping_variable), shape = !!sym(shape_variable),
          label = !!sym(label_column)
        ))
    }
  }

  # Calculate the percentage label for axis (e.g., 0.90 -> "top 10%")
  cv_percent_label <- paste0("top ", round((1 - cv_percentile) * 100, 0), "%")

  output <- base_plot +
    geom_point(size = 3) +
    xlab(paste("PC1 (", round(proportion_explained$X[["PC1"]] * 100, 0), "% of ", cv_percent_label, " CV)", sep = "")) +
    ylab(paste("PC2 (", round(proportion_explained$X[["PC2"]] * 100, 0), "% of ", cv_percent_label, " CV)", sep = "")) +
    labs(title = title) +
    theme(legend.title = element_blank()) +
    scale_x_continuous(labels = function(x) format(x, scientific = FALSE, digits = 3)) +
    scale_y_continuous(labels = function(x) format(x, scientific = FALSE, digits = 3))

  # Add explicit color and shape scales to support >6 discrete levels
  # Get color palette (supports many levels)
  categorical_colors <- getCategoricalColourPalette()

  # Add explicit color scale
  output <- output + scale_color_manual(values = categorical_colors)

  # Add explicit shape scale (15 distinct shapes - first 6 match ggplot2 defaults)
  if (!is.null(shape_variable)) {
    shape_values <- c(16, 17, 15, 3, 7, 8, 0, 1, 2, 4, 5, 6, 9, 10, 11)
    output <- output + scale_shape_manual(values = shape_values)
  }

  # Calculate axis limits based on the full range of data
  pc1_range <- range(temp_tbl$PC1, na.rm = TRUE)
  pc2_range <- range(temp_tbl$PC2, na.rm = TRUE)
  buffer_pc1 <- (pc1_range[2] - pc1_range[1]) * 0.05 # 5% buffer
  buffer_pc2 <- (pc2_range[2] - pc2_range[1]) * 0.05 # 5% buffer

  output <- output + coord_cartesian(
    xlim = c(pc1_range[1] - buffer_pc1, pc1_range[2] + buffer_pc1),
    ylim = c(pc2_range[1] - buffer_pc2, pc2_range[2] + buffer_pc2)
  )

  if (!is.null(label_column) && label_column != "") {
    output <- output + geom_text_repel(size = geom.text.size, show.legend = FALSE)
  }

  checkMem("After ggplot creation")
  message(sprintf("   DEBUG66 [plotPcaHelper] output object size: %s", format(object.size(output), units = "auto")))

  checkMem("Exit")
  message("+===========================================================================+")
  message("|  DEBUG66: Exiting plotPcaHelper                                           |")
  message("+===========================================================================+")

  # Strip captured environment to prevent memory bloat
  output$plot_env <- rlang::base_env()

  # class(output) <- "ggplot"
  output
}

# ----------------------------------------------------------------------------
# plotPcaListHelper
# ----------------------------------------------------------------------------
#' @export
plotPcaListHelper <- function(data,
                              design_matrix,
                              sample_id_column = "Sample_ID",
                              grouping_variables_list = c("group"),
                              label_column = NULL,
                              title, geom.text.size = 11, ncomp = 2,
                              cv_percentile = 0.90,
                              ...) {
  # --- START of modification ---
  if (nrow(data) > 1) {
    # 1. Filter by presence in samples
    # Keep features that have values in at least 5% of the samples.
    num_samples <- ncol(data)
    min_samples_present <- floor(0.05 * num_samples)
    data_abundant <- data[rowSums(!is.na(data)) >= min_samples_present, ]

    # 2. Filter by Coefficient of Variation (CV) on the already filtered data
    if (nrow(data_abundant) > 1) {
      cvs <- apply(data_abundant, 1, function(x) {
        if (all(is.na(x))) {
          return(NA)
        }
        s <- sd(x, na.rm = TRUE)
        m <- mean(x, na.rm = TRUE)
        if (is.na(s) || is.na(m) || abs(m) < 1e-6) {
          return(0)
        }
        s / m
      })

      # Find the CV threshold based on the specified percentile
      cv_threshold <- quantile(cvs, cv_percentile, na.rm = TRUE)

      if (!is.na(cv_threshold) && cv_threshold > 0) {
        data_filtered <- data_abundant[which(cvs >= cv_threshold), ]
        message(sprintf("   [plotPcaHelper] Filtering by CV. Rows before: %d, Rows after: %d", nrow(data), nrow(data_filtered)))
      } else {
        data_filtered <- data_abundant
      }
    } else {
      data_filtered <- data_abundant
    }
  } else {
    data_filtered <- data
  }
  # --- END of modification ---


  pca.res <- mixOmics::pca(t(as.matrix(data_filtered)), ncomp = ncomp)
  proportion_explained <- pca.res$prop_expl_var

  temp_tbl <- pca.res$variates$X |>
    as.data.frame() |>
    rownames_to_column(var = sample_id_column) |>
    left_join(design_matrix, by = sample_id_column)


  plotOneGgplotPca <- function(grouping_variable) {
    unique_groups <- temp_tbl |>
      distinct(!!sym(grouping_variable)) |>
      dplyr::pull(!!sym(grouping_variable))

    if (is.null(label_column) || label_column == "") {
      output <- temp_tbl |>
        ggplot(aes(PC1, PC2, col = !!sym(grouping_variable))) +
        geom_point() +
        xlab(paste("PC1 (", round(proportion_explained$X[["PC1"]] * 100, 0), "%)", sep = "")) +
        ylab(paste("PC2 (", round(proportion_explained$X[["PC2"]] * 100, 0), "%)", sep = "")) +
        labs(title = title) +
        theme(legend.title = element_blank())
    } else {
      output <- temp_tbl |>
        ggplot(aes(PC1, PC2, col = !!sym(grouping_variable), label = !!sym(label_column))) +
        geom_point() +
        geom_text_repel(size = geom.text.size, show.legend = FALSE) +
        xlab(paste("PC1 (", round(proportion_explained$X[["PC1"]] * 100, 0), "%)", sep = "")) +
        ylab(paste("PC2 (", round(proportion_explained$X[["PC2"]] * 100, 0), "%)", sep = "")) +
        labs(title = title) +
        theme(legend.title = element_blank())
    }

    # Strip captured environment to prevent memory bloat
    output$plot_env <- rlang::base_env()
    return(output)
  }

  output_list <- purrr::map(grouping_variables_list, plotOneGgplotPca)

  output_list
}

# ----------------------------------------------------------------------------
# plotPcaGgpairs
# ----------------------------------------------------------------------------
#' @export
plotPcaGgpairs <- function(
  data_matrix,
  design_matrix,
  grouping_variable,
  sample_id_column,
  ncomp = 2
) {
  pca.res <- mixOmics::pca(t(as.matrix(data_matrix)), ncomp = ncomp)


  pca_prop_explained_helper <- function(pca_obj, comp_idx) {
    proportion_explained <- pca.res$prop_expl_var

    pc_label <- paste0("PC", comp_idx)

    perc_label <- paste(paste0(pc_label, " ("), round(proportion_explained$X[[pc_label]] * 100, 0), "%)", sep = "")

    perc_label
  }

  pc_list <- purrr::map_chr(
    seq_len(ncomp),
    \(comp_idx){
      pca_prop_explained_helper(pca.res, comp_idx)
    }
  )

  pca_variates_x <- pca.res$variates$X

  colnames(pca_variates_x) <- pc_list

  pca_plot_ggpairs <- pca_variates_x |>
    as.data.frame() |>
    rownames_to_column(sample_id_column) |>
    left_join(design_matrix,
      by = join_by(!!sym(sample_id_column) == !!sym(sample_id_column))
    ) |>
    ggpairs(
      columns = pc_list, aes(colour = !!sym(grouping_variable), fill = !!sym(grouping_variable), alpha = 0.4),
      legend = 1
    )

  pca_plot_ggpairs
}

# ----------------------------------------------------------------------------
# plotRleHelper
# ----------------------------------------------------------------------------
#' @title Plot RLE
#' @export
#' @param Y  Rows = Samples, Columns = Proteins or Peptides
plotRleHelper <- function(Y, rowinfo = NULL, probs = c(
                            0.05, 0.25, 0.5, 0.75,
                            0.95
                          ), yaxis_limit = c(-0.5, 0.5)) {
  #  checks = check.ggplot()
  # if (checks) {
  rle <- t(apply(t(Y) - apply(Y, 2, function(x) {
    median(x, na.rm = TRUE)
  }), 2, function(x) {
    quantile(x, probs = probs, na.rm = TRUE)
  }))
  colnames(rle) <- c(
    "min", "lower", "middle", "upper",
    "max"
  )
  df <- cbind(data.frame(rle.x.factor = rownames(rle)), data.frame(rle))

  if (!is.null(rowinfo)) {
    rowinfo <- data.frame(rowinfo = rowinfo)
    df_temp <- cbind(df, rowinfo)

    my.x.factor.levels <- df_temp |>
      arrange(rowinfo) |>
      distinct(rle.x.factor) |>
      dplyr::pull(rle.x.factor)

    df <- df_temp |>
      mutate(rle.x.factor = factor(rle.x.factor,
        levels = my.x.factor.levels
      )) |>
      arrange(rowinfo)
  }

  rleplot <- ggplot(df, aes(x = .data[["rle.x.factor"]])) +
    geom_boxplot(
      aes(
        lower = .data[["lower"]],
        middle = .data[["middle"]],
        upper = .data[["upper"]],
        max = .data[["max"]],
        min = .data[["min"]]
      ),
      stat = "identity"
    ) +
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 90) # , axis.ticks.x = element_blank()
    ) +
    theme(axis.title.y = element_blank(), axis.text.y = element_text(size = rel(1.5))) +
    geom_hline(yintercept = 0)


  if (length(yaxis_limit) == 2) {
    rleplot <- rleplot +
      coord_cartesian(ylim = yaxis_limit)
  }


  if (!is.null(rowinfo)) {
    if (ncol(rowinfo) == 1) {
      rleplot <- rleplot + aes(fill = rowinfo) + labs(fill = "")
    }
  }

  # Strip captured environment to prevent memory bloat
  rleplot$plot_env <- rlang::base_env()
  return(rleplot)
  # }
  # else return(FALSE)
}

# ----------------------------------------------------------------------------
# getMaxMinBoxplot
# ----------------------------------------------------------------------------
#' @title Get Max and Min Boxplot
#' @export
#' @description Input a ggplot2 boxplot, return the maximum and minimum data point adjusted by the adjust_factor.
#' @param input_boxplot A ggplot2 boxplot object.
#' @param adjust_factor A numeric value to adjust the maximum and minimum data point.
getMaxMinBoxplot <- function(input_boxplot, adjust_factor = 0.05) {
  df_min <- min(input_boxplot$data$min, na.rm = TRUE)

  df_max <- max(input_boxplot$data$max, na.rm = TRUE)

  if (df_min > 0) {
    df_min <- df_min * (1 - adjust_factor)
  } else {
    df_min <- df_min * (1 + adjust_factor)
  }

  if (df_max > 0) {
    df_max <- df_max * (1 + adjust_factor)
  } else {
    df_max <- df_max * (1 - adjust_factor)
  }

  return(c(df_min, df_max))
}

# ----------------------------------------------------------------------------
# rlePcaPlotList
# ----------------------------------------------------------------------------
#' @export
rlePcaPlotList <- function(list_of_data_matrix, list_of_design_matrix,
                           sample_id_column = Sample_ID, grouping_variable = group, list_of_descriptions) {
  rle_list <- purrr::pmap(
    list(data_matrix = list_of_data_matrix, description = list_of_descriptions, design_matrix = list_of_design_matrix),
    function(data_matrix, description, design_matrix) {
      plotRleHelper(t(as.matrix(data_matrix)),
        rowinfo = design_matrix[colnames(data_matrix), as_name(enquo(grouping_variable))]
      ) +
        labs(title = description)
    }
  )

  pca_list <- purrr::pmap(
    list(data_matrix = list_of_data_matrix, description = list_of_descriptions, design_matrix = list_of_design_matrix),
    function(data_matrix, description, design_matrix) {
      plotPcaHelper(data_matrix,
        design_matrix = design_matrix,
        sample_id_column = sample_id_column,
        grouping_variable = grouping_variable,
        title = description, cex = 7
      )
    }
  )

  list_of_plots <- c(rle_list, pca_list)

  rle_pca_plots_arranged <- ggarrange(
    plotlist = list_of_plots, nrow = 2, ncol = length(list_of_descriptions),
    common.legend = FALSE, legend = "bottom", widths = 10, heights = 10
  )

  rle_pca_plots_arranged
}

