## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Create empty QC Grid
#' @export
setClass("GridPlotData",
  slots = list(
    pca_plots = "list",
    density_plots = "list",
    rle_plots = "list",
    pearson_plots = "list",
    cancor_plots = "list",
    limpa_plots = "list",
    pca_titles = "list",
    density_titles = "list",
    rle_titles = "list",
    pearson_titles = "list",
    cancor_titles = "list",
    limpa_titles = "list"
  )
)

#' @export
setMethod(
  "InitialiseGrid",
  signature(dummy = "ANY"),
  function(dummy = NULL) {
    new("GridPlotData",
      pca_plots = list(),
      density_plots = list(),
      rle_plots = list(),
      pearson_plots = list(),
      cancor_plots = list(),
      limpa_plots = list(),
      pca_titles = list(),
      density_titles = list(),
      rle_titles = list(),
      pearson_titles = list(),
      cancor_titles = list(),
      limpa_titles = list()
    )
  }
)

#' @export
setMethod(
  f = "createGridQC",
  signature = "GridPlotData",
  definition = function(theObject, pca_titles = NULL, density_titles = NULL, rle_titles = NULL, pearson_titles = NULL, cancor_titles = NULL, limpa_titles = NULL, ncol = 3, save_path = NULL, file_name = "pca_density_rle_pearson_corr_plots_merged", workflow_name = NULL) {
    # MEMORY OPTIMIZED: Added gc() calls and single-format save
    message("--- [createGridQC]: Entry ---")
    mem_before <- sum(gc()[, 2])
    message(sprintf("   [createGridQC] Memory Usage at entry: %.1f MB", mem_before))

    # Use stored titles if not provided as parameters
    pca_titles <- if (is.null(pca_titles)) theObject@pca_titles else pca_titles
    density_titles <- if (is.null(density_titles)) theObject@density_titles else density_titles
    rle_titles <- if (is.null(rle_titles)) theObject@rle_titles else rle_titles
    pearson_titles <- if (is.null(pearson_titles)) theObject@pearson_titles else pearson_titles
    cancor_titles <- if (is.null(cancor_titles)) theObject@cancor_titles else cancor_titles
    limpa_titles <- if (is.null(limpa_titles)) theObject@limpa_titles else limpa_titles

    createLabelPlot <- function(title) {
      # Option 1: Use xlim to expand the plot area and position text at left edge
      ggplot() +
        annotate("text", x = 0, y = 0.5, label = title, size = 5, hjust = 0) +
        xlim(0, 1) + # Explicitly set the x limits
        theme_void() +
        theme(
          plot.margin = margin(5, 5, 5, 5),
          panel.background = element_blank()
        )
    }

    # Create basic plots without titles
    createPcaPlot <- function(plot) {
      plot +
        xlim(-40, 45) + ylim(-30, 25) +
        theme(
          text = element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()
        )
    }

    createDensityPlot <- function(plot) {
      # Safety check: if plot is not a ggplot or patchwork, return it as is or wrapped
      if (!inherits(plot, "ggplot") && !inherits(plot, "patchwork")) {
        return(plot)
      }

      # For all plots, just apply the theme without adding title
      if (inherits(plot, "patchwork")) {
        plot &
          theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            text = element_text(size = 15)
          )
      } else {
        plot +
          theme(
            text = element_text(size = 15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank()
          )
      }
    }

    createRlePlot <- function(plot) {
      plot +
        theme(
          text = element_text(size = 15),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()
        )
    }

    createPearsonPlot <- function(plot) {
      plot +
        theme(text = element_text(size = 15))
    }

    createCancorPlot <- function(plot) {
      plot +
        theme(
          text = element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()
        )
    }

    # Create plots without titles - FIXED ORDER
    # Define the correct order: before_cyclic_loess, before_ruvIIIc, after_ruvIIIc
    plot_order <- c("pca_plot_before_cyclic_loess_group", "pca_plot_before_ruvIIIc_group", "pca_plot_after_ruvIIIc_group")
    density_order <- c("density_plot_before_cyclic_loess_group", "density_plot_before_ruvIIIc_group", "density_plot_after_ruvIIIc_group")
    rle_order <- c("rle_plot_before_cyclic_loess_group", "rle_plot_before_ruvIIIc_group", "rle_plot_after_ruvIIIc_group")
    pearson_order <- c("pearson_correlation_pair_before_cyclic_loess", "pearson_correlation_pair_before_ruvIIIc", "pearson_correlation_pair_after_ruvIIIc_group")
    cancor_order <- c("cancor_plot_before_cyclic_loess", "cancor_plot_before_ruvIIIc", "cancor_plot_after_ruvIIIc")

    # Extract plots in the correct order using lapply
    # MEMORY OPTIMIZATION: Process each plot type and run gc() between sections
    message("   [createGridQC] Processing PCA plots...")
    created_pca_plots <- lapply(plot_order, function(name) {
      if (!is.null(theObject@pca_plots[[name]])) createPcaPlot(theObject@pca_plots[[name]]) else NULL
    })
    created_pca_plots <- created_pca_plots[!sapply(created_pca_plots, is.null)]
    gc()

    message("   [createGridQC] Processing Density plots...")
    created_density_plots <- lapply(density_order, function(name) {
      if (!is.null(theObject@density_plots[[name]])) createDensityPlot(theObject@density_plots[[name]]) else NULL
    })
    created_density_plots <- created_density_plots[!sapply(created_density_plots, is.null)]
    gc()

    message("   [createGridQC] Processing RLE plots...")
    created_rle_plots <- lapply(rle_order, function(name) {
      if (!is.null(theObject@rle_plots[[name]])) createRlePlot(theObject@rle_plots[[name]]) else NULL
    })
    created_rle_plots <- created_rle_plots[!sapply(created_rle_plots, is.null)]
    gc()

    message("   [createGridQC] Processing Pearson plots...")
    created_pearson_plots <- lapply(pearson_order, function(name) {
      if (!is.null(theObject@pearson_plots[[name]])) createPearsonPlot(theObject@pearson_plots[[name]]) else NULL
    })
    created_pearson_plots <- created_pearson_plots[!sapply(created_pearson_plots, is.null)]
    gc()

    message("   [createGridQC] Processing Limpa plots...")
    # Limpa plots come as a list of 4 plots (dpc_curve, missing_comparison, intensity_distribution, summary)
    # We will wrap them into a single section
    created_limpa_plots <- list()
    if (length(theObject@limpa_plots) > 0) {
      # Helper: safely coerce any plot object to a ggplot-compatible object
      safe_as_ggplot <- function(p, fallback_label = "Plot unavailable") {
        tryCatch({
          if (inherits(p, "ggplot") || inherits(p, "patchwork")) {
            return(p)
          }
          # Try cowplot::ggdraw for grobs/gtables
          if (inherits(p, c("grob", "gtable", "gTree"))) {
            if (inherits(p, "gtable") && (length(p$heights) <= 1 || length(p$widths) <= 1)) {
              return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = fallback_label) + theme_void())
            }
            if (requireNamespace("cowplot", quietly = TRUE)) {
              return(cowplot::ggdraw() + cowplot::draw_grob(p))
            }
          }
          # Fallback: return a placeholder
          ggplot() + annotate("text", x = 0.5, y = 0.5, label = fallback_label) + theme_void()
        }, error = function(e) {
          message(sprintf("   [createGridQC] Warning: Could not convert limpa plot: %s", e$message))
          ggplot() + annotate("text", x = 0.5, y = 0.5, label = fallback_label) + theme_void()
        })
      }

      # Extract each component if it exists
      created_limpa_plots <- lapply(c("dpc_curve", "missing_comparison", "intensity_distribution", "summary"), function(name) {
        if (!is.null(theObject@limpa_plots[[name]])) {
           p <- theObject@limpa_plots[[name]]
           # Safely coerce to ggplot first
           p <- safe_as_ggplot(p, fallback_label = paste("Limpa", name, "unavailable"))
           # Apply theme for consistency
           createDensityPlot(p)
        } else NULL
      })
      created_limpa_plots <- created_limpa_plots[!sapply(created_limpa_plots, is.null)]
    }
    gc()

    message("   [createGridQC] Processing Cancor plots...")
    created_cancor_plots <- lapply(cancor_order, function(name) {
      if (!is.null(theObject@cancor_plots[[name]])) createCancorPlot(theObject@cancor_plots[[name]]) else NULL
    })
    # Don't filter out NULL plots to maintain column alignment
    gc()

    # Create label plots
    pca_labels <- lapply(pca_titles, createLabelPlot)
    density_labels <- lapply(density_titles, createLabelPlot)
    rle_labels <- lapply(rle_titles, createLabelPlot)
    pearson_labels <- lapply(pearson_titles, createLabelPlot)
    cancor_labels <- lapply(cancor_titles, createLabelPlot)
    limpa_labels <- lapply(limpa_titles, createLabelPlot)

    # Combine with labels above each row - modified to keep legends with their plots
    # GUARD: Only add label rows when titles exist (wrap_plots of empty list causes unit() crash)
    plot_sections <- list()
    height_values <- c()

    # Add PCA plots if they exist
    if (length(theObject@pca_plots) > 0 && length(created_pca_plots) > 0) {
      if (length(pca_labels) > 0) {
        plot_sections <- append(plot_sections, list(wrap_plots(pca_labels, ncol = ncol)))
        height_values <- c(height_values, 0.1)
      }
      plot_sections <- append(plot_sections, list(wrap_plots(created_pca_plots, ncol = ncol)))
      height_values <- c(height_values, 1)
    }

    # Add Density plots if they exist
    if (length(theObject@density_plots) > 0 && length(created_density_plots) > 0) {
      if (length(density_labels) > 0) {
        plot_sections <- append(plot_sections, list(wrap_plots(density_labels, ncol = ncol)))
        height_values <- c(height_values, 0.1)
      }
      plot_sections <- append(plot_sections, list(wrap_plots(created_density_plots, ncol = ncol)))
      height_values <- c(height_values, 1)
    }

    # Add RLE plots if they exist
    if (length(theObject@rle_plots) > 0 && length(created_rle_plots) > 0) {
      if (length(rle_labels) > 0) {
        plot_sections <- append(plot_sections, list(wrap_plots(rle_labels, ncol = ncol)))
        height_values <- c(height_values, 0.1)
      }
      plot_sections <- append(plot_sections, list(wrap_plots(created_rle_plots, ncol = ncol)))
      height_values <- c(height_values, 1)
    }

    # Add Pearson plots if they exist
    if (length(theObject@pearson_plots) > 0 && length(created_pearson_plots) > 0) {
      if (length(pearson_labels) > 0) {
        plot_sections <- append(plot_sections, list(wrap_plots(pearson_labels, ncol = ncol)))
        height_values <- c(height_values, 0.1)
      }
      plot_sections <- append(plot_sections, list(wrap_plots(created_pearson_plots, ncol = ncol)))
      height_values <- c(height_values, 1)
    }

    # Add Cancor plots if they exist (check for any non-NULL plots)
    if (length(theObject@cancor_plots) > 0 && any(!sapply(created_cancor_plots, is.null))) {
      # Replace NULL plots with empty plots to maintain column alignment
      cancor_plots_aligned <- lapply(created_cancor_plots, function(plot) {
        if (is.null(plot)) {
          ggplot() +
            theme_void() # Empty plot for NULL positions
        } else {
          plot
        }
      })

      if (length(cancor_labels) > 0) {
        plot_sections <- append(plot_sections, list(wrap_plots(cancor_labels, ncol = ncol)))
        height_values <- c(height_values, 0.1)
      }
      plot_sections <- append(plot_sections, list(wrap_plots(cancor_plots_aligned, ncol = ncol)))
      height_values <- c(height_values, 1)
    }

    # Add Limpa plots if they exist and are non-empty
    if (length(theObject@limpa_plots) > 0 && length(created_limpa_plots) > 0) {
      if (length(limpa_labels) > 0) {
        plot_sections <- append(plot_sections, list(wrap_plots(limpa_labels, ncol = ncol)))
        height_values <- c(height_values, 0.1)
      }
      plot_sections <- append(plot_sections, list(
        wrap_plots(created_limpa_plots, ncol = 2) # Limpa plots are 4, so 2x2 grid fits well
      ))
      height_values <- c(height_values, 2) # Double height for the 2x2 section
    }

    # MEMORY CLEANUP before combining plots
    message("   [createGridQC] Combining plot sections...")
    gc()

    # Create combined plot from sections
    combined_plot <- wrap_plots(plot_sections, ncol = 1) +
      plot_layout(heights = height_values)

    # Clear intermediate objects
    rm(created_pca_plots, created_density_plots, created_rle_plots, created_pearson_plots, created_cancor_plots)
    rm(pca_labels, density_labels, rle_labels, pearson_labels, cancor_labels)
    rm(plot_sections)
    gc()

    if (!is.null(save_path)) {
      # Calculate dynamic width based on number of columns
      plot_width <- 4 + (ncol * 3) # Base width + 3 units per column
      plot_height <- 4 + (length(height_values) * 2) # Base height + 2 units per row

      # MEMORY OPTIMIZATION: Save only PNG (removed PDF/SVG triple-save to reduce memory)
      message("   [createGridQC] Saving PNG only (memory optimization)...")
      tryCatch({
        ggsave(
          plot = combined_plot,
          filename = file.path(save_path, paste0(file_name, ".png")),
          width = plot_width,
          height = plot_height,
          dpi = 150 # Reduced DPI for memory efficiency
        )
        message(paste("   [createGridQC] Plot saved to:", save_path))
      }, error = function(e) {
        message(sprintf("   [createGridQC] WARNING: Failed to save plot: %s", e$message))
        message("   [createGridQC] The combined_plot object is returned but could not be rendered to PNG.")
        message("   [createGridQC] This may indicate a non-ggplot object in the plot assembly.")
      })
    }

    # Final memory cleanup
    mem_after <- sum(gc()[, 2])
    message(sprintf("   [createGridQC] Memory Usage at exit: %.1f MB (delta: %.1f MB)", mem_after, mem_after - mem_before))

    return(combined_plot)
  }
)
