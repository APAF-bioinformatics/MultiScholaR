# Keep the composite image-assembly helper cluster top-level so later waves can
# move this seam without reopening the module wrapper body.
buildMetabNormLabelPlot <- function(title) {
    ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0, y = 0.5, label = title, size = 5, hjust = 0) +
        ggplot2::xlim(0, 1) +
        ggplot2::theme_void() +
        ggplot2::theme(
            plot.margin = ggplot2::margin(5, 5, 5, 5),
            panel.background = ggplot2::element_blank()
        )
}

buildMetabNormTitlePlot <- function(title) {
    ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = title, size = 6, fontface = "bold", hjust = 0.5) +
        ggplot2::xlim(0, 1) +
        ggplot2::theme_void() +
        ggplot2::theme(
            plot.margin = ggplot2::margin(5, 5, 10, 5),
            panel.background = ggplot2::element_blank()
        )
}

loadMetabNormImageAsPlot <- function(
    file_path,
    fileExistsFn = file.exists,
    readPngFn = png::readPNG,
    rasterGrobFn = grid::rasterGrob,
    logWarnFn = logger::log_warn
) {
    if (is.na(file_path) || !fileExistsFn(file_path)) {
        return(ggplot2::ggplot() + ggplot2::theme_void())
    }

    tryCatch({
        img <- readPngFn(file_path)
        grob <- rasterGrobFn(img, interpolate = TRUE)
        ggplot2::ggplot() +
            ggplot2::annotation_custom(grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
            ggplot2::theme_void()
    }, error = function(e) {
        logWarnFn(sprintf("[generateCompositeFromFiles] Could not load image: %s", file_path))
        ggplot2::ggplot() + ggplot2::theme_void()
    })
}

generateMetabNormCompositeFromFiles <- function(
    plot_files,
    ncol = 3,
    row_labels = NULL,
    column_labels = NULL,
    logInfoFn = logger::log_info,
    logErrorFn = logger::log_error,
    warningFn = warning,
    requireNamespaceFn = requireNamespace,
    buildLabelPlotFn = buildMetabNormLabelPlot,
    buildTitlePlotFn = buildMetabNormTitlePlot,
    loadImageAsPlotFn = loadMetabNormImageAsPlot,
    wrapPlotsFn = patchwork::wrap_plots,
    plotLayoutFn = patchwork::plot_layout,
    combineLayoutFn = function(plot, layout) plot + layout,
    fileExistsFn = file.exists,
    gcFn = gc
) {
    logInfoFn(sprintf("[generateCompositeFromFiles] Generating composite from %d files...", length(plot_files)))

    if (!requireNamespaceFn("patchwork", quietly = TRUE)) {
        warningFn("patchwork package required for composite generation")
        return(NULL)
    }
    if (!requireNamespaceFn("ggplot2", quietly = TRUE)) {
        warningFn("ggplot2 package required for composite generation")
        return(NULL)
    }
    if (!requireNamespaceFn("png", quietly = TRUE)) {
        warningFn("png package required for composite generation")
        return(NULL)
    }

    tryCatch({
        n_files <- length(plot_files)
        n_plot_types <- n_files / ncol

        if (is.null(row_labels)) {
            all_labels <- letters[1:n_files]
            row_labels <- split(paste0(all_labels, ")"), rep(1:n_plot_types, each = ncol))
            names(row_labels) <- paste0("row", seq_len(n_plot_types))
        }

        plot_sections <- list()
        height_values <- c()

        if (!is.null(column_labels) && length(column_labels) == ncol) {
            title_plots <- lapply(column_labels, buildTitlePlotFn)
            plot_sections <- append(plot_sections, list(
                wrapPlotsFn(title_plots, ncol = ncol)
            ))
            height_values <- c(height_values, 0.2)
            logInfoFn("[generateCompositeFromFiles] Added column titles")
        }

        row_names <- names(row_labels)

        for (i in seq_along(row_names)) {
            row_name <- row_names[i]
            labels <- row_labels[[row_name]]

            start_idx <- (i - 1) * ncol + 1
            end_idx <- min(i * ncol, n_files)
            row_files <- plot_files[start_idx:end_idx]

            has_files <- any(!is.na(row_files) & vapply(
                row_files,
                function(path) !is.na(path) && fileExistsFn(path),
                logical(1)
            ))

            if (has_files) {
                label_plots <- lapply(labels, buildLabelPlotFn)
                image_plots <- lapply(row_files, loadImageAsPlotFn)

                plot_sections <- append(plot_sections, list(
                    wrapPlotsFn(label_plots, ncol = ncol),
                    wrapPlotsFn(image_plots, ncol = ncol)
                ))
                height_values <- c(height_values, 0.1, 1)

                logInfoFn(sprintf("[generateCompositeFromFiles] Added row: %s", row_name))
            } else {
                logInfoFn(sprintf("[generateCompositeFromFiles] Skipping empty row: %s", row_name))
            }
        }

        if (length(plot_sections) == 0) {
            warningFn("No valid plot sections to combine")
            return(NULL)
        }

        logInfoFn("[generateCompositeFromFiles] Combining plot sections...")
        combined_plot <- combineLayoutFn(
            wrapPlotsFn(plot_sections, ncol = 1),
            plotLayoutFn(heights = height_values)
        )

        plot_width <- 4 + (ncol * 3)
        plot_height <- 4 + (length(height_values) * 2)

        rm(plot_sections)
        gcFn()

        list(plot = combined_plot, width = plot_width, height = plot_height)
    }, error = function(e) {
        logErrorFn(paste("[generateCompositeFromFiles] Error:", e$message))
        NULL
    })
}

