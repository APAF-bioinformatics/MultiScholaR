registerLipidNormStaticQcImageOutputs <- function(output, renderQcImageForAssay) {
    output_specs <- list(
        pca_post_filter_assay1 = list(1, "pca", "pre_norm")
        , pca_post_norm_assay1 = list(1, "pca", "post_norm")
        , pca_ruv_corrected_assay1 = list(1, "pca", "ruv_corrected")
        , pca_post_filter_assay2 = list(2, "pca", "pre_norm")
        , pca_post_norm_assay2 = list(2, "pca", "post_norm")
        , pca_ruv_corrected_assay2 = list(2, "pca", "ruv_corrected")
        , density_post_filter_assay1 = list(1, "density", "pre_norm")
        , density_post_norm_assay1 = list(1, "density", "post_norm")
        , density_ruv_corrected_assay1 = list(1, "density", "ruv_corrected")
        , density_post_filter_assay2 = list(2, "density", "pre_norm")
        , density_post_norm_assay2 = list(2, "density", "post_norm")
        , density_ruv_corrected_assay2 = list(2, "density", "ruv_corrected")
        , rle_post_filter_assay1 = list(1, "rle", "pre_norm")
        , rle_post_norm_assay1 = list(1, "rle", "post_norm")
        , rle_ruv_corrected_assay1 = list(1, "rle", "ruv_corrected")
        , rle_post_filter_assay2 = list(2, "rle", "pre_norm")
        , rle_post_norm_assay2 = list(2, "rle", "post_norm")
        , rle_ruv_corrected_assay2 = list(2, "rle", "ruv_corrected")
        , correlation_post_filter_assay1 = list(1, "correlation", "pre_norm")
        , correlation_post_norm_assay1 = list(1, "correlation", "post_norm")
        , correlation_ruv_corrected_assay1 = list(1, "correlation", "ruv_corrected")
        , correlation_post_filter_assay2 = list(2, "correlation", "pre_norm")
        , correlation_post_norm_assay2 = list(2, "correlation", "post_norm")
        , correlation_ruv_corrected_assay2 = list(2, "correlation", "ruv_corrected")
    )

    for (output_id in names(output_specs)) {
        spec <- output_specs[[output_id]]
        output[[output_id]] <- renderQcImageForAssay(spec[[1]], spec[[2]], spec[[3]])
    }

    invisible(output)
}

registerLipidNormAssayLabelOutputs <- function(output, renderAssayLabel) {
    output_specs <- list(
        assay1_label_pca = 1
        , assay2_label_pca = 2
        , assay1_label_density = 1
        , assay2_label_density = 2
        , assay1_label_rle = 1
        , assay2_label_rle = 2
        , assay1_label_correlation = 1
        , assay2_label_correlation = 2
    )

    for (output_id in names(output_specs)) {
        output[[output_id]] <- renderAssayLabel(output_specs[[output_id]])
    }

    invisible(output)
}

registerLipidNormLogOutput <- function(output, renderNormLog) {
    output$norm_log <- renderNormLog()

    invisible(output)
}

registerLipidNormItsdSelectionOutput <- function(output, renderItsdSelectionUi) {
    output$itsd_selection_ui <- renderItsdSelectionUi()

    invisible(output)
}

registerLipidNormRuvQcOutput <- function(output, renderRuvQcUi) {
    output$ruv_qc_ui <- renderRuvQcUi()

    invisible(output)
}

registerLipidNormCorrelationFilterSummaryOutput <- function(output, renderCorrelationFilterSummary) {
    output$correlation_filter_summary <- renderCorrelationFilterSummary()

    invisible(output)
}

registerLipidNormFinalQcPlotOutput <- function(output, renderFinalQcPlot) {
    output$final_qc_plot <- renderFinalQcPlot()

    invisible(output)
}

buildLipidNormAssayLabelRenderer <- function(
    normData,
    renderTextFn = shiny::renderText
) {
    function(assaySlot) {
        renderTextFn({
            if (!is.null(normData$assay_names) && length(normData$assay_names) >= assaySlot) {
                paste("Assay:", normData$assay_names[[assaySlot]])
            } else {
                paste0("Assay ", assaySlot, ": (detecting...)")
            }
        })
    }
}

buildLipidNormQcImageRenderer <- function(
    normData,
    experimentPaths,
    renderImageFn = shiny::renderImage,
    filePathFn = file.path,
    fileExistsFn = file.exists,
    sanitizeAssayNameFn = function(assayName) {
        gsub("[^A-Za-z0-9]", "_", tolower(assayName))
    }
) {
    function(assaySlot, plotType, stagePrefix) {
        renderImageFn({
            normData$plot_refresh_trigger
            assay_names <- normData$assay_names

            if (is.null(assay_names) || length(assay_names) < assaySlot) {
                list(src = "", alt = "Assay not detected yet")
            } else {
                assay_name <- assay_names[[assaySlot]]
                safe_name <- sanitizeAssayNameFn(assay_name)
                filename <- paste0(safe_name, "_", stagePrefix, "_", plotType, ".png")

                qc_dir <- experimentPaths$lipid_qc_dir
                if (is.null(qc_dir)) {
                    list(src = "", alt = "QC directory not configured")
                } else {
                    img_path <- filePathFn(qc_dir, filename)

                    if (fileExistsFn(img_path)) {
                        list(
                            src = img_path
                            , contentType = "image/png"
                            , width = "100%"
                            , height = "auto"
                            , alt = paste(plotType, "-", assay_name)
                        )
                    } else {
                        list(src = "", alt = paste("Plot not generated yet:", filename))
                    }
                }
            }
        }, deleteFile = FALSE)
    }
}

buildLipidNormLogRenderer <- function(
    normData,
    renderTextFn = shiny::renderText
) {
    function() {
        renderTextFn({
            if (length(normData$normalization_log) == 0) {
                return("Normalization log will appear here as you apply steps...")
            }

            paste(normData$normalization_log, collapse = "\n")
        })
    }
}

buildLipidNormItsdSelectionUiRenderer <- function(
    normData,
    ns,
    renderUiFn = shiny::renderUI,
    reqFn = shiny::req,
    mapFn = purrr::map,
    wellPanelFn = shiny::wellPanel,
    headerFn = shiny::h5,
    dataTableOutputFn = DT::dataTableOutput,
    brFn = shiny::br,
    tagListFn = shiny::tagList
) {
    function() {
        renderUiFn({
            reqFn(normData$assay_names)

            assay_uis <- mapFn(normData$assay_names, \(assay_name) {
                safe_name <- gsub("[^A-Za-z0-9]", "_", tolower(assay_name))
                wellPanelFn(
                    headerFn(paste("Assay:", assay_name))
                    , dataTableOutputFn(ns(paste0("itsd_table_", safe_name)))
                    , brFn()
                )
            })

            tagListFn(assay_uis)
        })
    }
}

buildLipidNormRuvQcUiRenderer <- function(
    normData,
    ns,
    renderUiFn = shiny::renderUI,
    reqFn = shiny::req,
    mapFn = purrr::map,
    tagListFn = shiny::tagList,
    fluidRowFn = shiny::fluidRow,
    columnFn = shiny::column,
    headerFn = shiny::h5,
    subHeaderFn = shiny::h6,
    plotOutputFn = shiny::plotOutput,
    resizableFn = shinyjqui::jqui_resizable,
    wellPanelFn = shiny::wellPanel,
    verbatimTextOutputFn = shiny::verbatimTextOutput,
    brFn = shiny::br,
    dataTableOutputFn = DT::dataTableOutput,
    hrFn = shiny::hr
) {
    function() {
        renderUiFn({
            reqFn(normData$assay_names)

            assay_uis <- mapFn(normData$assay_names, \(assay_name) {
                safe_name <- gsub("[^A-Za-z0-9]", "_", tolower(assay_name))
                tagListFn(
                    fluidRowFn(
                        columnFn(
                            12,
                            headerFn(
                                paste("Assay:", assay_name),
                                style = "border-bottom: 1px solid #ddd; padding-bottom: 5px;"
                            )
                        )
                    ),
                    fluidRowFn(
                        columnFn(
                            8,
                            resizableFn(
                                plotOutputFn(ns(paste0("cancor_plot_", safe_name)), height = "400px")
                            )
                        ),
                        columnFn(
                            4,
                            wellPanelFn(
                                subHeaderFn("Optimization Summary"),
                                verbatimTextOutputFn(ns(paste0("ruv_summary_", safe_name))),
                                brFn(),
                                subHeaderFn("Results Table"),
                                dataTableOutputFn(ns(paste0("ruv_table_", safe_name)))
                            )
                        )
                    ),
                    hrFn()
                )
            })

            tagListFn(assay_uis)
        })
    }
}

buildLipidNormCorrelationFilterSummaryRenderer <- function(
    normData,
    renderTextFn = shiny::renderText
) {
    function() {
        renderTextFn({
            if (!normData$correlation_filtering_complete) {
                return("Apply correlation filter to see results...")
            }

            corr_results <- normData$correlation_results
            filtered_obj <- normData$correlation_filtered_obj
            original_obj <- if (!is.null(normData$ruv_corrected_obj)) {
                normData$ruv_corrected_obj
            } else {
                normData$post_norm_obj
            }

            if (is.null(corr_results) || length(corr_results) == 0) {
                return("No correlation results available.")
            }

            summary_lines <- c("=== Correlation Filtering Summary ===\n")

            for (assay_name in names(corr_results)) {
                assay_corr <- corr_results[[assay_name]]
                if (!is.null(assay_corr) && nrow(assay_corr) > 0) {
                    n_pairs <- nrow(assay_corr)
                    mean_corr <- round(mean(assay_corr$pearson_correlation, na.rm = TRUE), 3)
                    min_corr <- round(min(assay_corr$pearson_correlation, na.rm = TRUE), 3)
                    max_corr <- round(max(assay_corr$pearson_correlation, na.rm = TRUE), 3)

                    summary_lines <- c(summary_lines, sprintf(
                        "\n[%s]\n  Sample pairs: %d\n  Correlation: mean=%.3f, min=%.3f, max=%.3f",
                        assay_name, n_pairs, mean_corr, min_corr, max_corr
                    ))
                }
            }

            if (!is.null(original_obj) && !is.null(filtered_obj)) {
                original_samples <- nrow(original_obj@design_matrix)
                filtered_samples <- nrow(filtered_obj@design_matrix)
                removed <- original_samples - filtered_samples

                summary_lines <- c(summary_lines, sprintf(
                    "\n\n[Sample Filtering]\n  Original: %d samples\n  After filtering: %d samples\n  Removed: %d samples",
                    original_samples, filtered_samples, removed
                ))
            }

            paste(summary_lines, collapse = "")
        })
    }
}

buildLipidNormFinalQcPlotRenderer <- function(
    normData,
    getPlotAestheticsFn,
    renderPlotFn = shiny::renderPlot,
    reqFn = shiny::req,
    plotPcaFn = plotPca,
    wrapPlotsFn = patchwork::wrap_plots,
    emptyPlotFn = function(label, size) {
        ggplot2::ggplot() +
            ggplot2::annotate("text", x = 0.5, y = 0.5, label = label, size = size) +
            ggplot2::theme_void()
    }
) {
    function() {
        renderPlotFn({
            reqFn(normData$correlation_filtering_complete || normData$ruv_complete)

            current_s4 <- if (!is.null(normData$correlation_filtered_obj)) {
                normData$correlation_filtered_obj
            } else if (!is.null(normData$ruv_corrected_obj)) {
                normData$ruv_corrected_obj
            } else {
                normData$post_norm_obj
            }

            if (is.null(current_s4)) {
                return(emptyPlotFn("No data", 6))
            }

            aesthetics <- getPlotAestheticsFn()

            tryCatch({
                pca_plots <- plotPcaFn(
                    current_s4
                    , grouping_variable = aesthetics$color_var
                    , shape_variable = aesthetics$shape_var
                    , title = "Final QC - PCA"
                )

                if (is.list(pca_plots) && length(pca_plots) > 1) {
                    wrapPlotsFn(pca_plots, ncol = 1)
                } else if (is.list(pca_plots) && length(pca_plots) == 1) {
                    pca_plots[[1]]
                } else if (inherits(pca_plots, "ggplot")) {
                    pca_plots
                } else {
                    ggplot2::ggplot() + ggplot2::theme_void()
                }
            }, error = function(e) {
                emptyPlotFn(paste("Error:", e$message), 4)
            })
        })
    }
}

buildLipidNormAddLog <- function(
    normData,
    timeFn = Sys.time,
    formatTimeFn = format,
    sprintfFn = sprintf
) {
    function(message) {
        timestamp <- formatTimeFn(timeFn(), "%H:%M:%S")
        normData$normalization_log <- c(
            normData$normalization_log,
            sprintfFn("[%s] %s", timestamp, message)
        )

        invisible(normData$normalization_log)
    }
}

buildLipidNormPlotAestheticsGetter <- function(
    input,
    defaultVariable = "group"
) {
    function() {
        list(
            color_var = if (is.null(input$color_variable) || identical(input$color_variable, "")) {
                defaultVariable
            } else {
                input$color_variable
            },
            shape_var = if (is.null(input$shape_variable) || identical(input$shape_variable, "")) {
                defaultVariable
            } else {
                input$shape_variable
            }
        )
    }
}

buildLipidNormCompositeFromFilesGenerator <- function(
    requireNamespaceFn = requireNamespace,
    warningFn = warning,
    fileExistsFn = file.exists,
    readPngFn = png::readPNG,
    rasterGrobFn = grid::rasterGrob,
    wrapPlotsFn = patchwork::wrap_plots,
    plotLayoutFn = patchwork::plot_layout,
    logInfoFn = logger::log_info,
    logWarnFn = logger::log_warn,
    logErrorFn = logger::log_error
) {
    createLabelPlot <- function(title) {
        ggplot2::ggplot() +
            ggplot2::annotate("text", x = 0, y = 0.5, label = title, size = 5, hjust = 0) +
            ggplot2::xlim(0, 1) +
            ggplot2::theme_void() +
            ggplot2::theme(
                plot.margin = ggplot2::margin(5, 5, 5, 5)
                , panel.background = ggplot2::element_blank()
            )
    }

    createTitlePlot <- function(title) {
        ggplot2::ggplot() +
            ggplot2::annotate("text", x = 0.5, y = 0.5, label = title, size = 6, fontface = "bold", hjust = 0.5) +
            ggplot2::xlim(0, 1) +
            ggplot2::theme_void() +
            ggplot2::theme(
                plot.margin = ggplot2::margin(5, 5, 10, 5)
                , panel.background = ggplot2::element_blank()
            )
    }

    loadImageAsPlot <- function(filePath) {
        if (is.na(filePath) || !fileExistsFn(filePath)) {
            return(ggplot2::ggplot() + ggplot2::theme_void())
        }

        tryCatch({
            img <- readPngFn(filePath)
            grob <- rasterGrobFn(img, interpolate = TRUE)

            ggplot2::ggplot() +
                ggplot2::annotation_custom(grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
                ggplot2::theme_void()
        }, error = function(e) {
            logWarnFn(sprintf("[generateCompositeFromFiles] Could not load image: %s", filePath))
            ggplot2::ggplot() + ggplot2::theme_void()
        })
    }

    function(plot_files, ncol = 3, row_labels = NULL, column_labels = NULL) {
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
                all_labels <- letters[seq_len(n_files)]
                row_labels <- split(paste0(all_labels, ")"), rep(seq_len(n_plot_types), each = ncol))
                names(row_labels) <- paste0("row", seq_len(n_plot_types))
            }

            plot_sections <- list()
            height_values <- c()

            if (!is.null(column_labels) && length(column_labels) == ncol) {
                title_plots <- lapply(column_labels, createTitlePlot)
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

                has_files <- any(!is.na(row_files) & vapply(row_files, \(filePath) {
                    !is.na(filePath) && fileExistsFn(filePath)
                }, logical(1)))

                if (has_files) {
                    label_plots <- lapply(labels, createLabelPlot)
                    image_plots <- lapply(row_files, loadImageAsPlot)

                    plot_sections <- append(plot_sections, list(
                        wrapPlotsFn(label_plots, ncol = ncol)
                        , wrapPlotsFn(image_plots, ncol = ncol)
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
            combined_plot <- wrapPlotsFn(plot_sections, ncol = 1) +
                plotLayoutFn(heights = height_values)

            plot_width <- 4 + (ncol * 3)
            plot_height <- 4 + (length(height_values) * 2)

            rm(plot_sections)
            gc()

            list(plot = combined_plot, width = plot_width, height = plot_height)
        }, error = function(e) {
            logErrorFn(paste("[generateCompositeFromFiles] Error:", e$message))
            NULL
        })
    }
}

