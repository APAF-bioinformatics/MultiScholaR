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

# Keep the auto pre-normalization QC helper top-level so a later exact-source
# wave can move it without reopening the selected-tab observer body.
generateMetabNormPreNormalizationQc <- function(
    workflowData,
    experimentPaths,
    normData,
    getPlotAestheticsFn,
    addLogFn = function(entry) invisible(entry),
    reqFn = shiny::req,
    generateMetabQcPlotsFn = generateMetabQcPlots,
    logInfoFn = logger::log_info,
    logWarnFn = logger::log_warn,
    logErrorFn = logger::log_error
) {
    logInfoFn("=== GENERATING PRE-NORMALIZATION QC PLOTS ===")

    reqFn(workflowData$state_manager)
    current_s4 <- workflowData$state_manager$getState()

    if (is.null(current_s4)) {
        logWarnFn("No S4 object available for QC plot generation")
        return()
    }

    if (inherits(current_s4, "MetaboliteAssayData")) {
        detected_assays <- names(current_s4@metabolite_data)
        if (length(detected_assays) > 0 && is.null(normData$assay_names)) {
            normData$assay_names <- detected_assays
            logInfoFn(paste("Set assay names:", paste(detected_assays, collapse = ", ")))
        }
    }

    aesthetics <- getPlotAestheticsFn()

    tryCatch({
        generateMetabQcPlotsFn(
            theObject = current_s4
            , experiment_paths = experimentPaths
            , stage = "post_filter"
            , grouping_variable = aesthetics$color_var
            , shape_variable = aesthetics$shape_var
        )

        normData$plot_refresh_trigger <- normData$plot_refresh_trigger + 1
        normData$pre_norm_qc_generated <- TRUE
        logInfoFn("Pre-normalization QC plots generated successfully")

    }, error = function(e) {
        logErrorFn(paste("Error generating pre-normalization QC:", e$message))
        addLogFn(paste("Error generating Pre-QC:", e$message))
    })
}

# Keep plot-aesthetics resolution top-level so a later extraction wave can move
# it without reopening the module server wrapper.
getPlotAesthetics <- function(
    colorVariable = NULL,
    shapeVariable = NULL
) {
    list(
        color_var = if (is.null(colorVariable) || colorVariable == "") {
            "group"
        } else {
            colorVariable
        }
        , shape_var = if (is.null(shapeVariable) || shapeVariable == "") {
            "group"
        } else {
            shapeVariable
        }
    )
}

# Keep assay-label rendering top-level so a later extraction wave can move it
# without reopening the module server wrapper's static label bindings.
renderMetabNormAssayLabel <- function(
    assaySlot,
    getAssayNamesFn,
    renderTextFn = shiny::renderText
) {
    renderTextFn({
        assayNames <- getAssayNamesFn()

        if (!is.null(assayNames) && length(assayNames) >= assaySlot) {
            paste("Assay:", assayNames[[assaySlot]])
        } else {
            paste0("Assay ", assaySlot, ": (detecting...)")
        }
    })
}

# Keep QC image rendering top-level so a later extraction wave can move it
# without reopening the module server wrapper's static QC output bindings.
renderMetabNormQcImageForAssay <- function(
    assaySlot,
    plotType,
    stagePrefix,
    normData,
    qcDir,
    renderImageFn = shiny::renderImage,
    fileExistsFn = file.exists,
    filePathFn = file.path,
    sanitizeAssayNameFn = function(assayName) {
        gsub("[^A-Za-z0-9]", "_", tolower(assayName))
    }
) {
    renderImageFn({
        normData$plot_refresh_trigger
        assayNames <- normData$assay_names

        if (is.null(assayNames) || length(assayNames) < assaySlot) {
            return(list(src = "", alt = "Assay not detected yet"))
        }

        assayName <- assayNames[[assaySlot]]
        safeName <- sanitizeAssayNameFn(assayName)
        filename <- paste0(safeName, "_", stagePrefix, "_", plotType, ".png")

        if (is.null(qcDir)) {
            return(list(src = "", alt = "QC directory not configured"))
        }

        imgPath <- filePathFn(qcDir, filename)

        if (fileExistsFn(imgPath)) {
            list(
                src = imgPath
                , contentType = "image/png"
                , width = "100%"
                , height = "auto"
                , alt = paste(plotType, "-", assayName)
            )
        } else {
            list(src = "", alt = paste("Plot not generated yet:", filename))
        }
    }, deleteFile = FALSE)
}

# Keep normalization-log mutation top-level so a later extraction wave can move
# it without reopening the module server wrapper's log state shell.
appendMetabNormNormalizationLog <- function(
    normData,
    message,
    timestampFn = function() {
        format(Sys.time(), "%H:%M:%S")
    }
) {
    timestamp <- timestampFn()
    updatedLog <- c(
        normData$normalization_log
        , sprintf("[%s] %s", timestamp, message)
    )
    normData$normalization_log <- updatedLog

    invisible(updatedLog)
}

# Keep normalization-log rendering top-level so a later extraction wave can move
# it without reopening the module server wrapper's output binding shell.
renderMetabNormNormalizationLog <- function(
    normData,
    renderTextFn = shiny::renderText,
    emptyMessage = "Normalization log will appear here as you apply steps..."
) {
    renderTextFn({
        normalizationLog <- normData$normalization_log

        if (length(normalizationLog) == 0) {
            return(emptyMessage)
        }

        paste(normalizationLog, collapse = "\n")
    })
}

# Keep correlation-filter summary rendering top-level so a later extraction
# wave can move it without reopening the module server wrapper's output shell.
renderMetabNormCorrelationFilterSummary <- function(
    normData,
    renderTextFn = shiny::renderText,
    buildSummaryFn = buildMetabNormCorrelationFilterSummary,
    incompleteMessage = "Apply correlation filter to see results..."
) {
    renderTextFn({
        if (!normData$correlation_filtering_complete) {
            return(incompleteMessage)
        }

        corrResults <- normData$correlation_results
        filteredObject <- normData$correlation_filtered_obj
        originalObject <- if (!is.null(normData$ruv_corrected_obj)) {
            normData$ruv_corrected_obj
        } else {
            normData$post_norm_obj
        }

        buildSummaryFn(
            corrResults = corrResults
            , filteredObject = filteredObject
            , originalObject = originalObject
        )
    })
}

# Keep final-QC rendering top-level so a later extraction wave can move it
# without reopening the module server wrapper's output shell.
renderMetabNormFinalQcPlot <- function(
    normData,
    colorVariableFn = function() NULL,
    shapeVariableFn = function() NULL,
    renderPlotFn = shiny::renderPlot,
    reqFn = shiny::req,
    resolveRenderStateFn = resolveMetabNormFinalQcRenderState,
    getPlotAestheticsFn = getPlotAesthetics,
    buildPcaPlotFn = buildMetabNormFinalQcPcaPlot
) {
    renderPlotFn({
        reqFn(normData$correlation_filtering_complete || normData$ruv_complete)

        finalQcState <- resolveRenderStateFn(
            correlationFilteredObject = normData$correlation_filtered_obj
            , ruvCorrectedObject = normData$ruv_corrected_obj
            , postNormObject = normData$post_norm_obj
        )

        if (isTRUE(finalQcState$isFallback)) {
            return(finalQcState$plot)
        }

        aesthetics <- getPlotAestheticsFn(
            colorVariable = colorVariableFn()
            , shapeVariable = shapeVariableFn()
        )

        buildPcaPlotFn(
            sourceObject = finalQcState$sourceObject
            , colorVar = aesthetics$color_var
            , shapeVar = aesthetics$shape_var
        )
    })
}

# Keep the reset-normalization observer wrapper top-level so a later
# extraction wave can move it without reopening the module server wrapper.
runMetabNormResetNormalizationObserverWrapper <- function(
    workflowData,
    normData,
    addLogFn = function(entry) invisible(entry),
    showNotificationFn = shiny::showNotification,
    reqFn = shiny::req,
    runObserverShellFn = runMetabNormResetNormalizationObserverShell
) {
    runObserverShellFn(
        workflowData = workflowData
        , normData = normData
        , addLogFn = addLogFn
        , showNotificationFn = showNotificationFn
        , reqFn = reqFn
    )
}

# Keep the apply-correlation observer wrapper top-level so a later extraction
# wave can move it without reopening the module server wrapper.
runMetabNormApplyCorrelationObserverWrapper <- function(
    workflowData,
    input,
    normData,
    addLogFn = function(entry) invisible(entry),
    showNotificationFn = shiny::showNotification,
    removeNotificationFn = shiny::removeNotification,
    reqFn = shiny::req,
    runObserverShellFn = runMetabNormApplyCorrelationObserverShell,
    runObserverEntryFn = runMetabNormApplyCorrelationObserverEntry
) {
    runObserverShellFn(
        workflowData = workflowData
        , normData = normData
        , threshold = input$min_pearson_correlation_threshold
        , groupingVariable = input$ruv_grouping_variable
        , addLogFn = addLogFn
        , showNotificationFn = showNotificationFn
        , removeNotificationFn = removeNotificationFn
        , reqFn = reqFn
        , runObserverEntryFn = runObserverEntryFn
    )
}

# Keep the skip-correlation observer wrapper top-level so a later extraction
# wave can move it without reopening the module server wrapper.
runMetabNormSkipCorrelationObserverWrapper <- function(
    workflowData,
    normData,
    addLogFn = function(entry) invisible(entry),
    showNotificationFn = shiny::showNotification,
    reqFn = shiny::req,
    runObserverShellFn = runMetabNormSkipCorrelationObserverShell,
    runObserverEntryFn = runMetabNormSkipCorrelationObserverEntry
) {
    runObserverShellFn(
        workflowData = workflowData
        , normData = normData
        , addLogFn = addLogFn
        , showNotificationFn = showNotificationFn
        , reqFn = reqFn
        , runObserverEntryFn = runObserverEntryFn
    )
}

# Keep the export-session observer wrapper top-level so a later extraction
# wave can move it without reopening the module server wrapper.
runMetabNormExportSessionObserverWrapper <- function(
    workflowData,
    input,
    normData,
    experimentPaths,
    experimentLabel,
    addLogFn = function(entry) invisible(entry),
    logInfoFn = logger::log_info,
    reqFn = shiny::req,
    runObserverShellFn = runMetabNormExportSessionObserverShell,
    checkReadyFn = checkMetabNormExportSessionReady,
    dispatchExportSessionFn = dispatchMetabNormExportSession
) {
    inputValues <- list(
        export_session = input$export_session
        , norm_method = input$norm_method
        , ruv_mode = input$ruv_mode
        , apply_itsd = input$apply_itsd
        , itsd_aggregation = input$itsd_aggregation
        , log_offset = input$log_offset
        , min_pearson_correlation_threshold = input$min_pearson_correlation_threshold
        , ruv_grouping_variable = input$ruv_grouping_variable
    )

    runObserverShellFn(
        workflowData = workflowData
        , normData = normData
        , inputValues = inputValues
        , experimentPaths = experimentPaths
        , experimentLabel = experimentLabel
        , addLogFn = addLogFn
        , logInfoFn = logInfoFn
        , reqFn = reqFn
        , checkReadyFn = checkReadyFn
        , dispatchExportSessionFn = dispatchExportSessionFn
    )
}

