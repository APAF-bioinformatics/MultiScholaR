# Keep design-driven input-choice updates top-level so a later extraction wave
# can move them without reopening the module server wrapper.
updateMetabNormDesignDrivenChoices <- function(
    session,
    designMatrix,
    updateSelectInputFn = shiny::updateSelectInput
) {
    if (is.null(designMatrix)) {
        return(invisible(NULL))
    }

    designCols <- colnames(designMatrix)

    plotAvailableVars <- intersect(
        designCols,
        c("group", "factor1", "factor2", "batch", "technical_replicate_id", "sample_id")
    )

    if (length(plotAvailableVars) > 0) {
        defaultPlotVar <- if ("group" %in% plotAvailableVars) {
            "group"
        } else {
            plotAvailableVars[[1]]
        }

        updateSelectInputFn(
            session, "color_variable"
            , choices = plotAvailableVars
            , selected = defaultPlotVar
        )
        updateSelectInputFn(
            session, "shape_variable"
            , choices = plotAvailableVars
            , selected = defaultPlotVar
        )
    }

    ruvAvailableVars <- intersect(designCols, c("group", "factor1", "factor2", "batch"))

    if (length(ruvAvailableVars) > 0) {
        defaultRuvVar <- if ("group" %in% ruvAvailableVars) {
            "group"
        } else {
            ruvAvailableVars[[1]]
        }

        updateSelectInputFn(
            session, "ruv_grouping_variable"
            , choices = ruvAvailableVars
            , selected = defaultRuvVar
        )
    }

    invisible(NULL)
}

# Keep assay-name detection top-level so a later extraction wave can move the
# observer shell without reopening the module server wrapper.
initializeMetabNormAssayNames <- function(
    stateManager,
    normData,
    reqFn = shiny::req,
    getStateFn = function(manager) manager$getState(),
    logInfoFn = logger::log_info,
    logWarnFn = logger::log_warn
) {
    reqFn(stateManager)

    tryCatch({
        currentS4 <- getStateFn(stateManager)

        if (!inherits(currentS4, "MetaboliteAssayData")) {
            return(invisible(NULL))
        }

        detectedAssays <- names(currentS4@metabolite_data)
        normData$assay_names <- detectedAssays
        logInfoFn(paste("Detected assays:", paste(detectedAssays, collapse = ", ")))

        for (assayName in normData$assay_names) {
            if (!assayName %in% names(normData$itsd_selections)) {
                normData$itsd_selections[[assayName]] <- NULL
            }
        }

        invisible(detectedAssays)
    }, error = function(e) {
        logWarnFn(paste("Could not detect assay names:", e$message))
        invisible(NULL)
    })
}

# Keep the selected-tab auto pre-normalization QC shell top-level so a later
# extraction wave can move it without reopening the module server wrapper.
runMetabNormAutoPreNormalizationQcObserverShell <- function(
    selectedTab,
    workflowData,
    experimentPaths,
    normData,
    colorVariable = NULL,
    shapeVariable = NULL,
    addLogFn = function(entry) invisible(entry),
    reqFn = shiny::req,
    getStateFn = function(manager) manager$getState(),
    withProgressFn = shiny::withProgress,
    generatePreNormalizationQcFn = generateMetabNormPreNormalizationQc,
    getPlotAestheticsFn = getPlotAesthetics,
    generateMetabQcPlotsFn = generateMetabQcPlots,
    logInfoFn = logger::log_info,
    logWarnFn = logger::log_warn,
    logErrorFn = logger::log_error
) {
    if (is.null(selectedTab) || selectedTab != "norm") {
        return(invisible(NULL))
    }

    logInfoFn("Normalization tab selected - checking if pre-QC needed")

    if (isTRUE(normData$pre_norm_qc_generated)) {
        return(invisible(NULL))
    }

    reqFn(workflowData$state_manager)
    currentS4 <- tryCatch(getStateFn(workflowData$state_manager), error = function(e) NULL)

    if (is.null(currentS4) || !inherits(currentS4, "MetaboliteAssayData")) {
        return(invisible(NULL))
    }

    logInfoFn("Auto-triggering pre-normalization QC plots")

    withProgressFn(
        message = "Generating Pre-Normalization QC..."
        , value = 0.5
        , {
            generatePreNormalizationQcFn(
                workflowData = workflowData
                , experimentPaths = experimentPaths
                , normData = normData
                , getPlotAestheticsFn = function() {
                    getPlotAestheticsFn(
                        colorVariable = colorVariable
                        , shapeVariable = shapeVariable
                    )
                }
                , addLogFn = addLogFn
                , reqFn = reqFn
                , generateMetabQcPlotsFn = generateMetabQcPlotsFn
                , logInfoFn = logInfoFn
                , logWarnFn = logWarnFn
                , logErrorFn = logErrorFn
            )
        }
    )

    normData$pre_norm_qc_generated <- TRUE

    invisible(list(
        selectedTab = selectedTab
        , currentS4 = currentS4
        , preNormQcGenerated = normData$pre_norm_qc_generated
    ))
}

# Keep the dynamic ITSD selection tables render shell top-level so a later
# extraction wave can move it without reopening the module server wrapper.
renderMetabNormItsdSelectionUi <- function(
    normData,
    ns,
    renderUIFn = shiny::renderUI,
    reqFn = shiny::req,
    mapFn = purrr::map,
    wellPanelFn = shiny::wellPanel,
    h5Fn = shiny::h5,
    dataTableOutputFn = DT::dataTableOutput,
    brFn = shiny::br,
    tagListFn = shiny::tagList
) {
    renderUIFn({
        reqFn(normData$assay_names)

        assayUis <- mapFn(normData$assay_names, \(assayName) {
            safeName <- gsub("[^A-Za-z0-9]", "_", tolower(assayName))
            wellPanelFn(
                h5Fn(paste("Assay:", assayName))
                , dataTableOutputFn(ns(paste0("itsd_table_", safeName)))
                , brFn()
            )
        })

        tagListFn(assayUis)
    })
}

# Keep the per-assay ITSD selection table render shell top-level so a later
# extraction wave can move it without reopening the module server wrapper.
renderMetabNormItsdSelectionTable <- function(
    assayName,
    currentS4,
    renderDataTableFn = DT::renderDataTable,
    buildItsdSelectionTableFn = buildItsdSelectionTable,
    datatableFn = DT::datatable,
    formatStyleFn = DT::formatStyle,
    styleEqualFn = DT::styleEqual,
    formatRoundFn = DT::formatRound,
    datatableOptions = list(
        pageLength = 10,
        scrollX = TRUE,
        order = list(list(4, "desc"), list(3, "asc"))
    ),
    candidateHighlight = "#d4edda",
    roundedColumns = c("mean_intensity", "cv_percent"),
    digits = 2,
    filter = "top"
) {
    renderDataTableFn({
        assayData <- currentS4@metabolite_data[[assayName]]
        if (is.null(assayData)) {
            return(NULL)
        }

        metaboliteIdCol <- currentS4@metabolite_id_column
        annotationCol <- currentS4@annotation_id_column

        selectionTable <- buildItsdSelectionTableFn(
            assay_data = assayData,
            metabolite_id_col = metaboliteIdCol,
            annotation_cols = annotationCol
        )

        preselected <- which(selectionTable$is_candidate)

        datatableFn(
            selectionTable,
            selection = list(mode = "multiple", selected = preselected),
            filter = filter,
            options = datatableOptions,
            rownames = FALSE
        ) |>
            formatStyleFn(
                "is_candidate",
                backgroundColor = styleEqualFn(TRUE, candidateHighlight)
            ) |>
            formatRoundFn(columns = roundedColumns, digits = digits)
    })
}

# Keep the per-assay ITSD selection tracking registration top-level so a later
# extraction wave can move it without reopening the module server wrapper.
registerMetabNormItsdSelectionTracking <- function(
    normData,
    input,
    walkFn = purrr::walk,
    observeEventFn = shiny::observeEvent,
    selectionGetter = function(input, inputId) input[[inputId]],
    safeNameFn = function(assayName) gsub("[^A-Za-z0-9]", "_", tolower(assayName)),
    logInfoFn = logger::log_info
) {
    walkFn(normData$assay_names, \(assayName) {
        safeName <- safeNameFn(assayName)
        inputId <- paste0("itsd_table_", safeName, "_rows_selected")

        observeEventFn(selectionGetter(input, inputId), {
            selectedRows <- selectionGetter(input, inputId)
            normData$itsd_selections[[assayName]] <- selectedRows
            logInfoFn(paste(
                "ITSD selection updated for", assayName, ":"
                , length(selectedRows), "features selected"
            ))
        }, ignoreNULL = FALSE)
    })

    invisible(NULL)
}

# Keep the ITSD selection tracking observer shell top-level so a later
# extraction wave can move it without reopening the module server wrapper.
runMetabNormItsdSelectionTrackingObserverShell <- function(
    normData,
    input,
    reqFn = shiny::req,
    registerItsdSelectionTrackingFn = registerMetabNormItsdSelectionTracking
) {
    reqFn(normData$assay_names)

    registerItsdSelectionTrackingFn(
        normData = normData,
        input = input
    )

    invisible(NULL)
}

# Keep the per-assay ITSD selection table binding observer shell top-level so a
# later extraction wave can move it without reopening the module server wrapper.
runMetabNormItsdSelectionTableObserverShell <- function(
    normData,
    workflowData,
    output,
    reqFn = shiny::req,
    getStateFn = function(stateManager) stateManager$getState(),
    walkFn = purrr::walk,
    safeNameFn = function(assayName) gsub("[^A-Za-z0-9]", "_", tolower(assayName)),
    renderItsdSelectionTableFn = renderMetabNormItsdSelectionTable
) {
    reqFn(normData$assay_names)
    reqFn(workflowData$state_manager)

    currentS4 <- tryCatch({
        getStateFn(workflowData$state_manager)
    }, error = function(e) NULL)

    if (is.null(currentS4) || !inherits(currentS4, "MetaboliteAssayData")) {
        return(invisible(NULL))
    }

    walkFn(normData$assay_names, \(assayName) {
        safeName <- safeNameFn(assayName)
        outputId <- paste0("itsd_table_", safeName)

        output[[outputId]] <- renderItsdSelectionTableFn(
            assayName = assayName,
            currentS4 = currentS4
        )
    })

    invisible(NULL)
}

# Keep the dynamic RUV QC plots render shell top-level so a later extraction
# wave can move it without reopening the module server wrapper.
renderMetabNormRuvQcUi <- function(
    normData,
    ns,
    renderUIFn = shiny::renderUI,
    reqFn = shiny::req,
    mapFn = purrr::map,
    tagListFn = shiny::tagList,
    fluidRowFn = shiny::fluidRow,
    columnFn = shiny::column,
    h5Fn = shiny::h5,
    h6Fn = shiny::h6,
    jquiResizableFn = shinyjqui::jqui_resizable,
    plotOutputFn = shiny::plotOutput,
    wellPanelFn = shiny::wellPanel,
    verbatimTextOutputFn = shiny::verbatimTextOutput,
    brFn = shiny::br,
    dataTableOutputFn = DT::dataTableOutput,
    hrFn = shiny::hr
) {
    renderUIFn({
        reqFn(normData$assay_names)

        assayUis <- mapFn(normData$assay_names, \(assayName) {
            safeName <- gsub("[^A-Za-z0-9]", "_", tolower(assayName))

            tagListFn(
                fluidRowFn(
                    columnFn(
                        12,
                        h5Fn(
                            paste("Assay:", assayName),
                            style = "border-bottom: 1px solid #ddd; padding-bottom: 5px;"
                        )
                    )
                ),
                fluidRowFn(
                    columnFn(
                        8,
                        jquiResizableFn(
                            plotOutputFn(
                                ns(paste0("cancor_plot_", safeName)),
                                height = "400px"
                            )
                        )
                    ),
                    columnFn(
                        4,
                        wellPanelFn(
                            h6Fn("Optimization Summary"),
                            verbatimTextOutputFn(ns(paste0("ruv_summary_", safeName))),
                            brFn(),
                            h6Fn("Results Table"),
                            dataTableOutputFn(ns(paste0("ruv_table_", safeName)))
                        )
                    )
                ),
                hrFn()
            )
        })

        do.call(tagListFn, assayUis)
    })
}

# Keep the per-assay RUV cancor plot render shell top-level so a later
# extraction wave can move it without reopening the module server wrapper.
renderMetabNormRuvCancorPlot <- function(
    assayName,
    normData,
    renderPlotFn = shiny::renderPlot,
    buildFallbackPlotFn = function() {
        ggplot2::ggplot() +
            ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data", size = 6) +
            ggplot2::theme_void()
    }
) {
    renderPlotFn({
        result <- normData$ruv_optimization_results[[assayName]]
        if (!is.null(result) && isTRUE(result$success) && !is.null(result$cancor_plot)) {
            result$cancor_plot
        } else {
            buildFallbackPlotFn()
        }
    })
}

# Keep the per-assay RUV optimization summary render shell top-level so a later
# extraction wave can move it without reopening the module server wrapper.
renderMetabNormRuvOptimizationSummary <- function(
    assayName,
    normData,
    renderTextFn = shiny::renderText,
    pendingMessage = "Not yet computed"
) {
    renderTextFn({
        result <- normData$ruv_optimization_results[[assayName]]

        if (!is.null(result) && isTRUE(result$success)) {
            return(paste0(
                "Best k: ", result$best_k, "\n"
                , "Best %: ", result$best_percentage, "\n"
                , "Separation: ", round(result$separation_score, 4), "\n"
                , "Controls: ", sum(result$control_genes_index, na.rm = TRUE)
            ))
        }

        if (!is.null(result)) {
            return(paste0("Failed: ", result$error))
        }

        pendingMessage
    })
}

# Keep the per-assay RUV results table render shell top-level so a later
# extraction wave can move it without reopening the module server wrapper.
renderMetabNormRuvResultsTable <- function(
    assayName,
    normData,
    renderDataTableFn = DT::renderDataTable,
    datatableFn = DT::datatable,
    datatableOptions = list(pageLength = 5, dom = "t")
) {
    renderDataTableFn({
        result <- normData$ruv_optimization_results[[assayName]]

        if (!is.null(result) && isTRUE(result$success) && !is.null(result$optimization_results)) {
            return(datatableFn(
                result$optimization_results
                , options = datatableOptions
                , rownames = FALSE
            ))
        }

        NULL
    })
}

# Keep the per-assay RUV plot/summary/table binding observer shell top-level
# so a later extraction wave can move it without reopening the module server
# wrapper.
runMetabNormRuvBindingObserverShell <- function(
    normData,
    output,
    reqFn = shiny::req,
    walkFn = purrr::walk,
    safeNameFn = function(assayName) gsub("[^A-Za-z0-9]", "_", tolower(assayName)),
    renderRuvCancorPlotFn = renderMetabNormRuvCancorPlot,
    renderRuvOptimizationSummaryFn = renderMetabNormRuvOptimizationSummary,
    renderRuvResultsTableFn = renderMetabNormRuvResultsTable
) {
    reqFn(normData$assay_names)
    reqFn(length(normData$ruv_optimization_results) > 0)

    walkFn(normData$assay_names, \(assayName) {
        safeName <- safeNameFn(assayName)

        output[[paste0("cancor_plot_", safeName)]] <- renderRuvCancorPlotFn(
            assayName = assayName,
            normData = normData
        )

        output[[paste0("ruv_summary_", safeName)]] <- renderRuvOptimizationSummaryFn(
            assayName = assayName,
            normData = normData
        )

        output[[paste0("ruv_table_", safeName)]] <- renderRuvResultsTableFn(
            assayName = assayName,
            normData = normData
        )
    })

    invisible(NULL)
}

# Keep the static QC image output binding shell top-level so a later extraction
# wave can move it without reopening the module server wrapper.
runMetabNormQcImageBindingShell <- function(
    output,
    normData,
    qcDir,
    renderQcImageFn = renderMetabNormQcImageForAssay
) {
    output$pca_post_filter_assay1 <- renderQcImageFn(1, "pca", "pre_norm", normData, qcDir)
    output$pca_post_norm_assay1 <- renderQcImageFn(1, "pca", "post_norm", normData, qcDir)
    output$pca_ruv_corrected_assay1 <- renderQcImageFn(1, "pca", "ruv_corrected", normData, qcDir)
    output$pca_post_filter_assay2 <- renderQcImageFn(2, "pca", "pre_norm", normData, qcDir)
    output$pca_post_norm_assay2 <- renderQcImageFn(2, "pca", "post_norm", normData, qcDir)
    output$pca_ruv_corrected_assay2 <- renderQcImageFn(2, "pca", "ruv_corrected", normData, qcDir)

    output$density_post_filter_assay1 <- renderQcImageFn(1, "density", "pre_norm", normData, qcDir)
    output$density_post_norm_assay1 <- renderQcImageFn(1, "density", "post_norm", normData, qcDir)
    output$density_ruv_corrected_assay1 <- renderQcImageFn(1, "density", "ruv_corrected", normData, qcDir)
    output$density_post_filter_assay2 <- renderQcImageFn(2, "density", "pre_norm", normData, qcDir)
    output$density_post_norm_assay2 <- renderQcImageFn(2, "density", "post_norm", normData, qcDir)
    output$density_ruv_corrected_assay2 <- renderQcImageFn(2, "density", "ruv_corrected", normData, qcDir)

    output$rle_post_filter_assay1 <- renderQcImageFn(1, "rle", "pre_norm", normData, qcDir)
    output$rle_post_norm_assay1 <- renderQcImageFn(1, "rle", "post_norm", normData, qcDir)
    output$rle_ruv_corrected_assay1 <- renderQcImageFn(1, "rle", "ruv_corrected", normData, qcDir)
    output$rle_post_filter_assay2 <- renderQcImageFn(2, "rle", "pre_norm", normData, qcDir)
    output$rle_post_norm_assay2 <- renderQcImageFn(2, "rle", "post_norm", normData, qcDir)
    output$rle_ruv_corrected_assay2 <- renderQcImageFn(2, "rle", "ruv_corrected", normData, qcDir)

    output$correlation_post_filter_assay1 <- renderQcImageFn(1, "correlation", "pre_norm", normData, qcDir)
    output$correlation_post_norm_assay1 <- renderQcImageFn(1, "correlation", "post_norm", normData, qcDir)
    output$correlation_ruv_corrected_assay1 <- renderQcImageFn(1, "correlation", "ruv_corrected", normData, qcDir)
    output$correlation_post_filter_assay2 <- renderQcImageFn(2, "correlation", "pre_norm", normData, qcDir)
    output$correlation_post_norm_assay2 <- renderQcImageFn(2, "correlation", "post_norm", normData, qcDir)
    output$correlation_ruv_corrected_assay2 <- renderQcImageFn(2, "correlation", "ruv_corrected", normData, qcDir)

    invisible(NULL)
}

# Keep the static assay-label output binding shell top-level so a later
# extraction wave can move it without reopening the module server wrapper.
runMetabNormAssayLabelBindingShell <- function(
    output,
    getAssayNamesFn,
    renderAssayLabelFn = renderMetabNormAssayLabel
) {
    output$assay1_label_pca <- renderAssayLabelFn(1, getAssayNamesFn)
    output$assay2_label_pca <- renderAssayLabelFn(2, getAssayNamesFn)

    output$assay1_label_density <- renderAssayLabelFn(1, getAssayNamesFn)
    output$assay2_label_density <- renderAssayLabelFn(2, getAssayNamesFn)

    output$assay1_label_rle <- renderAssayLabelFn(1, getAssayNamesFn)
    output$assay2_label_rle <- renderAssayLabelFn(2, getAssayNamesFn)

    output$assay1_label_correlation <- renderAssayLabelFn(1, getAssayNamesFn)
    output$assay2_label_correlation <- renderAssayLabelFn(2, getAssayNamesFn)

    invisible(NULL)
}

