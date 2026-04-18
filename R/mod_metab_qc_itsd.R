# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# ============================================================================
# mod_metab_qc_itsd.R
# ============================================================================
# Purpose: Internal standard (ITSD) QC visualization Shiny module
#
# This module provides visualization of internal standard performance metrics
# including CV distribution and intensity trends across samples.
# ============================================================================

#' @title Internal Standard QC Visualization Module
#' @description A Shiny module for visualizing internal standard quality metrics.
#'              Displays IS detection, CV distributions, and intensity trends.
#' @name mod_metab_qc_itsd
NULL


analyzeMetabQcItsdData <- function(
    currentS4,
    inputPattern = NULL,
    getInternalStandardMetricsFn = getInternalStandardMetrics,
    getMetaboliteQuantDataFn = getMetaboliteQuantData,
    pivotLongerFn = tidyr::pivot_longer,
    allOfFn = dplyr::all_of,
    logInfoFn = logger::log_info
) {
    if (!inherits(currentS4, "MetaboliteAssayData")) {
        stop("Current state is not a MetaboliteAssayData object")
    }

    isPattern <- if (!is.null(inputPattern) && nzchar(inputPattern)) {
        inputPattern
    } else if (!is.na(currentS4@internal_standard_regex) && nzchar(currentS4@internal_standard_regex)) {
        currentS4@internal_standard_regex
    } else {
        NULL
    }

    if (is.null(isPattern) || !nzchar(isPattern)) {
        stop("No internal standard pattern provided. Please enter a regex pattern to identify IS.")
    }

    logInfoFn(paste("Analyzing internal standards with pattern:", isPattern))

    metaboliteIdCol <- currentS4@metabolite_id_column
    annotationIdCol <- currentS4@annotation_id_column
    sampleIdCol <- currentS4@sample_id
    assayList <- currentS4@metabolite_data
    assayNames <- names(assayList)
    colsToTry <- unique(c(
        metaboliteIdCol,
        annotationIdCol,
        "Metabolite name",
        "Name",
        "metabolite"
    ))

    metricsList <- list()
    longDataList <- list()

    for (i in seq_along(assayList)) {
        assayData <- assayList[[i]]
        assayName <- if (!is.null(assayNames)) assayNames[[i]] else paste0("Assay_", i)
        isMet <- data.frame(is_id = character(), mean_intensity = numeric(), cv = numeric())
        usedCol <- NULL

        for (col in colsToTry) {
            if (!is.null(col) && col %in% colnames(assayData)) {
                isMet <- getInternalStandardMetricsFn(
                    assay_data = assayData,
                    is_pattern = isPattern,
                    metabolite_id_col = col,
                    sample_id_col = sampleIdCol
                )

                if (nrow(isMet) > 0) {
                    usedCol <- col
                    logInfoFn(paste("Found IS in column:", col))
                    break
                }
            }
        }

        if (nrow(isMet) == 0) {
            next
        }

        isMet$assay <- assayName
        metricsList[[assayName]] <- isMet

        idColForFilter <- if (!is.null(usedCol)) usedCol else metaboliteIdCol
        quantInfo <- getMetaboliteQuantDataFn(assayData)
        sampleCols <- quantInfo$sample_names

        if (length(sampleCols) == 0 || !idColForFilter %in% colnames(assayData)) {
            next
        }

        isIds <- isMet$is_id
        isRows <- assayData[[idColForFilter]] %in% isIds

        if (!any(isRows)) {
            next
        }

        isData <- assayData[isRows, c(idColForFilter, sampleCols), drop = FALSE]
        isLong <- pivotLongerFn(
            isData,
            cols = allOfFn(sampleCols),
            names_to = "Sample",
            values_to = "Intensity"
        )
        isLong$assay <- assayName
        names(isLong)[1] <- "IS_ID"
        longDataList[[assayName]] <- isLong
    }

    if (length(metricsList) == 0) {
        firstAssayCols <- if (length(assayList) > 0) colnames(assayList[[1]]) else character()
        colsSearched <- colsToTry[colsToTry %in% firstAssayCols]
        stop(paste0(
            "No internal standards found. Pattern: ", isPattern,
            "\nColumns searched: ", paste(colsSearched, collapse = ", "),
            "\nTip: Check if your data has a column with metabolite names (e.g., 'Metabolite name')"
        ))
    }

    allMetrics <- do.call(rbind, metricsList)
    allLong <- if (length(longDataList) > 0) do.call(rbind, longDataList) else NULL
    nIsTotal <- nrow(allMetrics)
    medianCv <- median(allMetrics$cv, na.rm = TRUE)
    maxCv <- max(allMetrics$cv, na.rm = TRUE)
    resultText <- paste(
        "Internal Standard Analysis Complete",
        "====================================",
        sprintf("Pattern used: %s", isPattern),
        sprintf("Total IS detected: %d", nIsTotal),
        "",
        "CV Statistics:",
        sprintf("  Median CV: %.1f%%", medianCv),
        sprintf("  Max CV: %.1f%%", maxCv),
        sprintf("  IS with CV > 30%%: %d", sum(allMetrics$cv > 30, na.rm = TRUE)),
        "",
        "Per-Assay IS Counts:",
        paste(vapply(names(metricsList), function(name) {
            sprintf("  %s: %d internal standards", name, nrow(metricsList[[name]]))
        }, character(1)), collapse = "\n"),
        sep = "\n"
    )

    list(
        metrics = allMetrics,
        longData = allLong,
        resultText = resultText,
        nIsTotal = nIsTotal,
        pattern = isPattern,
        metricsByAssay = metricsList
    )
}

buildMetabQcItsdSummaryUi <- function(
    metrics,
    splitFn = split,
    medianFn = stats::median,
    paragraphFn = shiny::p,
    iconFn = shiny::icon,
    listTagFn = shiny::tags$ul,
    itemTagFn = shiny::tags$li,
    horizontalRuleFn = shiny::hr,
    smallTagFn = shiny::tags$small,
    tagListFn = shiny::tagList,
    sprintfFn = sprintf
) {
    if (is.null(metrics)) {
        return(paragraphFn(
            iconFn("info-circle"),
            " Click 'Analyze' to detect internal standards.",
            style = "color: #666;"
        ))
    }

    summaryByAssay <- splitFn(metrics, metrics$assay)

    summaryItems <- lapply(names(summaryByAssay), function(assayName) {
        assayMetrics <- summaryByAssay[[assayName]]
        nIs <- nrow(assayMetrics)
        medianCv <- medianFn(assayMetrics$cv, na.rm = TRUE)

        cvColor <- if (medianCv <= 15) "green" else if (medianCv <= 30) "orange" else "red"
        cvIcon <- if (medianCv <= 15) "check-circle" else if (medianCv <= 30) "exclamation-circle" else "times-circle"

        itemTagFn(
            iconFn(cvIcon, style = paste("color:", cvColor)),
            sprintfFn(" %s: %d IS (median CV: %.1f%%)", assayName, nIs, medianCv)
        )
    })

    tagListFn(
        listTagFn(summaryItems, style = "list-style: none; padding-left: 0;"),
        horizontalRuleFn(),
        smallTagFn(
            iconFn("info-circle"),
            " CV < 15%: Good | 15-30%: Acceptable | > 30%: Review",
            style = "color: #666;"
        )
    )
}

# Keep the visualization-tab assembly isolated so the wrapper can shed
# renderUI wiring without changing the public server entry point.
buildMetabQcItsdVizTabsUi <- function(
    metrics,
    ns,
    tabPanelFn = shiny::tabPanel,
    brFn = shiny::br,
    jquiResizableFn = shinyjqui::jqui_resizable,
    plotOutputFn = shiny::plotOutput,
    tabsetPanelFn = shiny::tabsetPanel
) {
    if (is.null(metrics)) {
        return(NULL)
    }

    tabList <- list(
        tabPanelFn(
            "CV Distribution",
            brFn(),
            jquiResizableFn(
                plotOutputFn(ns("cv_plot"), height = "500px")
            )
        ),
        tabPanelFn(
            "Intensity Trends",
            brFn(),
            jquiResizableFn(
                plotOutputFn(ns("intensity_plot"), height = "500px")
            )
        )
    )

    do.call(tabsetPanelFn, c(list(id = ns("is_viz_tabset")), tabList))
}

# Keep the CV lollipop-plot assembly isolated so the wrapper can shed
# renderPlot logic without changing the public server entry point.
buildMetabQcItsdCvPlot <- function(metrics) {
    metrics$is_id <- stats::reorder(metrics$is_id, metrics$cv)

    metrics$cv_status <- dplyr::case_when(
        metrics$cv <= 15 ~ "Good (<15%)"
        , metrics$cv <= 30 ~ "Acceptable (15-30%)"
        , TRUE ~ "Review (>30%)"
    )
    metrics$cv_status <- factor(
        metrics$cv_status,
        levels = c("Good (<15%)", "Acceptable (15-30%)", "Review (>30%)")
    )

    ggplot2::ggplot(metrics, ggplot2::aes(x = is_id, y = cv)) +
        ggplot2::geom_hline(
            yintercept = 15,
            linetype = "dashed",
            color = "#2ecc71",
            linewidth = 0.8,
            alpha = 0.7
        ) +
        ggplot2::geom_hline(
            yintercept = 30,
            linetype = "dashed",
            color = "#e67e22",
            linewidth = 0.8,
            alpha = 0.7
        ) +
        ggplot2::geom_segment(
            ggplot2::aes(xend = is_id, yend = 0),
            color = "grey60",
            linewidth = 0.8
        ) +
        ggplot2::geom_point(ggplot2::aes(color = cv_status), size = 4) +
        ggplot2::facet_wrap(~assay, scales = "free_y") +
        ggplot2::scale_color_manual(
            values = c(
                "Good (<15%)" = "#2ecc71",
                "Acceptable (15-30%)" = "#e67e22",
                "Review (>30%)" = "#e74c3c"
            ),
            name = "CV Status"
        ) +
        ggplot2::labs(
            title = "Internal Standard CV",
            subtitle = "Dashed lines: 15% (green) and 30% (orange) thresholds",
            x = NULL,
            y = "CV (%)"
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
            legend.position = "bottom",
            panel.grid.major.x = ggplot2::element_blank(),
            strip.text = ggplot2::element_text(face = "bold")
        ) +
        ggplot2::coord_flip()
}

# Keep the intensity-trend plot assembly isolated so the wrapper can shed
# renderPlot logic without changing the public server entry point.
buildMetabQcItsdIntensityPlot <- function(longData) {
    ggplot2::ggplot(
        longData,
        ggplot2::aes(x = Sample, y = Intensity, color = IS_ID, group = IS_ID)
    ) +
        ggplot2::geom_line(alpha = 0.7) +
        ggplot2::geom_point(size = 1, alpha = 0.5) +
        ggplot2::facet_wrap(~assay, scales = "free_y") +
        ggplot2::labs(
            title = "Internal Standard Intensity Across Samples",
            x = "Sample",
            y = "Intensity",
            color = "Internal Standard"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(
                angle = 90,
                hjust = 1,
                vjust = 0.5,
                size = 6
            ),
            legend.position = "bottom"
        )
}

# Keep the module server body isolated so the public wrapper can shed
# observer and render wiring without changing the entrypoint signature.
runMetabQcItsdServerBody <- function(
    input,
    output,
    session,
    workflowData,
    omicType = NULL,
    experimentLabel = NULL,
    reactiveValFn = shiny::reactiveVal,
    observeEventFn = shiny::observeEvent,
    reqFn = shiny::req,
    showNotificationFn = shiny::showNotification,
    removeNotificationFn = shiny::removeNotification,
    renderTextFn = shiny::renderText,
    renderUiFn = shiny::renderUI,
    renderPlotFn = shiny::renderPlot,
    analyzeMetabQcItsdDataFn = analyzeMetabQcItsdData,
    buildMetabQcItsdSummaryUiFn = buildMetabQcItsdSummaryUi,
    buildMetabQcItsdVizTabsUiFn = buildMetabQcItsdVizTabsUi,
    buildMetabQcItsdCvPlotFn = buildMetabQcItsdCvPlot,
    buildMetabQcItsdIntensityPlotFn = buildMetabQcItsdIntensityPlot,
    logInfoFn = logger::log_info,
    logErrorFn = logger::log_error
) {
    ns <- session$ns

    is_metrics <- reactiveValFn(NULL)
    is_data_long <- reactiveValFn(NULL)

    # Analyze internal standards
    observeEventFn(input$analyze_is, {
        reqFn(workflowData$state_manager)

        showNotificationFn(
            "Analyzing internal standards..."
            , id = "is_analysis_working"
            , duration = NULL
        )

        tryCatch({
            current_s4 <- workflowData$state_manager$getState()
            reqFn(current_s4)
            analysis <- analyzeMetabQcItsdDataFn(
                currentS4 = current_s4,
                inputPattern = input$is_pattern
            )
            is_metrics(analysis$metrics)
            is_data_long(analysis$longData)
            output$is_results <- renderTextFn(analysis$resultText)

            logInfoFn(paste("IS analysis complete:", analysis$nIsTotal, "standards found"))

            removeNotificationFn("is_analysis_working")
            showNotificationFn(
                sprintf("Found %d internal standards", analysis$nIsTotal)
                , type = "message"
            )

        }, error = function(e) {
            msg <- paste("Error analyzing internal standards:", e$message)
            logErrorFn(msg)
            showNotificationFn(msg, type = "error", duration = 10)
            removeNotificationFn("is_analysis_working")
        })
    })

    # Render IS summary
    output$is_summary <- renderUiFn({
        buildMetabQcItsdSummaryUiFn(metrics = is_metrics())
    })

    # Render IS visualization tabs
    output$is_viz_tabs <- renderUiFn({
        buildMetabQcItsdVizTabsUiFn(
            metrics = is_metrics(),
            ns = ns
        )
    })

    # Render CV distribution plot - lollipop chart
    output$cv_plot <- renderPlotFn({
        metrics <- is_metrics()
        reqFn(metrics)
        buildMetabQcItsdCvPlotFn(metrics = metrics)
    })

    # Render intensity trend plot
    output$intensity_plot <- renderPlotFn({
        long_data <- is_data_long()
        reqFn(long_data)
        buildMetabQcItsdIntensityPlotFn(longData = long_data)
    })

    invisible(NULL)
}

