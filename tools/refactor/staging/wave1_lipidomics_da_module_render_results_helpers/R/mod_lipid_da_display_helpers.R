buildLipidDaStatusText <- function(
    analysisComplete,
    daResultsList,
    currentS4Object
) {
    if (isTRUE(analysisComplete)) {
        significantCounts <- daResultsList$significant_counts

        if (!is.null(significantCounts)) {
            assayLines <- vapply(
                names(significantCounts),
                function(assayName) {
                    counts <- significantCounts[[assayName]]
                    sprintf(
                        "%s:\n  Up: %d | Down: %d | NS: %d",
                        assayName, counts$up, counts$down, counts$ns
                    )
                },
                FUN.VALUE = character(1)
            )

            return(paste(c("Analysis Complete", "", assayLines), collapse = "\n"))
        }

        return("Analysis complete.")
    }

    if (!is.null(currentS4Object)) {
        return("Data loaded. Ready to run analysis.")
    }

    "Waiting for data.\nClick 'Load Filtered Session' to begin."
}

buildLipidDaContrastsText <- function(
    contrastsTbl,
    fallbackFormatter = function(x) capture.output(print(x))
) {
    if (is.null(contrastsTbl) || nrow(contrastsTbl) == 0) {
        return("No contrasts defined.\nLoad a filtered session or define contrasts in the design tab.")
    }

    if ("friendly_names" %in% colnames(contrastsTbl)) {
        return(paste(contrastsTbl$friendly_names, collapse = "\n"))
    }

    if ("contrasts" %in% colnames(contrastsTbl)) {
        return(paste(contrastsTbl$contrasts, collapse = "\n"))
    }

    paste(fallbackFormatter(contrastsTbl), collapse = "\n")
}

buildLipidDaVolcanoGlimmaUi <- function(
    daResultsList,
    selectedContrast,
    selectedAssay,
    daQValThresh,
    widgetFactory = generateLipidDAVolcanoPlotGlimma,
    logError = logger::log_error
) {
    shiny::req(daResultsList, selectedContrast)

    if (is.null(selectedAssay) || identical(selectedAssay, "Combined")) {
        return(shiny::div(
            class = "alert alert-info",
            style = "margin: 20px;",
            shiny::icon("info-circle"),
            shiny::tags$strong(" Combined View: "),
            "Interactive Glimma plots require a single assay selection. ",
            "Please select a specific assay (e.g., LCMS_Pos) or uncheck 'Interactive Plot' to use the static volcano plot for combined view."
        ))
    }

    tryCatch(
        {
            widget <- widgetFactory(
                da_results_list = daResultsList,
                selected_contrast = selectedContrast,
                selected_assay = selectedAssay,
                da_q_val_thresh = daQValThresh
            )

            if (is.null(widget)) {
                return(shiny::div(
                    class = "alert alert-warning",
                    "Could not generate Glimma plot. Try the static plot option."
                ))
            }

            widget
        },
        error = function(e) {
            logError(sprintf("Glimma error: %s", e$message))
            shiny::div(
                class = "alert alert-danger",
                sprintf("Error generating plot: %s", e$message)
            )
        }
    )
}

buildLipidDaHeatmapManualSaveWarning <- function(analysisComplete) {
    if (!isTRUE(analysisComplete)) {
        return(NULL)
    }

    shiny::div(
        class = "alert alert-info",
        style = "margin-bottom: 15px; padding: 10px; border-left: 5px solid #17a2b8;",
        shiny::icon("info-circle"),
        shiny::tags$b(" Note: "),
        "Heatmaps are NOT saved automatically during analysis. Please navigate to the ",
        shiny::tags$b("Heatmap"),
        " tab and use the ",
        shiny::tags$b("Save Heatmap"),
        " button to save your specific clustering configuration and image."
    )
}

buildLipidDaVolcanoStaticPlot <- function(
    daResultsList,
    selectedContrast,
    selectedAssay,
    daQValThresh,
    lfcThreshold,
    plotFactory = generateLipidDAVolcanoStatic
) {
    shiny::req(daResultsList, selectedContrast)

    plotFactory(
        da_results_list = daResultsList,
        selected_contrast = selectedContrast,
        selected_assay = selectedAssay,
        da_q_val_thresh = daQValThresh,
        lfc_threshold = lfcThreshold,
        show_labels = TRUE,
        n_labels = 15
    )
}

buildLipidDaHeatmapRenderState <- function(
    daResultsList,
    selectedContrast,
    selectedAssay,
    topN,
    heatmapClustering,
    scaleData,
    clusteringMethod,
    distanceMethod,
    colorScheme,
    showLipidNames,
    daQValThresh,
    treeCutMethod,
    nClusters,
    cutHeight,
    minClusterSize,
    heatmapFactory = generateLipidDAHeatmap
) {
    shiny::req(daResultsList, selectedContrast)

    hm <- heatmapFactory(
        da_results_list = daResultsList,
        selected_contrast = selectedContrast,
        selected_assay = selectedAssay,
        top_n = topN,
        clustering_method = clusteringMethod,
        distance_method = distanceMethod,
        cluster_rows = heatmapClustering %in% c("both", "row"),
        cluster_cols = heatmapClustering %in% c("both", "column"),
        scale_data = scaleData,
        color_scheme = colorScheme,
        show_lipid_names = showLipidNames,
        da_q_val_thresh = daQValThresh,
        tree_cut_method = treeCutMethod,
        n_clusters = nClusters,
        cut_height = cutHeight,
        min_cluster_size = minClusterSize
    )

    if (is.null(hm)) {
        return(NULL)
    }

    if (is.list(hm) && "plot" %in% names(hm)) {
        return(list(
            plot = hm$plot,
            rowClusters = hm$row_clusters,
            colClusters = hm$col_clusters,
            currentHeatmapPlot = hm$plot
        ))
    }

    list(
        plot = hm,
        rowClusters = NULL,
        colClusters = NULL,
        currentHeatmapPlot = hm
    )
}

storeLipidDaHeatmapRenderState <- function(heatmapState, daData) {
    if (is.null(heatmapState)) {
        return(NULL)
    }

    daData$current_row_clusters <- heatmapState$rowClusters
    daData$current_col_clusters <- heatmapState$colClusters
    daData$current_heatmap_plot <- heatmapState$currentHeatmapPlot

    heatmapState$plot
}

buildLipidDaClusterSummaryText <- function(rowClusters) {
    if (is.null(rowClusters)) {
        return("No clusters defined. Enable clustering and tree cutting on the heatmap.")
    }

    uniqueClusters <- sort(unique(rowClusters))
    summaryLines <- c(
        sprintf("Total Clusters: %d", length(uniqueClusters)),
        "----------------------------------"
    )

    for (clusterId in uniqueClusters) {
        members <- names(rowClusters)[rowClusters == clusterId]
        memberLine <- paste(head(members, 20), collapse = ", ")

        if (length(members) > 20) {
            memberLine <- paste0(
                memberLine,
                sprintf(", ... and %d more", length(members) - 20)
            )
        }

        summaryLines <- c(
            summaryLines,
            "",
            sprintf("Cluster %d (%d lipids):", clusterId, length(members)),
            memberLine
        )
    }

    paste(summaryLines, collapse = "\n")
}

buildLipidDaClusterSummaryRenderText <- function(
    treeCutMethod,
    rowClusters,
    reqFn = shiny::req,
    buildClusterSummaryTextFn = buildLipidDaClusterSummaryText
) {
    reqFn(treeCutMethod != "none")

    buildClusterSummaryTextFn(rowClusters)
}

