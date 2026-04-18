filterMetabDaDisplayResults <- function(
    results,
    selectedContrast = NULL,
    selectedAssay = NULL,
    significanceFilter = "all",
    daQValThresh = 0.05,
    treatLfcCutoff = 0,
    maxRows = NULL
) {
    if (is.null(results) || nrow(results) == 0) {
        return(results)
    }

    if (!is.null(selectedContrast)) {
        comparisonMatches <- if ("comparison" %in% colnames(results)) {
            results$comparison == selectedContrast
        } else {
            rep(FALSE, nrow(results))
        }
        friendlyNameMatches <- if ("friendly_name" %in% colnames(results)) {
            results$friendly_name == selectedContrast
        } else {
            rep(FALSE, nrow(results))
        }

        results <- results[comparisonMatches | friendlyNameMatches, , drop = FALSE]
    }

    if (!is.null(selectedAssay) && selectedAssay != "All") {
        results <- results[results$assay == selectedAssay, , drop = FALSE]
    }

    if (significanceFilter == "significant") {
        results <- results[results$fdr_qvalue < daQValThresh, , drop = FALSE]
    } else if (significanceFilter == "up") {
        results <- results[
            results$fdr_qvalue < daQValThresh &
                results$logFC > treatLfcCutoff,
            ,
            drop = FALSE
        ]
    } else if (significanceFilter == "down") {
        results <- results[
            results$fdr_qvalue < daQValThresh &
                results$logFC < -treatLfcCutoff,
            ,
            drop = FALSE
        ]
    }

    if (!is.null(maxRows) && nrow(results) > maxRows) {
        results <- results[seq_len(maxRows), , drop = FALSE]
    }

    results
}

summarizeMetabDaDisplayResults <- function(results, daQValThresh = 0.05) {
    if (is.null(results) || nrow(results) == 0) {
        return(NULL)
    }

    total <- nrow(results)
    significant <- sum(results$significant != "NS", na.rm = TRUE)
    upRegulated <- sum(results$significant == "Up", na.rm = TRUE)
    downRegulated <- sum(results$significant == "Down", na.rm = TRUE)

    list(
        total = total,
        significant = significant,
        upRegulated = upRegulated,
        downRegulated = downRegulated,
        significantPct = 100 * significant / max(total, 1),
        qValueThreshold = daQValThresh
    )
}

buildMetabDaSummaryStatsText <- function(results, daQValThresh = 0.05) {
    if (is.null(results) || nrow(results) == 0) {
        return("No results available.")
    }

    summaryStats <- summarizeMetabDaDisplayResults(
        results = results,
        daQValThresh = daQValThresh
    )
    if (is.null(summaryStats)) {
        return("No results available.")
    }

    paste(
        c(
            sprintf("Total metabolites: %d", summaryStats$total),
            sprintf(
                "Significant (Q < %.3f): %d (%.1f%%)",
                summaryStats$qValueThreshold,
                summaryStats$significant,
                summaryStats$significantPct
            ),
            sprintf("  Up-regulated: %d", summaryStats$upRegulated),
            sprintf("  Down-regulated: %d", summaryStats$downRegulated)
        ),
        collapse = "\n"
    )
}

buildMetabDaSummaryStatsRenderOutput <- function(
    daResultsList,
    selectedContrast = NULL,
    selectedAssay = NULL,
    daQValThresh = 0.05,
    requireInputs = shiny::req,
    filterResults = filterMetabDaDisplayResults,
    buildSummaryText = buildMetabDaSummaryStatsText,
    writeOutput = cat
) {
    requireInputs(daResultsList)

    results <- filterResults(
        results = daResultsList$da_metabolites_long,
        selectedContrast = selectedContrast,
        selectedAssay = selectedAssay
    )

    writeOutput(buildSummaryText(
        results = results,
        daQValThresh = daQValThresh
    ))
}

shapeMetabDaTableDisplayResults <- function(results) {
    if (is.null(results)) {
        return(results)
    }

    displayCols <- c(
        "metabolite_id", "metabolite_name", "assay", "logFC",
        "raw_pvalue", "fdr_qvalue", "significant"
    )
    displayCols <- intersect(displayCols, colnames(results))

    results[, displayCols, drop = FALSE]
}

buildMetabDaResultsDatatable <- function(results) {
    if (is.null(results) || nrow(results) == 0) {
        return(NULL)
    }

    DT::datatable(
        results,
        options = list(
            pageLength = 25,
            scrollX = TRUE,
            dom = "Bfrtip",
            buttons = c("copy", "csv", "excel")
        ),
        extensions = "Buttons",
        rownames = FALSE,
        filter = "top"
    ) |>
        DT::formatRound(columns = c("logFC"), digits = 3) |>
        DT::formatRound(columns = c("raw_pvalue", "fdr_qvalue"), digits = 6) |>
        DT::formatStyle(
            "significant",
            backgroundColor = DT::styleEqual(
                c("Up", "Down", "NS"),
                c("#ffcccc", "#cce5ff", "white")
            )
        )
}

buildMetabDaResultsTableRenderOutput <- function(
    daResultsList,
    selectedContrast = NULL,
    selectedAssay = NULL,
    significanceFilter = "all",
    daQValThresh = 0.05,
    treatLfcCutoff = 0,
    maxRows = NULL,
    requireInputs = shiny::req,
    filterResults = filterMetabDaDisplayResults,
    shapeResults = shapeMetabDaTableDisplayResults,
    buildDatatable = buildMetabDaResultsDatatable
) {
    requireInputs(daResultsList)

    results <- filterResults(
        results = daResultsList$da_metabolites_long,
        selectedContrast = selectedContrast,
        selectedAssay = selectedAssay,
        significanceFilter = significanceFilter,
        daQValThresh = daQValThresh,
        treatLfcCutoff = treatLfcCutoff,
        maxRows = maxRows
    )
    if (is.null(results) || nrow(results) == 0) {
        return(NULL)
    }

    buildDatatable(shapeResults(results))
}

buildMetabDaClusterSummaryText <- function(clusters, maxMembers = 20) {
    if (is.null(clusters)) {
        return("No clusters defined. Enable clustering and tree cutting on the heatmap.")
    }

    clusterIds <- sort(unique(clusters))
    clusterSections <- vapply(
        clusterIds,
        function(clusterId) {
            members <- names(clusters)[clusters == clusterId]
            memberSummary <- paste(head(members, maxMembers), collapse = ", ")

            if (length(members) > maxMembers) {
                memberSummary <- paste0(
                    memberSummary,
                    sprintf(", ... and %d more", length(members) - maxMembers)
                )
            }

            paste0(
                "\nCluster ", clusterId, " (", length(members), " metabolites):\n",
                memberSummary,
                "\n"
            )
        },
        character(1)
    )

    paste0(
        "Total Clusters: ", length(clusterIds), "\n",
        "----------------------------------\n",
        paste(clusterSections, collapse = "")
    )
}

buildMetabDaClusterSummaryRenderOutput <- function(
    treeCutMethod,
    clusters,
    requireInputs = shiny::req,
    buildSummaryText = buildMetabDaClusterSummaryText,
    writeOutput = cat
) {
    requireInputs(treeCutMethod != "none")

    writeOutput(buildSummaryText(clusters))
}

buildMetabDaHeatmapManualSaveWarning <- function(analysisComplete) {
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

buildMetabDaHeatmapManualSaveWarningRenderOutput <- function(
    analysisComplete = FALSE,
    buildWarning = buildMetabDaHeatmapManualSaveWarning
) {
    buildWarning(analysisComplete)
}

buildMetabDaContrastsDisplayText <- function(contrastsTbl) {
    if (is.null(contrastsTbl) || nrow(contrastsTbl) == 0) {
        return("No contrasts defined.\nLoad a filtered session or define contrasts in the design tab.")
    }

    if ("friendly_names" %in% colnames(contrastsTbl)) {
        return(paste(contrastsTbl$friendly_names, collapse = "\n"))
    }

    if ("contrasts" %in% colnames(contrastsTbl)) {
        return(paste(contrastsTbl$contrasts, collapse = "\n"))
    }

    paste(capture.output(print(contrastsTbl)), collapse = "\n")
}

buildMetabDaContrastsRenderOutput <- function(
    contrastsTbl,
    buildDisplayText = buildMetabDaContrastsDisplayText,
    writeOutput = function(text, ...) cat(text)
) {
    writeOutput(buildDisplayText(contrastsTbl))
}

buildMetabDaStatusText <- function(
    analysisComplete = FALSE,
    currentS4Object = NULL,
    significantCounts = NULL
) {
    if (isTRUE(analysisComplete)) {
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
                character(1)
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

buildMetabDaStatusRenderOutput <- function(
    daResultsList = NULL,
    analysisComplete = FALSE,
    currentS4Object = NULL,
    buildStatusText = buildMetabDaStatusText,
    writeOutput = function(text, ...) cat(text)
) {
    significantCounts <- NULL
    if (!is.null(daResultsList)) {
        significantCounts <- daResultsList$significant_counts
    }

    writeOutput(
        buildStatusText(
            analysisComplete = analysisComplete,
            currentS4Object = currentS4Object,
            significantCounts = significantCounts
        )
    )
}

buildMetabDaGlimmaCombinedViewInfoBanner <- function(selectedAssay) {
    if (!is.null(selectedAssay) && selectedAssay != "Combined") {
        return(NULL)
    }

    shiny::div(
        class = "alert alert-info",
        style = "margin: 20px;",
        shiny::icon("info-circle"),
        shiny::tags$strong(" Combined View: "),
        "Interactive Glimma plots require a single assay selection. ",
        "Please select a specific assay (e.g., LCMS_Pos) or uncheck 'Interactive Plot' to use the static volcano plot for combined view."
    )
}

resolveMetabDaGlimmaWidgetOutput <- function(widget) {
    if (!is.null(widget)) {
        return(widget)
    }

    shiny::div(
        class = "alert alert-warning",
        "Could not generate Glimma plot. Try the static plot option."
    )
}

buildMetabDaGlimmaErrorBanner <- function(errorMessage) {
    shiny::div(
        class = "alert alert-danger",
        sprintf("Error generating plot: %s", errorMessage)
    )
}

buildMetabDaGlimmaRenderOutput <- function(
    daResultsList,
    selectedContrast,
    selectedAssay = NULL,
    daQValThresh = 0.05,
    requireInputs = shiny::req,
    buildCombinedViewInfo = buildMetabDaGlimmaCombinedViewInfoBanner,
    generatePlot = generateMetabDAVolcanoPlotGlimma,
    resolveWidgetOutput = resolveMetabDaGlimmaWidgetOutput,
    buildErrorBanner = buildMetabDaGlimmaErrorBanner,
    logError = logger::log_error
) {
    requireInputs(daResultsList, selectedContrast)

    combinedViewInfo <- buildCombinedViewInfo(selectedAssay)
    if (!is.null(combinedViewInfo)) {
        return(combinedViewInfo)
    }

    tryCatch(
        {
            widget <- generatePlot(
                da_results_list = daResultsList,
                selected_contrast = selectedContrast,
                selected_assay = selectedAssay,
                da_q_val_thresh = daQValThresh
            )

            resolveWidgetOutput(widget)
        },
        error = function(e) {
            logError(sprintf("Glimma error: %s", e$message))
            buildErrorBanner(e$message)
        }
    )
}

buildMetabDaStaticVolcanoPlot <- function(
    daResultsList,
    selectedContrast,
    selectedAssay = NULL,
    daQValThresh = 0.05,
    treatLfcCutoff = 0
) {
    generateMetabDAVolcanoStatic(
        da_results_list = daResultsList,
        selected_contrast = selectedContrast,
        selected_assay = selectedAssay,
        da_q_val_thresh = daQValThresh,
        lfc_threshold = treatLfcCutoff,
        show_labels = TRUE,
        n_labels = 15
    )
}

buildMetabDaStaticVolcanoRenderOutput <- function(
    daResultsList,
    selectedContrast,
    selectedAssay = NULL,
    daQValThresh = 0.05,
    treatLfcCutoff = 0,
    requireInputs = shiny::req,
    buildPlot = buildMetabDaStaticVolcanoPlot
) {
    requireInputs(daResultsList, selectedContrast)

    buildPlot(
        daResultsList = daResultsList,
        selectedContrast = selectedContrast,
        selectedAssay = selectedAssay,
        daQValThresh = daQValThresh,
        treatLfcCutoff = treatLfcCutoff
    )
}

buildMetabDaHeatmapPlotOutput <- function(
    daResultsList,
    selectedContrast,
    selectedAssay = NULL,
    topN = 50,
    clusteringMethod = "complete",
    distanceMethod = "euclidean",
    clusterRows = TRUE,
    clusterCols = TRUE,
    scaleData = "row",
    colorScheme = "RdBu",
    showMetaboliteNames = TRUE,
    daQValThresh = 0.05,
    treeCutMethod = "none",
    nClusters = NULL,
    cutHeight = NULL,
    minClusterSize = 1,
    daData = NULL
) {
    hm <- generateMetabDAHeatmap(
        da_results_list = daResultsList,
        selected_contrast = selectedContrast,
        selected_assay = selectedAssay,
        top_n = topN,
        clustering_method = clusteringMethod,
        distance_method = distanceMethod,
        cluster_rows = clusterRows,
        cluster_cols = clusterCols,
        scale_data = scaleData,
        color_scheme = colorScheme,
        show_metabolite_names = showMetaboliteNames,
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
        if (!is.null(daData)) {
            daData$current_row_clusters <- hm$row_clusters
            daData$current_col_clusters <- hm$col_clusters
            daData$current_heatmap_plot <- hm$plot
        }

        return(hm$plot)
    }

    hm
}

buildMetabDaHeatmapRenderOutput <- function(
    daResultsList,
    selectedContrast,
    selectedAssay = NULL,
    topN = 50,
    clusteringMethod = "complete",
    distanceMethod = "euclidean",
    heatmapClustering = "both",
    scaleData = "row",
    colorScheme = "RdBu",
    showMetaboliteNames = TRUE,
    daQValThresh = 0.05,
    treeCutMethod = "none",
    nClusters = NULL,
    cutHeight = NULL,
    minClusterSize = 1,
    daData = NULL,
    requireInputs = shiny::req,
    buildPlot = buildMetabDaHeatmapPlotOutput
) {
    requireInputs(daResultsList, selectedContrast)

    buildPlot(
        daResultsList = daResultsList,
        selectedContrast = selectedContrast,
        selectedAssay = selectedAssay,
        topN = topN,
        clusteringMethod = clusteringMethod,
        distanceMethod = distanceMethod,
        clusterRows = heatmapClustering %in% c("both", "row"),
        clusterCols = heatmapClustering %in% c("both", "column"),
        scaleData = scaleData,
        colorScheme = colorScheme,
        showMetaboliteNames = showMetaboliteNames,
        daQValThresh = daQValThresh,
        treeCutMethod = treeCutMethod,
        nClusters = nClusters,
        cutHeight = cutHeight,
        minClusterSize = minClusterSize,
        daData = daData
    )
}

