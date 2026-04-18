registerMetabDaOverviewOutputs <- function(
    output,
    workflowData,
    daData,
    renderPrint = shiny::renderPrint,
    buildContrastsOutput = buildMetabDaContrastsRenderOutput,
    buildStatusOutput = buildMetabDaStatusRenderOutput
) {
    output$contrasts_display <- renderPrint({
        buildContrastsOutput(
            contrastsTbl = workflowData$contrasts_tbl
        )
    })

    output$da_status <- renderPrint({
        buildStatusOutput(
            daResultsList = daData$da_results_list,
            analysisComplete = daData$analysis_complete,
            currentS4Object = daData$current_s4_object
        )
    })

    invisible(output)
}

registerMetabDaVisualizationOutputs <- function(
    input,
    output,
    daData,
    renderUi = shiny::renderUI,
    renderPlot = shiny::renderPlot,
    renderPrint = shiny::renderPrint,
    buildHeatmapWarningOutput = buildMetabDaHeatmapManualSaveWarningRenderOutput,
    buildGlimmaOutput = buildMetabDaGlimmaRenderOutput,
    buildStaticVolcanoOutput = buildMetabDaStaticVolcanoRenderOutput,
    buildHeatmapOutput = buildMetabDaHeatmapRenderOutput,
    buildClusterSummaryOutput = buildMetabDaClusterSummaryRenderOutput
) {
    output$heatmap_manual_save_warning <- renderUi({
        buildHeatmapWarningOutput(
            analysisComplete = daData$analysis_complete
        )
    })

    output$volcano_glimma <- renderUi({
        buildGlimmaOutput(
            daResultsList = daData$da_results_list,
            selectedContrast = input$volcano_contrast,
            selectedAssay = input$volcano_assay,
            daQValThresh = input$da_q_val_thresh
        )
    })

    output$volcano_static <- renderPlot({
        buildStaticVolcanoOutput(
            daResultsList = daData$da_results_list,
            selectedContrast = input$volcano_contrast,
            selectedAssay = input$volcano_assay,
            daQValThresh = input$da_q_val_thresh,
            treatLfcCutoff = input$treat_lfc_cutoff
        )
    })

    output$heatmap_plot <- renderPlot({
        buildHeatmapOutput(
            daResultsList = daData$da_results_list,
            selectedContrast = input$heatmap_contrast,
            selectedAssay = input$heatmap_assay,
            topN = input$heatmap_top_n,
            clusteringMethod = input$heatmap_cluster_method,
            distanceMethod = input$heatmap_distance_method,
            heatmapClustering = input$heatmap_clustering,
            scaleData = input$heatmap_scaling,
            colorScheme = input$heatmap_color_scheme,
            showMetaboliteNames = input$heatmap_show_labels,
            daQValThresh = input$da_q_val_thresh,
            treeCutMethod = input$heatmap_tree_cut_method,
            nClusters = input$heatmap_n_clusters,
            cutHeight = input$heatmap_cut_height,
            minClusterSize = input$heatmap_min_cluster_size,
            daData = daData
        )
    })

    output$cluster_summary <- renderPrint({
        buildClusterSummaryOutput(
            treeCutMethod = input$heatmap_tree_cut_method,
            clusters = daData$current_row_clusters
        )
    })

    invisible(output)
}

registerMetabDaResultsOutputs <- function(
    input,
    output,
    daData,
    renderPrint = shiny::renderPrint,
    renderDt = DT::renderDT,
    buildSummaryStatsOutput = buildMetabDaSummaryStatsRenderOutput,
    buildResultsTableOutput = buildMetabDaResultsTableRenderOutput,
    buildResultsDownloadHandler = buildMetabDaResultsDownloadHandler
) {
    output$da_summary_stats <- renderPrint({
        buildSummaryStatsOutput(
            daResultsList = daData$da_results_list,
            selectedContrast = input$table_contrast,
            selectedAssay = input$table_assay,
            daQValThresh = input$da_q_val_thresh
        )
    })

    output$da_results_table <- renderDt({
        buildResultsTableOutput(
            daResultsList = daData$da_results_list,
            selectedContrast = input$table_contrast,
            selectedAssay = input$table_assay,
            significanceFilter = input$table_significance,
            daQValThresh = input$da_q_val_thresh,
            treatLfcCutoff = input$treat_lfc_cutoff,
            maxRows = input$table_max_rows
        )
    })

    output$download_da_results <- buildResultsDownloadHandler(
        daData = daData
    )

    invisible(output)
}

registerMetabDaSaveHeatmapObserver <- function(
    input,
    daData,
    experimentPaths,
    observeEvent = shiny::observeEvent,
    runSaveHeatmapObserverEntry = runMetabDaSaveHeatmapObserverEntry
) {
    observer <- observeEvent(input$save_heatmap, {
        runSaveHeatmapObserverEntry(
            currentHeatmapPlot = daData$current_heatmap_plot,
            publicationGraphsDir = experimentPaths$publication_graphs_dir,
            currentRowClusters = daData$current_row_clusters,
            selectedContrast = input$heatmap_contrast,
            topN = input$heatmap_top_n,
            clusteringMethod = input$heatmap_cluster_method,
            distanceMethod = input$heatmap_distance_method,
            heatmapClustering = input$heatmap_clustering,
            scaling = input$heatmap_scaling,
            colorScheme = input$heatmap_color_scheme,
            treeCutMethod = input$heatmap_tree_cut_method,
            nClusters = input$heatmap_n_clusters,
            cutHeight = input$heatmap_cut_height,
            minClusterSize = input$heatmap_min_cluster_size
        )
    })

    invisible(observer)
}

registerMetabDaLoadFilteredSessionObserver <- function(
    input,
    workflowData,
    daData,
    session,
    experimentPaths,
    observeEvent = shiny::observeEvent,
    runLoadSessionObserverEntry = runMetabDaLoadSessionObserverEntry,
    logInfo = logger::log_info,
    messageFn = message,
    sysTime = Sys.time
) {
    observer <- observeEvent(input$load_filtered_session, {
        d66_log <- function(...) messageFn(sprintf("[D66] %s", paste0(...)))
        d66_log("--- ENTER load_filtered_session observer ---")
        d66_start <- sysTime()

        logInfo("=== LOAD FILTERED SESSION BUTTON CLICKED ===")

        runLoadSessionObserverEntry(
            experimentPaths = experimentPaths,
            workflowData = workflowData,
            daData = daData,
            session = session,
            debugLog = d66_log,
            startTime = d66_start
        )
    })

    invisible(observer)
}

registerMetabDaRunAnalysisObserver <- function(
    input,
    workflowData,
    daData,
    session,
    experimentPaths,
    observeEvent = shiny::observeEvent,
    runAnalysisObserverEntry = runMetabDaAnalysisObserverEntry,
    logInfo = logger::log_info
) {
    observer <- observeEvent(input$run_da_analysis, {
        logInfo("=== RUN DA ANALYSIS BUTTON CLICKED ===")

        runAnalysisObserverEntry(
            currentS4Object = daData$current_s4_object,
            workflowData = workflowData,
            formulaString = input$formula_string,
            daQValThresh = input$da_q_val_thresh,
            treatLfcCutoff = input$treat_lfc_cutoff,
            daData = daData,
            session = session,
            experimentPaths = experimentPaths
        )
    })

    invisible(observer)
}

registerMetabDaMainOutputsAndObservers <- function(
    input,
    output,
    workflowData,
    daData,
    session,
    experimentPaths,
    registerLoadFilteredSessionObserver = registerMetabDaLoadFilteredSessionObserver,
    registerRunAnalysisObserver = registerMetabDaRunAnalysisObserver,
    registerVisualizationOutputs = registerMetabDaVisualizationOutputs,
    registerSaveHeatmapObserver = registerMetabDaSaveHeatmapObserver,
    registerResultsOutputs = registerMetabDaResultsOutputs
) {
    registrations <- list(
        loadFilteredSessionObserver = registerLoadFilteredSessionObserver(
            input = input,
            workflowData = workflowData,
            daData = daData,
            session = session,
            experimentPaths = experimentPaths
        ),
        runAnalysisObserver = registerRunAnalysisObserver(
            input = input,
            workflowData = workflowData,
            daData = daData,
            session = session,
            experimentPaths = experimentPaths
        ),
        visualizationOutputs = registerVisualizationOutputs(
            input = input,
            output = output,
            daData = daData
        ),
        saveHeatmapObserver = registerSaveHeatmapObserver(
            input = input,
            daData = daData,
            experimentPaths = experimentPaths
        ),
        resultsOutputs = registerResultsOutputs(
            input = input,
            output = output,
            daData = daData
        )
    )

    invisible(registrations)
}

