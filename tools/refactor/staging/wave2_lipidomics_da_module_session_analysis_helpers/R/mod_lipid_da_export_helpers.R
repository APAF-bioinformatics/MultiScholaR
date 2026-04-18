bootstrapLipidDaSaveHeatmap <- function(
    currentHeatmapPlot,
    currentRowClusters,
    experimentPaths,
    heatmapContrast,
    heatmapTopN,
    heatmapClusterMethod,
    heatmapDistanceMethod,
    heatmapClustering,
    heatmapScaling,
    heatmapColorScheme,
    heatmapTreeCutMethod,
    heatmapNClusters,
    heatmapCutHeight,
    heatmapMinClusterSize,
    handleSaveHeatmapFn = handleLipidDaSaveHeatmap
) {
    handleSaveHeatmapFn(
        currentHeatmapPlot = currentHeatmapPlot,
        currentRowClusters = currentRowClusters,
        publicationGraphsDir = experimentPaths$publication_graphs_dir,
        heatmapContrast = heatmapContrast,
        heatmapTopN = heatmapTopN,
        heatmapClusterMethod = heatmapClusterMethod,
        heatmapDistanceMethod = heatmapDistanceMethod,
        heatmapClustering = heatmapClustering,
        heatmapScaling = heatmapScaling,
        heatmapColorScheme = heatmapColorScheme,
        heatmapTreeCutMethod = heatmapTreeCutMethod,
        heatmapNClusters = heatmapNClusters,
        heatmapCutHeight = heatmapCutHeight,
        heatmapMinClusterSize = heatmapMinClusterSize
    )
}

handleLipidDaSaveHeatmap <- function(
    currentHeatmapPlot,
    currentRowClusters,
    publicationGraphsDir,
    heatmapContrast,
    heatmapTopN,
    heatmapClusterMethod,
    heatmapDistanceMethod,
    heatmapClustering,
    heatmapScaling,
    heatmapColorScheme,
    heatmapTreeCutMethod,
    heatmapNClusters,
    heatmapCutHeight,
    heatmapMinClusterSize,
    saveHeatmapProducts = save_heatmap_products,
    logInfo = logger::log_info,
    notify = shiny::showNotification
) {
    shiny::req(currentHeatmapPlot, publicationGraphsDir)

    logInfo("Save Heatmap button clicked")

    params <- list(
        contrast = heatmapContrast,
        top_n = heatmapTopN,
        clustering_method = heatmapClusterMethod,
        distance_method = heatmapDistanceMethod,
        cluster_rows = heatmapClustering %in% c("both", "row"),
        cluster_cols = heatmapClustering %in% c("both", "column"),
        scaling = heatmapScaling,
        color_scheme = heatmapColorScheme,
        tree_cut_method = heatmapTreeCutMethod,
        n_clusters = heatmapNClusters,
        cut_height = heatmapCutHeight,
        min_cluster_size = heatmapMinClusterSize
    )

    prefix <- gsub("[^A-Za-z0-9_]", "_", paste0("lipid_", heatmapContrast))

    saveHeatmapProducts(
        heatmap_obj = currentHeatmapPlot,
        row_clusters = currentRowClusters,
        params_list = params,
        output_dir = publicationGraphsDir,
        file_prefix = prefix
    )

    notify(
        "Heatmap and cluster info saved to publication_graphs/Heatmap",
        type = "message",
        duration = 5
    )

    invisible(list(params = params, prefix = prefix))
}

