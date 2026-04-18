library(testthat)

source(test_path("..", "..", "R", "mod_lipid_da_display_helpers.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_da_results_helpers.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_da_session_helpers.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_da_analysis_helpers.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_da_export_helpers.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_da_ui.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_da_server.R"), local = environment())
source(test_path("..", "..", "R", "mod_lipid_da.R"), local = environment())

if (!methods::isClass("MockLipidDaArgsCarrier")) {
    methods::setClass("MockLipidDaArgsCarrier", slots = c(args = "list"))
}

newMockLipidDaArgsCarrier <- function(args = list()) {
    methods::new("MockLipidDaArgsCarrier", args = args)
}

test_that("initializeLipidDaServerState keeps the reactive-value defaults stable", {
    factory_call <- NULL

    state <- initializeLipidDaServerState(
        reactiveValuesFactory = function(...) {
            factory_call <<- list(...)
            structure(list(args = factory_call), class = "mock_lipid_da_state")
        }
    )

    expect_s3_class(state, "mock_lipid_da_state")
    expect_identical(names(factory_call), c(
        "da_results_list",
        "contrasts_available",
        "assays_available",
        "analysis_complete",
        "current_s4_object",
        "formula_from_s4",
        "current_row_clusters",
        "current_col_clusters",
        "current_heatmap_plot"
    ))
    expect_null(factory_call$da_results_list)
    expect_null(factory_call$contrasts_available)
    expect_null(factory_call$assays_available)
    expect_false(factory_call$analysis_complete)
    expect_null(factory_call$current_s4_object)
    expect_null(factory_call$formula_from_s4)
    expect_null(factory_call$current_row_clusters)
    expect_null(factory_call$current_col_clusters)
    expect_null(factory_call$current_heatmap_plot)
})

test_that("registerLipidDaPrimaryTextOutputs keeps the paired render registrations stable", {
    output <- new.env(parent = emptyenv())
    workflow_data <- new.env(parent = emptyenv())
    da_data <- new.env(parent = emptyenv())
    render_calls <- character()
    captured <- list()

    workflow_data$contrasts_tbl <- data.frame(
        friendly_names = c("B vs A", "C vs A"),
        stringsAsFactors = FALSE
    )
    da_data$analysis_complete <- TRUE
    da_data$da_results_list <- list(marker = "results")
    da_data$current_s4_object <- structure(list(marker = "s4"), class = "mock_s4")

    result <- registerLipidDaPrimaryTextOutputs(
        output = output,
        workflowData = workflow_data,
        daData = da_data,
        renderPrintFn = function(expr) {
            render_calls <<- c(render_calls, force(expr))
            structure(list(kind = "mock_render", index = length(render_calls)), class = "mock_render")
        },
        buildContrastsTextFn = function(contrastsTbl) {
            captured$contrastsTbl <<- contrastsTbl
            "contrasts-text"
        },
        buildStatusTextFn = function(analysisComplete, daResultsList, currentS4Object) {
            captured$statusArgs <<- list(
                analysisComplete = analysisComplete,
                daResultsList = daResultsList,
                currentS4Object = currentS4Object
            )
            "status-text"
        },
        catFn = function(...) paste0(..., collapse = "")
    )

    expect_identical(result, output)
    expect_identical(render_calls, c("contrasts-text", "status-text"))
    expect_identical(captured$contrastsTbl, workflow_data$contrasts_tbl)
    expect_true(captured$statusArgs$analysisComplete)
    expect_identical(captured$statusArgs$daResultsList, da_data$da_results_list)
    expect_identical(captured$statusArgs$currentS4Object, da_data$current_s4_object)
    expect_s3_class(output$contrasts_display, "mock_render")
    expect_s3_class(output$da_status, "mock_render")
})

test_that("registerLipidDaHeatmapWarningOutput keeps the warning render registration stable", {
    output <- new.env(parent = emptyenv())
    da_data <- new.env(parent = emptyenv())
    render_calls <- character()
    captured <- list()

    da_data$analysis_complete <- TRUE

    result <- registerLipidDaHeatmapWarningOutput(
        output = output,
        daData = da_data,
        renderUiFn = function(expr) {
            render_calls <<- c(render_calls, force(expr))
            structure(list(kind = "mock_render_ui"), class = "mock_render_ui")
        },
        buildWarningFn = function(analysisComplete) {
            captured$analysisComplete <<- analysisComplete
            "warning-ui"
        }
    )

    expect_identical(result, output)
    expect_identical(render_calls, "warning-ui")
    expect_true(captured$analysisComplete)
    expect_s3_class(output$heatmap_manual_save_warning, "mock_render_ui")
})

test_that("registerLipidDaVolcanoGlimmaOutput keeps the render registration handoff stable", {
    output <- new.env(parent = emptyenv())
    input <- list(
        volcano_contrast = "B_vs_A",
        volcano_assay = "Assay1",
        da_q_val_thresh = 0.1
    )
    da_data <- new.env(parent = emptyenv())
    render_calls <- character()
    captured <- list()

    da_data$da_results_list <- list(marker = TRUE)

    result <- registerLipidDaVolcanoGlimmaOutput(
        output = output,
        input = input,
        daData = da_data,
        renderUiFn = function(expr) {
            render_calls <<- c(render_calls, force(expr))
            structure(list(kind = "mock_render_ui"), class = "mock_render_ui")
        },
        buildVolcanoGlimmaUiFn = function(
            daResultsList,
            selectedContrast,
            selectedAssay,
            daQValThresh
        ) {
            captured$daResultsList <<- daResultsList
            captured$selectedContrast <<- selectedContrast
            captured$selectedAssay <<- selectedAssay
            captured$daQValThresh <<- daQValThresh
            "volcano-ui"
        }
    )

    expect_identical(result, output)
    expect_identical(render_calls, "volcano-ui")
    expect_identical(captured$daResultsList, da_data$da_results_list)
    expect_identical(captured$selectedContrast, "B_vs_A")
    expect_identical(captured$selectedAssay, "Assay1")
    expect_identical(captured$daQValThresh, 0.1)
    expect_s3_class(output$volcano_glimma, "mock_render_ui")
})

test_that("registerLipidDaVolcanoStaticOutput keeps the render registration handoff stable", {
    output <- new.env(parent = emptyenv())
    input <- list(
        volcano_contrast = "B_vs_A",
        volcano_assay = "Assay1",
        da_q_val_thresh = 0.1,
        treat_lfc_cutoff = 0.5
    )
    da_data <- new.env(parent = emptyenv())
    render_calls <- character()
    captured <- list()

    da_data$da_results_list <- list(marker = TRUE)

    result <- registerLipidDaVolcanoStaticOutput(
        output = output,
        input = input,
        daData = da_data,
        renderPlotFn = function(expr) {
            render_calls <<- c(render_calls, force(expr))
            structure(list(kind = "mock_render_plot"), class = "mock_render_plot")
        },
        buildVolcanoStaticPlotFn = function(
            daResultsList,
            selectedContrast,
            selectedAssay,
            daQValThresh,
            lfcThreshold
        ) {
            captured$daResultsList <<- daResultsList
            captured$selectedContrast <<- selectedContrast
            captured$selectedAssay <<- selectedAssay
            captured$daQValThresh <<- daQValThresh
            captured$lfcThreshold <<- lfcThreshold
            "volcano-static-plot"
        }
    )

    expect_identical(result, output)
    expect_identical(render_calls, "volcano-static-plot")
    expect_identical(captured$daResultsList, da_data$da_results_list)
    expect_identical(captured$selectedContrast, "B_vs_A")
    expect_identical(captured$selectedAssay, "Assay1")
    expect_identical(captured$daQValThresh, 0.1)
    expect_identical(captured$lfcThreshold, 0.5)
    expect_s3_class(output$volcano_static, "mock_render_plot")
})

test_that("registerLipidDaHeatmapPlotOutput keeps the render registration handoff stable", {
    output <- new.env(parent = emptyenv())
    input <- list(
        heatmap_contrast = "B_vs_A",
        heatmap_assay = "Assay1",
        heatmap_top_n = 75,
        heatmap_clustering = "column",
        heatmap_scaling = "both",
        heatmap_cluster_method = "complete",
        heatmap_distance_method = "pearson",
        heatmap_color_scheme = "plasma",
        heatmap_show_labels = TRUE,
        da_q_val_thresh = 0.1,
        heatmap_tree_cut_method = "dynamic",
        heatmap_n_clusters = 6,
        heatmap_cut_height = 12,
        heatmap_min_cluster_size = 7
    )
    da_data <- new.env(parent = emptyenv())
    render_values <- list()
    captured <- list()
    stored_plot <- structure(list(kind = "stored-plot"), class = "mock_plot")

    da_data$da_results_list <- list(marker = TRUE)

    result <- registerLipidDaHeatmapPlotOutput(
        output = output,
        input = input,
        daData = da_data,
        renderPlotFn = function(expr) {
            render_values$plot <<- force(expr)
            structure(list(kind = "mock_render_plot"), class = "mock_render_plot")
        },
        buildHeatmapRenderStateFn = function(
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
            minClusterSize
        ) {
            captured$buildArgs <<- list(
                daResultsList = daResultsList,
                selectedContrast = selectedContrast,
                selectedAssay = selectedAssay,
                topN = topN,
                heatmapClustering = heatmapClustering,
                scaleData = scaleData,
                clusteringMethod = clusteringMethod,
                distanceMethod = distanceMethod,
                colorScheme = colorScheme,
                showLipidNames = showLipidNames,
                daQValThresh = daQValThresh,
                treeCutMethod = treeCutMethod,
                nClusters = nClusters,
                cutHeight = cutHeight,
                minClusterSize = minClusterSize
            )
            list(
                plot = structure(list(kind = "rendered-plot"), class = "mock_plot"),
                rowClusters = c(LipidA = 1L, LipidB = 2L),
                colClusters = structure(list(kind = "col-clusters"), class = "mock_clusters"),
                currentHeatmapPlot = stored_plot
            )
        },
        storeHeatmapRenderStateFn = function(heatmapState, daData) {
            captured$storeArgs <<- list(
                heatmapState = heatmapState,
                daData = daData
            )
            stored_plot
        }
    )

    expect_identical(result, output)
    expect_identical(render_values$plot, stored_plot)
    expect_identical(captured$buildArgs$daResultsList, da_data$da_results_list)
    expect_identical(captured$buildArgs$selectedContrast, "B_vs_A")
    expect_identical(captured$buildArgs$selectedAssay, "Assay1")
    expect_identical(captured$buildArgs$topN, 75)
    expect_identical(captured$buildArgs$heatmapClustering, "column")
    expect_identical(captured$buildArgs$scaleData, "both")
    expect_identical(captured$buildArgs$clusteringMethod, "complete")
    expect_identical(captured$buildArgs$distanceMethod, "pearson")
    expect_identical(captured$buildArgs$colorScheme, "plasma")
    expect_true(captured$buildArgs$showLipidNames)
    expect_identical(captured$buildArgs$daQValThresh, 0.1)
    expect_identical(captured$buildArgs$treeCutMethod, "dynamic")
    expect_identical(captured$buildArgs$nClusters, 6)
    expect_identical(captured$buildArgs$cutHeight, 12)
    expect_identical(captured$buildArgs$minClusterSize, 7)
    expect_identical(captured$storeArgs$heatmapState$currentHeatmapPlot, stored_plot)
    expect_identical(captured$storeArgs$daData, da_data)
    expect_s3_class(output$heatmap_plot, "mock_render_plot")
})

test_that("registerLipidDaClusterSummaryOutput keeps the render registration handoff stable", {
    output <- new.env(parent = emptyenv())
    input <- list(heatmap_tree_cut_method = "dynamic")
    da_data <- new.env(parent = emptyenv())
    render_calls <- character()
    captured <- list()

    da_data$current_row_clusters <- c(LipidA = 1L, LipidB = 2L)

    result <- registerLipidDaClusterSummaryOutput(
        output = output,
        input = input,
        daData = da_data,
        renderPrintFn = function(expr) {
            render_calls <<- c(render_calls, force(expr))
            structure(list(kind = "mock_render"), class = "mock_render")
        },
        buildClusterSummaryRenderTextFn = function(treeCutMethod, rowClusters) {
            captured$treeCutMethod <<- treeCutMethod
            captured$rowClusters <<- rowClusters
            "cluster-summary-text"
        },
        catFn = function(...) paste0(..., collapse = "")
    )

    expect_identical(result, output)
    expect_identical(render_calls, "cluster-summary-text")
    expect_identical(captured$treeCutMethod, "dynamic")
    expect_identical(captured$rowClusters, da_data$current_row_clusters)
    expect_s3_class(output$cluster_summary, "mock_render")
})

test_that("registerLipidDaSaveHeatmapObserver keeps the observer shell handoff stable", {
    input <- list(
        save_heatmap = 4L,
        heatmap_contrast = "B_vs_A",
        heatmap_top_n = 75,
        heatmap_cluster_method = "complete",
        heatmap_distance_method = "pearson",
        heatmap_clustering = "column",
        heatmap_scaling = "both",
        heatmap_color_scheme = "plasma",
        heatmap_tree_cut_method = "dynamic",
        heatmap_n_clusters = 6,
        heatmap_cut_height = 12,
        heatmap_min_cluster_size = 7
    )
    da_data <- new.env(parent = emptyenv())
    observer_event <- NULL
    bootstrap_args <- NULL
    heatmap_plot <- structure(list(kind = "heatmap"), class = "mock_heatmap_plot")
    row_clusters <- c(LipidA = 1L, LipidB = 2L)

    da_data$current_heatmap_plot <- heatmap_plot
    da_data$current_row_clusters <- row_clusters

    registerLipidDaSaveHeatmapObserver(
        input = input,
        daData = da_data,
        experimentPaths = list(publication_graphs_dir = "/tmp/publication_graphs"),
        bootstrapSaveHeatmapFn = function(...) {
            bootstrap_args <<- list(...)
            "saved"
        },
        observeEventFn = function(event, handler) {
            observer_event <<- event
            force(handler)
            "registered"
        }
    )

    expect_identical(observer_event, 4L)
    expect_identical(bootstrap_args$currentHeatmapPlot, heatmap_plot)
    expect_identical(bootstrap_args$currentRowClusters, row_clusters)
    expect_identical(
        bootstrap_args$experimentPaths$publication_graphs_dir,
        "/tmp/publication_graphs"
    )
    expect_identical(bootstrap_args$heatmapContrast, "B_vs_A")
    expect_identical(bootstrap_args$heatmapTopN, 75)
    expect_identical(bootstrap_args$heatmapClusterMethod, "complete")
    expect_identical(bootstrap_args$heatmapDistanceMethod, "pearson")
    expect_identical(bootstrap_args$heatmapClustering, "column")
    expect_identical(bootstrap_args$heatmapScaling, "both")
    expect_identical(bootstrap_args$heatmapColorScheme, "plasma")
    expect_identical(bootstrap_args$heatmapTreeCutMethod, "dynamic")
    expect_identical(bootstrap_args$heatmapNClusters, 6)
    expect_identical(bootstrap_args$heatmapCutHeight, 12)
    expect_identical(bootstrap_args$heatmapMinClusterSize, 7)
})

test_that("registerLipidDaSummaryStatsOutput keeps the render registration handoff stable", {
    output <- new.env(parent = emptyenv())
    input <- list(
        table_contrast = "B vs A",
        table_assay = "Assay1",
        da_q_val_thresh = 0.05
    )
    da_data <- new.env(parent = emptyenv())
    render_calls <- character()
    captured <- NULL

    da_data$da_results_list <- list(
        da_lipids_long = data.frame(marker = "value", stringsAsFactors = FALSE)
    )

    result <- registerLipidDaSummaryStatsOutput(
        output = output,
        input = input,
        daData = da_data,
        renderPrintFn = function(expr) {
            render_calls <<- c(render_calls, force(expr))
            structure(list(kind = "mock_render"), class = "mock_render")
        },
        buildSummaryStatsRenderTextFn = function(...) {
            captured <<- list(...)
            "summary-stats-text"
        },
        catFn = function(...) paste0(..., collapse = "")
    )

    expect_identical(result, output)
    expect_identical(render_calls, "summary-stats-text")
    expect_identical(captured$daResultsList, da_data$da_results_list)
    expect_identical(captured$selectedContrast, "B vs A")
    expect_identical(captured$selectedAssay, "Assay1")
    expect_identical(captured$daQValThresh, 0.05)
    expect_s3_class(output$da_summary_stats, "mock_render")
})

test_that("registerLipidDaResultsTableOutput keeps the render registration handoff stable", {
    output <- new.env(parent = emptyenv())
    input <- list(
        table_contrast = "B vs A",
        table_assay = "Assay1",
        table_significance = "significant",
        da_q_val_thresh = 0.05,
        treat_lfc_cutoff = 0.5,
        table_max_rows = 250
    )
    da_data <- new.env(parent = emptyenv())
    render_values <- list()
    captured <- NULL

    da_data$da_results_list <- list(
        da_lipids_long = data.frame(marker = "value", stringsAsFactors = FALSE)
    )

    result <- registerLipidDaResultsTableOutput(
        output = output,
        input = input,
        daData = da_data,
        renderDtFn = function(expr) {
            render_values$widget <<- force(expr)
            structure(list(kind = "mock_render_dt"), class = "mock_render_dt")
        },
        buildResultsTableRenderWidgetFn = function(...) {
            captured <<- list(...)
            "results-table-widget"
        }
    )

    expect_identical(result, output)
    expect_identical(render_values$widget, "results-table-widget")
    expect_identical(captured$daResultsList, da_data$da_results_list)
    expect_identical(captured$selectedContrast, "B vs A")
    expect_identical(captured$selectedAssay, "Assay1")
    expect_identical(captured$tableSignificance, "significant")
    expect_identical(captured$daQValThresh, 0.05)
    expect_identical(captured$lfcThreshold, 0.5)
    expect_identical(captured$tableMaxRows, 250)
    expect_s3_class(output$da_results_table, "mock_render_dt")
})

test_that("registerLipidDaResultsDownloadOutput keeps the output registration handoff stable", {
    output <- new.env(parent = emptyenv())
    da_data <- new.env(parent = emptyenv())
    captured <- NULL

    da_data$da_results_list <- list(
        da_lipids_long = data.frame(marker = "value", stringsAsFactors = FALSE)
    )

    result <- registerLipidDaResultsDownloadOutput(
        output = output,
        daData = da_data,
        buildResultsDownloadOutputHandlerFn = function(...) {
            captured <<- list(...)
            structure(list(kind = "mock_download_handler"), class = "mock_download_handler")
        }
    )

    expect_identical(result, output)
    expect_identical(captured$daResultsList, da_data$da_results_list)
    expect_s3_class(output$download_da_results, "mock_download_handler")
})

test_that("buildLipidDaHeatmapManualSaveWarning stays hidden before analysis completes", {
    expect_null(buildLipidDaHeatmapManualSaveWarning(FALSE))
    expect_null(buildLipidDaHeatmapManualSaveWarning(NULL))
})

test_that("buildLipidDaHeatmapManualSaveWarning keeps the manual-save guidance stable", {
    rendered <- htmltools::renderTags(
        buildLipidDaHeatmapManualSaveWarning(TRUE)
    )$html

    expect_match(rendered, "alert alert-info", fixed = TRUE)
    expect_match(rendered, "Heatmaps are NOT saved automatically", fixed = TRUE)
    expect_match(rendered, "Save Heatmap", fixed = TRUE)
})

test_that("buildLipidDaVolcanoGlimmaUi keeps the combined-view guidance stable", {
    rendered <- htmltools::renderTags(
        buildLipidDaVolcanoGlimmaUi(
            daResultsList = list(marker = TRUE),
            selectedContrast = "B_vs_A",
            selectedAssay = "Combined",
            daQValThresh = 0.05,
            widgetFactory = function(...) stop("widget factory should not run"),
            logError = function(...) stop("logger should not run")
        )
    )$html

    expect_match(rendered, "alert alert-info", fixed = TRUE)
    expect_match(rendered, "Combined View:", fixed = TRUE)
    expect_match(rendered, "single assay selection", fixed = TRUE)
})

test_that("buildLipidDaVolcanoGlimmaUi delegates single-assay rendering through the lipid Glimma helper", {
    widget_calls <- list()
    widget <- shiny::div(class = "mock-widget", "ok")

    result <- buildLipidDaVolcanoGlimmaUi(
        daResultsList = list(marker = TRUE),
        selectedContrast = "B_vs_A",
        selectedAssay = "Assay1",
        daQValThresh = 0.1,
        widgetFactory = function(...) {
            widget_calls <<- list(...)
            widget
        },
        logError = function(...) stop("logger should not run")
    )

    expect_identical(result, widget)
    expect_identical(widget_calls$da_results_list$marker, TRUE)
    expect_identical(widget_calls$selected_contrast, "B_vs_A")
    expect_identical(widget_calls$selected_assay, "Assay1")
    expect_identical(widget_calls$da_q_val_thresh, 0.1)
})

test_that("buildLipidDaVolcanoGlimmaUi keeps warning and error fallbacks stable", {
    warning_html <- htmltools::renderTags(
        buildLipidDaVolcanoGlimmaUi(
            daResultsList = list(marker = TRUE),
            selectedContrast = "B_vs_A",
            selectedAssay = "Assay1",
            daQValThresh = 0.05,
            widgetFactory = function(...) NULL,
            logError = function(...) stop("logger should not run")
        )
    )$html

    logged_message <- NULL
    error_html <- htmltools::renderTags(
        buildLipidDaVolcanoGlimmaUi(
            daResultsList = list(marker = TRUE),
            selectedContrast = "B_vs_A",
            selectedAssay = "Assay1",
            daQValThresh = 0.05,
            widgetFactory = function(...) stop("boom"),
            logError = function(message) {
                logged_message <<- message
            }
        )
    )$html

    expect_match(warning_html, "alert alert-warning", fixed = TRUE)
    expect_match(warning_html, "Could not generate Glimma plot", fixed = TRUE)

    expect_match(error_html, "alert alert-danger", fixed = TRUE)
    expect_match(error_html, "Error generating plot: boom", fixed = TRUE)
    expect_identical(logged_message, "Glimma error: boom")
})

test_that("buildLipidDaVolcanoStaticPlot delegates static rendering through the lipid helper", {
    plot_calls <- list()
    plot_result <- structure(list(kind = "static-volcano"), class = "mock_plot")

    result <- buildLipidDaVolcanoStaticPlot(
        daResultsList = list(marker = TRUE),
        selectedContrast = "B_vs_A",
        selectedAssay = "Assay1",
        daQValThresh = 0.1,
        lfcThreshold = 0.5,
        plotFactory = function(...) {
            plot_calls <<- list(...)
            plot_result
        }
    )

    expect_identical(result, plot_result)
    expect_identical(plot_calls$da_results_list$marker, TRUE)
    expect_identical(plot_calls$selected_contrast, "B_vs_A")
    expect_identical(plot_calls$selected_assay, "Assay1")
    expect_identical(plot_calls$da_q_val_thresh, 0.1)
    expect_identical(plot_calls$lfc_threshold, 0.5)
    expect_identical(plot_calls$show_labels, TRUE)
    expect_identical(plot_calls$n_labels, 15)
})

test_that("buildLipidDaHeatmapRenderState delegates the heatmap render contract through the lipid helper", {
    heatmap_calls <- list()
    heatmap_plot <- structure(list(kind = "heatmap"), class = "mock_heatmap_plot")
    heatmap_result <- list(
        plot = heatmap_plot,
        row_clusters = c(LipidA = 1L, LipidB = 2L),
        col_clusters = structure(list(kind = "col-clusters"), class = "mock_clusters")
    )

    result <- buildLipidDaHeatmapRenderState(
        daResultsList = list(marker = TRUE),
        selectedContrast = "B_vs_A",
        selectedAssay = "Assay1",
        topN = 75,
        heatmapClustering = "column",
        scaleData = "both",
        clusteringMethod = "complete",
        distanceMethod = "pearson",
        colorScheme = "plasma",
        showLipidNames = TRUE,
        daQValThresh = 0.1,
        treeCutMethod = "dynamic",
        nClusters = 6,
        cutHeight = 12,
        minClusterSize = 7,
        heatmapFactory = function(...) {
            heatmap_calls <<- list(...)
            heatmap_result
        }
    )

    expect_identical(result$plot, heatmap_plot)
    expect_identical(result$rowClusters, heatmap_result$row_clusters)
    expect_identical(result$colClusters, heatmap_result$col_clusters)
    expect_identical(result$currentHeatmapPlot, heatmap_plot)

    expect_identical(heatmap_calls$da_results_list$marker, TRUE)
    expect_identical(heatmap_calls$selected_contrast, "B_vs_A")
    expect_identical(heatmap_calls$selected_assay, "Assay1")
    expect_identical(heatmap_calls$top_n, 75)
    expect_identical(heatmap_calls$clustering_method, "complete")
    expect_identical(heatmap_calls$distance_method, "pearson")
    expect_identical(heatmap_calls$cluster_rows, FALSE)
    expect_identical(heatmap_calls$cluster_cols, TRUE)
    expect_identical(heatmap_calls$scale_data, "both")
    expect_identical(heatmap_calls$color_scheme, "plasma")
    expect_identical(heatmap_calls$show_lipid_names, TRUE)
    expect_identical(heatmap_calls$da_q_val_thresh, 0.1)
    expect_identical(heatmap_calls$tree_cut_method, "dynamic")
    expect_identical(heatmap_calls$n_clusters, 6)
    expect_identical(heatmap_calls$cut_height, 12)
    expect_identical(heatmap_calls$min_cluster_size, 7)
})

test_that("storeLipidDaHeatmapRenderState keeps the NULL early return stable", {
    da_data <- new.env(parent = emptyenv())
    da_data$current_row_clusters <- c(LipidA = 1L)
    da_data$current_col_clusters <- "existing-col-clusters"
    da_data$current_heatmap_plot <- "existing-heatmap-plot"

    result <- storeLipidDaHeatmapRenderState(
        heatmapState = NULL,
        daData = da_data
    )

    expect_null(result)
    expect_identical(da_data$current_row_clusters, c(LipidA = 1L))
    expect_identical(da_data$current_col_clusters, "existing-col-clusters")
    expect_identical(da_data$current_heatmap_plot, "existing-heatmap-plot")
})

test_that("storeLipidDaHeatmapRenderState keeps the state-persistence handoff stable", {
    da_data <- new.env(parent = emptyenv())
    rendered_plot <- structure(list(kind = "rendered-plot"), class = "mock_plot")
    stored_plot <- structure(list(kind = "stored-plot"), class = "mock_plot")
    row_clusters <- c(LipidA = 1L, LipidB = 2L)
    col_clusters <- structure(list(kind = "col-clusters"), class = "mock_clusters")

    result <- storeLipidDaHeatmapRenderState(
        heatmapState = list(
            plot = rendered_plot,
            rowClusters = row_clusters,
            colClusters = col_clusters,
            currentHeatmapPlot = stored_plot
        ),
        daData = da_data
    )

    expect_identical(result, rendered_plot)
    expect_identical(da_data$current_row_clusters, row_clusters)
    expect_identical(da_data$current_col_clusters, col_clusters)
    expect_identical(da_data$current_heatmap_plot, stored_plot)
})

test_that("buildLipidDaClusterSummaryText keeps the empty-cluster guidance stable", {
    expect_identical(
        buildLipidDaClusterSummaryText(NULL),
        "No clusters defined. Enable clustering and tree cutting on the heatmap."
    )
})

test_that("buildLipidDaClusterSummaryText keeps the cluster summary formatting stable", {
    cluster_members <- c(sprintf("Lipid%02d", seq_len(22)), "LipidZ")
    cluster_assignments <- stats::setNames(c(rep(1L, 22), 2L), cluster_members)

    summary_text <- buildLipidDaClusterSummaryText(cluster_assignments)

    expect_match(summary_text, "Total Clusters: 2", fixed = TRUE)
    expect_match(summary_text, "Cluster 1 \\(22 lipids\\):")
    expect_match(summary_text, "Lipid01, Lipid02, Lipid03", fixed = TRUE)
    expect_match(summary_text, "\\.\\.\\. and 2 more")
    expect_match(summary_text, "Cluster 2 \\(1 lipids\\):")
    expect_match(summary_text, "LipidZ", fixed = TRUE)
})

test_that("buildLipidDaClusterSummaryRenderText keeps the render-shell handoff stable", {
    req_arg <- NULL
    build_arg <- NULL
    row_clusters <- c(LipidA = 1L, LipidB = 2L)

    result <- buildLipidDaClusterSummaryRenderText(
        treeCutMethod = "dynamic",
        rowClusters = row_clusters,
        reqFn = function(value) {
            req_arg <<- value
            value
        },
        buildClusterSummaryTextFn = function(value) {
            build_arg <<- value
            "cluster-summary-text"
        }
    )

    expect_identical(req_arg, TRUE)
    expect_identical(build_arg, row_clusters)
    expect_identical(result, "cluster-summary-text")
})

test_that("buildLipidDaStatusText keeps the waiting, ready, and plain-complete guidance stable", {
    expect_identical(
        buildLipidDaStatusText(
            analysisComplete = FALSE,
            daResultsList = NULL,
            currentS4Object = NULL
        ),
        "Waiting for data.\nClick 'Load Filtered Session' to begin."
    )

    expect_identical(
        buildLipidDaStatusText(
            analysisComplete = FALSE,
            daResultsList = NULL,
            currentS4Object = structure(list(marker = TRUE), class = "mock_s4")
        ),
        "Data loaded. Ready to run analysis."
    )

    expect_identical(
        buildLipidDaStatusText(
            analysisComplete = TRUE,
            daResultsList = list(),
            currentS4Object = NULL
        ),
        "Analysis complete."
    )
})

test_that("buildLipidDaStatusText keeps the significant-count status summary stable", {
    status_text <- buildLipidDaStatusText(
        analysisComplete = TRUE,
        daResultsList = list(
            significant_counts = list(
                LCMS_Pos = list(up = 3L, down = 1L, ns = 9L),
                LCMS_Neg = list(up = 0L, down = 2L, ns = 5L)
            )
        ),
        currentS4Object = NULL
    )

    expect_identical(
        status_text,
        paste(
            c(
                "Analysis Complete",
                "",
                "LCMS_Pos:\n  Up: 3 | Down: 1 | NS: 9",
                "LCMS_Neg:\n  Up: 0 | Down: 2 | NS: 5"
            ),
            collapse = "\n"
        )
    )
})

test_that("buildLipidDaContrastsText keeps the empty-contrasts guidance stable", {
    expect_identical(
        buildLipidDaContrastsText(NULL),
        "No contrasts defined.\nLoad a filtered session or define contrasts in the design tab."
    )

    expect_identical(
        buildLipidDaContrastsText(data.frame()),
        "No contrasts defined.\nLoad a filtered session or define contrasts in the design tab."
    )
})

test_that("buildLipidDaContrastsText prefers friendly contrast names when available", {
    expect_identical(
        buildLipidDaContrastsText(data.frame(
            friendly_names = c("B vs A", "C vs A"),
            contrasts = c("B_vs_A", "C_vs_A"),
            stringsAsFactors = FALSE
        )),
        "B vs A\nC vs A"
    )
})

test_that("buildLipidDaContrastsText falls back to raw contrast names when needed", {
    expect_identical(
        buildLipidDaContrastsText(data.frame(
            contrasts = c("B_vs_A", "C_vs_A"),
            stringsAsFactors = FALSE
        )),
        "B_vs_A\nC_vs_A"
    )
})

test_that("buildLipidDaContrastsText keeps the generic table fallback stable", {
    formatted_output <- buildLipidDaContrastsText(
        data.frame(other = "value", stringsAsFactors = FALSE),
        fallbackFormatter = function(x) {
            expect_identical(x$other, "value")
            c("  other", "1 value")
        }
    )

    expect_identical(formatted_output, "  other\n1 value")
})

test_that("restoreLipidDaContrastsFromSession keeps the empty-contrast contract stable", {
    workflow_data <- new.env(parent = emptyenv())
    da_data <- new.env(parent = emptyenv())
    update_calls <- list()
    logged_message <- NULL

    result <- restoreLipidDaContrastsFromSession(
        sessionData = list(contrasts_tbl = data.frame()),
        workflowData = workflow_data,
        daData = da_data,
        session = structure(list(kind = "mock-session"), class = "mock_session"),
        updateSelectInputFn = function(...) {
            update_calls[[length(update_calls) + 1]] <<- list(...)
        },
        logInfo = function(message) {
            logged_message <<- message
        }
    )

    expect_null(result)
    expect_length(update_calls, 0)
    expect_false(exists("contrasts_tbl", envir = workflow_data, inherits = FALSE))
    expect_false(exists("contrasts_available", envir = da_data, inherits = FALSE))
    expect_null(logged_message)
})

test_that("restoreLipidDaContrastsFromSession keeps the contrast-restore dropdown contract stable", {
    workflow_data <- new.env(parent = emptyenv())
    da_data <- new.env(parent = emptyenv())
    update_calls <- list()
    logged_message <- NULL
    mock_session <- structure(list(kind = "mock-session"), class = "mock_session")
    contrasts_tbl <- data.frame(
        friendly_names = c("B vs A", "C vs A"),
        contrasts = c("B_vs_A", "C_vs_A"),
        stringsAsFactors = FALSE
    )

    result <- restoreLipidDaContrastsFromSession(
        sessionData = list(contrasts_tbl = contrasts_tbl),
        workflowData = workflow_data,
        daData = da_data,
        session = mock_session,
        updateSelectInputFn = function(session, inputId, choices, selected) {
            update_calls[[length(update_calls) + 1]] <<- list(
                session = session,
                inputId = inputId,
                choices = choices,
                selected = selected
            )
        },
        logInfo = function(message) {
            logged_message <<- message
        }
    )

    expect_identical(result$contrastsTbl, contrasts_tbl)
    expect_identical(result$contrastChoices, c("B vs A", "C vs A"))
    expect_identical(workflow_data$contrasts_tbl, contrasts_tbl)
    expect_identical(da_data$contrasts_available, contrasts_tbl)
    expect_length(update_calls, 3)
    expect_identical(vapply(update_calls, `[[`, character(1), "inputId"), c(
        "volcano_contrast",
        "heatmap_contrast",
        "table_contrast"
    ))
    expect_identical(update_calls[[1]]$session, mock_session)
    expect_identical(update_calls[[1]]$choices, c("B vs A", "C vs A"))
    expect_identical(update_calls[[1]]$selected, "B vs A")
    expect_identical(update_calls[[2]]$choices, c("B vs A", "C vs A"))
    expect_identical(update_calls[[2]]$selected, "B vs A")
    expect_identical(update_calls[[3]]$choices, c("B vs A", "C vs A"))
    expect_identical(update_calls[[3]]$selected, "B vs A")
    expect_identical(logged_message, "   Restored 2 contrasts")
})

test_that("restoreLipidDaContrastsFromSession falls back to raw contrast names when needed", {
    workflow_data <- new.env(parent = emptyenv())
    da_data <- new.env(parent = emptyenv())
    update_calls <- list()

    result <- restoreLipidDaContrastsFromSession(
        sessionData = list(contrasts_tbl = data.frame(
            contrasts = c("B_vs_A", "C_vs_A"),
            stringsAsFactors = FALSE
        )),
        workflowData = workflow_data,
        daData = da_data,
        session = structure(list(kind = "mock-session"), class = "mock_session"),
        updateSelectInputFn = function(session, inputId, choices, selected) {
            update_calls[[length(update_calls) + 1]] <<- list(
                inputId = inputId,
                choices = choices,
                selected = selected
            )
        },
        logInfo = function(...) NULL
    )

    expect_identical(result$contrastChoices, c("B_vs_A", "C_vs_A"))
    expect_identical(update_calls[[1]]$choices, c("B_vs_A", "C_vs_A"))
    expect_identical(update_calls[[1]]$selected, "B_vs_A")
    expect_identical(da_data$contrasts_available$contrasts, c("B_vs_A", "C_vs_A"))
})

test_that("restoreLipidDaAssaysFromSession keeps the empty-assay contract stable", {
    da_data <- new.env(parent = emptyenv())
    update_calls <- list()
    debug_messages <- character()

    result <- restoreLipidDaAssaysFromSession(
        sessionData = list(assay_names = NULL),
        daData = da_data,
        session = structure(list(kind = "mock-session"), class = "mock_session"),
        updateSelectInputFn = function(...) {
            update_calls[[length(update_calls) + 1]] <<- list(...)
        },
        debugLog = function(message) {
            debug_messages <<- c(debug_messages, message)
        }
    )

    expect_null(result)
    expect_length(update_calls, 0)
    expect_false(exists("assays_available", envir = da_data, inherits = FALSE))
    expect_length(debug_messages, 0)
})

test_that("restoreLipidDaAssaysFromSession keeps the assay-restore dropdown contract stable", {
    da_data <- new.env(parent = emptyenv())
    update_calls <- list()
    debug_messages <- character()
    mock_session <- structure(list(kind = "mock-session"), class = "mock_session")

    result <- restoreLipidDaAssaysFromSession(
        sessionData = list(assay_names = c("Assay1", "Assay2")),
        daData = da_data,
        session = mock_session,
        updateSelectInputFn = function(session, inputId, choices) {
            update_calls[[length(update_calls) + 1]] <<- list(
                session = session,
                inputId = inputId,
                choices = choices
            )
        },
        debugLog = function(message) {
            debug_messages <<- c(debug_messages, message)
        }
    )

    expect_identical(result$assayNames, c("Assay1", "Assay2"))
    expect_identical(result$assayChoices, c("Combined", "Assay1", "Assay2"))
    expect_identical(result$tableAssayChoices, c("All", "Assay1", "Assay2"))
    expect_identical(da_data$assays_available, c("Assay1", "Assay2"))
    expect_length(update_calls, 3)
    expect_identical(vapply(update_calls, `[[`, character(1), "inputId"), c(
        "volcano_assay",
        "heatmap_assay",
        "table_assay"
    ))
    expect_identical(update_calls[[1]]$session, mock_session)
    expect_identical(update_calls[[1]]$choices, c("Combined", "Assay1", "Assay2"))
    expect_identical(update_calls[[2]]$choices, c("Combined", "Assay1", "Assay2"))
    expect_identical(update_calls[[3]]$choices, c("All", "Assay1", "Assay2"))
    expect_identical(debug_messages, c(
        "  STEP 8a: Updating volcano_assay dropdown...",
        "  STEP 8b: Updating heatmap_assay dropdown...",
        "  STEP 8c: Updating table_assay dropdown..."
    ))
})

test_that("restoreLipidDaFormulaFromSession keeps the missing-formula contract stable", {
    da_data <- new.env(parent = emptyenv())
    update_calls <- list()
    warn_message <- NULL

    result <- restoreLipidDaFormulaFromSession(
        sessionData = list(
            current_s4_object = newMockLipidDaArgsCarrier(args = list(
                daAnalysisParameters = list(formula_string = "")
            ))
        ),
        daData = da_data,
        session = structure(list(kind = "mock-session"), class = "mock_session"),
        updateTextAreaInputFn = function(...) {
            update_calls[[length(update_calls) + 1]] <<- list(...)
        },
        hasArgsMethodFn = function() TRUE,
        logWarn = function(message) {
            warn_message <<- message
        }
    )

    expect_null(result)
    expect_length(update_calls, 0)
    expect_false(exists("formula_from_s4", envir = da_data, inherits = FALSE))
    expect_null(warn_message)
})

test_that("restoreLipidDaFormulaFromSession keeps the formula-restore contract stable", {
    da_data <- new.env(parent = emptyenv())
    update_calls <- list()
    debug_messages <- character()
    mock_session <- structure(list(kind = "mock-session"), class = "mock_session")

    result <- restoreLipidDaFormulaFromSession(
        sessionData = list(
            current_s4_object = newMockLipidDaArgsCarrier(args = list(
                daAnalysisParameters = list(formula_string = "~ 0 + group")
            ))
        ),
        daData = da_data,
        session = mock_session,
        updateTextAreaInputFn = function(session, inputId, value) {
            update_calls[[length(update_calls) + 1]] <<- list(
                session = session,
                inputId = inputId,
                value = value
            )
        },
        debugLog = function(...) {
            debug_messages <<- c(debug_messages, paste0(..., collapse = ""))
        },
        hasArgsMethodFn = function() TRUE,
        logWarn = function(...) stop("warning logger should not run")
    )

    expect_identical(result$formula, "~ 0 + group")
    expect_identical(da_data$formula_from_s4, "~ 0 + group")
    expect_length(update_calls, 1)
    expect_identical(update_calls[[1]]$session, mock_session)
    expect_identical(update_calls[[1]]$inputId, "formula_string")
    expect_identical(update_calls[[1]]$value, "~ 0 + group")
    expect_identical(debug_messages, c(
        "  STEP 9a: Checking S4 @args slot...",
        "    hasSlot 'args': TRUE",
        "  STEP 9b: s4_args is NULL: FALSE",
        "    s4_args class: list",
        "    names(s4_args): daAnalysisParameters"
    ))
})

test_that("restoreLipidDaFormulaFromSession keeps the warning path stable on S4 extraction errors", {
    da_data <- new.env(parent = emptyenv())
    update_calls <- list()
    debug_messages <- character()
    warn_message <- NULL

    result <- restoreLipidDaFormulaFromSession(
        sessionData = list(current_s4_object = NULL),
        daData = da_data,
        session = structure(list(kind = "mock-session"), class = "mock_session"),
        updateTextAreaInputFn = function(...) {
            update_calls[[length(update_calls) + 1]] <<- list(...)
        },
        debugLog = function(...) {
            debug_messages <<- c(debug_messages, paste0(..., collapse = ""))
        },
        hasArgsMethodFn = function() FALSE,
        logWarn = function(message) {
            warn_message <<- message
        }
    )

    expect_null(result)
    expect_length(update_calls, 0)
    expect_false(exists("formula_from_s4", envir = da_data, inherits = FALSE))
    expect_identical(debug_messages[1:2], c(
        "  STEP 9a: Checking S4 @args slot...",
        "    hasSlot 'args': FALSE"
    ))
    expect_match(debug_messages[3], "^  STEP 9 ERROR: ")
    expect_match(warn_message, "^Could not extract formula from S4: ")
})

test_that("finalizeLipidDaSessionLoadSuccess keeps the success-notification contract stable", {
    removed_ids <- character()
    notification <- NULL
    debug_messages <- character()

    result <- finalizeLipidDaSessionLoadSuccess(
        removeNotificationFn = function(id) {
            removed_ids <<- c(removed_ids, id)
        },
        notify = function(message, type, duration) {
            notification <<- list(message = message, type = type, duration = duration)
        },
        debugLog = function(...) {
            debug_messages <<- c(debug_messages, paste0(..., collapse = ""))
        }
    )

    expect_identical(debug_messages, "  STEP 10: Removing notification and showing success...")
    expect_identical(removed_ids, "loading_session")
    expect_identical(notification$message, "Session loaded successfully!")
    expect_identical(notification$type, "message")
    expect_identical(notification$duration, 3)
    expect_identical(result$loadingNotificationId, "loading_session")
    expect_identical(result$successMessage, "Session loaded successfully!")
    expect_identical(result$successType, "message")
    expect_identical(result$successDuration, 3)
})

test_that("notifyLipidDaSessionSourceDirError keeps the preflight source-dir notification contract stable", {
    notification <- NULL
    debug_messages <- character()

    result <- notifyLipidDaSessionSourceDirError(
        notify = function(message, type, duration) {
            notification <<- list(message = message, type = type, duration = duration)
        },
        debugLog = function(...) {
            debug_messages <<- c(debug_messages, paste0(..., collapse = ""))
        }
    )

    expect_identical(debug_messages, "  ERROR: No valid source directory found")
    expect_identical(notification$message, "Could not find source directory for session data.")
    expect_identical(notification$type, "error")
    expect_identical(notification$duration, 5)
    expect_identical(result$loggedMessage, "  ERROR: No valid source directory found")
    expect_identical(result$notificationMessage, "Could not find source directory for session data.")
    expect_identical(result$notificationType, "error")
    expect_identical(result$notificationDuration, 5)
})

test_that("notifyLipidDaSessionFileMissing keeps the missing-session-file notification contract stable", {
    notification <- NULL

    result <- notifyLipidDaSessionFileMissing(
        sessionFile = "/tmp/lipid_filtered_session_data_latest.rds",
        notify = function(message, type, duration) {
            notification <<- list(message = message, type = type, duration = duration)
        }
    )

    expect_identical(result$sessionFile, "/tmp/lipid_filtered_session_data_latest.rds")
    expect_identical(
        result$notificationMessage,
        "Session file not found: /tmp/lipid_filtered_session_data_latest.rds"
    )
    expect_identical(result$notificationType, "error")
    expect_identical(result$notificationDuration, 5)
    expect_identical(
        notification$message,
        "Session file not found: /tmp/lipid_filtered_session_data_latest.rds"
    )
    expect_identical(notification$type, "error")
    expect_identical(notification$duration, 5)
})

test_that("resolveLipidDaSessionFile keeps the source-dir fallback contract stable", {
    debug_messages <- character()
    source_dir_notifications <- list()
    missing_file_notifications <- list()

    result <- resolveLipidDaSessionFile(
        experimentPaths = list(
            source_dir = "/tmp/missing-source",
            export_dir = "/tmp/export-dir"
        ),
        debugLog = function(...) {
            debug_messages <<- c(debug_messages, paste0(..., collapse = ""))
        },
        dirExistsFn = function(path) identical(path, "/tmp/export-dir"),
        fileExistsFn = function(path) identical(
            path,
            "/tmp/export-dir/lipid_filtered_session_data_latest.rds"
        ),
        notifySourceDirFn = function(...) {
            source_dir_notifications[[length(source_dir_notifications) + 1]] <<- list(...)
        },
        notifySessionFileMissingFn = function(...) {
            missing_file_notifications[[length(missing_file_notifications) + 1]] <<- list(...)
        }
    )

    expect_identical(result$sourceDir, "/tmp/export-dir")
    expect_identical(
        result$sessionFile,
        "/tmp/export-dir/lipid_filtered_session_data_latest.rds"
    )
    expect_length(source_dir_notifications, 0)
    expect_length(missing_file_notifications, 0)
    expect_identical(debug_messages, c(
        "  STEP 2: source_dir = /tmp/missing-source",
        "  BRANCH: source_dir NULL or not exists, trying export_dir",
        "  STEP 3: export_dir = /tmp/export-dir",
        "  STEP 4: session_file = /tmp/export-dir/lipid_filtered_session_data_latest.rds",
        "  STEP 4: file.exists = TRUE"
    ))
})

test_that("resolveLipidDaSessionFile keeps the missing-source-dir notification contract stable", {
    debug_messages <- character()
    source_dir_notification <- NULL
    missing_file_notifications <- list()

    result <- resolveLipidDaSessionFile(
        experimentPaths = list(
            source_dir = NULL,
            export_dir = "/tmp/missing-export"
        ),
        debugLog = function(...) {
            debug_messages <<- c(debug_messages, paste0(..., collapse = ""))
        },
        dirExistsFn = function(path) FALSE,
        fileExistsFn = function(...) stop("file existence should not be checked"),
        notifySourceDirFn = function(...) {
            source_dir_notification <<- list(...)
        },
        notifySessionFileMissingFn = function(...) {
            missing_file_notifications[[length(missing_file_notifications) + 1]] <<- list(...)
        }
    )

    expect_null(result)
    expect_type(source_dir_notification$debugLog, "closure")
    expect_identical(debug_messages, c(
        "  STEP 2: source_dir = NULL",
        "  BRANCH: source_dir NULL or not exists, trying export_dir",
        "  STEP 3: export_dir = /tmp/missing-export"
    ))
    expect_length(missing_file_notifications, 0)
})

test_that("resolveLipidDaSessionFile keeps the missing-session-file contract stable", {
    debug_messages <- character()
    source_dir_notifications <- list()
    missing_file_notification <- NULL

    result <- resolveLipidDaSessionFile(
        experimentPaths = list(
            source_dir = "/tmp/source-dir",
            export_dir = "/tmp/export-dir"
        ),
        debugLog = function(...) {
            debug_messages <<- c(debug_messages, paste0(..., collapse = ""))
        },
        dirExistsFn = function(path) identical(path, "/tmp/source-dir"),
        fileExistsFn = function(path) FALSE,
        notifySourceDirFn = function(...) {
            source_dir_notifications[[length(source_dir_notifications) + 1]] <<- list(...)
        },
        notifySessionFileMissingFn = function(...) {
            missing_file_notification <<- list(...)
        }
    )

    expect_null(result)
    expect_length(source_dir_notifications, 0)
    expect_identical(
        missing_file_notification$sessionFile,
        "/tmp/source-dir/lipid_filtered_session_data_latest.rds"
    )
    expect_identical(debug_messages, c(
        "  STEP 2: source_dir = /tmp/source-dir",
        "  STEP 4: session_file = /tmp/source-dir/lipid_filtered_session_data_latest.rds",
        "  STEP 4: file.exists = FALSE"
    ))
})

test_that("showLipidDaSessionLoadingNotification keeps the loading-notification contract stable", {
    notification <- NULL

    result <- showLipidDaSessionLoadingNotification(
        notify = function(message, id, duration) {
            notification <<- list(message = message, id = id, duration = duration)
        }
    )

    expect_identical(result$notificationMessage, "Loading filtered session...")
    expect_identical(result$notificationId, "loading_session")
    expect_null(result$notificationDuration)
    expect_identical(notification$message, "Loading filtered session...")
    expect_identical(notification$id, "loading_session")
    expect_null(notification$duration)
})

test_that("readLipidDaSessionData keeps the session-read logging contract stable", {
    debug_messages <- character()
    info_messages <- character()
    session_payload <- list(
        current_s4_object = structure(list(kind = "mock-s4"), class = "mock_s4"),
        contrasts_tbl = data.frame(contrast = "B_vs_A", stringsAsFactors = FALSE)
    )

    result <- readLipidDaSessionData(
        sessionFile = "/tmp/lipid_filtered_session_data_latest.rds",
        readRdsFn = function(path) {
            expect_identical(path, "/tmp/lipid_filtered_session_data_latest.rds")
            session_payload
        },
        debugLog = function(...) {
            debug_messages <<- c(debug_messages, paste0(..., collapse = ""))
        },
        logInfo = function(message) {
            info_messages <<- c(info_messages, message)
        }
    )

    expect_identical(result$sessionFile, "/tmp/lipid_filtered_session_data_latest.rds")
    expect_identical(result$sessionData, session_payload)
    expect_identical(result$sessionDataNames, c("current_s4_object", "contrasts_tbl"))
    expect_identical(debug_messages, c(
        "  STEP 5: Reading RDS file...",
        "  STEP 5: RDS loaded successfully",
        "    names(session_data): current_s4_object, contrasts_tbl"
    ))
    expect_identical(
        info_messages,
        "   Loaded session from: /tmp/lipid_filtered_session_data_latest.rds"
    )
})

test_that("restoreLipidDaCurrentS4FromSession keeps the empty-S4 restore contract stable", {
    da_data <- new.env(parent = emptyenv())
    workflow_data <- new.env(parent = emptyenv())
    save_state_calls <- list()
    debug_messages <- character()
    info_messages <- character()
    workflow_data$state_manager <- structure(list(
        saveState = function(...) {
            save_state_calls[[length(save_state_calls) + 1]] <<- list(...)
        }
    ), class = "mock_state_manager")

    result <- restoreLipidDaCurrentS4FromSession(
        sessionData = list(current_s4_object = NULL),
        daData = da_data,
        workflowData = workflow_data,
        sessionFile = "/tmp/lipid_filtered_session_data_latest.rds",
        debugLog = function(...) {
            debug_messages <<- c(debug_messages, paste0(..., collapse = ""))
        },
        logInfo = function(message) {
            info_messages <<- c(info_messages, message)
        }
    )

    expect_null(result)
    expect_false(exists("current_s4_object", envir = da_data, inherits = FALSE))
    expect_length(save_state_calls, 0)
    expect_length(debug_messages, 0)
    expect_length(info_messages, 0)
})

test_that("restoreLipidDaCurrentS4FromSession keeps the S4/state-manager restore contract stable", {
    da_data <- new.env(parent = emptyenv())
    workflow_data <- new.env(parent = emptyenv())
    debug_messages <- character()
    info_messages <- character()
    save_state_call <- NULL
    current_s4_object <- structure(list(kind = "mock-s4"), class = "mock_s4")
    workflow_data$state_manager <- structure(list(
        saveState = function(state_name, s4_data_object, config_object, description) {
            save_state_call <<- list(
                state_name = state_name,
                s4_data_object = s4_data_object,
                config_object = config_object,
                description = description
            )
        }
    ), class = "mock_state_manager")

    result <- restoreLipidDaCurrentS4FromSession(
        sessionData = list(
            current_s4_object = current_s4_object,
            r6_current_state_name = NULL
        ),
        daData = da_data,
        workflowData = workflow_data,
        sessionFile = "/tmp/lipid_filtered_session_data_latest.rds",
        debugLog = function(...) {
            debug_messages <<- c(debug_messages, paste0(..., collapse = ""))
        },
        logInfo = function(message) {
            info_messages <<- c(info_messages, message)
        }
    )

    expect_identical(result$currentS4Object, current_s4_object)
    expect_identical(result$stateName, "loaded_for_de")
    expect_identical(result$sessionFile, "/tmp/lipid_filtered_session_data_latest.rds")
    expect_identical(result$description, "Loaded from filtered session for DA analysis")
    expect_identical(da_data$current_s4_object, current_s4_object)
    expect_identical(save_state_call$state_name, "loaded_for_de")
    expect_identical(save_state_call$s4_data_object, current_s4_object)
    expect_identical(
        save_state_call$config_object,
        list(loaded_from = "/tmp/lipid_filtered_session_data_latest.rds")
    )
    expect_identical(
        save_state_call$description,
        "Loaded from filtered session for DA analysis"
    )
    expect_identical(debug_messages, c(
        "  STEP 6a: S4 object class: mock_s4",
        "  STEP 6b: Assigning to da_data$current_s4_object...",
        "  STEP 6c: Getting state_name...",
        "    r6_current_state_name: NULL",
        "  STEP 6d: Checking workflow_data$state_manager...",
        "    workflow_data is NULL: FALSE",
        "    workflow_data class: environment",
        "    state_manager exists: TRUE",
        "    state_manager class: mock_state_manager",
        "  STEP 6e: Calling workflow_data$state_manager$saveState...",
        "  STEP 6f: saveState completed"
    ))
    expect_identical(info_messages, "   Restored S4 object to state manager")
})

test_that("restoreLipidDaPostReadSessionState keeps the empty contrast/assay orchestration contract stable", {
    workflow_data <- new.env(parent = emptyenv())
    da_data <- new.env(parent = emptyenv())
    debug_messages <- character()
    contrast_calls <- list()
    assay_calls <- list()
    formula_calls <- list()
    formula_restore <- list(formula = NULL)

    result <- restoreLipidDaPostReadSessionState(
        sessionData = list(
            contrasts_tbl = NULL,
            assay_names = NULL,
            current_s4_object = NULL
        ),
        workflowData = workflow_data,
        daData = da_data,
        session = "mock-session",
        debugLog = function(...) {
            debug_messages <<- c(debug_messages, paste0(..., collapse = ""))
        },
        restoreContrastsFn = function(...) {
            contrast_calls[[length(contrast_calls) + 1]] <<- list(...)
            stop("contrast restore should not run")
        },
        restoreAssaysFn = function(...) {
            assay_calls[[length(assay_calls) + 1]] <<- list(...)
            stop("assay restore should not run")
        },
        restoreFormulaFn = function(...) {
            formula_calls[[length(formula_calls) + 1]] <<- list(...)
            formula_restore
        }
    )

    expect_null(result$contrastRestore)
    expect_null(result$assayRestore)
    expect_identical(result$formulaRestore, formula_restore)
    expect_length(contrast_calls, 0)
    expect_length(assay_calls, 0)
    expect_length(formula_calls, 1)
    expect_identical(formula_calls[[1]]$sessionData$current_s4_object, NULL)
    expect_identical(formula_calls[[1]]$daData, da_data)
    expect_identical(formula_calls[[1]]$session, "mock-session")
    expect_identical(debug_messages, c(
        "  STEP 7: Checking contrasts_tbl...",
        "    is.null(session_data$contrasts_tbl): TRUE",
        "  STEP 8: Checking assay_names...",
        "    is.null(session_data$assay_names): TRUE",
        "  STEP 9: Attempting to extract formula from S4 args..."
    ))
})

test_that("restoreLipidDaPostReadSessionState keeps the contrast/assay/formula restore shell stable", {
    workflow_data <- new.env(parent = emptyenv())
    da_data <- new.env(parent = emptyenv())
    debug_messages <- character()
    contrast_calls <- list()
    assay_calls <- list()
    formula_calls <- list()
    session_data <- list(
        contrasts_tbl = data.frame(
            friendly_names = c("B vs A", "C vs A"),
            contrasts = c("B_vs_A", "C_vs_A"),
            stringsAsFactors = FALSE
        ),
        assay_names = c("Assay1", "Assay2"),
        current_s4_object = newMockLipidDaArgsCarrier()
    )
    contrast_restore <- list(contrastChoices = c("B vs A", "C vs A"))
    assay_restore <- list(
        assayChoices = c("Combined", "Assay1", "Assay2"),
        tableAssayChoices = c("All", "Assay1", "Assay2")
    )
    formula_restore <- list(formula = "~ 0 + group")

    result <- restoreLipidDaPostReadSessionState(
        sessionData = session_data,
        workflowData = workflow_data,
        daData = da_data,
        session = "mock-session",
        debugLog = function(...) {
            debug_messages <<- c(debug_messages, paste0(..., collapse = ""))
        },
        restoreContrastsFn = function(...) {
            contrast_calls[[length(contrast_calls) + 1]] <<- list(...)
            contrast_restore
        },
        restoreAssaysFn = function(...) {
            assay_calls[[length(assay_calls) + 1]] <<- list(...)
            assay_restore
        },
        restoreFormulaFn = function(...) {
            formula_calls[[length(formula_calls) + 1]] <<- list(...)
            formula_restore
        }
    )

    expect_identical(result$contrastRestore, contrast_restore)
    expect_identical(result$assayRestore, assay_restore)
    expect_identical(result$formulaRestore, formula_restore)
    expect_length(contrast_calls, 1)
    expect_identical(contrast_calls[[1]]$sessionData, session_data)
    expect_identical(contrast_calls[[1]]$workflowData, workflow_data)
    expect_identical(contrast_calls[[1]]$daData, da_data)
    expect_identical(contrast_calls[[1]]$session, "mock-session")
    expect_length(assay_calls, 1)
    expect_identical(assay_calls[[1]]$sessionData, session_data)
    expect_identical(assay_calls[[1]]$daData, da_data)
    expect_identical(assay_calls[[1]]$session, "mock-session")
    expect_type(assay_calls[[1]]$debugLog, "closure")
    expect_length(formula_calls, 1)
    expect_identical(formula_calls[[1]]$sessionData, session_data)
    expect_identical(formula_calls[[1]]$daData, da_data)
    expect_identical(formula_calls[[1]]$session, "mock-session")
    expect_type(formula_calls[[1]]$debugLog, "closure")
    expect_identical(debug_messages, c(
        "  STEP 7: Checking contrasts_tbl...",
        "    is.null(session_data$contrasts_tbl): FALSE",
        "    nrow(contrasts_tbl): 2",
        "    colnames: friendly_names, contrasts",
        "  STEP 7a: Restoring contrast state and dropdowns...",
        "    contrast_choices: B vs A, C vs A",
        "  STEP 8: Checking assay_names...",
        "    is.null(session_data$assay_names): FALSE",
        "    assay_names: Assay1, Assay2",
        "    assay_choices: Combined, Assay1, Assay2",
        "    table_assay_choices: All, Assay1, Assay2",
        "  STEP 9: Attempting to extract formula from S4 args..."
    ))
})

test_that("loadLipidDaSessionFromFile keeps the success orchestration contract stable", {
    workflow_data <- new.env(parent = emptyenv())
    da_data <- new.env(parent = emptyenv())
    debug_messages <- character()
    call_sequence <- character()
    session_payload <- list(
        current_s4_object = newMockLipidDaArgsCarrier(args = list(formula = "~ 0 + group")),
        contrasts_tbl = data.frame(contrasts = "B_vs_A", stringsAsFactors = FALSE),
        assay_names = "Assay1"
    )
    show_loading_contract <- list(notificationId = "loading_session")
    session_load <- list(
        sessionFile = "/tmp/lipid_filtered_session_data_latest.rds",
        sessionData = session_payload
    )
    current_s4_restore <- list(stateName = "loaded_for_de")
    post_read_restore <- list(formulaRestore = list(formula = "~ 0 + group"))
    success_contract <- list(successMessage = "Session loaded successfully!")
    start_time <- as.POSIXct("2026-04-15 12:00:00", tz = "UTC")

    result <- loadLipidDaSessionFromFile(
        sessionFile = "/tmp/lipid_filtered_session_data_latest.rds",
        workflowData = workflow_data,
        daData = da_data,
        session = "mock-session",
        debugLog = function(...) {
            debug_messages <<- c(debug_messages, paste0(..., collapse = ""))
        },
        startTime = start_time,
        showLoadingFn = function() {
            call_sequence <<- c(call_sequence, "showLoading")
            show_loading_contract
        },
        readSessionFn = function(sessionFile, debugLog) {
            call_sequence <<- c(call_sequence, "readSession")
            expect_identical(sessionFile, "/tmp/lipid_filtered_session_data_latest.rds")
            expect_type(debugLog, "closure")
            session_load
        },
        restoreCurrentS4Fn = function(sessionData, daData, workflowData, sessionFile, debugLog) {
            call_sequence <<- c(call_sequence, "restoreCurrentS4")
            expect_identical(sessionData, session_payload)
            expect_identical(daData, da_data)
            expect_identical(workflowData, workflow_data)
            expect_identical(sessionFile, "/tmp/lipid_filtered_session_data_latest.rds")
            expect_type(debugLog, "closure")
            current_s4_restore
        },
        restorePostReadFn = function(sessionData, workflowData, daData, session, debugLog) {
            call_sequence <<- c(call_sequence, "restorePostRead")
            expect_identical(sessionData, session_payload)
            expect_identical(workflowData, workflow_data)
            expect_identical(daData, da_data)
            expect_identical(session, "mock-session")
            expect_type(debugLog, "closure")
            post_read_restore
        },
        finalizeSuccessFn = function(debugLog) {
            call_sequence <<- c(call_sequence, "finalizeSuccess")
            expect_type(debugLog, "closure")
            success_contract
        },
        finalizeErrorFn = function(...) {
            call_sequence <<- c(call_sequence, "finalizeError")
            stop("error finalizer should not run")
        },
        nowFn = function() {
            start_time + 2.5
        }
    )

    expect_identical(call_sequence, c(
        "showLoading",
        "readSession",
        "restoreCurrentS4",
        "restorePostRead",
        "finalizeSuccess"
    ))
    expect_identical(result$status, "success")
    expect_identical(result$showLoadingContract, show_loading_contract)
    expect_identical(result$sessionLoad, session_load)
    expect_identical(result$currentS4Restore, current_s4_restore)
    expect_identical(result$postReadRestore, post_read_restore)
    expect_identical(result$successContract, success_contract)
    expect_identical(as.numeric(result$elapsedSeconds, units = "secs"), 2.5)
    expect_identical(debug_messages, c(
        "  STEP 6: Checking current_s4_object...",
        "    is.null(session_data$current_s4_object): FALSE",
        "--- EXIT load_filtered_session (2.5s) --- SUCCESS"
    ))
})

test_that("loadLipidDaSessionFromFile keeps the fatal-error orchestration contract stable", {
    debug_messages <- character()
    call_sequence <- character()
    error_contract <- list(notificationMessage = "Error loading session: boom")

    result <- loadLipidDaSessionFromFile(
        sessionFile = "/tmp/lipid_filtered_session_data_latest.rds",
        workflowData = new.env(parent = emptyenv()),
        daData = new.env(parent = emptyenv()),
        session = "mock-session",
        debugLog = function(...) {
            debug_messages <<- c(debug_messages, paste0(..., collapse = ""))
        },
        showLoadingFn = function() {
            call_sequence <<- c(call_sequence, "showLoading")
            list(notificationId = "loading_session")
        },
        readSessionFn = function(...) {
            call_sequence <<- c(call_sequence, "readSession")
            stop(simpleError("boom", call = quote(readRDS(sessionFile))))
        },
        restoreCurrentS4Fn = function(...) {
            call_sequence <<- c(call_sequence, "restoreCurrentS4")
            stop("S4 restore should not run")
        },
        restorePostReadFn = function(...) {
            call_sequence <<- c(call_sequence, "restorePostRead")
            stop("post-read restore should not run")
        },
        finalizeSuccessFn = function(...) {
            call_sequence <<- c(call_sequence, "finalizeSuccess")
            stop("success finalizer should not run")
        },
        finalizeErrorFn = function(errorMessage) {
            call_sequence <<- c(call_sequence, "finalizeError")
            expect_identical(errorMessage, "boom")
            error_contract
        }
    )

    expect_identical(call_sequence, c("showLoading", "readSession", "finalizeError"))
    expect_identical(result$status, "error")
    expect_identical(result$errorMessage, "boom")
    expect_identical(result$errorCall, "readRDS(sessionFile)")
    expect_identical(result$errorContract, error_contract)
    expect_identical(debug_messages, c(
        "  FATAL ERROR: boom",
        "  ERROR CALL: readRDS(sessionFile)",
        "--- EXIT load_filtered_session --- FAILED"
    ))
})

test_that("handleLipidDaLoadFilteredSession keeps the missing-session-source contract stable", {
    workflow_data <- new.env(parent = emptyenv())
    da_data <- new.env(parent = emptyenv())
    debug_messages <- character()
    info_messages <- character()
    call_sequence <- character()

    result <- handleLipidDaLoadFilteredSession(
        experimentPaths = list(
            source_dir = "/tmp/source-dir",
            export_dir = "/tmp/export-dir"
        ),
        workflowData = workflow_data,
        daData = da_data,
        session = "mock-session",
        debugLog = function(...) {
            debug_messages <<- c(debug_messages, paste0(..., collapse = ""))
        },
        logInfo = function(message) {
            info_messages <<- c(info_messages, message)
        },
        resolveSessionSourceFn = function(experimentPaths, debugLog) {
            call_sequence <<- c(call_sequence, "resolveSessionSource")
            expect_identical(experimentPaths$source_dir, "/tmp/source-dir")
            expect_identical(experimentPaths$export_dir, "/tmp/export-dir")
            expect_type(debugLog, "closure")
            NULL
        },
        loadSessionFromFileFn = function(...) {
            call_sequence <<- c(call_sequence, "loadSessionFromFile")
            stop("session loader should not run")
        }
    )

    expect_null(result)
    expect_identical(call_sequence, "resolveSessionSource")
    expect_identical(info_messages, "=== LOAD FILTERED SESSION BUTTON CLICKED ===")
    expect_identical(debug_messages, c(
        "  STEP 1: Getting source_dir from experiment_paths",
        "    experiment_paths class: list",
        "    experiment_paths is NULL: FALSE"
    ))
})

test_that("handleLipidDaLoadFilteredSession keeps the session-source handoff contract stable", {
    workflow_data <- new.env(parent = emptyenv())
    da_data <- new.env(parent = emptyenv())
    debug_messages <- character()
    info_messages <- character()
    call_sequence <- character()
    session_source <- list(
        sourceDir = "/tmp/source-dir",
        sessionFile = "/tmp/source-dir/lipid_filtered_session_data_latest.rds"
    )
    load_result <- list(status = "success")
    start_time <- as.POSIXct("2026-04-15 12:00:00", tz = "UTC")

    result <- handleLipidDaLoadFilteredSession(
        experimentPaths = list(
            source_dir = "/tmp/source-dir",
            export_dir = "/tmp/export-dir"
        ),
        workflowData = workflow_data,
        daData = da_data,
        session = "mock-session",
        debugLog = function(...) {
            debug_messages <<- c(debug_messages, paste0(..., collapse = ""))
        },
        startTime = start_time,
        logInfo = function(message) {
            info_messages <<- c(info_messages, message)
        },
        resolveSessionSourceFn = function(experimentPaths, debugLog) {
            call_sequence <<- c(call_sequence, "resolveSessionSource")
            expect_identical(experimentPaths$source_dir, "/tmp/source-dir")
            expect_type(debugLog, "closure")
            session_source
        },
        loadSessionFromFileFn = function(sessionFile, workflowData, daData, session, debugLog, startTime) {
            call_sequence <<- c(call_sequence, "loadSessionFromFile")
            expect_identical(sessionFile, session_source$sessionFile)
            expect_identical(workflowData, workflow_data)
            expect_identical(daData, da_data)
            expect_identical(session, "mock-session")
            expect_type(debugLog, "closure")
            expect_identical(startTime, start_time)
            load_result
        }
    )

    expect_identical(call_sequence, c("resolveSessionSource", "loadSessionFromFile"))
    expect_identical(result$sessionSource, session_source)
    expect_identical(result$loadResult, load_result)
    expect_identical(info_messages, "=== LOAD FILTERED SESSION BUTTON CLICKED ===")
    expect_identical(debug_messages, c(
        "  STEP 1: Getting source_dir from experiment_paths",
        "    experiment_paths class: list",
        "    experiment_paths is NULL: FALSE"
    ))
})

test_that("bootstrapLipidDaLoadFilteredSession keeps the D66 bootstrap contract stable", {
    workflow_data <- new.env(parent = emptyenv())
    da_data <- new.env(parent = emptyenv())
    emitted_messages <- character()
    handled_call <- NULL
    start_time <- as.POSIXct("2026-04-15 13:45:00", tz = "UTC")
    handled_result <- list(status = "success")

    result <- bootstrapLipidDaLoadFilteredSession(
        experimentPaths = list(source_dir = "/tmp/source-dir"),
        workflowData = workflow_data,
        daData = da_data,
        session = "mock-session",
        messageFn = function(message) {
            emitted_messages <<- c(emitted_messages, message)
        },
        nowFn = function() {
            start_time
        },
        handleLoadFn = function(experimentPaths, workflowData, daData, session, debugLog, startTime) {
            handled_call <<- list(
                experimentPaths = experimentPaths,
                workflowData = workflowData,
                daData = daData,
                session = session,
                debugLog = debugLog,
                startTime = startTime
            )
            debugLog("nested message")
            handled_result
        }
    )

    expect_identical(result, handled_result)
    expect_identical(handled_call$experimentPaths, list(source_dir = "/tmp/source-dir"))
    expect_identical(handled_call$workflowData, workflow_data)
    expect_identical(handled_call$daData, da_data)
    expect_identical(handled_call$session, "mock-session")
    expect_type(handled_call$debugLog, "closure")
    expect_identical(handled_call$startTime, start_time)
    expect_identical(emitted_messages, c(
        "[D66] --- ENTER load_filtered_session observer ---",
        "[D66] nested message"
    ))
})

test_that("bootstrapLipidDaRunAnalysis keeps the observer bootstrap handoff stable", {
    da_data <- new.env(parent = emptyenv())
    workflow_data <- new.env(parent = emptyenv())
    analysis_context <- list(
        currentS4 = "mock-s4",
        contrastsTbl = data.frame(
            contrasts = "B_vs_A",
            stringsAsFactors = FALSE
        )
    )
    prepare_call <- NULL
    preflight_call <- NULL

    result <- bootstrapLipidDaRunAnalysis(
        daData = da_data,
        workflowData = workflow_data,
        formulaString = "~ condition",
        daQValThresh = 0.05,
        treatLfcCutoff = 1,
        session = "mock-session",
        experimentPaths = list(
            da_output_dir = "/tmp/da-output",
            publication_graphs_dir = "/tmp/publication-graphs"
        ),
        prepareAnalysisContextFn = function(daData, workflowData) {
            prepare_call <<- list(
                daData = daData,
                workflowData = workflowData
            )
            analysis_context
        },
        handlePreflightFn = function(...) {
            preflight_call <<- list(...)
            list(marker = "preflight-result")
        }
    )

    expect_identical(prepare_call$daData, da_data)
    expect_identical(prepare_call$workflowData, workflow_data)
    expect_identical(preflight_call$analysisContext, analysis_context)
    expect_identical(preflight_call$formulaString, "~ condition")
    expect_identical(preflight_call$daQValThresh, 0.05)
    expect_identical(preflight_call$treatLfcCutoff, 1)
    expect_identical(preflight_call$daData, da_data)
    expect_identical(preflight_call$workflowData, workflow_data)
    expect_identical(preflight_call$session, "mock-session")
    expect_identical(preflight_call$experimentPaths, list(
        da_output_dir = "/tmp/da-output",
        publication_graphs_dir = "/tmp/publication-graphs"
    ))
    expect_identical(result, list(marker = "preflight-result"))
})

test_that("prepareLipidDaRunAnalysisContext keeps the missing-data notification contract stable", {
    da_data <- new.env(parent = emptyenv())
    workflow_data <- new.env(parent = emptyenv())
    workflow_data$contrasts_tbl <- data.frame(
        contrasts = "B_vs_A",
        stringsAsFactors = FALSE
    )
    info_messages <- character()
    notifications <- list()
    resolve_calls <- 0L

    result <- prepareLipidDaRunAnalysisContext(
        daData = da_data,
        workflowData = workflow_data,
        notify = function(message, type = NULL, duration = NULL, id = NULL) {
            notifications[[length(notifications) + 1]] <<- list(
                message = message,
                type = type,
                duration = duration,
                id = id
            )
        },
        logInfo = function(message) {
            info_messages <<- c(info_messages, message)
        },
        resolveCurrentStateFn = function(workflowData) {
            resolve_calls <<- resolve_calls + 1L
            expect_identical(workflowData$contrasts_tbl$contrasts, "B_vs_A")
            NULL
        }
    )

    expect_null(result)
    expect_identical(resolve_calls, 1L)
    expect_identical(info_messages, "=== RUN DA ANALYSIS BUTTON CLICKED ===")
    expect_length(notifications, 1)
    expect_identical(notifications[[1]], list(
        message = "No lipidomics data loaded. Please load a filtered session first.",
        type = "error",
        duration = 5,
        id = NULL
    ))
})

test_that("prepareLipidDaRunAnalysisContext keeps the missing-contrasts notification contract stable", {
    da_data <- new.env(parent = emptyenv())
    da_data$current_s4_object <- structure(list(marker = TRUE), class = "LipidomicsAssayData")
    workflow_data <- new.env(parent = emptyenv())
    workflow_data$contrasts_tbl <- data.frame()
    notifications <- list()
    resolve_calls <- 0L

    result <- prepareLipidDaRunAnalysisContext(
        daData = da_data,
        workflowData = workflow_data,
        notify = function(message, type = NULL, duration = NULL, id = NULL) {
            notifications[[length(notifications) + 1]] <<- list(
                message = message,
                type = type,
                duration = duration,
                id = id
            )
        },
        logInfo = function(...) NULL,
        resolveCurrentStateFn = function(...) {
            resolve_calls <<- resolve_calls + 1L
            stop("resolver should not run when da_data already has an S4 object")
        }
    )

    expect_null(result)
    expect_identical(resolve_calls, 0L)
    expect_length(notifications, 1)
    expect_identical(notifications[[1]], list(
        message = "No contrasts defined. Please define contrasts in the design tab.",
        type = "error",
        duration = 5,
        id = NULL
    ))
})

test_that("prepareLipidDaRunAnalysisContext keeps the state-resolution and running-notification contract stable", {
    da_data <- new.env(parent = emptyenv())
    workflow_data <- new.env(parent = emptyenv())
    workflow_data$contrasts_tbl <- data.frame(
        contrasts = c("B_vs_A", "C_vs_A"),
        stringsAsFactors = FALSE
    )
    resolved_s4 <- structure(list(marker = TRUE), class = "LipidomicsAssayData")
    notifications <- list()
    info_messages <- character()
    resolve_calls <- 0L

    result <- prepareLipidDaRunAnalysisContext(
        daData = da_data,
        workflowData = workflow_data,
        notify = function(message, type = NULL, duration = NULL, id = NULL) {
            notifications[[length(notifications) + 1]] <<- list(
                message = message,
                type = type,
                duration = duration,
                id = id
            )
        },
        logInfo = function(message) {
            info_messages <<- c(info_messages, message)
        },
        resolveCurrentStateFn = function(workflowData) {
            resolve_calls <<- resolve_calls + 1L
            expect_identical(
                workflowData$contrasts_tbl$contrasts,
                c("B_vs_A", "C_vs_A")
            )
            resolved_s4
        }
    )

    expect_identical(resolve_calls, 1L)
    expect_identical(result$currentS4, resolved_s4)
    expect_identical(result$contrastsTbl, workflow_data$contrasts_tbl)
    expect_identical(result$runningNotification, list(
        message = "Running differential expression analysis...",
        id = "da_running",
        duration = NULL
    ))
    expect_identical(info_messages, "=== RUN DA ANALYSIS BUTTON CLICKED ===")
    expect_length(notifications, 1)
    expect_identical(notifications[[1]], list(
        message = "Running differential expression analysis...",
        type = NULL,
        duration = NULL,
        id = "da_running"
    ))
})

test_that("handleLipidDaRunAnalysisPreflight keeps the NULL early-return contract stable", {
    execute_called <- FALSE

    result <- handleLipidDaRunAnalysisPreflight(
        analysisContext = NULL,
        formulaString = "~ condition",
        daQValThresh = 0.05,
        treatLfcCutoff = 1,
        daData = "mock-da-data",
        workflowData = "mock-workflow-data",
        session = "mock-session",
        experimentPaths = list(da_output_dir = "/tmp/da-output"),
        executeRunAnalysisFn = function(...) {
            execute_called <<- TRUE
        }
    )

    expect_false(execute_called)
    expect_null(result)
})

test_that("handleLipidDaRunAnalysisPreflight keeps the execution handoff stable", {
    execute_call <- NULL

    returned <- handleLipidDaRunAnalysisPreflight(
        analysisContext = list(
            currentS4 = "mock-s4",
            contrastsTbl = data.frame(
                contrasts = "B_vs_A",
                stringsAsFactors = FALSE
            )
        ),
        formulaString = "~ condition",
        daQValThresh = 0.05,
        treatLfcCutoff = 1,
        daData = "mock-da-data",
        workflowData = "mock-workflow-data",
        session = "mock-session",
        experimentPaths = list(
            da_output_dir = "/tmp/da-output",
            publication_graphs_dir = "/tmp/publication-graphs"
        ),
        executeRunAnalysisFn = function(...) {
            execute_call <<- list(...)
            list(marker = "executed")
        }
    )

    expect_identical(execute_call$currentS4, "mock-s4")
    expect_identical(execute_call$contrastsTbl$contrasts, "B_vs_A")
    expect_identical(execute_call$formulaString, "~ condition")
    expect_identical(execute_call$daQValThresh, 0.05)
    expect_identical(execute_call$treatLfcCutoff, 1)
    expect_identical(execute_call$daData, "mock-da-data")
    expect_identical(execute_call$workflowData, "mock-workflow-data")
    expect_identical(execute_call$session, "mock-session")
    expect_identical(execute_call$experimentPaths, list(
        da_output_dir = "/tmp/da-output",
        publication_graphs_dir = "/tmp/publication-graphs"
    ))
    expect_identical(returned, list(marker = "executed"))
})

test_that("executeLipidDaRunAnalysis keeps the analysis-call and success-finalizer handoff stable", {
    run_call <- NULL
    success_call <- NULL
    results <- list(da_lipids_long = data.frame(marker = "ok", stringsAsFactors = FALSE))

    returned <- executeLipidDaRunAnalysis(
        currentS4 = "mock-s4",
        contrastsTbl = data.frame(
            contrasts = "B_vs_A",
            stringsAsFactors = FALSE
        ),
        formulaString = "~ condition",
        daQValThresh = 0.05,
        treatLfcCutoff = 1,
        daData = "mock-da-data",
        workflowData = "mock-workflow-data",
        session = "mock-session",
        experimentPaths = list(
            da_output_dir = "/tmp/da-output",
            publication_graphs_dir = "/tmp/publication-graphs"
        ),
        runLipidsDaFn = function(...) {
            run_call <<- list(...)
            results
        },
        successFinalizerFn = function(...) {
            success_call <<- list(...)
            list(marker = "success-finalizer")
        },
        errorFinalizerFn = function(...) {
            stop("error finalizer should not run")
        }
    )

    expect_identical(run_call$theObject, "mock-s4")
    expect_identical(run_call$contrasts_tbl$contrasts, "B_vs_A")
    expect_identical(run_call$formula_string, "~ condition")
    expect_identical(run_call$da_q_val_thresh, 0.05)
    expect_identical(run_call$treat_lfc_cutoff, 1)
    expect_true(run_call$eBayes_trend)
    expect_true(run_call$eBayes_robust)

    expect_identical(success_call$results, results)
    expect_identical(success_call$daData, "mock-da-data")
    expect_identical(success_call$workflowData, "mock-workflow-data")
    expect_identical(success_call$session, "mock-session")
    expect_identical(success_call$experimentPaths, list(
        da_output_dir = "/tmp/da-output",
        publication_graphs_dir = "/tmp/publication-graphs"
    ))
    expect_identical(success_call$daQValThresh, 0.05)
    expect_identical(success_call$treatLfcCutoff, 1)

    expect_identical(returned, list(
        status = "success",
        results = results,
        finalizerResult = list(marker = "success-finalizer")
    ))
})

test_that("executeLipidDaRunAnalysis keeps the error-finalizer handoff stable", {
    error_call <- NULL

    returned <- executeLipidDaRunAnalysis(
        currentS4 = "mock-s4",
        contrastsTbl = data.frame(
            contrasts = "B_vs_A",
            stringsAsFactors = FALSE
        ),
        formulaString = "~ condition",
        daQValThresh = 0.05,
        treatLfcCutoff = 1,
        daData = "mock-da-data",
        workflowData = "mock-workflow-data",
        session = "mock-session",
        experimentPaths = list(
            da_output_dir = "/tmp/da-output",
            publication_graphs_dir = "/tmp/publication-graphs"
        ),
        runLipidsDaFn = function(...) {
            stop("boom")
        },
        successFinalizerFn = function(...) {
            stop("success finalizer should not run")
        },
        errorFinalizerFn = function(...) {
            error_call <<- list(...)
            list(marker = "error-finalizer")
        }
    )

    expect_identical(error_call, list(errorMessage = "boom"))
    expect_identical(returned, list(
        status = "error",
        errorMessage = "boom",
        finalizerResult = list(marker = "error-finalizer")
    ))
})

test_that("finalizeLipidDaRunAnalysisSuccess keeps the success/update and disk-write contract stable", {
    da_data <- new.env(parent = emptyenv())
    workflow_data <- new.env(parent = emptyenv())
    workflow_data$tab_status <- list(import = "complete", differential_analysis = "pending")
    notifications <- list()
    removed_ids <- character()
    info_messages <- character()
    warn_messages <- character()
    update_calls <- list()
    write_call <- NULL
    results <- list(
        da_lipids_long = data.frame(
            friendly_name = c("B vs A", "B vs A", "C vs A"),
            comparison = c("B_vs_A", "B_vs_A", "C_vs_A"),
            assay = c("LCMS_Pos", "LCMS_Neg", "LCMS_Pos"),
            stringsAsFactors = FALSE
        )
    )

    result <- finalizeLipidDaRunAnalysisSuccess(
        results = results,
        daData = da_data,
        workflowData = workflow_data,
        session = "mock-session",
        experimentPaths = list(
            da_output_dir = "/tmp/da-output",
            publication_graphs_dir = "/tmp/publication-graphs"
        ),
        daQValThresh = 0.05,
        treatLfcCutoff = 1,
        removeNotificationFn = function(id) {
            removed_ids <<- c(removed_ids, id)
        },
        notify = function(message, type = NULL, duration = NULL, id = NULL) {
            notifications[[length(notifications) + 1]] <<- list(
                message = message,
                type = type,
                duration = duration,
                id = id
            )
        },
        updateSelectInputFn = function(session, inputId, choices, selected) {
            update_calls[[length(update_calls) + 1]] <<- list(
                session = session,
                inputId = inputId,
                choices = choices,
                selected = selected
            )
        },
        logInfo = function(message) {
            info_messages <<- c(info_messages, message)
        },
        logWarn = function(message) {
            warn_messages <<- c(warn_messages, message)
        },
        writeResultsFn = function(...) {
            write_call <<- list(...)
            TRUE
        }
    )

    expect_identical(da_data$da_results_list, results)
    expect_true(da_data$analysis_complete)
    expect_identical(workflow_data$tab_status, list(
        import = "complete",
        differential_analysis = "complete"
    ))
    expect_identical(removed_ids, "da_running")
    expect_identical(notifications, list(
        list(
            message = "Differential expression analysis complete!",
            type = "message",
            duration = 5,
            id = NULL
        ),
        list(
            message = "DA results saved to disk (tables, volcano plots, heatmaps)",
            type = "message",
            duration = 5,
            id = NULL
        )
    ))
    expect_identical(update_calls, list(
        list(
            session = "mock-session",
            inputId = "volcano_contrast",
            choices = c("B vs A", "C vs A"),
            selected = "B vs A"
        ),
        list(
            session = "mock-session",
            inputId = "heatmap_contrast",
            choices = c("B vs A", "C vs A"),
            selected = "B vs A"
        ),
        list(
            session = "mock-session",
            inputId = "table_contrast",
            choices = c("B vs A", "C vs A"),
            selected = "B vs A"
        ),
        list(
            session = "mock-session",
            inputId = "volcano_assay",
            choices = c("Combined", "LCMS_Pos", "LCMS_Neg"),
            selected = "Combined"
        ),
        list(
            session = "mock-session",
            inputId = "heatmap_assay",
            choices = c("Combined", "LCMS_Pos", "LCMS_Neg"),
            selected = "Combined"
        ),
        list(
            session = "mock-session",
            inputId = "table_assay",
            choices = c("All", "LCMS_Pos", "LCMS_Neg"),
            selected = "All"
        )
    ))
    expect_identical(info_messages, c(
        "   DA analysis completed successfully",
        "   Updated dropdowns: 2 contrasts, 2 assays",
        "   Writing DA results to disk...",
        "   da_output_dir = /tmp/da-output",
        "   publication_graphs_dir = /tmp/publication-graphs",
        "   All DA results written to disk successfully"
    ))
    expect_identical(warn_messages, character())
    expect_identical(write_call$da_results_list, results)
    expect_identical(write_call$da_output_dir, "/tmp/da-output")
    expect_identical(write_call$publication_graphs_dir, "/tmp/publication-graphs")
    expect_identical(write_call$da_q_val_thresh, 0.05)
    expect_identical(write_call$lfc_threshold, 1)
    expect_identical(write_call$heatmap_top_n, 50)
    expect_identical(write_call$heatmap_clustering, "both")
    expect_identical(write_call$heatmap_color_scheme, "RdBu")
    expect_identical(result$updatedStatus, workflow_data$tab_status)
    expect_identical(result$contrastChoices, c("B vs A", "C vs A"))
    expect_identical(result$assayChoices, c("Combined", "LCMS_Pos", "LCMS_Neg"))
    expect_identical(result$tableAssayChoices, c("All", "LCMS_Pos", "LCMS_Neg"))
    expect_identical(result$diskWrite, list(status = "completed", success = TRUE))
})

test_that("finalizeLipidDaRunAnalysisSuccess keeps the disk-write warning contract stable", {
    da_data <- new.env(parent = emptyenv())
    workflow_data <- new.env(parent = emptyenv())
    workflow_data$tab_status <- list()
    notifications <- list()
    removed_ids <- character()
    info_messages <- character()
    warn_messages <- character()
    update_calls <- list()

    result <- finalizeLipidDaRunAnalysisSuccess(
        results = list(),
        daData = da_data,
        workflowData = workflow_data,
        session = "mock-session",
        experimentPaths = list(
            da_output_dir = "/tmp/da-output",
            publication_graphs_dir = "/tmp/publication-graphs"
        ),
        daQValThresh = 0.05,
        treatLfcCutoff = 1,
        removeNotificationFn = function(id) {
            removed_ids <<- c(removed_ids, id)
        },
        notify = function(message, type = NULL, duration = NULL, id = NULL) {
            notifications[[length(notifications) + 1]] <<- list(
                message = message,
                type = type,
                duration = duration,
                id = id
            )
        },
        updateSelectInputFn = function(...) {
            update_calls[[length(update_calls) + 1]] <<- list(...)
        },
        logInfo = function(message) {
            info_messages <<- c(info_messages, message)
        },
        logWarn = function(message) {
            warn_messages <<- c(warn_messages, message)
        },
        writeResultsFn = function(...) {
            stop("disk full")
        }
    )

    expect_identical(da_data$da_results_list, list())
    expect_true(da_data$analysis_complete)
    expect_identical(workflow_data$tab_status, list(differential_analysis = "complete"))
    expect_identical(removed_ids, "da_running")
    expect_identical(update_calls, list())
    expect_identical(notifications, list(
        list(
            message = "Differential expression analysis complete!",
            type = "message",
            duration = 5,
            id = NULL
        ),
        list(
            message = "Warning: Could not save results to disk: disk full",
            type = "warning",
            duration = 8,
            id = NULL
        )
    ))
    expect_identical(info_messages, c(
        "   DA analysis completed successfully",
        "   Writing DA results to disk...",
        "   da_output_dir = /tmp/da-output",
        "   publication_graphs_dir = /tmp/publication-graphs"
    ))
    expect_identical(warn_messages, "   Could not write DA results to disk: disk full")
    expect_identical(result$updatedStatus, workflow_data$tab_status)
    expect_null(result$contrastChoices)
    expect_null(result$assayChoices)
    expect_null(result$tableAssayChoices)
    expect_identical(result$diskWrite, list(
        status = "error",
        success = FALSE,
        message = "disk full"
    ))
})

test_that("finalizeLipidDaRunAnalysisError keeps the error-finalizer contract stable", {
    removed_ids <- character()
    notifications <- list()
    logged_messages <- character()

    result <- finalizeLipidDaRunAnalysisError(
        errorMessage = "boom",
        removeNotificationFn = function(id) {
            removed_ids <<- c(removed_ids, id)
        },
        notify = function(message, type = NULL, duration = NULL, id = NULL) {
            notifications[[length(notifications) + 1]] <<- list(
                message = message,
                type = type,
                duration = duration,
                id = id
            )
        },
        logError = function(message) {
            logged_messages <<- c(logged_messages, message)
        }
    )

    expect_identical(logged_messages, "   DA analysis error: boom")
    expect_identical(removed_ids, "da_running")
    expect_identical(notifications, list(
        list(
            message = "Analysis error: boom",
            type = "error",
            duration = 10,
            id = NULL
        )
    ))
    expect_identical(result, list(
        errorMessage = "boom",
        loggedMessage = "   DA analysis error: boom",
        removedNotificationId = "da_running",
        notificationMessage = "Analysis error: boom",
        notificationType = "error",
        notificationDuration = 10
    ))
})

test_that("finalizeLipidDaSessionLoadError keeps the fatal-error notification contract stable", {
    removed_ids <- character()
    logged_messages <- character()
    notification <- NULL

    result <- finalizeLipidDaSessionLoadError(
        errorMessage = "boom",
        removeNotificationFn = function(id) {
            removed_ids <<- c(removed_ids, id)
        },
        notify = function(message, type, duration) {
            notification <<- list(message = message, type = type, duration = duration)
        },
        logError = function(message) {
            logged_messages <<- c(logged_messages, message)
        }
    )

    expect_identical(logged_messages, "   Error loading session: boom")
    expect_identical(removed_ids, "loading_session")
    expect_identical(notification$message, "Error loading session: boom")
    expect_identical(notification$type, "error")
    expect_identical(notification$duration, 10)
    expect_identical(result$loadingNotificationId, "loading_session")
    expect_identical(result$loggedMessage, "   Error loading session: boom")
    expect_identical(result$notificationMessage, "Error loading session: boom")
    expect_identical(result$notificationType, "error")
    expect_identical(result$notificationDuration, 10)
})

test_that("buildLipidDaSummaryStatsRenderText keeps the render-shell handoff stable", {
    req_arg <- NULL
    build_call <- NULL
    da_results <- data.frame(
        comparison = "B_vs_A",
        friendly_name = "B vs A",
        assay = "Assay1",
        significant = "Up",
        stringsAsFactors = FALSE
    )
    da_results_list <- list(
        da_lipids_long = da_results,
        marker = "results-list"
    )

    result <- buildLipidDaSummaryStatsRenderText(
        daResultsList = da_results_list,
        selectedContrast = "B vs A",
        selectedAssay = "Assay1",
        daQValThresh = 0.05,
        reqFn = function(value) {
            req_arg <<- value
            value
        },
        buildSummaryStatsTextFn = function(...) {
            build_call <<- list(...)
            "summary-stats-text"
        }
    )

    expect_identical(req_arg, da_results_list)
    expect_identical(build_call$daResults, da_results)
    expect_identical(build_call$selectedContrast, "B vs A")
    expect_identical(build_call$selectedAssay, "Assay1")
    expect_identical(build_call$daQValThresh, 0.05)
    expect_identical(result, "summary-stats-text")
})

test_that("buildLipidDaResultsTableRenderWidget keeps the render-shell handoff stable", {
    req_arg <- NULL
    build_call <- NULL
    da_results <- data.frame(
        lipid_id = "L1",
        comparison = "B_vs_A",
        friendly_name = "B vs A",
        assay = "Assay1",
        stringsAsFactors = FALSE
    )
    da_results_list <- list(
        da_lipids_long = da_results,
        marker = "results-list"
    )

    result <- buildLipidDaResultsTableRenderWidget(
        daResultsList = da_results_list,
        selectedContrast = "B vs A",
        selectedAssay = "Assay1",
        tableSignificance = "significant",
        daQValThresh = 0.05,
        lfcThreshold = 0.5,
        tableMaxRows = 25,
        reqFn = function(value) {
            req_arg <<- value
            value
        },
        buildResultsTableWidgetFn = function(...) {
            build_call <<- list(...)
            "results-table-widget"
        }
    )

    expect_identical(req_arg, da_results_list)
    expect_identical(build_call$daResults, da_results)
    expect_identical(build_call$selectedContrast, "B vs A")
    expect_identical(build_call$selectedAssay, "Assay1")
    expect_identical(build_call$tableSignificance, "significant")
    expect_identical(build_call$daQValThresh, 0.05)
    expect_identical(build_call$lfcThreshold, 0.5)
    expect_identical(build_call$tableMaxRows, 25)
    expect_identical(result, "results-table-widget")
})

test_that("buildLipidDaSummaryStatsText keeps the empty-results guidance stable", {
    expect_identical(
        buildLipidDaSummaryStatsText(
            daResults = NULL,
            selectedContrast = "B_vs_A",
            selectedAssay = "Assay1",
            daQValThresh = 0.05
        ),
        "No results available."
    )
})

test_that("buildLipidDaSummaryStatsText keeps the DA summary filtering and formatting stable", {
    da_results <- data.frame(
        comparison = c("B_vs_A", "B_vs_A", "C_vs_A", "B_vs_A"),
        friendly_name = c("B vs A", "B vs A", "C vs A", "B vs A"),
        assay = c("Assay1", "Assay1", "Assay1", "Assay2"),
        significant = c("Up", "NS", "Down", "Down"),
        stringsAsFactors = FALSE
    )

    summary_text <- buildLipidDaSummaryStatsText(
        daResults = da_results,
        selectedContrast = "B vs A",
        selectedAssay = "Assay1",
        daQValThresh = 0.05
    )

    expect_identical(
        summary_text,
        paste(
            c(
                "Total lipids: 2",
                "Significant (Q < 0.050): 1 (50.0%)",
                "  Up-regulated: 1",
                "  Down-regulated: 0"
            ),
            collapse = "\n"
        )
    )
})

test_that("buildLipidDaResultsTableWidget keeps the empty-results contract stable", {
    expect_null(
        buildLipidDaResultsTableWidget(
            daResults = NULL,
            selectedContrast = "B vs A",
            selectedAssay = "Assay1",
            tableSignificance = "all",
            daQValThresh = 0.05,
            lfcThreshold = 0.5,
            tableMaxRows = 25
        )
    )
})

test_that("buildLipidDaResultsTableWidget keeps filtering and DT formatting stable", {
    datatable_calls <- NULL
    style_equal_calls <- NULL

    table_widget <- buildLipidDaResultsTableWidget(
        daResults = data.frame(
            comparison = c("B_vs_A", "B_vs_A", "B_vs_A", "C_vs_A"),
            friendly_name = c("B vs A", "B vs A", "B vs A", "C vs A"),
            assay = c("Assay1", "Assay1", "Assay2", "Assay1"),
            lipid_id = c("L1", "L2", "L3", "L4"),
            lipid_name = c("Lipid 1", "Lipid 2", "Lipid 3", "Lipid 4"),
            logFC = c(1.2, 0.8, 1.5, -1.1),
            raw_pvalue = c(0.001, 0.002, 0.003, 0.004),
            fdr_qvalue = c(0.01, 0.03, 0.02, 0.01),
            significant = c("Up", "Up", "Up", "Down"),
            extra_column = c("drop1", "drop2", "drop3", "drop4"),
            stringsAsFactors = FALSE
        ),
        selectedContrast = "B vs A",
        selectedAssay = "Assay1",
        tableSignificance = "up",
        daQValThresh = 0.05,
        lfcThreshold = 0.5,
        tableMaxRows = 1,
        datatableFactory = function(data, options, extensions, rownames, filter) {
            datatable_calls <<- list(
                data = data,
                options = options,
                extensions = extensions,
                rownames = rownames,
                filter = filter
            )

            structure(
                list(data = data, roundCalls = list(), styleCall = NULL),
                class = "mock_dt"
            )
        },
        formatRoundFactory = function(widget, columns, digits) {
            widget$roundCalls <- c(
                widget$roundCalls,
                list(list(columns = columns, digits = digits))
            )
            widget
        },
        formatStyleFactory = function(widget, columns, backgroundColor) {
            widget$styleCall <- list(columns = columns, backgroundColor = backgroundColor)
            widget
        },
        styleEqualFactory = function(levels, values) {
            style_equal_calls <<- list(levels = levels, values = values)
            list(levels = levels, values = values)
        }
    )

    expect_s3_class(table_widget, "mock_dt")
    expect_identical(names(datatable_calls$data), c(
        "lipid_id", "lipid_name", "assay", "logFC",
        "raw_pvalue", "fdr_qvalue", "significant"
    ))
    expect_identical(nrow(datatable_calls$data), 1L)
    expect_identical(datatable_calls$data$lipid_id, "L1")
    expect_identical(datatable_calls$options$pageLength, 25)
    expect_true(isTRUE(datatable_calls$options$scrollX))
    expect_identical(datatable_calls$options$dom, "Bfrtip")
    expect_identical(datatable_calls$options$buttons, c("copy", "csv", "excel"))
    expect_identical(datatable_calls$extensions, "Buttons")
    expect_identical(datatable_calls$rownames, FALSE)
    expect_identical(datatable_calls$filter, "top")
    expect_length(table_widget$roundCalls, 2)
    expect_identical(
        table_widget$roundCalls[[1]],
        list(columns = "logFC", digits = 3)
    )
    expect_identical(
        table_widget$roundCalls[[2]],
        list(columns = c("raw_pvalue", "fdr_qvalue"), digits = 6)
    )
    expect_identical(table_widget$styleCall$columns, "significant")
    expect_identical(style_equal_calls$levels, c("Up", "Down", "NS"))
    expect_identical(style_equal_calls$values, c("#ffcccc", "#cce5ff", "white"))
    expect_identical(table_widget$styleCall$backgroundColor, style_equal_calls)
})

test_that("buildLipidDaResultsDownloadHandler keeps the DA results export contract stable", {
    handler_calls <- list()
    csv_call <- NULL

    handler <- buildLipidDaResultsDownloadHandler(
        daResultsReactive = function() {
            data.frame(
                lipid_id = "L1",
                logFC = 1.2,
                stringsAsFactors = FALSE
            )
        },
        downloadHandlerFactory = function(filename, content) {
            handler_calls <<- list(filename = filename, content = content)
            structure(list(filename = filename, content = content), class = "mock_download_handler")
        },
        currentDate = function() as.Date("2026-04-15"),
        csvWriter = function(data, file, row.names) {
            csv_call <<- list(data = data, file = file, row.names = row.names)
        },
        reqFn = function(value) {
            expect_s3_class(value, "data.frame")
            value
        }
    )

    expect_s3_class(handler, "mock_download_handler")
    expect_identical(handler_calls$filename(), "lipidomics_da_results_2026-04-15.csv")

    handler_calls$content("/tmp/lipid-da-results.csv")

    expect_identical(csv_call$file, "/tmp/lipid-da-results.csv")
    expect_identical(csv_call$row.names, FALSE)
    expect_identical(csv_call$data$lipid_id, "L1")
    expect_identical(csv_call$data$logFC, 1.2)
})

test_that("resolveLipidDaResultsDownloadData keeps the results forwarding stable", {
    da_results_list <- list(
        da_lipids_long = data.frame(
            lipid_id = c("L1", "L2"),
            logFC = c(1.2, -0.7),
            stringsAsFactors = FALSE
        ),
        other_payload = "ignored"
    )

    expect_identical(
        resolveLipidDaResultsDownloadData(da_results_list),
        da_results_list$da_lipids_long
    )
    expect_null(resolveLipidDaResultsDownloadData(NULL))
})

test_that("buildLipidDaResultsDownloadOutputHandler keeps the output-shell handoff stable", {
    resolve_arg <- NULL
    reactive_value <- data.frame(
        lipid_id = "L9",
        logFC = 3.1,
        stringsAsFactors = FALSE
    )
    handler_calls <- list()
    da_results_list <- list(da_lipids_long = data.frame(lipid_id = "unused"))

    handler <- buildLipidDaResultsDownloadOutputHandler(
        daResultsList = da_results_list,
        buildResultsDownloadHandlerFn = function(daResultsReactive) {
            handler_calls <<- list(daResultsReactive = daResultsReactive)
            structure(list(daResultsReactive = daResultsReactive), class = "mock_download_handler_shell")
        },
        resolveResultsDownloadDataFn = function(value) {
            resolve_arg <<- value
            reactive_value
        }
    )

    expect_s3_class(handler, "mock_download_handler_shell")
    expect_identical(resolve_arg, NULL)
    expect_identical(handler_calls$daResultsReactive(), reactive_value)
    expect_identical(resolve_arg, da_results_list)
})

test_that("bootstrapLipidDaSaveHeatmap keeps the observer bootstrap handoff stable", {
    handle_call <- NULL
    heatmap_plot <- structure(list(kind = "heatmap"), class = "mock_heatmap_plot")
    row_clusters <- c(LipidA = 1L, LipidB = 2L)

    result <- bootstrapLipidDaSaveHeatmap(
        currentHeatmapPlot = heatmap_plot,
        currentRowClusters = row_clusters,
        experimentPaths = list(
            publication_graphs_dir = "/tmp/publication_graphs"
        ),
        heatmapContrast = "Contrast A/B",
        heatmapTopN = 75,
        heatmapClusterMethod = "complete",
        heatmapDistanceMethod = "pearson",
        heatmapClustering = "column",
        heatmapScaling = "both",
        heatmapColorScheme = "plasma",
        heatmapTreeCutMethod = "dynamic",
        heatmapNClusters = 6,
        heatmapCutHeight = 12,
        heatmapMinClusterSize = 7,
        handleSaveHeatmapFn = function(...) {
            handle_call <<- list(...)
            list(marker = "save-result")
        }
    )

    expect_identical(handle_call$currentHeatmapPlot, heatmap_plot)
    expect_identical(handle_call$currentRowClusters, row_clusters)
    expect_identical(handle_call$publicationGraphsDir, "/tmp/publication_graphs")
    expect_identical(handle_call$heatmapContrast, "Contrast A/B")
    expect_identical(handle_call$heatmapTopN, 75)
    expect_identical(handle_call$heatmapClusterMethod, "complete")
    expect_identical(handle_call$heatmapDistanceMethod, "pearson")
    expect_identical(handle_call$heatmapClustering, "column")
    expect_identical(handle_call$heatmapScaling, "both")
    expect_identical(handle_call$heatmapColorScheme, "plasma")
    expect_identical(handle_call$heatmapTreeCutMethod, "dynamic")
    expect_identical(handle_call$heatmapNClusters, 6)
    expect_identical(handle_call$heatmapCutHeight, 12)
    expect_identical(handle_call$heatmapMinClusterSize, 7)
    expect_identical(result, list(marker = "save-result"))
})

test_that("handleLipidDaSaveHeatmap keeps the save-heatmap observer contract stable", {
    save_calls <- list()
    logged_message <- NULL
    notification <- NULL
    heatmap_plot <- structure(list(kind = "heatmap"), class = "mock_heatmap_plot")
    row_clusters <- c(LipidA = 1L, LipidB = 2L)

    result <- handleLipidDaSaveHeatmap(
        currentHeatmapPlot = heatmap_plot,
        currentRowClusters = row_clusters,
        publicationGraphsDir = "/tmp/publication_graphs",
        heatmapContrast = "Contrast A/B",
        heatmapTopN = 75,
        heatmapClusterMethod = "complete",
        heatmapDistanceMethod = "pearson",
        heatmapClustering = "column",
        heatmapScaling = "both",
        heatmapColorScheme = "plasma",
        heatmapTreeCutMethod = "dynamic",
        heatmapNClusters = 6,
        heatmapCutHeight = 12,
        heatmapMinClusterSize = 7,
        saveHeatmapProducts = function(...) {
            save_calls <<- list(...)
        },
        logInfo = function(message) {
            logged_message <<- message
        },
        notify = function(message, type, duration) {
            notification <<- list(message = message, type = type, duration = duration)
        }
    )

    expect_identical(logged_message, "Save Heatmap button clicked")
    expect_identical(save_calls$heatmap_obj, heatmap_plot)
    expect_identical(save_calls$row_clusters, row_clusters)
    expect_identical(save_calls$output_dir, "/tmp/publication_graphs")
    expect_identical(save_calls$file_prefix, "lipid_Contrast_A_B")
    expect_identical(save_calls$params_list$contrast, "Contrast A/B")
    expect_identical(save_calls$params_list$top_n, 75)
    expect_identical(save_calls$params_list$clustering_method, "complete")
    expect_identical(save_calls$params_list$distance_method, "pearson")
    expect_identical(save_calls$params_list$cluster_rows, FALSE)
    expect_identical(save_calls$params_list$cluster_cols, TRUE)
    expect_identical(save_calls$params_list$scaling, "both")
    expect_identical(save_calls$params_list$color_scheme, "plasma")
    expect_identical(save_calls$params_list$tree_cut_method, "dynamic")
    expect_identical(save_calls$params_list$n_clusters, 6)
    expect_identical(save_calls$params_list$cut_height, 12)
    expect_identical(save_calls$params_list$min_cluster_size, 7)
    expect_identical(notification$message, "Heatmap and cluster info saved to publication_graphs/Heatmap")
    expect_identical(notification$type, "message")
    expect_identical(notification$duration, 5)
    expect_identical(result$prefix, "lipid_Contrast_A_B")
    expect_identical(result$params, save_calls$params_list)
})
