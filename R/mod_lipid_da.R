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
# mod_lipid_da.R
# ============================================================================
# Purpose: Lipidomics differential abundance analysis Shiny module
#
# This module provides per-assay differential abundance analysis using limma,
# with interactive volcano plots (Glimma), heatmaps, and result tables.
# UI and functionality matches the proteomics DA module (mod_prot_da.R).
# ============================================================================

#' @title Lipidomics Differential Analysis Module
#' @description A Shiny module for performing per-assay differential abundance
#'              analysis on lipidomics data using limma.
#' @name mod_lipid_da
NULL


initializeLipidDaServerState <- function(
    reactiveValuesFactory = shiny::reactiveValues
) {
    reactiveValuesFactory(
        da_results_list = NULL,
        contrasts_available = NULL,
        assays_available = NULL,
        analysis_complete = FALSE,
        current_s4_object = NULL,
        formula_from_s4 = NULL,
        current_row_clusters = NULL,
        current_col_clusters = NULL,
        current_heatmap_plot = NULL
    )
}

registerLipidDaPrimaryTextOutputs <- function(
    output,
    workflowData,
    daData,
    renderPrintFn = shiny::renderPrint,
    buildContrastsTextFn = buildLipidDaContrastsText,
    buildStatusTextFn = buildLipidDaStatusText,
    catFn = cat
) {
    output$contrasts_display <- renderPrintFn({
        catFn(buildContrastsTextFn(workflowData$contrasts_tbl))
    })

    output$da_status <- renderPrintFn({
        catFn(buildStatusTextFn(
            analysisComplete = daData$analysis_complete,
            daResultsList = daData$da_results_list,
            currentS4Object = daData$current_s4_object
        ))
    })

    invisible(output)
}

registerLipidDaHeatmapWarningOutput <- function(
    output,
    daData,
    renderUiFn = shiny::renderUI,
    buildWarningFn = buildLipidDaHeatmapManualSaveWarning
) {
    output$heatmap_manual_save_warning <- renderUiFn({
        buildWarningFn(
            analysisComplete = daData$analysis_complete
        )
    })

    invisible(output)
}

registerLipidDaVolcanoGlimmaOutput <- function(
    output,
    input,
    daData,
    renderUiFn = shiny::renderUI,
    buildVolcanoGlimmaUiFn = buildLipidDaVolcanoGlimmaUi
) {
    output$volcano_glimma <- renderUiFn({
        buildVolcanoGlimmaUiFn(
            daResultsList = daData$da_results_list,
            selectedContrast = input$volcano_contrast,
            selectedAssay = input$volcano_assay,
            daQValThresh = input$da_q_val_thresh
        )
    })

    invisible(output)
}

registerLipidDaVolcanoStaticOutput <- function(
    output,
    input,
    daData,
    renderPlotFn = shiny::renderPlot,
    buildVolcanoStaticPlotFn = buildLipidDaVolcanoStaticPlot
) {
    output$volcano_static <- renderPlotFn({
        buildVolcanoStaticPlotFn(
            daResultsList = daData$da_results_list,
            selectedContrast = input$volcano_contrast,
            selectedAssay = input$volcano_assay,
            daQValThresh = input$da_q_val_thresh,
            lfcThreshold = input$treat_lfc_cutoff
        )
    })

    invisible(output)
}

registerLipidDaHeatmapPlotOutput <- function(
    output,
    input,
    daData,
    renderPlotFn = shiny::renderPlot,
    buildHeatmapRenderStateFn = buildLipidDaHeatmapRenderState,
    storeHeatmapRenderStateFn = storeLipidDaHeatmapRenderState
) {
    output$heatmap_plot <- renderPlotFn({
        heatmapState <- buildHeatmapRenderStateFn(
            daResultsList = daData$da_results_list,
            selectedContrast = input$heatmap_contrast,
            selectedAssay = input$heatmap_assay,
            topN = input$heatmap_top_n,
            heatmapClustering = input$heatmap_clustering,
            scaleData = input$heatmap_scaling,
            clusteringMethod = input$heatmap_cluster_method,
            distanceMethod = input$heatmap_distance_method,
            colorScheme = input$heatmap_color_scheme,
            showLipidNames = input$heatmap_show_labels,
            daQValThresh = input$da_q_val_thresh,
            treeCutMethod = input$heatmap_tree_cut_method,
            nClusters = input$heatmap_n_clusters,
            cutHeight = input$heatmap_cut_height,
            minClusterSize = input$heatmap_min_cluster_size
        )

        storeHeatmapRenderStateFn(
            heatmapState = heatmapState,
            daData = daData
        )
    })

    invisible(output)
}

registerLipidDaClusterSummaryOutput <- function(
    output,
    input,
    daData,
    renderPrintFn = shiny::renderPrint,
    buildClusterSummaryRenderTextFn = buildLipidDaClusterSummaryRenderText,
    catFn = cat
) {
    output$cluster_summary <- renderPrintFn({
        catFn(buildClusterSummaryRenderTextFn(
            treeCutMethod = input$heatmap_tree_cut_method,
            rowClusters = daData$current_row_clusters
        ))
    })

    invisible(output)
}

registerLipidDaSaveHeatmapObserver <- function(
    input,
    daData,
    experimentPaths,
    bootstrapSaveHeatmapFn = bootstrapLipidDaSaveHeatmap,
    observeEventFn = shiny::observeEvent
) {
    observeEventFn(input$save_heatmap, {
        bootstrapSaveHeatmapFn(
            currentHeatmapPlot = daData$current_heatmap_plot,
            currentRowClusters = daData$current_row_clusters,
            experimentPaths = experimentPaths,
            heatmapContrast = input$heatmap_contrast,
            heatmapTopN = input$heatmap_top_n,
            heatmapClusterMethod = input$heatmap_cluster_method,
            heatmapDistanceMethod = input$heatmap_distance_method,
            heatmapClustering = input$heatmap_clustering,
            heatmapScaling = input$heatmap_scaling,
            heatmapColorScheme = input$heatmap_color_scheme,
            heatmapTreeCutMethod = input$heatmap_tree_cut_method,
            heatmapNClusters = input$heatmap_n_clusters,
            heatmapCutHeight = input$heatmap_cut_height,
            heatmapMinClusterSize = input$heatmap_min_cluster_size
        )
    })
}

registerLipidDaSummaryStatsOutput <- function(
    output,
    input,
    daData,
    renderPrintFn = shiny::renderPrint,
    buildSummaryStatsRenderTextFn = buildLipidDaSummaryStatsRenderText,
    catFn = cat
) {
    output$da_summary_stats <- renderPrintFn({
        catFn(buildSummaryStatsRenderTextFn(
            daResultsList = daData$da_results_list,
            selectedContrast = input$table_contrast,
            selectedAssay = input$table_assay,
            daQValThresh = input$da_q_val_thresh
        ))
    })

    invisible(output)
}

registerLipidDaResultsTableOutput <- function(
    output,
    input,
    daData,
    renderDtFn = DT::renderDT,
    buildResultsTableRenderWidgetFn = buildLipidDaResultsTableRenderWidget
) {
    output$da_results_table <- renderDtFn({
        buildResultsTableRenderWidgetFn(
            daResultsList = daData$da_results_list,
            selectedContrast = input$table_contrast,
            selectedAssay = input$table_assay,
            tableSignificance = input$table_significance,
            daQValThresh = input$da_q_val_thresh,
            lfcThreshold = input$treat_lfc_cutoff,
            tableMaxRows = input$table_max_rows
        )
    })

    invisible(output)
}

registerLipidDaResultsDownloadOutput <- function(
    output,
    daData,
    buildResultsDownloadOutputHandlerFn = buildLipidDaResultsDownloadOutputHandler
) {
    output$download_da_results <- buildResultsDownloadOutputHandlerFn(
        daResultsList = daData$da_results_list
    )

    invisible(output)
}








































