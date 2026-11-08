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

# KIT-097: the transcriptomics DA tab, the lipid DA module under the transcriptomics profile: its outputs, observers and
# helpers with the transcriptomics collaborators. The analysis is the one KIT-095 records: edgeR glmQLFit(robust = TRUE),
# glmQLFTest, Benjamini-Hochberg FDR (shown as "FDR (BH)", never a q-value), the realised formula with the RUVs factors
# W_1..W_k, no TREAT; the fold-change field is a display threshold only. One pairwise result per contrast on one fit.

#' The model the DA tab shows in place of the lipid formula input.
#' @noRd
trxDaModelUi <- function(ns) {
    shiny::tagList(shiny::h5("Model"), shiny::verbatimTextOutput(ns("formula_display")),
        shiny::helpText("edgeR quasi-likelihood: estimateDisp and glmQLFit(robust = TRUE), then glmQLFTest per contrast on the",
            "one fit; Benjamini-Hochberg FDR. The formula is realised from the design: ~ 0 + group plus the RUVs factors",
            "W_1..W_k as covariates. No fold-change test is applied."))
}

#' The tab's status text: the lipid text once the analysis has run.
#' @noRd
trxDaStatusText <- function(analysisComplete, daResultsList, currentS4Object) {
    if (isTRUE(analysisComplete)) return(buildLipidDaStatusText(analysisComplete, daResultsList, currentS4Object))
    "Ready. The DA runs on the normalised (RUVs) state."
}

#' The interactive volcano's widget: the binding's Glimma widget in process, the worker's page when managed; the contrast
#' name mapped to its entry.
#' @noRd
trxDaWidget <- function(workflowData, runner = NULL) {
    widget <- if (is.null(runner)) trxVolcanoWidget else
        workflowOmicsVolcanoPageFactory(transcriptomicsDeclarationBinding, runner, workflowData$state_manager)
    function(da_results_list, selected_contrast, selected_assay, da_q_val_thresh) {
        widget(da_results_list, trxDaContrastEntry(workflowData, selected_contrast), selected_assay, da_q_val_thresh)
    }
}

#' The DA results with the head state read for one in-process display (the Glimma widget's annotation).
#' @noRd
trxDaWithSource <- function(results, workflowData) {
    if (is.null(results)) return(results)
    manager <- workflowData$state_manager
    results$theObject <- withWorkflowStateScope(manager, function() readTrxWorkflowState(manager))
    results
}

#' The in-process DA publication (the worker's `trx_da_publication`, run here).
#' @noRd
trxDaWriteResults <- function(workflowData, experimentPaths) {
    function(...) {
        display <- runWorkflowOmicsPublicationInProcess(transcriptomicsDeclarationBinding, workflowData$state_manager, "trx_da_publication", list(), experimentPaths)
        if (!isTRUE(display$saved)) stop(display$message %||% "the DA publication reported no success", call. = FALSE)
        TRUE
    }
}

#' The run button's analysis: the lipid bootstrap with the transcriptomics collaborators.
#' @noRd
runTrxDaAnalysis <- function(input, daData, workflowData, session, experimentPaths, profile) {
    finalize <- function(results, ...) {
        # The lipid selectors and filters match a contrast by `comparison` or `friendly_name`; here they are the same name.
        if (!is.null(results$da_genes_long)) results$da_genes_long$friendly_name <- results$da_genes_long$comparison
        finalizeLipidDaRunAnalysisSuccess(results = results, ..., longTable = profile$long_table, singleAssay = TRUE,
            writeResultsFn = trxDaWriteResults(workflowData, experimentPaths), withSourceFn = function(results, ...) results)
    }
    bootstrapLipidDaRunAnalysis(daData = daData, workflowData = workflowData, formulaString = NULL, daQValThresh = input$da_q_val_thresh,
        treatLfcCutoff = input$treat_lfc_cutoff, session = session, experimentPaths = experimentPaths,
        prepareAnalysisContextFn = function(daData, workflowData) prepareLipidDaRunAnalysisContext(daData, workflowData,
            resolveCurrentStateFn = trxDaReadMetadata, expectedClass = "TranscriptomicsAssayData",
            noDataMessage = "Run the normalisation first: the DA runs on the normalised (RUVs) state."),
        handlePreflightFn = function(...) handleLipidDaRunAnalysisPreflight(..., executeRunAnalysisFn = function(...) {
            executeLipidDaRunAnalysis(..., successFinalizerFn = finalize, profile = profile)
        }))
}

#' @title Transcriptomics DA Tab
#' @description The lipid DA module's server under the transcriptomics profile (`mod_lipid_da_server(profile =)`).
#' @return The tab's reactive values.
#' @name mod_trx_da
#' @noRd
mod_trx_da_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
    shiny::moduleServer(id, function(input, output, session) {
        profile <- trxOmicProfile()
        da_data <- initializeLipidDaServerState()
        registerLipidDaPrimaryTextOutputs(output, workflow_data, da_data, buildStatusTextFn = trxDaStatusText)
        output$formula_display <- shiny::renderText({
            workflow_data$tab_status
            trxDaModelText(workflow_data)
        })
        shiny::observeEvent(input$run_da_analysis, runTrxDaAnalysis(input, da_data, workflow_data, session, experiment_paths, profile))
        registerLipidDaVolcanoGlimmaOutput(output, input, da_data, workflow_data, profile = profile)
        registerLipidDaVolcanoStaticOutput(output, input, da_data,
            buildVolcanoStaticPlotFn = function(...) buildLipidDaVolcanoStaticPlot(..., plotFactory = trxDaStaticVolcano))
        long <- profile$long_table
        registerLipidDaSummaryStatsOutput(output, input, da_data, buildSummaryStatsRenderTextFn = function(...) {
            buildLipidDaSummaryStatsRenderText(..., longTable = long, buildSummaryStatsTextFn = function(...) {
                buildLipidDaSummaryStatsText(..., featureNoun = "genes", significanceLabel = "FDR (BH)")
            })
        })
        registerLipidDaResultsTableOutput(output, input, da_data, buildResultsTableRenderWidgetFn = function(...) {
            buildLipidDaResultsTableRenderWidget(..., longTable = long,
                buildResultsTableWidgetFn = function(...) buildLipidDaResultsTableWidget(..., columns = profile$da_columns))
        })
        registerLipidDaResultsDownloadOutput(output, da_data, buildResultsDownloadOutputHandlerFn = function(readResultsFn) {
            buildLipidDaResultsDownloadOutputHandler(readResultsFn = readResultsFn,
                resolveResultsDownloadDataFn = function(results) resolveLipidDaResultsDownloadData(results, long),
                buildResultsDownloadHandlerFn = function(daResultsReactive) buildLipidDaResultsDownloadHandler(daResultsReactive, filePrefix = "transcriptomics"))
        })
        da_data
    })
}
