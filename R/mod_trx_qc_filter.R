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

# KIT-096: the transcriptomics QC tab's expression-filter sub-tab. KIT-104 W6: it is the `trx_qc_filter` plugin of the
# shared QC shell (mod_omic_qc.R), beside the shared finalise; its figure and text functions are func_trx_qc_plots.R's.
# The filter is KIT-095's `trx_expression_filter` (filterByExpr on the TMM-scaled library sizes of the design's groups,
# then the zero-variance genes); the tab passes no library size of its own. Its metrics
# (genes before, after filterByExpr and after the zero-variance step; library size and genes detected per run) and the
# raw-logCPM PCA and RLE come from the worker on a managed session, or from the state read for one action in process.

# ---------------------------------------------------------------------------------------------------------------
# The sub-tab.
# ---------------------------------------------------------------------------------------------------------------

#' @title Transcriptomics Expression Filter Sub-tab
#' @param id Sub-module id (`filter`).
#' @name mod_trx_qc_filter
#' @noRd
mod_trx_qc_filter_ui <- function(id) {
    ns <- shiny::NS(id)
    shiny::tabPanel("Expression Filter", shiny::br(), shiny::fluidRow(
        shiny::column(4, shiny::wellPanel(
            shiny::h4("Gene Expression Filter"),
            shiny::p("Genes are kept by edgeR's filterByExpr on the TMM-scaled library sizes of the design's groups; genes",
                "with zero variance across the samples are then removed."),
            shiny::helpText("filterByExpr's defaults: min.count 10, min.total.count 15, large.n 10, min.prop 0.7."),
            shiny::hr(),
            shiny::div(testid(shiny::actionButton(ns("apply_filter"), "Apply Filter", class = "btn-primary", width = "48%"), "trx-qc-apply-filter"),
                shiny::actionButton(ns("revert_filter"), "Revert", class = "btn-warning", width = "48%", style = "margin-left: 4%")))),
        shiny::column(8, shiny::verbatimTextOutput(ns("filter_results")), DT::DTOutput(ns("sample_metrics")), shiny::br(),
            shinyjqui::jqui_resizable(shiny::plotOutput(ns("filter_plot"), height = "600px", width = "100%")),
            shiny::fluidRow(shiny::column(6, shiny::plotOutput(ns("raw_pca"), height = "400px")),
                shiny::column(6, shiny::plotOutput(ns("raw_rle"), height = "400px"))))))
}

#' Apply the filter: on the worker under a managed session (the commit, then the raw figures), or through the wrapper.
#' @param shown The sub-tab's `reactiveVal` of what it displays.
#' @noRd
runTrxQcExpressionFilter <- function(workflowData, shown, output) {
    manager <- workflowData$state_manager
    finish <- function(counts, metrics, stats) {
        shown(list(counts = counts, metrics = metrics, stats = stats))
        output$filter_results <- shiny::renderText(trxQcFilterText(counts))
        shiny::showNotification(sprintf("Expression filter applied: %d genes retained", counts$after_zero_variance), type = "message")
    }
    tryCatch({
        if (!identical(readWorkflowOmicsMetadata(transcriptomicsDeclarationBinding, manager)$class, "TranscriptomicsAssayData")) {
            stop("Current state is not a TranscriptomicsAssayData object")
        }
        runner <- workflowOmicsManagedActions("transcriptomics")
        if (!is.null(runner)) {
            held <- new.env(parent = emptyenv())
            runner$run(list(omicsManagedCommit(trxExpressionFilterStep(), workflowData$config_list, function(result) {
                held$counts <- result$value$counts
                held$metrics <- readWorkflowOmicsFilteringMetrics(manager, result)
            }), omicsManagedProduct("trx_qc_figures", list(stage = "raw"), function(result) {
                if (!is.null(result$display$stats)) held$stats <- resolveWorkflowStateValue(manager$.context, result$display$stats)
            })), done = function() finish(held$counts, held$metrics, held$stats))
            return(invisible(NULL))
        }
        run <- runTrxReferenceExpressionFilter(manager, workflowData$config_list)
        # KIT-075: the state is read for this action only; the display keeps numbers.
        object <- readTrxWorkflowState(manager, run$reference)
        metrics <- trxFilteringMetrics(object)
        stats <- trxQcRawStats(object)
        rm(object)
        finish(run$value$counts, metrics, stats)
    }, error = function(error) {
        message <- paste("Error applying the expression filter:", trxConditionText(error))
        logger::log_error(message)
        shiny::showNotification(message, type = "error", duration = 15)
    })
}

#' @rdname mod_trx_qc_filter
#' @param workflow_data,omic_type,experiment_label The QC tab's arguments.
#' @return The reactive of what the sub-tab displays (counts, metrics, raw statistics).
#' @noRd
mod_trx_qc_filter_server <- function(id, workflow_data, omic_type, experiment_label) {
    shiny::moduleServer(id, function(input, output, session) {
        shown <- shiny::reactiveVal(NULL)
        # KIT-101 F-17: the server starts after the design; a click made before it started is not run then.
        clickAfterStart <- qcClickAfterStart(input, "apply_filter")
        shiny::observeEvent(input$apply_filter, {
            if (!clickAfterStart()) return(invisible(NULL))
            shiny::req(workflow_data$state_manager)
            # KIT-101 F-11: before the import the session's owner is still the pending placeholder; refuse saying so.
            if (isDIASessionOwnerPending(workflow_data$state_manager)) {
                shiny::showNotification("Please import the count data before applying the expression filter.", type = "warning")
                return(invisible(NULL))
            }
            runTrxQcExpressionFilter(workflow_data, shown, output)
        })
        # The shared revert (mod_omic_qc.R): the previous state selected, the display cleared.
        registerOmicQcRevertObserver(input, output, "revert_filter", workflow_data, "filter_results",
            function() shown(NULL), "Reverted lipid intensity filter to")
        output$sample_metrics <- DT::renderDT({
            shiny::req(shown()$metrics)
            DT::datatable(trxQcMetricsTable(shown()$metrics), options = list(dom = "tp", pageLength = 12), rownames = FALSE)
        })
        output$filter_plot <- shiny::renderPlot({
            shiny::req(shown()$metrics)
            print(trxQcFilteringFigure(shown()$counts, shown()$metrics))
        })
        output$raw_pca <- shiny::renderPlot({
            shiny::req(shown()$stats)
            print(trxPcaPlot(shown()$stats$pca, "PCA: raw logCPM"))
        })
        output$raw_rle <- shiny::renderPlot({
            shiny::req(shown()$stats)
            print(trxRlePlot(shown()$stats$rle, "RLE: raw logCPM"))
        })
        shown
    })
}
