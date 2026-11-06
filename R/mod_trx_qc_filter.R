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

# KIT-096: the transcriptomics QC tab's expression-filter sub-tab and the collaborators the lipid QC finalise takes under
# the transcriptomics profile. The filter is KIT-095's `trx_expression_filter` (filterByExpr on the TMM-scaled library
# sizes of the design's groups, then the zero-variance genes); the tab passes no library size of its own. Its metrics
# (genes before, after filterByExpr and after the zero-variance step; library size and genes detected per run) and the
# raw-logCPM PCA and RLE come from the worker on a managed session, or from the state read for one action in process.

# ---------------------------------------------------------------------------------------------------------------
# What the sub-tab shows.
# ---------------------------------------------------------------------------------------------------------------

#' The raw-logCPM PCA and RLE statistics of a state (the `trx_qc_figures` operation's stage `raw`, in process).
#' @noRd
trxQcRawStats <- function(object) {
    x <- trxAssayInputs(object)
    design <- trxFigureDesign(x)
    raw <- trxTmmLogCpm(x$counts, x$design$group)$log_cpm_raw
    list(pca = trxPcaStats(raw, design), rle = trxRleStats(raw, design))
}

#' The filter's result text: the genes at each step, the zero-variance step named.
#' @noRd
trxQcFilterText <- function(counts) {
    paste("Expression Filter Applied", "=========================", "filterByExpr on TMM-scaled library sizes, then zero-variance genes removed",
        "", sprintf("Genes before filtering: %d", counts$before), sprintf("After filterByExpr: %d", counts$after_filter_by_expr),
        sprintf("After the Zero-variance step: %d", counts$after_zero_variance), "", "State saved as: 'trx_expression_filtered'", sep = "\n")
}

#' The per-run metrics table: library size and genes detected.
#' @noRd
trxQcMetricsTable <- function(metrics) {
    data.frame(Run = names(metrics$library_size), `Library size` = unname(metrics$library_size),
        `Genes detected` = unname(metrics$detected_genes), check.names = FALSE, stringsAsFactors = FALSE)
}

#' The filtering progress figure from the counts and the metrics: the genes at each step (when given), and per run the
#' library size and the genes detected.
#' @return A ggplot.
#' @noRd
trxQcFilteringFigure <- function(counts, metrics) {
    panel <- function(name, x, value) data.frame(panel = name, x = x, value = as.numeric(value), stringsAsFactors = FALSE)
    rows <- rbind(if (!is.null(counts)) panel("Genes retained by step", c("Before", "filterByExpr", "Zero variance"), unlist(counts)),
        panel("Library size per run", names(metrics$library_size), metrics$library_size),
        panel("Genes detected per run", names(metrics$detected_genes), metrics$detected_genes))
    rows$x <- factor(rows$x, levels = unique(rows$x))
    rows$panel <- factor(rows$panel, levels = unique(rows$panel))
    ggplot2::ggplot(rows, ggplot2::aes(x = .data$x, y = .data$value)) + ggplot2::geom_col(fill = "steelblue") +
        ggplot2::facet_wrap(ggplot2::vars(.data$panel), scales = "free", ncol = 1) + ggplot2::labs(x = NULL, y = NULL) +
        ggplot2::theme_bw() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

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
        shiny::observeEvent(input$apply_filter, {
            shiny::req(workflow_data$state_manager)
            # KIT-101 F-11: before the import the session's owner is still the pending placeholder; refuse saying so.
            if (isDIASessionOwnerPending(workflow_data$state_manager)) {
                shiny::showNotification("Please import the count data before applying the expression filter.", type = "warning")
                return(invisible(NULL))
            }
            runTrxQcExpressionFilter(workflow_data, shown, output)
        })
        # The lipid intensity sub-tab's revert: the previous state selected, the display cleared.
        registerLipidIntensityRevertObserver(input, workflow_data, output, filterStatsVal = shown, filterPlotVal = shown)
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

# ---------------------------------------------------------------------------------------------------------------
# The QC tab under the profile: its sub-tabs, and the finalise's collaborators.
# ---------------------------------------------------------------------------------------------------------------

#' The transcriptomics QC sub-tabs: the expression filter and the lipid finalise under the profile.
#' @noRd
trxQcTabsUi <- function(ns) {
    shiny::tabsetPanel(id = ns("trx_qc_tabs"), mod_trx_qc_filter_ui(ns("filter")),
        mod_lipid_qc_s4_ui(ns("s4_finalize"), profile = trxOmicProfile()))
}

#' Start the transcriptomics QC sub-modules (the lipid QC orchestrator's `initializeLipidQcSubmodules` seam).
#' @noRd
initializeTrxQcSubmodules <- function(workflow_data, omic_type, experiment_label) {
    mod_trx_qc_filter_server("filter", workflow_data, omic_type, experiment_label)
    mod_lipid_qc_s4_server("s4_finalize", workflow_data, omic_type, experiment_label, profile = trxOmicProfile())
}

#' The finalise's data summary from the state's metadata: the stored row count and the design; no count is read.
#' @noRd
registerTrxQcSummaryOutput <- function(output, workflowData) {
    output$data_summary <- shiny::renderUI({
        shiny::req(workflowData$state_manager)
        meta <- tryCatch(readWorkflowOmicsMetadata(transcriptomicsDeclarationBinding, workflowData$state_manager), error = function(e) NULL)
        if (is.null(meta) || !identical(meta$class, "TranscriptomicsAssayData")) {
            return(shiny::p(shiny::icon("exclamation-triangle", style = "color: orange;"), " No TranscriptomicsAssayData object available."))
        }
        row <- function(label, value) shiny::tags$tr(shiny::tags$td(shiny::strong(label)), shiny::tags$td(value))
        design <- meta$design_matrix
        shiny::tags$table(class = "table table-condensed", style = "margin-bottom: 0;", shiny::tags$tbody(
            row("Total Genes:", format(sum(meta$assay_rows), big.mark = ",")), row("Number of Samples:", nrow(design)),
            row("Experimental Groups:", length(unique(design[[meta$slots$group_id]]))), row("Gene ID Column:", meta$slots$gene_id_column),
            row("Sample ID Column:", meta$slots$sample_id)))
    })
    output
}

#' The finalise's results text (genes retained from the stored row count).
#' @noRd
trxQcFinalizeText <- function(currentS4, history) {
    paste("QC Finalization Complete", "========================", "", sprintf("Total genes retained: %d", sum(currentS4$assay_rows)),
        sprintf("Processing steps completed: %d", length(history)), "", "Processing History:",
        paste(sprintf("  %d. %s", seq_along(history), history), collapse = "\n"), "", "State saved as: 'trx_qc_complete'", "",
        "You can now proceed to the Normalization tab.", sep = "\n")
}

#' The collaborators `runLipidQcS4FinalizeWorkflow` takes under the transcriptomics profile: the state read as metadata, its
#' checks, the store step `trx_qc_finalize`, the filtering figure from the metrics and the results text.
#' @noRd
trxQcFinalizeCollaborators <- function() {
    list(
        getFinalizeStateFn = function(stateManager) readWorkflowOmicsMetadata(transcriptomicsDeclarationBinding, stateManager),
        validateFinalizeStateFn = function(currentS4) {
            if (is.null(currentS4) || !identical(currentS4$class, "TranscriptomicsAssayData")) stop("Current state is not a TranscriptomicsAssayData object")
            if (!length(currentS4$assay_rows) || any(currentS4$assay_rows == 0L) || !nrow(currentS4$design_matrix)) {
                stop("Current TranscriptomicsAssayData object has no genes or no design rows to finalize.")
            }
            currentS4
        },
        saveCompletedStateFn = function(stateManager, currentS4, configObject) {
            commitWorkflowOmicsStoreStep(stateManager, "trx_qc_finalize", "trx_qc_complete", configObject, "QC processing complete - ready for normalization")
        },
        updateTrackingPlotFn = function(currentS4, omicType, setFilterPlotFn, metrics = NULL) {
            plot <- tryCatch(trxQcFilteringFigure(NULL, metrics %||% trxFilteringMetrics(currentS4)), error = function(e) NULL)
            setFilterPlotFn(plot)
        },
        reportFinalizeSuccessFn = function(currentS4, history, output) {
            reportLipidQcS4FinalizeSuccess(currentS4, history, output, buildResultsTextFn = trxQcFinalizeText,
                logInfoFn = function(...) logger::log_info("Transcriptomics QC finalized successfully"))
        })
}
