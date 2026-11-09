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

# KIT-097: the transcriptomics Summary & Export tab, the lipid summary module under the transcriptomics profile. The
# profile adds the RUVIII export (KIT-095's `trx_ruviii_export`, a product after DA, accepted on the head; the DA head and
# results are unchanged) with its M choice: the group (the notebook) or the design's technical replicates, offered only
# when the design has them. The integration export is that RUVIII matrix (KIT-095's writer); the template is the
# package's `reports/transcriptomics/transcriptomics_report.rmd`.


#' The export panel above the lipid summary's panels.
#' @noRd
trxSummaryExportUi <- function(ns) {
    shiny::fluidRow(shiny::column(12, shiny::wellPanel(
        shiny::h4("RUVIII Integration Export"),
        shiny::helpText("RUVIII (ruv::RUVIII) on the TMM logCPM, with the RUVs inferred negative controls and k, for the",
            "integration export only: the DA results are not changed. M is the replicate structure: the group (the notebook),",
            "or the design's technical replicates when it has them."),
        shiny::uiOutput(ns("ruviii_m_ui")),
        testid(shiny::actionButton(ns("run_ruviii_export"), "Run RUVIII Export", class = "btn-primary", icon = shiny::icon("play")),
            "trx-summary-ruviii-export"),
        shiny::br(), shiny::br(), shiny::verbatimTextOutput(ns("ruviii_status")),
        shinyjqui::jqui_resizable(shiny::imageOutput(ns("ruviii_before_after"), height = "500px")))))
}


#' Run the export on the worker (the product, its acceptance on the head, its figures) or in process (the wrapper).
#' @noRd
runTrxRuviiiExport <- function(m_source, workflowData, experimentPaths, shown) {
    manager <- workflowData$state_manager
    failed <- function(error) {
        message <- paste("RUVIII export error:", trxConditionText(error))
        logger::log_error(message)
        shiny::showNotification(message, type = "error", duration = 15)
    }
    figures <- list(list("trx_qc_figures", list(stage = "ruviii")), list("trx_composite_figure", list(figure = "ruviii_before_and_after")))
    runner <- workflowOmicsManagedActions("transcriptomics")
    if (!is.null(runner)) {
        held <- new.env(parent = emptyenv())
        source <- getWorkflowStateReference(manager)
        ruvs <- withWorkflowStateScope(manager, function() trxStateRuvs(manager, source))
        export <- omicsManagedProduct("trx_ruviii_export", list(logcpm = ruvs$logcpm$entity_id, grid = ruvs$grid$entity_id, ruviii_m_source = m_source),
            function(result) {
                acceptTrxRuviiiExport(manager, source, list(reference = result$reference, value = result$value))
                held$value <- result$value
            })
        actions <- c(list(export), lapply(figures, function(f) omicsManagedProduct(f[[1L]], f[[2L]], function(result) invisible(NULL))))
        runner$run(actions, done = function() shown(list(text = trxRuviiiStatusText(held$value), at = Sys.time())), failed = failed)
        return(invisible(TRUE))
    }
    tryCatch({
        run <- runTrxReferenceRuviiiExport(manager, m_source)
        for (f in figures) runWorkflowOmicsPublicationInProcess(transcriptomicsDeclarationBinding, manager, f[[1L]], f[[2L]], experimentPaths)
        shown(list(text = trxRuviiiStatusText(run$value), at = Sys.time()))
    }, error = failed)
}

#' The export panel's server (the profile's `summary_server`).
#' @noRd
trxSummaryServer <- function(input, output, session, workflowData, experimentPaths) {
    shown <- shiny::reactiveVal(NULL)
    output$ruviii_m_ui <- shiny::renderUI({
        workflowData$tab_status
        replicates <- trxHeadHasTechnicalReplicates(workflowData$state_manager)
        shiny::tagList(shiny::radioButtons(session$ns("ruviii_m_source"), "Replicate structure M:", choices = trxRuviiiMChoices(replicates),
            selected = "group"), if (!replicates) shiny::helpText("Technical replicates: not available; the design has no technical",
            "replicates (trx_ruviii_replicates_invalid)."))
    })
    shiny::observeEvent(input$run_ruviii_export, {
        shiny::req(workflowData$state_manager)
        runTrxRuviiiExport(input$ruviii_m_source %||% "group", workflowData, experimentPaths, shown)
    })
    output$ruviii_status <- shiny::renderText(if (is.null(shown())) "Not run yet." else shown()$text)
    output$ruviii_before_after <- shiny::renderImage({
        shown()
        path <- file.path(experimentPaths$publication_graphs_dir %||% "", "ruviii_before_and_after.png")
        if (!is.null(shown()) && file.exists(path)) list(src = path, contentType = "image/png", width = "100%", height = "auto",
            alt = "RUVIII before and after") else list(src = "", alt = "RUVIII export not run yet")
    }, deleteFile = FALSE)
    invisible(shown)
}

