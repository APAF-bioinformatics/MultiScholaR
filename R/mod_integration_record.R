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

# KIT-142: the integration workflow's Record tab. It saves the integration's analysis record, its workbook, its study
# parameters and the report (`integrationWriteRecordBundle`) into the project's study report folder, and shows the
# draft methods paragraph read from the record. It opens once the results are saved (the shell's status chain), behind
# the integration gate like every integration tab.

#' The Record tab's UI.
#' @param id Module id.
#' @noRd
mod_integration_record_ui <- function(id) {
    ns <- shiny::NS(id)
    texts <- integrationRecordTexts()
    shiny::fluidRow(shiny::column(12,
        shiny::h3(texts[["record_heading"]]),
        shiny::p(texts[["record_intro"]]),
        testid(shiny::actionButton(ns("save"), texts[["record_button"]], class = "btn-primary", icon = shiny::icon("save")),
            "integration-record-save"),
        testid(shiny::uiOutput(ns("record_status")), "integration-record-status"),
        testid(shiny::uiOutput(ns("methods")), "integration-record-methods")))
}

#' The Record tab's server (KIT-142).
#' @param id Module id.
#' @param workflow_data The integration bus (`state_manager`, `integration_model`, `integration_record`, `tab_status`).
#' @param paths `project_dirs$integration` (`study_report_dir`, else `results_dir`).
#' @param texts `integrationGuiTexts()` (the refusals' words).
#' @param renderFn The report's renderer (`renderIntegrationReport`).
#' @param format The report's format.
#' @param withProgressFn The progress display around the save (`shiny::withProgress`).
#' @return The saved bundle (`integrationWriteRecordBundle`'s value) or NULL, as a reactive.
#' @noRd
mod_integration_record_server <- function(id, workflow_data, paths, texts = integrationGuiTexts(),
    renderFn = renderIntegrationReport, format = "docx", withProgressFn = shiny::withProgress) {
    shiny::moduleServer(id, function(input, output, session) {
        record_texts <- integrationRecordTexts()
        saved <- shiny::reactiveValues(bundle = NULL, error = NULL, saved_at = NULL)
        # A new model is a new analysis: its record is not saved yet.
        shiny::observeEvent(workflow_data$integration_model, {
            saved$bundle <- NULL
            saved$error <- NULL
            workflow_data$integration_record <- NULL
        }, ignoreNULL = FALSE, ignoreInit = TRUE)

        shiny::observeEvent(input$save, {
            saved$error <- NULL
            manager <- workflow_data$state_manager
            if (is.null(manager) || is.null(workflow_data$integration_model)) {
                saved$error <- list(sentence = record_texts[["record_not_ready"]], detail = NULL)
                return()
            }
            directory <- paths$study_report_dir %||% file.path(paths$results_dir, "Study_report")
            bundle <- tryCatch(withProgressFn(message = record_texts[["record_working"]], value = 0.5,
                integrationWriteRecordBundle(manager, directory, format = format, renderFn = renderFn)),
                error = function(error) {
                    saved$error <- integrationRefusalText(error, texts)
                    NULL
                })
            if (is.null(bundle)) return()
            saved$bundle <- bundle
            saved$saved_at <- Sys.time()
            workflow_data$integration_record <- list(directory = bundle$directory, files = bundle$files,
                methods = bundle$methods, saved_at = saved$saved_at)
            tab_status <- shiny::isolate(workflow_data$tab_status)
            tab_status$integration_record <- "complete"
            workflow_data$tab_status <- tab_status
        })

        output$record_status <- shiny::renderUI({
            if (!is.null(saved$error)) {
                return(integrationRefusalUi(sprintf(record_texts[["record_refused"]], saved$error$sentence), saved$error$detail,
                    texts = texts))
            }
            if (is.null(workflow_data$integration_model)) return(shiny::div(class = "alert alert-info", record_texts[["record_not_ready"]]))
            bundle <- saved$bundle
            if (is.null(bundle)) return(NULL)
            shiny::div(class = "alert alert-success",
                shiny::p(sprintf(record_texts[["record_done"]], format(saved$saved_at, "%Y-%m-%d %H:%M"), length(bundle$files),
                    bundle$directory)),
                do.call(shiny::tags$ul, lapply(basename(bundle$files), shiny::tags$li)))
        })
        output$methods <- shiny::renderUI({
            bundle <- saved$bundle
            if (is.null(bundle) || is.na(bundle$methods)) return(NULL)
            shiny::tagList(shiny::h4(record_texts[["record_methods_heading"]]), shiny::p(shiny::em(record_texts[["methods_draft_note"]])),
                shiny::p(bundle$methods))
        })
        shiny::reactive(saved$bundle)
    })
}
