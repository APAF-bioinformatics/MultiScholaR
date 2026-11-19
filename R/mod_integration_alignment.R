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

# KIT-139: the integration workflow's Sample Matching tab. The cross-omic sample table the model layer builds
# (`integrationAlignViews`, through `integrationAlignmentPreview`) on the imported data sets: one row per biological
# sample matched by the biological-sample column each omic saved, which data sets measure it, a sample missing from a
# data set kept as missing (D4); what happens to technical replicates (INT-2) and internal standards (INT-3) stated as
# facts with counts from the parameters as resolved (never set here); the refusals in words. Matching works only
# through the user's confirmation of this table (KIT-060 D5): the confirmed table is saved as the integration design
# and kept on the bus with the parameters and manifests it was confirmed on, for the Model tab (KIT-140).

#' The Sample Matching tab's UI.
#' @param id Module id.
#' @noRd
mod_integration_alignment_ui <- function(id) {
    ns <- shiny::NS(id)
    texts <- integrationGuiTexts()
    shiny::fluidRow(shiny::column(12,
        shiny::h3(texts[["align_heading"]]),
        shiny::p(texts[["align_intro"]]),
        testid(shiny::uiOutput(ns("alignment_status")), "integration-alignment-status"),
        testid(shiny::uiOutput(ns("facts")), "integration-alignment-facts"),
        shiny::h4(texts[["view_table_heading"]]),
        testid(DT::DTOutput(ns("view_table")), "integration-alignment-views"),
        shiny::h4(texts[["sample_table_heading"]]),
        testid(DT::DTOutput(ns("sample_table")), "integration-alignment-samples"),
        testid(shiny::actionButton(ns("confirm"), texts[["confirm_button"]], class = "btn-primary",
            icon = shiny::icon("check")), "integration-alignment-confirm"),
        testid(shiny::uiOutput(ns("confirm_status")), "integration-alignment-confirmed")))
}

#' The Sample Matching tab's server (KIT-139).
#' @param id Module id.
#' @param workflow_data The integration bus (`state_manager`, `integration_import`, `integration_alignment`).
#' @param paths `project_dirs$integration` (`multiomic_inputs_dir`, else `results_dir`).
#' @param texts `integrationGuiTexts()`.
#' @return The preview (or its refusal), as a reactive.
#' @noRd
mod_integration_alignment_server <- function(id, workflow_data, paths, texts = integrationGuiTexts()) {
    shiny::moduleServer(id, function(input, output, session) {
        confirm_error <- shiny::reactiveVal(NULL)
        preview <- shiny::reactive({
            summary <- workflow_data$integration_import
            manager <- workflow_data$state_manager
            if (is.null(summary) || is.null(manager)) return(NULL)
            tryCatch(integrationAlignmentPreview(integrationImportedObject(manager), summary$views),
                integration_refused = function(error) {
                    structure(list(refusal = integrationAlignmentRefusal(error, summary$views, texts)),
                        class = "integration_alignment_refused")
                })
        })
        usable <- function(value) !is.null(value) && !inherits(value, "integration_alignment_refused")
        # A new import is matched again: a confirmation is for the import it was made on.
        shiny::observeEvent(workflow_data$integration_import, {
            workflow_data$integration_alignment <- NULL
            confirm_error(NULL)
        }, ignoreInit = TRUE)

        output$alignment_status <- shiny::renderUI({
            value <- preview()
            if (is.null(value)) return(shiny::div(class = "alert alert-info", texts[["align_not_imported"]]))
            if (!usable(value)) return(shiny::div(class = "alert alert-danger", value$refusal))
            keys <- lapply(seq_len(nrow(value$keys)), function(i) {
                shiny::tags$li(sprintf(texts[["align_key"]], tools::toTitleCase(value$keys$omic[[i]]), value$keys$column[[i]]))
            })
            shiny::tagList(shiny::h4(texts[["align_key_heading"]]), do.call(shiny::tags$ul, keys),
                shiny::p(texts[["align_caution"]]),
                if (isTRUE(value$few_samples)) shiny::div(class = "alert alert-warning",
                    sprintf(texts[["few_samples"]], value$aligned$shared_samples, INTEGRATION_FEW_SAMPLES)))
        })
        output$facts <- shiny::renderUI({
            value <- preview()
            if (!usable(value)) return(NULL)
            shiny::tagList(shiny::h4(texts[["align_rules_heading"]]),
                do.call(shiny::tags$ul, lapply(integrationAlignmentFacts(value, texts), shiny::tags$li)),
                shiny::p(shiny::em(texts[["align_provisional"]])))
        })
        output$view_table <- DT::renderDT({
            value <- preview()
            shiny::req(usable(value))
            DT::datatable(integrationAlignmentViewTable(value, texts), rownames = FALSE, selection = "none",
                options = list(dom = "t", paging = FALSE, ordering = FALSE))
        })
        output$sample_table <- DT::renderDT({
            value <- preview()
            shiny::req(usable(value))
            DT::datatable(integrationAlignmentSampleTable(value, texts), rownames = FALSE, selection = "none",
                options = list(pageLength = 25))
        })
        output$confirm_status <- shiny::renderUI({
            error <- confirm_error()
            if (!is.null(error)) return(shiny::div(class = "alert alert-danger", error))
            confirmed <- workflow_data$integration_alignment
            if (is.null(confirmed)) return(shiny::p(texts[["not_confirmed"]]))
            shiny::div(class = "alert alert-success", sprintf(texts[["confirmed"]],
                format(confirmed$confirmed_at, "%Y-%m-%d %H:%M"), confirmed$file))
        })

        shiny::observeEvent(input$confirm, {
            confirm_error(NULL)
            summary <- workflow_data$integration_import
            manager <- workflow_data$state_manager
            if (is.null(summary) || is.null(manager)) {
                confirm_error(texts[["align_not_imported"]])
                return()
            }
            # Matched again at the confirmation, so that a data set changed on disk since the tab drew is refused.
            value <- tryCatch(integrationAlignmentPreview(integrationImportedObject(manager), summary$views),
                integration_refused = function(error) integrationAlignmentRefusal(error, summary$views, texts))
            if (is.character(value)) {
                confirm_error(value)
                return()
            }
            dir <- paths$multiomic_inputs_dir %||% file.path(paths$results_dir, "multiomic_inputs")
            file <- integrationSaveDesign(value$aligned$samples, dir)
            workflow_data$integration_alignment <- list(samples = value$aligned$samples,
                view_summary = value$aligned$view_summary, parameters = value$parameters,
                manifests = summary$manifests, file = file, confirmed_at = Sys.time())
            tab_status <- workflow_data$tab_status
            tab_status$integration_alignment <- "complete"
            workflow_data$tab_status <- tab_status
        })
        preview
    })
}
