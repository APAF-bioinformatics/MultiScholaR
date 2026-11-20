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

# KIT-063: the integration workflow's Results tab. It shows the standard results of the model the Model tab trained
# last (`workflow_data$integration_model`), read from the model's saved product (`integrationResults`): MOFA2's
# warnings in words, the variance explained, the factors against the groups with the claim limit beside them (KIT-060
# D10), the features with the largest weights, the factors' correlation; and saves them to the project on request
# (`integrationSaveResults`). It opens once a model is trained (the shell's status chain) and follows the latest model:
# a new model, or a new import, shows that model's results or none. Enrichment on the weights is KIT-141.

#' The Results tab's UI.
#' @param id Module id.
#' @noRd
mod_integration_results_ui <- function(id) {
    ns <- shiny::NS(id)
    texts <- integrationGuiTexts()
    shiny::fluidRow(shiny::column(12,
        shiny::h3(texts[["results_heading"]]),
        shiny::p(texts[["results_intro"]]),
        testid(shiny::uiOutput(ns("results_status")), "integration-results-status"),
        testid(shiny::uiOutput(ns("advisories")), "integration-results-advisories"),
        shiny::h4(texts[["variance_heading"]]),
        shiny::p(texts[["variance_intro"]]),
        testid(shiny::plotOutput(ns("variance_plot"), height = "420px"), "integration-results-variance-plot"),
        testid(shiny::tableOutput(ns("variance_table")), "integration-results-variance-table"),
        shiny::h5(texts[["variance_total_heading"]]),
        testid(shiny::plotOutput(ns("variance_total_plot"), height = "320px"), "integration-results-variance-total-plot"),
        shiny::h4(texts[["association_heading"]]),
        testid(shiny::uiOutput(ns("association_intro")), "integration-results-association-intro"),
        testid(shiny::uiOutput(ns("claim")), "integration-results-claim"),
        testid(shiny::plotOutput(ns("factors_plot"), height = "520px"), "integration-results-factors-plot"),
        testid(shiny::tableOutput(ns("association_table")), "integration-results-association-table"),
        shiny::h5(texts[["factor_values_heading"]]),
        testid(DT::DTOutput(ns("factor_values")), "integration-results-factor-values"),
        shiny::h4(texts[["weights_heading"]]),
        testid(shiny::uiOutput(ns("weights_intro")), "integration-results-weights-intro"),
        shiny::fluidRow(
            shiny::column(6, shiny::selectInput(ns("weights_factor"), texts[["weights_factor"]], choices = character())),
            shiny::column(6, shiny::selectInput(ns("weights_view"), texts[["weights_view"]], choices = character()))),
        testid(shiny::plotOutput(ns("weights_plot"), height = "420px"), "integration-results-weights-plot"),
        testid(shiny::tableOutput(ns("weights_table")), "integration-results-weights-table"),
        shiny::h4(texts[["correlation_heading"]]),
        shiny::p(texts[["correlation_intro"]]),
        testid(shiny::plotOutput(ns("correlation_plot"), height = "380px"), "integration-results-correlation-plot"),
        testid(shiny::tableOutput(ns("correlation_table")), "integration-results-correlation-table"),
        shiny::h4(texts[["save_heading"]]),
        shiny::p(shiny::em(texts[["results_provisional"]])),
        testid(shiny::actionButton(ns("save"), texts[["save_button"]], class = "btn-primary", icon = shiny::icon("save")),
            "integration-results-save"),
        testid(shiny::uiOutput(ns("save_status")), "integration-results-saved")))
}

#' The Results tab's server (KIT-063).
#' @param id Module id.
#' @param workflow_data The integration bus (`state_manager`, `integration_model`, `integration_results`, `tab_status`).
#' @param paths `project_dirs$integration` (`mofa_tables_dir`, `mofa_plots_dir`, `source_dir`).
#' @param texts `integrationGuiTexts()`.
#' @param parameters The Results parameters (absent ones take `integrationResultsParameterSpec()`'s defaults).
#' @return The results shown (`integrationResults`'s value, NULL, or a refusal), as a reactive.
#' @noRd
mod_integration_results_server <- function(id, workflow_data, paths, texts = integrationGuiTexts(), parameters = list()) {
    shiny::moduleServer(id, function(input, output, session) {
        saved <- shiny::reactiveValues(display = NULL, error = NULL)
        results <- shiny::reactive({
            model <- workflow_data$integration_model
            manager <- workflow_data$state_manager
            if (is.null(model) || is.null(manager)) return(NULL)
            tryCatch({
                value <- withWorkflowStateScope(manager, function()
                    hydrateSavedWorkflowEntity(manager$.context$store, model$reference)$value)
                integrationResults(value, integrationImportedObject(manager), parameters, texts)
            }, error = function(error) structure(list(refusal = integrationRefusalText(error, texts)),
                class = "integration_results_refused"))
        })
        usable <- function(value) !is.null(value) && !inherits(value, "integration_results_refused")
        # A new model (or none, after a new import): its results are not saved yet.
        shiny::observeEvent(workflow_data$integration_model, {
            saved$display <- NULL
            saved$error <- NULL
            workflow_data$integration_results <- NULL
        }, ignoreNULL = FALSE, ignoreInit = TRUE)
        shiny::observe({
            value <- results()
            if (!usable(value)) return()
            pick <- shiny::isolate(chosen(value))
            shiny::updateSelectInput(session, "weights_factor", choices = value$shown, selected = pick$factor)
            shiny::updateSelectInput(session, "weights_view", choices = value$views, selected = pick$view)
        })
        # The factor and data set the weights show: the chosen ones while they are the model's, else the first.
        chosen <- function(value) {
            factor <- input$weights_factor
            view <- input$weights_view
            list(factor = if (!is.null(factor) && factor %in% value$shown) factor else value$shown[[1L]],
                view = if (!is.null(view) && view %in% value$views) view else value$views[[1L]])
        }

        output$results_status <- shiny::renderUI({
            value <- results()
            if (is.null(value)) return(shiny::div(class = "alert alert-info", texts[["results_not_trained"]]))
            if (!usable(value)) {
                return(integrationRefusalUi(sprintf(texts[["results_refused"]], value$refusal$sentence), value$refusal$detail,
                    texts = texts))
            }
            model <- workflow_data$integration_model
            shiny::p(sprintf(texts[["results_model"]], format(model$trained_at, "%Y-%m-%d %H:%M")))
        })
        output$advisories <- shiny::renderUI({
            value <- results()
            if (!usable(value) || (!length(value$advisories) && !length(value$notes))) return(NULL)
            notes <- lapply(value$advisories, function(advisory) integrationRefusalUi(advisory$sentence, advisory$detail,
                class = "alert alert-warning", texts = texts))
            if (length(value$notes)) notes <- c(notes, list(shiny::div(class = "alert alert-warning",
                sprintf(texts[["few_samples"]], as.integer(value$shared_samples), INTEGRATION_FEW_SAMPLES))))
            shiny::tagList(shiny::h4(texts[["advisories_heading"]]), notes)
        })
        output$variance_plot <- shiny::renderPlot({
            value <- results()
            shiny::req(usable(value))
            integrationPlotVariance(value, texts)
        })
        output$variance_table <- shiny::renderTable({
            value <- results()
            shiny::req(usable(value))
            integrationResultsVarianceTable(value, texts)
        }, rownames = FALSE)
        output$variance_total_plot <- shiny::renderPlot({
            value <- results()
            shiny::req(usable(value))
            integrationPlotVarianceTotal(value, texts)
        })
        output$association_intro <- shiny::renderUI({
            value <- results()
            shiny::req(usable(value))
            shiny::p(sprintf(texts[["association_intro"]], texts[[paste0("padjust_", value$parameters$p_adjust_method)]]))
        })
        output$claim <- shiny::renderUI({
            value <- results()
            shiny::req(usable(value))
            lines <- integrationClaimLimitText(value, texts)
            sentences <- if (length(value$claim$group_corrected)) 2L else 1L
            shiny::div(class = "alert alert-secondary", shiny::h5(texts[["claim_heading"]]),
                lapply(lines[seq_len(sentences)], shiny::p), do.call(shiny::tags$ul, lapply(lines[-seq_len(sentences)], shiny::tags$li)))
        })
        output$factors_plot <- shiny::renderPlot({
            value <- results()
            shiny::req(usable(value))
            integrationPlotFactorsByGroup(value, texts)
        })
        output$association_table <- shiny::renderTable({
            value <- results()
            shiny::req(usable(value))
            integrationResultsAssociationTable(value, texts)
        }, rownames = FALSE)
        output$factor_values <- DT::renderDT({
            value <- results()
            shiny::req(usable(value))
            table <- value$factor_values
            for (name in value$factors) table[[name]] <- round(table[[name]], 3)
            names(table)[1:2] <- unname(texts[c("col_sample", "col_group")])
            DT::datatable(table, rownames = FALSE, selection = "none", options = list(pageLength = 25))
        })
        output$weights_intro <- shiny::renderUI({
            value <- results()
            shiny::req(usable(value))
            shiny::p(sprintf(texts[["weights_intro"]], value$parameters$top_features))
        })
        output$weights_plot <- shiny::renderPlot({
            value <- results()
            shiny::req(usable(value))
            pick <- chosen(value)
            integrationPlotWeights(value, pick$factor, pick$view, texts)
        })
        output$weights_table <- shiny::renderTable({
            value <- results()
            shiny::req(usable(value))
            pick <- chosen(value)
            integrationResultsWeightsTable(value, pick$factor, pick$view, texts)
        }, rownames = FALSE)
        output$correlation_plot <- shiny::renderPlot({
            value <- results()
            shiny::req(usable(value))
            integrationPlotFactorCorrelation(value, texts)
        })
        output$correlation_table <- shiny::renderTable({
            value <- results()
            shiny::req(usable(value))
            integrationResultsCorrelationTable(value, texts)
        }, rownames = FALSE)

        shiny::observeEvent(input$save, {
            saved$error <- NULL
            value <- results()
            model <- workflow_data$integration_model
            manager <- workflow_data$state_manager
            if (!usable(value) || is.null(model) || is.null(manager)) {
                saved$error <- list(sentence = texts[["results_not_trained"]], detail = NULL)
                return()
            }
            display <- tryCatch(integrationSaveResults(manager, model$reference, value, paths, texts),
                error = function(error) {
                    saved$error <- integrationRefusalText(error, texts)
                    NULL
                })
            if (is.null(display)) return()
            saved$display <- c(display, list(saved_at = Sys.time()))
            workflow_data$integration_results <- list(model = model$reference$entity_id, tables_dir = display$tables_dir,
                plots_dir = display$plots_dir, tables = display$tables, figures = display$figures, saved_at = saved$display$saved_at)
            tab_status <- shiny::isolate(workflow_data$tab_status)
            tab_status$integration_results <- "complete"
            workflow_data$tab_status <- tab_status
        })
        output$save_status <- shiny::renderUI({
            if (!is.null(saved$error)) {
                return(integrationRefusalUi(sprintf(texts[["save_refused"]], saved$error$sentence), saved$error$detail,
                    texts = texts))
            }
            display <- saved$display
            if (is.null(display)) return(NULL)
            shiny::div(class = "alert alert-success", sprintf(texts[["save_done"]], format(display$saved_at, "%Y-%m-%d %H:%M"),
                length(display$tables), display$tables_dir, length(display$figures), display$plots_dir))
        })
        results
    })
}
