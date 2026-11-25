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

# KIT-140: the integration workflow's Model tab. The MOFA2 settings built from the model's parameter table
# (`integrationMofaParameterSpec`: every default is the table's, each marked with its status), the software check
# (`integrationMofaCapability`), the training on the confirmed sample table (KIT-139) through KIT-062's model step, and
# what the trained model is in a few lines: whether training settled, the factors kept, the variance explained per data
# set, the saved model's checksum. Under a managed application the training is a managed-worker action (the omics'
# runner, `installWorkflowOmicsActions`, installed by the shell); without one it runs in process
# (`runIntegrationModel`). Both run the one declaration and accept the model on the head the same way. Each model is
# kept in the project's integration store; a run with other settings is a new model, one with the same settings reuses
# the saved one. Anything richer about the model (plots, weights, factors against the design) is KIT-063's.

#' The Model tab's UI. The settings are built here from the parameter table, so a default changed in the table changes
#' the screen with no edit to this file.
#' @param id Module id.
#' @noRd
mod_integration_model_ui <- function(id) {
    ns <- shiny::NS(id)
    texts <- integrationGuiTexts()
    shiny::fluidRow(shiny::column(12,
        shiny::h3(texts[["model_heading"]]),
        shiny::p(texts[["model_intro"]]),
        testid(shiny::uiOutput(ns("model_status")), "integration-model-status"),
        shiny::h4(texts[["readiness_heading"]]),
        testid(shiny::uiOutput(ns("readiness")), "integration-model-readiness"),
        testid(shiny::actionButton(ns("recheck"), texts[["readiness_recheck"]], icon = shiny::icon("sync")),
            "integration-model-recheck"),
        shiny::h4(texts[["settings_heading"]]),
        shiny::h5(texts[["settings_matching_heading"]]),
        testid(shiny::uiOutput(ns("matching_settings")), "integration-model-matching"),
        testid(shiny::div(integrationModelSettingsUi(ns, texts = texts)), "integration-model-settings"),
        testid(shiny::actionButton(ns("run"), texts[["run_button"]], class = "btn-primary", icon = shiny::icon("play")),
            "integration-model-run"),
        testid(shiny::uiOutput(ns("run_status")), "integration-model-run-status"),
        testid(shiny::uiOutput(ns("result")), "integration-model-result"),
        testid(shiny::tableOutput(ns("variance_table")), "integration-model-variance")))
}

#' The Model tab's server (KIT-140).
#' @param id Module id.
#' @param workflow_data The integration bus (`state_manager`, `integration_alignment`, `integration_model`, `tab_status`).
#' @param paths `project_dirs$integration` (unused; the shell's signature).
#' @param texts `integrationGuiTexts()`.
#' @param probe The software check (`integrationMofaCapability`; a test may pass its own).
#' @param trainer NULL for MOFA2, or a stand-in with `integrationTrainMofa`'s form (tests; in process only).
#' @param withProgressFn The progress display around an in-process training (`shiny::withProgress`).
#' @return The tab's last model (`list(value, reference, status, parameters, trained_at)`) as a reactive.
#' @noRd
mod_integration_model_server <- function(id, workflow_data, paths, texts = integrationGuiTexts(),
    probe = integrationMofaCapability, trainer = NULL, withProgressFn = shiny::withProgress) {
    shiny::moduleServer(id, function(input, output, session) {
        spec <- integrationMofaParameterSpec()
        method <- setdiff(names(spec), INTEGRATION_GUI_MATCHING_PARAMETERS)
        model <- shiny::reactiveValues(result = NULL, error = NULL, running = FALSE, references = character())
        settings_values <- shiny::reactive(lapply(stats::setNames(method, method), function(name) input[[paste0("param_", name)]]))
        # The software is checked once the tab can train (a confirmed sample table) and again on request, not at start.
        capability <- shiny::reactive({
            input$recheck
            if (is.null(workflow_data$integration_alignment)) return(NULL)
            probe()
        })
        # A new import is a new integration: its models are counted afresh and an earlier one is no longer shown.
        shiny::observeEvent(workflow_data$state_manager, {
            model$result <- NULL
            model$error <- NULL
            model$running <- FALSE
            model$references <- character()
            workflow_data$integration_model <- NULL
        }, ignoreInit = TRUE)

        deliver <- function(value, reference, status, parameters) {
            model$references <- unique(c(model$references, reference$entity_id))
            model$result <- list(value = value, reference = reference, status = status, parameters = parameters,
                trained_at = Sys.time())
            model$running <- FALSE
            workflow_data$integration_model <- list(reference = reference, parameters = parameters, status = status,
                trained_at = model$result$trained_at, models = model$references)
            tab_status <- shiny::isolate(workflow_data$tab_status)
            tab_status$integration_model <- "complete"
            workflow_data$tab_status <- tab_status
        }
        refuse <- function(error) {
            model$running <- FALSE
            model$error <- integrationRefusalText(error, texts)
        }
        # The few-samples warning (D7) was shown at Sample Matching; the trained model keeps it as a note.
        quiet <- function(expr) withCallingHandlers(expr, integration_few_samples = function(warning) invokeRestart("muffleWarning"))

        shiny::observeEvent(input$run, {
            model$error <- NULL
            confirmed <- workflow_data$integration_alignment
            manager <- workflow_data$state_manager
            if (is.null(confirmed) || is.null(manager)) {
                model$error <- list(sentence = texts[["model_not_confirmed"]], detail = NULL)
                return()
            }
            tryCatch({
                parameters <- integrationModelParameters(confirmed$parameters, integrationModelSettings(settings_values(), spec, texts))
                # KIT-142: the confirmation of the sample table is kept with the model, for the record.
                confirmation <- integrationConfirmationRecord(confirmed)
                runner <- workflowOmicsManagedActions("integration")
                if (is.null(runner)) {
                    run <- withProgressFn(message = texts[["run_working"]], value = 0.5,
                        quiet(runIntegrationModel(manager, parameters, probe = probe, trainer = trainer, confirmation = confirmation)))
                    export <- manager$.context$snapshot$exports[[INTEGRATION_MODEL_OPERATION]]
                    deliver(run$value, run$reference, run$status, export$parameters)
                    return(invisible(NULL))
                }
                checked <- integrationModelPrecheck(manager, parameters, probe)
                model$running <- TRUE
                # KIT-148: the worker is handed the sample match as its table in the store (`checked$request`).
                runner$run(list(omicsManagedProduct(INTEGRATION_MODEL_OPERATION, checked$request, function(result) {
                    acceptIntegrationModel(manager, checked$source, list(reference = result$reference, value = result$value),
                        checked$parameters, confirmation)
                    deliver(result$value, result$reference, result$display$status %||% "computed", checked$parameters)
                })), failed = refuse)
            }, error = refuse)
        })

        output$model_status <- shiny::renderUI({
            if (is.null(workflow_data$integration_alignment)) {
                return(shiny::div(class = "alert alert-info", texts[["model_not_confirmed"]]))
            }
            NULL
        })
        output$readiness <- shiny::renderUI({
            value <- capability()
            if (is.null(value)) return(NULL)
            if (isTRUE(value$available)) return(shiny::div(class = "alert alert-success", texts[["readiness_ok"]]))
            # KIT-240: the missing piece named, and how to get it; a check that does not name one gives the general sentence.
            key <- paste0("readiness_missing_", tolower(value$missing %||% ""))
            sentence <- if (key %in% names(texts)) texts[[key]] else texts[["refusal_mofa_unavailable"]]
            integrationRefusalUi(sentence, value$reason, class = "alert alert-warning", texts = texts)
        })
        output$matching_settings <- shiny::renderUI({
            confirmed <- workflow_data$integration_alignment
            if (is.null(confirmed)) return(NULL)
            do.call(shiny::tags$ul, lapply(integrationModelMatchingSettings(confirmed$parameters, spec, texts), shiny::tags$li))
        })
        output$run_status <- shiny::renderUI({
            error <- model$error
            if (!is.null(error)) {
                return(integrationRefusalUi(sprintf(texts[["run_refused"]], error$sentence), error$detail, texts = texts))
            }
            if (isTRUE(model$running)) return(shiny::div(class = "alert alert-info", texts[["run_working"]]))
            result <- model$result
            if (is.null(result)) return(shiny::p(texts[["run_none"]]))
            changed <- tryCatch({
                settings <- integrationModelSettings(settings_values(), spec, texts)
                !all(vapply(names(settings), function(name) identical(settings[[name]], result$parameters[[name]]), logical(1)))
            }, error = function(error) TRUE)
            shiny::div(class = "alert alert-success",
                shiny::p(sprintf(texts[["run_done"]], format(result$trained_at, "%Y-%m-%d %H:%M"))),
                if (identical(result$status, "hit")) shiny::p(texts[["run_reused"]]),
                shiny::p(sprintf(texts[["run_count"]], length(model$references))),
                if (changed) shiny::p(shiny::strong(texts[["run_settings_changed"]])))
        })
        output$result <- shiny::renderUI({
            result <- model$result
            if (is.null(result)) return(NULL)
            summary <- integrationModelSummary(result$value, result$parameters, texts)
            shiny::tagList(shiny::h4(texts[["result_heading"]]),
                shiny::p(class = if (summary$converged) NULL else "text-danger", summary$convergence),
                shiny::p(summary$factors), shiny::p(summary$samples),
                # KIT-239: where the limitations of the model are stated.
                shiny::p(texts[["model_limits_pointer"]]),
                if (length(result$value$notes)) shiny::div(class = "alert alert-warning",
                    sprintf(texts[["few_samples"]], as.integer(result$value$shared_samples), INTEGRATION_FEW_SAMPLES)),
                if (!is.null(summary$saved)) shiny::p(summary$saved),
                shiny::h5(texts[["result_variance_heading"]]))
        })
        output$variance_table <- shiny::renderTable({
            result <- model$result
            shiny::req(result)
            integrationModelSummary(result$value, result$parameters, texts)$variance
        }, rownames = FALSE)

        shiny::reactive(model$result)
    })
}
