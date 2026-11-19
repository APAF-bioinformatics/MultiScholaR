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

# KIT-138: the integration workflow's Import tab. No upload: the tab scans the project's integration folder
# (`<project>/integration/<omic>/`, KIT-060 D1) for the inputs the omics saved at their Summary tab, lists what it found
# per omic and per data set with the correction for unwanted variation each names, lists what cannot be used with the
# contract reader's message, and commits the selected omics through `importIntegrationInputs` (KIT-062).

#' The Import tab's UI.
#' @param id Module id.
#' @noRd
mod_integration_import_ui <- function(id) {
    ns <- shiny::NS(id)
    texts <- integrationGuiTexts()
    shiny::fluidRow(shiny::column(12,
        shiny::h3(texts[["import_heading"]]),
        shiny::p(texts[["import_intro"]]),
        shiny::uiOutput(ns("scan_summary")),
        testid(shiny::actionButton(ns("rescan"), texts[["import_rescan"]], icon = shiny::icon("sync")),
            "integration-import-rescan"),
        shiny::h4(texts[["omics_heading"]]),
        testid(DT::DTOutput(ns("omics_table")), "integration-import-omics"),
        shiny::uiOutput(ns("not_ready")),
        shiny::h4(texts[["views_heading"]]),
        testid(DT::DTOutput(ns("views_table")), "integration-import-views"),
        shiny::uiOutput(ns("view_choice")),
        testid(shiny::actionButton(ns("import"), texts[["import_button"]], class = "btn-primary",
            icon = shiny::icon("file-import")), "integration-import-commit"),
        shiny::uiOutput(ns("import_status"))))
}

#' The Import tab's server (KIT-138).
#' @param id Module id.
#' @param workflow_data The integration bus.
#' @param paths `project_dirs$integration` (`integration_dir`, `source_dir`).
#' @param texts `integrationGuiTexts()`.
#' @param newOwner Opens the integration store under the project's integration scripts folder (a test may pass its own).
#' @return The scan, as a reactive.
#' @noRd
mod_integration_import_server <- function(id, workflow_data, paths, texts = integrationGuiTexts(),
    newOwner = integrationSessionOwner) {
    shiny::moduleServer(id, function(input, output, session) {
        status <- shiny::reactiveVal(NULL)
        imported <- shiny::reactive({
            summary <- workflow_data$integration_import
            if (is.null(summary)) NULL else summary$manifests
        })
        scan <- shiny::eventReactive(list(input$rescan, imported()), {
            integrationScanProject(paths$integration_dir, imported(), texts)
        }, ignoreNULL = FALSE)
        shiny::observe(workflow_data$integration_scan <- scan())

        output$scan_summary <- shiny::renderUI({
            found <- scan()
            ready <- sum(found$omics$status == "ready")
            text <- if (nrow(found$omics)) sprintf(texts[["import_found"]], ready, nrow(found$omics) - ready) else
                texts[["import_none"]]
            shiny::tagList(shiny::p(sprintf(texts[["import_folder"]], found$dir)), shiny::p(text))
        })
        output$omics_table <- DT::renderDT({
            DT::datatable(integrationScanDisplay(scan(), texts), rownames = FALSE, selection = "none",
                options = list(dom = "t", paging = FALSE, ordering = FALSE))
        })
        output$not_ready <- shiny::renderUI({
            omics <- scan()$omics
            rows <- omics[omics$status %in% c("not_saved", "refused", "changed"), , drop = FALSE]
            if (!nrow(rows)) return(NULL)
            items <- lapply(seq_len(nrow(rows)), function(i) {
                shiny::tags$li(sprintf(texts[["refused_reason"]], tools::toTitleCase(rows$omic[[i]]), rows$reason[[i]]))
            })
            testid(shiny::div(class = "alert alert-warning", shiny::strong(texts[["refused_heading"]]),
                do.call(shiny::tags$ul, items)), "integration-import-not-ready")
        })
        output$views_table <- DT::renderDT({
            views <- scan()$views
            table <- data.frame(views$view, tools::toTitleCase(views$omic), views$features, views$runs,
                views$biological_samples, views$groups, stringsAsFactors = FALSE)
            names(table) <- unname(texts[c("col_view", "col_omic", "col_features", "col_runs", "col_biological_samples",
                "col_groups")])
            DT::datatable(table, rownames = FALSE, selection = "none",
                options = list(dom = "t", paging = FALSE, ordering = FALSE))
        })
        output$view_choice <- shiny::renderUI({
            views <- integrationScanReadyViews(scan())
            if (!nrow(views)) return(NULL)
            labels <- sprintf(texts[["view_choice_item"]], views$view, tools::toTitleCase(views$omic),
                as.integer(views$features), as.integer(views$biological_samples))
            chosen <- shiny::isolate(workflow_data$integration_import$views) %||% views$view
            testid(shiny::checkboxGroupInput(session$ns("views"), texts[["view_choice_label"]],
                choiceNames = labels, choiceValues = views$view, selected = intersect(chosen, views$view)),
                "integration-import-choice")
        })
        output$import_status <- shiny::renderUI({
            message <- status()
            if (is.null(message)) return(NULL)
            testid(shiny::div(class = paste("alert", if (identical(message$type, "ok")) "alert-success" else "alert-danger"),
                message$text), "integration-import-status")
        })

        shiny::observeEvent(input$import, {
            result <- tryCatch({
                request <- integrationImportRequest(scan(), input$views %||% character(), texts)
                manager <- shiny::isolate(workflow_data$state_manager)
                if (is.null(manager)) {
                    manager <- tryCatch(newOwner(paths), error = function(error) {
                        integrationRefuse(sprintf(texts[["store_unavailable"]], conditionMessage(error)),
                            c("integration_store_unavailable", "integration_gui_refusal"))
                    })
                    workflow_data$state_manager <- manager
                }
                integrationCommitImport(manager, request)
            }, integration_refused = function(error) error, integration_input_invalid = function(error) error)
            if (inherits(result, "error")) {
                # The tab's own refusals are sentences already; the import's are prefixed.
                text <- if (inherits(result, "integration_gui_refusal")) conditionMessage(result) else
                    sprintf(texts[["import_refused"]], conditionMessage(result))
                status(list(type = "error", text = text))
                return()
            }
            workflow_data$integration_import <- result
            workflow_data$tab_status <- integrationImportStatus(shiny::isolate(workflow_data$tab_status))
            status(list(type = "ok", text = sprintf(texts[["import_done"]], length(result$omics), length(result$views),
                result$runs, result$biological_samples)))
        })
        scan
    })
}

#' The tab statuses after a committed import: the import complete, the tab after it open and every later tab closed
#' (each later step is for the import it was made on; the shell's chain opens them again in order).
#' @param status The bus's `tab_status`.
#' @return The statuses.
#' @noRd
integrationImportStatus <- function(status) {
    keys <- names(status)
    at <- match("integration_import", keys)
    status[["integration_import"]] <- "complete"
    later <- keys[seq_along(keys) > at]
    if (length(later)) {
        status[[later[[1L]]]] <- "pending"
        for (key in later[-1L]) status[[key]] <- "disabled"
    }
    status
}

#' Open the integration store for an app session (KIT-138): a fresh folder under the project's integration scripts
#' folder, as the omics' session owners are opened (`createWorkflowSessionOwner`).
#' @param paths `project_dirs$integration`.
#' @return The store (`newIntegrationStateOwner`).
#' @noRd
integrationSessionOwner <- function(paths) {
    if (!is.character(paths$source_dir) || length(paths$source_dir) != 1L || !dir.exists(paths$source_dir)) {
        cli::cli_abort("the project has no integration scripts folder", class = "saved_entity_scope")
    }
    if (!workflowOwnerCapable()) {
        cli::cli_abort("the qs2, DuckDB and Arrow packages are needed", class = "saved_entity_unsupported")
    }
    root <- tempfile(".integration-native-", tmpdir = normalizePath(paths$source_dir, winslash = "/", mustWork = TRUE))
    newIntegrationStateOwner(root, "integration")
}
