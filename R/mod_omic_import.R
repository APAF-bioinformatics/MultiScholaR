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

# KIT-102 W4: the shared import frame. It renders the configured title and source panel in the markup the three import
# modules share today (shinyjs, a full-width row, a well panel under an h3), calls the panel's server in its own
# session, and offers the panels one completion helper. It adds no DOM id: the panel's inputs and outputs are in the
# frame's namespace, as today's are in the import module's.
#
# The source-panel plugin contract (the panels are KIT-102 phase C's: the feature table W5, proteomics W6):
#   ui(ns, config)                 the tags inside the well panel, after the title;
#   server(input, output, session, workflow_data, experiment_paths, volumes, config, complete)
#                                  called inside the frame's moduleServer; its value is the frame's;
#   complete(notify = TRUE)        the completion helper (`completeOmicImport`): the panel calls it after its commit.

#' The import frame's UI.
#' @param id Module id (the import tab's).
#' @param config An omic GUI configuration (`omicGuiConfig()`).
#' @noRd
mod_omic_import_ui <- function(id, config) {
    validateOmicGuiConfig(config)
    ns <- shiny::NS(id)
    panel <- omicImportSourcePanel(config$import$source_panel)
    shiny::tagList(
        shinyjs::useShinyjs(),
        shiny::fluidRow(shiny::column(12, shiny::wellPanel(
            shiny::h3(config$import$title),
            omicGuiFunction(panel$ui)(ns, config)))))
}

#' The import frame's server: the configured panel's server, given the completion helper.
#' @param id Module id.
#' @param workflow_data The workflow bus.
#' @param experiment_paths The experiment's project paths.
#' @param volumes shinyFiles volumes, passed to the panel.
#' @param config An omic GUI configuration (`omicGuiConfig()`).
#' @return What the panel's server returns.
#' @noRd
mod_omic_import_server <- function(id, workflow_data, experiment_paths, volumes = NULL, config) {
    validateOmicGuiConfig(config)
    panel <- omicImportSourcePanel(config$import$source_panel)
    shiny::moduleServer(id, function(input, output, session) {
        complete <- function(notify = TRUE) completeOmicImport(workflow_data, config, panel$design_status, notify)
        omicGuiFunction(panel$server)(input, output, session, workflow_data = workflow_data,
            experiment_paths = experiment_paths, volumes = volumes, config = config, complete = complete)
    })
}

#' The source panels by id: their UI and server (function names, looked up when called) and how their import opens
#' the Design step (`design_status`): `open_if_disabled` (the lipid rule, which metab joins: KIT-102's unification), or
#' `pending` (proteomics sets it unconditionally, completeProtImportSuccessState).
#' @param source_panel `config$import$source_panel`.
#' @noRd
omicImportSourcePanel <- function(source_panel) {
    panels <- list(
        feature_table = list(ui = "omicImportFeatureTableUi", server = "omicImportFeatureTableServer",
            design_status = "open_if_disabled"),
        proteomics = list(ui = "protImportPanelUi", server = "protImportPanelServer", design_status = "pending"))
    if (!is.character(source_panel) || length(source_panel) != 1L || !source_panel %in% names(panels)) {
        cli::cli_abort("{.val {source_panel}} is not an import source panel.", class = "omic_gui_config_invalid")
    }
    panels[[source_panel]]
}

#' The import's completion, in one `tab_status` write: the import tab complete and the step after it opened (by the
#' panel's `design_status`); then the configured success notification.
#' @param workflow_data The workflow bus.
#' @param config An omic GUI configuration: the first two tabs' status keys and `import$success_text`.
#' @param design_status `open_if_disabled` or `pending` (`omicImportSourcePanel`).
#' @param notify Whether to show the success notification.
#' @return The new statuses, invisibly.
#' @noRd
completeOmicImport <- function(workflow_data, config, design_status = "open_if_disabled", notify = TRUE,
                               showNotification = shiny::showNotification) {
    import_key <- config$shell$tabs[[1L]]$status_key
    design_key <- config$shell$tabs[[2L]]$status_key
    status <- workflow_data$tab_status
    status[[import_key]] <- "complete"
    if (identical(design_status, "pending") || identical(status[[design_key]], "disabled")) {
        status[[design_key]] <- "pending"
    }
    workflow_data$tab_status <- status
    if (isTRUE(notify)) showNotification(config$import$success_text, type = "message")
    invisible(status)
}
