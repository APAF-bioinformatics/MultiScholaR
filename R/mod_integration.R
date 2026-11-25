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

# KIT-138: the integration workflow's shell. app_server finds `mod_integration_ui/server` by name, as it finds an
# omic's, when Integration is selected on the home page. The shell composes the single-omic shell's helpers itself
# (the bus, the stepper, the tab chain; `integrationGuiConfig()$shell`), since `validateOmicGuiConfig` requires the
# single-omic groups (KIT-060 audit section 5.3). Behind the gate `MultiScholaR.integration_gui`, on unless set FALSE
# (KIT-240); set FALSE, the module is the "not yet implemented" panel app_server showed before it existed. Recorded divergence
# from conventions/R-shiny.md (module function names, fluidRow and pill tabset, no bslib): the app's existing module
# pattern and layout are kept, under the existing-code rule of conventions/R.md.

#' The integration workflow's UI: the stepper and the tabs (KIT-138).
#' @param id Module id (`integration_workflow` in the app).
#' @noRd
mod_integration_ui <- function(id) {
    if (!integrationGuiEnabled()) return(integrationGuiUnavailablePanel())
    config <- integrationGuiConfig()
    shell <- config$shell
    ns <- shiny::NS(id)
    panels <- lapply(shell$tabs, function(tab) {
        ui <- omicGuiFunction(paste0("mod_", tab$component, "_ui"))
        shiny::tabPanel(title = shiny::tagList(shiny::icon(tab$icon), tab$title), value = tab$value,
            shiny::div(`data-testid` = tab$testid, ui(ns(tab$module_id))))
    })
    tabset <- testid(do.call(shiny::tabsetPanel, c(list(id = ns(shell$tabset_id), type = "pills"), panels)),
        paste0(shell$testid, "-workflow-tabs"))
    shiny::tagList(shiny::fluidRow(shiny::column(12, shiny::uiOutput(ns("workflow_progress")))), tabset)
}

#' The integration workflow's server (KIT-138): the bus, the tab servers in tab order, `active_tab`, the stepper and the
#' status chain. The integration store is opened by the Import tab when an import is committed.
#' @param id Module id.
#' @param project_dirs Project directories by omic (`setupDirectories()`); `project_dirs$integration` is used.
#' @param omic_type `integration`.
#' @param experiment_label The experiment label.
#' @param volumes Unused (the import reads the project's own folder).
#' @return The bus, or NULL with the gate off or without the integration directories.
#' @noRd
mod_integration_server <- function(id, project_dirs, omic_type = "integration", experiment_label = NULL, volumes = NULL) {
    if (!integrationGuiEnabled()) return(invisible(NULL))
    config <- integrationGuiConfig()
    shiny::moduleServer(id, function(input, output, session) {
        shell <- config$shell
        workflow_data <- do.call(shiny::reactiveValues, omicWorkflowInitialBus(shell))
        session$userData$workflow_data <- workflow_data
        if (!omic_type %in% names(project_dirs)) {
            logger::log_error(paste("No directory information found for", omic_type))
            return()
        }
        paths <- project_dirs[[omic_type]]
        workflow_data$active_tab <- shiny::reactive(input[[shell$tabset_id]])
        # KIT-140: under a managed application the Model tab trains on the managed worker through the omics' runner
        # (NULL without one: `run_app()`, a test; the tab then trains in process).
        if (is.character(experiment_label) && length(experiment_label) == 1L && !is.na(experiment_label) &&
            nzchar(experiment_label)) {
            installWorkflowOmicsActions(session, workflow_data, omic_type, experiment_label, paths)
        }
        session$onSessionEnded(function() {
            manager <- shiny::isolate(workflow_data$state_manager)
            if (!is.null(manager)) closeWorkflowStateOwner(manager)
        })
        for (tab in shell$tabs) {
            server <- omicGuiFunction(paste0("mod_", tab$component, "_server"))
            server(tab$module_id, workflow_data, paths, config$texts)
        }
        steps <- omicWorkflowSteps(shell$tabs)
        output$workflow_progress <- shiny::renderUI(render_workflow_stepper(steps, workflow_data$tab_status))
        installOmicWorkflowChain(workflow_data, shell)
        workflow_data
    })
}
