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

# ============================================================================
# mod_metabolomics.R
# ============================================================================
# Purpose: Top-level metabolomics workflow orchestrator Shiny module
#
# This module coordinates all metabolomics sub-modules in a tabbed interface,
# managing the complete workflow from import to export.
# ============================================================================

#' @title Metabolomics Workflow Orchestrator Module
#' @description Top-level Shiny module that coordinates the complete metabolomics
#'              workflow including import, design, QC, normalization, differential
#'              analysis, and summary/export.
#' @name mod_metabolomics
NULL

#' @rdname mod_metabolomics
#' @export
#' @importFrom shiny NS tagList fluidRow column h2 tabsetPanel tabPanel uiOutput icon
mod_metabolomics_ui <- function(id) {
    # KIT-102: the shared shell (mod_omic_workflow.R) under the metabolomics configuration.
    mod_omic_workflow_ui(id, omicGuiConfig("metabolomics"))
}

#' @rdname mod_metabolomics
#' @param id Module ID
#' @param project_dirs List of project directories from setup
#' @param omic_type The omics type (should be "metabolomics")
#' @param experiment_label The experiment label for this analysis
#' @param volumes Volumes for shinyFiles file browser
#' @export
#' @importFrom shiny moduleServer reactiveValues reactive observeEvent renderUI req tags reactiveVal
#' @importFrom logger log_info log_error log_warn
mod_metabolomics_server <- function(id, project_dirs, omic_type, experiment_label, volumes = NULL) {
    # KIT-102: the shared shell (mod_omic_workflow.R) under the metabolomics configuration.
    mod_omic_workflow_server(id, project_dirs, omic_type, experiment_label, volumes, omicGuiConfig("metabolomics"))
}


# ============================================================================
# Standalone app launcher for testing
# ============================================================================

#' @title Launch Metabolomics Workflow App
#' @description Launches a standalone Shiny app for the metabolomics workflow.
#'              Useful for testing and development.
#'
#' @param base_dir Base directory for the project. Defaults to a temp directory.
#' @param ... Additional arguments passed to shiny::shinyApp
#'
#' @return A Shiny app object
#' @export
#'
#' @examples
#' \dontrun{
#' run_metabolomics_app()
#' }
run_metabolomics_app <- function(base_dir = NULL, ...) {
    # Set up default base directory
    if (is.null(base_dir)) {
        base_dir <- file.path(tempdir(), "metabolomics_test")
        if (!dir.exists(base_dir)) {
            dir.create(base_dir, recursive = TRUE)
        }
    }

    # Create mock project_dirs structure matching what app_server.R creates
    project_dirs <- list(
        metabolomics = list(
            base_dir = base_dir,
            data_dir = file.path(base_dir, "data"),
            results_dir = file.path(base_dir, "results"),
            source_dir = file.path(base_dir, "scripts")
        )
    )

    # Create directories if they don't exist
    lapply(project_dirs$metabolomics, function(d) {
        if (!dir.exists(d)) dir.create(d, recursive = TRUE)
    })

    ui <- shiny::fluidPage(
        shiny::tags$head(
            shiny::tags$style(shiny::HTML("
                body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; }
                .nav-pills > li > a { border-radius: 0; }
                .badge { padding: 8px 12px; }
            "))
        ),
        mod_metabolomics_ui("metabolomics_app")
    )

    server <- function(input, output, session) {
        # Get volumes for shinyFiles
        volumes <- NULL
        if (requireNamespace("shinyFiles", quietly = TRUE)) {
            volumes <- shinyFiles::getVolumes()()
        }

        mod_metabolomics_server(
            id = "metabolomics_app",
            project_dirs = project_dirs,
            omic_type = "metabolomics",
            experiment_label = "Metabolomics Test",
            volumes = volumes
        )
    }

    shiny::shinyApp(ui, server, ...)
}
