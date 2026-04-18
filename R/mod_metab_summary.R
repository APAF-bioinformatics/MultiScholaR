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
# mod_metab_summary.R
# ============================================================================
# Purpose: Metabolomics session summary and integration export Shiny module
#
# This module provides a summary of the processing session and enables
# export of results for multi-omics integration.
# ============================================================================

#' @title Metabolomics Session Summary Module
#' @description A Shiny module for displaying processing summary and enabling
#'              data export for integration workflows.
#' @name mod_metab_summary
NULL

#' @rdname mod_metab_summary
#' @export
#' @importFrom shiny NS tagList fluidPage h3 fluidRow column wellPanel h4 textInput textAreaInput actionButton icon br hr downloadButton verbatimTextOutput checkboxInput conditionalPanel textOutput
mod_metab_summary_ui <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::fluidPage(
            shiny::h3("Session Summary & Report Generation"),

            shiny::fluidRow(
                # Left column: Workflow Parameters
                shiny::column(6,
                    shiny::wellPanel(
                        shiny::h4("Workflow Parameters"),
                        shiny::textInput(ns("experiment_label"), "Experiment Label:",
                                        value = "", placeholder = "e.g., my_metabolomics_analysis"),
                        shiny::textAreaInput(ns("description"), "Description:",
                                            value = "Full metabolomics analysis workflow with normalization and DA",
                                            rows = 3, resize = "vertical"),
                        shiny::br(),
                        shiny::actionButton(ns("save_workflow_args"), "Save Workflow Arguments",
                                           class = "btn-primary", icon = shiny::icon("save"))
                    )
                ),

                # Right column: File Management
                shiny::column(6,
                    shiny::wellPanel(
                        shiny::h4("File Management"),
                        shiny::br(),
                        shiny::actionButton(ns("copy_to_publication"), "Copy to Publication Directory",
                                           class = "btn-info", icon = shiny::icon("copy")),
                        shiny::br(), shiny::br(),
                        shiny::verbatimTextOutput(ns("copy_status"))
                    )
                )
            ),

            shiny::fluidRow(
                # Left column: Report Generation
                shiny::column(6,
                    shiny::wellPanel(
                        shiny::h4("Report Generation"),
                        shiny::textOutput(ns("template_status")),
                        shiny::br(),
                        shiny::actionButton(ns("generate_report"), "Generate Report",
                                           class = "btn-success", icon = shiny::icon("file-pdf")),
                        shiny::br(), shiny::br(),
                        shiny::conditionalPanel(
                            condition = paste0("output['", ns("report_ready"), "']"),
                            shiny::downloadButton(ns("download_report"), "Download Report",
                                                 class = "btn-success")
                        )
                    )
                ),

                # Right column: GitHub Integration (Optional)
                shiny::column(6,
                    shiny::wellPanel(
                        shiny::h4("Version Control (Optional)"),
                        shiny::checkboxInput(ns("enable_github"), "Enable GitHub Push", FALSE),
                        shiny::conditionalPanel(
                            condition = paste0("input['", ns("enable_github"), "']"),
                            shiny::textInput(ns("github_org"), "GitHub Organization:", ""),
                            shiny::textInput(ns("github_email"), "GitHub Email:", ""),
                            shiny::textInput(ns("github_username"), "GitHub Username:", ""),
                            shiny::textInput(ns("project_id"), "Project ID:", ""),
                            shiny::br(),
                            shiny::actionButton(ns("push_to_github"), "Push to GitHub",
                                               class = "btn-warning", icon = shiny::icon("github"))
                        )
                    )
                )
            ),

            shiny::fluidRow(
                shiny::column(12,
                    shiny::wellPanel(
                        shiny::h4("Session Summary"),
                        shiny::verbatimTextOutput(ns("session_summary")),
                        shiny::br(),
                        shiny::actionButton(ns("export_session_state"), "Export Session State (.RDS)",
                                           class = "btn-secondary", icon = shiny::icon("download"))
                    )
                )
            )
        )
    )
}










#' @rdname mod_metab_summary
#' @export
#' @importFrom shiny moduleServer reactive reactiveValues observeEvent req renderText showNotification downloadHandler withProgress incProgress
#' @importFrom logger log_info log_error
#' @importFrom utils write.table
#' @importFrom purrr detect
#' @importFrom readr read_tsv
mod_metab_summary_server <- function(id, project_dirs, omic_type = "metabolomics", experiment_label = NULL, workflow_data = NULL) {
    shiny::moduleServer(id, function(input, output, session) {

        values <- setupMetabSummaryServerBootstrapState(
            session = session,
            experimentLabel = experiment_label
        )

        setupMetabSummaryTemplateStatusOutput(
            output = output,
            projectDirs = project_dirs,
            omicType = omic_type
        )

        registerMetabSummaryServerObservers(
            input = input,
            output = output,
            projectDirs = project_dirs,
            omicType = omic_type,
            workflowData = workflow_data,
            values = values
        )

        setupMetabSummaryBootstrapOutputs(output = output)

    })
}
