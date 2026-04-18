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
# mod_lipid_qc_duplicates.R
# ============================================================================
# Purpose: Lipid duplicate resolution Shiny module
#
# This module wraps the resolveDuplicateFeaturesByIntensity() function to
# provide an interactive UI for identifying and resolving duplicate lipids.
# ============================================================================

#' @title Lipid Duplicate Resolution Module
#' @description A Shiny module for identifying and resolving duplicate lipid features.
#'              Wraps the resolveDuplicateFeaturesByIntensity() function.
#' @name mod_lipid_qc_duplicates
NULL

#' @rdname mod_lipid_qc_duplicates
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 p hr div actionButton verbatimTextOutput uiOutput plotOutput
#' @importFrom DT DTOutput
#' @importFrom shinyjqui jqui_resizable
mod_lipid_qc_duplicates_ui <- function(id) {
    ns <- shiny::NS(id)

    shiny::tabPanel(
        "Duplicate Resolution"
        , shiny::br()
        , shiny::fluidRow(
            shiny::column(4
                , shiny::wellPanel(
                    shiny::h4("Duplicate Feature Resolution")
                    , shiny::p("Identify and resolve duplicate lipid IDs within each assay. Duplicates are resolved by keeping the feature with the highest average intensity across samples.")
                    , shiny::hr()

                    # Detection summary
                    , shiny::h5("Detection Summary")
                    , shiny::uiOutput(ns("duplicate_summary"))

                    , shiny::hr()
                    , shiny::div(
                        shiny::actionButton(
                            ns("detect_duplicates")
                            , "Detect Duplicates"
                            , class = "btn-info"
                            , width = "100%"
                        )
                    )
                    , shiny::br()
                    , shiny::div(
                        shiny::actionButton(
                            ns("resolve_duplicates")
                            , "Resolve Duplicates"
                            , class = "btn-primary"
                            , width = "48%"
                        )
                        , shiny::actionButton(
                            ns("revert_duplicates")
                            , "Revert"
                            , class = "btn-warning"
                            , width = "48%"
                            , style = "margin-left: 4%"
                        )
                    )
                )
            )
            , shiny::column(8
                , shiny::verbatimTextOutput(ns("resolution_results"))
                , shiny::br()
                # Per-assay duplicate tables
                , shiny::uiOutput(ns("duplicate_tables"))
                , shiny::br()
                # QC Progress Grid
                , shinyjqui::jqui_resizable(
                    shiny::plotOutput(ns("filter_plot"), height = "600px", width = "100%")
                )
            )
        )
    )
}






















#' @rdname mod_lipid_qc_duplicates
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification removeNotification renderText renderUI tabsetPanel tabPanel tags renderPlot
#' @importFrom DT renderDT datatable
#' @importFrom logger log_info log_error log_warn
#' @importFrom grid grid.draw
mod_lipid_qc_duplicates_server <- function(id, workflow_data, omic_type, experiment_label) {
    shiny::moduleServer(id, function(input, output, session) {
        runLipidDuplicateModuleServerShell(
            input = input
            , output = output
            , session = session
            , workflowData = workflow_data
            , omicType = omic_type
        )
    })
}
