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
# mod_metab_qc_s4.R
# ============================================================================
# Purpose: Final S4 state confirmation module for metabolomics QC
#
# This module provides a summary of QC processing and allows the user to
# finalize the QC step before proceeding to normalization.
# ============================================================================

#' @title Metabolomics QC Finalization Module
#' @description A Shiny module for finalizing the QC step and confirming the S4 state.
#'              Displays QC summary statistics and state history.
#' @name mod_metab_qc_s4
NULL

#' @rdname mod_metab_qc_s4
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 h5 p hr div actionButton verbatimTextOutput uiOutput icon tags plotOutput
#' @importFrom DT DTOutput
#' @importFrom shinyjqui jqui_resizable
mod_metab_qc_s4_ui <- function(id) {
    ns <- shiny::NS(id)
    
    shiny::tabPanel(
        "Finalize QC"
        , shiny::br()
        , shiny::fluidRow(
            shiny::column(6
                , shiny::wellPanel(
                    shiny::h4("QC Processing Summary")
                    , shiny::p("Review the QC processing steps before finalizing. You can revert to any previous state if needed.")
                    , shiny::hr()
                    
                    # State history
                    , shiny::h5("Processing History")
                    , shiny::uiOutput(ns("state_history"))
                    
                    , shiny::hr()
                    , shiny::div(
                        shiny::actionButton(
                            ns("finalize_qc")
                            , "Finalize QC"
                            , class = "btn-success"
                            , width = "100%"
                            , icon = shiny::icon("check")
                        )
                    )
                )
            )
            , shiny::column(6
                , shiny::wellPanel(
                    shiny::h4("Current Data Summary")
                    , shiny::hr()
                    , shiny::uiOutput(ns("data_summary"))
                )
                , shiny::wellPanel(
                    shiny::h4("Per-Assay Statistics")
                    , shiny::hr()
                    , DT::DTOutput(ns("assay_stats_table"))
                )
            )
        )
        , shiny::fluidRow(
            shiny::column(12
                , shiny::verbatimTextOutput(ns("finalize_results"))
                , shiny::br()
                # QC Progress Grid
                , shinyjqui::jqui_resizable(
                    shiny::plotOutput(ns("filter_plot"), height = "600px", width = "100%")
                )
            )
        )
    )
}























#' @rdname mod_metab_qc_s4
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification renderText renderUI observe tags renderPlot
#' @importFrom DT renderDT datatable
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
mod_metab_qc_s4_server <- function(id, workflow_data, omic_type, experiment_label) {
    shiny::moduleServer(id, function(input, output, session) {
        runMetabQcS4ServerBody(
            input = input,
            output = output,
            session = session,
            workflowData = workflow_data,
            omicType = omic_type,
            experimentLabel = experiment_label
        )
    })
}
