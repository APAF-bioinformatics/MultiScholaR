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
# mod_lipid_qc_intensity.R
# ============================================================================
# Purpose: Lipid intensity filtering Shiny module
#
# This module wraps the lipidIntensityFiltering() S4 method to provide
# an interactive UI for filtering lipids based on intensity thresholds.
# ============================================================================

#' @title Lipid Intensity Filter Module
#' @description A Shiny module for applying lipid intensity and missing value filters.
#'              Wraps the lipidIntensityFiltering() S4 method.
#' @name mod_lipid_qc_intensity
NULL

#' @rdname mod_lipid_qc_intensity
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 p hr numericInput sliderInput div actionButton verbatimTextOutput uiOutput
#' @importFrom shinyjqui jqui_resizable
mod_lipid_qc_intensity_ui <- function(id) {
    ns <- shiny::NS(id)
    
    shiny::tabPanel(
        "Intensity Filter"
        , shiny::br()
        , shiny::fluidRow(
            shiny::column(4
                , shiny::wellPanel(
                    shiny::h4("Lipid Intensity Filter Parameters")
                    , shiny::p("Filter lipids based on intensity thresholds. Lipids with low intensities across many samples will be removed.")
                    , shiny::hr()
                    
                    # Intensity cutoff percentile
                    , shiny::sliderInput(
                        ns("intensity_cutoff_percentile")
                        , "Intensity Cutoff Percentile (%)"
                        , min = 1
                        , max = 50
                        , value = 10
                        , step = 1
                    )
                    , shiny::helpText("Intensities below this percentile are considered 'low' (default: 10%)")
                    
                    # Proportion threshold
                    , shiny::sliderInput(
                        ns("proportion_below_cutoff")
                        , "Max Proportion of Samples Below Threshold"
                        , min = 0.1
                        , max = 1.0
                        , value = 0.5
                        , step = 0.05
                    )
                    , shiny::helpText("Remove lipids where this proportion of samples are below the intensity threshold (default: 0.5)")
                    
                    , shiny::hr()
                    , shiny::div(
                        shiny::actionButton(
                            ns("apply_filter")
                            , "Apply Filter"
                            , class = "btn-primary"
                            , width = "48%"
                        )
                        , shiny::actionButton(
                            ns("revert_filter")
                            , "Revert"
                            , class = "btn-warning"
                            , width = "48%"
                            , style = "margin-left: 4%"
                        )
                    )
                )
            )
            , shiny::column(8
                , shiny::verbatimTextOutput(ns("filter_results"))
                , shiny::br()
                # Per-assay results tabs
                , shiny::uiOutput(ns("assay_results_tabs"))
                , shiny::br()
                , shinyjqui::jqui_resizable(
                    shiny::plotOutput(ns("filter_plot"), height = "600px", width = "100%")
                )
            )
        )
    )
}





#' @rdname mod_lipid_qc_intensity
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification removeNotification renderText renderPlot renderUI tabsetPanel tabPanel
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
mod_lipid_qc_intensity_server <- function(id, workflow_data, omic_type, experiment_label) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns
        
        filter_plot <- shiny::reactiveVal(NULL)
        filter_stats <- shiny::reactiveVal(NULL)

        registerLipidIntensityApplyFilterObserver(
            input = input,
            workflowData = workflow_data,
            output = output,
            filterStatsVal = filter_stats,
            filterPlotVal = filter_plot,
            omicType = omic_type
        )
        
        registerLipidIntensityRevertObserver(
            input = input,
            workflowData = workflow_data,
            output = output,
            filterStatsVal = filter_stats,
            filterPlotVal = filter_plot
        )
        
        registerLipidIntensityAssayResultsOutput(
            output = output,
            filterStatsVal = filter_stats,
            ns = ns
        )

        registerLipidIntensityFilterPlotOutput(
            output = output,
            filterPlotVal = filter_plot
        )
    })
}
