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
# KIT-104 W5: the exported wrappers over the shared intensity sub-tab (mod_omic_qc_tabs.R) under
# the lipid configuration.
# ============================================================================

#' @title Lipid Intensity Filter Module
#' @description A Shiny module for applying lipid intensity and missing value filters.
#'              Wraps the lipidIntensityFiltering() S4 method.
#' @name mod_lipid_qc_intensity
#' @param id,workflow_data,omic_type,experiment_label Runtime inputs used by this function; see the usage section for accepted values.
NULL

#' @rdname mod_lipid_qc_intensity
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 p hr numericInput sliderInput div actionButton verbatimTextOutput uiOutput
#' @importFrom shinyjqui jqui_resizable
mod_lipid_qc_intensity_ui <- function(id) {
    mod_omic_qc_intensity_ui(id, omicGuiConfig("lipidomics"))
}

#' @rdname mod_lipid_qc_intensity
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification removeNotification renderText renderPlot renderUI tabsetPanel tabPanel
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
mod_lipid_qc_intensity_server <- function(id, workflow_data, omic_type, experiment_label) {
    mod_omic_qc_intensity_server(id, workflow_data, omic_type, experiment_label, omicGuiConfig("lipidomics"))
}

#' The lipid intensity sub-tab's Revert, which trx's expression filter borrows (mod_trx_qc_filter.R) until KIT-104 W6:
#' the shared revert (`registerOmicQcRevertObserver`) clearing the two values given.
#' @noRd
registerLipidIntensityRevertObserver <- function(input, workflowData, output, filterStatsVal, filterPlotVal) {
    registerOmicQcRevertObserver(input, output, "revert_filter", workflowData, "filter_results", function() {
        filterStatsVal(NULL)
        filterPlotVal(NULL)
    }, "Reverted lipid intensity filter to")
}
