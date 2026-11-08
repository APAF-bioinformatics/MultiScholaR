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
# KIT-104 W5: the exported wrappers over the shared duplicate sub-tab (mod_omic_qc_tabs.R) under
# the lipid configuration.
# ============================================================================

#' @title Lipid Duplicate Resolution Module
#' @description A Shiny module for identifying and resolving duplicate lipid features.
#'              Wraps the resolveDuplicateFeaturesByIntensity() function.
#' @name mod_lipid_qc_duplicates
#' @param id,workflow_data,omic_type,experiment_label Runtime inputs used by this function; see the usage section for accepted values.
NULL

#' @rdname mod_lipid_qc_duplicates
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 p hr div actionButton verbatimTextOutput uiOutput plotOutput
#' @importFrom DT DTOutput
#' @importFrom shinyjqui jqui_resizable
mod_lipid_qc_duplicates_ui <- function(id) {
    mod_omic_qc_duplicates_ui(id, omicGuiConfig("lipidomics"))
}

#' @rdname mod_lipid_qc_duplicates
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification removeNotification renderText renderUI tabsetPanel tabPanel tags renderPlot
#' @importFrom DT renderDT datatable
#' @importFrom logger log_info log_error log_warn
#' @importFrom grid grid.draw
mod_lipid_qc_duplicates_server <- function(id, workflow_data, omic_type, experiment_label) {
    mod_omic_qc_duplicates_server(id, workflow_data, omic_type, experiment_label, omicGuiConfig("lipidomics"))
}
