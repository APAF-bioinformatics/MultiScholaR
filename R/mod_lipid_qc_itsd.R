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
# mod_lipid_qc_itsd.R
# ============================================================================
# Purpose: Internal standard (ITSD) QC visualization Shiny module
#
# KIT-104 W5: the exported wrappers over the shared internal-standard sub-tab (mod_omic_qc_tabs.R)
# under the lipid configuration.
# ============================================================================

#' @title Internal Standard QC Visualization Module
#' @description A Shiny module for visualizing internal standard quality metrics.
#'              Displays IS detection, CV distributions, and intensity trends.
#' @name mod_lipid_qc_itsd
#' @param id,workflow_data,omic_type,experiment_label Runtime inputs used by this function; see the usage section for accepted values.
NULL

#' @rdname mod_lipid_qc_itsd
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 h5 p hr textInput actionButton verbatimTextOutput uiOutput plotOutput
#' @importFrom shinyjqui jqui_resizable
mod_lipid_qc_itsd_ui <- function(id) {
    mod_omic_qc_itsd_ui(id, omicGuiConfig("lipidomics"))
}

#' @rdname mod_lipid_qc_itsd
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification renderText renderUI renderPlot tabsetPanel tabPanel tags
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_point geom_line labs theme_minimal theme element_text coord_flip scale_fill_brewer facet_wrap geom_hline
#' @importFrom logger log_info log_error log_warn
mod_lipid_qc_itsd_server <- function(id, workflow_data, omic_type, experiment_label) {
    mod_omic_qc_itsd_server(id, workflow_data, omic_type, experiment_label, omicGuiConfig("lipidomics"))
}
