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
# mod_lipid_qc_s4.R
# ============================================================================
# Purpose: Final S4 state confirmation module for lipidomics QC
#
# KIT-104 W5: the exported wrappers over the shared finalise sub-tab (mod_omic_qc_tabs.R) under the
# lipid configuration, or the profile's omic (trx) under its own (KIT-104 W6 removed trx's own
# finalise path).
# ============================================================================

#' @title Lipidomics QC Finalization Module
#' @description A Shiny module for finalizing the QC step and confirming the S4 state.
#'              Displays QC summary statistics and state history.
#' @name mod_lipid_qc_s4
#' @param id,workflow_data,omic_type,experiment_label Runtime inputs used by this function; see the usage section for accepted values.
#' @param profile An omic profile (KIT-096): its `omic_type` selects that omic's configuration. `NULL`, the default, is lipidomics.
NULL

#' @rdname mod_lipid_qc_s4
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 h5 p hr div actionButton verbatimTextOutput uiOutput icon tags plotOutput
#' @importFrom DT DTOutput
#' @importFrom shinyjqui jqui_resizable
mod_lipid_qc_s4_ui <- function(id, profile = NULL) {
    mod_omic_qc_finalize_ui(id, omicGuiConfig(profile$omic_type %||% "lipidomics"))
}

#' @rdname mod_lipid_qc_s4
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification renderText renderUI observe tags renderPlot
#' @importFrom DT renderDT datatable
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
mod_lipid_qc_s4_server <- function(id, workflow_data, omic_type, experiment_label, profile = NULL) {
    mod_omic_qc_finalize_server(id, workflow_data, omic_type, experiment_label,
        omicGuiConfig(profile$omic_type %||% "lipidomics"))
}
