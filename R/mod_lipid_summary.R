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
# mod_lipid_summary.R
# ============================================================================
# Purpose: Lipidomics session summary and integration export Shiny module
#
# This module provides a summary of the processing session and enables
# export of results for multi-omics integration.
# ============================================================================

#' @title Lipidomics Session Summary Module
#' @description A Shiny module for displaying processing summary and enabling
#'              data export for integration workflows.
#' @name mod_lipid_summary
#' @param id,project_dirs,omic_type,experiment_label,workflow_data Runtime inputs used by this function; see the usage section for accepted values.
NULL

#' @rdname mod_lipid_summary
#' @export
#' @importFrom shiny NS tagList fluidPage h3 fluidRow column wellPanel h4 textInput textAreaInput actionButton icon br hr downloadButton verbatimTextOutput checkboxInput conditionalPanel textOutput
#' @param profile An omic profile (KIT-097: `trxOmicProfile()`): its export panel, labels, test ids, integration export and
#'   report template; `NULL`, the default, is lipidomics, unchanged.
mod_lipid_summary_ui <- function(id, profile = NULL) {
    # KIT-108 W3 (D2): the shared tab (mod_omic_summary.R) under the omic's `summary` group; the profile selects the omic.
    mod_omic_summary_ui(id, omicGuiConfig(profile$omic_type %||% "lipidomics")$summary)
}
#' @rdname mod_lipid_summary
#' @export
mod_lipid_summary_server <- function(id, project_dirs, omic_type = "lipidomics", experiment_label = NULL, workflow_data = NULL, profile = NULL) {
    mod_omic_summary_server(id, project_dirs, omic_type, experiment_label, workflow_data,
        omicGuiConfig(profile$omic_type %||% "lipidomics")$summary)
}
