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
# mod_metab_qc.R
# ============================================================================
# Purpose: Metabolomics QC orchestrator module
#
# KIT-104 W6: the exported wrappers over the shared QC shell (mod_omic_qc.R) under the metab
# configuration (`omicGuiConfig("metabolomics")`): its title, tabset, four shared sub-tabs, their
# start (the design's trigger, or the state's class) and the finalise's alignment check.
# ============================================================================

#' @title Metabolomics Quality Control Orchestrator Module
#' @description A Shiny module coordinator for the Quality Control step of the
#'              metabolomics workflow. Integrates sub-modules for intensity filtering,
#'              duplicate resolution, internal standard QC, and finalization.
#' @name mod_metab_qc
#' @param id,workflow_data,experiment_paths,omic_type,experiment_label,qc_trigger Runtime inputs used by this function; see the usage section for accepted values.
NULL

#' @rdname mod_metab_qc
#' @export
#' @importFrom shiny NS tagList fluidRow column wellPanel h3 uiOutput
mod_metab_qc_ui <- function(id) {
    mod_omic_qc_ui(id, omicGuiConfig("metabolomics"))
}

#' @rdname mod_metab_qc
#' @export
#' @importFrom shiny moduleServer reactive observeEvent req renderUI tabsetPanel tabPanel
#' @importFrom logger log_info log_error
mod_metab_qc_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label, qc_trigger = NULL) {
    mod_omic_qc_server(id, workflow_data, experiment_paths, omic_type, experiment_label, qc_trigger,
        config = omicGuiConfig("metabolomics"))
}
