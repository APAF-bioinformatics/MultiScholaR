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
# mod_lipid_qc.R
# ============================================================================
# Purpose: Lipidomics QC orchestrator module
#
# KIT-104 W5: the exported wrappers over the shared QC shell (mod_omic_qc.R) under the lipid
# configuration, or the profile's omic (trx) under its own (KIT-104 W6: its tab is wholly the shared
# shell's, with the filter plugin and the shared finalise).
# ============================================================================

#' @title Lipidomics Quality Control Orchestrator Module
#' @description A Shiny module coordinator for the Quality Control step of the
#'              lipidomics workflow. Integrates sub-modules for intensity filtering,
#'              duplicate resolution, internal standard QC, and finalization.
#' @name mod_lipid_qc
#' @param id,workflow_data,omic_type,experiment_label Runtime inputs used by this function; see the usage section for accepted values.
#' @param profile An omic profile (KIT-096): its `omic_type` selects that omic's configuration. `NULL`, the default, is lipidomics.
NULL

#' @rdname mod_lipid_qc
#' @export
#' @importFrom shiny NS tagList fluidRow column wellPanel h3 uiOutput
mod_lipid_qc_ui <- function(id, profile = NULL) {
    # KIT-104 W5: the shared shell under the lipid configuration, or the profile's omic (trx) under its own.
    mod_omic_qc_ui(id, omicGuiConfig(profile$omic_type %||% "lipidomics"))
}

#' @rdname mod_lipid_qc
#' @export
#' @importFrom shiny moduleServer reactive observeEvent req renderUI tabsetPanel tabPanel
#' @importFrom logger log_info log_error
initializeLipidQcSubmodules <- function(workflow_data, omic_type, experiment_label) {
    # KIT-104 W5: the four shared sub-tabs, started by the shell's start function.
    config <- omicGuiConfig("lipidomics")
    startOmicQcSubtabs(workflow_data, NULL, omic_type, experiment_label, config, config$shell$routes[[1L]])
}

mod_lipid_qc_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label, qc_trigger = NULL, profile = NULL) {
    # KIT-104 W5/W6: the shared shell under the lipid configuration, or the profile's omic (trx) under its own.
    mod_omic_qc_server(id, workflow_data, experiment_paths, omic_type, experiment_label, qc_trigger,
        config = omicGuiConfig(profile$omic_type %||% "lipidomics"))
}
