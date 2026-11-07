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
# mod_metab_design.R
# ============================================================================
# Purpose: the metabolomics design tab's exported names. KIT-103 W5: thin wrappers over the shared
# design tab (mod_omic_design.R) under the metab configuration (`omicGuiConfig("metabolomics")`); the
# import, commit, builder registration and previews are the shared frame's.
# ============================================================================

#' @title Metabolomics Design Matrix Applet Module
#'
#' @description A Shiny module that serves as the main host for the design
#' matrix creation workflow step. Embeds mod_metab_design_builder_ui/server
#' and handles Import Existing Design, S4 object creation, and state saving.
#' Architecture is 1:1 with mod_prot_design.R.
#'
#' @param id Module ID
#'
#' @name metabolomicsDesignMatrixAppletModule
#' @param output,workflowData,renderDt,renderText,req Runtime inputs used by this function; see the usage section for accepted values.
NULL

#' @rdname metabolomicsDesignMatrixAppletModule
#' @export
#' @importFrom shiny NS tagList wellPanel h3 h4 p conditionalPanel div icon tags HTML fluidRow column actionButton
#' @importFrom DT DTOutput renderDT
mod_metab_design_ui <- function(id) {
    # KIT-103 W5: the shared design tab (mod_omic_design.R) with the metab texts.
    mod_omic_design_ui(id, omicGuiConfig("metabolomics"))
}

#' @rdname metabolomicsDesignMatrixAppletModule
#' @export
#' @importFrom shiny moduleServer reactive observeEvent req renderUI showNotification removeNotification outputOptions renderText
#' @importFrom logger log_info log_error log_warn
#' @importFrom utils write.table
#' @importFrom vroom vroom
#' @importFrom shinyFiles shinyDirButton shinyDirChoose parseDirPath getVolumes
#' @importFrom jsonlite write_json read_json
registerMetabDesignPreviewOutputs <- function(
    output,
    workflowData,
    renderDt = DT::renderDT,
    renderText = shiny::renderText,
    req = shiny::req
) {
    # KIT-103 W5: the shared previews (registerOmicsDesignPreviewOutputs) under the metab
    # configuration. `renderDt`, `renderText` and `req` stay in the signature for its callers and
    # are not used: the previews are shiny's.
    registerOmicsDesignPreviewOutputs(output, workflowData, omicGuiConfig("metabolomics"))
}

#' The metab design tab's server: the shared design tab (mod_omic_design.R) under the metab
#' configuration (KIT-103 W5).
#' @noRd
mod_metab_design_server <- function(id, workflow_data, experiment_paths, volumes = NULL, qc_trigger = NULL) {
    mod_omic_design_server(id, workflow_data, experiment_paths, volumes = volumes, qc_trigger = qc_trigger,
        config = omicGuiConfig("metabolomics"))
}
