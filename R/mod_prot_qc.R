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

# KIT-104 W7: proteomics' QC tab is the shared QC shell (mod_omic_qc.R) under
# `omicGuiConfig("proteomics")`: the title, the tabsets by route, the "Peptide QC" and "Protein QC"
# groups whose plugins are mod_prot_qc_{peptide,protein} (looked up by name), the start on the
# design's trigger only, the route in the shell's local `qc_workflow_type`, and the completion on the
# route's final QC state through `completeProtQcWorkflowStatus`. The coordinator these wrappers
# replace is at 827f2328.

#' @title qualityControlAppletModule
#'
#' @description A Shiny module coordinator for the Quality Control step of the proteomics
#' workflow. This is now a lightweight coordinator that integrates extracted components
#' while maintaining backward compatibility.
#'
#' @name qualityControlAppletModule
#' @param id,workflowData,logInfoFn Runtime inputs used by this function; see the usage section for accepted values.
NULL

#' @rdname qualityControlAppletModule
#' @export
#' @importFrom shiny NS tagList fluidRow column wellPanel h3 uiOutput
mod_prot_qc_ui <- function(id) {
  mod_omic_qc_ui(id, omicGuiConfig("proteomics"))
}

#' @rdname qualityControlAppletModule
#' @export
#' @importFrom shiny moduleServer reactive observeEvent req renderUI showNotification removeNotification
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
completeProtQcWorkflowStatus <- function(
    workflowData,
    logInfoFn = logger::log_info
) {
  updatedStatus <- workflowData$tab_status
  if (is.null(updatedStatus)) {
    updatedStatus <- list()
  }

  updatedStatus$quality_control <- "complete"
  if (is.null(updatedStatus$normalization) || identical(updatedStatus$normalization, "disabled")) {
    updatedStatus$normalization <- "pending"
  }
  workflowData$tab_status <- updatedStatus

  logInfoFn("Quality Control complete. Normalization is ready.")
  invisible(updatedStatus)
}

mod_prot_qc_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label, qc_trigger = NULL) {
  mod_omic_qc_server(id, workflow_data, experiment_paths, omic_type, experiment_label, qc_trigger,
    config = omicGuiConfig("proteomics"))
}

#' One proteomics QC group's spec table on a route, from the configuration (KIT-104 D9): each tab on `route` as the
#' group plugins read it, `list(id, title, uiFn, serverFn)` plus `needsExperimentPaths = TRUE` where the tab takes the
#' experiment paths (the protein steps), in the configuration's order.
#' @param group_id The group's id (`peptide_qc`, `protein_qc`).
#' @param route The route whose tabs are listed.
#' @return A list of module specs (`getProtQc{Peptide,Protein}ModuleSpecs`).
#' @noRd
protQcModuleSpecs <- function(group_id, route) {
  groups <- omicGuiConfig("proteomics")$qc$groups
  group <- Filter(function(group) identical(group$id, group_id), groups)[[1L]]
  tabs <- Filter(function(tab) route %in% tab$routes, group$tabs)
  lapply(tabs, function(tab) {
    c(list(id = tab$id, title = tab$title, uiFn = paste0("mod_", tab$component, "_ui"),
      serverFn = paste0("mod_", tab$component, "_server")),
      if (isTRUE(tab$needs_paths)) list(needsExperimentPaths = TRUE))
  })
}
