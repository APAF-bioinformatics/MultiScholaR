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

#' @title Peptide QC Orchestrator Module
#'
#' @description Server orchestrator for the peptide-level quality control workflow.
#' This module reads the workflow type from the state manager and calls the
#' appropriate sub-modules for the given workflow (e.g., LFQ/DIA).
#'
#' @name mod_prot_qc_peptide
#' @param id,workflow_data,experiment_paths,omic_type,experiment_label Runtime inputs used by this function; see the usage section for accepted values.
#' @param workflowType The workflow type whose peptide step list the UI renders (`DIA` or `DIA_LIMPA`, KIT-046).
NULL

getProtQcPeptideModuleSpecs <- function(workflowType = "DIA") {
  # KIT-046: the DIA_LIMPA route's peptide stage (the script's checkpoints): the Q-value filter and the precursor
  # rollup are the shared declarations; the log2 matrix and the peptide-level accession cleanup are KIT-047's
  # declarations and modules. KIT-104 W7 (D9): the table is the configuration's peptide group; every route other
  # than DIA_LIMPA lists the DIA steps, as the literal table did (827f2328).
  protQcModuleSpecs("peptide_qc", if (identical(workflowType, "DIA_LIMPA")) "DIA_LIMPA" else "DIA")
}

buildProtQcPeptideTab <- function(moduleSpec, ns) {
  if (exists(moduleSpec$uiFn)) {
    get(moduleSpec$uiFn)(ns(moduleSpec$id))
  } else if (!is.null(moduleSpec$placeholder)) {
    shiny::tabPanel(moduleSpec$title, shiny::p(sprintf("%s: declared in %s.", moduleSpec$title, moduleSpec$placeholder)))
  } else {
    shiny::tabPanel(moduleSpec$title, "Module not loaded")
  }
}

runProtQcPeptideSubmodule <- function(moduleSpec, workflow_data, omic_type, experiment_label,
                                    experiment_paths = NULL) {
  if (exists(moduleSpec$serverFn)) {
    arguments <- list(moduleSpec$id, workflow_data, omic_type, experiment_label)
    if (!is.null(experiment_paths)) arguments$experiment_paths <- experiment_paths
    do.call(get(moduleSpec$serverFn), arguments)
  }
}

#' @rdname mod_prot_qc_peptide
#' @export
#' @importFrom shiny NS tagList tabsetPanel tabPanel
mod_prot_qc_peptide_ui <- function(id, workflowType = "DIA") {
  ns <- shiny::NS(id)

  tabList <- lapply(getProtQcPeptideModuleSpecs(workflowType), buildProtQcPeptideTab, ns = ns)
  do.call(shiny::tabsetPanel, c(list(id = ns("peptide_filter_tabs")), tabList))
}

#' @rdname mod_prot_qc_peptide
#' @export
#' @importFrom shiny moduleServer isolate
#' @importFrom logger log_info
mod_prot_qc_peptide_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    # Check workflow type to determine which modules to run
    workflow_type <- shiny::isolate(workflow_data$state_manager$workflow_type)
    
    logger::log_info(paste("Peptide QC orchestrator: workflow type is", workflow_type))
    
    # Conditional execution for DIA workflows which require peptide processing
    # Note: LFQ is now a protein-level workflow (like TMT) and bypasses peptide QC
    if (workflow_type %in% c("DIA", "DIA_LIMPA")) { # KIT-046: both DIA-NN routes are peptide-level
      logger::log_info(paste("Running", workflow_type, "peptide processing sub-modules."))

      for (moduleSpec in getProtQcPeptideModuleSpecs(workflow_type)) {
        runProtQcPeptideSubmodule(
          moduleSpec = moduleSpec,
          workflow_data = workflow_data,
          omic_type = omic_type,
          experiment_label = experiment_label,
          experiment_paths = experiment_paths
        )
      }
    } else {
      logger::log_info(paste("Skipping peptide processing for workflow type:", workflow_type))
      # For TMT or other protein-level workflows, this orchestrator does nothing.
    }
    
  })
}
