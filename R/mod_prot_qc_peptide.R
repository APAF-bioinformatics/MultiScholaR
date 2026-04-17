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
NULL

getProtQcPeptideModuleSpecs <- function() {
  list(
    list(
      id = "qvalue_filter",
      title = "Q-Value Filter",
      uiFn = "mod_prot_qc_peptide_qvalue_ui",
      serverFn = "mod_prot_qc_peptide_qvalue_server"
    ),
    list(
      id = "rollup",
      title = "Precursor Rollup",
      uiFn = "mod_prot_qc_peptide_rollup_ui",
      serverFn = "mod_prot_qc_peptide_rollup_server"
    ),
    list(
      id = "intensity_filter",
      title = "Intensity Filter",
      uiFn = "mod_prot_qc_peptide_intensity_ui",
      serverFn = "mod_prot_qc_peptide_intensity_server"
    ),
    list(
      id = "protein_peptide_filter",
      title = "Protein Peptides",
      uiFn = "mod_prot_qc_peptide_protein_ui",
      serverFn = "mod_prot_qc_peptide_protein_server"
    ),
    list(
      id = "sample_filter",
      title = "Sample Quality",
      uiFn = "mod_prot_qc_peptide_sample_ui",
      serverFn = "mod_prot_qc_peptide_sample_server"
    ),
    list(
      id = "replicate_filter",
      title = "Replicate Filter",
      uiFn = "mod_prot_qc_peptide_replicate_ui",
      serverFn = "mod_prot_qc_peptide_replicate_server"
    ),
    list(
      id = "imputation",
      title = "Imputation",
      uiFn = "mod_prot_qc_peptide_impute_ui",
      serverFn = "mod_prot_qc_peptide_impute_server"
    )
  )
}

buildProtQcPeptideTab <- function(moduleSpec, ns) {
  if (exists(moduleSpec$uiFn)) {
    get(moduleSpec$uiFn)(ns(moduleSpec$id))
  } else {
    shiny::tabPanel(moduleSpec$title, "Module not loaded")
  }
}

runProtQcPeptideSubmodule <- function(moduleSpec, workflow_data, omic_type, experiment_label) {
  if (exists(moduleSpec$serverFn)) {
    get(moduleSpec$serverFn)(moduleSpec$id, workflow_data, omic_type, experiment_label)
  }
}

#' @rdname mod_prot_qc_peptide
#' @export
#' @importFrom shiny NS tagList tabsetPanel tabPanel
mod_prot_qc_peptide_ui <- function(id) {
  ns <- shiny::NS(id)

  tabList <- lapply(getProtQcPeptideModuleSpecs(), buildProtQcPeptideTab, ns = ns)
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
    if (workflow_type == "DIA") {
      logger::log_info("Running DIA peptide processing sub-modules.")

      for (moduleSpec in getProtQcPeptideModuleSpecs()) {
        runProtQcPeptideSubmodule(
          moduleSpec = moduleSpec,
          workflow_data = workflow_data,
          omic_type = omic_type,
          experiment_label = experiment_label
        )
      }
    } else {
      logger::log_info(paste("Skipping peptide processing for workflow type:", workflow_type))
      # For TMT or other protein-level workflows, this orchestrator does nothing.
    }
    
  })
}
