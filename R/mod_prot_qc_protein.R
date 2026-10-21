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

#' @title Protein QC Orchestrator Module
#'
#' @description Server orchestrator for the protein-level quality control workflow.
#' This module reads the workflow type from the state manager and calls the
#' appropriate sub-modules for the given workflow.
#'
#' @name mod_prot_qc_protein
#' @param id,workflow_type,workflow_data,experiment_paths,omic_type,experiment_label Runtime inputs used by this function; see the usage section for accepted values.
NULL

getProtQcProteinRollupModuleSpec <- function() {
  list(
    id = "rollup",
    title = "IQ Protein Rollup",
    uiFn = "mod_prot_qc_protein_rollup_ui",
    serverFn = "mod_prot_qc_protein_rollup_server",
    needsExperimentPaths = TRUE
  )
}

getProtQcProteinCommonModuleSpecs <- function() {
  list(
    list(
      id = "cleanup",
      title = "Accession Cleanup",
      uiFn = "mod_prot_qc_protein_cleanup_ui",
      serverFn = "mod_prot_qc_protein_cleanup_server",
      needsExperimentPaths = TRUE
    ),
    list(
      id = "intensity_filter",
      title = "Protein Intensity Filter",
      uiFn = "mod_prot_qc_protein_intensity_ui",
      serverFn = "mod_prot_qc_protein_intensity_server",
      needsExperimentPaths = TRUE
    ),
    list(
      id = "duplicate_removal",
      title = "Duplicate Removal",
      uiFn = "mod_prot_qc_protein_dedup_ui",
      serverFn = "mod_prot_qc_protein_dedup_server",
      needsExperimentPaths = TRUE
    ),
    list(
      id = "replicate_filter",
      title = "Protein Replicate Filter",
      uiFn = "mod_prot_qc_protein_replicate_ui",
      serverFn = "mod_prot_qc_protein_replicate_server",
      needsExperimentPaths = TRUE
    )
  )
}

getProtQcProteinModuleSpecs <- function(workflowType = NULL) {
  # KIT-046: the DIA_LIMPA route's protein stage (the script's checkpoints): DPC quantification (KIT-047), the
  # optional between-sample normalisation and the peptide-count filter (KIT-041); placeholders fix the ids here.
  if (identical(workflowType, "DIA_LIMPA")) {
    return(list(
      list(id = "dpc_quant", title = "DPC Quantification", uiFn = "mod_prot_qc_protein_dpc_quant_ui", serverFn = "mod_prot_qc_protein_dpc_quant_server", needsExperimentPaths = TRUE, placeholder = "KIT-047"),
      list(id = "normalisation", title = "Between-Sample Normalisation", uiFn = "mod_prot_qc_protein_normalisation_ui", serverFn = "mod_prot_qc_protein_normalisation_server", needsExperimentPaths = TRUE, placeholder = "KIT-041"),
      list(id = "peptide_count", title = "Peptide Count Filter", uiFn = "mod_prot_qc_protein_peptide_count_ui", serverFn = "mod_prot_qc_protein_peptide_count_server", needsExperimentPaths = TRUE, placeholder = "KIT-041")
    ))
  }
  moduleSpecs <- getProtQcProteinCommonModuleSpecs()

  if (is.null(workflowType) || identical(workflowType, "DIA")) {
    moduleSpecs <- c(list(getProtQcProteinRollupModuleSpec()), moduleSpecs)
  }

  moduleSpecs
}

buildProtQcProteinTab <- function(moduleSpec, ns) {
  if (exists(moduleSpec$uiFn)) {
    get(moduleSpec$uiFn)(ns(moduleSpec$id))
  } else if (!is.null(moduleSpec$placeholder)) {
    shiny::tabPanel(moduleSpec$title, shiny::p(sprintf("%s: declared in %s.", moduleSpec$title, moduleSpec$placeholder)))
  } else {
    shiny::tabPanel(moduleSpec$title, "Module not loaded")
  }
}

runProtQcProteinSubmodule <- function(moduleSpec,
                                      workflow_data,
                                      experiment_paths,
                                      omic_type,
                                      experiment_label) {
  if (!exists(moduleSpec$serverFn)) {
    return(invisible(NULL))
  }

  # KIT-015: every protein module takes the experiment paths (its managed step's publication context) in the
  # third position, as the rollup and replicate modules always did.
  callArgs <- list(moduleSpec$id, workflow_data)
  if (isTRUE(moduleSpec$needsExperimentPaths)) {
    callArgs <- c(callArgs, list(experiment_paths))
  }
  callArgs <- c(callArgs, list(omic_type, experiment_label))

  do.call(get(moduleSpec$serverFn), callArgs)
}

#' @rdname mod_prot_qc_protein
#' @export
#' @importFrom shiny NS tagList tabsetPanel tabPanel
mod_prot_qc_protein_ui <- function(id, workflow_type = NULL) {
  ns <- shiny::NS(id)

  tabList <- lapply(
    getProtQcProteinModuleSpecs(workflowType = workflow_type),
    buildProtQcProteinTab,
    ns = ns
  )
  do.call(shiny::tabsetPanel, c(list(id = ns("protein_filter_tabs")), tabList))
}

#' @rdname mod_prot_qc_protein
#' @export
#' @importFrom shiny moduleServer isolate
#' @importFrom logger log_info
mod_prot_qc_protein_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    workflow_type <- shiny::isolate(workflow_data$state_manager$workflow_type)
    logger::log_info(paste("Protein QC orchestrator: workflow type is", workflow_type))

    if (workflow_type == "DIA") {
      logger::log_info("Running DIA protein rollup sub-module.")
    }
    logger::log_info("Running common protein processing sub-modules.")

    for (moduleSpec in getProtQcProteinModuleSpecs(workflowType = workflow_type)) {
      runProtQcProteinSubmodule(
        moduleSpec = moduleSpec,
        workflow_data = workflow_data,
        experiment_paths = experiment_paths,
        omic_type = omic_type,
        experiment_label = experiment_label
      )
    }
  })
}
