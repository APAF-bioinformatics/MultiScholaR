#' @title proteinQCAppletServer (Orchestrator)
#'
#' @description Server orchestrator for the protein-level quality control workflow.
#' This module reads the workflow type from the state manager and calls the
#' appropriate sub-modules for the given workflow.
#'
#' @param id Module ID
#' @param workflow_data A reactive values object to store workflow data
#' @param experiment_paths A list of paths for the current experiment
#' @param omic_type The omics type (e.g., "proteomics")
#' @param experiment_label The experiment label
#' @export
#' @import shiny
#' @importFrom logger log_info

# Source all sub-modules for the protein QC workflow
source("server/proteomics/qc/protein/protein_rollup_server.R")
source("server/proteomics/qc/protein/accession_cleanup_server.R")
source("server/proteomics/qc/protein/protein_intensity_filter_server.R")
source("server/proteomics/qc/protein/duplicate_removal_server.R")
source("server/proteomics/qc/protein/protein_replicate_filter_server.R")

proteinQCAppletServer <- function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    # Check workflow type to determine which modules to run
    workflow_type <- shiny::isolate(workflow_data$state_manager$workflow_type)
    
    log_info(paste("Protein QC orchestrator: workflow type is", workflow_type))
    
    # Conditional execution for LFQ/DIA workflows
    if (workflow_type %in% c("LFQ", "DIA")) {
      log_info("Running LFQ/DIA protein processing sub-modules.")
      
      # Call LFQ/DIA-specific rollup module
      protein_rollup_server("protein_rollup", workflow_data, experiment_paths, omic_type, experiment_label)
    }
    
    # These modules are common to all workflows (including TMT)
    log_info("Running common protein processing sub-modules.")
    accession_cleanup_server(input, output, session, workflow_data, omic_type, experiment_label)
    protein_intensity_filter_server(input, output, session, workflow_data, omic_type, experiment_label)
    duplicate_removal_server(input, output, session, workflow_data, omic_type, experiment_label)
    protein_replicate_filter_server(input, output, session, workflow_data, experiment_paths, omic_type, experiment_label)
    
  })
} 
