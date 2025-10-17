#' @title peptideQCAppletServer (Orchestrator)
#'
#' @description Server orchestrator for the peptide-level quality control workflow.
#' This module reads the workflow type from the state manager and calls the
#' appropriate sub-modules for the given workflow (e.g., LFQ/DIA).
#'
#' @param id Module ID
#' @param workflow_data A reactive values object to store workflow data
#' @param experiment_paths A list of paths for the current experiment
#' @param omic_type The omics type (e.g., "proteomics")
#' @param experiment_label The experiment label
#' @export
#' @import shiny
#' @importFrom logger log_info

# Source all sub-modules for the peptide QC workflow
source("server/proteomics/qc/peptide/qvalue_filter_server.R")
source("server/proteomics/qc/peptide/precursor_rollup_server.R")
source("server/proteomics/qc/peptide/intensity_filter_server.R")
source("server/proteomics/qc/peptide/protein_peptide_filter_server.R")
source("server/proteomics/qc/peptide/sample_filter_server.R")
source("server/proteomics/qc/peptide/replicate_filter_server.R")
source("server/proteomics/qc/peptide/imputation_server.R")

peptideQCAppletServer <- function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    # Check workflow type to determine which modules to run
    workflow_type <- shiny::isolate(workflow_data$state_manager$workflow_type)
    
    log_info(paste("Peptide QC orchestrator: workflow type is", workflow_type))
    
    # Conditional execution for LFQ/DIA workflows which require peptide processing
    if (workflow_type %in% c("LFQ", "DIA")) {
      log_info("Running LFQ/DIA peptide processing sub-modules.")
      
      # Call each of the extracted sub-modules
      qvalue_filter_server(input, output, session, workflow_data, omic_type, experiment_label)
      precursor_rollup_server(input, output, session, workflow_data, omic_type, experiment_label)
      intensity_filter_server(input, output, session, workflow_data, omic_type, experiment_label)
      protein_peptide_filter_server(input, output, session, workflow_data, omic_type, experiment_label)
      sample_filter_server(input, output, session, workflow_data, omic_type, experiment_label)
      replicate_filter_server(input, output, session, workflow_data, omic_type, experiment_label)
      imputation_server(input, output, session, workflow_data, omic_type, experiment_label)
      
        } else {
      log_info(paste("Skipping peptide processing for workflow type:", workflow_type))
      # For TMT or other protein-level workflows, this orchestrator does nothing.
    }
    
  })
} 
