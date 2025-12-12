#' @title Protein QC Orchestrator Module
#'
#' @description Server orchestrator for the protein-level quality control workflow.
#' This module reads the workflow type from the state manager and calls the
#' appropriate sub-modules for the given workflow.
#'
#' @name mod_prot_qc_protein
NULL

#' @rdname mod_prot_qc_protein
#' @export
#' @importFrom shiny NS tagList tabsetPanel tabPanel
mod_prot_qc_protein_ui <- function(id, workflow_type = NULL) {
  ns <- shiny::NS(id)
  
  # Build tab list conditionally based on workflow type
  tab_list <- list()
  
  # Step 1: IQ Protein Rollup - ONLY for DIA workflows (peptide-to-protein)
  # Note: LFQ is now protein-level (like TMT), so it doesn't need rollup
  if (is.null(workflow_type) || workflow_type == "DIA") {
    if (exists("mod_prot_qc_protein_rollup_ui")) {
      tab_list[[length(tab_list) + 1]] <- mod_prot_qc_protein_rollup_ui(ns("rollup"))
    } else {
      tab_list[[length(tab_list) + 1]] <- shiny::tabPanel("IQ Protein Rollup", "Module not loaded")
    }
  }
  
  # Step 2-5: Common tabs for ALL workflows
  if (exists("mod_prot_qc_protein_cleanup_ui")) {
    tab_list[[length(tab_list) + 1]] <- mod_prot_qc_protein_cleanup_ui(ns("cleanup"))
  } else {
    tab_list[[length(tab_list) + 1]] <- shiny::tabPanel("Accession Cleanup", "Module not loaded")
  }
  
  if (exists("mod_prot_qc_protein_intensity_ui")) {
    tab_list[[length(tab_list) + 1]] <- mod_prot_qc_protein_intensity_ui(ns("intensity_filter"))
  } else {
    tab_list[[length(tab_list) + 1]] <- shiny::tabPanel("Protein Intensity Filter", "Module not loaded")
  }
  
  if (exists("mod_prot_qc_protein_dedup_ui")) {
    tab_list[[length(tab_list) + 1]] <- mod_prot_qc_protein_dedup_ui(ns("duplicate_removal"))
  } else {
    tab_list[[length(tab_list) + 1]] <- shiny::tabPanel("Duplicate Removal", "Module not loaded")
  }
  
  if (exists("mod_prot_qc_protein_replicate_ui")) {
    tab_list[[length(tab_list) + 1]] <- mod_prot_qc_protein_replicate_ui(ns("replicate_filter"))
  } else {
    tab_list[[length(tab_list) + 1]] <- shiny::tabPanel("Protein Replicate Filter", "Module not loaded")
  }
  
  # Construct tabsetPanel
  do.call(shiny::tabsetPanel, c(list(id = ns("protein_filter_tabs")), tab_list))
}

#' @rdname mod_prot_qc_protein
#' @export
#' @importFrom shiny moduleServer isolate
#' @importFrom logger log_info
mod_prot_qc_protein_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    # Check workflow type to determine which modules to run
    workflow_type <- shiny::isolate(workflow_data$state_manager$workflow_type)
    
    logger::log_info(paste("Protein QC orchestrator: workflow type is", workflow_type))
    
    # Conditional execution for DIA workflows which require protein rollup
    # Note: LFQ is now protein-level from the start (like TMT), so no rollup needed
    if (workflow_type == "DIA") {
      logger::log_info("Running DIA protein rollup sub-module.")
      
      if (exists("mod_prot_qc_protein_rollup_server")) {
        mod_prot_qc_protein_rollup_server("rollup", workflow_data, experiment_paths, omic_type, experiment_label)
      }
    }
    
    # These modules are common to all workflows (including TMT and LFQ)
    logger::log_info("Running common protein processing sub-modules.")
    
    if (exists("mod_prot_qc_protein_cleanup_server")) {
      mod_prot_qc_protein_cleanup_server("cleanup", workflow_data, omic_type, experiment_label)
    }
    
    if (exists("mod_prot_qc_protein_intensity_server")) {
      mod_prot_qc_protein_intensity_server("intensity_filter", workflow_data, omic_type, experiment_label)
    }
    
    if (exists("mod_prot_qc_protein_dedup_server")) {
      mod_prot_qc_protein_dedup_server("duplicate_removal", workflow_data, omic_type, experiment_label)
    }
    
    if (exists("mod_prot_qc_protein_replicate_server")) {
      mod_prot_qc_protein_replicate_server("replicate_filter", workflow_data, experiment_paths, omic_type, experiment_label)
    }
    
  })
}

