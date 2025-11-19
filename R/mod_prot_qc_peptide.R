#' @title Peptide QC Orchestrator Module
#'
#' @description Server orchestrator for the peptide-level quality control workflow.
#' This module reads the workflow type from the state manager and calls the
#' appropriate sub-modules for the given workflow (e.g., LFQ/DIA).
#'
#' @name mod_prot_qc_peptide
#' @export
NULL

#' @rdname mod_prot_qc_peptide
#' @export
#' @import shiny
#' @import shinydashboard
mod_prot_qc_peptide_ui <- function(id) {
  ns <- NS(id)
  
  # Nested tabs for each peptide filtering step
  shiny::tabsetPanel(
    id = ns("peptide_filter_tabs"),
    
    # Use the UI functions from the sub-modules
    # Note: We assume these are available in the package namespace
    
    # Step 1: Q-Value Filter
    if (exists("mod_prot_qc_peptide_qvalue_ui")) {
      mod_prot_qc_peptide_qvalue_ui(ns("qvalue_filter"))
    } else {
      shiny::tabPanel("Q-Value Filter", "Module not loaded")
    },
    
    # Step 2: Precursor Rollup
    if (exists("mod_prot_qc_peptide_rollup_ui")) {
      mod_prot_qc_peptide_rollup_ui(ns("rollup"))
    } else {
      shiny::tabPanel("Precursor Rollup", "Module not loaded")
    },
    
    # Step 3: Intensity Filter
    if (exists("mod_prot_qc_peptide_intensity_ui")) {
      mod_prot_qc_peptide_intensity_ui(ns("intensity_filter"))
    } else {
      shiny::tabPanel("Intensity Filter", "Module not loaded")
    },
    
    # Step 4: Protein Peptide Count Filter
    if (exists("mod_prot_qc_peptide_protein_ui")) {
      mod_prot_qc_peptide_protein_ui(ns("protein_peptide_filter"))
    } else {
      shiny::tabPanel("Protein Peptides", "Module not loaded")
    },
    
    # Step 5: Sample Quality Filter
    if (exists("mod_prot_qc_peptide_sample_ui")) {
      mod_prot_qc_peptide_sample_ui(ns("sample_filter"))
    } else {
      shiny::tabPanel("Sample Quality", "Module not loaded")
    },
    
    # Step 6: Replicate Filter
    if (exists("mod_prot_qc_peptide_replicate_ui")) {
      mod_prot_qc_peptide_replicate_ui(ns("replicate_filter"))
    } else {
      shiny::tabPanel("Replicate Filter", "Module not loaded")
    },
    
    # Step 7: Missing Value Imputation
    if (exists("mod_prot_qc_peptide_impute_ui")) {
      mod_prot_qc_peptide_impute_ui(ns("imputation"))
    } else {
      shiny::tabPanel("Imputation", "Module not loaded")
    }
  )
}

#' @rdname mod_prot_qc_peptide
#' @export
#' @import shiny
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
      
      # Call each of the extracted sub-modules
      # Ensure function existence checks before calling
      
      if (exists("mod_prot_qc_peptide_qvalue_server")) {
        mod_prot_qc_peptide_qvalue_server("qvalue_filter", workflow_data, omic_type, experiment_label)
      }
      
      if (exists("mod_prot_qc_peptide_rollup_server")) {
        mod_prot_qc_peptide_rollup_server("rollup", workflow_data, omic_type, experiment_label)
      }
      
      if (exists("mod_prot_qc_peptide_intensity_server")) {
        mod_prot_qc_peptide_intensity_server("intensity_filter", workflow_data, omic_type, experiment_label)
      }
      
      if (exists("mod_prot_qc_peptide_protein_server")) {
        mod_prot_qc_peptide_protein_server("protein_peptide_filter", workflow_data, omic_type, experiment_label)
      }
      
      if (exists("mod_prot_qc_peptide_sample_server")) {
        mod_prot_qc_peptide_sample_server("sample_filter", workflow_data, omic_type, experiment_label)
      }
      
      if (exists("mod_prot_qc_peptide_replicate_server")) {
        mod_prot_qc_peptide_replicate_server("replicate_filter", workflow_data, omic_type, experiment_label)
      }
      
      if (exists("mod_prot_qc_peptide_impute_server")) {
        mod_prot_qc_peptide_impute_server("imputation", workflow_data, omic_type, experiment_label)
      }
      
    } else {
      logger::log_info(paste("Skipping peptide processing for workflow type:", workflow_type))
      # For TMT or other protein-level workflows, this orchestrator does nothing.
    }
    
  })
}

