# Import extracted components at top of file
# Note: These will need to be available in the package or sourced appropriately
# Since we are moving to R/, they should be loaded by package loading mechanisms
# or we might need to source them if they are not functions but scripts (which is bad practice)
# Given the plan, we should assume they will be refactored into functions in R/ as well.
# For now, we will assume the sub-modules are available as functions.

#' @title qualityControlAppletModule
#'
#' @description A Shiny module coordinator for the Quality Control step of the proteomics
#' workflow. This is now a lightweight coordinator that integrates extracted components
#' while maintaining backward compatibility.
#'
#' @name qualityControlAppletModule
#' @export
NULL

#' @rdname qualityControlAppletModule
#' @export
#' @import shiny
#' @import shinydashboard
mod_prot_qc_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::fluidRow(
    shiny::column(12,
      shiny::wellPanel(
        shiny::h3("Quality Control & Filtering"),
        # Replace the conditional logic with a single UI Output
        shiny::uiOutput(ns("dynamic_qc_tabs"))
      )
    )
  )
}

#' @rdname qualityControlAppletModule
#' @export
#' @import shiny
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
mod_prot_qc_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label, qc_trigger = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    output$dynamic_qc_tabs <- shiny::renderUI({
      message("=== DEBUG66: mod_prot_qc_server renderUI fired ===")
      
      # This logic will run once workflow_type is available
      shiny::req(workflow_data$state_manager$workflow_type)
      workflow_type <- workflow_data$state_manager$workflow_type
      message(sprintf("   DEBUG66: workflow_type = %s", workflow_type))
      
      if (workflow_type %in% c("TMT", "LFQ")) {
        # Check if S4 object already exists from design matrix step
        # Both TMT and LFQ are protein-level workflows that use "protein_s4_initial" state
        s4_exists <- tryCatch({
          s4_obj <- workflow_data$state_manager$getState("protein_s4_initial")
          exists <- !is.null(s4_obj)
          message(sprintf("   DEBUG66: %s - Checking for protein_s4_initial: exists = %s", workflow_type, exists))
          if (exists) {
            message(sprintf("   DEBUG66: %s - S4 object class = %s", workflow_type, class(s4_obj)))
          }
          exists
        }, error = function(e) {
          message(sprintf("   DEBUG66: %s - Error checking S4: %s", workflow_type, e$message))
          FALSE
        })
        
        message(sprintf("   DEBUG66: %s - s4_exists final = %s", workflow_type, s4_exists))
        
        if (s4_exists) {
          # S4 already created in design matrix step, skip to protein QC
          logger::log_info("QC UI: S4 object exists from design matrix. Skipping S4 creation UI.")
          message(sprintf("   DEBUG66: %s - Rendering single Protein QC tab", workflow_type))
          shiny::tabsetPanel(
            id = ns("qc_tabs_tmt"),
            # Using mod_prot_qc_protein_ui (assuming it will be available)
            shiny::tabPanel("Protein QC", mod_prot_qc_protein_ui(ns("protein_qc"), workflow_type = workflow_type))
          )
        } else {
          # Fallback: show S4 creation UI if somehow not created
          logger::log_info("QC UI: S4 object not found. Showing S4 creation UI.")
          message(sprintf("   DEBUG66: %s - Rendering S4 creation + Protein QC tabs", workflow_type))
          shiny::tabsetPanel(
            id = ns("qc_tabs_tmt"),
            # Using mod_prot_qc_protein_s4_ui and mod_prot_qc_protein_ui
            shiny::tabPanel("Initial Protein Processing", mod_prot_qc_protein_s4_ui(ns("create_s4_tmt"))),
            shiny::tabPanel("Protein QC", mod_prot_qc_protein_ui(ns("protein_qc"), workflow_type = workflow_type))
          )
        }
      } else { # DIA (peptide-level workflow)
        message(sprintf("   DEBUG66: %s - Rendering Peptide + Protein QC tabs", workflow_type))
        # Define and return UI for peptide-based workflows
        shiny::tabsetPanel(
          id = ns("qc_tabs_lfq"),
          # Using mod_prot_qc_peptide_ui and mod_prot_qc_protein_ui
          shiny::tabPanel("Peptide QC", mod_prot_qc_peptide_ui(ns("peptide_qc"))),
          shiny::tabPanel("Protein QC", mod_prot_qc_protein_ui(ns("protein_qc"), workflow_type = workflow_type))
        )
      }
    })
    
    # == Component Server Calls ==========================================
    
    # Call the appropriate sub-applet orchestrators based on workflow type
    # This now waits for the qc_trigger from the design matrix step.
    shiny::observeEvent(qc_trigger(), {
      message("=== DEBUG66: mod_prot_qc_server observeEvent(qc_trigger) FIRED ===")
      message(sprintf("   DEBUG66: qc_trigger() value = %s", qc_trigger()))
      
      shiny::req(qc_trigger() == TRUE) # Only run when trigger is explicitly set to TRUE
      message("   DEBUG66: qc_trigger is TRUE, proceeding with module activation")
      
      workflow_type <- shiny::isolate(workflow_data$state_manager$workflow_type)
      logger::log_info(paste("Main QC Applet: Routing for workflow type:", workflow_type))
      message(sprintf("   DEBUG66: Isolated workflow_type = %s", workflow_type))
      
      if (workflow_type == "DIA") {
        message(sprintf("   DEBUG66: %s workflow - calling peptide and protein servers", workflow_type))
        # The peptide applet orchestrator will handle all peptide-level steps
        # Using mod_prot_qc_peptide_server
        mod_prot_qc_peptide_server("peptide_qc", workflow_data, experiment_paths, omic_type, experiment_label)
        message("   DEBUG66: mod_prot_qc_peptide_server called")
        
        # The protein applet orchestrator will handle all protein-level steps
        # Using mod_prot_qc_protein_server
        mod_prot_qc_protein_server("protein_qc", workflow_data, experiment_paths, omic_type, experiment_label)
        message("   DEBUG66: mod_prot_qc_protein_server called")
        
      } else if (workflow_type %in% c("TMT", "LFQ")) {
        message(sprintf("   DEBUG66: %s workflow - checking S4 existence", workflow_type))
        # Check if S4 object already exists from design matrix step
        # Both TMT and LFQ are protein-level workflows that use "protein_s4_initial" state
        s4_exists <- tryCatch({
          s4_obj <- workflow_data$state_manager$getState("protein_s4_initial")
          exists <- !is.null(s4_obj)
          message(sprintf("   DEBUG66: %s observeEvent - S4 exists = %s", workflow_type, exists))
          exists
        }, error = function(e) {
          message(sprintf("   DEBUG66: %s observeEvent - Error checking S4: %s", workflow_type, e$message))
          FALSE
        })
        
        if (s4_exists) {
          # S4 already created in design matrix step, go straight to protein QC
          logger::log_info(sprintf("%s: S4 object exists from design matrix. Skipping S4 creation, starting protein QC.", workflow_type))
          message(sprintf("   DEBUG66: %s - S4 exists, calling mod_prot_qc_protein_server only", workflow_type))
          mod_prot_qc_protein_server("protein_qc", workflow_data, experiment_paths, omic_type, experiment_label)
          message(sprintf("   DEBUG66: %s - mod_prot_qc_protein_server called", workflow_type))
        } else {
          # S4 doesn't exist (shouldn't happen in normal flow), use the creation module
          logger::log_warn(sprintf("%s: S4 object not found. Using S4 creation module.", workflow_type))
          message(sprintf("   DEBUG66: %s - S4 MISSING, calling create_s4 + protein servers", workflow_type))
          # Using mod_prot_qc_protein_s4_server
          mod_prot_qc_protein_s4_server("create_s4_tmt", workflow_data, omic_type, experiment_label)
          message(sprintf("   DEBUG66: %s - mod_prot_qc_protein_s4_server called", workflow_type))
          mod_prot_qc_protein_server("protein_qc", workflow_data, experiment_paths, omic_type, experiment_label)
          message(sprintf("   DEBUG66: %s - mod_prot_qc_protein_server called", workflow_type))
        }
      }
      message("=== DEBUG66: mod_prot_qc_server observeEvent complete ===")
    }, ignoreNULL = TRUE, once = TRUE) # Run only once when triggered
    
    # Final step: update tab status when the final state is reached
    shiny::observeEvent(workflow_data$state_manager$states$protein_replicate_filtered, {
      # This observer is triggered when the final protein state is reached.
      # It can be used to update the UI or perform final actions.
      # For now, we'll just log it.
      logger::log_info("Final protein state reached. Quality Control module is ready.")
    })
    
  })
}

