#' @title Create Protein S4 Module
#'
#' @description UI component for creating a ProteinQuantitativeData S4 object
#' from protein-level data. This provides a user-facing button to trigger
#' the process for workflows like TMT.
#'
#' @name mod_prot_qc_protein_s4
NULL

#' @rdname mod_prot_qc_protein_s4
#' @export
#' @importFrom shiny NS tagList fluidPage wellPanel h4 p actionButton hr h5 verbatimTextOutput
mod_prot_qc_protein_s4_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tagList(
    shiny::fluidPage(
    shiny::wellPanel(
      shiny::h4("Finalise Protein Data"),
      shiny::p("This workflow starts with protein-level data. Click the button below to format the data into an S4 object for downstream analysis."),
      shiny::actionButton(ns("create_protein_s4"), "Create Protein S4 Object", class = "btn-primary"),
      shiny::hr(),
      shiny::h5("Results"),
      shiny::verbatimTextOutput(ns("s4_creation_results"))
    )
  )
  )
}

#' @rdname mod_prot_qc_protein_s4
#' @export
#' @importFrom shiny moduleServer observeEvent req showNotification removeNotification renderText
#' @importFrom logger log_info log_error
mod_prot_qc_protein_s4_server <- function(id, workflow_data, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    # This module is triggered by an action button in the UI.
    shiny::observeEvent(input$create_protein_s4, {
      shiny::req(
        workflow_data$data_cln,
        workflow_data$design_matrix,
        workflow_data$column_mapping,
        workflow_data$config_list
      )
      
      shiny::showNotification("Creating Protein S4 object from imported data...", id = "s4_creation_working", duration = NULL)
      
      tryCatch({
        
        logger::log_info("Protein S4 Creation: Starting process for protein-level workflow (e.g., TMT)")
        
        protein_id_col <- workflow_data$column_mapping$protein_col
        
        # Ensure the protein ID column exists
        if (!protein_id_col %in% names(workflow_data$data_cln)) {
          stop(paste("Protein ID column", protein_id_col, "not found in the data."))
        }
        
        # Create ProteinQuantitativeData S4 object directly from workflow data
        protein_obj <- ProteinQuantitativeData(
          protein_quant_table = workflow_data$data_cln,
          protein_id_column = protein_id_col,
          protein_id_table = workflow_data$data_cln |> dplyr::distinct(!!rlang::sym(protein_id_col)),
          design_matrix = workflow_data$design_matrix,
          sample_id = "Run",
          group_id = "group",
          technical_replicate_id = "replicates",
          args = workflow_data$config_list 
        )
        
        # Save the new S4 object to the state manager
        workflow_data$state_manager$saveState(
          state_name = "protein_s4_initial",
          s4_data_object = protein_obj,
          config_object = list(
            s4_class = "ProteinQuantitativeData",
            protein_id_column = protein_id_col
          ),
          description = "Created initial ProteinQuantitativeData S4 object from protein-level data."
        )
        
        # Generate summary
        protein_count <- protein_obj@protein_quant_table |>
          dplyr::distinct(!!rlang::sym(protein_id_col)) |>
          nrow()
        
        result_text <- paste(
          "Protein S4 Object Created Successfully\n",
          "======================================\n",
          sprintf("Proteins loaded: %d\n", protein_count),
          sprintf("S4 Class: %s\n", class(protein_obj)[1]),
          "State saved as: 'protein_s4_initial'\n",
          "\nReady for protein accession cleanup."
        )
        
        output$s4_creation_results <- shiny::renderText(result_text)
        
        logger::log_info("Protein S4 object creation from protein data completed successfully")
        shiny::removeNotification("s4_creation_working")
        shiny::showNotification("Protein S4 object created successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error creating Protein S4 object:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error", duration = 15)
        shiny::removeNotification("s4_creation_working")
      })
    })
    
    # Revert S4 Creation (if needed)
    # Note: The original code had this, although there might not be a previous state to revert to if this is the first step.
    # We'll keep it but make it revert to initial.
    shiny::observeEvent(input$revert_s4_creation, {
      tryCatch({
        # As a simple revert, we can go back to the 'initial' empty state
        workflow_data$state_manager$revertToState("initial")
        
        output$s4_creation_results <- shiny::renderText("Reverted to initial empty state. You may need to re-run previous steps.")
        
        logger::log_info("Reverted S4 object creation.")
        shiny::showNotification("Reverted successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error reverting:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error")
      })
    })
    
  })
}

