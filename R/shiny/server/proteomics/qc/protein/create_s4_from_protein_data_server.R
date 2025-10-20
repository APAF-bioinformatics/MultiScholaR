#' @title create_s4_from_protein_data_server
#'
#' @description Server module for creating a ProteinQuantitativeData S4 object directly
#' from protein-level data. This is used for workflows like TMT that bypass
#' the peptide-to-protein rollup stage.
#'
#' @param id Module ID
#' @param workflow_data A reactive values object to store workflow data.
#' @param omic_type The omics type (e.g., "proteomics").
#' @param experiment_label The experiment label.
#'
#' @export
#' @import shiny
#' @importFrom logger log_info log_error
create_s4_from_protein_data_server <- function(id, workflow_data, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    # This module is triggered by an action button in the UI.
    observeEvent(input$create_protein_s4, {
      shiny::req(
        workflow_data$data_cln,
        workflow_data$design_matrix,
        workflow_data$column_mapping,
        workflow_data$config_list
      )
      
      shiny::showNotification("Creating Protein S4 object from imported data...", id = "s4_creation_working", duration = NULL)
      
      tryCatch({
        
        log_info("Protein S4 Creation: Starting process for protein-level workflow (e.g., TMT)")
        
        protein_id_col <- workflow_data$column_mapping$protein_col
        
        # Ensure the protein ID column exists
        if (!protein_id_col %in% names(workflow_data$data_cln)) {
          stop(paste("Protein ID column", protein_id_col, "not found in the data."))
        }
        
        # Create ProteinQuantitativeData S4 object directly from workflow data
        protein_obj <- ProteinQuantitativeData(
          protein_quant_table = workflow_data$data_cln,
          protein_id_column = protein_id_col,
          protein_id_table = workflow_data$data_cln |> dplyr::distinct(!!sym(protein_id_col)),
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
          dplyr::distinct(!!sym(protein_id_col)) |>
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
    
    # Revert S4 Creation
    observeEvent(input$revert_s4_creation, {
      tryCatch({
        # There's no previous S4 state to revert to, so this action might clear the state
        # For now, we'll just inform the user. A more robust implementation could
        # revert workflow_data itself.
        
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
