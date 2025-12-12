#' @title Protein Peptide Count Filter Module
#'
#' @description A Shiny module for applying the protein peptide count filter.
#'
#' @name mod_prot_qc_peptide_protein
NULL

#' @rdname mod_prot_qc_peptide_protein
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 p hr numericInput helpText div actionButton verbatimTextOutput plotOutput
mod_prot_qc_peptide_protein_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tabPanel(
    "Protein Peptides",
    shiny::br(),
    shiny::fluidRow(
      shiny::column(4,
        shiny::wellPanel(
          shiny::h4("Minimum Peptides per Protein"),
          shiny::p("Keep proteins only if they have sufficient peptide evidence (two-peptide rule)."),
          shiny::hr(),
          
          shiny::numericInput(ns("min_peptides_per_protein"), 
            "Min Peptides per Protein", 
            value = 2, min = 1, max = 5, step = 1,
            width = "100%"
          ),
          shiny::helpText("Higher = requires more unique peptides per protein (default: 2)"),
          
          shiny::numericInput(ns("min_peptidoforms_per_protein"), 
            "Min Peptidoforms per Protein", 
            value = 2, min = 1, max = 5, step = 1,
            width = "100%"
          ),
          shiny::helpText("Higher = requires more peptide forms per protein (default: 2)"),
          
          shiny::hr(),
          shiny::div(
            shiny::actionButton(ns("apply_protein_peptide_filter"), "Apply Filter", 
              class = "btn-primary", width = "48%"),
            shiny::actionButton(ns("revert_protein_peptide"), "Revert", 
              class = "btn-warning", width = "48%", style = "margin-left: 4%")
          )
        )
      ),
      shiny::column(8,
        shiny::verbatimTextOutput(ns("protein_peptide_results")),
        shiny::br(),
        shinyjqui::jqui_resizable(
          shiny::plotOutput(ns("protein_peptide_plot"), height = "800px", width = "100%")
        )
      )
    )
  )
}

#' @rdname mod_prot_qc_peptide_protein
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification removeNotification renderText renderPlot
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
mod_prot_qc_peptide_protein_server <- function(id, workflow_data, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    protein_peptide_plot <- shiny::reactiveVal(NULL)
    
    # Step 4: Protein Peptide Count Filter (chunk 13)
    shiny::observeEvent(input$apply_protein_peptide_filter, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying protein peptide count filter...", id = "protein_peptide_working", duration = NULL)
      
      tryCatch({
        # Get current S4 object from the active state
        current_s4 <- workflow_data$state_manager$getState()
        shiny::req(current_s4)
        
        logger::log_info(sprintf("QC Step: Applying protein peptide count filter (min: %s)", input$min_peptides_per_protein))
        
        # ✅ FIXED: Use updateConfigParameter to sync S4 object AND global config_list
        # Use config.ini parameter names, not function parameter names
        current_s4 <- updateConfigParameter(
          theObject = current_s4,
          function_name = "filterMinNumPeptidesPerProtein",
          parameter_name = "peptides_per_protein_cutoff",
          new_value = input$min_peptides_per_protein
        )
        
        current_s4 <- updateConfigParameter(
          theObject = current_s4,
          function_name = "filterMinNumPeptidesPerProtein",
          parameter_name = "peptidoforms_per_protein_cutoff",
          new_value = input$min_peptidoforms_per_protein
        )
        
        # Apply S4 transformation (EXISTING S4 CODE - UNCHANGED)
        filtered_s4 <- filterMinNumPeptidesPerProtein(theObject = current_s4)
        
        # Track QC parameters in workflow_data
        if (is.null(workflow_data$qc_params)) {
          workflow_data$qc_params <- list()
        }
        if (is.null(workflow_data$qc_params$peptide_qc)) {
          workflow_data$qc_params$peptide_qc <- list()
        }
        
        workflow_data$qc_params$peptide_qc$protein_peptide_filter <- list(
          min_peptides_per_protein = input$min_peptides_per_protein,
          min_peptidoforms_per_protein = input$min_peptidoforms_per_protein,
          timestamp = Sys.time()
        )
        
        # Save new state in R6 manager
        workflow_data$state_manager$saveState(
          state_name = "protein_peptide_filtered",
          s4_data_object = filtered_s4,
          config_object = list(
            min_peptides_per_protein = input$min_peptides_per_protein,
            min_peptidoforms_per_protein = input$min_peptidoforms_per_protein
          ),
          description = "Applied minimum peptides per protein filter"
        )
        
        # Generate summary
        protein_count <- filtered_s4@peptide_data |>
          dplyr::distinct(Protein.Ids) |>
          nrow()
        
        result_text <- paste(
          "Protein Peptide Count Filter Applied Successfully\n",
          "===============================================\n",
          sprintf("Proteins remaining: %d\n", protein_count),
          sprintf("Min peptides per protein: %d\n", input$min_peptides_per_protein),
          sprintf("Min peptidoforms per protein: %d\n", input$min_peptidoforms_per_protein),
          "State saved as: 'protein_peptide_filtered'\n"
        )
        
        output$protein_peptide_results <- shiny::renderText(result_text)
        
        # Update filtering visualization and capture plot
        plot_grid <- updateProteinFiltering(
          data = filtered_s4@peptide_data,
          step_name = "5_protein_peptide_filtered",
          omic_type = omic_type,
          experiment_label = experiment_label,
          return_grid = TRUE,
          overwrite = TRUE
        )
        protein_peptide_plot(plot_grid)
        
        logger::log_info("Protein peptide count filter applied successfully")
        shiny::removeNotification("protein_peptide_working")
        shiny::showNotification("Protein peptide count filter applied successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error applying protein peptide count filter:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error", duration = 15)
        shiny::removeNotification("protein_peptide_working")
      })
    })
    
    # Revert Protein Peptide Filter
    shiny::observeEvent(input$revert_protein_peptide, {
      tryCatch({
        # Revert to the previous state in the history
        history <- workflow_data$state_manager$getHistory()
        if (length(history) > 1) {
          prev_state_name <- history[length(history) - 1]
          reverted_s4 <- workflow_data$state_manager$revertToState(prev_state_name)
          output$protein_peptide_results <- shiny::renderText(paste("Reverted to previous state:", prev_state_name))
          logger::log_info(paste("Reverted protein peptide filter to", prev_state_name))
          shiny::showNotification("Reverted successfully", type = "message")
        } else {
          stop("No previous state to revert to.")
        }
        
      }, error = function(e) {
        msg <- paste("Error reverting:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error")
      })
    })
    
    # Render protein peptide filter plot
    output$protein_peptide_plot <- shiny::renderPlot({
      shiny::req(protein_peptide_plot())
      grid::grid.draw(protein_peptide_plot())
    })
  })
}

