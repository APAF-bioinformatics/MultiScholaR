#' @title Protein Duplicate Removal Module
#'
#' @description A Shiny module for removing duplicate proteins.
#'
#' @name mod_prot_qc_protein_dedup
NULL

#' @rdname mod_prot_qc_protein_dedup
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 p hr selectInput helpText div actionButton verbatimTextOutput plotOutput
mod_prot_qc_protein_dedup_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tabPanel(
    "Duplicate Removal",
    shiny::br(),
    shiny::fluidRow(
      shiny::column(4,
        shiny::wellPanel(
          shiny::h4("Remove Duplicate Proteins"),
          shiny::p("Aggregate duplicate protein entries by taking the mean across matching proteins."),
          shiny::hr(),
          
          shiny::selectInput(ns("duplicate_aggregation_method"), 
            "Aggregation Method", 
            choices = c("mean", "median", "max"),
            selected = "mean",
            width = "100%"
          ),
          shiny::helpText("Method for combining duplicate protein measurements (default: mean)"),
          
          shiny::hr(),
          shiny::div(
            shiny::actionButton(ns("apply_duplicate_removal"), "Remove Duplicates", 
              class = "btn-primary", width = "48%"),
            shiny::actionButton(ns("revert_duplicate_removal"), "Revert", 
              class = "btn-warning", width = "48%", style = "margin-left: 4%")
          )
        )
      ),
      shiny::column(8,
        shiny::verbatimTextOutput(ns("duplicate_removal_results")),
        shiny::br(),
        shinyjqui::jqui_resizable(
          shiny::plotOutput(ns("duplicate_removal_plot"), height = "800px", width = "100%")
        )
      )
    )
  )
}

#' @rdname mod_prot_qc_protein_dedup
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification removeNotification renderText renderPlot
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
mod_prot_qc_protein_dedup_server <- function(id, workflow_data, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    duplicate_removal_plot <- shiny::reactiveVal(NULL)
    
    # Step 4: Duplicate Protein Removal (chunk 22)
    shiny::observeEvent(input$apply_duplicate_removal, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Removing duplicate proteins...", id = "duplicate_removal_working", duration = NULL)
      
      tryCatch({
        # Get current ProteinQuantitativeData S4 object from the active state
        current_s4 <- workflow_data$state_manager$getState()
        shiny::req(current_s4)
        
        logger::log_info("Protein Processing: Removing duplicate proteins using {input$duplicate_aggregation_method}")
        
        # Identify duplicates first
        duplicates <- current_s4@protein_quant_table |>
          dplyr::group_by(Protein.Ids) |>
          dplyr::filter(dplyr::n() > 1) |>
          dplyr::select(Protein.Ids) |>
          dplyr::distinct() |>
          dplyr::pull(Protein.Ids)
        
        # Apply duplicate removal
        current_s4@protein_quant_table <- current_s4@protein_quant_table |>
          dplyr::group_by(Protein.Ids) |>
          dplyr::summarise(
            dplyr::across(dplyr::matches("\\d+"), ~ get(input$duplicate_aggregation_method)(.x, na.rm = TRUE))
          ) |>
          dplyr::ungroup()
        
        # Save new state
        workflow_data$state_manager$saveState(
          state_name = "duplicates_removed",
          s4_data_object = current_s4,
          config_object = list(
            aggregation_method = input$duplicate_aggregation_method,
            duplicates_found = duplicates,
            num_duplicates = length(duplicates)
          ),
          description = "Removed duplicate proteins by aggregation"
        )
        
        # Generate summary
        protein_count <- current_s4@protein_quant_table |>
          dplyr::distinct(Protein.Ids) |>
          nrow()
        
        result_text <- paste(
          "Duplicate Protein Removal Completed Successfully\n",
          "===============================================\n",
          sprintf("Proteins remaining: %d\n", protein_count),
          sprintf("Duplicates found: %d\n", length(duplicates)),
          sprintf("Aggregation method: %s\n", input$duplicate_aggregation_method),
          "State saved as: 'duplicates_removed'\n"
        )
        
        output$duplicate_removal_results <- shiny::renderText(result_text)
        
        # Update filtering visualization and capture plot
        plot_grid <- updateProteinFiltering(
          data = current_s4@protein_quant_table,
          step_name = "12_duplicates_removed",
          omic_type = omic_type,
          experiment_label = experiment_label,
          return_grid = TRUE,
          overwrite = TRUE
        )
        duplicate_removal_plot(plot_grid)
        
        logger::log_info("Duplicate protein removal completed successfully")
        shiny::removeNotification("duplicate_removal_working")
        shiny::showNotification("Duplicate protein removal completed successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error removing duplicate proteins:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error", duration = 15)
        shiny::removeNotification("duplicate_removal_working")
      })
    })
    
    # Revert Duplicate Removal
    shiny::observeEvent(input$revert_duplicate_removal, {
      tryCatch({
        # Revert to the previous state in the history
        history <- workflow_data$state_manager$getHistory()
        if (length(history) > 1) {
          prev_state_name <- history[length(history) - 1]
          reverted_s4 <- workflow_data$state_manager$revertToState(prev_state_name)
          output$duplicate_removal_results <- shiny::renderText(paste("Reverted to previous state:", prev_state_name))
          logger::log_info(paste("Reverted duplicate removal to", prev_state_name))
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
    
    # Render duplicate removal plot
    output$duplicate_removal_plot <- shiny::renderPlot({
      shiny::req(duplicate_removal_plot())
      grid::grid.draw(duplicate_removal_plot())
    })
  })
}

