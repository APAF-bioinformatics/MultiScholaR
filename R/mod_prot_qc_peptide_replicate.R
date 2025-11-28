#' @title Peptide Replicate Filter Module
#'
#' @description A Shiny module for applying the replicate filter.
#'
#' @name mod_prot_qc_peptide_replicate
NULL

#' @rdname mod_prot_qc_peptide_replicate
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 p hr textInput helpText div actionButton verbatimTextOutput plotOutput
mod_prot_qc_peptide_replicate_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tabPanel(
    "Replicate Filter",
    shiny::br(),
    shiny::fluidRow(
      shiny::column(4,
        shiny::wellPanel(
          shiny::h4("Remove Single-Replicate Peptides"),
          shiny::p("Remove peptides that appear in only one replicate across all groups."),
          shiny::hr(),
          
          shiny::textInput(ns("replicate_group_column"), 
            "Replicate Group Column", 
            value = "replicates",
            width = "100%"
          ),
          shiny::helpText("Column name for grouping replicates (default: 'replicates')"),
          
          shiny::hr(),
          shiny::div(
            shiny::actionButton(ns("apply_replicate_filter"), "Apply Filter", 
              class = "btn-primary", width = "48%"),
            shiny::actionButton(ns("revert_replicate"), "Revert", 
              class = "btn-warning", width = "48%", style = "margin-left: 4%")
          )
        )
      ),
      shiny::column(8,
        shiny::verbatimTextOutput(ns("replicate_results")),
        shiny::br(),
        shinyjqui::jqui_resizable(
          shiny::plotOutput(ns("replicate_plot"), height = "800px", width = "100%")
        )
      )
    )
  )
}

#' @rdname mod_prot_qc_peptide_replicate
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification removeNotification renderText renderPlot
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
mod_prot_qc_peptide_replicate_server <- function(id, workflow_data, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    replicate_plot <- shiny::reactiveVal(NULL)
    
    # Step 6: Replicate Filter (chunk 15)
    shiny::observeEvent(input$apply_replicate_filter, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying replicate filter...", id = "replicate_working", duration = NULL)
      
      tryCatch({
        # Get current S4 object from the active state
        current_s4 <- workflow_data$state_manager$getState()
        shiny::req(current_s4)
        
        logger::log_info(sprintf("QC Step: Applying replicate filter (column: %s)", input$replicate_group_column))
        
        # Apply S4 transformation (EXISTING S4 CODE - UNCHANGED)
        # Note: This function takes the replicate column as a parameter
        filtered_s4 <- removePeptidesWithOnlyOneReplicate(
          current_s4,
          replicate_group_column = input$replicate_group_column
        )
        
        # Track QC parameters in workflow_data
        if (is.null(workflow_data$qc_params)) {
          workflow_data$qc_params <- list()
        }
        if (is.null(workflow_data$qc_params$peptide_qc)) {
          workflow_data$qc_params$peptide_qc <- list()
        }
        
        workflow_data$qc_params$peptide_qc$replicate_filter <- list(
          replicate_group_column = input$replicate_group_column,
          timestamp = Sys.time()
        )
        
        # Save new state in R6 manager
        workflow_data$state_manager$saveState(
          state_name = "replicate_filtered",
          s4_data_object = filtered_s4,
          config_object = list(
            replicate_group_column = input$replicate_group_column
          ),
          description = "Applied replicate filter (removed single-replicate peptides)"
        )
        
        # Generate summary
        protein_count <- filtered_s4@peptide_data |>
          dplyr::distinct(Protein.Ids) |>
          nrow()
        
        result_text <- paste(
          "Replicate Filter Applied Successfully\n",
          "====================================\n",
          sprintf("Proteins remaining: %d\n", protein_count),
          sprintf("Replicate group column: %s\n", input$replicate_group_column),
          "State saved as: 'replicate_filtered'\n"
        )
        
        output$replicate_results <- shiny::renderText(result_text)
        
        # Update filtering visualization and capture plot
        plot_grid <- updateProteinFiltering(
          data = filtered_s4@peptide_data,
          step_name = "7_replicate_filtered",
          omic_type = omic_type,
          experiment_label = experiment_label,
          return_grid = TRUE,
          overwrite = TRUE
        )
        replicate_plot(plot_grid)
        
        logger::log_info("Replicate filter applied successfully")
        shiny::removeNotification("replicate_working")
        shiny::showNotification("Replicate filter applied successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error applying replicate filter:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error", duration = 15)
        shiny::removeNotification("replicate_working")
      })
    })
    
    # Revert Replicate Filter
    shiny::observeEvent(input$revert_replicate, {
      tryCatch({
        # Revert to the previous state in the history
        history <- workflow_data$state_manager$getHistory()
        if (length(history) > 1) {
          prev_state_name <- history[length(history) - 1]
          reverted_s4 <- workflow_data$state_manager$revertToState(prev_state_name)
          output$replicate_results <- shiny::renderText(paste("Reverted to previous state:", prev_state_name))
          logger::log_info(paste("Reverted replicate filter to", prev_state_name))
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
    
    # Render replicate filter plot
    output$replicate_plot <- shiny::renderPlot({
      shiny::req(replicate_plot())
      grid::grid.draw(replicate_plot())
    })
  })
}

