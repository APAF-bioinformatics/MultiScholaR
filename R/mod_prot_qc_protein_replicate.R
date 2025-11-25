#' @title Protein Replicate Filter Module
#'
#' @description A Shiny module for applying the protein replicate filter.
#'
#' @name mod_prot_qc_protein_replicate
NULL

#' @rdname mod_prot_qc_protein_replicate
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 p hr textInput helpText numericInput div actionButton verbatimTextOutput plotOutput
mod_prot_qc_protein_replicate_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tagList(
    shiny::tabPanel(
    "Protein Replicate Filter",
    shiny::br(),
    shiny::fluidRow(
      shiny::column(4,
        shiny::wellPanel(
          shiny::h4("Remove Single-Replicate Proteins"),
          shiny::p("Remove proteins detected in only one replicate across experimental groups."),
          shiny::hr(),
          
          shiny::textInput(ns("protein_grouping_variable"), 
            "Grouping Variable", 
            value = "group",
            width = "100%"
          ),
          shiny::helpText("Column name for experimental grouping (default: 'group')"),
          
          shiny::numericInput(ns("parallel_cores"), 
            "Parallel Processing Cores", 
            value = 4, min = 1, max = 8, step = 1,
            width = "100%"
          ),
          shiny::helpText("Number of cores for parallel processing (default: 4)"),
          
          shiny::hr(),
          shiny::div(
            shiny::actionButton(ns("apply_protein_replicate_filter"), "Apply Filter", 
              class = "btn-primary", width = "48%"),
            shiny::actionButton(ns("revert_protein_replicate_filter"), "Revert", 
              class = "btn-warning", width = "48%", style = "margin-left: 4%")
          )
        )
      ),
      shiny::column(8,
        shiny::verbatimTextOutput(ns("protein_replicate_filter_results")),
        shiny::br(),
        shinyjqui::jqui_resizable(
          shiny::plotOutput(ns("protein_replicate_filter_plot"), height = "800px", width = "100%")
        )
      )
    )
  )
  )
}

#' @rdname mod_prot_qc_protein_replicate
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification removeNotification renderText renderPlot
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
mod_prot_qc_protein_replicate_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    protein_replicate_filter_plot <- shiny::reactiveVal(NULL)
    
    # Step 5: Protein Replicate Filter (chunk 23)
    shiny::observeEvent(input$apply_protein_replicate_filter, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying protein replicate filter...", id = "protein_replicate_filter_working", duration = NULL)
      
      tryCatch({
        # Get current ProteinQuantitativeData S4 object from the active state
        current_s4 <- workflow_data$state_manager$getState()
        shiny::req(current_s4)
        
        logger::log_info("Protein Processing: Applying protein replicate filter with {input$parallel_cores} cores")
        
        # Set up parallel processing
        # Note: The 'new_cluster' function must be available in the environment
        core_utilisation <- if(exists("new_cluster")) new_cluster(input$parallel_cores) else NULL
        
        # Apply S4 transformation (EXISTING S4 CODE - UNCHANGED)
        filtered_s4 <- removeProteinsWithOnlyOneReplicate(
          current_s4,
          core_utilisation,
          grouping_variable = input$protein_grouping_variable
        )
        
        # Save filtered data to file (as in original workflow)
        if (!is.null(experiment_paths$protein_qc_dir)) {
          output_file <- file.path(experiment_paths$protein_qc_dir, "remove_proteins_with_only_one_rep.tsv")
          vroom::vroom_write(filtered_s4@protein_quant_table, output_file)
        } else {
          output_file <- "remove_proteins_with_only_one_rep.tsv"
        }
        
        # Track QC parameters in workflow_data
        if (is.null(workflow_data$qc_params)) {
          workflow_data$qc_params <- list()
        }
        if (is.null(workflow_data$qc_params$protein_qc)) {
          workflow_data$qc_params$protein_qc <- list()
        }
        
        workflow_data$qc_params$protein_qc$replicate_filter <- list(
          grouping_variable = input$protein_grouping_variable,
          parallel_cores = input$parallel_cores,
          timestamp = Sys.time()
        )
        
        # Save QC parameters to file for persistence
        tryCatch({
          if (exists("experiment_paths") && !is.null(experiment_paths$source_dir)) {
            qc_params_file <- file.path(experiment_paths$source_dir, "qc_params.RDS")
            saveRDS(workflow_data$qc_params, qc_params_file)
            logger::log_info(sprintf("Saved QC parameters to: %s", qc_params_file))
          }
        }, error = function(e) {
          logger::log_warn(sprintf("Could not save QC parameters file: %s", e$message))
        })
        
        # Save new state
        workflow_data$state_manager$saveState(
          state_name = "protein_replicate_filtered",
          s4_data_object = filtered_s4,
          config_object = list(
            grouping_variable = input$protein_grouping_variable,
            parallel_cores = input$parallel_cores,
            output_file = output_file
          ),
          description = "Applied protein replicate filter (removed single-replicate proteins)"
        )
        
        # Generate summary
        protein_count <- filtered_s4@protein_quant_table |>
          dplyr::distinct(Protein.Ids) |>
          nrow()
        
        # Track protein count after QC filtering (final QC step before normalization)
        if (is.null(workflow_data$protein_counts)) {
          workflow_data$protein_counts <- list()
        }
        workflow_data$protein_counts$after_qc_filtering <- protein_count
        logger::log_info(sprintf("Tracked protein count after QC filtering: %d", protein_count))
        
        result_text <- paste(
          "Protein Replicate Filter Applied Successfully\n",
          "============================================\n",
          sprintf("Proteins remaining: %d\n", protein_count),
          sprintf("Grouping variable: %s\n", input$protein_grouping_variable),
          sprintf("Parallel cores used: %d\n", input$parallel_cores),
          sprintf("Output file: %s\n", basename(output_file)),
          "State saved as: 'protein_replicate_filtered'\n",
          "\nProtein filtering pipeline complete!"
        )
        
        output$protein_replicate_filter_results <- shiny::renderText(result_text)
        
        # Update filtering visualization and capture plot
        plot_grid <- updateProteinFiltering(
          data = filtered_s4@protein_quant_table,
          step_name = "13_protein_replicate_filtered",
          omic_type = omic_type,
          experiment_label = experiment_label,
          return_grid = TRUE,
          overwrite = TRUE
        )
        protein_replicate_filter_plot(plot_grid)
        
        logger::log_info("Protein replicate filter applied successfully")
        shiny::removeNotification("protein_replicate_filter_working")
        shiny::showNotification("Protein replicate filter applied successfully. Pipeline complete!", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error applying protein replicate filter:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error", duration = 15)
        shiny::removeNotification("protein_replicate_filter_working")
      })
    })
    
    # Revert Protein Replicate Filter
    shiny::observeEvent(input$revert_protein_replicate_filter, {
      tryCatch({
        # Revert to the previous state in the history
        history <- workflow_data$state_manager$getHistory()
        if (length(history) > 1) {
          prev_state_name <- history[length(history) - 1]
          reverted_s4 <- workflow_data$state_manager$revertToState(prev_state_name)
          output$protein_replicate_filter_results <- shiny::renderText(paste("Reverted to previous state:", prev_state_name))
          logger::log_info(paste("Reverted protein replicate filter to", prev_state_name))
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
    
    # Render protein replicate filter plot
    output$protein_replicate_filter_plot <- shiny::renderPlot({
      shiny::req(protein_replicate_filter_plot())
      grid::grid.draw(protein_replicate_filter_plot())
    })
  })
}

