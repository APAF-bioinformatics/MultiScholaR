#' @title Peptide Sample Filter Module
#'
#' @description A Shiny module for applying the sample quality filter.
#'
#' @name mod_prot_qc_peptide_sample
#' @export
NULL

#' @rdname mod_prot_qc_peptide_sample
#' @export
#' @import shiny
#' @import shinydashboard
mod_prot_qc_peptide_sample_ui <- function(id) {
  ns <- NS(id)
  
  shiny::tabPanel(
    "Sample Quality",
    shiny::br(),
    shiny::fluidRow(
      shiny::column(4,
        shiny::wellPanel(
          shiny::h4("Minimum Peptides per Sample"),
          shiny::p("Remove samples with insufficient peptide counts (poor sample performance)."),
          shiny::hr(),
          
          shiny::numericInput(ns("min_peptides_per_sample"), 
            "Min Peptides per Sample", 
            value = 500, min = 100, max = 5000, step = 100,
            width = "100%"
          ),
          shiny::helpText("Higher = stricter sample quality filter (default: 500)"),
          
          shiny::hr(),
          shiny::div(
            shiny::actionButton(ns("apply_sample_filter"), "Apply Filter", 
              class = "btn-primary", width = "48%"),
            shiny::actionButton(ns("revert_sample"), "Revert", 
              class = "btn-warning", width = "48%", style = "margin-left: 4%")
          )
        )
      ),
      shiny::column(8,
        shiny::verbatimTextOutput(ns("sample_results")),
        shiny::br(),
        shinyjqui::jqui_resizable(
          shiny::plotOutput(ns("sample_plot"), height = "800px", width = "100%")
        )
      )
    )
  )
}

#' @rdname mod_prot_qc_peptide_sample
#' @export
#' @import shiny
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
mod_prot_qc_peptide_sample_server <- function(id, workflow_data, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    sample_plot <- shiny::reactiveVal(NULL)
    
    # Step 5: Sample Quality Filter (chunk 14)
    shiny::observeEvent(input$apply_sample_filter, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying sample quality filter...", id = "sample_working", duration = NULL)
      
      tryCatch({
        # Get current S4 object from the active state
        current_s4 <- workflow_data$state_manager$getState()
        shiny::req(current_s4)
        
        logger::log_info(sprintf("QC Step: Applying sample quality filter (min: %s)", input$min_peptides_per_sample))
        
        # ✅ FIXED: Use updateConfigParameter to sync S4 object AND global config_list
        current_s4 <- updateConfigParameter(
          theObject = current_s4,
          function_name = "filterMinNumPeptidesPerSample",
          parameter_name = "peptides_per_sample_cutoff",
          new_value = input$min_peptides_per_sample
        )
        
        # Track samples before filtering
        samples_before <- current_s4@peptide_data |>
          dplyr::distinct(!!rlang::sym(current_s4@sample_id)) |>
          dplyr::pull(!!rlang::sym(current_s4@sample_id))
        
        # Apply S4 transformation (EXISTING S4 CODE - UNCHANGED)
        filtered_s4 <- filterMinNumPeptidesPerSample(theObject = current_s4)
        
        # Track samples after filtering
        samples_after <- filtered_s4@peptide_data |>
          dplyr::distinct(!!rlang::sym(filtered_s4@sample_id)) |>
          dplyr::pull(!!rlang::sym(filtered_s4@sample_id))
        
        # Identify removed samples
        samples_removed <- setdiff(samples_before, samples_after)
        samples_removed_count <- length(samples_removed)
        
        # Store removed samples info in S4 object @args for report generation
        if (is.null(filtered_s4@args$filterMinNumPeptidesPerSample)) {
          filtered_s4@args$filterMinNumPeptidesPerSample <- list()
        }
        filtered_s4@args$filterMinNumPeptidesPerSample$samples_removed <- samples_removed
        filtered_s4@args$filterMinNumPeptidesPerSample$samples_removed_count <- samples_removed_count
        filtered_s4@args$filterMinNumPeptidesPerSample$samples_before_count <- length(samples_before)
        filtered_s4@args$filterMinNumPeptidesPerSample$samples_after_count <- length(samples_after)
        
        # Track QC parameters in workflow_data
        if (is.null(workflow_data$qc_params)) {
          workflow_data$qc_params <- list()
        }
        if (is.null(workflow_data$qc_params$peptide_qc)) {
          workflow_data$qc_params$peptide_qc <- list()
        }
        
        workflow_data$qc_params$peptide_qc$sample_filter <- list(
          min_peptides_per_sample = input$min_peptides_per_sample,
          samples_removed = samples_removed,
          samples_removed_count = samples_removed_count,
          samples_before_count = length(samples_before),
          samples_after_count = length(samples_after),
          timestamp = Sys.time()
        )
        
        # Save new state in R6 manager
        workflow_data$state_manager$saveState(
          state_name = "sample_filtered",
          s4_data_object = filtered_s4,
          config_object = list(
            min_peptides_per_sample = input$min_peptides_per_sample,
            samples_removed = samples_removed,
            samples_removed_count = samples_removed_count
          ),
          description = "Applied minimum peptides per sample filter"
        )
        
        # Generate summary
        protein_count <- filtered_s4@peptide_data |>
          dplyr::distinct(Protein.Ids) |>
          nrow()
        run_count <- filtered_s4@peptide_data |>
          dplyr::distinct(Run) |>
          nrow()
        
        result_text <- paste(
          "Sample Quality Filter Applied Successfully\n",
          "=========================================\n",
          sprintf("Proteins remaining: %d\n", protein_count),
          sprintf("Samples remaining: %d\n", run_count),
          sprintf("Samples removed: %d\n", samples_removed_count),
          sprintf("Min peptides per sample: %d\n", input$min_peptides_per_sample),
          "State saved as: 'sample_filtered'\n"
        )
        
        # Add removed sample names if any
        if (samples_removed_count > 0) {
          result_text <- paste0(
            result_text,
            "\nRemoved samples:\n",
            paste(samples_removed, collapse = ", ")
          )
        }
        
        output$sample_results <- shiny::renderText(result_text)
        
        # Update filtering visualization and capture plot
        plot_grid <- updateProteinFiltering(
          data = filtered_s4@peptide_data,
          step_name = "6_sample_filtered",
          omic_type = omic_type,
          experiment_label = experiment_label,
          return_grid = TRUE,
          overwrite = TRUE
        )
        sample_plot(plot_grid)
        
        logger::log_info("Sample quality filter applied successfully")
        shiny::removeNotification("sample_working")
        shiny::showNotification("Sample quality filter applied successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error applying sample quality filter:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error", duration = 15)
        shiny::removeNotification("sample_working")
      })
    })
    
    # Revert Sample Filter
    shiny::observeEvent(input$revert_sample, {
      tryCatch({
        # Revert to the previous state in the history
        history <- workflow_data$state_manager$getHistory()
        if (length(history) > 1) {
          prev_state_name <- history[length(history) - 1]
          reverted_s4 <- workflow_data$state_manager$revertToState(prev_state_name)
          output$sample_results <- shiny::renderText(paste("Reverted to previous state:", prev_state_name))
          logger::log_info(paste("Reverted sample filter to", prev_state_name))
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
    
    # Render sample filter plot
    output$sample_plot <- shiny::renderPlot({
      shiny::req(sample_plot())
      grid::grid.draw(sample_plot())
    })
  })
}

