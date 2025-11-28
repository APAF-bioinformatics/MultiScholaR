#' @title Protein Accession Cleanup Module
#'
#' @description A Shiny module for performing protein accession cleanup.
#'
#' @name mod_prot_qc_protein_cleanup
NULL

#' @rdname mod_prot_qc_protein_cleanup
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column wellPanel h4 p hr textInput helpText div icon strong textOutput verbatimTextOutput selectInput actionButton plotOutput
mod_prot_qc_protein_cleanup_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tabPanel(
    "Accession Cleanup",
    shiny::br(),
    shiny::fluidRow(
      shiny::column(4,
        shiny::wellPanel(
          shiny::h4("Best Protein Accession Selection"),
          shiny::p("Choose the best protein accession for ambiguous mappings using FASTA metadata."),
          shiny::hr(),
          
          shiny::textInput(ns("delimiter"), 
            "Protein ID Delimiter", 
            value = ";",
            width = "100%"
          ),
          shiny::helpText("Character separating multiple protein IDs (default: ';')"),
          
          # Detected delimiters info box
          shiny::div(
            style = "background-color: #e7f3ff; padding: 10px; border-left: 4px solid #2196F3; border-radius: 4px; margin-top: 10px;",
            shiny::icon("info-circle", style = "color: #2196F3;"),
            shiny::strong(" Detected Delimiters in Your Data:"),
            shiny::br(),
            shiny::textOutput(ns("detected_delimiters"), inline = TRUE)
          ),
          
          # Sample protein IDs display
          shiny::div(
            style = "background-color: #f5f5f5; padding: 10px; border-left: 4px solid #9E9E9E; border-radius: 4px; margin-top: 10px; font-family: monospace; font-size: 11px;",
            shiny::icon("list", style = "color: #9E9E9E;"),
            shiny::strong(" Sample Protein IDs:"),
            shiny::br(),
            shiny::verbatimTextOutput(ns("sample_protein_ids"))
          ),
          
          shiny::selectInput(ns("aggregation_method"), 
            "Aggregation Method", 
            choices = c("sum", "mean", "median"),
            selected = "mean",
            width = "100%"
          ),
          shiny::helpText("Method for combining duplicate protein measurements (S4 default: mean)"),
          
          shiny::hr(),
          shiny::div(
            shiny::actionButton(ns("apply_accession_cleanup"), "Apply Cleanup", 
              class = "btn-primary", width = "48%"),
            shiny::actionButton(ns("revert_accession_cleanup"), "Revert", 
              class = "btn-warning", width = "48%", style = "margin-left: 4%")
          )
        )
      ),
      shiny::column(8,
        shiny::verbatimTextOutput(ns("accession_cleanup_results")),
        shiny::br(),
        shinyjqui::jqui_resizable(
          shiny::plotOutput(ns("accession_cleanup_plot"), height = "800px", width = "100%")
        )
      )
    )
  )
}

#' @rdname mod_prot_qc_protein_cleanup
#' @export
#' @importFrom shiny moduleServer reactiveVal observeEvent req showNotification removeNotification renderText renderPlot textOutput
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
mod_prot_qc_protein_cleanup_server <- function(id, workflow_data, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    accession_cleanup_plot <- shiny::reactiveVal(NULL)
    
    # Step 2: Protein Accession Cleanup (chunk 19)
    shiny::observeEvent(input$apply_accession_cleanup, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying protein accession cleanup...", id = "accession_cleanup_working", duration = NULL)
      
      tryCatch({
        # Get current ProteinQuantitativeData S4 object from the active state
        current_s4 <- workflow_data$state_manager$getState()
        shiny::req(current_s4)
        
        logger::log_info("Protein Processing: Applying accession cleanup with delimiter: {input$delimiter}")
        
        # Count proteins before cleanup
        proteins_before <- current_s4@protein_quant_table |>
          dplyr::distinct(Protein.Ids) |>
          nrow()
        
        # Check if aa_seq_tbl_final exists in workflow_data (from FASTA processing)
        if (exists("aa_seq_tbl_final", envir = .GlobalEnv)) {
          aa_seq_tbl_final <- get("aa_seq_tbl_final", envir = .GlobalEnv)
          aa_seq_tbl_final <- aa_seq_tbl_final |>
            dplyr::rename(uniprot_acc = database_id)
          
          # Apply S4 transformation (EXISTING S4 CODE - UNCHANGED)
          cleaned_s4 <- chooseBestProteinAccession(
            theObject = current_s4,
            delim = input$delimiter,
            seqinr_obj = aa_seq_tbl_final,
            seqinr_accession_column = "uniprot_acc",
            replace_zero_with_na = TRUE,
            aggregation_method = input$aggregation_method
          )
          
          cleanup_applied <- TRUE
        } else {
          # No FASTA data available, skip cleanup
          cleaned_s4 <- current_s4
          cleanup_applied <- FALSE
        }
        
        # Count proteins after cleanup
        proteins_after <- cleaned_s4@protein_quant_table |>
          dplyr::distinct(Protein.Ids) |>
          nrow()
        
        # Track cleanup results in workflow_data
        workflow_data$accession_cleanup_results <- list(
          cleanup_applied = cleanup_applied,
          delimiter_used = input$delimiter,
          aggregation_method = input$aggregation_method,
          proteins_before = proteins_before,
          proteins_after = proteins_after,
          had_full_metadata = if (!is.null(workflow_data$fasta_metadata)) {
            workflow_data$fasta_metadata$has_protein_evidence && 
            workflow_data$fasta_metadata$has_gene_names
          } else {
            FALSE
          },
          timestamp = Sys.time()
        )
        
        # Also store in qc_params for consistency
        if (is.null(workflow_data$qc_params)) {
          workflow_data$qc_params <- list()
        }
        if (is.null(workflow_data$qc_params$protein_qc)) {
          workflow_data$qc_params$protein_qc <- list()
        }
        workflow_data$qc_params$protein_qc$accession_cleanup <- workflow_data$accession_cleanup_results
        
        # Save accession cleanup results to file for persistence
        tryCatch({
          if (exists("experiment_paths") && !is.null(experiment_paths$source_dir)) {
            accession_cleanup_file <- file.path(experiment_paths$source_dir, "accession_cleanup_results.RDS")
            saveRDS(workflow_data$accession_cleanup_results, accession_cleanup_file)
            logger::log_info(sprintf("Saved accession cleanup results to: %s", accession_cleanup_file))
          }
        }, error = function(e) {
          logger::log_warn(sprintf("Could not save accession cleanup results file: %s", e$message))
        })
        
        logger::log_info(sprintf("Accession cleanup results tracked: %d -> %d proteins", 
                                proteins_before, proteins_after))
        
        # Save new state
        workflow_data$state_manager$saveState(
          state_name = "protein_accession_cleaned",
          s4_data_object = cleaned_s4,
          config_object = list(
            delimiter = input$delimiter,
            aggregation_method = input$aggregation_method,
            cleanup_applied = cleanup_applied
          ),
          description = if (cleanup_applied) "Applied protein accession cleanup" else "Skipped accession cleanup (no FASTA data)"
        )
        
        # Generate summary
        protein_count <- cleaned_s4@protein_quant_table |>
          dplyr::distinct(Protein.Ids) |>
          nrow()
        
        result_text <- if (cleanup_applied) {
          paste(
            "Protein Accession Cleanup Applied Successfully\n",
            "==============================================\n",
            sprintf("Proteins remaining: %d\n", protein_count),
            sprintf("Delimiter: %s\n", input$delimiter),
            sprintf("Aggregation method: %s\n", input$aggregation_method),
            "State saved as: 'protein_accession_cleaned'\n"
          )
        } else {
          paste(
            "Protein Accession Cleanup Skipped\n",
            "=================================\n",
            sprintf("Proteins remaining: %d\n", protein_count),
            "Reason: No FASTA data available for cleanup\n",
            "State saved as: 'protein_accession_cleaned'\n"
          )
        }
        
        output$accession_cleanup_results <- shiny::renderText(result_text)
        
        # Update filtering visualization and capture plot
        plot_grid <- updateProteinFiltering(
          data = cleaned_s4@protein_quant_table,
          step_name = "10_protein_accession_cleaned",
          omic_type = omic_type,
          experiment_label = experiment_label,
          return_grid = TRUE,
          overwrite = TRUE
        )
        accession_cleanup_plot(plot_grid)
        
        logger::log_info("Protein accession cleanup completed")
        shiny::removeNotification("accession_cleanup_working")
        shiny::showNotification("Protein accession cleanup completed", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error in protein accession cleanup:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error", duration = 15)
        shiny::removeNotification("accession_cleanup_working")
      })
    })
    
    # Revert Accession Cleanup
    shiny::observeEvent(input$revert_accession_cleanup, {
      tryCatch({
        # Revert to the previous state in the history
        history <- workflow_data$state_manager$getHistory()
        if (length(history) > 1) {
          prev_state_name <- history[length(history) - 1]
          reverted_s4 <- workflow_data$state_manager$revertToState(prev_state_name)
          output$accession_cleanup_results <- shiny::renderText(paste("Reverted to previous state:", prev_state_name))
          logger::log_info(paste("Reverted accession cleanup to", prev_state_name))
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
    
    # Render accession cleanup plot
    output$accession_cleanup_plot <- shiny::renderPlot({
      shiny::req(accession_cleanup_plot())
      grid::grid.draw(accession_cleanup_plot())
    })
  })
}

