#' @title qualityControlAppletModule
#'
#' @description A Shiny module for the Quality Control step of the proteomics
#' workflow. It contains sub-tabs for different QC operations.
#'
#' @name qualityControlAppletModule
NULL

#' @rdname qualityControlAppletModule
#' @export
#' @import shiny
#' @import shinydashboard
qualityControlAppletUI <- function(id) {
  ns <- NS(id)
  
  shiny::fluidRow(
    shiny::column(12,
      shiny::wellPanel(
        shiny::h3("Quality Control"),
        shiny::p("This section provides tools for assessing and filtering the data based on quality metrics."),
        
        shiny::tabsetPanel(
          id = ns("qc_tabs"),
          
          # == Raw Data QC Sub-Tab ========================================
          shiny::tabPanel(
            "Raw Data Overview",
            icon = shiny::icon("images"),
            shiny::br(),
            shiny::fluidRow(
              shiny::column(3,
                shiny::wellPanel(
                  shiny::h4("Raw Data Visualization"),
                  shiny::p("This generates a plot summarizing the number of proteins identified at various stages of the raw data processing."),
                  shiny::actionButton(ns("run_raw_data_qc"), "Generate Plot", class = "btn-primary", width = "100%"),
                  shiny::hr(),
                  shiny::h5("Session Reset"),
                  shiny::p("Reset the analysis back to the raw data state. This will clear all filtering steps."),
                  shiny::actionButton(ns("revert_to_raw"), "Revert to Raw Data", 
                    class = "btn-warning", width = "100%")
                )
              ),
              shiny::column(9,
                                    shinyjqui::jqui_resizable(
                      shiny::plotOutput(ns("raw_data_plot"), height = "800px", width = "100%")
                    )
              )
            )
          ),
          
          # == Peptide Filtering Sub-Tab =================================
          shiny::tabPanel(
            "Peptide Filtering",
            icon = shiny::icon("filter"),
            shiny::br(),
            
            # Nested tabs for each filtering step
            shiny::tabsetPanel(
              id = ns("peptide_filter_tabs"),
              
              # Step 1: Q-Value Filter (chunk 10)
              shiny::tabPanel(
                "Q-Value Filter",
                shiny::br(),
                shiny::fluidRow(
                  shiny::column(4,
                    shiny::wellPanel(
                      shiny::h4("Q-Value & Proteotypic Peptide Filter"),
                      shiny::p("Filter peptides based on statistical confidence (q-value) and ensure only proteotypic peptides are retained."),
                      shiny::hr(),
                      
                      shiny::numericInput(ns("qvalue_threshold"), 
                        "Q-Value Threshold", 
                        value = 0.01, min = 0.001, max = 0.1, step = 0.001,
                        width = "100%"
                      ),
                      shiny::helpText("Lower = stricter peptide ID confidence (default: 0.01)"),
                      
                      shiny::numericInput(ns("global_qvalue_threshold"), 
                        "Global Q-Value Threshold", 
                        value = 0.01, min = 0.001, max = 0.1, step = 0.001,
                        width = "100%"
                      ),
                      shiny::helpText("Lower = stricter protein group ID confidence (default: 0.01)"),
                      
                      shiny::checkboxInput(ns("proteotypic_only"), 
                        "Keep only proteotypic peptides", 
                        value = TRUE
                      ),
                      shiny::helpText("TRUE = unique peptides only, FALSE = allow shared peptides (default: TRUE)"),
                      
                      shiny::hr(),
                      shiny::div(
                        shiny::actionButton(ns("apply_qvalue_filter"), "Apply Filter", 
                          class = "btn-primary", width = "48%"),
                        shiny::actionButton(ns("revert_qvalue"), "Revert", 
                          class = "btn-warning", width = "48%", style = "margin-left: 4%")
                      )
                    )
                  ),
                  shiny::column(8,
                    shiny::verbatimTextOutput(ns("qvalue_results")),
                    shiny::br(),
                    shinyjqui::jqui_resizable(
                      shiny::plotOutput(ns("qvalue_plot"), height = "800px", width = "100%")
                    )
                  )
                )
              ),
              
              # Step 2: Precursor Rollup (chunk 11)
              shiny::tabPanel(
                "Precursor Rollup",
                shiny::br(),
                shiny::fluidRow(
                  shiny::column(4,
                    shiny::wellPanel(
                      shiny::h4("Precursor to Peptide Rollup"),
                      shiny::p("Aggregate intensity measurements from multiple precursor ions to peptide level."),
                      shiny::hr(),
                      shiny::p("This step has no configurable parameters - it uses statistical methods to combine precursor measurements."),
                      shiny::hr(),
                      shiny::div(
                        shiny::actionButton(ns("apply_rollup"), "Apply Rollup", 
                          class = "btn-primary", width = "48%"),
                        shiny::actionButton(ns("revert_rollup"), "Revert", 
                          class = "btn-warning", width = "48%", style = "margin-left: 4%")
                      )
                    )
                  ),
                  shiny::column(8,
                    shiny::verbatimTextOutput(ns("rollup_results")),
                    shiny::br(),
                    shinyjqui::jqui_resizable(
                      shiny::plotOutput(ns("rollup_plot"), height = "800px", width = "100%")
                    )
                  )
                )
              ),
              
              # Step 3: Intensity Filter (chunk 12)
              shiny::tabPanel(
                "Intensity Filter",
                shiny::br(),
                shiny::fluidRow(
                  shiny::column(4,
                    shiny::wellPanel(
                      shiny::h4("Peptide Intensity Filter"),
                      shiny::p("Remove peptides with low intensity values or detected in too few samples."),
                      shiny::hr(),
                      
                      shiny::numericInput(ns("intensity_cutoff_percentile"), 
                        "Intensity Cutoff Percentile (%)", 
                        value = 1, min = 0.1, max = 10, step = 0.1,
                        width = "100%"
                      ),
                      shiny::helpText("Higher = stricter intensity threshold (default: 1%)"),
                      
                      shiny::numericInput(ns("proportion_samples_below_cutoff"), 
                        "Max Proportion Below Cutoff", 
                        value = 0.8, min = 0.1, max = 1.0, step = 0.1,
                        width = "100%"
                      ),
                      shiny::helpText("Higher = more relaxed filtering (default: 0.8)"),
                      
                      shiny::hr(),
                      shiny::div(
                        shiny::actionButton(ns("apply_intensity_filter"), "Apply Filter", 
                          class = "btn-primary", width = "48%"),
                        shiny::actionButton(ns("revert_intensity"), "Revert", 
                          class = "btn-warning", width = "48%", style = "margin-left: 4%")
                      )
                    )
                  ),
                  shiny::column(8,
                    shiny::verbatimTextOutput(ns("intensity_results")),
                    shiny::br(),
                    shinyjqui::jqui_resizable(
                      shiny::plotOutput(ns("intensity_plot"), height = "800px", width = "100%")
                    )
                  )
                )
              ),
              
              # Step 4: Protein Peptide Count Filter (chunk 13)
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
              ),
              
              # Step 5: Sample Quality Filter (chunk 14)
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
              ),
              
              # Step 6: Replicate Filter (chunk 15)
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
              ),
              
              # Step 7: Missing Value Imputation (chunk 16)
              shiny::tabPanel(
                "Imputation",
                shiny::br(),
                shiny::fluidRow(
                  shiny::column(4,
                    shiny::wellPanel(
                      shiny::h4("Missing Value Imputation"),
                      shiny::p("Impute missing values using technical replicate averages."),
                      shiny::hr(),
                      
                      shiny::numericInput(ns("proportion_missing_values"), 
                        "Max Proportion Missing", 
                        value = 0.8, min = 0.1, max = 1.0, step = 0.1,
                        width = "100%"
                      ),
                      shiny::helpText("Lower = imputes more readily, higher = more stringent (default: 0.8)"),
                      
                      shiny::hr(),
                      shiny::div(
                        shiny::actionButton(ns("apply_imputation"), "Apply Imputation", 
                          class = "btn-primary", width = "48%"),
                        shiny::actionButton(ns("revert_imputation"), "Revert", 
                          class = "btn-warning", width = "48%", style = "margin-left: 4%")
                      )
                    )
                  ),
                  shiny::column(8,
                    shiny::verbatimTextOutput(ns("imputation_results")),
                    shiny::br(),
                    shinyjqui::jqui_resizable(
                      shiny::plotOutput(ns("imputation_plot"), height = "800px", width = "100%")
                    )
                  )
                )
              )
            )
          ),
          
          # == Protein Filtering Sub-Tab =================================
          shiny::tabPanel(
            "Protein Filtering",
            icon = shiny::icon("filter"),
            shiny::br(),
            
            # Nested tabs for each protein processing step
            shiny::tabsetPanel(
              id = ns("protein_filter_tabs"),
              
              # Step 1: IQ Protein Rollup (chunk 17)
              shiny::tabPanel(
                "IQ Protein Rollup",
                shiny::br(),
                shiny::fluidRow(
                  shiny::column(4,
                    shiny::wellPanel(
                      shiny::h4("Peptide-to-Protein Rollup"),
                      shiny::p("Aggregate peptide-level data to protein-level quantification using the IQ algorithm (MaxLFQ)."),
                      shiny::hr(),
                      shiny::p("This step uses the IQ tool to implement the MaxLFQ algorithm for protein quantification, then automatically creates a ProteinQuantitativeData S4 object."),
                      shiny::hr(),
                      shiny::div(
                        shiny::actionButton(ns("apply_iq_rollup"), "Run IQ Rollup & Create S4 Object", 
                          class = "btn-primary", width = "48%"),
                        shiny::actionButton(ns("revert_iq_rollup"), "Revert", 
                          class = "btn-warning", width = "48%", style = "margin-left: 4%")
                      )
                    )
                  ),
                  shiny::column(8,
                    shiny::verbatimTextOutput(ns("iq_rollup_results")),
                    shiny::br(),
                    shinyjqui::jqui_resizable(
                      shiny::plotOutput(ns("iq_rollup_plot"), height = "800px", width = "100%")
                    )
                  )
                )
              ),
              
              # Step 2: Protein Accession Cleanup (chunk 19)
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
              ),
              
              # Step 4: Protein Intensity Filter (chunks 20+21 combined)
              shiny::tabPanel(
                "Protein Intensity Filter",
                shiny::br(),
                shiny::fluidRow(
                  shiny::column(4,
                    shiny::wellPanel(
                      shiny::h4("Missing Value Parameters & Intensity Filter"),
                      shiny::p("Configure missing value thresholds and filter proteins on intensity and missing values."),
                      shiny::hr(),
                      
                      # Missing value parameters (chunk 20)
                      shiny::h5("Missing Value Parameters"),
                      shiny::numericInput(ns("min_reps_per_group"), 
                        "Min Replicates per Group", 
                        value = 2, min = 1, max = 10, step = 1,
                        width = "100%"
                      ),
                      shiny::helpText("Minimum valid measurements required within each group (default: 2)"),
                      
                      shiny::numericInput(ns("min_groups"), 
                        "Min Groups Required", 
                        value = 2, min = 1, max = 10, step = 1,
                        width = "100%"
                      ),
                      shiny::helpText("Minimum groups that must meet the replicate threshold (default: 2)"),
                      
                      shiny::hr(),
                      
                      # Calculated percentages (read-only display)
                      shiny::h5("Calculated Thresholds"),
                      shiny::p(style = "background-color: #f0f0f0; padding: 10px; border-radius: 5px;",
                        shiny::strong("These values are automatically calculated from your replicate settings above:"),
                        shiny::br(),
                        shiny::textOutput(ns("calculated_groupwise_percent"), inline = TRUE),
                        shiny::br(),
                        shiny::textOutput(ns("calculated_max_groups_percent"), inline = TRUE)
                      ),
                      
                      shiny::numericInput(ns("proteins_intensity_cutoff_percentile"), 
                        "Protein Intensity Cutoff Percentile (%)", 
                        value = 1, min = 0.1, max = 10, step = 0.1,
                        width = "100%"
                      ),
                      shiny::helpText("Intensity threshold percentile (default: 1%)"),
                      
                      shiny::hr(),
                      shiny::div(
                        shiny::actionButton(ns("apply_protein_intensity_filter"), "Apply Filter", 
                          class = "btn-primary", width = "48%"),
                        shiny::actionButton(ns("revert_protein_intensity_filter"), "Revert", 
                          class = "btn-warning", width = "48%", style = "margin-left: 4%")
                      )
                    )
                  ),
                  shiny::column(8,
                    shiny::verbatimTextOutput(ns("protein_intensity_filter_results")),
                    shiny::br(),
                    shinyjqui::jqui_resizable(
                      shiny::plotOutput(ns("protein_intensity_filter_plot"), height = "800px", width = "100%")
                    )
                  )
                )
              ),
              
              # Step 5: Duplicate Protein Removal (chunk 22)
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
              ),
              
              # Step 6: Protein Replicate Filter (chunk 23)
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
          )
        )
      )
    )
  )
}

#' @rdname qualityControlAppletModule
#' @export
#' @import shiny
#' @importFrom logger log_info log_error
#' @importFrom grid grid.draw
qualityControlAppletServer <- function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
  shiny::moduleServer(id, function(input, output, session) {
    
    # Create global project_dirs object that updateProteinFiltering expects
    # This maps our experiment_paths to the expected global structure
    if (!exists("project_dirs", envir = .GlobalEnv)) {
      project_dirs <- list()
      project_dirs[[paste0(omic_type, "_", experiment_label)]] <- experiment_paths
      assign("project_dirs", project_dirs, envir = .GlobalEnv)
      logger::log_info("Created global project_dirs object for updateProteinFiltering compatibility")
    }
    
    # Reactive values to store plots from each filtering step
    raw_qc_plot <- reactiveVal(NULL)
    qvalue_plot <- reactiveVal(NULL)
    rollup_plot <- reactiveVal(NULL)
    intensity_plot <- reactiveVal(NULL)
    protein_peptide_plot <- reactiveVal(NULL)
    sample_plot <- reactiveVal(NULL)
    replicate_plot <- reactiveVal(NULL)
    imputation_plot <- reactiveVal(NULL)
    
    # Protein filtering plots
    iq_rollup_plot <- reactiveVal(NULL)
    accession_cleanup_plot <- reactiveVal(NULL)
    protein_intensity_filter_plot <- reactiveVal(NULL)
    duplicate_removal_plot <- reactiveVal(NULL)
    protein_replicate_filter_plot <- reactiveVal(NULL)
    
    # == Raw Data QC Logic ==================================================
    
    observeEvent(input$run_raw_data_qc, {
      shiny::req(workflow_data$state_manager)
      
      # Show a notification while running
      shiny::showNotification("Generating raw data QC plot...", id = "raw_qc_plot_working", duration = NULL)
      
      tryCatch({
        # 1. Retrieve the initial S4 object from the state manager
        logger::log_info("QC Step: Retrieving 'raw_data_s4' object from state manager.")
        s4_object <- workflow_data$state_manager$getState("raw_data_s4")
        
        shiny::req(s4_object)
        
        # 2. Run the updateProteinFiltering function
        # This function generates a plot and saves it, but we also want to capture it.
        logger::log_info("QC Step: Running updateProteinFiltering for '1_Raw Data'.")
        
        # The function needs 'data_cln' from the global env, which is not ideal.
        # For now, we get it from workflow_data.
        # A better long-term solution would be to refactor updateProteinFiltering to accept the S4 object.
        plot_grid <- updateProteinFiltering(
          data = s4_object, # Pass the S4 object directly
          step_name = "1_Raw Data",
          omic_type = omic_type,
          experiment_label = experiment_label,
          return_grid = TRUE,
          overwrite = TRUE
        )
        
        raw_qc_plot(plot_grid)
        
        logger::log_info("QC Step: Raw data QC plot generated successfully.")
        shiny::removeNotification("raw_qc_plot_working")
        shiny::showNotification("Plot generated.", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error generating raw data QC plot:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error", duration = 15)
        shiny::removeNotification("raw_qc_plot_working")
      })
    })
    
    # Render the plot
    output$raw_data_plot <- shiny::renderPlot({
      shiny::req(raw_qc_plot())
      # The return object is a grob, so it needs to be drawn
      grid::grid.draw(raw_qc_plot())
    })
    
    # Revert to Raw Data functionality
    observeEvent(input$revert_to_raw, {
      tryCatch({
        # Revert to raw data state and clear all subsequent states
        reverted_s4 <- workflow_data$state_manager$revertToState("raw_data_s4")
        
        # Clear all result outputs
        output$qvalue_results <- shiny::renderText("Reset to raw data state")
        output$rollup_results <- shiny::renderText("Reset to raw data state")
        output$intensity_results <- shiny::renderText("Reset to raw data state")
        output$protein_peptide_results <- shiny::renderText("Reset to raw data state")
        output$sample_results <- shiny::renderText("Reset to raw data state")
        output$replicate_results <- shiny::renderText("Reset to raw data state")
        output$imputation_results <- shiny::renderText("Reset to raw data state")
        
        logger::log_info("Session reset to raw data state - all filtering steps cleared")
        shiny::showNotification("Session reset to raw data state. All filtering steps have been cleared.", 
          type = "message", duration = 5)
        
      }, error = function(e) {
        msg <- paste("Error resetting session:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error")
      })
    })
    
    # == Peptide Filtering Logic ===============================================
    
    # Step 1: Q-Value Filter (chunk 10)
    observeEvent(input$apply_qvalue_filter, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying Q-value filter...", id = "qvalue_working", duration = NULL)
      
      tryCatch({
        # Get current S4 object from state manager
        current_s4 <- workflow_data$state_manager$getState("raw_data_s4")
        shiny::req(current_s4)
        
        logger::log_info("QC Step: Applying Q-value filter with thresholds {input$qvalue_threshold}, {input$global_qvalue_threshold}")
        
        # ✅ FIXED: Use updateConfigParameter to sync S4 object AND global config_list
        current_s4 <- updateConfigParameter(
          theObject = current_s4,
          function_name = "srlQvalueProteotypicPeptideClean",
          parameter_name = "qvalue_threshold",
          new_value = input$qvalue_threshold
        )
        
        current_s4 <- updateConfigParameter(
          theObject = current_s4,
          function_name = "srlQvalueProteotypicPeptideClean",
          parameter_name = "global_qvalue_threshold",
          new_value = input$global_qvalue_threshold
        )
        
        current_s4 <- updateConfigParameter(
          theObject = current_s4,
          function_name = "srlQvalueProteotypicPeptideClean",
          parameter_name = "choose_only_proteotypic_peptide",
          new_value = as.numeric(input$proteotypic_only)
        )
        
        # Apply S4 transformation (EXISTING S4 CODE - UNCHANGED)
        filtered_s4 <- srlQvalueProteotypicPeptideClean(theObject = current_s4)
        
        # Save new state in R6 manager
        workflow_data$state_manager$saveState(
          state_name = "qvalue_filtered",
          s4_data_object = filtered_s4,
          config_object = list(
            qvalue_threshold = input$qvalue_threshold,
            global_qvalue_threshold = input$global_qvalue_threshold,
            proteotypic_only = input$proteotypic_only
          ),
          description = "Applied Q-value and proteotypic peptide filter"
        )
        
        # Generate filtering summary
        protein_count <- filtered_s4@peptide_data |>
          dplyr::distinct(Protein.Ids) |>
          nrow()
        
        result_text <- paste(
          "Q-Value Filter Applied Successfully\n",
          "================================\n",
          sprintf("Proteins remaining: %d\n", protein_count),
          sprintf("Q-value threshold: %g\n", input$qvalue_threshold),
          sprintf("Global Q-value threshold: %g\n", input$global_qvalue_threshold),
          sprintf("Proteotypic only: %s\n", input$proteotypic_only),
          sprintf("State saved as: 'qvalue_filtered'\n")
        )
        
        output$qvalue_results <- shiny::renderText(result_text)
        
        # Update filtering visualization and capture plot
        plot_grid <- updateProteinFiltering(
          data = filtered_s4@peptide_data,
          step_name = "2_qval_filtered",
          omic_type = omic_type,
          experiment_label = experiment_label,
          return_grid = TRUE,
          overwrite = TRUE
        )
        qvalue_plot(plot_grid)
        
        logger::log_info("Q-value filter applied successfully")
        shiny::removeNotification("qvalue_working")
        shiny::showNotification("Q-value filter applied successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error applying Q-value filter:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error", duration = 15)
        shiny::removeNotification("qvalue_working")
      })
    })
    
    # Revert Q-Value Filter
    observeEvent(input$revert_qvalue, {
      tryCatch({
        # Revert to raw data state
        reverted_s4 <- workflow_data$state_manager$revertToState("raw_data_s4")
        
        output$qvalue_results <- shiny::renderText("Reverted to raw data state")
        
        logger::log_info("Reverted to raw data state")
        shiny::showNotification("Reverted to raw data state", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error reverting:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error")
      })
    })
    
    # Step 2: Precursor Rollup (chunk 11)
    observeEvent(input$apply_rollup, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying precursor rollup...", id = "rollup_working", duration = NULL)
      
      tryCatch({
        # Get current S4 object (should be qvalue_filtered or raw_data_s4)
        current_state <- workflow_data$state_manager$current_state
        current_s4 <- workflow_data$state_manager$getState(current_state)
        shiny::req(current_s4)
        
        logger::log_info("QC Step: Applying precursor rollup")
        
        # Apply S4 transformation (EXISTING S4 CODE - UNCHANGED)
        rolled_up_s4 <- rollUpPrecursorToPeptide(current_s4)
        
        # Save new state in R6 manager
        workflow_data$state_manager$saveState(
          state_name = "precursor_rollup",
          s4_data_object = rolled_up_s4,
          config_object = list(),
          description = "Applied precursor to peptide rollup"
        )
        
        # Generate summary
        protein_count <- rolled_up_s4@peptide_data |>
          dplyr::distinct(Protein.Ids) |>
          nrow()
        
        result_text <- paste(
          "Precursor Rollup Applied Successfully\n",
          "====================================\n",
          sprintf("Proteins remaining: %d\n", protein_count),
          "State saved as: 'precursor_rollup'\n"
        )
        
        output$rollup_results <- shiny::renderText(result_text)
        
        # Update filtering visualization and capture plot
        plot_grid <- updateProteinFiltering(
          data = rolled_up_s4@peptide_data,
          step_name = "3_precursor_rollup",
          omic_type = omic_type,
          experiment_label = experiment_label,
          return_grid = TRUE,
          overwrite = TRUE
        )
        rollup_plot(plot_grid)
        
        logger::log_info("Precursor rollup applied successfully")
        shiny::removeNotification("rollup_working")
        shiny::showNotification("Precursor rollup applied successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error applying precursor rollup:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error", duration = 15)
        shiny::removeNotification("rollup_working")
      })
    })
    
    # Revert Precursor Rollup
    observeEvent(input$revert_rollup, {
      tryCatch({
        # Revert to previous state (qvalue_filtered or raw_data_s4)
        history <- workflow_data$state_manager$getHistory()
        if ("qvalue_filtered" %in% history) {
          reverted_s4 <- workflow_data$state_manager$revertToState("qvalue_filtered")
          output$rollup_results <- shiny::renderText("Reverted to Q-value filtered state")
        } else {
          reverted_s4 <- workflow_data$state_manager$revertToState("raw_data_s4")
          output$rollup_results <- shiny::renderText("Reverted to raw data state")
        }
        
        logger::log_info("Reverted precursor rollup")
        shiny::showNotification("Reverted successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error reverting:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error")
      })
    })
    
    # Step 3: Intensity Filter (chunk 12)
    observeEvent(input$apply_intensity_filter, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying intensity filter...", id = "intensity_working", duration = NULL)
      
      tryCatch({
        # Get current S4 object
        current_state <- workflow_data$state_manager$current_state
        current_s4 <- workflow_data$state_manager$getState(current_state)
        shiny::req(current_s4)
        
        logger::log_info("QC Step: Applying intensity filter with cutoff {input$intensity_cutoff_percentile}%")
        
        # ✅ FIXED: Use updateConfigParameter to sync S4 object AND global config_list
        current_s4 <- updateConfigParameter(
          theObject = current_s4,
          function_name = "peptideIntensityFiltering",
          parameter_name = "peptides_intensity_cutoff_percentile",
          new_value = input$intensity_cutoff_percentile
        )
        
        current_s4 <- updateConfigParameter(
          theObject = current_s4,
          function_name = "peptideIntensityFiltering",
          parameter_name = "peptides_proportion_of_samples_below_cutoff",
          new_value = input$proportion_samples_below_cutoff
        )
        
        # Apply S4 transformation (EXISTING S4 CODE - UNCHANGED)
        filtered_s4 <- peptideIntensityFiltering(theObject = current_s4)
        
        # Save new state in R6 manager
        workflow_data$state_manager$saveState(
          state_name = "intensity_filtered",
          s4_data_object = filtered_s4,
          config_object = list(
            intensity_cutoff_percentile = input$intensity_cutoff_percentile,
            proportion_samples_below_cutoff = input$proportion_samples_below_cutoff
          ),
          description = "Applied peptide intensity filter"
        )
        
        # Generate summary
        protein_count <- filtered_s4@peptide_data |>
          dplyr::distinct(Protein.Ids) |>
          nrow()
        
        result_text <- paste(
          "Intensity Filter Applied Successfully\n",
          "====================================\n",
          sprintf("Proteins remaining: %d\n", protein_count),
          sprintf("Intensity cutoff percentile: %.1f%%\n", input$intensity_cutoff_percentile),
          sprintf("Max proportion below cutoff: %.1f\n", input$proportion_samples_below_cutoff),
          "State saved as: 'intensity_filtered'\n"
        )
        
        output$intensity_results <- shiny::renderText(result_text)
        
        # Update filtering visualization and capture plot
        plot_grid <- updateProteinFiltering(
          data = filtered_s4@peptide_data,
          step_name = "4_intensity_filtered",
          omic_type = omic_type,
          experiment_label = experiment_label,
          return_grid = TRUE,
          overwrite = TRUE
        )
        intensity_plot(plot_grid)
        
        logger::log_info("Intensity filter applied successfully")
        shiny::removeNotification("intensity_working")
        shiny::showNotification("Intensity filter applied successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error applying intensity filter:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error", duration = 15)
        shiny::removeNotification("intensity_working")
      })
    })
    
    # Revert Intensity Filter
    observeEvent(input$revert_intensity, {
      tryCatch({
        # Revert to previous state
        history <- workflow_data$state_manager$getHistory()
        prev_states <- c("precursor_rollup", "qvalue_filtered", "raw_data_s4")
        prev_state <- intersect(rev(history), prev_states)[1]
        
        reverted_s4 <- workflow_data$state_manager$revertToState(prev_state)
        output$intensity_results <- shiny::renderText(paste("Reverted to", prev_state, "state"))
        
        logger::log_info("Reverted intensity filter to {prev_state}")
        shiny::showNotification("Reverted successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error reverting:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error")
      })
    })
    
    # Step 4: Protein Peptide Count Filter (chunk 13)
    observeEvent(input$apply_protein_peptide_filter, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying protein peptide count filter...", id = "protein_peptide_working", duration = NULL)
      
      tryCatch({
        # Get current S4 object
        current_state <- workflow_data$state_manager$current_state
        current_s4 <- workflow_data$state_manager$getState(current_state)
        shiny::req(current_s4)
        
        logger::log_info("QC Step: Applying protein peptide count filter (min: {input$min_peptides_per_protein})")
        
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
    observeEvent(input$revert_protein_peptide, {
      tryCatch({
        # Revert to previous state
        history <- workflow_data$state_manager$getHistory()
        prev_states <- c("intensity_filtered", "precursor_rollup", "qvalue_filtered", "raw_data_s4")
        prev_state <- intersect(rev(history), prev_states)[1]
        
        reverted_s4 <- workflow_data$state_manager$revertToState(prev_state)
        output$protein_peptide_results <- shiny::renderText(paste("Reverted to", prev_state, "state"))
        
        logger::log_info("Reverted protein peptide filter to {prev_state}")
        shiny::showNotification("Reverted successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error reverting:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error")
      })
    })
    
    # Step 5: Sample Quality Filter (chunk 14)
    observeEvent(input$apply_sample_filter, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying sample quality filter...", id = "sample_working", duration = NULL)
      
      tryCatch({
        # Get current S4 object
        current_state <- workflow_data$state_manager$current_state
        current_s4 <- workflow_data$state_manager$getState(current_state)
        shiny::req(current_s4)
        
        logger::log_info("QC Step: Applying sample quality filter (min: {input$min_peptides_per_sample})")
        
        # ✅ FIXED: Use updateConfigParameter to sync S4 object AND global config_list
        current_s4 <- updateConfigParameter(
          theObject = current_s4,
          function_name = "filterMinNumPeptidesPerSample",
          parameter_name = "peptides_per_sample_cutoff",
          new_value = input$min_peptides_per_sample
        )
        
        # Apply S4 transformation (EXISTING S4 CODE - UNCHANGED)
        filtered_s4 <- filterMinNumPeptidesPerSample(theObject = current_s4)
        
        # Save new state in R6 manager
        workflow_data$state_manager$saveState(
          state_name = "sample_filtered",
          s4_data_object = filtered_s4,
          config_object = list(
            min_peptides_per_sample = input$min_peptides_per_sample
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
          sprintf("Min peptides per sample: %d\n", input$min_peptides_per_sample),
          "State saved as: 'sample_filtered'\n"
        )
        
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
    observeEvent(input$revert_sample, {
      tryCatch({
        # Revert to previous state
        history <- workflow_data$state_manager$getHistory()
        prev_states <- c("protein_peptide_filtered", "intensity_filtered", "precursor_rollup", "qvalue_filtered", "raw_data_s4")
        prev_state <- intersect(rev(history), prev_states)[1]
        
        reverted_s4 <- workflow_data$state_manager$revertToState(prev_state)
        output$sample_results <- shiny::renderText(paste("Reverted to", prev_state, "state"))
        
        logger::log_info("Reverted sample filter to {prev_state}")
        shiny::showNotification("Reverted successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error reverting:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error")
      })
    })
    
    # Step 6: Replicate Filter (chunk 15)
    observeEvent(input$apply_replicate_filter, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying replicate filter...", id = "replicate_working", duration = NULL)
      
      tryCatch({
        # Get current S4 object
        current_state <- workflow_data$state_manager$current_state
        current_s4 <- workflow_data$state_manager$getState(current_state)
        shiny::req(current_s4)
        
        logger::log_info("QC Step: Applying replicate filter (column: {input$replicate_group_column})")
        
        # Apply S4 transformation (EXISTING S4 CODE - UNCHANGED)
        # Note: This function takes the replicate column as a parameter
        filtered_s4 <- removePeptidesWithOnlyOneReplicate(
          current_s4,
          replicate_group_column = input$replicate_group_column
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
    observeEvent(input$revert_replicate, {
      tryCatch({
        # Revert to previous state
        history <- workflow_data$state_manager$getHistory()
        prev_states <- c("sample_filtered", "protein_peptide_filtered", "intensity_filtered", "precursor_rollup", "qvalue_filtered", "raw_data_s4")
        prev_state <- intersect(rev(history), prev_states)[1]
        
        reverted_s4 <- workflow_data$state_manager$revertToState(prev_state)
        output$replicate_results <- shiny::renderText(paste("Reverted to", prev_state, "state"))
        
        logger::log_info("Reverted replicate filter to {prev_state}")
        shiny::showNotification("Reverted successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error reverting:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error")
      })
    })
    
    # Step 7: Missing Value Imputation (chunk 16)
    observeEvent(input$apply_imputation, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying missing value imputation...", id = "imputation_working", duration = NULL)
      
      tryCatch({
        # Get current S4 object
        current_state <- workflow_data$state_manager$current_state
        current_s4 <- workflow_data$state_manager$getState(current_state)
        shiny::req(current_s4)
        
        logger::log_info("QC Step: Applying missing value imputation (proportion: {input$proportion_missing_values})")
        
        # ✅ FIXED: Use updateConfigParameter to sync S4 object AND global config_list
        current_s4 <- updateConfigParameter(
          theObject = current_s4,
          function_name = "peptideMissingValueImputation",
          parameter_name = "proportion_missing_values",
          new_value = input$proportion_missing_values
        )
        
        # Apply S4 transformation (EXISTING S4 CODE - UNCHANGED)
        imputed_s4 <- peptideMissingValueImputation(theObject = current_s4)
        
        # Save new state in R6 manager
        workflow_data$state_manager$saveState(
          state_name = "imputed",
          s4_data_object = imputed_s4,
          config_object = list(
            proportion_missing_values = input$proportion_missing_values
          ),
          description = "Applied missing value imputation using technical replicates"
        )
        
        # Generate summary
        protein_count <- imputed_s4@peptide_data |>
          dplyr::distinct(Protein.Ids) |>
          nrow()
        
        result_text <- paste(
          "Missing Value Imputation Applied Successfully\n",
          "============================================\n",
          sprintf("Proteins remaining: %d\n", protein_count),
          sprintf("Max proportion missing: %.1f\n", input$proportion_missing_values),
          "State saved as: 'imputed'\n",
          "\nReady for peptide-to-protein rollup step."
        )
        
        output$imputation_results <- shiny::renderText(result_text)
        
        # Update filtering visualization and capture plot
        plot_grid <- updateProteinFiltering(
          data = imputed_s4@peptide_data,
          step_name = "8_imputed",
          omic_type = omic_type,
          experiment_label = experiment_label,
          return_grid = TRUE,
          overwrite = TRUE
        )
        imputation_plot(plot_grid)
        
        logger::log_info("Missing value imputation applied successfully")
        shiny::removeNotification("imputation_working")
        shiny::showNotification("Missing value imputation applied successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error applying missing value imputation:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error", duration = 15)
        shiny::removeNotification("imputation_working")
      })
    })
    
    # Revert Imputation
    observeEvent(input$revert_imputation, {
      tryCatch({
        # Revert to previous state
        history <- workflow_data$state_manager$getHistory()
        prev_states <- c("replicate_filtered", "sample_filtered", "protein_peptide_filtered", "intensity_filtered", "precursor_rollup", "qvalue_filtered", "raw_data_s4")
        prev_state <- intersect(rev(history), prev_states)[1]
        
        reverted_s4 <- workflow_data$state_manager$revertToState(prev_state)
        output$imputation_results <- shiny::renderText(paste("Reverted to", prev_state, "state"))
        
        logger::log_info("Reverted imputation to {prev_state}")
        shiny::showNotification("Reverted successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error reverting:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error")
      })
    })
    
    # == Real-time Calculated Values Display ================================
    
    # Show calculated percentages based on replicate settings
    output$calculated_groupwise_percent <- shiny::renderText({
      shiny::req(input$min_reps_per_group, input$min_groups)
      
      # Get design matrix if available
      current_state <- workflow_data$state_manager$current_state
      if (!is.null(current_state)) {
        current_s4 <- workflow_data$state_manager$getState(current_state)
        if (!is.null(current_s4) && !is.null(current_s4@design_matrix)) {
          design_matrix <- current_s4@design_matrix
          
          # Calculate the same way as updateMissingValueParameters
          reps_per_group_tbl <- design_matrix |>
            dplyr::group_by(group) |>
            dplyr::summarise(n_reps = n()) |>
            dplyr::ungroup()
          
          group_thresholds <- reps_per_group_tbl |>
            dplyr::mutate(
              adjusted_min_reps = pmin(n_reps, input$min_reps_per_group),
              max_missing = n_reps - adjusted_min_reps,
              missing_percent = round((max_missing / n_reps) * 100, 3)
            )
          
          groupwise_cutoff <- max(group_thresholds$missing_percent)
          return(sprintf("Groupwise Percentage Cutoff: %.3f%% (calculated)", groupwise_cutoff))
        }
      }
      
      return("Groupwise Percentage Cutoff: [Design matrix needed]")
    })
    
    output$calculated_max_groups_percent <- shiny::renderText({
      shiny::req(input$min_reps_per_group, input$min_groups)
      
      # Get design matrix if available
      current_state <- workflow_data$state_manager$current_state
      if (!is.null(current_state)) {
        current_s4 <- workflow_data$state_manager$getState(current_state)
        if (!is.null(current_s4) && !is.null(current_s4@design_matrix)) {
          design_matrix <- current_s4@design_matrix
          
          # Calculate the same way as updateMissingValueParameters
          total_groups <- design_matrix |> dplyr::distinct(group) |> nrow()
          max_failing_groups <- total_groups - input$min_groups
          max_groups_cutoff <- round((max_failing_groups / total_groups) * 100, 3)
          
          return(sprintf("Max Groups Percentage Cutoff: %.3f%% (calculated)", max_groups_cutoff))
        }
      }
      
      return("Max Groups Percentage Cutoff: [Design matrix needed]")
    })

    # == Plot Rendering Functions ==============================================
    
    # Render Q-value filter plot
    output$qvalue_plot <- shiny::renderPlot({
      shiny::req(qvalue_plot())
      grid::grid.draw(qvalue_plot())
    })
    
    # Render rollup plot
    output$rollup_plot <- shiny::renderPlot({
      shiny::req(rollup_plot())
      grid::grid.draw(rollup_plot())
    })
    
    # Render intensity filter plot
    output$intensity_plot <- shiny::renderPlot({
      shiny::req(intensity_plot())
      grid::grid.draw(intensity_plot())
    })
    
    # Render protein peptide filter plot
    output$protein_peptide_plot <- shiny::renderPlot({
      shiny::req(protein_peptide_plot())
      grid::grid.draw(protein_peptide_plot())
    })
    
    # Render sample filter plot
    output$sample_plot <- shiny::renderPlot({
      shiny::req(sample_plot())
      grid::grid.draw(sample_plot())
    })
    
    # Render replicate filter plot
    output$replicate_plot <- shiny::renderPlot({
      shiny::req(replicate_plot())
      grid::grid.draw(replicate_plot())
    })
    
    # Render imputation plot
    output$imputation_plot <- shiny::renderPlot({
      shiny::req(imputation_plot())
      grid::grid.draw(imputation_plot())
    })
    
    # == Protein Filtering Logic ===========================================
    
    # Step 1: IQ Protein Rollup (chunk 17)
    observeEvent(input$apply_iq_rollup, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Running IQ protein rollup & creating S4 object...", id = "iq_rollup_working", duration = NULL)
      
      tryCatch({
        # Get the final peptide S4 object (should be 'imputed' state)
        current_state <- workflow_data$state_manager$current_state
        peptide_s4 <- workflow_data$state_manager$getState(current_state)
        shiny::req(peptide_s4)
        
        logger::log_info("Protein Processing: Starting IQ rollup from peptide state")
        
        # Save peptide data to TSV file for IQ input
        peptide_values_imputed_file <- file.path(
          experiment_paths$peptide_qc_dir,
          "peptide_values_imputed.tsv"
        )
        
        # Prepare data for IQ (ensure correct column names and format)
        peptide_data_for_iq <- peptide_s4@peptide_data |>
          dplyr::mutate(
            Q.Value = 0.0009,
            PG.Q.Value = 0.009,
            Peptide.Imputed = ifelse(is.na(Peptide.Imputed), 0, Peptide.Imputed)
          )
        
        vroom::vroom_write(peptide_data_for_iq, peptide_values_imputed_file)
        
        # Run IQ processing
        iq_output_file <- file.path(experiment_paths$protein_qc_dir, "iq_output_file.txt")
        
        iq::process_long_format(
          peptide_values_imputed_file,
          output_filename = iq_output_file,
          sample_id = "Run",
          primary_id = "Protein.Ids",
          secondary_id = "Stripped.Sequence",
          intensity_col = "Peptide.Imputed",
          filter_double_less = c("Q.Value" = "0.01", "PG.Q.Value" = "0.01"),
          normalization = "none"  # Critical - no normalization at this stage
        )
        
        # Wait for IQ output file to be available
        max_wait <- 30  # Maximum 30 seconds
        wait_count <- 0
        while (!file.exists(iq_output_file) && wait_count < max_wait) {
          Sys.sleep(1)
          wait_count <- wait_count + 1
        }
        
        if (!file.exists(iq_output_file)) {
          stop("IQ output file not created within timeout period")
        }
        
        # Read IQ output
        protein_log2_quant <- vroom::vroom(iq_output_file)
        
        logger::log_info("Protein Processing: Creating ProteinQuantitativeData S4 object")
        
        # Create ProteinQuantitativeData S4 object
        protein_obj <- ProteinQuantitativeData(
          protein_quant_table = protein_log2_quant,
          protein_id_column = "Protein.Ids",
          protein_id_table = protein_log2_quant |> dplyr::distinct(Protein.Ids),
          design_matrix = peptide_s4@design_matrix,
          sample_id = "Run",
          group_id = "group",
          technical_replicate_id = "replicates",
          args = peptide_s4@args
        )
        
        # Save S4 object as new state (combined rollup + S4 creation)
        workflow_data$state_manager$saveState(
          state_name = "protein_s4_created",
          s4_data_object = protein_obj,
          config_object = list(
            iq_output_file = iq_output_file,
            peptide_input_file = peptide_values_imputed_file,
            s4_class = "ProteinQuantitativeData",
            protein_id_column = "Protein.Ids"
          ),
          description = "IQ protein rollup completed and ProteinQuantitativeData S4 object created"
        )
        
        # Generate summary
        protein_count <- protein_obj@protein_quant_table |>
          dplyr::distinct(Protein.Ids) |>
          nrow()
        
        result_text <- paste(
          "IQ Protein Rollup & S4 Object Creation Completed Successfully\n",
          "============================================================\n",
          sprintf("Proteins quantified: %d\n", protein_count),
          sprintf("Samples: %d\n", ncol(protein_obj@protein_quant_table) - 1),
          sprintf("Algorithm: MaxLFQ (via IQ tool)\n"),
          sprintf("S4 Class: %s\n", class(protein_obj)[1]),
          sprintf("Design matrix: %s\n", paste(colnames(protein_obj@design_matrix), collapse = ", ")),
          sprintf("Output file: %s\n", basename(iq_output_file)),
          "State saved as: 'protein_s4_created'\n",
          "\nReady for protein accession cleanup."
        )
        
        output$iq_rollup_results <- shiny::renderText(result_text)
        
        # Update filtering visualization and capture plot
        plot_grid <- updateProteinFiltering(
          data = protein_obj@protein_quant_table,
          step_name = "9_protein_s4_created",
          omic_type = omic_type,
          experiment_label = experiment_label,
          return_grid = TRUE,
          overwrite = TRUE
        )
        iq_rollup_plot(plot_grid)
        
        logger::log_info("IQ protein rollup and S4 object creation completed successfully")
        shiny::removeNotification("iq_rollup_working")
        shiny::showNotification("IQ protein rollup & S4 object creation completed successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error in IQ protein rollup & S4 creation:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error", duration = 15)
        shiny::removeNotification("iq_rollup_working")
      })
    })
    
    # Revert IQ Rollup
    observeEvent(input$revert_iq_rollup, {
      tryCatch({
        # Revert to final peptide state (should be 'imputed')
        history <- workflow_data$state_manager$getHistory()
        peptide_states <- c("imputed", "replicate_filtered", "sample_filtered", "protein_peptide_filtered", 
                           "intensity_filtered", "precursor_rollup", "qvalue_filtered", "raw_data_s4")
        prev_state <- intersect(rev(history), peptide_states)[1]
        
        reverted_s4 <- workflow_data$state_manager$revertToState(prev_state)
        output$iq_rollup_results <- shiny::renderText(paste("Reverted to", prev_state, "state"))
        
        logger::log_info("Reverted IQ rollup to {prev_state}")
        shiny::showNotification("Reverted successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error reverting:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error")
      })
    })
    
    # Step 2: Protein Accession Cleanup (chunk 19)
    observeEvent(input$apply_accession_cleanup, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying protein accession cleanup...", id = "accession_cleanup_working", duration = NULL)
      
      tryCatch({
        # Get current ProteinQuantitativeData S4 object
        current_s4 <- workflow_data$state_manager$getState("protein_s4_created")
        shiny::req(current_s4)
        
        logger::log_info("Protein Processing: Applying accession cleanup with delimiter: {input$delimiter}")
        
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
    observeEvent(input$revert_accession_cleanup, {
      tryCatch({
        # Revert to protein S4 created state
        reverted_s4 <- workflow_data$state_manager$revertToState("protein_s4_created")
        output$accession_cleanup_results <- shiny::renderText("Reverted to protein S4 created state")
        
        logger::log_info("Reverted accession cleanup to protein_s4_created")
        shiny::showNotification("Reverted successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error reverting:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error")
      })
    })
    
    # Step 4: Protein Intensity Filter (chunks 20+21 combined)
    observeEvent(input$apply_protein_intensity_filter, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying protein intensity filter...", id = "protein_intensity_filter_working", duration = NULL)
      
      tryCatch({
        # Get current ProteinQuantitativeData S4 object
        current_state <- workflow_data$state_manager$current_state
        current_s4 <- workflow_data$state_manager$getState(current_state)
        shiny::req(current_s4)
        
        logger::log_info("Protein Processing: Applying missing value parameters and intensity filter")
        
        # First: Update missing value parameters (chunk 20)
        config_list <- current_s4@args
        config_list <- updateMissingValueParameters(
          current_s4@design_matrix,
          config_list,
          min_reps_per_group = input$min_reps_per_group,
          min_groups = input$min_groups
        )
        current_s4@args <- config_list
        
        # ✅ FIXED: Use ONLY the calculated percentages from updateMissingValueParameters
        # The replicate-based inputs (min_reps_per_group, min_groups) automatically calculate the percentages
        # No need to override with UI percentage inputs - they would conflict!
        
        # Only update the intensity cutoff percentile since it's not calculated by updateMissingValueParameters
        current_s4 <- updateConfigParameter(
          theObject = current_s4,
          function_name = "removeRowsWithMissingValuesPercent",
          parameter_name = "proteins_intensity_cutoff_percentile",
          new_value = input$proteins_intensity_cutoff_percentile
        )
        
        # Apply S4 transformation (EXISTING S4 CODE - UNCHANGED)
        filtered_s4 <- removeRowsWithMissingValuesPercent(current_s4)
        
        # Save new state
        workflow_data$state_manager$saveState(
          state_name = "protein_intensity_filtered",
          s4_data_object = filtered_s4,
          config_object = list(
            min_reps_per_group = input$min_reps_per_group,
            min_groups = input$min_groups,
            groupwise_percentage_cutoff = config_list$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff,
            max_groups_percentage_cutoff = config_list$removeRowsWithMissingValuesPercent$max_groups_percentage_cutoff,
            proteins_intensity_cutoff_percentile = input$proteins_intensity_cutoff_percentile
          ),
          description = "Applied missing value parameters and protein intensity filter"
        )
        
        # Generate summary
        protein_count <- filtered_s4@protein_quant_table |>
          dplyr::distinct(Protein.Ids) |>
          nrow()
        
        result_text <- paste(
          "Protein Intensity Filter Applied Successfully\n",
          "============================================\n",
          sprintf("Proteins remaining: %d\n", protein_count),
          sprintf("Min replicates per group: %d\n", input$min_reps_per_group),
          sprintf("Min groups required: %d\n", input$min_groups),
          sprintf("Groupwise %% cutoff: %.3f%% (calculated)\n", config_list$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff),
          sprintf("Max groups %% cutoff: %.3f%% (calculated)\n", config_list$removeRowsWithMissingValuesPercent$max_groups_percentage_cutoff),
          sprintf("Intensity cutoff percentile: %.1f%%\n", input$proteins_intensity_cutoff_percentile),
          "State saved as: 'protein_intensity_filtered'\n"
        )
        
        output$protein_intensity_filter_results <- shiny::renderText(result_text)
        
        # Update filtering visualization and capture plot
        plot_grid <- updateProteinFiltering(
          data = filtered_s4@protein_quant_table,
          step_name = "11_protein_intensity_filtered",
          omic_type = omic_type,
          experiment_label = experiment_label,
          return_grid = TRUE,
          overwrite = TRUE
        )
        protein_intensity_filter_plot(plot_grid)
        
        logger::log_info("Protein intensity filter applied successfully")
        shiny::removeNotification("protein_intensity_filter_working")
        shiny::showNotification("Protein intensity filter applied successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error applying protein intensity filter:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error", duration = 15)
        shiny::removeNotification("protein_intensity_filter_working")
      })
    })
    
    # Revert Protein Intensity Filter
    observeEvent(input$revert_protein_intensity_filter, {
      tryCatch({
        # Revert to previous state
        history <- workflow_data$state_manager$getHistory()
        prev_states <- c("protein_accession_cleaned", "protein_s4_created")
        prev_state <- intersect(rev(history), prev_states)[1]
        
        reverted_s4 <- workflow_data$state_manager$revertToState(prev_state)
        output$protein_intensity_filter_results <- shiny::renderText(paste("Reverted to", prev_state, "state"))
        
        logger::log_info("Reverted protein intensity filter to {prev_state}")
        shiny::showNotification("Reverted successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error reverting:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error")
      })
    })
    
    # Step 5: Duplicate Protein Removal (chunk 22)
    observeEvent(input$apply_duplicate_removal, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Removing duplicate proteins...", id = "duplicate_removal_working", duration = NULL)
      
      tryCatch({
        # Get current ProteinQuantitativeData S4 object
        current_state <- workflow_data$state_manager$current_state
        current_s4 <- workflow_data$state_manager$getState(current_state)
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
    observeEvent(input$revert_duplicate_removal, {
      tryCatch({
        # Revert to previous state
        history <- workflow_data$state_manager$getHistory()
        prev_states <- c("protein_intensity_filtered", "protein_accession_cleaned", "protein_s4_created")
        prev_state <- intersect(rev(history), prev_states)[1]
        
        reverted_s4 <- workflow_data$state_manager$revertToState(prev_state)
        output$duplicate_removal_results <- shiny::renderText(paste("Reverted to", prev_state, "state"))
        
        logger::log_info("Reverted duplicate removal to {prev_state}")
        shiny::showNotification("Reverted successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error reverting:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error")
      })
    })
    
    # Step 6: Protein Replicate Filter (chunk 23)
    observeEvent(input$apply_protein_replicate_filter, {
      shiny::req(workflow_data$state_manager)
      
      shiny::showNotification("Applying protein replicate filter...", id = "protein_replicate_filter_working", duration = NULL)
      
      tryCatch({
        # Get current ProteinQuantitativeData S4 object
        current_state <- workflow_data$state_manager$current_state
        current_s4 <- workflow_data$state_manager$getState(current_state)
        shiny::req(current_s4)
        
        logger::log_info("Protein Processing: Applying protein replicate filter with {input$parallel_cores} cores")
        
        # Set up parallel processing
        core_utilisation <- new_cluster(input$parallel_cores)
        
        # Apply S4 transformation (EXISTING S4 CODE - UNCHANGED)
        filtered_s4 <- removeProteinsWithOnlyOneReplicate(
          current_s4,
          core_utilisation,
          grouping_variable = input$protein_grouping_variable
        )
        
        # Save filtered data to file (as in original workflow)
        output_file <- file.path(experiment_paths$protein_qc_dir, "remove_proteins_with_only_one_rep.tsv")
        vroom::vroom_write(filtered_s4@protein_quant_table, output_file)
        
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
    observeEvent(input$revert_protein_replicate_filter, {
      tryCatch({
        # Revert to previous state
        history <- workflow_data$state_manager$getHistory()
        prev_states <- c("duplicates_removed", "protein_intensity_filtered", "protein_accession_cleaned", "protein_s4_created")
        prev_state <- intersect(rev(history), prev_states)[1]
        
        reverted_s4 <- workflow_data$state_manager$revertToState(prev_state)
        output$protein_replicate_filter_results <- shiny::renderText(paste("Reverted to", prev_state, "state"))
        
        logger::log_info("Reverted protein replicate filter to {prev_state}")
        shiny::showNotification("Reverted successfully", type = "message")
        
      }, error = function(e) {
        msg <- paste("Error reverting:", e$message)
        logger::log_error(msg)
        shiny::showNotification(msg, type = "error")
      })
    })
    
    # == Protein Plot Rendering Functions ==================================
    
    # Render IQ rollup plot
    output$iq_rollup_plot <- shiny::renderPlot({
      shiny::req(iq_rollup_plot())
      grid::grid.draw(iq_rollup_plot())
    })
    
    # Render accession cleanup plot
    output$accession_cleanup_plot <- shiny::renderPlot({
      shiny::req(accession_cleanup_plot())
      grid::grid.draw(accession_cleanup_plot())
    })
    
    # Render protein intensity filter plot
    output$protein_intensity_filter_plot <- shiny::renderPlot({
      shiny::req(protein_intensity_filter_plot())
      grid::grid.draw(protein_intensity_filter_plot())
    })
    
    # Render duplicate removal plot
    output$duplicate_removal_plot <- shiny::renderPlot({
      shiny::req(duplicate_removal_plot())
      grid::grid.draw(duplicate_removal_plot())
    })
    
    # Render protein replicate filter plot
    output$protein_replicate_filter_plot <- shiny::renderPlot({
      shiny::req(protein_replicate_filter_plot())
      grid::grid.draw(protein_replicate_filter_plot())
    })
  })
} 