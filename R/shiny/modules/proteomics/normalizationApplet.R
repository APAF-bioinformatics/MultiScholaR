#' @title normalizationAppletModule
#'
#' @description A Shiny module for the Normalization step of the proteomics
#' workflow. Handles normalization methods, RUV-III batch correction, and
#' correlation-based sample filtering before preparing data for DE analysis.
#'
#' @name normalizationAppletModule
NULL

#' @rdname normalizationAppletModule
#' @export
#' @import shiny
#' @import shinydashboard
normalizationAppletUI <- function(id) {
  ns <- NS(id)
  
  shiny::fluidRow(
    shiny::column(3,
      shiny::wellPanel(
        shiny::h4("Normalization Options"),
        
        # Normalization method selection
        shiny::selectInput(
          ns("norm_method"),
          "Normalization Method:",
          choices = list(
            "Cyclic Loess" = "cyclicloess",
            "Quantile" = "quantile", 
            "Scale (Median Absolute Values)" = "scale",
            "None" = "none"
          ),
          selected = "cyclicloess",
          width = "100%"
        ),
        shiny::helpText("Method for between-sample normalization (default: Cyclic Loess)"),
        
        shiny::hr(),
        
        # Plot aesthetics controls
        shiny::h4("Plot Aesthetics"),
        shiny::selectInput(
          ns("color_variable"),
          "Color by:",
          choices = c("group", "factor1", "factor2"),
          selected = "group",
          width = "100%"
        ),
        shiny::helpText("Variable to use for plot coloring"),
        
        shiny::selectInput(
          ns("shape_variable"),
          "Shape by:",
          choices = c("group", "factor1", "factor2"),
          selected = "group",
          width = "100%"
        ),
        shiny::helpText("Variable to use for plot shapes"),
        
        shiny::hr(),
        
        # RUV-III options
        shiny::h4("RUV-III Batch Correction"),
        shiny::radioButtons(
          ns("ruv_mode"),
          "RUV Parameter Tuning:",
          choices = list(
            "Automatic (recommended)" = "automatic",
            "Manual (advanced users)" = "manual"
          ),
          selected = "automatic",
          width = "100%"
        ),
        
        # Manual RUV controls (shown when manual mode selected)
        shiny::conditionalPanel(
          condition = "input.ruv_mode == 'manual'",
          ns = ns,
          
          shiny::sliderInput(
            ns("ruv_percentage"),
            "% Proteins as Negative Controls:",
            min = 1,
            max = 20,
            value = 5,
            step = 1,
            width = "100%"
          ),
          shiny::helpText("Percentage of most stable proteins used as negative controls"),
          
          shiny::numericInput(
            ns("ruv_k"),
            "Number of Factors (k):",
            value = NULL,
            min = 1,
            max = 10,
            width = "100%"
          ),
          shiny::helpText("Will be auto-populated from canonical correlation analysis"),
          
          shiny::br()
        ),
        
        # Normalization action button
        shiny::actionButton(
          ns("run_normalization"),
          "Run Normalization & RUV",
          class = "btn-primary",
          width = "100%"
        ),
        
        shiny::br(),
        shiny::br(),
        
        # Reset button
        shiny::actionButton(
          ns("reset_normalization"),
          "Reset to Pre-Normalization",
          class = "btn-warning",
          width = "100%"
        )
      )
    ),
    
    # QC plot panels with 3-column layout structure and correlation filtering tab
    shiny::column(9,
      shiny::tabsetPanel(
        id = ns("norm_qc_tabs"),
        
        # PCA Tab
        shiny::tabPanel(
          "PCA",
          icon = shiny::icon("project-diagram"),
          shiny::br(),
          shiny::fluidRow(
            shiny::column(4,
              shiny::h5("Post-Filtering", style = "text-align: center;"),
              shinyjqui::jqui_resizable(
                shiny::plotOutput(ns("pca_post_filtering"), height = "400px")
              )
            ),
            shiny::column(4,
              shiny::h5("Post-Normalization", style = "text-align: center;"),
              shinyjqui::jqui_resizable(
                shiny::plotOutput(ns("pca_post_normalization"), height = "400px")
              )
            ),
            shiny::column(4,
              shiny::h5("RUV-Corrected", style = "text-align: center;"),
              shinyjqui::jqui_resizable(
                shiny::plotOutput(ns("pca_ruv_corrected"), height = "400px")
              )
            )
          )
        ),
        
        # Density Tab
        shiny::tabPanel(
          "Density",
          icon = shiny::icon("chart-area"),
          shiny::br(),
          shiny::fluidRow(
            shiny::column(4,
              shiny::h5("Post-Filtering", style = "text-align: center;"),
              shinyjqui::jqui_resizable(
                shiny::plotOutput(ns("density_post_filtering"), height = "400px")
              )
            ),
            shiny::column(4,
              shiny::h5("Post-Normalization", style = "text-align: center;"),
              shinyjqui::jqui_resizable(
                shiny::plotOutput(ns("density_post_normalization"), height = "400px")
              )
            ),
            shiny::column(4,
              shiny::h5("RUV-Corrected", style = "text-align: center;"),
              shinyjqui::jqui_resizable(
                shiny::plotOutput(ns("density_ruv_corrected"), height = "400px")
              )
            )
          )
        ),
        
        # RLE Tab
        shiny::tabPanel(
          "RLE",
          icon = shiny::icon("chart-line"),
          shiny::br(),
          shiny::fluidRow(
            shiny::column(4,
              shiny::h5("Post-Filtering", style = "text-align: center;"),
              shinyjqui::jqui_resizable(
                shiny::plotOutput(ns("rle_post_filtering"), height = "400px")
              )
            ),
            shiny::column(4,
              shiny::h5("Post-Normalization", style = "text-align: center;"),
              shinyjqui::jqui_resizable(
                shiny::plotOutput(ns("rle_post_normalization"), height = "400px")
              )
            ),
            shiny::column(4,
              shiny::h5("RUV-Corrected", style = "text-align: center;"),
              shinyjqui::jqui_resizable(
                shiny::plotOutput(ns("rle_ruv_corrected"), height = "400px")
              )
            )
          )
        ),
        
        # Correlation Tab
        shiny::tabPanel(
          "Correlation",
          icon = shiny::icon("th"),
          shiny::br(),
          shiny::fluidRow(
            shiny::column(4,
              shiny::h5("Post-Filtering", style = "text-align: center;"),
              shinyjqui::jqui_resizable(
                shiny::plotOutput(ns("correlation_post_filtering"), height = "400px")
              )
            ),
            shiny::column(4,
              shiny::h5("Post-Normalization", style = "text-align: center;"),
              shinyjqui::jqui_resizable(
                shiny::plotOutput(ns("correlation_post_normalization"), height = "400px")
              )
            ),
            shiny::column(4,
              shiny::h5("RUV-Corrected", style = "text-align: center;"),
              shinyjqui::jqui_resizable(
                shiny::plotOutput(ns("correlation_ruv_corrected"), height = "400px")
              )
            )
          )
        ),
        
        # NEW: Post-Normalisation QC Tab (Filtering Summary) - BEFORE Correlation Filtering
        shiny::tabPanel(
          "Post-Normalisation QC",
          icon = shiny::icon("chart-bar"),
          shiny::br(),
          shiny::fluidRow(
            shiny::column(12,
              shiny::h5("Filtering Progression Summary", style = "text-align: center;"),
              shiny::helpText("Summary of protein filtering steps through normalization and RUV correction. This step removes proteins with excessive missing values after RUV processing."),
              shinyjqui::jqui_resizable(
                shiny::plotOutput(ns("post_norm_filtering_summary"), height = "600px")
              ),
              shiny::br(),
              shiny::verbatimTextOutput(ns("filtering_summary_text"))
            )
          )
        ),
        
        # NEW: Correlation Filtering Tab (Chunk 28) - FINAL step before DE
        shiny::tabPanel(
          "Correlation Filtering",
          icon = shiny::icon("filter"),
          shiny::br(),
          shiny::fluidRow(
            shiny::column(4,
              shiny::wellPanel(
                shiny::h5("Sample Correlation Filtering"),
                shiny::helpText("FINAL STEP: Filter samples based on their correlation with other samples in the same group. Low correlation indicates potential technical issues or sample quality problems. This is the last filtering step before differential expression analysis."),
                
                shiny::sliderInput(
                  ns("min_pearson_correlation_threshold"),
                  "Min. Pearson Correlation Threshold:",
                  min = 0.0,
                  max = 1.0,
                  value = 0.5,
                  step = 0.05,
                  width = "100%"
                ),
                shiny::helpText("Samples with correlation below this threshold will be removed from the analysis"),
                
                shiny::br(),
                
                shiny::actionButton(
                  ns("apply_correlation_filter"),
                  "Apply Correlation Filter",
                  class = "btn-primary",
                  width = "100%",
                  icon = shiny::icon("play")
                ),
                
                shiny::br(),
                shiny::br(),
                
                shiny::h5("Filter Results"),
                shiny::verbatimTextOutput(ns("correlation_filter_summary"))
              )
            ),
            shiny::column(8,
              shiny::h5("Final Filtered Data QC", style = "text-align: center;"),
              shiny::helpText("Quality control plots for the final dataset after correlation filtering"),
              shinyjqui::jqui_resizable(
                shiny::plotOutput(ns("final_qc_plot"), height = "500px")
              )
            )
          )
        )
      )
    )
  )
}

#' @rdname normalizationAppletModule 
#' @export
normalizationAppletServer <- function(id, workflow_data, experiment_paths, omic_type, experiment_label, selected_tab = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    
    message("=== NORMALIZATION MODULE SERVER STARTED ===")
    message(sprintf("Module ID: %s", id))
    message(sprintf("workflow_data is NULL: %s", is.null(workflow_data)))
    if (!is.null(workflow_data$state_manager)) {
      message(sprintf("Current state at module start: %s", workflow_data$state_manager$current_state))
    }
    
    # Initialize reactive values for normalization state
    norm_data <- reactiveValues(
      pre_norm_qc_generated = FALSE,
      normalization_complete = FALSE,
      ruv_complete = FALSE,
      correlation_filtering_complete = FALSE,
      
      # QC plot objects - 3 columns: post-filtering, post-normalization, ruv-corrected
      qc_plots = list(
        post_filtering = list(
          pca = NULL,
          density = NULL, 
          rle = NULL,
          correlation = NULL
        ),
        post_normalization = list(
          pca = NULL,
          density = NULL,
          rle = NULL, 
          correlation = NULL
        ),
        ruv_corrected = list(
          pca = NULL,
          density = NULL,
          rle = NULL,
          correlation = NULL
        )
      ),
      
      # Normalization results
      normalized_protein_obj = NULL,
      ruv_normalized_obj = NULL,
      correlation_filtered_obj = NULL,
      best_k = NULL,
      control_genes_index = NULL,
      
      # Correlation filtering results
      correlation_vector = NULL,
      correlation_threshold = NULL,
      final_qc_plot = NULL,
      final_filtering_plot = NULL,
      
      # Post-normalisation filtering summary
      post_norm_filtering_plot = NULL,
      filtering_summary_text = NULL
    )
    
    # Helper function to get current plot aesthetics
    getPlotAesthetics <- function() {
      list(
        color_var = if(is.null(input$color_variable) || input$color_variable == "") "group" else input$color_variable,
        shape_var = if(is.null(input$shape_variable) || input$shape_variable == "") "group" else input$shape_variable
      )
    }
    
    # Helper function to generate pre-normalization QC plots
    generatePreNormalizationQc <- function() {
      message("=== GENERATING PRE-NORMALIZATION QC PLOTS ===")
      message("Step 1: Getting current S4 object...")
      
      # Get current S4 object using CORRECT R6 methods
      shiny::req(workflow_data$state_manager)
      current_state <- workflow_data$state_manager$current_state
      message(sprintf("Current state: '%s'", current_state))
      
      current_s4 <- workflow_data$state_manager$getState(current_state)
      message(sprintf("S4 object is NULL: %s", is.null(current_s4)))
      
      if (!is.null(current_s4)) {
        message(sprintf("S4 object class: %s", class(current_s4)))
      }
      
      if (is.null(current_s4)) {
        stop("No S4 object available for QC plot generation")
      }
      
      message("Step 2: S4 object retrieved successfully")
      
      # Get current plot aesthetics
      aesthetics <- getPlotAesthetics()
      
      # Generate QC plots and store in norm_data$qc_plots$post_filtering
      
      # PCA plot
      norm_data$qc_plots$post_filtering$pca <- plotPca(
        current_s4,
        grouping_variable = aesthetics$color_var,
        label_column = "",
        shape_variable = aesthetics$shape_var, 
        title = "Pre-Normalization PCA",
        font_size = 8
      )
      
      # RLE plot  
      norm_data$qc_plots$post_filtering$rle <- plotRle(
        current_s4,
        group = aesthetics$color_var,
        yaxis_limit = c(-6, 6)
      )
      
      # Density plot
      norm_data$qc_plots$post_filtering$density <- plotDensity(
        norm_data$qc_plots$post_filtering$pca,
        grouping_variable = aesthetics$color_var
      )
      
      # Correlation plot
      norm_data$qc_plots$post_filtering$correlation <- plotPearson(
        current_s4,
        tech_rep_remove_regex = "pool",
        correlation_group = aesthetics$color_var
      )
      
      message("Pre-normalization QC plots generated successfully")
    }
    
    # Helper function to generate post-normalization QC plots
    generatePostNormalizationQc <- function(normalized_s4) {
      aesthetics <- getPlotAesthetics()
      
      norm_data$qc_plots$post_normalization$pca <- plotPca(
        normalized_s4,
        grouping_variable = aesthetics$color_var,
        label_column = "",
        shape_variable = aesthetics$shape_var,
        title = "Post-Normalization PCA", 
        font_size = 8
      )
      
      norm_data$qc_plots$post_normalization$rle <- plotRle(
        normalized_s4,
        group = aesthetics$color_var,
        yaxis_limit = c(-6, 6)
      )
      
      norm_data$qc_plots$post_normalization$density <- plotDensity(
        norm_data$qc_plots$post_normalization$pca,
        grouping_variable = aesthetics$color_var
      )
      
      norm_data$qc_plots$post_normalization$correlation <- plotPearson(
        normalized_s4,
        tech_rep_remove_regex = "pool",
        correlation_group = aesthetics$color_var
      )
    }
    
    # Helper function to generate RUV-corrected QC plots
    generateRuvCorrectedQc <- function(ruv_corrected_s4) {
      aesthetics <- getPlotAesthetics()
      
      norm_data$qc_plots$ruv_corrected$pca <- plotPca(
        ruv_corrected_s4,
        grouping_variable = aesthetics$color_var,
        label_column = "",
        shape_variable = aesthetics$shape_var,
        title = "RUV-Corrected PCA",
        font_size = 8
      )
      
      norm_data$qc_plots$ruv_corrected$rle <- plotRle(
        ruv_corrected_s4,
        group = aesthetics$color_var, 
        yaxis_limit = c(-6, 6)
      )
      
      norm_data$qc_plots$ruv_corrected$density <- plotDensity(
        norm_data$qc_plots$ruv_corrected$pca,
        grouping_variable = aesthetics$color_var
      )
      
      norm_data$qc_plots$ruv_corrected$correlation <- plotPearson(
        ruv_corrected_s4,
        tech_rep_remove_regex = "pool",
        correlation_group = aesthetics$color_var
      )
    }
    
    # Update plot aesthetic choices based on design matrix
    observe({
      if (!is.null(workflow_data$design_matrix)) {
        design_cols <- colnames(workflow_data$design_matrix)
        
        # Filter to common experimental variables
        available_vars <- intersect(design_cols, c("group", "factor1", "factor2", "technical_replicate_id", "sample_id"))
        
        if (length(available_vars) > 0) {
          # Update color variable choices
          shiny::updateSelectInput(session, "color_variable",
            choices = available_vars,
            selected = if("group" %in% available_vars) "group" else available_vars[1]
          )
          
          # Update shape variable choices  
          shiny::updateSelectInput(session, "shape_variable",
            choices = available_vars,
            selected = if("group" %in% available_vars) "group" else available_vars[1]
          )
        }
      }
    })
    
    # Auto-trigger pre-normalization QC when normalization tab is clicked
    # Observe tab selection passed from parent workflow
    if (!is.null(selected_tab)) {
      observeEvent(selected_tab(), {
        
        # Only trigger if normalization tab is selected
        if (!is.null(selected_tab()) && selected_tab() == "normalization") {
          
          message("=== NORMALIZATION TAB CLICKED ===")
          message(sprintf("workflow_data$state_manager is NULL: %s", is.null(workflow_data$state_manager)))
          
          if (!is.null(workflow_data$state_manager)) {
            current_state <- workflow_data$state_manager$current_state
            
            # Define valid states where this tab can be active
            valid_states_for_norm_tab <- c("protein_replicate_filtered", "ruv_corrected", "correlation_filtered")
            
            message(sprintf("Current state: '%s'", current_state))
            message(sprintf("Target trigger state: 'protein_replicate_filtered'"))
            message(sprintf("Pre-norm QC generated: %s", norm_data$pre_norm_qc_generated))
            
            # Auto-trigger only fires if we've just completed QC and haven't run pre-norm QC yet
            if (current_state == "protein_replicate_filtered" && !norm_data$pre_norm_qc_generated) {
              
              message("*** AUTO-TRIGGERING PRE-NORMALIZATION QC (chunk 24) ***")
              
              tryCatch({
                # Generate pre-normalization QC plots (chunk 24 equivalent)
                generatePreNormalizationQc()
                norm_data$pre_norm_qc_generated <- TRUE
                
                message("*** PRE-NORMALIZATION QC COMPLETED SUCCESSFULLY ***")
                
              }, error = function(e) {
                message(paste("*** ERROR generating pre-normalization QC:", e$message))
                shiny::showNotification(
                  paste("Error generating pre-normalization QC:", e$message),
                  type = "error",
                  duration = 10
                )
              })
              
            } else if (current_state %in% valid_states_for_norm_tab) {
              # If we are in a later, valid state (e.g., correlation_filtered),
              # or have already run the QC, we don't need to do anything.
              message(sprintf("State is '%s'. Skipping auto-trigger as it's not applicable or already done.", current_state))
              
            } else {
              # This case handles when the user has not completed the required prior steps
              message(sprintf("*** State '%s' is not valid for normalization. User needs to complete QC. ***", current_state))
              shiny::showNotification(
                "Please complete the Quality Control filtering steps before accessing normalization.",
                type = "warning",
                duration = 5
              )
            }
          } else {
            message("*** workflow_data$state_manager is NULL - cannot check state ***")
          }
        }
      }, ignoreInit = TRUE)
    }
    
    # Regenerate plots when aesthetic variables change
    observeEvent(c(input$color_variable, input$shape_variable), {
      # Only regenerate if we have plots already
      if (norm_data$pre_norm_qc_generated) {
        message("Regenerating pre-normalization QC plots with new aesthetics...")
        tryCatch({
          generatePreNormalizationQc()
        }, error = function(e) {
          message(paste("Error regenerating pre-normalization QC:", e$message))
        })
      }
      
      # Also regenerate post-normalization plots if they exist
      if (norm_data$normalization_complete && !is.null(norm_data$normalized_protein_obj)) {
        message("Regenerating post-normalization QC plots with new aesthetics...")
        tryCatch({
          generatePostNormalizationQc(norm_data$normalized_protein_obj)
        }, error = function(e) {
          message(paste("Error regenerating post-normalization QC:", e$message))
        })
      }
      
      # Also regenerate RUV-corrected plots if they exist
      if (norm_data$ruv_complete && !is.null(norm_data$ruv_normalized_obj)) {
        message("Regenerating RUV-corrected QC plots with new aesthetics...")
        tryCatch({
          generateRuvCorrectedQc(norm_data$ruv_normalized_obj)
        }, error = function(e) {
          message(paste("Error regenerating RUV-corrected QC:", e$message))
        })
      }
    })
    
    # Normalization button logic
    observeEvent(input$run_normalization, {
      message("=== NORMALIZATION BUTTON CLICKED ===")
      message("Starting normalization workflow...")
      
      tryCatch({
        # Get current S4 object using CORRECT R6 methods
        shiny::req(workflow_data$state_manager)
        current_state <- workflow_data$state_manager$current_state
        current_s4 <- workflow_data$state_manager$getState(current_state)
        
        if (is.null(current_s4)) {
          stop("No S4 object available for normalization")
        }
        
        # Update progress
        shiny::withProgress(message = "Running normalization...", value = 0, {
          
          # Step 1: Between-samples normalization (chunk 25 equivalent)
          shiny::incProgress(0.2, detail = "Normalizing between samples...")
          message("*** STEP 1: Starting between-samples normalization ***")
          
          tryCatch({
            # Update config with user selection
            if (exists("config_list", envir = .GlobalEnv)) {
              config_list$normaliseBetweenSamples$method <- input$norm_method
              assign("config_list", config_list, envir = .GlobalEnv)
            }
            
            normalized_s4 <- normaliseBetweenSamples(current_s4, normalisation_method = input$norm_method)
            norm_data$normalized_protein_obj <- normalized_s4
            message("*** STEP 1: Between-samples normalization completed ***")
            
          }, error = function(e) {
            stop(paste("Step 1 (normalization) error:", e$message))
          })
          
          # Generate post-normalization QC plots
          shiny::incProgress(0.2, detail = "Generating post-normalization QC plots...")
          message("*** STEP 2: Generating post-normalization QC plots ***")
          
          tryCatch({
            generatePostNormalizationQc(normalized_s4)
            norm_data$normalization_complete <- TRUE
            message("*** STEP 2: Post-normalization QC completed ***")
            
          }, error = function(e) {
            stop(paste("Step 2 (post-norm QC) error:", e$message))
          })
          
          # Step 3: RUV-III correction (chunks 26-27 equivalent)
          shiny::incProgress(0.2, detail = "Determining RUV parameters...")
          message("*** STEP 3: Starting RUV parameter determination ***")
          
          tryCatch({
            # Get RUV parameters based on mode
            if (input$ruv_mode == "automatic") {
              # Use default parameters (placeholder for automatic tuning)
              percentage_as_neg_ctrl <- 5
              ruv_k <- 3  # Default, will be optimized
            } else {
              percentage_as_neg_ctrl <- input$ruv_percentage
              ruv_k <- if(is.null(input$ruv_k) || is.na(input$ruv_k)) 3 else input$ruv_k
            }
            
            # Get negative control proteins
            message("*** STEP 3a: Getting negative control proteins ***")
            control_genes_index <- getNegCtrlProtAnova(
              normalized_s4,
              ruv_grouping_variable = getPlotAesthetics()$color_var,
              percentage_as_neg_ctrl = percentage_as_neg_ctrl,
              ruv_qval_cutoff = 0.05,
              ruv_fdr_method = "BH"
            )
            norm_data$control_genes_index <- control_genes_index
            message("*** STEP 3a: Negative control proteins completed ***")
            
          }, error = function(e) {
            stop(paste("Step 3a (negative controls) error:", e$message))
          })
          
          # Determine best k using canonical correlation
          tryCatch({
            if (input$ruv_mode == "automatic") {
              message("*** STEP 3b: Determining best k using canonical correlation ***")
              cancor_result <- ruvCancor(
                normalized_s4,
                ctrl = control_genes_index,
                num_components_to_impute = 2,
                ruv_grouping_variable = getPlotAesthetics()$color_var
              )
              # Find best k (simplified - take the one with highest canonical correlation)
              best_k <- which.max(cancor_result) + 1  # Add 1 since index starts at 0
              norm_data$best_k <- best_k
              message("*** STEP 3b: Canonical correlation completed ***")
            } else {
              best_k <- ruv_k
              norm_data$best_k <- best_k
              message("*** STEP 3b: Using manual k value ***")
            }
            
          }, error = function(e) {
            stop(paste("Step 3b (canonical correlation) error:", e$message))
          })
          
          # Apply RUV-III correction
          shiny::incProgress(0.2, detail = "Applying RUV-III batch correction...")
          message("*** STEP 4: Applying RUV-III correction ***")
          
          tryCatch({
            ruv_corrected_s4 <- ruvIII_C_Varying(
              normalized_s4,
              ruv_grouping_variable = getPlotAesthetics()$color_var,
              ruv_number_k = best_k,
              ctrl = control_genes_index
            )
            message("*** STEP 4: RUV-III correction completed ***")
            
          }, error = function(e) {
            stop(paste("Step 4 (RUV-III correction) error:", e$message))
          })
          
          # Clean up any proteins with excessive missing values after RUV
          message("*** STEP 5: Cleaning up missing values ***")
          
          tryCatch({
            ruv_corrected_s4_clean <- removeRowsWithMissingValuesPercent(
              theObject = ruv_corrected_s4
            )
            
            # Log protein count after RUV filtering (matching RMarkdown chunk 27)
            ruvfilt_protein_count <- ruv_corrected_s4_clean@protein_quant_table |>
              dplyr::distinct(Protein.Ids) |>
              nrow()
            message(sprintf("Number of distinct proteins remaining after RUV normalization and filtering: %d", ruvfilt_protein_count))
            
            # Update protein filtering tracking (matching RMarkdown chunk 27)
            filtering_plot <- updateProteinFiltering(
              data = ruv_corrected_s4_clean@protein_quant_table,
              step_name = "11_RUV_filtered",
              omic_type = omic_type,
              experiment_label = experiment_label,
              return_grid = TRUE,
              overwrite = TRUE
            )
            
            # Store the filtering summary plot for display
            norm_data$post_norm_filtering_plot <- filtering_plot
            
            # Generate filtering summary text
            norm_data$filtering_summary_text <- sprintf(
              "Filtering Summary through Normalization & RUV:\n\n• Pre-normalization proteins: [Previous step]\n• Post-RUV filtering: %d proteins\n• Proteins removed by RUV: [Calculated from difference]\n\nRUV Parameters:\n• Normalization method: %s\n• RUV mode: %s\n• RUV k value: %d\n• Negative control %%: %.1f",
              ruvfilt_protein_count,
              input$norm_method,
              input$ruv_mode,
              best_k,
              percentage_as_neg_ctrl
            )
            
            norm_data$ruv_normalized_obj <- ruv_corrected_s4_clean
            message("*** STEP 5: Missing values cleanup completed ***")
            
            # Save post-normalization state to R6 state manager (includes RUV + missing value filtering)
            workflow_data$state_manager$saveState(
              state_name = "ruv_corrected",
              s4_data_object = ruv_corrected_s4_clean,
              config_object = list(
                norm_method = input$norm_method,
                ruv_mode = input$ruv_mode,
                ruv_k = best_k,
                percentage_as_neg_ctrl = percentage_as_neg_ctrl
              ),
              description = "Post-normalization complete: RUV-III correction and missing value cleanup completed"
            )
            
          }, error = function(e) {
            stop(paste("Step 5 (missing values cleanup) error:", e$message))
          })
          
          # Generate RUV-corrected QC plots
          shiny::incProgress(0.2, detail = "Generating RUV-corrected QC plots...")
          message("*** STEP 6: Generating RUV-corrected QC plots ***")
          
          tryCatch({
            generateRuvCorrectedQc(ruv_corrected_s4_clean)
            norm_data$ruv_complete <- TRUE
            message("*** STEP 6: RUV-corrected QC plots completed ***")
            
          }, error = function(e) {
            stop(paste("Step 6 (RUV-corrected QC) error:", e$message))
          })
          
          # CRITICAL FIX: Only enable correlation filtering tab, NOT differential expression
          message("*** STEP 7: Enabling correlation filtering step ***")
          # Remove premature DE enablement - correlation filtering must happen first
          message("*** STEP 7: Normalization and RUV workflow completed - ready for correlation filtering ***")
          
        })
        
        shiny::showNotification(
          "Normalization and RUV correction completed! Check the 'Post-Normalisation QC' tab for filtering summary, then proceed to 'Correlation Filtering' tab for the final step.",
          type = "success",
          duration = 10
        )
        
        message("Normalization workflow completed successfully")
        
      }, error = function(e) {
        message(paste("Error in normalization workflow:", e$message))
        shiny::showNotification(
          paste("Error in normalization:", e$message),
          type = "error",
          duration = 10
        )
      })
    })
    
    # NEW: Correlation Filtering Button Logic (Chunk 28)
    observeEvent(input$apply_correlation_filter, {
      message("=== CORRELATION FILTERING BUTTON CLICKED ===")
      
      tryCatch({
        # Get RUV-corrected object (should exist after normalization)
        ruv_s4 <- norm_data$ruv_normalized_obj
        if (is.null(ruv_s4)) {
          stop("RUV correction must be completed before correlation filtering")
        }
        
        shiny::withProgress(message = "Applying correlation filter...", value = 0, {
          
          # Step 1: Calculate correlation vector (chunk 28)
          shiny::incProgress(0.3, detail = "Calculating sample correlations...")
          message("*** CORRELATION STEP 1: Calculating sample correlations ***")
          
          correlation_vec <- pearsonCorForSamplePairs(
            ruv_s4,
            tech_rep_remove_regex = "pool",
            correlation_group = getPlotAesthetics()$color_var
          )
          norm_data$correlation_vector <- correlation_vec
          norm_data$correlation_threshold <- input$min_pearson_correlation_threshold
          message("*** CORRELATION STEP 1: Sample correlations calculated ***")
          
          # Step 2: Apply correlation threshold filter (chunk 28)
          shiny::incProgress(0.4, detail = "Filtering low-correlation samples...")
          message("*** CORRELATION STEP 2: Applying correlation threshold filter ***")
          
          final_s4_for_de <- filterSamplesByProteinCorrelationThreshold(
            ruv_s4,
            pearson_correlation_per_pair = correlation_vec,
            min_pearson_correlation_threshold = input$min_pearson_correlation_threshold
          )
          norm_data$correlation_filtered_obj <- final_s4_for_de
          message("*** CORRELATION STEP 2: Correlation filtering applied ***")
          
          # Step 3: Update protein filtering tracking (chunk 28)
          shiny::incProgress(0.2, detail = "Updating tracking...")
          message("*** CORRELATION STEP 3: Updating protein filtering tracking ***")
          
          # Fix Issue 3: Capture the filtering plot for final QC display
          final_filtering_plot <- updateProteinFiltering(
            data = final_s4_for_de@protein_quant_table,
            step_name = "12_correlation_filtered",
            omic_type = omic_type,
            experiment_label = experiment_label,
            return_grid = TRUE,
            overwrite = TRUE
          )
          
          # Store the complete filtering progression for final QC display
          norm_data$final_filtering_plot <- final_filtering_plot
          
          # Generate final QC plot
          tryCatch({
            norm_data$final_qc_plot <- plotPca(
              final_s4_for_de,
              grouping_variable = getPlotAesthetics()$color_var,
              label_column = "",
              shape_variable = getPlotAesthetics()$shape_var,
              title = "Final Correlation-Filtered Data",
              font_size = 8
            )
          }, error = function(e) {
            message(paste("Error generating final QC plot:", e$message))
          })
          
          # Step 4: Save final results as TSV and RDS (chunk 28)
          shiny::incProgress(0.1, detail = "Saving results...")
          message("*** CORRELATION STEP 4: Saving final results ***")
          
          tryCatch({
            # Save TSV file
            if (!is.null(experiment_paths) && "protein_qc_dir" %in% names(experiment_paths)) {
              vroom::vroom_write(
                final_s4_for_de@protein_quant_table,
                file.path(experiment_paths$protein_qc_dir, "ruv_normalised_results_cln_with_replicates.tsv")
              )
              
              # Save RDS file
              saveRDS(
                final_s4_for_de,
                file.path(experiment_paths$protein_qc_dir, "ruv_normalised_results_cln_with_replicates.RDS")
              )
            }
          }, error = function(e) {
            message(paste("Warning: Could not save files:", e$message))
          })
          
        })
        
        # CRITICAL: THIS is where the final object becomes ready for DE
        workflow_data$ruv_normalised_for_de_analysis_obj <- final_s4_for_de
        
        # Save to R6 state manager with new state
        workflow_data$state_manager$saveState(
          state_name = "correlation_filtered",
          s4_data_object = final_s4_for_de,
          config_object = list(min_pearson_correlation_threshold = input$min_pearson_correlation_threshold),
          description = "Applied final sample correlation filter (chunk 28)"
        )
        
        # DEBUG66: CRITICAL STATE UPDATE TRIGGER
        cat("--- Entering STATE UPDATE TRIGGER setting ---\n")
        cat("   STATE UPDATE Step: Setting workflow_data$state_update_trigger...\n")
        old_trigger_value <- workflow_data$state_update_trigger
        new_trigger_value <- Sys.time()
        workflow_data$state_update_trigger <- new_trigger_value
        cat(sprintf("   STATE UPDATE Step: Old trigger value = %s\n", old_trigger_value))
        cat(sprintf("   STATE UPDATE Step: New trigger value = %s\n", new_trigger_value))
        cat("   STATE UPDATE Step: state_update_trigger SET SUCCESSFULLY\n")
        
        # NOW enable differential expression tab
        cat("   STATE UPDATE Step: Setting tab status...\n")
        workflow_data$tab_status$normalization <- "complete"
        workflow_data$tab_status$differential_expression <- "pending"
        cat(sprintf("   STATE UPDATE Step: normalization status = %s\n", workflow_data$tab_status$normalization))
        cat(sprintf("   STATE UPDATE Step: differential_expression status = %s\n", workflow_data$tab_status$differential_expression))
        cat("--- Exiting STATE UPDATE TRIGGER setting ---\n")
        
        # Update summary display
        correlation_summary <- sprintf(
          "Correlation filtering completed successfully!\n\nThreshold: %.2f\nProteins remaining: %d\nSamples remaining: %d\n\nReady for differential expression analysis.",
          input$min_pearson_correlation_threshold,
          length(unique(final_s4_for_de@protein_quant_table$Protein.Ids)),
          # Fix Issue 2: Use correct column counting for samples 
          # Sample columns are all columns except the protein ID column
          length(setdiff(colnames(final_s4_for_de@protein_quant_table), final_s4_for_de@protein_id_column))
        )
        output$correlation_filter_summary <- shiny::renderText(correlation_summary)
        
        norm_data$correlation_filtering_complete <- TRUE
        
        shiny::showNotification(
          "Correlation filtering completed! Ready for differential expression analysis.",
          type = "success",
          duration = 5
        )
        
        message("=== CORRELATION FILTERING COMPLETED SUCCESSFULLY ===")
        
      }, error = function(e) {
        message(paste("Error in correlation filtering:", e$message))
        shiny::showNotification(
          paste("Error in correlation filtering:", e$message),
          type = "error",
          duration = 10
        )
      })
    })
    
    # Render QC plots - Post-filtering column
    output$pca_post_filtering <- shiny::renderPlot({
      if (!is.null(norm_data$qc_plots$post_filtering$pca)) {
        norm_data$qc_plots$post_filtering$pca
      } else {
        plot.new()
        text(0.5, 0.5, "Pre-normalization QC not yet generated", cex = 1.2)
      }
    })
    
    output$density_post_filtering <- shiny::renderPlot({
      if (!is.null(norm_data$qc_plots$post_filtering$density)) {
        norm_data$qc_plots$post_filtering$density
      } else {
        plot.new()
        text(0.5, 0.5, "Pre-normalization QC not yet generated", cex = 1.2)
      }
    })
    
    output$rle_post_filtering <- shiny::renderPlot({
      if (!is.null(norm_data$qc_plots$post_filtering$rle)) {
        norm_data$qc_plots$post_filtering$rle
      } else {
        plot.new()
        text(0.5, 0.5, "Pre-normalization QC not yet generated", cex = 1.2)
      }
    })
    
    output$correlation_post_filtering <- shiny::renderPlot({
      if (!is.null(norm_data$qc_plots$post_filtering$correlation)) {
        norm_data$qc_plots$post_filtering$correlation
      } else {
        plot.new()
        text(0.5, 0.5, "Pre-normalization QC not yet generated", cex = 1.2)
      }
    })
    
    # Render QC plots - Post-normalization column
    output$pca_post_normalization <- shiny::renderPlot({
      if (!is.null(norm_data$qc_plots$post_normalization$pca)) {
        norm_data$qc_plots$post_normalization$pca
      } else {
        plot.new()
        text(0.5, 0.5, "Run normalization to generate plots", cex = 1.2)
      }
    })
    
    output$density_post_normalization <- shiny::renderPlot({
      if (!is.null(norm_data$qc_plots$post_normalization$density)) {
        norm_data$qc_plots$post_normalization$density
      } else {
        plot.new()
        text(0.5, 0.5, "Run normalization to generate plots", cex = 1.2)
      }
    })
    
    output$rle_post_normalization <- shiny::renderPlot({
      if (!is.null(norm_data$qc_plots$post_normalization$rle)) {
        norm_data$qc_plots$post_normalization$rle
      } else {
        plot.new()
        text(0.5, 0.5, "Run normalization to generate plots", cex = 1.2)
      }
    })
    
    output$correlation_post_normalization <- shiny::renderPlot({
      if (!is.null(norm_data$qc_plots$post_normalization$correlation)) {
        norm_data$qc_plots$post_normalization$correlation
      } else {
        plot.new()
        text(0.5, 0.5, "Run normalization to generate plots", cex = 1.2)
      }
    })
    
    # Render QC plots - RUV-corrected column
    output$pca_ruv_corrected <- shiny::renderPlot({
      if (!is.null(norm_data$qc_plots$ruv_corrected$pca)) {
        norm_data$qc_plots$ruv_corrected$pca
      } else {
        plot.new()
        text(0.5, 0.5, "Run RUV correction to generate plots", cex = 1.2)
      }
    })
    
    output$density_ruv_corrected <- shiny::renderPlot({
      if (!is.null(norm_data$qc_plots$ruv_corrected$density)) {
        norm_data$qc_plots$ruv_corrected$density
      } else {
        plot.new()
        text(0.5, 0.5, "Run RUV correction to generate plots", cex = 1.2)
      }
    })
    
    output$rle_ruv_corrected <- shiny::renderPlot({
      if (!is.null(norm_data$qc_plots$ruv_corrected$rle)) {
        norm_data$qc_plots$ruv_corrected$rle
      } else {
        plot.new()
        text(0.5, 0.5, "Run RUV correction to generate plots", cex = 1.2)
      }
    })
    
    output$correlation_ruv_corrected <- shiny::renderPlot({
      if (!is.null(norm_data$qc_plots$ruv_corrected$correlation)) {
        norm_data$qc_plots$ruv_corrected$correlation
      } else {
        plot.new()
        text(0.5, 0.5, "Run RUV correction to generate plots", cex = 1.2)
      }
    })
    
    # Reset normalization button logic
    observeEvent(input$reset_normalization, {
      message("Resetting normalization...")
      
      # Reset normalization state
      norm_data$normalization_complete <- FALSE
      norm_data$ruv_complete <- FALSE
      norm_data$correlation_filtering_complete <- FALSE
      norm_data$normalized_protein_obj <- NULL
      norm_data$ruv_normalized_obj <- NULL
      norm_data$correlation_filtered_obj <- NULL
      norm_data$best_k <- NULL
      norm_data$control_genes_index <- NULL
      norm_data$correlation_vector <- NULL
      norm_data$correlation_threshold <- NULL
      norm_data$final_qc_plot <- NULL
      norm_data$final_filtering_plot <- NULL
      norm_data$post_norm_filtering_plot <- NULL
      norm_data$filtering_summary_text <- NULL
      
      # Clear post-normalization, RUV, and correlation plots
      norm_data$qc_plots$post_normalization <- list(pca = NULL, density = NULL, rle = NULL, correlation = NULL)
      norm_data$qc_plots$ruv_corrected <- list(pca = NULL, density = NULL, rle = NULL, correlation = NULL)
      
      # Reset workflow data
      workflow_data$ruv_normalised_for_de_analysis_obj <- NULL
      workflow_data$tab_status$normalization <- "pending"
      workflow_data$tab_status$differential_expression <- "disabled"
      
      # Clear correlation filter summary and filtering summary
      output$correlation_filter_summary <- shiny::renderText("No correlation filtering applied yet")
      output$filtering_summary_text <- shiny::renderText("Filtering summary will be available after normalization and RUV correction.")
      
      shiny::showNotification(
        "Normalization has been reset to pre-normalization state",
        type = "warning",
        duration = 3
      )
      
      message("Normalization reset completed")
    })
    
    # NEW: Render post-normalisation filtering summary plot
    output$post_norm_filtering_summary <- shiny::renderPlot({
      if (!is.null(norm_data$post_norm_filtering_plot)) {
        # Fix Issue 1: Add proper grob rendering like qualityControlApplet
        grid::grid.draw(norm_data$post_norm_filtering_plot)
      } else {
        plot.new()
        text(0.5, 0.5, "Complete normalization and RUV correction\nto generate filtering summary", cex = 1.2)
      }
    })
    
    # NEW: Render filtering summary text
    output$filtering_summary_text <- shiny::renderText({
      if (!is.null(norm_data$filtering_summary_text)) {
        norm_data$filtering_summary_text
      } else {
        "Filtering summary will be available after normalization and RUV correction."
      }
    })
    
    # NEW: Render final QC plot after correlation filtering
    output$final_qc_plot <- shiny::renderPlot({
      # Fix Issue 3: Show both PCA and complete filtering progression
      if (!is.null(norm_data$final_qc_plot) && !is.null(norm_data$final_filtering_plot)) {
        # Create a combined layout: PCA on top, filtering progression below
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 1, heights = c(0.4, 0.6))))
        
        # Top: PCA plot
        grid::pushViewport(grid::viewport(layout.pos.row = 1))
        grid::grid.draw(ggplotGrob(norm_data$final_qc_plot))
        grid::popViewport()
        
        # Bottom: Complete filtering progression 
        grid::pushViewport(grid::viewport(layout.pos.row = 2))
        grid::grid.draw(norm_data$final_filtering_plot)
        grid::popViewport()
        
        grid::popViewport()
      } else if (!is.null(norm_data$final_filtering_plot)) {
        # Show just filtering progression if PCA failed
        grid::grid.draw(norm_data$final_filtering_plot)
      } else if (!is.null(norm_data$final_qc_plot)) {
        # Show just PCA if filtering progression failed
        norm_data$final_qc_plot
      } else {
        plot.new()
        text(0.5, 0.5, "Apply correlation filter to generate final QC plot", cex = 1.2)
      }
    })
    
    # Return normalization data for potential use by parent module
    return(norm_data)
  })
} 