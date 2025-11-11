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
          choices = c("group", "factor1", "factor2", "batch"),
          selected = "group",
          width = "100%"
        ),
        shiny::helpText("Variable to use for plot coloring"),
        
        shiny::selectInput(
          ns("shape_variable"),
          "Shape by:",
          choices = c("group", "factor1", "factor2", "batch"),
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
        
        # RUV grouping variable selection
        shiny::selectInput(
          ns("ruv_grouping_variable"),
          "RUV Grouping Variable:",
          choices = c("group", "factor1", "factor2", "batch"),
          selected = "group",
          width = "100%"
        ),
        shiny::helpText("Experimental factor to preserve during batch correction. RUV-III will remove unwanted variation while preserving differences between groups defined by this variable. For TMT workflows with batches, select 'batch' to remove batch effects."),
        
        shiny::hr(),
        
        # Automatic RUV controls (shown when automatic mode selected)
        shiny::conditionalPanel(
          condition = "input.ruv_mode == 'automatic'",
          ns = ns,
          
          shiny::h5("Automatic Optimization Parameters"),
          shiny::helpText("These parameters control the automatic optimization of RUV-III negative control selection. The function tests different percentages and finds the optimal balance between separation quality and k value."),
          
          # Percentage range control
          shiny::fluidRow(
            shiny::column(6,
              shiny::numericInput(
                ns("auto_percentage_min"),
                "Min %:",
                value = 1,
                min = 1,
                max = 50,
                width = "100%"
              )
            ),
            shiny::column(6,
              shiny::numericInput(
                ns("auto_percentage_max"), 
                "Max %:",
                value = 20,
                min = 1,
                max = 50,
                width = "100%"
              )
            )
          ),
          shiny::helpText("Range of percentages to test (1-50%). Start with 1-20% for most datasets. Larger datasets may benefit from higher ranges."),
          
          shiny::selectInput(
            ns("separation_metric"),
            "Separation Quality Metric:",
            choices = list(
              "Maximum Difference (recommended)" = "max_difference",
              "Mean Difference" = "mean_difference", 
              "Area Under Curve" = "auc",
              "Weighted Difference" = "weighted_difference"
            ),
            selected = "max_difference",
            width = "100%"
          ),
          shiny::helpText("Method for evaluating separation between All and Control groups in canonical correlation plots. Max difference is most robust for typical datasets."),
          
          shiny::sliderInput(
            ns("k_penalty_weight"),
            "K Value Penalty Weight:",
            min = 0.1,
            max = 0.9,
            value = 0.5,
            step = 0.1,
            width = "100%"
          ),
          shiny::helpText("Controls how strongly high k values are penalized (0.1 = lenient, 0.9 = strict). Higher values favor lower k to preserve biological signal."),
          
          shiny::checkboxInput(
            ns("adaptive_k_penalty"),
            "Adaptive K Penalty (sample size aware)",
            value = TRUE,
            width = "100%"
          ),
          shiny::helpText("Automatically adjusts penalty thresholds based on dataset size. Recommended for most analyses."),
          
          shiny::br()
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
        ),
        
        shiny::br(),
        shiny::br(),
        
        # Export session button
        shiny::actionButton(
          ns("export_filtered_session"),
          "Export Filtered Data Session",
          class = "btn-success",
          width = "100%",
          icon = shiny::icon("download")
        ),
        shiny::helpText("Save filtered data and contrasts for later DE analysis")
      )
    ),
    
    # QC plot panels with 3-column layout structure and correlation filtering tab
    shiny::column(9,
      shiny::tabsetPanel(
        id = ns("norm_qc_tabs"),
        
        # NEW: RUV QC Tab (first tab)
        shiny::tabPanel(
          "RUV QC",
          icon = shiny::icon("chart-line"),
          shiny::br(),
          shiny::fluidRow(
            shiny::column(12,
              shiny::h5("RUV Parameter Optimization Results", style = "text-align: center;"),
              shiny::helpText("Canonical correlation plot and optimization summary for RUV-III parameter selection. This plot shows the separation between 'All' and 'Control' protein groups across different k values."),
              
              shiny::fluidRow(
                shiny::column(8,
                  shinyjqui::jqui_resizable(
                    shiny::plotOutput(ns("ruv_canonical_correlation_plot"), height = "500px")
                  )
                ),
                shiny::column(4,
                  shiny::wellPanel(
                    shiny::h5("Optimization Summary"),
                    shiny::verbatimTextOutput(ns("ruv_optimization_summary")),
                    
                    shiny::br(),
                    
                    shiny::h5("Optimization Results"),
                    shiny::helpText("Detailed results for all tested percentages:"),
                    DT::dataTableOutput(ns("ruv_optimization_table"))
                  )
                )
              )
            )
          )
        ),
        
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
      
      # Initialize S4 QC composite figure object
      QC_composite_figure = NULL,
      
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
      filtering_summary_text = NULL,
      
      # RUV optimization results
      ruv_optimization_result = NULL
    )
    
    # Update plot aesthetic choices based on design matrix
    observe({
      if (!is.null(workflow_data$design_matrix)) {
        design_cols <- colnames(workflow_data$design_matrix)
        
        # Filter to common experimental variables for plot aesthetics (include batch)
        plot_available_vars <- intersect(design_cols, c("group", "factor1", "factor2", "batch", "technical_replicate_id", "sample_id"))
        
        if (length(plot_available_vars) > 0) {
          # Update color variable choices
          shiny::updateSelectInput(session, "color_variable",
            choices = plot_available_vars,
            selected = if("group" %in% plot_available_vars) "group" else plot_available_vars[1]
          )
          
          # Update shape variable choices  
          shiny::updateSelectInput(session, "shape_variable",
            choices = plot_available_vars,
            selected = if("group" %in% plot_available_vars) "group" else plot_available_vars[1]
          )
        }
        
        # Filter to experimental grouping variables for RUV (exclude technical IDs, include batch)
        ruv_available_vars <- intersect(design_cols, c("group", "factor1", "factor2", "batch"))
        
        if (length(ruv_available_vars) > 0) {
          # If batch exists and has values, suggest it as default for RUV
          default_ruv <- if("batch" %in% ruv_available_vars && 
                            any(!is.na(workflow_data$design_matrix$batch))) {
            message("*** RUV: Batch column detected with values - suggesting 'batch' as RUV variable to set to shapes ***")
            "group"
          } else if("group" %in% ruv_available_vars) {
            "group"
          } else {
            ruv_available_vars[1]
          }
          
          # Update RUV grouping variable choices
          shiny::updateSelectInput(session, "ruv_grouping_variable",
            choices = ruv_available_vars,
            selected = default_ruv
          )
          
          message(sprintf("*** RUV: Updated grouping variable choices to: %s ***", paste(ruv_available_vars, collapse = ", ")))
        }
      }
    })
    
    # Helper function to get current plot aesthetics
    getPlotAesthetics <- function() {
      list(
        color_var = if(is.null(input$color_variable) || input$color_variable == "") "group" else input$color_variable,
        shape_var = if(is.null(input$shape_variable) || input$shape_variable == "") "group" else input$shape_variable
      )
    }
    
    # Helper function to get current RUV grouping variable
    getRuvGroupingVariable <- function() {
      if(is.null(input$ruv_grouping_variable) || input$ruv_grouping_variable == "") {
        "group"  # Default fallback
      } else {
        input$ruv_grouping_variable
      }
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
      
      # PCA plot (Step 1/4)
      message("*** PRE-NORM QC: Generating PCA plot (Step 1/4) ***")
      pca_start_time <- Sys.time()
      pca_plot <- plotPca(
        current_s4,
        grouping_variable = aesthetics$color_var,
        label_column = "",
        shape_variable = aesthetics$shape_var, 
                    title = "",
        font_size = 8
      )
      norm_data$qc_plots$post_filtering$pca <- pca_plot
      pca_elapsed <- as.numeric(difftime(Sys.time(), pca_start_time, units = "secs"))
      message(sprintf("*** PRE-NORM QC: PCA plot completed in %.1f seconds ***", pca_elapsed))
      
      # RLE plot (Step 2/4)
      message("*** PRE-NORM QC: Generating RLE plot (Step 2/4) ***")
      rle_start_time <- Sys.time()
      rle_plot <- plotRle(
        current_s4,
        group = aesthetics$color_var,
        yaxis_limit = c(-6, 6)
      )
      norm_data$qc_plots$post_filtering$rle <- rle_plot
      rle_elapsed <- as.numeric(difftime(Sys.time(), rle_start_time, units = "secs"))
      message(sprintf("*** PRE-NORM QC: RLE plot completed in %.1f seconds ***", rle_elapsed))
      
      # Density plot (Step 3/4)
      message("*** PRE-NORM QC: Generating Density plot (Step 3/4) ***")
      density_start_time <- Sys.time()
      density_plot <- plotDensity(
        pca_plot,
        grouping_variable = aesthetics$color_var
      )
      norm_data$qc_plots$post_filtering$density <- density_plot
      density_elapsed <- as.numeric(difftime(Sys.time(), density_start_time, units = "secs"))
      message(sprintf("*** PRE-NORM QC: Density plot completed in %.1f seconds ***", density_elapsed))
      
      # Correlation plot (Step 4/4) - This is the slow one
      message("*** PRE-NORM QC: Generating Pearson correlation plot (Step 4/4) ***")
      num_samples <- length(setdiff(colnames(current_s4@protein_quant_table), current_s4@protein_id_column))
      estimated_pairs <- choose(num_samples, 2)
      message(sprintf("*** PRE-NORM QC: Sample count = %d, Expected pairs ≈ %d ***", num_samples, estimated_pairs))
      message("*** PRE-NORM QC: This may take 30-90 seconds for pairwise correlation calculations ***")
      
      pearson_start_time <- Sys.time()
      
      # Check if we're in a Shiny context and wrap with progress indicator
      if (shiny::isRunning()) {
        correlation_plot <- shiny::withProgress(
          message = "Generating Pearson correlation plot...", 
          detail = "Calculating pairwise correlations (30-90s)", 
          value = 0.5, {
            result <- plotPearson(
              current_s4,
              tech_rep_remove_regex = "pool",
              correlation_group = aesthetics$color_var
            )
            shiny::incProgress(0.5)
            result
          }
        )
      } else {
        correlation_plot <- plotPearson(
          current_s4,
          tech_rep_remove_regex = "pool",
          correlation_group = aesthetics$color_var
        )
      }
      
      norm_data$qc_plots$post_filtering$correlation <- correlation_plot
      pearson_elapsed <- as.numeric(difftime(Sys.time(), pearson_start_time, units = "secs"))
      message(sprintf("*** PRE-NORM QC: Pearson correlation completed in %.1f seconds ***", pearson_elapsed))
      
      # ONLY populate S4 object if it exists (during normalization workflow)
      if (!is.null(norm_data$QC_composite_figure)) {
        norm_data$QC_composite_figure@pca_plots$pca_plot_before_cyclic_loess_group <- pca_plot
        norm_data$QC_composite_figure@rle_plots$rle_plot_before_cyclic_loess_group <- rle_plot
        norm_data$QC_composite_figure@density_plots$density_plot_before_cyclic_loess_group <- density_plot
        norm_data$QC_composite_figure@pearson_plots$pearson_correlation_pair_before_cyclic_loess <- correlation_plot
      }
      
      message("Pre-normalization QC plots generated successfully")
    }
    
    # Helper function to generate post-normalization QC plots
    generatePostNormalizationQc <- function(normalized_s4) {
      aesthetics <- getPlotAesthetics()
      
      pca_plot <- plotPca(
        normalized_s4,
        grouping_variable = aesthetics$color_var,
        label_column = "",
        shape_variable = aesthetics$shape_var,
                    title = "", 
        font_size = 8
      )
      norm_data$qc_plots$post_normalization$pca <- pca_plot
      
      rle_plot <- plotRle(
        normalized_s4,
        group = aesthetics$color_var,
        yaxis_limit = c(-6, 6)
      )
      norm_data$qc_plots$post_normalization$rle <- rle_plot
      
      density_plot <- plotDensity(
        pca_plot,
        grouping_variable = aesthetics$color_var
      )
      norm_data$qc_plots$post_normalization$density <- density_plot
      
      correlation_plot <- plotPearson(
        normalized_s4,
        tech_rep_remove_regex = "pool",
        correlation_group = aesthetics$color_var
      )
      norm_data$qc_plots$post_normalization$correlation <- correlation_plot
      
      # Populate S4 object (should exist during normalization workflow)
      if (!is.null(norm_data$QC_composite_figure)) {
        norm_data$QC_composite_figure@pca_plots$pca_plot_before_ruvIIIc_group <- pca_plot
        norm_data$QC_composite_figure@rle_plots$rle_plot_before_ruvIIIc_group <- rle_plot
        norm_data$QC_composite_figure@density_plots$density_plot_before_ruvIIIc_group <- density_plot
        norm_data$QC_composite_figure@pearson_plots$pearson_correlation_pair_before_ruvIIIc <- correlation_plot
      }
    }
    
    # Helper function to generate RUV-corrected QC plots
    generateRuvCorrectedQc <- function(ruv_corrected_s4) {
      aesthetics <- getPlotAesthetics()
      
      pca_plot <- plotPca(
        ruv_corrected_s4,
        grouping_variable = aesthetics$color_var,
        label_column = "",
        shape_variable = aesthetics$shape_var,
                    title = "",
        font_size = 8
      )
      norm_data$qc_plots$ruv_corrected$pca <- pca_plot
      
      rle_plot <- plotRle(
        ruv_corrected_s4,
        group = aesthetics$color_var, 
        yaxis_limit = c(-6, 6)
      )
      norm_data$qc_plots$ruv_corrected$rle <- rle_plot
      
      density_plot <- plotDensity(
        pca_plot,
        grouping_variable = aesthetics$color_var
      )
      norm_data$qc_plots$ruv_corrected$density <- density_plot
      
      correlation_plot <- plotPearson(
        ruv_corrected_s4,
        tech_rep_remove_regex = "pool",
        correlation_group = aesthetics$color_var
      )
      norm_data$qc_plots$ruv_corrected$correlation <- correlation_plot
      
      # Populate S4 object (should exist during normalization workflow)
      if (!is.null(norm_data$QC_composite_figure)) {
        norm_data$QC_composite_figure@pca_plots$pca_plot_after_ruvIIIc_group <- pca_plot
        norm_data$QC_composite_figure@rle_plots$rle_plot_after_ruvIIIc_group <- rle_plot
        norm_data$QC_composite_figure@density_plots$density_plot_after_ruvIIIc_group <- density_plot
        norm_data$QC_composite_figure@pearson_plots$pearson_correlation_pair_after_ruvIIIc_group <- correlation_plot
      }
    }
    
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
        
        # Initialize the S4 QC composite figure object at the START
        message("*** INITIALIZING S4 QC COMPOSITE FIGURE ***")
        norm_data$QC_composite_figure <- InitialiseGrid()
        
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
            
            # Save post-normalization protein matrix (before RUV)
            tryCatch({
              if (!is.null(experiment_paths$protein_qc_dir)) {
                vroom::vroom_write(
                  normalized_s4@protein_quant_table,
                  file.path(experiment_paths$protein_qc_dir, "normalized_protein_matrix_pre_ruv.tsv")
                )
                message("*** STEP 1b: Saved post-normalization protein matrix to protein_qc_dir ***")
              }
            }, error = function(e) {
              message(paste("Warning: Could not save post-normalization matrix:", e$message))
            })
            
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
              message("*** STEP 3a: Running automatic RUV optimization ***")
              
              # Validate automatic parameters
              percentage_min <- if(is.null(input$auto_percentage_min) || is.na(input$auto_percentage_min)) 1 else input$auto_percentage_min
              percentage_max <- if(is.null(input$auto_percentage_max) || is.na(input$auto_percentage_max)) 20 else input$auto_percentage_max
              
              if (percentage_min >= percentage_max) {
                stop("Minimum percentage must be less than maximum percentage")
              }
              
              percentage_range <- seq(percentage_min, percentage_max, by = 1)
              
              # Run findBestNegCtrlPercentage optimization
              shiny::withProgress(message = "Optimizing RUV parameters...", detail = "Testing percentage range...", value = 0.5, {
                
                # 🔧 DEBUG: Log all parameters before calling optimization
                message("=== AUTOMATIC RUV OPTIMIZATION DEBUG ===")
                message(sprintf("*** DEBUG: percentage_range = %s ***", paste(percentage_range, collapse = ", ")))
                message(sprintf("*** DEBUG: separation_metric = %s ***", input$separation_metric))
                message(sprintf("*** DEBUG: k_penalty_weight = %.2f ***", input$k_penalty_weight))
                message(sprintf("*** DEBUG: adaptive_k_penalty = %s ***", input$adaptive_k_penalty))
                message(sprintf("*** DEBUG: ruv_grouping_variable = %s ***", getRuvGroupingVariable()))
                message(sprintf("*** DEBUG: normalized_s4 class = %s ***", class(normalized_s4)))
                message(sprintf("*** DEBUG: normalized_s4 dims = %d x %d ***", nrow(normalized_s4@protein_quant_table), ncol(normalized_s4@protein_quant_table)))
                
                optimization_result <- findBestNegCtrlPercentage(
                  normalised_protein_matrix_obj = normalized_s4,
                  percentage_range = percentage_range,
                  separation_metric = input$separation_metric,
                  k_penalty_weight = input$k_penalty_weight,
                  adaptive_k_penalty = input$adaptive_k_penalty,
                  ruv_grouping_variable = getRuvGroupingVariable(),
                  verbose = TRUE
                )
                
                # 🔧 DEBUG: Log optimization results
                message("*** DEBUG: Optimization completed ***")
                message(sprintf("*** DEBUG: optimization_result class = %s ***", class(optimization_result)))
                message(sprintf("*** DEBUG: best_percentage = %.1f ***", optimization_result$best_percentage))
                message(sprintf("*** DEBUG: best_k = %d ***", optimization_result$best_k))
                message(sprintf("*** DEBUG: best_composite_score = %.6f ***", optimization_result$best_composite_score))
                message(sprintf("*** DEBUG: best_separation_score = %.6f ***", optimization_result$best_separation_score))
                if (!is.null(optimization_result$optimization_results)) {
                  message("*** DEBUG: First few optimization results: ***")
                  print(head(optimization_result$optimization_results, 3))
                } else {
                  message("*** DEBUG: optimization_results is NULL! ***")
                }
                
                # 🔧 DETAILED DEBUG: Check optimization_results structure
                message("*** DEBUG: Checking optimization_results structure ***")
                message(sprintf("*** DEBUG: optimization_results class = %s ***", class(optimization_result$optimization_results)))
                if (is.data.frame(optimization_result$optimization_results)) {
                  message(sprintf("*** DEBUG: optimization_results dims = %d x %d ***", 
                                 nrow(optimization_result$optimization_results), 
                                 ncol(optimization_result$optimization_results)))
                  message(sprintf("*** DEBUG: optimization_results columns = %s ***", 
                                 paste(colnames(optimization_result$optimization_results), collapse = ", ")))
                  
                  # Show full results table for debugging
                  message("*** DEBUG: FULL optimization_results table: ***")
                  print(optimization_result$optimization_results)
                } else {
                  message("*** DEBUG: optimization_results is not a data.frame ***")
                }
                
                # 🔧 DEBUG: Test if verbose is working by calling with message override
                message("*** DEBUG: Testing verbose logging override ***")
                message("*** DEBUG: This should help us see if the function is running at all ***")
                
                            # Store optimization results for display
            norm_data$ruv_optimization_result <- optimization_result
            
            # ✅ NEW: Store RUV optimization results in workflow_data for session summary
            workflow_data$ruv_optimization_result <- optimization_result
            cat("*** STEP 3a: Stored RUV optimization results in workflow_data ***\n")
            
            # ✅ NEW: Save RUV optimization results to file for persistence
            if (!is.null(experiment_paths$source_dir)) {
              tryCatch({
                ruv_file <- file.path(experiment_paths$source_dir, "ruv_optimization_results.RDS")
                saveRDS(optimization_result, ruv_file)
                cat(sprintf("*** STEP 3a: Saved RUV optimization results to: %s ***\n", ruv_file))
              }, error = function(e) {
                cat(sprintf("*** STEP 3a: Warning - could not save RUV results file: %s ***\n", e$message))
              })
            }
            
            # Extract optimized parameters
            percentage_as_neg_ctrl <- optimization_result$best_percentage
            ruv_k <- optimization_result$best_k
            control_genes_index <- optimization_result$best_control_genes_index
            
            message(sprintf("*** STEP 3a: Automatic optimization completed - Best %%: %.1f, Best k: %d ***", 
                           percentage_as_neg_ctrl, ruv_k))
              })
              
              # ✅ AUDIT TRAIL: Log automatic RUV parameters to global config_list
              tryCatch({
                if (exists("config_list", envir = .GlobalEnv)) {
                  config_list <- get("config_list", envir = .GlobalEnv)
                  config_list <- updateRuvParameters(config_list, ruv_k, control_genes_index, percentage_as_neg_ctrl)
                  assign("config_list", config_list, envir = .GlobalEnv)
                  message("*** AUDIT: Logged automatic RUV parameters to config_list ***")
                } else {
                  message("*** WARNING: config_list not found in global environment - audit trail not updated ***")
                }
              }, error = function(e) {
                message(paste("*** WARNING: Could not update config_list audit trail:", e$message, "***"))
              })
              
            } else {
              # Manual mode - use existing logic
              message("*** STEP 3a: Using manual RUV parameters ***")
              percentage_as_neg_ctrl <- input$ruv_percentage
              ruv_k <- if(is.null(input$ruv_k) || is.na(input$ruv_k)) 3 else input$ruv_k
              
              # Get negative control proteins for manual parameters
              control_genes_index <- getNegCtrlProtAnova(
                normalized_s4,
                ruv_grouping_variable = getRuvGroupingVariable(),
                percentage_as_neg_ctrl = percentage_as_neg_ctrl,
                ruv_qval_cutoff = 0.05,
                ruv_fdr_method = "BH"
              )
              
              # ✅ AUDIT TRAIL: Log manual RUV parameters to global config_list
              tryCatch({
                if (exists("config_list", envir = .GlobalEnv)) {
                  config_list <- get("config_list", envir = .GlobalEnv)
                  config_list <- updateRuvParameters(config_list, ruv_k, control_genes_index, percentage_as_neg_ctrl)
                  assign("config_list", config_list, envir = .GlobalEnv)
                  message("*** AUDIT: Logged manual RUV parameters to config_list ***")
                } else {
                  message("*** WARNING: config_list not found in global environment - audit trail not updated ***")
                }
              }, error = function(e) {
                message(paste("*** WARNING: Could not update config_list audit trail:", e$message, "***"))
              })
              
              # For manual mode, generate canonical correlation plot for display
              tryCatch({
                cancor_plot <- ruvCancor(
                  normalized_s4,
                  ctrl = control_genes_index,
                  num_components_to_impute = 2,
                  ruv_grouping_variable = getRuvGroupingVariable()
                )
                
                              # Store manual results in same format as automatic for consistent display
              manual_ruv_result <- list(
                best_percentage = percentage_as_neg_ctrl,
                best_k = ruv_k,
                best_control_genes_index = control_genes_index,
                best_cancor_plot = cancor_plot,
                optimization_results = data.frame(
                  percentage = percentage_as_neg_ctrl,
                  separation_score = NA,
                  best_k = ruv_k,
                  composite_score = NA,
                  num_controls = sum(control_genes_index, na.rm = TRUE),
                  valid_plot = TRUE
                ),
                separation_metric_used = "manual",
                k_penalty_weight = NA,
                adaptive_k_penalty_used = FALSE
              )
              
              norm_data$ruv_optimization_result <- manual_ruv_result
              
              # ✅ NEW: Store manual RUV results in workflow_data for session summary
              workflow_data$ruv_optimization_result <- manual_ruv_result
              cat("*** STEP 3a: Stored manual RUV results in workflow_data ***\n")
              
              # ✅ NEW: Save manual RUV results to file for persistence
              if (!is.null(experiment_paths$source_dir)) {
                tryCatch({
                  ruv_file <- file.path(experiment_paths$source_dir, "ruv_optimization_results.RDS")
                  saveRDS(manual_ruv_result, ruv_file)
                  cat(sprintf("*** STEP 3a: Saved manual RUV results to: %s ***\n", ruv_file))
                }, error = function(e) {
                  cat(sprintf("*** STEP 3a: Warning - could not save manual RUV results file: %s ***\n", e$message))
                })
              }
              }, error = function(e) {
                message(paste("Warning: Could not generate canonical correlation plot for manual mode:", e$message))
              })
            }
            
            norm_data$control_genes_index <- control_genes_index
            norm_data$best_k <- ruv_k
            message("*** STEP 3: RUV parameter determination completed ***")
            
          }, error = function(e) {
            stop(paste("Step 3 (RUV parameter determination) error:", e$message))
          })
          
          # Apply RUV-III correction
          shiny::incProgress(0.2, detail = "Applying RUV-III batch correction...")
          message("*** STEP 4: Applying RUV-III correction ***")
          
          tryCatch({
            # Get RUV parameters from norm_data to ensure accessibility
            ruv_k <- norm_data$best_k
            control_genes_index <- norm_data$control_genes_index
            
            ruv_corrected_s4 <- ruvIII_C_Varying(
              normalized_s4,
              ruv_grouping_variable = getRuvGroupingVariable(),
              ruv_number_k = ruv_k,
              ctrl = control_genes_index
            )
            
            # IMMEDIATELY store the RUV result for Step 6 access
            norm_data$ruv_normalized_obj <- ruv_corrected_s4
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
            # Wrap in tryCatch to prevent workflow failure
            message("*** STEP 5: About to call updateProteinFiltering ***")
            message(sprintf("*** STEP 5: omic_type = %s, experiment_label = %s ***", 
                          ifelse(is.null(omic_type), "NULL", omic_type),
                          ifelse(is.null(experiment_label), "NULL", experiment_label)))
            
            # SKIP updateProteinFiltering for now to test if this is the problem
            if (FALSE) {  # TEMPORARILY DISABLED
              tryCatch({
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
                message("*** STEP 5: Protein filtering tracking updated successfully ***")
                
              }, error = function(e) {
                message(paste("*** WARNING: updateProteinFiltering failed:", e$message, "***"))
                message("*** STEP 5: Continuing without filtering plot update ***")
                norm_data$post_norm_filtering_plot <- NULL
              })
            } else {
              message("*** STEP 5: SKIPPING updateProteinFiltering to test workflow ***")
              norm_data$post_norm_filtering_plot <- NULL
            }
            
            # Get RUV parameters BEFORE using them
            best_percentage <- if (!is.null(norm_data$ruv_optimization_result)) {
              norm_data$ruv_optimization_result$best_percentage
            } else {
              input$ruv_percentage  # Fallback for manual mode
            }
            
            # Generate filtering summary text
            norm_data$filtering_summary_text <- sprintf(
              "Filtering Summary through Normalization & RUV:\n\n• Pre-normalization proteins: [Previous step]\n• Post-RUV filtering: %d proteins\n• Proteins removed by RUV: [Calculated from difference]\n\nRUV Parameters:\n• Normalization method: %s\n• RUV mode: %s\n• RUV k value: %d\n• Negative control %%: %.1f",
              ruvfilt_protein_count,
              input$norm_method,
              input$ruv_mode,
              norm_data$best_k,
              best_percentage
            )
            
            norm_data$ruv_normalized_obj <- ruv_corrected_s4_clean
            message("*** STEP 5: Missing values cleanup completed ***")
            
            # Save post-normalization state to R6 state manager (includes RUV + missing value filtering)
            message("*** STEP 5: Saving state to R6 state manager ***")
            
            tryCatch({
            workflow_data$state_manager$saveState(
              state_name = "ruv_corrected",
              s4_data_object = ruv_corrected_s4_clean,
              config_object = list(
                norm_method = input$norm_method,
                ruv_mode = input$ruv_mode,
                  ruv_k = norm_data$best_k,
                  percentage_as_neg_ctrl = best_percentage
              ),
              description = "Post-normalization complete: RUV-III correction and missing value cleanup completed"
            )
              message("*** STEP 5: State saved successfully ***")
            }, error = function(e) {
              message(paste("*** WARNING: Could not save state to R6 manager:", e$message, "***"))
              message("*** STEP 5: Continuing without state save (Step 6 will still proceed) ***")
            })
            
          }, error = function(e) {
            message(paste("*** ERROR in Step 5:", e$message, "***"))
            message("*** STEP 5: Continuing to Step 6 despite error ***")
          })
          
          # STEP 6 - ALWAYS GENERATE COMPOSITE QC (MOVED OUTSIDE STEP 5 TRY-CATCH)
          # Generate RUV-corrected QC plots
          shiny::incProgress(0.2, detail = "Generating RUV-corrected QC plots...")
          message("*** STEP 6: STARTING RUV-corrected QC plot generation ***")
          
          # Check if we have the cleaned object from Step 5
          if (exists("ruv_corrected_s4_clean") && !is.null(ruv_corrected_s4_clean)) {
            message(sprintf("*** STEP 6: ruv_corrected_s4_clean is available ***"))
          } else if (!is.null(norm_data$ruv_normalized_obj)) {
            # Fallback to the object stored in norm_data if Step 5 failed
            ruv_corrected_s4_clean <- norm_data$ruv_normalized_obj
            message("*** STEP 6: Using fallback ruv_normalized_obj from norm_data ***")
          } else {
            message("*** STEP 6: WARNING - No RUV corrected object available ***")
          }
          
          message(sprintf("*** STEP 6: experiment_paths$protein_qc_dir: %s ***", experiment_paths$protein_qc_dir))
          
          # STEP 6A: Try to generate RUV QC plots (but don't fail if this breaks)
          tryCatch({
            generateRuvCorrectedQc(ruv_corrected_s4_clean)
            message("*** STEP 6A: RUV-corrected QC plots completed ***")
          }, error = function(e) {
            message(paste("Warning: generateRuvCorrectedQc failed:", e$message))
            message("*** STEP 6A: Continuing without RUV QC plots ***")
          })
          
          # STEP 6B: ALWAYS generate composite QC figure (critical for copyToResultsSummary)
          message(sprintf("*** STEP 6B: norm_data$qc_plots$post_filtering$pca is NULL: %s ***", is.null(norm_data$qc_plots$post_filtering$pca)))
          message(sprintf("*** STEP 6B: norm_data$qc_plots$post_normalization$pca is NULL: %s ***", is.null(norm_data$qc_plots$post_normalization$pca)))
          message(sprintf("*** STEP 6B: norm_data$qc_plots$ruv_corrected$pca is NULL: %s ***", is.null(norm_data$qc_plots$ruv_corrected$pca)))
          
          # ENSURE S4 object has all plots populated before calling createGridQC
          if (!is.null(norm_data$QC_composite_figure)) {
            message("*** STEP 6B: Ensuring S4 object has all plots populated ***")
            
            # Populate pre-norm plots if missing
            if (is.null(norm_data$QC_composite_figure@pca_plots$pca_plot_before_cyclic_loess_group) && 
                !is.null(norm_data$qc_plots$post_filtering$pca)) {
              norm_data$QC_composite_figure@pca_plots$pca_plot_before_cyclic_loess_group <- norm_data$qc_plots$post_filtering$pca
              norm_data$QC_composite_figure@rle_plots$rle_plot_before_cyclic_loess_group <- norm_data$qc_plots$post_filtering$rle
              norm_data$QC_composite_figure@density_plots$density_plot_before_cyclic_loess_group <- norm_data$qc_plots$post_filtering$density
              norm_data$QC_composite_figure@pearson_plots$pearson_correlation_pair_before_cyclic_loess <- norm_data$qc_plots$post_filtering$correlation
            }
            
            # Populate post-norm plots if missing
            if (is.null(norm_data$QC_composite_figure@pca_plots$pca_plot_before_ruvIIIc_group) && 
                !is.null(norm_data$qc_plots$post_normalization$pca)) {
              norm_data$QC_composite_figure@pca_plots$pca_plot_before_ruvIIIc_group <- norm_data$qc_plots$post_normalization$pca
              norm_data$QC_composite_figure@rle_plots$rle_plot_before_ruvIIIc_group <- norm_data$qc_plots$post_normalization$rle
              norm_data$QC_composite_figure@density_plots$density_plot_before_ruvIIIc_group <- norm_data$qc_plots$post_normalization$density
              norm_data$QC_composite_figure@pearson_plots$pearson_correlation_pair_before_ruvIIIc <- norm_data$qc_plots$post_normalization$correlation
            }
            
            # Populate RUV plots if missing
            if (is.null(norm_data$QC_composite_figure@pca_plots$pca_plot_after_ruvIIIc_group) && 
                !is.null(norm_data$qc_plots$ruv_corrected$pca)) {
              norm_data$QC_composite_figure@pca_plots$pca_plot_after_ruvIIIc_group <- norm_data$qc_plots$ruv_corrected$pca
              norm_data$QC_composite_figure@rle_plots$rle_plot_after_ruvIIIc_group <- norm_data$qc_plots$ruv_corrected$rle
              norm_data$QC_composite_figure@density_plots$density_plot_after_ruvIIIc_group <- norm_data$qc_plots$ruv_corrected$density
              norm_data$QC_composite_figure@pearson_plots$pearson_correlation_pair_after_ruvIIIc_group <- norm_data$qc_plots$ruv_corrected$correlation
            }
            
            # Add cancor plots for the three stages
            if (!is.null(norm_data$ruv_optimization_result) && !is.null(norm_data$ruv_optimization_result$best_cancor_plot)) {
              # For simplicity, use the same cancor plot for all three stages (or generate different ones if available)
              norm_data$QC_composite_figure@cancor_plots$cancor_plot_before_cyclic_loess <- NULL  # No cancor before normalization
              norm_data$QC_composite_figure@cancor_plots$cancor_plot_before_ruvIIIc <- NULL      # No cancor before RUV
              norm_data$QC_composite_figure@cancor_plots$cancor_plot_after_ruvIIIc <- norm_data$ruv_optimization_result$best_cancor_plot
            }
          }
          
          tryCatch({
            # Use the S4 QC_composite_figure object that we've been populating
            message("*** STEP 6B: Using S4 QC_composite_figure object ***")
            message(sprintf("*** STEP 6B: QC_composite_figure class: %s ***", class(norm_data$QC_composite_figure)))
            
            # Generate composite QC figure and save to protein_qc_dir
            message("*** STEP 6B: About to call createGridQC ***")
            
            pca_ruv_rle_correlation_merged <- createGridQC(
              norm_data$QC_composite_figure,
              pca_titles = c("a) Pre-normalization", "b) Normalized data", "c) RUV-corrected data"),
              density_titles = c("d)", "e)", "f)"),
              rle_titles = c("g)", "h)", "i)"),
              pearson_titles = c("j)", "k)", "l)"),
              cancor_titles = c("", "", "m)"),  # Empty for first two columns, "m)" for RUV column
              ncol = 3,
              save_path = experiment_paths$protein_qc_dir,
              file_name = "composite_QC_figure"
            )
            
            message("*** STEP 6B: Composite QC figure saved to protein_qc_dir ***")
            
          }, error = function(e) {
            message(paste("Warning: Could not generate composite QC figure:", e$message))
            message("*** STEP 6B: This may cause issues with session summary file copying ***")
          })
          
          # Mark normalization as complete regardless of QC plot generation
          norm_data$ruv_complete <- TRUE
          message("*** STEP 6: RUV-corrected workflow completed ***")
          
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
            correlation_group = getRuvGroupingVariable()
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
            aesthetics <- getPlotAesthetics()
            norm_data$final_qc_plot <- plotPca(
              final_s4_for_de,
              grouping_variable = aesthetics$color_var,
              label_column = "",
              shape_variable = aesthetics$shape_var,
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
    
    # NEW: Render RUV QC outputs
    output$ruv_canonical_correlation_plot <- shiny::renderPlot({
      if (!is.null(norm_data$ruv_optimization_result) && !is.null(norm_data$ruv_optimization_result$best_cancor_plot)) {
        norm_data$ruv_optimization_result$best_cancor_plot
      } else {
        plot.new()
        text(0.5, 0.5, "Run normalization to generate RUV canonical correlation plot", cex = 1.2)
      }
    })
    
    output$ruv_optimization_summary <- shiny::renderText({
      if (!is.null(norm_data$ruv_optimization_result)) {
        result <- norm_data$ruv_optimization_result
        
        if (input$ruv_mode == "automatic") {
          sprintf(
            "Automatic RUV Optimization Results:\n\n• Best percentage: %.1f%%\n• Best k value: %d\n• Separation score: %.4f\n• Composite score: %.4f\n• Control genes: %d\n• RUV grouping variable: %s\n• Separation metric: %s\n• K penalty weight: %.1f\n• Adaptive penalty: %s\n• Sample size: %s",
            result$best_percentage,
            result$best_k,
            ifelse(is.null(result$best_separation_score), 0, result$best_separation_score),
            ifelse(is.null(result$best_composite_score), 0, result$best_composite_score),
            sum(result$best_control_genes_index, na.rm = TRUE),
            getRuvGroupingVariable(),
            result$separation_metric_used,
            ifelse(is.null(result$k_penalty_weight), 0, result$k_penalty_weight),
            ifelse(is.null(result$adaptive_k_penalty_used), FALSE, result$adaptive_k_penalty_used),
            ifelse(is.null(result$sample_size), "N/A", result$sample_size)
          )
        } else {
          sprintf(
            "Manual RUV Parameters:\n\n• Selected percentage: %.1f%%\n• Selected k value: %d\n• Control genes: %d\n• RUV grouping variable: %s\n• Mode: Manual selection",
            result$best_percentage,
            result$best_k,
            sum(result$best_control_genes_index, na.rm = TRUE),
            getRuvGroupingVariable()
          )
        }
      } else {
        "Run normalization to see RUV optimization results"
      }
    })
    
    output$ruv_optimization_table <- DT::renderDataTable({
      if (!is.null(norm_data$ruv_optimization_result) && !is.null(norm_data$ruv_optimization_result$optimization_results)) {
        
        # Format the optimization results table
        results_df <- norm_data$ruv_optimization_result$optimization_results
        
        # Round numeric columns for display
        if ("separation_score" %in% colnames(results_df)) {
          results_df$separation_score <- round(results_df$separation_score, 4)
        }
        if ("composite_score" %in% colnames(results_df)) {
          results_df$composite_score <- round(results_df$composite_score, 4)
        }
        
        # Highlight the best result
        best_percentage <- norm_data$ruv_optimization_result$best_percentage
        
        DT::datatable(
          results_df,
          options = list(
            pageLength = 10,
            scrollY = "300px",
            scrollCollapse = TRUE,
            dom = 't'  # Only show table, no search/pagination controls
          ),
          rownames = FALSE
        ) |>
          DT::formatStyle(
            "percentage",
            target = "row",
            backgroundColor = DT::styleEqual(best_percentage, "#e6f3ff")
          )
      } else {
        DT::datatable(data.frame(Message = "Run normalization to see optimization results"))
      }
    })
    
    # Reset normalization button logic
    observeEvent(input$reset_normalization, {
      message("Resetting normalization...")
      
      tryCatch({
        # ✅ CRITICAL FIX: Revert R6 state manager to pre-normalization state
        if (!is.null(workflow_data$state_manager)) {
          # Use smart revert pattern like QC tabs - find the actual previous state
          history <- workflow_data$state_manager$getHistory()
          # Define possible pre-normalization states in reverse chronological order
          pre_norm_states <- c("protein_replicate_filtered", "imputed", "replicate_filtered", "sample_filtered", "protein_peptide_filtered", "intensity_filtered", "precursor_rollup", "qvalue_filtered", "raw_data_s4")
          # Find the most recent state that actually exists
          prev_state <- intersect(rev(history), pre_norm_states)[1]
          
          if (!is.null(prev_state)) {
            reverted_s4 <- workflow_data$state_manager$revertToState(prev_state)
            message(sprintf("*** RESET: Reverted R6 state manager to '%s' (actual previous state) ***", prev_state))
          } else {
            message("*** WARNING: No valid pre-normalization state found in history ***")
          }
        } else {
          message("*** WARNING: workflow_data$state_manager is NULL - cannot revert state ***")
        }
        
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
        norm_data$ruv_optimization_result <- NULL  # Reset RUV optimization results
        
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
        output$ruv_optimization_summary <- shiny::renderText("Run normalization to see RUV optimization results")
        
        # Create dynamic notification message based on actual revert state
        notification_msg <- if (!is.null(prev_state)) {
          sprintf("Normalization has been reset to pre-normalization state (%s)", prev_state)
        } else {
          "Normalization has been reset (no valid previous state found)"
        }
        
        shiny::showNotification(
          notification_msg,
          type = "warning",
          duration = 5
        )
        
        message("*** RESET: Normalization reset completed successfully ***")
        
      }, error = function(e) {
        message(paste("*** ERROR in normalization reset:", e$message, "***"))
        shiny::showNotification(
          paste("Error resetting normalization:", e$message),
          type = "error",
          duration = 10
        )
      })
    })
    
    # Export filtered session button logic
    observeEvent(input$export_filtered_session, {
      message("=== EXPORT FILTERED SESSION BUTTON CLICKED ===")
      
      # Check if correlation filtering is complete
      if (!norm_data$correlation_filtering_complete || is.null(norm_data$correlation_filtered_obj)) {
        shiny::showNotification(
          "Please complete correlation filtering before exporting session data.",
          type = "warning",
          duration = 5
        )
        return()
      }
      
      tryCatch({
        # Get the source directory for saving
        source_dir <- experiment_paths$source_dir
        if (is.null(source_dir) || !dir.exists(source_dir)) {
          stop("Could not find the source directory to save session data.")
        }
        
        shiny::withProgress(message = "Exporting filtered session data...", value = 0, {
          
          # Step 1: Gather all essential data for DE analysis
          shiny::incProgress(0.2, detail = "Gathering data...")
          
                     # Get the current R6 state manager state and complete state info
           current_state_name <- workflow_data$state_manager$current_state
           current_s4_object <- workflow_data$state_manager$getState(current_state_name)
           
           # Export the complete R6 states structure for proper restoration
           r6_complete_states <- workflow_data$state_manager$states
           r6_state_history <- workflow_data$state_manager$state_history
           
           session_data <- list(
             # Complete R6 State Manager Information for proper restoration
             r6_current_state_name = current_state_name,
             r6_complete_states = r6_complete_states,
             r6_state_history = r6_state_history,
             
             # Current S4 object from R6 state manager (for quick access)
             current_s4_object = current_s4_object,
             
             # Backup S4 object from normalization module (should be same)
             correlation_filtered_s4 = norm_data$correlation_filtered_obj,
             
             # Contrasts table with all required columns
             contrasts_tbl = workflow_data$contrasts_tbl,
             
             # Design matrix
             design_matrix = workflow_data$design_matrix,
             
             # Configuration
             config_list = workflow_data$config_list,
             
             # Additional metadata for verification
             export_timestamp = Sys.time(),
            normalization_method = input$norm_method,
            ruv_mode = input$ruv_mode,
            ruv_k = norm_data$best_k,
            correlation_threshold = norm_data$correlation_threshold,
            
            # NEW: Store workflow type for report template selection
            workflow_type = if (!is.null(workflow_data$config_list) && 
                                !is.null(workflow_data$config_list$globalParameters) &&
                                !is.null(workflow_data$config_list$globalParameters$workflow_type)) {
              workflow_data$config_list$globalParameters$workflow_type
            } else if (!is.null(current_s4_object@args$globalParameters$workflow_type)) {
              current_s4_object@args$globalParameters$workflow_type
            } else {
              "DIA"  # Default fallback for backward compatibility
            },
            
            # NEW: FASTA metadata for complete audit trail
            fasta_metadata = workflow_data$fasta_metadata,
            
            # NEW: Accession cleanup results
            accession_cleanup_results = workflow_data$accession_cleanup_results,
            
            # NEW: Complete RUV optimization results (not just ruv_k)
            ruv_optimization_result = workflow_data$ruv_optimization_result,
            
            # NEW: QC parameters from all QC steps
            qc_params = workflow_data$qc_params,
            
            # Data dimensions for verification
             final_protein_count = length(unique(current_s4_object@protein_quant_table$Protein.Ids)),
             final_sample_count = length(setdiff(colnames(current_s4_object@protein_quant_table), current_s4_object@protein_id_column))
           )
          
          message("*** EXPORT: Gathered session data successfully ***")
          message(sprintf("*** EXPORT: Final protein count: %d ***", session_data$final_protein_count))
          message(sprintf("*** EXPORT: Final sample count: %d ***", session_data$final_sample_count))
          message(sprintf("*** EXPORT: Contrasts available: %d ***", ifelse(is.null(session_data$contrasts_tbl), 0, nrow(session_data$contrasts_tbl))))
          
          # Step 2: Save to RDS file with timestamp
          shiny::incProgress(0.4, detail = "Saving to file...")
          
          timestamp_str <- format(Sys.time(), "%Y%m%d_%H%M%S")
          session_filename <- sprintf("filtered_session_data_%s.rds", timestamp_str)
          session_filepath <- file.path(source_dir, session_filename)
          
          saveRDS(session_data, session_filepath)
          message(sprintf("*** EXPORT: Session data saved to: %s ***", session_filepath))
          
          # Step 3: Also save a "latest" version for easy access
          shiny::incProgress(0.2, detail = "Creating latest version...")
          
          latest_filename <- "filtered_session_data_latest.rds"
          latest_filepath <- file.path(source_dir, latest_filename)
          
          saveRDS(session_data, latest_filepath)
          message(sprintf("*** EXPORT: Latest version saved to: %s ***", latest_filepath))
          
          # Step 3.5: Save individual metadata files for redundancy and easy access
          shiny::incProgress(0.1, detail = "Saving metadata files...")
          
          tryCatch({
            # Save accession cleanup results
            if (!is.null(session_data$accession_cleanup_results)) {
              accession_cleanup_file <- file.path(source_dir, "accession_cleanup_results.RDS")
              saveRDS(session_data$accession_cleanup_results, accession_cleanup_file)
              message(sprintf("*** EXPORT: Saved accession_cleanup_results.RDS ***"))
            }
            
            # Save QC parameters
            if (!is.null(session_data$qc_params)) {
              qc_params_file <- file.path(source_dir, "qc_params.RDS")
              saveRDS(session_data$qc_params, qc_params_file)
              message(sprintf("*** EXPORT: Saved qc_params.RDS ***"))
            }
            
            # RUV optimization and FASTA metadata already saved in earlier steps, but verify/update
            if (!is.null(session_data$fasta_metadata)) {
              fasta_metadata_file <- file.path(source_dir, "fasta_metadata.RDS")
              if (!file.exists(fasta_metadata_file)) {
                saveRDS(session_data$fasta_metadata, fasta_metadata_file)
                message(sprintf("*** EXPORT: Saved fasta_metadata.RDS ***"))
              }
            }
            
            if (!is.null(session_data$ruv_optimization_result)) {
              ruv_file <- file.path(source_dir, "ruv_optimization_results.RDS")
              if (!file.exists(ruv_file)) {
                saveRDS(session_data$ruv_optimization_result, ruv_file)
                message(sprintf("*** EXPORT: Saved ruv_optimization_results.RDS ***"))
              }
            }
            
          }, error = function(e) {
            message(sprintf("*** WARNING: Some metadata files could not be saved: %s ***", e$message))
          })
          
          # Step 4: Create a summary file for user reference
          shiny::incProgress(0.1, detail = "Creating summary...")
          
          summary_content <- sprintf(
            "Filtered Session Data Export Summary\n=====================================\n\nExport Timestamp: %s\nSession File: %s\n\nData Summary:\n- Proteins: %d\n- Samples: %d\n- Contrasts: %d\n- Normalization: %s\n- RUV Mode: %s\n- RUV K: %d\n- Correlation Threshold: %.3f\n\nContrasts:\n%s\n\nThis data is ready for differential expression analysis.\nUse 'Load Filtered Session' in the DE tab to import.\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            session_filename,
            session_data$final_protein_count,
            session_data$final_sample_count,
            ifelse(is.null(session_data$contrasts_tbl), 0, nrow(session_data$contrasts_tbl)),
            session_data$normalization_method,
            session_data$ruv_mode,
            ifelse(is.null(session_data$ruv_k), "NA", session_data$ruv_k),
            ifelse(is.null(session_data$correlation_threshold), "NA", session_data$correlation_threshold),
            if (!is.null(session_data$contrasts_tbl)) paste(session_data$contrasts_tbl$friendly_names, collapse = "\n") else "None"
          )
          
          summary_filepath <- file.path(source_dir, "filtered_session_summary.txt")
          writeLines(summary_content, summary_filepath)
          
        })
        
        shiny::showNotification(
          sprintf("Filtered session data exported successfully!\nSaved as: %s\nSee summary file for details.", session_filename),
          type = "success",
          duration = 10
        )
        
        message("=== EXPORT FILTERED SESSION COMPLETED SUCCESSFULLY ===")
        
      }, error = function(e) {
        message(paste("*** ERROR in session export:", e$message, "***"))
        shiny::showNotification(
          paste("Error exporting session:", e$message),
          type = "error",
          duration = 10
        )
      })
    })
    
    # Return normalization data for potential use by parent module
    return(norm_data)
  })
} 