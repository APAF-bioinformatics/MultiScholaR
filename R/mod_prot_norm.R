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
#' @importFrom shiny NS tagList fluidRow column wellPanel h4 selectInput helpText hr numericInput checkboxInput actionButton icon br sliderInput h5 verbatimTextOutput plotOutput tabsetPanel tabPanel conditionalPanel
mod_prot_norm_ui <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::tagList(
    shiny::wellPanel(
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
            "Manual (advanced users)" = "manual",
            "Skip RUV (troublesome datasets)" = "skip"
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
        
        # Skip RUV mode information (shown when skip mode selected)
        shiny::conditionalPanel(
          condition = "input.ruv_mode == 'skip'",
          ns = ns,
          
          shiny::div(
            style = "background-color: #fff3cd; border: 1px solid #ffc107; padding: 10px; border-radius: 4px;",
            shiny::h5("RUV-III Will Be Skipped", style = "color: #856404; margin-top: 0;"),
            shiny::helpText(
              "Batch correction will not be applied. Use this option for datasets where RUV-III is problematic or unnecessary.",
              style = "margin-bottom: 5px;"
            ),
            shiny::helpText(
              "Alternative approach: Consider including batch as a covariate in your differential expression model if batch effects are present.",
              style = "margin-bottom: 0; font-style: italic;"
            )
          ),
          
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
                
                shiny::actionButton(
                  ns("skip_correlation_filter"),
                  "Skip Filtering & Proceed to DE",
                  class = "btn-info",
                  width = "100%",
                  icon = shiny::icon("forward")
                ),
                shiny::helpText("Use this option to bypass sample filtering and proceed directly to Differential Expression with the RUV-corrected data."),
                
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
  )
  )
}

#' @rdname normalizationAppletModule 
#' @export
mod_prot_norm_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label, selected_tab = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    
    message("=== NORMALIZATION MODULE SERVER STARTED ===")
    message(sprintf("Module ID: %s", id))
    message(sprintf("workflow_data is NULL: %s", is.null(workflow_data)))
    if (!is.null(workflow_data$state_manager)) {
      message(sprintf("Current state at module start: %s", workflow_data$state_manager$current_state))
    }
    
    # MEMORY MONITORING: Helper function to check memory usage and warn if high
    checkMemoryUsage <- function(threshold_gb = 8, context = "") {
      mem_info <- gc()
      mem_used_mb <- sum(mem_info[,2])  # Used memory in MB
      mem_used_gb <- mem_used_mb / 1024
      
      if (mem_used_gb > threshold_gb) {
        warning(sprintf("*** HIGH MEMORY WARNING [%s]: %.1f GB used (threshold: %.1f GB) ***", 
                       context, mem_used_gb, threshold_gb))
        message(sprintf("*** HIGH MEMORY WARNING [%s]: %.1f GB used ***", context, mem_used_gb))
      } else {
        message(sprintf("*** MEMORY CHECK [%s]: %.1f GB used ***", context, mem_used_gb))
      }
      
      invisible(mem_used_gb)
    }
    
    # Initialize reactive values for normalization state
    norm_data <- shiny::reactiveValues(
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
    shiny::observe({
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
    
    # Helper function to generate composite QC figure from saved images
    # TRUE TO SOURCE: Mirrors createGridQC structure with row labels and patchwork layout
    generateCompositeFromFiles <- function(plot_files, output_path, ncol = 3, row_labels = NULL, column_labels = NULL) {
      message(sprintf("   [generateCompositeFromFiles] Generating composite from %d files...", length(plot_files)))
      
      if (!requireNamespace("patchwork", quietly = TRUE)) {
        warning("patchwork package required for composite generation")
        return(NULL)
      }
      if (!requireNamespace("ggplot2", quietly = TRUE)) {
        warning("ggplot2 package required for composite generation")
        return(NULL)
      }
      if (!requireNamespace("png", quietly = TRUE)) {
        warning("png package required for composite generation")
        return(NULL)
      }
      
      # --- Helper: Create label plot (true to createGridQC source) ---
      createLabelPlot <- function(title) {
        ggplot2::ggplot() + 
          ggplot2::annotate("text", x = 0, y = 0.5, label = title, size = 5, hjust = 0) +
          ggplot2::xlim(0, 1) +
          ggplot2::theme_void() +
          ggplot2::theme(
            plot.margin = ggplot2::margin(5, 5, 5, 5)
            , panel.background = ggplot2::element_blank()
          )
      }

      # --- Helper: Create column title plot ---
      createTitlePlot <- function(title) {
        ggplot2::ggplot() + 
          ggplot2::annotate("text", x = 0.5, y = 0.5, label = title, size = 6, fontface = "bold", hjust = 0.5) +
          ggplot2::xlim(0, 1) +
          ggplot2::theme_void() +
          ggplot2::theme(
            plot.margin = ggplot2::margin(5, 5, 10, 5) # Extra bottom margin
            , panel.background = ggplot2::element_blank()
          )
      }
      
      # --- Helper: Load image as ggplot ---
      loadImageAsPlot <- function(file_path) {
        if (is.na(file_path) || !file.exists(file_path)) {
          # Return empty plot for missing files (maintains grid alignment)
          return(ggplot2::ggplot() + ggplot2::theme_void())
        }
        tryCatch({
          img <- png::readPNG(file_path)
          g <- grid::rasterGrob(img, interpolate = TRUE)
          ggplot2::ggplot() +
            ggplot2::annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
            ggplot2::theme_void()
        }, error = function(e) {
          message(sprintf("   [generateCompositeFromFiles] Could not load image: %s", file_path))
          ggplot2::ggplot() + ggplot2::theme_void()
        })
      }
      
      tryCatch({
        # Determine number of plot types (rows) based on file count
        n_files <- length(plot_files)
        n_plot_types <- n_files / ncol
        
        # Default row labels if not provided
        if (is.null(row_labels)) {
          # Generate default labels: a), b), c), ...
          all_labels <- letters[1:n_files]
          row_labels <- split(paste0(all_labels, ")"), rep(1:n_plot_types, each = ncol))
          names(row_labels) <- paste0("row", seq_len(n_plot_types))
        }
        
        # Build plot sections (labels row + images row for each plot type)
        plot_sections <- list()
        height_values <- c()

        # Add Column Titles if provided
        if (!is.null(column_labels)) {
           if (length(column_labels) != ncol) {
             warning("Number of column labels does not match ncol")
           } else {
             # Create a blank plot for the row label column (top-left corner)
             blank_plot <- ggplot2::ggplot() + ggplot2::theme_void()
             
             # Create title plots
             title_plots <- lapply(column_labels, createTitlePlot)
             
             # Prepend the blank plot to align with the grid (labels column + data columns)
             # Wait, the structure is:
             # Row 1: [Label a] [Plot 1] [Label b] [Plot 2] ... 
             # No, structure is:
             # Row Labels: [a] [b] [c]
             # Plots:      [P1] [P2] [P3]
             
             # So for column titles, we just need a row of titles matching the columns.
             # But the labels (a, b, c) are also in columns.
             # So it should be: [Title 1] [Title 2] [Title 3]
             
             plot_sections <- append(plot_sections, list(
               patchwork::wrap_plots(title_plots, ncol = ncol)
             ))
             height_values <- c(height_values, 0.2) # Height for title row
             message("   [generateCompositeFromFiles] Added column titles")
           }
        }
        
        row_names <- names(row_labels)
        
        for (i in seq_along(row_names)) {
          row_name <- row_names[i]
          labels <- row_labels[[row_name]]
          
          # Get file indices for this row
          start_idx <- (i - 1) * ncol + 1
          end_idx <- min(i * ncol, n_files)
          row_files <- plot_files[start_idx:end_idx]
          
          # Check if any non-NA files exist for this row
          has_files <- any(!is.na(row_files) & sapply(row_files, function(f) !is.na(f) && file.exists(f)))
          
          if (has_files || row_name == "cancor") {
            # Create label plots for this row
            label_plots <- lapply(labels, createLabelPlot)
            
            # Load image plots for this row
            image_plots <- lapply(row_files, loadImageAsPlot)
            
            # Add label row and image row to sections
            plot_sections <- append(plot_sections, list(
              patchwork::wrap_plots(label_plots, ncol = ncol)
              , patchwork::wrap_plots(image_plots, ncol = ncol)
            ))
            height_values <- c(height_values, 0.1, 1)
            
            message(sprintf("   [generateCompositeFromFiles] Added row: %s", row_name))
          } else {
            message(sprintf("   [generateCompositeFromFiles] Skipping empty row: %s", row_name))
          }
        }
        
        if (length(plot_sections) == 0) {
          warning("No valid plot sections to combine")
          return(NULL)
        }
        
        # Combine all sections vertically (true to createGridQC source)
        message("   [generateCompositeFromFiles] Combining plot sections...")
        combined_plot <- patchwork::wrap_plots(plot_sections, ncol = 1) +
          patchwork::plot_layout(heights = height_values)
        
        # Calculate dimensions based on content
        plot_width <- 4 + (ncol * 3)  # Base width + 3 units per column
        plot_height <- 4 + (length(height_values) * 2)  # Base height + 2 units per row
        
        # Save composite
        ggplot2::ggsave(
          output_path
          , combined_plot
          , width = plot_width
          , height = plot_height
          , dpi = 150
          , limitsize = FALSE
        )
        
        # Clear memory
        rm(plot_sections, combined_plot)
        gc()
        
        message(sprintf("   [generateCompositeFromFiles] Composite saved to %s", output_path))
        return(output_path)
        
      }, error = function(e) {
        message(paste("   [generateCompositeFromFiles] Error:", e$message))
        return(NULL)
      })
    }

    # Helper function to generate pre-normalization QC plots
    # MEMORY OPTIMIZED: Save to disk immediately, clear from memory
    generatePreNormalizationQc <- function() {
      message("=== GENERATING PRE-NORMALIZATION QC PLOTS ===")
      message("Step 1: Getting current S4 object...")
      
      # Get QC directory for saving plots
      qc_dir <- experiment_paths$protein_qc_dir
      
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
      
      # Initialize plot paths storage if needed
      if (is.null(norm_data$qc_plot_paths)) {
        norm_data$qc_plot_paths <- list(
          post_filtering = list(),
          post_normalization = list(),
          ruv_corrected = list()
        )
      }
      
      # Generate QC plots - SAVE TO DISK IMMEDIATELY, clear memory
      
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
      
      if (!is.null(qc_dir)) {
        pca_path <- file.path(qc_dir, "pre_norm_pca.png")
        ggplot2::ggsave(pca_path, pca_plot, width = 8, height = 6, dpi = 150)
        norm_data$qc_plot_paths$post_filtering$pca <- pca_path
      }
      
      # CLEAR MEMORY
      rm(pca_plot)
      gc()
      
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
      
      if (!is.null(qc_dir)) {
        rle_path <- file.path(qc_dir, "pre_norm_rle.png")
        ggplot2::ggsave(rle_path, rle_plot, width = 10, height = 6, dpi = 150)
        norm_data$qc_plot_paths$post_filtering$rle <- rle_path
      }
      
      # CLEAR MEMORY
      rm(rle_plot)
      gc()
      
      rle_elapsed <- as.numeric(difftime(Sys.time(), rle_start_time, units = "secs"))
      message(sprintf("*** PRE-NORM QC: RLE plot completed in %.1f seconds ***", rle_elapsed))
      
      # Density plot (Step 3/4)
      message("*** PRE-NORM QC: Generating Density plot (Step 3/4) ***")
      # Need PCA plot again for density plot input? 
      # plotPcaBox takes the PCA plot object. 
      # Since we deleted it, we need to regenerate it or redesign plotPcaBox to take data.
      # plotPcaBox in func_general_plotting.R takes a ggplot object.
      # We have to regenerate PCA plot briefly.
      
      density_start_time <- Sys.time()
      
      # REGENERATE minimal PCA for density plot (unfortunate but saves long-term memory)
      pca_plot_temp <- plotPca(
        current_s4,
        grouping_variable = aesthetics$color_var,
        label_column = "",
        shape_variable = aesthetics$shape_var, 
        title = "",
        font_size = 8
      )
      
      density_plot <- plotPcaBox(
        pca_plot_temp,
        grouping_variable = aesthetics$color_var,
        show_legend = TRUE
      )
      
      if (!is.null(qc_dir)) {
        density_path <- file.path(qc_dir, "pre_norm_density.png")
        ggplot2::ggsave(density_path, density_plot, width = 8, height = 6, dpi = 150)
        norm_data$qc_plot_paths$post_filtering$density <- density_path
      }
      
      # CLEAR MEMORY
      rm(pca_plot_temp, density_plot)
      gc()
      
      density_elapsed <- as.numeric(difftime(Sys.time(), density_start_time, units = "secs"))
      message(sprintf("*** PRE-NORM QC: Density plot completed in %.1f seconds ***", density_elapsed))
      
      # Correlation plot (Step 4/4)
      message("*** PRE-NORM QC: Generating Pearson correlation plot (Step 4/4) ***")
      num_samples <- length(setdiff(colnames(current_s4@protein_quant_table), current_s4@protein_id_column))
      estimated_pairs <- choose(num_samples, 2)
      message(sprintf("*** PRE-NORM QC: Sample count = %d, Expected pairs ≈ %d ***", num_samples, estimated_pairs))
      
      pearson_start_time <- Sys.time()
      
      # Check if we're in a Shiny context and wrap with progress indicator
      if (shiny::isRunning()) {
        correlation_plot <- shiny::withProgress(
          message = "Generating Pearson correlation plot...", 
          detail = "Calculating pairwise correlations...", 
          value = 0.5, {
            result <- plotPearson(
              current_s4,
              correlation_group = aesthetics$color_var,
              exclude_pool_samples = TRUE
            )
            shiny::incProgress(0.5)
            result
          }
        )
      } else {
        correlation_plot <- plotPearson(
          current_s4,
          correlation_group = aesthetics$color_var,
          exclude_pool_samples = TRUE
        )
      }
      
      if (!is.null(qc_dir)) {
        corr_path <- file.path(qc_dir, "pre_norm_correlation.png")
        ggplot2::ggsave(corr_path, correlation_plot, width = 10, height = 8, dpi = 150)
        norm_data$qc_plot_paths$post_filtering$correlation <- corr_path
      }
      
      # CLEAR MEMORY
      rm(correlation_plot)
      gc()
      
      pearson_elapsed <- as.numeric(difftime(Sys.time(), pearson_start_time, units = "secs"))
      message(sprintf("*** PRE-NORM QC: Pearson correlation completed in %.1f seconds ***", pearson_elapsed))
      
      # NO S4 POPULATION - We are dropping QC_composite_figure logic
      
      # MEMORY CLEANUP
      message("*** PRE-NORM QC: Running garbage collection ***")
      gc()
      
      message("Pre-normalization QC plots generated successfully")
    }
    
    # Helper function to generate post-normalization QC plots
    # MEMORY OPTIMIZED: Save to disk immediately, clear from memory
    generatePostNormalizationQc <- function(normalized_s4) {
      message("=== GENERATING POST-NORMALIZATION QC PLOTS ===")
      
      # Get QC directory for saving plots
      qc_dir <- experiment_paths$protein_qc_dir
      
      aesthetics <- getPlotAesthetics()
      
      # Initialize plot paths storage if needed
      if (is.null(norm_data$qc_plot_paths)) {
        norm_data$qc_plot_paths <- list(
          post_filtering = list(),
          post_normalization = list(),
          ruv_corrected = list()
        )
      }
      
      # PCA plot
      message("*** POST-NORM QC: Generating PCA plot ***")
      pca_plot <- plotPca(
        normalized_s4,
        grouping_variable = aesthetics$color_var,
        label_column = "",
        shape_variable = aesthetics$shape_var,
                    title = "", 
        font_size = 8
      )
      
      if (!is.null(qc_dir)) {
        pca_path <- file.path(qc_dir, "post_norm_pca.png")
        ggplot2::ggsave(pca_path, pca_plot, width = 8, height = 6, dpi = 150)
        norm_data$qc_plot_paths$post_normalization$pca <- pca_path
      }
      
      # CLEAR MEMORY
      rm(pca_plot)
      gc()
      
      # RLE plot
      message("*** POST-NORM QC: Generating RLE plot ***")
      rle_plot <- plotRle(
        normalized_s4,
        group = aesthetics$color_var,
        yaxis_limit = c(-6, 6)
      )
      
      if (!is.null(qc_dir)) {
        rle_path <- file.path(qc_dir, "post_norm_rle.png")
        ggplot2::ggsave(rle_path, rle_plot, width = 10, height = 6, dpi = 150)
        norm_data$qc_plot_paths$post_normalization$rle <- rle_path
      }
      
      # CLEAR MEMORY
      rm(rle_plot)
      gc()
      
      # Density plot
      message("*** POST-NORM QC: Generating Density plot ***")
      
      # Regenerate minimal PCA for density plot
      pca_plot_temp <- plotPca(
        normalized_s4,
        grouping_variable = aesthetics$color_var,
        label_column = "",
        shape_variable = aesthetics$shape_var,
        title = "", 
        font_size = 8
      )
      
      density_plot <- plotPcaBox(
        pca_plot_temp,
        grouping_variable = aesthetics$color_var,
        show_legend = TRUE
      )
      
      if (!is.null(qc_dir)) {
        density_path <- file.path(qc_dir, "post_norm_density.png")
        ggplot2::ggsave(density_path, density_plot, width = 8, height = 6, dpi = 150)
        norm_data$qc_plot_paths$post_normalization$density <- density_path
      }
      
      # CLEAR MEMORY
      rm(pca_plot_temp, density_plot)
      gc()
      
      # Correlation plot - memory intensive
      message("*** POST-NORM QC: Generating Pearson correlation plot ***")
      correlation_plot <- plotPearson(
        normalized_s4,
        correlation_group = aesthetics$color_var,
        exclude_pool_samples = TRUE
      )
      
      if (!is.null(qc_dir)) {
        corr_path <- file.path(qc_dir, "post_norm_correlation.png")
        ggplot2::ggsave(corr_path, correlation_plot, width = 10, height = 8, dpi = 150)
        norm_data$qc_plot_paths$post_normalization$correlation <- corr_path
      }
      
      # CLEAR MEMORY
      rm(correlation_plot)
      gc()
      
      # NO S4 POPULATION
      
      # MEMORY CLEANUP
      message("*** POST-NORM QC: Running garbage collection ***")
      gc()
      
      message("Post-normalization QC plots generated successfully")
    }
    
    # Helper function to generate RUV-corrected QC plots
    # MEMORY OPTIMIZED: Save to disk immediately, clear from memory
    generateRuvCorrectedQc <- function(ruv_corrected_s4) {
      message("=== GENERATING RUV-CORRECTED QC PLOTS ===")
      
      # Get QC directory for saving plots
      qc_dir <- experiment_paths$protein_qc_dir
      
      aesthetics <- getPlotAesthetics()
      
      # Initialize plot paths storage if needed
      if (is.null(norm_data$qc_plot_paths)) {
        norm_data$qc_plot_paths <- list(
          post_filtering = list(),
          post_normalization = list(),
          ruv_corrected = list()
        )
      }
      
      # PCA plot
      message("*** RUV QC: Generating PCA plot ***")
      pca_plot <- plotPca(
        ruv_corrected_s4,
        grouping_variable = aesthetics$color_var,
        label_column = "",
        shape_variable = aesthetics$shape_var,
                    title = "",
        font_size = 8
      )
      
      if (!is.null(qc_dir)) {
        pca_path <- file.path(qc_dir, "ruv_corrected_pca.png")
        ggplot2::ggsave(pca_path, pca_plot, width = 8, height = 6, dpi = 150)
        norm_data$qc_plot_paths$ruv_corrected$pca <- pca_path
      }
      
      # CLEAR MEMORY
      rm(pca_plot)
      gc()
      
      # RLE plot
      message("*** RUV QC: Generating RLE plot ***")
      rle_plot <- plotRle(
        ruv_corrected_s4,
        group = aesthetics$color_var, 
        yaxis_limit = c(-6, 6)
      )
      
      if (!is.null(qc_dir)) {
        rle_path <- file.path(qc_dir, "ruv_corrected_rle.png")
        ggplot2::ggsave(rle_path, rle_plot, width = 10, height = 6, dpi = 150)
        norm_data$qc_plot_paths$ruv_corrected$rle <- rle_path
      }
      
      # CLEAR MEMORY
      rm(rle_plot)
      gc()
      
      # Density plot
      message("*** RUV QC: Generating Density plot ***")
      
      # Regenerate minimal PCA
      pca_plot_temp <- plotPca(
        ruv_corrected_s4,
        grouping_variable = aesthetics$color_var,
        label_column = "",
        shape_variable = aesthetics$shape_var,
        title = "",
        font_size = 8
      )
      
      density_plot <- plotPcaBox(
        pca_plot_temp,
        grouping_variable = aesthetics$color_var,
        show_legend = TRUE
      )
      
      if (!is.null(qc_dir)) {
        density_path <- file.path(qc_dir, "ruv_corrected_density.png")
        ggplot2::ggsave(density_path, density_plot, width = 8, height = 6, dpi = 150)
        norm_data$qc_plot_paths$ruv_corrected$density <- density_path
      }
      
      # CLEAR MEMORY
      rm(pca_plot_temp, density_plot)
      gc()
      
      # Correlation plot - memory intensive
      message("*** RUV QC: Generating Pearson correlation plot ***")
      correlation_plot <- plotPearson(
        ruv_corrected_s4,
        correlation_group = aesthetics$color_var,
        exclude_pool_samples = TRUE
      )
      
      if (!is.null(qc_dir)) {
        corr_path <- file.path(qc_dir, "ruv_corrected_correlation.png")
        ggplot2::ggsave(corr_path, correlation_plot, width = 10, height = 8, dpi = 150)
        norm_data$qc_plot_paths$ruv_corrected$correlation <- corr_path
      }
      
      # CLEAR MEMORY
      rm(correlation_plot)
      gc()
      
      # NO S4 POPULATION
      
      # MEMORY CLEANUP
      message("*** RUV QC: Running garbage collection ***")
      gc()
      
      message("RUV-corrected QC plots generated successfully")
    }
    
    # Auto-trigger pre-normalization QC when normalization tab is clicked
    # Observe tab selection passed from parent workflow
    if (!is.null(selected_tab)) {
      shiny::observeEvent(selected_tab(), {
        
        # Only trigger if normalization tab is selected
        if (!is.null(selected_tab()) && selected_tab() == "normalization") {
          
          message("=== NORMALIZATION TAB CLICKED ===")
          message(sprintf("workflow_data$state_manager is NULL: %s", is.null(workflow_data$state_manager)))
          
          if (!is.null(workflow_data$state_manager)) {
            current_state <- workflow_data$state_manager$current_state
            
            # Define valid states where this tab can be active
            valid_states_for_norm_tab <- c("protein_replicate_filtered", "normalized", "ruv_corrected", "correlation_filtered")
            
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
    shiny::observeEvent(c(input$color_variable, input$shape_variable), {
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
    shiny::observeEvent(input$run_normalization, {
      message("=== NORMALIZATION BUTTON CLICKED ===")
      message("Starting normalization workflow...")
      
      # MEMORY MONITORING: Check memory at start of workflow
      checkMemoryUsage(threshold_gb = 8, context = "Normalization Start")
      
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
          
          # CHECK IF RUV SHOULD BE SKIPPED
          message(sprintf("*** DEBUG: input$ruv_mode value = '%s' ***", input$ruv_mode))
          message(sprintf("*** DEBUG: Checking if '%s' == 'skip': %s ***", input$ruv_mode, input$ruv_mode == "skip"))
          
          if (input$ruv_mode == "skip") {
            message("*** RUV MODE: SKIP - Bypassing RUV-III correction ***")
            
            # Use normalized data directly as final object
            norm_data$ruv_normalized_obj <- normalized_s4
            norm_data$ruv_complete <- TRUE
            norm_data$best_k <- NA
            norm_data$control_genes_index <- NA
            norm_data$ruv_optimization_result <- list(
              best_percentage = NA,
              best_k = NA,
              best_control_genes_index = NA,
              best_cancor_plot = NULL,
              ruv_skipped = TRUE,
              skip_reason = "User selected skip due to dataset constraints"
            )
            
            # ✅ CRITICAL: Store in workflow_data so session summary can find it
            workflow_data$ruv_optimization_result <- norm_data$ruv_optimization_result
            message("*** RUV SKIP: Stored skip result in workflow_data for session summary ***")
            
            # ✅ CRITICAL: Save to file to overwrite any old RUV results
            if (!is.null(experiment_paths$source_dir)) {
              tryCatch({
                ruv_file <- file.path(experiment_paths$source_dir, "ruv_optimization_results.RDS")
                saveRDS(norm_data$ruv_optimization_result, ruv_file)
                message(sprintf("*** RUV SKIP: Saved skip result to file: %s (overwrites old results) ***", ruv_file))
              }, error = function(e) {
                message(sprintf("*** RUV SKIP: Warning - could not save skip result file: %s ***", e$message))
              })
            }
            
            message("*** RUV SKIP: Using normalized data directly for correlation filtering ***")
            
            # Save state as "normalized" instead of "ruv_corrected"
            message("*** RUV SKIP: Saving state to R6 state manager ***")
            tryCatch({
              workflow_data$state_manager$saveState(
                state_name = "normalized",
                s4_data_object = normalized_s4,
                config_object = list(
                  norm_method = input$norm_method,
                  ruv_mode = "skip",
                  ruv_applied = FALSE,
                  ruv_k = NA,
                  percentage_as_neg_ctrl = NA
                ),
                description = "Post-normalization complete: RUV-III skipped by user"
              )
              message("*** RUV SKIP: State saved successfully ***")
            }, error = function(e) {
              message(paste("*** WARNING: Could not save state to R6 manager:", e$message, "***"))
            })
            
            # Skip to Step 6B for composite QC figure (2-column layout)
            message("*** RUV SKIP: Proceeding to QC figure generation (2-column layout) ***")
            
          } else {
            # Normal RUV workflow - Steps 3-5
            
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
            
            # Track protein count after RUV filtering
            if (is.null(workflow_data$protein_counts)) {
              workflow_data$protein_counts <- list()
            }
            workflow_data$protein_counts$after_ruv_filtering <- ruvfilt_protein_count
            message(sprintf("*** STEP 5: Tracked protein count after RUV: %d ***", ruvfilt_protein_count))
            
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
          
          } # End of else block (normal RUV workflow)
          
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
          # Skip this step if RUV was skipped
          if (input$ruv_mode != "skip") {
            tryCatch({
              generateRuvCorrectedQc(ruv_corrected_s4_clean)
              message("*** STEP 6A: RUV-corrected QC plots completed ***")
            }, error = function(e) {
              message(paste("Warning: generateRuvCorrectedQc failed:", e$message))
              message("*** STEP 6A: Continuing without RUV QC plots ***")
            })
            
            # Save cancor plot to disk if available
            tryCatch({
              if (!is.null(norm_data$ruv_optimization_result$best_cancor_plot)) {
                qc_dir <- experiment_paths$protein_qc_dir
                if (!is.null(qc_dir)) {
                  cancor_path <- file.path(qc_dir, "ruv_corrected_cancor.png")
                  ggplot2::ggsave(
                    cancor_path
                    , norm_data$ruv_optimization_result$best_cancor_plot
                    , width = 8
                    , height = 6
                    , dpi = 150
                  )
                  norm_data$qc_plot_paths$ruv_corrected$cancor <- cancor_path
                  message("*** STEP 6A: Cancor plot saved to disk ***")
                }
              } else {
                message("*** STEP 6A: No cancor plot available to save ***")
              }
            }, error = function(e) {
              message(paste("Warning: Could not save cancor plot:", e$message))
            })
          } else {
            message("*** STEP 6A: Skipping RUV-corrected QC plots (RUV was skipped) ***")
          }
          
          # STEP 6B: Generate composite QC figure from saved images
          message("*** STEP 6B: Generating composite QC figure from saved images ***")
          
          qc_dir <- experiment_paths$protein_qc_dir
          
          # Define file lists based on RUV mode
          # Structure: files grouped by plot type (row) then by stage (column)
          if (input$ruv_mode == "skip") {
             ncol_composite <- 2
             # 2 columns: Pre-Norm, Post-Norm
             # Order by row: PCA, Density, RLE, Correlation (no Cancor for skip mode)
             file_names <- c(
               "pre_norm_pca.png", "post_norm_pca.png"
               , "pre_norm_density.png", "post_norm_density.png"
               , "pre_norm_rle.png", "post_norm_rle.png"
               , "pre_norm_correlation.png", "post_norm_correlation.png"
             )
             # Row labels for 2-column mode
             row_labels <- list(
               pca = c("a)", "b)")
               , density = c("c)", "d)")
               , rle = c("e)", "f)")
               , correlation = c("g)", "h)")
             )
             column_labels <- c("Pre-Normalisation", "Post-Normalisation")
          } else {
             ncol_composite <- 3
             # 3 columns: Pre, Post, RUV
             # Order by row: PCA, Density, RLE, Correlation, Cancor
             file_names <- c(
               "pre_norm_pca.png", "post_norm_pca.png", "ruv_corrected_pca.png"
               , "pre_norm_density.png", "post_norm_density.png", "ruv_corrected_density.png"
               , "pre_norm_rle.png", "post_norm_rle.png", "ruv_corrected_rle.png"
               , "pre_norm_correlation.png", "post_norm_correlation.png", "ruv_corrected_correlation.png"
               , NA, NA, "ruv_corrected_cancor.png"  # Cancor only exists for RUV column
             )
             # Row labels for 3-column mode
             row_labels <- list(
               pca = c("a)", "b)", "c)")
               , density = c("d)", "e)", "f)")
               , rle = c("g)", "h)", "i)")
               , correlation = c("j)", "k)", "l)")
               , cancor = c("", "", "m)")  # Only RUV column has cancor
             )
             column_labels <- c("Pre-Normalisation", "Post-Normalisation", "RUV-Corrected")
          }
          
          # Full paths
          if (!is.null(qc_dir)) {
            # Build full paths, handling NA placeholders for empty slots
            plot_files <- sapply(file_names, function(fn) {
              if (is.na(fn)) NA else file.path(qc_dir, fn)
            })
            
            # Output path
            composite_path <- file.path(qc_dir, "composite_QC_figure.png")
            
            # Generate composite with titles
            tryCatch({
              generateCompositeFromFiles(
                plot_files
                , composite_path
                , ncol = ncol_composite
                , row_labels = row_labels
                , column_labels = column_labels
              )
              message("*** STEP 6B: Composite QC figure saved to protein_qc_dir ***")
            }, error = function(e) {
              message(paste("Warning: Could not generate composite QC figure:", e$message))
            })
          }
          
          # Clear QC_composite_figure if it exists (legacy cleanup)
          norm_data$QC_composite_figure <- NULL
          
          # MEMORY OPTIMIZATION: Clear individual QC plots from memory since they're saved to disk
          # Keep only what's needed for UI display (the plots are also saved as individual PNGs)
          message("*** STEP 6: Clearing redundant plot objects from memory ***")
          # Note: We keep norm_data$qc_plots for UI rendering, but clear the S4 composite
          gc()
          
          # Mark normalization as complete regardless of QC plot generation
          norm_data$ruv_complete <- TRUE
          message("*** STEP 6: RUV-corrected workflow completed ***")
          
          # CRITICAL FIX: Only enable correlation filtering tab, NOT differential expression
          message("*** STEP 7: Enabling correlation filtering step ***")
          # Remove premature DE enablement - correlation filtering must happen first
          message("*** STEP 7: Normalization and RUV workflow completed - ready for correlation filtering ***")
          
        })
        
        # Update notification based on RUV mode
        notification_msg <- if (input$ruv_mode == "skip") {
          "Normalization completed (RUV skipped)! Check the 'Post-Normalisation QC' tab for filtering summary, then proceed to 'Correlation Filtering' tab for the final step."
        } else {
          "Normalization and RUV correction completed! Check the 'Post-Normalisation QC' tab for filtering summary, then proceed to 'Correlation Filtering' tab for the final step."
        }
        
        shiny::showNotification(
          notification_msg,
          type = "message",
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
    shiny::observeEvent(input$apply_correlation_filter, {
      message("=== CORRELATION FILTERING BUTTON CLICKED (DEBUG66 ACTIVE) ===")
      
      # MEMORY MONITORING: Check memory at start of correlation filtering
      checkMemoryUsage(threshold_gb = 8, context = "Correlation Filtering Start")
      
      tryCatch({
        # Get RUV-corrected object (should exist after normalization)
        ruv_s4 <- norm_data$ruv_normalized_obj
        if (is.null(ruv_s4)) {
          stop("RUV correction must be completed before correlation filtering")
        }
        
        message(sprintf("--- DEBUG66 [mod_prot_norm]: Starting Correlation Filter Flow ---"))
        message(sprintf("   [mod_prot_norm] Input Object (ruv_s4) Dimensions: %d proteins x %d samples", 
                        nrow(ruv_s4@protein_quant_table), 
                        ncol(ruv_s4@protein_quant_table) - 1)) # Approx
        
        shiny::withProgress(message = "Applying correlation filter...", value = 0, {
          
          # Step 1: Calculate correlation vector (chunk 28)
          shiny::incProgress(0.3, detail = "Calculating sample correlations...")
          message("*** CORRELATION STEP 1: Calculating sample correlations ***")
          
          message("   [mod_prot_norm] Calling pearsonCorForSamplePairs...")
          start_time_step1 <- Sys.time()
          correlation_vec <- pearsonCorForSamplePairs(
            ruv_s4,
            tech_rep_remove_regex = "pool",
            correlation_group = getRuvGroupingVariable()
          )
          end_time_step1 <- Sys.time()
          message(sprintf("   [mod_prot_norm] pearsonCorForSamplePairs returned. Duration: %.2f secs", 
                          as.numeric(difftime(end_time_step1, start_time_step1, units = "secs"))))
          message(sprintf("   [mod_prot_norm] Correlation Vector Size: %d rows", nrow(correlation_vec)))
          
          norm_data$correlation_vector <- correlation_vec
          norm_data$correlation_threshold <- input$min_pearson_correlation_threshold
          message("*** CORRELATION STEP 1: Sample correlations calculated ***")
          
          # Step 2: Apply correlation threshold filter (chunk 28)
          shiny::incProgress(0.4, detail = "Filtering low-correlation samples...")
          message("*** CORRELATION STEP 2: Applying correlation threshold filter ***")
          
          message("   [mod_prot_norm] Calling filterSamplesByProteinCorrelationThreshold...")
          start_time_step2 <- Sys.time()
          final_s4_for_de <- filterSamplesByProteinCorrelationThreshold(
            ruv_s4,
            pearson_correlation_per_pair = correlation_vec,
            min_pearson_correlation_threshold = input$min_pearson_correlation_threshold
          )
          end_time_step2 <- Sys.time()
          message(sprintf("   [mod_prot_norm] filterSamplesByProteinCorrelationThreshold returned. Duration: %.2f secs", 
                          as.numeric(difftime(end_time_step2, start_time_step2, units = "secs"))))
          
          norm_data$correlation_filtered_obj <- final_s4_for_de
          message("*** CORRELATION STEP 2: Correlation filtering applied ***")
          
          # MEMORY CLEANUP: Force garbage collection after filtering
          message("*** CORRELATION STEP 2: Running garbage collection ***")
          gc()
          
          # Step 3: Update protein filtering tracking (chunk 28)
          shiny::incProgress(0.2, detail = "Updating tracking...")
          message("*** CORRELATION STEP 3: Updating protein filtering tracking ***")
          
          # Fix Issue 3: Capture the filtering plot for final QC display
          # Wrap in tryCatch to prevent hangs and allow workflow to continue
          final_filtering_plot <- tryCatch({
            message("   [mod_prot_norm] Calling updateProteinFiltering...")
            result <- updateProteinFiltering(
              data = final_s4_for_de@protein_quant_table,
              step_name = "12_correlation_filtered",
              omic_type = omic_type,
              experiment_label = experiment_label,
              return_grid = TRUE,
              overwrite = TRUE
            )
            message("   [mod_prot_norm] updateProteinFiltering returned.")
            result
          }, error = function(e) {
            message(paste("*** WARNING: updateProteinFiltering failed:", e$message, "***"))
            message("*** CORRELATION STEP 3: Continuing without filtering plot update ***")
            return(NULL)
          })
          
          # Store the complete filtering progression for final QC display (if available)
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
          
          # MEMORY CLEANUP: Force garbage collection before file saves
          message("*** CORRELATION STEP 4: Running garbage collection before saves ***")
          gc()
          
          tryCatch({
            # Determine if RUV was skipped for file naming
            ruv_was_skipped <- isTRUE(norm_data$ruv_optimization_result$ruv_skipped) || 
                              isTRUE(workflow_data$ruv_optimization_result$ruv_skipped)
            
            # Use conditional file names based on RUV status
            if (ruv_was_skipped) {
              tsv_filename <- "normalised_results_cln_with_replicates.tsv"
              rds_filename <- "normalised_results_cln_with_replicates.RDS"
              message("*** CORRELATION STEP 4: RUV was skipped - using 'normalised' file names ***")
            } else {
              tsv_filename <- "ruv_normalised_results_cln_with_replicates.tsv"
              rds_filename <- "ruv_normalised_results_cln_with_replicates.RDS"
              message("*** CORRELATION STEP 4: RUV was applied - using 'ruv_normalised' file names ***")
            }
            
            # Save TSV file
            if (!is.null(experiment_paths) && "protein_qc_dir" %in% names(experiment_paths)) {
              # Use readr::write_tsv instead of vroom::vroom_write to avoid potential crashes with large datasets on Windows
              readr::write_tsv(
                final_s4_for_de@protein_quant_table,
                file.path(experiment_paths$protein_qc_dir, tsv_filename)
              )
              
              # Save RDS file
              saveRDS(
                final_s4_for_de,
                file.path(experiment_paths$protein_qc_dir, rds_filename)
              )
              
              message(sprintf("*** CORRELATION STEP 4: Saved files: %s, %s ***", tsv_filename, rds_filename))
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
        
        # Calculate final protein count
        final_protein_count <- length(unique(final_s4_for_de@protein_quant_table$Protein.Ids))
        
        # Track protein count after correlation filtering (final count for DE)
        if (is.null(workflow_data$protein_counts)) {
          workflow_data$protein_counts <- list()
        }
        workflow_data$protein_counts$final_for_de <- final_protein_count
        message(sprintf("*** CORRELATION: Tracked final protein count for DE: %d ***", final_protein_count))
        
        # Update summary display
        correlation_summary <- sprintf(
          "Correlation filtering completed successfully!\n\nThreshold: %.2f\nProteins remaining: %d\nSamples remaining: %d\n\nReady for differential expression analysis.",
          input$min_pearson_correlation_threshold,
          final_protein_count,
          # Fix Issue 2: Use correct column counting for samples 
          # Sample columns are all columns except the protein ID column
          length(setdiff(colnames(final_s4_for_de@protein_quant_table), final_s4_for_de@protein_id_column))
        )
        output$correlation_filter_summary <- shiny::renderText(correlation_summary)
        
        norm_data$correlation_filtering_complete <- TRUE
        
        shiny::showNotification(
          "Correlation filtering completed! Ready for differential expression analysis.",
          type = "message",
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
    
    # NEW: Skip Correlation Filtering Button Logic
    shiny::observeEvent(input$skip_correlation_filter, {
      message("=== SKIP CORRELATION FILTERING BUTTON CLICKED (DEBUG66 ACTIVE) ===")
      
      # MEMORY MONITORING
      checkMemoryUsage(threshold_gb = 8, context = "Skip Correlation Filtering Start")
      
      tryCatch({
        # Get RUV-corrected object (should exist after normalization)
        ruv_s4 <- norm_data$ruv_normalized_obj
        if (is.null(ruv_s4)) {
          stop("RUV correction must be completed before proceeding")
        }
        
        message(sprintf("--- DEBUG66 [mod_prot_norm]: Skipping Correlation Filter Flow ---"))
        message(sprintf("   [mod_prot_norm] Input Object (ruv_s4) Dimensions: %d proteins x %d samples", 
                        nrow(ruv_s4@protein_quant_table), 
                        ncol(ruv_s4@protein_quant_table) - 1))
        
        shiny::withProgress(message = "Skipping filter and saving...", value = 0, {
          
          # Step 1: Bypass filtering
          shiny::incProgress(0.3, detail = "Bypassing filter...")
          message("*** SKIP CORRELATION: Bypassing sample filtering ***")
          
          # Use the RUV object directly as the final object
          final_s4_for_de <- ruv_s4
          
          # No correlation vector needed when skipping
          norm_data$correlation_vector <- NULL
          norm_data$correlation_threshold <- NULL
          norm_data$correlation_filtered_obj <- final_s4_for_de
          
          message("*** SKIP CORRELATION: Object passed through without filtering ***")
          
          # MEMORY CLEANUP
          gc()
          
          # Step 2: Update protein filtering tracking
          shiny::incProgress(0.2, detail = "Updating tracking...")
          message("*** SKIP CORRELATION: Updating protein filtering tracking ***")
          
          # Capture the filtering plot for final QC display
          final_filtering_plot <- tryCatch({
            message("   [mod_prot_norm] Calling updateProteinFiltering (skip mode)...")
            result <- updateProteinFiltering(
              data = final_s4_for_de@protein_quant_table,
              step_name = "12_correlation_filtered",
              omic_type = omic_type,
              experiment_label = experiment_label,
              return_grid = TRUE,
              overwrite = TRUE
            )
            message("   [mod_prot_norm] updateProteinFiltering returned.")
            result
          }, error = function(e) {
            message(paste("*** WARNING: updateProteinFiltering failed:", e$message, "***"))
            return(NULL)
          })
          
          # Store the complete filtering progression
          norm_data$final_filtering_plot <- final_filtering_plot
          
          # Generate final QC plot (same as standard flow)
          tryCatch({
            aesthetics <- getPlotAesthetics()
            norm_data$final_qc_plot <- plotPca(
              final_s4_for_de,
              grouping_variable = aesthetics$color_var,
              label_column = "",
              shape_variable = aesthetics$shape_var,
              title = "Final Data (Correlation Filter Skipped)",
              font_size = 8
            )
          }, error = function(e) {
            message(paste("Error generating final QC plot:", e$message))
          })
          
          # Step 3: Save final results as TSV and RDS (Same as Apply button)
          shiny::incProgress(0.2, detail = "Saving results...")
          message("*** SKIP CORRELATION: Saving final results ***")
          
          # MEMORY CLEANUP
          gc()
          
          tryCatch({
            # Determine if RUV was skipped for file naming
            ruv_was_skipped <- isTRUE(norm_data$ruv_optimization_result$ruv_skipped) || 
                              isTRUE(workflow_data$ruv_optimization_result$ruv_skipped)
            
            # Use conditional file names based on RUV status
            if (ruv_was_skipped) {
              tsv_filename <- "normalised_results_cln_with_replicates.tsv"
              rds_filename <- "normalised_results_cln_with_replicates.RDS"
              message("*** SKIP CORRELATION: RUV was skipped - using 'normalised' file names ***")
            } else {
              tsv_filename <- "ruv_normalised_results_cln_with_replicates.tsv"
              rds_filename <- "ruv_normalised_results_cln_with_replicates.RDS"
              message("*** SKIP CORRELATION: RUV was applied - using 'ruv_normalised' file names ***")
            }
            
            # Save TSV file
            if (!is.null(experiment_paths) && "protein_qc_dir" %in% names(experiment_paths)) {
              # Use readr::write_tsv instead of vroom::vroom_write
              readr::write_tsv(
                final_s4_for_de@protein_quant_table,
                file.path(experiment_paths$protein_qc_dir, tsv_filename)
              )
              
              # Save RDS file
              saveRDS(
                final_s4_for_de,
                file.path(experiment_paths$protein_qc_dir, rds_filename)
              )
              
              message(sprintf("*** SKIP CORRELATION: Saved files: %s, %s ***", tsv_filename, rds_filename))
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
          config_object = list(min_pearson_correlation_threshold = 0, skipped = TRUE),
          description = "Skipped final sample correlation filter"
        )
        
        # DEBUG66: CRITICAL STATE UPDATE TRIGGER
        cat("--- Entering STATE UPDATE TRIGGER setting (SKIP MODE) ---\n")
        cat("   STATE UPDATE Step: Setting workflow_data$state_update_trigger...\n")
        old_trigger_value <- workflow_data$state_update_trigger
        new_trigger_value <- Sys.time()
        workflow_data$state_update_trigger <- new_trigger_value
        cat("   STATE UPDATE Step: state_update_trigger SET SUCCESSFULLY\n")
        
        # NOW enable differential expression tab
        cat("   STATE UPDATE Step: Setting tab status...\n")
        workflow_data$tab_status$normalization <- "complete"
        workflow_data$tab_status$differential_expression <- "pending"
        cat(sprintf("   STATE UPDATE Step: normalization status = %s\n", workflow_data$tab_status$normalization))
        cat("--- Exiting STATE UPDATE TRIGGER setting (SKIP MODE) ---\n")
        
        # Calculate final protein count
        final_protein_count <- length(unique(final_s4_for_de@protein_quant_table$Protein.Ids))
        
        # Track protein count
        if (is.null(workflow_data$protein_counts)) {
          workflow_data$protein_counts <- list()
        }
        workflow_data$protein_counts$final_for_de <- final_protein_count
        message(sprintf("*** SKIP CORRELATION: Tracked final protein count for DE: %d ***", final_protein_count))
        
        # Update summary display
        correlation_summary <- sprintf(
          "Correlation filtering SKIPPED.\n\nAll samples retained.\nProteins remaining: %d\nSamples remaining: %d\n\nReady for differential expression analysis.",
          final_protein_count,
          length(setdiff(colnames(final_s4_for_de@protein_quant_table), final_s4_for_de@protein_id_column))
        )
        output$correlation_filter_summary <- shiny::renderText(correlation_summary)
        
        norm_data$correlation_filtering_complete <- TRUE
        
        shiny::showNotification(
          "Correlation filtering skipped! Ready for differential expression analysis.",
          type = "message",
          duration = 5
        )
        
        message("=== SKIP CORRELATION FILTERING COMPLETED SUCCESSFULLY ===")
        
      }, error = function(e) {
        message(paste("Error in skipping correlation filtering:", e$message))
        shiny::showNotification(
          paste("Error in skipping correlation filtering:", e$message),
          type = "error",
          duration = 10
        )
      })
    })
    
    # Render QC plots - Post-filtering column
    # MEMORY OPTIMIZATION: Use renderImage to load saved PNGs instead of keeping ggplot objects in memory
    
    # Define a helper to render images safely
    render_qc_image <- function(filename, alt_text) {
      shiny::renderImage({
        # Define path
        if (!is.null(experiment_paths$protein_qc_dir)) {
           img_path <- file.path(experiment_paths$protein_qc_dir, filename)
        } else {
           img_path <- ""
        }
        
        # Check if file exists
        if (img_path != "" && file.exists(img_path)) {
          list(src = img_path,
               contentType = 'image/png',
               width = "100%",
               height = "auto",
               alt = alt_text)
        } else {
          # Return a transparent 1x1 pixel or a placeholder if file doesn't exist
          # For now, we return a list that will result in a broken image icon or we can use a placeholder
          # Better: Return a list with a flag to let the UI know, or use a default placeholder
          list(src = "", alt = "Plot not generated yet")
        }
      }, deleteFile = FALSE)
    }

    output$pca_post_filtering <- render_qc_image("pre_norm_pca.png", "PCA Post-Filtering")
    output$density_post_filtering <- render_qc_image("pre_norm_density.png", "Density Post-Filtering")
    output$rle_post_filtering <- render_qc_image("pre_norm_rle.png", "RLE Post-Filtering")
    output$correlation_post_filtering <- render_qc_image("pre_norm_correlation.png", "Correlation Post-Filtering")
    
    output$pca_post_normalization <- render_qc_image("post_norm_pca.png", "PCA Post-Normalization")
    output$density_post_normalization <- render_qc_image("post_norm_density.png", "Density Post-Normalization")
    output$rle_post_normalization <- render_qc_image("post_norm_rle.png", "RLE Post-Normalization")
    output$correlation_post_normalization <- render_qc_image("post_norm_correlation.png", "Correlation Post-Normalization")

    output$pca_ruv_corrected <- render_qc_image("ruv_corrected_pca.png", "PCA RUV Corrected")
    output$density_ruv_corrected <- render_qc_image("ruv_corrected_density.png", "Density RUV Corrected")
    output$rle_ruv_corrected <- render_qc_image("ruv_corrected_rle.png", "RLE RUV Corrected")
    output$correlation_ruv_corrected <- render_qc_image("ruv_corrected_correlation.png", "Correlation RUV Corrected")

    
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
      message("--- DEBUG66 [final_qc_plot]: Rendering ---")
      
      # Fix Issue 3: Show both PCA and complete filtering progression
      if (!is.null(norm_data$final_qc_plot) && !is.null(norm_data$final_filtering_plot)) {
        message("   [final_qc_plot] Drawing PCA + Filtering Progression...")
        
        # Create a combined layout: PCA on top, filtering progression below
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 1, heights = c(0.4, 0.6))))
        
        # Top: PCA plot
        grid::pushViewport(grid::viewport(layout.pos.row = 1))
        grid::grid.draw(ggplotGrob(norm_data$final_qc_plot))
        grid::popViewport()
        
        # Bottom: Complete filtering progression 
        grid::pushViewport(grid::viewport(layout.pos.row = 2))
        
        # CAREFUL: norm_data$final_filtering_plot might be a grid already (from arrangeGrob)
        # or a ggplot object. grid::grid.draw handles both, but let's be explicit.
        grid::grid.draw(norm_data$final_filtering_plot)
        grid::popViewport()
        
        grid::popViewport()
        message("   [final_qc_plot] Render complete.")
        
      } else if (!is.null(norm_data$final_filtering_plot)) {
        # Show just filtering progression if PCA failed
        message("   [final_qc_plot] Drawing Filtering Progression ONLY...")
        grid::grid.draw(norm_data$final_filtering_plot)
        message("   [final_qc_plot] Render complete.")
        
      } else if (!is.null(norm_data$final_qc_plot)) {
        # Show just PCA if filtering progression failed
        message("   [final_qc_plot] Drawing PCA ONLY...")
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
    shiny::observeEvent(input$reset_normalization, {
      message("Resetting normalization...")
      
      tryCatch({
        # ✅ CRITICAL FIX: Revert R6 state manager to pre-normalization state
        if (!is.null(workflow_data$state_manager)) {
          # Use smart revert pattern like QC tabs - find the actual previous state
          history <- workflow_data$state_manager$getHistory()
          # Define possible pre-normalization states in reverse chronological order
          # Include both peptide-level (raw_data_s4) and protein-level (protein_s4_initial) initial states
          pre_norm_states <- c("protein_replicate_filtered", "imputed", "replicate_filtered", "sample_filtered", "protein_peptide_filtered", "intensity_filtered", "precursor_rollup", "qvalue_filtered", "raw_data_s4", "protein_s4_initial")
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
    shiny::observeEvent(input$export_filtered_session, {
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
           
           # MEMORY FIX: Do NOT export the full R6 history or complete states list
           # These objects contain copies of every previous step's S4 object (10GB+ each)
           # Instead, we only save what is strictly needed to restore the CURRENT state
           
           session_data <- list(
             # Minimal R6 State Info for restoration
             r6_current_state_name = current_state_name,
             # We do NOT save r6_complete_states or r6_state_history to prevent 50GB+ files
             
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
            ruv_applied = (input$ruv_mode != "skip"),
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
            
            # NEW: Protein counts through workflow stages
            protein_counts = workflow_data$protein_counts,
            
            # ✅ NEW: Mixed species FASTA analysis metadata for enrichment filtering
            mixed_species_analysis = workflow_data$mixed_species_analysis,
            
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
            
            # ✅ NEW: Save mixed species analysis metadata for enrichment filtering
            if (!is.null(session_data$mixed_species_analysis)) {
              mixed_species_file <- file.path(source_dir, "mixed_species_analysis.RDS")
              saveRDS(session_data$mixed_species_analysis, mixed_species_file)
              message(sprintf("*** EXPORT: Saved mixed_species_analysis.RDS (enabled: %s) ***", 
                             isTRUE(session_data$mixed_species_analysis$enabled)))
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
          type = "message",
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

