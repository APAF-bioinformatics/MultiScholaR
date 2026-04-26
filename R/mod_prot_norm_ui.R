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
              "Weighted Difference (deprecated)" = "weighted_difference"
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

