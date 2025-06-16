# Proteomics Workflow UI Module

#' Proteomics Workflow UI
#' 
#' Main UI for the proteomics workflow interface
#' 
#' @param id Module ID
#' 
#' @importFrom shiny NS tagList fluidRow column wellPanel h4 uiOutput tabsetPanel
#' @importFrom shiny tabPanel icon br conditionalPanel div DTOutput hr
#' @importFrom shiny actionButton sliderInput numericInput plotOutput radioButtons
#' @importFrom shiny textAreaInput checkboxGroupInput textInput downloadButton
#' @importFrom shiny verbatimTextOutput h3 p
#' @export
proteomicsWorkflowUi <- function(id) {
  message(sprintf("--- Entering proteomicsWorkflowUi ---"))
  message(sprintf("   proteomicsWorkflowUi Arg: id = %s", id))
  
  ns <- NS(id)
  message(sprintf("   proteomicsWorkflowUi Step: NS created. Type: %s, Class: %s", 
                  typeof(ns), class(ns)))
  
  message("   proteomicsWorkflowUi Step: Building tagList...")
  
  # Create workflow progress section
  message("   proteomicsWorkflowUi Step: Creating workflow progress section...")
  workflow_progress_section <- shiny::fluidRow(
    shiny::column(12,
      shiny::wellPanel(
        shiny::h4("Workflow Progress"),
        shiny::uiOutput(ns("workflow_progress"))
      )
    )
  )
  
  # Create setup import tab
  message("   proteomicsWorkflowUi Step: Calling setupImportUi...")
  message(sprintf("   proteomicsWorkflowUi: Checking if setupImportUi exists: %s", exists("setupImportUi")))
  if (exists("setupImportUi")) {
    message(sprintf("   proteomicsWorkflowUi: setupImportUi type: %s, class: %s", 
                    typeof(setupImportUi), class(setupImportUi)))
  }
  message(sprintf("   proteomicsWorkflowUi: About to call setupImportUi with: %s", ns("setup_import")))
  
  # Try to catch the error more precisely
  tryCatch({
    setup_import_content <- setupImportUi(ns("setup_import"))
    message(sprintf("   proteomicsWorkflowUi: setupImportUi returned successfully"))
    message(sprintf("   proteomicsWorkflowUi: Return type: %s, class: %s", 
                    typeof(setup_import_content), class(setup_import_content)))
    message("   proteomicsWorkflowUi: Return structure:")
    utils::str(setup_import_content, max.level = 2)
  }, error = function(e) {
    message(sprintf("   proteomicsWorkflowUi: ERROR calling setupImportUi: %s", e$message))
    message("   proteomicsWorkflowUi: Error traceback:")
    print(traceback())
    stop(e)
  })
  
  # Create design matrix tab
  message("   proteomicsWorkflowUi Step: Calling designMatrixAppletUI...")
  design_matrix_content <- designMatrixAppletUI(ns("design_matrix"))
  
  # Create QC tab content
  message("   proteomicsWorkflowUi Step: Creating QC tab content...")
  qc_tab_content <- shiny::fluidRow(
          shiny::column(3,
            shiny::wellPanel(
              shiny::h4("QC Parameters"),
              shiny::sliderInput(
                ns("q_value_threshold"),
                "Q-value Threshold:",
                min = 0.001,
                max = 0.1,
                value = 0.01,
                step = 0.001
              ),
              shiny::sliderInput(
                ns("intensity_percentile"),
                "Intensity Cutoff Percentile:",
                min = 0,
                max = 10,
                value = 1,
                step = 0.5
              ),
              shiny::numericInput(
                ns("min_peptides_per_protein"),
                "Min Peptides per Protein:",
                value = 2,
                min = 1,
                max = 10
              ),
              shiny::numericInput(
                ns("min_peptides_per_sample"),
                "Min Peptides per Sample:",
                value = 1000,
                min = 100,
                max = 10000,
                step = 100
              ),
              shiny::br(),
              shiny::actionButton(
                ns("run_qc"),
                "Run QC Filters",
                class = "btn-primary",
                width = "100%"
              )
            )
          ),
          shiny::column(9,
            shiny::tabsetPanel(
              shiny::tabPanel(
                "Filtering Summary",
                shiny::br(),
                shiny::plotOutput(ns("filtering_summary_plot"), height = "600px")
              ),
              shiny::tabPanel(
                "QC Metrics",
                shiny::br(),
                shiny::fluidRow(
                  shiny::column(6, shiny::plotOutput(ns("rle_plot"))),
                  shiny::column(6, shiny::plotOutput(ns("pca_plot")))
                ),
                shiny::fluidRow(
                  shiny::column(6, shiny::plotOutput(ns("density_plot"))),
                  shiny::column(6, shiny::plotOutput(ns("correlation_plot")))
                )
              ),
              shiny::tabPanel(
                "Filtered Data",
                shiny::br(),
                DT::DTOutput(ns("filtered_data_preview"))
              )
            )
          )
        )
  
  # Now build the complete tagList
  message("   proteomicsWorkflowUi Step: Creating tabsetPanel...")
  result <- shiny::tagList(
    # Workflow progress indicator
    workflow_progress_section,
    
    # Main tabset for workflow steps
    shiny::tabsetPanel(
      id = ns("workflow_tabs"),
      type = "pills",
      
      # Tab 1: Setup and Data Import
      shiny::tabPanel(
        "Setup & Import",
        value = "setup",
        icon = shiny::icon("upload"),
        shiny::br(),
        setup_import_content
      ),
      
      # Tab 2: Design Matrix
      shiny::tabPanel(
        "Design Matrix",
        value = "design",
        icon = shiny::icon("table"),
        shiny::br(),
        design_matrix_content
      ),
      
      # Tab 3: Quality Control
      shiny::tabPanel(
        "Quality Control",
        value = "qc",
        icon = shiny::icon("chart-line"),
        shiny::br(),
        qualityControlAppletUI(ns("quality_control"))
      ),
      
      # Tab 4: Normalization
      shiny::tabPanel(
        "Normalization",
        value = "normalization",
        icon = shiny::icon("balance-scale"),
        shiny::br(),
        shiny::fluidRow(
            shiny::column(3,
              shiny::wellPanel(
                shiny::h4("QC Parameters"),
                shiny::sliderInput(
                  ns("q_value_threshold"),
                  "Q-value Threshold:",
                  min = 0.001,
                  max = 0.1,
                  value = 0.01,
                  step = 0.001
                ),
                shiny::sliderInput(
                  ns("intensity_percentile"),
                  "Intensity Cutoff Percentile:",
                  min = 0,
                  max = 10,
                  value = 1,
                  step = 0.5
                ),
                shiny::numericInput(
                  ns("min_peptides_per_protein"),
                  "Min Peptides per Protein:",
                  value = 2,
                  min = 1,
                  max = 10
                ),
                shiny::numericInput(
                  ns("min_peptides_per_sample"),
                  "Min Peptides per Sample:",
                  value = 1000,
                  min = 100,
                  max = 10000,
                  step = 100
                ),
                shiny::br(),
                shiny::actionButton(
                  ns("run_qc"),
                  "Run QC Filters",
                  class = "btn-primary",
                  width = "100%"
                )
              )
            ),
            shiny::column(9,
              shiny::tabsetPanel(
                shiny::tabPanel(
                  "Filtering Summary",
                  shiny::br(),
                  shiny::plotOutput(ns("filtering_summary_plot"), height = "600px")
                ),
                shiny::tabPanel(
                  "QC Metrics",
                  shiny::br(),
                  shiny::fluidRow(
                    shiny::column(6, shiny::plotOutput(ns("rle_plot"))),
                    shiny::column(6, shiny::plotOutput(ns("pca_plot")))
                  ),
                  shiny::fluidRow(
                    shiny::column(6, shiny::plotOutput(ns("density_plot"))),
                    shiny::column(6, shiny::plotOutput(ns("correlation_plot")))
                  )
                ),
                shiny::tabPanel(
                  "Filtered Data",
                  shiny::br(),
                  DT::DTOutput(ns("filtered_data_preview"))
                )
              )
            )
          )
      ),
      
      # Tab 4: Normalization
      shiny::tabPanel(
        "Normalization",
        value = "normalization",
        icon = shiny::icon("balance-scale"),
        shiny::br(),
        shiny::fluidRow(
          shiny::column(3,
            shiny::wellPanel(
              shiny::h4("Normalization Options"),
              shiny::radioButtons(
                ns("norm_method"),
                "Normalization Method:",
                choices = list(
                  "Cyclic Loess" = "cyclicloess",
                  "Quantile" = "quantile",
                  "RLE" = "rle",
                  "TMM" = "TMM"
                ),
                selected = "cyclicloess"
              ),
              shiny::hr(),
              shiny::h4("RUV-III Parameters"),
              shiny::sliderInput(
                ns("ruv_percentage"),
                "% Proteins as Negative Controls:",
                min = 1,
                max = 20,
                value = 5,
                step = 1
              ),
              numericInput(
                ns("ruv_k"),
                "Number of Factors (k):",
                value = NULL,
                min = 1,
                max = 10
              ),
              br(),
              actionButton(
                ns("run_normalization"),
                "Run Normalization",
                class = "btn-primary",
                width = "100%"
              )
            )
          ),
          column(9,
            tabsetPanel(
              tabPanel(
                "Canonical Correlation",
                br(),
                plotOutput(ns("cancor_plot"), height = "500px"),
                br(),
                verbatimTextOutput(ns("best_k_suggestion"))
              ),
              tabPanel(
                "Before/After Comparison",
                br(),
                fluidRow(
                  column(6, 
                    h4("Before Normalization"),
                    plotOutput(ns("before_norm_plots"))
                  ),
                  column(6,
                    h4("After Normalization"),
                    plotOutput(ns("after_norm_plots"))
                  )
                )
              ),
              tabPanel(
                "Normalized Data",
                br(),
                DTOutput(ns("normalized_data_preview"))
              )
            )
          )
        )
      ),
      
      # Tab 5: Differential Expression
      tabPanel(
        "Differential Expression",
        value = "de_analysis",
        icon = shiny::icon("dna"),
        br(),
        fluidRow(
          column(3,
            wellPanel(
              h4("DE Analysis Settings"),
              textAreaInput(
                ns("formula_string"),
                "Model Formula:",
                value = "~ 0 + group",
                height = "80px"
              ),
              hr(),
              h4("Contrasts"),
              uiOutput(ns("contrast_builder")),
              actionButton(
                ns("add_contrast"),
                "Add Contrast",
                icon = shiny::icon("plus"),
                class = "btn-sm"
              ),
              hr(),
              sliderInput(
                ns("de_qval_threshold"),
                "Q-value Threshold:",
                min = 0.01,
                max = 0.1,
                value = 0.05,
                step = 0.01
              ),
              sliderInput(
                ns("treat_lfc"),
                "Min Log2 Fold Change:",
                min = 0,
                max = 2,
                value = 0,
                step = 0.1
              ),
              br(),
              actionButton(
                ns("run_de_analysis"),
                "Run DE Analysis",
                class = "btn-primary",
                width = "100%"
              )
            )
          ),
          column(9,
            tabsetPanel(
              tabPanel(
                "Results Summary",
                br(),
                DTOutput(ns("de_summary_table"))
              ),
              tabPanel(
                "Volcano Plots",
                br(),
                uiOutput(ns("volcano_plot_selector")),
                plotOutput(ns("volcano_plot"), height = "600px")
              ),
              tabPanel(
                "DE Tables",
                br(),
                uiOutput(ns("de_table_selector")),
                br(),
                downloadButton(ns("download_de_results"), "Download Results"),
                br(), br(),
                DTOutput(ns("de_results_table"))
              )
            )
          )
        )
      ),
      
      # Tab 6: Enrichment Analysis
      tabPanel(
        "Enrichment Analysis",
        value = "enrichment",
        icon = shiny::icon("sitemap"),
        br(),
        fluidRow(
          column(3,
            wellPanel(
              h4("Enrichment Settings"),
              radioButtons(
                ns("enrichment_method"),
                "Enrichment Method:",
                choices = list(
                  "g:Profiler" = "gprofiler",
                  "clusterProfiler" = "clusterprofiler"
                ),
                selected = "gprofiler"
              ),
              checkboxGroupInput(
                ns("enrichment_sources"),
                "Data Sources:",
                choices = list(
                  "GO:BP" = "GO:BP",
                  "GO:MF" = "GO:MF",
                  "GO:CC" = "GO:CC",
                  "KEGG" = "KEGG",
                  "Reactome" = "REAC"
                ),
                selected = c("GO:BP", "GO:MF", "GO:CC")
              ),
              sliderInput(
                ns("enrichment_qval"),
                "Q-value Threshold:",
                min = 0.01,
                max = 0.1,
                value = 0.05,
                step = 0.01
              ),
              br(),
              actionButton(
                ns("run_enrichment"),
                "Run Enrichment",
                class = "btn-primary",
                width = "100%"
              )
            )
          ),
          column(9,
            tabsetPanel(
              tabPanel(
                "Enrichment Results",
                br(),
                uiOutput(ns("enrichment_contrast_selector")),
                br(),
                DTOutput(ns("enrichment_results_table"))
              ),
              tabPanel(
                "Enrichment Plots",
                br(),
                plotOutput(ns("enrichment_plot"), height = "700px")
              ),
              tabPanel(
                "Network View",
                br(),
                plotOutput(ns("enrichment_network"), height = "700px")
              )
            )
          )
        )
      ),
      
      # Tab 7: Report Generation
      tabPanel(
        "Report",
        value = "report",
        icon = shiny::icon("file-alt"),
        br(),
        fluidRow(
          column(6,
            wellPanel(
              h4("Report Settings"),
              textInput(
                ns("report_title"),
                "Report Title:",
                value = "Proteomics Analysis Report"
              ),
              textInput(
                ns("report_author"),
                "Author:",
                value = ""
              ),
              checkboxGroupInput(
                ns("report_sections"),
                "Include Sections:",
                choices = list(
                  "Data Summary" = "summary",
                  "Quality Control" = "qc",
                  "Normalization" = "norm",
                  "Differential Expression" = "de",
                  "Enrichment Analysis" = "enrichment",
                  "Methods" = "methods"
                ),
                selected = c("summary", "qc", "norm", "de", "enrichment", "methods")
              ),
              br(),
              actionButton(
                ns("generate_report"),
                "Generate Report",
                class = "btn-primary",
                icon = shiny::icon("file-export"),
                width = "100%"
              )
            )
          ),
          column(6,
            wellPanel(
              h4("Export Options"),
              downloadButton(
                ns("download_processed_data"),
                "Download Processed Data",
                class = "btn-info",
                style = "width: 100%; margin-bottom: 10px;"
              ),
              downloadButton(
                ns("download_session"),
                "Download R Session",
                class = "btn-info",
                style = "width: 100%; margin-bottom: 10px;"
              ),
              downloadButton(
                ns("download_workflow_params"),
                "Download Workflow Parameters",
                class = "btn-info",
                style = "width: 100%;"
              )
            )
          )
        )
      )
    )
  )
  
  message(sprintf("   proteomicsWorkflowUi Step: tagList created. Type: %s, Class: %s", 
                  typeof(result), class(result)))
  message(sprintf("--- Exiting proteomicsWorkflowUi ---"))
  
  return(result)
} 