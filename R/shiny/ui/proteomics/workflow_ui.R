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
          
          # QC plot panels with 3-column layout structure
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
              )
            )
          )
        )
      ),
      
      # Tab 5: Differential Expression
      shiny::tabPanel(
        "Differential Expression",
        value = "de",
        icon = shiny::icon("chart-bar"),
        shiny::br(),
        shiny::fluidRow(
          shiny::column(3,
            shiny::wellPanel(
              shiny::h4("DE Analysis Settings"),
              shiny::textAreaInput(
                ns("formula_string"),
                "Model Formula:",
                value = "~ 0 + group",
                height = "60px"
              ),
              shiny::helpText("Linear model formula for DE analysis"),
              shiny::br(),
              shiny::numericInput(
                ns("de_q_val_thresh"),
                "Q-value Threshold:",
                value = 0.05,
                min = 0.001,
                max = 0.2,
                step = 0.005
              ),
              shiny::numericInput(
                ns("treat_lfc_cutoff"),
                "Log Fold-Change Cutoff:",
                value = 0,
                min = 0,
                max = 2,
                step = 0.1
              ),
              shiny::br(),
              shiny::actionButton(
                ns("run_de_analysis"),
                "Run DE Analysis",
                class = "btn-primary",
                width = "100%"
              ),
              shiny::br(),
              shiny::br(),
              shiny::downloadButton(
                ns("download_de_results"),
                "Download Results",
                class = "btn-success",
                width = "100%"
              )
            )
          ),
          shiny::column(9,
            shiny::tabsetPanel(
              shiny::tabPanel(
                "Volcano Plot",
                shiny::br(),
                shiny::uiOutput(ns("contrast_selector")),
                shiny::br(),
                shiny::plotOutput(ns("volcano_plot"), height = "600px")
              ),
              shiny::tabPanel(
                "Results Table",
                shiny::br(),
                shiny::uiOutput(ns("results_contrast_selector")),
                shiny::br(),
                DT::DTOutput(ns("de_results_table"))
              ),
              shiny::tabPanel(
                "MA Plot",
                shiny::br(),
                shiny::uiOutput(ns("ma_contrast_selector")),
                shiny::br(),
                shiny::plotOutput(ns("ma_plot"), height = "600px")
              )
            )
          )
        )
      ),
      
      # Tab 6: Enrichment Analysis
      shiny::tabPanel(
        "Enrichment Analysis",
        value = "enrichment",
        icon = shiny::icon("network-wired"),
        shiny::br(),
        shiny::fluidRow(
          shiny::column(3,
            shiny::wellPanel(
              shiny::h4("Enrichment Settings"),
              shiny::selectInput(
                ns("enrichment_method"),
                "Method:",
                choices = list(
                  "g:Profiler" = "gprofiler",
                  "clusterProfiler" = "clusterprofiler"
                ),
                selected = "gprofiler"
              ),
              shiny::checkboxGroupInput(
                ns("enrichment_databases"),
                "Databases:",
                choices = list(
                  "GO Biological Process" = "GO:BP",
                  "GO Molecular Function" = "GO:MF", 
                  "GO Cellular Component" = "GO:CC",
                  "KEGG Pathways" = "KEGG",
                  "Reactome" = "REAC"
                ),
                selected = c("GO:BP", "KEGG", "REAC")
              ),
              shiny::numericInput(
                ns("enrichment_q_threshold"),
                "Q-value Threshold:",
                value = 0.05,
                min = 0.001,
                max = 0.2,
                step = 0.005
              ),
              shiny::br(),
              shiny::actionButton(
                ns("run_enrichment"),
                "Run Enrichment",
                class = "btn-primary",
                width = "100%"
              )
            )
          ),
          shiny::column(9,
            shiny::tabsetPanel(
              shiny::tabPanel(
                "Enrichment Results",
                shiny::br(),
                shiny::uiOutput(ns("enrichment_contrast_selector")),
                shiny::br(),
                DT::DTOutput(ns("enrichment_results_table"))
              ),
              shiny::tabPanel(
                "Enrichment Plots",
                shiny::br(),
                shiny::plotOutput(ns("enrichment_plot"), height = "700px")
              ),
              shiny::tabPanel(
                "Network View",
                shiny::br(),
                shiny::plotOutput(ns("enrichment_network"), height = "700px")
              )
            )
          )
        )
      ),
      
      # Tab 7: Report Generation
      shiny::tabPanel(
        "Report",
        value = "report",
        icon = shiny::icon("file-alt"),
        shiny::br(),
        shiny::fluidRow(
          shiny::column(6,
            shiny::wellPanel(
              shiny::h4("Report Settings"),
              shiny::textInput(
                ns("report_title"),
                "Report Title:",
                value = "Proteomics Analysis Report"
              ),
              shiny::textInput(
                ns("report_author"),
                "Author:",
                value = ""
              ),
              shiny::checkboxGroupInput(
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
              shiny::br(),
              shiny::actionButton(
                ns("generate_report"),
                "Generate Report",
                class = "btn-primary",
                icon = shiny::icon("file-export"),
                width = "100%"
              )
            )
          ),
          shiny::column(6,
            shiny::wellPanel(
              shiny::h4("Export Options"),
              shiny::checkboxGroupInput(
                ns("export_formats"),
                "Report Formats:",
                choices = list(
                  "HTML" = "html",
                  "PDF" = "pdf",
                  "Word" = "docx"
                ),
                selected = c("html", "pdf")
              ),
              shiny::br(),
              shiny::h5("Data Exports"),
              shiny::downloadButton(
                ns("download_processed_data"),
                "Download Processed Data",
                class = "btn-info",
                width = "100%"
              ),
              shiny::br(),
              shiny::br(),
              shiny::downloadButton(
                ns("download_session"),
                "Download Session File",
                class = "btn-info", 
                width = "100%"
              ),
              shiny::br(),
              shiny::br(),
              shiny::downloadButton(
                ns("download_parameters"),
                "Download Parameters",
                class = "btn-info",
                width = "100%"
              )
            )
          )
        )
      )
    )
  )
  
  message("   proteomicsWorkflowUi Step: tagList created successfully")
  return(result)
} 