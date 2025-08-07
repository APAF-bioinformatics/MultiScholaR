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
  }, error = function(e) {
    message(sprintf("   proteomicsWorkflowUi: ERROR calling setupImportUi: %s", e$message))
    stop(e)
  })
  
  # Create design matrix tab
  message("   proteomicsWorkflowUi Step: Calling designMatrixAppletUI...")
  tryCatch({
    design_matrix_content <- designMatrixAppletUI(ns("design_matrix"))
    message(sprintf("   proteomicsWorkflowUi: designMatrixAppletUI returned successfully"))
  }, error = function(e) {
    message(sprintf("   proteomicsWorkflowUi: ERROR calling designMatrixAppletUI: %s", e$message))
    stop(e)
  })
  
  # Test other UI functions
  message("   proteomicsWorkflowUi Step: Testing other UI functions...")
  tryCatch({
    qc_content <- qualityControlAppletUI(ns("quality_control"))
    message(sprintf("   proteomicsWorkflowUi: qualityControlAppletUI returned successfully"))
  }, error = function(e) {
    message(sprintf("   proteomicsWorkflowUi: ERROR calling qualityControlAppletUI: %s", e$message))
    stop(e)
  })
  
  tryCatch({
    norm_content <- normalizationAppletUI(ns("normalization"))
    message(sprintf("   proteomicsWorkflowUi: normalizationAppletUI returned successfully"))
  }, error = function(e) {
    message(sprintf("   proteomicsWorkflowUi: ERROR calling normalizationAppletUI: %s", e$message))
    stop(e)
  })
  
  tryCatch({
    de_content <- differentialExpressionAppletUI(ns("differential_expression"))
    message(sprintf("   proteomicsWorkflowUi: differentialExpressionAppletUI returned successfully"))
  }, error = function(e) {
    message(sprintf("   proteomicsWorkflowUi: ERROR calling differentialExpressionAppletUI: %s", e$message))
    stop(e)
  })
  
  tryCatch({
    enrich_content <- enrichmentAnalysisAppletUI(ns("enrichment_analysis"))
    message(sprintf("   proteomicsWorkflowUi: enrichmentAnalysisAppletUI returned successfully"))
  }, error = function(e) {
    message(sprintf("   proteomicsWorkflowUi: ERROR calling enrichmentAnalysisAppletUI: %s", e$message))
    stop(e)
  })
  
  tryCatch({
    session_content <- sessionSummaryUI(ns("session_summary"))
    message(sprintf("   proteomicsWorkflowUi: sessionSummaryUI returned successfully"))
  }, error = function(e) {
    message(sprintf("   proteomicsWorkflowUi: ERROR calling sessionSummaryUI: %s", e$message))
    stop(e)
  })

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
        qc_content
      ),
      
      # Tab 4: Normalization
      shiny::tabPanel(
        "Normalization",
        value = "normalization",
        icon = shiny::icon("balance-scale"),
        shiny::br(),
        norm_content
      ),
      
      # Tab 5: Differential Expression
      shiny::tabPanel(
        "Differential Expression",
        value = "de",
        icon = shiny::icon("chart-bar"),
        shiny::br(),
        de_content
      ),
      
      # Tab 6: Enrichment Analysis
      shiny::tabPanel(
        "Enrichment Analysis",
        value = "enrichment",
        icon = shiny::icon("network-wired"),
        shiny::br(),
        enrich_content
      ),
      
      # Tab 7: Session Summary & Report Generation
      shiny::tabPanel(
        "Session Summary & Report",
        value = "session_summary",
        icon = shiny::icon("file-export"),
        shiny::br(),
        session_content
      )
    )
  )
  
  message("   proteomicsWorkflowUi Step: tagList created successfully")
  return(result)
} 