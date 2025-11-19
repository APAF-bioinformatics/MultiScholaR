#' Proteomics Workflow UI Module
#'
#' @description Main UI for the proteomics workflow interface
#'
#' @param id Module ID
#'
#' @importFrom shiny NS tagList fluidRow column wellPanel h4 uiOutput tabsetPanel
#' @importFrom shiny tabPanel icon br conditionalPanel div hr
#' @importFrom shiny actionButton sliderInput numericInput plotOutput radioButtons
#' @importFrom shiny textAreaInput checkboxGroupInput textInput downloadButton
#' @importFrom shiny verbatimTextOutput h3 p
#' @importFrom DT DTOutput
#' @export
mod_proteomics_ui <- function(id) {
  ns <- NS(id)
  
  # Create workflow progress section
  workflow_progress_section <- shiny::fluidRow(
    shiny::column(12,
      shiny::wellPanel(
        shiny::h4("Workflow Progress"),
        shiny::uiOutput(ns("workflow_progress"))
      )
    )
  )
  
  # Create setup import tab
  # Using mod_prot_import_ui
  if (exists("mod_prot_import_ui")) {
      setup_import_content <- mod_prot_import_ui(ns("setup_import"))
  } else {
      setup_import_content <- shiny::div("Import module not loaded")
  }
  
  # Create design matrix tab
  # Using mod_prot_design_ui
  if (exists("mod_prot_design_ui")) {
    design_matrix_content <- mod_prot_design_ui(ns("design_matrix"))
  } else {
    design_matrix_content <- shiny::div("Design matrix module not loaded")
  }

  # Create QC tab
  # Using mod_prot_qc_ui
  if (exists("mod_prot_qc_ui")) {
    qc_content <- mod_prot_qc_ui(ns("quality_control"))
  } else {
    qc_content <- shiny::div("QC module not loaded")
  }
  
  # Create Normalization tab
  # Using mod_prot_norm_ui
  if (exists("mod_prot_norm_ui")) {
    norm_content <- mod_prot_norm_ui(ns("normalization"))
  } else {
    norm_content <- shiny::div("Normalization module not loaded")
  }
  
  # Create DE tab
  # Using mod_prot_de_ui
  if (exists("mod_prot_de_ui")) {
    de_content <- mod_prot_de_ui(ns("differential_expression"))
  } else {
    de_content <- shiny::div("DE module not loaded")
  }
  
  # Create Enrichment tab
  # Using mod_prot_enrich_ui
  if (exists("mod_prot_enrich_ui")) {
    enrich_content <- mod_prot_enrich_ui(ns("enrichment_analysis"))
  } else {
    enrich_content <- shiny::div("Enrichment module not loaded")
  }
  
  # Create Session Summary tab
  # Using mod_prot_summary_ui
  if (exists("mod_prot_summary_ui")) {
    session_content <- mod_prot_summary_ui(ns("session_summary"))
  } else {
    session_content <- shiny::div("Session Summary module not loaded")
  }

  # Now build the complete tagList
  shiny::tagList(
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
}

#' Proteomics Workflow Server
#' 
#' @description Server logic for the proteomics workflow interface
#' 
#' @param id Module ID
#' @param project_dirs Project directory structure from setupDirectories
#' @param omic_type The omics type (should be "proteomics")
#' @param experiment_label The experiment label for this analysis
#' @param volumes Volumes for shinyFiles
#' 
#' @importFrom shiny moduleServer reactive reactiveValues observeEvent req reactiveVal
#' @importFrom logger log_info log_error
#' @export
mod_proteomics_server <- function(id, project_dirs, omic_type, experiment_label, volumes = NULL) {
  moduleServer(id, function(input, output, session) {
    
    # Initialize reactive values to share data between tabs
    workflow_data <- reactiveValues(
      # Data objects
      data_tbl = NULL,
      fasta_file_path = NULL,
      config_list = NULL,
      taxon_id = NULL,
      organism_name = NULL,
      design_matrix = NULL,
      data_cln = NULL,
      contrasts_tbl = NULL,
      
      # R6 state manager for S4 objects
      state_manager = WorkflowState$new(),
      
      # Legacy S4 object slots (can be deprecated later)
      peptide_data = NULL,
      protein_log2_quant = NULL,
      protein_data = NULL,
      ruv_normalised_for_de_analysis_obj = NULL,
      de_analysis_results_list = NULL,
      uniprot_dat_cln = NULL,
      enrichment_results = NULL,
      
      # Status tracking
      tab_status = list(
        setup_import = "pending",
        design_matrix = "disabled", 
        quality_control = "disabled",
        normalization = "disabled",
        differential_expression = "disabled",
        enrichment_analysis = "disabled",
        session_summary = "disabled"
      ),
      
      # Initialize state update trigger
      state_update_trigger = NULL,
      
      # Processing logs
      processing_log = list()
    )
    
    # Reactive trigger for initializing QC modules
    qc_trigger <- reactiveVal(NULL)
    
    # Get paths for this proteomics experiment
    paths_key <- omic_type
    
    if (!paths_key %in% names(project_dirs)) {
      log_error(paste("No directory information found for", paths_key))
      return()
    }
    
    experiment_paths <- project_dirs[[paths_key]]
    
    # Load aa_seq_tbl_final from scripts directory if resuming session
    if (!is.null(experiment_paths) && !is.null(experiment_paths$source_dir)) {
      aa_seq_file_path <- file.path(experiment_paths$source_dir, "aa_seq_tbl_final.RDS")
      if (file.exists(aa_seq_file_path)) {
        log_info("Loading existing aa_seq_tbl_final from scripts directory for session resumption")
        tryCatch({
          aa_seq_tbl_final <- readRDS(aa_seq_file_path)
          workflow_data$aa_seq_tbl_final <- aa_seq_tbl_final
          assign("aa_seq_tbl_final", aa_seq_tbl_final, envir = .GlobalEnv)
          log_info(sprintf("Successfully loaded aa_seq_tbl_final with %d sequences", nrow(aa_seq_tbl_final)))
        }, error = function(e) {
          log_warn(paste("Error loading aa_seq_tbl_final:", e$message))
        })
      }
    }
    
    # Tab 1: Setup & Import
    # Using mod_prot_import_server
    if (exists("mod_prot_import_server")) {
      mod_prot_import_server("setup_import", workflow_data, experiment_paths, volumes)
    }
    
    # Tab 2: Design Matrix
    # Using mod_prot_design_server
    if (exists("mod_prot_design_server")) {
      mod_prot_design_server("design_matrix", workflow_data, experiment_paths, volumes, qc_trigger = qc_trigger)
    }
    
    # Tab 3: Quality Control
    # Using mod_prot_qc_server
    if (exists("mod_prot_qc_server")) {
      mod_prot_qc_server("quality_control", workflow_data, experiment_paths, omic_type, experiment_label, qc_trigger = qc_trigger)
    }
    
    # Create reactive for selected tab to pass to normalization module
    selected_tab <- reactive({
      input$workflow_tabs
    })
    
    # Tab 4: Normalization
    # Using mod_prot_norm_server
    if (exists("mod_prot_norm_server")) {
      mod_prot_norm_server("normalization", workflow_data, experiment_paths, omic_type, experiment_label, selected_tab)
    }
    
    # Tab 5: Differential Expression
    # Using mod_prot_de_server
    if (exists("mod_prot_de_server")) {
      mod_prot_de_server("differential_expression", workflow_data, experiment_paths, omic_type, experiment_label, selected_tab)
    }
    
    # Tab 6: Enrichment Analysis
    # Using mod_prot_enrich_server
    if (exists("mod_prot_enrich_server")) {
      mod_prot_enrich_server("enrichment_analysis", workflow_data, experiment_paths, omic_type, experiment_label, selected_tab)
    }
    
    # Tab 7: Session Summary & Report Generation
    # Using mod_prot_summary_server
    if (exists("mod_prot_summary_server")) {
      mod_prot_summary_server("session_summary", project_dirs, omic_type, experiment_label, workflow_data)
    }
    
    observeEvent(workflow_data$tab_status$design_matrix, {
      # Enable QC tab after design matrix is complete
      if (workflow_data$tab_status$design_matrix == "complete") {
        workflow_type <- shiny::isolate(workflow_data$state_manager$workflow_type)
        log_info(paste("Workflow server: Design matrix complete. Workflow type:", workflow_type))
        
        if (workflow_type %in% c("TMT", "LFQ")) {
          # For TMT and LFQ (protein-level workflows), bypass QC and go straight to Normalization
          log_info(sprintf("%s workflow detected, bypassing QC tab.", workflow_type))
          workflow_data$tab_status$quality_control <- "complete"
          workflow_data$tab_status$normalization <- "pending"
        } else {
          # For DIA (peptide-level workflow), proceed to QC tab as normal
          workflow_data$tab_status$quality_control <- "pending"
        }
      }
    }, ignoreNULL = TRUE)
    
    observeEvent(workflow_data$tab_status$quality_control, {
      # Enable normalization tab after QC is complete
      if (workflow_data$tab_status$quality_control == "complete") {
        workflow_data$tab_status$normalization <- "pending"
      }
    }, ignoreNULL = TRUE)
    
    observeEvent(workflow_data$tab_status$normalization, {
      # Enable DE tab after normalization is complete
      if (workflow_data$tab_status$normalization == "complete") {
        workflow_data$tab_status$differential_expression <- "pending"
      }
    }, ignoreNULL = TRUE)
    
    observeEvent(workflow_data$tab_status$differential_expression, {
      # Enable enrichment analysis tab after DE is complete
      if (workflow_data$tab_status$differential_expression == "complete") {
        workflow_data$tab_status$enrichment_analysis <- "pending"
      }
    }, ignoreNULL = TRUE)
    
    observeEvent(workflow_data$tab_status$enrichment_analysis, {
      # Enable session summary after enrichment is complete
      if (workflow_data$tab_status$enrichment_analysis == "complete") {
        workflow_data$tab_status$session_summary <- "pending"
      }
    }, ignoreNULL = TRUE)
    
    # Return workflow data for potential use by parent module
    return(workflow_data)
  })
}

