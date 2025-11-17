#' Proteomics Workflow Server
#' 
#' Server logic for the proteomics workflow interface
#' 
#' @param id Module ID
#' @param project_dirs Project directory structure from setupDirectories
#' @param omic_type The omics type (should be "proteomics")
#' @param experiment_label The experiment label for this analysis
#' 
#' @importFrom shiny moduleServer reactive reactiveValues observeEvent req
#' @importFrom logger log_info log_error
#' @export
proteomicsWorkflowServer <- function(id, project_dirs, omic_type, experiment_label, volumes = NULL) {
  message(sprintf("--- Entering proteomicsWorkflowServer ---"))
  message(sprintf("   proteomicsWorkflowServer Arg: id = %s", id))
  message(sprintf("   proteomicsWorkflowServer Arg: omic_type = %s", omic_type))
  message(sprintf("   proteomicsWorkflowServer Arg: experiment_label = %s", experiment_label))
  message(sprintf("   proteomicsWorkflowServer Arg: project_dirs type = %s, class = %s", 
                  typeof(project_dirs), class(project_dirs)))
  message(sprintf("   proteomicsWorkflowServer Arg: volumes is NULL = %s", is.null(volumes)))
  if (!is.null(volumes)) {
    message(sprintf("   proteomicsWorkflowServer Arg: volumes type = %s, class = %s", 
                    typeof(volumes), paste(class(volumes), collapse = ", ")))
  }
  
  moduleServer(id, function(input, output, session) {
    message(sprintf("   proteomicsWorkflowServer Step: Inside moduleServer function"))
    
    # Initialize reactive values to share data between tabs
    message("   proteomicsWorkflowServer Step: Creating workflow_data reactiveValues...")
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
      
      # DEBUG66: Initialize state update trigger
      state_update_trigger = NULL,
      
      # Processing logs
      processing_log = list()
    )
    
    # Reactive trigger for initializing QC modules
    qc_trigger <- reactiveVal(NULL)
    
    cat("   proteomicsWorkflowServer Step: workflow_data created successfully\n")
    cat(sprintf("   proteomicsWorkflowServer Step: state_update_trigger initialized = %s\n", !is.null(workflow_data$state_update_trigger)))
    cat(sprintf("   proteomicsWorkflowServer Step: state_manager initialized = %s\n", !is.null(workflow_data$state_manager)))
    
    # Get paths for this proteomics experiment
    # The key in project_dirs is just the omic type (no experiment label)
    paths_key <- omic_type
    log_info(paste("Looking for project_dirs key:", paths_key))
    log_info(paste("Available keys:", paste(names(project_dirs), collapse = ', ')))
    
    message(sprintf("   proteomicsWorkflowServer Step: Checking project_dirs for key: %s", paths_key))
    message(sprintf("   proteomicsWorkflowServer Step: project_dirs names: %s", 
                    paste(names(project_dirs), collapse = ", ")))
    
    if (!paths_key %in% names(project_dirs)) {
      message(sprintf("   proteomicsWorkflowServer ERROR: No directory information found for %s", paths_key))
      log_error(paste("No directory information found for", paths_key))
      return()
    }
    
    experiment_paths <- project_dirs[[paths_key]]
    message(sprintf("   proteomicsWorkflowServer Step: Retrieved experiment_paths. Type: %s", 
                    typeof(experiment_paths)))
    
    # Debug what's in experiment_paths
    if (!is.null(experiment_paths)) {
      message(sprintf("   proteomicsWorkflowServer Step: experiment_paths is list: %s", is.list(experiment_paths)))
      message(sprintf("   proteomicsWorkflowServer Step: experiment_paths names: %s", 
                      paste(names(experiment_paths), collapse = ", ")))
      if ("results_dir" %in% names(experiment_paths)) {
        message(sprintf("   proteomicsWorkflowServer Step: results_dir = %s", experiment_paths$results_dir))
      } else {
        message("   proteomicsWorkflowServer Step: results_dir NOT found in experiment_paths")
      }
    } else {
      message("   proteomicsWorkflowServer Step: experiment_paths is NULL!")
    }
    
    # ✅ FIXED: Load aa_seq_tbl_final from scripts directory if resuming session
    if (!is.null(experiment_paths) && !is.null(experiment_paths$source_dir)) {
      aa_seq_file_path <- file.path(experiment_paths$source_dir, "aa_seq_tbl_final.RDS")
      if (file.exists(aa_seq_file_path)) {
        message("   proteomicsWorkflowServer Step: Loading existing aa_seq_tbl_final from scripts directory for session resumption")
        log_info("Loading existing aa_seq_tbl_final from scripts directory for session resumption")
        tryCatch({
          aa_seq_tbl_final <- readRDS(aa_seq_file_path)
          workflow_data$aa_seq_tbl_final <- aa_seq_tbl_final
          assign("aa_seq_tbl_final", aa_seq_tbl_final, envir = .GlobalEnv)
          log_info(sprintf("Successfully loaded aa_seq_tbl_final with %d sequences", nrow(aa_seq_tbl_final)))
        }, error = function(e) {
          log_warn(paste("Error loading aa_seq_tbl_final:", e$message))
        })
      } else {
        message("   proteomicsWorkflowServer Step: No existing aa_seq_tbl_final found in scripts directory")
      }
    }
    
    # Tab 1: Setup & Import
    message("   proteomicsWorkflowServer Step: Calling setupImportServer...")
    message(sprintf("   proteomicsWorkflowServer Step: Passing volumes to setupImportServer. is.null = %s, type = %s", 
                    is.null(volumes), typeof(volumes)))
    setupImportServer("setup_import", workflow_data, experiment_paths, volumes)
    
    # Tab 2: Design Matrix
    message("   proteomicsWorkflowServer Step: Calling designMatrixAppletServer...")
    designMatrixAppletServer("design_matrix", workflow_data, experiment_paths, volumes, qc_trigger = qc_trigger)
    
    # Tab 3: Quality Control
    qualityControlAppletServer("quality_control", workflow_data, experiment_paths, omic_type, experiment_label, qc_trigger = qc_trigger)
    
    # Create reactive for selected tab to pass to normalization module
    selected_tab <- reactive({
      input$workflow_tabs
    })
    
    # Tab 4: Normalization
    message("   proteomicsWorkflowServer Step: Calling normalizationAppletServer...")
    normalizationAppletServer("normalization", workflow_data, experiment_paths, omic_type, experiment_label, selected_tab)
    
    # Tab 5: Differential Expression
    message("   proteomicsWorkflowServer Step: Calling differentialExpressionAppletServer...")
    differentialExpressionAppletServer("differential_expression", workflow_data, experiment_paths, omic_type, experiment_label, selected_tab)
    
    # Tab 6: Enrichment Analysis
    message("   proteomicsWorkflowServer Step: Calling enrichmentAnalysisAppletServer...")
    enrichmentAnalysisAppletServer("enrichment_analysis", workflow_data, experiment_paths, omic_type, experiment_label, selected_tab)
    
    # Tab 7: Session Summary & Report Generation
    message("   proteomicsWorkflowServer Step: Calling sessionSummaryServer...")
    sessionSummaryServer("session_summary", project_dirs, omic_type, experiment_label, workflow_data)
    
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
    message("   proteomicsWorkflowServer Step: Returning workflow_data")
    return(workflow_data)
  })
} 