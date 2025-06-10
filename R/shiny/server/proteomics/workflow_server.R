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
proteomicsWorkflowServer <- function(id, project_dirs, omic_type, experiment_label) {
  message(sprintf("--- Entering proteomicsWorkflowServer ---"))
  message(sprintf("   proteomicsWorkflowServer Arg: id = %s", id))
  message(sprintf("   proteomicsWorkflowServer Arg: omic_type = %s", omic_type))
  message(sprintf("   proteomicsWorkflowServer Arg: experiment_label = %s", experiment_label))
  message(sprintf("   proteomicsWorkflowServer Arg: project_dirs type = %s, class = %s", 
                  typeof(project_dirs), class(project_dirs)))
  
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
        report_generation = "disabled"
      ),
      
      # Processing logs
      processing_log = list()
    )
    message("   proteomicsWorkflowServer Step: workflow_data created successfully")
    
    # Get paths for this proteomics experiment
    paths_key <- paste0(omic_type, "_", experiment_label)
    log_info("Looking for project_dirs key: {paths_key}")
    log_info("Available keys: {paste(names(project_dirs), collapse = ', ')}")
    
    message(sprintf("   proteomicsWorkflowServer Step: Checking project_dirs for key: %s", paths_key))
    message(sprintf("   proteomicsWorkflowServer Step: project_dirs names: %s", 
                    paste(names(project_dirs), collapse = ", ")))
    
    if (!paths_key %in% names(project_dirs)) {
      message(sprintf("   proteomicsWorkflowServer ERROR: No directory information found for %s", paths_key))
      log_error("No directory information found for {paths_key}")
      return()
    }
    
    experiment_paths <- project_dirs[[paths_key]]
    message(sprintf("   proteomicsWorkflowServer Step: Retrieved experiment_paths. Type: %s", 
                    typeof(experiment_paths)))
    
    # Tab 1: Setup & Import
    message("   proteomicsWorkflowServer Step: Calling setupImportServer...")
    setupImportServer("setup_import", workflow_data, experiment_paths)
    
    # Tab 2: Design Matrix
    message("   proteomicsWorkflowServer Step: Calling designMatrixServer...")
    designMatrixServer("design_matrix", workflow_data, experiment_paths, omic_type, experiment_label)
    
    # Tab 3: Quality Control
    # qualityControlServer("quality_control", workflow_data, experiment_paths)
    
    # Tab 4: Normalization
    # normalizationServer("normalization", workflow_data, experiment_paths)
    
    # Tab 5: Differential Expression
    # differentialExpressionServer("differential_expression", workflow_data, experiment_paths)
    
    # Tab 6: Enrichment Analysis
    # enrichmentAnalysisServer("enrichment_analysis", workflow_data, experiment_paths)
    
    # Tab 7: Report Generation
    # reportGenerationServer("report_generation", workflow_data, experiment_paths)
    
    # Enable tabs based on completion status
    observeEvent(workflow_data$tab_status, {
      # Enable design matrix tab after setup is complete
      if (workflow_data$tab_status$setup_import == "complete") {
        workflow_data$tab_status$design_matrix <- "pending"
      }
      
      # Enable QC tab after design matrix is complete
      if (workflow_data$tab_status$design_matrix == "complete") {
        workflow_data$tab_status$quality_control <- "pending"
      }
      
      # And so on for other tabs...
    })
    
    # Return workflow data for potential use by parent module
    message("   proteomicsWorkflowServer Step: Returning workflow_data")
    return(workflow_data)
  })
} 