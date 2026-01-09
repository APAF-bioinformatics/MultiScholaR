# ============================================================================
# mod_metabolomics.R
# ============================================================================
# Purpose: Top-level metabolomics workflow orchestrator Shiny module
#
# This module coordinates all metabolomics sub-modules in a tabbed interface,
# managing the complete workflow from import to export.
# ============================================================================

#' @title Metabolomics Workflow Orchestrator Module
#' @description Top-level Shiny module that coordinates the complete metabolomics
#'              workflow including import, design, QC, normalization, differential
#'              analysis, and summary/export.
#' @name mod_metabolomics
NULL

#' @rdname mod_metabolomics
#' @export
#' @importFrom shiny NS tagList fluidRow column h2 tabsetPanel tabPanel uiOutput icon
mod_metabolomics_ui <- function(id) {
    ns <- shiny::NS(id)
    
    shiny::tagList(
        shiny::fluidRow(
            shiny::column(12
                , shiny::h2(
                    shiny::icon("flask")
                    , " Metabolomics Workflow"
                    , style = "margin-bottom: 20px;"
                )
                
                # Workflow progress indicator
                , shiny::uiOutput(ns("workflow_progress"))
                
                # Main workflow tabs
                , shiny::tabsetPanel(
                    id = ns("metabolomics_tabs")
                    , type = "pills"
                    
                    # Tab 1: Setup & Import
                    , shiny::tabPanel(
                        title = shiny::tagList(
                            shiny::icon("file-import")
                            , " Setup & Import"
                        )
                        , value = "import"
                        , mod_metab_import_ui(ns("import"))
                    )
                    
                    # Tab 2: Design Matrix
                    , shiny::tabPanel(
                        title = shiny::tagList(
                            shiny::icon("th")
                            , " Design Matrix"
                        )
                        , value = "design"
                        , mod_metab_design_ui(ns("design"))
                    )
                    
                    # Tab 3: Quality Control
                    , shiny::tabPanel(
                        title = shiny::tagList(
                            shiny::icon("check-double")
                            , " Quality Control"
                        )
                        , value = "qc"
                        , mod_metab_qc_ui(ns("qc"))
                    )
                    
                    # Tab 4: Normalization
                    , shiny::tabPanel(
                        title = shiny::tagList(
                            shiny::icon("balance-scale")
                            , " Normalization"
                        )
                        , value = "norm"
                        , mod_metab_norm_ui(ns("norm"))
                    )
                    
                    # Tab 5: Differential Analysis
                    , shiny::tabPanel(
                        title = shiny::tagList(
                            shiny::icon("chart-line")
                            , " Differential Analysis"
                        )
                        , value = "de"
                        , mod_metab_de_ui(ns("de"))
                    )
                    
                    # Tab 6: Summary & Export
                    , shiny::tabPanel(
                        title = shiny::tagList(
                            shiny::icon("download")
                            , " Summary & Export"
                        )
                        , value = "summary"
                        , mod_metab_summary_ui(ns("summary"))
                    )
                )
            )
        )
    )
}

#' @rdname mod_metabolomics
#' @param id Module ID
#' @param project_dirs List of project directories from setup
#' @param omic_type The omics type (should be "metabolomics")
#' @param experiment_label The experiment label for this analysis
#' @param volumes Volumes for shinyFiles file browser
#' @export
#' @importFrom shiny moduleServer reactiveValues reactive observeEvent renderUI req tags reactiveVal
#' @importFrom logger log_info log_error log_warn
mod_metabolomics_server <- function(id, project_dirs, omic_type, experiment_label, volumes = NULL) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns
        
        # Initialize reactive values to share data between tabs
        workflow_data <- shiny::reactiveValues(
            # Data storage
            data_tbl = NULL
            , data_format = NULL
            , data_type = "metabolite"
            
            # Column mapping from import
            , column_mapping = NULL
            
            # Config
            , config_list = list()
            
            # R6 state manager for S4 objects
            , state_manager = WorkflowState$new()
            
            # S4 objects for metabolomics
            , metabolite_assay_data = NULL
            
            # Contrasts for DE
            , contrasts = list()
            , contrasts_tbl = NULL
            
            # Tab status tracking
            , tab_status = list(
                setup_import = "pending"
                , design_matrix = "disabled"
                , quality_control = "disabled"
                , normalization = "disabled"
                , differential_analysis = "disabled"
                , session_summary = "disabled"
            )
            
            # Initialize state update trigger
            , state_update_trigger = NULL
            
            # Processing log
            , processing_log = list()
        )
        
        # Get paths for this metabolomics experiment
        paths_key <- omic_type
        
        if (!paths_key %in% names(project_dirs)) {
            logger::log_error(paste("No directory information found for", paths_key))
            return()
        }
        
        experiment_paths <- project_dirs[[paths_key]]
        
        logger::log_info(sprintf("Metabolomics module initialized with paths for: %s", paths_key))
        
        # Reactive trigger for initializing QC modules (like proteomics pattern)
        qc_trigger <- shiny::reactiveVal(NULL)
        
        # Initialize all sub-module servers
        logger::log_info("Initializing metabolomics workflow modules")
        
        # Import module
        mod_metab_import_server(
            "import"
            , workflow_data
            , experiment_paths
            , volumes
        )
        
        # Design module
        mod_metab_design_server(
            "design"
            , workflow_data
            , experiment_paths
            , volumes = volumes
            , qc_trigger = qc_trigger
        )
        
        # QC module
        mod_metab_qc_server(
            "qc"
            , workflow_data
            , experiment_paths
            , omic_type
            , experiment_label
            , qc_trigger
        )
        
        # Normalization module
        mod_metab_norm_server(
            "norm"
            , workflow_data
            , experiment_paths
            , omic_type
            , experiment_label
            , selected_tab = shiny::reactive(input$metabolomics_tabs)
        )
        
        # DE module
        mod_metab_de_server(
            "de"
            , workflow_data
            , experiment_paths
            , omic_type
            , experiment_label
        )
        
        # Summary module
        mod_metab_summary_server(
            "summary"
            , workflow_data
            , experiment_paths
            , omic_type
            , experiment_label
        )
        
        logger::log_info("Metabolomics workflow modules initialized")
        
        # Workflow progress indicator using shared stepper component
        output$workflow_progress <- shiny::renderUI({
            # Define metabolomics workflow steps (6 steps)
            steps <- list(
                list(name = "Import", key = "setup_import", icon = "file-import")
                , list(name = "Design", key = "design_matrix", icon = "th")
                , list(name = "QC", key = "quality_control", icon = "check-double")
                , list(name = "Normalize", key = "normalization", icon = "balance-scale")
                , list(name = "DE", key = "differential_analysis", icon = "chart-line")
                , list(name = "Summary", key = "session_summary", icon = "download")
            )
            
            render_workflow_stepper(steps, workflow_data$tab_status)
        })
        
        # --- Workflow State Observers ---
        # Enable QC tab after design matrix is complete
        shiny::observeEvent(workflow_data$tab_status$design_matrix, {
            if (workflow_data$tab_status$design_matrix == "complete") {
                logger::log_info("Metabolomics: Design matrix complete, enabling QC tab")
                workflow_data$tab_status$quality_control <- "pending"
            }
        }, ignoreNULL = TRUE)
        
        # Enable normalization tab after QC is complete
        shiny::observeEvent(workflow_data$tab_status$quality_control, {
            if (workflow_data$tab_status$quality_control == "complete") {
                logger::log_info("Metabolomics: QC complete, enabling Normalization tab")
                workflow_data$tab_status$normalization <- "pending"
            }
        }, ignoreNULL = TRUE)
        
        # Enable DE tab after normalization is complete
        shiny::observeEvent(workflow_data$tab_status$normalization, {
            if (workflow_data$tab_status$normalization == "complete") {
                logger::log_info("Metabolomics: Normalization complete, enabling DE tab")
                workflow_data$tab_status$differential_analysis <- "pending"
            }
        }, ignoreNULL = TRUE)
        
        # Enable summary tab after DE is complete
        shiny::observeEvent(workflow_data$tab_status$differential_analysis, {
            if (workflow_data$tab_status$differential_analysis == "complete") {
                logger::log_info("Metabolomics: DE complete, enabling Summary tab")
                workflow_data$tab_status$session_summary <- "pending"
            }
        }, ignoreNULL = TRUE)
        
        # Return workflow_data for external access if needed
        return(workflow_data)
    })
}


# ============================================================================
# Standalone app launcher for testing
# ============================================================================

#' @title Launch Metabolomics Workflow App
#' @description Launches a standalone Shiny app for the metabolomics workflow.
#'              Useful for testing and development.
#'
#' @param base_dir Base directory for the project. Defaults to a temp directory.
#' @param ... Additional arguments passed to shiny::shinyApp
#'
#' @return A Shiny app object
#' @export
#'
#' @examples
#' \dontrun{
#' run_metabolomics_app()
#' }
run_metabolomics_app <- function(base_dir = NULL, ...) {
    # Set up default base directory
    if (is.null(base_dir)) {
        base_dir <- file.path(tempdir(), "metabolomics_test")
        if (!dir.exists(base_dir)) {
            dir.create(base_dir, recursive = TRUE)
        }
    }
    
    # Create mock project_dirs structure matching what app_server.R creates
    project_dirs <- list(
        metabolomics = list(
            base_dir = base_dir
            , data_dir = file.path(base_dir, "data")
            , results_dir = file.path(base_dir, "results")
            , source_dir = file.path(base_dir, "scripts")
        )
    )
    
    # Create directories if they don't exist
    lapply(project_dirs$metabolomics, function(d) {
        if (!dir.exists(d)) dir.create(d, recursive = TRUE)
    })
    
    ui <- shiny::fluidPage(
        shiny::tags$head(
            shiny::tags$style(shiny::HTML("
                body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; }
                .nav-pills > li > a { border-radius: 0; }
                .badge { padding: 8px 12px; }
            "))
        )
        , mod_metabolomics_ui("metabolomics_app")
    )
    
    server <- function(input, output, session) {
        # Get volumes for shinyFiles
        volumes <- NULL
        if (requireNamespace("shinyFiles", quietly = TRUE)) {
            volumes <- shinyFiles::getVolumes()()
        }
        
        mod_metabolomics_server(
            id = "metabolomics_app"
            , project_dirs = project_dirs
            , omic_type = "metabolomics"
            , experiment_label = "Metabolomics Test"
            , volumes = volumes
        )
    }
    
    shiny::shinyApp(ui, server, ...)
}

