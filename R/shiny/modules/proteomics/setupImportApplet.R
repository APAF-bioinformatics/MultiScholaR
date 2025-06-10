#' Setup & Import Applet UI
#' 
#' UI for data import and initial setup
#' 
#' @param id Module ID
#' 
#' @importFrom shiny NS fileInput numericInput textInput actionButton uiOutput
#' @importFrom shiny fluidRow column wellPanel h3 h4 tags hr
#' @export
setupImportUi <- function(id) {
  ns <- shiny::NS(id)
  
  shiny::fluidRow(
    shiny::column(12,
      shiny::wellPanel(
        shiny::h3("Setup & Import Data"),
        
        # Data import section
        shiny::fluidRow(
          shiny::column(6,
            shiny::h4("Step 1: Import DIA-NN Report"),
            shiny::fileInput(ns("diann_report"), 
                     "Select DIA-NN report.tsv file:",
                     accept = c(".tsv", ".txt", ".tab")),
            
            shiny::h4("Step 2: Import FASTA File"),
            shiny::fileInput(ns("fasta_file"), 
                     "Select FASTA file:",
                     accept = c(".fasta", ".fa", ".faa")),
            
            shiny::h4("Step 3: Organism Information"),
            shiny::numericInput(ns("taxon_id"), 
                        "Taxonomy ID:", 
                        value = 9606,  # Human default
                        min = 1),
            shiny::textInput(ns("organism_name"), 
                     "Organism Name:", 
                     value = "Homo sapiens")
          ),
          
          shiny::column(6,
            shiny::h4("Optional: ID Mapping Files"),
            shiny::fileInput(ns("uniprot_mapping"), 
                     "UniProt ID mapping file (optional):",
                     accept = c(".tsv", ".txt", ".tab")),
            shiny::fileInput(ns("uniparc_mapping"), 
                     "UniParc ID mapping file (optional):",
                     accept = c(".tsv", ".txt", ".tab")),
            
            shiny::h4("Configuration"),
            shiny::fileInput(ns("config_file"), 
                     "Configuration file (optional):",
                     accept = c(".ini", ".cfg")),
            shiny::tags$small("If not provided, default configuration will be used")
          )
        ),
        
        shiny::hr(),
        
        # Process button and status
        shiny::fluidRow(
          shiny::column(12,
            shiny::actionButton(ns("process_data"), 
                        "Process Imported Data", 
                        class = "btn-success",
                        width = "100%"),
            shiny::br(),
            shiny::br(),
            shiny::uiOutput(ns("import_status")),
            shiny::uiOutput(ns("data_summary"))
          )
        )
      )
    )
  )
}

#' Setup & Import Applet Server
#' 
#' Server logic for data import and initial setup
#' 
#' @param id Module ID
#' @param workflow_data Reactive values object to store workflow data
#' @param experiment_paths List of paths for this experiment
#' 
#' @importFrom shiny moduleServer reactive reactiveValues observeEvent req
#' @importFrom shiny renderUI showNotification removeNotification
#' @importFrom logger log_info log_error log_warn
#' @importFrom vroom vroom
#' @importFrom ini read.ini
#' @export
setupImportServer <- function(id, workflow_data, experiment_paths) {
  shiny::moduleServer(id, function(input, output, session) {
    
    # Local reactive values
    local_data <- shiny::reactiveValues(
      processing = FALSE
    )
    
    # Process imported data
    shiny::observeEvent(input$process_data, {
      shiny::req(input$diann_report)
      shiny::req(input$fasta_file)
      
      local_data$processing <- TRUE
      
      shiny::showNotification("Processing imported data...", 
                      id = "processing_notification",
                      duration = NULL)
      
      tryCatch({
        # Read DIA-NN report
        log_info("Reading DIA-NN report from {input$diann_report$datapath}")
        data_tbl <- vroom::vroom(
          input$diann_report$datapath,
          show_col_types = FALSE,
          progress = TRUE
        )
        
        # Store in workflow data
        workflow_data$data_tbl <- data_tbl
        
        # Read FASTA file path (we don't parse it here, just store the path)
        workflow_data$fasta_file_path <- input$fasta_file$datapath
        
        # Set organism info
        workflow_data$taxon_id <- input$taxon_id
        workflow_data$organism_name <- input$organism_name
        
        # Load configuration
        if (!is.null(input$config_file)) {
          log_info("Reading configuration from {input$config_file$datapath}")
          config_list <- ini::read.ini(input$config_file$datapath)
        } else {
          log_info("Using default configuration")
          config_list <- getDefaultProteomicsConfig()
        }
        workflow_data$config_list <- config_list
        
        # Read optional mapping files
        if (!is.null(input$uniprot_mapping)) {
          log_info("Reading UniProt mapping file")
          workflow_data$uniprot_mapping <- vroom::vroom(
            input$uniprot_mapping$datapath,
            show_col_types = FALSE
          )
        }
        
        if (!is.null(input$uniparc_mapping)) {
          log_info("Reading UniParc mapping file")
          workflow_data$uniparc_mapping <- vroom::vroom(
            input$uniparc_mapping$datapath,
            show_col_types = FALSE
          )
        }
        
        # Update processing log
        workflow_data$processing_log$setup_import <- list(
          timestamp = Sys.time(),
          diann_file = input$diann_report$name,
          fasta_file = input$fasta_file$name,
          taxon_id = input$taxon_id,
          organism = input$organism_name,
          n_rows = nrow(data_tbl),
          n_runs = length(unique(data_tbl$Run)),
          n_proteins = length(unique(data_tbl$Protein.Group)),
          n_peptides = length(unique(data_tbl$Stripped.Sequence))
        )
        
        # Mark this tab as complete
        workflow_data$tab_status$setup_import <- "complete"
        
        shiny::removeNotification("processing_notification")
        shiny::showNotification("Data import successful!", type = "success")
        
        local_data$processing <- FALSE
        
      }, error = function(e) {
        log_error("Error during data import: {e$message}")
        shiny::removeNotification("processing_notification")
        shiny::showNotification(paste("Error:", e$message), type = "error")
        local_data$processing <- FALSE
      })
    })
    
    # Render import status
    output$import_status <- shiny::renderUI({
      if (local_data$processing) {
        return(shiny::tags$div(
          class = "alert alert-info",
          shiny::tags$strong("Processing..."),
          " Please wait while data is being imported and validated."
        ))
      }
      
      if (!is.null(workflow_data$data_tbl)) {
        return(shiny::tags$div(
          class = "alert alert-success",
          shiny::tags$strong("✓ Data imported successfully!")
        ))
      }
      
      return(NULL)
    })
    
    # Render data summary
    output$data_summary <- shiny::renderUI({
      shiny::req(workflow_data$data_tbl)
      
      summary_info <- workflow_data$processing_log$setup_import
      
      shiny::tags$div(
        class = "well",
        shiny::h4("Data Summary"),
        shiny::tags$ul(
          shiny::tags$li(paste("Total rows:", format(summary_info$n_rows, big.mark = ","))),
          shiny::tags$li(paste("Number of runs:", summary_info$n_runs)),
          shiny::tags$li(paste("Number of protein groups:", format(summary_info$n_proteins, big.mark = ","))),
          shiny::tags$li(paste("Number of peptides:", format(summary_info$n_peptides, big.mark = ","))),
          shiny::tags$li(paste("Organism:", summary_info$organism, "(taxon:", summary_info$taxon_id, ")"))
        )
      )
    })
    
  })
}

#' Get Default Proteomics Configuration
#' 
#' Returns default configuration settings for proteomics analysis
#' 
#' @return List with default configuration parameters
#' @export
getDefaultProteomicsConfig <- function() {
  list(
    generalParameters = list(
      min_peptides_per_protein = 2,
      min_peptides_per_sample = 2,
      q_value_threshold = 0.01,
      intensity_threshold = 0
    ),
    deAnalysisParameters = list(
      formula_string = "~ 0 + group",
      q_value_threshold = 0.05,
      log2_fc_threshold = 1
    ),
    normalizationParameters = list(
      normalisation_method = "cyclicloess"
    ),
    ruvParameters = list(
      percentage_as_neg_ctrl = 33,
      ruv_k = NULL
    )
  )
} 