#' Setup & Import Applet UI
#' 
#' UI for data import and initial setup with multi-format support
#' 
#' @param id Module ID
#' 
#' @importFrom shiny NS tagList fileInput numericInput textInput actionButton uiOutput
#' @importFrom shiny fluidRow column wellPanel h3 h4 tags hr radioButtons verbatimTextOutput helpText div br conditionalPanel icon
#' @importFrom shinyjs useShinyjs disable enable
#' @export
mod_prot_import_ui <- function(id) {
  ns <- shiny::NS(id)
  
  # Check if shinyFiles is available
  use_shiny_files <- requireNamespace("shinyFiles", quietly = TRUE)
  
  shiny::tagList(
    shinyjs::useShinyjs(),  # Enable shinyjs for dynamic UI control
    shiny::fluidRow(
    shiny::column(12,
      shiny::wellPanel(
        shiny::h3("Setup & Import Data"),
        
        # Data import section
        shiny::fluidRow(
          shiny::column(6,
            shiny::h4("Step 1: Import Searched Data File"),
            if (use_shiny_files) {
              shiny::tagList(
                shinyFiles::shinyFilesButton(
                  ns("search_results"), 
                  "Select proteomics search results file",
                  "Select File",
                  multiple = FALSE,
                  icon = shiny::icon("file")
                ),
                shiny::br(),
                shiny::verbatimTextOutput(ns("search_results_path"), placeholder = TRUE)
              )
            } else {
              shiny::fileInput(ns("search_results_standard"), 
                       "Select proteomics search results file:",
                       accept = c(".tsv", ".txt", ".tab", ".csv", ".xlsx", ".zip"))
            },
            
            # Format detection output
            shiny::uiOutput(ns("format_detection")),
            
            # Manual format override option
            shiny::radioButtons(ns("format_override"),
                        "Override detected format:",
                        choices = c("Auto-detect" = "auto",
                                  "DIA-NN" = "diann",
                                  "Spectronaut DIA" = "spectronaut",
                                  "FragPipe LFQ" = "fragpipe",
                                  "MaxQuant LFQ" = "maxquant",
                                  "Proteome Discoverer TMT" = "pd_tmt"),
                        selected = "auto",
                        inline = TRUE),
            
            shiny::h4("Step 2: Import FASTA File"),
            if (use_shiny_files) {
              shiny::tagList(
                shinyFiles::shinyFilesButton(
                  ns("fasta_file"), 
                  "Select FASTA file",
                  "Select File",
                  multiple = FALSE,
                  icon = shiny::icon("file")
                ),
                shiny::br(),
                shiny::verbatimTextOutput(ns("fasta_file_path"), placeholder = TRUE)
              )
            } else {
              shiny::fileInput(ns("fasta_file_standard"), 
                       "Select FASTA file:",
                       accept = c(".fasta", ".fa", ".faa"))
            },
            
            # Mixed species FASTA option
            shiny::checkboxInput(ns("mixed_species_fasta"),
                         "Mixed species FASTA (analyze organism distribution)",
                         value = FALSE),
            shiny::helpText("Check this if your FASTA contains proteins from multiple species (e.g., spiked-in standards, contaminants)"),
            
            shiny::h4("Step 3: Organism Information"),
            shiny::div(
              id = ns("organism_inputs_container"),
              shiny::numericInput(ns("taxon_id"), 
                          "Taxonomy ID:", 
                          value = 9606,  # Human default
                          min = 1),
              shiny::textInput(ns("organism_name"), 
                       "Organism Name:", 
                       value = "Homo sapiens")
            ),
            # Help text shown when mixed species is checked
            shiny::conditionalPanel(
              condition = paste0("input['", ns("mixed_species_fasta"), "'] == true"),
              shiny::helpText(
                shiny::icon("info-circle"),
                " Organism will be selected from the analysis popup after processing."
              ),
              style = "color: #666; font-style: italic;"
            )
          ),
          
          shiny::column(6,
            shiny::h4("Optional: ID Mapping Files"),
            if (use_shiny_files) {
              shiny::tagList(
                shinyFiles::shinyFilesButton(
                  ns("uniprot_mapping"), 
                  "UniProt ID mapping file (optional)",
                  "Select File",
                  multiple = FALSE,
                  icon = shiny::icon("file")
                ),
                shiny::verbatimTextOutput(ns("uniprot_mapping_path"), placeholder = TRUE),
                shiny::br(),
                
                shinyFiles::shinyFilesButton(
                  ns("uniparc_mapping"), 
                  "UniParc ID mapping file (optional)",
                  "Select File",
                  multiple = FALSE,
                  icon = shiny::icon("file")
                ),
                shiny::verbatimTextOutput(ns("uniparc_mapping_path"), placeholder = TRUE)
              )
            } else {
              shiny::tagList(
                shiny::fileInput(ns("uniprot_mapping_standard"), 
                         "UniProt ID mapping file (optional):",
                         accept = c(".tsv", ".txt", ".tab")),
                shiny::fileInput(ns("uniparc_mapping_standard"), 
                         "UniParc ID mapping file (optional):",
                         accept = c(".tsv", ".txt", ".tab"))
              )
            },
            
            shiny::h4("Configuration"),
            if (use_shiny_files) {
              shiny::tagList(
                shinyFiles::shinyFilesButton(
                  ns("config_file"), 
                  "Configuration file (optional)",
                  "Select File",
                  multiple = FALSE,
                  icon = shiny::icon("file")
                ),
                shiny::verbatimTextOutput(ns("config_file_path"), placeholder = TRUE)
              )
            } else {
              shiny::fileInput(ns("config_file_standard"), 
                       "Configuration file (optional):",
                       accept = c(".ini", ".cfg"))
            },
            shiny::tags$small("If not provided, default configuration will be used"),
            
            # Format-specific options will appear here
            shiny::uiOutput(ns("format_specific_options"))
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
  )
}

#' Setup & Import Applet Server
#' 
#' Server logic for data import and initial setup with multi-format support
#' 
#' @param id Module ID
#' @param workflow_data Reactive values object to store workflow data
#' @param experiment_paths List of paths for this experiment
#' 
#' @importFrom shiny moduleServer reactive reactiveValues observeEvent req
#' @importFrom shiny renderUI showNotification removeNotification showModal modalDialog removeModal updateNumericInput updateTextInput
#' @importFrom shinyjs html
#' @importFrom logger log_info log_error log_warn
#' @importFrom vroom vroom
#' @importFrom ini read.ini
#' @export
mod_prot_import_server <- function(id, workflow_data, experiment_paths, volumes = NULL) {
  message(sprintf("--- Entering mod_prot_import_server ---"))
  message(sprintf("   mod_prot_import_server Arg: id = %s", id))
  
  shiny::moduleServer(id, function(input, output, session) {
    message(sprintf("   mod_prot_import_server Step: Inside moduleServer function"))
    
    # Local reactive values
    local_data <- shiny::reactiveValues(
      processing = FALSE,
      detected_format = NULL,
      format_confidence = NULL,
      # Store file paths for shinyFiles
      search_results_file = NULL,
      fasta_file_path = NULL,
      uniprot_mapping_file = NULL,
      uniparc_mapping_file = NULL,
      config_file_path = NULL,
      # Mixed species FASTA analysis
      organism_distribution = NULL,
      organism_mapping = NULL,
      pending_import_completion = FALSE,
      waiting_for_organism_selection = FALSE
    )
    
    # Check if shinyFiles is available
    use_shiny_files <- requireNamespace("shinyFiles", quietly = TRUE)
    
    # Set up shinyFiles if available
    if (use_shiny_files) {
      # If volumes not passed in, create them
      if (is.null(volumes)) {
        volumes <- shinyFiles::getVolumes()()  # Call the function to get the actual volumes list
      }
      
      # Search results file
      shinyFiles::shinyFileChoose(
        input, 
        "search_results", 
        roots = volumes, 
        session = session,
        filetypes = c("tsv", "txt", "tab", "csv", "xlsx", "zip")
      )
      
      observeEvent(input$search_results, {
        if (!is.null(input$search_results)) {
          if (!is.integer(input$search_results)) {
            tryCatch({
              file_info <- shinyFiles::parseFilePaths(volumes, input$search_results)
              if (nrow(file_info) > 0) {
                local_data$search_results_file <- as.character(file_info$datapath[1])
                output$search_results_path <- renderText(local_data$search_results_file)
              }
            }, error = function(e) {
              message(sprintf("   mod_prot_import_server ERROR in parseFilePaths: %s", e$message))
            })
          }
        }
      })
      
      # FASTA file
      shinyFiles::shinyFileChoose(input, "fasta_file", 
                                 roots = volumes, 
                                 session = session,
                                 filetypes = c("fasta", "fa", "faa"))
      
      observeEvent(input$fasta_file, {
        if (!is.null(input$fasta_file) && !is.integer(input$fasta_file)) {
          file_info <- shinyFiles::parseFilePaths(volumes, input$fasta_file)
          if (nrow(file_info) > 0) {
            local_data$fasta_file_path <- as.character(file_info$datapath[1])
            output$fasta_file_path <- renderText(local_data$fasta_file_path)
          }
        }
      })
      
      # UniProt mapping
      shinyFiles::shinyFileChoose(input, "uniprot_mapping", 
                                 roots = volumes, 
                                 session = session,
                                 filetypes = c("tsv", "txt", "tab"))
      
      observeEvent(input$uniprot_mapping, {
        if (!is.null(input$uniprot_mapping) && !is.integer(input$uniprot_mapping)) {
          file_info <- shinyFiles::parseFilePaths(volumes, input$uniprot_mapping)
          if (nrow(file_info) > 0) {
            local_data$uniprot_mapping_file <- as.character(file_info$datapath[1])
            output$uniprot_mapping_path <- renderText(local_data$uniprot_mapping_file)
          }
        }
      })
      
      # UniParc mapping
      shinyFiles::shinyFileChoose(input, "uniparc_mapping", 
                                 roots = volumes, 
                                 session = session,
                                 filetypes = c("tsv", "txt", "tab"))
      
      observeEvent(input$uniparc_mapping, {
        if (!is.null(input$uniparc_mapping) && !is.integer(input$uniparc_mapping)) {
          file_info <- shinyFiles::parseFilePaths(volumes, input$uniparc_mapping)
          if (nrow(file_info) > 0) {
            local_data$uniparc_mapping_file <- as.character(file_info$datapath[1])
            output$uniparc_mapping_path <- renderText(local_data$uniparc_mapping_file)
          }
        }
      })
      
      # Config file
      shinyFiles::shinyFileChoose(input, "config_file", 
                                 roots = volumes, 
                                 session = session,
                                 filetypes = c("ini", "cfg"))
      
      observeEvent(input$config_file, {
        if (!is.null(input$config_file) && !is.integer(input$config_file)) {
          file_info <- shinyFiles::parseFilePaths(volumes, input$config_file)
          if (nrow(file_info) > 0) {
            local_data$config_file_path <- as.character(file_info$datapath[1])
            output$config_file_path <- renderText(local_data$config_file_path)
          }
        }
      })
    }
    
    # --- Observer: Toggle Step 3 inputs when mixed species is checked ---
    shiny::observeEvent(input$mixed_species_fasta, {
      if (isTRUE(input$mixed_species_fasta)) {
        # Disable the organism inputs when mixed species is checked
        shinyjs::disable("taxon_id")
        shinyjs::disable("organism_name")
        message("[mod_prot_import] Mixed species FASTA enabled - organism inputs disabled")
      } else {
        # Re-enable the organism inputs when unchecked
        shinyjs::enable("taxon_id")
        shinyjs::enable("organism_name")
        message("[mod_prot_import] Mixed species FASTA disabled - organism inputs enabled")
      }
    }, ignoreInit = FALSE)
    
    # Get the search results file path (works for both shinyFiles and standard)
    get_search_results_path <- reactive({
      if (use_shiny_files) {
        local_data$search_results_file
      } else {
        if (!is.null(input$search_results_standard)) {
          input$search_results_standard$datapath
        } else {
          NULL
        }
      }
    })
    
    # Get the FASTA file path
    get_fasta_path <- reactive({
      if (use_shiny_files) {
        local_data$fasta_file_path
      } else {
        if (!is.null(input$fasta_file_standard)) {
          input$fasta_file_standard$datapath
        } else {
          NULL
        }
      }
    })
    
    # Detect file format when file is uploaded
    shiny::observeEvent(get_search_results_path(), {
      file_path <- get_search_results_path()
      shiny::req(file_path)
      
      tryCatch({
        headers <- NULL
        
        # NEW LOGIC: Handle ZIP files differently
        if (tolower(tools::file_ext(file_path)) == "zip") {
          # List files inside the zip without extracting all
          zip_contents <- unzip(file_path, list = TRUE)
          # Find the first data file (xlsx, tsv, csv)
          first_data_file <- zip_contents$Name[grep("\\.(xlsx|tsv|csv)$", zip_contents$Name, ignore.case = TRUE)[1]]
          
          if (is.na(first_data_file)) {
            stop("No data files (.xlsx, .tsv, .csv) found inside the ZIP archive.")
          }
          
          # Read just the header row from the first data file within the zip
          if (tolower(tools::file_ext(first_data_file)) == "xlsx") {
             if (!requireNamespace("readxl", quietly = TRUE)) stop("Package 'readxl' needed for .xlsx files.")
             # readxl needs to read the file from a temp extraction
             temp_dir <- tempfile()
             dir.create(temp_dir)
             unzip(file_path, files = first_data_file, exdir = temp_dir, junkpaths = TRUE)
             unzipped_file_path <- file.path(temp_dir, basename(first_data_file))
             headers <- names(readxl::read_excel(unzipped_file_path, n_max = 0))
             unlink(temp_dir, recursive = TRUE)
          } else {
            con <- unz(file_path, first_data_file)
            headers <- strsplit(readLines(con, n = 1), "\t|,")[[1]]
            close(con)
          }
          
        } else {
          # Original logic for non-zip files
          preview_lines <- readLines(file_path, n = 10)
          headers <- strsplit(preview_lines[1], "\t|,")[[1]]
        }
        
        if (is.null(headers)) {
            stop("Could not read headers from the provided file.")
        }
        
        # Get filename for detection
        filename <- if (use_shiny_files) {
          basename(file_path)
        } else {
          input$search_results_standard$name
        }
        
        # Detect format based on column headers and filename
        format_info <- detectProteomicsFormat(
          headers = headers,
          filename = filename
        )
        
        local_data$detected_format <- format_info$format
        local_data$format_confidence <- format_info$confidence
        
        log_info(sprintf("Detected format: %s (confidence: %s)", format_info$format, format_info$confidence))
        
      }, error = function(e) {
        log_error(paste("Error detecting file format:", e$message))
        local_data$detected_format <- "unknown"
        local_data$format_confidence <- 0
      })
    })
    
    # Get active format (either detected or overridden)
    active_format <- shiny::reactive({
      if (input$format_override == "auto") {
        return(local_data$detected_format)
      } else {
        return(input$format_override)
      }
    })
    
    # Render format detection output
    output$format_detection <- shiny::renderUI({
      shiny::req(local_data$detected_format)
      
      confidence_color <- if (local_data$format_confidence >= 0.8) {
        "success"
      } else if (local_data$format_confidence >= 0.5) {
        "warning"
      } else {
        "danger"
      }
      
      format_display <- switch(local_data$detected_format,
        "diann" = "DIA-NN",
        "spectronaut" = "Spectronaut DIA",
        "fragpipe" = "FragPipe LFQ",
        "maxquant" = "MaxQuant LFQ",
        "pd_tmt" = "Proteome Discoverer TMT",
        "unknown" = "Unknown format"
      )
      
      shiny::tags$div(
        class = paste("alert", paste0("alert-", confidence_color)),
        shiny::tags$strong("Detected format: "),
        format_display,
        shiny::tags$br(),
        shiny::tags$small(
          paste0("Confidence: ", round(local_data$format_confidence * 100), "%")
        )
      )
    })
    
    # Render format-specific options
    output$format_specific_options <- shiny::renderUI({
      format <- active_format()
      shiny::req(format)
      
      ns <- session$ns
      
      switch(format,
        "diann" = shiny::tagList(
          shiny::h5("DIA-NN Specific Options"),
          shiny::checkboxInput(ns("diann_use_precursor_norm"),
                       "Use Precursor.Normalised values",
                       value = TRUE)
        ),
        "spectronaut" = shiny::tagList(
          shiny::h5("Spectronaut Specific Options"),
          shiny::radioButtons(ns("spectronaut_quantity"),
                      "Quantity type:",
                      choices = c("PG.Quantity" = "pg",
                                "PEP.Quantity" = "pep"),
                      selected = "pg")
        ),
        "fragpipe" = shiny::tagList(
          shiny::h5("FragPipe LFQ Options"),
          shiny::checkboxInput(ns("fragpipe_use_maxlfq"),
                       "Use MaxLFQ intensities",
                       value = FALSE)
        ),
        "maxquant" = shiny::tagList(
          shiny::h5("MaxQuant LFQ Options"),
          shiny::checkboxInput(ns("maxquant_use_lfq"),
                       "Use LFQ intensities",
                       value = TRUE),
          shiny::checkboxInput(ns("maxquant_filter_contaminants"),
                       "Filter contaminants",
                       value = TRUE)
        ),
        shiny::tags$div(
          class = "alert alert-warning",
          "Format-specific options not available"
        )
      )
    })
    
    # Helper function to update processing modal status
    update_processing_status <- function(status_text) {
      message(sprintf("[mod_prot_import] %s", status_text))
      shinyjs::html("processing_status_text", status_text)
    }
    
    # Process imported data
    shiny::observeEvent(input$process_data, {
      # Get file paths
      search_results_path <- get_search_results_path()
      fasta_path <- get_fasta_path()
      
      shiny::req(search_results_path)
      shiny::req(fasta_path)
      
      local_data$processing <- TRUE
      format <- active_format()
      
      ns <- session$ns
      
      message("========================================")
      message("[mod_prot_import] Starting data import process")
      message("========================================")
      message(sprintf("[mod_prot_import] Search results: %s", search_results_path))
      message(sprintf("[mod_prot_import] FASTA file: %s", fasta_path))
      message(sprintf("[mod_prot_import] Detected format: %s", format))
      
      # Show processing modal with spinner
      shiny::showModal(shiny::modalDialog(
        title = shiny::tagList(
          shiny::icon("cog", class = "fa-spin")
          , " Processing Data"
        )
        , shiny::tags$div(
            style = "text-align: center; padding: 20px;"
            , shiny::tags$div(
                style = "font-size: 48px; color: #3c8dbc; margin-bottom: 20px;"
                , shiny::icon("spinner", class = "fa-spin fa-pulse")
              )
            , shiny::tags$div(
                id = ns("processing_status_text")
                , style = "font-size: 16px; color: #333;"
                , "Initializing..."
              )
            , shiny::tags$div(
                style = "margin-top: 15px; font-size: 12px; color: #999;"
                , "Please wait while your data is being processed."
              )
          )
        , footer = NULL
        , size = "s"
        , easyClose = FALSE
      ))
      
      tryCatch({
        # Read data based on format
        update_processing_status(sprintf("Reading %s data...", toupper(format)))
        log_info(sprintf("Reading %s data from %s", format, search_results_path))
        
        # Import data with error handling
        data_import_result <- tryCatch({
          switch(format,
            "diann" = importDIANNData(
              filepath = search_results_path,
              use_precursor_norm = input$diann_use_precursor_norm %||% TRUE
            ),
            "spectronaut" = importSpectronautData(
              filepath = search_results_path,
              quantity_type = input$spectronaut_quantity %||% "pg"
            ),
            "fragpipe" = importFragPipeData(
              filepath = search_results_path,
              use_maxlfq = input$fragpipe_use_maxlfq %||% TRUE
            ),
            "maxquant" = importMaxQuantData(
              filepath = search_results_path,
              use_lfq = input$maxquant_use_lfq %||% TRUE,
              filter_contaminants = input$maxquant_filter_contaminants %||% TRUE
            ),
            "pd_tmt" = importProteomeDiscovererTMTData(
              filepath = search_results_path
            ),
            stop("Unsupported format: ", format)
          )
        }, error = function(e) {
          log_error(paste("Error in data import function:", e$message))
          stop("Failed to import data: ", e$message)
        })
        
        # Validate import result
        if (is.null(data_import_result) || is.null(data_import_result$data)) {
          stop("Data import returned NULL or empty data")
        }
        
        log_info(sprintf("Data imported successfully. Rows: %d", nrow(data_import_result$data)))
        message(sprintf("[mod_prot_import] Data imported: %d rows", nrow(data_import_result$data)))
        
        update_processing_status("Storing imported data...")
        
        # Store in workflow data
        workflow_data$data_tbl <- data_import_result$data
        workflow_data$data_format <- format
        workflow_data$data_type <- data_import_result$data_type  # "peptide" or "protein"
        workflow_data$column_mapping <- data_import_result$column_mapping
        
        # Initialize data_cln as a copy of data_tbl (will be modified by design matrix)
        workflow_data$data_cln <- data_import_result$data
        
        # Set workflow type based on detected format
        if (format == "pd_tmt") {
          workflow_data$state_manager$setWorkflowType("TMT")
        } else if (format == "diann") {
          workflow_data$state_manager$setWorkflowType("DIA")
        } else {
          workflow_data$state_manager$setWorkflowType("LFQ")
        }
        
        # Read FASTA file path (we don't parse it here, just store the path)
        workflow_data$fasta_file_path <- fasta_path
        
        # Process FASTA file to create aa_seq_tbl_final
        update_processing_status("Processing FASTA file...")
        log_info("Processing FASTA file...")
        
        # Get mapping files if provided
        uniprot_mapping <- if (!is.null(workflow_data$uniprot_mapping)) {
          workflow_data$uniprot_mapping
        } else {
          NULL
        }
        
        uniparc_mapping <- if (!is.null(workflow_data$uniparc_mapping)) {
          workflow_data$uniparc_mapping
        } else {
          NULL
        }
        
        # Create fasta meta file path in experiment directory
        log_info("Checking experiment_paths...")
        if (!is.null(experiment_paths) && !is.null(experiment_paths$results_dir)) {
          fasta_meta_file <- file.path(
            experiment_paths$results_dir, 
            "cache", 
            "aa_seq_tbl.RDS"
          )
          cache_dir <- file.path(experiment_paths$results_dir, "cache")
        } else {
          log_warn("experiment_paths not properly initialized, using temp directory for cache")
          cache_dir <- file.path(tempdir(), "proteomics_cache")
          fasta_meta_file <- file.path(cache_dir, "aa_seq_tbl.RDS")
        }
        
        # Create cache directory if it doesn't exist
        if (!dir.exists(cache_dir)) {
          dir.create(cache_dir, recursive = TRUE)
        }
        
        # Process FASTA file
        tryCatch({
          fasta_result <- processFastaFile(
            fasta_file_path = fasta_path,
            uniprot_search_results = uniprot_mapping,
            uniparc_search_results = uniparc_mapping,
            fasta_meta_file = fasta_meta_file,
            organism_name = input$organism_name
          )
          
          # Extract data and metadata from result
          aa_seq_tbl_final <- fasta_result$aa_seq_tbl_final
          fasta_metadata <- fasta_result$fasta_metadata
          
          workflow_data$aa_seq_tbl_final <- aa_seq_tbl_final
          workflow_data$fasta_metadata <- fasta_metadata
          
          # ✅ FIXED: Store in global environment for chooseBestProteinAccession to find
          assign("aa_seq_tbl_final", aa_seq_tbl_final, envir = .GlobalEnv)
          
          # ✅ FIXED: Save to /scripts directory for session persistence
          if (!is.null(experiment_paths) && !is.null(experiment_paths$source_dir)) {
            scripts_aa_seq_path <- file.path(experiment_paths$source_dir, "aa_seq_tbl_final.RDS")
            saveRDS(aa_seq_tbl_final, scripts_aa_seq_path)
            log_info(sprintf("Saved aa_seq_tbl_final to scripts directory: %s", scripts_aa_seq_path))
            
            # Save FASTA metadata
            fasta_metadata_path <- file.path(experiment_paths$source_dir, "fasta_metadata.RDS")
            saveRDS(fasta_metadata, fasta_metadata_path)
            log_info(sprintf("Saved FASTA metadata to scripts directory: %s", fasta_metadata_path))
            log_info(sprintf("FASTA Format: %s, Sequences: %d", fasta_metadata$fasta_format, fasta_metadata$num_sequences))
          }
          
          log_info(sprintf("FASTA file processed successfully. Found %d sequences", nrow(aa_seq_tbl_final)))
          
        }, error = function(e) {
          log_warn(paste("Error processing FASTA file:", e$message))
          log_warn("Continuing without FASTA processing - protein ID conversion will be skipped")
          workflow_data$aa_seq_tbl_final <- NULL
          workflow_data$fasta_metadata <- NULL
        })
        
        # --- Mixed Species FASTA Analysis ---
        # If user checked mixed species option, analyze organism distribution
        if (isTRUE(input$mixed_species_fasta) && !is.null(workflow_data$aa_seq_tbl_final)) {
          tryCatch({
            update_processing_status("Analyzing organism distribution...")
            message("[mod_prot_import] Mixed species FASTA - analyzing organism distribution")
            log_info("Mixed species FASTA detected - analyzing organism distribution...")
            
            # Extract organism info from FASTA
            organism_mapping <- extractOrganismsFromFasta(fasta_path)
            local_data$organism_mapping <- organism_mapping
            
            # Get protein IDs from imported data
            protein_col <- data_import_result$column_mapping$protein_col
            if (!is.null(protein_col) && protein_col %in% names(data_import_result$data)) {
              protein_ids <- unique(data_import_result$data[[protein_col]])
              
              # Analyze distribution
              organism_dist <- analyzeOrganismDistribution(protein_ids, organism_mapping)
              local_data$organism_distribution <- organism_dist
              workflow_data$organism_distribution <- organism_dist
              
              log_info(sprintf("Found %d organisms in data", nrow(organism_dist)))
              
              # Show modal for organism selection
              ns <- session$ns
              
              # Build radio button choices from distribution
              # Format: "Organism Name (Taxon: XXXXX) - XX.X%"
              valid_organisms <- organism_dist |>
                dplyr::filter(!is.na(taxon_id) & organism_name != "[Unmatched/Unknown]")
              
              if (nrow(valid_organisms) > 0) {
                # Set flag to indicate we're waiting for organism selection
                local_data$waiting_for_organism_selection <- TRUE
                
                choices <- stats::setNames(
                  as.character(valid_organisms$taxon_id)
                  , paste0(
                      valid_organisms$organism_name
                      , " (Taxon: ", valid_organisms$taxon_id, ") - "
                      , valid_organisms$protein_count, " proteins ("
                      , valid_organisms$percentage, "%)"
                    )
                )
                
                # Pre-select the organism with highest count
                selected_taxon <- as.character(valid_organisms$taxon_id[1])
                
                shiny::showModal(shiny::modalDialog(
                  title = "Select Primary Organism"
                  , size = "l"
                  , shiny::tags$div(
                      shiny::tags$p(
                        shiny::tags$strong("Multiple organisms detected in your FASTA database.")
                        , " The distribution of matched proteins is shown below."
                        , " Please select the primary organism for this analysis."
                      )
                      , shiny::tags$hr()
                      , shiny::tags$h5("Organism Distribution")
                      , DT::DTOutput(ns("organism_dist_table"))
                      , shiny::tags$hr()
                      , shiny::radioButtons(
                          ns("selected_organism")
                          , "Select primary organism:"
                          , choices = choices
                          , selected = selected_taxon
                        )
                      , shiny::checkboxInput(
                          ns("filter_to_organism")
                          , "Filter data to keep only proteins from selected organism"
                          , value = FALSE
                        )
                      , shiny::helpText(
                          "If checked, proteins from other organisms will be removed from the analysis."
                        )
                    )
                  , footer = shiny::tagList(
                      shiny::modalButton("Cancel")
                      , shiny::actionButton(ns("confirm_organism"), "Confirm Selection", class = "btn-primary")
                    )
                ))
              } else {
                log_warn("No valid organisms with taxon IDs found in FASTA")
                shiny::showNotification(
                  "Could not extract organism information from FASTA headers. Using default organism settings."
                  , type = "warning"
                )
              }
            }
          }, error = function(e) {
            log_warn(paste("Error analyzing organism distribution:", e$message))
            shiny::showNotification(
              paste("Could not analyze organism distribution:", e$message)
              , type = "warning"
            )
          })
        }
        
        # Load configuration
        update_processing_status("Loading configuration...")
        message("[mod_prot_import] Loading configuration")
        
        config_path <- if (use_shiny_files) {
          local_data$config_file_path
        } else {
          if (!is.null(input$config_file_standard)) input$config_file_standard$datapath else NULL
        }
        
        if (!is.null(config_path)) {
          log_info(paste("Reading configuration from", config_path))
          config_list <- readConfigFile(file = config_path)
        } else {
          log_info("Using default configuration")
          default_config_path <- file.path(experiment_paths$source_dir, "config.ini")
          if (file.exists(default_config_path)) {
            config_list <- readConfigFile(file = default_config_path)
          } else {
            log_info("config.ini not found in project. Retrieving default config.")
            tryCatch({
              # 1. Try to get from package installation first (system.file)
              pkg_config <- system.file("config", "config.ini", package = "MultiScholaR")
              
              if (file.exists(pkg_config) && pkg_config != "") {
                log_info(paste("Found default config.ini in package:", pkg_config))
                file.copy(pkg_config, default_config_path)
                log_info(paste("Default config.ini copied to:", default_config_path))
                shiny::showNotification("Copied default config.ini to scripts directory.", type = "message")
                config_list <- readConfigFile(file = default_config_path)
              } else {
                # 2. Fallback to GitHub download
                log_info("Default config.ini not found in package, downloading from GitHub...")
                default_config_url <- "https://raw.githubusercontent.com/APAF-bioinformatics/MultiScholaR/main/inst/config/config.ini"
                download.file(default_config_url, destfile = default_config_path, quiet = TRUE)
                log_info(paste("Default config.ini downloaded to:", default_config_path))
                shiny::showNotification("Downloaded default config.ini to scripts directory.", type = "message")
                config_list <- readConfigFile(file = default_config_path)
              }
            }, error = function(e_download) {
              msg <- paste("Failed to retrieve default config.ini:", e_download$message)
              log_error(msg)
              shiny::showNotification(msg, type = "warning", duration = 10)
              log_info("Using minimal fallback configuration")
              config_list <- getDefaultProteomicsConfig()
            })
          }
        }
        workflow_data$config_list <- config_list

        # Create global config_list for compatibility
        assign("config_list", config_list, envir = .GlobalEnv)
        log_info("Created global config_list for compatibility")
        
        # Read optional mapping files
        uniprot_path <- if (use_shiny_files) {
          local_data$uniprot_mapping_file
        } else {
          if (!is.null(input$uniprot_mapping_standard)) input$uniprot_mapping_standard$datapath else NULL
        }
        
        uniparc_path <- if (use_shiny_files) {
          local_data$uniparc_mapping_file
        } else {
          if (!is.null(input$uniparc_mapping_standard)) input$uniparc_mapping_standard$datapath else NULL
        }
        
        if (!is.null(uniprot_path)) {
          log_info(paste("Reading UniProt mapping from", uniprot_path))
          uniprot_mapping <- tryCatch({
            vroom::vroom(uniprot_path, show_col_types = FALSE)
          }, error = function(e) {
            log_error(paste("Failed to read UniProt mapping file:", e$message))
            stop("Failed to read UniProt mapping file: ", e$message)
          })
          workflow_data$uniprot_mapping <- uniprot_mapping
          log_info(sprintf("Read %d rows from UniProt mapping file", nrow(uniprot_mapping)))
        }
        
        if (!is.null(uniparc_path)) {
          log_info(paste("Reading UniParc mapping from", uniparc_path))
          uniparc_mapping <- tryCatch({
            vroom::vroom(uniparc_path, show_col_types = FALSE)
          }, error = function(e) {
            log_error(paste("Failed to read UniParc mapping file:", e$message))
            stop("Failed to read UniParc mapping file: ", e$message)
          })
          workflow_data$uniparc_mapping <- uniparc_mapping
          log_info(sprintf("Read %d rows from UniParc mapping file", nrow(uniparc_mapping)))
        }
        
        # Get filename for logging
        search_filename <- if (use_shiny_files) {
          basename(search_results_path)
        } else {
          input$search_results_standard$name
        }
        
        fasta_filename <- if (use_shiny_files) {
          basename(fasta_path)
        } else {
          input$fasta_file_standard$name
        }
        
        # Update processing log
        log_info("Creating processing log...")
        
        # Safely get column names with defaults
        run_col <- data_import_result$column_mapping$run_col
        protein_col <- data_import_result$column_mapping$protein_col
        peptide_col <- data_import_result$column_mapping$peptide_col
        
        # Calculate statistics with error handling
        n_runs <- if (!is.null(run_col) && run_col %in% names(data_import_result$data)) {
          length(unique(data_import_result$data[[run_col]]))
        } else {
          NA
        }
        
        n_proteins <- if (!is.null(protein_col) && protein_col %in% names(data_import_result$data)) {
          length(unique(data_import_result$data[[protein_col]]))
        } else {
          NA
        }
        
        n_peptides <- if (!is.null(peptide_col) && peptide_col %in% names(data_import_result$data)) {
          length(unique(data_import_result$data[[peptide_col]]))
        } else {
          NA
        }
        
        # Store taxon_id and organism at top level for use by other modules
        workflow_data$taxon_id <- input$taxon_id
        workflow_data$organism_name <- input$organism_name
        
        workflow_data$processing_log$setup_import <- list(
          timestamp = Sys.time(),
          search_file = search_filename,
          detected_format = format,
          data_type = data_import_result$data_type,
          fasta_file = fasta_filename,
          taxon_id = input$taxon_id,
          organism = input$organism_name,
          n_rows = nrow(data_import_result$data),
          n_runs = n_runs,
          n_proteins = n_proteins,
          n_peptides = n_peptides
        )
        
        log_info("Processing log created successfully")
        
        update_processing_status("Finalizing import...")
        message("[mod_prot_import] Finalizing import")
        
        # Mark this tab as complete
        workflow_data$tab_status$setup_import <- "complete"
        
        message("========================================")
        message("[mod_prot_import] Data import completed successfully!")
        message(sprintf("[mod_prot_import] Rows: %d, Proteins: %s, Format: %s", 
                       nrow(data_import_result$data), n_proteins, format))
        message("========================================")
        
        # Only close modal and show success if we're NOT waiting for organism selection
        # (The organism selection modal replaces the processing modal)
        if (!isTRUE(local_data$waiting_for_organism_selection)) {
          shiny::removeModal()
          shiny::showNotification("Data import successful!", type = "message")
        } else {
          message("[mod_prot_import] Waiting for organism selection from user...")
        }
        
        local_data$processing <- FALSE
        
      }, error = function(e) {
        message("========================================")
        message(sprintf("[mod_prot_import] ERROR: %s", e$message))
        message("========================================")
        
        log_error(paste("Error during data import:", e$message))
        shiny::removeModal()
        shiny::showNotification(paste("Error:", e$message), type = "error", duration = 10)
        local_data$processing <- FALSE
        
        # Clean up workflow_data on error to prevent inconsistent state
        workflow_data$data_tbl <- NULL
        workflow_data$data_format <- NULL
        workflow_data$data_type <- NULL
        workflow_data$column_mapping <- NULL
        workflow_data$data_cln <- NULL
        workflow_data$fasta_file_path <- NULL
        workflow_data$aa_seq_tbl_final <- NULL
        workflow_data$config_list <- NULL
        workflow_data$processing_log$setup_import <- NULL
        workflow_data$tab_status$setup_import <- "incomplete"
      })
    })
    
    # --- Organism Distribution Table Output ---
    output$organism_dist_table <- DT::renderDT({
      shiny::req(local_data$organism_distribution)
      
      dist_data <- local_data$organism_distribution |>
        dplyr::select(
          `Organism` = organism_name
          , `Taxon ID` = taxon_id
          , `Proteins` = protein_count
          , `%` = percentage
        )
      
      DT::datatable(
        dist_data
        , options = list(
            pageLength = 10
            , searching = FALSE
            , lengthChange = FALSE
            , ordering = FALSE
            , info = FALSE
          )
        , rownames = FALSE
        , class = "compact stripe"
      ) |>
        DT::formatStyle(
          columns = c("Organism", "Taxon ID", "Proteins", "%")
          , fontSize = "12px"
        )
    })
    
    # --- Observer for Organism Selection Confirmation ---
    shiny::observeEvent(input$confirm_organism, {
      shiny::req(input$selected_organism)
      shiny::req(local_data$organism_distribution)
      
      selected_taxon <- as.integer(input$selected_organism)
      
      # Find the selected organism info
      selected_org <- local_data$organism_distribution |>
        dplyr::filter(taxon_id == selected_taxon)
      
      if (nrow(selected_org) > 0) {
        new_organism_name <- selected_org$organism_name[1]
        new_taxon_id <- selected_org$taxon_id[1]
        
        log_info(sprintf(
          "User selected organism: %s (Taxon: %d)"
          , new_organism_name
          , new_taxon_id
        ))
        
        # Update workflow_data with selected organism
        workflow_data$taxon_id <- new_taxon_id
        workflow_data$organism_name <- new_organism_name
        
        # Update the input fields to reflect selection
        shiny::updateNumericInput(session, "taxon_id", value = new_taxon_id)
        shiny::updateTextInput(session, "organism_name", value = new_organism_name)
        
        # Update processing log if it exists
        if (!is.null(workflow_data$processing_log$setup_import)) {
          workflow_data$processing_log$setup_import$taxon_id <- new_taxon_id
          workflow_data$processing_log$setup_import$organism <- new_organism_name
          workflow_data$processing_log$setup_import$mixed_species_selection <- list(
            selected_organism = new_organism_name
            , selected_taxon_id = new_taxon_id
            , filter_applied = input$filter_to_organism
            , organism_distribution = local_data$organism_distribution
          )
        }
        
        # Filter data if requested
        if (isTRUE(input$filter_to_organism)) {
          log_info("Filtering data to selected organism only...")
          
          # Get protein accessions for the selected organism from FASTA
          organism_proteins <- local_data$organism_mapping |>
            dplyr::filter(taxon_id == selected_taxon) |>
            dplyr::pull(uniprot_acc)
          
          # Clean accessions for matching (remove isoform numbers)
          clean_acc <- function(acc) sub("-\\d+$", "", acc)
          organism_proteins_clean <- unique(clean_acc(organism_proteins))
          
          # Filter the data tables
          protein_col <- workflow_data$column_mapping$protein_col
          
          if (!is.null(protein_col) && !is.null(workflow_data$data_tbl)) {
            original_count <- nrow(workflow_data$data_tbl)
            
            # Filter function that checks if any protein in a group matches
            filter_by_organism <- function(protein_ids) {
              # Split protein groups
              ids <- unlist(strsplit(as.character(protein_ids), ";"))
              ids_clean <- clean_acc(trimws(ids))
              any(ids_clean %in% organism_proteins_clean)
            }
            
            # Apply filter to data_tbl
            keep_rows <- sapply(
              workflow_data$data_tbl[[protein_col]]
              , filter_by_organism
            )
            workflow_data$data_tbl <- workflow_data$data_tbl[keep_rows, ]
            
            # Apply filter to data_cln if it exists
            if (!is.null(workflow_data$data_cln)) {
              keep_rows_cln <- sapply(
                workflow_data$data_cln[[protein_col]]
                , filter_by_organism
              )
              workflow_data$data_cln <- workflow_data$data_cln[keep_rows_cln, ]
            }
            
            filtered_count <- nrow(workflow_data$data_tbl)
            removed_count <- original_count - filtered_count
            
            log_info(sprintf(
              "Filtered data: kept %d rows, removed %d rows (%.1f%%)"
              , filtered_count
              , removed_count
              , (removed_count / original_count) * 100
            ))
            
            shiny::showNotification(
              sprintf(
                "Filtered to %s: kept %d rows, removed %d rows"
                , new_organism_name
                , filtered_count
                , removed_count
              )
              , type = "message"
              , duration = 5
            )
          }
        }
        
        # Reset the waiting flag
        local_data$waiting_for_organism_selection <- FALSE
        
        # Close the modal
        shiny::removeModal()
        
        shiny::showNotification(
          sprintf("Primary organism set to: %s (Taxon: %d)", new_organism_name, new_taxon_id)
          , type = "message"
        )
        
        # Show the import success notification that was deferred
        shiny::showNotification("Data import successful!", type = "message")
        
      } else {
        shiny::showNotification("Could not find selected organism", type = "error")
      }
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
          shiny::tags$strong("✓ Data imported successfully!"),
          shiny::tags$br(),
          paste("Format:", workflow_data$data_format, "| Type:", workflow_data$data_type)
        ))
      }
      
      return(NULL)
    })
    
    # Render data summary
    output$data_summary <- shiny::renderUI({
      shiny::req(workflow_data$data_tbl)
      shiny::req(workflow_data$processing_log$setup_import)
      
      summary_info <- workflow_data$processing_log$setup_import
      
      # Build the list items, only including peptides if available
      list_items <- list(
        shiny::tags$li(paste("Data format:", toupper(summary_info$detected_format))),
        shiny::tags$li(paste("Data type:", summary_info$data_type)),
        shiny::tags$li(paste("Total rows:", format(summary_info$n_rows, big.mark = ","))),
        shiny::tags$li(paste("Number of runs:", summary_info$n_runs)),
        shiny::tags$li(paste("Number of protein groups:", format(summary_info$n_proteins, big.mark = ",")))
      )
      
      # Add peptides info only if available
      if (!is.null(summary_info$n_peptides) && !is.na(summary_info$n_peptides)) {
        list_items <- append(list_items, list(
          shiny::tags$li(paste("Number of peptides:", format(summary_info$n_peptides, big.mark = ",")))
        ))
      }
      
      # Add organism info
      list_items <- append(list_items, list(
        shiny::tags$li(paste("Organism:", summary_info$organism, "(taxon:", summary_info$taxon_id, ")"))
      ))
      
      shiny::tags$div(
        class = "well",
        shiny::h4("Data Summary"),
        do.call(shiny::tags$ul, list_items)
      )
    })
    
  })
}

#' Detect Proteomics Data Format
#' 
#' Detects the format of proteomics search results based on headers and filename
#' 
#' @param headers Character vector of column headers
#' @param filename Name of the file
#' @param preview_lines First few lines of the file
#' 
#' @return List with format and confidence score
#' @export
detectProteomicsFormat <- function(headers, filename, preview_lines = NULL) {
  # Convert to lowercase for comparison
  headers_lower <- tolower(headers)
  filename_lower <- tolower(filename)
  
  # DIA-NN detection
  diann_score <- 0
  diann_markers <- c("protein.group", "protein.ids", "protein.names", 
                     "precursor.id", "modified.sequence", "stripped.sequence",
                     "precursor.charge", "q.value", "pg.q.value", "run")
  diann_found <- sum(diann_markers %in% headers_lower)
  diann_score <- diann_found / length(diann_markers)
  
  # Spectronaut detection
  spectronaut_score <- 0
  spectronaut_markers <- c("pg.proteingroups", "pg.proteinaccessions", "pg.genes",
                          "eg.precursorid", "eg.modifiedsequence", "fg.charge",
                          "eg.qvalue", "pg.qvalue", "r.filename", "pg.quantity")
  spectronaut_found <- sum(spectronaut_markers %in% headers_lower)
  spectronaut_score <- spectronaut_found / length(spectronaut_markers)
  
  # FragPipe detection
  fragpipe_score <- 0
  
  # Check for key markers (flexible matching)
  has_protein_id <- any(grepl("^protein\\s+id$|^protein\\.id$|^protein_id$", headers_lower))
  has_protein <- "protein" %in% headers_lower
  has_gene <- "gene" %in% headers_lower
  has_description <- any(grepl("description", headers_lower))
  has_spectral_count <- any(grepl("spectral.*count", headers_lower))
  
  # Strong indicator: columns ending with "Intensity" (case-insensitive)
  intensity_cols <- sum(grepl("intensity$", headers_lower))
  has_intensity_cols <- intensity_cols > 0
  
  # Check for MaxLFQ Intensity columns specifically
  has_maxlfq <- any(grepl("maxlfq.*intensity$", headers_lower))
  
  # Count basic markers found
  basic_markers_found <- sum(c(has_protein_id, has_protein, has_gene, has_description, has_spectral_count))
  
  # Weighted scoring similar to TMT detection
  # Basic markers worth 40% (8% each)
  basic_score <- (basic_markers_found / 5) * 0.4
  
  # Intensity columns presence is worth 40% (strong indicator)
  intensity_score <- if (has_intensity_cols) {
    # Bonus if multiple intensity columns (typical of FragPipe)
    min(0.4, 0.3 + (min(intensity_cols, 10) / 10) * 0.1)
  } else {
    0
  }
  
  # MaxLFQ presence is worth 20% (very specific to FragPipe)
  maxlfq_score <- if (has_maxlfq) 0.2 else 0
  
  fragpipe_score <- basic_score + intensity_score + maxlfq_score
  
  # Filename bonus (cap at 1.0)
  if (grepl("fragpipe|msfragger", filename_lower)) {
    fragpipe_score <- min(1.0, fragpipe_score + 0.1)
  }
  
  # MaxQuant detection
  maxquant_score <- 0
  maxquant_markers <- c("proteins", "majority.protein.ids", "protein.names",
                       "gene.names", "peptide.counts", "unique.peptides",
                       "intensity", "lfq.intensity", "ms.ms.count")
  maxquant_found <- sum(maxquant_markers %in% headers_lower)
  maxquant_score <- maxquant_found / length(maxquant_markers)
  if (grepl("proteingroups", filename_lower)) maxquant_score <- maxquant_score + 0.3
  
  # PD-TMT detection
  pd_tmt_score <- 0
  pd_tmt_markers <- c("protein fdr confidence", "master", "accession", "exp. q-value", "sum pep score")
  pd_tmt_found <- sum(pd_tmt_markers %in% headers_lower)
  abundance_found <- any(grepl("^abundance:", headers_lower))
  
  # Weighted scoring to prevent exceeding 100%
  # Header markers are worth 60% total (12% each)
  header_score <- (pd_tmt_found / length(pd_tmt_markers)) * 0.6
  # Abundance column presence is worth 40%
  abundance_score <- if (abundance_found) 0.4 else 0
  
  pd_tmt_score <- header_score + abundance_score
  
  # Determine best match
  scores <- c(diann = diann_score, 
              spectronaut = spectronaut_score,
              fragpipe = fragpipe_score,
              maxquant = maxquant_score,
              pd_tmt = pd_tmt_score)
  
  best_format <- names(which.max(scores))
  best_score <- max(scores)
  
  if (best_score < 0.3) {
    best_format <- "unknown"
  }
  
  return(list(
    format = best_format,
    confidence = best_score,
    all_scores = scores
  ))
}

#' Import DIA-NN Data
#' 
#' @param filepath Path to DIA-NN report file
#' @param use_precursor_norm Whether to use Precursor.Normalised values
#' @return List with data, data_type, and column_mapping
#' @export
importDIANNData <- function(filepath, use_precursor_norm = TRUE) {
  log_info(paste("Starting DIA-NN import from:", filepath))
  
  # Check file exists
  if (!file.exists(filepath)) {
    stop("File not found: ", filepath)
  }
  
  # Read data
  data <- tryCatch({
    vroom::vroom(filepath, show_col_types = FALSE)
  }, error = function(e) {
    stop("Failed to read file: ", e$message)
  })
  
  log_info(sprintf("Read %d rows and %d columns", nrow(data), ncol(data)))
  
  # Check for required columns
  required_cols <- c("Protein.Group", "Stripped.Sequence", "Run")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    log_error(paste("Missing required columns:", paste(missing_cols, collapse = ', ')))
    log_error(paste("Available columns:", paste(head(names(data), 20), collapse = ', ')))
    stop("Missing required DIA-NN columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Convert Precursor.Normalised if needed
  if (use_precursor_norm && "Precursor.Normalised" %in% names(data)) {
    log_info("Converting Precursor.Normalised to numeric")
    data$Precursor.Normalised <- as.numeric(data$Precursor.Normalised)
  }
  
  # Determine quantity column
  quantity_col <- if (use_precursor_norm && "Precursor.Normalised" %in% names(data)) {
    "Precursor.Normalised"
  } else if ("Precursor.Quantity" %in% names(data)) {
    "Precursor.Quantity"
  } else {
    stop("No suitable quantity column found (Precursor.Normalised or Precursor.Quantity)")
  }
  
  log_info(paste("Using quantity column:", quantity_col))
  
  return(list(
    data = data,
    data_type = "peptide",
    column_mapping = list(
      protein_col = "Protein.Group",
      peptide_col = "Stripped.Sequence",
      run_col = "Run",
      quantity_col = quantity_col,
      qvalue_col = "Q.Value",
      pg_qvalue_col = "PG.Q.Value"
    )
  ))
}

#' Import Spectronaut Data
#' 
#' @param filepath Path to Spectronaut report file
#' @param quantity_type Either "pg" for protein or "pep" for peptide quantities
#' @return List with data, data_type, and column_mapping
#' @export
importSpectronautData <- function(filepath, quantity_type = "pg") {
  data <- vroom::vroom(filepath, show_col_types = FALSE)
  
  data_type <- if (quantity_type == "pg") "protein" else "peptide"
  
  return(list(
    data = data,
    data_type = data_type,
    column_mapping = list(
      protein_col = "PG.ProteinGroups",
      peptide_col = if (quantity_type == "pep") "EG.ModifiedSequence" else NULL,
      run_col = "R.FileName",
      quantity_col = if (quantity_type == "pg") "PG.Quantity" else "PEP.Quantity",
      qvalue_col = if (quantity_type == "pg") "PG.Qvalue" else "EG.Qvalue"
    )
  ))
}

#' Import FragPipe Data
#' 
#' @param filepath Path to FragPipe output file
#' @param use_maxlfq Whether to use MaxLFQ intensities
#' @return List with data, data_type, and column_mapping
#' @export
importFragPipeData <- function(filepath, use_maxlfq = TRUE) {
  log_info(paste("Starting FragPipe LFQ import from:", filepath))
  
  # Check file exists
  if (!file.exists(filepath)) {
    stop("File not found: ", filepath)
  }
  
  # Read data
  data <- tryCatch({
    vroom::vroom(filepath, show_col_types = FALSE)
  }, error = function(e) {
    stop("Failed to read file: ", e$message)
  })
  
  log_info(sprintf("Read %d rows and %d columns", nrow(data), ncol(data)))
  
  # Find Protein ID column (case-insensitive, handle variations)
  protein_id_candidates <- c("Protein ID", "Protein.ID", "Protein_ID", "Protein", "protein id", "protein.id", "protein_id")
  protein_id_col <- NULL
  
  for (candidate in protein_id_candidates) {
    if (candidate %in% names(data)) {
      protein_id_col <- candidate
      break
    }
  }
  
  # Case-insensitive search if exact match not found
  if (is.null(protein_id_col)) {
    names_lower <- tolower(names(data))
    for (candidate in tolower(protein_id_candidates)) {
      idx <- which(names_lower == candidate)
      if (length(idx) > 0) {
        protein_id_col <- names(data)[idx[1]]
        break
      }
    }
  }
  
  if (is.null(protein_id_col)) {
    log_error(paste("Protein ID column not found. Available columns:", paste(head(names(data), 10), collapse = ", ")))
    stop("Required 'Protein ID' column not found in FragPipe file. Available columns: ", paste(head(names(data), 10), collapse = ", "))
  }
  
  log_info(sprintf("Found Protein ID column: %s", protein_id_col))
  
  # Find all columns ending with "Intensity" (case-insensitive)
  all_intensity_cols <- grep("Intensity$", names(data), value = TRUE, ignore.case = TRUE)
  
  if (length(all_intensity_cols) == 0) {
    log_error("No columns ending with 'Intensity' found in FragPipe file.")
    stop("No intensity columns found. Expected columns ending with 'Intensity' (e.g., 'WLP530_1 Intensity', 'WLP530_1 MaxLFQ Intensity')")
  }
  
  log_info(sprintf("Found %d columns ending with 'Intensity'", length(all_intensity_cols)))
  
  # Separate MaxLFQ and regular Intensity columns
  maxlfq_cols <- grep("MaxLFQ.*Intensity$", all_intensity_cols, value = TRUE, ignore.case = TRUE)
  regular_intensity_cols <- setdiff(all_intensity_cols, maxlfq_cols)
  
  # Select which intensity columns to use
  if (use_maxlfq && length(maxlfq_cols) > 0) {
    intensity_cols <- maxlfq_cols
    log_info(sprintf("Using MaxLFQ Intensity columns (%d columns)", length(intensity_cols)))
  } else {
    intensity_cols <- regular_intensity_cols
    log_info(sprintf("Using regular Intensity columns (%d columns)", length(intensity_cols)))
  }
  
  if (length(intensity_cols) == 0) {
    stop("No suitable intensity columns found. Check use_maxlfq parameter and file format.")
  }
  
  # Extract sample names by removing " Intensity" or " MaxLFQ Intensity" suffix
  sample_names <- gsub("\\s+(MaxLFQ\\s+)?Intensity$", "", intensity_cols, ignore.case = TRUE)
  
  # Select only the protein ID column and intensity columns for pivoting
  cols_to_keep <- c(protein_id_col, intensity_cols)
  data_subset <- data |> dplyr::select(dplyr::all_of(cols_to_keep))
  
  # Convert to long format
  log_info("Converting data from wide to long format...")
  long_data <- tryCatch({
    data_subset |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(intensity_cols),
        names_to = "Run",
        values_to = "Intensity"
      ) |>
      # Clean up Run column names (remove " Intensity" or " MaxLFQ Intensity" suffix)
      dplyr::mutate(
        Run = gsub("\\s+(MaxLFQ\\s+)?Intensity$", "", Run, ignore.case = TRUE)
      ) |>
      # Ensure Intensity is numeric
      dplyr::mutate(Intensity = as.numeric(Intensity)) |>
      # Rename protein ID column to standardized name
      dplyr::rename(Protein.Ids = !!rlang::sym(protein_id_col))
  }, error = function(e) {
    log_error(paste("Error converting to long format:", e$message))
    stop("Failed to convert data to long format: ", e$message)
  })
  
  log_info(sprintf("Converted to long format: %d rows", nrow(long_data)))
  log_info(sprintf("Unique proteins: %d, Unique samples: %d", 
                   length(unique(long_data$Protein.Ids)),
                   length(unique(long_data$Run))))
  
  return(list(
    data = long_data,
    data_type = "protein",
    column_mapping = list(
      protein_col = "Protein.Ids",
      run_col = "Run",
      quantity_col = "Intensity",
      qvalue_col = NULL  # FragPipe doesn't typically include q-values
    )
  ))
}

#' Import MaxQuant Data
#' 
#' @param filepath Path to MaxQuant proteinGroups.txt file
#' @param use_lfq Whether to use LFQ intensities
#' @param filter_contaminants Whether to filter contaminants
#' @return List with data, data_type, and column_mapping
#' @export
importMaxQuantData <- function(filepath, use_lfq = TRUE, filter_contaminants = TRUE) {
  data <- vroom::vroom(filepath, show_col_types = FALSE)
  
  # Filter contaminants and reverse hits if requested
  if (filter_contaminants) {
    if ("Potential.contaminant" %in% names(data)) {
      data <- data[data$Potential.contaminant != "+", ]
    }
    if ("Reverse" %in% names(data)) {
      data <- data[data$Reverse != "+", ]
    }
  }
  
  # Find intensity columns
  intensity_cols <- grep("^Intensity\\.", names(data), value = TRUE)
  lfq_cols <- grep("^LFQ\\.intensity\\.", names(data), value = TRUE)
  
  quantity_cols <- if (use_lfq && length(lfq_cols) > 0) {
    lfq_cols
  } else {
    intensity_cols
  }
  
  return(list(
    data = data,
    data_type = "protein",  # MaxQuant proteinGroups is protein-level
    column_mapping = list(
      protein_col = "Majority.protein.IDs",
      peptide_col = NULL,
      run_col = NULL,  # MaxQuant has wide format
      quantity_cols = quantity_cols,  # Multiple columns for wide format
      qvalue_col = "Q.value"  # If Perseus was used
    )
  ))
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

# Helper function for null-coalescing
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

