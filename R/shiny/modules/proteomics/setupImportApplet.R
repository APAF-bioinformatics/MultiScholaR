#' Setup & Import Applet UI
#' 
#' UI for data import and initial setup with multi-format support
#' 
#' @param id Module ID
#' 
#' @importFrom shiny NS fileInput numericInput textInput actionButton uiOutput
#' @importFrom shiny fluidRow column wellPanel h3 h4 tags hr radioButtons
#' @export
setupImportUi <- function(id) {
  ns <- shiny::NS(id)
  
  # Check if shinyFiles is available
  use_shiny_files <- requireNamespace("shinyFiles", quietly = TRUE)
  
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
                       accept = c(".tsv", ".txt", ".tab", ".csv"))
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
                                  "MaxQuant LFQ" = "maxquant"),
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
#' @importFrom shiny renderUI showNotification removeNotification
#' @importFrom logger log_info log_error log_warn
#' @importFrom vroom vroom
#' @importFrom ini read.ini
#' @export
setupImportServer <- function(id, workflow_data, experiment_paths, volumes = NULL) {
  message(sprintf("--- Entering setupImportServer ---"))
  message(sprintf("   setupImportServer Arg: id = %s", id))
  message(sprintf("   setupImportServer Arg: volumes is NULL = %s", is.null(volumes)))
  if (!is.null(volumes)) {
    message(sprintf("   setupImportServer Arg: volumes type = %s, class = %s", 
                    typeof(volumes), paste(class(volumes), collapse = ", ")))
  }
  
  shiny::moduleServer(id, function(input, output, session) {
    message(sprintf("   setupImportServer Step: Inside moduleServer function"))
    
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
      config_file_path = NULL
    )
    
    # Check if shinyFiles is available
    message("   setupImportServer Step: Checking if shinyFiles is available...")
    use_shiny_files <- requireNamespace("shinyFiles", quietly = TRUE)
    message(sprintf("   setupImportServer Step: shinyFiles available = %s", use_shiny_files))
    
    # Set up shinyFiles if available
    if (use_shiny_files) {
      message("   setupImportServer Step: Setting up shinyFiles...")
      
      # If volumes not passed in, create them
      if (is.null(volumes)) {
        message("   setupImportServer Step: volumes is NULL, creating new volumes...")
        volumes <- shinyFiles::getVolumes()
        message(sprintf("   setupImportServer Step: Created volumes. Type = %s, class = %s", 
                        typeof(volumes), paste(class(volumes), collapse = ", ")))
        
        # Try to inspect volumes
        message("   setupImportServer Step: Attempting to inspect volumes...")
        tryCatch({
          if (is.function(volumes)) {
            message("   setupImportServer Step: volumes is a FUNCTION")
            # Try calling it
            vol_result <- volumes()
            message(sprintf("   setupImportServer Step: Called volumes(). Result type = %s, class = %s", 
                            typeof(vol_result), paste(class(vol_result), collapse = ", ")))
            message(sprintf("   setupImportServer Step: volumes() result length = %d", length(vol_result)))
            if (length(vol_result) > 0) {
              message(sprintf("   setupImportServer Step: First few volume names: %s", 
                              paste(head(names(vol_result), 3), collapse = ", ")))
            }
          } else {
            message("   setupImportServer Step: volumes is NOT a function")
            message(sprintf("   setupImportServer Step: volumes length = %d", length(volumes)))
          }
        }, error = function(e) {
          message(sprintf("   setupImportServer ERROR inspecting volumes: %s", e$message))
        })
      } else {
        message("   setupImportServer Step: volumes passed in from parent")
        message(sprintf("   setupImportServer Step: Passed volumes type = %s, class = %s", 
                        typeof(volumes), paste(class(volumes), collapse = ", ")))
      }
      
      # Search results file
      message("   setupImportServer Step: Setting up search_results file chooser...")
      message(sprintf("   setupImportServer Step: About to call shinyFileChoose with volumes type = %s", 
                      typeof(volumes)))
      
      shinyFiles::shinyFileChoose(
        input, 
        "search_results", 
        roots = volumes, 
        session = session,
        filetypes = c("tsv", "txt", "tab", "csv")
      )
      
      message("   setupImportServer Step: shinyFileChoose for search_results completed")
      
      observeEvent(input$search_results, {
        message("   setupImportServer Step: search_results button clicked")
        message(sprintf("   setupImportServer Step: input$search_results is NULL = %s", 
                        is.null(input$search_results)))
        
        if (!is.null(input$search_results)) {
          message(sprintf("   setupImportServer Step: input$search_results type = %s, class = %s", 
                          typeof(input$search_results), paste(class(input$search_results), collapse = ", ")))
          message(sprintf("   setupImportServer Step: is.integer(input$search_results) = %s", 
                          is.integer(input$search_results)))
          
          if (!is.integer(input$search_results)) {
            message("   setupImportServer Step: Attempting to parse file paths...")
            message(sprintf("   setupImportServer Step: volumes type for parsing = %s", typeof(volumes)))
            
            tryCatch({
              file_info <- shinyFiles::parseFilePaths(volumes, input$search_results)
              message(sprintf("   setupImportServer Step: parseFilePaths returned. Type = %s, class = %s", 
                              typeof(file_info), paste(class(file_info), collapse = ", ")))
              message(sprintf("   setupImportServer Step: file_info nrow = %d", nrow(file_info)))
              
              if (nrow(file_info) > 0) {
                message("   setupImportServer Step: file_info has rows, extracting datapath...")
                local_data$search_results_file <- as.character(file_info$datapath[1])
                message(sprintf("   setupImportServer Step: Set search_results_file = %s", 
                                local_data$search_results_file))
                output$search_results_path <- renderText(local_data$search_results_file)
              } else {
                message("   setupImportServer Step: file_info has 0 rows")
              }
            }, error = function(e) {
              message(sprintf("   setupImportServer ERROR in parseFilePaths: %s", e$message))
            })
          } else {
            message("   setupImportServer Step: input$search_results is integer (cancelled)")
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
        # Read first few lines to detect format
        preview_lines <- readLines(file_path, n = 10)
        headers <- strsplit(preview_lines[1], "\t|,")[[1]]
        
        # Get filename for detection
        filename <- if (use_shiny_files) {
          basename(file_path)
        } else {
          input$search_results_standard$name
        }
        
        # Detect format based on column headers and filename
        format_info <- detectProteomicsFormat(
          headers = headers,
          filename = filename,
          preview_lines = preview_lines
        )
        
        local_data$detected_format <- format_info$format
        local_data$format_confidence <- format_info$confidence
        
        log_info("Detected format: {format_info$format} (confidence: {format_info$confidence})")
        
      }, error = function(e) {
        log_error("Error detecting file format: {e$message}")
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
                       value = TRUE)
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
    
    # Process imported data
    shiny::observeEvent(input$process_data, {
      # Get file paths
      search_results_path <- get_search_results_path()
      fasta_path <- get_fasta_path()
      
      shiny::req(search_results_path)
      shiny::req(fasta_path)
      
      local_data$processing <- TRUE
      format <- active_format()
      
      shiny::showNotification("Processing imported data...", 
                      id = "processing_notification",
                      duration = NULL)
      
              tryCatch({
          # Read data based on format
          log_info("Reading {format} data from {search_results_path}")
        
                  data_import_result <- switch(format,
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
          stop("Unsupported format: ", format)
        )
        
        # Store in workflow data
        workflow_data$data_tbl <- data_import_result$data
        workflow_data$data_format <- format
        workflow_data$data_type <- data_import_result$data_type  # "peptide" or "protein"
        workflow_data$column_mapping <- data_import_result$column_mapping
        
        # Read FASTA file path (we don't parse it here, just store the path)
        workflow_data$fasta_file_path <- fasta_path
        
        # Set organism info
        workflow_data$taxon_id <- input$taxon_id
        workflow_data$organism_name <- input$organism_name
        
        # Load configuration
        config_path <- if (use_shiny_files) {
          local_data$config_file_path
        } else {
          if (!is.null(input$config_file_standard)) input$config_file_standard$datapath else NULL
        }
        
        if (!is.null(config_path)) {
          log_info("Reading configuration from {config_path}")
          config_list <- ini::read.ini(config_path)
        } else {
          log_info("Using default configuration")
          config_list <- getDefaultProteomicsConfig()
        }
        workflow_data$config_list <- config_list
        
        # Read optional mapping files
        uniprot_path <- if (use_shiny_files) {
          local_data$uniprot_mapping_file
        } else {
          if (!is.null(input$uniprot_mapping_standard)) input$uniprot_mapping_standard$datapath else NULL
        }
        
        if (!is.null(uniprot_path)) {
          log_info("Reading UniProt mapping file")
          workflow_data$uniprot_mapping <- vroom::vroom(
            uniprot_path,
            show_col_types = FALSE
          )
        }
        
        uniparc_path <- if (use_shiny_files) {
          local_data$uniparc_mapping_file
        } else {
          if (!is.null(input$uniparc_mapping_standard)) input$uniparc_mapping_standard$datapath else NULL
        }
        
        if (!is.null(uniparc_path)) {
          log_info("Reading UniParc mapping file")
          workflow_data$uniparc_mapping <- vroom::vroom(
            uniparc_path,
            show_col_types = FALSE
          )
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
        workflow_data$processing_log$setup_import <- list(
          timestamp = Sys.time(),
          search_file = search_filename,
          detected_format = format,
          data_type = data_import_result$data_type,
          fasta_file = fasta_filename,
          taxon_id = input$taxon_id,
          organism = input$organism_name,
          n_rows = nrow(data_import_result$data),
          n_runs = length(unique(data_import_result$data[[data_import_result$column_mapping$run_col]])),
          n_proteins = length(unique(data_import_result$data[[data_import_result$column_mapping$protein_col]])),
          n_peptides = if (!is.null(data_import_result$column_mapping$peptide_col)) {
            length(unique(data_import_result$data[[data_import_result$column_mapping$peptide_col]]))
          } else {
            NA
          }
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
      
      summary_info <- workflow_data$processing_log$setup_import
      
      shiny::tags$div(
        class = "well",
        shiny::h4("Data Summary"),
        shiny::tags$ul(
          shiny::tags$li(paste("Data format:", toupper(summary_info$detected_format))),
          shiny::tags$li(paste("Data type:", summary_info$data_type)),
          shiny::tags$li(paste("Total rows:", format(summary_info$n_rows, big.mark = ","))),
          shiny::tags$li(paste("Number of runs:", summary_info$n_runs)),
          shiny::tags$li(paste("Number of protein groups:", format(summary_info$n_proteins, big.mark = ","))),
          if (!is.na(summary_info$n_peptides)) {
            shiny::tags$li(paste("Number of peptides:", format(summary_info$n_peptides, big.mark = ",")))
          },
          shiny::tags$li(paste("Organism:", summary_info$organism, "(taxon:", summary_info$taxon_id, ")"))
        )
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
  fragpipe_markers <- c("protein", "peptide", "modified.peptide", "charge",
                       "gene", "protein.description", "intensity", "spectral.count")
  fragpipe_found <- sum(fragpipe_markers %in% headers_lower)
  fragpipe_score <- fragpipe_found / length(fragpipe_markers)
  if (grepl("fragpipe|msfragger", filename_lower)) fragpipe_score <- fragpipe_score + 0.3
  
  # MaxQuant detection
  maxquant_score <- 0
  maxquant_markers <- c("proteins", "majority.protein.ids", "protein.names",
                       "gene.names", "peptide.counts", "unique.peptides",
                       "intensity", "lfq.intensity", "ms.ms.count")
  maxquant_found <- sum(maxquant_markers %in% headers_lower)
  maxquant_score <- maxquant_found / length(maxquant_markers)
  if (grepl("proteingroups", filename_lower)) maxquant_score <- maxquant_score + 0.3
  
  # Determine best match
  scores <- c(diann = diann_score, 
              spectronaut = spectronaut_score,
              fragpipe = fragpipe_score,
              maxquant = maxquant_score)
  
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
  data <- vroom::vroom(filepath, show_col_types = FALSE)
  
  # Convert Precursor.Normalised if needed
  if (use_precursor_norm && "Precursor.Normalised" %in% names(data)) {
    data$Precursor.Normalised <- as.numeric(data$Precursor.Normalised)
  }
  
  return(list(
    data = data,
    data_type = "peptide",
    column_mapping = list(
      protein_col = "Protein.Group",
      peptide_col = "Stripped.Sequence",
      run_col = "Run",
      quantity_col = if (use_precursor_norm) "Precursor.Normalised" else "Precursor.Quantity",
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
  data <- vroom::vroom(filepath, show_col_types = FALSE)
  
  # Determine if this is protein or peptide level data
  data_type <- if ("Peptide" %in% names(data)) "peptide" else "protein"
  
  # Find intensity columns
  intensity_cols <- grep("Intensity", names(data), value = TRUE)
  maxlfq_cols <- grep("MaxLFQ", names(data), value = TRUE)
  
  quantity_cols <- if (use_maxlfq && length(maxlfq_cols) > 0) {
    maxlfq_cols
  } else {
    intensity_cols
  }
  
  return(list(
    data = data,
    data_type = data_type,
    column_mapping = list(
      protein_col = "Protein",
      peptide_col = if (data_type == "peptide") "Peptide" else NULL,
      run_col = NULL,  # FragPipe has wide format with sample columns
      quantity_cols = quantity_cols,  # Multiple columns for wide format
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