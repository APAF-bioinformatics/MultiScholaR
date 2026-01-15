# ============================================================================
# mod_metab_import.R
# ============================================================================
# Purpose: Metabolomics data import Shiny module
#
# This module provides multi-assay file import with vendor format detection
# and dynamic column mapping for metabolomics data.
# ============================================================================

#' @title Metabolomics Import Module
#' @description A Shiny module for importing metabolomics data with vendor format
#'              detection and dynamic column mapping.
#' @name mod_metab_import
NULL

#' @rdname mod_metab_import
#' @export
#' @importFrom shiny NS tagList fluidRow column wellPanel h3 h4 h5 p hr br radioButtons selectInput textInput actionButton uiOutput verbatimTextOutput icon tags conditionalPanel helpText div
#' @importFrom shinyjs useShinyjs disable enable
mod_metab_import_ui <- function(id) {
    ns <- shiny::NS(id)
    
    # Check if shinyFiles is available
    use_shiny_files <- requireNamespace("shinyFiles", quietly = TRUE)
    
    shiny::tagList(
        shinyjs::useShinyjs()
        , shiny::fluidRow(
            shiny::column(12
                , shiny::wellPanel(
                    shiny::h3("Metabolomics Data Import")
                    
                    , shiny::fluidRow(
                        # Left column: File import
                        shiny::column(6
                            , shiny::h4("Step 1: Select Vendor Format")
                            , shiny::radioButtons(
                                ns("vendor_format")
                                , NULL
                                , choices = c(
                                    "MS-DIAL" = "msdial"
                                    , "Progenesis QI" = "progenesis"
                                    , "XCMS" = "xcms"
                                    , "Compound Discoverer" = "compound_discoverer"
                                    , "Other/Custom" = "custom"
                                )
                                , selected = "msdial"
                                , inline = TRUE
                            )
                            
                            , shiny::hr()
                            , shiny::h4("Step 2: Import Data Files")
                            , shiny::p("Import one or more assay files (e.g., LCMS_Pos, LCMS_Neg, GCMS)")
                            
                            # Assay 1
                            , shiny::wellPanel(
                                style = "background-color: #f8f9fa;"
                                , shiny::fluidRow(
                                    shiny::column(4
                                        , shiny::textInput(
                                            ns("assay1_name")
                                            , "Assay Name"
                                            , value = "LCMS_Pos"
                                        )
                                    )
                                    , shiny::column(8
                                        , if (use_shiny_files) {
                                            shiny::tagList(
                                                shinyFiles::shinyFilesButton(
                                                    ns("assay1_file")
                                                    , "Select File"
                                                    , "Choose Data File"
                                                    , multiple = FALSE
                                                    , icon = shiny::icon("file")
                                                )
                                                , shiny::br()
                                                , shiny::verbatimTextOutput(ns("assay1_path"), placeholder = TRUE)
                                            )
                                        } else {
                                            shiny::fileInput(
                                                ns("assay1_file_std")
                                                , NULL
                                                , accept = c(".tsv", ".txt", ".csv", ".xlsx")
                                            )
                                        }
                                    )
                                )
                            )
                            
                            # Assay 2 (optional)
                            , shiny::wellPanel(
                                style = "background-color: #f8f9fa;"
                                , shiny::fluidRow(
                                    shiny::column(4
                                        , shiny::textInput(
                                            ns("assay2_name")
                                            , "Assay Name (Optional)"
                                            , value = ""
                                            , placeholder = "e.g., LCMS_Neg"
                                        )
                                    )
                                    , shiny::column(8
                                        , if (use_shiny_files) {
                                            shiny::tagList(
                                                shinyFiles::shinyFilesButton(
                                                    ns("assay2_file")
                                                    , "Select File"
                                                    , "Choose Data File"
                                                    , multiple = FALSE
                                                    , icon = shiny::icon("file")
                                                )
                                                , shiny::br()
                                                , shiny::verbatimTextOutput(ns("assay2_path"), placeholder = TRUE)
                                            )
                                        } else {
                                            shiny::fileInput(
                                                ns("assay2_file_std")
                                                , NULL
                                                , accept = c(".tsv", ".txt", ".csv", ".xlsx")
                                            )
                                        }
                                    )
                                )
                            )
                        )
                        
                        # Right column: Column mapping
                        , shiny::column(6
                            , shiny::h4("Step 3: Column Mapping")
                            , shiny::p("Configure column mappings for your data. Dropdowns are auto-populated based on detected format.")
                            
                            , shiny::conditionalPanel(
                                condition = paste0("output['", ns("file_loaded"), "']")
                                , shiny::wellPanel(
                                    # Format detection status (hidden for custom)
                                    shiny::conditionalPanel(
                                        condition = paste0("input['", ns("vendor_format"), "'] != 'custom'")
                                        , shiny::uiOutput(ns("format_detection_status"))
                                        , shiny::hr()
                                    )
                                    
                                    # Custom format instructions
                                    , shiny::conditionalPanel(
                                        condition = paste0("input['", ns("vendor_format"), "'] == 'custom'")
                                        , shiny::div(
                                            class = "alert alert-info"
                                            , shiny::icon("edit")
                                            , shiny::strong(" Custom Format")
                                            , shiny::br()
                                            , "Type the exact column names from your file below."
                                        )
                                        , shiny::hr()
                                    )
                                    
                                    # Metabolite ID column - Dropdown mode (non-custom)
                                    , shiny::conditionalPanel(
                                        condition = paste0("input['", ns("vendor_format"), "'] != 'custom'")
                                        , shiny::fluidRow(
                                            shiny::column(6
                                                , shiny::selectInput(
                                                    ns("metabolite_id_col")
                                                    , "Metabolite ID Column"
                                                    , choices = NULL
                                                )
                                            )
                                            , shiny::column(6
                                                , shiny::uiOutput(ns("metabolite_id_status"))
                                            )
                                        )
                                    )
                                    
                                    # Metabolite ID column - Text input mode (custom)
                                    , shiny::conditionalPanel(
                                        condition = paste0("input['", ns("vendor_format"), "'] == 'custom'")
                                        , shiny::fluidRow(
                                            shiny::column(6
                                                , shiny::textInput(
                                                    ns("metabolite_id_col_custom")
                                                    , "Metabolite ID Column"
                                                    , value = ""
                                                    , placeholder = "e.g., Compound_ID"
                                                )
                                            )
                                            , shiny::column(6
                                                , shiny::uiOutput(ns("metabolite_id_status_custom"))
                                            )
                                        )
                                    )
                                    
                                    # Annotation column - Dropdown mode (non-custom)
                                    , shiny::conditionalPanel(
                                        condition = paste0("input['", ns("vendor_format"), "'] != 'custom'")
                                        , shiny::fluidRow(
                                            shiny::column(6
                                                , shiny::selectInput(
                                                    ns("annotation_col")
                                                    , "Annotation Column"
                                                    , choices = NULL
                                                )
                                            )
                                            , shiny::column(6
                                                , shiny::uiOutput(ns("annotation_status"))
                                            )
                                        )
                                    )
                                    
                                    # Annotation column - Text input mode (custom)
                                    , shiny::conditionalPanel(
                                        condition = paste0("input['", ns("vendor_format"), "'] == 'custom'")
                                        , shiny::fluidRow(
                                            shiny::column(6
                                                , shiny::textInput(
                                                    ns("annotation_col_custom")
                                                    , "Annotation Column (optional)"
                                                    , value = ""
                                                    , placeholder = "e.g., Metabolite_Name"
                                                )
                                            )
                                            , shiny::column(6
                                                , shiny::uiOutput(ns("annotation_status_custom"))
                                            )
                                        )
                                    )
                                    
                                    # Sample columns pattern (custom mode only)
                                    , shiny::conditionalPanel(
                                        condition = paste0("input['", ns("vendor_format"), "'] == 'custom'")
                                        , shiny::textInput(
                                            ns("sample_cols_pattern")
                                            , "Sample Column Pattern (Regex)"
                                            , value = ""
                                            , placeholder = "e.g., ^Sample_|_Rep[0-9]+$"
                                        )
                                        , shiny::helpText("Leave blank to auto-detect numeric columns as samples")
                                    )
                                    
                                    # Internal standard pattern
                                    , shiny::textInput(
                                        ns("is_pattern")
                                        , "Internal Standard Pattern (Regex)"
                                        , value = ""
                                        , placeholder = "e.g., ^IS_|_d[0-9]+$|ISTD"
                                    )
                                    , shiny::helpText("Regular expression to identify internal standards")
                                    
                                    , shiny::hr()
                                    
                                    # Sample column detection
                                    , shiny::h5("Detected Sample Columns")
                                    , shiny::verbatimTextOutput(ns("sample_columns_display"))
                                    
                                    # Available columns helper (custom mode)
                                    , shiny::conditionalPanel(
                                        condition = paste0("input['", ns("vendor_format"), "'] == 'custom'")
                                        , shiny::hr()
                                        , shiny::h5("Available Columns in File")
                                        , shiny::verbatimTextOutput(ns("available_columns_display"))
                                    )
                                )
                                
                                , shiny::hr()
                                
                                # Validation summary
                                , shiny::h4("Step 4: Validation")
                                , shiny::uiOutput(ns("validation_summary"))
                            )
                            
                            , shiny::conditionalPanel(
                                condition = paste0("!output['", ns("file_loaded"), "']")
                                , shiny::div(
                                    class = "alert alert-info"
                                    , shiny::icon("info-circle")
                                    , " Import a data file to configure column mappings."
                                )
                            )
                        )
                    )
                    
                    , shiny::hr()
                    
                    # Process button
                    , shiny::fluidRow(
                        shiny::column(12
                            , shiny::actionButton(
                                ns("process_import")
                                , "Process Imported Data"
                                , class = "btn-success"
                                , width = "100%"
                                , icon = shiny::icon("check")
                            )
                            , shiny::br()
                            , shiny::br()
                            , shiny::uiOutput(ns("import_status"))
                        )
                    )
                )
            )
        )
    )
}

#' @rdname mod_metab_import
#' @export
#' @importFrom shiny moduleServer reactiveValues reactive observeEvent req renderUI renderText showNotification removeNotification updateSelectInput outputOptions
#' @importFrom shinyFiles shinyFileChoose parseFilePaths getVolumes
#' @importFrom logger log_info log_error log_warn
mod_metab_import_server <- function(id, workflow_data, experiment_paths, volumes = NULL) {
    message("--- Entering mod_metab_import_server ---")
    message(sprintf("   mod_metab_import_server: volumes param is NULL = %s", is.null(volumes)))
    
    shiny::moduleServer(id, function(input, output, session) {
        message("   mod_metab_import_server: Inside moduleServer function")
        
        # Check if shinyFiles is available
        use_shiny_files <- requireNamespace("shinyFiles", quietly = TRUE)
        message(sprintf("   mod_metab_import_server: shinyFiles available = %s", use_shiny_files))
        
        # Local reactive values
        local_data <- shiny::reactiveValues(
            assay1_file = NULL
            , assay1_data = NULL
            , assay1_import_result = NULL
            , assay2_file = NULL
            , assay2_data = NULL
            , assay2_import_result = NULL
            , detected_format = NULL
            , format_confidence = NULL
            , all_headers = NULL
        )
        
        # Set up shinyFiles if available
        if (use_shiny_files) {
            if (is.null(volumes)) {
                message("   mod_metab_import_server: volumes is NULL, creating from getVolumes()")
                volumes <- shinyFiles::getVolumes()()
            }
            message(sprintf("   mod_metab_import_server: volumes type = %s, length = %d", 
                          typeof(volumes), length(volumes)))
            if (length(volumes) > 0) {
                message(sprintf("   mod_metab_import_server: volumes names = %s", 
                              paste(names(volumes), collapse = ", ")))
            } else {
                message("   mod_metab_import_server: WARNING - volumes is empty!")
            }
            
            # Assay 1 file chooser
            shinyFiles::shinyFileChoose(
                input
                , "assay1_file"
                , roots = volumes
                , session = session
                , filetypes = c("tsv", "txt", "csv", "xlsx")
            )
            
            shiny::observeEvent(input$assay1_file, {
                if (!is.null(input$assay1_file) && !is.integer(input$assay1_file)) {
                    tryCatch({
                        file_info <- shinyFiles::parseFilePaths(volumes, input$assay1_file)
                        if (nrow(file_info) > 0) {
                            local_data$assay1_file <- as.character(file_info$datapath[1])
                            output$assay1_path <- shiny::renderText(local_data$assay1_file)
                            
                            # Import and detect format
                            import_data()
                        }
                    }, error = function(e) {
                        logger::log_error(paste("Error parsing file path:", e$message))
                    })
                }
            })
            
            # Assay 2 file chooser
            shinyFiles::shinyFileChoose(
                input
                , "assay2_file"
                , roots = volumes
                , session = session
                , filetypes = c("tsv", "txt", "csv", "xlsx")
            )
            
            shiny::observeEvent(input$assay2_file, {
                if (!is.null(input$assay2_file) && !is.integer(input$assay2_file)) {
                    tryCatch({
                        file_info <- shinyFiles::parseFilePaths(volumes, input$assay2_file)
                        if (nrow(file_info) > 0) {
                            local_data$assay2_file <- as.character(file_info$datapath[1])
                            output$assay2_path <- shiny::renderText(local_data$assay2_file)
                        }
                    }, error = function(e) {
                        logger::log_error(paste("Error parsing file path:", e$message))
                    })
                }
            })
        }
        
        # Import data and detect format
        import_data <- function() {
            shiny::req(local_data$assay1_file)
            
            tryCatch({
                # Read headers first for format detection
                headers <- tryCatch({
                    con <- file(local_data$assay1_file, "r")
                    first_line <- readLines(con, n = 1)
                    close(con)

                    # Detect delimiter
                    raw_headers <- if (grepl(",", first_line)) {
                        strsplit(first_line, ",")[[1]]
                    } else {
                        strsplit(first_line, "\t")[[1]]
                    }

                    # Strip surrounding quotes from headers (CSV quoting)
                    gsub('^"|"$', '', raw_headers)
                }, error = function(e) {
                    character(0)
                })
                
                if (length(headers) == 0) {
                    stop("Could not read headers from file")
                }
                
                local_data$all_headers <- headers
                
                # Detect format
                format_info <- detectMetabolomicsFormat(
                    headers = headers
                    , filename = basename(local_data$assay1_file)
                )
                
                local_data$detected_format <- format_info$format
                local_data$format_confidence <- format_info$confidence
                
                # Import data based on format
                import_result <- switch(format_info$format
                    , "msdial" = importMSDIALData(local_data$assay1_file)
                    , importMSDIALData(local_data$assay1_file)  # Default to MS-DIAL parser
                )
                
                local_data$assay1_import_result <- import_result
                local_data$assay1_data <- import_result$data
                
                # Update column mapping dropdowns
                shiny::updateSelectInput(
                    session
                    , "metabolite_id_col"
                    , choices = headers
                    , selected = import_result$detected_columns$metabolite_id
                )
                
                shiny::updateSelectInput(
                    session
                    , "annotation_col"
                    , choices = c("(None)" = "", headers)
                    , selected = import_result$detected_columns$annotation
                )
                
                # Set IS pattern from defaults
                if (!is.null(import_result$is_pattern) && !is.na(import_result$is_pattern)) {
                    shiny::updateTextInput(session, "is_pattern", value = import_result$is_pattern)
                }
                
                logger::log_info(sprintf(
                    "Imported assay 1: %d rows, %d columns, format: %s"
                    , nrow(import_result$data)
                    , ncol(import_result$data)
                    , format_info$format
                ))
                
            }, error = function(e) {
                logger::log_error(paste("Error importing data:", e$message))
                shiny::showNotification(paste("Error importing data:", e$message), type = "error")
            })
        }
        
        # Check if file is loaded
        output$file_loaded <- shiny::reactive({
            !is.null(local_data$assay1_data)
        })
        shiny::outputOptions(output, "file_loaded", suspendWhenHidden = FALSE)
        
        # Format detection status
        output$format_detection_status <- shiny::renderUI({
            shiny::req(local_data$detected_format)
            
            confidence_pct <- round(local_data$format_confidence * 100)
            color_class <- if (confidence_pct >= 70) "success" else if (confidence_pct >= 40) "warning" else "danger"
            
            format_display <- switch(local_data$detected_format
                , "msdial" = "MS-DIAL"
                , "progenesis" = "Progenesis QI"
                , "xcms" = "XCMS"
                , "compound_discoverer" = "Compound Discoverer"
                , "Unknown"
            )
            
            shiny::tags$div(
                class = paste("alert", paste0("alert-", color_class))
                , shiny::tags$strong("Detected format: ")
                , format_display
                , shiny::tags$br()
                , shiny::tags$small(sprintf("Confidence: %d%%", confidence_pct))
            )
        })
        
        # Metabolite ID validation status
        output$metabolite_id_status <- shiny::renderUI({
            shiny::req(local_data$assay1_data, input$metabolite_id_col)
            
            if (input$metabolite_id_col %in% names(local_data$assay1_data)) {
                n_unique <- length(unique(local_data$assay1_data[[input$metabolite_id_col]]))
                shiny::tags$span(
                    shiny::icon("check-circle", style = "color: green;")
                    , sprintf(" %d unique IDs", n_unique)
                )
            } else {
                shiny::tags$span(
                    shiny::icon("times-circle", style = "color: red;")
                    , " Column not found"
                )
            }
        })
        
        # Annotation validation status
        output$annotation_status <- shiny::renderUI({
            shiny::req(local_data$assay1_data)
            
            if (!is.null(input$annotation_col) && nzchar(input$annotation_col)) {
                if (input$annotation_col %in% names(local_data$assay1_data)) {
                    shiny::tags$span(
                        shiny::icon("check-circle", style = "color: green;")
                        , " OK"
                    )
                } else {
                    shiny::tags$span(
                        shiny::icon("times-circle", style = "color: red;")
                        , " Column not found"
                    )
                }
            } else {
                shiny::tags$span(
                    shiny::icon("minus-circle", style = "color: gray;")
                    , " Optional"
                )
            }
        })
        
        # Display detected sample columns
        output$sample_columns_display <- shiny::renderText({
            shiny::req(local_data$assay1_import_result)
            
            sample_cols <- local_data$assay1_import_result$sample_columns
            if (length(sample_cols) > 10) {
                paste(
                    paste(head(sample_cols, 10), collapse = ", ")
                    , sprintf("... and %d more", length(sample_cols) - 10)
                )
            } else {
                paste(sample_cols, collapse = ", ")
            }
        })
        
        # Display available columns (for custom mode)
        output$available_columns_display <- shiny::renderText({
            shiny::req(local_data$all_headers)
            paste(local_data$all_headers, collapse = ", ")
        })
        
        # Custom mode: Metabolite ID validation status
        output$metabolite_id_status_custom <- shiny::renderUI({
            shiny::req(local_data$assay1_data)
            
            col_name <- input$metabolite_id_col_custom
            
            if (is.null(col_name) || !nzchar(col_name)) {
                return(shiny::tags$span(
                    shiny::icon("question-circle", style = "color: gray;")
                    , " Enter column name"
                ))
            }
            
            # Case-insensitive matching
            headers_lower <- tolower(names(local_data$assay1_data))
            col_lower <- tolower(col_name)
            
            if (col_lower %in% headers_lower) {
                actual_col <- names(local_data$assay1_data)[which(headers_lower == col_lower)[1]]
                n_unique <- length(unique(local_data$assay1_data[[actual_col]]))
                shiny::tags$span(
                    shiny::icon("check-circle", style = "color: green;")
                    , sprintf(" Found: %d unique IDs", n_unique)
                )
            } else {
                shiny::tags$span(
                    shiny::icon("times-circle", style = "color: red;")
                    , " Column not found"
                )
            }
        })
        
        # Custom mode: Annotation validation status
        output$annotation_status_custom <- shiny::renderUI({
            shiny::req(local_data$assay1_data)
            
            col_name <- input$annotation_col_custom
            
            if (is.null(col_name) || !nzchar(col_name)) {
                return(shiny::tags$span(
                    shiny::icon("minus-circle", style = "color: gray;")
                    , " Optional"
                ))
            }
            
            # Case-insensitive matching
            headers_lower <- tolower(names(local_data$assay1_data))
            col_lower <- tolower(col_name)
            
            if (col_lower %in% headers_lower) {
                shiny::tags$span(
                    shiny::icon("check-circle", style = "color: green;")
                    , " Found"
                )
            } else {
                shiny::tags$span(
                    shiny::icon("times-circle", style = "color: red;")
                    , " Column not found"
                )
            }
        })
        
        # Helper to get the effective metabolite ID column (dropdown or custom)
        get_metabolite_id_col <- shiny::reactive({
            if (input$vendor_format == "custom") {
                col_name <- input$metabolite_id_col_custom
                if (!is.null(col_name) && nzchar(col_name)) {
                    # Find actual column name (case-insensitive)
                    headers <- names(local_data$assay1_data)
                    headers_lower <- tolower(headers)
                    col_lower <- tolower(col_name)
                    if (col_lower %in% headers_lower) {
                        return(headers[which(headers_lower == col_lower)[1]])
                    }
                }
                return(col_name)
            } else {
                return(input$metabolite_id_col)
            }
        })
        
        # Helper to get the effective annotation column
        get_annotation_col <- shiny::reactive({
            if (input$vendor_format == "custom") {
                col_name <- input$annotation_col_custom
                if (!is.null(col_name) && nzchar(col_name)) {
                    headers <- names(local_data$assay1_data)
                    headers_lower <- tolower(headers)
                    col_lower <- tolower(col_name)
                    if (col_lower %in% headers_lower) {
                        return(headers[which(headers_lower == col_lower)[1]])
                    }
                }
                return(col_name)
            } else {
                return(input$annotation_col)
            }
        })
        
        # Helper to get sample columns (with custom pattern support)
        get_sample_columns <- shiny::reactive({
            shiny::req(local_data$assay1_data)
            
            if (input$vendor_format == "custom" && nzchar(input$sample_cols_pattern)) {
                # Use regex pattern to match columns
                pattern <- input$sample_cols_pattern
                all_cols <- names(local_data$assay1_data)
                matched <- all_cols[grepl(pattern, all_cols, ignore.case = TRUE)]
                
                if (length(matched) > 0) {
                    return(matched)
                }
            }
            
            # Default: use detected sample columns from import
            if (!is.null(local_data$assay1_import_result)) {
                return(local_data$assay1_import_result$sample_columns)
            }
            
            # Fallback: detect numeric columns
            return(names(local_data$assay1_data)[sapply(local_data$assay1_data, is.numeric)])
        })
        
        # Validation summary
        output$validation_summary <- shiny::renderUI({
            shiny::req(local_data$assay1_data)
            
            metabolite_col <- get_metabolite_id_col()
            sample_cols <- get_sample_columns()
            
            shiny::req(metabolite_col)
            
            validation <- validateColumnMapping(
                data = local_data$assay1_data
                , metabolite_id_column = metabolite_col
                , sample_columns = sample_cols
            )
            
            if (validation$valid) {
                shiny::tagList(
                    shiny::tags$div(
                        class = "alert alert-success"
                        , shiny::icon("check-circle")
                        , shiny::tags$strong(" Validation Passed")
                    )
                    , shiny::tags$ul(
                        shiny::tags$li(sprintf("Metabolites: %d", validation$summary$n_metabolites))
                        , shiny::tags$li(sprintf("Samples: %d", validation$summary$n_samples))
                        , shiny::tags$li(sprintf("Missing values: %.1f%%", validation$summary$pct_missing))
                    )
                    , if (length(validation$warnings) > 0) {
                        shiny::tags$div(
                            class = "alert alert-warning"
                            , shiny::icon("exclamation-triangle")
                            , " Warnings:"
                            , shiny::tags$ul(
                                lapply(validation$warnings, shiny::tags$li)
                            )
                        )
                    }
                )
            } else {
                shiny::tags$div(
                    class = "alert alert-danger"
                    , shiny::icon("times-circle")
                    , shiny::tags$strong(" Validation Failed")
                    , shiny::tags$ul(
                        lapply(validation$errors, shiny::tags$li)
                    )
                )
            }
        })
        
        # Process import
        shiny::observeEvent(input$process_import, {
            shiny::req(local_data$assay1_data)
            
            metabolite_col <- get_metabolite_id_col()
            annotation_col <- get_annotation_col()
            sample_cols <- get_sample_columns()
            
            shiny::req(metabolite_col)
            
            shiny::showNotification(
                "Processing imported data..."
                , id = "metab_import_working"
                , duration = NULL
            )
            
            tryCatch({
                # Build assay list
                assay_list <- list()
                assay_list[[input$assay1_name]] <- local_data$assay1_data
                
                # Add second assay if provided
                if (!is.null(local_data$assay2_file) && nzchar(input$assay2_name)) {
                    assay2_import <- importMSDIALData(local_data$assay2_file)
                    assay_list[[input$assay2_name]] <- assay2_import$data
                }
                
                # Store in workflow_data
                workflow_data$data_tbl <- assay_list
                workflow_data$data_format <- if (input$vendor_format == "custom") "custom" else local_data$detected_format
                workflow_data$data_type <- "metabolite"
                
                workflow_data$column_mapping <- list(
                    metabolite_id_col = metabolite_col
                    , annotation_col = if (!is.null(annotation_col) && nzchar(annotation_col)) annotation_col else NULL
                    , sample_columns = sample_cols
                    , is_pattern = if (nzchar(input$is_pattern)) input$is_pattern else NA_character_
                )
                
                # Initialize state manager workflow type
                # Workflow type determines module set, not vendor format
                if (!is.null(workflow_data$state_manager)) {
                    workflow_data$state_manager$setWorkflowType("metabolomics_standard")
                    logger::log_info("Workflow type set to: metabolomics_standard")
                }
                
                # Update processing log
                workflow_data$processing_log$setup_import <- list(
                    timestamp = Sys.time()
                    , n_assays = length(assay_list)
                    , assay_names = names(assay_list)
                    , detected_format = workflow_data$data_format
                    , n_metabolites = sapply(assay_list, function(a) {
                        if (metabolite_col %in% names(a)) {
                            length(unique(a[[metabolite_col]]))
                        } else {
                            nrow(a)
                        }
                    })
                    , n_samples = length(sample_cols)
                )
                
                # Mark tab as complete - must replace entire list to trigger reactivity
                updated_status <- workflow_data$tab_status
                updated_status$setup_import <- "complete"
                workflow_data$tab_status <- updated_status
                
                logger::log_info(sprintf(
                    "Metabolomics import complete: %d assays, %d total metabolites"
                    , length(assay_list)
                    , sum(sapply(assay_list, nrow))
                ))
                
                shiny::removeNotification("metab_import_working")
                shiny::showNotification("Data imported successfully!", type = "message")
                
            }, error = function(e) {
                logger::log_error(paste("Error processing import:", e$message))
                shiny::removeNotification("metab_import_working")
                shiny::showNotification(paste("Error:", e$message), type = "error", duration = 10)
            })
        })
        
        # Import status display
        output$import_status <- shiny::renderUI({
            if (!is.null(workflow_data$tab_status$setup_import) && 
                workflow_data$tab_status$setup_import == "complete") {
                
                log_info <- workflow_data$processing_log$setup_import
                
                shiny::tags$div(
                    class = "alert alert-success"
                    , shiny::icon("check-circle")
                    , shiny::tags$strong(" Import Complete")
                    , shiny::tags$br()
                    , sprintf(
                        "Format: %s | Assays: %d | Samples: %d"
                        , toupper(log_info$detected_format)
                        , log_info$n_assays
                        , log_info$n_samples
                    )
                )
            } else {
                NULL
            }
        })
    })
}

