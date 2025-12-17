# ============================================================================
# mod_metab_design.R
# ============================================================================
# Purpose: Metabolomics design matrix builder Shiny module
#
# This module handles the experimental design setup for metabolomics data,
# including sample metadata assignment and contrast specification.
# ============================================================================

#' @title Metabolomics Design Matrix Module
#' @description A Shiny module for building and validating the experimental design
#'              matrix for metabolomics workflows.
#' @name mod_metab_design
NULL

#' @rdname mod_metab_design
#' @export
#' @importFrom shiny NS tagList fluidRow column wellPanel h3 h4 h5 p hr br radioButtons selectInput textInput actionButton uiOutput verbatimTextOutput icon tags fileInput numericInput checkboxInput
#' @importFrom DT DTOutput
mod_metab_design_ui <- function(id) {
    ns <- shiny::NS(id)
    
    shiny::tagList(
        shiny::fluidRow(
            shiny::column(12
                , shiny::wellPanel(
                    shiny::h3("Experimental Design Setup")
                    
                    , shiny::fluidRow(
                        # Left column: Design matrix input
                        shiny::column(6
                            , shiny::h4("Step 1: Import Design Matrix")
                            , shiny::radioButtons(
                                ns("design_source")
                                , "Design matrix source:"
                                , choices = c(
                                    "Upload file" = "upload"
                                    , "Generate from sample names" = "generate"
                                )
                                , selected = "upload"
                            )
                            
                            , shiny::conditionalPanel(
                                condition = paste0("input['", ns("design_source"), "'] == 'upload'")
                                , shiny::fileInput(
                                    ns("design_file")
                                    , "Upload Design Matrix"
                                    , accept = c(".csv", ".tsv", ".txt", ".xlsx")
                                )
                            )
                            
                            , shiny::conditionalPanel(
                                condition = paste0("input['", ns("design_source"), "'] == 'generate'")
                                , shiny::textInput(
                                    ns("group_pattern")
                                    , "Group extraction pattern (regex)"
                                    , value = ""
                                    , placeholder = "e.g., ^([A-Za-z]+)_ to extract prefix"
                                )
                                , shiny::actionButton(
                                    ns("generate_design")
                                    , "Generate Design"
                                    , class = "btn-info"
                                )
                            )
                            
                            , shiny::hr()
                            
                            , shiny::h4("Step 2: Column Mapping")
                            , shiny::selectInput(
                                ns("sample_id_col")
                                , "Sample ID Column"
                                , choices = NULL
                            )
                            , shiny::selectInput(
                                ns("group_col")
                                , "Group Column"
                                , choices = NULL
                            )
                            , shiny::selectInput(
                                ns("replicate_col")
                                , "Technical Replicate Column (optional)"
                                , choices = NULL
                            )
                            
                            , shiny::hr()
                            
                            , shiny::h4("Step 3: Contrast Specification")
                            , shiny::p("Define comparisons for differential analysis.")
                            , shiny::textInput(
                                ns("contrast_formula")
                                , "Contrast Formula"
                                , value = ""
                                , placeholder = "e.g., Treatment-Control"
                            )
                            , shiny::helpText("Use group names separated by minus sign. Multiple contrasts: comma-separated.")
                            , shiny::actionButton(
                                ns("add_contrast")
                                , "Add Contrast"
                                , class = "btn-secondary"
                            )
                            , shiny::br()
                            , shiny::br()
                            , shiny::uiOutput(ns("contrast_list"))
                        )
                        
                        # Right column: Preview and validation
                        , shiny::column(6
                            , shiny::h4("Design Matrix Preview")
                            , DT::DTOutput(ns("design_preview"))
                            
                            , shiny::hr()
                            
                            , shiny::h4("Validation Summary")
                            , shiny::uiOutput(ns("validation_summary"))
                            
                            , shiny::hr()
                            
                            , shiny::h4("Sample-Assay Matching")
                            , shiny::verbatimTextOutput(ns("sample_match_summary"))
                        )
                    )
                    
                    , shiny::hr()
                    
                    # Create S4 object button
                    , shiny::fluidRow(
                        shiny::column(12
                            , shiny::actionButton(
                                ns("create_s4")
                                , "Create MetaboliteAssayData Object"
                                , class = "btn-success"
                                , width = "100%"
                                , icon = shiny::icon("database")
                            )
                            , shiny::br()
                            , shiny::br()
                            , shiny::uiOutput(ns("s4_status"))
                        )
                    )
                )
            )
        )
    )
}

#' @rdname mod_metab_design
#' @export
#' @importFrom shiny moduleServer reactiveValues reactive observeEvent req renderUI renderText showNotification removeNotification updateSelectInput
#' @importFrom DT renderDT datatable
#' @importFrom logger log_info log_error log_warn
mod_metab_design_server <- function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns
        
        # Local reactive values
        local_data <- shiny::reactiveValues(
            design_matrix = NULL
            , contrasts = list()
            , sample_match_info = NULL
        )
        
        # Handle design file upload
        shiny::observeEvent(input$design_file, {
            shiny::req(input$design_file)
            
            tryCatch({
                file_path <- input$design_file$datapath
                file_name <- input$design_file$name
                
                # Determine format and read
                if (grepl("\\.xlsx$", file_name, ignore.case = TRUE)) {
                    design_df <- readxl::read_excel(file_path)
                } else if (grepl("\\.csv$", file_name, ignore.case = TRUE)) {
                    design_df <- vroom::vroom(file_path, delim = ",", show_col_types = FALSE)
                } else {
                    design_df <- vroom::vroom(file_path, delim = "\t", show_col_types = FALSE)
                }
                
                design_df <- as.data.frame(design_df)
                local_data$design_matrix <- design_df
                
                # Update column dropdowns
                col_names <- names(design_df)
                
                shiny::updateSelectInput(
                    session
                    , "sample_id_col"
                    , choices = col_names
                    , selected = findMatchingColumn(col_names, c("Sample", "SampleID", "Sample_ID", "sample_id", "sample"))
                )
                
                shiny::updateSelectInput(
                    session
                    , "group_col"
                    , choices = col_names
                    , selected = findMatchingColumn(col_names, c("Group", "group", "Condition", "condition", "Treatment", "treatment"))
                )
                
                shiny::updateSelectInput(
                    session
                    , "replicate_col"
                    , choices = c("(None)" = "", col_names)
                    , selected = ""
                )
                
                logger::log_info(sprintf("Design matrix loaded: %d samples, %d columns", nrow(design_df), ncol(design_df)))
                
            }, error = function(e) {
                logger::log_error(paste("Error reading design file:", e$message))
                shiny::showNotification(paste("Error reading design file:", e$message), type = "error")
            })
        })
        
        # Generate design from sample names
        shiny::observeEvent(input$generate_design, {
            shiny::req(workflow_data$column_mapping)
            
            tryCatch({
                sample_cols <- workflow_data$column_mapping$sample_columns
                
                if (length(sample_cols) == 0) {
                    stop("No sample columns available. Import data first.")
                }
                
                # Create basic design
                design_df <- data.frame(
                    Sample_ID = sample_cols
                    , stringsAsFactors = FALSE
                )
                
                # Try to extract group from pattern
                if (nzchar(input$group_pattern)) {
                    matches <- regmatches(sample_cols, regexpr(input$group_pattern, sample_cols, perl = TRUE))
                    design_df$Group <- ifelse(nchar(matches) > 0, matches, "Unknown")
                } else {
                    # Default: assign all to one group
                    design_df$Group <- "Sample"
                }
                
                local_data$design_matrix <- design_df
                
                # Update dropdowns
                shiny::updateSelectInput(session, "sample_id_col", choices = names(design_df), selected = "Sample_ID")
                shiny::updateSelectInput(session, "group_col", choices = names(design_df), selected = "Group")
                shiny::updateSelectInput(session, "replicate_col", choices = c("(None)" = "", names(design_df)), selected = "")
                
                logger::log_info(sprintf("Generated design matrix for %d samples", nrow(design_df)))
                shiny::showNotification("Design matrix generated", type = "message")
                
            }, error = function(e) {
                logger::log_error(paste("Error generating design:", e$message))
                shiny::showNotification(paste("Error:", e$message), type = "error")
            })
        })
        
        # Render design preview
        output$design_preview <- DT::renderDT({
            shiny::req(local_data$design_matrix)
            
            DT::datatable(
                local_data$design_matrix
                , options = list(
                    pageLength = 10
                    , scrollX = TRUE
                    , dom = 'frtip'
                )
                , rownames = FALSE
                , class = "compact stripe"
            )
        })
        
        # Add contrast
        shiny::observeEvent(input$add_contrast, {
            shiny::req(input$contrast_formula)
            
            # Parse contrast formula
            contrast_str <- trimws(input$contrast_formula)
            if (nzchar(contrast_str)) {
                # Split on comma for multiple contrasts
                new_contrasts <- strsplit(contrast_str, ",")[[1]]
                new_contrasts <- trimws(new_contrasts)
                
                # Add to list
                for (contrast in new_contrasts) {
                    if (nzchar(contrast) && !contrast %in% local_data$contrasts) {
                        local_data$contrasts <- c(local_data$contrasts, contrast)
                    }
                }
                
                # Clear input
                shiny::updateTextInput(session, "contrast_formula", value = "")
            }
        })
        
        # Render contrast list
        output$contrast_list <- shiny::renderUI({
            contrasts <- local_data$contrasts
            
            if (length(contrasts) == 0) {
                return(shiny::p(shiny::icon("info-circle"), " No contrasts defined yet.", style = "color: #666;"))
            }
            
            contrast_items <- lapply(seq_along(contrasts), function(i) {
                shiny::tags$div(
                    style = "margin-bottom: 5px;"
                    , shiny::tags$span(
                        class = "badge badge-primary"
                        , style = "margin-right: 5px;"
                        , contrasts[[i]]
                    )
                    , shiny::actionLink(
                        ns(paste0("remove_contrast_", i))
                        , shiny::icon("times")
                        , style = "color: red;"
                    )
                )
            })
            
            shiny::tagList(
                shiny::h5("Defined Contrasts:")
                , contrast_items
            )
        })
        
        # Validation summary
        output$validation_summary <- shiny::renderUI({
            shiny::req(local_data$design_matrix, input$sample_id_col, input$group_col)
            
            dm <- local_data$design_matrix
            errors <- character(0)
            warnings <- character(0)
            info <- character(0)
            
            # Check sample ID column
            if (!input$sample_id_col %in% names(dm)) {
                errors <- c(errors, "Sample ID column not found")
            } else {
                n_samples <- nrow(dm)
                n_unique <- length(unique(dm[[input$sample_id_col]]))
                info <- c(info, sprintf("Samples: %d", n_samples))
                
                if (n_unique != n_samples) {
                    warnings <- c(warnings, sprintf("Duplicate sample IDs: %d duplicates", n_samples - n_unique))
                }
            }
            
            # Check group column
            if (!input$group_col %in% names(dm)) {
                errors <- c(errors, "Group column not found")
            } else {
                groups <- unique(dm[[input$group_col]])
                info <- c(info, sprintf("Groups: %s", paste(groups, collapse = ", ")))
                
                # Check group sizes
                group_sizes <- table(dm[[input$group_col]])
                min_size <- min(group_sizes)
                if (min_size < 3) {
                    warnings <- c(warnings, sprintf("Small group(s) detected: min size = %d", min_size))
                }
            }
            
            # Check sample matching with data
            if (!is.null(workflow_data$column_mapping)) {
                data_samples <- workflow_data$column_mapping$sample_columns
                design_samples <- dm[[input$sample_id_col]]
                
                matched <- sum(design_samples %in% data_samples)
                unmatched_design <- sum(!design_samples %in% data_samples)
                unmatched_data <- sum(!data_samples %in% design_samples)
                
                info <- c(info, sprintf("Matched samples: %d/%d", matched, length(data_samples)))
                
                if (unmatched_design > 0) {
                    warnings <- c(warnings, sprintf("%d design samples not in data", unmatched_design))
                }
                if (unmatched_data > 0) {
                    warnings <- c(warnings, sprintf("%d data samples not in design", unmatched_data))
                }
                
                local_data$sample_match_info <- list(
                    matched = matched
                    , unmatched_design = unmatched_design
                    , unmatched_data = unmatched_data
                    , total_data = length(data_samples)
                )
            }
            
            # Build output
            if (length(errors) > 0) {
                return(shiny::tags$div(
                    class = "alert alert-danger"
                    , shiny::icon("times-circle")
                    , shiny::tags$strong(" Errors:")
                    , shiny::tags$ul(lapply(errors, shiny::tags$li))
                ))
            }
            
            shiny::tagList(
                if (length(warnings) > 0) {
                    shiny::tags$div(
                        class = "alert alert-warning"
                        , shiny::icon("exclamation-triangle")
                        , " Warnings:"
                        , shiny::tags$ul(lapply(warnings, shiny::tags$li))
                    )
                }
                , shiny::tags$div(
                    class = "alert alert-success"
                    , shiny::icon("check-circle")
                    , shiny::tags$strong(" Design Valid")
                    , shiny::tags$ul(lapply(info, shiny::tags$li))
                )
            )
        })
        
        # Sample match summary
        output$sample_match_summary <- shiny::renderText({
            match_info <- local_data$sample_match_info
            
            if (is.null(match_info)) {
                return("Import data first to check sample matching.")
            }
            
            paste(
                sprintf("Matched: %d / %d samples", match_info$matched, match_info$total_data)
                , sprintf("Design-only: %d", match_info$unmatched_design)
                , sprintf("Data-only: %d", match_info$unmatched_data)
                , sep = "\n"
            )
        })
        
        # Create S4 object
        shiny::observeEvent(input$create_s4, {
            shiny::req(
                local_data$design_matrix
                , input$sample_id_col
                , input$group_col
                , workflow_data$data_tbl
                , workflow_data$column_mapping
            )
            
            shiny::showNotification(
                "Creating MetaboliteAssayData object..."
                , id = "metab_s4_working"
                , duration = NULL
            )
            
            tryCatch({
                dm <- local_data$design_matrix
                col_map <- workflow_data$column_mapping
                assay_list <- workflow_data$data_tbl
                
                # Create S4 object
                s4_obj <- createMetaboliteAssayData(
                    metabolite_data = assay_list
                    , design_matrix = dm
                    , metabolite_id_column = col_map$metabolite_id_col
                    , annotation_id_column = if (!is.null(col_map$annotation_col) && nzchar(col_map$annotation_col)) {
                        col_map$annotation_col
                    } else {
                        NA_character_
                    }
                    , sample_id = input$sample_id_col
                    , group_id = input$group_col
                    , technical_replicate_id = if (nzchar(input$replicate_col)) input$replicate_col else NA_character_
                    , database_identifier_type = "Unknown"
                    , internal_standard_regex = col_map$is_pattern
                    , args = workflow_data$config_list
                )
                
                # Initialize state manager if needed
                if (is.null(workflow_data$state_manager)) {
                    workflow_data$state_manager <- WorkflowState$new("metabolomics")
                }
                
                # Save initial state
                workflow_data$state_manager$saveState(
                    state_name = "metab_raw_data_s4"
                    , s4_data_object = s4_obj
                    , config_object = workflow_data$config_list
                    , description = "Initial MetaboliteAssayData object created"
                )
                
                # Store contrasts
                workflow_data$contrasts <- local_data$contrasts
                
                # Update tab status
                workflow_data$tab_status$design_matrix <- "complete"
                
                # Update processing log
                workflow_data$processing_log$design_matrix <- list(
                    timestamp = Sys.time()
                    , n_samples = nrow(dm)
                    , n_groups = length(unique(dm[[input$group_col]]))
                    , groups = unique(dm[[input$group_col]])
                    , n_contrasts = length(local_data$contrasts)
                    , contrasts = local_data$contrasts
                )
                
                logger::log_info("MetaboliteAssayData S4 object created successfully")
                
                shiny::removeNotification("metab_s4_working")
                shiny::showNotification(
                    "MetaboliteAssayData object created! Proceed to QC."
                    , type = "message"
                )
                
            }, error = function(e) {
                logger::log_error(paste("Error creating S4 object:", e$message))
                shiny::removeNotification("metab_s4_working")
                shiny::showNotification(paste("Error:", e$message), type = "error", duration = 10)
            })
        })
        
        # S4 status display
        output$s4_status <- shiny::renderUI({
            if (!is.null(workflow_data$tab_status$design_matrix) && 
                workflow_data$tab_status$design_matrix == "complete") {
                
                shiny::tags$div(
                    class = "alert alert-success"
                    , shiny::icon("check-circle")
                    , shiny::tags$strong(" MetaboliteAssayData Object Created")
                    , shiny::tags$br()
                    , "You can now proceed to the Quality Control tab."
                )
            } else {
                NULL
            }
        })
    })
}

