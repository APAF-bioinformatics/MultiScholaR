# ============================================================================
# mod_metab_design.R
# ============================================================================
# Purpose: Metabolomics design matrix host module - 1:1 with proteomics
#
# This module embeds mod_metab_design_builder.R and handles orchestration,
# including Import Existing Design, S4 object creation, and state management.
# ============================================================================

#' @title Metabolomics Design Matrix Applet Module
#'
#' @description A Shiny module that serves as the main host for the design
#' matrix creation workflow step. Embeds mod_metab_design_builder_ui/server
#' and handles Import Existing Design, S4 object creation, and state saving.
#' Architecture is 1:1 with mod_prot_design.R.
#'
#' @param id Module ID
#' @param workflow_data A reactive values object to store workflow data.
#' @param experiment_paths A list of paths for the current experiment.
#' @param volumes A list of volumes for shinyFiles (optional).
#' @param qc_trigger A reactive trigger for QC execution (optional).
#'
#' @name metabolomicsDesignMatrixAppletModule
NULL

#' @rdname metabolomicsDesignMatrixAppletModule
#' @export
#' @importFrom shiny NS tagList wellPanel h3 h4 p conditionalPanel div icon tags HTML fluidRow column actionButton
#' @importFrom DT DTOutput renderDT
mod_metab_design_ui <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::fluidRow(
            shiny::column(12
                , shiny::wellPanel(
                    shiny::fluidRow(
                        shiny::column(8
                            , shiny::h3("Design Matrix Builder")
                        )
                        , shiny::column(4
                            , shiny::actionButton(ns("show_import_modal"), "Import Existing Design"
                                , icon = shiny::icon("folder-open"), class = "btn-info pull-right")
                        )
                    )
                    , shiny::p("Use the tools below to define your experimental groups and contrasts for differential analysis. All assays (Pos/Neg modes) share the same design matrix.")

                    # Multi-assay info banner
                    , shiny::div(
                        class = "alert alert-info"
                        , style = "margin-top: 10px;"
                        , shiny::icon("info-circle")
                        , shiny::tags$strong(" Multi-Assay Workflow: ")
                        , "Design applies uniformly across all assay modes (LCMS_Pos, LCMS_Neg, GCMS, etc.). "
                        , "Sample assignments and contrasts are shared."
                    )

                    # Conditional panel: show builder when data is available
                    , shiny::conditionalPanel(
                        condition = paste0("output['", ns("data_available"), "']")
                        , if (exists("mod_metab_design_builder_ui")) {
                            mod_metab_design_builder_ui(ns("builder"))
                        } else {
                            shiny::div("Design builder module not loaded")
                        }

                        # Saved results preview
                        , shiny::hr()
                        , shiny::h3("Saved Results Preview")
                        , shiny::p("This section shows the design matrix and contrasts saved to the workflow.")
                        , shiny::conditionalPanel(
                            condition = paste0("output['", ns("design_matrix_exists"), "']")
                            , shiny::wellPanel(
                                shiny::h4("Current Design Matrix")
                                , DT::DTOutput(ns("design_matrix_preview"))
                                , shiny::br()
                                , shiny::h4("Defined Contrasts")
                                , DT::DTOutput(ns("contrasts_preview"))
                                , shiny::br()
                                , shiny::h4("Assays Included")
                                , shiny::verbatimTextOutput(ns("assays_preview"))
                            )
                        )
                    )

                    # Conditional panel: message if no data
                    , shiny::conditionalPanel(
                        condition = paste0("!output['", ns("data_available"), "']")
                        , shiny::div(
                            class = "alert alert-info"
                            , shiny::icon("info-circle")
                            , " Please complete the 'Import' step first. The builder will appear once data is available."
                        )
                    )
                )
            )
        )
    )
}

#' @rdname metabolomicsDesignMatrixAppletModule
#' @export
#' @importFrom shiny moduleServer reactive observeEvent req renderUI showNotification removeNotification outputOptions renderText
#' @importFrom logger log_info log_error log_warn
#' @importFrom utils write.table
#' @importFrom vroom vroom
#' @importFrom shinyFiles shinyDirButton shinyDirChoose parseDirPath getVolumes
#' @importFrom jsonlite write_json read_json
mod_metab_design_server <- function(id, workflow_data, experiment_paths, volumes = NULL, qc_trigger = NULL) {
    message(sprintf("--- Entering mod_metab_design_server ---"))
    message(sprintf("   mod_metab_design_server Arg: id = %s", id))

    shiny::moduleServer(id, function(input, output, session) {
        message(sprintf("   mod_metab_design_server Step: Inside moduleServer function"))

        # == Setup shinyFiles =======================================================
        resolved_volumes <- shiny::isolate({
            base_volumes <- if (is.function(volumes)) {
                volumes()
            } else if (is.null(volumes)) {
                shinyFiles::getVolumes()()
            } else {
                volumes
            }

            if (!is.null(experiment_paths) && !is.null(experiment_paths$base_dir) &&
                dir.exists(experiment_paths$base_dir)) {
                enhanced_volumes <- c("Project Base Dir" = experiment_paths$base_dir, base_volumes)
                logger::log_info(paste("Added base_dir to volumes:", experiment_paths$base_dir))
                enhanced_volumes
            } else {
                base_volumes
            }
        })

        shinyFiles::shinyDirChoose(input, "import_dir", roots = resolved_volumes, session = session)

        # == Modal Logic for Import =================================================

        shiny::observeEvent(input$show_import_modal, {
            ns <- session$ns
            shiny::showModal(shiny::modalDialog(
                title = "Import Existing Design Matrix"
                , shiny::p("Select the folder containing 'design_matrix.tab', assay data files, and optionally 'contrast_strings.tab'.")
                , shiny::helpText("Required files: design_matrix.tab, assay_manifest.txt, data_cln_*.tab files")
                , shinyFiles::shinyDirButton(ns("import_dir"), "Select Folder", "Choose a directory")
                , shiny::verbatimTextOutput(ns("import_dir_path"), placeholder = TRUE)
                , footer = shiny::tagList(
                    shiny::modalButton("Cancel")
                    , shiny::actionButton(ns("confirm_import"), "Import", class = "btn-primary")
                )
            ))
        })

        # Display selected directory path in modal
        output$import_dir_path <- shiny::renderText({
            shiny::req(input$import_dir)
            shinyFiles::parseDirPath(resolved_volumes, input$import_dir)
        })

        # == Handle Import Confirmation =============================================

        shiny::observeEvent(input$confirm_import, {
            shiny::req(input$import_dir)

            import_path <- shinyFiles::parseDirPath(resolved_volumes, input$import_dir)
            shiny::req(import_path)

            shiny::removeModal()

            # Define file paths
            design_file <- file.path(import_path, "design_matrix.tab")
            contrast_file <- file.path(import_path, "contrast_strings.tab")
            manifest_file <- file.path(import_path, "assay_manifest.txt")
            col_map_file <- file.path(import_path, "column_mapping.json")
            config_file <- file.path(import_path, "config.ini")

            # Validate required files exist
            if (!file.exists(design_file)) {
                msg <- "Import failed: 'design_matrix.tab' not found in the selected directory."
                logger::log_error(msg)
                shiny::showNotification(msg, type = "error", duration = 10)
                return()
            }

            if (!file.exists(manifest_file)) {
                msg <- "Import failed: 'assay_manifest.txt' not found. Cannot determine which assay files to load."
                logger::log_error(msg)
                shiny::showNotification(msg, type = "error", duration = 10)
                return()
            }

            shiny::showNotification("Importing design files...", id = "importing_design", duration = NULL)

            tryCatch({
                # DEBUG66: Entry point
                message("DEBUG66: ========================================")
                message("DEBUG66: STARTING IMPORT EXISTING DESIGN")
                message(sprintf("DEBUG66: import_path = %s", import_path))
                message(sprintf("DEBUG66: config_file exists = %s", file.exists(config_file)))
                message(sprintf("DEBUG66: design_file exists = %s", file.exists(design_file)))
                message(sprintf("DEBUG66: manifest_file exists = %s", file.exists(manifest_file)))
                message("DEBUG66: ========================================")
                
                # --- 1. Load config.ini ---
                if (file.exists(config_file)) {
                    logger::log_info("Loading config.ini from import directory.")
                    message("DEBUG66: About to call readConfigFile()...")
                    workflow_data$config_list <- readConfigFile(file = config_file)
                    message("DEBUG66: readConfigFile() completed successfully")
                    message(sprintf("DEBUG66: config_list names: %s", paste(names(workflow_data$config_list), collapse = ", ")))
                    assign("config_list", workflow_data$config_list, envir = .GlobalEnv)
                    message("DEBUG66: config assigned to global env")
                    logger::log_info("Loaded config.ini and assigned to global environment.")
                } else {
                    # Try to load from experiment_paths or use default
                    default_config <- file.path(experiment_paths$source_dir, "config.ini")
                    if (file.exists(default_config)) {
                        workflow_data$config_list <- readConfigFile(file = default_config)
                        assign("config_list", workflow_data$config_list, envir = .GlobalEnv)
                        logger::log_info("Loaded config.ini from source_dir.")
                    } else {
                        logger::log_warn("No config.ini found. Using empty config.")
                        workflow_data$config_list <- list()
                    }
                }

                # --- 2. Load design matrix ---
                logger::log_info("Loading design_matrix.tab")
                imported_design <- vroom::vroom(design_file, show_col_types = FALSE)
                
                # DEBUG66: Log design matrix details
                message(sprintf("DEBUG66: Loaded design_matrix with %d rows, %d cols", 
                    nrow(imported_design), ncol(imported_design)))
                message(sprintf("DEBUG66: design_matrix columns: %s", 
                    paste(names(imported_design), collapse = ", ")))
                message(sprintf("DEBUG66: design_matrix$Run values: %s%s", 
                    paste(head(imported_design$Run, 5), collapse = ", "), 
                    if(nrow(imported_design) > 5) "..." else ""))

                # --- 3. Load assay manifest and assay data files ---
                logger::log_info("Loading assay_manifest.txt")
                assay_names <- readLines(manifest_file)
                assay_names <- assay_names[nzchar(assay_names)]  # Remove empty lines

                if (length(assay_names) == 0) {
                    stop("assay_manifest.txt is empty. No assays to load.")
                }

                logger::log_info(sprintf("Found %d assays in manifest: %s", length(assay_names), paste(assay_names, collapse = ", ")))

                # Load each assay data file
                assay_list <- list()
                for (assay_name in assay_names) {
                    assay_file <- file.path(import_path, paste0("data_cln_", assay_name, ".tab"))
                    if (!file.exists(assay_file)) {
                        stop(sprintf("Assay data file not found: %s", assay_file))
                    }
                    logger::log_info(sprintf("Loading assay: %s", assay_name))
                    assay_list[[assay_name]] <- vroom::vroom(assay_file, show_col_types = FALSE)
                }
                
                # DEBUG66: Log assay list details
                message(sprintf("DEBUG66: Loaded %d assays: %s", 
                    length(assay_list), paste(names(assay_list), collapse = ", ")))
                for (a_name in names(assay_list)) {
                    message(sprintf("DEBUG66: Assay '%s': %d rows, %d cols, cols: %s", 
                        a_name, nrow(assay_list[[a_name]]), ncol(assay_list[[a_name]]),
                        paste(head(names(assay_list[[a_name]]), 5), collapse = ", ")))
                }

                # --- 4. Load column mapping ---
                if (file.exists(col_map_file)) {
                    logger::log_info("Loading column_mapping.json")
                    col_map <- jsonlite::read_json(col_map_file, simplifyVector = TRUE)
                    workflow_data$column_mapping <- col_map
                    
                    # DEBUG66: Log column mapping from file
                    message(sprintf("DEBUG66: Loaded column_mapping from JSON"))
                    message(sprintf("DEBUG66: metabolite_id_col = '%s'", col_map$metabolite_id_col))
                    message(sprintf("DEBUG66: annotation_col = '%s'", col_map$annotation_col))
                    message(sprintf("DEBUG66: sample_columns count = %d", length(col_map$sample_columns)))
                    message(sprintf("DEBUG66: is_pattern = '%s'", col_map$is_pattern))
                } else {
                    # Infer column mapping from data structure
                    logger::log_warn("column_mapping.json not found. Inferring from data structure.")
                    first_assay <- assay_list[[1]]
                    all_cols <- names(first_assay)

                    # Try to detect metabolite ID column
                    id_candidates <- c("Peak ID", "Alignment ID", "Alignment.ID", "peak id", "alignment id")
                    metabolite_id_col <- NULL
                    for (candidate in id_candidates) {
                        if (candidate %in% all_cols) {
                            metabolite_id_col <- candidate
                            break
                        }
                    }
                    if (is.null(metabolite_id_col)) {
                        metabolite_id_col <- all_cols[1]  # Fallback to first column
                    }

                    # Try to detect annotation column
                    annot_candidates <- c("Metabolite name", "Metabolite.name", "Name", "name")
                    annotation_col <- NULL
                    for (candidate in annot_candidates) {
                        if (candidate %in% all_cols) {
                            annotation_col <- candidate
                            break
                        }
                    }

                    # Sample columns are those matching design matrix Run values
                    sample_cols <- intersect(imported_design$Run, all_cols)

                    col_map <- list(
                        metabolite_id_col = metabolite_id_col
                        , annotation_col = annotation_col
                        , sample_columns = sample_cols
                        , is_pattern = NA_character_
                    )
                    workflow_data$column_mapping <- col_map
                    logger::log_info(sprintf("Inferred column mapping: metabolite_id=%s, annotation=%s, %d samples",
                        metabolite_id_col, annotation_col, length(sample_cols)))
                }

                # --- 5. Load contrasts ---
                imported_contrasts <- if (file.exists(contrast_file)) {
                    contrast_strings <- readLines(contrast_file)

                    if (length(contrast_strings) > 0) {
                        contrast_info <- lapply(contrast_strings, function(contrast_string) {
                            clean_string <- gsub("^group", "", contrast_string)
                            clean_string <- gsub("-group", "-", clean_string)
                            friendly_name <- gsub("-", "_vs_", clean_string)
                            full_format <- paste0(friendly_name, "=", contrast_string)

                            list(
                                contrast_string = contrast_string
                                , friendly_name = friendly_name
                                , full_format = full_format
                            )
                        })

                        data.frame(
                            contrasts = sapply(contrast_info, function(x) x$contrast_string)
                            , friendly_names = sapply(contrast_info, function(x) x$friendly_name)
                            , full_format = sapply(contrast_info, function(x) x$full_format)
                            , stringsAsFactors = FALSE
                        )
                    } else {
                        NULL
                    }
                } else {
                    logger::log_info("No contrast_strings.tab found.")
                    NULL
                }
                
                # DEBUG66: Log contrasts
                if (!is.null(imported_contrasts)) {
                    message(sprintf("DEBUG66: Loaded %d contrasts", nrow(imported_contrasts)))
                    message(sprintf("DEBUG66: Contrast strings: %s", 
                        paste(imported_contrasts$contrasts, collapse = ", ")))
                } else {
                    message("DEBUG66: No contrasts loaded (imported_contrasts is NULL)")
                }

                # --- 6. Update workflow_data ---
                workflow_data$design_matrix <- imported_design
                workflow_data$contrasts_tbl <- imported_contrasts
                workflow_data$data_tbl <- assay_list
                workflow_data$data_cln <- assay_list

                if (!is.null(imported_contrasts)) {
                    assign("contrasts_tbl", imported_contrasts, envir = .GlobalEnv)
                    logger::log_info("Saved contrasts_tbl to global environment.")
                }

                # Create tech_rep_group column
                workflow_data$design_matrix <- workflow_data$design_matrix |>
                    dplyr::mutate(tech_rep_group = paste(group, replicates, sep = "_"))

                # --- 7. Create S4 Object ---
                col_map <- workflow_data$column_mapping

                logger::log_info("Creating MetaboliteAssayData S4 object from imported data")
                
                # DEBUG66: Log S4 creation parameters
                message("DEBUG66: === S4 Object Creation Parameters ===")
                message(sprintf("DEBUG66: assay_list names: %s", paste(names(assay_list), collapse = ", ")))
                message(sprintf("DEBUG66: design_matrix dims: %d x %d", 
                    nrow(workflow_data$design_matrix), ncol(workflow_data$design_matrix)))
                message(sprintf("DEBUG66: design_matrix$Run: %s", 
                    paste(head(workflow_data$design_matrix$Run, 5), collapse = ", ")))
                message(sprintf("DEBUG66: col_map is NULL: %s", is.null(col_map)))
                if (!is.null(col_map)) {
                    message(sprintf("DEBUG66: col_map$metabolite_id_col: '%s'", col_map$metabolite_id_col))
                    message(sprintf("DEBUG66: col_map$annotation_col: '%s'", col_map$annotation_col))
                }
                message("DEBUG66: Calling createMetaboliteAssayData()...")

                s4_obj <- createMetaboliteAssayData(
                    metabolite_data = assay_list
                    , design_matrix = workflow_data$design_matrix
                    , metabolite_id_column = col_map$metabolite_id_col
                    , annotation_id_column = if (!is.null(col_map$annotation_col) && !is.na(col_map$annotation_col) && nzchar(col_map$annotation_col)) {
                        col_map$annotation_col
                    } else {
                        NA_character_
                    }
                    , sample_id = "Run"
                    , group_id = "group"
                    , technical_replicate_id = "tech_rep_group"
                    , database_identifier_type = "Unknown"
                    , internal_standard_regex = if (!is.null(col_map$is_pattern) && !is.na(col_map$is_pattern)) {
                        col_map$is_pattern
                    } else {
                        NA_character_
                    }
                    , args = workflow_data$config_list
                )
                
                # DEBUG66: S4 creation successful
                message("DEBUG66: createMetaboliteAssayData() completed successfully")
                message(sprintf("DEBUG66: S4 object class: %s", class(s4_obj)[1]))

                # --- 8. Initialize state manager and save ---
                message("DEBUG66: Initializing state manager...")
                if (is.null(workflow_data$state_manager)) {
                    workflow_data$state_manager <- WorkflowState$new("metabolomics")
                    message("DEBUG66: Created new WorkflowState")
                } else {
                    message("DEBUG66: Using existing WorkflowState")
                }

                message("DEBUG66: Saving state to state_manager...")
                workflow_data$state_manager$saveState(
                    state_name = "metab_raw_data_s4"
                    , s4_data_object = s4_obj
                    , config_object = workflow_data$config_list
                    , description = "MetaboliteAssayData S4 object created from imported design"
                )
                message("DEBUG66: State saved successfully")

                # Initialize QC progress tracking with raw data baseline
                tryCatch({
                    updateMetaboliteFiltering(
                        theObject = s4_obj
                        , step_name = "1_Raw_Data"
                        , omics_type = "metabolomics"
                        , return_grid = FALSE
                        , overwrite = TRUE
                    )
                    logger::log_info("Initialized QC progress tracking with raw data baseline")
                }, error = function(e) {
                    logger::log_warn(paste("Could not initialize QC progress:", e$message))
                })

                logger::log_info("S4 object saved to state manager as 'metab_raw_data_s4'")
                logger::log_info("Import complete - user can proceed to QC")

                # --- 9. Trigger QC ---
                if (!is.null(qc_trigger)) {
                    qc_trigger(TRUE)
                }

                # Must replace entire list to trigger reactivity
                updated_status <- workflow_data$tab_status
                updated_status$design_matrix <- "complete"
                workflow_data$tab_status <- updated_status

                shiny::removeNotification("importing_design")
                shiny::showNotification(
                    sprintf("Design imported successfully! Loaded %d assays with %d samples.",
                        length(assay_list), nrow(imported_design))
                    , type = "message"
                )

            }, error = function(e) {
                msg <- paste("Error during import:", e$message)
                logger::log_error(msg)
                shiny::showNotification(msg, type = "error", duration = 15)
                shiny::removeNotification("importing_design")
            })
        })

        # == Reactivity Checks ======================================================

        output$data_available <- shiny::reactive({
            !is.null(workflow_data$data_tbl) && !is.null(workflow_data$config_list)
        })
        shiny::outputOptions(output, "data_available", suspendWhenHidden = FALSE)

        output$design_matrix_exists <- shiny::reactive({
            !is.null(workflow_data$design_matrix)
        })
        shiny::outputOptions(output, "design_matrix_exists", suspendWhenHidden = FALSE)

        # == Module Integration =====================================================

        if (exists("mod_metab_design_builder_server")) {
            builder_results_rv <- mod_metab_design_builder_server(
                "builder"
                , data_tbl = shiny::reactive(workflow_data$data_tbl)
                , config_list = shiny::reactive(workflow_data$config_list)
                , column_mapping = shiny::reactive(workflow_data$column_mapping)
                , existing_design_matrix = shiny::reactive(workflow_data$design_matrix)
                , existing_contrasts = shiny::reactive(workflow_data$contrasts_tbl)
            )
        } else {
            builder_results_rv <- shiny::reactiveVal(NULL)
        }

        # == Handle Builder Results =================================================

        shiny::observeEvent(builder_results_rv(), {
            results <- builder_results_rv()
            shiny::req(results)

            shiny::showModal(shiny::modalDialog(
                title = "Processing Design Matrix"
                , shiny::div(
                    style = "text-align: center; padding: 20px;"
                    , shiny::icon("spinner", class = "fa-spin fa-3x")
                    , shiny::br()
                    , shiny::br()
                    , shiny::p("Saving design matrix and preparing data...")
                )
                , footer = NULL
                , easyClose = FALSE
            ))

            logger::log_info("Received results from metabolomics design builder. Saving to workflow and disk.")

            # 1. Update workflow_data
            workflow_data$design_matrix <- results$design_matrix
            workflow_data$data_cln <- results$data_cln
            workflow_data$contrasts_tbl <- results$contrasts_tbl
            workflow_data$config_list <- results$config_list

            if (!is.null(results$contrasts_tbl)) {
                assign("contrasts_tbl", results$contrasts_tbl, envir = .GlobalEnv)
                logger::log_info("Updated contrasts_tbl in global environment.")
            }

            assign("config_list", workflow_data$config_list, envir = .GlobalEnv)
            logger::log_info("Updated global config_list.")

            # 2. Get source directory for saving files
            source_dir <- experiment_paths$source_dir
            if (is.null(source_dir) || !dir.exists(source_dir)) {
                msg <- "Could not find source directory to save files."
                logger::log_error(msg)
                shiny::showNotification(msg, type = "error", duration = 15)
                shiny::removeModal()
                return()
            }

            # 3. Save files
            tryCatch({
                # --- Save Design Matrix ---
                design_matrix_path <- file.path(source_dir, "design_matrix.tab")
                logger::log_info(paste("Writing design matrix to:", design_matrix_path))
                utils::write.table(results$design_matrix, file = design_matrix_path
                    , sep = "\t", row.names = FALSE, quote = FALSE)

                # --- Save Contrasts ---
                if (!is.null(results$contrasts_tbl) && nrow(results$contrasts_tbl) > 0) {
                    contrast_path <- file.path(source_dir, "contrast_strings.tab")
                    logger::log_info(paste("Writing contrasts to:", contrast_path))
                    writeLines(results$contrasts_tbl$contrasts, contrast_path)
                }

                # --- Save each assay's cleaned data (CRITICAL for import to work) ---
                assay_names <- names(results$data_cln)
                for (assay_name in assay_names) {
                    assay_path <- file.path(source_dir, paste0("data_cln_", assay_name, ".tab"))
                    logger::log_info(paste("Writing assay data to:", assay_path))
                    utils::write.table(results$data_cln[[assay_name]], file = assay_path
                        , sep = "\t", row.names = FALSE, quote = FALSE)
                }

                # --- Save assay manifest ---
                manifest_path <- file.path(source_dir, "assay_manifest.txt")
                writeLines(assay_names, manifest_path)
                logger::log_info(sprintf("Saved assay manifest with %d assays: %s",
                    length(assay_names), paste(assay_names, collapse = ", ")))

                # --- Save column mapping as JSON ---
                col_map <- workflow_data$column_mapping
                if (!is.null(col_map)) {
                    col_map_path <- file.path(source_dir, "column_mapping.json")
                    jsonlite::write_json(col_map, col_map_path, auto_unbox = TRUE)
                    logger::log_info("Saved column_mapping.json")
                }

                # --- Save config.ini ---
                if (!is.null(workflow_data$config_list)) {
                    config_path <- file.path(source_dir, "config.ini")
                    logger::log_info(paste("Writing config.ini to:", config_path))
                    tryCatch({
                        ini::write.ini(workflow_data$config_list, config_path)
                        logger::log_info("Saved config.ini.")
                    }, error = function(e) {
                        logger::log_warn(paste("Could not save config.ini:", e$message))
                    })
                }

                # --- 4. Create S4 Object ---
                s4_obj <- createMetaboliteAssayData(
                    metabolite_data = results$data_cln
                    , design_matrix = results$design_matrix
                    , metabolite_id_column = col_map$metabolite_id_col
                    , annotation_id_column = if (!is.null(col_map$annotation_col) && !is.na(col_map$annotation_col) && nzchar(col_map$annotation_col)) {
                        col_map$annotation_col
                    } else {
                        NA_character_
                    }
                    , sample_id = "Run"
                    , group_id = "group"
                    , technical_replicate_id = "tech_rep_group"
                    , database_identifier_type = "Unknown"
                    , internal_standard_regex = if (!is.null(col_map$is_pattern) && !is.na(col_map$is_pattern)) {
                        col_map$is_pattern
                    } else {
                        NA_character_
                    }
                    , args = results$config_list
                )

                # Initialize state manager if needed
                if (is.null(workflow_data$state_manager)) {
                    workflow_data$state_manager <- WorkflowState$new("metabolomics")
                }

                logger::log_info("Saving MetaboliteAssayData S4 object to state manager as 'metab_raw_data_s4'")
                workflow_data$state_manager$saveState(
                    state_name = "metab_raw_data_s4"
                    , s4_data_object = s4_obj
                    , config_object = results$config_list
                    , description = "Initial MetaboliteAssayData S4 object created after design matrix"
                )

                # --- 5. Trigger QC ---
                if (!is.null(qc_trigger)) {
                    qc_trigger(TRUE)
                    logger::log_info("QC trigger set to TRUE")
                }

                # --- 6. Update tab status ---
                # Must replace entire list to trigger reactivity
                updated_status <- workflow_data$tab_status
                updated_status$design_matrix <- "complete"
                workflow_data$tab_status <- updated_status

                shiny::removeModal()
                shiny::showNotification("Design matrix and contrasts saved successfully!", type = "message")

                logger::log_info(sprintf("Design save complete. Files saved to: %s", source_dir))
                logger::log_info(sprintf("Saved: design_matrix.tab, %d assay files, assay_manifest.txt, column_mapping.json, config.ini",
                    length(assay_names)))

            }, error = function(e) {
                msg <- paste("Error saving design matrix results:", e$message)
                logger::log_error(msg)
                shiny::removeModal()
                shiny::showNotification(msg, type = "error", duration = 15)
            })

        }, ignoreNULL = TRUE)

        # == Previews of Saved Data =================================================

        output$design_matrix_preview <- DT::renderDT({
            shiny::req(workflow_data$design_matrix)
            workflow_data$design_matrix
        }, options = list(pageLength = 5, scrollX = TRUE))

        output$contrasts_preview <- DT::renderDT({
            shiny::req(workflow_data$contrasts_tbl)
            workflow_data$contrasts_tbl
        }, options = list(pageLength = 5, scrollX = TRUE))

        output$assays_preview <- shiny::renderText({
            shiny::req(workflow_data$data_tbl)
            assay_names <- names(workflow_data$data_tbl)
            paste("Included assays:", paste(assay_names, collapse = ", "))
        })
    })
}
