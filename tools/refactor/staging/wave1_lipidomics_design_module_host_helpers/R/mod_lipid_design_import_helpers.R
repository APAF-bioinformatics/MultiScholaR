initializeLipidDesignImportBootstrap <- function(
    input,
    session,
    experimentPaths,
    volumes = NULL,
    dirChooseFn = shinyFiles::shinyDirChoose,
    isolateFn = shiny::isolate,
    getVolumesFn = shinyFiles::getVolumes,
    dirExistsFn = dir.exists,
    logInfo = logger::log_info
) {
    resolvedVolumes <- isolateFn({
        baseVolumes <- if (is.function(volumes)) {
            volumes()
        } else if (is.null(volumes)) {
            getVolumesFn()()
        } else {
            volumes
        }

        if (!is.null(experimentPaths) &&
            !is.null(experimentPaths$base_dir) &&
            dirExistsFn(experimentPaths$base_dir)) {
            enhancedVolumes <- c("Project Base Dir" = experimentPaths$base_dir, baseVolumes)
            logInfo(paste("Added base_dir to volumes:", experimentPaths$base_dir))
            enhancedVolumes
        } else {
            baseVolumes
        }
    })

    dirChooseFn(input, "import_dir", roots = resolvedVolumes, session = session)

    list(
        resolvedVolumes = resolvedVolumes
    )
}

buildLipidDesignImportModal <- function(
    ns,
    modalDialogFn = shiny::modalDialog,
    paragraphFn = shiny::p,
    helpTextFn = shiny::helpText,
    dirButtonFn = shinyFiles::shinyDirButton,
    verbatimTextOutputFn = shiny::verbatimTextOutput,
    tagListFn = shiny::tagList,
    modalButtonFn = shiny::modalButton,
    actionButtonFn = shiny::actionButton
) {
    modalDialogFn(
        title = "Import Existing Design Matrix",
        paragraphFn("Select the folder containing 'design_matrix.tab', assay data files, and optionally 'contrast_strings.tab'."),
        helpTextFn("Required files: design_matrix.tab, assay_manifest.txt, data_cln_*.tab files"),
        dirButtonFn(ns("import_dir"), "Select Folder", "Choose a directory"),
        verbatimTextOutputFn(ns("import_dir_path"), placeholder = TRUE),
        footer = tagListFn(
            modalButtonFn("Cancel"),
            actionButtonFn(ns("confirm_import"), "Import", class = "btn-primary")
        )
    )
}

registerLipidDesignImportModalShell <- function(
    input,
    output,
    session,
    resolvedVolumes,
    observeEventFn = shiny::observeEvent,
    showModalFn = shiny::showModal,
    renderTextFn = shiny::renderText,
    reqFn = shiny::req,
    parseDirPathFn = shinyFiles::parseDirPath,
    buildImportModal = buildLipidDesignImportModal
) {
    observeEventFn(input$show_import_modal, {
        showModalFn(buildImportModal(session$ns))
    })

    output$import_dir_path <- renderTextFn({
        reqFn(input$import_dir)
        parseDirPathFn(resolvedVolumes, input$import_dir)
    })

    invisible(output)
}

runLipidDesignImportConfirmationShell <- function(
    workflowData,
    experimentPaths,
    importPath,
    qcTrigger = NULL,
    removeModalFn = shiny::removeModal,
    removeNotificationFn = shiny::removeNotification,
    showNotificationFn = shiny::showNotification,
    fileExistsFn = file.exists,
    readConfigFileFn = readConfigFile,
    vroomFn = vroom::vroom,
    readLinesFn = readLines,
    readJsonFn = jsonlite::read_json,
    assignFn = assign,
    createLipidomicsAssayDataFn = createLipidomicsAssayData,
    workflowStateClass = WorkflowState,
    updateLipidFilteringFn = updateLipidFiltering
) {
    removeModalFn()

    design_file <- file.path(importPath, "design_matrix.tab")
    contrast_file <- file.path(importPath, "contrast_strings.tab")
    manifest_file <- file.path(importPath, "assay_manifest.txt")
    col_map_file <- file.path(importPath, "column_mapping.json")
    config_file <- file.path(importPath, "config.ini")

    if (!fileExistsFn(design_file)) {
        msg <- "Import failed: 'design_matrix.tab' not found in the selected directory."
        logger::log_error(msg)
        showNotificationFn(msg, type = "error", duration = 10)
        return(invisible(NULL))
    }

    if (!fileExistsFn(manifest_file)) {
        msg <- "Import failed: 'assay_manifest.txt' not found. Cannot determine which assay files to load."
        logger::log_error(msg)
        showNotificationFn(msg, type = "error", duration = 10)
        return(invisible(NULL))
    }

    showNotificationFn("Importing design files...", id = "importing_design", duration = NULL)

    tryCatch({
        # DEBUG66: Entry point
        message("DEBUG66: ========================================")
        message("DEBUG66: STARTING IMPORT EXISTING DESIGN")
        message(sprintf("DEBUG66: import_path = %s", importPath))
        message(sprintf("DEBUG66: config_file exists = %s", fileExistsFn(config_file)))
        message(sprintf("DEBUG66: design_file exists = %s", fileExistsFn(design_file)))
        message(sprintf("DEBUG66: manifest_file exists = %s", fileExistsFn(manifest_file)))
        message("DEBUG66: ========================================")

        # --- 1. Load config.ini ---
        if (fileExistsFn(config_file)) {
            logger::log_info("Loading config.ini from import directory.")
            message("DEBUG66: About to call readConfigFile()...")
            workflowData$config_list <- readConfigFileFn(file = config_file)
            message("DEBUG66: readConfigFile() completed successfully")
            message(sprintf("DEBUG66: config_list names: %s", paste(names(workflowData$config_list), collapse = ", ")))
            assignFn("config_list", workflowData$config_list, envir = .GlobalEnv)
            message("DEBUG66: config assigned to global env")
            logger::log_info("Loaded config.ini and assigned to global environment.")
        } else {
            default_config <- file.path(experimentPaths$source_dir, "config.ini")
            if (fileExistsFn(default_config)) {
                workflowData$config_list <- readConfigFileFn(file = default_config)
                assignFn("config_list", workflowData$config_list, envir = .GlobalEnv)
                logger::log_info("Loaded config.ini from source_dir.")
            } else {
                logger::log_warn("No config.ini found. Using empty config.")
                workflowData$config_list <- list()
            }
        }

        # --- 2. Load design matrix ---
        logger::log_info("Loading design_matrix.tab")
        imported_design <- vroomFn(design_file, show_col_types = FALSE)

        # DEBUG66: Log design matrix details
        message(sprintf("DEBUG66: Loaded design_matrix with %d rows, %d cols",
            nrow(imported_design), ncol(imported_design)))
        message(sprintf("DEBUG66: design_matrix columns: %s",
            paste(names(imported_design), collapse = ", ")))
        message(sprintf("DEBUG66: design_matrix$Run values: %s%s",
            paste(head(imported_design$Run, 5), collapse = ", "),
            if (nrow(imported_design) > 5) "..." else ""))

        # --- 3. Load assay manifest and assay data files ---
        logger::log_info("Loading assay_manifest.txt")
        assay_names <- readLinesFn(manifest_file)
        assay_names <- assay_names[nzchar(assay_names)]

        if (length(assay_names) == 0) {
            stop("assay_manifest.txt is empty. No assays to load.")
        }

        logger::log_info(sprintf("Found %d assays in manifest: %s", length(assay_names), paste(assay_names, collapse = ", ")))

        assay_list <- list()
        for (assay_name in assay_names) {
            assay_file <- file.path(importPath, paste0("data_cln_", assay_name, ".tab"))
            if (!fileExistsFn(assay_file)) {
                stop(sprintf("Assay data file not found: %s", assay_file))
            }
            logger::log_info(sprintf("Loading assay: %s", assay_name))
            assay_list[[assay_name]] <- vroomFn(assay_file, show_col_types = FALSE)
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
        if (fileExistsFn(col_map_file)) {
            logger::log_info("Loading column_mapping.json")
            col_map <- readJsonFn(col_map_file, simplifyVector = TRUE)
            workflowData$column_mapping <- col_map

            # DEBUG66: Log column mapping from file
            message(sprintf("DEBUG66: Loaded column_mapping from JSON"))
            message(sprintf("DEBUG66: lipid_id_col = '%s'", col_map$lipid_id_col))
            message(sprintf("DEBUG66: annotation_col = '%s'", col_map$annotation_col))
            message(sprintf("DEBUG66: sample_columns count = %d", length(col_map$sample_columns)))
            message(sprintf("DEBUG66: is_pattern = '%s'", col_map$is_pattern))
        } else {
            logger::log_warn("column_mapping.json not found. Inferring from data structure.")
            first_assay <- assay_list[[1]]
            all_cols <- names(first_assay)

            id_candidates <- c("Peak ID", "Alignment ID", "Alignment.ID", "peak id", "alignment id")
            lipid_id_col <- NULL
            for (candidate in id_candidates) {
                if (candidate %in% all_cols) {
                    lipid_id_col <- candidate
                    break
                }
            }
            if (is.null(lipid_id_col)) {
                lipid_id_col <- all_cols[1]
            }

            annot_candidates <- c("Lipid name", "Lipid.name", "Name", "name")
            annotation_col <- NULL
            for (candidate in annot_candidates) {
                if (candidate %in% all_cols) {
                    annotation_col <- candidate
                    break
                }
            }

            sample_cols <- intersect(imported_design$Run, all_cols)

            col_map <- list(
                lipid_id_col = lipid_id_col,
                annotation_col = annotation_col,
                sample_columns = sample_cols,
                is_pattern = NA_character_
            )
            workflowData$column_mapping <- col_map
            logger::log_info(sprintf("Inferred column mapping: lipid_id=%s, annotation=%s, %d samples",
                lipid_id_col, annotation_col, length(sample_cols)))
        }

        # --- 5. Load contrasts ---
        imported_contrasts <- if (fileExistsFn(contrast_file)) {
            contrast_strings <- readLinesFn(contrast_file)

            if (length(contrast_strings) > 0) {
                contrast_info <- lapply(contrast_strings, function(contrast_string) {
                    clean_string <- gsub("^group", "", contrast_string)
                    clean_string <- gsub("-group", "-", clean_string)
                    friendly_name <- gsub("-", "_vs_", clean_string)
                    full_format <- paste0(friendly_name, "=", contrast_string)

                    list(
                        contrast_string = contrast_string,
                        friendly_name = friendly_name,
                        full_format = full_format
                    )
                })

                data.frame(
                    contrasts = sapply(contrast_info, function(x) x$contrast_string),
                    friendly_names = sapply(contrast_info, function(x) x$friendly_name),
                    full_format = sapply(contrast_info, function(x) x$full_format),
                    stringsAsFactors = FALSE
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
        workflowData$design_matrix <- imported_design
        workflowData$contrasts_tbl <- imported_contrasts
        workflowData$data_tbl <- assay_list
        workflowData$data_cln <- assay_list

        if (!is.null(imported_contrasts)) {
            assignFn("contrasts_tbl", imported_contrasts, envir = .GlobalEnv)
            logger::log_info("Saved contrasts_tbl to global environment.")
        }

        workflowData$design_matrix <- workflowData$design_matrix |>
            dplyr::mutate(tech_rep_group = paste(group, replicates, sep = "_"))

        # --- 7. Create S4 Object ---
        col_map <- workflowData$column_mapping

        logger::log_info("Creating LipidomicsAssayData S4 object from imported data")

        # DEBUG66: Log S4 creation parameters
        message("DEBUG66: === S4 Object Creation Parameters ===")
        message(sprintf("DEBUG66: assay_list names: %s", paste(names(assay_list), collapse = ", ")))
        message(sprintf("DEBUG66: design_matrix dims: %d x %d",
            nrow(workflowData$design_matrix), ncol(workflowData$design_matrix)))
        message(sprintf("DEBUG66: design_matrix$Run: %s",
            paste(head(workflowData$design_matrix$Run, 5), collapse = ", ")))
        message(sprintf("DEBUG66: col_map is NULL: %s", is.null(col_map)))
        if (!is.null(col_map)) {
            message(sprintf("DEBUG66: col_map$lipid_id_col: '%s'", col_map$lipid_id_col))
            message(sprintf("DEBUG66: col_map$annotation_col: '%s'", col_map$annotation_col))
        }
        message("DEBUG66: Calling createLipidomicsAssayData()...")

        s4_obj <- createLipidomicsAssayDataFn(
            lipid_data = assay_list,
            design_matrix = workflowData$design_matrix,
            lipid_id_column = col_map$lipid_id_col,
            annotation_id_column = if (!is.null(col_map$annotation_col) && !is.na(col_map$annotation_col) && nzchar(col_map$annotation_col)) {
                col_map$annotation_col
            } else {
                NA_character_
            },
            sample_id = "Run",
            group_id = "group",
            technical_replicate_id = "tech_rep_group",
            database_identifier_type = "Unknown",
            internal_standard_regex = if (!is.null(col_map$is_pattern) && !is.na(col_map$is_pattern)) {
                col_map$is_pattern
            } else {
                NA_character_
            },
            args = workflowData$config_list
        )

        # DEBUG66: S4 creation successful
        message("DEBUG66: createLipidomicsAssayData() completed successfully")
        message(sprintf("DEBUG66: S4 object class: %s", class(s4_obj)[1]))

        # --- 8. Initialize state manager and save ---
        message("DEBUG66: Initializing state manager...")
        if (is.null(workflowData$state_manager)) {
            workflowData$state_manager <- workflowStateClass$new("lipidomics")
            message("DEBUG66: Created new WorkflowState")
        } else {
            message("DEBUG66: Using existing WorkflowState")
        }

        message("DEBUG66: Saving state to state_manager...")
        workflowData$state_manager$saveState(
            state_name = "lipid_raw_data_s4",
            s4_data_object = s4_obj,
            config_object = workflowData$config_list,
            description = "LipidomicsAssayData S4 object created from imported design"
        )
        message("DEBUG66: State saved successfully")

        tryCatch({
            updateLipidFilteringFn(
                theObject = s4_obj,
                step_name = "1_Raw_Data",
                omics_type = "lipidomics",
                return_grid = FALSE,
                overwrite = TRUE
            )
            logger::log_info("Initialized QC progress tracking with raw data baseline")
        }, error = function(e) {
            logger::log_warn(paste("Could not initialize QC progress:", e$message))
        })

        logger::log_info("S4 object saved to state manager as 'lipid_raw_data_s4'")
        logger::log_info("Import complete - user can proceed to QC")

        if (!is.null(qcTrigger)) {
            qcTrigger(TRUE)
        }

        updated_status <- workflowData$tab_status
        updated_status$design_matrix <- "complete"
        workflowData$tab_status <- updated_status

        removeNotificationFn("importing_design")
        showNotificationFn(
            sprintf("Design imported successfully! Loaded %d assays with %d samples.",
                length(assay_list), nrow(imported_design)),
            type = "message"
        )
    }, error = function(e) {
        msg <- paste("Error during import:", e$message)
        logger::log_error(msg)
        showNotificationFn(msg, type = "error", duration = 15)
        removeNotificationFn("importing_design")
    })
}

registerLipidDesignImportConfirmationObserver <- function(
    input,
    resolvedVolumes,
    workflowData,
    experimentPaths,
    qcTrigger = NULL,
    observeEventFn = shiny::observeEvent,
    reqFn = shiny::req,
    parseDirPathFn = shinyFiles::parseDirPath,
    runImportConfirmationShell = runLipidDesignImportConfirmationShell
) {
    observeEventFn(input$confirm_import, {
        reqFn(input$import_dir)

        importPath <- parseDirPathFn(resolvedVolumes, input$import_dir)
        reqFn(importPath)

        runImportConfirmationShell(
            workflowData = workflowData,
            experimentPaths = experimentPaths,
            importPath = importPath,
            qcTrigger = qcTrigger
        )
    })
}

