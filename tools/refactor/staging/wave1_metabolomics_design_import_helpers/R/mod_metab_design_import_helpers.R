initializeMetabDesignImportBootstrap <- function(
    input,
    session,
    experimentPaths,
    volumes = NULL,
    isolateFn = shiny::isolate,
    getVolumesFn = shinyFiles::getVolumes,
    dirExistsFn = dir.exists,
    dirChooseFn = shinyFiles::shinyDirChoose,
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

    resolvedVolumes
}

registerMetabDesignImportModalShell <- function(
    input,
    output,
    session,
    resolvedVolumes,
    observeEventFn = shiny::observeEvent,
    showModalFn = shiny::showModal,
    modalDialogFn = shiny::modalDialog,
    paragraphFn = shiny::p,
    helpTextFn = shiny::helpText,
    dirButtonFn = shinyFiles::shinyDirButton,
    verbatimTextOutputFn = shiny::verbatimTextOutput,
    tagListFn = shiny::tagList,
    modalButtonFn = shiny::modalButton,
    actionButtonFn = shiny::actionButton,
    renderTextFn = shiny::renderText,
    reqFn = shiny::req,
    parseDirPathFn = shinyFiles::parseDirPath
) {
    observeEventFn(input$show_import_modal, {
        ns <- session$ns
        showModalFn(modalDialogFn(
            title = "Import Existing Design Matrix"
            , paragraphFn("Select the folder containing 'design_matrix.tab', assay data files, and optionally 'contrast_strings.tab'.")
            , helpTextFn("Required files: design_matrix.tab, assay_manifest.txt, data_cln_*.tab files")
            , dirButtonFn(ns("import_dir"), "Select Folder", "Choose a directory")
            , verbatimTextOutputFn(ns("import_dir_path"), placeholder = TRUE)
            , footer = tagListFn(
                modalButtonFn("Cancel")
                , actionButtonFn(ns("confirm_import"), "Import", class = "btn-primary")
            )
        ))
    })

    output$import_dir_path <- renderTextFn({
        reqFn(input$import_dir)
        parseDirPathFn(resolvedVolumes, input$import_dir)
    })

    invisible(output)
}

resolveMetabDesignImportPreflight <- function(
    input,
    resolvedVolumes,
    reqFn = shiny::req,
    parseDirPathFn = shinyFiles::parseDirPath,
    fileExistsFn = file.exists
) {
    reqFn(input$import_dir)

    importPath <- parseDirPathFn(resolvedVolumes, input$import_dir)
    reqFn(importPath)

    designFile <- file.path(importPath, "design_matrix.tab")
    contrastFile <- file.path(importPath, "contrast_strings.tab")
    manifestFile <- file.path(importPath, "assay_manifest.txt")
    colMapFile <- file.path(importPath, "column_mapping.json")
    configFile <- file.path(importPath, "config.ini")

    if (!fileExistsFn(designFile)) {
        return(list(
            ok = FALSE,
            errorMessage = "Import failed: 'design_matrix.tab' not found in the selected directory."
        ))
    }

    if (!fileExistsFn(manifestFile)) {
        return(list(
            ok = FALSE,
            errorMessage = "Import failed: 'assay_manifest.txt' not found. Cannot determine which assay files to load."
        ))
    }

    list(
        ok = TRUE,
        importPath = importPath,
        designFile = designFile,
        contrastFile = contrastFile,
        manifestFile = manifestFile,
        colMapFile = colMapFile,
        configFile = configFile
    )
}

hydrateMetabDesignImportArtifacts <- function(
    workflowData,
    experimentPaths,
    importPath,
    designFile,
    manifestFile,
    configFile,
    readConfigFn = readConfigFile,
    readTabularFn = vroom::vroom,
    readLinesFn = readLines,
    fileExistsFn = file.exists,
    assignFn = assign,
    logInfo = logger::log_info,
    logWarn = logger::log_warn,
    messageFn = message
) {
    messageFn("DEBUG66: ========================================")
    messageFn("DEBUG66: STARTING IMPORT EXISTING DESIGN")
    messageFn(sprintf("DEBUG66: import_path = %s", importPath))
    messageFn(sprintf("DEBUG66: config_file exists = %s", fileExistsFn(configFile)))
    messageFn(sprintf("DEBUG66: design_file exists = %s", fileExistsFn(designFile)))
    messageFn(sprintf("DEBUG66: manifest_file exists = %s", fileExistsFn(manifestFile)))
    messageFn("DEBUG66: ========================================")

    if (fileExistsFn(configFile)) {
        logInfo("Loading config.ini from import directory.")
        messageFn("DEBUG66: About to call readConfigFile()...")
        workflowData$config_list <- readConfigFn(file = configFile)
        messageFn("DEBUG66: readConfigFile() completed successfully")
        messageFn(sprintf("DEBUG66: config_list names: %s", paste(names(workflowData$config_list), collapse = ", ")))
        assignFn("config_list", workflowData$config_list, envir = .GlobalEnv)
        messageFn("DEBUG66: config assigned to global env")
        logInfo("Loaded config.ini and assigned to global environment.")
    } else {
        defaultConfig <- file.path(experimentPaths$source_dir, "config.ini")
        if (fileExistsFn(defaultConfig)) {
            workflowData$config_list <- readConfigFn(file = defaultConfig)
            assignFn("config_list", workflowData$config_list, envir = .GlobalEnv)
            logInfo("Loaded config.ini from source_dir.")
        } else {
            logWarn("No config.ini found. Using empty config.")
            workflowData$config_list <- list()
        }
    }

    logInfo("Loading design_matrix.tab")
    importedDesign <- readTabularFn(designFile, show_col_types = FALSE)

    messageFn(sprintf("DEBUG66: Loaded design_matrix with %d rows, %d cols",
        nrow(importedDesign), ncol(importedDesign)))
    messageFn(sprintf("DEBUG66: design_matrix columns: %s",
        paste(names(importedDesign), collapse = ", ")))
    messageFn(sprintf("DEBUG66: design_matrix$Run values: %s%s",
        paste(head(importedDesign$Run, 5), collapse = ", "),
        if (nrow(importedDesign) > 5) "..." else ""))

    logInfo("Loading assay_manifest.txt")
    assayNames <- readLinesFn(manifestFile)
    assayNames <- assayNames[nzchar(assayNames)]

    if (length(assayNames) == 0) {
        stop("assay_manifest.txt is empty. No assays to load.")
    }

    logInfo(sprintf("Found %d assays in manifest: %s", length(assayNames), paste(assayNames, collapse = ", ")))

    assayList <- list()
    for (assayName in assayNames) {
        assayFile <- file.path(importPath, paste0("data_cln_", assayName, ".tab"))
        if (!fileExistsFn(assayFile)) {
            stop(sprintf("Assay data file not found: %s", assayFile))
        }
        logInfo(sprintf("Loading assay: %s", assayName))
        assayList[[assayName]] <- readTabularFn(assayFile, show_col_types = FALSE)
    }

    messageFn(sprintf("DEBUG66: Loaded %d assays: %s",
        length(assayList), paste(names(assayList), collapse = ", ")))
    for (assayName in names(assayList)) {
        messageFn(sprintf("DEBUG66: Assay '%s': %d rows, %d cols, cols: %s",
            assayName, nrow(assayList[[assayName]]), ncol(assayList[[assayName]]),
            paste(head(names(assayList[[assayName]]), 5), collapse = ", ")))
    }

    list(
        importedDesign = importedDesign,
        assayNames = assayNames,
        assayList = assayList
    )
}

hydrateMetabDesignImportMetadata <- function(
    workflowData,
    importedDesign,
    assayList,
    colMapFile,
    contrastFile,
    fileExistsFn = file.exists,
    readJsonFn = function(path) jsonlite::read_json(path, simplifyVector = TRUE),
    readLinesFn = readLines,
    logInfo = logger::log_info,
    logWarn = logger::log_warn,
    messageFn = message
) {
    if (fileExistsFn(colMapFile)) {
        logInfo("Loading column_mapping.json")
        colMap <- readJsonFn(colMapFile)
        workflowData$column_mapping <- colMap

        messageFn("DEBUG66: Loaded column_mapping from JSON")
        messageFn(sprintf("DEBUG66: metabolite_id_col = '%s'", colMap$metabolite_id_col))
        messageFn(sprintf("DEBUG66: annotation_col = '%s'", colMap$annotation_col))
        messageFn(sprintf("DEBUG66: sample_columns count = %d", length(colMap$sample_columns)))
        messageFn(sprintf("DEBUG66: is_pattern = '%s'", colMap$is_pattern))
    } else {
        logWarn("column_mapping.json not found. Inferring from data structure.")
        firstAssay <- assayList[[1]]
        allCols <- names(firstAssay)

        idCandidates <- c("Peak ID", "Alignment ID", "Alignment.ID", "peak id", "alignment id")
        metaboliteIdCol <- NULL
        for (candidate in idCandidates) {
            if (candidate %in% allCols) {
                metaboliteIdCol <- candidate
                break
            }
        }
        if (is.null(metaboliteIdCol)) {
            metaboliteIdCol <- allCols[1]
        }

        annotCandidates <- c("Metabolite name", "Metabolite.name", "Name", "name")
        annotationCol <- NULL
        for (candidate in annotCandidates) {
            if (candidate %in% allCols) {
                annotationCol <- candidate
                break
            }
        }

        sampleCols <- intersect(importedDesign$Run, allCols)

        colMap <- list(
            metabolite_id_col = metaboliteIdCol,
            annotation_col = annotationCol,
            sample_columns = sampleCols,
            is_pattern = NA_character_
        )
        workflowData$column_mapping <- colMap
        logInfo(sprintf(
            "Inferred column mapping: metabolite_id=%s, annotation=%s, %d samples",
            metaboliteIdCol,
            annotationCol,
            length(sampleCols)
        ))
    }

    importedContrasts <- if (fileExistsFn(contrastFile)) {
        contrastStrings <- readLinesFn(contrastFile)

        if (length(contrastStrings) > 0) {
            contrastInfo <- lapply(contrastStrings, function(contrastString) {
                cleanString <- gsub("^group", "", contrastString)
                cleanString <- gsub("-group", "-", cleanString)
                friendlyName <- gsub("-", "_vs_", cleanString)
                fullFormat <- paste0(friendlyName, "=", contrastString)

                list(
                    contrast_string = contrastString,
                    friendly_name = friendlyName,
                    full_format = fullFormat
                )
            })

            data.frame(
                contrasts = sapply(contrastInfo, function(x) x$contrast_string),
                friendly_names = sapply(contrastInfo, function(x) x$friendly_name),
                full_format = sapply(contrastInfo, function(x) x$full_format),
                stringsAsFactors = FALSE
            )
        } else {
            NULL
        }
    } else {
        logInfo("No contrast_strings.tab found.")
        NULL
    }

    if (!is.null(importedContrasts)) {
        messageFn(sprintf("DEBUG66: Loaded %d contrasts", nrow(importedContrasts)))
        messageFn(sprintf(
            "DEBUG66: Contrast strings: %s",
            paste(importedContrasts$contrasts, collapse = ", ")
        ))
    } else {
        messageFn("DEBUG66: No contrasts loaded (imported_contrasts is NULL)")
    }

    list(
        columnMapping = workflowData$column_mapping,
        importedContrasts = importedContrasts
    )
}

hydrateMetabDesignImportWorkflowState <- function(
    workflowData,
    importedDesign,
    assayList,
    importedContrasts,
    assignFn = assign,
    mutateFn = dplyr::mutate,
    logInfo = logger::log_info
) {
    workflowData$design_matrix <- importedDesign
    workflowData$contrasts_tbl <- importedContrasts
    workflowData$data_tbl <- assayList
    workflowData$data_cln <- assayList

    if (!is.null(importedContrasts)) {
        assignFn("contrasts_tbl", importedContrasts, envir = .GlobalEnv)
        logInfo("Saved contrasts_tbl to global environment.")
    }

    workflowData$design_matrix <- mutateFn(
        workflowData$design_matrix,
        tech_rep_group = paste(group, replicates, sep = "_")
    )

    list(
        designMatrix = workflowData$design_matrix,
        assayList = workflowData$data_tbl,
        importedContrasts = workflowData$contrasts_tbl
    )
}

createMetabDesignImportedS4Object <- function(
    workflowData,
    assayList,
    colMap,
    createS4Fn = createMetaboliteAssayData,
    logInfo = logger::log_info,
    messageFn = message
) {
    logInfo("Creating MetaboliteAssayData S4 object from imported data")

    messageFn("DEBUG66: === S4 Object Creation Parameters ===")
    messageFn(sprintf("DEBUG66: assay_list names: %s", paste(names(assayList), collapse = ", ")))
    messageFn(sprintf(
        "DEBUG66: design_matrix dims: %d x %d",
        nrow(workflowData$design_matrix),
        ncol(workflowData$design_matrix)
    ))
    messageFn(sprintf(
        "DEBUG66: design_matrix$Run: %s",
        paste(head(workflowData$design_matrix$Run, 5), collapse = ", ")
    ))
    messageFn(sprintf("DEBUG66: col_map is NULL: %s", is.null(colMap)))
    if (!is.null(colMap)) {
        messageFn(sprintf("DEBUG66: col_map$metabolite_id_col: '%s'", colMap$metabolite_id_col))
        messageFn(sprintf("DEBUG66: col_map$annotation_col: '%s'", colMap$annotation_col))
    }
    messageFn("DEBUG66: Calling createMetaboliteAssayData()...")

    s4Obj <- createS4Fn(
        metabolite_data = assayList,
        design_matrix = workflowData$design_matrix,
        metabolite_id_column = colMap$metabolite_id_col,
        annotation_id_column = if (!is.null(colMap$annotation_col) && !is.na(colMap$annotation_col) && nzchar(colMap$annotation_col)) {
            colMap$annotation_col
        } else {
            NA_character_
        },
        sample_id = "Run",
        group_id = "group",
        technical_replicate_id = "tech_rep_group",
        database_identifier_type = "Unknown",
        internal_standard_regex = if (!is.null(colMap$is_pattern) && !is.na(colMap$is_pattern)) {
            colMap$is_pattern
        } else {
            NA_character_
        },
        args = workflowData$config_list
    )

    messageFn("DEBUG66: createMetaboliteAssayData() completed successfully")
    messageFn(sprintf("DEBUG66: S4 object class: %s", class(s4Obj)[1]))

    s4Obj
}

saveMetabDesignImportedS4State <- function(
    workflowData,
    s4Object,
    stateManagerFactory = function(omicsType) WorkflowState$new(omicsType),
    messageFn = message
) {
    messageFn("DEBUG66: Initializing state manager...")
    if (is.null(workflowData$state_manager)) {
        workflowData$state_manager <- stateManagerFactory("metabolomics")
        messageFn("DEBUG66: Created new WorkflowState")
    } else {
        messageFn("DEBUG66: Using existing WorkflowState")
    }

    messageFn("DEBUG66: Saving state to state_manager...")
    workflowData$state_manager$saveState(
        state_name = "metab_raw_data_s4"
        , s4_data_object = s4Object
        , config_object = workflowData$config_list
        , description = "MetaboliteAssayData S4 object created from imported design"
    )
    messageFn("DEBUG66: State saved successfully")

    workflowData$state_manager
}

initializeMetabDesignImportedQcBaseline <- function(
    s4Object,
    updateFilteringFn = updateMetaboliteFiltering,
    logInfo = logger::log_info,
    logWarn = logger::log_warn
) {
    tryCatch({
        updateFilteringFn(
            theObject = s4Object,
            step_name = "1_Raw_Data",
            omics_type = "metabolomics",
            return_grid = FALSE,
            overwrite = TRUE
        )
        logInfo("Initialized QC progress tracking with raw data baseline")
        TRUE
    }, error = function(e) {
        logWarn(paste("Could not initialize QC progress:", e$message))
        FALSE
    })
}

completeMetabDesignImportedPostCheckpoint <- function(
    workflowData,
    assayNames,
    importedDesign,
    qcTrigger = NULL,
    successNotificationId = "importing_design",
    logInfo = logger::log_info,
    removeNotification = shiny::removeNotification,
    showNotification = shiny::showNotification
) {
    logInfo("S4 object saved to state manager as 'metab_raw_data_s4'")
    logInfo("Import complete - user can proceed to QC")

    if (!is.null(qcTrigger)) {
        qcTrigger(TRUE)
    }

    updatedStatus <- workflowData$tab_status
    updatedStatus$design_matrix <- "complete"
    workflowData$tab_status <- updatedStatus

    if (!is.null(successNotificationId)) {
        removeNotification(successNotificationId)
    }

    successMessage <- sprintf(
        "Design imported successfully! Loaded %d assays with %d samples.",
        length(assayNames),
        nrow(importedDesign)
    )
    showNotification(successMessage, type = "message")

    invisible(list(
        tabStatus = workflowData$tab_status,
        successMessage = successMessage
    ))
}

registerMetabDesignImportObserverShell <- function(
    input,
    resolvedVolumes,
    workflowData,
    experimentPaths,
    qcTrigger = NULL,
    observeEventFn = shiny::observeEvent,
    resolveImportPreflight = resolveMetabDesignImportPreflight,
    hydrateImportArtifacts = hydrateMetabDesignImportArtifacts,
    hydrateImportMetadata = hydrateMetabDesignImportMetadata,
    hydrateImportWorkflowState = hydrateMetabDesignImportWorkflowState,
    createImportedS4Object = createMetabDesignImportedS4Object,
    saveImportedS4State = saveMetabDesignImportedS4State,
    initializeImportedQcBaseline = initializeMetabDesignImportedQcBaseline,
    completeImportedPostCheckpoint = completeMetabDesignImportedPostCheckpoint,
    removeModal = shiny::removeModal,
    showNotification = shiny::showNotification,
    removeNotification = shiny::removeNotification,
    logError = logger::log_error
) {
    observeEventFn(input$confirm_import, {
        importPreflight <- resolveImportPreflight(
            input = input,
            resolvedVolumes = resolvedVolumes
        )

        if (!isTRUE(importPreflight$ok)) {
            msg <- importPreflight$errorMessage
            logError(msg)
            showNotification(msg, type = "error", duration = 10)
            return(invisible(NULL))
        }

        removeModal()

        showNotification("Importing design files...", id = "importing_design", duration = NULL)

        tryCatch({
            importedArtifacts <- hydrateImportArtifacts(
                workflowData = workflowData,
                experimentPaths = experimentPaths,
                importPath = importPreflight$importPath,
                designFile = importPreflight$designFile,
                manifestFile = importPreflight$manifestFile,
                configFile = importPreflight$configFile
            )

            importedDesign <- importedArtifacts$importedDesign
            assayNames <- importedArtifacts$assayNames
            assayList <- importedArtifacts$assayList

            importedMetadata <- hydrateImportMetadata(
                workflowData = workflowData,
                importedDesign = importedDesign,
                assayList = assayList,
                colMapFile = importPreflight$colMapFile,
                contrastFile = importPreflight$contrastFile
            )

            hydrateImportWorkflowState(
                workflowData = workflowData,
                importedDesign = importedDesign,
                assayList = assayList,
                importedContrasts = importedMetadata$importedContrasts
            )

            s4Object <- createImportedS4Object(
                workflowData = workflowData,
                assayList = assayList,
                colMap = importedMetadata$columnMapping
            )

            saveImportedS4State(
                workflowData = workflowData,
                s4Object = s4Object
            )

            initializeImportedQcBaseline(
                s4Object = s4Object
            )

            completeImportedPostCheckpoint(
                workflowData = workflowData,
                assayNames = assayNames,
                importedDesign = importedDesign,
                qcTrigger = qcTrigger
            )
        }, error = function(e) {
            msg <- paste("Error during import:", e$message)
            logError(msg)
            showNotification(msg, type = "error", duration = 15)
            removeNotification("importing_design")
            invisible(NULL)
        })
    })
}

