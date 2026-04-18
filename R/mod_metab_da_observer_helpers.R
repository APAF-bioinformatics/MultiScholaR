runMetabDaSaveHeatmapObserverShell <- function(
    currentHeatmapPlot,
    currentRowClusters = NULL,
    publicationGraphsDir,
    selectedContrast,
    topN = 50,
    clusteringMethod = "complete",
    distanceMethod = "euclidean",
    heatmapClustering = "both",
    scaling = "row",
    colorScheme = "RdBu",
    treeCutMethod = "none",
    nClusters = NULL,
    cutHeight = NULL,
    minClusterSize = 1,
    saveHeatmapProducts = save_heatmap_products,
    showNotification = shiny::showNotification,
    logInfo = logger::log_info
) {
    logInfo("Save Heatmap button clicked")

    params <- list(
        contrast = selectedContrast,
        top_n = topN,
        clustering_method = clusteringMethod,
        distance_method = distanceMethod,
        cluster_rows = heatmapClustering %in% c("both", "row"),
        cluster_cols = heatmapClustering %in% c("both", "column"),
        scaling = scaling,
        color_scheme = colorScheme,
        tree_cut_method = treeCutMethod,
        n_clusters = nClusters,
        cut_height = cutHeight,
        min_cluster_size = minClusterSize
    )

    filePrefix <- paste0("metab_", selectedContrast)
    filePrefix <- gsub("[^A-Za-z0-9_]", "_", filePrefix)

    saveHeatmapProducts(
        heatmap_obj = currentHeatmapPlot,
        row_clusters = currentRowClusters,
        params_list = params,
        output_dir = publicationGraphsDir,
        file_prefix = filePrefix
    )

    showNotification(
        "Heatmap and cluster info saved to publication_graphs/Heatmap",
        type = "message",
        duration = 5
    )

    invisible(filePrefix)
}

runMetabDaSaveHeatmapObserverEntry <- function(
    currentHeatmapPlot,
    publicationGraphsDir,
    currentRowClusters = NULL,
    selectedContrast,
    topN = 50,
    clusteringMethod = "complete",
    distanceMethod = "euclidean",
    heatmapClustering = "both",
    scaling = "row",
    colorScheme = "RdBu",
    treeCutMethod = "none",
    nClusters = NULL,
    cutHeight = NULL,
    minClusterSize = 1,
    requireInputs = shiny::req,
    saveHeatmapShell = runMetabDaSaveHeatmapObserverShell
) {
    requireInputs(currentHeatmapPlot, publicationGraphsDir)

    saveState <- saveHeatmapShell(
        currentHeatmapPlot = currentHeatmapPlot,
        currentRowClusters = currentRowClusters,
        publicationGraphsDir = publicationGraphsDir,
        selectedContrast = selectedContrast,
        topN = topN,
        clusteringMethod = clusteringMethod,
        distanceMethod = distanceMethod,
        heatmapClustering = heatmapClustering,
        scaling = scaling,
        colorScheme = colorScheme,
        treeCutMethod = treeCutMethod,
        nClusters = nClusters,
        cutHeight = cutHeight,
        minClusterSize = minClusterSize
    )

    invisible(list(
        status = "success",
        saveState = saveState
    ))
}

updateMetabDaResultsSelectorInputs <- function(
    daResultsLong,
    session,
    updateSelectInput = shiny::updateSelectInput,
    logInfo = logger::log_info
) {
    if (is.null(daResultsLong) || nrow(daResultsLong) == 0) {
        return(invisible(NULL))
    }

    contrastChoices <- character()
    if ("friendly_name" %in% colnames(daResultsLong)) {
        contrastChoices <- unique(daResultsLong$friendly_name)
    }
    if (length(contrastChoices) == 0 && "comparison" %in% colnames(daResultsLong)) {
        contrastChoices <- unique(daResultsLong$comparison)
    }

    assayValues <- character()
    if ("assay" %in% colnames(daResultsLong)) {
        assayValues <- unique(daResultsLong$assay)
    }

    if (length(contrastChoices) > 0) {
        updateSelectInput(session, "volcano_contrast",
            choices = contrastChoices, selected = contrastChoices[1]
        )
        updateSelectInput(session, "heatmap_contrast",
            choices = contrastChoices, selected = contrastChoices[1]
        )
        updateSelectInput(session, "table_contrast",
            choices = contrastChoices, selected = contrastChoices[1]
        )
    }

    if (length(assayValues) > 0) {
        assayChoices <- c("Combined", assayValues)
        updateSelectInput(session, "volcano_assay",
            choices = assayChoices, selected = "Combined"
        )
        updateSelectInput(session, "heatmap_assay",
            choices = assayChoices, selected = "Combined"
        )

        tableAssayChoices <- c("All", assayValues)
        updateSelectInput(session, "table_assay",
            choices = tableAssayChoices, selected = "All"
        )
    } else {
        assayChoices <- "Combined"
        tableAssayChoices <- "All"
    }

    logInfo(sprintf(
        "   Updated dropdowns: %d contrasts, %d assays",
        length(contrastChoices), length(assayChoices) - 1
    ))

    invisible(list(
        contrastChoices = contrastChoices,
        assayChoices = assayChoices,
        tableAssayChoices = tableAssayChoices
    ))
}

writeMetabDaResultsArtifacts <- function(
    results,
    experimentPaths,
    daQValThresh = 0.05,
    treatLfcCutoff = 0,
    outputResults = outputMetabDaResultsAllContrasts,
    showNotification = shiny::showNotification,
    logInfo = logger::log_info,
    logWarn = logger::log_warn
) {
    daOutputDir <- experimentPaths$da_output_dir
    publicationGraphsDir <- experimentPaths$publication_graphs_dir
    formatPathValue <- function(pathValue) {
        if (is.null(pathValue)) {
            return("NULL")
        }

        pathValue
    }

    logInfo("   Writing DA results to disk...")
    logInfo(sprintf("   da_output_dir = %s", formatPathValue(daOutputDir)))
    logInfo(sprintf(
        "   publication_graphs_dir = %s",
        formatPathValue(publicationGraphsDir)
    ))

    if (is.null(daOutputDir) || is.null(publicationGraphsDir)) {
        logWarn("   Output directories not configured, skipping file output")
        return(invisible(list(
            status = "skipped",
            daOutputDir = daOutputDir,
            publicationGraphsDir = publicationGraphsDir
        )))
    }

    tryCatch(
        {
            success <- outputResults(
                da_results_list = results,
                da_output_dir = daOutputDir,
                publication_graphs_dir = publicationGraphsDir,
                da_q_val_thresh = daQValThresh,
                lfc_threshold = treatLfcCutoff,
                heatmap_top_n = 50,
                heatmap_clustering = "both",
                heatmap_color_scheme = "RdBu"
            )

            if (isTRUE(success)) {
                logInfo("   All DA results written to disk successfully")
                showNotification(
                    "DA results saved to disk (tables, volcano plots, heatmaps)",
                    type = "message",
                    duration = 5
                )
            }

            invisible(list(
                status = if (isTRUE(success)) "written" else "attempted",
                success = success,
                daOutputDir = daOutputDir,
                publicationGraphsDir = publicationGraphsDir
            ))
        },
        error = function(e) {
            logWarn(paste("   Could not write DA results to disk:", e$message))
            showNotification(
                paste("Warning: Could not save results to disk:", e$message),
                type = "warning",
                duration = 8
            )

            invisible(list(
                status = "error",
                errorMessage = e$message,
                daOutputDir = daOutputDir,
                publicationGraphsDir = publicationGraphsDir
            ))
        }
    )
}

writeMetabDaResultsDownloadCsv <- function(
    results,
    file,
    req = shiny::req,
    writeCsv = utils::write.csv
) {
    req(results)
    writeCsv(results, file, row.names = FALSE)

    invisible(file)
}

buildMetabDaResultsDownloadHandler <- function(
    daData,
    downloadHandler = shiny::downloadHandler,
    buildFilename = function() {
        paste0("metabolomics_da_results_", Sys.Date(), ".csv")
    },
    writeDownloadCsv = writeMetabDaResultsDownloadCsv
) {
    downloadHandler(
        filename = function() {
            buildFilename()
        },
        content = function(file) {
            results <- daData$da_results_list$da_metabolites_long
            writeDownloadCsv(
                results = results,
                file = file
            )
        }
    )
}

resolveMetabDaAnalysisInputs <- function(
    currentS4Object = NULL,
    workflowData,
    assayDataClass = "MetaboliteAssayData",
    getState = NULL
) {
    currentS4 <- currentS4Object
    if (is.null(currentS4) && is.function(getState)) {
        currentS4 <- getState()
    }

    if (is.null(currentS4) || !inherits(currentS4, assayDataClass)) {
        return(list(
            ok = FALSE,
            currentS4 = currentS4,
            contrastsTbl = NULL,
            errorMessage = "No metabolomics data loaded. Please load a filtered session first."
        ))
    }

    contrastsTbl <- workflowData$contrasts_tbl
    if (is.null(contrastsTbl) || nrow(contrastsTbl) == 0) {
        return(list(
            ok = FALSE,
            currentS4 = currentS4,
            contrastsTbl = contrastsTbl,
            errorMessage = "No contrasts defined. Please define contrasts in the design tab."
        ))
    }

    list(
        ok = TRUE,
        currentS4 = currentS4,
        contrastsTbl = contrastsTbl,
        errorMessage = NULL
    )
}

runMetabDaAnalysisObserverEntry <- function(
    currentS4Object = NULL,
    workflowData,
    formulaString,
    daQValThresh,
    treatLfcCutoff,
    daData,
    session,
    experimentPaths,
    resolveAnalysisInputs = resolveMetabDaAnalysisInputs,
    runAnalysisShell = runMetabDaAnalysisObserverShell,
    showNotification = shiny::showNotification
) {
    analysisInputs <- resolveAnalysisInputs(
        currentS4Object = currentS4Object,
        workflowData = workflowData,
        getState = if (!is.null(workflowData$state_manager)) {
            workflowData$state_manager$getState
        } else {
            NULL
        }
    )

    if (!isTRUE(analysisInputs$ok)) {
        showNotification(
            analysisInputs$errorMessage,
            type = "error",
            duration = 5
        )

        return(invisible(list(
            status = "error",
            stage = "resolve_analysis_inputs",
            analysisInputs = analysisInputs,
            errorMessage = analysisInputs$errorMessage
        )))
    }

    shellState <- runAnalysisShell(
        currentS4 = analysisInputs$currentS4,
        contrastsTbl = analysisInputs$contrastsTbl,
        formulaString = formulaString,
        daQValThresh = daQValThresh,
        treatLfcCutoff = treatLfcCutoff,
        workflowData = workflowData,
        daData = daData,
        session = session,
        experimentPaths = experimentPaths
    )

    invisible(list(
        status = shellState$status,
        analysisInputs = analysisInputs,
        shellState = shellState
    ))
}

resolveMetabDaLoadSessionFile <- function(
    experimentPaths,
    sessionFilename = "metab_filtered_session_data_latest.rds",
    dirExists = dir.exists,
    fileExists = file.exists,
    filePath = file.path,
    debugLog = NULL
) {
    if (is.null(debugLog)) {
        debugLog <- function(...) invisible(NULL)
    }

    debugLog("  STEP 1: Getting source_dir from experiment_paths")
    debugLog("    experiment_paths class: ", class(experimentPaths)[1])
    debugLog("    experiment_paths is NULL: ", is.null(experimentPaths))

    sourceDir <- experimentPaths$source_dir
    debugLog("  STEP 2: source_dir = ", if (is.null(sourceDir)) "NULL" else sourceDir)

    if (is.null(sourceDir) || !dirExists(sourceDir)) {
        debugLog("  BRANCH: source_dir NULL or not exists, trying export_dir")
        sourceDir <- experimentPaths$export_dir
        debugLog("  STEP 3: export_dir = ", if (is.null(sourceDir)) "NULL" else sourceDir)
    }

    if (is.null(sourceDir) || !dirExists(sourceDir)) {
        debugLog("  ERROR: No valid source directory found")
        return(list(
            ok = FALSE,
            sourceDir = sourceDir,
            sessionFile = NULL,
            errorMessage = "Could not find source directory for session data."
        ))
    }

    sessionFile <- filePath(sourceDir, sessionFilename)
    debugLog("  STEP 4: session_file = ", sessionFile)
    debugLog("  STEP 4: file.exists = ", fileExists(sessionFile))

    if (!fileExists(sessionFile)) {
        return(list(
            ok = FALSE,
            sourceDir = sourceDir,
            sessionFile = sessionFile,
            errorMessage = sprintf("Session file not found: %s", sessionFile)
        ))
    }

    list(
        ok = TRUE,
        sourceDir = sourceDir,
        sessionFile = sessionFile,
        errorMessage = NULL
    )
}

runMetabDaLoadSessionObserverEntry <- function(
    experimentPaths,
    workflowData,
    daData,
    session,
    resolveSessionFile = resolveMetabDaLoadSessionFile,
    loadSessionShell = runMetabDaLoadSessionObserverShell,
    showNotification = shiny::showNotification,
    debugLog = NULL,
    startTime = NULL
) {
    if (is.null(debugLog)) {
        debugLog <- function(...) invisible(NULL)
    }

    sessionResolution <- resolveSessionFile(
        experimentPaths = experimentPaths,
        debugLog = debugLog
    )

    if (!isTRUE(sessionResolution$ok)) {
        showNotification(
            sessionResolution$errorMessage,
            type = "error",
            duration = 5
        )

        return(invisible(list(
            status = "error",
            stage = "resolve_session_file",
            sessionResolution = sessionResolution,
            errorMessage = sessionResolution$errorMessage
        )))
    }

    loadState <- loadSessionShell(
        sessionFile = sessionResolution$sessionFile,
        workflowData = workflowData,
        daData = daData,
        session = session,
        debugLog = debugLog,
        startTime = startTime
    )

    invisible(list(
        status = loadState$status,
        sessionResolution = sessionResolution,
        loadState = loadState
    ))
}

runMetabDaLoadSessionObserverShell <- function(
    sessionFile,
    workflowData,
    daData,
    session,
    readSessionData = readRDS,
    restoreState = restoreMetabDaLoadedSessionState,
    showNotification = shiny::showNotification,
    removeNotification = shiny::removeNotification,
    logInfo = logger::log_info,
    logError = logger::log_error,
    debugLog = NULL,
    startTime = NULL,
    timeNow = Sys.time
) {
    if (is.null(debugLog)) {
        debugLog <- function(...) invisible(NULL)
    }

    showNotification(
        "Loading filtered session...",
        id = "loading_session",
        duration = NULL
    )

    tryCatch(
        {
            debugLog("  STEP 5: Reading RDS file...")
            sessionData <- readSessionData(sessionFile)
            debugLog("  STEP 5: RDS loaded successfully")
            debugLog("    names(session_data): ", paste(names(sessionData), collapse = ", "))

            logInfo(sprintf("   Loaded session from: %s", sessionFile))

            restoreStateResult <- restoreState(
                sessionData = sessionData,
                sessionFile = sessionFile,
                workflowData = workflowData,
                daData = daData,
                session = session,
                debugLog = debugLog
            )

            debugLog("  STEP 10: Removing notification and showing success...")
            removeNotification("loading_session")
            showNotification("Session loaded successfully!", type = "message", duration = 3)

            if (!is.null(startTime)) {
                debugLog(
                    "--- EXIT load_filtered_session (",
                    round(difftime(timeNow(), startTime, units = "secs"), 3), "s) --- SUCCESS"
                )
            }

            invisible(list(
                status = "success",
                sessionData = sessionData,
                restoreState = restoreStateResult
            ))
        },
        error = function(e) {
            debugLog("  FATAL ERROR: ", e$message)
            debugLog("  ERROR CALL: ", deparse(e$call))
            debugLog("--- EXIT load_filtered_session --- FAILED")
            logError(sprintf("   Error loading session: %s", e$message))
            removeNotification("loading_session")
            showNotification(
                sprintf("Error loading session: %s", e$message),
                type = "error",
                duration = 10
            )

            invisible(list(
                status = "error",
                errorMessage = e$message
            ))
        }
    )
}

runMetabDaAnalysisObserverShell <- function(
    currentS4,
    contrastsTbl,
    formulaString,
    daQValThresh,
    treatLfcCutoff,
    workflowData,
    daData,
    session,
    experimentPaths,
    runAnalysis = runMetabolitesDA,
    updateSelectors = updateMetabDaResultsSelectorInputs,
    writeArtifacts = writeMetabDaResultsArtifacts,
    showNotification = shiny::showNotification,
    removeNotification = shiny::removeNotification,
    logInfo = logger::log_info,
    logError = logger::log_error
) {
    showNotification(
        "Running differential expression analysis...",
        id = "da_running",
        duration = NULL
    )

    tryCatch(
        {
            results <- runAnalysis(
                theObject = currentS4,
                contrasts_tbl = contrastsTbl,
                formula_string = formulaString,
                da_q_val_thresh = daQValThresh,
                treat_lfc_cutoff = treatLfcCutoff,
                eBayes_trend = TRUE,
                eBayes_robust = TRUE
            )

            daData$da_results_list <- results
            daData$analysis_complete <- TRUE

            updatedStatus <- workflowData$tab_status
            updatedStatus$differential_analysis <- "complete"
            workflowData$tab_status <- updatedStatus

            removeNotification("da_running")
            showNotification(
                "Differential expression analysis complete!",
                type = "message",
                duration = 5
            )

            logInfo("   DA analysis completed successfully")

            selectorState <- updateSelectors(
                daResultsLong = results$da_metabolites_long,
                session = session
            )

            diskState <- writeArtifacts(
                results = results,
                experimentPaths = experimentPaths,
                daQValThresh = daQValThresh,
                treatLfcCutoff = treatLfcCutoff
            )

            invisible(list(
                status = "success",
                results = results,
                selectorState = selectorState,
                diskState = diskState
            ))
        },
        error = function(e) {
            logError(sprintf("   DA analysis error: %s", e$message))
            removeNotification("da_running")
            showNotification(
                sprintf("Analysis error: %s", e$message),
                type = "error",
                duration = 10
            )

            invisible(list(
                status = "error",
                errorMessage = e$message
            ))
        }
    )
}

restoreMetabDaLoadedSessionState <- function(
    sessionData,
    sessionFile,
    workflowData,
    daData,
    session,
    updateSelectInput = shiny::updateSelectInput,
    updateTextAreaInput = shiny::updateTextAreaInput,
    logInfo = logger::log_info,
    logWarn = logger::log_warn,
    debugLog = NULL
) {
    if (is.null(debugLog)) {
        debugLog <- function(...) invisible(NULL)
    }

    stateName <- NULL
    contrastChoices <- character()
    assayChoices <- character()
    formulaValue <- NULL

    debugLog("  STEP 6: Checking current_s4_object...")
    debugLog(
        "    is.null(sessionData$current_s4_object): ",
        is.null(sessionData$current_s4_object)
    )

    if (!is.null(sessionData$current_s4_object)) {
        debugLog("  STEP 6a: S4 object class: ", class(sessionData$current_s4_object)[1])
        debugLog("  STEP 6b: Assigning to da_data$current_s4_object...")
        daData$current_s4_object <- sessionData$current_s4_object

        debugLog("  STEP 6c: Getting state_name...")
        stateName <- sessionData$r6_current_state_name
        debugLog("    r6_current_state_name: ", if (is.null(stateName)) "NULL" else stateName)
        if (is.null(stateName)) {
            stateName <- "loaded_for_de"
        }

        debugLog("  STEP 6d: Checking workflow_data$state_manager...")
        debugLog("    workflow_data is NULL: ", is.null(workflowData))
        debugLog("    workflow_data class: ", class(workflowData)[1])
        debugLog("    state_manager exists: ", !is.null(workflowData$state_manager))
        if (!is.null(workflowData$state_manager)) {
            debugLog("    state_manager class: ", class(workflowData$state_manager)[1])
        }

        debugLog("  STEP 6e: Calling workflow_data$state_manager$saveState...")
        workflowData$state_manager$saveState(
            state_name = stateName,
            s4_data_object = sessionData$current_s4_object,
            config_object = list(loaded_from = sessionFile),
            description = "Loaded from filtered session for DA analysis"
        )
        debugLog("  STEP 6f: saveState completed")
        logInfo("   Restored S4 object to state manager")
    }

    debugLog("  STEP 7: Checking contrasts_tbl...")
    debugLog("    is.null(sessionData$contrasts_tbl): ", is.null(sessionData$contrasts_tbl))
    if (!is.null(sessionData$contrasts_tbl)) {
        debugLog("    nrow(contrasts_tbl): ", nrow(sessionData$contrasts_tbl))
        debugLog("    colnames: ", paste(colnames(sessionData$contrasts_tbl), collapse = ", "))
    }

    if (!is.null(sessionData$contrasts_tbl) && nrow(sessionData$contrasts_tbl) > 0) {
        debugLog("  STEP 7a: Assigning contrasts to workflow_data...")
        workflowData$contrasts_tbl <- sessionData$contrasts_tbl
        daData$contrasts_available <- sessionData$contrasts_tbl

        debugLog("  STEP 7b: Getting contrast_choices...")
        if ("friendly_names" %in% colnames(sessionData$contrasts_tbl)) {
            contrastChoices <- sessionData$contrasts_tbl$friendly_names
        } else {
            contrastChoices <- sessionData$contrasts_tbl$contrasts
        }
        debugLog("    contrast_choices: ", paste(contrastChoices, collapse = ", "))

        debugLog("  STEP 7c: Updating volcano_contrast dropdown...")
        updateSelectInput(session, "volcano_contrast",
            choices = contrastChoices, selected = contrastChoices[1]
        )
        debugLog("  STEP 7d: Updating heatmap_contrast dropdown...")
        updateSelectInput(session, "heatmap_contrast",
            choices = contrastChoices, selected = contrastChoices[1]
        )
        debugLog("  STEP 7e: Updating table_contrast dropdown...")
        updateSelectInput(session, "table_contrast",
            choices = contrastChoices, selected = contrastChoices[1]
        )

        logInfo(sprintf("   Restored %d contrasts", nrow(sessionData$contrasts_tbl)))
    }

    debugLog("  STEP 8: Checking assay_names...")
    debugLog("    is.null(sessionData$assay_names): ", is.null(sessionData$assay_names))
    if (!is.null(sessionData$assay_names)) {
        debugLog("    assay_names: ", paste(sessionData$assay_names, collapse = ", "))
    }

    if (!is.null(sessionData$assay_names)) {
        daData$assays_available <- sessionData$assay_names
        assayChoices <- c("Combined", sessionData$assay_names)
        debugLog("  STEP 8a: Updating volcano_assay dropdown...")
        updateSelectInput(session, "volcano_assay", choices = assayChoices)
        debugLog("  STEP 8b: Updating heatmap_assay dropdown...")
        updateSelectInput(session, "heatmap_assay", choices = assayChoices)
        debugLog("  STEP 8c: Updating table_assay dropdown...")
        updateSelectInput(session, "table_assay", choices = c("All", sessionData$assay_names))
    }

    debugLog("  STEP 9: Attempting to extract formula from S4 args...")
    tryCatch(
        {
            debugLog("  STEP 9a: Checking S4 @args slot...")
            debugLog("    hasSlot 'args': ", methods::hasMethod("@", "MetaboliteAssayData"))
            s4_args <- sessionData$current_s4_object@args
            debugLog("  STEP 9b: s4_args is NULL: ", is.null(s4_args))
            if (!is.null(s4_args)) {
                debugLog("    s4_args class: ", class(s4_args)[1])
                debugLog("    names(s4_args): ", paste(names(s4_args), collapse = ", "))
            }

            if (!is.null(s4_args) && !is.null(s4_args$daAnalysisParameters)) {
                formulaValue <- s4_args$daAnalysisParameters$formula_string
                if (!is.null(formulaValue) && nzchar(formulaValue)) {
                    daData$formula_from_s4 <- formulaValue
                    updateTextAreaInput(session, "formula_string", value = formulaValue)
                }
            }
        },
        error = function(e) {
            debugLog("  STEP 9 ERROR: ", e$message)
            logWarn(sprintf("Could not extract formula from S4: %s", e$message))
        }
    )

    invisible(list(
        stateName = stateName,
        contrastChoices = contrastChoices,
        assayChoices = assayChoices,
        formulaValue = formulaValue
    ))
}

