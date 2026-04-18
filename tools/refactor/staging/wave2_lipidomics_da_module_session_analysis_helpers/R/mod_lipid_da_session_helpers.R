restoreLipidDaContrastsFromSession <- function(
    sessionData,
    workflowData,
    daData,
    session,
    updateSelectInputFn = shiny::updateSelectInput,
    logInfo = logger::log_info
) {
    contrastsTbl <- sessionData$contrasts_tbl

    if (is.null(contrastsTbl) || nrow(contrastsTbl) == 0) {
        return(NULL)
    }

    workflowData$contrasts_tbl <- contrastsTbl
    daData$contrasts_available <- contrastsTbl

    contrastChoices <- if ("friendly_names" %in% colnames(contrastsTbl)) {
        contrastsTbl$friendly_names
    } else {
        contrastsTbl$contrasts
    }

    updateSelectInputFn(session, "volcano_contrast",
        choices = contrastChoices, selected = contrastChoices[1]
    )
    updateSelectInputFn(session, "heatmap_contrast",
        choices = contrastChoices, selected = contrastChoices[1]
    )
    updateSelectInputFn(session, "table_contrast",
        choices = contrastChoices, selected = contrastChoices[1]
    )

    logInfo(sprintf("   Restored %d contrasts", nrow(contrastsTbl)))

    list(
        contrastsTbl = contrastsTbl,
        contrastChoices = contrastChoices
    )
}

restoreLipidDaAssaysFromSession <- function(
    sessionData,
    daData,
    session,
    updateSelectInputFn = shiny::updateSelectInput,
    debugLog = NULL
) {
    assayNames <- sessionData$assay_names

    if (is.null(assayNames)) {
        return(NULL)
    }

    daData$assays_available <- assayNames

    assayChoices <- c("Combined", assayNames)
    tableAssayChoices <- c("All", assayNames)

    if (!is.null(debugLog)) {
        debugLog("  STEP 8a: Updating volcano_assay dropdown...")
    }
    updateSelectInputFn(session, "volcano_assay", choices = assayChoices)

    if (!is.null(debugLog)) {
        debugLog("  STEP 8b: Updating heatmap_assay dropdown...")
    }
    updateSelectInputFn(session, "heatmap_assay", choices = assayChoices)

    if (!is.null(debugLog)) {
        debugLog("  STEP 8c: Updating table_assay dropdown...")
    }
    updateSelectInputFn(session, "table_assay", choices = tableAssayChoices)

    list(
        assayNames = assayNames,
        assayChoices = assayChoices,
        tableAssayChoices = tableAssayChoices
    )
}

restoreLipidDaFormulaFromSession <- function(
    sessionData,
    daData,
    session,
    updateTextAreaInputFn = shiny::updateTextAreaInput,
    debugLog = NULL,
    hasArgsMethodFn = function() methods::hasMethod("@", "LipidomicsAssayData"),
    logWarn = logger::log_warn
) {
    formulaVal <- NULL

    tryCatch(
        {
            if (!is.null(debugLog)) {
                debugLog("  STEP 9a: Checking S4 @args slot...")
                debugLog("    hasSlot 'args': ", hasArgsMethodFn())
            }

            s4Args <- sessionData$current_s4_object@args

            if (!is.null(debugLog)) {
                debugLog("  STEP 9b: s4_args is NULL: ", is.null(s4Args))
                if (!is.null(s4Args)) {
                    debugLog("    s4_args class: ", class(s4Args)[1])
                    debugLog("    names(s4_args): ", paste(names(s4Args), collapse = ", "))
                }
            }

            if (!is.null(s4Args) && !is.null(s4Args$daAnalysisParameters)) {
                formulaVal <- s4Args$daAnalysisParameters$formula_string
                if (!is.null(formulaVal) && nzchar(formulaVal)) {
                    daData$formula_from_s4 <- formulaVal
                    updateTextAreaInputFn(session, "formula_string", value = formulaVal)
                }
            }
        },
        error = function(e) {
            if (!is.null(debugLog)) {
                debugLog("  STEP 9 ERROR: ", e$message)
            }
            logWarn(sprintf("Could not extract formula from S4: %s", e$message))
        }
    )

    if (is.null(formulaVal) || !nzchar(formulaVal)) {
        return(NULL)
    }

    list(formula = formulaVal)
}

notifyLipidDaSessionSourceDirError <- function(
    notify = shiny::showNotification,
    debugLog = NULL
) {
    contract <- list(
        loggedMessage = "  ERROR: No valid source directory found",
        notificationMessage = "Could not find source directory for session data.",
        notificationType = "error",
        notificationDuration = 5
    )

    if (!is.null(debugLog)) {
        debugLog(contract$loggedMessage)
    }

    notify(
        contract$notificationMessage,
        type = contract$notificationType,
        duration = contract$notificationDuration
    )

    contract
}

notifyLipidDaSessionFileMissing <- function(
    sessionFile,
    notify = shiny::showNotification
) {
    contract <- list(
        sessionFile = sessionFile,
        notificationMessage = sprintf("Session file not found: %s", sessionFile),
        notificationType = "error",
        notificationDuration = 5
    )

    notify(
        contract$notificationMessage,
        type = contract$notificationType,
        duration = contract$notificationDuration
    )

    contract
}

resolveLipidDaSessionFile <- function(
    experimentPaths,
    debugLog = NULL,
    dirExistsFn = dir.exists,
    fileExistsFn = file.exists,
    filePathFn = file.path,
    notifySourceDirFn = notifyLipidDaSessionSourceDirError,
    notifySessionFileMissingFn = notifyLipidDaSessionFileMissing
) {
    sourceDir <- experimentPaths$source_dir

    if (!is.null(debugLog)) {
        debugLog("  STEP 2: source_dir = ", if (is.null(sourceDir)) "NULL" else sourceDir)
    }

    if (is.null(sourceDir) || !dirExistsFn(sourceDir)) {
        if (!is.null(debugLog)) {
            debugLog("  BRANCH: source_dir NULL or not exists, trying export_dir")
        }

        sourceDir <- experimentPaths$export_dir

        if (!is.null(debugLog)) {
            debugLog("  STEP 3: export_dir = ", if (is.null(sourceDir)) "NULL" else sourceDir)
        }
    }

    if (is.null(sourceDir) || !dirExistsFn(sourceDir)) {
        notifySourceDirFn(debugLog = debugLog)
        return(NULL)
    }

    sessionFile <- filePathFn(sourceDir, "lipid_filtered_session_data_latest.rds")
    sessionFileExists <- fileExistsFn(sessionFile)

    if (!is.null(debugLog)) {
        debugLog("  STEP 4: session_file = ", sessionFile)
        debugLog("  STEP 4: file.exists = ", sessionFileExists)
    }

    if (!sessionFileExists) {
        notifySessionFileMissingFn(sessionFile = sessionFile)
        return(NULL)
    }

    list(
        sourceDir = sourceDir,
        sessionFile = sessionFile
    )
}

showLipidDaSessionLoadingNotification <- function(
    notify = shiny::showNotification
) {
    contract <- list(
        notificationMessage = "Loading filtered session...",
        notificationId = "loading_session",
        notificationDuration = NULL
    )

    notify(
        contract$notificationMessage,
        id = contract$notificationId,
        duration = contract$notificationDuration
    )

    contract
}

readLipidDaSessionData <- function(
    sessionFile,
    readRdsFn = readRDS,
    debugLog = NULL,
    logInfo = logger::log_info
) {
    if (!is.null(debugLog)) {
        debugLog("  STEP 5: Reading RDS file...")
    }

    sessionData <- readRdsFn(sessionFile)
    sessionDataNames <- names(sessionData)

    if (!is.null(debugLog)) {
        debugLog("  STEP 5: RDS loaded successfully")
        debugLog("    names(session_data): ", paste(sessionDataNames, collapse = ", "))
    }

    logInfo(sprintf("   Loaded session from: %s", sessionFile))

    list(
        sessionFile = sessionFile,
        sessionData = sessionData,
        sessionDataNames = sessionDataNames
    )
}

restoreLipidDaCurrentS4FromSession <- function(
    sessionData,
    daData,
    workflowData,
    sessionFile,
    debugLog = NULL,
    logInfo = logger::log_info
) {
    currentS4Object <- sessionData$current_s4_object
    if (is.null(currentS4Object)) {
        return(NULL)
    }

    if (!is.null(debugLog)) {
        debugLog("  STEP 6a: S4 object class: ", class(currentS4Object)[1])
        debugLog("  STEP 6b: Assigning to da_data$current_s4_object...")
    }
    daData$current_s4_object <- currentS4Object

    if (!is.null(debugLog)) {
        debugLog("  STEP 6c: Getting state_name...")
    }
    stateName <- sessionData$r6_current_state_name
    if (!is.null(debugLog)) {
        debugLog("    r6_current_state_name: ", if (is.null(stateName)) "NULL" else stateName)
    }
    if (is.null(stateName)) {
        stateName <- "loaded_for_de"
    }

    if (!is.null(debugLog)) {
        debugLog("  STEP 6d: Checking workflow_data$state_manager...")
        debugLog("    workflow_data is NULL: ", is.null(workflowData))
        debugLog("    workflow_data class: ", class(workflowData)[1])
        debugLog("    state_manager exists: ", !is.null(workflowData$state_manager))
        if (!is.null(workflowData$state_manager)) {
            debugLog("    state_manager class: ", class(workflowData$state_manager)[1])
        }
        debugLog("  STEP 6e: Calling workflow_data$state_manager$saveState...")
    }

    workflowData$state_manager$saveState(
        state_name = stateName,
        s4_data_object = currentS4Object,
        config_object = list(loaded_from = sessionFile),
        description = "Loaded from filtered session for DA analysis"
    )

    if (!is.null(debugLog)) {
        debugLog("  STEP 6f: saveState completed")
    }
    logInfo("   Restored S4 object to state manager")

    list(
        currentS4Object = currentS4Object,
        stateName = stateName,
        sessionFile = sessionFile,
        description = "Loaded from filtered session for DA analysis"
    )
}

restoreLipidDaPostReadSessionState <- function(
    sessionData,
    workflowData,
    daData,
    session,
    debugLog = NULL,
    restoreContrastsFn = restoreLipidDaContrastsFromSession,
    restoreAssaysFn = restoreLipidDaAssaysFromSession,
    restoreFormulaFn = restoreLipidDaFormulaFromSession
) {
    contrastRestore <- NULL
    assayRestore <- NULL
    formulaRestore <- NULL

    if (!is.null(debugLog)) {
        debugLog("  STEP 7: Checking contrasts_tbl...")
        debugLog("    is.null(session_data$contrasts_tbl): ", is.null(sessionData$contrasts_tbl))
        if (!is.null(sessionData$contrasts_tbl)) {
            debugLog("    nrow(contrasts_tbl): ", nrow(sessionData$contrasts_tbl))
            debugLog("    colnames: ", paste(colnames(sessionData$contrasts_tbl), collapse = ", "))
        }
    }

    if (!is.null(sessionData$contrasts_tbl) && nrow(sessionData$contrasts_tbl) > 0) {
        if (!is.null(debugLog)) {
            debugLog("  STEP 7a: Restoring contrast state and dropdowns...")
        }
        contrastRestore <- restoreContrastsFn(
            sessionData = sessionData,
            workflowData = workflowData,
            daData = daData,
            session = session
        )
        if (!is.null(debugLog)) {
            debugLog(
                "    contrast_choices: ",
                paste(contrastRestore$contrastChoices, collapse = ", ")
            )
        }
    }

    if (!is.null(debugLog)) {
        debugLog("  STEP 8: Checking assay_names...")
        debugLog("    is.null(session_data$assay_names): ", is.null(sessionData$assay_names))
        if (!is.null(sessionData$assay_names)) {
            debugLog("    assay_names: ", paste(sessionData$assay_names, collapse = ", "))
        }
    }

    if (!is.null(sessionData$assay_names)) {
        assayRestore <- restoreAssaysFn(
            sessionData = sessionData,
            daData = daData,
            session = session,
            debugLog = debugLog
        )
        if (!is.null(debugLog)) {
            debugLog(
                "    assay_choices: ",
                paste(assayRestore$assayChoices, collapse = ", ")
            )
            debugLog(
                "    table_assay_choices: ",
                paste(assayRestore$tableAssayChoices, collapse = ", ")
            )
        }
    }

    if (!is.null(debugLog)) {
        debugLog("  STEP 9: Attempting to extract formula from S4 args...")
    }
    formulaRestore <- restoreFormulaFn(
        sessionData = sessionData,
        daData = daData,
        session = session,
        debugLog = debugLog
    )

    list(
        contrastRestore = contrastRestore,
        assayRestore = assayRestore,
        formulaRestore = formulaRestore
    )
}

loadLipidDaSessionFromFile <- function(
    sessionFile,
    workflowData,
    daData,
    session,
    debugLog = NULL,
    startTime = NULL,
    showLoadingFn = showLipidDaSessionLoadingNotification,
    readSessionFn = readLipidDaSessionData,
    restoreCurrentS4Fn = restoreLipidDaCurrentS4FromSession,
    restorePostReadFn = restoreLipidDaPostReadSessionState,
    finalizeSuccessFn = finalizeLipidDaSessionLoadSuccess,
    finalizeErrorFn = finalizeLipidDaSessionLoadError,
    nowFn = Sys.time,
    deparseCallFn = deparse
) {
    tryCatch(
        {
            showLoadingContract <- showLoadingFn()

            sessionLoad <- readSessionFn(
                sessionFile = sessionFile,
                debugLog = debugLog
            )
            sessionData <- sessionLoad$sessionData

            if (!is.null(debugLog)) {
                debugLog("  STEP 6: Checking current_s4_object...")
                debugLog("    is.null(session_data$current_s4_object): ", is.null(sessionData$current_s4_object))
            }

            currentS4Restore <- NULL
            if (!is.null(sessionData$current_s4_object)) {
                currentS4Restore <- restoreCurrentS4Fn(
                    sessionData = sessionData,
                    daData = daData,
                    workflowData = workflowData,
                    sessionFile = sessionFile,
                    debugLog = debugLog
                )
            }

            postReadRestore <- restorePostReadFn(
                sessionData = sessionData,
                workflowData = workflowData,
                daData = daData,
                session = session,
                debugLog = debugLog
            )

            successContract <- finalizeSuccessFn(debugLog = debugLog)
            elapsedSeconds <- NULL

            if (!is.null(startTime)) {
                elapsedSeconds <- round(difftime(nowFn(), startTime, units = "secs"), 3)
            }

            if (!is.null(debugLog)) {
                if (is.null(elapsedSeconds)) {
                    debugLog("--- EXIT load_filtered_session --- SUCCESS")
                } else {
                    debugLog(
                        "--- EXIT load_filtered_session (",
                        elapsedSeconds, "s) --- SUCCESS"
                    )
                }
            }

            list(
                status = "success",
                showLoadingContract = showLoadingContract,
                sessionLoad = sessionLoad,
                currentS4Restore = currentS4Restore,
                postReadRestore = postReadRestore,
                successContract = successContract,
                elapsedSeconds = elapsedSeconds
            )
        },
        error = function(e) {
            errorCall <- deparseCallFn(e$call)

            if (!is.null(debugLog)) {
                debugLog("  FATAL ERROR: ", e$message)
                debugLog("  ERROR CALL: ", errorCall)
                debugLog("--- EXIT load_filtered_session --- FAILED")
            }

            errorContract <- finalizeErrorFn(errorMessage = e$message)

            list(
                status = "error",
                errorMessage = e$message,
                errorCall = errorCall,
                errorContract = errorContract
            )
        }
    )
}

finalizeLipidDaSessionLoadSuccess <- function(
    removeNotificationFn = shiny::removeNotification,
    notify = shiny::showNotification,
    debugLog = NULL
) {
    contract <- list(
        loadingNotificationId = "loading_session",
        successMessage = "Session loaded successfully!",
        successType = "message",
        successDuration = 3
    )

    if (!is.null(debugLog)) {
        debugLog("  STEP 10: Removing notification and showing success...")
    }

    removeNotificationFn(contract$loadingNotificationId)
    notify(
        contract$successMessage,
        type = contract$successType,
        duration = contract$successDuration
    )

    contract
}

finalizeLipidDaSessionLoadError <- function(
    errorMessage,
    removeNotificationFn = shiny::removeNotification,
    notify = shiny::showNotification,
    logError = logger::log_error
) {
    contract <- list(
        loadingNotificationId = "loading_session",
        loggedMessage = sprintf("   Error loading session: %s", errorMessage),
        notificationMessage = sprintf("Error loading session: %s", errorMessage),
        notificationType = "error",
        notificationDuration = 10
    )

    logError(contract$loggedMessage)
    removeNotificationFn(contract$loadingNotificationId)
    notify(
        contract$notificationMessage,
        type = contract$notificationType,
        duration = contract$notificationDuration
    )

    contract
}

bootstrapLipidDaLoadFilteredSession <- function(
    experimentPaths,
    workflowData,
    daData,
    session,
    messageFn = message,
    nowFn = Sys.time,
    handleLoadFn = handleLipidDaLoadFilteredSession
) {
    d66Log <- function(...) {
        messageFn(sprintf("[D66] %s", paste0(...)))
    }

    d66Log("--- ENTER load_filtered_session observer ---")
    d66Start <- nowFn()

    handleLoadFn(
        experimentPaths = experimentPaths,
        workflowData = workflowData,
        daData = daData,
        session = session,
        debugLog = d66Log,
        startTime = d66Start
    )
}

handleLipidDaLoadFilteredSession <- function(
    experimentPaths,
    workflowData,
    daData,
    session,
    debugLog = NULL,
    startTime = NULL,
    logInfo = logger::log_info,
    resolveSessionSourceFn = resolveLipidDaSessionFile,
    loadSessionFromFileFn = loadLipidDaSessionFromFile
) {
    logInfo("=== LOAD FILTERED SESSION BUTTON CLICKED ===")

    if (!is.null(debugLog)) {
        debugLog("  STEP 1: Getting source_dir from experiment_paths")
        debugLog("    experiment_paths class: ", class(experimentPaths)[1])
        debugLog("    experiment_paths is NULL: ", is.null(experimentPaths))
    }

    sessionSource <- resolveSessionSourceFn(
        experimentPaths = experimentPaths,
        debugLog = debugLog
    )

    if (is.null(sessionSource)) {
        return(NULL)
    }

    loadResult <- loadSessionFromFileFn(
        sessionFile = sessionSource$sessionFile,
        workflowData = workflowData,
        daData = daData,
        session = session,
        debugLog = debugLog,
        startTime = startTime
    )

    list(
        sessionSource = sessionSource,
        loadResult = loadResult
    )
}

