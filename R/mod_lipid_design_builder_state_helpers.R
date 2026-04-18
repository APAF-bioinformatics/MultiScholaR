buildLipidDesignResetConfirmationBody <- function(
    resetScope,
    htmlFn = shiny::HTML
) {
    htmlFn(paste0(
        "<p>This will revert <strong>", resetScope, "</strong> to their initial state.</p>"
    ))
}

buildLipidDesignResetConfirmationModal <- function(
    resetScope,
    nsFn,
    modalDialogFn = shiny::modalDialog,
    tagListFn = shiny::tagList,
    modalButtonFn = shiny::modalButton,
    actionButtonFn = shiny::actionButton,
    buildResetConfirmationBodyFn = buildLipidDesignResetConfirmationBody
) {
    modalDialogFn(
        title = "Confirm Reset"
        , buildResetConfirmationBodyFn(resetScope)
        , footer = tagListFn(
            modalButtonFn("Cancel")
            , actionButtonFn(nsFn("confirm_reset"), "Reset", class = "btn-danger")
        )
        , easyClose = TRUE
    )
}

registerLipidDesignResetRequestShells <- function(
    input,
    session,
    observeEventFn = shiny::observeEvent,
    showModalFn = shiny::showModal,
    buildResetConfirmationModalFn = buildLipidDesignResetConfirmationModal
) {
    observeEventFn(input$reset_changes, {
        showModalFn(buildResetConfirmationModalFn(
            resetScope = input$reset_scope
            , nsFn = session$ns
        ))
    })

    invisible(NULL)
}

applyLipidDesignResetState <- function(
    scope,
    initialState,
    designMatrix,
    dataClnReactive,
    sampleNamesReactive,
    removedSamples,
    factors,
    groups,
    contrasts,
    updateFormulaFn = function(value) invisible(value)
) {
    if (scope == "all" || scope == "sample_names") {
        designMatrix(initialState$design_matrix)
        dataClnReactive(initialState$data_cln)
        sampleNamesReactive(initialState$sample_names)
    }
    if (scope == "all" || scope == "removed_samples") {
        removedSamples(character(0))
    }
    if (scope == "all" || scope == "factors") {
        factors(initialState$factors)
        currentMatrix <- designMatrix()
        currentMatrix$factor1 <- NA_character_
        currentMatrix$factor2 <- NA_character_
        currentMatrix$factor3 <- NA_character_
        currentMatrix$group <- NA_character_
        currentMatrix$tech_reps <- NA_integer_
        designMatrix(currentMatrix)
        groups(initialState$groups)
    }
    if (scope == "all" || scope == "contrasts") {
        contrasts(initialState$contrasts)
    }
    if (scope == "all" || scope == "formula") {
        updateFormulaFn(initialState$formula)
    }

    invisible(NULL)
}

runLipidDesignResetConfirmationShell <- function(
    scope,
    initialState,
    designMatrix,
    dataClnReactive,
    sampleNamesReactive,
    removedSamples,
    factors,
    groups,
    contrasts,
    session,
    applyResetStateFn = applyLipidDesignResetState,
    removeModalFn = shiny::removeModal,
    showNotificationFn = shiny::showNotification,
    updateTextInputFn = shiny::updateTextInput
) {
    applyResetStateFn(
        scope = scope,
        initialState = initialState,
        designMatrix = designMatrix,
        dataClnReactive = dataClnReactive,
        sampleNamesReactive = sampleNamesReactive,
        removedSamples = removedSamples,
        factors = factors,
        groups = groups,
        contrasts = contrasts,
        updateFormulaFn = function(value) {
            updateTextInputFn(session, "formula_string", value = value)
        }
    )

    removeModalFn()
    showNotificationFn(paste("Reset of", scope, "completed."), type = "message")

    invisible(NULL)
}

registerLipidDesignResetConfirmationShells <- function(
    input,
    initialState,
    designMatrix,
    dataClnReactive,
    sampleNamesReactive,
    removedSamples,
    factors,
    groups,
    contrasts,
    session,
    observeEventFn = shiny::observeEvent,
    runResetConfirmationShellFn = runLipidDesignResetConfirmationShell
) {
    observeEventFn(input$confirm_reset, {
        runResetConfirmationShellFn(
            scope = input$reset_scope,
            initialState = initialState(),
            designMatrix = designMatrix,
            dataClnReactive = dataClnReactive,
            sampleNamesReactive = sampleNamesReactive,
            removedSamples = removedSamples,
            factors = factors,
            groups = groups,
            contrasts = contrasts,
            session = session
        )
    })

    invisible(NULL)
}

getLipidDesignSampleColumns <- function(
    assayList,
    colMap,
    isNumericFn = is.numeric,
    greplFn = grepl,
    messageFn = message
) {
    messageFn("DEBUG66: get_sample_columns() called")
    messageFn(sprintf(
        "DEBUG66: assay_list is NULL: %s, length: %d",
        is.null(assayList),
        if (is.null(assayList)) 0 else length(assayList)
    ))
    messageFn(sprintf("DEBUG66: col_map is NULL: %s", is.null(colMap)))

    if (is.null(assayList) || length(assayList) == 0) {
        messageFn("DEBUG66: Returning empty character(0) - assay_list is NULL or empty")
        return(character(0))
    }

    firstAssay <- assayList[[1]]
    if (!is.data.frame(firstAssay)) {
        messageFn(sprintf(
            "DEBUG66: WARNING - first_assay is not a data frame, class: %s",
            class(firstAssay)[1]
        ))
        return(character(0))
    }

    allCols <- names(firstAssay)
    messageFn(sprintf("DEBUG66: first_assay has %d columns", length(allCols)))

    if (is.null(colMap)) {
        messageFn("DEBUG66: col_map is NULL, using all columns as sample candidates")
        numericCols <- names(firstAssay)[vapply(firstAssay, isNumericFn, logical(1))]
        messageFn(sprintf("DEBUG66: Found %d numeric columns", length(numericCols)))
        return(numericCols)
    }

    excludeCols <- c(
        colMap$lipid_id_col,
        colMap$annotation_col
    )
    excludeCols <- excludeCols[!is.na(excludeCols) & nzchar(excludeCols)]
    messageFn(sprintf(
        "DEBUG66: Excluding columns: %s",
        paste(excludeCols, collapse = ", ")
    ))

    sampleCols <- setdiff(allCols, excludeCols)
    nonSamplePatterns <- c(
        "^Alignment.ID$",
        "^Average.Rt",
        "^Average.Mz",
        "^Lipid.name$",
        "^Adduct.type$",
        "^Post.curation.result$",
        "^Fill.%$",
        "^MS/MS.spectrum$",
        "^Reference.RT$",
        "^Reference.m/z$",
        "^Formula$",
        "^Ontology$",
        "^INCHIKEY$",
        "^SMILES$",
        "^Annotation.tag",
        "^RT.matched$",
        "^m/z.matched$",
        "^MS/MS.matched$",
        "^Comment$",
        "^Manually.modified$",
        "^Total.score$",
        "^RT.similarity$",
        "^Dot.product$",
        "^Reverse.dot.product$",
        "^Fragment.presence",
        "^S/N.average$",
        "^Spectrum.reference",
        "^MS1.isotopic.spectrum$",
        "^MS/MS.spectrum$"
    )

    for (pattern in nonSamplePatterns) {
        sampleCols <- sampleCols[!greplFn(pattern, sampleCols, ignore.case = TRUE)]
    }

    if (!is.null(colMap$sample_columns) && length(colMap$sample_columns) > 0) {
        sampleCols <- colMap$sample_columns
        messageFn(sprintf(
            "DEBUG66: Using sample_columns from col_map: %d columns",
            length(sampleCols)
        ))
    }

    messageFn(sprintf("DEBUG66: Returning %d sample columns", length(sampleCols)))
    sampleCols
}

buildLipidDesignContrastState <- function(existingContrasts) {
    if (is.null(existingContrasts) || !is.data.frame(existingContrasts) || nrow(existingContrasts) == 0) {
        return(data.frame(
            contrast_name = character(),
            numerator = character(),
            denominator = character(),
            stringsAsFactors = FALSE
        ))
    }

    parsedContrasts <- lapply(existingContrasts$contrasts, function(contrastString) {
        parts <- strsplit(gsub("^group", "", contrastString), "-group")[[1]]
        if (length(parts) == 2) {
            return(list(
                contrast_name = paste0(parts[1], ".vs.", parts[2]),
                numerator = parts[1],
                denominator = parts[2]
            ))
        }

        list(
            contrast_name = contrastString,
            numerator = contrastString,
            denominator = ""
        )
    })

    data.frame(
        contrast_name = vapply(parsedContrasts, `[[`, character(1), "contrast_name"),
        numerator = vapply(parsedContrasts, `[[`, character(1), "numerator"),
        denominator = vapply(parsedContrasts, `[[`, character(1), "denominator"),
        stringsAsFactors = FALSE
    )
}

buildLipidDesignInitialState <- function(
    assayList,
    configList,
    colMap,
    existingDesignMatrix = NULL,
    existingContrasts = NULL,
    getSampleColumnsFn = getLipidDesignSampleColumns,
    buildContrastStateFn = buildLipidDesignContrastState,
    mixedsortFn = gtools::mixedsort,
    logErrorFn = logger::log_error,
    messageFn = message
) {
    sampleNames <- getSampleColumnsFn(assayList = assayList, colMap = colMap)

    if (length(sampleNames) == 0) {
        logErrorFn("No sample columns detected in lipidomics data")
        return(NULL)
    }

    formulaString <- configList[["deAnalysisParameters"]][["formula_string"]]
    resolvedExistingDesign <- NULL
    resolvedExistingContrasts <- NULL
    useExisting <- FALSE

    if (!is.null(existingDesignMatrix)) {
        resolvedExistingDesign <- tryCatch(
            if (is.function(existingDesignMatrix)) existingDesignMatrix() else existingDesignMatrix,
            error = function(error) NULL
        )

        if (!is.null(resolvedExistingDesign) &&
            is.data.frame(resolvedExistingDesign) &&
            nrow(resolvedExistingDesign) > 0 &&
            "group" %in% names(resolvedExistingDesign)) {
            nonNaGroups <- resolvedExistingDesign$group[
                !is.na(resolvedExistingDesign$group) &
                    resolvedExistingDesign$group != ""
            ]

            if (length(nonNaGroups) > 0) {
                useExisting <- TRUE
                messageFn(sprintf(
                    "DEBUG66: Using existing design matrix with %d assigned samples",
                    length(nonNaGroups)
                ))
            }
        }
    }

    if (!is.null(existingContrasts)) {
        resolvedExistingContrasts <- tryCatch(
            if (is.function(existingContrasts)) existingContrasts() else existingContrasts,
            error = function(error) NULL
        )
    }

    if (useExisting && !is.null(resolvedExistingDesign)) {
        messageFn("DEBUG66: initial_state using IMPORTED design matrix")

        groupsFromDesign <- unique(
            resolvedExistingDesign$group[
                !is.na(resolvedExistingDesign$group) &
                    resolvedExistingDesign$group != ""
            ]
        )
        factorsFromDesign <- character(0)

        if ("factor1" %in% names(resolvedExistingDesign)) {
            factorsFromDesign <- c(
                factorsFromDesign,
                unique(
                    resolvedExistingDesign$factor1[
                        !is.na(resolvedExistingDesign$factor1) &
                            resolvedExistingDesign$factor1 != ""
                    ]
                )
            )
        }

        if ("factor2" %in% names(resolvedExistingDesign)) {
            factorsFromDesign <- c(
                factorsFromDesign,
                unique(
                    resolvedExistingDesign$factor2[
                        !is.na(resolvedExistingDesign$factor2) &
                            resolvedExistingDesign$factor2 != ""
                    ]
                )
            )
        }

        return(list(
            design_matrix = resolvedExistingDesign,
            data_cln = assayList,
            sample_names = sampleNames,
            groups = groupsFromDesign,
            factors = unique(factorsFromDesign),
            formula = formulaString,
            contrasts = buildContrastStateFn(resolvedExistingContrasts)
        ))
    }

    messageFn("DEBUG66: initial_state creating FRESH design matrix")

    designMatrixRaw <- tibble::tibble(
        Run = mixedsortFn(sampleNames),
        group = NA_character_,
        factor1 = NA_character_,
        factor2 = NA_character_,
        factor3 = NA_character_,
        batch = NA_character_,
        replicates = NA_integer_,
        tech_reps = NA_integer_
    )

    list(
        design_matrix = designMatrixRaw,
        data_cln = assayList,
        sample_names = sampleNames,
        groups = unique(
            designMatrixRaw$group[!is.na(designMatrixRaw$group) & designMatrixRaw$group != ""]
        ),
        factors = character(0),
        formula = formulaString,
        contrasts = buildContrastStateFn(NULL)
    )
}

runLipidDesignInitialStateShell <- function(
    state,
    session,
    designMatrix,
    dataClnReactive,
    sampleNamesReactive,
    groups,
    factors,
    contrasts,
    removedSamples,
    updateSelectizeInputFn = shiny::updateSelectizeInput,
    updateTextInputFn = shiny::updateTextInput
) {
    designMatrix(state$design_matrix)
    dataClnReactive(state$data_cln)
    sampleNamesReactive(state$sample_names)
    groups(state$groups)
    factors(state$factors)
    contrasts(state$contrasts)
    removedSamples(character(0))

    sortedRuns <- state$design_matrix$Run
    updateSelectizeInputFn(session, "sample_to_rename", choices = sortedRuns, selected = "")
    updateSelectizeInputFn(session, "selected_runs", choices = sortedRuns, selected = "")
    updateSelectizeInputFn(session, "samples_to_transform", choices = sortedRuns, selected = "")
    updateSelectizeInputFn(session, "tech_rep_samples", choices = sortedRuns, selected = "")
    updateSelectizeInputFn(session, "samples_to_remove", choices = sortedRuns, selected = "")
    updateTextInputFn(session, "formula_string", value = state$formula)

    invisible(state)
}

registerLipidDesignInitialStateShells <- function(
    dataTbl,
    input,
    initialState,
    session,
    designMatrix,
    dataClnReactive,
    sampleNamesReactive,
    groups,
    factors,
    contrasts,
    removedSamples,
    observeFn = shiny::observe,
    bindEventFn = shiny::bindEvent,
    reqFn = shiny::req,
    runInitialStateShellFn = runLipidDesignInitialStateShell
) {
    initializationObserver <- observeFn({
        state <- initialState()
        reqFn(state)

        runInitialStateShellFn(
            state = state,
            session = session,
            designMatrix = designMatrix,
            dataClnReactive = dataClnReactive,
            sampleNamesReactive = sampleNamesReactive,
            groups = groups,
            factors = factors,
            contrasts = contrasts,
            removedSamples = removedSamples
        )
    })

    bindEventFn(initializationObserver, dataTbl(), input$reset_changes)

    invisible(NULL)
}

buildLipidDesignSaveResultsDesignMatrix <- function(
    designMatrix,
    currentRemovedSamples = character(0)
) {
    designMatrix |>
        dplyr::filter(!is.na(group) & group != "") |>
        dplyr::filter(!Run %in% currentRemovedSamples) |>
        dplyr::mutate(tech_rep_group = paste(group, replicates, sep = "_"))
}

buildLipidDesignSaveResultsDataList <- function(
    dataCln,
    assignedSamples,
    colMap
) {
    lapply(dataCln, function(assayDf) {
        keepCols <- c(
            colMap$lipid_id_col
            , colMap$annotation_col
        )
        keepCols <- keepCols[!is.na(keepCols) & nzchar(keepCols)]
        keepCols <- keepCols[keepCols %in% names(assayDf)]

        sampleColsInAssay <- intersect(assignedSamples, names(assayDf))
        finalCols <- unique(c(keepCols, sampleColsInAssay))

        assayDf[, finalCols, drop = FALSE]
    })
}

buildLipidDesignSaveResultsContrastsTable <- function(
    contrastData,
    formulaString
) {
    if (is.null(contrastData) || nrow(contrastData) == 0) {
        return(NULL)
    }

    formulaUsesGroupPrefix <- grepl("~ *0 *\\+ *group", formulaString)

    contrastInfo <- lapply(seq_len(nrow(contrastData)), function(index) {
        friendlyName <- gsub("\\.", "_", contrastData$contrast_name[index])

        contrastString <- if (formulaUsesGroupPrefix) {
            paste0(
                "group", contrastData$numerator[index]
                , "-group", contrastData$denominator[index]
            )
        } else {
            paste0(
                contrastData$numerator[index]
                , "-", contrastData$denominator[index]
            )
        }

        list(
            friendly_name = friendlyName
            , contrast_string = contrastString
            , full_format = paste0(friendlyName, "=", contrastString)
        )
    })

    data.frame(
        contrasts = vapply(contrastInfo, `[[`, character(1), "contrast_string")
        , friendly_names = vapply(contrastInfo, `[[`, character(1), "friendly_name")
        , full_format = vapply(contrastInfo, `[[`, character(1), "full_format")
        , stringsAsFactors = FALSE
    )
}

buildLipidDesignSaveResultsPayload <- function(
    designMatrix,
    currentRemovedSamples = character(0),
    dataCln,
    colMap,
    contrastData,
    configList,
    formulaString,
    buildSaveResultsDesignMatrixFn = buildLipidDesignSaveResultsDesignMatrix,
    buildSaveResultsDataListFn = buildLipidDesignSaveResultsDataList,
    buildSaveResultsContrastsTableFn = buildLipidDesignSaveResultsContrastsTable
) {
    designMatrixFinal <- buildSaveResultsDesignMatrixFn(
        designMatrix = designMatrix,
        currentRemovedSamples = currentRemovedSamples
    )

    if (nrow(designMatrixFinal) == 0) {
        return(NULL)
    }

    configListFinal <- configList
    configListFinal[["deAnalysisParameters"]][["formula_string"]] <- formulaString

    list(
        design_matrix = designMatrixFinal
        , data_cln = buildSaveResultsDataListFn(
            dataCln = dataCln,
            assignedSamples = designMatrixFinal$Run,
            colMap = colMap
        )
        , contrasts_tbl = buildSaveResultsContrastsTableFn(
            contrastData = contrastData,
            formulaString = formulaString
        )
        , config_list = configListFinal
    )
}

runLipidDesignSaveResultsShell <- function(
    designMatrix,
    currentRemovedSamples = character(0),
    dataCln,
    colMap,
    contrastData,
    configList,
    formulaString,
    resultSetter,
    buildSaveResultsPayloadFn = buildLipidDesignSaveResultsPayload,
    showNotificationFn = shiny::showNotification
) {
    finalResult <- buildSaveResultsPayloadFn(
        designMatrix = designMatrix,
        currentRemovedSamples = currentRemovedSamples,
        dataCln = dataCln,
        colMap = colMap,
        contrastData = contrastData,
        configList = configList,
        formulaString = formulaString
    )

    if (is.null(finalResult)) {
        showNotificationFn("No samples have been assigned to groups.", type = "warning")
        return(NULL)
    }

    resultSetter(finalResult)
    showNotificationFn("Design saved successfully!", type = "message", duration = 5)

    finalResult
}

registerLipidDesignSaveResultsShells <- function(
    input,
    designMatrix,
    removedSamples,
    dataClnReactive,
    columnMapping,
    contrasts,
    configList,
    resultRv,
    observeEventFn = shiny::observeEvent,
    runSaveResultsShellFn = runLipidDesignSaveResultsShell
) {
    observeEventFn(input$save_results, {
        runSaveResultsShellFn(
            designMatrix = designMatrix(),
            currentRemovedSamples = removedSamples(),
            dataCln = dataClnReactive(),
            colMap = columnMapping(),
            contrastData = contrasts(),
            configList = configList(),
            formulaString = input$formula_string,
            resultSetter = resultRv
        )
    })

    invisible(NULL)
}

