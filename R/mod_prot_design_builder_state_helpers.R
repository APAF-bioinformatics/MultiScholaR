buildProtDesignSaveResultsContrastsTable <- function(
    contrastData,
    formulaString
) {
    if (is.null(contrastData) || nrow(contrastData) == 0) {
        return(NULL)
    }

    formulaUsesGroupPrefix <- grepl("~ *0 *\\+ *group", formulaString)

    contrastInfo <- lapply(seq_len(nrow(contrastData)), function(i) {
        friendlyName <- gsub("\\.", "_", contrastData$contrast_name[i])
        friendlyName <- gsub("\\.vs\\.", "_minus_", friendlyName)

        contrastString <- if (formulaUsesGroupPrefix) {
            paste0(
                "group", contrastData$numerator[i],
                "-group", contrastData$denominator[i]
            )
        } else {
            paste0(
                contrastData$numerator[i],
                "-", contrastData$denominator[i]
            )
        }

        list(
            friendly_name = friendlyName,
            contrast_string = contrastString,
            full_format = paste0(friendlyName, "=", contrastString)
        )
    })

    data.frame(
        contrasts = vapply(contrastInfo, `[[`, character(1), "contrast_string"),
        friendly_names = vapply(contrastInfo, `[[`, character(1), "friendly_name"),
        full_format = vapply(contrastInfo, `[[`, character(1), "full_format"),
        stringsAsFactors = FALSE
    )
}

buildProtDesignSaveResultsPayload <- function(
    designMatrix,
    currentRemovedSamples = character(0),
    dataCln,
    contrastData,
    configList,
    formulaString
) {
    designMatrixFinal <- designMatrix |>
        dplyr::filter(!is.na(group) & group != "") |>
        dplyr::filter(!Run %in% currentRemovedSamples) |>
        # group + replicates uniquely identifies a biological sample for
        # downstream technical-replicate grouping.
        dplyr::mutate(tech_rep_group = paste(group, replicates, sep = "_"))

    if (nrow(designMatrixFinal) == 0) {
        return(NULL)
    }

    assignedRuns <- designMatrixFinal$Run
    dataClnFinal <- dataCln |>
        dplyr::filter(Run %in% assignedRuns)

    configListFinal <- configList
    configListFinal[["deAnalysisParameters"]][["formula_string"]] <- formulaString

    list(
        design_matrix = designMatrixFinal,
        data_cln = dataClnFinal,
        contrasts_tbl = buildProtDesignSaveResultsContrastsTable(
            contrastData = contrastData,
            formulaString = formulaString
        ),
        config_list = configListFinal
    )
}

runProtDesignSaveResultsObserverShell <- function(
    designMatrix,
    currentRemovedSamples = character(0),
    dataCln,
    contrastData,
    configList,
    formulaString,
    resultSetter,
    buildSaveResultsPayload = buildProtDesignSaveResultsPayload,
    showNotification = shiny::showNotification
) {
    finalResult <- buildSaveResultsPayload(
        designMatrix = designMatrix,
        currentRemovedSamples = currentRemovedSamples,
        dataCln = dataCln,
        contrastData = contrastData,
        configList = configList,
        formulaString = formulaString
    )

    if (is.null(finalResult)) {
        showNotification(
            "No samples have been assigned to groups. Please assign metadata before saving.",
            type = "warning"
        )
        return(NULL)
    }

    resultSetter(finalResult)

    showNotification(
        "Design saved successfully. You can close this builder.",
        type = "message",
        duration = 5
    )

    finalResult
}

registerProtDesignSaveResultsObserver <- function(
    input,
    designMatrix,
    removedSamples,
    dataCln,
    contrastData,
    configList,
    resultRv,
    observeEventFn = shiny::observeEvent,
    runSaveResultsObserverShell = runProtDesignSaveResultsObserverShell
) {
    observeEventFn(input$save_results, {
        runSaveResultsObserverShell(
            designMatrix = designMatrix(),
            currentRemovedSamples = removedSamples(),
            dataCln = dataCln(),
            contrastData = contrastData(),
            configList = configList(),
            formulaString = input$formula_string,
            resultSetter = resultRv
        )
    })
}

showProtDesignResetConfirmationModal <- function(
    resetScope,
    nsFn = identity,
    showModalFn = shiny::showModal,
    modalDialogFn = shiny::modalDialog,
    htmlFn = shiny::HTML,
    tagListFn = shiny::tagList,
    modalButtonFn = shiny::modalButton,
    actionButtonFn = shiny::actionButton
) {
    showModalFn(modalDialogFn(
        title = "Confirm Reset",
        htmlFn(paste0(
            "<p>This will revert <strong>", resetScope,
            "</strong> to their initial state. This action cannot be undone.</p>"
        )),
        footer = tagListFn(
            modalButtonFn("Cancel"),
            actionButtonFn(nsFn("confirm_reset"), "Reset", class = "btn-danger")
        ),
        easyClose = TRUE
    ))
}

registerProtDesignResetRequestObserver <- function(
    input,
    session,
    observeEventFn = shiny::observeEvent,
    showResetConfirmationModalFn = showProtDesignResetConfirmationModal
) {
    observeEventFn(input$reset_changes, {
        showResetConfirmationModalFn(
            resetScope = input$reset_scope,
            nsFn = session$ns
        )
    })
}

runProtDesignResetConfirmationObserverShell <- function(
    scope,
    initialState,
    designMatrix,
    dataClnReactive,
    removedSamples,
    factors,
    groups,
    contrasts,
    session,
    applyResetStateFn = applyProtDesignResetState,
    removeModalFn = shiny::removeModal,
    showNotificationFn = shiny::showNotification,
    updateTextInputFn = shiny::updateTextInput
) {
    applyResetStateFn(
        scope = scope,
        initialState = initialState,
        designMatrixGetter = designMatrix,
        designMatrixSetter = designMatrix,
        dataClnSetter = dataClnReactive,
        currentRemovedSamples = removedSamples(),
        removedSamplesSetter = removedSamples,
        factorsSetter = factors,
        groupsSetter = groups,
        contrastsSetter = contrasts,
        updateFormulaFn = function(value) {
            updateTextInputFn(session, "formula_string", value = value)
        }
    )

    removeModalFn()
    showNotificationFn(paste("Reset of", scope, "completed."), type = "message")

    invisible(NULL)
}

registerProtDesignResetConfirmationObserver <- function(
    input,
    initialState,
    designMatrix,
    dataClnReactive,
    removedSamples,
    factors,
    groups,
    contrasts,
    session,
    observeEventFn = shiny::observeEvent,
    runResetConfirmationObserverShell = runProtDesignResetConfirmationObserverShell
) {
    observeEventFn(input$confirm_reset, {
        runResetConfirmationObserverShell(
            scope = input$reset_scope,
            initialState = initialState(),
            designMatrix = designMatrix,
            dataClnReactive = dataClnReactive,
            removedSamples = removedSamples,
            factors = factors,
            groups = groups,
            contrasts = contrasts,
            session = session
        )
    })
}

applyProtDesignResetState <- function(
    scope,
    initialState,
    designMatrixGetter,
    designMatrixSetter,
    dataClnSetter,
    currentRemovedSamples,
    removedSamplesSetter,
    factorsSetter,
    groupsSetter,
    contrastsSetter,
    updateFormulaFn = function(value) invisible(value),
    applyRemovedSamplesUpdateFn = applyProtDesignRemovedSamplesUpdate
) {
    if (scope == "all" || scope == "sample_names") {
        designMatrixSetter(initialState$design_matrix)
        dataClnSetter(initialState$data_cln)
    }
    if (scope == "all" || scope == "removed_samples") {
        removedSamplesSetter(applyRemovedSamplesUpdateFn(
            currentRemovedSamples = currentRemovedSamples,
            reset = TRUE
        ))
    }
    if (scope == "all" || scope == "factors") {
        factorsSetter(initialState$factors)

        currentMatrix <- designMatrixGetter()
        currentMatrix$factor1 <- NA_character_
        currentMatrix$factor2 <- NA_character_
        currentMatrix$factor3 <- NA_character_
        currentMatrix$group <- NA_character_
        currentMatrix$tech_reps <- NA_integer_

        # Preserve import-time columns such as batch while clearing assignments.
        designMatrixSetter(currentMatrix)
        groupsSetter(initialState$groups)
    }
    if (scope == "all" || scope == "contrasts") {
        contrastsSetter(initialState$contrasts)
    }
    if (scope == "all" || scope == "formula") {
        updateFormulaFn(initialState$formula)
    }

    invisible(NULL)
}

applyProtDesignRemovedSamplesUpdate <- function(
    currentRemovedSamples,
    selectedSamples = character(0),
    reset = FALSE
) {
    if (reset) {
        return(character(0))
    }

    unique(c(currentRemovedSamples, selectedSamples))
}

buildProtDesignInitialState <- function(
    dataTbl,
    configList,
    columnMapping
) {
    dataCln <- dataTbl
    quantityColumns <- unique(Filter(
        Negate(is.null),
        list(
            columnMapping$norm_quantity_col,
            columnMapping$raw_quantity_col,
            columnMapping$quantity_col
        )
    ))

    for (columnName in quantityColumns) {
        if (columnName %in% names(dataCln)) {
            dataCln <- dataCln |>
                dplyr::mutate(
                    !!rlang::sym(columnName) := as.numeric(!!rlang::sym(columnName))
                )
        }
    }

    designMatrixRaw <- tibble::tibble(
        Run = dataCln |>
            dplyr::distinct(Run) |>
            dplyr::select(Run) |>
            dplyr::pull(Run) |>
            gtools::mixedsort(),
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
        data_cln = dataCln,
        groups = unique(
            designMatrixRaw$group[
                !is.na(designMatrixRaw$group) & designMatrixRaw$group != ""
            ]
        ),
        factors = character(0),
        formula = configList[["deAnalysisParameters"]][["formula_string"]],
        contrasts = data.frame(
            contrast_name = character(),
            numerator = character(),
            denominator = character(),
            stringsAsFactors = FALSE
        )
    )
}

