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
    manager <- getDIAActiveManager()
    if (inherits(manager, "DIAStateAdapter") && inherits(dataCln, "spec_tbl_df") &&
        length(getDIAImportMetadata(manager)$numeric_columns)) {
        # The legacy numeric mutate also removes reader-only metadata. Defer
        # that metadata-only reconstruction to this required table consumer;
        # the imported reference retains its validated parser diagnostics.
        dataCln <- dplyr::mutate(dataCln)
    }
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

    if (inherits(getDIAActiveManager(), "DIAStateAdapter")) {
        finalResult$data_reference <- saveDIAOwnedValue(finalResult$data_cln, "design-assigned")
        finalResult["data_cln"] <- list(NULL)
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
        withDIAStateScope(getDIAActiveManager(), function() {
        runSaveResultsObserverShell(
            designMatrix = designMatrix(),
            currentRemovedSamples = removedSamples(),
            dataCln = getDIAReactiveReference(dataCln),
            contrastData = contrastData(),
            configList = configList(),
            formulaString = input$formula_string,
            resultSetter = resultRv
        )
        })
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
        dataClnSetter(initialState$data_reference %||% initialState$data_cln)
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
