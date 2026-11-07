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

#' The proteomics builder's reset confirmation (its text adds "This action cannot be undone."); the
#' `dia_long` adapter's `reset_modal` (KIT-103 W6). The shared builder shows it.
#' @noRd
protDesignResetModal <- function(resetScope, ns) {
    shiny::modalDialog(
        title = "Confirm Reset",
        shiny::HTML(paste0(
            "<p>This will revert <strong>", resetScope,
            "</strong> to their initial state. This action cannot be undone.</p>"
        )),
        footer = shiny::tagList(
            shiny::modalButton("Cancel"),
            shiny::actionButton(ns("confirm_reset"), "Reset", class = "btn-danger")
        ),
        easyClose = TRUE
    )
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
