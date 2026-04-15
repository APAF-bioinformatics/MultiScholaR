formatLipidDesignTechRepSummary <- function(designMatrix) {
    techRepGroups <- designMatrix |>
        dplyr::filter(!is.na(tech_reps)) |>
        dplyr::group_by(group, replicates) |>
        dplyr::summarise(
            samples = paste(Run, collapse = ", ")
            , techRepNumbers = paste(tech_reps, collapse = ", ")
            , .groups = "drop"
        )

    if (nrow(techRepGroups) == 0) {
        return("No technical replicates assigned yet.")
    }

    outputText <- techRepGroups |>
        dplyr::mutate(
            formatted = sprintf(
                "Group: %s, Biological Replicate: %s\n  Samples: %s\n  Technical Replicates: %s"
                , group, replicates, samples, techRepNumbers
            )
        ) |>
        dplyr::pull(formatted)

    paste(outputText, collapse = "\n\n")
}

formatLipidDesignRemovedSamplesDisplay <- function(removedSamples) {
    if (length(removedSamples) == 0) {
        return("No samples have been removed.")
    }

    paste0(
        "Removed ", length(removedSamples), " sample(s):\n"
        , paste(gtools::mixedsort(removedSamples), collapse = "\n")
    )
}

formatLipidDesignContrastFactorsInfo <- function(formulaString) {
    formulaUsesGroupPrefix <- grepl("~ *0 *\\+ *group", formulaString)
    if (formulaUsesGroupPrefix) {
        "Note: Contrasts will use 'group' prefix based on formula: ~ 0 + group"
    } else {
        "Note: Contrasts will use group names as-is"
    }
}

registerLipidDesignSummaryOutputShells <- function(
    output,
    designMatrix,
    removedSamples,
    formulaString,
    renderTextFn = shiny::renderText,
    reqFn = shiny::req,
    formatTechRepSummaryFn = formatLipidDesignTechRepSummary,
    formatRemovedSamplesDisplayFn = formatLipidDesignRemovedSamplesDisplay,
    formatContrastFactorsInfoFn = formatLipidDesignContrastFactorsInfo
) {
    output$tech_rep_summary <- renderTextFn({
        currentDesignMatrix <- designMatrix()
        reqFn(currentDesignMatrix)

        formatTechRepSummaryFn(currentDesignMatrix)
    })

    output$removed_samples_display <- renderTextFn({
        formatRemovedSamplesDisplayFn(removedSamples())
    })

    output$contrast_factors_info <- renderTextFn({
        formatContrastFactorsInfoFn(formulaString())
    })

    output
}

formatLipidDesignAvailableFactorsDisplay <- function(currentFactors) {
    if (length(currentFactors) == 0) {
        return("No factors defined yet (use the 'Factors' tab).")
    }

    paste(currentFactors, collapse = ", ")
}

formatLipidDesignDefinedContrastLines <- function(contrastData, formulaString) {
    if (is.null(contrastData) || nrow(contrastData) == 0) {
        return(character(0))
    }

    formulaUsesGroupPrefix <- grepl("~ *0 *\\+ *group", formulaString)

    vapply(seq_len(nrow(contrastData)), function(index) {
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

        paste0(friendlyName, "=", contrastString)
    }, character(1))
}

formatLipidDesignRangePreview <- function(
    samplesToTransform,
    rangeStart,
    rangeEnd,
    extractExperimentFn = extract_experiment
) {
    firstSample <- samplesToTransform[1]

    tryCatch({
        previewResult <- extractExperimentFn(
            firstSample
            , mode = "range"
            , start = rangeStart
            , end = rangeEnd
        )
        paste0("\"", firstSample, "\" -> \"", previewResult, "\"")
    }, error = function(error) {
        paste("Error:", error$message)
    })
}

buildLipidDesignReplicateInputLabel <- function(selectedRuns) {
    paste("Starting replicate number for", length(selectedRuns), "selected samples:")
}

registerLipidDesignAdjacentOutputShells <- function(
    output,
    session,
    factors,
    contrasts,
    formulaString,
    samplesToTransform,
    rangeStart,
    rangeEnd,
    selectedRuns,
    renderUiFn = shiny::renderUI,
    renderTextFn = shiny::renderText,
    reqFn = shiny::req,
    paragraphFn = shiny::p,
    tagListFn = shiny::tagList,
    codeTagFn = shiny::tags$code,
    extractExperimentFn = extract_experiment,
    numericInputFn = shiny::numericInput,
    nsFn = session$ns,
    formatAvailableFactorsDisplayFn = formatLipidDesignAvailableFactorsDisplay,
    formatDefinedContrastLinesFn = formatLipidDesignDefinedContrastLines,
    formatRangePreviewFn = formatLipidDesignRangePreview,
    buildReplicateInputLabelFn = buildLipidDesignReplicateInputLabel
) {
    output$available_factors_display <- renderUiFn({
        paragraphFn(formatAvailableFactorsDisplayFn(factors()))
    })

    output$defined_contrasts_display <- renderUiFn({
        contrastLines <- formatDefinedContrastLinesFn(contrasts(), formulaString())
        if (length(contrastLines) == 0) {
            return(paragraphFn("No contrasts defined yet."))
        }

        tagListFn(lapply(contrastLines, function(fullFormat) {
            paragraphFn(codeTagFn(fullFormat))
        }))
    })

    output$range_preview <- renderTextFn({
        currentSamplesToTransform <- samplesToTransform()
        currentRangeStart <- rangeStart()
        currentRangeEnd <- rangeEnd()
        reqFn(currentSamplesToTransform)
        reqFn(currentRangeStart, currentRangeEnd)

        formatRangePreviewFn(
            currentSamplesToTransform,
            currentRangeStart,
            currentRangeEnd,
            extractExperimentFn = extractExperimentFn
        )
    })

    output$replicate_inputs <- renderUiFn({
        currentSelectedRuns <- selectedRuns()
        reqFn(currentSelectedRuns)

        numericInputFn(
            nsFn("replicate_start")
            , buildReplicateInputLabelFn(currentSelectedRuns)
            , value = 1
            , min = 1
        )
    })

    output
}

runLipidDesignSampleSelectionInputShell <- function(
    input,
    session,
    designMatrix,
    removedSamples,
    reqFn = shiny::req,
    mixedsortFn = gtools::mixedsort,
    isolateFn = shiny::isolate,
    updateSelectizeInputFn = shiny::updateSelectizeInput
) {
    currentDesignMatrix <- designMatrix()
    reqFn(currentDesignMatrix)

    allRuns <- mixedsortFn(currentDesignMatrix$Run)

    currentlyRemoved <- removedSamples()
    availableRuns <- allRuns[!allRuns %in% currentlyRemoved]

    isolateFn({
        selectedRename <- input$sample_to_rename
        selectedMeta <- input$selected_runs
        selectedTransform <- input$samples_to_transform
        selectedTech <- input$tech_rep_samples
        selectedRemove <- input$samples_to_remove

        selectedRename <- if (!is.null(selectedRename) && selectedRename %in% availableRuns) {
            selectedRename
        } else {
            ""
        }
        selectedMeta <- selectedMeta[selectedMeta %in% availableRuns]
        selectedTransform <- selectedTransform[selectedTransform %in% availableRuns]
        selectedTech <- selectedTech[selectedTech %in% availableRuns]
        selectedRemove <- selectedRemove[selectedRemove %in% availableRuns]
    })

    updateSelectizeInputFn(session, "sample_to_rename", choices = availableRuns, selected = selectedRename)
    updateSelectizeInputFn(session, "selected_runs", choices = availableRuns, selected = selectedMeta)
    updateSelectizeInputFn(session, "samples_to_transform", choices = availableRuns, selected = selectedTransform)
    updateSelectizeInputFn(session, "tech_rep_samples", choices = availableRuns, selected = selectedTech)
    updateSelectizeInputFn(session, "samples_to_remove", choices = availableRuns, selected = selectedRemove)

    invisible(NULL)
}

registerLipidDesignSampleSelectionInputShells <- function(
    input,
    session,
    designMatrix,
    removedSamples,
    observeFn = shiny::observe,
    runSampleSelectionInputShellFn = runLipidDesignSampleSelectionInputShell
) {
    observeFn({
        runSampleSelectionInputShellFn(
            input = input,
            session = session,
            designMatrix = designMatrix,
            removedSamples = removedSamples
        )
    })

    invisible(NULL)
}

runLipidDesignFactorGroupDropdownShell <- function(
    input,
    session,
    factors,
    groups,
    isolateFn = shiny::isolate,
    updateSelectInputFn = shiny::updateSelectInput
) {
    isolateFn({
        factor1Selection <- input$factor1_select
        factor2Selection <- input$factor2_select
        factor3Selection <- input$factor3_select
        contrastGroup1 <- input$contrast_group1
        contrastGroup2 <- input$contrast_group2
    })

    currentFactors <- factors()
    currentGroups <- groups()

    updateSelectInputFn(session, "factor1_select", choices = c("", currentFactors), selected = factor1Selection)
    updateSelectInputFn(session, "factor2_select", choices = c("", currentFactors), selected = factor2Selection)
    updateSelectInputFn(session, "factor3_select", choices = c("", currentFactors), selected = factor3Selection)
    updateSelectInputFn(session, "contrast_group1", choices = c("", currentGroups), selected = contrastGroup1)
    updateSelectInputFn(session, "contrast_group2", choices = c("", currentGroups), selected = contrastGroup2)

    invisible(NULL)
}

registerLipidDesignFactorGroupDropdownShells <- function(
    input,
    session,
    factors,
    groups,
    observeFn = shiny::observe,
    runFactorGroupDropdownShellFn = runLipidDesignFactorGroupDropdownShell
) {
    observeFn({
        runFactorGroupDropdownShellFn(
            input = input,
            session = session,
            factors = factors,
            groups = groups
        )
    })

    invisible(NULL)
}

buildLipidDesignActiveDataTable <- function(
    designMatrix,
    currentRemovedSamples = character(0)
) {
    designMatrix |>
        dplyr::filter(!Run %in% currentRemovedSamples)
}

registerLipidDesignDataTableShells <- function(
    output,
    proxyDataTable,
    designMatrix,
    removedSamples,
    renderDtFn = DT::renderDT,
    observeFn = shiny::observe,
    reqFn = shiny::req,
    replaceDataFn = DT::replaceData,
    buildActiveDataTableFn = buildLipidDesignActiveDataTable
) {
    output$data_table <- renderDtFn({
        currentDesignMatrix <- designMatrix()
        reqFn(currentDesignMatrix)

        buildActiveDataTableFn(
            designMatrix = currentDesignMatrix,
            currentRemovedSamples = removedSamples()
        )
    }
    , selection = "none"
    , options = list(pageLength = 10, scrollX = TRUE, server = FALSE)
    )

    observeFn({
        reqFn(proxyDataTable)
        currentDesignMatrix <- designMatrix()
        reqFn(currentDesignMatrix)

        replaceDataFn(
            proxyDataTable,
            buildActiveDataTableFn(
                designMatrix = currentDesignMatrix,
                currentRemovedSamples = removedSamples()
            ),
            resetPaging = FALSE
        )
    })

    invisible(NULL)
}

