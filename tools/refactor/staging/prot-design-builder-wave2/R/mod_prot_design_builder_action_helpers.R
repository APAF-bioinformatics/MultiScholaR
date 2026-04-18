buildProtDesignReplicateInputs <- function(selectedRuns, nsFn = identity, numericInputFn = shiny::numericInput) {
    numericInputFn(
        nsFn("replicate_start"),
        paste("Starting replicate number for", length(selectedRuns), "selected runs:"),
        value = 1,
        min = 1
    )
}

formatProtDesignRangePreview <- function(
    selectedSamples,
    rangeStart,
    rangeEnd,
    extractExperimentFn = extract_experiment
) {
    firstSample <- selectedSamples[[1]]

    tryCatch({
        previewResult <- extractExperimentFn(
            firstSample,
            mode = "range",
            start = rangeStart,
            end = rangeEnd
        )
        paste0("\"", firstSample, "\" -> \"", previewResult, "\"")
    }, error = function(e) {
        paste("Error:", e$message)
    })
}

transformProtDesignSampleNames <- function(
    selectedSamples,
    transformMode,
    rangeStart = NULL,
    rangeEnd = NULL,
    extractExperimentFn = extract_experiment,
    reqFn = shiny::req
) {
    vapply(selectedSamples, function(sampleName) {
        if (transformMode == "range") {
            reqFn(rangeStart, rangeEnd)
            return(extractExperimentFn(
                sampleName,
                mode = "range",
                start = rangeStart,
                end = rangeEnd
            ))
        }

        if (transformMode == "before_underscore") {
            return(extractExperimentFn(sampleName, mode = "start"))
        }

        if (transformMode == "after_underscore") {
            return(extractExperimentFn(sampleName, mode = "end"))
        }

        stop("Unsupported transform mode: ", transformMode, call. = FALSE)
    }, character(1), USE.NAMES = FALSE)
}

applyProtDesignBulkRenameUpdates <- function(
    designMatrix,
    dataCln,
    selectedSamples,
    newNames
) {
    updatedDesignMatrix <- designMatrix
    updatedDataCln <- dataCln

    for (i in seq_along(selectedSamples)) {
        originalName <- selectedSamples[[i]]
        newName <- newNames[[i]]
        updatedDesignMatrix$Run[updatedDesignMatrix$Run == originalName] <- newName
        updatedDataCln$Run[updatedDataCln$Run == originalName] <- newName
    }

    list(
        designMatrix = updatedDesignMatrix,
        dataCln = updatedDataCln
    )
}

applyProtDesignSingleRenameUpdate <- function(
    designMatrix,
    dataCln,
    originalName,
    newName
) {
    updatedDesignMatrix <- designMatrix
    updatedDataCln <- dataCln

    updatedDesignMatrix$Run[updatedDesignMatrix$Run == originalName] <- newName
    updatedDataCln$Run[updatedDataCln$Run == originalName] <- newName

    list(
        designMatrix = updatedDesignMatrix,
        dataCln = updatedDataCln
    )
}

registerProtDesignRenameSampleObserver <- function(
    input,
    designMatrix,
    dataCln,
    session,
    observeEventFn = shiny::observeEvent,
    reqFn = shiny::req,
    applySingleRenameUpdateFn = applyProtDesignSingleRenameUpdate,
    updateTextInputFn = shiny::updateTextInput
) {
    observeEventFn(input$rename_sample, {
        reqFn(input$sample_to_rename, input$new_sample_name)

        renameUpdate <- applySingleRenameUpdateFn(
            designMatrix = designMatrix(),
            dataCln = dataCln(),
            originalName = input$sample_to_rename,
            newName = input$new_sample_name
        )

        designMatrix(renameUpdate$designMatrix)
        dataCln(renameUpdate$dataCln)
        updateTextInputFn(session, "new_sample_name", value = "")
    })
}

registerProtDesignBulkRenameObserver <- function(
    input,
    designMatrix,
    dataCln,
    observeEventFn = shiny::observeEvent,
    reqFn = shiny::req,
    transformSampleNamesFn = transformProtDesignSampleNames,
    applyBulkRenameUpdatesFn = applyProtDesignBulkRenameUpdates
) {
    observeEventFn(input$bulk_rename, {
        reqFn(input$samples_to_transform)

        renameUpdates <- applyBulkRenameUpdatesFn(
            designMatrix = designMatrix(),
            dataCln = dataCln(),
            selectedSamples = input$samples_to_transform,
            newNames = transformSampleNamesFn(
                selectedSamples = input$samples_to_transform,
                transformMode = input$transform_mode,
                rangeStart = input$range_start,
                rangeEnd = input$range_end
            )
        )

        designMatrix(renameUpdates$designMatrix)
        dataCln(renameUpdates$dataCln)
    })
}

applyProtDesignFactorAppendReset <- function(
    currentFactors,
    newFactorInput
) {
    newFactorName <- trimws(newFactorInput)
    updatedFactors <- currentFactors

    if (newFactorName != "" && !newFactorName %in% currentFactors) {
        updatedFactors <- c(currentFactors, newFactorName)
    }

    list(
        factors = updatedFactors,
        newFactorValue = ""
    )
}

registerProtDesignAddFactorObserver <- function(
    input,
    factors,
    session,
    observeEventFn = shiny::observeEvent,
    reqFn = shiny::req,
    applyFactorAppendResetFn = applyProtDesignFactorAppendReset,
    updateTextInputFn = shiny::updateTextInput
) {
    observeEventFn(input$add_factor, {
        reqFn(input$new_factor)

        factorUpdate <- applyFactorAppendResetFn(
            currentFactors = factors(),
            newFactorInput = input$new_factor
        )

        factors(factorUpdate$factors)
        updateTextInputFn(session, "new_factor", value = factorUpdate$newFactorValue)
    })
}

registerProtDesignAssignMetadataObserver <- function(
    input,
    designMatrix,
    groups,
    observeEventFn = shiny::observeEvent,
    reqFn = shiny::req,
    seqFn = seq.int
) {
    observeEventFn(input$assign_metadata, {
        reqFn(input$selected_runs, input$factor1_select)

        currentMatrix <- designMatrix()

        replicateNumbers <- if (!is.null(input$replicate_start) && !is.na(input$replicate_start)) {
            seqFn(input$replicate_start, length.out = length(input$selected_runs))
        } else {
            NA_integer_
        }

        selectedIndices <- which(currentMatrix$Run %in% input$selected_runs)
        currentMatrix$factor1[selectedIndices] <- input$factor1_select
        currentMatrix$factor2[selectedIndices] <- input$factor2_select
        currentMatrix$factor3[selectedIndices] <- if (!is.null(input$factor3_select) && input$factor3_select != "") {
            input$factor3_select
        } else {
            NA_character_
        }
        currentMatrix$replicates[selectedIndices] <- replicateNumbers

        factorParts <- c()
        if (!is.null(input$factor1_select) && input$factor1_select != "") {
            factorParts <- c(factorParts, input$factor1_select)
        }
        if (!is.null(input$factor2_select) && input$factor2_select != "") {
            factorParts <- c(factorParts, input$factor2_select)
        }
        if (!is.null(input$factor3_select) && input$factor3_select != "") {
            factorParts <- c(factorParts, input$factor3_select)
        }

        groupName <- if (length(factorParts) > 0) {
            paste(factorParts, collapse = "_")
        } else {
            NA_character_
        }
        currentMatrix$group[selectedIndices] <- groupName

        designMatrix(currentMatrix)

        uniqueGroups <- unique(currentMatrix$group[!is.na(currentMatrix$group) & currentMatrix$group != ""])
        groups(uniqueGroups)
    })
}

registerProtDesignAssignTechRepsObserver <- function(
    input,
    designMatrix,
    observeEventFn = shiny::observeEvent,
    reqFn = shiny::req,
    showNotificationFn = shiny::showNotification
) {
    observeEventFn(input$assign_tech_reps, {
        reqFn(input$tech_rep_samples)
        currentMatrix <- designMatrix()

        if (length(input$tech_rep_samples) < 2) {
            showNotificationFn(
                "Please select at least two samples to assign as technical replicates.",
                type = "warning"
            )
            return()
        }

        selectedIndices <- which(currentMatrix$Run %in% input$tech_rep_samples)
        selectedGroups <- unique(currentMatrix$group[selectedIndices])
        if (length(selectedGroups) > 1 || any(is.na(selectedGroups))) {
            showNotificationFn(
                "All selected samples must belong to the same group. Please assign metadata first.",
                type = "warning"
            )
            return()
        }

        if (input$tech_rep_assignment_mode == "lowest") {
            baseReplicateNumber <- min(currentMatrix$replicates[selectedIndices], na.rm = TRUE)
        } else if (input$tech_rep_assignment_mode == "first") {
            baseReplicateNumber <- currentMatrix$replicates[currentMatrix$Run == input$tech_rep_samples[1]]
        } else if (input$tech_rep_assignment_mode == "manual") {
            baseReplicateNumber <- input$manual_replicate_number
        }

        sampleGroup <- selectedGroups[1]
        originalReplicateNumbers <- unique(currentMatrix$replicates[selectedIndices])
        numReplicatesConsolidated <- length(originalReplicateNumbers) - 1
        highestConsolidatedReplicate <- max(originalReplicateNumbers, na.rm = TRUE)

        samplesToAdjust <- which(
            currentMatrix$group == sampleGroup &
            !is.na(currentMatrix$group) &
            currentMatrix$replicates > highestConsolidatedReplicate &
            !currentMatrix$Run %in% input$tech_rep_samples
        )

        if (length(samplesToAdjust) > 0) {
            currentMatrix$replicates[samplesToAdjust] <-
                currentMatrix$replicates[samplesToAdjust] - numReplicatesConsolidated
        }

        currentMatrix$replicates[selectedIndices] <- baseReplicateNumber

        techRepMapping <- setNames(seq_along(input$tech_rep_samples), input$tech_rep_samples)
        currentMatrix$tech_reps[selectedIndices] <- techRepMapping[currentMatrix$Run[selectedIndices]]

        designMatrix(currentMatrix)

        showNotificationFn(
            paste(
                "Assigned", length(input$tech_rep_samples),
                "samples as technical replicates with biological replicate number",
                baseReplicateNumber, "and adjusted subsequent replicate numbers"
            ),
            type = "message"
        )
    })
}

applyProtDesignContrastAppend <- function(
    currentContrasts,
    group1,
    group2
) {
    if (
        is.null(group1) || is.null(group2) ||
        group1 == "" || group2 == "" ||
        identical(group1, group2)
    ) {
        return(currentContrasts)
    }

    contrastName <- paste0(group1, ".vs.", group2)
    if (!is.null(currentContrasts) && contrastName %in% currentContrasts$contrast_name) {
        return(currentContrasts)
    }

    newContrast <- data.frame(
        contrast_name = contrastName,
        numerator = group1,
        denominator = group2,
        stringsAsFactors = FALSE
    )

    if (is.null(currentContrasts) || nrow(currentContrasts) == 0) {
        return(newContrast)
    }

    rbind(currentContrasts, newContrast)
}

registerProtDesignAddContrastObserver <- function(
    input,
    contrasts,
    observeEventFn = shiny::observeEvent,
    reqFn = shiny::req,
    appendContrastFn = applyProtDesignContrastAppend
) {
    observeEventFn(input$add_contrast, {
        reqFn(input$contrast_group1, input$contrast_group2)
        updatedContrasts <- appendContrastFn(
            currentContrasts = contrasts(),
            group1 = input$contrast_group1,
            group2 = input$contrast_group2
        )
        contrasts(updatedContrasts)
    })
}

registerProtDesignRemoveSamplesObserver <- function(
    input,
    removedSamples,
    session,
    observeEventFn = shiny::observeEvent,
    reqFn = shiny::req,
    applyRemovedSamplesUpdateFn = applyProtDesignRemovedSamplesUpdate,
    updateSelectizeInputFn = shiny::updateSelectizeInput,
    showNotificationFn = shiny::showNotification
) {
    observeEventFn(input$remove_samples, {
        reqFn(input$samples_to_remove)

        samplesToRemove <- input$samples_to_remove
        removedSamples(applyRemovedSamplesUpdateFn(
            currentRemovedSamples = removedSamples(),
            selectedSamples = samplesToRemove
        ))

        updateSelectizeInputFn(session, "samples_to_remove", selected = "")
        showNotificationFn(
            paste("Removed", length(samplesToRemove), "sample(s) from analysis."),
            type = "message"
        )
    })
}

