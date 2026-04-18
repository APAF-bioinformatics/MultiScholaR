renameLipidDesignDesignMatrixRuns <- function(designMatrix, renameMap) {
    updatedDesignMatrix <- designMatrix

    for (oldName in names(renameMap)) {
        updatedDesignMatrix$Run[updatedDesignMatrix$Run == oldName] <- unname(renameMap[[oldName]])
    }

    updatedDesignMatrix
}

renameLipidDesignAssayColumns <- function(assayList, renameMap) {
    updatedAssayList <- assayList

    for (assayName in names(updatedAssayList)) {
        assayColumns <- names(updatedAssayList[[assayName]])
        renamedColumns <- renameMap[assayColumns]
        columnsToRename <- !is.na(renamedColumns)
        assayColumns[columnsToRename] <- unname(renamedColumns[columnsToRename])
        names(updatedAssayList[[assayName]]) <- assayColumns
    }

    updatedAssayList
}

renameLipidDesignTrackedSamples <- function(sampleNames, renameMap) {
    updatedSampleNames <- sampleNames
    renamedSamples <- renameMap[updatedSampleNames]
    samplesToRename <- !is.na(renamedSamples)
    updatedSampleNames[samplesToRename] <- unname(renamedSamples[samplesToRename])

    updatedSampleNames
}

applyLipidDesignSampleRenameMap <- function(
    designMatrix,
    assayList,
    sampleNames,
    renameMap,
    renameDesignMatrixRunsFn = renameLipidDesignDesignMatrixRuns,
    renameAssayColumnsFn = renameLipidDesignAssayColumns,
    renameTrackedSamplesFn = renameLipidDesignTrackedSamples
) {
    list(
        designMatrix = renameDesignMatrixRunsFn(designMatrix, renameMap)
        , assayList = renameAssayColumnsFn(assayList, renameMap)
        , sampleNames = renameTrackedSamplesFn(sampleNames, renameMap)
    )
}

buildLipidDesignBulkRenameMap <- function(
    samplesToTransform,
    transformMode,
    rangeStart = NULL,
    rangeEnd = NULL,
    extractExperimentFn = extract_experiment,
    reqFn = shiny::req
) {
    renamedSamples <- vapply(samplesToTransform, function(sampleName) {
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

        stop("Unsupported transform mode.")
    }, character(1))

    stats::setNames(renamedSamples, samplesToTransform)
}

registerLipidDesignSampleRenameShells <- function(
    input,
    session,
    designMatrix,
    dataClnReactive,
    sampleNamesReactive,
    observeEventFn = shiny::observeEvent,
    reqFn = shiny::req,
    updateTextInputFn = shiny::updateTextInput,
    buildBulkRenameMapFn = buildLipidDesignBulkRenameMap,
    applySampleRenameMapFn = applyLipidDesignSampleRenameMap
) {
    observeEventFn(input$rename_sample, {
        reqFn(input$sample_to_rename, input$new_sample_name)

        renamedState <- applySampleRenameMapFn(
            designMatrix = designMatrix()
            , assayList = dataClnReactive()
            , sampleNames = sampleNamesReactive()
            , renameMap = stats::setNames(input$new_sample_name, input$sample_to_rename)
        )

        designMatrix(renamedState$designMatrix)
        dataClnReactive(renamedState$assayList)
        sampleNamesReactive(renamedState$sampleNames)

        updateTextInputFn(session, "new_sample_name", value = "")
    })

    observeEventFn(input$bulk_rename, {
        reqFn(input$samples_to_transform)

        renameMap <- buildBulkRenameMapFn(
            samplesToTransform = input$samples_to_transform
            , transformMode = input$transform_mode
            , rangeStart = input$range_start
            , rangeEnd = input$range_end
        )

        renamedState <- applySampleRenameMapFn(
            designMatrix = designMatrix()
            , assayList = dataClnReactive()
            , sampleNames = sampleNamesReactive()
            , renameMap = renameMap
        )

        designMatrix(renamedState$designMatrix)
        dataClnReactive(renamedState$assayList)
        sampleNamesReactive(renamedState$sampleNames)
    })

    invisible(NULL)
}

appendLipidDesignFactorName <- function(currentFactors, newFactorName) {
    trimmedFactorName <- trimws(newFactorName)

    if (trimmedFactorName == "" || trimmedFactorName %in% currentFactors) {
        return(currentFactors)
    }

    c(currentFactors, trimmedFactorName)
}

buildLipidDesignMetadataReplicateNumbers <- function(selectedRuns, replicateStart) {
    if (!is.null(replicateStart) && !is.na(replicateStart)) {
        return(seq(replicateStart, length.out = length(selectedRuns)))
    }

    NA_integer_
}

buildLipidDesignMetadataGroupName <- function(
    factor1Selection,
    factor2Selection = NULL,
    factor3Selection = NULL
) {
    factorParts <- c()

    if (!is.null(factor1Selection) && factor1Selection != "") {
        factorParts <- c(factorParts, factor1Selection)
    }
    if (!is.null(factor2Selection) && factor2Selection != "") {
        factorParts <- c(factorParts, factor2Selection)
    }
    if (!is.null(factor3Selection) && factor3Selection != "") {
        factorParts <- c(factorParts, factor3Selection)
    }

    if (length(factorParts) == 0) {
        return(NA_character_)
    }

    paste(factorParts, collapse = "_")
}

listLipidDesignAssignedGroups <- function(designMatrix) {
    unique(designMatrix$group[!is.na(designMatrix$group) & designMatrix$group != ""])
}

applyLipidDesignMetadataAssignment <- function(
    designMatrix,
    selectedRuns,
    factor1Selection,
    factor2Selection = NULL,
    factor3Selection = NULL,
    replicateStart = NULL,
    buildMetadataReplicateNumbersFn = buildLipidDesignMetadataReplicateNumbers,
    buildMetadataGroupNameFn = buildLipidDesignMetadataGroupName
) {
    updatedDesignMatrix <- designMatrix
    selectedIndices <- which(updatedDesignMatrix$Run %in% selectedRuns)
    replicateNumbers <- buildMetadataReplicateNumbersFn(selectedRuns, replicateStart)

    updatedDesignMatrix$factor1[selectedIndices] <- factor1Selection
    updatedDesignMatrix$factor2[selectedIndices] <- factor2Selection
    updatedDesignMatrix$factor3[selectedIndices] <- if (!is.null(factor3Selection) && factor3Selection != "") {
        factor3Selection
    } else {
        NA_character_
    }
    updatedDesignMatrix$replicates[selectedIndices] <- replicateNumbers
    updatedDesignMatrix$group[selectedIndices] <- buildMetadataGroupNameFn(
        factor1Selection = factor1Selection
        , factor2Selection = factor2Selection
        , factor3Selection = factor3Selection
    )

    updatedDesignMatrix
}

registerLipidDesignMetadataAssignmentShells <- function(
    input,
    session,
    designMatrix,
    factors,
    groups,
    observeEventFn = shiny::observeEvent,
    reqFn = shiny::req,
    updateTextInputFn = shiny::updateTextInput,
    appendFactorNameFn = appendLipidDesignFactorName,
    applyMetadataAssignmentFn = applyLipidDesignMetadataAssignment,
    listAssignedGroupsFn = listLipidDesignAssignedGroups
) {
    observeEventFn(input$add_factor, {
        reqFn(input$new_factor)
        factors(appendFactorNameFn(factors(), input$new_factor))
        updateTextInputFn(session, "new_factor", value = "")
    })

    observeEventFn(input$assign_metadata, {
        reqFn(input$selected_runs, input$factor1_select)

        updatedDesignMatrix <- applyMetadataAssignmentFn(
            designMatrix = designMatrix()
            , selectedRuns = input$selected_runs
            , factor1Selection = input$factor1_select
            , factor2Selection = input$factor2_select
            , factor3Selection = input$factor3_select
            , replicateStart = input$replicate_start
        )

        designMatrix(updatedDesignMatrix)
        groups(listAssignedGroupsFn(updatedDesignMatrix))
    })

    invisible(NULL)
}

applyLipidDesignTechRepAssignment <- function(
    designMatrix,
    techRepSamples,
    assignmentMode,
    manualReplicateNumber = NULL
) {
    updatedDesignMatrix <- designMatrix

    if (length(techRepSamples) < 2) {
        return(list(
            ok = FALSE
            , type = "warning"
            , message = "Please select at least two samples."
            , designMatrix = updatedDesignMatrix
        ))
    }

    selectedIndices <- which(updatedDesignMatrix$Run %in% techRepSamples)
    selectedGroups <- unique(updatedDesignMatrix$group[selectedIndices])

    if (length(selectedGroups) > 1 || any(is.na(selectedGroups))) {
        return(list(
            ok = FALSE
            , type = "warning"
            , message = "All selected samples must belong to the same group."
            , designMatrix = updatedDesignMatrix
        ))
    }

    if (assignmentMode == "lowest") {
        baseReplicateNumber <- min(updatedDesignMatrix$replicates[selectedIndices], na.rm = TRUE)
    } else if (assignmentMode == "first") {
        baseReplicateNumber <- updatedDesignMatrix$replicates[
            updatedDesignMatrix$Run == techRepSamples[1]
        ]
    } else if (assignmentMode == "manual") {
        baseReplicateNumber <- manualReplicateNumber
    } else {
        stop("Unsupported technical replicate assignment mode.")
    }

    sampleGroup <- selectedGroups[1]
    originalReplicateNumbers <- unique(updatedDesignMatrix$replicates[selectedIndices])
    numReplicatesConsolidated <- length(originalReplicateNumbers) - 1
    highestConsolidatedReplicate <- max(originalReplicateNumbers, na.rm = TRUE)

    samplesToAdjust <- which(
        updatedDesignMatrix$group == sampleGroup &
        !is.na(updatedDesignMatrix$group) &
        updatedDesignMatrix$replicates > highestConsolidatedReplicate &
        !updatedDesignMatrix$Run %in% techRepSamples
    )

    if (length(samplesToAdjust) > 0) {
        updatedDesignMatrix$replicates[samplesToAdjust] <-
            updatedDesignMatrix$replicates[samplesToAdjust] - numReplicatesConsolidated
    }

    updatedDesignMatrix$replicates[selectedIndices] <- baseReplicateNumber
    techRepMapping <- stats::setNames(seq_along(techRepSamples), techRepSamples)
    updatedDesignMatrix$tech_reps[selectedIndices] <-
        techRepMapping[updatedDesignMatrix$Run[selectedIndices]]

    list(
        ok = TRUE
        , type = "message"
        , message = paste("Assigned", length(techRepSamples), "samples as technical replicates")
        , designMatrix = updatedDesignMatrix
    )
}

registerLipidDesignTechRepAssignmentShells <- function(
    input,
    designMatrix,
    observeEventFn = shiny::observeEvent,
    reqFn = shiny::req,
    showNotificationFn = shiny::showNotification,
    applyTechRepAssignmentFn = applyLipidDesignTechRepAssignment
) {
    observeEventFn(input$assign_tech_reps, {
        reqFn(input$tech_rep_samples)

        assignmentResult <- applyTechRepAssignmentFn(
            designMatrix = designMatrix()
            , techRepSamples = input$tech_rep_samples
            , assignmentMode = input$tech_rep_assignment_mode
            , manualReplicateNumber = input$manual_replicate_number
        )

        if (!isTRUE(assignmentResult$ok)) {
            showNotificationFn(assignmentResult$message, type = assignmentResult$type)
            return()
        }

        designMatrix(assignmentResult$designMatrix)
        showNotificationFn(assignmentResult$message, type = assignmentResult$type)
    })

    invisible(NULL)
}

appendLipidDesignContrast <- function(currentContrasts, group1, group2) {
    updatedContrasts <- currentContrasts

    if (is.null(group1) || is.null(group2) || group1 == "" || group2 == "" || group1 == group2) {
        return(updatedContrasts)
    }

    contrastName <- paste0(group1, ".vs.", group2)

    if (contrastName %in% updatedContrasts$contrast_name) {
        return(updatedContrasts)
    }

    newContrast <- data.frame(
        contrast_name = contrastName
        , numerator = group1
        , denominator = group2
        , stringsAsFactors = FALSE
    )

    rbind(updatedContrasts, newContrast)
}

registerLipidDesignContrastManagementShells <- function(
    input,
    contrasts,
    observeEventFn = shiny::observeEvent,
    reqFn = shiny::req,
    appendContrastFn = appendLipidDesignContrast
) {
    observeEventFn(input$add_contrast, {
        reqFn(input$contrast_group1, input$contrast_group2)

        contrasts(appendContrastFn(
            currentContrasts = contrasts()
            , group1 = input$contrast_group1
            , group2 = input$contrast_group2
        ))
    })

    invisible(NULL)
}

appendLipidDesignRemovedSamples <- function(currentRemovedSamples, selectedSamples) {
    unique(c(currentRemovedSamples, selectedSamples))
}

registerLipidDesignSampleRemovalShells <- function(
    input,
    removedSamples,
    session,
    observeEventFn = shiny::observeEvent,
    reqFn = shiny::req,
    updateSelectizeInputFn = shiny::updateSelectizeInput,
    showNotificationFn = shiny::showNotification,
    appendRemovedSamplesFn = appendLipidDesignRemovedSamples
) {
    observeEventFn(input$remove_samples, {
        reqFn(input$samples_to_remove)

        samplesToRemove <- input$samples_to_remove
        removedSamples(appendRemovedSamplesFn(
            currentRemovedSamples = removedSamples()
            , selectedSamples = samplesToRemove
        ))

        updateSelectizeInputFn(session, "samples_to_remove", selected = "")
        showNotificationFn(
            paste("Removed", length(samplesToRemove), "sample(s) from analysis.")
            , type = "message"
        )
    })

    invisible(NULL)
}

