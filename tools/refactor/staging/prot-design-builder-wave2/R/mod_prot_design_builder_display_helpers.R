formatProtDesignTechRepSummary <- function(designMatrix) {
    techRepGroups <- designMatrix |>
        dplyr::filter(!is.na(tech_reps)) |>
        dplyr::group_by(group, replicates) |>
        dplyr::summarise(
            samples = paste(Run, collapse = ", "),
            techRepNumbers = paste(tech_reps, collapse = ", "),
            .groups = "drop"
        )

    if (nrow(techRepGroups) == 0) {
        return("No technical replicates assigned yet.")
    }

    outputText <- techRepGroups |>
        dplyr::mutate(
            formatted = sprintf(
                "Group: %s, Biological Replicate: %s\n  Samples: %s\n  Technical Replicates: %s",
                group,
                replicates,
                samples,
                techRepNumbers
            )
        ) |>
        dplyr::pull(formatted)

    paste(outputText, collapse = "\n\n")
}

registerProtDesignTechRepSummaryOutput <- function(
    output,
    designMatrix,
    renderTextFn = shiny::renderText,
    reqFn = shiny::req,
    formatTechRepSummaryFn = formatProtDesignTechRepSummary
) {
    output$tech_rep_summary <- renderTextFn({
        dm <- designMatrix()
        reqFn(dm)

        formatTechRepSummaryFn(dm)
    })
}

registerProtDesignRemovedSamplesDisplayOutput <- function(
    output,
    removedSamples,
    renderTextFn = shiny::renderText,
    formatRemovedSamplesDisplayFn = formatProtDesignRemovedSamplesDisplay
) {
    output$removed_samples_display <- renderTextFn({
        formatRemovedSamplesDisplayFn(removedSamples())
    })
}

registerProtDesignContrastFactorsInfoOutput <- function(
    output,
    formulaString,
    renderTextFn = shiny::renderText,
    formatContrastFactorsInfoFn = formatProtDesignContrastFactorsInfo
) {
    output$contrast_factors_info <- renderTextFn({
        formatContrastFactorsInfoFn(formulaString())
    })
}

registerProtDesignDataTableOutput <- function(
    output,
    designMatrix,
    removedSamples,
    renderDTFn = DT::renderDT,
    reqFn = shiny::req,
    filterDesignMatrixFn = function(currentDesignMatrix, currentlyRemoved) {
        currentDesignMatrix |> dplyr::filter(!Run %in% currentlyRemoved)
    }
) {
    output$data_table <- renderDTFn({
            reqFn(designMatrix())
            currentDesignMatrix <- designMatrix()
            currentlyRemoved <- removedSamples()

            filterDesignMatrixFn(currentDesignMatrix, currentlyRemoved)
        },
        selection = "none",
        options = list(pageLength = 10, scrollX = TRUE, server = FALSE)
    )
}

registerProtDesignDataTableProxyRefreshObserver <- function(
    proxyDataTable,
    designMatrix,
    removedSamples,
    observeFn = shiny::observe,
    reqFn = shiny::req,
    filterDesignMatrixFn = function(currentDesignMatrix, currentlyRemoved) {
        currentDesignMatrix |> dplyr::filter(!Run %in% currentlyRemoved)
    },
    replaceDataFn = DT::replaceData
) {
    observeFn({
        reqFn(proxyDataTable)
        reqFn(designMatrix())

        currentDesignMatrix <- designMatrix()
        currentlyRemoved <- removedSamples()
        filteredDesignMatrix <- filterDesignMatrixFn(
            currentDesignMatrix,
            currentlyRemoved
        )

        replaceDataFn(
            proxyDataTable,
            filteredDesignMatrix,
            resetPaging = FALSE
        )
    })
}

registerProtDesignSampleSelectionSyncObserver <- function(
    input,
    session,
    designMatrix,
    removedSamples,
    observeFn = shiny::observe,
    reqFn = shiny::req,
    sortRunsFn = gtools::mixedsort,
    isolateFn = shiny::isolate,
    updateSelectizeInputFn = shiny::updateSelectizeInput
) {
    observeFn({
        reqFn(designMatrix())
        allRuns <- sortRunsFn(designMatrix()$Run)

        currentlyRemoved <- removedSamples()
        availableRuns <- allRuns[!allRuns %in% currentlyRemoved]

        preservedSelections <- isolateFn({
            list(
                rename = input$sample_to_rename,
                meta = input$selected_runs,
                transform = input$samples_to_transform,
                tech = input$tech_rep_samples,
                remove = input$samples_to_remove
            )
        })

        selectedRename <- if (
            !is.null(preservedSelections$rename) &&
            preservedSelections$rename %in% availableRuns
        ) {
            preservedSelections$rename
        } else {
            ""
        }

        selectedMeta <- preservedSelections$meta[
            preservedSelections$meta %in% availableRuns
        ]
        selectedTransform <- preservedSelections$transform[
            preservedSelections$transform %in% availableRuns
        ]
        selectedTech <- preservedSelections$tech[
            preservedSelections$tech %in% availableRuns
        ]
        selectedRemove <- preservedSelections$remove[
            preservedSelections$remove %in% availableRuns
        ]

        updateSelectizeInputFn(
            session,
            "sample_to_rename",
            choices = availableRuns,
            selected = selectedRename
        )
        updateSelectizeInputFn(
            session,
            "selected_runs",
            choices = availableRuns,
            selected = selectedMeta
        )
        updateSelectizeInputFn(
            session,
            "samples_to_transform",
            choices = availableRuns,
            selected = selectedTransform
        )
        updateSelectizeInputFn(
            session,
            "tech_rep_samples",
            choices = availableRuns,
            selected = selectedTech
        )
        updateSelectizeInputFn(
            session,
            "samples_to_remove",
            choices = availableRuns,
            selected = selectedRemove
        )
    })
}

registerProtDesignFactorGroupSyncObserver <- function(
    input,
    session,
    factors,
    groups,
    observeFn = shiny::observe,
    isolateFn = shiny::isolate,
    updateSelectInputFn = shiny::updateSelectInput
) {
    observeFn({
        preservedSelections <- isolateFn({
            list(
                factor1 = input$factor1_select,
                factor2 = input$factor2_select,
                factor3 = input$factor3_select,
                group1 = input$contrast_group1,
                group2 = input$contrast_group2
            )
        })

        currentFactors <- factors()
        currentGroups <- groups()

        updateSelectInputFn(
            session,
            "factor1_select",
            choices = c("", currentFactors),
            selected = preservedSelections$factor1
        )
        updateSelectInputFn(
            session,
            "factor2_select",
            choices = c("", currentFactors),
            selected = preservedSelections$factor2
        )
        updateSelectInputFn(
            session,
            "factor3_select",
            choices = c("", currentFactors),
            selected = preservedSelections$factor3
        )
        updateSelectInputFn(
            session,
            "contrast_group1",
            choices = c("", currentGroups),
            selected = preservedSelections$group1
        )
        updateSelectInputFn(
            session,
            "contrast_group2",
            choices = c("", currentGroups),
            selected = preservedSelections$group2
        )
    })
}

registerProtDesignInputSyncObserverShells <- function(
    input,
    session,
    designMatrix,
    removedSamples,
    factors,
    groups,
    registerSampleSelectionSyncObserver = registerProtDesignSampleSelectionSyncObserver,
    registerFactorGroupSyncObserver = registerProtDesignFactorGroupSyncObserver
) {
    registerSampleSelectionSyncObserver(
        input = input,
        session = session,
        designMatrix = designMatrix,
        removedSamples = removedSamples
    )

    registerFactorGroupSyncObserver(
        input = input,
        session = session,
        factors = factors,
        groups = groups
    )

    invisible(NULL)
}

registerProtDesignReplicateInputsOutput <- function(
    output,
    input,
    nsFn,
    renderUIFn = shiny::renderUI,
    reqFn = shiny::req,
    buildReplicateInputsFn = buildProtDesignReplicateInputs
) {
    output$replicate_inputs <- renderUIFn({
        reqFn(input$selected_runs)

        buildReplicateInputsFn(
            selectedRuns = input$selected_runs,
            nsFn = nsFn
        )
    })
}

registerProtDesignRangePreviewOutput <- function(
    output,
    input,
    renderTextFn = shiny::renderText,
    reqFn = shiny::req,
    formatRangePreviewFn = formatProtDesignRangePreview
) {
    output$range_preview <- renderTextFn({
        reqFn(input$samples_to_transform)
        reqFn(input$range_start, input$range_end)

        formatRangePreviewFn(
            selectedSamples = input$samples_to_transform,
            rangeStart = input$range_start,
            rangeEnd = input$range_end
        )
    })
}

registerProtDesignAvailableFactorsDisplayOutput <- function(
    output,
    factors,
    renderUIFn = shiny::renderUI,
    buildAvailableFactorsDisplayFn = buildProtDesignAvailableFactorsDisplay
) {
    output$available_factors_display <- renderUIFn({
        buildAvailableFactorsDisplayFn(factors())
    })
}

registerProtDesignDefinedContrastsDisplayOutput <- function(
    output,
    contrasts,
    formulaString,
    renderUIFn = shiny::renderUI,
    buildDefinedContrastsDisplayFn = buildProtDesignDefinedContrastsDisplay
) {
    output$defined_contrasts_display <- renderUIFn({
        buildDefinedContrastsDisplayFn(
            contrastData = contrasts(),
            formulaString = formulaString()
        )
    })
}

buildProtDesignDefinedContrastsDisplay <- function(contrastData, formulaString) {
    if (is.null(contrastData) || nrow(contrastData) == 0) {
        return(shiny::p("No contrasts defined yet."))
    }

    formulaUsesGroupPrefix <- grepl("~ *0 *\\+ *group", formulaString)

    contrastInfo <- vapply(seq_len(nrow(contrastData)), function(i) {
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

        paste0(friendlyName, "=", contrastString)
    }, character(1))

    shiny::tagList(
        lapply(contrastInfo, function(fullFormat) {
            shiny::p(shiny::tags$code(fullFormat))
        })
    )
}

buildProtDesignAvailableFactorsDisplay <- function(currentFactors) {
    if (length(currentFactors) == 0) {
        return(shiny::p("No factors defined yet (use the 'Factors' tab)."))
    }

    shiny::p(paste(currentFactors, collapse = ", "))
}

formatProtDesignRemovedSamplesDisplay <- function(removedSamples) {
    if (length(removedSamples) == 0) {
        return("No samples have been removed.")
    }

    paste0(
        "Removed ", length(removedSamples), " sample(s):\n",
        paste(gtools::mixedsort(removedSamples), collapse = "\n")
    )
}

formatProtDesignContrastFactorsInfo <- function(formulaString) {
    formulaUsesGroupPrefix <- grepl("~ *0 *\\+ *group", formulaString)
    if (formulaUsesGroupPrefix) {
        "Note: Contrasts will use 'group' prefix (e.g., groupGA_Control-groupGA_Elevated)\nbased on current formula: ~ 0 + group"
    } else {
        "Note: Contrasts will use group names as-is (e.g., GA_Control-GA_Elevated)"
    }
}

