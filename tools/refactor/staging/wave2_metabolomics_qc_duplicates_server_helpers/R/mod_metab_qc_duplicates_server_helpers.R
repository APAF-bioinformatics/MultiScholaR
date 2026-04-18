resolveMetabDuplicateAssayData <- function(
    assayList
    , metaboliteIdCol
    , assayNames = names(assayList)
    , resolveDuplicateFeaturesByIntensityFn = resolveDuplicateFeaturesByIntensity
    , logWarnFn = logger::log_warn
) {
    statsList <- list()

    resolvedAssayList <- lapply(seq_along(assayList), function(i) {
        assayData <- assayList[[i]]
        assayName <- if (!is.null(assayNames)) assayNames[i] else paste0("Assay_", i)

        originalRows <- nrow(assayData)

        sampleCols <- names(assayData)[sapply(assayData, is.numeric)]

        if (length(sampleCols) == 0) {
            logWarnFn(paste("No numeric columns found in assay:", assayName))
            statsList[[assayName]] <<- list(
                original = originalRows
                , resolved = originalRows
                , removed = 0
            )
            return(assayData)
        }

        resolvedAssay <- resolveDuplicateFeaturesByIntensityFn(
            assay_tibble = assayData
            , id_col = metaboliteIdCol
            , sample_cols = sampleCols
        )

        resolvedRows <- nrow(resolvedAssay)

        statsList[[assayName]] <<- list(
            original = originalRows
            , resolved = resolvedRows
            , removed = originalRows - resolvedRows
        )

        resolvedAssay
    })

    names(resolvedAssayList) <- assayNames

    list(
        resolvedAssayList = resolvedAssayList
        , statsList = statsList
    )
}

buildMetabDuplicateResolutionSummary <- function(
    statsList
    , stateName = "metab_duplicates_resolved"
) {
    perAssayText <- if (length(statsList) == 0) {
        "  No assay statistics available"
    } else {
        paste(vapply(names(statsList), function(assayName) {
            stats <- statsList[[assayName]]
            sprintf(
                "  %s: %d -> %d rows (removed %d duplicates)"
                , assayName
                , stats$original
                , stats$resolved
                , stats$removed
            )
        }, character(1)), collapse = "\n")
    }

    totalRemoved <- if (length(statsList) == 0) {
        0
    } else {
        sum(vapply(statsList, `[[`, numeric(1), "removed"))
    }

    list(
        totalRemoved = totalRemoved
        , resultText = paste(
            "Duplicate Resolution Complete"
            , "=============================="
            , "Strategy: Keep feature with highest average intensity"
            , ""
            , "Per-Assay Results:"
            , perAssayText
            , ""
            , sprintf("Total duplicate rows removed: %d", totalRemoved)
            , ""
            , sprintf("State saved as: '%s'", stateName)
            , sep = "\n"
        )
    )
}

detectMetabDuplicateFeatures <- function(
    stateManager
    , duplicateFinderFn = findMetabDuplicateFeatureIDs
    , reqFn = shiny::req
    , inheritsFn = inherits
    , expectedClass = "MetaboliteAssayData"
) {
    reqFn(stateManager)

    currentS4 <- stateManager$getState()
    reqFn(currentS4)

    if (!inheritsFn(currentS4, expectedClass)) {
        stop(sprintf("Current state is not a %s object", expectedClass))
    }

    duplicatesList <- duplicateFinderFn(currentS4)
    totalDuplicates <- sum(vapply(duplicatesList, function(dupDf) {
        if (is.null(dupDf)) {
            return(0L)
        }

        nrow(dupDf)
    }, integer(1)))

    list(
        duplicatesList = duplicatesList
        , totalDuplicates = totalDuplicates
    )
}

revertMetabDuplicateResolution <- function(
    stateManager
    , reqFn = shiny::req
    , historyGetterFn = function(manager) manager$getHistory()
    , revertStateFn = function(manager, stateName) manager$revertToState(stateName)
) {
    reqFn(stateManager)

    history <- historyGetterFn(stateManager)
    if (length(history) <= 1) {
        stop("No previous state to revert to.")
    }

    previousStateName <- history[[length(history) - 1L]]
    revertStateFn(stateManager, previousStateName)

    list(
        previousStateName = previousStateName
        , resultText = paste("Reverted to previous state:", previousStateName)
    )
}

