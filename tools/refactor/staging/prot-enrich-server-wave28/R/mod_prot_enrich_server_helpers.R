resolveProtEnrichCurrentS4Object <- function(workflowData, daResultsList) {
  currentS4 <- NULL
  sourceLabel <- NULL

  if (!is.null(workflowData$state_manager)) {
    currentState <- workflowData$state_manager$current_state
    currentS4 <- workflowData$state_manager$getState(currentState)
    if (!is.null(currentS4)) {
      sourceLabel <- "state_manager"
    }
  }

  if (is.null(currentS4) && !is.null(daResultsList) && length(daResultsList) > 0) {
    firstResult <- daResultsList[[1]]
    firstObject <- tryCatch(firstResult$theObject, error = function(e) NULL)
    if (!is.null(firstObject)) {
      currentS4 <- firstObject
      sourceLabel <- "da_results_first_result"
    }
  }

  if (is.null(currentS4) && !is.null(daResultsList)) {
    combinedObject <- tryCatch(daResultsList$theObject, error = function(e) NULL)
    if (!is.null(combinedObject)) {
      currentS4 <- combinedObject
      sourceLabel <- "da_results_combined"
    }
  }

  list(
    currentS4 = currentS4,
    source = sourceLabel
  )
}

buildProtEnrichContrastChoices <- function(daResultsList, contrastsTbl = NULL) {
  contrastNames <- names(daResultsList)
  if (is.null(contrastNames)) {
    contrastNames <- character()
  }

  if (!is.null(contrastsTbl) && "friendly_names" %in% names(contrastsTbl)) {
    friendlyNames <- contrastsTbl$friendly_names
    return(list(
      contrastsAvailable = friendlyNames,
      contrastChoices = setNames(friendlyNames, friendlyNames),
      source = "friendly_names"
    ))
  }

  list(
    contrastsAvailable = contrastNames,
    contrastChoices = setNames(contrastNames, contrastNames),
    source = "raw_names"
  )
}

resolveProtEnrichRawContrastName <- function(selectedContrast, contrastsTbl = NULL) {
  rawContrastName <- selectedContrast
  sourceLabel <- "selected_contrast"

  if (!is.null(contrastsTbl) &&
      all(c("friendly_names", "contrasts") %in% names(contrastsTbl))) {
    matchingIdx <- which(contrastsTbl$friendly_names == selectedContrast)
    if (length(matchingIdx) > 0) {
      rawContrastName <- contrastsTbl$contrasts[matchingIdx[1]]
      sourceLabel <- "friendly_name"
    }
  }

  list(
    rawContrastName = rawContrastName,
    source = sourceLabel
  )
}

createProtEnrichRawContrastNameReactive <- function(input,
                                                    globalEnv = .GlobalEnv,
                                                    reactiveFn = shiny::reactive,
                                                    reqFn = shiny::req,
                                                    existsFn = exists,
                                                    getFn = get,
                                                    resolveRawContrastNameFn = resolveProtEnrichRawContrastName) {
  reactiveFn({
    reqFn(input$selected_contrast)

    contrastsTbl <- if (existsFn("contrasts_tbl", envir = globalEnv)) {
      getFn("contrasts_tbl", envir = globalEnv)
    } else {
      NULL
    }

    resolveRawContrastNameFn(input$selected_contrast, contrastsTbl)$rawContrastName
  })
}

resolveProtEnrichSelectedContrastResults <- function(selectedContrast, allEnrichmentResults, contrastsTbl = NULL) {
  resolvedContrast <- resolveProtEnrichRawContrastName(selectedContrast, contrastsTbl)
  availableContrasts <- names(allEnrichmentResults)
  if (is.null(availableContrasts)) {
    availableContrasts <- character()
  }

  contrastResults <- NULL
  if (!is.null(allEnrichmentResults) &&
      resolvedContrast$rawContrastName %in% availableContrasts) {
    contrastResults <- allEnrichmentResults[[resolvedContrast$rawContrastName]]
  }

  list(
    rawContrastName = resolvedContrast$rawContrastName,
    source = resolvedContrast$source,
    found = !is.null(contrastResults),
    availableContrasts = availableContrasts,
    gprofilerResults = if (!is.null(contrastResults)) contrastResults$gprofiler_results else NULL,
    clusterprofilerResults = if (!is.null(contrastResults)) contrastResults$clusterprofiler_results else NULL,
    stringdbResults = if (!is.null(contrastResults)) contrastResults$stringdb_results else NULL
  )
}

resolveProtEnrichSelectedDaResults <- function(selectedContrast, daResultsData, contrastsTbl = NULL) {
  availableKeys <- names(daResultsData)
  if (is.null(availableKeys)) {
    availableKeys <- character()
  }

  resolvedContrast <- resolveProtEnrichRawContrastName(selectedContrast, contrastsTbl)
  rawContrastName <- if (identical(resolvedContrast$source, "friendly_name")) {
    resolvedContrast$rawContrastName
  } else {
    NULL
  }
  selectedDaResults <- NULL
  sourceLabel <- NULL

  if (!is.null(rawContrastName)) {
    selectedDaResults <- daResultsData[[rawContrastName]]
    if (!is.null(selectedDaResults)) {
      sourceLabel <- "friendly_name"
    }
  }

  if (is.null(selectedDaResults)) {
    contrastParts <- stringr::str_split(selectedContrast, "_vs_")[[1]]
    if (length(contrastParts) == 2) {
      part1 <- contrastParts[1]
      part2 <- contrastParts[2]

      for (key in availableKeys) {
        if (stringr::str_detect(key, part1) && stringr::str_detect(key, part2)) {
          selectedDaResults <- daResultsData[[key]]
          rawContrastName <- key
          sourceLabel <- "fuzzy_match"
          break
        }
      }
    }
  }

  if (is.null(selectedDaResults)) {
    selectedDaResults <- daResultsData[[selectedContrast]]
    if (!is.null(selectedDaResults)) {
      rawContrastName <- selectedContrast
      sourceLabel <- "direct_key"
    }
  }

  list(
    selectedDaResults = selectedDaResults,
    rawContrastName = rawContrastName,
    source = sourceLabel,
    availableKeys = availableKeys,
    mappedRawContrastName = if (identical(resolvedContrast$source, "friendly_name")) {
      resolvedContrast$rawContrastName
    } else {
      NULL
    }
  )
}

