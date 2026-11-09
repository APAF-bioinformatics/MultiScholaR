resolveProtEnrichCurrentS4Object <- function(workflowData, daResultsList) {
  if (isDIAResultManifest(daResultsList)) return(list(
    currentS4 = daResultsList$source, source = "component_references"))
  currentS4 <- NULL
  sourceLabel <- NULL

  if (!is.null(workflowData$state_manager)) {
    currentState <- workflowData$state_manager$current_state
    currentS4 <- readDIAStateObject(workflowData$state_manager, currentState)
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

createProtEnrichRawContrastNameReactive <- function(input,
                                                    workflowData = NULL,
                                                    globalEnv = NULL,
                                                    reactiveFn = shiny::reactive,
                                                    reqFn = shiny::req,
                                                    existsFn = exists,
                                                    getFn = get,
                                                    resolveRawContrastNameFn = resolveProtEnrichRawContrastName) {
  reactiveFn({
    reqFn(input$selected_contrast)

    contrastsTbl <- resolveProtEnrichContrastsTable(workflowData, globalEnv, existsFn, getFn)

    resolveRawContrastNameFn(input$selected_contrast, contrastsTbl)$rawContrastName
  })
}
