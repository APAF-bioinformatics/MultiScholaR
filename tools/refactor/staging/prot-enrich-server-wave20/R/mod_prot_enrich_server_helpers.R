updateProtEnrichStateManagerUiParams <- function(workflowData,
                                                 storedUiParams,
                                                 currentDataStates = c(
                                                   "correlation_filtered",
                                                   "ruv_corrected",
                                                   "protein_replicate_filtered"
                                                 ),
                                                 detectFn = purrr::detect,
                                                 catFn = cat) {
  if (!isTRUE(storedUiParams)) {
    return(list(
      attempted = FALSE,
      currentDataState = NULL,
      availableStates = NULL,
      warning = NULL,
      updated = FALSE
    ))
  }

  catFn("   ENRICHMENT Step: Updating R6 state with enrichment UI parameters\n")

  currentDataState <- NULL
  availableStates <- NULL
  warningMessage <- NULL

  tryCatch({
    availableStates <- workflowData$state_manager$getHistory()
    currentDataState <- detectFn(currentDataStates, ~ .x %in% availableStates)

    if (!is.null(currentDataState)) {
      catFn("   ENRICHMENT Step: Skipping state update (updateState method verification needed)\n")
    }
  }, error = function(e) {
    warningMessage <<- e$message
    catFn(sprintf(
      "   ENRICHMENT Step: Warning - could not update state with UI parameters: %s\n",
      e$message
    ))
  })

  list(
    attempted = TRUE,
    currentDataState = currentDataState,
    availableStates = availableStates,
    warning = warningMessage,
    updated = FALSE
  )
}

saveProtEnrichCompletedState <- function(workflowData,
                                         enrichmentResults,
                                         selectedContrast,
                                         methodInfo,
                                         upCutoff,
                                         downCutoff,
                                         qCutoff,
                                         organismTaxid,
                                         pathwayDir,
                                         catFn = cat) {
  catFn("   ENRICHMENT Step: Saving results to R6 state manager...\n")

  warningMessage <- NULL
  saved <- FALSE

  tryCatch({
    workflowData$state_manager$saveState(
      state_name = "enrichment_completed",
      s4_data_object = enrichmentResults,
      config_object = list(
        selected_contrast = selectedContrast,
        analysis_method = methodInfo$method,
        organism_supported = methodInfo$supported,
        up_cutoff = upCutoff,
        down_cutoff = downCutoff,
        q_cutoff = qCutoff,
        organism_taxid = organismTaxid,
        pathway_dir = pathwayDir
      ),
      description = paste(
        "Enrichment analysis completed using",
        methodInfo$method,
        "for contrast:",
        selectedContrast
      )
    )

    saved <- TRUE
    catFn("   ENRICHMENT Step: Successfully saved state to R6 state manager\n")
  }, error = function(e) {
    warningMessage <<- e$message
    catFn(sprintf(
      "   ENRICHMENT Step: Warning saving to R6 state manager: %s\n",
      e$message
    ))
  })

  list(
    attempted = TRUE,
    saved = saved,
    warning = warningMessage
  )
}

completeProtEnrichTabStatus <- function(workflowData,
                                        tabName = "enrichment_analysis",
                                        status = "complete") {
  updatedStatus <- workflowData$tab_status
  updatedStatus[[tabName]] <- status
  workflowData$tab_status <- updatedStatus

  list(
    tabName = tabName,
    status = status,
    tabStatus = updatedStatus
  )
}

