registerProtEnrichRunObserver <- function(input,
                                          enrichmentData,
                                          workflowData,
                                          experimentPaths,
                                          currentAnalysisMethodFn,
                                          observeEventFn = shiny::observeEvent,
                                          catFn = cat,
                                          runObserverPreflightFn = runProtEnrichObserverPreflight) {
  observeEventFn(input$run_enrichment_analysis, {
    catFn("=== STARTING ENRICHMENT ANALYSIS ===\n")

    runObserverPreflightFn(
      input = input,
      enrichmentData = enrichmentData,
      workflowData = workflowData,
      experimentPaths = experimentPaths,
      currentAnalysisMethodFn = currentAnalysisMethodFn
    )
  })
}

registerProtEnrichSelectedTabObserver <- function(selectedTabFn,
                                                  workflowData,
                                                  enrichmentData,
                                                  session,
                                                  observeEventFn = shiny::observeEvent,
                                                  resolveCurrentS4ObjectFn = resolveProtEnrichCurrentS4Object,
                                                  buildContrastChoicesFn = buildProtEnrichContrastChoices,
                                                  existsFn = exists,
                                                  getFn = get,
                                                  globalEnv = .GlobalEnv,
                                                  updateSelectInputFn = shiny::updateSelectInput,
                                                  showNotificationFn = shiny::showNotification,
                                                  catFn = cat) {
  observeEventFn(selectedTabFn(), {
    selectedTab <- selectedTabFn()
    observerState <- list(
      selectedTab = selectedTab,
      initialized = FALSE,
      reason = NULL,
      currentState = NULL,
      source = NULL,
      contrastChoices = NULL,
      errorMessage = NULL
    )

    catFn("--- Entering enrichment tab selection observer ---\n")
    catFn(sprintf("   enrichment_tab_observer Step: Selected tab = %s\n", selectedTab))

    if (!is.null(selectedTab) && selectedTab == "enrichment") {
      catFn("=== ENRICHMENT TAB CLICKED ===\n")
      catFn(sprintf(
        "   ENRICHMENT TAB Step: workflow_data$state_manager is NULL = %s\n",
        is.null(workflowData$state_manager)
      ))

      if (!is.null(workflowData$state_manager)) {
        currentState <- workflowData$state_manager$current_state
        validStatesForEnrichment <- c(
          "correlation_filtered",
          "normalized",
          "ruv_corrected",
          "protein_replicate_filtered"
        )
        observerState$currentState <- currentState

        catFn(sprintf("   ENRICHMENT TAB Step: Current state = '%s'\n", currentState))
        catFn(sprintf(
          "   ENRICHMENT TAB Step: Valid states for enrichment = %s\n",
          paste(validStatesForEnrichment, collapse = ", ")
        ))
        catFn(sprintf(
          "   ENRICHMENT TAB Step: DE results available = %s\n",
          !is.null(workflowData$da_analysis_results_list)
        ))

        if (currentState %in% validStatesForEnrichment &&
            !is.null(workflowData$da_analysis_results_list)) {
          catFn("*** AUTO-TRIGGERING ENRICHMENT INITIALIZATION (DE results found) ***\n")

          tryCatch({
            catFn("   ENRICHMENT TAB Step: Getting S4 object and DE results from workflow_data...\n")

            daResultsList <- workflowData$da_analysis_results_list
            resolvedContext <- resolveCurrentS4ObjectFn(workflowData, daResultsList)
            currentS4 <- resolvedContext$currentS4
            observerState$source <- resolvedContext$source

            if (!is.null(currentS4)) {
              catFn(sprintf(
                "   ENRICHMENT TAB Step: Got S4 from %s (class: %s)\n",
                resolvedContext$source,
                class(currentS4)
              ))
            }

            if (!is.null(currentS4) && !is.null(daResultsList)) {
              catFn(sprintf(
                "   ENRICHMENT TAB Step: S4 object retrieved, class = %s\n",
                class(currentS4)
              ))
              catFn(sprintf(
                "   ENRICHMENT TAB Step: DE results available for %d contrasts\n",
                length(daResultsList)
              ))

              enrichmentData$current_s4_object <- currentS4
              enrichmentData$da_results_data <- daResultsList

              contrastNames <- names(daResultsList)
              catFn(sprintf(
                "   ENRICHMENT TAB Step: Available contrast names: %s\n",
                paste(contrastNames, collapse = ", ")
              ))

              contrastsTbl <- if (existsFn("contrasts_tbl", envir = globalEnv)) {
                getFn("contrasts_tbl", envir = globalEnv)
              } else {
                NULL
              }

              contrastConfig <- buildContrastChoicesFn(daResultsList, contrastsTbl)
              enrichmentData$contrasts_available <- contrastConfig$contrastsAvailable
              observerState$contrastChoices <- contrastConfig$contrastChoices

              if (identical(contrastConfig$source, "friendly_names")) {
                catFn(sprintf(
                  "   ENRICHMENT TAB Step: Using friendly names: %s\n",
                  paste(contrastConfig$contrastsAvailable, collapse = ", ")
                ))
              }

              updateSelectInputFn(
                session,
                "selected_contrast",
                choices = contrastConfig$contrastChoices
              )

              observerState$initialized <- TRUE
              observerState$reason <- "initialized"
            } else {
              catFn("   ENRICHMENT TAB Step: S4 object or DE results are NULL\n")
              observerState$reason <- "missing_current_s4_or_results"
            }

            catFn("*** ENRICHMENT INITIALIZATION COMPLETED SUCCESSFULLY ***\n")
          }, error = function(e) {
            observerState$reason <<- "initialization_error"
            observerState$errorMessage <<- e$message
            catFn(paste("*** ERROR in enrichment initialization:", e$message, "\n"))
            showNotificationFn(
              paste("Error initializing enrichment analysis:", e$message),
              type = "error",
              duration = 10
            )
          })
        } else if (is.null(workflowData$da_analysis_results_list)) {
          observerState$reason <- "missing_da_results"
          catFn("*** No DE analysis results found. User needs to complete differential expression analysis first. ***\n")
          showNotificationFn(
            "Please complete the differential expression analysis before accessing enrichment analysis.",
            type = "warning",
            duration = 5
          )
        } else {
          observerState$reason <- "invalid_state"
          catFn(sprintf(
            "*** State '%s' is not valid for enrichment analysis. User needs to complete DE analysis. ***\n",
            currentState
          ))
          showNotificationFn(
            "Please complete the differential expression analysis before accessing enrichment analysis.",
            type = "warning",
            duration = 5
          )
        }
      } else {
        observerState$reason <- "missing_state_manager"
        catFn("*** workflow_data$state_manager is NULL - cannot check state ***\n")
      }
    } else {
      observerState$reason <- "tab_not_enrichment"
      catFn(sprintf(
        "   enrichment_tab_observer Step: Tab '%s' is not enrichment tab, ignoring\n",
        selectedTab
      ))
    }

    catFn("--- Exiting enrichment tab selection observer ---\n")
    observerState
  }, ignoreInit = TRUE)
}

registerProtEnrichDaResultsObserver <- function(workflowData,
                                                enrichmentData,
                                                session,
                                                observeEventFn = shiny::observeEvent,
                                                resolveCurrentS4ObjectFn = resolveProtEnrichCurrentS4Object,
                                                buildContrastChoicesFn = buildProtEnrichContrastChoices,
                                                existsFn = exists,
                                                getFn = get,
                                                globalEnv = .GlobalEnv,
                                                updateSelectInputFn = shiny::updateSelectInput,
                                                catFn = cat) {
  observeEventFn(workflowData$da_analysis_results_list, {
    catFn("*** ENRICHMENT: DE results detected - updating contrasts ***\n")

    daResultsList <- workflowData$da_analysis_results_list
    observerState <- list(
      hasResults = FALSE,
      currentS4Stored = FALSE,
      source = NULL,
      contrastChoices = NULL,
      reason = NULL
    )

    if (!is.null(daResultsList) && length(daResultsList) > 0) {
      observerState$hasResults <- TRUE
      resolvedContext <- resolveCurrentS4ObjectFn(workflowData, daResultsList)
      currentS4 <- resolvedContext$currentS4
      observerState$source <- resolvedContext$source

      if (!is.null(currentS4)) {
        enrichmentData$current_s4_object <- currentS4
        observerState$currentS4Stored <- TRUE
        catFn(sprintf(
          "   ENRICHMENT DE Observer: Got S4 from %s (class: %s)\n",
          resolvedContext$source,
          class(currentS4)
        ))
        catFn("   ENRICHMENT DE Observer: S4 object stored successfully\n")
      } else {
        catFn("   ENRICHMENT DE Observer: WARNING - Could not retrieve S4 object from any source\n")
      }

      enrichmentData$da_results_data <- daResultsList

      contrastNames <- names(daResultsList)
      catFn(sprintf(
        "   ENRICHMENT DE Observer: Available contrast names: %s\n",
        paste(contrastNames, collapse = ", ")
      ))

      contrastsTbl <- if (existsFn("contrasts_tbl", envir = globalEnv)) {
        getFn("contrasts_tbl", envir = globalEnv)
      } else {
        NULL
      }

      contrastConfig <- buildContrastChoicesFn(daResultsList, contrastsTbl)
      enrichmentData$contrasts_available <- contrastConfig$contrastsAvailable
      observerState$contrastChoices <- contrastConfig$contrastChoices

      if (identical(contrastConfig$source, "friendly_names")) {
        catFn(sprintf(
          "   ENRICHMENT DE Observer: Using friendly names: %s\n",
          paste(contrastConfig$contrastsAvailable, collapse = ", ")
        ))
      }

      updateSelectInputFn(
        session,
        "selected_contrast",
        choices = contrastConfig$contrastChoices
      )

      observerState$reason <- "updated"
      catFn(sprintf(
        "*** ENRICHMENT: Updated contrasts dropdown with %d contrasts ***\n",
        length(contrastConfig$contrastChoices)
      ))
    } else {
      observerState$reason <- "missing_or_empty_da_results"
    }

    observerState
  }, ignoreInit = TRUE, ignoreNULL = TRUE)
}

registerProtEnrichSelectedContrastObserver <- function(input,
                                                       enrichmentData,
                                                       observeFn = shiny::observe,
                                                       reqFn = shiny::req,
                                                       resolveSelectedContrastResultsFn = resolveProtEnrichSelectedContrastResults,
                                                       existsFn = exists,
                                                       getFn = get,
                                                       globalEnv = .GlobalEnv,
                                                       rowCountFn = nrow,
                                                       catFn = cat) {
  selectedContrastObserver <- observeFn({
    reqFn(input$selected_contrast)
    reqFn(enrichmentData$analysis_complete)
    reqFn(enrichmentData$all_enrichment_results)

    catFn(sprintf(
      "*** CONTRAST CHANGED: User selected '%s' ***\n",
      input$selected_contrast
    ))

    contrastsTbl <- if (existsFn("contrasts_tbl", envir = globalEnv)) {
      getFn("contrasts_tbl", envir = globalEnv)
    } else {
      NULL
    }

    contrastState <- resolveSelectedContrastResultsFn(
      input$selected_contrast,
      enrichmentData$all_enrichment_results,
      contrastsTbl
    )

    if (identical(contrastState$source, "friendly_name")) {
      catFn(sprintf(
        "*** CONTRAST MAPPING: '%s' -> '%s' ***\n",
        input$selected_contrast,
        contrastState$rawContrastName
      ))
    }

    if (contrastState$found) {
      catFn(sprintf(
        "*** UPDATING RESULTS: Found results for contrast '%s' ***\n",
        contrastState$rawContrastName
      ))

      enrichmentData$gprofiler_results <- contrastState$gprofilerResults
      enrichmentData$clusterprofiler_results <- contrastState$clusterprofilerResults
      enrichmentData$stringdb_results <- contrastState$stringdbResults

      if (!is.null(contrastState$gprofilerResults)) {
        catFn(sprintf(
          "*** UPDATED: %d gprofiler2 results ***\n",
          rowCountFn(contrastState$gprofilerResults)
        ))
      }
      if (!is.null(contrastState$clusterprofilerResults)) {
        catFn(sprintf(
          "*** UPDATED: %d clusterProfileR results ***\n",
          rowCountFn(contrastState$clusterprofilerResults)
        ))
      }
      if (!is.null(contrastState$stringdbResults)) {
        catFn(sprintf(
          "*** UPDATED: %d stringDB results ***\n",
          rowCountFn(contrastState$stringdbResults)
        ))
      }
    } else {
      catFn(sprintf(
        "*** WARNING: No results found for contrast '%s' ***\n",
        contrastState$rawContrastName
      ))
      catFn(sprintf(
        "*** AVAILABLE CONTRASTS: %s ***\n",
        paste(contrastState$availableContrasts, collapse = ", ")
      ))

      enrichmentData$gprofiler_results <- NULL
      enrichmentData$clusterprofiler_results <- NULL
      enrichmentData$stringdb_results <- NULL
    }

    contrastState
  })

  selectedContrastObserver
}

