registerProtEnrichPlotOutputs <- function(output,
                                          input,
                                          enrichmentData,
                                          rawContrastNameFn,
                                          renderGprofilerPlotFn = renderProtEnrichGprofilerPlot,
                                          renderClusterProfilerPlotFn = renderProtEnrichClusterProfilerPlot) {
  gprofilerPlot <- renderGprofilerPlotFn(
    analysisComplete = enrichmentData$analysis_complete,
    enrichmentResultsFull = enrichmentData$enrichment_results_full,
    rawContrast = rawContrastNameFn(),
    directionFilter = input$gprofiler_direction_filter
  )

  clusterprofilerPlot <- renderClusterProfilerPlotFn(
    analysisComplete = enrichmentData$analysis_complete,
    enrichmentResultsFull = enrichmentData$enrichment_results_full,
    rawContrast = rawContrastNameFn(),
    directionFilter = input$clusterprofiler_direction_filter
  )

  output$gprofiler_plot <- gprofilerPlot
  output$clusterprofiler_plot <- clusterprofilerPlot

  list(
    gprofilerPlot = gprofilerPlot,
    clusterprofilerPlot = clusterprofilerPlot
  )
}

registerProtEnrichStringDbPlotOutput <- function(output,
                                                 renderStringDbPlotFn = renderProtEnrichStringDbPlot) {
  stringdbPlot <- renderStringDbPlotFn()
  output$stringdb_plot <- stringdbPlot
  stringdbPlot
}

registerProtEnrichResultsDownloadHandler <- function(output,
                                                     input,
                                                     enrichmentData,
                                                     currentAnalysisMethodFn,
                                                     downloadHandlerFn = shiny::downloadHandler,
                                                     buildFilenameFn = buildProtEnrichResultsDownloadFilename,
                                                     writeArchiveFn = writeProtEnrichResultsDownloadArchive) {
  downloadHandler <- downloadHandlerFn(
    filename = function() {
      buildFilenameFn(input$selected_contrast)
    },
    content = function(file) {
      writeArchiveFn(
        file = file,
        selectedContrast = input$selected_contrast,
        methodInfo = currentAnalysisMethodFn(),
        organismTaxid = input$organism_taxid,
        upCutoff = input$up_cutoff,
        downCutoff = input$down_cutoff,
        qCutoff = input$q_cutoff,
        gprofilerResults = enrichmentData$gprofiler_results,
        clusterprofilerResults = enrichmentData$clusterprofiler_results,
        stringdbResults = enrichmentData$stringdb_results
      )
    }
  )

  output$download_enrichment_results <- downloadHandler
  downloadHandler
}

registerProtEnrichAnalysisMethodDisplayOutput <- function(output,
                                                          currentAnalysisMethodFn,
                                                          renderTextFn = shiny::renderText,
                                                          formatAnalysisMethodTextFn = formatProtEnrichAnalysisMethodText) {
  analysisMethodDisplay <- renderTextFn({
    formatAnalysisMethodTextFn(currentAnalysisMethodFn())
  })

  output$analysis_method_display <- analysisMethodDisplay
  analysisMethodDisplay
}

registerProtEnrichTaxonIdObserver <- function(workflowData,
                                              session,
                                              observeEventFn = shiny::observeEvent,
                                              updateTextInputFn = shiny::updateTextInput) {
  observeEventFn(
    eventExpr = workflowData$taxon_id,
    handlerExpr = {
    taxonId <- workflowData$taxon_id
    observerState <- list(
      taxonId = taxonId,
      updated = FALSE,
      reason = NULL
    )

    if (!is.null(taxonId)) {
      updateTextInputFn(session, "organism_taxid", value = taxonId)
      observerState$updated <- TRUE
      observerState$reason <- "updated"
    } else {
      observerState$reason <- "missing_taxon_id"
    }

    observerState
  },
    ignoreInit = TRUE,
    ignoreNULL = TRUE
  )
}

registerProtEnrichMixedSpeciesObserver <- function(workflowData,
                                                   session,
                                                   observeEventFn = shiny::observeEvent,
                                                   updateCheckboxInputFn = shiny::updateCheckboxInput,
                                                   showNotificationFn = shiny::showNotification,
                                                   catFn = cat) {
  observeEventFn(
    eventExpr = workflowData$mixed_species_analysis,
    handlerExpr = {
    mixedSpeciesAnalysis <- workflowData$mixed_species_analysis
    observerState <- list(
      mixedSpeciesAnalysis = mixedSpeciesAnalysis,
      filterEnabled = FALSE,
      notificationShown = FALSE,
      reason = NULL
    )

    if (!is.null(mixedSpeciesAnalysis) &&
        isTRUE(mixedSpeciesAnalysis$enabled)) {
      updateCheckboxInputFn(session, "enable_organism_filter", value = TRUE)
      catFn("*** ENRICHMENT: Auto-enabled organism filter (mixed species FASTA detected at import) ***\n")

      showNotificationFn(
        sprintf(
          "Multi-species FASTA detected. Filtering to %s enabled.",
          mixedSpeciesAnalysis$selected_organism
        ),
        type = "message",
        duration = 5
      )

      observerState$filterEnabled <- TRUE
      observerState$notificationShown <- TRUE
      observerState$reason <- "enabled_filter"
    } else {
      observerState$reason <- "mixed_species_disabled"
    }

    observerState
  },
    ignoreInit = TRUE,
    ignoreNULL = TRUE
  )
}

registerProtEnrichContrastsDisplayOutput <- function(output,
                                                     enrichmentData,
                                                     renderTextFn = shiny::renderText,
                                                     formatContrastsTextFn = formatProtEnrichContrastsText) {
  contrastsDisplay <- renderTextFn({
    formatContrastsTextFn(enrichmentData$contrasts_available)
  })

  output$contrasts_display <- contrastsDisplay
  contrastsDisplay
}

registerProtEnrichGprofilerResultsTableOutput <- function(output,
                                                          input,
                                                          enrichmentData,
                                                          renderResultsTableFn = renderProtEnrichGprofilerResultsTable) {
  gprofilerResultsTable <- renderResultsTableFn(
    gprofilerResults = enrichmentData$gprofiler_results,
    directionFilter = input$gprofiler_direction_filter
  )

  output$gprofiler_results_table <- gprofilerResultsTable
  gprofilerResultsTable
}

registerProtEnrichGprofilerSummaryOutput <- function(output,
                                                     input,
                                                     enrichmentData,
                                                     renderTextFn = shiny::renderText,
                                                     formatSummaryTextFn = formatProtEnrichGprofilerSummaryText) {
  gprofilerSummaryStats <- renderTextFn({
    formatSummaryTextFn(
      gprofilerResults = enrichmentData$gprofiler_results,
      directionFilter = input$gprofiler_direction_filter
    )
  })

  output$gprofiler_summary_stats <- gprofilerSummaryStats
  gprofilerSummaryStats
}

registerProtEnrichClusterProfilerResultsTableOutput <- function(output,
                                                                input,
                                                                enrichmentData,
                                                                renderResultsTableFn = renderProtEnrichClusterProfilerResultsTable) {
  clusterprofilerResultsTable <- renderResultsTableFn(
    clusterprofilerResults = enrichmentData$clusterprofiler_results,
    directionFilter = input$clusterprofiler_direction_filter
  )

  output$clusterprofiler_results_table <- clusterprofilerResultsTable
  clusterprofilerResultsTable
}

registerProtEnrichClusterProfilerSummaryOutput <- function(output,
                                                           input,
                                                           enrichmentData,
                                                           renderTextFn = shiny::renderText,
                                                           formatSummaryTextFn = formatProtEnrichClusterProfilerSummaryText) {
  clusterprofilerSummaryStats <- renderTextFn({
    formatSummaryTextFn(
      clusterprofilerResults = enrichmentData$clusterprofiler_results,
      directionFilter = input$clusterprofiler_direction_filter
    )
  })

  output$clusterprofiler_summary_stats <- clusterprofilerSummaryStats
  clusterprofilerSummaryStats
}

registerProtEnrichStringDbResultsTableOutput <- function(output,
                                                         input,
                                                         enrichmentData,
                                                         renderResultsTableFn = renderProtEnrichStringDbResultsTable) {
  stringdbResultsTable <- renderResultsTableFn(
    stringdbResults = enrichmentData$stringdb_results,
    filterSignificant = input$stringdb_filter_significant,
    enrichmentPValThresh = input$enrichment_p_val_thresh,
    maxResults = input$stringdb_max_results
  )

  output$stringdb_results_table <- stringdbResultsTable
  stringdbResultsTable
}

registerProtEnrichStringDbSummaryOutput <- function(output,
                                                    enrichmentData,
                                                    renderTextFn = shiny::renderText,
                                                    formatSummaryTextFn = formatProtEnrichStringDbSummaryText) {
  stringdbSummaryStats <- renderTextFn({
    formatSummaryTextFn(
      stringdbResults = enrichmentData$stringdb_results
    )
  })

  output$stringdb_summary_stats <- stringdbSummaryStats
  stringdbSummaryStats
}

registerProtEnrichStatusOutput <- function(output,
                                           input,
                                           enrichmentData,
                                           currentAnalysisMethodFn,
                                           renderTextFn = shiny::renderText,
                                           formatStatusTextFn = formatProtEnrichStatusText) {
  enrichmentStatus <- renderTextFn({
    methodInfo <- NULL
    if (enrichmentData$analysis_complete) {
      methodInfo <- currentAnalysisMethodFn()
    }

    formatStatusTextFn(
      analysisComplete = enrichmentData$analysis_complete,
      methodInfo = methodInfo,
      selectedContrast = input$selected_contrast,
      upCutoff = input$up_cutoff,
      downCutoff = input$down_cutoff,
      qCutoff = input$q_cutoff,
      gprofilerResults = enrichmentData$gprofiler_results,
      clusterprofilerResults = enrichmentData$clusterprofiler_results,
      stringdbResults = enrichmentData$stringdb_results
    )
  })

  output$enrichment_status <- enrichmentStatus
  enrichmentStatus
}

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
                                                  initialiseArtifactDataFn = initialiseProtDiaEnrichArtifactData,
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
        currentState <- workflowStateCurrentName(workflowData$state_manager)
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
            if (protDiaEnrichArtifactEligible(workflowData, "queries") &&
                isProtDiaDaArtifactIndex(daResultsList)) {
              artifactState <- initialiseArtifactDataFn(
                workflowData,
                enrichmentData
              )
              observerState$source <- "artifact_index"
              observerState$contrastChoices <-
                artifactState$contrastConfig$contrastChoices
              updateSelectInputFn(
                session,
                "selected_contrast",
                choices = observerState$contrastChoices
              )
              observerState$initialized <- TRUE
              observerState$reason <- "artifact_initialized"
              return(observerState)
            }
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
                                                initialiseArtifactDataFn = initialiseProtDiaEnrichArtifactData,
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
      if (protDiaEnrichArtifactEligible(workflowData, "queries") &&
          isProtDiaDaArtifactIndex(daResultsList)) {
        artifactState <- initialiseArtifactDataFn(
          workflowData,
          enrichmentData
        )
        observerState$currentS4Stored <- !is.null(artifactState$currentS4)
        observerState$source <- "artifact_index"
        observerState$contrastChoices <-
          artifactState$contrastConfig$contrastChoices
        updateSelectInputFn(
          session,
          "selected_contrast",
          choices = observerState$contrastChoices
        )
        observerState$reason <- "artifact_updated"
        return(observerState)
      }
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
