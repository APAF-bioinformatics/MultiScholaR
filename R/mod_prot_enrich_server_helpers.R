setupProtEnrichReactiveValues <- function(createReactiveValuesFn = createProtEnrichReactiveValues,
                                          catFn = cat) {
  enrichmentData <- createReactiveValuesFn()
  catFn("   mod_prot_enrich_server Step: Reactive values initialized\n")

  list(
    enrichmentData = enrichmentData,
    reason = "created"
  )
}

setupProtEnrichSupportedOrganismsReactive <- function(createSupportedOrganismsReactiveFn = createProtEnrichSupportedOrganismsReactive) {
  supportedOrganisms <- createSupportedOrganismsReactiveFn()

  list(
    supportedOrganisms = supportedOrganisms,
    reason = "created"
  )
}

setupProtEnrichCurrentAnalysisMethodReactive <- function(input,
                                                         enrichmentData,
                                                         supportedOrganismsFn,
                                                         createCurrentAnalysisMethodReactiveFn = createProtEnrichCurrentAnalysisMethodReactive) {
  currentAnalysisMethod <- createCurrentAnalysisMethodReactiveFn(
    input = input,
    enrichmentData = enrichmentData,
    supportedOrganismsFn = supportedOrganismsFn
  )

  list(
    currentAnalysisMethod = currentAnalysisMethod,
    reason = "created"
  )
}

setupProtEnrichAnalysisMethodBootstrap <- function(input,
                                                   enrichmentData,
                                                   setupSupportedOrganismsReactiveFn = setupProtEnrichSupportedOrganismsReactive,
                                                   setupCurrentAnalysisMethodReactiveFn = setupProtEnrichCurrentAnalysisMethodReactive) {
  supportedOrganismsSetup <- setupSupportedOrganismsReactiveFn()
  currentAnalysisMethodSetup <- setupCurrentAnalysisMethodReactiveFn(
    input = input,
    enrichmentData = enrichmentData,
    supportedOrganismsFn = supportedOrganismsSetup$supportedOrganisms
  )

  list(
    supportedOrganisms = supportedOrganismsSetup$supportedOrganisms,
    currentAnalysisMethod = currentAnalysisMethodSetup$currentAnalysisMethod,
    reason = "created"
  )
}

setupProtEnrichRawContrastNameReactive <- function(input,
                                                   createRawContrastNameReactiveFn = createProtEnrichRawContrastNameReactive) {
  rawContrastName <- createRawContrastNameReactiveFn(input = input)

  list(
    rawContrastName = rawContrastName,
    reason = "created"
  )
}

setupProtEnrichDaResultsObserverRegistration <- function(workflowData,
                                                         enrichmentData,
                                                         session,
                                                         registerDaResultsObserverFn = registerProtEnrichDaResultsObserver) {
  registration <- registerDaResultsObserverFn(
    workflowData = workflowData,
    enrichmentData = enrichmentData,
    session = session
  )

  list(
    registration = registration,
    reason = "registered"
  )
}

setupProtEnrichSelectedContrastObserverRegistration <- function(input,
                                                                enrichmentData,
                                                                registerSelectedContrastObserverFn = registerProtEnrichSelectedContrastObserver) {
  registration <- registerSelectedContrastObserverFn(
    input = input,
    enrichmentData = enrichmentData
  )

  list(
    registration = registration,
    reason = "registered"
  )
}

setupProtEnrichContrastsDisplayOutputRegistration <- function(output,
                                                              enrichmentData,
                                                              registerContrastsDisplayOutputFn = registerProtEnrichContrastsDisplayOutput) {
  registration <- registerContrastsDisplayOutputFn(
    output = output,
    enrichmentData = enrichmentData
  )

  list(
    registration = registration,
    reason = "registered"
  )
}

setupProtEnrichStatusOutputRegistration <- function(output,
                                                    input,
                                                    enrichmentData,
                                                    currentAnalysisMethodFn,
                                                    registerStatusOutputFn = registerProtEnrichStatusOutput) {
  registration <- registerStatusOutputFn(
    output = output,
    input = input,
    enrichmentData = enrichmentData,
    currentAnalysisMethodFn = currentAnalysisMethodFn
  )

  list(
    registration = registration,
    reason = "registered"
  )
}

setupProtEnrichDisplayStatusOutputBootstrap <- function(output,
                                                        input,
                                                        enrichmentData,
                                                        currentAnalysisMethodFn,
                                                        setupAnalysisMethodDisplayOutputRegistrationFn = setupProtEnrichAnalysisMethodDisplayOutputRegistration,
                                                        setupContrastsDisplayOutputRegistrationFn = setupProtEnrichContrastsDisplayOutputRegistration,
                                                        setupStatusOutputRegistrationFn = setupProtEnrichStatusOutputRegistration) {
  analysisMethodDisplayOutputRegistration <- setupAnalysisMethodDisplayOutputRegistrationFn(
    output = output,
    currentAnalysisMethodFn = currentAnalysisMethodFn
  )
  contrastsDisplayOutputRegistration <- setupContrastsDisplayOutputRegistrationFn(
    output = output,
    enrichmentData = enrichmentData
  )
  statusOutputRegistration <- setupStatusOutputRegistrationFn(
    output = output,
    input = input,
    enrichmentData = enrichmentData,
    currentAnalysisMethodFn = currentAnalysisMethodFn
  )

  list(
    analysisMethodDisplayOutputRegistration = analysisMethodDisplayOutputRegistration$registration,
    contrastsDisplayOutputRegistration = contrastsDisplayOutputRegistration$registration,
    statusOutputRegistration = statusOutputRegistration$registration,
    reason = "registered"
  )
}

setupProtEnrichRunObserverRegistration <- function(input,
                                                   enrichmentData,
                                                   workflowData,
                                                   experimentPaths,
                                                   currentAnalysisMethodFn,
                                                   runObserverPreflightFn = runProtEnrichObserverPreflight,
                                                   registerRunObserverFn = registerProtEnrichRunObserver) {
  registration <- registerRunObserverFn(
    input = input,
    enrichmentData = enrichmentData,
    workflowData = workflowData,
    experimentPaths = experimentPaths,
    currentAnalysisMethodFn = currentAnalysisMethodFn,
    runObserverPreflightFn = runObserverPreflightFn
  )

  list(
    registration = registration,
    reason = "registered"
  )
}

setupProtEnrichGprofilerResultsTableOutputRegistration <- function(output,
                                                                   input,
                                                                   enrichmentData,
                                                                   registerGprofilerResultsTableOutputFn = registerProtEnrichGprofilerResultsTableOutput) {
  registration <- registerGprofilerResultsTableOutputFn(
    output = output,
    input = input,
    enrichmentData = enrichmentData
  )

  list(
    registration = registration,
    reason = "registered"
  )
}

setupProtEnrichGprofilerSummaryOutputRegistration <- function(output,
                                                              input,
                                                              enrichmentData,
                                                              registerGprofilerSummaryOutputFn = registerProtEnrichGprofilerSummaryOutput) {
  registration <- registerGprofilerSummaryOutputFn(
    output = output,
    input = input,
    enrichmentData = enrichmentData
  )

  list(
    registration = registration,
    reason = "registered"
  )
}

setupProtEnrichClusterProfilerResultsTableOutputRegistration <- function(output,
                                                                         input,
                                                                         enrichmentData,
                                                                         registerClusterProfilerResultsTableOutputFn = registerProtEnrichClusterProfilerResultsTableOutput) {
  registration <- registerClusterProfilerResultsTableOutputFn(
    output = output,
    input = input,
    enrichmentData = enrichmentData
  )

  list(
    registration = registration,
    reason = "registered"
  )
}

setupProtEnrichClusterProfilerSummaryOutputRegistration <- function(output,
                                                                    input,
                                                                    enrichmentData,
                                                                    registerClusterProfilerSummaryOutputFn = registerProtEnrichClusterProfilerSummaryOutput) {
  registration <- registerClusterProfilerSummaryOutputFn(
    output = output,
    input = input,
    enrichmentData = enrichmentData
  )

  list(
    registration = registration,
    reason = "registered"
  )
}

setupProtEnrichStringDbResultsTableOutputRegistration <- function(output,
                                                                  input,
                                                                  enrichmentData,
                                                                  registerStringDbResultsTableOutputFn = registerProtEnrichStringDbResultsTableOutput) {
  registration <- registerStringDbResultsTableOutputFn(
    output = output,
    input = input,
    enrichmentData = enrichmentData
  )

  list(
    registration = registration,
    reason = "registered"
  )
}

setupProtEnrichStringDbSummaryOutputRegistration <- function(output,
                                                             enrichmentData,
                                                             registerStringDbSummaryOutputFn = registerProtEnrichStringDbSummaryOutput) {
  registration <- registerStringDbSummaryOutputFn(
    output = output,
    enrichmentData = enrichmentData
  )

  list(
    registration = registration,
    reason = "registered"
  )
}

setupProtEnrichResultsSummaryOutputBootstrap <- function(output,
                                                         input,
                                                         enrichmentData,
                                                         setupGprofilerResultsTableOutputRegistrationFn = setupProtEnrichGprofilerResultsTableOutputRegistration,
                                                         setupGprofilerSummaryOutputRegistrationFn = setupProtEnrichGprofilerSummaryOutputRegistration,
                                                         setupClusterProfilerResultsTableOutputRegistrationFn = setupProtEnrichClusterProfilerResultsTableOutputRegistration,
                                                         setupClusterProfilerSummaryOutputRegistrationFn = setupProtEnrichClusterProfilerSummaryOutputRegistration,
                                                         setupStringDbResultsTableOutputRegistrationFn = setupProtEnrichStringDbResultsTableOutputRegistration,
                                                         setupStringDbSummaryOutputRegistrationFn = setupProtEnrichStringDbSummaryOutputRegistration) {
  gprofilerResultsTableOutputRegistration <- setupGprofilerResultsTableOutputRegistrationFn(
    output = output,
    input = input,
    enrichmentData = enrichmentData
  )
  gprofilerSummaryOutputRegistration <- setupGprofilerSummaryOutputRegistrationFn(
    output = output,
    input = input,
    enrichmentData = enrichmentData
  )
  clusterProfilerResultsTableOutputRegistration <- setupClusterProfilerResultsTableOutputRegistrationFn(
    output = output,
    input = input,
    enrichmentData = enrichmentData
  )
  clusterProfilerSummaryOutputRegistration <- setupClusterProfilerSummaryOutputRegistrationFn(
    output = output,
    input = input,
    enrichmentData = enrichmentData
  )
  stringDbResultsTableOutputRegistration <- setupStringDbResultsTableOutputRegistrationFn(
    output = output,
    input = input,
    enrichmentData = enrichmentData
  )
  stringDbSummaryOutputRegistration <- setupStringDbSummaryOutputRegistrationFn(
    output = output,
    enrichmentData = enrichmentData
  )

  list(
    gprofilerResultsTableOutputRegistration = gprofilerResultsTableOutputRegistration$registration,
    gprofilerSummaryOutputRegistration = gprofilerSummaryOutputRegistration$registration,
    clusterProfilerResultsTableOutputRegistration = clusterProfilerResultsTableOutputRegistration$registration,
    clusterProfilerSummaryOutputRegistration = clusterProfilerSummaryOutputRegistration$registration,
    stringDbResultsTableOutputRegistration = stringDbResultsTableOutputRegistration$registration,
    stringDbSummaryOutputRegistration = stringDbSummaryOutputRegistration$registration,
    reason = "registered"
  )
}

setupProtEnrichPlotOutputsRegistration <- function(output,
                                                   input,
                                                   enrichmentData,
                                                   rawContrastNameFn,
                                                   registerPlotOutputsFn = registerProtEnrichPlotOutputs) {
  registration <- registerPlotOutputsFn(
    output = output,
    input = input,
    enrichmentData = enrichmentData,
    rawContrastNameFn = rawContrastNameFn
  )

  list(
    registration = registration,
    reason = "registered"
  )
}

setupProtEnrichStringDbPlotOutputRegistration <- function(output,
                                                          registerStringDbPlotOutputFn = registerProtEnrichStringDbPlotOutput) {
  registration <- registerStringDbPlotOutputFn(
    output = output
  )

  list(
    registration = registration,
    reason = "registered"
  )
}

setupProtEnrichPlotOutputBootstrap <- function(output,
                                               input,
                                               enrichmentData,
                                               setupRawContrastNameReactiveFn = setupProtEnrichRawContrastNameReactive,
                                               setupPlotOutputsRegistrationFn = setupProtEnrichPlotOutputsRegistration,
                                               setupStringDbPlotOutputRegistrationFn = setupProtEnrichStringDbPlotOutputRegistration) {
  rawContrastNameSetup <- setupRawContrastNameReactiveFn(input = input)
  plotOutputsRegistration <- setupPlotOutputsRegistrationFn(
    output = output,
    input = input,
    enrichmentData = enrichmentData,
    rawContrastNameFn = rawContrastNameSetup$rawContrastName
  )
  stringDbPlotOutputRegistration <- setupStringDbPlotOutputRegistrationFn(
    output = output
  )

  list(
    rawContrastName = rawContrastNameSetup$rawContrastName,
    plotOutputsRegistration = plotOutputsRegistration$registration,
    stringDbPlotOutputRegistration = stringDbPlotOutputRegistration$registration,
    reason = "registered"
  )
}

setupProtEnrichRunOutputDownloadBootstrap <- function(output,
                                                      input,
                                                      enrichmentData,
                                                      workflowData,
                                                      experimentPaths,
                                                      currentAnalysisMethodFn,
                                                      runObserverPreflightFn = runProtEnrichObserverPreflight,
                                                      setupRunObserverRegistrationFn = setupProtEnrichRunObserverRegistration,
                                                      setupResultsSummaryOutputBootstrapFn = setupProtEnrichResultsSummaryOutputBootstrap,
                                                      setupPlotOutputBootstrapFn = setupProtEnrichPlotOutputBootstrap,
                                                      setupResultsDownloadHandlerRegistrationFn = setupProtEnrichResultsDownloadHandlerRegistration) {
  runObserverRegistration <- setupRunObserverRegistrationFn(
    input = input,
    enrichmentData = enrichmentData,
    workflowData = workflowData,
    experimentPaths = experimentPaths,
    currentAnalysisMethodFn = currentAnalysisMethodFn,
    runObserverPreflightFn = runObserverPreflightFn
  )
  resultsSummaryOutputBootstrap <- setupResultsSummaryOutputBootstrapFn(
    output = output,
    input = input,
    enrichmentData = enrichmentData
  )
  plotOutputBootstrap <- setupPlotOutputBootstrapFn(
    output = output,
    input = input,
    enrichmentData = enrichmentData
  )
  resultsDownloadHandlerRegistration <- setupResultsDownloadHandlerRegistrationFn(
    output = output,
    input = input,
    enrichmentData = enrichmentData,
    currentAnalysisMethodFn = currentAnalysisMethodFn
  )

  list(
    runObserverRegistration = runObserverRegistration$registration,
    gprofilerResultsTableOutputRegistration = resultsSummaryOutputBootstrap$gprofilerResultsTableOutputRegistration,
    gprofilerSummaryOutputRegistration = resultsSummaryOutputBootstrap$gprofilerSummaryOutputRegistration,
    clusterProfilerResultsTableOutputRegistration = resultsSummaryOutputBootstrap$clusterProfilerResultsTableOutputRegistration,
    clusterProfilerSummaryOutputRegistration = resultsSummaryOutputBootstrap$clusterProfilerSummaryOutputRegistration,
    stringDbResultsTableOutputRegistration = resultsSummaryOutputBootstrap$stringDbResultsTableOutputRegistration,
    stringDbSummaryOutputRegistration = resultsSummaryOutputBootstrap$stringDbSummaryOutputRegistration,
    rawContrastName = plotOutputBootstrap$rawContrastName,
    plotOutputsRegistration = plotOutputBootstrap$plotOutputsRegistration,
    stringDbPlotOutputRegistration = plotOutputBootstrap$stringDbPlotOutputRegistration,
    resultsDownloadHandlerRegistration = resultsDownloadHandlerRegistration$registration,
    reason = "registered"
  )
}

setupProtEnrichAnalysisMethodDisplayOutputRegistration <- function(output,
                                                                   currentAnalysisMethodFn,
                                                                   registerAnalysisMethodDisplayOutputFn = registerProtEnrichAnalysisMethodDisplayOutput) {
  registration <- registerAnalysisMethodDisplayOutputFn(
    output = output,
    currentAnalysisMethodFn = currentAnalysisMethodFn
  )

  list(
    registration = registration,
    reason = "registered"
  )
}

setupProtEnrichResultsDownloadHandlerRegistration <- function(output,
                                                              input,
                                                              enrichmentData,
                                                              currentAnalysisMethodFn,
                                                              registerResultsDownloadHandlerFn = registerProtEnrichResultsDownloadHandler) {
  registration <- registerResultsDownloadHandlerFn(
    output = output,
    input = input,
    enrichmentData = enrichmentData,
    currentAnalysisMethodFn = currentAnalysisMethodFn
  )

  list(
    registration = registration,
    reason = "registered"
  )
}

setupProtEnrichSelectedTabObserverRegistration <- function(selectedTabFn,
                                                           workflowData,
                                                           enrichmentData,
                                                           session,
                                                           registerSelectedTabObserverFn = registerProtEnrichSelectedTabObserver,
                                                           catFn = cat) {
  registrationState <- list(
    selectedTabProvided = !is.null(selectedTabFn),
    registration = NULL,
    reason = NULL
  )

  if (!is.null(selectedTabFn)) {
    catFn("   mod_prot_enrich_server Step: Setting up tab selection observer\n")
    registrationState$registration <- registerSelectedTabObserverFn(
      selectedTabFn = selectedTabFn,
      workflowData = workflowData,
      enrichmentData = enrichmentData,
      session = session
    )
    registrationState$reason <- "selected_tab_provided"
  } else {
    catFn("   mod_prot_enrich_server Step: No selected_tab parameter provided - tab selection observer NOT set up\n")
    registrationState$reason <- "selected_tab_missing"
  }

  registrationState
}

setupProtEnrichTaxonIdObserverRegistration <- function(workflowData,
                                                       session,
                                                       registerTaxonIdObserverFn = registerProtEnrichTaxonIdObserver) {
  registration <- registerTaxonIdObserverFn(
    workflowData = workflowData,
    session = session
  )

  list(
    registration = registration,
    reason = "registered"
  )
}

setupProtEnrichMixedSpeciesObserverRegistration <- function(workflowData,
                                                            session,
                                                            registerMixedSpeciesObserverFn = registerProtEnrichMixedSpeciesObserver) {
  registration <- registerMixedSpeciesObserverFn(
    workflowData = workflowData,
    session = session
  )

  list(
    registration = registration,
    reason = "registered"
  )
}

setupProtEnrichObserverRegistrationBootstrap <- function(selectedTabFn,
                                                         input,
                                                         workflowData,
                                                         enrichmentData,
                                                         session,
                                                         setupTaxonIdObserverRegistrationFn = setupProtEnrichTaxonIdObserverRegistration,
                                                         setupMixedSpeciesObserverRegistrationFn = setupProtEnrichMixedSpeciesObserverRegistration,
                                                         setupSelectedContrastObserverRegistrationFn = setupProtEnrichSelectedContrastObserverRegistration,
                                                         setupSelectedTabObserverRegistrationFn = setupProtEnrichSelectedTabObserverRegistration,
                                                         setupDaResultsObserverRegistrationFn = setupProtEnrichDaResultsObserverRegistration) {
  taxonIdObserverRegistration <- setupTaxonIdObserverRegistrationFn(
    workflowData = workflowData,
    session = session
  )
  mixedSpeciesObserverRegistration <- setupMixedSpeciesObserverRegistrationFn(
    workflowData = workflowData,
    session = session
  )
  selectedContrastObserverRegistration <- setupSelectedContrastObserverRegistrationFn(
    input = input,
    enrichmentData = enrichmentData
  )
  selectedTabObserverRegistration <- setupSelectedTabObserverRegistrationFn(
    selectedTabFn = selectedTabFn,
    workflowData = workflowData,
    enrichmentData = enrichmentData,
    session = session
  )
  daResultsObserverRegistration <- setupDaResultsObserverRegistrationFn(
    workflowData = workflowData,
    enrichmentData = enrichmentData,
    session = session
  )

  list(
    taxonIdObserverRegistration = taxonIdObserverRegistration$registration,
    mixedSpeciesObserverRegistration = mixedSpeciesObserverRegistration$registration,
    selectedContrastObserverRegistration = selectedContrastObserverRegistration$registration,
    selectedTabObserverRegistration = selectedTabObserverRegistration$registration,
    selectedTabProvided = selectedTabObserverRegistration$selectedTabProvided,
    selectedTabReason = selectedTabObserverRegistration$reason,
    daResultsObserverRegistration = daResultsObserverRegistration$registration,
    reason = "registered"
  )
}

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

runProtEnrichObserverShell <- function(selectedContrast,
                                       runAnalysisBody,
                                       withProgressFn = shiny::withProgress,
                                       finalizeFailureFn = finalizeProtEnrichObserverFailure,
                                       finalizeObserverRunFn = finalizeProtEnrichObserverRun,
                                       progressMessage = "Running enrichment analysis...",
                                       progressValue = 0) {
  analysisCompleted <- FALSE
  failureResult <- NULL
  successResult <- NULL

  tryCatch({
    withProgressFn(message = progressMessage, value = progressValue, {
      runAnalysisBody()
    })
    analysisCompleted <- TRUE
  }, error = function(e) {
    failureResult <<- finalizeFailureFn(errorMessage = e$message)
  })

  if (isTRUE(analysisCompleted)) {
    successResult <- finalizeObserverRunFn(
      completed = TRUE,
      selectedContrast = selectedContrast
    )
  }

  list(
    completed = analysisCompleted,
    selectedContrast = selectedContrast,
    failureResult = failureResult,
    successResult = successResult
  )
}

finalizeProtEnrichObserverFailure <- function(errorMessage,
                                              reportAnalysisErrorFn = reportProtEnrichAnalysisError,
                                              removeWorkingNotificationFn = removeProtEnrichWorkingNotification) {
  reportResult <- reportAnalysisErrorFn(errorMessage = errorMessage)
  cleanupResult <- removeWorkingNotificationFn()

  list(
    errorMessage = errorMessage,
    reportResult = reportResult,
    cleanupResult = cleanupResult
  )
}

reportProtEnrichCompletion <- function(selectedContrast,
                                       notifyCompletionFn = notifyProtEnrichCompletion,
                                       logCompletionFn = logProtEnrichCompletion) {
  notificationResult <- notifyCompletionFn(selectedContrast = selectedContrast)
  logResult <- logCompletionFn()

  list(
    selectedContrast = selectedContrast,
    notificationResult = notificationResult,
    logResult = logResult
  )
}

finalizeProtEnrichObserverRun <- function(completed,
                                          selectedContrast = NULL,
                                          reportCompletionFn = reportProtEnrichCompletion,
                                          removeWorkingNotificationFn = removeProtEnrichWorkingNotification) {
  reportResult <- NULL

  if (isTRUE(completed)) {
    reportResult <- reportCompletionFn(selectedContrast = selectedContrast)
  }

  cleanupResult <- removeWorkingNotificationFn()

  list(
    completed = isTRUE(completed),
    selectedContrast = selectedContrast,
    reportResult = reportResult,
    cleanupResult = cleanupResult
  )
}

logProtEnrichCompletion <- function(message = "=== ENRICHMENT ANALYSIS COMPLETED ===\n",
                                    catFn = cat) {
  catFn(message)

  list(message = message)
}

removeProtEnrichWorkingNotification <- function(notificationId = "enrichment_working",
                                                removeNotificationFn = shiny::removeNotification) {
  removeNotificationFn(notificationId)

  list(notificationId = notificationId)
}

handoffProtEnrichObserverRun <- function(input,
                                         enrichmentData,
                                         workflowData,
                                         experimentPaths,
                                         currentAnalysisMethodFn,
                                         showNotificationFn = shiny::showNotification,
                                         runObserverShellFn = runProtEnrichObserverShell,
                                         runAnalysisBodyFn = runProtEnrichAnalysisBody) {
  selectedContrast <- input$selected_contrast
  notificationMessage <- "Running enrichment analysis..."
  notificationId <- "enrichment_working"
  notificationDuration <- NULL

  showNotificationFn(
    notificationMessage,
    id = notificationId,
    duration = notificationDuration
  )

  shellResult <- runObserverShellFn(
    selectedContrast = selectedContrast,
    runAnalysisBody = function() {
      runAnalysisBodyFn(
        input = input,
        enrichmentData = enrichmentData,
        workflowData = workflowData,
        experimentPaths = experimentPaths,
        currentAnalysisMethodFn = currentAnalysisMethodFn
      )
    }
  )

  list(
    selectedContrast = selectedContrast,
    notificationMessage = notificationMessage,
    notificationId = notificationId,
    notificationDuration = notificationDuration,
    shellResult = shellResult
  )
}

runProtEnrichObserverPreflight <- function(input,
                                           enrichmentData,
                                           workflowData,
                                           experimentPaths,
                                           currentAnalysisMethodFn,
                                           reqFn = shiny::req,
                                           handoffObserverRunFn = handoffProtEnrichObserverRun) {
  selectedContrast <- input$selected_contrast
  daResultsData <- enrichmentData$da_results_data

  reqFn(selectedContrast, daResultsData)

  handoffResult <- handoffObserverRunFn(
    input = input,
    enrichmentData = enrichmentData,
    workflowData = workflowData,
    experimentPaths = experimentPaths,
    currentAnalysisMethodFn = currentAnalysisMethodFn
  )

  list(
    selectedContrast = selectedContrast,
    daResultsData = daResultsData,
    handoffResult = handoffResult
  )
}

buildProtEnrichSupportedOrganisms <- function() {
  tibble::tribble(
    ~taxid,     ~id,             ~name,
    "9606",     "hsapiens",      "Homo sapiens",
    "10090",    "mmusculus",     "Mus musculus",
    "10116",    "rnorvegicus",   "Rattus norvegicus",
    "7227",     "dmelanogaster", "Drosophila melanogaster",
    "6239",     "celegans",      "Caenorhabditis elegans",
    "4932",     "scerevisiae",   "Saccharomyces cerevisiae",
    "3702",     "athaliana",     "Arabidopsis thaliana",
    "7955",     "drerio",        "Danio rerio",
    "9031",     "ggallus",       "Gallus gallus",
    "9823",     "sscrofa",       "Sus scrofa",
    "9913",     "btaurus",       "Bos taurus",
    "9544",     "mmulatta",      "Macaca mulatta",
    "9598",     "ptroglodytes",  "Pan troglodytes"
  )
}

createProtEnrichSupportedOrganismsReactive <- function(reactiveFn = shiny::reactive,
                                                       buildSupportedOrganismsFn = buildProtEnrichSupportedOrganisms) {
  reactiveFn({
    buildSupportedOrganismsFn()
  })
}

createProtEnrichReactiveValues <- function(reactiveValuesFn = shiny::reactiveValues) {
  reactiveValuesFn(
    enrichment_results = NULL,
    contrasts_available = NULL,
    analysis_complete = FALSE,
    current_s4_object = NULL,
    da_results_data = NULL,
    gprofiler_results = NULL,
    clusterprofiler_results = NULL,
    stringdb_results = NULL,
    analysis_method = NULL,
    organism_supported = NULL,
    all_enrichment_results = list(),
    current_contrast_results = list(),
    enrichment_plots = list()
  )
}

resolveProtEnrichAnalysisMethod <- function(organismTaxid, supportedOrganisms) {
  isSupported <- organismTaxid %in% supportedOrganisms$taxid

  if (isTRUE(isSupported)) {
    speciesInfo <- supportedOrganisms |>
      dplyr::filter(taxid == organismTaxid)

    methodInfo <- list(
      method = "gprofiler2",
      supported = TRUE,
      species_id = speciesInfo$id[1],
      species_name = speciesInfo$name[1],
      description = paste("gprofiler2 analysis for", speciesInfo$name[1])
    )
  } else {
    methodInfo <- list(
      method = "clusterprofiler",
      supported = FALSE,
      species_id = NULL,
      species_name = paste("Taxon ID", organismTaxid),
      description = paste("clusterProfileR analysis with custom GO annotations for taxon", organismTaxid)
    )
  }

  list(
    methodInfo = methodInfo,
    analysisMethod = methodInfo$method,
    organismSupported = methodInfo$supported
  )
}

createProtEnrichCurrentAnalysisMethodReactive <- function(input,
                                                          enrichmentData,
                                                          supportedOrganismsFn,
                                                          reactiveFn = shiny::reactive,
                                                          reqFn = shiny::req,
                                                          resolveAnalysisMethodFn = resolveProtEnrichAnalysisMethod) {
  reactiveFn({
    reqFn(input$organism_taxid)

    methodState <- resolveAnalysisMethodFn(
      organismTaxid = input$organism_taxid,
      supportedOrganisms = supportedOrganismsFn()
    )

    enrichmentData$analysis_method <- methodState$analysisMethod
    enrichmentData$organism_supported <- methodState$organismSupported

    methodState$methodInfo
  })
}

formatProtEnrichStatusText <- function(analysisComplete,
                                       methodInfo = NULL,
                                       selectedContrast = NULL,
                                       upCutoff = NULL,
                                       downCutoff = NULL,
                                       qCutoff = NULL,
                                       gprofilerResults = NULL,
                                       clusterprofilerResults = NULL,
                                       stringdbResults = NULL) {
  if (isTRUE(analysisComplete)) {
    gprofilerCount <- if (!is.null(gprofilerResults)) nrow(gprofilerResults) else 0
    clusterprofilerCount <- if (!is.null(clusterprofilerResults)) nrow(clusterprofilerResults) else 0
    stringdbCount <- if (!is.null(stringdbResults)) nrow(stringdbResults) else 0

    return(paste(
      "[OK] Analysis Complete\n",
      sprintf("Method: %s\n", methodInfo$method),
      sprintf("Contrast: %s\n", selectedContrast),
      sprintf("Up log2FC cutoff: %.1f\n", upCutoff),
      sprintf("Down log2FC cutoff: %.1f\n", downCutoff),
      sprintf("Q-value cutoff: %.3f\n", qCutoff),
      sprintf("Organism: %s\n", methodInfo$species_name),
      "",
      "Results Available:",
      sprintf("* gprofiler2: %d terms", gprofilerCount),
      sprintf("* clusterProfileR: %d terms", clusterprofilerCount),
      sprintf("* STRING-DB: %d networks", stringdbCount),
      "",
      "[OK] Results saved to workflow state",
      sep = "\n"
    ))
  }

  paste(
    "[WAITING] Ready for analysis\n",
    "",
    "Steps:",
    "1. Select contrast from DE results",
    "2. Set log fold change cutoffs",
    "3. Set Q-value cutoff (significance threshold)",
    "4. Click 'Run Enrichment Analysis'",
    "",
    "Method automatically determined by organism.",
    sep = "\n"
  )
}

formatProtEnrichAnalysisMethodText <- function(methodInfo) {
  if (isTRUE(methodInfo$supported)) {
    return(paste(
      "[OK] SUPPORTED ORGANISM\n",
      sprintf("Method: %s\n", methodInfo$method),
      sprintf("Species: %s\n", methodInfo$species_name),
      "All enrichment methods available"
    ))
  }

  paste(
    "[WARNING] CUSTOM ORGANISM\n",
    sprintf("Method: %s\n", methodInfo$method),
    sprintf("Organism: %s\n", methodInfo$species_name),
    "Using UniProt GO annotations"
  )
}

formatProtEnrichContrastsText <- function(contrastsAvailable) {
  if (!is.null(contrastsAvailable)) {
    return(paste(contrastsAvailable, collapse = "\n"))
  }

  "No contrasts available.\nComplete differential expression\nanalysis first."
}

formatProtEnrichGprofilerSummaryText <- function(gprofilerResults,
                                                 directionFilter = "all") {
  if (is.null(gprofilerResults) || nrow(gprofilerResults) == 0) {
    return("No gprofiler2 results available.")
  }

  tryCatch({
    if (!identical(directionFilter, "all") &&
        "directionality" %in% names(gprofilerResults)) {
      directionValue <- if (identical(directionFilter, "up")) "positive" else "negative"
      filteredResults <- gprofilerResults |> dplyr::filter(directionality == directionValue)
      displayedCount <- nrow(filteredResults)

      if (identical(directionFilter, "up")) {
        messageText <- sprintf("Showing %d up-regulated pathways", displayedCount)
      } else {
        messageText <- sprintf("Showing %d down-regulated pathways", displayedCount)
      }
    } else {
      totalTerms <- nrow(gprofilerResults)
      positiveTerms <- sum(gprofilerResults$directionality == "positive", na.rm = TRUE)
      negativeTerms <- sum(gprofilerResults$directionality == "negative", na.rm = TRUE)

      messageText <- paste(
        sprintf("Total enrichment terms: %d", totalTerms),
        sprintf("Up-regulated pathways: %d", positiveTerms),
        sprintf("Down-regulated pathways: %d", negativeTerms),
        sep = "\n"
      )
    }

    paste(
      messageText,
      "",
      "Results displayed in table below.",
      sep = "\n"
    )
  }, error = function(e) {
    paste("Error calculating statistics:", e$message)
  })
}

formatProtEnrichClusterProfilerSummaryText <- function(clusterprofilerResults,
                                                       directionFilter = "all") {
  if (is.null(clusterprofilerResults) || nrow(clusterprofilerResults) == 0) {
    return("No clusterProfileR results available.")
  }

  tryCatch({
    if (!identical(directionFilter, "all") &&
        "directionality" %in% names(clusterprofilerResults)) {
      filteredResults <- clusterprofilerResults |>
        dplyr::filter(directionality == directionFilter)
      displayedCount <- nrow(filteredResults)

      if (identical(directionFilter, "up")) {
        messageText <- sprintf("Showing %d up-regulated GO terms", displayedCount)
      } else {
        messageText <- sprintf("Showing %d down-regulated GO terms", displayedCount)
      }
    } else {
      totalTerms <- nrow(clusterprofilerResults)
      upTerms <- sum(clusterprofilerResults$directionality == "up", na.rm = TRUE)
      downTerms <- sum(clusterprofilerResults$directionality == "down", na.rm = TRUE)

      messageText <- paste(
        sprintf("Total GO terms: %d", totalTerms),
        sprintf("Up-regulated: %d", upTerms),
        sprintf("Down-regulated: %d", downTerms),
        sep = "\n"
      )
    }

    paste(
      messageText,
      "",
      "Results displayed in table below.",
      sep = "\n"
    )
  }, error = function(e) {
    paste("Error calculating statistics:", e$message)
  })
}

formatProtEnrichStringDbSummaryText <- function(stringdbResults) {
  if (is.null(stringdbResults) || nrow(stringdbResults) == 0) {
    return("STRING-DB analysis not yet implemented.\nThis will show network enrichment statistics.")
  }

  paste(
    "STRING-DB Network Analysis",
    "Status: Implementation pending",
    "",
    "Features planned:",
    "* Protein-protein interaction networks",
    "* Functional cluster identification",
    "* Network topology analysis",
    "* Interactive network visualization",
    sep = "\n"
  )
}

renderProtEnrichGprofilerResultsTable <- function(gprofilerResults,
                                                  directionFilter = "all",
                                                  renderDtFn = DT::renderDT,
                                                  datatableFn = DT::datatable,
                                                  formatRoundFn = DT::formatRound,
                                                  catFn = cat) {
  renderDtFn({
    if (is.null(gprofilerResults) || nrow(gprofilerResults) == 0) {
      return(datatableFn(data.frame(
        Message = "No gprofiler2 results available. Run analysis first."
      )))
    }

    tryCatch({
      currentResults <- gprofilerResults

      if (!identical(directionFilter, "all") &&
          "directionality" %in% names(currentResults)) {
        directionValue <- if (identical(directionFilter, "up")) "positive" else "negative"
        currentResults <- currentResults |>
          dplyr::filter(directionality == directionValue)
      }

      datatableFn(
        currentResults,
        options = list(
          pageLength = 25,
          scrollX = TRUE,
          dom = "Bfrtip",
          buttons = c("copy", "csv", "excel")
        ),
        extensions = "Buttons"
      ) |>
        formatRoundFn(
          columns = intersect(c("pvalue", "p.adjust", "qvalue"), names(currentResults)),
          digits = 4
        )
    }, error = function(e) {
      catFn(paste("*** ERROR in gprofiler2 results table:", e$message, "\n"))
      datatableFn(data.frame(Message = paste("Error:", e$message)))
    })
  })
}

renderProtEnrichGprofilerPlot <- function(analysisComplete,
                                          enrichmentResultsFull,
                                          rawContrast,
                                          directionFilter = "all",
                                          renderPlotlyFn = plotly::renderPlotly,
                                          plotLyFn = plotly::plot_ly,
                                          addTextFn = plotly::add_text) {
  renderPlotlyFn({
    buildPlaceholderPlot <- function(messageText) {
      plotLyFn() |>
        addTextFn(
          x = 0.5,
          y = 0.5,
          text = messageText,
          showlegend = FALSE
        )
    }

    if (!analysisComplete) {
      return(buildPlaceholderPlot("Run enrichment analysis first"))
    }

    if (is.null(enrichmentResultsFull)) {
      return(buildPlaceholderPlot("No enrichment results"))
    }

    tryCatch({
      plots <- enrichmentResultsFull@enrichment_plotly[[rawContrast]]

      if (identical(directionFilter, "up") && !is.null(plots$up)) {
        return(plots$up)
      }

      if (identical(directionFilter, "down") && !is.null(plots$down)) {
        return(plots$down)
      }

      if (identical(directionFilter, "all")) {
        if (!is.null(plots$up)) {
          return(plots$up)
        }

        if (!is.null(plots$down)) {
          return(plots$down)
        }
      }

      if (identical(directionFilter, "up")) {
        return(buildPlaceholderPlot("No up-regulated enrichment data"))
      }

      if (identical(directionFilter, "down")) {
        return(buildPlaceholderPlot("No down-regulated enrichment data"))
      }

      buildPlaceholderPlot("No plot data for this contrast")
    }, error = function(e) {
      buildPlaceholderPlot(paste("Plot error:", e$message))
    })
  })
}

renderProtEnrichClusterProfilerPlot <- function(analysisComplete,
                                                enrichmentResultsFull,
                                                rawContrast,
                                                directionFilter = "all",
                                                renderPlotlyFn = plotly::renderPlotly,
                                                plotLyFn = plotly::plot_ly,
                                                addTextFn = plotly::add_text) {
  renderPlotlyFn({
    buildPlaceholderPlot <- function(messageText) {
      plotLyFn() |>
        addTextFn(
          x = 0.5,
          y = 0.5,
          text = messageText,
          showlegend = FALSE
        )
    }

    if (!analysisComplete) {
      return(buildPlaceholderPlot("Run enrichment analysis first"))
    }

    if (is.null(enrichmentResultsFull)) {
      return(buildPlaceholderPlot("No enrichment results"))
    }

    tryCatch({
      plots <- enrichmentResultsFull@enrichment_plotly[[rawContrast]]

      if (identical(directionFilter, "up") && !is.null(plots$up)) {
        return(plots$up)
      }

      if (identical(directionFilter, "down") && !is.null(plots$down)) {
        return(plots$down)
      }

      if (identical(directionFilter, "all")) {
        if (!is.null(plots$up)) {
          return(plots$up)
        }

        if (!is.null(plots$down)) {
          return(plots$down)
        }
      }

      if (identical(directionFilter, "up")) {
        return(buildPlaceholderPlot("No up-regulated enrichment data"))
      }

      if (identical(directionFilter, "down")) {
        return(buildPlaceholderPlot("No down-regulated enrichment data"))
      }

      buildPlaceholderPlot("No plot data for this contrast")
    }, error = function(e) {
      buildPlaceholderPlot(paste("Plot error:", e$message))
    })
  })
}

renderProtEnrichClusterProfilerResultsTable <- function(clusterprofilerResults,
                                                        directionFilter = "all",
                                                        renderDtFn = DT::renderDT,
                                                        datatableFn = DT::datatable,
                                                        formatRoundFn = DT::formatRound,
                                                        catFn = cat) {
  renderDtFn({
    if (is.null(clusterprofilerResults) || nrow(clusterprofilerResults) == 0) {
      return(datatableFn(data.frame(
        Message = "No clusterProfileR results available. Run analysis first."
      )))
    }

    tryCatch({
      currentResults <- clusterprofilerResults

      if (!identical(directionFilter, "all") &&
          "directionality" %in% names(currentResults)) {
        currentResults <- currentResults |>
          dplyr::filter(directionality == directionFilter)
      }

      datatableFn(
        currentResults,
        options = list(
          pageLength = 25,
          scrollX = TRUE,
          dom = "Bfrtip",
          buttons = c("copy", "csv", "excel")
        ),
        extensions = "Buttons"
      ) |>
        formatRoundFn(
          columns = intersect(c("pvalue", "p.adjust", "qvalue"), names(currentResults)),
          digits = 4
        )
    }, error = function(e) {
      catFn(sprintf("*** ERROR in clusterProfileR results table: %s ***\n", e$message))
      datatableFn(data.frame(Message = paste("Error:", e$message)))
    })
  })
}

renderProtEnrichStringDbResultsTable <- function(stringdbResults,
                                                 filterSignificant = FALSE,
                                                 enrichmentPValThresh = NULL,
                                                 maxResults = Inf,
                                                 renderDtFn = DT::renderDT,
                                                 datatableFn = DT::datatable) {
  renderDtFn({
    if (is.null(stringdbResults) || nrow(stringdbResults) == 0) {
      return(datatableFn(data.frame(
        Message = "STRING-DB analysis not yet implemented.",
        Note = "This tab will show protein-protein interaction network enrichment results.",
        Status = "Coming soon..."
      )))
    }

    tryCatch({
      currentResults <- stringdbResults

      if (isTRUE(filterSignificant)) {
        currentResults <- currentResults |>
          dplyr::filter(p_value < enrichmentPValThresh)
      }

      if (nrow(currentResults) > maxResults) {
        currentResults <- currentResults |>
          dplyr::slice_head(n = maxResults)
      }

      datatableFn(
        currentResults,
        options = list(
          pageLength = 25,
          scrollX = TRUE,
          dom = "Bfrtip",
          buttons = c("copy", "csv", "excel")
        ),
        extensions = "Buttons"
      )
    }, error = function(e) {
      datatableFn(data.frame(Message = paste("STRING-DB Error:", e$message)))
    })
  })
}

renderProtEnrichStringDbPlot <- function(renderPlotlyFn = plotly::renderPlotly,
                                         plotLyFn = plotly::plot_ly,
                                         addTextFn = plotly::add_text) {
  renderPlotlyFn({
    plotLyFn() |>
      addTextFn(
        x = 0.5,
        y = 0.5,
        text = "STRING-DB coming soon",
        showlegend = FALSE
      )
  })
}

buildProtEnrichResultsDownloadFilename <- function(selectedContrast,
                                                   date = Sys.Date(),
                                                   prefix = "Enrichment_results") {
  contrastSafe <- gsub("[^A-Za-z0-9_.-]", "_", selectedContrast)
  paste0(prefix, "_", contrastSafe, "_", date, ".zip")
}

writeProtEnrichResultsDownloadArchive <- function(file,
                                                  selectedContrast,
                                                  methodInfo,
                                                  organismTaxid,
                                                  upCutoff,
                                                  downCutoff,
                                                  qCutoff,
                                                  gprofilerResults = NULL,
                                                  clusterprofilerResults = NULL,
                                                  stringdbResults = NULL,
                                                  tempDir = tempdir(),
                                                  writeTsvFn = readr::write_tsv,
                                                  writeLinesFn = writeLines,
                                                  zipFn = utils::zip,
                                                  sysTimeFn = Sys.time) {
  tryCatch({
    filesToZip <- character()

    if (!is.null(gprofilerResults) && nrow(gprofilerResults) > 0) {
      gprofilerFile <- file.path(tempDir, "gprofiler2_results.tsv")
      writeTsvFn(gprofilerResults, gprofilerFile)
      filesToZip <- c(filesToZip, gprofilerFile)
    }

    if (!is.null(clusterprofilerResults) && nrow(clusterprofilerResults) > 0) {
      clusterprofilerFile <- file.path(tempDir, "clusterProfileR_results.tsv")
      writeTsvFn(clusterprofilerResults, clusterprofilerFile)
      filesToZip <- c(filesToZip, clusterprofilerFile)
    }

    if (!is.null(stringdbResults) && nrow(stringdbResults) > 0) {
      stringdbFile <- file.path(tempDir, "stringdb_results.tsv")
      writeTsvFn(stringdbResults, stringdbFile)
      filesToZip <- c(filesToZip, stringdbFile)
    }

    summaryFile <- file.path(tempDir, "enrichment_analysis_summary.txt")
    summaryContent <- paste(
      "# MultiScholaR Enrichment Analysis Results",
      paste("Date:", sysTimeFn()),
      paste("Contrast:", selectedContrast),
      paste("Analysis Method:", methodInfo$method),
      paste("Organism:", methodInfo$species_name),
      paste("Taxonomy ID:", organismTaxid),
      "",
      "## Analysis Parameters:",
      paste("Up log2FC Cutoff:", upCutoff),
      paste("Down log2FC Cutoff:", downCutoff),
      paste("Q-value Cutoff:", qCutoff),
      "",
      "## Results Summary:",
      paste("gprofiler2 terms:", if (!is.null(gprofilerResults)) nrow(gprofilerResults) else 0),
      paste("clusterProfileR terms:", if (!is.null(clusterprofilerResults)) nrow(clusterprofilerResults) else 0),
      paste("STRING-DB networks:", if (!is.null(stringdbResults)) nrow(stringdbResults) else 0),
      "",
      "## Files Included:",
      if (length(filesToZip) > 0) paste("*", basename(filesToZip), collapse = "\n") else "* No result files (analysis may have failed)",
      "",
      "Generated by MultiScholaR Enrichment Analysis Module",
      sep = "\n"
    )

    writeLinesFn(summaryContent, summaryFile)
    filesToZip <- c(filesToZip, summaryFile)

    if (length(filesToZip) > 0) {
      zipFn(zipfile = file, files = filesToZip, flags = "-j")
    } else {
      noteFile <- file.path(tempDir, "no_results.txt")
      writeLinesFn("No enrichment results available to download.", noteFile)
      zipFn(zipfile = file, files = noteFile, flags = "-j")
      filesToZip <- noteFile
    }

    invisible(list(
      status = "ok",
      files = filesToZip,
      zipfile = file
    ))
  }, error = function(e) {
    errorFile <- file.path(tempDir, "download_error.txt")
    writeLinesFn(paste("Error creating download:", e$message), errorFile)
    zipFn(zipfile = file, files = errorFile, flags = "-j")

    invisible(list(
      status = "error",
      files = errorFile,
      zipfile = file,
      error = e$message
    ))
  })
}

buildProtEnrichAnalysisResultsPayload <- function(gprofilerResults,
                                                  clusterprofilerResults,
                                                  stringdbResults,
                                                  fullEnrichmentResults,
                                                  selectedContrast,
                                                  analysisMethod,
                                                  upCutoff,
                                                  downCutoff,
                                                  qCutoff,
                                                  organismTaxid) {
  list(
    gprofiler_results = gprofilerResults,
    clusterprofiler_results = clusterprofilerResults,
    stringdb_results = stringdbResults,
    full_enrichment_results = fullEnrichmentResults,
    selected_contrast = selectedContrast,
    analysis_method = analysisMethod,
    parameters = list(
      up_cutoff = upCutoff,
      down_cutoff = downCutoff,
      q_cutoff = qCutoff,
      organism_taxid = organismTaxid
    )
  )
}

propagateProtEnrichResultsArgs <- function(enrichmentResults,
                                           currentS4Object,
                                           selectedContrast,
                                           methodInfo,
                                           upCutoff,
                                           downCutoff,
                                           qCutoff,
                                           organismTaxid,
                                           pathwayDir,
                                           timeFn = Sys.time,
                                           catFn = cat) {
  catFn("   ENRICHMENT Step: Copying @args from original data object...\n")

  dataHasArgs <- tryCatch({
    !is.null(currentS4Object) && !is.null(currentS4Object@args)
  }, error = function(e) {
    FALSE
  })

  resultsHasArgs <- tryCatch({
    !is.null(enrichmentResults@args)
  }, error = function(e) {
    FALSE
  })

  if (dataHasArgs) {
    if (resultsHasArgs) {
      enrichmentResults@args <- currentS4Object@args

      if (is.null(enrichmentResults@args$enrichmentAnalysis)) {
        enrichmentResults@args$enrichmentAnalysis <- list()
      }

      enrichmentResults@args$enrichmentAnalysis <- list(
        selected_contrast = selectedContrast,
        analysis_method = methodInfo$method,
        organism_supported = methodInfo$supported,
        up_cutoff = upCutoff,
        down_cutoff = downCutoff,
        q_cutoff = qCutoff,
        organism_taxid = organismTaxid,
        pathway_dir = pathwayDir
      )

      catFn("   ENRICHMENT Step: Storing UI parameters in @args\n")
      if (is.null(enrichmentResults@args$enrichmentAnalysisUI)) {
        enrichmentResults@args$enrichmentAnalysisUI <- list()
      }

      enrichmentResults@args$enrichmentAnalysisUI <- list(
        up_log2fc_cutoff = upCutoff,
        down_log2fc_cutoff = downCutoff,
        q_value_cutoff = qCutoff,
        organism_taxon_id = organismTaxid,
        analysis_method = methodInfo$method,
        organism_name = methodInfo$species_name,
        organism_supported = methodInfo$supported,
        selected_contrast = selectedContrast,
        timestamp = timeFn()
      )

      catFn("   ENRICHMENT Step: Successfully copied and updated @args\n")
    } else {
      catFn("   ENRICHMENT Step: EnrichmentResults doesn't have @args slot\n")
    }
  } else {
    catFn("   ENRICHMENT Step: Original data object doesn't have @args to copy\n")
  }

  list(
    enrichmentResults = enrichmentResults,
    dataHasArgs = dataHasArgs,
    resultsHasArgs = resultsHasArgs,
    copiedArgs = isTRUE(dataHasArgs) && isTRUE(resultsHasArgs)
  )
}

propagateProtEnrichUiParams <- function(currentS4Object,
                                        workflowData,
                                        selectedContrast,
                                        methodInfo,
                                        upCutoff,
                                        downCutoff,
                                        qCutoff,
                                        organismTaxid,
                                        timeFn = Sys.time,
                                        catFn = cat) {
  dataHasArgs <- tryCatch({
    !is.null(currentS4Object) && !is.null(currentS4Object@args)
  }, error = function(e) {
    FALSE
  })

  if (dataHasArgs) {
    catFn("   ENRICHMENT Step: Storing UI parameters in original data object @args\n")
    if (is.null(currentS4Object@args$enrichmentAnalysisUI)) {
      currentS4Object@args$enrichmentAnalysisUI <- list()
    }

    currentS4Object@args$enrichmentAnalysisUI <- list(
      up_log2fc_cutoff = upCutoff,
      down_log2fc_cutoff = downCutoff,
      q_value_cutoff = qCutoff,
      organism_taxon_id = organismTaxid,
      analysis_method = methodInfo$method,
      organism_name = methodInfo$species_name,
      organism_supported = methodInfo$supported,
      selected_contrast = selectedContrast,
      timestamp = timeFn()
    )

    workflowData$enrichment_ui_params <- list(
      up_log2fc_cutoff = upCutoff,
      down_log2fc_cutoff = downCutoff,
      q_value_cutoff = qCutoff,
      organism_selected = organismTaxid,
      database_source = methodInfo$method,
      organism_name = methodInfo$species_name,
      organism_supported = methodInfo$supported,
      selected_contrast = selectedContrast,
      timestamp = timeFn()
    )
    catFn("   ENRICHMENT Step: Stored UI parameters in workflow_data for sessionSummary\n")
  }

  list(
    currentS4Object = currentS4Object,
    dataHasArgs = dataHasArgs,
    storedUiParams = isTRUE(dataHasArgs),
    workflowUiParams = if (isTRUE(dataHasArgs)) workflowData$enrichment_ui_params else NULL
  )
}

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

completeProtEnrichProgress <- function(value = 1.0,
                                       detail = "Complete!",
                                       incProgressFn = shiny::incProgress) {
  incProgressFn(value, detail = detail)

  list(
    value = value,
    detail = detail
  )
}

notifyProtEnrichCompletion <- function(selectedContrast,
                                       type = "message",
                                       duration = 5,
                                       showNotificationFn = shiny::showNotification) {
  notificationMessage <- paste(
    "Enrichment analysis completed successfully for contrast:",
    selectedContrast
  )

  showNotificationFn(
    notificationMessage,
    type = type,
    duration = duration
  )

  list(
    message = notificationMessage,
    type = type,
    duration = duration,
    selectedContrast = selectedContrast
  )
}

notifyProtEnrichAnalysisError <- function(errorMessage,
                                          type = "error",
                                          duration = 10,
                                          showNotificationFn = shiny::showNotification) {
  notificationMessage <- sprintf("Error in enrichment analysis: %s", errorMessage)

  showNotificationFn(
    notificationMessage,
    type = type,
    duration = duration
  )

  list(
    message = notificationMessage,
    type = type,
    duration = duration,
    errorMessage = errorMessage
  )
}

logProtEnrichAnalysisError <- function(errorMessage,
                                       template = "*** ERROR in enrichment analysis: %s ***\n",
                                       catFn = cat) {
  message <- sprintf(template, errorMessage)
  catFn(message)

  list(
    message = message,
    errorMessage = errorMessage
  )
}

messageProtEnrichAnalysisError <- function(errorMessage,
                                           template = "*** ERROR in enrichment analysis: %s",
                                           messageFn = message) {
  formattedMessage <- sprintf(template, errorMessage)
  messageFn(formattedMessage)

  list(
    message = formattedMessage,
    errorMessage = errorMessage
  )
}

reportProtEnrichAnalysisError <- function(errorMessage,
                                          messageErrorFn = messageProtEnrichAnalysisError,
                                          logErrorFn = logProtEnrichAnalysisError,
                                          notifyErrorFn = notifyProtEnrichAnalysisError) {
  messageResult <- messageErrorFn(errorMessage = errorMessage)
  logResult <- logErrorFn(errorMessage = errorMessage)
  notificationResult <- notifyErrorFn(errorMessage = errorMessage)

  list(
    errorMessage = errorMessage,
    messageResult = messageResult,
    logResult = logResult,
    notificationResult = notificationResult
  )
}
captureProtEnrichPostProcessResults <- function(selectedContrast,
                                                enrichmentResults,
                                                enrichmentData,
                                                contrastsTbl,
                                                methodInfo,
                                                buildAllContrastResultsFn = buildProtEnrichAllContrastResults,
                                                resolveSelectedContrastResultsFn = resolveProtEnrichSelectedContrastResults,
                                                catFn = cat) {
  if (!is.null(enrichmentResults) && !is.null(enrichmentResults@enrichment_data)) {
    allContrastResults <- buildAllContrastResultsFn(
      enrichmentResults = enrichmentResults,
      methodInfo = methodInfo
    )

    enrichmentData$all_enrichment_results <- allContrastResults

    initialContrastState <- resolveSelectedContrastResultsFn(
      selectedContrast,
      allContrastResults,
      contrastsTbl
    )

    if (initialContrastState$found) {
      enrichmentData$gprofiler_results <- initialContrastState$gprofilerResults
      enrichmentData$clusterprofiler_results <- initialContrastState$clusterprofilerResults
      enrichmentData$stringdb_results <- initialContrastState$stringdbResults
      catFn(sprintf(
        "   ENRICHMENT Step: Set initial display to contrast '%s'\n",
        initialContrastState$rawContrastName
      ))
    }

    return(list(
      allContrastResults = allContrastResults,
      initialContrastState = initialContrastState,
      hasResults = TRUE
    ))
  }

  catFn("   ENRICHMENT Step: processEnrichments returned NULL or empty results\n")

  list(
    allContrastResults = NULL,
    initialContrastState = NULL,
    hasResults = FALSE
  )
}

persistProtEnrichAnalysisResults <- function(input,
                                             enrichmentData,
                                             workflowData,
                                             selectedContrast,
                                             enrichmentResults,
                                             methodInfo,
                                             pathwayDir,
                                             buildAnalysisResultsPayloadFn = buildProtEnrichAnalysisResultsPayload,
                                             propagateResultsArgsFn = propagateProtEnrichResultsArgs,
                                             propagateUiParamsFn = propagateProtEnrichUiParams,
                                             updateStateManagerUiParamsFn = updateProtEnrichStateManagerUiParams,
                                             saveCompletedStateFn = saveProtEnrichCompletedState,
                                             completeTabStatusFn = completeProtEnrichTabStatus,
                                             completeProgressFn = completeProtEnrichProgress,
                                             incProgressFn = shiny::incProgress) {
  incProgressFn(0.8, detail = "Storing results...")

  enrichmentData$enrichment_results_full <- enrichmentResults
  enrichmentData$analysis_complete <- TRUE

  workflowData$enrichment_analysis_results <- buildAnalysisResultsPayloadFn(
    gprofilerResults = enrichmentData$gprofiler_results,
    clusterprofilerResults = enrichmentData$clusterprofiler_results,
    stringdbResults = enrichmentData$stringdb_results,
    fullEnrichmentResults = enrichmentData$enrichment_results_full,
    selectedContrast = selectedContrast,
    analysisMethod = methodInfo$method,
    upCutoff = input$up_cutoff,
    downCutoff = input$down_cutoff,
    qCutoff = input$q_cutoff,
    organismTaxid = input$organism_taxid
  )

  argsPropagation <- propagateResultsArgsFn(
    enrichmentResults = enrichmentResults,
    currentS4Object = enrichmentData$current_s4_object,
    selectedContrast = selectedContrast,
    methodInfo = methodInfo,
    upCutoff = input$up_cutoff,
    downCutoff = input$down_cutoff,
    qCutoff = input$q_cutoff,
    organismTaxid = input$organism_taxid,
    pathwayDir = pathwayDir
  )
  enrichmentResults <- argsPropagation$enrichmentResults

  uiParamsPropagation <- propagateUiParamsFn(
    currentS4Object = enrichmentData$current_s4_object,
    workflowData = workflowData,
    selectedContrast = selectedContrast,
    methodInfo = methodInfo,
    upCutoff = input$up_cutoff,
    downCutoff = input$down_cutoff,
    qCutoff = input$q_cutoff,
    organismTaxid = input$organism_taxid
  )
  enrichmentData$current_s4_object <- uiParamsPropagation$currentS4Object

  updateStateManagerUiParamsFn(
    workflowData = workflowData,
    storedUiParams = uiParamsPropagation$storedUiParams
  )

  saveCompletedStateFn(
    workflowData = workflowData,
    enrichmentResults = enrichmentResults,
    selectedContrast = selectedContrast,
    methodInfo = methodInfo,
    upCutoff = input$up_cutoff,
    downCutoff = input$down_cutoff,
    qCutoff = input$q_cutoff,
    organismTaxid = input$organism_taxid,
    pathwayDir = pathwayDir
  )

  completeTabStatusFn(workflowData = workflowData)
  completeProgressFn()

  list(
    enrichmentResults = enrichmentResults,
    analysisComplete = enrichmentData$analysis_complete,
    currentS4Object = enrichmentData$current_s4_object
  )
}

finalizeProtEnrichAnalysisBodyResults <- function(selectedContrast,
                                                  rawContrastName,
                                                  organismFilterApplied,
                                                  filterStats,
                                                  enrichmentResults,
                                                  enrichmentData,
                                                  workflowData,
                                                  input,
                                                  methodInfo,
                                                  contrastsTbl,
                                                  pathwayDir,
                                                  buildAllContrastResultsFn = buildProtEnrichAllContrastResults,
                                                  resolveSelectedContrastResultsFn = resolveProtEnrichSelectedContrastResults,
                                                  capturePostProcessResultsFn = captureProtEnrichPostProcessResults,
                                                  persistAnalysisResultsFn = persistProtEnrichAnalysisResults,
                                                  buildAnalysisResultsPayloadFn = buildProtEnrichAnalysisResultsPayload,
                                                  propagateResultsArgsFn = propagateProtEnrichResultsArgs,
                                                  propagateUiParamsFn = propagateProtEnrichUiParams,
                                                  updateStateManagerUiParamsFn = updateProtEnrichStateManagerUiParams,
                                                  saveCompletedStateFn = saveProtEnrichCompletedState,
                                                  completeTabStatusFn = completeProtEnrichTabStatus,
                                                  completeProgressFn = completeProtEnrichProgress,
                                                  incProgressFn = shiny::incProgress,
                                                  catFn = cat) {
  capturePostProcessResultsFn(
    selectedContrast = selectedContrast,
    enrichmentResults = enrichmentResults,
    enrichmentData = enrichmentData,
    contrastsTbl = contrastsTbl,
    methodInfo = methodInfo,
    buildAllContrastResultsFn = buildAllContrastResultsFn,
    resolveSelectedContrastResultsFn = resolveSelectedContrastResultsFn,
    catFn = catFn
  )

  persistenceResult <- persistAnalysisResultsFn(
    input = input,
    enrichmentData = enrichmentData,
    workflowData = workflowData,
    selectedContrast = selectedContrast,
    enrichmentResults = enrichmentResults,
    methodInfo = methodInfo,
    pathwayDir = pathwayDir,
    buildAnalysisResultsPayloadFn = buildAnalysisResultsPayloadFn,
    propagateResultsArgsFn = propagateResultsArgsFn,
    propagateUiParamsFn = propagateUiParamsFn,
    updateStateManagerUiParamsFn = updateStateManagerUiParamsFn,
    saveCompletedStateFn = saveCompletedStateFn,
    completeTabStatusFn = completeTabStatusFn,
    completeProgressFn = completeProgressFn,
    incProgressFn = incProgressFn
  )
  enrichmentResults <- persistenceResult$enrichmentResults

  list(
    selectedContrast = selectedContrast,
    rawContrastName = rawContrastName,
    analysisMethod = methodInfo$method,
    analysisComplete = persistenceResult$analysisComplete,
    organismFilterApplied = organismFilterApplied,
    filterStats = filterStats,
    enrichmentResults = enrichmentResults
  )
}

resolveProtEnrichAnalysisInputColumns <- function(methodInfo,
                                                  daResultsForEnrichment,
                                                  currentS4Object = NULL) {
  idColumn <- tryCatch({
    if (!is.null(currentS4Object)) {
      currentS4Object@protein_id_column
    } else {
      "uniprot_acc"
    }
  }, error = function(e) "uniprot_acc")
  sourceLabel <- if (!is.null(currentS4Object)) "s4_object" else "default"
  geneNameOverrideApplied <- FALSE

  firstDaData <- tryCatch(daResultsForEnrichment@da_data[[1]], error = function(e) NULL)
  if (!is.null(methodInfo) &&
      identical(methodInfo$method, "gprofiler2") &&
      !is.null(firstDaData) &&
      "gene_name" %in% names(firstDaData)) {
    idColumn <- "gene_name"
    sourceLabel <- "gprofiler_gene_name_override"
    geneNameOverrideApplied <- TRUE
  }

  list(
    idColumn = idColumn,
    source = sourceLabel,
    geneNameOverrideApplied = geneNameOverrideApplied
  )
}

buildProtEnrichProcessEnrichmentsArgs <- function(daResultsForEnrichment,
                                                  organismTaxid,
                                                  upCutoff,
                                                  downCutoff,
                                                  qCutoff,
                                                  pathwayDir,
                                                  goAnnotations,
                                                  proteinIdColumn,
                                                  contrastNames,
                                                  correctionMethod,
                                                  excludeIea = FALSE) {
  taxonId <- as.numeric(organismTaxid)

  list(
    checkpointArgs = list(
      da_results_s4 = daResultsForEnrichment,
      taxon_id = taxonId,
      up_cutoff = upCutoff,
      down_cutoff = downCutoff,
      q_cutoff = qCutoff,
      pathway_dir = pathwayDir,
      go_annotations = goAnnotations,
      exclude_iea = excludeIea,
      protein_id_column = proteinIdColumn,
      contrast_names = contrastNames,
      correction_method = correctionMethod
    ),
    processArgs = list(
      da_results = daResultsForEnrichment,
      taxon_id = taxonId,
      up_cutoff = upCutoff,
      down_cutoff = downCutoff,
      q_cutoff = qCutoff,
      pathway_dir = pathwayDir,
      go_annotations = goAnnotations,
      exclude_iea = excludeIea,
      protein_id_column = proteinIdColumn,
      contrast_names = contrastNames,
      correction_method = correctionMethod
    )
  )
}

prepareProtEnrichProcessExecution <- function(input,
                                              enrichmentData,
                                              daResultsForEnrichment,
                                              pathwayDir,
                                              goAnnotations,
                                              currentAnalysisMethodFn,
                                              resolveAnalysisInputColumnsFn = resolveProtEnrichAnalysisInputColumns,
                                              buildProcessEnrichmentsArgsFn = buildProtEnrichProcessEnrichmentsArgs,
                                              catFn = cat) {
  methodInfo <- currentAnalysisMethodFn()
  catFn(sprintf("   ENRICHMENT Step: Using analysis method: %s\n", methodInfo$method))

  inputColumnConfig <- resolveAnalysisInputColumnsFn(
    methodInfo = methodInfo,
    daResultsForEnrichment = daResultsForEnrichment,
    currentS4Object = enrichmentData$current_s4_object
  )

  enrichmentArgs <- buildProcessEnrichmentsArgsFn(
    daResultsForEnrichment = daResultsForEnrichment,
    organismTaxid = input$organism_taxid,
    upCutoff = input$up_cutoff,
    downCutoff = input$down_cutoff,
    qCutoff = input$q_cutoff,
    pathwayDir = pathwayDir,
    goAnnotations = goAnnotations,
    proteinIdColumn = inputColumnConfig$idColumn,
    contrastNames = names(enrichmentData$da_results_data),
    correctionMethod = input$correction_method
  )

  list(
    methodInfo = methodInfo,
    inputColumnConfig = inputColumnConfig,
    enrichmentArgs = enrichmentArgs
  )
}

executeProtEnrichProcessEnrichments <- function(enrichmentArgs,
                                                upCutoff,
                                                downCutoff,
                                                qCutoff,
                                                captureCheckpointFn = .capture_checkpoint,
                                                processEnrichmentsFn = processEnrichments,
                                                catFn = cat) {
  captureCheckpointFn(enrichmentArgs$checkpointArgs, "cp10", "enrichment_input")

  enrichmentResults <- do.call(processEnrichmentsFn, enrichmentArgs$processArgs)

  catFn(sprintf(
    "   ENRICHMENT Step: processEnrichments completed with up_cutoff: %.1f, down_cutoff: %.1f, q_cutoff: %.3f\n",
    upCutoff,
    downCutoff,
    qCutoff
  ))

  enrichmentResults
}

buildProtEnrichAllContrastResults <- function(enrichmentResults,
                                              methodInfo,
                                              isEnrichResultFn = methods::is,
                                              catFn = cat) {
  catFn("   ENRICHMENT Step: Processing results for ALL contrasts\n")

  allContrastResults <- list()

  for (rawContrastName in names(enrichmentResults@enrichment_data)) {
    catFn(sprintf("   ENRICHMENT Step: Processing contrast '%s'\n", rawContrastName))

    contrastEnrichment <- enrichmentResults@enrichment_data[[rawContrastName]]
    contrastResults <- list(
      gprofiler_results = NULL,
      clusterprofiler_results = NULL,
      stringdb_results = NULL
    )

    if (identical(methodInfo$method, "gprofiler2")) {
      gprofilerResults <- data.frame()

      if (!is.null(contrastEnrichment$up) || !is.null(contrastEnrichment$down)) {
        if (!is.null(contrastEnrichment$up)) {
          upResults <- contrastEnrichment$up
          if (!is.null(upResults) &&
              !is.null(upResults$result) &&
              length(upResults$result) > 0 &&
              nrow(upResults$result) > 0) {
            upDf <- upResults$result
            upDf$directionality <- "positive"
            upDf$analysis_method <- "gprofiler2"
            gprofilerResults <- rbind(gprofilerResults, upDf)
          }
        }

        if (!is.null(contrastEnrichment$down)) {
          downResults <- contrastEnrichment$down
          if (!is.null(downResults) &&
              !is.null(downResults$result) &&
              length(downResults$result) > 0 &&
              nrow(downResults$result) > 0) {
            downDf <- downResults$result
            downDf$directionality <- "negative"
            downDf$analysis_method <- "gprofiler2"
            gprofilerResults <- rbind(gprofilerResults, downDf)
          }
        }

        if (nrow(gprofilerResults) > 0) {
          if ("term_name" %in% names(gprofilerResults)) {
            gprofilerResults$Description <- gprofilerResults$term_name
          }
          if ("p_value" %in% names(gprofilerResults)) {
            gprofilerResults$pvalue <- gprofilerResults$p_value
            gprofilerResults$p.adjust <- gprofilerResults$p_value
            gprofilerResults$qvalue <- gprofilerResults$p_value
          }
          if ("term_size" %in% names(gprofilerResults)) {
            gprofilerResults$Count <- gprofilerResults$term_size
          }
          if ("source" %in% names(gprofilerResults)) {
            gprofilerResults$data_source <- gprofilerResults$source
          }
        }
      }

      contrastResults$gprofiler_results <- gprofilerResults
      catFn(sprintf(
        "   ENRICHMENT Step: Contrast '%s' - %d gprofiler2 results\n",
        rawContrastName,
        nrow(gprofilerResults)
      ))
    } else if (identical(methodInfo$method, "clusterprofiler")) {
      clusterprofilerResults <- data.frame()

      if (!is.null(contrastEnrichment$up) || !is.null(contrastEnrichment$down)) {
        if (!is.null(contrastEnrichment$up)) {
          upResults <- contrastEnrichment$up
          if (!is.null(upResults) &&
              isTRUE(isEnrichResultFn(upResults, "enrichResult")) &&
              nrow(upResults@result) > 0) {
            upDf <- upResults@result
            upDf$directionality <- "up"
            upDf$analysis_method <- "clusterprofiler"
            clusterprofilerResults <- rbind(clusterprofilerResults, upDf)
          }
        }

        if (!is.null(contrastEnrichment$down)) {
          downResults <- contrastEnrichment$down
          if (!is.null(downResults) &&
              isTRUE(isEnrichResultFn(downResults, "enrichResult")) &&
              nrow(downResults@result) > 0) {
            downDf <- downResults@result
            downDf$directionality <- "down"
            downDf$analysis_method <- "clusterprofiler"
            clusterprofilerResults <- rbind(clusterprofilerResults, downDf)
          }
        }
      }

      contrastResults$clusterprofiler_results <- clusterprofilerResults
      catFn(sprintf(
        "   ENRICHMENT Step: Contrast '%s' - %d clusterProfileR results\n",
        rawContrastName,
        nrow(clusterprofilerResults)
      ))
    }

    allContrastResults[[rawContrastName]] <- contrastResults
  }

  allContrastResults
}

prepareProtEnrichAnalysisBodySetup <- function(selectedContrast,
                                               input,
                                               enrichmentData,
                                               workflowData,
                                               experimentPaths,
                                               resolveSelectedDaResultsFn = resolveProtEnrichSelectedDaResults,
                                               resolveRunDependenciesFn = resolveProtEnrichRunDependencies,
                                               resolveOutputDirectoriesFn = resolveProtEnrichOutputDirectories,
                                               createDaResultsForEnrichmentFn = createDAResultsForEnrichment,
                                               resolveUniprotAnnotationsFn = resolveProtEnrichUniprotAnnotations,
                                               resolveAnnotationMatchingFn = resolveProtEnrichAnnotationMatching,
                                               resolveOrganismMappingFn = resolveProtEnrichOrganismMapping,
                                               applyOrganismFilterFn = applyProtEnrichOrganismFilter,
                                               persistOrganismFilterMetadataFn = persistProtEnrichOrganismFilterMetadata,
                                               showNotificationFn = shiny::showNotification,
                                               globalEnv = .GlobalEnv,
                                               catFn = cat) {
  catFn(sprintf("   ENRICHMENT Step: Selected contrast (friendly name) = %s\n", selectedContrast))
  catFn(sprintf(
    "   ENRICHMENT Step: Available DE results: %s\n",
    paste(names(enrichmentData$da_results_data), collapse = ", ")
  ))

  contrastsTbl <- if (exists("contrasts_tbl", envir = globalEnv)) {
    get("contrasts_tbl", envir = globalEnv)
  } else {
    NULL
  }

  selectedDaConfig <- resolveSelectedDaResultsFn(
    selectedContrast = selectedContrast,
    daResultsData = enrichmentData$da_results_data,
    contrastsTbl = contrastsTbl
  )
  selectedDaResults <- selectedDaConfig$selectedDaResults
  rawContrastName <- selectedDaConfig$rawContrastName
  availableKeys <- selectedDaConfig$availableKeys

  if (!is.null(selectedDaConfig$mappedRawContrastName)) {
    catFn(sprintf(
      "   ENRICHMENT Step: Mapped friendly name '%s' to raw name '%s'\n",
      selectedContrast,
      selectedDaConfig$mappedRawContrastName
    ))
  }

  if (identical(selectedDaConfig$source, "friendly_name")) {
    catFn(sprintf("   ENRICHMENT Step: Found DE results for raw contrast %s\n", rawContrastName))
  } else if (!is.null(selectedDaConfig$mappedRawContrastName)) {
    catFn(sprintf(
      "   ENRICHMENT Step: No DE results found for raw contrast '%s'\n",
      selectedDaConfig$mappedRawContrastName
    ))
  }

  if (is.null(selectedDaResults)) {
    catFn(sprintf(
      "   ENRICHMENT Step: Available DE result keys: %s\n",
      paste(availableKeys, collapse = ", ")
    ))
    catFn(sprintf("   ENRICHMENT Step: Looking for friendly name: %s\n", selectedContrast))
  }

  if (identical(selectedDaConfig$source, "fuzzy_match")) {
    catFn(sprintf(
      "   ENRICHMENT Step: Found matching DE results using fuzzy matching: %s\n",
      rawContrastName
    ))
  }

  if (identical(selectedDaConfig$source, "direct_key")) {
    catFn(sprintf(
      "   ENRICHMENT Step: Found DE results using direct key lookup: %s\n",
      selectedContrast
    ))
  }

  if (is.null(selectedDaResults)) {
    stop(sprintf(
      "Could not find DE results for contrast '%s'. Available contrasts: %s",
      selectedContrast,
      paste(availableKeys, collapse = ", ")
    ))
  }

  hadCurrentS4 <- !is.null(enrichmentData$current_s4_object)
  dependencyConfig <- resolveRunDependenciesFn(
    currentS4Object = enrichmentData$current_s4_object,
    daResultsData = enrichmentData$da_results_data,
    workflowData = workflowData,
    contrastsTbl = contrastsTbl
  )
  contrastsTbl <- dependencyConfig$contrastsTbl
  designMatrix <- dependencyConfig$designMatrix
  enrichmentData$current_s4_object <- dependencyConfig$currentS4

  if (!is.null(contrastsTbl)) {
    catFn("   ENRICHMENT Step: Found contrasts_tbl in global environment\n")
  }

  if (!hadCurrentS4 && !is.null(enrichmentData$current_s4_object)) {
    catFn("   ENRICHMENT Step: S4 object is NULL, trying to retrieve it...\n")

    if (identical(dependencyConfig$s4Source, "da_results_first_result")) {
      catFn(sprintf(
        "   ENRICHMENT Step: Got S4 from DE results theObject (class: %s)\n",
        class(enrichmentData$current_s4_object)
      ))
    }

    if (identical(dependencyConfig$s4Source, "state_manager")) {
      catFn(sprintf(
        "   ENRICHMENT Step: Got S4 from state manager (class: %s)\n",
        class(enrichmentData$current_s4_object)
      ))
    }
  }

  if (!is.null(dependencyConfig$designMatrixError)) {
    catFn(sprintf(
      "   ENRICHMENT Step: Could not access design_matrix slot: %s\n",
      dependencyConfig$designMatrixError
    ))
  }

  if (identical(dependencyConfig$designMatrixSource, "s4_object")) {
    catFn("   ENRICHMENT Step: Found design_matrix in S4 object\n")
  }

  if (identical(dependencyConfig$designMatrixSource, "global_environment")) {
    catFn("   ENRICHMENT Step: Found design_matrix in global environment\n")
  }

  if (is.null(contrastsTbl) || is.null(designMatrix)) {
    stop("Missing contrasts_tbl or design_matrix required for enrichment analysis. Please ensure DE analysis was completed successfully.")
  }

  pathConfig <- resolveOutputDirectoriesFn(
    experimentPaths = experimentPaths
  )
  daOutputDir <- pathConfig$daOutputDir
  pathwayDir <- pathConfig$pathwayDir

  if (identical(pathConfig$daOutputDirSource, "experiment_paths")) {
    catFn(sprintf("   ENRICHMENT Step: Using da_output_dir from experiment_paths: %s\n", daOutputDir))
  } else {
    catFn(sprintf("   ENRICHMENT Step: Falling back to da_proteins folder: %s\n", daOutputDir))
  }

  if (identical(pathConfig$pathwayDirSource, "experiment_paths")) {
    catFn(sprintf("   ENRICHMENT Step: Using pathway_dir from experiment_paths: %s\n", pathwayDir))
  } else {
    catFn(sprintf("   ENRICHMENT Step: Created fallback pathway directory: %s\n", pathwayDir))
  }

  daResultsForEnrichment <- createDaResultsForEnrichmentFn(
    contrasts_tbl = contrastsTbl,
    design_matrix = designMatrix,
    da_output_dir = daOutputDir
  )

  catFn("   ENRICHMENT Step: S4 object created successfully\n")

  uniprotConfig <- resolveUniprotAnnotationsFn(
    workflowData = workflowData,
    experimentPaths = experimentPaths,
    currentS4Object = enrichmentData$current_s4_object,
    organismTaxid = input$organism_taxid
  )
  uniprotDatCln <- uniprotConfig$uniprotDatCln

  annotationMatchConfig <- resolveAnnotationMatchingFn(
    uniprotDatCln = uniprotDatCln,
    daResultsForEnrichment = daResultsForEnrichment,
    currentS4Object = enrichmentData$current_s4_object
  )

  if (!is.null(annotationMatchConfig$annotationMatchResults)) {
    enrichmentData$annotation_match_results <- annotationMatchConfig$annotationMatchResults
  }

  organismFilterApplied <- FALSE
  filterStats <- list(proteins_before = 0, proteins_after = 0, proteins_removed = 0)

  if (isTRUE(input$enable_organism_filter)) {
    catFn("*** ENRICHMENT Step: Multi-species filtering ENABLED ***\n")

    targetTaxon <- as.character(input$organism_taxid)
    organismMappingConfig <- resolveOrganismMappingFn(
      workflowData = workflowData,
      uniprotDatCln = uniprotDatCln,
      targetTaxon = targetTaxon
    )
    organismMapping <- organismMappingConfig$organismMapping

    if (!is.null(organismMapping) && nrow(organismMapping) > 0) {
      organismFilterConfig <- applyOrganismFilterFn(
        daResultsForEnrichment = daResultsForEnrichment,
        organismMapping = organismMapping,
        targetTaxon = targetTaxon,
        currentS4Object = enrichmentData$current_s4_object
      )
      daResultsForEnrichment <- organismFilterConfig$daResultsForEnrichment
      organismFilterApplied <- organismFilterConfig$filterApplied
      filterStats <- organismFilterConfig$filterStats
    } else {
      catFn("   ENRICHMENT Step: WARNING - No organism mapping available, filtering skipped\n")
      showNotificationFn(
        "Multi-species filtering enabled but no organism mapping available. Proceeding without filtering.",
        type = "warning",
        duration = 8
      )
    }
  } else {
    catFn("   ENRICHMENT Step: Multi-species filtering DISABLED\n")
  }

  persistOrganismFilterMetadataFn(
    workflowData = workflowData,
    organismFilterEnabled = input$enable_organism_filter,
    organismFilterApplied = organismFilterApplied,
    targetTaxonId = input$organism_taxid,
    filterStats = filterStats
  )

  list(
    rawContrastName = rawContrastName,
    contrastsTbl = contrastsTbl,
    pathwayDir = pathwayDir,
    daResultsForEnrichment = daResultsForEnrichment,
    goAnnotations = uniprotDatCln,
    organismFilterApplied = organismFilterApplied,
    filterStats = filterStats
  )
}

runProtEnrichAnalysisBody <- function(input,
                                      enrichmentData,
                                      workflowData,
                                      experimentPaths,
                                      currentAnalysisMethodFn,
                                      prepareAnalysisSetupFn = prepareProtEnrichAnalysisBodySetup,
                                      captureCheckpointFn = .capture_checkpoint,
                                      resolveSelectedDaResultsFn = resolveProtEnrichSelectedDaResults,
                                      resolveRunDependenciesFn = resolveProtEnrichRunDependencies,
                                      resolveOutputDirectoriesFn = resolveProtEnrichOutputDirectories,
                                      createDaResultsForEnrichmentFn = createDAResultsForEnrichment,
                                      resolveUniprotAnnotationsFn = resolveProtEnrichUniprotAnnotations,
                                      resolveAnnotationMatchingFn = resolveProtEnrichAnnotationMatching,
                                      resolveOrganismMappingFn = resolveProtEnrichOrganismMapping,
                                      applyOrganismFilterFn = applyProtEnrichOrganismFilter,
                                      persistOrganismFilterMetadataFn = persistProtEnrichOrganismFilterMetadata,
                                      resolveAnalysisInputColumnsFn = resolveProtEnrichAnalysisInputColumns,
                                      buildProcessEnrichmentsArgsFn = buildProtEnrichProcessEnrichmentsArgs,
                                      prepareProcessExecutionFn = prepareProtEnrichProcessExecution,
                                      processEnrichmentsFn = processEnrichments,
                                      executeProcessEnrichmentsFn = executeProtEnrichProcessEnrichments,
                                      buildAllContrastResultsFn = buildProtEnrichAllContrastResults,
                                      resolveSelectedContrastResultsFn = resolveProtEnrichSelectedContrastResults,
                                      capturePostProcessResultsFn = captureProtEnrichPostProcessResults,
                                      persistAnalysisResultsFn = persistProtEnrichAnalysisResults,
                                      buildAnalysisResultsPayloadFn = buildProtEnrichAnalysisResultsPayload,
                                      propagateResultsArgsFn = propagateProtEnrichResultsArgs,
                                      propagateUiParamsFn = propagateProtEnrichUiParams,
                                      updateStateManagerUiParamsFn = updateProtEnrichStateManagerUiParams,
                                      saveCompletedStateFn = saveProtEnrichCompletedState,
                                      completeTabStatusFn = completeProtEnrichTabStatus,
                                      completeProgressFn = completeProtEnrichProgress,
                                      finalizeAnalysisResultsFn = finalizeProtEnrichAnalysisBodyResults,
                                      showNotificationFn = shiny::showNotification,
                                      incProgressFn = shiny::incProgress,
                                      globalEnv = .GlobalEnv,
                                      catFn = cat) {
  selectedContrast <- input$selected_contrast

  incProgressFn(0.2, detail = "Transforming DE data...")

  incProgressFn(0.3, detail = "Creating DE results S4 object...")

  catFn("   ENRICHMENT Step: Creating DA results S4 object using createDAResultsForEnrichment\n")

  setupConfig <- prepareAnalysisSetupFn(
    selectedContrast = selectedContrast,
    input = input,
    enrichmentData = enrichmentData,
    workflowData = workflowData,
    experimentPaths = experimentPaths,
    resolveSelectedDaResultsFn = resolveSelectedDaResultsFn,
    resolveRunDependenciesFn = resolveRunDependenciesFn,
    resolveOutputDirectoriesFn = resolveOutputDirectoriesFn,
    createDaResultsForEnrichmentFn = createDaResultsForEnrichmentFn,
    resolveUniprotAnnotationsFn = resolveUniprotAnnotationsFn,
    resolveAnnotationMatchingFn = resolveAnnotationMatchingFn,
    resolveOrganismMappingFn = resolveOrganismMappingFn,
    applyOrganismFilterFn = applyOrganismFilterFn,
    persistOrganismFilterMetadataFn = persistOrganismFilterMetadataFn,
    showNotificationFn = showNotificationFn,
    globalEnv = globalEnv,
    catFn = catFn
  )
  rawContrastName <- setupConfig$rawContrastName
  contrastsTbl <- setupConfig$contrastsTbl
  pathwayDir <- setupConfig$pathwayDir
  daResultsForEnrichment <- setupConfig$daResultsForEnrichment
  organismFilterApplied <- setupConfig$organismFilterApplied
  filterStats <- setupConfig$filterStats

  incProgressFn(0.5, detail = "Running enrichment analysis...")

  catFn("   ENRICHMENT Step: Running processEnrichments\n")

  processExecutionConfig <- prepareProcessExecutionFn(
    input = input,
    enrichmentData = enrichmentData,
    daResultsForEnrichment = daResultsForEnrichment,
    pathwayDir = pathwayDir,
    goAnnotations = setupConfig$goAnnotations,
    currentAnalysisMethodFn = currentAnalysisMethodFn,
    resolveAnalysisInputColumnsFn = resolveAnalysisInputColumnsFn,
    buildProcessEnrichmentsArgsFn = buildProcessEnrichmentsArgsFn,
    catFn = catFn
  )
  methodInfo <- processExecutionConfig$methodInfo
  enrichmentArgs <- processExecutionConfig$enrichmentArgs

  enrichmentResults <- executeProcessEnrichmentsFn(
    enrichmentArgs = enrichmentArgs,
    upCutoff = input$up_cutoff,
    downCutoff = input$down_cutoff,
    qCutoff = input$q_cutoff,
    captureCheckpointFn = captureCheckpointFn,
    processEnrichmentsFn = processEnrichmentsFn,
    catFn = catFn
  )

  finalizeAnalysisResultsFn(
    selectedContrast = selectedContrast,
    rawContrastName = rawContrastName,
    organismFilterApplied = organismFilterApplied,
    filterStats = filterStats,
    enrichmentResults = enrichmentResults,
    enrichmentData = enrichmentData,
    workflowData = workflowData,
    input = input,
    methodInfo = methodInfo,
    contrastsTbl = contrastsTbl,
    pathwayDir = pathwayDir,
    buildAllContrastResultsFn = buildAllContrastResultsFn,
    resolveSelectedContrastResultsFn = resolveSelectedContrastResultsFn,
    capturePostProcessResultsFn = capturePostProcessResultsFn,
    persistAnalysisResultsFn = persistAnalysisResultsFn,
    buildAnalysisResultsPayloadFn = buildAnalysisResultsPayloadFn,
    propagateResultsArgsFn = propagateResultsArgsFn,
    propagateUiParamsFn = propagateUiParamsFn,
    updateStateManagerUiParamsFn = updateStateManagerUiParamsFn,
    saveCompletedStateFn = saveCompletedStateFn,
    completeTabStatusFn = completeTabStatusFn,
    completeProgressFn = completeProgressFn,
    incProgressFn = incProgressFn,
    catFn = catFn
  )
}
resolveProtEnrichOrganismMapping <- function(workflowData,
                                             uniprotDatCln,
                                             targetTaxon,
                                             catFn = cat) {
  organismMapping <- NULL
  mappingSource <- NULL
  accessionColumn <- NULL
  taxonColumn <- NULL
  warningMessage <- NULL
  availableTaxonIds <- character()

  if (!is.null(workflowData$mixed_species_analysis) &&
      !is.null(workflowData$mixed_species_analysis$organism_mapping)) {
    organismMapping <- workflowData$mixed_species_analysis$organism_mapping
    mappingSource <- "workflow_data"
    catFn("   ENRICHMENT Step: Using organism_mapping from import module\n")
  }

  if (is.null(organismMapping) && !is.null(uniprotDatCln)) {
    catFn(sprintf(
      "   ENRICHMENT Step: Checking uniprot_dat_cln columns: %s\n",
      paste(names(uniprotDatCln), collapse = ", ")
    ))

    organismNameToTaxid <- c(
      "Homo sapiens" = "9606", "Human" = "9606",
      "Mus musculus" = "10090", "Mouse" = "10090",
      "Rattus norvegicus" = "10116", "Rat" = "10116",
      "Drosophila melanogaster" = "7227", "Fruit fly" = "7227",
      "Caenorhabditis elegans" = "6239",
      "Saccharomyces cerevisiae" = "4932", "Yeast" = "4932",
      "Arabidopsis thaliana" = "3702",
      "Danio rerio" = "7955", "Zebrafish" = "7955",
      "Gallus gallus" = "9031", "Chicken" = "9031",
      "Sus scrofa" = "9823", "Pig" = "9823",
      "Bos taurus" = "9913", "Bovine" = "9913", "Cow" = "9913"
    )
    possibleAccCols <- c(
      "Entry", "entry", "UniProt_Acc", "uniprot_acc", "Accession", "accession", "protein_id"
    )

    for (col in possibleAccCols) {
      if (col %in% names(uniprotDatCln)) {
        accessionColumn <- col
        break
      }
    }

    if ("Organism" %in% names(uniprotDatCln) && !is.null(accessionColumn)) {
      catFn("   ENRICHMENT Step: Found 'Organism' column with organism names - mapping to taxon IDs\n")

      mapOrganismToTaxid <- function(orgName) {
        if (is.na(orgName) || orgName == "") {
          return(NA_character_)
        }

        for (name in names(organismNameToTaxid)) {
          if (grepl(name, orgName, ignore.case = TRUE)) {
            return(unname(organismNameToTaxid[[name]]))
          }
        }

        NA_character_
      }

      organismMapping <- tryCatch({
        uniprotDatCln |>
          dplyr::select(uniprot_acc = dplyr::all_of(accessionColumn), organism_name = Organism) |>
          dplyr::mutate(
            taxon_id = unname(vapply(organism_name, mapOrganismToTaxid, character(1)))
          ) |>
          dplyr::select(uniprot_acc, taxon_id)
      }, error = function(e) {
        warningMessage <<- e$message
        catFn(sprintf("   ENRICHMENT Step: Error creating organism_mapping from names: %s\n", e$message))
        NULL
      })

      if (!is.null(organismMapping)) {
        mappingSource <- "organism_names"
        availableTaxonIds <- unique(organismMapping$taxon_id)
        catFn(sprintf(
          "   ENRICHMENT Step: Created organism_mapping from organism names (%d entries)\n",
          nrow(organismMapping)
        ))
        catFn(sprintf(
          "   ENRICHMENT Step: Unique taxon IDs mapped: %s\n",
          paste(availableTaxonIds[!is.na(availableTaxonIds)], collapse = ", ")
        ))

        if (targetTaxon %in% availableTaxonIds) {
          targetCount <- sum(organismMapping$taxon_id == targetTaxon, na.rm = TRUE)
          catFn(sprintf(
            "   ENRICHMENT Step: Found %d proteins matching target taxon %s\n",
            targetCount,
            targetTaxon
          ))
        } else {
          catFn(sprintf(
            "   ENRICHMENT Step: WARNING - Target taxon %s not found in mapped taxon IDs!\n",
            targetTaxon
          ))
          catFn(sprintf(
            "   ENRICHMENT Step: Available taxon IDs: %s\n",
            paste(availableTaxonIds[!is.na(availableTaxonIds)], collapse = ", ")
          ))
        }
      }
    } else {
      possibleTaxonCols <- c(
        "Organism (ID)", "organism_id", "Organism_ID", "taxon_id", "Taxon_ID",
        "Taxonomy ID", "taxonomy_id", "NCBI_TaxID", "ncbi_taxid", "OX"
      )

      for (col in possibleTaxonCols) {
        if (col %in% names(uniprotDatCln)) {
          taxonColumn <- col
          catFn(sprintf("   ENRICHMENT Step: Found taxon column: %s\n", taxonColumn))
          break
        }
      }

      if (!is.null(taxonColumn) && !is.null(accessionColumn)) {
        organismMapping <- tryCatch({
          uniprotDatCln |>
            dplyr::select(
              uniprot_acc = dplyr::all_of(accessionColumn),
              taxon_raw = dplyr::all_of(taxonColumn)
            ) |>
            dplyr::mutate(
              taxon_id = dplyr::case_when(
                grepl("ID=", taxon_raw) ~ stringr::str_extract(taxon_raw, "(?<=ID=)\\d+"),
                grepl("^\\d+$", as.character(taxon_raw)) ~ as.character(taxon_raw),
                TRUE ~ as.character(taxon_raw)
              )
            ) |>
            dplyr::select(uniprot_acc, taxon_id)
        }, error = function(e) {
          warningMessage <<- e$message
          catFn(sprintf("   ENRICHMENT Step: Error creating organism_mapping: %s\n", e$message))
          NULL
        })

        if (!is.null(organismMapping)) {
          mappingSource <- "taxon_column"
          catFn(sprintf(
            "   ENRICHMENT Step: Created organism_mapping from taxon column (%d entries)\n",
            nrow(organismMapping)
          ))
        }
      } else if (!is.null(accessionColumn)) {
        catFn("   ENRICHMENT Step: No organism column found - creating single-species mapping\n")
        organismMapping <- uniprotDatCln |>
          dplyr::select(uniprot_acc = dplyr::all_of(accessionColumn)) |>
          dplyr::mutate(taxon_id = targetTaxon)
        mappingSource <- "single_species_fallback"
        catFn(sprintf(
          "   ENRICHMENT Step: Created single-species organism_mapping (%d entries, all assigned to taxon %s)\n",
          nrow(organismMapping),
          targetTaxon
        ))
      }
    }
  }

  if (length(availableTaxonIds) == 0 &&
      !is.null(organismMapping) &&
      "taxon_id" %in% names(organismMapping)) {
    availableTaxonIds <- unique(organismMapping$taxon_id)
  }

  list(
    organismMapping = organismMapping,
    source = mappingSource,
    accessionColumn = accessionColumn,
    taxonColumn = taxonColumn,
    availableTaxonIds = availableTaxonIds[!is.na(availableTaxonIds)],
    warning = warningMessage
  )
}

applyProtEnrichOrganismFilter <- function(daResultsForEnrichment,
                                          organismMapping,
                                          targetTaxon,
                                          currentS4Object = NULL,
                                          normalizeUniprotFn = normalizeUniprotAccession,
                                          cleanAccFn = clean_acc,
                                          catFn = cat) {
  filterApplied <- FALSE
  proteinsBefore <- 0
  proteinsAfter <- 0
  proteinsRemoved <- 0
  proteinIdCol <- tryCatch({
    if (!is.null(currentS4Object)) {
      currentS4Object@protein_id_column
    } else {
      "uniprot_acc"
    }
  }, error = function(e) "uniprot_acc")

  targetProteins <- organismMapping |>
    dplyr::filter(taxon_id == targetTaxon) |>
    dplyr::pull(uniprot_acc) |>
    unique()
  targetProteinsClean <- unique(normalizeUniprotFn(targetProteins, remove_isoform = TRUE))

  catFn(sprintf(
    "   ENRICHMENT Step: Found %d proteins for target taxon %s\n",
    length(targetProteins),
    targetTaxon
  ))

  if (!is.null(daResultsForEnrichment@da_data) &&
      length(daResultsForEnrichment@da_data) > 0) {
    filteredDeData <- lapply(names(daResultsForEnrichment@da_data), function(contrastName) {
      contrastData <- daResultsForEnrichment@da_data[[contrastName]]

      if (!is.null(contrastData) && proteinIdCol %in% names(contrastData)) {
        originalCount <- nrow(contrastData)

        keepRows <- vapply(contrastData[[proteinIdCol]], function(proteinIds) {
          ids <- unlist(strsplit(as.character(proteinIds), ";"))
          idsClean <- cleanAccFn(trimws(ids))
          any(idsClean %in% targetProteinsClean) || any(ids %in% targetProteins)
        }, logical(1))

        filteredData <- contrastData[keepRows, , drop = FALSE]
        filteredCount <- nrow(filteredData)

        catFn(sprintf(
          "   ENRICHMENT Step: Contrast '%s': %d -> %d proteins (removed %d non-target organism)\n",
          contrastName,
          originalCount,
          filteredCount,
          originalCount - filteredCount
        ))

        proteinsBefore <<- proteinsBefore + originalCount
        proteinsAfter <<- proteinsAfter + filteredCount
        proteinsRemoved <<- proteinsRemoved + (originalCount - filteredCount)

        return(filteredData)
      }

      contrastData
    })
    names(filteredDeData) <- names(daResultsForEnrichment@da_data)

    daResultsForEnrichment@da_data <- filteredDeData
    filterApplied <- TRUE

    catFn(sprintf(
      "*** ENRICHMENT Step: Organism filtering complete - kept %d/%d proteins (%.1f%%) ***\n",
      proteinsAfter,
      proteinsBefore,
      (proteinsAfter / max(proteinsBefore, 1)) * 100
    ))
  }

  filterStats <- list(
    proteins_before = proteinsBefore,
    proteins_after = proteinsAfter,
    proteins_removed = proteinsRemoved
  )

  list(
    daResultsForEnrichment = daResultsForEnrichment,
    filterApplied = filterApplied,
    filterStats = filterStats,
    proteinIdColumn = proteinIdCol,
    targetProteins = targetProteins,
    targetProteinsClean = targetProteinsClean
  )
}

persistProtEnrichOrganismFilterMetadata <- function(workflowData,
                                                    organismFilterEnabled,
                                                    organismFilterApplied,
                                                    targetTaxonId,
                                                    filterStats,
                                                    timeFn = Sys.time) {
  metadata <- list(
    enabled = isTRUE(organismFilterEnabled),
    filter_applied = organismFilterApplied,
    target_taxon_id = targetTaxonId,
    proteins_before = filterStats$proteins_before,
    proteins_after = filterStats$proteins_after,
    proteins_removed = filterStats$proteins_removed,
    timestamp = timeFn()
  )

  workflowData$enrichment_organism_filter <- metadata

  metadata
}

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

resolveProtEnrichRunDependencies <- function(currentS4Object,
                                             daResultsData,
                                             workflowData,
                                             contrastsTbl = NULL,
                                             globalEnv = .GlobalEnv) {
  resolvedCurrentS4 <- currentS4Object
  s4Source <- if (!is.null(resolvedCurrentS4)) "current_s4_object" else NULL

  if (is.null(resolvedCurrentS4) &&
      !is.null(daResultsData) &&
      length(daResultsData) > 0) {
    firstResult <- daResultsData[[1]]
    firstObject <- tryCatch(firstResult$theObject, error = function(e) NULL)
    if (!is.null(firstObject)) {
      resolvedCurrentS4 <- firstObject
      s4Source <- "da_results_first_result"
    }
  }

  if (is.null(resolvedCurrentS4) && !is.null(workflowData$state_manager)) {
    currentState <- workflowData$state_manager$current_state
    stateObject <- workflowData$state_manager$getState(currentState)
    if (!is.null(stateObject)) {
      resolvedCurrentS4 <- stateObject
      s4Source <- "state_manager"
    }
  }

  designMatrix <- NULL
  designMatrixSource <- NULL
  designMatrixError <- NULL

  if (!is.null(resolvedCurrentS4)) {
    tryCatch({
      if (!is.null(resolvedCurrentS4@design_matrix)) {
        designMatrix <- resolvedCurrentS4@design_matrix
        designMatrixSource <- "s4_object"
      }
    }, error = function(e) {
      designMatrixError <<- e$message
    })
  }

  if (is.null(designMatrix) && exists("design_matrix", envir = globalEnv)) {
    designMatrix <- get("design_matrix", envir = globalEnv)
    designMatrixSource <- "global_environment"
  }

  list(
    contrastsTbl = contrastsTbl,
    currentS4 = resolvedCurrentS4,
    s4Source = s4Source,
    designMatrix = designMatrix,
    designMatrixSource = designMatrixSource,
    designMatrixError = designMatrixError
  )
}

resolveProtEnrichOutputDirectories <- function(experimentPaths,
                                               dirExistsFn = dir.exists,
                                               dirCreateFn = dir.create,
                                               filePathFn = file.path) {
  daOutputDir <- NULL
  daOutputDirSource <- NULL
  daOutputDirCreated <- FALSE

  if (!is.null(experimentPaths$da_output_dir) &&
      dirExistsFn(experimentPaths$da_output_dir)) {
    daOutputDir <- experimentPaths$da_output_dir
    daOutputDirSource <- "experiment_paths"
  } else {
    daOutputDir <- filePathFn(experimentPaths$results_dir, "da_proteins")
    daOutputDirSource <- "results_fallback"
    if (!dirExistsFn(daOutputDir)) {
      dirCreateFn(daOutputDir, recursive = TRUE)
      daOutputDirCreated <- TRUE
    }
  }

  pathwayDir <- NULL
  pathwayDirSource <- NULL
  pathwayDirCreated <- FALSE

  if (!is.null(experimentPaths$pathway_dir) &&
      dirExistsFn(experimentPaths$pathway_dir)) {
    pathwayDir <- experimentPaths$pathway_dir
    pathwayDirSource <- "experiment_paths"
  } else {
    pathwayDir <- filePathFn(experimentPaths$results_dir, "pathway_enrichment")
    pathwayDirSource <- "results_fallback"
    if (!dirExistsFn(pathwayDir)) {
      dirCreateFn(pathwayDir, recursive = TRUE)
      pathwayDirCreated <- TRUE
    }
  }

  list(
    daOutputDir = daOutputDir,
    daOutputDirSource = daOutputDirSource,
    daOutputDirCreated = daOutputDirCreated,
    pathwayDir = pathwayDir,
    pathwayDirSource = pathwayDirSource,
    pathwayDirCreated = pathwayDirCreated
  )
}

resolveProtEnrichUniprotAnnotations <- function(workflowData,
                                                experimentPaths,
                                                currentS4Object = NULL,
                                                organismTaxid = NULL,
                                                globalEnv = .GlobalEnv,
                                                fileExistsFn = file.exists,
                                                readRdsFn = readRDS,
                                                filePathFn = file.path,
                                                dirExistsFn = dir.exists,
                                                dirCreateFn = dir.create,
                                                getUniprotAnnotationsFn = getUniprotAnnotations,
                                                catFn = cat) {
  uniprotDatCln <- NULL
  annotationSource <- NULL
  sourcePath <- NULL
  cacheDir <- NULL
  cacheDirCreated <- FALSE
  loadError <- NULL
  creationError <- NULL

  if (exists("uniprot_dat_cln", envir = globalEnv)) {
    uniprotDatCln <- get("uniprot_dat_cln", envir = globalEnv)
    annotationSource <- "global_environment"
    catFn("   ENRICHMENT Step: Found uniprot_dat_cln in global environment\n")
  } else if (!is.null(workflowData$uniprot_dat_cln)) {
    uniprotDatCln <- workflowData$uniprot_dat_cln
    assign("uniprot_dat_cln", uniprotDatCln, envir = globalEnv)
    annotationSource <- "workflow_data"
    catFn("   ENRICHMENT Step: Found uniprot_dat_cln in workflow_data\n")
  } else {
    catFn("   ENRICHMENT Step: No uniprot_dat_cln found - checking source directory\n")

    sourcePath <- filePathFn(experimentPaths$source_dir, "uniprot_dat_cln.RDS")
    if (fileExistsFn(sourcePath)) {
      catFn(sprintf("   ENRICHMENT Step: Found uniprot_dat_cln.RDS at %s\n", sourcePath))
      tryCatch({
        uniprotDatCln <- readRdsFn(sourcePath)
        workflowData$uniprot_dat_cln <- uniprotDatCln
        assign("uniprot_dat_cln", uniprotDatCln, envir = globalEnv)
        annotationSource <- "source_directory"
        catFn(sprintf(
          "   ENRICHMENT Step: Successfully loaded %d UniProt annotations from source directory\n",
          nrow(uniprotDatCln)
        ))
      }, error = function(e) {
        loadError <<- e$message
        uniprotDatCln <<- NULL
        catFn(sprintf(
          "   ENRICHMENT Step: Error loading UniProt from source directory: %s\n",
          e$message
        ))
      })
    } else {
      catFn(sprintf(
        "   ENRICHMENT Step: No uniprot_dat_cln.RDS found at %s\n",
        sourcePath
      ))
    }
  }

  if (is.null(uniprotDatCln)) {
    catFn("   ENRICHMENT Step: Attempting to create UniProt annotations on-the-fly\n")

    tryCatch({
      cacheDir <- filePathFn(experimentPaths$results_dir, "cache", "uniprot_annotations")
      if (!dirExistsFn(cacheDir)) {
        dirCreateFn(cacheDir, recursive = TRUE)
        cacheDirCreated <- TRUE
      }

      if (!is.null(currentS4Object) && !is.null(currentS4Object@protein_quant_table)) {
        uniprotDatCln <- getUniprotAnnotationsFn(
          input_tbl = currentS4Object@protein_quant_table,
          cache_dir = cacheDir,
          taxon_id = as.numeric(organismTaxid)
        )

        workflowData$uniprot_dat_cln <- uniprotDatCln
        assign("uniprot_dat_cln", uniprotDatCln, envir = globalEnv)
        annotationSource <- "generated"
        catFn("   ENRICHMENT Step: Successfully created uniprot_dat_cln on-the-fly\n")
      } else {
        catFn("   ENRICHMENT Step: No protein table available for annotation creation\n")
      }
    }, error = function(e) {
      creationError <<- e$message
      catFn(sprintf("   ENRICHMENT Step: Error creating UniProt annotations: %s\n", e$message))
    })
  }

  list(
    uniprotDatCln = uniprotDatCln,
    source = annotationSource,
    sourcePath = sourcePath,
    cacheDir = cacheDir,
    cacheDirCreated = cacheDirCreated,
    loadError = loadError,
    creationError = creationError
  )
}

resolveProtEnrichAnnotationMatching <- function(uniprotDatCln,
                                                daResultsForEnrichment,
                                                currentS4Object,
                                                matchAnnotationsFn = matchAnnotations,
                                                catFn = cat) {
  if (is.null(uniprotDatCln) ||
      is.null(daResultsForEnrichment) ||
      is.null(currentS4Object)) {
    return(list(
      attempted = FALSE,
      proteinIdColumn = NULL,
      annotationMatchResults = NULL,
      matchRate = NULL,
      warning = NULL
    ))
  }

  catFn("   ENRICHMENT Step: Attempting UniProt annotation matching\n")

  proteinIdCol <- tryCatch(
    currentS4Object@protein_id_column,
    error = function(e) "uniprot_acc"
  )
  annotationMatchResults <- NULL
  matchRate <- NULL
  warningMessage <- NULL

  tryCatch({
    annotationMatchResults <- matchAnnotationsFn(
      da_results_s4 = daResultsForEnrichment,
      uniprot_annotations = uniprotDatCln,
      protein_id_column = proteinIdCol,
      uniprot_id_column = "Entry",
      gene_names_column = "gene_names"
    )

    matchRate <- annotationMatchResults$match_statistics$match_rate
    catFn(sprintf(
      "   ENRICHMENT Step: Annotation matching completed - %d%% match rate\n",
      matchRate
    ))
  }, error = function(e) {
    warningMessage <<- e$message
    annotationMatchResults <<- NULL
    catFn(sprintf("   ENRICHMENT Step: Warning in annotation matching: %s\n", e$message))
    catFn("   ENRICHMENT Step: Continuing with enrichment analysis...\n")
  })

  list(
    attempted = TRUE,
    proteinIdColumn = proteinIdCol,
    annotationMatchResults = annotationMatchResults,
    matchRate = matchRate,
    warning = warningMessage
  )
}

