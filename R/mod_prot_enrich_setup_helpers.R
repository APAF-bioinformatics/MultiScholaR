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
