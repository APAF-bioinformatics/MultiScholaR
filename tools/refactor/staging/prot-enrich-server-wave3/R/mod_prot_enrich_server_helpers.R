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

