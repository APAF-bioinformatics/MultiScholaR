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

