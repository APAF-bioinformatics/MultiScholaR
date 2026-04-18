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

