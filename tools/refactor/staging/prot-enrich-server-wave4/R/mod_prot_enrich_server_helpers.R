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

