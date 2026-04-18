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

setupProtEnrichRawContrastNameReactive <- function(input,
                                                   createRawContrastNameReactiveFn = createProtEnrichRawContrastNameReactive) {
  rawContrastName <- createRawContrastNameReactiveFn(input = input)

  list(
    rawContrastName = rawContrastName,
    reason = "created"
  )
}

