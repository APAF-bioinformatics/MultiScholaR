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

