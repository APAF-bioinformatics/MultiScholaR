if (!methods::isClass("ProtEnrichDeterministicResult")) {
  methods::setClass(
    "ProtEnrichDeterministicResult",
    slots = c(result = "data.frame")
  )
}

buildProtEnrichDeterministicUniprotAnnotations <- function(currentS4Object = NULL,
                                                           organismTaxid = "9606") {
  proteinIds <- tryCatch({
    idCol <- currentS4Object@protein_id_column
    unique(as.character(currentS4Object@protein_quant_table[[idCol]]))
  }, error = function(e) character())

  proteinIds <- proteinIds[!is.na(proteinIds) & nzchar(proteinIds)]
  if (length(proteinIds) == 0L) {
    proteinIds <- c("P04637", "P10600", "Q07817", "P55957", "O14757")
  }

  data.frame(
    Entry = proteinIds,
    gene_names = paste0("GENE", seq_along(proteinIds)),
    Organism = if (identical(as.character(organismTaxid), "9606")) "Homo sapiens" else "Synthetic test organism",
    go_id_go_biological_process = "GO:0006915; GO:0007049",
    go_term_go_biological_process = "apoptotic process; cell cycle",
    go_id_go_molecular_function = "GO:0005515",
    go_term_go_molecular_function = "protein binding",
    go_id_go_cellular_compartment = "GO:0005739",
    go_term_go_cellular_compartment = "mitochondrion",
    stringsAsFactors = FALSE
  )
}

buildProtEnrichDeterministicGprofilerResult <- function(direction = "positive") {
  result <- data.frame(
    source = c("GO:BP", "GO:BP", "GO:CC", "KEGG", "REAC"),
    term_name = c(
      "E2E apoptotic process",
      "E2E cell cycle",
      "E2E mitochondrion",
      "E2E cell-cycle pathway",
      "E2E hemostasis"
    ),
    term_id = c("GO:0006915", "GO:0007049", "GO:0005739", "KEGG:04110", "REAC:R-HSA-109582"),
    term_size = c(1392L, 1780L, 892L, 131L, 418L),
    query_size = rep(10L, 5),
    intersection_size = c(4L, 5L, 3L, 3L, 2L),
    p_value = c(0.0012, 0.0034, 0.0085, 0.0127, 0.0243),
    precision = c(0.4, 0.5, 0.3, 0.3, 0.2),
    recall = c(0.00287, 0.00281, 0.00336, 0.0229, 0.00478),
    effective_domain_size = rep(20717L, 5),
    source_order = seq_len(5),
    parents = I(list(c("GO:0008219"), c("GO:0022402"), c("GO:0043231"), character(), character())),
    significant = TRUE,
    query = "query_1",
    intersection = c(
      "P04637/P10600/Q07817/P55957",
      "P04637/P10600/Q07817/P55957/O14757",
      "P04637/P10600/Q07817",
      "P04637/P10600/Q07817",
      "P04637/P10600"
    ),
    analysis_method = "gprofiler2",
    deterministic_backend = "gprofiler2_test_double",
    directionality = direction,
    stringsAsFactors = FALSE
  )

  list(
    result = result,
    meta = list(query_metadata = list(user_threshold = 0.05))
  )
}

buildProtEnrichDeterministicClusterProfilerResult <- function(direction = "up") {
  result <- data.frame(
    ID = c("GO:0006915", "GO:0007049", "GO:0005739"),
    Description = c("E2E apoptotic process", "E2E cell cycle", "E2E mitochondrion"),
    GeneRatio = c("4/10", "5/10", "3/10"),
    BgRatio = c("1392/20717", "1780/20717", "892/20717"),
    pvalue = c(0.0012, 0.0034, 0.0085),
    p.adjust = c(0.0060, 0.0085, 0.0142),
    qvalue = c(0.0058, 0.0083, 0.0139),
    geneID = c(
      "P04637/P10600/Q07817/P55957",
      "P04637/P10600/Q07817/P55957/O14757",
      "P04637/P10600/Q07817"
    ),
    Count = c(4L, 5L, 3L),
    analysis_method = "clusterprofiler",
    deterministic_backend = "clusterprofiler_test_double",
    directionality = direction,
    stringsAsFactors = FALSE
  )

  methods::new("ProtEnrichDeterministicResult", result = result)
}

writeProtEnrichDeterministicResultTable <- function(resultTable,
                                                    pathwayDir,
                                                    contrast,
                                                    direction,
                                                    writeTsvFn = readr::write_tsv,
                                                    dirCreateFn = dir.create) {
  dirCreateFn(pathwayDir, recursive = TRUE, showWarnings = FALSE)
  resultTable$parents <- if ("parents" %in% names(resultTable)) {
    vapply(resultTable$parents, paste, character(1), collapse = ", ")
  } else {
    NULL
  }
  writeTsvFn(
    resultTable,
    file.path(pathwayDir, paste0(contrast, "_", direction, "_enrichment_results.tsv"))
  )
}

buildProtEnrichDeterministicPlot <- function(backend, direction) {
  plotly::plot_ly(
    x = c("GO:0006915", "GO:0007049", "GO:0005739"),
    y = c(4, 5, 3),
    type = "bar",
    name = paste(backend, direction)
  )
}

runProtEnrichDeterministicProcessEnrichments <- function(da_results,
                                                        taxon_id,
                                                        up_cutoff = 0,
                                                        down_cutoff = 0,
                                                        q_cutoff = 0.05,
                                                        pathway_dir,
                                                        go_annotations = NULL,
                                                        exclude_iea = FALSE,
                                                        protein_id_column = "Protein.IDs",
                                                        contrast_names = NULL,
                                                        correction_method = "gSCS") {
  supportedOrganisms <- buildProtEnrichSupportedOrganisms()
  backend <- if (as.character(taxon_id) %in% supportedOrganisms$taxid) {
    "gprofiler2"
  } else {
    "clusterprofiler"
  }

  internalNames <- names(da_results@da_data)
  if (is.null(internalNames) || length(internalNames) == 0L) {
    internalNames <- "E2E_contrast"
  }
  outputNames <- contrast_names
  if (is.null(outputNames) || length(outputNames) != length(internalNames)) {
    outputNames <- internalNames
  }

  enrichmentResults <- createEnrichmentResults(da_results@contrasts)
  enrichmentData <- vector("list", length(outputNames))
  names(enrichmentData) <- outputNames
  enrichmentPlots <- vector("list", length(outputNames))
  names(enrichmentPlots) <- outputNames
  enrichmentPlotly <- vector("list", length(outputNames))
  names(enrichmentPlotly) <- outputNames
  enrichmentSummaries <- vector("list", length(outputNames))
  names(enrichmentSummaries) <- outputNames

  for (contrast in outputNames) {
    if (identical(backend, "gprofiler2")) {
      upResult <- buildProtEnrichDeterministicGprofilerResult("positive")
      downResult <- buildProtEnrichDeterministicGprofilerResult("negative")
      writeProtEnrichDeterministicResultTable(upResult$result, pathway_dir, contrast, "up")
      writeProtEnrichDeterministicResultTable(downResult$result, pathway_dir, contrast, "down")
    } else {
      upResult <- buildProtEnrichDeterministicClusterProfilerResult("up")
      downResult <- buildProtEnrichDeterministicClusterProfilerResult("down")
      writeProtEnrichDeterministicResultTable(upResult@result, pathway_dir, contrast, "up")
      writeProtEnrichDeterministicResultTable(downResult@result, pathway_dir, contrast, "down")
    }

    enrichmentData[[contrast]] <- list(up = upResult, down = downResult)
    enrichmentPlotly[[contrast]] <- list(
      up = buildProtEnrichDeterministicPlot(backend, "up"),
      down = buildProtEnrichDeterministicPlot(backend, "down")
    )
    enrichmentPlots[[contrast]] <- list(up = NULL, down = NULL)
    enrichmentSummaries[[contrast]] <- list(
      up = list(backend = backend, direction = "up", deterministic = TRUE),
      down = list(backend = backend, direction = "down", deterministic = TRUE)
    )
  }

  enrichmentResults@enrichment_data <- enrichmentData
  enrichmentResults@enrichment_plots <- enrichmentPlots
  enrichmentResults@enrichment_plotly <- enrichmentPlotly
  enrichmentResults@enrichment_summaries <- enrichmentSummaries

  enrichmentResults
}

resolveProtEnrichProcessEnrichmentsFn <- function(processEnrichmentsFn = processEnrichments,
                                                  deterministicProcessFn = runProtEnrichDeterministicProcessEnrichments,
                                                  isTestModeFn = is_test_mode) {
  if (isTRUE(isTestModeFn())) {
    return(deterministicProcessFn)
  }

  processEnrichmentsFn
}

resolveProtEnrichAnalysisMethod <- function(organismTaxid, supportedOrganisms) {
  organismTaxid <- normaliseProtEnrichTaxonId(organismTaxid)
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
