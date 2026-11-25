# KIT-107 W2 (D3, D4, D5): the enrichment core.
#   - Four entries: `runWorkflowEnrichment` (the run on an accepted DA manifest), `restoreWorkflowEnrichment` (the
#     completed state's results after an Open), `saveWorkflowEnrichmentState` (the completed state, with the step's
#     provenance, F-16) and `replayWorkflowEnrichmentStep` (the analysis record's enrichment step, F-16).
#   - The 19 helpers they share with the enrichment module, moved unchanged with their names from
#     mod_prot_enrich_{analysis,execution,input,observer,state}_helpers.R at the start commit; the module calls them by
#     name. Their `theObject`, state-manager and `globalEnv` fallbacks are kept unchanged (D7).
# The proteomics values the entries read from the `enrichment` group of the omic configuration (the double, the organism
# filter, the annotation resolver and matcher) are in utils_workflow_dia_enrichment.R. The core calls no Shiny UI
# function: it notifies and reports progress through the caller's callbacks (`resolveProtEnrichContrastsTable` keeps its
# `shiny::isolate` read of the bus).

resolveProtEnrichAnalysisInputColumns <- function(methodInfo,
                                                  daResultsForEnrichment,
                                                  currentS4Object = NULL,
                                                  id_choice = NULL) {
  # KIT-167 (8.3 item 3): the reference identifiers when chosen; "imported" or NULL keeps today's resolution.
  if (identical(id_choice, "reference")) {
    return(list(idColumn = "reference_accession_used", source = "reference_mapping", geneNameOverrideApplied = FALSE))
  }
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
                                                  excludeIea = FALSE,
                                                  goPropagation = NULL,
                                                  significanceRule = NULL) {
  validatedParams <- validateProtEnrichProcessParameters(
    organismTaxid = organismTaxid,
    upCutoff = upCutoff,
    downCutoff = downCutoff,
    qCutoff = qCutoff,
    correctionMethod = correctionMethod
  )
  taxonId <- as.numeric(validatedParams$organism_taxid)

  list(
    checkpointArgs = list(
      da_results_s4 = daResultsForEnrichment,
      taxon_id = taxonId,
      up_cutoff = validatedParams$up_cutoff,
      down_cutoff = validatedParams$down_cutoff,
      q_cutoff = validatedParams$q_cutoff,
      pathway_dir = pathwayDir,
      go_annotations = goAnnotations,
      exclude_iea = excludeIea,
      protein_id_column = proteinIdColumn,
      contrast_names = contrastNames,
      correction_method = validatedParams$correction_method
    ),
    processArgs = c(list(
      da_results = daResultsForEnrichment,
      taxon_id = taxonId,
      up_cutoff = validatedParams$up_cutoff,
      down_cutoff = validatedParams$down_cutoff,
      q_cutoff = validatedParams$q_cutoff,
      pathway_dir = pathwayDir,
      go_annotations = goAnnotations,
      exclude_iea = excludeIea,
      protein_id_column = proteinIdColumn,
      contrast_names = contrastNames,
      correction_method = validatedParams$correction_method
    ),
    # KIT-173: the GO parent-term rule reaches the kernel only when given (the clusterProfiler route); KIT-203: so
    # does the significance rule.
    if (!is.null(goPropagation)) list(go_propagation = goPropagation),
    if (!is.null(significanceRule)) list(significance_rule = significanceRule))
  )
}

buildProtEnrichAllContrastResults <- function(enrichmentResults,
                                              methodInfo,
                                              isEnrichResultFn = methods::is,
                                              catFn = cat) {

  allContrastResults <- list()

  for (rawContrastName in names(enrichmentResults@enrichment_data)) {

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
    } else if (identical(methodInfo$method, "clusterprofiler")) {
      clusterprofilerResults <- data.frame()

      if (!is.null(contrastEnrichment$up) || !is.null(contrastEnrichment$down)) {
        if (!is.null(contrastEnrichment$up)) {
          upResults <- contrastEnrichment$up
          upHasResultSlot <- isS4(upResults) && "result" %in% methods::slotNames(upResults)
          if (!is.null(upResults) &&
              (isTRUE(isEnrichResultFn(upResults, "enrichResult")) || upHasResultSlot) &&
              nrow(upResults@result) > 0) {
            upDf <- upResults@result
            upDf$directionality <- "up"
            upDf$analysis_method <- "clusterprofiler"
            clusterprofilerResults <- rbind(clusterprofilerResults, upDf)
          }
        }

        if (!is.null(contrastEnrichment$down)) {
          downResults <- contrastEnrichment$down
          downHasResultSlot <- isS4(downResults) && "result" %in% methods::slotNames(downResults)
          if (!is.null(downResults) &&
              (isTRUE(isEnrichResultFn(downResults, "enrichResult")) || downHasResultSlot) &&
              nrow(downResults@result) > 0) {
            downDf <- downResults@result
            downDf$directionality <- "down"
            downDf$analysis_method <- "clusterprofiler"
            clusterprofilerResults <- rbind(clusterprofilerResults, downDf)
          }
        }
      }

      contrastResults$clusterprofiler_results <- clusterprofilerResults
    }

    allContrastResults[[rawContrastName]] <- contrastResults
  }

  allContrastResults
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

#' Resolve the contrasts table the enrichment step works from.
#' @param workflowData Proteomics workflow bus; its `contrasts_tbl` is the owned payload (KIT-013).
#' @param globalEnv Optional lookup environment kept as a seam for scripts and tests; NULL (the
#'   application default) means no environment lookup.
#' @param existsFn,getFn Injectable lookup functions for the optional environment.
#' @return The contrasts table or NULL.
#' @noRd
resolveProtEnrichContrastsTable <- function(workflowData, globalEnv = NULL, existsFn = exists, getFn = get) {
  contrastsTbl <- if (is.null(workflowData)) NULL else shiny::isolate(workflowData$contrasts_tbl)
  if (is.null(contrastsTbl) && !is.null(globalEnv) && existsFn("contrasts_tbl", envir = globalEnv)) {
    contrastsTbl <- getFn("contrasts_tbl", envir = globalEnv)
  }
  contrastsTbl
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
    stringdbResults = if (!is.null(contrastResults)) contrastResults$stringdb_results else NULL,
    gprofilerResultsTable = if (!is.null(contrastResults)) contrastResults$gprofiler_results_table else NULL,
    clusterprofilerResultsTable = if (!is.null(contrastResults)) contrastResults$clusterprofiler_results_table else NULL
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
                                             globalEnv = NULL) {
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
    stateObject <- readDIAStateObject(workflowData$state_manager, currentState)
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

  if (is.null(designMatrix) && !is.null(globalEnv) && exists("design_matrix", envir = globalEnv)) {
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

normaliseProtEnrichTaxonId <- function(organismTaxid) {
  if (is.null(organismTaxid) || length(organismTaxid) != 1L || is.na(organismTaxid)) {
    stop("organism_taxid must be a single positive integer NCBI taxonomy ID.", call. = FALSE)
  }

  taxonId <- trimws(as.character(organismTaxid))
  if (!nzchar(taxonId) || !grepl("^[0-9]+$", taxonId)) {
    stop("organism_taxid must be a single positive integer NCBI taxonomy ID.", call. = FALSE)
  }

  taxonId <- sub("^0+", "", taxonId)
  if (!nzchar(taxonId) || identical(taxonId, "0")) {
    stop("organism_taxid must be a single positive integer NCBI taxonomy ID.", call. = FALSE)
  }

  taxonId
}

validateProtEnrichProcessParameters <- function(organismTaxid,
                                                upCutoff,
                                                downCutoff,
                                                qCutoff,
                                                correctionMethod = NULL) {
  taxonId <- normaliseProtEnrichTaxonId(organismTaxid)
  numericParams <- list(
    up_cutoff = upCutoff,
    down_cutoff = downCutoff,
    q_cutoff = qCutoff
  )
  numericParams <- lapply(names(numericParams), function(name) {
    value <- numericParams[[name]]
    if (is.null(value) || length(value) != 1L || is.na(value)) {
      stop(sprintf("%s must be a single finite numeric value.", name), call. = FALSE)
    }
    value <- suppressWarnings(as.numeric(value))
    if (!is.finite(value)) {
      stop(sprintf("%s must be a single finite numeric value.", name), call. = FALSE)
    }
    value
  }) |>
    stats::setNames(names(numericParams))

  if (numericParams$q_cutoff < 0 || numericParams$q_cutoff > 1) {
    stop("q_cutoff must be between 0 and 1.", call. = FALSE)
  }

  if (!is.null(correctionMethod) &&
      (length(correctionMethod) != 1L || is.na(correctionMethod) || !nzchar(trimws(as.character(correctionMethod))))) {
    stop("correction_method must be a single non-empty string.", call. = FALSE)
  }

  list(
    organism_taxid = taxonId,
    up_cutoff = numericParams$up_cutoff,
    down_cutoff = numericParams$down_cutoff,
    q_cutoff = numericParams$q_cutoff,
    correction_method = if (is.null(correctionMethod)) NULL else trimws(as.character(correctionMethod))
  )
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

    }
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
                                                   "ruv_covariates",
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


  currentDataState <- NULL
  availableStates <- NULL
  warningMessage <- NULL

  tryCatch({
    availableStates <- diaStateHistory(workflowData$state_manager)
    currentDataState <- detectFn(currentDataStates, ~ .x %in% availableStates)

  }, error = function(e) {
    warningMessage <<- e$message
  })

  list(
    attempted = TRUE,
    currentDataState = currentDataState,
    availableStates = availableStates,
    warning = warningMessage,
    updated = FALSE
  )
}

# ------------------------------------------------------------------------------------------------
# The entries (KIT-107 D3, D4, D5).
# ------------------------------------------------------------------------------------------------

#' The multiple-testing correction an enrichment run's backend applies (GUI audit R3, A107-03)
#'
#' The kernel's clusterProfiler branch calls `clusterProfiler::enricher(pAdjustMethod = "BH")` whatever correction was
#' selected (func_general_enrichment_orchestration.R); its g:Profiler branch passes the selection to `gost`. The
#' kernel takes the clusterProfiler branch exactly when the taxon is not a supported g:Profiler organism, which is
#' when `method$method` is `clusterprofiler` (the two organism lists are the same thirteen taxa).
#' @param method The method information (`resolveProtEnrichAnalysisMethod`).
#' @param requested The selected correction.
#' @return The correction applied: `"BH"` on the clusterProfiler backend, else `requested`.
#' @noRd
workflowEnrichmentEffectiveCorrection <- function(method, requested) {
    if (identical(method$method, "clusterprofiler")) "BH" else requested
}

#' The annotation problems an enrichment run reports (KIT-175 item 5)
#'
#' One notice line and one `annotation_error` outcome from the resolved table's `annotation_outcome` attribute (an
#' empty or incomplete retrieval, wherever the table came from) and the resolver's `loadError` and `creationError`.
#' When the table is empty or missing, the notice is the first sentence of the design step's notice for the problem
#' (the run is then refused).
#' @param annotations The resolved annotation table, or NULL.
#' @param resolved The resolver's list (`loadError`, `creationError`, `creationReason`, `creationDetail`).
#' @return `list(notice, outcome)`, or NULL when nothing occurred.
#' @noRd
workflowEnrichmentAnnotationProblem <- function(annotations, resolved) {
    outcome <- attr(annotations, "annotation_outcome")
    empty <- is.null(annotations) || nrow(annotations) == 0L
    load <- if (!is.null(resolved$loadError)) redactUniprotDetail(resolved$loadError)
    notices <- character()
    outcomes <- character()
    if (!is.null(load)) {
        notices <- c(notices, if (empty) sprintf("Protein annotation could not be read from the project file (%s).", load)
            else sprintf("The saved protein annotation could not be read (%s); it was retrieved again.", load))
        outcomes <- c(outcomes, sprintf("annotation file could not be read (%s)", load))
    }
    if (!is.null(resolved$creationError)) {
        notices <- c(notices, uniprotAnnotationSentence(list(reason = resolved$creationReason,
            detail = resolved$creationDetail)))
    }
    if (!is.null(outcome)) {
        if (identical(outcome$reason, "incomplete")) {
            notices <- c(notices, sprintf(paste("Enrichment ran with incomplete protein annotation: %d of %d proteins",
                "have no annotation and were left out."), outcome$accessions_not_looked_up, outcome$accessions_requested))
            outcomes <- c(outcomes, uniprotRecordSummary("incomplete", outcome))
        } else {
            notices <- c(notices, uniprotAnnotationSentence(outcome))
        }
    }
    if (!length(notices)) return(NULL)
    list(notice = paste(notices, collapse = " "), outcome = paste(outcomes, collapse = "; "))
}

#' The proteins without a gene name per comparison on the g:Profiler route (KIT-208, VAL-147 (3))
#'
#' Tested rows of each comparison's DA table (`enrichmentTestedRows`, the background in force, KIT-205) whose identifier
#' (the route's `gene_name`) is NA or empty; they leave the lists and the background g:Profiler is sent.
#' @param da_results The run's `da_results_for_enrichment` (after the organism filter, when it applied).
#' @param contrasts_tbl The contrasts table (display names).
#' @param id_column The identifier column sent to g:Profiler.
#' @return One string, "<comparison>: <unnamed> of <tested rows>" per comparison joined by "; ", or NULL when every
#'   tested protein of every comparison has a name.
#' @noRd
workflowEnrichmentUnnamedProteins <- function(da_results, contrasts_tbl, id_column) {
    tables <- lapply(Filter(Negate(is.null), da_results@da_data), function(table) enrichmentTestedRows(table)$rows)
    unnamed <- vapply(tables, function(table) {
        ids <- table[[id_column]]
        sum(is.na(ids) | !nzchar(as.character(ids)))
    }, 0L)
    if (!any(unnamed > 0L)) return(NULL)
    paste(sprintf("%s: %d of %d", enrichmentContrastDisplayNames(contrasts_tbl, names(tables)), unnamed,
        vapply(tables, nrow, 0L)), collapse = "; ")
}

#' Run the enrichment on an accepted DA manifest (KIT-107 D3)
#'
#' The module's manifest path at the start commit (`runProtEnrichAnalysisBody`, its setup, finalise and persist steps,
#' mod_prot_enrich_analysis_helpers.R:257-487, :528-603; mod_prot_enrich_state_helpers.R:371-515), in the same order,
#' without the tab's values: the caller writes them from the returned list. The manifest is read and never replaced,
#' so the caller keeps the manifest it holds (F107-01). The bus writes are today's, in today's order: `uniprot_dat_cln`
#' when the annotations are loaded, `enrichment_organism_filter`, `enrichment_analysis_results`,
#' `enrichment_ui_params`, the completed state with its step provenance (F-16, `saveWorkflowEnrichmentState`), then
#' `tab_status[[config$status_key]]`. Progress 0.2, 0.3, 0.5, 0.8 and 1.0 as today.
#' GUI audit R3, A107-08: when the state's commit fails, the failure is notified as a warning and `tab_status` is left
#' as it was; the returned `saved` says so (`saved$saved`, `saved$warning`).
#' @param workflow The bus: an environment, or reactive values read inside an observer or `shiny::isolate()`.
#' @param config The omic's `enrichment` group (`omicGuiConfig(<omic>)$enrichment`).
#' @param manifest The accepted DA result manifest.
#' @param parameters The run's inputs: `selected_contrast`, `organism_taxid`, `up_cutoff`, `down_cutoff`, `q_cutoff`,
#'   `correction_method`, `enable_organism_filter`, `go_propagation` (KIT-173; absent takes the shipped value,
#'   `goPropagationShipped`; used on the clusterProfiler route only), and `significance_rule` (KIT-203; absent takes
#'   `enrichmentSignificanceShipped`; clusterProfiler route only, where every run records the rule it applied as the
#'   outcome `significance_rule` and notifies it).
#' @param paths Experiment paths (`da_output_dir`, `pathway_dir`, `results_dir`, `source_dir`).
#' @param process_fn The kernel; NULL selects `config$kernel$select()` when called (the double in test mode).
#' @param notify `function(message, type, duration)`: the organism filter's warning, and on the clusterProfiler route
#'   the notice of the significance rule applied (KIT-203), with the sentence that the selected correction was not used
#'   when the selection is not `fdr` or `BH` (A107-03; then a warning); on the g:Profiler route a warning naming the
#'   lists not tested because none of their proteins has a gene name (KIT-208); NULL notifies nobody.
#' @param progress `function(value, detail)`; NULL reports nothing.
#' @param refuse_missing_da TRUE refuses a run whose DA table is missing for a contrast before anything is written
#'   (class `workflow_enrichment_input`; F107-03). FALSE is the start commit's behaviour: the run continues on the
#'   missing table, with the reader's warning.
#' @return `list(results, method, contrasts, selected, raw_contrast, contrasts_tbl, pathway_dir, annotation_match,
#'   organism_filter, service, reference, parameters, saved)`: the kernel's `EnrichmentResults` as the tab holds it
#'   (the completed state carries it with the source's `@args`); the method information; every contrast's display
#'   results (with their saved tables); the selected contrast's (`resolveProtEnrichSelectedContrastResults`); the raw
#'   contrast name; the contrasts table; the pathway directory; the annotation match (NULL when none); the
#'   organism-filter metadata written to the bus; the service identity; the product reference; the declared values
#'   the product's identity used (9, and the organism filter's 2 when it is on; P-28); and the state save's result.
#'   KIT-208: a g:Profiler comparison without a named protein refuses the run (class
#'   `workflow_enrichment_no_identifiers`, from the kernel) before anything is saved; the step's outcomes gain
#'   `enrichment_lists_not_tested` when a list was not tested and `enrichment_unnamed_proteins` when a comparison has
#'   proteins without a gene name. KIT-215 (VAL-170): on both routes; a GO comparison with lists and no annotated
#'   tested protein or no term of 10 to 500 of them, and a test that raised an error, refuse the run (classes
#'   `workflow_enrichment_no_identifiers`, `workflow_enrichment_test_error`); the list also carries `not_tested` and
#'   `no_list` (`diaEnrichmentNoIdentifiers`) when there are some, and the state's config commits them.
#' @noRd
runWorkflowEnrichment <- function(workflow, config, manifest, parameters, paths, process_fn = NULL, notify = NULL,
    progress = NULL, refuse_missing_da = FALSE) {
    notify <- notify %||% function(message, type, duration) invisible(NULL)
    progress <- progress %||% function(value, detail) invisible(NULL)
    if (!isTRUE(config$da_source$manifest(manifest))) {
        rlang::abort("Enrichment needs an accepted DA result manifest.", class = "workflow_enrichment_input")
    }
    manager <- workflow$state_manager
    # KIT-016: the source state's design and id table, hydrated without the quant table, for this run only.
    source_object <- readDIAExecutionComponents(manager, manifest$source, fields = config$da_source$components)
    selected <- parameters$selected_contrast
    progress(0.2, "Transforming DE data...")
    progress(0.3, "Creating DE results S4 object...")
    contrasts_tbl <- resolveProtEnrichContrastsTable(workflow)
    selection <- resolveProtEnrichSelectedDaResults(selected, manifest$results, contrasts_tbl)
    if (is.null(selection$selectedDaResults)) {
        rlang::abort(sprintf("Could not find DE results for contrast '%s'. Available contrasts: %s", selected,
            paste(selection$availableKeys, collapse = ", ")), class = "workflow_enrichment_input")
    }
    dependencies <- resolveProtEnrichRunDependencies(source_object, manifest$results, workflow, contrasts_tbl)
    contrasts_tbl <- dependencies$contrastsTbl
    if (is.null(contrasts_tbl) || is.null(dependencies$designMatrix)) {
        rlang::abort(paste("Missing contrasts_tbl or design_matrix required for enrichment analysis.",
            "Please ensure DE analysis was completed successfully."), class = "workflow_enrichment_input")
    }
    directories <- resolveProtEnrichOutputDirectories(paths)
    da_results <- config$da_source$read(contrasts_tbl = contrasts_tbl, design_matrix = dependencies$designMatrix,
        da_output_dir = directories$daOutputDir)
    missing <- names(Filter(is.null, da_results@da_data))
    if (isTRUE(refuse_missing_da) && length(missing)) {
        rlang::abort(sprintf("The DA results table is missing for contrast %s.", paste(missing, collapse = ", ")),
            class = "workflow_enrichment_input")
    }
    # KIT-119 (A107-06): each contrast's DA table as the run read it, before the organism filter.
    da_inputs <- lapply(da_results@da_data, function(table) diaAnnotationDigest(table) %||% NA_character_)
    # KIT-167 (8.3): with the reference identifiers chosen (gate on, a mapping on the head), each contrast's table gains
    # the group's `reference_accession_used` and unmapped groups leave the query lists and the background; the tables
    # on disk are unchanged.
    reference <- NULL
    if (identical(parameters$id_source, "reference")) {
        if (!referenceMappingEnabled()) {
            rlang::abort("Orthologue identifiers need the reference mapping tool, which is switched off.",
                class = "workflow_enrichment_input")
        }
        stored <- readReferenceMapping(manager)
        if (is.null(stored)) {
            rlang::abort("Orthologue identifiers need a reference mapping saved in this project; there is none.",
                class = "workflow_enrichment_input")
        }
        attached <- attachReferenceMapping(da_results, stored$mapping, source_object@protein_id_column,
            q_cutoff = parameters$q_cutoff, up_cutoff = parameters$up_cutoff, down_cutoff = parameters$down_cutoff)
        da_results <- attached$da_results
        reference <- list(stored = stored, counts = attached$counts,
            mapping_sha256 = stored$record$value[stored$record$field == "mapping_sha256"])
        # The tab's 2 x 2 counts and the unmapped changed groups, per contrast (5.6); written only on this route.
        workflow$enrichment_reference_mapping <- attached$counts
    }
    # KIT-175 (items 5, 5a): the annotation's problems are reported, and a run without annotation is refused before
    # any service is called and before the GO outcome step (D1, VAL-30; KIT-173 audit I04). The conditions of the
    # step's own retrieval are seen, not muffled.
    retrieval <- NULL
    keep <- function(w) retrieval <<- w
    # KIT-167 (8.3 item 3, exemption): with the reference identifiers the resolver returns the reference's table (the
    # shipped GO capture, or one fetch for a supplied reference) and KIT-175's refusal below tests that table.
    resolved <- withCallingHandlers(if (is.null(reference)) {
        config$annotation$resolve(workflowData = workflow, experimentPaths = paths, currentS4Object = source_object,
            organismTaxid = parameters$organism_taxid)
    } else {
        reference$annotation <- referenceAnnotationTable(reference$stored$mapping, reference$stored$record,
            project_dir = paths$source_dir)
        config$annotation$resolve(workflowData = workflow, experimentPaths = paths, currentS4Object = source_object,
            organismTaxid = parameters$organism_taxid, referenceAnnotation = reference$annotation$table)
    }, prot_annotation_uniprot_empty = keep, prot_annotation_uniprot_incomplete = keep)
    annotations <- resolved$uniprotDatCln
    if (!is.null(retrieval) && !is.null(annotations)) {
        attr(annotations, "annotation_outcome") <- uniprotOutcomeFromCondition(retrieval)
        # KIT-167: the project's annotation slot never holds a reference's table.
        if (is.null(reference)) workflow$uniprot_dat_cln <- annotations
    }
    annotation_problem <- workflowEnrichmentAnnotationProblem(annotations, resolved)
    if (if (is.null(annotations)) !is_test_mode() else nrow(annotations) == 0L) {
        if (!is.null(annotation_problem)) notify(annotation_problem$notice, "warning", NULL)
        rlang::abort(paste("Enrichment needs protein annotation and this project has none. Save the design again to",
            "retrieve it, then run the later steps again."), class = "workflow_enrichment_annotation")
    }
    if (!is.null(annotation_problem)) notify(annotation_problem$notice, "warning", NULL)
    match_object <- source_object
    if (!is.null(reference)) match_object@protein_id_column <- "reference_accession_used"
    matching <- config$annotation$match(uniprotDatCln = annotations, daResultsForEnrichment = da_results,
        currentS4Object = match_object)
    filter_applied <- FALSE
    filter_stats <- list(proteins_before = 0, proteins_after = 0, proteins_removed = 0)
    # KIT-167 (8.3 item 2, D05): with the reference identifiers the organism filter is not applied (the join keeps only
    # the strain's entries) and none of its names is declared.
    filter_enabled <- isTRUE(parameters$enable_organism_filter) && !is.null(config$organism_filter) &&
        is.null(reference)
    if (filter_enabled) {
        target <- as.character(parameters$organism_taxid)
        mapping <- config$organism_filter$mapping(workflowData = workflow, uniprotDatCln = annotations,
            targetTaxon = target)$organismMapping
        if (!is.null(mapping) && nrow(mapping) > 0) {
            filtered <- config$organism_filter$apply(daResultsForEnrichment = da_results, organismMapping = mapping,
                targetTaxon = target, currentS4Object = source_object)
            da_results <- filtered$daResultsForEnrichment
            filter_applied <- filtered$filterApplied
            filter_stats <- filtered$filterStats
        } else {
            notify("Multi-species filtering enabled but no organism mapping available. Proceeding without filtering.",
                "warning", 8)
        }
    }
    organism_filter <- persistProtEnrichOrganismFilterMetadata(workflowData = workflow,
        organismFilterEnabled = isTRUE(parameters$enable_organism_filter) && is.null(reference),
        organismFilterApplied = filter_applied,
        targetTaxonId = parameters$organism_taxid, filterStats = filter_stats)
    progress(0.5, "Running enrichment analysis...")
    # `[[` (exact): `$organism` would partially match `organism_filter`.
    method <- config[["organism"]]$method(organismTaxid = parameters$organism_taxid,
        supportedOrganisms = config[["organism"]]$supported())$methodInfo
    columns <- resolveProtEnrichAnalysisInputColumns(methodInfo = method, daResultsForEnrichment = da_results,
        currentS4Object = source_object, id_choice = if (!is.null(reference)) "reference")
    # KIT-173: the GO parent-term rule, on the clusterProfiler route only (the g:Profiler route's arguments, declared
    # names and record are unchanged).
    go_propagation <- if (identical(method$method, "clusterprofiler"))
        checkGoPropagation(parameters$go_propagation %||% goPropagationShipped)
    # KIT-203: the significance rule, on the clusterProfiler route only.
    significance_rule <- if (identical(method$method, "clusterprofiler"))
        checkEnrichmentSignificanceRule(parameters$significance_rule %||% enrichmentSignificanceShipped)
    arguments <- buildProtEnrichProcessEnrichmentsArgs(daResultsForEnrichment = da_results,
        organismTaxid = parameters$organism_taxid, upCutoff = parameters$up_cutoff, downCutoff = parameters$down_cutoff,
        qCutoff = parameters$q_cutoff, pathwayDir = directories$pathwayDir, goAnnotations = annotations,
        proteinIdColumn = columns$idColumn, contrastNames = names(manifest$results),
        correctionMethod = parameters$correction_method, excludeIea = config$exclude_iea, goPropagation = go_propagation,
        significanceRule = significance_rule)
    kernel <- process_fn %||% config$kernel$select()
    .capture_checkpoint(arguments$checkpointArgs, "cp10", "enrichment_input")
    # GUI audit R3, P-28: the organism filter is declared when it is on, and what it did is recorded as the step's
    # outcomes (provenance only). A run with the filter off declares and records what it did before.
    run_arguments <- arguments$processArgs
    outcomes <- NULL
    # A107-03: the correction the backend applies is the declared one; the kernel is still called with the selection.
    correction <- workflowEnrichmentEffectiveCorrection(method, parameters$correction_method)
    if (!identical(correction, parameters$correction_method)) {
        run_arguments$correction_method_effective <- correction
        outcomes <- list(correction_method_requested = parameters$correction_method)
    }
    # KIT-203 (VAL-86, OD10): every clusterProfiler run is told which rule calls a term significant; the sentence that
    # the selected correction was not used follows when the selection is not g:Profiler's `fdr` (Benjamini-Hochberg
    # under another name) or BH, and makes the notice a warning.
    if (!is.null(significance_rule)) {
        unused <- !isTRUE(parameters$correction_method %in% c("fdr", "BH"))
        notify(paste0(ENRICHMENT_SIGNIFICANCE_RULES$notice[ENRICHMENT_SIGNIFICANCE_RULES$rule == significance_rule],
            if (unused) sprintf(" The selected correction, %s, applies to g:Profiler only and was not used.",
                parameters$correction_method)), if (unused) "warning" else "message", 8)
    }
    if (!is.null(reference)) {
        run_arguments$id_source <- "reference"
        run_arguments$reference_mapping_sha256 <- reference$mapping_sha256
        run_arguments$go_table_sha256 <- reference$annotation$go_table_sha256
        outcomes <- c(outcomes, list(reference_mapping_groups = reference$counts))
    }
    if (filter_enabled) {
        run_arguments$organism_filter_enabled <- TRUE
        run_arguments$organism_filter_taxon <- as.character(parameters$organism_taxid)
        outcomes <- c(outcomes, list(organism_filter_applied = isTRUE(filter_applied),
            organism_filter_rows_before = filter_stats$proteins_before,
            organism_filter_rows_after = filter_stats$proteins_after))
    }
    # KIT-119 (A107-06): the identities of the inputs the run read, recorded with the step's outcomes (provenance
    # only): the sha256 of the annotation table and of each contrast's DA table (`diaAnnotationDigest`: a table's
    # displayable columns as text, so the table in memory and a read-back of it agree). NA for no table or an empty one.
    outcomes <- c(outcomes, list(annotation_sha256 = diaAnnotationDigest(annotations) %||% NA_character_,
        da_inputs_sha256 = da_inputs))
    # KIT-175 (item 5): present only when the annotation had a problem on a run that proceeds.
    if (!is.null(annotation_problem) && nzchar(annotation_problem$outcome)) {
        outcomes$annotation_error <- annotation_problem$outcome
    }
    if (!is.null(go_propagation) && !identical(go_propagation, "none")) {
        outcomes <- c(outcomes, workflowEnrichmentGoOutcomes(annotations, go_propagation, da_results, columns$idColumn))
    }
    # KIT-203 (OD10): the rule applied, provenance only (not identity), last so positional reads keep their entries.
    if (!is.null(significance_rule)) outcomes <- c(outcomes, list(significance_rule = significance_rule))
    declared <- config$run(workflow, manifest, run_arguments, process_fn = kernel)
    # KIT-208 (VAL-147 (1), SD-1 (a)): a g:Profiler list whose proteins all lack a gene name was not sent; the run says
    # so once and records it, so that its empty result is not read as "no term enriched". Absent otherwise.
    # KIT-215 (VAL-170): a GO list none of whose proteins has an annotation the analysis tests likewise, in its own
    # sentence; each route has one reason, so the route names it.
    not_tested <- declared$not_tested
    if (NROW(not_tested) > 0L) {
        go <- identical(method$method, "clusterprofiler")
        notify(paste(if (go) sprintf(paste("Comparison %s, %s: none of its %d proteins has a Gene Ontology annotation",
            "that this analysis tests, so this list was not tested (its empty result does not mean that no term is",
            "enriched)."), not_tested$comparison, not_tested$direction, not_tested$proteins) else
            sprintf(paste("Comparison %s, %s: none of its %d proteins has a gene name, so this list was not",
            "sent to g:Profiler and was not tested (its empty result does not mean that no term is enriched)."),
            not_tested$comparison, not_tested$direction, not_tested$proteins), collapse = " "), "warning", NULL)
        outcomes$enrichment_lists_not_tested <- paste(not_tested$comparison, not_tested$direction, sep = ": ",
            collapse = "; ")
    }
    # KIT-208 (VAL-147 (3), SD-3 (c)): on the g:Profiler route, the proteins of each comparison without a gene name,
    # which leave its lists and its background; recorded only when there are some.
    if (identical(method$method, "gprofiler2")) {
        outcomes$enrichment_unnamed_proteins <- workflowEnrichmentUnnamedProteins(da_results, contrasts_tbl,
            columns$idColumn)
    }
    results <- declared$results
    contrast_results <- buildProtEnrichAllContrastResults(enrichmentResults = results, methodInfo = method)
    # KIT-018: each contrast's result tables are also saved display tables (content-addressed).
    for (contrast in names(contrast_results)) {
        contrast_results[[contrast]]$gprofiler_results_table <- saveDIADisplayTable(manager,
            contrast_results[[contrast]]$gprofiler_results)
        contrast_results[[contrast]]$clusterprofiler_results_table <- saveDIADisplayTable(manager,
            contrast_results[[contrast]]$clusterprofiler_results)
    }
    shown <- resolveProtEnrichSelectedContrastResults(selected, contrast_results, contrasts_tbl)
    progress(0.8, "Storing results...")
    workflow$enrichment_analysis_results <- buildProtEnrichAnalysisResultsPayload(
        gprofilerResults = if (shown$found) shown$gprofilerResults,
        clusterprofilerResults = if (shown$found) shown$clusterprofilerResults,
        stringdbResults = if (shown$found) shown$stringdbResults, fullEnrichmentResults = results,
        selectedContrast = selected, analysisMethod = method$method, upCutoff = parameters$up_cutoff,
        downCutoff = parameters$down_cutoff, qCutoff = parameters$q_cutoff, organismTaxid = parameters$organism_taxid)
    stamped <- propagateProtEnrichResultsArgs(enrichmentResults = results, currentS4Object = source_object,
        selectedContrast = selected, methodInfo = method, upCutoff = parameters$up_cutoff,
        downCutoff = parameters$down_cutoff, qCutoff = parameters$q_cutoff, organismTaxid = parameters$organism_taxid,
        pathwayDir = directories$pathwayDir)$enrichmentResults
    ui_parameters <- propagateProtEnrichUiParams(currentS4Object = source_object, workflowData = workflow,
        selectedContrast = selected, methodInfo = method, upCutoff = parameters$up_cutoff,
        downCutoff = parameters$down_cutoff, qCutoff = parameters$q_cutoff, organismTaxid = parameters$organism_taxid)
    updateProtEnrichStateManagerUiParams(workflowData = workflow, storedUiParams = ui_parameters$storedUiParams)
    saved <- saveWorkflowEnrichmentState(workflow, config, stamped, declared$parameters, method, parameters,
        directories$pathwayDir, outcomes = outcomes, not_tested = declared$not_tested, no_list = declared$no_list,
        source = manifest$source)
    if (isTRUE(saved$saved)) {
        status <- workflow$tab_status
        status[[config$status_key]] <- "complete"
        workflow$tab_status <- status
    } else {
        # A107-08: the commit's failure was swallowed and the run reported complete.
        notify(sprintf("The enrichment results were not saved to the workflow state: %s",
            saved$warning %||% "unknown error"), "warning", 10)
    }
    progress(1.0, "Complete!")
    out <- list(results = results, method = method, contrasts = contrast_results, selected = shown,
        raw_contrast = selection$rawContrastName, contrasts_tbl = contrasts_tbl, pathway_dir = directories$pathwayDir,
        annotation_match = matching$annotationMatchResults, organism_filter = organism_filter, service = declared$service,
        reference = declared$reference, parameters = declared$parameters, saved = saved)
    # KIT-215: the lists not tested and the directions with no list, for the tab; absent when there are none.
    out$not_tested <- declared$not_tested
    out$no_list <- declared$no_list
    out
}

#' The GO parent-term outcomes of a propagated clusterProfiler run (KIT-173)
#'
#' Computed from the annotation table and the gene sets the kernel builds from it (`buildUniprotGoGeneSets()`, the
#' same call `processEnrichments` makes), not from the result object, whose class does not change. Provenance only.
#' @param annotations The annotation table the run passes to the kernel (its `go_ids_not_in_go_db` attribute, set by
#'   `getUniprotAnnotations`, when present).
#' @param go_propagation The run's propagating value.
#' @param da_results The DA results the kernel receives.
#' @param id_column The protein id column of the DA tables.
#' @return Named list: `go_db_version`, `go_source_date`, `go_relations` (plain words), `go_ids_not_in_go_db` (by
#'   class), `go_ids_not_in_go_db_total`, `go_annotation_go_db_version`, `go_annotation_go_source_date` (the release
#'   in force when the table was built), `go_pairs_direct`, `go_pairs_propagated`, `go_root_only_proteins` and
#'   `go_unreviewed_share`; `"not recorded"` where the table does not say.
#' @noRd
workflowEnrichmentGoOutcomes <- function(annotations, go_propagation, da_results, id_column) {
    sets <- buildUniprotGoGeneSets(annotations, go_propagation)
    lost <- attr(annotations, "go_ids_not_in_go_db", exact = TRUE)
    da_ids <- unique(unlist(lapply(da_results@da_data, function(table)
        if (is.data.frame(table) && id_column %in% names(table)) as.character(table[[id_column]])), use.names = FALSE))
    pairs <- attr(sets, "pairs_propagated")
    unreviewed <- attr(sets, "unreviewed_pairs", exact = TRUE)
    plain <- c(isa = "is a", `part of` = "part of", regulates = "regulates",
        `positively regulates` = "positively regulates", `negatively regulates` = "negatively regulates")
    list(go_db_version = attr(sets, "go_db_version"), go_source_date = attr(sets, "go_source_date"),
        go_relations = as.list(unname(plain[attr(sets, "relations")])),
        go_ids_not_in_go_db = if (is.null(lost)) "not recorded" else lost[c("obsolete", "secondary", "unknown")],
        go_ids_not_in_go_db_total = if (is.null(lost)) "not recorded" else lost$total,
        go_annotation_go_db_version = lost$go_db_version %||% "not recorded",
        go_annotation_go_source_date = lost$go_source_date %||% "not recorded",
        go_pairs_direct = attr(sets, "pairs_direct"), go_pairs_propagated = pairs,
        go_root_only_proteins = sum(da_ids %in% attr(sets, "root_only_proteins")),
        go_unreviewed_share = if (is.null(unreviewed)) "not recorded" else if (pairs > 0) unreviewed / pairs else 0)
}

#' Save the completed enrichment state with its step provenance (KIT-107 D3, D4; F-16)
#'
#' Today's `saveProtEnrichCompletedState` (mod_prot_enrich_state_helpers.R:210-262 at the start commit): the same state
#' name, config object, description and failure warning. Immediately before the commit it sets the step's provenance as
#' the protein steps set theirs (utils_workflow_prot_declarations.R:1158-1164): the registry entry's operation, its
#' default variant, its declared parameter names and fixed values, and the declared values the product's identity used.
#' GUI audit R3, P-28: the names are those of the declared values the run used (the organism filter's two only when it
#' is on), and `outcomes` is committed in the state's config, where the analysis record reads a step's outcomes.
#' A failed commit clears the pending provenance, so the next commit is not attributed to the enrichment.
#' @param workflow The bus (its `state_manager`).
#' @param config The omic's `enrichment` group (`operation`, `state_name`).
#' @param results The `EnrichmentResults` to commit.
#' @param declared The declared values (`runDIAReferenceEnrichment`'s `parameters`).
#' @param method The method information (`method`, `supported`).
#' @param parameters The run's inputs (`selected_contrast`, cutoffs, `organism_taxid`).
#' @param pathway_dir The pathway directory.
#' @param outcomes Named list of provenance-only values of the run (the organism filter's applied flag and row
#'   counts; the selected correction when the backend applied another, A107-03; the sha256 of the annotation table and
#'   of each contrast's DA table the run read, KIT-119 A107-06), or NULL for none.
#' @param not_tested,no_list KIT-215: the lists not tested and the directions with no list
#'   (`diaEnrichmentNoIdentifiers`), committed in the config beside `outcomes` only when they have rows (display
#'   only: neither identity nor record).
#' @param source KIT-216 (VAL-175): the DA manifest's source, the table the enrichment was run on. On a reference-state
#'   owner the result is committed as a product of that table (`commitWorkflowEnrichmentProduct`): the head and the
#'   history stay the table. Another owner commits the state as before.
#' @return `list(attempted, saved, warning)`.
#' @noRd
saveWorkflowEnrichmentState <- function(workflow, config, results, declared, method, parameters, pathway_dir,
    outcomes = NULL, not_tested = NULL, no_list = NULL, source = NULL) {
    manager <- workflow$state_manager
    warning_message <- NULL
    saved <- FALSE
    tryCatch({
        results <- portableWorkflowEnrichmentResults(results) # F107-05
        entry <- workflowOperationEntry(config$operation()$name)
        # KIT-167: the reference identifiers' three values are recorded after the declaration's names when the run
        # used them (a run without them records exactly what it did before), so that its replay checks them.
        setWorkflowStepProvenance(manager, entry$name, entry$default_variant,
            c(intersect(entry$parameters, names(declared)), intersect(referenceEnrichmentNames, names(declared))),
            entry$declaration$fixed, values = declared)
        config_object <- c(list(
                selected_contrast = parameters$selected_contrast,
                analysis_method = method$method,
                organism_supported = method$supported,
                up_cutoff = parameters$up_cutoff,
                down_cutoff = parameters$down_cutoff,
                q_cutoff = parameters$q_cutoff,
                organism_taxid = parameters$organism_taxid,
                pathway_dir = pathway_dir),
                if (length(outcomes)) list(outcomes = outcomes),
                if (NROW(not_tested) > 0L) list(not_tested = not_tested),
                if (NROW(no_list) > 0L) list(no_list = no_list))
        description <- paste("Enrichment analysis completed using", method$method, "for contrast:",
            parameters$selected_contrast)
        if (usesDIAReferenceStates(manager)) {
            commitWorkflowEnrichmentProduct(manager, config$state_name, checkpointDIAExecutionResult(manager, results),
                config_object, description, source)
        } else {
            saveDIAOwnerState(manager, state_name = config$state_name, s4_data_object = results,
                config_object = config_object, description = description)
        }
        saved <- TRUE
    }, error = function(e) {
        clearWorkflowStepProvenance(manager)
        warning_message <<- e$message
    })
    list(attempted = TRUE, saved = saved, warning = warning_message)
}

#' The completed enrichment's results made portable for the commit (KIT-107 F107-05)
#'
#' The live kernel's static plots (`processEnrichments`: ggplot2, then `plotly::ggplotly` and `ggplot2::ggsave`) keep
#' the environments they were made in: a ggplot2 layer records its caller's environment, and building the plot leaves
#' ggproto method closures whose `self` (and a scale transform's `reverse`) are unforced promises. Through the kernel's
#' frame they can also reach the DA tables as readr read them, with a live parser pointer. The commit's portability check
#' (`checkSavedWorkflowPayload`) refuses both, so a run on the live kernel never saved its state. Here every lazy binding
#' reachable from the results is forced (it then holds the value it would have read) and every readr table's parser
#' diagnostics are materialised (`prepareDIAParserTable`), in place, before the commit. No value, method or result
#' changes. A binding that cannot be forced is left as it is, and the check refuses it as before.
#' @param results The `EnrichmentResults` to commit.
#' @return `results`, with the reachable readr tables' diagnostics materialised.
#' @noRd
portableWorkflowEnrichmentResults <- function(results) {
    seen <- new.env(parent = emptyenv())
    visit <- function(value) {
        if (inherits(value, "spec_tbl_df")) return(prepareDIAParserTable(value))
        if (is.environment(value)) {
            if (isNamespace(value) || identical(value, globalenv()) || identical(value, baseenv()) ||
                identical(value, emptyenv())) return(value)
            id <- rlang::obj_address(value)
            if (exists(id, envir = seen, inherits = FALSE)) return(value)
            assign(id, TRUE, envir = seen)
            for (name in ls(value, all.names = TRUE)) {
                if (bindingIsActive(name, value)) next
                if (rlang::env_binding_are_lazy(value, name)) {
                    tryCatch(get(name, envir = value, inherits = FALSE), error = function(e) NULL)
                    if (rlang::env_binding_are_lazy(value, name)) next
                }
                # An empty argument binding (a plotting closure's) holds no value to visit.
                if (rlang::is_missing(.subset2(value, name))) next
                current <- .subset2(value, name)
                updated <- visit(current)
                if (!identical(updated, current) && !bindingIsLocked(name, value)) assign(name, updated, envir = value)
            }
            visit(parent.env(value))
            return(value)
        }
        if (is.function(value)) {
            if (!is.primitive(value)) visit(environment(value))
            return(value)
        }
        if (isS4(value) && length(methods::slotNames(value))) {
            for (slot in methods::slotNames(value)) {
                current <- methods::slot(value, slot)
                updated <- visit(current)
                if (!identical(updated, current)) methods::slot(value, slot, check = FALSE) <- updated
            }
        } else if (is.list(value) && !is.pairlist(value)) {
            classes <- oldClass(value)
            for (i in seq_along(value)) {
                current <- .subset2(value, i)
                updated <- visit(current)
                if (!identical(updated, current)) {
                    oldClass(value) <- NULL
                    value[i] <- list(updated)
                    oldClass(value) <- classes
                }
            }
        }
        others <- setdiff(names(attributes(value)), c("names", "class", "row.names",
            if (isS4(value)) methods::slotNames(value)))
        for (name in others) {
            current <- attr(value, name, exact = TRUE)
            updated <- visit(current)
            if (!identical(updated, current)) attr(value, name) <- updated
        }
        value
    }
    visit(results)
}

#' The completed enrichment's results after an Open (KIT-107 D3; KIT-024)
#'
#' Today's restore body (mod_prot_enrich_state_helpers.R:588-607 at the start commit) without its writes to the tab:
#' the completed state's `EnrichmentResults` rebuilt into the per-contrast results exactly as the run builds them, each
#' contrast's tables saved as display tables.
#' @param manager The session owner.
#' @param config The omic's `enrichment` group (`state_name`).
#' @return NULL when the state is not in the history or is not an `EnrichmentResults`; otherwise
#'   `list(all_enrichment_results, enrichment_results_full, method)`, with KIT-215's `not_tested` and `no_list` from the
#'   state's config when it carries them.
#' @noRd
restoreWorkflowEnrichment <- function(manager, config) {
    if (is.null(manager)) return(NULL)
    # KIT-216 (VAL-175): on a reference-state owner the enrichment is the current product of the head's table (a store
    # adopted from an archive written before it is converted on adoption); another owner reads it from its history.
    # r1.4 (audit F01, F02): an old-shape store whose head is the enrichment state itself (not converted; in the application
    # only test-built stores) is read from the history as before; an enrichment below a table head never is.
    if (usesDIAReferenceStates(manager) && !identical(manager$current_state, config$state_name)) {
        if (is.null(diaEnrichmentProduct(manager))) return(NULL)
    } else if (!config$state_name %in% unlist(diaStateHistory(manager))) {
        return(NULL)
    }
    results <- readDIAStateObject(manager, config$state_name)
    if (!methods::is(results, "EnrichmentResults")) return(NULL)
    method <- list(method = diaEnrichmentMethodOf(results))
    contrast_results <- buildProtEnrichAllContrastResults(enrichmentResults = results, methodInfo = method)
    for (contrast in names(contrast_results)) {
        contrast_results[[contrast]]$gprofiler_results_table <- saveDIADisplayTable(manager,
            contrast_results[[contrast]]$gprofiler_results)
        contrast_results[[contrast]]$clusterprofiler_results_table <- saveDIADisplayTable(manager,
            contrast_results[[contrast]]$clusterprofiler_results)
    }
    out <- list(all_enrichment_results = contrast_results, enrichment_results_full = results, method = method)
    # KIT-215: the record's lists not tested and directions with no list (absent from a record without them).
    # A state stored as a whole legacy object has no component metadata and carries neither (as the record writer
    # reads it, `buildWorkflowAnalysisRecord`).
    reference <- getWorkflowStateReference(manager, config$state_name)
    if (!is.null(reference)) {
        recorded <- tryCatch(readComponentMetadata(manager, reference)$config,
            saved_entity_unsupported = function(e) NULL)
        out$not_tested <- recorded$not_tested
        out$no_list <- recorded$no_list
    }
    out
}

#' Replay a recorded enrichment step after the replayed DA (KIT-107 D5; F-16)
#'
#' 1. The replayed manifest's DA tables are published into `<scratch>/enrichment/{da,graphs}` through the DA
#'    publication body the GUI uses (`produceDIADaPublication`, assets not saved), with the annotation table the GUI's
#'    DA would have had: the prepared project's `<source_dir>/uniprot_dat_cln.RDS`, else none.
#' 2. `runWorkflowEnrichment` runs on a bus holding the owner and the record's DA contrasts, with the record's taxon,
#'    cutoffs and correction, the organism filter on when the record declares it (P-28), the first contrast's friendly
#'    name selected (a display choice) and the kernel the GUI would select. Its annotation resolver takes the GUI's
#'    order with an empty bus: the prepared project's `uniprot_dat_cln.RDS`, then the test-mode table, then the
#'    network.
#' 3. The run commits the completed state with the same provenance (D4). The declared values it produced must equal
#'    the record's 9 (`protein_id_column` and `contrast_names` among them); a difference refuses the replay, naming
#'    the parameter. `results` is compared by its contrasts: its values are the DA tables' entity ids, which are
#'    store-local (a saved entity's id hashes its payload's file name, `saveWorkflowTable`), so the replay store's ids
#'    never equal the origin's; the tables' content is the replayed DA's (proposed F107-04).
#' 4. A recorded organism filter (P-28) is applied with the mapping the replay can resolve: the prepared project's
#'    annotation table (the import's organism mapping is not in the record). What it did must equal the record's
#'    outcomes (applied, rows before and after); a difference refuses the replay.
#'
#' Limits: a live g:Profiler result reproduces only when the service's database version is the recorded one (the
#' record's `external_resources`; the replay warns when it is not); the annotation table is the prepared project's and is not hash-checked.
#' @param manager The replay's owner.
#' @param step The record's enrichment step.
#' @param parameters The step's declared values, decoded (`workflowReplayParameters`).
#' @param manifest The replayed DA's manifest.
#' @param contrasts The record's DA contrasts table.
#' @param paths Experiment paths of the prepared project.
#' @param scratch Directory for the publication and the pathway files.
#' @param recorded_resources The record's `external_resources`; a service the replay queried that is not among them
#'   (by name and release) gives a warning of class `workflow_replay_resource`.
#' @return The run's result (`runWorkflowEnrichment`).
#' @noRd
replayWorkflowEnrichmentStep <- function(manager, step, parameters, manifest, contrasts, paths, scratch,
    recorded_resources = NULL) {
    config <- omicGuiConfig("proteomics")$enrichment
    staged <- list(da_output_dir = file.path(scratch, "enrichment", "da"),
        publication_graphs_dir = file.path(scratch, "enrichment", "graphs"))
    pathway_dir <- file.path(scratch, "enrichment", "pathway")
    for (directory in c(staged, pathway_dir)) dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    annotation_file <- if (is.null(paths$source_dir)) "" else file.path(paths$source_dir, "uniprot_dat_cln.RDS")
    annotations <- if (nzchar(annotation_file) && file.exists(annotation_file)) readRDS(annotation_file) else data.frame()
    produceDIADaPublication(manager, manifest, annotations, manifest$parameters$qvalue, staged, save_assets = FALSE)
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    workflow$contrasts_tbl <- contrasts
    workflow$tab_status <- list()
    run_paths <- utils::modifyList(as.list(paths), list(da_output_dir = staged$da_output_dir, pathway_dir = pathway_dir))
    run_parameters <- list(selected_contrast = contrasts$friendly_names[[1L]],
        organism_taxid = format(parameters$taxon_id, scientific = FALSE, trim = TRUE), up_cutoff = parameters$up_cutoff,
        down_cutoff = parameters$down_cutoff, q_cutoff = parameters$q_cutoff,
        correction_method = parameters$correction_method,
        enable_organism_filter = isTRUE(parameters$organism_filter_enabled),
        # KIT-173: a record made before parent terms, or without them, replays without them.
        go_propagation = parameters$go_propagation %||% "none",
        # KIT-203 (OD6): a record made before the rule, or with the q-value, replays with the q-value.
        significance_rule = parameters$significance_rule %||% "qvalue")
    # KIT-167: a step recorded with the reference identifiers replays with them, the tool switched on for this step
    # only; its recorded mapping and GO-table SHA-256 are then compared below like every declared value.
    reference_ids <- identical(parameters$id_source, "reference")
    if (reference_ids) run_parameters$id_source <- "reference"
    if (reference_ids) {
        previous <- options(MultiScholaR.reference_mapping = TRUE)
        on.exit(options(previous), add = TRUE)
    }
    run <- runWorkflowEnrichment(workflow, config, manifest, run_parameters, run_paths)
    if (reference_ids) options(previous)
    # The record holds each value through its JSON encoding; the produced values are compared in the same form.
    as_recorded <- function(value) decodeWorkflowRecordValue(jsonlite::fromJSON(jsonlite::toJSON(
        encodeWorkflowRecordValue(value), auto_unbox = TRUE, digits = I(17), null = "null", na = "null"),
        simplifyVector = FALSE))
    produced <- run$parameters
    produced$results <- names(produced$results)
    parameters$results <- names(parameters$results)
    for (name in names(parameters)) {
        if (!identical(as_recorded(produced[[name]]), parameters[[name]])) {
            cli::cli_abort("replayed step {.val {step$state}} produced a different {.val {name}} from the record",
                class = "workflow_replay_invalid")
        }
    }
    if (isTRUE(parameters$organism_filter_enabled)) {
        produced_filter <- list(organism_filter_applied = isTRUE(run$organism_filter$filter_applied),
            organism_filter_rows_before = run$organism_filter$proteins_before,
            organism_filter_rows_after = run$organism_filter$proteins_after)
        for (name in names(produced_filter)) {
            if (!isTRUE(all.equal(produced_filter[[name]], step$outcomes[[name]]))) {
                cli::cli_abort("replayed step {.val {step$state}} produced a different {.val {name}} from the record",
                    class = "workflow_replay_invalid")
            }
        }
    }
    # KIT-115 (A107-07): the declared values do not name the service's release, and a result depends on it.
    as_named <- function(service, release) sub(":$", "", paste(service, release %||% "", sep = ":"))
    recorded <- vapply(recorded_resources %||% list(), function(resource) as_named(resource$service, resource$release), "")
    unrecorded <- setdiff(sub(":$", "", setdiff(run$service, "unavailable")), recorded)
    if (length(unrecorded)) {
        cli::cli_warn(c("replayed step {.val {step$state}} queried {.val {unrecorded}}, which the record does not name",
            i = "the record names {.val {recorded}}; an enrichment result depends on the release of its service"),
            class = "workflow_replay_resource")
    }
    # KIT-173: propagated sets depend on the GO release the record names.
    if (!is.null(step$outcomes$go_db_version) || !is.null(step$outcomes$go_source_date)) {
        release <- goDbRelease()
        recorded_go <- c(as.character(step$outcomes$go_db_version %||% NA), as.character(step$outcomes$go_source_date %||% NA))
        if (!identical(recorded_go, c(release$version, release$date))) {
            cli::cli_warn(c(paste("replayed step {.val {step$state}} propagated GO terms with GO.db {.val {release$version}}",
                "({.val {release$date}}), not the recorded {.val {recorded_go[[1L]]}} ({.val {recorded_go[[2L]]}})"),
                i = "parent-term gene sets depend on the GO release"), class = "workflow_replay_resource")
        }
    }
    run
}
