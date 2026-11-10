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
    processArgs = list(
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
    )
  )
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

  catFn("   ENRICHMENT Step: Updating R6 state with enrichment UI parameters\n")

  currentDataState <- NULL
  availableStates <- NULL
  warningMessage <- NULL

  tryCatch({
    availableStates <- diaStateHistory(workflowData$state_manager)
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

# ------------------------------------------------------------------------------------------------
# The entries (KIT-107 D3, D4, D5).
# ------------------------------------------------------------------------------------------------

#' Run the enrichment on an accepted DA manifest (KIT-107 D3)
#'
#' The module's manifest path at the start commit (`runProtEnrichAnalysisBody`, its setup, finalise and persist steps,
#' mod_prot_enrich_analysis_helpers.R:257-487, :528-603; mod_prot_enrich_state_helpers.R:371-515), in the same order,
#' without the tab's values: the caller writes them from the returned list. The manifest is read and never replaced,
#' so the caller keeps the manifest it holds (F107-01). The bus writes are today's, in today's order: `uniprot_dat_cln`
#' when the annotations are loaded, `enrichment_organism_filter`, `enrichment_analysis_results`,
#' `enrichment_ui_params`, the completed state with its step provenance (F-16, `saveWorkflowEnrichmentState`), then
#' `tab_status[[config$status_key]]`. Progress 0.2, 0.3, 0.5, 0.8 and 1.0 as today.
#' @param workflow The bus: an environment, or reactive values read inside an observer or `shiny::isolate()`.
#' @param config The omic's `enrichment` group (`omicGuiConfig(<omic>)$enrichment`).
#' @param manifest The accepted DA result manifest.
#' @param parameters The run's inputs: `selected_contrast`, `organism_taxid`, `up_cutoff`, `down_cutoff`, `q_cutoff`,
#'   `correction_method`, `enable_organism_filter`.
#' @param paths Experiment paths (`da_output_dir`, `pathway_dir`, `results_dir`, `source_dir`).
#' @param process_fn The kernel; NULL selects `config$kernel$select()` when called (the double in test mode).
#' @param notify `function(message, type, duration)`: the organism filter's warning; NULL notifies nobody.
#' @param progress `function(value, detail)`; NULL reports nothing.
#' @param refuse_missing_da TRUE refuses a run whose DA table is missing for a contrast before anything is written
#'   (class `workflow_enrichment_input`; F107-03). FALSE is the start commit's behaviour: the run continues on the
#'   missing table, with the reader's warning.
#' @return `list(results, method, contrasts, selected, raw_contrast, contrasts_tbl, pathway_dir, annotation_match,
#'   organism_filter, service, reference, parameters, saved)`: the kernel's `EnrichmentResults` as the tab holds it
#'   (the completed state carries it with the source's `@args`); the method information; every contrast's display
#'   results (with their saved tables); the selected contrast's (`resolveProtEnrichSelectedContrastResults`); the raw
#'   contrast name; the contrasts table; the pathway directory; the annotation match (NULL when none); the
#'   organism-filter metadata written to the bus; the service identity; the product reference; the 9 declared values
#'   the product's identity used; and the state save's result.
#' KIT-113 D3: the body is the helpers below, each a contiguous run of the start commit's statements, called in that
#' order; `requestWorkflowEnrichment`, `computeWorkflowEnrichment` and `acceptWorkflowEnrichment` are the same helpers
#' split at the worker.
#' @noRd
runWorkflowEnrichment <- function(workflow, config, manifest, parameters, paths, process_fn = NULL, notify = NULL,
    progress = NULL, refuse_missing_da = FALSE) {
    notify <- notify %||% function(message, type, duration) invisible(NULL)
    progress <- progress %||% function(value, detail) invisible(NULL)
    context <- workflowEnrichmentContext(workflow, config, manifest, parameters, paths, progress)
    da_results <- readWorkflowEnrichmentDa(config, context, refuse_missing_da)
    annotations <- config$annotation$resolve(workflowData = workflow, experimentPaths = paths,
        currentS4Object = context$source_object, organismTaxid = parameters$organism_taxid)$uniprotDatCln
    filtered <- workflowEnrichmentFilter(workflow, config, parameters, da_results, annotations, context$source_object,
        notify)
    progress(0.5, "Running enrichment analysis...")
    # `[[` (exact): `$organism` would partially match `organism_filter`.
    method <- config[["organism"]]$method(organismTaxid = parameters$organism_taxid,
        supportedOrganisms = config[["organism"]]$supported())$methodInfo
    declared <- workflowEnrichmentRun(config, workflow, manifest, filtered$da_results, annotations, method,
        context$source_object, parameters, context$directories, process_fn)
    publishWorkflowEnrichment(workflow, config, context, declared, method, parameters,
        filtered$matching$annotationMatchResults, filtered$organism_filter, progress)
}

#' The run's context (KIT-113 D3; the start of `runWorkflowEnrichment`): the manifest check, the source's design and
#' identifier table, progress 0.2 and 0.3, the contrasts table and the selection, the run's dependencies and the output
#' directories. It reads no assay and no DA table.
#' @inheritParams runWorkflowEnrichment
#' @return `list(manager, source_object, selected, selection, contrasts_tbl, dependencies, directories)`.
#' @noRd
workflowEnrichmentContext <- function(workflow, config, manifest, parameters, paths, progress) {
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
    list(manager = manager, source_object = source_object, selected = selected, selection = selection,
        contrasts_tbl = contrasts_tbl, dependencies = dependencies, directories = directories)
}

#' The DA tables of every contrast, with the missing-table refusal (KIT-113 D3; F107-03).
#' @param context `workflowEnrichmentContext`'s value (`contrasts_tbl`, `dependencies$designMatrix`,
#'   `directories$daOutputDir`).
#' @inheritParams runWorkflowEnrichment
#' @return The DA results for enrichment.
#' @noRd
readWorkflowEnrichmentDa <- function(config, context, refuse_missing_da) {
    da_results <- config$da_source$read(contrasts_tbl = context$contrasts_tbl,
        design_matrix = context$dependencies$designMatrix, da_output_dir = context$directories$daOutputDir)
    missing <- names(Filter(is.null, da_results@da_data))
    if (isTRUE(refuse_missing_da) && length(missing)) {
        rlang::abort(sprintf("The DA results table is missing for contrast %s.", paste(missing, collapse = ", ")),
            class = "workflow_enrichment_input")
    }
    da_results
}

#' The annotation match and the organism filter (KIT-113 D3): the match, the mapping with its warning, the filter's
#' application and the bus's filter metadata, in today's order. The managed halves call the three parts apart.
#' @param da_results The DA results for enrichment.
#' @param annotations The annotation table (`uniprot_dat_cln`), or NULL.
#' @param source_object The source's design and identifier table.
#' @inheritParams runWorkflowEnrichment
#' @return `list(da_results, matching, organism_filter)`: the DA results the kernel takes, the match, and the filter
#'   metadata written to the bus.
#' @noRd
workflowEnrichmentFilter <- function(workflow, config, parameters, da_results, annotations, source_object, notify) {
    matching <- config$annotation$match(uniprotDatCln = annotations, daResultsForEnrichment = da_results,
        currentS4Object = source_object)
    mapping <- workflowEnrichmentMapping(workflow, config, parameters, annotations, notify)
    filtered <- applyWorkflowEnrichmentFilter(config, parameters, da_results, mapping, source_object)
    list(da_results = filtered$da_results, matching = matching,
        organism_filter = persistWorkflowEnrichmentFilter(workflow, parameters, filtered))
}

#' The organism mapping from the bus when the filter is on (KIT-113 D3), else NULL; the warning when it is on and no
#' mapping is found.
#' @inheritParams workflowEnrichmentFilter
#' @return The non-empty mapping, or NULL.
#' @noRd
workflowEnrichmentMapping <- function(workflow, config, parameters, annotations, notify) {
    if (!isTRUE(parameters$enable_organism_filter) || is.null(config$organism_filter)) return(NULL)
    mapping <- config$organism_filter$mapping(workflowData = workflow, uniprotDatCln = annotations,
        targetTaxon = as.character(parameters$organism_taxid))$organismMapping
    if (!is.null(mapping) && nrow(mapping) > 0) return(mapping)
    notify("Multi-species filtering enabled but no organism mapping available. Proceeding without filtering.",
        "warning", 8)
    NULL
}

#' The organism filter's application to the DA tables (KIT-113 D3); nothing applied without a mapping.
#' @param mapping `workflowEnrichmentMapping`'s value.
#' @inheritParams workflowEnrichmentFilter
#' @return `list(da_results, applied, stats)`.
#' @noRd
applyWorkflowEnrichmentFilter <- function(config, parameters, da_results, mapping, source_object) {
    if (is.null(mapping)) {
        return(list(da_results = da_results, applied = FALSE,
            stats = list(proteins_before = 0, proteins_after = 0, proteins_removed = 0)))
    }
    filtered <- config$organism_filter$apply(daResultsForEnrichment = da_results, organismMapping = mapping,
        targetTaxon = as.character(parameters$organism_taxid), currentS4Object = source_object)
    list(da_results = filtered$daResultsForEnrichment, applied = filtered$filterApplied, stats = filtered$filterStats)
}

#' The bus's organism-filter metadata (KIT-113 D3).
#' @param filtered `applyWorkflowEnrichmentFilter`'s value (`applied`, `stats`).
#' @inheritParams runWorkflowEnrichment
#' @return The metadata written to the bus.
#' @noRd
persistWorkflowEnrichmentFilter <- function(workflow, parameters, filtered) {
    persistProtEnrichOrganismFilterMetadata(workflowData = workflow,
        organismFilterEnabled = parameters$enable_organism_filter, organismFilterApplied = filtered$applied,
        targetTaxonId = parameters$organism_taxid, filterStats = filtered$stats)
}

#' The kernel's run (KIT-113 D3): the identifier column, the kernel's arguments (validated), the kernel, the
#' option-gated debug capture, then the declaration's run and its product.
#' @param method The method information.
#' @param directories The output directories (`daOutputDir`, `pathwayDir`).
#' @param kernel The kernel; NULL selects `config$kernel$select()`.
#' @param run `function(workflow, manifest, arguments, process_fn)`: `config$run` in process; the worker's runs
#'   `computeDIAReferenceEnrichment`.
#' @inheritParams workflowEnrichmentFilter
#' @return The run's `list(results, service, reference, parameters)`.
#' @noRd
workflowEnrichmentRun <- function(config, workflow, manifest, da_results, annotations, method, source_object, parameters,
    directories, kernel = NULL, run = config$run) {
    columns <- resolveProtEnrichAnalysisInputColumns(methodInfo = method, daResultsForEnrichment = da_results,
        currentS4Object = source_object)
    arguments <- buildProtEnrichProcessEnrichmentsArgs(daResultsForEnrichment = da_results,
        organismTaxid = parameters$organism_taxid, upCutoff = parameters$up_cutoff, downCutoff = parameters$down_cutoff,
        qCutoff = parameters$q_cutoff, pathwayDir = directories$pathwayDir, goAnnotations = annotations,
        proteinIdColumn = columns$idColumn, contrastNames = names(manifest$results),
        correctionMethod = parameters$correction_method, excludeIea = config$exclude_iea)
    kernel <- kernel %||% config$kernel$select()
    .capture_checkpoint(arguments$checkpointArgs, "cp10", "enrichment_input")
    run(workflow, manifest, arguments$processArgs, process_fn = kernel)
}

#' The display and the commit (KIT-113 D3; the end of `runWorkflowEnrichment`): each contrast's display tables, the
#' bus's results payload, the results' `@args`, the UI parameters, the state commit with its provenance, `tab_status`,
#' progress 0.8 and 1.0.
#' @param context `workflowEnrichmentContext`'s value.
#' @param declared The run's `list(results, service, reference, parameters)`.
#' @param matching The annotation match results, or NULL.
#' @param organism_filter The filter metadata written to the bus.
#' @inheritParams runWorkflowEnrichment
#' @return `runWorkflowEnrichment`'s list.
#' @noRd
publishWorkflowEnrichment <- function(workflow, config, context, declared, method, parameters, matching,
    organism_filter, progress) {
    manager <- context$manager
    selected <- context$selected
    contrasts_tbl <- context$contrasts_tbl
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
    stamped <- propagateProtEnrichResultsArgs(enrichmentResults = results, currentS4Object = context$source_object,
        selectedContrast = selected, methodInfo = method, upCutoff = parameters$up_cutoff,
        downCutoff = parameters$down_cutoff, qCutoff = parameters$q_cutoff, organismTaxid = parameters$organism_taxid,
        pathwayDir = context$directories$pathwayDir)$enrichmentResults
    ui_parameters <- propagateProtEnrichUiParams(currentS4Object = context$source_object, workflowData = workflow,
        selectedContrast = selected, methodInfo = method, upCutoff = parameters$up_cutoff,
        downCutoff = parameters$down_cutoff, qCutoff = parameters$q_cutoff, organismTaxid = parameters$organism_taxid)
    updateProtEnrichStateManagerUiParams(workflowData = workflow, storedUiParams = ui_parameters$storedUiParams)
    saved <- saveWorkflowEnrichmentState(workflow, config, stamped, declared$parameters, method, parameters,
        context$directories$pathwayDir)
    status <- workflow$tab_status
    status[[config$status_key]] <- "complete"
    workflow$tab_status <- status
    progress(1.0, "Complete!")
    list(results = results, method = method, contrasts = contrast_results, selected = shown,
        raw_contrast = context$selection$rawContrastName, contrasts_tbl = contrasts_tbl,
        pathway_dir = context$directories$pathwayDir, annotation_match = matching, organism_filter = organism_filter,
        service = declared$service, reference = declared$reference, parameters = declared$parameters, saved = saved)
}

# ------------------------------------------------------------------------------------------------
# The managed halves (KIT-113 D3): the request and the accept in the controller, the compute on the worker's `omics`
# role (`computeWorkflowDiaEnrichmentAction`, utils_workflow_dia_da.R). Each is the in-process helpers above, split at
# the worker; the values that cross are saved `rds` state values (W1's P1).
# ------------------------------------------------------------------------------------------------

#' The kernels a request may name; each is looked up by name when called, never captured (KIT-105 W4).
#' @noRd
workflowEnrichmentKernels <- c("runProtEnrichDeterministicProcessEnrichments", "processEnrichments")

#' The selected kernel's name (KIT-113 D3); any other function is refused.
#' @param kernel The function `config$kernel$select()` returned.
#' @return One of `workflowEnrichmentKernels`.
#' @noRd
workflowEnrichmentKernelName <- function(kernel) {
    ns <- asNamespace("MultiScholaR")
    for (name in workflowEnrichmentKernels) if (identical(kernel, get(name, envir = ns))) return(name)
    rlang::abort("The enrichment kernel is neither the deterministic test double nor processEnrichments.",
        class = "workflow_enrichment_input")
}

#' The kernel a request names (KIT-113 D3), looked up when called; any other name is refused.
#' @param name One of `workflowEnrichmentKernels`.
#' @return The kernel.
#' @noRd
workflowEnrichmentKernel <- function(name) {
    if (!(is.character(name) && length(name) == 1L && name %in% workflowEnrichmentKernels)) {
        rlang::abort("The enrichment request names an unknown kernel.", class = "workflow_enrichment_input")
    }
    get(name, envir = asNamespace("MultiScholaR"))
}

#' Save a transfer value as an `rds` state value (KIT-113 D3, W1's P1); NULL stays NULL.
#' @noRd
saveWorkflowEnrichmentValue <- function(manager, value, label) {
    if (is.null(value)) return(NULL)
    saveWorkflowStateValue(manager$.context, value, label, codec = "rds")
}

#' Read a transfer value saved by `saveWorkflowEnrichmentValue`; NULL stays NULL.
#' @noRd
readWorkflowEnrichmentValue <- function(manager, reference) {
    if (is.null(reference)) return(NULL)
    resolveWorkflowStateValue(manager$.context, reference)
}

#' The managed request half (KIT-113 D3), in the controller: the context, the annotation resolve, the mapping and its
#' warning, the method and the kernel's name. The manifest must be the head's `da_results`. No DA table or assay is
#' read.
#' @inheritParams runWorkflowEnrichment
#' @return `list(context, parameters)`: the context the accept half takes (with the method and the run's inputs), and
#'   the request's parameters (the manifest's reference, the contrasts table, the run's inputs, the two directories,
#'   the annotation and mapping references, the method, the kernel's name).
#' @noRd
requestWorkflowEnrichment <- function(workflow, config, manifest, parameters, paths, notify = NULL, progress = NULL) {
    notify <- notify %||% function(message, type, duration) invisible(NULL)
    progress <- progress %||% function(value, detail) invisible(NULL)
    context <- workflowEnrichmentContext(workflow, config, manifest, parameters, paths, progress)
    annotations <- config$annotation$resolve(workflowData = workflow, experimentPaths = paths,
        currentS4Object = context$source_object, organismTaxid = parameters$organism_taxid)$uniprotDatCln
    mapping <- workflowEnrichmentMapping(workflow, config, parameters, annotations, notify)
    progress(0.5, "Running enrichment analysis...")
    method <- config[["organism"]]$method(organismTaxid = parameters$organism_taxid,
        supportedOrganisms = config[["organism"]]$supported())$methodInfo
    kernel <- workflowEnrichmentKernelName(config$kernel$select())
    manager <- context$manager
    held <- manager$.context$snapshot$da_results
    if (is.null(held) || !identical(resolveWorkflowStateValue(manager$.context, held), manifest)) {
        rlang::abort("The DA results the tab holds are not the current state's.", class = "workflow_enrichment_input")
    }
    list(context = c(context, list(method = method, parameters = parameters)),
        parameters = list(manifest = held, contrasts_tbl = context$contrasts_tbl, inputs = parameters,
            directories = context$directories[c("daOutputDir", "pathwayDir")],
            annotations = saveWorkflowEnrichmentValue(manager, annotations, "enrichment-annotations"),
            mapping = saveWorkflowEnrichmentValue(manager, mapping, "enrichment-organism-mapping"),
            method = method, kernel = kernel))
}

#' The managed compute half (KIT-113 D3), worker-safe: the source's components, the DA tables with the refusal, the
#' match, the filter's application and the run with its product. The kernel writes into the pathway directory, as in
#' process. The results are made portable (F107-05) and saved as an `rds` state value.
#' @param manager The worker's reopened owner.
#' @param config The proteomics `enrichment` group.
#' @param parameters `requestWorkflowEnrichment`'s parameters.
#' @return `list(reference, results, service, parameters, matching, filter)`.
#' @noRd
computeWorkflowEnrichment <- function(manager, config, parameters) {
    manifest <- resolveWorkflowStateValue(manager$.context, parameters$manifest)
    workflow <- list(state_manager = manager)
    source_object <- readDIAExecutionComponents(manager, manifest$source, fields = config$da_source$components)
    dependencies <- resolveProtEnrichRunDependencies(source_object, manifest$results, workflow, parameters$contrasts_tbl)
    da_results <- readWorkflowEnrichmentDa(config, list(contrasts_tbl = dependencies$contrastsTbl,
        dependencies = dependencies, directories = parameters$directories), refuse_missing_da = TRUE)
    annotations <- readWorkflowEnrichmentValue(manager, parameters$annotations)
    matching <- config$annotation$match(uniprotDatCln = annotations, daResultsForEnrichment = da_results,
        currentS4Object = source_object)
    filtered <- applyWorkflowEnrichmentFilter(config, parameters$inputs, da_results,
        readWorkflowEnrichmentValue(manager, parameters$mapping), source_object)
    declared <- workflowEnrichmentRun(config, workflow, manifest, filtered$da_results, annotations, parameters$method,
        source_object, parameters$inputs, parameters$directories, workflowEnrichmentKernel(parameters$kernel),
        run = function(workflow, manifest, arguments, process_fn)
            computeDIAReferenceEnrichment(workflow$state_manager, manifest, arguments, process_fn))
    list(reference = declared$reference,
        results = saveWorkflowEnrichmentValue(manager, portableWorkflowEnrichmentResults(declared$results),
            "enrichment-results"),
        service = declared$service, parameters = declared$parameters,
        # F113-02: the match results scale with the DA proteins, so they travel by reference, as the results do.
        matching = saveWorkflowEnrichmentValue(manager, matching$annotationMatchResults, "enrichment-annotation-match"),
        filter = filtered[c("applied", "stats")])
}

#' The managed accept half (KIT-113 D3), in the controller, inside the owner's scope: the completion's provenance
#' cleared, the results hydrated, the external-resource records, the bus's filter metadata, then the display and the
#' commit as in process.
#' @param context `requestWorkflowEnrichment`'s context.
#' @param delivered The managed result (`readManagedOmicsResult`): `reference` and the handler's `display`.
#' @inheritParams runWorkflowEnrichment
#' @return `runWorkflowEnrichment`'s list.
#' @noRd
acceptWorkflowEnrichment <- function(workflow, config, context, delivered, progress = NULL) {
    progress <- progress %||% function(value, detail) invisible(NULL)
    manager <- context$manager
    withDIAStateScope(manager, function() {
        clearWorkflowStepProvenance(manager)
        display <- delivered$display
        declared <- list(results = readWorkflowEnrichmentValue(manager, display$results), service = display$service,
            reference = delivered$reference, parameters = display$parameters)
        recordDIAEnrichmentResources(manager, declared$service, declared$parameters$results,
            declared$parameters$taxon_id)
        organism_filter <- persistWorkflowEnrichmentFilter(workflow, context$parameters, display$filter)
        publishWorkflowEnrichment(workflow, config, context, declared, context$method, context$parameters,
            readWorkflowEnrichmentValue(manager, display$matching), organism_filter, progress)
    })
}

#' Save the completed enrichment state with its step provenance (KIT-107 D3, D4; F-16)
#'
#' Today's `saveProtEnrichCompletedState` (mod_prot_enrich_state_helpers.R:210-262 at the start commit): the same state
#' name, config object, description and failure warning. Immediately before the commit it sets the step's provenance as
#' the protein steps set theirs (utils_workflow_prot_declarations.R:1158-1164): the registry entry's operation, its
#' default variant, its declared parameter names and fixed values, and the declared values the product's identity used.
#' A failed commit clears the pending provenance, so the next commit is not attributed to the enrichment.
#' @param workflow The bus (its `state_manager`).
#' @param config The omic's `enrichment` group (`operation`, `state_name`).
#' @param results The `EnrichmentResults` to commit.
#' @param declared The declared values (`runDIAReferenceEnrichment`'s `parameters`).
#' @param method The method information (`method`, `supported`).
#' @param parameters The run's inputs (`selected_contrast`, cutoffs, `organism_taxid`).
#' @param pathway_dir The pathway directory.
#' @return `list(attempted, saved, warning)`.
#' @noRd
saveWorkflowEnrichmentState <- function(workflow, config, results, declared, method, parameters, pathway_dir) {
    cat("   ENRICHMENT Step: Saving results to R6 state manager...\n")
    manager <- workflow$state_manager
    warning_message <- NULL
    saved <- FALSE
    tryCatch({
        results <- portableWorkflowEnrichmentResults(results) # F107-05
        entry <- workflowOperationEntry(config$operation()$name)
        setWorkflowStepProvenance(manager, entry$name, entry$default_variant, entry$parameters, entry$declaration$fixed,
            values = declared)
        saveDIAOwnerState(manager,
            state_name = config$state_name,
            s4_data_object = results,
            config_object = list(
                selected_contrast = parameters$selected_contrast,
                analysis_method = method$method,
                organism_supported = method$supported,
                up_cutoff = parameters$up_cutoff,
                down_cutoff = parameters$down_cutoff,
                q_cutoff = parameters$q_cutoff,
                organism_taxid = parameters$organism_taxid,
                pathway_dir = pathway_dir),
            description = paste("Enrichment analysis completed using", method$method, "for contrast:",
                parameters$selected_contrast))
        saved <- TRUE
        cat("   ENRICHMENT Step: Successfully saved state to R6 state manager\n")
    }, error = function(e) {
        clearWorkflowStepProvenance(manager)
        warning_message <<- e$message
        cat(sprintf("   ENRICHMENT Step: Warning saving to R6 state manager: %s\n", e$message))
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
#'   `list(all_enrichment_results, enrichment_results_full, method)`.
#' @noRd
restoreWorkflowEnrichment <- function(manager, config) {
    if (is.null(manager) || !config$state_name %in% unlist(diaStateHistory(manager))) return(NULL)
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
    list(all_enrichment_results = contrast_results, enrichment_results_full = results, method = method)
}

#' Replay a recorded enrichment step after the replayed DA (KIT-107 D5; F-16)
#'
#' 1. The replayed manifest's DA tables are published into `<scratch>/enrichment/{da,graphs}` through the DA
#'    publication body the GUI uses (`produceDIADaPublication`, assets not saved), with the annotation table the GUI's
#'    DA would have had: the prepared project's `<source_dir>/uniprot_dat_cln.RDS`, else none.
#' 2. `runWorkflowEnrichment` runs on a bus holding the owner and the record's DA contrasts, with the record's taxon,
#'    cutoffs and correction, the organism filter off, the first contrast's friendly name selected (a display choice)
#'    and the kernel the GUI would select. Its annotation resolver takes the GUI's order with an empty bus: the
#'    prepared project's `uniprot_dat_cln.RDS`, then the test-mode table, then the network.
#' 3. The run commits the completed state with the same provenance (D4). The declared values it produced must equal
#'    the record's 9 (`protein_id_column` and `contrast_names` among them); a difference refuses the replay, naming
#'    the parameter. `results` is compared by its contrasts: its values are the DA tables' entity ids, which are
#'    store-local (a saved entity's id hashes its payload's file name, `saveWorkflowTable`), so the replay store's ids
#'    never equal the origin's; the tables' content is the replayed DA's (proposed F107-04).
#'
#' Limits: a live g:Profiler result reproduces only when the service's database version is the recorded one (the
#' record's `external_resources`); the annotation table is the prepared project's and is not hash-checked; the organism
#' filter is not in the record (P-28), so the step replays unfiltered.
#' @param manager The replay's owner.
#' @param step The record's enrichment step.
#' @param parameters The step's declared values, decoded (`workflowReplayParameters`).
#' @param manifest The replayed DA's manifest.
#' @param contrasts The record's DA contrasts table.
#' @param paths Experiment paths of the prepared project.
#' @param scratch Directory for the publication and the pathway files.
#' @return The run's result (`runWorkflowEnrichment`).
#' @noRd
replayWorkflowEnrichmentStep <- function(manager, step, parameters, manifest, contrasts, paths, scratch) {
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
    run <- runWorkflowEnrichment(workflow, config, manifest, list(selected_contrast = contrasts$friendly_names[[1L]],
        organism_taxid = format(parameters$taxon_id, scientific = FALSE, trim = TRUE), up_cutoff = parameters$up_cutoff,
        down_cutoff = parameters$down_cutoff, q_cutoff = parameters$q_cutoff,
        correction_method = parameters$correction_method, enable_organism_filter = FALSE), run_paths)
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
    run
}
