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
  taxonId <- as.numeric(organismTaxid)

  list(
    checkpointArgs = list(
      da_results_s4 = daResultsForEnrichment,
      taxon_id = taxonId,
      up_cutoff = upCutoff,
      down_cutoff = downCutoff,
      q_cutoff = qCutoff,
      pathway_dir = pathwayDir,
      go_annotations = goAnnotations,
      exclude_iea = excludeIea,
      protein_id_column = proteinIdColumn,
      contrast_names = contrastNames,
      correction_method = correctionMethod
    ),
    processArgs = list(
      da_results = daResultsForEnrichment,
      taxon_id = taxonId,
      up_cutoff = upCutoff,
      down_cutoff = downCutoff,
      q_cutoff = qCutoff,
      pathway_dir = pathwayDir,
      go_annotations = goAnnotations,
      exclude_iea = excludeIea,
      protein_id_column = proteinIdColumn,
      contrast_names = contrastNames,
      correction_method = correctionMethod
    )
  )
}

prepareProtEnrichProcessExecution <- function(input,
                                              enrichmentData,
                                              daResultsForEnrichment,
                                              pathwayDir,
                                              goAnnotations,
                                              currentAnalysisMethodFn,
                                              resolveAnalysisInputColumnsFn = resolveProtEnrichAnalysisInputColumns,
                                              buildProcessEnrichmentsArgsFn = buildProtEnrichProcessEnrichmentsArgs,
                                              catFn = cat) {
  methodInfo <- currentAnalysisMethodFn()
  catFn(sprintf("   ENRICHMENT Step: Using analysis method: %s\n", methodInfo$method))

  inputColumnConfig <- resolveAnalysisInputColumnsFn(
    methodInfo = methodInfo,
    daResultsForEnrichment = daResultsForEnrichment,
    currentS4Object = enrichmentData$current_s4_object
  )

  enrichmentArgs <- buildProcessEnrichmentsArgsFn(
    daResultsForEnrichment = daResultsForEnrichment,
    organismTaxid = input$organism_taxid,
    upCutoff = input$up_cutoff,
    downCutoff = input$down_cutoff,
    qCutoff = input$q_cutoff,
    pathwayDir = pathwayDir,
    goAnnotations = goAnnotations,
    proteinIdColumn = inputColumnConfig$idColumn,
    contrastNames = names(enrichmentData$da_results_data),
    correctionMethod = input$correction_method
  )

  list(
    methodInfo = methodInfo,
    inputColumnConfig = inputColumnConfig,
    enrichmentArgs = enrichmentArgs
  )
}

executeProtEnrichProcessEnrichments <- function(enrichmentArgs,
                                                upCutoff,
                                                downCutoff,
                                                qCutoff,
                                                captureCheckpointFn = .capture_checkpoint,
                                                processEnrichmentsFn = processEnrichments,
                                                catFn = cat) {
  captureCheckpointFn(enrichmentArgs$checkpointArgs, "cp10", "enrichment_input")

  enrichmentResults <- do.call(processEnrichmentsFn, enrichmentArgs$processArgs)

  catFn(sprintf(
    "   ENRICHMENT Step: processEnrichments completed with up_cutoff: %.1f, down_cutoff: %.1f, q_cutoff: %.3f\n",
    upCutoff,
    downCutoff,
    qCutoff
  ))

  enrichmentResults
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
          if (!is.null(upResults) &&
              isTRUE(isEnrichResultFn(upResults, "enrichResult")) &&
              nrow(upResults@result) > 0) {
            upDf <- upResults@result
            upDf$directionality <- "up"
            upDf$analysis_method <- "clusterprofiler"
            clusterprofilerResults <- rbind(clusterprofilerResults, upDf)
          }
        }

        if (!is.null(contrastEnrichment$down)) {
          downResults <- contrastEnrichment$down
          if (!is.null(downResults) &&
              isTRUE(isEnrichResultFn(downResults, "enrichResult")) &&
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

