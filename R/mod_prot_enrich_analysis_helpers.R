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

prepareProtEnrichAnalysisBodySetup <- function(selectedContrast,
                                               input,
                                               enrichmentData,
                                               workflowData,
                                               experimentPaths,
                                               resolveSelectedDaResultsFn = resolveProtEnrichSelectedDaResults,
                                               resolveRunDependenciesFn = resolveProtEnrichRunDependencies,
                                               resolveOutputDirectoriesFn = resolveProtEnrichOutputDirectories,
                                               createDaResultsForEnrichmentFn = createDAResultsForEnrichment,
                                               resolveUniprotAnnotationsFn = resolveProtEnrichUniprotAnnotations,
                                               resolveAnnotationMatchingFn = resolveProtEnrichAnnotationMatching,
                                               resolveOrganismMappingFn = resolveProtEnrichOrganismMapping,
                                               applyOrganismFilterFn = applyProtEnrichOrganismFilter,
                                               persistOrganismFilterMetadataFn = persistProtEnrichOrganismFilterMetadata,
                                               showNotificationFn = shiny::showNotification,
                                               globalEnv = .GlobalEnv,
                                               catFn = cat) {
  catFn(sprintf("   ENRICHMENT Step: Selected contrast (friendly name) = %s\n", selectedContrast))
  catFn(sprintf(
    "   ENRICHMENT Step: Available DE results: %s\n",
    paste(names(enrichmentData$da_results_data), collapse = ", ")
  ))

  contrastsTbl <- if (exists("contrasts_tbl", envir = globalEnv)) {
    get("contrasts_tbl", envir = globalEnv)
  } else {
    NULL
  }

  selectedDaConfig <- resolveSelectedDaResultsFn(
    selectedContrast = selectedContrast,
    daResultsData = enrichmentData$da_results_data,
    contrastsTbl = contrastsTbl
  )
  selectedDaResults <- selectedDaConfig$selectedDaResults
  rawContrastName <- selectedDaConfig$rawContrastName
  availableKeys <- selectedDaConfig$availableKeys

  if (!is.null(selectedDaConfig$mappedRawContrastName)) {
    catFn(sprintf(
      "   ENRICHMENT Step: Mapped friendly name '%s' to raw name '%s'\n",
      selectedContrast,
      selectedDaConfig$mappedRawContrastName
    ))
  }

  if (identical(selectedDaConfig$source, "friendly_name")) {
    catFn(sprintf("   ENRICHMENT Step: Found DE results for raw contrast %s\n", rawContrastName))
  } else if (!is.null(selectedDaConfig$mappedRawContrastName)) {
    catFn(sprintf(
      "   ENRICHMENT Step: No DE results found for raw contrast '%s'\n",
      selectedDaConfig$mappedRawContrastName
    ))
  }

  if (is.null(selectedDaResults)) {
    catFn(sprintf(
      "   ENRICHMENT Step: Available DE result keys: %s\n",
      paste(availableKeys, collapse = ", ")
    ))
    catFn(sprintf("   ENRICHMENT Step: Looking for friendly name: %s\n", selectedContrast))
  }

  if (identical(selectedDaConfig$source, "fuzzy_match")) {
    catFn(sprintf(
      "   ENRICHMENT Step: Found matching DE results using fuzzy matching: %s\n",
      rawContrastName
    ))
  }

  if (identical(selectedDaConfig$source, "direct_key")) {
    catFn(sprintf(
      "   ENRICHMENT Step: Found DE results using direct key lookup: %s\n",
      selectedContrast
    ))
  }

  if (is.null(selectedDaResults)) {
    stop(sprintf(
      "Could not find DE results for contrast '%s'. Available contrasts: %s",
      selectedContrast,
      paste(availableKeys, collapse = ", ")
    ))
  }

  hadCurrentS4 <- !is.null(enrichmentData$current_s4_object)
  dependencyConfig <- resolveRunDependenciesFn(
    currentS4Object = enrichmentData$current_s4_object,
    daResultsData = enrichmentData$da_results_data,
    workflowData = workflowData,
    contrastsTbl = contrastsTbl
  )
  contrastsTbl <- dependencyConfig$contrastsTbl
  designMatrix <- dependencyConfig$designMatrix
  enrichmentData$current_s4_object <- dependencyConfig$currentS4

  if (!is.null(contrastsTbl)) {
    catFn("   ENRICHMENT Step: Found contrasts_tbl in global environment\n")
  }

  if (!hadCurrentS4 && !is.null(enrichmentData$current_s4_object)) {
    catFn("   ENRICHMENT Step: S4 object is NULL, trying to retrieve it...\n")

    if (identical(dependencyConfig$s4Source, "da_results_first_result")) {
      catFn(sprintf(
        "   ENRICHMENT Step: Got S4 from DE results theObject (class: %s)\n",
        class(enrichmentData$current_s4_object)
      ))
    }

    if (identical(dependencyConfig$s4Source, "state_manager")) {
      catFn(sprintf(
        "   ENRICHMENT Step: Got S4 from state manager (class: %s)\n",
        class(enrichmentData$current_s4_object)
      ))
    }
  }

  if (!is.null(dependencyConfig$designMatrixError)) {
    catFn(sprintf(
      "   ENRICHMENT Step: Could not access design_matrix slot: %s\n",
      dependencyConfig$designMatrixError
    ))
  }

  if (identical(dependencyConfig$designMatrixSource, "s4_object")) {
    catFn("   ENRICHMENT Step: Found design_matrix in S4 object\n")
  }

  if (identical(dependencyConfig$designMatrixSource, "global_environment")) {
    catFn("   ENRICHMENT Step: Found design_matrix in global environment\n")
  }

  if (is.null(contrastsTbl) || is.null(designMatrix)) {
    stop("Missing contrasts_tbl or design_matrix required for enrichment analysis. Please ensure DE analysis was completed successfully.")
  }

  pathConfig <- resolveOutputDirectoriesFn(
    experimentPaths = experimentPaths
  )
  daOutputDir <- pathConfig$daOutputDir
  pathwayDir <- pathConfig$pathwayDir

  if (identical(pathConfig$daOutputDirSource, "experiment_paths")) {
    catFn(sprintf("   ENRICHMENT Step: Using da_output_dir from experiment_paths: %s\n", daOutputDir))
  } else {
    catFn(sprintf("   ENRICHMENT Step: Falling back to da_proteins folder: %s\n", daOutputDir))
  }

  if (identical(pathConfig$pathwayDirSource, "experiment_paths")) {
    catFn(sprintf("   ENRICHMENT Step: Using pathway_dir from experiment_paths: %s\n", pathwayDir))
  } else {
    catFn(sprintf("   ENRICHMENT Step: Created fallback pathway directory: %s\n", pathwayDir))
  }

  daResultsForEnrichment <- createDaResultsForEnrichmentFn(
    contrasts_tbl = contrastsTbl,
    design_matrix = designMatrix,
    da_output_dir = daOutputDir
  )

  catFn("   ENRICHMENT Step: S4 object created successfully\n")

  uniprotConfig <- resolveUniprotAnnotationsFn(
    workflowData = workflowData,
    experimentPaths = experimentPaths,
    currentS4Object = enrichmentData$current_s4_object,
    organismTaxid = input$organism_taxid
  )
  uniprotDatCln <- uniprotConfig$uniprotDatCln

  annotationMatchConfig <- resolveAnnotationMatchingFn(
    uniprotDatCln = uniprotDatCln,
    daResultsForEnrichment = daResultsForEnrichment,
    currentS4Object = enrichmentData$current_s4_object
  )

  if (!is.null(annotationMatchConfig$annotationMatchResults)) {
    enrichmentData$annotation_match_results <- annotationMatchConfig$annotationMatchResults
  }

  organismFilterApplied <- FALSE
  filterStats <- list(proteins_before = 0, proteins_after = 0, proteins_removed = 0)

  if (isTRUE(input$enable_organism_filter)) {
    catFn("*** ENRICHMENT Step: Multi-species filtering ENABLED ***\n")

    targetTaxon <- as.character(input$organism_taxid)
    organismMappingConfig <- resolveOrganismMappingFn(
      workflowData = workflowData,
      uniprotDatCln = uniprotDatCln,
      targetTaxon = targetTaxon
    )
    organismMapping <- organismMappingConfig$organismMapping

    if (!is.null(organismMapping) && nrow(organismMapping) > 0) {
      organismFilterConfig <- applyOrganismFilterFn(
        daResultsForEnrichment = daResultsForEnrichment,
        organismMapping = organismMapping,
        targetTaxon = targetTaxon,
        currentS4Object = enrichmentData$current_s4_object
      )
      daResultsForEnrichment <- organismFilterConfig$daResultsForEnrichment
      organismFilterApplied <- organismFilterConfig$filterApplied
      filterStats <- organismFilterConfig$filterStats
    } else {
      catFn("   ENRICHMENT Step: WARNING - No organism mapping available, filtering skipped\n")
      showNotificationFn(
        "Multi-species filtering enabled but no organism mapping available. Proceeding without filtering.",
        type = "warning",
        duration = 8
      )
    }
  } else {
    catFn("   ENRICHMENT Step: Multi-species filtering DISABLED\n")
  }

  persistOrganismFilterMetadataFn(
    workflowData = workflowData,
    organismFilterEnabled = input$enable_organism_filter,
    organismFilterApplied = organismFilterApplied,
    targetTaxonId = input$organism_taxid,
    filterStats = filterStats
  )

  list(
    rawContrastName = rawContrastName,
    contrastsTbl = contrastsTbl,
    pathwayDir = pathwayDir,
    daResultsForEnrichment = daResultsForEnrichment,
    goAnnotations = uniprotDatCln,
    organismFilterApplied = organismFilterApplied,
    filterStats = filterStats
  )
}

runProtEnrichAnalysisBody <- function(input,
                                      enrichmentData,
                                      workflowData,
                                      experimentPaths,
                                      currentAnalysisMethodFn,
                                      prepareAnalysisSetupFn = prepareProtEnrichAnalysisBodySetup,
                                      captureCheckpointFn = .capture_checkpoint,
                                      resolveSelectedDaResultsFn = resolveProtEnrichSelectedDaResults,
                                      resolveRunDependenciesFn = resolveProtEnrichRunDependencies,
                                      resolveOutputDirectoriesFn = resolveProtEnrichOutputDirectories,
                                      createDaResultsForEnrichmentFn = createDAResultsForEnrichment,
                                      resolveUniprotAnnotationsFn = resolveProtEnrichUniprotAnnotations,
                                      resolveAnnotationMatchingFn = resolveProtEnrichAnnotationMatching,
                                      resolveOrganismMappingFn = resolveProtEnrichOrganismMapping,
                                      applyOrganismFilterFn = applyProtEnrichOrganismFilter,
                                      persistOrganismFilterMetadataFn = persistProtEnrichOrganismFilterMetadata,
                                      resolveAnalysisInputColumnsFn = resolveProtEnrichAnalysisInputColumns,
                                      buildProcessEnrichmentsArgsFn = buildProtEnrichProcessEnrichmentsArgs,
                                      prepareProcessExecutionFn = prepareProtEnrichProcessExecution,
                                      processEnrichmentsFn = processEnrichments,
                                      resolveProcessEnrichmentsFn = resolveProtEnrichProcessEnrichmentsFn,
                                      executeProcessEnrichmentsFn = executeProtEnrichProcessEnrichments,
                                      buildAllContrastResultsFn = buildProtEnrichAllContrastResults,
                                      resolveSelectedContrastResultsFn = resolveProtEnrichSelectedContrastResults,
                                      capturePostProcessResultsFn = captureProtEnrichPostProcessResults,
                                      persistAnalysisResultsFn = persistProtEnrichAnalysisResults,
                                      buildAnalysisResultsPayloadFn = buildProtEnrichAnalysisResultsPayload,
                                      propagateResultsArgsFn = propagateProtEnrichResultsArgs,
                                      propagateUiParamsFn = propagateProtEnrichUiParams,
                                      updateStateManagerUiParamsFn = updateProtEnrichStateManagerUiParams,
                                      saveCompletedStateFn = saveProtEnrichCompletedState,
                                      completeTabStatusFn = completeProtEnrichTabStatus,
                                      completeProgressFn = completeProtEnrichProgress,
                                      finalizeAnalysisResultsFn = finalizeProtEnrichAnalysisBodyResults,
                                      showNotificationFn = shiny::showNotification,
                                      incProgressFn = shiny::incProgress,
                                      globalEnv = .GlobalEnv,
                                      catFn = cat) {
  selectedContrast <- input$selected_contrast

  incProgressFn(0.2, detail = "Transforming DE data...")

  incProgressFn(0.3, detail = "Creating DE results S4 object...")

  catFn("   ENRICHMENT Step: Creating DA results S4 object using createDAResultsForEnrichment\n")

  setupConfig <- prepareAnalysisSetupFn(
    selectedContrast = selectedContrast,
    input = input,
    enrichmentData = enrichmentData,
    workflowData = workflowData,
    experimentPaths = experimentPaths,
    resolveSelectedDaResultsFn = resolveSelectedDaResultsFn,
    resolveRunDependenciesFn = resolveRunDependenciesFn,
    resolveOutputDirectoriesFn = resolveOutputDirectoriesFn,
    createDaResultsForEnrichmentFn = createDaResultsForEnrichmentFn,
    resolveUniprotAnnotationsFn = resolveUniprotAnnotationsFn,
    resolveAnnotationMatchingFn = resolveAnnotationMatchingFn,
    resolveOrganismMappingFn = resolveOrganismMappingFn,
    applyOrganismFilterFn = applyOrganismFilterFn,
    persistOrganismFilterMetadataFn = persistOrganismFilterMetadataFn,
    showNotificationFn = showNotificationFn,
    globalEnv = globalEnv,
    catFn = catFn
  )
  rawContrastName <- setupConfig$rawContrastName
  contrastsTbl <- setupConfig$contrastsTbl
  pathwayDir <- setupConfig$pathwayDir
  daResultsForEnrichment <- setupConfig$daResultsForEnrichment
  organismFilterApplied <- setupConfig$organismFilterApplied
  filterStats <- setupConfig$filterStats

  incProgressFn(0.5, detail = "Running enrichment analysis...")

  catFn("   ENRICHMENT Step: Running processEnrichments\n")

  processExecutionConfig <- prepareProcessExecutionFn(
    input = input,
    enrichmentData = enrichmentData,
    daResultsForEnrichment = daResultsForEnrichment,
    pathwayDir = pathwayDir,
    goAnnotations = setupConfig$goAnnotations,
    currentAnalysisMethodFn = currentAnalysisMethodFn,
    resolveAnalysisInputColumnsFn = resolveAnalysisInputColumnsFn,
    buildProcessEnrichmentsArgsFn = buildProcessEnrichmentsArgsFn,
    catFn = catFn
  )
  methodInfo <- processExecutionConfig$methodInfo
  enrichmentArgs <- processExecutionConfig$enrichmentArgs
  selectedProcessEnrichmentsFn <- resolveProcessEnrichmentsFn(
    processEnrichmentsFn = processEnrichmentsFn
  )

  enrichmentResults <- executeProcessEnrichmentsFn(
    enrichmentArgs = enrichmentArgs,
    upCutoff = input$up_cutoff,
    downCutoff = input$down_cutoff,
    qCutoff = input$q_cutoff,
    captureCheckpointFn = captureCheckpointFn,
    processEnrichmentsFn = selectedProcessEnrichmentsFn,
    catFn = catFn
  )

  finalizeAnalysisResultsFn(
    selectedContrast = selectedContrast,
    rawContrastName = rawContrastName,
    organismFilterApplied = organismFilterApplied,
    filterStats = filterStats,
    enrichmentResults = enrichmentResults,
    enrichmentData = enrichmentData,
    workflowData = workflowData,
    input = input,
    methodInfo = methodInfo,
    contrastsTbl = contrastsTbl,
    pathwayDir = pathwayDir,
    buildAllContrastResultsFn = buildAllContrastResultsFn,
    resolveSelectedContrastResultsFn = resolveSelectedContrastResultsFn,
    capturePostProcessResultsFn = capturePostProcessResultsFn,
    persistAnalysisResultsFn = persistAnalysisResultsFn,
    buildAnalysisResultsPayloadFn = buildAnalysisResultsPayloadFn,
    propagateResultsArgsFn = propagateResultsArgsFn,
    propagateUiParamsFn = propagateUiParamsFn,
    updateStateManagerUiParamsFn = updateStateManagerUiParamsFn,
    saveCompletedStateFn = saveCompletedStateFn,
    completeTabStatusFn = completeTabStatusFn,
    completeProgressFn = completeProgressFn,
    incProgressFn = incProgressFn,
    catFn = catFn
  )
}
