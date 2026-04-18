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

