resolveProtEnrichOrganismMapping <- function(workflowData,
                                             uniprotDatCln,
                                             targetTaxon,
                                             catFn = cat) {
  organismMapping <- NULL
  mappingSource <- NULL
  accessionColumn <- NULL
  taxonColumn <- NULL
  warningMessage <- NULL
  availableTaxonIds <- character()

  if (!is.null(workflowData$mixed_species_analysis) &&
      !is.null(workflowData$mixed_species_analysis$organism_mapping)) {
    organismMapping <- workflowData$mixed_species_analysis$organism_mapping
    mappingSource <- "workflow_data"
    catFn("   ENRICHMENT Step: Using organism_mapping from import module\n")
  }

  if (is.null(organismMapping) && !is.null(uniprotDatCln)) {
    catFn(sprintf(
      "   ENRICHMENT Step: Checking uniprot_dat_cln columns: %s\n",
      paste(names(uniprotDatCln), collapse = ", ")
    ))

    organismNameToTaxid <- c(
      "Homo sapiens" = "9606", "Human" = "9606",
      "Mus musculus" = "10090", "Mouse" = "10090",
      "Rattus norvegicus" = "10116", "Rat" = "10116",
      "Drosophila melanogaster" = "7227", "Fruit fly" = "7227",
      "Caenorhabditis elegans" = "6239",
      "Saccharomyces cerevisiae" = "4932", "Yeast" = "4932",
      "Arabidopsis thaliana" = "3702",
      "Danio rerio" = "7955", "Zebrafish" = "7955",
      "Gallus gallus" = "9031", "Chicken" = "9031",
      "Sus scrofa" = "9823", "Pig" = "9823",
      "Bos taurus" = "9913", "Bovine" = "9913", "Cow" = "9913"
    )
    possibleAccCols <- c(
      "Entry", "entry", "UniProt_Acc", "uniprot_acc", "Accession", "accession", "protein_id"
    )

    for (col in possibleAccCols) {
      if (col %in% names(uniprotDatCln)) {
        accessionColumn <- col
        break
      }
    }

    if ("Organism" %in% names(uniprotDatCln) && !is.null(accessionColumn)) {
      catFn("   ENRICHMENT Step: Found 'Organism' column with organism names - mapping to taxon IDs\n")

      mapOrganismToTaxid <- function(orgName) {
        if (is.na(orgName) || orgName == "") {
          return(NA_character_)
        }

        for (name in names(organismNameToTaxid)) {
          if (grepl(name, orgName, ignore.case = TRUE)) {
            return(unname(organismNameToTaxid[[name]]))
          }
        }

        NA_character_
      }

      organismMapping <- tryCatch({
        uniprotDatCln |>
          dplyr::select(uniprot_acc = dplyr::all_of(accessionColumn), organism_name = Organism) |>
          dplyr::mutate(
            taxon_id = unname(vapply(organism_name, mapOrganismToTaxid, character(1)))
          ) |>
          dplyr::select(uniprot_acc, taxon_id)
      }, error = function(e) {
        warningMessage <<- e$message
        catFn(sprintf("   ENRICHMENT Step: Error creating organism_mapping from names: %s\n", e$message))
        NULL
      })

      if (!is.null(organismMapping)) {
        mappingSource <- "organism_names"
        availableTaxonIds <- unique(organismMapping$taxon_id)
        catFn(sprintf(
          "   ENRICHMENT Step: Created organism_mapping from organism names (%d entries)\n",
          nrow(organismMapping)
        ))
        catFn(sprintf(
          "   ENRICHMENT Step: Unique taxon IDs mapped: %s\n",
          paste(availableTaxonIds[!is.na(availableTaxonIds)], collapse = ", ")
        ))

        if (targetTaxon %in% availableTaxonIds) {
          targetCount <- sum(organismMapping$taxon_id == targetTaxon, na.rm = TRUE)
          catFn(sprintf(
            "   ENRICHMENT Step: Found %d proteins matching target taxon %s\n",
            targetCount,
            targetTaxon
          ))
        } else {
          catFn(sprintf(
            "   ENRICHMENT Step: WARNING - Target taxon %s not found in mapped taxon IDs!\n",
            targetTaxon
          ))
          catFn(sprintf(
            "   ENRICHMENT Step: Available taxon IDs: %s\n",
            paste(availableTaxonIds[!is.na(availableTaxonIds)], collapse = ", ")
          ))
        }
      }
    } else {
      possibleTaxonCols <- c(
        "Organism (ID)", "organism_id", "Organism_ID", "taxon_id", "Taxon_ID",
        "Taxonomy ID", "taxonomy_id", "NCBI_TaxID", "ncbi_taxid", "OX"
      )

      for (col in possibleTaxonCols) {
        if (col %in% names(uniprotDatCln)) {
          taxonColumn <- col
          catFn(sprintf("   ENRICHMENT Step: Found taxon column: %s\n", taxonColumn))
          break
        }
      }

      if (!is.null(taxonColumn) && !is.null(accessionColumn)) {
        organismMapping <- tryCatch({
          uniprotDatCln |>
            dplyr::select(
              uniprot_acc = dplyr::all_of(accessionColumn),
              taxon_raw = dplyr::all_of(taxonColumn)
            ) |>
            dplyr::mutate(
              taxon_id = dplyr::case_when(
                grepl("ID=", taxon_raw) ~ stringr::str_extract(taxon_raw, "(?<=ID=)\\d+"),
                grepl("^\\d+$", as.character(taxon_raw)) ~ as.character(taxon_raw),
                TRUE ~ as.character(taxon_raw)
              )
            ) |>
            dplyr::select(uniprot_acc, taxon_id)
        }, error = function(e) {
          warningMessage <<- e$message
          catFn(sprintf("   ENRICHMENT Step: Error creating organism_mapping: %s\n", e$message))
          NULL
        })

        if (!is.null(organismMapping)) {
          mappingSource <- "taxon_column"
          catFn(sprintf(
            "   ENRICHMENT Step: Created organism_mapping from taxon column (%d entries)\n",
            nrow(organismMapping)
          ))
        }
      } else if (!is.null(accessionColumn)) {
        catFn("   ENRICHMENT Step: No organism column found - creating single-species mapping\n")
        organismMapping <- uniprotDatCln |>
          dplyr::select(uniprot_acc = dplyr::all_of(accessionColumn)) |>
          dplyr::mutate(taxon_id = targetTaxon)
        mappingSource <- "single_species_fallback"
        catFn(sprintf(
          "   ENRICHMENT Step: Created single-species organism_mapping (%d entries, all assigned to taxon %s)\n",
          nrow(organismMapping),
          targetTaxon
        ))
      }
    }
  }

  if (length(availableTaxonIds) == 0 &&
      !is.null(organismMapping) &&
      "taxon_id" %in% names(organismMapping)) {
    availableTaxonIds <- unique(organismMapping$taxon_id)
  }

  list(
    organismMapping = organismMapping,
    source = mappingSource,
    accessionColumn = accessionColumn,
    taxonColumn = taxonColumn,
    availableTaxonIds = availableTaxonIds[!is.na(availableTaxonIds)],
    warning = warningMessage
  )
}

applyProtEnrichOrganismFilter <- function(daResultsForEnrichment,
                                          organismMapping,
                                          targetTaxon,
                                          currentS4Object = NULL,
                                          normalizeUniprotFn = normalizeUniprotAccession,
                                          cleanAccFn = clean_acc,
                                          catFn = cat) {
  filterApplied <- FALSE
  proteinsBefore <- 0
  proteinsAfter <- 0
  proteinsRemoved <- 0
  proteinIdCol <- tryCatch({
    if (!is.null(currentS4Object)) {
      currentS4Object@protein_id_column
    } else {
      "uniprot_acc"
    }
  }, error = function(e) "uniprot_acc")

  targetProteins <- organismMapping |>
    dplyr::filter(taxon_id == targetTaxon) |>
    dplyr::pull(uniprot_acc) |>
    unique()
  targetProteinsClean <- unique(normalizeUniprotFn(targetProteins, remove_isoform = TRUE))

  catFn(sprintf(
    "   ENRICHMENT Step: Found %d proteins for target taxon %s\n",
    length(targetProteins),
    targetTaxon
  ))

  if (!is.null(daResultsForEnrichment@da_data) &&
      length(daResultsForEnrichment@da_data) > 0) {
    filteredDeData <- lapply(names(daResultsForEnrichment@da_data), function(contrastName) {
      contrastData <- daResultsForEnrichment@da_data[[contrastName]]

      if (!is.null(contrastData) && proteinIdCol %in% names(contrastData)) {
        originalCount <- nrow(contrastData)

        keepRows <- vapply(contrastData[[proteinIdCol]], function(proteinIds) {
          ids <- unlist(strsplit(as.character(proteinIds), ";"))
          idsClean <- cleanAccFn(trimws(ids))
          any(idsClean %in% targetProteinsClean) || any(ids %in% targetProteins)
        }, logical(1))

        filteredData <- contrastData[keepRows, , drop = FALSE]
        filteredCount <- nrow(filteredData)

        catFn(sprintf(
          "   ENRICHMENT Step: Contrast '%s': %d -> %d proteins (removed %d non-target organism)\n",
          contrastName,
          originalCount,
          filteredCount,
          originalCount - filteredCount
        ))

        proteinsBefore <<- proteinsBefore + originalCount
        proteinsAfter <<- proteinsAfter + filteredCount
        proteinsRemoved <<- proteinsRemoved + (originalCount - filteredCount)

        return(filteredData)
      }

      contrastData
    })
    names(filteredDeData) <- names(daResultsForEnrichment@da_data)

    daResultsForEnrichment@da_data <- filteredDeData
    filterApplied <- TRUE

    catFn(sprintf(
      "*** ENRICHMENT Step: Organism filtering complete - kept %d/%d proteins (%.1f%%) ***\n",
      proteinsAfter,
      proteinsBefore,
      (proteinsAfter / max(proteinsBefore, 1)) * 100
    ))
  }

  filterStats <- list(
    proteins_before = proteinsBefore,
    proteins_after = proteinsAfter,
    proteins_removed = proteinsRemoved
  )

  list(
    daResultsForEnrichment = daResultsForEnrichment,
    filterApplied = filterApplied,
    filterStats = filterStats,
    proteinIdColumn = proteinIdCol,
    targetProteins = targetProteins,
    targetProteinsClean = targetProteinsClean
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

resolveProtEnrichCurrentS4Object <- function(workflowData, daResultsList) {
  currentS4 <- NULL
  sourceLabel <- NULL

  if (!is.null(workflowData$state_manager)) {
    currentState <- workflowData$state_manager$current_state
    currentS4 <- workflowData$state_manager$getState(currentState)
    if (!is.null(currentS4)) {
      sourceLabel <- "state_manager"
    }
  }

  if (is.null(currentS4) && !is.null(daResultsList) && length(daResultsList) > 0) {
    firstResult <- daResultsList[[1]]
    firstObject <- tryCatch(firstResult$theObject, error = function(e) NULL)
    if (!is.null(firstObject)) {
      currentS4 <- firstObject
      sourceLabel <- "da_results_first_result"
    }
  }

  if (is.null(currentS4) && !is.null(daResultsList)) {
    combinedObject <- tryCatch(daResultsList$theObject, error = function(e) NULL)
    if (!is.null(combinedObject)) {
      currentS4 <- combinedObject
      sourceLabel <- "da_results_combined"
    }
  }

  list(
    currentS4 = currentS4,
    source = sourceLabel
  )
}

buildProtEnrichContrastChoices <- function(daResultsList, contrastsTbl = NULL) {
  contrastNames <- names(daResultsList)
  if (is.null(contrastNames)) {
    contrastNames <- character()
  }

  if (!is.null(contrastsTbl) && "friendly_names" %in% names(contrastsTbl)) {
    friendlyNames <- contrastsTbl$friendly_names
    return(list(
      contrastsAvailable = friendlyNames,
      contrastChoices = setNames(friendlyNames, friendlyNames),
      source = "friendly_names"
    ))
  }

  list(
    contrastsAvailable = contrastNames,
    contrastChoices = setNames(contrastNames, contrastNames),
    source = "raw_names"
  )
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

createProtEnrichRawContrastNameReactive <- function(input,
                                                    globalEnv = .GlobalEnv,
                                                    reactiveFn = shiny::reactive,
                                                    reqFn = shiny::req,
                                                    existsFn = exists,
                                                    getFn = get,
                                                    resolveRawContrastNameFn = resolveProtEnrichRawContrastName) {
  reactiveFn({
    reqFn(input$selected_contrast)

    contrastsTbl <- if (existsFn("contrasts_tbl", envir = globalEnv)) {
      getFn("contrasts_tbl", envir = globalEnv)
    } else {
      NULL
    }

    resolveRawContrastNameFn(input$selected_contrast, contrastsTbl)$rawContrastName
  })
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
    stringdbResults = if (!is.null(contrastResults)) contrastResults$stringdb_results else NULL
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
                                             globalEnv = .GlobalEnv) {
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
    stateObject <- workflowData$state_manager$getState(currentState)
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

  if (is.null(designMatrix) && exists("design_matrix", envir = globalEnv)) {
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

resolveProtEnrichUniprotAnnotations <- function(workflowData,
                                                experimentPaths,
                                                currentS4Object = NULL,
                                                organismTaxid = NULL,
                                                globalEnv = .GlobalEnv,
                                                fileExistsFn = file.exists,
                                                readRdsFn = readRDS,
                                                filePathFn = file.path,
                                                dirExistsFn = dir.exists,
                                                dirCreateFn = dir.create,
                                                getUniprotAnnotationsFn = getUniprotAnnotations,
                                                catFn = cat) {
  uniprotDatCln <- NULL
  annotationSource <- NULL
  sourcePath <- NULL
  cacheDir <- NULL
  cacheDirCreated <- FALSE
  loadError <- NULL
  creationError <- NULL

  if (exists("uniprot_dat_cln", envir = globalEnv)) {
    uniprotDatCln <- get("uniprot_dat_cln", envir = globalEnv)
    annotationSource <- "global_environment"
    catFn("   ENRICHMENT Step: Found uniprot_dat_cln in global environment\n")
  } else if (!is.null(workflowData$uniprot_dat_cln)) {
    uniprotDatCln <- workflowData$uniprot_dat_cln
    assign("uniprot_dat_cln", uniprotDatCln, envir = globalEnv)
    annotationSource <- "workflow_data"
    catFn("   ENRICHMENT Step: Found uniprot_dat_cln in workflow_data\n")
  } else {
    catFn("   ENRICHMENT Step: No uniprot_dat_cln found - checking source directory\n")

    sourcePath <- filePathFn(experimentPaths$source_dir, "uniprot_dat_cln.RDS")
    if (fileExistsFn(sourcePath)) {
      catFn(sprintf("   ENRICHMENT Step: Found uniprot_dat_cln.RDS at %s\n", sourcePath))
      tryCatch({
        uniprotDatCln <- readRdsFn(sourcePath)
        workflowData$uniprot_dat_cln <- uniprotDatCln
        assign("uniprot_dat_cln", uniprotDatCln, envir = globalEnv)
        annotationSource <- "source_directory"
        catFn(sprintf(
          "   ENRICHMENT Step: Successfully loaded %d UniProt annotations from source directory\n",
          nrow(uniprotDatCln)
        ))
      }, error = function(e) {
        loadError <<- e$message
        uniprotDatCln <<- NULL
        catFn(sprintf(
          "   ENRICHMENT Step: Error loading UniProt from source directory: %s\n",
          e$message
        ))
      })
    } else {
      catFn(sprintf(
        "   ENRICHMENT Step: No uniprot_dat_cln.RDS found at %s\n",
        sourcePath
      ))
    }
  }

  if (is.null(uniprotDatCln) && is_test_mode()) {
    uniprotDatCln <- buildProtEnrichDeterministicUniprotAnnotations(
      currentS4Object = currentS4Object,
      organismTaxid = organismTaxid
    )
    workflowData$uniprot_dat_cln <- uniprotDatCln
    assign("uniprot_dat_cln", uniprotDatCln, envir = globalEnv)
    annotationSource <- "test_mode_deterministic"
    catFn(sprintf(
      "   ENRICHMENT Step: Created %d deterministic test-mode UniProt annotations\n",
      nrow(uniprotDatCln)
    ))
  }

  if (is.null(uniprotDatCln)) {
    catFn("   ENRICHMENT Step: Attempting to create UniProt annotations on-the-fly\n")

    tryCatch({
      cacheDir <- filePathFn(experimentPaths$results_dir, "cache", "uniprot_annotations")
      if (!dirExistsFn(cacheDir)) {
        dirCreateFn(cacheDir, recursive = TRUE)
        cacheDirCreated <- TRUE
      }

      if (!is.null(currentS4Object) && !is.null(currentS4Object@protein_quant_table)) {
        uniprotDatCln <- getUniprotAnnotationsFn(
          input_tbl = currentS4Object@protein_quant_table,
          cache_dir = cacheDir,
          taxon_id = as.numeric(organismTaxid)
        )

        workflowData$uniprot_dat_cln <- uniprotDatCln
        assign("uniprot_dat_cln", uniprotDatCln, envir = globalEnv)
        annotationSource <- "generated"
        catFn("   ENRICHMENT Step: Successfully created uniprot_dat_cln on-the-fly\n")
      } else {
        catFn("   ENRICHMENT Step: No protein table available for annotation creation\n")
      }
    }, error = function(e) {
      creationError <<- e$message
      catFn(sprintf("   ENRICHMENT Step: Error creating UniProt annotations: %s\n", e$message))
    })
  }

  list(
    uniprotDatCln = uniprotDatCln,
    source = annotationSource,
    sourcePath = sourcePath,
    cacheDir = cacheDir,
    cacheDirCreated = cacheDirCreated,
    loadError = loadError,
    creationError = creationError
  )
}

resolveProtEnrichAnnotationMatching <- function(uniprotDatCln,
                                                daResultsForEnrichment,
                                                currentS4Object,
                                                matchAnnotationsFn = matchAnnotations,
                                                catFn = cat) {
  if (is.null(uniprotDatCln) ||
      is.null(daResultsForEnrichment) ||
      is.null(currentS4Object)) {
    return(list(
      attempted = FALSE,
      proteinIdColumn = NULL,
      annotationMatchResults = NULL,
      matchRate = NULL,
      warning = NULL
    ))
  }

  catFn("   ENRICHMENT Step: Attempting UniProt annotation matching\n")

  proteinIdCol <- tryCatch(
    currentS4Object@protein_id_column,
    error = function(e) "uniprot_acc"
  )
  annotationMatchResults <- NULL
  matchRate <- NULL
  warningMessage <- NULL

  tryCatch({
    annotationMatchResults <- matchAnnotationsFn(
      da_results_s4 = daResultsForEnrichment,
      uniprot_annotations = uniprotDatCln,
      protein_id_column = proteinIdCol,
      uniprot_id_column = "Entry",
      gene_names_column = "gene_names"
    )

    matchRate <- annotationMatchResults$match_statistics$match_rate
    catFn(sprintf(
      "   ENRICHMENT Step: Annotation matching completed - %d%% match rate\n",
      matchRate
    ))
  }, error = function(e) {
    warningMessage <<- e$message
    annotationMatchResults <<- NULL
    catFn(sprintf("   ENRICHMENT Step: Warning in annotation matching: %s\n", e$message))
    catFn("   ENRICHMENT Step: Continuing with enrichment analysis...\n")
  })

  list(
    attempted = TRUE,
    proteinIdColumn = proteinIdCol,
    annotationMatchResults = annotationMatchResults,
    matchRate = matchRate,
    warning = warningMessage
  )
}
