# KIT-107 W2 (D3): the proteomics (DIA) half of the enrichment core, moved unchanged with their names from the
# enrichment module's helpers at the start commit (mod_prot_enrich_execution_helpers.R:1-194,
# mod_prot_enrich_input_helpers.R:1-274, :599-774). They are proteomics' values for the `enrichment` group of the omic
# configuration (utils_omic_gui_config.R): the deterministic test double (`kernel$test_double`), the organism filter
# (`organism_filter$mapping`, `$apply`) and the UniProt annotation resolver and matcher (`annotation$resolve`, `$match`).
# The omic-neutral entries that call them are in utils_workflow_enrichment.R.

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

resolveProtEnrichUniprotAnnotations <- function(workflowData,
                                                experimentPaths,
                                                currentS4Object = NULL,
                                                organismTaxid = NULL,
                                                globalEnv = NULL,
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

  shareAnnotations <- function(value) {
    if (!is.null(globalEnv)) assign("uniprot_dat_cln", value, envir = globalEnv)
    invisible(value)
  }
  if (!is.null(workflowData$uniprot_dat_cln)) {
    uniprotDatCln <- workflowData$uniprot_dat_cln
    shareAnnotations(uniprotDatCln)
    annotationSource <- "workflow_data"
    catFn("   ENRICHMENT Step: Found uniprot_dat_cln in workflow_data\n")
  } else if (!is.null(globalEnv) && exists("uniprot_dat_cln", envir = globalEnv)) {
    # Script/test seam only: the application passes no environment (KIT-013).
    uniprotDatCln <- get("uniprot_dat_cln", envir = globalEnv)
    annotationSource <- "global_environment"
    catFn("   ENRICHMENT Step: Found uniprot_dat_cln in workflow_data\n")
  } else {
    catFn("   ENRICHMENT Step: No uniprot_dat_cln found - checking source directory\n")

    sourcePath <- filePathFn(experimentPaths$source_dir, "uniprot_dat_cln.RDS")
    if (fileExistsFn(sourcePath)) {
      catFn(sprintf("   ENRICHMENT Step: Found uniprot_dat_cln.RDS at %s\n", sourcePath))
      tryCatch({
        uniprotDatCln <- readRdsFn(sourcePath)
        workflowData$uniprot_dat_cln <- uniprotDatCln
        shareAnnotations(uniprotDatCln)
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
    shareAnnotations(uniprotDatCln)
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
        shareAnnotations(uniprotDatCln)
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
