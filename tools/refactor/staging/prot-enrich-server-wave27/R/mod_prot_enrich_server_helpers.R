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

