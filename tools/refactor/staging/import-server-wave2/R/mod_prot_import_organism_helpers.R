prepareProtImportOrganismSelectionChoices <- function(organismDistribution) {
  valid_organisms <- organismDistribution |>
    dplyr::filter(!is.na(taxon_id) & organism_name != "[Unmatched/Unknown]")

  if (nrow(valid_organisms) == 0) {
    return(NULL)
  }

  choices <- stats::setNames(
    as.character(valid_organisms$taxon_id)
    , paste0(
        valid_organisms$organism_name
        , " (Taxon: ", valid_organisms$taxon_id, ") - "
        , valid_organisms$protein_count, " proteins ("
        , valid_organisms$percentage, "%)"
      )
  )

  list(
    validOrganisms = valid_organisms
    , choices = choices
    , selectedTaxon = as.character(valid_organisms$taxon_id[1])
  )
}

buildProtImportOrganismSelectionModal <- function(
    ns
    , organismDistribution
    , prepareChoices = prepareProtImportOrganismSelectionChoices
) {
  selectionChoices <- prepareChoices(organismDistribution)

  if (is.null(selectionChoices)) {
    return(NULL)
  }

  shiny::modalDialog(
    title = "Select Primary Organism"
    , size = "l"
    , shiny::tags$div(
        shiny::tags$p(
          shiny::tags$strong("Multiple organisms detected in your FASTA database.")
          , " The distribution of matched proteins is shown below."
          , " Please select the primary organism for this analysis."
        )
        , shiny::tags$hr()
        , shiny::tags$h5("Organism Distribution")
        , DT::DTOutput(ns("organism_dist_table"))
        , shiny::tags$hr()
        , shiny::radioButtons(
            ns("selected_organism")
            , "Select primary organism:"
            , choices = selectionChoices$choices
            , selected = selectionChoices$selectedTaxon
          )
        , shiny::checkboxInput(
            ns("filter_to_organism")
            , "Filter data to keep only proteins from selected organism"
            , value = FALSE
          )
        , shiny::helpText(
            "If checked, proteins from other organisms will be removed from the analysis."
          )
      )
    , footer = shiny::tagList(
        shiny::modalButton("Cancel")
        , shiny::actionButton(ns("confirm_organism"), "Confirm Selection", class = "btn-primary")
      )
  )
}

analyzeProtImportMixedSpeciesData <- function(
    workflowData
    , localData
    , dataImportResult
    , fastaPath
    , session
    , updateStatus = updateProtImportProcessingStatus
    , messageFn = message
    , logInfo = logger::log_info
    , logWarn = logger::log_warn
    , showNotification = shiny::showNotification
    , extractOrganisms = extractOrganismsFromFasta
    , analyzeDistribution = analyzeOrganismDistribution
    , buildSelectionModal = buildProtImportOrganismSelectionModal
    , showModal = shiny::showModal
) {
  tryCatch({
    updateStatus("Analyzing organism distribution...")
    messageFn("[mod_prot_import] Mixed species FASTA - analyzing organism distribution")
    logInfo("Mixed species FASTA detected - analyzing organism distribution...")

    organism_mapping <- extractOrganisms(fastaPath)
    localData$organism_mapping <- organism_mapping

    protein_col <- dataImportResult$column_mapping$protein_col
    if (is.null(protein_col) || !(protein_col %in% names(dataImportResult$data))) {
      return(FALSE)
    }

    protein_ids <- unique(dataImportResult$data[[protein_col]])
    organism_dist <- analyzeDistribution(protein_ids, organism_mapping)
    localData$organism_distribution <- organism_dist
    workflowData$organism_distribution <- organism_dist

    logInfo(sprintf("Found %d organisms in data", nrow(organism_dist)))

    organism_modal <- buildSelectionModal(
      ns = session$ns
      , organismDistribution = organism_dist
    )

    if (is.null(organism_modal)) {
      logWarn("No valid organisms with taxon IDs found in FASTA")
      showNotification(
        "Could not extract organism information from FASTA headers. Using default organism settings."
        , type = "warning"
      )
      return(FALSE)
    }

    localData$waiting_for_organism_selection <- TRUE
    showModal(organism_modal)

    TRUE
  }, error = function(e) {
    logWarn(paste("Error analyzing organism distribution:", e$message))
    showNotification(
      paste("Could not analyze organism distribution:", e$message)
      , type = "warning"
    )
    FALSE
  })
}

buildProtImportOrganismDistributionTable <- function(organismDistribution) {
  dist_data <- organismDistribution |>
    dplyr::select(
      `Organism` = organism_name
      , `Taxon ID` = taxon_id
      , `Proteins` = protein_count
      , `%` = percentage
    )

  DT::datatable(
    dist_data
    , options = list(
        pageLength = 10
        , searching = FALSE
        , lengthChange = FALSE
        , ordering = FALSE
        , info = FALSE
      )
    , rownames = FALSE
    , class = "compact stripe"
  ) |>
    DT::formatStyle(
      columns = c("Organism", "Taxon ID", "Proteins", "%")
      , fontSize = "12px"
    )
}

filterProtImportDataToOrganism <- function(
    workflowData
    , organismMapping
    , selectedTaxon
    , organismName
    , normalizeUniprot = normalizeUniprotAccession
    , cleanAcc = clean_acc
    , logInfo = logger::log_info
    , showNotification = shiny::showNotification
) {
  logInfo("Filtering data to selected organism only...")

  organism_proteins <- organismMapping |>
    dplyr::filter(taxon_id == selectedTaxon) |>
    dplyr::pull(uniprot_acc)

  organism_proteins_clean <- unique(normalizeUniprot(organism_proteins, remove_isoform = TRUE))
  protein_col <- workflowData$column_mapping$protein_col

  if (is.null(protein_col) || is.null(workflowData$data_tbl)) {
    return(invisible(NULL))
  }

  original_count <- nrow(workflowData$data_tbl)

  filterByOrganism <- function(proteinIds) {
    ids <- unlist(strsplit(as.character(proteinIds), ";"))
    ids_clean <- cleanAcc(trimws(ids))
    any(ids_clean %in% organism_proteins_clean)
  }

  keep_rows <- sapply(workflowData$data_tbl[[protein_col]], filterByOrganism)
  workflowData$data_tbl <- workflowData$data_tbl[keep_rows, ]

  if (!is.null(workflowData$data_cln)) {
    keep_rows_cln <- sapply(workflowData$data_cln[[protein_col]], filterByOrganism)
    workflowData$data_cln <- workflowData$data_cln[keep_rows_cln, ]
  }

  filtered_count <- nrow(workflowData$data_tbl)
  removed_count <- original_count - filtered_count

  logInfo(sprintf(
    "Filtered data: kept %d rows, removed %d rows (%.1f%%)"
    , filtered_count
    , removed_count
    , (removed_count / original_count) * 100
  ))

  showNotification(
    sprintf(
      "Filtered to %s: kept %d rows, removed %d rows"
      , organismName
      , filtered_count
      , removed_count
    )
    , type = "message"
    , duration = 5
  )

  invisible(list(
    filtered_count = filtered_count
    , removed_count = removed_count
  ))
}

confirmProtImportOrganismSelection <- function(
    selectedTaxon
    , filterToOrganism
    , workflowData
    , localData
    , session
    , updateNumericInput = shiny::updateNumericInput
    , updateTextInput = shiny::updateTextInput
    , showNotification = shiny::showNotification
    , removeModal = shiny::removeModal
    , now = Sys.time
    , logInfo = logger::log_info
    , normalizeUniprot = normalizeUniprotAccession
    , cleanAcc = clean_acc
) {
  selected_org <- localData$organism_distribution |>
    dplyr::filter(taxon_id == selectedTaxon)

  if (nrow(selected_org) == 0) {
    showNotification("Could not find selected organism", type = "error")
    return(FALSE)
  }

  new_organism_name <- selected_org$organism_name[1]
  new_taxon_id <- selected_org$taxon_id[1]

  logInfo(sprintf(
    "User selected organism: %s (Taxon: %d)"
    , new_organism_name
    , new_taxon_id
  ))

  workflowData$taxon_id <- new_taxon_id
  workflowData$organism_name <- new_organism_name

  updateNumericInput(session, "taxon_id", value = new_taxon_id)
  updateTextInput(session, "organism_name", value = new_organism_name)

  if (!is.null(workflowData$processing_log$setup_import)) {
    workflowData$processing_log$setup_import$taxon_id <- new_taxon_id
    workflowData$processing_log$setup_import$organism <- new_organism_name
    workflowData$processing_log$setup_import$mixed_species_selection <- list(
      selected_organism = new_organism_name
      , selected_taxon_id = new_taxon_id
      , filter_applied = filterToOrganism
      , organism_distribution = localData$organism_distribution
    )
  }

  workflowData$mixed_species_analysis <- list(
    enabled = TRUE
    , selected_taxon_id = new_taxon_id
    , selected_organism = new_organism_name
    , organism_distribution = localData$organism_distribution
    , organism_mapping = localData$organism_mapping
    , filter_applied_at_import = filterToOrganism
    , timestamp = now()
  )

  logInfo(sprintf(
    "Stored mixed species analysis metadata: taxon=%d, organism=%s, filter_at_import=%s"
    , new_taxon_id
    , new_organism_name
    , filterToOrganism
  ))

  if (isTRUE(filterToOrganism)) {
    filterProtImportDataToOrganism(
      workflowData = workflowData
      , organismMapping = localData$organism_mapping
      , selectedTaxon = selectedTaxon
      , organismName = new_organism_name
      , normalizeUniprot = normalizeUniprot
      , cleanAcc = cleanAcc
      , logInfo = logInfo
      , showNotification = showNotification
    )
  }

  localData$waiting_for_organism_selection <- FALSE
  removeModal()

  showNotification(
    sprintf("Primary organism set to: %s (Taxon: %d)", new_organism_name, new_taxon_id)
    , type = "message"
  )

  showNotification("Data import successful!", type = "message")

  TRUE
}

