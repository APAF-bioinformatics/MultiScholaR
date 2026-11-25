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
        # KIT-265 (VAL-291 D-2): the single-organism limit, stated where the organism is chosen (draft, VAL-41).
        , shiny::helpText(PROT_IMPORT_ORGANISM_LIMIT_TEXT)
        , shiny::checkboxInput(
            ns("filter_to_organism")
            , "Filter data to keep only proteins from selected organism"
            , value = FALSE
          )
        , shiny::helpText(
            "If checked, proteins from other organisms will be removed from the analysis."
          )
      )
    # KIT-250: Cancel is the module's own button (a modalButton runs no server code); its observer clears the import.
    , footer = shiny::tagList(
        shiny::actionButton(ns("cancel_organism"), "Cancel")
        , shiny::actionButton(ns("confirm_organism"), "Confirm Selection", class = "btn-primary")
      )
  )
}

PROT_IMPORT_ORGANISM_LIMIT_TEXT <- paste0(
  "Only the organism you choose here is looked up in UniProt and tested for enrichment; the other organisms' "
  , "proteins are not, and they stay in the normalisation, RUV and the differential analysis (and in its correction "
  , "for multiple testing) unless you tick 'Filter data to keep only proteins from selected organism'."
)

#' KIT-265: with 'Mixed species FASTA' ticked and no organism window, the import never keeps an organism nobody chose
#' without saying so. A refusal is a classed error that the import's own error path (`runProtImportProcessingSafely`)
#' shows as "Error: <message>" and uses to clear what the import wrote. `rlang::abort`, not `cli::cli_abort`: the
#' messages embed an error's text and the form's name, which must not be interpolated.
#' @noRd
abortProtImportOrganismUnavailable <- function(message, parent = NULL) {
  rlang::abort(message, class = "prot_import_organism_unavailable", parent = parent)
}

PROT_IMPORT_UNTICK_ADVICE <- "untick 'Mixed species FASTA', enter the organism yourself and press 'Process Imported Data' again."

#' The persistent notice of route (i): no organism with a taxonomy ID could be matched to the data's proteins, so no
#' window is shown and the form's organism is used. Names that organism, or says the form has none.
#' @noRd
buildProtImportNoOrganismNotice <- function(organismName = NULL, taxonId = NULL) {
  has_name <- length(organismName) == 1L && !is.na(organismName) && nzchar(trimws(organismName))
  has_taxon <- length(taxonId) == 1L && !is.na(taxonId)

  used <- if (has_name && has_taxon) {
    sprintf(
      "This analysis will use the organism entered in the form: %s (taxon %s). "
      , trimws(organismName)
      , format(taxonId, scientific = FALSE)
    )
  } else {
    paste0(
      "The form has no organism (its name or Taxonomy ID is empty), so UniProt annotation and enrichment will have "
      , "no organism to use. "
    )
  }

  paste0(
    "No organism with a taxonomy ID could be matched to the proteins in the search results (the FASTA headers name "
    , "no organism or taxonomy ID for them, or the proteins are not in the FASTA file), so no organism window is "
    , "shown. "
    , used
    , "To use another organism, untick 'Mixed species FASTA', enter it and press 'Process Imported Data' again."
  )
}

#' Whether the FASTA names organisms (an entry with both an organism and a taxon ID) and none of the data's proteins
#' is in it: route (i)(c), the usual sign of a FASTA other than the one the search used. Accessions are split at ";",
#' trimmed, blanks dropped and isoform suffixes removed on both sides, as `analyzeOrganismDistribution` does. A data
#' accession found under an entry without organism fields is not (c): that FASTA does contain the protein.
#' @noRd
protImportFastaLacksDataProteins <- function(organismMapping, proteinIds) {
  needed <- c("uniprot_acc", "organism_name", "taxon_id")
  if (!is.data.frame(organismMapping) || !all(needed %in% names(organismMapping))) {
    return(FALSE)
  }
  if (!any(!is.na(organismMapping$organism_name) & !is.na(organismMapping$taxon_id))) {
    return(FALSE)
  }

  cleanAccessions <- function(x) {
    accessions <- trimws(unlist(strsplit(as.character(x), ";", fixed = TRUE)))
    accessions <- accessions[!is.na(accessions) & accessions != ""]
    unique(normalizeUniprotAccession(accessions, remove_isoform = TRUE))
  }

  data_accessions <- cleanAccessions(proteinIds)
  if (length(data_accessions) == 0L) {
    return(FALSE)
  }

  !any(data_accessions %in% cleanAccessions(organismMapping$uniprot_acc))
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
    , organismName = NULL
    , taxonId = NULL
) {
  tryCatch({
    updateStatus("Analyzing organism distribution...")
    messageFn("[mod_prot_import] Mixed species FASTA - analyzing organism distribution")
    logInfo("Mixed species FASTA detected - analyzing organism distribution...")

    organism_mapping <- extractOrganisms(fastaPath)
    localData$organism_mapping <- organism_mapping

    protein_col <- dataImportResult$column_mapping$protein_col
    # KIT-265: without the protein column the organisms cannot be matched to the data: the import is refused (it
    # kept the form's organism silently). The owner path's metadata carries `protein_ids`, NULL when the column
    # was absent from the table.
    owner_path <- !is.null(dataImportResult$import_summary)
    if (is.null(protein_col) || (!owner_path && !(protein_col %in% names(dataImportResult$data))) ||
        (owner_path && is.null(dataImportResult$protein_ids))) {
      column_text <- if (is.null(protein_col)) "protein column" else sprintf("'%s' column", protein_col)
      abortProtImportOrganismUnavailable(paste0(
        "The import was stopped: the search results have no ", column_text, ", so their proteins cannot be matched "
        , "to the organisms in the FASTA file. Nothing was imported. Export the search results with that column and "
        , "import them again."
      ))
    }

    protein_ids <- if (owner_path) dataImportResult$protein_ids else unique(dataImportResult$data[[protein_col]])
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
      # KIT-265 (VAL-304, route (i)(c)): the FASTA names organisms and none of the data's proteins is in it.
      if (protImportFastaLacksDataProteins(organism_mapping, protein_ids)) {
        abortProtImportOrganismUnavailable(paste0(
          "The import was stopped: none of the proteins in the search results could be found in the FASTA file, so "
          , "their organisms cannot be counted. Nothing was imported. Check that this is the FASTA file the search "
          , "used, or ", PROT_IMPORT_UNTICK_ADVICE
        ))
      }
      # KIT-265: routes (i)(a) and (i)(b): the import continues on the form's organism, and the notice says which,
      # until closed (the one notice that an unconfirmed organism is in use; it shares Cancel's slot).
      showNotification(
        buildProtImportNoOrganismNotice(organismName, taxonId)
        , type = "warning"
        , duration = NULL
        , id = PROT_IMPORT_ORGANISM_CANCELLED_ID
      )
      return(FALSE)
    }

    localData$waiting_for_organism_selection <- TRUE
    showModal(organism_modal)

    TRUE
  }, error = function(e) {
    logWarn(paste("Error analyzing organism distribution:", e$message))
    # KIT-265: a refusal raised above passes through unchanged; any other error is route (ii) and refuses the import.
    if (inherits(e, "prot_import_organism_unavailable")) {
      stop(e)
    }
    abortProtImportOrganismUnavailable(
      paste0(
        "The import was stopped: the organisms in the FASTA file could not be counted (", e$message, "). The "
        , "organism is chosen from that count, so nothing was imported. Check the FASTA file, or "
        , PROT_IMPORT_UNTICK_ADVICE
      )
      , parent = e
    )
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
    , complete = function(notify) completeOmicImport(workflowData, omicGuiConfig("proteomics"),
        omicImportSourcePanel("proteomics")$design_status, notify)
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

  # KIT-250: with 'Mixed species FASTA' the import completes here, on the confirmed organism, after the optional filter.
  complete(notify = FALSE)

  localData$waiting_for_organism_selection <- FALSE
  localData$organism_unconfirmed <- FALSE
  removeModal()

  showNotification(
    sprintf("Primary organism set to: %s (Taxon: %d)", new_organism_name, new_taxon_id)
    , type = "message"
  )

  showNotification("Data import successful!", type = "message")

  TRUE
}

#' KIT-250: Cancel in the organism dialog. The import is not finished: it is cleared as a failed import is
#' (`clearProtImportWorkflowState`), so no design can be built or saved until an organism is confirmed or the option
#' is unticked, and the user is told so in a notification that stays until closed.
#' @return The new statuses, invisibly.
#' @noRd
cancelProtImportOrganismSelection <- function(
    workflowData
    , localData
    , removeModal = shiny::removeModal
    , showNotification = shiny::showNotification
    , logInfo = logger::log_info
) {
  localData$waiting_for_organism_selection <- FALSE
  localData$organism_unconfirmed <- TRUE
  removeModal()
  updated_status <- clearProtImportWorkflowState(workflowData)
  logInfo("Organism selection cancelled: the import was cleared (no organism confirmed)")
  showNotification(PROT_IMPORT_ORGANISM_CANCELLED_TEXT, type = "warning", duration = NULL,
    id = PROT_IMPORT_ORGANISM_CANCELLED_ID)

  invisible(updated_status)
}
