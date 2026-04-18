buildProtDesignStateCheckpoint <- function(
    workflowData,
    workflowType,
    actionLabel,
    validateColumnMapping = FALSE
) {
  workflowData$design_matrix <- workflowData$design_matrix |>
    dplyr::mutate(tech_rep_group = paste(group, replicates, sep = "_"))
  logger::log_info("Created tech_rep_group column for proper technical replicate identification")

  if (workflowType == "DIA") {
    s4Object <- PeptideQuantitativeDataDiann(
      peptide_data = workflowData$data_cln,
      design_matrix = workflowData$design_matrix,
      technical_replicate_id = "tech_rep_group",
      args = workflowData$config_list
    )
    stateName <- "raw_data_s4"
    description <- "Initial Peptide S4 object created after design matrix."
  } else if (workflowType %in% c("TMT", "LFQ")) {
    logger::log_info(sprintf("%s workflow: Reshaping protein data from long to wide format for S4 creation.", workflowType))

    if (validateColumnMapping) {
      if (is.null(workflowData$column_mapping$protein_col) ||
          is.null(workflowData$column_mapping$run_col) ||
          is.null(workflowData$column_mapping$quantity_col)) {
        stop(sprintf("Missing required column mappings for %s data reshape. Cannot proceed.", workflowType))
      }

      requiredCols <- c(
        workflowData$column_mapping$protein_col,
        workflowData$column_mapping$run_col,
        workflowData$column_mapping$quantity_col
      )
      missingCols <- requiredCols[!requiredCols %in% names(workflowData$data_cln)]
      if (length(missingCols) > 0) {
        stop(paste("Missing required columns in data_cln:", paste(missingCols, collapse = ", ")))
      }
    }

    proteinQuantTableWide <- workflowData$data_cln |>
      tidyr::pivot_wider(
        id_cols = !!sym(workflowData$column_mapping$protein_col),
        names_from = !!sym(workflowData$column_mapping$run_col),
        values_from = !!sym(workflowData$column_mapping$quantity_col)
      )

    logger::log_info(sprintf(
      "%s %s: Pivoted data from %d rows (long) to %d rows (wide)",
      workflowType,
      actionLabel,
      nrow(workflowData$data_cln),
      nrow(proteinQuantTableWide)
    ))
    logger::log_info(sprintf(
      "%s %s: Wide format has %d protein rows, %d sample columns",
      workflowType,
      actionLabel,
      nrow(proteinQuantTableWide),
      ncol(proteinQuantTableWide) - 1
    ))

    logger::log_info(sprintf("%s %s: Applying log2 transformation to abundance values...", workflowType, actionLabel))

    proteinIdCol <- workflowData$column_mapping$protein_col
    proteinQuantTableWide <- proteinQuantTableWide |>
      dplyr::mutate(
        dplyr::across(
          -!!sym(proteinIdCol),
          ~ log2(.x + 1)
        )
      )

    logger::log_info(sprintf("%s %s: Log2 transformation completed", workflowType, actionLabel))

    s4Object <- ProteinQuantitativeData(
      protein_quant_table = proteinQuantTableWide,
      protein_id_column = workflowData$column_mapping$protein_col,
      protein_id_table = proteinQuantTableWide |> dplyr::distinct(!!sym(workflowData$column_mapping$protein_col)),
      design_matrix = workflowData$design_matrix,
      sample_id = "Run",
      group_id = "group",
      technical_replicate_id = "tech_rep_group",
      args = workflowData$config_list
    )
    stateName <- "protein_s4_initial"
    description <- sprintf("Initial Protein S4 object created from %s data after design matrix.", workflowType)
  } else {
    stop("Unknown workflow_type in designMatrixAppletServer: ", workflowType)
  }

  logger::log_info(paste("Saving S4 object to R6 state manager as:", stateName))
  workflowData$state_manager$saveState(
    state_name = stateName,
    s4_data_object = s4Object,
    config_object = workflowData$config_list,
    description = description
  )

  .capture_checkpoint(s4Object, "cp04", "design_matrix")

  stateName
}

completeProtDesignPostCheckpoint <- function(
    workflowData,
    experimentPaths,
    session,
    qcTrigger = NULL,
    successMessage,
    successNotificationId = NULL,
    debugQcTrigger = FALSE
) {
  log_info("Design Matrix complete. Triggering UniProt annotation.")

  shiny::showModal(shiny::modalDialog(
    title = "Retrieving UniProt Annotations",
    shiny::p("Please wait while we retrieve protein annotations from UniProt..."),
    shiny::div(id = "uniprot_progress_text", "Initializing..."),
    shiny::tags$div(class = "progress",
      shiny::tags$div(
        id = "uniprot_progress_bar",
        class = "progress-bar progress-bar-striped active",
        role = "progressbar",
        style = "width: 0%",
        "0%"
      )
    ),
    footer = NULL,
    easyClose = FALSE
  ))

  tryCatch({
    proteinColumn <- workflowData$column_mapping$protein_col
    if (is.null(proteinColumn)) {
      stop("Protein column not found in column_mapping. Cannot get annotations.")
    }

    cacheDir <- if (!is.null(experimentPaths) && !is.null(experimentPaths$results_dir)) {
      file.path(experimentPaths$results_dir, "cache")
    } else {
      file.path(tempdir(), "proteomics_cache")
    }

    uniprotCacheDir <- file.path(cacheDir, "uniprot_annotations")
    if (!dir.exists(uniprotCacheDir)) {
      dir.create(uniprotCacheDir, recursive = TRUE)
    }

    progressUpdater <- function(current, total) {
      tryCatch({
        percent <- round((current / total) * 100)
        session$sendCustomMessage("updateUniprotProgress", list(
          percent = percent,
          text = sprintf("Processing chunk %d of %d (%d%%)", current, total, percent)
        ))
      }, error = function(e) {
        message(paste("Progress update failed:", e$message))
      })
    }

    uniprotDatCln <- getUniprotAnnotationsFull(
      data_tbl = workflowData$data_cln,
      protein_id_column = proteinColumn,
      cache_dir = uniprotCacheDir,
      taxon_id = workflowData$taxon_id,
      progress_callback = progressUpdater
    )

    workflowData$uniprot_dat_cln <- uniprotDatCln
    assign("uniprot_dat_cln", uniprotDatCln, envir = .GlobalEnv)

    if (!is.null(experimentPaths) && !is.null(experimentPaths$source_dir)) {
      scriptsUniprotPath <- file.path(experimentPaths$source_dir, "uniprot_dat_cln.RDS")
      saveRDS(uniprotDatCln, scriptsUniprotPath)
      log_info(sprintf("Saved uniprot_dat_cln to scripts directory: %s", scriptsUniprotPath))
    }
    log_info(sprintf("UniProt annotations retrieved successfully. Found %d annotations", nrow(uniprotDatCln)))
    shiny::removeModal()
    shiny::showNotification("UniProt annotations retrieved successfully.", type = "message")
  }, error = function(e) {
    log_warn(paste("Error getting UniProt annotations:", e$message))
    workflowData$uniprot_dat_cln <- NULL
    shiny::removeModal()
    shiny::showNotification(
      paste("Warning: Could not retrieve UniProt annotations:", e$message),
      type = "warning",
      duration = 8
    )
  })

  if (debugQcTrigger) {
    message("=== DEBUG66: designMatrixApplet Save Design - About to set qc_trigger ===")
    message(sprintf("   DEBUG66: qc_trigger is NULL = %s", is.null(qcTrigger)))
  }

  if (!is.null(qcTrigger)) {
    if (debugQcTrigger) {
      message("   DEBUG66: Setting qc_trigger(TRUE)")
    }
    qcTrigger(TRUE)
    if (debugQcTrigger) {
      message(sprintf("   DEBUG66: qc_trigger set. Current value = %s", qcTrigger()))
    }
  } else if (debugQcTrigger) {
    message("   DEBUG66: qc_trigger is NULL, cannot set")
  }

  if (debugQcTrigger) {
    message("=== DEBUG66: designMatrixApplet - qc_trigger setting complete ===")
  }

  updatedStatus <- workflowData$tab_status
  updatedStatus$design_matrix <- "complete"
  workflowData$tab_status <- updatedStatus

  if (!is.null(successNotificationId)) {
    shiny::removeNotification(successNotificationId)
  }
  shiny::showNotification(successMessage, type = "message")
}

