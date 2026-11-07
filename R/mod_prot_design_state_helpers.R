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

    # KIT-103 W3: the computation is the core's annotation step.
    uniprotDatCln <- retrieveProtDesignUniprotAnnotations(
      manager = workflowData$state_manager,
      dataCln = workflowData$data_cln,
      proteinColumn = workflowData$column_mapping$protein_col,
      taxonId = workflowData$taxon_id,
      resultsDir = experimentPaths$results_dir,
      sourceDir = experimentPaths$source_dir,
      progressCallback = progressUpdater
    )

    workflowData$uniprot_dat_cln <- uniprotDatCln
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
