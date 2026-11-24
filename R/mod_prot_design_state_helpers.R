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
    shiny::removeModal()
    # KIT-175: a complete retrieval removes an earlier annotation notice; an empty or incomplete one shows the notice
    # of its outcome until it is closed, under a fixed id so that a new attempt replaces it.
    outcome <- attr(uniprotDatCln, "annotation_outcome")
    if (nrow(uniprotDatCln) > 0 && is.null(outcome)) {
      log_info(sprintf("UniProt annotations retrieved successfully. Found %d annotations", nrow(uniprotDatCln)))
      shiny::removeNotification("prot_uniprot_annotation")
      shiny::showNotification("UniProt annotations retrieved successfully.", type = "message")
    } else {
      notice <- uniprotAnnotationNotice(outcome)
      log_warn(gsub("[{}]", "", notice))
      shiny::showNotification(notice, type = "warning", duration = NULL, id = "prot_uniprot_annotation")
    }
  }, error = function(e) {
    log_warn(paste("Error getting UniProt annotations:", e$message))
    workflowData$uniprot_dat_cln <- NULL
    shiny::removeModal()
    outcome <- list(reason = e$reason %||% "not_attempted", stage = e$stage,
      detail = e$detail %||% redactUniprotDetail(conditionMessage(e)))
    shiny::showNotification(uniprotAnnotationNotice(outcome), type = "warning", duration = NULL,
      id = "prot_uniprot_annotation")
  })

  if (!is.null(qcTrigger)) {
    qcTrigger(TRUE)
  }


  updatedStatus <- workflowData$tab_status
  updatedStatus$design_matrix <- "complete"
  workflowData$tab_status <- updatedStatus

  if (!is.null(successNotificationId)) {
    shiny::removeNotification(successNotificationId)
  }
  shiny::showNotification(successMessage, type = "message")
}
