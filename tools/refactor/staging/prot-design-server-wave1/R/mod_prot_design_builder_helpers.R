persistProtDesignBuilderArtifacts <- function(results, workflowData, sourceDir) {
  designMatrixPath <- file.path(sourceDir, "design_matrix.tab")
  logger::log_info(paste("Writing design matrix to:", designMatrixPath))
  utils::write.table(results$design_matrix, file = designMatrixPath, sep = "\t", row.names = FALSE, quote = FALSE)

  dataClnPath <- file.path(sourceDir, "data_cln.tab")
  logger::log_info(paste("Writing cleaned data to:", dataClnPath))
  utils::write.table(results$data_cln, file = dataClnPath, sep = "\t", row.names = FALSE, quote = FALSE)

  if (!is.null(results$contrasts_tbl) && nrow(results$contrasts_tbl) > 0) {
    contrastPath <- file.path(sourceDir, "contrast_strings.tab")
    logger::log_info(paste("Writing contrasts to:", contrastPath))
    writeLines(results$contrasts_tbl$contrasts, contrastPath)

    assign("contrasts_tbl", results$contrasts_tbl, envir = .GlobalEnv)
    logger::log_info("Saved contrasts_tbl to global environment for DE analysis.")
  }

  manifestPath <- file.path(sourceDir, "manifest.json")
  manifestData <- list(
    data_path = "data_cln.tab",
    design_matrix_path = "design_matrix.tab",
    contrast_strings_path = "contrast_strings.tab"
  )
  jsonlite::write_json(manifestData, manifestPath, auto_unbox = TRUE, pretty = TRUE)
  logger::log_info(paste("Saved manifest.json to:", manifestPath))

  if (!is.null(workflowData$config_list)) {
    if (!is.null(workflowData$state_manager$workflow_type)) {
      workflowData$config_list$globalParameters$workflow_type <- workflowData$state_manager$workflow_type
      logger::log_info(paste("Added workflow_type to config.ini:", workflowData$state_manager$workflow_type))
    }

    configPath <- file.path(sourceDir, "config.ini")
    logger::log_info(paste("Writing config.ini to:", configPath))
    tryCatch({
      ini::write.ini(workflowData$config_list, configPath)
      logger::log_info("Saved config.ini for future import/export.")
    }, error = function(e) {
      logger::log_warn(paste("Could not save config.ini:", e$message))
      logger::log_warn("Import will use default config if this design is imported.")
    })
  } else {
    logger::log_warn("config_list is NULL, cannot save config.ini.")
  }
}

hydrateProtDesignBuilderResults <- function(results, workflowData) {
  workflowData$design_matrix <- results$design_matrix
  workflowData$data_cln <- results$data_cln
  workflowData$contrasts_tbl <- results$contrasts_tbl
  workflowData$config_list <- results$config_list

  if (!is.null(results$contrasts_tbl)) {
    assign("contrasts_tbl", results$contrasts_tbl, envir = .GlobalEnv)
    logger::log_info("Updated contrasts_tbl in global environment from design builder.")
  }

  assign("config_list", workflowData$config_list, envir = .GlobalEnv)
  logger::log_info("Updated global config_list for updateConfigParameter compatibility")
}

runProtDesignBuilderSaveFlow <- function(
    results,
    workflowData,
    experimentPaths,
    session,
    qcTrigger = NULL
) {
  sourceDir <- experimentPaths$source_dir
  if (is.null(sourceDir) || !dir.exists(sourceDir)) {
    msg <- "Could not find the source directory to save files. Check experiment paths."
    logger::log_error(msg)
    shiny::showNotification(msg, type = "error", duration = 15)
    return(FALSE)
  }

  tryCatch({
    persistProtDesignBuilderArtifacts(
      results = results,
      workflowData = workflowData,
      sourceDir = sourceDir
    )

    workflowType <- shiny::isolate(workflowData$state_manager$workflow_type)
    buildProtDesignStateCheckpoint(
      workflowData = workflowData,
      workflowType = workflowType,
      actionLabel = "Save Design"
    )
    completeProtDesignPostCheckpoint(
      workflowData = workflowData,
      experimentPaths = experimentPaths,
      session = session,
      qcTrigger = qcTrigger,
      successMessage = "Design matrix and contrasts saved successfully!",
      debugQcTrigger = TRUE
    )

    TRUE
  }, error = function(e) {
    msg <- paste("Error saving design matrix results:", e$message)
    logger::log_error(msg)
    shiny::showNotification(msg, type = "error", duration = 15)
    FALSE
  })
}

runProtDesignBuilderObserverShell <- function(
    results,
    workflowData,
    experimentPaths,
    session,
    qcTrigger = NULL,
    showProcessingModal = NULL,
    hydrateBuilderResults = hydrateProtDesignBuilderResults,
    runBuilderSaveFlow = runProtDesignBuilderSaveFlow,
    logInfo = logger::log_info
) {
  if (is.null(showProcessingModal)) {
    showProcessingModal <- function() {
      shiny::showModal(shiny::modalDialog(
        title = "Processing Design Matrix",
        shiny::div(
          style = "text-align: center; padding: 20px;",
          shiny::icon("spinner", class = "fa-spin fa-3x"),
          shiny::br(),
          shiny::br(),
          shiny::p("Saving design matrix and preparing data..."),
          shiny::p(shiny::tags$small("This may take a moment for large datasets."))
        ),
        footer = NULL,
        easyClose = FALSE
      ))
    }
  }

  showProcessingModal()
  logInfo("Received results from design matrix builder module. Now saving to workflow and disk.")

  hydrateBuilderResults(
    results = results,
    workflowData = workflowData
  )

  runBuilderSaveFlow(
    results = results,
    workflowData = workflowData,
    experimentPaths = experimentPaths,
    session = session,
    qcTrigger = qcTrigger
  )
}

registerProtDesignPreviewOutputs <- function(
    output,
    workflowData,
    reactiveFn = shiny::reactive,
    outputOptionsFn = shiny::outputOptions,
    renderDT = DT::renderDT,
    reqFn = shiny::req
) {
  output$data_available <- reactiveFn({
    !is.null(workflowData$data_tbl) && !is.null(workflowData$config_list)
  })
  outputOptionsFn(output, "data_available", suspendWhenHidden = FALSE)

  output$design_matrix_exists <- reactiveFn({
    !is.null(workflowData$design_matrix)
  })
  outputOptionsFn(output, "design_matrix_exists", suspendWhenHidden = FALSE)

  output$design_matrix_preview <- renderDT({
    reqFn(workflowData$design_matrix)
    workflowData$design_matrix
  }, options = list(pageLength = 5, scrollX = TRUE))

  output$contrasts_preview <- renderDT({
    reqFn(workflowData$contrasts_tbl)
    workflowData$contrasts_tbl
  }, options = list(pageLength = 5, scrollX = TRUE))

  invisible(output)
}

registerProtDesignBuilderModule <- function(
    workflowData,
    moduleId = "builder",
    builderServerExists = exists("mod_prot_design_builder_server"),
    builderServerFn = NULL,
    reactiveFn = shiny::reactive,
    reactiveValFn = shiny::reactiveVal
) {
  if (isTRUE(builderServerExists)) {
    if (is.null(builderServerFn)) {
      builderServerFn <- mod_prot_design_builder_server
    }

    return(builderServerFn(
      moduleId,
      data_tbl = reactiveFn(workflowData$data_tbl),
      config_list = reactiveFn(workflowData$config_list),
      column_mapping = reactiveFn(workflowData$column_mapping)
    ))
  }

  reactiveValFn(NULL)
}

registerProtDesignBuilderResultsObserver <- function(
    builderResultsRv,
    workflowData,
    experimentPaths,
    session,
    qcTrigger = NULL,
    observeEventFn = shiny::observeEvent,
    reqFn = shiny::req,
    runBuilderObserverShell = runProtDesignBuilderObserverShell
) {
  observeEventFn(builderResultsRv(), {
    results <- builderResultsRv()
    reqFn(results)

    runBuilderObserverShell(
      results = results,
      workflowData = workflowData,
      experimentPaths = experimentPaths,
      session = session,
      qcTrigger = qcTrigger
    )
  }, ignoreNULL = TRUE)
}

