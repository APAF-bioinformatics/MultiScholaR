persistProtDesignBuilderArtifacts <- function(results, workflowData, sourceDir) {
  designMatrixPath <- file.path(sourceDir, "design_matrix.tab")
  logger::log_info(paste("Writing design matrix to:", designMatrixPath))
  utils::write.table(results$design_matrix, file = designMatrixPath, sep = "\t", row.names = FALSE, quote = FALSE)

  dataClnPath <- file.path(sourceDir, "data_cln.tab")
  logger::log_info(paste("Writing cleaned data to:", dataClnPath))
  writeDIADesignTable(workflowData$state_manager, results$data_reference, dataClnPath)

  if (!is.null(results$contrasts_tbl) && nrow(results$contrasts_tbl) > 0) {
    contrastPath <- file.path(sourceDir, "contrast_strings.tab")
    logger::log_info(paste("Writing contrasts to:", contrastPath))
    writeLines(results$contrasts_tbl$contrasts, contrastPath)
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
  workflowData <- setDIAOwnedPayload(workflowData, "data_cln",
    results$data_reference %||% results$data_cln)
  workflowData$contrasts_tbl <- results$contrasts_tbl
  workflowData$config_list <- results$config_list
}

runProtDesignBuilderSaveFlow <- function(
    results,
    workflowData,
    experimentPaths,
    session,
    qcTrigger = NULL,
    cut = FALSE
) {
  sourceDir <- experimentPaths$source_dir
  if (is.null(sourceDir) || !dir.exists(sourceDir)) {
    msg <- "Could not find the source directory to save files. Check experiment paths."
    logger::log_error(msg)
    shiny::showNotification(msg, type = "error", duration = 15)
    return(FALSE)
  }

  # KIT-118 (P-10): the contrasts are checked against the design's model matrix before anything is
  # written, as lipid and metab do at save (`validateOmicsDesignSave`): a contrast naming a term the
  # formula does not produce is refused here, not later at differential abundance. The sample
  # alignment half of that preflight does not apply: the builder lists the data's own runs.
  preflight <- validateOmicsDesignContrasts(
    designMatrix = results$design_matrix, contrastsTbl = results$contrasts_tbl,
    rules = list(label = "Proteomics"),
    formulaString = resolveOmicsDesignFormula(results$config_list), allowEmpty = TRUE
  )
  if (!isTRUE(preflight$valid)) {
    msg <- paste("Invalid proteomics design:", paste(preflight$errors, collapse = "; "))
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
    workflowData$prot_design_log2_counts <- NULL
    stateName <- buildProtDesignStateCheckpoint(
      workflowData = workflowData,
      workflowType = workflowType,
      actionLabel = "Save Design"
    )
    # KIT-247 (VAL-239, OD-1 (D)): a changed design confirmed over later steps sets them aside.
    rememberOmicsDesignCommit(workflowData$state_manager, workflowData$design_matrix, workflowData$data_cln,
      workflowData$config_list)
    didCut <- isTRUE(cut) && cutOmicsDesignHistory(workflowData$state_manager, stateName)
    # KIT-177: the values the design step set to missing, named and counted.
    notifyProtDesignLog2Counts(workflowData)
    completeProtDesignPostCheckpoint(
      workflowData = workflowData,
      experimentPaths = experimentPaths,
      session = session,
      qcTrigger = qcTrigger,
      successMessage = "Design matrix and contrasts saved successfully!",
      debugQcTrigger = TRUE
    )
    if (didCut) completeOmicsDesignCut(workflowData, "proteomics", experimentPaths)

    manager <- workflowData$state_manager
    reference <- getDIAStateReference(manager)
    requestDIAReclamation(manager, reference, paste0("design-", reference$entity_id), "design")
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

  commit <- function(cut) {
    showProcessingModal()
    logInfo("Received results from design matrix builder module. Now saving to workflow and disk.")

    hydrateBuilderResults(
      results = results,
      workflowData = workflowData
    )

    if (isTRUE(cut)) {
      runBuilderSaveFlow(results = results, workflowData = workflowData, experimentPaths = experimentPaths,
        session = session, qcTrigger = qcTrigger, cut = TRUE)
    } else {
      runBuilderSaveFlow(
        results = results,
        workflowData = workflowData,
        experimentPaths = experimentPaths,
        session = session,
        qcTrigger = qcTrigger
      )
    }
  }
  # KIT-247 (VAL-239, OD-1 (D)): compared before the bus is hydrated; a changed design over later steps asks first.
  manager <- workflowData$state_manager
  if (!inherits(manager, "DIAStateAdapter")) return(commit(FALSE))
  stateName <- protDesignStateName(manager)
  verdict <- compareOmicsDesignResult(manager,
    list(design_matrix = results$design_matrix, data = results$data_reference %||% results$data_cln), stateName,
    omicsDesignCompared(manager, workflowData))
  later <- length(omicsDesignLaterStates(manager, stateName)) > 0L
  if (!later || !verdict %in% c("unchanged", "changed")) return(commit(FALSE))
  # KIT-277 (stage 1): a formula-only change asks and is carried through on Save.
  formula <- if (identical(verdict, "unchanged")) omicsDesignFormulaChange(manager, workflowData, results$config_list)
  runOmicsDesignResaveGate(session, verdict, later, commit = commit, unchanged = function() {
    applyOmicsDesignUnchanged(workflowData, results, "differential_expression", quiet = !is.null(formula))
    sourceDir <- experimentPaths$source_dir
    if (!is.null(sourceDir) && dir.exists(sourceDir)) {
      tryCatch(persistProtDesignBuilderArtifacts(results = results, workflowData = workflowData, sourceDir = sourceDir),
        error = function(e) logger::log_warn(paste("Design files not written:", e$message)))
    }
    if (!is.null(formula)) {
      tryCatch(shiny::showNotification(carryOmicsDesignFormula(workflowData, "proteomics", formula[["to"]],
        "differential_expression"), type = "message"), error = function(e) {
          msg <- paste("Error saving design matrix results:", conditionMessage(e))
          logger::log_error(msg)
          shiny::showNotification(msg, type = "error", duration = 15)
        })
    }
  }, limited = if (!is.null(formula)) list(formula = formula) else character())
}

#' KIT-247. The proteomics design state of the owner's route (`buildProtDesignStateCheckpoint`).
#' @noRd
protDesignStateName <- function(manager) {
  if (isTRUE(shiny::isolate(manager$workflow_type) %in% c("DIA", "DIA_LIMPA"))) "raw_data_s4" else "protein_s4_initial"
}

# KIT-103 W6: the saved-results previews are the shared design tab's (the status outputs are the
# frame's); kept by name as `registerProtDesignServerShells`' default preview registration.
registerProtDesignPreviewOutputs <- function(output, workflowData) {
  registerOmicsDesignPreviewOutputs(output, workflowData, omicGuiConfig("proteomics"))
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

    data_reference <- reactiveFn(workflowData$data_tbl)
    data_provider <- function() resolveDIAOwnedValue(data_reference())
    attr(data_provider, "dia_reference") <- data_reference
    return(builderServerFn(
      moduleId,
      data_tbl = data_provider,
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
  observeEventFn(getDIAReactiveReference(builderResultsRv), {
    withDIAStateScope(workflowData$state_manager, function() {
    results <- builderResultsRv()
    reqFn(results)

    runBuilderObserverShell(
      results = results,
      workflowData = workflowData,
      experimentPaths = experimentPaths,
      session = session,
      qcTrigger = qcTrigger
    )
    })
  }, ignoreNULL = TRUE)
}
