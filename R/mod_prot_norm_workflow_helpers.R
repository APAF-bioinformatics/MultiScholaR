#' Record the chosen between-sample normalisation method on the workflow bus configuration.
#' @param workflowData Proteomics workflow bus (`reactiveValues` or environment) or NULL.
#' @param normMethod Method name chosen in the GUI.
#' @return NULL invisibly; a bus without a configuration is left untouched.
#' @noRd
recordProtNormMethod <- function(workflowData, normMethod) {
  if (is.null(workflowData)) return(invisible(NULL))
  config_list <- shiny::isolate(workflowData$config_list)
  if (is.null(config_list)) return(invisible(NULL))
  config_list$normaliseBetweenSamples$method <- normMethod
  workflowData$config_list <- config_list
  invisible(NULL)
}

runProtNormBetweenSamplesStep <- function(
  currentS4,
  normMethod,
  normData,
  proteinQcDir = NULL,
  workflowData = NULL,
  normaliseBetweenSamplesFn = normaliseBetweenSamples,
  captureCheckpointFn = .capture_checkpoint,
  saveMatrixFn = vroom::vroom_write,
  messageFn = message
) {
  messageFn("*** STEP 1: Starting between-samples normalization ***")

  recordProtNormMethod(workflowData, normMethod)

  normalized_s4 <- normaliseBetweenSamplesFn(currentS4, normalisation_method = normMethod)
  normalized_s4 <- inheritDIAExecutionSource(normalized_s4, getDIAActiveManager(),
    attr(currentS4, "dia_component_source", exact = TRUE))
  normData <- setDIAOwnedPayload(normData, "normalized_protein_obj", normalized_s4)

  captureCheckpointFn(normalized_s4, "cp05", "normalised")

  messageFn("*** STEP 1: Between-samples normalization completed ***")

  tryCatch({
    if (!is.null(proteinQcDir)) {
      saveMatrixFn(
        normalized_s4@protein_quant_table,
        file.path(proteinQcDir, "normalized_protein_matrix_pre_ruv.tsv")
      )
      messageFn("*** STEP 1b: Saved post-normalization protein matrix to protein_qc_dir ***")
    }
  }, error = function(e) {
    messageFn(paste("Warning: Could not save post-normalization matrix:", e$message))
  })

  normalized_s4
}

#' The RUV mode choices of the normalisation tab by workflow type (KIT-041): the DIA_LIMPA route computes the
#' RUV-III factors as covariates after the correlation filter (or skips RUV); every other type keeps the
#' automatic / manual / skip correction modes the UI ships.
#' @noRd
protNormRuvModeChoices <- function(workflowType = NULL) {
  if (identical(workflowType, "DIA_LIMPA")) {
    # KIT-056: the DIA route's form as a second mode; the same fit, k and controls, the matrix corrected instead.
    return(list("RUV-III factors as covariates (recommended)" = "covariates", "RUV-III correction of the quantified matrix" = "correction", "Skip RUV" = "skip"))
  }
  list("Automatic (recommended)" = "automatic", "Manual (advanced users)" = "manual", "Skip RUV (troublesome datasets)" = "skip")
}

#' The RUV-covariates action of the normalisation tab on the DIA_LIMPA route (KIT-041/053): after the correlation
#' sub-step (applied or skipped), the `apply_ruv_covariates` button runs the `protein_ruv_covariates` declaration
#' through `runDIAReferenceRuvCovariates` (or records a skip when the tab's RUV mode is `skip`) and writes the
#' tab's `ruv_covariates_summary`.
#' @noRd
runProtNormRuvCovariatesObserver <- function(input, output, workflowData, normData, experimentPaths, omicType, experimentLabel,
                                             generateRuvCorrectedQcFn, getRuvGroupingVariableFn,
                                             renderTextFn = shiny::renderText,
                                             showNotificationFn = shiny::showNotification, messageFn = message) {
  manager <- workflowData$state_manager
  if (!identical(tryCatch(manager$workflow_type, error = function(e) NULL), "DIA_LIMPA")) return(invisible(NULL))
  tryCatch({
    if (!isTRUE(normData$correlation_filtering_complete)) stop("apply or skip the correlation filter first")
    # KIT-129 fourth pass item 2: with the gate on, the RUV design section's choice gives the mode and the rules.
    choice <- omicNormRuvDesignChoice(input, omicNormRuvRunTable("proteomics", "DIA_LIMPA"))
    mode <- protNormDiaRuvDesignMode(choice) %||% input$ruv_mode %||% "covariates"
    if (identical(mode, "skip")) {
      # KIT-129 fourth pass item 3: under the gate the design none is what was done.
      omicNormSetRuvDesignResult(normData, if (identical(choice$ruv_design, "none")) omicNormRuvDesignResult("none", NA, NA, NA, list()))
      messageFn("*** RUV covariates skipped by the user (DIA with limpa) ***")
      normData$ruv_covariates_complete <- TRUE
      text <- "RUV-III covariates skipped by the user; the DA model carries no W terms."
      output$ruv_covariates_summary <- renderTextFn(text)
      showNotificationFn("RUV-III covariates skipped. Ready for differential expression analysis.", type = "message", duration = 10)
      return(invisible(NULL))
    }
    if (!mode %in% c("covariates", "correction")) stop("the RUV mode must be covariates, correction or skip")
    settings <- c(list(ruv_mode = mode, color_variable = input$color_variable, shape_variable = input$shape_variable),
      protNormDiaRuvDesignRules(choice))
    # KIT-174 (a), VAL-29: the RUV grouping warning (on this route declared technical replicates
    # are the replicate sets under every design).
    omicNormWarnRuvGrouping(workflowData, getRuvGroupingVariableFn(), mode, choice$ruv_design,
      technical = omicNormRuvHeadTechnicalReplicates(manager),
      showNotificationFn = showNotificationFn)
    finish <- function(accepted) finishProtNormRuvCovariates(accepted, mode, output, workflowData, normData, renderTextFn,
      showNotificationFn, messageFn)
    # KIT-111 D3: on a managed session the declaration runs on the worker and the delivery finishes the tab.
    runner <- workflowOmicsManagedActions("proteomics")
    if (!is.null(runner)) return(dispatchProtNormRuvCovariates(runner, settings, workflowData, normData, experimentPaths,
      generateRuvCorrectedQcFn, getRuvGroupingVariableFn, finish, function(e) reportProtNormRuvCovariatesError(e, output,
        renderTextFn, showNotificationFn, messageFn)))
    accepted <- runDIAReferenceRuvCovariates(settings, workflowData, normData, experimentPaths, omicType, experimentLabel,
      generateRuvCorrectedQcFn, getRuvGroupingVariableFn)
    finish(accepted)
  }, error = function(e) reportProtNormRuvCovariatesError(e, output, renderTextFn, showNotificationFn, messageFn))
}

#' The DIA with limpa RUV mode of the RUV design section's choice (KIT-129 fourth pass item 2): the legacy designs'
#' covariates and correction modes, `skip` for none, the covariates mode for the new designs; NULL without a choice.
#' @param choice `omicNormRuvDesignChoice`'s value, or NULL.
#' @noRd
protNormDiaRuvDesignMode <- function(choice) {
  if (is.null(choice)) return(NULL)
  switch(choice$ruv_design, legacy_correction = "correction", none = "skip", "covariates")
}

#' The DIA with limpa declared rules of the RUV design section's choice (KIT-129 fourth pass item 2), by the names
#' `diaRuvCovariatesRequestParameters` reads: the design's control rule (all proteins for the recommended design and the
#' technical-replicate design, which the route computes on the declared units either way), the k rule, the typed k, the
#' supplied list and the seed. An empty list without a new design: today's request.
#' @noRd
protNormDiaRuvDesignRules <- function(choice) {
  if (is.null(choice) || isTRUE(choice$legacy) || identical(choice$ruv_design, "none")) return(list())
  rule <- c(groups_all_features = "all_features", technical_replicates = "all_features", supplied_controls = "supplied_list",
    first_pass_nulls = "first_pass_drop_top_25")[[choice$ruv_design]]
  list(ruv_control_rule = rule, ruv_k_rule = choice$ruv_k_rule, ruv_k = choice$ruv_k, ruv_control_ids = choice$ruv_control_ids,
    ruv_k_seed = choice$ruv_k_seed, ruv_design = choice$ruv_design)
}

#' The classic normalisation's settings under the RUV design section's choice (KIT-129 fourth pass item 2): the new
#' design, its k rule and seed and its supplied list (`diaNormalizationRequest` sends them), the step's mode from the rule (the typed k under
#' `manual`, the estimated number or the separation rule under `automatic`), and the design `none` with the skip mode.
#' A legacy design, or no choice, keeps today's settings.
#' @param settings `resolveProtNormReferenceSettings`' list.
#' @param choice `omicNormRuvDesignChoice`'s value, or NULL.
#' @return `settings`.
#' @noRd
applyProtNormRuvDesign <- function(settings, choice) {
  if (is.null(choice) || isTRUE(choice$legacy)) return(settings)
  if (identical(choice$ruv_design, "none")) {
    settings$ruv_mode <- "skip"
    settings$ruv_design <- "none"
    return(settings)
  }
  settings$ruv_mode <- if (identical(choice$ruv_k_rule, "user")) "manual" else "automatic"
  if (!is.null(choice$ruv_k)) settings$ruv_k <- choice$ruv_k
  settings$ruv_design <- choice$ruv_design
  settings$ruv_k_rule <- choice$ruv_k_rule
  settings$ruv_k_seed <- choice$ruv_k_seed
  # KIT-129 fifth pass item 3: the supplied list (`diaNormalizationRequest` sends it under that design only).
  settings$ruv_control_ids <- choice$ruv_control_ids
  settings
}

#' The RUV covariates' completion on the tab (KIT-041/053/056), after the runner or the managed delivery: the RUV audit
#' trail, the summary text and the notice.
#' @param accepted The accepted reference.
#' @param mode `covariates` or `correction`.
#' @return `accepted`, invisibly.
#' @noRd
finishProtNormRuvCovariates <- function(accepted, mode, output, workflowData, normData, renderTextFn, showNotificationFn,
                                        messageFn) {
  updateProtNormRuvAuditTrail(ruvK = normData$best_k, controlGenesIndex = normData$control_genes_index,
    percentageAsNegCtrl = normData$ruv_optimization_result$best_percentage, modeLabel = mode,
    workflowData = workflowData, messageFn = messageFn)
  result <- normData$ruv_optimization_result
  controls <- paste0("negative controls ", result$best_percentage, "% (", sum(normData$control_genes_index, na.rm = TRUE), " proteins); replicate structure ", result$replicate_structure)
  text <- if (identical(mode, "correction")) paste0("RUV-III correction applied to the quantified matrix: k = ", normData$best_k, "; ", controls,
    ". The corrected matrix and the DPC EList carry Y - W alpha; the DPC standard errors are kept as quantified; no W terms enter the DA model.\n",
    "State saved as: 'ruv_corrected'") else paste0("RUV-III covariates appended to the design: k = ", normData$best_k, " (W", if (normData$best_k > 1) paste0("1..W", normData$best_k) else "1",
    "); ", controls, ". The quantified matrix is unchanged; the adjusted matrix is shown in the RUV-corrected QC plots only.\n",
    "State saved as: 'ruv_covariates'")
  output$ruv_covariates_summary <- renderTextFn(text)
  showNotificationFn(if (identical(mode, "correction")) sprintf("RUV-III correction applied (k = %s). Ready for differential expression analysis.", normData$best_k) else
    sprintf("RUV-III factors appended as covariates (k = %s). Ready for differential expression analysis.", normData$best_k), type = "message", duration = 10)
  invisible(accepted)
}

#' The RUV covariates' error on the tab: the summary text and the notice.
#' @noRd
reportProtNormRuvCovariatesError <- function(e, output, renderTextFn, showNotificationFn, messageFn) {
  message <- paste("Error computing the RUV-III covariates:", e$message)
  messageFn(message)
  output$ruv_covariates_summary <- renderTextFn(message)
  showNotificationFn(message, type = "error", duration = 10)
  invisible(NULL)
}

#' One managed action of the proteomics normalisation tab (KIT-111 D3): the operation on its named `input` (a re-run's
#' source is not the head). A state-producing action names its state; the worker's commit description replaces the
#' context's config and description.
#' @noRd
protNormManagedAction <- function(operation, parameters, input, stateName, deliver) {
  list(operation = operation, parameters = parameters, input = input, deliver = deliver,
    context = if (!is.null(stateName)) list(state_name = stateName, config = list(), description = operation))
}

#' One managed QC-figures action of the tab (KIT-122, audit E1/F3): the worker saves the QC plot products of `source`
#' and, for a stage's figures, their images (`normalization-qc-figures`), so that the tab's generators and publication
#' halves, run after its delivery, find them saved and neither hydrate the state nor draw.
#' @param source The saved state reference the figures are of.
#' @param aesthetics `list(color_var, shape_var)`.
#' @param figures `stage`, `final_pca` or `composite` (`diaQcFiguresParameters`).
#' @param deliver `function(result)`, run in the owner's scope when the worker has delivered, after the session's
#'   store has indexed what the worker saved (`indexWorkflowHeadProducts`).
#' @param composite For `composite`: `list(product, ruv_mode)` (`diaQcFiguresParameters`).
#' @noRd
protNormFiguresAction <- function(source, aesthetics, figures = "stage", deliver = function(result) invisible(NULL),
                                  composite = NULL) {
  entry <- proteomicsDiaBinding$norm$figures
  protNormManagedAction(entry$operation, entry$parameters(aesthetics, figures, composite), source, NULL, function(result) {
    entry$accept(getDIAActiveManager())
    deliver(result)
  })
}

#' The figures actions of several sources, in their order; `deliver` is the last one's.
#' @noRd
protNormFiguresActions <- function(sources, aesthetics, deliver) {
  lapply(seq_along(sources), function(i) protNormFiguresAction(sources[[i]], aesthetics,
    deliver = if (i == length(sources)) deliver else function(result) invisible(NULL)))
}

#' The tab's dispatcher of stage figures that no button's sequence owns (the pre-normalisation QC on tab entry, an
#' aesthetics change; KIT-122): NULL on a session with no runner, where the tab draws in process as before.
#' @param workflowData The tab's bus.
#' @param getPlotAestheticsFn The tab's aesthetics getter.
#' @return NULL, or `function(sources, then, failed)`, which gives FALSE when nothing was dispatched because another
#'   module's action holds the application (the caller then draws in process): `sources(manager)` gives the saved references, read in the
#'   owner's scope when the actions start (a request that waited reads what the action before it left); `then()` runs in
#'   the last delivery; `failed(error)` hears a failure.
#' @noRd
protNormFiguresDispatcher <- function(workflowData, getPlotAestheticsFn) {
  runner <- workflowOmicsManagedActions("proteomics")
  if (is.null(runner)) return(NULL)
  manager <- workflowData$state_manager
  function(sources, then, failed) {
    # Another module's action (a QC step still completing when the tab is opened) holds the application: the runner
    # cannot queue behind it, so nothing is dispatched and the caller draws in process, as an unmanaged session does.
    application <- workflowSessionApplication(shiny::getDefaultReactiveDomain())
    if (!is.null(application$active) && is.null(shiny::isolate(runner$step$active()))) return(invisible(FALSE))
    tryCatch(runner$run(function() shiny::isolate({
      references <- Filter(isSavedWorkflowReference, withDIAStateScope(manager, function() sources(manager)))
      references <- references[!duplicated(vapply(references, function(reference) reference$entity_id, ""))]
      if (!length(references)) stop("there is no saved state to draw the QC figures of")
      protNormFiguresActions(references, getPlotAestheticsFn(), function(result) then())
    }), failed = failed), error = failed)
    invisible(TRUE)
  }
}

#' The managed DIA_LIMPA pass-through (KIT-122): the worker saves the QC figures of the source the pass-through will
#' designate (`diaLimpaPassThroughSource`), then the pass-through runs in the delivery, as it does in process.
#' @param settings The run's settings (its aesthetics).
#' @param passThrough `function()`, the in-process pass-through with its notice; `failed` `function(error)`.
#' @param experimentPaths The omic's project paths: KIT-257 removes the earlier run's own files when the action starts.
#' @return TRUE, invisibly, once dispatched.
#' @noRd
dispatchProtNormPassThrough <- function(runner, settings, workflowData, normData, passThrough, failed, experimentPaths = NULL) {
  manager <- workflowData$state_manager
  held <- new.env(parent = emptyenv())
  runner$run(function() shiny::isolate({
    source <- withDIAStateScope(manager, function() {
      assertDIAExportReady(manager, allow_action = "normalization")
      source <- diaLimpaPassThroughSource(manager, normData)
      # AUDIT-KIT-257-impl I01 (KIT-244 I07): only a press that designates the head makes the earlier files stale.
      held$removed <- if (identical(source$entity_id, getDIAStateReference(manager)$entity_id))
        removeProtNormOwnFiles(experimentPaths, "normalization")
      source
    })
    list(protNormFiguresAction(source, omicNormPlotAestheticsGetter(settings)(), deliver = function(result) passThrough()))
  }), failed = function(error) failed(protNormTagRemovedFiles(error, held$removed)))
  invisible(TRUE)
}

#' The managed RUV covariates (DIA_LIMPA, KIT-111 D3): the gate, the source and the parameters as
#' `runDIAReferenceRuvCovariates` resolves them, the declaration on the worker, then the publication half and `finish`.
#' KIT-122: between the two, the worker saves the QC figures of the adjusted matrix the publication half draws.
#' @param settings `ruv_mode`, `color_variable`, `shape_variable`.
#' @param ruvQc,grouping The module's RUV QC callback and grouping function.
#' @param finish `function(accepted)`; `failed` `function(error)`.
#' @return TRUE, invisibly, once dispatched.
#' @noRd
dispatchProtNormRuvCovariates <- function(runner, settings, workflowData, normData, experimentPaths, ruvQc, grouping,
                                          finish, failed) {
  manager <- workflowData$state_manager
  group <- grouping()
  entry <- proteomicsDiaBinding$norm$ruv_covariates
  state <- if (identical(settings$ruv_mode, "correction")) "ruv_corrected" else "ruv_covariates"
  held <- new.env(parent = emptyenv())
  publish <- function(result) {
    prior_plot <- ggplot2::get_last_plot()
    on.exit(ggplot2::set_last_plot(prior_plot), add = TRUE)
    # Published before the tab's completion reads it (an argument would be a promise).
    accepted <- withDIAStateScope(manager, function() entry$publish(manager, settings, held$result[c("value", "reference")],
      workflowData, normData, experimentPaths, ruvQc))
    finish(accepted)
  }
  figures <- function() runner$run(list(protNormFiguresAction(held$result$value$adjusted,
    omicNormPlotAestheticsGetter(settings)(), deliver = publish)), failed = failed)
  # GUI audit A111-01: the gate and the source are read when the action starts, so a click that waited behind a
  # correlation re-apply runs on that re-apply's result, as in process. The step's observer does not depend on the reads.
  runner$run(function() shiny::isolate({
    request <- withDIAStateScope(manager, function() {
      assertDIAExportReady(manager, allow_action = "ruv-covariates")
      list(source = diaRuvCovariatesSource(manager, normData),
        parameters = diaRuvCovariatesRequestParameters(workflowData$config_list, group, settings))
    })
    list(protNormManagedAction(entry$operation, request$parameters, request$source, state,
      function(result) held$result <- result))
  }), done = figures, failed = failed)
  invisible(TRUE)
}

buildProtNormCompletionNotification <- function(ruvMode) {
  if (ruvMode == "covariates") {
    return("Normalisation stage ready (DIA with limpa). Apply or skip the correlation filter; the RUV-III covariates are computed after it.")
  }
  if (ruvMode == "correction") { # KIT-056
    return("Normalisation stage ready (DIA with limpa). Apply or skip the correlation filter; the RUV-III correction of the quantified matrix runs after it.")
  }
  if (ruvMode == "skip") {
    return("Normalization completed (RUV skipped)! Check the 'Post-Normalisation QC' tab for filtering summary, then proceed to 'Correlation Filtering' tab for the final step.")
  }

  "Normalization and RUV correction completed! Check the 'Post-Normalisation QC' tab for filtering summary, then proceed to 'Correlation Filtering' tab for the final step."
}

notifyProtNormNormalizationPrereqWarning <- function(
  currentState,
  showNotificationFn = shiny::showNotification,
  messageFn = message
) {
  messageFn(sprintf(
    "*** State '%s' is not valid for normalization. User needs to complete QC. ***",
    currentState
  ))
  showNotificationFn(
    "Please complete the Quality Control filtering steps before accessing normalization.",
    type = "warning",
    duration = 5
  )

  invisible(FALSE)
}

handleProtNormPreQcGenerationError <- function(
  error,
  showNotificationFn = shiny::showNotification,
  messageFn = message
) {
  messageFn(paste("*** ERROR generating pre-normalization QC:", error$message))
  showNotificationFn(
    paste("Error generating pre-normalization QC:", error$message),
    type = "error",
    duration = 10
  )

  invisible(FALSE)
}

runProtNormTabEntryWorkflow <- function(
  selectedTab,
  workflowData,
  normData,
  generatePreNormalizationQcFn,
  shouldAutoGeneratePreQcFn = shouldProtNormAutoGeneratePreQc,
  notifyInvalidStateFn = notifyProtNormNormalizationPrereqWarning,
  handlePreQcErrorFn = handleProtNormPreQcGenerationError,
  messageFn = message,
  dispatchFiguresFn = NULL
) {
  withDIAStateScope(getDIAActiveManager(), function() {
  if (is.null(selectedTab) || selectedTab != "normalization") {
    return(invisible(FALSE))
  }

  messageFn("=== NORMALIZATION TAB CLICKED ===")
  messageFn(sprintf(
    "workflow_data$state_manager is NULL: %s",
    is.null(workflowData$state_manager)
  ))

  if (is.null(workflowData$state_manager)) {
    messageFn("*** workflow_data$state_manager is NULL - cannot check state ***")
    return(invisible(FALSE))
  }

  current_state <- workflowData$state_manager$current_state

  messageFn(sprintf("Current state: '%s'", current_state))
  messageFn(sprintf("Target trigger state: 'protein_replicate_filtered'"))
  messageFn(sprintf("Pre-norm QC generated: %s", normData$pre_norm_qc_generated))

  if (!shouldAutoGeneratePreQcFn(
    selectedTab,
    current_state,
    normData$pre_norm_qc_generated
  )) {
    # KIT-101 F-10: the gate is also FALSE when the pre-normalisation QC was already generated, and on a route whose
    # last QC state is not a trigger state; neither means QC is incomplete, so warn only while it is.
    if (shouldAutoGeneratePreQcFn(selectedTab, current_state, FALSE) ||
        identical(workflowData$tab_status$quality_control, "complete")) {
      return(invisible(FALSE))
    }
    # KIT-116 F116-01: on DIA_LIMPA the protein QC's last state and the RUV covariates after it are reached only once
    # the QC has finished, so they need no status (a restored session, or the DA tab's load-session, may not have one).
    if (identical(workflowData$state_manager$workflow_type, "DIA_LIMPA") &&
        current_state %in% c("peptide_count_filtered", "ruv_covariates")) {
      return(invisible(FALSE))
    }
    return(notifyInvalidStateFn(currentState = current_state))
  }

  messageFn("*** AUTO-TRIGGERING PRE-NORMALIZATION QC (chunk 24) ***")

  tryCatch({
    # KIT-122: on a managed session the worker saves the figures first and the generator runs in its delivery. The
    # QC is marked generated when dispatched, as the shared frame's pre-QC is, and unmarked if the worker fails; the
    # image slots are refreshed when the files exist (in process they are written before the slots first render).
    if (is.function(dispatchFiguresFn)) {
      normData$pre_norm_qc_generated <- TRUE
      dispatched <- dispatchFiguresFn(
        sources = function(manager) list(resolveProtNormQcStateObject(manager, messageFn = messageFn)$currentS4),
        then = function() {
          generatePreNormalizationQcFn()
          normData$plot_refresh_trigger <- normData$plot_refresh_trigger + 1
        },
        failed = function(e) {
          normData$pre_norm_qc_generated <- FALSE
          handlePreQcErrorFn(error = e)
        })
      if (!isFALSE(dispatched)) return(TRUE)
    }
    generatePreNormalizationQcFn()
    normData$pre_norm_qc_generated <- TRUE

    messageFn("*** PRE-NORMALIZATION QC COMPLETED SUCCESSFULLY ***")
    messageFn(sprintf("State is '%s' and PCA already generated. skipping.", current_state))

    TRUE
  }, error = function(e) {
    handlePreQcErrorFn(error = e)
  })
  })
}

handleProtNormNormalizationError <- function(
  error,
  showNotificationFn = shiny::showNotification,
  messageFn = message
) {
  messageFn(paste("Error in normalization workflow:", error$message))
  showNotificationFn(
    paste("Error in normalization:", error$message),
    type = "error",
    duration = 10
  )

  invisible(FALSE)
}

# The plain settings the reference normalisation declaration reads (KIT-004): the module edge is the
# only place GUI inputs are read.
resolveProtNormReferenceSettings <- function(input) {
  fields <- c("norm_method", "ruv_mode", "auto_percentage_min", "auto_percentage_max",
    "separation_metric", "k_penalty_weight", "adaptive_k_penalty", "ruv_percentage", "ruv_k",
    "color_variable", "shape_variable")
  stats::setNames(lapply(fields, function(field) input[[field]]), fields)
}

runProtNormNormalizationWorkflow <- function(
  input,
  workflowData,
  normData,
  experimentPaths,
  omicType,
  experimentLabel,
  generatePostNormalizationQcFn,
  generateRuvCorrectedQcFn,
  getRuvGroupingVariableFn,
  buildCompletionNotificationFn = buildProtNormCompletionNotification,
  showNotificationFn = shiny::showNotification,
  handleErrorFn = handleProtNormNormalizationError,
  messageFn = message,
  generatePreNormalizationQcFn = NULL
) {
  # AUDIT-KIT-257-impl I05: a run that stopped after removing its earlier files says so in its error notice.
  fail <- function(e) handleErrorFn(error = protNormRemovedFilesNotice(e))
  return(tryCatch({
    settings <- resolveProtNormReferenceSettings(input)
    if (identical(workflowData$state_manager$workflow_type, "DIA_LIMPA")) {
      # KIT-182 (VAL-91 D1): a table to which quality control was applied after normalisation (a session saved before
      # the fix) is not designated.
      withDIAStateScope(workflowData$state_manager, function()
        assertDIALimpaDesignationNotNormalised(workflowData$state_manager, normData))
      # KIT-041: the route's normalisation ran in the protein QC stage; RUV covariates follow the correlation step.
      # KIT-259 (B2): after the completion notice, one warning names the figures that could not be drawn.
      pass_through <- function() withProtNormFigureNotice({
        completed <- runDIAReferenceLimpaNormalisationPassThrough(settings, workflowData, normData, generatePostNormalizationQcFn)
        showNotificationFn(buildCompletionNotificationFn(settings$ruv_mode %||% "covariates"), type = "message", duration = 10)
        completed
      }, showNotificationFn)
      # KIT-122: on a managed session the worker saves the source's QC figures before the pass-through draws them.
      runner <- workflowOmicsManagedActions("proteomics")
      if (is.null(runner)) {
        # KIT-257 (VAL-250): the earlier run's own files go after the route's refusals, before the pass-through (whose
        # definition is pinned, test-prot-norm-core.R); managed, the dispatched action removes them when it starts.
        removed <- withDIAStateScope(workflowData$state_manager, function() {
          manager <- workflowData$state_manager
          assertDIAExportReady(manager, allow_action = "normalization")
          # AUDIT-KIT-257-impl I01 (KIT-244 I07): only a press that designates the head makes the earlier correlation's
          # table and RUV figures stale; pressed again on the tab's correlation or RUV result, the pass-through keeps the
          # input it designated before (P-35, KIT-244), so that result and its files stay current.
          if (identical(diaLimpaPassThroughSource(manager, normData)$entity_id, getDIAStateReference(manager)$entity_id))
            removeProtNormOwnFiles(experimentPaths, "normalization")
        })
        return(withProtNormRemovedFiles(removed, pass_through()))
      }
      return(dispatchProtNormPassThrough(runner, settings, workflowData, normData, pass_through, fail, experimentPaths))
    }
    # KIT-182: the request runs synchronously here on both paths; when the run reads another table than the previous
    # run, its classed message is shown as the tab's notice.
    withCallingHandlers({
      # KIT-129 fourth pass item 2: with the gate on, the RUV design section's choice (the classic routes' table).
      route <- workflowData$state_manager$workflow_type
      choice <- omicNormRuvDesignChoice(input, omicNormRuvRunTable("proteomics", route))
      settings <- applyProtNormRuvDesign(settings, choice)
      # KIT-174 (a), VAL-29: the RUV grouping warning, once the mode and the design are known.
      omicNormWarnRuvGrouping(workflowData, getRuvGroupingVariableFn(), settings$ruv_mode,
        choice$ruv_design, classic = TRUE, showNotificationFn = showNotificationFn,
        technical = omicNormRuvHeadTechnicalReplicates(workflowData$state_manager))
      recordProtNormMethod(workflowData, settings$norm_method)
      finish <- function(completed) finishProtNormNormalization(completed, settings, workflowData, normData,
        buildCompletionNotificationFn, showNotificationFn, messageFn)
      # KIT-111 D3: on a managed session the declaration runs on the worker and the delivery finishes the tab.
      runner <- workflowOmicsManagedActions("proteomics")
      if (!is.null(runner)) dispatchProtNormNormalization(runner, settings, workflowData, normData, experimentPaths,
        omicType, experimentLabel, generatePostNormalizationQcFn, generateRuvCorrectedQcFn, getRuvGroupingVariableFn,
        finish, fail, preQc = generatePreNormalizationQcFn) else
        # KIT-259 (B2): after the completion notice, one warning names the figures that could not be drawn.
        withProtNormFigureNotice(finish(runDIAReferenceNormalization(settings, workflowData, normData, experimentPaths,
          omicType, experimentLabel, generatePostNormalizationQcFn, generateRuvCorrectedQcFn,
          getRuvGroupingVariableFn, pre_qc = generatePreNormalizationQcFn)), showNotificationFn)
    }, prot_norm_input_current = function(condition) {
      showNotificationFn(conditionMessage(condition), type = "message", duration = 10)
      invokeRestart("muffleMessage")
    })
  }, error = fail))
}

#' The normalisation's completion on the tab, after the runner or the managed delivery: the RUV audit trail (unless
#' skipped) and the completion notice.
#' @param completed The publication half's value.
#' @return `completed`.
#' @noRd
finishProtNormNormalization <- function(completed, settings, workflowData, normData, buildCompletionNotificationFn,
                                        showNotificationFn, messageFn) {
  if (isTRUE(completed) && !identical(settings$ruv_mode, "skip")) {
    updateProtNormRuvAuditTrail(
      ruvK = normData$best_k,
      controlGenesIndex = normData$control_genes_index,
      percentageAsNegCtrl = if (identical(settings$ruv_mode, "automatic"))
        normData$ruv_optimization_result$best_percentage else settings$ruv_percentage,
      modeLabel = settings$ruv_mode,
      workflowData = workflowData,
      messageFn = messageFn
    )
  }
  showNotificationFn(buildCompletionNotificationFn(settings$ruv_mode), type = "message", duration = 10)
  warnProtNormScaleNegativeLogs(completed, settings, workflowData, normData, showNotificationFn, messageFn)
  completed
}

#' KIT-177 D6: `scale` divides each sample by the median of its absolute log2 values, so on plain log2 values below 0
#' its result depends on the table's unit. When `scale` ran and the normalisation's input protein table holds log2
#' values below 0: one step-log line and one warning notification with their count. `scale` itself is unchanged.
#' @return The count, invisibly (NA when not counted).
#' @noRd
warnProtNormScaleNegativeLogs <- function(completed, settings, workflowData, normData, showNotificationFn, messageFn) {
  if (!isTRUE(completed) || !identical(settings$norm_method, "scale")) return(invisible(NA_integer_))
  manager <- shiny::isolate(workflowData$state_manager)
  reference <- shiny::isolate(normData$normalization_input_reference)
  if (is.null(manager) || is.null(reference)) return(invisible(NA_integer_))
  count <- withDIAStateScope(manager, function() {
    object <- readDIAReferenceProteinObject(manager, reference)
    table <- object@protein_quant_table
    values <- table[setdiff(names(table), object@protein_id_column)]
    sum(vapply(values, function(x) if (is.numeric(x)) sum(x < 0, na.rm = TRUE) else 0L, integer(1)))
  })
  if (count > 0L) {
    texts <- omicGuiProtLog2Texts()
    messageFn(sprintf(texts[["scale_log"]], count))
    showNotificationFn(sprintf(texts[["scale_warning"]], count), type = "warning", duration = NULL)
  }
  invisible(count)
}

#' The managed normalisation (KIT-111 D3): the gate and the request (`diaNormalizationRequest`, which records the
#' source) in the controller, `normalization-ruv` on the worker on that source, then the publication half and `finish`.
#' KIT-122: between the two, the worker saves the stage figures the publication half draws: the normalised state's and,
#' unless RUV was skipped, the final's; then it builds the composite QC figure from them and saves its files, which the
#' publication half writes (`writeDIACompositeFigureFiles`).
#' @param postQc,ruvQc,grouping The module's QC callbacks and grouping function.
#' @param finish `function(completed)`; `failed` `function(error)`.
#' @param preQc KIT-182 (F-e): the tab's pre-normalisation figure generator. When the request's source differs from the
#'   table the figures were drawn from (`redraw_pre`), the actions start with the source's figures, whose delivery calls
#'   it; NULL dispatches no such action, as before.
#' @return TRUE, invisibly, once dispatched.
#' @noRd
dispatchProtNormNormalization <- function(runner, settings, workflowData, normData, experimentPaths, omicType,
                                          experimentLabel, postQc, ruvQc, grouping, finish, failed, preQc = NULL) {
  manager <- workflowData$state_manager
  request <- withDIAStateScope(manager, function() {
    assertDIAExportReady(manager, allow_action = "normalization")
    diaNormalizationRequest(manager, settings, normData, grouping())
  })
  entry <- proteomicsDiaBinding$norm$run
  state <- if (identical(settings$ruv_mode, "skip")) "normalized" else "ruv_corrected"
  held <- new.env(parent = emptyenv())
  # AUDIT-KIT-257-impl I05: a failure after the removal carries the removed files, for the tab's notice.
  stopped <- function(error) failed(protNormTagRemovedFiles(error, held$removed))
  publish <- function(result) {
    prior_plot <- ggplot2::get_last_plot()
    on.exit(ggplot2::set_last_plot(prior_plot), add = TRUE)
    # KIT-259 (B2): after the completion notice, one warning names the figures that could not be drawn.
    withProtNormFigureNotice({
      completed <- withDIAStateScope(manager, function() entry$publish(manager, settings, held$result[c("value", "reference")],
        workflowData, normData, experimentPaths, omicType, experimentLabel, postQc, ruvQc))
      finish(completed)
    })
  }
  figures <- function() {
    saved <- held$result$value
    aesthetics <- omicNormPlotAestheticsGetter(settings)()
    sources <- c(list(saved$normalized), if (!identical(settings$ruv_mode, "skip")) list(saved$final))
    composite <- protNormFiguresAction(saved$final, aesthetics, "composite", publish,
      list(product = held$result$reference$entity_id, ruv_mode = settings$ruv_mode))
    runner$run(c(protNormFiguresActions(sources, aesthetics, function(result) invisible(NULL)), list(composite)),
      failed = stopped)
  }
  actions <- list(protNormManagedAction(entry$operation, request$parameters, request$source, state,
    function(result) held$result <- result))
  if (isTRUE(request$redraw_pre) && is.function(preQc)) actions <- c(list(protNormFiguresAction(request$source,
    omicNormPlotAestheticsGetter(settings)(), deliver = function(result) preQc())), actions)
  # KIT-257 (VAL-250, R1-01): the earlier run's own files are removed when the actions start, not at the click: a click
  # that waited behind another action (a correlation still publishing its table) removes what that action wrote.
  runner$run(function() {
    held$removed <- removeProtNormOwnFiles(experimentPaths, "normalization")
    actions
  }, done = figures, failed = stopped)
  invisible(TRUE)
}

runProtNormApplyCorrelationObserver <- function(
  input,
  output,
  workflowData,
  normData,
  experimentPaths,
  omicType,
  experimentLabel,
  getRuvGroupingVariableFn,
  getPlotAestheticsFn,
  checkMemoryUsageFn,
  resolveCorrelationInputObjectFn = resolveProtNormCorrelationInputObject,
  runApplyCorrelationWorkflowFn = runProtNormApplyCorrelationWorkflow,
  completeCorrelationWorkflowFn = completeProtNormCorrelationWorkflow,
  handleCorrelationErrorFn = handleProtNormCorrelationError,
  messageFn = message,
  generateRuvCorrectedQcFn = NULL
) {
  failed <- function(e) handleCorrelationErrorFn(error = protNormRemovedFilesNotice(e),
    logPrefix = "Error in correlation filtering:", notificationPrefix = "Error in correlation filtering:")
  # KIT-111 D3: on a managed session the pairs and the selection run on the worker and the delivery finishes the tab.
  runner <- workflowOmicsManagedActions("proteomics")
  threshold <- input$min_pearson_correlation_threshold
  if (!is.null(runner)) return(tryCatch(dispatchProtNormCorrelation(runner, workflowData, normData, experimentPaths,
    threshold, omicType, experimentLabel, getRuvGroupingVariableFn, getPlotAestheticsFn,
    function(final) completeCorrelationWorkflowFn(
      finalS4ForDe = final,
      workflowData = workflowData,
      output = output,
      normData = normData,
      correlationThreshold = threshold,
      skipped = FALSE,
      successNotification = "Correlation filtering completed! Ready for differential expression analysis.",
      completionMessage = "Correlation workflow completed",
      messagePrefix = "Correlation"
    ), failed), error = failed))
  return(tryCatch(withDIAStateScope(getDIAActiveManager(), function() {
    final <- runDIAReferenceCorrelation(workflowData, normData, experimentPaths,
      input$min_pearson_correlation_threshold, FALSE, omicType, experimentLabel, getRuvGroupingVariableFn, getPlotAestheticsFn)
    completeCorrelationWorkflowFn(
      finalS4ForDe = final,
      workflowData = workflowData,
      output = output,
      normData = normData,
      correlationThreshold = input$min_pearson_correlation_threshold,
      skipped = FALSE,
      successNotification = "Correlation filtering completed! Ready for differential expression analysis.",
      completionMessage = "Correlation workflow completed",
      messagePrefix = "Correlation"
    )
    # KIT-041/053: on the DIA_LIMPA route the RUV-III covariates are the tab's own action after the filter
    # (`apply_ruv_covariates`), the script's order; nothing runs here.
    final
  }), error = function(e) handleCorrelationErrorFn(error = protNormRemovedFilesNotice(e),
      logPrefix = "Error in correlation filtering:", notificationPrefix = "Error in correlation filtering:")))
}

#' The managed correlation filter (KIT-111 D3): the gate and the source (the normalisation's final, as
#' `runDIAReferenceCorrelation` reads it), then `pearson-pairs` and `pearson-selection` on that source on the worker, the
#' selection naming the pairs product; then the publication half and `finish`.
#' @param threshold The filter's threshold.
#' @param grouping,aesthetics The module's grouping and aesthetics functions.
#' @param finish `function(final)`; `failed` `function(error)`.
#' @return TRUE, invisibly, once dispatched.
#' @noRd
dispatchProtNormCorrelation <- function(runner, workflowData, normData, experimentPaths, threshold, omicType,
                                        experimentLabel, grouping, aesthetics, finish, failed) {
  manager <- workflowData$state_manager
  group <- grouping()
  entry <- proteomicsDiaBinding$norm$correlation
  held <- new.env(parent = emptyenv())
  stopped <- function(error) failed(protNormTagRemovedFiles(error, held$removed))
  # KIT-122: after the selection the worker saves the PCA product of its checkpoint, which the publication half draws.
  figures <- function() runner$run(list(protNormFiguresAction(held$final, aesthetics(), "final_pca", function(result) {
    run <- list(pairs = held$pairs$reference, final = held$final)
    final <- withDIAStateScope(manager, function() entry$publish(manager, run, workflowData, normData, experimentPaths,
      threshold, FALSE, omicType, experimentLabel, aesthetics))
    finish(final)
  })), failed = stopped)
  select <- function() runner$run(list(protNormManagedAction(entry$selection,
    list(pairs = held$pairs$reference$entity_id, threshold = threshold), held$source, "correlation_filtered",
    function(result) held$final <- result$value)), done = figures, failed = stopped)
  # GUI audit A111-01: the gate and the source are read when the pairs start, so a click that waited behind a
  # normalisation re-run runs on that re-run's final, as in process. The step's observer does not depend on the reads.
  runner$run(function() shiny::isolate({
    held$source <- withDIAStateScope(manager, function() {
      assertDIAExportReady(manager, allow_action = "correlation")
      source <- diaCorrelationSource(manager, normData)
      # KIT-257 (VAL-250): the earlier results table (either name) goes before the pairs start.
      held$removed <- removeProtNormOwnFiles(experimentPaths, "correlation")
      readDIAStateMetadata(manager, source)
      source
    })
    list(protNormManagedAction(entry$pairs, diaPearsonPairsParameters(group), held$source, NULL,
      function(result) held$pairs <- result))
  }), done = select, failed = stopped)
  invisible(TRUE)
}

runProtNormSkipCorrelationObserver <- function(
  output,
  workflowData,
  normData,
  experimentPaths,
  omicType,
  experimentLabel,
  getPlotAestheticsFn,
  checkMemoryUsageFn,
  resolveCorrelationInputObjectFn = resolveProtNormCorrelationInputObject,
  runSkipCorrelationWorkflowFn = runProtNormSkipCorrelationWorkflow,
  completeCorrelationWorkflowFn = completeProtNormCorrelationWorkflow,
  handleCorrelationErrorFn = handleProtNormCorrelationError,
  messageFn = message,
  input = NULL,
  getRuvGroupingVariableFn = NULL,
  generateRuvCorrectedQcFn = NULL
) {
  return(tryCatch(withDIAStateScope(getDIAActiveManager(), function() {
    final <- runDIAReferenceCorrelation(workflowData, normData, experimentPaths,
      0, TRUE, omicType, experimentLabel, function() NULL, getPlotAestheticsFn)
    completeCorrelationWorkflowFn(
      finalS4ForDe = final,
      workflowData = workflowData,
      output = output,
      normData = normData,
      correlationThreshold = 0,
      skipped = TRUE,
      successNotification = "Correlation filtering skipped! Ready for differential expression analysis.",
      completionMessage = "Correlation workflow completed",
      messagePrefix = "Correlation"
    )
    # KIT-041/053: the RUV-III covariates follow the (skipped) filter through the tab's own action.
    final
  }), error = function(e) handleCorrelationErrorFn(error = protNormRemovedFilesNotice(e),
      logPrefix = "Error in correlation filtering:", notificationPrefix = "Error in correlation filtering:")))
}

notifyProtNormExportSessionPrereqWarning <- function(
  showNotificationFn = shiny::showNotification,
  messageFn = message
) {
  messageFn("Export skipped because correlation filtering is not complete.")
  showNotificationFn(
    "Please complete correlation filtering before exporting session data.",
    type = "warning",
    duration = 5
  )

  invisible(FALSE)
}

runProtNormExportObserver <- function(
  input,
  workflowData,
  normData,
  experimentPaths,
  canExportFilteredSessionFn = canProtNormExportFilteredSession,
  notifyExportPrereqFn = notifyProtNormExportSessionPrereqWarning,
  resolveExportSourceDirFn = resolveProtNormExportSourceDir,
  runExportSessionWorkflowFn = runProtNormExportSessionWorkflow,
  handleExportErrorFn = handleProtNormExportError,
  showNotificationFn = shiny::showNotification,
  messageFn = message
) {
  withDIAStateScope(getDIAActiveManager(), function() {
  messageFn("=== EXPORT FILTERED SESSION BUTTON CLICKED ===")

  if (!canExportFilteredSessionFn(
    correlationFilteringComplete = normData$correlation_filtering_complete,
    correlationFilteredObj = getDIAOwnedInput(normData, "correlation_filtered_obj")
  )) {
    return(notifyExportPrereqFn())
  }

  tryCatch({
    # KIT-206: refused, before anything is written or left pending, when the current table is not the tab's last
    # result. Here and not in the workflow: KIT-182's F1 (h) exports a stale head through the workflow on purpose.
    assertProtNormExportHeadCurrent(if (is.list(workflowData) || is.environment(workflowData))
      workflowData$state_manager, normData$correlation_filtered_obj)
    source_dir <- resolveExportSourceDirFn(experimentPaths)

    export_result <- runExportSessionWorkflowFn(
      workflowData = workflowData,
      normData = normData,
      input = input,
      sourceDir = source_dir
    )

    showNotificationFn(
      sprintf(
        "Filtered session data exported successfully!\nSaved as: %s\nSee summary file for details.",
        export_result$exportArtifacts$sessionFilename
      ),
      type = "message",
      duration = 10
    )

    messageFn("=== EXPORT FILTERED SESSION COMPLETED SUCCESSFULLY ===")

    export_result
  }, error = function(e) {
    handleExportErrorFn(e)
  })
  })
}

runProtNormResetObserver <- function(
  workflowData,
  normData,
  output,
  runResetWorkflowFn = runProtNormResetWorkflow,
  handleResetErrorFn = handleProtNormResetError,
  messageFn = message
) {
  withDIAStateScope(getDIAActiveManager(), function() {
  messageFn("Resetting normalization...")

  tryCatch({
    runResetWorkflowFn(
      workflowData = workflowData,
      normData = normData,
      output = output
    )
  }, error = function(e) {
    handleResetErrorFn(e)
  })
  })
}
