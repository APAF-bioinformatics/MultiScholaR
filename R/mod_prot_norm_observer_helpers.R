registerProtNormServerObservers <- function(
  input,
  output,
  session,
  selectedTab = NULL,
  workflowData,
  normData,
  experimentPaths,
  omicType,
  experimentLabel,
  generatePreNormalizationQcFn,
  generatePostNormalizationQcFn,
  generateRuvCorrectedQcFn,
  getPlotAestheticsFn,
  getRuvGroupingVariableFn,
  checkMemoryUsageFn = checkProtNormMemoryUsage,
  updateDesignDrivenChoicesFn = updateProtNormDesignDrivenChoices,
  runTabEntryWorkflowFn = runProtNormTabEntryWorkflow,
  regenerateQcForAestheticChangeFn = regenerateProtNormQcForAestheticChange,
  runNormalizationWorkflowFn = runProtNormNormalizationWorkflow,
  runApplyCorrelationObserverFn = runProtNormApplyCorrelationObserver,
  runSkipCorrelationObserverFn = runProtNormSkipCorrelationObserver,
  runResetObserverFn = runProtNormResetObserver,
  runExportObserverFn = runProtNormExportObserver,
  observeFn = shiny::observe,
  observeEventFn = shiny::observeEvent,
  updateRadioButtonsFn = shiny::updateRadioButtons,
  renderUIFn = shiny::renderUI,
  runRuvCovariatesObserverFn = runProtNormRuvCovariatesObserver
) {
  observeFn({
    updateDesignDrivenChoicesFn(session, workflowData$design_matrix)
  })

  if (!is.null(selectedTab)) {
    observeEventFn(selectedTab(), {
      runTabEntryWorkflowFn(
        selectedTab = selectedTab(),
        workflowData = workflowData,
        normData = normData,
        generatePreNormalizationQcFn = generatePreNormalizationQcFn
      )
    }, ignoreInit = TRUE)
  }

  observeEventFn(c(input$color_variable, input$shape_variable), {
    regenerateQcForAestheticChangeFn(
      normData = normData,
      generatePreNormalizationQcFn = generatePreNormalizationQcFn,
      generatePostNormalizationQcFn = generatePostNormalizationQcFn,
      generateRuvCorrectedQcFn = generateRuvCorrectedQcFn
    )
  })

  observeEventFn(input$run_normalization, {
    runNormalizationWorkflowFn(
      input = input,
      workflowData = workflowData,
      normData = normData,
      experimentPaths = experimentPaths,
      omicType = omicType,
      experimentLabel = experimentLabel,
      generatePostNormalizationQcFn = generatePostNormalizationQcFn,
      generateRuvCorrectedQcFn = generateRuvCorrectedQcFn,
      getRuvGroupingVariableFn = getRuvGroupingVariableFn
    )
  })

  # KIT-041: on the DIA_LIMPA route the RUV-III factors are covariates computed after the correlation filter (the
  # limpa script's order); the tab's RUV mode offers `covariates` or `skip`, set here once the type is known.
  observeFn({
    workflowData$state_update_trigger
    manager <- workflowData$state_manager
    if (identical(tryCatch(manager$workflow_type, error = function(e) NULL), "DIA_LIMPA") &&
        !(shiny::isolate(input$ruv_mode) %||% "") %in% c("covariates", "correction", "skip")) { # KIT-056: correction
      updateRadioButtonsFn(session, "ruv_mode", choices = protNormRuvModeChoices("DIA_LIMPA"), selected = "covariates")
    }
  })

  # KIT-053: the RUV covariates panel exists on the DIA_LIMPA route only (rendered once the type is known).
  output$ruv_covariates_panel <- renderUIFn({
    workflowData$state_update_trigger
    manager <- workflowData$state_manager
    if (!identical(tryCatch(manager$workflow_type, error = function(e) NULL), "DIA_LIMPA")) return(NULL)
    ns <- session$ns
    shiny::tagList(
      shiny::br(),
      shiny::h5("RUV-III Factors as Covariates (DIA with limpa)"),
      shiny::helpText("After the correlation filter the negative-control sweep chooses k and the controls. Covariates: the W factors are appended to the design for the DA model and the quantified matrix is not corrected. Correction: the quantified matrix and the DPC EList are corrected (Y - W alpha) and the DA model carries no W terms."),
      shiny::actionButton(ns("apply_ruv_covariates"), "Compute RUV-III Covariates", class = "btn-primary", width = "100%", icon = shiny::icon("play")),
      shiny::br(), shiny::br(),
      shiny::verbatimTextOutput(ns("ruv_covariates_summary")))
  })

  observeEventFn(input$apply_ruv_covariates, {
    runRuvCovariatesObserverFn(
      input = input,
      output = output,
      workflowData = workflowData,
      normData = normData,
      experimentPaths = experimentPaths,
      omicType = omicType,
      experimentLabel = experimentLabel,
      generateRuvCorrectedQcFn = generateRuvCorrectedQcFn,
      getRuvGroupingVariableFn = getRuvGroupingVariableFn
    )
  })

  observeEventFn(input$apply_correlation_filter, {
    runApplyCorrelationObserverFn(
      input = input,
      output = output,
      workflowData = workflowData,
      normData = normData,
      experimentPaths = experimentPaths,
      omicType = omicType,
      experimentLabel = experimentLabel,
      getRuvGroupingVariableFn = getRuvGroupingVariableFn,
      getPlotAestheticsFn = getPlotAestheticsFn,
      checkMemoryUsageFn = checkMemoryUsageFn,
      generateRuvCorrectedQcFn = generateRuvCorrectedQcFn
    )
  })

  observeEventFn(input$skip_correlation_filter, {
    runSkipCorrelationObserverFn(
      output = output,
      workflowData = workflowData,
      normData = normData,
      experimentPaths = experimentPaths,
      omicType = omicType,
      experimentLabel = experimentLabel,
      getPlotAestheticsFn = getPlotAestheticsFn,
      checkMemoryUsageFn = checkMemoryUsageFn,
      input = input,
      getRuvGroupingVariableFn = getRuvGroupingVariableFn,
      generateRuvCorrectedQcFn = generateRuvCorrectedQcFn
    )
  })

  observeEventFn(input$reset_normalization, {
    runResetObserverFn(
      workflowData = workflowData,
      normData = normData,
      output = output,
      ruvMode = input$ruv_mode,
      groupingVariable = getRuvGroupingVariableFn()
    )
  })

  observeEventFn(input$export_filtered_session, {
    runExportObserverFn(
      input = input,
      workflowData = workflowData,
      normData = normData,
      experimentPaths = experimentPaths
    )
  })

  invisible(TRUE)
}

#' KIT-120: proteomics' registrations on the shared normalisation frame (`omicGuiConfig("proteomics")$norm$server`;
#' mod_omic_norm.R calls it with the frame's state and callbacks). It is the body of the tab's own server shell
#' (`runProtNormModuleServerShell`, removed by KIT-120 W5) on the frame's state: the tab's own fields, the managed actions, the S4 download, the restore, the QC generators, the
#' observers (`registerProtNormServerObservers`) and the outputs (`registerProtNormRenderOutputs`). The frame registers
#' the design's choices and the stage grid's images itself, so neither is registered here; the RUV QC tab's content
#' (the frame's `ruv_qc_ui`) is proteomics' plot, summary and table under their own ids.
#' @param input,output,session The frame's module session.
#' @param context The frame's `workflow_data`, `experiment_paths`, `omic_type`, `experiment_label`, `selected_tab`,
#'   `norm_data` and `plot_aesthetics`.
#' @return The tab's state, invisibly.
#' @noRd
protNormFrameServer <- function(input, output, session, context) {
  workflowData <- context$workflow_data
  experimentPaths <- context$experiment_paths
  omicType <- context$omic_type
  experimentLabel <- context$experiment_label
  normData <- context$norm_data
  fields <- protNormStateFields()
  for (field in names(fields)) normData[[field]] <- fields[[field]]
  # The tab's own state has no QC plot paths until its first figures (the frame's starts as an empty list).
  normData$qc_plot_paths <- NULL
  if (inherits(session, c("ShinySession", "session_proxy", "MockShinySession"))) {
    installWorkflowOmicsActions(session, workflowData, "proteomics", experimentLabel, experimentPaths)
  }
  output$export_s4_compatibility <- protNormS4DownloadHandler(workflowData, normData, experimentPaths, omicType, experimentLabel)
  shiny::observeEvent(workflowData$native_session_restored, {
    restoreProtNormTabState(session, workflowData, normData)
  }, ignoreInit = TRUE)
  getPlotAesthetics <- context$plot_aesthetics
  getRuvGroupingVariable <- function() getProtNormRuvGroupingVariable(input$ruv_grouping_variable)
  generators <- protNormQcGenerators(workflowData, normData, experimentPaths, getPlotAesthetics)
  registerProtNormServerObservers(
    input = input, output = output, session = session, selectedTab = context$selected_tab,
    workflowData = workflowData, normData = normData, experimentPaths = experimentPaths,
    omicType = omicType, experimentLabel = experimentLabel,
    generatePreNormalizationQcFn = generators$pre, generatePostNormalizationQcFn = generators$post,
    generateRuvCorrectedQcFn = generators$ruv,
    getPlotAestheticsFn = getPlotAesthetics, getRuvGroupingVariableFn = getRuvGroupingVariable,
    # The frame registers the design's choices (omicNormDesignChoices: the same choices and selections).
    updateDesignDrivenChoicesFn = function(session, designMatrix) invisible(NULL)
  )
  ns <- session$ns
  output$ruv_qc_ui <- shiny::renderUI(shiny::fluidRow(
    shiny::column(8, shinyjqui::jqui_resizable(shiny::plotOutput(ns("ruv_canonical_correlation_plot"), height = "500px"))),
    shiny::column(4, shiny::wellPanel(
      shiny::h5("Optimization Summary"), shiny::verbatimTextOutput(ns("ruv_optimization_summary")), shiny::br(),
      shiny::h5("Optimization Results"), shiny::helpText("Detailed results for all tested percentages:"),
      DT::dataTableOutput(ns("ruv_optimization_table"))))))
  registerProtNormRenderOutputs(output = output, normData = normData, ruvMode = input$ruv_mode,
    groupingVariable = getRuvGroupingVariable())
  invisible(normData)
}

#' The S4 compatibility download of proteomics' normalisation tab (`export_s4_compatibility`), for its own shell and
#' for the shared frame (KIT-120).
#' @noRd
protNormS4DownloadHandler <- function(workflowData, normData, experimentPaths, omicType, experimentLabel) {
  shiny::downloadHandler(
    filename = function() paste0(omicType, "_", experimentLabel, "_final_s4.RDS"),
    content = function(file) {
      names <- resolveProtNormCorrelationResultFilenames(normData$ruv_optimization_result,
        workflowData$ruv_optimization_result)
      canonical <- file.path(experimentPaths$protein_qc_dir, names$rdsFilename)
      exportDIAS4Compatibility(workflowData$state_manager, canonical)
      if (!file.copy(canonical, file, overwrite = TRUE)) stop("Could not deliver compatibility export")
      if (!is.null(experimentPaths$integration_dir)) {
        destination <- file.path(experimentPaths$integration_dir,
          paste0(omicType, "_", experimentLabel, "_final_s4.RDS"))
        if (!file.copy(canonical, destination, overwrite = TRUE)) stop("Could not publish compatibility Integration file")
      }
    }
  )
}

#' A restored native session's normalisation tab state (the `native_session_restored` observer's body), for the tab's
#' own shell and for the shared frame (KIT-120).
#' @noRd
restoreProtNormTabState <- function(session, workflowData, normData) {
  manager <- workflowData$state_manager
  withDIAStateScope(manager, function() {
    context <- manager$.context
    extra <- context$native_extra_references %||% list()
    metadata <- if (!is.null(context$native_metadata_reference))
      resolveDIAStateValue(context, context$native_metadata_reference) else list()
    normData$normalization_input_reference <- extra$normalization_input_reference %||%
      getDIAStateReference(manager, "protein_replicate_filtered")
    normData$normalized_protein_obj <- extra$normalized_protein_obj %||% getDIAStateReference(manager, "normalized")
    normData$ruv_normalized_obj <- extra$ruv_normalized_obj %||%
      getDIAStateReference(manager, "ruv_corrected") %||% getDIAStateReference(manager, "normalized")
    normData$correlation_filtered_obj <- getDIAStateReference(manager, "correlation_filtered")
    # GUI audit P-35 (review): a restored session's RUV covariates start from its own correlation result.
    normData$ruv_covariates_input_reference <- NULL
    active_history <- diaStateHistory(manager)
    normalized_active <- any(c("normalized", "ruv_corrected", "correlation_filtered") %in% active_history)
    normData$normalization_complete <- normalized_active && !is.null(normData$normalized_protein_obj)
    normData$ruv_complete <- normalized_active && !is.null(normData$ruv_normalized_obj)
    normData$correlation_filtering_complete <- "correlation_filtered" %in% active_history &&
      !is.null(normData$correlation_filtered_obj)
    if (!normalized_active) normData$normalized_protein_obj <- normData$ruv_normalized_obj <- NULL
    if (!"correlation_filtered" %in% active_history) normData$correlation_filtered_obj <- NULL
    ruv <- metadata$ruv_optimization_result
    if (!is.null(ruv$best_cancor_plot) && !inherits(ruv$best_cancor_plot, "grob"))
      ruv$best_cancor_plot <- ggplot2::ggplotGrob(stripProtNormQcPlotOwners(ruv$best_cancor_plot))
    normData$ruv_optimization_result <- ruv
    normData$best_k <- metadata$ruv_k
    normData$control_genes_index <- ruv$best_control_genes_index
    normData$correlation_threshold <- metadata$correlation_threshold
    normData$pre_norm_qc_generated <- FALSE
    normData$plot_refresh_trigger <- normData$plot_refresh_trigger + 1L
    if (!is.null(metadata$normalization_method)) shiny::updateSelectInput(session, "norm_method",
      selected = metadata$normalization_method)
    invisible(NULL)
  })
}

#' The tab's three QC generators on its state (pre-normalisation, post-normalisation, RUV-corrected), as the tab's own
#' shell builds them; for the shared frame (KIT-120).
#' @noRd
protNormQcGenerators <- function(workflowData, normData, experimentPaths, getPlotAesthetics) {
  list(
    pre = function() {
      normData$qc_plot_paths <- generateProtNormPreNormalizationQcArtifacts(
        stateManager = workflowData$state_manager,
        qcDir = experimentPaths$protein_qc_dir,
        aesthetics = getPlotAesthetics(),
        qcPlotPaths = normData$qc_plot_paths
      )
    },
    post = function(normalizedS4) {
      normData$qc_plot_paths <- generateProtNormPostNormalizationQcArtifacts(
        normalizedS4 = normalizedS4,
        qcDir = normData$publication_paths$protein_qc_dir %||% experimentPaths$protein_qc_dir,
        aesthetics = getPlotAesthetics(),
        qcPlotPaths = normData$qc_plot_paths
      )
      normData$plot_refresh_trigger <- normData$plot_refresh_trigger + 1
    },
    ruv = function(ruvCorrectedS4) {
      normData$qc_plot_paths <- generateProtNormRuvCorrectedQcArtifacts(
        ruvCorrectedS4 = ruvCorrectedS4,
        qcDir = normData$publication_paths$protein_qc_dir %||% experimentPaths$protein_qc_dir,
        aesthetics = getPlotAesthetics(),
        qcPlotPaths = normData$qc_plot_paths
      )
      normData$plot_refresh_trigger <- normData$plot_refresh_trigger + 1
    }
  )
}
