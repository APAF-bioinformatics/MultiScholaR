prepareProtNormNormalizationRun <- function(
  stateManager,
  normData,
  reqFn = shiny::req,
  initialiseGridFn = InitialiseGrid,
  messageFn = message
) {
  reqFn(stateManager)

  current_state <- workflowStateCurrentName(stateManager)
  current_s4 <- stateManager$getState(current_state)

  if (is.null(current_s4)) {
    stop("No S4 object available for normalization")
  }

  messageFn("*** INITIALIZING S4 QC COMPOSITE FIGURE ***")
  normData$QC_composite_figure <- initialiseGridFn()

  list(currentState = current_state, currentS4 = current_s4)
}

runProtNormBetweenSamplesStep <- function(
  currentS4,
  normMethod,
  normData,
  workflowData = NULL,
  proteinQcDir = NULL,
  normaliseBetweenSamplesFn = normaliseBetweenSamples,
  captureCheckpointFn = .capture_checkpoint,
  existsFn = exists,
  getFn = get,
  assignFn = assign,
  saveMatrixFn = vroom::vroom_write,
  messageFn = message
) {
  messageFn("*** STEP 1: Starting between-samples normalization ***")

  if (protNormWorkflowDataAvailable(workflowData)) {
    config_list <- workflowData$config_list
    if (!is.list(config_list)) config_list <- list()
    if (!is.list(config_list$normaliseBetweenSamples)) {
      config_list$normaliseBetweenSamples <- list()
    }
    config_list$normaliseBetweenSamples$method <- normMethod
    workflowData$config_list <- config_list
  } else if (existsFn("config_list", envir = .GlobalEnv)) {
    config_list <- getFn("config_list", envir = .GlobalEnv)
    config_list$normaliseBetweenSamples$method <- normMethod
    assignFn("config_list", config_list, envir = .GlobalEnv)
  }

  normalized_s4 <- normaliseBetweenSamplesFn(currentS4, normalisation_method = normMethod)
  normData$normalized_protein_obj <- normalized_s4

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

persistProtNormNormalizedState <- function(
  currentS4,
  normalizedS4,
  input,
  workflowData,
  normData,
  experimentPaths,
  saveStateFn = saveProtNormState
) {
  if (!methods::is(normalizedS4, "ProteinQuantitativeData")) return(normalizedS4)
  if (!protDiaNormWorkflowIsDia(workflowData, workflowData$state_manager) &&
      is.null(protNonDiaNormDescriptor(workflowData))) {
    return(normalizedS4)
  }
  parameters <- list(
    normalization_method = input$norm_method,
    ruv_mode = input$ruv_mode,
    ruv_grouping_variable = workflowData$normalization_context$ruv_grouping_variable,
    skip_reason = if (identical(input$ruv_mode, "skip")) {
      buildProtNormSkippedRuvResult()$skip_reason
    } else {
      NULL
    },
    normalized_matrix_file = "normalized_protein_matrix_pre_ruv.tsv"
  )
  saved <- saveStateFn(
    workflow_data = workflowData,
    state_manager = workflowData$state_manager,
    before = currentS4,
    after = normalizedS4,
    stage_id = "normalization",
    state_name = "normalized",
    config_object = parameters,
    description = "Applied between-samples protein normalization",
    parameters = parameters,
    status = "applied",
    transformation_type = "normalization"
  )
  settleProtNormArtifactState(
    workflowData,
    normData,
    "normalization",
    "normalized",
    saved
  )
}

runProtNormPostNormalizationQcStep <- function(
  normalizedS4,
  normData,
  generatePostNormalizationQcFn,
  messageFn = message
) {
  messageFn("*** STEP 2: Generating post-normalization QC plots ***")
  generatePostNormalizationQcFn(normalizedS4)
  normData$normalization_complete <- TRUE
  messageFn("*** STEP 2: Post-normalization QC completed ***")

  invisible(NULL)
}

buildProtNormSkippedRuvResult <- function() {
  list(
    best_percentage = NA,
    best_k = NA,
    best_control_genes_index = NA,
    best_cancor_plot = NULL,
    ruv_skipped = TRUE,
    skip_reason = "User selected skip due to dataset constraints"
  )
}

applyProtNormSkippedRuvState <- function(
  normalizedS4,
  normMethod,
  normData,
  workflowData,
  sourceDir = NULL,
  skipResult = buildProtNormSkippedRuvResult(),
  saveRdsFn = saveRDS,
  messageFn = message
) {
  messageFn("*** RUV MODE: SKIP - Bypassing RUV-III correction ***")

  normData$ruv_normalized_obj <- normalizedS4
  normData$ruv_complete <- TRUE
  normData$best_k <- NA
  normData$control_genes_index <- NA
  normData$ruv_optimization_result <- skipResult

  workflowData$ruv_optimization_result <- skipResult
  updateProtNormWorkflowRuvContext(
    workflow_data = workflowData,
    mode = "skip",
    grouping_variable = workflowData$normalization_context$ruv_grouping_variable,
    percentage = NA_real_,
    k = NA_integer_,
    controls = NA,
    optimization_result = skipResult,
    input = list()
  )
  messageFn("*** RUV SKIP: Stored skip result in workflow_data for session summary ***")

  if (!is.null(sourceDir)) {
    tryCatch({
      ruv_file <- file.path(sourceDir, "ruv_optimization_results.RDS")
      saveRdsFn(skipResult, ruv_file)
      messageFn(sprintf("*** RUV SKIP: Saved skip result to file: %s (overwrites old results) ***", ruv_file))
    }, error = function(e) {
      messageFn(sprintf("*** RUV SKIP: Warning - could not save skip result file: %s ***", e$message))
    })
  }

  messageFn("*** RUV SKIP: Using normalized data directly for correlation filtering ***")
  messageFn("*** RUV SKIP: Saving state to R6 state manager ***")

  current_state <- tryCatch(
    workflowStateCurrentName(workflowData$state_manager),
    error = \(error) NULL
  )
  if (!identical(current_state, "normalized")) {
    tryCatch({
      workflowData$state_manager$saveState(
      state_name = "normalized",
      s4_data_object = normalizedS4,
      config_object = list(
        norm_method = normMethod,
        ruv_mode = "skip",
        ruv_applied = FALSE,
        ruv_k = NA,
        percentage_as_neg_ctrl = NA
      ),
      description = "Post-normalization complete: RUV-III skipped by user"
    )
      messageFn("*** RUV SKIP: State saved successfully ***")
    }, error = function(e) {
      messageFn(paste("*** WARNING: Could not save state to R6 manager:", e$message, "***"))
    })
  }

  messageFn("*** RUV SKIP: Proceeding to QC figure generation (2-column layout) ***")

  invisible(skipResult)
}

resolveProtNormStep6QcObject <- function(
  step5Object = NULL,
  normData,
  messageFn = message
) {
  if (!is.null(step5Object)) {
    messageFn("*** STEP 6: ruv_corrected_s4_clean is available ***")
    return(step5Object)
  }

  if (!is.null(normData$ruv_normalized_obj)) {
    messageFn("*** STEP 6: Using fallback ruv_normalized_obj from norm_data ***")
    return(normData$ruv_normalized_obj)
  }

  messageFn("*** STEP 6: WARNING - No RUV corrected object available ***")
  NULL
}

runProtNormStep6RuvQc <- function(
  ruvMode,
  step6Object,
  normData,
  qcDir,
  generateRuvCorrectedQcFn,
  ggsaveFn = ggplot2::ggsave,
  initPathsFn = initializeProtNormQcPlotPaths,
  messageFn = message
) {
  messageFn(sprintf("*** STEP 6: experiment_paths$protein_qc_dir: %s ***", qcDir))

  if (ruvMode == "skip") {
    messageFn("*** STEP 6A: Skipping RUV-corrected QC plots (RUV was skipped) ***")
    return(invisible(NULL))
  }

  tryCatch({
    generateRuvCorrectedQcFn(step6Object)
    messageFn("*** STEP 6A: RUV-corrected QC plots completed ***")
  }, error = function(e) {
    messageFn(paste("Warning: generateRuvCorrectedQc failed:", e$message))
    messageFn("*** STEP 6A: Continuing without RUV QC plots ***")
  })

  tryCatch({
    if (!is.null(normData$ruv_optimization_result$best_cancor_plot)) {
      if (!is.null(qcDir)) {
        cancor_path <- file.path(qcDir, "ruv_corrected_cancor.png")
        ggsaveFn(
          cancor_path,
          normData$ruv_optimization_result$best_cancor_plot,
          width = 8,
          height = 6,
          dpi = 150
        )
        normData$qc_plot_paths <- initPathsFn(normData$qc_plot_paths)
        normData$qc_plot_paths$ruv_corrected$cancor <- cancor_path
        messageFn("*** STEP 6A: Cancor plot saved to disk ***")
      }
    } else {
      messageFn("*** STEP 6A: No cancor plot available to save ***")
    }
  }, error = function(e) {
    messageFn(paste("Warning: Could not save cancor plot:", e$message))
  })

  invisible(NULL)
}

resolveProtNormCompositeFigureInputs <- function(ruvMode, qcDir) {
  if (ruvMode == "skip") {
    ncol_composite <- 2
    file_names <- c(
      "pre_norm_pca.png", "post_norm_pca.png",
      "pre_norm_density.png", "post_norm_density.png",
      "pre_norm_rle.png", "post_norm_rle.png",
      "pre_norm_correlation.png", "post_norm_correlation.png"
    )
    row_labels <- list(
      pca = c("a)", "b)"),
      density = c("c)", "d)"),
      rle = c("e)", "f)"),
      correlation = c("g)", "h)")
    )
    column_labels <- c("Pre-Normalisation", "Post-Normalisation")
  } else {
    ncol_composite <- 3
    file_names <- c(
      "pre_norm_pca.png", "post_norm_pca.png", "ruv_corrected_pca.png",
      "pre_norm_density.png", "post_norm_density.png", "ruv_corrected_density.png",
      "pre_norm_rle.png", "post_norm_rle.png", "ruv_corrected_rle.png",
      "pre_norm_correlation.png", "post_norm_correlation.png", "ruv_corrected_correlation.png",
      NA, NA, "ruv_corrected_cancor.png"
    )
    row_labels <- list(
      pca = c("a)", "b)", "c)"),
      density = c("d)", "e)", "f)"),
      rle = c("g)", "h)", "i)"),
      correlation = c("j)", "k)", "l)"),
      cancor = c("", "", "m)")
    )
    column_labels <- c("Pre-Normalisation", "Post-Normalisation", "RUV-Corrected")
  }

  plot_files <- if (!is.null(qcDir)) {
    vapply(
      file_names,
      function(fn) {
        if (is.na(fn)) NA_character_ else file.path(qcDir, fn)
      },
      character(1)
    )
  } else {
    NULL
  }

  list(
    ncol = ncol_composite,
    plotFiles = plot_files,
    rowLabels = row_labels,
    columnLabels = column_labels
  )
}

generateProtNormCompositeQcFigure <- function(
  ruvMode,
  qcDir,
  omicType,
  resolveInputsFn = resolveProtNormCompositeFigureInputs,
  generateCompositeFn = generateProtNormCompositeFromFiles,
  savePlotFn = savePlot,
  messageFn = message
) {
  messageFn("*** STEP 6B: Generating composite QC figure from saved images ***")

  if (is.null(qcDir)) {
    return(invisible(NULL))
  }

  composite_inputs <- resolveInputsFn(ruvMode, qcDir)

  tryCatch({
    composite_res <- generateCompositeFn(
      plotFiles = composite_inputs$plotFiles,
      ncol = composite_inputs$ncol,
      rowLabels = composite_inputs$rowLabels,
      columnLabels = composite_inputs$columnLabels
    )

    if (!is.null(composite_res)) {
      savePlotFn(
        composite_res$plot,
        base_path = qcDir,
        plot_name = paste0(omicType, "_composite_QC_figure"),
        width = composite_res$width,
        height = composite_res$height,
        dpi = 150,
        limitsize = FALSE
      )
      messageFn("*** STEP 6B: Composite QC figure saved to protein_qc_dir ***")
    }
  }, error = function(e) {
    messageFn(paste("Warning: Could not generate composite QC figure:", e$message))
  })

  invisible(NULL)
}

finalizeProtNormWorkflowState <- function(
  normData,
  gcFn = gc,
  messageFn = message
) {
  normData$QC_composite_figure <- NULL
  releaseProtNormArtifactStageObjects(normData)
  messageFn("*** STEP 6: Clearing redundant plot objects from memory ***")
  gcFn()

  normData$ruv_complete <- TRUE
  messageFn("*** STEP 6: RUV-corrected workflow completed ***")
  messageFn("*** STEP 7: Enabling correlation filtering step ***")
  messageFn("*** STEP 7: Normalization and RUV workflow completed - ready for correlation filtering ***")

  invisible(NULL)
}

buildProtNormCompletionNotification <- function(ruvMode) {
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
  messageFn = message
) {
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

  current_state <- workflowStateCurrentName(workflowData$state_manager)

  messageFn(sprintf("Current state: '%s'", current_state))
  messageFn(sprintf("Target trigger state: 'protein_replicate_filtered'"))
  messageFn(sprintf("Pre-norm QC generated: %s", normData$pre_norm_qc_generated))

  if (!shouldAutoGeneratePreQcFn(
    selectedTab,
    current_state,
    normData$pre_norm_qc_generated
  )) {
    return(notifyInvalidStateFn(currentState = current_state))
  }

  messageFn("*** AUTO-TRIGGERING PRE-NORMALIZATION QC (chunk 24) ***")

  tryCatch({
    generatePreNormalizationQcFn()
    normData$pre_norm_qc_generated <- TRUE

    messageFn("*** PRE-NORMALIZATION QC COMPLETED SUCCESSFULLY ***")
    messageFn(sprintf("State is '%s' and PCA already generated. skipping.", current_state))

    TRUE
  }, error = function(e) {
    handlePreQcErrorFn(error = e)
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

runProtNormNormalizationWorkflow <- function(
  input,
  workflowData,
  normData,
  experimentPaths,
  omicType,
  experimentLabel,
  checkMemoryUsageFn,
  generatePostNormalizationQcFn,
  generateRuvCorrectedQcFn,
  getRuvGroupingVariableFn,
  initializeWorkflowContextFn = initializeProtNormWorkflowContext,
  persistNormalizedStateFn = persistProtNormNormalizedState,
  prepareNormalizationRunFn = prepareProtNormNormalizationRun,
  runBetweenSamplesStepFn = runProtNormBetweenSamplesStep,
  runPostNormalizationQcStepFn = runProtNormPostNormalizationQcStep,
  applySkippedRuvStateFn = applyProtNormSkippedRuvState,
  resolveRuvParametersFn = resolveProtNormRuvParameters,
  applyRuvCorrectionStepFn = applyProtNormRuvCorrectionStep,
  finalizeRuvCleanupStepFn = finalizeProtNormRuvCleanupStep,
  resolveStep6QcObjectFn = resolveProtNormStep6QcObject,
  runStep6RuvQcFn = runProtNormStep6RuvQc,
  generateCompositeQcFigureFn = generateProtNormCompositeQcFigure,
  finalizeWorkflowStateFn = finalizeProtNormWorkflowState,
  buildCompletionNotificationFn = buildProtNormCompletionNotification,
  withProgressFn = shiny::withProgress,
  incProgressFn = shiny::incProgress,
  showNotificationFn = shiny::showNotification,
  handleErrorFn = handleProtNormNormalizationError,
  messageFn = message
) {
  messageFn("=== NORMALIZATION BUTTON CLICKED ===")
  messageFn("Starting normalization workflow...")

  checkMemoryUsageFn(threshold_gb = 8, context = "Normalization Start")

  tryCatch({
    initializeWorkflowContextFn(
      workflow_data = workflowData,
      input = input,
      experiment_paths = experimentPaths,
      omic_type = omicType,
      experiment_label = experimentLabel,
      grouping_variable = getRuvGroupingVariableFn()
    )
    run_context <- prepareNormalizationRunFn(
      stateManager = workflowData$state_manager,
      normData = normData
    )
    current_s4 <- run_context$currentS4
    ruv_corrected_s4_clean <- NULL

    withProgressFn(message = "Running normalization...", value = 0, {
      incProgressFn(0.2, detail = "Normalizing between samples...")

      normalized_s4 <- tryCatch({
        between_args <- list(
          currentS4 = current_s4,
          normMethod = input$norm_method,
          normData = normData,
          proteinQcDir = experimentPaths$protein_qc_dir
        )
        supported <- names(formals(runBetweenSamplesStepFn))
        if ("workflowData" %in% supported || "..." %in% supported) {
          between_args$workflowData <- workflowData
        }
        do.call(runBetweenSamplesStepFn, between_args)
      }, error = function(e) {
        stop(paste("Step 1 (normalization) error:", e$message))
      })

      normalized_s4 <- persistNormalizedStateFn(
        currentS4 = current_s4,
        normalizedS4 = normalized_s4,
        input = input,
        workflowData = workflowData,
        normData = normData,
        experimentPaths = experimentPaths
      )

      incProgressFn(0.2, detail = "Generating post-normalization QC plots...")

      tryCatch({
        runPostNormalizationQcStepFn(
          normalizedS4 = normalized_s4,
          normData = normData,
          generatePostNormalizationQcFn = generatePostNormalizationQcFn
        )
      }, error = function(e) {
        stop(paste("Step 2 (post-norm QC) error:", e$message))
      })

      messageFn(sprintf("*** DEBUG: input$ruv_mode value = '%s' ***", input$ruv_mode))
      messageFn(sprintf(
        "*** DEBUG: Checking if '%s' == 'skip': %s ***",
        input$ruv_mode,
        input$ruv_mode == "skip"
      ))

      if (input$ruv_mode == "skip") {
        applySkippedRuvStateFn(
          normalizedS4 = normalized_s4,
          normMethod = input$norm_method,
          normData = normData,
          workflowData = workflowData,
          sourceDir = experimentPaths$source_dir
        )
      } else {
        incProgressFn(0.2, detail = "Determining RUV parameters...")

        tryCatch({
          resolveRuvParametersFn(
            normalizedS4 = normalized_s4,
            input = input,
            normData = normData,
            workflowData = workflowData,
            sourceDir = experimentPaths$source_dir,
            getRuvGroupingVariableFn = getRuvGroupingVariableFn
          )
        }, error = function(e) {
          stop(paste("Step 3 (RUV parameter determination) error:", e$message))
        })

        incProgressFn(0.2, detail = "Applying RUV-III batch correction...")

        ruv_corrected_s4 <- tryCatch({
          applyRuvCorrectionStepFn(
            normalizedS4 = normalized_s4,
            normData = normData,
            getRuvGroupingVariableFn = getRuvGroupingVariableFn
          )
        }, error = function(e) {
          stop(paste("Step 4 (RUV-III correction) error:", e$message))
        })

        tryCatch({
          ruv_corrected_s4_clean <- finalizeRuvCleanupStepFn(
            ruvCorrectedS4 = ruv_corrected_s4,
            input = input,
            normData = normData,
            workflowData = workflowData,
            omicType = omicType,
            experimentLabel = experimentLabel
          )
        }, error = function(e) {
          messageFn(paste("*** ERROR in Step 5:", e$message, "***"))
          messageFn("*** STEP 5: Continuing to Step 6 despite error ***")
        })
      }

      incProgressFn(0.2, detail = "Generating RUV-corrected QC plots...")
      messageFn("*** STEP 6: STARTING RUV-corrected QC plot generation ***")

      step6_object <- resolveStep6QcObjectFn(
        step5Object = ruv_corrected_s4_clean,
        normData = normData
      )

      runStep6RuvQcFn(
        ruvMode = input$ruv_mode,
        step6Object = step6_object,
        normData = normData,
        qcDir = experimentPaths$protein_qc_dir,
        generateRuvCorrectedQcFn = generateRuvCorrectedQcFn
      )

      generateCompositeQcFigureFn(
        ruvMode = input$ruv_mode,
        qcDir = experimentPaths$protein_qc_dir,
        omicType = omicType
      )

      finalizeWorkflowStateFn(normData)
    })

    showNotificationFn(
      buildCompletionNotificationFn(input$ruv_mode),
      type = "message",
      duration = 10
    )

    messageFn("Normalization workflow completed successfully")

    TRUE
  }, error = function(e) {
    handleErrorFn(error = e)
  })
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
  messageFn = message
) {
  messageFn("=== CORRELATION FILTERING BUTTON CLICKED (DEBUG66 ACTIVE) ===")

  checkMemoryUsageFn(threshold_gb = 8, context = "Correlation Filtering Start")

  tryCatch({
    correlation_input <- resolveProtNormStateObject(
      workflow_data = workflowData,
      norm_data = normData,
      state_names = c("ruv_corrected", "normalized"),
      legacy_object = normData$ruv_normalized_obj,
      stage_id = "ruv_correction"
    )
    ruv_s4 <- resolveCorrelationInputObjectFn(
      ruvNormalizedObj = correlation_input,
      startMessage = "Starting Correlation Filter Flow",
      missingObjectMessage = "RUV correction must be completed before correlation filtering"
    )
    final_s4_for_de <- runApplyCorrelationWorkflowFn(
      ruvS4 = ruv_s4,
      correlationThreshold = input$min_pearson_correlation_threshold,
      normData = normData,
      experimentPaths = experimentPaths,
      workflowData = workflowData,
      omicType = omicType,
      experimentLabel = experimentLabel,
      getRuvGroupingVariableFn = getRuvGroupingVariableFn,
      getPlotAestheticsFn = getPlotAestheticsFn
    )

    completeCorrelationWorkflowFn(
      finalS4ForDe = final_s4_for_de,
      workflowData = workflowData,
      output = output,
      normData = normData,
      correlationThreshold = input$min_pearson_correlation_threshold,
      skipped = FALSE,
      successNotification = "Correlation filtering completed! Ready for differential expression analysis.",
      completionMessage = "=== CORRELATION FILTERING COMPLETED SUCCESSFULLY ===",
      messagePrefix = "*** CORRELATION"
    )

    final_s4_for_de
  }, error = function(e) {
    handleCorrelationErrorFn(
      error = e,
      logPrefix = "Error in correlation filtering:",
      notificationPrefix = "Error in correlation filtering:"
    )
  })
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
  messageFn = message
) {
  messageFn("=== SKIP CORRELATION FILTERING BUTTON CLICKED (DEBUG66 ACTIVE) ===")

  checkMemoryUsageFn(threshold_gb = 8, context = "Skip Correlation Filtering Start")

  tryCatch({
    correlation_input <- resolveProtNormStateObject(
      workflow_data = workflowData,
      norm_data = normData,
      state_names = c("ruv_corrected", "normalized"),
      legacy_object = normData$ruv_normalized_obj,
      stage_id = "ruv_correction"
    )
    ruv_s4 <- resolveCorrelationInputObjectFn(
      ruvNormalizedObj = correlation_input,
      startMessage = "Skipping Correlation Filter Flow",
      missingObjectMessage = "RUV correction must be completed before proceeding"
    )
    final_s4_for_de <- runSkipCorrelationWorkflowFn(
      ruvS4 = ruv_s4,
      normData = normData,
      experimentPaths = experimentPaths,
      workflowData = workflowData,
      omicType = omicType,
      experimentLabel = experimentLabel,
      getPlotAestheticsFn = getPlotAestheticsFn
    )

    completeCorrelationWorkflowFn(
      finalS4ForDe = final_s4_for_de,
      workflowData = workflowData,
      output = output,
      normData = normData,
      correlationThreshold = 0,
      skipped = TRUE,
      successNotification = "Correlation filtering skipped! Ready for differential expression analysis.",
      completionMessage = "=== SKIP CORRELATION FILTERING COMPLETED SUCCESSFULLY ===",
      messagePrefix = "*** SKIP CORRELATION"
    )

    final_s4_for_de
  }, error = function(e) {
    handleCorrelationErrorFn(
      error = e,
      logPrefix = "Error in skipping correlation filtering:",
      notificationPrefix = "Error in skipping correlation filtering:"
    )
  })
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
  messageFn("=== EXPORT FILTERED SESSION BUTTON CLICKED ===")

  export_object <- normData$correlation_filtered_obj
  if (is.null(export_object)) export_object <- workflowData$final_for_da_ref
  if (!canExportFilteredSessionFn(
    correlationFilteringComplete = normData$correlation_filtering_complete,
    correlationFilteredObj = export_object
  )) {
    return(notifyExportPrereqFn())
  }

  tryCatch({
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
}

runProtNormResetObserver <- function(
  workflowData,
  normData,
  output,
  ruvMode,
  groupingVariable,
  runResetWorkflowFn = runProtNormResetWorkflow,
  handleResetErrorFn = handleProtNormResetError,
  messageFn = message
) {
  messageFn("Resetting normalization...")

  tryCatch({
    runResetWorkflowFn(
      workflowData = workflowData,
      normData = normData,
      output = output,
      ruvMode = ruvMode,
      groupingVariable = groupingVariable
    )
  }, error = function(e) {
    handleResetErrorFn(e)
  })
}
