renderProtNormQcImage <- function(
  filename,
  altText,
  normData,
  proteinQcDir,
  renderImageFn = shiny::renderImage,
  buildQcImagePayloadFn = buildProtNormQcImagePayload
) {
  renderImageFn({
    normData$plot_refresh_trigger

    buildQcImagePayloadFn(proteinQcDir, filename, altText)
  }, deleteFile = FALSE)
}

registerProtNormQcImageOutputs <- function(
  output,
  normData,
  proteinQcDir,
  renderQcImageFn = renderProtNormQcImage
) {
  output$pca_post_filtering <- renderQcImageFn(
    filename = "pre_norm_pca.png",
    altText = "PCA Post-Filtering",
    normData = normData,
    proteinQcDir = proteinQcDir
  )
  output$density_post_filtering <- renderQcImageFn(
    filename = "pre_norm_density.png",
    altText = "Density Post-Filtering",
    normData = normData,
    proteinQcDir = proteinQcDir
  )
  output$rle_post_filtering <- renderQcImageFn(
    filename = "pre_norm_rle.png",
    altText = "RLE Post-Filtering",
    normData = normData,
    proteinQcDir = proteinQcDir
  )
  output$correlation_post_filtering <- renderQcImageFn(
    filename = "pre_norm_correlation.png",
    altText = "Correlation Post-Filtering",
    normData = normData,
    proteinQcDir = proteinQcDir
  )
  output$pca_post_normalization <- renderQcImageFn(
    filename = "post_norm_pca.png",
    altText = "PCA Post-Normalization",
    normData = normData,
    proteinQcDir = proteinQcDir
  )
  output$density_post_normalization <- renderQcImageFn(
    filename = "post_norm_density.png",
    altText = "Density Post-Normalization",
    normData = normData,
    proteinQcDir = proteinQcDir
  )
  output$rle_post_normalization <- renderQcImageFn(
    filename = "post_norm_rle.png",
    altText = "RLE Post-Normalization",
    normData = normData,
    proteinQcDir = proteinQcDir
  )
  output$correlation_post_normalization <- renderQcImageFn(
    filename = "post_norm_correlation.png",
    altText = "Correlation Post-Normalization",
    normData = normData,
    proteinQcDir = proteinQcDir
  )
  output$pca_ruv_corrected <- renderQcImageFn(
    filename = "ruv_corrected_pca.png",
    altText = "PCA RUV Corrected",
    normData = normData,
    proteinQcDir = proteinQcDir
  )
  output$density_ruv_corrected <- renderQcImageFn(
    filename = "ruv_corrected_density.png",
    altText = "Density RUV Corrected",
    normData = normData,
    proteinQcDir = proteinQcDir
  )
  output$rle_ruv_corrected <- renderQcImageFn(
    filename = "ruv_corrected_rle.png",
    altText = "RLE RUV Corrected",
    normData = normData,
    proteinQcDir = proteinQcDir
  )
  output$correlation_ruv_corrected <- renderQcImageFn(
    filename = "ruv_corrected_correlation.png",
    altText = "Correlation RUV Corrected",
    normData = normData,
    proteinQcDir = proteinQcDir
  )

  invisible(output)
}

renderProtNormPostNormFilteringSummary <- function(
  normData,
  renderPlotFn = shiny::renderPlot,
  gridDrawFn = grid::grid.draw,
  plotNewFn = graphics::plot.new,
  textFn = graphics::text
) {
  renderPlotFn({
    if (!is.null(normData$post_norm_filtering_plot)) {
      gridDrawFn(normData$post_norm_filtering_plot)
    } else {
      plotNewFn()
      textFn(
        0.5,
        0.5,
        "Complete normalization and RUV correction\nto generate filtering summary",
        cex = 1.2
      )
    }
  })
}

checkProtNormMemoryUsage <- function(
  threshold_gb = 8,
  context = "",
  gcFn = gc,
  warningFn = warning,
  messageFn = message
) {
  mem_info <- gcFn()
  mem_used_mb <- sum(mem_info[, 2])
  mem_used_gb <- mem_used_mb / 1024

  if (mem_used_gb > threshold_gb) {
    warningFn(sprintf(
      "*** HIGH MEMORY WARNING [%s]: %.1f GB used (threshold: %.1f GB) ***",
      context,
      mem_used_gb,
      threshold_gb
    ))
    messageFn(sprintf(
      "*** HIGH MEMORY WARNING [%s]: %.1f GB used ***",
      context,
      mem_used_gb
    ))
  } else {
    messageFn(sprintf(
      "*** MEMORY CHECK [%s]: %.1f GB used ***",
      context,
      mem_used_gb
    ))
  }

  invisible(mem_used_gb)
}

createProtNormReactiveState <- function(
  reactiveValuesFn = shiny::reactiveValues
) {
  reactiveValuesFn(
    pre_norm_qc_generated = FALSE,
    normalization_complete = FALSE,
    ruv_complete = FALSE,
    correlation_filtering_complete = FALSE,
    QC_composite_figure = NULL,
    qc_plots = list(
      post_filtering = list(
        pca = NULL,
        density = NULL,
        rle = NULL,
        correlation = NULL
      ),
      post_normalization = list(
        pca = NULL,
        density = NULL,
        rle = NULL,
        correlation = NULL
      ),
      ruv_corrected = list(
        pca = NULL,
        density = NULL,
        rle = NULL,
        correlation = NULL
      )
    ),
    normalized_protein_obj = NULL,
    ruv_normalized_obj = NULL,
    correlation_filtered_obj = NULL,
    best_k = NULL,
    control_genes_index = NULL,
    correlation_vector = NULL,
    correlation_threshold = NULL,
    final_qc_plot = NULL,
    final_filtering_plot = NULL,
    post_norm_filtering_plot = NULL,
    filtering_summary_text = NULL,
    plot_refresh_trigger = 0,
    ruv_optimization_result = NULL
  )
}

renderProtNormFilteringSummaryText <- function(
  normData,
  renderTextFn = shiny::renderText,
  getSummaryTextFn = getProtNormFilteringSummaryText
) {
  renderTextFn({
    getSummaryTextFn(normData$filtering_summary_text)
  })
}

renderProtNormFinalQcPlot <- function(
  normData,
  renderPlotFn = shiny::renderPlot,
  resolveFinalQcRenderStateFn = resolveProtNormFinalQcRenderState,
  gridNewpageFn = grid::grid.newpage,
  pushViewportFn = grid::pushViewport,
  popViewportFn = grid::popViewport,
  viewportFn = grid::viewport,
  gridLayoutFn = grid::grid.layout,
  gridDrawFn = grid::grid.draw,
  ggplotGrobFn = ggplotGrob,
  plotNewFn = graphics::plot.new,
  textFn = graphics::text,
  messageFn = message
) {
  renderPlotFn({
    messageFn("--- DEBUG66 [final_qc_plot]: Rendering ---")

    final_qc_state <- resolveFinalQcRenderStateFn(
      finalQcPlot = normData$final_qc_plot,
      finalFilteringPlot = normData$final_filtering_plot
    )

    if (final_qc_state$mode == "combined") {
      messageFn("   [final_qc_plot] Drawing PCA + Filtering Progression...")
      gridNewpageFn()
      pushViewportFn(
        viewportFn(layout = gridLayoutFn(2, 1, heights = c(0.4, 0.6)))
      )

      pushViewportFn(viewportFn(layout.pos.row = 1))
      gridDrawFn(ggplotGrobFn(final_qc_state$finalQcPlot))
      popViewportFn()

      pushViewportFn(viewportFn(layout.pos.row = 2))
      gridDrawFn(final_qc_state$finalFilteringPlot)
      popViewportFn()

      popViewportFn()
      messageFn("   [final_qc_plot] Render complete.")
    } else if (final_qc_state$mode == "filtering_only") {
      messageFn("   [final_qc_plot] Drawing Filtering Progression ONLY...")
      gridDrawFn(final_qc_state$finalFilteringPlot)
      messageFn("   [final_qc_plot] Render complete.")
    } else if (final_qc_state$mode == "pca_only") {
      messageFn("   [final_qc_plot] Drawing PCA ONLY...")
      final_qc_state$finalQcPlot
    } else {
      plotNewFn()
      textFn(0.5, 0.5, "Apply correlation filter to generate final QC plot", cex = 1.2)
    }
  })
}

renderProtNormRuvCanonicalCorrelationPlot <- function(
  normData,
  renderPlotFn = shiny::renderPlot,
  getCanonicalCorrelationPlotFn = getProtNormRuvCanonicalCorrelationPlot,
  plotNewFn = graphics::plot.new,
  textFn = graphics::text
) {
  renderPlotFn({
    cancor_plot <- getCanonicalCorrelationPlotFn(normData$ruv_optimization_result)

    if (!is.null(cancor_plot)) {
      cancor_plot
    } else {
      plotNewFn()
      textFn(0.5, 0.5, "Run normalization to generate RUV canonical correlation plot", cex = 1.2)
    }
  })
}

renderProtNormRuvOptimizationSummary <- function(
  normData,
  ruvMode,
  groupingVariable,
  renderTextFn = shiny::renderText,
  buildSummaryFn = buildProtNormRuvOptimizationSummary
) {
  renderTextFn({
    buildSummaryFn(
      ruvOptimizationResult = normData$ruv_optimization_result,
      ruvMode = ruvMode,
      groupingVariable = groupingVariable
    )
  })
}

renderProtNormRuvOptimizationTable <- function(
  normData,
  renderDataTableFn = DT::renderDataTable,
  prepareOptimizationResultsTableFn = prepareProtNormOptimizationResultsTable,
  datatableFn = DT::datatable,
  formatStyleFn = DT::formatStyle,
  styleEqualFn = DT::styleEqual
) {
  renderDataTableFn({
    optimization_table <- prepareOptimizationResultsTableFn(normData$ruv_optimization_result)

    datatable_output <- datatableFn(
      optimization_table$results,
      options = list(
        pageLength = 10,
        scrollY = "300px",
        scrollCollapse = TRUE,
        dom = "t"
      ),
      rownames = FALSE
    )

    if (!optimization_table$hasResults) {
      return(datatable_output)
    }

    formatStyleFn(
      datatable_output,
      "percentage",
      target = "row",
      backgroundColor = styleEqualFn(optimization_table$bestPercentage, "#e6f3ff")
    )
  })
}

registerProtNormRenderOutputs <- function(
  output,
  normData,
  ruvMode,
  groupingVariable,
  renderPostNormFilteringSummaryFn = renderProtNormPostNormFilteringSummary,
  renderFilteringSummaryTextFn = renderProtNormFilteringSummaryText,
  renderFinalQcPlotFn = renderProtNormFinalQcPlot,
  renderRuvCanonicalCorrelationPlotFn = renderProtNormRuvCanonicalCorrelationPlot,
  renderRuvOptimizationSummaryFn = renderProtNormRuvOptimizationSummary,
  renderRuvOptimizationTableFn = renderProtNormRuvOptimizationTable
) {
  output$post_norm_filtering_summary <- renderPostNormFilteringSummaryFn(
    normData = normData
  )
  output$filtering_summary_text <- renderFilteringSummaryTextFn(
    normData = normData
  )
  output$final_qc_plot <- renderFinalQcPlotFn(
    normData = normData
  )
  output$ruv_canonical_correlation_plot <- renderRuvCanonicalCorrelationPlotFn(
    normData = normData
  )
  output$ruv_optimization_summary <- renderRuvOptimizationSummaryFn(
    normData = normData,
    ruvMode = ruvMode,
    groupingVariable = groupingVariable
  )
  output$ruv_optimization_table <- renderRuvOptimizationTableFn(
    normData = normData
  )

  invisible(output)
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

  current_state <- workflowData$state_manager$current_state

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
    run_context <- prepareNormalizationRunFn(
      stateManager = workflowData$state_manager,
      normData = normData
    )
    current_s4 <- run_context$currentS4
    ruv_corrected_s4_clean <- NULL

    withProgressFn(message = "Running normalization...", value = 0, {
      incProgressFn(0.2, detail = "Normalizing between samples...")

      normalized_s4 <- tryCatch({
        runBetweenSamplesStepFn(
          currentS4 = current_s4,
          normMethod = input$norm_method,
          normData = normData,
          proteinQcDir = experimentPaths$protein_qc_dir
        )
      }, error = function(e) {
        stop(paste("Step 1 (normalization) error:", e$message))
      })

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
    ruv_s4 <- resolveCorrelationInputObjectFn(
      ruvNormalizedObj = normData$ruv_normalized_obj,
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
    ruv_s4 <- resolveCorrelationInputObjectFn(
      ruvNormalizedObj = normData$ruv_normalized_obj,
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

  if (!canExportFilteredSessionFn(
    correlationFilteringComplete = normData$correlation_filtering_complete,
    correlationFilteredObj = normData$correlation_filtered_obj
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
  observeEventFn = shiny::observeEvent
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
      checkMemoryUsageFn = checkMemoryUsageFn,
      generatePostNormalizationQcFn = generatePostNormalizationQcFn,
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
      checkMemoryUsageFn = checkMemoryUsageFn
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
      checkMemoryUsageFn = checkMemoryUsageFn
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

runProtNormModuleServerShell <- function(
  input,
  output,
  session,
  id,
  workflowData,
  experimentPaths,
  omicType,
  experimentLabel,
  selectedTab = NULL,
  messageFn = message,
  createReactiveStateFn = createProtNormReactiveState,
  getPlotAestheticsFn = getProtNormPlotAesthetics,
  getRuvGroupingVariableFn = getProtNormRuvGroupingVariable,
  generatePreNormalizationQcArtifactsFn = generateProtNormPreNormalizationQcArtifacts,
  generatePostNormalizationQcArtifactsFn = generateProtNormPostNormalizationQcArtifacts,
  generateRuvCorrectedQcArtifactsFn = generateProtNormRuvCorrectedQcArtifacts,
  registerServerObserversFn = registerProtNormServerObservers,
  registerQcImageOutputsFn = registerProtNormQcImageOutputs,
  registerRenderOutputsFn = registerProtNormRenderOutputs,
  checkMemoryUsageFn = checkProtNormMemoryUsage
) {
  messageFn("=== NORMALIZATION MODULE SERVER STARTED ===")
  messageFn(sprintf("Module ID: %s", id))
  messageFn(sprintf("workflow_data is NULL: %s", is.null(workflowData)))
  if (!is.null(workflowData$state_manager)) {
    messageFn(sprintf("Current state at module start: %s", workflowData$state_manager$current_state))
  }

  normData <- createReactiveStateFn()

  getPlotAesthetics <- function() {
    getPlotAestheticsFn(input$color_variable, input$shape_variable)
  }

  getRuvGroupingVariable <- function() {
    getRuvGroupingVariableFn(input$ruv_grouping_variable)
  }

  generatePreNormalizationQc <- function() {
    normData$qc_plot_paths <- generatePreNormalizationQcArtifactsFn(
      stateManager = workflowData$state_manager,
      qcDir = experimentPaths$protein_qc_dir,
      aesthetics = getPlotAesthetics(),
      qcPlotPaths = normData$qc_plot_paths
    )
  }

  generatePostNormalizationQc <- function(normalizedS4) {
    normData$qc_plot_paths <- generatePostNormalizationQcArtifactsFn(
      normalizedS4 = normalizedS4,
      qcDir = experimentPaths$protein_qc_dir,
      aesthetics = getPlotAesthetics(),
      qcPlotPaths = normData$qc_plot_paths
    )

    normData$plot_refresh_trigger <- normData$plot_refresh_trigger + 1
  }

  generateRuvCorrectedQc <- function(ruvCorrectedS4) {
    normData$qc_plot_paths <- generateRuvCorrectedQcArtifactsFn(
      ruvCorrectedS4 = ruvCorrectedS4,
      qcDir = experimentPaths$protein_qc_dir,
      aesthetics = getPlotAesthetics(),
      qcPlotPaths = normData$qc_plot_paths
    )

    normData$plot_refresh_trigger <- normData$plot_refresh_trigger + 1
  }

  registerServerObserversFn(
    input = input,
    output = output,
    session = session,
    selectedTab = selectedTab,
    workflowData = workflowData,
    normData = normData,
    experimentPaths = experimentPaths,
    omicType = omicType,
    experimentLabel = experimentLabel,
    generatePreNormalizationQcFn = generatePreNormalizationQc,
    generatePostNormalizationQcFn = generatePostNormalizationQc,
    generateRuvCorrectedQcFn = generateRuvCorrectedQc,
    getPlotAestheticsFn = getPlotAesthetics,
    getRuvGroupingVariableFn = getRuvGroupingVariable,
    checkMemoryUsageFn = checkMemoryUsageFn
  )

  registerQcImageOutputsFn(
    output = output,
    normData = normData,
    proteinQcDir = experimentPaths$protein_qc_dir
  )

  registerRenderOutputsFn(
    output = output,
    normData = normData,
    ruvMode = input$ruv_mode,
    groupingVariable = getRuvGroupingVariable()
  )

  invisible(normData)
}

runProtNormModuleServerEntryShell <- function(
  id,
  workflowData,
  experimentPaths,
  omicType,
  experimentLabel,
  selectedTab = NULL,
  moduleServerFn = shiny::moduleServer,
  runModuleServerShellFn = runProtNormModuleServerShell
) {
  moduleServerFn(id, function(input, output, session) {
    runModuleServerShellFn(
      input = input,
      output = output,
      session = session,
      id = id,
      workflowData = workflowData,
      experimentPaths = experimentPaths,
      omicType = omicType,
      experimentLabel = experimentLabel,
      selectedTab = selectedTab
    )
  })
}

runProtNormModuleServerPublicWrapper <- function(
  id,
  workflow_data,
  experiment_paths,
  omic_type,
  experiment_label,
  selected_tab = NULL,
  runModuleServerEntryShellFn = runProtNormModuleServerEntryShell
) {
  runModuleServerEntryShellFn(
    id = id,
    workflowData = workflow_data,
    experimentPaths = experiment_paths,
    omicType = omic_type,
    experimentLabel = experiment_label,
    selectedTab = selected_tab
  )
}
