prepareProtNormNormalizationRun <- function(
  stateManager,
  normData,
  reqFn = shiny::req,
  initialiseGridFn = InitialiseGrid,
  messageFn = message
) {
  reqFn(stateManager)

  current_state <- stateManager$current_state
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

  if (existsFn("config_list", envir = .GlobalEnv)) {
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

