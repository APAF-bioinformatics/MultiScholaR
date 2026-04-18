buildProtNormManualRuvResult <- function(
  percentageAsNegCtrl,
  ruvK,
  controlGenesIndex,
  cancorPlot
) {
  list(
    best_percentage = percentageAsNegCtrl,
    best_k = ruvK,
    best_control_genes_index = controlGenesIndex,
    best_cancor_plot = cancorPlot,
    optimization_results = data.frame(
      percentage = percentageAsNegCtrl,
      separation_score = NA,
      best_k = ruvK,
      composite_score = NA,
      num_controls = sum(controlGenesIndex, na.rm = TRUE),
      valid_plot = TRUE
    ),
    separation_metric_used = "manual",
    k_penalty_weight = NA,
    adaptive_k_penalty_used = FALSE
  )
}

updateProtNormRuvAuditTrail <- function(
  ruvK,
  controlGenesIndex,
  percentageAsNegCtrl,
  modeLabel,
  existsFn = exists,
  getFn = get,
  assignFn = assign,
  updateRuvParametersFn = updateRuvParameters,
  messageFn = message
) {
  tryCatch({
    if (existsFn("config_list", envir = .GlobalEnv)) {
      config_list <- getFn("config_list", envir = .GlobalEnv)
      config_list <- updateRuvParametersFn(
        config_list,
        ruvK,
        controlGenesIndex,
        percentageAsNegCtrl
      )
      assignFn("config_list", config_list, envir = .GlobalEnv)
      messageFn(sprintf("*** AUDIT: Logged %s RUV parameters to config_list ***", modeLabel))
    } else {
      messageFn("*** WARNING: config_list not found in global environment - audit trail not updated ***")
    }
  }, error = function(e) {
    messageFn(paste("*** WARNING: Could not update config_list audit trail:", e$message, "***"))
  })

  invisible(NULL)
}

persistProtNormRuvResult <- function(
  ruvResult,
  workflowData,
  sourceDir = NULL,
  resultLabel,
  saveRdsFn = saveRDS,
  catFn = cat
) {
  workflowData$ruv_optimization_result <- ruvResult
  catFn(sprintf("*** STEP 3a: Stored %s in workflow_data ***\n", resultLabel))

  if (!is.null(sourceDir)) {
    tryCatch({
      ruv_file <- file.path(sourceDir, "ruv_optimization_results.RDS")
      saveRdsFn(ruvResult, ruv_file)
      catFn(sprintf("*** STEP 3a: Saved %s to: %s ***\n", resultLabel, ruv_file))
    }, error = function(e) {
      catFn(sprintf("*** STEP 3a: Warning - could not save %s file: %s ***\n", resultLabel, e$message))
    })
  }

  invisible(ruvResult)
}

resolveProtNormRuvParameters <- function(
  normalizedS4,
  input,
  normData,
  workflowData,
  sourceDir = NULL,
  getRuvGroupingVariableFn,
  withProgressFn = shiny::withProgress,
  findBestNegCtrlPercentageFn = findBestNegCtrlPercentage,
  getNegCtrlProtAnovaFn = getNegCtrlProtAnova,
  ruvCancorFn = ruvCancor,
  buildManualRuvResultFn = buildProtNormManualRuvResult,
  updateAuditTrailFn = updateProtNormRuvAuditTrail,
  persistRuvResultFn = persistProtNormRuvResult,
  messageFn = message
) {
  messageFn("*** STEP 3: Starting RUV parameter determination ***")

  if (input$ruv_mode == "automatic") {
    messageFn("*** STEP 3a: Running automatic RUV optimization ***")

    percentage_min <- if (is.null(input$auto_percentage_min) || is.na(input$auto_percentage_min)) 1 else input$auto_percentage_min
    percentage_max <- if (is.null(input$auto_percentage_max) || is.na(input$auto_percentage_max)) 20 else input$auto_percentage_max

    if (percentage_min >= percentage_max) {
      stop("Minimum percentage must be less than maximum percentage")
    }

    percentage_range <- seq(percentage_min, percentage_max, by = 1)

    optimization_result <- withProgressFn(
      message = "Optimizing RUV parameters...",
      detail = "Testing percentage range...",
      value = 0.5,
      expr = {
        messageFn("=== AUTOMATIC RUV OPTIMIZATION DEBUG ===")
        messageFn(sprintf("*** DEBUG: percentage_range = %s ***", paste(percentage_range, collapse = ", ")))
        messageFn(sprintf("*** DEBUG: separation_metric = %s ***", input$separation_metric))
        messageFn(sprintf("*** DEBUG: k_penalty_weight = %.2f ***", input$k_penalty_weight))
        messageFn(sprintf("*** DEBUG: adaptive_k_penalty = %s ***", input$adaptive_k_penalty))
        messageFn(sprintf("*** DEBUG: ruv_grouping_variable = %s ***", getRuvGroupingVariableFn()))
        messageFn(sprintf("*** DEBUG: normalized_s4 class = %s ***", class(normalizedS4)))
        messageFn(sprintf("*** DEBUG: normalized_s4 dims = %d x %d ***", nrow(normalizedS4@protein_quant_table), ncol(normalizedS4@protein_quant_table)))

        optimization_result <- findBestNegCtrlPercentageFn(
          normalised_protein_matrix_obj = normalizedS4,
          percentage_range = percentage_range,
          separation_metric = input$separation_metric,
          k_penalty_weight = input$k_penalty_weight,
          adaptive_k_penalty = input$adaptive_k_penalty,
          ruv_grouping_variable = getRuvGroupingVariableFn(),
          verbose = TRUE
        )

        messageFn("*** DEBUG: Optimization completed ***")
        messageFn(sprintf("*** DEBUG: optimization_result class = %s ***", class(optimization_result)))
        messageFn(sprintf("*** DEBUG: best_percentage = %.1f ***", optimization_result$best_percentage))
        messageFn(sprintf("*** DEBUG: best_k = %d ***", optimization_result$best_k))
        messageFn(sprintf("*** DEBUG: best_composite_score = %.6f ***", optimization_result$best_composite_score))
        messageFn(sprintf("*** DEBUG: best_separation_score = %.6f ***", optimization_result$best_separation_score))
        if (!is.null(optimization_result$optimization_results)) {
          messageFn("*** DEBUG: First few optimization results: ***")
          print(head(optimization_result$optimization_results, 3))
        } else {
          messageFn("*** DEBUG: optimization_results is NULL! ***")
        }

        messageFn("*** DEBUG: Checking optimization_results structure ***")
        messageFn(sprintf("*** DEBUG: optimization_results class = %s ***", class(optimization_result$optimization_results)))
        if (is.data.frame(optimization_result$optimization_results)) {
          messageFn(sprintf(
            "*** DEBUG: optimization_results dims = %d x %d ***",
            nrow(optimization_result$optimization_results),
            ncol(optimization_result$optimization_results)
          ))
          messageFn(sprintf(
            "*** DEBUG: optimization_results columns = %s ***",
            paste(colnames(optimization_result$optimization_results), collapse = ", ")
          ))
          messageFn("*** DEBUG: FULL optimization_results table: ***")
          print(optimization_result$optimization_results)
        } else {
          messageFn("*** DEBUG: optimization_results is not a data.frame ***")
        }

        messageFn("*** DEBUG: Testing verbose logging override ***")
        messageFn("*** DEBUG: This should help us see if the function is running at all ***")

        optimization_result
      }
    )

    normData$ruv_optimization_result <- optimization_result
    persistRuvResultFn(
      ruvResult = optimization_result,
      workflowData = workflowData,
      sourceDir = sourceDir,
      resultLabel = "RUV optimization results"
    )

    percentage_as_neg_ctrl <- optimization_result$best_percentage
    ruv_k <- optimization_result$best_k
    control_genes_index <- optimization_result$best_control_genes_index

    messageFn(sprintf(
      "*** STEP 3a: Automatic optimization completed - Best %%: %.1f, Best k: %d ***",
      percentage_as_neg_ctrl,
      ruv_k
    ))

    updateAuditTrailFn(
      ruvK = ruv_k,
      controlGenesIndex = control_genes_index,
      percentageAsNegCtrl = percentage_as_neg_ctrl,
      modeLabel = "automatic"
    )
  } else {
    messageFn("*** STEP 3a: Using manual RUV parameters ***")
    percentage_as_neg_ctrl <- input$ruv_percentage
    ruv_k <- if (is.null(input$ruv_k) || is.na(input$ruv_k)) 3 else input$ruv_k

    control_genes_index <- getNegCtrlProtAnovaFn(
      normalizedS4,
      ruv_grouping_variable = getRuvGroupingVariableFn(),
      percentage_as_neg_ctrl = percentage_as_neg_ctrl,
      ruv_qval_cutoff = 0.05,
      ruv_fdr_method = "BH"
    )

    updateAuditTrailFn(
      ruvK = ruv_k,
      controlGenesIndex = control_genes_index,
      percentageAsNegCtrl = percentage_as_neg_ctrl,
      modeLabel = "manual"
    )

    tryCatch({
      cancor_plot <- ruvCancorFn(
        normalizedS4,
        ctrl = control_genes_index,
        num_components_to_impute = 2,
        ruv_grouping_variable = getRuvGroupingVariableFn()
      )

      manual_ruv_result <- buildManualRuvResultFn(
        percentageAsNegCtrl = percentage_as_neg_ctrl,
        ruvK = ruv_k,
        controlGenesIndex = control_genes_index,
        cancorPlot = cancor_plot
      )

      normData$ruv_optimization_result <- manual_ruv_result
      persistRuvResultFn(
        ruvResult = manual_ruv_result,
        workflowData = workflowData,
        sourceDir = sourceDir,
        resultLabel = "manual RUV results"
      )
    }, error = function(e) {
      messageFn(paste("Warning: Could not generate canonical correlation plot for manual mode:", e$message))
    })
  }

  normData$control_genes_index <- control_genes_index
  normData$best_k <- ruv_k
  messageFn("*** STEP 3: RUV parameter determination completed ***")

  list(
    percentageAsNegCtrl = percentage_as_neg_ctrl,
    ruvK = ruv_k,
    controlGenesIndex = control_genes_index
  )
}

applyProtNormRuvCorrectionStep <- function(
  normalizedS4,
  normData,
  getRuvGroupingVariableFn,
  ruvIII_C_VaryingFn = ruvIII_C_Varying,
  captureCheckpointFn = .capture_checkpoint,
  messageFn = message
) {
  messageFn("*** STEP 4: Applying RUV-III correction ***")

  ruv_k <- normData$best_k
  control_genes_index <- normData$control_genes_index

  ruv_corrected_s4 <- ruvIII_C_VaryingFn(
    normalizedS4,
    ruv_grouping_variable = getRuvGroupingVariableFn(),
    ruv_number_k = ruv_k,
    ctrl = control_genes_index
  )

  captureCheckpointFn(ruv_corrected_s4, "cp06", "ruv_corrected")
  normData$ruv_normalized_obj <- ruv_corrected_s4
  messageFn("*** STEP 4: RUV-III correction completed ***")

  ruv_corrected_s4
}

finalizeProtNormRuvCleanupStep <- function(
  ruvCorrectedS4,
  input,
  normData,
  workflowData,
  omicType,
  experimentLabel,
  removeRowsWithMissingValuesPercentFn = removeRowsWithMissingValuesPercent,
  countDistinctProteinsFn = function(s4Object) {
    s4Object@protein_quant_table |>
      dplyr::distinct(Protein.Ids) |>
      nrow()
  },
  updateProteinFilteringFn = updateProteinFiltering,
  messageFn = message
) {
  messageFn("*** STEP 5: Cleaning up missing values ***")

  ruv_corrected_s4_clean <- removeRowsWithMissingValuesPercentFn(
    theObject = ruvCorrectedS4
  )

  ruvfilt_protein_count <- countDistinctProteinsFn(ruv_corrected_s4_clean)
  messageFn(sprintf(
    "Number of distinct proteins remaining after RUV normalization and filtering: %d",
    ruvfilt_protein_count
  ))

  if (is.null(workflowData$protein_counts)) {
    workflowData$protein_counts <- list()
  }
  workflowData$protein_counts$after_ruv_filtering <- ruvfilt_protein_count
  messageFn(sprintf("*** STEP 5: Tracked protein count after RUV: %d ***", ruvfilt_protein_count))

  messageFn("*** STEP 5: About to call updateProteinFiltering ***")
  messageFn(sprintf(
    "*** STEP 5: omic_type = %s, experiment_label = %s ***",
    ifelse(is.null(omicType), "NULL", omicType),
    ifelse(is.null(experimentLabel), "NULL", experimentLabel)
  ))

  tryCatch({
    filtering_plot <- updateProteinFilteringFn(
      data = ruv_corrected_s4_clean@protein_quant_table,
      step_name = "11_RUV_filtered",
      omic_type = omicType,
      experiment_label = experimentLabel,
      return_grid = TRUE,
      overwrite = TRUE
    )

    normData$post_norm_filtering_plot <- filtering_plot
    messageFn("*** STEP 5: Protein filtering tracking updated successfully ***")
  }, error = function(e) {
    messageFn(paste("*** WARNING: updateProteinFiltering failed:", e$message, "***"))
    messageFn("*** STEP 5: Continuing without filtering plot update ***")
    normData$post_norm_filtering_plot <- NULL
  })

  best_percentage <- if (!is.null(normData$ruv_optimization_result)) {
    normData$ruv_optimization_result$best_percentage
  } else {
    input$ruv_percentage
  }

  normData$filtering_summary_text <- sprintf(
    "Filtering Summary through Normalization & RUV:\n\n* Pre-normalization proteins: [Previous step]\n* Post-RUV filtering: %d proteins\n* Proteins removed by RUV: [Calculated from difference]\n\nRUV Parameters:\n* Normalization method: %s\n* RUV mode: %s\n* RUV k value: %d\n* Negative control %%: %.1f",
    ruvfilt_protein_count,
    input$norm_method,
    input$ruv_mode,
    normData$best_k,
    best_percentage
  )

  normData$ruv_normalized_obj <- ruv_corrected_s4_clean
  messageFn("*** STEP 5: Missing values cleanup completed ***")
  messageFn("*** STEP 5: Saving state to R6 state manager ***")

  tryCatch({
    workflowData$state_manager$saveState(
      state_name = "ruv_corrected",
      s4_data_object = ruv_corrected_s4_clean,
      config_object = list(
        norm_method = input$norm_method,
        ruv_mode = input$ruv_mode,
        ruv_k = normData$best_k,
        percentage_as_neg_ctrl = best_percentage
      ),
      description = "Post-normalization complete: RUV-III correction and missing value cleanup completed"
    )
    messageFn("*** STEP 5: State saved successfully ***")
  }, error = function(e) {
    messageFn(paste("*** WARNING: Could not save state to R6 manager:", e$message, "***"))
    messageFn("*** STEP 5: Continuing without state save (Step 6 will still proceed) ***")
  })

  ruv_corrected_s4_clean
}

