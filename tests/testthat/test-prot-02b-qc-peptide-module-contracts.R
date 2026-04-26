# fidelity-coverage-compare: shared
library(testthat)
library(shiny)

skipIfMissingProtQcPeptideRefactorHelpers <- function() {
  package_namespace <- asNamespace("MultiScholaR")
  required_symbols <- c(
    "getProtQcPeptideModuleSpecs",
    "buildProtQcPeptideTab",
    "runProtQcPeptideSubmodule"
  )
  missing <- required_symbols[
    !vapply(
      required_symbols,
      exists,
      logical(1),
      where = package_namespace,
      mode = "function",
      inherits = FALSE
    )
  ]
  if (length(missing) > 0) {
    testthat::skip(
      sprintf(
        "Target-only prot QC peptide helper(s) not present: %s",
        paste(missing, collapse = ", ")
      )
    )
  }
}

skipIfMissingProtQcPeptideRefactorHelpers()

makeFunctionWithOverrides <- function(fun, replacements) {
  funOverride <- fun
  environment(funOverride) <- list2env(replacements, parent = environment(fun))
  funOverride
}

test_that("mod_prot_qc_peptide_ui renders each namespaced submodule tab when available", {
  observed_ids <- list()
  ui_fn_names <- vapply(getProtQcPeptideModuleSpecs(), `[[`, character(1), "uiFn")
  overrides <- list(
    exists = function(x, where = -1, inherits = TRUE) {
      if (x %in% ui_fn_names) {
        TRUE
      } else {
        base::exists(x, where = where, inherits = inherits)
      }
    },
    get = function(x, pos = -1, envir = as.environment(pos), mode = "any", inherits = TRUE) {
      if (x %in% ui_fn_names) {
        function(id) {
          observed_ids[[x]] <<- id
          shiny::tabPanel(
            title = paste("ready", x),
            shiny::div(`data-module-id` = id, x)
          )
        }
      } else {
        base::get(x, pos = pos, envir = envir, mode = mode, inherits = inherits)
      }
    }
  )
  ui_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_peptide_ui,
    c(
      overrides,
      list(
        getProtQcPeptideModuleSpecs = getProtQcPeptideModuleSpecs,
        buildProtQcPeptideTab = makeFunctionWithOverrides(buildProtQcPeptideTab, overrides)
      )
    )
  )

  html <- htmltools::renderTags(ui_under_test("peptide-qc"))$html

  expect_match(html, "id=\"peptide-qc-peptide_filter_tabs\"", fixed = TRUE)
  expect_false(grepl("Module not loaded", html, fixed = TRUE))
  expect_identical(observed_ids$mod_prot_qc_peptide_qvalue_ui, "peptide-qc-qvalue_filter")
  expect_identical(observed_ids$mod_prot_qc_peptide_rollup_ui, "peptide-qc-rollup")
  expect_identical(observed_ids$mod_prot_qc_peptide_intensity_ui, "peptide-qc-intensity_filter")
  expect_identical(observed_ids$mod_prot_qc_peptide_protein_ui, "peptide-qc-protein_peptide_filter")
  expect_identical(observed_ids$mod_prot_qc_peptide_sample_ui, "peptide-qc-sample_filter")
  expect_identical(observed_ids$mod_prot_qc_peptide_replicate_ui, "peptide-qc-replicate_filter")
  expect_identical(observed_ids$mod_prot_qc_peptide_impute_ui, "peptide-qc-imputation")
})

test_that("mod_prot_qc_peptide_ui keeps fallback labels when submodules are unavailable", {
  ui_fn_names <- vapply(getProtQcPeptideModuleSpecs(), `[[`, character(1), "uiFn")
  overrides <- list(
    exists = function(x, where = -1, inherits = TRUE) {
      if (x %in% ui_fn_names) {
        FALSE
      } else {
        base::exists(x, where = where, inherits = inherits)
      }
    }
  )
  ui_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_peptide_ui,
    c(
      overrides,
      list(
        getProtQcPeptideModuleSpecs = getProtQcPeptideModuleSpecs,
        buildProtQcPeptideTab = makeFunctionWithOverrides(buildProtQcPeptideTab, overrides)
      )
    )
  )

  html <- htmltools::renderTags(ui_under_test("peptide-qc"))$html
  fallback_count <- lengths(regmatches(html, gregexpr("Module not loaded", html, fixed = TRUE)))

  expect_identical(fallback_count, 7L)
  expect_match(html, "Q-Value Filter", fixed = TRUE)
  expect_match(html, "Precursor Rollup", fixed = TRUE)
  expect_match(html, "Intensity Filter", fixed = TRUE)
  expect_match(html, "Protein Peptides", fixed = TRUE)
  expect_match(html, "Sample Quality", fixed = TRUE)
  expect_match(html, "Replicate Filter", fixed = TRUE)
  expect_match(html, "Imputation", fixed = TRUE)
})

test_that("mod_prot_qc_peptide_server forwards DIA orchestration to available submodules in order", {
  call_log <- character()
  call_args <- list()
  server_fn_names <- vapply(getProtQcPeptideModuleSpecs(), `[[`, character(1), "serverFn")
  workflow_data <- shiny::reactiveValues(state_manager = list(workflow_type = "DIA"))
  overrides <- list(
    exists = function(x, where = -1, inherits = TRUE) {
      if (x %in% server_fn_names) {
        TRUE
      } else {
        base::exists(x, where = where, inherits = inherits)
      }
    },
    get = function(x, pos = -1, envir = as.environment(pos), mode = "any", inherits = TRUE) {
      if (x %in% server_fn_names) {
        function(id, workflowData, omicType, experimentLabel) {
          call_log <<- c(call_log, x)
          call_args[[x]] <<- list(
            id = id,
            workflowData = workflowData,
            omicType = omicType,
            experimentLabel = experimentLabel
          )
          invisible(NULL)
        }
      } else {
        base::get(x, pos = pos, envir = envir, mode = mode, inherits = inherits)
      }
    }
  )
  server_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_peptide_server,
    c(
      overrides,
      list(
        getProtQcPeptideModuleSpecs = getProtQcPeptideModuleSpecs,
        runProtQcPeptideSubmodule = makeFunctionWithOverrides(runProtQcPeptideSubmodule, overrides)
      )
    )
  )

  testServer(
    server_under_test,
    args = list(
      workflow_data = workflow_data,
      experiment_paths = list(results_dir = tempdir()),
      omic_type = "proteomics",
      experiment_label = "DIA Experiment"
    ),
    {
      expect_true(TRUE)
    }
  )

  expect_identical(call_log, unname(server_fn_names))
  expect_identical(call_args$mod_prot_qc_peptide_qvalue_server$id, "qvalue_filter")
  expect_identical(call_args$mod_prot_qc_peptide_rollup_server$id, "rollup")
  expect_identical(call_args$mod_prot_qc_peptide_intensity_server$id, "intensity_filter")
  expect_identical(call_args$mod_prot_qc_peptide_protein_server$id, "protein_peptide_filter")
  expect_identical(call_args$mod_prot_qc_peptide_sample_server$id, "sample_filter")
  expect_identical(call_args$mod_prot_qc_peptide_replicate_server$id, "replicate_filter")
  expect_identical(call_args$mod_prot_qc_peptide_impute_server$id, "imputation")
  expect_identical(call_args$mod_prot_qc_peptide_qvalue_server$workflowData, workflow_data)
  expect_identical(call_args$mod_prot_qc_peptide_qvalue_server$omicType, "proteomics")
  expect_identical(call_args$mod_prot_qc_peptide_qvalue_server$experimentLabel, "DIA Experiment")
})

test_that("mod_prot_qc_peptide_server skips peptide submodules for non-DIA workflows", {
  call_log <- character()
  server_fn_names <- vapply(getProtQcPeptideModuleSpecs(), `[[`, character(1), "serverFn")
  workflow_data <- shiny::reactiveValues(state_manager = list(workflow_type = "LFQ"))
  overrides <- list(
    exists = function(x, where = -1, inherits = TRUE) {
      if (x %in% server_fn_names) {
        TRUE
      } else {
        base::exists(x, where = where, inherits = inherits)
      }
    },
    get = function(x, pos = -1, envir = as.environment(pos), mode = "any", inherits = TRUE) {
      if (x %in% server_fn_names) {
        function(id, workflowData, omicType, experimentLabel) {
          call_log <<- c(call_log, x)
          invisible(NULL)
        }
      } else {
        base::get(x, pos = pos, envir = envir, mode = mode, inherits = inherits)
      }
    }
  )
  server_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_peptide_server,
    c(
      overrides,
      list(
        getProtQcPeptideModuleSpecs = getProtQcPeptideModuleSpecs,
        runProtQcPeptideSubmodule = makeFunctionWithOverrides(runProtQcPeptideSubmodule, overrides)
      )
    )
  )

  testServer(
    server_under_test,
    args = list(
      workflow_data = workflow_data,
      experiment_paths = list(results_dir = tempdir()),
      omic_type = "proteomics",
      experiment_label = "LFQ Experiment"
    ),
    {
      expect_true(TRUE)
    }
  )

  expect_identical(call_log, character())
})

test_that("buildPeptideIntensityThresholdPreview formats both preview strings from one temporary state update", {
  if (!methods::isClass("FakePeptideIntensityPreview")) {
    methods::setClass("FakePeptideIntensityPreview", slots = c(args = "list"))
  }

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  current_state <- list(marker = "current_state")
  workflow_data$state_manager$getState <- function() {
    current_state
  }

  captured <- new.env(parent = emptyenv())
  captured$args <- NULL

  preview <- buildPeptideIntensityThresholdPreview(
    workflowData = workflow_data,
    minRepsPerGroup = 2,
    minGroups = 3,
    updateMissingValueParametersFn = function(theObject,
                                              min_reps_per_group,
                                              min_groups,
                                              function_name,
                                              grouping_variable) {
      captured$args <- list(
        theObject = theObject,
        min_reps_per_group = min_reps_per_group,
        min_groups = min_groups,
        function_name = function_name,
        grouping_variable = grouping_variable
      )

      methods::new(
        "FakePeptideIntensityPreview",
        args = list(
          peptideIntensityFiltering = list(
            groupwise_percentage_cutoff = 12.345,
            max_groups_percentage_cutoff = 67.89
          )
        )
      )
    }
  )

  expect_identical(captured$args$theObject, current_state)
  expect_identical(captured$args$min_reps_per_group, 2)
  expect_identical(captured$args$min_groups, 3)
  expect_identical(captured$args$function_name, "peptideIntensityFiltering")
  expect_identical(captured$args$grouping_variable, "group")
  expect_identical(preview$groupwiseText, "Groupwise % cutoff: 12.345%")
  expect_identical(preview$maxGroupsText, "Max groups % cutoff: 67.890%")
})

test_that("runPeptideIntensityRevertStep reverts to the previous state in history", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$getHistory <- function() {
    c("raw_data_s4", "qvalue_filtered", "intensity_filtered")
  }

  reverted_to <- NULL
  workflow_data$state_manager$revertToState <- function(state) {
    reverted_to <<- state
    "reverted_object"
  }

  result <- runPeptideIntensityRevertStep(workflowData = workflow_data)

  expect_equal(result$previousState, "qvalue_filtered")
  expect_equal(result$revertedS4, "reverted_object")
  expect_equal(result$resultText, "Reverted to previous state: qvalue_filtered")
  expect_equal(reverted_to, "qvalue_filtered")
})

test_that("runPeptideIntensityRevertStep errors when there is no prior state", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$getHistory <- function() {
    "raw_data_s4"
  }

  expect_error(
    runPeptideIntensityRevertStep(workflowData = workflow_data),
    "No previous state to revert to.",
    fixed = TRUE
  )
})

test_that("runPeptideIntensityRevertObserver delegates the revert workflow and success notification", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data_ref <- workflow_data
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()

  completed <- runPeptideIntensityRevertObserver(
    workflowData = workflow_data,
    output = output,
    runRevertStepFn = function(workflowData) {
      captured$calls <- c(captured$calls, "run")
      expect_identical(workflowData, workflow_data_ref)
      list(
        previousState = "qvalue_filtered",
        resultText = "Reverted to previous state: qvalue_filtered"
      )
    },
    renderTextFn = function(value) {
      captured$calls <- c(captured$calls, "render")
      value
    },
    showNotificationFn = function(...) {
      captured$calls <- c(captured$calls, "show")
      captured$notifications[[length(captured$notifications) + 1]] <- list(...)
    },
    logErrorFn = function(...) {
      stop("error logging should not run on the happy path")
    }
  )

  expect_identical(captured$calls, c("run", "render", "show"))
  expect_identical(
    output$intensity_results,
    "Reverted to previous state: qvalue_filtered"
  )
  expect_identical(captured$notifications[[1]][[1]], "Reverted successfully")
  expect_identical(captured$notifications[[1]]$type, "message")
  expect_identical(completed$status, "success")
  expect_identical(
    completed$revertResult$resultText,
    "Reverted to previous state: qvalue_filtered"
  )
})

test_that("runPeptideIntensityRevertObserver reports revert errors without rendering results", {
  workflow_data <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()

  completed <- runPeptideIntensityRevertObserver(
    workflowData = workflow_data,
    output = output,
    runRevertStepFn = function(...) {
      captured$calls <- c(captured$calls, "run")
      stop("mock revert failure")
    },
    renderTextFn = function(...) {
      stop("rendering should not happen on the error path")
    },
    showNotificationFn = function(...) {
      captured$calls <- c(captured$calls, "show")
      captured$notifications[[length(captured$notifications) + 1]] <- list(...)
    },
    logErrorFn = function(message) {
      captured$calls <- c(captured$calls, "log_error")
      captured$error <- message
    }
  )

  expect_identical(captured$calls, c("run", "log_error", "show"))
  expect_identical(captured$error, "Error reverting: mock revert failure")
  expect_identical(
    captured$notifications[[1]][[1]],
    "Error reverting: mock revert failure"
  )
  expect_identical(captured$notifications[[1]]$type, "error")
  expect_false(exists("intensity_results", envir = output, inherits = FALSE))
  expect_identical(completed$status, "error")
  expect_identical(completed$errorMessage, "Error reverting: mock revert failure")
})

test_that("runPeptideIntensityApplyStep applies flexible-mode updates and saves the filtered state", {
  if (!methods::isClass("FakePeptideIntensityState")) {
    methods::setClass(
      "FakePeptideIntensityState",
      slots = c(args = "list", peptide_data = "data.frame")
    )
  }

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$getState <- function() {
    methods::new(
      "FakePeptideIntensityState",
      args = list(
        peptideIntensityFiltering = list(
          groupwise_percentage_cutoff = 1,
          max_groups_percentage_cutoff = 2,
          peptides_intensity_cutoff_percentile = 0.5
        )
      ),
      peptide_data = data.frame(Protein.Ids = c("P0"))
    )
  }

  captured <- new.env(parent = emptyenv())
  captured$configUpdates <- list()
  captured$saveState <- NULL
  workflow_data$state_manager$saveState <- function(...) {
    captured$saveState <- list(...)
  }

  filtered_s4 <- methods::new(
    "FakePeptideIntensityState",
    args = list(
      peptideIntensityFiltering = list(
        groupwise_percentage_cutoff = 12.345,
        max_groups_percentage_cutoff = 67.89,
        peptides_intensity_cutoff_percentile = 1.5
      )
    ),
    peptide_data = data.frame(Protein.Ids = c("P1", "P1", "P2"))
  )
  timestamp <- as.POSIXct("2026-04-16 12:34:56", tz = "UTC")

  result <- runPeptideIntensityApplyStep(
    workflowData = workflow_data,
    useStrictMode = FALSE,
    minRepsPerGroup = 2,
    minGroups = 3,
    intensityCutoffPercentile = 1.5,
    updateConfigParameterFn = function(theObject, function_name, parameter_name, new_value) {
      captured$configUpdates[[length(captured$configUpdates) + 1]] <<- list(
        function_name = function_name,
        parameter_name = parameter_name,
        new_value = new_value
      )
      updated_args <- theObject@args
      updated_args$peptideIntensityFiltering[[parameter_name]] <- new_value
      methods::new(
        "FakePeptideIntensityState",
        args = updated_args,
        peptide_data = theObject@peptide_data
      )
    },
    updateMissingValueParametersFn = function(theObject,
                                              min_reps_per_group,
                                              min_groups,
                                              function_name,
                                              grouping_variable) {
      captured$missingUpdate <<- list(
        theObject = theObject,
        min_reps_per_group = min_reps_per_group,
        min_groups = min_groups,
        function_name = function_name,
        grouping_variable = grouping_variable
      )
      updated_args <- theObject@args
      updated_args$peptideIntensityFiltering$groupwise_percentage_cutoff <- 12.345
      updated_args$peptideIntensityFiltering$max_groups_percentage_cutoff <- 67.89
      methods::new(
        "FakePeptideIntensityState",
        args = updated_args,
        peptide_data = theObject@peptide_data
      )
    },
    peptideIntensityFilteringFn = function(theObject) {
      captured$filterInput <<- theObject
      filtered_s4
    },
    logInfoFn = function(message) {
      captured$info <- message
    },
    nowFn = function() timestamp
  )

  expect_identical(captured$info, "Peptide Processing: Using FLEXIBLE MODE")
  expect_identical(captured$missingUpdate$min_reps_per_group, 2)
  expect_identical(captured$missingUpdate$min_groups, 3)
  expect_identical(captured$missingUpdate$function_name, "peptideIntensityFiltering")
  expect_identical(captured$missingUpdate$grouping_variable, "group")
  expect_length(captured$configUpdates, 1)
  expect_identical(
    captured$configUpdates[[1]]$parameter_name,
    "peptides_intensity_cutoff_percentile"
  )
  expect_equal(captured$configUpdates[[1]]$new_value, 1.5)
  expect_identical(captured$filterInput@args$peptideIntensityFiltering$groupwise_percentage_cutoff, 12.345)
  expect_identical(captured$filterInput@args$peptideIntensityFiltering$max_groups_percentage_cutoff, 67.89)
  expect_equal(captured$filterInput@args$peptideIntensityFiltering$peptides_intensity_cutoff_percentile, 1.5)
  expect_false(workflow_data$qc_params$peptide_qc$intensity_filter$strict_mode)
  expect_identical(workflow_data$qc_params$peptide_qc$intensity_filter$min_reps_per_group, 2)
  expect_identical(workflow_data$qc_params$peptide_qc$intensity_filter$min_groups, 3)
  expect_equal(workflow_data$qc_params$peptide_qc$intensity_filter$intensity_cutoff_percentile, 1.5)
  expect_identical(workflow_data$qc_params$peptide_qc$intensity_filter$timestamp, timestamp)
  expect_identical(captured$saveState$state_name, "intensity_filtered")
  expect_identical(captured$saveState$s4_data_object, filtered_s4)
  expect_false(captured$saveState$config_object$strict_mode)
  expect_equal(captured$saveState$config_object$intensity_cutoff_percentile, 1.5)
  expect_identical(captured$saveState$description, "Applied FLEXIBLE peptide intensity filter")
  expect_identical(result$filteredS4, filtered_s4)
  expect_match(result$resultText, "Mode: FLEXIBLE", fixed = TRUE)
  expect_match(result$resultText, "Proteins remaining: 2", fixed = TRUE)
  expect_match(result$resultText, "Intensity cutoff percentile: 1.5%", fixed = TRUE)
  expect_match(result$resultText, "Groupwise % cutoff: 12.345%", fixed = TRUE)
  expect_match(result$resultText, "Max groups % cutoff: 67.890%", fixed = TRUE)
})

test_that("runPeptideIntensityApplyStep applies strict-mode zero cutoffs before filtering", {
  if (!methods::isClass("FakePeptideIntensityState")) {
    methods::setClass(
      "FakePeptideIntensityState",
      slots = c(args = "list", peptide_data = "data.frame")
    )
  }

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$getState <- function() {
    methods::new(
      "FakePeptideIntensityState",
      args = list(
        peptideIntensityFiltering = list(
          groupwise_percentage_cutoff = 55,
          max_groups_percentage_cutoff = 66,
          peptides_intensity_cutoff_percentile = 0.5
        )
      ),
      peptide_data = data.frame(Protein.Ids = c("P0"))
    )
  }

  captured <- new.env(parent = emptyenv())
  captured$configUpdates <- list()
  captured$saveState <- NULL
  workflow_data$state_manager$saveState <- function(...) {
    captured$saveState <- list(...)
  }

  filtered_s4 <- methods::new(
    "FakePeptideIntensityState",
    args = list(
      peptideIntensityFiltering = list(
        groupwise_percentage_cutoff = 0,
        max_groups_percentage_cutoff = 0,
        peptides_intensity_cutoff_percentile = 2.5
      )
    ),
    peptide_data = data.frame(Protein.Ids = c("P1", "P2", "P3"))
  )
  timestamp <- as.POSIXct("2026-04-16 13:00:00", tz = "UTC")

  result <- runPeptideIntensityApplyStep(
    workflowData = workflow_data,
    useStrictMode = TRUE,
    minRepsPerGroup = 4,
    minGroups = 5,
    intensityCutoffPercentile = 2.5,
    updateConfigParameterFn = function(theObject, function_name, parameter_name, new_value) {
      captured$configUpdates[[length(captured$configUpdates) + 1]] <<- list(
        function_name = function_name,
        parameter_name = parameter_name,
        new_value = new_value
      )
      updated_args <- theObject@args
      updated_args$peptideIntensityFiltering[[parameter_name]] <- new_value
      methods::new(
        "FakePeptideIntensityState",
        args = updated_args,
        peptide_data = theObject@peptide_data
      )
    },
    updateMissingValueParametersFn = function(...) {
      stop("flexible missing-value updates should not run in strict mode")
    },
    peptideIntensityFilteringFn = function(theObject) {
      captured$filterInput <<- theObject
      filtered_s4
    },
    logInfoFn = function(message) {
      captured$info <- message
    },
    nowFn = function() timestamp
  )

  expect_identical(captured$info, "Peptide Processing: Using STRICT MODE")
  expect_identical(
    vapply(captured$configUpdates, `[[`, character(1), "parameter_name"),
    c(
      "groupwise_percentage_cutoff",
      "max_groups_percentage_cutoff",
      "peptides_intensity_cutoff_percentile"
    )
  )
  expect_equal(
    vapply(captured$configUpdates, `[[`, numeric(1), "new_value"),
    c(0, 0, 2.5)
  )
  expect_identical(captured$filterInput@args$peptideIntensityFiltering$groupwise_percentage_cutoff, 0)
  expect_identical(captured$filterInput@args$peptideIntensityFiltering$max_groups_percentage_cutoff, 0)
  expect_equal(captured$filterInput@args$peptideIntensityFiltering$peptides_intensity_cutoff_percentile, 2.5)
  expect_true(workflow_data$qc_params$peptide_qc$intensity_filter$strict_mode)
  expect_true(is.na(workflow_data$qc_params$peptide_qc$intensity_filter$min_reps_per_group))
  expect_true(is.na(workflow_data$qc_params$peptide_qc$intensity_filter$min_groups))
  expect_equal(workflow_data$qc_params$peptide_qc$intensity_filter$intensity_cutoff_percentile, 2.5)
  expect_identical(workflow_data$qc_params$peptide_qc$intensity_filter$timestamp, timestamp)
  expect_identical(captured$saveState$description, "Applied STRICT peptide intensity filter")
  expect_true(captured$saveState$config_object$strict_mode)
  expect_identical(result$filteredS4, filtered_s4)
  expect_match(result$resultText, "Mode: STRICT", fixed = TRUE)
  expect_match(result$resultText, "Proteins remaining: 3", fixed = TRUE)
  expect_match(result$resultText, "Groupwise % cutoff: 0.000%", fixed = TRUE)
  expect_match(result$resultText, "Max groups % cutoff: 0.000%", fixed = TRUE)
})

test_that("updatePeptideIntensityOutputs refreshes result text and plot grid", {
  if (!methods::isClass("FakePeptideIntensityState")) {
    methods::setClass(
      "FakePeptideIntensityState",
      slots = c(args = "list", peptide_data = "data.frame")
    )
  }

  intensity_state <- methods::new(
    "FakePeptideIntensityState",
    args = list(peptideIntensityFiltering = list()),
    peptide_data = data.frame(Protein.Ids = c("P1", "P2"))
  )
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$args <- NULL
  captured$plot <- NULL
  intensity_plot <- function(value) {
    captured$plot <- value
    invisible(value)
  }

  plot_grid <- updatePeptideIntensityOutputs(
    output = output,
    intensityPlot = intensity_plot,
    intensityResult = list(
      filteredS4 = intensity_state,
      resultText = "Intensity completed"
    ),
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    renderTextFn = function(value) value,
    updateProteinFilteringFn = function(data, step_name, omic_type, experiment_label, return_grid, overwrite) {
      captured$args <- list(
        data = data,
        step_name = step_name,
        omic_type = omic_type,
        experiment_label = experiment_label,
        return_grid = return_grid,
        overwrite = overwrite
      )
      "plot_grid"
    }
  )

  expect_identical(output$intensity_results, "Intensity completed")
  expect_identical(plot_grid, "plot_grid")
  expect_identical(captured$plot, "plot_grid")
  expect_identical(captured$args$data, intensity_state@peptide_data)
  expect_identical(captured$args$step_name, "4_intensity_filtered")
  expect_identical(captured$args$omic_type, "proteomics")
  expect_identical(captured$args$experiment_label, "DIA Experiment")
  expect_true(captured$args$return_grid)
  expect_true(captured$args$overwrite)
})

test_that("runPeptideIntensityApplyObserver delegates the apply workflow and notifications", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data_ref <- workflow_data
  output <- new.env(parent = emptyenv())
  output_ref <- output
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()
  captured$removed <- character()
  intensity_plot <- function(...) "plot_reactive"

  completed <- runPeptideIntensityApplyObserver(
    workflowData = workflow_data,
    useStrictMode = TRUE,
    minRepsPerGroup = 2,
    minGroups = 3,
    intensityCutoffPercentile = 1.5,
    output = output,
    intensityPlot = intensity_plot,
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    runApplyStepFn = function(workflowData, useStrictMode, minRepsPerGroup, minGroups, intensityCutoffPercentile) {
      captured$calls <- c(captured$calls, "run")
      expect_identical(workflowData, workflow_data_ref)
      expect_true(useStrictMode)
      expect_identical(minRepsPerGroup, 2)
      expect_identical(minGroups, 3)
      expect_equal(intensityCutoffPercentile, 1.5)
      list(
        filteredS4 = "filtered_s4",
        resultText = "Intensity completed"
      )
    },
    updateOutputsFn = function(output, intensityPlot, intensityResult, omicType, experimentLabel) {
      captured$calls <- c(captured$calls, "update")
      expect_identical(output, output_ref)
      expect_identical(intensityPlot, intensity_plot)
      expect_identical(intensityResult$filteredS4, "filtered_s4")
      expect_identical(intensityResult$resultText, "Intensity completed")
      expect_identical(omicType, "proteomics")
      expect_identical(experimentLabel, "DIA Experiment")
      "plot_grid"
    },
    showNotificationFn = function(...) {
      captured$calls <- c(captured$calls, "show")
      captured$notifications[[length(captured$notifications) + 1]] <- list(...)
    },
    removeNotificationFn = function(id) {
      captured$calls <- c(captured$calls, "remove")
      captured$removed <- c(captured$removed, id)
    },
    logErrorFn = function(...) {
      stop("error logging should not run on the happy path")
    }
  )

  expect_identical(
    captured$calls,
    c("show", "run", "update", "remove", "show")
  )
  expect_identical(captured$notifications[[1]][[1]], "Applying intensity filter...")
  expect_identical(captured$notifications[[1]]$id, "intensity_working")
  expect_null(captured$notifications[[1]]$duration)
  expect_identical(
    captured$notifications[[2]][[1]],
    "Intensity filter applied successfully"
  )
  expect_identical(captured$notifications[[2]]$type, "message")
  expect_identical(captured$removed, "intensity_working")
  expect_identical(completed$status, "success")
  expect_identical(completed$intensityResult$filteredS4, "filtered_s4")
  expect_identical(completed$plotGrid, "plot_grid")
})

test_that("runPeptideIntensityApplyObserver reports apply errors and clears the working notification", {
  workflow_data <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()
  captured$removed <- character()
  intensity_plot <- function(...) "plot_reactive"

  completed <- runPeptideIntensityApplyObserver(
    workflowData = workflow_data,
    useStrictMode = FALSE,
    minRepsPerGroup = 2,
    minGroups = 3,
    intensityCutoffPercentile = 1.5,
    output = output,
    intensityPlot = intensity_plot,
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    runApplyStepFn = function(...) {
      captured$calls <- c(captured$calls, "run")
      stop("mock failure")
    },
    updateOutputsFn = function(...) {
      stop("output refresh should not run after an intensity error")
    },
    showNotificationFn = function(...) {
      captured$calls <- c(captured$calls, "show")
      captured$notifications[[length(captured$notifications) + 1]] <- list(...)
    },
    removeNotificationFn = function(id) {
      captured$calls <- c(captured$calls, "remove")
      captured$removed <- c(captured$removed, id)
    },
    logErrorFn = function(message) {
      captured$calls <- c(captured$calls, "log_error")
      captured$error <- message
    }
  )

  expect_identical(
    captured$calls,
    c("show", "run", "log_error", "show", "remove")
  )
  expect_identical(captured$notifications[[1]][[1]], "Applying intensity filter...")
  expect_identical(captured$notifications[[1]]$id, "intensity_working")
  expect_identical(
    captured$notifications[[2]][[1]],
    "Error applying intensity filter: mock failure"
  )
  expect_identical(captured$notifications[[2]]$type, "error")
  expect_identical(captured$notifications[[2]]$duration, 15)
  expect_identical(captured$removed, "intensity_working")
  expect_identical(captured$error, "Error applying intensity filter: mock failure")
  expect_identical(completed$status, "error")
  expect_identical(completed$errorMessage, "Error applying intensity filter: mock failure")
})

test_that("runPeptideQvalueApplyStep applies thresholds, records diagnostics, and saves the filtered state", {
  if (!methods::isClass("FakePeptideQvalueState")) {
    methods::setClass(
      "FakePeptideQvalueState",
      slots = c(args = "list", peptide_data = "data.frame")
    )
  }

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$getState <- function() {
    methods::new(
      "FakePeptideQvalueState",
      args = list(
        srlQvalueProteotypicPeptideClean = list(
          input_matrix_column_ids = c(" sample1", "", "sample3 "),
          qvalue_threshold = 0.05,
          global_qvalue_threshold = 0.05,
          choose_only_proteotypic_peptide = 0
        )
      ),
      peptide_data = data.frame(
        Protein.Ids = c("P0", "P0"),
        Intensity = c(1, 2)
      )
    )
  }

  captured <- new.env(parent = emptyenv())
  captured$configUpdates <- list()
  captured$info <- character()
  captured$warn <- character()
  captured$saveState <- NULL
  workflow_data$state_manager$saveState <- function(...) {
    captured$saveState <- list(...)
  }

  filtered_s4 <- methods::new(
    "FakePeptideQvalueState",
    args = list(
      srlQvalueProteotypicPeptideClean = list(
        qvalue_threshold = 0.01,
        global_qvalue_threshold = 0.02,
        choose_only_proteotypic_peptide = 1
      )
    ),
    peptide_data = data.frame(
      Protein.Ids = c("P1", "P1", "P2"),
      Intensity = c(10, 20, 30)
    )
  )
  timestamp <- as.POSIXct("2026-04-16 16:30:00", tz = "UTC")

  result <- runPeptideQvalueApplyStep(
    workflowData = workflow_data,
    qvalueThreshold = 0.01,
    globalQvalueThreshold = 0.02,
    proteotypicOnly = TRUE,
    updateConfigParameterFn = function(theObject, function_name, parameter_name, new_value) {
      captured$configUpdates[[length(captured$configUpdates) + 1]] <<- list(
        function_name = function_name,
        parameter_name = parameter_name,
        new_value = new_value
      )
      updated_args <- theObject@args
      updated_args$srlQvalueProteotypicPeptideClean[[parameter_name]] <- new_value
      methods::new(
        "FakePeptideQvalueState",
        args = updated_args,
        peptide_data = theObject@peptide_data
      )
    },
    qvalueFilterFn = function(theObject) {
      captured$filterInput <<- theObject
      filtered_s4
    },
    logInfoFn = function(message) {
      captured$info <- c(captured$info, message)
    },
    logWarnFn = function(message) {
      captured$warn <- c(captured$warn, message)
    },
    nowFn = function() timestamp
  )

  expect_true(any(grepl("S4 class = FakePeptideQvalueState", captured$info, fixed = TRUE)))
  expect_true(any(grepl("input_matrix_column_ids length = 3", captured$info, fixed = TRUE)))
  expect_true(any(grepl("thresholds 0.01, 0.02", captured$info, fixed = TRUE)))
  expect_identical(
    captured$warn,
    c(
      "Q-value filter: input_matrix_column_ids contains values with leading/trailing whitespace!",
      "Q-value filter: input_matrix_column_ids contains empty strings!"
    )
  )
  expect_length(captured$configUpdates, 3)
  expect_identical(captured$configUpdates[[1]]$parameter_name, "qvalue_threshold")
  expect_identical(captured$configUpdates[[1]]$new_value, 0.01)
  expect_identical(captured$configUpdates[[2]]$parameter_name, "global_qvalue_threshold")
  expect_identical(captured$configUpdates[[2]]$new_value, 0.02)
  expect_identical(captured$configUpdates[[3]]$parameter_name, "choose_only_proteotypic_peptide")
  expect_identical(captured$configUpdates[[3]]$new_value, 1)
  expect_identical(
    captured$filterInput@args$srlQvalueProteotypicPeptideClean$qvalue_threshold,
    0.01
  )
  expect_identical(
    captured$filterInput@args$srlQvalueProteotypicPeptideClean$global_qvalue_threshold,
    0.02
  )
  expect_identical(
    captured$filterInput@args$srlQvalueProteotypicPeptideClean$choose_only_proteotypic_peptide,
    1
  )
  expect_identical(workflow_data$qc_params$peptide_qc$qvalue_filter$qvalue_threshold, 0.01)
  expect_identical(workflow_data$qc_params$peptide_qc$qvalue_filter$global_qvalue_threshold, 0.02)
  expect_identical(workflow_data$qc_params$peptide_qc$qvalue_filter$proteotypic_only, TRUE)
  expect_identical(workflow_data$qc_params$peptide_qc$qvalue_filter$timestamp, timestamp)
  expect_identical(captured$saveState$state_name, "qvalue_filtered")
  expect_identical(captured$saveState$s4_data_object, filtered_s4)
  expect_identical(captured$saveState$config_object$qvalue_threshold, 0.01)
  expect_identical(captured$saveState$config_object$global_qvalue_threshold, 0.02)
  expect_identical(captured$saveState$config_object$proteotypic_only, TRUE)
  expect_identical(
    captured$saveState$description,
    "Applied Q-value and proteotypic peptide filter"
  )
  expect_identical(result$filteredS4, filtered_s4)
  expect_match(result$resultText, "Proteins remaining: 2", fixed = TRUE)
  expect_match(result$resultText, "Q-value threshold: 0.01", fixed = TRUE)
  expect_match(result$resultText, "Global Q-value threshold: 0.02", fixed = TRUE)
  expect_match(result$resultText, "Proteotypic only: TRUE", fixed = TRUE)
})

test_that("updatePeptideQvalueOutputs refreshes result text and plot grid", {
  if (!methods::isClass("FakePeptideQvalueState")) {
    methods::setClass(
      "FakePeptideQvalueState",
      slots = c(args = "list", peptide_data = "data.frame")
    )
  }

  qvalue_state <- methods::new(
    "FakePeptideQvalueState",
    args = list(srlQvalueProteotypicPeptideClean = list()),
    peptide_data = data.frame(Protein.Ids = c("P1", "P2"))
  )
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$args <- NULL
  captured$plot <- NULL
  qvalue_plot <- function(value) {
    captured$plot <- value
    invisible(value)
  }

  plot_grid <- updatePeptideQvalueOutputs(
    output = output,
    qvaluePlot = qvalue_plot,
    qvalueResult = list(
      filteredS4 = qvalue_state,
      resultText = "Q-value completed"
    ),
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    renderTextFn = function(value) value,
    updateProteinFilteringFn = function(data, step_name, omic_type, experiment_label, return_grid, overwrite) {
      captured$args <- list(
        data = data,
        step_name = step_name,
        omic_type = omic_type,
        experiment_label = experiment_label,
        return_grid = return_grid,
        overwrite = overwrite
      )
      "plot_grid"
    }
  )

  expect_identical(output$qvalue_results, "Q-value completed")
  expect_identical(plot_grid, "plot_grid")
  expect_identical(captured$plot, "plot_grid")
  expect_identical(captured$args$data, qvalue_state@peptide_data)
  expect_identical(captured$args$step_name, "2_qval_filtered")
  expect_identical(captured$args$omic_type, "proteomics")
  expect_identical(captured$args$experiment_label, "DIA Experiment")
  expect_true(captured$args$return_grid)
  expect_true(captured$args$overwrite)
})

test_that("runPeptideQvalueApplyObserver delegates the apply workflow and notifications", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data_ref <- workflow_data
  output <- new.env(parent = emptyenv())
  output_ref <- output
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()
  captured$removed <- character()
  qvalue_plot <- function(...) "plot_reactive"

  completed <- runPeptideQvalueApplyObserver(
    workflowData = workflow_data,
    qvalueThreshold = 0.01,
    globalQvalueThreshold = 0.02,
    proteotypicOnly = TRUE,
    output = output,
    qvaluePlot = qvalue_plot,
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    runApplyStepFn = function(workflowData, qvalueThreshold, globalQvalueThreshold, proteotypicOnly) {
      captured$calls <- c(captured$calls, "run")
      expect_identical(workflowData, workflow_data_ref)
      expect_identical(qvalueThreshold, 0.01)
      expect_identical(globalQvalueThreshold, 0.02)
      expect_identical(proteotypicOnly, TRUE)
      list(
        filteredS4 = "filtered_s4",
        resultText = "Q-value completed"
      )
    },
    updateOutputsFn = function(output, qvaluePlot, qvalueResult, omicType, experimentLabel) {
      captured$calls <- c(captured$calls, "update")
      expect_identical(output, output_ref)
      expect_identical(qvaluePlot, qvalue_plot)
      expect_identical(qvalueResult$filteredS4, "filtered_s4")
      expect_identical(qvalueResult$resultText, "Q-value completed")
      expect_identical(omicType, "proteomics")
      expect_identical(experimentLabel, "DIA Experiment")
      "plot_grid"
    },
    showNotificationFn = function(...) {
      captured$calls <- c(captured$calls, "show")
      captured$notifications[[length(captured$notifications) + 1]] <- list(...)
    },
    removeNotificationFn = function(id) {
      captured$calls <- c(captured$calls, "remove")
      captured$removed <- c(captured$removed, id)
    },
    logInfoFn = function(message) {
      captured$calls <- c(captured$calls, "log_info")
      captured$info <- message
    },
    logErrorFn = function(...) {
      stop("error logging should not run on the happy path")
    }
  )

  expect_identical(
    captured$calls,
    c("show", "run", "update", "log_info", "remove", "show")
  )
  expect_identical(captured$notifications[[1]][[1]], "Applying Q-value filter...")
  expect_identical(captured$notifications[[1]]$id, "qvalue_working")
  expect_null(captured$notifications[[1]]$duration)
  expect_identical(
    captured$notifications[[2]][[1]],
    "Q-value filter applied successfully"
  )
  expect_identical(captured$notifications[[2]]$type, "message")
  expect_identical(captured$removed, "qvalue_working")
  expect_identical(captured$info, "Q-value filter applied successfully")
  expect_identical(completed$status, "success")
  expect_identical(completed$qvalueResult$filteredS4, "filtered_s4")
  expect_identical(completed$plotGrid, "plot_grid")
})

test_that("runPeptideQvalueApplyObserver reports apply errors and clears the working notification", {
  workflow_data <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()
  captured$removed <- character()
  qvalue_plot <- function(...) "plot_reactive"

  completed <- runPeptideQvalueApplyObserver(
    workflowData = workflow_data,
    qvalueThreshold = 0.01,
    globalQvalueThreshold = 0.02,
    proteotypicOnly = TRUE,
    output = output,
    qvaluePlot = qvalue_plot,
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    runApplyStepFn = function(...) {
      captured$calls <- c(captured$calls, "run")
      stop("mock failure")
    },
    updateOutputsFn = function(...) {
      stop("output refresh should not run after a qvalue error")
    },
    showNotificationFn = function(...) {
      captured$calls <- c(captured$calls, "show")
      captured$notifications[[length(captured$notifications) + 1]] <- list(...)
    },
    removeNotificationFn = function(id) {
      captured$calls <- c(captured$calls, "remove")
      captured$removed <- c(captured$removed, id)
    },
    logInfoFn = function(...) {
      stop("success logging should not run on the error path")
    },
    logErrorFn = function(message) {
      captured$calls <- c(captured$calls, "log_error")
      captured$error <- message
    }
  )

  expect_identical(
    captured$calls,
    c("show", "run", "log_error", "show", "remove")
  )
  expect_identical(captured$notifications[[1]][[1]], "Applying Q-value filter...")
  expect_identical(captured$notifications[[1]]$id, "qvalue_working")
  expect_identical(
    captured$notifications[[2]][[1]],
    "Error applying Q-value filter: mock failure"
  )
  expect_identical(captured$notifications[[2]]$type, "error")
  expect_identical(captured$notifications[[2]]$duration, 15)
  expect_identical(captured$removed, "qvalue_working")
  expect_identical(captured$error, "Error applying Q-value filter: mock failure")
  expect_identical(completed$status, "error")
  expect_identical(
    completed$errorMessage,
    "Error applying Q-value filter: mock failure"
  )
})

test_that("runPeptideQvalueRevertStep reverts to the raw data state when present in history", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$getHistory <- function() {
    c("raw_data_s4", "qvalue_filtered")
  }

  captured <- new.env(parent = emptyenv())
  captured$revertedTo <- NULL
  captured$info <- NULL
  workflow_data$state_manager$revertToState <- function(state) {
    captured$revertedTo <- state
    "reverted_object"
  }

  result <- runPeptideQvalueRevertStep(
    workflowData = workflow_data,
    logInfoFn = function(message) {
      captured$info <- message
    }
  )

  expect_identical(captured$revertedTo, "raw_data_s4")
  expect_identical(captured$info, "Reverted to raw data state")
  expect_identical(result$revertedS4, "reverted_object")
  expect_identical(result$revertStateName, "raw_data_s4")
  expect_identical(result$resultText, "Reverted to raw data state")
})

test_that("runPeptideQvalueRevertStep errors when the raw data state is missing", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$getHistory <- function() {
    c("qvalue_filtered", "protein_peptide_filtered")
  }

  expect_error(
    runPeptideQvalueRevertStep(workflowData = workflow_data),
    "Cannot revert: 'raw_data_s4' state not found in history.",
    fixed = TRUE
  )
})

test_that("runPeptideQvalueRevertObserver delegates the revert workflow and success notification", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data_ref <- workflow_data
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()

  completed <- runPeptideQvalueRevertObserver(
    workflowData = workflow_data,
    output = output,
    runRevertStepFn = function(workflowData) {
      captured$calls <- c(captured$calls, "run")
      expect_identical(workflowData, workflow_data_ref)
      list(
        revertStateName = "raw_data_s4",
        resultText = "Reverted to raw data state"
      )
    },
    renderTextFn = function(value) {
      captured$calls <- c(captured$calls, "render")
      value
    },
    showNotificationFn = function(...) {
      captured$calls <- c(captured$calls, "show")
      captured$notifications[[length(captured$notifications) + 1]] <- list(...)
    },
    logErrorFn = function(...) {
      stop("error logging should not run on the happy path")
    }
  )

  expect_identical(captured$calls, c("run", "render", "show"))
  expect_identical(output$qvalue_results, "Reverted to raw data state")
  expect_identical(captured$notifications[[1]][[1]], "Reverted to raw data state")
  expect_identical(captured$notifications[[1]]$type, "message")
  expect_identical(completed$status, "success")
  expect_identical(
    completed$revertResult$resultText,
    "Reverted to raw data state"
  )
})

test_that("runPeptideQvalueRevertObserver reports revert errors without rendering results", {
  workflow_data <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()

  completed <- runPeptideQvalueRevertObserver(
    workflowData = workflow_data,
    output = output,
    runRevertStepFn = function(...) {
      captured$calls <- c(captured$calls, "run")
      stop("mock revert failure")
    },
    renderTextFn = function(...) {
      stop("rendering should not happen on the error path")
    },
    showNotificationFn = function(...) {
      captured$calls <- c(captured$calls, "show")
      captured$notifications[[length(captured$notifications) + 1]] <- list(...)
    },
    logErrorFn = function(message) {
      captured$calls <- c(captured$calls, "log_error")
      captured$error <- message
    }
  )

  expect_identical(captured$calls, c("run", "log_error", "show"))
  expect_identical(
    captured$notifications[[1]][[1]],
    "Error reverting: mock revert failure"
  )
  expect_identical(captured$notifications[[1]]$type, "error")
  expect_identical(captured$error, "Error reverting: mock revert failure")
  expect_identical(completed$status, "error")
  expect_identical(
    completed$errorMessage,
    "Error reverting: mock revert failure"
  )
})

test_that("mod_prot_qc_peptide_qvalue_server wires the apply observer through the helper seam", {
  workflow_data <- shiny::reactiveValues(state_manager = list(kind = "state_manager"))
  captured <- new.env(parent = emptyenv())
  captured$apply <- NULL

  server_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_peptide_qvalue_server,
    list(
      runPeptideQvalueApplyObserver = function(workflowData,
                                               qvalueThreshold,
                                               globalQvalueThreshold,
                                               proteotypicOnly,
                                               output,
                                               qvaluePlot,
                                               omicType,
                                               experimentLabel) {
        captured$apply <<- list(
          workflowData = workflowData,
          qvalueThreshold = qvalueThreshold,
          globalQvalueThreshold = globalQvalueThreshold,
          proteotypicOnly = proteotypicOnly,
          output = output,
          qvaluePlot = qvaluePlot,
          omicType = omicType,
          experimentLabel = experimentLabel
        )
        invisible(NULL)
      }
    )
  )

  testServer(
    server_under_test,
    args = list(
      workflow_data = workflow_data,
      omic_type = "proteomics",
      experiment_label = "DIA Experiment"
    ),
    {
      session$setInputs(qvalue_threshold = 0.01)
      session$setInputs(global_qvalue_threshold = 0.02)
      session$setInputs(proteotypic_only = TRUE)
      session$flushReact()
      session$setInputs(apply_qvalue_filter = 1)
      session$flushReact()

      expect_false(is.null(captured$apply))
      expect_identical(captured$apply$workflowData, workflow_data)
      expect_identical(captured$apply$qvalueThreshold, 0.01)
      expect_identical(captured$apply$globalQvalueThreshold, 0.02)
      expect_identical(captured$apply$proteotypicOnly, TRUE)
      expect_identical(captured$apply$output, output)
      expect_identical(captured$apply$omicType, "proteomics")
      expect_identical(captured$apply$experimentLabel, "DIA Experiment")
      expect_type(captured$apply$qvaluePlot, "closure")
    }
  )
})

test_that("mod_prot_qc_peptide_qvalue_server wires the revert observer through the helper seam", {
  workflow_data <- shiny::reactiveValues(state_manager = list(kind = "state_manager"))
  captured <- new.env(parent = emptyenv())
  captured$revert <- NULL

  server_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_peptide_qvalue_server,
    list(
      runPeptideQvalueRevertObserver = function(workflowData, output) {
        captured$revert <<- list(
          workflowData = workflowData,
          output = output
        )
        invisible(NULL)
      }
    )
  )

  testServer(
    server_under_test,
    args = list(
      workflow_data = workflow_data,
      omic_type = "proteomics",
      experiment_label = "DIA Experiment"
    ),
    {
      session$setInputs(revert_qvalue = 1)
      session$flushReact()

      expect_false(is.null(captured$revert))
      expect_identical(captured$revert$workflowData, workflow_data)
      expect_identical(captured$revert$output, output)
    }
  )
})

test_that("setupPeptideIntensityServerBootstrap wires the remaining preview, observer, and plot bindings through one seam", {
  input <- new.env(parent = emptyenv())
  input$min_reps_per_group <- 2
  input$min_groups <- 3
  input$use_strict_mode <- TRUE
  input$intensity_cutoff_percentile <- 1.5
  output <- new.env(parent = emptyenv())
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- list(kind = "state_manager")
  captured <- new.env(parent = emptyenv())
  captured$previewCalls <- 0L
  captured$previewArgs <- NULL
  captured$observeEvents <- list()
  captured$apply <- NULL
  captured$revert <- NULL
  captured$drawn <- NULL

  intensity_plot <- function(value) {
    if (missing(value)) {
      return("plot_grid")
    }

    stop("setupPeptideIntensityServerBootstrap should only read from the reactive")
  }

  reactive_fn <- function(expr) {
    expr_quo <- substitute(expr)
    expr_env <- parent.frame()
    cache <- new.env(parent = emptyenv())
    cache$initialized <- FALSE

    function() {
      if (!cache$initialized) {
        cache$value <- eval(expr_quo, expr_env)
        cache$initialized <- TRUE
      }

      cache$value
    }
  }

  render_text_fn <- function(expr) {
    expr_quo <- substitute(expr)
    expr_env <- parent.frame()
    eval(expr_quo, expr_env)
  }

  observe_event_fn <- function(eventExpr, handlerExpr) {
    event_quo <- substitute(eventExpr)
    handler_quo <- substitute(handlerExpr)
    handler_env <- parent.frame()

    captured$observeEvents[[length(captured$observeEvents) + 1L]] <<- list(
      event = function() eval(event_quo, handler_env),
      run = function() eval(handler_quo, handler_env)
    )

    invisible(NULL)
  }

  render_plot_fn <- function(expr) {
    expr_quo <- substitute(expr)
    expr_env <- parent.frame()
    eval(expr_quo, expr_env)
    "rendered_plot"
  }

  bootstrap <- setupPeptideIntensityServerBootstrap(
    input = input,
    output = output,
    workflowData = workflow_data,
    intensityPlot = intensity_plot,
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    buildThresholdPreviewFn = function(workflowData, minRepsPerGroup, minGroups) {
      captured$previewCalls <- captured$previewCalls + 1L
      captured$previewArgs <- list(
        workflowData = workflowData,
        minRepsPerGroup = minRepsPerGroup,
        minGroups = minGroups
      )
      list(
        groupwiseText = "Groupwise % cutoff: 12.345%",
        maxGroupsText = "Max groups % cutoff: 67.890%"
      )
    },
    runApplyObserverFn = function(workflowData,
                                  useStrictMode,
                                  minRepsPerGroup,
                                  minGroups,
                                  intensityCutoffPercentile,
                                  output,
                                  intensityPlot,
                                  omicType,
                                  experimentLabel) {
      captured$apply <- list(
        workflowData = workflowData,
        useStrictMode = useStrictMode,
        minRepsPerGroup = minRepsPerGroup,
        minGroups = minGroups,
        intensityCutoffPercentile = intensityCutoffPercentile,
        output = output,
        intensityPlot = intensityPlot,
        omicType = omicType,
        experimentLabel = experimentLabel
      )
      invisible(NULL)
    },
    runRevertObserverFn = function(workflowData, output) {
      captured$revert <- list(
        workflowData = workflowData,
        output = output
      )
      invisible(NULL)
    },
    reactiveFn = reactive_fn,
    renderTextFn = render_text_fn,
    observeEventFn = observe_event_fn,
    renderPlotFn = render_plot_fn,
    reqFn = function(value) value,
    gridDrawFn = function(value) {
      captured$drawn <- value
      invisible(NULL)
    }
  )

  expect_identical(bootstrap$reason, "wired")
  expect_identical(bootstrap$thresholdPreview()$groupwiseText, "Groupwise % cutoff: 12.345%")
  expect_identical(output$calculated_groupwise_percent, "Groupwise % cutoff: 12.345%")
  expect_identical(output$calculated_max_groups_percent, "Max groups % cutoff: 67.890%")
  expect_identical(output$intensity_plot, "rendered_plot")
  expect_identical(captured$drawn, "plot_grid")
  expect_identical(captured$previewCalls, 1L)
  expect_identical(captured$previewArgs$workflowData, workflow_data)
  expect_identical(captured$previewArgs$minRepsPerGroup, 2)
  expect_identical(captured$previewArgs$minGroups, 3)
  expect_length(captured$observeEvents, 2L)

  captured$observeEvents[[1]]$run()

  expect_identical(captured$apply$workflowData, workflow_data)
  expect_true(captured$apply$useStrictMode)
  expect_identical(captured$apply$minRepsPerGroup, 2)
  expect_identical(captured$apply$minGroups, 3)
  expect_equal(captured$apply$intensityCutoffPercentile, 1.5)
  expect_identical(captured$apply$output, output)
  expect_identical(captured$apply$intensityPlot, intensity_plot)
  expect_identical(captured$apply$omicType, "proteomics")
  expect_identical(captured$apply$experimentLabel, "DIA Experiment")

  captured$observeEvents[[2]]$run()

  expect_identical(captured$revert$workflowData, workflow_data)
  expect_identical(captured$revert$output, output)
})

test_that("mod_prot_qc_peptide_intensity_server delegates the remaining server wiring through the bootstrap seam", {
  workflow_data <- shiny::reactiveValues(state_manager = list(kind = "state_manager"))
  captured <- new.env(parent = emptyenv())
  captured$args <- NULL

  server_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_peptide_intensity_server,
    list(
      setupPeptideIntensityServerBootstrap = function(input,
                                                      output,
                                                      workflowData,
                                                      intensityPlot,
                                                      omicType,
                                                      experimentLabel) {
        captured$args <<- list(
          input = input,
          output = output,
          workflowData = workflowData,
          intensityPlot = intensityPlot,
          omicType = omicType,
          experimentLabel = experimentLabel
        )
        invisible(list(reason = "wired"))
      }
    )
  )

  testServer(
    server_under_test,
    args = list(
      workflow_data = workflow_data,
      omic_type = "proteomics",
      experiment_label = "DIA Experiment"
    ),
    {
      expect_false(is.null(captured$args))
      expect_identical(captured$args$input, input)
      expect_identical(captured$args$output, output)
      expect_identical(captured$args$workflowData, workflow_data)
      expect_true(is.function(captured$args$intensityPlot))
      expect_null(captured$args$intensityPlot())
      expect_identical(captured$args$omicType, "proteomics")
      expect_identical(captured$args$experimentLabel, "DIA Experiment")
    }
  )
})

test_that("runPeptideImputationRevertStep reverts to the previous state in history", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$getHistory <- function() {
    c("raw_data_s4", "sample_filtered", "imputed")
  }

  reverted_to <- NULL
  workflow_data$state_manager$revertToState <- function(state) {
    reverted_to <<- state
    "reverted_object"
  }

  result <- runPeptideImputationRevertStep(workflow_data = workflow_data)

  expect_equal(result$previousState, "sample_filtered")
  expect_equal(result$revertedS4, "reverted_object")
  expect_equal(result$resultText, "Reverted to previous state: sample_filtered")
  expect_equal(reverted_to, "sample_filtered")
})

test_that("runPeptideImputationRevertStep errors when there is no prior state", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$getHistory <- function() {
    "raw_data_s4"
  }

  expect_error(
    runPeptideImputationRevertStep(workflow_data = workflow_data),
    "No previous state to revert to.",
    fixed = TRUE
  )
})

test_that("updatePeptideImputationOutputs refreshes result text and plot grid", {
  if (!methods::isClass("FakePeptideImputationState")) {
    methods::setClass("FakePeptideImputationState", slots = c(peptide_data = "data.frame"))
  }

  fake_state <- methods::new(
    "FakePeptideImputationState",
    peptide_data = data.frame(Protein.Ids = c("P1", "P2"))
  )
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$plot <- NULL
  captured$args <- NULL
  imputation_plot <- function(value) {
    if (missing(value)) {
      return(captured$plot)
    }

    captured$plot <- value
    invisible(value)
  }

  helper_under_test <- makeFunctionWithOverrides(
    updatePeptideImputationOutputs,
    list(
      renderText = function(value) value,
      updateProteinFiltering = function(data, step_name, omic_type, experiment_label, return_grid, overwrite) {
        captured$args <- list(
          data = data,
          step_name = step_name,
          omic_type = omic_type,
          experiment_label = experiment_label,
          return_grid = return_grid,
          overwrite = overwrite
        )
        "plot_grid"
      }
    )
  )

  helper_under_test(
    output = output,
    imputationPlot = imputation_plot,
    imputationResult = list(
      imputedS4 = fake_state,
      resultText = "Imputation completed"
    ),
    omicType = "proteomics",
    experimentLabel = "DIA Experiment"
  )

  expect_identical(output$imputation_results, "Imputation completed")
  expect_identical(captured$plot, "plot_grid")
  expect_identical(captured$args$data, fake_state@peptide_data)
  expect_identical(captured$args$step_name, "8_imputed")
  expect_identical(captured$args$omic_type, "proteomics")
  expect_identical(captured$args$experiment_label, "DIA Experiment")
  expect_true(captured$args$return_grid)
  expect_true(captured$args$overwrite)
})

test_that("runPeptideImputationApplyObserver delegates the apply workflow and notifications", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data_ref <- workflow_data
  output <- new.env(parent = emptyenv())
  output_ref <- output
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()
  captured$removed <- character()
  imputation_plot <- function(...) "plot_reactive"

  completed <- runPeptideImputationApplyObserver(
    workflow_data = workflow_data,
    proportionMissingValues = 0.4,
    output = output,
    imputationPlot = imputation_plot,
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    runImputationStepFn = function(workflow_data, proportionMissingValues) {
      captured$calls <- c(captured$calls, "run")
      expect_identical(workflow_data, workflow_data_ref)
      expect_equal(proportionMissingValues, 0.4)
      list(
        imputedS4 = "imputed_s4",
        resultText = "Imputation completed"
      )
    },
    updateOutputsFn = function(output, imputationPlot, imputationResult, omicType, experimentLabel) {
      captured$calls <- c(captured$calls, "update")
      expect_identical(output, output_ref)
      expect_identical(imputationPlot, imputation_plot)
      expect_identical(imputationResult$imputedS4, "imputed_s4")
      expect_identical(imputationResult$resultText, "Imputation completed")
      expect_identical(omicType, "proteomics")
      expect_identical(experimentLabel, "DIA Experiment")
      "plot_grid"
    },
    showNotificationFn = function(...) {
      captured$calls <- c(captured$calls, "show")
      captured$notifications[[length(captured$notifications) + 1]] <- list(...)
    },
    removeNotificationFn = function(id) {
      captured$calls <- c(captured$calls, "remove")
      captured$removed <- c(captured$removed, id)
    },
    logInfoFn = function(message) {
      captured$calls <- c(captured$calls, "log_info")
      captured$info <- message
    },
    logErrorFn = function(...) {
      stop("error logging should not run on the happy path")
    }
  )

  expect_identical(
    captured$calls,
    c("show", "run", "update", "log_info", "remove", "show")
  )
  expect_identical(captured$notifications[[1]][[1]], "Applying missing value imputation...")
  expect_identical(captured$notifications[[1]]$id, "imputation_working")
  expect_null(captured$notifications[[1]]$duration)
  expect_identical(
    captured$notifications[[2]][[1]],
    "Missing value imputation applied successfully"
  )
  expect_identical(captured$notifications[[2]]$type, "message")
  expect_identical(captured$removed, "imputation_working")
  expect_identical(captured$info, "Missing value imputation applied successfully")
  expect_identical(completed$status, "success")
  expect_identical(completed$imputationResult$imputedS4, "imputed_s4")
  expect_identical(completed$plotGrid, "plot_grid")
})

test_that("runPeptideImputationApplyObserver reports apply errors and clears the working notification", {
  workflow_data <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()
  captured$removed <- character()
  imputation_plot <- function(...) "plot_reactive"

  completed <- runPeptideImputationApplyObserver(
    workflow_data = workflow_data,
    proportionMissingValues = 0.4,
    output = output,
    imputationPlot = imputation_plot,
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    runImputationStepFn = function(...) {
      captured$calls <- c(captured$calls, "run")
      stop("mock failure")
    },
    updateOutputsFn = function(...) {
      stop("output refresh should not run after an imputation error")
    },
    showNotificationFn = function(...) {
      captured$calls <- c(captured$calls, "show")
      captured$notifications[[length(captured$notifications) + 1]] <- list(...)
    },
    removeNotificationFn = function(id) {
      captured$calls <- c(captured$calls, "remove")
      captured$removed <- c(captured$removed, id)
    },
    logInfoFn = function(...) {
      stop("success logging should not run on the error path")
    },
    logErrorFn = function(message) {
      captured$calls <- c(captured$calls, "log_error")
      captured$error <- message
    }
  )

  expect_identical(
    captured$calls,
    c("show", "run", "log_error", "show", "remove")
  )
  expect_identical(captured$notifications[[1]][[1]], "Applying missing value imputation...")
  expect_identical(captured$notifications[[1]]$id, "imputation_working")
  expect_identical(
    captured$notifications[[2]][[1]],
    "Error applying missing value imputation: mock failure"
  )
  expect_identical(captured$notifications[[2]]$type, "error")
  expect_identical(captured$notifications[[2]]$duration, 15)
  expect_identical(captured$removed, "imputation_working")
  expect_identical(captured$error, "Error applying missing value imputation: mock failure")
  expect_identical(completed$status, "error")
  expect_identical(completed$errorMessage, "Error applying missing value imputation: mock failure")
})

test_that("runPeptideImputationRevertObserver delegates the revert workflow and success notification", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data_ref <- workflow_data
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()

  completed <- runPeptideImputationRevertObserver(
    workflow_data = workflow_data,
    output = output,
    runRevertStepFn = function(workflow_data) {
      captured$calls <- c(captured$calls, "run")
      expect_identical(workflow_data, workflow_data_ref)
      list(
        previousState = "sample_filtered",
        resultText = "Reverted to previous state: sample_filtered"
      )
    },
    renderTextFn = function(value) {
      captured$calls <- c(captured$calls, "render")
      value
    },
    showNotificationFn = function(...) {
      captured$calls <- c(captured$calls, "show")
      captured$notifications[[length(captured$notifications) + 1]] <- list(...)
    },
    logErrorFn = function(...) {
      stop("error logging should not run on the happy path")
    }
  )

  expect_identical(captured$calls, c("run", "render", "show"))
  expect_identical(output$imputation_results, "Reverted to previous state: sample_filtered")
  expect_identical(captured$notifications[[1]][[1]], "Reverted successfully")
  expect_identical(captured$notifications[[1]]$type, "message")
  expect_identical(completed$status, "success")
  expect_identical(completed$revertResult$previousState, "sample_filtered")
})

test_that("runPeptideImputationRevertObserver reports revert errors without rendering results", {
  workflow_data <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()

  completed <- runPeptideImputationRevertObserver(
    workflow_data = workflow_data,
    output = output,
    runRevertStepFn = function(...) {
      captured$calls <- c(captured$calls, "run")
      stop("mock revert failure")
    },
    renderTextFn = function(...) {
      stop("renderText should not run on the error path")
    },
    showNotificationFn = function(...) {
      captured$calls <- c(captured$calls, "show")
      captured$notifications[[length(captured$notifications) + 1]] <- list(...)
    },
    logErrorFn = function(message) {
      captured$calls <- c(captured$calls, "log_error")
      captured$error <- message
    }
  )

  expect_identical(captured$calls, c("run", "log_error", "show"))
  expect_false(exists("imputation_results", envir = output, inherits = FALSE))
  expect_identical(captured$notifications[[1]][[1]], "Error reverting: mock revert failure")
  expect_identical(captured$notifications[[1]]$type, "error")
  expect_identical(captured$error, "Error reverting: mock revert failure")
  expect_identical(completed$status, "error")
  expect_identical(completed$errorMessage, "Error reverting: mock revert failure")
})

test_that("bindPeptideImputationPlot binds the rendered plot output to the stored grid", {
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$drawn <- NULL

  imputation_plot <- function(value) {
    if (missing(value)) {
      return("plot_grid")
    }

    stop("bindPeptideImputationPlot should only read from the reactive")
  }

  helper_under_test <- makeFunctionWithOverrides(
    bindPeptideImputationPlot,
    list(
      renderPlot = function(expr) {
        force(expr)
        "rendered_plot"
      },
      req = function(value) value,
      grid.draw = function(value) {
        captured$drawn <- value
        invisible(NULL)
      }
    )
  )

  helper_under_test(output = output, imputationPlot = imputation_plot)

  expect_identical(output$imputation_plot, "rendered_plot")
  expect_identical(captured$drawn, "plot_grid")
})

test_that("mod_prot_qc_peptide_impute_server wires apply and revert observers through the helper seams", {
  workflow_data <- shiny::reactiveValues(state_manager = list(kind = "state_manager"))
  captured <- new.env(parent = emptyenv())
  captured$apply <- NULL
  captured$revert <- NULL
  captured$bind <- NULL

  server_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_peptide_impute_server,
    list(
      runPeptideImputationApplyObserver = function(workflow_data,
                                                   proportionMissingValues,
                                                   output,
                                                   imputationPlot,
                                                   omicType,
                                                   experimentLabel) {
        captured$apply <<- list(
          workflow_data = workflow_data,
          proportionMissingValues = proportionMissingValues,
          output = output,
          imputationPlot = imputationPlot,
          omicType = omicType,
          experimentLabel = experimentLabel
        )
        invisible(NULL)
      },
      runPeptideImputationRevertObserver = function(workflow_data, output) {
        captured$revert <<- list(
          workflow_data = workflow_data,
          output = output
        )
        invisible(NULL)
      },
      bindPeptideImputationPlot = function(output, imputationPlot) {
        captured$bind <<- list(
          output = output,
          imputationPlot = imputationPlot
        )
        invisible(NULL)
      }
    )
  )

  testServer(
    server_under_test,
    args = list(
      workflow_data = workflow_data,
      omic_type = "proteomics",
      experiment_label = "DIA Experiment"
    ),
    {
      expect_false(is.null(captured$bind))
      expect_identical(captured$bind$output, output)

      session$setInputs(proportion_missing_values = 0.3)
      session$flushReact()
      session$setInputs(apply_imputation = 1)
      session$flushReact()

      expect_false(is.null(captured$apply))
      expect_identical(captured$apply$workflow_data, workflow_data)
      expect_identical(captured$apply$output, output)
      expect_identical(captured$apply$imputationPlot, captured$bind$imputationPlot)
      expect_identical(captured$apply$omicType, "proteomics")
      expect_identical(captured$apply$experimentLabel, "DIA Experiment")
      expect_equal(captured$apply$proportionMissingValues, 0.3)

      session$setInputs(revert_imputation = 1)
      session$flushReact()

      expect_false(is.null(captured$revert))
      expect_identical(captured$revert$workflow_data, workflow_data)
      expect_identical(captured$revert$output, output)
    }
  )
})

test_that("runProteinPeptideApplyStep applies both cutoffs and saves the filtered state", {
  if (!methods::isClass("FakeProteinPeptideState")) {
    methods::setClass(
      "FakeProteinPeptideState",
      slots = c(args = "list", peptide_data = "data.frame")
    )
  }

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$getState <- function() {
    methods::new(
      "FakeProteinPeptideState",
      args = list(
        filterMinNumPeptidesPerProtein = list(
          peptides_per_protein_cutoff = 1,
          peptidoforms_per_protein_cutoff = 1
        )
      ),
      peptide_data = data.frame(Protein.Ids = c("P0"))
    )
  }

  captured <- new.env(parent = emptyenv())
  captured$configUpdates <- list()
  captured$saveState <- NULL
  workflow_data$state_manager$saveState <- function(...) {
    captured$saveState <- list(...)
  }

  filtered_s4 <- methods::new(
    "FakeProteinPeptideState",
    args = list(
      filterMinNumPeptidesPerProtein = list(
        peptides_per_protein_cutoff = 3,
        peptidoforms_per_protein_cutoff = 2
      )
    ),
    peptide_data = data.frame(Protein.Ids = c("P1", "P1", "P2"))
  )
  timestamp <- as.POSIXct("2026-04-16 15:00:00", tz = "UTC")

  result <- runProteinPeptideApplyStep(
    workflowData = workflow_data,
    minPeptidesPerProtein = 3,
    minPeptidoformsPerProtein = 2,
    updateConfigParameterFn = function(theObject, function_name, parameter_name, new_value) {
      captured$configUpdates[[length(captured$configUpdates) + 1]] <<- list(
        function_name = function_name,
        parameter_name = parameter_name,
        new_value = new_value
      )
      updated_args <- theObject@args
      updated_args$filterMinNumPeptidesPerProtein[[parameter_name]] <- new_value
      methods::new(
        "FakeProteinPeptideState",
        args = updated_args,
        peptide_data = theObject@peptide_data
      )
    },
    filterMinNumPeptidesPerProteinFn = function(theObject) {
      captured$filterInput <<- theObject
      filtered_s4
    },
    logInfoFn = function(message) {
      captured$info <- message
    },
    nowFn = function() timestamp
  )

  expect_identical(
    captured$info,
    "QC Step: Applying protein peptide count filter (min: 3)"
  )
  expect_length(captured$configUpdates, 2)
  expect_identical(
    captured$configUpdates[[1]]$parameter_name,
    "peptides_per_protein_cutoff"
  )
  expect_identical(captured$configUpdates[[1]]$new_value, 3)
  expect_identical(
    captured$configUpdates[[2]]$parameter_name,
    "peptidoforms_per_protein_cutoff"
  )
  expect_identical(captured$configUpdates[[2]]$new_value, 2)
  expect_identical(
    captured$filterInput@args$filterMinNumPeptidesPerProtein$peptides_per_protein_cutoff,
    3
  )
  expect_identical(
    captured$filterInput@args$filterMinNumPeptidesPerProtein$peptidoforms_per_protein_cutoff,
    2
  )
  expect_identical(workflow_data$qc_params$peptide_qc$protein_peptide_filter$min_peptides_per_protein, 3)
  expect_identical(workflow_data$qc_params$peptide_qc$protein_peptide_filter$min_peptidoforms_per_protein, 2)
  expect_identical(workflow_data$qc_params$peptide_qc$protein_peptide_filter$timestamp, timestamp)
  expect_identical(captured$saveState$state_name, "protein_peptide_filtered")
  expect_identical(captured$saveState$s4_data_object, filtered_s4)
  expect_identical(captured$saveState$config_object$min_peptides_per_protein, 3)
  expect_identical(captured$saveState$config_object$min_peptidoforms_per_protein, 2)
  expect_identical(
    captured$saveState$description,
    "Applied minimum peptides per protein filter"
  )
  expect_identical(result$filteredS4, filtered_s4)
  expect_match(result$resultText, "Proteins remaining: 2", fixed = TRUE)
  expect_match(result$resultText, "Min peptides per protein: 3", fixed = TRUE)
  expect_match(result$resultText, "Min peptidoforms per protein: 2", fixed = TRUE)
})

test_that("updateProteinPeptideOutputs refreshes result text and plot grid", {
  if (!methods::isClass("FakeProteinPeptideState")) {
    methods::setClass(
      "FakeProteinPeptideState",
      slots = c(args = "list", peptide_data = "data.frame")
    )
  }

  protein_state <- methods::new(
    "FakeProteinPeptideState",
    args = list(filterMinNumPeptidesPerProtein = list()),
    peptide_data = data.frame(Protein.Ids = c("P1", "P2"))
  )
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$args <- NULL
  captured$plot <- NULL
  protein_peptide_plot <- function(value) {
    captured$plot <- value
    invisible(value)
  }

  plot_grid <- updateProteinPeptideOutputs(
    output = output,
    proteinPeptidePlot = protein_peptide_plot,
    proteinPeptideResult = list(
      filteredS4 = protein_state,
      resultText = "Protein peptide completed"
    ),
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    renderTextFn = function(value) value,
    updateProteinFilteringFn = function(data, step_name, omic_type, experiment_label, return_grid, overwrite) {
      captured$args <- list(
        data = data,
        step_name = step_name,
        omic_type = omic_type,
        experiment_label = experiment_label,
        return_grid = return_grid,
        overwrite = overwrite
      )
      "plot_grid"
    }
  )

  expect_identical(output$protein_peptida_results, "Protein peptide completed")
  expect_identical(plot_grid, "plot_grid")
  expect_identical(captured$plot, "plot_grid")
  expect_identical(captured$args$data, protein_state@peptide_data)
  expect_identical(captured$args$step_name, "5_protein_peptide_filtered")
  expect_identical(captured$args$omic_type, "proteomics")
  expect_identical(captured$args$experiment_label, "DIA Experiment")
  expect_true(captured$args$return_grid)
  expect_true(captured$args$overwrite)
})

test_that("runProteinPeptideApplyObserver delegates the apply workflow and notifications", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data_ref <- workflow_data
  output <- new.env(parent = emptyenv())
  output_ref <- output
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()
  captured$removed <- character()
  protein_peptide_plot <- function(...) "plot_reactive"

  completed <- runProteinPeptideApplyObserver(
    workflowData = workflow_data,
    minPeptidesPerProtein = 3,
    minPeptidoformsPerProtein = 2,
    output = output,
    proteinPeptidePlot = protein_peptide_plot,
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    runApplyStepFn = function(workflowData, minPeptidesPerProtein, minPeptidoformsPerProtein) {
      captured$calls <- c(captured$calls, "run")
      expect_identical(workflowData, workflow_data_ref)
      expect_identical(minPeptidesPerProtein, 3)
      expect_identical(minPeptidoformsPerProtein, 2)
      list(
        filteredS4 = "filtered_s4",
        resultText = "Protein peptide completed"
      )
    },
    updateOutputsFn = function(output, proteinPeptidePlot, proteinPeptideResult, omicType, experimentLabel) {
      captured$calls <- c(captured$calls, "update")
      expect_identical(output, output_ref)
      expect_identical(proteinPeptidePlot, protein_peptide_plot)
      expect_identical(proteinPeptideResult$filteredS4, "filtered_s4")
      expect_identical(proteinPeptideResult$resultText, "Protein peptide completed")
      expect_identical(omicType, "proteomics")
      expect_identical(experimentLabel, "DIA Experiment")
      "plot_grid"
    },
    showNotificationFn = function(...) {
      captured$calls <- c(captured$calls, "show")
      captured$notifications[[length(captured$notifications) + 1]] <- list(...)
    },
    removeNotificationFn = function(id) {
      captured$calls <- c(captured$calls, "remove")
      captured$removed <- c(captured$removed, id)
    },
    logInfoFn = function(message) {
      captured$calls <- c(captured$calls, "log_info")
      captured$info <- message
    },
    logErrorFn = function(...) {
      stop("error logging should not run on the happy path")
    }
  )

  expect_identical(
    captured$calls,
    c("show", "run", "update", "log_info", "remove", "show")
  )
  expect_identical(
    captured$notifications[[1]][[1]],
    "Applying protein peptide count filter..."
  )
  expect_identical(captured$notifications[[1]]$id, "protein_peptide_working")
  expect_null(captured$notifications[[1]]$duration)
  expect_identical(
    captured$notifications[[2]][[1]],
    "Protein peptide count filter applied successfully"
  )
  expect_identical(captured$notifications[[2]]$type, "message")
  expect_identical(captured$removed, "protein_peptide_working")
  expect_identical(
    captured$info,
    "Protein peptide count filter applied successfully"
  )
  expect_identical(completed$status, "success")
  expect_identical(completed$proteinPeptideResult$filteredS4, "filtered_s4")
  expect_identical(completed$plotGrid, "plot_grid")
})

test_that("runProteinPeptideApplyObserver reports apply errors and clears the working notification", {
  workflow_data <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()
  captured$removed <- character()
  protein_peptide_plot <- function(...) "plot_reactive"

  completed <- runProteinPeptideApplyObserver(
    workflowData = workflow_data,
    minPeptidesPerProtein = 3,
    minPeptidoformsPerProtein = 2,
    output = output,
    proteinPeptidePlot = protein_peptide_plot,
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    runApplyStepFn = function(...) {
      captured$calls <- c(captured$calls, "run")
      stop("mock failure")
    },
    updateOutputsFn = function(...) {
      stop("output refresh should not run after a protein peptide error")
    },
    showNotificationFn = function(...) {
      captured$calls <- c(captured$calls, "show")
      captured$notifications[[length(captured$notifications) + 1]] <- list(...)
    },
    removeNotificationFn = function(id) {
      captured$calls <- c(captured$calls, "remove")
      captured$removed <- c(captured$removed, id)
    },
    logInfoFn = function(...) {
      stop("success logging should not run on the error path")
    },
    logErrorFn = function(message) {
      captured$calls <- c(captured$calls, "log_error")
      captured$error <- message
    }
  )

  expect_identical(
    captured$calls,
    c("show", "run", "log_error", "show", "remove")
  )
  expect_identical(
    captured$notifications[[1]][[1]],
    "Applying protein peptide count filter..."
  )
  expect_identical(captured$notifications[[1]]$id, "protein_peptide_working")
  expect_identical(
    captured$notifications[[2]][[1]],
    "Error applying protein peptide count filter: mock failure"
  )
  expect_identical(captured$notifications[[2]]$type, "error")
  expect_identical(captured$notifications[[2]]$duration, 15)
  expect_identical(captured$removed, "protein_peptide_working")
  expect_identical(
    captured$error,
    "Error applying protein peptide count filter: mock failure"
  )
  expect_identical(completed$status, "error")
  expect_identical(
    completed$errorMessage,
    "Error applying protein peptide count filter: mock failure"
  )
})

test_that("runProteinPeptideRevertStep reverts to the previous state in history", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$getHistory <- function() {
    c("raw_data_s4", "qvalue_filtered", "protein_peptide_filtered")
  }

  reverted_to <- NULL
  workflow_data$state_manager$revertToState <- function(state) {
    reverted_to <<- state
    "reverted_object"
  }

  result <- runProteinPeptideRevertStep(workflowData = workflow_data)

  expect_identical(result$previousState, "qvalue_filtered")
  expect_identical(result$revertedS4, "reverted_object")
  expect_identical(result$resultText, "Reverted to previous state: qvalue_filtered")
  expect_identical(reverted_to, "qvalue_filtered")
})

test_that("runProteinPeptideRevertStep errors when there is no prior state", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$getHistory <- function() {
    "raw_data_s4"
  }

  expect_error(
    runProteinPeptideRevertStep(workflowData = workflow_data),
    "No previous state to revert to.",
    fixed = TRUE
  )
})

test_that("runPeptideReplicateRevertStep reverts to the previous state in history", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$getHistory <- function() {
    c("raw_data_s4", "sample_filtered", "replicate_filtered")
  }

  reverted_to <- NULL
  workflow_data$state_manager$revertToState <- function(state) {
    reverted_to <<- state
    "reverted_object"
  }

  result <- runPeptideReplicateRevertStep(workflowData = workflow_data)

  expect_identical(result$previousState, "sample_filtered")
  expect_identical(result$revertedS4, "reverted_object")
  expect_identical(result$resultText, "Reverted to previous state: sample_filtered")
  expect_identical(reverted_to, "sample_filtered")
})

test_that("runPeptideReplicateRevertStep errors when there is no prior state", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$getHistory <- function() {
    "raw_data_s4"
  }

  expect_error(
    runPeptideReplicateRevertStep(workflowData = workflow_data),
    "No previous state to revert to.",
    fixed = TRUE
  )
})

test_that("runPeptideReplicateApplyStep saves the filtered state and QC parameters", {
  if (!methods::isClass("FakePeptideReplicateState")) {
    methods::setClass(
      "FakePeptideReplicateState",
      slots = c(peptide_data = "data.frame")
    )
  }

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  current_state <- methods::new(
    "FakePeptideReplicateState",
    peptide_data = data.frame(Protein.Ids = c("P0"))
  )
  workflow_data$state_manager$getState <- function() {
    current_state
  }

  captured <- new.env(parent = emptyenv())
  captured$saveState <- NULL
  workflow_data$state_manager$saveState <- function(...) {
    captured$saveState <- list(...)
  }

  filtered_s4 <- methods::new(
    "FakePeptideReplicateState",
    peptide_data = data.frame(Protein.Ids = c("P1", "P1", "P2"))
  )
  timestamp <- as.POSIXct("2026-04-16 15:30:00", tz = "UTC")

  result <- runPeptideReplicateApplyStep(
    workflowData = workflow_data,
    replicateGroupColumn = "replicates",
    removePeptidesWithOnlyOneReplicateFn = function(theObject, replicate_group_column) {
      captured$filterInput <<- list(
        theObject = theObject,
        replicate_group_column = replicate_group_column
      )
      filtered_s4
    },
    logInfoFn = function(message) {
      captured$info <- message
    },
    nowFn = function() timestamp
  )

  expect_identical(
    captured$info,
    "QC Step: Applying replicate filter (column: replicates)"
  )
  expect_identical(captured$filterInput$theObject, current_state)
  expect_identical(captured$filterInput$replicate_group_column, "replicates")
  expect_identical(
    workflow_data$qc_params$peptide_qc$replicate_filter$replicate_group_column,
    "replicates"
  )
  expect_identical(
    workflow_data$qc_params$peptide_qc$replicate_filter$timestamp,
    timestamp
  )
  expect_identical(captured$saveState$state_name, "replicate_filtered")
  expect_identical(captured$saveState$s4_data_object, filtered_s4)
  expect_identical(
    captured$saveState$config_object$replicate_group_column,
    "replicates"
  )
  expect_identical(
    captured$saveState$description,
    "Applied replicate filter (removed single-replicate peptides)"
  )
  expect_identical(result$filteredS4, filtered_s4)
  expect_match(result$resultText, "Proteins remaining: 2", fixed = TRUE)
  expect_match(result$resultText, "Replicate group column: replicates", fixed = TRUE)
})

test_that("updatePeptideReplicateOutputs refreshes result text and plot grid", {
  if (!methods::isClass("FakePeptideReplicateState")) {
    methods::setClass(
      "FakePeptideReplicateState",
      slots = c(peptide_data = "data.frame")
    )
  }

  replicate_state <- methods::new(
    "FakePeptideReplicateState",
    peptide_data = data.frame(Protein.Ids = c("P1", "P2"))
  )
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$args <- NULL
  captured$plot <- NULL
  replicate_plot <- function(value) {
    captured$plot <- value
    invisible(value)
  }

  plot_grid <- updatePeptideReplicateOutputs(
    output = output,
    replicatePlot = replicate_plot,
    replicateResult = list(
      filteredS4 = replicate_state,
      resultText = "Replicate completed"
    ),
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    renderTextFn = function(value) value,
    updateProteinFilteringFn = function(data, step_name, omic_type, experiment_label, return_grid, overwrite) {
      captured$args <- list(
        data = data,
        step_name = step_name,
        omic_type = omic_type,
        experiment_label = experiment_label,
        return_grid = return_grid,
        overwrite = overwrite
      )
      "plot_grid"
    }
  )

  expect_identical(output$replicate_results, "Replicate completed")
  expect_identical(plot_grid, "plot_grid")
  expect_identical(captured$plot, "plot_grid")
  expect_identical(captured$args$data, replicate_state@peptide_data)
  expect_identical(captured$args$step_name, "7_replicate_filtered")
  expect_identical(captured$args$omic_type, "proteomics")
  expect_identical(captured$args$experiment_label, "DIA Experiment")
  expect_true(captured$args$return_grid)
  expect_true(captured$args$overwrite)
})

test_that("runPeptideReplicateApplyObserver delegates the apply workflow and notifications", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data_ref <- workflow_data
  output <- new.env(parent = emptyenv())
  output_ref <- output
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()
  captured$removed <- character()
  replicate_plot <- function(...) "plot_reactive"

  completed <- runPeptideReplicateApplyObserver(
    workflowData = workflow_data,
    replicateGroupColumn = "replicates",
    output = output,
    replicatePlot = replicate_plot,
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    runApplyStepFn = function(workflowData, replicateGroupColumn) {
      captured$calls <- c(captured$calls, "run")
      expect_identical(workflowData, workflow_data_ref)
      expect_identical(replicateGroupColumn, "replicates")
      list(
        filteredS4 = "filtered_s4",
        resultText = "Replicate completed"
      )
    },
    updateOutputsFn = function(output, replicatePlot, replicateResult, omicType, experimentLabel) {
      captured$calls <- c(captured$calls, "update")
      expect_identical(output, output_ref)
      expect_identical(replicatePlot, replicate_plot)
      expect_identical(replicateResult$filteredS4, "filtered_s4")
      expect_identical(replicateResult$resultText, "Replicate completed")
      expect_identical(omicType, "proteomics")
      expect_identical(experimentLabel, "DIA Experiment")
      "plot_grid"
    },
    showNotificationFn = function(...) {
      captured$calls <- c(captured$calls, "show")
      captured$notifications[[length(captured$notifications) + 1]] <- list(...)
    },
    removeNotificationFn = function(id) {
      captured$calls <- c(captured$calls, "remove")
      captured$removed <- c(captured$removed, id)
    },
    logInfoFn = function(message) {
      captured$calls <- c(captured$calls, "log_info")
      captured$info <- message
    },
    logErrorFn = function(...) {
      stop("error logging should not run on the happy path")
    }
  )

  expect_identical(
    captured$calls,
    c("show", "run", "update", "log_info", "remove", "show")
  )
  expect_identical(captured$notifications[[1]][[1]], "Applying replicate filter...")
  expect_identical(captured$notifications[[1]]$id, "replicate_working")
  expect_null(captured$notifications[[1]]$duration)
  expect_identical(
    captured$notifications[[2]][[1]],
    "Replicate filter applied successfully"
  )
  expect_identical(captured$notifications[[2]]$type, "message")
  expect_identical(captured$removed, "replicate_working")
  expect_identical(captured$info, "Replicate filter applied successfully")
  expect_identical(completed$status, "success")
  expect_identical(completed$replicateResult$filteredS4, "filtered_s4")
  expect_identical(completed$plotGrid, "plot_grid")
})

test_that("runPeptideReplicateApplyObserver reports apply errors and clears the working notification", {
  workflow_data <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()
  captured$removed <- character()
  replicate_plot <- function(...) "plot_reactive"

  completed <- runPeptideReplicateApplyObserver(
    workflowData = workflow_data,
    replicateGroupColumn = "replicates",
    output = output,
    replicatePlot = replicate_plot,
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    runApplyStepFn = function(...) {
      captured$calls <- c(captured$calls, "run")
      stop("mock failure")
    },
    updateOutputsFn = function(...) {
      stop("output refresh should not run after a replicate error")
    },
    showNotificationFn = function(...) {
      captured$calls <- c(captured$calls, "show")
      captured$notifications[[length(captured$notifications) + 1]] <- list(...)
    },
    removeNotificationFn = function(id) {
      captured$calls <- c(captured$calls, "remove")
      captured$removed <- c(captured$removed, id)
    },
    logInfoFn = function(...) {
      stop("success logging should not run on the error path")
    },
    logErrorFn = function(message) {
      captured$calls <- c(captured$calls, "log_error")
      captured$error <- message
    }
  )

  expect_identical(
    captured$calls,
    c("show", "run", "log_error", "show", "remove")
  )
  expect_identical(captured$notifications[[1]][[1]], "Applying replicate filter...")
  expect_identical(captured$notifications[[1]]$id, "replicate_working")
  expect_identical(
    captured$notifications[[2]][[1]],
    "Error applying replicate filter: mock failure"
  )
  expect_identical(captured$notifications[[2]]$type, "error")
  expect_identical(captured$notifications[[2]]$duration, 15)
  expect_identical(captured$removed, "replicate_working")
  expect_identical(
    captured$error,
    "Error applying replicate filter: mock failure"
  )
  expect_identical(completed$status, "error")
  expect_identical(
    completed$errorMessage,
    "Error applying replicate filter: mock failure"
  )
})

test_that("runPeptideReplicateRevertObserver delegates the revert workflow and success notification", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data_ref <- workflow_data
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()

  completed <- runPeptideReplicateRevertObserver(
    workflowData = workflow_data,
    output = output,
    runRevertStepFn = function(workflowData) {
      captured$calls <- c(captured$calls, "run")
      expect_identical(workflowData, workflow_data_ref)
      list(
        previousState = "sample_filtered",
        resultText = "Reverted to previous state: sample_filtered"
      )
    },
    renderTextFn = function(value) {
      captured$calls <- c(captured$calls, "render")
      value
    },
    showNotificationFn = function(...) {
      captured$calls <- c(captured$calls, "show")
      captured$notifications[[length(captured$notifications) + 1]] <- list(...)
    },
    logErrorFn = function(...) {
      stop("error logging should not run on the happy path")
    }
  )

  expect_identical(captured$calls, c("run", "render", "show"))
  expect_identical(
    output$replicate_results,
    "Reverted to previous state: sample_filtered"
  )
  expect_identical(captured$notifications[[1]][[1]], "Reverted successfully")
  expect_identical(captured$notifications[[1]]$type, "message")
  expect_identical(completed$status, "success")
  expect_identical(
    completed$revertResult$resultText,
    "Reverted to previous state: sample_filtered"
  )
})

test_that("runPeptideReplicateRevertObserver reports revert errors without rendering results", {
  workflow_data <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()

  completed <- runPeptideReplicateRevertObserver(
    workflowData = workflow_data,
    output = output,
    runRevertStepFn = function(...) {
      captured$calls <- c(captured$calls, "run")
      stop("mock revert failure")
    },
    renderTextFn = function(...) {
      stop("renderText should not run on the error path")
    },
    showNotificationFn = function(...) {
      captured$calls <- c(captured$calls, "show")
      captured$notifications[[length(captured$notifications) + 1]] <- list(...)
    },
    logErrorFn = function(message) {
      captured$calls <- c(captured$calls, "log_error")
      captured$error <- message
    }
  )

  expect_identical(captured$calls, c("run", "log_error", "show"))
  expect_identical(captured$error, "Error reverting: mock revert failure")
  expect_identical(
    captured$notifications[[1]][[1]],
    "Error reverting: mock revert failure"
  )
  expect_identical(captured$notifications[[1]]$type, "error")
  expect_false(exists("replicate_results", envir = output, inherits = FALSE))
  expect_identical(completed$status, "error")
  expect_identical(completed$errorMessage, "Error reverting: mock revert failure")
})

test_that("runProteinPeptideRevertObserver delegates the revert workflow and success notification", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data_ref <- workflow_data
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()

  completed <- runProteinPeptideRevertObserver(
    workflowData = workflow_data,
    output = output,
    runRevertStepFn = function(workflowData) {
      captured$calls <- c(captured$calls, "run")
      expect_identical(workflowData, workflow_data_ref)
      list(
        previousState = "qvalue_filtered",
        resultText = "Reverted to previous state: qvalue_filtered"
      )
    },
    renderTextFn = function(value) {
      captured$calls <- c(captured$calls, "render")
      value
    },
    showNotificationFn = function(...) {
      captured$calls <- c(captured$calls, "show")
      captured$notifications[[length(captured$notifications) + 1]] <- list(...)
    },
    logErrorFn = function(...) {
      stop("error logging should not run on the happy path")
    }
  )

  expect_identical(captured$calls, c("run", "render", "show"))
  expect_identical(
    output$protein_peptida_results,
    "Reverted to previous state: qvalue_filtered"
  )
  expect_identical(captured$notifications[[1]][[1]], "Reverted successfully")
  expect_identical(captured$notifications[[1]]$type, "message")
  expect_identical(completed$status, "success")
  expect_identical(
    completed$revertResult$resultText,
    "Reverted to previous state: qvalue_filtered"
  )
})

test_that("runProteinPeptideRevertObserver reports revert errors without rendering results", {
  workflow_data <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()

  completed <- runProteinPeptideRevertObserver(
    workflowData = workflow_data,
    output = output,
    runRevertStepFn = function(...) {
      captured$calls <- c(captured$calls, "run")
      stop("mock revert failure")
    },
    renderTextFn = function(...) {
      stop("renderText should not run on the error path")
    },
    showNotificationFn = function(...) {
      captured$calls <- c(captured$calls, "show")
      captured$notifications[[length(captured$notifications) + 1]] <- list(...)
    },
    logErrorFn = function(message) {
      captured$calls <- c(captured$calls, "log_error")
      captured$error <- message
    }
  )

  expect_identical(captured$calls, c("run", "log_error", "show"))
  expect_identical(captured$error, "Error reverting: mock revert failure")
  expect_identical(
    captured$notifications[[1]][[1]],
    "Error reverting: mock revert failure"
  )
  expect_identical(captured$notifications[[1]]$type, "error")
  expect_false(exists("protein_peptida_results", envir = output, inherits = FALSE))
  expect_identical(completed$status, "error")
  expect_identical(completed$errorMessage, "Error reverting: mock revert failure")
})

test_that("mod_prot_qc_peptide_protein_server wires the apply observer through the helper seam", {
  workflow_data <- shiny::reactiveValues(state_manager = list(kind = "state_manager"))
  captured <- new.env(parent = emptyenv())
  captured$apply <- NULL

  server_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_peptide_protein_server,
    list(
      runProteinPeptideApplyObserver = function(workflowData,
                                                minPeptidesPerProtein,
                                                minPeptidoformsPerProtein,
                                                output,
                                                proteinPeptidePlot,
                                                omicType,
                                                experimentLabel) {
        captured$apply <<- list(
          workflowData = workflowData,
          minPeptidesPerProtein = minPeptidesPerProtein,
          minPeptidoformsPerProtein = minPeptidoformsPerProtein,
          output = output,
          proteinPeptidePlot = proteinPeptidePlot,
          omicType = omicType,
          experimentLabel = experimentLabel
        )
        invisible(NULL)
      }
    )
  )

  testServer(
    server_under_test,
    args = list(
      workflow_data = workflow_data,
      omic_type = "proteomics",
      experiment_label = "DIA Experiment"
    ),
    {
      session$setInputs(min_peptides_per_protein = 3)
      session$setInputs(min_peptidoforms_per_protein = 2)
      session$flushReact()
      session$setInputs(apply_protein_peptide_filter = 1)
      session$flushReact()

      expect_false(is.null(captured$apply))
      expect_identical(captured$apply$workflowData, workflow_data)
      expect_identical(captured$apply$minPeptidesPerProtein, 3)
      expect_identical(captured$apply$minPeptidoformsPerProtein, 2)
      expect_identical(captured$apply$output, output)
      expect_identical(captured$apply$omicType, "proteomics")
      expect_identical(captured$apply$experimentLabel, "DIA Experiment")
      expect_type(captured$apply$proteinPeptidePlot, "closure")
    }
  )
})

test_that("mod_prot_qc_peptide_replicate_server wires the apply observer through the helper seam", {
  workflow_data <- shiny::reactiveValues(state_manager = list(kind = "state_manager"))
  captured <- new.env(parent = emptyenv())
  captured$apply <- NULL

  server_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_peptide_replicate_server,
    list(
      runPeptideReplicateApplyObserver = function(workflowData,
                                                  replicateGroupColumn,
                                                  output,
                                                  replicatePlot,
                                                  omicType,
                                                  experimentLabel) {
        captured$apply <<- list(
          workflowData = workflowData,
          replicateGroupColumn = replicateGroupColumn,
          output = output,
          replicatePlot = replicatePlot,
          omicType = omicType,
          experimentLabel = experimentLabel
        )
        invisible(NULL)
      }
    )
  )

  testServer(
    server_under_test,
    args = list(
      workflow_data = workflow_data,
      omic_type = "proteomics",
      experiment_label = "DIA Experiment"
    ),
    {
      session$setInputs(replicate_group_column = "replicates")
      session$flushReact()
      session$setInputs(apply_replicate_filter = 1)
      session$flushReact()

      expect_false(is.null(captured$apply))
      expect_identical(captured$apply$workflowData, workflow_data)
      expect_identical(captured$apply$replicateGroupColumn, "replicates")
      expect_identical(captured$apply$output, output)
      expect_identical(captured$apply$omicType, "proteomics")
      expect_identical(captured$apply$experimentLabel, "DIA Experiment")
      expect_type(captured$apply$replicatePlot, "closure")
    }
  )
})

test_that("mod_prot_qc_peptide_protein_server wires the revert observer through the helper seam", {
  workflow_data <- shiny::reactiveValues(state_manager = list(kind = "state_manager"))
  captured <- new.env(parent = emptyenv())
  captured$revert <- NULL

  server_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_peptide_protein_server,
    list(
      runProteinPeptideRevertObserver = function(workflowData, output) {
        captured$revert <<- list(
          workflowData = workflowData,
          output = output
        )
        invisible(NULL)
      }
    )
  )

  testServer(
    server_under_test,
    args = list(
      workflow_data = workflow_data,
      omic_type = "proteomics",
      experiment_label = "DIA Experiment"
    ),
    {
      session$setInputs(revert_protein_peptide = 1)
      session$flushReact()

      expect_false(is.null(captured$revert))
      expect_identical(captured$revert$workflowData, workflow_data)
      expect_identical(captured$revert$output, output)
    }
  )
})

test_that("mod_prot_qc_peptide_replicate_server wires the revert observer through the helper seam", {
  workflow_data <- shiny::reactiveValues(state_manager = list(kind = "state_manager"))
  captured <- new.env(parent = emptyenv())
  captured$revert <- NULL

  server_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_peptide_replicate_server,
    list(
      runPeptideReplicateRevertObserver = function(workflowData, output) {
        captured$revert <<- list(
          workflowData = workflowData,
          output = output
        )
        invisible(NULL)
      }
    )
  )

  testServer(
    server_under_test,
    args = list(
      workflow_data = workflow_data,
      omic_type = "proteomics",
      experiment_label = "DIA Experiment"
    ),
    {
      session$setInputs(revert_replicate = 1)
      session$flushReact()

      expect_false(is.null(captured$revert))
      expect_identical(captured$revert$workflowData, workflow_data)
      expect_identical(captured$revert$output, output)
    }
  )
})

test_that("runPeptideRollupRevertStep reverts to the previous state in history", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$getHistory <- function() {
    c("raw_data_s4", "qvalue_filtered", "precursor_rollup")
  }

  reverted_to <- NULL
  captured <- new.env(parent = emptyenv())
  workflow_data$state_manager$revertToState <- function(state) {
    reverted_to <<- state
    "reverted_object"
  }

  result <- runPeptideRollupRevertStep(
    workflowData = workflow_data,
    logInfoFn = function(message) {
      captured$info <- message
    }
  )

  expect_identical(result$previousState, "qvalue_filtered")
  expect_identical(result$revertedS4, "reverted_object")
  expect_identical(result$resultText, "Reverted to previous state: qvalue_filtered")
  expect_identical(reverted_to, "qvalue_filtered")
  expect_identical(captured$info, "Reverted precursor rollup to qvalue_filtered")
})

test_that("runPeptideRollupRevertStep errors when there is no prior state", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$getHistory <- function() {
    "precursor_rollup"
  }

  expect_error(
    runPeptideRollupRevertStep(workflowData = workflow_data),
    "No previous state to revert to.",
    fixed = TRUE
  )
})

test_that("runPeptideRollupApplyStep saves the rolled-up state and QC parameters", {
  if (!methods::isClass("FakePeptideRollupState")) {
    methods::setClass("FakePeptideRollupState", slots = c(peptide_data = "data.frame"))
  }

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  current_state <- methods::new(
    "FakePeptideRollupState",
    peptide_data = data.frame(Protein.Ids = c("P1", "P1", "P2"))
  )
  rolled_up_state <- methods::new(
    "FakePeptideRollupState",
    peptide_data = data.frame(Protein.Ids = c("P1", "P2", "P2"))
  )
  workflow_data$state_manager$getState <- function() {
    current_state
  }

  captured <- new.env(parent = emptyenv())
  captured$saveState <- NULL
  workflow_data$state_manager$saveState <- function(state_name, s4_data_object, config_object, description) {
    captured$saveState <- list(
      state_name = state_name,
      s4_data_object = s4_data_object,
      config_object = config_object,
      description = description
    )
    invisible(NULL)
  }

  timestamp <- as.POSIXct("2026-04-16 12:00:00", tz = "UTC")
  result <- runPeptideRollupApplyStep(
    workflowData = workflow_data,
    rollupFn = function(theObject) {
      expect_identical(theObject, current_state)
      rolled_up_state
    },
    logInfoFn = function(message) {
      captured$info <- message
    },
    nowFn = function() timestamp
  )

  expect_identical(captured$info, "QC Step: Applying precursor rollup")
  expect_true(workflow_data$qc_params$peptide_qc$precursor_rollup$applied)
  expect_identical(workflow_data$qc_params$peptide_qc$precursor_rollup$timestamp, timestamp)
  expect_identical(captured$saveState$state_name, "precursor_rollup")
  expect_identical(captured$saveState$s4_data_object, rolled_up_state)
  expect_identical(captured$saveState$config_object, list())
  expect_identical(
    captured$saveState$description,
    "Applied precursor to peptide rollup"
  )
  expect_identical(result$rolledUpS4, rolled_up_state)
  expect_match(result$resultText, "Proteins remaining: 2", fixed = TRUE)
  expect_match(result$resultText, "State saved as: 'precursor_rollup'", fixed = TRUE)
})

test_that("updatePeptideRollupOutputs refreshes result text and plot grid", {
  if (!methods::isClass("FakePeptideRollupOutputState")) {
    methods::setClass("FakePeptideRollupOutputState", slots = c(peptide_data = "data.frame"))
  }

  rollup_state <- methods::new(
    "FakePeptideRollupOutputState",
    peptide_data = data.frame(Protein.Ids = c("P1", "P2"))
  )
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$args <- NULL
  captured$plot <- NULL
  rollup_plot <- function(value) {
    captured$plot <- value
    invisible(value)
  }

  plot_grid <- updatePeptideRollupOutputs(
    output = output,
    rollupPlot = rollup_plot,
    rollupResult = list(
      rolledUpS4 = rollup_state,
      resultText = "Rollup completed"
    ),
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    renderTextFn = function(value) value,
    updateProteinFilteringFn = function(data, step_name, omic_type, experiment_label, return_grid, overwrite) {
      captured$args <- list(
        data = data,
        step_name = step_name,
        omic_type = omic_type,
        experiment_label = experiment_label,
        return_grid = return_grid,
        overwrite = overwrite
      )
      "plot_grid"
    }
  )

  expect_identical(output$rollup_results, "Rollup completed")
  expect_identical(plot_grid, "plot_grid")
  expect_identical(captured$plot, "plot_grid")
  expect_identical(captured$args$data, rollup_state@peptide_data)
  expect_identical(captured$args$step_name, "3_precursor_rollup")
  expect_identical(captured$args$omic_type, "proteomics")
  expect_identical(captured$args$experiment_label, "DIA Experiment")
  expect_true(captured$args$return_grid)
  expect_true(captured$args$overwrite)
})

test_that("runPeptideRollupApplyObserver delegates the apply workflow and notifications", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data_ref <- workflow_data
  output <- new.env(parent = emptyenv())
  output_ref <- output
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()
  captured$removed <- character()
  rollup_plot <- function(...) "plot_reactive"

  completed <- runPeptideRollupApplyObserver(
    workflowData = workflow_data,
    output = output,
    rollupPlot = rollup_plot,
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    runApplyStepFn = function(workflowData) {
      captured$calls <- c(captured$calls, "run")
      expect_identical(workflowData, workflow_data_ref)
      list(
        rolledUpS4 = "rolled_up_s4",
        resultText = "Rollup completed"
      )
    },
    updateOutputsFn = function(output, rollupPlot, rollupResult, omicType, experimentLabel) {
      captured$calls <- c(captured$calls, "update")
      expect_identical(output, output_ref)
      expect_identical(rollupPlot, rollup_plot)
      expect_identical(rollupResult$rolledUpS4, "rolled_up_s4")
      expect_identical(rollupResult$resultText, "Rollup completed")
      expect_identical(omicType, "proteomics")
      expect_identical(experimentLabel, "DIA Experiment")
      "plot_grid"
    },
    showNotificationFn = function(...) {
      captured$calls <- c(captured$calls, "show")
      captured$notifications[[length(captured$notifications) + 1]] <- list(...)
    },
    removeNotificationFn = function(id) {
      captured$calls <- c(captured$calls, "remove")
      captured$removed <- c(captured$removed, id)
    },
    logInfoFn = function(message) {
      captured$calls <- c(captured$calls, "log_info")
      captured$info <- message
    },
    logErrorFn = function(...) {
      stop("error logging should not run on the happy path")
    }
  )

  expect_identical(
    captured$calls,
    c("show", "run", "update", "log_info", "remove", "show")
  )
  expect_identical(captured$notifications[[1]][[1]], "Applying precursor rollup...")
  expect_identical(captured$notifications[[1]]$id, "rollup_working")
  expect_null(captured$notifications[[1]]$duration)
  expect_identical(
    captured$notifications[[2]][[1]],
    "Precursor rollup applied successfully"
  )
  expect_identical(captured$notifications[[2]]$type, "message")
  expect_identical(captured$removed, "rollup_working")
  expect_identical(captured$info, "Precursor rollup applied successfully")
  expect_identical(completed$status, "success")
  expect_identical(completed$rollupResult$rolledUpS4, "rolled_up_s4")
  expect_identical(completed$plotGrid, "plot_grid")
})

test_that("runPeptideRollupApplyObserver reports apply errors and clears the working notification", {
  workflow_data <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()
  captured$removed <- character()
  rollup_plot <- function(...) "plot_reactive"

  completed <- runPeptideRollupApplyObserver(
    workflowData = workflow_data,
    output = output,
    rollupPlot = rollup_plot,
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    runApplyStepFn = function(...) {
      captured$calls <- c(captured$calls, "run")
      stop("mock failure")
    },
    updateOutputsFn = function(...) {
      stop("output refresh should not run after a rollup error")
    },
    showNotificationFn = function(...) {
      captured$calls <- c(captured$calls, "show")
      captured$notifications[[length(captured$notifications) + 1]] <- list(...)
    },
    removeNotificationFn = function(id) {
      captured$calls <- c(captured$calls, "remove")
      captured$removed <- c(captured$removed, id)
    },
    logInfoFn = function(...) {
      stop("success logging should not run on the error path")
    },
    logErrorFn = function(message) {
      captured$calls <- c(captured$calls, "log_error")
      captured$error <- message
    }
  )

  expect_identical(
    captured$calls,
    c("show", "run", "log_error", "show", "remove")
  )
  expect_identical(captured$notifications[[1]][[1]], "Applying precursor rollup...")
  expect_identical(captured$notifications[[1]]$id, "rollup_working")
  expect_identical(
    captured$notifications[[2]][[1]],
    "Error applying precursor rollup: mock failure"
  )
  expect_identical(captured$notifications[[2]]$type, "error")
  expect_identical(captured$notifications[[2]]$duration, 15)
  expect_identical(captured$removed, "rollup_working")
  expect_identical(captured$error, "Error applying precursor rollup: mock failure")
  expect_identical(completed$status, "error")
  expect_identical(
    completed$errorMessage,
    "Error applying precursor rollup: mock failure"
  )
})

test_that("runPeptideRollupRevertObserver delegates the revert workflow and success notification", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data_ref <- workflow_data
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()

  completed <- runPeptideRollupRevertObserver(
    workflowData = workflow_data,
    output = output,
    runRevertStepFn = function(workflowData) {
      captured$calls <- c(captured$calls, "run")
      expect_identical(workflowData, workflow_data_ref)
      list(
        previousState = "qvalue_filtered",
        resultText = "Reverted to previous state: qvalue_filtered"
      )
    },
    renderTextFn = function(value) {
      captured$calls <- c(captured$calls, "render")
      value
    },
    showNotificationFn = function(...) {
      captured$calls <- c(captured$calls, "show")
      captured$notifications[[length(captured$notifications) + 1]] <- list(...)
    },
    logErrorFn = function(...) {
      stop("error logging should not run on the happy path")
    }
  )

  expect_identical(captured$calls, c("run", "render", "show"))
  expect_identical(
    output$rollup_results,
    "Reverted to previous state: qvalue_filtered"
  )
  expect_identical(captured$notifications[[1]][[1]], "Reverted successfully")
  expect_identical(captured$notifications[[1]]$type, "message")
  expect_identical(completed$status, "success")
  expect_identical(
    completed$revertResult$resultText,
    "Reverted to previous state: qvalue_filtered"
  )
})

test_that("runPeptideRollupRevertObserver reports revert errors without rendering results", {
  workflow_data <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()

  completed <- runPeptideRollupRevertObserver(
    workflowData = workflow_data,
    output = output,
    runRevertStepFn = function(...) {
      captured$calls <- c(captured$calls, "run")
      stop("mock revert failure")
    },
    renderTextFn = function(...) {
      stop("renderText should not run on the error path")
    },
    showNotificationFn = function(...) {
      captured$calls <- c(captured$calls, "show")
      captured$notifications[[length(captured$notifications) + 1]] <- list(...)
    },
    logErrorFn = function(message) {
      captured$calls <- c(captured$calls, "log_error")
      captured$error <- message
    }
  )

  expect_identical(captured$calls, c("run", "log_error", "show"))
  expect_identical(captured$error, "Error reverting: mock revert failure")
  expect_identical(
    captured$notifications[[1]][[1]],
    "Error reverting: mock revert failure"
  )
  expect_identical(captured$notifications[[1]]$type, "error")
  expect_false(exists("rollup_results", envir = output, inherits = FALSE))
  expect_identical(completed$status, "error")
  expect_identical(completed$errorMessage, "Error reverting: mock revert failure")
})

test_that("mod_prot_qc_peptide_rollup_server wires the apply observer through the helper seam", {
  workflow_data <- shiny::reactiveValues(state_manager = list(kind = "state_manager"))
  captured <- new.env(parent = emptyenv())
  captured$apply <- NULL

  server_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_peptide_rollup_server,
    list(
      runPeptideRollupApplyObserver = function(workflowData,
                                               output,
                                               rollupPlot,
                                               omicType,
                                               experimentLabel) {
        captured$apply <<- list(
          workflowData = workflowData,
          output = output,
          rollupPlot = rollupPlot,
          omicType = omicType,
          experimentLabel = experimentLabel
        )
        invisible(NULL)
      }
    )
  )

  testServer(
    server_under_test,
    args = list(
      workflow_data = workflow_data,
      omic_type = "proteomics",
      experiment_label = "DIA Experiment"
    ),
    {
      session$setInputs(apply_rollup = 1)
      session$flushReact()

      expect_false(is.null(captured$apply))
      expect_identical(captured$apply$workflowData, workflow_data)
      expect_identical(captured$apply$output, output)
      expect_identical(captured$apply$omicType, "proteomics")
      expect_identical(captured$apply$experimentLabel, "DIA Experiment")
      expect_type(captured$apply$rollupPlot, "closure")
    }
  )
})

test_that("mod_prot_qc_peptide_rollup_server wires the revert observer through the helper seam", {
  workflow_data <- shiny::reactiveValues(state_manager = list(kind = "state_manager"))
  captured <- new.env(parent = emptyenv())
  captured$revert <- NULL

  server_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_peptide_rollup_server,
    list(
      runPeptideRollupRevertObserver = function(workflowData, output) {
        captured$revert <<- list(
          workflowData = workflowData,
          output = output
        )
        invisible(NULL)
      }
    )
  )

  testServer(
    server_under_test,
    args = list(
      workflow_data = workflow_data,
      omic_type = "proteomics",
      experiment_label = "DIA Experiment"
    ),
    {
      session$setInputs(revert_rollup = 1)
      session$flushReact()

      expect_false(is.null(captured$revert))
      expect_identical(captured$revert$workflowData, workflow_data)
      expect_identical(captured$revert$output, output)
    }
  )
})

test_that("runPeptideSampleApplyStep saves the filtered state, removed samples, and QC parameters", {
  if (!methods::isClass("FakePeptideSampleState")) {
    methods::setClass(
      "FakePeptideSampleState",
      slots = c(args = "list", peptide_data = "data.frame", sample_id = "character")
    )
  }

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  current_state <- methods::new(
    "FakePeptideSampleState",
    args = list(filterMinNumPeptidesPerSample = list()),
    peptide_data = data.frame(
      Protein.Ids = c("P1", "P1", "P2", "P3"),
      Run = c("S1", "S1", "S2", "S3")
    ),
    sample_id = "Run"
  )
  workflow_data$state_manager$getState <- function() {
    current_state
  }

  captured <- new.env(parent = emptyenv())
  captured$configUpdate <- NULL
  captured$saveState <- NULL
  workflow_data$state_manager$saveState <- function(...) {
    captured$saveState <- list(...)
  }

  filtered_s4 <- methods::new(
    "FakePeptideSampleState",
    args = list(filterMinNumPeptidesPerSample = list()),
    peptide_data = data.frame(
      Protein.Ids = c("P1", "P1", "P3"),
      Run = c("S1", "S1", "S3")
    ),
    sample_id = "Run"
  )
  timestamp <- as.POSIXct("2026-04-16 16:05:00", tz = "UTC")

  result <- runPeptideSampleApplyStep(
    workflowData = workflow_data,
    minPeptidesPerSample = 500,
    updateConfigParameterFn = function(theObject, function_name, parameter_name, new_value) {
      captured$configUpdate <<- list(
        theObject = theObject,
        function_name = function_name,
        parameter_name = parameter_name,
        new_value = new_value
      )
      theObject
    },
    filterMinNumPeptidesPerSampleFn = function(theObject) {
      captured$filterInput <<- theObject
      filtered_s4
    },
    logInfoFn = function(message) {
      captured$info <- message
    },
    nowFn = function() timestamp
  )

  expect_identical(
    captured$info,
    "QC Step: Applying sample quality filter (min: 500)"
  )
  expect_identical(captured$configUpdate$theObject, current_state)
  expect_identical(
    captured$configUpdate$function_name,
    "filterMinNumPeptidesPerSample"
  )
  expect_identical(
    captured$configUpdate$parameter_name,
    "peptides_per_sample_cutoff"
  )
  expect_identical(captured$configUpdate$new_value, 500)
  expect_identical(captured$filterInput, current_state)
  expect_identical(
    result$filteredS4@args$filterMinNumPeptidesPerSample$samples_removed,
    "S2"
  )
  expect_identical(
    result$filteredS4@args$filterMinNumPeptidesPerSample$samples_removed_count,
    1L
  )
  expect_identical(
    result$filteredS4@args$filterMinNumPeptidesPerSample$samples_before_count,
    3L
  )
  expect_identical(
    result$filteredS4@args$filterMinNumPeptidesPerSample$samples_after_count,
    2L
  )
  expect_identical(
    workflow_data$qc_params$peptide_qc$sample_filter$min_peptides_per_sample,
    500
  )
  expect_identical(
    workflow_data$qc_params$peptide_qc$sample_filter$samples_removed,
    "S2"
  )
  expect_identical(
    workflow_data$qc_params$peptide_qc$sample_filter$samples_removed_count,
    1L
  )
  expect_identical(
    workflow_data$qc_params$peptide_qc$sample_filter$samples_before_count,
    3L
  )
  expect_identical(
    workflow_data$qc_params$peptide_qc$sample_filter$samples_after_count,
    2L
  )
  expect_identical(
    workflow_data$qc_params$peptide_qc$sample_filter$timestamp,
    timestamp
  )
  expect_identical(captured$saveState$state_name, "sample_filtered")
  expect_identical(captured$saveState$s4_data_object, result$filteredS4)
  expect_identical(captured$saveState$config_object$min_peptides_per_sample, 500)
  expect_identical(captured$saveState$config_object$samples_removed, "S2")
  expect_identical(captured$saveState$config_object$samples_removed_count, 1L)
  expect_identical(
    captured$saveState$description,
    "Applied minimum peptides per sample filter"
  )
  expect_match(result$resultText, "Proteins remaining: 2", fixed = TRUE)
  expect_match(result$resultText, "Samples remaining: 2", fixed = TRUE)
  expect_match(result$resultText, "Samples removed: 1", fixed = TRUE)
  expect_match(result$resultText, "Min peptides per sample: 500", fixed = TRUE)
  expect_match(result$resultText, "Removed samples:\nS2", fixed = TRUE)
})

test_that("updatePeptideSampleOutputs refreshes result text and plot grid", {
  if (!methods::isClass("FakePeptideSampleOutputState")) {
    methods::setClass(
      "FakePeptideSampleOutputState",
      slots = c(peptide_data = "data.frame")
    )
  }

  sample_state <- methods::new(
    "FakePeptideSampleOutputState",
    peptide_data = data.frame(Protein.Ids = c("P1", "P2"), Run = c("S1", "S3"))
  )
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$args <- NULL
  captured$plot <- NULL
  sample_plot <- function(value) {
    captured$plot <- value
    invisible(value)
  }

  plot_grid <- updatePeptideSampleOutputs(
    output = output,
    samplePlot = sample_plot,
    sampleResult = list(
      filteredS4 = sample_state,
      resultText = "Sample completed"
    ),
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    renderTextFn = function(value) value,
    updateProteinFilteringFn = function(data, step_name, omic_type, experiment_label, return_grid, overwrite) {
      captured$args <- list(
        data = data,
        step_name = step_name,
        omic_type = omic_type,
        experiment_label = experiment_label,
        return_grid = return_grid,
        overwrite = overwrite
      )
      "plot_grid"
    }
  )

  expect_identical(output$sample_results, "Sample completed")
  expect_identical(plot_grid, "plot_grid")
  expect_identical(captured$plot, "plot_grid")
  expect_identical(captured$args$data, sample_state@peptide_data)
  expect_identical(captured$args$step_name, "6_sample_filtered")
  expect_identical(captured$args$omic_type, "proteomics")
  expect_identical(captured$args$experiment_label, "DIA Experiment")
  expect_true(captured$args$return_grid)
  expect_true(captured$args$overwrite)
})

test_that("runPeptideSampleApplyObserver delegates the apply workflow and notifications", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data_ref <- workflow_data
  output <- new.env(parent = emptyenv())
  output_ref <- output
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()
  captured$removed <- character()
  sample_plot <- function(...) "plot_reactive"

  completed <- runPeptideSampleApplyObserver(
    workflowData = workflow_data,
    minPeptidesPerSample = 500,
    output = output,
    samplePlot = sample_plot,
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    runApplyStepFn = function(workflowData, minPeptidesPerSample) {
      captured$calls <- c(captured$calls, "run")
      expect_identical(workflowData, workflow_data_ref)
      expect_identical(minPeptidesPerSample, 500)
      list(
        filteredS4 = "filtered_s4",
        resultText = "Sample completed"
      )
    },
    updateOutputsFn = function(output, samplePlot, sampleResult, omicType, experimentLabel) {
      captured$calls <- c(captured$calls, "update")
      expect_identical(output, output_ref)
      expect_identical(samplePlot, sample_plot)
      expect_identical(sampleResult$filteredS4, "filtered_s4")
      expect_identical(sampleResult$resultText, "Sample completed")
      expect_identical(omicType, "proteomics")
      expect_identical(experimentLabel, "DIA Experiment")
      "plot_grid"
    },
    showNotificationFn = function(...) {
      captured$calls <- c(captured$calls, "show")
      captured$notifications[[length(captured$notifications) + 1]] <- list(...)
    },
    removeNotificationFn = function(id) {
      captured$calls <- c(captured$calls, "remove")
      captured$removed <- c(captured$removed, id)
    },
    logInfoFn = function(message) {
      captured$calls <- c(captured$calls, "log_info")
      captured$info <- message
    },
    logErrorFn = function(...) {
      stop("error logging should not run on the happy path")
    }
  )

  expect_identical(
    captured$calls,
    c("show", "run", "update", "log_info", "remove", "show")
  )
  expect_identical(
    captured$notifications[[1]][[1]],
    "Applying sample quality filter..."
  )
  expect_identical(captured$notifications[[1]]$id, "sample_working")
  expect_null(captured$notifications[[1]]$duration)
  expect_identical(
    captured$notifications[[2]][[1]],
    "Sample quality filter applied successfully"
  )
  expect_identical(captured$notifications[[2]]$type, "message")
  expect_identical(captured$removed, "sample_working")
  expect_identical(captured$info, "Sample quality filter applied successfully")
  expect_identical(completed$status, "success")
  expect_identical(completed$sampleResult$filteredS4, "filtered_s4")
  expect_identical(completed$plotGrid, "plot_grid")
})

test_that("runPeptideSampleApplyObserver reports apply errors and clears the working notification", {
  workflow_data <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()
  captured$removed <- character()
  sample_plot <- function(...) "plot_reactive"

  completed <- runPeptideSampleApplyObserver(
    workflowData = workflow_data,
    minPeptidesPerSample = 500,
    output = output,
    samplePlot = sample_plot,
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    runApplyStepFn = function(...) {
      captured$calls <- c(captured$calls, "run")
      stop("mock failure")
    },
    updateOutputsFn = function(...) {
      stop("output refresh should not run after a sample error")
    },
    showNotificationFn = function(...) {
      captured$calls <- c(captured$calls, "show")
      captured$notifications[[length(captured$notifications) + 1]] <- list(...)
    },
    removeNotificationFn = function(id) {
      captured$calls <- c(captured$calls, "remove")
      captured$removed <- c(captured$removed, id)
    },
    logInfoFn = function(...) {
      stop("success logging should not run on the error path")
    },
    logErrorFn = function(message) {
      captured$calls <- c(captured$calls, "log_error")
      captured$error <- message
    }
  )

  expect_identical(
    captured$calls,
    c("show", "run", "log_error", "show", "remove")
  )
  expect_identical(
    captured$notifications[[1]][[1]],
    "Applying sample quality filter..."
  )
  expect_identical(captured$notifications[[1]]$id, "sample_working")
  expect_identical(
    captured$notifications[[2]][[1]],
    "Error applying sample quality filter: mock failure"
  )
  expect_identical(captured$notifications[[2]]$type, "error")
  expect_identical(captured$notifications[[2]]$duration, 15)
  expect_identical(captured$removed, "sample_working")
  expect_identical(
    captured$error,
    "Error applying sample quality filter: mock failure"
  )
  expect_identical(completed$status, "error")
  expect_identical(
    completed$errorMessage,
    "Error applying sample quality filter: mock failure"
  )
})

test_that("mod_prot_qc_peptide_sample_server wires the apply observer through the helper seam", {
  workflow_data <- shiny::reactiveValues(state_manager = list(kind = "state_manager"))
  captured <- new.env(parent = emptyenv())
  captured$apply <- NULL

  server_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_peptide_sample_server,
    list(
      runPeptideSampleApplyObserver = function(workflowData,
                                               minPeptidesPerSample,
                                               output,
                                               samplePlot,
                                               omicType,
                                               experimentLabel) {
        captured$apply <<- list(
          workflowData = workflowData,
          minPeptidesPerSample = minPeptidesPerSample,
          output = output,
          samplePlot = samplePlot,
          omicType = omicType,
          experimentLabel = experimentLabel
        )
        invisible(NULL)
      }
    )
  )

  testServer(
    server_under_test,
    args = list(
      workflow_data = workflow_data,
      omic_type = "proteomics",
      experiment_label = "DIA Experiment"
    ),
    {
      session$setInputs(min_peptides_per_sample = 700)
      session$flushReact()
      session$setInputs(apply_sample_filter = 1)
      session$flushReact()

      expect_false(is.null(captured$apply))
      expect_identical(captured$apply$workflowData, workflow_data)
      expect_identical(captured$apply$minPeptidesPerSample, 700)
      expect_identical(captured$apply$output, output)
      expect_identical(captured$apply$omicType, "proteomics")
      expect_identical(captured$apply$experimentLabel, "DIA Experiment")
      expect_type(captured$apply$samplePlot, "closure")
    }
  )
})

test_that("runPeptideSampleRevertStep reverts to the previous state in history", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$getHistory <- function() {
    c("raw_data_s4", "sample_filtered", "replicate_filtered")
  }

  reverted_to <- NULL
  workflow_data$state_manager$revertToState <- function(state) {
    reverted_to <<- state
    "reverted_object"
  }

  result <- runPeptideSampleRevertStep(workflowData = workflow_data)

  expect_identical(result$previousState, "sample_filtered")
  expect_identical(result$revertedS4, "reverted_object")
  expect_identical(result$resultText, "Reverted to previous state: sample_filtered")
  expect_identical(reverted_to, "sample_filtered")
})

test_that("runPeptideSampleRevertStep errors when there is no prior state", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$getHistory <- function() {
    "raw_data_s4"
  }

  expect_error(
    runPeptideSampleRevertStep(workflowData = workflow_data),
    "No previous state to revert to.",
    fixed = TRUE
  )
})

test_that("runPeptideSampleRevertObserver delegates the revert workflow and success notification", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data_ref <- workflow_data
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()

  completed <- runPeptideSampleRevertObserver(
    workflowData = workflow_data,
    output = output,
    runRevertStepFn = function(workflowData) {
      captured$calls <- c(captured$calls, "run")
      expect_identical(workflowData, workflow_data_ref)
      list(
        previousState = "sample_filtered",
        resultText = "Reverted to previous state: sample_filtered"
      )
    },
    renderTextFn = function(value) {
      captured$calls <- c(captured$calls, "render")
      value
    },
    showNotificationFn = function(...) {
      captured$calls <- c(captured$calls, "show")
      captured$notifications[[length(captured$notifications) + 1]] <- list(...)
    },
    logErrorFn = function(...) {
      stop("error logging should not run on the happy path")
    }
  )

  expect_identical(captured$calls, c("run", "render", "show"))
  expect_identical(
    output$sample_results,
    "Reverted to previous state: sample_filtered"
  )
  expect_identical(captured$notifications[[1]][[1]], "Reverted successfully")
  expect_identical(captured$notifications[[1]]$type, "message")
  expect_identical(completed$status, "success")
  expect_identical(
    completed$revertResult$resultText,
    "Reverted to previous state: sample_filtered"
  )
})

test_that("runPeptideSampleRevertObserver reports revert errors without rendering results", {
  workflow_data <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()

  completed <- runPeptideSampleRevertObserver(
    workflowData = workflow_data,
    output = output,
    runRevertStepFn = function(...) {
      captured$calls <- c(captured$calls, "run")
      stop("mock revert failure")
    },
    renderTextFn = function(...) {
      stop("renderText should not run on the error path")
    },
    showNotificationFn = function(...) {
      captured$calls <- c(captured$calls, "show")
      captured$notifications[[length(captured$notifications) + 1]] <- list(...)
    },
    logErrorFn = function(message) {
      captured$calls <- c(captured$calls, "log_error")
      captured$error <- message
    }
  )

  expect_identical(captured$calls, c("run", "log_error", "show"))
  expect_identical(captured$error, "Error reverting: mock revert failure")
  expect_identical(
    captured$notifications[[1]][[1]],
    "Error reverting: mock revert failure"
  )
  expect_identical(captured$notifications[[1]]$type, "error")
  expect_false(exists("sample_results", envir = output, inherits = FALSE))
  expect_identical(completed$status, "error")
  expect_identical(completed$errorMessage, "Error reverting: mock revert failure")
})

test_that("mod_prot_qc_peptide_sample_server wires the revert observer through the helper seam", {
  workflow_data <- shiny::reactiveValues(state_manager = list(kind = "state_manager"))
  captured <- new.env(parent = emptyenv())
  captured$revert <- NULL

  server_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_peptide_sample_server,
    list(
      runPeptideSampleRevertObserver = function(workflowData, output) {
        captured$revert <<- list(
          workflowData = workflowData,
          output = output
        )
        invisible(NULL)
      }
    )
  )

  testServer(
    server_under_test,
    args = list(
      workflow_data = workflow_data,
      omic_type = "proteomics",
      experiment_label = "DIA Experiment"
    ),
    {
      session$setInputs(revert_sample = 1)
      session$flushReact()

      expect_false(is.null(captured$revert))
      expect_identical(captured$revert$workflowData, workflow_data)
      expect_identical(captured$revert$output, output)
    }
  )
})
