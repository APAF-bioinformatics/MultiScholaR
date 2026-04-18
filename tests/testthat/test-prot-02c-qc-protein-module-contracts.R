library(testthat)
library(shiny)

makeFunctionWithOverrides <- function(fun, replacements) {
  funOverride <- fun
  environment(funOverride) <- list2env(replacements, parent = environment(fun))
  funOverride
}

test_that("mod_prot_qc_protein_ui renders DIA submodule tabs in order when available", {
  observed_ids <- list()
  ui_fn_names <- vapply(getProtQcProteinModuleSpecs("DIA"), `[[`, character(1), "uiFn")
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
    mod_prot_qc_protein_ui,
    c(
      overrides,
      list(
        getProtQcProteinModuleSpecs = getProtQcProteinModuleSpecs,
        buildProtQcProteinTab = makeFunctionWithOverrides(buildProtQcProteinTab, overrides)
      )
    )
  )

  html <- htmltools::renderTags(ui_under_test("protein-qc", workflow_type = "DIA"))$html

  expect_match(html, "id=\"protein-qc-protein_filter_tabs\"", fixed = TRUE)
  expect_false(grepl("Module not loaded", html, fixed = TRUE))
  expect_identical(observed_ids$mod_prot_qc_protein_rollup_ui, "protein-qc-rollup")
  expect_identical(observed_ids$mod_prot_qc_protein_cleanup_ui, "protein-qc-cleanup")
  expect_identical(observed_ids$mod_prot_qc_protein_intensity_ui, "protein-qc-intensity_filter")
  expect_identical(observed_ids$mod_prot_qc_protein_dedup_ui, "protein-qc-duplicate_removal")
  expect_identical(observed_ids$mod_prot_qc_protein_replicate_ui, "protein-qc-replicate_filter")
})

test_that("mod_prot_qc_protein_ui omits DIA rollup for LFQ workflows", {
  observed_ids <- list()
  ui_fn_names <- vapply(getProtQcProteinModuleSpecs("LFQ"), `[[`, character(1), "uiFn")
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
    mod_prot_qc_protein_ui,
    c(
      overrides,
      list(
        getProtQcProteinModuleSpecs = getProtQcProteinModuleSpecs,
        buildProtQcProteinTab = makeFunctionWithOverrides(buildProtQcProteinTab, overrides)
      )
    )
  )

  html <- htmltools::renderTags(ui_under_test("protein-qc", workflow_type = "LFQ"))$html

  expect_false(grepl("protein-qc-rollup", html, fixed = TRUE))
  expect_null(observed_ids$mod_prot_qc_protein_rollup_ui)
  expect_identical(observed_ids$mod_prot_qc_protein_cleanup_ui, "protein-qc-cleanup")
  expect_identical(observed_ids$mod_prot_qc_protein_intensity_ui, "protein-qc-intensity_filter")
  expect_identical(observed_ids$mod_prot_qc_protein_dedup_ui, "protein-qc-duplicate_removal")
  expect_identical(observed_ids$mod_prot_qc_protein_replicate_ui, "protein-qc-replicate_filter")
})

test_that("mod_prot_qc_protein_ui keeps fallback labels when submodules are unavailable", {
  ui_fn_names <- vapply(getProtQcProteinModuleSpecs("DIA"), `[[`, character(1), "uiFn")
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
    mod_prot_qc_protein_ui,
    c(
      overrides,
      list(
        getProtQcProteinModuleSpecs = getProtQcProteinModuleSpecs,
        buildProtQcProteinTab = makeFunctionWithOverrides(buildProtQcProteinTab, overrides)
      )
    )
  )

  html <- htmltools::renderTags(ui_under_test("protein-qc", workflow_type = "DIA"))$html
  fallback_count <- lengths(regmatches(html, gregexpr("Module not loaded", html, fixed = TRUE)))

  expect_identical(fallback_count, 5L)
  expect_match(html, "IQ Protein Rollup", fixed = TRUE)
  expect_match(html, "Accession Cleanup", fixed = TRUE)
  expect_match(html, "Protein Intensity Filter", fixed = TRUE)
  expect_match(html, "Duplicate Removal", fixed = TRUE)
  expect_match(html, "Protein Replicate Filter", fixed = TRUE)
})

test_that("mod_prot_qc_protein_server forwards DIA orchestration to available submodules in order", {
  call_log <- character()
  call_args <- list()
  server_fn_names <- vapply(getProtQcProteinModuleSpecs("DIA"), `[[`, character(1), "serverFn")
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
        if (x %in% c("mod_prot_qc_protein_rollup_server", "mod_prot_qc_protein_replicate_server")) {
          function(id, workflow_data, experiment_paths, omic_type, experiment_label) {
            call_log <<- c(call_log, x)
            call_args[[x]] <<- list(
              id = id,
              workflowData = workflow_data,
              experimentPaths = experiment_paths,
              omicType = omic_type,
              experimentLabel = experiment_label
            )
            invisible(NULL)
          }
        } else {
          function(id, workflow_data, omic_type, experiment_label) {
            call_log <<- c(call_log, x)
            call_args[[x]] <<- list(
              id = id,
              workflowData = workflow_data,
              omicType = omic_type,
              experimentLabel = experiment_label
            )
            invisible(NULL)
          }
        }
      } else {
        base::get(x, pos = pos, envir = envir, mode = mode, inherits = inherits)
      }
    }
  )
  server_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_protein_server,
    c(
      overrides,
      list(
        getProtQcProteinModuleSpecs = getProtQcProteinModuleSpecs,
        runProtQcProteinSubmodule = makeFunctionWithOverrides(runProtQcProteinSubmodule, overrides)
      )
    )
  )

  testServer(
    server_under_test,
    args = list(
      workflow_data = workflow_data,
      experiment_paths = list(results_dir = tempdir(), protein_qc_dir = tempdir()),
      omic_type = "proteomics",
      experiment_label = "DIA Experiment"
    ),
    {
      expect_true(TRUE)
    }
  )

  expect_identical(call_log, unname(server_fn_names))
  expect_identical(call_args$mod_prot_qc_protein_rollup_server$id, "rollup")
  expect_identical(call_args$mod_prot_qc_protein_cleanup_server$id, "cleanup")
  expect_identical(call_args$mod_prot_qc_protein_intensity_server$id, "intensity_filter")
  expect_identical(call_args$mod_prot_qc_protein_dedup_server$id, "duplicate_removal")
  expect_identical(call_args$mod_prot_qc_protein_replicate_server$id, "replicate_filter")
  expect_identical(call_args$mod_prot_qc_protein_rollup_server$workflowData, workflow_data)
  expect_identical(call_args$mod_prot_qc_protein_rollup_server$experimentPaths$protein_qc_dir, tempdir())
  expect_identical(call_args$mod_prot_qc_protein_cleanup_server$workflowData, workflow_data)
  expect_identical(call_args$mod_prot_qc_protein_cleanup_server$omicType, "proteomics")
  expect_identical(call_args$mod_prot_qc_protein_cleanup_server$experimentLabel, "DIA Experiment")
  expect_identical(call_args$mod_prot_qc_protein_replicate_server$experimentPaths$protein_qc_dir, tempdir())
})

test_that("mod_prot_qc_protein_server skips rollup for non-DIA workflows and keeps common fan-out", {
  call_log <- character()
  server_fn_names <- vapply(getProtQcProteinModuleSpecs("LFQ"), `[[`, character(1), "serverFn")
  workflow_data <- shiny::reactiveValues(state_manager = list(workflow_type = "LFQ"))
  overrides <- list(
    exists = function(x, where = -1, inherits = TRUE) {
      if (x %in% c(server_fn_names, "mod_prot_qc_protein_rollup_server")) {
        TRUE
      } else {
        base::exists(x, where = where, inherits = inherits)
      }
    },
    get = function(x, pos = -1, envir = as.environment(pos), mode = "any", inherits = TRUE) {
      if (x %in% c(server_fn_names, "mod_prot_qc_protein_rollup_server")) {
        if (x %in% c("mod_prot_qc_protein_rollup_server", "mod_prot_qc_protein_replicate_server")) {
          function(...) {
            call_log <<- c(call_log, x)
            invisible(NULL)
          }
        } else {
          function(...) {
            call_log <<- c(call_log, x)
            invisible(NULL)
          }
        }
      } else {
        base::get(x, pos = pos, envir = envir, mode = mode, inherits = inherits)
      }
    }
  )
  server_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_protein_server,
    c(
      overrides,
      list(
        getProtQcProteinModuleSpecs = getProtQcProteinModuleSpecs,
        runProtQcProteinSubmodule = makeFunctionWithOverrides(runProtQcProteinSubmodule, overrides)
      )
    )
  )

  testServer(
    server_under_test,
    args = list(
      workflow_data = workflow_data,
      experiment_paths = list(results_dir = tempdir(), protein_qc_dir = tempdir()),
      omic_type = "proteomics",
      experiment_label = "LFQ Experiment"
    ),
    {
      expect_true(TRUE)
    }
  )

  expect_identical(call_log, unname(server_fn_names))
  expect_false("mod_prot_qc_protein_rollup_server" %in% call_log)
})

test_that("runProteinAccessionCleanupStep applies cleanup, tracks workflow state, and saves results", {
  if (!methods::isClass("FakeProteinCleanupState")) {
    methods::setClass(
      "FakeProteinCleanupState",
      slots = c(protein_quant_table = "data.frame")
    )
  }

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$fasta_metadata <- list(
    has_protein_evidence = TRUE,
    has_gene_names = TRUE
  )
  workflow_data$state_manager$getState <- function() {
    methods::new(
      "FakeProteinCleanupState",
      protein_quant_table = data.frame(Protein.Ids = c("P1", "P1", "P2"))
    )
  }

  captured <- new.env(parent = emptyenv())
  captured$info <- character()
  captured$saveState <- NULL
  workflow_data$state_manager$saveState <- function(...) {
    captured$saveState <- list(...)
  }

  cleaned_s4 <- methods::new(
    "FakeProteinCleanupState",
    protein_quant_table = data.frame(Protein.Ids = c("P1", "P1"))
  )
  timestamp <- as.POSIXct("2026-04-16 16:30:00", tz = "UTC")
  source_dir <- tempdir()

  result <- runProteinAccessionCleanupStep(
    workflowData = workflow_data,
    delimiter = ";",
    aggregationMethod = "median",
    chooseBestProteinAccessionFn = function(theObject,
                                            delim,
                                            seqinr_obj,
                                            seqinr_accession_column,
                                            replace_zero_with_na,
                                            aggregation_method) {
      captured$chooseArgs <- list(
        theObject = theObject,
        delim = delim,
        seqinr_obj = seqinr_obj,
        seqinr_accession_column = seqinr_accession_column,
        replace_zero_with_na = replace_zero_with_na,
        aggregation_method = aggregation_method
      )
      cleaned_s4
    },
    nowFn = function() timestamp,
    logInfoFn = function(message) {
      captured$info <- c(captured$info, message)
    },
    logWarnFn = function(message) {
      captured$warn <- message
    },
    saveRdsFn = function(object, path) {
      captured$saved <- list(object = object, path = path)
    },
    existsFn = function(x, envir = parent.frame(), inherits = TRUE) {
      x %in% c("aa_seq_tbl_final", "experiment_paths")
    },
    getFn = function(x, envir = parent.frame(), inherits = TRUE) {
      switch(
        x,
        aa_seq_tbl_final = data.frame(
          database_id = c("P1", "P2"),
          sequence = c("AAA", "BBB")
        ),
        experiment_paths = list(source_dir = source_dir)
      )
    }
  )

  expect_identical(
    captured$info[[1]],
    "Protein Processing: Applying accession cleanup with delimiter: ;"
  )
  expect_identical(captured$chooseArgs$theObject@protein_quant_table$Protein.Ids, c("P1", "P1", "P2"))
  expect_identical(captured$chooseArgs$delim, ";")
  expect_identical(captured$chooseArgs$seqinr_accession_column, "uniprot_acc")
  expect_true(captured$chooseArgs$replace_zero_with_na)
  expect_identical(captured$chooseArgs$aggregation_method, "median")
  expect_true("uniprot_acc" %in% names(captured$chooseArgs$seqinr_obj))
  expect_null(captured$warn)
  expect_identical(workflow_data$accession_cleanup_results$cleanup_applied, TRUE)
  expect_identical(workflow_data$accession_cleanup_results$delimiter_used, ";")
  expect_identical(workflow_data$accession_cleanup_results$aggregation_method, "median")
  expect_identical(workflow_data$accession_cleanup_results$proteins_before, 2L)
  expect_identical(workflow_data$accession_cleanup_results$proteins_after, 1L)
  expect_identical(workflow_data$accession_cleanup_results$had_full_metadata, TRUE)
  expect_identical(workflow_data$accession_cleanup_results$timestamp, timestamp)
  expect_identical(workflow_data$qc_params$protein_qc$accession_cleanup, workflow_data$accession_cleanup_results)
  expect_identical(
    captured$saved$path,
    file.path(source_dir, "accession_cleanup_results.RDS")
  )
  expect_identical(captured$saved$object, workflow_data$accession_cleanup_results)
  expect_identical(captured$saveState$state_name, "protein_accession_cleaned")
  expect_identical(captured$saveState$s4_data_object, cleaned_s4)
  expect_identical(captured$saveState$config_object$delimiter, ";")
  expect_identical(captured$saveState$config_object$aggregation_method, "median")
  expect_identical(captured$saveState$config_object$cleanup_applied, TRUE)
  expect_identical(captured$saveState$description, "Applied protein accession cleanup")
  expect_identical(result$cleanedS4, cleaned_s4)
  expect_identical(result$cleanupApplied, TRUE)
  expect_match(result$resultText, "Proteins remaining: 1", fixed = TRUE)
  expect_match(result$resultText, "Aggregation method: median", fixed = TRUE)
})

test_that("updateProteinAccessionCleanupOutputs refreshes result text and plot grid", {
  if (!methods::isClass("FakeProteinCleanupState")) {
    methods::setClass(
      "FakeProteinCleanupState",
      slots = c(protein_quant_table = "data.frame")
    )
  }

  cleanup_state <- methods::new(
    "FakeProteinCleanupState",
    protein_quant_table = data.frame(Protein.Ids = c("P1", "P2"))
  )
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$plot <- NULL
  accession_cleanup_plot <- function(value) {
    captured$plot <- value
    invisible(value)
  }

  plot_grid <- updateProteinAccessionCleanupOutputs(
    output = output,
    accessionCleanupPlot = accession_cleanup_plot,
    cleanupResult = list(
      cleanedS4 = cleanup_state,
      resultText = "cleanup complete"
    ),
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    renderTextFn = function(value) value,
    updateProteinFilteringFn = function(data,
                                        step_name,
                                        omic_type,
                                        experiment_label,
                                        return_grid,
                                        overwrite) {
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

  expect_identical(output$accession_cleanup_results, "cleanup complete")
  expect_identical(plot_grid, "plot_grid")
  expect_identical(captured$plot, "plot_grid")
  expect_identical(captured$args$data, cleanup_state@protein_quant_table)
  expect_identical(captured$args$step_name, "10_protein_accession_cleaned")
  expect_identical(captured$args$omic_type, "proteomics")
  expect_identical(captured$args$experiment_label, "DIA Experiment")
  expect_true(captured$args$return_grid)
  expect_true(captured$args$overwrite)
})

test_that("runProteinAccessionCleanupApplyObserver delegates apply workflow and notifications", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data_ref <- workflow_data
  output <- new.env(parent = emptyenv())
  output_ref <- output
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()
  captured$removed <- character()
  accession_cleanup_plot <- function(...) "plot_reactive"

  completed <- runProteinAccessionCleanupApplyObserver(
    workflowData = workflow_data,
    delimiter = ";",
    aggregationMethod = "mean",
    output = output,
    accessionCleanupPlot = accession_cleanup_plot,
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    runApplyStepFn = function(workflowData, delimiter, aggregationMethod) {
      captured$calls <- c(captured$calls, "run")
      expect_identical(workflowData, workflow_data_ref)
      expect_identical(delimiter, ";")
      expect_identical(aggregationMethod, "mean")
      list(
        cleanedS4 = "cleaned_s4",
        cleanupApplied = TRUE,
        resultText = "cleanup complete"
      )
    },
    updateOutputsFn = function(output, accessionCleanupPlot, cleanupResult, omicType, experimentLabel) {
      captured$calls <- c(captured$calls, "update")
      expect_identical(output, output_ref)
      expect_identical(accessionCleanupPlot, accession_cleanup_plot)
      expect_identical(cleanupResult$cleanedS4, "cleaned_s4")
      expect_identical(cleanupResult$resultText, "cleanup complete")
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
    "Applying protein accession cleanup..."
  )
  expect_identical(captured$notifications[[1]]$id, "accession_cleanup_working")
  expect_null(captured$notifications[[1]]$duration)
  expect_identical(
    captured$notifications[[2]][[1]],
    "Protein accession cleanup completed"
  )
  expect_identical(captured$notifications[[2]]$type, "message")
  expect_identical(captured$removed, "accession_cleanup_working")
  expect_identical(captured$info, "Protein accession cleanup completed")
  expect_identical(completed$status, "success")
  expect_identical(completed$cleanupResult$cleanedS4, "cleaned_s4")
  expect_identical(completed$plotGrid, "plot_grid")
})

test_that("runProteinAccessionCleanupApplyObserver reports apply errors and clears the working notification", {
  workflow_data <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()
  captured$removed <- character()
  accession_cleanup_plot <- function(...) "plot_reactive"

  completed <- runProteinAccessionCleanupApplyObserver(
    workflowData = workflow_data,
    delimiter = ";",
    aggregationMethod = "mean",
    output = output,
    accessionCleanupPlot = accession_cleanup_plot,
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    runApplyStepFn = function(...) {
      captured$calls <- c(captured$calls, "run")
      stop("mock failure")
    },
    updateOutputsFn = function(...) {
      stop("output refresh should not run after an accession cleanup error")
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
    "Applying protein accession cleanup..."
  )
  expect_identical(captured$notifications[[1]]$id, "accession_cleanup_working")
  expect_identical(
    captured$notifications[[2]][[1]],
    "Error in protein accession cleanup: mock failure"
  )
  expect_identical(captured$notifications[[2]]$type, "error")
  expect_identical(captured$notifications[[2]]$duration, 15)
  expect_identical(captured$removed, "accession_cleanup_working")
  expect_identical(captured$error, "Error in protein accession cleanup: mock failure")
  expect_identical(completed$status, "error")
  expect_identical(
    completed$errorMessage,
    "Error in protein accession cleanup: mock failure"
  )
})

test_that("runProteinAccessionCleanupRevertStep reverts to the previous state in history", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$getHistory <- function() c("imported", "protein_accession_cleaned")
  workflow_data$state_manager$revertToState <- function(state_name) {
    expect_identical(state_name, "imported")
    "reverted_s4"
  }

  result <- runProteinAccessionCleanupRevertStep(workflowData = workflow_data)

  expect_identical(result$previousState, "imported")
  expect_identical(result$revertedS4, "reverted_s4")
  expect_identical(result$resultText, "Reverted to previous state: imported")
})

test_that("runProteinAccessionCleanupRevertStep errors when there is no prior state", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$getHistory <- function() "protein_accession_cleaned"

  expect_error(
    runProteinAccessionCleanupRevertStep(workflowData = workflow_data),
    "No previous state to revert to\\."
  )
})

test_that("runProteinAccessionCleanupRevertObserver updates output and reports success", {
  workflow_data <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  output_ref <- output
  workflow_data_ref <- workflow_data
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()

  completed <- runProteinAccessionCleanupRevertObserver(
    workflowData = workflow_data,
    output = output,
    runRevertStepFn = function(workflowData) {
      captured$calls <- c(captured$calls, "run")
      expect_identical(workflowData, workflow_data_ref)

      list(
        previousState = "imported",
        revertedS4 = "reverted_s4",
        resultText = "Reverted to previous state: imported"
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
    logInfoFn = function(message) {
      captured$calls <- c(captured$calls, "log_info")
      captured$info <- message
    },
    logErrorFn = function(...) {
      stop("error logging should not run on the happy path")
    }
  )

  expect_identical(output, output_ref)
  expect_identical(
    captured$calls,
    c("run", "render", "log_info", "show")
  )
  expect_identical(output$accession_cleanup_results, "Reverted to previous state: imported")
  expect_identical(captured$info, "Reverted accession cleanup to imported")
  expect_identical(captured$notifications[[1]][[1]], "Reverted successfully")
  expect_identical(captured$notifications[[1]]$type, "message")
  expect_identical(completed$status, "success")
  expect_identical(completed$revertResult$previousState, "imported")
})

test_that("runProteinAccessionCleanupRevertObserver reports revert errors", {
  workflow_data <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()

  completed <- runProteinAccessionCleanupRevertObserver(
    workflowData = workflow_data,
    output = output,
    runRevertStepFn = function(...) {
      captured$calls <- c(captured$calls, "run")
      stop("mock revert failure")
    },
    renderTextFn = function(...) {
      stop("rendering should not run after a revert failure")
    },
    showNotificationFn = function(...) {
      captured$calls <- c(captured$calls, "show")
      captured$notifications[[length(captured$notifications) + 1]] <- list(...)
    },
    logInfoFn = function(...) {
      stop("success logging should not run on the error path")
    },
    logErrorFn = function(message) {
      captured$calls <- c(captured$calls, "log_error")
      captured$error <- message
    }
  )

  expect_identical(captured$calls, c("run", "log_error", "show"))
  expect_identical(captured$error, "Error reverting: mock revert failure")
  expect_identical(captured$notifications[[1]][[1]], "Error reverting: mock revert failure")
  expect_identical(captured$notifications[[1]]$type, "error")
  expect_identical(completed$status, "error")
  expect_identical(completed$errorMessage, "Error reverting: mock revert failure")
})

test_that("bindProteinAccessionCleanupPlot binds the rendered plot output to the stored grid", {
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$drawn <- NULL

  accession_cleanup_plot <- function(value) {
    if (missing(value)) {
      return("plot_grid")
    }

    stop("bindProteinAccessionCleanupPlot should only read from the reactive")
  }

  helper_under_test <- makeFunctionWithOverrides(
    bindProteinAccessionCleanupPlot,
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

  helper_under_test(output = output, accessionCleanupPlot = accession_cleanup_plot)

  expect_identical(output$accession_cleanup_plot, "rendered_plot")
  expect_identical(captured$drawn, "plot_grid")
})

test_that("mod_prot_qc_protein_cleanup_server wires the apply observer through the helper seam", {
  workflow_data <- shiny::reactiveValues(state_manager = list(kind = "state_manager"))
  captured <- new.env(parent = emptyenv())
  captured$apply <- NULL

  server_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_protein_cleanup_server,
    list(
      runProteinAccessionCleanupApplyObserver = function(workflowData,
                                                         delimiter,
                                                         aggregationMethod,
                                                         output,
                                                         accessionCleanupPlot,
                                                         omicType,
                                                         experimentLabel) {
        captured$apply <<- list(
          workflowData = workflowData,
          delimiter = delimiter,
          aggregationMethod = aggregationMethod,
          output = output,
          accessionCleanupPlot = accessionCleanupPlot,
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
      session$setInputs(delimiter = ";")
      session$setInputs(aggregation_method = "median")
      session$flushReact()
      session$setInputs(apply_accession_cleanup = 1)
      session$flushReact()

      expect_false(is.null(captured$apply))
      expect_identical(captured$apply$workflowData, workflow_data)
      expect_identical(captured$apply$delimiter, ";")
      expect_identical(captured$apply$aggregationMethod, "median")
      expect_identical(captured$apply$output, output)
      expect_type(captured$apply$accessionCleanupPlot, "closure")
      expect_identical(captured$apply$omicType, "proteomics")
      expect_identical(captured$apply$experimentLabel, "DIA Experiment")
    }
  )
})

test_that("mod_prot_qc_protein_cleanup_server wires the render binding through the helper seam", {
  workflow_data <- shiny::reactiveValues(state_manager = list(kind = "state_manager"))
  captured <- new.env(parent = emptyenv())
  captured$bind <- NULL

  server_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_protein_cleanup_server,
    list(
      bindProteinAccessionCleanupPlot = function(output, accessionCleanupPlot) {
        captured$bind <<- list(
          output = output,
          accessionCleanupPlot = accessionCleanupPlot
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
      expect_type(captured$bind$accessionCleanupPlot, "closure")
    }
  )
})

test_that("mod_prot_qc_protein_cleanup_server wires the revert observer through the helper seam", {
  workflow_data <- shiny::reactiveValues(state_manager = list(kind = "state_manager"))
  captured <- new.env(parent = emptyenv())
  captured$revert <- NULL

  server_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_protein_cleanup_server,
    list(
      runProteinAccessionCleanupRevertObserver = function(workflowData, output) {
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
      session$setInputs(revert_accession_cleanup = 1)
      session$flushReact()

      expect_false(is.null(captured$revert))
      expect_identical(captured$revert$workflowData, workflow_data)
      expect_identical(captured$revert$output, output)
    }
  )
})

test_that("runProteinDuplicateRemovalStep aggregates duplicates and saves state", {
  if (!methods::isClass("FakeProteinDedupState")) {
    methods::setClass(
      "FakeProteinDedupState",
      slots = c(protein_quant_table = "data.frame")
    )
  }

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$getState <- function() {
    methods::new(
      "FakeProteinDedupState",
      protein_quant_table = data.frame(
        Protein.Ids = c("P1", "P1", "P2"),
        `1` = c(1, 3, 5),
        `2` = c(2, 4, 6),
        check.names = FALSE
      )
    )
  }

  captured <- new.env(parent = emptyenv())
  captured$saveState <- NULL
  captured$aggregationMethod <- NULL
  captured$info <- NULL
  workflow_data$state_manager$saveState <- function(...) {
    captured$saveState <- list(...)
  }

  result <- runProteinDuplicateRemovalStep(
    workflowData = workflow_data,
    aggregationMethod = "mean",
    aggregationResolverFn = function(name) {
      captured$aggregationMethod <- name
      base::mean
    },
    logInfoFn = function(message) {
      captured$info <- message
    }
  )

  expect_identical(captured$aggregationMethod, "mean")
  expect_identical(
    captured$info,
    "Protein Processing: Removing duplicate proteins using mean"
  )
  expect_identical(result$duplicates, "P1")
  expect_identical(
    as.data.frame(result$deduplicatedS4@protein_quant_table),
    data.frame(
      Protein.Ids = c("P1", "P2"),
      `1` = c(2, 5),
      `2` = c(3, 6),
      check.names = FALSE
    )
  )
  expect_identical(captured$saveState$state_name, "duplicates_removed")
  expect_identical(captured$saveState$config_object$aggregation_method, "mean")
  expect_identical(captured$saveState$config_object$duplicates_found, "P1")
  expect_identical(captured$saveState$config_object$num_duplicates, 1L)
  expect_identical(
    captured$saveState$description,
    "Removed duplicate proteins by aggregation"
  )
  expect_match(result$resultText, "Proteins remaining: 2", fixed = TRUE)
  expect_match(result$resultText, "Duplicates found: 1", fixed = TRUE)
  expect_match(result$resultText, "Aggregation method: mean", fixed = TRUE)
})

test_that("updateProteinDuplicateRemovalOutputs refreshes results text and plot grid", {
  if (!methods::isClass("FakeProteinDedupState")) {
    methods::setClass(
      "FakeProteinDedupState",
      slots = c(protein_quant_table = "data.frame")
    )
  }

  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$plot <- NULL
  captured$args <- NULL

  duplicate_removal_plot <- function(value) {
    captured$plot <- value
  }

  duplicate_result <- list(
    resultText = "duplicate removal complete",
    deduplicatedS4 = methods::new(
      "FakeProteinDedupState",
      protein_quant_table = data.frame(
        Protein.Ids = c("P1", "P2"),
        `1` = c(2, 5),
        check.names = FALSE
      )
    )
  )

  plot_grid <- updateProteinDuplicateRemovalOutputs(
    output = output,
    duplicateRemovalPlot = duplicate_removal_plot,
    duplicateRemovalResult = duplicate_result,
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    renderTextFn = function(value) value,
    updateProteinFilteringFn = function(data,
                                        step_name,
                                        omic_type,
                                        experiment_label,
                                        return_grid,
                                        overwrite) {
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

  expect_identical(output$duplicate_removal_results, "duplicate removal complete")
  expect_identical(plot_grid, "plot_grid")
  expect_identical(captured$plot, "plot_grid")
  expect_identical(captured$args$data, duplicate_result$deduplicatedS4@protein_quant_table)
  expect_identical(captured$args$step_name, "12_duplicates_removed")
  expect_identical(captured$args$omic_type, "proteomics")
  expect_identical(captured$args$experiment_label, "DIA Experiment")
  expect_true(captured$args$return_grid)
  expect_true(captured$args$overwrite)
})

test_that("runProteinDuplicateRemovalApplyObserver delegates apply workflow and notifications", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data_ref <- workflow_data
  output <- new.env(parent = emptyenv())
  output_ref <- output
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()
  captured$removed <- character()
  duplicate_removal_plot <- function(...) "plot_reactive"

  completed <- runProteinDuplicateRemovalApplyObserver(
    workflowData = workflow_data,
    aggregationMethod = "median",
    output = output,
    duplicateRemovalPlot = duplicate_removal_plot,
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    runApplyStepFn = function(workflowData, aggregationMethod) {
      captured$calls <- c(captured$calls, "run")
      expect_identical(workflowData, workflow_data_ref)
      expect_identical(aggregationMethod, "median")
      list(
        deduplicatedS4 = "deduplicated_s4",
        duplicates = "P1",
        resultText = "duplicate removal complete"
      )
    },
    updateOutputsFn = function(output,
                               duplicateRemovalPlot,
                               duplicateRemovalResult,
                               omicType,
                               experimentLabel) {
      captured$calls <- c(captured$calls, "update")
      expect_identical(output, output_ref)
      expect_identical(duplicateRemovalPlot, duplicate_removal_plot)
      expect_identical(duplicateRemovalResult$deduplicatedS4, "deduplicated_s4")
      expect_identical(duplicateRemovalResult$resultText, "duplicate removal complete")
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
    "Removing duplicate proteins..."
  )
  expect_identical(captured$notifications[[1]]$id, "duplicate_removal_working")
  expect_null(captured$notifications[[1]]$duration)
  expect_identical(
    captured$notifications[[2]][[1]],
    "Duplicate protein removal completed successfully"
  )
  expect_identical(captured$notifications[[2]]$type, "message")
  expect_identical(captured$removed, "duplicate_removal_working")
  expect_identical(
    captured$info,
    "Duplicate protein removal completed successfully"
  )
  expect_identical(completed$status, "success")
  expect_identical(completed$duplicateRemovalResult$deduplicatedS4, "deduplicated_s4")
  expect_identical(completed$plotGrid, "plot_grid")
})

test_that("runProteinDuplicateRemovalApplyObserver reports apply errors and clears the working notification", {
  workflow_data <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()
  captured$removed <- character()
  duplicate_removal_plot <- function(...) "plot_reactive"

  completed <- runProteinDuplicateRemovalApplyObserver(
    workflowData = workflow_data,
    aggregationMethod = "median",
    output = output,
    duplicateRemovalPlot = duplicate_removal_plot,
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    runApplyStepFn = function(...) {
      captured$calls <- c(captured$calls, "run")
      stop("mock failure")
    },
    updateOutputsFn = function(...) {
      stop("output refresh should not run after a duplicate removal error")
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
    "Removing duplicate proteins..."
  )
  expect_identical(captured$notifications[[1]]$id, "duplicate_removal_working")
  expect_identical(
    captured$notifications[[2]][[1]],
    "Error removing duplicate proteins: mock failure"
  )
  expect_identical(captured$notifications[[2]]$type, "error")
  expect_identical(captured$notifications[[2]]$duration, 15)
  expect_identical(captured$removed, "duplicate_removal_working")
  expect_identical(
    captured$error,
    "Error removing duplicate proteins: mock failure"
  )
  expect_identical(completed$status, "error")
  expect_identical(
    completed$errorMessage,
    "Error removing duplicate proteins: mock failure"
  )
})

test_that("mod_prot_qc_protein_dedup_server wires the apply observer through the helper seam", {
  workflow_data <- shiny::reactiveValues(state_manager = list(kind = "state_manager"))
  captured <- new.env(parent = emptyenv())
  captured$apply <- NULL

  server_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_protein_dedup_server,
    list(
      runProteinDuplicateRemovalApplyObserver = function(workflowData,
                                                         aggregationMethod,
                                                         output,
                                                         duplicateRemovalPlot,
                                                         omicType,
                                                         experimentLabel) {
        captured$apply <<- list(
          workflowData = workflowData,
          aggregationMethod = aggregationMethod,
          output = output,
          duplicateRemovalPlot = duplicateRemovalPlot,
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
      session$setInputs(duplicate_aggregation_method = "median")
      session$flushReact()
      session$setInputs(apply_duplicate_removal = 1)
      session$flushReact()

      expect_false(is.null(captured$apply))
      expect_identical(captured$apply$workflowData, workflow_data)
      expect_identical(captured$apply$aggregationMethod, "median")
      expect_identical(captured$apply$output, output)
      expect_type(captured$apply$duplicateRemovalPlot, "closure")
      expect_identical(captured$apply$omicType, "proteomics")
      expect_identical(captured$apply$experimentLabel, "DIA Experiment")
    }
  )
})

test_that("bindProteinDuplicateRemovalPlot binds the rendered plot output to the stored grid", {
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$drawn <- NULL

  duplicate_removal_plot <- function(value) {
    if (missing(value)) {
      return("plot_grid")
    }

    stop("bindProteinDuplicateRemovalPlot should only read from the reactive")
  }

  helper_under_test <- makeFunctionWithOverrides(
    bindProteinDuplicateRemovalPlot,
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

  helper_under_test(output = output, duplicateRemovalPlot = duplicate_removal_plot)

  expect_identical(output$duplicate_removal_plot, "rendered_plot")
  expect_identical(captured$drawn, "plot_grid")
})

test_that("mod_prot_qc_protein_dedup_server wires the render binding through the helper seam", {
  workflow_data <- shiny::reactiveValues(state_manager = list(kind = "state_manager"))
  captured <- new.env(parent = emptyenv())
  captured$bind <- NULL

  server_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_protein_dedup_server,
    list(
      bindProteinDuplicateRemovalPlot = function(output, duplicateRemovalPlot) {
        captured$bind <<- list(
          output = output,
          duplicateRemovalPlot = duplicateRemovalPlot
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
      expect_type(captured$bind$duplicateRemovalPlot, "closure")
    }
  )
})

test_that("runProteinDuplicateRemovalRevertStep reverts to the previous state in history", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$getHistory <- function() c("imported", "duplicates_removed")
  workflow_data$state_manager$revertToState <- function(state_name) {
    expect_identical(state_name, "imported")
    "reverted_s4"
  }

  result <- runProteinDuplicateRemovalRevertStep(workflowData = workflow_data)

  expect_identical(result$previousState, "imported")
  expect_identical(result$revertedS4, "reverted_s4")
  expect_identical(result$resultText, "Reverted to previous state: imported")
})

test_that("runProteinDuplicateRemovalRevertStep errors when there is no prior state", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$getHistory <- function() "duplicates_removed"

  expect_error(
    runProteinDuplicateRemovalRevertStep(workflowData = workflow_data),
    "No previous state to revert to\\."
  )
})

test_that("runProteinDuplicateRemovalRevertObserver updates output and reports success", {
  workflow_data <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  output_ref <- output
  workflow_data_ref <- workflow_data
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()

  completed <- runProteinDuplicateRemovalRevertObserver(
    workflowData = workflow_data,
    output = output,
    runRevertStepFn = function(workflowData) {
      captured$calls <- c(captured$calls, "run")
      expect_identical(workflowData, workflow_data_ref)

      list(
        previousState = "imported",
        revertedS4 = "reverted_s4",
        resultText = "Reverted to previous state: imported"
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
    logInfoFn = function(message) {
      captured$calls <- c(captured$calls, "log_info")
      captured$info <- message
    },
    logErrorFn = function(...) {
      stop("error logging should not run on the happy path")
    }
  )

  expect_identical(output, output_ref)
  expect_identical(captured$calls, c("run", "render", "log_info", "show"))
  expect_identical(output$duplicate_removal_results, "Reverted to previous state: imported")
  expect_identical(captured$info, "Reverted duplicate removal to imported")
  expect_identical(captured$notifications[[1]][[1]], "Reverted successfully")
  expect_identical(captured$notifications[[1]]$type, "message")
  expect_identical(completed$status, "success")
  expect_identical(completed$revertResult$previousState, "imported")
})

test_that("runProteinDuplicateRemovalRevertObserver reports revert errors", {
  workflow_data <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()

  completed <- runProteinDuplicateRemovalRevertObserver(
    workflowData = workflow_data,
    output = output,
    runRevertStepFn = function(...) {
      captured$calls <- c(captured$calls, "run")
      stop("mock revert failure")
    },
    renderTextFn = function(...) {
      stop("rendering should not run after a revert failure")
    },
    showNotificationFn = function(...) {
      captured$calls <- c(captured$calls, "show")
      captured$notifications[[length(captured$notifications) + 1]] <- list(...)
    },
    logInfoFn = function(...) {
      stop("success logging should not run on the error path")
    },
    logErrorFn = function(message) {
      captured$calls <- c(captured$calls, "log_error")
      captured$error <- message
    }
  )

  expect_identical(captured$calls, c("run", "log_error", "show"))
  expect_identical(captured$error, "Error reverting: mock revert failure")
  expect_identical(captured$notifications[[1]][[1]], "Error reverting: mock revert failure")
  expect_identical(captured$notifications[[1]]$type, "error")
  expect_identical(completed$status, "error")
  expect_identical(completed$errorMessage, "Error reverting: mock revert failure")
})

test_that("mod_prot_qc_protein_dedup_server wires the revert observer through the helper seam", {
  workflow_data <- shiny::reactiveValues(state_manager = list(kind = "state_manager"))
  captured <- new.env(parent = emptyenv())
  captured$revert <- NULL

  server_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_protein_dedup_server,
    list(
      runProteinDuplicateRemovalRevertObserver = function(workflowData, output) {
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
      session$setInputs(revert_duplicate_removal = 1)
      session$flushReact()

      expect_false(is.null(captured$revert))
      expect_identical(captured$revert$workflowData, workflow_data)
      expect_identical(captured$revert$output, output)
    }
  )
})

test_that("runProteinIntensityFilterApplyStep applies flexible-mode updates and saves the filtered state", {
  if (!methods::isClass("FakeProteinIntensityState")) {
    methods::setClass(
      "FakeProteinIntensityState",
      slots = c(protein_quant_table = "data.frame", args = "list")
    )
  }

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  original_state <- methods::new(
    "FakeProteinIntensityState",
    protein_quant_table = data.frame(Protein.Ids = c("P1", "P2", "P3")),
    args = list(
      removeRowsWithMissingValuesPercent = list(
        groupwise_percentage_cutoff = 25,
        max_groups_percentage_cutoff = 50,
        proteins_intensity_cutoff_percentile = 1
      )
    )
  )
  filtered_state <- methods::new(
    "FakeProteinIntensityState",
    protein_quant_table = data.frame(Protein.Ids = c("P1", "P2")),
    args = original_state@args
  )
  workflow_data$state_manager$getState <- function() original_state

  captured <- new.env(parent = emptyenv())
  captured$saveState <- NULL
  captured$info <- NULL
  captured$missingUpdate <- NULL
  captured$filterInput <- NULL
  timestamp <- as.POSIXct("2026-04-16 12:00:00", tz = "UTC")

  workflow_data$state_manager$saveState <- function(...) {
    captured$saveState <- list(...)
  }

  result <- runProteinIntensityFilterApplyStep(
    workflowData = workflow_data,
    useStrictMode = FALSE,
    minRepsPerGroup = 3,
    minGroups = 2,
    intensityCutoffPercentile = 1.5,
    updateConfigParameterFn = function(theObject, function_name, parameter_name, new_value) {
      updated_args <- theObject@args
      updated_args$removeRowsWithMissingValuesPercent[[parameter_name]] <- new_value
      methods::new(
        "FakeProteinIntensityState",
        protein_quant_table = theObject@protein_quant_table,
        args = updated_args
      )
    },
    updateMissingValueParametersFn = function(theObject, min_reps_per_group, min_groups) {
      captured$missingUpdate <- list(
        min_reps_per_group = min_reps_per_group,
        min_groups = min_groups
      )
      updated_args <- theObject@args
      updated_args$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff <- 12.345
      updated_args$removeRowsWithMissingValuesPercent$max_groups_percentage_cutoff <- 67.89
      methods::new(
        "FakeProteinIntensityState",
        protein_quant_table = theObject@protein_quant_table,
        args = updated_args
      )
    },
    removeRowsWithMissingValuesPercentFn = function(theObject) {
      captured$filterInput <- theObject
      filtered_state
    },
    logInfoFn = function(message) {
      captured$info <- message
    },
    nowFn = function() timestamp
  )

  expect_identical(
    captured$info,
    "Protein Processing: Using FLEXIBLE MODE (adaptive thresholds)"
  )
  expect_identical(captured$missingUpdate$min_reps_per_group, 3)
  expect_identical(captured$missingUpdate$min_groups, 2)
  expect_identical(
    captured$filterInput@args$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff,
    12.345
  )
  expect_identical(
    captured$filterInput@args$removeRowsWithMissingValuesPercent$max_groups_percentage_cutoff,
    67.89
  )
  expect_equal(
    captured$filterInput@args$removeRowsWithMissingValuesPercent$proteins_intensity_cutoff_percentile,
    1.5
  )
  expect_identical(workflow_data$qc_params$protein_qc$intensity_filter$strict_mode, FALSE)
  expect_identical(workflow_data$qc_params$protein_qc$intensity_filter$min_reps_per_group, 3)
  expect_identical(workflow_data$qc_params$protein_qc$intensity_filter$min_groups, 2)
  expect_identical(
    workflow_data$qc_params$protein_qc$intensity_filter$groupwise_percentage_cutoff,
    12.345
  )
  expect_identical(
    workflow_data$qc_params$protein_qc$intensity_filter$max_groups_percentage_cutoff,
    67.89
  )
  expect_equal(
    workflow_data$qc_params$protein_qc$intensity_filter$proteins_intensity_cutoff_percentile,
    1.5
  )
  expect_identical(workflow_data$qc_params$protein_qc$intensity_filter$timestamp, timestamp)
  expect_identical(captured$saveState$state_name, "protein_intensity_filtered")
  expect_identical(captured$saveState$s4_data_object, filtered_state)
  expect_identical(captured$saveState$config_object$strict_mode, FALSE)
  expect_identical(captured$saveState$config_object$min_reps_per_group, 3)
  expect_identical(captured$saveState$config_object$min_groups, 2)
  expect_identical(captured$saveState$config_object$groupwise_percentage_cutoff, 12.345)
  expect_identical(captured$saveState$config_object$max_groups_percentage_cutoff, 67.89)
  expect_equal(captured$saveState$config_object$proteins_intensity_cutoff_percentile, 1.5)
  expect_identical(
    captured$saveState$description,
    "Applied FLEXIBLE protein intensity filter (adaptive thresholds)"
  )
  expect_identical(result$filteredS4, filtered_state)
  expect_match(result$resultText, "Mode: FLEXIBLE \\(Adaptive Thresholds\\)")
  expect_match(result$resultText, "Proteins remaining: 2", fixed = TRUE)
  expect_match(result$resultText, "Min replicates per group: 3", fixed = TRUE)
  expect_match(result$resultText, "Min groups required: 2", fixed = TRUE)
})

test_that("runProteinIntensityFilterApplyStep applies strict-mode zero cutoffs before filtering", {
  if (!methods::isClass("FakeProteinIntensityState")) {
    methods::setClass(
      "FakeProteinIntensityState",
      slots = c(protein_quant_table = "data.frame", args = "list")
    )
  }

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  original_state <- methods::new(
    "FakeProteinIntensityState",
    protein_quant_table = data.frame(Protein.Ids = c("P1", "P2", "P3")),
    args = list(
      removeRowsWithMissingValuesPercent = list(
        groupwise_percentage_cutoff = 25,
        max_groups_percentage_cutoff = 50,
        proteins_intensity_cutoff_percentile = 1
      )
    )
  )
  filtered_state <- methods::new(
    "FakeProteinIntensityState",
    protein_quant_table = data.frame(Protein.Ids = c("P1")),
    args = original_state@args
  )
  workflow_data$state_manager$getState <- function() original_state

  captured <- new.env(parent = emptyenv())
  captured$saveState <- NULL
  captured$info <- NULL
  captured$filterInput <- NULL
  workflow_data$state_manager$saveState <- function(...) {
    captured$saveState <- list(...)
  }

  result <- runProteinIntensityFilterApplyStep(
    workflowData = workflow_data,
    useStrictMode = TRUE,
    minRepsPerGroup = 3,
    minGroups = 2,
    intensityCutoffPercentile = 2.5,
    updateConfigParameterFn = function(theObject, function_name, parameter_name, new_value) {
      updated_args <- theObject@args
      updated_args$removeRowsWithMissingValuesPercent[[parameter_name]] <- new_value
      methods::new(
        "FakeProteinIntensityState",
        protein_quant_table = theObject@protein_quant_table,
        args = updated_args
      )
    },
    updateMissingValueParametersFn = function(...) {
      stop("flexible-mode updates should not run in strict mode")
    },
    removeRowsWithMissingValuesPercentFn = function(theObject) {
      captured$filterInput <- theObject
      filtered_state
    },
    logInfoFn = function(message) {
      captured$info <- message
    }
  )

  expect_identical(
    captured$info,
    "Protein Processing: Using STRICT MODE (no missing values allowed)"
  )
  expect_identical(
    captured$filterInput@args$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff,
    0
  )
  expect_identical(
    captured$filterInput@args$removeRowsWithMissingValuesPercent$max_groups_percentage_cutoff,
    0
  )
  expect_equal(
    captured$filterInput@args$removeRowsWithMissingValuesPercent$proteins_intensity_cutoff_percentile,
    2.5
  )
  expect_identical(workflow_data$qc_params$protein_qc$intensity_filter$strict_mode, TRUE)
  expect_true(is.na(workflow_data$qc_params$protein_qc$intensity_filter$min_reps_per_group))
  expect_true(is.na(workflow_data$qc_params$protein_qc$intensity_filter$min_groups))
  expect_identical(captured$saveState$state_name, "protein_intensity_filtered")
  expect_identical(
    captured$saveState$description,
    "Applied STRICT protein intensity filter (no missing values)"
  )
  expect_identical(result$filteredS4, filtered_state)
  expect_match(result$resultText, "Mode: STRICT \\(No Missing Values\\)")
  expect_match(
    result$resultText,
    "Groupwise % cutoff: 0.000% (strict - no missing allowed)",
    fixed = TRUE
  )
})

test_that("updateProteinIntensityFilterOutputs refreshes result text and plot grid", {
  if (!methods::isClass("FakeProteinIntensityState")) {
    methods::setClass(
      "FakeProteinIntensityState",
      slots = c(protein_quant_table = "data.frame", args = "list")
    )
  }

  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$plot <- NULL
  captured$args <- NULL

  protein_intensity_filter_plot <- function(value) {
    captured$plot <- value
  }

  plot_grid <- updateProteinIntensityFilterOutputs(
    output = output,
    proteinIntensityFilterPlot = protein_intensity_filter_plot,
    intensityFilterResult = list(
      resultText = "protein intensity complete",
      filteredS4 = methods::new(
        "FakeProteinIntensityState",
        protein_quant_table = data.frame(Protein.Ids = c("P1", "P2")),
        args = list(removeRowsWithMissingValuesPercent = list())
      )
    ),
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    renderTextFn = function(value) value,
    updateProteinFilteringFn = function(data,
                                        step_name,
                                        omic_type,
                                        experiment_label,
                                        return_grid,
                                        overwrite) {
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

  expect_identical(output$protein_intensity_filter_results, "protein intensity complete")
  expect_identical(plot_grid, "plot_grid")
  expect_identical(captured$plot, "plot_grid")
  expect_identical(captured$args$data$Protein.Ids, c("P1", "P2"))
  expect_identical(captured$args$step_name, "11_protein_intensity_filtered")
  expect_identical(captured$args$omic_type, "proteomics")
  expect_identical(captured$args$experiment_label, "DIA Experiment")
  expect_true(captured$args$return_grid)
  expect_true(captured$args$overwrite)
})

test_that("runProteinIntensityFilterApplyObserver delegates the apply workflow and notifications", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data_ref <- workflow_data
  output <- new.env(parent = emptyenv())
  output_ref <- output
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()
  captured$removed <- character()
  protein_intensity_filter_plot <- function(...) "plot_reactive"

  completed <- runProteinIntensityFilterApplyObserver(
    workflowData = workflow_data,
    useStrictMode = FALSE,
    minRepsPerGroup = 3,
    minGroups = 2,
    intensityCutoffPercentile = 1.5,
    output = output,
    proteinIntensityFilterPlot = protein_intensity_filter_plot,
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    runApplyStepFn = function(workflowData,
                              useStrictMode,
                              minRepsPerGroup,
                              minGroups,
                              intensityCutoffPercentile) {
      captured$calls <- c(captured$calls, "run")
      expect_identical(workflowData, workflow_data_ref)
      expect_identical(useStrictMode, FALSE)
      expect_identical(minRepsPerGroup, 3)
      expect_identical(minGroups, 2)
      expect_equal(intensityCutoffPercentile, 1.5)
      list(
        filteredS4 = "filtered_s4",
        resultText = "protein intensity complete"
      )
    },
    updateOutputsFn = function(output,
                               proteinIntensityFilterPlot,
                               intensityFilterResult,
                               omicType,
                               experimentLabel) {
      captured$calls <- c(captured$calls, "update")
      expect_identical(output, output_ref)
      expect_identical(proteinIntensityFilterPlot, protein_intensity_filter_plot)
      expect_identical(intensityFilterResult$filteredS4, "filtered_s4")
      expect_identical(intensityFilterResult$resultText, "protein intensity complete")
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
    "Applying protein intensity filter..."
  )
  expect_identical(captured$notifications[[1]]$id, "protein_intensity_filter_working")
  expect_null(captured$notifications[[1]]$duration)
  expect_identical(
    captured$notifications[[2]][[1]],
    "Protein intensity filter applied successfully"
  )
  expect_identical(captured$notifications[[2]]$type, "message")
  expect_identical(captured$removed, "protein_intensity_filter_working")
  expect_identical(captured$info, "Protein intensity filter applied successfully")
  expect_identical(completed$status, "success")
  expect_identical(completed$intensityFilterResult$filteredS4, "filtered_s4")
  expect_identical(completed$plotGrid, "plot_grid")
})

test_that("runProteinIntensityFilterApplyObserver reports apply errors and clears the working notification", {
  workflow_data <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()
  captured$removed <- character()
  protein_intensity_filter_plot <- function(...) "plot_reactive"

  completed <- runProteinIntensityFilterApplyObserver(
    workflowData = workflow_data,
    useStrictMode = TRUE,
    minRepsPerGroup = 3,
    minGroups = 2,
    intensityCutoffPercentile = 2.5,
    output = output,
    proteinIntensityFilterPlot = protein_intensity_filter_plot,
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    runApplyStepFn = function(...) {
      captured$calls <- c(captured$calls, "run")
      stop("mock failure")
    },
    updateOutputsFn = function(...) {
      stop("output refresh should not run after a protein intensity error")
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
    "Applying protein intensity filter..."
  )
  expect_identical(captured$notifications[[1]]$id, "protein_intensity_filter_working")
  expect_identical(
    captured$notifications[[2]][[1]],
    "Error applying protein intensity filter: mock failure"
  )
  expect_identical(captured$notifications[[2]]$type, "error")
  expect_identical(captured$notifications[[2]]$duration, 15)
  expect_identical(captured$removed, "protein_intensity_filter_working")
  expect_identical(captured$error, "Error applying protein intensity filter: mock failure")
  expect_identical(completed$status, "error")
  expect_identical(
    completed$errorMessage,
    "Error applying protein intensity filter: mock failure"
  )
})

test_that("runProteinIntensityFilterRevertStep reverts to the previous state in history", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$getHistory <- function() c("imported", "protein_intensity_filtered")
  workflow_data$state_manager$revertToState <- function(state_name) {
    expect_identical(state_name, "imported")
    "reverted_s4"
  }

  result <- runProteinIntensityFilterRevertStep(workflowData = workflow_data)

  expect_identical(result$previousState, "imported")
  expect_identical(result$revertedS4, "reverted_s4")
  expect_identical(result$resultText, "Reverted to previous state: imported")
})

test_that("runProteinIntensityFilterRevertStep errors when there is no prior state", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$getHistory <- function() "protein_intensity_filtered"

  expect_error(
    runProteinIntensityFilterRevertStep(workflowData = workflow_data),
    "No previous state to revert to\\."
  )
})

test_that("runProteinIntensityFilterRevertObserver updates output and reports success", {
  workflow_data <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  output_ref <- output
  workflow_data_ref <- workflow_data
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()

  completed <- runProteinIntensityFilterRevertObserver(
    workflowData = workflow_data,
    output = output,
    runRevertStepFn = function(workflowData) {
      captured$calls <- c(captured$calls, "run")
      expect_identical(workflowData, workflow_data_ref)

      list(
        previousState = "imported",
        revertedS4 = "reverted_s4",
        resultText = "Reverted to previous state: imported"
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
    logInfoFn = function(message) {
      captured$calls <- c(captured$calls, "log_info")
      captured$info <- message
    },
    logErrorFn = function(...) {
      stop("error logging should not run on the happy path")
    }
  )

  expect_identical(output, output_ref)
  expect_identical(captured$calls, c("run", "render", "log_info", "show"))
  expect_identical(
    output$protein_intensity_filter_results,
    "Reverted to previous state: imported"
  )
  expect_identical(captured$info, "Reverted protein intensity filter to imported")
  expect_identical(captured$notifications[[1]][[1]], "Reverted successfully")
  expect_identical(captured$notifications[[1]]$type, "message")
  expect_identical(completed$status, "success")
  expect_identical(completed$revertResult$previousState, "imported")
})

test_that("runProteinIntensityFilterRevertObserver reports revert errors", {
  workflow_data <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()

  completed <- runProteinIntensityFilterRevertObserver(
    workflowData = workflow_data,
    output = output,
    runRevertStepFn = function(...) {
      captured$calls <- c(captured$calls, "run")
      stop("mock revert failure")
    },
    renderTextFn = function(...) {
      stop("rendering should not run after a revert failure")
    },
    showNotificationFn = function(...) {
      captured$calls <- c(captured$calls, "show")
      captured$notifications[[length(captured$notifications) + 1]] <- list(...)
    },
    logInfoFn = function(...) {
      stop("success logging should not run on the error path")
    },
    logErrorFn = function(message) {
      captured$calls <- c(captured$calls, "log_error")
      captured$error <- message
    }
  )

  expect_identical(captured$calls, c("run", "log_error", "show"))
  expect_identical(captured$error, "Error reverting: mock revert failure")
  expect_identical(captured$notifications[[1]][[1]], "Error reverting: mock revert failure")
  expect_identical(captured$notifications[[1]]$type, "error")
  expect_identical(completed$status, "error")
  expect_identical(completed$errorMessage, "Error reverting: mock revert failure")
})

test_that("bindProteinIntensityFilterPlot binds the rendered plot output to the stored grid", {
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$drawn <- NULL

  protein_intensity_filter_plot <- function(value) {
    if (missing(value)) {
      return("plot_grid")
    }

    stop("bindProteinIntensityFilterPlot should only read from the reactive")
  }

  helper_under_test <- makeFunctionWithOverrides(
    bindProteinIntensityFilterPlot,
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

  helper_under_test(
    output = output,
    proteinIntensityFilterPlot = protein_intensity_filter_plot
  )

  expect_identical(output$protein_intensity_filter_plot, "rendered_plot")
  expect_identical(captured$drawn, "plot_grid")
})

test_that("mod_prot_qc_protein_intensity_server wires the apply observer through the helper seam", {
  workflow_data <- shiny::reactiveValues(state_manager = list(kind = "state_manager"))
  captured <- new.env(parent = emptyenv())
  captured$apply <- NULL

  server_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_protein_intensity_server,
    list(
      runProteinIntensityFilterApplyObserver = function(workflowData,
                                                        useStrictMode,
                                                        minRepsPerGroup,
                                                        minGroups,
                                                        intensityCutoffPercentile,
                                                        output,
                                                        proteinIntensityFilterPlot,
                                                        omicType,
                                                        experimentLabel) {
        captured$apply <<- list(
          workflowData = workflowData,
          useStrictMode = useStrictMode,
          minRepsPerGroup = minRepsPerGroup,
          minGroups = minGroups,
          intensityCutoffPercentile = intensityCutoffPercentile,
          output = output,
          proteinIntensityFilterPlot = proteinIntensityFilterPlot,
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
      session$setInputs(use_strict_mode = TRUE)
      session$setInputs(min_reps_per_group = 3)
      session$setInputs(min_groups = 2)
      session$setInputs(proteins_intensity_cutoff_percentile = 1.5)
      session$flushReact()
      session$setInputs(apply_protein_intensity_filter = 1)
      session$flushReact()

      expect_false(is.null(captured$apply))
      expect_identical(captured$apply$workflowData, workflow_data)
      expect_identical(captured$apply$useStrictMode, TRUE)
      expect_identical(captured$apply$minRepsPerGroup, 3)
      expect_identical(captured$apply$minGroups, 2)
      expect_equal(captured$apply$intensityCutoffPercentile, 1.5)
      expect_identical(captured$apply$output, output)
      expect_type(captured$apply$proteinIntensityFilterPlot, "closure")
      expect_identical(captured$apply$omicType, "proteomics")
      expect_identical(captured$apply$experimentLabel, "DIA Experiment")
    }
  )
})

test_that("mod_prot_qc_protein_intensity_server wires the render binding through the helper seam", {
  workflow_data <- shiny::reactiveValues(state_manager = list(kind = "state_manager"))
  captured <- new.env(parent = emptyenv())
  captured$bind <- NULL

  server_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_protein_intensity_server,
    list(
      bindProteinIntensityFilterPlot = function(output, proteinIntensityFilterPlot) {
        captured$bind <<- list(
          output = output,
          proteinIntensityFilterPlot = proteinIntensityFilterPlot
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
      expect_type(captured$bind$proteinIntensityFilterPlot, "closure")
    }
  )
})

test_that("mod_prot_qc_protein_intensity_server wires the revert observer through the helper seam", {
  workflow_data <- shiny::reactiveValues(state_manager = list(kind = "state_manager"))
  captured <- new.env(parent = emptyenv())
  captured$revert <- NULL

  server_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_protein_intensity_server,
    list(
      runProteinIntensityFilterRevertObserver = function(workflowData, output) {
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
      session$setInputs(revert_protein_intensity_filter = 1)
      session$flushReact()

      expect_false(is.null(captured$revert))
      expect_identical(captured$revert$workflowData, workflow_data)
      expect_identical(captured$revert$output, output)
    }
  )
})

test_that("runProteinReplicateFilterApplyStep saves the filtered state and QC parameters", {
  if (!methods::isClass("FakeProteinReplicateState")) {
    methods::setClass(
      "FakeProteinReplicateState",
      slots = c(protein_quant_table = "data.frame")
    )
  }

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  current_state <- methods::new(
    "FakeProteinReplicateState",
    protein_quant_table = data.frame(Protein.Ids = c("P0", "P1", "P2"))
  )
  filtered_state <- methods::new(
    "FakeProteinReplicateState",
    protein_quant_table = data.frame(Protein.Ids = c("P1", "P1", "P2"))
  )
  workflow_data$state_manager$getState <- function() current_state

  captured <- new.env(parent = emptyenv())
  captured$saveState <- NULL
  captured$info <- character()
  captured$written <- NULL
  captured$cluster <- NULL
  captured$filter <- NULL
  captured$saved <- NULL
  workflow_data$state_manager$saveState <- function(...) {
    captured$saveState <- list(...)
  }
  timestamp <- as.POSIXct("2026-04-16 17:15:00", tz = "UTC")
  source_dir <- tempdir()
  protein_qc_dir <- tempdir()

  result <- runProteinReplicateFilterApplyStep(
    workflowData = workflow_data,
    experimentPaths = list(
      protein_qc_dir = protein_qc_dir,
      source_dir = source_dir
    ),
    groupingVariable = "condition",
    parallelCores = 3,
    removeProteinsWithOnlyOneReplicateFn = function(theObject,
                                                    core_utilisation,
                                                    grouping_variable) {
      captured$filter <- list(
        theObject = theObject,
        core_utilisation = core_utilisation,
        grouping_variable = grouping_variable
      )
      filtered_state
    },
    writeTsvFn = function(data, path) {
      captured$written <- list(data = data, path = path)
      invisible(path)
    },
    saveRdsFn = function(object, path) {
      captured$saved <- list(object = object, path = path)
      invisible(path)
    },
    logInfoFn = function(message) {
      captured$info <- c(captured$info, message)
    },
    logWarnFn = function(...) {
      stop("warning path should not run on the happy path")
    },
    newClusterFn = function(cores) {
      captured$cluster <- cores
      paste("cluster", cores)
    },
    nowFn = function() timestamp
  )

  expect_identical(
    captured$info[[1]],
    "Protein Processing: Applying protein replicate filter with 3 cores"
  )
  expect_identical(captured$cluster, 3)
  expect_identical(captured$filter$theObject, current_state)
  expect_identical(captured$filter$core_utilisation, "cluster 3")
  expect_identical(captured$filter$grouping_variable, "condition")
  expect_identical(
    captured$written$path,
    file.path(protein_qc_dir, "remove_proteins_with_only_one_rep.tsv")
  )
  expect_identical(captured$written$data, filtered_state@protein_quant_table)
  expect_identical(
    workflow_data$qc_params$protein_qc$replicate_filter$grouping_variable,
    "condition"
  )
  expect_identical(
    workflow_data$qc_params$protein_qc$replicate_filter$parallel_cores,
    3
  )
  expect_identical(
    workflow_data$qc_params$protein_qc$replicate_filter$timestamp,
    timestamp
  )
  expect_identical(
    captured$saved$path,
    file.path(source_dir, "qc_params.RDS")
  )
  expect_identical(captured$saved$object, workflow_data$qc_params)
  expect_identical(captured$saveState$state_name, "protein_replicate_filtered")
  expect_identical(captured$saveState$s4_data_object, filtered_state)
  expect_identical(captured$saveState$config_object$grouping_variable, "condition")
  expect_identical(captured$saveState$config_object$parallel_cores, 3)
  expect_identical(
    captured$saveState$config_object$output_file,
    file.path(protein_qc_dir, "remove_proteins_with_only_one_rep.tsv")
  )
  expect_identical(
    captured$saveState$description,
    "Applied protein replicate filter (removed single-replicate proteins)"
  )
  expect_identical(workflow_data$protein_counts$after_qc_filtering, 2L)
  expect_identical(result$filteredS4, filtered_state)
  expect_identical(
    result$outputFile,
    file.path(protein_qc_dir, "remove_proteins_with_only_one_rep.tsv")
  )
  expect_identical(result$proteinCount, 2L)
  expect_match(result$resultText, "Proteins remaining: 2", fixed = TRUE)
  expect_match(result$resultText, "Grouping variable: condition", fixed = TRUE)
  expect_match(result$resultText, "Parallel cores used: 3", fixed = TRUE)
})

test_that("updateProteinReplicateFilterOutputs refreshes result text and plot grid", {
  if (!methods::isClass("FakeProteinReplicateState")) {
    methods::setClass(
      "FakeProteinReplicateState",
      slots = c(protein_quant_table = "data.frame")
    )
  }

  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$plot <- NULL
  captured$args <- NULL

  protein_replicate_filter_plot <- function(value) {
    captured$plot <- value
    invisible(value)
  }

  plot_grid <- updateProteinReplicateFilterOutputs(
    output = output,
    proteinReplicateFilterPlot = protein_replicate_filter_plot,
    replicateFilterResult = list(
      resultText = "protein replicate complete",
      filteredS4 = methods::new(
        "FakeProteinReplicateState",
        protein_quant_table = data.frame(Protein.Ids = c("P1", "P2"))
      )
    ),
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    renderTextFn = function(value) value,
    updateProteinFilteringFn = function(data,
                                        step_name,
                                        omic_type,
                                        experiment_label,
                                        return_grid,
                                        overwrite) {
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

  expect_identical(output$protein_replicate_filter_results, "protein replicate complete")
  expect_identical(plot_grid, "plot_grid")
  expect_identical(captured$plot, "plot_grid")
  expect_identical(captured$args$data$Protein.Ids, c("P1", "P2"))
  expect_identical(captured$args$step_name, "13_protein_replicate_filtered")
  expect_identical(captured$args$omic_type, "proteomics")
  expect_identical(captured$args$experiment_label, "DIA Experiment")
  expect_true(captured$args$return_grid)
  expect_true(captured$args$overwrite)
})

test_that("runProteinReplicateFilterApplyObserver delegates the apply workflow and notifications", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data_ref <- workflow_data
  output <- new.env(parent = emptyenv())
  output_ref <- output
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()
  captured$removed <- character()
  protein_replicate_filter_plot <- function(...) "plot_reactive"
  experiment_paths <- list(protein_qc_dir = tempdir(), source_dir = tempdir())

  completed <- runProteinReplicateFilterApplyObserver(
    workflowData = workflow_data,
    experimentPaths = experiment_paths,
    groupingVariable = "condition",
    parallelCores = 4,
    output = output,
    proteinReplicateFilterPlot = protein_replicate_filter_plot,
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    runApplyStepFn = function(workflowData,
                              experimentPaths,
                              groupingVariable,
                              parallelCores) {
      captured$calls <- c(captured$calls, "run")
      expect_identical(workflowData, workflow_data_ref)
      expect_identical(experimentPaths, experiment_paths)
      expect_identical(groupingVariable, "condition")
      expect_identical(parallelCores, 4)
      list(
        filteredS4 = "filtered_s4",
        resultText = "protein replicate complete"
      )
    },
    updateOutputsFn = function(output,
                               proteinReplicateFilterPlot,
                               replicateFilterResult,
                               omicType,
                               experimentLabel) {
      captured$calls <- c(captured$calls, "update")
      expect_identical(output, output_ref)
      expect_identical(proteinReplicateFilterPlot, protein_replicate_filter_plot)
      expect_identical(replicateFilterResult$filteredS4, "filtered_s4")
      expect_identical(replicateFilterResult$resultText, "protein replicate complete")
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
    "Applying protein replicate filter..."
  )
  expect_identical(captured$notifications[[1]]$id, "protein_replicate_filter_working")
  expect_null(captured$notifications[[1]]$duration)
  expect_identical(
    captured$notifications[[2]][[1]],
    "Protein replicate filter applied successfully. Pipeline complete!"
  )
  expect_identical(captured$notifications[[2]]$type, "message")
  expect_identical(captured$removed, "protein_replicate_filter_working")
  expect_identical(captured$info, "Protein replicate filter applied successfully")
  expect_identical(completed$status, "success")
  expect_identical(completed$replicateFilterResult$filteredS4, "filtered_s4")
  expect_identical(completed$plotGrid, "plot_grid")
})

test_that("runProteinReplicateFilterApplyObserver reports apply errors and clears the working notification", {
  workflow_data <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()
  captured$removed <- character()
  protein_replicate_filter_plot <- function(...) "plot_reactive"

  completed <- runProteinReplicateFilterApplyObserver(
    workflowData = workflow_data,
    experimentPaths = list(protein_qc_dir = tempdir()),
    groupingVariable = "condition",
    parallelCores = 2,
    output = output,
    proteinReplicateFilterPlot = protein_replicate_filter_plot,
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    runApplyStepFn = function(...) {
      captured$calls <- c(captured$calls, "run")
      stop("mock failure")
    },
    updateOutputsFn = function(...) {
      stop("output refresh should not run after a protein replicate error")
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
    "Applying protein replicate filter..."
  )
  expect_identical(captured$notifications[[1]]$id, "protein_replicate_filter_working")
  expect_identical(
    captured$notifications[[2]][[1]],
    "Error applying protein replicate filter: mock failure"
  )
  expect_identical(captured$notifications[[2]]$type, "error")
  expect_identical(captured$notifications[[2]]$duration, 15)
  expect_identical(captured$removed, "protein_replicate_filter_working")
  expect_identical(captured$error, "Error applying protein replicate filter: mock failure")
  expect_identical(completed$status, "error")
  expect_identical(
    completed$errorMessage,
    "Error applying protein replicate filter: mock failure"
  )
})

test_that("runProteinReplicateFilterRevertStep reverts to the previous state in history", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$getHistory <- function() {
    c("sample_filtered", "protein_replicate_filtered")
  }
  workflow_data$state_manager$revertToState <- function(state_name) {
    expect_identical(state_name, "sample_filtered")
    "reverted_s4"
  }

  result <- runProteinReplicateFilterRevertStep(workflowData = workflow_data)

  expect_identical(result$previousState, "sample_filtered")
  expect_identical(result$revertedS4, "reverted_s4")
  expect_identical(result$resultText, "Reverted to previous state: sample_filtered")
})

test_that("runProteinReplicateFilterRevertStep errors when there is no prior state", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$getHistory <- function() "protein_replicate_filtered"

  expect_error(
    runProteinReplicateFilterRevertStep(workflowData = workflow_data),
    "No previous state to revert to\\."
  )
})

test_that("runProteinReplicateFilterRevertObserver updates output and reports success", {
  workflow_data <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  output_ref <- output
  workflow_data_ref <- workflow_data
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()

  completed <- runProteinReplicateFilterRevertObserver(
    workflowData = workflow_data,
    output = output,
    runRevertStepFn = function(workflowData) {
      captured$calls <- c(captured$calls, "run")
      expect_identical(workflowData, workflow_data_ref)

      list(
        previousState = "sample_filtered",
        revertedS4 = "reverted_s4",
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
    logInfoFn = function(message) {
      captured$calls <- c(captured$calls, "log_info")
      captured$info <- message
    },
    logErrorFn = function(...) {
      stop("error logging should not run on the happy path")
    }
  )

  expect_identical(output, output_ref)
  expect_identical(captured$calls, c("run", "render", "log_info", "show"))
  expect_identical(
    output$protein_replicate_filter_results,
    "Reverted to previous state: sample_filtered"
  )
  expect_identical(captured$info, "Reverted protein replicate filter to sample_filtered")
  expect_identical(captured$notifications[[1]][[1]], "Reverted successfully")
  expect_identical(captured$notifications[[1]]$type, "message")
  expect_identical(completed$status, "success")
  expect_identical(completed$revertResult$previousState, "sample_filtered")
})

test_that("runProteinReplicateFilterRevertObserver reports revert errors", {
  workflow_data <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()

  completed <- runProteinReplicateFilterRevertObserver(
    workflowData = workflow_data,
    output = output,
    runRevertStepFn = function(...) {
      captured$calls <- c(captured$calls, "run")
      stop("mock revert failure")
    },
    renderTextFn = function(...) {
      stop("rendering should not run after a revert failure")
    },
    showNotificationFn = function(...) {
      captured$calls <- c(captured$calls, "show")
      captured$notifications[[length(captured$notifications) + 1]] <- list(...)
    },
    logInfoFn = function(...) {
      stop("success logging should not run on the error path")
    },
    logErrorFn = function(message) {
      captured$calls <- c(captured$calls, "log_error")
      captured$error <- message
    }
  )

  expect_identical(captured$calls, c("run", "log_error", "show"))
  expect_identical(captured$error, "Error reverting: mock revert failure")
  expect_identical(captured$notifications[[1]][[1]], "Error reverting: mock revert failure")
  expect_identical(captured$notifications[[1]]$type, "error")
  expect_identical(completed$status, "error")
  expect_identical(completed$errorMessage, "Error reverting: mock revert failure")
})

test_that("bindProteinReplicateFilterPlot binds the rendered plot output to the stored grid", {
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$drawn <- NULL

  protein_replicate_filter_plot <- function(value) {
    if (missing(value)) {
      return("plot_grid")
    }

    stop("bindProteinReplicateFilterPlot should only read from the reactive")
  }

  helper_under_test <- makeFunctionWithOverrides(
    bindProteinReplicateFilterPlot,
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

  helper_under_test(
    output = output,
    proteinReplicateFilterPlot = protein_replicate_filter_plot
  )

  expect_identical(output$protein_replicate_filter_plot, "rendered_plot")
  expect_identical(captured$drawn, "plot_grid")
})

test_that("mod_prot_qc_protein_replicate_server wires the apply observer through the helper seam", {
  workflow_data <- shiny::reactiveValues(state_manager = list(kind = "state_manager"))
  experiment_paths <- list(protein_qc_dir = tempdir(), source_dir = tempdir())
  captured <- new.env(parent = emptyenv())
  captured$apply <- NULL

  server_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_protein_replicate_server,
    list(
      runProteinReplicateFilterApplyObserver = function(workflowData,
                                                        experimentPaths,
                                                        groupingVariable,
                                                        parallelCores,
                                                        output,
                                                        proteinReplicateFilterPlot,
                                                        omicType,
                                                        experimentLabel) {
        captured$apply <<- list(
          workflowData = workflowData,
          experimentPaths = experimentPaths,
          groupingVariable = groupingVariable,
          parallelCores = parallelCores,
          output = output,
          proteinReplicateFilterPlot = proteinReplicateFilterPlot,
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
      experiment_paths = experiment_paths,
      omic_type = "proteomics",
      experiment_label = "DIA Experiment"
    ),
    {
      session$setInputs(protein_grouping_variable = "condition")
      session$setInputs(parallel_cores = 5)
      session$flushReact()
      session$setInputs(apply_protein_replicate_filter = 1)
      session$flushReact()

      expect_false(is.null(captured$apply))
      expect_identical(captured$apply$workflowData, workflow_data)
      expect_identical(captured$apply$experimentPaths, experiment_paths)
      expect_identical(captured$apply$groupingVariable, "condition")
      expect_identical(captured$apply$parallelCores, 5)
      expect_identical(captured$apply$output, output)
      expect_type(captured$apply$proteinReplicateFilterPlot, "closure")
      expect_identical(captured$apply$omicType, "proteomics")
      expect_identical(captured$apply$experimentLabel, "DIA Experiment")
    }
  )
})

test_that("mod_prot_qc_protein_replicate_server wires the render binding through the helper seam", {
  workflow_data <- shiny::reactiveValues(state_manager = list(kind = "state_manager"))
  experiment_paths <- list(protein_qc_dir = tempdir(), source_dir = tempdir())
  captured <- new.env(parent = emptyenv())
  captured$bind <- NULL

  server_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_protein_replicate_server,
    list(
      bindProteinReplicateFilterPlot = function(output, proteinReplicateFilterPlot) {
        captured$bind <<- list(
          output = output,
          proteinReplicateFilterPlot = proteinReplicateFilterPlot
        )
        invisible(NULL)
      }
    )
  )

  testServer(
    server_under_test,
    args = list(
      workflow_data = workflow_data,
      experiment_paths = experiment_paths,
      omic_type = "proteomics",
      experiment_label = "DIA Experiment"
    ),
    {
      expect_false(is.null(captured$bind))
      expect_identical(captured$bind$output, output)
      expect_type(captured$bind$proteinReplicateFilterPlot, "closure")
    }
  )
})

test_that("mod_prot_qc_protein_replicate_server wires the revert observer through the helper seam", {
  workflow_data <- shiny::reactiveValues(state_manager = list(kind = "state_manager"))
  experiment_paths <- list(protein_qc_dir = tempdir(), source_dir = tempdir())
  captured <- new.env(parent = emptyenv())
  captured$revert <- NULL

  server_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_protein_replicate_server,
    list(
      runProteinReplicateFilterRevertObserver = function(workflowData, output) {
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
      experiment_paths = experiment_paths,
      omic_type = "proteomics",
      experiment_label = "DIA Experiment"
    ),
    {
      session$setInputs(revert_protein_replicate_filter = 1)
      session$flushReact()

      expect_false(is.null(captured$revert))
      expect_identical(captured$revert$workflowData, workflow_data)
      expect_identical(captured$revert$output, output)
    }
  )
})

test_that("runProteinIqRollupApplyStep saves the protein S4 state and restores sample ids", {
  if (!methods::isClass("FakeProteinIqInputState")) {
    methods::setClass(
      "FakeProteinIqInputState",
      slots = c(
        peptide_data = "data.frame",
        design_matrix = "data.frame",
        args = "list"
      )
    )
  }
  if (!methods::isClass("FakeProteinIqOutputState")) {
    methods::setClass(
      "FakeProteinIqOutputState",
      slots = c(
        protein_quant_table = "data.frame",
        design_matrix = "data.frame"
      )
    )
  }

  peptide_state <- methods::new(
    "FakeProteinIqInputState",
    peptide_data = data.frame(
      Protein.Ids = c("P1", "P1"),
      Stripped.Sequence = c("PEP1", "PEP1"),
      Run = c("sample_a", "sample_b"),
      Peptide.Imputed = c(10, NA_real_),
      stringsAsFactors = FALSE
    ),
    design_matrix = data.frame(
      Run = c("sample_a", "sample_b"),
      group = c("G1", "G2"),
      replicates = c(1, 1),
      stringsAsFactors = FALSE
    ),
    args = list(source = "test")
  )

  peptide_dir <- tempfile("protein-iq-peptide-")
  protein_dir <- tempfile("protein-iq-output-")
  dir.create(peptide_dir)
  dir.create(protein_dir)

  captured <- new.env(parent = emptyenv())
  captured$info <- character()
  captured$notifications <- list()

  state_manager <- new.env(parent = emptyenv())
  state_manager$current_state <- "imputed"
  state_manager$getState <- function(stateName) {
    expect_identical(stateName, "imputed")
    peptide_state
  }
  state_manager$saveState <- function(state_name, s4_data_object, config_object, description) {
    captured$saveState <- list(
      state_name = state_name,
      s4_data_object = s4_data_object,
      config_object = config_object,
      description = description
    )
    invisible(NULL)
  }

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- state_manager

  result <- runProteinIqRollupApplyStep(
    workflowData = workflow_data,
    experimentPaths = list(
      peptide_qc_dir = peptide_dir,
      protein_qc_dir = protein_dir
    ),
    createProteinDataFn = function(protein_quant_table,
                                   protein_id_column,
                                   protein_id_table,
                                   design_matrix,
                                   sample_id,
                                   group_id,
                                   technical_replicate_id,
                                   args) {
      captured$create <- list(
        protein_quant_table = protein_quant_table,
        protein_id_column = protein_id_column,
        protein_id_table = protein_id_table,
        design_matrix = design_matrix,
        sample_id = sample_id,
        group_id = group_id,
        technical_replicate_id = technical_replicate_id,
        args = args
      )
      methods::new(
        "FakeProteinIqOutputState",
        protein_quant_table = protein_quant_table,
        design_matrix = design_matrix
      )
    },
    writeTsvFn = function(data, path) {
      captured$write <- list(data = data, path = path)
      invisible(NULL)
    },
    processLongFormatFn = function(input_filename,
                                   output_filename,
                                   sample_id,
                                   primary_id,
                                   secondary_id,
                                   intensity_col,
                                   filter_double_less,
                                   normalization) {
      captured$process <- list(
        input_filename = input_filename,
        output_filename = output_filename,
        sample_id = sample_id,
        primary_id = primary_id,
        secondary_id = secondary_id,
        intensity_col = intensity_col,
        filter_double_less = filter_double_less,
        normalization = normalization
      )
      file.create(output_filename)
      invisible(NULL)
    },
    readTsvFn = function(path, .name_repair) {
      captured$read <- list(path = path, name_repair = .name_repair)
      data.frame(
        Protein.Ids = "P1",
        S_001 = 42,
        check.names = FALSE
      )
    },
    captureCheckpointFn = function(object, checkpointId, checkpointLabel) {
      captured$checkpoint <- list(
        object = object,
        checkpointId = checkpointId,
        checkpointLabel = checkpointLabel
      )
      invisible(NULL)
    },
    showNotificationFn = function(...) {
      captured$notifications[[length(captured$notifications) + 1]] <- list(...)
    },
    logInfoFn = function(message) {
      captured$info <- c(captured$info, message)
    },
    logWarnFn = function(message) {
      captured$warn <- message
    },
    sleepFn = function(...) {
      stop("sleep should not run when the mock IQ output already exists")
    }
  )

  expect_identical(captured$write$path, file.path(peptide_dir, "peptide_values_imputed.tsv"))
  expect_identical(captured$write$data$Run, c("S_001", "S_002"))
  expect_identical(captured$process$sample_id, "Run")
  expect_identical(captured$process$primary_id, "Protein.Ids")
  expect_identical(captured$process$secondary_id, "Stripped.Sequence")
  expect_identical(captured$process$intensity_col, "Peptide.Imputed")
  expect_identical(captured$process$normalization, "none")
  expect_identical(captured$read$name_repair, "minimal")
  expect_identical(colnames(captured$create$protein_quant_table), c("Protein.Ids", "sample_a"))
  expect_identical(captured$create$protein_id_column, "Protein.Ids")
  expect_identical(captured$create$design_matrix$Run, "sample_a")
  expect_identical(captured$create$sample_id, "Run")
  expect_identical(captured$create$group_id, "group")
  expect_identical(captured$create$technical_replicate_id, "replicates")
  expect_identical(captured$create$args, list(source = "test"))
  expect_identical(captured$saveState$state_name, "protein_s4_created")
  expect_identical(captured$saveState$s4_data_object, result$proteinObj)
  expect_identical(
    captured$saveState$config_object$iq_output_file,
    file.path(protein_dir, "iq_output_file.txt")
  )
  expect_identical(captured$saveState$config_object$s4_class, "ProteinQuantitativeData")
  expect_identical(captured$checkpoint$object, result$proteinObj)
  expect_identical(captured$checkpoint$checkpointId, "cp03")
  expect_identical(captured$checkpoint$checkpointLabel, "rolled_up_protein")
  expect_identical(
    captured$warn,
    "Protein Processing: 1 samples dropped during IQ rollup: sample_b"
  )
  expect_match(
    captured$notifications[[1]][[1]],
    "Warning: The following samples were dropped during protein rollup due to low peptide quality: sample_b",
    fixed = TRUE
  )
  expect_identical(captured$notifications[[1]]$type, "warning")
  expect_identical(captured$notifications[[1]]$duration, 10)
  expect_identical(
    captured$info,
    c(
      "Protein Processing: Starting IQ rollup from peptide state",
      "Protein Processing: Creating ProteinQuantitativeData S4 object"
    )
  )
  expect_match(result$resultText, "Proteins quantified: 1", fixed = TRUE)
  expect_match(result$resultText, "Samples: 1", fixed = TRUE)
  expect_match(result$resultText, "Output file: iq_output_file.txt", fixed = TRUE)
})

test_that("updateProteinIqRollupOutputs refreshes result text and plot grid", {
  if (!methods::isClass("FakeProteinIqOutputState")) {
    methods::setClass(
      "FakeProteinIqOutputState",
      slots = c(
        protein_quant_table = "data.frame",
        design_matrix = "data.frame"
      )
    )
  }

  protein_state <- methods::new(
    "FakeProteinIqOutputState",
    protein_quant_table = data.frame(Protein.Ids = c("P1", "P2")),
    design_matrix = data.frame(Run = c("S1", "S2"))
  )
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$plot <- NULL

  plot_grid <- updateProteinIqRollupOutputs(
    output = output,
    iqRollupPlot = function(value) {
      captured$plot <- value
      invisible(value)
    },
    iqRollupResult = list(
      proteinObj = protein_state,
      resultText = "IQ rollup complete"
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

  expect_identical(output$iq_rollup_results, "IQ rollup complete")
  expect_identical(plot_grid, "plot_grid")
  expect_identical(captured$plot, "plot_grid")
  expect_identical(captured$args$data, protein_state@protein_quant_table)
  expect_identical(captured$args$step_name, "9_protein_s4_created")
  expect_identical(captured$args$omic_type, "proteomics")
  expect_identical(captured$args$experiment_label, "DIA Experiment")
  expect_true(captured$args$return_grid)
  expect_true(captured$args$overwrite)
})

test_that("runProteinIqRollupApplyObserver delegates the apply workflow and notifications", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data_ref <- workflow_data
  output <- new.env(parent = emptyenv())
  output_ref <- output
  experiment_paths <- list(peptide_qc_dir = tempdir(), protein_qc_dir = tempdir())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()
  captured$removed <- character()
  iq_rollup_plot <- function(...) "plot_reactive"

  completed <- runProteinIqRollupApplyObserver(
    workflowData = workflow_data,
    experimentPaths = experiment_paths,
    output = output,
    iqRollupPlot = iq_rollup_plot,
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    runApplyStepFn = function(workflowData, experimentPaths) {
      captured$calls <- c(captured$calls, "run")
      expect_identical(workflowData, workflow_data_ref)
      expect_identical(experimentPaths, experiment_paths)
      list(
        proteinObj = "protein_s4",
        resultText = "IQ rollup complete"
      )
    },
    updateOutputsFn = function(output, iqRollupPlot, iqRollupResult, omicType, experimentLabel) {
      captured$calls <- c(captured$calls, "update")
      expect_identical(output, output_ref)
      expect_identical(iqRollupPlot, iq_rollup_plot)
      expect_identical(iqRollupResult$proteinObj, "protein_s4")
      expect_identical(iqRollupResult$resultText, "IQ rollup complete")
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
    "Running IQ protein rollup & creating S4 object..."
  )
  expect_identical(captured$notifications[[1]]$id, "iq_rollup_working")
  expect_null(captured$notifications[[1]]$duration)
  expect_identical(
    captured$notifications[[2]][[1]],
    "IQ protein rollup & S4 object creation completed successfully"
  )
  expect_identical(captured$notifications[[2]]$type, "message")
  expect_identical(captured$removed, "iq_rollup_working")
  expect_identical(captured$info, "IQ protein rollup and S4 object creation completed successfully")
  expect_identical(completed$status, "success")
  expect_identical(completed$iqRollupResult$proteinObj, "protein_s4")
  expect_identical(completed$plotGrid, "plot_grid")
})

test_that("runProteinIqRollupApplyObserver reports apply errors and clears the working notification", {
  workflow_data <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  experiment_paths <- list(peptide_qc_dir = tempdir(), protein_qc_dir = tempdir())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()
  captured$removed <- character()

  completed <- runProteinIqRollupApplyObserver(
    workflowData = workflow_data,
    experimentPaths = experiment_paths,
    output = output,
    iqRollupPlot = function(...) "plot_reactive",
    omicType = "proteomics",
    experimentLabel = "DIA Experiment",
    runApplyStepFn = function(...) {
      captured$calls <- c(captured$calls, "run")
      stop("mock failure")
    },
    updateOutputsFn = function(...) {
      stop("output refresh should not run after an IQ rollup error")
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
    "Running IQ protein rollup & creating S4 object..."
  )
  expect_identical(captured$notifications[[1]]$id, "iq_rollup_working")
  expect_identical(
    captured$notifications[[2]][[1]],
    "Error in IQ protein rollup & S4 creation: mock failure"
  )
  expect_identical(captured$notifications[[2]]$type, "error")
  expect_identical(captured$notifications[[2]]$duration, 15)
  expect_identical(captured$removed, "iq_rollup_working")
  expect_identical(captured$error, "Error in IQ protein rollup & S4 creation: mock failure")
  expect_identical(completed$status, "error")
  expect_identical(
    completed$errorMessage,
    "Error in IQ protein rollup & S4 creation: mock failure"
  )
})

test_that("mod_prot_qc_protein_rollup_server wires the apply observer through the helper seam", {
  workflow_data <- shiny::reactiveValues(state_manager = list(kind = "state_manager"))
  experiment_paths <- list(peptide_qc_dir = tempdir(), protein_qc_dir = tempdir())
  captured <- new.env(parent = emptyenv())
  captured$apply <- NULL

  server_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_protein_rollup_server,
    list(
      runProteinIqRollupApplyObserver = function(workflowData,
                                                 experimentPaths,
                                                 output,
                                                 iqRollupPlot,
                                                 omicType,
                                                 experimentLabel) {
        captured$apply <<- list(
          workflowData = workflowData,
          experimentPaths = experimentPaths,
          output = output,
          iqRollupPlot = iqRollupPlot,
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
      experiment_paths = experiment_paths,
      omic_type = "proteomics",
      experiment_label = "DIA Experiment"
    ),
    {
      session$setInputs(apply_iq_rollup = 1)
      session$flushReact()

      expect_false(is.null(captured$apply))
      expect_identical(captured$apply$workflowData, workflow_data)
      expect_identical(captured$apply$experimentPaths, experiment_paths)
      expect_identical(captured$apply$output, output)
      expect_type(captured$apply$iqRollupPlot, "closure")
      expect_identical(captured$apply$omicType, "proteomics")
      expect_identical(captured$apply$experimentLabel, "DIA Experiment")
    }
  )
})

test_that("runProteinIqRollupRevertStep reverts to the latest peptide state in history", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$getHistory <- function() {
    c("raw_data_s4", "qvalue_filtered", "protein_s4_created")
  }
  workflow_data$state_manager$revertToState <- function(state_name) {
    expect_identical(state_name, "qvalue_filtered")
    "reverted_s4"
  }

  result <- runProteinIqRollupRevertStep(workflowData = workflow_data)

  expect_identical(result$previousState, "qvalue_filtered")
  expect_identical(result$revertedS4, "reverted_s4")
  expect_identical(result$resultText, "Reverted to qvalue_filtered state")
})

test_that("runProteinIqRollupRevertStep errors when there is no peptide state to revert to", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$state_manager <- new.env(parent = emptyenv())
  workflow_data$state_manager$getHistory <- function() "protein_s4_created"

  expect_error(
    runProteinIqRollupRevertStep(workflowData = workflow_data),
    "No previous peptide state to revert to\\."
  )
})

test_that("runProteinIqRollupRevertObserver updates output and reports success", {
  workflow_data <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  workflow_data_ref <- workflow_data
  output_ref <- output
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()

  completed <- runProteinIqRollupRevertObserver(
    workflowData = workflow_data,
    output = output,
    runRevertStepFn = function(workflowData) {
      captured$calls <- c(captured$calls, "run")
      expect_identical(workflowData, workflow_data_ref)

      list(
        previousState = "qvalue_filtered",
        revertedS4 = "reverted_s4",
        resultText = "Reverted to qvalue_filtered state"
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
    logInfoFn = function(message) {
      captured$calls <- c(captured$calls, "log_info")
      captured$info <- message
    },
    logErrorFn = function(...) {
      stop("error logging should not run on the happy path")
    }
  )

  expect_identical(output, output_ref)
  expect_identical(captured$calls, c("run", "render", "log_info", "show"))
  expect_identical(output$iq_rollup_results, "Reverted to qvalue_filtered state")
  expect_identical(captured$info, "Reverted IQ rollup to qvalue_filtered")
  expect_identical(captured$notifications[[1]][[1]], "Reverted successfully")
  expect_identical(captured$notifications[[1]]$type, "message")
  expect_identical(completed$status, "success")
  expect_identical(completed$revertResult$previousState, "qvalue_filtered")
})

test_that("runProteinIqRollupRevertObserver reports revert errors", {
  workflow_data <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()

  completed <- runProteinIqRollupRevertObserver(
    workflowData = workflow_data,
    output = output,
    runRevertStepFn = function(...) {
      captured$calls <- c(captured$calls, "run")
      stop("mock revert failure")
    },
    renderTextFn = function(...) {
      stop("rendering should not run after a revert failure")
    },
    showNotificationFn = function(...) {
      captured$calls <- c(captured$calls, "show")
      captured$notifications[[length(captured$notifications) + 1]] <- list(...)
    },
    logInfoFn = function(...) {
      stop("success logging should not run on the error path")
    },
    logErrorFn = function(message) {
      captured$calls <- c(captured$calls, "log_error")
      captured$error <- message
    }
  )

  expect_identical(captured$calls, c("run", "log_error", "show"))
  expect_identical(captured$error, "Error reverting: mock revert failure")
  expect_identical(captured$notifications[[1]][[1]], "Error reverting: mock revert failure")
  expect_identical(captured$notifications[[1]]$type, "error")
  expect_identical(completed$status, "error")
  expect_identical(completed$errorMessage, "Error reverting: mock revert failure")
})

test_that("bindProteinIqRollupPlot binds the rendered plot output to the stored grid", {
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$drawn <- NULL

  iq_rollup_plot <- function(value) {
    if (missing(value)) {
      return("plot_grid")
    }

    stop("bindProteinIqRollupPlot should only read from the reactive")
  }

  helper_under_test <- makeFunctionWithOverrides(
    bindProteinIqRollupPlot,
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

  helper_under_test(output = output, iqRollupPlot = iq_rollup_plot)

  expect_identical(output$iq_rollup_plot, "rendered_plot")
  expect_identical(captured$drawn, "plot_grid")
})

test_that("mod_prot_qc_protein_rollup_server wires the render binding through the helper seam", {
  workflow_data <- shiny::reactiveValues(state_manager = list(kind = "state_manager"))
  experiment_paths <- list(peptide_qc_dir = tempdir(), protein_qc_dir = tempdir())
  captured <- new.env(parent = emptyenv())
  captured$bind <- NULL

  server_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_protein_rollup_server,
    list(
      bindProteinIqRollupPlot = function(output, iqRollupPlot) {
        captured$bind <<- list(
          output = output,
          iqRollupPlot = iqRollupPlot
        )
        invisible(NULL)
      }
    )
  )

  testServer(
    server_under_test,
    args = list(
      workflow_data = workflow_data,
      experiment_paths = experiment_paths,
      omic_type = "proteomics",
      experiment_label = "DIA Experiment"
    ),
    {
      expect_false(is.null(captured$bind))
      expect_identical(captured$bind$output, output)
      expect_type(captured$bind$iqRollupPlot, "closure")
    }
  )
})

test_that("mod_prot_qc_protein_rollup_server wires the revert observer through the helper seam", {
  workflow_data <- shiny::reactiveValues(state_manager = list(kind = "state_manager"))
  experiment_paths <- list(peptide_qc_dir = tempdir(), protein_qc_dir = tempdir())
  captured <- new.env(parent = emptyenv())
  captured$revert <- NULL

  server_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_protein_rollup_server,
    list(
      runProteinIqRollupRevertObserver = function(workflowData, output) {
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
      experiment_paths = experiment_paths,
      omic_type = "proteomics",
      experiment_label = "DIA Experiment"
    ),
    {
      session$setInputs(revert_iq_rollup = 1)
      session$flushReact()

      expect_false(is.null(captured$revert))
      expect_identical(captured$revert$workflowData, workflow_data)
      expect_identical(captured$revert$output, output)
    }
  )
})

test_that("runProteinS4CreationStep creates the initial protein S4 state and result text", {
  if (!methods::isClass("FakeProteinS4State")) {
    methods::setClass(
      "FakeProteinS4State",
      slots = c(protein_quant_table = "data.frame")
    )
  }

  workflow_data <- new.env(parent = emptyenv())
  workflow_data$data_cln <- data.frame(
    Protein.Ids = c("P1", "P1", "P2"),
    Run = c("S1", "S2", "S1"),
    group = c("A", "A", "B"),
    replicates = c("R1", "R2", "R1"),
    stringsAsFactors = FALSE
  )
  workflow_data$design_matrix <- data.frame(
    Run = c("S1", "S2"),
    group = c("A", "A"),
    replicates = c("R1", "R2"),
    stringsAsFactors = FALSE
  )
  workflow_data$column_mapping <- list(protein_col = "Protein.Ids")
  workflow_data$config_list <- list(normalize = FALSE)
  workflow_data$state_manager <- new.env(parent = emptyenv())

  captured <- new.env(parent = emptyenv())
  captured$saveState <- NULL
  captured$logInfo <- NULL
  workflow_data$state_manager$saveState <- function(...) {
    captured$saveState <- list(...)
  }

  result <- runProteinS4CreationStep(
    workflowData = workflow_data,
    proteinQuantitativeDataFn = function(...) {
      captured$constructorArgs <- list(...)
      methods::new(
        "FakeProteinS4State",
        protein_quant_table = workflow_data$data_cln
      )
    },
    logInfoFn = function(message) {
      captured$logInfo <- message
    }
  )

  expect_identical(
    captured$logInfo,
    "Protein S4 Creation: Starting process for protein-level workflow (e.g., TMT)"
  )
  expect_identical(captured$constructorArgs$protein_quant_table, workflow_data$data_cln)
  expect_identical(captured$constructorArgs$protein_id_column, "Protein.Ids")
  expect_identical(captured$constructorArgs$protein_id_table$Protein.Ids, c("P1", "P2"))
  expect_identical(captured$constructorArgs$design_matrix, workflow_data$design_matrix)
  expect_identical(captured$constructorArgs$sample_id, "Run")
  expect_identical(captured$constructorArgs$group_id, "group")
  expect_identical(captured$constructorArgs$technical_replicate_id, "replicates")
  expect_identical(captured$constructorArgs$args, workflow_data$config_list)
  expect_identical(captured$saveState$state_name, "protein_s4_initial")
  expect_identical(captured$saveState$s4_data_object, result$proteinObj)
  expect_identical(captured$saveState$config_object$s4_class, "ProteinQuantitativeData")
  expect_identical(captured$saveState$config_object$protein_id_column, "Protein.Ids")
  expect_identical(
    captured$saveState$description,
    "Created initial ProteinQuantitativeData S4 object from protein-level data."
  )
  expect_identical(result$proteinCount, 2L)
  expect_identical(result$proteinIdCol, "Protein.Ids")
  expect_match(result$resultText, "Proteins loaded: 2", fixed = TRUE)
  expect_match(result$resultText, "S4 Class: FakeProteinS4State", fixed = TRUE)
  expect_match(result$resultText, "State saved as: 'protein_s4_initial'", fixed = TRUE)
})

test_that("runProteinS4CreationStep errors when the configured protein id column is missing", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data$data_cln <- data.frame(Run = "S1", stringsAsFactors = FALSE)
  workflow_data$design_matrix <- data.frame(Run = "S1", stringsAsFactors = FALSE)
  workflow_data$column_mapping <- list(protein_col = "Protein.Ids")
  workflow_data$config_list <- list()
  workflow_data$state_manager <- new.env(parent = emptyenv())

  expect_error(
    runProteinS4CreationStep(
      workflowData = workflow_data,
      proteinQuantitativeDataFn = function(...) stop("constructor should not run")
    ),
    "Protein ID column Protein\\.Ids not found in the data\\."
  )
})

test_that("runProteinS4CreationObserver delegates the create workflow and notifications", {
  workflow_data <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  output_ref <- output
  workflow_data_ref <- workflow_data
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()
  captured$removed <- character()

  completed <- runProteinS4CreationObserver(
    workflowData = workflow_data,
    output = output,
    runCreationStepFn = function(workflowData) {
      captured$calls <- c(captured$calls, "run")
      expect_identical(workflowData, workflow_data_ref)
      list(
        proteinObj = "protein_s4",
        proteinCount = 2L,
        proteinIdCol = "Protein.Ids",
        resultText = "creation complete"
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

  expect_identical(output, output_ref)
  expect_identical(
    captured$calls,
    c("show", "run", "render", "log_info", "remove", "show")
  )
  expect_identical(
    captured$notifications[[1]][[1]],
    "Creating Protein S4 object from imported data..."
  )
  expect_identical(captured$notifications[[1]]$id, "s4_creation_working")
  expect_null(captured$notifications[[1]]$duration)
  expect_identical(output$s4_creation_results, "creation complete")
  expect_identical(
    captured$notifications[[2]][[1]],
    "Protein S4 object created successfully"
  )
  expect_identical(captured$notifications[[2]]$type, "message")
  expect_identical(captured$removed, "s4_creation_working")
  expect_identical(
    captured$info,
    "Protein S4 object creation from protein data completed successfully"
  )
  expect_identical(completed$status, "success")
  expect_identical(completed$creationResult$proteinObj, "protein_s4")
})

test_that("runProteinS4CreationObserver reports create errors and clears the working notification", {
  workflow_data <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()
  captured$removed <- character()

  completed <- runProteinS4CreationObserver(
    workflowData = workflow_data,
    output = output,
    runCreationStepFn = function(...) {
      captured$calls <- c(captured$calls, "run")
      stop("mock failure")
    },
    renderTextFn = function(...) {
      stop("rendering should not run after a create failure")
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
    "Creating Protein S4 object from imported data..."
  )
  expect_identical(captured$notifications[[1]]$id, "s4_creation_working")
  expect_identical(
    captured$notifications[[2]][[1]],
    "Error creating Protein S4 object: mock failure"
  )
  expect_identical(captured$notifications[[2]]$type, "error")
  expect_identical(captured$notifications[[2]]$duration, 15)
  expect_identical(captured$removed, "s4_creation_working")
  expect_identical(
    captured$error,
    "Error creating Protein S4 object: mock failure"
  )
  expect_identical(completed$status, "error")
  expect_identical(
    completed$errorMessage,
    "Error creating Protein S4 object: mock failure"
  )
})

test_that("mod_prot_qc_protein_s4_server wires the create observer through the helper seam", {
  workflow_data <- shiny::reactiveValues(state_manager = list(kind = "state_manager"))
  captured <- new.env(parent = emptyenv())
  captured$create <- NULL

  server_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_protein_s4_server,
    list(
      runProteinS4CreationObserver = function(workflowData, output) {
        captured$create <<- list(
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
      experiment_label = "TMT Experiment"
    ),
    {
      session$setInputs(create_protein_s4 = 1)
      session$flushReact()

      expect_false(is.null(captured$create))
      expect_identical(captured$create$workflowData, workflow_data)
      expect_identical(captured$create$output, output)
    }
  )
})

test_that("runProteinS4RevertStep reverts to the initial state and returns result text", {
  workflow_data <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$stateName <- NULL
  captured$info <- NULL
  workflow_data$state_manager <- list(
    revertToState = function(stateName) {
      captured$stateName <- stateName
      "reverted_s4"
    }
  )

  result <- runProteinS4RevertStep(
    workflowData = workflow_data,
    logInfoFn = function(message) {
      captured$info <- message
    }
  )

  expect_identical(captured$stateName, "initial")
  expect_identical(captured$info, "Reverted S4 object creation.")
  expect_identical(result$revertedState, "reverted_s4")
  expect_identical(result$stateName, "initial")
  expect_identical(
    result$resultText,
    "Reverted to initial empty state. You may need to re-run previous steps."
  )
})

test_that("runProteinS4RevertObserver delegates the revert workflow and success notification", {
  workflow_data <- new.env(parent = emptyenv())
  workflow_data_ref <- workflow_data
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()

  completed <- runProteinS4RevertObserver(
    workflowData = workflow_data,
    output = output,
    runRevertStepFn = function(workflowData) {
      captured$calls <- c(captured$calls, "run")
      expect_identical(workflowData, workflow_data_ref)
      list(
        revertedState = "reverted_s4",
        stateName = "initial",
        resultText = "Reverted to initial empty state. You may need to re-run previous steps."
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
    output$s4_creation_results,
    "Reverted to initial empty state. You may need to re-run previous steps."
  )
  expect_identical(captured$notifications[[1]][[1]], "Reverted successfully")
  expect_identical(captured$notifications[[1]]$type, "message")
  expect_identical(completed$status, "success")
  expect_identical(completed$revertResult$stateName, "initial")
  expect_identical(completed$revertResult$revertedState, "reverted_s4")
})

test_that("runProteinS4RevertObserver reports revert errors", {
  workflow_data <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  captured <- new.env(parent = emptyenv())
  captured$calls <- character()
  captured$notifications <- list()

  completed <- runProteinS4RevertObserver(
    workflowData = workflow_data,
    output = output,
    runRevertStepFn = function(...) {
      captured$calls <- c(captured$calls, "run")
      stop("mock revert failure")
    },
    renderTextFn = function(...) {
      stop("rendering should not run after a revert failure")
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
  expect_identical(captured$notifications[[1]][[1]], "Error reverting: mock revert failure")
  expect_identical(captured$notifications[[1]]$type, "error")
  expect_identical(completed$status, "error")
  expect_identical(completed$errorMessage, "Error reverting: mock revert failure")
})

test_that("mod_prot_qc_protein_s4_server wires the revert observer through the helper seam", {
  workflow_data <- shiny::reactiveValues(state_manager = list(kind = "state_manager"))
  captured <- new.env(parent = emptyenv())
  captured$revert <- NULL

  server_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_protein_s4_server,
    list(
      runProteinS4RevertObserver = function(workflowData, output) {
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
      experiment_label = "TMT Experiment"
    ),
    {
      session$setInputs(revert_s4_creation = 1)
      session$flushReact()

      expect_false(is.null(captured$revert))
      expect_identical(captured$revert$workflowData, workflow_data)
      expect_identical(captured$revert$output, output)
    }
  )
})

test_that("mod_prot_qc_protein_s4_ui preserves the public create contract", {
  html <- htmltools::renderTags(mod_prot_qc_protein_s4_ui("protein-s4"))$html

  expect_match(html, "Finalise Protein Data", fixed = TRUE)
  expect_match(html, "Create Protein S4 Object", fixed = TRUE)
  expect_match(html, "id=\"protein-s4-create_protein_s4\"", fixed = TRUE)
  expect_match(html, "id=\"protein-s4-s4_creation_results\"", fixed = TRUE)
})

test_that("mod_prot_qc_protein_s4_server keeps both helper seams wired through the thin wrapper", {
  workflow_data <- shiny::reactiveValues(state_manager = list(kind = "state_manager"))
  captured <- new.env(parent = emptyenv())
  captured$create <- NULL
  captured$revert <- NULL

  server_under_test <- makeFunctionWithOverrides(
    mod_prot_qc_protein_s4_server,
    list(
      runProteinS4CreationObserver = function(workflowData, output) {
        captured$create <<- list(
          workflowData = workflowData,
          output = output
        )
        invisible(NULL)
      },
      runProteinS4RevertObserver = function(workflowData, output) {
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
      experiment_label = "TMT Experiment"
    ),
    {
      session$setInputs(create_protein_s4 = 1)
      session$flushReact()
      session$setInputs(revert_s4_creation = 1)
      session$flushReact()

      expect_false(is.null(captured$create))
      expect_false(is.null(captured$revert))
      expect_identical(captured$create$workflowData, workflow_data)
      expect_identical(captured$revert$workflowData, workflow_data)
      expect_identical(captured$create$output, output)
      expect_identical(captured$revert$output, output)
    }
  )
})
