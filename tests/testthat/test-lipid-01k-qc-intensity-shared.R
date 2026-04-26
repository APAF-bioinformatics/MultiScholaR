# fidelity-coverage-compare: shared
library(testthat)

getLipidIntensityServer <- function() {
  package_ns <- asNamespace("MultiScholaR")
  if (exists("mod_lipid_qc_intensity_server", envir = package_ns, inherits = FALSE)) {
    return(get("mod_lipid_qc_intensity_server", envir = package_ns, inherits = FALSE))
  }

  mod_lipid_qc_intensity_server
}

makeSharedLipidIntensityState <- function(lipid_data = NULL) {
  if (is.null(lipid_data)) {
    lipid_data <- list(
      AssayA = data.frame(
        lipid_id = c("L1", "L1", "L2"),
        S1 = c(12, 9, 7),
        S2 = c(15, 10, 8),
        stringsAsFactors = FALSE
      ),
      AssayB = data.frame(
        lipid_id = c("F1", "F2"),
        S1 = c(4, 2),
        S2 = c(5, 3),
        stringsAsFactors = FALSE
      )
    )
  }

  methods::new(
    "LipidomicsAssayData",
    lipid_data = lipid_data,
    lipid_id_column = "lipid_id",
    annotation_id_column = "annotation",
    database_identifier_type = "mock",
    internal_standard_regex = NA_character_,
    design_matrix = data.frame(
      Sample_ID = c("S1", "S2"),
      group = c("A", "B"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Sample_ID",
    group_id = "group",
    technical_replicate_id = NA_character_,
    args = list()
  )
}

makeSharedFilteredLipidIntensityState <- function() {
  makeSharedLipidIntensityState(
    lipid_data = list(
      AssayA = data.frame(
        lipid_id = "L1",
        S1 = 12,
        S2 = 15,
        stringsAsFactors = FALSE
      ),
      AssayB = data.frame(
        lipid_id = "F1",
        S1 = 4,
        S2 = 5,
        stringsAsFactors = FALSE
      )
    )
  )
}

newSharedLipidIntensityCapture <- function() {
  captured <- new.env(parent = emptyenv())
  captured$notifications <- list()
  captured$removed_notifications <- character()
  captured$log_info <- character()
  captured$log_warn <- character()
  captured$log_error <- character()
  captured
}

makeSharedLipidIntensityWorkflow <- function(current_state,
                                             captured,
                                             history = c("raw_lipid_s4", "lipid_intensity_filtered")) {
  state_manager <- list(
    getState = function() current_state,
    saveState = function(...) {
      captured$save_state <- list(...)
      invisible(NULL)
    },
    getHistory = function() history,
    revertToState = function(state_name) {
      captured$reverted_state <- state_name
      invisible(state_name)
    }
  )

  list(state_manager = state_manager, config_list = list(cache = "cfg"))
}

localSharedLipidIntensityNamespaceBinding <- function(env, name, value, .local_envir = parent.frame()) {
  had_binding <- exists(name, envir = env, inherits = FALSE)
  old_value <- if (had_binding) get(name, envir = env, inherits = FALSE) else NULL
  was_locked <- had_binding && bindingIsLocked(name, env)

  if (was_locked) {
    unlockBinding(name, env)
  }
  assign(name, value, envir = env)
  if (was_locked) {
    lockBinding(name, env)
  }

  withr::defer({
    if (exists(name, envir = env, inherits = FALSE) && bindingIsLocked(name, env)) {
      unlockBinding(name, env)
    }
    if (had_binding) {
      assign(name, old_value, envir = env)
    } else if (exists(name, envir = env, inherits = FALSE)) {
      rm(list = name, envir = env)
    }
    if (was_locked && exists(name, envir = env, inherits = FALSE)) {
      lockBinding(name, env)
    }
  }, envir = .local_envir)
}

withSharedLipidIntensityPackageMocks <- function(server_env,
                                                 captured,
                                                 filtered_state,
                                                 plot_result = structure(list(id = "qc-grid"), class = "grob"),
                                                 filter_error = NULL,
                                                 plot_error = NULL) {
  mock_frame <- parent.frame()

  localSharedLipidIntensityNamespaceBinding(
    server_env,
    "lipidIntensityFiltering",
    function(theObject,
             lipids_intensity_cutoff_percentile,
             lipids_proportion_of_samples_below_cutoff) {
      captured$filter_args <- list(
        theObject = theObject,
        percentile = lipids_intensity_cutoff_percentile,
        proportion = lipids_proportion_of_samples_below_cutoff
      )
      if (!is.null(filter_error)) {
        stop(filter_error)
      }
      filtered_state
    },
    mock_frame
  )
  localSharedLipidIntensityNamespaceBinding(
    server_env,
    "updateLipidFiltering",
    function(theObject, step_name, omics_type, return_grid, overwrite) {
      captured$plot_args <- list(
        theObject = theObject,
        step_name = step_name,
        omics_type = omics_type,
        return_grid = return_grid,
        overwrite = overwrite
      )
      if (!is.null(plot_error)) {
        stop(plot_error)
      }
      plot_result
    },
    mock_frame
  )
}

withSharedLipidIntensityUiMocks <- function(captured,
                                            input_values,
                                            evaluate_plot = TRUE) {
  mock_frame <- parent.frame()

  testthat::local_mocked_bindings(
    moduleServer = function(id, module) {
      output <- new.env(parent = emptyenv())
      module(
        input_values,
        output,
        list(ns = function(value) paste(id, value, sep = "-"))
      )
      captured$output <- output
      captured$module_id <- id
      invisible(NULL)
    },
    reactiveVal = function(value = NULL) {
      stored <- value
      function(new_value) {
        if (missing(new_value)) {
          stored
        } else {
          stored <<- new_value
          invisible(NULL)
        }
      }
    },
    observeEvent = function(eventExpr, handlerExpr, ...) {
      event_value <- eval(substitute(eventExpr), parent.frame())
      if (!is.null(event_value) && !identical(event_value, FALSE) && !identical(event_value, 0L)) {
        eval(substitute(handlerExpr), parent.frame())
      }
      invisible(NULL)
    },
    renderText = function(expr) {
      eval(substitute(expr), parent.frame())
    },
    renderUI = function(expr) {
      eval(substitute(expr), parent.frame())
    },
    renderPlot = function(expr) {
      if (isTRUE(evaluate_plot)) {
        eval(substitute(expr), parent.frame())
      }
      "rendered-plot"
    },
    req = function(...) {
      values <- list(...)
      for (value in values) {
        if (is.null(value) || identical(value, FALSE)) {
          stop("required value missing", call. = FALSE)
        }
      }
      values[[1]]
    },
    showNotification = function(message, ...) {
      captured$notifications[[length(captured$notifications) + 1L]] <- c(
        list(message = message),
        list(...)
      )
      invisible(NULL)
    },
    removeNotification = function(id, ...) {
      captured$removed_notifications <- c(captured$removed_notifications, id)
      invisible(NULL)
    },
    .package = "shiny",
    .env = mock_frame
  )

  testthat::local_mocked_bindings(
    log_info = function(message, ...) {
      captured$log_info <- c(captured$log_info, message)
      invisible(NULL)
    },
    log_warn = function(message, ...) {
      captured$log_warn <- c(captured$log_warn, message)
      invisible(NULL)
    },
    log_error = function(message, ...) {
      captured$log_error <- c(captured$log_error, message)
      invisible(NULL)
    },
    .package = "logger",
    .env = mock_frame
  )

  testthat::local_mocked_bindings(
    grid.draw = function(value, ...) {
      captured$drawn_plot <- value
      invisible(value)
    },
    .package = "grid",
    .env = mock_frame
  )
}

test_that("lipidomics intensity module preserves apply-filter success behavior", {
  captured <- newSharedLipidIntensityCapture()
  server_fn <- getLipidIntensityServer()
  server_env <- environment(server_fn)
  current_state <- makeSharedLipidIntensityState()
  filtered_state <- makeSharedFilteredLipidIntensityState()
  plot_result <- structure(list(id = "qc-grid"), class = "grob")

  withSharedLipidIntensityPackageMocks(server_env, captured, filtered_state, plot_result)
  withSharedLipidIntensityUiMocks(
    captured,
    input_values = list(
      apply_filter = TRUE,
      revert_filter = FALSE,
      intensity_cutoff_percentile = 10L,
      proportion_below_cutoff = 0.5
    )
  )

  workflow_data <- makeSharedLipidIntensityWorkflow(current_state, captured)
  server_fn(
    id = "intensity",
    workflow_data = workflow_data,
    omic_type = "lipidomics",
    experiment_label = "Lipid Experiment"
  )

  expect_identical(captured$module_id, "intensity")
  expect_identical(captured$filter_args$theObject, current_state)
  expect_identical(captured$filter_args$percentile, 10L)
  expect_identical(captured$filter_args$proportion, 0.5)
  expect_identical(captured$plot_args$theObject, filtered_state)
  expect_identical(captured$plot_args$step_name, "2_Intensity_Filtered")
  expect_identical(captured$plot_args$omics_type, "lipidomics")
  expect_true(captured$plot_args$return_grid)
  expect_true(captured$plot_args$overwrite)
  expect_identical(captured$save_state$state_name, "lipid_intensity_filtered")
  expect_identical(captured$save_state$s4_data_object, filtered_state)
  expect_identical(captured$save_state$config_object, list(cache = "cfg"))
  expect_match(captured$save_state$description, "percentile: 10%", fixed = TRUE)
  expect_match(captured$output$filter_results, "Lipid Intensity Filter Applied Successfully", fixed = TRUE)
  expect_match(captured$output$filter_results, "Total: 4 -> 2 lipids (removed 2)", fixed = TRUE)
  expect_match(
    as.character(htmltools::renderTags(captured$output$assay_results_tabs)$html),
    "intensity-assay_stats_tabs",
    fixed = TRUE
  )
  expect_identical(captured$output$filter_plot, "rendered-plot")
  expect_identical(captured$drawn_plot, plot_result)
  expect_identical(captured$removed_notifications, "lipid_intensity_filter_working")
  expect_true(any(grepl("Lipid intensity filter applied: removed 2 lipids", captured$log_info, fixed = TRUE)))
  expect_true(any(vapply(
    captured$notifications,
    function(notification) identical(notification$type, "message"),
    logical(1)
  )))
})

test_that("lipidomics intensity module preserves plot warning behavior", {
  captured <- newSharedLipidIntensityCapture()
  server_fn <- getLipidIntensityServer()
  server_env <- environment(server_fn)

  withSharedLipidIntensityPackageMocks(
    server_env,
    captured,
    filtered_state = makeSharedFilteredLipidIntensityState(),
    plot_error = "plot failed"
  )
  withSharedLipidIntensityUiMocks(
    captured,
    input_values = list(
      apply_filter = TRUE,
      revert_filter = FALSE,
      intensity_cutoff_percentile = 20L,
      proportion_below_cutoff = 0.6
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedLipidIntensityWorkflow(makeSharedLipidIntensityState(), captured)
  server_fn(
    id = "intensity",
    workflow_data = workflow_data,
    omic_type = "lipidomics",
    experiment_label = "Lipid Experiment"
  )

  expect_true(any(grepl("Could not generate QC plot: plot failed", captured$log_warn, fixed = TRUE)))
  expect_identical(captured$save_state$state_name, "lipid_intensity_filtered")
  expect_identical(captured$removed_notifications, "lipid_intensity_filter_working")
})

test_that("lipidomics intensity module preserves apply-filter error behavior", {
  captured <- newSharedLipidIntensityCapture()
  server_fn <- getLipidIntensityServer()
  server_env <- environment(server_fn)

  withSharedLipidIntensityPackageMocks(
    server_env,
    captured,
    filtered_state = makeSharedFilteredLipidIntensityState()
  )
  withSharedLipidIntensityUiMocks(
    captured,
    input_values = list(
      apply_filter = TRUE,
      revert_filter = FALSE,
      intensity_cutoff_percentile = 10L,
      proportion_below_cutoff = 0.5
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedLipidIntensityWorkflow(list(not = "s4"), captured)
  server_fn(
    id = "intensity",
    workflow_data = workflow_data,
    omic_type = "lipidomics",
    experiment_label = "Lipid Experiment"
  )

  expect_false(exists("save_state", envir = captured, inherits = FALSE))
  expect_false(exists("plot_args", envir = captured, inherits = FALSE))
  expect_true(any(grepl("Current state is not a LipidomicsAssayData object", captured$log_error, fixed = TRUE)))
  expect_true(any(vapply(
    captured$notifications,
    function(notification) {
      identical(notification$type, "error") &&
        grepl("Current state is not a LipidomicsAssayData object", notification$message, fixed = TRUE)
    },
    logical(1)
  )))
  expect_identical(captured$removed_notifications, "lipid_intensity_filter_working")
})

test_that("lipidomics intensity module preserves ggplot output dispatch", {
  captured <- newSharedLipidIntensityCapture()
  server_fn <- getLipidIntensityServer()
  server_env <- environment(server_fn)
  ggplot_result <- ggplot2::ggplot(
    data.frame(x = 1, y = 1),
    ggplot2::aes(x = x, y = y)
  ) + ggplot2::geom_point()

  withSharedLipidIntensityPackageMocks(
    server_env,
    captured,
    filtered_state = makeSharedFilteredLipidIntensityState(),
    plot_result = ggplot_result
  )
  withSharedLipidIntensityUiMocks(
    captured,
    input_values = list(
      apply_filter = TRUE,
      revert_filter = FALSE,
      intensity_cutoff_percentile = 10L,
      proportion_below_cutoff = 0.5
    )
  )

  workflow_data <- makeSharedLipidIntensityWorkflow(makeSharedLipidIntensityState(), captured)
  server_fn(
    id = "intensity",
    workflow_data = workflow_data,
    omic_type = "lipidomics",
    experiment_label = "Lipid Experiment"
  )

  expect_identical(captured$output$filter_plot, "rendered-plot")
})

test_that("lipidomics intensity module preserves revert success behavior", {
  captured <- newSharedLipidIntensityCapture()
  server_fn <- getLipidIntensityServer()

  withSharedLipidIntensityUiMocks(
    captured,
    input_values = list(
      apply_filter = FALSE,
      revert_filter = TRUE,
      intensity_cutoff_percentile = 10L,
      proportion_below_cutoff = 0.5
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedLipidIntensityWorkflow(
    makeSharedLipidIntensityState(),
    captured,
    history = c("raw_lipid_s4", "lipid_intensity_filtered", "lipid_qc")
  )
  server_fn(
    id = "intensity",
    workflow_data = workflow_data,
    omic_type = "lipidomics",
    experiment_label = "Lipid Experiment"
  )

  expect_identical(captured$reverted_state, "lipid_intensity_filtered")
  expect_identical(captured$output$filter_results, "Reverted to previous state: lipid_intensity_filtered")
  expect_true(any(grepl("Reverted lipid intensity filter to lipid_intensity_filtered", captured$log_info, fixed = TRUE)))
  expect_true(any(vapply(
    captured$notifications,
    function(notification) identical(notification$type, "message"),
    logical(1)
  )))
})

test_that("lipidomics intensity module preserves revert error behavior", {
  captured <- newSharedLipidIntensityCapture()
  server_fn <- getLipidIntensityServer()

  withSharedLipidIntensityUiMocks(
    captured,
    input_values = list(
      apply_filter = FALSE,
      revert_filter = TRUE,
      intensity_cutoff_percentile = 10L,
      proportion_below_cutoff = 0.5
    ),
    evaluate_plot = FALSE
  )

  workflow_data <- makeSharedLipidIntensityWorkflow(
    makeSharedLipidIntensityState(),
    captured,
    history = "raw_lipid_s4"
  )
  server_fn(
    id = "intensity",
    workflow_data = workflow_data,
    omic_type = "lipidomics",
    experiment_label = "Lipid Experiment"
  )

  expect_false(exists("reverted_state", envir = captured, inherits = FALSE))
  expect_true(any(grepl("No previous state to revert to.", captured$log_error, fixed = TRUE)))
  expect_true(any(vapply(
    captured$notifications,
    function(notification) {
      identical(notification$type, "error") &&
        grepl("No previous state to revert to.", notification$message, fixed = TRUE)
    },
    logical(1)
  )))
})
