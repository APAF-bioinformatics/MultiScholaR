library(testthat)

repo_root <- normalizePath(file.path("..", ".."), mustWork = TRUE)

skipIfMissingProtNormSplitFiles <- function() {
  required_paths <- c(
    "R/mod_prot_norm_server_helpers.R",
    "R/mod_prot_norm_server.R"
  )
  missing <- required_paths[!file.exists(file.path(repo_root, required_paths))]
  if (length(missing) > 0) {
    testthat::skip(
      sprintf(
        "Target-only prot norm split file(s) not present: %s",
        paste(basename(missing), collapse = ", ")
      )
    )
  }
}

skipIfMissingProtNormSplitFiles()

loadProtNormModuleEnv <- function() {
  module_env <- new.env(parent = globalenv())
  sys.source(file.path(repo_root, "R", "mod_prot_norm_server_helpers.R"), envir = module_env)
  sys.source(file.path(repo_root, "R", "mod_prot_norm_server.R"), envir = module_env)
  module_env
}

test_that("runProtNormModuleServerShell preserves normalization entry-shell orchestration", {
  module_env <- loadProtNormModuleEnv()
  captured <- new.env(parent = emptyenv())
  output <- new.env(parent = emptyenv())
  norm_state <- new.env(parent = emptyenv())
  norm_state$qc_plot_paths <- list()
  norm_state$plot_refresh_trigger <- 0

  input <- list(
    color_variable = "batch",
    shape_variable = "group",
    ruv_grouping_variable = "group",
    ruv_mode = "automatic"
  )
  workflow_data <- list(
    state_manager = list(current_state = "post_filter")
  )
  experiment_paths <- list(protein_qc_dir = tempdir())

  result <- module_env$runProtNormModuleServerShell(
    input = input,
    output = output,
    session = "prot-session",
    id = "norm",
    workflowData = workflow_data,
    experimentPaths = experiment_paths,
    omicType = "proteomics",
    experimentLabel = "Proteomics",
    selectedTab = function() "norm",
    messageFn = function(...) invisible(NULL),
    createReactiveStateFn = function() {
      captured$create_reactive_state <- TRUE
      norm_state
    },
    getPlotAestheticsFn = function(colorVariable, shapeVariable) {
      captured$plot_aesthetics_request <- list(
        colorVariable = colorVariable,
        shapeVariable = shapeVariable
      )
      list(color_var = colorVariable, shape_var = shapeVariable)
    },
    getRuvGroupingVariableFn = function(groupingVariable) {
      captured$ruv_grouping_request <- groupingVariable
      paste0("resolved:", groupingVariable)
    },
    generatePreNormalizationQcArtifactsFn = function(...) {
      captured$pre_qc <- list(...)
      list(pre = "qc-paths")
    },
    generatePostNormalizationQcArtifactsFn = function(...) {
      captured$post_qc <- list(...)
      list(post = "qc-paths")
    },
    generateRuvCorrectedQcArtifactsFn = function(...) {
      captured$ruv_qc <- list(...)
      list(ruv = "qc-paths")
    },
    registerServerObserversFn = function(...) {
      captured$register_observers <- list(...)
      invisible(NULL)
    },
    registerQcImageOutputsFn = function(...) {
      captured$register_qc_images <- list(...)
      invisible(NULL)
    },
    registerRenderOutputsFn = function(...) {
      captured$register_render_outputs <- list(...)
      invisible(NULL)
    },
    checkMemoryUsageFn = function(...) "memory-ok"
  )

  expect_identical(result, norm_state)
  expect_true(isTRUE(captured$create_reactive_state))
  expect_identical(captured$register_qc_images$normData, norm_state)
  expect_identical(captured$register_qc_images$proteinQcDir, experiment_paths$protein_qc_dir)
  expect_identical(captured$register_render_outputs$normData, norm_state)
  expect_identical(captured$register_render_outputs$ruvMode, "automatic")
  expect_identical(captured$register_render_outputs$groupingVariable, "resolved:group")
  expect_identical(captured$ruv_grouping_request, "group")

  observers_args <- captured$register_observers
  expect_identical(observers_args$workflowData, workflow_data)
  expect_identical(observers_args$experimentPaths, experiment_paths)
  expect_identical(observers_args$omicType, "proteomics")
  expect_identical(observers_args$experimentLabel, "Proteomics")
  expect_true(is.function(observers_args$generatePreNormalizationQcFn))
  expect_true(is.function(observers_args$generatePostNormalizationQcFn))
  expect_true(is.function(observers_args$generateRuvCorrectedQcFn))
  expect_true(is.function(observers_args$getPlotAestheticsFn))
  expect_true(is.function(observers_args$getRuvGroupingVariableFn))

  expect_identical(
    observers_args$getPlotAestheticsFn(),
    list(color_var = "batch", shape_var = "group")
  )
  expect_identical(captured$plot_aesthetics_request$colorVariable, "batch")
  expect_identical(captured$plot_aesthetics_request$shapeVariable, "group")

  observers_args$generatePreNormalizationQcFn()
  expect_identical(captured$pre_qc$stateManager, workflow_data$state_manager)
  expect_identical(captured$pre_qc$qcDir, experiment_paths$protein_qc_dir)
  expect_identical(captured$pre_qc$aesthetics, list(color_var = "batch", shape_var = "group"))
  expect_identical(norm_state$qc_plot_paths, list(pre = "qc-paths"))

  observers_args$generatePostNormalizationQcFn("normalized-s4")
  expect_identical(captured$post_qc$normalizedS4, "normalized-s4")
  expect_identical(norm_state$qc_plot_paths, list(post = "qc-paths"))
  expect_identical(norm_state$plot_refresh_trigger, 1)

  observers_args$generateRuvCorrectedQcFn("ruv-s4")
  expect_identical(captured$ruv_qc$ruvCorrectedS4, "ruv-s4")
  expect_identical(norm_state$qc_plot_paths, list(ruv = "qc-paths"))
  expect_identical(norm_state$plot_refresh_trigger, 2)
})

test_that("runProtNormModuleServerEntryShell preserves moduleServer handoff", {
  module_env <- loadProtNormModuleEnv()
  captured <- new.env(parent = emptyenv())

  result <- module_env$runProtNormModuleServerEntryShell(
    id = "norm",
    workflowData = "workflow-state",
    experimentPaths = "experiment-paths",
    omicType = "proteomics",
    experimentLabel = "Proteomics",
    selectedTab = function() "norm",
    moduleServerFn = function(id, module) {
      captured$id <- id
      captured$module <- module
      "module-server-state"
    },
    runModuleServerShellFn = function(...) {
      captured$shell_args <- list(...)
      "module-shell-state"
    }
  )

  expect_identical(result, "module-server-state")
  expect_identical(captured$id, "norm")
  expect_true(is.function(captured$module))

  captured$module(
    input = list(),
    output = new.env(parent = emptyenv()),
    session = "prot-session"
  )

  expect_identical(captured$shell_args$id, "norm")
  expect_identical(captured$shell_args$workflowData, "workflow-state")
  expect_identical(captured$shell_args$experimentPaths, "experiment-paths")
  expect_identical(captured$shell_args$omicType, "proteomics")
  expect_identical(captured$shell_args$experimentLabel, "Proteomics")
  expect_identical(captured$shell_args$selectedTab(), "norm")
})

test_that("runProtNormModuleServerPublicWrapper preserves breadcrumb public-wrapper handoff", {
  module_env <- loadProtNormModuleEnv()
  captured <- new.env(parent = emptyenv())

  result <- module_env$runProtNormModuleServerPublicWrapper(
    id = "norm",
    workflow_data = "workflow-state",
    experiment_paths = "experiment-paths",
    omic_type = "proteomics",
    experiment_label = "Proteomics",
    selected_tab = function() "norm",
    runModuleServerEntryShellFn = function(...) {
      captured$args <- list(...)
      "entry-shell-state"
    }
  )

  expect_identical(result, "entry-shell-state")
  expect_identical(captured$args$id, "norm")
  expect_identical(captured$args$workflowData, "workflow-state")
  expect_identical(captured$args$experimentPaths, "experiment-paths")
  expect_identical(captured$args$omicType, "proteomics")
  expect_identical(captured$args$experimentLabel, "Proteomics")
  expect_identical(captured$args$selectedTab(), "norm")
})

test_that("mod_prot_norm_server routes through the top-level public seam", {
  module_env <- loadProtNormModuleEnv()
  captured <- new.env(parent = emptyenv())

  module_env$runProtNormModuleServerPublicWrapper <- function(...) {
    captured$args <- list(...)
    "public-wrapper-state"
  }

  result <- module_env$mod_prot_norm_server(
    id = "norm",
    workflow_data = "workflow-state",
    experiment_paths = "experiment-paths",
    omic_type = "proteomics",
    experiment_label = "Proteomics",
    selected_tab = function() "norm"
  )

  expect_identical(result, "public-wrapper-state")
  expect_identical(captured$args$id, "norm")
  expect_identical(captured$args$workflow_data, "workflow-state")
  expect_identical(captured$args$experiment_paths, "experiment-paths")
  expect_identical(captured$args$omic_type, "proteomics")
  expect_identical(captured$args$experiment_label, "Proteomics")
  expect_identical(captured$args$selected_tab(), "norm")
})
