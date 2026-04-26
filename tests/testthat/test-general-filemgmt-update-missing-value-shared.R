# fidelity-coverage-compare: shared
library(testthat)

getSharedUpdateMissingValueParameters <- function() {
  package_ns <- asNamespace("MultiScholaR")
  if (exists("updateMissingValueParameters", envir = package_ns, inherits = FALSE)) {
    return(get("updateMissingValueParameters", envir = package_ns, inherits = FALSE))
  }

  updateMissingValueParameters
}

if (!methods::isClass("SharedMissingValueParameterObject")) {
  methods::setClass(
    "SharedMissingValueParameterObject",
    slots = c(design_matrix = "data.frame", args = "list")
  )
}

if (!methods::isClass("SharedMissingValueNoArgsObject")) {
  methods::setClass(
    "SharedMissingValueNoArgsObject",
    slots = c(design_matrix = "data.frame")
  )
}

makeSharedMissingValueObject <- function(design_matrix, args = list()) {
  methods::new(
    "SharedMissingValueParameterObject",
    design_matrix = design_matrix,
    args = args
  )
}

makeSharedMissingValueDesign <- function() {
  data.frame(
    sample = paste0("sample_", seq_len(5)),
    group = c("A", "A", "B", "B", "B"),
    stringsAsFactors = FALSE
  )
}

captureSharedMissingValueMessages <- function(expr) {
  messages <- character()
  value <- withCallingHandlers(
    expr,
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  list(value = value, messages = messages)
}

runSharedUpdateMissingValueParameters <- function(fun,
                                                  design_matrix = makeSharedMissingValueDesign(),
                                                  config_list = list(),
                                                  args = list(),
                                                  min_reps_per_group = 2,
                                                  min_groups = 2,
                                                  ...) {
  contract_env <- new.env(parent = emptyenv())
  contract_env$config_list <- config_list

  captured <- captureSharedMissingValueMessages(
    fun(
      theObject = makeSharedMissingValueObject(design_matrix, args = args),
      min_reps_per_group = min_reps_per_group,
      min_groups = min_groups,
      config_list_name = "config_list",
      env = contract_env,
      ...
    )
  )

  list(
    object = captured$value,
    messages = captured$messages,
    config_list = contract_env$config_list
  )
}

test_that("updateMissingValueParameters synchronizes default missingness thresholds", {
  update_fn <- getSharedUpdateMissingValueParameters()
  result <- runSharedUpdateMissingValueParameters(
    update_fn,
    config_list = list(removeRowsWithMissingValuesPercent = list(existing = TRUE))
  )

  config_params <- result$config_list$removeRowsWithMissingValuesPercent
  object_params <- result$object@args$removeRowsWithMissingValuesPercent

  expect_true(config_params$existing)
  expect_equal(config_params$groupwise_percentage_cutoff, 33.333)
  expect_equal(config_params$max_groups_percentage_cutoff, 0)
  expect_equal(config_params$proteins_intensity_cutoff_percentile, 1)
  expect_identical(object_params, config_params[names(object_params)])
  expect_true(any(grepl("Updated missing value parameters", result$messages, fixed = TRUE)))
  expect_true(any(grepl("Group A: 2 out of 2 replicates required", result$messages, fixed = TRUE)))
  expect_true(any(grepl("Group B: 2 out of 3 replicates required", result$messages, fixed = TRUE)))
  expect_true(any(grepl("S4 object @args and global config_list", result$messages, fixed = TRUE)))
})

test_that("updateMissingValueParameters initializes absent config and args entries", {
  update_fn <- getSharedUpdateMissingValueParameters()
  result <- runSharedUpdateMissingValueParameters(
    update_fn,
    config_list = list(),
    args = list(removeRowsWithMissingValuesPercent = list(existing_arg = TRUE))
  )

  expect_equal(
    result$config_list$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff,
    33.333
  )
  expect_equal(
    result$object@args$removeRowsWithMissingValuesPercent$groupwise_percentage_cutoff,
    33.333
  )
  expect_true(result$object@args$removeRowsWithMissingValuesPercent$existing_arg)
})

test_that("updateMissingValueParameters preserves validation failures", {
  update_fn <- getSharedUpdateMissingValueParameters()
  design_matrix <- makeSharedMissingValueDesign()
  contract_env <- new.env(parent = emptyenv())
  contract_env$config_list <- list()

  expect_error(
    update_fn(
      theObject = list(),
      config_list_name = "config_list",
      env = contract_env
    ),
    "'theObject' must be an S4 object.",
    fixed = TRUE
  )

  expect_error(
    update_fn(
      theObject = methods::new("SharedMissingValueNoArgsObject", design_matrix = design_matrix),
      config_list_name = "config_list",
      env = contract_env
    ),
    "@args",
    fixed = TRUE
  )

  expect_error(
    update_fn(
      theObject = makeSharedMissingValueObject(design_matrix),
      config_list_name = "missing_config",
      env = contract_env
    ),
    "Global config list 'missing_config' not found in the specified environment.",
    fixed = TRUE
  )

  expect_error(
    runSharedUpdateMissingValueParameters(update_fn, min_groups = 3),
    "min_groups cannot be larger than total number of groups",
    fixed = TRUE
  )
})
