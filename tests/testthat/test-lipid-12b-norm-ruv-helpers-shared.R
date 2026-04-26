# fidelity-coverage-compare: shared
library(testthat)

makeSharedLipidRuvObject <- function(include_negative = TRUE) {
  lipid_data <- list(
    `Positive Mode` = data.frame(
      database_identifier = c("L1", "L2"),
      lipid_identification = c("A1", "A2"),
      S1 = c(10, 20),
      S2 = c(30, 40),
      stringsAsFactors = FALSE
    )
  )

  if (isTRUE(include_negative)) {
    lipid_data$`Negative Mode` <- data.frame(
      database_identifier = c("L3", "L4"),
      lipid_identification = c("B1", "B2"),
      S1 = c(15, 25),
      S2 = c(35, 45),
      stringsAsFactors = FALSE
    )
  }

  methods::new(
    "LipidomicsAssayData",
    lipid_data = lipid_data,
    lipid_id_column = "database_identifier",
    annotation_id_column = "lipid_identification",
    database_identifier_type = "MockDB",
    internal_standard_regex = "",
    design_matrix = data.frame(
      Sample_ID = c("S1", "S2"),
      group = c("A", "B"),
      batch = c("B1", "B2"),
      stringsAsFactors = FALSE
    ),
    sample_id = "Sample_ID",
    group_id = "group",
    technical_replicate_id = NA_character_,
    args = list()
  )
}

localSharedLipidNormBinding <- function(env, name, value, .local_envir = parent.frame()) {
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

test_that("runLipidPerAssayRuvOptimization covers automatic mode and extractor helpers", {
  lipid_object <- makeSharedLipidRuvObject()
  package_ns <- asNamespace("MultiScholaR")
  captured <- new.env(parent = emptyenv())
  captured$ctrl_calls <- list()
  captured$cancor_calls <- list()

  localSharedLipidNormBinding(
    package_ns,
    "getNegCtrlMetabAnova",
    function(theObject, ruv_grouping_variable, percentage_as_neg_ctrl) {
      captured$ctrl_calls[[length(captured$ctrl_calls) + 1L]] <<- list(
        grouping = ruv_grouping_variable,
        percentage = percentage_as_neg_ctrl
      )
      list(
        `Positive Mode` = c(TRUE, FALSE, TRUE),
        `Negative Mode` = c(FALSE, TRUE, TRUE)
      )
    }
  )
  localSharedLipidNormBinding(
    package_ns,
    "ruvCancor",
    function(theObject, ctrl, ruv_grouping_variable) {
      captured$cancor_calls[[length(captured$cancor_calls) + 1L]] <<- list(
        ctrl = ctrl,
        grouping = ruv_grouping_variable
      )
      list(
        `Positive Mode` = list(best_k = 2L, score = 0.91),
        `Negative Mode` = list(best_k = 3L, score = 0.87)
      )
    }
  )
  localSharedLipidNormBinding(package_ns, "findBestK", function(cancor_plot) cancor_plot$best_k)
  localSharedLipidNormBinding(
    package_ns,
    "calculateSeparationScore",
    function(cancor_plot, metric) {
      expect_identical(metric, "max_difference")
      cancor_plot$score
    }
  )

  results <- runLipidPerAssayRuvOptimization(
    lipid_object,
    ruv_mode = "automatic",
    params = list(
      percentage_max = 15,
      ruv_grouping_variable = "missing_group"
    )
  )

  expect_identical(names(results), c("Positive Mode", "Negative Mode"))
  expect_true(all(vapply(results, function(x) isTRUE(x$success), logical(1))))
  expect_identical(results$`Positive Mode`$best_k, 2L)
  expect_identical(results$`Negative Mode`$best_k, 3L)
  expect_identical(results$`Positive Mode`$best_percentage, 15)
  expect_identical(results$`Negative Mode`$best_percentage, 15)
  expect_identical(results$`Positive Mode`$control_genes_index, c(TRUE, FALSE, TRUE))
  expect_identical(results$`Negative Mode`$control_genes_index, c(FALSE, TRUE, TRUE))
  expect_identical(captured$ctrl_calls[[1]]$grouping, "group")
  expect_identical(captured$cancor_calls[[1]]$grouping, "group")

  best_k <- extractLipidBestKPerAssay(results)
  ctrl_idx <- extractLipidCtrlPerAssay(results)
  combined <- buildLipidCombinedRuvTable(results)

  expect_identical(best_k$`Positive Mode`, 2L)
  expect_identical(best_k$`Negative Mode`, 3L)
  expect_identical(ctrl_idx$`Positive Mode`, c(TRUE, FALSE, TRUE))
  expect_identical(ctrl_idx$`Negative Mode`, c(FALSE, TRUE, TRUE))
  expect_identical(combined$Status, c("Success", "Success"))
  expect_identical(combined$Num_Controls, c(2L, 2L))
})

test_that("runLipidPerAssayRuvOptimization covers manual mode", {
  lipid_object <- makeSharedLipidRuvObject(include_negative = FALSE)
  package_ns <- asNamespace("MultiScholaR")
  captured <- new.env(parent = emptyenv())

  localSharedLipidNormBinding(
    package_ns,
    "getNegCtrlMetabAnova",
    function(theObject, ruv_grouping_variable, percentage_as_neg_ctrl) {
      captured$percentage <<- percentage_as_neg_ctrl
      captured$grouping <<- ruv_grouping_variable
      c(TRUE, FALSE, TRUE, FALSE)
    }
  )
  localSharedLipidNormBinding(
    package_ns,
    "ruvCancor",
    function(theObject, ctrl, ruv_grouping_variable) {
      captured$ctrl <<- ctrl
      list(best_k = 99L, score = 0.5)
    }
  )

  results <- runLipidPerAssayRuvOptimization(
    lipid_object,
    ruv_mode = "manual",
    params = list(
      ruv_grouping_variable = "batch",
      manual_k = 4L,
      manual_percentage = 12
    )
  )

  expect_identical(results$`Positive Mode`$best_k, 4L)
  expect_identical(results$`Positive Mode`$best_percentage, 12)
  expect_identical(results$`Positive Mode`$separation_score, NA_real_)
  expect_null(results$`Positive Mode`$optimization_results)
  expect_identical(captured$percentage, 12)
  expect_identical(captured$grouping, "batch")
  expect_identical(captured$ctrl, c(TRUE, FALSE, TRUE, FALSE))
})

test_that("runLipidPerAssayRuvOptimization records failures and table helpers surface them", {
  lipid_object <- makeSharedLipidRuvObject(include_negative = FALSE)
  package_ns <- asNamespace("MultiScholaR")

  localSharedLipidNormBinding(
    package_ns,
    "getNegCtrlMetabAnova",
    function(...) {
      stop("mock automatic failure")
    }
  )

  results <- runLipidPerAssayRuvOptimization(
    lipid_object,
    ruv_mode = "automatic",
    params = list(ruv_grouping_variable = "group")
  )

  expect_false(results$`Positive Mode`$success)
  expect_match(results$`Positive Mode`$error, "mock automatic failure", fixed = TRUE)
  expect_true(is.na(extractLipidBestKPerAssay(results)$`Positive Mode`))
  expect_null(extractLipidCtrlPerAssay(results)$`Positive Mode`)

  combined <- buildLipidCombinedRuvTable(results)
  expect_match(combined$Status[[1]], "Failed: mock automatic failure", fixed = TRUE)
  expect_true(is.na(combined$Best_K[[1]]))
})
