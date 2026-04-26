# fidelity-coverage-compare: shared
library(testthat)

localBinding <- function(env, name, value, .local_envir = parent.frame()) {
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

makeCancorPlot <- function(differences = c(0.1, 0.5, 0.2)) {
  plot_data <- data.frame(
    featureset = rep(c("Control", "All"), each = 3),
    cc = c(0.2, 0.3, 0.4, 0.2 + differences),
    K = rep(1:3, times = 2)
  )
  ggplot2::ggplot(plot_data, ggplot2::aes(K, cc, colour = featureset)) +
    ggplot2::geom_point()
}

makeProteinRuvObject <- function() {
  methods::new(
    "ProteinQuantitativeData",
    protein_quant_table = data.frame(
      Protein.Ids = c("P1", "P2", "P3"),
      S1 = c(1, 2, 3),
      S2 = c(2, 3, 4),
      S3 = c(3, 4, 5),
      S4 = c(4, 5, 6),
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    protein_id_column = "Protein.Ids",
    protein_id_table = data.frame(Protein.Ids = c("P1", "P2", "P3"), stringsAsFactors = FALSE),
    design_matrix = data.frame(
      Run = c("S1", "S2", "S3", "S4"),
      group = c("A", "A", "B", "B"),
      replicate = c("R1", "R2", "R1", "R2"),
      stringsAsFactors = FALSE,
      row.names = c("S1", "S2", "S3", "S4")
    ),
    sample_id = "Run",
    group_id = "group",
    technical_replicate_id = "replicate",
    args = list()
  )
}

test_that("protein RUV helper functions preserve K and scoring calculations", {
  plot <- makeCancorPlot()

  expect_identical(findBestK(plot), 2L)

  scored_list <- suppressWarnings(findBestKForAssayList(list(Assay1 = plot, Assay2 = 42)))
  unnamed_list <- suppressWarnings(findBestKForAssayList(list(plot)))

  expect_identical(scored_list$Assay1, 2L)
  expect_true(is.na(scored_list$Assay2))
  expect_identical(names(unnamed_list), "Plot_1")

  expect_equal(calculateSeparationScore(plot, "max_difference"), 0.4)
  expect_equal(calculateSeparationScore(plot, "mean_difference"), mean(c(0.1, 0.4, 0)))
  expect_gt(calculateSeparationScore(plot, "auc"), 0)
  expect_gt(calculateSeparationScore(plot, "weighted_difference"), 0)
  expect_true(is.na(calculateSeparationScore(plot, "not-a-metric")))
  expect_true(is.na(calculateSeparationScore(list(), "max_difference")))

  expect_equal(calculateCompositeScore(1, 1, 0.5, 3), 1)
  expect_lt(calculateCompositeScore(1, 5, 0.5, 3), 1)
  expect_identical(calculateCompositeScore(NA_real_, 2, 0.5, 3), NA_real_)
  expect_identical(calculateCompositeScore(-1, 2, 0.5, 3), 0)

  expect_identical(calculateAdaptiveMaxK(10), 2L)
  expect_identical(calculateAdaptiveMaxK(20), 3L)
  expect_identical(calculateAdaptiveMaxK(60), 4L)
  expect_identical(calculateAdaptiveMaxK(100), 5L)
})

test_that("protein RUV helper utilities preserve parameter, replicate, and scaling outputs", {
  updated <- updateRuvParameters(
    config_list = list(ruvParameters = list()),
    best_k = 3L,
    control_genes_index = c(TRUE, FALSE, TRUE),
    percentage_as_neg_ctrl = 10
  )

  replicate_matrix <- getRuvIIIReplicateMatrixHelper(
    design_matrix = data.frame(
      Run = c("S1", "S2", "S3"),
      group = c("A", "A", "B"),
      stringsAsFactors = FALSE
    ),
    sample_id_column = Run,
    grouping_variable = group
  )

  control_index <- getNegCtrlProtAnovaHelper(
    data_matrix = structure(
      c(
        10, 11, 12, 13, 14,
        20, 21, 22, 23, 24,
        30, 31, 32, 33, 34
      ),
      dim = c(3, 5),
      dimnames = list(
        c("P1", "P2", "P3"),
        c("Pool1", "A1", "A2", "B1", "B2")
      )
    ),
    design_matrix = data.frame(
      group = c("Pool", "A", "A", "B", "B"),
      stringsAsFactors = FALSE,
      row.names = c("Pool1", "A1", "A2", "B1", "B2")
    ),
    grouping_variable = "group",
    percentage_as_neg_ctrl = 50,
    ruv_qval_cutoff = 0.05,
    ruv_fdr_method = "BH",
    exclude_pool_samples = TRUE
  )

  extracted <- extractRuvResults(list(
    one = list(results = data.frame(value = 1)),
    two = list(results = data.frame(value = 2))
  ))

  scaled <- scaleCenterAndFillMissing(matrix(
    c(1, 2, NA, 4),
    nrow = 2,
    dimnames = list(c("P1", "P2"), c("S1", "S2"))
  ))

  expect_identical(updated$ruvParameters$best_k, 3L)
  expect_identical(updated$ruvParameters$num_neg_ctrl, 3L)
  expect_identical(updated$ruvParameters$percentage_as_neg_ctrl, 10)
  expect_identical(colnames(replicate_matrix), c("A", "B"))
  expect_equal(replicate_matrix["S1", "A"], 1)
  expect_identical(length(control_index), 3L)
  expect_identical(names(control_index), c("P1", "P2", "P3"))
  expect_identical(names(extracted), c("one", "two"))
  expect_false(any(is.na(scaled)))
})

test_that("protein RUV percentage optimization preserves ranking and validation branches", {
  helper_env <- environment(findBestNegCtrlPercentage)
  localBinding(helper_env, "getNegCtrlProtAnova", function(
    normalised_protein_matrix_obj,
    ruv_grouping_variable,
    percentage_as_neg_ctrl,
    ruv_qval_cutoff,
    ruv_fdr_method
  ) {
    structure(
      c(
        rep(TRUE, if (percentage_as_neg_ctrl == 1) 3 else if (percentage_as_neg_ctrl == 5) 6 else 8),
        rep(FALSE, if (percentage_as_neg_ctrl == 1) 7 else if (percentage_as_neg_ctrl == 5) 4 else 2)
      ),
      percentage = percentage_as_neg_ctrl
    )
  })
  localBinding(helper_env, "ruvCancor", function(
    normalised_protein_matrix_obj,
    ctrl,
    num_components_to_impute,
    ruv_grouping_variable
  ) {
    list(percentage = attr(ctrl, "percentage"))
  })
  localBinding(helper_env, "calculateSeparationScore", function(cancorplot, metric = "max_difference") {
    cancorplot$percentage / 10
  })
  localBinding(helper_env, "findBestK", function(cancorplot) {
    if (identical(cancorplot$percentage, 5)) 2L else 3L
  })

  result <- findBestNegCtrlPercentage(
    normalised_protein_matrix_obj = makeProteinRuvObject(),
    percentage_range = c(1, 5, 10),
    adaptive_k_penalty = TRUE,
    verbose = FALSE
  )

  expect_identical(result$best_percentage, 5)
  expect_identical(result$best_k, 2L)
  expect_identical(sum(result$best_control_genes_index), 6L)
  expect_equal(result$best_separation_score, 0.5)
  expect_gt(result$best_composite_score, 0)
  expect_identical(result$max_acceptable_k, 2L)
  expect_true(is.data.frame(result$optimization_results))
  expect_identical(result$sample_size, 4)

  expect_error(
    findBestNegCtrlPercentage(
      normalised_protein_matrix_obj = list(),
      percentage_range = 1:2,
      verbose = FALSE
    ),
    "normalised_protein_matrix_obj must be a ProteinQuantitativeData object",
    fixed = TRUE
  )

  expect_error(
    findBestNegCtrlPercentage(
      normalised_protein_matrix_obj = makeProteinRuvObject(),
      percentage_range = c(0, 5),
      verbose = FALSE
    ),
    "percentage_range must contain values between 0 and 100",
    fixed = TRUE
  )
})
