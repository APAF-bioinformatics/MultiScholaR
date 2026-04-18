library(testthat)

repo_root <- normalizePath(file.path("..", ".."), mustWork = TRUE)

loadSelectedFunctions <- function(paths, symbols, env) {
  for (path in paths[file.exists(paths)]) {
    exprs <- parse(file = path, keep.source = TRUE)

    for (expr in exprs) {
      is_assignment <- is.call(expr) &&
        length(expr) >= 3 &&
        as.character(expr[[1]]) %in% c("<-", "=")

      if (!is_assignment || !is.symbol(expr[[2]])) {
        next
      }

      symbol_name <- as.character(expr[[2]])
      if (symbol_name %in% symbols) {
        eval(expr, envir = env)
      }
    }
  }
}

loadSelectedFunctions(
  paths = c(
    file.path(repo_root, "R", "func_general_helpers.R"),
    file.path(repo_root, "R", "func_metab_norm_support_helpers.R"),
    file.path(repo_root, "R", "func_metab_norm.R")
  ),
  symbols = c(
    "%||%",
    "extractBestKPerAssay",
    "extractCtrlPerAssay",
    "buildCombinedRuvTable",
    "buildNormConfig",
    "buildItsdSelectionTable"
  ),
  env = environment()
)

test_that("metabolomics normalization support helpers preserve the current RUV summary contracts", {
  ruv_results <- list(
    Plasma = list(
      success = TRUE,
      best_k = 3L,
      best_percentage = 12.5,
      control_genes_index = c(TRUE, FALSE, TRUE, FALSE),
      separation_score = 0.45678
    ),
    Urine = list(
      success = FALSE,
      best_k = NA_integer_,
      best_percentage = NA_real_,
      control_genes_index = NULL,
      separation_score = NA_real_,
      error = "no controls available"
    )
  )

  best_k <- extractBestKPerAssay(ruv_results)
  ctrl <- extractCtrlPerAssay(ruv_results)
  combined <- buildCombinedRuvTable(ruv_results)

  expect_identical(best_k, list(Plasma = 3L, Urine = NA_integer_))
  expect_identical(ctrl$Plasma, c(TRUE, FALSE, TRUE, FALSE))
  expect_null(ctrl$Urine)

  expect_identical(combined$Assay, c("Plasma", "Urine"))
  expect_identical(combined$Best_K, c(3L, NA_integer_))
  expect_equal(combined$Best_Percentage, c(12.5, NA_real_))
  expect_equal(combined$Separation_Score, c(0.4568, NA_real_))
  expect_identical(combined$Num_Controls, c(2L, NA_integer_))
  expect_identical(combined$Status[[1]], "Success")
  expect_identical(combined$Status[[2]], "Failed: no controls available")
})

test_that("metabolomics normalization support helpers preserve current normalization config defaults", {
  default_config <- buildNormConfig(list())

  expect_false(default_config$itsd$applied)
  expect_identical(default_config$itsd$method, "median")
  expect_identical(default_config$log2$offset, 1)
  expect_identical(default_config$normalization$method, "cyclicloess")
  expect_identical(default_config$ruv$mode, "skip")
  expect_null(default_config$ruv$grouping_variable)
  expect_null(default_config$ruv$auto_percentage_min)
  expect_null(default_config$ruv$manual_k)
  expect_s3_class(default_config$timestamp, "POSIXct")

  explicit_config <- buildNormConfig(list(
    apply_itsd = TRUE,
    itsd_method = "mean",
    log_offset = 5,
    norm_method = "quantile",
    ruv_mode = "manual",
    ruv_grouping_variable = "Condition",
    auto_percentage_min = 2,
    auto_percentage_max = 15,
    separation_metric = "mean_difference",
    k_penalty_weight = 0.75,
    adaptive_k_penalty = FALSE,
    ruv_k = 4,
    ruv_percentage = 20
  ))

  expect_true(explicit_config$itsd$applied)
  expect_identical(explicit_config$itsd$method, "mean")
  expect_identical(explicit_config$log2$offset, 5)
  expect_identical(explicit_config$normalization$method, "quantile")
  expect_identical(explicit_config$ruv$mode, "manual")
  expect_identical(explicit_config$ruv$grouping_variable, "Condition")
  expect_identical(explicit_config$ruv$auto_percentage_min, 2)
  expect_identical(explicit_config$ruv$auto_percentage_max, 15)
  expect_identical(explicit_config$ruv$separation_metric, "mean_difference")
  expect_identical(explicit_config$ruv$k_penalty_weight, 0.75)
  expect_identical(explicit_config$ruv$adaptive_k_penalty, FALSE)
  expect_identical(explicit_config$ruv$manual_k, 4)
  expect_identical(explicit_config$ruv$manual_percentage, 20)
  expect_s3_class(explicit_config$timestamp, "POSIXct")
})

test_that("metabolomics normalization support helpers preserve current ITSD selection table contracts", {
  assay_data <- data.frame(
    feature_code = c("M1", "M2", "M3"),
    metabolite = c("Lactic acid", "ITSD-d5", "Unknown"),
    annotation = c("reference", "internal standard", "rt_peak"),
    sample_a = c(10, 100, 50),
    sample_b = c(12, 98, 55),
    batch = c("B1", "B1", "B2"),
    stringsAsFactors = FALSE
  )

  selection_table <- buildItsdSelectionTable(
    assay_data = assay_data,
    metabolite_id_col = "feature_code",
    annotation_cols = c("metabolite", "annotation")
  )

  expect_identical(
    names(selection_table),
    c("feature_id", "annotation", "mean_intensity", "cv_percent", "is_candidate")
  )
  expect_identical(selection_table$feature_id, c("M2", "M3", "M1"))
  expect_identical(
    selection_table$annotation,
    c("ITSD-d5 | internal standard", "Unknown | rt_peak", "Lactic acid | reference")
  )
  expect_equal(selection_table$mean_intensity, c(99, 52.5, 11))
  expect_equal(
    selection_table$cv_percent,
    c(1.428499, 6.73435, 12.85649),
    tolerance = 1e-5
  )
  expect_identical(selection_table$is_candidate, c(TRUE, FALSE, FALSE))
})

test_that("metabolomics normalization support helpers preserve empty-sample ITSD fallback", {
  assay_data <- data.frame(
    feature_code = c("M1", "M2"),
    metabolite = c("A", "B"),
    batch = c("B1", "B2"),
    stringsAsFactors = FALSE
  )

  expect_warning(
    selection_table <- buildItsdSelectionTable(
      assay_data = assay_data,
      metabolite_id_col = "feature_code"
    ),
    "No sample columns found in assay_data"
  )

  expect_identical(
    selection_table,
    data.frame(
      feature_id = character(0),
      annotation = character(0),
      mean_intensity = numeric(0),
      cv_percent = numeric(0),
      is_candidate = logical(0)
    )
  )
})
