# fidelity-coverage-compare: shared
library(testthat)

test_that("metabolomics QC variability helpers preserve grouped CV and internal standard metrics", {
  cv_assay_data <- data.frame(
    metabolite_id = c("M1", "M2"),
    Sample1 = c(10, 5),
    Sample2 = c(20, 10),
    Sample3 = c(30, 0),
    check.names = FALSE
  )
  design_matrix <- data.frame(
    Run = c("Sample1", "Sample2", "Sample3"),
    Group = c("A", "A", "B"),
    stringsAsFactors = FALSE
  )

  suppressMessages(
    cv_metrics <- calculateMetaboliteCVs(
      assay_data = cv_assay_data,
      design_matrix = design_matrix,
      group_id_col = "Group",
      replicate_id_col = NULL,
      sample_id_col = "Run",
      metabolite_id_col = "metabolite_id",
      sample_columns = c("Sample1", "Sample2", "Sample3")
    )
  )

  cv_metrics <- cv_metrics[order(cv_metrics$metabolite_id, cv_metrics$group), , drop = FALSE]
  expect_equal(cv_metrics$metabolite_id, c("M1", "M1", "M2"))
  expect_equal(cv_metrics$group, c("A", "B", "A"))
  expect_equal(cv_metrics$n_samples, c(2L, 1L, 2L))
  expect_equal(
    cv_metrics$cv[c(1, 3)],
    c(47.1404520791032, 47.1404520791032),
    tolerance = 1e-9
  )
  expect_true(is.na(cv_metrics$cv[2]))

  is_assay_data <- data.frame(
    metabolite_id = c("Met1", "IS_Caffeine", "Caffeine-d3"),
    Sample1 = c(10, 100, 50),
    Sample2 = c(20, 110, 70),
    Sample3 = c(30, 120, 90),
    check.names = FALSE
  )

  expect_message(
    is_metrics <- getInternalStandardMetrics(
      assay_data = is_assay_data,
      is_pattern = "^Z$",
      metabolite_id_col = "metabolite_id",
      sample_id_col = "Run",
      sample_columns = c("Sample1", "Sample2", "Sample3")
    ),
    "Using fallback patterns, found 2 candidates"
  )

  is_metrics <- is_metrics[order(is_metrics$is_id), , drop = FALSE]
  expect_equal(is_metrics$is_id, c("Caffeine-d3", "IS_Caffeine"))
  expect_equal(is_metrics$mean_intensity, c(70, 110))
  expect_equal(
    is_metrics$cv,
    c(28.5714285714286, 9.09090909090909),
    tolerance = 1e-9
  )
})

test_that("metabolomics QC variability helpers preserve current empty-result and guardrail branches", {
  cv_assay_data <- data.frame(
    metabolite_id = c("M1", "M2"),
    Sample1 = c(10, 5),
    Sample2 = c(20, 10),
    check.names = FALSE
  )

  expect_warning(
    missing_design_result <- calculateMetaboliteCVs(
      assay_data = cv_assay_data,
      design_matrix = data.frame(Run = c("Sample1", "Sample2")),
      group_id_col = "Group",
      replicate_id_col = NULL,
      sample_id_col = "Run",
      metabolite_id_col = "metabolite_id",
      sample_columns = c("Sample1", "Sample2")
    ),
    "Design matrix missing required columns"
  )
  expect_equal(
    missing_design_result,
    data.frame(metabolite_id = character(), group = character(), cv = numeric())
  )

  expect_warning(
    no_quant_result <- calculateMetaboliteCVs(
      assay_data = data.frame(metabolite_id = c("M1", "M2"), annotation = c("a", "b")),
      design_matrix = data.frame(Run = c("Sample1", "Sample2"), Group = c("A", "A")),
      group_id_col = "Group",
      replicate_id_col = NULL,
      sample_id_col = "Run",
      metabolite_id_col = "metabolite_id",
      sample_columns = character()
    ),
    "No sample columns or no data rows found for CV calculation"
  )
  expect_equal(
    no_quant_result,
    data.frame(metabolite_id = character(), group = character(), cv = numeric())
  )

  expect_equal(
    getInternalStandardMetrics(
      assay_data = cv_assay_data,
      is_pattern = "",
      metabolite_id_col = "metabolite_id",
      sample_id_col = "Run",
      sample_columns = c("Sample1", "Sample2")
    ),
    data.frame(is_id = character(), mean_intensity = numeric(), cv = numeric())
  )

  expect_warning(
    missing_id_result <- getInternalStandardMetrics(
      assay_data = data.frame(Sample1 = 1, Sample2 = 2),
      is_pattern = "IS",
      metabolite_id_col = "metabolite_id",
      sample_id_col = "Run",
      sample_columns = c("Sample1", "Sample2")
    ),
    "Metabolite ID column"
  )
  expect_equal(
    missing_id_result,
    data.frame(is_id = character(), mean_intensity = numeric(), cv = numeric())
  )
})
