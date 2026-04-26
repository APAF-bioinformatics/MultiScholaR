library(testthat)

# fidelity-coverage-compare: shared

newLipidLogTransformObject <- function(assay_data = NULL, args = list()) {
  if (is.null(assay_data)) {
    assay_data <- tibble::tibble(
      Name = c("L1", "L2"),
      Annotation = c("Alpha", "Beta"),
      Class = c("C1", "C2"),
      S1 = c(0, -3),
      S2 = c(6, 14)
    )
  }

  createLipidomicsAssayData(
    lipid_data = list(LCMS_Pos = assay_data),
    design_matrix = data.frame(
      Run = c("S1", "S2"),
      group = c("A", "B"),
      stringsAsFactors = FALSE
    ),
    lipid_id_column = "Name",
    annotation_id_column = "Annotation",
    sample_id = "Run",
    group_id = "group",
    args = args
  )
}

test_that("lipidomics S4 log-transform method preserves current assay transformation contract", {
  transformed <- logTransformAssays(newLipidLogTransformObject(), offset = 2)
  transformed_assay <- transformed@lipid_data$LCMS_Pos

  expect_s3_class(transformed_assay, "tbl_df")
  expect_equal(transformed_assay$S1, c(1, 1))
  expect_equal(transformed_assay$S2, log2(c(8, 16)))
  expect_identical(transformed_assay$Name, c("L1", "L2"))
  expect_identical(transformed_assay$Annotation, c("Alpha", "Beta"))
  expect_identical(transformed_assay$Class, c("C1", "C2"))
  expect_true(transformed@args$log_transformed)
  expect_identical(transformed@args$log_transform_offset, 2)
})

test_that("lipidomics S4 log-transform method rejects non-positive offsets", {
  expect_error(
    logTransformAssays(newLipidLogTransformObject(), offset = 0),
    "`offset` must be a single positive numeric value.",
    fixed = TRUE
  )
})

test_that("lipidomics S4 log-transform method preserves empty assay objects", {
  lipid_object <- createLipidomicsAssayData(
    lipid_data = list(),
    design_matrix = data.frame(
      Run = c("S1", "S2"),
      group = c("A", "B"),
      stringsAsFactors = FALSE
    ),
    lipid_id_column = "Name",
    annotation_id_column = "Annotation",
    sample_id = "Run",
    group_id = "group",
    args = list(existing = TRUE)
  )

  transformed <- suppressWarnings(logTransformAssays(lipid_object, offset = 1))

  expect_length(transformed@lipid_data, 0)
  expect_identical(transformed@args, list(existing = TRUE))
})

test_that("lipidomics S4 log-transform method handles unnamed non-tibble assays", {
  assay_data <- data.frame(
    Name = c("L1", "L2"),
    Annotation = c("Alpha", "Beta"),
    S1 = c(1, 3),
    S2 = c(7, 15),
    stringsAsFactors = FALSE
  )
  lipid_object <- newLipidLogTransformObject(assay_data = assay_data)
  lipid_object@lipid_data <- unname(lipid_object@lipid_data)

  transformed <- suppressWarnings(logTransformAssays(lipid_object, offset = 1))
  transformed_assay <- transformed@lipid_data$Assay_1

  expect_s3_class(transformed_assay, "tbl_df")
  expect_equal(transformed_assay$S1, c(1, 2))
  expect_equal(transformed_assay$S2, c(3, 4))
})

test_that("lipidomics S4 log-transform method skips assays without required sample metadata", {
  lipid_object <- newLipidLogTransformObject()
  lipid_object@design_matrix <- data.frame(
    Run = c("S3", "S4"),
    group = c("A", "B"),
    stringsAsFactors = FALSE
  )

  transformed <- suppressWarnings(logTransformAssays(lipid_object, offset = 1))

  expect_identical(transformed@lipid_data$LCMS_Pos, lipid_object@lipid_data$LCMS_Pos)
})

test_that("lipidomics S4 log-transform method preserves current nonnumeric coercion behavior", {
  assay_data <- tibble::tibble(
    Name = c("L1", "L2"),
    Annotation = c("Alpha", "Beta"),
    S1 = list(list(1), list(3)),
    S2 = c(1, 3)
  )
  lipid_object <- newLipidLogTransformObject(assay_data = assay_data)

  transformed <- suppressWarnings(logTransformAssays(lipid_object, offset = 1))
  transformed_assay <- transformed@lipid_data$LCMS_Pos

  expect_true(all(is.na(transformed_assay$S1)))
  expect_equal(transformed_assay$S2, c(1, 2))
  expect_identical(transformed_assay$Name, assay_data$Name)
  expect_identical(transformed_assay$Annotation, assay_data$Annotation)
})
