# fidelity-coverage-compare: shared
library(testthat)

writeMaxQuantDelimited <- function(path, data) {
  utils::write.table(data, path, sep = "\t", row.names = FALSE, quote = FALSE)
}

test_that("importMaxQuantData reshapes proteinGroups LFQ columns to long protein data", {
  input_path <- tempfile(fileext = ".txt")
  writeMaxQuantDelimited(
    input_path,
    data.frame(
      Protein.IDs = c("P1", "P2", "P3"),
      Gene.names = c("G1", "G2", "G3"),
      LFQ.intensity.WT_1 = c(10, 20, 30),
      LFQ.intensity.KO_1 = c(40, 50, 60),
      Potential.contaminant = c("", "+", ""),
      Reverse = c("", "", "+"),
      check.names = FALSE
    )
  )

  imported <- suppressWarnings(importMaxQuantData(input_path))

  expect_identical(imported$data_type, "protein")
  expect_identical(imported$column_mapping$protein_col, "Protein.Ids")
  expect_identical(imported$column_mapping$run_col, "Run")
  expect_identical(imported$column_mapping$quantity_col, "Intensity")
  expect_null(imported$column_mapping$peptide_col)
  expect_equal(nrow(imported$data), 2L)
  expect_identical(unique(imported$data$Protein.Ids), "P1")
  expect_setequal(unique(imported$data$Run), c("WT_1", "KO_1"))
  expect_equal(
    imported$data$Intensity[imported$data$Protein.Ids == "P1" & imported$data$Run == "KO_1"],
    40
  )
})

test_that("importMaxQuantData supports Majority protein IDs and raw intensity fallback", {
  input_path <- tempfile(fileext = ".txt")
  writeMaxQuantDelimited(
    input_path,
    data.frame(
      Majority.protein.IDs = c("P1", "P2"),
      Intensity.WT_1 = c(11, 21),
      Intensity.KO_1 = c(12, 22),
      check.names = FALSE
    )
  )

  imported <- suppressWarnings(importMaxQuantData(input_path, use_lfq = FALSE))

  expect_equal(nrow(imported$data), 4L)
  expect_setequal(unique(imported$data$Protein.Ids), c("P1", "P2"))
  expect_setequal(unique(imported$data$Run), c("WT_1", "KO_1"))
  expect_true(all(is.numeric(imported$data$Intensity)))
})

test_that("importMaxQuantData reports malformed proteinGroups exports clearly", {
  no_intensity <- tempfile(fileext = ".txt")
  no_protein <- tempfile(fileext = ".txt")

  writeMaxQuantDelimited(
    no_intensity,
    data.frame(
      Protein.IDs = "P1",
      Gene.names = "G1",
      check.names = FALSE
    )
  )
  writeMaxQuantDelimited(
    no_protein,
    data.frame(
      Gene.names = "G1",
      LFQ.intensity.WT_1 = 10,
      check.names = FALSE
    )
  )

  expect_error(
    importMaxQuantData(no_intensity),
    "No MaxQuant intensity columns found",
    fixed = TRUE
  )
  expect_error(
    importMaxQuantData(no_protein),
    "Required MaxQuant protein identifier column not found",
    fixed = TRUE
  )
})
