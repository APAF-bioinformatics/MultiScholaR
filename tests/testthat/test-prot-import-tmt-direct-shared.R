# fidelity-coverage-compare: shared
library(testthat)

makeFunctionWithOverrides <- function(fun, replacements) {
  fun_override <- fun
  environment(fun_override) <- list2env(replacements, parent = environment(fun))
  fun_override
}

writePdDelimited <- function(path, data) {
  if (grepl("\\.tsv$", path, ignore.case = TRUE)) {
    utils::write.table(data, path, sep = "\t", row.names = FALSE, quote = FALSE)
  } else {
    utils::write.csv(data, path, row.names = FALSE, quote = FALSE)
  }
}

test_that("importProteomeDiscovererTMTData reshapes a single delimited export", {
  input_path <- tempfile(fileext = ".tsv")
  writePdDelimited(
    input_path,
    data.frame(
      Accession = c("P1", "P2"),
      Description = c("protein 1", "protein 2"),
      `Abundance: F1: 126, Sample_A` = c(10, 20),
      `Abundance: F2: 127N, Sample_B` = c(11, 21),
      check.names = FALSE
    )
  )

  imported <- suppressMessages(suppressWarnings(importProteomeDiscovererTMTData(input_path)))

  expect_identical(imported$data_type, "protein")
  expect_identical(imported$column_mapping$protein_col, "Protein.Ids")
  expect_identical(imported$column_mapping$run_col, "Run")
  expect_identical(imported$column_mapping$quantity_col, "Abundance")
  expect_null(imported$column_mapping$batch_col)
  expect_equal(nrow(imported$data), 4L)
  expect_identical(sort(unique(imported$data$Protein.Ids)), c("P1", "P2"))
  expect_identical(sort(unique(imported$data$Run)), c("126_Sample_A", "127N_Sample_B"))
  expect_equal(
    imported$data$Abundance[imported$data$Protein.Ids == "P2" & imported$data$Run == "127N_Sample_B"],
    21
  )
})

test_that("importProteomeDiscovererTMTData covers ZIP batching and multi-file combine paths", {
  file_one <- tempfile(fileext = ".tsv")
  file_two <- tempfile(fileext = ".tsv")
  zip_path <- tempfile(fileext = ".zip")

  writePdDelimited(
    file_one,
    data.frame(
      Accession = "P1",
      `Abundance: F1: 126, Sample_A` = 10,
      `Abundance: F2: 127N, Sample_B` = 11,
      check.names = FALSE
    )
  )
  writePdDelimited(
    file_two,
    data.frame(
      Accession = "P2",
      `Abundance: F3: 128C, Sample_C` = 12,
      `Abundance: F4: 129N, Sample_D` = 13,
      check.names = FALSE
    )
  )

  import_fun <- makeFunctionWithOverrides(
    importProteomeDiscovererTMTData,
    list(
      unzip = function(zipfile, exdir) {
        invisible(NULL)
      },
      list.files = function(path, pattern, full.names, recursive, ignore.case) {
        c(file_one, file_two)
      }
    )
  )

  imported <- suppressMessages(suppressWarnings(import_fun(zip_path)))

  expect_identical(imported$data_type, "protein")
  expect_identical(imported$column_mapping$batch_col, "Batch")
  expect_equal(nrow(imported$data), 4L)
  expect_identical(sort(unique(imported$data$Batch)), c("b1", "b2"))
  expect_true(all(startsWith(imported$data$Run, paste0(imported$data$Batch, "_"))))
  expect_equal(
    imported$data$Abundance[imported$data$Protein.Ids == "P2" & imported$data$Run == "b2_129N_Sample_D"],
    13
  )
})

test_that("importProteomeDiscovererTMTData reports malformed exports clearly", {
  missing_accession <- tempfile(fileext = ".tsv")
  no_abundance <- tempfile(fileext = ".tsv")
  zip_path <- tempfile(fileext = ".zip")

  writePdDelimited(
    missing_accession,
    data.frame(
      Protein = "P1",
      `Abundance: F1: 126, Sample_A` = 10,
      check.names = FALSE
    )
  )
  writePdDelimited(
    no_abundance,
    data.frame(
      Accession = "P1",
      Description = "protein 1",
      check.names = FALSE
    )
  )

  import_fun <- makeFunctionWithOverrides(
    importProteomeDiscovererTMTData,
    list(
      unzip = function(zipfile, exdir) {
        invisible(NULL)
      },
      list.files = function(path, pattern, full.names, recursive, ignore.case) {
        character()
      }
    )
  )

  expect_error(
    importProteomeDiscovererTMTData(missing_accession),
    "Required column 'Accession' not found",
    fixed = TRUE
  )
  expect_error(
    importProteomeDiscovererTMTData(no_abundance),
    "No 'Abundance: ' columns found",
    fixed = TRUE
  )
  expect_error(
    import_fun(zip_path),
    "No data files (.xlsx, .csv, .tsv) found in the ZIP archive.",
    fixed = TRUE
  )
})
