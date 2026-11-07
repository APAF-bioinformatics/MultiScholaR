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

test_that("importMaxQuantData drops sparse contaminant and reverse marks that a sampled type guess would miss", {
  # KIT-021 corpus check: vroom types a column from a sample of rows across the file; a "+" column with a few
  # marks among thousands of rows was read as logical with every mark lost, so decoys and contaminants survived.
  path <- withr::local_tempfile(fileext = ".txt")
  n <- 5000L
  table <- data.frame(
    Protein.IDs = sprintf("P%05d", seq_len(n)),
    Gene.names = sprintf("G%d", seq_len(n)),
    LFQ.intensity.WT_1 = seq_len(n) + 0.5,
    LFQ.intensity.KO_1 = seq_len(n) + 1.5,
    Potential.contaminant = "",
    Reverse = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  table$Potential.contaminant[[1234L]] <- "+"
  table$Reverse[[4321L]] <- "+"
  table$Protein.IDs[[1234L]] <- "CON__P01234"
  table$Protein.IDs[[4321L]] <- "REV__P04321"
  utils::write.table(table, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  imported <- importMaxQuantData(path, use_lfq = TRUE, filter_contaminants = TRUE)
  ids <- unique(imported$data$Protein.Ids)
  expect_length(ids, n - 2L)
  expect_false(any(grepl("^CON__|^REV__", ids)))
  kept <- importMaxQuantData(path, use_lfq = TRUE, filter_contaminants = FALSE)
  expect_length(unique(kept$data$Protein.Ids), n)
})

test_that("importMaxQuantData reads genuine MaxQuant output with spaced column names as it reads the dotted fixtures", {
  # KIT-054: a real proteinGroups.txt names its columns "LFQ intensity <experiment>", "Potential contaminant",
  # "Majority protein IDs"; the package's fixtures carry the dotted form. Both must give the same long table.
  dotted <- withr::local_tempfile(fileext = ".txt")
  spaced <- withr::local_tempfile(fileext = ".txt")
  n <- 300L
  table <- data.frame(
    Protein.IDs = sprintf("P%05d", seq_len(n)),
    Majority.protein.IDs = sprintf("P%05d", seq_len(n)),
    Gene.names = sprintf("G%d", seq_len(n)),
    LFQ.intensity.EA2_S1 = seq_len(n) + 0.5,
    LFQ.intensity.FB1_S2 = seq_len(n) + 1.5,
    Potential.contaminant = "",
    Reverse = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  table$Potential.contaminant[[7L]] <- "+"
  table$Reverse[[9L]] <- "+"
  utils::write.table(table, dotted, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  real <- table
  names(real) <- c("Protein IDs", "Majority protein IDs", "Gene names", "LFQ intensity EA2_S1", "LFQ intensity FB1_S2",
    "Potential contaminant", "Reverse")
  utils::write.table(real, spaced, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  from_dotted <- importMaxQuantData(dotted, use_lfq = TRUE, filter_contaminants = TRUE)
  from_spaced <- importMaxQuantData(spaced, use_lfq = TRUE, filter_contaminants = TRUE)
  expect_identical(as.data.frame(from_spaced$data), as.data.frame(from_dotted$data))
  expect_identical(from_spaced$column_mapping, from_dotted$column_mapping)
  expect_identical(from_spaced$data_type, from_dotted$data_type)
  expect_equal(length(unique(from_spaced$data$Protein.Ids)), n - 2L)
  expect_setequal(unique(from_spaced$data$Run), c("EA2_S1", "FB1_S2"))
  expect_identical(
    normaliseMaxQuantColumnNames(c("LFQ intensity A B", "Intensity A B", "Potential contaminant", "Reverse", "Other col")),
    c("LFQ.intensity.A B", "Intensity.A B", "Potential.contaminant", "Reverse", "Other col")
  )
})

test_that("KIT-101 F-18: the MaxQuant and FragPipe readers read a 0 intensity (not quantified) as NA", {
  # MaxQuant and FragPipe write 0 where a protein was not quantified in a sample. The reader turns that 0 into NA,
  # as an empty field already reads; a 0 kept as a value became log2(0 + 1) = 0 at the design step.
  maxquant <- tempfile(fileext = ".txt")
  writeLines(c(
    paste("Protein IDs", "Gene names", "LFQ intensity WT_1", "LFQ intensity WT_2", "LFQ intensity KO_1",
      "LFQ intensity KO_2", "Potential contaminant", "Reverse", sep = "\t"),
    paste("P00001", "G1", "1000", "0", "1200", "1100", "", "", sep = "\t"),
    paste("P00002", "G2", "2000", "", "2100", "1900", "", "", sep = "\t"),
    paste("P00003", "G3", "3000", "3100", "2900", "3050", "", "", sep = "\t")), maxquant)
  fragpipe <- tempfile(fileext = ".tsv")
  writeLines(c(
    paste("Protein", "Protein ID", "Gene", "WT_1 MaxLFQ Intensity", "WT_2 MaxLFQ Intensity", "KO_1 MaxLFQ Intensity",
      "KO_2 MaxLFQ Intensity", sep = "\t"),
    paste("sp|P00001|A", "P00001", "G1", "1000", "0", "1200", "1100", sep = "\t"),
    paste("sp|P00002|B", "P00002", "G2", "2000", "", "2100", "1900", sep = "\t"),
    paste("sp|P00003|C", "P00003", "G3", "3000", "3100", "2900", "3050", sep = "\t")), fragpipe)

  readers <- list(maxquant = suppressMessages(importMaxQuantData(maxquant)),
    fragpipe = suppressMessages(importFragPipeData(fragpipe)))
  for (imported in readers) {
    data <- as.data.frame(imported$data)
    cell <- function(protein, run) data$Intensity[data$Protein.Ids == protein & data$Run == run]
    expect_true(is.na(cell("P00001", "WT_2")))
    expect_true(is.na(cell("P00002", "WT_2")))
    expect_false(any(data$Intensity == 0, na.rm = TRUE))
    expect_identical(sum(is.na(data$Intensity)), 2L)
    expect_equal(cell("P00001", "WT_1"), 1000)
    expect_equal(cell("P00003", "KO_2"), 3050)
  }
})
