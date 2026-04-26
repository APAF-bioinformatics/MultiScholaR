# fidelity-coverage-compare: shared
library(testthat)

makeFunctionWithOverrides <- function(fun, replacements) {
  fun_override <- fun
  environment(fun_override) <- list2env(replacements, parent = environment(fun))
  fun_override
}

processAndFilterData <- get(
  "processAndFilterData",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
saveResults <- get(
  "saveResults",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
updateProteinIDs <- get(
  "updateProteinIDs",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
cleanMaxQuantProteins <- get(
  "cleanMaxQuantProteins",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)

test_that("annotation cleaning helpers preserve ID replacement behavior", {
  protein_data <- data.frame(
    Protein.Ids = c("WP_123456.1 sample", "NP_777777.1", "UNCHANGED"),
    abundance = c(1, 2, 3),
    stringsAsFactors = FALSE
  )
  aa_seq_tbl_final <- data.frame(
    accession = c("WP_123456.1", "NP_777777.1"),
    database_id = c("UPI000001", NA_character_),
    ncbi_refseq = c("WP_123456.1", "NP_777777.1"),
    stringsAsFactors = FALSE
  )

  updated <- updateProteinIDs(protein_data, aa_seq_tbl_final)
  untouched <- updateProteinIDs(protein_data, aa_seq_tbl_final["accession"])

  expect_identical(updated$Protein.Ids, c("UPI000001", "NP_777777.1", "UNCHANGED"))
  expect_identical(updated$abundance, protein_data$abundance)
  expect_identical(untouched, protein_data)
})

test_that("annotation processing helpers preserve filtering, optional contaminant removal, renaming, and result writing", {
  skip_if_not_installed("vroom")

  evidence_tbl <- data.frame(
    maxquant_row_id = c(1, 2, 3),
    protein_ids = c("P1;CON__BAD", "P2;P2_ALT", "REV__P3"),
    razor_unique_peptides = c("3;3", "4;4", "5;5"),
    unique_peptides = c("3;3", "4;4", "5;5"),
    reverse = c(NA_character_, NA_character_, "+"),
    potential_contaminant = c(NA_character_, NA_character_, NA_character_),
    reporter_intensity_corrected_1 = c(100, 200, 300),
    reporter_intensity_corrected_2 = c(110, 210, 310),
    stringsAsFactors = FALSE
  )

  default_args <- list(
    remove_more_peptides = TRUE,
    razor_unique_peptides_group_thresh = 2,
    unique_peptides_group_thresh = 2,
    group_pattern = ""
  )

  choose_best <- function(input_tbl,
                          acc_detail_tab,
                          accessions_column,
                          row_id_column,
                          group_id,
                          delim = ":") {
    data.frame(
      maxquant_row_id = 2,
      uniprot_acc = "BEST_P2",
      stringsAsFactors = FALSE
    )
  }

  process_default <- makeFunctionWithOverrides(
    processAndFilterData,
    list(chooseBestProteinAccessionHelper = choose_best)
  )

  default_result <- process_default(
    evidence_tbl = evidence_tbl,
    args = default_args,
    razor_unique_peptides_group_col = "razor_unique_peptides",
    unique_peptides_group_col = "unique_peptides",
    column_pattern = "reporter_intensity_corrected_\\d+",
    aa_seq_tbl = data.frame(accession = "P2", stringsAsFactors = FALSE),
    extract_replicate_group = "reporter_intensity_corrected_(\\d+)"
  )

  expect_identical(colnames(default_result$evidence_tbl_filt), c("uniprot_acc", "1", "2"))
  expect_identical(default_result$evidence_tbl_filt$uniprot_acc, "BEST_P2")
  expect_equal(
    as.numeric(default_result$num_proteins_remaining),
    c(3, 1, 1)
  )
  expect_true("protein_ids" %in% names(default_result$accession_gene_name_tbl_record))

  grouped_tbl <- data.frame(
    maxquant_row_id = 1,
    protein_ids = "P4;P4_ALT",
    razor_unique_peptides_condition = "3;3",
    unique_peptides_condition = "3;3",
    reverse = NA_character_,
    potential_contaminant = NA_character_,
    reporter_intensity_corrected_1_condition = 400,
    reporter_intensity_corrected_2_condition = 410,
    stringsAsFactors = FALSE
  )
  grouped_args <- default_args
  grouped_args$remove_more_peptides <- FALSE
  grouped_args$group_pattern <- "condition"

  grouped_result <- process_default(
    evidence_tbl = grouped_tbl,
    args = grouped_args,
    razor_unique_peptides_group_col = "razor_unique_peptides_condition",
    unique_peptides_group_col = "unique_peptides_condition",
    column_pattern = "reporter_intensity_corrected_\\d+_condition",
    aa_seq_tbl = data.frame(accession = "P4", stringsAsFactors = FALSE),
    extract_replicate_group = "reporter_intensity_corrected_(\\d+)_(condition)"
  )

  expect_identical(
    colnames(grouped_result$evidence_tbl_filt),
    c("uniprot_acc", "1_CONDITION", "2_CONDITION")
  )

  output_dir <- tempfile("annotation-processing-results-")
  dir.create(output_dir, recursive = TRUE)
  withr::defer(unlink(output_dir, recursive = TRUE, force = TRUE))

  saveResults(
    filtered_data = list(
      evidence_tbl_filt = default_result$evidence_tbl_filt,
      accession_gene_name_tbl_record = default_result$accession_gene_name_tbl_record,
      num_proteins_remaining = default_result$num_proteins_remaining
    ),
    args = list(
      output_dir = output_dir,
      output_counts_file = "counts_table_cleaned.tab",
      accession_record_file = "accession_record.tab"
    )
  )

  expect_true(file.exists(file.path(output_dir, "counts_table_cleaned.tab")))
  expect_true(file.exists(file.path(output_dir, "accession_record.tab")))
  expect_true(file.exists(file.path(output_dir, "number_of_proteins_remaining_after_each_filtering_step.tab")))
  expect_true(file.exists(file.path(output_dir, "sample_names.tab")))
})

test_that("cleanMaxQuantProteins preserves file validation and orchestration with cached and parsed FASTA metadata", {
  skip_if_not_installed("vroom")
  skip_if_not_installed("janitor")

  expect_error(
    cleanMaxQuantProteins(
      fasta_file = tempfile("missing-fasta-"),
      raw_counts_file = tempfile("missing-counts-")
    ),
    "Missing required files",
    fixed = TRUE
  )

  root <- tempfile("clean-maxquant-")
  dir.create(root, recursive = TRUE)
  withr::defer(unlink(root, recursive = TRUE, force = TRUE))

  fasta_file <- file.path(root, "input.fasta")
  raw_counts_file <- file.path(root, "proteinGroups.txt")
  writeLines(">sp|P1|TEST\nAAAA", fasta_file)
  writeLines(
    c(
      "id\tReporter intensity corrected 1",
      "1\t100"
    ),
    raw_counts_file
  )

  captured <- new.env(parent = emptyenv())
  cached_tbl <- data.frame(
    accession = "P1",
    database_id = "DB1",
    ncbi_refseq = "WP_1",
    stringsAsFactors = FALSE
  )

  clean_helper <- makeFunctionWithOverrides(
    cleanMaxQuantProteins,
    list(
      tic = function() NULL,
      toc = function(quiet = TRUE) list(tic = 1, toc = 3),
      logReset = function() NULL,
      addHandler = function(...) NULL,
      writeToConsole = function(...) NULL,
      writeToFile = function(...) NULL,
      loglevels = c(DEBUG = 1L, INFO = 2L, ERROR = 3L),
      setLevel = function(level) {
        captured$level <- level
        invisible(NULL)
      },
      loginfo = function(...) invisible(NULL),
      parseFastaFile = function(fasta_file) {
        captured$parsed_fasta <- fasta_file
        cached_tbl
      },
      processAndFilterData = function(evidence_tbl,
                                      args,
                                      razor_unique_peptides_group_col,
                                      unique_peptides_group_col,
                                      column_pattern,
                                      aa_seq_tbl,
                                      extract_replicate_group,
                                      delim = ":") {
        captured$orchestration <- list(
          evidence_tbl = evidence_tbl,
          razor_unique_peptides_group_col = razor_unique_peptides_group_col,
          unique_peptides_group_col = unique_peptides_group_col,
          column_pattern = column_pattern,
          extract_replicate_group = extract_replicate_group,
          aa_seq_tbl = aa_seq_tbl,
          output_dir = args$output_dir
        )
        list(
          evidence_tbl_filt = data.frame(uniprot_acc = "P1", SAMPLE_1 = 100, stringsAsFactors = FALSE),
          accession_gene_name_tbl_record = data.frame(maxquant_row_id = 1, protein_ids = "P1", stringsAsFactors = FALSE),
          num_proteins_remaining = c("raw" = 1)
        )
      },
      saveResults = function(filtered_data, args) {
        captured$saved <- list(filtered_data = filtered_data, args = args)
        writeLines("saved", file.path(args$output_dir, "saved.flag"))
      },
      capture.output = function(...) "sessionInfo"
    )
  )

  parsed_result <- clean_helper(
    fasta_file = fasta_file,
    raw_counts_file = raw_counts_file,
    output_dir = file.path(root, "results"),
    tmp_dir = file.path(root, "cache")
  )

  expect_true(file.exists(file.path(root, "results", "saved.flag")))
  expect_true(file.exists(file.path(root, "results", "sessionInfo.txt")))
  expect_identical(captured$parsed_fasta, fasta_file)
  expect_identical(captured$orchestration$column_pattern, "reporter_intensity_corrected_\\d+")
  expect_identical(captured$orchestration$extract_replicate_group, "reporter_intensity_corrected_(\\d+)")
  expect_identical(captured$orchestration$razor_unique_peptides_group_col, "razor_unique_peptides")
  expect_identical(captured$orchestration$unique_peptides_group_col, "unique_peptides")
  expect_identical(captured$orchestration$aa_seq_tbl, cached_tbl)
  expect_identical(parsed_result$evidence_tbl_filt$uniprot_acc, "P1")
  expect_true(file.exists(file.path(root, "cache", "aa_seq_tbl.RDS")))

  saveRDS(
    data.frame(accession = "CACHED", database_id = "DB2", ncbi_refseq = "WP_2", stringsAsFactors = FALSE),
    file.path(root, "cache", "aa_seq_tbl.RDS")
  )
  captured$parsed_fasta <- NULL

  cached_result <- clean_helper(
    fasta_file = fasta_file,
    raw_counts_file = raw_counts_file,
    output_dir = file.path(root, "results_cached"),
    tmp_dir = file.path(root, "cache")
  )

  expect_null(captured$parsed_fasta)
  expect_identical(captured$orchestration$aa_seq_tbl$accession, "CACHED")
  expect_identical(cached_result$evidence_tbl_filt$uniprot_acc, "P1")
})
