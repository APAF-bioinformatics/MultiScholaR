# fidelity-coverage-compare: shared
library(testthat)

makeFunctionWithOverrides <- function(fun, replacements) {
  fun_override <- fun
  environment(fun_override) <- list2env(replacements, parent = environment(fun))
  fun_override
}

newPeptideSupportObject <- function(peptide_matrix, design_matrix = NULL) {
  if (is.null(design_matrix)) {
    design_matrix <- data.frame(
      Run = colnames(peptide_matrix),
      group = rep(c("A", "B"), length.out = ncol(peptide_matrix)),
      stringsAsFactors = FALSE
    )
  }

  peptide_rows <- expand.grid(
    peptide_index = seq_len(nrow(peptide_matrix)),
    Run = colnames(peptide_matrix),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  peptide_rows$Protein.Ids <- sub("%.*$", "", rownames(peptide_matrix)[peptide_rows$peptide_index])
  peptide_rows$Stripped.Sequence <- sub("^.*%", "", rownames(peptide_matrix)[peptide_rows$peptide_index])
  peptide_rows$Precursor.Normalised <- as.numeric(peptide_matrix[cbind(
    peptide_rows$peptide_index,
    match(peptide_rows$Run, colnames(peptide_matrix))
  )])
  peptide_rows$Q.Value <- 0.001
  peptide_rows$Global.Q.Value <- 0.001

  methods::new(
    "PeptideQuantitativeData",
    peptide_data = peptide_rows[
      c("Protein.Ids", "Stripped.Sequence", "Run", "Precursor.Normalised", "Q.Value", "Global.Q.Value")
    ],
    peptide_matrix = peptide_matrix,
    protein_id_column = "Protein.Ids",
    peptide_sequence_column = "Stripped.Sequence",
    q_value_column = "Q.Value",
    global_q_value_column = "Global.Q.Value",
    raw_quantity_column = "Precursor.Normalised",
    norm_quantity_column = "Precursor.Normalised",
    is_logged_data = TRUE,
    design_matrix = design_matrix,
    sample_id = "Run",
    group_id = "group",
    technical_replicate_id = NA_character_,
    args = list()
  )
}

test_that("peptide support helpers preserve NA summaries and cohort filtering branches", {
  peptide_object <- newPeptideSupportObject(
    peptide_matrix = matrix(
      c(
        1, NA, 3, NA,
        4, 5, NA, 7,
        NA, NA, 8, 9
      ),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(c("P1%pep1", "P1%pep2", "P2%pep3"), c("S1", "S2", "S3", "S4"))
    )
  )

  summary_verbose <- capture.output(
    na_results <- checkPeptideNAPercentages(peptide_object, verbose = TRUE)
  )
  expect_true(any(grepl("Missing Value Analysis", summary_verbose, fixed = TRUE)))
  expect_equal(na_results$total_na_percent, 41.66667, tolerance = 1e-5)
  expect_true(all(c("sample", "na_count", "na_percentage", "group") %in% names(na_results$per_sample_na)))
  expect_true(all(c("group", "mean_na_percentage") %in% names(na_results$per_group_na)))

  expect_error(checkPeptideNAPercentages(list(), verbose = FALSE), "PeptideQuantitativeData", fixed = TRUE)
  expect_error(
    checkPeptideNAPercentages(
      newPeptideSupportObject(
        peptide_matrix = matrix(
          c(1, 2, 3, 4),
          nrow = 2,
          dimnames = list(c("P1%pep1", "P2%pep2"), c("S1", "S2"))
        ),
        design_matrix = data.frame(Run = "S1", group = "A", stringsAsFactors = FALSE)
      ),
      verbose = FALSE
    ),
    "doesn't match",
    fixed = TRUE
  )

  input_table <- data.frame(
    Run = c("S1", "S2", "S3", "S4"),
    Protein.Ids = c("P1", "P1", "P2", "P2"),
    Stripped.Sequence = c("pep1", "pep1", "pep2", "pep2"),
    stringsAsFactors = FALSE
  )
  metadata_table <- data.frame(
    ms_filename = c("S1", "S2", "S3", "S4"),
    general_sample_info = c("HEK_control", "HEK_control", "Cohort_01", "Cohort_02"),
    stringsAsFactors = FALSE
  )

  remove_parallel <- makeFunctionWithOverrides(
    removePeptidesOnlyInHek293,
    list(
      partition = function(x, cores) x,
      collect = function(x) x
    )
  )

  direct_removed <- removePeptidesOnlyInHek293(
    input_table = input_table,
    metadata_table = metadata_table,
    input_table_sample_id_column = Run,
    sample_id_tbl_sample_id_column = ms_filename,
    protein_id_column = Protein.Ids,
    peptide_sequence_column = Stripped.Sequence,
    hek_string = "HEK",
    general_sample_info = general_sample_info,
    core_utilisation = NA
  )
  parallel_removed <- remove_parallel(
    input_table = input_table,
    metadata_table = metadata_table,
    input_table_sample_id_column = Run,
    sample_id_tbl_sample_id_column = ms_filename,
    protein_id_column = Protein.Ids,
    peptide_sequence_column = Stripped.Sequence,
    hek_string = "HEK",
    general_sample_info = general_sample_info,
    core_utilisation = 2
  )

  expect_identical(direct_removed$Protein.Ids, "P2")
  expect_identical(parallel_removed$Protein.Ids, "P2")
})

test_that("peptide support helpers preserve object comparisons and sample correlations", {
  object_a <- newPeptideSupportObject(
    peptide_matrix = matrix(
      c(
        1, 2, 3,
        4, 5, 6,
        7, 8, 9
      ),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(c("P1%pep1", "P1%pep2", "P2%pep3"), c("S1", "S2", "S3"))
    )
  )
  object_b <- newPeptideSupportObject(
    peptide_matrix = matrix(
      c(
        1, 2, 3,
        4, 5, 6,
        10, 11, 12
      ),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(c("P1%pep1", "P3%pep4", "P2%pep3"), c("S1", "S2", "S4"))
    ),
    design_matrix = data.frame(
      Run = c("S1", "S2", "S4"),
      group = c("A", "A", "B"),
      stringsAsFactors = FALSE
    )
  )

  comparison_tbl <- compareTwoPeptideDataObjects(object_a, object_b)
  summary_list <- summarisePeptideObject(object_a)
  correlation_tbl <- calculatePeptidePearsonCorrelation(
    temp_obj = list(
      data_table = data.frame(
        peptide_id = c("pep1", "pep2", "pep3"),
        S1 = c(1, 2, 3),
        S2 = c(1, 2, 3),
        TechPool = c(9, 9, 9),
        S3 = c(3, 2, 1),
        check.names = FALSE,
        stringsAsFactors = FALSE
      ),
      id_column = "peptide_id",
      design_matrix = data.frame(
        Run = c("S1", "S2", "TechPool", "S3"),
        group = c("A", "A", "QC", "B"),
        stringsAsFactors = FALSE
      ),
      sample_id = "Run"
    ),
    tech_rep_remove_regex = "Tech",
    correlation_group = "group"
  )

  expect_true(all(c("Levels", "in_a_not_b", "intersect_a_and_b", "in_b_not_a") %in% names(comparison_tbl)))
  expect_identical(summary_list$num_peptides, 3L)
  expect_identical(summary_list$num_proteins, 2L)
  expect_identical(summary_list$num_samples, 3L)
  expect_true(all(c("sample1", "sample2", "pearson_correlation") %in% names(correlation_tbl)))
  expect_true(all(correlation_tbl$sample1 != "TechPool"))
  expect_true(all(correlation_tbl$sample2 != "TechPool"))
})
