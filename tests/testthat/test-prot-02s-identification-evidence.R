library(testthat)

makeIdentificationEvidenceFixture <- function() {
  data.frame(
    Run = c(
      "S1", "S2", "S3",
      "S1", "S2",
      "S1", "S2",
      "S1", "S2"
    ),
    Precursor.Id = paste0("precursor_", seq_len(9)),
    Protein.Group = c(
      rep("P_REPEAT", 3),
      rep("P_MULTIFORM", 2),
      rep("P_TWO_PEPTIDES", 2),
      rep("P_BAD_Q", 2)
    ),
    Protein.Ids = c(
      rep("ID_REPEAT", 3),
      rep("ID_MULTIFORM", 2),
      "ID_TWO_A", "ID_TWO_B",
      rep("ID_BAD_Q", 2)
    ),
    Stripped.Sequence = c(
      rep("SAME_SEQUENCE", 3),
      rep("MODIFIABLE", 2),
      "PEPTIDE_A", "PEPTIDE_B",
      rep("BAD_Q_SEQUENCE", 2)
    ),
    Modified.Sequence = c(
      rep("SAME_SEQUENCE", 3),
      "MODIFIABLE", "M(UniMod:35)ODIFIABLE",
      "PEPTIDE_A", "PEPTIDE_B",
      "BAD_Q_SEQUENCE", "BAD_Q(UniMod:35)_SEQUENCE"
    ),
    Precursor.Charge = 2L,
    Precursor.Quantity = seq(100, 900, by = 100),
    Precursor.Normalised = seq(10, 90, by = 10),
    Q.Value = c(rep(0.001, 8), 0.02),
    Global.Q.Value = c(rep(0.001, 8), 0.02),
    Global.PG.Q.Value = c(rep(0.001, 8), 0.02),
    Proteotypic = 1,
    stringsAsFactors = FALSE
  )
}

identificationEvidenceColumns <- function() {
  c(
    "Run",
    "Precursor.Id",
    "Protein.Group",
    "Protein.Ids",
    "Stripped.Sequence",
    "Modified.Sequence",
    "Precursor.Charge",
    "Precursor.Quantity",
    "Precursor.Normalised"
  )
}

test_that("the default protein evidence rule is one distinct peptide and two distinct peptidoforms", {
  q_valid <- srlQvalueProteotypicPeptideCleanHelper(
    input_table = makeIdentificationEvidenceFixture(),
    input_matrix_column_ids = identificationEvidenceColumns(),
    protein_id_column = Protein.Group,
    peptide_sequence_column = Stripped.Sequence,
    modified_peptide_sequence_column = Modified.Sequence,
    qvalue_threshold = 0.01,
    global_qvalue_threshold = 0.01,
    choose_only_proteotypic_peptide = 1
  )

  evidence <- q_valid |>
    dplyr::distinct(
      Protein.Group,
      identification_peptide_count,
      identification_peptidoform_count
    ) |>
    dplyr::arrange(Protein.Group)

  expect_equal(
    evidence,
    data.frame(
      Protein.Group = c(
        "P_BAD_Q",
        "P_MULTIFORM",
        "P_REPEAT",
        "P_TWO_PEPTIDES"
      ),
      identification_peptide_count = c(1L, 1L, 1L, 2L),
      identification_peptidoform_count = c(1L, 2L, 1L, 2L),
      stringsAsFactors = FALSE
    )
  )

  filtered <- filterMinNumPeptidesPerProteinHelper(
    input_table = q_valid,
    protein_id_column = Protein.Group,
    peptide_sequence_column = Stripped.Sequence,
    modified_peptide_sequence_column = Modified.Sequence,
    core_utilisation = NA
  )

  expect_setequal(
    unique(filtered$Protein.Group),
    c("P_MULTIFORM", "P_TWO_PEPTIDES")
  )
  expect_true(all(filtered$peptides_for_protein_count >= 1L))
  expect_true(all(filtered$peptidoforms_for_protein_count >= 2L))
})

test_that("identification evidence is invariant to repeated runs and later row loss", {
  input <- makeIdentificationEvidenceFixture()
  duplicated_input <- dplyr::bind_rows(input, input)

  clean <- function(data) {
    srlQvalueProteotypicPeptideCleanHelper(
      input_table = data,
      input_matrix_column_ids = identificationEvidenceColumns(),
      protein_id_column = Protein.Group,
      qvalue_threshold = 0.01,
      global_qvalue_threshold = 0.01,
      choose_only_proteotypic_peptide = 1
    )
  }

  evidence <- clean(input) |>
    dplyr::distinct(
      Protein.Group,
      identification_peptide_count,
      identification_peptidoform_count
    )
  duplicated_evidence <- clean(duplicated_input) |>
    dplyr::distinct(
      Protein.Group,
      identification_peptide_count,
      identification_peptidoform_count
    )

  expect_equal(
    dplyr::arrange(evidence, Protein.Group),
    dplyr::arrange(duplicated_evidence, Protein.Group)
  )

  one_surviving_row <- clean(input) |>
    dplyr::filter(Protein.Group %in% c("P_MULTIFORM", "P_TWO_PEPTIDES")) |>
    dplyr::group_by(Protein.Group) |>
    dplyr::slice_head(n = 1L) |>
    dplyr::ungroup()

  filtered <- filterMinNumPeptidesPerProteinHelper(
    input_table = one_surviving_row,
    protein_id_column = Protein.Group,
    core_utilisation = NA
  )

  expect_setequal(
    unique(filtered$Protein.Group),
    c("P_MULTIFORM", "P_TWO_PEPTIDES")
  )
  expect_identical(
    filtered$peptides_for_protein_count[filtered$Protein.Group == "P_TWO_PEPTIDES"],
    2L
  )
})

test_that("exact evidence honours configured identity columns and missing forms", {
  input <- data.frame(
    protein = rep("P_CUSTOM", 3),
    peptide = c("PEPTIDE_A", "PEPTIDE_A", "PEPTIDE_B"),
    form = c("PEPTIDE_A", "PEPTIDE_A[mod]", NA_character_),
    stringsAsFactors = FALSE
  )

  filtered <- filterMinNumPeptidesPerProteinHelper(
    input_table = input,
    num_peptides_per_protein_thresh = 2,
    num_peptidoforms_per_protein_thresh = 2,
    protein_id_column = "protein",
    peptide_sequence_column = "peptide",
    modified_peptide_sequence_column = "form",
    core_utilisation = NA
  )

  expect_identical(unique(filtered$peptides_for_protein_count), 2L)
  expect_identical(unique(filtered$peptidoforms_for_protein_count), 2L)
})

test_that("DIA-NN objects use Protein.Group while retaining Protein.Ids provenance", {
  input <- makeIdentificationEvidenceFixture()
  design <- data.frame(
    Run = unique(input$Run),
    group = "A",
    replicates = seq_along(unique(input$Run)),
    stringsAsFactors = FALSE
  )

  object <- PeptideQuantitativeDataDiann(
    peptide_data = input,
    design_matrix = design,
    args = list()
  )

  expect_identical(object@protein_id_column, "Protein.Group")

  q_valid <- srlQvalueProteotypicPeptideClean(
    object,
    qvalue_threshold = 0.01,
    global_qvalue_threshold = 0.01,
    choose_only_proteotypic_peptide = 1,
    input_matrix_column_ids = identificationEvidenceColumns()
  )
  filtered <- filterMinNumPeptidesPerProtein(q_valid)
  rolled <- rollUpPrecursorToPeptide(filtered, core_utilisation = NA)

  expect_setequal(
    unique(filtered@peptide_data$Protein.Group),
    c("P_MULTIFORM", "P_TWO_PEPTIDES")
  )
  expect_true("Protein.Ids" %in% names(rolled@peptide_data))
  expect_true(all(c(
    "identification_peptide_count",
    "identification_peptidoform_count"
  ) %in% names(rolled@peptide_data)))
  expect_identical(
    unique(rolled@peptide_data$identification_peptide_count[
      rolled@peptide_data$Protein.Group == "P_TWO_PEPTIDES"
    ]),
    2L
  )
})

test_that("config-only historical protein evidence aliases execute exactly", {
  input <- makeIdentificationEvidenceFixture()
  design <- data.frame(
    Run = unique(input$Run),
    group = "A",
    replicates = seq_along(unique(input$Run)),
    stringsAsFactors = FALSE
  )
  object <- PeptideQuantitativeDataDiann(input, design, args = list())
  object <- srlQvalueProteotypicPeptideClean(
    object,
    input_matrix_column_ids = identificationEvidenceColumns()
  )
  object@args$filterMinNumPeptidesPerProtein <- list(
    peptides_per_protein_cutoff = 2,
    peptidoforms_per_protein_cutoff = 2,
    core_utilisation = NA
  )

  filtered <- filterMinNumPeptidesPerProtein(object)

  expect_identical(unique(filtered@peptide_data$Protein.Group), "P_TWO_PEPTIDES")
  expect_identical(
    filtered@args$filterMinNumPeptidesPerProtein$num_peptides_per_protein_thresh,
    2
  )
  expect_identical(
    filtered@args$filterMinNumPeptidesPerProtein$num_peptidoforms_per_protein_thresh,
    2
  )
})
