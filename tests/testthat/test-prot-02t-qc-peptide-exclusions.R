library(testthat)

# Frozen pre-optimization classifier from VEqJKP. Only the classifier is copied;
# its unchanged parsing/manifest helpers resolve from the package. This reference
# must not be updated to follow the optimized implementation.
classifyPeptideExclusionsReference <- function(
    input_table,
    protein_id_column = "Protein.Group",
    protein_ids_column = "Protein.Ids",
    contaminant_manifest = NULL,
    exclude_decoys = TRUE,
    exclude_contaminants = TRUE) {
  if (!is.data.frame(input_table)) {
    stop("input_table must be a data frame.", call. = FALSE)
  }
  if (length(protein_id_column) != 1L || !protein_id_column %in% names(input_table)) {
    stop("Active protein_id_column is missing from the DIA-NN table.", call. = FALSE)
  }
  exclude_decoys <- .asPeptideExclusionFlag(exclude_decoys, "exclude_decoys")
  exclude_contaminants <- .asPeptideExclusionFlag(
    exclude_contaminants,
    "exclude_contaminants"
  )

  manifest <- readPeptideContaminantManifest(contaminant_manifest)
  provenance_columns <- unique(c(
    protein_id_column,
    if (protein_ids_column %in% names(input_table)) protein_ids_column else character()
  ))
  combined_ids <- if (nrow(input_table)) {
    apply(
      input_table[provenance_columns],
      1L,
      function(values) paste(values[!is.na(values) & nzchar(trimws(values))], collapse = ";")
    )
  } else {
    character()
  }

  decoy_columns <- .findPeptideExclusionColumns(
    names(input_table),
    c("Decoy", "Reverse", "Is.Decoy", "IsDecoy")
  )
  contaminant_columns <- .findPeptideExclusionColumns(
    names(input_table),
    c(
      "Contaminant", "Potential.contaminant", "Potential.Contaminant",
      "Is.Contaminant", "IsContaminant"
    )
  )

  any_truthy <- function(columns) {
    if (!length(columns)) {
      return(rep(FALSE, nrow(input_table)))
    }
    Reduce(`|`, lapply(columns, function(column) {
      .peptideExclusionTruthy(input_table[[column]])
    }))
  }
  explicit_decoy <- any_truthy(decoy_columns)
  explicit_contaminant <- any_truthy(contaminant_columns)
  decoy_tag <- grepl(
    "(^|[;,[:space:]])(?:REV__|DECOY_)",
    combined_ids,
    ignore.case = TRUE,
    perl = TRUE
  )
  contaminant_tag <- grepl(
    "(^|[;,[:space:]])(?:CON__|CONTAMINANT_|cRAP-)",
    combined_ids,
    ignore.case = TRUE,
    perl = TRUE
  )

  manifest_id_text <- if (protein_ids_column %in% names(input_table)) {
    input_table[[protein_ids_column]]
  } else {
    input_table[[protein_id_column]]
  }
  tokenised <- .tokenizeProteinAccessions(manifest_id_text)
  manifest_match_count <- integer(nrow(input_table))
  manifest_token_count <- lengths(tokenised)
  if (length(manifest$accessions)) {
    manifest_match_count <- vapply(tokenised, function(tokens) {
      if (!length(tokens)) {
        return(0L)
      }
      canonical_tokens <- sub("-[0-9]+$", "", tokens)
      sum(tokens %in% manifest$accessions | canonical_tokens %in% manifest$accessions)
    }, integer(1))
  }
  manifest_any <- manifest_match_count > 0L
  manifest_all <- manifest_any & manifest_match_count == manifest_token_count
  manifest_partial <- manifest_any & !manifest_all

  is_decoy <- explicit_decoy | decoy_tag
  is_contaminant <- explicit_contaminant | contaminant_tag | manifest_all
  exclusion_reason <- dplyr::case_when(
    explicit_decoy ~ "explicit_decoy_column",
    decoy_tag ~ "decoy_identifier_tag",
    explicit_contaminant ~ "explicit_contaminant_column",
    contaminant_tag ~ "contaminant_identifier_tag",
    manifest_all ~ "all_group_accessions_in_contaminant_manifest",
    manifest_partial ~ "mixed_group_partial_contaminant_manifest_match",
    TRUE ~ NA_character_
  )
  classification_source <- dplyr::case_when(
    explicit_decoy ~ paste(decoy_columns, collapse = ","),
    decoy_tag ~ "identifier_tag",
    explicit_contaminant ~ paste(contaminant_columns, collapse = ","),
    contaminant_tag ~ "identifier_tag",
    manifest_any ~ manifest$source,
    TRUE ~ "no_reliable_signal"
  )
  exclude_from_biological_analysis <-
    (is_decoy & exclude_decoys) |
    (is_contaminant & exclude_contaminants)

  annotated <- input_table
  annotated$is_decoy <- is_decoy
  annotated$is_contaminant <- is_contaminant
  annotated$contaminant_manifest_partial_match <- manifest_partial
  annotated$exclude_from_biological_analysis <- exclude_from_biological_analysis
  annotated$exclusion_reason <- exclusion_reason
  annotated$classification_source <- classification_source

  identity_columns <- unique(c(
    ".source_row_id", "Run", "Precursor.Id", protein_id_column,
    if (protein_ids_column %in% names(annotated)) protein_ids_column else character(),
    "Stripped.Sequence", "Modified.Sequence", "is_decoy", "is_contaminant",
    "exclusion_reason", "classification_source"
  ))
  identity_columns <- identity_columns[identity_columns %in% names(annotated)]
  exclusion_ledger <- annotated[exclude_from_biological_analysis, identity_columns, drop = FALSE]

  summary <- data.frame(
    input_rows = nrow(input_table),
    decoy_rows = sum(is_decoy),
    contaminant_rows = sum(is_contaminant),
    partial_manifest_group_rows = sum(manifest_partial),
    excluded_rows = sum(exclude_from_biological_analysis),
    biological_rows = sum(!exclude_from_biological_analysis),
    manifest_accessions = length(manifest$accessions),
    classification_status = if (length(manifest$accessions) ||
        length(decoy_columns) || length(contaminant_columns) ||
        any(decoy_tag) || any(contaminant_tag)) {
      "classified_from_available_signals"
    } else {
      "not_classified_no_reliable_signal"
    },
    stringsAsFactors = FALSE
  )

  list(
    raw_data = input_table,
    classified_data = annotated,
    analysis_data = annotated[!exclude_from_biological_analysis, , drop = FALSE],
    exclusion_ledger = exclusion_ledger,
    summary = summary,
    manifest = manifest
  )
}

test_that("identifier reuse is exactly equivalent to the frozen classifier", {
    skip("quarantined Q-10 (docs/qa/workflow-execution/kit/quarantined-tests.md)")
  input <- data.frame(
    Protein.Group = rep(c("g1", "REV__g2", "g3", NA, "", "  ", "x;CON__y", "g8"), 5),
    Protein.Ids = rep(c("P1", "sp|P2|name", "P02769;P1", "P02769-2", "", NA,
      " P02769 , p02769-2 ", "tr|p4|name"), 5),
    Decoy = rep(c("false", "true", NA, "0", "yes"), 8),
    Potential.contaminant = rep(c(FALSE, TRUE), 20),
    Run = rep(c("s1", "s2"), 20), .source_row_id = seq_len(40),
    quantity = rep(c(NA, NaN, Inf, -Inf, 0), 8), stringsAsFactors = FALSE)
  rownames(input) <- paste0("row-", rev(seq_len(nrow(input))))
  attr(input, "annotation") <- list(group = ordered(c("B", "A"), levels = c("A", "B", "unused")))
  factored <- input
  factored$Protein.Ids <- factor(factored$Protein.Ids, levels = c(unique(na.omit(input$Protein.Ids)), "unused"))
  named <- input
  names(named$Protein.Ids) <- paste0("accession-row-", seq_len(nrow(input)))
  numeric_group <- input
  numeric_group$Protein.Group <- rep(c(1.1, 2, NA, 1e10), 10)
  fallback <- input
  fallback$Protein.Ids <- NULL
  manifest_file <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("accession\tversion", "P02769\tv1", "P4\tv1"), manifest_file)
  for (data in list(input, factored, named, numeric_group, fallback, tibble::as_tibble(input))) {
    before <- serialize(data, NULL)
    for (manifest in list(NULL, c("P02769", "P4"), manifest_file)) {
      for (decoys in c(FALSE, TRUE)) for (contaminants in c(FALSE, TRUE)) {
        args <- list(input_table = data, contaminant_manifest = manifest,
          exclude_decoys = decoys, exclude_contaminants = contaminants)
        expect_identical(do.call(classifyPeptideBiologicalExclusions, args),
          do.call(classifyPeptideExclusionsReference, args))
      }
    }
    expect_identical(serialize(data, NULL), before)
  }
  # The preserved implementation rejects zero-row input during tokenization.
  # This performance repair does not silently change that public error behavior.
  expect_error(classifyPeptideExclusionsReference(input[FALSE, ]), "non-character argument")
  expect_error(classifyPeptideBiologicalExclusions(input[FALSE, ]), "non-character argument")
})

test_that("accession tokenization is bounded by distinct texts within each action", {
  input <- data.frame(Protein.Group = rep("g", 1000),
    Protein.Ids = rep(c("P1", "P02769;P1", NA, ""), 250))
  counts <- integer()
  original <- .tokenizeProteinAccessions
  testthat::local_mocked_bindings(.tokenizeProteinAccessions = function(text) {
    counts <<- c(counts, length(text))
    original(text)
  }, .package = "MultiScholaR")
  classifyPeptideBiologicalExclusions(input)
  classifyPeptideBiologicalExclusions(input)
  expect_identical(counts, c(3L, 3L))
})

makePeptideExclusionFixture <- function() {
  data.frame(
    .source_row_id = 1:7,
    Run = paste0("S", 1:7),
    Precursor.Id = paste0("prec", 1:7),
    Protein.Group = paste0("GROUP", 1:7),
    Protein.Ids = c(
      "P_TARGET",
      "P_DECOY_COLUMN",
      "CON__P_TAGGED",
      "P02769",
      "P02769;P_TARGET",
      "REV__P_DECOY;P_TARGET",
      "XP02769_SUFFIX"
    ),
    Stripped.Sequence = paste0("PEP", 1:7),
    Modified.Sequence = paste0("PEP", 1:7),
    Decoy = c(0, 1, 0, 0, 0, 0, 0),
    stringsAsFactors = FALSE
  )
}

test_that("DIA-NN decoys and contaminants produce an immutable analysis split", {
  input <- makePeptideExclusionFixture()
  result <- classifyPeptideBiologicalExclusions(
    input,
    contaminant_manifest = c("P02769")
  )

  expect_identical(result$raw_data, input)
  expect_identical(
    result$exclusion_ledger$.source_row_id,
    c(2L, 3L, 4L, 6L)
  )
  expect_identical(
    result$analysis_data$.source_row_id,
    c(1L, 5L, 7L)
  )
  expect_true(result$classified_data$contaminant_manifest_partial_match[[5]])
  expect_false(result$classified_data$is_contaminant[[5]])
  expect_false(result$classified_data$is_contaminant[[7]])
  expect_identical(
    result$classified_data$exclusion_reason[[5]],
    "mixed_group_partial_contaminant_manifest_match"
  )
  expect_identical(result$summary$decoy_rows, 2L)
  expect_identical(result$summary$contaminant_rows, 2L)
  expect_identical(result$summary$excluded_rows, 4L)
})

test_that("contaminant inclusion is explicit without including decoys", {
  result <- classifyPeptideBiologicalExclusions(
    makePeptideExclusionFixture(),
    contaminant_manifest = c("P02769"),
    exclude_contaminants = FALSE
  )

  expect_identical(result$exclusion_ledger$.source_row_id, c(2L, 6L))
  expect_true(all(c(3L, 4L) %in% result$analysis_data$.source_row_id))
  expect_true(result$classified_data$is_contaminant[[3]])
  expect_true(result$classified_data$is_contaminant[[4]])
})

test_that("standard reports without exclusion signals record an explicit no-op", {
  input <- data.frame(
    Protein.Group = "P1",
    Protein.Ids = "P1",
    stringsAsFactors = FALSE
  )
  result <- classifyPeptideBiologicalExclusions(input)

  expect_identical(result$analysis_data$Protein.Group, "P1")
  expect_identical(nrow(result$exclusion_ledger), 0L)
  expect_identical(
    result$summary$classification_status,
    "not_classified_no_reliable_signal"
  )
})

test_that("user manifest files are normalised and checksummed", {
  manifest_file <- tempfile(fileext = ".tsv")
  writeLines(c("accession\tversion", "sp|P02769|ALBU_BOVIN\t2012.01.01"), manifest_file)

  manifest <- readPeptideContaminantManifest(manifest_file)

  expect_identical(manifest$accessions, "P02769")
  expect_identical(manifest$version, "2012.01.01")
  expect_true(nzchar(manifest$checksum))
})

test_that("biological exclusions occur before frozen identification evidence", {
  input <- data.frame(
    Run = c("S1", "S2", "S1", "S2"),
    Precursor.Id = paste0("prec", 1:4),
    Protein.Group = c("P_TARGET", "P_TARGET", "P_CONT", "P_CONT"),
    Protein.Ids = c("P_TARGET", "P_TARGET", "P02769", "P02769"),
    Stripped.Sequence = c("PEP_A", "PEP_B", "CONT_A", "CONT_B"),
    Modified.Sequence = c("PEP_A", "PEP_B", "CONT_A", "CONT_B"),
    Precursor.Charge = 2L,
    Precursor.Quantity = 100,
    Precursor.Normalised = 10,
    Q.Value = 0.001,
    Global.Q.Value = 0.001,
    Global.PG.Q.Value = 0.001,
    Proteotypic = 1,
    stringsAsFactors = FALSE
  )

  filtered <- srlQvalueProteotypicPeptideCleanHelper(
    input,
    input_matrix_column_ids = c(
      "Run", "Precursor.Id", "Protein.Group", "Protein.Ids",
      "Stripped.Sequence", "Modified.Sequence", "Precursor.Charge",
      "Precursor.Quantity", "Precursor.Normalised"
    ),
    protein_id_column = Protein.Group,
    contaminant_manifest = "P02769",
    return_exclusion_result = TRUE
  )

  expect_identical(unique(filtered$data$Protein.Group), "P_TARGET")
  expect_identical(unique(filtered$data$identification_peptide_count), 2L)
  expect_identical(unique(filtered$data$identification_peptidoform_count), 2L)
  expect_identical(unique(filtered$exclusion_ledger$Protein.Group), "P_CONT")
  expect_identical(filtered$exclusion_summary$excluded_rows, 2L)
})

test_that("DIA-NN S4 filtering stores the exclusion ledger and manifest provenance", {
  input <- data.frame(
    Run = c("S1", "S2"),
    Precursor.Id = c("prec1", "prec2"),
    Protein.Group = c("P_TARGET", "P_CONT"),
    Protein.Ids = c("P_TARGET", "P02769"),
    Stripped.Sequence = c("PEP_A", "CONT_A"),
    Modified.Sequence = c("PEP_A", "CONT_A"),
    Precursor.Charge = 2L,
    Precursor.Quantity = 100,
    Precursor.Normalised = 10,
    Q.Value = 0.001,
    Global.Q.Value = 0.001,
    Global.PG.Q.Value = 0.001,
    Proteotypic = 1,
    stringsAsFactors = FALSE
  )
  design <- data.frame(
    Run = c("S1", "S2"),
    group = "A",
    replicates = c("R1", "R2"),
    stringsAsFactors = FALSE
  )
  object <- PeptideQuantitativeDataDiann(input, design, args = list())

  filtered <- srlQvalueProteotypicPeptideClean(
    object,
    input_matrix_column_ids = c(
      "Run", "Precursor.Id", "Protein.Group", "Protein.Ids",
      "Stripped.Sequence", "Modified.Sequence", "Precursor.Charge",
      "Precursor.Quantity", "Precursor.Normalised"
    ),
    contaminant_manifest = "P02769"
  )
  metadata <- filtered@args$srlQvalueProteotypicPeptideClean

  expect_identical(filtered@peptide_data$Protein.Group, "P_TARGET")
  expect_identical(metadata$biological_exclusion_summary$excluded_rows, 1L)
  expect_identical(metadata$biological_exclusion_ledger$Protein.Group, "P_CONT")
  expect_identical(metadata$contaminant_manifest_provenance$accessions, "P02769")
})

test_that("an explicit contaminant column takes precedence over missing manifests", {
  input <- data.frame(
    Protein.Group = c("P1", "P2"),
    Protein.Ids = c("P1", "P2"),
    Potential.contaminant = c("", "+"),
    stringsAsFactors = FALSE
  )
  result <- classifyPeptideBiologicalExclusions(input)

  expect_identical(result$analysis_data$Protein.Group, "P1")
  expect_identical(
    result$exclusion_ledger$exclusion_reason,
    "explicit_contaminant_column"
  )
})
