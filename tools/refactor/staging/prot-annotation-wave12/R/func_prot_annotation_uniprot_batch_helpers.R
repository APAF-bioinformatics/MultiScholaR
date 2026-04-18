# ----------------------------------------------------------------------------
# subsetQuery
# ----------------------------------------------------------------------------
# Filter for a batch and run analysis on that batch of uniprot accession keys only.
subsetQuery <- function(data, subset, accessions_col_name, uniprot_handle, uniprot_columns = c("EXISTENCE", "SCORE", "REVIEWED", "GENENAME", "PROTEIN-NAMES", "LENGTH"),
                        uniprot_keytype = "UniProtKB") {


  # print(subset)
  my_keys <- data |>
    dplyr::filter(round == subset) |>
    dplyr::pull({ { accessions_col_name } })

   print(head( my_keys) )

  # print(uniprot_keytype)
  # Learning in progress

  output <- UniProt.ws::select(uniprot_handle,
                     keys = my_keys,
                     columns = uniprot_columns,
                     keytype = uniprot_keytype)

  print(head(output))
}

# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# prepareUniprotBatchInput (Internal)
# ----------------------------------------------------------------------------
# Internal helper to handle the common logic for batch query preparation
.prepareUniprotBatchInput <- function(input_tbl, id_column, delim = ";", batch_size = 100, clean_isoform = TRUE) {
  patterns <- getUniprotRegexPatterns()
  
  # Ensure we have the column name as a string for sym/!!
  col_label <- rlang::as_label(rlang::enquo(id_column))
  
  res <- input_tbl |>
    dplyr::select({{ id_column }}) |>
    tidyr::separate_longer_delim({{ id_column }}, delim = delim) |>
    dplyr::mutate({{ id_column }} := trimws({{ id_column }})) |>
    dplyr::distinct({{ id_column }}) |>
    dplyr::filter(grepl(patterns$entry, {{ id_column }})) |>
    dplyr::arrange({{ id_column }})
    
  if (clean_isoform) {
    res <- res |>
      dplyr::mutate({{ id_column }} := cleanIsoformNumber({{ id_column }})) |>
      dplyr::distinct({{ id_column }})
  }
    
  res <- res |>
    dplyr::mutate(round = ceiling(dplyr::row_number() / batch_size))
    
  return(res)
}

# ----------------------------------------------------------------------------
# batchQueryEvidenceHelper
# ----------------------------------------------------------------------------
# The UniProt.ws::select function limits the number of keys queried to 100. This gives a batch number for it to be queried in batches.
batchQueryEvidenceHelper <- function(uniprot_acc_tbl, uniprot_acc_column) {
  # [OK] REFACTORED: Use centralized batch preparation logic
  .prepareUniprotBatchInput(uniprot_acc_tbl, {{ uniprot_acc_column }}, delim = ";", batch_size = 100, clean_isoform = TRUE)
}

# ----------------------------------------------------------------------------
# batchQueryEvidence
# ----------------------------------------------------------------------------
## Run evidence collection online, giving a table of keys (uniprot_acc_tbl) and the column name (uniprot_acc_column)
#'@export
batchQueryEvidence <- function(uniprot_acc_tbl, uniprot_acc_column, uniprot_handle,
                               uniprot_columns = c("EXISTENCE", "SCORE", "REVIEWED", "GENENAME", "PROTEIN-NAMES", "LENGTH"),
                               uniprot_keytype = "UniProtKB") {

  all_uniprot_acc <- batchQueryEvidenceHelper(uniprot_acc_tbl,
                                              { { uniprot_acc_column } })

  partial_subset_query <- partial(subsetQuery,
                                  data = all_uniprot_acc,
                                  accessions_col_name = { { uniprot_acc_column } },
                                  uniprot_handle = uniprot_handle,
                                  uniprot_columns = uniprot_columns,
                                  uniprot_keytype = uniprot_keytype)

  rounds_list <- all_uniprot_acc |>
    dplyr::distinct(round) |>
    dplyr::arrange(round) |>
    dplyr::pull(round)

  all_uniprot_evidence <- purrr::map(rounds_list, \(x){ partial_subset_query(subset = x) }) |>
    dplyr::bind_rows()

  return(all_uniprot_evidence)
}

# ----------------------------------------------------------------------------
# batchQueryEvidenceHelperGeneId
# ----------------------------------------------------------------------------
# The UniProt.ws::select function limits the number of keys queried to 100. This gives a batch number for it to be queried in batches.
batchQueryEvidenceHelperGeneId <- function(input_tbl, gene_id_column, delim = ":") {
  # [OK] REFACTORED: Use centralized batch preparation logic
  # NOTE: Batch size was 25 in original, keeping it 25 here as per original implementation preference
  .prepareUniprotBatchInput(input_tbl, {{ gene_id_column }}, delim = delim, batch_size = 25, clean_isoform = FALSE)
}

# ----------------------------------------------------------------------------
# batchQueryEvidenceGeneId
# ----------------------------------------------------------------------------
#' @export
batchQueryEvidenceGeneId <- function(input_tbl, gene_id_column, uniprot_handle, uniprot_keytype = "UniProtKB",
                                     uniprot_columns = c("EXISTENCE", "SCORE", "REVIEWED", "GENENAME", "PROTEIN-NAMES", "LENGTH")) {


  all_gene_id <- batchQueryEvidenceHelperGeneId(input_tbl, {{ gene_id_column }})

  partial_subset_query <- partial(subsetQuery,
                                  data = all_gene_id,
                                  accessions_col_name = {{ gene_id_column }},
                                  uniprot_handle = uniprot_handle,
                                  uniprot_columns = uniprot_columns,
                                  uniprot_keytype = uniprot_keytype)

  rounds_list <- all_gene_id |>
    dplyr::distinct(round) |>
    dplyr::arrange(round) |>
    dplyr::pull(round)

  all_uniprot_evidence <- purrr::map(rounds_list, \(x) { partial_subset_query(subset = x) }) |>
    dplyr::bind_rows()

  return(all_uniprot_evidence)
}

