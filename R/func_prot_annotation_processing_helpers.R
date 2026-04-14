# ----------------------------------------------------------------------------
# processAndFilterData
# ----------------------------------------------------------------------------
#' Helper function to process and filter data
#' @noRd
#' @export
processAndFilterData <- function(
    evidence_tbl,
    args,
    razor_unique_peptides_group_col,
    unique_peptides_group_col,
    column_pattern,
    aa_seq_tbl,
    extract_replicate_group,
    delim = ":"
) {
  # Initialize tracking of protein numbers
  num_proteins_remaining <- numeric(3)
  names(num_proteins_remaining) <- c(
    "Number of proteins in raw unfiltered file",
    "Number of proteins after removing reverse decoy and contaminant proteins",
    paste0(
      "Number of proteins after removing proteins with no. of razor + unique peptides < ",
      args$razor_unique_peptides_group_thresh,
      " and no. of unique peptides < ",
      args$unique_peptides_group_thresh
    )
  )

  # Filter and process data
  select_columns <- evidence_tbl %>%
    dplyr::select(
      maxquant_row_id,
      protein_ids,
      !!rlang::sym(razor_unique_peptides_group_col),
      !!rlang::sym(unique_peptides_group_col),
      reverse,
      potential_contaminant,
      matches(column_pattern)
    )

  num_proteins_remaining[1] <- nrow(select_columns)

  remove_reverse_and_contaminant <- select_columns  %>%
    dplyr::filter( is.na(reverse) &
                     is.na(potential_contaminant)) %>%
    dplyr::filter( !str_detect(protein_ids, "^CON__") &
                     !str_detect(protein_ids, "^REV__") )

  remove_reverse_and_contaminant_more_hits <- remove_reverse_and_contaminant

  # Remove reverse decoy peptides and contaminant peptides even if it is not the first ranked Protein IDs (e.g. it is lower down in the list of protein IDs)
  if( args$remove_more_peptides == TRUE) {
    remove_reverse_and_contaminant_more_hits <- remove_reverse_and_contaminant  %>%
      dplyr::filter( is.na(reverse) &
                       is.na(potential_contaminant)) %>%
      dplyr::filter( !str_detect(protein_ids, "CON__") &
                       !str_detect(protein_ids, "REV__") )
  }

  # Record the number of proteins after removing reverse decoy and contaminant proteins
  # The numbers will be saved into the file 'number_of_proteins_remaining_after_each_filtering_step.tab'
  num_proteins_remaining[2] <- nrow(remove_reverse_and_contaminant_more_hits)

  helper_unnest_unique_and_razor_peptides <- remove_reverse_and_contaminant_more_hits %>%
    dplyr::mutate(protein_ids = str_split(protein_ids, ";")) %>%
    dplyr::mutate(!!rlang::sym(razor_unique_peptides_group_col) := str_split(!!rlang::sym(razor_unique_peptides_group_col), ";")) %>%
    dplyr::mutate(!!rlang::sym(unique_peptides_group_col) := str_split(!!rlang::sym(unique_peptides_group_col), ";")) %>%
    unnest(cols = c(protein_ids,
                    !!rlang::sym(razor_unique_peptides_group_col),
                    !!rlang::sym(unique_peptides_group_col)))


  evidence_tbl_cleaned <- helper_unnest_unique_and_razor_peptides %>%
    dplyr::filter(!!rlang::sym(razor_unique_peptides_group_col) >= args$razor_unique_peptides_group_thresh &
                    !!rlang::sym(unique_peptides_group_col) >= args$unique_peptides_group_thresh)


  ## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  accession_gene_name_tbl <- chooseBestProteinAccessionHelper(input_tbl = evidence_tbl_cleaned,
                                                        acc_detail_tab = aa_seq_tbl,
                                                        accessions_column = protein_ids,
                                                        row_id_column = "uniprot_acc",
                                                        group_id = maxquant_row_id,
                                                        delim = delim)



   print( accession_gene_name_tbl|>
    dplyr::filter (str_detect( uniprot_acc, "A0A024R1R8")) )



  accession_gene_name_tbl_record <- accession_gene_name_tbl %>%
    left_join(evidence_tbl %>% dplyr::select(maxquant_row_id, protein_ids), by = c("maxquant_row_id"))


  evidence_tbl_filt <- evidence_tbl_cleaned |>
    inner_join(accession_gene_name_tbl |>
                 dplyr::select(maxquant_row_id, uniprot_acc), by = "maxquant_row_id") |>
    dplyr::select(uniprot_acc, matches(column_pattern), -contains(c("razor", "unique"))) |>
    distinct()

  # Record the number of proteins after removing proteins with low no. of razor + unique peptides and low no. of unique peptides
  num_proteins_remaining[3] <- nrow( evidence_tbl_filt)

  # Record the number of proteins remaining after each filtering step into the file 'number_of_proteins_remaining_after_each_filtering_step.tab'
  num_proteins_remaining_tbl <- data.frame( step=names( num_proteins_remaining), num_proteins_remaining=num_proteins_remaining)

  #TODO: This part need improvement. There is potential for bugs.
  extraction_pattern <- "\\1"
  if (args$group_pattern != "") {
    extraction_pattern <- "\\1_\\2"
  }

  colnames(evidence_tbl_filt) <- str_replace_all(colnames(evidence_tbl_filt), tolower(extract_replicate_group), extraction_pattern) %>%
    toupper( ) %>%
    str_replace_all( "UNIPROT_ACC", "uniprot_acc")


  return(list(
    evidence_tbl_filt = evidence_tbl_filt,
    num_proteins_remaining = num_proteins_remaining,
    accession_gene_name_tbl_record = accession_gene_name_tbl_record
  ))
}

# ----------------------------------------------------------------------------
# saveResults
# ----------------------------------------------------------------------------
#' Helper function to save results
#' @noRd
saveResults <- function(filtered_data, args) {
  # Save cleaned counts
  vroom::vroom_write(
    filtered_data$evidence_tbl_filt,
    file.path(args$output_dir, args$output_counts_file)
  )

  # Save accession records
  vroom::vroom_write(
    filtered_data$accession_gene_name_tbl_record,
    file.path(args$output_dir, args$accession_record_file)
  )

  # Save protein numbers
  vroom::vroom_write(
    data.frame(
      step = names(filtered_data$num_proteins_remaining),
      num_proteins_remaining = filtered_data$num_proteins_remaining
    ),
    file.path(args$output_dir, "number_of_proteins_remaining_after_each_filtering_step.tab")
  )

  # Save sample names
  sample_names <- colnames(filtered_data$evidence_tbl_filt)[-1]
  vroom::vroom_write(
    data.frame(sample_names = t(t(sample_names))),
    file.path(args$output_dir, "sample_names.tab")
  )
}

