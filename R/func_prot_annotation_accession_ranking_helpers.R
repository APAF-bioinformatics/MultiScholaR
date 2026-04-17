# ----------------------------------------------------------------------------
# chooseBestProteinAccessionHelper
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@title Choose the best UniProt accession for a protein group
#'@description From a list of UniProt accessions, choose the best accession to use based on the UniProt score for quality of annotation for the protein entries
#'@param input_tbl Contain the following columns, 'group_id' which is the Id for each protein group, 'accessions_column' which is the column with the accession of the protein
#'@param acc_detail_tabl The out table from running the function 'parseFastaFile'
#'@param accessions_column The name of the column with the list of protein accessions, separated by ';' semi-colon. No need to quote the name as we are using tidyverse programming quosure.
#'@param group_id The name of the column with the group ID for each protein group. No need to quote the name as we are using tidyverse programming quosure.
#' @returns A table with the following columns:
#'  maxquant_row_id: Row ID
#'  num_gene_names: Number of gene names associated with this row ID
#'  gene_names: The gene names
#'  uniprot_acc: List of uniprot accessions, but with the list ordered by the best one to less useful one to use
#'  is_unique: Is the protein group assined to a unique UniProt accession or multiple UniProt accessions
#'@export
chooseBestProteinAccessionHelper <- function(input_tbl
                                             , acc_detail_tab
                                             , accessions_column
                                             , row_id_column = "uniprot_acc"
                                             , group_id
                                             , delim= ":") {
  
  cat("\n\n+===========================================================================+\n")
  cat("|      ENTERING chooseBestProteinAccessionHelper (DEBUG66)              |\n")
  cat("+===========================================================================+\n\n")
  
  cat(">>> HELPER STEP 1: INPUTS <<<\n")
  cat(sprintf("   Input 'delim' parameter: '%s'\n", delim))
  cat(sprintf("   Input table dimensions: %d rows x %d cols\n", nrow(input_tbl), ncol(input_tbl)))
  accessions_col_name <- rlang::as_name(rlang::enquo(accessions_column))
  if (accessions_col_name %in% names(input_tbl)) {
    cat(sprintf("   First 5 accession values to split: %s\n", paste(head(input_tbl[[accessions_col_name]], 5), collapse = ", ")))
  }
  cat(sprintf("   row_id_column: %s\n", row_id_column))
  cat("\n")

  cat(">>> HELPER STEP 2: SPLITTING ACCESSIONS (str_split with delim) <<<\n")
  patterns <- getUniprotRegexPatterns()
  
  resolve_acc_temp <- input_tbl |>
    dplyr::select( { { group_id } }, { { accessions_column } }) |>
    mutate(row_id_column_with_isoform = str_split({ { accessions_column } }, delim)) |>
    unnest( row_id_column_with_isoform ) |>
    mutate( !!sym(row_id_column) := cleanIsoformNumber( row_id_column_with_isoform)) |>
    dplyr::filter( !str_detect(!!sym(row_id_column), patterns$decoy)) |>
    dplyr::filter( !str_detect(!!sym(row_id_column), patterns$contaminant))
  
  cat(sprintf("   After splitting and cleaning: %d rows\n", nrow(resolve_acc_temp)))
  if (nrow(resolve_acc_temp) > 0) {
    cat(sprintf("   First 5 split/cleaned accessions: %s\n", paste(head(resolve_acc_temp[[row_id_column]], 5), collapse = ", ")))
  }
  cat("\n")

  cat(">>> HELPER STEP 3: JOINING WITH FASTA DETAIL TABLE <<<\n")
  
  # Join with FASTA data
  resolve_acc_joined <- resolve_acc_temp |>
    left_join( acc_detail_tab ,
               by = join_by( !!sym(row_id_column) == !!sym(row_id_column) ),
               copy = TRUE,
               keep = NULL)
  
  # Detect which columns are available and add missing ones with defaults
  cat("   Detecting available FASTA metadata columns...\n")
  available_cols <- colnames(resolve_acc_joined)
  
  if (!"gene_name" %in% available_cols) {
    cat("   WARNING: 'gene_name' column not found in FASTA data. Using NA as default.\n")
    resolve_acc_joined <- resolve_acc_joined |> mutate(gene_name = NA_character_)
  }
  
  if (!"cleaned_acc" %in% available_cols) {
    cat("   WARNING: 'cleaned_acc' column not found. Using row_id_column as fallback.\n")
    resolve_acc_joined <- resolve_acc_joined |> mutate(cleaned_acc = !!sym(row_id_column))
  }
  
  if (!"protein_evidence" %in% available_cols) {
    cat("   WARNING: 'protein_evidence' column not found. Using neutral value (3) for sorting.\n")
    resolve_acc_joined <- resolve_acc_joined |> mutate(protein_evidence = 3)
  }
  
  if (!"status" %in% available_cols) {
    cat("   WARNING: 'status' column not found. Using 'unknown' as default.\n")
    resolve_acc_joined <- resolve_acc_joined |> mutate(status = "unknown")
  }
  
  if (!"is_isoform" %in% available_cols) {
    cat("   WARNING: 'is_isoform' column not found. Using 'Canonical' as default.\n")
    resolve_acc_joined <- resolve_acc_joined |> mutate(is_isoform = "Canonical")
  }
  
  if (!"isoform_num" %in% available_cols) {
    cat("   WARNING: 'isoform_num' column not found. Using 0 as default.\n")
    resolve_acc_joined <- resolve_acc_joined |> mutate(isoform_num = 0)
  }
  
  if (!"seq_length" %in% available_cols) {
    cat("   WARNING: 'seq_length' column not found. Using NA as default.\n")
    resolve_acc_joined <- resolve_acc_joined |> mutate(seq_length = NA_integer_)
  }
  
  # Now select and arrange using the complete set of columns
  resolve_acc_helper <- resolve_acc_joined |>
    dplyr::select( { { group_id } }, one_of(c(row_id_column, "gene_name", "cleaned_acc",
                                              "protein_evidence", "status", "is_isoform", "isoform_num", "seq_length", "annotation_score"))) |>
    distinct() |>
    arrange( { { group_id } }, protein_evidence, status, is_isoform, desc(seq_length), isoform_num)
  
  cat(sprintf("   After joining with FASTA data: %d rows\n", nrow(resolve_acc_helper)))
  cat("\n")

  cat(">>> HELPER STEP 4: SCORING AND RANKING ISOFORMS <<<\n")
  score_isoforms <- resolve_acc_helper |>
    mutate(gene_name = ifelse(is.na(gene_name) | gene_name == "", "NA", gene_name)) |>
    mutate(annotation_score = if_else(is.na(annotation_score), 0, annotation_score)) |>
    group_by({ { group_id } }, gene_name) |>
    arrange( { { group_id } }, desc(annotation_score), protein_evidence,
             status, is_isoform, desc(seq_length), isoform_num, cleaned_acc) |>
    mutate( ranking = row_number()) |>
    ungroup()
  
  cat(sprintf("   After scoring: %d rows\n", nrow(score_isoforms)))
  cat("\n")

  cat(">>> HELPER STEP 5: SELECTING BEST ACCESSION PER GENE & RECOMBINING <<<\n")
  cat(sprintf("   [WARNING]  CRITICAL: Using 'delim' for recombining: '%s'\n", delim))
  ## For each gene name find the uniprot_acc with the lowest rankinG
  group_gene_names_and_uniprot_accs <- score_isoforms |>
    distinct( { { group_id } }, gene_name, ranking) |>
    dplyr::filter(ranking == 1) |>
    left_join(score_isoforms |>
                dplyr::select({ { group_id } }, ranking, gene_name, !!sym(row_id_column), protein_evidence, annotation_score),
              by = join_by( {{ group_id }} == {{ group_id }}
                            , ranking == ranking
                            , gene_name == gene_name)) |>
    dplyr::select(-ranking) |>
    group_by({ { group_id } }) |>
    mutate(annotation_score = if_else(is.na(annotation_score), 0, annotation_score)) |>
    arrange( {{group_id}}, desc(annotation_score), protein_evidence) |>
    summarise(num_gene_names = n(),
              gene_names = paste(gene_name, collapse = delim),
              !!sym(row_id_column) := paste(!!sym(row_id_column), collapse = delim)) |>
    ungroup() |>
    mutate(is_unique = case_when(num_gene_names == 1 ~ "Unique",
                                 TRUE ~ "Multimapped"))
  
  cat(sprintf("   Final output: %d rows\n", nrow(group_gene_names_and_uniprot_accs)))
  if (nrow(group_gene_names_and_uniprot_accs) > 0) {
    cat(sprintf("   First 5 recombined accessions: %s\n", paste(head(group_gene_names_and_uniprot_accs[[row_id_column]], 5), collapse = ", ")))
  }
  cat("\n")
  cat("+===========================================================================+\n")
  cat("|      EXITING chooseBestProteinAccessionHelper (DEBUG66)               |\n")
  cat("+===========================================================================+\n\n")

  return(group_gene_names_and_uniprot_accs)

}

# ----------------------------------------------------------------------------
# rankProteinAccessionHelper
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@title Rank the UniProt accessions for a protein group
#'@description From a list of UniProt accessions, rank the accession to use based on the UniProt score for quality of annotation for the protein entries
#'@param input_tbl Contain the following columns, 'group_id' which is the Id for each protein group, 'accessions_column' which is the column with the accession of the protein
#'@param acc_detail_tabl The out table from running the function 'parseFastaFile'
#'@param accessions_column The name of the column with the list of protein accessions, separated by ';' semi-colon. No need to quote the name as we are using tidyverse programming quosure.
#'@param group_id The name of the column with the group ID for each protein group. No need to quote the name as we are using tidyverse programming quosure.
#' @returns A table with the following columns:
#'  maxquant_row_id: Row ID
#'  num_gene_names: Number of gene names associated with this row ID
#'  gene_names: The gene names
#'  uniprot_acc: List of uniprot accessions, but with the list ordered by the best one to less useful one to use
#'  is_unique: Is the protein group assined to a unique UniProt accession or multiple UniProt accessions
#'@export
rankProteinAccessionHelper <- function(input_tbl
                                       , acc_detail_tab
                                       , accessions_column
                                       , row_id_column = "uniprot_acc"
                                       , group_id
                                       , delim= ";") {
  
  cat("\n\n+===========================================================================+\n")
  cat("|      ENTERING rankProteinAccessionHelper (DEBUG66)                    |\n")
  cat("+===========================================================================+\n\n")
  
  cat(">>> RANK HELPER STEP 1: INPUTS <<<\n")
  cat(sprintf("   Input 'delim' parameter: '%s'\n", delim))
  cat(sprintf("   Input table dimensions: %d rows x %d cols\n", nrow(input_tbl), ncol(input_tbl)))
  accessions_col_name <- rlang::as_name(rlang::enquo(accessions_column))
  if (accessions_col_name %in% names(input_tbl)) {
    cat(sprintf("   First 5 accession lists to split: %s\n", paste(head(input_tbl[[accessions_col_name]], 5), collapse = " | ")))
  }
  cat(sprintf("   row_id_column: %s\n", row_id_column))
  cat("\n")

  cat(">>> RANK HELPER STEP 2: SPLITTING ACCESSION LISTS (str_split with delim) <<<\n")
  
  # Split and join with FASTA data
  resolve_acc_joined <- input_tbl |>
    dplyr::select( { { group_id } }, { { accessions_column } }) |>
    mutate( !!sym(row_id_column) := str_split({ { accessions_column } }, delim)) |>
    unnest( !!sym(row_id_column)) |>
    mutate( cleaned_acc = cleanIsoformNumber(row_id_column))   |>
    left_join( acc_detail_tab ,
               by = join_by( cleaned_acc == !!sym(row_id_column) ),
               copy = TRUE,
               keep = NULL)
  
  # Detect which columns are available and add missing ones with defaults
  cat("   Detecting available FASTA metadata columns...\n")
  available_cols <- colnames(resolve_acc_joined)
  
  if (!"gene_name" %in% available_cols) {
    cat("   WARNING: 'gene_name' column not found in FASTA data. Using NA as default.\n")
    resolve_acc_joined <- resolve_acc_joined |> mutate(gene_name = NA_character_)
  }
  
  if (!"protein_evidence" %in% available_cols) {
    cat("   WARNING: 'protein_evidence' column not found. Using neutral value (3) for sorting.\n")
    resolve_acc_joined <- resolve_acc_joined |> mutate(protein_evidence = 3)
  }
  
  if (!"status" %in% available_cols) {
    cat("   WARNING: 'status' column not found. Using 'unknown' as default.\n")
    resolve_acc_joined <- resolve_acc_joined |> mutate(status = "unknown")
  }
  
  if (!"is_isoform" %in% available_cols) {
    cat("   WARNING: 'is_isoform' column not found. Using 'Canonical' as default.\n")
    resolve_acc_joined <- resolve_acc_joined |> mutate(is_isoform = "Canonical")
  }
  
  if (!"isoform_num" %in% available_cols) {
    cat("   WARNING: 'isoform_num' column not found. Using 0 as default.\n")
    resolve_acc_joined <- resolve_acc_joined |> mutate(isoform_num = 0)
  }
  
  if (!"seq_length" %in% available_cols) {
    cat("   WARNING: 'seq_length' column not found. Using NA as default.\n")
    resolve_acc_joined <- resolve_acc_joined |> mutate(seq_length = NA_integer_)
  }
  
  # Now select and arrange using the complete set of columns
  resolve_acc_helper <- resolve_acc_joined |>
    dplyr::select( { { group_id } }, one_of(c(row_id_column, "gene_name", "cleaned_acc",
                                              "protein_evidence", "status", "is_isoform", "isoform_num", "seq_length", "annotation_score"))) |>
    distinct() |>
    arrange( { { group_id } }, protein_evidence, status, is_isoform, desc(seq_length), isoform_num)
  
  cat(sprintf("   After splitting and joining with FASTA: %d rows\n", nrow(resolve_acc_helper)))
  if (nrow(resolve_acc_helper) > 0) {
    cat(sprintf("   First 5 split accessions: %s\n", paste(head(resolve_acc_helper[[row_id_column]], 5), collapse = ", ")))
  }
  cat("\n")

  cat(">>> RANK HELPER STEP 3: SCORING AND RANKING <<<\n")
  score_isoforms <- resolve_acc_helper |>
    mutate(gene_name = ifelse(is.na(gene_name) | gene_name == "", "NA", gene_name)) |>
    mutate(annotation_score = if_else(is.na(annotation_score), 0, annotation_score)) |>
    group_by({ { group_id } }, gene_name) |>
    arrange( { { group_id } }, desc(annotation_score), protein_evidence,
             status, is_isoform, desc(seq_length), isoform_num, cleaned_acc) |>
    mutate( ranking = row_number()) |>
    ungroup()
  
  cat(sprintf("   After scoring: %d rows\n", nrow(score_isoforms)))
  cat("\n")

  cat(">>> RANK HELPER STEP 4: SELECTING BEST & RECOMBINING WITH DELIM <<<\n")
  cat(sprintf("   [WARNING]  CRITICAL: Using 'delim' for recombining: '%s'\n", delim))
  ## For each gene name find the uniprot_acc with the lowest rankinG
  group_gene_names_and_uniprot_accs <- score_isoforms |>
    distinct( { { group_id } }, gene_name, ranking) |>
    left_join(score_isoforms |>
                dplyr::select({ { group_id } }, ranking, gene_name, !!sym(row_id_column), annotation_score),
              by = join_by( {{ group_id }} == {{ group_id }}
                            , ranking == ranking
                            , gene_name == gene_name)) |>

    dplyr::select(-ranking) |>
    group_by({ { group_id } }) |>
    summarise(num_gene_names = n(),
              gene_names = paste(gene_name, collapse = delim),
              !!sym(row_id_column) := paste(!!sym(row_id_column), collapse = delim)) |>
    ungroup() |>
    mutate(is_unique = case_when(num_gene_names == 1 ~ "Unique",
                                 TRUE ~ "Multimapped"))
  
  cat(sprintf("   Final output: %d rows\n", nrow(group_gene_names_and_uniprot_accs)))
  if (nrow(group_gene_names_and_uniprot_accs) > 0) {
    cat(sprintf("   First 5 recombined accessions: %s\n", paste(head(group_gene_names_and_uniprot_accs[[row_id_column]], 5), collapse = ", ")))
  }
  cat("\n")
  cat("+===========================================================================+\n")
  cat("|      EXITING rankProteinAccessionHelper (DEBUG66)                     |\n")
  cat("+===========================================================================+\n\n")

  return(group_gene_names_and_uniprot_accs)

}

