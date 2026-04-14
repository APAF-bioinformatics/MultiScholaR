# ----------------------------------------------------------------------------
# updateProteinIDs
# ----------------------------------------------------------------------------
#' @title Update Protein IDs based on Sequence Table
#' @description Matches Protein.Ids in the data with accessions in the sequence table 
#' and replaces them with standardized IDs (e.g., database_id or ncbi_refseq).
#' @param protein_data Data frame containing a 'Protein.Ids' column.
#' @param aa_seq_tbl_final Data frame containing accessions and standardized IDs.
#' @return Updated protein_data data frame.
#' @export
updateProteinIDs <- function(protein_data, aa_seq_tbl_final) {
  # Check if ncbi_refseq column exists in aa_seq_tbl_final
  if (!"ncbi_refseq" %in% colnames(aa_seq_tbl_final)) {
    message("No ncbi_refseq column found in aa_seq_tbl_final. Returning original data unchanged.")
    return(protein_data)
  }

  # Generic NCBI protein ID patterns - escaped special characters
  ncbi_patterns <- c(
    "WP_\\d+\\.?\\d*",                     # WP_123456789.1
    "[A-Z]{2}_\\d+\\.?\\d*",              # NP_123456.1, XP_123456.1
    "[A-Z]{3}\\d+\\.?\\d*",               # ABC12345.1
    "\\w+\\.\\d+_prot_\\w+_\\d+"          # Assembly specific patterns like NZ_LR130543.1_prot_ABC_123
  )

  pattern <- paste0("(", paste(ncbi_patterns, collapse = "|"), ")")

  protein_data <- protein_data |>
    dplyr::mutate(matching_id = stringr::str_extract(Protein.Ids, pattern))

  lookup_table <- aa_seq_tbl_final |>
    dplyr::mutate(matching_id = stringr::str_extract(accession, pattern)) |>
    dplyr::select(matching_id, database_id, ncbi_refseq)

  updated_protein_data <- protein_data |>
    dplyr::left_join(lookup_table, by = "matching_id") |>
    dplyr::mutate(Protein.Ids_new = coalesce(database_id, ncbi_refseq, Protein.Ids)) |>
    dplyr::select(-matching_id, -database_id, -ncbi_refseq)

  changes <- sum(updated_protein_data$Protein.Ids_new != updated_protein_data$Protein.Ids)
  cat("Number of Protein.Ids that would be updated:", changes, "\n")

  if (changes > 0) {
    cat("\nSample of changes:\n")
    changed <- which(updated_protein_data$Protein.Ids_new != updated_protein_data$Protein.Ids)
    sample_changes <- head(changed, 5)
    for (i in sample_changes) {
      cat("Old:", updated_protein_data$Protein.Ids[i], "-> New:", updated_protein_data$Protein.Ids_new[i], "\n")
    }
  }

  # Replace old Protein.Ids with new ones
  updated_protein_data$Protein.Ids <- updated_protein_data$Protein.Ids_new
  updated_protein_data$Protein.Ids_new <- NULL

  return(updated_protein_data)
}

# ----------------------------------------------------------------------------
# cleanMaxQuantProteins
# ----------------------------------------------------------------------------
#' Clean MaxQuant Protein Data
#'
#' This function processes and cleans protein data from MaxQuant output,
#' filtering based on peptide counts and removing contaminants.
#'
#' @param fasta_file Path to input FASTA file
#' @param raw_counts_file Path to MaxQuant proteinGroups.txt file
#' @param output_counts_file Name of cleaned counts table output file
#' @param accession_record_file Name of cleaned accession to protein group mapping file
#' @param column_pattern Pattern to match intensity columns (e.g., "Reporter intensity corrected")
#' @param group_pattern Pattern to identify experimental groups (default: "")
#' @param razor_unique_peptides_group_thresh Threshold for razor + unique peptides (default: 0)
#' @param unique_peptides_group_thresh Threshold for unique peptides (default: 1)
#' @param fasta_meta_file Name of FASTA metadata RDS file (default: "aa_seq_tbl.RDS")
#' @param output_dir Directory for results (default: "results/proteomics/clean_proteins")
#' @param tmp_dir Directory for temporary files (default: "cache")
#' @param log_file Name of log file (default: "output.log")
#' @param debug Enable debug output (default: FALSE)
#' @param silent Only print critical information (default: FALSE)
#' @param no_backup Deactivate backup of previous run (default: FALSE)
#' @return List containing cleaned data and statistics
#' @import vroom magrittr knitr rlang optparse janitor tictoc configr logging
#' @export
#'
cleanMaxQuantProteins <- function(
    fasta_file,
    raw_counts_file,
    output_counts_file = "counts_table_cleaned.tab",
    accession_record_file = "cleaned_accession_to_protein_group.tab",
    column_pattern = "Reporter intensity corrected",
    group_pattern = "",
    razor_unique_peptides_group_thresh = 0,
    unique_peptides_group_thresh = 1,
    fasta_meta_file = "aa_seq_tbl.RDS",
    output_dir = "results/proteomics/clean_proteins",
    tmp_dir = "cache",
    log_file = "output.log",
    debug = FALSE,
    silent = FALSE,
    no_backup = FALSE
) {
  tic()

  # Initialize argument list with direct parameters
  args <- list(
    fasta_file = fasta_file,
    raw_counts_file = raw_counts_file,
    output_counts_file = output_counts_file,
    accession_record_file = accession_record_file,
    column_pattern = column_pattern,
    group_pattern = group_pattern,
    razor_unique_peptides_group_thresh = razor_unique_peptides_group_thresh,
    unique_peptides_group_thresh = unique_peptides_group_thresh,
    fasta_meta_file = fasta_meta_file,
    output_dir = output_dir,
    tmp_dir = tmp_dir,
    log_file = log_file,
    debug = debug,
    silent = silent,
    no_backup = no_backup
  )

  # Create directories
  if (!dir.exists(args$output_dir)) {
    dir.create(args$output_dir, recursive = TRUE)
  }
  if (!dir.exists(args$tmp_dir)) {
    dir.create(args$tmp_dir, recursive = TRUE)
  }

  # Configure logging
  logReset()
  addHandler(writeToConsole)
  addHandler(writeToFile, file = file.path(args$output_dir, args$log_file))

  level <- ifelse(args$debug, loglevels["DEBUG"], loglevels["INFO"])
  setLevel(level = ifelse(args$silent, loglevels["ERROR"], level))

  # Log start of processing
  loginfo("Starting protein data cleaning")

  # Validate required files
  required_files <- c(args$fasta_file, args$raw_counts_file)
  missing_files <- required_files[!file.exists(required_files)]
  if (length(missing_files) > 0) {
    stop("Missing required files: ", paste(missing_files, collapse = ", "))
  }

  # Set default values for pattern suffixes
  args$pattern_suffix <- "_\\d+"
  args$extract_patt_suffix <- "_(\\d+)"
  args$remove_more_peptides <- FALSE

  # Read counts file
  loginfo("Reading the counts file")
  dat_tbl <- vroom::vroom(args$raw_counts_file)

  # Clean counts table header
  loginfo("Cleaning counts table header")
  dat_cln <- janitor::clean_names(dat_tbl)
  colnames(dat_cln) <- str_replace(colnames(dat_cln), "_i_ds", "_ids")

  # Prepare regular expressions
  pattern_suffix <- args$pattern_suffix
  if (args$group_pattern != "") {
    pattern_suffix <- paste(args$pattern_suffix, tolower(args$group_pattern), sep = "_")
  }

  extract_patt_suffix <- args$extract_patt_suffix
  if (args$group_pattern != "") {
    extract_patt_suffix <- paste0(args$extract_patt_suffix, "_(", tolower(args$group_pattern), ")")
  }

  column_pattern <- tolower(paste0(make_clean_names(args$column_pattern), pattern_suffix))
  extract_replicate_group <- tolower(paste0(make_clean_names(args$column_pattern), extract_patt_suffix))

  # Prepare peptide count columns
  razor_unique_peptides_group_col <- "razor_unique_peptides"
  unique_peptides_group_col <- "unique_peptides"

  if (args$group_pattern != "") {
    razor_unique_peptides_group_col <- paste0("razor_unique_peptides_", tolower(args$group_pattern))
    unique_peptides_group_col <- paste0("unique_peptides_", tolower(args$group_pattern))
  }

  # Process FASTA file
  fasta_meta_file <- file.path(args$tmp_dir, args$fasta_meta_file)
  loginfo("Processing FASTA file")

  if (file.exists(fasta_meta_file)) {
    aa_seq_tbl <- readRDS(fasta_meta_file)
  } else {
    aa_seq_tbl <- parseFastaFile(args$fasta_file)

    saveRDS(aa_seq_tbl, fasta_meta_file)
  }

  # Process and filter data
  evidence_tbl <- dat_cln %>%
    mutate(maxquant_row_id = id)

  # Filter and clean data
  loginfo("Identify best UniProt accession per entry, extract sample number and simplify column header")

  filtered_data <- processAndFilterData(
    evidence_tbl,
    args,
    razor_unique_peptides_group_col,
    unique_peptides_group_col,
    column_pattern,
    aa_seq_tbl,
    extract_replicate_group
  )

  # Save results
  saveResults(filtered_data, args)

  # Log completion and session info
  te <- toc(quiet = TRUE)
  loginfo("%f sec elapsed", te$toc - te$tic)
  writeLines(capture.output(sessionInfo()), file.path(args$output_dir, "sessionInfo.txt"))

  return(filtered_data)
}

