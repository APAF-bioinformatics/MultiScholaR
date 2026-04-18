# ----------------------------------------------------------------------------
# extractOrganismsFromFasta
# ----------------------------------------------------------------------------
#' @title Extract Organism Information from FASTA File
#' 
#' @description Parses FASTA file headers to extract organism name (OS) and 
#' taxonomy ID (OX) for each protein entry. Used for analyzing mixed-species
#' FASTA databases.
#' 
#' @param fasta_file_path Path to the FASTA file
#' 
#' @return A tibble with columns: uniprot_acc, organism_name, taxon_id
#' 
#' @importFrom seqinr read.fasta
#' @importFrom dplyr tibble mutate
#' @importFrom stringr str_detect str_replace_all str_extract
#' @importFrom logger log_info log_warn
#' @export
extractOrganismsFromFasta <- function(fasta_file_path) {
  log_info(paste("Extracting organism information from FASTA:", fasta_file_path))
  
  # Check file exists

  if (!file.exists(fasta_file_path)) {
    stop("FASTA file not found: ", fasta_file_path)
  }
  
  # Read FASTA file
  aa_seqinr <- tryCatch({
    seqinr::read.fasta(
      file = fasta_file_path
      , seqtype = "AA"
      , whole.header = TRUE
      , as.string = TRUE
    )
  }, error = function(e) {
    stop("Failed to read FASTA file: ", e$message)
  })
  
  headers <- names(aa_seqinr)
  total_entries <- length(headers)
  
  log_info(sprintf("Processing %d FASTA entries for organism information...", total_entries))
  
  # Helper function to extract field value from FASTA header
  # Replicates getFastaFields logic inline for efficiency
  extractField <- function(header_string, field_name) {
    pattern <- paste0(field_name, "=")
    if (!stringr::str_detect(header_string, pattern)) {
      return(NA_character_)
    }
    # Extract the field value (everything after = until next field or end)
    extracted <- stringr::str_replace_all(
      header_string
      , paste0("(.*)", field_name, "=(.*?)(\\s..=.*|$)")
      , "\\2"
    )
    return(trimws(extracted))
  }
  
  # Parse each header
  result_list <- vector("list", length(headers))
  
  for (i in seq_along(headers)) {
    header <- headers[i]
    
    # Extract accession from standard UniProt format: >sp|ACCESSION|ID or >tr|ACCESSION|ID
    parts <- strsplit(header, "\\|", fixed = FALSE)[[1]]
    
    if (length(parts) >= 2) {
      accession <- parts[2]
    } else {
      # Fallback: use first word
      accession <- strsplit(header, "\\s+")[[1]][1]
      accession <- sub("^>", "", accession)
    }
    
    # Extract organism fields
    organism_name <- extractField(header, "OS")
    taxon_id_str <- extractField(header, "OX")
    taxon_id <- suppressWarnings(as.integer(taxon_id_str))
    
    result_list[[i]] <- list(
      uniprot_acc = accession
      , organism_name = organism_name
      , taxon_id = taxon_id
    )
  }
  
  # Bind into tibble
  organism_tbl <- dplyr::bind_rows(result_list)
  
  # Count organisms found
  n_with_organism <- sum(!is.na(organism_tbl$organism_name))
  n_with_taxon <- sum(!is.na(organism_tbl$taxon_id))
  
  log_info(sprintf(
    "Extracted organism info: %d/%d have organism name, %d/%d have taxon ID"
    , n_with_organism, total_entries
    , n_with_taxon, total_entries
  ))
  
  return(organism_tbl)
}

# ----------------------------------------------------------------------------
# analyzeOrganismDistribution
# ----------------------------------------------------------------------------
#' @title Analyze Organism Distribution in Proteomics Data
#' 
#' @description Analyzes the distribution of proteins across different organisms
#' by matching protein IDs from search results against FASTA organism annotations.
#' Used to identify the primary organism in mixed-species databases.
#' 
#' @param protein_ids Character vector of protein IDs from search results
#' @param organism_mapping Tibble from extractOrganismsFromFasta() with 
#'   columns: uniprot_acc, organism_name, taxon_id
#' 
#' @return A tibble with columns:
#'   \item{organism_name}{Name of the organism}
#'   \item{taxon_id}{NCBI taxonomy ID}
#'   \item{protein_count}{Number of proteins matched to this organism}
#'   \item{percentage}{Percentage of total proteins}
#'   \item{rank}{Rank by protein count (1 = highest)}
#' 
#' @importFrom dplyr tibble filter group_by summarise arrange mutate n desc row_number left_join
#' @importFrom stringr str_split
#' @importFrom logger log_info
#' @export
analyzeOrganismDistribution <- function(protein_ids, organism_mapping) {
  log_info(sprintf("Analyzing organism distribution for %d protein IDs", length(protein_ids)))
  
  # Handle protein groups (semicolon-separated IDs)
  # Extract all unique accessions from protein groups
  all_accessions <- unique(unlist(stringr::str_split(protein_ids, ";")))
  all_accessions <- trimws(all_accessions)
  all_accessions <- all_accessions[all_accessions != ""]
  
  log_info(sprintf("Expanded to %d unique accessions", length(all_accessions)))
  
  # Clean accessions (remove isoform numbers for matching)
  clean_accession <- function(acc) {
    # [OK] REFACTORED: Use centralized UniProt normalization
    normalizeUniprotAccession(acc, remove_isoform = TRUE)
  }
  
  # Create lookup table with cleaned accessions
  organism_mapping_clean <- organism_mapping |>
    dplyr::mutate(cleaned_acc = clean_accession(uniprot_acc))
  
  # Match accessions to organisms
  accession_tbl <- dplyr::tibble(
    accession = all_accessions
    , cleaned_acc = clean_accession(all_accessions)
  )
  
  # Join with organism mapping
  matched <- accession_tbl |>
    dplyr::left_join(
      organism_mapping_clean |>
        dplyr::select(cleaned_acc, organism_name, taxon_id) |>
        dplyr::distinct()
      , by = "cleaned_acc"
    )
  
  n_matched <- sum(!is.na(matched$organism_name))
  n_unmatched <- sum(is.na(matched$organism_name))
  
  log_info(sprintf("Matched %d accessions, %d unmatched", n_matched, n_unmatched))
  
  # Calculate distribution by organism
  distribution <- matched |>
    dplyr::filter(!is.na(organism_name)) |>
    dplyr::group_by(organism_name, taxon_id) |>
    dplyr::summarise(
      protein_count = dplyr::n()
      , .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(protein_count)) |>
    dplyr::mutate(
      percentage = round(protein_count / sum(protein_count) * 100, 2)
      , rank = dplyr::row_number()
    )
  
  # Add unmatched as a separate row if any
  if (n_unmatched > 0) {
    unmatched_row <- dplyr::tibble(
      organism_name = "[Unmatched/Unknown]"
      , taxon_id = NA_integer_
      , protein_count = n_unmatched
      , percentage = round(n_unmatched / length(all_accessions) * 100, 2)
      , rank = nrow(distribution) + 1
    )
    distribution <- dplyr::bind_rows(distribution, unmatched_row)
  }
  
  log_info(sprintf("Found %d distinct organisms in data", nrow(distribution) - ifelse(n_unmatched > 0, 1, 0)))
  
  return(distribution)
}

