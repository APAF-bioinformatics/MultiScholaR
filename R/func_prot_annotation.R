# ============================================================================
# func_prot_annotation.R
# ============================================================================
# Purpose: Annotation and protein ID management functions
# 
# This file contains functions for protein annotation, UniProt data retrieval,
# protein ID cleaning and selection, FASTA file processing, and phosphoproteomics
# annotation. Functions in this file are used across proteomics workflows.
#
# Functions to extract here:
# - UniProt annotation functions (getUniprotAnnotations, etc.)
# - Protein ID cleaning functions (.cleanProteinIds, cleanIsoformNumber, etc.)
# - Best accession selection functions (chooseBestProteinAccession, etc.)
# - FASTA processing functions (parseFastaFile, processFastaFile, etc.)
# - Phosphoproteomics annotation functions (addPhosphositesPositionsString, etc.)
# - ID conversion functions (convertEnsemblToUniprot, etc.)
# - Additional annotation helper functions
#
# Dependencies:
# - UniProt.ws, seqinr, GO.db
# - func_general_helpers.R (for utility functions)
# ============================================================================

# TODO: Extract the following functions from their current locations:

# === UniProt Annotation Functions ===

# Function 1: getUniprotAnnotations()
# Current location: R/annotation.R
# Description: Gets UniProt annotations
# getUniprotAnnotations <- function(...) {
#   # Extract from R/annotation.R
# }

# Function 2: getUniprotAnnotationsFull()
# Current location: R/annotation.R
# Description: Gets full UniProt annotations
# getUniprotAnnotationsFull <- function(...) {
#   # Extract from R/annotation.R
# }

# Function 3: getUniProtAnnotation()
# Current location: R/annotation.R
# Description: Gets UniProt annotation for input table
# getUniProtAnnotation <- function(...) {
#   # Extract from R/annotation.R
# }

# Function 4: download_uniprot_data()
# Current location: R/annotation.R
# Description: Downloads UniProt data
# download_uniprot_data <- function(...) {
#   # Extract from R/annotation.R
# }

# Function 5: directUniprotDownload()
# Current location: R/annotation.R
# Description: Direct UniProt download
# directUniprotDownload <- function(...) {
#   # Extract from R/annotation.R
# }

# Function 6: buildAnnotationIdToAnnotationNameDictionary()
# Current location: R/enrichment_functions.R
# Description: Builds annotation ID to name dictionary
# buildAnnotationIdToAnnotationNameDictionary <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 7: buildOneProteinToAnnotationList()
# Current location: R/enrichment_functions.R
# Description: Builds one protein to annotation list
# buildOneProteinToAnnotationList <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 8: convertIdToAnnotation()
# Current location: R/enrichment_functions.R
# Description: Converts ID to annotation
# convertIdToAnnotation <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 9: convertKeyToAttribute()
# Current location: R/helper_functions.R
# Description: Converts key to attribute (also in helpers)
# convertKeyToAttribute <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 10: createIdToAttributeHash()
# Current location: R/helper_functions.R
# Description: Creates ID to attribute hash (also in helpers)
# createIdToAttributeHash <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 11: convertProteinAccToGeneSymbol()
# Current location: R/enrichment_functions.R
# Description: Converts protein accession to gene symbol
# convertProteinAccToGeneSymbol <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 12: getUniprotAccToGeneSymbolDictionary()
# Current location: R/enrichment_functions.R
# Description: Gets UniProt accession to gene symbol dictionary
# getUniprotAccToGeneSymbolDictionary <- function(...) {
#   # Extract from R/enrichment_functions.R
# }

# Function 13: convertEnsemblToUniprot()
# Current location: R/annotation.R
# Description: Converts Ensembl IDs to UniProt
# convertEnsemblToUniprot <- function(...) {
#   # Extract from R/annotation.R
# }

# Function 14: detectEnsemblIds()
# Current location: R/annotation.R
# Description: Detects Ensembl IDs
# detectEnsemblIds <- function(...) {
#   # Extract from R/annotation.R
# }

# Function 15: matchAnnotations()
# Current location: R/annotation.R
# Description: Matches annotations
# matchAnnotations <- function(...) {
#   # Extract from R/annotation.R
# }

# Function 16: convertDpcDEToStandardFormat()
# Current location: R/limpa_functions.R
# Description: Converts DPC DE to standard format
# convertDpcDEToStandardFormat <- function(...) {
#   # Extract from R/limpa_functions.R
# }

# === FASTA Processing Functions ===

# Function 17: getFastaFields()
# Current location: R/get_best_accession_helper.R
# Description: Gets FASTA fields
# getFastaFields <- function(...) {
#   # Extract from R/get_best_accession_helper.R
# }

# Function 18: parseFastaFile()
# Current location: R/get_best_accession_helper.R
# Description: Parses FASTA file
# parseFastaFile <- function(...) {
#   # Extract from R/get_best_accession_helper.R
# }

# Function 19: parseFastaObject()
# Current location: R/get_best_accession_helper.R
# Description: Parses FASTA object
# parseFastaObject <- function(...) {
#   # Extract from R/get_best_accession_helper.R
# }

# Function 20: processFastaFile()
# Current location: R/get_best_accession_helper.R
# Description: Processes FASTA file
# processFastaFile <- function(...) {
#   # Extract from R/get_best_accession_helper.R
# }

# === Protein ID Selection Functions ===

# Function 21: chooseBestProteinAccession()
# Current location: R/proteinVsSamplesS4Objects.R, R/get_best_accession_helper.R
# Type: S4 method (exportMethods) and helper
# Description: Chooses best protein accession
# setMethod(f = "chooseBestProteinAccession", ...) {
#   # Extract from R/proteinVsSamplesS4Objects.R
# }

# Function 22: chooseBestProteinAccessionHelper()
# Current location: R/get_best_accession_helper.R
# Description: Helper for choosing best protein accession
# chooseBestProteinAccessionHelper <- function(...) {
#   # Extract from R/get_best_accession_helper.R
# }

# Function 23: chooseBestProteinAccession_s3()
# Current location: R/helper_functions.R
# Description: S3 version of choose best protein accession
# chooseBestProteinAccession_s3 <- function(...) {
#   # Extract from R/helper_functions.R
# }

# Function 24: rankProteinAccessionHelper()
# Current location: R/get_best_accession_helper.R
# Description: Ranks protein accessions
# rankProteinAccessionHelper <- function(...) {
#   # Extract from R/get_best_accession_helper.R
# }

# === Protein ID Cleaning Functions ===

# Function 25: .cleanProteinIds()
# Current location: R/get_best_accession_helper.R
# Description: Cleans protein IDs (internal function)
# .cleanProteinIds <- function(...) {
#   # Extract from R/get_best_accession_helper.R
# }

# Function 26: cleanIsoformNumber()
# Current location: R/get_best_accession_helper.R
# Description: Cleans isoform number from accession
# cleanIsoformNumber <- function(...) {
#   # Extract from R/get_best_accession_helper.R
# }

# Function 27: cleanMaxQuantProteins()
# Current location: R/get_best_accession_helper.R
# Description: Cleans MaxQuant protein IDs
# cleanMaxQuantProteins <- function(...) {
#   # Extract from R/get_best_accession_helper.R
# }

# Function 28: updateProteinIDs()
# Current location: R/get_best_accession_helper.R
# Description: Updates protein IDs
# updateProteinIDs <- function(...) {
#   # Extract from R/get_best_accession_helper.R
# }

# === Phosphoproteomics Annotation Functions ===

# Function 29: getBestPosition()
# Current location: R/phosphoproteomics_helper.R
# Description: Gets best position for phosphosite
# getBestPosition <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 30: getBestPositionFutureMap()
# Current location: R/phosphoproteomics_helper.R
# Description: Gets best position using future map
# getBestPositionFutureMap <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 31: getMaxProb()
# Current location: R/phosphoproteomics_helper.R
# Description: Gets maximum probability
# getMaxProb <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 32: getMaxProbFutureMap()
# Current location: R/phosphoproteomics_helper.R
# Description: Gets max probability using future map
# getMaxProbFutureMap <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 33: getPosString()
# Current location: R/phosphoproteomics_helper.R
# Description: Gets position string
# getPosString <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 34: getXMerString()
# Current location: R/phosphoproteomics_helper.R
# Description: Gets X-mer string
# getXMerString <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 35: getXMersList()
# Current location: R/phosphoproteomics_helper.R
# Description: Gets list of X-mers
# getXMersList <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 36: addXMerStrings()
# Current location: R/phosphoproteomics_helper.R
# Description: Adds X-mer strings
# addXMerStrings <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 37: addPeptideStartAndEnd()
# Current location: R/phosphoproteomics_helper.R
# Description: Adds peptide start and end positions
# addPeptideStartAndEnd <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 38: addPhosphositesPositionsString()
# Current location: R/phosphoproteomics_helper.R
# Description: Adds phosphosite position strings
# addPhosphositesPositionsString <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 39: chooseBestPhosphositeAccession()
# Current location: R/get_best_accession_helper.R
# Description: Chooses best phosphosite accession
# chooseBestPhosphositeAccession <- function(...) {
#   # Extract from R/get_best_accession_helper.R
# }

# Function 40: getUniprotAccRankFromSitesId()
# Current location: R/phosphoproteomics_helper.R
# Description: Gets UniProt accession rank from sites ID
# getUniprotAccRankFromSitesId <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 41: formatPhosphositePosition()
# Current location: R/phosphoproteomics_helper.R
# Description: Formats phosphosite position
# formatPhosphositePosition <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 42: allPhosphositesPivotLonger()
# Current location: R/phosphoproteomics_helper.R
# Description: Pivots all phosphosites to long format
# allPhosphositesPivotLonger <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 43: allPhosphositesPivotWider()
# Current location: R/phosphoproteomics_helper.R
# Description: Pivots all phosphosites to wide format
# allPhosphositesPivotWider <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 44: uniquePhosphositesSummariseLongList()
# Current location: R/phosphoproteomics_helper.R
# Description: Summarizes unique phosphosites in long list
# uniquePhosphositesSummariseLongList <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 45: uniquePhosphositesSummariseWideList()
# Current location: R/phosphoproteomics_helper.R
# Description: Summarizes unique phosphosites in wide list
# uniquePhosphositesSummariseWideList <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 46: processMultisiteEvidence()
# Current location: R/phosphoproteomics_helper.R
# Description: Processes multisite evidence
# processMultisiteEvidence <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 47: addColumnsToEvidenceTbl()
# Current location: R/phosphoproteomics_helper.R
# Description: Adds columns to evidence table
# addColumnsToEvidenceTbl <- function(...) {
#   # Extract from R/phosphoproteomics_helper.R
# }

# Function 48: batchQueryEvidence()
# Current location: R/annotation.R
# Description: Batch queries UniProt evidence
# batchQueryEvidence <- function(...) {
#   # Extract from R/annotation.R
# }

# Function 49: batchQueryEvidenceGeneId()
# Current location: R/annotation.R
# Description: Batch queries evidence with gene ID
# batchQueryEvidenceGeneId <- function(...) {
#   # Extract from R/annotation.R
# }


# ----------------------------------------------------------------------------
# cleanIsoformNumber
# ----------------------------------------------------------------------------
#'@export
cleanIsoformNumber <- function(string) {
  # "Q8K4R4-2"
  str_replace(string, "-\\d+$", "")

}

# ----------------------------------------------------------------------------
# .cleanProteinIds
# ----------------------------------------------------------------------------
#' Clean Protein IDs
#' @description Remove isoform numbers from protein IDs (wrapper for cleanIsoformNumber)
#' @param string Character vector of protein IDs
#' @return Cleaned protein IDs with isoform numbers removed
#' @importFrom purrr partial compact imap map map_int map_chr reduce
#' @importFrom stringr str_match str_remove str_length str_extract str_split str_detect str_replace
#' @importFrom dplyr filter select mutate bind_rows arrange distinct left_join anti_join rowwise ungroup pull rename rename_with coalesce row_number
#' @importFrom seqinr read.fasta
#' @importFrom vroom vroom vroom_write
#' @importFrom janitor clean_names
#' @importFrom httr GET timeout status_code content
#' @export
.cleanProteinIds <- function(string) {
  cleanIsoformNumber(string)
}


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
# batchQueryEvidenceHelper
# ----------------------------------------------------------------------------
# The UniProt.ws::select function limits the number of keys queried to 100. This gives a batch number for it to be queried in batches.
batchQueryEvidenceHelper <- function(uniprot_acc_tbl, uniprot_acc_column) {

  all_uniprot_acc <- uniprot_acc_tbl |>
    dplyr::select({ { uniprot_acc_column } }) |>
    mutate(Proteins = str_split({ { uniprot_acc_column } }, ";")) |>
    unnest(Proteins) |>
    distinct() |>
    arrange(Proteins) |>
    mutate(Proteins = cleanIsoformNumber(Proteins)) |>
    dplyr::mutate(round = ceiling(row_number() / 100))  ## 100 is the maximum number of queries at one time
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
batchQueryEvidenceHelperGeneId <- function(input_tbl, gene_id_column, delim =":") {

  all_uniprot_acc <- input_tbl |>
    dplyr::select({ { gene_id_column } }) |>
    separate_longer_delim({ { gene_id_column } }, delim= delim ) |>
    dplyr::arrange({ { gene_id_column } }) |>
    dplyr::mutate(round = ceiling(dplyr::row_number() / 25))  ## 25 is the maximum number of queries at one time

  return(all_uniprot_acc)
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


# ----------------------------------------------------------------------------
# goIdToTerm
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' Convert a list of Gene Ontology IDs to their respective human readable name (e.g. GO term).
#' @param go_string A string consisting of a list of Gene Ontology ID, separated by a delimiter
#' @param goterms Output from running \code{goterms <- Term(GOTERM)} from the GO.db library.
#' @param gotypes Output from running \code{gotypes <- Ontology(GOTERM)} from the GO.db library.
#' @return A table with three columns. go_biological_process, go_celluar_compartment, and go_molecular_function. Each column is a list of gene ontology terms, separated by '; '.
#' @export
#' @examples
#' go_string <- "GO:0016021; GO:0030659; GO:0031410; GO:0035915; GO:0042742; GO:0045087; GO:0045335; GO:0050829; GO:0050830"
#' go_id_to_term(go_string)
goIdToTerm <- function(go_string, sep = "; ", goterms, gotypes) {

  if (!is.na(go_string)) {
    go_string_tbl <- data.frame(go_id = go_string) |>
      separate_rows(go_id, sep = sep) |>
      mutate(go_term = purrr::map_chr(go_id, function(x) { if (x %in% names(goterms)) { return(goterms[[x]]) }; return(NA) })) |>
      mutate(go_type = purrr::map_chr(go_id, function(x) { if (x %in% names(gotypes)) { return(gotypes[[x]]) }; return(NA) })) |>
      group_by(go_type) |>
      summarise(go_term = paste(go_term, collapse = "; ")) |>
      ungroup() |>
      mutate(go_type = case_when(go_type == "BP" ~ "go_biological_process",
                                 go_type == "CC" ~ "go_cellular_compartment",
                                 go_type == "MF" ~ "go_molecular_function")) |>
      pivot_wider(names_from = go_type,
                  values_from = go_term)

    return(go_string_tbl)
  }
  return(NA)
}


# ----------------------------------------------------------------------------
# uniprotGoIdToTerm
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title Convert UniProt Accession to GO Term
#' @param uniprot_dat  a table with uniprot accessions and a column with GO-ID
#' @param uniprot_id_column The name of the column with the uniprot accession, as a tidyverse header format, not a string
#' @param go_id_column The name of the column with the GO-ID, as a tidyverse header format, not a string
#' @param goterms Output from running \code{goterms <- Term(GOTERM)} from the GO.db library.
#' @param gotypes Output from running \code{gotypes <- Ontology(GOTERM)} from the GO.db library.
#' @return A table with three columns. go_biological_process, go_celluar_compartment, and go_molecular_function. Each column is a list of gene ontology terms, separated by '; '.
#' @export
uniprotGoIdToTerm <- function(uniprot_dat, uniprot_id_column = UNIPROTKB
                              , go_id_column = `GO-IDs`,  sep = "; "
                              , goterms = AnnotationDbi::Term(GO.db::GOTERM)
                              , gotypes = AnnotationDbi::Ontology(GO.db::GOTERM)) {

  uniprot_acc_to_go_id <- uniprot_dat |>
    dplyr::distinct({{uniprot_id_column}}, {{go_id_column}}) |>
    separate_rows({{go_id_column}}, sep = sep) |>
    dplyr::distinct({{uniprot_id_column}}, {{go_id_column}}) |>
    dplyr::filter(!is.na({{go_id_column}}))

  go_term_temp <- uniprot_acc_to_go_id |>
    dplyr::distinct({{go_id_column}}) |>
    mutate(go_term = purrr::map_chr({{go_id_column}}, function(x) { if (x %in% names(goterms)) { return(goterms[[x]]) }; return(NA) })) |>
    mutate(go_type = purrr::map_chr({{go_id_column}}, function(x) { if (x %in% names(gotypes)) { return(gotypes[[x]]) }; return(NA) })) |>
    mutate(go_type = case_when(go_type == "BP" ~ "go_biological_process",
                               go_type == "CC" ~ "go_cellular_compartment",
                               go_type == "MF" ~ "go_molecular_function"))

  uniprot_acc_to_go_term <- uniprot_acc_to_go_id |>
    left_join(go_term_temp, by = join_by({{go_id_column}} == {{go_id_column}})) |>
    dplyr::filter(!is.na(go_term)) |>
    group_by({{uniprot_id_column}}, go_type) |>
    summarise(go_id = paste({{go_id_column}}, collapse = "; ")
              , go_term = paste(go_term, collapse = "; ")
              , .groups = 'drop' ) |>
    pivot_wider(id_cols = {{uniprot_id_column}},
                names_from = go_type,
                values_from = c(go_term, go_id))


  output_uniprot_dat <- uniprot_dat |>
    left_join(uniprot_acc_to_go_term, by = join_by({{uniprot_id_column}} == {{uniprot_id_column}})) |>
    dplyr::select( -{{go_id_column}})

  return(output_uniprot_dat)

}


# ----------------------------------------------------------------------------
# getUniprotAnnotations
# ----------------------------------------------------------------------------
#' Download and Process UniProt Annotations
#' 
#' @description
#' Downloads protein information from UniProt for a list of protein IDs,
#' processes the results including Gene Ontology annotations, and caches
#' the result for future use.
#'
#' @param input_tbl Data frame containing protein IDs in a column named 'Protein.Ids'
#' @param cache_dir Directory path for caching the results
#' @param taxon_id Taxonomic identifier for the organism (e.g., 9606 for human)
#' @param force_download Logical; if TRUE, forces new download even if cache exists
#' @param batch_size Number of protein IDs to query in each batch
#' @param timeout Timeout in seconds for the download operation
#' @param api_delay Sleep time in seconds between API calls
#'
#' @return A data frame containing UniProt annotations and GO terms
#'
#' @export
getUniprotAnnotations <- function(input_tbl, 
                                 cache_dir, 
                                 taxon_id,
                                 force_download = FALSE,
                                 batch_size = 25,
                                 timeout = 600,
                                 api_delay = 1,
                                 progress_callback = NULL) {
  
  message("=== DEBUG66: Entering getUniprotAnnotations ===")
  message(sprintf("   Cache dir: %s", cache_dir))
  message(sprintf("   Taxon ID: %d", taxon_id))
  message(sprintf("   Input table rows: %d", nrow(input_tbl)))
  message(sprintf("   Force download: %s", force_download))
  
  # Ensure cache directory exists
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }
  
  # Define cache file paths
  cache_file <- file.path(cache_dir, "uniprot_annotations.RDS")
  raw_results_file <- file.path(cache_dir, "uniprot_results.tsv")
  
  message("   DEBUG66: Checking cache...")
  message(sprintf("   DEBUG66: Cache file: %s", cache_file))
  message(sprintf("   DEBUG66: Cache exists: %s", file.exists(cache_file)))
  
  # Check if cache exists and should be used
  if (!force_download && file.exists(cache_file)) {
    message("   DEBUG66: Cache hit, loading from cache")
    message("Loading cached UniProt annotations...")
    return(readRDS(cache_file))
  }
  
  # Download annotations if needed
  message("   DEBUG66: Cache miss, calling directUniprotDownload...")
  message("Fetching UniProt annotations...")
  annotations <- directUniprotDownload(
    input_tbl = input_tbl,
    output_path = raw_results_file,
    taxon_id = taxon_id,
    batch_size = batch_size,
    timeout = timeout,
    api_delay = api_delay,
    progress_callback = progress_callback
  )
  
  # Process annotations or create empty table if download failed
  if (!is.null(annotations) && nrow(annotations) > 0) {
    message("Processing GO terms...")
    processed_annotations <- annotations |>
      uniprotGoIdToTerm(
        uniprot_id_column = Entry,
        go_id_column = Gene.Ontology.IDs,
        sep = "; "
      )
    
    # Standardize column names
    uniprot_dat_cln <- standardizeUniprotColumns(processed_annotations)
    
    # Save to cache
    saveRDS(uniprot_dat_cln, cache_file)
    message("UniProt annotations saved to cache.")
  } else {
    warning("Failed to retrieve UniProt annotations. Using empty table.")
    uniprot_dat_cln <- createEmptyUniprotTable()
    saveRDS(uniprot_dat_cln, cache_file)
  }
  
  return(uniprot_dat_cln)
}


# ----------------------------------------------------------------------------
# directUniprotDownload
# ----------------------------------------------------------------------------
#' Download Protein Data Directly from UniProt REST API
#'
#' @description
#' Downloads protein information from UniProt REST API for a list of protein IDs.
#' Processes proteins in batches to avoid overwhelming the API.
#'
#' @param input_tbl Data frame containing protein IDs in a column named 'Protein.Ids'
#' @param output_path File path to save the raw results
#' @param taxon_id Taxonomic identifier for the organism
#' @param batch_size Number of protein IDs to query in each batch
#' @param timeout Timeout in seconds for the download operation
#' @param api_delay Sleep time in seconds between API calls
#'
#' @return A data frame containing the raw UniProt results
#'
#' @export
directUniprotDownload <- function(input_tbl, 
                                 output_path, 
                                 taxon_id, 
                                 batch_size = 25,
                                 timeout = 600, 
                                 api_delay = 1,
                                 progress_callback = NULL) {
  
  message("=== DEBUG66: Entering directUniprotDownload ===")
  message(sprintf("   Input: %d unique protein IDs", nrow(input_tbl)))
  message(sprintf("   Taxon ID: %d", taxon_id))
  message(sprintf("   Batch size: %d", batch_size))
  message(sprintf("   Output path: %s", output_path))
  
  # Set a longer timeout
  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout))
  options(timeout = timeout)
  
  # Extract unique protein IDs
  protein_ids <- unique(input_tbl$Protein.Ids)
  message(paste("Found", length(protein_ids), "unique protein IDs to query"))
  message("   DEBUG66: First 10 protein IDs:")
  print(head(protein_ids, 10))
  
  # Split into batches
  chunks <- split(protein_ids, ceiling(seq_along(protein_ids)/batch_size))
  message(paste("Split into", length(chunks), "chunks for processing"))
  
  # Function to process one chunk
  process_chunk <- function(chunk, chunk_idx, total_chunks) {
    # Convert chunk_idx to integer (purrr::imap passes names as characters)
    chunk_idx <- as.integer(chunk_idx)
    
    message(sprintf("=== DEBUG66: Entering process_chunk %d/%d ===", chunk_idx, total_chunks))
    message(paste("Processing chunk", chunk_idx, "of", total_chunks, "with", length(chunk), "IDs"))
    message("   Chunk protein IDs:")
    print(chunk)
    
    # Create query for this batch
    message("   DEBUG66: Constructing query string...")
    query <- paste0("(", paste(chunk, collapse=" OR "), ") AND organism_id:", taxon_id)
    message(sprintf("   DEBUG66: Query string length: %d characters", nchar(query)))
    message(sprintf("   DEBUG66: Query string: %s", substr(query, 1, 500)))  # First 500 chars
    
    # Use httr to download
    response <- httr::GET(
      url = "https://rest.uniprot.org/uniprotkb/search",
      query = list(
        query = query,
        format = "tsv",
        fields = "accession,id,protein_name,gene_names,organism_name,length,go_id,reviewed,protein_existence,annotation_score"
      ),
      httr::timeout(30)
    )
    
    # Be nice to the API
    Sys.sleep(api_delay)
    
    # Update progress every 5 chunks
    if (!is.null(progress_callback) && chunk_idx %% 5 == 0) {
      progress_callback(chunk_idx, total_chunks)
    }
    
    # Check if successful
    if (httr::status_code(response) == 200) {
      message("   DEBUG66: Response status 200, processing content...")
      tryCatch({
        message("   DEBUG66: Extracting content from response...")
        content <- httr::content(response, "text", encoding = "UTF-8")
        message(sprintf("   DEBUG66: Content retrieved, length: %d chars", nchar(content)))
        message(sprintf("   DEBUG66: First 200 chars of content: %s", substr(content, 1, 200)))
        
        temp_file <- tempfile(fileext = ".tsv")
        message(sprintf("   DEBUG66: Writing to temp file: %s", temp_file))
        writeLines(content, temp_file)
        
        message("   DEBUG66: Reading TSV from temp file...")
        chunk_result <- suppressWarnings(
          read.delim(temp_file, sep="\t", quote="", stringsAsFactors=FALSE)
        )
        
        message(sprintf("   DEBUG66: Parsed TSV, rows: %d, cols: %d", nrow(chunk_result), ncol(chunk_result)))
        if (ncol(chunk_result) > 0) {
          message("   DEBUG66: Column names:")
          print(names(chunk_result))
        }
        
        if (nrow(chunk_result) > 0) {
          message(paste("  Found", nrow(chunk_result), "results"))
          message("   DEBUG66: Returning chunk result")
          return(chunk_result)
        } else {
          message("  API returned 200 but result has 0 rows")
          return(NULL)
        }
      }, error = function(e) {
        message("   DEBUG66: ERROR processing API response")
        message(paste("      Error message:", e$message))
        message(paste("      Error class:", paste(class(e), collapse = ", ")))
        message(paste("      Error occurred at chunk", chunk_idx))
        print(str(e))
        return(NULL)
      })
    } else {
      message(paste("  Request failed with status", httr::status_code(response)))
      message("   DEBUG66: Non-200 status, returning NULL")
      return(NULL)
    }
    
    return(NULL)
  }
  
  # Process all chunks using imap (provides both value and index)
  total_chunks <- length(chunks)
  
  message(sprintf("=== DEBUG66: About to process %d chunks ===", total_chunks))
  message("   First chunk protein IDs:")
  print(head(chunks[[1]], 10))
  
  # Wrap purrr::imap in tryCatch to catch errors with full context
  results <- tryCatch({
    message("   DEBUG66: Starting purrr::imap...")
    res <- purrr::imap(chunks, ~ process_chunk(.x, .y, total_chunks)) |>
      purrr::compact() # Remove NULL results
    message(sprintf("=== DEBUG66: purrr::imap completed. Got %d results ===", length(res)))
    res
  }, error = function(e) {
    message("=== DEBUG66: ERROR in purrr::imap ===")
    message(sprintf("   Error message: %s", e$message))
    message(sprintf("   Error class: %s", paste(class(e), collapse = ", ")))
    message("   Full error structure:")
    print(str(e))
    stop(e)
  })
  
  # Send final progress update to show 100% completion
  if (!is.null(progress_callback)) {
    progress_callback(total_chunks, total_chunks)
  }
  
  # Combine results
  if (length(results) > 0) {
    tryCatch({
      all_results <- purrr::reduce(results, rbind)
      
      # Standardize column names for downstream processing
      names(all_results) <- gsub(" ", ".", names(all_results))
      
      # Add From column needed for downstream processing
      all_results$From <- all_results$Entry
      
      # Write to file
      write.table(all_results, output_path, sep="\t", quote=FALSE, row.names=FALSE)
      message(paste("Successfully retrieved", nrow(all_results), "entries from UniProt"))
      return(all_results)
    }, error = function(e) {
      message(paste("Error combining results from chunks:", e$message))
      message(paste("Number of successful chunks:", length(results)))
      if (length(results) > 0) {
        message(paste("First chunk columns:", paste(names(results[[1]]), collapse = ", ")))
      }
      return(NULL)
    })
  } else {
    message("Failed to retrieve any data from UniProt")
    return(NULL)
  }
}


# ----------------------------------------------------------------------------
# standardizeUniprotColumns
# ----------------------------------------------------------------------------
#' Standardize UniProt Column Names
#'
#' @description
#' Standardizes column names from UniProt results for downstream processing.
#' Handles missing columns gracefully.
#'
#' @param df Data frame with UniProt results
#'
#' @return Data frame with standardized column names
#'
#' @keywords internal
standardizeUniprotColumns <- function(df) {
  # Handle Protein existence column
  if ("Protein.existence" %in% colnames(df)) {
    df <- df |> dplyr::rename(Protein_existence = "Protein.existence")
  } else {
    df$Protein_existence <- NA_character_
  }
  
  # Handle Protein names column
  if ("Protein.names" %in% colnames(df)) {
    df <- df |> dplyr::rename(Protein_names = "Protein.names")
  } else {
    df$Protein_names <- NA_character_
  }
  
  # Handle Gene Names column (different versions of API may use different names)
  gene_names_col <- grep("Gene\\.Names", colnames(df), value = TRUE)
  if (length(gene_names_col) > 0) {
    df <- df |> dplyr::rename_with(~"gene_names", .cols = matches("Gene.Names"))
  } else {
    df$gene_names <- NA_character_
  }

  # Handle Annotation Score column - ensure it's always present with fallback to 0
  if ("Annotation.score" %in% colnames(df)) {
    df <- df |> dplyr::rename(annotation_score = "Annotation.score")
  } else if (!"annotation_score" %in% colnames(df)) {
    df$annotation_score <- 0
  }

  # Ensure annotation_score is numeric and handle any NA values
  df$annotation_score <- as.numeric(df$annotation_score)
  df$annotation_score <- ifelse(is.na(df$annotation_score), 0, df$annotation_score)

  return(df)
}


# ----------------------------------------------------------------------------
# createEmptyUniprotTable
# ----------------------------------------------------------------------------
#' Create Empty UniProt Table
#'
#' @description
#' Creates an empty table with standard UniProt columns when download fails.
#'
#' @return Empty data frame with standard UniProt columns
#'
#' @keywords internal
createEmptyUniprotTable <- function() {
  data.frame(
    Entry = character(0),
    From = character(0),
    "gene_names" = character(0),
    "Protein_existence" = character(0),
    "Protein_names" = character(0),
    "annotation_score" = numeric(0),
    stringsAsFactors = FALSE
  )
}


# ----------------------------------------------------------------------------
# getUniProtAnnotation
# ----------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
getUniProtAnnotation <- function(   input_table, taxonomy_id  =9606, protein_id_column = "Protein.Ids",  protein_id_delim=":", output_dir = "." ) {


  uniprot_dat <- NA
  uniprot_acc_tbl <- input_table |>
    mutate( uniprot_acc_copy = !!sym(protein_id_column) )  |>
    separate_rows(uniprot_acc_copy, sep=protein_id_delim )  |>
    mutate( join_uniprot_acc = cleanIsoformNumber(uniprot_acc_copy))  |>
    dplyr::distinct( !!sym(protein_id_column) , join_uniprot_acc)  |>
    group_by( !!sym(protein_id_column) )  |>
    mutate( acc_order_id = row_number())  |>
    ungroup()

  uniprot_file<-file.path( output_dir, "uniprot_dat.rds")
  if( ! file.exists( uniprot_file )) {

    up <- UniProt.ws(taxId= taxonomy_id)
    list_of_sp_columns <- c("EXISTENCE"
                            , "SCORE"
                            , "REVIEWED"
                            , "GENENAME"
                            , "PROTEIN-NAMES"
                            , "LENGTH"
                            , "ENSEMBL"
                            , "GO-ID"
                            , "KEYWORDS"

                            ,"protein_existence"
                            ,"annotation_score"#?
                            ,"reviewed"
                            ,"gene_names"
                            ,"protein_name"
                            ,"length"
                            ,"xref_ensembl"
                            , "go_id"
                            , "keyword"
    )
    up_cls<-unlist(columns(up))
    list_intersect<-intersect(list_of_sp_columns,up_cls)
    if(length(setdiff( list_of_sp_columns,list_intersect)) > 0)
    {
      print(paste("UniProt fields not found:", setdiff( list_of_sp_columns,list_intersect),sep=", "))
    }

    my_keytype <- "UniProtKB"
    if( "UNIPROTKB" %in% keytypes(up) ) {
      my_keytype <- "UNIPROTKB"
    }

    uniprot_dat <- batchQueryEvidence(uniprot_acc_tbl, join_uniprot_acc, uniprot_handle=up,
                                      uniprot_columns = list_intersect, uniprot_keytype=my_keytype)

    if( my_keytype == "UniProtKB") {
      uniprot_dat <- uniprot_dat %>%
        dplyr::select(-From) %>%
        dplyr::rename( UNIPROTKB = "Entry",
                       EXISTENCE = "Protein.existence",
                       SCORE = "Annotation",
                       REVIEWED = "Reviewed",
                       GENENAME = "Gene.Names",
                       `PROTEIN-NAMES` = "Protein.names",
                       LENGTH = "Length",
                       ENSEMBL = "Ensembl",
                       `GO-ID` = "Gene.Ontology.IDs",
                       KEYWORDS   = "Keywords")


    }


    ## Merge with Gene Ontology terms.
    goterms <- Term(GOTERM)
    gotypes <- Ontology(GOTERM)


    uniprot_dat_cln <- uniprotGoIdToTerm(uniprot_dat, sep="; ", goterms, gotypes  )


    uniprot_dat_multiple_acc <- uniprot_acc_tbl %>%
      left_join( uniprot_dat_cln, by=c("join_uniprot_acc" = "UNIPROTKB") ) %>%
      arrange( !!sym(protein_id_column) , acc_order_id) %>%
      group_by(!!sym(protein_id_column)  ) %>%
      summarise( across( .cols=setdiff( colnames( uniprot_dat_cln), "UNIPROTKB")   , ~paste(., collapse=":"))   ) %>%
      ungroup() %>%
      dplyr::rename( UNIPROT_GENENAME = "GENENAME")

    saveRDS( uniprot_dat_multiple_acc, uniprot_file)

    return( uniprot_dat_multiple_acc)

  } else {
    uniprot_dat <- readRDS(uniprot_file)

    return( uniprot_dat)

  }

}


# ----------------------------------------------------------------------------
# matchAnnotations
# ----------------------------------------------------------------------------
#' Match Protein Annotations to DE Results
#'
#' @description Matches protein identifiers from differential expression results
#' with UniProt annotations, enriching the DE results with gene names and other
#' annotation information.
#'
#' @param de_results_s4 S4 object containing DE results from differential expression analysis
#' @param uniprot_annotations Data frame with UniProt annotations (typically uniprot_dat_cln)
#' @param protein_id_column Character string specifying the protein ID column name in DE results
#' @param uniprot_id_column Character string specifying the UniProt ID column in annotations (default: "Entry")
#' @param gene_names_column Character string specifying the gene names column in annotations (default: "gene_names")
#'
#' @return List containing:
#' \itemize{
#'   \item annotated_de_results: DE results enriched with UniProt annotations
#'   \item match_statistics: Summary of matching success rates
#'   \item unmatched_proteins: Proteins that could not be matched
#' }
#'
#' @details
#' This function performs fuzzy matching between protein identifiers in DE results
#' and UniProt entries. It handles common protein ID formats including:
#' \itemize{
#'   \item UniProt accessions (P12345)
#'   \item UniProt accessions with isoform suffixes (P12345-1)
#'   \item Protein group IDs separated by semicolons
#' }
#'
#' The function first attempts exact matching, then falls back to partial matching
#' by removing isoform suffixes and taking the first protein ID from groups.
#'
#' @examples
#' \dontrun{
#' # Match DE results with UniProt annotations
#' matched_results <- matchAnnotations(
#'   de_results_s4 = de_results_obj,
#'   uniprot_annotations = uniprot_dat_cln,
#'   protein_id_column = "uniprot_acc"
#' )
#' 
#' # Check matching statistics
#' print(matched_results$match_statistics)
#' }
#'
#' @export
matchAnnotations <- function(de_results_s4,
                            uniprot_annotations,
                            protein_id_column = NULL,
                            uniprot_id_column = "Entry",
                            gene_names_column = "gene_names") {
  
  log_info("=== Starting protein annotation matching ===")
  
  # Validate inputs
  if (is.null(de_results_s4)) {
    stop("de_results_s4 cannot be NULL")
  }
  
  if (is.null(uniprot_annotations) || nrow(uniprot_annotations) == 0) {
    stop("uniprot_annotations cannot be NULL or empty")
  }
  
  # Get protein ID column from S4 object if not specified
  if (is.null(protein_id_column)) {
    if (!is.null(de_results_s4@protein_id_column)) {
      protein_id_column <- de_results_s4@protein_id_column
    } else {
      stop("protein_id_column must be specified or available in S4 object @protein_id_column")
    }
  }
  
  log_info(sprintf("Using protein ID column: %s", protein_id_column))
  log_info(sprintf("UniProt annotations: %d entries", nrow(uniprot_annotations)))
  log_info(sprintf("UniProt ID column: %s", uniprot_id_column))
  log_info(sprintf("Gene names column: %s", gene_names_column))
  
  # Validate required columns exist
  if (!uniprot_id_column %in% names(uniprot_annotations)) {
    stop(sprintf("UniProt ID column '%s' not found in annotations", uniprot_id_column))
  }
  
  if (!gene_names_column %in% names(uniprot_annotations)) {
    log_warn(sprintf("Gene names column '%s' not found in annotations", gene_names_column))
    gene_names_column <- NULL
  }
  
  # Extract protein IDs from DE results
  # This could be from protein_quant_table or protein_id_table depending on S4 structure
  # OR from de_data for de_results_for_enrichment objects
  protein_ids <- NULL
  
  # ✅ ENHANCED: Handle de_results_for_enrichment objects
  if (inherits(de_results_s4, "de_results_for_enrichment")) {
    log_info("Detected de_results_for_enrichment S4 object")
    
    # Extract protein IDs from all DE data contrasts
    if (!is.null(de_results_s4@de_data) && length(de_results_s4@de_data) > 0) {
      # Get protein IDs from the first available contrast (they should all have the same proteins)
      available_contrasts <- names(de_results_s4@de_data)
      log_info(sprintf("Available DE contrasts: %s", paste(available_contrasts, collapse = ", ")))
      
      for (contrast_name in available_contrasts) {
        contrast_data <- de_results_s4@de_data[[contrast_name]]
        if (!is.null(contrast_data) && nrow(contrast_data) > 0 && protein_id_column %in% names(contrast_data)) {
          protein_ids <- unique(contrast_data[[protein_id_column]])
          log_info(sprintf("Extracted protein IDs from DE contrast: %s", contrast_name))
          break  # Use first available contrast with valid data
        }
      }
      
      if (is.null(protein_ids)) {
        stop(sprintf("Protein ID column '%s' not found in any DE data contrasts. Available columns in first contrast: %s", 
                     protein_id_column, 
                     if(length(available_contrasts) > 0 && !is.null(de_results_s4@de_data[[available_contrasts[1]]])) {
                       paste(names(de_results_s4@de_data[[available_contrasts[1]]]), collapse = ", ")
                     } else {
                       "No valid data found"
                     }))
      }
    } else {
      stop("No DE data found in de_results_for_enrichment object @de_data slot")
    }
    
  } else {
    # ✅ PRESERVED: Original logic for ProteinQuantitativeData and similar objects
    if (!is.null(de_results_s4@protein_quant_table) && protein_id_column %in% names(de_results_s4@protein_quant_table)) {
      protein_ids <- unique(de_results_s4@protein_quant_table[[protein_id_column]])
      log_info("Extracted protein IDs from @protein_quant_table")
    } else if (!is.null(de_results_s4@protein_id_table) && protein_id_column %in% names(de_results_s4@protein_id_table)) {
      protein_ids <- unique(de_results_s4@protein_id_table[[protein_id_column]])
      log_info("Extracted protein IDs from @protein_id_table")
    } else {
      stop(sprintf("Protein ID column '%s' not found in DE results S4 object. Object class: %s", 
                   protein_id_column, class(de_results_s4)[1]))
    }
  }
  
  log_info(sprintf("Found %d unique protein IDs in DE results", length(protein_ids)))
  
  # Clean protein IDs for matching
  # Remove isoform suffixes (-1, -2, etc.) and take first protein from groups
  cleaned_protein_ids <- protein_ids |>
    purrr::map_chr(\(x) {
      # Split by semicolon and take first protein
      first_protein <- stringr::str_split(x, ";")[[1]][1]
      # Remove isoform suffix
      stringr::str_replace(first_protein, "-\\d+$", "")
    })
  
  # Create matching table
  protein_mapping <- data.frame(
    original_id = protein_ids,
    cleaned_id = cleaned_protein_ids,
    stringsAsFactors = FALSE
  )
  
  log_info("Performing exact matching...")
  
  # Perform exact matching first
  exact_matches <- protein_mapping |>
    dplyr::left_join(
      uniprot_annotations,
      by = setNames(uniprot_id_column, "cleaned_id")
    ) |>
    dplyr::filter(!is.na(.data[[uniprot_id_column]]))
  
  log_info(sprintf("Exact matches found: %d", nrow(exact_matches)))
  
  # For unmatched proteins, try partial matching strategies
  unmatched_proteins <- protein_mapping |>
    dplyr::anti_join(exact_matches, by = "original_id")
  
  log_info(sprintf("Proteins requiring fuzzy matching: %d", nrow(unmatched_proteins)))
  
  fuzzy_matches <- data.frame()
  
  if (nrow(unmatched_proteins) > 0) {
    log_info("Performing fuzzy matching for remaining proteins...")
    
    # Strategy 1: Remove version numbers and try matching
    version_cleaned <- unmatched_proteins |>
      dplyr::mutate(
        version_cleaned = stringr::str_replace(.data$cleaned_id, "\\.\\d+$", "")
      ) |>
      dplyr::left_join(
        uniprot_annotations,
        by = setNames(uniprot_id_column, "version_cleaned")
      ) |>
      dplyr::filter(!is.na(.data[[uniprot_id_column]]))
    
    if (nrow(version_cleaned) > 0) {
      fuzzy_matches <- dplyr::bind_rows(fuzzy_matches, version_cleaned)
      log_info(sprintf("Version-based matches found: %d", nrow(version_cleaned)))
    }
    
    # Strategy 2: Substring matching for remaining proteins
    still_unmatched <- unmatched_proteins |>
      dplyr::anti_join(fuzzy_matches, by = "original_id")
    
    if (nrow(still_unmatched) > 0) {
      log_info("Attempting substring matching for remaining proteins...")
      
      # This is computationally expensive, so limit to reasonable number
      if (nrow(still_unmatched) <= 1000) {
        substring_matches <- still_unmatched |>
          dplyr::rowwise() |>
          dplyr::mutate(
            matched_uniprot = {
              # Look for partial matches in UniProt IDs
              matches <- stringr::str_detect(uniprot_annotations[[uniprot_id_column]], .data$cleaned_id)
              if (any(matches)) {
                uniprot_annotations[[uniprot_id_column]][matches][1]  # Take first match
              } else {
                NA_character_
              }
            }
          ) |>
          dplyr::ungroup() |>
          dplyr::filter(!is.na(.data$matched_uniprot)) |>
          dplyr::left_join(
            uniprot_annotations,
            by = setNames(uniprot_id_column, "matched_uniprot")
          )
        
        if (nrow(substring_matches) > 0) {
          fuzzy_matches <- dplyr::bind_rows(fuzzy_matches, substring_matches)
          log_info(sprintf("Substring matches found: %d", nrow(substring_matches)))
        }
      } else {
        log_warn(sprintf("Skipping substring matching for %d proteins (too many for efficient processing)", nrow(still_unmatched)))
      }
    }
  }
  
  # Combine all matches
  all_matches <- dplyr::bind_rows(exact_matches, fuzzy_matches)
  
  # Create final annotated results
  annotated_results <- protein_mapping |>
    dplyr::left_join(all_matches |> dplyr::select(-cleaned_id), by = "original_id")
  
  # Add gene names if available
  if (!is.null(gene_names_column)) {
    annotated_results <- annotated_results |>
      dplyr::mutate(
        gene_name = if (gene_names_column %in% names(annotated_results)) {
          purrr::map_chr(.data[[gene_names_column]], \(x) {
            if (is.na(x) || x == "") {
              NA_character_
            } else {
              # Extract first gene name from semicolon/space separated list
              stringr::str_split(x, "[; ]")[[1]][1]
            }
          })
        } else {
          NA_character_
        }
      )
  }
  
  # Calculate matching statistics
  total_proteins <- length(protein_ids)
  matched_proteins <- sum(!is.na(annotated_results[[uniprot_id_column]]))
  exact_match_count <- nrow(exact_matches)
  fuzzy_match_count <- nrow(fuzzy_matches)
  unmatched_count <- total_proteins - matched_proteins
  
  match_statistics <- list(
    total_proteins = total_proteins,
    matched_proteins = matched_proteins,
    exact_matches = exact_match_count,
    fuzzy_matches = fuzzy_match_count,
    unmatched_proteins = unmatched_count,
    match_rate = round(matched_proteins / total_proteins * 100, 2),
    exact_match_rate = round(exact_match_count / total_proteins * 100, 2),
    fuzzy_match_rate = round(fuzzy_match_count / total_proteins * 100, 2)
  )
  
  # Get list of unmatched proteins for troubleshooting
  unmatched_protein_list <- annotated_results |>
    dplyr::filter(is.na(.data[[uniprot_id_column]])) |>
    dplyr::pull(.data$original_id)
  
  log_info("=== Annotation matching completed ===")
  log_info(sprintf("Total proteins: %d", match_statistics$total_proteins))
  log_info(sprintf("Matched proteins: %d (%.1f%%)", match_statistics$matched_proteins, match_statistics$match_rate))
  log_info(sprintf("Exact matches: %d (%.1f%%)", match_statistics$exact_matches, match_statistics$exact_match_rate))
  log_info(sprintf("Fuzzy matches: %d (%.1f%%)", match_statistics$fuzzy_matches, match_statistics$fuzzy_match_rate))
  log_info(sprintf("Unmatched proteins: %d (%.1f%%)", match_statistics$unmatched_proteins, 100 - match_statistics$match_rate))
  
  return(list(
    annotated_de_results = annotated_results,
    match_statistics = match_statistics,
    unmatched_proteins = unmatched_protein_list
  ))
}


# ----------------------------------------------------------------------------
# detectEnsemblIds
# ----------------------------------------------------------------------------
#' Detect if Protein IDs are Ensembl Format
#'
#' @description Detects whether a vector of protein IDs contains Ensembl protein IDs
#' based on characteristic Ensembl naming patterns (e.g., ENSP, ENSOARP, etc.)
#'
#' @param protein_ids Character vector of protein identifiers
#'
#' @return List with:
#' \itemize{
#'   \item is_ensembl: Logical, TRUE if Ensembl IDs detected
#'   \item ensembl_prefix: Character, the detected Ensembl prefix (e.g., "ENSP", "ENSOARP")
#'   \item detection_rate: Numeric, proportion of IDs matching Ensembl pattern
#' }
#'
#' @details
#' Ensembl protein IDs follow the pattern: `ENS[SPECIES]P[0-9]+`
#' Examples: ENSP00000123456 (human), ENSOARP00020006936 (sheep), ENSMUSP00000123456 (mouse)
#'
#' @export
detectEnsemblIds <- function(protein_ids) {
  # Remove NAs and empty strings
  valid_ids <- protein_ids[!is.na(protein_ids) & protein_ids != ""]
  
  if (length(valid_ids) == 0) {
    return(list(is_ensembl = FALSE, ensembl_prefix = NA_character_, detection_rate = 0))
  }
  
  # Pattern for Ensembl protein IDs: ENS followed by optional species code, then P, then numbers
  # Can have version suffix like .1, .2
  ensembl_pattern <- "^ENS[A-Z]*P[0-9]+(\\.[0-9]+)?$"
  
  matches <- grepl(ensembl_pattern, valid_ids)
  detection_rate <- sum(matches) / length(valid_ids)
  
  # Consider it Ensembl if >80% of IDs match the pattern
  is_ensembl <- detection_rate > 0.8
  
  # Extract common prefix if Ensembl
  ensembl_prefix <- NA_character_
  if (is_ensembl && length(valid_ids) > 0) {
    # Extract the prefix (everything before the numbers)
    first_match <- valid_ids[matches][1]
    prefix_match <- regmatches(first_match, regexpr("^ENS[A-Z]*P", first_match))
    if (length(prefix_match) > 0) {
      ensembl_prefix <- prefix_match[1]
    }
  }
  
  return(list(
    is_ensembl = is_ensembl,
    ensembl_prefix = ensembl_prefix,
    detection_rate = detection_rate
  ))
}


# ----------------------------------------------------------------------------
# taxonIdToGprofilerOrganism
# ----------------------------------------------------------------------------
#' Map NCBI Taxon ID to gprofiler2 Organism Code
#'
#' @description Converts NCBI taxonomy IDs to gprofiler2 organism codes
#' for use with gprofiler2::gconvert() function
#'
#' @param taxon_id Numeric NCBI taxonomy identifier
#'
#' @return Character string with gprofiler2 organism code
#'
#' @details
#' Supported organisms include common model organisms and livestock species.
#' If the taxon_id is not in the mapping, the function will stop with an error
#' listing supported organisms.
#'
#' @examples
#' \dontrun{
#' taxonIdToGprofilerOrganism(9606)  # Returns "hsapiens"
#' taxonIdToGprofilerOrganism(9940)  # Returns "oaries"
#' }
#'
#' @export
taxonIdToGprofilerOrganism <- function(taxon_id) {
  # Comprehensive mapping of common organisms
  taxon_mapping <- c(
    "9606" = "hsapiens",        # Human
    "10090" = "mmusculus",       # Mouse
    "10116" = "rnorvegicus",     # Rat
    "9940" = "oaries",           # Sheep
    "9913" = "btaurus",          # Cow
    "9031" = "ggallus",          # Chicken
    "9823" = "sscrofa",          # Pig
    "7955" = "drerio",           # Zebrafish
    "6239" = "celegans",         # C. elegans
    "7227" = "dmelanogaster",    # Fruit fly
    "4932" = "scerevisiae",      # Yeast
    "3702" = "athaliana",        # Arabidopsis
    "9615" = "cfamiliaris",      # Dog
    "9796" = "ecaballus",        # Horse
    "9544" = "mmulatta"          # Rhesus macaque
  )
  
  taxon_str <- as.character(taxon_id)
  
  if (taxon_str %in% names(taxon_mapping)) {
    return(taxon_mapping[[taxon_str]])
  } else {
    supported_organisms <- paste(
      sprintf("%s (%s)", names(taxon_mapping), taxon_mapping),
      collapse = "\n  "
    )
    stop(sprintf(
      "Taxon ID %s not supported for gprofiler2 conversion.\nSupported organisms:\n  %s\n\nFalling back to direct UniProt query.",
      taxon_id,
      supported_organisms
    ))
  }
}


# ----------------------------------------------------------------------------
# convertEnsemblToUniprot
# ----------------------------------------------------------------------------
#' Convert Ensembl Protein IDs to UniProt IDs
#'
#' @description Uses gprofiler2::gconvert() to convert Ensembl protein IDs
#' to UniProt/SwissProt accessions. Automatically strips version suffixes
#' (e.g., .1, .2) before conversion and maps results back to original IDs.
#'
#' @param ensembl_ids Character vector of Ensembl protein IDs (with or without version suffixes)
#' @param organism_code Character string with gprofiler2 organism code (e.g., "hsapiens", "oaries")
#'
#' @return Data frame with columns:
#' \itemize{
#'   \item original_id: Original Ensembl protein ID (with version suffix if provided)
#'   \item converted_id: Converted UniProt accession (NA if no match)
#'   \item conversion_status: "success" or "failed"
#'   \item n_targets: Number of UniProt IDs found for this Ensembl ID
#' }
#'
#' @details
#' This function uses gprofiler2's gene identifier conversion service.
#' 
#' **Version Suffix Handling:**
#' - Ensembl IDs with version suffixes (e.g., ENSOARP00020006936.1) are automatically 
#'   stripped before conversion
#' - Results are mapped back to original IDs with versions preserved
#' 
#' **Multiple Matches:**
#' - Some Ensembl IDs may map to multiple UniProt IDs
#' - In such cases, the first match is returned (gprofiler2 orders by relevance)
#'
#' @examples
#' \dontrun{
#' # Human IDs
#' convertEnsemblToUniprot(
#'   c("ENSP00000269305", "ENSP00000256078"),
#'   "hsapiens"
#' )
#' 
#' # Sheep IDs with version suffixes
#' convertEnsemblToUniprot(
#'   c("ENSOARP00020006936.1", "ENSOARP00020015587.2"),
#'   "oaries"
#' )
#' }
#'
#' @export
convertEnsemblToUniprot <- function(ensembl_ids, organism_code) {
  # Check if gprofiler2 is available
  if (!requireNamespace("gprofiler2", quietly = TRUE)) {
    stop("Package 'gprofiler2' is required for Ensembl ID conversion but is not installed.\nPlease install it with: install.packages('gprofiler2')")
  }
  
  cat(sprintf("Converting %d Ensembl IDs to UniProt using gprofiler2...\n", length(ensembl_ids)))
  
  # ✅ FIX: Strip version suffixes (.1, .2, etc.) before conversion
  # gprofiler2 doesn't recognize Ensembl IDs with version suffixes
  cat("Stripping version suffixes from Ensembl IDs (e.g., .1, .2)...\n")
  
  # Create mapping: original_with_version → stripped_version
  version_mapping <- data.frame(
    original_with_version = ensembl_ids,
    stripped_version = gsub("\\.\\d+$", "", ensembl_ids),  # Remove .X at end
    stringsAsFactors = FALSE
  )
  
  # Get unique stripped IDs (multiple versioned IDs may map to same stripped ID)
  unique_stripped <- unique(version_mapping$stripped_version)
  
  cat(sprintf("Stripped versions: %d original IDs → %d unique stripped IDs\n", 
              length(ensembl_ids), length(unique_stripped)))
  cat(sprintf("Example: %s → %s\n", 
              version_mapping$original_with_version[1], 
              version_mapping$stripped_version[1]))
  
  # Call gprofiler2::gconvert with STRIPPED IDs
  tryCatch({
    cat("Calling gprofiler2::gconvert()...\n")
    conversion_result <- gprofiler2::gconvert(
      query = unique_stripped,
      organism = organism_code,
      target = "UNIPROTSPTREMBL",  # UniProt Swiss-Prot + TrEMBL
      mthreshold = Inf,             # Return all matches
      filter_na = FALSE             # Keep non-matching IDs
    )
    
    cat(sprintf("✓ gprofiler2 returned %d conversion results\n", 
                if(is.null(conversion_result)) 0 else nrow(conversion_result)))
    
    if (!is.null(conversion_result) && nrow(conversion_result) > 0) {
      cat("First 3 conversions from gprofiler2:\n")
      print(head(conversion_result[, c("input", "target", "name")], 3))
    }
    
    # Create mapping table for STRIPPED → UniProt
    stripped_to_uniprot <- data.frame(
      stripped_version = unique_stripped,
      converted_id = NA_character_,
      conversion_status = "failed",
      n_targets = 0,
      stringsAsFactors = FALSE
    )
    
    if (!is.null(conversion_result) && nrow(conversion_result) > 0) {
      # For each stripped ID, get the first matching UniProt ID
      for (i in seq_len(nrow(stripped_to_uniprot))) {
        stripped <- stripped_to_uniprot$stripped_version[i]
        
        # Find matches in conversion result
        matches <- conversion_result[conversion_result$input == stripped, ]
        
        if (nrow(matches) > 0) {
          # Take first match (gprofiler2 orders by relevance)
          stripped_to_uniprot$converted_id[i] <- matches$target[1]
          stripped_to_uniprot$conversion_status[i] <- "success"
          stripped_to_uniprot$n_targets[i] <- nrow(matches)
        }
      }
    }
    
    # ✅ FIX: Map back to ORIGINAL IDs (with versions)
    # Join: original_with_version → stripped_version → converted_id
    final_mapping <- merge(
      version_mapping,
      stripped_to_uniprot,
      by = "stripped_version",
      all.x = TRUE
    )
    
    # Rename for consistency with expected output
    final_mapping <- final_mapping[, c("original_with_version", "converted_id", "conversion_status", "n_targets")]
    names(final_mapping)[1] <- "original_id"
    
    success_count <- sum(final_mapping$conversion_status == "success", na.rm = TRUE)
    cat(sprintf("✓ Successfully converted %d/%d Ensembl IDs to UniProt IDs (%.1f%%)\n", 
                success_count, length(ensembl_ids), 
                100 * success_count / length(ensembl_ids)))
    
    if (success_count < length(ensembl_ids)) {
      failed_count <- length(ensembl_ids) - success_count
      cat(sprintf("⚠ Warning: %d Ensembl IDs could not be converted (will use original ID with NA annotations)\n",
                  failed_count))
      if (failed_count <= 10) {
        cat("Failed IDs:\n")
        print(final_mapping$original_id[final_mapping$conversion_status != "success"])
      }
    }
    
    return(final_mapping)
    
  }, error = function(e) {
    cat(sprintf("✗ ERROR in gprofiler2::gconvert: %s\n", e$message))
    warning(sprintf("gprofiler2::gconvert failed: %s", e$message))
    
    # Return empty mapping on error
    return(data.frame(
      original_id = ensembl_ids,
      converted_id = NA_character_,
      conversion_status = "error",
      n_targets = 0,
      stringsAsFactors = FALSE
    ))
  })
}


# ----------------------------------------------------------------------------
# getUniprotAnnotationsFull
# ----------------------------------------------------------------------------
#' Get UniProt Annotations from Actual Data Proteins
#'
#' @description Creates a comprehensive UniProt annotation lookup table using
#' protein IDs from the actual imported data. This function extracts unique protein
#' accessions by splitting protein groups and retrieves their annotations efficiently.
#' 
#' Automatically detects and converts Ensembl protein IDs to UniProt IDs using gprofiler2
#' before querying the UniProt API.
#'
#' @param data_tbl Data frame containing the actual imported quantitative data
#' @param protein_id_column Character string specifying the protein ID column in data_tbl
#' @param cache_dir Character string path to cache directory for storing annotations
#' @param taxon_id Numeric taxonomy ID for the organism
#' @param chunk_size Numeric size of chunks for batch processing (default: 25)
#'
#' @return Data frame with comprehensive UniProt annotations for proteins in the dataset
#'
#' @details
#' This function:
#' \itemize{
#'   \item Extracts unique protein groups from actual quantitative data
#'   \item Splits protein groups by semicolon to get individual protein IDs
#'   \item Detects if IDs are Ensembl format and converts to UniProt if needed
#'   \item Creates master list of unique individual proteins
#'   \item Uses efficient caching and retrieval logic for annotations
#'   \item Does NOT modify the original quantitative data
#' }
#'
#' This approach is much more efficient than annotating all FASTA proteins because:
#' \itemize{
#'   \item Only annotates proteins actually present in the dataset
#'   \item Handles protein groups properly (A;B;C → A, B, C)
#'   \item Significantly reduces API calls and processing time
#'   \item Supports both UniProt and Ensembl protein IDs
#' }
#'
#' @examples
#' \dontrun{
#' # Create efficient UniProt annotations from actual data
#' uniprot_lookup <- getUniprotAnnotationsFull(
#'   data_tbl = imported_data,
#'   protein_id_column = "Protein.Group",
#'   cache_dir = "cache/uniprot_annotations",
#'   taxon_id = 9606
#' )
#' }
#'
#' @export
.extractProteinIdFromHeader <- function(id_string) {
  # 1. Try to match standard UniProt/SwissProt header: sp|ACCESSION|...
  # or tr|ACCESSION|...
  if (grepl("^(?:sp|tr)\\|[^|]+\\|", id_string)) {
    return(sub("^(?:sp|tr)\\|([^|]+)\\|.*", "\\1", id_string))
  }
  
  # 2. Try to match generic pipe format: generic|ACCESSION|...
  if (grepl("^[^|]+\\|[^|]+\\|", id_string)) {
    return(sub("^[^|]+\\|([^|]+)\\|.*", "\\1", id_string))
  }
  
  # 3. Take the first whitespace-separated token
  first_token <- strsplit(id_string, "\\s+")[[1]][1]
  
  # 4. Handle Colon-separated prefixes (e.g. TREMBL:A2A5Y0, ENSEMBL:ENSP00000...)
  # Only if the token contains a colon
  if (grepl(":", first_token)) {
    # Split by colon and take the last part
    # Use tail(..., 1) to get the last element
    parts <- strsplit(first_token, ":")[[1]]
    return(tail(parts, 1))
  }
  
  return(first_token)
}

#' @export
getUniprotAnnotationsFull <- function(data_tbl,
                                     protein_id_column,
                                     cache_dir,
                                     taxon_id = 9606,
                                     chunk_size = 25,
                                     progress_callback = NULL) {
  
  message("=== DEBUG66: Entering getUniprotAnnotationsFull ===")
  message(sprintf("   Data rows: %d", nrow(data_tbl)))
  message(sprintf("   Protein column: %s", protein_id_column))
  message(sprintf("   Taxon ID: %d", taxon_id))
  
  cat("=== Starting getUniprotAnnotationsFull ===\n")
  
  # Validate inputs
  if (is.null(data_tbl) || !is.data.frame(data_tbl)) {
    stop("data_tbl must be a non-null data frame")
  }
  
  if (nrow(data_tbl) == 0) {
    stop("data_tbl cannot be empty")
  }
  
  if (!protein_id_column %in% names(data_tbl)) {
    stop(sprintf("Protein ID column '%s' not found. Available columns: %s", 
                 protein_id_column, paste(names(data_tbl), collapse = ", ")))
  }
  
  cat(sprintf("Processing quantitative data with %d rows for UniProt annotation\n", nrow(data_tbl)))
  cat(sprintf("Using protein ID column: %s\n", protein_id_column))
  cat(sprintf("Target organism taxon ID: %d\n", taxon_id))
  
  # Extract unique protein groups from actual data
  cat("Extracting unique protein groups from actual data...\n")
  
  # Get unique protein groups from the imported data
  raw_protein_groups <- unique(data_tbl[[protein_id_column]])
  
  cat(sprintf("Found %d unique protein groups in data\n", length(raw_protein_groups)))
  cat(sprintf("First 5 protein groups: %s\n", paste(head(raw_protein_groups, 5), collapse = ", ")))
  
  # Split protein groups by semicolon to get individual proteins
  cat("Splitting protein groups to individual proteins...\n")
  
  all_individual_proteins <- character()
  
  for (group in raw_protein_groups) {
    if (!is.na(group) && group != "") {
      # Split by semicolon
      individual_proteins <- strsplit(group, ";")[[1]]
      # Clean up (remove whitespace, empty strings)
      individual_proteins <- trimws(individual_proteins)
      individual_proteins <- individual_proteins[individual_proteins != "" & !is.na(individual_proteins)]
      # Add to master list
      all_individual_proteins <- c(all_individual_proteins, individual_proteins)
    }
  }
  
  # Get unique individual proteins
  unique_proteins <- unique(all_individual_proteins)
  
  # Clean up IDs (extract from FASTA headers if present)
  cleaned_ids <- vapply(unique_proteins, .extractProteinIdFromHeader, FUN.VALUE = character(1), USE.NAMES = FALSE)

  # Remove potential isoform suffixes for broader matching
  cleaned_proteins <- unique(gsub("-\\d+$", "", cleaned_ids))
  
  cat(sprintf("Expanded %d protein groups into %d unique individual proteins\n", 
                   length(raw_protein_groups), length(cleaned_proteins)))
  cat(sprintf("First 10 individual proteins: %s\n", paste(head(cleaned_proteins, 10), collapse = ", ")))
  
  message(sprintf("   DEBUG66: Extracted %d unique proteins", length(cleaned_proteins)))
  message("   DEBUG66: First 10 proteins:")
  print(head(cleaned_proteins, 10))
  
  # === ENSEMBL ID CONVERSION LOGIC ===
  # Detect if protein IDs are Ensembl format and convert to UniProt if necessary
  ensembl_detection <- detectEnsemblIds(cleaned_proteins)
  ensembl_mapping <- NULL  # Initialize mapping variable
  proteins_to_query <- cleaned_proteins  # Default: use original proteins
  
  if (ensembl_detection$is_ensembl) {
    cat("\n=== ENSEMBL ID CONVERSION ===\n")
    cat(sprintf("Ensembl protein IDs detected (%.1f%% match rate)\n", 
                ensembl_detection$detection_rate * 100))
    cat(sprintf("Ensembl prefix: %s\n", ensembl_detection$ensembl_prefix))
    cat("Converting to UniProt IDs using gprofiler2...\n")
    
    tryCatch({
      # Map taxon ID to gprofiler2 organism
      organism_code <- taxonIdToGprofilerOrganism(taxon_id)
      cat(sprintf("Organism code for gprofiler2: %s\n", organism_code))
      
      # Convert Ensembl → UniProt
      conversion_result <- convertEnsemblToUniprot(cleaned_proteins, organism_code)
      
      # Filter to successfully converted IDs for UniProt query
      successfully_converted <- conversion_result[
        conversion_result$conversion_status == "success" & 
        !is.na(conversion_result$converted_id), 
      ]
      
      if (nrow(successfully_converted) > 0) {
        # ✅ CRITICAL: Extract UNIPROT IDs (not Ensembl IDs!) from conversion result
        proteins_to_query <- successfully_converted$converted_id
        ensembl_mapping <- conversion_result  # Save for later merge
        
        cat(sprintf("\n✓ Conversion summary:\n"))
        cat(sprintf("  - Successfully converted: %d/%d (%.1f%%)\n", 
                    nrow(successfully_converted), 
                    nrow(conversion_result),
                    100 * nrow(successfully_converted) / nrow(conversion_result)))
        cat(sprintf("  - Failed conversions: %d (will have NA annotations)\n", 
                    nrow(conversion_result) - nrow(successfully_converted)))
        cat(sprintf("  - UniProt IDs extracted for API query: %d\n", length(proteins_to_query)))
        cat(sprintf("\n✓ CRITICAL CHECK - First 5 IDs to be submitted to UniProt API:\n"))
        cat(sprintf("  Type: %s\n", if(any(grepl("^ENS", head(proteins_to_query, 5)))) "❌ ENSEMBL (WRONG!)" else "✓ UniProt (CORRECT)"))
        print(head(proteins_to_query, 5))
        cat("=== END ENSEMBL CONVERSION ===\n\n")
      } else {
        cat("✗ ERROR: gprofiler2 conversion returned 0 results!\n")
        warning("No Ensembl IDs could be converted. Falling back to direct query with original IDs.")
        proteins_to_query <- cleaned_proteins
        ensembl_mapping <- NULL
      }
      
    }, error = function(e) {
      warning(sprintf("Ensembl conversion failed: %s. Falling back to direct UniProt query.", e$message))
      cat(sprintf("Error details: %s\n", e$message))
      proteins_to_query <- cleaned_proteins
      ensembl_mapping <- NULL
    })
  } else {
    cat("UniProt IDs detected (or non-Ensembl format). Proceeding with direct query.\n")
  }
  # === END ENSEMBL CONVERSION ===
  
  # Create a temporary data frame for getUniProtAnnotation
  # This function expects protein IDs in a standard format
  # Use proteins_to_query which may be converted UniProt IDs or original IDs
  temp_protein_table <- data.frame(
    Protein.Ids = proteins_to_query,
    stringsAsFactors = FALSE
  )
  
  cat("\n=== FINAL CHECK BEFORE UniProt API SUBMISSION ===\n")
  cat(sprintf("Temporary table dimensions: %d rows x %d cols\n", nrow(temp_protein_table), ncol(temp_protein_table)))
  cat(sprintf("First 5 IDs in temp_protein_table$Protein.Ids:\n"))
  print(head(temp_protein_table$Protein.Ids, 5))
  cat(sprintf("ID Type Check: %s\n", 
              if(any(grepl("^ENS", head(temp_protein_table$Protein.Ids, 10)))) {
                "❌ STILL ENSEMBL IDS - BUG IN CONVERSION LOGIC!"
              } else {
                "✓ UniProt IDs - Ready for API"
              }))
  cat("=== END FINAL CHECK ===\n\n")
  
  # Create cache directory if needed
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
    cat(sprintf("Created cache directory: %s\n", cache_dir))
  }
  
  # Call getUniprotAnnotations with our optimized protein list
  cat("Calling getUniprotAnnotations for comprehensive annotation retrieval...\n")
  message("   DEBUG66: About to call getUniprotAnnotations...")
  message(sprintf("   DEBUG66: temp_protein_table has %d rows", nrow(temp_protein_table)))
  
  tryCatch({
    uniprot_annotations <- getUniprotAnnotations(
      input_tbl = temp_protein_table,
      cache_dir = cache_dir,
      taxon_id = taxon_id,
      progress_callback = progress_callback
    )
    
    cat(sprintf("Successfully retrieved UniProt annotations for %d proteins\n", nrow(uniprot_annotations)))
    
    # === MERGE BACK ORIGINAL ENSEMBL IDs IF CONVERSION WAS PERFORMED ===
    if (!is.null(ensembl_mapping)) {
      cat("\n=== MERGING ENSEMBL IDs WITH UNIPROT ANNOTATIONS ===\n")
      cat("Preserving original Ensembl IDs alongside UniProt annotations...\n")
      
      # The uniprot_annotations table has Entry column with UniProt IDs
      # We need to map these back to original Ensembl IDs
      
      # Add Original_ID column to preserve Ensembl IDs
      if ("Entry" %in% names(uniprot_annotations)) {
        # Create a lookup from converted_id to original_id
        ensembl_lookup <- ensembl_mapping[, c("original_id", "converted_id")]
        names(ensembl_lookup) <- c("Original_Ensembl_ID", "Entry")
        
        # Merge annotations with Ensembl mapping
        uniprot_annotations <- merge(
          uniprot_annotations,
          ensembl_lookup,
          by = "Entry",
          all.x = TRUE
        )
        
        # For any annotations without Ensembl mapping, use Entry as Original_ID
        if ("Original_Ensembl_ID" %in% names(uniprot_annotations)) {
          uniprot_annotations$Original_Ensembl_ID <- ifelse(
            is.na(uniprot_annotations$Original_Ensembl_ID),
            uniprot_annotations$Entry,
            uniprot_annotations$Original_Ensembl_ID
          )
        }
        
        cat(sprintf("Added Original_Ensembl_ID column to preserve Ensembl identifiers\n"))
        cat(sprintf("Final annotation table: %d rows, %d columns\n", 
                    nrow(uniprot_annotations), ncol(uniprot_annotations)))
      } else {
        warning("Entry column not found in annotations. Cannot map back to Ensembl IDs.")
      }
      
      # Add metadata flag
      attr(uniprot_annotations, "ensembl_conversion_applied") <- TRUE
      attr(uniprot_annotations, "ensembl_conversion_success_rate") <- 
        sum(ensembl_mapping$conversion_status == "success") / nrow(ensembl_mapping)
      
      cat("=== ENSEMBL ID MERGING COMPLETE ===\n\n")
    } else {
      # No Ensembl conversion was performed
      attr(uniprot_annotations, "ensembl_conversion_applied") <- FALSE
    }
    # === END ENSEMBL ID MERGING ===
    
    # Add some metadata about the annotation process
    attr(uniprot_annotations, "annotation_source") <- "DATA_OPTIMIZED"
    attr(uniprot_annotations, "original_protein_groups") <- length(raw_protein_groups)
    attr(uniprot_annotations, "unique_proteins_processed") <- length(cleaned_proteins)
    attr(uniprot_annotations, "taxon_id") <- taxon_id
    attr(uniprot_annotations, "processing_timestamp") <- Sys.time()
    
    cat("=== getUniprotAnnotationsFull completed successfully ===\n")
    
    return(uniprot_annotations)
    
  }, error = function(e) {
    cat("=== DEBUG66: ERROR in getUniprotAnnotationsFull ===\n")
    cat(sprintf("Error message: %s\n", e$message))
    cat(sprintf("Error class: %s\n", paste(class(e), collapse = ", ")))
    cat("Full error details:\n")
    print(str(e))
    cat(sprintf("Error in getUniProtAnnotation: %s\n", e$message))
    
    # Return a minimal data frame if annotation fails
    cat("Returning minimal annotation table due to error\n")
    
    minimal_annotations <- data.frame(
      Protein.Ids = cleaned_proteins,
      UNIPROT_GENENAME = NA_character_,
      stringsAsFactors = FALSE
    )
    
    return(minimal_annotations)
  })
}


# ----------------------------------------------------------------------------
# getFastaFields
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @export
getFastaFields <- function(string, pattern) {

  field_found <- str_detect( {{string}}, paste0(pattern, "="))

  extract_data <- NA_character_
  if( field_found ) {

    extract_data <- str_replace_all( {{string}},
                                     paste0("(.*)",
                                            pattern,
                                            "=(.*?)(\\s..=.*|$)"), "\\2")
  }

  case_when(  field_found ~ extract_data,
              TRUE ~ NA_character_  )
}


# ----------------------------------------------------------------------------
# parseFastaObject
# ----------------------------------------------------------------------------
#' Parse FASTA object from seqinr
#' @description parse_fasta_object: Parse FASTA headers
#' @param aa_seq AAStringSet object, output from running seqinr
#' @return A table containing the protein evidence, isoform number, uniprot accession without isoform number in the uniprot_acc column, gene name
#' @export
parseFastaObject <- function(aa_seq ) {

  accession_tab <-  data.frame( header=names(aa_seq)) %>%
    separate( header, into=c("db", "uniprot_acc", "description"), sep="\\|") %>%
    mutate( uniprot_id = str_replace( description, "(.*?)\\s(.*)", "\\1" ) ) %>%
    mutate( OS = purrr::map_chr(description, ~getFastaFields(., "OS")))  %>%
    mutate( OX = purrr::map_int(description, ~as.integer(getFastaFields(., "OX")))) %>%
    mutate( GN = purrr::map_chr(description, ~getFastaFields(., "GN"))) %>%
    mutate( GN = ifelse( is.na(GN), "", GN)) %>%
    mutate( PE = purrr::map_int(description, ~as.integer(getFastaFields(., "PE")))) %>%
    mutate( SV = purrr::map_int(description, ~as.integer(getFastaFields(., "SV")))) %>%
    dplyr::select(-description) %>%
    dplyr::rename( species = "OS",
                   tax_id = "OX",
                   gene_name = "GN",
                   protein_evidence = "PE",
                   sequence_version = "SV")

  acc_detail_tab <- accession_tab %>%
    mutate( is_isoform = case_when( str_detect( uniprot_acc, "-\\d+") ~ "Isoform",
                                    TRUE ~ "Canonical") ) %>%
    mutate (isoform_num = case_when ( is_isoform == "Isoform" ~ str_replace_all( uniprot_acc,
                                                                                 "(.*)(-)(\\d{1,})",
                                                                                 "\\3") %>%
                                        as.numeric,
                                      is_isoform == "Canonical" ~ 0,
                                      TRUE ~ NA_real_ ) ) %>%
    mutate( cleaned_acc = cleanIsoformNumber(uniprot_acc)) %>%
    mutate( protein_evidence  = factor(protein_evidence, levels =1:5 )) %>%
    mutate( status = factor( db, levels =c( "sp", "tr"), labels=c("reviewed", "unreviewed"))) %>%
    mutate( is_isoform = factor(is_isoform, levels =c("Canonical", "Isoform")))

  return(acc_detail_tab )

}


# ----------------------------------------------------------------------------
# parseFastaFile
# ----------------------------------------------------------------------------
#'Parse the headers of a Uniprot FASTA file and extract the headers and sequences into a data frame
#'Use seqinr object instead as it seems to be a lot faster to run substring
#' @param path to input faster file with header format described in https://www.uniprot.org/help/fasta-headers
#' @return A table containing the following columns:
#' db  sp for Swiss-Prot, tr for TrEMBL
#' uniprot_acc Uniprot Accession
#' uniprot_id  Uniprot ID
#' species     Species
#' tax_id      Taxonomy ID
#' gene_name   Gene symbol
#' protein_evidence 1 to 5, the lower the value, the more evidence that supports the existence of this protein
#' sequence_version Sequence version
#' is_isoform  Is it a protein isoform (not the canonical form)
#' isoform_num     Isoform number.
#' cleaned_acc Cleaned accession without isoform number.
#' status  Reviewed or unreviewed.
#' seq     Amino acid sequence.
#' seq_length      Sequence length (integer).
#' @export
parseFastaFile <- function(fasta_file) {

  aa_seqinr <-  read.fasta( file = fasta_file,
                            seqtype="AA",
                            whole.header	=TRUE,
                            as.string=TRUE)

  acc_detail_tab <- parseFastaObject(aa_seqinr)

  names(aa_seqinr) <- str_match( names(aa_seqinr), "(sp|tr)\\|(.+?)\\|(.*)\\s+" )[,3]

  aa_seq_tbl <- acc_detail_tab %>%
    mutate(seq = map_chr( aa_seqinr, 1)) %>%
    mutate(seq_length = purrr::map_int(seq, str_length) )
}


# ----------------------------------------------------------------------------
# chooseBestPhosphositeAccession
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @export
chooseBestPhosphositeAccession <- function(input_tbl, acc_detail_tab, accessions_column, group_id) {

  # Join with FASTA data
  resolve_acc_joined <- input_tbl %>%
    dplyr::select( {{group_id}}, {{accessions_column}}, cleaned_peptide) %>%
    mutate( uniprot_acc = str_split( {{accessions_column}}, ";") ) %>%
    unnest( uniprot_acc )   %>%
    mutate( cleaned_acc = cleanIsoformNumber(uniprot_acc)) %>%
    left_join( acc_detail_tab,
               by=c("uniprot_acc" = "uniprot_acc",
                    "cleaned_acc" = "cleaned_acc") )
  
  # Detect which columns are available and add missing ones with defaults
  available_cols <- colnames(resolve_acc_joined)
  
  if (!"gene_name" %in% available_cols) {
    resolve_acc_joined <- resolve_acc_joined %>% mutate(gene_name = NA_character_)
  }
  
  if (!"protein_evidence" %in% available_cols) {
    resolve_acc_joined <- resolve_acc_joined %>% mutate(protein_evidence = 3)
  }
  
  if (!"status" %in% available_cols) {
    resolve_acc_joined <- resolve_acc_joined %>% mutate(status = "unknown")
  }
  
  if (!"is_isoform" %in% available_cols) {
    resolve_acc_joined <- resolve_acc_joined %>% mutate(is_isoform = "Canonical")
  }
  
  if (!"isoform_num" %in% available_cols) {
    resolve_acc_joined <- resolve_acc_joined %>% mutate(isoform_num = 0)
  }
  
  if (!"seq_length" %in% available_cols) {
    resolve_acc_joined <- resolve_acc_joined %>% mutate(seq_length = NA_integer_)
  }
  
  if (!"seq" %in% available_cols) {
    resolve_acc_joined <- resolve_acc_joined %>% mutate(seq = NA_character_)
  }
  
  # Now select, filter, and arrange using the complete set of columns
  resolve_acc_helper <- resolve_acc_joined %>%
    ## Just a sanity check that the peptide is actually in the sequence (skip if seq not available)
    {if ("seq" %in% colnames(.) && all(!is.na(.$seq))) filter(., str_detect( seq, cleaned_peptide  )) else .} %>%
    dplyr::select({{group_id}}, one_of(c( "uniprot_acc", "gene_name", "cleaned_acc",
                                          "protein_evidence", "status", "is_isoform", "isoform_num", "seq_length"  ))) %>%
    distinct %>%
    mutate(annotation_score = if_else(is.na(annotation_score), 0, annotation_score)) %>%
    arrange( {{group_id}}, desc(annotation_score), protein_evidence, status, is_isoform, desc(seq_length), isoform_num )

  # print( colnames(head(resolve_acc_helper)) )


  score_isoforms <- resolve_acc_helper %>%
    mutate( gene_name = ifelse( is.na(gene_name) | gene_name == "", "NA", gene_name)) %>%
    mutate(annotation_score = if_else(is.na(annotation_score), 0, annotation_score)) %>%
    group_by( {{group_id}},  gene_name ) %>%
    arrange( {{group_id}},  desc(annotation_score), protein_evidence,
             status, is_isoform, desc(seq_length), isoform_num, cleaned_acc )  %>%
    mutate(ranking = row_number()) %>%
    ungroup


  # print( colnames(head(score_isoforms)) )

  ## For each gene name find the uniprot_acc with the lowest ranking
  group_gene_names_and_uniprot_accs <- score_isoforms %>%
    distinct( {{group_id}}, gene_name, ranking ) %>%
    dplyr::filter( ranking == 1) %>%
    left_join( score_isoforms %>%
                 dplyr::select( {{group_id}}, ranking, gene_name, uniprot_acc),
               by = join_by( {{group_id}} == {{group_id}}
                             , ranking == ranking
                             , gene_name == gene_name ) )   %>%
    dplyr::select(-ranking)

  # %>%
  #   group_by({{group_id}}) %>%
  #   summarise( num_gene_names = n(),
  #              gene_names = paste( gene_name, collapse=":"),
  #              uniprot_acc = paste( uniprot_acc, collapse=":")) %>%
  #   ungroup() %>%
  #   mutate( is_unique = case_when( num_gene_names == 1 ~ "Unique",
  #                                  TRUE ~ "Multimapped"))


  return( group_gene_names_and_uniprot_accs )

}


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
  
  cat("\n\n╔═══════════════════════════════════════════════════════════════════════════╗\n")
  cat("║      ENTERING chooseBestProteinAccessionHelper (DEBUG66)              ║\n")
  cat("╚═══════════════════════════════════════════════════════════════════════════╝\n\n")
  
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
  resolve_acc_temp <- input_tbl |>
    dplyr::select( { { group_id } }, { { accessions_column } }) |>
    mutate(row_id_column_with_isoform = str_split({ { accessions_column } }, delim)) |>
    unnest( row_id_column_with_isoform ) |>
    mutate( !!sym(row_id_column) := cleanIsoformNumber( row_id_column_with_isoform)) |>
    dplyr::filter( !str_detect(!!sym(row_id_column), "REV__")) |>
    dplyr::filter( !str_detect(!!sym(row_id_column), "CON__"))
  
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
  cat(sprintf("   ⚠️  CRITICAL: Using 'delim' for recombining: '%s'\n", delim))
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
  cat("╔═══════════════════════════════════════════════════════════════════════════╗\n")
  cat("║      EXITING chooseBestProteinAccessionHelper (DEBUG66)               ║\n")
  cat("╚═══════════════════════════════════════════════════════════════════════════╝\n\n")

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
  
  cat("\n\n╔═══════════════════════════════════════════════════════════════════════════╗\n")
  cat("║      ENTERING rankProteinAccessionHelper (DEBUG66)                    ║\n")
  cat("╚═══════════════════════════════════════════════════════════════════════════╝\n\n")
  
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
  cat(sprintf("   ⚠️  CRITICAL: Using 'delim' for recombining: '%s'\n", delim))
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
  cat("╔═══════════════════════════════════════════════════════════════════════════╗\n")
  cat("║      EXITING rankProteinAccessionHelper (DEBUG66)                     ║\n")
  cat("╚═══════════════════════════════════════════════════════════════════════════╝\n\n")

  return(group_gene_names_and_uniprot_accs)

}


# ----------------------------------------------------------------------------
# processFastaFile
# ----------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export
processFastaFile <- function(fasta_file_path, uniprot_search_results = NULL, uniparc_search_results = NULL, fasta_meta_file, organism_name) {
  # Properly suppress all vroom messages
  withr::local_options(list(
    vroom.show_col_types = FALSE,
    vroom.show_progress = FALSE
  ))

  startsWith <- function(x, prefix) {
    substr(x, 1, nchar(prefix)) == prefix
  }

  parseFastaFileStandard <- function(fasta_file) {
    message("Reading FASTA file with seqinr...")
    utils::flush.console()

    aa_seqinr <- seqinr::read.fasta(file = fasta_file, seqtype = "AA",
                                    whole.header = TRUE, as.string = TRUE)
    headers <- names(aa_seqinr)
    total_entries <- length(headers)

    message(sprintf("\nProcessing %d FASTA entries...", total_entries))
    utils::flush.console()

    # Create a text progress bar
    pb <- utils::txtProgressBar(min = 0, max = total_entries, style = 3, width = 50)

    parsed_headers <- vector("list", length(headers))

    for(i in seq_along(headers)) {
      header <- headers[i]
      parsed_headers[[i]] <- {
        parts <- strsplit(substr(header, 2, nchar(header)), " ", fixed = TRUE)[[1]]
        id_parts <- strsplit(parts[1], "|", fixed = TRUE)[[1]]

        # Extract protein evidence level
        protein_evidence <- stringr::str_extract(header, "PE=[0-9]") |>
          stringr::str_extract("[0-9]") |>
          as.integer()

        # Determine status based on entry type
        status <- if(startsWith(header, ">sp|")) "reviewed" else "unreviewed"

        # Extract gene name (GN=)
        gene_name <- stringr::str_extract(header, "GN=\\S+") |>
          stringr::str_remove("GN=")

        # For entries without isoforms, set defaults
        is_isoform <- FALSE
        isoform_num <- 0L
        cleaned_acc <- id_parts[2]

        list(
          accession = id_parts[2],
          database_id = id_parts[2],
          cleaned_acc = cleaned_acc,
          gene_name = gene_name,
          protein = paste(parts[-1], collapse = " "),
          attributes = paste(parts[-1], collapse = " "),
          protein_evidence = protein_evidence,
          status = status,
          is_isoform = is_isoform,
          isoform_num = isoform_num,
          annotation_score = NA_real_
        )
      }

      # Update progress bar every 100 entries
      if(i %% 100 == 0 || i == total_entries) {
        utils::setTxtProgressBar(pb, i)
      }
    }

    close(pb)

    message("\nBinding rows and creating final table...")
    utils::flush.console()

    acc_detail_tab <- dplyr::bind_rows(parsed_headers)
    aa_seq_tbl <- acc_detail_tab |>
      dplyr::mutate(
        seq = purrr::map_chr(aa_seqinr, 1),
        seq_length = stringr::str_length(seq),
        description = headers
      )

    return(aa_seq_tbl)
  }

  parseFastaHeader <- function(header) {
    parts <- strsplit(substr(header, 2, nchar(header)), " ", fixed = TRUE)[[1]]
    id_parts <- strsplit(parts[1], "|", fixed = TRUE)[[1]]
    accession <- id_parts[2]
    attributes <- paste(parts[-1], collapse = " ")
    locus_tag <- stringr::str_extract(attributes, "(?<=\\[locus_tag=)[^\\]]+")
    protein <- stringr::str_extract(attributes, "(?<=\\[protein=)[^\\]]+")
    ncbi_refseq <- stringr::str_extract(attributes, "(?<=\\[protein_id=)WP_[^\\]]+")
    list(
      accession = accession,
      protein_id = locus_tag,
      protein = protein,
      ncbi_refseq = ncbi_refseq,
      attributes = attributes,
      annotation_score = NA_real_
    )
  }

  parseFastaFileNonStandard <- function(fasta_file) {
    message("Reading FASTA file with seqinr...")
    utils::flush.console()

    aa_seqinr <- seqinr::read.fasta(file = fasta_file, seqtype = "AA",
                                    whole.header = TRUE, as.string = TRUE)
    headers <- names(aa_seqinr)
    total_entries <- length(headers)

    message(sprintf("\nProcessing %d non-standard FASTA entries...", total_entries))
    utils::flush.console()

    # Create a text progress bar
    pb <- utils::txtProgressBar(min = 0, max = total_entries, style = 3, width = 50)

    parsed_headers <- vector("list", length(headers))

    for(i in seq_along(headers)) {
      header <- headers[i]
      parsed_headers[[i]] <- parseFastaHeader(header)

      # Update progress bar every 100 entries
      if(i %% 100 == 0 || i == total_entries) {
        utils::setTxtProgressBar(pb, i)
      }
    }

    close(pb)

    message("\nBinding rows and creating final table...")
    utils::flush.console()

    acc_detail_tab <- dplyr::bind_rows(parsed_headers)
    aa_seq_tbl <- acc_detail_tab |>
      dplyr::mutate(
        seq = purrr::map_chr(aa_seqinr, 1),
        seq_length = stringr::str_length(seq),
        description = headers
      )

    return(aa_seq_tbl)
  }

  matchAndUpdateDataFrames <- function(aa_seq_tbl, uniprot_search_results, uniparc_search_results, organism_name) {
    message("Matching and updating dataframes...")
    flush.console()

    uniprot_filtered <- uniprot_search_results |>
      dplyr::filter(Organism == organism_name) |>
      dplyr::select("ncbi_refseq", "uniprot_id")

    uniparc_prepared <- uniparc_search_results |>
      dplyr::select(ncbi_refseq, uniparc_id = uniprot_id)

    aa_seq_tbl_updated <- aa_seq_tbl |>
      dplyr::left_join(uniprot_filtered, by = "ncbi_refseq") |>
      dplyr::left_join(uniparc_prepared, by = "ncbi_refseq") |>
      dplyr::mutate(database_id = dplyr::coalesce(uniprot_id, uniparc_id)) |>
      dplyr::select(-uniprot_id, -uniparc_id)

    return(aa_seq_tbl_updated)
  }

  message("Reading FASTA file...")
  flush.console()

  suppressMessages({
    fasta_file_raw <- vroom::vroom(fasta_file_path, delim = "\n", col_names = FALSE, progress = FALSE)
  })
  first_line <- fasta_file_raw$X1[1]

  if (startsWith(first_line, ">sp|") || startsWith(first_line, ">tr|")) {
    message("Processing standard UniProt FASTA format...")
    flush.console()
    aa_seq_tbl <- parseFastaFileStandard(fasta_file_path)
    
    # Create metadata for standard UniProt format
    fasta_metadata <- list(
      fasta_format = "standard_uniprot",
      available_columns = colnames(aa_seq_tbl),
      has_protein_evidence = "protein_evidence" %in% colnames(aa_seq_tbl),
      has_gene_names = "gene_name" %in% colnames(aa_seq_tbl),
      has_isoform_info = "is_isoform" %in% colnames(aa_seq_tbl),
      has_status_info = "status" %in% colnames(aa_seq_tbl),
      num_sequences = nrow(aa_seq_tbl),
      processing_timestamp = Sys.time()
    )
    
    message("Saving results...")
    flush.console()
    saveRDS(aa_seq_tbl, fasta_meta_file)
    
    # Return list with data and metadata
    return(list(
      aa_seq_tbl_final = aa_seq_tbl,
      fasta_metadata = fasta_metadata
    ))
  } else {
    message("Processing non-standard FASTA format...")
    flush.console()
    aa_seq_tbl <- parseFastaFileNonStandard(fasta_file_path)

    if (!is.null(uniprot_search_results) && !is.null(uniparc_search_results)) {
      aa_seq_tbl_final <- matchAndUpdateDataFrames(aa_seq_tbl, uniprot_search_results, uniparc_search_results, organism_name)
    } else {
      aa_seq_tbl_final <- aa_seq_tbl |>
        dplyr::mutate(database_id = NA_character_)
    }
    
    # Create metadata for non-standard format
    fasta_metadata <- list(
      fasta_format = "non_standard",
      available_columns = colnames(aa_seq_tbl_final),
      has_protein_evidence = "protein_evidence" %in% colnames(aa_seq_tbl_final),
      has_gene_names = "gene_name" %in% colnames(aa_seq_tbl_final),
      has_isoform_info = "is_isoform" %in% colnames(aa_seq_tbl_final),
      has_status_info = "status" %in% colnames(aa_seq_tbl_final),
      num_sequences = nrow(aa_seq_tbl_final),
      processing_timestamp = Sys.time()
    )

    message("Writing results...")
    flush.console()

    vroom::vroom_write(aa_seq_tbl_final,
                       file = "aa_seq_tbl.tsv",
                       delim = "\t",
                       na = "",
                       quote = "none",
                       progress = FALSE)

    saveRDS(aa_seq_tbl_final, fasta_meta_file)
    
    # Return list with data and metadata
    return(list(
      aa_seq_tbl_final = aa_seq_tbl_final,
      fasta_metadata = fasta_metadata
    ))
  }
}


# ----------------------------------------------------------------------------
# updateProteinIDs
# ----------------------------------------------------------------------------
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

