# ----------------------------------------------------------------------------
# getUniprotRegexPatterns
# ----------------------------------------------------------------------------
#' @title Get UniProt Regex Patterns
#' @description Returns a named list of regex patterns used for UniProt ID management.
#' @return A list with patterns: isoform, decoy, contaminant, delimiter
#' @export
getUniprotRegexPatterns <- function() {
  list(
    isoform = "-\\d+$",
    decoy = "REV__",
    contaminant = "CON__",
    delimiter = " |;|:|\\|",
    entry = "^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-\\d+)?$"
  )
}

# ----------------------------------------------------------------------------
# normalizeUniprotAccession
# ----------------------------------------------------------------------------
#' @title Normalize UniProt Accession
#' @description Cleans UniProt accessions by removing isoforms, decoys, and contaminants.
#' @param string Character vector of protein accessions.
#' @param remove_isoform Logical, whether to remove isoform suffixes (default: TRUE).
#' @param remove_decoy Logical, whether to remove decoy prefixes (default: FALSE).
#' @param remove_contaminant Logical, whether to remove contaminant prefixes (default: FALSE).
#' @return Character vector of normalized accessions.
#' @export
normalizeUniprotAccession <- function(string, 
                                     remove_isoform = TRUE, 
                                     remove_decoy = FALSE, 
                                     remove_contaminant = FALSE) {
  if (is.null(string) || length(string) == 0) return(string)
  
  patterns <- getUniprotRegexPatterns()
  res <- string
  
  if (remove_decoy) {
    res <- str_replace(res, patterns$decoy, "")
  }
  
  if (remove_contaminant) {
    res <- str_replace(res, patterns$contaminant, "")
  }
  
  if (remove_isoform) {
    res <- str_replace(res, patterns$isoform, "")
  }
  
  return(res)
}

# ----------------------------------------------------------------------------
# cleanIsoformNumber
# ----------------------------------------------------------------------------
#' @title Clean Isoform Number
#' @description Remove isoform numbers from protein IDs.
#' @param string Character vector of protein IDs.
#' @return Cleaned protein IDs with isoform numbers removed.
#' @export
cleanIsoformNumber <- function(string) {
  # Wrapper for the centralized normalization function
  normalizeUniprotAccession(string, remove_isoform = TRUE)
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

