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
  
  # [OK] FIX: Strip version suffixes (.1, .2, etc.) before conversion
  # gprofiler2 doesn't recognize Ensembl IDs with version suffixes
  cat("Stripping version suffixes from Ensembl IDs (e.g., .1, .2)...\n")
  
  # Create mapping: original_with_version -> stripped_version
  version_mapping <- data.frame(
    original_with_version = ensembl_ids,
    stripped_version = gsub("\\.\\d+$", "", ensembl_ids),  # Remove .X at end
    stringsAsFactors = FALSE
  )
  
  # Get unique stripped IDs (multiple versioned IDs may map to same stripped ID)
  unique_stripped <- unique(version_mapping$stripped_version)
  
  cat(sprintf("Stripped versions: %d original IDs -> %d unique stripped IDs\n", 
              length(ensembl_ids), length(unique_stripped)))
  cat(sprintf("Example: %s -> %s\n", 
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
    
    cat(sprintf("[OK] gprofiler2 returned %d conversion results\n", 
                if(is.null(conversion_result)) 0 else nrow(conversion_result)))
    
    if (!is.null(conversion_result) && nrow(conversion_result) > 0) {
      cat("First 3 conversions from gprofiler2:\n")
      print(head(conversion_result[, c("input", "target", "name")], 3))
    }
    
    # Create mapping table for STRIPPED -> UniProt
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
    
    # [OK] FIX: Map back to ORIGINAL IDs (with versions)
    # Join: original_with_version -> stripped_version -> converted_id
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
    cat(sprintf("[OK] Successfully converted %d/%d Ensembl IDs to UniProt IDs (%.1f%%)\n", 
                success_count, length(ensembl_ids), 
                100 * success_count / length(ensembl_ids)))
    
    if (success_count < length(ensembl_ids)) {
      failed_count <- length(ensembl_ids) - success_count
      cat(sprintf("[WARNING] Warning: %d Ensembl IDs could not be converted (will use original ID with NA annotations)\n",
                  failed_count))
      if (failed_count <= 10) {
        cat("Failed IDs:\n")
        print(final_mapping$original_id[final_mapping$conversion_status != "success"])
      }
    }
    
    return(final_mapping)
    
  }, error = function(e) {
    cat(sprintf("[FAIL] ERROR in gprofiler2::gconvert: %s\n", e$message))
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

