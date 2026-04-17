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
#'   \item Handles protein groups properly (A;B;C -> A, B, C)
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
# ----------------------------------------------------------------------------
# getUniprotAnnotationsFull
# ----------------------------------------------------------------------------
#' @title Get Full UniProt Annotations
#' @description Retrieves comprehensive protein information from UniProt, including Cross-references, 
#' Gene Ontology, and functional data.
#' @param data_tbl Data frame containing protein IDs.
#' @param protein_id_column Name of the column containing protein IDs.
#' @param cache_dir Directory of the UniProt cache.
#' @param taxon_id Taxonomic ID of the organism.
#' @return Data frame with UniProt annotations.
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
  # [OK] REFACTORED: Use centralized UniProt normalization
  cleaned_proteins <- unique(normalizeUniprotAccession(cleaned_ids, remove_isoform = TRUE))
  
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
      
      # Convert Ensembl -> UniProt
      conversion_result <- convertEnsemblToUniprot(cleaned_proteins, organism_code)
      
      # Filter to successfully converted IDs for UniProt query
      successfully_converted <- conversion_result[
        conversion_result$conversion_status == "success" & 
        !is.na(conversion_result$converted_id), 
      ]
      
      if (nrow(successfully_converted) > 0) {
        # [OK] CRITICAL: Extract UNIPROT IDs (not Ensembl IDs!) from conversion result
        proteins_to_query <- successfully_converted$converted_id
        ensembl_mapping <- conversion_result  # Save for later merge
        
        cat(sprintf("\n[OK] Conversion summary:\n"))
        cat(sprintf("  - Successfully converted: %d/%d (%.1f%%)\n", 
                    nrow(successfully_converted), 
                    nrow(conversion_result),
                    100 * nrow(successfully_converted) / nrow(conversion_result)))
        cat(sprintf("  - Failed conversions: %d (will have NA annotations)\n", 
                    nrow(conversion_result) - nrow(successfully_converted)))
        cat(sprintf("  - UniProt IDs extracted for API query: %d\n", length(proteins_to_query)))
        cat(sprintf("\n[OK] CRITICAL CHECK - First 5 IDs to be submitted to UniProt API:\n"))
        cat(sprintf("  Type: %s\n", if(any(grepl("^ENS", head(proteins_to_query, 5)))) "[FAIL] ENSEMBL (WRONG!)" else "[OK] UniProt (CORRECT)"))
        print(head(proteins_to_query, 5))
        cat("=== END ENSEMBL CONVERSION ===\n\n")
      } else {
        cat("[FAIL] ERROR: gprofiler2 conversion returned 0 results!\n")
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
                "[FAIL] STILL ENSEMBL IDS - BUG IN CONVERSION LOGIC!"
              } else {
                "[OK] UniProt IDs - Ready for API"
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

