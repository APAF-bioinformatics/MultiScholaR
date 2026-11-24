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
#' @param id_string Runtime inputs used by this function; see the usage section for accepted values.
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
#' @details
#' When the retrieval fails, a warning of class `prot_annotation_uniprot_empty` with its `reason` (`request_failed`,
#' `processing_failed`, or `id_conversion_failed` when the Ensembl conversion failed) is signalled and an empty
#' annotation table is returned. Requests carry a time-out and up to three attempts (see `directUniprotDownload`).
#' @export
#' @param chunk_size,progress_callback Runtime inputs used by this function; see the usage section for accepted values.
getUniprotAnnotationsFull <- function(data_tbl,
                                     protein_id_column,
                                     cache_dir,
                                     taxon_id = 9606,
                                     chunk_size = 25,
                                     progress_callback = NULL) {
  
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
  
  
  # === ENSEMBL ID CONVERSION LOGIC ===
  # Detect if protein IDs are Ensembl format and convert to UniProt if necessary
  ensembl_detection <- detectEnsemblIds(cleaned_proteins)
  ensembl_mapping <- NULL  # Initialize mapping variable
  proteins_to_query <- cleaned_proteins  # Default: use original proteins
  conversion_failed <- NULL  # KIT-175: the reason the Ensembl conversion failed, when it did
  
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
        cat("=== END ENSEMBL CONVERSION ===\n\n")
      } else {
        cat("[FAIL] ERROR: gprofiler2 conversion returned 0 results!\n")
        warning("No Ensembl IDs could be converted. Falling back to direct query with original IDs.")
        conversion_failed <- "no identifier could be converted"
        proteins_to_query <- cleaned_proteins
        ensembl_mapping <- NULL
      }
      
    }, error = function(e) {
      warning(sprintf("Ensembl conversion failed: %s. Falling back to direct UniProt query.", e$message))
      conversion_failed <<- uniprotConditionDetail(e)
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
  
  tryCatch({
    # KIT-175: when the Ensembl conversion failed, the query found no UniProt accession; that outcome is recorded and
    # signalled as `id_conversion_failed` instead (the only condition muffled on the way).
    uniprot_annotations <- withCallingHandlers(
      getUniprotAnnotations(
        input_tbl = temp_protein_table,
        cache_dir = cache_dir,
        taxon_id = taxon_id,
        progress_callback = progress_callback
      ),
      prot_annotation_uniprot_empty = function(w) {
        if (!is.null(conversion_failed) && identical(w$reason, "no_uniprot_accessions")) {
          recordUniprotOutcome(NA_character_, "not_queried",
            "not queried: the Ensembl identifiers could not be converted to UniProt accessions")
          signalUniprotEmpty("id_conversion_failed", w$requested, detail = conversion_failed, taxon_id = taxon_id)
          invokeRestart("muffleWarning")
        }
      }
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
        go_ids_lost <- attr(uniprot_annotations, "go_ids_not_in_go_db", exact = TRUE)
        uniprot_annotations <- merge(
          uniprot_annotations,
          ensembl_lookup,
          by = "Entry",
          all.x = TRUE
        )
        # KIT-173: merge() drops the table's attributes; the lost-GO-id count is put back.
        attr(uniprot_annotations, "go_ids_not_in_go_db") <- go_ids_lost
        
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
    # KIT-175: an error already recorded (a transport error, a processing failure) is signalled with its reason; any
    # other error here (the Ensembl merge-back) is recorded as `processing_failed` and signalled. The empty annotation
    # table getUniprotAnnotations returns when its download fails (KIT-101 F-28) is returned; it is not cached.
    if (!inherits(e, "prot_annotation_uniprot_recorded")) {
      recordUniprotOutcome(NA_character_, "failed", "annotation failed after retrieval: the answer could not be processed")
      e <- markUniprotRecorded(e, "processing_failed", list(detail = uniprotConditionDetail(e),
        requested = length(proteins_to_query), taxon_id = taxon_id))
    }
    message(sprintf("UniProt annotation ended without a table: %s", e$reason))
    counts <- unclass(e)[intersect(c("requests_total", "requests_failed", "requests_not_sent", "requests_truncated",
      "accessions_requested", "accessions_not_looked_up", "hits_beyond_rows", "statuses", "stage"), names(e))]
    signalUniprotEmpty(e$reason, e$requested %||% length(proteins_to_query), detail = e$detail %||% NA_character_,
      taxon_id = taxon_id, counts = counts)
    return(createEmptyUniprotTable())
  })
}

#' The design step's UniProt annotation (KIT-103 W3), without its GUI.
#' @description The cache directory, the retrieval and the scripts-directory copy that
#'   `completeProtDesignPostCheckpoint` ran inline; its modal, progress messages and notifications stay there.
#' @param manager Current owner (the workflow's state manager).
#' @param dataCln Accepted assay table reference.
#' @param proteinColumn Declared protein identifier column; NULL is refused with today's message.
#' @param taxonId Taxonomic ID, passed through unchanged.
#' @param resultsDir Experiment results directory, or NULL for the session's temporary cache.
#' @param sourceDir Experiment scripts directory, or NULL to skip the `uniprot_dat_cln.RDS` copy.
#' @param progressCallback Optional `function(current, total)`.
#' @return The UniProt annotation table.
#' @noRd
retrieveProtDesignUniprotAnnotations <- function(manager, dataCln, proteinColumn, taxonId,
    resultsDir = NULL, sourceDir = NULL, progressCallback = NULL) {
    # KIT-175: an error before any request is recorded as `not_attempted` and re-raised, message unchanged.
    notAttempted <- function(e, detail) {
        recordUniprotOutcome(NA_character_, "not_queried", sprintf("not queried: annotation was not attempted (%s)", detail))
        stop(markUniprotRecorded(e, "not_attempted", list(detail = detail)))
    }
    if (is.null(proteinColumn)) {
        notAttempted(simpleError("Protein column not found in column_mapping. Cannot get annotations."),
            "no protein identifier column is selected")
    }

    cacheDir <- if (!is.null(resultsDir)) {
        file.path(resultsDir, "cache")
    } else {
        file.path(tempdir(), "proteomics_cache")
    }

    uniprotCacheDir <- file.path(cacheDir, "uniprot_annotations")
    if (!dir.exists(uniprotCacheDir)) {
        dir.create(uniprotCacheDir, recursive = TRUE)
    }

    # KIT-175: the classed conditions are kept, never muffled (they reach every caller); the last one becomes the
    # table's `annotation_outcome` attribute, which travels with it to the bus, the project file and the session.
    rowsBefore <- uniprotRecordRows()
    lastCondition <- NULL
    keep <- function(w) lastCondition <<- w
    uniprotDatCln <- tryCatch(
        withCallingHandlers(
            getUniprotAnnotationsFull(
                data_tbl = readDIAAnnotationIdentifiers(manager, dataCln, proteinColumn),
                protein_id_column = proteinColumn,
                cache_dir = uniprotCacheDir,
                taxon_id = taxonId,
                progress_callback = progressCallback
            ),
            prot_annotation_uniprot_empty = keep,
            prot_annotation_uniprot_incomplete = keep
        ),
        error = function(e) {
            if (inherits(e, "prot_annotation_uniprot_recorded")) stop(e)
            detail <- if (inherits(e, "saved_entity_limit")) {
                "the data have more than 10,000 distinct protein groups, the most the annotation step reads"
            } else if (identical(conditionMessage(e), "data_tbl cannot be empty")) {
                "the data have no protein identifiers"
            } else {
                uniprotConditionDetail(e)
            }
            notAttempted(e, detail)
        }
    )
    if (!is.null(lastCondition)) {
        attr(uniprotDatCln, "annotation_outcome") <- uniprotOutcomeFromCondition(lastCondition)
    }

    # KIT-115 (A101-09): an empty table is a failed retrieval (F-28) and is not written as the project's annotations.
    if (!is.null(sourceDir) && nrow(uniprotDatCln) > 0) {
        scriptsUniprotPath <- file.path(sourceDir, "uniprot_dat_cln.RDS")
        tryCatch({
            saveRDS(uniprotDatCln, scriptsUniprotPath)
            log_info(sprintf("Saved uniprot_dat_cln to scripts directory: %s", scriptsUniprotPath))
        }, error = function(e) {
            # KIT-175 (C04): recorded only after a retrieval that wrote a row in this call (none on a cache hit).
            detail <- uniprotConditionDetail(e)
            if (uniprotRecordRows() > rowsBefore) {
                recordUniprotOutcome(NA_character_, "failed",
                    sprintf("annotation retrieved but not saved to the project file (%s)", detail))
            }
            stop(markUniprotRecorded(e, "processing_failed", list(detail = detail, stage = "save")))
        })
    }

    uniprotDatCln
}
