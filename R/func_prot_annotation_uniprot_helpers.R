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
#' @details
#' A cache file with zero rows, or one that cannot be read, is not used: it is removed and the query runs. An empty
#' or incomplete answer (the conditions of `directUniprotDownload`, which pass through) is never written to the
#' cache. An error in the processing of the answer or in the cache write is recorded once and re-raised with the
#' class `prot_annotation_uniprot_recorded` added (field `reason = "processing_failed"`).
#'
#' @return A data frame containing UniProt annotations and GO terms. A table built here carries the attribute
#'   `go_ids_not_in_go_db`: the GO ids of the download that the installed GO.db does not hold (and that the table
#'   therefore omits), counted by class (obsolete, secondary, unknown), with up to 20 examples and the GO.db release.
#'
#' @export
#' @param progress_callback Runtime inputs used by this function; see the usage section for accepted values.
#' @param label NULL (the project's own annotation), or a text appended to the `query_summary` of every analysis-record
#'   row the retrieval writes (KIT-167: `"reference annotation: <reference> (<UPID> or supplied)"`).
getUniprotAnnotations <- function(input_tbl, 
                                 cache_dir, 
                                 taxon_id,
                                 force_download = FALSE,
                                 batch_size = 25,
                                 timeout = 600,
                                 api_delay = 1,
                                 progress_callback = NULL,
                                 label = NULL) {
  
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
  
  
  # Check if cache exists and should be used. KIT-175: a zero-row or unreadable cache file is a miss and is removed
  # (a successful retrieval writes it again), so an empty answer cached by an earlier version is not returned.
  if (!force_download && file.exists(cache_file)) {
    cached <- tryCatch(readRDS(cache_file), error = function(e) NULL)
    if (is.data.frame(cached) && nrow(cached) > 0) {
      message("   DEBUG66: Cache hit, loading from cache")
      message("Loading cached UniProt annotations...")
      return(cached)
    }
    message("   The UniProt cache file is empty or unreadable; it is ignored and removed.")
    unlink(cache_file)
  }
  
  # Download annotations if needed
  message("   DEBUG66: Cache miss, calling directUniprotDownload...")
  message("Fetching UniProt annotations...")
  # KIT-175: an empty answer is signalled by directUniprotDownload (class `prot_annotation_uniprot_empty`) and an
  # incomplete one (class `prot_annotation_uniprot_incomplete`) is seen here without being muffled; neither is cached.
  # An error in the processing or the cache write is recorded once and re-raised with its reason.
  incomplete <- FALSE
  uniprot_dat_cln <- tryCatch({
    annotations <- withCallingHandlers(
      directUniprotDownload(
        input_tbl = input_tbl,
        output_path = raw_results_file,
        taxon_id = taxon_id,
        batch_size = batch_size,
        timeout = timeout,
        api_delay = api_delay,
        progress_callback = progress_callback,
        label = label
      ),
      prot_annotation_uniprot_incomplete = function(w) incomplete <<- TRUE
    )

    # Process annotations or create empty table if download failed
    if (!is.null(annotations) && nrow(annotations) > 0) {
      message("Processing GO terms...")
      # KIT-173: the raw GO ids the installed GO.db does not hold are dropped by the conversion below; they are counted
      # here, by class, with the GO.db release in force, and the count travels with the table.
      go_ids_lost <- goIdsNotInGoDb(unlist(strsplit(as.character(annotations$Gene.Ontology.IDs), "; ", fixed = TRUE)))
      processed_annotations <- annotations |>
        uniprotGoIdToTerm(
          uniprot_id_column = Entry,
          go_id_column = Gene.Ontology.IDs,
          sep = "; "
        )

      # Standardize column names
      processed <- standardizeUniprotColumns(processed_annotations)
      attr(processed, "go_ids_not_in_go_db") <- go_ids_lost

      if (incomplete) {
        message("UniProt annotations are incomplete; they are not saved to the cache.")
      } else {
        # Save to cache
        saveRDS(processed, cache_file)
        message("UniProt annotations saved to cache.")
      }
      processed
    } else {
      createEmptyUniprotTable()
    }
  }, error = function(e) {
    if (!inherits(e, "prot_annotation_uniprot_recorded")) {
      recordUniprotOutcome(NA_character_, "failed",
        uniprotLabelSummary("annotation failed after retrieval: the answer could not be processed", label))
      e <- markUniprotRecorded(e, "processing_failed", list(detail = uniprotConditionDetail(e),
        requested = length(unique(input_tbl$Protein.Ids)), taxon_id = taxon_id))
    }
    stop(e)
  })
  
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
#' @details
#' The query terms are field-qualified (`accession:<acc>`; D8, VAL-92). Each request (one chunk) has a time-out of 30
#' seconds per attempt and up to three attempts; a transport failure and the statuses 429, 500, 502, 503 and 504 are
#' tried again (`Retry-After` honoured up to 10 seconds). After three consecutive requests refused with such a status
#' the remaining requests are not sent. A request fails when its final status is not 200 or its answer is not a
#' UniProt table; an answer is cut short when a requested accession (isoform suffixes are removed before the query;
#' the counts are of distinct canonical accessions) is absent and the answer announces a next page
#' (`Link: rel="next"`); its rows are kept and its absent accessions count as not looked up. `X-Total-Results` is a
#' count only. When no request returns rows a warning of class `prot_annotation_uniprot_empty` is signalled (field
#' `reason`: `service_status`, `invalid_response`, `no_rows` or `no_uniprot_accessions`) and NULL is returned; when
#' some requests fail or are cut short and others return rows, a warning of class
#' `prot_annotation_uniprot_incomplete` with the request counts is signalled and the rows are returned. A transport
#' failure is re-raised with the class `prot_annotation_uniprot_recorded` added (field `reason = "request_failed"`).
#' Under a workflow owner the outcome is written once to the analysis record.
#'
#' @return A data frame containing the raw UniProt results
#'
#' @export
#' @param progress_callback Runtime inputs used by this function; see the usage section for accepted values.
#' @param label NULL, or a text appended to the `query_summary` of every analysis-record row written (KIT-167).
directUniprotDownload <- function(input_tbl, 
                                 output_path, 
                                 taxon_id, 
                                 batch_size = 25,
                                 timeout = 600, 
                                 api_delay = 1,
                                 progress_callback = NULL,
                                 label = NULL) {
  # KIT-167: a reference annotation's record rows name the reference (default NULL keeps the text).
  recordUniprotOutcome <- function(release, release_kind, summary) {
    get("recordUniprotOutcome", envir = topenv())(release, release_kind, uniprotLabelSummary(summary, label))
  }
  
  message(sprintf("   Input: %d unique protein IDs", nrow(input_tbl)))
  message(sprintf("   Taxon ID: %d", taxon_id))
  message(sprintf("   Batch size: %d", batch_size))
  message(sprintf("   Output path: %s", output_path))
  
  # Set a longer timeout
  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout))
  options(timeout = timeout)
  
  # Extract unique protein IDs
  original_protein_ids <- unique(input_tbl$Protein.Ids)
  message(paste("Found", length(original_protein_ids), "unique protein IDs from input"))
  
  # Official UniProt accession regex (with optional isoform suffix)
  uniprot_regex <- getUniprotRegexPatterns()$entry
  
  # Filter IDs matching the pattern
  protein_ids <- original_protein_ids[grepl(uniprot_regex, original_protein_ids)]
  
  n_skipped <- length(original_protein_ids) - length(protein_ids)
  if (n_skipped > 0) {
    message(sprintf("   Filtered out %d IDs that do not match UniProt accession pattern", n_skipped))
    skipped_ids <- setdiff(original_protein_ids, protein_ids)
    message(sprintf("   Example skipped IDs: %s", paste(head(skipped_ids, 5), collapse = ", ")))
  }
  # KIT-175 (I01): an isoform is asked for by its canonical accession. A field-qualified isoform term
  # (`accession:P05067-4`) makes UniProt return isoform entries for every accession of the request; the bare query
  # of earlier versions returned the canonical entry. The design step already strips isoforms.
  protein_ids <- unique(normalizeUniprotAccession(protein_ids, remove_isoform = TRUE))
  
  if (length(protein_ids) == 0) {
    message("   No valid UniProt IDs found to query. Skipping download.")
    # KIT-175: no request is made; the outcome is recorded and signalled.
    requested <- length(original_protein_ids)
    recordUniprotOutcome(NA_character_, "not_queried", if (requested == 0L) {
      "not queried: this step had no protein identifiers (annotation is retrieved at the design step)"
    } else {
      sprintf("not queried: none of %d identifiers is a UniProt accession", requested)
    })
    signalUniprotEmpty("no_uniprot_accessions", requested, taxon_id = taxon_id)
    return(NULL)
  }
  
  message(paste("Querying", length(protein_ids), "valid UniProt IDs"))
  
  # Split into batches
  chunks <- split(protein_ids, ceiling(seq_along(protein_ids)/batch_size))
  message(paste("Split into", length(chunks), "chunks for processing"))
  
  # Function to process one chunk. KIT-175: each request is one `req_perform()` of up to three attempts
  # (`buildUniprotSearchRequest`); its outcome is returned as a list and the retrieval's outcome is decided once all
  # requests are made. After three consecutive requests that each ended with a retried status, the remaining requests
  # are not sent (outage stop).
  consecutive_refused <- 0L
  process_chunk <- function(chunk, chunk_idx, total_chunks) {
    # Convert chunk_idx to integer (purrr::imap passes names as characters)
    chunk_idx <- as.integer(chunk_idx)
    outcome <- list(kind = "failed", size = length(chunk), status = NA_integer_, truncated = FALSE,
      not_looked_up = length(chunk), hits_beyond_rows = NA_integer_, release = NA_character_, rows = NULL)
    if (consecutive_refused >= 3L) {
      outcome$kind <- "not_sent"
      return(outcome)
    }
    
    message(paste("Processing chunk", chunk_idx, "of", total_chunks, "with", length(chunk), "IDs"))
    message("   Chunk protein IDs:")
    print(chunk)
    
    # Create query for this batch. KIT-175 (D8, VAL-92): field-qualified terms, so that a term matches that accession
    # only and a request of at most 25 terms has at most 25 hits.
    query <- paste0("(", paste0("accession:", chunk, collapse = " OR "), ") AND organism_id:", taxon_id)
    
    response <- httr2::req_perform(buildUniprotSearchRequest(query))
    
    # Be nice to the API
    Sys.sleep(api_delay)
    
    # Update progress every 5 chunks
    if (!is.null(progress_callback) && chunk_idx %% 5 == 0) {
      progress_callback(chunk_idx, total_chunks)
    }
    
    # Check if successful
    outcome$status <- as.integer(httr2::resp_status(response))
    if (outcome$status != 200L) {
      message(paste("  Request failed with status", outcome$status))
      consecutive_refused <<- if (outcome$status %in% UNIPROT_TRANSIENT_STATUSES) consecutive_refused + 1L else 0L
      return(outcome)
    }
    consecutive_refused <<- 0L
    # KIT-026: the release header of a successful response, when the service sends one.
    outcome$release <- httr2::resp_header(response, "x-uniprot-release") %||% NA_character_
    chunk_result <- tryCatch({
      content <- httr2::resp_body_string(response, encoding = "UTF-8")
      
      temp_file <- tempfile(fileext = ".tsv")
      writeLines(content, temp_file)
      
      suppressWarnings(
        read.delim(temp_file, sep="\t", quote="", stringsAsFactors=FALSE)
      )
    }, error = function(e) {
      message(paste("      Error message:", e$message))
      message(paste("      Error class:", paste(class(e), collapse = ", ")))
      message(paste("      Error occurred at chunk", chunk_idx))
      NULL
    })
    if (is.null(chunk_result) || !"Entry" %in% names(chunk_result)) {
      message("  API returned 200 but the answer is not a UniProt table")
      return(outcome)
    }
    # KIT-175 (item 1a): `X-Total-Results` is a count only. An answer is cut short when a requested accession (isoform
    # suffix removed: UniProt returns the canonical entry) is absent and the answer announces a next page.
    total <- suppressWarnings(as.integer(httr2::resp_header(response, "x-total-results") %||% NA_character_))
    if (!is.na(total)) outcome$hits_beyond_rows <- max(total - nrow(chunk_result), 0L)
    absent <- !sub("-[0-9]+$", "", chunk) %in% as.character(chunk_result$Entry)
    link <- httr2::resp_header(response, "link") %||% ""
    outcome$truncated <- any(absent) && grepl('rel="next"', link, fixed = TRUE)
    outcome$not_looked_up <- if (outcome$truncated) sum(absent) else 0L
    if (outcome$truncated) {
      message(sprintf("  API returned %d rows and a next page with %d requested IDs missing", nrow(chunk_result),
        sum(absent)))
    }
    if (nrow(chunk_result) > 0) {
      message(paste("  Found", nrow(chunk_result), "results"))
      outcome$kind <- "rows"
      outcome$rows <- chunk_result
    } else {
      message("  API returned 200 but result has 0 rows")
      outcome$kind <- "empty"
    }
    outcome
  }
  
  # Process all chunks using imap (provides both value and index)
  total_chunks <- length(chunks)
  accessions <- length(protein_ids)
  
  message("   First chunk protein IDs:")
  
  # Wrap purrr::imap in tryCatch to catch errors with full context. KIT-175: a transport error (after the request's
  # attempts) ends the retrieval; it is recorded once and re-raised with the class `prot_annotation_uniprot_recorded`.
  results <- tryCatch({
    purrr::imap(chunks, ~ process_chunk(.x, .y, total_chunks))
  }, error = function(e) {
    message(sprintf("   Error message: %s", e$message))
    message(sprintf("   Error class: %s", paste(class(e), collapse = ", ")))
    if (isUniprotTransportError(e)) {
      recordUniprotOutcome(NA_character_, "failed",
        sprintf("annotation failed: UniProt could not be reached; 0 of %d accessions annotated", accessions))
      e <- markUniprotRecorded(e, "request_failed", list(detail = uniprotConditionDetail(e), requested = accessions,
        requests_total = total_chunks, taxon_id = taxon_id))
    }
    stop(e)
  })
  
  # Send final progress update to show 100% completion
  if (!is.null(progress_callback)) {
    progress_callback(total_chunks, total_chunks)
  }
  
  # KIT-175: the requests counted by outcome (items 1, 1a); a request not sent counts as failed.
  kinds <- vapply(results, function(outcome) outcome$kind, character(1))
  statuses <- vapply(results, function(outcome) as.integer(outcome$status), integer(1))
  truncated <- vapply(results, function(outcome) isTRUE(outcome$truncated), logical(1))
  not_looked_up <- vapply(results, function(outcome) as.integer(outcome$not_looked_up), integer(1))
  hits <- vapply(results, function(outcome) as.integer(outcome$hits_beyond_rows), integer(1))
  failed <- kinds %in% c("failed", "not_sent")
  refused <- failed & !is.na(statuses) & statuses != 200L
  releases <- stats::na.omit(vapply(results, function(outcome) as.character(outcome$release), character(1)))
  release <- if (length(releases)) releases[[1L]] else NA_character_
  counts <- list(
    requests_total = length(results),
    requests_failed = sum(failed),
    requests_not_sent = sum(kinds == "not_sent"),
    requests_truncated = sum(truncated & !failed),
    accessions_requested = accessions,
    accessions_not_looked_up = sum(not_looked_up[failed | truncated]),
    hits_beyond_rows = sum(hits, na.rm = TRUE),
    statuses = sort(unique(statuses[refused]))
  )
  
  if (!any(kinds == "rows")) {
    reason <- if (any(refused)) "service_status" else if (any(failed | truncated)) "invalid_response" else "no_rows"
    message("Failed to retrieve any data from UniProt")
    if (identical(reason, "no_rows")) {
      recordUniprotOutcome(release, "header",
        sprintf("no entry found for %d accessions under organism %s", accessions, taxon_id))
    } else {
      recordUniprotOutcome(NA_character_, "failed", uniprotRecordSummary(reason, counts))
    }
    signalUniprotEmpty(reason, accessions, taxon_id = taxon_id, counts = counts)
    return(NULL)
  }
  
  # Combine results
  all_results <- tryCatch({
    combined <- purrr::reduce(lapply(results[kinds == "rows"], function(outcome) outcome$rows), rbind)
    
    # Standardize column names for downstream processing
    names(combined) <- gsub(" ", ".", names(combined))
    
    # Add From column needed for downstream processing
    combined$From <- combined$Entry
    combined
  }, error = function(e) {
    message(paste("Error combining results from chunks:", e$message))
    message(paste("Number of successful chunks:", sum(kinds == "rows")))
    NULL
  })
  if (is.null(all_results)) {
    recordUniprotOutcome(NA_character_, "failed", uniprotRecordSummary("invalid_response", counts))
    signalUniprotEmpty("invalid_response", accessions, taxon_id = taxon_id, counts = counts)
    return(NULL)
  }
  
  # Write to file (KIT-175: outside the combine handler, so that its failure is a processing failure)
  write.table(all_results, output_path, sep="\t", quote=FALSE, row.names=FALSE)
  message(paste("Successfully retrieved", nrow(all_results), "entries from UniProt"))
  
  if (counts$requests_failed > 0 || counts$requests_truncated > 0) {
    recordUniprotOutcome(release, "header", uniprotRecordSummary("incomplete", counts))
    signalUniprotIncomplete(accessions, taxon_id = taxon_id, counts = counts)
  } else {
    recordUniprotOutcome(release, "header", sprintf("uniprotkb search, organism_id %s", taxon_id))
  }
  all_results
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
  requireOptionalPackage("UniProt.ws", "UniProt annotation") # KIT-064 (roar BUG-24)
  requireOptionalPackage("AnnotationDbi", "UniProt annotation") # KIT-064 (roar BUG-24)
  requireOptionalPackage("GO.db", "UniProt annotation") # KIT-064 (roar BUG-24)


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

    up <- UniProt.ws::UniProt.ws(taxId = taxonomy_id)
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
    up_cls <- unlist(AnnotationDbi::columns(up))
    list_intersect<-intersect(list_of_sp_columns,up_cls)
    if(length(setdiff( list_of_sp_columns,list_intersect)) > 0)
    {
      print(paste("UniProt fields not found:", setdiff( list_of_sp_columns,list_intersect),sep=", "))
    }

    my_keytype <- "UniProtKB"
    if ("UNIPROTKB" %in% AnnotationDbi::keytypes(up)) {
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
    goterms <- AnnotationDbi::Term(GO.db::GOTERM)
    gotypes <- AnnotationDbi::Ontology(GO.db::GOTERM)


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
# KIT-175: the UniProt request, and the outcome of a retrieval that ends empty or incomplete
# ----------------------------------------------------------------------------

#' The response statuses after which a UniProt request is tried again (KIT-175, D5).
#' @noRd
UNIPROT_TRANSIENT_STATUSES <- c(429L, 500L, 502L, 503L, 504L)

#' The UniProt search request for one chunk (KIT-175, D5, VAL-30)
#'
#' A time-out of 30 s per attempt, up to three attempts (a transport failure or one of the statuses in
#' `UNIPROT_TRANSIENT_STATUSES` is tried again; `Retry-After` is honoured up to 10 s, otherwise httr2's backoff),
#' the user agent `MultiScholaR`, and no HTTP error raised: the caller reads the status.
#' @param query The UniProt query string.
#' @return An httr2 request.
#' @noRd
buildUniprotSearchRequest <- function(query) {
  httr2::request("https://rest.uniprot.org/uniprotkb/search") |>
    httr2::req_url_query(
      query = query,
      format = "tsv",
      fields = "accession,id,protein_name,gene_names,organism_name,length,go_id,reviewed,protein_existence,annotation_score"
    ) |>
    httr2::req_user_agent("MultiScholaR") |>
    httr2::req_timeout(30) |>
    httr2::req_retry(
      max_tries = 3L,
      retry_on_failure = TRUE,
      is_transient = function(resp) httr2::resp_status(resp) %in% UNIPROT_TRANSIENT_STATUSES,
      after = function(resp) {
        wait <- httr2::resp_retry_after(resp)
        if (is.na(wait)) NA else min(wait, 10)
      }
    ) |>
    httr2::req_error(is_error = function(response) FALSE)
}

#' A condition's text as it may be shown or recorded: cut at the first `?` (a URL's query holds accessions) and with
#' every token matching the UniProt accession pattern replaced by `<accession>` (licence guard, ENG-07).
#' @param text Character.
#' @return One string.
#' @noRd
redactUniprotDetail <- function(text) {
  text <- sub("\\?.*$", "", paste(as.character(text), collapse = " "))
  gsub("\\b([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-[0-9]+)?\\b", "<accession>",
    text, perl = TRUE)
}

#' The redacted message of a condition's innermost parent (purrr wraps a request error as `purrr_error_indexed`,
#' whose own message is the index; httr2 wraps the curl error as `httr2_failure`).
#' @param cnd A condition.
#' @return One string.
#' @noRd
uniprotConditionDetail <- function(cnd) {
  while (inherits(cnd$parent, "condition")) cnd <- cnd$parent
  redactUniprotDetail(conditionMessage(cnd))
}

#' Whether an error is a transport failure of an httr2 request (a condition of class `httr2_failure` in its chain).
#' @param cnd A condition.
#' @return TRUE or FALSE.
#' @noRd
isUniprotTransportError <- function(cnd) {
  while (inherits(cnd, "condition")) {
    if (inherits(cnd, "httr2_failure")) return(TRUE)
    cnd <- cnd$parent
  }
  FALSE
}

#' The same condition object with the class `prot_annotation_uniprot_recorded` prepended and its fields set: its
#' message and existing classes are unchanged. The class says that the outcome is already in the record.
#' @param cnd A condition.
#' @param reason The reason (`request_failed`, `not_attempted` or `processing_failed`).
#' @param fields Further named fields (`detail`, counts).
#' @return The condition.
#' @noRd
markUniprotRecorded <- function(cnd, reason, fields = list()) {
  for (name in names(fields)) cnd[[name]] <- fields[[name]]
  cnd$reason <- reason
  class(cnd) <- unique(c("prot_annotation_uniprot_recorded", class(cnd)))
  cnd
}

#' Record the outcome of a UniProt retrieval in the analysis record (`external_resources`, service `uniprot`).
#' @param release The release header, or NA.
#' @param release_kind `header`, `failed` or `not_queried`.
#' @param summary The query summary (counts, never identifiers).
#' @return NULL invisibly.
#' @noRd
recordUniprotOutcome <- function(release, release_kind, summary) {
  manager <- getDIAActiveManager()
  if (!usesDIAReferenceStates(manager)) return(invisible(NULL))
  recordWorkflowExternalResource(manager, "uniprot", release, release_kind, query_summary = summary)
}

#' The number of rows of the active owner's `external_resources` (0 without an owner).
#' @noRd
uniprotRecordRows <- function() {
  manager <- getDIAActiveManager()
  if (!usesDIAReferenceStates(manager)) return(0L)
  length(manager$.context$snapshot$external_resources)
}

#' The statuses of a retrieval's refused requests as text, or "" when there are none.
#' @noRd
uniprotStatusText <- function(statuses) {
  if (length(statuses)) paste(statuses, collapse = ", ") else ""
}

#' A record row's `query_summary` with a reference annotation's label appended (KIT-167 8.3 item 3, KIT-175 B15)
#' @param summary KIT-175's text.
#' @param label NULL (today's text), or the label.
#' @return One string.
#' @noRd
uniprotLabelSummary <- function(summary, label = NULL) {
  if (is.null(label)) summary else paste0(summary, "; ", label)
}

#' The record's `query_summary` of a retrieval that failed or was incomplete (KIT-175 section 12).
#' @param reason `service_status`, `invalid_response` or `incomplete`.
#' @param counts The request counts.
#' @return One string.
#' @noRd
uniprotRecordSummary <- function(reason, counts) {
  n <- counts$requests_total %||% 0L
  r <- counts$requests_not_sent %||% 0L
  k <- (counts$requests_failed %||% 0L) - r
  t <- counts$requests_truncated %||% 0L
  status <- uniprotStatusText(counts$statuses)
  not_sent <- if (r > 0) sprintf("; %d requests not sent", r) else ""
  switch(reason,
    service_status = sprintf(paste("annotation failed: UniProt refused %d of %d requests (status %s) and returned no",
      "entry for the others; 0 of %d accessions annotated%s"), k, n, status, counts$accessions_requested, not_sent),
    invalid_response = sprintf("annotation failed: the answer was not a complete UniProt table; 0 of %d accessions annotated",
      counts$accessions_requested),
    incomplete = paste0(
      sprintf("%d of %d requests failed%s; %d of %d accessions not looked up", k, n,
        if (nzchar(status)) sprintf(" (status %s)", status) else "", counts$accessions_not_looked_up,
        counts$accessions_requested),
      if (t > 0) sprintf(" (%d answers cut short with a requested accession missing)", t) else "",
      not_sent))
}

#' The reason of an empty annotation in words, for the condition's message.
#' @noRd
uniprotReasonWords <- function(reason) {
  switch(reason,
    request_failed = "UniProt could not be reached",
    service_status = "UniProt answered with an error",
    invalid_response = "the answer was not a complete UniProt table",
    no_rows = "UniProt has no entry for the accessions",
    no_uniprot_accessions = "no identifier is a UniProt accession",
    id_conversion_failed = "the Ensembl identifiers could not be converted to UniProt accessions",
    processing_failed = "the answer could not be processed",
    not_attempted = "annotation was not attempted",
    "the reason is not known")
}

#' Signal that a UniProt retrieval ended with zero rows (KIT-175 item 1): a warning of class
#' `prot_annotation_uniprot_empty` with the fields `reason`, `requested`, `detail`, `taxon_id` and, where requests
#' were made, the request counts.
#' @param reason The reason.
#' @param requested The number of identifiers asked about.
#' @param detail The underlying text (redacted), or NA.
#' @param taxon_id The organism queried.
#' @param counts The request counts (`requests_total`, `requests_failed`, `requests_not_sent`, `requests_truncated`,
#'   `hits_beyond_rows`, `statuses`, ...), or an empty list.
#' @return NULL invisibly.
#' @noRd
signalUniprotEmpty <- function(reason, requested, detail = NA_character_, taxon_id = NA, counts = list()) {
  text <- sprintf("UniProt annotation is empty: %s (0 of %s accessions annotated)", uniprotReasonWords(reason),
    format(requested))
  do.call(rlang::warn, c(list(message = text, class = "prot_annotation_uniprot_empty", reason = reason,
    requested = requested, detail = detail, taxon_id = taxon_id), counts))
  invisible(NULL)
}

#' Signal that a UniProt retrieval is incomplete (KIT-175 item 1a): a warning of class
#' `prot_annotation_uniprot_incomplete` with the request counts.
#' @param requested The number of accessions asked about.
#' @param taxon_id The organism queried.
#' @param counts The request counts.
#' @return NULL invisibly.
#' @noRd
signalUniprotIncomplete <- function(requested, taxon_id = NA, counts = list()) {
  text <- sprintf("UniProt annotation is incomplete: %d of %d requests failed, %d of %d accessions not looked up",
    counts$requests_failed, counts$requests_total, counts$accessions_not_looked_up, counts$accessions_requested)
  do.call(rlang::warn, c(list(message = text, class = "prot_annotation_uniprot_incomplete", reason = "incomplete",
    requested = requested, taxon_id = taxon_id), counts))
  invisible(NULL)
}

#' The outcome of an annotation retrieval as a plain named list (the attribute `annotation_outcome`), built from the
#' classed condition: not the condition object, so the session store's payload check accepts it.
#' @param cnd A condition of class `prot_annotation_uniprot_empty`, `prot_annotation_uniprot_incomplete` or
#'   `prot_annotation_uniprot_recorded`, or NULL.
#' @return A named list, or NULL.
#' @noRd
uniprotOutcomeFromCondition <- function(cnd) {
  if (is.null(cnd)) return(NULL)
  fields <- c("reason", "stage", "requested", "detail", "taxon_id", "requests_total", "requests_failed",
    "requests_not_sent", "requests_truncated", "accessions_requested", "accessions_not_looked_up", "hits_beyond_rows",
    "statuses")
  values <- unclass(cnd)[intersect(fields, names(cnd))]
  values <- values[!vapply(values, is.null, logical(1))]
  c(list(class = if (inherits(cnd, "prot_annotation_uniprot_incomplete")) "prot_annotation_uniprot_incomplete" else
    "prot_annotation_uniprot_empty"), values)
}

#' The first sentence of the annotation notice for an outcome (KIT-175 section 12; wording VAL-41, the owner's to edit).
#' @param outcome The outcome list (`uniprotOutcomeFromCondition`), or NULL for "no reason known".
#' @return One string.
#' @noRd
uniprotAnnotationSentence <- function(outcome) {
  reason <- outcome$reason %||% ""
  detail <- outcome$detail %||% NA_character_
  detail <- if (length(detail) && !is.na(detail[[1L]])) detail[[1L]] else "no detail"
  n <- outcome$requests_total %||% 0L
  r <- outcome$requests_not_sent %||% 0L
  k <- (outcome$requests_failed %||% 0L) - r
  t <- outcome$requests_truncated %||% 0L
  status <- uniprotStatusText(outcome$statuses)
  switch(reason,
    request_failed = sprintf("Protein annotation was not retrieved: UniProt could not be reached (%s).", detail),
    service_status = paste0(
      sprintf("Protein annotation was not retrieved: UniProt answered with an error (status %s) for %d of %d requests.",
        status, k, n),
      if (r > 0) sprintf(" UniProt stopped answering, so the remaining %d requests were not sent.", r) else ""),
    no_rows = sprintf(paste("Protein annotation is empty: UniProt has no entry for any of the %s accessions under",
      "organism %s. Check the organism (taxon ID)."), format(outcome$requested), format(outcome$taxon_id)),
    no_uniprot_accessions = if (identical(as.integer(outcome$requested %||% NA), 0L)) {
      paste("Protein annotation could not be requested at this step: it has no protein identifiers to look up.",
        "Save the design again to retrieve the annotation.")
    } else {
      sprintf("Protein annotation was not requested: none of the %s protein identifiers is a UniProt accession.",
        format(outcome$requested))
    },
    incomplete = {
      clauses <- c(
        if (k > 0) sprintf("%d of %d requests to UniProt failed%s", k, n,
          if (nzchar(status)) sprintf(" (status %s)", status) else ""),
        if (r > 0) sprintf("%d more were not sent after UniProt stopped answering", r),
        if (t > 0) sprintf("%d answers were cut short with a requested protein missing", t))
      lead <- if (k == 0 && r == 0) {
        sprintf("%d of %d answers from UniProt were cut short with a requested protein missing", t, n)
      } else {
        paste(clauses, collapse = "; ")
      }
      sprintf("Protein annotation is incomplete: %s, so %d of %d proteins were not looked up.", lead,
        outcome$accessions_not_looked_up, outcome$accessions_requested)
    },
    not_attempted = sprintf("Protein annotation was not attempted: %s.", detail),
    invalid_response = paste("Protein annotation was not retrieved: the answer received was not UniProt's complete",
      "table (for example a network sign-in page, or an answer cut short)."),
    id_conversion_failed = sprintf(
      "Protein annotation was not retrieved: the Ensembl identifiers could not be converted to UniProt accessions (%s).",
      detail),
    processing_failed = if (identical(outcome$stage, "save")) {
      sprintf("Protein annotation was retrieved but could not be saved to the project (%s).", detail)
    } else {
      sprintf("Protein annotation was retrieved but could not be processed (%s).", detail)
    },
    "Protein annotation is empty.")
}

#' The design step's annotation notice for an outcome: the first sentence, the consequence and, where a retry can
#' cure it, the retry advice (KIT-175 section 12).
#' @param outcome The outcome list, or NULL for "no reason known".
#' @return One string.
#' @noRd
uniprotAnnotationNotice <- function(outcome) {
  reason <- outcome$reason %||% ""
  consequence <- if (identical(reason, "incomplete")) {
    "The analysis can continue. The proteins not looked up have no gene names and are left out of enrichment."
  } else {
    paste("The analysis can continue. Without annotation the results tables have no gene names, plots are labelled",
      "with accessions, and enrichment cannot be run.")
  }
  retry <- if (reason %in% c("request_failed", "service_status", "invalid_response", "incomplete", "processing_failed")) {
    " To try again, save the design again."
  } else {
    ""
  }
  paste0(uniprotAnnotationSentence(outcome), " ", consequence, retry)
}
