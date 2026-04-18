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
  
  if (length(protein_ids) == 0) {
    message("   No valid UniProt IDs found to query. Skipping download.")
    return(NULL)
  }
  
  message(paste("Querying", length(protein_ids), "valid UniProt IDs"))
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

