#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#'@export

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

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

#' Get UniProt Annotations from Actual Data Proteins
#'
#' @description Creates a comprehensive UniProt annotation lookup table using
#' protein IDs from the actual imported data. This function extracts unique protein
#' accessions by splitting protein groups and retrieves their annotations efficiently.
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
getUniprotAnnotationsFull <- function(data_tbl,
                                     protein_id_column,
                                     cache_dir,
                                     taxon_id = 9606,
                                     chunk_size = 25) {
  
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
  
  # Remove potential isoform suffixes for broader matching
  cleaned_proteins <- unique(gsub("-\\d+$", "", unique_proteins))
  
  cat(sprintf("Expanded %d protein groups into %d unique individual proteins\n", 
                   length(raw_protein_groups), length(cleaned_proteins)))
  cat(sprintf("First 10 individual proteins: %s\n", paste(head(cleaned_proteins, 10), collapse = ", ")))
  
  # Create a temporary data frame for getUniProtAnnotation
  # This function expects protein IDs in a standard format
  temp_protein_table <- data.frame(
    Protein.Ids = cleaned_proteins,
    stringsAsFactors = FALSE
  )
  
  cat("Created temporary protein table for UniProt annotation retrieval\n")
  cat(sprintf("Temporary table dimensions: %d rows x %d cols\n", nrow(temp_protein_table), ncol(temp_protein_table)))
  
  # Create cache directory if needed
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
    cat(sprintf("Created cache directory: %s\n", cache_dir))
  }
  
  # Call getUniprotAnnotations with our optimized protein list
  cat("Calling getUniprotAnnotations for comprehensive annotation retrieval...\n")
  
  tryCatch({
    uniprot_annotations <- getUniprotAnnotations(
      input_tbl = temp_protein_table,
      cache_dir = cache_dir,
      taxon_id = taxon_id
    )
    
    cat(sprintf("Successfully retrieved UniProt annotations for %d proteins\n", nrow(uniprot_annotations)))
    
    # Add some metadata about the annotation process
    attr(uniprot_annotations, "annotation_source") <- "DATA_OPTIMIZED"
    attr(uniprot_annotations, "original_protein_groups") <- length(raw_protein_groups)
    attr(uniprot_annotations, "unique_proteins_processed") <- length(cleaned_proteins)
    attr(uniprot_annotations, "taxon_id") <- taxon_id
    attr(uniprot_annotations, "processing_timestamp") <- Sys.time()
    
    cat("=== getUniprotAnnotationsFull completed successfully ===\n")
    
    return(uniprot_annotations)
    
  }, error = function(e) {
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




