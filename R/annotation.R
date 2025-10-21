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
#' Ensembl protein IDs follow the pattern: ENS[SPECIES]P[0-9]+
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
  
  # Remove potential isoform suffixes for broader matching
  cleaned_proteins <- unique(gsub("-\\d+$", "", unique_proteins))
  
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




