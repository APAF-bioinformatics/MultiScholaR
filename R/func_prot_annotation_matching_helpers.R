# ----------------------------------------------------------------------------
# matchAnnotations
# ----------------------------------------------------------------------------
#' Match Protein Annotations to DE Results
#'
#' @description Matches protein identifiers from differential expression results
#' with UniProt annotations, enriching the DE results with gene names and other
#' annotation information.
#'
#' @param da_results_s4 S4 object containing DE results from differential expression analysis
#' @param uniprot_annotations Data frame with UniProt annotations (typically uniprot_dat_cln)
#' @param protein_id_column Character string specifying the protein ID column name in DE results
#' @param uniprot_id_column Character string specifying the UniProt ID column in annotations (default: "Entry")
#' @param gene_names_column Character string specifying the gene names column in annotations (default: "gene_names")
#'
#' @return List containing:
#' \itemize{
#'   \item annotated_da_results: DE results enriched with UniProt annotations
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
#'   da_results_s4 = da_results_obj,
#'   uniprot_annotations = uniprot_dat_cln,
#'   protein_id_column = "uniprot_acc"
#' )
#' 
#' # Check matching statistics
#' print(matched_results$match_statistics)
#' }
#'
# findFuzzyAnnotationMatches
# ----------------------------------------------------------------------------
findFuzzyAnnotationMatches <- function(unmatched_proteins,
                                       uniprot_annotations,
                                       uniprot_id_column) {
  fuzzy_matches <- data.frame()

  if (nrow(unmatched_proteins) == 0) {
    return(fuzzy_matches)
  }

  log_info("Performing fuzzy matching for remaining proteins...")

  version_cleaned <- unmatched_proteins |>
    dplyr::mutate(
      version_cleaned = stringr::str_replace(.data$cleaned_id, "\\.\\d+$", "")
    ) |>
    dplyr::inner_join(
      uniprot_annotations,
      by = setNames(uniprot_id_column, "version_cleaned")
    ) |>
    dplyr::mutate(!!uniprot_id_column := .data$version_cleaned)

  if (nrow(version_cleaned) > 0) {
    fuzzy_matches <- dplyr::bind_rows(fuzzy_matches, version_cleaned)
    log_info(sprintf("Version-based matches found: %d", nrow(version_cleaned)))
  }

  still_unmatched <- unmatched_proteins |>
    dplyr::anti_join(fuzzy_matches, by = "original_id")

  if (nrow(still_unmatched) == 0) {
    return(fuzzy_matches)
  }

  log_info("Attempting substring matching for remaining proteins...")

  if (nrow(still_unmatched) > 1000) {
    log_warn(sprintf("Skipping substring matching for %d proteins (too many for efficient processing)", nrow(still_unmatched)))
    return(fuzzy_matches)
  }

  substring_matches <- still_unmatched |>
    dplyr::rowwise() |>
    dplyr::mutate(
      matched_uniprot = {
        matches <- stringr::str_detect(uniprot_annotations[[uniprot_id_column]], .data$cleaned_id)
        if (any(matches)) {
          uniprot_annotations[[uniprot_id_column]][matches][1]
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

  fuzzy_matches
}

# ----------------------------------------------------------------------------
# matchAnnotations
# ----------------------------------------------------------------------------
#' @export
matchAnnotations <- function(da_results_s4,
                            uniprot_annotations,
                            protein_id_column = NULL,
                            uniprot_id_column = "Entry",
                            gene_names_column = "gene_names") {
  
  log_info("=== Starting protein annotation matching ===")
  
  # Validate inputs
  if (is.null(da_results_s4)) {
    stop("da_results_s4 cannot be NULL")
  }
  
  if (is.null(uniprot_annotations) || nrow(uniprot_annotations) == 0) {
    stop("uniprot_annotations cannot be NULL or empty")
  }
  
  # Get protein ID column from S4 object if not specified
  if (is.null(protein_id_column)) {
    if (!is.null(da_results_s4@protein_id_column)) {
      protein_id_column <- da_results_s4@protein_id_column
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
  # OR from de_data for da_results_for_enrichment objects
  protein_ids <- NULL
  
  # [OK] ENHANCED: Handle da_results_for_enrichment objects
  if (inherits(da_results_s4, "da_results_for_enrichment")) {
    log_info("Detected da_results_for_enrichment S4 object")
    
    # Extract protein IDs from all DE data contrasts
    if (!is.null(da_results_s4@da_data) && length(da_results_s4@da_data) > 0) {
      # Get protein IDs from the first available contrast (they should all have the same proteins)
      available_contrasts <- names(da_results_s4@da_data)
      log_info(sprintf("Available DE contrasts: %s", paste(available_contrasts, collapse = ", ")))
      
      for (contrast_name in available_contrasts) {
        contrast_data <- da_results_s4@da_data[[contrast_name]]
        if (!is.null(contrast_data) && nrow(contrast_data) > 0 && protein_id_column %in% names(contrast_data)) {
          protein_ids <- unique(contrast_data[[protein_id_column]])
          log_info(sprintf("Extracted protein IDs from DE contrast: %s", contrast_name))
          break  # Use first available contrast with valid data
        }
      }
      
      if (is.null(protein_ids)) {
        stop(sprintf("Protein ID column '%s' not found in any DE data contrasts. Available columns in first contrast: %s", 
                     protein_id_column, 
                     if(length(available_contrasts) > 0 && !is.null(da_results_s4@da_data[[available_contrasts[1]]])) {
                       paste(names(da_results_s4@da_data[[available_contrasts[1]]]), collapse = ", ")
                     } else {
                       "No valid data found"
                     }))
      }
    } else {
      stop("No DE data found in da_results_for_enrichment object @da_data slot")
    }
    
  } else {
    # [OK] PRESERVED: Original logic for ProteinQuantitativeData and similar objects
    if (!is.null(da_results_s4@protein_quant_table) && protein_id_column %in% names(da_results_s4@protein_quant_table)) {
      protein_ids <- unique(da_results_s4@protein_quant_table[[protein_id_column]])
      log_info("Extracted protein IDs from @protein_quant_table")
    } else if (!is.null(da_results_s4@protein_id_table) && protein_id_column %in% names(da_results_s4@protein_id_table)) {
      protein_ids <- unique(da_results_s4@protein_id_table[[protein_id_column]])
      log_info("Extracted protein IDs from @protein_id_table")
    } else {
      stop(sprintf("Protein ID column '%s' not found in DE results S4 object. Object class: %s", 
                   protein_id_column, class(da_results_s4)[1]))
    }
  }
  
  log_info(sprintf("Found %d unique protein IDs in DE results", length(protein_ids)))
  
  # Clean protein IDs for matching
  # Remove isoform suffixes (-1, -2, etc.) and take first protein from groups
  cleaned_protein_ids <- protein_ids |>
    purrr::map_chr(\(x) {
      # Split by semicolon and take first protein
      first_protein <- stringr::str_split(x, ";")[[1]][1]
      # [OK] REFACTORED: Use centralized UniProt normalization
      normalizeUniprotAccession(first_protein, remove_isoform = TRUE)
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
    dplyr::inner_join(
      uniprot_annotations,
      by = setNames(uniprot_id_column, "cleaned_id")
    ) |>
    dplyr::mutate(!!uniprot_id_column := .data$cleaned_id)
  
  log_info(sprintf("Exact matches found: %d", nrow(exact_matches)))
  
  # For unmatched proteins, try partial matching strategies
  unmatched_proteins <- protein_mapping |>
    dplyr::anti_join(exact_matches, by = "original_id")
  
  log_info(sprintf("Proteins requiring fuzzy matching: %d", nrow(unmatched_proteins)))
  
  fuzzy_matches <- findFuzzyAnnotationMatches(
    unmatched_proteins = unmatched_proteins,
    uniprot_annotations = uniprot_annotations,
    uniprot_id_column = uniprot_id_column
  )
  
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
    annotated_da_results = annotated_results,
    match_statistics = match_statistics,
    unmatched_proteins = unmatched_protein_list
  ))
}

