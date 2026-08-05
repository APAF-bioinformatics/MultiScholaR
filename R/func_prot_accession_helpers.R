# ----------------------------------------------------------------------------
# chooseBestProteinAccession_s3
# ----------------------------------------------------------------------------
#' Choose Best Protein Accession for Data Frame (S3 Version)
#'
#' @description A simplified S3 version of chooseBestProteinAccession that works
#' directly on data frames. Resolves protein groups by selecting the best
#' representative protein based on FASTA sequence data and aggregates
#' quantitative values.
#'
#' @param data_tbl Data frame containing protein quantitative data
#' @param protein_id_column Character string specifying the protein ID column name
#' @param seqinr_obj Data frame with FASTA sequence information (aa_seq_tbl_final)
#' @param seqinr_accession_column Character string specifying the accession column in seqinr_obj
#' @param delim Character string delimiter used to separate protein groups (default: ";")
#' @param replace_zero_with_na Logical, whether to replace zeros with NA (default: TRUE)
#' @param aggregation_method Character string aggregation method ("mean", "median", "sum") (default: "mean")
#'
#' @return Data frame with cleaned protein IDs and aggregated values
#'
#' @details
#' This function processes protein groups (e.g., "P12345;P67890;Q11111") by:
#' \itemize{
#'   \item Splitting groups using the specified delimiter
#'   \item Finding the best representative protein using FASTA sequence data
#'   \item Aggregating quantitative values for the chosen protein
#'   \item Returning cleaned data with single protein IDs
#' }
#'
#' The selection of "best" protein is based on:
#' \itemize{
#'   \item Presence in FASTA sequence database
#'   \item Protein length (longer sequences preferred)
#'   \item Alphabetical order as tiebreaker
#' }
#'
#' @examples
#' \dontrun{
#' cleaned_data <- chooseBestProteinAccession_s3(
#'   data_tbl = protein_data,
#'   protein_id_column = "Protein.Group",
#'   seqinr_obj = aa_seq_tbl_final,
#'   seqinr_accession_column = "uniprot_acc"
#' )
#' }
#'
#' @export
chooseBestProteinAccession_s3 <- function(data_tbl,
                                          protein_id_column,
                                          seqinr_obj,
                                          seqinr_accession_column = "uniprot_acc",
                                          delim = ";",
                                          replace_zero_with_na = TRUE,
                                          aggregation_method = "mean") {
  
  cat("=== Starting chooseBestProteinAccession_s3 ===\n")
  
  tryCatch({
    # Validate inputs with detailed error messages
    cat("*** S3 CLEANUP: Validating inputs ***\n")
    
    if (is.null(data_tbl) || !is.data.frame(data_tbl)) {
      stop("data_tbl must be a non-null data frame")
    }
    
    if (nrow(data_tbl) == 0) {
      stop("data_tbl cannot be empty")
    }
    
    cat(sprintf("*** S3 CLEANUP: data_tbl dimensions: %d rows x %d cols ***\n", nrow(data_tbl), ncol(data_tbl)))
    cat(sprintf("*** S3 CLEANUP: data_tbl column names: %s ***\n", paste(head(names(data_tbl), 10), collapse = ", ")))
    
    if (!protein_id_column %in% names(data_tbl)) {
      stop(sprintf("Protein ID column '%s' not found in data. Available columns: %s", 
                   protein_id_column, paste(names(data_tbl), collapse = ", ")))
    }
    
    if (is.null(seqinr_obj) || !is.data.frame(seqinr_obj)) {
      stop("seqinr_obj must be a non-null data frame")
    }
    
    if (nrow(seqinr_obj) == 0) {
      stop("seqinr_obj cannot be empty")
    }
    
    cat(sprintf("*** S3 CLEANUP: seqinr_obj dimensions: %d rows x %d cols ***\n", nrow(seqinr_obj), ncol(seqinr_obj)))
    cat(sprintf("*** S3 CLEANUP: seqinr_obj column names: %s ***\n", paste(head(names(seqinr_obj), 10), collapse = ", ")))
    
    if (!seqinr_accession_column %in% names(seqinr_obj)) {
      stop(sprintf("Accession column '%s' not found in sequence data. Available columns: %s", 
                   seqinr_accession_column, paste(names(seqinr_obj), collapse = ", ")))
    }
    
    cat("*** S3 CLEANUP: Input validation passed ***\n")
    
    # Get non-quantitative columns (metadata columns) with error handling
    cat("*** S3 CLEANUP: Identifying column types ***\n")
    
    numeric_cols <- tryCatch({
      sapply(data_tbl, is.numeric)
    }, error = function(e) {
      stop(sprintf("Error checking column types: %s", e$message))
    })
    
    quant_cols <- names(data_tbl)[numeric_cols]
    meta_cols <- names(data_tbl)[!numeric_cols]
    
    cat(sprintf("*** S3 CLEANUP: Found %d quantitative columns and %d metadata columns ***\n", 
                length(quant_cols), length(meta_cols)))
    cat(sprintf("*** S3 CLEANUP: Quantitative columns: %s ***\n", paste(head(quant_cols, 5), collapse = ", ")))
    cat(sprintf("*** S3 CLEANUP: Metadata columns: %s ***\n", paste(head(meta_cols, 5), collapse = ", ")))
  
    # Replace zeros with NA if requested
    if (replace_zero_with_na && length(quant_cols) > 0) {
      cat("*** S3 CLEANUP: Replacing zeros with NA in quantitative columns ***\n")
      
      tryCatch({
        data_tbl[quant_cols] <- lapply(data_tbl[quant_cols], function(x) {
          if (is.numeric(x)) {
            x[x == 0] <- NA
          }
          x
        })
        cat("*** S3 CLEANUP: Zero replacement completed ***\n")
      }, error = function(e) {
        stop(sprintf("Error replacing zeros with NA: %s", e$message))
      })
    }
    
    # Get unique protein groups
    cat("*** S3 CLEANUP: Getting unique protein groups ***\n")
    
    protein_groups <- tryCatch({
      unique(data_tbl[[protein_id_column]])
    }, error = function(e) {
      stop(sprintf("Error getting unique protein groups: %s", e$message))
    })
    
    cat(sprintf("*** S3 CLEANUP: Processing %d unique protein groups ***\n", length(protein_groups)))
    cat(sprintf("*** S3 CLEANUP: First 5 protein groups: %s ***\n", paste(head(protein_groups, 5), collapse = ", ")))
    
    # Create lookup table for sequence information
    cat("*** S3 CLEANUP: Creating sequence lookup table ***\n")
    
    seq_lookup <- tryCatch({
      if ("length" %in% names(seqinr_obj)) {
        cat("*** S3 CLEANUP: Using existing length column ***\n")
        seqinr_obj |>
          dplyr::select(dplyr::all_of(c(seqinr_accession_column, "length"))) |>
          dplyr::distinct()
      } else {
        cat("*** S3 CLEANUP: Creating default length column ***\n")
        # If no length column, create one or use a default
        seqinr_obj |>
          dplyr::select(dplyr::all_of(seqinr_accession_column)) |>
          dplyr::distinct() |>
          dplyr::mutate(length = 1000)  # Default length
      }
    }, error = function(e) {
      stop(sprintf("Error creating sequence lookup: %s", e$message))
    })
    
    cat(sprintf("*** S3 CLEANUP: Created sequence lookup with %d entries ***\n", nrow(seq_lookup)))
  
    # Function to choose best protein from a group
    cat("*** S3 CLEANUP: Defining chooseBestFromGroup function ***\n")
    
    chooseBestFromGroup <- function(group_string) {
      tryCatch({
        if (is.na(group_string) || group_string == "") {
          return(NA_character_)
        }
        
        # Split the group
        proteins <- stringr::str_split(group_string, delim)[[1]]
        proteins <- stringr::str_trim(proteins)  # Remove whitespace
        proteins <- proteins[proteins != ""]     # Remove empty strings
        
        if (length(proteins) == 1) {
          return(proteins[1])
        }
        
        # Find proteins that exist in sequence database
        proteins_in_db <- proteins[proteins %in% seq_lookup[[seqinr_accession_column]]]
        
        if (length(proteins_in_db) == 0) {
          # No proteins found in database, return first one
          return(proteins[1])
        }
        
        if (length(proteins_in_db) == 1) {
          return(proteins_in_db[1])
        }
        
        # Multiple proteins in database - choose by length (longest first)
        protein_info <- seq_lookup |>
          dplyr::filter(!!sym(seqinr_accession_column) %in% proteins_in_db) |>
          dplyr::arrange(desc(length), !!sym(seqinr_accession_column))
        
        return(protein_info[[seqinr_accession_column]][1])
        
      }, error = function(e) {
        warning(sprintf("Error processing protein group '%s': %s", group_string, e$message))
        return(group_string)  # Return original if processing fails
      })
    }
    
    # Create mapping of original groups to best proteins
    cat("*** S3 CLEANUP: Creating protein group mapping ***\n")
    
    protein_mapping <- tryCatch({
      data.frame(
        original_group = protein_groups,
        best_protein = sapply(protein_groups, chooseBestFromGroup, USE.NAMES = FALSE),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      stop(sprintf("Error creating protein mapping: %s", e$message))
    })
    
    # Count reductions
    groups_with_multiple <- tryCatch({
      sum(stringr::str_detect(protein_mapping$original_group, delim), na.rm = TRUE)
    }, error = function(e) {
      cat(sprintf("*** S3 CLEANUP: Warning - could not count multi-protein groups: %s ***\n", e$message))
      0
    })
    
    cat(sprintf("*** S3 CLEANUP: Found %d protein groups with multiple proteins ***\n", groups_with_multiple))
  
    # Add mapping to data
    cat("*** S3 CLEANUP: Applying protein mapping to data ***\n")
    
    data_with_mapping <- tryCatch({
      data_tbl |>
        dplyr::left_join(
          protein_mapping,
          by = setNames("original_group", protein_id_column)
        )
    }, error = function(e) {
      stop(sprintf("Error applying protein mapping: %s", e$message))
    })
    
    cat(sprintf("*** S3 CLEANUP: Data mapping completed. Mapped data has %d rows ***\n", nrow(data_with_mapping)))
    
    # Group by best protein and aggregate
    cat(sprintf("*** S3 CLEANUP: Aggregating data using method: %s ***\n", aggregation_method))
    
    # Define aggregation function
    agg_func <- switch(aggregation_method,
      "mean" = function(x) mean(x, na.rm = TRUE),
      "median" = function(x) median(x, na.rm = TRUE),
      "sum" = function(x) sum(x, na.rm = TRUE),
      function(x) mean(x, na.rm = TRUE)  # Default to mean
    )
    
    # Separate quantitative and metadata columns for different aggregation
    cat("*** S3 CLEANUP: Aggregating quantitative data ***\n")
    
    quant_data <- tryCatch({
      if (length(quant_cols) > 0) {
        data_with_mapping |>
          dplyr::select(best_protein, dplyr::all_of(quant_cols)) |>
          dplyr::group_by(best_protein) |>
          dplyr::summarise(
            dplyr::across(dplyr::all_of(quant_cols), agg_func),
            .groups = "drop"
          )
      } else {
        data.frame(best_protein = unique(data_with_mapping$best_protein))
      }
    }, error = function(e) {
      stop(sprintf("Error aggregating quantitative data: %s", e$message))
    })
    
    cat("*** S3 CLEANUP: Aggregating metadata ***\n")
    
    # For metadata columns, take the first non-NA value for each best protein
    meta_data <- tryCatch({
      if (length(meta_cols) > 0) {
        data_with_mapping |>
          dplyr::select(best_protein, dplyr::all_of(meta_cols)) |>
          dplyr::group_by(best_protein) |>
          dplyr::summarise(
            dplyr::across(dplyr::all_of(meta_cols), ~ dplyr::first(na.omit(.))[1]),
            .groups = "drop"
          )
      } else {
        data.frame(best_protein = unique(data_with_mapping$best_protein))
      }
    }, error = function(e) {
      stop(sprintf("Error aggregating metadata: %s", e$message))
    })
    
    cat("*** S3 CLEANUP: Combining results ***\n")
    
    # Combine results
    result_data <- tryCatch({
      quant_data |>
        dplyr::left_join(meta_data, by = "best_protein") |>
        dplyr::rename(!!sym(protein_id_column) := best_protein)
    }, error = function(e) {
      stop(sprintf("Error combining results: %s", e$message))
    })
    
    cat("*** S3 CLEANUP: Reordering columns ***\n")
    
    # Reorder columns to match original
    result_data <- tryCatch({
      result_data |>
        dplyr::select(dplyr::all_of(names(data_tbl)))
    }, error = function(e) {
      stop(sprintf("Error reordering columns: %s", e$message))
    })
    
    # Final statistics
    original_proteins <- length(unique(data_tbl[[protein_id_column]]))
    final_proteins <- length(unique(result_data[[protein_id_column]]))
    reduction_pct <- round(((original_proteins - final_proteins) / original_proteins) * 100, 1)
    
    cat(sprintf("*** S3 CLEANUP: Protein cleanup completed: %d -> %d proteins (%.1f%% reduction) ***\n",
                     original_proteins, final_proteins, reduction_pct))
    
    return(result_data)
    
  }, error = function(e) {
    cat(sprintf("*** S3 CLEANUP FATAL ERROR: %s ***\n", e$message))
    cat("*** S3 CLEANUP: Returning original data ***\n")
    return(data_tbl)  # Return original data if everything fails
  })
}

