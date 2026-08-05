# runOneStringDbRankEnrichment
# ----------------------------------------------------------------------------
#' Run STRING DB Rank Enrichment Analysis
#'
#' @description
#' Performs STRING DB enrichment analysis on ranked data. This function processes
#' differential expression results by calculating a score from log fold change and
#' FDR values, extracts protein identifiers, and submits them to STRING DB for
#' functional enrichment analysis. Results are automatically saved to files.
#'
#' @param input_table A data frame containing differential expression results with
#'   columns: `log2FC` (log fold change), `fdr_qvalue` (adjusted p-value), and
#'   `Protein.Ids` (protein identifiers, potentially with isoform information).
#' @param result_label Character string: A label for the analysis results, used
#'   in output file names.
#' @param pathway_dir Character string: The pathway directory where results will be saved.
#'   All outputs are saved to pathway_dir/string_db/ subdirectory.
#' @param api_key Character string: Your personal STRING API key.
#' @param species Character or numeric: NCBI/STRING species identifier. Default is 9606 (Homo sapiens).
#' @param ge_fdr Numeric: FDR threshold for gene expression enrichment. Default is 0.05.
#' @param ge_enrichment_rank_direction Integer: Direction for enrichment rank.
#'   (-1, 0, or 1). Default is -1.
#' @param polling_interval_seconds Numeric: Seconds to wait between polling attempts.
#'   Default is 10.
#' @param max_polling_attempts Numeric: Maximum number of polling attempts. Default is 30.
#'
#' @return A data frame containing the STRING DB enrichment results.
#'
#' @details
#' This function saves enrichment results to the specified pathway directory:
#' - `{pathway_dir}/string_db/{result_label}_string_enrichment_page_url.txt` - URLs for STRING DB results
#' - `{pathway_dir}/string_db/{result_label}_string_enrichment_results.tab` - Tab-delimited enrichment results
#' - `{pathway_dir}/string_db/{result_label}_string_enrichment_graph.png` - Network graph image
#'
#' @importFrom dplyr mutate relocate arrange desc
#' @importFrom stringr str_split
#' @importFrom purrr map_chr
#' @importFrom readr write_lines
#' @importFrom vroom vroom_write
#'
#' @examples
#' \dontrun{
#' # Assume 'da_results' is a data frame with log2FC, fdr_qvalue, and Protein.Ids columns
#' enrichment_results <- runOneStringDbRankEnrichment(
#'   input_table = da_results,
#'   result_label = "treatment_vs_control",
#'   pathway_dir = "path/to/pathway_enrichment",
#'   api_key = "YOUR_API_KEY"
#' )
#' }
#'
#' @export
runOneStringDbRankEnrichment <- function( input_table
                                          ,  result_label
                                          , pathway_dir
                                          , api_key = NULL
                                          , species = "9606"
                                          , ge_fdr = 0.05
                                          , ge_enrichment_rank_direction = -1
                                          , polling_interval_seconds = 10
                                          , max_polling_attempts = 30) {

  stringdb_input_table <-  input_table |>
    mutate( score = sign(log2FC) * -log10(fdr_qvalue)) |>
    relocate(score, .after="log2FC") |>
    arrange(desc(log2FC)) |>
    mutate( protein_id = purrr::map_chr(Protein.Ids, ~str_split(.x, ":")[[1]][1])) |>
    relocate(protein_id, .after="Protein.Ids")


  parsed_response <- submitStringDBEnrichment (input_data_frame = stringdb_input_table ,
                                               identifier_column_name = "protein_id",
                                               value_column_name = "score",
                                               caller_identity = result_label,
                                               api_key = api_key,
                                               species = species,
                                               ge_fdr = ge_fdr,
                                               ge_enrichment_rank_direction = ge_enrichment_rank_direction)


  output_tbl <- retrieveStringDBEnrichmentResults( submission_info = parsed_response,
                                                   polling_interval_seconds = polling_interval_seconds,
                                                   max_polling_attempts = max_polling_attempts)

  enrichment_dir <- file.path(pathway_dir, "string_db")
  dir.create(enrichment_dir, showWarnings = FALSE, recursive = TRUE)
  
  write_lines(c("page_url", output_tbl$page_url
                , "download_url" , output_tbl$download_url
                , "graph_url" , output_tbl$graph_url)
              , file.path(enrichment_dir, paste0(result_label, "_string_enrichment_page_url.txt")))

  vroom::vroom_write( output_tbl$enrichment_data
                      , file = file.path(enrichment_dir
                                          , paste0(result_label, "_string_enrichment_results.tab")))

  writeBin(output_tbl$graph_image_content
           , file.path(enrichment_dir, paste0(result_label, "_string_enrichment_graph.png")))

  return(output_tbl$enrichment_data)

}


# ----------------------------------------------------------------------------
# runOneStringDbRankEnrichmentMofa
# ----------------------------------------------------------------------------
#' Run STRING DB Rank Enrichment Analysis for MOFA
#'
#' @description
#' Performs STRING DB enrichment analysis for MOFA (Multi-Omics Factor Analysis)
#' results. This function takes pre-processed data with identifier and value columns
#' and submits them to STRING DB for functional enrichment analysis. Results are
#' automatically saved to the specified results directory.
#'
#' @param input_table A data frame containing the data to analyze with identifier
#'   and value columns as specified by the column name parameters.
#' @param identifier_column_name Character string: The name of the column in
#'   `input_table` that contains the protein/gene identifiers. Default is "protein_id".
#' @param value_column_name Character string: The name of the column in `input_table`
#'   that contains the numerical values (e.g., scores, weights) associated with
#'   each identifier. Default is "score".
#' @param result_label Character string: A label for the analysis results, used
#'   in output file names.
#' @param results_dir Character string: The directory path where results should be saved.
#' @param api_key Character string: Your personal STRING API key.
#' @param species Character or numeric: NCBI/STRING species identifier. Default is 9606 (Homo sapiens).
#' @param ge_fdr Numeric: FDR threshold for gene expression enrichment. Default is 0.05.
#' @param ge_enrichment_rank_direction Integer: Direction for enrichment rank.
#'   (-1, 0, or 1). Default is -1.
#' @param polling_interval_seconds Numeric: Seconds to wait between polling attempts.
#'   Default is 10.
#' @param max_polling_attempts Numeric: Maximum number of polling attempts. Default is 30.
#'
#' @return A data frame containing the STRING DB enrichment results.
#'
#' @details
#' This function saves results directly to the specified `results_dir`:
#' - URLs for STRING DB results page, download, and graph
#' - Tab-delimited enrichment results file
#' - PNG image of the enrichment network graph
#'
#' Unlike `runOneStringDbRankEnrichment`, this function expects the input data to
#' already be properly formatted with identifier and value columns.
#'
#' @importFrom readr write_lines
#' @importFrom vroom vroom_write
#'
#' @examples
#' \dontrun{
#' # Assume 'mofa_results' is a data frame with protein_id and score columns
#' enrichment_results <- runOneStringDbRankEnrichmentMofa(
#'   input_table = mofa_results,
#'   identifier_column_name = "protein_id",
#'   value_column_name = "score",
#'   result_label = "MOFA_factor_1",
#'   results_dir = "/path/to/results",
#'   api_key = "YOUR_API_KEY"
#' )
#' }
#'
#' @export
runOneStringDbRankEnrichmentMofa <- function( input_table
                                              ,   identifier_column_name = "protein_id"
                                              ,   value_column_name = "score"
                                              ,  result_label
                                              , results_dir
                                              , api_key = NULL
                                              , species = "9606"
                                              , ge_fdr = 0.05
                                              , ge_enrichment_rank_direction = -1
                                              , polling_interval_seconds = 10
                                              , max_polling_attempts = 30) {



  parsed_response <- submitStringDBEnrichment (input_data_frame = input_table ,
                                               identifier_column_name = identifier_column_name,
                                               value_column_name = value_column_name,
                                               caller_identity = result_label,
                                               api_key = api_key,
                                               species = species,
                                               ge_fdr = ge_fdr,
                                               ge_enrichment_rank_direction = ge_enrichment_rank_direction)


  output_tbl <- retrieveStringDBEnrichmentResults( submission_info = parsed_response,
                                                   polling_interval_seconds = polling_interval_seconds,
                                                   max_polling_attempts = max_polling_attempts)

  write_lines(c("page_url", output_tbl$page_url
                , "download_url" , output_tbl$download_url
                , "graph_url" , output_tbl$graph_url)
              , file.path( results_dir,  paste0( result_label, "_string_enrichment_page_url.txt") ))

  vroom::vroom_write( output_tbl$enrichment_data
                      , file = file.path( results_dir

                                          , paste0( result_label, "_string_enrichment_results.tab") ))

  dir.create( file.path( results_dir), showWarnings = TRUE, recursive = TRUE)
  writeBin(output_tbl$graph_image_content
           , file.path( results_dir , paste0( result_label, "string_enrichment_graph.png") ))

  return(output_tbl$enrichment_data)

}


# ----------------------------------------------------------------------------
# runStringDbEnrichmentFromDAResults
# ----------------------------------------------------------------------------
#' Run STRING DB Enrichment Analysis from DE Results S4 Object
#'
#' @description
#' This function extracts differential expression data from a da_results_for_enrichment S4 object,
#' applies the specified ranking method, and performs STRING DB enrichment analysis using the
#' existing runOneStringDbRankEnrichmentMofa function.
#'
#' @param da_results_for_enrichment An S4 object of class da_results_for_enrichment containing
#'   differential expression results across multiple contrasts.
#' @param contrast_name Character string. Name of the specific contrast to analyze.
#'   Must match one of the names in da_results_for_enrichment@da_data.
#'   If NULL, will use the first contrast. Default: NULL.
#' @param ranking_method Character string. Method for ranking proteins. Options:
#'   - "fdr_qvalue": Rank by FDR q-value (ascending, most significant first)
#'   - "log2fc": Rank by log2 fold change (descending, highest FC first)  
#'   - "combined_score": Use sign(log2FC) * (-log10(fdr_qvalue)) for ranking
#'   - "none": No ranking applied (proteins in original order)
#'   Default: "combined_score".
#' @param identifier_column Character string. Name of the column containing protein identifiers.
#'   Default: "Protein.Ids".
#' @param filter_significant Logical. Whether to filter to only significant proteins (fdr_qvalue < 0.05).
#'   Default: FALSE (include all proteins).
#' @param fdr_threshold Numeric. FDR threshold for filtering significant proteins when filter_significant=TRUE.
#'   Default: 0.05.
#' @param result_label Character string. A label used for naming output files.
#'   If NULL, will use the contrast name. Default: NULL.
#' @param results_dir Character string. The path to the directory where enrichment
#'   results will be saved. Default: "string_enrichment_results".
#' @param api_key Character string. Your personal STRING API key.
#'   Default: NULL.
#' @param species Character string. NCBI/STRING species identifier.
#'   Default: "9606" (Homo sapiens).
#' @param ge_fdr Numeric. FDR threshold for gene expression enrichment.
#'   Default: 0.05.
#' @param ge_enrichment_rank_direction Integer. Direction for enrichment rank
#'   (-1, 0, or 1). Default: -1.
#' @param polling_interval_seconds Numeric. Seconds to wait between polling attempts.
#'   Default: 10.
#' @param max_polling_attempts Numeric. Maximum polling attempts before timing out.
#'   Default: 30.
#'
#' @return A data frame containing the enrichment results from STRING DB.
#'   Returns NULL if the process fails.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming you have a da_results_for_enrichment object
#' enrichment_results <- runStringDbEnrichmentFromDAResults(
#'   da_results_for_enrichment = my_da_results,
#'   contrast_name = "T2.minus.MSO=groupT2-groupMSO",
#'   ranking_method = "combined_score",
#'   filter_significant = FALSE,
#'   result_label = "T2_vs_MSO_enrichment",
#'   results_dir = "string_enrichment_output",
#'   api_key = "your_api_key",
#'   species = "9606"
#' )
#' }
runStringDbEnrichmentFromDAResults <- function(da_results_for_enrichment,
                                             contrast_name = NULL,
                                             ranking_method = "combined_score",
                                             identifier_column = "Protein.Ids",
                                             filter_significant = FALSE,
                                             fdr_threshold = 0.05,
                                             result_label = NULL,
                                             results_dir = "string_enrichment_results",
                                             api_key = NULL,
                                             species = "9606",
                                             ge_fdr = 0.05,
                                             ge_enrichment_rank_direction = -1,
                                             polling_interval_seconds = 10,
                                             max_polling_attempts = 30) {
  
  # Load required packages
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required but not installed.")
  }
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("Package 'stringr' is required but not installed.")
  }
  
  # Validate input S4 object
  if (!inherits(da_results_for_enrichment, "da_results_for_enrichment")) {
    stop("Input must be an S4 object of class 'da_results_for_enrichment'")
  }
  
  # Get available contrasts
  available_contrasts <- names(da_results_for_enrichment@da_data)
  
  if (length(available_contrasts) == 0) {
    stop("No contrasts found in da_results_for_enrichment@da_data")
  }
  
  # Select contrast
  if (is.null(contrast_name)) {
    contrast_name <- available_contrasts[1]
    message(paste("No contrast specified. Using:", contrast_name))
  } else if (!contrast_name %in% available_contrasts) {
    stop(paste("Contrast", contrast_name, "not found. Available contrasts:", 
               paste(available_contrasts, collapse = ", ")))
  }
  
  # Extract data for the specified contrast
  de_data <- da_results_for_enrichment@da_data[[contrast_name]]
  
  if (is.null(de_data) || nrow(de_data) == 0) {
    stop(paste("No data found for contrast:", contrast_name))
  }
  
  # Validate required columns
  required_cols <- c(identifier_column, "fdr_qvalue", "log2FC")
  missing_cols <- required_cols[!required_cols %in% colnames(de_data)]
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # Validate ranking method
  valid_methods <- c("fdr_qvalue", "log2fc", "combined_score", "none")
  if (!ranking_method %in% valid_methods) {
    stop(paste("Invalid ranking_method. Must be one of:", paste(valid_methods, collapse = ", ")))
  }
  
  # Filter significant proteins if requested
  if (filter_significant) {
    initial_count <- nrow(de_data)
    de_data <- de_data |>
      dplyr::filter(fdr_qvalue < fdr_threshold)
    final_count <- nrow(de_data)
    message(paste("Filtered from", initial_count, "to", final_count, 
                  "proteins using FDR threshold of", fdr_threshold))
    
    if (nrow(de_data) == 0) {
      stop("No proteins remain after filtering for significance")
    }
  }
  
  # Clean protein IDs (remove anything after ":")  
  de_data_processed <- de_data |>
    dplyr::mutate(
      protein_id = purrr::map_chr(!!sym(identifier_column), 
                                  ~stringr::str_split(.x, ":")[[1]][1])
    ) |>
    dplyr::filter(!is.na(protein_id), protein_id != "")
  
  # Apply ranking method and create score column
  if (ranking_method == "fdr_qvalue") {
    # Rank by FDR (ascending - most significant first)
    # Use negative log10 for proper ordering in STRING DB
    de_data_processed <- de_data_processed |>
      dplyr::arrange(fdr_qvalue) |>
      dplyr::mutate(score = -log10(fdr_qvalue + 1e-10))  # Add small value to avoid log(0)
    
  } else if (ranking_method == "log2fc") {
    # Rank by log2FC (descending - highest FC first)
    de_data_processed <- de_data_processed |>
      dplyr::arrange(desc(abs(log2FC))) |>
      dplyr::mutate(score = log2FC)
    
  } else if (ranking_method == "combined_score") {
    # Use sign(log2FC) * (-log10(fdr_qvalue))
    de_data_processed <- de_data_processed |>
      dplyr::mutate(score = sign(log2FC) * (-log10(fdr_qvalue + 1e-10))) |>
      dplyr::arrange(desc(abs(score)))
    
  } else if (ranking_method == "none") {
    # No ranking - use original order with a neutral score
    de_data_processed <- de_data_processed |>
      dplyr::mutate(score = 1)
  }
  
  # Prepare input table for STRING DB
  string_input_table <- de_data_processed |>
    dplyr::select(protein_id, score) |>
    dplyr::filter(!is.na(score), !is.infinite(score))
  
  if (nrow(string_input_table) == 0) {
    stop("No valid protein-score pairs remain after processing")
  }
  
  # Set result label
  if (is.null(result_label)) {
    result_label <- paste0(contrast_name, "_", ranking_method)
  }
  
  message(paste("Submitting", nrow(string_input_table), "proteins to STRING DB"))
  message(paste("Ranking method:", ranking_method))
  message(paste("Result label:", result_label))
  
  # Call the existing MOFA function
  enrichment_results <- runOneStringDbRankEnrichmentMofa(
    input_table = string_input_table,
    identifier_column_name = "protein_id",
    value_column_name = "score",
    result_label = result_label,
    results_dir = results_dir,
    api_key = api_key,
    species = species,
    ge_fdr = ge_fdr,
    ge_enrichment_rank_direction = ge_enrichment_rank_direction,
    polling_interval_seconds = polling_interval_seconds,
    max_polling_attempts = max_polling_attempts
  )
  
  return(enrichment_results)
}


# ----------------------------------------------------------------------------
# runStringDbEnrichmentFromDAResultsMultiple
# ----------------------------------------------------------------------------
#' Run STRING DB Enrichment Analysis for Multiple Contrasts from DE Results
#'
#' @description
#' This function runs STRING DB enrichment analysis for all or selected contrasts
#' from a da_results_for_enrichment S4 object.
#'
#' @param da_results_for_enrichment An S4 object of class da_results_for_enrichment.
#' @param contrast_names Character vector. Names of contrasts to analyze. 
#'   If NULL, analyzes all available contrasts. Default: NULL.
#' @param ranking_method Character string. Same options as runStringDbEnrichmentFromDAResults.
#'   Default: "combined_score".
#' @param ... Additional arguments passed to runStringDbEnrichmentFromDAResults.
#'
#' @return A named list of enrichment results, one for each contrast.
#'
#' @export
runStringDbEnrichmentFromDAResultsMultiple <- function(da_results_for_enrichment,
                                                      contrast_names = NULL,
                                                      ranking_method = "combined_score",
                                                      ...) {
  
  # Get available contrasts
  available_contrasts <- names(da_results_for_enrichment@da_data)
  
  if (is.null(contrast_names)) {
    contrast_names <- available_contrasts
  } else {
    # Validate contrast names
    invalid_contrasts <- contrast_names[!contrast_names %in% available_contrasts]
    if (length(invalid_contrasts) > 0) {
      stop(paste("Invalid contrast names:", paste(invalid_contrasts, collapse = ", ")))
    }
  }
  
  # Run enrichment for each contrast
  results_list <- purrr::map(contrast_names, function(contrast) {
    message(paste("Processing contrast:", contrast))
    
    tryCatch({
      runStringDbEnrichmentFromDAResults(
        da_results_for_enrichment = da_results_for_enrichment,
        contrast_name = contrast,
        ranking_method = ranking_method,
        result_label = paste0(contrast, "_", ranking_method),
        ...
      )
    }, error = function(e) {
      message(paste("Error processing contrast", contrast, ":", e$message))
      return(NULL)
    })
  })
  
  # Name the results list
  names(results_list) <- contrast_names
  
  return(results_list)
}


# ----------------------------------------------------------------------------

