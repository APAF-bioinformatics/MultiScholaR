# printStringDbFunctionalEnrichmentBarGraph
# ----------------------------------------------------------------------------
printStringDbFunctionalEnrichmentBarGraph <- function (input_table, word_limit = 10)
{
  plot_data <- input_table |>
    group_by(comparison, category, termDescription) |>
    arrange( desc( enrichmentScore), falseDiscoveryRate ) |>
    #summarise(enrichmentScore = max(enrichmentScore), falseDiscoveryRate =min(falseDiscoveryRate), direction = first(direction), genesMapped = max(genesMapped)) |>
    dplyr::slice(1) |>
    ungroup() |>
    mutate(termDescriptionAbbrev = sapply(strsplit(as.character(termDescription), " "), function(x) {
      if (length(x) > word_limit) {
        paste(c(head(x, word_limit), "..."), collapse = " ")
      } else {
        paste(x, collapse = " ")
      }
    })) |>
    #mutate(facet_group = interaction(category, comparison)) |>
    # Ensure 'direction' is a factor with all possible levels
    mutate(direction = factor(direction, levels = c("top", "bottom", "both ends")))

  output_group_enrichment_table <- ggplot(plot_data,
                                          aes(y =  reorder_within(  termDescriptionAbbrev, enrichmentScore, list( category )  ) , x = enrichmentScore )) +
    geom_bar(aes(fill=direction),stat = "identity",  width = 0.1) +
    scale_fill_manual(
      name = "Direction",
      values = c("top" = "red", "bottom" = "blue", "both ends" = "grey"),
      drop = FALSE
    ) +
    geom_point(aes(y = reorder_within( termDescriptionAbbrev, enrichmentScore, list( category  )) ,
                   x = enrichmentScore, colour = -log10(falseDiscoveryRate),
                   size = (genesMapped))) + theme(strip.text.y = element_text(angle = 0)) +
    facet_grid( category ~ comparison, scales = "free_y",
                               space = "free") +
    scale_y_reordered() +
    ylab("Term Description") +
    xlab("Enrichment Score")
  output_group_enrichment_table
}


# ----------------------------------------------------------------------------
# getStringDbSpecies
# ----------------------------------------------------------------------------
#' Get Available Species from STRING DB
#'
#' @description
#' This function queries the STRING DB API to retrieve all available species
#' and their corresponding species identifiers. This is useful for finding
#' the correct species code to use in STRING DB enrichment analyses.
#'
#' @param search_term Character string. Optional search term to filter species
#'   by name (case-insensitive partial matching). If NULL, returns all species.
#'   Default: NULL.
#' @param api_key Character string. Your STRING DB API key. Default: NULL.
#'
#' @return A data frame containing available species with columns:
#'   - species_id: STRING DB species identifier (e.g., "9606", "STRG0A62HCE")
#'   - official_name: Official species name
#'   - compact_name: Compact species name
#'   - taxon_id: NCBI taxonomy ID
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Get all available species
#' all_species <- getStringDbSpecies()
#' 
#' # Search for human species
#' human_species <- getStringDbSpecies("sapiens")
#' 
#' # Search for mouse species  
#' mouse_species <- getStringDbSpecies("musculus")
#' 
#' # Search for a specific genus
#' zingiber_species <- getStringDbSpecies("zingiber")
#' }
getStringDbSpecies <- function(search_term = NULL, api_key = NULL) {
  
  # Load required packages
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package 'httr' is required but not installed.")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required but not installed.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required but not installed.")
  }
  
  # STRING DB API endpoint for species
  api_url <- "https://version-12-0.string-db.org/api/json/get_string_ids"
  
  # Prepare parameters
  params <- list(format = "json")
  if (!is.null(api_key)) {
    params$api_key <- api_key
  }
  
  message("Querying STRING DB for available species...")
  
  # Make API request
  response <- tryCatch({
    httr::GET(url = "https://version-12-0.string-db.org/api/json/species", 
              query = params)
  }, error = function(e) {
    stop(paste("Failed to connect to STRING DB API:", e$message))
  })
  
  # Check if request was successful
  if (httr::http_error(response)) {
    stop(paste("STRING DB API request failed with status:", httr::status_code(response)))
  }
  
  # Parse JSON response
  content_text <- httr::content(response, "text", encoding = "UTF-8")
  species_data <- tryCatch({
    jsonlite::fromJSON(content_text, flatten = TRUE)
  }, error = function(e) {
    stop(paste("Failed to parse JSON response from STRING DB:", e$message))
  })
  
  # Convert to data frame if it's not already
  if (!is.data.frame(species_data)) {
    stop("Unexpected response format from STRING DB species API")
  }
  
  # Standardize column names (STRING DB API might have different column names)
  # Try to identify the key columns
  col_mapping <- list()
  
  for (col in colnames(species_data)) {
    lower_col <- tolower(col)
    if (grepl("species.*id|taxon.*id", lower_col)) {
      col_mapping$species_id <- col
    } else if (grepl("official.*name|species.*name", lower_col)) {
      col_mapping$official_name <- col
    } else if (grepl("compact.*name|short.*name", lower_col)) {
      col_mapping$compact_name <- col
    } else if (grepl("^taxon$|ncbi.*id", lower_col)) {
      col_mapping$taxon_id <- col
    }
  }
  
  # Create standardized data frame
  result_df <- species_data
  
  # Rename columns if mappings were found
  for (new_name in names(col_mapping)) {
    old_name <- col_mapping[[new_name]]
    if (old_name %in% colnames(result_df)) {
      result_df <- result_df |>
        dplyr::rename(!!new_name := !!old_name)
    }
  }
  
  # Ensure we have at least species_id and one name column
  if (!"species_id" %in% colnames(result_df)) {
    # Try to find any ID column
    id_cols <- colnames(result_df)[grepl("id", colnames(result_df), ignore.case = TRUE)]
    if (length(id_cols) > 0) {
      result_df <- result_df |>
        dplyr::rename(species_id = !!id_cols[1])
    } else {
      stop("Could not identify species ID column in STRING DB response")
    }
  }
  
  # Apply search filter if provided
  if (!is.null(search_term)) {
    message(paste("Filtering species with search term:", search_term))
    
    # Create a search pattern (case-insensitive)
    search_pattern <- paste0("(?i)", search_term)
    
    # Search across all text columns
    text_cols <- colnames(result_df)[sapply(result_df, is.character)]
    
    # Create a combined search field
    if (length(text_cols) > 0) {
      result_df <- result_df |>
        dplyr::mutate(
          search_field = paste(!!!syms(text_cols), sep = " ", collapse = " ")
        ) |>
        dplyr::filter(grepl(search_pattern, search_field)) |>
        dplyr::select(-search_field)
    }
  }
  
  # Sort by species_id for consistent output
  result_df <- result_df |>
    dplyr::arrange(species_id)
  
  message(paste("Found", nrow(result_df), "species matching your criteria"))
  
  return(result_df)
}


# ----------------------------------------------------------------------------
# searchStringDbSpecies
# ----------------------------------------------------------------------------
#' Search STRING DB Species by Name
#'
#' @description
#' A convenience function to search for STRING DB species by name.
#' This is a wrapper around getStringDbSpecies with better error handling
#' and formatted output for common use cases.
#'
#' @param species_name Character string. The species name to search for
#'   (case-insensitive partial matching).
#' @param api_key Character string. Your STRING DB API key. Default: NULL.
#' @param show_top_n Integer. Maximum number of results to display. Default: 10.
#'
#' @return A data frame with the most relevant species matches.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Search for human
#' searchStringDbSpecies("human")
#' searchStringDbSpecies("homo sapiens")
#' 
#' # Search for mouse
#' searchStringDbSpecies("mouse")
#' searchStringDbSpecies("mus musculus")
#' 
#' # Search for your organism
#' searchStringDbSpecies("zingiber")
#' searchStringDbSpecies("ginger")
#' }
searchStringDbSpecies <- function(species_name, api_key = NULL, show_top_n = 10) {
  
  if (missing(species_name) || is.null(species_name) || species_name == "") {
    stop("Please provide a species name to search for")
  }
  
  # Get species data
  species_df <- getStringDbSpecies(search_term = species_name, api_key = api_key)
  
  if (nrow(species_df) == 0) {
    message(paste("No species found matching:", species_name))
    message("Try using a more general search term or check the spelling")
    return(tibble::tibble())
  }
  
  # Limit results
  if (nrow(species_df) > show_top_n) {
    message(paste("Showing top", show_top_n, "results out of", nrow(species_df), "total matches"))
    species_df <- species_df |> dplyr::slice_head(n = show_top_n)
  }
  
  # Display helpful information
  message("\nSpecies search results:")
  message("Use the 'species_id' column value as the 'species' parameter in STRING DB functions")
  
  return(species_df)
}


# ----------------------------------------------------------------------------
# runStringDbEnrichmentAllContrasts
# ----------------------------------------------------------------------------
#' Run STRING DB Enrichment for All Contrasts
#'
#' @description
#' Executes STRING DB enrichment analysis for all contrasts in a differential expression
#' analysis results list. Results are cached to avoid redundant API calls. If a cached
#' result file exists, it will be loaded instead of re-running the analysis.
#'
#' @param da_analysis_results_list A named list containing differential expression results.
#'   Each element should contain a `da_proteins_long` data frame with columns:
#'   `log2FC`, `fdr_qvalue`, and `Protein.Ids`.
#' @param project_dirs A list containing project directory structure. Must include
#'   an element named `paste0(omic_type, "_", experiment_label)` with a `results_dir` component.
#' @param omic_type Character string: The omics type (e.g., "proteomics", "transcriptomics").
#' @param experiment_label Character string: The experiment label for directory construction.
#' @param api_key Character string: Your personal STRING API key. Default is "bjcR4Px5rByQ".
#' @param species Character or numeric: NCBI/STRING species identifier. Default is "9606" (Homo sapiens).
#'   This parameter can be set to the taxon_id variable defined in your environment.
#' @param ge_fdr Numeric: FDR threshold for gene expression enrichment. Default is 0.05.
#' @param ge_enrichment_rank_direction Integer: Direction for enrichment rank
#'   (-1, 0, or 1). Default is -1.
#' @param polling_interval_seconds Numeric: Seconds to wait between polling attempts.
#'   Default is 10.
#' @param max_polling_attempts Numeric: Maximum polling attempts. Default is 30.
#' @param force_refresh Logical: If TRUE, forces re-running the analysis even if
#'   cached results exist. Default is FALSE.
#' @param comparison_name_transform Function: Optional function to transform contrast names
#'   for display in the comparison column. Default is NULL (uses simple transformation).
#'
#' @return A data frame containing enrichment results for all contrasts with columns:
#'   - `comparison`: The contrast name (transformed)
#'   - All columns returned by STRING DB enrichment analysis
#'
#' @details
#' This function:
#' 1. Checks for cached results in `{pathway_dir}/string_db/all_enrichment_results.rds`
#' 2. If cache exists and `force_refresh = FALSE`, loads and returns cached results
#' 3. If cache doesn't exist or `force_refresh = TRUE`, runs enrichment for all contrasts
#' 4. Saves results to both RDS (for caching) and TSV (for inspection) formats
#'
#' Results are saved to the pathway enrichment directory structure defined in `setupDirectories()`.
#' The directory structure is automatically created if it doesn't exist.
#'
#' @importFrom dplyr bind_rows
#' @importFrom purrr map map_chr
#' @importFrom stringr str_split str_split_i
#' @importFrom vroom vroom_write
#'
#' @examples
#' \dontrun{
#' # Run enrichment for all contrasts
#' all_enrichment_results <- runStringDbEnrichmentAllContrasts(
#'   da_analysis_results_list = my_da_results,
#'   project_dirs = project_dirs,
#'   omic_type = "proteomics",
#'   experiment_label = "standard",
#'   api_key = "YOUR_API_KEY",
#'   species = taxon_id  # Use taxon_id from your environment
#' )
#'
#' # Force refresh of cached results
#' refreshed_results <- runStringDbEnrichmentAllContrasts(
#'   da_analysis_results_list = my_da_results,
#'   project_dirs = project_dirs,
#'   omic_type = "proteomics",
#'   experiment_label = "standard",
#'   api_key = "YOUR_API_KEY",
#'   species = taxon_id,
#'   force_refresh = TRUE
#' )
#' }
#'
#' @export
runStringDbEnrichmentAllContrasts <- function(da_analysis_results_list = NULL,
                                             project_dirs,
                                             omic_type,
                                             experiment_label,
                                             api_key = "bjcR4Px5rByQ",
                                             species = "9606",
                                             ge_fdr = 0.05,
                                             ge_enrichment_rank_direction = -1,
                                             polling_interval_seconds = 10,
                                             max_polling_attempts = 30,
                                             force_refresh = FALSE,
                                             comparison_name_transform = NULL,
                                             de_analysis_results_list = NULL) {
  
  # Handle alias for backward compatibility
  if (is.null(da_analysis_results_list)) {
    if (!is.null(de_analysis_results_list)) {
      da_analysis_results_list <- de_analysis_results_list
    } else {
      stop("Either da_analysis_results_list or de_analysis_results_list must be provided")
    }
  }
  
  # Load required packages
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required but not installed.")
  }
  if (!requireNamespace("purrr", quietly = TRUE)) {
    stop("Package 'purrr' is required but not installed.")
  }
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("Package 'stringr' is required but not installed.")
  }
  if (!requireNamespace("vroom", quietly = TRUE)) {
    stop("Package 'vroom' is required but not installed.")
  }
  
  # Validate inputs
  if (!is.list(da_analysis_results_list) || length(da_analysis_results_list) == 0) {
    stop("da_analysis_results_list must be a non-empty list")
  }
  if (!is.list(project_dirs)) {
    stop("project_dirs must be a list")
  }
  if (is.null(api_key) || api_key == "") {
    stop("api_key is required for STRING DB enrichment")
  }
  
  # Construct directory paths
  dir_key <- paste0(omic_type, "_", experiment_label)
  if (!dir_key %in% names(project_dirs)) {
    stop(paste("project_dirs does not contain element:", dir_key))
  }
  
  # Use the pathway directory defined in project structure
  pathway_dir <- project_dirs[[dir_key]]$pathway_dir
  enrichment_dir <- file.path(pathway_dir, "string_db")
  enrichment_cache_file <- file.path(enrichment_dir, "all_enrichment_results.rds")
  enrichment_tsv_file <- file.path(enrichment_dir, "all_enrichment_results.tsv")
  
  # Check if cache exists and should be used
  if (file.exists(enrichment_cache_file) && !force_refresh) {
    message("Loading cached enrichment results from: ", enrichment_cache_file)
    all_enrichment_results <- readRDS(enrichment_cache_file)
    message("Loaded ", nrow(all_enrichment_results), " rows from cached enrichment results.")
    return(all_enrichment_results)
  }
  
  # Run enrichment for all contrasts
  message("Running STRING DB enrichment for all contrasts...")
  message("Number of contrasts to process: ", length(da_analysis_results_list))
  
  list_of_contrasts <- names(da_analysis_results_list)
  
  output_group_tables_list <- purrr::map(list_of_contrasts, function(contrast_name) {
    message("Processing contrast: ", contrast_name)
    
    input_table <- da_analysis_results_list[[contrast_name]]$da_proteins_long
    if (is.null(input_table)) {
      input_table <- da_analysis_results_list[[contrast_name]]$de_proteins_long
    }
    
    if (is.null(input_table)) {
      message("Warning: No DA/DE results table found for contrast: ", contrast_name)
      return(NULL)
    }
    
    result_label <- stringr::str_split_i(contrast_name, "=", 1)
    
    # Run enrichment for this contrast
    tryCatch({
      runOneStringDbRankEnrichment(
        input_table = input_table,
        result_label = result_label,
        pathway_dir = pathway_dir,
        api_key = api_key,
        species = species,
        ge_fdr = ge_fdr,
        ge_enrichment_rank_direction = ge_enrichment_rank_direction,
        polling_interval_seconds = polling_interval_seconds,
        max_polling_attempts = max_polling_attempts
      )
    }, error = function(e) {
      message("Error processing contrast ", contrast_name, ": ", e$message)
      return(NULL)
    })
  })
  
  # Name the list elements
  if (is.null(comparison_name_transform)) {
    # Default transformation: extract first part before underscore and add "(M-C)"
    names(output_group_tables_list) <- purrr::map_chr(list_of_contrasts, function(x) {
      paste(stringr::str_split(x, "_")[[1]][1], "(M-C)")
    })
  } else {
    names(output_group_tables_list) <- purrr::map_chr(list_of_contrasts, comparison_name_transform)
  }
  
  # Combine all results
  output_group_main_table <- dplyr::bind_rows(output_group_tables_list, .id = "comparison")
  
  # Create directory if it doesn't exist
  dir.create(enrichment_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Save results
  message("Saving enrichment results to: ", enrichment_cache_file)
  saveRDS(output_group_main_table, file = enrichment_cache_file)
  vroom::vroom_write(output_group_main_table, enrichment_tsv_file)
  
  message("Enrichment analysis complete. Processed ", nrow(output_group_main_table), " total results.")
  
  return(output_group_main_table)
}


# ----------------------------------------------------------------------------
# plotStringDbEnrichmentResults
# ----------------------------------------------------------------------------
#' Plot STRING DB Functional Enrichment Results
#'
#' @description
#' Loads cached STRING DB enrichment results, filters to top terms per category,
#' and generates publication-quality plots organized by functional category groups.
#' Creates separate plots for Gene Ontology terms, pathways, protein domains, and other categories.
#'
#' @param project_dirs A list containing project directory structure. Must include
#'   an element named `paste0(omic_type, "_", experiment_label)` with a `results_dir` component.
#' @param omic_type Character string: The omics type (e.g., "proteomics", "transcriptomics").
#' @param experiment_label Character string: The experiment label for directory construction.
#' @param enrichment_results Data frame: Optional pre-loaded enrichment results.
#'   If NULL (default), will attempt to load from cache file.
#' @param top_n_terms Integer: Number of top terms to include per category and comparison.
#'   Default is 5.
#' @param word_limit Integer: Maximum number of words in term descriptions before truncation.
#'   Default is 7.
#' @param plot_width Numeric: Width of output plots in inches. Default is 16.
#' @param plot_height Numeric: Height of output plots in inches. Default is 12.
#' @param plot_dpi Numeric: DPI resolution for PNG output. Default is 300.
#' @param save_plots Logical: Whether to save plots to files. Default is TRUE.
#' @param return_plots Logical: Whether to return plot objects in a list. Default is TRUE.
#' @param print_plots Logical: Whether to print plots to graphics device. Default is TRUE.
#'
#' @return A list containing:
#'   - `filtered_results`: The filtered enrichment results data frame
#'   - `plots`: A named list of ggplot objects (if `return_plots = TRUE`):
#'     - `go_terms`: Plot of Gene Ontology terms
#'     - `pathways`: Plot of pathway enrichment (KEGG, Reactome, WikiPathways, STRING clusters)
#'     - `protein_domains`: Plot of protein domain enrichment (InterPro, Pfam, SMART)
#'     - `others`: Plot of other enrichment categories
#'   Returns NULL if cached results file is not found.
#'
#' @details
#' This function:
#' 1. Loads cached enrichment results from `{pathway_dir}/string_db/all_enrichment_results.rds`
#' 2. Filters to the top N terms per category and comparison based on enrichmentScore and falseDiscoveryRate
#' 3. Saves filtered results to TSV file
#' 4. Generates four separate plots organized by category type:
#'    - **GO Terms**: Gene Ontology (Biological Process, Molecular Function, Cellular Component)
#'    - **Pathways**: KEGG, Reactome, WikiPathways, STRING clusters
#'    - **Protein Domains**: InterPro, Pfam, SMART
#'    - **Others**: All remaining categories
#' 5. Saves plots as both PNG and PDF formats
#'
#' Results are loaded from and saved to the pathway enrichment directory structure defined in `setupDirectories()`.
#' All plots use the `printStringDbFunctionalEnrichmentBarGraph` function for consistent
#' visualization with faceting by category and comparison.
#'
#' @importFrom dplyr filter inner_join arrange group_by ungroup mutate distinct
#' @importFrom stringr str_detect
#' @importFrom vroom vroom_write
#' @importFrom ggplot2 ggsave
#' @importFrom purrr walk
#'
#' @examples
#' \dontrun{
#' # Plot enrichment results with default settings
#' plot_results <- plotStringDbEnrichmentResults(
#'   project_dirs = project_dirs,
#'   omic_type = "proteomics",
#'   experiment_label = "standard"
#' )
#'
#' # Access individual plots
#' print(plot_results$plots$go_terms)
#' print(plot_results$plots$pathways)
#'
#' # Use custom filtering and plotting parameters
#' custom_plots <- plotStringDbEnrichmentResults(
#'   project_dirs = project_dirs,
#'   omic_type = "proteomics",
#'   experiment_label = "standard",
#'   top_n_terms = 10,
#'   word_limit = 10,
#'   plot_width = 20,
#'   plot_height = 14,
#'   save_plots = FALSE,
#'   print_plots = FALSE
#' )
#'
#' # Use pre-loaded enrichment results
#' my_results <- readRDS("path/to/enrichment_results.rds")
#' plots <- plotStringDbEnrichmentResults(
#'   enrichment_results = my_results,
#'   project_dirs = project_dirs,
#'   omic_type = "proteomics",
#'   experiment_label = "standard"
#' )
#' }
#'
#' @export
plotStringDbEnrichmentResults <- function(project_dirs,
                                         omic_type,
                                         experiment_label,
                                         enrichment_results = NULL,
                                         top_n_terms = 5,
                                         word_limit = 7,
                                         plot_width = 16,
                                         plot_height = 12,
                                         plot_dpi = 300,
                                         fdr_threshold = 0.05,
                                         save_plots = TRUE,
                                         return_plots = TRUE,
                                         print_plots = TRUE) {
  
  # Load required packages
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required but not installed.")
  }
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("Package 'stringr' is required but not installed.")
  }
  if (!requireNamespace("vroom", quietly = TRUE)) {
    stop("Package 'vroom' is required but not installed.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but not installed.")
  }
  if (!requireNamespace("purrr", quietly = TRUE)) {
    stop("Package 'purrr' is required but not installed.")
  }
  
  # Validate inputs
  if (!is.list(project_dirs)) {
    stop("project_dirs must be a list")
  }
  
  # Construct directory paths
  dir_key <- paste0(omic_type, "_", experiment_label)
  if (!dir_key %in% names(project_dirs)) {
    stop(paste("project_dirs does not contain element:", dir_key))
  }
  
  # Use the pathway directory defined in project structure
  pathway_dir <- project_dirs[[dir_key]]$pathway_dir
  enrichment_dir <- file.path(pathway_dir, "string_db")
  enrichment_cache_file <- file.path(enrichment_dir, "all_enrichment_results.rds")
  
  # Load enrichment results if not provided
  if (is.null(enrichment_results)) {
    if (!file.exists(enrichment_cache_file)) {
      message("Enrichment results file not found at: ", enrichment_cache_file)
      message("Please run runStringDbEnrichmentAllContrasts() first.")
      return(NULL)
    }
    
    message("Loading enrichment results from: ", enrichment_cache_file)
    all_enrichment_results <- readRDS(enrichment_cache_file)
    message("Loaded ", nrow(all_enrichment_results), " rows from cached enrichment results.")
  } else {
    all_enrichment_results <- enrichment_results
    message("Using provided enrichment results with ", nrow(all_enrichment_results), " rows.")
  }
  
  # Check for empty results
  cat(sprintf("DEBUG_STR: Checking results. Rows: %d, Cols: %d\n", nrow(all_enrichment_results), ncol(all_enrichment_results)))
  cat(sprintf("DEBUG_STR: Columns: %s\n", paste(colnames(all_enrichment_results), collapse = ", ")))

  if (nrow(all_enrichment_results) == 0) {
    message("No enrichment results found (0 rows). Skipping plotting.")
    return(NULL)
  }

  # Check for required columns
  required_cols <- c("comparison", "category", "enrichmentScore", "falseDiscoveryRate", "termID")
  missing_cols <- setdiff(required_cols, colnames(all_enrichment_results))
  if (length(missing_cols) > 0) {
    warning("Enrichment results missing required columns: ", paste(missing_cols, collapse = ", "))
    return(NULL)
  }

  # Filter to top N terms per category and comparison
  message("Filtering to top ", top_n_terms, " terms per category and comparison...")
  
  sorted_enrichment_results <- all_enrichment_results |>
    dplyr::group_by(comparison, category) |>
    dplyr::arrange(comparison, category, desc(enrichmentScore), falseDiscoveryRate) |>
    dplyr::mutate(rank = dplyr::row_number()) |>
    dplyr::ungroup()
  
  included_functional_category_and_term <- sorted_enrichment_results |>
    dplyr::filter(rank <= top_n_terms & falseDiscoveryRate < fdr_threshold) |>
    dplyr::distinct(category, termID)
  
  filtered_enrichment_results <- sorted_enrichment_results |>
    dplyr::inner_join(included_functional_category_and_term, by = c("category", "termID")) |>
    dplyr::arrange(comparison, category, desc(enrichmentScore), falseDiscoveryRate)
  
  message("Filtered results contain ", nrow(filtered_enrichment_results), " rows.")
  
  # Save filtered results
  filtered_results_file <- file.path(enrichment_dir, "filtered_enrichment_results.tsv")
  message("Saving filtered results to: ", filtered_results_file)
  vroom::vroom_write(filtered_enrichment_results, filtered_results_file)
  
  # Define category groups for splitting plots
  fig1_pattern <- "GO Component|GO Function|GO Process|Gene Ontology"
  fig2_categories <- c("STRING clusters", "Reactome", "KEGG", "WikiPathways")
  fig3_categories <- c("InterPro", "Pfam", "SMART")
  
  # Filter data for each figure
  message("Organizing results by category groups...")
  
  results_fig1 <- filtered_enrichment_results |>
    dplyr::filter(stringr::str_detect(category, fig1_pattern))
  
  results_fig2 <- filtered_enrichment_results |>
    dplyr::filter(category %in% fig2_categories)
  
  results_fig3 <- filtered_enrichment_results |>
    dplyr::filter(category %in% fig3_categories)
  
  results_fig4 <- filtered_enrichment_results |>
    dplyr::filter(!stringr::str_detect(category, fig1_pattern) &
                    !category %in% fig2_categories &
                    !category %in% fig3_categories)
  
  message("Category distribution:")
  message("  GO terms: ", nrow(results_fig1), " rows")
  message("  Pathways: ", nrow(results_fig2), " rows")
  message("  Protein domains: ", nrow(results_fig3), " rows")
  message("  Others: ", nrow(results_fig4), " rows")
  
  # Helper function to generate, save, and optionally print each plot
  plotAndSave <- function(data, filename_suffix, word_limit_val, width, height, dpi) {
    if (nrow(data) == 0) {
      message("No enrichment results to plot for ", filename_suffix)
      return(NULL)
    }
    
    message("Generating enrichment plot for ", filename_suffix, "...")
    
    enrichment_plot <- printStringDbFunctionalEnrichmentBarGraph(data, word_limit = word_limit_val)
    
    if (save_plots) {
      output_filename <- file.path(enrichment_dir, paste0("enrichment_summary_", filename_suffix))
      message("Saving plot to: ", output_filename, ".{png,pdf}")
      
      purrr::walk(c(".png", ".pdf"), function(ext) {
        ggplot2::ggsave(
          plot = enrichment_plot,
          filename = paste0(output_filename, ext),
          width = width,
          height = height,
          dpi = dpi
        )
      })
    }
    
    if (print_plots) {
      print(enrichment_plot)
    }
    
    return(enrichment_plot)
  }
  
  # Generate all plots
  message("\n=== Generating Plots ===")
  
  plot_go_terms <- plotAndSave(
    data = results_fig1,
    filename_suffix = "go_terms",
    word_limit_val = word_limit,
    width = plot_width,
    height = plot_height,
    dpi = plot_dpi
  )
  
  plot_pathways <- plotAndSave(
    data = results_fig2,
    filename_suffix = "pathways",
    word_limit_val = word_limit,
    width = plot_width,
    height = plot_height,
    dpi = plot_dpi
  )
  
  plot_protein_domains <- plotAndSave(
    data = results_fig3,
    filename_suffix = "protein_domains",
    word_limit_val = word_limit,
    width = plot_width,
    height = plot_height,
    dpi = plot_dpi
  )
  
  plot_others <- plotAndSave(
    data = results_fig4,
    filename_suffix = "others",
    word_limit_val = word_limit,
    width = plot_width,
    height = plot_height,
    dpi = plot_dpi
  )
  
  message("\n=== Plotting Complete ===")
  
  # Prepare return list
  result_list <- list(
    filtered_results = filtered_enrichment_results
  )
  
  if (return_plots) {
    result_list$plots <- list(
      go_terms = plot_go_terms,
      pathways = plot_pathways,
      protein_domains = plot_protein_domains,
      others = plot_others
    )
  }
  
  return(result_list)
}


# ----------------------------------------------------------------------------

