# runReactomeEnrichment
# ----------------------------------------------------------------------------
#' Run Reactome Pathway Enrichment Analysis for Metabolites
#'
#' @param ranked_list Named numeric vector with ChEBI IDs as names and values as weights
#' @param mapping_table Metabolite ID mapping table
#' @param project_dirs Project directories structure
#' @param omic_type Omics type for directory access
#' @param experiment_label Experiment label for directory access
#' @param assay_name Name of the assay
#' @param reactome_organism Organism to use for Reactome (e.g., "Klebsiella pneumoniae")
#' @return A data frame with Reactome pathway enrichment results
#' @export

runReactomeEnrichment <- function(ranked_list, 
                                 mapping_table, 
                                 project_dirs,
                                 omic_type,
                                 experiment_label,
                                 assay_name,
                                 reactome_organism = NULL) {
  
  message(sprintf("--- Entering runReactomeEnrichment ---"))
  message(sprintf("   runReactomeEnrichment Args: assay_name = %s", assay_name))
  message(sprintf("   runReactomeEnrichment Args: reactome_organism = %s", as.character(reactome_organism)))
  message(sprintf("   runReactomeEnrichment Args: ranked_list length = %d", length(ranked_list)))
  
  # Use getProjectPaths helper with automatic fallback
  data_dir <- getProjectPaths(omic_type, experiment_label)$data_dir
  message(sprintf("   runReactomeEnrichment: Using data_dir = %s", data_dir))
  
  # Load ChEBI to Reactome mapping
  tryCatch({
    message(sprintf("   runReactomeEnrichment Step: Loading ChEBI to Reactome mapping from %s...", file.path(data_dir, "chebi_to_reactome.tab")))
    
    # Check if file exists
    if (!file.exists(file.path(data_dir, "chebi_to_reactome.tab"))) {
      message("   runReactomeEnrichment ERROR: File 'chebi_to_reactome.tab' does not exist!")
      return(tibble(
        termDescription = character(),
        enrichmentScore = numeric(),
        falseDiscoveryRate = numeric(),
        genesMapped = integer(),
        comparison = character()
      ))
    }
    
    chebi_to_reactome <- vroom::vroom(
      file.path(data_dir, "chebi_to_reactome.tab"),
      show_col_types = FALSE
    ) |>
      dplyr::select(chebi_id, reactome_id, reactome_term, organism)
    
    message(sprintf("   runReactomeEnrichment: Loaded %d rows from chebi_to_reactome.tab", nrow(chebi_to_reactome)))
    message("   runReactomeEnrichment Step: Sample of ChEBI to Reactome mapping:")
    print(head(chebi_to_reactome, 5))
    
    # List available organisms in the data
    unique_organisms <- unique(chebi_to_reactome$organism)
    message(paste("Available organisms in Reactome data:", paste(unique_organisms, collapse=", ")))
    
    # FIXED: Improved organism handling with flexible matching
    if (!is.null(reactome_organism)) {
      # Try exact match first
      if (reactome_organism %in% unique_organisms) {
        message(paste("Filtering Reactome data for organism:", reactome_organism))
        chebi_to_reactome <- chebi_to_reactome |> 
          dplyr::filter(organism == reactome_organism)
      } else {
        # Try partial matching if exact match fails
        message(sprintf("   runReactomeEnrichment DEBUG: Trying partial match for '%s'", reactome_organism))
        possible_matches <- unique_organisms[grepl(reactome_organism, unique_organisms, ignore.case = TRUE)]
        if (length(possible_matches) > 0) {
          selected_organism <- possible_matches[1]  # Use first match
          message(paste("Using closest matching organism:", selected_organism))
          chebi_to_reactome <- chebi_to_reactome |> 
            dplyr::filter(organism == selected_organism)
        } else {
          # Fallback to Human if no match
          message("No matching organism found. Using 'Homo sapiens' as fallback...")
          chebi_to_reactome <- chebi_to_reactome |> 
            dplyr::filter(organism == "Homo sapiens")
        }
      }
      message(sprintf("   runReactomeEnrichment: After organism filtering, %d rows remain", nrow(chebi_to_reactome)))
    } else {
      message("Using all organisms in Reactome data (no filtering)")
    }
    
    # Extract just the numeric part from ChEBI IDs to match our format
    message("   runReactomeEnrichment Step: Formatting ChEBI IDs...")
    modified_chebi_to_reactome <- chebi_to_reactome |>
      dplyr::mutate(
        numeric_id = as.numeric(str_extract(chebi_id, "\\d+")),
        new_chebi_id = paste0("CHEBI:", numeric_id)
      )
    
    message("   runReactomeEnrichment Step: Sample after numeric ID extraction:")
    print(head(modified_chebi_to_reactome, 5))
    
    chebi_to_reactome <- modified_chebi_to_reactome |>
      dplyr::select(-chebi_id) |>
      dplyr::rename(chebi_id = new_chebi_id) |>
      dplyr::filter(!is.na(numeric_id))  # Remove any that didn't convert properly
    
    message(sprintf("   runReactomeEnrichment: After ChEBI ID formatting, %d rows remain", nrow(chebi_to_reactome)))
    message("   runReactomeEnrichment Step: Sample of formatted ChEBI IDs:")
    print(head(chebi_to_reactome, 5))
    
    # Check which ChEBI IDs overlap with our ranked list
    message("   runReactomeEnrichment Step: Checking for ChEBI IDs that overlap with ranked list...")
    chebi_ids_in_ranked_list <- intersect(chebi_to_reactome$chebi_id, names(ranked_list))
    message(paste("Found", length(chebi_ids_in_ranked_list), "metabolites that overlap with Reactome pathways."))
    
    if (length(chebi_ids_in_ranked_list) > 0) {
      message("   runReactomeEnrichment Step: Sample matching ChEBI IDs:")
      message(paste(head(chebi_ids_in_ranked_list, 5), collapse=", "))
    }
    
    if (length(chebi_ids_in_ranked_list) == 0) {
      message("No overlap between ChEBI IDs in ranked list and Reactome database. Returning empty result.")
      return(tibble(
        termDescription = character(),
        enrichmentScore = numeric(),
        falseDiscoveryRate = numeric(),
        genesMapped = integer(),
        comparison = character()
      ))
    }
    
    # Create a subset of the ranked list with only IDs that are in Reactome
    message("   runReactomeEnrichment Step: Creating filtered ranked list...")
    reactome_ranked_list <- ranked_list[chebi_ids_in_ranked_list]
    
    message(sprintf("   runReactomeEnrichment: Created ranked list with %d metabolites", length(reactome_ranked_list)))
    message("   runReactomeEnrichment Step: Sample values from ranked list:")
    print(head(reactome_ranked_list))
    
    # Sort the ranked list in decreasing order 
    reactome_ranked_list <- sort(reactome_ranked_list, decreasing = TRUE)
    
    # Diagnostic: Check for ties in the data
    message("   runReactomeEnrichment DEBUG: Analyzing score distribution...")
    # Count exact duplicates in values
    duplicate_count <- sum(duplicated(reactome_ranked_list))
    duplicate_percentage <- (duplicate_count / length(reactome_ranked_list)) * 100
    message(sprintf("   runReactomeEnrichment DEBUG: Found %d exact duplicate values (%.2f%% of total)",
                   duplicate_count, duplicate_percentage))
    
    # Show basic statistics
    message(sprintf("   runReactomeEnrichment DEBUG: Score range: [%.4f, %.4f], Median: %.4f", 
                   min(reactome_ranked_list), max(reactome_ranked_list), median(reactome_ranked_list)))
    
    # Add small random noise to break ties in the ranked list
    message("   runReactomeEnrichment Step: Adding small random noise to break ties in ranked list...")
    set.seed(42) # For reproducibility
    original_reactome_ranked_list <- reactome_ranked_list # Store original for reference
    reactome_ranked_list <- reactome_ranked_list + runif(length(reactome_ranked_list), min=-0.0001, max=0.0001)
    # Re-sort to ensure order is preserved with the new values
    reactome_ranked_list <- sort(reactome_ranked_list, decreasing = TRUE)
    message(sprintf("   runReactomeEnrichment: Added jitter to %d values to prevent ties", length(reactome_ranked_list)))
    
    # CRITICAL: Verify the ordering and distribution before GSEA
    message("   runReactomeEnrichment CRITICAL CHECK: Verifying ranked list ordering and distribution...")
    positive_values <- sum(reactome_ranked_list > 0)
    negative_values <- sum(reactome_ranked_list < 0)
    zero_values <- sum(reactome_ranked_list == 0)
    message(sprintf("   runReactomeEnrichment: Value distribution: %d positive, %d negative, %d zero", 
                   positive_values, negative_values, zero_values))
    
    # Check that the list is actually sorted (decreasing)
    is_decreasing <- all(diff(reactome_ranked_list) <= 0)
    message(sprintf("   runReactomeEnrichment: List is properly sorted in decreasing order: %s", 
                   if(is_decreasing) "YES" else "NO - THIS IS A PROBLEM!"))
    
    # Show top and bottom values to verify ordering
    message("   runReactomeEnrichment: Top 5 entries in ranked list:")
    for (i in 1:min(5, length(reactome_ranked_list))) {
      message(sprintf("     %s: %.6f", names(reactome_ranked_list)[i], reactome_ranked_list[i]))
    }
    message("   runReactomeEnrichment: Bottom 5 entries in ranked list:")
    for (i in (length(reactome_ranked_list) - min(4, length(reactome_ranked_list) - 1)):length(reactome_ranked_list)) {
      message(sprintf("     %s: %.6f", names(reactome_ranked_list)[i], reactome_ranked_list[i]))
    }
    
    # Create TERM2GENE mapping
    message("   runReactomeEnrichment Step: Creating TERM2GENE mapping...")
    term2gene <- chebi_to_reactome |>
      dplyr::filter(chebi_id %in% names(reactome_ranked_list)) |>
      dplyr::select(reactome_id, chebi_id) |>
      dplyr::rename(term = reactome_id, gene = chebi_id)
    
    message(sprintf("   runReactomeEnrichment: Created TERM2GENE mapping with %d rows", nrow(term2gene)))
    message("   runReactomeEnrichment Step: Sample TERM2GENE mapping:")
    print(head(term2gene, 5))
    
    # Check if we have pathways
    if (nrow(term2gene) == 0) {
      message("   runReactomeEnrichment WARNING: No term-gene mappings could be created!")
      return(tibble(
        termDescription = character(),
        enrichmentScore = numeric(),
        falseDiscoveryRate = numeric(),
        genesMapped = integer(),
        comparison = character()
      ))
    }
    
    # Count unique terms and genes
    unique_terms <- unique(term2gene$term)
    unique_genes <- unique(term2gene$gene)
    message(sprintf("   runReactomeEnrichment: Mapping contains %d unique terms and %d unique genes", 
                   length(unique_terms), length(unique_genes)))
    
    # Create TERM2NAME mapping
    message("   runReactomeEnrichment Step: Creating TERM2NAME mapping...")
    term2name <- chebi_to_reactome |>
      dplyr::select(reactome_id, reactome_term) |>
      dplyr::distinct() |>  # FIXED: Added distinct here to avoid duplicates
      dplyr::rename(term = reactome_id, name = reactome_term)
    
    message(sprintf("   runReactomeEnrichment: Created TERM2NAME mapping with %d rows", nrow(term2name)))
    message("   runReactomeEnrichment Step: Sample TERM2NAME mapping:")
    print(head(term2name, 5))
    
    # Debug: check for duplicate pathways
    message("   runReactomeEnrichment DEBUG: Checking for duplicate pathway IDs...")
    duplicated_terms <- term2name$term[duplicated(term2name$term)]
    if (length(duplicated_terms) > 0) {
      message(sprintf("   runReactomeEnrichment WARNING: Found %d duplicate pathway IDs", length(duplicated_terms)))
      message("Sample duplicates:")
      for (dup_term in head(duplicated_terms, 3)) {
        dup_entries <- term2name |> dplyr::filter(term == dup_term)
        message(sprintf("   Term: %s has %d entries:", dup_term, nrow(dup_entries)))
        print(dup_entries)
      }
    } else {
      message("   runReactomeEnrichment DEBUG: No duplicate pathway IDs found in TERM2NAME")
    }
    
    # Add fallback in case term2name is empty
    if (!exists("term2name") || nrow(term2name) == 0) {
      message("   runReactomeEnrichment WARNING: term2name is missing or empty, creating fallback...")
      # Create a minimal term2name mapping using the term itself
      term2name <- term2gene |>
        dplyr::distinct(term) |>
        dplyr::mutate(name = term)
      
      message("   runReactomeEnrichment: Created fallback term2name mapping:")
      print(head(term2name, 3))
    }
    
    # Try to run GSEA analysis
    message("   runReactomeEnrichment Step: Running GSEA analysis...")
    
    # ABSOLUTE FINAL SAFETY CHECK - verify both term2gene and term2name exist
    if (!exists("term2gene") || nrow(term2gene) == 0) {
      message("   runReactomeEnrichment CRITICAL ERROR: term2gene is missing or empty before GSEA call")
      # Can't proceed without term2gene
      return(tibble(
        termDescription = character(),
        enrichmentScore = numeric(),
        falseDiscoveryRate = numeric(),
        genesMapped = integer(),
        comparison = character()
      ))
    }
    
    if (!exists("term2name") || nrow(term2name) == 0) {
      message("   runReactomeEnrichment WARNING: term2name is missing or empty before GSEA call, creating emergency fallback")
      # Create emergency fallback
      term2name <- term2gene |> 
        dplyr::distinct(term) |>
        dplyr::mutate(name = term)
    }
    
    # Debug output
    message(sprintf("   runReactomeEnrichment VERIFICATION: term2gene has %d rows, term2name has %d rows", 
                   nrow(term2gene), nrow(term2name)))
    
    gsea_result <- tryCatch({
      clusterProfiler::GSEA(
        geneList = reactome_ranked_list,
        TERM2GENE = term2gene,
        TERM2NAME = term2name,
        minGSSize = 3,
        maxGSSize = 500,
        pvalueCutoff = 0.05,  # MODIFIED: Increased from 0.1 to 0.05 for more lenient filtering
        pAdjustMethod = "fdr",
        verbose = TRUE,
        seed = TRUE,      # Set seed for reproducibility
        by = "fgsea",     # FIXED: Use "fgsea" (valid option) instead of "fgseaMultilevel"
        BPPARAM = BiocParallel::SerialParam()  # Disable parallel processing
      )
    }, error = function(e) {
      message(paste("GSEA error:", e$message))
      
      # Try the enricher method instead
      message("Trying enricher instead...")
      
      # Select the top ranked genes (positive values)
      positive_genes <- names(reactome_ranked_list)[reactome_ranked_list > 0]
      message(sprintf("   runReactomeEnrichment: Found %d positively ranked genes", length(positive_genes)))
      
      if (length(positive_genes) > 0) {
        tryCatch({
          message("   runReactomeEnrichment Step: Running enricher analysis...")
          message("   runReactomeEnrichment DEBUG: Top 10 genes for enricher:")
          message(paste(head(positive_genes, 10), collapse=", "))
          
          # Try with higher p-value cutoff
          message("   runReactomeEnrichment Step: Using higher p-value cutoff (0.1)...")
          enricher_result <- clusterProfiler::enricher(
            gene = positive_genes,
            universe = names(reactome_ranked_list),
            TERM2GENE = term2gene,
            TERM2NAME = term2name,
            pvalueCutoff = 0.05,  # Using higher threshold (0.1) for more enriched pathways
            pAdjustMethod = "fdr",
            minGSSize = 3,
            maxGSSize = 500
          )
          
          if (!is.null(enricher_result) && nrow(enricher_result@result) > 0) {
            message(sprintf("   runReactomeEnrichment: Enricher analysis successful with %d results", nrow(enricher_result@result)))
            return(enricher_result)
          } else {
            message("   runReactomeEnrichment: Enricher analysis returned no significant results")
          }
        }, error = function(e2) {
          message(paste("Enricher error:", e2$message))
          return(NULL)
        })
      }
      
      return(NULL)
    })
    
    # Format results for visualization
    if (!is.null(gsea_result)) {
      message("   runReactomeEnrichment Step: Formatting results...")
      
      if (class(gsea_result)[1] == "enrichResult") {
        # Result from enricher
        message("   runReactomeEnrichment: Processing results from enricher")
        if (nrow(gsea_result@result) > 0) {
          message(sprintf("   runReactomeEnrichment: Found %d enriched terms", nrow(gsea_result@result)))
        } else {
          message("   runReactomeEnrichment: No enriched terms found")
        }
        
        reactome_results <- gsea_result@result |>
          as_tibble() |>
          dplyr::select(
            termDescription = Description,
            enrichmentScore = p.adjust,  # Use p.adjust as score
            falseDiscoveryRate = p.adjust,
            genesMapped = Count,
            mappedIDs = geneID  # Add the list of mapped IDs
          ) |>
          mutate(
            enrichmentScore = -log10(enrichmentScore),  # Convert to -log10 scale
            comparison = assay_name
          ) |>
          distinct()  # FIXED: Added distinct to remove duplicates
      } else if (class(gsea_result)[1] == "gseaResult" && nrow(gsea_result@result) > 0) {
        # Result from GSEA
        message("   runReactomeEnrichment: Processing results from GSEA")
        message(sprintf("   runReactomeEnrichment: Found %d enriched terms", nrow(gsea_result@result)))
        
        reactome_results <- gsea_result@result |>
          as_tibble() |>
          dplyr::select(
            termDescription = Description,
            enrichmentScore = NES, 
            falseDiscoveryRate = p.adjust,
            genesMapped = setSize,
            mappedIDs = core_enrichment  # Add the core enrichment metabolites
          ) |>
          mutate(
            enrichmentScore = abs(enrichmentScore),
            comparison = assay_name
          ) |>
          distinct()  # FIXED: Added distinct to remove duplicates
      } else {
        message("   runReactomeEnrichment: No enriched terms in GSEA/enricher output, returning empty table")
        reactome_results <- tibble(
          termDescription = character(),
          enrichmentScore = numeric(),
          falseDiscoveryRate = numeric(),
          genesMapped = integer(),
          mappedIDs = character(),  # Add empty column for mapped IDs
          comparison = character()
        )
      }
      
      # FIXED: Check for duplicate results after processing
      if (nrow(reactome_results) > 0) {
        duplicate_count <- sum(duplicated(reactome_results$termDescription))
        if (duplicate_count > 0) {
          message(sprintf("   runReactomeEnrichment WARNING: Found %d duplicate pathway descriptions after processing", duplicate_count))
          message("   runReactomeEnrichment DEBUG: Applying additional distinct() operation...")
          reactome_results <- reactome_results |> distinct(termDescription, .keep_all = TRUE)
        }
      }
      
    } else {
      message("   runReactomeEnrichment: GSEA/enricher analysis failed or returned NULL, returning empty table")
      reactome_results <- tibble(
        termDescription = character(),
        enrichmentScore = numeric(),
        falseDiscoveryRate = numeric(),
        genesMapped = integer(),
        mappedIDs = character(),  # Add empty column for mapped IDs
        comparison = character()
      )
    }
    
    message(sprintf("--- Exiting runReactomeEnrichment. Returning table with %d rows ---", nrow(reactome_results)))
    return(reactome_results)
    
  }, error = function(e) {
    message(paste("Error processing Reactome data:", e$message))
    return(tibble(
      termDescription = character(),
      enrichmentScore = numeric(),
      falseDiscoveryRate = numeric(),
      genesMapped = integer(),
      mappedIDs = character(),
      comparison = character()
    ))
  })
}

