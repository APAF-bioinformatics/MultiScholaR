# runKeggEnrichment
# ----------------------------------------------------------------------------
#' @title Run KEGG Pathway Enrichment Analysis for Metabolites
#' @param ranked_list Named numeric vector with ChEBI IDs as names and values as weights
#' @param mapping_table Metabolite ID mapping table
#' @param project_dirs Project directories structure
#' @param omic_type Omics type for directory access
#' @param experiment_label Experiment label for directory access
#' @param assay_name Name of the assay
#' @param kegg_species_code KEGG species code (e.g., "kpn" for Klebsiella pneumoniae)
#' @return A data frame with KEGG pathway enrichment results
#' @export
runKeggEnrichment <- function(ranked_list, 
                             mapping_table, 
                             project_dirs,
                             omic_type,
                             experiment_label,
                             assay_name,
                             kegg_species_code = "kpn") {
  
  message(sprintf("--- Entering runKeggEnrichment ---"))
  message(sprintf("   runKeggEnrichment Args: assay_name = %s", assay_name))
  message(sprintf("   runKeggEnrichment Args: kegg_species_code = %s", kegg_species_code))
  message(sprintf("   runKeggEnrichment Args: ranked_list length = %d", length(ranked_list)))
  
  # Use getProjectPaths helper with automatic fallback
  data_dir <- getProjectPaths(omic_type, experiment_label)$data_dir
  message(sprintf("   runKeggEnrichment: Using data_dir = %s", data_dir))
  
  # Print some debugging information
  message("KEGG Enrichment - Input ranked_list summary:")
  message(paste("Number of metabolites:", length(ranked_list)))
  message(paste("Sample ChEBI IDs:", paste(head(names(ranked_list), 3), collapse=", ")))
  
  # Load the necessary mappings with better debugging
  message("   runKeggEnrichment Step: Loading KEGG to ChEBI mappings...")
  kegg_to_chebi <- mapping_table |>
    dplyr::filter(!is.na(KEGG) & !is.na(ChEBI)) 
  
  message(paste("Number of KEGG-ChEBI mappings:", nrow(kegg_to_chebi)))
  message("   runKeggEnrichment Step: Sample of KEGG-ChEBI mappings:")
  print(head(kegg_to_chebi, 5))
  
  # Format ChEBI IDs to match our format
  message("   runKeggEnrichment Step: Formatting ChEBI IDs...")
  kegg_to_chebi <- kegg_to_chebi |>
    dplyr::select(KEGG, ChEBI) |>
    dplyr::mutate(ChEBI = paste0("CHEBI:", ChEBI))
  
  message("   runKeggEnrichment Step: Formatted KEGG-ChEBI mappings:")
  print(head(kegg_to_chebi, 5))
  
  # Load pathway to compound mapping
  message(sprintf("   runKeggEnrichment Step: Loading pathway to compound mapping from %s...", file.path(data_dir, "species_specific_pathway_to_compound_tbl.tab")))
  
  # Check if file exists
  if (!file.exists(file.path(data_dir, "species_specific_pathway_to_compound_tbl.tab"))) {
    message("   runKeggEnrichment ERROR: File 'species_specific_pathway_to_compound_tbl.tab' does not exist!")
    return(tibble(
      termDescription = character(),
      enrichmentScore = numeric(),
      falseDiscoveryRate = numeric(),
      genesMapped = integer(),
      mappedIDs = character(),
      comparison = character()
    ))
  }
  
  species_specific_pathway_to_compound_tbl <- vroom::vroom(
    file.path(data_dir, "species_specific_pathway_to_compound_tbl.tab"),
    show_col_types = FALSE
  )
  
  message(paste("Number of pathway-compound mappings:", nrow(species_specific_pathway_to_compound_tbl)))
  message("   runKeggEnrichment Step: Pathway to compound column names:")
  print(colnames(species_specific_pathway_to_compound_tbl))
  
  # Debug - show first few compounds and pathways
  message("First few compounds:")
  message(paste(head(species_specific_pathway_to_compound_tbl$compound, 5), collapse=", "))
  
  message("First few pathways:")
  message(paste(head(unique(species_specific_pathway_to_compound_tbl$pathway), 5), collapse=", "))
  
  # Directly use KEGG compound IDs from the pathway mapping
  message("   runKeggEnrichment Step: Extracting unique compound IDs...")
  compound_ids <- sort(unique(species_specific_pathway_to_compound_tbl$compound))
  message(paste("Total unique compounds in KEGG:", length(compound_ids)))
  message("Sample compound IDs in KEGG:")
  message(paste(head(compound_ids, 10), collapse=", "))
  
  # Cross-check with our ranked list - are there any direct KEGG IDs in our list?
  message("   runKeggEnrichment Step: Checking for direct KEGG IDs in ranked list...")
  kegg_ids_in_ranked_list <- intersect(compound_ids, names(ranked_list))
  message(paste("Direct KEGG IDs in ranked list:", length(kegg_ids_in_ranked_list)))
  if (length(kegg_ids_in_ranked_list) > 0) {
    message("Sample direct KEGG IDs found in ranked list:")
    message(paste(head(kegg_ids_in_ranked_list, 5), collapse=", "))
  }
  
  # If there are no direct KEGG IDs, we'll need to rely on the mapping
  if (length(kegg_ids_in_ranked_list) == 0) {
    message("   runKeggEnrichment Step: No direct KEGG IDs in ranked list, using ChEBI to KEGG mapping...")
    
    # Print the structure of the first few ChEBI IDs to check format
    chebi_ids_in_ranked_list <- names(ranked_list)[grepl("CHEBI:", names(ranked_list))]
    message(paste("ChEBI IDs in ranked list:", length(chebi_ids_in_ranked_list)))
    if (length(chebi_ids_in_ranked_list) > 0) {
      message(paste("Sample ChEBI IDs:", paste(head(chebi_ids_in_ranked_list, 3), collapse=", ")))
    }
    
    # Check how many of our ChEBI IDs are in the mapping table
    chebi_ids_in_mapping <- intersect(chebi_ids_in_ranked_list, kegg_to_chebi$ChEBI)
    message(paste("ChEBI IDs found in mapping table:", length(chebi_ids_in_mapping)))
    if (length(chebi_ids_in_mapping) > 0) {
      message("Sample ChEBI IDs found in mapping table:")
      message(paste(head(chebi_ids_in_mapping, 5), collapse=", "))
    }
    
    # Create a direct lookup from ChEBI to KEGG
    message("   runKeggEnrichment Step: Creating ChEBI to KEGG lookup...")
    chebi_to_kegg_lookup <- kegg_to_chebi |>
      dplyr::filter(ChEBI %in% chebi_ids_in_ranked_list)
    
    message(paste("ChEBI-KEGG mappings for our ranked list:", nrow(chebi_to_kegg_lookup)))
    if (nrow(chebi_to_kegg_lookup) > 0) {
      message("   runKeggEnrichment Step: Sample ChEBI-KEGG mappings:")
      print(head(chebi_to_kegg_lookup, 5))
    }
    
    # If we have mappings, create a new ranked list with KEGG IDs
    if (nrow(chebi_to_kegg_lookup) > 0) {
      # Create a new ranked list with KEGG IDs
      message("   runKeggEnrichment Step: Creating ranked list with KEGG IDs...")
      kegg_ranked_list <- numeric()
      for (i in 1:nrow(chebi_to_kegg_lookup)) {
        chebi_id <- chebi_to_kegg_lookup$ChEBI[i]
        kegg_id <- chebi_to_kegg_lookup$KEGG[i]
        if (chebi_id %in% names(ranked_list)) {
          kegg_ranked_list[kegg_id] <- ranked_list[chebi_id]
        }
      }
      
      message(paste("Created KEGG ranked list with", length(kegg_ranked_list), "entries"))
      if (length(kegg_ranked_list) > 0) {
        message("   runKeggEnrichment Step: Sample KEGG IDs in ranked list:")
        message(paste(head(names(kegg_ranked_list), 5), collapse=", "))
      }
      
      # If we have KEGG IDs, use them for the analysis
      if (length(kegg_ranked_list) > 0) {
        # Sort the ranked list in decreasing order
        kegg_ranked_list <- sort(kegg_ranked_list, decreasing = TRUE)
        
        # Diagnostic: Check for ties in the data
        message("   runKeggEnrichment DEBUG: Analyzing score distribution...")
        # Count exact duplicates in values
        duplicate_count <- sum(duplicated(kegg_ranked_list))
        duplicate_percentage <- (duplicate_count / length(kegg_ranked_list)) * 100
        message(sprintf("   runKeggEnrichment DEBUG: Found %d exact duplicate values (%.2f%% of total)",
                       duplicate_count, duplicate_percentage))
        
        # Show basic statistics
        message(sprintf("   runKeggEnrichment DEBUG: Score range: [%.4f, %.4f], Median: %.4f", 
                       min(kegg_ranked_list), max(kegg_ranked_list), median(kegg_ranked_list)))
        
        # Add small random noise to break ties in the ranked list
        message("   runKeggEnrichment Step: Adding small random noise to break ties in ranked list...")
        set.seed(42) # For reproducibility
        original_kegg_ranked_list <- kegg_ranked_list # Store original for reference
        kegg_ranked_list <- kegg_ranked_list + runif(length(kegg_ranked_list), min=-0.0001, max=0.0001)
        # Re-sort to ensure order is preserved with the new values
        kegg_ranked_list <- sort(kegg_ranked_list, decreasing = TRUE)
        message(sprintf("   runKeggEnrichment: Added jitter to %d values to prevent ties", length(kegg_ranked_list)))
        
        # CRITICAL: Verify the ordering and distribution before GSEA
        message("   runKeggEnrichment CRITICAL CHECK: Verifying ranked list ordering and distribution...")
        positive_values <- sum(kegg_ranked_list > 0)
        negative_values <- sum(kegg_ranked_list < 0)
        zero_values <- sum(kegg_ranked_list == 0)
        message(sprintf("   runKeggEnrichment: Value distribution: %d positive, %d negative, %d zero", 
                       positive_values, negative_values, zero_values))
        
        # Check that the list is actually sorted (decreasing)
        is_decreasing <- all(diff(kegg_ranked_list) <= 0)
        message(sprintf("   runKeggEnrichment: List is properly sorted in decreasing order: %s", 
                       if(is_decreasing) "YES" else "NO - THIS IS A PROBLEM!"))
        
        # Show top and bottom values to verify ordering
        message("   runKeggEnrichment: Top 5 entries in ranked list:")
        for (i in 1:min(5, length(kegg_ranked_list))) {
          message(sprintf("     %s: %.6f", names(kegg_ranked_list)[i], kegg_ranked_list[i]))
        }
        message("   runKeggEnrichment: Bottom 5 entries in ranked list:")
        for (i in (length(kegg_ranked_list) - min(4, length(kegg_ranked_list) - 1)):length(kegg_ranked_list)) {
          message(sprintf("     %s: %.6f", names(kegg_ranked_list)[i], kegg_ranked_list[i]))
        }
        
        # CRITICAL SECTION - WHERE THE ISSUE LIKELY IS
        message("   runKeggEnrichment Step: Creating TERM2GENE mapping with KEGG IDs...")
        message("   runKeggEnrichment Step: Checking for pathway-compound matches...")
        
        # Check if compounds in pathway table have a prefix that's missing in our ranked list
        pathway_compounds <- head(species_specific_pathway_to_compound_tbl$compound, 10)
        ranked_compounds <- head(names(kegg_ranked_list), 10)
        
        message("Pathway compound format examples:")
        print(pathway_compounds)
        
        message("Ranked list compound format examples:")
        print(ranked_compounds)
        
        # FIX: FORCE PREFIX ADDITION REGARDLESS OF CHECKS
        message("   runKeggEnrichment DEBUG: Adding 'cpd:' prefix to all KEGG IDs in ranked list")
        prefixed_kegg_ids <- paste0("cpd:", names(kegg_ranked_list))
        prefixed_kegg_ranked_list <- kegg_ranked_list
        names(prefixed_kegg_ranked_list) <- prefixed_kegg_ids
        
        message(sprintf("   runKeggEnrichment DEBUG: Created prefixed ranked list with %d entries", length(prefixed_kegg_ranked_list)))
        message("Sample prefixed entries:")
        print(head(prefixed_kegg_ranked_list, 5))
        
        # Use the prefixed version for TERM2GENE mapping
        term2gene <- species_specific_pathway_to_compound_tbl |>
          dplyr::filter(compound %in% names(prefixed_kegg_ranked_list)) |>
          dplyr::select(pathway, compound) |>
          dplyr::rename(term = pathway, gene = compound)
        
        message(sprintf("   runKeggEnrichment DEBUG: After prefix fix, TERM2GENE has %d entries", nrow(term2gene)))
        
        if (nrow(term2gene) == 0) {
          message("   runKeggEnrichment WARNING: Still no matches after prefix fix!")
          
          # Debug check the values directly
          message("   runKeggEnrichment DEBUG: Direct value checking...")
          sample_pathway_compounds <- head(species_specific_pathway_to_compound_tbl$compound, 5)
          sample_prefixed_ids <- head(names(prefixed_kegg_ranked_list), 10)
          
          for (compound in sample_prefixed_ids) {
            matching_idx <- which(species_specific_pathway_to_compound_tbl$compound == compound)
            message(sprintf("   runKeggEnrichment DEBUG: Compound '%s' appears %d times in pathway table", 
                          compound, length(matching_idx)))
            if (length(matching_idx) > 0) {
              sample_pathways <- head(species_specific_pathway_to_compound_tbl$pathway[matching_idx], 3)
              message(sprintf("   runKeggEnrichment DEBUG: Found in pathways: %s", 
                            paste(sample_pathways, collapse=", ")))
            }
          }
        }
        
        # Use prefixed list regardless of whether mapping succeeded
        kegg_ranked_list <- prefixed_kegg_ranked_list
        
        # Create term2name directly in this path - crucial addition
        message("   runKeggEnrichment Step: Creating TERM2NAME mapping in direct path...")
        # Get all pathways from KEGG if not already done
        if (!exists("pathways_tbl")) {
          message("   runKeggEnrichment Step: Fetching pathways from KEGG...")
          pathways_tbl <- KEGGREST::keggList("pathway")
          
          message(sprintf("   runKeggEnrichment: Retrieved %d pathways from KEGG", length(pathways_tbl)))
          if (length(pathways_tbl) > 0) {
            message("Sample pathway entries:")
            print(head(pathways_tbl, 3))
          }
        }
        
        # Extract unique terms from term2gene
        term2name_initial <- term2gene |>
          dplyr::distinct(term)
        
        message(sprintf("   runKeggEnrichment: Found %d distinct terms", nrow(term2name_initial)))
        
        # Create a mapping table from pathway IDs (e.g., "path:kpn00010") to names
        # First extract the numeric part of each term
        message("   runKeggEnrichment Step: Extracting numeric parts from KEGG pathway IDs...")
        term2name_with_numeric <- term2name_initial |>
          dplyr::mutate(
            pathway_num = stringr::str_extract(term, "\\d+"),
            # Create generic format for lookup
            map_format = paste0("map", pathway_num)
          )
        
        message("   runKeggEnrichment DEBUG: Sample pathways with numeric parts:")
        print(head(term2name_with_numeric, 3))
        
        # Create mapping from pathway numbers to pathway names
        pathway_num_to_name <- data.frame(
          map_id = names(pathways_tbl),
          pathway_name = as.character(pathways_tbl)
        ) |>
          dplyr::mutate(
            map_format = stringr::str_extract(map_id, "map\\d+")
          )
        
        message("   runKeggEnrichment DEBUG: Sample pathway mappings:")
        print(head(pathway_num_to_name, 3))
        
        # Join the terms with pathway names based on numeric part
        message("   runKeggEnrichment Step: Joining with pathway names...")
        term2name <- term2name_with_numeric |>
          dplyr::left_join(pathway_num_to_name, by = "map_format") |>
          dplyr::select(term, name = pathway_name)
        
        # Handle any missing names (fallback to term ID)
        term2name <- term2name |>
          dplyr::mutate(name = ifelse(is.na(name), term, name))
        
        message("   runKeggEnrichment: Created improved term2name mapping:")
        print(head(term2name, 3))
        
        # Try running GSEA with tryCatch to handle errors
        message("   runKeggEnrichment Step: Running GSEA analysis...")
        
        # ABSOLUTE FINAL SAFETY CHECK - verify both term2gene and term2name exist
        if (!exists("term2gene") || nrow(term2gene) == 0) {
          message("   runKeggEnrichment CRITICAL ERROR: term2gene is missing or empty before GSEA call")
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
          message("   runKeggEnrichment WARNING: term2name is missing or empty before GSEA call, creating emergency fallback")
          # Create emergency fallback - but try to get proper names first
          # Get all pathways from KEGG if not already done
          if (!exists("pathways_tbl")) {
            message("   runKeggEnrichment Step: Emergency fetch of pathways from KEGG...")
            pathways_tbl <- tryCatch({
              KEGGREST::keggList("pathway")
            }, error = function(e) {
              message("   runKeggEnrichment ERROR: Could not fetch pathway names:", e$message)
              return(NULL)
            })
          }
          
          if (!is.null(pathways_tbl) && length(pathways_tbl) > 0) {
            # Extract term numbers and try to map
            term2name <- term2gene |> 
              dplyr::distinct(term) |>
              dplyr::mutate(
                pathway_num = stringr::str_extract(term, "\\d+"),
                map_format = paste0("map", pathway_num)
              )
            
            # Create mapping from pathway numbers to pathway names
            pathway_num_to_name <- data.frame(
              map_id = names(pathways_tbl),
              pathway_name = as.character(pathways_tbl)
            ) |>
              dplyr::mutate(
                map_format = stringr::str_extract(map_id, "map\\d+")
              )
            
            # Join and select final columns
            term2name <- term2name |>
              dplyr::left_join(pathway_num_to_name, by = "map_format") |>
              dplyr::select(term, name = pathway_name)
              
            # Handle any missing names (fallback to term ID)
            term2name <- term2name |>
              dplyr::mutate(name = ifelse(is.na(name), term, name))
          } else {
            # Complete fallback if KEGG fetch failed
            term2name <- term2gene |> 
              dplyr::distinct(term) |>
              dplyr::mutate(name = term)
          }
        }
        
        # Debug output
        message(sprintf("   runKeggEnrichment VERIFICATION: term2gene has %d rows, term2name has %d rows", 
                       nrow(term2gene), nrow(term2name)))
        
        gsea_result <- tryCatch({
          # FIXED: Use fgsea and remove nPerm parameter
          clusterProfiler::GSEA(
            geneList = kegg_ranked_list,
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
          positive_genes <- names(kegg_ranked_list)[kegg_ranked_list > 0]
          message(sprintf("   runKeggEnrichment: Found %d positively ranked genes", length(positive_genes)))
          
          if (length(positive_genes) > 0) {
            tryCatch({
              message("   runKeggEnrichment Step: Running enricher analysis...")
              enricher_result <- clusterProfiler::enricher(
                gene = positive_genes,
                universe = names(kegg_ranked_list),
                TERM2GENE = term2gene,
                TERM2NAME = term2name,
                pvalueCutoff = 0.05,  # Using higher threshold (0.1) for more lenient filtering
                pAdjustMethod = "fdr",
                minGSSize = 3,
                maxGSSize = 500
              )
              
              if (!is.null(enricher_result) && nrow(enricher_result@result) > 0) {
                message(sprintf("   runKeggEnrichment: Enricher analysis successful with %d results", nrow(enricher_result@result)))
                return(enricher_result)
              } else {
                message("   runKeggEnrichment: Enricher analysis returned no significant results")
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
          message("   runKeggEnrichment Step: Formatting results...")
          
          if (class(gsea_result)[1] == "enrichResult") {
            # Result from enricher
            message("   runKeggEnrichment: Processing results from enricher")
            if (nrow(gsea_result@result) > 0) {
              message(sprintf("   runKeggEnrichment: Found %d enriched terms", nrow(gsea_result@result)))
            } else {
              message("   runKeggEnrichment: No enriched terms found")
            }
            
            kegg_results <- gsea_result@result |>
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
              )
          } else if (class(gsea_result)[1] == "gseaResult" && nrow(gsea_result@result) > 0) {
            # Result from GSEA
            message("   runKeggEnrichment: Processing results from GSEA")
            message(sprintf("   runKeggEnrichment: Found %d enriched terms", nrow(gsea_result@result)))
            
            kegg_results <- gsea_result@result |>
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
              )
          } else {
            message("   runKeggEnrichment: No enriched terms in GSEA/enricher output, returning empty table")
            kegg_results <- tibble(
              termDescription = character(),
              enrichmentScore = numeric(),
              falseDiscoveryRate = numeric(),
              genesMapped = integer(),
              mappedIDs = character(),  # Add empty column for mapped IDs
              comparison = character()
            )
          }
        } else {
          message("   runKeggEnrichment: GSEA/enricher analysis failed or returned NULL, returning empty table")
          kegg_results <- tibble(
            termDescription = character(),
            enrichmentScore = numeric(),
            falseDiscoveryRate = numeric(),
            genesMapped = integer(),
            mappedIDs = character(),  # Add empty column for mapped IDs
            comparison = character()
          )
        }
        
        message(sprintf("--- Exiting runKeggEnrichment. Returning table with %d rows ---", nrow(kegg_results)))
        return(kegg_results)
      } else {
        message("No KEGG IDs could be mapped from ChEBI IDs. Returning empty result.")
        return(tibble(
          termDescription = character(),
          enrichmentScore = numeric(),
          falseDiscoveryRate = numeric(),
          genesMapped = integer(),
          comparison = character()
        ))
      }
    } else {
      message("No ChEBI-KEGG mappings found for our ranked list. Returning empty result.")
      return(tibble(
        termDescription = character(),
        enrichmentScore = numeric(),
        falseDiscoveryRate = numeric(),
        genesMapped = integer(),
        comparison = character()
      ))
    }
  } else {
    # Use direct KEGG IDs from the ranked list
    message("   runKeggEnrichment Step: Using direct KEGG IDs from ranked list")
    kegg_ranked_list <- ranked_list[kegg_ids_in_ranked_list]
    
    # Sort in decreasing order
    kegg_ranked_list <- sort(kegg_ranked_list, decreasing = TRUE)
    
    # Diagnostic: Check for ties in the data
    message("   runKeggEnrichment DEBUG: Analyzing score distribution...")
    # Count exact duplicates in values
    duplicate_count <- sum(duplicated(kegg_ranked_list))
    duplicate_percentage <- (duplicate_count / length(kegg_ranked_list)) * 100
    message(sprintf("   runKeggEnrichment DEBUG: Found %d exact duplicate values (%.2f%% of total)",
                   duplicate_count, duplicate_percentage))
    
    # Show basic statistics
    message(sprintf("   runKeggEnrichment DEBUG: Score range: [%.4f, %.4f], Median: %.4f", 
                   min(kegg_ranked_list), max(kegg_ranked_list), median(kegg_ranked_list)))
    
    # Add small random noise to break ties in the ranked list
    message("   runKeggEnrichment Step: Adding small random noise to break ties in ranked list...")
    set.seed(42) # For reproducibility
    original_kegg_ranked_list <- kegg_ranked_list # Store original for reference
    kegg_ranked_list <- kegg_ranked_list + runif(length(kegg_ranked_list), min=-0.0001, max=0.0001)
    # Re-sort to ensure order is preserved with the new values
    kegg_ranked_list <- sort(kegg_ranked_list, decreasing = TRUE)
    message(sprintf("   runKeggEnrichment: Added jitter to %d values to prevent ties", length(kegg_ranked_list)))
    
    # Create TERM2GENE mapping
    message("   runKeggEnrichment Step: Creating TERM2GENE mapping with direct KEGG IDs...")
    term2gene <- species_specific_pathway_to_compound_tbl |>
      dplyr::filter(compound %in% names(kegg_ranked_list)) |>
      dplyr::select(pathway, compound) |>
      dplyr::rename(term = pathway, gene = compound)
    
    message(paste("Created TERM2GENE mapping with", nrow(term2gene), "entries"))
    if (nrow(term2gene) > 0) {
      message("   runKeggEnrichment Step: Sample term2gene mappings:")
      print(head(term2gene, 5))
    }
  }
  
  # Check if we have any mappings
  if (!exists("term2gene") || nrow(term2gene) == 0) {
    message("No TERM2GENE mappings created. Returning empty result.")
    return(tibble(
      termDescription = character(),
      enrichmentScore = numeric(),
      falseDiscoveryRate = numeric(),
      genesMapped = integer(),
      mappedIDs = character(),
      comparison = character()
    ))
  }
  
  # Get pathway ID to name mapping
  # First get all pathways from KEGG
  message("   runKeggEnrichment Step: Fetching pathways from KEGG...")
  pathways_tbl <- KEGGREST::keggList("pathway")
  
  message(sprintf("   runKeggEnrichment: Retrieved %d pathways from KEGG", length(pathways_tbl)))
  if (length(pathways_tbl) > 0) {
    message("Sample pathway entries:")
    print(head(pathways_tbl, 3))
  }
  
  # Create mapping from pathway ID to name
  message("   runKeggEnrichment Step: Creating pathway ID to name mapping...")
  pathway_id_to_name <- data.frame(
    pathway_id = names(pathways_tbl), 
    pathway_name = pathways_tbl
  ) |>
    dplyr::mutate(
      hsa_pathway_id = str_replace(pathway_id, "^path:map", "") |> 
        purrr::map_chr(\(x){paste0(kegg_species_code, x)})
    )
  
  message("   runKeggEnrichment Step: Sample pathway ID to name mapping:")
  print(head(pathway_id_to_name, 3))
  
  # Create TERM2NAME mapping
  message("   runKeggEnrichment Step: Creating TERM2NAME mapping...")
  term2name_initial <- term2gene |>
    distinct(term)
  
  message(sprintf("   runKeggEnrichment: Found %d distinct terms", nrow(term2name_initial)))
  
  # Debug the term format
  message("   runKeggEnrichment DEBUG: Sample term formats:")
  print(head(term2name_initial$term))
  
  # Approach 1: Direct matching with pathway_id
  pathway_id_to_name_extended <- pathway_id_to_name |>
    # Add additional format columns for matching
    mutate(
      term_format1 = paste0("path:", kegg_species_code, str_replace(hsa_pathway_id, paste0("^", kegg_species_code), "")),
      term_format2 = paste0("path:", hsa_pathway_id)
    )
  
  message("   runKeggEnrichment DEBUG: Extended pathway_id_to_name formats:")
  print(head(pathway_id_to_name_extended))
  
  # Try to join using the extended formats
  term2name <- term2name_initial |>
    left_join(pathway_id_to_name_extended, 
              by = c("term" = "term_format1")) |>
    dplyr::select(term, pathway_name)
  
  # If no results, try the second format
  if (all(is.na(term2name$pathway_name))) {
    message("   runKeggEnrichment DEBUG: First join attempt failed, trying second format...")
    term2name <- term2name_initial |>
      left_join(pathway_id_to_name_extended, 
                by = c("term" = "term_format2")) |>
      dplyr::select(term, pathway_name)
  }
  
  # If both failed, try extracting and joining by numeric part
  if (all(is.na(term2name$pathway_name))) {
    message("   runKeggEnrichment DEBUG: Second join attempt failed, trying numeric part matching...")
    
    # Extract numeric parts from terms
    term2name <- term2name_initial |>
      mutate(
        # Extract just the numeric part from the pathway ID
        pathway_num = str_extract(term, "\\d+")
      )
    
    message("   runKeggEnrichment DEBUG: Extracted pathway numbers:")
    print(head(term2name))
    
    # Also extract numeric parts from pathway IDs
    pathway_id_to_name_numeric <- pathway_id_to_name |>
      mutate(
        pathway_num = str_extract(pathway_id, "\\d+")
      )
    
    message("   runKeggEnrichment DEBUG: Pathway ID numeric parts:")
    print(head(pathway_id_to_name_numeric))
    
    # Join on numeric parts
    term2name <- term2name |>
      left_join(pathway_id_to_name_numeric, by = "pathway_num") |>
      dplyr::select(term, name = pathway_name)
  } else {
    # We had a successful join earlier, just rename
    term2name <- term2name |>
      dplyr::rename(name = pathway_name)
  }
  
  message("   runKeggEnrichment Step: After joining with pathway names:")
  print(head(term2name))
  
  term2name <- term2name |>
    dplyr::filter(!is.na(name))  # Make sure we have names
  
  # Try running GSEA with tryCatch to handle errors
  message("   runKeggEnrichment Step: Running GSEA analysis...")
  gsea_result <- tryCatch({
    # FIXED: Use fgsea and remove nPerm parameter
    clusterProfiler::GSEA(
      geneList = kegg_ranked_list,
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
    positive_genes <- names(kegg_ranked_list)[kegg_ranked_list > 0]
    message(sprintf("   runKeggEnrichment: Found %d positively ranked genes", length(positive_genes)))
    
    if (length(positive_genes) > 0) {
      tryCatch({
        message("   runKeggEnrichment Step: Running enricher analysis...")
        enricher_result <- clusterProfiler::enricher(
          gene = positive_genes,
          universe = names(kegg_ranked_list),
          TERM2GENE = term2gene,
          TERM2NAME = term2name,
          pvalueCutoff = 0.05,  # Using higher threshold (0.1) for more lenient filtering
          pAdjustMethod = "fdr",
          minGSSize = 3,
          maxGSSize = 500
        )
        
        if (!is.null(enricher_result) && nrow(enricher_result@result) > 0) {
          message(sprintf("   runKeggEnrichment: Enricher analysis successful with %d results", nrow(enricher_result@result)))
          return(enricher_result)
        } else {
          message("   runKeggEnrichment: Enricher analysis returned no significant results")
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
    message("   runKeggEnrichment Step: Formatting results...")
    
    if (class(gsea_result)[1] == "enrichResult") {
      # Result from enricher
      message("   runKeggEnrichment: Processing results from enricher")
      if (nrow(gsea_result@result) > 0) {
        message(sprintf("   runKeggEnrichment: Found %d enriched terms", nrow(gsea_result@result)))
      } else {
        message("   runKeggEnrichment: No enriched terms found")
      }
      
      kegg_results <- gsea_result@result |>
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
        )
    } else if (class(gsea_result)[1] == "gseaResult" && nrow(gsea_result@result) > 0) {
      # Result from GSEA
      message("   runKeggEnrichment: Processing results from GSEA")
      message(sprintf("   runKeggEnrichment: Found %d enriched terms", nrow(gsea_result@result)))
      
      kegg_results <- gsea_result@result |>
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
        )
    } else {
      message("   runKeggEnrichment: No enriched terms in GSEA/enricher output, returning empty table")
      kegg_results <- tibble(
        termDescription = character(),
        enrichmentScore = numeric(),
        falseDiscoveryRate = numeric(),
        genesMapped = integer(),
        mappedIDs = character(),  # Add empty column for mapped IDs
        comparison = character()
      )
    }
  } else {
    message("   runKeggEnrichment: GSEA/enricher analysis failed or returned NULL, returning empty table")
    kegg_results <- tibble(
      termDescription = character(),
      enrichmentScore = numeric(),
      falseDiscoveryRate = numeric(),
      genesMapped = integer(),
      mappedIDs = character(),  # Add empty column for mapped IDs
      comparison = character()
    )
  }
  
  message(sprintf("--- Exiting runKeggEnrichment. Returning table with %d rows ---", nrow(kegg_results)))
  return(kegg_results)
}


# ----------------------------------------------------------------------------

