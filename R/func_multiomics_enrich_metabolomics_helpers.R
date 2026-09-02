# runMetabolomicsEnrichmentAnalysis
# ----------------------------------------------------------------------------
#' Run Metabolomics Enrichment Analysis for MOFA Factors
#'
#' @param weights Data frame containing MOFA weights
#' @param metabolomics_obj Metabolomics S4 object containing ChEBI IDs
#' @param mapping_table The metabolite ID mapping table from AnnotationHub
#' @param project_dirs Project directories structure
#' @param omic_type Omics type for directory access
#' @param experiment_label Experiment label for directory access
#' @param assay_name Name of the assay ("metabolome_lc" or "metabolome_gc")
#' @param kegg_species_code KEGG species code (e.g., "kpn" for Klebsiella pneumoniae)
#' @param reactome_organism Organism name to use for Reactome filtering (e.g., "Homo sapiens")
#' @return A data frame with enrichment results formatted for visualization
#' @export
runMetabolomicsEnrichmentAnalysis <- function(weights, 
                                             metabolomics_obj,
                                             mapping_table,
                                             project_dirs,
                                             omic_type,
                                             experiment_label,
                                             assay_name,
                                             kegg_species_code = "kpn",
                                             reactome_organism = NULL) {
  
  message(sprintf("--- Entering runMetabolomicsEnrichmentAnalysis ---"))
  message(sprintf("   runMetabolomicsEnrichmentAnalysis Args: assay_name = %s", assay_name))
  message(sprintf("   runMetabolomicsEnrichmentAnalysis Args: kegg_species_code = %s", kegg_species_code))
  message(sprintf("   runMetabolomicsEnrichmentAnalysis Args: reactome_organism = %s", as.character(reactome_organism)))
  
  # Use getProjectPaths helper with automatic fallback
  data_dir <- getProjectPaths(omic_type, experiment_label)$data_dir
  results_dir <- getProjectPaths(omic_type, experiment_label)$integration_enrichment_plots_dir
  
  # Create directory structure if it doesn't exist
  dir.create(file.path(results_dir, assay_name), recursive = TRUE, showWarnings = FALSE)
  
  # 1. Extract metabolite weights from MOFA for the specific assay
  message("   runMetabolomicsEnrichmentAnalysis Step: Extracting metabolite weights for assay...")
  metabolite_weights <- weights |>
    dplyr::filter(view == assay_name & factor == "Factor1") |>
    mutate(feature = str_replace_all(feature, paste0("_", assay_name), ""))
  
  message(sprintf("   runMetabolomicsEnrichmentAnalysis: Found %d metabolite weights", nrow(metabolite_weights)))
  message("   runMetabolomicsEnrichmentAnalysis: Sample weights:")
  print(head(metabolite_weights, 5))
  
  # 2. Get ChEBI IDs from metabolomics object
  message("   runMetabolomicsEnrichmentAnalysis Step: Getting ChEBI IDs from metabolomics object...")
  assay_index <- if(assay_name == "metabolome_lc") 1 else 2
  
  # Extract ChEBI IDs and corresponding metabolite names
  chebi_mapping <- metabolomics_obj@metabolite_data[[assay_index]] |>
    dplyr::filter(!stringr::str_detect(database_identifier, "^ITSD") & 
                 !stringr::str_detect(metabolite_identification, "^ITSD")) |>
    dplyr::select(metabolite, database_identifier) |>
    # Extract just the ChEBI ID number from the database_identifier column
    dplyr::mutate(chebi_id = stringr::str_extract(database_identifier, "CHEBI:\\d+"))
  
  message(sprintf("   runMetabolomicsEnrichmentAnalysis: Extracted %d ChEBI IDs", nrow(chebi_mapping)))
  message("   runMetabolomicsEnrichmentAnalysis: Sample ChEBI mappings:")
  print(head(chebi_mapping, 5))
  
  # Check missing ChEBI IDs
  missing_chebi_count <- sum(is.na(chebi_mapping$chebi_id))
  if (missing_chebi_count > 0) {
    message(sprintf("   runMetabolomicsEnrichmentAnalysis WARNING: %d metabolites missing ChEBI IDs", missing_chebi_count))
    message("Sample entries with missing ChEBI IDs:")
    print(head(chebi_mapping[is.na(chebi_mapping$chebi_id),], 5))
  }
  
  # 3. Join metabolite weights with ChEBI IDs
  message("   runMetabolomicsEnrichmentAnalysis Step: Joining weights with ChEBI IDs...")
  metabolite_weights_with_ids <- metabolite_weights |>
    left_join(chebi_mapping, by = c("feature" = "metabolite"))
  
  message(sprintf("   runMetabolomicsEnrichmentAnalysis: Joined table has %d rows", nrow(metabolite_weights_with_ids)))
  message("   runMetabolomicsEnrichmentAnalysis: Sample joined data:")
  print(head(metabolite_weights_with_ids, 5))
  
  # 4. Filter out entries without ChEBI IDs
  message("   runMetabolomicsEnrichmentAnalysis Step: Filtering out entries without ChEBI IDs...")
  metabolite_weights_with_ids <- metabolite_weights_with_ids |>
    dplyr::filter(!is.na(chebi_id))
  
  message(sprintf("   runMetabolomicsEnrichmentAnalysis: After filtering, %d entries remain", nrow(metabolite_weights_with_ids)))
  
  # Create ranked gene list for GSEA
  message("   runMetabolomicsEnrichmentAnalysis Step: Creating ranked list for GSEA...")
  ranked_list <- metabolite_weights_with_ids |>
    dplyr::arrange(desc(value)) |>
    dplyr::pull(value, name = chebi_id)
  
  message(sprintf("   runMetabolomicsEnrichmentAnalysis: Created ranked list with %d entries", length(ranked_list)))
  message("   runMetabolomicsEnrichmentAnalysis: Top metabolites in ranked list:")
  print(head(ranked_list, 5))
  
  # Run KEGG pathway enrichment analysis
  message("   runMetabolomicsEnrichmentAnalysis Step: Running KEGG pathway enrichment...")
  kegg_results <- runKeggEnrichment(
    ranked_list = ranked_list, 
    mapping_table = mapping_table,
    project_dirs = project_dirs,
    omic_type = omic_type,
    experiment_label = experiment_label,
    assay_name = assay_name,
    kegg_species_code = kegg_species_code
  )
  
  message(sprintf("   runMetabolomicsEnrichmentAnalysis: KEGG analysis returned %d results", nrow(kegg_results)))
  
  # Run Reactome pathway enrichment analysis  
  message("   runMetabolomicsEnrichmentAnalysis Step: Running Reactome pathway enrichment...")
  reactome_results <- runReactomeEnrichment(
    ranked_list = ranked_list, 
    mapping_table = mapping_table,
    project_dirs = project_dirs,
    omic_type = omic_type,
    experiment_label = experiment_label,
    assay_name = assay_name,
    reactome_organism = reactome_organism
  )
  
  message(sprintf("   runMetabolomicsEnrichmentAnalysis: Reactome analysis returned %d results", nrow(reactome_results)))
  
  # Combine results and format for visualization
  message("   runMetabolomicsEnrichmentAnalysis Step: Combining KEGG and Reactome results...")
  combined_results <- bind_rows(
    kegg_results |> mutate(category = "KEGG"),
    reactome_results |> mutate(category = "Reactome")
  ) |> 
  distinct()  # Ensure no duplicates after combining
  
  message(sprintf("   runMetabolomicsEnrichmentAnalysis: Combined %d results total", nrow(combined_results)))
  
  # Save the combined results
  message(sprintf("   runMetabolomicsEnrichmentAnalysis Step: Saving results to %s...", 
                 file.path(results_dir, paste0(assay_name, "_enrichment_results.tab"))))
  vroom::vroom_write(combined_results, 
                    file.path(results_dir, paste0(assay_name, "_enrichment_results.tab")))
  
  message("--- Exiting runMetabolomicsEnrichmentAnalysis ---")
  return(combined_results)
}

multiomicsFirstLookup <- function(keys, values) {
  keep <- !is.na(keys) & nzchar(keys) & !duplicated(keys)
  entries <- as.list(as.character(values[keep]))
  names(entries) <- keys[keep]
  list2env(entries, envir = new.env(hash = TRUE, parent = emptyenv()))
}

multiomicsLookupValue <- function(lookup, key) {
  get0(
    key,
    envir = lookup,
    mode = "character",
    inherits = FALSE,
    ifnotfound = NA_character_
  )
}

multiomicsPathwayNameLookups <- function(all_names_mapping, kegg_to_chebi) {
  list(
    kegg_to_chebi = multiomicsFirstLookup(
      kegg_to_chebi$kegg_id,
      kegg_to_chebi$chebi_id
    ),
    assay_chebi = multiomicsFirstLookup(
      paste(all_names_mapping$assay, all_names_mapping$chebi_id, sep = "\r"),
      all_names_mapping$metabolite
    ),
    chebi = multiomicsFirstLookup(
      all_names_mapping$chebi_id,
      all_names_mapping$metabolite
    )
  )
}

multiomicsMapPathwayIds <- function(id_string, assay_type, lookups) {
  if (is.null(id_string) || !length(id_string) ||
      isTRUE(is.na(id_string[[1L]])) || identical(id_string[[1L]], "")) {
    return(NA_character_)
  }
  ids <- unlist(strsplit(as.character(id_string), split = ",|/"))
  names_vec <- vapply(ids, \(id) {
    if (is.na(id) || !nzchar(id)) return(NA_character_)
    chebi_id <- if (startsWith(id, "cpd:")) {
      multiomicsLookupValue(lookups$kegg_to_chebi, id)
    } else if (startsWith(id, "CHEBI:")) {
      id
    } else {
      return(id)
    }
    if (is.na(chebi_id)) return(id)
    assay_name <- multiomicsLookupValue(
      lookups$assay_chebi,
      paste(assay_type, chebi_id, sep = "\r")
    )
    if (!is.na(assay_name)) return(assay_name)
    fallback <- multiomicsLookupValue(lookups$chebi, chebi_id)
    if (is.na(fallback)) id else fallback
  }, character(1))
  paste(names_vec[!is.na(names_vec)], collapse = ", ")
}


# ----------------------------------------------------------------------------
# runMetabolomicsPathwayEnrichment
# ----------------------------------------------------------------------------
#' Run Metabolomics Pathway Enrichment for Both Assays
#'
#' @param weights Data frame containing MOFA weights
#' @param metabolomics_obj Metabolomics S4 object containing ChEBI IDs
#' @param mapping_table The metabolite ID mapping table from AnnotationHub
#' @param project_dirs Project directories structure
#' @param omic_type Omics type for directory access
#' @param experiment_label Experiment label for directory access
#' @param kegg_species_code KEGG species code (e.g., "kpn" for Klebsiella pneumoniae)
#' @param reactome_organism Organism name to use for Reactome filtering (optional)
#' @return A combined data frame with enrichment results for both assays
#' @export
runMetabolomicsPathwayEnrichment <- function(weights, 
                                            metabolomics_obj, 
                                            mapping_table, 
                                            project_dirs,
                                            omic_type,
                                            experiment_label,
                                            kegg_species_code = "kpn",
                                            reactome_organism = NULL) {
  
  message(sprintf("--- Entering runMetabolomicsPathwayEnrichment ---"))
  message(sprintf("   runMetabolomicsPathwayEnrichment Args: kegg_species_code = %s", kegg_species_code))
  message(sprintf("   runMetabolomicsPathwayEnrichment Args: reactome_organism = %s", as.character(reactome_organism)))
  
  # Use getProjectPaths helper with automatic fallback
  results_dir <- getProjectPaths(omic_type, experiment_label)$integration_enrichment_plots_dir
  
  # Inspect the structure of input objects
  message("   runMetabolomicsPathwayEnrichment Step: Inspecting metabolomics object structure...")
  message(sprintf("   runMetabolomicsPathwayEnrichment: Metabolomics object has %d assays", length(metabolomics_obj@metabolite_data)))
  message(sprintf("   runMetabolomicsPathwayEnrichment: Assay 1 (LC) has %d metabolites", nrow(metabolomics_obj@metabolite_data[[1]])))
  message(sprintf("   runMetabolomicsPathwayEnrichment: Assay 2 (GC) has %d metabolites", nrow(metabolomics_obj@metabolite_data[[2]])))
  
  message("   runMetabolomicsPathwayEnrichment Step: Inspecting weights data frame...")
  message(sprintf("   runMetabolomicsPathwayEnrichment: Weights data frame has %d entries", nrow(weights)))
  metabolome_lc_weights <- sum(weights$view == "metabolome_lc")
  metabolome_gc_weights <- sum(weights$view == "metabolome_gc")
  message(sprintf("   runMetabolomicsPathwayEnrichment: Found %d LC-MS weights and %d GC-MS weights", 
                 metabolome_lc_weights, metabolome_gc_weights))
  
  message("   runMetabolomicsPathwayEnrichment Step: Inspecting mapping table...")
  message(sprintf("   runMetabolomicsPathwayEnrichment: Mapping table has %d entries", nrow(mapping_table)))
  kegg_mappings <- sum(!is.na(mapping_table$KEGG))
  chebi_mappings <- sum(!is.na(mapping_table$ChEBI))
  message(sprintf("   runMetabolomicsPathwayEnrichment: Found %d KEGG IDs and %d ChEBI IDs in mapping table", 
                 kegg_mappings, chebi_mappings))
  
  # Run enrichment for LC-MS metabolomics
  message("   runMetabolomicsPathwayEnrichment Step: Running LC-MS metabolomics enrichment...")
  lc_results <- runMetabolomicsEnrichmentAnalysis(
    weights = weights,
    metabolomics_obj = metabolomics_obj,
    mapping_table = mapping_table,
    project_dirs = project_dirs,
    omic_type = omic_type,
    experiment_label = experiment_label,
    assay_name = "metabolome_lc",
    kegg_species_code = kegg_species_code,
    reactome_organism = reactome_organism
  )
  
  message(sprintf("   runMetabolomicsPathwayEnrichment: LC-MS analysis returned %d results", nrow(lc_results)))
  
  # Run enrichment for GC-MS metabolomics
  message("   runMetabolomicsPathwayEnrichment Step: Running GC-MS metabolomics enrichment...")
  gc_results <- runMetabolomicsEnrichmentAnalysis(
    weights = weights,
    metabolomics_obj = metabolomics_obj,
    mapping_table = mapping_table,
    project_dirs = project_dirs,
    omic_type = omic_type,
    experiment_label = experiment_label,
    assay_name = "metabolome_gc",
    kegg_species_code = kegg_species_code,
    reactome_organism = reactome_organism
  )
  
  message(sprintf("   runMetabolomicsPathwayEnrichment: GC-MS analysis returned %d results", nrow(gc_results)))
  
  # Combine results
  message("   runMetabolomicsPathwayEnrichment Step: Combining LC-MS and GC-MS results...")
  combined_results <- bind_rows(
    lc_results |> mutate(assay = "LC-MS"),
    gc_results |> mutate(assay = "GC-MS")
  ) |>
  mutate(comparison = "Metabolome") |> # Set overall comparison label
  distinct()  # Ensure no duplicates
  
  message(sprintf("   runMetabolomicsPathwayEnrichment: Combined results has %d entries", nrow(combined_results)))
  
  # Add original metabolite names
  message("   runMetabolomicsPathwayEnrichment Step: Adding original metabolite names...")
  
  # Create a mapping from KEGG IDs to metabolite names
  message("   runMetabolomicsPathwayEnrichment Step: Creating ID to name mapping tables...")
  
  # Extract all metabolite names and IDs from metabolomics object
  lc_names_mapping <- metabolomics_obj@metabolite_data[[1]] |>
    dplyr::select(metabolite, database_identifier) |>
    dplyr::filter(!is.na(database_identifier)) |>
    dplyr::mutate(chebi_id = stringr::str_extract(database_identifier, "CHEBI:\\d+"))
    
  gc_names_mapping <- metabolomics_obj@metabolite_data[[2]] |>
    dplyr::select(metabolite, database_identifier) |>
    dplyr::filter(!is.na(database_identifier)) |>
    dplyr::mutate(chebi_id = stringr::str_extract(database_identifier, "CHEBI:\\d+"))
    
  # Combine mappings from both assays
  all_names_mapping <- bind_rows(
    lc_names_mapping |> mutate(assay = "LC-MS"),
    gc_names_mapping |> mutate(assay = "GC-MS")
  )
  
  message(sprintf("   runMetabolomicsPathwayEnrichment: Created mapping table with %d entries", nrow(all_names_mapping)))
  
  # Also get KEGG to ChEBI mapping for translating KEGG IDs
  kegg_to_chebi <- mapping_table |>
    dplyr::filter(!is.na(KEGG) & !is.na(ChEBI)) |>
    dplyr::select(KEGG, ChEBI) |>
    dplyr::mutate(
      kegg_id = paste0("cpd:", KEGG),
      chebi_id = paste0("CHEBI:", ChEBI)
    )
  
  # Apply the mapping function to each row
  message("   runMetabolomicsPathwayEnrichment Step: Mapping IDs to names...")

  if (nrow(combined_results) == 0) {
    message("   runMetabolomicsPathwayEnrichment: No results to map. Adding empty mappedNames column.")
    combined_results$mappedNames <- character(0)
  } else {
    if (!"mappedIDs" %in% colnames(combined_results)) {
      message("   runMetabolomicsPathwayEnrichment WARNING: mappedIDs column missing! Adding empty mappedNames.")
      combined_results$mappedNames <- rep(NA_character_, nrow(combined_results))
    } else {
      lookups <- multiomicsPathwayNameLookups(
        all_names_mapping,
        kegg_to_chebi
      )
      combined_results$mappedNames <- purrr::map2_chr(
        combined_results$mappedIDs,
        combined_results$assay,
        \(id_string, assay_type) multiomicsMapPathwayIds(
          id_string,
          assay_type,
          lookups
        )
      )
      combined_results <- tibble::as_tibble(combined_results)
    }
  }
  
  # Save combined results
  message(sprintf("   runMetabolomicsPathwayEnrichment Step: Saving combined results to %s...", 
                 file.path(results_dir, "combined_metabolomics_enrichment_results.tab")))
  vroom::vroom_write(combined_results, 
                    file.path(results_dir, "combined_metabolomics_enrichment_results.tab"))
  
  message("--- Exiting runMetabolomicsPathwayEnrichment ---")
  return(combined_results)
}


# ----------------------------------------------------------------------------
