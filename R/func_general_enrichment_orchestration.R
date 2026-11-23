# stripEnrichmentIdSuffix
# ----------------------------------------------------------------------------
#' An enrichment id without its colon suffix (KIT-119)
#'
#' The one rule `processEnrichments` applies to the ids it sends to g:Profiler, the query and the custom background
#' alike: everything from the first colon is dropped (`P12345:iso1` becomes `P12345`). An id without a colon is
#' unchanged.
#' @param ids Character vector of ids.
#' @return `ids` without their colon suffixes.
#' @noRd
stripEnrichmentIdSuffix <- function(ids) {
  stringr::str_remove(ids, ":.*")
}

#' The proteins of each `geneID` cell in alphabetical order (KIT-173, VAL-61)
#'
#' enrichit's `ora_cpp` writes a cell in the iteration order of a hash set filled from `Rcpp::unique()`, whose order
#' follows the strings' memory addresses, so the same analysis listed the same proteins in another order in another R
#' process. Sorting fixes the order; counts, sizes, ratios, p-values and the sets are unchanged.
#' @param result An `enrichResult`; anything else (NULL, a stand-in without `geneID`) is returned as it is.
#' @return `result` with each `geneID` cell's proteins sorted (byte order, `method = "radix"`).
#' @noRd
sortEnrichmentGeneIds <- function(result) {
  if (!isS4(result) || !"result" %in% methods::slotNames(result) || !is.data.frame(result@result) ||
      !nrow(result@result) || !"geneID" %in% names(result@result)) return(result)
  result@result$geneID <- vapply(strsplit(as.character(result@result$geneID), "/", fixed = TRUE),
    function(proteins) paste(sort(proteins, method = "radix"), collapse = "/"), character(1))
  result
}

# processEnrichments
# ----------------------------------------------------------------------------
#' Process Enrichments
#'
#' @param da_results S4 object containing differential expression results
#' @param taxon_id NCBI taxonomy ID for the organism
#' @param up_cutoff Log2 fold change cutoff for up-regulated proteins (default: 0)
#' @param down_cutoff Log2 fold change cutoff for down-regulated proteins (default: 0)
#' @param q_cutoff FDR q-value threshold for enrichment significance (default: 0.05)
#' @param pathway_dir Directory for saving pathway results
#' @param go_annotations UniProt GO annotations (required for unsupported organisms)
#' @param exclude_iea Whether to exclude IEA (Inferred Electronic Annotation) terms (default: FALSE). It applies to
#'   the 13 organisms served by g:Profiler only: the UniProt annotations used for any other organism carry no evidence
#'   codes, so for those `TRUE` gives a warning (class `multischolar_enrichment_iea_not_excludable`) and electronic
#'   annotations are included.
#' @param protein_id_column Name of the protein ID column (default: "Protein.IDs")
#' @param contrast_names Vector of contrast names for output labeling
#' @param correction_method Method for FDR correction (default: "gSCS")
#' @param go_propagation For organisms outside the 13 served by g:Profiler, how the Gene Ontology sets built from
#'   `go_annotations` take parent terms (KIT-173): `"is_a_part_of_regulates"` (the default) adds to each term the
#'   proteins annotated to its more specific terms over the relations is a, part of, regulates, positively regulates
#'   and negatively regulates, as g:Profiler and GO.db's ancestor tables do (a biological-process term then also holds
#'   proteins annotated as regulating that process); `"is_a_part_of"` follows is a and part of only, the relations the
#'   GO Consortium calls safe for grouping annotations; `"none"` uses the direct annotations only (the sets before
#'   KIT-173). On a propagating value the three aspect roots are not tested, so proteins annotated only to them leave
#'   the background, and the GO release is that of the installed GO.db. Not used on the g:Profiler route.
#'
#' @return S4 EnrichmentResults object containing enrichment data and summaries. Its plot slots are left empty: the
#'   plots are written to `pathway_dir` and rebuilt from the enrichment data when shown.
#'
#' @export
processEnrichments <- function(da_results,
                               taxon_id,
                               up_cutoff = 0,
                               down_cutoff = 0,
                               q_cutoff = 0.05,
                               pathway_dir,
                               go_annotations = NULL,
                               exclude_iea = FALSE,
                               protein_id_column = "Protein.IDs",
                               contrast_names = NULL,
                               correction_method = "gSCS",
                               go_propagation = goPropagationShipped) {

  message("--- RUNNING processEnrichments VERSION [Timestamp: ", Sys.time(), "] ---")

  # Validate exclude_iea parameter
  if (is.null(exclude_iea)) {
    stop("exclude_iea must be explicitly set to TRUE or FALSE")
  }

  if (!is.logical(exclude_iea)) {
    stop("exclude_iea must be a logical value (TRUE or FALSE)")
  }

  # Common model organisms lookup
  supported_organisms <- tibble::tribble(
    ~taxid,     ~id,            ~name,
    "9606",     "hsapiens",     "Homo sapiens",
    "10090",    "mmusculus",    "Mus musculus",
    "10116",    "rnorvegicus",  "Rattus norvegicus",
    "7227",     "dmelanogaster", "Drosophila melanogaster",
    "6239",     "celegans",     "Caenorhabditis elegans",
    "4932",     "scerevisiae",  "Saccharomyces cerevisiae",
    "3702",     "athaliana",    "Arabidopsis thaliana",
    "7955",     "drerio",       "Danio rerio",
    "9031",     "ggallus",      "Gallus gallus",
    "9823",     "sscrofa",      "Sus scrofa",
    "9913",     "btaurus",      "Bos taurus",
    "9544",     "mmulatta",     "Macaca mulatta",
    "9598",     "ptroglodytes", "Pan troglodytes"
  )

  is_supported <- as.character(taxon_id) %in% supported_organisms$taxid

  if(is_supported) {
    message(sprintf("Taxon ID %s found in supported organisms. Proceeding with gprofiler2 analysis...", taxon_id))

    # Convert taxon_id to species
    species <- supported_organisms |>
      dplyr::filter(.data$taxid == as.character(taxon_id)) |>
      dplyr::pull(.data$id)

    enrichment_results <- createEnrichmentResults(da_results@contrasts)

    # Process each contrast
    results <- da_results@da_data |>
      purrr::map(function(de_data) {
        tryCatch({
        if(is.null(de_data)) {
          warning("No DE data found for contrast")
          return(NULL)
        }

        # Split data into up/down regulated
        message(sprintf("Total proteins before filtering: %d", nrow(de_data)))
        
        # [OK] DIAGNOSTIC: Check protein_id_column and data structure
        
        if (!protein_id_column %in% names(de_data)) {
          stop(sprintf("ERROR: Column '%s' not found in DE data. Available columns: %s", 
                      protein_id_column, paste(names(de_data), collapse = ", ")))
        }
        
        # Check if the column has valid data
        protein_col_data <- de_data[[protein_id_column]]

         
         na_count_before <- sum(is.na(de_data[[protein_id_column]]))
         
         passes_q_filter <- de_data |> dplyr::filter(.data$fdr_qvalue < q_cutoff)
         if (nrow(passes_q_filter) > 0) {
           na_in_sig <- sum(is.na(passes_q_filter[[protein_id_column]]))
         }

         subset_sig <- de_data |>
            dplyr::mutate(
              !!rlang::sym(protein_id_column) := stripEnrichmentIdSuffix(.data[[protein_id_column]])
            ) |>
           dplyr::filter(.data$fdr_qvalue < q_cutoff)
         
         na_count_after_q <- sum(is.na(subset_sig[[protein_id_column]]))

         if(nrow(subset_sig) > 0) {
            # Also check for the literal "Protein.Ids" to be sure
            if(any(subset_sig[[protein_id_column]] == "Protein.Ids", na.rm = TRUE)) {
              message("      WARNING: Literal 'Protein.Ids' found in the protein ID column after modification!")
            }
         }

        message(sprintf("Proteins passing q-value cutoff (%.3f): %d", q_cutoff, nrow(subset_sig)))
        
        # No longer need subset_for_enrichment, as subset_sig is already filtered
        message(sprintf("Proteins available for enrichment: %d", nrow(subset_sig)))
        
        up_matrix <- subset_sig |>
          dplyr::filter(.data$log2FC > up_cutoff)
        message(sprintf("Up-regulated proteins (log2FC > %g): %d", up_cutoff, nrow(up_matrix)))
        
        if (nrow(up_matrix) > 0) {
          na_in_up <- sum(is.na(up_matrix[[protein_id_column]]))
        }
        
        down_matrix <- subset_sig |>
          dplyr::filter(.data$log2FC < -down_cutoff)
        message(sprintf("Down-regulated proteins (log2FC < -%g): %d", down_cutoff, nrow(down_matrix)))
        
        if (nrow(down_matrix) > 0) {
          na_in_down <- sum(is.na(down_matrix[[protein_id_column]]))
        }
        
        # Get background IDs from the full de_data for this contrast. KIT-119: stripped as the query ids are
        # (they were not, so a query id with a colon suffix was missing from its own background).
        custom_bg <- stripEnrichmentIdSuffix(de_data[[protein_id_column]]) |>
          unique()

        message(sprintf("Using %d unique proteins as background for enrichment analysis", length(custom_bg)))
        
        # Process up and down regulated genes
        list(
          up = tryCatch({
            if(nrow(up_matrix) > 0) {
              # [OK] DEBUG 66: Log call to perform_enrichment for up-regulated
              message(sprintf("      About to analyze %d up-regulated proteins", nrow(up_matrix)))
              
              protein_col <- protein_id_column
              message(sprintf("      Using protein_id_column: %s", protein_col))
              message(sprintf("      Using species: %s", species))
              message(sprintf("      Using threshold: %s", q_cutoff))
              message(sprintf("      Using exclude_iea: %s", exclude_iea))
              
              up_result <- perform_enrichment(
                data_subset = up_matrix,
                species = species,
                threshold = q_cutoff,
                sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"),
                domain_scope = "custom",
                custom_bg = custom_bg,
                exclude_iea = exclude_iea,
                protein_id_column = protein_col,
                correction_method = correction_method
              )
              
              up_result
            } else {
              NULL
            }
          }, error = function(e) {
            warning(sprintf("Error processing up-regulated genes: %s", e$message))
            message("Debug info for up-regulated genes:")
            message("Number of proteins: ", nrow(up_matrix))
            NULL
          }),

          down = tryCatch({
            if(nrow(down_matrix) > 0) {
              # [OK] DEBUG 66: Log call to perform_enrichment for down-regulated
              message(sprintf("      About to analyze %d down-regulated proteins", nrow(down_matrix)))
              
              protein_col <- protein_id_column
              message(sprintf("      Using protein_id_column: %s", protein_col))
              message(sprintf("      Using species: %s", species))
              message(sprintf("      Using threshold: %s", q_cutoff))
              message(sprintf("      Using exclude_iea: %s", exclude_iea))
              
              down_result <- perform_enrichment(
                data_subset = down_matrix,
                species = species,
                threshold = q_cutoff,
                sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"),
                domain_scope = "custom",
                custom_bg = custom_bg,
                exclude_iea = exclude_iea,
                protein_id_column = protein_col,
                correction_method = correction_method
              )
              
              down_result
            } else {
              NULL
            }
          }, error = function(e) {
            warning(sprintf("Error processing down-regulated genes: %s", e$message))
            message("Debug info for down-regulated genes:")
            message("Number of proteins: ", nrow(down_matrix))
            NULL
          })
        )
        }, error = function(e) {
          message(sprintf("*** ERROR in contrast processing: %s", e$message))
          message(sprintf("*** ERROR details: %s", class(e)))
          message(sprintf("*** Call stack: %s", paste(deparse(sys.calls()), collapse = " -> ")))
          return(NULL)
        })
      })

    # Store enrichment results
    enrichment_results@enrichment_data <- results

    # Generate both static and interactive plots for the saved files. KIT-119: they are not stored; the
    # `enrichment_plots` and `enrichment_plotly` slots stay empty and the tab rebuilds a plot from
    # `enrichment_data` with the same builder (`buildEnrichmentPlot`).
    plot_results <- purrr::map(names(results), function(contrast) {
      list(
        up = if(!is.null(results[[contrast]]$up)) {
          generate_enrichment_plots(results[[contrast]]$up, contrast, "up", pathway_dir)
        } else NULL,

        down = if(!is.null(results[[contrast]]$down)) {
          generate_enrichment_plots(results[[contrast]]$down, contrast, "down", pathway_dir)
        } else NULL
      )
    }) |>
      purrr::set_names(names(results))

    # Generate and store summaries
    enrichment_results@enrichment_summaries <- purrr::map(names(results), function(contrast) {
      list(
        up = if(!is.null(results[[contrast]]$up)) summarize_enrichment(results[[contrast]]$up) else NULL,
        down = if(!is.null(results[[contrast]]$down)) summarize_enrichment(results[[contrast]]$down) else NULL
      )
    }) |>
      purrr::set_names(names(results))

    # Explicitly ensure the names of plot_results match contrast_names_to_use
    # This guards against issues if names were lost or mismatched during creation
    if (!identical(names(plot_results), contrast_names)) {
       # Check if the number of elements still matches before trying to assign names
       if(length(plot_results) == length(contrast_names)) {
         names(plot_results) <- contrast_names
       } else {
         stop("Critical error: Number of plot results does not match the number of contrast names.")
       }
    }

    # Save plots using the desired contrast names for files
    purrr::iwalk(contrast_names, function(contrast, i) {
      plots <- plot_results[[contrast]]

      # Double-check if plots object is NULL, which might happen if naming failed
      if(is.null(plots)) {
          warning(paste("Could not find plot data for contrast:", contrast, "- Skipping save for this contrast."))
          return(NULL) # Skip to the next iteration
      }

      # Simple sanitization (should be redundant now)
      sanitized_contrast <- gsub("[^A-Za-z0-9_.-]", "_", contrast)

      message(paste("Loop", i, "- Saving plots for contrast:", contrast, " (Sanitized:", sanitized_contrast, ")")) # Debug

      purrr::walk(c("up", "down"), function(direction) {
        # Check if the plot component exists and is not NULL
        if (!is.null(plots[[direction]]) && !is.null(plots[[direction]]$interactive)) {

          file_name <- paste0(sanitized_contrast, "_", direction, "_enrichment_plot.html")
          file_path <- file.path(pathway_dir, file_name)
          
          # Add PNG file path
          png_file_name <- paste0(sanitized_contrast, "_", direction, "_enrichment_plot.png")
          png_file_path <- file.path(pathway_dir, png_file_name)

          message(paste("  Attempting to save:", file_path)) # Debug

          # Ensure the directory exists before saving
          target_dir_for_file <- dirname(file_path)
          if (!dir.exists(target_dir_for_file)) {
             message(paste("  Creating directory:", target_dir_for_file)) # Debug
             dir.create(target_dir_for_file, recursive = TRUE)
          }

          tryCatch({
            # Save the interactive Plotly HTML file
            saveSelfContainedWidget(plots[[direction]]$interactive, file_path) # KIT-036: no pandoc pass
            message(paste("  Successfully saved:", file_path)) # Debug success
            
            # Save the static ggplot as PNG
            ggplot2::ggsave(
              filename = png_file_path,
              plot = plots[[direction]]$static,
              width = 10, 
              height = 8,
              dpi = 300,
              bg = "white"
            )
            message(paste("  Successfully saved:", png_file_path)) # Debug success
          }, error = function(e) {
            # Print more detailed error context
            warning(paste("  ERROR saving plots for contrast:", contrast,
                          "direction:", direction,
                          "path:", file_path,
                          "- Error message:", e$message))
          })
        } else {
           message(paste("  Skipping save for direction:", direction, "- Plot component is NULL or not interactive."))
        }
      })
    })

    return(enrichment_results)

  } else {
    if(is.null(go_annotations)) {
      stop("For unsupported organisms, GO annotations must be provided")
    }

    message(sprintf("Using custom GO annotations for taxon ID %s", taxon_id))
    checkGoPropagation(go_propagation) # KIT-173
    if (isTRUE(exclude_iea)) {
      cli::cli_warn(paste("Electronic annotations cannot be excluded for this organism: the UniProt annotations used",
        "here carry no evidence codes. They are included."), class = "multischolar_enrichment_iea_not_excludable")
    }

    # Prepare GO term mappings (KIT-173: with parent terms unless `go_propagation` is "none")
    gene_sets <- buildUniprotGoGeneSets(go_annotations, go_propagation)
    all_terms <- gene_sets$all_terms
    term2gene <- gene_sets$term2gene
    term2name <- gene_sets$term2name

    # Get the internal long names (only used initially if short names aren't provided or needed for mapping)
    internal_contrast_names <- names(da_results@da_data)

    # Determine which names to use (Prefer explicitly passed short names)
    if (is.null(contrast_names)) {
      warning("Explicit contrast_names not provided, using internal names from da_results@da_data which might be long or contain invalid characters.")
      contrast_names_to_use <- internal_contrast_names
    } else {
      if(length(contrast_names) != length(internal_contrast_names)) {
        stop("Length of provided 'contrast_names' does not match the number of contrasts in 'da_results@da_data'.")
      }
      contrast_names_to_use <- contrast_names # Use the short names
    }
    print(contrast_names_to_use)

    # Initialize lists
    results_list_long_names <- list() # Temporary list to hold results with original names
    enrichment_results <- createEnrichmentResults(da_results@contrasts)

    # --- Step 1: Process enrichment using internal loop mapped to short names ---
    message("--- Starting Enrichment Processing Loop ---")
    results_list_long_names <- purrr::map2(contrast_names_to_use, internal_contrast_names, function(short_name, long_name) {
        message(paste("Processing contrast:", short_name, "(original:", long_name, ")"))

        de_data <- da_results@da_data[[long_name]] # Access input data using long name

        if(is.null(de_data)) {
          warning(paste("No DE data found for internal contrast:", long_name))
          # Return NULL for both up and down to keep list structure aligned
          return(list(up = NULL, down = NULL))
        }

        # Split data into up/down regulated
        message(sprintf("Total proteins before filtering: %d", nrow(de_data)))

        subset_sig <- de_data |>
          dplyr::filter(.data$fdr_qvalue < q_cutoff)

        message(sprintf("Proteins passing FDR cutoff (%g): %d", q_cutoff, nrow(subset_sig)))

        up_genes <- subset_sig |>
          dplyr::filter(.data$log2FC > up_cutoff) |>
          dplyr::pull({{protein_id_column}})

        down_genes <- subset_sig |>
          dplyr::filter(.data$log2FC < -down_cutoff) |>
          dplyr::pull({{protein_id_column}})

        # Background genes
        background_IDs <- unique(de_data |> dplyr::pull({{protein_id_column}}))

        # Perform enrichment for up-regulated genes
        up_enrich <- tryCatch({
          if(length(up_genes) > 0) {
            sortEnrichmentGeneIds(clusterProfiler::enricher(
              gene = up_genes,
              universe = background_IDs,
              TERM2GENE = term2gene,
              TERM2NAME = term2name,
              pvalueCutoff = q_cutoff,
              pAdjustMethod = "BH"
            )) # KIT-173 (VAL-61): geneID order fixed
          } else NULL
        }, error = function(e) {
          warning(sprintf("Error processing up-regulated genes: %s", e$message))
          NULL
        })

        # Perform enrichment for down-regulated genes
        down_enrich <- tryCatch({
          if(length(down_genes) > 0) {
            sortEnrichmentGeneIds(clusterProfiler::enricher(
              gene = down_genes,
              universe = background_IDs,
              TERM2GENE = term2gene,
              TERM2NAME = term2name,
              pvalueCutoff = q_cutoff,
              pAdjustMethod = "BH"
            )) # KIT-173 (VAL-61): geneID order fixed
          } else NULL
        }, error = function(e) {
          warning(sprintf("Error processing down-regulated genes: %s", e$message))
          NULL
        })

        # Return results for this contrast
        list(up = up_enrich, down = down_enrich)
    }) |> purrr::set_names(internal_contrast_names)
    message("--- Finished Enrichment Processing Loop ---")

    # --- Step 2: Assign results to the S4 object with SHORT names ---
    if(length(results_list_long_names) == length(contrast_names_to_use)) {
        names(results_list_long_names) <- contrast_names_to_use # Rename the temporary list
        enrichment_results@enrichment_data <- results_list_long_names # Assign renamed list
    } else {
        stop("Mismatch between number of processed results and expected contrast names.")
    }
    # Now enrichment_results@enrichment_data uses SHORT names

    # Create GO term mappings once (moved outside the plotting function). KIT-173: with parent terms, the names come
    # from the gene sets' GO.db lookup, which covers every tested term.
    go_term_map <- if (!identical(go_propagation, "none")) {
      data.frame(ID = term2name$TERM, term = term2name$NAME, stringsAsFactors = FALSE)
    } else dplyr::bind_rows(
      go_annotations |>
        tidyr::separate_rows(.data$go_id_go_biological_process, .data$go_term_go_biological_process, sep = "; ") |>
        dplyr::select(.data$go_id_go_biological_process, .data$go_term_go_biological_process) |>
        dplyr::rename(ID = .data$go_id_go_biological_process, term = .data$go_term_go_biological_process),

      go_annotations |>
        tidyr::separate_rows(.data$go_id_go_molecular_function, .data$go_term_go_molecular_function, sep = "; ") |>
        dplyr::select(.data$go_id_go_molecular_function, .data$go_term_go_molecular_function) |>
        dplyr::rename(ID = .data$go_id_go_molecular_function, term = .data$go_term_go_molecular_function),

      go_annotations |>
        tidyr::separate_rows(.data$go_id_go_cellular_compartment, .data$go_term_go_cellular_compartment, sep = "; ") |>
        dplyr::select(.data$go_id_go_cellular_compartment, .data$go_term_go_cellular_compartment) |>
        dplyr::rename(ID = .data$go_id_go_cellular_compartment, term = .data$go_term_go_cellular_compartment)
    ) |>
      dplyr::distinct()

    # Create category mapping once
    go_category_map <- all_terms |>
      dplyr::distinct(.data$TERM, .data$ONTOLOGY) |>
      dplyr::mutate(
        source = dplyr::case_when(
          .data$ONTOLOGY == "BP" ~ "GO:BP",
          .data$ONTOLOGY == "CC" ~ "GO:CC",
          .data$ONTOLOGY == "MF" ~ "GO:MF"
        )
      )

    # --- Step 3: Generate Plots using SHORT names ---
    message("--- Starting Plot Generation Loop ---")
    plot_results_list <- purrr::map(contrast_names_to_use, function(contrast) {
        message(paste("Generating plots for contrast:", contrast))
        # Access enrichment data using SHORT name
        current_enrich_data <- enrichment_results@enrichment_data[[contrast]]

        # Process up and down directions
        direction_results <- purrr::map(c("up", "down"), function(direction) {
            result_data <- current_enrich_data[[direction]] # Access up/down list

            if(!is.null(result_data) && nrow(result_data@result) > 0) {
                message(sprintf("Processing %s-regulated genes for contrast %s", direction, contrast))

                # Prepare data for plotting and tables
                plot_data <- result_data@result |>
                  dplyr::left_join(go_category_map, by = c("ID" = "TERM")) |>
                  dplyr::left_join(go_term_map, by = "ID") |>
                  dplyr::mutate(
                    source = dplyr::coalesce(.data$source, "Other"),
                    source = factor(.data$source, levels = c("GO:BP", "GO:CC", "GO:MF", "Other")),
                    neg_log10_q = -log10(.data$qvalue),  # Using qvalue directly from clusterProfiler output
                    gene_count = .data$Count,
                    significant = .data$qvalue < q_cutoff  # Add significance flag based on q_cutoff
                  ) |>
                  dplyr::mutate(
                    term = dplyr::coalesce(.data$term, .data$Description)
                  )

                # Save results table
                readr::write_tsv(
                  plot_data,
                  file.path(pathway_dir,
                            paste0(contrast, "_", direction, "_enrichment_results.tsv"))
                )

                # KIT-119: the stored result carries its terms' GO category and name and the cutoff, and the plot
                # is built from that stored result by the builder the tab uses at render time.
                result_data <- setEnrichmentPlotInputs(result_data, plot_data, q_cutoff)
                c(buildEnrichmentPlot(result_data, contrast, direction), list(result = result_data))
            } else {
                NULL # Return NULL if no results
            }
        }) |> purrr::set_names(c("up", "down"))
        
        # Return plot components for this contrast
        direction_results
    }) |> purrr::set_names(contrast_names_to_use)
    message("--- Finished Plot Generation Loop ---")

    # --- Step 4: Store what the plots are rebuilt from (KIT-119) ---
    # No plot is stored: the `enrichment_plots` and `enrichment_plotly` slots stay empty. Each stored result
    # carries its plot inputs, and the tab rebuilds the plot from it (`buildEnrichmentPlot`).
    for (contrast in contrast_names_to_use) {
      for (direction in c("up", "down")) {
        stored_result <- plot_results_list[[contrast]][[direction]]$result
        if (!is.null(stored_result)) {
          enrichment_results@enrichment_data[[contrast]][[direction]] <- stored_result
        }
      }
    }

    # --- Step 5: Save interactive plots using SHORT names ---
    message("--- Starting Final Save Loop ---")
    purrr::walk(contrast_names_to_use, function(contrast) {
        # Access plot data using SHORT name
        plots <- plot_results_list[[contrast]]

        if(is.null(plots)) {
            # This warning should NOT trigger if steps above worked
            warning(paste("Could not find plot data for contrast:", contrast, "- Skipping save. (This indicates an earlier issue)"))
            return(NULL)
        }

        # Sanitize (optional if short names are clean)
        sanitized_contrast <- gsub("[^A-Za-z0-9_.-]", "_", contrast)
        message(paste("Saving plots for contrast:", contrast))

        # Process both directions
        purrr::walk(c("up", "down"), function(direction) {
            # Check if interactive plot exists
            if (!is.null(plots[[direction]]$interactive)) { # Access up/down list
                interactive_plot <- plots[[direction]]$interactive # The actual plotly object
                static_plot <- plots[[direction]]$static # The static ggplot

                file_name <- paste0(sanitized_contrast, "_", direction, "_enrichment_plot.html")
                file_path <- file.path(pathway_dir, file_name)
                
                # Add PNG file path
                png_file_name <- paste0(sanitized_contrast, "_", direction, "_enrichment_plot.png")
                png_file_path <- file.path(pathway_dir, png_file_name)
                
                message(paste("  Attempting to save:", file_path))

                # Ensure directory exists
                target_dir_for_file <- dirname(file_path)
                if (!dir.exists(target_dir_for_file)) {
                   dir.create(target_dir_for_file, recursive = TRUE)
                }

                # Save the widget
                tryCatch({
                  # Save HTML interactive plot
                  saveSelfContainedWidget(interactive_plot, file_path) # KIT-036: no pandoc pass
                  message(paste("  Successfully saved:", file_path))
                  
                  # Save static ggplot as PNG
                  if (!is.null(static_plot)) {
                    ggplot2::ggsave(
                      filename = png_file_path,
                      plot = static_plot,
                      width = 10, 
                      height = 8,
                      dpi = 300,
                      bg = "white"
                    )
                    message(paste("  Successfully saved:", png_file_path))
                  }
                }, error = function(e) {
                  warning(paste("  ERROR saving plots:", file_path, "-", e$message))
                })
            } else {
                message(paste("  Skipping save for direction:", direction, "- Plot component is NULL."))
            }
        })
    })
    message("--- Finished Final Save Loop ---")

    return(enrichment_results)
  }
}


# ----------------------------------------------------------------------------
# getEnrichmentResult
# ----------------------------------------------------------------------------
# Helper function to access results
#' @export
getEnrichmentResult <- function(enrichment_results, contrast, direction) {
  if(!contrast %in% names(enrichment_results@enrichment_data)) {
    stop("Contrast not found")
  }
  if(!direction %in% c("up", "down")) {
    stop("Direction must be 'up' or 'down'")
  }
  enrichment_results@enrichment_data[[contrast]][[direction]]
}


# ----------------------------------------------------------------------------
# getEnrichmentPlotly
# ----------------------------------------------------------------------------
# Helper function to access plotly objects
#' @export
getEnrichmentPlotly <- function(enrichment_results, contrast, direction) {
  # KIT-119: a run no longer stores its plots; the plot is rebuilt from the stored result unless one is stored.
  if(!contrast %in% c(names(enrichment_results@enrichment_plotly), names(enrichment_results@enrichment_data))) {
    stop("Contrast not found")
  }
  if(!direction %in% c("up", "down")) {
    stop("Direction must be 'up' or 'down'")
  }
  resolveEnrichmentPlotly(enrichment_results, contrast, direction)
}


# ----------------------------------------------------------------------------
# getEnrichmentSummary
# ----------------------------------------------------------------------------
# Helper function to get summary
#' @export

getEnrichmentSummary <- function(enrichment_results) {
  summaries <- purrr::map_df(names(enrichment_results@enrichment_summaries), function(contrast) {
    summary <- enrichment_results@enrichment_summaries[[contrast]]

    data.frame(
      contrast = contrast,
      up_total = if(!is.null(summary$up)) summary$up$total else 0,
      down_total = if(!is.null(summary$down)) summary$down$total else 0,
      up_GO_BP = if(!is.null(summary$up)) summary$up$GO_BP else 0,
      down_GO_BP = if(!is.null(summary$down)) summary$down$GO_BP else 0,
      up_GO_CC = if(!is.null(summary$up)) summary$up$GO_CC else 0,
      down_GO_CC = if(!is.null(summary$down)) summary$down$GO_CC else 0,
      up_GO_MF = if(!is.null(summary$up)) summary$up$GO_MF else 0,
      down_GO_MF = if(!is.null(summary$down)) summary$down$GO_MF else 0,
      up_KEGG = if(!is.null(summary$up)) summary$up$KEGG else 0,
      down_KEGG = if(!is.null(summary$down)) summary$down$KEGG else 0,
      up_REAC = if(!is.null(summary$up)) summary$up$REAC else 0,
      down_REAC = if(!is.null(summary$down)) summary$down$REAC else 0
    )
  })

  return(summaries)
}
