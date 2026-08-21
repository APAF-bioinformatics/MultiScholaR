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
#' @param exclude_iea Whether to exclude IEA (Inferred Electronic Annotation) terms (default: FALSE)
#' @param protein_id_column Name of the protein ID column (default: "Protein.IDs")
#' @param contrast_names Vector of contrast names for output labeling
#' @param correction_method Method for FDR correction (default: "gSCS")
#' @param execution_context Optional artifact provenance execution context.
#' @param gost_fn Injectable gprofiler2-compatible service function.
#' @param enricher_fn Injectable clusterProfiler-compatible enrichment function.
#'
#' @return S4 EnrichmentResults object containing enrichment data, plots, and summaries
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
                               execution_context = NULL,
                               gost_fn = gprofiler2::gost,
                               enricher_fn = clusterProfiler::enricher) {

  message("--- RUNNING processEnrichments VERSION [Timestamp: ", Sys.time(), "] ---")

  # Validate exclude_iea parameter
  if (is.null(exclude_iea)) {
    stop("exclude_iea must be explicitly set to TRUE or FALSE")
  }

  if (!is.logical(exclude_iea)) {
    stop("exclude_iea must be a logical value (TRUE or FALSE)")
  }

  organism <- resolveEnrichmentOrganism(taxon_id)
  is_supported <- organism$supported

  if(is_supported) {
    message(sprintf("Taxon ID %s found in supported organisms. Proceeding with gprofiler2 analysis...", taxon_id))

    species <- organism$species

    enrichment_results <- createEnrichmentResults(da_results@contrasts)

    # Process each contrast
    results <- da_results@da_data |>
      purrr::imap(function(de_data, contrast_key) {
        tryCatch({
        if(is.null(de_data)) {
          warning("No DE data found for contrast")
          return(NULL)
        }

        # Split data into up/down regulated
        message(sprintf("Total proteins before filtering: %d", nrow(de_data)))
        
        # [OK] DIAGNOSTIC: Check protein_id_column and data structure
        message(sprintf("DEBUG: protein_id_column = '%s'", protein_id_column))
        message(sprintf("DEBUG: Available columns in de_data: %s", paste(names(de_data), collapse = ", ")))
        
        if (!protein_id_column %in% names(de_data)) {
          stop(sprintf("ERROR: Column '%s' not found in DE data. Available columns: %s", 
                      protein_id_column, paste(names(de_data), collapse = ", ")))
        }
        
        # Check if the column has valid data
        protein_col_data <- de_data[[protein_id_column]]
        message(sprintf("DEBUG: First 3 values in protein_id_column: %s", 
                       paste(head(protein_col_data, 3), collapse = ", ")))
        message(sprintf("DEBUG: Column class: %s", class(protein_col_data)))

                 message("   processEnrichments Step: About to apply protein ID splitting and FDR filtering...")
         message(sprintf("      Data State (de_data Before Modify): Dims=%d rows, %d cols. First 5 IDs: %s", nrow(de_data), ncol(de_data), paste(head(de_data[[protein_id_column]], 5), collapse=", ")))
         
         # [OK] DEBUG66: Check NA values BEFORE filtering
         na_count_before <- sum(is.na(de_data[[protein_id_column]]))
         message(sprintf("      DEBUG66: NA values in %s column BEFORE filtering: %d out of %d (%.1f%%)", 
                        protein_id_column, na_count_before, nrow(de_data), (na_count_before/nrow(de_data))*100))
         
         # [OK] DEBUG66: Check how many pass q-value filter (before NA filter)
         passes_q_filter <- de_data |> dplyr::filter(.data$fdr_qvalue < q_cutoff)
         message(sprintf("      DEBUG66: Proteins passing q-value filter (< %.3f): %d", q_cutoff, nrow(passes_q_filter)))
         if (nrow(passes_q_filter) > 0) {
           na_in_sig <- sum(is.na(passes_q_filter[[protein_id_column]]))
           message(sprintf("      DEBUG66: NA values in %s among q-significant proteins: %d out of %d (%.1f%%)",
                          protein_id_column, na_in_sig, nrow(passes_q_filter), (na_in_sig/nrow(passes_q_filter))*100))
           message(sprintf("      DEBUG66: First 10 %s values in q-significant proteins: %s",
                          protein_id_column, paste(head(passes_q_filter[[protein_id_column]], 10), collapse = ", ")))
         }

         subset_sig <- de_data |>
            dplyr::mutate(
              !!rlang::sym(protein_id_column) := stringr::str_remove(.data[[protein_id_column]], ":.*")
            ) |>
           dplyr::filter(.data$fdr_qvalue < q_cutoff)
         
         # [OK] DEBUG66: Check if NAs are being filtered somewhere
         na_count_after_q <- sum(is.na(subset_sig[[protein_id_column]]))
         message(sprintf("      DEBUG66: NA values in %s AFTER q-filter (before any NA removal): %d out of %d",
                        protein_id_column, na_count_after_q, nrow(subset_sig)))

         message("   processEnrichments Step: Protein ID splitting and filtering complete.")
         message(sprintf("      Data State (subset_sig After Modify): Dims=%d rows, %d cols.", nrow(subset_sig), ncol(subset_sig)))
         if(nrow(subset_sig) > 0) {
            message(sprintf("      Data State (subset_sig After Modify): First 5 IDs: %s", paste(head(subset_sig[[protein_id_column]], 5), collapse=", ")))
            # Also check for the literal "Protein.Ids" to be sure
            if(any(subset_sig[[protein_id_column]] == "Protein.Ids", na.rm = TRUE)) {
              message("      WARNING: Literal 'Protein.Ids' found in the protein ID column after modification!")
            }
         }

        message(sprintf("Proteins passing q-value cutoff (%.3f): %d", q_cutoff, nrow(subset_sig)))
        
        # [OK] DEBUG 66: Extensive logging for subset_sig
        message("      Data State (subset_sig) Structure:")
        utils::str(subset_sig)
        if (nrow(subset_sig) > 0) {
          message("      Data State (subset_sig) Head:")
          print(head(subset_sig))
          message("      Data State (subset_sig fdr_qvalue range):")
          print(range(subset_sig$fdr_qvalue, na.rm = TRUE))
          message("      Data State (subset_sig log2FC range):")
          print(range(subset_sig$log2FC, na.rm = TRUE))
        }

        # No longer need subset_for_enrichment, as subset_sig is already filtered
        message(sprintf("Proteins available for enrichment: %d", nrow(subset_sig)))
        
        # [OK] DEBUG 66: Check subset_for_enrichment
        if (nrow(subset_sig) == 0) {
          message("      Data State (subset_sig): NO PROTEINS PASS ENRICHMENT CUTOFF!")
        }

        up_matrix <- subset_sig |>
          dplyr::filter(.data$log2FC > up_cutoff)
        message(sprintf("Up-regulated proteins (log2FC > %g): %d", up_cutoff, nrow(up_matrix)))
        
        # [OK] DEBUG66: Check NA values in up_matrix
        if (nrow(up_matrix) > 0) {
          na_in_up <- sum(is.na(up_matrix[[protein_id_column]]))
          message(sprintf("      DEBUG66: NA values in %s for UP-regulated proteins: %d out of %d (%.1f%%)",
                         protein_id_column, na_in_up, nrow(up_matrix), (na_in_up/nrow(up_matrix))*100))
          message(sprintf("      DEBUG66: First 10 UP gene names: %s",
                         paste(head(up_matrix[[protein_id_column]], 10), collapse = ", ")))
        }
        
        # [OK] DEBUG 66: Check up_matrix details
        if (nrow(up_matrix) > 0) {
          message("      Data State (up_matrix) Structure:")
          utils::str(up_matrix)
          message("      Data State (up_matrix) Head:")
          print(head(up_matrix))
          message("      Data State (up_matrix log2FC values):")
          print(summary(up_matrix$log2FC))
        } else {
          message("      Data State (up_matrix): NO UP-REGULATED PROTEINS FOUND!")
          message(sprintf("      Debug: up_cutoff = %g, max log2FC in subset_sig = %g", 
                         up_cutoff, if(nrow(subset_sig) > 0) max(subset_sig$log2FC, na.rm = TRUE) else NA))
        }

        down_matrix <- subset_sig |>
          dplyr::filter(.data$log2FC < -down_cutoff)
        message(sprintf("Down-regulated proteins (log2FC < -%g): %d", down_cutoff, nrow(down_matrix)))
        
        # [OK] DEBUG66: Check NA values in down_matrix
        if (nrow(down_matrix) > 0) {
          na_in_down <- sum(is.na(down_matrix[[protein_id_column]]))
          message(sprintf("      DEBUG66: NA values in %s for DOWN-regulated proteins: %d out of %d (%.1f%%)",
                         protein_id_column, na_in_down, nrow(down_matrix), (na_in_down/nrow(down_matrix))*100))
          message(sprintf("      DEBUG66: First 10 DOWN gene names: %s",
                         paste(head(down_matrix[[protein_id_column]], 10), collapse = ", ")))
        }
        
        # [OK] DEBUG 66: Check down_matrix details
        if (nrow(down_matrix) > 0) {
          message("      Data State (down_matrix) Structure:")
          utils::str(down_matrix)
          message("      Data State (down_matrix) Head:")
          print(head(down_matrix))
          message("      Data State (down_matrix log2FC values):")
          print(summary(down_matrix$log2FC))
        } else {
          message("      Data State (down_matrix): NO DOWN-REGULATED PROTEINS FOUND!")
          message(sprintf("      Debug: down_cutoff = %g, min log2FC in subset_sig = %g", 
                         down_cutoff, if(nrow(subset_sig) > 0) min(subset_sig$log2FC, na.rm = TRUE) else NA))
        }

        # Get background IDs from the full de_data for this contrast
        custom_bg <- de_data[[protein_id_column]] |>
          unique()

        message(sprintf("Using %d unique proteins as background for enrichment analysis", length(custom_bg)))
        
        # [OK] DEBUG 66: Check background details
        message("      Data State (custom_bg): First 5 background proteins:")
        message(paste(head(custom_bg, 5), collapse = ", "))

        run_gprofiler_direction <- function(data_subset, direction) {
          enrichment_args <- list(
            data_subset = data_subset,
            species = species,
            threshold = q_cutoff,
            sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"),
            domain_scope = "custom",
            custom_bg = custom_bg,
            exclude_iea = exclude_iea,
            protein_id_column = protein_id_column,
            correction_method = correction_method
          )
          if (isProtDiaEnrichExecutionContext(execution_context)) {
            enrichment_args$execution_context <- execution_context
            enrichment_args$request_context <- protDiaEnrichGprofilerContext(
              contrast_key,
              direction,
              protein_id_column
            )
            enrichment_args$gost_fn <- gost_fn
          }
          do.call(perform_enrichment, enrichment_args)
        }

        # Process up and down regulated genes
        list(
          up = tryCatch({
            if(nrow(up_matrix) > 0) {
              # [OK] DEBUG 66: Log call to perform_enrichment for up-regulated
              message("   processEnrichments Step: CALLING perform_enrichment for UP-REGULATED proteins...")
              message(sprintf("      About to analyze %d up-regulated proteins", nrow(up_matrix)))
              
              protein_col <- protein_id_column
              message(sprintf("      Using protein_id_column: %s", protein_col))
              message(sprintf("      Using species: %s", species))
              message(sprintf("      Using threshold: %s", q_cutoff))
              message(sprintf("      Using exclude_iea: %s", exclude_iea))
              
              up_result <- run_gprofiler_direction(up_matrix, "up")
              
              # [OK] DEBUG 66: Log result from perform_enrichment
              message("   processEnrichments Step: UP-REGULATED enrichment completed")
              message("      Data State (up_result): Checking perform_enrichment result...")
              
              if (is.null(up_result)) {
                message("      Data State (up_result): Result is NULL")
              } else {
                message("      Data State (up_result) Structure:")
                utils::str(up_result)
                
                if ("result" %in% names(up_result) && !is.null(up_result$result)) {
                  message(sprintf("      Data State (up_result$result): Found %d enriched terms", nrow(up_result$result)))
                  if (nrow(up_result$result) > 0) {
                    message("      Data State (up_result$result) Head:")
                    print(head(up_result$result))
                  }
                } else {
                  message("      Data State (up_result): No 'result' component or result is NULL")
                }
              }
              
              up_result
            } else {
              message("   processEnrichments Step: SKIPPING UP-REGULATED enrichment (no proteins)")
              if (isProtDiaEnrichExecutionContext(execution_context)) {
                protDiaEnrichRecordGprofilerSkip(
                  execution_context,
                  protDiaEnrichGprofilerContext(
                    contrast_key, "up", protein_id_column
                  ),
                  custom_bg, species, q_cutoff, correction_method, exclude_iea
                )
              }
              NULL
            }
          }, error = function(e) {
            message(sprintf("   processEnrichments Step: ERROR in up-regulated enrichment: %s", e$message))
            warning(sprintf("Error processing up-regulated genes: %s", e$message))
            message("Debug info for up-regulated genes:")
            message("Number of proteins: ", nrow(up_matrix))
            NULL
          }),

          down = tryCatch({
            if(nrow(down_matrix) > 0) {
              # [OK] DEBUG 66: Log call to perform_enrichment for down-regulated
              message("   processEnrichments Step: CALLING perform_enrichment for DOWN-REGULATED proteins...")
              message(sprintf("      About to analyze %d down-regulated proteins", nrow(down_matrix)))
              
              protein_col <- protein_id_column
              message(sprintf("      Using protein_id_column: %s", protein_col))
              message(sprintf("      Using species: %s", species))
              message(sprintf("      Using threshold: %s", q_cutoff))
              message(sprintf("      Using exclude_iea: %s", exclude_iea))
              
              down_result <- run_gprofiler_direction(down_matrix, "down")
              
              # [OK] DEBUG 66: Log result from perform_enrichment
              message("   processEnrichments Step: DOWN-REGULATED enrichment completed")
              message("      Data State (down_result): Checking perform_enrichment result...")
              
              if (is.null(down_result)) {
                message("      Data State (down_result): Result is NULL")
              } else {
                message("      Data State (down_result) Structure:")
                utils::str(down_result)
                
                if ("result" %in% names(down_result) && !is.null(down_result$result)) {
                  message(sprintf("      Data State (down_result$result): Found %d enriched terms", nrow(down_result$result)))
                  if (nrow(down_result$result) > 0) {
                    message("      Data State (down_result$result) Head:")
                    print(head(down_result$result))
                  }
                } else {
                  message("      Data State (down_result): No 'result' component or result is NULL")
                }
              }
              
              down_result
            } else {
              message("   processEnrichments Step: SKIPPING DOWN-REGULATED enrichment (no proteins)")
              if (isProtDiaEnrichExecutionContext(execution_context)) {
                protDiaEnrichRecordGprofilerSkip(
                  execution_context,
                  protDiaEnrichGprofilerContext(
                    contrast_key, "down", protein_id_column
                  ),
                  custom_bg, species, q_cutoff, correction_method, exclude_iea
                )
              }
              NULL
            }
          }, error = function(e) {
            message(sprintf("   processEnrichments Step: ERROR in down-regulated enrichment: %s", e$message))
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

    # Generate and store both static and interactive plots
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

    # Store static plots
    enrichment_results@enrichment_plots <- purrr::map(plot_results, function(x) {
      list(
        up = if(!is.null(x$up)) x$up$static else NULL,
        down = if(!is.null(x$down)) x$down$static else NULL
      )
    })

    # Store interactive plotly objects
    enrichment_results@enrichment_plotly <- purrr::map(plot_results, function(x) {
      list(
        up = if(!is.null(x$up)) x$up$interactive else NULL,
        down = if(!is.null(x$down)) x$down$interactive else NULL
      )
    })

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
       message("DEBUG: Mismatch detected or names missing in plot_results. Reassigning names.")
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
            htmlwidgets::saveWidget(
              plots[[direction]]$interactive,
              file_path,
              selfcontained = TRUE
            )
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

    go_mappings <- prepareClusterProfilerGoMappings(go_annotations)
    all_terms <- go_mappings$all_terms
    term2gene <- go_mappings$term2gene
    term2name <- go_mappings$term2name

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
    message("DEBUG: Using the following contrast names for processing and output:")
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

        execute_cluster_request <- function(genes, direction) {
          protDiaEnrichExecuteClusterRequest(
            execution_context, genes, background_IDs, short_name, direction,
            taxon_id, q_cutoff, exclude_iea, protein_id_column, term2gene,
            term2name, enricher_fn
          )
        }

        # Perform enrichment for up-regulated genes
        up_enrich <- tryCatch({
          execute_cluster_request(up_genes, "up")
        }, error = function(e) {
          warning(sprintf("Error processing up-regulated genes: %s", e$message))
          NULL
        })

        # Perform enrichment for down-regulated genes
        down_enrich <- tryCatch({
          execute_cluster_request(down_genes, "down")
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
        message("DEBUG: Assigned enrichment data to S4 object with short names.")
    } else {
        stop("Mismatch between number of processed results and expected contrast names.")
    }
    # Now enrichment_results@enrichment_data uses SHORT names

    # Create GO term mappings once (moved outside the plotting function)
    go_term_map <- dplyr::bind_rows(
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

                buildClusterProfilerEnrichmentPlots(
                  plot_data,
                  contrast,
                  direction,
                  q_cutoff
                )
            } else {
                NULL # Return NULL if no results
            }
        }) |> purrr::set_names(c("up", "down"))
        
        # Return plot components for this contrast
        direction_results
    }) |> purrr::set_names(contrast_names_to_use)
    message("--- Finished Plot Generation Loop ---")

    # --- Step 4: Store plots in S4 object using SHORT names ---
    enrichment_results@enrichment_plots <- purrr::map(plot_results_list, function(x) {
      list(up = x$up$static, down = x$down$static) # Handle NULLs gracefully if needed
    })
    enrichment_results@enrichment_plotly <- purrr::map(plot_results_list, function(x) {
      list(up = x$up$interactive, down = x$down$interactive) # Handle NULLs
    })
    message("DEBUG: Assigned plot data to S4 object slots.")
    message("DEBUG: Names in enrichment_results@enrichment_plotly:")
    print(names(enrichment_results@enrichment_plotly))

    # --- Step 5: Save interactive plots using SHORT names ---
    message("--- Starting Final Save Loop ---")
    purrr::walk(contrast_names_to_use, function(contrast) {
        # Access plot data using SHORT name (now stored with short names)
        plots <- enrichment_results@enrichment_plotly[[contrast]]

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
            if (!is.null(plots[[direction]])) { # Access up/down list
                interactive_plot <- plots[[direction]] # The actual plotly object
                static_plot <- enrichment_results@enrichment_plots[[contrast]][[direction]] # The static ggplot

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
                  htmlwidgets::saveWidget(
                    interactive_plot, # Pass the plot object directly
                    file_path,
                    selfcontained = TRUE
                  )
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
  if(!contrast %in% names(enrichment_results@enrichment_plotly)) {
    stop("Contrast not found")
  }
  if(!direction %in% c("up", "down")) {
    stop("Direction must be 'up' or 'down'")
  }
  enrichment_results@enrichment_plotly[[contrast]][[direction]]
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
