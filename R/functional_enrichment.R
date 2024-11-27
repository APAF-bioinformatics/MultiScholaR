#' Create DE Results for Enrichment
' 
' @param contrasts_tbl A tibble containing contrast information
' @param design_matrix A data frame containing the design matrix
' @param de_output_dir Directory containing DE results files
' @return An S4 object of class de_results_for_enrichment
' @export
reateDEResultsForEnrichment <- function(contrasts_tbl, design_matrix, de_output_dir) {
   # Define the S4 class if it doesn't exist
   if (!isClass("de_results_for_enrichment")) {
       setClass("de_results_for_enrichment",
               slots = list(
                   contrasts = "tbl_df",
                   de_data = "list",
                   design_matrix = "data.frame"
               ))
   }
   
   # Helper function to format contrast filenames
   format_contrast_filename <- function(contrast_string) {
       contrast_name <- stringr::str_split(contrast_string, "=")[[1]][1] |>
           stringr::str_replace_all("\\.", "_")
       
       paste0("de_proteins_", contrast_name, "_long_annot.tsv")
   }
   
   # Create and populate the S4 object
   de_results <- new("de_results_for_enrichment")
   de_results@contrasts <- contrasts_tbl
   de_results@design_matrix <- design_matrix
   de_results@de_data <- contrasts_tbl$contrasts |>
       purrr::set_names() |>
       purrr::map(function(contrast) {
           filename <- format_contrast_filename(contrast)
           filepath <- file.path(de_output_dir, filename)
           
           if (!file.exists(filepath)) {
               warning("File not found: ", filepath)
               return(NULL)
           }
           
           readr::read_tsv(filepath, show_col_types = FALSE)
       })
   
   return(de_results)


# S4 class definition
setClass("EnrichmentResults",
         slots = list(
           contrasts = "tbl_df",
           enrichment_data = "list",
           enrichment_plots = "list",        # gostplot objects
           enrichment_plotly = "list",       # interactive plotly objects
           enrichment_summaries = "list"
         ))

# Constructor function
createEnrichmentResults <- function(contrasts_tbl) {
  new("EnrichmentResults",
      contrasts = contrasts_tbl,
      enrichment_data = list(),
      enrichment_plots = list(),
      enrichment_plotly = list(),
      enrichment_summaries = list())
}

perform_enrichment <- function(data_subset, species, threshold, sources, domain_scope, background_IDs, max_retries = 5, wait_time = 5) {
  if (nrow(data_subset) == 0) {
    return(NULL)
  }
  
  # Clean data before enrichment
  if (any(is.na(data_subset$uniprot_acc))) {
    warning("NA values found in uniprot_acc column")
    data_subset <- data_subset |> filter(!is.na(uniprot_acc))
  }

  if (any(is.na(background_IDs))) {
    warning("NA values found in background IDs")
    background_IDs <- background_IDs[!is.na(background_IDs)]
  }
  
  result <- NULL
  attempt <- 1
  
  while (is.null(result) && attempt <= max_retries) {
    tryCatch({
      result <- gprofiler2::gost(
        query = data_subset$uniprot_acc,
        organism = species,
        ordered_query = FALSE,
        sources = sources,
        user_threshold = threshold,
        domain_scope = domain_scope,
        correction_method = "gSCS",
        evcodes = TRUE,
        custom_bg = background_IDs,
        significant = FALSE
      )
      
      # If no significant results, return NULL immediately without retrying
      if (is.null(result$result)) {
        message("No significant results found. Moving to next analysis.")
        return(NULL)
      }
      
    }, error = function(e) {
      # Retry only for connection/timeout errors
      if (grepl("408", e$message) || 
          grepl("Could not resolve host", e$message) ||
          grepl("Connection refused", e$message)) {
        message(paste("Attempt", attempt, "failed with connection error. Retrying in", wait_time, "seconds..."))
        Sys.sleep(wait_time)
        result <- NULL  # Ensure the loop continues
      } else if (grepl("Please make sure that the organism is correct", e$message)) {
        warning("Organism check failed. Please verify species identifier.")
        return(NULL)
      } else {
        warning(paste("Error in enrichment analysis:", e$message))
        return(NULL)
      }
    })
    
    attempt <- attempt + 1
  }
  
  if (is.null(result) && attempt > max_retries) {
    message("Failed to get a valid response after ", max_retries, " attempts. Returning NULL.")
    return(NULL)
  }
  
  return(result)
}

# Plot generation function
generate_enrichment_plots <- function(enrichment_result, contrast, direction, pathway_dir) {
  if (is.null(enrichment_result) || nrow(enrichment_result$result) == 0) {
    return(list(
      static = NULL,
      interactive = NULL
    ))
  }
  
  # Generate static gostplot
  static_plot <- gostplot(
    enrichment_result,
    capped = TRUE,
    interactive = FALSE,
    pal = c(`GO:MF` = "#dc3912", `GO:BP` = "#ff9900", `GO:CC` = "#109618", 
            KEGG = "#dd4477", REAC = "#3366cc")
  )
  
  # Convert to plotly
  interactive_plot <- plotly::ggplotly(static_plot)
  
  # Save interactive plot
  plot_file <- file.path(pathway_dir, 
                        paste0(contrast, "_", direction, "_enrichment_plot.html"), 
                        fsep = "/")
  lib_dir <- file.path(pathway_dir, 
                      paste0(contrast, "_", direction, "_libs"), 
                      fsep = "/")
  dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)
  
  tryCatch({
    htmlwidgets::saveWidget(interactive_plot, 
                          file = plot_file, 
                          selfcontained = FALSE, 
                          libdir = lib_dir)
  }, error = function(e) {
    warning(sprintf("Error saving plot for %s_%s: %s", contrast, direction, e$message))
  })
  
  # Save results table
  result_table <- enrichment_result$result
  result_table$parents <- sapply(result_table$parents, paste, collapse = ", ")
  write.table(result_table, 
              file = file.path(pathway_dir, 
                             paste0(contrast, "_", direction, "_enrichment_results.tsv")),
              sep = "\t", 
              row.names = FALSE, 
              quote = FALSE)
  
  return(list(
    static = static_plot,
    interactive = interactive_plot
  ))
}

# Summary function
summarize_enrichment <- function(enrichment_result) {
  if (is.null(enrichment_result) || length(enrichment_result$result) == 0) {
    return(data.frame(
      total = 0,
      GO_BP = 0,
      GO_CC = 0,
      GO_MF = 0,
      KEGG = 0,
      REAC = 0
    ))
  }
  
  result_df <- enrichment_result$result
  
  data.frame(
    total = nrow(result_df),
    GO_BP = sum(result_df$source == "GO:BP"),
    GO_CC = sum(result_df$source == "GO:CC"),
    GO_MF = sum(result_df$source == "GO:MF"),
    KEGG = sum(result_df$source == "KEGG"),
    REAC = sum(result_df$source == "REAC")
  )
}

processEnrichments <- function(de_results, 
                             taxon_id,           
                             up_cutoff = 0, 
                             down_cutoff = 0, 
                             q_cutoff = 0.05, 
                             pathway_dir) {
  
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
      filter(taxid == as.character(taxon_id)) |>
      pull(id)
    
    enrichment_results <- createEnrichmentResults(de_results@contrasts)
    
    # Process each contrast
    results <- de_results@de_data |>
      purrr::map(function(de_data) {
        if(is.null(de_data)) {
          warning("No DE data found for contrast")
          return(NULL)
        }
        
        # Split data into up/down regulated
        subset_sig <- de_data |>
          filter(fdr_qvalue < q_cutoff)
        
        up_matrix <- subset_sig |>
          filter(log2FC > up_cutoff)
        
        down_matrix <- subset_sig |>
          filter(log2FC < -down_cutoff)
        
        # Get background IDs - using uniprot_acc instead of best_uniprot_acc
        background_IDs <- unique(de_data$uniprot_acc)
        
        # Process up and down regulated genes
        list(
          up = tryCatch({
            if(nrow(up_matrix) > 0) {
              perform_enrichment(
                data_subset = up_matrix,
                species = species,
                threshold = q_cutoff,
                sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"),
                domain_scope = "custom",
                background_IDs = background_IDs
              )
            } else NULL
          }, error = function(e) {
            warning(sprintf("Error processing up-regulated genes: %s", e$message))
            NULL
          }),
          
          down = tryCatch({
            if(nrow(down_matrix) > 0) {
              perform_enrichment(
                data_subset = down_matrix,
                species = species,
                threshold = q_cutoff,
                sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"),
                domain_scope = "custom",
                background_IDs = background_IDs
              )
            } else NULL
          }, error = function(e) {
            warning(sprintf("Error processing down-regulated genes: %s", e$message))
            NULL
          })
        )
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
    
    return(enrichment_results)
    
  } else {
    stop(sprintf("Taxon ID %s not found in gprofiler2 supported organisms. 
                 ClusterProfiler implementation pending. 
                 Current supported organisms in gprofiler2: %s", 
                 taxon_id,
                 paste(unique(supported_organisms$taxon_id), collapse = ", ")))
  }
}

# Helper function to access results
getEnrichmentResult <- function(enrichment_results, contrast, direction) {
  if(!contrast %in% names(enrichment_results@enrichment_data)) {
    stop("Contrast not found")
  }
  if(!direction %in% c("up", "down")) {
    stop("Direction must be 'up' or 'down'")
  }
  enrichment_results@enrichment_data[[contrast]][[direction]]
}

# Helper function to access plotly objects
getEnrichmentPlotly <- function(enrichment_results, contrast, direction) {
  if(!contrast %in% names(enrichment_results@enrichment_plotly)) {
    stop("Contrast not found")
  }
  if(!direction %in% c("up", "down")) {
    stop("Direction must be 'up' or 'down'")
  }
  enrichment_results@enrichment_plotly[[contrast]][[direction]]
}

# Helper function to get summary
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
