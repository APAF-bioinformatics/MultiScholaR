# Enrichment constructor functions removed - migrated to R/func_general_s4_objects.R
# ----------------------------------------------------------------------------

resolveEnrichmentOrganism <- function(taxon_id) {
  supported <- tibble::tribble(
    ~taxid, ~id, ~name,
    "9606", "hsapiens", "Homo sapiens",
    "10090", "mmusculus", "Mus musculus",
    "10116", "rnorvegicus", "Rattus norvegicus",
    "7227", "dmelanogaster", "Drosophila melanogaster",
    "6239", "celegans", "Caenorhabditis elegans",
    "4932", "scerevisiae", "Saccharomyces cerevisiae",
    "3702", "athaliana", "Arabidopsis thaliana",
    "7955", "drerio", "Danio rerio",
    "9031", "ggallus", "Gallus gallus",
    "9823", "sscrofa", "Sus scrofa",
    "9913", "btaurus", "Bos taurus",
    "9544", "mmulatta", "Macaca mulatta",
    "9598", "ptroglodytes", "Pan troglodytes"
  )
  match <- supported |> dplyr::filter(.data$taxid == as.character(taxon_id))
  list(
    supported = nrow(match) == 1L,
    species = if (nrow(match) == 1L) match$id[[1L]] else NULL
  )
}

prepareClusterProfilerGoMappings <- function(go_annotations) {
  term_table <- function(id_column, ontology) {
    go_annotations |>
      dplyr::filter(!is.na(.data[[id_column]])) |>
      tidyr::separate_rows(dplyr::all_of(id_column), sep = "; ") |>
      dplyr::select(.data$Entry, dplyr::all_of(id_column)) |>
      dplyr::rename(TERM = dplyr::all_of(id_column)) |>
      dplyr::mutate(ONTOLOGY = ontology)
  }
  all_terms <- rbind(
    term_table("go_id_go_biological_process", "BP"),
    term_table("go_id_go_molecular_function", "MF"),
    term_table("go_id_go_cellular_compartment", "CC")
  )
  term2gene <- all_terms |>
    dplyr::select(.data$TERM, .data$Entry) |>
    dplyr::distinct()
  term2name <- data.frame(
    TERM = unique(all_terms$TERM),
    NAME = purrr::map_chr(unique(all_terms$TERM), function(term) {
      tryCatch(AnnotationDbi::Term(GO.db::GOTERM[[term]]), error = function(e) term)
    })
  )
  list(all_terms = all_terms, term2gene = term2gene, term2name = term2name)
}



# ----------------------------------------------------------------------------
# perform_enrichment
# ----------------------------------------------------------------------------
#' Run gprofiler2 Enrichment for One Direction
#'
#' @param data_subset Differential-abundance rows selected for one direction.
#' @param species gprofiler2 organism identifier.
#' @param threshold Adjusted significance threshold sent to gprofiler2.
#' @param sources Annotation sources sent to gprofiler2.
#' @param domain_scope gprofiler2 background-domain mode.
#' @param custom_bg Identifier background sent to gprofiler2.
#' @param exclude_iea Whether to exclude electronically inferred annotations.
#' @param max_retries Maximum number of service attempts.
#' @param wait_time Seconds between service attempts.
#' @param protein_id_column Identifier column in `data_subset`.
#' @param correction_method Multiple-testing correction method.
#' @param execution_context Optional artifact provenance execution context.
#' @param request_context Optional contrast, direction, and mapping provenance.
#' @param gost_fn Injectable gprofiler2-compatible service function.
#'
#' @return A gprofiler2 response list, or `NULL` for no result in memory mode.
#' @export
perform_enrichment <- function(data_subset,
                               species,
                               threshold,
                               sources,
                               domain_scope,
                               custom_bg,
                               exclude_iea = FALSE,
                               max_retries = 5,
                               wait_time = 5,
                               protein_id_column,
                               correction_method = "gSCS",
                               execution_context = NULL,
                               request_context = list(),
                               gost_fn = gprofiler2::gost) {
  
  message("--- Entering perform_enrichment ---")
  message(sprintf("   perform_enrichment Arg: species = %s", species))
  message(sprintf("   perform_enrichment Arg: threshold = %s", threshold))
  message(sprintf("   perform_enrichment Arg: sources = %s", paste(sources, collapse = ", ")))
  message(sprintf("   perform_enrichment Arg: domain_scope = %s", domain_scope))
  message(sprintf("   perform_enrichment Arg: exclude_iea = %s", exclude_iea))
  message(sprintf("   perform_enrichment Arg: protein_id_column = %s", protein_id_column))
  message(sprintf("   perform_enrichment Arg: max_retries = %s", max_retries))
  
  message("      Data State (data_subset): Checking input data...")
  message(sprintf("      Data State (data_subset): Dims = %d rows, %d cols", nrow(data_subset), ncol(data_subset)))
  message("      Data State (data_subset) Structure:")
  utils::str(data_subset)
  message("      Data State (data_subset) Head:")
  print(head(data_subset))
  
  if (nrow(data_subset) == 0) {
    if (isProtDiaEnrichExecutionContext(execution_context)) {
      request <- protDiaEnrichCanonicalRequest(
        backend = "gprofiler2",
        contrast = request_context$contrast,
        direction = request_context$direction,
        identifiers = character(),
        background = custom_bg,
        parameters = list(
          organism = species,
          ordered_query = FALSE,
          sources = as.list(sources),
          user_threshold = as.double(threshold),
          correction_method = correction_method,
          exclude_iea = isTRUE(exclude_iea),
          evcodes = TRUE,
          domain_scope = domain_scope,
          significant = TRUE,
          highlight = TRUE
        ),
        mapping = request_context$mapping
      )
      protDiaEnrichRecordSkippedRequest(execution_context, request)
    }
    message("   perform_enrichment Step: Data subset is empty, returning NULL")
    message("--- Exiting perform_enrichment. Returning: NULL ---")
    return(NULL)
  }

  message("   perform_enrichment Step: Extracting protein IDs from data...")
  # Clean data before enrichment
  protein_ids <- data_subset[[protein_id_column]]
  
  message(sprintf("      Data State (protein_ids): Length = %d", length(protein_ids)))
  message(sprintf("      Data State (protein_ids): Class = %s", class(protein_ids)))
  message(sprintf("      Data State (protein_ids): First 5 values = %s", paste(head(protein_ids, 5), collapse = ", ")))
  
  # [OK] DEBUG66: Comprehensive NA analysis
  na_count_initial <- sum(is.na(protein_ids))
  empty_count_initial <- sum(protein_ids == "", na.rm = TRUE)
  valid_count_initial <- length(protein_ids) - na_count_initial - empty_count_initial
  
  message(sprintf("      +==================================================================+"))
  message(sprintf("      | DEBUG66: PROTEIN ID AVAILABILITY ANALYSIS                       |"))
  message(sprintf("      +==================================================================+"))
  message(sprintf("      | Total proteins in data_subset:     %5d (100.0%%)                 |", length(protein_ids)))
  message(sprintf("      | NA values in %s column:   %5d (%5.1f%%)                 |", 
                 protein_id_column, na_count_initial, (na_count_initial/length(protein_ids))*100))
  message(sprintf("      | Empty string values:                %5d (%5.1f%%)                 |", 
                 empty_count_initial, (empty_count_initial/length(protein_ids))*100))
  message(sprintf("      | Valid (non-NA, non-empty) values:  %5d (%5.1f%%)                 |", 
                 valid_count_initial, (valid_count_initial/length(protein_ids))*100))
  message(sprintf("      +==================================================================+"))

  if (any(is.na(protein_ids))) {
    na_count <- sum(is.na(protein_ids))
    message(sprintf("   [WARNING]  WARNING: Filtering out %d proteins with NA %s values!", na_count, protein_id_column))
    message(sprintf("   [WARNING]  %d proteins will be EXCLUDED from enrichment analysis!", na_count))
    
    if (na_count > (length(protein_ids) * 0.5)) {
      message("   +===========================================================================+")
      message("   | [WARNING]  CRITICAL WARNING: > 50% of proteins have NA gene names!              |")
      message("   | Consider using 'Protein.Ids' instead of 'gene_name' for enrichment!      |")
      message("   | OR ensure gene names are properly annotated in your data.                |")
      message("   +===========================================================================+")
    }
    
    warning(paste("NA values found in", protein_id_column, "column"))
    data_subset <- data_subset |> dplyr::filter(!is.na(.data[[protein_id_column]]))
    protein_ids <- data_subset[[protein_id_column]]
    message(sprintf("      Data State (protein_ids after NA filter): Length = %d", length(protein_ids)))
  }

  protein_ids <- unique(protein_ids[!is.na(protein_ids) & protein_ids != ""])
  
  message(sprintf("      DEBUG66: Final protein IDs to submit to gprofiler2: %d unique IDs", length(protein_ids)))
  if (length(protein_ids) > 0 && length(protein_ids) <= 20) {
    message(sprintf("      DEBUG66: All protein IDs being submitted: %s", paste(protein_ids, collapse = ", ")))
  } else if (length(protein_ids) > 20) {
    message(sprintf("      DEBUG66: First 20 protein IDs being submitted: %s", paste(head(protein_ids, 20), collapse = ", ")))
  }

  message("   perform_enrichment Step: Checking custom background...")
  message(sprintf("      Data State (custom_bg): Length = %d", length(custom_bg)))
  message(sprintf("      Data State (custom_bg): Class = %s", class(custom_bg)))
  message(sprintf("      Data State (custom_bg): First 5 values = %s", paste(head(custom_bg, 5), collapse = ", ")))

  if (any(is.na(custom_bg))) {
    na_bg_count <- sum(is.na(custom_bg))
    message(sprintf("   perform_enrichment Step: Found %d NA values in custom background, filtering...", na_bg_count))
    warning("NA values found in custom background IDs")
    custom_bg <- custom_bg[!is.na(custom_bg)]
    message(sprintf("      Data State (custom_bg after NA filter): Length = %d", length(custom_bg)))
  }

  custom_bg <- unique(custom_bg)

  if (isProtDiaEnrichExecutionContext(execution_context)) {
    request <- protDiaEnrichCanonicalRequest(
      backend = "gprofiler2",
      contrast = request_context$contrast,
      direction = request_context$direction,
      identifiers = protein_ids,
      background = custom_bg,
      parameters = list(
        organism = species,
        ordered_query = FALSE,
        sources = as.list(sources),
        user_threshold = as.double(threshold),
        correction_method = correction_method,
        exclude_iea = isTRUE(exclude_iea),
        evcodes = TRUE,
        domain_scope = domain_scope,
        significant = TRUE,
        highlight = TRUE
      ),
      mapping = request_context$mapping
    )
    if (length(protein_ids) == 0L) {
      protDiaEnrichRecordSkippedRequest(execution_context, request)
      return(NULL)
    }
    return(protDiaEnrichExecuteRequest(
      execution_context,
      request,
      call = function() {
        gost_fn(
          query = request$identifiers,
          organism = species,
          ordered_query = FALSE,
          sources = sources,
          user_threshold = threshold,
          correction_method = correction_method,
          exclude_iea = exclude_iea,
          evcodes = TRUE,
          domain_scope = domain_scope,
          custom_bg = request$background,
          significant = TRUE,
          highlight = TRUE
        )
      },
      max_retries = max_retries,
      wait_time = wait_time
    ))
  }

  result <- NULL
  attempt <- 1

  message(sprintf("   perform_enrichment Step: Starting gprofiler2 gost() retry loop (max %d attempts)...", max_retries))

  while (is.null(result) && attempt <= max_retries) {
    message(sprintf("   perform_enrichment Attempt %d: Calling gprofiler2::gost()...", attempt))
    
    tryCatch({
      message("   perform_enrichment Step: About to call gprofiler2::gost() with parameters:")
      message(sprintf("      gost query: %d protein IDs", length(protein_ids)))
      message(sprintf("      gost organism: %s", species))
      message(sprintf("      gost sources: %s", paste(sources, collapse = ", ")))
      message(sprintf("      gost user_threshold: %s", threshold))
      message(sprintf("      gost domain_scope: %s", domain_scope))
      message(sprintf("      gost custom_bg: %d background IDs", length(custom_bg)))
      message(sprintf("      gost exclude_iea: %s", exclude_iea))
      
      result <- gost_fn(
        query = protein_ids,
        organism = species,
        ordered_query = FALSE,
        sources = sources,
        user_threshold = threshold,
        correction_method = correction_method,
        exclude_iea = exclude_iea,
        evcodes = TRUE,
        domain_scope = domain_scope,
        custom_bg = custom_bg,
        significant = TRUE,
        highlight = TRUE
      )
      
      message("   perform_enrichment Step: gprofiler2::gost() completed successfully")
      message("      Data State (gost result): Checking gost result structure...")
      
      if (is.null(result)) {
        message("      Data State (gost result): Result is NULL")
      } else {
        message("      Data State (gost result) Structure:")
        utils::str(result)
        
        if ("result" %in% names(result)) {
          if (is.null(result$result)) {
            message("      Data State (gost result$result): Is NULL")
          } else {
            message(sprintf("      Data State (gost result$result): Dims = %d rows, %d cols", nrow(result$result), ncol(result$result)))
            message("      Data State (gost result$result) Column names:")
            print(names(result$result))
            if (nrow(result$result) > 0) {
              message("      Data State (gost result$result) Head:")
              print(head(result$result))
            } else {
              message("      Data State (gost result$result): NO ROWS - Empty results!")
            }
          }
        } else {
          message("      Data State (gost result): No 'result' component found")
        }
      }
      
    }, error = function(e) {
      message(sprintf("   perform_enrichment Attempt %d: ERROR in gprofiler2::gost(): %s", attempt, e$message))
      message(sprintf("   perform_enrichment Attempt %d: Will wait %d seconds before retry...", attempt, wait_time))
      Sys.sleep(wait_time)
    })
    
    attempt <- attempt + 1
  }

  if (is.null(result)) {
    message("   perform_enrichment Step: All retry attempts failed, returning NULL")
    message("--- Exiting perform_enrichment. Returning: NULL ---")
    return(NULL)
  }

  message("   perform_enrichment Step: Successfully obtained gprofiler2 result")
  message("--- Exiting perform_enrichment. Returning: gost result object ---")
  return(result)
}


# ----------------------------------------------------------------------------
# generate_enrichment_plots
# ----------------------------------------------------------------------------
buildGprofilerEnrichmentPlots <- function(enrichment_result, contrast, direction) {
  empty <- is.null(enrichment_result) ||
    is.null(enrichment_result$result) ||
    nrow(enrichment_result$result) == 0L
  if (empty) {
    return(list(static = NULL, interactive = NULL))
  }

  # Extract the significance threshold from the result object metadata
  significance_threshold <- enrichment_result$meta$query_metadata$user_threshold

  # Prepare data for plotting
  plot_data <- enrichment_result$result |>
    dplyr::mutate(
      neg_log10_p = -log10(.data$p_value),
      # Ensure 'source' is a factor with a consistent level order for plotting
      source = factor(
        .data$source,
        levels = c("GO:BP", "GO:CC", "GO:MF", "KEGG", "REAC")
      )
    ) |>
    # Drop any levels that are not actually present in the data to avoid empty spaces on the plot
    dplyr::mutate(source = forcats::fct_drop(.data$source))

  # Identify the top significant terms to add labels for
  top_terms <- plot_data |>
    dplyr::group_by(.data$source) |>
    dplyr::arrange(.data$p_value) |>
    dplyr::slice_head(n = 3) |>
    dplyr::ungroup()

  # Create the static ggplot object
  static_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data$source,
      y = .data$neg_log10_p,
      text = paste0(
        "Term: ", .data$term_name, "\n",
        "ID: ", .data$term_id, "\n",
        "P-value: ", signif(.data$p_value, 3), "\n",
        "Genes: ", .data$intersection_size
      )
    )
  ) +
    # Add a line for the significance threshold
    ggplot2::geom_hline(
      yintercept = -log10(significance_threshold),
      linetype = "dashed",
      color = "red"
    ) +
    ggplot2::geom_jitter(
      ggplot2::aes(color = .data$source, size = .data$term_size),
      alpha = 0.7,
      width = 0.2
    ) +
    ggrepel::geom_text_repel(
      data = top_terms,
      ggplot2::aes(label = .data$term_name),
      size = 3, max.overlaps = 15, box.padding = 0.5, point.padding = 0.3, force = 5
    ) +
    ggplot2::scale_color_manual(
      values = c(`GO:BP` = "#ff9900", `GO:CC` = "#109618", `GO:MF` = "#dc3912",
                 KEGG = "#dd4477", REAC = "#3366cc"),
      name = "Source", drop = FALSE
    ) +
    ggplot2::scale_size_continuous(name = "Term Size", range = c(3, 10)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10),
      axis.title = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      legend.position = "right"
    ) +
    ggplot2::labs(
      title = paste0("Enrichment for ", contrast, " (", direction, "-regulated)"),
      x = "Annotation Source", y = "-log10(p-value)"
    ) +
    ggplot2::guides(color = "none")

  # Convert to an interactive plotly object
  interactive_plot <- tryCatch(
    plotly::ggplotly(static_plot, tooltip = "text"),
    error = function(error) {
      warning(paste(
        "Failed to convert ggplot to plotly for",
        contrast,
        direction,
        ":",
        conditionMessage(error)
      ))
      return(NULL)
    }
  )

  list(
    static = static_plot,
    interactive = interactive_plot
  )
}

# Plot generation function
#' @export
generate_enrichment_plots <- function(enrichment_result, contrast, direction, pathway_dir) {
  plots <- buildGprofilerEnrichmentPlots(
    enrichment_result,
    contrast,
    direction
  )
  if (is.null(plots$static)) {
    return(plots)
  }

  # Save the results data table
  tryCatch({
    result_table <- enrichment_result$result
    result_table$parents <- sapply(result_table$parents, paste, collapse = ", ")
    readr::write_tsv(result_table,
      file = file.path(pathway_dir, paste0(contrast, "_", direction, "_enrichment_results.tsv"))
    )
  }, error = function(e) {
    warning(paste(
      "Failed to write enrichment results table for",
      contrast,
      direction,
      ":",
      e$message
    ))
  })

  plots
}

buildClusterProfilerEnrichmentPlots <- function(plot_data,
                                                contrast,
                                                direction,
                                                q_cutoff) {
  if (is.null(plot_data) || nrow(plot_data) == 0L) {
    return(list(static = NULL, interactive = NULL))
  }
  plot_data$source <- factor(
    plot_data$source,
    levels = c("GO:BP", "GO:CC", "GO:MF", "Other")
  )
  top_terms <- plot_data |>
    dplyr::filter(.data$qvalue < q_cutoff) |>
    dplyr::arrange(.data$qvalue) |>
    dplyr::slice_head(n = 5)
  static_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data$source,
      y = .data$neg_log10_q,
      text = paste0(
        "Term: ", .data$term, "\n",
        "ID: ", .data$ID, "\n",
        "Genes: ", .data$Count, "\n",
        "Gene Ratio: ", .data$GeneRatio, "\n",
        "Background Ratio: ", .data$BgRatio, "\n",
        "Q-value: ", signif(.data$qvalue, 3)
      )
    )
  ) +
    ggplot2::geom_hline(
      yintercept = -log10(q_cutoff),
      linetype = "dashed",
      color = "darkgrey"
    ) +
    ggplot2::geom_jitter(
      ggplot2::aes(
        size = .data$gene_count,
        color = -log10(.data$qvalue)
      ),
      alpha = 0.7,
      width = 0.2
    ) +
    ggrepel::geom_text_repel(
      data = top_terms,
      ggplot2::aes(label = .data$term),
      size = 3,
      max.overlaps = 15,
      box.padding = 0.5,
      point.padding = 0.3,
      force = 5
    ) +
    ggplot2::scale_color_gradient(
      low = "#FED976",
      high = "#800026",
      name = "-log10(q-value)"
    ) +
    ggplot2::scale_size_continuous(name = "Gene Count", range = c(3, 12)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 10, angle = 0),
      axis.text.y = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 8)
    ) +
    ggplot2::labs(
      title = paste0(
        contrast,
        " ",
        tools::toTitleCase(direction),
        "-regulated"
      ),
      x = "GO Category",
      y = "-log10(q-value)"
    )
  list(
    static = static_plot,
    interactive = plotly::ggplotly(static_plot, tooltip = "text")
  )
}


# ----------------------------------------------------------------------------
# summarize_enrichment
# ----------------------------------------------------------------------------
# Summary function
#' @export
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


# ----------------------------------------------------------------------------
