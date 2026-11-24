# Enrichment constructor functions removed - migrated to R/func_general_s4_objects.R
# ----------------------------------------------------------------------------



# ----------------------------------------------------------------------------
# perform_enrichment
# ----------------------------------------------------------------------------
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
                               correction_method = "gSCS") {
  
  
  
  if (nrow(data_subset) == 0) {
    return(NULL)
  }

  # Clean data before enrichment
  protein_ids <- data_subset[[protein_id_column]]
  
  
  na_count_initial <- sum(is.na(protein_ids))
  empty_count_initial <- sum(protein_ids == "", na.rm = TRUE)
  valid_count_initial <- length(protein_ids) - na_count_initial - empty_count_initial
  

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
  }

  protein_ids <- unique(protein_ids[!is.na(protein_ids) & protein_ids != ""])
  


  if (any(is.na(custom_bg))) {
    na_bg_count <- sum(is.na(custom_bg))
    warning("NA values found in custom background IDs")
    custom_bg <- custom_bg[!is.na(custom_bg)]
  }

  custom_bg <- unique(custom_bg)

  result <- NULL
  attempt <- 1
  # GUI audit R3, A107-04: whether g:Profiler answered any attempt, and the last error when it did not. `gost` returns
  # NULL for "no significant terms", so a NULL result alone does not say which of the two happened.
  answered <- FALSE
  last_error <- NULL


  while (is.null(result) && attempt <= max_retries) {
    message(sprintf("   perform_enrichment Attempt %d: Calling gprofiler2::gost()...", attempt))
    
    tryCatch({
      message(sprintf("      gost query: %d protein IDs", length(protein_ids)))
      message(sprintf("      gost organism: %s", species))
      message(sprintf("      gost sources: %s", paste(sources, collapse = ", ")))
      message(sprintf("      gost user_threshold: %s", threshold))
      message(sprintf("      gost domain_scope: %s", domain_scope))
      message(sprintf("      gost custom_bg: %d background IDs", length(custom_bg)))
      message(sprintf("      gost exclude_iea: %s", exclude_iea))
      
      result <- gprofiler2::gost(
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
      answered <- TRUE
      
      
    }, error = function(e) {
      last_error <<- e$message
      message(sprintf("   perform_enrichment Attempt %d: ERROR in gprofiler2::gost(): %s", attempt, e$message))
      message(sprintf("   perform_enrichment Attempt %d: Will wait %d seconds before retry...", attempt, wait_time))
      Sys.sleep(wait_time)
    })
    
    attempt <- attempt + 1
  }

  if (is.null(result)) {
    # A107-04: the return value stays NULL for every caller. When no attempt was answered, a classed condition is
    # signalled first; it does nothing unless a caller handles `multischolar_enrichment_service_failure`
    # (`diaEnrichmentOperation`'s kernel does, and fails the run instead of saving "no terms").
    if (!answered) {
      signalCondition(structure(class = c("multischolar_enrichment_service_failure", "condition"),
        list(message = sprintf("g:Profiler did not answer after %d attempts: %s", max_retries,
          last_error %||% "no attempt was made"), call = NULL)))
    }
    return(NULL)
  }

  return(result)
}


# ----------------------------------------------------------------------------
# Enrichment plots, built from the stored result (KIT-119)
# ----------------------------------------------------------------------------
# The seed of the plots' point jitter and label placement: two builds from the same result are the same plot.
ENRICHMENT_PLOT_SEED <- 119L

#' The static enrichment plot of one contrast and direction, from its stored result (KIT-119)
#'
#' The one place the enrichment plots are drawn. The run writes its HTML and PNG files from it and the tab draws
#' from it at render time, so the saved file and the tab show the same plot, and no plot is kept in the
#' `EnrichmentResults` plot slots. The jitter and the label placement are seeded (`ENRICHMENT_PLOT_SEED`).
#' @param enrichment_result One entry of `EnrichmentResults@enrichment_data[[contrast]][[direction]]`: a gprofiler2
#'   `gost` result list (`result`, `meta`), or a clusterProfiler `enrichResult` (see `setEnrichmentPlotInputs`).
#' @param contrast Contrast name, for the title.
#' @param direction `"up"` or `"down"`, for the title.
#' @return A ggplot, or NULL when the result is NULL or has no terms.
#' @noRd
buildEnrichmentGgplot <- function(enrichment_result, contrast, direction) {
  if (isS4(enrichment_result)) {
    return(buildClusterProfilerEnrichmentGgplot(enrichment_result, contrast, direction))
  }
  # Defensive check for empty or NULL results
  if (is.null(enrichment_result) || is.null(enrichment_result$result) || nrow(enrichment_result$result) == 0) {
    return(NULL)
  }

  # Extract the significance threshold from the result object metadata
  significance_threshold <- enrichment_result$meta$query_metadata$user_threshold

  # Prepare data for plotting
  plot_data <- enrichment_result$result %>%
    dplyr::mutate(
      neg_log10_p = -log10(p_value),
      # Ensure 'source' is a factor with a consistent level order for plotting
      source = factor(source, levels = c("GO:BP", "GO:CC", "GO:MF", "KEGG", "REAC"))
    ) %>%
    # Drop any levels that are not actually present in the data to avoid empty spaces on the plot
    dplyr::mutate(source = forcats::fct_drop(source))

  # Identify the top significant terms to add labels for
  top_terms <- plot_data %>%
    dplyr::group_by(source) %>%
    dplyr::arrange(p_value) %>%
    dplyr::slice_head(n = 3) %>%
    dplyr::ungroup()

  # Create the static ggplot object
  ggplot2::ggplot(plot_data, ggplot2::aes(x = source, y = neg_log10_p, text = paste0(
    "Term: ", term_name, "\n",
    "ID: ", term_id, "\n",
    "P-value: ", signif(p_value, 3), "\n",
    "Genes: ", intersection_size
  ))) +
    # Add a line for the significance threshold
    ggplot2::geom_hline(yintercept = -log10(significance_threshold), linetype = "dashed", color = "red") +
    ggplot2::geom_jitter(ggplot2::aes(color = source, size = term_size), alpha = 0.7,
      position = ggplot2::position_jitter(width = 0.2, seed = ENRICHMENT_PLOT_SEED)) +
    ggrepel::geom_text_repel(
      data = top_terms,
      ggplot2::aes(label = term_name),
      size = 3, max.overlaps = 15, box.padding = 0.5, point.padding = 0.3, force = 5,
      seed = ENRICHMENT_PLOT_SEED
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
}

#' Carry what a clusterProfiler result's plot needs beside its table (KIT-119)
#'
#' The plot needs each term's GO category and display name, which come from the run's GO annotation maps, and the
#' run's q-value cutoff. They are carried as the attribute `enrichment_plot_inputs` of the stored `enrichResult`, so
#' the plot can be rebuilt from the stored result alone. The `@result` table is not changed: the tab's results table
#' and its downloads show that table whole.
#' @param enrichment_result The `enrichResult`.
#' @param plot_data The run's plot table: `@result` joined to the GO maps, with `ID`, `source` and `term`.
#' @param q_cutoff The run's q-value cutoff.
#' @param significance_rule KIT-203: `"bh_whole_family"` on a run flagged by the whole-family value, carried so the
#'   plot is drawn from `p_adjust_family`; NULL (a `"qvalue"` run) adds nothing, as before.
#' @return `enrichment_result` with the attribute: `list(go_terms, q_cutoff)`, `go_terms` the distinct `ID`,
#'   `source` and `term` of `plot_data`, in its order, and `significance_rule` when given.
#' @noRd
setEnrichmentPlotInputs <- function(enrichment_result, plot_data, q_cutoff, significance_rule = NULL) {
  go_terms <- data.frame(ID = plot_data$ID, source = as.character(plot_data$source), term = plot_data$term,
    stringsAsFactors = FALSE)
  inputs <- list(go_terms = unique(go_terms), q_cutoff = q_cutoff)
  if (!is.null(significance_rule)) inputs$significance_rule <- significance_rule
  attr(enrichment_result, "enrichment_plot_inputs") <- inputs
  enrichment_result
}

#' The static plot of a clusterProfiler result (KIT-119: `processEnrichments`' plot, moved here)
#' @param enrichment_result The stored `enrichResult`, with `enrichment_plot_inputs` when the run set them. Without
#'   them (a result stored before KIT-119) every term is `Other` and named by its `Description`, and the threshold is
#'   the result's `pvalueCutoff`, which the run sets to its q-value cutoff. KIT-203: inputs naming
#'   `significance_rule = "bh_whole_family"` draw the whole-family value (`p_adjust_family`) in place of the q-value;
#'   any other stored result (every result stored before KIT-203, and every `"qvalue"` run) is drawn as before.
#' @param contrast,direction As `buildEnrichmentGgplot`.
#' @return A ggplot, or NULL when the result has no terms.
#' @noRd
buildClusterProfilerEnrichmentGgplot <- function(enrichment_result, contrast, direction) {
  if (nrow(enrichment_result@result) == 0) {
    return(NULL)
  }
  inputs <- attr(enrichment_result, "enrichment_plot_inputs", exact = TRUE)
  q_cutoff <- inputs$q_cutoff %||% enrichment_result@pvalueCutoff
  go_terms <- inputs$go_terms %||%
    data.frame(ID = character(), source = character(), term = character(), stringsAsFactors = FALSE)
  # KIT-203: the value drawn and its names, by the run's rule (the q-value's are the earlier versions').
  family <- identical(inputs$significance_rule, "bh_whole_family")
  value <- if (family) "p_adjust_family" else "qvalue"
  value_tooltip <- if (family) "Adjusted p (all terms): " else "Q-value: "
  value_legend <- if (family) "-log10(adjusted p-value)" else "-log10(q-value)"
  value_axis <- if (family) "-log10(adjusted p-value, all terms of 10 to 500 proteins)" else "-log10(q-value)"

  plot_data <- enrichment_result@result |>
    dplyr::left_join(go_terms, by = "ID") |>
    dplyr::mutate(
      source = dplyr::coalesce(.data$source, "Other"),
      source = factor(.data$source, levels = c("GO:BP", "GO:CC", "GO:MF", "Other")),
      neg_log10_q = -log10(.data[[value]]),  # the q-value from clusterProfiler, or the whole-family value
      gene_count = .data$Count,
      term = dplyr::coalesce(.data$term, .data$Description)
    )

  # Generate static plot with the cutoff's threshold line
  # First identify top significant terms to label
  top_terms <- plot_data %>%
    dplyr::filter(.data[[value]] < q_cutoff) %>%
    dplyr::arrange(.data[[value]]) %>%
    dplyr::slice_head(n = 5) # Label top 5 most significant terms

  ggplot2::ggplot(plot_data,
                  ggplot2::aes(x = .data$source,
                               y = .data$neg_log10_q,
                               text = paste0(
                                 "Term: ", .data$term, "\n",
                                 "ID: ", .data$ID, "\n",
                                 "Genes: ", .data$Count, "\n",
                                 "Gene Ratio: ", .data$GeneRatio, "\n",
                                 "Background Ratio: ", .data$BgRatio, "\n",
                                 value_tooltip, signif(.data[[value]], 3)
                               ))) +
    ggplot2::geom_hline(yintercept = -log10(q_cutoff),
                        linetype = "dashed",
                        color = "darkgrey") +
    ggplot2::geom_jitter(ggplot2::aes(size = .data$gene_count,
                                      color = -log10(.data[[value]])),
                         alpha = 0.7,
                         position = ggplot2::position_jitter(width = 0.2, seed = ENRICHMENT_PLOT_SEED)) +
    # Add labels for significant terms
    ggrepel::geom_text_repel(
      data = top_terms,
      ggplot2::aes(label = .data$term),
      size = 3,
      max.overlaps = 15,
      box.padding = 0.5,
      point.padding = 0.3,
      force = 5,
      seed = ENRICHMENT_PLOT_SEED
    ) +
    ggplot2::scale_color_gradient(low = "#FED976",
                                  high = "#800026",
                                  name = value_legend) +
    ggplot2::scale_size_continuous(name = "Gene Count",
                                   range = c(3, 12)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 10, angle = 0),
      axis.text.y = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.title = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 8)
    ) +
    ggplot2::labs(
      title = paste0(contrast, " ", tools::toTitleCase(direction), "-regulated"),
      x = "GO Category",
      y = value_axis
    )
}

#' The static and interactive enrichment plot of one contrast and direction, from its stored result (KIT-119)
#' @param enrichment_result,contrast,direction As `buildEnrichmentGgplot`.
#' @return `list(static, interactive)`: the ggplot (`buildEnrichmentGgplot`) and its plotly conversion; both NULL
#'   when the result is NULL or has no terms.
#' @noRd
buildEnrichmentPlot <- function(enrichment_result, contrast, direction) {
  static <- buildEnrichmentGgplot(enrichment_result, contrast, direction)
  if (is.null(static)) {
    return(list(static = NULL, interactive = NULL))
  }
  list(static = static, interactive = plotly::ggplotly(static, tooltip = "text"))
}

#' The interactive plot the enrichment tab shows for one contrast and direction (KIT-119)
#'
#' A plot the results already hold (a session saved before KIT-119, or the deterministic test kernel's) is shown as
#' stored; otherwise the plot is built from the stored result (`buildEnrichmentPlot`), as the run built its saved file.
#' @param enrichment_results `EnrichmentResults`.
#' @param contrast Contrast name, as the results' lists are named.
#' @param direction `"up"` or `"down"`.
#' @return A plotly object, or NULL when the contrast and direction have no result with terms.
#' @noRd
resolveEnrichmentPlotly <- function(enrichment_results, contrast, direction) {
  stored <- enrichment_results@enrichment_plotly[[contrast]][[direction]]
  if (!is.null(stored)) {
    return(stored)
  }
  buildEnrichmentPlot(enrichment_results@enrichment_data[[contrast]][[direction]], contrast, direction)$interactive
}


# ----------------------------------------------------------------------------
# generate_enrichment_plots
# ----------------------------------------------------------------------------
# Plot generation function
#' @export
generate_enrichment_plots <- function(enrichment_result, contrast, direction, pathway_dir) {
  # KIT-119: the plot is drawn by the builder the tab uses at render time.
  static_plot <- buildEnrichmentGgplot(enrichment_result, contrast, direction)
  if (is.null(static_plot)) {
    return(list(static = NULL, interactive = NULL))
  }

  # Convert to an interactive plotly object
  interactive_plot <- tryCatch({
      plotly::ggplotly(static_plot, tooltip = "text")
  }, error = function(e) {
      warning(paste("Failed to convert ggplot to plotly for", contrast, direction, ":", e$message))
      return(NULL)
  })

  # Save the results data table
  tryCatch({
    result_table <- enrichment_result$result
    result_table$parents <- sapply(result_table$parents, paste, collapse = ", ")
    readr::write_tsv(result_table,
      file = file.path(pathway_dir, paste0(contrast, "_", direction, "_enrichment_results.tsv"))
    )
  }, error = function(e) {
    warning(paste("Failed to write enrichment results table for", contrast, direction, ":", e$message))
  })
  
  return(list(
    static = static_plot,
    interactive = interactive_plot
  ))
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

