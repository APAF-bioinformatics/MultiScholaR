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

#' The shipped significance rule of the clusterProfiler route (KIT-203; VAL-123, VAL-128)
#' @noRd
enrichmentSignificanceShipped <- "bh_whole_family"

#' The significance rules of the clusterProfiler route and what the user is shown of each (KIT-203 5.4)
#'
#' One row per rule: `label` (the tab's selector), `note` (shown beside the selector), `notice` (the run's notice),
#' `methods` (the report paragraph's rule-specific sentences; `{cutoff}` is the run's cutoff) and `source` (where the
#' figures come from; options-k203.R regenerates every figure). Proposed wording (VAL-41): the owner edits before release.
#' @noRd
ENRICHMENT_SIGNIFICANCE_RULES <- data.frame(
  rule = c("bh_whole_family", "qvalue"),
  label = c("Benjamini-Hochberg over every term of 10 to 500 tested proteins (default)",
    "Storey q-value over the terms with at least one listed protein (as in earlier versions)"),
  note = c(paste("Each term's p-value is adjusted together with every Gene Ontology term (biological process, molecular",
    "function and cellular component together) that has 10 to 500 of your tested proteins, including terms with none",
    "of the listed proteins (they count as p = 1). In a check on simulated protein lists at a cutoff of 0.05 (reviewed",
    "E. coli and yeast annotations; lists of 25, 126 and 400 proteins drawn at random; with and without parent terms;",
    "2,000 random lists per setting), lists with no real signal had at least one significant term in 1.0 % to 2.5 % of",
    "cases (nominal 5 %). On lists of 126 proteins built around one process, the average share of a list's significant",
    "terms that were false (0 when nothing was called) was 0.1 % to 1.7 % (counted over all significant terms of those",
    "lists together, at most 3.0 %), and the process's own term was found in 74.8 % to 100.0 % of the lists in which it",
    "was tested: up to 5.8 percentage points less often than with the q-value with parent terms and up to 8.8 points",
    "less often without them, never more often. Not simulated: down lists, other organisms, and lists of proteins that",
    "change together as in real data. Its error guarantee assumes that the terms' tests are independent or tend to",
    "agree; overlapping Gene Ontology terms do not ensure this, so the check measured it instead, in the settings named",
    "only."),
    paste("Terms with none of your listed proteins are left out before the adjustment. This breaks the q-value's",
    "assumption (p-values of terms with no real change spread evenly between 0 and 1), so the q-value can estimate that",
    "most terms are real. In the same check, lists with no real signal had at least one significant term in 1.7 % to",
    "42.9 % of cases: clearly above 7.5 %, the tolerance allowed for a nominal 5 % (lower confidence limit above it), in",
    "6 of 18 settings and above it without certainty in 2 more; with parent terms the clear cases were lists of 25",
    "proteins (E. coli full proteome 22.9 %, its 2,000-protein subset 16.7 %). Up to 23.6 % of such lists had 10 or more",
    "significant terms. On lists of 126 proteins built around one process, the average share of a list's significant",
    "terms that were false (0 when nothing was called) was up to 35.7 %, clearly above 7.5 % in 3 of 18 settings, all",
    "without parent terms; counted over all significant terms of those lists together it was up to 82.1 %. On lists of",
    "25 with parent terms the average was 10.7 %. It found the process's own term up to 5.8 percentage points more often",
    "than the default with parent terms and up to 8.8 points more often without.")),
  notice = c(paste("For this organism, terms are called significant by the Benjamini-Hochberg adjusted p-value over every",
    "term of 10 to 500 tested proteins."),
    paste("For this organism, terms are called significant by the Storey q-value over the terms with at least one listed",
    "protein (Benjamini-Hochberg adjusted p-values over the same terms are also given).")),
  methods = c(paste("A term was called significant when its Benjamini-Hochberg adjusted p-value over that family was below",
    "{cutoff}; Storey's q-value and the Benjamini-Hochberg value over the terms with at least one listed protein are",
    "also given in the results table. In simulations through this route at a cutoff of 0.05 (lists drawn at random from",
    "the reviewed UniProt proteomes of E. coli K-12, a 2,000-protein subset of it, and S. cerevisiae, the last run",
    "through this route as a stand-in; lists of 25, 126 and 400 proteins; with parent terms over the five relations and",
    "without them; 2,000 random lists per setting), 1.0 % to 2.5 % of random lists had any significant term, and on",
    "lists of 126 built around one process the mean false discovery proportion was at most 1.7 % (in a design where the",
    "other listed proteins avoid that process, which tends to lower false calls). Not simulated: down lists, other",
    "cutoffs, lists of proteins that change together as in real data, other organisms or annotation sources, and parent",
    "terms over is a and part of only. Transcriptomics tests sets of 5 to 500 genes; this route 10 to 500 proteins."),
    paste("Because terms without a list protein are left out before the adjustment, the estimated share of true null",
    "terms falls far below 1 on many random lists; on simulated random lists at a cutoff of 0.05 the share of lists",
    "with any significant term was 1.7 % to 42.9 % depending on the proteome, list size and parent terms (failing a",
    "7.5 % line in 6 of 18 settings, two of them with parent terms, both with lists of 25 proteins), against 1.6 % to",
    "8.2 % for Benjamini-Hochberg over the same terms.")),
  source = c(paste("KIT-189 simulation, score run 45311934798, results/options.tsv and classifications.tsv; the paired",
    "comparison with the q-value and the share over all significant terms computed after scoring (results.md section 5)"),
    paste("KIT-189 simulation, score run 45311934798, results/options.tsv and classifications.tsv; the paired comparison",
    "with the default, the share over all significant terms and the lists of 25 computed after scoring (results.md",
    "section 5)")),
  stringsAsFactors = FALSE
)

#' Check a significance rule of the clusterProfiler route (KIT-203)
#' @param significance_rule One of `ENRICHMENT_SIGNIFICANCE_RULES$rule`.
#' @return `significance_rule`, invisibly; otherwise an error of class `enrichment_significance_rule_invalid`.
#' @noRd
checkEnrichmentSignificanceRule <- function(significance_rule) {
  if (!is.character(significance_rule) || length(significance_rule) != 1L || is.na(significance_rule) ||
      !significance_rule %in% ENRICHMENT_SIGNIFICANCE_RULES$rule) {
    cli::cli_abort(c("{.arg significance_rule} must be one of {.val {ENRICHMENT_SIGNIFICANCE_RULES$rule}}",
      x = "got {.val {format(significance_rule)}}"), class = "enrichment_significance_rule_invalid")
  }
  invisible(significance_rule)
}

#' Benjamini-Hochberg over the whole size-filtered family of one enricher call (KIT-203 5.1)
#'
#' enrichit keeps in `@geneSets` every term of 10 to 500 members in the universe (m terms) and in `@result` only the m'
#' terms with at least one list protein. The family value adjusts each kept term's p-value together with the m - m'
#' terms of count 0, which enter with p = 1 (the exact p-value of a count of 0):
#' `p.adjust(c(p, rep(1, m - m')), "BH")[seq_len(m')]`, in `@result`'s own row order. `@geneSets` is read directly: no
#' accessor of the class returns the size-filtered family (a recorded divergence from the accessor rule).
#' @param result An `enrichResult`, or NULL (no list in that direction, or `enricher` returned NULL).
#' @param significance_rule The run's rule; only `"bh_whole_family"` adds the columns.
#' @return `result` with `p_adjust_family` and `family_size` (m) appended to `@result`; NULL, a result with no rows and
#'   any result under `"qvalue"` unchanged. A result without `pvalue`, with an NA p-value, without `geneSets`, or with
#'   m < m' is refused (class `enrichment_family_unavailable`).
#' @noRd
addEnrichmentFamilyAdjustment <- function(result, significance_rule) {
  if (is.null(result) || !identical(significance_rule, "bh_whole_family")) return(result)
  refuse <- function(why) cli::cli_abort(paste("the whole-family adjustment cannot be formed:", why),
    class = "enrichment_family_unavailable")
  if (!isS4(result) || !"result" %in% methods::slotNames(result) || !is.data.frame(result@result)) {
    refuse("the result has no result table")
  }
  if (!nrow(result@result)) return(result)
  if (!"geneSets" %in% methods::slotNames(result)) refuse("the result has no gene sets")
  p <- result@result$pvalue
  if (!is.numeric(p) || anyNA(p)) refuse("the result has no p-value for every term")
  m <- length(result@geneSets)
  if (m < length(p)) refuse(sprintf("%d gene sets for %d tested terms", m, length(p)))
  result@result$p_adjust_family <- stats::p.adjust(c(p, rep(1, m - length(p))), method = "BH")[seq_along(p)]
  result@result$family_size <- m
  result
}

#' The rows of one comparison's DA table that were tested (KIT-205; VAL-129 (2))
#'
#' A protein without a p-value in a comparison (no value in a compared group) can never be in that comparison's lists,
#' so it is not in their background (universe) either. Tested means a finite `fdr_qvalue` and, where the table has the
#' column, a finite `raw_pvalue`; on the package's DA tables the two coincide. A table with no untested row is returned
#' as it is. Under the earlier rule (the option of `workflowReplayRuleMarkers$enrichment_background` holding its
#' earlier token, set only by the record replay of a step made before this rule) every row is kept.
#' @param de_data One comparison's DA table.
#' @return `list(rows, left_out)`: the rows the lists and the background are formed from, and how many were left out.
#' @noRd
enrichmentTestedRows <- function(de_data) {
  marker <- workflowReplayRuleMarkers$enrichment_background
  if (identical(getOption(marker$option), marker$earlier)) return(list(rows = de_data, left_out = 0L))
  tested <- is.finite(de_data[["fdr_qvalue"]])
  if ("raw_pvalue" %in% names(de_data)) tested <- tested & is.finite(de_data[["raw_pvalue"]])
  if (all(tested)) return(list(rows = de_data, left_out = 0L))
  list(rows = de_data[tested, , drop = FALSE], left_out = sum(!tested))
}

# processEnrichments
# ----------------------------------------------------------------------------
#' Process Enrichments
#'
#' @param da_results S4 object containing differential expression results
#' @param taxon_id NCBI taxonomy ID for the organism
#' @param up_cutoff Log2 fold change cutoff for up-regulated proteins (default: 0)
#' @param down_cutoff Log2 fold change cutoff for down-regulated proteins (default: 0)
#' @param q_cutoff Cutoff for the differential-abundance q-value that forms the lists and for the term's significance
#'   value (default: 0.05)
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
#' @param significance_rule For organisms outside the 13 served by g:Profiler, how a term is called significant
#'   (KIT-203): `"bh_whole_family"` (the default) by its Benjamini-Hochberg adjusted p-value over every term of 10 to
#'   500 tested proteins, terms without a listed protein entering with p = 1 (columns `p_adjust_family` and
#'   `family_size` added to the result); `"qvalue"` by Storey's q-value over the terms with at least one listed protein,
#'   as in earlier versions (in simulations it called terms on up to 42.9 % of random protein lists). Not used on the
#'   g:Profiler route.
#'
#' @details The background (universe) of each comparison is the proteins tested in that comparison (KIT-205, VAL-129
#'   (2)): before the up and down lists and the background are formed, the comparison's table is reduced to its rows
#'   with a finite `fdr_qvalue` and, where the table has a `raw_pvalue` column, a finite `raw_pvalue` (a table without
#'   `raw_pvalue` is reduced by `fdr_qvalue` alone). A protein without a p-value could never be in a list, so it is not
#'   counted in the background either. One background per comparison serves both lists, on both routes; the number of
#'   rows left out is written to the log. Analyses recorded before this rule replay with every row.
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
                               go_propagation = goPropagationShipped,
                               significance_rule = enrichmentSignificanceShipped) {

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

        # KIT-205 (VAL-129 (2)): the lists and the background are formed from the rows tested in this comparison.
        reduced <- enrichmentTestedRows(de_data)
        de_data <- reduced$rows
        
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

        message(sprintf("Using %d unique proteins as background for enrichment analysis (%d rows without a p-value left out)",
          length(custom_bg), reduced$left_out))
        
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
    checkEnrichmentSignificanceRule(significance_rule) # KIT-203
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

        # KIT-205 (VAL-129 (2)): the lists and the background are formed from the rows tested in this comparison.
        reduced <- enrichmentTestedRows(de_data)
        de_data <- reduced$rows

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
        message(sprintf("Using %d unique proteins as background for enrichment analysis (%d rows without a p-value left out)",
          length(background_IDs), reduced$left_out))

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

        # KIT-203: the whole-family value, formed after each call and outside the calls' tryCatch, so a result it
        # cannot be formed for is refused, not turned into "no result".
        up_enrich <- addEnrichmentFamilyAdjustment(up_enrich, significance_rule)
        down_enrich <- addEnrichmentFamilyAdjustment(down_enrich, significance_rule)

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
    family <- identical(significance_rule, "bh_whole_family") # KIT-203
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
                    # KIT-203: the flag by the run's rule (the q-value's is the earlier versions' flag)
                    significant = if (family) .data$p_adjust_family < q_cutoff else .data$qvalue < q_cutoff
                  ) |>
                  dplyr::mutate(
                    term = dplyr::coalesce(.data$term, .data$Description)
                  )
                if (family) {
                  plot_data$significance_rule <- significance_rule
                  plot_data <- dplyr::relocate(plot_data, "significance_rule", .after = "significant")
                }

                # Save results table
                readr::write_tsv(
                  plot_data,
                  file.path(pathway_dir,
                            paste0(contrast, "_", direction, "_enrichment_results.tsv"))
                )

                # KIT-119: the stored result carries its terms' GO category and name and the cutoff, and the plot
                # is built from that stored result by the builder the tab uses at render time.
                result_data <- setEnrichmentPlotInputs(result_data, plot_data, q_cutoff,
                  significance_rule = if (family) significance_rule)
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
