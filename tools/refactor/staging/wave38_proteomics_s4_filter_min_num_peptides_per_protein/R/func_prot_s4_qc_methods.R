#' @export
setMethod(
  f = "filterMinNumPeptidesPerProtein",
  signature = "ProteinQuantitativeData",
  definition = function(theObject, ...) {
    # Extract specific parameters from ...
    args <- list(...)
    num_peptides_per_protein_thresh <- args$num_peptides_per_protein_thresh
    num_peptidoforms_per_protein_thresh <- args$num_peptidoforms_per_protein_thresh
    verbose <- args$verbose

    # --- Parameter validation and defaults ---
    num_peptides_per_protein_thresh <- checkParamsObjectFunctionSimplify(
      theObject,
      "num_peptides_per_protein_thresh",
      1
    )

    num_peptidoforms_per_protein_thresh <- checkParamsObjectFunctionSimplify(
      theObject,
      "num_peptidoforms_per_protein_thresh",
      2
    )

    verbose <- checkParamsObjectFunctionSimplify(theObject, "verbose", TRUE)

    # Update parameters in object
    theObject <- updateParamInObject(theObject, "num_peptides_per_protein_thresh")
    theObject <- updateParamInObject(theObject, "num_peptidoforms_per_protein_thresh")
    theObject <- updateParamInObject(theObject, "verbose")

    if (verbose) {
      log_info("Starting protein filtering based on peptide and peptidoform evidence...")
      log_info("Minimum unique peptides per protein: {num_peptides_per_protein_thresh}")
      log_info("Minimum total peptidoforms per protein: {num_peptidoforms_per_protein_thresh}")
    }

    # Get the peptide summary table (which is now guaranteed to be in sync)
    peptide_summary <- theObject@args$limpa_dpc_quant_results$peptide_counts_per_protein
    if (is.null(peptide_summary)) {
      stop("Could not find the peptide summary table. Please run chooseBestProteinAccession first.")
    }

    # --- Perform the filtering ---
    protein_quant_table <- theObject@protein_quant_table
    protein_id_column <- theObject@protein_id_column
    proteins_before <- nrow(protein_quant_table)

    protein_ids_to_keep <- peptide_summary |>
      dplyr::filter(peptide_count >= num_peptides_per_protein_thresh & peptidoform_count >= num_peptidoforms_per_protein_thresh) |>
      dplyr::pull(!!sym(protein_id_column))

    filtered_protein_table <- protein_quant_table |>
      dplyr::filter(!!sym(protein_id_column) %in% protein_ids_to_keep)

    proteins_after <- nrow(filtered_protein_table)

    if (verbose) {
      log_info("Proteins before filtering: {proteins_before}")
      log_info("Proteins after filtering: {proteins_after}")
      log_info("Proteins removed: {proteins_before - proteins_after}")
      if (proteins_before > 0) {
        log_info("Retention rate: {round(100 * proteins_after / proteins_before, 1)}%")
      }
    }

    # Update the main data table
    theObject@protein_quant_table <- filtered_protein_table

    # Also filter the EList for consistency
    if (!is.null(theObject@args$limpa_dpc_quant_results$quantified_elist)) {
      original_elist <- theObject@args$limpa_dpc_quant_results$quantified_elist
      indices_to_keep <- which(original_elist$genes$protein.id %in% protein_ids_to_keep)
      theObject@args$limpa_dpc_quant_results$quantified_elist <- original_elist[indices_to_keep, ]
    }

    return(theObject)
  }
)

