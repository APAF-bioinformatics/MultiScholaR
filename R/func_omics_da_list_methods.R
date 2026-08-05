# ----------------------------------------------------------------------------
# Shared differential-abundance list dispatch
# ----------------------------------------------------------------------------
.resolveOmicsDaListDomain <- function(objects_list, metabolomics_class, lipidomics_class) {
  if (!is.list(objects_list)) {
    stop("Differential-abundance list methods require a list input.", call. = FALSE)
  }
  if (!length(objects_list)) {
    return("lipidomics")
  }

  is_metabolomics <- vapply(
    objects_list,
    inherits,
    logical(1L),
    what = metabolomics_class
  )
  is_lipidomics <- vapply(
    objects_list,
    inherits,
    logical(1L),
    what = lipidomics_class
  )

  if (all(is_metabolomics)) {
    return("metabolomics")
  }
  if (all(is_lipidomics)) {
    return("lipidomics")
  }

  stop(
    sprintf(
      "List elements must all inherit from either %s or %s.",
      metabolomics_class,
      lipidomics_class
    ),
    call. = FALSE
  )
}

#' @export
setMethod(
  f = "differentialAbundanceAnalysis",
  signature = "list",
  definition = function(
    theObject,
    contrasts_tbl = NULL,
    formula_string = NULL,
    group_id = NULL,
    da_q_val_thresh = NULL,
    treat_lfc_cutoff = NULL,
    eBayes_trend = NULL,
    eBayes_robust = NULL,
    args_group_pattern = NULL,
    args_row_id = NULL,
    qvalue_column = NULL,
    raw_pvalue_column = NULL
  ) {
    domain <- .resolveOmicsDaListDomain(
      theObject,
      "MetaboliteAssayData",
      "LipidomicsAssayData"
    )
    helper <- if (identical(domain, "metabolomics")) {
      .differentialAbundanceAnalysisMetabolomicsList
    } else {
      .differentialAbundanceAnalysisLipidomicsList
    }

    helper(
      theObject,
      contrasts_tbl = contrasts_tbl,
      formula_string = formula_string,
      group_id = group_id,
      da_q_val_thresh = da_q_val_thresh,
      treat_lfc_cutoff = treat_lfc_cutoff,
      eBayes_trend = eBayes_trend,
      eBayes_robust = eBayes_robust,
      args_group_pattern = args_group_pattern,
      args_row_id = args_row_id,
      qvalue_column = qvalue_column,
      raw_pvalue_column = raw_pvalue_column
    )
  }
)

#' @export
setMethod(
  f = "plotNumSigDiffExpBarPlot",
  signature = "list",
  definition = function(objectsList) {
    domain <- .resolveOmicsDaListDomain(
      objectsList,
      "MetabolomicsDifferentialAbundanceResults",
      "LipidomicsDifferentialAbundanceResults"
    )
    helper <- if (identical(domain, "metabolomics")) {
      .plotNumSigDiffExpBarPlotMetabolomicsList
    } else {
      .plotNumSigDiffExpBarPlotLipidomicsList
    }
    helper(objectsList)
  }
)

#' @export
setMethod(
  f = "plotVolcanoS4",
  signature = "list",
  definition = function(
    objectsList,
    da_q_val_thresh = 0.05,
    qvalue_column = "fdr_qvalue",
    log2fc_column = "logFC"
  ) {
    domain <- .resolveOmicsDaListDomain(
      objectsList,
      "MetabolomicsDifferentialAbundanceResults",
      "LipidomicsDifferentialAbundanceResults"
    )
    helper <- if (identical(domain, "metabolomics")) {
      .plotVolcanoS4MetabolomicsList
    } else {
      .plotVolcanoS4LipidomicsList
    }
    helper(
      objectsList,
      da_q_val_thresh = da_q_val_thresh,
      qvalue_column = qvalue_column,
      log2fc_column = log2fc_column
    )
  }
)

#' @export
setMethod(
  f = "plotInteractiveVolcano",
  signature = "list",
  definition = function(objectsList, anno_list = NULL) {
    domain <- .resolveOmicsDaListDomain(
      objectsList,
      "MetabolomicsDifferentialAbundanceResults",
      "LipidomicsDifferentialAbundanceResults"
    )
    helper <- if (identical(domain, "metabolomics")) {
      .plotInteractiveVolcanoMetabolomicsList
    } else {
      .plotInteractiveVolcanoLipidomicsList
    }
    helper(objectsList, anno_list = anno_list)
  }
)

# ----------------------------------------------------------------------------
# getDaResultsLongFormat
# ----------------------------------------------------------------------------
# Get the differential abundance results in wide format
.resolveDaAccessorSource <- function(object, assay_name = NULL) {
  source_object <- object@theObject
  source_slots <- methods::slotNames(source_object)

  if ("lipid_data" %in% source_slots && "lipid_id_column" %in% source_slots) {
    counts_data_slot <- source_object@lipid_data
    id_col_name <- source_object@lipid_id_column
  } else if ("metabolite_data" %in% source_slots && "metabolite_id_column" %in% source_slots) {
    counts_data_slot <- source_object@metabolite_data
    id_col_name <- source_object@metabolite_id_column
  } else {
    stop(
      sprintf(
        "Unsupported DA result source object class for result accessors: %s",
        paste(class(source_object), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  counts_table_to_use <- if (is.list(counts_data_slot) && !is.data.frame(counts_data_slot)) {
    slot_names <- names(counts_data_slot)
    if (!is.null(assay_name) && nzchar(assay_name) &&
        !is.null(slot_names) && assay_name %in% slot_names) {
      counts_data_slot[[assay_name]]
    } else {
      counts_data_slot[[1]]
    }
  } else {
    counts_data_slot
  }

  list(
    counts_table = counts_table_to_use,
    id_col_name = id_col_name
  )
}

#' @export
setMethod(
  f = "getDaResultsLongFormat",
  signature = "list",
  definition = function(objectsList) {
    return_object_list <- purrr::imap(objectsList, function(object, assay_name) {
      accessor_source <- .resolveDaAccessorSource(object, assay_name)
      counts_table_to_use <- accessor_source$counts_table
      id_col_name <- accessor_source$id_col_name

      # Bind the list of data frames into a single tidy data frame
      tidy_results <- object@contrasts_results_table |>
        dplyr::bind_rows(.id = "comparison") |>
        dplyr::mutate(comparision_short = str_split_i(comparison, "=", 1))

      # Pivot the tidy data frame to a wide format using the correct ID column
      long_results <- tidy_results |>
        # Join with original counts using the correct ID column.
        dplyr::left_join(counts_table_to_use, by = join_by(!!sym(id_col_name) == !!sym(id_col_name))) |>
        dplyr::distinct()

      print(head(long_results))

      # Assign to the correct slot and return the object
      object@results_table_long <- long_results
      return(object)
    })

    return(return_object_list)
  }
)

# ----------------------------------------------------------------------------
# getDaResultsWideFormat
# ----------------------------------------------------------------------------
#' @export
setMethod(
  f = "getDaResultsWideFormat",
  signature = "list",
  definition = function(
    objectsList,
    qvalue_column = "fdr_qvalue",
    raw_pvalue_column = "raw_pvalue",
    log2fc_column = "logFC"
  ) {
    return_object_list <- purrr::imap(objectsList, function(object, assay_name) {
      accessor_source <- .resolveDaAccessorSource(object, assay_name)
      counts_table_to_use <- accessor_source$counts_table
      id_col_name <- accessor_source$id_col_name

      # Bind the list of data frames into a single tidy data frame
      tidy_results <- object@contrasts_results_table |>
        dplyr::bind_rows(.id = "comparison") |>
        dplyr::mutate(comparision_short = str_split_i(comparison, "=", 1))

      # Pivot the tidy data frame to a wide format using the correct ID column
      wida_results <- tidy_results |>
        tidyr::pivot_wider(
          id_cols = c(!!sym(id_col_name)),
          names_from = c(comparision_short),
          names_sep = ":",
          values_from = c(!!sym(log2fc_column), !!sym(qvalue_column), !!sym(raw_pvalue_column))
        ) |>
        # Join with original counts using the correct ID column.
        dplyr::left_join(counts_table_to_use, by = join_by(!!sym(id_col_name) == !!sym(id_col_name))) |>
        dplyr::arrange(dplyr::across(matches(qvalue_column))) |>
        dplyr::distinct()

      # Assign to the correct slot and return the object
      object@results_table_wide <- wida_results
      return(object)
    })

    return(return_object_list)
  }
)
