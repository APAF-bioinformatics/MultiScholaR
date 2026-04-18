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
        # Validate that all objects in the list are MetaboliteAssayData
        objectsList <- theObject;
            if (!all(purrr::map_lgl(objectsList, ~ inherits(.x, "MetaboliteAssayData")))) {
            stop("All objects in objectsList must be of class MetaboliteAssayData")
        }

        # Run DE analysis and explicitly set names
        results_list <- purrr::map(
            objectsList,
            \(obj) {
                differentialAbundanceAnalysisHelper(obj,
                    contrasts_tbl = contrasts_tbl,
                    formula_string = formula_string,
                    group_id = group_id,
                    da_q_val_thresh = da_q_val_thresh,
                    treat_lfc_cutoff = treat_lfc_cutoff,
                    eBayes_trend = eBayes_trend,
                    eBayes_robust = eBayes_robust,
                    args_group_pattern = args_group_pattern
                )
            }
        )

        # Set names if the input list had names
        if (!is.null(names(objectsList))) {
            names(results_list) <- names(objectsList)
        }

        return(results_list)
    }
)

