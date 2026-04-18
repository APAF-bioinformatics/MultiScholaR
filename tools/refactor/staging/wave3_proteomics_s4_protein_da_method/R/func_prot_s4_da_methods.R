#' @export
setMethod(
  f = "differentialAbundanceAnalysis",
  signature = "ProteinQuantitativeData",
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
    qvalue_column = "fdr_qvalue",
    raw_pvalue_column = "raw_pvalue"
  ) {
    # IMMEDIATE ERROR CATCH - Check if we even get here
    message("*** ENTERING differentialAbundanceAnalysis METHOD ***")
    message(sprintf("*** METHOD SIGNATURE MATCHED: ProteinQuantitativeData ***"))

    # Try to catch the index error immediately
    tryCatch(
      {
        message("*** Testing parameter access ***")
        if (!is.null(contrasts_tbl)) {
          test <- contrasts_tbl[[1]]
          message("*** Parameter access successful ***")
        }
      },
      error = function(e) {
        message(sprintf("*** IMMEDIATE ERROR: %s ***", e$message))
        message(sprintf("*** ERROR CLASS: %s ***", class(e)))
        stop(e)
      }
    )

    message("--- Entering differentialAbundanceAnalysis ---")
    message(sprintf("   differentialAbundanceAnalysis: theObject class = %s", class(theObject)))
    message(sprintf("   differentialAbundanceAnalysis: contrasts_tbl provided = %s", !is.null(contrasts_tbl)))
    if (!is.null(contrasts_tbl)) {
      message(sprintf("   differentialAbundanceAnalysis: contrasts_tbl dims = %d x %d", nrow(contrasts_tbl), ncol(contrasts_tbl)))
      message(sprintf("   differentialAbundanceAnalysis: contrasts_tbl content = %s", paste(contrasts_tbl[[1]], collapse = ", ")))
    }

    # Wrap the helper function call in tryCatch to get better error info
    message("   differentialAbundanceAnalysis: About to call differentialAbundanceAnalysisHelper...")

    results_list <- tryCatch(
      {
        differentialAbundanceAnalysisHelper(theObject,
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
      },
      error = function(e) {
        # CRITICAL FIX: Use paste() for logger calls in error handlers to avoid interpolation bug
        message(paste("   differentialAbundanceAnalysis ERROR in helper function:", e$message))
        message(paste("   differentialAbundanceAnalysis ERROR call stack:", capture.output(traceback())))
        stop(e)
      }
    )

    message("   differentialAbundanceAnalysis: Helper function completed successfully!")
    message("--- Exiting differentialAbundanceAnalysis ---")
    return(results_list)
  }
)

