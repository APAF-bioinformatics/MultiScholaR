# KIT-093: the notebook's expression filter (dev/distributable/plans/transcriptomics_workflow.rmd:166-175) as a kernel.
# The edgeR calls are the notebook's, unchanged.

#' The notebook's expression filter.
#' @param counts Numeric genes-by-runs count matrix in design order (whole numbers; checked).
#' @param group The design's group factor, in the same order.
#' @return A list: `effective_lib_size` (lib.size x TMM norm.factors of the unfiltered counts, the sizes
#'   `filterByExpr.DGEList` uses), `min_sample_size` and `median_lib_size` (as `filterByExpr.default` computes them,
#'   recorded for a partitioned variant; the mask always comes from the edgeR call), `keep_expr` (the
#'   `edgeR::filterByExpr` mask), `keep_var` (the zero-variance mask on the expression-filtered counts) and
#'   `filtered_counts`. Refuses `trx_filter_empty` when no gene survives.
#' @noRd
trxExpressionFilter <- function(counts, group) {
    validateTrxCounts(counts)
    requireOptionalPackage("edgeR", "the transcriptomics expression filter")
    dge_raw <- edgeR::DGEList(counts = counts, group = group)
    dge_raw <- edgeR::calcNormFactors(dge_raw, method = "TMM")
    keep_expr <- edgeR::filterByExpr(dge_raw, group = group)
    # drop = FALSE: a single surviving gene stays a matrix (the notebook would fail there); identical otherwise.
    filtered <- counts[keep_expr, , drop = FALSE]
    gene_vars <- apply(filtered, 1, stats::var)
    keep_var <- gene_vars > 0
    filtered <- filtered[keep_var, , drop = FALSE]
    if (!length(keep_var) || !any(keep_var)) {
        trxAbort("no gene survives the expression and zero-variance filters", c("trx_filter_empty", "trx_error"))
    }
    effective <- dge_raw$samples$lib.size * dge_raw$samples$norm.factors
    list(effective_lib_size = effective, min_sample_size = trxFilterMinSampleSize(group),
        median_lib_size = stats::median(effective), keep_expr = keep_expr, keep_var = keep_var, filtered_counts = filtered)
}

#' filterByExpr.default's minimum group size, by edgeR's own rule and defaults (large.n, min.prop).
#' @noRd
trxFilterMinSampleSize <- function(group) {
    f <- formals(utils::getS3method("filterByExpr", "default", envir = asNamespace("edgeR")))
    n <- tabulate(as.factor(group))
    m <- min(n[n > 0L])
    if (m > f$large.n) m <- f$large.n + (m - f$large.n) * f$min.prop
    m
}
