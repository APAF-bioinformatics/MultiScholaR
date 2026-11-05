# KIT-093: the transcriptomics whole-number count policy. Counts are stored as double (as vroom reads them); only
# integer-valued, finite, non-negative values pass, and every sample must have a non-zero library.

#' Raise a classed transcriptomics refusal.
#' @param message The message.
#' @param class The specific class, followed by its parent classes.
#' @return Does not return.
#' @noRd
trxAbort <- function(message, class) {
    rlang::abort(message, class = class)
}

#' Validate transcriptomics counts against the whole-number policy.
#' Runs before any statistics: in the class validity and at the entry of every kernel that takes a count matrix.
#' Callable on its own so that an import can show the refusal by class.
#' @param counts A numeric genes-by-samples matrix (or a data frame of numeric sample columns).
#' @return `counts`, unchanged, when it passes. Otherwise a classed error inheriting `trx_counts_invalid`:
#'   `trx_counts_na` (NA or NaN), `trx_counts_nonfinite` (`Inf`, `-Inf`), `trx_counts_negative`,
#'   `trx_counts_fractional` (not whole) or `trx_counts_zero_library` (a sample summing to zero).
#' @noRd
validateTrxCounts <- function(counts) {
    x <- as.matrix(counts)
    refuse <- function(fault, what, cells) {
        at <- arrayInd(which(cells)[[1L]], dim(x))
        trxAbort(sprintf("%s count (%d cell(s); first: gene %s, sample %s)", what, sum(cells),
            rownames(x)[at[[1L]]] %||% at[[1L]], colnames(x)[at[[2L]]] %||% at[[2L]]), c(fault, "trx_counts_invalid", "trx_error"))
    }
    if (!is.numeric(x)) trxAbort("counts must be numeric", c("trx_counts_invalid", "trx_error"))
    if (anyNA(x)) refuse("trx_counts_na", "NA or NaN", is.na(x))
    if (any(is.infinite(x))) refuse("trx_counts_nonfinite", "non-finite", is.infinite(x))
    if (any(x < 0)) refuse("trx_counts_negative", "negative", x < 0)
    if (any(x != round(x))) refuse("trx_counts_fractional", "non-whole", x != round(x))
    empty <- colSums(x) == 0
    if (any(empty)) {
        trxAbort(sprintf("sample(s) with a zero library: %s", paste(colnames(x)[empty] %||% which(empty), collapse = ", ")),
            c("trx_counts_zero_library", "trx_counts_invalid", "trx_error"))
    }
    counts
}
