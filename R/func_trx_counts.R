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

# KIT-102 W2: the transcriptomics import's preview and count check, moved into the core from
# mod_transcriptomics.R (`trxImportPreview`, `checkTrxImportCounts`, :82-111). The import itself is
# `importOmicsTables()` (utils_workflow_omics_import.R) with `checkTrxImportCounts` as its `check`.

#' The import's preview of a count table, in `previewOmicsImport()`'s shape: read as the notebook reads
#' it (vroom), the id and name columns proposed (`GeneID`, `GeneName` when present), every other column
#' a sample (`trxImportPreview`).
#' @param path The count table.
#' @param import The omic's `config$import`; not used (the count table has one format and reader).
#' @noRd
previewTrxCountsImport <- function(path, import = NULL) {
    data <- vroom::vroom(path, show_col_types = FALSE, progress = FALSE)
    # The reader's lazy problems pointer is not a value the session can hold by reference.
    attr(data, "problems") <- NULL
    headers <- names(data)
    id <- if ("GeneID" %in% headers) "GeneID" else headers[[1L]]
    name <- if ("GeneName" %in% headers) "GeneName" else setdiff(headers, id)[[1L]]
    list(headers = headers, detected_format = "trx_counts", format_confidence = 1,
        import_result = list(data = data, sample_columns = setdiff(headers, c(id, name))), assay_data = data,
        updates = list(id = list(choices = headers, selected = id),
            annotation = list(choices = headers, selected = name), is_pattern = NULL))
}

#' KIT-093's count policy at import on every column but the id and name columns: the `check` plugin of
#' the transcriptomics import. A refusal keeps KIT-093's class and names it in its message, which the
#' tab shows.
#' @param assayList The assembled assays (one count table).
#' @return `assayList`, invisibly and unchanged, when it passes.
#' @noRd
checkTrxImportCounts <- function(assayList, idCol, annotationCol) {
    table <- assayList[[1L]]
    if (!omicsImportScalarString(annotationCol)) {
        trxAbort("a gene name column is required (trx_counts_invalid)", c("trx_counts_invalid", "trx_error"))
    }
    absent <- setdiff(c(idCol, annotationCol), names(table))
    if (length(absent)) {
        trxAbort(sprintf("the count table lacks the column(s) %s (trx_counts_invalid)", paste(absent, collapse = ", ")),
            c("trx_counts_invalid", "trx_error"))
    }
    withCallingHandlers(validateTrxCounts(as.data.frame(table)[setdiff(names(table), c(idCol, annotationCol))]),
        trx_counts_invalid = \(error) trxAbort(sprintf("%s (%s)", conditionMessage(error), class(error)[[1L]]),
            setdiff(class(error), c("rlang_error", "error", "condition"))))
    invisible(assayList)
}
