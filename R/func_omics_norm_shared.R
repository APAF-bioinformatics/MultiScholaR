# Shared between-sample normalisation rules of the lipidomics and metabolomics S4 methods.

#' The fewest features a sample must have measured for cyclic loess to be fitted (KIT-145).
#'
#' Decision M17 (user, 2026-11-18; docs/qa/workflow-execution/kit/kit-124/design-review-20261118.md, section 6).
#' 10 is a judgement about when a fitted curve means something, not a measurement. The numerical minimum is 4: with 2
#' or 3 features `limma::normalizeCyclicLoess` is a line through two points, its robustness passes weight points by
#' rounding-level residuals, and the result differs between machines (KIT-134, Status 2026-11-18).
#' @noRd
OMICS_CYCLIC_LOESS_MIN_FEATURES <- 10L

#' Whether the cyclic loess minimum is enforced (KIT-145).
#'
#' Off by default until the single re-record: the committed metabolomics entry-lane fixtures keep 3 features per assay
#' and their stored references were frozen on today's behaviour. The switch commit changes the default below to `TRUE`.
#' The option is read in the process that runs the method; it does not reach a managed worker, which takes the default.
#' @return `TRUE` when `options(MultiScholaR.cyclic_loess_minimum = TRUE)`.
#' @noRd
omicsCyclicLoessMinimumActive <- function() {
    isTRUE(getOption("MultiScholaR.cyclic_loess_minimum", FALSE))
}

#' Refuse cyclic loess on an assay with too few features (KIT-145).
#'
#' "Features at the point of the fit": `limma::normalizeCyclicLoess` (method "fast") fits one curve per sample through
#' the features that have a finite value in that sample (`limma::loessFit` drops the others). The count is therefore
#' taken per sample on the matrix handed to limma (after the filters, ITSD and log2 steps, non-finite values as
#' missing), and the sample with the fewest decides. A zero is a value and counts; a missing value does not.
#'
#' @param assay_matrix Numeric matrix, features by samples, as it would be passed to `limma::normalizeCyclicLoess`.
#' @param assay_name Name of the assay, for the message.
#' @param feature_noun Plural noun for the omic's features ("metabolites", "lipids").
#' @return `TRUE`, invisibly; an error of class `omics_cyclic_loess_too_few_features` (fields `assay`, `sample`,
#'   `fewest`, `minimum`) when a sample has fewer measured features than `OMICS_CYCLIC_LOESS_MIN_FEATURES`.
#' @noRd
checkOmicsCyclicLoessMinimum <- function(assay_matrix, assay_name, feature_noun) {
    if (!omicsCyclicLoessMinimumActive()) {
        return(invisible(TRUE))
    }
    measured <- colSums(is.finite(assay_matrix))
    fewest <- as.integer(min(measured))
    minimum <- OMICS_CYCLIC_LOESS_MIN_FEATURES
    if (fewest >= minimum) {
        return(invisible(TRUE))
    }
    sample_name <- colnames(assay_matrix)[[which.min(measured)]]
    # rlang::abort, not cli::cli_abort: cli wraps the text at the console width, and the tab shows it in a notification.
    rlang::abort(
        c(
            sprintf(
                paste(
                    "Cyclic loess normalisation was not run: assay \"%s\" has only %d %s with a measured value in",
                    "sample \"%s\", and the method needs at least %d in every sample."
                ),
                assay_name, fewest, feature_noun, sample_name, minimum
            ),
            i = sprintf(
                paste(
                    "Cyclic loess fits a curve through the %s of each sample. With so few points the curve is not",
                    "meaningful, and its result can differ from one computer to another."
                ),
                feature_noun
            ),
            i = paste(
                "Choose another between-sample method under 'Normalization Method' (Quantile, Scale (Median Absolute Values)",
                "or None) and run the normalisation again."
            )
        ),
        class = "omics_cyclic_loess_too_few_features",
        assay = assay_name,
        sample = sample_name,
        fewest = fewest,
        minimum = minimum,
        call = NULL
    )
}

# ---------------------------------------------------------------------------------------------------------------------
# KIT-242 (VAL-242): between-sample normalisation of an assay that holds a sample with no value, a method failure and a
# malformed assay. Older records replay under the rule they were made with (`workflowReplayRuleMarkers`).
# ---------------------------------------------------------------------------------------------------------------------

#' Whether the lipidomics and metabolomics between-sample normalisation runs under KIT-242's rule (VAL-242 OD-1 (c)):
#' TRUE unless the replay of a record made before sets the marker's earlier token ("warn": the whole-assay skip with a
#' warning, as then).
#' @noRd
omicsNormalisationCurrentRule <- function() {
    marker <- workflowReplayRuleMarkers$normalisation_fallback
    !identical(getOption(marker$option), marker$earlier)
}

#' Stop a between-sample normalisation (VAL-242: a method failure stops; a malformed assay stops).
#' @param assay The assay's name.
#' @param reason `normaliser_error`, `malformed_assay` or `invalid_method`.
#' @param method The method.
#' @param text The message.
#' @param parent The method's own error, for `normaliser_error`.
#' @return Does not return: an error of class `omics_normalisation_failed` (fields `assay`, `sample`, `method`,
#'   `reason`).
#' @noRd
abortOmicsNormalisation <- function(assay, reason, method, text, parent = NULL) {
    rlang::abort(text, class = "omics_normalisation_failed", assay = assay, sample = NULL, method = method, reason = reason,
        parent = parent, call = NULL)
}

#' Normalise one assay's matrix between samples under KIT-242's rule (VAL-242 OD-1 (c)): a sample whose column holds no
#' finite value is left out of the method (and out of KIT-145's count, item viii) and put back as missing, so the other
#' samples are normalised exactly as the method normalises the assay without it; a failure of the method stops. With no
#' such sample the method receives the assay's matrix unchanged and its result is returned as it is.
#' @param assay_matrix Numeric matrix, features by samples, non-finite values already missing.
#' @param assay_name The assay's name.
#' @param method `cyclicloess`, `quantile`, `scale` or `none`.
#' @param feature_noun Plural noun for the omic's features, for KIT-145's refusal.
#' @return `list(matrix, excluded)`: the normalised matrix and the samples left out (each signalled by a warning of class
#'   `omics_normalisation_sample_excluded`, fields `assay`, `sample`, `method`).
#' @noRd
normaliseOmicsAssayMatrix <- function(assay_matrix, assay_name, method, feature_noun) {
    if (identical(method, "none")) {
        message("   Normalization method is 'none'. Skipping application.")
        return(list(matrix = assay_matrix, excluded = character()))
    }
    normaliser <- switch(method, cyclicloess = limma::normalizeCyclicLoess, quantile = limma::normalizeQuantiles,
        scale = limma::normalizeMedianAbsValues, NULL)
    if (is.null(normaliser)) {
        abortOmicsNormalisation(assay_name, "invalid_method", method,
            sprintf("Between-sample normalisation was not run: '%s' is not a method (cyclicloess, quantile, scale or none).", method))
    }
    empty <- colnames(assay_matrix)[colSums(!is.na(assay_matrix)) == 0]
    kept <- setdiff(colnames(assay_matrix), empty)
    if (!length(kept)) {
        abortOmicsNormalisation(assay_name, "malformed_assay", method,
            sprintf("Between-sample normalisation was not run: assay '%s' has no value in any sample.", assay_name))
    }
    input <- if (length(empty)) assay_matrix[, kept, drop = FALSE] else assay_matrix
    if (identical(method, "cyclicloess")) checkOmicsCyclicLoessMinimum(input, assay_name, feature_noun)
    message(sprintf("   Applying %s normalization...", method))
    fitted <- tryCatch(normaliser(input), error = function(e) {
        abortOmicsNormalisation(assay_name, "normaliser_error", method, sprintf(paste("Between-sample normalisation (%s) failed on",
            "assay '%s': %s. The step was not completed; no assay was left unnormalised."), method, assay_name, conditionMessage(e)),
            parent = e)
    })
    if (!length(empty)) return(list(matrix = fitted, excluded = character()))
    out <- assay_matrix
    out[, kept] <- fitted
    out[, empty] <- NA_real_
    rlang::warn(sprintf(paste("Assay '%s': samples %s have no value in this assay. The assay was normalised (%s) without",
        "them, and they are left empty in this assay only."), assay_name, paste(empty, collapse = ", "), method),
        class = "omics_normalisation_sample_excluded", assay = assay_name, sample = empty, method = method)
    list(matrix = out, excluded = empty)
}

#' Remove the samples excluded from every assay (KIT-242, item 1a: a sample with no value in any assay leaves the design,
#' with a notice).
#' @param assays Named list of assay tables.
#' @param design The design matrix.
#' @param sample_id The design's sample column.
#' @param excluded Named list, per assay, of the samples excluded there.
#' @return `list(assays, design, removed)`; each removal is signalled by a warning of class
#'   `omics_normalisation_sample_removed` (field `sample`).
#' @noRd
omicsRemoveSamplesEmptyEverywhere <- function(assays, design, sample_id, excluded) {
    if (!length(assays)) return(list(assays = assays, design = design, removed = character()))
    removed <- Reduce(intersect, lapply(names(assays), function(assay) excluded[[assay]] %||% character()))
    if (!length(removed)) return(list(assays = assays, design = design, removed = character()))
    assays <- lapply(assays, function(table) table[, setdiff(names(table), removed), drop = FALSE])
    design <- design[!as.character(design[[sample_id]]) %in% removed, , drop = FALSE]
    rlang::warn(sprintf("Samples %s have no value in any assay and were removed from the design.", paste(removed, collapse = ", ")),
        class = "omics_normalisation_sample_removed", sample = removed)
    list(assays = assays, design = design, removed = removed)
}
