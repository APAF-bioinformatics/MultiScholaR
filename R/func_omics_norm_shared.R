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
