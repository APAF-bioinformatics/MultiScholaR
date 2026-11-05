#' NIPALS components of a matrix with missing values, without a dense feature-by-feature matrix
#'
#' NIPALS (Wold, 1966) estimates each component's scores and loadings by alternating regressions over the observed
#' cells only. Each regression divides by the sum of squares over a column's or row's observed cells, which needs one
#' number per column or row, so the cost is one pass over the observed-cell mask. At 200 samples by 18,000 features
#' this saves about 5 GB and most of the time against `mixOmics::nipals`, whose result this equals bit for bit on the
#' package's reference BLAS (KIT-087; proof in `docs/qa/workflow-execution/kit/kit-087/`).
#'
#' @param X Numeric matrix with at least one missing value, samples in rows and features in columns.
#' @param ncomp Number of components.
#' @param max_iter,tol Iteration limit and convergence threshold for each component.
#' @return `list(scores, loadings, singular_values)`, the scores scaled to unit length.
#' @noRd
nipalsComponents <- function(X, ncomp, max_iter, tol) {
    missing_cells <- is.na(X)
    if (any(colSums(!missing_cells) == 0L) || any(rowSums(!missing_cells) == 0L)) {
        stop("Some rows or columns are entirely missing: remove them first.", call. = FALSE)
    }

    observed <- (!missing_cells) * 1
    observed_by_feature <- t(observed)
    n_features <- ncol(X)
    scores <- matrix(0, nrow(X), ncomp)
    loadings <- matrix(0, n_features, ncomp)
    singular_values <- numeric(ncomp)
    residual <- X
    start_column <- which.min(colSums(missing_cells))

    for (component in seq_len(ncomp)) {
        score <- residual[, start_column]
        score[is.na(score)] <- 0
        filled <- residual
        filled[missing_cells] <- 0
        previous <- 1 / rep(sqrt(n_features), n_features)
        iteration <- 1L
        change <- 1

        while (change > tol && iteration <= max_iter) {
            loading <- crossprod(filled, score) / drop(crossprod(score^2, observed))
            loading <- loading / drop(sqrt(crossprod(loading)))
            score <- (filled %*% loading) / drop(crossprod(loading^2, observed_by_feature))
            change <- sum((loading - previous)^2, na.rm = TRUE)
            previous <- loading
            iteration <- iteration + 1L
        }
        if (iteration > max_iter) {
            message(sprintf("NIPALS reached the iteration limit on component %d.", component))
        }

        residual <- residual - score %*% t(loading)
        loadings[, component] <- loading
        singular_values[[component]] <- sqrt(sum(score * score, na.rm = TRUE))
        scores[, component] <- score / singular_values[[component]]
    }

    list(scores = scores, loadings = loadings, singular_values = singular_values)
}

#' Impute missing values from a NIPALS reconstruction
#'
#' @param X Numeric matrix, samples in rows and features in columns, as the RUV canonical-correlation input is.
#' @param ncomp,max_iter,tol As `nipalsComponents`.
#' @return `X` with its missing values replaced from the `ncomp`-component reconstruction; observed values unchanged.
#' @noRd
imputeNipalsMatrix <- function(X, ncomp = 2L, max_iter = 500L, tol = 1e-06) {
    if (!anyNA(X)) {
        return(X)
    }
    # One component is refused because `mixOmics::impute.nipals` fails on it, so callers see the behaviour they always had.
    if (ncomp < 2L || ncomp > min(dim(X))) {
        stop("`ncomp` must be between 2 and the smaller dimension of `X`.", call. = FALSE)
    }
    components <- nipalsComponents(X, ncomp, max_iter, tol)
    reconstruction <- components$scores %*% diag(components$singular_values) %*% t(components$loadings)
    missing_cells <- is.na(X)
    X[missing_cells] <- reconstruction[missing_cells]
    X
}

#' Principal components for the QC plots
#'
#' Complete data goes to `mixOmics::pca` as before. With missing values `mixOmics::pca` runs its NIPALS, whose dense
#' matrix stopped a QC figure at the worker deadline at 200 samples by 18,000 features (KIT-086), so this computes the
#' same result from `nipalsComponents`, with the fields the plots read.
#' @param X Numeric matrix, samples in rows and features in columns.
#' @param ncomp Number of components, 2 as in `mixOmics::pca`.
#' @return The `mixOmics::pca` object for complete data; otherwise a list with its `X`, `ncomp`, `sdev`, `loadings`,
#'   `variates`, `prop_expl_var`, `var.tot` and `cum.var`.
#' @noRd
pcaWithMissing <- function(X, ncomp = 2L) {
    if (!anyNA(X)) {
        return(mixOmics::pca(X, ncomp = ncomp))
    }
    if (any(is.infinite(X))) {
        stop("`X` must hold finite values or missing values.", call. = FALSE)
    }
    ncomp <- round(ncomp)
    if (ncomp < 1L || ncomp > min(dim(X))) {
        stop("Use a smaller `ncomp`.", call. = FALSE)
    }

    centred <- scale(X, center = TRUE, scale = FALSE)
    components <- nipalsComponents(centred, ncomp, 500L, 1e-09)
    denominator <- max(1, nrow(X) - 1)
    pc_names <- paste0("PC", seq_len(ncomp))

    filled <- centred
    filled[is.na(filled)] <- 0
    loadings <- components$loadings
    variates <- filled %*% loadings
    dimnames(loadings) <- list(colnames(centred), pc_names)
    dimnames(variates) <- list(rownames(centred), pc_names)
    sdev <- stats::setNames(components$singular_values / sqrt(denominator), pc_names)

    var_tot <- sum(centred^2, na.rm = TRUE) / denominator
    proportion <- sdev^2 / var_tot
    list(X = centred, ncomp = ncomp, sdev = sdev, loadings = list(X = loadings), variates = list(X = variates),
        prop_expl_var = list(X = proportion), var.tot = var_tot, cum.var = cumsum(proportion))
}
