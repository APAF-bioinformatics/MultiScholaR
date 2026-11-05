#' Impute missing values by NIPALS principal components, without a dense feature-by-feature matrix
#'
#' NIPALS (Wold, 1966) estimates each component's scores and loadings by alternating regressions over the observed
#' cells only. Each regression divides by the sum of squares over a column's or row's observed cells, which needs one
#' number per column or row, so the cost is one pass over the observed-cell mask. At 200 samples by 18,000 features
#' this saves about 5 GB and most of the time against `mixOmics::impute.nipals`, whose result this equals bit for bit
#' on the package's reference BLAS (KIT-087; proof in `docs/qa/workflow-execution/kit/kit-087/`).
#'
#' @param X Numeric matrix, samples in rows and features in columns, as the RUV canonical-correlation input is.
#' @param ncomp Number of components.
#' @param max_iter,tol Iteration limit and convergence threshold for each component.
#' @return `X` with its missing values replaced from the `ncomp`-component reconstruction; observed values unchanged.
#' @noRd
imputeNipalsMatrix <- function(X, ncomp = 2L, max_iter = 500L, tol = 1e-06) {
    missing_cells <- is.na(X)
    if (!any(missing_cells)) {
        return(X)
    }
    if (any(colSums(!missing_cells) == 0L) || any(rowSums(!missing_cells) == 0L)) {
        stop("Some rows or columns are entirely missing: remove them before imputing.", call. = FALSE)
    }
    # One component is refused because `mixOmics::impute.nipals` fails on it, so callers see the behaviour they always had.
    if (ncomp < 2L || ncomp > min(dim(X))) {
        stop("`ncomp` must be between 2 and the smaller dimension of `X`.", call. = FALSE)
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

    reconstruction <- scores %*% diag(singular_values) %*% t(loadings)
    X[missing_cells] <- reconstruction[missing_cells]
    X
}
