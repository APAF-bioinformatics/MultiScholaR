# NIPALS principal components, imputation and PCA for matrices with missing values.
#
# Method: H. Wold (1966), "Estimation of principal components and related models by iterative least squares", in
# P. R. Krishnaiah (ed.), Multivariate Analysis, Academic Press, pp. 391-420.
#
# Written in a clean room from the KIT-090 behaviour specification and the published method, checked against
# mixOmics::impute.nipals and mixOmics::pca as black boxes only. Every sum inside the iteration is a BLAS product
# (crossprod or %*%), so the results match the reference bit for bit. No features-by-features or samples-by-samples
# matrix is built: the largest intermediates have the input's own dimensions.

#' NIPALS components over the observed cells of a matrix
#'
#' @param X Numeric matrix, samples in rows and features in columns, `NA` for missing cells.
#' @param ncomp Number of components.
#' @param max_iter Iteration limit per component.
#' @param tol Convergence threshold on the squared change of the loading vector.
#' @return `list(scores, loadings, singular_values)`: an `n x ncomp` matrix of unit-length scores, a `p x ncomp`
#'   matrix of unit-length loadings and the `ncomp` singular values.
#' @noRd
nipalsComponents <- function(X, ncomp, max_iter, tol) {
  observed <- !is.na(X)
  if (any(colSums(observed) == 0L) || any(rowSums(observed) == 0L)) {
    stop("Some rows or columns are entirely missing: remove them first.", call. = FALSE)
  }
  n <- nrow(X)
  p <- ncol(X)
  observed_01 <- observed * 1
  start <- which.min(colSums(!observed))
  scores <- matrix(0, n, ncomp)
  loadings <- matrix(0, p, ncomp)
  singular_values <- numeric(ncomp)
  residual <- X

  for (h in seq_len(ncomp)) {
    t <- residual[, start]
    t[is.na(t)] <- 0
    residual_0 <- residual
    residual_0[!observed] <- 0
    previous <- rep(1 / sqrt(p), p)
    k <- 1L
    change <- 1
    while (change > tol && k <= max_iter) {
      w <- drop(crossprod(residual_0, t)) / drop(crossprod(observed_01, t^2))
      w <- w / drop(sqrt(crossprod(w)))
      t <- drop(residual_0 %*% w) / drop(observed_01 %*% w^2)
      change <- sum((w - previous)^2, na.rm = TRUE)
      previous <- w
      k <- k + 1L
    }
    if (k > max_iter) message("NIPALS reached the iteration limit on component ", h, ".")
    residual <- residual - tcrossprod(t, w)
    loadings[, h] <- w
    singular_values[h] <- sqrt(sum(t^2, na.rm = TRUE))
    scores[, h] <- t / singular_values[h]
  }
  list(scores = scores, loadings = loadings, singular_values = singular_values)
}

#' Fill missing values from a rank-`ncomp` NIPALS reconstruction
#'
#' @param X Numeric matrix, samples in rows and features in columns, `NA` for missing cells.
#' @param ncomp Number of components, between 2 and `min(dim(X))`.
#' @param max_iter Iteration limit per component.
#' @param tol Convergence threshold.
#' @return `X` with its missing cells filled; observed cells, dimnames and attributes unchanged. Without missing values,
#'   `X` itself.
#' @noRd
imputeNipalsMatrix <- function(X, ncomp = 2L, max_iter = 500L, tol = 1e-06) {
  if (!anyNA(X)) {
    return(X)
  }
  if (ncomp < 2 || ncomp > min(dim(X))) {
    stop("`ncomp` must be between 2 and the smaller dimension of `X`.", call. = FALSE)
  }
  fit <- nipalsComponents(X, ncomp, max_iter, tol)
  reconstruction <- fit$scores %*% diag(fit$singular_values) %*% t(fit$loadings)
  missing <- is.na(X)
  X[missing] <- reconstruction[missing]
  X
}

#' Principal components for the QC plots, tolerating missing values
#'
#' @param X Numeric matrix, samples in rows and features in columns, `NA` for missing cells.
#' @param ncomp Number of components.
#' @return Without missing values, `mixOmics::pca(X, ncomp = ncomp)`. Otherwise a list with the fields of that result
#'   the plots read: `X`, `ncomp`, `sdev`, `loadings`, `variates`, `prop_expl_var`, `var.tot`, `cum.var`.
#' @noRd
pcaWithMissing <- function(X, ncomp = 2L) {
  if (!anyNA(X)) {
    return(mixOmics::pca(X, ncomp = ncomp))
  }
  if (any(is.infinite(X))) {
    stop("`X` must hold finite values or missing values.", call. = FALSE)
  }
  ncomp <- round(ncomp)
  if (ncomp < 1 || ncomp > min(dim(X))) {
    stop("Use a smaller `ncomp`.", call. = FALSE)
  }
  centred <- scale(X, center = TRUE, scale = FALSE)
  fit <- nipalsComponents(centred, ncomp, 500L, 1e-09)
  denominator <- max(1, nrow(X) - 1)
  components <- paste0("PC", seq_len(ncomp))

  loadings <- fit$loadings
  dimnames(loadings) <- list(colnames(centred), components)
  centred_0 <- centred
  centred_0[is.na(centred_0)] <- 0
  variates <- centred_0 %*% loadings
  dimnames(variates) <- list(rownames(centred), components)
  sdev <- fit$singular_values / sqrt(denominator)
  names(sdev) <- components
  var_tot <- sum(centred^2, na.rm = TRUE) / denominator
  proportion <- sdev^2 / var_tot

  list(X = centred, ncomp = ncomp, sdev = sdev, loadings = list(X = loadings), variates = list(X = variates),
    prop_expl_var = list(X = proportion), var.tot = var_tot, cum.var = cumsum(proportion))
}
