# Our NIPALS imputation returns exactly what mixOmics::impute.nipals returns, without its dense feature-by-feature matrix.
library(testthat)

nipalsMatrix <- function(n_samples, n_features, seed, missing) {
  set.seed(seed)
  values <- matrix(rnorm(n_samples * n_features, 15, 2), n_samples, n_features,
    dimnames = list(sprintf("S%03d", seq_len(n_samples)), sprintf("F%05d", seq_len(n_features)))) +
    rep(rnorm(n_samples, 0, 0.4), times = n_features)
  missing(values)
}

clustered <- function(values) { # the lowest-abundance 30% of features, in 15% of their samples
  sparse <- order(colMeans(values))[seq_len(round(ncol(values) * 0.3))]
  for (feature in sparse) values[sample.int(nrow(values), max(1L, round(nrow(values) * 0.15))), feature] <- NA_real_
  values
}

uniform <- function(values) {
  values[sample.int(length(values), floor(length(values) * 0.05))] <- NA_real_
  values
}

reference <- function(X, ...) suppressWarnings(suppressMessages(mixOmics::impute.nipals(X, ...)))

test_that("imputation is bit-identical to mixOmics across patterns, sizes and components", {
  for (case in list(
    list(24L, 300L, clustered, 2L), list(24L, 300L, uniform, 2L),
    list(24L, 300L, clustered, 5L), list(200L, 1000L, clustered, 2L), list(200L, 1000L, uniform, 2L),
    list(300L, 60L, clustered, 2L), list(6L, 5L, uniform, 2L))) {
    X <- nipalsMatrix(case[[1L]], case[[2L]], seed = case[[2L]] + case[[4L]], missing = case[[3L]])
    expect_true(anyNA(X))
    expect_identical(MultiScholaR:::imputeNipalsMatrix(X, ncomp = case[[4L]]), reference(X, ncomp = case[[4L]]))
  }
})

test_that("observed values are returned unchanged and every gap is filled", {
  X <- nipalsMatrix(24L, 300L, 7L, clustered)
  imputed <- MultiScholaR:::imputeNipalsMatrix(X, ncomp = 2L)
  expect_identical(imputed[!is.na(X)], X[!is.na(X)])
  expect_false(anyNA(imputed))
  expect_identical(dimnames(imputed), dimnames(X))
})

test_that("a run stopped by the iteration limit is identical too", {
  X <- nipalsMatrix(24L, 300L, 8L, uniform)
  expect_identical(suppressMessages(MultiScholaR:::imputeNipalsMatrix(X, ncomp = 2L, max_iter = 2L)), reference(X, ncomp = 2L, max.iter = 2L))
  expect_identical(MultiScholaR:::imputeNipalsMatrix(X, ncomp = 2L, tol = 1e-12), reference(X, ncomp = 2L, tol = 1e-12))
})

test_that("complete data comes back untouched, and impossible input is refused as mixOmics refuses it", {
  X <- nipalsMatrix(24L, 300L, 9L, identity)
  expect_identical(MultiScholaR:::imputeNipalsMatrix(X, ncomp = 2L), X)
  one_component <- nipalsMatrix(24L, 300L, 13L, clustered)
  expect_error(reference(one_component, ncomp = 1L))
  expect_error(MultiScholaR:::imputeNipalsMatrix(one_component, ncomp = 1L), "between 2")
  empty_row <- nipalsMatrix(24L, 300L, 10L, uniform); empty_row[3L, ] <- NA_real_
  expect_error(reference(empty_row, ncomp = 2L))
  expect_error(MultiScholaR:::imputeNipalsMatrix(empty_row, ncomp = 2L), "entirely missing")
  empty_column <- nipalsMatrix(24L, 300L, 11L, uniform); empty_column[, 3L] <- NA_real_
  expect_error(MultiScholaR:::imputeNipalsMatrix(empty_column, ncomp = 2L), "entirely missing")
})

test_that("no feature-by-feature matrix is built", {
  X <- nipalsMatrix(50L, 4000L, 12L, clustered)
  gc(reset = TRUE)
  MultiScholaR:::imputeNipalsMatrix(X, ncomp = 2L)
  used_mb <- sum(gc()[, 6L])
  gc(reset = TRUE)
  reference(X, ncomp = 2L)
  reference_mb <- sum(gc()[, 6L])
  # A 4000 x 4000 double matrix is about 122 MB: mixOmics carries it, and we must not.
  expect_gt(reference_mb - used_mb, 100)
})

referencePca <- function(X, ...) suppressWarnings(suppressMessages(mixOmics::pca(X, ...)))

expectSamePca <- function(ours, theirs) {
  expect_identical(ours$variates$X, theirs$variates$X)
  expect_identical(ours$prop_expl_var$X, theirs$prop_expl_var$X)
  expect_identical(ours$sdev, theirs$sdev)
  expect_identical(ours$loadings$X, theirs$loadings$X)
  expect_identical(ours$var.tot, theirs$var.tot)
  expect_identical(ours$cum.var, theirs$cum.var)
}

test_that("PCA with missing values is bit-identical to mixOmics in everything the plots read", {
  for (case in list(list(24L, 300L, clustered, 2L), list(24L, 300L, uniform, 2L), list(24L, 300L, clustered, 1L),
    list(24L, 300L, clustered, 5L), list(200L, 1000L, clustered, 2L), list(300L, 60L, uniform, 2L))) {
    X <- nipalsMatrix(case[[1L]], case[[2L]], seed = 100L + case[[2L]] + case[[4L]], missing = case[[3L]])
    expectSamePca(MultiScholaR:::pcaWithMissing(X, ncomp = case[[4L]]), referencePca(X, ncomp = case[[4L]]))
  }
  X <- nipalsMatrix(24L, 300L, 31L, clustered) # the sites that take the default number of components
  expectSamePca(MultiScholaR:::pcaWithMissing(X), referencePca(X))
})

test_that("PCA on complete data goes to mixOmics untouched", {
  X <- nipalsMatrix(24L, 300L, 32L, identity)
  calls <- 0L
  original <- mixOmics::pca
  local_mocked_bindings(pca = function(...) { calls <<- calls + 1L; original(...) }, .package = "mixOmics")
  ours <- MultiScholaR:::pcaWithMissing(X, ncomp = 2L)
  expect_identical(calls, 1L)
  expect_s3_class(ours, "pca")
})

test_that("PCA with missing values builds no feature-by-feature matrix, and refuses what mixOmics refuses", {
  X <- nipalsMatrix(50L, 4000L, 33L, clustered)
  gc(reset = TRUE)
  MultiScholaR:::pcaWithMissing(X, ncomp = 2L)
  used_mb <- sum(gc()[, 6L])
  gc(reset = TRUE)
  referencePca(X, ncomp = 2L)
  expect_gt(sum(gc()[, 6L]) - used_mb, 100)
  empty_row <- nipalsMatrix(24L, 300L, 34L, uniform); empty_row[3L, ] <- NA_real_
  expect_error(referencePca(empty_row, ncomp = 2L))
  expect_error(MultiScholaR:::pcaWithMissing(empty_row, ncomp = 2L), "entirely missing")
  too_many <- nipalsMatrix(6L, 5L, 35L, uniform)
  expect_error(referencePca(too_many, ncomp = 9L), "smaller")
  expect_error(MultiScholaR:::pcaWithMissing(too_many, ncomp = 9L), "smaller")
  infinite <- nipalsMatrix(24L, 300L, 36L, uniform); infinite[1L, 1L] <- Inf
  expect_error(referencePca(infinite, ncomp = 2L))
  expect_error(MultiScholaR:::pcaWithMissing(infinite, ncomp = 2L), "finite")
})
