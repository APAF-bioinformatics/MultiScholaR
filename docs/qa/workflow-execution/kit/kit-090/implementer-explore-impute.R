# Black-box comparison of candidate imputation variants with mixOmics::impute.nipals outputs.
source("load-candidate.R"); load_candidate()
nipalsMatrix <- function(n_samples, n_features, seed, missing) {
  set.seed(seed)
  values <- matrix(rnorm(n_samples * n_features, 15, 2), n_samples, n_features,
    dimnames = list(sprintf("S%03d", seq_len(n_samples)), sprintf("F%05d", seq_len(n_features)))) +
    rep(rnorm(n_samples, 0, 0.4), times = n_features)
  missing(values)
}
clustered <- function(values) {
  sparse <- order(colMeans(values))[seq_len(round(ncol(values) * 0.3))]
  for (feature in sparse) values[sample.int(nrow(values), max(1L, round(nrow(values) * 0.15))), feature] <- NA_real_
  values
}
uniform <- function(values) {
  values[sample.int(length(values), floor(length(values) * 0.05))] <- NA_real_
  values
}
cases <- list(list(24L, 300L, clustered, 2L), list(24L, 300L, uniform, 2L), list(24L, 300L, clustered, 5L),
  list(300L, 60L, clustered, 2L), list(6L, 5L, uniform, 2L))
data <- lapply(cases, function(case) {
  X <- nipalsMatrix(case[[1L]], case[[2L]], seed = case[[2L]] + case[[4L]], missing = case[[3L]])
  list(X = X, ncomp = case[[4L]], ref = suppressWarnings(suppressMessages(mixOmics::impute.nipals(X, ncomp = case[[4L]]))))
})
grid <- expand.grid(wnum = c("blas", "cs"), wden = c("blas", "cs"), norm = c("sum", "blas"),
  tnum = c("blas", "rs"), tden = c("blas", "rs"), stringsAsFactors = FALSE)
if (nzchar(Sys.getenv("GRID_FIX"))) eval(parse(text = Sys.getenv("GRID_FIX")))
for (g in seq_len(nrow(grid))) {
  do.call(options, setNames(as.list(grid[g, ]), paste0("nip.", names(grid))))
  res <- vapply(data, function(d) {
    ours <- suppressMessages(MultiScholaR:::imputeNipalsMatrix(d$X, ncomp = d$ncomp))
    miss <- is.na(d$X)
    sprintf("%d/%d %.1e", sum(ours[miss] != d$ref[miss]), sum(miss), max(abs(ours[miss] - d$ref[miss])))
  }, "")
  cat(paste(grid[g, ], collapse = " "), "|", paste(res, collapse = " | "), "\n")
}
