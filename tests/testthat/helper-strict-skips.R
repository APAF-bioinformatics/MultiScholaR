# KIT-101: the strict-skip sentinel. The file-level guards that skip a whole test file when named `R/` files are
# absent (13 guards in 12 files) or named internals have left the namespace (8 guards in 8 files), and the three
# partial loads that drop absent `R/` files without a skip, call these helpers instead of `testthat::skip()` and
# `paths[file.exists(paths)]`. By default (and under `R CMD check`) they behave exactly as before. With
# `MULTISCHOLAR_STRICT_SKIPS=1` the guard raises a classed error, so the file fails instead of disappearing from the
# count, and a partial load that would drop a file fails the same way.

#' @return TRUE when `MULTISCHOLAR_STRICT_SKIPS` is `"1"`.
#' @noRd
strictSkipsEnabled <- function() identical(Sys.getenv("MULTISCHOLAR_STRICT_SKIPS"), "1")

#' A file-level guard's skip: `testthat::skip(message)`, or under the sentinel a `multischolar_strict_skip` error.
#' @param message The guard's skip reason, unchanged.
#' @noRd
strictSkip <- function(message) {
  if (strictSkipsEnabled()) {
    rlang::abort(paste("strict skip (MULTISCHOLAR_STRICT_SKIPS=1):", message), class = "multischolar_strict_skip")
  }
  testthat::skip(message)
}

#' A partial load's file list: the paths that exist, as `paths[file.exists(paths)]`; under the sentinel an absent path
#' is a `multischolar_strict_skip` error naming it.
#' @param paths The `R/` files the test file sources.
#' @noRd
strictSkipPresent <- function(paths) {
  present <- file.exists(paths)
  if (!all(present) && strictSkipsEnabled()) {
    rlang::abort(paste("strict skip (MULTISCHOLAR_STRICT_SKIPS=1): partial load dropped",
      paste(basename(paths[!present]), collapse = ", ")), class = "multischolar_strict_skip")
  }
  paths[present]
}
