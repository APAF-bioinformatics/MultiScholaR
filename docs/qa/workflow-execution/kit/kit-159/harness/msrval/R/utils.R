# Small shared helpers (internal).

#' Raise a classed harness condition (KIT-159 7.12)
#'
#' @param msg Message (cli format).
#' @param class One of the error classes of section 7.12.
#' @param .envir Environment for cli interpolation.
#' @return Never returns.
#' @noRd
valAbort <- function(msg, class, .envir = parent.frame()) {
    cli::cli_abort(msg, class = c(class, "val_error"), .envir = .envir)
}

#' The fixed random-number generator of every resampling (KIT-159 4.4)
#'
#' @param seed Integer seed.
#' @param code Code to evaluate.
#' @return The value of `code`.
#' @noRd
valWithSeed <- function(seed, code) {
    withr::with_seed(seed, code, .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion",
                     .rng_sample_kind = "Rejection")
}

#' Rows of a scored table that enter a metric: scored truth, key status `ok`
#' @noRd
valScoredOnly <- function(scored) {
    keep <- scored$truth_status %in% "scored" & scored$key_status %in% "ok"
    scored[keep, , drop = FALSE]
}

#' Two-sided critical value of the t distribution at a confidence level
#' @noRd
valTCritical <- function(conf, df) {
    stats::qt(1 - (1 - conf) / 2, df)
}
