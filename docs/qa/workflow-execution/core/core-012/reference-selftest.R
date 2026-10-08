# CORE-012 exact-reference self-test. Usage: Rscript reference-selftest.R <isolated library> <receipt path>
# The reference library is test-only: never on the production path, never in DESCRIPTION.
args <- commandArgs(TRUE)
lib <- normalizePath(args[[1L]], mustWork = TRUE)
suppressPackageStartupMessages({
    library(gmp, lib.loc = lib)
    library(Rmpfr, lib.loc = lib)
})
# Exact sum: every finite double is an exact rational; add the rationals exactly.
exactSum <- function(x) Reduce(`+`, gmp::as.bigq(x), gmp::as.bigq(0))
relativeDeviation <- function(approx, exact) {
    e <- Rmpfr::mpfr(exact, precBits = 256L)
    a <- Rmpfr::mpfr(gmp::as.bigq(approx), precBits = 256L)
    as.numeric(abs(a - e) / abs(e))
}
checks <- list(
    # R's long-double sum keeps 1e16 + 1 but loses the 1 beyond its 64-bit mantissa; the reference keeps it.
    cancellation = identical(as.character(exactSum(c(1e20, 1, -1e20))), "1") && sum(c(1e20, 1, -1e20)) != 1,
    tenth_is_exact_rational = identical(as.character(gmp::as.bigq(0.1)), "3602879701896397/36028797018963968"),
    order_independent = identical(exactSum(c(0.1, 0.2, 0.3)), exactSum(c(0.3, 0.1, 0.2))),
    resolves_one_ulp = {
        x <- 1 + 2^-52
        d <- relativeDeviation(1, gmp::as.bigq(x))
        d > 2^-53 && d < 2^-51
    },
    tolerance_scale = relativeDeviation(1 + 2^-50, gmp::as.bigq(1)) == 2^-50,
    zero_exact = identical(as.character(exactSum(c(-2.5, 2.5))), "0"))
receipt <- list(schema = 1L, ticket = "CORE-012", created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    library = lib, r_version = R.version.string,
    packages = list(gmp = as.character(packageVersion("gmp", lib.loc = lib)),
        Rmpfr = as.character(packageVersion("Rmpfr", lib.loc = lib))),
    system = list(gmp = trimws(system2("pkg-config", c("--modversion", "gmp"), stdout = TRUE)),
        mpfr = trimws(system2("pkg-config", c("--modversion", "mpfr"), stdout = TRUE))),
    method = "exact rational sum of the doubles (gmp::as.bigq); relative deviation in 256-bit MPFR",
    self_test = lapply(checks, isTRUE), self_test_passed = all(vapply(checks, isTRUE, logical(1))),
    source = "CRAN https://cloud.r-project.org, installed 2026-10-08 into the isolated library")
jsonlite::write_json(receipt, args[[2L]], auto_unbox = TRUE, pretty = TRUE)
cat("self_test_passed:", receipt$self_test_passed, "\n"); print(unlist(receipt$self_test))
if (!receipt$self_test_passed) quit(status = 1L)
