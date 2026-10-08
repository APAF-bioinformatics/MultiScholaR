#' Refuse with a classed error when an optional (Suggests) package is absent.
#' @param package Package name.
#' @param purpose Short description of what needs it, for the message.
#' @return TRUE invisibly.
#' @noRd
requireOptionalPackage <- function(package, purpose) {
    if (!optionalPackageAvailable(package)) {
        cli::cli_abort("{.pkg {package}} is required for {purpose} but is not installed.",
            class = "multischolar_package_missing")
    }
    invisible(TRUE)
}

#' Availability check kept separate so tests can simulate an absent package.
#' @noRd
optionalPackageAvailable <- function(package) requireNamespace(package, quietly = TRUE)
