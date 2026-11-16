# KIT-130 packaging spike: helpers shared by the check scripts. Sourced, not a package.

# The folder the bundle is expected to be running from, given by the caller of the check.
expectedBundle <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) < 1L) {
        stop("usage: Rscript --vanilla <check>.R <bundle directory>")
    }
    normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
}

# TRUE for each path that lies inside the bundle (long names, case-insensitive on Windows).
insideBundle <- function(paths, bundle) {
    canon <- \(p) {
        p <- normalizePath(p, winslash = "/", mustWork = FALSE)
        if (.Platform$OS.type == "windows") tolower(p) else p
    }
    startsWith(canon(paths), paste0(canon(bundle), "/"))
}

# The programs of these names found on the PATH that are not the bundle's own. (R for Windows
# puts its own bin directory on the PATH of its process, so finding the bundled one is expected.)
foreignOnPath <- function(programs, bundle) {
    found <- Sys.which(programs)
    found <- found[nzchar(found)]
    unname(found[!insideBundle(found, bundle)])
}

# TRUE when a TCP connection to the outside can be opened from this R process.
canReachNetwork <- function(host = "cran.r-project.org", port = 80L) {
    tryCatch(
        {
            con <- socketConnection(host, port, open = "r+", blocking = TRUE, timeout = 10)
            close(con)
            TRUE
        },
        error = \(e) FALSE,
        warning = \(w) FALSE
    )
}

checkPassed <- function(name) {
    cat(sprintf("KIT-130 CHECK PASSED: %s\n", name))
}
