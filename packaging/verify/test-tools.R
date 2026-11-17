# KIT-134: testthat and what the bundle lacks for it, from the bundle's snapshot date, into a
# separate library for packaging/verify/run-lanes.R (the bundle's own library is not written to).
# Needs the network; run with the bundle's Rscript.
#
# Usage: Rscript --vanilla test-tools.R <library to create> <snapshot date>

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) stop("usage: Rscript --vanilla test-tools.R <library to create> <snapshot date>")
lib <- arguments[[1L]]
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
windows <- identical(.Platform$OS.type, "windows")
# Posit Package Manager serves Linux binaries by the R user agent.
options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(),
    paste(getRversion(), R.version[["platform"]], R.version[["arch"]], R.version[["os"]])))
repos <- sprintf(if (windows) "https://packagemanager.posit.co/cran/%s" else "https://packagemanager.posit.co/cran/__linux__/jammy/%s",
    arguments[[2L]])
available <- utils::available.packages(repos = repos)
wanted <- unique(c("testthat", unlist(tools::package_dependencies("testthat", db = available,
    which = c("Depends", "Imports", "LinkingTo"), recursive = TRUE))))
needed <- setdiff(wanted, rownames(utils::installed.packages(lib.loc = .libPaths())))
cat("testthat needs, not in the bundle:", paste(needed, collapse = ", "), "\n")
utils::install.packages(needed, lib = lib, repos = repos, dependencies = FALSE, type = if (windows) "binary" else getOption("pkgType"))
.libPaths(c(.libPaths(), lib))
library(testthat)
cat("testthat", as.character(utils::packageVersion("testthat")), "from", find.package("testthat"), "\n")
