# run-tests: packages the snapshot library lacks (params extra_packages), from the same snapshot as the library:
# CRAN and Bioconductor 3.23 jammy binaries of Posit Package Manager at the lockfile's snapshot date (as the KIT-124
# workflow installs sva), with their hard dependencies, into the library (.Library).
#
# Usage: Rscript --vanilla extra-packages.R <lockfile> <package>...

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) < 2L) stop("usage: Rscript --vanilla extra-packages.R <lockfile> <package>...")
snapshot <- jsonlite::fromJSON(arguments[[1L]])$Snapshot$Date
packages <- arguments[-1L]
options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(),
    paste(getRversion(), R.version[["platform"]], R.version[["arch"]], R.version[["os"]])))
ppm <- "https://packagemanager.posit.co"
bioc <- sprintf("%s/bioconductor/__linux__/jammy/%s/packages/3.23", ppm, snapshot)
repos <- c(CRAN = sprintf("%s/cran/__linux__/jammy/%s", ppm, snapshot), BioCsoft = paste0(bioc, "/bioc"),
    BioCann = paste0(bioc, "/data/annotation"))
available <- utils::available.packages(repos = repos)
unknown <- setdiff(packages, rownames(available))
if (length(unknown)) {
    cat(sprintf("::error::not in the snapshot %s (CRAN, Bioconductor 3.23): %s\n", snapshot, paste(unknown, collapse = ", ")))
    quit(save = "no", status = 1L)
}
wanted <- unique(c(packages, unlist(tools::package_dependencies(packages, db = available,
    which = c("Depends", "Imports", "LinkingTo"), recursive = TRUE))))
needed <- setdiff(wanted, rownames(utils::installed.packages(lib.loc = .Library)))
cat(sprintf("snapshot %s; not in the library: %s\n", snapshot, if (length(needed)) paste(needed, collapse = ", ") else "none"))
if (length(needed)) utils::install.packages(needed, lib = .Library, repos = repos, dependencies = FALSE)
missing <- setdiff(packages, rownames(utils::installed.packages(lib.loc = .Library)))
if (length(missing)) {
    cat(sprintf("::error::not installed: %s\n", paste(missing, collapse = ", ")))
    quit(save = "no", status = 1L)
}
for (p in packages) cat(p, as.character(utils::packageVersion(p, lib.loc = .Library)), "\n")
