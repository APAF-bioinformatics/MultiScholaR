# KIT-130 packaging spike, build time only: install the three hard packages and shiny (for the
# hello-world app) into the library of the R that runs this script. Nothing here runs in the
# shipped application. MultiScholaR itself is deliberately not installed.

BIOC_VERSION <- "3.23" # the Bioconductor release for R 4.6 (bioconductor.org/config.yaml, 2026-11-17)
RUVIIIC_SOURCE <- "https://cran.r-project.org/src/contrib/Archive/RUVIIIC/RUVIIIC_1.0.19.tar.gz"

lib <- .Library
is_windows <- .Platform$OS.type == "windows"
options(timeout = 900, Ncpus = 4L)

if (is_windows) {
    bioc_root <- sprintf("https://bioconductor.org/packages/%s", BIOC_VERSION)
    cran <- "https://cloud.r-project.org"
    pkg_type <- "binary"
} else {
    # Posit Public Package Manager serves binaries for jammy when R identifies itself.
    bioc_root <- sprintf(
        "https://packagemanager.posit.co/bioconductor/__linux__/jammy/latest/packages/%s",
        BIOC_VERSION
    )
    cran <- "https://packagemanager.posit.co/cran/__linux__/jammy/latest"
    pkg_type <- "source"
    options(HTTPUserAgent = sprintf(
        "R/%s R (%s)",
        getRversion(),
        paste(getRversion(), R.version[["platform"]], R.version[["arch"]], R.version[["os"]])
    ))
}
repos <- c(
    BioCsoft = paste0(bioc_root, "/bioc"),
    BioCann = paste0(bioc_root, "/data/annotation"),
    BioCexp = paste0(bioc_root, "/data/experiment"),
    CRAN = cran
)

# RUVIIIC was archived on CRAN in 2023, so its dependencies are named here and it is compiled
# from the archived source on both platforms.
from_repos <- c("limpa", "MOFA2", "shiny", "RSpectra", "progress", "Rcpp", "RcppEigen", "RcppProgress")
utils::install.packages(from_repos, lib = lib, repos = repos, type = pkg_type)
utils::install.packages(RUVIIIC_SOURCE, lib = lib, repos = NULL, type = "source")

installed <- utils::installed.packages(lib.loc = lib)
needed <- c("limpa", "RUVIIIC", "MOFA2", "shiny", "reticulate")
missing <- setdiff(needed, rownames(installed))
cat(sprintf("library: %s\npackages in the library: %d\n", lib, nrow(installed)))
print(installed[intersect(needed, rownames(installed)), c("Version", "Built"), drop = FALSE])
if (length(missing) > 0L) {
    stop("not installed: ", paste(missing, collapse = ", "))
}
