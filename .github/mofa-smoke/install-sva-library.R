# KIT-151 d (branch run-mofa-evidence-k151d), task ruv only: sva, which the package's estimated count calls
# (sva::num.sv; R/func_general_ruv_factor_count.R at e8a2ac67), installed with its dependencies into the library of the
# R that runs this script, from the same dated snapshot as install-mofa-library.R (the snapshot the bundle installs
# from). sva must have the version the package's lockfile at e8a2ac67 records (3.60.0; the lockfile is copied to
# docs/qa/workflow-execution/kit/kit-151d/ref/lock-e8a2ac67.json), otherwise the run stops; every other package of the
# closure is compared with the lockfile and the comparison is written out (a difference there is reported, not fatal:
# the lockfile was taken from a local library and the snapshot holds 365 of 394 repository packages at the lock).
# No MultiScholaR, no internal package, no token.
#
# Usage: Rscript --vanilla install-sva-library.R <checkout> <output directory>   (jsonlite on the library path)

PPM <- "https://packagemanager.posit.co"
BIOC_VERSION <- "3.23"
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 2L)
checkout <- normalizePath(arguments[[1L]], winslash = "/", mustWork = TRUE)
out_dir <- arguments[[2L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
lib <- .Library
options(timeout = 1800, Ncpus = 4L, warn = 1)

lock <- jsonlite::fromJSON(file.path(checkout, "docs", "qa", "workflow-execution", "kit", "kit-151d", "ref",
    "lock-e8a2ac67.json"), simplifyVector = FALSE)
snapshot <- lock$Snapshot$Date
bioc_dated <- sprintf("%s/bioconductor/%s/packages/%s", PPM, snapshot, BIOC_VERSION)
source_repos <- c(CRAN = sprintf("%s/cran/%s", PPM, snapshot),
    BioCsoft = paste0(bioc_dated, "/bioc"), BioCann = paste0(bioc_dated, "/data/annotation"))
bioc_root <- sprintf("%s/bioconductor/__linux__/jammy/%s/packages/%s", PPM, snapshot, BIOC_VERSION)
install_repos <- c(CRAN = sprintf("%s/cran/__linux__/jammy/%s", PPM, snapshot),
    BioCsoft = paste0(bioc_root, "/bioc"), BioCann = paste0(bioc_root, "/data/annotation"))
options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(),
    paste(getRversion(), R.version[["platform"]], R.version[["arch"]], R.version[["os"]])))

db <- utils::available.packages(repos = source_repos, type = "source")
closure <- unique(c("sva", unlist(tools::package_dependencies("sva", db = db,
    which = c("Depends", "Imports", "LinkingTo"), recursive = TRUE), use.names = FALSE)))
base_packages <- rownames(utils::installed.packages(lib.loc = .Library, priority = "base"))
closure <- sort(setdiff(closure, base_packages))
present <- rownames(utils::installed.packages(lib.loc = lib))
needed <- setdiff(closure, present)
cat(sprintf("snapshot %s; sva closure %d packages; %d already in the library; installing %s\n", snapshot,
    length(closure), length(closure) - length(needed), paste(needed, collapse = ", ")))
if (length(needed)) utils::install.packages(needed, lib = lib, repos = install_repos, type = "source",
    dependencies = FALSE, INSTALL_opts = "--no-test-load")

installed <- utils::installed.packages(lib.loc = lib)
missing <- setdiff(closure, rownames(installed))
if (length(missing)) stop("not installed: ", paste(missing, collapse = ", "))
locked <- vapply(closure, function(p) {
    entry <- lock$Packages[[p]]
    if (is.null(entry) || is.null(entry$Version)) NA_character_ else as.character(entry$Version)
}, character(1))
versions <- data.frame(package = closure, installed = unname(installed[closure, "Version"]),
    lockfile_e8a2ac67 = unname(locked), stringsAsFactors = FALSE)
versions$same <- !is.na(versions$lockfile_e8a2ac67) & versions$installed == versions$lockfile_e8a2ac67
utils::write.table(versions, file.path(out_dir, "sva-library-versions.tsv"), sep = "\t", quote = FALSE,
    row.names = FALSE)
print(versions, row.names = FALSE)
if (!isTRUE(versions$same[versions$package == "sva"])) stop("sva differs from the lockfile at e8a2ac67")
if (!all(versions$same)) cat(sprintf("note: %d of %d closure packages differ from the lockfile (see the table)\n",
    sum(!versions$same), nrow(versions)))
cat("sva library ready\n")
