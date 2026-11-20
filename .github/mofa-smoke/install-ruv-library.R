# KIT-151 c (branch run-mofa-evidence-k151c), task 3 only: the RUV-III functions the package calls, installed into the
# library of the R that runs this script, at the bundle's versions. RUVIIIC (RUVIII_C_Varying: proteomics, lipidomics,
# metabolomics) from the bundle's GitHub source cran/RUVIIIC@c534151f; ruv (ruv::RUVIII: transcriptomics) and the
# dependencies of both from the KIT-134 snapshot, as install-mofa-library.R installs MOFA2's closure. No MultiScholaR,
# no internal package, no token. Any version that differs from the bundle's (kit-134/p2/linux/bundle-sources.tsv)
# stops the run.
#
# Usage: Rscript --vanilla install-ruv-library.R <checkout> <output directory>   (jsonlite on the library path)

PPM <- "https://packagemanager.posit.co"
RUVIIIC_SHA <- "c534151f6b2265b36294ff630f4dd13c2e86a9fa"
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 2L)
checkout <- normalizePath(arguments[[1L]], winslash = "/", mustWork = TRUE)
out_dir <- arguments[[2L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
lib <- .Library
options(timeout = 1800, Ncpus = 4L, warn = 1)

lockfile <- jsonlite::fromJSON(file.path(checkout, "packaging", "multischolar.lock.json"), simplifyVector = FALSE)
snapshot <- lockfile$Snapshot$Date
source_repos <- c(CRAN = sprintf("%s/cran/%s", PPM, snapshot))
install_repos <- c(CRAN = sprintf("%s/cran/__linux__/jammy/%s", PPM, snapshot))
options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(),
    paste(getRversion(), R.version[["platform"]], R.version[["arch"]], R.version[["os"]])))

db <- utils::available.packages(repos = source_repos, type = "source")
roots <- c("ruv", "RSpectra", "progress", "Rcpp", "RcppEigen", "RcppProgress")
closure <- unique(c(roots, unlist(tools::package_dependencies(roots, db = db,
    which = c("Depends", "Imports", "LinkingTo"), recursive = TRUE), use.names = FALSE)))
base_packages <- rownames(utils::installed.packages(lib.loc = .Library, priority = "base"))
closure <- sort(setdiff(closure, base_packages))
present <- rownames(utils::installed.packages(lib.loc = lib))
needed <- setdiff(closure, present)
cat(sprintf("closure %d packages; %d already in the library; installing %s\n", length(closure),
    length(closure) - length(needed), paste(needed, collapse = ", ")))
if (length(needed)) utils::install.packages(needed, lib = lib, repos = install_repos, type = "source",
    dependencies = FALSE, INSTALL_opts = "--no-test-load")

tarball <- file.path(tempdir(), "RUVIIIC.tar.gz")
utils::download.file(sprintf("https://github.com/cran/RUVIIIC/archive/%s.tar.gz", RUVIIIC_SHA), tarball, mode = "wb",
    quiet = TRUE)
status <- system2(file.path(R.home("bin"), "R"), c("CMD", "INSTALL", "--no-test-load", paste0("--library=", lib),
    shQuote(tarball)))
if (!identical(as.integer(status), 0L)) stop("RUVIIIC did not install")

all_packages <- c(closure, "RUVIIIC")
installed <- utils::installed.packages(lib.loc = lib)
missing <- setdiff(all_packages, rownames(installed))
if (length(missing)) stop("not installed: ", paste(missing, collapse = ", "))
reference <- utils::read.delim(file.path(checkout, "docs", "qa", "workflow-execution", "kit", "kit-134", "p2", "linux",
    "bundle-sources.tsv"), colClasses = "character")
versions <- data.frame(package = all_packages, installed = unname(installed[all_packages, "Version"]),
    bundle = reference$version[match(all_packages, reference$package)], stringsAsFactors = FALSE)
versions$same <- !is.na(versions$bundle) & versions$installed == versions$bundle
utils::write.table(versions, file.path(out_dir, "ruv-library-versions.tsv"), sep = "\t", quote = FALSE,
    row.names = FALSE)
print(versions, row.names = FALSE)
if (!all(versions$same)) stop("installed versions differ from the bundle's")
cat("RUV library ready\n")
