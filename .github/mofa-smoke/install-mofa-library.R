# MOFA2 smoke lane (run-mofa-smoke; KIT-060 risks 2 and 6, KIT-134): install MOFA2 and its recursive
# R dependencies into the library of the R that runs this script, from the KIT-134 snapshot and the
# same repositories packaging/r/install-library.R installs the bundle from, without the internal
# packages and without MultiScholaR (the smoke does not load either). Then check every installed
# package's version against the bundle's recorded sources for this platform
# (docs/qa/workflow-execution/kit/kit-134/p2/<platform>/bundle-sources.tsv, the same commit as the
# lockfile) and MOFA2's mofapy2 pins against packaging/python/requirements.txt (the same check as
# install-library.R). Any difference stops the run.
#
# Usage: Rscript --vanilla install-mofa-library.R <plan|install> <checkout> <output directory>
# jsonlite must be on the library path (the build-tool library).

PPM <- "https://packagemanager.posit.co"
BIOC_VERSION <- "3.23"
ROOTS <- c("MOFA2", "jsonlite") # jsonlite: the smoke writes its result with it (the bundle ships it)

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 3L || !arguments[[1L]] %in% c("plan", "install")) {
    stop("usage: Rscript --vanilla install-mofa-library.R <plan|install> <checkout> <output directory>")
}
mode <- arguments[[1L]]
checkout <- normalizePath(arguments[[2L]], winslash = "/", mustWork = TRUE)
out_dir <- arguments[[3L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
is_windows <- .Platform$OS.type == "windows"
lib <- .Library
options(timeout = 1800, Ncpus = 4L, warn = 1)

# Repositories: as packaging/r/install-library.R (lines 80-107 on kit-134 3c81b922).
lockfile <- jsonlite::fromJSON(file.path(checkout, "packaging", "multischolar.lock.json"), simplifyVector = FALSE)
snapshot <- lockfile$Snapshot$Date
windows_bioc <- lockfile$Snapshot$BioconductorWindows
bioc_dated <- sprintf("%s/bioconductor/%s/packages/%s", PPM, snapshot, BIOC_VERSION)
source_repos <- c(CRAN = sprintf("%s/cran/%s", PPM, snapshot),
    BioCsoft = paste0(bioc_dated, "/bioc"), BioCann = paste0(bioc_dated, "/data/annotation"))
if (is_windows) {
    bioc_root <- if (identical(windows_bioc, "snapshot")) bioc_dated else
        sprintf("https://bioconductor.org/packages/%s", BIOC_VERSION)
    install_repos <- c(CRAN = sprintf("%s/cran/%s", PPM, snapshot),
        BioCsoft = paste0(bioc_root, "/bioc"), BioCann = paste0(bioc_root, "/data/annotation"))
} else {
    bioc_root <- sprintf("%s/bioconductor/__linux__/jammy/%s/packages/%s", PPM, snapshot, BIOC_VERSION)
    install_repos <- c(CRAN = sprintf("%s/cran/__linux__/jammy/%s", PPM, snapshot),
        BioCsoft = paste0(bioc_root, "/bioc"), BioCann = paste0(bioc_root, "/data/annotation"))
    options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(),
        paste(getRversion(), R.version[["platform"]], R.version[["arch"]], R.version[["os"]])))
}
cat(sprintf("R %s; library %s; snapshot %s; mode %s\n", getRversion(), lib, snapshot, mode))
print(install_repos)

db <- utils::available.packages(repos = source_repos, type = "source")
base_packages <- rownames(utils::installed.packages(lib.loc = .Library, priority = "base"))
closure <- unique(c(ROOTS, unlist(tools::package_dependencies(ROOTS, db = db,
    which = c("Depends", "Imports", "LinkingTo"), recursive = TRUE), use.names = FALSE)))
closure <- sort(setdiff(closure, base_packages))
unresolved <- setdiff(closure, rownames(db))
if (length(unresolved)) stop("not in the snapshot: ", paste(unresolved, collapse = ", "))
# install-library.R builds these differently (GitHub, source, lockfile version): none may be in the closure.
special <- intersect(closure, c("RUVIIIC", "Glimma", "GO.db", "nanoarrow", "limpa"))
if (length(special)) stop("the closure holds packages the bundle builds another way: ", paste(special, collapse = ", "))
repository <- ifelse(startsWith(db[closure, "Repository"], source_repos[["CRAN"]]), "CRAN", "Bioconductor")
utils::write.table(data.frame(package = closure, repository = repository), file.path(out_dir, "plan-packages.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("closure of %s: %d packages\n", paste(ROOTS, collapse = ", "), length(closure)))
if (mode == "plan") quit(save = "no", status = 0L)

if (is_windows) {
    binaries <- utils::available.packages(repos = install_repos, type = "win.binary")
    as_binary <- intersect(closure, rownames(binaries))
    utils::install.packages(as_binary, lib = lib, repos = install_repos, type = "binary", dependencies = FALSE,
        INSTALL_opts = "--no-test-load")
    rest <- setdiff(closure, as_binary)
    if (length(rest)) utils::install.packages(rest, lib = lib, repos = install_repos, type = "source",
        dependencies = FALSE, INSTALL_opts = "--no-test-load")
} else {
    utils::install.packages(closure, lib = lib, repos = install_repos, type = "source", dependencies = FALSE,
        INSTALL_opts = "--no-test-load")
}

installed <- utils::installed.packages(lib.loc = lib)
missing <- setdiff(closure, rownames(installed))
if (length(missing)) stop("not installed: ", paste(missing, collapse = ", "))

# The bundle's versions for this platform.
platform <- if (is_windows) "windows" else "linux"
reference <- utils::read.delim(file.path(checkout, "docs", "qa", "workflow-execution", "kit", "kit-134", "p2",
    platform, "bundle-sources.tsv"), colClasses = "character")
versions <- data.frame(package = closure, installed = unname(installed[closure, "Version"]),
    bundle = reference$version[match(closure, reference$package)], stringsAsFactors = FALSE)
versions$same <- !is.na(versions$bundle) & versions$installed == versions$bundle
utils::write.table(versions, file.path(out_dir, "mofa-library-versions.tsv"), sep = "\t", quote = FALSE,
    row.names = FALSE)
cat(sprintf("against the bundle (%s): %d of %d at the bundle's version\n", platform, sum(versions$same),
    nrow(versions)))
if (!all(versions$same)) {
    print(versions[!versions$same, ], row.names = FALSE)
    stop("installed versions differ from the bundle's")
}

# MOFA2's mofapy2 pins against requirements.txt: as packaging/r/install-library.R lines 357-371.
requirements <- readLines(file.path(checkout, "packaging", "python", "requirements.txt"), warn = FALSE)
pins <- regmatches(requirements, regexpr("^[A-Za-z0-9_.-]+==[^ \\\\]+", requirements))
pins <- stats::setNames(sub("^.*==", "", pins), tolower(sub("==.*$", "", pins)))
mofa <- asNamespace("MOFA2")
dependency_pins <- strsplit(as.character(get0(".mofapy2_dependencies", envir = mofa, inherits = FALSE)), "=+")
mofa_pins <- c(mofapy2 = as.character(get(".mofapy2_version", envir = mofa)),
    stats::setNames(vapply(dependency_pins, function(x) x[[2L]], character(1)),
        tolower(vapply(dependency_pins, function(x) x[[1L]], character(1)))))
mofa_pins <- mofa_pins[names(mofa_pins) != "python"]
disagree <- names(mofa_pins)[is.na(pins[names(mofa_pins)]) | pins[names(mofa_pins)] != mofa_pins]
cat(sprintf("MOFA2 %s pins: %s\n", utils::packageVersion("MOFA2", lib.loc = lib),
    paste(names(mofa_pins), mofa_pins, sep = "==", collapse = ", ")))
if (length(disagree)) stop("requirements.txt disagrees with MOFA2's pins: ", paste(disagree, collapse = ", "))
cat("MOFA2 library ready\n")
