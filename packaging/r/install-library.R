# KIT-134 P2 (revision 2.0, section 2.0.3): install the bundle's R library into the library of the R
# that runs this script. Build time only; nothing here runs in the shipped application.
#
# Usage:
#   Rscript --vanilla install-library.R plan    <checkout> <output directory>
#   Rscript --vanilla install-library.R install <checkout> <output directory>
#
# <checkout> is the MultiScholaR source tree. Its DESCRIPTION is the root of the closure; the
# lockfile (packaging/multischolar.lock.json, its Snapshot block) and the exclusion list
# (packaging/r/library-exclude.txt) are read from it. jsonlite must be on the library path (a
# build-tool library, not the bundle's).
#
# plan    resolves the closure against the snapshot's source indexes and writes plan-packages.tsv
#         (every package to install and its repository, for the system requirements query);
#         installs nothing.
# install installs, in this order: (1) the snapshot's packages, as binaries where the snapshot
#         serves them; (2) the three source builds (GO.db from the snapshot, RUVIIIC and the Glimma
#         fork from GitHub by commit) and nanoarrow and limpa at the lockfile's versions (U4);
#         (2b) with MSR_INTERNAL set, the two internal packages (fastduckdb, fastlimpadpc) at
#         DESCRIPTION's Remotes commits, with the read token in GITHUB_PAT; (3) MultiScholaR from the
#         checkout. Writes bundle-sources.tsv (where each package came from) and, with MSR_INTERNAL,
#         internal-packages.tsv (each internal package: installed, or the failure).
#
# The closure: DESCRIPTION's Depends, Imports, LinkingTo and Suggests, less the internal packages,
# less the Suggests that were never installed on the laptop (the lockfile's NotInstalled), plus MOFA2
# (D7), less the exclusion list; then recursively Depends, Imports and LinkingTo.

PPM <- "https://packagemanager.posit.co"
BIOC_VERSION <- "3.23"
INTERNAL <- c("fastduckdb", "fastlimpadpc")
ADDED <- "MOFA2" # D7: carried, although it was never installed on the laptop
FROM_GITHUB <- list(
    RUVIIIC = list(user = "cran", repo = "RUVIIIC", sha = "c534151f6b2265b36294ff630f4dd13c2e86a9fa"),
    Glimma = list(user = "APAF-bioinformatics", repo = "GlimmaV2", sha = "743af0944957ce2515b17b38d329506edd06ca5a")
)
FROM_SOURCE <- "GO.db" # annotation data: source only; nothing to compile
# U4 (KIT-134 status 2026-11-17): the snapshot has nanoarrow 0.8.0-1 (DESCRIPTION declares >= 0.9.0)
# and limpa 1.4.0 (the lockfile has 1.4.2); both are built from source at the lockfile's version, on
# both platforms. The current file first, then the repository's archive.
AT_LOCK <- list(
    nanoarrow = list(version = "0.9.0", urls = c("https://cran.r-project.org/src/contrib/nanoarrow_0.9.0.tar.gz",
        "https://cran.r-project.org/src/contrib/Archive/nanoarrow/nanoarrow_0.9.0.tar.gz")),
    limpa = list(version = "1.4.2", urls = c("https://bioconductor.org/packages/3.23/bioc/src/contrib/limpa_1.4.2.tar.gz",
        "https://bioconductor.org/packages/3.23/bioc/src/contrib/Archive/limpa/limpa_1.4.2.tar.gz"))
)
# The internal packages, named explicitly (a GitHub reference names its repository, not the package
# inside it) at the commits DESCRIPTION's Remotes pin on v0.9-staging (formerly performance-v2) since 3ec13e95: fastduckdb 0.2.0
# at fast_duckdb e724f90b (the package code of 2694b972; Cargo goes offline only with vendored crates)
# and fastlimpadpc at fast_limpa_dpc 0ad8c527 (branch portable-x86-64-v3: 87a90b89, which declares Imports limpa and
# methods, built for x86-64-v3 instead of the building CPU, KIT-157). The checkout's
# Remotes must name the same commits, or the build stops. Installed only when MSR_INTERNAL is
# "required" (a failure stops the build) or "attempt" (a failure is recorded and the build goes on:
# Windows, D4), with the read token in GITHUB_PAT. While the repositories are internal their sources,
# and their R CMD INSTALL output (compiler messages quote source lines), stay under the runner's
# temporary folder: never in the output directory, never printed; only the INTERNAL result line is.
INTERNAL_REMOTES <- list(
    fastduckdb = list(user = "APAF-bioinformatics", repo = "fast_duckdb", subdir = "fastduckdb",
        sha = "e724f90bdbdf1073635b7112c80f8dd406f2c79b"),
    fastlimpadpc = list(user = "APAF-bioinformatics", repo = "fast_limpa_dpc", subdir = "",
        sha = "0ad8c5271da97605af3a1c8048d7b10fdfc0745b")
)
internal_mode <- Sys.getenv("MSR_INTERNAL", "")
if (!internal_mode %in% c("", "required", "attempt")) stop("MSR_INTERNAL must be empty, required or attempt")
# The token is read once and removed from the environment, so no configure script, Makevars or crate
# build script this process starts can see it; it is used only in the internal downloads' headers.
internal_token <- Sys.getenv("GITHUB_PAT")
Sys.unsetenv("GITHUB_PAT")

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 3L || !arguments[[1L]] %in% c("plan", "install")) {
    stop("usage: Rscript --vanilla install-library.R <plan|install> <checkout> <output directory>")
}
mode <- arguments[[1L]]
checkout <- normalizePath(arguments[[2L]], winslash = "/", mustWork = TRUE)
out_dir <- arguments[[3L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
is_windows <- .Platform$OS.type == "windows"
lib <- .Library
options(timeout = 1800, Ncpus = 4L, warn = 1)

lockfile <- jsonlite::fromJSON(file.path(checkout, "packaging", "multischolar.lock.json"), simplifyVector = FALSE)
snapshot <- lockfile$Snapshot$Date
if (is.null(snapshot) || !grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", snapshot)) {
    stop("the lockfile has no Snapshot$Date; run packaging/r/snapshot-date.R first")
}
windows_bioc <- lockfile$Snapshot$BioconductorWindows
excluded <- readLines(file.path(checkout, "packaging", "r", "library-exclude.txt"), warn = FALSE)
excluded <- trimws(excluded[!grepl("^\\s*(#|$)", excluded)])
cat(sprintf("R %s at %s; library %s\nsnapshot %s; mode %s\n", getRversion(), R.home(), lib, snapshot, mode))

# The snapshot's source indexes: the same versions on both platforms, used to resolve the closure.
bioc_dated <- sprintf("%s/bioconductor/%s/packages/%s", PPM, snapshot, BIOC_VERSION)
source_repos <- c(CRAN = sprintf("%s/cran/%s", PPM, snapshot),
    BioCsoft = paste0(bioc_dated, "/bioc"), BioCann = paste0(bioc_dated, "/data/annotation"))
# The repositories installed from.
if (is_windows) {
    bioc_root <- if (identical(windows_bioc, "snapshot")) bioc_dated else
        sprintf("https://bioconductor.org/packages/%s", BIOC_VERSION)
    install_repos <- c(CRAN = sprintf("%s/cran/%s", PPM, snapshot),
        BioCsoft = paste0(bioc_root, "/bioc"), BioCann = paste0(bioc_root, "/data/annotation"))
} else {
    # Posit Package Manager serves jammy binaries when R names itself in the user agent.
    bioc_root <- sprintf("%s/bioconductor/__linux__/jammy/%s/packages/%s", PPM, snapshot, BIOC_VERSION)
    install_repos <- c(CRAN = sprintf("%s/cran/__linux__/jammy/%s", PPM, snapshot),
        BioCsoft = paste0(bioc_root, "/bioc"), BioCann = paste0(bioc_root, "/data/annotation"))
    options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(),
        paste(getRversion(), R.version[["platform"]], R.version[["arch"]], R.version[["os"]])))
}
cat("installing from:\n")
print(install_repos)

# Package names from DESCRIPTION dependency fields, without version requirements.
dependencyNames <- function(fields) {
    fields <- fields[!is.na(fields)]
    if (!length(fields)) return(character())
    names <- trimws(sub("\\(.*$", "", unlist(strsplit(gsub("\n", " ", fields), ","))))
    setdiff(names[nzchar(names)], "R")
}

# The two GitHub packages, unpacked into a folder named after the package, with their provenance
# written into DESCRIPTION so the installed copy records it (RemoteType, RemoteSha, ...).
github_dir <- file.path(out_dir, "github-sources")
dir.create(github_dir, showWarnings = FALSE)
github <- lapply(names(FROM_GITHUB), function(package) {
    spec <- FROM_GITHUB[[package]]
    tarball <- file.path(github_dir, paste0(package, "-", spec$sha, ".tar.gz"))
    if (!file.exists(tarball)) {
        utils::download.file(sprintf("https://github.com/%s/%s/archive/%s.tar.gz", spec$user, spec$repo, spec$sha),
            tarball, mode = "wb", quiet = TRUE)
    }
    unpacked <- file.path(github_dir, "unpacked", package)
    unlink(unpacked, recursive = TRUE)
    staging <- file.path(github_dir, "staging")
    unlink(staging, recursive = TRUE)
    utils::untar(tarball, exdir = staging)
    top <- list.dirs(staging, recursive = FALSE)
    stopifnot(length(top) == 1L)
    dir.create(dirname(unpacked), showWarnings = FALSE)
    stopifnot(file.rename(top, unpacked))
    description <- read.dcf(file.path(unpacked, "DESCRIPTION"))
    if (!identical(unname(description[1L, "Package"]), package)) {
        stop(sprintf("%s/%s@%s holds package %s, not %s", spec$user, spec$repo, spec$sha,
            description[1L, "Package"], package))
    }
    provenance <- c(RemoteType = "github", RemoteHost = "api.github.com", RemoteUsername = spec$user,
        RemoteRepo = spec$repo, RemoteSha = spec$sha, RemoteRef = spec$sha)
    keep <- setdiff(colnames(description), names(provenance))
    description <- cbind(description[, keep, drop = FALSE], t(provenance))
    write.dcf(description, file.path(unpacked, "DESCRIPTION"))
    list(package = package, dir = unpacked, description = description, spec = spec)
})
names(github) <- names(FROM_GITHUB)

# The resolution database: the snapshot's source indexes, with the GitHub packages' own
# DESCRIPTIONs (the Glimma fork replaces Bioconductor's Glimma; RUVIIIC is not in any index).
db <- utils::available.packages(repos = source_repos, type = "source")
for (entry in github) {
    row <- matrix(NA_character_, 1L, ncol(db), dimnames = list(entry$package, colnames(db)))
    shared <- intersect(colnames(db), colnames(entry$description))
    row[1L, shared] <- entry$description[1L, shared]
    row[1L, "Repository"] <- sprintf("https://github.com/%s/%s", entry$spec$user, entry$spec$repo)
    db <- rbind(db[rownames(db) != entry$package, , drop = FALSE], row)
}

description <- read.dcf(file.path(checkout, "DESCRIPTION"))
# Remotes in either form: "user/repo[/subdir]@sha" or "package=github::user/repo[/subdir]@sha". When the
# internal packages are installed, the repository, subfolder and commit must be the ones DESCRIPTION
# names; any difference stops the build.
remotes <- trimws(unlist(strsplit(gsub("\n", " ", description[1L, "Remotes"]), ",")))
remote_paths <- sub("@.*$", "", sub("^[A-Za-z0-9.]+=(github::)?", "", remotes))
remote_shas <- ifelse(grepl("@", remotes), sub("^.*@", "", remotes), "")
for (package in if (nzchar(internal_mode)) names(INTERNAL_REMOTES) else character()) {
    spec <- INTERNAL_REMOTES[[package]]
    path <- sprintf("%s/%s%s", spec$user, spec$repo, if (nzchar(spec$subdir)) paste0("/", spec$subdir) else "")
    if (!path %in% remote_paths) stop(sprintf("%s: %s is not in DESCRIPTION's Remotes", package, path))
    named <- remote_shas[match(path, remote_paths)]
    if (!identical(named, spec$sha)) {
        stop(sprintf("%s: DESCRIPTION's Remotes pin %s@%s, this script %s; they must agree", package, path, named, spec$sha))
    }
}
declared <-dependencyNames(description[1L, intersect(c("Depends", "Imports", "LinkingTo", "Suggests"),
    colnames(description))])
# An internal package that DESCRIPTION declares outside Suggests (fastduckdb: Imports since 2026-11-18, the
# user's decision) must be installed before MultiScholaR. A build without the read token stops here,
# before anything is installed, instead of at MultiScholaR's own install after every other package.
internal_required <- intersect(INTERNAL, dependencyNames(description[1L, intersect(c("Depends", "Imports", "LinkingTo"),
    colnames(description))]))
if (identical(mode, "install") && length(internal_required) && !nzchar(internal_mode)) {
    stop(sprintf("DESCRIPTION requires %s (an internal repository): this build needs MSR_INTERNAL and the read token",
        paste(internal_required, collapse = ", ")))
}
never_installed <- setdiff(unlist(lockfile$Lockfile$NotInstalled), ADDED)
roots <- setdiff(unique(c(declared, ADDED)), c(INTERNAL, never_installed, excluded))
base_packages <- rownames(utils::installed.packages(lib.loc = .Library, priority = "base"))
closure <- unique(c(roots, unlist(tools::package_dependencies(roots, db = db,
    which = c("Depends", "Imports", "LinkingTo"), recursive = TRUE), use.names = FALSE)))
closure <- sort(setdiff(closure, base_packages))

needs_excluded <- intersect(closure, excluded)
if (length(needs_excluded)) {
    direct <- tools::package_dependencies(setdiff(closure, excluded), db = db,
        which = c("Depends", "Imports", "LinkingTo"))
    for (tool in needs_excluded) {
        cat(sprintf("%s is needed by: %s\n", tool,
            paste(names(Filter(function(d) tool %in% d, direct)), collapse = ", ")))
    }
    stop("shipped packages need excluded tools: ", paste(needs_excluded, collapse = ", "))
}
unresolved <- setdiff(closure, rownames(db))
if (length(unresolved)) stop("not in the snapshot or the GitHub sources: ", paste(unresolved, collapse = ", "))

locked <- Filter(function(p) isTRUE(p$Installed), lockfile$Packages)
locked_hard <- vapply(Filter(function(p) identical(p$Closure, "hard") && !identical(p$Source, "base R"), locked),
    function(p) p$Package, character(1))
missing_hard <- setdiff(locked_hard, closure)
if (length(missing_hard)) stop("hard-closure packages of the lockfile not in the plan: ", paste(missing_hard, collapse = ", "))
locked_names <- vapply(locked, function(p) p$Package, character(1))
cat(sprintf(paste0("closure: %d packages (roots %d); excluded tools %d; the lockfile's installed packages",
    " not in the plan: %d; planned packages not in the lockfile: %d\n"),
    length(closure), length(roots), length(excluded),
    length(setdiff(locked_names, c(closure, base_packages, INTERNAL))), length(setdiff(closure, locked_names))))
cat("lockfile packages left out:", paste(sort(setdiff(locked_names, c(closure, base_packages, INTERNAL))), collapse = ", "), "\n")
cat("planned, not in the lockfile:", paste(sort(setdiff(closure, locked_names)), collapse = ", "), "\n")
repository <- ifelse(closure %in% names(FROM_GITHUB), "GitHub",
    ifelse(startsWith(db[closure, "Repository"], source_repos[["CRAN"]]), "CRAN", "Bioconductor"))
utils::write.table(data.frame(package = closure, repository = repository), file.path(out_dir, "plan-packages.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)
if (mode == "plan") quit(save = "no", status = 0L)

# (1) The snapshot's packages.
from_repositories <- setdiff(closure, c(names(FROM_GITHUB), FROM_SOURCE))
requested <- stats::setNames(rep("source", length(from_repositories)), from_repositories)
install <- function(packages, type, repos = install_repos) {
    if (!length(packages)) return(invisible())
    cat(sprintf("\n== installing %d packages (type %s)\n", length(packages), type))
    utils::install.packages(packages, lib = lib, repos = repos, type = type, dependencies = FALSE,
        INSTALL_opts = "--no-test-load")
}
if (is_windows) {
    binaries <- utils::available.packages(repos = install_repos, type = "win.binary")
    as_binary <- intersect(from_repositories, rownames(binaries))
    requested[as_binary] <- "binary"
    install(as_binary, "binary")
    install(setdiff(from_repositories, as_binary), "source")
} else {
    install(from_repositories, "source")
}

# (2) The source builds.
install(FROM_SOURCE, "source", repos = install_repos[["BioCann"]])
requested[FROM_SOURCE] <- "source"
r_command <- file.path(R.home("bin"), if (is_windows) "R.exe" else "R")
rInstall <- function(dir, log = "") {
    if (!nzchar(log)) cat(sprintf("\n== R CMD INSTALL %s\n", dir))
    quote <- function(x) shQuote(x, type = if (is_windows) "cmd" else "sh")
    status <- system2(r_command, c("CMD", "INSTALL", "--no-test-load", paste0("--library=", quote(lib)), quote(dir)),
        stdout = log, stderr = log)
    if (status != 0L) stop("R CMD INSTALL failed for ", dir)
}
for (entry in github) {
    rInstall(entry$dir)
    requested[entry$package] <- "source"
}
at_lock_origin <- character()
for (package in names(AT_LOCK)) {
    spec <- AT_LOCK[[package]]
    tarball <- file.path(tempdir(), sprintf("%s_%s.tar.gz", package, spec$version))
    for (url in spec$urls) {
        fetched <- tryCatch(utils::download.file(url, tarball, mode = "wb", quiet = TRUE) == 0L, error = function(e) FALSE,
            warning = function(w) FALSE)
        if (fetched) break
    }
    if (!fetched) stop(sprintf("%s %s: no source at %s", package, spec$version, paste(spec$urls, collapse = " or ")))
    cat(sprintf("%s %s from %s, sha256 %s\n", package, spec$version, url, digest::digest(file = tarball, algo = "sha256")))
    rInstall(tarball)
    if (!identical(as.character(utils::packageVersion(package, lib.loc = lib)), spec$version)) {
        stop(sprintf("%s: %s installed, %s expected", package, utils::packageVersion(package, lib.loc = lib), spec$version))
    }
    requested[package] <- "source"
    at_lock_origin[package] <- sprintf("source %s (lockfile version, U4)", url)
}

# (2b) The internal packages, after duckdb and before MultiScholaR, whose configure records fastduckdb.
internalSource <- function(package, spec) {
    if (!nzchar(internal_token)) stop("GITHUB_PAT is empty: the internal packages need the read token")
    # A short folder: Cargo's build-script paths under R's session folder and the tarball's
    # <owner>-<repo>-<sha> folder reached Windows' 260-character limit (run 44015367782, LNK1104).
    dir <- file.path(internal_root, package)
    unlink(dir, recursive = TRUE)
    dir.create(dir, recursive = TRUE)
    tarball <- file.path(dir, "source.tar.gz")
    utils::download.file(sprintf("https://api.github.com/repos/%s/%s/tarball/%s", spec$user, spec$repo, spec$sha),
        tarball, mode = "wb", quiet = TRUE, method = "libcurl",
        headers = c(Authorization = paste("Bearer", internal_token), Accept = "application/vnd.github+json"))
    utils::untar(tarball, exdir = file.path(dir, "u"))
    unlink(tarball)
    top <- list.dirs(file.path(dir, "u"), recursive = FALSE)
    stopifnot(length(top) == 1L, file.rename(top, file.path(dir, "s")))
    top <- file.path(dir, "s")
    root <- if (nzchar(spec$subdir)) file.path(top, spec$subdir) else top
    package_description <- read.dcf(file.path(root, "DESCRIPTION"))
    if (!identical(unname(package_description[1L, "Package"]), package)) {
        stop(sprintf("%s/%s@%s holds package %s, not %s", spec$user, spec$repo, spec$sha,
            package_description[1L, "Package"], package))
    }
    provenance <- c(RemoteType = "github", RemoteHost = "api.github.com", RemoteUsername = spec$user,
        RemoteRepo = spec$repo, RemoteSubdir = spec$subdir, RemoteSha = spec$sha, RemoteRef = spec$sha)
    provenance <- provenance[nzchar(provenance)]
    keep <- setdiff(colnames(package_description), names(provenance))
    write.dcf(cbind(package_description[, keep, drop = FALSE], t(provenance)), file.path(root, "DESCRIPTION"))
    root
}
internal_rows <- list()
# Outside the output directory (uploaded): the runner's temporary folder, or R's session folder.
internal_root <- file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "i")
internal_logs <- file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "internal-build-logs")
if (nzchar(internal_mode)) {
    dir.create(internal_logs, showWarnings = FALSE)
    for (package in names(INTERNAL_REMOTES)) {
        spec <- INTERNAL_REMOTES[[package]]
        started <- Sys.time()
        result <- tryCatch({
            rInstall(internalSource(package, spec), log = file.path(internal_logs, paste0(package, ".log")))
            "installed"
        }, error = function(e) paste("FAILED:", gsub("[\t\n]+", " ", conditionMessage(e))))
        version <- if (identical(result, "installed")) as.character(utils::packageVersion(package, lib.loc = lib)) else ""
        cat(sprintf("INTERNAL %s %s@%s: %s %s (%.0f s)\n", package, spec$repo, substr(spec$sha, 1L, 8L), result, version,
            as.numeric(difftime(Sys.time(), started, units = "secs"))))
        internal_rows[[package]] <- data.frame(package = package, commit = spec$sha, result = result, version = version,
            seconds = round(as.numeric(difftime(Sys.time(), started, units = "secs"))), stringsAsFactors = FALSE)
        if (identical(result, "installed")) requested[package] <- "source"
        else if (identical(internal_mode, "required")) stop(sprintf("internal package %s: %s", package, result))
    }
    utils::write.table(do.call(rbind, internal_rows), file.path(out_dir, "internal-packages.tsv"), sep = "\t",
        quote = FALSE, row.names = FALSE)
    unlink(internal_root, recursive = TRUE)
}

# Linux bundle build (MSR_BUNDLE_RUNPATH=1, packaging/linux/build.sh): the run path to the bundle's
# syslib/ goes into duckdb's, duckplyr's and fastduckdb's libraries now. MultiScholaR's install
# records their sha256 (src/build_dia_rsum.R), and the DIA rollup and the direct writer refuse a
# build whose binaries changed afterwards; assemble.sh writes the run path of everything else.
if (!is_windows && identical(Sys.getenv("MSR_BUNDLE_RUNPATH"), "1")) {
    hashed <- file.path(lib, intersect(c("duckdb", "duckplyr", "fastduckdb"), rownames(utils::installed.packages(lib.loc = lib))))
    if (any(c("duckdb", "duckplyr", "fastduckdb") %in% loadedNamespaces())) stop("a package whose libraries are about to be rewritten is loaded")
    status <- system2("bash", shQuote(c(file.path(checkout, "packaging", "linux", "runpath.sh"),
        file.path(dirname(R.home()), "syslib"), hashed)))
    if (status != 0L) stop("the run path could not be written into ", paste(basename(hashed), collapse = ", "))
}

# (3) MultiScholaR from the checkout.
rInstall(checkout)

# KIT-176 (D4, VAL-31): the bundle lane's package commit (MSR_PACKAGE_COMMIT, packaging/package-commit.sh) is written
# into the installed package as Config/MultiScholaR/SourceCommit, in DESCRIPTION and Meta/package.rds (which
# packageDescription and so workflowSoftwareIdentity() read), as tools/profiling/run_workflow_kit.sh does, and read
# back. Unset (the run-tests lane's snapshot-library action), nothing is written.
package_commit <- Sys.getenv("MSR_PACKAGE_COMMIT")
if (nzchar(package_commit)) {
    if (!grepl("^[0-9a-f]{40}$", package_commit)) {
        stop(sprintf("MSR_PACKAGE_COMMIT must be the 40-character package commit, not '%s'", package_commit))
    }
    installed_package <- file.path(lib, "MultiScholaR")
    meta <- file.path(installed_package, "Meta", "package.rds")
    package_meta <- readRDS(meta)
    package_meta$DESCRIPTION["Config/MultiScholaR/SourceCommit"] <- package_commit
    saveRDS(package_meta, meta)
    cat("Config/MultiScholaR/SourceCommit: ", package_commit, "\n", sep = "", file = file.path(installed_package, "DESCRIPTION"),
        append = TRUE)
    written <- c(DESCRIPTION = unname(read.dcf(file.path(installed_package, "DESCRIPTION"),
        fields = "Config/MultiScholaR/SourceCommit")[1L, 1L]),
        package.rds = unname(utils::packageDescription("MultiScholaR", lib.loc = lib,
            fields = "Config/MultiScholaR/SourceCommit")))
    if (!all(written %in% package_commit)) {
        stop(sprintf("Config/MultiScholaR/SourceCommit reads back as %s, not %s",
            paste(names(written), written, sep = " ", collapse = ", "), package_commit))
    }
    cat(sprintf("MultiScholaR: Config/MultiScholaR/SourceCommit %s written into the installed package\n", package_commit))
}

installed <- utils::installed.packages(lib.loc = lib, fields = c("RemoteSha", "RemoteUsername", "RemoteRepo"))
missing <- setdiff(c(closure, "MultiScholaR"), rownames(installed))
if (length(missing)) stop("not installed: ", paste(missing, collapse = ", "))

# MOFA2's run_mofa() checks the mofapy2 it finds against the pin it carries; the bundled Python is
# installed from packaging/python/requirements.txt, so the two must agree.
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
cat(sprintf("MOFA2 %s pins: %s\n", utils::packageVersion("MOFA2"),
    paste(names(mofa_pins), mofa_pins, sep = "==", collapse = ", ")))
if (length(disagree)) stop("requirements.txt disagrees with MOFA2's pins: ", paste(disagree, collapse = ", "))

commit <- Sys.getenv("GITHUB_SHA", unset = NA_character_)
origin <- vapply(names(requested), function(package) {
    if (package %in% names(FROM_GITHUB)) {
        spec <- FROM_GITHUB[[package]]
        return(sprintf("GitHub %s/%s@%s", spec$user, spec$repo, spec$sha))
    }
    if (package %in% names(INTERNAL_REMOTES)) {
        spec <- INTERNAL_REMOTES[[package]]
        return(sprintf("GitHub %s/%s%s@%s (internal, read token)", spec$user, spec$repo,
            if (nzchar(spec$subdir)) paste0("/", spec$subdir) else "", spec$sha))
    }
    if (package %in% names(at_lock_origin)) return(at_lock_origin[[package]])
    repository <- db[package, "Repository"]
    if (startsWith(repository, source_repos[["CRAN"]])) "CRAN" else if (is_windows && !identical(windows_bioc, "snapshot"))
        "bioconductor.org" else "Bioconductor"
}, character(1))
sources <- data.frame(package = c(names(requested), "MultiScholaR"),
    version = unname(installed[c(names(requested), "MultiScholaR"), "Version"]),
    origin = c(unname(origin), sprintf("checkout %s", commit)),
    requested = c(unname(requested), "source"),
    snapshot = snapshot, stringsAsFactors = FALSE)
utils::write.table(sources, file.path(out_dir, "bundle-sources.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\ninstalled %d packages and MultiScholaR %s into %s\n", length(closure),
    installed["MultiScholaR", "Version"], lib))
