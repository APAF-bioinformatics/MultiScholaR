# KIT-157 (branch run-k157-rc only): fastlimpadpc at the previous pin into a separate library (see old-limpa.sh).
# Usage: Rscript --vanilla old-limpa.R <library>
# Prints one line: K157OLD fastlimpadpc 87a90b89: installed|FAILED ... (seconds). Exit status 0 either way (the comparison
# step reports a missing build); nothing else is printed. As packaging/r/install-library.R: the token from GITHUB_PAT,
# read once and unset; the download in a short folder under RUNNER_TEMP; R CMD INSTALL output in a log there.
SHA <- "87a90b894b2057744fc9c315a4fb7d7db3ce3af7"
lib <- commandArgs(trailingOnly = TRUE)[[1L]]
token <- Sys.getenv("GITHUB_PAT")
Sys.unsetenv("GITHUB_PAT")
root <- file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "k157-i")
log <- file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "k157-old-build.log")
started <- Sys.time()
result <- tryCatch({
    if (!nzchar(token)) stop("GITHUB_PAT is empty")
    unlink(root, recursive = TRUE)
    dir.create(root, recursive = TRUE)
    dir.create(lib, recursive = TRUE, showWarnings = FALSE)
    tarball <- file.path(root, "source.tar.gz")
    utils::download.file(sprintf("https://api.github.com/repos/APAF-bioinformatics/fast_limpa_dpc/tarball/%s", SHA),
        tarball, mode = "wb", quiet = TRUE, method = "libcurl",
        headers = c(Authorization = paste("Bearer", token), Accept = "application/vnd.github+json"))
    utils::untar(tarball, exdir = file.path(root, "u"))
    unlink(tarball)
    top <- list.dirs(file.path(root, "u"), recursive = FALSE)
    stopifnot(length(top) == 1L, file.rename(top, file.path(root, "s")))
    source <- file.path(root, "s")
    description <- read.dcf(file.path(source, "DESCRIPTION"))
    stopifnot(identical(unname(description[1L, "Package"]), "fastlimpadpc"))
    provenance <- c(RemoteType = "github", RemoteHost = "api.github.com", RemoteUsername = "APAF-bioinformatics",
        RemoteRepo = "fast_limpa_dpc", RemoteSha = SHA, RemoteRef = SHA)
    keep <- setdiff(colnames(description), names(provenance))
    write.dcf(cbind(description[, keep, drop = FALSE], t(provenance)), file.path(source, "DESCRIPTION"))
    status <- system2(file.path(R.home("bin"), "R"), c("CMD", "INSTALL", "--no-test-load", paste0("--library=", shQuote(lib)),
        shQuote(source)), stdout = log, stderr = log)
    if (status != 0L) stop("R CMD INSTALL failed (exit ", status, "; log under RUNNER_TEMP, not printed)")
    "installed"
}, error = function(e) paste("FAILED:", gsub("[\t\n]+", " ", conditionMessage(e))))
unlink(root, recursive = TRUE)
cat(sprintf("K157OLD fastlimpadpc %s: %s (%.0f s)\n", substr(SHA, 1L, 8L), result,
    as.numeric(difftime(Sys.time(), started, units = "secs"))))
