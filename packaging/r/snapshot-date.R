# KIT-134 P2 (revision 2.0, section 2.0.3): choose the Posit Package Manager snapshot date for the
# bundle. Build time only; run once, on a runner.
#
# Usage: Rscript --vanilla snapshot-date.R <lockfile> <first date> <last date> <output directory>
#
# Every calendar date from <first date> to <last date> (Posit Package Manager serves the latest
# snapshot on or before any date) is scored against the lockfile's CRAN and Bioconductor packages,
# using that date's source indexes (CRAN, Bioconductor software and annotation). The rule, in order:
#   1. duckdb is 1.5.5 (pinned by the package's code; a date without it is not eligible);
#   2. the most hard-closure packages at their locked version;
#   3. the most packages at their locked version overall;
#   4. the latest date on a tie.
# For the chosen date it then asks the snapshot for its Windows binary indexes (CRAN, and
# Bioconductor software as KIT-130 asked for jammy) and, for every lockfile package, whether a
# Windows binary exists and at which version.
#
# Writes into <output directory>: snapshot-scores.tsv (one row per date), snapshot-packages.tsv
# (the chosen date, one row per lockfile package), snapshot-windows-probe.tsv and
# snapshot-choice.txt (key: value lines). Needs jsonlite on the library path (a build-tool library).

PPM <- "https://packagemanager.posit.co"
BIOC_VERSION <- "3.23"
R_MINOR <- "4.6"
PINNED <- c(duckdb = "1.5.5")

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 4L) {
    stop("usage: Rscript --vanilla snapshot-date.R <lockfile> <first date> <last date> <output directory>")
}
lockfile <- jsonlite::fromJSON(arguments[[1L]], simplifyVector = FALSE)
dates <- seq(as.Date(arguments[[2L]]), as.Date(arguments[[3L]]), by = "day")
out_dir <- arguments[[4L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cache <- file.path(tempdir(), "indexes")
dir.create(cache, showWarnings = FALSE)
options(timeout = 300)

# The lockfile packages that come from a repository: CRAN and Bioconductor, installed on the laptop.
locked <- do.call(rbind, lapply(lockfile$Packages, function(p) {
    if (!isTRUE(p$Installed) || !p$Source %in% c("CRAN", "Bioconductor")) return(NULL)
    data.frame(package = p$Package, locked = p$Version, source = p$Source,
        hard = identical(p$Closure, "hard"), stringsAsFactors = FALSE)
}))
cat(sprintf("lockfile packages from a repository: %d (CRAN %d, Bioconductor %d; hard closure %d)\n",
    nrow(locked), sum(locked$source == "CRAN"), sum(locked$source == "Bioconductor"), sum(locked$hard)))

indexUrls <- function(date) {
    bioc <- sprintf("%s/bioconductor/%s/packages/%s", PPM, date, BIOC_VERSION)
    c(cran = sprintf("%s/cran/%s/src/contrib/PACKAGES.gz", PPM, date),
        bioc = sprintf("%s/bioc/src/contrib/PACKAGES.gz", bioc),
        annotation = sprintf("%s/data/annotation/src/contrib/PACKAGES.gz", bioc))
}

# Download many index files at once (libcurl), then retry any that failed one at a time.
fetch <- function(urls, files) {
    todo <- !file.exists(files)
    for (batch in split(which(todo), ceiling(seq_along(which(todo)) / 12))) {
        try(suppressWarnings(utils::download.file(urls[batch], files[batch], method = "libcurl",
            quiet = TRUE, mode = "wb")), silent = TRUE)
    }
    for (i in which(!file.exists(files) | file.size(files) == 0)) {
        unlink(files[[i]])
        try(suppressWarnings(utils::download.file(urls[[i]], files[[i]], method = "libcurl",
            quiet = TRUE, mode = "wb")), silent = TRUE)
    }
    file.exists(files) & file.size(files) > 0
}

# Package -> version from one index (the highest version when a package is listed twice).
readIndex <- function(file, gz = TRUE) {
    if (!file.exists(file)) return(NULL)
    connection <- if (gz) gzfile(file) else file(file)
    db <- tryCatch(read.dcf(connection, fields = c("Package", "Version")), error = function(e) NULL)
    if (is.null(db) || !nrow(db)) return(NULL)
    db <- db[order(db[, "Package"], package_version(db[, "Version"]), decreasing = TRUE), , drop = FALSE]
    db <- db[!duplicated(db[, "Package"]), , drop = FALSE]
    stats::setNames(db[, "Version"], db[, "Package"])
}

urls <- unlist(lapply(as.character(dates), indexUrls))
files <- file.path(cache, sprintf("%s-%s.gz", rep(as.character(dates), each = 3L), names(urls)))
cat(sprintf("fetching %d index files for %d dates\n", length(urls), length(dates)))
started <- Sys.time()
fetched <- fetch(urls, files)
cat(sprintf("fetched %d of %d in %.0f s\n", sum(fetched), length(urls),
    as.numeric(difftime(Sys.time(), started, units = "secs"))))

snapshotVersions <- function(date) {
    parts <- file.path(cache, sprintf("%s-%s.gz", date, c("cran", "bioc", "annotation")))
    cran <- readIndex(parts[[1L]])
    bioc <- c(readIndex(parts[[2L]]), readIndex(parts[[3L]]))
    if (is.null(cran) || is.null(bioc)) return(NULL)
    version <- ifelse(locked$source == "CRAN", cran[locked$package], bioc[locked$package])
    list(version = unname(version), cran_md5 = unname(tools::md5sum(parts[[1L]])),
        bioc_md5 = unname(tools::md5sum(parts[[2L]])))
}

scores <- do.call(rbind, lapply(as.character(dates), function(date) {
    snapshot <- snapshotVersions(date)
    if (is.null(snapshot)) {
        return(data.frame(date = date, indexes = "missing", duckdb = NA_character_, eligible = FALSE,
            hard_same = NA_integer_, all_same = NA_integer_, present = NA_integer_,
            cran_index_md5 = NA_character_, bioc_index_md5 = NA_character_))
    }
    same <- !is.na(snapshot$version) & snapshot$version == locked$locked
    duckdb <- snapshot$version[locked$package == "duckdb"]
    data.frame(date = date, indexes = "ok", duckdb = duckdb,
        eligible = identical(duckdb, PINNED[["duckdb"]]),
        hard_same = sum(same & locked$hard), all_same = sum(same),
        present = sum(!is.na(snapshot$version)),
        cran_index_md5 = snapshot$cran_md5, bioc_index_md5 = snapshot$bioc_md5)
}))
utils::write.table(scores, file.path(out_dir, "snapshot-scores.tsv"), sep = "\t", quote = FALSE,
    row.names = FALSE, na = "")

eligible <- scores[scores$eligible, , drop = FALSE]
cat(sprintf("\ndates scored: %d; with indexes: %d; eligible (duckdb %s): %d\n", nrow(scores),
    sum(scores$indexes == "ok"), PINNED[["duckdb"]], nrow(eligible)))
if (!nrow(eligible)) {
    stop("no date has duckdb ", PINNED[["duckdb"]], ": build duckdb from CRAN's archive source (revision 2.0, 2.0.3)")
}
cat(sprintf("eligible window: %s to %s\n", min(eligible$date), max(eligible$date)))
ranked <- eligible[order(eligible$hard_same, eligible$all_same, eligible$date, decreasing = TRUE,
    method = "radix"), , drop = FALSE]
chosen <- ranked$date[[1L]]
cat("\nbest eligible dates (hard closure, then overall, then latest):\n")
print(utils::head(ranked[, c("date", "duckdb", "hard_same", "all_same", "present")], 15L), row.names = FALSE)
cat(sprintf("\nCHOSEN SNAPSHOT DATE: %s\n", chosen))

# The chosen date, one row per lockfile package.
snapshot <- snapshotVersions(chosen)
packages <- data.frame(locked, snapshot = snapshot$version,
    same = !is.na(snapshot$version) & snapshot$version == locked$locked)
packages$relation <- ifelse(is.na(packages$snapshot), "absent",
    ifelse(packages$same, "same", ifelse(package_version(ifelse(is.na(packages$snapshot), "0", packages$snapshot)) >
        package_version(packages$locked), "newer", "older")))
utils::write.table(packages, file.path(out_dir, "snapshot-packages.tsv"), sep = "\t", quote = FALSE,
    row.names = FALSE, na = "")
cat("\nat the chosen date, by relation to the lock (hard closure / all):\n")
print(table(relation = packages$relation, hard = packages$hard))
cat("\nhard-closure packages not at their locked version:\n")
print(packages[packages$hard & !packages$same, c("package", "source", "locked", "snapshot")], row.names = FALSE)

# Windows binary indexes at the chosen date.
windows <- c(cran = sprintf("%s/cran/%s/bin/windows/contrib/%s/PACKAGES", PPM, chosen, R_MINOR),
    bioc = sprintf("%s/bioconductor/%s/packages/%s/bioc/bin/windows/contrib/%s/PACKAGES", PPM, chosen,
        BIOC_VERSION, R_MINOR))
windows_files <- file.path(cache, paste0("windows-", names(windows)))
windows_fetched <- fetch(windows, windows_files)
cran_win <- readIndex(windows_files[[1L]], gz = FALSE)
bioc_win <- readIndex(windows_files[[2L]], gz = FALSE)
probe <- data.frame(packages[, c("package", "source", "hard", "locked", "snapshot")],
    windows_binary = unname(ifelse(packages$source == "CRAN", cran_win[packages$package], bioc_win[packages$package])))
probe$windows_same_as_snapshot <- !is.na(probe$windows_binary) & probe$windows_binary == probe$snapshot
utils::write.table(probe, file.path(out_dir, "snapshot-windows-probe.tsv"), sep = "\t", quote = FALSE,
    row.names = FALSE, na = "")
# One Bioconductor binary actually served: the index alone does not show that the file is there.
bioc_file_ok <- FALSE
if (!is.null(bioc_win) && "limma" %in% names(bioc_win)) {
    target <- file.path(cache, "limma-windows.zip")
    url <- sprintf("%s/bioconductor/%s/packages/%s/bioc/bin/windows/contrib/%s/limma_%s.zip", PPM, chosen,
        BIOC_VERSION, R_MINOR, bioc_win[["limma"]])
    bioc_file_ok <- isTRUE(fetch(url, target)) && identical(readBin(target, "raw", 2L), charToRaw("PK"))
}
bioc_rows <- probe$source == "Bioconductor"
bioc_windows <- if (!is.null(bioc_win) && bioc_file_ok) "snapshot" else "bioconductor.org"
cat(sprintf(paste0("\nWindows binaries at %s: CRAN index %s (%d of %d CRAN packages, %d at the snapshot's version);",
    " Bioconductor index %s (%d of %d, %d at the snapshot's version); limma zip served: %s\n"),
    chosen, if (is.null(cran_win)) "missing" else "present",
    sum(!is.na(probe$windows_binary[!bioc_rows])), sum(!bioc_rows), sum(probe$windows_same_as_snapshot[!bioc_rows]),
    if (is.null(bioc_win)) "missing" else "present",
    sum(!is.na(probe$windows_binary[bioc_rows])), sum(bioc_rows), sum(probe$windows_same_as_snapshot[bioc_rows]),
    bioc_file_ok))
cat("packages without a Windows binary at the snapshot's version:\n")
print(probe[!probe$windows_same_as_snapshot, c("package", "source", "hard", "snapshot", "windows_binary")],
    row.names = FALSE)
cat(sprintf("WINDOWS BIOCONDUCTOR SOURCE: %s\n", bioc_windows))

writeLines(c(
    sprintf("date: %s", chosen),
    sprintf("rule: duckdb %s; most hard-closure packages at the lock; most packages at the lock; latest", PINNED[["duckdb"]]),
    sprintf("scored: %s to %s (%d dates, %d with indexes)", min(scores$date), max(scores$date), nrow(scores),
        sum(scores$indexes == "ok")),
    sprintf("eligible window: %s to %s (%d dates)", min(eligible$date), max(eligible$date), nrow(eligible)),
    sprintf("hard closure at the lock: %d of %d", ranked$hard_same[[1L]], sum(locked$hard)),
    sprintf("all at the lock: %d of %d", ranked$all_same[[1L]], nrow(locked)),
    sprintf("windows bioconductor: %s", bioc_windows)
), file.path(out_dir, "snapshot-choice.txt"))
