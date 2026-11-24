# KIT-167 step 1 (VAL-38): capture of the four UniProt reference proteomes shipped with the reference-mapping
# tool and of the UniParc proteomes of the study authors' three references, at the UniProt release current at
# capture (KIT-167 revision 1.9, section 4, "Pinning"). Runs on a GitHub-hosted runner only, as a run-tests
# `scripts` request; UniProt is never called from the owner's machine.
#
# Usage (from the package root): Rscript docs/qa/workflow-execution/kit/kit-167/capture-references.R
# Writes only into capture-out/ (relative to the working directory; no setwd()):
#   <UPID>.fasta            per reference, the FTP per-proteome FASTA, decompressed, bytes unmodified
#   <UPID>_xref.tsv         per reference, UniProtKB stream (cross-reference fields): header line verbatim,
#                           data lines sorted by the first field (C collation)
#   <UPID>_go.tsv           per reference, UniProtKB stream (the fields getUniprotAnnotations requests), as above
#   <UPID>_uniparc.fasta    per authors' reference, UniParc stream, records sorted by UPI (C collation)
#   hashes.tsv              the fifteen data files: file, bytes, sha256, entries (sorted by file, C collation)
#   release.tsv             url, source, release, release_date, used (sorted by url, C collation; no time stamps)
#   checks.tsv              check, reference, value, result; the job's start and end time (UTC)
#   citation.tsv            UniProt's paper as the help page lists it, `confirmed` or `not_confirmed`
#   citation-page.txt       the help page's text (not committed; read instead of the live page)
#   scan-precount.tsv       per file and scan pattern, the matching lines and the entries' accessions (never the
#                           pattern or the matched text)
# Exits 1 when a check fails. Every download: up to three tries, 30 s apart, 300 s timeout each; a download still
# failing is the check `download` (an instrumentation fault). Packages: base R, curl, digest, jsonlite (in the lock).

OUT_DIR <- "capture-out"
REFERENCES <- c("UP000000625", "UP000000265", "UP000008816", "UP000000750")
AUTHORS_REFERENCES <- c("UP000001734", "UP000000530", "UP000002436")
FTP_RELEASE <- "https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
FTP_PROTEOMES <- paste0(FTP_RELEASE, "knowledgebase/reference_proteomes/")
REST <- "https://rest.uniprot.org/"
XREF_FIELDS <- "accession,gene_primary,gene_oln,xref_string,xref_kegg"
GO_FIELDS <- paste0(
    "accession,id,protein_name,gene_names,organism_name,length,go_id,reviewed,",
    "protein_existence,annotation_score"
)
CITATION_URL <- paste0(REST, "help/publications")
CITATION_TITLE <- "UniProt: the Universal Protein Knowledgebase in ([0-9]{4})"
DOI_PATTERN <- "10\\.[0-9]+/[^[:space:]\"'<>()\\[\\]{}]+"
MAX_BYTES <- 5e6
TRIES <- 3L
PAUSE_SECONDS <- 30
TIMEOUT_SECONDS <- 300

# Verbatim copy of PATTERNS in .github/run-tests/scan.R at run-tests commit
# 8caa68dc5b31883f30738e1eea205f04eade6f06 (unchanged at run-tests 198d73c3). Used only to count matching lines;
# neither the patterns nor the matched text is written or printed.
SCAN_COMMIT <- "8caa68dc5b31883f30738e1eea205f04eade6f06"
PATTERNS <- c(
    rust_location = "\\.rs:[0-9]", rust_tree = "src/rust|\\.cargo|Cargo\\.(toml|lock)", rust_build = "rustc|cargo |panicked at",
    rust_messages = "Compiling [a-z_]+ v[0-9]|error\\[E[0-9]+\\]|warning: unused",
    c_compiler = "(^|[ /])(gcc|g\\+\\+|clang|cc1|ld)( |:)|-fPIC|-std=gnu", make = "make(\\[[0-9]+\\])?: \\*\\*\\*",
    c_location = "\\.(c|cc|cpp|h|hpp|o):[0-9]+", internal_repository = "fast_limpa_dpc|fast_duckdb|internal-build-logs",
    token = "gh[pousr]_[A-Za-z0-9]{20,}|github_pat_|Bearer "
)

# --- Records ---

newRun <- function() {
    run <- new.env(parent = emptyenv())
    run$checks <- list()
    run$releases <- list()
    run$citation <- NULL
    run$downloads <- 0L
    run$start <- utcNow()
    run
}

utcNow <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

oneLine <- function(x) {
    x <- as.character(x)
    x[is.na(x) | !nzchar(x)] <- NA_character_
    gsub("[[:space:]]+", " ", x)
}

recordCheck <- function(run, check, reference, value, pass) {
    result <- if (isTRUE(pass)) "pass" else "fail"
    run$checks[[length(run$checks) + 1L]] <- data.frame(
        check = check, reference = reference, value = oneLine(value), result = result
    )
    cat(sprintf("check\t%s\t%s\t%s\t%s\n", check, reference, oneLine(value), result))
    invisible(isTRUE(pass))
}

# A failed check that ends the job: recorded, then a classed condition stops the capture.
stopCapture <- function(run, check, reference, value) {
    recordCheck(run, check, reference, value, FALSE)
    stop(structure(
        class = c("capture_stop", "error", "condition"),
        list(message = sprintf("check %s failed (%s)", check, reference), call = NULL)
    ))
}

recordRelease <- function(run, url, source, release, release_date = NA, used = NA) {
    run$releases[[length(run$releases) + 1L]] <- data.frame(
        url = url, source = source, release = oneLine(release),
        release_date = oneLine(release_date), used = used
    )
}

writeTsv <- function(table, path) {
    table[] <- lapply(table, \(x) ifelse(is.na(x), "NA", as.character(x)))
    lines <- c(
        paste(names(table), collapse = "\t"),
        if (nrow(table)) do.call(paste, c(unname(as.list(table)), sep = "\t"))
    )
    con <- file(path, open = "wb")
    on.exit(close(con))
    writeLines(enc2utf8(lines), con, useBytes = TRUE)
}

# --- Downloads ---

# One URL to a file: up to TRIES attempts, PAUSE_SECONDS apart, each limited to TIMEOUT_SECONDS. A 404 is an
# answer, not a transient failure, only where the caller says so (`absent_ok`, RELEASE.metalink).
fetchUrl <- function(run, url, dest, absent_ok = FALSE) {
    outcome <- "not tried"
    for (attempt in seq_len(TRIES)) {
        if (attempt > 1L) Sys.sleep(PAUSE_SECONDS)
        handle <- curl::new_handle(timeout = TIMEOUT_SECONDS, connecttimeout = 60)
        reply <- tryCatch(curl::curl_fetch_disk(url, dest, handle = handle), error = \(e) e)
        if (inherits(reply, "error")) {
            outcome <- conditionMessage(reply)
            cat(sprintf("download try %d of %d failed: %s (%s)\n", attempt, TRIES, url, outcome))
            next
        }
        if (reply$status_code == 200L) {
            run$downloads <- run$downloads + 1L
            return(list(path = dest, headers = curl::parse_headers_list(reply$headers)))
        }
        if (reply$status_code == 404L && absent_ok) return(NULL)
        outcome <- sprintf("HTTP %d", reply$status_code)
        cat(sprintf("download try %d of %d failed: %s (%s)\n", attempt, TRIES, url, outcome))
    }
    stopCapture(run, "download", url, sprintf("%d tries; last: %s", TRIES, outcome))
}

fetchLines <- function(run, url, absent_ok = FALSE) {
    dest <- tempfile()
    reply <- fetchUrl(run, url, dest, absent_ok)
    if (is.null(reply)) return(NULL)
    reply$lines <- readLines(dest, warn = FALSE, encoding = "UTF-8")
    reply
}

releaseHeader <- function(headers) {
    value <- headers[["x-uniprot-release"]]
    if (is.null(value)) NA_character_ else trimws(value[[1L]])
}

listingEntries <- function(run, url) {
    html <- paste(fetchLines(run, url)$lines, collapse = "\n")
    hrefs <- regmatches(html, gregexpr("href=\"[^\"]+\"", html))[[1L]]
    hrefs <- sub("^href=\"", "", sub("\"$", "", hrefs))
    hrefs[!grepl("^([?/]|\\.\\.|[a-z]+:)", hrefs)]
}

# --- The FTP release ---

readFtpRelease <- function(run) {
    metalink_url <- paste0(FTP_RELEASE, "RELEASE.metalink")
    reldate_url <- paste0(FTP_RELEASE, "knowledgebase/complete/reldate.txt")
    metalink <- fetchLines(run, metalink_url, absent_ok = TRUE)
    text <- if (is.null(metalink)) "" else paste(metalink$lines, collapse = "\n")
    release <- regmatches(text, regexpr("<version>[[:space:]]*[0-9]{4}_[0-9]{2}", text))
    release <- if (length(release)) sub("^<version>[[:space:]]*", "", release) else NA_character_
    date <- regmatches(text, regexpr("<(pubdate|published|releasedate|date)>[^<]+<", text))
    date <- if (length(date)) sub("<$", "", sub("^<[a-z]+>", "", date)) else NA_character_
    found <- if (is.null(metalink)) "absent (HTTP 404)" else "present"
    recordRelease(run, metalink_url, paste("ftp RELEASE.metalink,", found), release, date, !is.na(release))
    if (!is.na(release)) return(recordCheckRelease(run, release, "RELEASE.metalink"))
    reldate <- fetchLines(run, reldate_url, absent_ok = TRUE)
    first <- if (is.null(reldate) || !length(reldate$lines)) "" else reldate$lines[[1L]]
    pattern <- "^UniProt Knowledgebase Release ([0-9]{4}_[0-9]{2}) consists of:"
    release <- if (grepl(pattern, first)) sub(paste0(pattern, ".*$"), "\\1", first) else NA_character_
    date <- if (is.null(reldate)) character() else
        regmatches(reldate$lines, regexpr("Release [0-9]{4}_[0-9]{2} of [^[:space:]]+", reldate$lines))
    date <- if (length(date)) sub("^.* of ", "", date[[1L]]) else NA_character_
    found <- if (is.null(reldate)) "absent (HTTP 404)" else "present"
    recordRelease(run, reldate_url, paste("ftp reldate.txt,", found), release, date, !is.na(release))
    if (is.na(release)) {
        stopCapture(run, "ftp_release", "all", "neither RELEASE.metalink nor reldate.txt names a release")
    }
    recordCheckRelease(run, release, "reldate.txt")
}

recordCheckRelease <- function(run, release, source) {
    recordCheck(run, "ftp_release", "all", sprintf("%s from %s", release, source), TRUE)
    release
}

# --- Reference FASTA (FTP) ---

# Each UPID's `<UPID>_<digits>.fasta.gz`, found by listing the kingdom folders and the proteome's own folder.
findReferenceFastas <- function(run) {
    folders <- grep("/$", listingEntries(run, FTP_PROTEOMES), value = TRUE)
    found <- setNames(vector("list", length(REFERENCES)), REFERENCES)
    for (folder in folders) {
        entries <- listingEntries(run, paste0(FTP_PROTEOMES, folder))
        for (upid in REFERENCES) {
            file_pattern <- sprintf("^%s_[0-9]+\\.fasta\\.gz$", upid)
            direct <- grep(file_pattern, entries, value = TRUE)
            if (length(direct)) found[[upid]] <- c(found[[upid]], paste0(FTP_PROTEOMES, folder, direct))
            if (paste0(upid, "/") %in% entries) {
                inner_url <- paste0(FTP_PROTEOMES, folder, upid, "/")
                inner <- grep(file_pattern, listingEntries(run, inner_url), value = TRUE)
                if (length(inner)) found[[upid]] <- c(found[[upid]], paste0(inner_url, inner))
            }
        }
    }
    found
}

readGzBytes <- function(path) {
    con <- gzfile(path, open = "rb")
    on.exit(close(con))
    chunks <- list()
    repeat {
        chunk <- readBin(con, "raw", n = 1e6)
        if (!length(chunk)) break
        chunks[[length(chunks) + 1L]] <- chunk
    }
    unlist(chunks, use.names = FALSE)
}

captureReferenceFasta <- function(run, upid, urls) {
    value <- if (length(urls)) paste(basename(urls), collapse = ";") else "none listed"
    if (!recordCheck(run, "fasta_file", upid, value, length(urls) == 1L)) return(NULL)
    archive <- tempfile(fileext = ".fasta.gz")
    fetchUrl(run, urls[[1L]], archive)
    target <- file.path(OUT_DIR, paste0(upid, ".fasta"))
    writeBin(readGzBytes(archive), target)
    unlink(archive)
    target
}

fastaAccessions <- function(path) {
    headers <- grep("^>", readLines(path, warn = FALSE), value = TRUE)
    vapply(strsplit(sub("^>", "", headers), "|", fixed = TRUE), \(x) x[[min(2L, length(x))]], "")
}

# --- UniProtKB tables and UniParc proteomes (REST) ---

captureTable <- function(run, upid, kind, fields) {
    url <- sprintf("%suniprotkb/stream?query=proteome:%s&format=tsv&fields=%s", REST, upid, fields)
    reply <- fetchLines(run, url)
    release <- releaseHeader(reply$headers)
    recordRelease(run, url, paste("uniprotkb", kind, upid), release)
    body <- reply$lines[-1L]
    body <- body[nzchar(body)]
    body <- body[order(sub("\t.*$", "", body), method = "radix")]
    target <- file.path(OUT_DIR, sprintf("%s_%s.tsv", upid, kind))
    writeLinesBytes(c(reply$lines[[1L]], body), target)
    list(path = target, release = release, accessions = sub("\t.*$", "", body))
}

captureUniparc <- function(run, upid) {
    url <- sprintf("%suniparc/stream?query=upid:%s&format=fasta", REST, upid)
    reply <- fetchLines(run, url)
    recordRelease(run, url, paste("uniparc", upid), releaseHeader(reply$headers))
    lines <- reply$lines[nzchar(reply$lines)]
    starts <- grep("^>", lines)
    records <- split(lines, cumsum(seq_along(lines) %in% starts))
    if (length(starts) && starts[[1L]] > 1L) records <- records[-1L]
    upis <- vapply(records, \(r) sub("^>([^[:space:]]+).*$", "\\1", r[[1L]]), "")
    target <- file.path(OUT_DIR, sprintf("%s_uniparc.fasta", upid))
    sorted <- as.character(unlist(records[order(upis, method = "radix")], use.names = FALSE))
    writeLinesBytes(sorted, target)
    recordCheck(run, "uniparc_records", upid, sprintf("%d records", length(starts)), length(starts) >= 1L)
    target
}

writeLinesBytes <- function(lines, path) {
    con <- file(path, open = "wb")
    on.exit(close(con))
    writeLines(lines, con, useBytes = TRUE)
}

compareAccessions <- function(run, upid, fasta, table) {
    fasta_ids <- sort(fastaAccessions(fasta), method = "radix")
    table_ids <- sort(table$accessions, method = "radix")
    value <- sprintf(
        "fasta %d; table %d; only in fasta %d; only in table %d",
        length(fasta_ids), length(table_ids),
        length(setdiff(fasta_ids, table_ids)), length(setdiff(table_ids, fasta_ids))
    )
    check <- paste0("accessions_", sub("^.*_([a-z]+)\\.tsv$", "\\1", table$path))
    recordCheck(run, check, upid, value, identical(fasta_ids, table_ids))
}

# --- Citation ---

captureCitation <- function(run) {
    page <- fetchUrl(run, CITATION_URL, tempfile())
    raw_text <- readChar(page$path, file.size(page$path), useBytes = TRUE)
    if (!length(raw_text)) raw_text <- ""
    parsed <- tryCatch(jsonlite::fromJSON(raw_text, simplifyVector = FALSE), error = \(e) NULL)
    is_json <- is.list(parsed) && is.character(parsed$content)
    text <- if (is_json) parsed$content else raw_text
    Encoding(text) <- "UTF-8"
    writeLinesBytes(text, file.path(OUT_DIR, "citation-page.txt"))
    cited <- pairCitation(text)
    cited$url <- CITATION_URL
    cited$page_format <- if (is_json) "json content field" else "text"
    cited$page_sha256 <- digest::digest(file = page$path, algo = "sha256")
    run$citation <- cited
    cat(sprintf("citation\t%s\t%s\t%s\n", cited$status, cited$title, cited$doi))
}

# A title of CITATION_TITLE's form followed, before the next such title, by a DOI; the highest year paired wins.
pairCitation <- function(text) {
    starts <- gregexpr(CITATION_TITLE, text, ignore.case = TRUE, perl = TRUE)[[1L]]
    empty <- data.frame(
        status = "not_confirmed", authors = NA, title = NA, journal = NA, volume = NA, pages = NA,
        year = NA, doi = NA, listed_text = NA
    )
    if (starts[[1L]] == -1L) return(empty)
    ends <- starts + attr(starts, "match.length") - 1L
    limits <- c(starts[-1L] - 1L, nchar(text))
    pairs <- lapply(seq_along(starts), \(i) {
        segment <- substr(text, ends[[i]] + 1L, limits[[i]])
        doi <- regmatches(segment, regexpr(DOI_PATTERN, segment, perl = TRUE))
        if (!length(doi)) return(NULL)
        doi <- sub("[.,;:]+$", "", doi)
        title <- substr(text, starts[[i]], ends[[i]])
        doi_end <- ends[[i]] + regexpr(doi, segment, fixed = TRUE) + nchar(doi) - 1L
        listed <- substr(text, starts[[i]], doi_end)
        before <- substr(text, max(1L, starts[[i]] - 300L), starts[[i]] - 1L)
        citationFields(title, doi, segment, listed, before)
    })
    pairs <- do.call(rbind, pairs)
    if (is.null(pairs)) return(empty)
    pairs[order(-as.integer(pairs$title_year))[[1L]], setdiff(names(pairs), "title_year")]
}

citationFields <- function(title, doi, segment, listed, before) {
    # journal, volume, optional issue, first and last page (a hyphen or an en dash between the pages)
    reference <- paste0(
        "([A-Z][A-Za-z.& ]+?)\\.?[[:space:]]+([0-9]{1,3})(\\([0-9]+\\))?:[[:space:]]*",
        "([A-Z]?[0-9]+)[[:space:]]*[-\\x{2013}]+[[:space:]]*([A-Z]?[0-9]+)"
    )
    found <- regmatches(segment, regexec(reference, segment, perl = TRUE))[[1L]]
    year <- regmatches(segment, regexec("\\((20[0-9]{2})\\)", segment))[[1L]]
    data.frame(
        status = "confirmed",
        authors = if (grepl("UniProt Consortium", before, ignore.case = TRUE)) "The UniProt Consortium" else NA,
        title = oneLine(gsub("[*_]", "", title)),
        journal = if (length(found)) trimws(found[[2L]]) else NA,
        volume = if (length(found)) found[[3L]] else NA,
        pages = if (length(found)) paste0(found[[5L]], "-", found[[6L]]) else NA,
        year = if (length(year)) year[[2L]] else NA,
        doi = doi,
        listed_text = oneLine(listed),
        title_year = sub("^.* ", "", title)
    )
}

# --- Hashes, size check, scan pre-count ---

dataFiles <- function() {
    files <- list.files(OUT_DIR, pattern = "^UP[0-9]+(\\.fasta|_xref\\.tsv|_go\\.tsv|_uniparc\\.fasta)$")
    sort(files, method = "radix")
}

hashTable <- function() {
    files <- dataFiles()
    paths <- file.path(OUT_DIR, files)
    entries <- vapply(paths, \(p) {
        lines <- readLines(p, warn = FALSE)
        if (grepl("\\.fasta$", p)) sum(startsWith(lines, ">")) else length(lines) - 1L
    }, 1L)
    data.frame(
        file = files, bytes = as.integer(file.size(paths)),
        sha256 = vapply(paths, \(p) digest::digest(file = p, algo = "sha256"), "", USE.NAMES = FALSE),
        entries = unname(entries)
    )
}

checkSizes <- function(run) {
    files <- list.files(OUT_DIR)
    for (file in sort(files, method = "radix")) {
        bytes <- as.integer(file.size(file.path(OUT_DIR, file)))
        recordCheck(run, "file_size", file, sprintf("%d bytes", bytes), bytes < MAX_BYTES)
    }
}

# The accession a line belongs to: the header's accession for a FASTA, the first field for a table.
lineAccessions <- function(file, lines) {
    if (grepl("\\.fasta$", file)) {
        headers <- startsWith(lines, ">")
        ids <- sub("^>([^[:space:]]+).*$", "\\1", lines[headers])
        ids <- vapply(strsplit(ids, "|", fixed = TRUE), \(x) x[[min(2L, length(x))]], "")
        return(c("before first header", ids)[cumsum(headers) + 1L])
    }
    if (grepl("^UP[0-9]+_(xref|go)\\.tsv$", file)) return(c("header line", sub("\t.*$", "", lines[-1L])))
    rep(NA_character_, length(lines))
}

scanPrecount <- function() {
    files <- sort(setdiff(list.files(OUT_DIR), "scan-precount.tsv"), method = "radix")
    rows <- lapply(files, \(file) {
        lines <- readLines(file.path(OUT_DIR, file), warn = FALSE)
        accessions <- lineAccessions(file, lines)
        do.call(rbind, lapply(names(PATTERNS), \(name) {
            hit <- grepl(PATTERNS[[name]], lines, useBytes = TRUE)
            ids <- unique(accessions[hit])
            data.frame(
                file = file, pattern = name, matching_lines = sum(hit),
                accessions = if (any(hit)) paste(ids, collapse = ";") else NA
            )
        }))
    })
    table <- do.call(rbind, rows)
    writeTsv(table, file.path(OUT_DIR, "scan-precount.tsv"))
    flagged <- table[table$matching_lines > 0L, c("file", "pattern", "matching_lines")]
    cat(sprintf("scan pre-count (scan.R at %s): %d file and pattern pairs match\n", SCAN_COMMIT, nrow(flagged)))
    if (nrow(flagged)) {
        cat(sprintf("precount\t%s\t%s\t%d\n", flagged$file, flagged$pattern, flagged$matching_lines), sep = "")
    }
}

# --- The capture ---

capture <- function(run) {
    release <- readFtpRelease(run)
    fastas <- findReferenceFastas(run)
    for (upid in REFERENCES) {
        fasta <- captureReferenceFasta(run, upid, fastas[[upid]])
        for (kind in c("xref", "go")) {
            table <- captureTable(run, upid, kind, if (kind == "xref") XREF_FIELDS else GO_FIELDS)
            recordCheck(
                run, paste0("release_uniprotkb_", kind), upid,
                sprintf("reply %s; ftp %s", table$release, release), identical(table$release, release)
            )
            if (!is.null(fasta)) compareAccessions(run, upid, fasta, table)
        }
    }
    for (upid in AUTHORS_REFERENCES) captureUniparc(run, upid)
    captureCitation(run)
}

writeRecords <- function(run) {
    hashes <- hashTable()
    writeTsv(hashes, file.path(OUT_DIR, "hashes.tsv"))
    cat(sprintf("hashes\t%s\t%s\t%s\t%s\n", hashes$file, hashes$bytes, hashes$sha256, hashes$entries), sep = "")
    releases <- do.call(rbind, run$releases)
    if (!is.null(releases)) {
        releases <- releases[order(releases$url, method = "radix"), ]
        writeTsv(releases, file.path(OUT_DIR, "release.tsv"))
        cat(sprintf("release\t%s\t%s\t%s\t%s\t%s\n", releases$url, releases$source, releases$release,
            releases$release_date, releases$used), sep = "")
    }
    if (!is.null(run$citation)) writeTsv(run$citation, file.path(OUT_DIR, "citation.tsv"))
    checkSizes(run)
    failed_downloads <- vapply(run$checks, \(row) row$check == "download" && row$result == "fail", TRUE)
    if (!any(failed_downloads)) {
        recordCheck(run, "download", "all", sprintf("%d downloads completed", run$downloads), TRUE)
    }
    times <- data.frame(
        check = c("job_start_utc", "job_end_utc"), reference = "all",
        value = c(run$start, utcNow()), result = "recorded"
    )
    checks <- rbind(do.call(rbind, run$checks), times)
    writeTsv(checks, file.path(OUT_DIR, "checks.tsv"))
    checks
}

main <- function() {
    if (dir.exists(OUT_DIR)) stop(sprintf("%s exists in %s: not overwritten", OUT_DIR, getwd()))
    dir.create(OUT_DIR)
    run <- newRun()
    cat(sprintf("KIT-167 capture; start %s (UTC)\n", run$start))
    tryCatch(
        capture(run),
        capture_stop = \(e) cat(sprintf("capture stopped: %s\n", conditionMessage(e))),
        error = \(e) recordCheck(run, "script_error", "all", conditionMessage(e), FALSE)
    )
    checks <- tryCatch(writeRecords(run), error = \(e) {
        recordCheck(run, "script_error", "records", conditionMessage(e), FALSE)
        checks <- do.call(rbind, run$checks)
        writeTsv(checks, file.path(OUT_DIR, "checks.tsv"))
        checks
    })
    scanPrecount()
    failed <- unique(checks$check[checks$result == "fail"])
    cat(sprintf("checks failed: %d (%s)\n", length(failed), paste(failed, collapse = ", ")))
    quit(save = "no", status = if (length(failed)) 1L else 0L)
}

main()
