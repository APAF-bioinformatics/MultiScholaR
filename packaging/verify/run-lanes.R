# KIT-134 P5 step 3 (the token build): GUI-entry lanes of tests/testthat/test-gui-entry-<omic>.R,
# in process, against the bundle's installed MultiScholaR and library, with the test tools in a
# separate library after the bundle's. Build-time verification, not package code.
#
# Usage: Rscript --vanilla run-lanes.R <checkout> <test-tool library> <output directory> <lane id>...
# Run with the bundle's launcher (bin/Rscript), so R, its home and its library are the bundle's.
#
# Each lane runs as the test file runs it (testthat::test_file, the installed package, the lane's
# test_that block only), against its frozen oracle (tests/testdata/gui-entry/<lane>.json, recorded
# on the laptop). The lane's record (MULTISCHOLAR_GUI_ENTRY_RECORD_DIR) is then compared with the
# oracle step by step and key by key, and every difference is written with its kind:
#   discrete  a head, history, tab status, notification, refusal check or text value: compared
#             exactly under revision 1.2's rule;
#   digest    a sha256 the record holds instead of values (state contents rounded to 10 significant
#             digits, tables, widgets, long text, artifacts): the rule's tolerance (relative 1e-6,
#             absolute 1e-9) cannot be applied to a digest; the values are needed on both machines.
# Writes lane-summary.tsv, lane-differences.tsv and lane-failures.tsv (text only). Exits 1 when a
# lane does not pass, after writing them; 2 when the lanes could not be run.
#
# With MSR_REPORT_CHECK=1 (the installer stage; W1 in outcome.md) the question is another one: does a
# Summary report render to .docx with every Chrome-type browser hidden? The caller hides them
# (none on the PATH, CHROMOTE_CHROME at a path that does not exist); this script refuses to start
# (2) if one can still be found, prints every rmarkdown render with its output's size, and reads
# each lane's record at the step summary_report: the report's .docx among the artifacts and the
# application's own "Report generated successfully!". It then exits 0 when at least one lane
# rendered its report and 1 when none did; the lanes' comparison with their oracles is still
# written and printed, and does not decide the exit status in this mode (that is P5's).

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) < 4L) stop("usage: Rscript --vanilla run-lanes.R <checkout> <test-tool library> <output directory> <lane id>...")
checkout <- normalizePath(arguments[[1L]], mustWork = TRUE)
tools_lib <- normalizePath(arguments[[2L]], mustWork = TRUE)
out_dir <- arguments[[3L]]
lanes <- arguments[-(1:3)]
dir.create(file.path(out_dir, "records"), recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)

bundle_libs <- .libPaths()
.libPaths(c(bundle_libs, tools_lib))
shadowing <- intersect(rownames(utils::installed.packages(tools_lib)), rownames(utils::installed.packages(bundle_libs)))
cat(sprintf("library paths: %s\n", paste(.libPaths(), collapse = " | ")))
if (length(shadowing)) {
    cat("the test-tool library holds packages the bundle also has:", paste(shadowing, collapse = ", "), "\n")
    quit(save = "no", status = 2L)
}
cat(sprintf("test-tool library: %s\n", paste(sort(rownames(utils::installed.packages(tools_lib))), collapse = ", ")))
suppressPackageStartupMessages(library(MultiScholaR))
versions <- vapply(c("MultiScholaR", "testthat", "fastduckdb", "fastlimpadpc", "limpa", "duckdb"), function(package) {
    path <- find.package(package, quiet = TRUE)
    if (length(path)) sprintf("%s %s (%s)", package, utils::packageVersion(package), dirname(path)) else paste(package, "absent")
}, character(1))
cat(versions, sep = "\n")
cat(sprintf("locale: %s\n", Sys.getlocale()))
reached <- tryCatch({
    close(socketConnection("cran.r-project.org", 80L, open = "r+", blocking = TRUE, timeout = 10))
    TRUE
}, error = function(e) FALSE, warning = function(w) FALSE)
cat(sprintf("network: %s\n", if (reached) "OPEN (cran.r-project.org reached)" else "closed (cran.r-project.org not reached)"))
cat(sprintf("DPC engine resolved by default: %s\n", tryCatch(MultiScholaR:::resolveLimpaDpcEngine(), error = conditionMessage)))

report_check <- identical(Sys.getenv("MSR_REPORT_CHECK"), "1")
if (report_check) {
    browsers <- c("google-chrome", "google-chrome-stable", "chrome", "chromium", "chromium-browser", "chrome-headless-shell",
        "microsoft-edge", "microsoft-edge-stable", "msedge")
    on_path <- Sys.which(browsers)
    chromote_chrome <- Sys.getenv("CHROMOTE_CHROME")
    found <- tryCatch(as.character(chromote::find_chrome()), error = function(e) paste("error:", conditionMessage(e)))
    cat(sprintf("BROWSERS PATH: %s\n", Sys.getenv("PATH")))
    cat(sprintf("BROWSERS on the PATH: %s\n", if (any(nzchar(on_path))) paste(on_path[nzchar(on_path)], collapse = ", ") else "none"))
    cat(sprintf("BROWSERS CHROMOTE_CHROME: '%s' (exists: %s); chromote::find_chrome(): '%s'\n", chromote_chrome,
        file.exists(chromote_chrome), paste(found, collapse = " ")))
    if (any(nzchar(on_path)) || !nzchar(chromote_chrome) || file.exists(chromote_chrome) || !identical(found, chromote_chrome)) {
        cat("a browser can still be found: the report check does not start\n")
        quit(save = "no", status = 2L)
    }
    # Every rmarkdown render of the lanes, with what it wrote (the lane's folder is removed when the lane ends).
    traced <- tryCatch({
        suppressMessages(trace("render", where = asNamespace("rmarkdown"), print = FALSE, exit = quote(try({
            written <- returnValue()
            if (is.character(written) && length(written) == 1L && file.exists(written)) {
                cat(sprintf("RENDER %s: %.0f bytes, zip container %s, %s\n", basename(written), file.size(written),
                    identical(readBin(written, "raw", 2L), charToRaw("PK")), format(Sys.time(), "%H:%M:%S")), file = stderr())
            } else cat("RENDER returned no file\n", file = stderr())
        }, silent = TRUE))))
        TRUE
    }, error = function(e) {
        cat("rmarkdown::render is not traced:", conditionMessage(e), "\n")
        FALSE
    })
    cat(sprintf("rmarkdown::render traced: %s\n", traced))
}

test_files <- list.files(file.path(checkout, "tests", "testthat"), pattern = "^test-gui-entry-.*\\.R$", full.names = TRUE)
blocks <- unlist(lapply(test_files, function(file) {
    found <- regmatches(readLines(file), regexpr("test_that\\(\"KIT-101 entry lane [^\"]+\"", readLines(file)))
    stats::setNames(gsub("^test_that\\(\"|\"$", "", found), rep(file, length(found)))
}))
Sys.setenv(MULTISCHOLAR_GUI_ENTRY_RECORD_DIR = file.path(out_dir, "records"), NOT_CRAN = "true")
truncate <- function(x, n = 160L) {
    x <- gsub("[\t\r\n]+", " ", paste(as.character(unlist(x)), collapse = " "))
    if (nchar(x) > n) paste0(substr(x, 1L, n), "...") else x
}

#' The differences between a lane's recorded steps and its oracle's, one row per differing key.
compareLane <- function(lane, recorded, frozen) {
    isDigest <- function(value) is.character(value) && length(value) == 1L &&
        grepl("^((<table>|<widget>|<text>|<ui>) )?[0-9a-f]{64}$", value)
    rows <- list()
    for (step in union(names(frozen), names(recorded))) {
        ours <- recorded[[step]]
        theirs <- frozen[[step]]
        if (is.null(ours) || is.null(theirs)) {
            rows[[length(rows) + 1L]] <- data.frame(lane = lane, step = step, field = "<step>", key = "", kind = "discrete",
                runner = if (is.null(ours)) "<absent>" else "<present>", oracle = if (is.null(theirs)) "<absent>" else "<present>",
                stringsAsFactors = FALSE)
            next
        }
        for (field in union(names(ours), names(theirs))) {
            a <- ours[[field]]
            b <- theirs[[field]]
            if (identical(a, b)) next
            keyed <- is.list(a) && is.list(b) && !is.null(names(a)) && !is.null(names(b))
            for (key in if (keyed) union(names(a), names(b)) else "") {
                x <- if (keyed) a[[key]] else a
                y <- if (keyed) b[[key]] else b
                if (identical(x, y)) next
                rows[[length(rows) + 1L]] <- data.frame(lane = lane, step = step, field = field, key = truncate(key, 200L),
                    kind = if (isDigest(x) && isDigest(y)) "digest" else "discrete",
                    runner = if (is.null(x)) "<absent>" else truncate(x), oracle = if (is.null(y)) "<absent>" else truncate(y),
                    stringsAsFactors = FALSE)
            }
        }
    }
    do.call(rbind, rows)
}

# KIT-134 open item (run 44241609555): prot_lfq_fragpipe's da_proteins_groupKO-groupWT_long_annot.tsv differed from
# its oracle by digest on Windows. Each copy the recorder reads (helper-gui-entry.R reads every .tsv with readLines) is
# kept in the output folder, which never leaves the job, and only fingerprints are printed (FPRINT lines): the bytes,
# line endings, the header, per column its kind, a digest of its text and, for numbers, digests of the values rounded
# to 15, 12, 10, 8 and 6 significant digits, two sample cells, so two platforms' tables can be compared from the logs.
# The recorder's digest (`digest(lines, serialize = FALSE)`) hashes the first line only, so the oracle's value is
# compared with the header's.
capture_name <- "da_proteins_groupKO-groupWT_long_annot.tsv"
capture <- new.env()
capture$lane <- ""
capture$n <- 0L
dir.create(file.path(out_dir, "captured"), showWarnings = FALSE)
captureTable <- function(con) {
    if (is.character(con) && length(con) == 1L && identical(basename(con), capture_name) && file.exists(con)) {
        capture$n <- capture$n + 1L
        file.copy(con, file.path(out_dir, "captured", sprintf("%s-%02d.tsv", capture$lane, capture$n)))
    }
    invisible()
}
invisible(suppressMessages(trace("readLines", tracer = quote(captureTable(con)), where = baseenv(), print = FALSE)))

fingerprintTables <- function(lane, oracle) {
    files <- sort(list.files(file.path(out_dir, "captured"), pattern = paste0("^", lane, "-[0-9]+\\.tsv$"), full.names = TRUE))
    if (!length(files)) return(invisible())
    stored <- unique(unlist(lapply(oracle$steps, function(step) step$artifacts[grepl(capture_name, names(step$artifacts), fixed = TRUE)])))
    hashes <- vapply(files, function(file) digest::digest(file = file, algo = "sha256"), "")
    cat(sprintf("FPRINT %s copies read by the recorder: %d; distinct contents: %d; the oracle's digests: %s\n", lane, length(files),
        length(unique(hashes)), paste(substr(stored, 1L, 16L), collapse = ", ")))
    short <- function(x) substr(x, 1L, 16L)
    for (file in files[!duplicated(hashes)]) {
        label <- sprintf("%s %s", lane, sub("\\.tsv$", "", basename(file)))
        raw <- readBin(file, "raw", file.size(file))
        lf <- which(raw == as.raw(10L))
        crlf <- sum(lf > 1L & raw[pmax(lf - 1L, 1L)] == as.raw(13L))
        lines <- readLines(file, warn = FALSE)
        header <- digest::digest(lines, algo = "sha256", serialize = FALSE)
        cat(sprintf("FPRINT %s: bytes %d, sha256 %s; lines %d, LF %d, CRLF %d, UTF-8 BOM %s; first-line digest %s (the oracle's: %s)\n",
            label, length(raw), short(hashes[[file]]), length(lines), length(lf), crlf,
            identical(raw[1:3], as.raw(c(0xef, 0xbb, 0xbf))), short(header), if (header %in% stored) "equal" else "differs"))
        cat(sprintf("FPRINT %s header: %s\n", label, substr(lines[[1L]], 1L, 1500L)))
        table <- utils::read.delim(file, colClasses = "character", check.names = FALSE, na.strings = character(), quote = "\"",
            comment.char = "")
        cat(sprintf("FPRINT %s shape: %d rows, %d columns; first column order %s, sorted %s; row-text digest %s\n", label, nrow(table),
            ncol(table), short(digest::digest(table[[1L]])), short(digest::digest(sort(table[[1L]]))),
            short(digest::digest(unname(as.list(table))))))
        for (j in seq_along(table)) {
            text <- table[[j]]
            value <- suppressWarnings(as.numeric(text))
            numeric <- any(!is.na(value)) && all(!is.na(value) | text %in% c("", "NA", "NaN", "Inf", "-Inf"))
            rounded <- if (numeric) vapply(c(15L, 12L, 10L, 8L, 6L), function(s)
                short(digest::digest(formatC(signif(value, s), digits = s, format = "g"))), "") else character()
            cat(sprintf("FPRINT %s col %d %s: %s; text %s%s; sample %s\n", label, j, names(table)[[j]],
                if (numeric) sprintf("numeric (%d NA)", sum(is.na(value))) else "text", short(digest::digest(text)),
                if (numeric) paste0("; rounded 15/12/10/8/6 ", paste(rounded, collapse = "/")) else "",
                paste(substr(utils::head(text, 2L), 1L, 40L), collapse = " | ")))
        }
    }
    invisible()
}

summaries <- list()
reports <- list()
failures <- list()
differences <- list()
for (lane in lanes) {
    block <- blocks[startsWith(blocks, sprintf("KIT-101 entry lane %s:", lane))]
    if (length(block) != 1L) {
        cat("no single test_that block for lane", lane, "\n")
        quit(save = "no", status = 2L)
    }
    test_file <- names(block)
    block <- unname(block)
    if (report_check) reports[[lane]] <- FALSE
    unlink(file.path(out_dir, "records", paste0(lane, ".json")))
    capture$lane <- lane
    capture$n <- 0L
    started <- Sys.time()
    reporter <- testthat::ListReporter$new()
    outcome <- tryCatch({
        testthat::test_file(test_file, desc = block, reporter = reporter, stop_on_failure = FALSE,
            package = "MultiScholaR", load_package = "installed")
        ""
    }, error = function(e) conditionMessage(e))
    expectations <- unlist(lapply(reporter$get_results(), function(test) test$results), recursive = FALSE)
    classes <- vapply(expectations, function(x) class(x)[[1L]], character(1))
    for (x in expectations[!classes %in% "expectation_success"]) {
        failures[[length(failures) + 1L]] <- data.frame(lane = lane, class = class(x)[[1L]],
            message = truncate(conditionMessage(x), 2000L), stringsAsFactors = FALSE)
    }
    record_path <- file.path(out_dir, "records", paste0(lane, ".json"))
    oracle <- jsonlite::read_json(file.path(checkout, "tests", "testdata", "gui-entry", paste0(lane, ".json")), simplifyVector = FALSE)
    lane_differences <- NULL
    if (file.exists(record_path) && report_check) {
        step <- jsonlite::read_json(record_path, simplifyVector = FALSE)$steps$summary_report
        documents <- grep("_report_[^/]*\\.docx$", names(step$artifacts), value = TRUE)
        notes <- vapply(step$notifications, function(note) paste0(note$type, ": ", note$message), character(1))
        reports[[lane]] <- length(documents) > 0L && any(grepl("Report generated successfully", notes, fixed = TRUE))
        cat(sprintf("REPORT %s: %s; .docx: %s; notifications: %s\n", lane, if (reports[[lane]]) "rendered" else "NOT rendered",
            if (length(documents)) paste(documents, collapse = ", ") else "none", truncate(paste(notes, collapse = " | "), 600L)))
    }
    if (file.exists(record_path)) {
        lane_differences <- tryCatch(compareLane(lane, jsonlite::read_json(record_path, simplifyVector = FALSE)$steps, oracle$steps),
            error = function(e) data.frame(lane = lane, step = "<comparison failed>", field = "", key = "", kind = "discrete",
                runner = truncate(conditionMessage(e)), oracle = "", stringsAsFactors = FALSE))
    }
    tryCatch(fingerprintTables(lane, oracle), error = function(e) cat(sprintf("FPRINT %s failed: %s\n", lane, conditionMessage(e))))
    differences <- c(differences, list(lane_differences))
    kinds <- if (is.null(lane_differences)) character() else lane_differences$kind
    summaries[[lane]] <- data.frame(lane = lane,
        passed = sum(classes == "expectation_success"), failed = sum(classes == "expectation_failure"),
        errors = sum(classes == "expectation_error"), skipped = sum(classes == "expectation_skip"),
        record = file.exists(record_path), oracle_steps = length(oracle$steps),
        steps_differing = if (is.null(lane_differences)) 0L else length(unique(lane_differences$step)),
        discrete_differences = sum(kinds == "discrete"), digest_differences = sum(kinds == "digest"),
        seconds = round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1),
        run_error = truncate(outcome, 300L), stringsAsFactors = FALSE)
    cat(sprintf("LANE %s: pass %d fail %d error %d skip %d; record %s; differing steps %d (discrete %d, digest %d); %.0f s\n",
        lane, summaries[[lane]]$passed, summaries[[lane]]$failed, summaries[[lane]]$errors, summaries[[lane]]$skipped,
        summaries[[lane]]$record, summaries[[lane]]$steps_differing, summaries[[lane]]$discrete_differences,
        summaries[[lane]]$digest_differences, summaries[[lane]]$seconds))
}

summary <- do.call(rbind, summaries)
utils::write.table(summary, file.path(out_dir, "lane-summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(if (length(failures)) do.call(rbind, failures) else data.frame(lane = character(), class = character(), message = character()),
    file.path(out_dir, "lane-failures.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
differences <- do.call(rbind, differences)
utils::write.table(if (is.null(differences)) data.frame(lane = character(), step = character(), field = character(), key = character(),
    kind = character(), runner = character(), oracle = character()) else differences,
    file.path(out_dir, "lane-differences.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
for (row in seq_len(nrow(summary))) if (summary$failed[[row]] + summary$errors[[row]] > 0L || !summary$record[[row]]) {
    cat(sprintf("first problems of %s:\n", summary$lane[[row]]))
    messages <- Filter(function(x) identical(x$lane, summary$lane[[row]]), failures)
    for (x in utils::head(messages, 5L)) cat("  ", x$class, ":", substr(x$message, 1L, 600L), "\n")
}
passed <- all(summary$failed == 0L & summary$errors == 0L & summary$passed > 0L)
cat(if (passed) "KIT-134 LANES PASSED" else "KIT-134 LANES FAILED", "\n")
if (report_check) {
    rendered <- names(reports)[vapply(reports, isTRUE, logical(1))]
    cat(if (length(rendered)) sprintf("KIT-134 REPORT CHECK PASSED: rendered to .docx with no browser to be found (%s)\n",
        paste(rendered, collapse = ", ")) else "KIT-134 REPORT CHECK FAILED: no lane rendered its report\n")
    quit(save = "no", status = if (length(rendered)) 0L else 1L)
}
quit(save = "no", status = if (passed) 0L else 1L)
