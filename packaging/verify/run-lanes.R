# KIT-134 P5 step 3 (the token build): GUI-entry lanes of tests/testthat/test-gui-entry-proteomics.R,
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

test_file <- file.path(checkout, "tests", "testthat", "test-gui-entry-proteomics.R")
blocks <- regmatches(readLines(test_file), regexpr("test_that\\(\"KIT-101 entry lane [^\"]+\"", readLines(test_file)))
blocks <- gsub("^test_that\\(\"|\"$", "", blocks)
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

summaries <- list()
failures <- list()
differences <- list()
for (lane in lanes) {
    block <- blocks[startsWith(blocks, sprintf("KIT-101 entry lane %s:", lane))]
    if (length(block) != 1L) {
        cat("no single test_that block for lane", lane, "\n")
        quit(save = "no", status = 2L)
    }
    unlink(file.path(out_dir, "records", paste0(lane, ".json")))
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
    if (file.exists(record_path)) {
        lane_differences <- tryCatch(compareLane(lane, jsonlite::read_json(record_path, simplifyVector = FALSE)$steps, oracle$steps),
            error = function(e) data.frame(lane = lane, step = "<comparison failed>", field = "", key = "", kind = "discrete",
                runner = truncate(conditionMessage(e)), oracle = "", stringsAsFactors = FALSE))
    }
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
quit(save = "no", status = if (passed) 0L else 1L)
