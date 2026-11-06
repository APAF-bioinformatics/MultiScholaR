# KIT-101: the per-file test baseline. It extends the series' one-process-per-file runner (kit/scripts-078/run1.R):
# each test file runs in its own R process, and for each file it records the static test_that blocks, the blocks run,
# the expectations passed, failed and skipped, the errors, a file-level skip and its reason (testthat reports a skip
# raised outside test_that() to the reporter only, so the list reporter alone shows such a file as empty), the hidden
# blocks, and every skip reason classified by regex. The run header records the git SHA, the R version, the library
# paths, the load mode and the sha256 of the sorted R/ file list.
#
# Load modes:
#   load_all  pkgload::load_all(<package source>, export_all = TRUE, helpers = TRUE), as run1.R;
#   check     the installed package (library(MultiScholaR) from R_LIBS) and test_file(package = "MultiScholaR",
#             load_package = "installed"), on a copy of tests/ inside the R CMD check directory, as test_check runs it.
#
# Usage:
#   Rscript tools/testing/test-baseline.R file <mode> <package source> <test file> <out dir>
#   Rscript tools/testing/test-baseline.R suite <mode> <package source> <tests/testthat dir> <out dir> [parallel = 6] [per-file timeout s = 1800]
#   Rscript tools/testing/test-baseline.R aggregate <mode> <package source> <tests/testthat dir> <out dir>
# `suite` writes <out dir>/rows/<file>.tsv, <out dir>/skips/<file>.tsv and <out dir>/logs/<file>.log, then aggregates
# into <out dir>/test-baseline.tsv, <out dir>/skip-reasons.tsv and <out dir>/header.tsv.

baselineSkipCategory <- function(message) {
    rules <- c(
        binding_guard = "extracted helper\\(s\\) not present|helper\\(s\\) not present|public binding|binding\\(s\\) not present",
        file_guard = "file\\(s\\) not present|Target-only|duplicate\\(s\\) still present|split file",
        quarantine = "quarantin|\\bQ-[0-9]+",
        environment = "Linux|not installed|installed|capable|browser|chromote|On CRAN|not available|interactive",
        fixture = "LFS pointer|[Ss]napshot|not found|corpus generator|source tree|fixture")
    for (category in names(rules)) if (grepl(rules[[category]], message)) return(category)
    "other"
}

baselineStaticBlocks <- function(path) {
    lines <- readLines(path, warn = FALSE)
    sum(grepl("^\\s*(testthat::)?test_that\\(", lines))
}

baselineFileRun <- function(mode, package, path, out_dir) {
    started <- Sys.time()
    if (identical(mode, "load_all")) {
        suppressMessages(pkgload::load_all(package, quiet = TRUE, export_all = TRUE, helpers = TRUE))
    } else if (identical(mode, "check")) {
        suppressMessages(library(MultiScholaR))
    } else stop("unknown load mode: ", mode)
    # Every result the file raises, inside or outside test_that(): the list reporter drops a top-level skip.
    Log <- R6::R6Class("BaselineLog", inherit = testthat::Reporter, public = list(entries = list(),
        add_result = function(context, test, result) {
            self$entries[[length(self$entries) + 1L]] <- list(test = if (is.null(test)) NA_character_ else test,
                class = class(result)[[1L]], message = gsub("[\t\r\n]+", " ", conditionMessage(result)))
        }))
    log <- Log$new()
    list_reporter <- testthat::ListReporter$new()
    reporter <- testthat::MultiReporter$new(list(list_reporter, log))
    outcome <- tryCatch({
        if (identical(mode, "check")) {
            testthat::test_file(path, reporter = reporter, stop_on_failure = FALSE, package = "MultiScholaR", load_package = "installed")
        } else testthat::test_file(path, reporter = reporter, stop_on_failure = FALSE)
        NULL
    }, error = function(error) conditionMessage(error))
    df <- as.data.frame(list_reporter$get_results())
    entries <- log$entries
    tests <- vapply(entries, `[[`, "", "test")
    classes <- vapply(entries, `[[`, "", "class")
    messages <- vapply(entries, `[[`, "", "message")
    top <- is.na(tests)
    file_skip <- any(top & classes == "expectation_skip")
    file_error <- any(top & classes %in% c("expectation_error", "expectation_failure"))
    file_skip_reason <- if (file_skip) sub("^Reason: ", "", messages[top & classes == "expectation_skip"][[1L]]) else ""
    file_error_message <- if (file_error) messages[top & classes %in% c("expectation_error", "expectation_failure")][[1L]] else outcome %||% ""
    static <- baselineStaticBlocks(path)
    blocks_run <- length(unique(tests[!top]))
    skips <- classes == "expectation_skip"
    reasons <- sub("^Reason: ", "", messages[skips])
    categories <- vapply(reasons, baselineSkipCategory, "", USE.NAMES = FALSE)
    counts <- table(factor(categories, levels = c("file_guard", "binding_guard", "environment", "quarantine", "fixture", "other")))
    row <- data.frame(file = basename(path), mode = mode, static_blocks = static, blocks_run = blocks_run,
        hidden_blocks = if (file_skip || file_error || !is.null(outcome)) max(static - blocks_run, 0L) else 0L,
        expectations = sum(df$nb), passed = sum(df$passed), failed = sum(df$failed),
        error_blocks = sum(df$error), skipped_blocks = sum(df$skipped), warnings = sum(df$warning),
        skip_results = sum(skips), file_skip = file_skip, file_skip_reason = file_skip_reason,
        file_skip_category = if (file_skip) baselineSkipCategory(file_skip_reason) else "",
        file_error = file_error || !is.null(outcome), file_error_message = substr(file_error_message, 1L, 300L),
        skip_file_guard = counts[["file_guard"]], skip_binding_guard = counts[["binding_guard"]],
        skip_environment = counts[["environment"]], skip_quarantine = counts[["quarantine"]],
        skip_fixture = counts[["fixture"]], skip_other = counts[["other"]],
        seconds = round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1), status = "complete",
        stringsAsFactors = FALSE)
    name <- sub("[.]R$", "", basename(path))
    utils::write.table(row, file.path(out_dir, "rows", paste0(name, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
    skip_rows <- data.frame(file = rep(basename(path), sum(skips)), mode = rep(mode, sum(skips)),
        test = ifelse(is.na(tests[skips]), "<file level>", tests[skips]), category = categories,
        reason = substr(reasons, 1L, 300L), stringsAsFactors = FALSE)
    utils::write.table(skip_rows, file.path(out_dir, "skips", paste0(name, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("COUNTS %s pass %d fail %d error %d skip %d file_skip %s\n", basename(path), row$passed, row$failed,
        row$error_blocks, row$skipped_blocks, file_skip))
    invisible(row)
}

baselineHeader <- function(mode, package, out_dir) {
    sha <- tryCatch(system2("git", c("-C", package, "rev-parse", "HEAD"), stdout = TRUE), error = function(e) NA_character_)
    r_files <- sort(list.files(file.path(package, "R"), pattern = "[.][Rr]$"))
    header <- data.frame(field = c("sha", "r_version", "load_mode", "r_files", "r_file_list_sha256", "lib_paths", "not_cran",
        "strict_skips", "generated"), value = c(sha[[1L]], R.version.string, mode, length(r_files),
        digest::digest(paste(r_files, collapse = "\n"), algo = "sha256", serialize = FALSE),
        paste(.libPaths(), collapse = ":"), Sys.getenv("NOT_CRAN"), Sys.getenv("MULTISCHOLAR_STRICT_SKIPS"),
        format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")), stringsAsFactors = FALSE)
    utils::write.table(header, file.path(out_dir, "header.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    header
}

baselineAggregate <- function(mode, package, tests_dir, out_dir) {
    files <- sort(list.files(tests_dir, pattern = "^test-.*[.]R$", full.names = TRUE))
    rows <- lapply(files, function(path) {
        name <- sub("[.]R$", "", basename(path))
        row_path <- file.path(out_dir, "rows", paste0(name, ".tsv"))
        if (file.exists(row_path)) return(utils::read.delim(row_path, stringsAsFactors = FALSE, quote = "", na.strings = character()))
        # No row: the process was killed (per-file timeout or the run's ceiling) or crashed before writing.
        data.frame(file = basename(path), mode = mode, static_blocks = baselineStaticBlocks(path), blocks_run = 0L,
            hidden_blocks = baselineStaticBlocks(path), expectations = 0L, passed = 0L, failed = 0L, error_blocks = 0L,
            skipped_blocks = 0L, warnings = 0L, skip_results = 0L, file_skip = FALSE, file_skip_reason = "",
            file_skip_category = "", file_error = TRUE, file_error_message = "no result row (timeout, ceiling or crash; see logs)",
            skip_file_guard = 0L, skip_binding_guard = 0L, skip_environment = 0L, skip_quarantine = 0L, skip_fixture = 0L,
            skip_other = 0L, seconds = NA_real_, status = "no_result", stringsAsFactors = FALSE)
    })
    table <- do.call(rbind, rows)
    utils::write.table(table, file.path(out_dir, "test-baseline.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    skip_files <- file.path(out_dir, "skips", paste0(sub("[.]R$", "", basename(files)), ".tsv"))
    skips <- do.call(rbind, lapply(skip_files[file.exists(skip_files)], function(path)
        utils::read.delim(path, stringsAsFactors = FALSE, quote = "", na.strings = character())))
    if (!is.null(skips)) utils::write.table(skips, file.path(out_dir, "skip-reasons.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    baselineHeader(mode, package, out_dir)
    cat(sprintf("TOTAL files %d pass %d fail %d error %d skip %d file_skips %d file_errors %d no_result %d hidden_blocks %d\n",
        nrow(table), sum(table$passed), sum(table$failed), sum(table$error_blocks), sum(table$skipped_blocks),
        sum(as.logical(table$file_skip)), sum(as.logical(table$file_error)), sum(table$status == "no_result"),
        sum(table$hidden_blocks)))
    invisible(table)
}

baselineSuite <- function(script, mode, package, tests_dir, out_dir, parallel = 6L, timeout = 1800L) {
    for (sub_dir in c("rows", "skips", "logs")) dir.create(file.path(out_dir, sub_dir), recursive = TRUE, showWarnings = FALSE)
    files <- sort(list.files(tests_dir, pattern = "^test-.*[.]R$", full.names = TRUE))
    invisible(parallel::mclapply(files, function(path) {
        log <- file.path(out_dir, "logs", paste0(sub("[.]R$", "", basename(path)), ".log"))
        system2("timeout", c(as.character(timeout), file.path(R.home("bin"), "Rscript"), shQuote(script), "file", mode,
            shQuote(package), shQuote(path), shQuote(out_dir)), stdout = log, stderr = log)
    }, mc.cores = parallel, mc.preschedule = FALSE))
    baselineAggregate(mode, package, tests_dir, out_dir)
}

if (sys.nframe() == 0L) {
    `%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
    args <- commandArgs(TRUE)
    if (length(args) < 5L) stop("usage: test-baseline.R <file|suite|aggregate> <load_all|check> <package source> <test file|tests dir> <out dir> [parallel] [timeout]")
    command <- args[[1L]]; mode <- args[[2L]]; package <- normalizePath(args[[3L]]); target <- normalizePath(args[[4L]])
    out_dir <- args[[5L]]
    for (sub_dir in c("rows", "skips", "logs")) dir.create(file.path(out_dir, sub_dir), recursive = TRUE, showWarnings = FALSE)
    out_dir <- normalizePath(out_dir)
    script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]))
    switch(command,
        file = baselineFileRun(mode, package, target, out_dir),
        suite = baselineSuite(script, mode, package, target, out_dir, if (length(args) > 5L) as.integer(args[[6L]]) else 6L,
            if (length(args) > 6L) as.integer(args[[7L]]) else 1800L),
        aggregate = baselineAggregate(mode, package, target, out_dir),
        stop("unknown command: ", command))
}
