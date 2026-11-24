# KIT-192 script S1: the transcriptomics DA export, base against head (ticket section 9 "S1", acceptance criterion A2).
#
# Usage, from the package root of the fix commit on a run-tests runner (full git history):
#   Rscript docs/qa/workflow-execution/kit/kit-192/export-routes-control.R <base commit> [<out dir>]
#   <out dir>  default docs/qa/workflow-execution/kit/kit-192/s1 (made absolute)
#
# The mechanism of KIT-186's ruv-routes-control.R (that file is not edited). Each side (base: a detached git worktree of
# <base commit>; head: this tree) runs in a fresh Rscript process with its working directory at its tree's root:
# pkgload::load_all('.', export_all = TRUE, helpers = FALSE), then the test helpers sourced into the global environment,
# and each route inside testthat::test_that(). Routes, each on a fresh module_ci_trx_project(module_ci_trx_corpus()),
# driven through the real normalisation and DA modules (module_ci_trx_tabs_journey):
#   1 gate off, automatic, in process          (legacy grid)
#   2 gate off, automatic, managed             (legacy grid)
#   3 gate off, manual, percentage 10, k 2     (legacy grid)
#   4 gate on, none, in process
#   5 gate on, none, managed
#   6 gate on, groups_all_features, k rule user, k 2, in process
#   7 gate on, first_pass_nulls, k rule user, k 1, managed
# Per route and side: errors; head; k; SHA-256 of the DA long table (da_genes_long); the notifications (masked by
# guiEntryMask, joined by " || "); the sorted relative paths of every file under da_output_dir and publication_graphs_dir
# with a size > 0 flag (a folder name of the form YYYYMMDD_HHMMSS masked); SHA-256 of every .csv and .tsv under
# da_output_dir; SHA-256 of ruv_optimization_summary.txt with its Timestamp: line removed, or "absent". PNG and PDF
# bytes are not compared (the PDF carries its creation time).
# Written: <out>/export-routes-<side>.tsv by each side, <out>/export-routes.tsv (both sides) and
# <out>/export-routes-compare.tsv (route, field, base, head, equal, criterion) here. Synthetic data only.
# Exit status non-zero when a side's process fails or a route fails to run on either side.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1L, grepl("^[0-9a-f]{40}$", args[[1L]]))
base_commit <- args[[1L]]
out <- if (length(args) >= 2L) args[[2L]] else "docs/qa/workflow-execution/kit/kit-192/s1"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
head_root <- normalizePath(".")
head_commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
cat("KIT-192 S1: head", head_commit, "base", base_commit, "out", out, "\n")
rscript <- file.path(R.home("bin"), "Rscript")

child <- file.path(tempdir(), "kit192-s1-child.R")
writeLines(c(
    "args <- commandArgs(trailingOnly = TRUE)",
    "setwd(args[[1L]])",
    "side <- args[[2L]]",
    "out <- args[[3L]]",
    "suppressPackageStartupMessages(pkgload::load_all('.', quiet = TRUE, export_all = TRUE, helpers = FALSE))",
    "suppressMessages(testthat::source_test_helpers('tests/testthat', env = globalenv()))",
    "Sys.setenv(MULTISCHOLAR_TRX_MANIFEST_NOTE = '1')",
    "sha <- function(x) digest::digest(serialize(x, NULL, version = 3L), algo = 'sha256', serialize = FALSE)",
    "flat <- function(x) gsub('[\\t\\r\\n]+', ' ', paste(as.character(x), collapse = ' | '))",
    "fileSha <- function(path) digest::digest(file = path, algo = 'sha256')",
    "listing <- function(paths) {",
    "    entries <- unlist(lapply(c(da = paths$da_output_dir, graphs = paths$publication_graphs_dir), function(dir) {",
    "        files <- sort(list.files(dir, recursive = TRUE, all.files = TRUE))",
    "        if (!length(files)) return(character())",
    "        sizes <- as.integer(file.size(file.path(dir, files)) > 0)",
    "        # A directory named by its creation time (setupDirectories' timestamped folders) is masked.",
    "        paste0(basename(dir), '/', gsub('[0-9]{8}_[0-9]{6}', '<time>', files), ':', sizes)",
    "    }), use.names = FALSE)",
    "    paste(entries, collapse = ';')",
    "}",
    "tableDigests <- function(dir) {",
    "    files <- sort(list.files(dir, pattern = '\\\\.(csv|tsv)$', recursive = TRUE))",
    "    paste(paste0(files, '=', vapply(file.path(dir, files), fileSha, character(1))), collapse = ';')",
    "}",
    "summaryDigest <- function(dir) {",
    "    path <- file.path(dir, 'ruv_optimization_summary.txt')",
    "    if (!file.exists(path)) return('absent')",
    "    lines <- readLines(path, warn = FALSE)",
    "    digest::digest(paste(lines[!startsWith(lines, 'Timestamp:')], collapse = '\\n'), algo = 'sha256', serialize = FALSE)",
    "}",
    "routes <- list(",
    "    `1` = list(gate = FALSE, managed = FALSE, norm = list(ruv_mode = 'automatic')),",
    "    `2` = list(gate = FALSE, managed = TRUE, norm = list(ruv_mode = 'automatic')),",
    "    `3` = list(gate = FALSE, managed = FALSE, norm = list(ruv_mode = 'manual', ruv_percentage = 10, ruv_k = 2)),",
    "    `4` = list(gate = TRUE, managed = FALSE, norm = list(ruv_mode = 'automatic', ruv_design = 'none')),",
    "    `5` = list(gate = TRUE, managed = TRUE, norm = list(ruv_mode = 'automatic', ruv_design = 'none')),",
    "    `6` = list(gate = TRUE, managed = FALSE, norm = list(ruv_mode = 'automatic',",
    "        ruv_design = 'groups_all_features', ruv_k_rule = 'user', ruv_k_user = 2)),",
    "    `7` = list(gate = TRUE, managed = TRUE, norm = list(ruv_mode = 'automatic',",
    "        ruv_design = 'first_pass_nulls', ruv_k_rule = 'user', ruv_k_user = 1)))",
    "rows <- list()",
    "failed <- FALSE",
    "for (route in names(routes)) {",
    "    r <- routes[[route]]",
    "    row <- NULL",
    "    ok <- tryCatch(testthat::test_that(paste('KIT-192 S1', side, 'route', route), {",
    "        withr::local_options(MultiScholaR.ruv_design_panel = r$gate)",
    "        fixture <- module_ci_trx_project(module_ci_trx_corpus())",
    "        if (isTRUE(r$managed)) {",
    "            probe <- new.env()",
    "            omicsManagedChainMocks(fixture, probe)",
    "        }",
    "        journey <- module_ci_trx_tabs_journey(fixture, managed = r$managed, norm = r$norm)",
    "        owner <- fixture$owner",
    "        roots <- tempdir()",
    "        record <- buildWorkflowAnalysisRecord(owner)",
    "        steps <- Filter(function(step) identical(step$operation, 'trx_ruvs'), record$steps)",
    "        ruv <- if (length(steps)) steps[[length(steps)]]$outcomes else list()",
    "        k <- ruv[['ruv_k', exact = TRUE]] %||% ruv[['k', exact = TRUE]] %||% NA",
    "        notes <- vapply(journey$notes, function(note) paste0('[', note$type %||% '', '] ', note$message), character(1))",
    "        row <<- data.frame(side = side, route = route, errors = flat(journey$errors), head = owner$current_state %||% '',",
    "            k = as.character(k),",
    "            table_sha256 = if (is.null(journey$da$da_genes_long)) '<no table>' else sha(journey$da$da_genes_long),",
    "            notifications = flat(guiEntryMask(paste(notes, collapse = ' || '), roots)),",
    "            files = listing(fixture$paths), table_digests = tableDigests(fixture$paths$da_output_dir),",
    "            summary_sha256 = summaryDigest(fixture$paths$da_output_dir), stringsAsFactors = FALSE)",
    "        testthat::succeed()",
    "    }), error = function(e) {",
    "        cat('route', route, 'error', class(e)[[1L]], conditionMessage(e), '\\n')",
    "        FALSE",
    "    })",
    "    if (is.null(row)) row <- data.frame(side = side, route = route, errors = '<route failed>', head = '', k = '',",
    "        table_sha256 = '', notifications = '', files = '', table_digests = '', summary_sha256 = '',",
    "        stringsAsFactors = FALSE)",
    "    cat('KIT-192 S1', side, 'route', route, if (isTRUE(ok) && !nzchar(row$errors)) 'ok' else 'FAILED', '\\n')",
    "    if (!isTRUE(ok) || identical(row$errors, '<route failed>')) failed <- TRUE",
    "    rows[[route]] <- row",
    "}",
    "utils::write.table(do.call(rbind, rows), file.path(out, paste0('export-routes-', side, '.tsv')), sep = '\\t',",
    "    quote = FALSE, row.names = FALSE)",
    "quit(status = if (failed) 1L else 0L)"), child)

runSide <- function(tree, side) {
    console <- suppressWarnings(system2(rscript, c(shQuote(child), shQuote(tree), side, shQuote(out)), stdout = TRUE,
        stderr = TRUE))
    cat(grep("^KIT-192 S1|^route ", console, value = TRUE), sep = "\n")
    cat(utils::tail(console, 15L), sep = "\n")
    status <- attr(console, "status")
    if (is.null(status)) 0L else status
}

base_tree <- file.path(tempdir(), "kit192-s1-base")
invisible(system2("git", c("worktree", "add", "--detach", shQuote(base_tree), base_commit), stdout = TRUE, stderr = TRUE))
stopifnot(dir.exists(file.path(base_tree, "tests", "testthat")))
base_status <- runSide(base_tree, "base")
head_status <- runSide(head_root, "head")

readSide <- function(side) {
    path <- file.path(out, paste0("export-routes-", side, ".tsv"))
    if (!file.exists(path)) return(NULL)
    utils::read.delim(path, colClasses = "character", quote = "", na.strings = character(), comment.char = "")
}
base <- readSide("base")
head <- readSide("head")
stopifnot(!is.null(base), !is.null(head))
utils::write.table(rbind(base, head), file.path(out, "export-routes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# A2 (ticket section 8). Routes 1 to 3 (legacy grid): every field identical. Routes 4 to 7 (no grid): (i) the base's last
# notification is the save warning and the head's the success message; (ii) the lists without that last element are
# identical; (iii) every other field identical, and the grid summary "absent" on both sides.
saved <- "[message] DA results saved to disk (tables, volcano plots, heatmaps)"
warned <- "[warning] Warning: Could not save results to disk: non-numeric argument to mathematical function"
fields <- c("errors", "head", "k", "table_sha256", "notifications", "files", "table_digests", "summary_sha256")
row <- function(route, field, b, h, equal, criterion) {
    data.frame(route = route, field = field, base = b, head = h, equal = equal, criterion = criterion, stringsAsFactors = FALSE)
}
last <- function(x) if (length(x)) x[[length(x)]] else ""
compare <- list()
for (route in base$route) {
    b <- base[base$route == route, , drop = FALSE]
    h <- head[head$route == route, , drop = FALSE]
    legacy <- route %in% c("1", "2", "3")
    for (field in fields) {
        if (!legacy && identical(field, "notifications")) {
            bn <- strsplit(b[[field]], " || ", fixed = TRUE)[[1L]]
            hn <- strsplit(h[[field]], " || ", fixed = TRUE)[[1L]]
            compare[[length(compare) + 1L]] <- row(route, "notifications_last", last(bn), last(hn),
                identical(last(bn), warned) && identical(last(hn), saved),
                "A2 (i): base last = save warning, head last = success message")
            compare[[length(compare) + 1L]] <- row(route, "notifications_rest", paste(utils::head(bn, -1L), collapse = " || "),
                paste(utils::head(hn, -1L), collapse = " || "), identical(utils::head(bn, -1L), utils::head(hn, -1L)),
                "A2 (ii): identical without the last element")
        } else if (!legacy && identical(field, "summary_sha256")) {
            compare[[length(compare) + 1L]] <- row(route, field, b[[field]], h[[field]],
                identical(b[[field]], "absent") && identical(h[[field]], "absent"), "A2 (iii): absent on both sides (D1 (a))")
        } else {
            compare[[length(compare) + 1L]] <- row(route, field, b[[field]], h[[field]], identical(b[[field]], h[[field]]),
                if (legacy) "A2 legacy: identical" else "A2 (iii): identical")
        }
    }
}
compare <- do.call(rbind, compare)
utils::write.table(compare, file.path(out, "export-routes-compare.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("KIT-192 S1: A2 rows", nrow(compare), "held", sum(compare$equal), "failing", paste(unique(paste0(compare$route[!compare$equal],
    ":", compare$field[!compare$equal])), collapse = ","), "\n")
cat("KIT-192 S1: base status", base_status, "head status", head_status, "\n")
quit(status = if (base_status != 0L || head_status != 0L) 1L else 0L)
