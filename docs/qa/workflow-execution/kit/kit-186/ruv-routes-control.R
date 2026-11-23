# KIT-186 script S1: RUV routes unchanged, base against head (ticket section 9 "S1", acceptance criterion A2).
#
# Usage, from the package root of the fix commit on a run-tests runner (full git history):
#   Rscript docs/qa/workflow-execution/kit/kit-186/ruv-routes-control.R <base commit> [<out dir>]
#   <out dir>  default docs/qa/workflow-execution/kit/kit-186/s1 (made absolute)
#
# Each side (base: a git worktree of <base commit>; head: this tree) runs in a fresh Rscript process with its working
# directory at its tree's root: pkgload::load_all(export_all = TRUE), the test helpers sourced into the global
# environment, and each route inside testthat::test_that() (C collation, as the test files; KIT-180 R4). Routes, each on
# a fresh module_ci_trx_project(module_ci_trx_corpus()), driven through the real normalisation, DA and summary modules
# (module_ci_trx_tabs_journey):
#   1 gate off, automatic, in process, with the summary steps (the report's text)
#   2 gate off, automatic, managed (omicsManagedChainMocks)
#   3 gate off, manual, percentage 10, k 2
#   4 gate on, groups_all_features, k rule user, k 2
#   5 gate on, groups_all_features, k rule num_sv_be (if k resolves to 0 the route is a k = 0 route, not part of A2)
#   6 gate on, first_pass_nulls, k rule user, k 1
# Per route and side: head, k, W_ columns, SHA-256 of the DA long table (da_genes_long), the record's formula_used, the
# model text, the notifications (masked by guiEntryMask) and, on route 1, the SHA-256 of the report's extracted text
# (white space collapsed to one space, then masked). Written: <out>/ruv-routes-<side>.tsv by each side,
# <out>/ruv-routes.tsv (both sides) and <out>/ruv-routes-compare.tsv (route, field, base, head, equal, a2) here, and
# <out>/report-text-<side>.txt (route 1's masked text). Synthetic data only.
# Exit status non-zero when a side's process fails or a route fails to run on either side.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1L, grepl("^[0-9a-f]{40}$", args[[1L]]))
base_commit <- args[[1L]]
out <- if (length(args) >= 2L) args[[2L]] else "docs/qa/workflow-execution/kit/kit-186/s1"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
head_root <- normalizePath(".")
head_commit <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
cat("KIT-186 S1: head", head_commit, "base", base_commit, "out", out, "\n")
rscript <- file.path(R.home("bin"), "Rscript")

child <- file.path(tempdir(), "kit186-s1-child.R")
writeLines(c(
    "args <- commandArgs(trailingOnly = TRUE)",
    "setwd(args[[1L]])",
    "side <- args[[2L]]",
    "out <- args[[3L]]",
    "suppressPackageStartupMessages(pkgload::load_all('.', quiet = TRUE, export_all = TRUE, helpers = FALSE))",
    "suppressMessages(testthat::source_test_helpers('tests/testthat', env = globalenv()))",
    "Sys.setenv(MULTISCHOLAR_TRX_MANIFEST_NOTE = '1')",
    "# MCI-032's summary steps and docx text (tests/testthat/test-module-ci-trx-summary.R:13-36), not helpers.",
    "summarySteps <- function(session, output, out) {",
    "    session$setInputs(`summary-experiment_label` = 'kit097', `summary-description` = 'KIT-097 synthetic transcriptomics')",
    "    session$setInputs(`summary-save_workflow_args` = 1L)",
    "    session$flushReact()",
    "    session$setInputs(`summary-copy_to_publication` = 1L)",
    "    session$flushReact()",
    "    session$setInputs(`summary-generate_report` = 1L)",
    "    session$flushReact()",
    "}",
    "docxText <- function(path) {",
    "    xml <- utils::unzip(path, 'word/document.xml', exdir = tempfile())",
    "    gsub('<[^>]+>', ' ', paste(readLines(xml, warn = FALSE), collapse = ' '))",
    "}",
    "sha <- function(x) digest::digest(serialize(x, NULL, version = 3L), algo = 'sha256', serialize = FALSE)",
    "flat <- function(x) gsub('[\\t\\r\\n]+', ' ', paste(as.character(x), collapse = ' | '))",
    "routes <- list(",
    "    `1` = list(gate = FALSE, managed = FALSE, summary = TRUE, norm = list(ruv_mode = 'automatic')),",
    "    `2` = list(gate = FALSE, managed = TRUE, summary = FALSE, norm = list(ruv_mode = 'automatic')),",
    "    `3` = list(gate = FALSE, managed = FALSE, summary = FALSE, norm = list(ruv_mode = 'manual', ruv_percentage = 10, ruv_k = 2)),",
    "    `4` = list(gate = TRUE, managed = FALSE, summary = FALSE, norm = list(ruv_mode = 'automatic',",
    "        ruv_design = 'groups_all_features', ruv_k_rule = 'user', ruv_k_user = 2)),",
    "    `5` = list(gate = TRUE, managed = FALSE, summary = FALSE, norm = list(ruv_mode = 'automatic',",
    "        ruv_design = 'groups_all_features', ruv_k_rule = 'num_sv_be')),",
    "    `6` = list(gate = TRUE, managed = FALSE, summary = FALSE, norm = list(ruv_mode = 'automatic',",
    "        ruv_design = 'first_pass_nulls', ruv_k_rule = 'user', ruv_k_user = 1)))",
    "rows <- list()",
    "failed <- FALSE",
    "for (route in names(routes)) {",
    "    r <- routes[[route]]",
    "    row <- NULL",
    "    ok <- tryCatch(testthat::test_that(paste('KIT-186 S1', side, 'route', route), {",
    "        withr::local_options(MultiScholaR.ruv_design_panel = r$gate)",
    "        fixture <- module_ci_trx_project(module_ci_trx_corpus())",
    "        if (isTRUE(r$summary)) {",
    "            assign('project_dirs', fixture$project_dirs, envir = globalenv())",
    "            withr::defer(if (exists('project_dirs', envir = globalenv())) rm('project_dirs', envir = globalenv()))",
    "        }",
    "        if (isTRUE(r$managed)) {",
    "            probe <- new.env()",
    "            omicsManagedChainMocks(fixture, probe)",
    "        }",
    "        journey <- module_ci_trx_tabs_journey(fixture, managed = r$managed, norm = r$norm,",
    "            summary = if (isTRUE(r$summary)) summarySteps)",
    "        owner <- fixture$owner",
    "        roots <- tempdir()",
    "        record <- buildWorkflowAnalysisRecord(owner)",
    "        steps <- Filter(function(step) identical(step$operation, 'trx_ruvs'), record$steps)",
    "        ruv <- if (length(steps)) steps[[length(steps)]]$outcomes else list()",
    "        design <- names(readWorkflowOmicsMetadata(transcriptomicsDeclarationBinding, owner)$design_matrix)",
    "        notes <- vapply(journey$notes, function(note) paste0('[', note$type %||% '', '] ', note$message), character(1))",
    "        report <- ''",
    "        if (isTRUE(r$summary)) {",
    "            docx <- list.files(fixture$paths$results_summary_dir, pattern = '^transcriptomics_report.*\\\\.docx$', full.names = TRUE)",
    "            text <- if (length(docx)) docxText(docx[[1L]]) else '<no report>'",
    "            text <- guiEntryMask(gsub('\\\\s+', ' ', text), roots)",
    "            writeLines(text, file.path(out, paste0('report-text-', side, '.txt')))",
    "            report <- digest::digest(text, algo = 'sha256', serialize = FALSE)",
    "        }",
    "        row <<- data.frame(side = side, route = route, errors = flat(journey$errors), head = owner$current_state %||% '',",
    "            k = as.character(ruv$k %||% NA), w_columns = paste(design[startsWith(design, 'W_')], collapse = ','),",
    "            table_sha256 = if (is.null(journey$da$da_genes_long)) '<no table>' else sha(journey$da$da_genes_long),",
    "            formula_used = record$analyses$differential_abundance$model$formula_used %||% '<none>',",
    "            model_text = flat(trxDaModelText(fixture$workflow)),",
    "            notifications = flat(guiEntryMask(paste(notes, collapse = ' || '), roots)),",
    "            report_sha256 = report, stringsAsFactors = FALSE)",
    "        testthat::succeed()",
    "    }), error = function(e) {",
    "        cat('route', route, 'error', class(e)[[1L]], conditionMessage(e), '\\n')",
    "        FALSE",
    "    })",
    "    if (is.null(row)) row <- data.frame(side = side, route = route, errors = '<route failed>', head = '', k = '',",
    "        w_columns = '', table_sha256 = '', formula_used = '', model_text = '', notifications = '', report_sha256 = '',",
    "        stringsAsFactors = FALSE)",
    "    cat('KIT-186 S1', side, 'route', route, if (isTRUE(ok) && !nzchar(row$errors)) 'ok' else 'FAILED', '\\n')",
    "    if (!isTRUE(ok) || identical(row$errors, '<route failed>')) failed <- TRUE",
    "    rows[[route]] <- row",
    "}",
    "utils::write.table(do.call(rbind, rows), file.path(out, paste0('ruv-routes-', side, '.tsv')), sep = '\\t', quote = FALSE,",
    "    row.names = FALSE)",
    "quit(status = if (failed) 1L else 0L)"), child)

runSide <- function(tree, side) {
    console <- suppressWarnings(system2(rscript, c(shQuote(child), shQuote(tree), side, shQuote(out)), stdout = TRUE,
        stderr = TRUE))
    cat(grep("^KIT-186 S1|^route ", console, value = TRUE), sep = "\n")
    cat(utils::tail(console, 15L), sep = "\n")
    status <- attr(console, "status")
    if (is.null(status)) 0L else status
}

base_tree <- file.path(tempdir(), "kit186-s1-base")
invisible(system2("git", c("worktree", "add", "--detach", shQuote(base_tree), base_commit), stdout = TRUE, stderr = TRUE))
stopifnot(dir.exists(file.path(base_tree, "tests", "testthat")))
base_status <- runSide(base_tree, "base")
head_status <- runSide(head_root, "head")

readSide <- function(side) {
    path <- file.path(out, paste0("ruv-routes-", side, ".tsv"))
    if (!file.exists(path)) return(NULL)
    utils::read.delim(path, colClasses = "character", quote = "", na.strings = character(), comment.char = "")
}
base <- readSide("base")
head <- readSide("head")
stopifnot(!is.null(base), !is.null(head))
utils::write.table(rbind(base, head), file.path(out, "ruv-routes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

fields <- c("errors", "head", "k", "w_columns", "table_sha256", "formula_used", "model_text", "notifications", "report_sha256")
# A2's fields (ticket section 8). `errors`, `head` and `w_columns` are recorded beside them: the base tree's journey helper
# records the grid summary that a declared design does not define as an error (the R1 instrument defect, fixed in the
# head tree's helper), so `errors` differs on routes 4 to 6 by that entry alone.
a2_fields <- c("table_sha256", "formula_used", "k", "model_text", "notifications", "report_sha256")
compare <- list()
for (route in base$route) {
    b <- base[base$route == route, , drop = FALSE]
    h <- head[head$route == route, , drop = FALSE]
    # A route whose k resolves to 0 (route 5 only can) is a k = 0 route: recorded, not part of A2.
    k_zero <- identical(b$k, "0") || identical(h$k, "0")
    for (field in fields) {
        compare[[length(compare) + 1L]] <- data.frame(route = route, field = field, base = b[[field]], head = h[[field]],
            equal = identical(b[[field]], h[[field]]), a2 = !k_zero && field %in% a2_fields, stringsAsFactors = FALSE)
    }
}
compare <- do.call(rbind, compare)
utils::write.table(compare, file.path(out, "ruv-routes-compare.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
a2 <- compare[compare$a2, , drop = FALSE]
cat("KIT-186 S1: A2 rows", nrow(a2), "equal", sum(a2$equal), "routes", paste(unique(a2$route), collapse = ","), "\n")
cat("KIT-186 S1: base status", base_status, "head status", head_status, "\n")
quit(status = if (base_status != 0L || head_status != 0L) 1L else 0L)
