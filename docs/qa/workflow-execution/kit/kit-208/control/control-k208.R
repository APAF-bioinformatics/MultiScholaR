# KIT-208 control (after stop rule 2 in R2, run 45370085901): what moved behind the stored digest of
# `da_output/da_proteins_groupB-groupA_long_annot.tsv` in test-omic-da-module.R's "(b) in process: prot_load".
#
# Usage (run-tests scripts job, from the package root, full git history):
#   Rscript docs/qa/workflow-execution/kit/kit-208/control/control-k208.R <base commit> <head commit> <out dir>
# Internal: Rscript --vanilla control-k208.R child <tree> <side dir>
#
# Three sides, one process each, in detached worktrees under $RUNNER_TEMP/k208: base (the tests-only commit), head1
# and head2 (the fix commit, two separate processes). Each child loads its tree as the test jobs do (pkgload::load_all
# with helpers; the helpers sourced again as testthat::test_file does, working directory tests/testthat, TESTTHAT set),
# sources test-omic-da-module.R into an environment whose `test_that` does nothing (so only its definitions are made),
# rebinds `guiEntryArtifacts` there to a wrapper that calls the helper's own and also writes the masked text of the
# DA table each time a step's files are captured, then runs omicDaW1Case("prot_load") inside testthat::test_that().
# The digest recorded is the one guiEntryArtifacts itself returns (the value the test compares with the stored one).
# Writes, under <out dir>: <side>/da-table-step<k>.tsv (masked text, as digested), <side>/digests.tsv, and in the
# parent comparison.tsv, columns.tsv and status.tsv. No file of the package is edited; nothing is re-frozen.

STORED_PROT_LOAD <- "180c330673a89c3d58f1e830639d9011ddfc16f96aa8ee5ebdcbfc489d85c7c0" # test-omic-da-module.R:3343
STORED_PROT_MAIN <- "fbe4b6a921c2abe20abdb68f3fae2a34f3a45e70062415d3da8acebed58e73c6" # :3464
TABLE_KEY <- "da_output/da_proteins_groupB-groupA_long_annot.tsv"
CHILD_TIMEOUT_SECONDS <- 1800

childMain <- function(tree, side_dir) {
    dir.create(side_dir, recursive = TRUE, showWarnings = FALSE)
    setwd(tree)
    suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = TRUE))
    suppressPackageStartupMessages(library(testthat))
    test_dir <- normalizePath("tests/testthat")
    helper_env <- new.env(parent = globalenv())
    test_env <- new.env(parent = helper_env)
    assign("test_that", function(...) invisible(NULL), envir = test_env)
    withr::with_dir(test_dir, withr::with_envvar(c(TESTTHAT = "true"), {
        suppressMessages(testthat::source_test_helpers(test_dir, env = helper_env))
        sys.source(file.path(test_dir, "test-omic-da-module.R"), envir = test_env)
    }))
    original <- get("guiEntryArtifacts", envir = helper_env, mode = "function")
    mask <- get("guiEntryMask", envir = helper_env, mode = "function")
    calls <- 0L
    assign("guiEntryArtifacts", function(base_dir, roots) {
        out <- original(base_dir, roots)
        calls <<- calls + 1L
        path <- file.path(base_dir, TABLE_KEY)
        if (file.exists(path)) {
            writeLines(mask(readLines(path, warn = FALSE), roots), file.path(side_dir, sprintf("da-table-step%d.tsv", calls)))
            cat(paste(calls, out[[TABLE_KEY]] %||% NA_character_, sep = "\t"), "\n", sep = "",
                file = file.path(side_dir, "digests.tsv"), append = TRUE)
        }
        out
    }, envir = test_env)
    writeLines(paste("capture", "sha256", sep = "\t"), file.path(side_dir, "digests.tsv"))
    case_fn <- get("omicDaW1Case", envir = test_env, mode = "function")
    ok <- FALSE
    withr::with_dir(test_dir, withr::with_envvar(c(TESTTHAT = "true"), {
        testthat::test_that("KIT-208 control prot_load", {
            case_fn("prot_load")
            ok <<- TRUE
            testthat::succeed()
        })
    }))
    if (!ok) stop("the prot_load case did not complete")
    0L
}

readTable <- function(path) utils::read.delim(path, sep = "\t", quote = "", check.names = FALSE,
    colClasses = "character", na.strings = character(), stringsAsFactors = FALSE)

lastTable <- function(side_dir) {
    digests <- utils::read.delim(file.path(side_dir, "digests.tsv"), colClasses = "character")
    last <- digests[nrow(digests), , drop = FALSE]
    list(path = file.path(side_dir, sprintf("da-table-step%s.tsv", last$capture)), digests = digests)
}

parentMain <- function(base, head, out) {
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    out <- normalizePath(out)
    script <- normalizePath("docs/qa/workflow-execution/kit/kit-208/control/control-k208.R")
    work <- file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "k208")
    dir.create(work, recursive = TRUE, showWarnings = FALSE)
    status <- file.path(out, "status.tsv")
    writeLines(paste("side", "commit", "exit", "seconds", sep = "\t"), status)
    sides <- c(base = base, head1 = head, head2 = head)
    for (side in names(sides)) {
        tree <- file.path(work, side)
        stopifnot(system2("git", c("worktree", "add", "--detach", tree, sides[[side]])) == 0L)
        started <- proc.time()[["elapsed"]]
        code <- system2("timeout", c(CHILD_TIMEOUT_SECONDS, "Rscript", "--vanilla", script, "child", tree,
            file.path(out, side)), stdout = file.path(work, paste0(side, ".log")), stderr = file.path(work, paste0(side, ".log")))
        cat(paste(side, sides[[side]], code, round(proc.time()[["elapsed"]] - started, 1), sep = "\t"), "\n",
            sep = "", file = status, append = TRUE)
        if (code != 0L) {
            writeLines(utils::tail(readLines(file.path(work, paste0(side, ".log")), warn = FALSE), 40L),
                file.path(out, paste0(side, "-log-tail.txt")))
            stop(sprintf("side %s exited %s", side, code))
        }
    }
    b <- lastTable(file.path(out, "base"))
    h1 <- lastTable(file.path(out, "head1"))
    h2 <- lastTable(file.path(out, "head2"))
    tb <- readTable(b$path)
    th <- readTable(h1$path)
    common <- intersect(names(tb), names(th))
    other <- setdiff(common, c("gene_name", "gene_names"))
    same_rows <- nrow(tb) == nrow(th) && identical(tb[[1L]], th[[1L]])
    blank <- function(x) if (is.null(x)) NA_integer_ else sum(!nzchar(x) | x == "NA")
    differing <- other[!vapply(other, function(col) identical(tb[[col]], th[[col]]), NA)]
    rows <- list(
        c("precondition: base digest equals the stored prot_load digest (:3343)", b$digests$sha256[nrow(b$digests)],
            STORED_PROT_LOAD, identical(b$digests$sha256[nrow(b$digests)], STORED_PROT_LOAD)),
        c("head digest", h1$digests$sha256[nrow(h1$digests)], "", NA),
        c("head digest equals the stored prot_main digest (:3464)", h1$digests$sha256[nrow(h1$digests)],
            STORED_PROT_MAIN, identical(h1$digests$sha256[nrow(h1$digests)], STORED_PROT_MAIN)),
        c("head1 and head2 dumps byte-identical", unname(tools::md5sum(h1$path)), unname(tools::md5sum(h2$path)),
            identical(unname(tools::md5sum(h1$path)), unname(tools::md5sum(h2$path)))),
        c("head1 and head2 digest sequences identical", paste(h1$digests$sha256, collapse = ";"),
            paste(h2$digests$sha256, collapse = ";"), identical(h1$digests$sha256, h2$digests$sha256)),
        c("row count", nrow(tb), nrow(th), nrow(tb) == nrow(th)),
        c("row order (first column, in order)", names(tb)[[1L]], names(th)[[1L]], same_rows),
        c("columns only on base", paste(setdiff(names(tb), names(th)), collapse = ";"), "", NA),
        c("columns only on head", "", paste(setdiff(names(th), names(tb)), collapse = ";"), NA),
        c("common columns other than gene_name and gene_names", length(other), length(other), NA),
        c("of those, byte-identical", length(other) - length(differing), length(other), length(differing) == 0L),
        c("of those, differing", paste(differing, collapse = ";"), "", length(differing) == 0L),
        c("gene_name blank rows", blank(tb$gene_name), blank(th$gene_name), NA),
        c("gene_names column count", sum(names(tb) == "gene_names"), sum(names(th) == "gene_names"), NA),
        c("head gene_name equals first word of head gene_names", "", "",
            if (is.null(th$gene_names)) NA else identical(th$gene_name,
                vapply(strsplit(th$gene_names, " |:"), function(x) if (length(x)) x[[1L]] else "", ""))))
    comparison <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
    names(comparison) <- c("check", "base", "head", "holds")
    utils::write.table(comparison, file.path(out, "comparison.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    cols <- union(names(tb), names(th))
    utils::write.table(data.frame(column = cols, on_base = cols %in% names(tb), on_head = cols %in% names(th),
        identical = vapply(cols, function(col) if (col %in% common) identical(tb[[col]], th[[col]]) else NA, NA)),
        file.path(out, "columns.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    print(comparison)
    invisible(0L)
}

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) >= 1L && identical(arguments[[1L]], "child")) {
    quit(status = childMain(arguments[[2L]], arguments[[3L]]))
}
if (length(arguments) != 3L) stop("usage: control-k208.R <base commit> <head commit> <out dir>")
parentMain(arguments[[1L]], arguments[[2L]], arguments[[3L]])
