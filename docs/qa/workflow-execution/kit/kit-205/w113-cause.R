# KIT-205 W113 cause control (ticket 6.3, A7, stop rule 5; KIT-187's provenance-cause.R form, VAL-104). Evidence code
# under docs/qa: no package code. Runs on a GitHub-hosted runner (run-tests `scripts`) from the package root of the
# commit under test (head), after head is installed.
#
# Usage: Rscript docs/qa/workflow-execution/kit/kit-205/w113-cause.R <out dir> <base commit>
#
# Base (the base commit's tree, git archive) and head each run, in a fresh Rscript with pkgload::load_all(<tree>,
# export_all = TRUE, helpers = TRUE), the one test "(KIT-113 W1) each W1 (b) session's ..." of a copy of that tree's
# tests/testthat/test-omic-enrich-module.R. The copy adds one step after each case's D6 capture that keeps the step
# provenance of `enrichment_completed` and the analysis record's enrichment step unhashed (committed_at, build and
# source_commit removed; store-local ids and the fixture's root masked as the capture masks them), stored by case. The
# copy is removed afterwards. The two sides are then compared leaf by leaf (each leaf by the SHA-256 of its
# serialisation). Declared expectation: the provenance differs only by `fixed/enrichment_background` = "tested" on
# head; the record step only by one parameter entry on head (name `enrichment_background`, source "fixed", provenance
# "fixed", value "tested"). Anything else is reported.
# Writes <out dir>/w113-cause.tsv (case, field, leaves_equal, marker_leaves, other_leaves), w113-detail.tsv (every
# leaf that differs), w113-verdict.txt and the side logs. Exit status 0 (a difference is a result).

DESC <- "(KIT-113 W1) each W1 (b) session's committed record, results, product, provenance and analysis record step"
MARKER_PATH <- "fixed/enrichment_background"
CASES_WITH_STATE <- c("gprofiler", "clusterprofiler", "kernel_gprofiler", "kernel_clusterprofiler",
    "organism_filter/no_mapping", "organism_filter/mapped", "restore")

substituteOnce <- function(text, old, new) {
    hits <- gregexpr(old, text, fixed = TRUE)[[1L]]
    if (length(hits) != 1L || hits[[1L]] < 0L) stop("anchor not found exactly once: ", substr(old, 1L, 80L))
    paste0(substr(text, 1L, hits[[1L]] - 1L), new, substr(text, hits[[1L]] + nchar(old), nchar(text)))
}

RAW <- paste(
    "k205Raw <- function(ctx) {",
    "    manager <- shiny::isolate(ctx$bus()$state_manager)",
    "    snapshot <- withWorkflowStateScope(manager, function() manager$.context$snapshot)",
    "    provenance <- snapshot$step_provenance[[\"enrichment_completed\"]]",
    "    if (!is.null(provenance)) provenance <- omicEnrichW113Mask(provenance[setdiff(names(provenance),",
    "        c(\"committed_at\", \"build\", \"source_commit\"))], ctx$roots)",
    "    steps <- tryCatch(buildWorkflowAnalysisRecord(manager)$steps, error = function(e) NULL)",
    "    step <- Filter(function(s) identical(s$operation, \"enrichment\"), steps %||% list())",
    "    record <- if (length(step)) omicEnrichW113Mask(step[[length(step)]][c(\"state\", \"operation\", \"variant\",",
    "        \"step_label\", \"description\", \"parameters\")], ctx$roots) else NULL",
    "    list(provenance = provenance, record = record)",
    "}", "", sep = "\n")

runSide <- function(tree, side, work) {
    setwd(tree)
    suppressPackageStartupMessages(pkgload::load_all(".", export_all = TRUE, helpers = TRUE, quiet = TRUE))
    source_file <- file.path("tests", "testthat", "test-omic-enrich-module.R")
    text <- paste(readLines(source_file, warn = FALSE), collapse = "\n")
    text <- substituteOnce(text, 'list(kit113_d6 = function(ctx) ctx$keep("d6", omicEnrichW113Fields(ctx)))), live = spec$live)',
        'list(kit113_d6 = function(ctx) ctx$keep("d6", omicEnrichW113Fields(ctx)), k205_raw = function(ctx) ctx$keep("raw", k205Raw(ctx)))), live = spec$live)')
    text <- substituteOnce(text, "            capture <- seen$kit113_d6$d6\n",
        "            capture <- seen$kit113_d6$d6\n            assign(key, seen$k205_raw$raw, envir = get(\"k205_store\", envir = globalenv()))\n")
    text <- substituteOnce(text, paste0('test_that("', DESC, '", {'), paste0(RAW, 'test_that("', DESC, '", {'))
    copy <- file.path("tests", "testthat", "test-k205-w113-dump.R")
    writeLines(text, copy)
    on.exit(unlink(copy), add = TRUE)
    assign("k205_store", new.env(), envir = globalenv())
    result <- testthat::test_file(copy, desc = DESC, reporter = "summary", stop_on_failure = FALSE)
    cat(sprintf("side %s: kept %d cases\n", side, length(ls(k205_store))))
    print(as.data.frame(result)[, c("test", "nb", "failed", "error")])
    saveRDS(as.list(k205_store), file.path(work, paste0(side, ".rds")))
}

# Leaves of a nested value as path -> SHA-256 of the serialised leaf (a list's unnamed elements by position).
flatten <- function(x, path = character()) {
    if (is.list(x) && !is.object(x) && length(x)) {
        keys <- names(x) %||% rep("", length(x))
        keys[!nzchar(keys)] <- paste0("[", which(!nzchar(keys)), "]")
        out <- character()
        for (i in seq_along(x)) out <- c(out, flatten(x[[i]], c(path, keys[[i]])))
        return(out)
    }
    stats::setNames(digest::digest(x, algo = "sha256"), paste(path, collapse = "/"))
}

args <- commandArgs(TRUE)
if (length(args) && identical(args[[1L]], "--side")) {
    runSide(args[[2L]], args[[3L]], args[[4L]])
    quit(status = 0L)
}
stopifnot(length(args) == 2L)
dir.create(args[[1L]], showWarnings = FALSE, recursive = TRUE)
out <- normalizePath(args[[1L]], mustWork = TRUE)
base <- args[[2L]]
head_root <- normalizePath(".")
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]))
work <- file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "k205-w113")
base_root <- file.path(work, "base-src")
dir.create(base_root, recursive = TRUE, showWarnings = FALSE)
archive_exit <- system(sprintf("git archive %s | tar -x -C %s", shQuote(base), shQuote(base_root)))
rscript <- file.path(R.home("bin"), "Rscript")
exits <- c()
for (side in c("base", "head")) {
    tree <- if (side == "base") base_root else head_root
    exits[[side]] <- system2(rscript, c(shQuote(script), "--side", shQuote(tree), side, shQuote(work)),
        stdout = file.path(out, paste0("side-", side, ".log.txt")), stderr = file.path(out, paste0("side-", side, ".log.txt")),
        env = "NOT_CRAN=true")
}
if (!all(file.exists(file.path(work, c("base.rds", "head.rds"))))) {
    writeLines(c(sprintf("archive %d; side exits %s", archive_exit, paste(names(exits), exits, collapse = ", ")),
        "VERDICT: INSTRUMENT (a side wrote no result)"), file.path(out, "w113-verdict.txt"))
    quit(status = 0L)
}
b <- readRDS(file.path(work, "base.rds"))
h <- readRDS(file.path(work, "head.rds"))
suppressPackageStartupMessages(library(MultiScholaR))
marker_value <- MultiScholaR:::encodeWorkflowRecordValue("tested")
marker_sha <- digest::digest("tested", algo = "sha256")
null_sha <- digest::digest(NULL, algo = "sha256")
isMarkerEntry <- function(p) is.list(p) && identical(p$name, "enrichment_background") && identical(p$source, "fixed") &&
    identical(p$provenance, "fixed") && identical(p$value, marker_value)

summary_rows <- list()
detail_rows <- list()
detail <- function(case, field, path, bv, hv, status) detail_rows[[length(detail_rows) + 1L]] <<- data.frame(case = case,
    field = field, path = path, base = substr(bv, 1L, 16L), head = substr(hv, 1L, 16L), status = status,
    stringsAsFactors = FALSE)
cases <- sort(union(names(b), names(h)))
for (case in cases) {
    # provenance, leaf by leaf
    bl <- if (is.null(b[[case]]$provenance)) character() else flatten(b[[case]]$provenance)
    hl <- if (is.null(h[[case]]$provenance)) character() else flatten(h[[case]]$provenance)
    equal <- 0L
    marker <- 0L
    other <- 0L
    for (path in sort(union(names(bl), names(hl)))) {
        bv <- if (path %in% names(bl)) bl[[path]] else NA_character_
        hv <- if (path %in% names(hl)) hl[[path]] else NA_character_
        # A base entry may hold `fixed` as NULL where head holds the marker list.
        if (identical(bv, hv)) equal <- equal + 1L else if ((is.na(bv) && identical(path, MARKER_PATH) &&
            identical(hv, marker_sha)) || (identical(path, "fixed") && identical(bv, null_sha) && is.na(hv))) {
            marker <- marker + 1L
            detail(case, "provenance", path, bv, hv, "marker_on_head")
        } else {
            other <- other + 1L
            detail(case, "provenance", path, bv, hv, "UNEXPLAINED")
        }
    }
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(case = case, field = "provenance", leaves_equal = equal,
        marker_leaves = marker, other_leaves = other, stringsAsFactors = FALSE)
    # record: head with its one marker entry taken out must equal base leaf by leaf
    head_record <- h[[case]]$record
    marker <- 0L
    if (is.list(head_record) && length(head_record$parameters)) {
        hit <- vapply(head_record$parameters, isMarkerEntry, logical(1))
        if (sum(hit) == 1L) {
            marker <- length(flatten(head_record$parameters[hit]))
            head_record$parameters <- head_record$parameters[!hit]
            detail(case, "record", "parameters/<enrichment_background entry>", NA_character_, "fixed;fixed;tested",
                "marker_on_head")
        }
    }
    bl <- if (is.null(b[[case]]$record)) character() else flatten(b[[case]]$record)
    hl <- if (is.null(head_record)) character() else flatten(head_record)
    equal <- 0L
    other <- 0L
    for (path in sort(union(names(bl), names(hl)))) {
        bv <- if (path %in% names(bl)) bl[[path]] else NA_character_
        hv <- if (path %in% names(hl)) hl[[path]] else NA_character_
        if (identical(bv, hv)) equal <- equal + 1L else {
            other <- other + 1L
            detail(case, "record", path, bv, hv, "UNEXPLAINED")
        }
    }
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(case = case, field = "record", leaves_equal = equal,
        marker_leaves = marker, other_leaves = other, stringsAsFactors = FALSE)
}
summary <- do.call(rbind, summary_rows)
details <- if (length(detail_rows)) do.call(rbind, detail_rows) else data.frame(case = character(), field = character(),
    path = character(), base = character(), head = character(), status = character())
utils::write.table(summary, file.path(out, "w113-cause.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(details, file.path(out, "w113-detail.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
with_state <- summary[summary$case %in% CASES_WITH_STATE, , drop = FALSE]
missing <- setdiff(CASES_WITH_STATE, intersect(names(b), names(h)))
ok <- !sum(summary$other_leaves) && !length(missing) && nrow(with_state) == 2L * length(CASES_WITH_STATE) &&
    all(with_state$marker_leaves > 0L)
verdict <- c(sprintf("base %s, head %s; archive %d; side exits %s", base, system2("git", c("rev-parse", "HEAD"),
        stdout = TRUE), archive_exit, paste(names(exits), exits, collapse = ", ")),
    sprintf("cases: base %d, head %d; cases with a completed state missing: %s", length(b), length(h),
        if (length(missing)) paste(missing, collapse = ", ") else "none"),
    sprintf("provenance: %d leaves equal, %d marker leaves, %d other", sum(summary$leaves_equal[summary$field == "provenance"]),
        sum(summary$marker_leaves[summary$field == "provenance"]), sum(summary$other_leaves[summary$field == "provenance"])),
    sprintf("record: %d leaves equal (head without its marker entry), %d marker leaves, %d other",
        sum(summary$leaves_equal[summary$field == "record"]), sum(summary$marker_leaves[summary$field == "record"]),
        sum(summary$other_leaves[summary$field == "record"])),
    if (ok) "VERDICT: only the marker differs (provenance and record)" else
        "VERDICT: STOP (an unexplained difference, a missing case or no marker)")
writeLines(verdict, file.path(out, "w113-verdict.txt"))
cat(verdict, sep = "\n")
print(summary)
