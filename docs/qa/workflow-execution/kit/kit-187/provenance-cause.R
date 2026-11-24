# KIT-187 (VAL-104): the cause of the moved KIT-111 provenance and record digests in tests/testthat/test-prot-norm-core.R.
# Base and head each run the three KIT-111 tests whose values moved (W1 (a), W2 (W1 (b)), W2 (R7)) and keep each
# case's `provenance` and `record` fields unhashed, as the test captures them (`protNormW1Provenance`: committed_at,
# build and source_commit removed; `protNormW1RecordSteps`: the analysis record's steps without ids, times and
# software; the project root masked). The two sides are then compared leaf by leaf. Declared expectation: provenance
# differs only by the field `fixed$missing_value_comparison` (value "counts") present on head; the record differs only
# by one parameter entry per normalisation step (`operation` "normalization-ruv") present on head, name
# `missing_value_comparison`, source "fixed", provenance "fixed", value "counts". Anything else is reported and the
# script exits with status 1. (Revision 2 adds the record half; the provenance comparison is unchanged.)
#
# Usage, from the package root of the bound commit on a run-tests runner (full git history):
#   Rscript docs/qa/workflow-execution/kit/kit-187/provenance-cause.R <base commit> <out dir>
# Writes <out dir>/provenance-compare.tsv (case, state, path, base, head, status) and provenance-verdict.txt.
# Each side runs in a fresh Rscript with pkgload::load_all(<tree>, export_all = TRUE, helpers = TRUE); the test file
# is copied to a temporary test file in that side's tests/testthat with two substitutions (the capture keeps the
# provenance instead of its digest; the two expectation helpers store it under the case's name) and removed after.

DESCS <- c(w1 = "KIT-111 W1 (a): the DIA normalisation and correlation clicks give the recorded states, products, provenance, files and display",
    w2 = "KIT-111 W2 (W1 (b)): the DIA normalisation and correlation clicks through the worker give W1 (a)'s values",
    r7 = "KIT-111 W2 (R7): a managed publication that fails is completed by the repeated click, a hit that commits nothing")
MARKER_PATH <- "fixed/missing_value_comparison"

substituteOnce <- function(text, old, new) {
    hits <- gregexpr(old, text, fixed = TRUE)[[1L]]
    if (length(hits) < 1L || hits[[1L]] < 0L) stop("anchor not found: ", old)
    at <- hits[[1L]]
    paste0(substr(text, 1L, at - 1L), new, substr(text, at + nchar(old), nchar(text)))
}

runSide <- function(tree, side, work) {
    setwd(tree)
    suppressPackageStartupMessages(pkgload::load_all(".", export_all = TRUE, helpers = TRUE, quiet = TRUE))
    source_file <- file.path("tests", "testthat", "test-prot-norm-core.R")
    text <- paste(readLines(source_file, warn = FALSE), collapse = "\n")
    text <- substituteOnce(text, "provenance = protNormW1Sha(protNormW1Provenance(snapshot, root)),",
        "provenance = protNormW1Provenance(snapshot, root),")
    text <- substituteOnce(text, "record = protNormW1Sha(protNormW1RecordSteps(manager, root)),",
        "record = protNormW1RecordSteps(manager, root),")
    text <- substituteOnce(text, "    frozen <- protNormW1Frozen()[[case]]\n    expect_false(is.null(frozen), info = case)\n    for (field in union(",
        "    assign(paste0(get('k187_prefix', envir = globalenv()), '/', case), capture[c('provenance', 'record')], envir = get('k187_store', envir = globalenv()))\n    frozen <- protNormW1Frozen()[[case]]\n    expect_false(is.null(frozen), info = case)\n    for (field in union(")
    text <- substituteOnce(text, "    frozen <- protNormW1Frozen()[[case]]\n    expect_false(is.null(frozen), info = case)\n    for (field in setdiff(",
        "    assign(paste0(get('k187_prefix', envir = globalenv()), '/', case), capture[c('provenance', 'record')], envir = get('k187_store', envir = globalenv()))\n    frozen <- protNormW1Frozen()[[case]]\n    expect_false(is.null(frozen), info = case)\n    for (field in setdiff(")
    copy <- file.path("tests", "testthat", "test-k187-provenance-dump.R")
    writeLines(text, copy)
    on.exit(unlink(copy), add = TRUE)
    assign("k187_store", new.env(), envir = globalenv())
    for (prefix in names(DESCS)) {
        assign("k187_prefix", prefix, envir = globalenv())
        invisible(utils::capture.output(suppressMessages(testthat::test_file(copy, desc = DESCS[[prefix]],
            reporter = "silent", stop_on_failure = FALSE))))
        cat(sprintf("side %s: %s kept %d cases so far\n", side, prefix, length(ls(k187_store))))
    }
    saveRDS(as.list(k187_store), file.path(work, paste0(side, ".rds")))
}

# Leaves of a nested value as path -> SHA-256 of the serialised leaf (a list's unnamed elements by position).
flatten <- function(x, path = character()) {
    if (is.list(x) && length(x)) {
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
stopifnot(length(args) == 2L, grepl("^[0-9a-f]{40}$", args[[1L]]))
out <- args[[2L]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)
head_root <- normalizePath(".")
script_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1L]]))
base_root <- file.path(tempdir(), "k187-prov-base")
invisible(system2("git", c("worktree", "add", "--detach", shQuote(base_root), args[[1L]]), stdout = TRUE, stderr = TRUE))
stopifnot(dir.exists(file.path(base_root, "R")))
work <- file.path(tempdir(), "k187-prov-work")
dir.create(work, showWarnings = FALSE)
rscript <- file.path(R.home("bin"), "Rscript")
for (side in c("base", "head")) {
    tree <- if (side == "base") base_root else head_root
    status <- system2(rscript, c(shQuote(script_file), "--side", shQuote(tree), side, shQuote(work)))
    if (!identical(as.integer(status), 0L)) stop("side ", side, " failed with status ", status)
}
base <- readRDS(file.path(work, "base.rds"))
head <- readRDS(file.path(work, "head.rds"))
marker_sha <- digest::digest("counts", algo = "sha256")
# The record half: on head, take out each normalisation step's marker entry of the declared form, counting its leaves.
isMarkerEntry <- function(p) is.list(p) && identical(p$name, "missing_value_comparison") && identical(p$source, "fixed") &&
    identical(p$provenance, "fixed") && identical(decodeValue(p$value), "counts")
decodeValue <- function(v) if (is.list(v) && identical(v$type, "character")) unlist(v$value) else v
stripRecordMarkers <- function(steps) {
    removed <- 0L
    leaves <- 0L
    if (!is.list(steps)) return(list(steps = steps, removed = 0L, leaves = 0L))
    for (i in seq_along(steps)) {
        step <- steps[[i]]
        if (!is.list(step) || !identical(step$operation, "normalization-ruv") || !length(step$parameters)) next
        hit <- vapply(step$parameters, isMarkerEntry, logical(1L))
        if (sum(hit) > 1L) next
        if (any(hit)) {
            leaves <- leaves + length(flatten(step$parameters[hit]))
            step$parameters <- step$parameters[!hit]
            removed <- removed + 1L
            steps[[i]] <- step
        }
    }
    list(steps = steps, removed = removed, leaves = leaves)
}
rows <- list()
unexplained <- c(provenance = 0L, record = 0L)
markers <- c(provenance = 0L, record = 0L)
record_same <- 0L
cases <- sort(union(names(base), names(head)))
for (case in cases) {
    # provenance, leaf by leaf (as revision 1)
    b <- if (is.null(base[[case]]$provenance)) character() else flatten(base[[case]]$provenance)
    h <- if (is.null(head[[case]]$provenance)) character() else flatten(head[[case]]$provenance)
    for (path in sort(union(names(b), names(h)))) {
        bv <- if (path %in% names(b)) b[[path]] else NA_character_
        hv <- if (path %in% names(h)) h[[path]] else NA_character_
        status <- if (identical(bv, hv)) "same" else if (is.na(bv) && endsWith(path, MARKER_PATH) &&
            identical(hv, marker_sha)) "marker_on_head" else "UNEXPLAINED"
        if (identical(status, "UNEXPLAINED")) unexplained[["provenance"]] <- unexplained[["provenance"]] + 1L
        if (identical(status, "marker_on_head")) markers[["provenance"]] <- markers[["provenance"]] + 1L
        rows[[length(rows) + 1L]] <- data.frame(field = "provenance", case = case, path = path, base = substr(bv, 1L, 16L),
            head = substr(hv, 1L, 16L), status = status, stringsAsFactors = FALSE)
    }
    # record: head with its marker entries taken out must equal base leaf by leaf
    stripped <- stripRecordMarkers(head[[case]]$record)
    markers[["record"]] <- markers[["record"]] + stripped$leaves
    rows[[length(rows) + 1L]] <- data.frame(field = "record", case = case, path = "<marker entries removed on head>",
        base = NA_character_, head = as.character(stripped$removed), status = sprintf("marker_entries_%d_leaves_%d",
            stripped$removed, stripped$leaves), stringsAsFactors = FALSE)
    if (!stripped$removed) unexplained[["record"]] <- unexplained[["record"]] + 1L
    b <- if (is.null(base[[case]]$record)) character() else flatten(base[[case]]$record)
    h <- if (is.null(stripped$steps)) character() else flatten(stripped$steps)
    for (path in sort(union(names(b), names(h)))) {
        bv <- if (path %in% names(b)) b[[path]] else NA_character_
        hv <- if (path %in% names(h)) h[[path]] else NA_character_
        status <- if (identical(bv, hv)) "same" else "UNEXPLAINED"
        if (identical(status, "same")) record_same <- record_same + 1L else {
            unexplained[["record"]] <- unexplained[["record"]] + 1L
            rows[[length(rows) + 1L]] <- data.frame(field = "record", case = case, path = path, base = substr(bv, 1L, 16L),
                head = substr(hv, 1L, 16L), status = status, stringsAsFactors = FALSE)
        }
    }
}
table <- do.call(rbind, rows)
utils::write.table(table, file.path(out, "provenance-compare.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
expected_cases <- as.vector(outer(c("w1", "w2"), c("dia_normalization_skip", "dia_normalization_skip_repeat",
    "dia_normalization_manual", "dia_correlation_apply", "dia_correlation_reapply", "dia_correlation_skip",
    "dia_normalization_refused", "dia_normalization_retry"), paste, sep = "/"))
missing_cases <- setdiff(c(expected_cases, "r7/dia_normalization_skip"), intersect(names(base), names(head)))
ok <- !sum(unexplained) && !length(missing_cases) && all(markers > 0L)
verdict <- c(sprintf("base %s head %s", args[[1L]], system2("git", c("rev-parse", "HEAD"), stdout = TRUE)),
    sprintf("cases compared: %d (base %d, head %d); missing expected cases: %s", length(cases), length(base), length(head),
        if (length(missing_cases)) paste(missing_cases, collapse = ", ") else "none"),
    sprintf("provenance leaves: %d same, %d marker_on_head, %d UNEXPLAINED", sum(table$field == "provenance" &
        table$status == "same"), markers[["provenance"]], unexplained[["provenance"]]),
    sprintf("record leaves: %d same (head after removing its marker entries), %d in marker entries on head, %d UNEXPLAINED (a case with no marker entry counts as one)",
        record_same, markers[["record"]], unexplained[["record"]]),
    if (ok) "VERDICT: only the marker differs (provenance and record)" else
        "VERDICT: STOP (an unexplained difference, a missing case or no marker)")
writeLines(verdict, file.path(out, "provenance-verdict.txt"))
cat(verdict, sep = "\n")
if (!ok) quit(status = 1L)

