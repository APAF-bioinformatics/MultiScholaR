# The licence guard: an allow-list that covers logs (KIT-159 section 5.2; VAL-47).

#' Columns of the fixed evidence tables (section 5.7)
#' @noRd
EVIDENCE_COLUMNS <- list(
    cells.tsv = c("unit", "run", "arm", "exit", "stop_reason", "step", "class", "message_sha256",
                  "fallback", "options_at_start", "ruv_k", "ruv_controls", "seconds", "chunk"),
    calls.tsv = c("unit", "arm", "contrast", "call_statistic", "fallback_count"),
    metrics.tsv = c("ticket", "dataset", "truth_id", "unit", "contrast", "arm", "assay", "view",
                    "class", "call_rule", "metric", "statistic_class", "value", "method", "form",
                    "scope", "ci_low", "ci_high", "achieved_level", "n_units", "unit_kind", "B",
                    "seed", "na_share", "status", "primary", "family", "blinded"),
    comparisons.tsv = c("arm", "comparator", "truth_id", "metric", "view", "call_rule", "diff",
                        "method", "scope", "ci_low", "ci_high", "margin", "conf", "outcome",
                        "no_worse", "outcome_family", "n_refused", "regime", "flags", "member",
                        "component", "not_sensitive", "amended", "headline", "blinded"),
    coverage.tsv = c("truth_id", "unit", "arm", "contrast", "class", "keys", "quantified",
                     "testable", "called_own", "called_bh", "completeness"),
    attribution.tsv = c("role", "kind", "licence", "redistribution", "citation"),
    determinism.tsv = c("table", "sha256_job1", "sha256_job2", "identical", "sha256_10_job1",
                        "sha256_10_job2", "max_abs_diff", "max_rel_diff", "cause"),
    ledger = c("step", "action", "input_id", "value", "source", "location")
)

#' File names admitted without a column schema (their content is fixed by the harness)
#' @noRd
EVIDENCE_FREE_NAMES <- c(
    "provenance.json", "results.md", "tables.sha256", "not-committed.sha256", "outputs.sha256",
    "inputs.tsv", "library.tsv", "coverage-table.tsv", "coverage-detail.tsv", "coverage-run.txt",
    "ticket-coverage.tsv", "harness-tests.tsv", "harness-failures.txt", "ticket-tests.tsv",
    "ticket-failures.txt", "code-hashes.tsv", "fsim.tsv", "cost-probe.tsv"
)
RESULT_TYPES <- c("label", "count", "share", "quantile", "sd", "n", "abs_diff", "hash", "statistic")
CALL_COUNT_METRICS <- c("TP", "FP", "calls", "keys_called", "called")
MAX_FILE_BYTES <- 5 * 1024^2
ARCHIVE_PATTERN <- "\\.(zip|gz|tgz|tar|bz2|xz|7z|rar|rds|rdata|so|o|a|dll)$"

#' Is every `kind = data` input in a lineage redistributable?
#'
#' A program's licence never makes an output private (VAL-47).
#'
#' @param manifest The input manifest.
#' @param lineage Roles in the output's lineage.
#' @return `TRUE` when every data role in the lineage is `permitted`.
#' @export
valPublishable <- function(manifest, lineage) {
    if (!all(lineage %in% manifest$role)) {
        valAbort("a lineage names an undeclared role", "val_spec_invalid")
    }
    rows <- manifest[manifest$role %in% lineage & manifest$kind == "data", , drop = FALSE]
    all(rows$redistribution == "permitted")
}

#' Must logs be private (any `kind = data` input with `redistribution = not_permitted`)?
#'
#' @param manifest The input manifest.
#' @return `TRUE` or `FALSE`.
#' @export
valPrivateLogs <- function(manifest) {
    any(manifest$kind == "data" & manifest$redistribution == "not_permitted")
}

#' Print one fixed-template status line (the only console text under a private log)
#'
#' Values are reduced to `[A-Za-z0-9._:-]`, so no data text can pass. With `VAL_TEMPLATE_FILE`
#' set (a child of [valRunAsChild()]) the line goes to that file.
#'
#' @param phase,dataset,units,arm,status,class,count Template fields.
#' @return The line, invisibly.
#' @export
valStatus <- function(phase = NULL, dataset = NULL, units = NULL, arm = NULL, status = NULL,
                      class = NULL, count = NULL) {
    fields <- list(phase = phase, dataset = dataset, units = units, arm = arm, status = status,
                   class = class, count = count)
    fields <- fields[!vapply(fields, is.null, logical(1))]
    clean <- vapply(fields, function(v) gsub("[^A-Za-z0-9._:-]", "?", as.character(v)[[1L]]),
                    character(1))
    line <- paste("[msrval]", paste0(names(clean), "=", clean, collapse = " "))
    target <- Sys.getenv("VAL_TEMPLATE_FILE")
    if (nzchar(target)) {
        cat(line, "\n", file = target, sep = "", append = TRUE)
    } else {
        cat(line, "\n", sep = "")
    }
    invisible(line)
}

#' Re-run a script as a child whose stdout and stderr go to private files (A27)
#'
#' `sink()` does not divert text written to the file descriptors by compiled code or by a tool
#' started with `system2`; a child process with both descriptors on files does. Only the template
#' lines the child writes through [valStatus()] are returned and printed.
#'
#' @param script R script.
#' @param args Its arguments.
#' @param private_dir `$RUNNER_TEMP/val-private`.
#' @return The template lines, with attribute `status` (the child's exit status).
#' @export
valRunAsChild <- function(script, args, private_dir) {
    logs <- file.path(private_dir, "logs")
    dir.create(logs, recursive = TRUE, showWarnings = FALSE)
    tpl <- file.path(private_dir, "template.txt")
    file.create(tpl)
    # The shell opens both private files in append mode, so a child that reopens /dev/stderr
    # (a fresh descriptor at offset 0) cannot make a later write overwrite earlier text.
    res <- processx::run("sh", c("-c", 'exec "$0" "$@" >> "$VAL_OUT" 2>> "$VAL_ERR"',
                                 file.path(R.home("bin"), "Rscript"), "--vanilla", script, args),
                         env = c("current", VAL_TEMPLATE_FILE = tpl,
                                 VAL_OUT = file.path(logs, "child.out"),
                                 VAL_ERR = file.path(logs, "child.err")),
                         error_on_status = FALSE)
    lines <- readLines(tpl, warn = FALSE)
    if (length(lines) > 0L) cat(lines, sep = "\n")
    structure(lines, status = res$status)
}

#' Expand a declared output path's placeholders (section 5.2, C77)
#'
#' `{unit}` and `{arm}` are built in; any other name must be declared in `placeholders`; a declared
#' placeholder the path does not use, or a value outside `^[A-Za-z0-9._-]+$`, is refused.
#'
#' @param output One declared output (`path`, `placeholders`).
#' @param arms The declaration's arm ids.
#' @param units The data set's unit ids.
#' @return The expanded paths.
#' @export
valExpandOutputPaths <- function(output, arms = character(), units = character()) {
    path <- output$path
    used <- unique(regmatches(path, gregexpr("\\{[A-Za-z0-9_]+\\}", path))[[1L]])
    used <- gsub("[{}]", "", used)
    declared <- output$placeholders %||% list()
    values <- list(unit = units, arm = arms)
    for (nm in names(declared)) values[[nm]] <- as.character(unlist(declared[[nm]]))
    if (!all(used %in% c("unit", "arm", names(declared)))) {
        valAbort("output {path}: a placeholder is neither built in nor declared",
                 "val_spec_invalid")
    }
    if (!all(names(declared) %in% used)) {
        valAbort("output {path}: a declared placeholder is not used", "val_spec_invalid")
    }
    vals <- values[used]
    if (any(!grepl("^[A-Za-z0-9._-]+$", unlist(vals)))) {
        valAbort("output {path}: a placeholder value is not of the allowed form",
                 "val_spec_invalid")
    }
    if (length(used) == 0L) return(path)
    grid <- expand.grid(vals, stringsAsFactors = FALSE)
    vapply(seq_len(nrow(grid)), function(i) {
        p <- path
        for (nm in used) p <- gsub(paste0("{", nm, "}"), grid[[nm]][[i]], p, fixed = TRUE)
        p
    }, character(1))
}

#' Check one result table: typed columns, declared labels, unique label rows, a row bound, and no
#' universe key among the labels
#' @noRd
valCheckResultTable <- function(path, name, columns, labels, max_rows, universe = NULL) {
    forbid <- function(what) valAbort("{name}: {what}", "val_publish_forbidden")
    if (!all(columns %in% RESULT_TYPES)) forbid("a column type outside the result types")
    tab <- utils::read.delim(path, colClasses = "character", na.strings = character(),
                             check.names = FALSE, stringsAsFactors = FALSE)
    if (!identical(names(tab), names(columns))) forbid("columns differ from the declaration")
    if (nrow(tab) > min(max_rows, 1000L)) forbid("more rows than declared")
    label_cols <- names(columns)[columns == "label"]
    for (col in label_cols) {
        if (!all(tab[[col]] %in% labels[[col]])) forbid("a label outside the declared vocabulary")
        if (!is.null(universe) && any(tab[[col]] %in% universe)) forbid("a label is a feature key")
    }
    if (length(label_cols) > 0L && anyDuplicated(tab[, label_cols, drop = FALSE])) {
        forbid("label rows are not unique")
    }
    for (col in names(columns)[columns == "hash"]) {
        if (!all(grepl("^[0-9a-f]{64}$", tab[[col]]))) forbid("a hash column holds other text")
    }
    for (col in names(columns)[!columns %in% c("label", "hash")]) {
        v <- tab[[col]][nzchar(tab[[col]]) & tab[[col]] != "NA"]
        if (anyNA(suppressWarnings(as.numeric(v)))) forbid("a typed column holds text")
    }
    invisible(TRUE)
}

#' Header line of a tab-separated file
#' @noRd
valHeader <- function(path) {
    con <- if (grepl("\\.gz$", path)) gzfile(path, "r") else file(path, "r")
    on.exit(close(con))
    strsplit(readLines(con, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1L]]
}

#' The allow-list check of an evidence folder (section 5.2)
#'
#' Only the file names of section 5.7 and the declared outputs, each with exactly its columns;
#' `class` values only from the declared classes; row counts within their bounds; no file whose
#' SHA-256 equals a `not_permitted` input; result tables by their typed columns; keyed tables and
#' copies only where every data source in their lineage permits redistribution; `attribution.tsv`
#' listing every input; no file over 5 MB; no archive except `scored.tsv.gz` (never in phase
#' `lifecycle`); in phase `lifecycle` no call count of a primary unit.
#'
#' @param dir Evidence folder (`val-out/` or a committed run folder).
#' @param manifest The input manifest.
#' @param outputs The declaration's `outputs`.
#' @param arms,units Arm and unit ids for placeholder expansion.
#' @param universe Keys and feature identifiers of `universe.tsv` (for `not_permitted` lineages).
#' @param classes Declared truth classes.
#' @param row_bounds Named list: file name = maximum rows.
#' @param phase The phase.
#' @param scored_lineage Data roles in the lineage of `scored.tsv.gz` (default: every data role).
#' @param primary_units Units whose call counts stay private in the lifecycle.
#' @param forbidden_sha256 Further hashes never to leave (prepared inputs of `not_permitted` data).
#' @return `TRUE`, or a `val_publish_forbidden` refusal.
#' @export
valCheckEvidence <- function(dir, manifest, outputs = list(), arms = character(),
                             units = character(), universe = NULL, classes = NULL,
                             row_bounds = list(), phase = "counted", scored_lineage = NULL,
                             primary_units = character(), forbidden_sha256 = character()) {
    forbid <- function(what) valAbort("evidence: {what}", "val_publish_forbidden")
    closed <- manifest$kind == "data" & manifest$redistribution == "not_permitted"
    forbidden <- c(manifest$sha256[closed], forbidden_sha256)
    copies <- manifest[manifest$commit %in% TRUE, , drop = FALSE]
    declared <- list()
    for (o in outputs) {
        for (p in valExpandOutputPaths(o, arms = arms, units = units)) declared[[p]] <- o
    }
    att <- file.path(dir, "attribution.tsv")
    if (!file.exists(att)) forbid("attribution.tsv is missing")
    att_tab <- utils::read.delim(att, colClasses = "character", stringsAsFactors = FALSE)
    if (!setequal(att_tab$role, manifest$role)) forbid("attribution.tsv does not list every input")
    files <- list.files(dir, recursive = TRUE, all.files = TRUE)
    for (f in files) {
        path <- file.path(dir, f)
        base <- basename(f)
        if (file.size(path) > MAX_FILE_BYTES) forbid(paste(f, "is over 5 MB"))
        if (valSha256(path) %in% forbidden) forbid(paste(f, "copies a not-permitted input"))
        if (identical(f, "scored.tsv.gz")) {
            lineage <- scored_lineage %||% manifest$role[manifest$kind == "data"]
            if (phase == "lifecycle") forbid("scored.tsv.gz never leaves a lifecycle job")
            if (!valPublishable(manifest, lineage)) {
                forbid("scored.tsv.gz has a not-permitted source")
            }
            if (!identical(valHeader(path), SCORED_COLUMNS)) forbid("scored.tsv.gz columns")
            next
        }
        if (grepl(ARCHIVE_PATTERN, f, ignore.case = TRUE)) {
            forbid(paste(f, "is an archive or binary"))
        }
        if (!is.null(declared[[f]])) {
            valCheckDeclaredOutput(path, f, declared[[f]], manifest, universe, phase)
            next
        }
        if (startsWith(f, "inputs/") && base %in% copies$file) {
            row <- copies[copies$file == base, ][1L, ]
            if (!identical(valSha256(path), row$sha256)) {
                forbid(paste(f, "is not the declared copy"))
            }
            next
        }
        schema <- if (grepl("^ledger-.*\\.tsv$", base)) "ledger" else f
        if (!is.null(EVIDENCE_COLUMNS[[schema]])) {
            valCheckFixedTable(path, f, schema, classes, row_bounds, phase, primary_units)
            next
        }
        if (f %in% EVIDENCE_FREE_NAMES || grepl("^members-.*\\.tsv$", f)) next
        forbid(paste(f, "is not on the allow-list"))
    }
    TRUE
}

#' Check a fixed evidence table: columns, classes, row bound, lifecycle call counts
#' @noRd
valCheckFixedTable <- function(path, f, schema, classes, row_bounds, phase, primary_units) {
    forbid <- function(what) valAbort("evidence: {f} {what}", "val_publish_forbidden")
    tab <- utils::read.delim(path, colClasses = "character", na.strings = character(),
                             check.names = FALSE, stringsAsFactors = FALSE)
    if (!identical(names(tab), EVIDENCE_COLUMNS[[schema]])) forbid("has other columns")
    if (!is.null(classes) && "class" %in% names(tab) && schema != "cells.tsv") {
        if (!all(tab$class %in% c(classes, "", "NA", "all"))) forbid("has an undeclared class")
    }
    bound <- row_bounds[[f]]
    if (!is.null(bound) && nrow(tab) > bound) forbid("has more rows than its bound")
    if (phase == "lifecycle" && schema == "metrics.tsv" &&
        any(tab$metric %in% CALL_COUNT_METRICS & tab$unit %in% primary_units)) {
        forbid("holds a call count of a primary unit in the lifecycle")
    }
    invisible(TRUE)
}

#' Check a declared output: a keyed table only under a permitted lineage (never in the lifecycle);
#' a result table by its typed columns
#' @noRd
valCheckDeclaredOutput <- function(path, f, o, manifest, universe, phase) {
    lineage <- unlist(o$lineage)
    publishable <- valPublishable(manifest, lineage)
    cols <- vapply(o$columns, function(x) x$type, character(1))
    names(cols) <- vapply(o$columns, function(x) x$name, character(1))
    if (!identical(o$key %||% "none", "none")) {
        if (!publishable || phase == "lifecycle") {
            valAbort("evidence: {f} is a keyed table of a not-permitted or blind source",
                     "val_publish_forbidden")
        }
        if (!identical(valHeader(path), names(cols))) {
            valAbort("evidence: {f} has other columns", "val_publish_forbidden")
        }
        return(invisible(TRUE))
    }
    check_universe <- if (publishable) NULL else universe
    if (!is.null(check_universe)) {
        vals <- unlist(o$placeholders %||% list())
        if (any(vals %in% check_universe)) {
            valAbort("evidence: {f} names a feature key in its path", "val_publish_forbidden")
        }
    }
    labels <- lapply(o$labels %||% list(), function(x) as.character(unlist(x)))
    valCheckResultTable(path, f, cols, labels, o$max_rows %||% 1000L, check_universe)
}
