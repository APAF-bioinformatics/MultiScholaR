# KIT-268 G0 (observation only; no expectation; no stored file touched).
# Usage: Rscript g0-observe.R <mode>   (mode = "control" or a tests/testthat file name); run from the package root.
# One process per mode (the runner's per-file convention): pkgload::load_all(export_all, helpers), then
# trace(completeOmicDaRun) in the namespace and in the attached copy; each call prints one row to <mode>-calls.tsv and one
# row per (assay, key) of results$qvalue_warnings to <mode>-keys.tsv, with the test name, the lane and step (the GUI-entry
# harness's frame variables) and the contrasts table's `contrasts` and `friendly_names`.
mode <- commandArgs(trailingOnly = TRUE)[[1L]]
pkg <- normalizePath(".")
out <- file.path(pkg, "kit268-g0") # absolute: testthat::test_file changes the working directory
dir.create(out, showWarnings = FALSE)
tsv <- function(name) file.path(out, paste0(mode, "-", name, ".tsv"))
if (!identical(mode, "control")) {
  control <- file.path(out, "control-keys.tsv")
  n <- if (file.exists(control)) length(readLines(control)) - 1L else 0L
  if (n < 1L) {
    cat("INSTRUMENT BLIND: the positive control printed no fallback comparison; not reading", mode, "\n")
    quit(status = 3L)
  }
}
suppressMessages(pkgload::load_all(pkg, quiet = TRUE, export_all = TRUE, helpers = TRUE))

.g0 <- new.env()
.g0$seq <- 0L
.g0$file <- mode
.g0$test <- NA_character_
assign(".g0", .g0, envir = globalenv())
flat <- function(x) paste(x, collapse = "|")
cat("seq\tfile\ttest\tlane\tstep\tn_assays\tn_keys\tfailed_contrasts\tcontrasts\tfriendly_names\n", file = tsv("calls"))
cat("seq\tfile\ttest\tlane\tstep\tassay\tkey\texact_match\texact_friendly\ttext_before_eq\trule_result\trule_result_in_friendly\n",
  file = tsv("keys"))
assign("g0_record", function(results, workflowData, frames) {
  tryCatch({
    g <- get(".g0", envir = globalenv())
    g$seq <- g$seq + 1L
    lane <- NA_character_
    step <- NA_character_
    for (fr in rev(frames)) {
      if (is.na(step) && exists("step", envir = fr, inherits = FALSE)) {
        s <- get("step", envir = fr, inherits = FALSE)
        if (is.list(s) && is.character(s$id) && length(s$id) == 1L) step <- s$id
      }
      if (is.na(lane) && exists("lane_id", envir = fr, inherits = FALSE)) {
        l <- get("lane_id", envir = fr, inherits = FALSE)
        if (is.character(l) && length(l) == 1L) lane <- l
      }
    }
    tbl <- tryCatch(shiny::isolate(workflowData$contrasts_tbl), error = function(e) NULL)
    ctr <- if (is.data.frame(tbl) && "contrasts" %in% names(tbl)) as.character(tbl$contrasts) else character()
    fri <- if (is.data.frame(tbl) && "friendly_names" %in% names(tbl)) as.character(tbl$friendly_names) else character()
    qw <- results$qvalue_warnings
    rows <- list()
    for (assay in names(qw)) {
      item <- qw[[assay]]
      keys <- if (is.list(item)) names(item) else assay
      for (key in keys) rows[[length(rows) + 1L]] <- c(if (is.list(item)) assay else NA_character_, key)
    }
    cat(g$seq, g$file, g$test, lane, step, length(qw), length(rows), flat(names(results$failed_contrasts)),
      flat(ctr), flat(fri), sep = "\t", file = tsv("calls"), append = TRUE)
    cat("\n", file = tsv("calls"), append = TRUE)
    for (row in rows) {
      key <- row[[2L]]
      exact <- key %in% ctr
      exact_friendly <- if (exact) fri[match(key, ctr)] else NA_character_
      before <- sub("=.*$", "", key)
      rule <- if (exact) exact_friendly else if (grepl("=", key, fixed = TRUE)) before else key
      cat(g$seq, g$file, g$test, lane, step, row[[1L]], key, exact, exact_friendly, before, rule, rule %in% fri,
        sep = "\t", file = tsv("keys"), append = TRUE)
      cat("\n", file = tsv("keys"), append = TRUE)
    }
  }, error = function(e) cat("g0_record error:", conditionMessage(e), "\n", file = tsv("errors"), append = TRUE))
  invisible(NULL)
}, envir = globalenv())
# tsv() is looked up from the tracer's frame: make it and `mode` visible there.
assign("tsv", tsv, envir = globalenv())
assign("mode", mode, envir = globalenv())
tracer <- quote(g0_record(results, workflowData, sys.frames()))
for (where in list(asNamespace("MultiScholaR"), as.environment("package:MultiScholaR"))) {
  suppressMessages(trace("completeOmicDaRun", tracer = tracer, print = FALSE, where = where))
  cat("traced in", environmentName(where), ":", is(get("completeOmicDaRun", envir = where), "functionWithTrace"), "\n")
}

G0Reporter <- R6::R6Class("G0Reporter", inherit = testthat::Reporter, public = list(
  start_test = function(context, test) { assign("test", gsub("[\t\r\n]+", " ", test), envir = get(".g0", envir = globalenv())); invisible() }))
path <- if (identical(mode, "control")) file.path("docs/qa/workflow-execution/kit/kit-268/g0/g0-control.R") else
  file.path("tests/testthat", mode)
list_reporter <- testthat::ListReporter$new()
reporter <- testthat::MultiReporter$new(list(testthat::SummaryReporter$new(), list_reporter, G0Reporter$new()))
testthat::test_file(normalizePath(path), reporter = reporter, stop_on_failure = FALSE)
df <- as.data.frame(list_reporter$get_results())
utils::write.table(df[, intersect(c("file", "test", "nb", "passed", "failed", "skipped", "error", "real"), names(df))],
  tsv("tests"), sep = "\t", quote = FALSE, row.names = FALSE)
calls <- utils::read.delim(tsv("calls"), stringsAsFactors = FALSE, quote = "")
keys <- utils::read.delim(tsv("keys"), stringsAsFactors = FALSE, quote = "")
cat("\n=== G0", mode, ": completeOmicDaRun calls:", nrow(calls), "; fallback keys:", nrow(keys), "\n")
if (nrow(keys)) print(keys[, c("seq", "test", "lane", "step", "assay", "key", "rule_result", "rule_result_in_friendly")], row.names = FALSE)
if (identical(mode, "control")) {
  cat("POSITIVE CONTROL fallback comparisons:", if (nrow(keys)) paste(keys$assay, keys$key, sep = ":", collapse = "; ") else "NONE", "\n")
  if (!nrow(keys)) quit(status = 2L)
}
