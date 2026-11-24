# KIT-203 X: the mutation run (ticket section 11). Evidence code under docs/qa: no package code. Runs on a
# GitHub-hosted runner (run-tests `scripts`) from the package root of the commit under test (head).
#
# Usage: Rscript docs/qa/workflow-execution/kit/kit-203/mutate-k203.R <out dir> [<case> ...]
#   cases default: none m1 m2 m3 m4 (X); X2 (implementation audit K203-I03): none m5 m6 m7
#
# The head tree (git archive HEAD) is copied once per case: "none" (the control, unmutated) and the four mutations
# (m1) to (m4). Each mutation is an exact text replacement in R/ that must match exactly once (otherwise the case is
# an instrument defect: `applied` FALSE, no test run). Per case, a separate process runs
# pkgload::load_all(<copy>, export_all = TRUE, helpers = TRUE) and testthat::test_file() on the copy's
# tests/testthat/test-enrichment-significance.R, writing per test its passed, failed and error counts.
# The ticket's must-fail guards: S3 under (m1) and (m3); S1 under (m3); S4 under (m4); S14 under (m2). X2: S7 under
# (m5) (base's correction notice restored), S8 under (m6) (the plot never draws the family value), S9 under (m7) (one
# report template never states the family rule). An edit may span lines; it is matched in the file's whole text.
# Writes <out dir>/x-tests.tsv (case, test, passed, failed, errors, skipped), x-cases.tsv (case, file, applied) and
# x-verdict.tsv (case, guard, must_fail, failed_or_errored, pass). Exit status 0 (a guard that does not fail is a
# result).

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 1L)
# The out dir is created before it is made absolute (normalizePath leaves a missing path relative, and
# test_file runs from tests/testthat; R2's OD5 dumps failed on that, VAL-54).
dir.create(args[[1L]], showWarnings = FALSE, recursive = TRUE)
out <- normalizePath(args[[1L]], mustWork = TRUE)
work <- file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "k203-x")
dir.create(work, recursive = TRUE, showWarnings = FALSE)
rscript <- file.path(R.home("bin"), "Rscript")

orch <- file.path("R", "func_general_enrichment_orchestration.R")
flow <- file.path("R", "utils_workflow_enrichment.R")
plot <- file.path("R", "func_general_enrichment_engine_helpers.R")
lfq <- file.path("inst", "reports", "proteomics", "LFQ_report.rmd")
cases <- list(
    none = list(),
    m1 = list(list(file = orch, from = 'enrichmentSignificanceShipped <- "bh_whole_family"',
        to = 'enrichmentSignificanceShipped <- "qvalue"')),
    m2 = list(list(file = flow, from = 'significance_rule = parameters$significance_rule %||% "qvalue"), run_paths)',
        to = 'significance_rule = parameters$significance_rule %||% enrichmentSignificanceShipped), run_paths)')),
    m3 = list(list(file = orch, from = "m <- length(result@geneSets)", to = "m <- nrow(result@result)")),
    m4 = list(
        list(file = orch, from = 'if (is.null(result) || !identical(significance_rule, "bh_whole_family")) return(result)',
            to = "if (is.null(result)) return(result)"),
        list(file = orch, from = "significance_rule = if (family) significance_rule)",
            to = "significance_rule = significance_rule)")),
    m5 = list(list(file = flow, from = paste(
        "    if (!is.null(significance_rule)) {",
        "        unused <- !isTRUE(parameters$correction_method %in% c(\"fdr\", \"BH\"))",
        "        notify(paste0(ENRICHMENT_SIGNIFICANCE_RULES$notice[ENRICHMENT_SIGNIFICANCE_RULES$rule == significance_rule],",
        "            if (unused) sprintf(\" The selected correction, %s, applies to g:Profiler only and was not used.\",",
        "                parameters$correction_method)), if (unused) \"warning\" else \"message\", 8)",
        "    }", sep = "\n"),
        to = paste(
        "    if (!is.null(significance_rule) && !isTRUE(parameters$correction_method %in% c(\"fdr\", \"BH\"))) {",
        "        notify(sprintf(paste(\"clusterProfiler adjusts p-values with %s (Benjamini-Hochberg). The selected\",",
        "            \"correction, %s, applies to g:Profiler only and was not used.\"), correction,",
        "            parameters$correction_method), \"warning\", 8)",
        "    }", sep = "\n"))),
    m6 = list(list(file = plot, from = 'family <- identical(inputs$significance_rule, "bh_whole_family")',
        to = "family <- FALSE")),
    m7 = list(list(file = lfq,
        from = 'family_rule <- identical(numbers$enrichment_significance_rule %||% NA, "bh_whole_family")',
        to = "family_rule <- FALSE")))
selected <- if (length(args) > 1L) args[-1L] else c("none", "m1", "m2", "m3", "m4")
stopifnot(all(selected %in% names(cases)))
cases <- cases[selected]
guards <- data.frame(case = c("m1", "m3", "m3", "m4", "m2", "m5", "m6", "m7"),
    guard = c("KIT-203 S3", "KIT-203 S3", "KIT-203 S1", "KIT-203 S4", "KIT-203 S14", "KIT-203 S7", "KIT-203 S8",
        "KIT-203 S9"), stringsAsFactors = FALSE)
guards <- guards[guards$case %in% selected, , drop = FALSE]

case_rows <- list()
test_rows <- list()
for (case in names(cases)) {
    tree <- file.path(work, case)
    dir.create(tree, showWarnings = FALSE)
    system(sprintf("git archive HEAD | tar -x -C %s", shQuote(tree)))
    applied <- TRUE
    for (edit in cases[[case]]) {
        path <- file.path(tree, edit$file)
        text <- paste(readLines(path, warn = FALSE), collapse = "\n")
        ok <- sum(gregexpr(edit$from, text, fixed = TRUE)[[1L]] > 0L) == 1L
        if (ok) writeLines(sub(edit$from, edit$to, text, fixed = TRUE), path)
        case_rows[[length(case_rows) + 1L]] <- data.frame(case = case, file = edit$file, applied = ok)
        applied <- applied && ok
    }
    if (!length(cases[[case]])) case_rows[[length(case_rows) + 1L]] <- data.frame(case = case, file = "", applied = TRUE)
    if (!applied) next
    counts <- file.path(work, paste0(case, ".tsv"))
    code <- paste0("suppressMessages(pkgload::load_all(", deparse(tree), ", export_all = TRUE, helpers = TRUE, ",
        "quiet = TRUE)); res <- testthat::test_file(file.path(", deparse(tree), ", 'tests', 'testthat', ",
        "'test-enrichment-significance.R'), reporter = 'silent', stop_on_failure = FALSE); ",
        "df <- as.data.frame(res); utils::write.table(data.frame(test = df$test, passed = df$passed, ",
        "failed = df$failed, errors = as.integer(df$error), skipped = as.integer(df$skipped)), ", deparse(counts),
        ", sep = '\\t', quote = FALSE, row.names = FALSE)")
    log <- file.path(out, paste0("x-", case, ".log.txt"))
    system2(rscript, c("-e", shQuote(code)), stdout = log, stderr = log, env = "NOT_CRAN=true")
    if (file.exists(counts)) {
        x <- utils::read.delim(counts, stringsAsFactors = FALSE, quote = "")
        test_rows[[case]] <- cbind(case = case, x)
    }
}
tests <- do.call(rbind, test_rows)
write <- function(x, name) utils::write.table(x, file.path(out, name), sep = "\t", quote = FALSE, row.names = FALSE)
write(do.call(rbind, case_rows), "x-cases.tsv")
write(tests, "x-tests.tsv")
verdict <- do.call(rbind, lapply(c(seq_len(nrow(guards))), function(i) {
    rows <- tests[tests$case == guards$case[[i]] & startsWith(tests$test, paste0(guards$guard[[i]], ":")), ]
    bad <- sum(rows$failed) + sum(rows$errors)
    data.frame(case = guards$case[[i]], guard = guards$guard[[i]], must_fail = TRUE,
        failed_or_errored = bad, pass = nrow(rows) == 1L && bad > 0L)
}))
control <- tests[tests$case == "none", ]
verdict <- rbind(verdict, data.frame(case = "none", guard = "every test", must_fail = FALSE,
    failed_or_errored = sum(control$failed) + sum(control$errors),
    pass = nrow(control) > 0L && sum(control$failed) + sum(control$errors) == 0L))
write(verdict, "x-verdict.tsv")
print(verdict)
