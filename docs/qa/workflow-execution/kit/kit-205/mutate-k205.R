# KIT-205 X: the mutation run (ticket 10.1). Evidence code under docs/qa: no package code. Runs on a GitHub-hosted runner
# (run-tests `scripts`) from the package root of the commit under test (head). The form of KIT-203's mutate-k203.R.
#
# Usage: Rscript docs/qa/workflow-execution/kit/kit-205/mutate-k205.R <out dir>
#
# The head tree (git archive HEAD) is copied once per case: "none" (the control, unmutated) and the three mutations.
# Each mutation is an exact text replacement in R/ that must match exactly once (otherwise the case is an instrument
# defect: `applied` FALSE, no test run). Per case, a separate process runs pkgload::load_all(<copy>, export_all = TRUE,
# helpers = TRUE) and testthat::test_file() on the copy's tests/testthat/test-kit205-enrichment-background.R, writing per
# test its passed, failed and error counts.
#   (m1) the reduction pooled over every comparison of the run (a protein is kept only where it is tested in all);
#   (m2) the record replay's wrapper around the enrichment step removed;
#   (m3) the reduction on a finite log2FC instead of the q- and p-values.
# Must fail: B2 under (m1), B6 (b, c) under (m2), B5 under (m3); the control passes every test.
# Writes <out dir>/x-tests.tsv, x-cases.tsv and x-verdict.tsv. Exit status 0 (a guard that does not fail is a result).

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
dir.create(args[[1L]], showWarnings = FALSE, recursive = TRUE)
out <- normalizePath(args[[1L]], mustWork = TRUE)
work <- file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "k205-x")
dir.create(work, recursive = TRUE, showWarnings = FALSE)
rscript <- file.path(R.home("bin"), "Rscript")

orch <- file.path("R", "func_general_enrichment_orchestration.R")
record <- file.path("R", "utils_workflow_analysis_record.R")
cases <- list(
    none = list(),
    m1 = list(list(file = orch, from = paste(
        "        reduced <- enrichmentTestedRows(de_data)",
        "        de_data <- reduced$rows",
        "        ",
        "        # Check if the column has valid data", sep = "\n"),
        to = paste(
        "        pooled <- Reduce(intersect, lapply(da_results@da_data, function(d)",
        "          as.character(enrichmentTestedRows(d)$rows[[protein_id_column]])))",
        "        reduced <- list(rows = de_data[de_data[[protein_id_column]] %in% pooled, , drop = FALSE], left_out = 0L)",
        "        de_data <- reduced$rows",
        "",
        "        # Check if the column has valid data", sep = "\n"))),
    m2 = list(list(file = record,
        from = 'withr::with_options(workflowReplayRuleOptions(step, workflowOperationEntry("enrichment")$declaration),',
        to = "identity(")),
    m3 = list(list(file = orch, from = paste(
        "  tested <- is.finite(de_data[[\"fdr_qvalue\"]])",
        "  if (\"raw_pvalue\" %in% names(de_data)) tested <- tested & is.finite(de_data[[\"raw_pvalue\"]])", sep = "\n"),
        to = "  tested <- is.finite(de_data[[\"log2FC\"]])")))
guards <- data.frame(case = c("m1", "m2", "m3"), guard = c("KIT-205 B2", "KIT-205 B6 (b, c)", "KIT-205 B5"),
    stringsAsFactors = FALSE)

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
        "'test-kit205-enrichment-background.R'), reporter = 'silent', stop_on_failure = FALSE); ",
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
verdict <- do.call(rbind, lapply(seq_len(nrow(guards)), function(i) {
    rows <- tests[tests$case == guards$case[[i]] & startsWith(tests$test, paste0(guards$guard[[i]], ":")), ]
    bad <- sum(rows$failed) + sum(rows$errors)
    data.frame(case = guards$case[[i]], guard = guards$guard[[i]], must_fail = TRUE, failed_or_errored = bad,
        pass = nrow(rows) == 1L && bad > 0L)
}))
control <- tests[tests$case == "none", ]
verdict <- rbind(verdict, data.frame(case = "none", guard = "every test", must_fail = FALSE,
    failed_or_errored = sum(control$failed) + sum(control$errors),
    pass = nrow(control) > 0L && sum(control$failed) + sum(control$errors) == 0L))
write(verdict, "x-verdict.tsv")
print(verdict)
