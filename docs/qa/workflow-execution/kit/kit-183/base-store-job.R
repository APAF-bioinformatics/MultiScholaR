# KIT-183 (VAL-64 (2)): the base store built and used inside one runner job. The store is made by the base commit's
# code and read by the head's tests T5, T8 and T15 in the same job; it never leaves the job (only result tables are
# written to <out dir>). Then X: each of three mutations of the head is loaded and the same tests run, which must fail.
#
# Usage (from the head package root, full git history): Rscript base-store-job.R <base commit> <out dir>
# Writes <out dir>/guards.tsv (variant, test, passed, failed, skipped, error) and <out dir>/verdict.txt.

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 2L)
base <- arguments[[1L]]
out_dir <- arguments[[2L]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)
head <- normalizePath(".")
work <- file.path(tempdir(), "k183-guards")
dir.create(work)
rscript <- file.path(R.home("bin"), "Rscript")

# 1. The store, by the base commit's code (its own tree, loaded in its own process).
base_tree <- file.path(work, "base")
stopifnot(system2("git", c("worktree", "add", "--detach", base_tree, base)) == 0L)
store_out <- file.path(work, "base-store")
capture <- file.path(head, "docs/qa/workflow-execution/kit/kit-183/capture-base-store.R")
code <- sprintf("suppressMessages(pkgload::load_all('%s', quiet = TRUE, export_all = TRUE, helpers = FALSE)); source('%s')",
    base_tree, capture)
status <- system2(rscript, c("--vanilla", "-e", shQuote(code), shQuote(store_out)), stdout = file.path(out_dir,
    "capture.log"), stderr = file.path(out_dir, "capture.log"))
stopifnot(status == 0L, dir.exists(file.path(store_out, "store", ".workflow-entities")))

# 2. The guards on a tree: T5, T8, T15 of test-integration-enrichment-tests.R, with K183_BASE_STORE set.
runGuards <- function(tree, variant) {
    result <- file.path(work, paste0(variant, ".tsv"))
    child <- sprintf(paste0(
        "suppressMessages(pkgload::load_all('%s', quiet = TRUE, export_all = TRUE, helpers = TRUE)); ",
        "r <- testthat::test_file('%s', reporter = testthat::ListReporter$new(), stop_on_failure = FALSE); ",
        "d <- as.data.frame(r); d <- d[grepl('KIT-183 T(5|8|15):', d$test), ]; ",
        "utils::write.table(data.frame(variant = '%s', test = substr(d$test, 1, 12), passed = d$nb - d$failed - d$skipped, ",
        "failed = d$failed, skipped = d$skipped, error = d$error), '%s', sep = '\\t', quote = FALSE, row.names = FALSE)"),
        tree, file.path(tree, "tests/testthat/test-integration-enrichment-tests.R"), variant, result)
    env <- paste0("K183_BASE_STORE=", file.path(store_out, "store"))
    system2(rscript, c("--vanilla", "-e", shQuote(child)), env = env, stdout = file.path(out_dir, paste0(variant, ".log")),
        stderr = file.path(out_dir, paste0(variant, ".log")))
    if (file.exists(result)) utils::read.delim(result, stringsAsFactors = FALSE) else
        data.frame(variant = variant, test = NA, passed = NA, failed = NA, skipped = NA, error = TRUE)
}

# 3. The mutations (ticket 9 A13): each on its own copy of the head tree.
mutations <- list(
    m1_no_restore = list(file = "R/utils_integration_enrichment.R",
        from = 'parameters$test <- if (identical(database, "string")) NULL else settings$test %||% "gsea"',
        to = 'parameters$test <- NULL'),
    m2_test_declared = list(file = "R/utils_integration_enrichment.R",
        from = 'parameters = c("ranking", "database", "id_column", "service_response", "mapping", "organism", "string_species",',
        to = 'parameters = c("ranking", "test", "database", "id_column", "service_response", "mapping", "organism", "string_species",',
        also_from = 'min_set_size = p$min_set_size, max_set_size = p$max_set_size, seed = p$seed)',
        also_to = 'min_set_size = p$min_set_size, max_set_size = p$max_set_size, seed = p$seed, test = p$test)'),
    m3_gsea_default = list(file = "R/utils_integration_enrichment.R",
        from = 'integrationEnrichmentGsea <- function(ranked, term2gene, term2name, seed, settings = integrationGseaSettings()) {',
        to = 'integrationEnrichmentGsea <- function(ranked, term2gene, term2name, seed, settings = integrationGseaSettings(min_size = 5L)) {'))
mutate <- function(name, m) {
    tree <- file.path(work, name)
    stopifnot(system2("git", c("worktree", "add", "--detach", tree, "HEAD")) == 0L)
    path <- file.path(tree, m$file)
    text <- readLines(path, warn = FALSE)
    replace <- function(text, from, to) {
        hit <- which(grepl(from, text, fixed = TRUE))
        stopifnot(length(hit) == 1L)
        text[[hit]] <- sub(from, to, text[[hit]], fixed = TRUE)
        text
    }
    text <- replace(text, m$from, m$to)
    if (!is.null(m$also_from)) text <- replace(text, m$also_from, m$also_to)
    writeLines(text, path)
    tree
}

rows <- list(runGuards(head, "head"))
for (name in names(mutations)) rows[[length(rows) + 1L]] <- runGuards(mutate(name, mutations[[name]]), name)
guards <- do.call(rbind, rows)
utils::write.table(guards, file.path(out_dir, "guards.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# A test id may cover several test_that blocks (T8 restated: (a), then (b) and (c)): every block must pass.
ok <- function(variant, test) {
    d <- guards[guards$variant == variant & guards$test == test, , drop = FALSE]
    nrow(d) >= 1L && all(d$failed == 0L) && !any(as.logical(d$skipped)) && !any(as.logical(d$error)) && all(d$passed > 0L)
}
verdict <- c(
    sprintf("head: T5 %s, T8 %s, T15 %s", ok("head", "KIT-183 T5: "), ok("head", "KIT-183 T8: "), ok("head", "KIT-183 T15:")),
    sprintf("m1_no_restore: T5 fails %s", !ok("m1_no_restore", "KIT-183 T5: ")),
    sprintf("m2_test_declared: T8 fails %s (T8 restated: (a) and (b) must fail)", !ok("m2_test_declared", "KIT-183 T8: ")),
    sprintf("m3_gsea_default: T15 fails %s", !ok("m3_gsea_default", "KIT-183 T15:")))
writeLines(verdict, file.path(out_dir, "verdict.txt"))
cat(verdict, sep = "\n")
for (tree in c(base_tree, file.path(work, names(mutations)))) system2("git", c("worktree", "remove", "--force", tree))
