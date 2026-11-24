# KIT-184 (VAL-64 (2)): the base stores built and used inside one runner job. The stores and the inputs are made by the
# base commit's code (capture-base-store.R, loaded with the base tree in its own process) and read by the head's tests
# in the same job (tests/testthat/test-integration-removed-factors.R with K184_BASE_STORE set: K3's base parts, K7 (i)
# (b), (c) and (iii), K9 (a), K10's second case); nothing leaves the job but the result tables. With `mutations` (run X,
# ticket 9 A16), each of four mutations of the head is then loaded on its own copy of the tree and the same file run,
# and the guards named below must fail.
#
# Usage (from the head package root, full git history): Rscript base-store-job.R <base commit> <out dir> [mutations]
# Writes <out dir>/guards.tsv (variant, test, passed, failed, skipped, error), <out dir>/verdict.txt, the capture's
# provenance.tsv and one log per variant.

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) %in% 2:3)
base <- arguments[[1L]]
out_dir <- arguments[[2L]]
with_mutations <- length(arguments) == 3L && identical(arguments[[3L]], "mutations")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)
head <- normalizePath(".")
work <- file.path(tempdir(), "k184-guards")
dir.create(work)
rscript <- file.path(R.home("bin"), "Rscript")

# 1. The stores, by the base commit's code (its own tree, loaded in its own process, its tree the working directory).
base_tree <- file.path(work, "base")
stopifnot(system2("git", c("worktree", "add", "--detach", base_tree, base)) == 0L)
store_out <- file.path(work, "base-store")
capture <- file.path(head, "docs/qa/workflow-execution/kit/kit-184/capture-base-store.R")
code <- sprintf(paste0("setwd('%s'); suppressMessages(pkgload::load_all('.', quiet = TRUE, export_all = TRUE, ",
    "helpers = FALSE)); source('%s')"), base_tree, capture)
status <- system2(rscript, c("--vanilla", "-e", shQuote(code), shQuote(store_out)), stdout = file.path(out_dir,
    "capture.log"), stderr = file.path(out_dir, "capture.log"))
stopifnot(status == 0L, dir.exists(file.path(store_out, "corrected", "store", ".workflow-entities")),
    dir.exists(file.path(store_out, "default", "store", ".workflow-entities")))
file.copy(file.path(store_out, "provenance.tsv"), file.path(out_dir, "provenance.tsv"))

# 2. The new test file on a tree, with K184_BASE_STORE set: one row per test_that block.
runFile <- function(tree, variant) {
    result <- file.path(work, paste0(variant, ".tsv"))
    child <- sprintf(paste0(
        "suppressMessages(pkgload::load_all('%s', quiet = TRUE, export_all = TRUE, helpers = TRUE)); ",
        "r <- testthat::test_file('%s', reporter = testthat::ListReporter$new(), stop_on_failure = FALSE); ",
        "d <- as.data.frame(r); ",
        "utils::write.table(data.frame(variant = '%s', test = d$test, passed = d$nb - d$failed - d$skipped, ",
        "failed = d$failed, skipped = d$skipped, error = d$error), '%s', sep = '\\t', quote = FALSE, row.names = FALSE)"),
        tree, file.path(tree, "tests/testthat/test-integration-removed-factors.R"), variant, result)
    env <- c(paste0("K184_BASE_STORE=", store_out), "NOT_CRAN=true")
    system2(rscript, c("--vanilla", "-e", shQuote(child)), env = env, stdout = file.path(out_dir, paste0(variant, ".log")),
        stderr = file.path(out_dir, paste0(variant, ".log")))
    if (file.exists(result)) utils::read.delim(result, stringsAsFactors = FALSE, quote = "") else
        data.frame(variant = variant, test = NA, passed = NA, failed = NA, skipped = NA, error = TRUE)
}

# 3. The mutations (ticket 9 A16, section 11 run X): each on its own copy of the head tree, one line replaced.
mutations <- list(
    # K3's guard: an (empty) `removed_factors` element added to every omic's entry of the inputs state.
    m1_always_element = list(file = "R/utils_workflow_integration_model.R",
        from = "            if (is.null(removed)) next",
        to = "            if (is.null(removed)) { record[[omic]][\"removed_factors\"] <- list(record[[omic]]$removed_factors); next }"),
    # K7 (ii): covariates whenever the input has the factors, whatever the model was trained on.
    m2_covariates_default = list(file = "R/utils_integration_enrichment.R",
        from = "    adjust <- if (is.null(given)) identical(model_parameters$input_matrix, \"ruv_corrected\") else given",
        to = paste0("    adjust <- if (is.null(given)) !is.null(object@inputs[[integrationViewOmics(object)[[p$view]]]]",
            "$removed_factors[[p$view]]) else given")),
    # K9: the replay's restore dropped (an entry without the setting takes the default for its model).
    m3_no_restore = list(file = "R/utils_integration_enrichment.R",
        from = "            settings$adjust_ruv_factors %||% FALSE",
        to = "            settings$adjust_ruv_factors"),
    # K7 (i) (a) and (iii): the setting put among the declared parameters.
    m4_setting_declared = list(file = "R/utils_integration_enrichment.R",
        from = "            \"test_settings\", \"values_slot\", \"design_covariates\"),",
        to = "            \"test_settings\", \"values_slot\", \"design_covariates\", \"adjust_ruv_factors\"),",
        also_from = "                design_covariates = p$design_covariates %||% \"none\")",
        also_to = "                design_covariates = p$design_covariates %||% \"none\", adjust_ruv_factors = p$adjust_ruv_factors)"))
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

rows <- list(runFile(head, "head"))
if (with_mutations) {
    for (name in names(mutations)) rows[[length(rows) + 1L]] <- runFile(mutate(name, mutations[[name]]), name)
}
guards <- do.call(rbind, rows)
utils::write.table(guards, file.path(out_dir, "guards.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# A test id may cover several test_that blocks: every block whose name starts with the id must run and pass.
ok <- function(variant, id) {
    d <- guards[guards$variant == variant & startsWith(as.character(guards$test), id), , drop = FALSE]
    nrow(d) >= 1L && all(d$failed == 0L) && !any(as.logical(d$skipped)) && !any(as.logical(d$error)) && all(d$passed > 0L)
}
ids <- c(k3_base = "KIT-184 K3 (base store)", k3_import = "KIT-184 K3: the import", k7_ia = "KIT-184 K7 (i) (a)",
    k7_ib = "KIT-184 K7 (i) (b)", k7_ii = "KIT-184 K7 (ii), K16", k7_iii_k9a = "KIT-184 K7 (iii), K9 (a)")
head_rows <- guards[guards$variant == "head", , drop = FALSE]
verdict <- c(
    sprintf("head: %d test blocks, %d with a failure or an error, %d skipped", nrow(head_rows),
        sum(head_rows$failed > 0L | as.logical(head_rows$error)), sum(as.logical(head_rows$skipped))),
    sprintf("head: %s passes %s", ids, vapply(ids, function(id) ok("head", id), logical(1))))
if (with_mutations) {
    verdict <- c(verdict,
        sprintf("m1_always_element: K3 (base store) fails %s; K3 (the import) fails %s", !ok("m1_always_element", ids[["k3_base"]]),
            !ok("m1_always_element", ids[["k3_import"]])),
        sprintf("m2_covariates_default: K7 (ii) fails %s", !ok("m2_covariates_default", ids[["k7_ii"]])),
        sprintf("m3_no_restore: K9 (a) fails %s; the replay of an entry without the setting (in K7 (ii), K16) fails %s",
            !ok("m3_no_restore", ids[["k7_iii_k9a"]]), !ok("m3_no_restore", ids[["k7_ii"]])),
        sprintf("m4_setting_declared: K7 (i) (a) fails %s; K7 (iii) fails %s", !ok("m4_setting_declared", ids[["k7_ia"]]),
            !ok("m4_setting_declared", ids[["k7_iii_k9a"]])))
}
writeLines(verdict, file.path(out_dir, "verdict.txt"))
cat(verdict, sep = "\n")
for (tree in c(base_tree, if (with_mutations) file.path(work, names(mutations)))) {
    system2("git", c("worktree", "remove", "--force", tree))
}
