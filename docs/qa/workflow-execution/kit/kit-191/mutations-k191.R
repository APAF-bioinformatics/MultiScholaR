# KIT-191 X (ticket 9 A12, 11): the guards shown failing on deliberate mutations of the head. Each mutation is applied
# to its own copy of the head tree (a git worktree), loaded in a fresh process, and the guards U5, U7 and U8 of
# tests/testthat/test-integration-enrichment-options.R run: (m1) camera's canonical text changed: U7 must fail; (m2)
# the replay drops the recorded test: U8 (in "U5, U8") must fail; (m3) the permutation's text without the seed: U5 (in
# "U5, U8") must fail. The unmutated head is the control (all pass). AUDIT-KIT-191-impl I07: with a base commit given,
# m2 also replays the branch-point store's old entries (base-head.R run from the m2 tree), which must then fail.
#
# Usage (from the head package root, full git history): Rscript mutations-k191.R <out dir> [<base commit>]
# Writes <out dir>/guards.tsv, verdict.txt and each variant's log.

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) %in% 1:2)
out_dir <- arguments[[1L]]
base <- if (length(arguments) == 2L) arguments[[2L]] else NA_character_
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)
head <- normalizePath(".")
work <- file.path(tempdir(), "k191-mutations")
dir.create(work)
rscript <- file.path(R.home("bin"), "Rscript")

runGuards <- function(tree, variant) {
    result <- file.path(work, paste0(variant, ".tsv"))
    child <- sprintf(paste0(
        "suppressMessages(pkgload::load_all('%s', quiet = TRUE, export_all = TRUE, helpers = TRUE)); ",
        "r <- testthat::test_file('%s', reporter = testthat::ListReporter$new(), stop_on_failure = FALSE); ",
        "d <- as.data.frame(r); d <- d[grepl('KIT-191 U(5|7)', d$test), ]; ",
        "utils::write.table(data.frame(variant = '%s', test = substr(d$test, 1, 16), passed = d$nb - d$failed - d$skipped, ",
        "failed = d$failed, skipped = d$skipped, error = d$error), '%s', sep = '\\t', quote = FALSE, row.names = FALSE)"),
        tree, file.path(tree, "tests/testthat/test-integration-enrichment-options.R"), variant, result)
    system2(rscript, c("--vanilla", "-e", shQuote(child)), stdout = file.path(out_dir, paste0(variant, ".log")),
        stderr = file.path(out_dir, paste0(variant, ".log")))
    if (file.exists(result)) utils::read.delim(result, stringsAsFactors = FALSE) else
        data.frame(variant = variant, test = NA, passed = NA, failed = NA, skipped = NA, error = TRUE)
}

mutations <- list(
    m1_camera_text = list(file = "R/utils_integration_enrichment.R",
        from = 'camera = paste0("inter.gene.cor=NA;allow.neg.cor=TRUE;use.ranks=FALSE;',
        to = 'camera = paste0("inter.gene.cor=NA;allow.neg.cor=TRUE;use.ranks=TRUE;'),
    m2_replay_drops_test = list(file = "R/utils_integration_enrichment.R",
        from = 'parameters$test <- if (identical(database, "string")) NULL else settings$test %||% "gsea"',
        to = 'parameters$test <- NULL'),
    m3_permutation_text_no_seed = list(file = "R/utils_integration_enrichment.R",
        from = '"seed=%d;",',
        to = '"",'))
mutate <- function(name, m) {
    tree <- file.path(work, name)
    stopifnot(system2("git", c("worktree", "add", "--detach", tree, "HEAD")) == 0L)
    path <- file.path(tree, m$file)
    text <- readLines(path, warn = FALSE)
    hit <- which(grepl(m$from, text, fixed = TRUE))
    stopifnot(length(hit) == 1L)
    text[[hit]] <- sub(m$from, m$to, text[[hit]], fixed = TRUE)
    writeLines(text, path)
    tree
}

rows <- list(runGuards(head, "head"))
trees <- list()
for (name in names(mutations)) {
    trees[[name]] <- mutate(name, mutations[[name]])
    rows[[length(rows) + 1L]] <- runGuards(trees[[name]], name)
}
# The old entries under m2: the branch-point store replayed by the mutated code (its own process tree, as base-head.R).
old_entries <- NA
if (!is.na(base)) {
    old_out <- file.path(out_dir, "m2-old-entries")
    tree <- trees[["m2_replay_drops_test"]]
    withr::with_dir(tree, system2(rscript, c("--vanilla", "docs/qa/workflow-execution/kit/kit-191/base-head.R", base,
        old_out), stdout = file.path(out_dir, "m2-old-entries.log"), stderr = file.path(out_dir, "m2-old-entries.log")))
    verdict <- file.path(old_out, "verdict.txt")
    replay_ok <- file.exists(verdict) && any(grepl("^replay_terms_identical\tTRUE$", readLines(verdict))) &&
        any(grepl("^replay_tests_kept\tTRUE$", readLines(verdict)))
    old_entries <- !replay_ok
}
guards <- do.call(rbind, rows)
utils::write.table(guards, file.path(out_dir, "guards.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
failing <- function(variant, pattern) {
    rows <- guards[guards$variant == variant & grepl(pattern, guards$test), , drop = FALSE]
    nrow(rows) > 0L && any(rows$failed > 0L | rows$error)
}
control <- guards[guards$variant == "head", , drop = FALSE]
checks <- c(head_passes = nrow(control) >= 2L && all(control$failed == 0L & !control$error),
    m1_fails_U7 = failing("m1_camera_text", "U7"),
    m2_fails_U8 = failing("m2_replay_drops_test", "U5, U8"),
    m3_fails_U5 = failing("m3_permutation_text_no_seed", "U5, U8"),
    m2_fails_old_entries = if (is.na(base)) NA else old_entries)
writeLines(sprintf("%s\t%s", names(checks), checks), file.path(out_dir, "verdict.txt"))
print(guards)
print(checks)
stopifnot(all(checks, na.rm = TRUE))
