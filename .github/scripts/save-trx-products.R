# KIT-134 B-01: save the values behind the transcriptomics oracles that are stored only as digests
# (tests/testthat/helper-workflow-trx-declarations.R against docs/qa/workflow-execution/kit/kit-040), so that
# packaging/verify/compare-results.R can apply KIT-134's rule (revision 2.0.4) against the laptop's KIT-040 products.
# It repeats what the tests do (MCI-030.1, MCI-031.1, MCI-032.1 on trx_edge_12) through the same helpers, without
# their expectations; no test and no oracle is changed. Run from the package root.
#
# Usage: Rscript --vanilla .github/scripts/save-trx-products.R <output directory> [--corpus-only]
#
# Writes trx-corpus.rds: the regenerated corpus (each file's bytes and sha256, the manifest's sha256, the manifest, the
# counts table, truth.json parsed). Unless --corpus-only (cheap; the laptop side): trx-products.rds (top_genes,
# best_res, best_control_genes, ruviii_newY, the normalised TSV parsed, each journey's errors) and
# trx-products-sha256.tsv (each product's sha256 as the tests compute it, beside the KIT-040 fingerprint).

arguments <- commandArgs(trailingOnly = TRUE)
if (!length(arguments)) stop("usage: Rscript --vanilla .github/scripts/save-trx-products.R <output directory> [--corpus-only]")
out_dir <- arguments[[1L]]
corpus_only <- "--corpus-only" %in% arguments
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
corpus <- "trx_edge_12"

# The helper's constants, read from it so that they cannot drift.
helper_file <- "tests/testthat/helper-workflow-trx-declarations.R"
constants <- new.env()
for (expression in parse(helper_file)) {
    if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
        as.character(expression[[2L]]) %in% c("TRX_DECL_STORE_SHA256", "TRX_DECL_PER_GROUP")) eval(expression, constants)
}
store_path <- file.path("docs/qa/workflow-execution/kit/kit-040", sprintf("reference-%s.rds", corpus))
stopifnot(identical(digest::digest(file = store_path, algo = "sha256"), constants$TRX_DECL_STORE_SHA256[[corpus]]))
store <- readRDS(store_path)

# The corpus, as trxDeclReference() writes it.
corpus_dir <- file.path(tempdir(), "trx-corpus")
generator <- new.env(parent = globalenv())
sys.source("tools/profiling/kit_omics_corpus.R", generator)
generator$writeKitTrxCorpus(corpus_dir, generator$kitOmicsRecipe("trx_screen", constants$TRX_DECL_PER_GROUP[[corpus]]))
files <- sort(list.files(corpus_dir, recursive = TRUE))
corpus_value <- list(
    files = data.frame(file = files, bytes = unname(file.size(file.path(corpus_dir, files))),
        sha256 = unname(vapply(file.path(corpus_dir, files), function(f) digest::digest(file = f, algo = "sha256"), character(1)))),
    manifest_sha256 = digest::digest(file = file.path(corpus_dir, "manifest.json"), algo = "sha256"),
    manifest = jsonlite::read_json(file.path(corpus_dir, "manifest.json"), simplifyVector = TRUE),
    counts = utils::read.delim(file.path(corpus_dir, "counts.tsv"), check.names = FALSE, stringsAsFactors = FALSE),
    truth = jsonlite::read_json(file.path(corpus_dir, "truth.json"), simplifyVector = TRUE))
saveRDS(corpus_value, file.path(out_dir, "trx-corpus.rds"), version = 3)
cat(sprintf("corpus manifest sha256 %s (KIT-040 %s)\n", corpus_value$manifest_sha256, store$corpus_manifest_sha256))
print(corpus_value$files, row.names = FALSE)
if (corpus_only) quit(save = "no", status = 0L)

# The products, through the tests' own helpers. load_all(export_all = TRUE) attaches every package object; the helpers
# go into an environment under the global one (a helper defines an S4 class, which a locked namespace refuses).
pkgload::load_all(".", export_all = TRUE, helpers = FALSE, attach_testthat = FALSE, quiet = TRUE)
env <- new.env(parent = globalenv())
testthat::source_test_helpers("tests/testthat", env = env)
# The summary test's two local helpers (module_ci_trx_ruviii_value, module_ci_trx_summary_steps), taken from it.
for (expression in parse("tests/testthat/test-module-ci-trx-summary.R")) {
    if (is.call(expression) && identical(expression[[1L]], as.name("<-"))) eval(expression, env)
}
collect <- function() {
    # MCI-031.1: the DA tab's first contrast's topTags table.
    fixture <- module_ci_trx_project(corpus_dir)
    out <- module_ci_trx_tabs_journey(fixture)
    owner <- fixture$owner
    manifest <- withWorkflowStateScope(owner, function() resolveWorkflowStateValue(owner$.context, owner$.context$snapshot$da_results))
    top_genes <- trxDeclValue(owner, manifest$results[[1L]]$products$da_contrast)$assays[[1L]]$results$top_table
    # MCI-030.1: the automatic grid's winner (its own journey without DA, as the test).
    fixture_norm <- module_ci_trx_project(corpus_dir)
    out_norm <- module_ci_trx_tabs_journey(fixture_norm, da = FALSE)
    owner_norm <- fixture_norm$owner
    grid <- withWorkflowStateScope(owner_norm, function()
        trxProductValue(owner_norm, trxStateRuvs(owner_norm, getWorkflowStateReference(owner_norm))$grid))
    # MCI-032.1: the RUVIII export's newY and the integration TSV written from it.
    fixture_summary <- module_ci_trx_project(corpus_dir)
    assign("project_dirs", fixture_summary$project_dirs, envir = globalenv())
    out_summary <- module_ci_trx_tabs_journey(fixture_summary, ruviii = "group", summary = module_ci_trx_summary_steps)
    value <- module_ci_trx_ruviii_value(fixture_summary$owner)
    tsv <- file.path(fixture_summary$paths$integration_dir, "transcriptomics_normalized_data.tsv")
    list(top_genes = top_genes, best_control_genes = grid$control_genes,
        best_res = grid$best[c("percentage", "n_controls", "best_k", "separation_score", "k_penalty", "composite_score",
            "ctl_idx", "control_genes")],
        ruviii_newY = value$newY,
        normalized_tsv = if (file.exists(tsv)) list(sha256 = digest::digest(file = tsv, algo = "sha256"),
            table = as.data.frame(vroom::vroom(tsv, show_col_types = FALSE))) else NULL,
        errors = list(da = out$errors, norm = out_norm$errors, summary = out_summary$errors))
}
environment(collect) <- env
# The helpers' own expectations are reported, not raised (this is not a test run).
products <- withCallingHandlers(collect(), expectation = function(e) {
    if (inherits(e, c("expectation_failure", "expectation_error"))) message("helper expectation: ", conditionMessage(e))
    invokeRestart("continue_test")
})
saveRDS(products, file.path(out_dir, "trx-products.rds"), version = 3)

sha <- function(x) digest::digest(serialize(x, NULL, version = 3L), algo = "sha256", serialize = FALSE)
prints <- store$fingerprints
names <- c("top_genes", "best_control_genes", "best_res", "ruviii_newY")
digests <- data.frame(product = c(names, "transcriptomics_normalized_data.tsv"),
    runner_sha256 = c(vapply(products[names], sha, character(1)), products$normalized_tsv$sha256 %||% NA_character_),
    kit040_sha256 = c(prints$sha256[match(names, prints$product)],
        store$stable_outputs$sha256[store$stable_outputs$path == "results/transcriptomics/transcriptomics_normalized_data.tsv"]))
digests$equal <- digests$runner_sha256 == digests$kit040_sha256
utils::write.table(digests, file.path(out_dir, "trx-products-sha256.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(digests, row.names = FALSE)
cat("journey errors:", if (length(unlist(products$errors))) paste(unlist(products$errors), collapse = "; ") else "none", "\n")
