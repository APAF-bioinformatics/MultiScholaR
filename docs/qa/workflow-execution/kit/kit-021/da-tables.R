# KIT-021 acceptance criterion 2 (owner side): the DA tables the oracle run saved on the kit_021 owner path, read back
# from the oracle project's store through the DA module's own reader and hashed with the canonical digest of
# control-check.R. Run with the bound build (kit_021) on R_LIBS and the installed package:
#   R_LIBS=<proof>/core/browser-smoke.<kit_021>/library Rscript docs/qa/workflow-execution/kit/kit-021/da-tables.R lfq_300
suppressMessages(library(MultiScholaR))
args <- commandArgs(TRUE)
corpus_key <- if (length(args)) args[[1L]] else "lfq_300"
repo <- normalizePath(".")
ns <- asNamespace("MultiScholaR")
api <- function(name) get(name, ns)
`%||%` <- function(a, b) if (is.null(a)) b else a
table_hash <- function(table) {
    table <- as.data.frame(table)
    table <- table[do.call(order, unname(as.list(table[, seq_len(min(2L, ncol(table))), drop = FALSE]))), order(names(table)), drop = FALSE]
    rownames(table) <- NULL
    digest::digest(lapply(table, function(column) if (is.numeric(column)) signif(column, 12) else as.character(column)), algo = "sha256")
}
# KIT-054: the evidence directory (oracles in, tables out) may be another ticket's (KIT_CHECK_DIR).
recorded <- readRDS(file.path(repo, Sys.getenv("KIT_CHECK_DIR", "docs/qa/workflow-execution/kit/kit-021"), paste0("oracles-", corpus_key, ".rds")))
root <- recorded$store
if (identical(basename(root), ".workflow-entities")) root <- dirname(root)
stopifnot(dir.exists(root), file.exists(file.path(root, ".workflow-entities", "store.json")))
identity <- api("readSavedStoreIdentity")(file.path(root, ".workflow-entities", "store.json"))
manager <- api("newDIAStateAdapter")(api("newWorkflowEntityStore")(root, identity$project_id, identity$workflow_id), reopen = TRUE)
on.exit(api("closeDIAStateAdapter")(manager), add = TRUE)
out <- api("withDIAStateScope")(manager, function() {
    manifest <- api("resolveDIAStateValue")(manager$.context, manager$.context$snapshot$da_results)
    stopifnot(identical(manifest$schema, "dia.da_results.v1"))
    tables <- lapply(manifest$results, function(r) lapply(
        api("readDIAContrastOutput")(r, fields = c("da_proteins_long", "da_proteins_wide", "contrasts_results_table")), table_hash))
    rows <- lapply(manifest$results, function(r) {
        long <- api("readDIAContrastOutput")(r, fields = "da_proteins_long")$da_proteins_long
        nrow(long)
    })
    list(corpus = corpus_key, corpus_sha256 = recorded$corpus_sha256, store = root, head_sha256 = recorded$head_sha256,
        build = recorded$build, contrasts = manifest$contrasts, model = manifest$model, table_hashes = tables, rows = rows)
})
out$library <- list(path = find.package("MultiScholaR"), source_commit = utils::packageDescription("MultiScholaR")$`Config/MultiScholaR/SourceCommit` %||% NA_character_)
jsonlite::write_json(out, file.path(repo, Sys.getenv("KIT_CHECK_DIR", "docs/qa/workflow-execution/kit/kit-021"), paste0("da-tables-", corpus_key, ".json")),
    auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = I(17))
cat("da tables", corpus_key, "| contrasts", paste(names(out$table_hashes), collapse = ","), "| rows", paste(unlist(out$rows), collapse = ","), "\n")
