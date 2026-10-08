# CORE-014 premise probe (2026-10-08, read-only; not a counted run).
# Classifies a store's entities against the closure Save would archive, by descriptor kind and identity.
# Usage: Rscript store-closure-probe.R <project>/.workflow-entities   (from the repository root)
.libPaths(c(Sys.getenv("CORE_LIBRARY"), .libPaths()))
ns <- asNamespace("MultiScholaR"); store_dir <- commandArgs(TRUE)[1]
project <- dirname(store_dir); id <- jsonlite::read_json(file.path(store_dir, "store.json"))
store <- ns$newWorkflowEntityStore(project, id$project_id, id$workflow_id, max_cache_bytes = 0)
manager <- ns$newDIAStateAdapter(store, reopen = TRUE)
ctx <- manager$.context
refs <- ns$expandWorkflowTableReferences(store, ns$collectDIANativeReferences(list(ctx$snapshot, ctx$head)))
live <- unique(vapply(refs, function(r) basename(dirname(r$payload$relative_path)), character(1)))
ents <- list.dirs(file.path(store_dir, "entities"), recursive = FALSE)
size <- vapply(ents, function(d) sum(file.size(list.files(d, recursive = TRUE, full.names = TRUE))), numeric(1))
dead <- !basename(ents) %in% live
kind <- sub("-[0-9a-f]+$", "", basename(ents))
cat(sprintf("RESULT entities=%d closure=%d dead=%d dead_bytes=%.3fGB of %.3fGB states=%s\n", length(ents), sum(!dead), sum(dead), sum(size[dead])/1e9, sum(size)/1e9, paste(names(ctx$snapshot$states), collapse=",")))
print(tapply(size[dead]/1e6, kind[dead], function(x) sprintf("%d x, %.1f MB", length(x), sum(x))))
ns$closeDIAStateAdapter(manager)
desc <- lapply(ents, function(d) readRDS(file.path(d, "descriptor.rds")))
tag <- vapply(desc, function(x) paste(x$kind, x$identity$parameters, x$identity$selectors, sep = "/"), character(1))
cat("DEAD by kind/parameters/selectors:\n"); print(table(tag[dead]))
cat("LIVE by kind/parameters/selectors:\n"); print(table(tag[!dead]))
heads <- grepl("/state-head/", tag)
cat(sprintf("RESULT2 head_snapshots total=%d dead=%d dead_bytes=%.3fGB; other dead=%d bytes=%.3fGB\n", sum(heads), sum(heads & dead), sum(size[heads & dead])/1e9, sum(dead & !heads), sum(size[dead & !heads])/1e9))
