# KIT-077 acceptance on a journey member: every figure file against the baseline member's paths; every figure file's bytes
# against the store's rendered assets (is it really written from the store?); every asset rooted; each QC product's plot
# data against the baseline hash. Paths, counts and hashes only.
# Usage: R_LIBS=<library> Rscript docs/qa/workflow-execution/kit/kit-077/figure-compare.R <member dir> <baseline rds>
suppressMessages(library(MultiScholaR))
args <- commandArgs(TRUE); member <- args[[1L]]; baseline <- readRDS(args[[2L]])
ns <- asNamespace("MultiScholaR"); api <- function(name) get(name, ns)
h <- new.env(parent = globalenv()); sys.source("tools/profiling/workflow_kit.R", h)
work <- tempfile("kit077-compare-"); dir.create(work); on.exit(unlink(work, recursive = TRUE), add = TRUE)
stopifnot(system2("cp", c("-a", "--reflink=auto", file.path(member, "project"), work)) == 0L)
project <- file.path(work, "project")
store_dir <- h$kitFindStore(project)
identity <- api("readSavedStoreIdentity")(file.path(store_dir, "store.json"))
manager <- api("newDIAStateAdapter")(api("newWorkflowEntityStore")(dirname(store_dir), identity$project_id, identity$workflow_id), reopen = TRUE)
store <- manager$.context$store
files <- list.files(project, pattern = "\\.(png|pdf|html|svg)$", recursive = TRUE)
files <- files[!grepl("^\\.|/\\.", files)]
missing <- setdiff(baseline$files$path, files)
cat(sprintf("files: %d now, %d in the baseline; missing from the baseline's paths: %d; new: %d; empty: %d\n", length(files), nrow(baseline$files), length(missing),
    length(setdiff(files, baseline$files$path)), sum(file.size(file.path(project, files)) == 0)))
for (m in utils::head(missing, 20L)) cat("  MISSING", m, "\n")
new_kinds <- sub("^.*/", "", setdiff(files, baseline$files$path))
cat("  new files, by ending:", paste(names(table(sub("^.*_", "", new_kinds))), table(sub("^.*_", "", new_kinds)), collapse = "; "), "\n")
assets <- api("findSavedWorkflowEntities")(store, list(kind = "rendered_asset"), max_matches = 2000L)$matches
asset_hash <- vapply(assets, function(r) digest::digest(api("hydrateSavedWorkflowEntity")(store, r)$bytes, serialize = FALSE, algo = "sha256"), "")
roots <- names(manager$.context$snapshot$products)
asset_ids <- vapply(assets, function(r) r$entity_id, "")
cat(sprintf("assets: %d; rooted in the head: %d\n", length(assets), sum(asset_ids %in% roots)))
file_hash <- vapply(file.path(project, files), function(f) digest::digest(file = f, algo = "sha256"), "", USE.NAMES = FALSE)
outside <- files[!file_hash %in% asset_hash]
cat(sprintf("figure files whose bytes are a stored asset: %d of %d; not stored: %d\n", sum(file_hash %in% asset_hash), length(files), length(outside)))
for (f in outside) cat("  NOT STORED", f, "\n")
products <- api("findSavedWorkflowEntities")(store, list(kind = "scientific_product"), max_matches = 2000L)$matches
values <- Filter(function(v) is.list(v) && !is.null(v$figure_type), lapply(products, function(r) api("hydrateSavedWorkflowEntity")(store, r)$value))
cat("products by type:", paste(names(table(vapply(values, `[[`, "", "figure_type"))), table(vapply(values, `[[`, "", "figure_type")), collapse = "; "), "\n")
stage_state <- function(stage) paste0(baseline$omic, c(post_filter = "_qc_complete", post_norm = "_normalized", ruv_corrected = "_ruv_corrected")[[stage]])
compared <- 0L; different <- character()
for (v in values) if (identical(v$schema, "omics_qc_plot_product_v1")) {
    expected <- baseline$plot_data[[stage_state(v$parameters$stage)]][[v$figure_type]][[v$assay]]
    if (is.null(expected)) { different <- c(different, paste("no baseline:", v$figure)); next }
    compared <- compared + 1L
    if (!identical(digest::digest(v$data, algo = "sha256"), expected)) different <- c(different, v$figure)
}
cat(sprintf("plot data compared with the baseline: %d; different: %d\n", compared, length(different)))
for (d in different) cat("  DIFFERENT", d, "\n")
api("closeDIAStateAdapter")(manager)
cat(if (length(missing) || length(different)) "RESULT differs\n" else "RESULT every baseline file present; plot data equal\n")
