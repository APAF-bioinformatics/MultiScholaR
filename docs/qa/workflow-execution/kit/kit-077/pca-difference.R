# KIT-077: why a QC product's PCA data differs from the baseline hash. On the OLD build (kit_076), the same state and the
# parameters the app used (read from the product): does plotPca give the product's data? Hashes and column names only.
# Usage: R_LIBS=<kit_076 library> Rscript docs/qa/workflow-execution/kit/kit-077/pca-difference.R <member dir> <lipid|metab>
suppressMessages(library(MultiScholaR))
args <- commandArgs(TRUE); member <- args[[1L]]; omic <- args[[2L]]
ns <- asNamespace("MultiScholaR"); api <- function(name) get(name, ns)
h <- new.env(parent = globalenv()); sys.source("tools/profiling/workflow_kit.R", h)
work <- tempfile("kit077-pca-"); dir.create(work); on.exit(unlink(work, recursive = TRUE), add = TRUE)
stopifnot(system2("cp", c("-a", "--reflink=auto", file.path(member, "project"), work)) == 0L)
store_dir <- h$kitFindStore(file.path(work, "project"))
identity <- api("readSavedStoreIdentity")(file.path(store_dir, "store.json"))
manager <- api("newDIAStateAdapter")(api("newWorkflowEntityStore")(dirname(store_dir), identity$project_id, identity$workflow_id), reopen = TRUE)
store <- manager$.context$store
products <- api("findSavedWorkflowEntities")(store, list(kind = "scientific_product"), max_matches = 2000L)$matches
values <- Filter(function(v) is.list(v) && identical(v$figure_type, "pca"), lapply(products, function(r) api("hydrateSavedWorkflowEntity")(store, r)$value))
state_of <- c(post_filter = "_qc_complete", post_norm = "_normalized", ruv_corrected = "_ruv_corrected")
same <- 0L
for (v in values) {
    object <- api("withWorkflowStateScope")(manager, function() manager$getState(paste0(omic, state_of[[v$parameters$stage]])))
    p <- v$parameters
    plots <- suppressWarnings(suppressMessages(plotPca(object, grouping_variable = p$grouping_variable, label_column = "", shape_variable = p$shape_variable, title = "", font_size = p$font_size)))
    match <- identical(digest::digest(as.data.frame(plots[[v$assay]]$data), algo = "sha256"), digest::digest(v$data, algo = "sha256"))
    same <- same + match
    bare <- suppressWarnings(suppressMessages(plotPca(object, grouping_variable = object@group_id, label_column = "", shape_variable = NULL, title = "", font_size = 8)))
    cat(sprintf("%-32s app parameters: grouping %s, shape %s; identical on the old build: %s; columns the bare call lacks: [%s]; extra: [%s]\n", v$figure, format(p$grouping_variable), format(p$shape_variable), match,
        paste(setdiff(names(v$data), names(bare[[v$assay]]$data)), collapse = ", "), paste(setdiff(names(bare[[v$assay]]$data), names(v$data)), collapse = ", ")))
}
cat(sprintf("RESULT %d of %d PCA products identical to plotPca on the old build with the app's parameters\n", same, length(values)))
api("closeDIAStateAdapter")(manager)
