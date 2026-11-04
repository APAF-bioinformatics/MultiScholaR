# KIT-082 item 8: what `readWorkflowOmicsDaResults` keeps on the DA tab, on a screen journey member, beside the assay it
# comes from. Sizes only. Usage: R_LIBS=<library> Rscript docs/qa/workflow-execution/kit/kit-082/da-tab-read-size.R <member dir>...
suppressMessages(library(MultiScholaR)); ns <- asNamespace("MultiScholaR"); api <- function(name) get(name, ns)
h <- new.env(parent = globalenv()); sys.source("tools/profiling/workflow_kit.R", h)
mb <- function(x) round(as.numeric(utils::object.size(x)) / 2^20, 2)
for (member in commandArgs(TRUE)) {
    work <- tempfile("kit082-da-"); dir.create(work)
    stopifnot(system2("cp", c("-a", "--reflink=auto", file.path(member, "project"), work)) == 0L)
    store_dir <- h$kitFindStore(file.path(work, "project")); identity <- api("readSavedStoreIdentity")(file.path(store_dir, "store.json"))
    manager <- api("newDIAStateAdapter")(api("newWorkflowEntityStore")(dirname(store_dir), identity$project_id, identity$workflow_id), reopen = TRUE)
    store <- manager$.context$store
    products <- api("findSavedWorkflowEntities")(store, list(kind = "scientific_product"), max_matches = 2000L)$matches
    values <- Filter(function(v) is.list(v) && is.list(v$assays) && all(vapply(v$assays, function(a) is.list(a) && !is.null(a$table), TRUE)),
        lapply(products, function(r) tryCatch(api("hydrateSavedWorkflowEntity")(store, r)$value, error = function(e) NULL)))
    cat(sprintf("%s: %s, %d contrast products, %.2f MiB together\n", basename(dirname(member)), manager$.context$omic_type, length(values), mb(values)))
    for (v in values) for (a in names(v$assays)) cat(sprintf("  %s: table %d rows x %d columns, %.2f MiB; whole entry %.2f MiB\n", a,
        nrow(v$assays[[a]]$table), ncol(v$assays[[a]]$table), mb(v$assays[[a]]$table), mb(v$assays[[a]])))
    api("closeDIAStateAdapter")(manager); unlink(work, recursive = TRUE)
}
