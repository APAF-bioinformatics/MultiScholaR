# KIT-082 item 8: what the DA tab held from `readWorkflowOmicsDaResults`, and what it holds now, on a screen journey member. Sizes
# only. An earlier version measured other products through a partial match on `$table`; its 0.54 and 0.80 MiB are wrong. Usage: R_LIBS=<library> Rscript docs/qa/workflow-execution/kit/kit-082/da-tab-read-size.R <member dir>...
suppressMessages(library(MultiScholaR)); ns <- asNamespace("MultiScholaR"); api <- function(name) get(name, ns)
h <- new.env(parent = globalenv()); sys.source("tools/profiling/workflow_kit.R", h)
mb <- function(x) round(as.numeric(utils::object.size(x)) / 2^20, 2)
for (member in commandArgs(TRUE)) {
    work <- tempfile("kit082-da-"); dir.create(work)
    stopifnot(system2("cp", c("-a", "--reflink=auto", file.path(member, "project"), work)) == 0L)
    store_dir <- h$kitFindStore(file.path(work, "project")); identity <- api("readSavedStoreIdentity")(file.path(store_dir, "store.json"))
    manager <- api("newDIAStateAdapter")(api("newWorkflowEntityStore")(dirname(store_dir), identity$project_id, identity$workflow_id), reopen = TRUE)
    store <- manager$.context$store
    binding <- api("omicsDeclarationBindings")()[[manager$.context$omic_type]]
    manifest <- api("resolveWorkflowStateValue")(manager$.context, manager$.context$snapshot$da_results)
    # What the tab held before revision 1.2, and what it holds since.
    whole <- api("readWorkflowOmicsDaResults")(binding, manager, manifest)
    tab <- api("omicsDaTabResults")(binding, whole)
    rows <- whole[[binding$long_table]]
    cat(sprintf("%s: %s, %d contrasts, %d rows x %d columns; whole read %.2f MiB; the tab's rows %.2f MiB, %d columns\n", basename(dirname(member)),
        manager$.context$omic_type, length(manifest$results), nrow(rows), ncol(rows), mb(whole), mb(tab), ncol(tab[[binding$long_table]])))
    api("closeDIAStateAdapter")(manager); unlink(work, recursive = TRUE)
}
