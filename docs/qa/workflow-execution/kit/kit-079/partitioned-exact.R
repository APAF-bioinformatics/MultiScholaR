# KIT-079 scientific contract on a real screen journey: each partitioned variant, run on a copy of a finished `whole` member's
# project from the state before the step with the step's recorded parameters, against the state that member committed: every
# assay table's values and order, the design, the recorded arguments, the dimensions and the step's value. Hashes only.
# Usage: R_LIBS=<library> Rscript docs/qa/workflow-execution/kit/kit-079/partitioned-exact.R <member dir> <lipid|metab>
suppressMessages(library(MultiScholaR))
args <- commandArgs(TRUE); member <- args[[1L]]; prefix <- args[[2L]]
ns <- asNamespace("MultiScholaR"); api <- function(name) get(name, ns)
h <- new.env(parent = globalenv()); sys.source("tools/profiling/workflow_kit.R", h)
`%||%` <- function(a, b) if (is.null(a)) b else a
work <- tempfile("kit079-exact-"); dir.create(work); on.exit(unlink(work, recursive = TRUE), add = TRUE)
stopifnot(system2("cp", c("-a", "--reflink=auto", file.path(member, "project"), work)) == 0L)
store_dir <- h$kitFindStore(file.path(work, "project"))
identity <- api("readSavedStoreIdentity")(file.path(store_dir, "store.json"))
manager <- api("newDIAStateAdapter")(api("newWorkflowEntityStore")(dirname(store_dir), identity$project_id, identity$workflow_id), reopen = TRUE)
store <- manager$.context$store
snapshot <- manager$.context$snapshot
history <- unlist(snapshot$state_history)
binding <- api("omicsDeclarationBindings")()[[if (prefix == "lipid") "lipidomics" else "metabolomics"]]
hash_tables <- function(tables) lapply(tables, function(table) digest::digest(as.data.frame(api("materializeWorkflowTable")(manager$.context, table)), algo = "sha256"))
same <- TRUE
for (step in c("intensity_filter", "duplicate_resolution", "log2")) {
  operation <- paste0(prefix, "_", step)
  state <- Filter(function(name) identical(snapshot$step_provenance[[name]]$operation, operation), names(snapshot$step_provenance))[[1L]]
  source <- snapshot$states[[history[[match(state, history) - 1L]]]]
  entry <- api("workflowOperationEntry")(operation)
  committed <- api("readWorkflowStateMetadata")(manager, snapshot$states[[state]])$data
  parameters <- snapshot$step_provenance[[state]]$values %||% list()
  started <- proc.time()[["elapsed"]]
  value <- suppressWarnings(suppressMessages(api("withWorkflowStateScope")(manager, function() {
    identity_now <- api("workflowComponentScientificIdentity")(manager, source, operation, parameters[entry$declaration$parameters])
    operation_store <- api("newWorkflowOperationStore")(manager, entry$declaration, source, identity_now)
    entry$declaration$partitioned(parameters, operation_store)
  })))
  seconds <- proc.time()[["elapsed"]] - started
  result <- api("readWorkflowStateMetadata")(manager, value$final)$data
  whole_value <- api("hydrateSavedWorkflowEntity")(store, committed$products[[paste0(operation, "_science")]])$value
  checks <- c(assays = identical(hash_tables(result$components[[binding$data_slot]]), hash_tables(committed$components[[binding$data_slot]])),
    design = identical(hash_tables(result$components["design_matrix"]), hash_tables(committed$components["design_matrix"])),
    arguments = identical(result$slots$args, committed$slots$args), dimensions = identical(result$dimensions, committed$dimensions),
    value = identical(value[setdiff(names(value), "final")], whole_value[setdiff(names(whole_value), "final")]))
  rows <- vapply(committed$components[[binding$data_slot]], function(table) as.integer(table$metadata$shape[[1L]]), integer(1))
  cat(sprintf("%s: %s; rows after %s; %.1f s\n", operation, if (all(checks)) "EXACT" else paste("DIFFERS:", paste(names(checks)[!checks], collapse = ", ")),
    paste(rows, collapse = "+"), seconds))
  same <- same && all(checks)
}
api("closeDIAStateAdapter")(manager)
cat(if (same) "RESULT exact\n" else "RESULT differs\n")
