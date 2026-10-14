# KIT-017 screen run, candidate-1 (kit_024): the RUV grid's server-side page request. Shiny serves the widget's page
# requests with no reactive domain and no open scope; the reader resolved the owner at request time and found none.
# Run on the installed library under test against the failed member's project (read-only: the display table save is
# content-addressed and the state is not touched). Usage: Rscript paged-render-probe.R <member dir> <label>
suppressMessages(library(MultiScholaR))
args <- commandArgs(TRUE); member <- args[[1L]]; label <- args[[2L]]
project <- file.path(member, "project", "dia_baseline_smoke")
source_dir <- file.path(project, "scripts", "proteomics")
roots <- list.files(source_dir, "^[.]dia-native-[0-9a-f]+$", all.files = TRUE, include.dirs = TRUE, full.names = TRUE)
stopifnot(length(roots) == 1L)
identity <- MultiScholaR:::readSavedStoreIdentity(file.path(roots, ".workflow-entities", "store.json"))
manager <- MultiScholaR:::newDIAStateAdapter(MultiScholaR:::newWorkflowEntityStore(roots, identity$project_id, identity$workflow_id, max_cache_bytes = 0), reopen = TRUE)
manager$.context$session_format <- "native_v1"; manager$.context$reference_states <- TRUE
reference <- MultiScholaR:::withDIAStateScope(manager, function() {
    result <- MultiScholaR:::resolveDIARuvOptimizationResult(manager)
    MultiScholaR:::saveDIADisplayTable(manager, result$optimization_results)
})
stopifnot(MultiScholaR:::isWorkflowTableReference(reference))
outside <- tryCatch(list(ok = TRUE, rows = nrow(MultiScholaR:::readDIADisplayTablePage(reference, limit = 10L)$data)),
    error = function(e) list(ok = FALSE, class = class(e)[[1L]], message = conditionMessage(e)))
explicit <- tryCatch(list(ok = TRUE, rows = nrow(MultiScholaR:::readDIADisplayTablePage(reference, limit = 10L, manager = manager)$data)),
    error = function(e) list(ok = FALSE, class = class(e)[[1L]], message = conditionMessage(e)))
record <- list(schema = "kit.paged_render_probe.v1", ticket = "KIT-017", label = label, date = format(Sys.Date()),
    library = dirname(system.file(package = "MultiScholaR")), build = MultiScholaR:::workflowExecutionBuildIdentity(),
    rows = reference$metadata$shape[[1L]], request_without_owner = outside, request_with_render_owner = explicit)
jsonlite::write_json(record, paste0("paged-render-probe-", label, ".json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
cat(jsonlite::toJSON(record[c("label", "rows", "request_without_owner", "request_with_render_owner")], auto_unbox = TRUE), "\n")
