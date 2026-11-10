# KIT-113 W1 probe P4 (R3) at W0's 16,000 x 300 x 3 (scratch, not frozen). Run on the start commit's R code:
#   Rscript docs/qa/workflow-execution/kit/kit-113/w1/probe-p4.R <package source>
# W0's fixture (`buildCase` and `newBus`, read from kit-113/w0/enrichment-heap.R's definitions: a synthetic DIA owner,
# the in-process DA and its publication), then the enrichment's request (probe-lib.R's request half, the annotation
# table as an `rds` transfer value) and the worker's value (the emulated handler through the worker's entry point and
# the controller's acceptance), with the live `processEnrichments` over the e2e doubles and with the deterministic
# double. Prints the request's parameters' and the value's sizes against the limits (65,536 and 8,388,608 bytes,
# measured by `object.size` as `validateWorkflowActionParameters` and `acceptWorkflowActionResult` measure them).
# Network: the lane's closed-port proxy; the doubles replace gprofiler2::gost and clusterProfiler::enricher.

args <- commandArgs(TRUE)
tree <- normalizePath(args[[1L]])
started <- proc.time()[["elapsed"]]
Sys.setenv(http_proxy = "http://127.0.0.1:9", https_proxy = "http://127.0.0.1:9", HTTP_PROXY = "http://127.0.0.1:9",
    HTTPS_PROXY = "http://127.0.0.1:9", no_proxy = "", NO_PROXY = "")
setwd(tree)
suppressMessages(pkgload::load_all(tree, quiet = TRUE, export_all = TRUE, helpers = TRUE))
suppressPackageStartupMessages(library(testthat))
source(file.path(tree, "docs", "qa", "workflow-execution", "kit", "kit-113", "w1", "probe-lib.R"))
loadKit113TestHelpers(tree)
install_enrichment_doubles(.local_envir = globalenv())
kit113Mocks(globalenv())
registerKit113EnrichmentEntry()
mib <- function(bytes) bytes / 2^20
for (expression in parse(file.path(tree, "docs", "qa", "workflow-execution", "kit", "kit-113", "w0", "enrichment-heap.R"),
    keep.source = FALSE)) {
    if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
        deparse(expression[[2L]]) %in% c("buildCase", "newBus")) eval(expression, globalenv())
}
ns <- asNamespace("MultiScholaR")

root <- normalizePath(tempfile("kit113-p4-"), mustWork = FALSE)
dir.create(root)
case <- buildCase(16000L, 300L, 3L, root)
say("P4 setup done at %.0f s", proc.time()[["elapsed"]] - started)
controller <- omicsWorkerController(list(owner = case$manager))
for (kernel in c("processEnrichments", "runProtEnrichDeterministicProcessEnrichments")) {
    workflow <- newBus(case)
    config <- omicsEnrichCoreConfig()
    config$kernel$select <- local({
        name <- kernel
        function() get(name, envir = ns)
    })
    parameters <- list(selected_contrast = case$contrasts$friendly_names[[1L]], organism_taxid = "9606", up_cutoff = 0,
        down_cutoff = 0, q_cutoff = 0.05, correction_method = "gSCS", enable_organism_filter = FALSE)
    t0 <- proc.time()[["elapsed"]]
    half <- kit113Request(workflow, config, case$manifest, parameters, case$paths, NULL, NULL, codec = "rds")
    sizes <- kit113Sizes(half$parameters)
    say("P4 16000 x 300 x 3, kernel %s: request parameters serialized %d B, object.size %d B (limit 65536): %s",
        kernel, sizes[["serialized"]], sizes[["object_size"]], if (sizes[["object_size"]] <= 65536) "within" else "OVER")
    for (field in names(half$parameters)) say("    parameters$%-14s object.size %8d B", field,
        as.integer(object.size(half$parameters[[field]])))
    context <- list(state_name = config$state_name, config = list(), description = "enrichment")
    dispatched <- withWorkflowStateScope(case$manager, function() kit113Dispatch(case$manager, controller, case$paths,
        "enrichment", half$parameters, case$manifest$source, context))
    value <- dispatched$message$result
    sizes <- kit113Sizes(value)
    say("P4 16000 x 300 x 3, kernel %s: value serialized %d B, object.size %d B (limit 8388608): %s; acceptance %s; %.0f s",
        kernel, sizes[["serialized"]], sizes[["object_size"]], if (sizes[["object_size"]] <= 8388608) "within" else "OVER",
        dispatched$receipt$status, proc.time()[["elapsed"]] - t0)
    for (field in names(value)) say("    value$%-12s object.size %8d B", field, as.integer(object.size(value[[field]])))
    for (field in names(value$display)) say("    value$display$%-10s object.size %8d B", field,
        as.integer(object.size(value$display[[field]])))
    results <- resolveWorkflowStateValue(case$manager$.context, value$display$results)
    say("    the results' rds state value: %.2f MiB on disk (%s)", mib(value$display$results$payload$bytes),
        value$display$results$payload$codec)
    rm(results)
    invisible(gc())
}
say("P4 total %.0f s", proc.time()[["elapsed"]] - started)
closeDIAStateAdapter(case$manager)
unlink(root, recursive = TRUE)
