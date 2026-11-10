# KIT-113 W1 probe P5 (R5, R6; scratch, not frozen). Run on the start commit's R code:
#   Rscript docs/qa/workflow-execution/kit/kit-113/w1/probe-p5.R <package source>
# After a managed DA (KIT-112 W3's fixture, `omicDaK112ProtSetup(managed = TRUE)` from test-omic-da-module.R: MCI-009's
# object on an on-disk DIA owner, the DA tab's run dispatched through the omics runner, `DA-contrast-output` then
# `da_publication`, through the worker's entry point and the controller's acceptance in this process), the enrichment's
# real request (probe-lib.R's request half on the bus and the manifest the DA left) goes through `runWorkflowAction`'s
# guards in its order, then the worker's entry point, the acceptance (its value checks) and the omics completion;
# first Run with the head as input, second Run with the manifest's source as input while the head is
# `enrichment_completed`. The start commit's registry entry is tried first; then D3's entry is registered in this
# process only.

args <- commandArgs(TRUE)
tree <- normalizePath(args[[1L]])
setwd(tree)
suppressMessages(pkgload::load_all(tree, quiet = TRUE, export_all = TRUE, helpers = TRUE))
library(testthat)
source(file.path(tree, "docs", "qa", "workflow-execution", "kit", "kit-113", "w1", "probe-lib.R"))
loadKit113TestHelpers(tree)
for (expression in parse(file.path(tree, "tests", "testthat", "test-omic-da-module.R"), keep.source = FALSE)) {
    if (is.call(expression) && identical(expression[[1L]], as.name("<-"))) eval(expression, globalenv())
}

report <- new.env()
enrich <- function(n, setup) function(ctx) {
    withr::local_envvar(MULTISCHOLAR_TEST_MODE = "true", http_proxy = "http://127.0.0.1:9", https_proxy = "http://127.0.0.1:9",
        HTTP_PROXY = "http://127.0.0.1:9", HTTPS_PROXY = "http://127.0.0.1:9", no_proxy = "", NO_PROXY = "")
    shiny::isolate({
        bus <- ctx$bus()
        manager <- bus$state_manager
        controller <- manager$.context$store$execution_owner
        for (field in c("enrichment_analysis_results", "enrichment_ui_params", "enrichment_organism_filter",
            "mixed_species_analysis", "taxon_id")) if (!field %in% names(bus)) bus[[field]] <- NULL
        paths <- c(manager$.context$experiment_paths, list(results_dir = file.path(setup$roots, "results"),
            pathway_dir = file.path(setup$roots, "results", "pathway_enrichment")))
        dir.create(paths$pathway_dir, recursive = TRUE, showWarnings = FALSE)
        config <- omicsEnrichCoreConfig()
        manifest <- bus$da_analysis_results_list
        inputs <- omicsEnrichCoreInputs("9606")
        snapshot_before <- withWorkflowStateScope(manager, function() manager$.context$snapshot)
        head <- getWorkflowStateReference(manager)
        out <- list(head = paste(tail(unlist(snapshot_before$state_history), 1L)),
            input_is_head = identical(manifest$source, head),
            manifest_source_is_head_entity = identical(manifest$source$entity_id, head$entity_id),
            history_before = paste(unlist(snapshot_before$state_history), collapse = ","),
            manifest_is_head_snapshot = identical(resolveWorkflowStateValue(manager$.context, snapshot_before$da_results), manifest))
        half <- kit113Request(bus, config, manifest, inputs, paths, NULL, NULL, codec = "rds")
        input <- list(state = manifest$source)
        context <- list(state_name = config$state_name, config = list(), description = "enrichment")
        if (n == 1L) out$guards_start_registry <- kit113Guards(manager, "enrichment", input, half$parameters, context)
        registerKit113EnrichmentEntry()
        out$guards <- kit113Guards(manager, "enrichment", input, half$parameters, context)
        out$sizes <- kit113Sizes(half$parameters)
        dispatched <- tryCatch(kit113Dispatch(manager, controller, paths, "enrichment", half$parameters, manifest$source,
            context), error = function(e) e)
        if (inherits(dispatched, "error")) {
            out$dispatch_error <- conditionMessage(dispatched)
            report[[paste0("run", n)]] <- out
            return(invisible())
        }
        value <- dispatched$message$result
        snapshot_after <- withWorkflowStateScope(manager, function() manager$.context$snapshot)
        out$value_checks <- kit113ValueChecks(manager$.context$store, value)
        out$value_sizes <- kit113Sizes(value)
        out$receipt <- dispatched$receipt$status
        out$completion <- c(
            history_unchanged = identical(snapshot_after$state_history, snapshot_before$state_history),
            head_state_unchanged = identical(getWorkflowStateReference(manager), head),
            product_rooted = value$reference$entity_id %in% names(snapshot_after$products),
            product_roots = length(snapshot_after$product_dependencies[[value$reference$entity_id]]),
            action_outputs_unchanged = identical(snapshot_after$action_outputs, snapshot_before$action_outputs),
            external_resources_unchanged = identical(snapshot_after$external_resources, snapshot_before$external_resources),
            staged_publication_files = length(list.files(file.path(dispatched$export$request_dir, "omics-publication"),
                recursive = TRUE)),
            pending_provenance = if (is.null(manager$.context$pending_provenance)) "<none>" else
                paste(manager$.context$pending_provenance$operation, "values:",
                    paste(names(manager$.context$pending_provenance$values), collapse = "|")))
        accepted <- tryCatch(kit113Accept(bus, config, half$context, value, inputs, NULL), error = function(e) e)
        if (inherits(accepted, "error")) out$accept_error <- conditionMessage(accepted) else {
            snapshot_done <- withWorkflowStateScope(manager, function() manager$.context$snapshot)
            out$accepted <- c(saved = isTRUE(accepted$saved$saved),
                history = paste(unlist(snapshot_done$state_history), collapse = ","),
                provenance_operation = snapshot_done$step_provenance[["enrichment_completed"]]$operation %||% "<none>",
                provenance_values = paste(names(snapshot_done$step_provenance[["enrichment_completed"]]$values), collapse = "|"),
                pending_after = if (is.null(manager$.context$pending_provenance)) "<none>" else "set",
                resources = length(snapshot_done$external_resources))
        }
        report[[paste0("run", n)]] <- out
    })
}

kit113Mocks(globalenv())
drive <- omicDaK112Drive(function(setup) c(omicDaK112Steady(setup), list(enrich1 = enrich(1L, setup),
    enrich2 = enrich(2L, setup))), managed = TRUE)
say("P5 the managed DA: operations dispatched %s", paste(drive$counts$operations, collapse = ", "))
for (step in names(drive$captures)) {
    errors <- drive$captures[[step]]$errors
    if (length(errors)) say("P5 step %s errors: %s", step, paste(errors, collapse = " | "))
}
for (n in c("run1", "run2")) {
    out <- report[[n]]
    say("")
    say("== P5 %s: head %s; history %s; input = manifest$source; input is the head: %s; the held manifest is the head snapshot's: %s",
        n, out$head, out$history_before, out$input_is_head, out$manifest_is_head_snapshot)
    if (!is.null(out$guards_start_registry)) {
        say("  the start commit's registry entry (in_process):")
        for (g in names(out$guards_start_registry)) say("    %-58s %s", g, out$guards_start_registry[[g]])
    }
    say("  D3's entry (omics, named_input):")
    for (g in names(out$guards)) say("    %-58s %s", g, out$guards[[g]])
    say("  request parameters: serialized %d B, object.size %d B (limit 65536)", out$sizes[["serialized"]],
        out$sizes[["object_size"]])
    if (!is.null(out$dispatch_error)) { say("  DISPATCH ERROR: %s", out$dispatch_error); next }
    say("  acceptance receipt: %s", out$receipt)
    for (g in names(out$value_checks)) say("    value %-52s %s", g, out$value_checks[[g]])
    say("  value: serialized %d B, object.size %d B (limit 8388608)", out$value_sizes[["serialized"]], out$value_sizes[["object_size"]])
    for (g in names(out$completion)) say("    completion %-32s %s", g, out$completion[[g]])
    if (!is.null(out$accept_error)) say("  ACCEPT ERROR: %s", out$accept_error) else
        for (g in names(out$accepted)) say("    accept half %-31s %s", g, out$accepted[[g]])
}
