# KIT-113 W1 probes P1 (the portable results), P2, P4 (lane fixtures) and P5's guards on W1 (b)'s fixture (scratch,
# not frozen). Run on the start commit's R code:
#   Rscript docs/qa/workflow-execution/kit/kit-113/w1/probe-managed.R <package source> <transfer codec: rds|table> [cases]
# Each case drives two sessions on W1 (b)'s fixture (test-omic-enrich-module.R: the DA through the module, the tab,
# the contrast): in process, `runWorkflowEnrichment` twice (as test-omics-enrich-core.R runs it); managed, D3's halves
# emulated by probe-lib.R twice, through the worker's entry point and the controller's acceptance in this process
# (the in-process test worker, `omicsWorkerController`), the first Run with the head as input, the second with the
# manifest's source as input while the head is `enrichment_completed`. Prints what it finds.

args <- commandArgs(TRUE)
tree <- normalizePath(args[[1L]])
codec <- if (length(args) > 1L) args[[2L]] else "rds"
cases <- if (length(args) > 2L) strsplit(args[[3L]], ",")[[1L]] else
    c("gprofiler", "clusterprofiler", "kernel_gprofiler", "kernel_clusterprofiler")
setwd(tree)
suppressMessages(pkgload::load_all(tree, quiet = TRUE, export_all = TRUE, helpers = TRUE))
library(testthat)
source(file.path(tree, "docs", "qa", "workflow-execution", "kit", "kit-113", "w1", "probe-lib.R"))
loadKit113TestHelpers(tree)
kit113Mocks(globalenv())
start <- registerKit113EnrichmentEntry()
say("registry (this process only): enrichment compute=%s handler=%s completion=%s named_input=%s", start$compute,
    start$handler, start$completion, start$named_input)

results_fields <- function(x) list(enrichment_data = x@enrichment_data, enrichment_summaries = x@enrichment_summaries,
    contrasts = x@contrasts, plots = lapply(x@enrichment_plots, function(p) lapply(p, class)))

#' One session: W1 (b)'s steps to `select`, then two runs by `run(ctx, inputs, n)`, each followed by its captures.
drive <- function(case, run) {
    spec <- omicsEnrichCoreCases[[case]]
    w1 <- omicEnrichW1CaseSpec(case, environment())[[1L]]
    inputs <- omicsEnrichCoreInputs(spec$taxon)
    before <- intersect(c("start", "da", "open", "select"), names(omicEnrichW1Frozen[[case]]))
    step <- function(n) function(ctx) {
        run(ctx, inputs, n)
        ctx$keep("d6", omicEnrichW113Fields(ctx))
        ctx$keep("pathway", kit113PathwayFiles(ctx$paths$pathway_dir, ctx$roots))
    }
    omicEnrichW1Drive(w1$setup, c(w1$steps[before], list(run = step(1L), again = step(2L))), live = spec$live)
}

notify <- function(message, type, duration) shiny::showNotification(message, type = type, duration = duration)
progress <- function(value, detail) shiny::incProgress(value, detail = detail)

for (case in cases) {
    config <- omicsEnrichCoreConfig()
    held <- new.env()
    in_process <- drive(case, function(ctx, inputs, n) shiny::isolate({
        bus <- ctx$bus()
        out <- runWorkflowEnrichment(bus, config, bus$da_analysis_results_list, inputs, ctx$paths, notify = notify,
            progress = progress, refuse_missing_da = TRUE)
        held[[paste0("out", n)]] <- out
        ctx$keep("fields", results_fields(out$results))
        ctx$keep("state", omicEnrichW1Digest(withWorkflowStateScope(bus$state_manager, function()
            bus$state_manager$getState("enrichment_completed")), ctx$roots))
    }))
    managed <- drive(case, function(ctx, inputs, n) shiny::isolate({
        bus <- ctx$bus()
        manager <- bus$state_manager
        manifest <- bus$da_analysis_results_list
        if (is.null(held$controller)) held$controller <- omicsWorkerController(list(owner = manager))
        half <- kit113Request(bus, config, manifest, inputs, ctx$paths, notify, progress, codec = codec)
        input <- list(state = manifest$source)
        context <- list(state_name = config$state_name, config = list(), description = "enrichment")
        head_before <- getWorkflowStateReference(manager)
        guards <- kit113Guards(manager, "enrichment", input, half$parameters, context)
        snapshot_before <- manager$.context$snapshot
        dispatched <- kit113Dispatch(manager, held$controller, ctx$paths, "enrichment", half$parameters,
            manifest$source, context)
        snapshot_after <- withWorkflowStateScope(manager, function() manager$.context$snapshot)
        value <- dispatched$message$result
        completion <- list(
            receipt = dispatched$receipt$status,
            history_unchanged = identical(snapshot_after$state_history, snapshot_before$state_history),
            head_state_unchanged = identical(getWorkflowStateReference(manager), head_before),
            product_rooted = value$reference$entity_id %in% names(snapshot_after$products),
            product_roots = length(snapshot_after$product_dependencies[[value$reference$entity_id]]),
            action_outputs_unchanged = identical(snapshot_after$action_outputs, snapshot_before$action_outputs),
            staged = length(list.files(file.path(dispatched$export$request_dir, "omics-publication"), recursive = TRUE)),
            pending_provenance = if (is.null(manager$.context$pending_provenance)) "<none>" else
                paste(manager$.context$pending_provenance$operation, paste(names(manager$.context$pending_provenance$values),
                    collapse = "|")))
        out <- kit113Accept(bus, config, half$context, value, inputs, progress)
        held[[paste0("out", n)]] <- out
        ctx$keep("guards", guards)
        ctx$keep("checks", kit113ValueChecks(manager$.context$store, value))
        ctx$keep("completion", completion)
        ctx$keep("sizes", list(parameters = kit113Sizes(half$parameters), value = kit113Sizes(value)))
        ctx$keep("fields", results_fields(out$results))
        ctx$keep("state", omicEnrichW1Digest(withWorkflowStateScope(manager, function()
            manager$getState("enrichment_completed")), ctx$roots))
        ctx$keep("head_at_dispatch", if (n == 1L) "head = the manifest's source" else
            paste("head =", paste(tail(unlist(snapshot_before$state_history), 1L))))
    }))
    say("")
    say("== %s (transfer codec %s) ==", case, codec)
    for (n in c("run", "again")) {
        a <- in_process[[n]]
        b <- managed[[n]]
        say("-- %s: managed errors: %s; in-process errors: %s", n, paste(b$errors, collapse = " | "), paste(a$errors, collapse = " | "))
        say("   head at dispatch: %s", b$head_at_dispatch)
        for (guard in names(b$guards)) say("   guard %-58s %s", guard, b$guards[[guard]])
        for (check in names(b$checks)) say("   value %-58s %s", check, b$checks[[check]])
        for (field in names(b$completion)) say("   completion %-30s %s", field, format(b$completion[[field]]))
        say("   sizes: parameters serialized %d B, object.size %d B (limit 65536); value serialized %d B, object.size %d B (limit 8388608)",
            b$sizes$parameters[["serialized"]], b$sizes$parameters[["object_size"]], b$sizes$value[["serialized"]],
            b$sizes$value[["object_size"]])
        # P1 (R1): the results the controller hydrated against the live object, field by field; the committed state.
        for (field in names(a$fields)) say("   P1 results$%-22s identical: %s", field, identical(a$fields[[field]], b$fields[[field]]))
        say("   P1 committed state digest: in process %s, managed %s, identical %s", a$state, b$state, identical(a$state, b$state))
        if (identical(n, "run")) {
            frozen <- omicEnrichW1Frozen[[case]]$run$states[["enrichment_completed"]]
            say("   P1 committed state digest = W1 (b)'s frozen run step (%s): in process %s, managed %s", frozen,
                identical(a$state, frozen), identical(b$state, frozen))
        }
        # P2: the pathway directory's files.
        say("   P2 pathway files: in process %d, managed %d; names identical %s; digests identical %s", length(a$pathway),
            length(b$pathway), identical(names(a$pathway), names(b$pathway)), identical(a$pathway, b$pathway))
        if (!identical(a$pathway, b$pathway)) {
            for (file in union(names(a$pathway), names(b$pathway))) if (!identical(a$pathway[file], b$pathway[file]))
                say("      DIFF %s: in process %s; managed %s", file, a$pathway[file], b$pathway[file])
        }
        # D6's other fields: the run step's (core) fields and the KIT-113 W1 block's fields.
        for (field in setdiff(omicsEnrichCoreFields, c("files", "errors"))) {
            same <- identical(a[[field]], b[[field]])
            say("   D6 %-12s identical: %s", field, same)
            if (!same) {
                say("      in process: %s", paste(utils::capture.output(str(a[[field]])), collapse = " / "))
                say("      managed:    %s", paste(utils::capture.output(str(b[[field]])), collapse = " / "))
            }
        }
        files <- function(x) if (is.character(x)) x[!startsWith(names(x), "<store entries>")] else x
        say("   D6 files (store entries excluded) identical: %s", identical(files(a$files), files(b$files)))
        store <- function(x) if (is.character(x)) x[startsWith(names(x), "<store entries>")] else x
        if (!identical(store(a$files), store(b$files))) {
            sa <- store(a$files); sb <- store(b$files)
            for (key in union(names(sa), names(sb))) if (!identical(sa[key], sb[key]))
                say("      store entries %s: in process %s, managed %s", sub("^<store entries>", "", key), sa[key], sb[key])
        }
        for (field in names(b$d6)) {
            same <- identical(a$d6[[field]], b$d6[[field]])
            say("   D6 (W1 block) %-10s identical: %s", field, same)
            if (!same) {
                say("      in process: %s", paste(utils::capture.output(str(a$d6[[field]])), collapse = " / "))
                say("      managed:    %s", paste(utils::capture.output(str(b$d6[[field]])), collapse = " / "))
            }
        }
    }
}
say("")
say("done")
