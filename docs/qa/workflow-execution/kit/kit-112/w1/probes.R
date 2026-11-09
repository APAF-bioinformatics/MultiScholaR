# KIT-112 W1 probes P1, P2 and P4 (scratch, not frozen). Run on the start commit's R code:
#   Rscript docs/qa/workflow-execution/kit/kit-112/w1/probes.R <package source> <P1|P2|P4>
# The owners are W1 (a)'s builders (`savedDaW1Owner`, `savedDaW1LimpaOwner`), read from test-workflow-saved-da.R's
# function definitions only. Prints what it finds; `kit-112/w1/probes.txt` records the outcome.

args <- commandArgs(TRUE)
tree <- normalizePath(args[[1L]])
probe <- args[[2L]]
setwd(tree)
suppressMessages(pkgload::load_all(tree, quiet = TRUE, export_all = TRUE, helpers = TRUE))
library(testthat)
for (expression in parse(file.path(tree, "tests", "testthat", "test-workflow-saved-da.R"))) {
    if (is.call(expression) && identical(expression[[1L]], as.name("<-")) && grepl("^savedDaW1", deparse(expression[[2L]])))
        eval(expression, globalenv())
}
quiet <- function(code) {
    value <- NULL
    utils::capture.output(value <- suppressWarnings(suppressMessages(code)))
    value
}
say <- function(...) cat(sprintf(...), "\n", sep = "")

#' The in-process parameters of one contrast, as `runDIAReferenceDaAnalysis` builds them (utils_workflow_dia_da.R:397-399).
contrastParameters <- function(manager, source, contrasts, i, formula = "~ 0 + group", qvalue = 0.05, lfc = 0) {
    spec <- readDIAStateMetadata(manager, source)$data
    c(list(formula = formula, qvalue = qvalue, lfc = lfc, weights = NA),
        list(contrast = contrasts$full_format[[i]], annotations = spec$components$protein_id_table$entity_id,
            name = contrasts$contrasts[[i]], comparison = contrasts$friendly_names[[i]]))
}

#' The compute half's prototype (D3): the flags set and restored, the declaration on `source`.
computePrototype <- function(manager, source, parameters, route) {
    context <- manager$.context
    previous <- list(context$da_reference_execution, context$da_source_reference)
    context$da_reference_execution <- TRUE
    context$da_source_reference <- source
    on.exit({ context$da_reference_execution <- previous[[1L]]; context$da_source_reference <- previous[[2L]] }, add = TRUE)
    before <- context$da_base_fit_reference
    run <- runSavedOperation(manager, diaContrastOperation(route), source, parameters)
    after <- context$da_base_fit_reference
    list(run = run, base_fit = if (!identical(before, after)) after)
}

#' A worker's owner on the controller's project: a new store handle and a bare reopen, as `computeWorkflowOmicsAction`.
workerOwner <- function(manager) {
    store <- manager$.context$store
    worker_store <- newWorkflowEntityStore(dirname(store$directory), store$identity$project_id, store$identity$workflow_id)
    owner <- newWorkflowStateOwner(worker_store, reopen = TRUE)
    owner$.context$deferred_roots <- list()
    owner$.context$deferred_owned <- list()
    owner
}

closeWorker <- function(owner) {
    store <- owner$.context$store
    closeWorkflowStateOwner(owner)
    closeWorkflowEntityStore(store)
}

#' Descriptors on disk without head snapshots and their ids: kind, layer, identity (ids masked) and factors.
entitySummary <- function(manager) {
    records <- savedDaW1Descriptors(manager)
    lapply(records, function(record) list(kind = record$kind, layer = record$metadata$factors$layer,
        identity = savedDaW1Mask(record$identity, "<no-root>"), factors = record$metadata$factors))
}

fitSummary <- function(manager) {
    Filter(function(record) identical(record$kind, "da_fit") || identical(record$layer, "full_results"), entitySummary(manager))
}

if (identical(probe, "P1")) {
    # R2: the delivered value (a contrast product hydrated by the controller, as `readManagedOmicsResult` hydrates it)
    # against the in-process `run$value`, on a computed run and on a hit; on DIA and on DIA_LIMPA.
    for (route in c("DIA", "DIA_LIMPA")) {
        fixture <- if (route == "DIA") savedDaW1Owner("k112-p1") else savedDaW1LimpaOwner()
        manager <- fixture$workflow$state_manager
        contrasts <- savedDaW1Contrasts(2L)
        source <- getDIAStateReference(manager)
        withDIAStateScope(manager, function() {
            for (i in 1:2) {
                parameters <- contrastParameters(manager, source, contrasts, i)
                computed <- quiet(computePrototype(manager, source, parameters, diaActiveRoute(manager)))
                hydrated <- hydrateSavedWorkflowEntity(manager$.context$store, computed$run$reference)
                hit <- quiet(computePrototype(manager, source, parameters, diaActiveRoute(manager)))
                say("P1 %s contrast %d: status %s then %s; identical(hydrated$value, run$value) = %s; identical(hit$value, computed$value) = %s; identical(hydrated envelope, run$value) = %s",
                    route, i, computed$run$status, hit$run$status, identical(hydrated$value, computed$run$value),
                    identical(hit$run$value, computed$run$value), identical(hydrated, computed$run$value))
                if (!identical(hydrated$value, computed$run$value)) {
                    print(all.equal(hydrated$value, computed$run$value))
                }
            }
        })
        # The whole in-process runner, then each contrast product reread from a fresh store handle (the controller's).
        result <- quiet(runDIAReferenceDaAnalysis(list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0),
            fixture$da, fixture$workflow, fixture$paths, contrasts, NULL, publish = FALSE))
        products <- manager$.context$snapshot$products
        for (name in names(result$results)) {
            value <- result$results[[name]]
            match <- Filter(function(reference) {
                envelope <- tryCatch(hydrateSavedWorkflowEntity(manager$.context$store, reference), error = function(e) NULL)
                identical(envelope$product_type, "dia-contrast-references") && identical(envelope$value$contrast, name)
            }, products)
            say("P1 %s manifest result %s: %d head product(s) for it; identical(hydrated$value, manifest result) = %s", route, name,
                length(match), length(match) == 1L && identical(hydrateSavedWorkflowEntity(manager$.context$store, match[[1L]])$value, value))
        }
    }
}

if (identical(probe, "P2")) {
    # R3 and D3's route: a worker's reopened owner, the compute prototype in its scope, against the in-process runner on a
    # twin owner. Entities by kind/layer/identity (ids masked) and factors; and the reopened owner's workflow type.
    for (route in c("DIA", "DIA_LIMPA")) {
        twin <- function() if (route == "DIA") savedDaW1Owner(paste0("k112-p2-", tolower(route)), envir = globalenv()) else
            savedDaW1LimpaOwner(envir = globalenv())
        inprocess <- twin()
        worker_fixture <- twin()
        contrasts <- savedDaW1Contrasts(2L)
        settings <- list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0)
        before_in <- length(savedDaW1Descriptors(inprocess$workflow$state_manager))
        manifest <- quiet(runDIAReferenceDaAnalysis(settings, inprocess$da, inprocess$workflow, inprocess$paths, contrasts, NULL,
            publish = FALSE))
        manager <- worker_fixture$workflow$state_manager
        before_worker <- length(savedDaW1Descriptors(manager))
        source <- getDIAStateReference(manager)
        worker <- workerOwner(manager)
        say("P2 %s: session workflow_type %s, reopened %s; route session %s, reopened %s; reopened reference_states %s",
            route, manager$workflow_type, worker$workflow_type %||% "NULL", diaActiveRoute(manager), diaActiveRoute(worker),
            usesDIAReferenceStates(worker))
        delivered <- list()
        for (i in seq_len(nrow(contrasts))) {
            parameters <- contrastParameters(manager, source, contrasts, i)
            out <- quiet(withWorkflowStateScope(worker, function() {
                inner <- getDIAActiveManager()
                stopifnot(identical(inner$.context$store$directory, worker$.context$store$directory))
                computePrototype(worker, source, parameters, diaActiveRoute(manager))
            }))
            delivered[[i]] <- out
            say("P2 %s contrast %d on the worker: status %s; base fit delivered %s; deferred roots %d", route, i, out$run$status,
                !is.null(out$base_fit), length(worker$.context$deferred_roots))
        }
        closeWorker(worker)
        a <- fitSummary(inprocess$workflow$state_manager)
        b <- fitSummary(manager)
        key <- function(x) vapply(x, function(e) paste(e$kind, e$layer %||% "", digest::digest(e$identity), digest::digest(e$factors)), "")
        say("P2 %s fits: in process %d [%s]; worker %d [%s]; same identities and factors (ids masked): %s", route, length(a),
            paste(sort(vapply(a, function(e) paste0(e$kind, ":", e$layer %||% ""), "")), collapse = ", "), length(b),
            paste(sort(vapply(b, function(e) paste0(e$kind, ":", e$layer %||% ""), "")), collapse = ", "),
            identical(sort(key(a)), sort(key(b))))
        for (e in a) if (identical(e$kind, "da_fit")) { say("P2 %s in-process base-fit factors:", route); str(e$factors) }
        for (e in b) if (identical(e$kind, "da_fit")) { say("P2 %s worker base-fit factors:", route); str(e$factors) }
        added_in <- length(savedDaW1Descriptors(inprocess$workflow$state_manager)) - before_in
        added_worker <- length(savedDaW1Descriptors(manager)) - before_worker
        kinds <- function(m) table(vapply(savedDaW1Descriptors(m), function(r) paste0(r$kind, ":", r$metadata$factors$layer %||%
            (if (identical(r$kind, "checkpoint")) r$identity$parameters else "")), ""))
        say("P2 %s entities added (no head snapshots): in process %d (incl. the arguments' commit and the manifest), worker %d", route,
            added_in, added_worker)
        print(kinds(inprocess$workflow$state_manager))
        print(kinds(manager))
        # The base-fit factors the controller would adopt (D3 R1): the delivered reference's against the in-process one's.
        delivered_fit <- Filter(Negate(is.null), lapply(delivered, `[[`, "base_fit"))
        in_fit <- inprocess$workflow$state_manager$.context$da_base_fit_reference
        say("P2 %s base-fit reference: in process set %s; delivered %d; factors identical %s; manifest model unset fields: %s", route,
            !is.null(in_fit), length(delivered_fit),
            if (length(delivered_fit)) identical(delivered_fit[[1L]]$metadata$factors, in_fit$metadata$factors) else NA,
            paste(names(manifest$model)[vapply(manifest$model, is.null, NA)], collapse = ","))
        # The contrast values: the worker's (on its own store) against the in-process manifest's, store-local ids masked.
        memo_a <- new.env(); memo_b <- new.env()
        for (i in seq_len(nrow(contrasts))) {
            name <- contrasts$contrasts[[i]]
            ca <- savedDaW1Canonical(manifest$results[[name]], inprocess$workflow$state_manager, inprocess$root, memo_a)
            cb <- savedDaW1Canonical(delivered[[i]]$run$value, manager, worker_fixture$root, memo_b)
            say("P2 %s contrast %s value canonical identical (source excluded: the in-process source is the arguments' commit): %s",
                route, name, identical(ca[setdiff(names(ca), "source")], cb[setdiff(names(cb), "source")]))
        }
    }
}

if (identical(probe, "P4")) {
    # F112-01: under a mocked managed session at the start commit, the run plugin shows neither the BH-fallback warning nor
    # the completion notice. The plugin's server on a DIA owner whose DA falls back to BH (the small recipe: q-value
    # estimation fails), once unmanaged and once with `workflowSessionApplication` mocked to an application and the
    # publication step's `apply` recorded.
    notices <- function(managed) {
        fixture <- savedDaW1Owner(paste0("k112-p4-", managed))
        manager <- fixture$workflow$state_manager
        contrasts <- savedDaW1Contrasts(1L)
        fixture$workflow$contrasts_tbl <- contrasts
        fixture$da$contrasts_available <- contrasts$contrasts
        shown <- character()
        applied <- 0L
        testthat::with_mocked_bindings(testthat::with_mocked_bindings({
            shiny::testServer(function(input, output, session) {
                session$userData$dia_capacity_manager <- manager
                da_server_run_analysis_handler(input, output, session, session$ns, fixture$da, fixture$workflow, fixture$paths)
            }, {
                session$setInputs(formula_string = "~ 0 + group", da_q_val_thresh = 0.05, treat_lfc_cutoff = 0)
                quiet(session$setInputs(run_da_analysis = 1L))
                session$flushReact()
            })
        }, workflowSessionApplication = function(session) if (managed) list(application = "mocked") else NULL,
        installWorkflowManagedStep = function(...) list(apply = function(...) applied <<- applied + 1L,
            active = function() FALSE, pending = function() FALSE),
        .package = "MultiScholaR"),
        showNotification = function(ui, action = NULL, duration = 5, closeButton = TRUE, id = NULL,
            type = c("default", "message", "warning", "error"), session = NULL) {
            shown <<- c(shown, paste0("show [", type[[1L]], "] ", if (is.null(id)) "" else paste0("(", id, ") "),
                paste(as.character(ui), collapse = "")))
            invisible(id)
        }, removeNotification = function(id, session = NULL) {
            shown <<- c(shown, paste0("remove (", id, ")"))
            invisible(id)
        }, .package = "shiny")
        snapshot <- withDIAStateScope(manager, function() manager$.context$snapshot)
        manifest <- if (isSavedWorkflowReference(snapshot$da_results)) resolveWorkflowStateValue(manager$.context, snapshot$da_results)
        warnings <- if (!is.null(manifest)) unlist(lapply(manifest$results, function(result) result$qvalue_warnings), recursive = FALSE)
        list(shown = shown, applied = applied, manifest = !is.null(manifest), warnings = length(warnings))
    }
    for (managed in c(FALSE, TRUE)) {
        out <- notices(managed)
        say("P4 managed = %s: publication step applied %d time(s); manifest on the head %s, with %d q-value warning(s) (BH fallback); notices:",
            managed, out$applied, out$manifest, out$warnings)
        for (line in out$shown) say("    %s", line)
    }
}
