# KIT-107 W2: the enrichment core (D3, D4), written before it.
#   - each moved function is in the package namespace and its source is its new file (utils_workflow_enrichment.R,
#     utils_workflow_dia_enrichment.R);
#   - the equivalence: `runWorkflowEnrichment` on W1 (b)'s fixture (test-omic-enrich-module.R, its helpers loaded below)
#     gives the start commit's run step for the two backends and the two kernel cases: the kernel's arguments (and the
#     fake services' query and background), the product, the state digests, the history, `tab_status`, the external
#     resources, the bus values, the project's files and store counts, the display tables and the progress steps.
#     The core side compares with W1 (b)'s frozen run step and stays live. The comparison side runs the start commit's
#     `runProtEnrichAnalysisBody` (mod_prot_enrich_analysis_helpers.R) on the same fixture and compares it with the core
#     side: one block, quarantined when W3 deletes that helper (Q-166);
#   - `restoreWorkflowEnrichment` gives W1 (b)'s `restore` capture;
#   - `refuse_missing_da` refuses a run whose DA table is missing before anything is written (F107-03's guard; W3 wires
#     it into the tab);
#   - F-16 (D6): test 4, a failed commit leaves no pending provenance; test 1, the module-level record after a Run (it
#     fails until W3 puts the tab on the core: the start commit's module still commits a `legacy` step).

# W1 (b)'s helpers and frozen captures: every top-level assignment of test-omic-enrich-module.R (its blocks are not run).
for (omicsEnrichCoreW1Expression in parse(testthat::test_path("test-omic-enrich-module.R"), keep.source = FALSE)) {
    if (is.call(omicsEnrichCoreW1Expression) && identical(omicsEnrichCoreW1Expression[[1L]], as.name("<-"))) {
        eval(omicsEnrichCoreW1Expression, environment())
    }
}
rm(omicsEnrichCoreW1Expression)

#' The proteomics `enrichment` group with its function fields (D1). Until the W2 merge wires them in
#' utils_omic_gui_config.R they are optional; this fills an absent one with the value the merge wires, and changes
#' nothing once they are wired.
omicsEnrichCoreConfig <- function(config = omicGuiConfig("proteomics")$enrichment) {
    config[["organism"]] <- config[["organism"]] %||% list(supported = buildProtEnrichSupportedOrganisms,
        method = resolveProtEnrichAnalysisMethod)
    config$annotation <- config$annotation %||% list(resolve = resolveProtEnrichUniprotAnnotations,
        match = resolveProtEnrichAnnotationMatching)
    config$kernel <- config$kernel %||% list(select = resolveProtEnrichProcessEnrichmentsFn,
        test_double = runProtEnrichDeterministicProcessEnrichments)
    config$organism_filter$mapping <- config$organism_filter$mapping %||% resolveProtEnrichOrganismMapping
    config$organism_filter$apply <- config$organism_filter$apply %||% applyProtEnrichOrganismFilter
    config
}

test_that("(W2) the moved functions are in the namespace, each in its new file", {
    core <- c("resolveProtEnrichAnalysisInputColumns", "buildProtEnrichProcessEnrichmentsArgs",
        "buildProtEnrichAllContrastResults", "resolveProtEnrichProcessEnrichmentsFn", "resolveProtEnrichAnalysisMethod",
        "persistProtEnrichOrganismFilterMetadata", "resolveProtEnrichRawContrastName", "resolveProtEnrichContrastsTable",
        "resolveProtEnrichSelectedContrastResults", "resolveProtEnrichSelectedDaResults", "resolveProtEnrichRunDependencies",
        "resolveProtEnrichOutputDirectories", "buildProtEnrichSupportedOrganisms", "normaliseProtEnrichTaxonId",
        "validateProtEnrichProcessParameters", "buildProtEnrichAnalysisResultsPayload", "propagateProtEnrichResultsArgs",
        "propagateProtEnrichUiParams", "updateProtEnrichStateManagerUiParams",
        "runWorkflowEnrichment", "restoreWorkflowEnrichment", "saveWorkflowEnrichmentState", "replayWorkflowEnrichmentStep")
    dia <- c("buildProtEnrichDeterministicUniprotAnnotations", "buildProtEnrichDeterministicGprofilerResult",
        "buildProtEnrichDeterministicClusterProfilerResult", "writeProtEnrichDeterministicResultTable",
        "buildProtEnrichDeterministicPlot", "runProtEnrichDeterministicProcessEnrichments", "resolveProtEnrichOrganismMapping",
        "applyProtEnrichOrganismFilter", "resolveProtEnrichUniprotAnnotations", "resolveProtEnrichAnnotationMatching")
    moved <- c(stats::setNames(rep("utils_workflow_enrichment.R", length(core)), core),
        stats::setNames(rep("utils_workflow_dia_enrichment.R", length(dia)), dia))
    expect_length(moved, 33L)
    ns <- asNamespace("MultiScholaR")
    for (name in names(moved)) {
        expect_true(exists(name, envir = ns, inherits = FALSE), info = name)
        if (!exists(name, envir = ns, inherits = FALSE)) next
        fn <- get(name, envir = ns, inherits = FALSE)
        expect_identical(environmentName(environment(fn)), "MultiScholaR", info = name)
        # Source files are known under load_all only; the installed package keeps no srcref (KIT-104 ca35c462's rule).
        if (!is.null(utils::getSrcref(fn))) expect_identical(utils::getSrcFilename(fn), moved[[name]], info = name)
    }
    expect_true(methods::isClass("ProtEnrichDeterministicResult"))
    if (!is.null(utils::getSrcref(get("runWorkflowEnrichment", envir = ns)))) {
        source <- readLines(testthat::test_path("..", "..", "R", "utils_workflow_dia_enrichment.R"))
        expect_true(any(grepl('methods::isClass("ProtEnrichDeterministicResult")', source, fixed = TRUE)))
        # The core calls no Shiny UI function (D3): the entries notify and report progress through callbacks.
        core_source <- readLines(testthat::test_path("..", "..", "R", "utils_workflow_enrichment.R"))
        expect_false(any(grepl("shiny::(showNotification|incProgress|withProgress|removeNotification|update[A-Z])",
            core_source)))
    }
})

# ------------------------------------------------------------------------------------------------
# The equivalence (D3): W1 (b)'s fixture, the DA through the module, then the run.
# ------------------------------------------------------------------------------------------------

#' The run's inputs as the tab holds them after W1 (b)'s select step (the UI's defaults for the rest).
omicsEnrichCoreInputs <- function(taxon) list(selected_contrast = "B_vs_A", organism_taxid = taxon, up_cutoff = 0,
    down_cutoff = 0, q_cutoff = 0.05, correction_method = "gSCS", enable_organism_filter = FALSE)

#' The four cases: W1 (b)'s case, its taxon and whether the kernel is the live one.
omicsEnrichCoreCases <- list(gprofiler = list(taxon = "9606", live = FALSE),
    clusterprofiler = list(taxon = "999999", live = FALSE), kernel_gprofiler = list(taxon = "9606", live = TRUE),
    kernel_clusterprofiler = list(taxon = "999999", live = TRUE))

#' The fields of a run step that the run itself produces (the tab's outputs, notices, selector updates and returned
#' values are the frame's, W3).
omicsEnrichCoreFields <- c("kernel", "services", "states", "history", "tab_status", "resources", "bus", "files",
    "errors")

omicsEnrichCoreRuns <- new.env(parent = emptyenv())

#' One session of W1 (b)'s case up to its select step, then `run(ctx, inputs)` as the run step.
omicsEnrichCoreDrive <- function(case, run, env = parent.frame()) {
    spec <- omicsEnrichCoreCases[[case]]
    w1 <- omicEnrichW1CaseSpec(case, env)[[1L]]
    inputs <- omicsEnrichCoreInputs(spec$taxon)
    # The steps before W1 (b)'s run step (a backend case starts with `start`), so each step's "<unchanged>" fields
    # are relative to the same previous step.
    before <- intersect(c("start", "da", "open", "select"), names(omicEnrichW1Frozen[[case]]))
    steps <- c(w1$steps[before], list(run = function(ctx) run(ctx, inputs)))
    omicEnrichW1Drive(w1$setup, steps, live = spec$live)
}

#' The core side: `runWorkflowEnrichment` on the bus and the manifest the DA published, with the frame's callbacks;
#' the product and display tables kept from its returned list.
omicsEnrichCoreRun <- function(case, env = parent.frame()) {
    if (is.null(omicsEnrichCoreRuns[[case]])) {
        config <- omicsEnrichCoreConfig()
        omicsEnrichCoreRuns[[case]] <- omicsEnrichCoreDrive(case, function(ctx, inputs) shiny::isolate({
            bus <- ctx$bus()
            manifest <- bus$da_analysis_results_list
            out <- runWorkflowEnrichment(bus, config, manifest, inputs, ctx$paths,
                notify = function(message, type, duration) shiny::showNotification(message, type = type, duration = duration),
                progress = function(value, detail) shiny::incProgress(value, detail = detail))
            ctx$keep("run_product", omicEnrichW1Product(out$reference, bus$state_manager, manifest))
            ctx$keep("run_display", omicEnrichW1Display(out$contrasts, bus$state_manager))
            ctx$keep("result", c(names = paste(names(out), collapse = ","), method = out$method$method,
                raw_contrast = out$raw_contrast, selected = out$selected$found, saved = out$saved$saved,
                saved_warning = out$saved$warning %||% "",
                manifest_unchanged = identical(bus$da_analysis_results_list, manifest)))
        }), env = env)
    }
    omicsEnrichCoreRuns[[case]]
}

for (case in names(omicsEnrichCoreCases)) {
    local({
        case <- case
        test_that(sprintf("(W2, core side) runWorkflowEnrichment gives W1 (b)'s %s run step", case), {
            seen <- omicsEnrichCoreRun(case)
            frozen <- omicEnrichW1Frozen[[case]]
            expect_identical(names(seen), c(intersect(c("start", "da", "open", "select"), names(frozen)), "run"))
            for (step in setdiff(names(seen), "run")) {
                expect_identical(seen[[step]][omicsEnrichCoreFields], frozen[[step]][omicsEnrichCoreFields],
                    info = paste(case, step))
            }
            run <- seen$run
            for (field in omicsEnrichCoreFields) {
                expect_identical(run[[field]], frozen$run[[field]], info = paste(case, "run", field))
            }
            expect_identical(run$run_product, frozen$run$product, info = case)
            expect_identical(run$run_display, frozen$run$display, info = case)
            # The progress steps (the frame's withProgress start line is the frame's).
            expect_identical(run$progress, frozen$run$progress[-1L], info = case)
            expect_identical(run$notes, character(), info = case)
            expect_identical(run$result[["names"]], paste(c("results", "method", "contrasts", "selected", "raw_contrast",
                "contrasts_tbl", "pathway_dir", "annotation_match", "organism_filter", "service", "reference", "parameters",
                "saved"), collapse = ","))
            expect_identical(run$result[["method"]], if (identical(omicsEnrichCoreCases[[case]]$taxon, "9606"))
                "gprofiler2" else "clusterprofiler")
            expect_identical(run$result[["raw_contrast"]], "groupB-groupA")
            expect_identical(run$result[["selected"]], "TRUE")
            # The state commits when it did at the start commit (W1 (b): the live kernel cases' commit fails today and is
            # reported as a warning; F-16 then leaves no provenance pending).
            expect_identical(run$result[["saved"]], if (identical(frozen$run$history, "<unchanged>")) "FALSE" else "TRUE",
                info = paste(case, run$result[["saved_warning"]]))
            # F107-01: the run never replaces the manifest the caller holds.
            expect_identical(run$result[["manifest_unchanged"]], "TRUE")
        })
    })
}

# Comparison side: quarantined when W3 deletes mod_prot_enrich_analysis_helpers.R (`runProtEnrichAnalysisBody`), Q-166.
test_that("(W2, comparison side) runWorkflowEnrichment gives the start commit's runProtEnrichAnalysisBody on the same fixture", {
    testthat::skip("quarantined: Q-166 (KIT-107 W3) the comparison side runs the start commit's `runProtEnrichAnalysisBody` (mod_prot_enrich_analysis_helpers.R, deleted by W3); the core side's blocks above stay live and compare `runWorkflowEnrichment` with W1 (b)'s frozen run steps, and W1 (b) drives the tab on the core")
    for (case in names(omicsEnrichCoreCases)) {
        core <- omicsEnrichCoreRun(case)
        body <- omicsEnrichCoreDrive(case, function(ctx, inputs) shiny::isolate({
            bus <- ctx$bus()
            tab <- ctx$returned()
            runProtEnrichAnalysisBody(input = inputs, enrichmentData = tab, workflowData = bus, experimentPaths = ctx$paths,
                currentAnalysisMethodFn = function() resolveProtEnrichAnalysisMethod(inputs$organism_taxid,
                    buildProtEnrichSupportedOrganisms())$methodInfo)
            ctx$keep("run_product", omicEnrichW1Product(tab$enrichment_product, bus$state_manager, bus$da_analysis_results_list))
            ctx$keep("run_display", omicEnrichW1Display(tab$all_enrichment_results, bus$state_manager))
        }))
        for (field in c(omicsEnrichCoreFields, "run_product", "run_display", "progress", "notes")) {
            expect_identical(core$run[[field]], body$run[[field]], info = paste(case, field))
        }
    }
})

test_that("(W2) restoreWorkflowEnrichment gives W1 (b)'s restore capture", {
    config <- omicsEnrichCoreConfig()
    w1 <- omicEnrichW1CaseSpec("restore", environment())[[1L]]
    steps <- w1$steps
    steps$restore <- function(ctx) {
        manager <- shiny::isolate(ctx$bus()$state_manager)
        out <- restoreWorkflowEnrichment(manager, config)
        ctx$keep("restore", c(restored = !is.null(out),
            enrichment_results_full = omicEnrichW1Digest(out$enrichment_results_full, ctx$roots),
            omicEnrichW1Display(out$all_enrichment_results, manager)))
        ctx$keep("method", out$method$method)
        ctx$keep("absent", is.null(restoreWorkflowEnrichment(manager, utils::modifyList(config,
            list(state_name = "not_a_state")))))
    }
    seen <- omicEnrichW1Drive(w1$setup, steps)
    frozen <- omicEnrichW1Frozen$restore$restore$restore
    expect_identical(seen$restore$restore, frozen[setdiff(names(frozen), c("again", "analysis_complete"))])
    expect_identical(seen$restore$method, "gprofiler2")
    expect_true(seen$restore$absent)
    expect_null(restoreWorkflowEnrichment(NULL, config))
})

test_that("(W2) refuse_missing_da refuses a run whose DA table is missing, before anything is written (F107-03's guard)", {
    config <- omicsEnrichCoreConfig()
    setup <- omicEnrichW1Setup(env = environment())
    steps <- list(da = omicEnrichW1Da, open = omicEnrichW1Open, refuse = function(ctx) shiny::isolate({
        unlink(file.path(ctx$paths$da_output_dir, "da_proteins_groupB-groupA_long_annot.tsv"))
        bus <- ctx$bus()
        manager <- bus$state_manager
        before <- list(history = diaStateHistory(manager), products = length(manager$.context$snapshot$products))
        error <- tryCatch(runWorkflowEnrichment(bus, config, bus$da_analysis_results_list, omicsEnrichCoreInputs("9606"),
            ctx$paths, refuse_missing_da = TRUE), error = function(e) e)
        ctx$keep("refused", c(class = paste(intersect(class(error), "workflow_enrichment_input"), collapse = ""),
            message = conditionMessage(error),
            history = identical(diaStateHistory(manager), before$history),
            products = identical(length(manager$.context$snapshot$products), before$products),
            filter = is.null(bus$enrichment_organism_filter), annotations = is.null(bus$uniprot_dat_cln),
            results = is.null(bus$enrichment_analysis_results)))
    }))
    seen <- suppressWarnings(omicEnrichW1Drive(setup, steps))
    expect_identical(seen$refuse$refused, c(class = "workflow_enrichment_input",
        message = "The DA results table is missing for contrast groupB-groupA.", history = "TRUE", products = "TRUE",
        filter = "TRUE", annotations = "TRUE", results = "TRUE"))
    expect_identical(seen$refuse$history, "<unchanged>")
})

# ------------------------------------------------------------------------------------------------
# F-16 (D4, D6).
# ------------------------------------------------------------------------------------------------

test_that("(F-16, D6 test 4) the enrichment state records its provenance, and a failed commit leaves none pending", {
    root <- withr::local_tempdir()
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "kit107-provenance", "DIA"), WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$reference_states <- TRUE
    manager$saveState("correlation_filtered", makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small")),
        list(), "fixture")
    config <- omicsEnrichCoreConfig()
    results <- runProtEnrichDeterministicProcessEnrichments(da_results = module_ci_prot_enrich_da_results(), taxon_id = 9606,
        up_cutoff = 0, down_cutoff = 0, q_cutoff = 0.05, pathway_dir = withr::local_tempdir(), go_annotations = NULL,
        exclude_iea = FALSE, protein_id_column = "uniprot_acc", contrast_names = "groupB-groupA", correction_method = "gSCS")
    declared <- list(taxon_id = 9606, up_cutoff = 0, down_cutoff = 0, q_cutoff = 0.05, exclude_iea = FALSE,
        correction_method = "gSCS", protein_id_column = "uniprot_acc", contrast_names = "groupB-groupA",
        results = c(`groupB-groupA` = strrep("a", 64L)))
    method <- list(method = "gprofiler2", supported = TRUE)
    parameters <- list(selected_contrast = "B_vs_A", organism_taxid = "9606", up_cutoff = 0, down_cutoff = 0, q_cutoff = 0.05)
    workflow <- list2env(list(state_manager = manager))
    # A failed commit: the provenance set before it is cleared, the next commit is not attributed to the enrichment.
    save <- function() {
        out <- NULL
        utils::capture.output(out <- saveWorkflowEnrichmentState(workflow, config, results, declared, method, parameters,
            "pathway"))
        out
    }
    failed <- (function() {
        testthat::local_mocked_bindings(saveDIAOwnerState = function(manager, ...) stop("disk full"))
        save()
    })()
    expect_identical(failed, list(attempted = TRUE, saved = FALSE, warning = "disk full"))
    expect_null(manager$.context$pending_provenance)
    manager$saveState("next_state", makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small")), list(), "next")
    expect_identical(manager$.context$snapshot$step_provenance$next_state$operation, "legacy")
    # The commit: the registry entry's operation, variant and names, and the declared values, on the head.
    saved <- save()
    expect_identical(saved, list(attempted = TRUE, saved = TRUE, warning = NULL))
    expect_null(manager$.context$pending_provenance)
    expect_identical(manager$current_state, "enrichment_completed")
    entry <- manager$.context$snapshot$step_provenance$enrichment_completed
    expect_identical(entry$operation, "enrichment")
    expect_identical(entry$variant, "whole")
    expect_identical(entry$parameter_names, diaEnrichmentOperation()$parameters)
    expect_identical(entry$values, declared)
    expect_null(entry$fixed)
    # The state's config and description as today's save writes them (a result-class state: KIT-101 F-15's read).
    metadata <- readComponentMetadata(manager, manager$.context$snapshot$states$enrichment_completed)
    expect_identical(metadata$config, list(selected_contrast = "B_vs_A", analysis_method = "gprofiler2",
        organism_supported = TRUE, up_cutoff = 0, down_cutoff = 0, q_cutoff = 0.05, organism_taxid = "9606",
        pathway_dir = "pathway"))
    expect_identical(metadata$description, "Enrichment analysis completed using gprofiler2 for contrast: B_vs_A")
})

test_that("(F-16, D6 test 1) after a Run the record's enrichment step carries its operation, variant and 9 declared parameters (passes from W3)", {
    w1 <- omicEnrichW1CaseSpec("gprofiler", environment())[[1L]]
    steps <- c(w1$steps[c("da", "open", "select", "run")], list(record = function(ctx) {
        manager <- shiny::isolate(ctx$bus()$state_manager)
        record <- jsonlite::fromJSON(jsonlite::toJSON(buildWorkflowAnalysisRecord(manager), auto_unbox = TRUE,
            digits = I(17), null = "null", na = "null"), simplifyVector = FALSE)
        ctx$keep("record", list(step = record$steps[[length(record$steps)]],
            incomplete = validateWorkflowAnalysisRecord(record)$incomplete))
    }))
    seen <- omicEnrichW1Drive(w1$setup, steps)
    step <- seen$record$record$step
    expect_identical(step$state, "enrichment_completed")
    expect_identical(step$operation, "enrichment")
    expect_identical(step$variant, "whole")
    parameters <- step$parameters
    expect_identical(vapply(parameters, `[[`, "", "name"), diaEnrichmentOperation()$parameters)
    expect_identical(vapply(parameters, `[[`, "", "source"), rep("provenance", 9L))
    expect_identical(vapply(parameters, `[[`, "", "provenance"), c("input", "no_shipped_default", "no_shipped_default",
        "no_shipped_default", "no_shipped_default", "no_shipped_default", "input", "input", "input"))
    expect_false(any(grepl("enrichment_completed", seen$record$record$incomplete, fixed = TRUE)))
})

# ------------------------------------------------------------------------------------------------
# KIT-113 W2: the enrichment on the managed worker, at core level (D3, D6). The frame is today's (W3 changes it); its
# call of the core is replaced by the managed halves as the runner runs them: `requestWorkflowEnrichment` in the
# controller, one `enrichment` action through the worker's entry point (`workflowWorkerMain`) on the request's named
# input (the DA manifest's source), the controller's acceptance and omics completion, then `acceptWorkflowEnrichment`.
# Each of KIT-107's eight W1 (b) cases then runs its frozen steps and must give its frozen in-process values and the
# KIT-113 W1 block's, apart from D6's recorded exclusions, each asserted as it is.
# ------------------------------------------------------------------------------------------------

#' One managed action, in this process (KIT-111 W2's `protNormW2Dispatch`, test-prot-norm-core.R, reproduced here:
#' test files do not see each other's locals): the request with its named input through the worker's entry point, the
#' acceptance and what the tab reads. A worker failure is raised with the worker's message; a refused acceptance
#' returns its status and the worker's value only.
omicsEnrichW2Dispatch <- function(manager, controller, paths, parameters, input, context) {
    store <- manager$.context$store
    operation <- "enrichment"
    id <- paste0("wfx-", substr(digest::digest(list(tempfile(), operation), algo = "sha256"), 1L, 32L))
    export <- workflowExportPaths(store$directory, id)
    token <- list(request_id = id, incarnation = store$identity$incarnation, expected_head = workflowExecutionHead(store)$entity_id,
        generation = controller$control$generation, controller_id = controller$control$controller_id)
    writeWorkflowExecutionMessage(file.path(export$request_dir, "publication-context.rds"),
        list(experiment_label = "core", experiment_paths = paths, omic_type = "proteomics"))
    writeWorkflowExportJournal(export, list(schema = "multischolar.workflow_action.v1", request_id = id, status = "admitting",
        expected_head = token$expected_head, token_hash = digest::digest(token, algo = "sha256")))
    request <- list(schema = 1L, build = workflowExecutionBuildIdentity(), operation = operation,
        input_references = list(state = input), parameters = parameters, token = token, owner = store$identity,
        project_root = dirname(store$directory), directory = export$request_dir, config = list(mode = "adapters"),
        worker_cap = 12 * 1024^3)
    path <- file.path(export$request_dir, "request.rds")
    writeWorkflowExecutionMessage(path, request)
    withr::with_envvar(c(MULTISCHOLAR_WORKFLOW_WORKER = "1"), suppressWarnings(suppressMessages(utils::capture.output(try(
        workflowWorkerMain(c(path, digest::digest(file = path, algo = "sha256"))), silent = TRUE)))))
    message <- readRDS(file.path(export$request_dir, "worker-result.rds"))
    if (!identical(message$status, "complete")) stop(message$error, call. = FALSE)
    controller$job <- list(request = list(token = token, operation = operation, config = list(mode = "adapters"),
        parameters = parameters), paths = export, completion_context = validateWorkflowOmicsCompletionContext(operation, context))
    receipt <- acceptWorkflowActionResult(controller, message[c("schema", "status", "token", "result")])
    if (!identical(receipt$status, "complete")) return(list(acceptance = receipt$status, result = message$result))
    accepted <- resumeWorkflowAction(controller, id)
    c(readManagedOmicsResult(manager, accepted, NULL), list(acceptance = receipt$status, result = message$result))
}

#' The build identity fixed and the service teardown replaced (KIT-111 W2's `protNormW2Mocks`), for the life of `env`.
omicsEnrichW2Mocks <- function(env = parent.frame()) {
    testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
        finishWorkflowExecutionJob = function(owner, status) { owner$job <- NULL; list(status = status) },
        .package = "MultiScholaR", .env = env)
}

#' The managed path in place of the frame's core call (`runWorkflowEnrichment`'s signature): the request half, one
#' dispatch on a controller kept per owner (the in-process test worker), then the accept half. `log` keeps each
#' dispatch: its input and the head's state when it was made, the acceptance, the provenance the completion left
#' pending, and the sizes the registry and the acceptance measure (`object.size`).
omicsEnrichW2Managed <- function(log) {
    function(workflow, config, manifest, parameters, paths, process_fn = NULL, notify = NULL, progress = NULL,
        refuse_missing_da = FALSE) {
        manager <- workflow$state_manager
        half <- requestWorkflowEnrichment(workflow, config, manifest, parameters, paths, notify, progress)
        key <- manager$.context$store$directory
        if (is.null(log$controllers[[key]])) log$controllers[[key]] <- omicsWorkerController(list(owner = manager))
        n <- length(log$dispatches) + 1L
        log$dispatches[[n]] <- list(input = manifest$source$entity_id, head = getWorkflowStateReference(manager)$entity_id,
            head_state = manager$current_state, parameters = as.numeric(object.size(half$parameters)), acceptance = "failed")
        delivered <- omicsEnrichW2Dispatch(manager, log$controllers[[key]], paths, half$parameters, manifest$source,
            list(state_name = config$state_name, config = list(), description = "enrichment"))
        log$dispatches[[n]]$acceptance <- delivered$acceptance
        log$dispatches[[n]]$value <- as.numeric(object.size(delivered$result))
        if (!identical(delivered$acceptance, "complete")) stop("the acceptance refused the worker's value", call. = FALSE)
        log$dispatches[[n]]$pending <- manager$.context$pending_provenance$operation %||% "<none>"
        out <- acceptWorkflowEnrichment(workflow, config, half$context, delivered, progress)
        log$dispatches[[n]]$pending_after <- manager$.context$pending_provenance$operation %||% "<none>"
        out
    }
}

omicsEnrichW2Runs <- new.env(parent = emptyenv())

#' One W1 (b) case on the managed path, run once per process: each session's frozen steps, then the KIT-113 W1 block's
#' capture (`kit113_d6`). Returns the captures by step (a session's steps prefixed as W1 (b) prefixes them) and the
#' dispatches.
omicsEnrichW2Case <- function(case) {
    if (is.null(omicsEnrichW2Runs[[case]])) {
        log <- new.env(parent = emptyenv())
        log$controllers <- list()
        log$dispatches <- list()
        omicsEnrichW2Mocks()
        testthat::local_mocked_bindings(runWorkflowEnrichment = omicsEnrichW2Managed(log), .package = "MultiScholaR")
        sessions <- omicEnrichW1CaseSpec(case, environment())
        captures <- list()
        for (session_name in names(sessions) %||% "") {
            spec <- if (nzchar(session_name)) sessions[[session_name]] else sessions[[1L]]
            seen <- omicEnrichW1Drive(spec$setup, c(spec$steps,
                list(kit113_d6 = function(ctx) ctx$keep("d6", omicEnrichW113Fields(ctx)))), live = spec$live)
            if (nzchar(session_name)) names(seen) <- paste0(session_name, "/", names(seen))
            captures <- c(captures, seen)
        }
        omicsEnrichW2Runs[[case]] <- list(captures = captures, dispatches = log$dispatches)
    }
    omicsEnrichW2Runs[[case]]
}

#' A field as a step holds it, "<unchanged>" resolved to the last step of the same session that set it.
omicsEnrichW2Held <- function(captures, step, field) {
    steps <- names(captures)
    i <- match(step, steps)
    while (i > 1L && identical(captures[[steps[[i]]]][[field]], "<unchanged>")) i <- i - 1L
    captures[[steps[[i]]]][[field]]
}

#' The store's `dia-state` entries (the saved state values, the transfer values among them) apart from the rest of the
#' project's files.
omicsEnrichW2StateEntries <- "<store entries>..workflow-entities/entities/dia-state-<ID>"

#' The transfer values one managed Run that completes adds to the store's `dia-state` entries (D6's first recorded
#' exclusion): the annotation table and the results.
omicsEnrichW2Transfers <- 2L

#' D6's recorded exclusions per case: `transfers`, the `dia-state` entries each step's managed Run adds (a Run the
#' worker refuses keeps the request's annotation table and mapping); `progress`, the steps whose progress gains the
#' request half's 0.5 before the worker refuses (D3, the progress before dispatch); `annotations`, the steps whose
#' request half writes `uniprot_dat_cln` before the worker refuses a missing DA table; `no_filter`, the steps whose
#' failure (on the worker, or the request half's taxon refusal: F113-03) leaves `enrichment_organism_filter` unwritten.
omicsEnrichW2Exclusions <- function(case) {
    s <- omicsEnrichW2Transfers
    switch(case,
        gprofiler = , clusterprofiler = , kernel_gprofiler = , kernel_clusterprofiler = , restore = list(transfers = c(run = s)),
        organism_filter = list(transfers = c(`no_mapping/run` = s, `mapped/run` = 2L), progress = "mapped/run"),
        refusals = list(transfers = c(`pre_da/run_taxon_abc` = 0L, `q_cutoff/run_q_cutoff` = 1L,
            `table_removed/run_table_removed` = 1L), progress = "table_removed/run_table_removed",
            annotations = "table_removed/run_table_removed",
            no_filter = c("pre_da/run_taxon_abc", "q_cutoff/run_q_cutoff")),
        prefill = list(transfers = integer()))
}

#' The bus capture with `enrichment_organism_filter` unwritten.
omicsEnrichW2NoFilter <- function(bus) {
    at <- which(startsWith(names(bus), "enrichment_organism_filter$"))
    c(bus[seq_len(at[[1L]] - 1L)], `enrichment_organism_filter$` = "<NULL>", bus[-seq_len(max(at))])
}

#' Compare one managed case with W1 (b)'s frozen captures: every field of every step identical, apart from D6's
#' recorded exclusions, each asserted as it is.
expectOmicsEnrichW2Case <- function(case) {
    seen <- omicsEnrichW2Case(case)$captures
    frozen <- omicEnrichW1Frozen[[case]]
    excluded <- omicsEnrichW2Exclusions(case)
    expect_identical(setdiff(names(seen), grep("kit113_d6$", names(seen), value = TRUE)), names(frozen))
    added <- 0L
    session <- ""
    for (step in names(frozen)) {
        prefix <- if (grepl("/", step, fixed = TRUE)) sub("/.*$", "", step) else ""
        if (!identical(prefix, session)) {
            session <- prefix
            added <- 0L
        }
        if (step %in% names(excluded$transfers)) added <- added + excluded$transfers[[step]]
        for (field in setdiff(union(names(frozen[[step]]), names(seen[[step]])), "files")) {
            expected <- frozen[[step]][[field]]
            if (identical(field, "progress") && step %in% excluded$progress) expected <- c(expected, "0.5 Running enrichment analysis...")
            if (identical(field, "bus") && step %in% excluded$annotations) {
                expected <- omicsEnrichW2Held(frozen, step, "bus")
                expected[["uniprot_dat_cln"]] <- "5x9 13df57c377de46dc"
            }
            if (identical(field, "bus") && step %in% excluded$no_filter) expected <- omicsEnrichW2NoFilter(expected)
            expect_identical(seen[[step]][[field]], expected, info = paste(case, step, field))
        }
        files <- omicsEnrichW2Held(seen, step, "files")
        files_frozen <- omicsEnrichW2Held(frozen, step, "files")
        others <- function(x) x[sort(setdiff(names(x), omicsEnrichW2StateEntries))]
        expect_identical(others(files), others(files_frozen), info = paste(case, step, "files"))
        expect_identical(as.integer(files[[omicsEnrichW2StateEntries]]) - as.integer(files_frozen[[omicsEnrichW2StateEntries]]),
            added, info = paste(case, step, "transfer values"))
    }
}

for (case in c("gprofiler", "clusterprofiler", "kernel_gprofiler", "kernel_clusterprofiler", "organism_filter", "refusals",
    "restore", "prefill")) {
    local({
        name <- case
        test_that(sprintf("(KIT-113 W2) managed at core level: %s gives W1 (b)'s frozen values, apart from D6's exclusions", name), {
            skip_if_not(workflowOwnerCapable())
            expectOmicsEnrichW2Case(name)
            run <- omicsEnrichW2Case(name)
            # The KIT-113 W1 block's fields after each session (the committed record, results, product, provenance and
            # the analysis record's step).
            for (key in grep("kit113_d6$", names(run$captures), value = TRUE)) {
                session <- sub("/?kit113_d6$", "", key)
                frozen <- if (nzchar(session)) paste0(name, "/", session) else name
                expect_identical(run$captures[[key]]$d6, omicEnrichW113Frozen[[frozen]], info = frozen)
            }
            # Each dispatch names the DA manifest's source; the completion leaves the request's provenance pending and
            # the accept half replaces it with the commit's; the request and the value are within their limits (R3).
            for (dispatch in run$dispatches) {
                expect_identical(dispatch$input, dispatch$head, info = name)
                expect_lt(dispatch$parameters, 65536)
                if (identical(dispatch$acceptance, "complete")) {
                    expect_lte(dispatch$value, 8388608)
                    expect_identical(c(dispatch$pending, dispatch$pending_after), c("enrichment", "<none>"), info = name)
                }
            }
            expect_identical(vapply(run$dispatches, `[[`, "", "acceptance"), switch(name,
                organism_filter = c("complete", "failed"), refusals = c("failed", "failed"), prefill = character(),
                "complete"), info = name)
        })
    })
}

#' W1 (b)'s gprofiler session to its select step, then `steps` and the KIT-113 W1 block's capture, on the managed path.
#' `log$fail` TRUE makes the worker's compute half fail. Returns the captures and the dispatches.
omicsEnrichW2Drive <- function(steps, log = new.env(parent = emptyenv())) {
    log$controllers <- list()
    log$dispatches <- list()
    omicsEnrichW2Mocks()
    compute <- get("computeWorkflowEnrichment", envir = asNamespace("MultiScholaR"))
    testthat::local_mocked_bindings(runWorkflowEnrichment = omicsEnrichW2Managed(log),
        computeWorkflowEnrichment = function(manager, config, parameters) {
            if (isTRUE(log$fail)) stop("forced worker failure", call. = FALSE)
            compute(manager, config, parameters)
        }, .package = "MultiScholaR")
    w1 <- omicEnrichW1CaseSpec("gprofiler", environment())[[1L]]
    seen <- omicEnrichW1Drive(w1$setup, c(w1$steps[c("da", "open", "select")], steps,
        list(kit113_d6 = function(ctx) ctx$keep("d6", omicEnrichW113Fields(ctx)))))
    list(captures = seen, dispatches = log$dispatches)
}

test_that("(KIT-113 W2, F107-01) a second managed Run names the DA manifest's source while the head is enrichment_completed, and gives the first Run's values", {
    skip_if_not(workflowOwnerCapable())
    run <- omicsEnrichW2Drive(list(run = omicEnrichW1Run(1L), again = omicEnrichW1Run(2L)))
    seen <- run$captures
    expect_identical(seen$again$errors, character())
    expect_identical(seen$again$notes, seen$run$notes)
    expect_identical(seen$again$progress, seen$run$progress)
    for (field in c("history", "states", "product", "display")) expect_identical(seen$again[[field]], "<unchanged>", info = field)
    expect_identical(seen$again$values, seen$run$values)
    expect_gt(as.integer(seen$again$files[[omicsEnrichW2StateEntries]]),
        as.integer(seen$run$files[[omicsEnrichW2StateEntries]]))
    first <- run$dispatches[[1L]]
    second <- run$dispatches[[2L]]
    expect_identical(c(first$head_state, second$head_state), c("correlation_filtered", "enrichment_completed"))
    expect_identical(second$input, first$input)
    expect_false(identical(second$input, second$head))
    expect_identical(c(first$acceptance, second$acceptance), c("complete", "complete"))
    expect_identical(second$pending_after, "<none>")
    fields <- c("history", "results", "provenance", "record")
    expect_identical(seen$kit113_d6$d6[fields], omicEnrichW113Frozen$gprofiler[fields])
    # The second commit's record is the first's, made on the first (as a second Run in process makes it).
    state <- omicEnrichW113Frozen$gprofiler$state
    state[["previous"]] <- "enrichment_completed"
    expect_identical(seen$kit113_d6$d6$state, state)
})

test_that("(KIT-113 W2) a forced worker failure commits, records and publishes nothing; then a managed Run gives W1 (b)'s run step", {
    skip_if_not(workflowOwnerCapable())
    log <- new.env(parent = emptyenv())
    fail <- function(ctx) {
        log$fail <- TRUE
        omicEnrichW1Run(1L)(ctx)
        log$fail <- FALSE
        ctx$keep("pending", shiny::isolate(ctx$bus()$state_manager)$.context$pending_provenance)
    }
    run <- omicsEnrichW2Drive(list(fail = fail, run = omicEnrichW1Run(2L)), log)
    seen <- run$captures
    failed <- seen$fail
    expect_identical(failed$notes, c("[default NULL enrichment_working] Running enrichment analysis...",
        "[error 10 -] Error in enrichment analysis: forced worker failure", "remove enrichment_working"))
    for (field in c("tab_status", "history", "states", "product", "resources", "display", "returned")) {
        expect_identical(failed[[field]], "<unchanged>", info = field)
    }
    expect_identical(failed$kernel, character())
    expect_null(failed$pending)
    # D6: the request half has written the annotation table; the filter metadata is the accept half's.
    expect_identical(failed$bus, c(enrichment_analysis_results = "<NULL>", `enrichment_ui_params$` = "<NULL>",
        `enrichment_organism_filter$` = "<NULL>", uniprot_dat_cln = "5x9 13df57c377de46dc"))
    expect_identical(as.integer(failed$files[[omicsEnrichW2StateEntries]]) -
        as.integer(omicsEnrichW2Held(seen, "select", "files")[[omicsEnrichW2StateEntries]]), 1L)
    expect_identical(vapply(run$dispatches, `[[`, "", "acceptance"), c("failed", "complete"))
    frozen <- omicEnrichW1Frozen$gprofiler
    for (field in c("values", "notes", "progress", "tab_status", "history", "states", "product", "resources", "kernel",
        "services", "bus", "display", "returned", "errors")) {
        expect_identical(omicsEnrichW2Held(seen, "run", field), omicsEnrichW2Held(frozen, "run", field), info = field)
    }
    expect_identical(seen$kit113_d6$d6, omicEnrichW113Frozen$gprofiler)
})

test_that("(KIT-113 W2) the kernel travels by name: the double and processEnrichments map both ways, looked up when called; any other is refused", {
    ns <- asNamespace("MultiScholaR")
    for (name in c("runProtEnrichDeterministicProcessEnrichments", "processEnrichments")) {
        expect_identical(workflowEnrichmentKernelName(get(name, envir = ns)), name)
        expect_identical(workflowEnrichmentKernel(name), get(name, envir = ns))
    }
    config <- omicGuiConfig("proteomics")$enrichment
    expect_identical(workflowEnrichmentKernelName(resolveProtEnrichProcessEnrichmentsFn(isTestModeFn = function() TRUE)),
        "runProtEnrichDeterministicProcessEnrichments")
    expect_identical(workflowEnrichmentKernelName(resolveProtEnrichProcessEnrichmentsFn(isTestModeFn = function() FALSE)),
        "processEnrichments")
    # A recorder in the double's binding (KIT-107 W1 (b)'s) is what the selection returns and what the name finds; the
    # double a configuration captured before it no longer matches.
    local({
        recorder <- function(...) "recorder"
        testthat::local_mocked_bindings(runProtEnrichDeterministicProcessEnrichments = recorder, .package = "MultiScholaR")
        selected <- withr::with_options(list(multischolar.test_mode = TRUE), config$kernel$select())
        expect_identical(selected, recorder)
        expect_identical(workflowEnrichmentKernelName(selected), "runProtEnrichDeterministicProcessEnrichments")
        expect_identical(workflowEnrichmentKernel("runProtEnrichDeterministicProcessEnrichments"), recorder)
        expect_error(workflowEnrichmentKernelName(config$kernel$test_double), class = "workflow_enrichment_input")
    })
    expect_error(workflowEnrichmentKernelName(function(...) NULL), class = "workflow_enrichment_input")
    expect_error(workflowEnrichmentKernelName(base::identity), class = "workflow_enrichment_input")
    for (name in list("identity", "base::eval", NULL, NA_character_, c("processEnrichments", "processEnrichments"), 1L)) {
        expect_error(workflowEnrichmentKernel(name), class = "workflow_enrichment_input")
    }
})
