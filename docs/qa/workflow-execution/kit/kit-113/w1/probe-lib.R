# KIT-113 W1 probes: the shared emulation (scratch, not frozen). Sourced by probe-p1.R, probe-p2-p5.R and probe-p4.R
# after `pkgload::load_all(<tree>, export_all = TRUE, helpers = TRUE)` with the working directory at <tree>.
#
# The managed halves do not exist on the start commit (132e5445). They are emulated here with today's statements of
# `runWorkflowEnrichment` (R/utils_workflow_enrichment.R:822-925) and `runDIAReferenceEnrichment`
# (R/utils_workflow_dia_da.R:656-685), grouped as KIT-113 D3 says:
#   request half (controller): the context (:824-845), the annotation resolve (:853-854), the mapping and its warning,
#     the method (:877-880), the kernel's name; the request's parameters hold references, never a DA table;
#   handler (worker, `computeWorkflowDiaEnrichmentAction`, defined in the global environment, which the registry's
#     `get(entry$handler, mode = "function")` reaches from the namespace): the source's components, the DA tables with the
#     refusal (:846-852), the match (:855-856), the filter's application, the kernel's arguments and the run (:881-890,
#     `runDIAReferenceEnrichment` whose records are a no-op on a worker), the results made portable and saved as an `rds`
#     state value;
#   accept half (controller, in the owner's scope): the completion's provenance cleared, the results hydrated, the
#     external-resource records (:677-683), the bus's filter metadata (:874-876), the display and the commit (:891-925).
# The registry entry is replaced in this process only (`registerWorkflowOperation(replace = TRUE)`), as D3 registers it.
# Nothing under R/ is edited.

say <- function(...) cat(sprintf(...), "\n", sep = "")

#' The W1 (b) helpers and frozen captures (every top-level assignment of the test file), as test-omics-enrich-core.R
#' loads them; and that file's own helpers.
loadKit113TestHelpers <- function(tree) {
    for (file in c("test-omic-enrich-module.R", "test-omics-enrich-core.R")) {
        for (expression in parse(file.path(tree, "tests", "testthat", file), keep.source = FALSE)) {
            if (is.call(expression) && identical(expression[[1L]], as.name("<-"))) eval(expression, globalenv())
        }
    }
}

#' D3's registry entry for `enrichment`, in this process only.
registerKit113EnrichmentEntry <- function() {
    declaration <- diaEnrichmentOperation()
    entry <- list(name = declaration$name, input = "checkpoint", input_class = declaration$input_class,
        compute = "omics", handler = "computeWorkflowDiaEnrichmentAction", validate = "validateWorkflowOmicsParameters",
        resolve = NULL, parameter_limit = 65536, parameters = declaration$parameters %||% character(),
        completion = "omics", publication = "omics", step_label = NULL, omic_type = "proteomics",
        variants = "whole", default_variant = "whole", diagnostic_only = FALSE, inclusion_absence = FALSE,
        named_input = TRUE, declaration = declaration, section = NULL)
    registerWorkflowOperation(entry, replace = TRUE)
    invisible(workflowActionSpec("enrichment"))
}

#' KIT-111 W2's mocks (`protNormW2Mocks`, test-prot-norm-core.R:1302-1306): the build identity fixed and the service
#' teardown replaced (no systemd unit exists in this process), for the life of `env`.
kit113Mocks <- function(env = parent.frame()) {
    testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
        finishWorkflowExecutionJob = function(owner, status) { owner$job <- NULL; list(status = status) },
        .package = "MultiScholaR", .env = env)
}

#' Save a transfer value: a KIT-019 display table (`table`) or an `rds` state value (`rds`); NULL stays NULL.
kit113Save <- function(manager, value, codec) {
    if (is.null(value)) return(NULL)
    if (identical(codec, "table")) return(saveDIADisplayTable(manager, value))
    saveWorkflowStateValue(manager$.context, value, "kit113-transfer", codec = "rds")
}

kit113Load <- function(manager, reference) {
    if (is.null(reference)) return(NULL)
    if (identical(reference$payload$codec, "parquet_duckdb")) return(materializeWorkflowTable(manager$.context, reference))
    resolveWorkflowStateValue(manager$.context, reference)
}

#' The request half (controller). Returns the context for the accept half and the request's parameters.
kit113Request <- function(workflow, config, manifest, parameters, paths, notify, progress, codec = "rds") {
    notify <- notify %||% function(message, type, duration) invisible(NULL)
    progress <- progress %||% function(value, detail) invisible(NULL)
    if (!isTRUE(config$da_source$manifest(manifest))) {
        rlang::abort("Enrichment needs an accepted DA result manifest.", class = "workflow_enrichment_input")
    }
    manager <- workflow$state_manager
    source_object <- readDIAExecutionComponents(manager, manifest$source, fields = config$da_source$components)
    selected <- parameters$selected_contrast
    progress(0.2, "Transforming DE data...")
    progress(0.3, "Creating DE results S4 object...")
    contrasts_tbl <- resolveProtEnrichContrastsTable(workflow)
    selection <- resolveProtEnrichSelectedDaResults(selected, manifest$results, contrasts_tbl)
    if (is.null(selection$selectedDaResults)) {
        rlang::abort(sprintf("Could not find DE results for contrast '%s'. Available contrasts: %s", selected,
            paste(selection$availableKeys, collapse = ", ")), class = "workflow_enrichment_input")
    }
    dependencies <- resolveProtEnrichRunDependencies(source_object, manifest$results, workflow, contrasts_tbl)
    contrasts_tbl <- dependencies$contrastsTbl
    if (is.null(contrasts_tbl) || is.null(dependencies$designMatrix)) {
        rlang::abort(paste("Missing contrasts_tbl or design_matrix required for enrichment analysis.",
            "Please ensure DE analysis was completed successfully."), class = "workflow_enrichment_input")
    }
    directories <- resolveProtEnrichOutputDirectories(paths)
    annotations <- config$annotation$resolve(workflowData = workflow, experimentPaths = paths,
        currentS4Object = source_object, organismTaxid = parameters$organism_taxid)$uniprotDatCln
    mapping <- NULL
    if (isTRUE(parameters$enable_organism_filter) && !is.null(config$organism_filter)) {
        target <- as.character(parameters$organism_taxid)
        mapping <- config$organism_filter$mapping(workflowData = workflow, uniprotDatCln = annotations,
            targetTaxon = target)$organismMapping
        if (is.null(mapping) || nrow(mapping) == 0) {
            notify("Multi-species filtering enabled but no organism mapping available. Proceeding without filtering.",
                "warning", 8)
            mapping <- NULL
        }
    }
    progress(0.5, "Running enrichment analysis...")
    method <- config[["organism"]]$method(organismTaxid = parameters$organism_taxid,
        supportedOrganisms = config[["organism"]]$supported())$methodInfo
    kernel <- config$kernel$select()
    ns <- asNamespace("MultiScholaR")
    kernel_name <- if (identical(kernel, get("runProtEnrichDeterministicProcessEnrichments", envir = ns))) {
        "runProtEnrichDeterministicProcessEnrichments"
    } else if (identical(kernel, get("processEnrichments", envir = ns))) {
        "processEnrichments"
    } else rlang::abort("unknown enrichment kernel", class = "workflow_enrichment_input")
    head_manifest <- manager$.context$snapshot$da_results
    if (!identical(resolveWorkflowStateValue(manager$.context, head_manifest), manifest))
        rlang::abort("the held manifest is not the head's", class = "workflow_enrichment_input")
    request <- list(manifest = head_manifest, contrasts_tbl = contrasts_tbl, inputs = parameters,
        directories = list(daOutputDir = directories$daOutputDir, pathwayDir = directories$pathwayDir),
        annotations = kit113Save(manager, annotations, codec), mapping = kit113Save(manager, mapping, codec),
        method = method, kernel = kernel_name)
    list(context = list(source_object = source_object, selected = selected, selection = selection,
        contrasts_tbl = contrasts_tbl, directories = directories, annotations = annotations, method = method),
        parameters = request)
}

#' The worker handler (D3), defined where the registry's name lookup reaches it.
computeWorkflowDiaEnrichmentAction <- function(manager, request, variant) {
    p <- request$parameters
    manifest <- resolveWorkflowStateValue(manager$.context, p$manifest)
    if (!identical(manifest$source$entity_id, request$input_references$state$entity_id))
        cli::cli_abort("enrichment input is not the manifest's source", class = "workflow_execution_stale")
    config <- omicsEnrichCoreConfig(omicGuiConfig("proteomics")$enrichment)
    workflow <- list2env(list(state_manager = manager))
    source_object <- readDIAExecutionComponents(manager, manifest$source, fields = config$da_source$components)
    dependencies <- resolveProtEnrichRunDependencies(source_object, manifest$results, workflow, p$contrasts_tbl)
    da_results <- config$da_source$read(contrasts_tbl = dependencies$contrastsTbl,
        design_matrix = dependencies$designMatrix, da_output_dir = p$directories$daOutputDir)
    missing <- names(Filter(is.null, da_results@da_data))
    if (length(missing)) {
        rlang::abort(sprintf("The DA results table is missing for contrast %s.", paste(missing, collapse = ", ")),
            class = "workflow_enrichment_input")
    }
    annotations <- kit113Load(manager, p$annotations)
    matching <- config$annotation$match(uniprotDatCln = annotations, daResultsForEnrichment = da_results,
        currentS4Object = source_object)
    filter_applied <- FALSE
    filter_stats <- list(proteins_before = 0, proteins_after = 0, proteins_removed = 0)
    mapping <- kit113Load(manager, p$mapping)
    if (!is.null(mapping)) {
        filtered <- config$organism_filter$apply(daResultsForEnrichment = da_results, organismMapping = mapping,
            targetTaxon = as.character(p$inputs$organism_taxid), currentS4Object = source_object)
        da_results <- filtered$daResultsForEnrichment
        filter_applied <- filtered$filterApplied
        filter_stats <- filtered$filterStats
    }
    columns <- resolveProtEnrichAnalysisInputColumns(methodInfo = p$method, daResultsForEnrichment = da_results,
        currentS4Object = source_object)
    arguments <- buildProtEnrichProcessEnrichmentsArgs(daResultsForEnrichment = da_results,
        organismTaxid = p$inputs$organism_taxid, upCutoff = p$inputs$up_cutoff, downCutoff = p$inputs$down_cutoff,
        qCutoff = p$inputs$q_cutoff, pathwayDir = p$directories$pathwayDir, goAnnotations = annotations,
        proteinIdColumn = columns$idColumn, contrastNames = names(manifest$results),
        correctionMethod = p$inputs$correction_method, excludeIea = config$exclude_iea)
    kernel <- get(p$kernel, envir = asNamespace("MultiScholaR"))
    declared <- runDIAReferenceEnrichment(workflow, manifest, arguments$processArgs, process_fn = kernel)
    results <- saveWorkflowStateValue(manager$.context, portableWorkflowEnrichmentResults(declared$results),
        "kit113-results", codec = "rds")
    list(checkpoint = readWorkflowStateMetadata(manager, request$input_references$state)$data,
        reference = declared$reference, commit = NULL,
        display = list(results = results, service = declared$service, parameters = declared$parameters,
            matching = matching$annotationMatchResults,
            filter = list(applied = filter_applied, stats = filter_stats), contrasts = length(manifest$results),
            taxon = arguments$processArgs$taxon_id))
}

#' One managed action through the worker's entry point and the controller's acceptance, in this process, with the
#' request's named input (KIT-111 W2's `protNormW2Dispatch`, test-prot-norm-core.R:1176-1200).
kit113Dispatch <- function(manager, controller, paths, operation, parameters, input, context) {
    store <- manager$.context$store
    id <- paste0("wfx-", substr(digest::digest(list(tempfile(), operation), algo = "sha256"), 1L, 32L))
    export <- workflowExportPaths(store$directory, id)
    token <- list(request_id = id, incarnation = store$identity$incarnation, expected_head = workflowExecutionHead(store)$entity_id,
        generation = controller$control$generation, controller_id = controller$control$controller_id)
    writeWorkflowExecutionMessage(file.path(export$request_dir, "publication-context.rds"),
        list(experiment_label = "kit113", experiment_paths = paths, omic_type = "proteomics"))
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
    list(message = message, receipt = receipt, export = export, request = request)
}

#' The accept half (controller), returning today's list.
kit113Accept <- function(workflow, config, context, delivered, parameters, progress) {
    progress <- progress %||% function(value, detail) invisible(NULL)
    manager <- workflow$state_manager
    withDIAStateScope(manager, function() {
        clearWorkflowStepProvenance(manager)
        display <- delivered$display
        results <- resolveWorkflowStateValue(manager$.context, display$results)
        for (service in setdiff(display$service, "unavailable")) {
            parts <- strsplit(service, ":", fixed = TRUE)[[1L]]
            recordWorkflowExternalResource(manager, parts[[1L]], paste(parts[-1L], collapse = ":"),
                if (identical(parts[[1L]], "clusterProfiler")) "package_version" else "database_version",
                query_summary = sprintf("enrichment of %d contrasts, taxon %s", display$contrasts, display$taxon %||% NA))
        }
        organism_filter <- persistProtEnrichOrganismFilterMetadata(workflowData = workflow,
            organismFilterEnabled = parameters$enable_organism_filter, organismFilterApplied = display$filter$applied,
            targetTaxonId = parameters$organism_taxid, filterStats = display$filter$stats)
        selected <- context$selected
        method <- context$method
        contrasts_tbl <- context$contrasts_tbl
        source_object <- context$source_object
        directories <- context$directories
        contrast_results <- buildProtEnrichAllContrastResults(enrichmentResults = results, methodInfo = method)
        for (contrast in names(contrast_results)) {
            contrast_results[[contrast]]$gprofiler_results_table <- saveDIADisplayTable(manager,
                contrast_results[[contrast]]$gprofiler_results)
            contrast_results[[contrast]]$clusterprofiler_results_table <- saveDIADisplayTable(manager,
                contrast_results[[contrast]]$clusterprofiler_results)
        }
        shown <- resolveProtEnrichSelectedContrastResults(selected, contrast_results, contrasts_tbl)
        progress(0.8, "Storing results...")
        workflow$enrichment_analysis_results <- buildProtEnrichAnalysisResultsPayload(
            gprofilerResults = if (shown$found) shown$gprofilerResults,
            clusterprofilerResults = if (shown$found) shown$clusterprofilerResults,
            stringdbResults = if (shown$found) shown$stringdbResults, fullEnrichmentResults = results,
            selectedContrast = selected, analysisMethod = method$method, upCutoff = parameters$up_cutoff,
            downCutoff = parameters$down_cutoff, qCutoff = parameters$q_cutoff, organismTaxid = parameters$organism_taxid)
        stamped <- propagateProtEnrichResultsArgs(enrichmentResults = results, currentS4Object = source_object,
            selectedContrast = selected, methodInfo = method, upCutoff = parameters$up_cutoff,
            downCutoff = parameters$down_cutoff, qCutoff = parameters$q_cutoff, organismTaxid = parameters$organism_taxid,
            pathwayDir = directories$pathwayDir)$enrichmentResults
        ui_parameters <- propagateProtEnrichUiParams(currentS4Object = source_object, workflowData = workflow,
            selectedContrast = selected, methodInfo = method, upCutoff = parameters$up_cutoff,
            downCutoff = parameters$down_cutoff, qCutoff = parameters$q_cutoff, organismTaxid = parameters$organism_taxid)
        updateProtEnrichStateManagerUiParams(workflowData = workflow, storedUiParams = ui_parameters$storedUiParams)
        saved <- saveWorkflowEnrichmentState(workflow, config, stamped, display$parameters, method, parameters,
            directories$pathwayDir)
        status <- workflow$tab_status
        status[[config$status_key]] <- "complete"
        workflow$tab_status <- status
        progress(1.0, "Complete!")
        list(results = results, method = method, contrasts = contrast_results, selected = shown,
            raw_contrast = context$selection$rawContrastName, contrasts_tbl = contrasts_tbl,
            pathway_dir = directories$pathwayDir, annotation_match = display$matching, organism_filter = organism_filter,
            service = display$service, reference = delivered$reference, parameters = display$parameters, saved = saved)
    })
}

#' The pathway directory's files: name, sha256 of the bytes, and for a table the sha256 of its lines with the fixture's
#' root masked (as the entry suite digests them).
kit113PathwayFiles <- function(dir, root) {
    files <- sort(list.files(dir, recursive = TRUE, all.files = TRUE, no.. = TRUE))
    vapply(files, function(file) {
        path <- file.path(dir, file)
        raw <- digest::digest(file = path, algo = "sha256")
        if (grepl("\\.(tsv|tab|csv)$", file, ignore.case = TRUE)) {
            masked <- digest::digest(guiEntryMask(readLines(path, warn = FALSE), root), algo = "sha256", serialize = FALSE)
            paste0("raw ", substr(raw, 1L, 16L), " masked ", substr(masked, 1L, 16L))
        } else paste0("raw ", substr(raw, 1L, 16L))
    }, "")
}

#' Serialized and in-memory sizes of a value.
kit113Sizes <- function(value) c(serialized = length(serialize(value, NULL)), object_size = as.numeric(object.size(value)))

#' `runWorkflowAction`'s guards (R/utils_workflow_execution.R:1515-1600), in its order, on a real request; the
#' application's own checks (its liveness, the busy check, the host headroom, the admission lock, the systemd units) are
#' not run. Each guard's outcome is "ok" or the refusal's class and message.
kit113Guards <- function(manager, operation, input_references, parameters, completion_context) {
    policy <- workflowExecutionPolicy()
    store <- manager$.context$store
    run <- function(code) tryCatch({ force(code); "ok" }, error = function(e)
        paste0("REFUSED [", paste(setdiff(class(e), c("rlang_error", "error", "condition")), collapse = "/"), "] ",
            conditionMessage(e)))
    c(`:1536 assertDIAExportReady` = run(assertDIAExportReady(manager)),
        `:1540 assertWorkflowOperationDispatchable` = run(assertWorkflowOperationDispatchable(operation)),
        `:1541 validateWorkflowActionParameters (resolved = FALSE)` = run(validateWorkflowActionParameters(operation,
            parameters, policy, resolved = FALSE)),
        `:1542-1544 input reference shape` = run(if (!is.list(input_references) ||
            !identical(names(input_references), "state") || !isSavedWorkflowReference(input_references$state))
            stop("unsupported input reference")),
        `:1565 validateSavedWorkflowReference (payload = FALSE)` = run(validateSavedWorkflowReference(store,
            input_references$state, payload = FALSE)),
        `:1566 assertWorkflowActionInputSelected` = run(assertWorkflowActionInputSelected(operation, manager,
            input_references)),
        `:1567 validateWorkflowActionInput` = run(validateWorkflowActionInput(operation, manager, input_references)),
        `:1568 resolveWorkflowActionParameters` = run(resolveWorkflowActionParameters(operation, manager,
            input_references$state, parameters, policy)),
        `:1586-1587 validateWorkflowOmicsCompletionContext` = run(validateWorkflowOmicsCompletionContext(operation,
            completion_context)))
}

#' The acceptance's checks on a worker value (R/utils_workflow_execution.R:1771-1774), each separately.
kit113ValueChecks <- function(store, value) {
    c(`:1771 !hasMutableWorkflowEntityOwner(value)` = !hasMutableWorkflowEntityOwner(value),
        `:1771 isDIAComponentCheckpoint(value$checkpoint)` = isDIAComponentCheckpoint(value$checkpoint),
        `:1772 isSavedWorkflowReference(value$reference)` = isSavedWorkflowReference(value$reference),
        `:1772 object.size(value) <= 8388608` = as.numeric(object.size(value)) <= 8388608,
        `:1774 validateSavedWorkflowReference(payload = TRUE)` = isTRUE(tryCatch({
            validateSavedWorkflowReference(store, value$reference, payload = TRUE)
            TRUE
        }, error = function(e) FALSE)))
}
