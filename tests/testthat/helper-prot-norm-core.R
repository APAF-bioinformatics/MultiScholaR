# KIT-122: the fixtures test-prot-norm-core.R and test-prot-norm-worker-figures.R share, moved unchanged from
# test-prot-norm-core.R when its worker-figures drives moved to their own file (the file had reached the suite's
# per-file ceiling): the native DIA and DIA_LIMPA owners, the managed session of the tab (each action through the
# worker's entry point in this process), the notices' recorder and the tab's clicks.

protNormCoreQuiet <- function(code) suppressWarnings(suppressMessages(code))

# A native DIA owner (the iq route) on the "small" protein recipe of test-workflow-checkpoint-dia-state.R, with
# the paths and the plain-environment workflow the GUI's actions read.
protNormCoreOwner <- function(envir = parent.frame()) {
    withr::local_options(multischolar.dia_native_sessions = TRUE, multischolar.capture_test_checkpoints = FALSE,
        .local_envir = envir)
    root <- withr::local_tempdir(.local_envir = envir)
    paths <- list(source_dir = file.path(root, "scripts", "proteomics"),
        protein_qc_dir = file.path(root, "results", "proteomics", "protein_qc"),
        peptide_qc_dir = file.path(root, "results", "proteomics", "peptide_qc"))
    paths$publication_graphs_dir <- file.path(root, "results", "proteomics", "publication_graphs")
    paths$time_dir <- file.path(paths$publication_graphs_dir, "filtering_qc")
    for (path in unlist(paths)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    session <- new.env(parent = emptyenv())
    session$userData <- new.env(parent = emptyenv())
    session$onFlush <- session$onFlushed <- session$onSessionEnded <- function(...) function() invisible(NULL)
    original <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- WorkflowState$new(audit_enabled = FALSE)
    workflow$design_matrix <- original@design_matrix
    workflow$config_list <- original@args
    initializeDIANativeImport(workflow, paths, session, "diann", "iq")
    manager <- workflow$state_manager
    withr::defer(closeDIAStateAdapter(manager), envir = envir)
    manager$saveState("protein_replicate_filtered", original, list(), "source")
    list(root = root, paths = paths, workflow = workflow, manager = manager, original = original)
}

protNormW1Sha <- function(x) digest::digest(x, algo = "sha256")

#' Store-local ids and the fixture's temporary root masked in every character value of plain lists.
protNormW1Mask <- function(x, root) {
    if (is.character(x)) {
        x[!is.na(x) & grepl("^[a-f0-9]{64}$", x)] <- "<id>"
        return(gsub(root, "<root>", x, fixed = TRUE))
    }
    if (is.list(x) && !is.object(x)) {
        x[] <- lapply(x, protNormW1Mask, root = root)
        if (!is.null(names(x))) names(x) <- protNormW1Mask(names(x), root)
    }
    x
}

#' A DIA_LIMPA owner (the KIT-053 record fixture of test-workflow-analysis-record.R): the "small" recipe rolled up and
#' quantified by the limpa DPC engine, then the protein QC stage's between-sample normalisation (none) and peptide
#' count filter, committed as the route commits them, with protNormCoreOwner's paths.
protNormW1LimpaOwner <- function(envir = parent.frame()) {
    withr::local_options(multischolar.dia_native_sessions = TRUE, multischolar.capture_test_checkpoints = FALSE,
        .local_envir = envir)
    root <- withr::local_tempdir(.local_envir = envir)
    paths <- list(source_dir = file.path(root, "scripts", "proteomics"),
        protein_qc_dir = file.path(root, "results", "proteomics", "protein_qc"),
        peptide_qc_dir = file.path(root, "results", "proteomics", "peptide_qc"))
    paths$publication_graphs_dir <- file.path(root, "results", "proteomics", "publication_graphs")
    paths$time_dir <- file.path(paths$publication_graphs_dir, "filtering_qc")
    for (path in unlist(paths)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(paths$source_dir, "kit111-limpa", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager), envir = envir)
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType("DIA_LIMPA")
    manager$.context$experiment_paths <- paths
    rolled <- rollUpPrecursorToPeptide(makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small")), core_utilisation = NA)
    quantified <- protNormCoreQuiet(proteinMissingValueImputationLimpa(log2TransformPeptideMatrix(calcPeptideMatrix(rolled)),
        dpc_slope = 0.8, verbose = FALSE, chunk = 1000, engine = "limpa"))
    setDIAStepProvenance(manager, "design_builder")
    manager$saveState("protein_dpc_quant", quantified, list(), "fixture")
    workflow <- list2env(list(state_manager = manager, qc_params = list(protein_qc = list()), config_list = list()))
    protNormCoreQuiet(withDIAStateScope(manager, function() {
        v1 <- computeDIAReferenceProteinNormalisation(manager, getDIAStateReference(manager), list(norm_method = "none"))
        setDIAStepProvenance(manager, "protein_between_sample_normalisation", "whole", "norm_method", NULL,
            values = list(norm_method = "none"))
        finishDIAReferenceProteinStep(workflow, v1, publish = FALSE)
        v2 <- computeDIAReferenceProteinPeptideCountFilter(manager, getDIAStateReference(manager),
            list(num_peptides_per_protein_thresh = 1, num_peptidoforms_per_protein_thresh = 2))
        setDIAStepProvenance(manager, "protein_peptide_count_filter", "whole",
            c("num_peptides_per_protein_thresh", "num_peptidoforms_per_protein_thresh"), NULL, values = v2$config[1:2])
        finishDIAReferenceProteinStep(workflow, v2, publish = FALSE)
    }))
    list(root = root, paths = paths, workflow = workflow, manager = manager)
}

#' One managed action, in this process: the frozen request (with its named input) through the worker's entry point,
#' the acceptance and what the tab reads. A worker failure is raised with the worker's message.
protNormW2Dispatch <- function(fixture, controller, operation, parameters, input, context = NULL) {
    store <- fixture$manager$.context$store
    id <- paste0("wfx-", substr(digest::digest(list(tempfile(), operation), algo = "sha256"), 1L, 32L))
    paths <- workflowExportPaths(store$directory, id)
    token <- list(request_id = id, incarnation = store$identity$incarnation, expected_head = workflowExecutionHead(store)$entity_id,
        generation = controller$control$generation, controller_id = controller$control$controller_id)
    writeWorkflowExecutionMessage(file.path(paths$request_dir, "publication-context.rds"),
        list(experiment_label = "core", experiment_paths = fixture$paths, omic_type = "proteomics"))
    writeWorkflowExportJournal(paths, list(schema = "multischolar.workflow_action.v1", request_id = id, status = "admitting",
        expected_head = token$expected_head, token_hash = digest::digest(token, algo = "sha256")))
    request <- list(schema = 1L, build = workflowExecutionBuildIdentity(), operation = operation,
        input_references = list(state = input), parameters = parameters, token = token, owner = store$identity,
        project_root = dirname(store$directory), directory = paths$request_dir, config = list(mode = "adapters"),
        worker_cap = 12 * 1024^3)
    path <- file.path(paths$request_dir, "request.rds")
    writeWorkflowExecutionMessage(path, request)
    withr::with_envvar(c(MULTISCHOLAR_WORKFLOW_WORKER = "1"), suppressWarnings(suppressMessages(utils::capture.output(try(
        workflowWorkerMain(c(path, digest::digest(file = path, algo = "sha256"))), silent = TRUE)))))
    message <- readRDS(file.path(paths$request_dir, "worker-result.rds"))
    if (!identical(message$status, "complete")) stop(message$error, call. = FALSE)
    controller$job <- list(request = list(token = token, operation = operation, config = list(mode = "adapters"),
        parameters = parameters), paths = paths, completion_context = validateWorkflowOmicsCompletionContext(operation, context))
    receipt <- acceptWorkflowActionResult(controller, message[c("schema", "status", "token", "result")])
    accepted <- resumeWorkflowAction(controller, id)
    c(readManagedOmicsResult(fixture$manager, accepted, NULL), list(acceptance = receipt$status, result = message$result))
}

protNormW3Inputs <- function(mode, ...) utils::modifyList(list(norm_method = "quantile", ruv_mode = mode,
    ruv_percentage = 50, ruv_k = 1, color_variable = "group", ruv_grouping_variable = "group",
    min_pearson_correlation_threshold = 0.5), list(...))

#' The notices shown (`[type duration] text`) and the ids removed, recorded in `log`.
protNormW3Notices <- function(log, root, env = parent.frame()) {
    notify <- function(ui, action = NULL, duration = 5, closeButton = TRUE, id = NULL,
        type = c("default", "message", "warning", "error"), session = NULL) {
        log$notes <- c(log$notes, sprintf("[%s %s] %s", type[[1L]], if (is.null(duration)) "NULL" else format(duration),
            protNormW1Mask(paste(as.character(ui), collapse = " "), root)))
        invisible(id %||% "notification")
    }
    remove <- function(id, session = NULL) {
        log$removed <- c(log$removed, id)
        invisible(id)
    }
    testthat::local_mocked_bindings(showNotification = notify, removeNotification = remove, .package = "shiny", .env = env)
    testthat::local_mocked_bindings(showNotification = notify, removeNotification = remove, .package = "MultiScholaR", .env = env)
}

#' The managed session: the chain mocks, with each action's request run on its named input and recorded in `probe`
#' (`actions`: operation and input entity id; `values`: the worker's values). `probe$hold` keeps the action running.
protNormW3Managed <- function(fixture, probe, env = parent.frame()) {
    controller <- omicsManagedChainMocks(list(owner = fixture$manager, paths = fixture$paths), probe, env = env)
    testthat::local_mocked_bindings(
        runWorkflowAction = function(operation, input_references, parameters, owner, reservation, completion_context = NULL) {
            probe$actions <- c(probe$actions, list(list(operation = operation, input = input_references$state$entity_id)))
            out <- protNormW2Dispatch(fixture, controller, operation, parameters, input_references$state, completion_context)
            probe$values <- c(probe$values, list(out$result))
            probe$handoff <- list(accepted = list(operation = operation, reference = out$reference,
                state_reference = out$stateReference, display = out$display))
            paste0("worker-", operation)
        },
        pollWorkflowAction = function(owner) list(status = if (isTRUE(probe$hold)) "running" else "complete"),
        .package = "MultiScholaR", .env = env)
    controller
}

#' A click on `button` (after setting `inputs`), counted per drive.
protNormW3Clicker <- function() {
    counts <- new.env(parent = emptyenv())
    function(button, ...) {
        inputs <- list(...)
        function(session) {
            if (length(inputs)) do.call(session$setInputs, inputs)
            counts[[button]] <- (counts[[button]] %||% 0) + 1
            do.call(session$setInputs, stats::setNames(list(counts[[button]]), button))
        }
    }
}
