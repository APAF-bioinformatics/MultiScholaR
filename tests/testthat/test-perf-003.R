# PERF-102 focused correctness. Source loading here is explicitly not installed-runtime proof.
perf102Code <- function(science = FALSE) {
    parent <- if (science) asNamespace("MultiScholaR") else baseenv()
    code <- new.env(parent = parent)
    sys.source(testthat::test_path("..", "..", "R/utils_workflow_execution.R"), code)
    if (science) for (name in c("newWorkflowEntityStore", "newDIAStateAdapter",
        "closeDIAStateAdapter", "WorkflowState", "getDIAStateReference",
        "materializeDIATable", "hasMutableWorkflowEntityOwner",
        "workflowExportPaths", "writeWorkflowExportJournal"))
        code[[name]] <- get(name, envir = parent, inherits = TRUE)
    code
}

perf102Evidence <- function(code, goal = "benefit") {
    keys <- c("operation", "input", "parameters", "runtime", "boundary", "policy", "outputs")
    bindings <- stats::setNames(as.list(rep(strrep("a", 64), length(keys))), keys)
    bindings$policy <- digest::digest(code$workflowPerfQualificationPolicy(), algo = "sha256")
    control <- list(status = "complete", peak_bytes = 6 * 1024^3, elapsed_seconds = 100,
        swap_bytes = 0, spill_bytes = 0, science = TRUE, outputs = TRUE,
        measurement_valid = TRUE, integrity = TRUE, guard = TRUE, admission = TRUE)
    candidate <- control
    candidate$peak_bytes <- 3 * 1024^3
    candidate$elapsed_seconds <- 115
    pair <- list(control = control, candidate = candidate)
    evidence <- list(schema = 2L, policy = code$workflowPerfQualificationPolicy()$id, goal = goal, bindings = bindings,
        domain_supported = TRUE, runtime_eligible = TRUE, open_holds = character(),
        release_holds = character(), populations = list(primary = list(pair, pair),
            wide = list(pair, pair)), arms = list(native = list(control, control),
            physical = list(candidate, candidate), selection = list(candidate, candidate)))
    plan <- list(schema = 1L, goal = goal,
        finalist = if (goal == "benefit") "physical" else NULL,
        screen_decision = if (goal == "benefit") digest::digest("accepted-screen", algo = "sha256") else NULL,
        members = list())
    stamp <- function(member, population, role, rep) {
        arm <- if (role == "control") "native" else if (role == "candidate") "physical" else role
        slot <- paste(population, role, rep, sep = ".")
        identity <- list(run_id = slot, population = population, arm = arm, repetition = rep,
            input = digest::digest(population, algo = "sha256"), parameters = bindings$parameters,
            runtime = digest::digest(paste(arm, "runtime"), algo = "sha256"),
            build = digest::digest(paste(arm, "build"), algo = "sha256"),
            boundary = bindings$boundary, output_contract = bindings$outputs)
        plan$members[[slot]] <<- identity
        c(member, identity, list(receipt_id = digest::digest(paste(slot, "receipt"), algo = "sha256")))
    }
    if (goal == "screen_feasibility") {
        for (arm in names(evidence$arms)) for (rep in seq_len(2L))
            evidence$arms[[arm]][[rep]] <- stamp(evidence$arms[[arm]][[rep]], "screen", arm, rep)
    } else {
        for (population in names(evidence$populations)) for (rep in seq_len(2L))
            for (arm in c("control", "candidate")) evidence$populations[[population]][[rep]][[arm]] <-
                stamp(evidence$populations[[population]][[rep]][[arm]], population, arm, rep)
    }
    evidence$execution_plan <- plan
    evidence$bindings$execution_plan <- digest::digest(plan, algo = "sha256")
    evidence
}

test_that("historical policy replay is unchanged and PERF requires both populations", {
    code <- perf102Code()
    manifest <- jsonlite::read_json(testthat::test_path("..", "testdata/workflow-execution/qualification.json"))
    keys <- c("operation", "input", "parameters", "runtime", "boundary", "policy", "outputs")
    bindings <- stats::setNames(as.list(rep(strrep("a", 64), length(keys))), keys)
    for (fixture in manifest$decision_fixtures$cases) {
        fixture$schema <- 1L
        fixture$bindings <- bindings
        fixture$open_holds <- as.character(unlist(fixture$open_holds, use.names = FALSE))
        expect_identical(code$evaluateWorkflowQualification(fixture, bindings)$status,
            fixture$expected, info = fixture$id)
    }
    historical <- Filter(function(x) identical(x$expected, "benefit") &&
        identical(x$goal, "benefit"), manifest$decision_fixtures$cases)[[1L]]
    historical$schema <- 1L
    historical$bindings <- bindings
    historical$open_holds <- character()
    for (i in 1:2) {
        historical$pairs[[i]]$candidate$peak_bytes <- 0.85 * historical$pairs[[i]]$control$peak_bytes
        historical$pairs[[i]]$candidate$elapsed_seconds <- 1.15 * historical$pairs[[i]]$control$elapsed_seconds
    }
    expect_identical(code$evaluateWorkflowQualification(historical, bindings)$status, "benefit")
    evidence <- perf102Evidence(code)
    malformed <- evidence
    malformed$bindings <- c(malformed$bindings, stats::setNames(list(strrep("a", 64)), NA_character_))
    expect_identical(code$evaluateWorkflowQualification(malformed, malformed$bindings)$status, "held")
    evaluate <- function(x) code$evaluateWorkflowQualification(x, evidence$bindings)
    expect_identical(evaluate(evidence)$status, "benefit")
    for (population in c("primary", "wide")) for (i in 1:2) {
        x <- evidence
        x$populations[[population]][[i]]$candidate$peak_bytes <- 0.85 * 6 * 1024^3
        expect_identical(evaluate(x)$status, "inconclusive")
        x$populations[[population]][[i]]$candidate$peak_bytes <- 3 * 1024^3
        x$populations[[population]][[i]]$candidate$elapsed_seconds <- 121
        expect_identical(evaluate(x)$status, "inconclusive")
    }
    x <- evidence
    x$populations$wide <- NULL
    expect_identical(evaluate(x)$status, "held")
    x <- evidence
    x$populations$primary[[1]]$control$status <- "memory_guard"
    x$populations$primary[[1]]$control$peak_bytes <- NA_real_
    expect_identical(evaluate(x)$status, "held")
    expect_null(evaluate(x)$ratios)
    for (counter in c(NA_real_, NaN, Inf, 0, -1)) {
        x <- evidence
        x$populations$wide[[2]]$control$peak_bytes <- counter
        expect_identical(evaluate(x)$status, "held")
    }
    for (failure in c("science", "outputs", "measurement_valid", "integrity", "guard", "admission")) {
        x <- evidence
        x$populations$wide[[2]]$candidate[[failure]] <- FALSE
        expect_identical(evaluate(x)$status, "rejected")
    }
    for (field in c("swap_bytes", "spill_bytes")) {
        x <- evidence
        x$populations$primary[[2]]$candidate[[field]] <- 1
        expect_identical(evaluate(x)$status, "rejected")
    }
    x <- evidence
    x$release_holds <- "unsigned runtime"
    expect_identical(evaluate(x)$status, "held")
    expect_identical(evaluate(x)$comparison_status, "benefit")
    expect_false(evaluate(x)$promotion_authorized)
    x <- evidence
    x$bindings$policy <- strrep("b", 64)
    expect_identical(evaluate(x)$status, "held")
    x <- evidence
    x$schema <- 3L
    expect_identical(evaluate(x)$status, "held")
})

test_that("PERF screening ranks feasibility without imposing full-size benefit", {
    code <- perf102Code()
    evidence <- perf102Evidence(code, "screen_feasibility")
    for (arm in c("physical", "selection")) for (i in 1:2) {
        evidence$arms[[arm]][[i]]$peak_bytes <- 5.8 * 1024^3
        evidence$arms[[arm]][[i]]$elapsed_seconds <- 140
    }
    evaluate <- function(x) code$evaluateWorkflowQualification(x, evidence$bindings)
    expect_identical(evaluate(evidence)$finalist, "physical")
    x <- evidence
    x$arms$selection[[1]]$elapsed_seconds <- 130
    expect_identical(evaluate(x)$finalist, "selection")
    x <- evidence
    x$arms$selection[[1]]$peak_bytes <- 5 * 1024^3
    x$arms$selection[[2]]$peak_bytes <- 5 * 1024^3
    expect_identical(evaluate(x)$finalist, "selection")
    for (failure in c("science", "outputs", "measurement_valid", "guard", "admission")) {
        x <- evidence
        x$arms$physical[[2]][[failure]] <- FALSE
        expect_identical(evaluate(x)$finalist, "selection")
        x$arms$selection[[1]][[failure]] <- FALSE
        expect_identical(evaluate(x)$status, "no_finalist")
    }
    x <- evidence
    x$release_holds <- "private API"
    x$runtime_eligible <- FALSE
    expect_identical(evaluate(x)$status, "screen_finalist")
    expect_false(evaluate(x)$promotion_authorized)
    x <- evidence
    x$arms$physical[[1]]$peak_bytes <- 15.5 * 1024^3
    expect_identical(evaluate(x)$finalist, "selection")
    expect_identical(code$evaluateWorkflowQualification(evidence, evidence$bindings,
        list(memory_max = 17 * 1024^3, deadline = 300))$status, "held")
})

test_that("registered sample parameters preserve native precedence without assay reads", {
    code <- perf102Code(TRUE)
    root <- normalizePath(testthat::test_path("..", ".."))
    sys.source(file.path(root, "tests/testthat/helper-workflow-checkpoint-fixtures.R"), code)
    recipe <- code$readWorkflowCheckpointRecipe("small", root)
    recipe$n_proteins <- 8L
    recipe$n_samples <- 8L
    object <- code$makeWorkflowCheckpointPeptide(recipe)
    object@args$filterMinNumPeptidesPerSample <- list(peptides_per_sample_cutoff = 2L,
        inclusion_list = c(" a ", "a", ""))
    before <- digest::digest(object, algo = "sha256")
    # The seam supplies real native prototypes, proving that no payload field is requested.
    fields <- NULL
    code$readDIAExecutionComponents <- function(manager, reference, fields) {
        stopifnot(identical(fields, character()))
        prototype <- object
        prototype@peptide_data <- prototype@peptide_data[0, , drop = FALSE]
        prototype@design_matrix <- prototype@design_matrix[0, , drop = FALSE]
        prototype
    }
    policy <- code$workflowExecutionPolicy(TRUE)
    for (parameters in list(list(), list(peptides_per_sample_cutoff = NULL, inclusion_list = NULL),
        list(peptides_per_sample_cutoff = 3L, inclusion_list = character()),
        list(peptides_per_sample_cutoff = 1, inclusion_list = c(" a ", "a", "")))) {
        expected <- do.call(get("filterMinNumPeptidesPerSample", asNamespace("MultiScholaR")),
            c(list(theObject = object, core_utilisation = NA), parameters))
        frozen <- code$resolveWorkflowActionParameters("minimum_peptides_per_sample", NULL,
            NULL, parameters, policy)
        expected_args <- expected@args$filterMinNumPeptidesPerSample
        expect_identical(frozen, expected_args[c("peptides_per_sample_cutoff", "inclusion_list")])
        parameters$peptides_per_sample_cutoff <- 999L
        expect_false(identical(frozen$peptides_per_sample_cutoff, 999L))
    }
    object@args$filterMinNumPeptidesPerSample <- NULL
    frozen <- code$resolveWorkflowActionParameters("minimum_peptides_per_sample", NULL, NULL, list(), policy)
    expect_equal(frozen$peptides_per_sample_cutoff, 500)
    expect_identical(frozen$inclusion_list, character())
    expect_error(code$resolveWorkflowActionParameters("minimum_peptides_per_sample", NULL, NULL,
        list(min_peptides = 3L), policy), "unknown parameter")
    for (invalid in list(0, -1, NA_real_, Inf, 1.5)) expect_error(
        code$resolveWorkflowActionParameters("minimum_peptides_per_sample", NULL, NULL,
            list(peptides_per_sample_cutoff = invalid), policy), "positive integer")
    expect_error(code$validateWorkflowActionParameters("minimum_peptides_per_sample",
        list(peptides_per_sample_cutoff = 2L, inclusion_list = new.env()), policy), "unsupported execution parameters")
    # Restore the deliberately edited fixture args and check scientific input bytes.
    object@args$filterMinNumPeptidesPerSample <- list(peptides_per_sample_cutoff = 2L,
        inclusion_list = c(" a ", "a", ""))
    expect_identical(digest::digest(object, algo = "sha256"), before)
})

test_that("finite input and completion handlers preserve controller-only ownership", {
    code <- perf102Code(TRUE)
    expect_error(code$workflowActionSpec("arbitrary_R_function"), class = "workflow_execution_unsupported")
    expect_error(code$workflowActionSpec(function() NULL), class = "workflow_execution_unsupported")
    record <- list(class = "PeptideQuantitativeData", components = list(peptide_data = list(), design_matrix = list()),
        slots = list(sample_id = "Run", protein_id_column = "Protein", peptide_sequence_column = "Peptide"),
        dimensions = list(columns = c("Run", "Protein", "Peptide")))
    code$isSavedWorkflowReference <- function(x) identical(x, "owned-reference")
    code$isDIAComponentCheckpoint <- function(x) is.list(x) && !is.null(x$components)
    code$readDIAStateMetadata <- function(...) list(data = record)
    expect_null(code$validateWorkflowActionInput("minimum_peptides_per_sample", NULL, list(state = "owned-reference")))
    record$slots$sample_id <- "missing"
    expect_error(code$validateWorkflowActionInput("minimum_peptides_per_sample", NULL,
        list(state = "owned-reference")), "roles")
    record$slots$sample_id <- "Run"
    record$class <- "ProteinQuantitativeData"
    expect_error(code$validateWorkflowActionInput("minimum_peptides_per_sample", NULL,
        list(state = "owned-reference")), "roles")
    expect_error(code$validateWorkflowActionInput("minimum_peptides_per_sample", NULL,
        list(file = "owned-reference")), "unsupported input")
    owner <- list(accepting = "token", manager = "controller-manager")
    observed <- NULL
    code$finishDIAReferencePeptideStep <- function(workflow, ...) {
        observed <<- list(manager = workflow$state_manager, arguments = list(...)); "accepted"
    }
    value <- list(checkpoint = list(), state_name = "sample_filtered", description = "filter",
        config = list(), section = "filterMinNumPeptidesPerSample", counts = NULL,
        summary = list(inclusion_entries_absent_from_data = "absent"))
    expect_identical(code$completeWorkflowAction(owner, "minimum_peptides_per_sample", value), "accepted")
    expect_identical(observed$manager, "controller-manager")
    expect_identical(observed$arguments$observed_inclusion_absence, "absent")
    expect_error(code$completeWorkflowAction(list(accepting = NULL), "minimum_peptides_per_sample", value), "lease")
    withr::with_envvar(c(MULTISCHOLAR_WORKFLOW_WORKER = "1"), expect_error(
        code$completeWorkflowAction(owner, "minimum_peptides_per_sample", value), "workers cannot"))
})

test_that("registered dispatch freezes inputs and preserves acceptance and output-only recovery", {
    code <- perf102Code(TRUE)
    root <- normalizePath(testthat::test_path("..", ".."))
    sys.source(file.path(root, "tests/testthat/helper-workflow-checkpoint-fixtures.R"), code)
    recipe <- code$readWorkflowCheckpointRecipe("small", root)
    recipe$n_proteins <- recipe$n_samples <- 8L
    object <- code$makeWorkflowCheckpointPeptide(recipe)
    object@args$filterMinNumPeptidesPerSample <- list(peptides_per_sample_cutoff = 3L,
        inclusion_list = character())
    object@peptide_data[[object@raw_quantity_column]][1:3] <- c(0, NA_real_, NaN)
    last_run <- tail(unique(object@peptide_data[[object@sample_id]]), 1L)
    keep <- object@peptide_data[[object@sample_id]] != last_run
    keep[which(!keep)[1L]] <- TRUE
    object@peptide_data <- object@peptide_data[keep, , drop = FALSE]
    temporary <- withr::local_tempdir()
    project <- file.path(temporary, "project")
    dir.create(project)
    store <- code$newWorkflowEntityStore(project, "perf102", "DIA")
    manager <- code$newDIAStateAdapter(store, code$WorkflowState$new(audit_enabled = FALSE))
    on.exit(code$closeDIAStateAdapter(manager), add = TRUE)
    manager$.context$reference_states <- TRUE
    manager$saveState("protein_peptide_filtered", object, list(), "PERF-102 small native fixture")
    app <- new.env(parent = emptyenv())
    app$config <- list(root = temporary, policy = code$workflowExecutionPolicy(TRUE), stem = "fixture")
    app$controllers <- list()
    app$active <- NULL
    app$cgroup <- "seam-only-no-resource-claim"
    # Only the OS launch/accounting seam is substituted. Native science, disk
    # references, leases, head acceptance and required publication remain real.
    code$checkWorkflowApplication <- function(application) invisible(TRUE)
    code$openWorkflowPeakWindow <- function(application) {
        con <- file(tempfile(tmpdir = temporary), "w+")
        writeLines("0", con)
        flush(con)
        con
    }
    code$workflowCgroupAccounting <- function(group) list()
    code$workflowHostAvailable <- function() 8 * 1024^3
    launches <- 0L
    code$workflowLaunchService <- function(...) {
        launches <<- launches + 1L
        args <- list(...)
        directory <- readRDS(args[[4L]])$directory
        saveRDS(list(max_observation_interval = 1, reason = "simulated OS monitor"),
            file.path(directory, "monitor-result.rds"))
        invisible(NULL)
    }
    code$workflowSystemd <- function(...) invisible(NULL)
    code$workflowServiceProperties <- function(...) list(available = TRUE, Result = "success")
    code$workflowServiceEmpty <- function(...) TRUE
    code$workflowMemoryCounter <- function(...) 0
    owner <- code$newWorkflowExecutionController(manager, app)
    on.exit(code$closeWorkflowExecutionController(owner), add = TRUE, after = FALSE)
    input <- list(state = code$getDIAStateReference(manager))
    before <- code$workflowExecutionHead(store)
    rows_before <- store$table_rows
    supplied <- list(peptides_per_sample_cutoff = NULL)
    id <- code$runWorkflowAction("minimum_peptides_per_sample", input, supplied, owner)
    job <- owner$job
    expect_identical(launches, 1L)
    expect_identical(store$table_rows, rows_before)
    expect_identical(code$workflowExecutionHead(store), before)
    supplied$peptides_per_sample_cutoff <- 999L
    request <- readRDS(job$worker_request)
    expect_identical(request$parameters$peptides_per_sample_cutoff, 3L)
    expect_identical(request$token$parameters, digest::digest(request$parameters, algo = "sha256"))
    expect_false(code$hasMutableWorkflowEntityOwner(request))
    request$worker_cap <- 12 * 1024^3
    saveRDS(request, job$worker_request)
    withr::with_envvar(c(MULTISCHOLAR_WORKFLOW_WORKER = "1"),
        code$workflowWorkerMain(c(job$worker_request,
            digest::digest(file = job$worker_request, algo = "sha256"))))
    worker_file <- file.path(job$paths$request_dir, "worker-result.rds")
    message <- readRDS(worker_file)
    expect_identical(message$status, "complete")
    value <- message$result
    expect_identical(code$workflowExecutionHead(store), before)
    native <- get("filterMinNumPeptidesPerSample", asNamespace("MultiScholaR"))(
        object, peptides_per_sample_cutoff = 3L, inclusion_list = character(), core_utilisation = NA)
    actual <- code$materializeDIATable(manager$.context, value$reference)
    expect_identical(actual, native@peptide_data)
    args <- native@args$filterMinNumPeptidesPerSample
    expect_identical(value$summary, args$filter_summary)
    for (name in c("support_table", "removal_ledger")) expect_identical(
        code$materializeDIATable(manager$.context,
            value$checkpoint$support$filterMinNumPeptidesPerSample[[name]]), args[[name]])
    worker_hash <- digest::digest(file = worker_file, algo = "sha256")
    withr::local_options(multischolar.workflow_execution_boundary = function(stage) {
        if (stage == "publication") stop("injected publication failure")
    })
    done <- code$acceptWorkflowActionResult(owner, message)
    expect_identical(done$status, "publication_pending")
    expect_true(done$cleanup$worker_stopped && done$cleanup$monitor_stopped)
    expect_false(done$measurement_valid)
    expect_error(isOpen(job$peak), "invalid connection")
    accepted <- code$getDIAStateReference(manager)
    options(multischolar.workflow_execution_boundary = NULL)
    result <- code$resumeWorkflowAction(owner, id)
    expect_identical(result$state_reference, accepted)
    expect_identical(code$getDIAStateReference(manager), accepted)
    expect_identical(digest::digest(file = worker_file, algo = "sha256"), worker_hash)
    expect_identical(code$acceptWorkflowActionResult(owner, message)$status, "duplicate")
    expect_identical(launches, 1L)
    paths <- list.files(job$paths$request_dir, pattern = "[.](png|pdf)$", recursive = TRUE, full.names = TRUE)
    paths <- paths[!grepl("/figures/", paths)]
    expect_length(paths, 12L)
    expect_true(all(file.size(paths) > 0))
    # The unchanged installed reader can read this accepted physical checkpoint.
    old_reader <- get("readDIAExecutionComponents", asNamespace("MultiScholaR"))
    expect_identical(old_reader(manager, accepted, fields = "peptide_data")@peptide_data, native@peptide_data)
    expect_identical(old_reader(manager, accepted, fields = "design_matrix")@design_matrix, native@design_matrix)
    expect_equal(nrow(native@design_matrix), 7L)
    expect_identical(store$open_resources, 0L)
    fake_sequence <- 0L
    fake_job <- function() {
        x <- new.env(parent = emptyenv())
        fake_sequence <<- fake_sequence + 1L
        x$request <- request
        x$request$token$request_id <- paste0("wfx-", substr(digest::digest(fake_sequence), 1L, 32L))
        x$request$token$expected_head <- code$workflowExecutionHead(store)$entity_id
        x$request$token$generation <- owner$control$generation
        x$paths <- code$workflowExportPaths(store$directory, x$request$token$request_id)
        code$writeWorkflowExportJournal(x$paths, list(schema = "multischolar.workflow_action.v1",
            status = "running", request_id = x$request$token$request_id))
        saveRDS(list(max_observation_interval = 1), file.path(x$paths$request_dir, "monitor-result.rds"))
        x$peak <- code$openWorkflowPeakWindow(app)
        x$admission <- filelock::lock(file.path(temporary, ".heavy-action.lock"), timeout = 0)
        stopifnot(!is.null(x$admission))
        x$accounting <- list()
        x$started <- proc.time()[["elapsed"]]
        x$stage <- "running"
        x$monitor_unit <- "monitor"
        x$worker_unit <- "worker"
        x
    }
    head <- code$workflowExecutionHead(store)
    owner$job <- fake_job()
    stale <- list(token = owner$job$request$token, result = value)
    stale$token$parameters <- strrep("f", 64)
    expect_identical(code$acceptWorkflowActionResult(owner, stale)$status, "stale")
    expect_identical(code$workflowExecutionHead(store), head)
    owner$job <- fake_job()
    code$workflowServiceProperties <- function(unit) if (unit == "monitor")
        list(available = TRUE, SubState = "running") else
        list(available = TRUE, SubState = "failed", Result = "exit-code")
    expect_identical(code$pollWorkflowAction(owner)$status, "worker_failed")
    expect_identical(code$workflowExecutionHead(store), head)
    owner$job <- fake_job()
    generation <- owner$control$generation
    code$cancelWorkflowAction(owner, "reverted")
    expect_gt(owner$control$generation, generation)
    expect_null(owner$job)
    expect_identical(code$workflowExecutionHead(store), head)
    code$closeWorkflowExecutionController(owner)
    expect_null(store$execution_owner)
    expect_error(code$resumeWorkflowAction(owner, id), class = "workflow_execution_stale")
})
