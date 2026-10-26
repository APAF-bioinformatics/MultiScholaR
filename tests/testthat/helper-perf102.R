# PERF-102 audit-repair fixtures shared by test-perf-003.R and test-perf-102-audit-repairs.R (KIT-034, quarantine Q-05:
# they lived at the top of test-perf-003.R, so the PERF-102 file could not run on its own).
perf102Code <- function(science = FALSE) {
    parent <- if (science) asNamespace("MultiScholaR") else baseenv()
    code <- new.env(parent = parent)
    sys.source(testthat::test_path("..", "..", "R/utils_workflow_execution.R"), code)
    # KIT-024 added the in-process dispatch refusal (in utils_workflow_operations.R); the sourced copy gets a pass-through.
    code$assertWorkflowOperationDispatchable <- function(operation) invisible(NULL)
    # KIT-007 moved the operation facts into the registry (utils_workflow_operations.R); the sourced copy of the
    # execution file gets a minimal spec, as it already gets a stubbed validator in the tests.
    code$workflowActionSpec <- utils::getFromNamespace("workflowActionSpec", "MultiScholaR")
    # KIT-045: the execution file calls the workflow names; the DIA names stay for the tests' direct calls.
    if (science) for (name in c("newWorkflowEntityStore", "newDIAStateAdapter", "newWorkflowStateOwner",
        "closeDIAStateAdapter", "closeWorkflowStateOwner", "WorkflowState", "getDIAStateReference", "getWorkflowStateReference",
        "materializeDIATable", "materializeWorkflowTable", "hasMutableWorkflowEntityOwner",
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

