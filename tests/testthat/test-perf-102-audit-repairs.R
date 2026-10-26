# Small audit-repair fixtures. No benchmark, GUI or managed worker is launched.
test_that("PERF member identities bind the frozen candidate and distinct executions", {
    code <- perf102Code()
    base <- perf102Evidence(code)
    evaluate <- function(x) code$evaluateWorkflowQualification(x, base$bindings)
    expect_identical(evaluate(base)$status, "benefit")
    for (field in c("input", "parameters", "runtime", "build", "boundary", "output_contract", "arm", "run_id")) {
        x <- base
        x$populations$wide[[2L]]$candidate[[field]] <- if (field == "arm") "selection" else strrep("f", 64)
        expect_identical(evaluate(x)$status, "held", info = field)
        expect_null(evaluate(x)$ratios)
    }
    x <- base
    x$populations$wide[[2L]]$candidate$receipt_id <- x$populations$primary[[1L]]$candidate$receipt_id
    expect_identical(evaluate(x)$status, "held")
    x <- base
    x$populations$primary[[2L]] <- x$populations$primary[[1L]]
    expect_identical(evaluate(x)$status, "held")
    x <- base
    x$execution_plan$finalist <- "selection"
    x$bindings$execution_plan <- digest::digest(x$execution_plan, algo = "sha256")
    expect_identical(evaluate(x)$status, "held")
    # Even a newly supplied independent plan must obey the fixed pilot structure.
    expect_identical(code$evaluateWorkflowQualification(x, x$bindings)$status, "held")
    x <- base
    x$execution_plan$screen_decision <- NULL
    x$bindings$execution_plan <- digest::digest(x$execution_plan, algo = "sha256")
    expect_identical(code$evaluateWorkflowQualification(x, x$bindings)$status, "held")
    x <- base
    for (rep in 1:2) for (role in c("control", "candidate")) {
        slot <- paste("wide", role, rep, sep = ".")
        x$execution_plan$members[[slot]]$input <- base$populations$primary[[1L]]$control$input
        x$populations$wide[[rep]][[role]]$input <- base$populations$primary[[1L]]$control$input
    }
    x$bindings$execution_plan <- digest::digest(x$execution_plan, algo = "sha256")
    expect_identical(code$evaluateWorkflowQualification(x, x$bindings)$status, "held")
    screen <- perf102Evidence(code, "screen_feasibility")
    expect_identical(code$evaluateWorkflowQualification(screen, screen$bindings)$status, "screen_finalist")
    changed <- screen
    changed$arms$selection[[2L]]$build <- strrep("f", 64)
    expect_identical(code$evaluateWorkflowQualification(changed, screen$bindings)$status, "held")
})

test_that("PERF collection validation rejects missing extra and malformed names", {
    code <- perf102Code()
    corrupt <- function(value, kind) {
        if (kind == "extra") value <- c(value, list(extra = value[[1L]]))
        if (kind == "NA") value <- c(value, stats::setNames(list(value[[1L]]), NA_character_))
        if (kind == "blank") value <- c(value, stats::setNames(list(value[[1L]]), ""))
        if (kind == "duplicate") names(value)[2L] <- names(value)[1L]
        if (kind == "missing") value <- value[-1L]
        if (kind == "unnamed") names(value) <- NULL
        value
    }
    for (goal in c("benefit", "screen_feasibility")) for (kind in
        c("extra", "NA", "blank", "duplicate", "missing", "unnamed")) {
        base <- perf102Evidence(code, goal)
        x <- base
        field <- if (goal == "benefit") "populations" else "arms"
        x[[field]] <- corrupt(x[[field]], kind)
        expect_identical(code$evaluateWorkflowQualification(x, base$bindings)$status, "held",
            info = paste(goal, kind))
        x <- base
        x$execution_plan$members <- corrupt(x$execution_plan$members, kind)
        x$bindings$execution_plan <- digest::digest(x$execution_plan, algo = "sha256")
        expect_identical(code$evaluateWorkflowQualification(x, x$bindings)$status, "held")
        if (goal == "benefit") {
            x <- base
            x$populations$primary[[1L]] <- corrupt(x$populations$primary[[1L]], kind)
            expect_identical(code$evaluateWorkflowQualification(x, base$bindings)$status, "held")
        }
    }
})

test_that("action refusal releases actual counters and locks before ownership transfer", {
    code <- perf102Code()
    root <- withr::local_tempdir()
    peak <- NULL
    reference <- list(entity_id = strrep("a", 64))
    owner <- new.env(parent = emptyenv())
    class(owner) <- "WorkflowExecutionController"
    owner$closed <- FALSE
    owner$control <- list(generation = 1L, controller_id = "test-controller")
    owner$manager <- list(.context = list(store = list(directory = root,
        identity = list(incarnation = "test-incarnation"))))
    owner$application <- new.env(parent = emptyenv())
    owner$application$config <- list(root = root, policy = list(startup_headroom = 1), stem = "fixture")
    code$openWorkflowPeakWindow <- function(app) {
        peak <<- file(file.path(root, "peak"), "w+")
        peak
    }
    code$checkWorkflowApplication <- code$assertDIAExportReady <- function(...) invisible(NULL)
    code$validateWorkflowActionParameters <- function(...) invisible(NULL)
    code$isSavedWorkflowReference <- function(...) TRUE
    code$workflowHostAvailable <- function() 2
    code$workflowCgroupAccounting <- function(...) list()
    code$validateSavedWorkflowReference <- function(...) invisible(NULL)
    code$getDIAStateReference <- code$getWorkflowStateReference <- function(...) reference # KIT-045: both names
    code$workflowExecutionHead <- function(...) reference
    code$workflowExecutionBuildIdentity <- function(...) strrep("b", 64)
    code$workflowExportPaths <- function(directory, id) {
        path <- file.path(directory, id)
        dir.create(path)
        list(request_dir = path)
    }
    code$writeWorkflowExportJournal <- code$writeWorkflowExecutionControl <- function(...) invisible(NULL)
    code$workflowLaunchService <- function(...) invisible(NULL)
    failure <- ""
    code$validateWorkflowActionInput <- function(...) {
        if (failure == "roles") stop("role refusal")
    }
    code$resolveWorkflowActionParameters <- function(...) {
        if (failure == "parameters") stop("parameter refusal")
        list()
    }
    code$writeWorkflowExecutionMessage <- function(path, value) {
        if (failure == "serialization") stop("serialization refusal")
        saveRDS(value, path)
    }
    lock_available <- function() callr::r(function(path) {
        lock <- filelock::lock(path, timeout = 0)
        if (is.null(lock)) return(FALSE)
        filelock::unlock(lock)
        TRUE
    }, args = list(file.path(root, ".heavy-action.lock")))
    for (stage in c("roles", "parameters", "serialization")) {
        failure <- stage
        expect_error(code$runWorkflowAction("minimum_peptides_per_sample",
            list(state = reference), list(), owner), "refusal")
        expect_error(isOpen(peak), "invalid connection")
        expect_true(lock_available())
        expect_null(owner$job)
        expect_null(owner$application$active)
    }
    failure <- ""
    id <- code$runWorkflowAction("minimum_peptides_per_sample", list(state = reference), list(), owner)
    expect_identical(owner$job$request$token$request_id, id)
    expect_true(isOpen(peak))
    expect_false(lock_available())
    # The real post-transfer teardown is exercised in the registered disk fixture.
    close(owner$job$peak)
    filelock::unlock(owner$job$admission)
    expect_true(lock_available())
})
