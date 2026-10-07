# Metadata-only runner regressions; never launch a scientific job or service.
perf106RunnerCode <- function() {
    code <- new.env(parent = baseenv())
    sys.source(testthat::test_path("..", "..", "tools/profiling/workflow_perf106.R"), code)
    code
}

test_that("verification requires bound successful process and terminal completion", {
    code <- perf106RunnerCode()
    root <- withr::local_tempdir()
    owned <- file.path(root, "operation")
    run <- file.path(root, "browser-smoke.Abc123")
    dir.create(owned)
    dir.create(run)
    unit <- "multischolar-dia003-browser-smoke.Abc123.service"
    invocation <- strrep("a", 32)
    operation <- file.path(owned, "operation-receipt.rds")
    saveRDS(list(status = "complete"), operation)
    value <- list(status = "complete", run_dir = run,
        operation_sha256 = digest::digest(file = operation, algo = "sha256"))
    saveRDS(value, file.path(owned, "verification.rds"))
    writeLines(jsonlite::toJSON(list(unit = unit, invocation_id = invocation), auto_unbox = TRUE),
        file.path(owned, "verification-launch.log"))
    terminal <- list(status = "service_tree_interval_only", service_succeeded = TRUE,
        service_result = "success", unit = unit, invocation_id = invocation)
    write_terminal <- function(x) jsonlite::write_json(x, file.path(run, "terminal-interval.json"),
        auto_unbox = TRUE)
    write_terminal(terminal)
    writeLines(c("unit=multischolar-dia003-browser-smoke.Abc123", "interval_exit=0", "stop_exit=0"),
        file.path(run, "controller-status.txt"))
    read <- function(status = 0L) code$readPerf106Verification(owned, root, status)
    expect_true(read()$accepted)
    for (status in list(1L, 74L, 124L, NA_integer_, NULL, integer()))
        expect_false(read(status)$accepted)
    for (field in c("status", "service_succeeded", "service_result", "unit", "invocation_id")) {
        bad <- terminal
        bad[[field]] <- if (field == "service_succeeded") FALSE else "wrong"
        write_terminal(bad)
        expect_false(read()$accepted)
    }
    write_terminal(terminal)
    writeLines("reason=job_early_stop", file.path(run, "guard-stop.txt"))
    expect_false(read()$accepted)
    unlink(file.path(run, "guard-stop.txt"))
    saveRDS(list(status = "different"), operation)
    expect_false(read()$accepted)
    saveRDS(list(status = "complete"), operation)
    writeLines(c("unit=multischolar-dia003-browser-smoke.Abc123", "interval_exit=0", "stop_exit=1"),
        file.path(run, "controller-status.txt"))
    expect_false(read()$accepted)
    unlink(file.path(run, "terminal-interval.json"))
    expect_false(read()$accepted)
})

test_that("member observations reset before an immediate guard can stop the loop", {
    lines <- readLines(testthat::test_path("..", "..", "tools/profiling/workflow_perf106_repair.R"))
    initialization <- grep('properties <- list(available = FALSE, reason = "service not yet observed")',
        lines, fixed = TRUE)
    loop <- grep("for (case in cases) {", lines, fixed = TRUE)[[1L]]
    stopifnot(length(initialization) == 1L)
    expect_identical(initialization, loop + 1L)
    receipt_line <- grep("receipt$cases[[case$id]] <- list(passed", lines, fixed = TRUE)
    expression <- parse(text = paste(lines[receipt_line:(receipt_line + 1L)], collapse = "\n"))
    state <- new.env(parent = baseenv())
    state$receipt <- list(cases = list())
    state$passed <- FALSE
    state$result <- state$tests <- NULL
    for (id in c("first", "second")) {
        state$case <- list(id = id)
        state$properties <- list(available = TRUE, unit = "previous-member")
        state$outer_guard <- "previous-reason"
        eval(parse(text = lines[initialization:(initialization + 2L)]), state)
        expect_null(state$outer_guard)
        state$outer_guard <- "memory_guard"
        eval(expression, state)
        expect_false(state$receipt$cases[[id]]$properties$available)
        expect_null(state$receipt$cases[[id]]$properties$unit)
    }
    status_line <- grep('saveRDS(list(status = if (exists("outer_guard")', lines, fixed = TRUE)
    stopifnot(length(status_line) == 1L)
    state$owned <- withr::local_tempdir()
    state$member <- list(id = "guard-fixture")
    state$error <- simpleError("terminal observation failed")
    eval(parse(text = paste(lines[status_line:(status_line + 2L)], collapse = "\n")), state)
    expect_identical(readRDS(file.path(state$owned, "operation-receipt.rds"))$status, "memory_guard")
})
