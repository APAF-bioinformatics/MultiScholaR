test_that("passive diagnostic boundaries never invoke their heap collector", {
    api <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling",
        "workflow_checkpoint_gui_diagnostic.R"), api)
    api$gc <- function(...) stop("diagnostic collection must stay disabled")
    root <- withr::local_tempdir()
    context <- api$startCheckpointGuiDiagnostic(root, collect_heap = FALSE)
    withr::defer(api$stopCheckpointGuiDiagnostic(context))
    expect_null(api$checkpointGuiCompletedHeap(context, "completed"))
    fixture <- new.env(parent = baseenv())
    fixture$work <- function(fail) {
        if (fail) stop("fixture failure")
        sum(numeric(100000L))
    }
    api$installCheckpointGuiTraces(context, fixture, "work")
    expect_identical(fixture$work(FALSE), 0)
    expect_error(fixture$work(TRUE), "fixture failure")
    records <- lapply(readLines(file.path(root, "application-phases.jsonl")), jsonlite::fromJSON)
    stages <- vapply(records, `[[`, character(1), "stage")
    expect_identical(sum(stages == "exit:work"), 2L)
    expect_false(any(startsWith(stages, "heap:")))
    expect_false(records[[1L]]$detail$completed_action_gc)
})

test_that("the actual diagnostic browser wrapper starts in passive mode", {
    skip("quarantined Q-14 (docs/qa/workflow-execution/kit/quarantined-tests.md)")
    api <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling",
        "workflow_checkpoint_dia_driver.R"), api)
    library <- Sys.getenv("MULTISCHOLAR_DIA008_LIBRARY")
    expect_true(dir.exists(library))
    api$checkpointInstalledApp <- function(library_path, run_dir, capacity_treatment,
        session_format = "legacy") {
        expect_identical(session_format, "legacy")
        launch <- function() {
            loadNamespace("MultiScholaR", lib.loc = fixed_library)
            shiny::shinyApp(shiny::fluidPage(shiny::actionButton("go", "Run"),
                shiny::textOutput("result")), function(input, output, session) {
                shiny::observeEvent(input$shutdown, shiny::stopApp(), ignoreInit = TRUE)
                output$result <- shiny::renderText({
                    shiny::req(input$go > 0)
                    MultiScholaR:::workflowExportBoundary("passive_fixture_complete")
                    as.character(sum(numeric(100000L)))
                })
            })
        }
        environment(launch) <- list2env(list(fixed_library = library_path), parent = baseenv())
        launch
    }
    root <- withr::local_tempdir()
    # Preserve the actual stages even when an assertion or shutdown fails.
    withr::defer({
        receipt <- Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT")
        if (nzchar(receipt)) {
            for (name in c("application-phases.jsonl", "application-console.log")) {
                path <- file.path(root, name)
                if (file.exists(path)) file.copy(path,
                    file.path(dirname(receipt), paste0("passive-browser-", name)))
            }
        }
    })
    helper <- normalizePath(test_path("..", "..", "tools", "profiling",
        "workflow_checkpoint_gui_diagnostic.R"))
    result <- api$newCheckpointBrowser(library, root, profile_application = TRUE,
        capacity_treatment = "qs2_duckdb", diagnostic_source = helper)
    driver <- result$driver
    withr::defer(driver$stop())
    driver$click("go")
    driver$wait_for_value(output = "result")
    expect_identical(driver$get_value(output = "result"), "0")
    # AppDriver$stop may terminate the R process after its interrupt timeout.
    # Exercise ordinary Shiny teardown explicitly before transport cleanup.
    driver$run_js("Shiny.setInputValue('shutdown', 1, {priority: 'event'});")
    deadline <- unname(proc.time()[["elapsed"]]) + 5
    repeat {
        phases <- readLines(file.path(root, "application-phases.jsonl"), warn = FALSE)
        if (any(grepl('"stage":"diagnostic:stop"', phases, fixed = TRUE)) ||
            unname(proc.time()[["elapsed"]]) >= deadline) break
        Sys.sleep(0.05)
    }
    driver$stop()
    records <- lapply(readLines(file.path(root, "application-phases.jsonl")), jsonlite::fromJSON)
    stages <- vapply(records, `[[`, character(1), "stage")
    expect_false(records[[1L]]$detail$completed_action_gc)
    expect_true("enter:workflowExportBoundary" %in% stages, info = paste(stages, collapse = ", "))
    expect_true("diagnostic:stop" %in% stages, info = paste(stages, collapse = ", "))
    expect_false(any(startsWith(stages, "heap:")))
})


test_that("allocator release precedes the recorded memory observation", {
    api <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "workflow_checkpoint_gui_diagnostic.R"), api)
    called <- FALSE
    context <- new.env(parent = emptyenv())
    context$allocator_snapshot <- function(trim) {
        called <<- TRUE
        c(trim_return = as.integer(trim))
    }
    api$checkpointGuiDiagnosticEvent <- function(context, stage, detail) {
        expect_true(called)
        expect_identical(detail$glibc$trim_return, 1L)
    }
    api$checkpointAllocatorBoundary(context, "ordered", trim = TRUE)
})


test_that("recorded release comparison preserves exact corresponding output bytes", {
    library <- Sys.getenv("MULTISCHOLAR_DIA008_LIBRARY")
    skip_if(!nzchar(library), "explicit installed evidence check")
    proof_root <- dirname(dirname(library))
    member <- file.path(proof_root, "browser-smoke.0dj9B1")
    reference <- file.path(proof_root, "browser-smoke.6EqMOx", "anonymous-memory-composite")
    api <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "workflow_checkpoint_gui_diagnostic.R"), api)
    before <- digest::digest(file = file.path(member, "allocator-release-comparison.json"), algo = "sha256")
    verified <- api$verifyCheckpointAllocatorOutputs(member, reference)
    expect_identical(verified$status, "passed_exact_corresponding_outputs")
    expect_false(verified$iteration_rds_bytes_identical)
    expect_identical(digest::digest(file = file.path(member, "allocator-release-comparison.json"),
        algo = "sha256"), before)
    verified$source_member <- member
    verified$original_receipt_sha256 <- before
    verified$original_member_status <- "failed_single_render_RDS_oracle"
    jsonlite::write_json(verified, file.path(dirname(Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT")),
        "allocator-output-verification.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
})
