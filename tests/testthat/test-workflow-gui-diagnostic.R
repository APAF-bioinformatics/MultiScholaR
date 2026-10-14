test_that("GUI diagnostics persist memory boundaries and failure without allocation tracing", {
    source <- normalizePath(test_path("..", "..", "tools", "profiling",
        "workflow_checkpoint_gui_diagnostic.R"))
    root <- withr::local_tempdir()
    result <- callr::r(function(source, root) {
        api <- new.env(parent = baseenv())
        sys.source(source, api)
        run <- function(fail) {
            directory <- file.path(root, if (fail) "failure" else "success")
            dir.create(directory)
            original_stderr <- sink.number(type = "message")
            context <- api$startCheckpointGuiDiagnostic(directory)
            on.exit(api$stopCheckpointGuiDiagnostic(context), add = TRUE)
            stopifnot(identical(sink.number(type = "message"), original_stderr))
            fixture <- new.env(parent = baseenv())
            fixture$work <- function(fail) {
                x <- numeric(2000000)
                started <- proc.time()[[1L]]
                while (proc.time()[[1L]] - started < 0.08) sum(x)
                if (fail) stop("intentional diagnostic fixture failure")
                sum(x)
            }
            api$installCheckpointGuiTraces(context, fixture, "work")
            fixture$work(fail)
        }
        value <- run(FALSE)
        error <- tryCatch(run(TRUE), error = conditionMessage)
        list(value = value, error = error,
            observer_cleared = is.null(getOption("multischolar.gui_diagnostic_observer")))
    }, args = list(source, root))
    expect_identical(result$value, 0)
    expect_match(result$error, "intentional diagnostic fixture failure")
    expect_true(result$observer_cleared)
    for (case in c("success", "failure")) {
        directory <- file.path(root, case)
        phases <- lapply(readLines(file.path(directory, "application-phases.jsonl")), jsonlite::fromJSON)
        stages <- vapply(phases, `[[`, character(1), "stage")
        expect_true(all(c("diagnostic:start", "enter:work", "exit:work", "diagnostic:stop") %in% stages))
        expect_gt(phases[[2L]]$memory$VmRSS, 0)
        expect_false(file.exists(file.path(directory, "application.Rprofmem")))
        expect_false(file.exists(file.path(directory, "application.Rprof")))
        expect_lt(file.size(file.path(directory, "application-phases.jsonl")), 65536)
    }
})

test_that("completed-action collection distinguishes retained and released live heap", {
    api <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "workflow_checkpoint_gui_diagnostic.R"), api)
    context <- api$startCheckpointGuiDiagnostic(withr::local_tempdir(), collect_heap = TRUE)
    withr::defer(api$stopCheckpointGuiDiagnostic(context))
    owner <- new.env(parent = emptyenv())
    before <- api$checkpointGuiCompletedHeap(context, "before")
    owner$payload <- numeric(2000000L)
    held <- api$checkpointGuiCompletedHeap(context, "held")
    expect_gt(held - before, 12 * 1024^2)
    frame <- list2env(list(normData = owner), parent = emptyenv())
    metadata <- api$checkpointGuiFrameOwners(frame)
    expect_identical(metadata$normData$payload$class, "numeric")
    expect_false("value" %in% names(metadata$normData$payload))
    rm("payload", envir = owner)
    released <- api$checkpointGuiCompletedHeap(context, "released")
    expect_gt(held - released, 12 * 1024^2)
    expect_identical(ls(context)[vapply(as.list(context), \(x) identical(x, owner), logical(1))], character())
})

test_that("small protein worker probe records payload and namespace cost without changing results", {
    skip_if(!identical(Sys.getenv("MULTISCHOLAR_DIA025_WORKER_PROBE"), "true"),
        "explicit bounded attribution batch only")
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    local <- removeProteinsWithOnlyOneReplicate(object, NA, "group")
    pool <- multidplyr::new_cluster(4L)
    withr::defer(stopDIAWorkerCluster(pool))
    observations <- list()
    ns <- asNamespace("multidplyr")
    previous <- getOption("multischolar.worker_payload_probe")
    withr::defer(options(multischolar.worker_payload_probe = previous))
    options(multischolar.worker_payload_probe = function(frame) {
        code <- get("code", frame, inherits = FALSE)
        observations[[length(observations) + 1L]] <<- list(
            serialized_code_bytes = length(serialize(code, NULL, version = 3)))
    })
    base::trace("cluster_call", exit = quote(getOption("multischolar.worker_payload_probe")(environment())),
        where = ns, print = FALSE)
    withr::defer(base::untrace("cluster_call", where = ns))
    parallel <- removeProteinsWithOnlyOneReplicate(object, pool, "group")
    workers <- multidplyr::cluster_call(pool, list(pid = Sys.getpid(),
        namespaces = sort(loadedNamespaces()), status = readLines("/proc/self/status")))
    expect_identical(parallel@protein_quant_table, local@protein_quant_table)
    expect_identical(parallel@design_matrix, local@design_matrix)
    receipt <- Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT")
    jsonlite::write_json(list(table_bytes = as.numeric(object.size(object@protein_quant_table)),
        observations = observations, workers = workers, performance_acceptance = FALSE),
        file.path(dirname(receipt), "worker-attribution.json"), auto_unbox = TRUE, pretty = TRUE)
})

test_that("a real browser survives a long blocking action with diagnostic polling", {
    source <- normalizePath(test_path("..", "..", "tools", "profiling",
        "workflow_checkpoint_gui_diagnostic.R"))
    root <- withr::local_tempdir()
    launch <- function() {
        api <- new.env(parent = baseenv())
        sys.source(fixed_source, api)
        context <- api$startCheckpointGuiDiagnostic(fixed_root)
        shiny::onStop(function() api$stopCheckpointGuiDiagnostic(context), session = NULL)
        shiny::shinyApp(shiny::fluidPage(shiny::actionButton("go", "Run"),
            shiny::textOutput("result")), function(input, output, session) {
            output$result <- shiny::renderText({
                shiny::req(input$go > 0)
                Sys.sleep(190)
                value <- sum(numeric(2000000))
                api$checkpointGuiDiagnosticEvent(context, "fixture:action_complete")
                as.character(value)
            })
        })
    }
    environment(launch) <- list2env(list(fixed_source = source, fixed_root = root), parent = baseenv())
    driver <- shinytest2::AppDriver$new(app_dir = launch, name = "dia-diagnostic-startup",
        load_timeout = 20000, timeout = 230000, check_names = FALSE)
    withr::defer(driver$stop())
    harness <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "workflow_checkpoint_dia_driver.R"), harness)
    withr::local_options(multischolar.dia_memory_audit_directory = root)
    driver$get_js("document.getElementById('go').click();true")
    harness$checkpointWaitForJs(driver, "document.getElementById('result').textContent==='0'",
        timeout = 230000)
    expect_identical(driver$get_value(output = "result"), "0")
    driver$stop()
    phases <- lapply(readLines(file.path(root, "application-phases.jsonl")), jsonlite::fromJSON)
    stages <- vapply(phases, `[[`, character(1), "stage")
    expect_true(all(c("diagnostic:start", "fixture:action_complete", "diagnostic:stop") %in% stages))
})

test_that("the installed journey trace list resolves and can be removed", {
    skip("quarantined Q-20 (docs/qa/workflow-execution/kit/quarantined-tests.md)")
    api <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "workflow_checkpoint_gui_diagnostic.R"), api)
    root <- withr::local_tempdir()
    context <- api$startCheckpointGuiDiagnostic(root)
    withr::defer(api$stopCheckpointGuiDiagnostic(context))
    api$installCheckpointGuiTraces(context, asNamespace("MultiScholaR"))
    expect_true("resolveDIAStateValue" %in% context$traced)
    expect_true("buildProtNormCorrelationPlot" %in% context$traced)
    expect_true("runTestsContrasts" %in% context$traced)
    expect_true(all(c("hydrateSavedWorkflowEntity", "runDIAOwnedSessionExport",
        "workflowExportBoundary", "stageWorkflowLegacyExport",
        "publishWorkflowLegacyExport") %in% context$traced))
    expect_true(all(c(".peptideQcCanonicalise", "runProteinIqRollupApplyStep",
        "generateProtNormCompositeFromFiles", "runProtSummaryPublicationCopy") %in% context$traced))
    metadata <- api$checkpointGuiOwnerSizes(list(table = data.frame(x = 1:3), env = new.env()))
    expect_identical(metadata$table$shape, c(3L, 1L))
    expect_true(is.na(metadata$table$bytes))
    expect_true(is.na(metadata$env$bytes))
})

test_that("export memory traces identify actual hydration and serialization including failure", {
    skip("quarantined Q-20 (docs/qa/workflow-execution/kit/quarantined-tests.md)")
    api <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "workflow_checkpoint_gui_diagnostic.R"), api)
    root <- withr::local_tempdir()
    context <- api$startCheckpointGuiDiagnostic(root)
    withr::defer(api$stopCheckpointGuiDiagnostic(context))
    namespace <- asNamespace("MultiScholaR")
    api$installCheckpointGuiTraces(context, namespace)
    store <- newWorkflowEntityStore(withr::local_tempdir(), "trace-export", "DIA")
    manager <- newDIAStateAdapter(store, WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    manager$saveState("raw", object, list(), "fixture")
    states <- withDIAStateScope(manager, function()
        namespace$materializeDIAStates(manager$.context, share_immutable = TRUE))
    paths <- workflowExportPaths(withr::local_tempdir(), "trace-success")
    namespace$stageWorkflowLegacyExport(list(r6_complete_states = states), paths,
        "trace-success", "20260922_000000")
    namespace$publishWorkflowLegacyExport(paths)
    failed_paths <- workflowExportPaths(withr::local_tempdir(), "trace-failure")
    expect_error(namespace$stageWorkflowLegacyExport(list(r6_complete_states = states), failed_paths,
        "trace-failure", "20260922_000001", serializeFn = function(...) stop("fixture failure")),
        "fixture failure")
    records <- lapply(readLines(file.path(root, "application-phases.jsonl")), jsonlite::fromJSON)
    stages <- vapply(records, `[[`, character(1), "stage")
    expect_true(all(c("enter:materializeDIAStates", "exit:materializeDIAStates",
        "enter:hydrateSavedWorkflowEntity", "exit:hydrateSavedWorkflowEntity",
        "enter:stageWorkflowLegacyExport", "exit:stageWorkflowLegacyExport",
        "enter:publishWorkflowLegacyExport", "exit:publishWorkflowLegacyExport") %in% stages))
    hydration <- records[[which(stages == "enter:hydrateSavedWorkflowEntity")[[1L]]]]
    expect_type(hydration$detail$label, "character")
    expect_gt(hydration$detail$payload_bytes, 0)
    expect_false(any(c("object", "value", "data") %in% names(hydration$detail)))
    boundaries <- vapply(records[stages == "enter:workflowExportBoundary"],
        \(record) record$detail$boundary, character(1))
    expect_identical(sum(boundaries == "before_serialization"), 2L)
    expect_identical(sum(boundaries == "after_serialization"), 1L)
    expect_identical(sum(stages == "exit:stageWorkflowLegacyExport"), 2L)
    expect_false("manager" %in% ls(context))
    expect_lt(file.size(file.path(root, "application-phases.jsonl")), 65536)
})

test_that("memory diagnostics cap their own records and report lost coverage", {
    api <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "workflow_checkpoint_gui_diagnostic.R"), api)
    root <- withr::local_tempdir()
    context <- api$startCheckpointGuiDiagnostic(root)
    withr::defer(api$stopCheckpointGuiDiagnostic(context))
    expect_identical(context$max_log_bytes, 64 * 1024^2)
    context$max_log_bytes <- 1
    api$checkpointGuiDiagnosticEvent(context, "fixture:over_limit")
    path <- file.path(root, "application-phases.jsonl")
    size <- file.size(path)
    api$checkpointGuiDiagnosticEvent(context, "fixture:must_not_grow")
    expect_identical(file.size(path), size)
    records <- lapply(readLines(path), jsonlite::fromJSON)
    expect_identical(utils::tail(records, 1L)[[1L]]$stage, "diagnostic:recording_limit")
    expect_false(utils::tail(records, 1L)[[1L]]$detail$coverage_complete)
})

test_that("the one-hour allowance is restricted to explicit memory diagnostics", {
    driver <- new.env(parent = baseenv())
    runner <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "workflow_checkpoint_dia_driver.R"), driver)
    sys.source(test_path("..", "..", "tools", "profiling", "run_workflow_checkpoint_proof.R"), runner)
    withr::local_options(multischolar.dia_memory_audit_seconds = NULL)
    expect_identical(driver$checkpointWaitBudget(elapsed_seconds = 899), 1000L)
    expect_error(driver$checkpointWaitBudget(elapsed_seconds = 900), "exhausted")
    options(multischolar.dia_memory_audit_seconds = 3600)
    expect_identical(driver$checkpointWaitBudget(elapsed_seconds = 3599), 1000L)
    expect_identical(driver$checkpointWaitBudget(control = TRUE, elapsed_seconds = 1000), 30000L)
    expect_error(driver$checkpointWaitBudget(elapsed_seconds = 3600), "exhausted")
    path <- file.path(withr::local_tempdir(), "request.nul")
    fields <- c("dia-controller-request-v1", "multischolar-dia003-browser-smoke.ABC123.service",
        "browser-baseline", "fixture-library", "17179869184", "0", "3600",
        "Rscript", "--vanilla", "runner.R", "--internal-browser-baseline", "fixture-library",
        "fixture-output", "report.tsv", "proteins.fasta", "design.tsv", "capacity:qs2_duckdb:diagnostic")
    writeBin(fields, path)
    expect_identical(runner$readCheckpointControllerRequest(path)$runtime_max_seconds, 3600)
    fields[[length(fields)]] <- "capacity:qs2_duckdb"
    writeBin(fields, path)
    expect_error(runner$readCheckpointControllerRequest(path), "invalid controller request")
    fields[[7L]] <- "900"
    writeBin(fields, path)
    expect_identical(runner$readCheckpointControllerRequest(path)$runtime_max_seconds, 900)
})

test_that("GUI request boundaries record labels without retaining input values", {
    driver <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "workflow_checkpoint_dia_driver.R"), driver)
    root <- withr::local_tempdir()
    withr::local_options(multischolar.dia_memory_audit_directory = root)
    browser <- list(wait_for_js = function(...) NULL, get_js = function(script, ...) {
        if (startsWith(script, "JSON.stringify")) '{"ok":true,"connected":true}' else NULL
    })
    driver$checkpointClickAndWaitText(browser, "load_session", "fixture", "loaded")
    records <- lapply(readLines(file.path(root, "driver-phases.jsonl")), jsonlite::fromJSON)
    expect_identical(vapply(records, `[[`, character(1), "stage"),
        c("enter:action:load_session", "exit:action:load_session"))
    expect_true(all(vapply(records, \(record) identical(names(record), c("time", "pid", "stage")), logical(1))))
})

test_that("diagnostic and ordinary waits reject dead sessions", {
    api <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "workflow_checkpoint_dia_driver.R"), api)
    withr::local_options(multischolar.dia_memory_audit_directory = tempdir())
    browser <- list(get_js = function(...) NULL)
    expect_error(api$checkpointWaitForJs(browser, "false", timeout = 100), "no response")
    browser$get_js <- function(...) '{"ok":true,"connected":false}'
    expect_error(api$checkpointWaitForJs(browser, "true", timeout = 100), "disconnected")
    options(multischolar.dia_memory_audit_directory = NULL)
    browser$wait_for_js <- function(script, timeout) list(script, timeout)
    result <- api$checkpointWaitForJs(browser, "fixture", timeout = 123)
    expect_match(result[[1L]], "Shiny session disconnected", fixed = TRUE)
    expect_match(result[[1L]], 'eval)("fixture")', fixed = TRUE)
    expect_identical(result[[2L]], 123)
})

test_that("process attribution samples real PID counters separately from the cgroup peak", {
    shell <- readLines(test_path("..", "..", "tools", "profiling", "run_workflow_checkpoint_proof.sh"))
    preamble <- shell[seq_len(which(startsWith(shell, "proof_script_dir=")) - 1L)]
    root <- withr::local_tempdir()
    group <- file.path(root, "group")
    dir.create(group)
    writeLines(as.character(Sys.getpid()), file.path(group, "cgroup.procs"))
    script <- file.path(root, "sample.sh")
    output <- file.path(root, "samples.tsv")
    writeLines(c(preamble, 'checkpointProcessSample "$1" fixture "$2"'), script)
    expect_identical(system2("bash", c(shQuote(script), shQuote(group), shQuote(output))), 0L)
    sample <- utils::read.delim(output, header = FALSE)
    expect_equal(sample[[2L]], Sys.getpid())
    expect_gt(sample[[5L]], 0)
    expect_gt(sample[[9L]], 0)
})


test_that("diagnostic SQL counts dispatch once and expose bounded Arrow pool totals", {
    api <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "workflow_checkpoint_gui_diagnostic.R"), api)
    con <- DBI::dbConnect(duckdb::duckdb())
    withr::defer(DBI::dbDisconnect(con, shutdown = TRUE))
    root <- withr::local_tempdir()
    context <- api$startCheckpointGuiDiagnostic(root)
    withr::defer(api$stopCheckpointGuiDiagnostic(context))
    fixture <- new.env(parent = baseenv())
    fixture$work <- function() NULL
    api$installCheckpointGuiTraces(context, fixture, "work")
    before <- context$sql_queries
    DBI::dbGetQuery(con, "SELECT 1")
    DBI::dbExecute(con, "CREATE TEMP TABLE diagnostics (id INTEGER)")
    expect_identical(context$sql_queries - before, 2)
    arrow::default_memory_pool()
    api$checkpointGuiDiagnosticEvent(context, "pool:observed")
    record <- jsonlite::fromJSON(tail(readLines(file.path(root, "application-phases.jsonl")), 1L))
    expect_gte(record$arrow_pool$allocated_bytes, 0)
    expect_gte(record$arrow_pool$peak_bytes, record$arrow_pool$allocated_bytes)
    expect_type(record$arrow_pool$backend, "character")
})


test_that("retained QC figure replay attributes live heap and allocator residency", {
    skip_if(!identical(Sys.getenv("MULTISCHOLAR_DIA025_WORKER_PROBE"), "true"),
        "explicit guarded attribution only")
    library_path <- Sys.getenv("MULTISCHOLAR_DIA008_LIBRARY")
    proof_root <- dirname(dirname(library_path))
    member <- file.path(proof_root, "browser-smoke.FKZX69")
    source <- file.path(member, "dia_baseline_smoke", "results", "proteomics",
        "publication_graphs", "filtering_qc", "20260925_131212", ".workflow-export",
        "qvalue-a277baf727685fa6432f382a889f20fdbb42cbcd", "display.rds")
    expect_true(file.exists(source))
    receipt_root <- dirname(Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT"))
    diagnostic_source <- normalizePath(test_path("..", "..", "tools", "profiling",
        "workflow_checkpoint_gui_diagnostic.R"))
    result <- callr::r(function(source, output, diagnostic_source) {
        api <- new.env(parent = baseenv())
        sys.source(diagnostic_source, api)
        Rcpp::cppFunction(code = '
Rcpp::NumericVector allocatorSnapshot(bool trim) {
    int released = trim ? malloc_trim(0) : -1;
    struct mallinfo2 info = mallinfo2();
    return Rcpp::NumericVector::create(
        Rcpp::Named("arena_bytes") = (double) info.arena,
        Rcpp::Named("in_use_arena_bytes") = (double) info.uordblks,
        Rcpp::Named("free_arena_bytes") = (double) info.fordblks,
        Rcpp::Named("mmap_bytes") = (double) info.hblkhd,
        Rcpp::Named("trim_return") = released);
}', includes = "#include <malloc.h>", env = environment())
        directory <- file.path(output, "anonymous-memory-figure")
        dir.create(directory)
        context <- api$startCheckpointGuiDiagnostic(directory, collect_heap = TRUE)
        on.exit(api$stopCheckpointGuiDiagnostic(context), add = TRUE)
        observations <- list()
        observe <- function(label, collect = FALSE, trim = FALSE) {
            live <- if (collect) api$checkpointGuiCompletedHeap(context, label) else NULL
            native <- as.list(allocatorSnapshot(trim))
            status <- readLines("/proc/self/status", warn = FALSE)
            memory <- lapply(c("VmRSS", "RssAnon"), function(field) {
                value <- grep(paste0("^", field, ":"), status, value = TRUE)
                as.numeric(strsplit(trimws(sub("^[^:]+:", "", value)), " +")[[1L]][[1L]]) * 1024
            })
            names(memory) <- c("rss_bytes", "anonymous_bytes")
            row <- c(list(label = label, live_R_bytes = live), memory, native)
            observations[[length(observations) + 1L]] <<- row
            api$checkpointGuiDiagnosticEvent(context, "allocator:observation", row)
        }
        loadNamespace("ggplot2")
        loadNamespace("ragg")
        before_hash <- digest::digest(file = source, algo = "sha256")
        display <- readRDS(source)
        observe("grid_loaded", collect = TRUE)
        observe("baseline_trim", collect = TRUE, trim = TRUE)
        for (format in c("png", "pdf")) {
            ggplot2::ggsave(file.path(directory, paste0("retained-grid.", format)),
                plot = display$grid, width = 15, height = 18, dpi = 300)
            observe(paste0("after_", format))
            observe(paste0("after_", format, "_gc"), collect = TRUE)
        }
        display <- NULL
        observe("grid_released_gc", collect = TRUE)
        observe("after_allocator_trim", collect = TRUE, trim = TRUE)
        result <- list(source = source, source_sha256 = before_hash,
            unchanged_source = identical(before_hash, digest::digest(file = source, algo = "sha256")),
            observations = observations, production_change = FALSE, timing_eligible = FALSE,
            limitation = "glibc arena statistics are not a census of other allocators; saved-grid replay is not the full GUI")
        jsonlite::write_json(result, file.path(directory, "attribution.json"),
            auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
        result
    }, args = list(source, receipt_root, diagnostic_source), libpath = .libPaths())
    expect_true(result$unchanged_source)
    expect_length(result$observations, 8L)
    expect_true(all(vapply(result$observations, function(x)
        x$arena_bytes >= x$in_use_arena_bytes && x$free_arena_bytes >= 0, logical(1))))
})


test_that("real allocator boundaries survive failure without retaining inputs", {
    api <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "workflow_checkpoint_gui_diagnostic.R"), api)
    root <- withr::local_tempdir()
    context <- api$startCheckpointGuiDiagnostic(root)
    withr::defer(api$stopCheckpointGuiDiagnostic(context))
    api$enableCheckpointAllocatorProbe(context, arena_details = TRUE)
    api$checkpointAllocatorBoundary(context, "before", collect = TRUE)
    expect_error(local({
        value <- numeric(2000000L)
        api$checkpointAllocatorBoundary(context, "held", collect = TRUE)
        stop("intentional allocator fixture failure")
    }), "intentional allocator fixture failure")
    api$checkpointAllocatorBoundary(context, "released", collect = TRUE)
    api$checkpointAllocatorBoundary(context, "trimmed", collect = TRUE, trim = TRUE)
    records <- lapply(readLines(file.path(root, "application-phases.jsonl")), jsonlite::fromJSON)
    records <- Filter(function(x) identical(x$stage, "allocator:boundary"), records)
    expect_length(records, 4L)
    expect_gt(records[[2L]]$detail$live_R_bytes - records[[1L]]$detail$live_R_bytes, 12000000)
    expect_gt(records[[2L]]$detail$live_R_bytes - records[[3L]]$detail$live_R_bytes, 12000000)
    expect_true(records[[4L]]$detail$trimmed)
    expect_gte(records[[4L]]$detail$glibc$free_arena_bytes, 0)
    expect_gt(records[[4L]]$detail$arenas$count, 0)
    expect_gt(records[[4L]]$detail$mappings$main_heap$rss_bytes, 0)
    expect_gte(records[[4L]]$detail$glibc$main_top_free_bytes, 0)
})


test_that("retained composite compares completed-boundary release and immediate reuse", {
    skip_if(!identical(Sys.getenv("MULTISCHOLAR_DIA025_WORKER_PROBE"), "true"),
        "explicit guarded attribution only")
    library_path <- Sys.getenv("MULTISCHOLAR_DIA008_LIBRARY")
    proof_root <- dirname(dirname(library_path))
    member <- file.path(proof_root, "browser-smoke.FKZX69")
    source <- list.files(file.path(member, "dia_baseline_smoke", "results"),
        pattern = "^pre_norm_pca[.]png$", recursive = TRUE, full.names = TRUE)
    expect_length(source, 1L)
    receipt_root <- dirname(Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT"))
    diagnostic_source <- normalizePath(test_path("..", "..", "tools", "profiling",
        "workflow_checkpoint_gui_diagnostic.R"))
    reference <- file.path(proof_root, "browser-smoke.6EqMOx", "anonymous-memory-composite")
    worker <- function(source, output, diagnostic_source, reference, trim_between) {
        api <- new.env(parent = baseenv())
        sys.source(diagnostic_source, api)
        ns <- loadNamespace("MultiScholaR")
        loadNamespace("patchwork")
        loadNamespace("png")
        loadNamespace("ragg")
        arm <- if (trim_between) "trim" else "control"
        directory <- file.path(output, paste0("allocator-release-cost-", arm))
        dir.create(directory)
        context <- api$startCheckpointGuiDiagnostic(directory)
        on.exit(api$stopCheckpointGuiDiagnostic(context), add = TRUE)
        api$enableCheckpointAllocatorProbe(context)
        inputs <- ns$resolveProtNormCompositeFigureInputs("automatic", dirname(source))
        paths <- inputs$plotFiles[!is.na(inputs$plotFiles)]
        hashes <- vapply(paths, function(path) digest::digest(file = path, algo = "sha256"), character(1))
        ggplot2::set_last_plot(NULL)
        api$checkpointAllocatorBoundary(context, "baseline", collect = TRUE, trim = TRUE)
        render_times <- vector("list", 2L)
        output_hashes <- vector("list", 2L)
        for (iteration in seq_len(2L)) {
            target <- file.path(directory, paste0("render-", iteration))
            dir.create(target)
            started <- proc.time()
            ns$generateProtNormCompositeQcFigure("automatic", target, "proteomics",
                resolveInputsFn = function(...) inputs)
            render_times[[iteration]] <- as.list(unclass(proc.time() - started))
            stopifnot(is.null(ggplot2::get_last_plot()))
            outputs <- file.path(target, paste0("proteomics_composite_QC_figure.", c("png", "pdf", "rds")))
            stopifnot(all(file.exists(outputs)), all(file.size(outputs) > 0))
            output_hashes[[iteration]] <- vapply(outputs[c(1L, 3L)], function(path)
                digest::digest(file = path, algo = "sha256"), character(1))
            api$checkpointAllocatorBoundary(context, paste0("after_render_", iteration), collect = TRUE)
            if (iteration == 1L) {
                started <- proc.time()
                native <- context$allocator_snapshot(trim_between)
                release_time <- as.list(unclass(proc.time() - started))
                api$checkpointAllocatorBoundary(context, "after_between_action_operation")
                release_return <- unname(native[["trim_return"]])
                native <- NULL
            }
        }
        expected <- vapply(file.path(reference, paste0("proteomics_composite_QC_figure.",
            c("png", "rds"))), function(path) digest::digest(file = path, algo = "sha256"), character(1))
        unchanged_outputs <- identical(unname(output_hashes[[1L]]), unname(expected)) &&
            all(vapply(output_hashes, function(value)
                identical(unname(value[[1L]]), unname(expected[[1L]])), logical(1)))
        result <- list(arm = arm, render_times = render_times,
            between_action_time = release_time, release_return = release_return,
            unchanged_outputs = unchanged_outputs,
            source_hashes = as.list(hashes), unchanged_sources = identical(hashes,
                vapply(paths, function(path) digest::digest(file = path, algo = "sha256"), character(1))),
            production_change = FALSE, performance_acceptance = FALSE,
            limitation = "one sequential diagnostic pair; second render measures immediate same-workload reuse, not the next GUI scientific action")
        jsonlite::write_json(result, file.path(directory, "comparison.json"),
            auto_unbox = TRUE, pretty = TRUE, digits = NA)
        result
    }
    results <- lapply(c(FALSE, TRUE), function(trim_between) callr::r(worker,
        args = list(source, receipt_root, diagnostic_source, reference, trim_between),
        libpath = .libPaths()))
    expect_true(all(vapply(results, `[[`, logical(1), "unchanged_sources")))
    expect_true(all(vapply(results, `[[`, logical(1), "unchanged_outputs")))
    expect_identical(vapply(results, `[[`, character(1), "arm"), c("control", "trim"))
    api <- new.env(parent = baseenv())
    sys.source(diagnostic_source, api)
    verification <- api$verifyCheckpointAllocatorOutputs(receipt_root, reference)
    jsonlite::write_json(verification, file.path(receipt_root, "allocator-output-verification.json"),
        auto_unbox = TRUE, pretty = TRUE, digits = NA)
    jsonlite::write_json(results, file.path(receipt_root, "allocator-release-comparison.json"),
        auto_unbox = TRUE, pretty = TRUE, digits = NA)
})
