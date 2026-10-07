test_that("PERF-108 attribution flushes real query phases on completion and failure", {
    root <- Sys.getenv("MULTISCHOLAR_PERF108_ATTRIBUTION_ROOT")
    skip_if(!nzchar(root), "explicit bounded attribution preflight")
    repo <- normalizePath(Sys.getenv("MULTISCHOLAR_PERF108_SOURCE"))
    helpers <- new.env(parent = asNamespace("MultiScholaR"))
    sys.source(file.path(repo, "tools/profiling/workflow_perf106.R"), helpers)
    sys.source(file.path(repo, "tools/profiling/workflow_perf108.R"), helpers)
    fixture <- paste(
        'local({ ns <- asNamespace("MultiScholaR");',
        'con <- fastduckdb::open_query_fast(); class(con) <- "workflow_fast_query";',
        'on.exit(fastduckdb::close_query_fast(con));',
        'phase <- getOption("multischolar.dia_rollup_phase");',
        'query <- "SELECT i%10000 AS k, sum(i) AS value, row_number() OVER (ORDER BY i%10000) AS n FROM range(100000) t(i) GROUP BY k ORDER BY k";',
        'for (cohort in 1:3) phase("ordered_result_stage", function() {',
        'get("workflowQueryExecute", ns)(con, paste0("COPY (", query, ") TO ",',
        'get("workflowQueryQuoteString", ns)(con, file.path(d,paste0("result-",cohort,".parquet"))), " (FORMAT PARQUET)")) });',
        'phase("fixture_allocation", function() { x <- rep.int(as.raw(1L),33554432L); Sys.sleep(0.2); length(x) });',
        'if (readRDS(commandArgs(TRUE)[[1L]])$fail) phase("fixture_failure", function() stop("owned attribution failure"));',
        'stopifnot(fastduckdb::query_data_fast(con,"SELECT count(*) n FROM read_parquet(?)", params=list(file.path(d,"result-3.parquet")))$n==10000); })'
    )
    expression <- helpers$perf108AttributionExpression(
        file.path(repo, "tools/profiling/workflow_perf108.R"), fixture)
    results <- list()
    for (fail in c(FALSE, TRUE)) {
        owned <- file.path(root, if (fail) "failure" else "completion")
        stopifnot(!dir.exists(owned), dir.create(owned, recursive = TRUE))
        request <- file.path(owned, "request.rds")
        saveRDS(list(fail = fail), request)
        status <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
            c("--vanilla", "-e", shQuote(expression), shQuote(request)),
            stdout = file.path(owned, "stdout.log"), stderr = file.path(owned, "stderr.log")))
        expect_identical(status == 0L, !fail)
        start <- readRDS(file.path(owned, "profile-start.rds"))
        finish <- readRDS(file.path(owned, "profile-finished.rds"))
        events <- utils::read.delim(file.path(owned, "phase-events.tsv"))
        expect_identical(finish$error, fail)
        expect_true(all(events$uptime >= start$start & events$uptime <= finish$finished))
        expect_identical(events$event[events$phase == "ordered_result_stage"], rep(c("start", "complete"), 3L))
        if (fail) expect_identical(events$event[events$phase == "fixture_failure"], c("start", "failed"))
        expect_identical(tail(events$event, 1L), "finished")
        profiles <- list.files(owned, "^query-profile-.*[.]json$", full.names = TRUE)
        expect_length(profiles, 3L)
        for (path in profiles) {
            profile <- jsonlite::read_json(path)
            expect_match(profile$query_name, "^COPY")
            expect_true(profile$rows_returned > 0)
        }
        expect_true(file.size(file.path(owned, "worker-stack.out")) > 0)
        allocations <- readLines(file.path(owned, "worker-allocations.out"))
        sizes <- suppressWarnings(as.numeric(sub(" .*", "", allocations)))
        expect_true(any(sizes >= 33554432, na.rm = TRUE))
        results[[basename(owned)]] <- list(exit_status = status, start = start, finish = finish,
            phase_events = nrow(events), copy_profile = profiles[[1L]])
    }
    saveRDS(results, file.path(root, "cases.rds"))
})
