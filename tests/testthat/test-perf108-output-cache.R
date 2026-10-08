test_that("staged-file cache ownership preserves ordered pages and failure isolation", {
    skip_if_not_installed("fastduckdb")
    role <- "candidate"  # CORE-015: the staged reader always clears the native cache on entry and exit.
    root <- tempfile("cache-fixture-")
    dir.create(root)
    on.exit(unlink(root, recursive = TRUE))
    path <- file.path(root, "staged.parquet")
    n <- 20005L
    original <- data.frame(.dia_output_order = seq_len(n),
        value = as.double(seq_len(n)), label = rep(c("a", "\u03b2", NA_character_), length.out = n))
    arrow::write_parquet(original[n:1L, ], path, chunk_size = 10240L)
    clear <- fastduckdb::clear_parquet_cache
    clears <- 0L
    testthat::local_mocked_bindings(clear_parquet_cache = function() {
        clears <<- clears + 1L
        clear()
    }, .package = "fastduckdb")
    parts <- list()
    consume <- function(part) {
        expect_lte(nrow(part), 10000L)
        expect_lte(as.numeric(object.size(part)), 67108864)
        parts[[length(parts) + 1L]] <<- part
    }
    expect_identical(streamDIARelationalResultFile(path, consume, rows = n), 3L)
    actual <- do.call(rbind, parts)
    rownames(actual) <- NULL
    expect_identical(actual, original[c("value", "label")])
    expect_identical(clears, if (role == "candidate") 2L else 3L)

    # Different projections/ranges must not become cached query results.
    expect_identical(fastduckdb::read_parquet_fast(path, columns = "label", n_threads = 8L,
        lazy_strings = FALSE, filter = list(.dia_output_order = c(n, n)))$label, original$label[[n]])
    replacement <- file.path(root, "replacement.parquet")
    updated <- original
    updated$value <- -updated$value
    arrow::write_parquet(updated, replacement, chunk_size = 10240L)
    stopifnot(file.rename(replacement, path))
    expect_identical(fastduckdb::read_parquet_fast(path, columns = "value", n_threads = 8L,
        lazy_strings = FALSE, filter = list(.dia_output_order = c(1, 2)))$value, c(-1, -2))

    clears <- 0L
    expect_error(streamDIARelationalResultFile(path,
        function(part) stop("owned cache consumer failure"), rows = n), "owned cache consumer failure")
    expect_identical(clears, if (role == "candidate") 2L else 1L)
    parts <- list()
    streamDIARelationalResultFile(path, consume, rows = n)
    actual <- do.call(rbind, parts)
    rownames(actual) <- NULL
    expect_identical(actual, updated[c("value", "label")])
    writeBin(as.raw(c(1L, 2L, 3L)), path)
    clears <- 0L
    expect_error(streamDIARelationalResultFile(path, consume, rows = n))
    expect_identical(clears, if (role == "candidate") 2L else 1L)
    arrow::write_parquet(updated, path, chunk_size = 10240L)
    expect_identical(fastduckdb::read_parquet_fast(path, columns = "value", n_threads = 8L,
        filter = list(.dia_output_order = c(n, n)))$value, -as.double(n))
    clears <- 0L
    expect_identical(streamDIARelationalResultFile(path,
        function(part) stop("empty input must not be consumed"), rows = 0), 0L)
    expect_identical(clears, if (role == "candidate") 2L else 0L)
    clear()
    expect_length(list.files(root, pattern = "replacement"), 0L)
})

test_that("cache qualification refuses invalid controls and preserves the selection rules", {
    skip_if(!Sys.getenv("MULTISCHOLAR_PERF108_CACHE_ROLE") %in% c("current", "candidate"),
        "explicit installed cache qualification")
    helpers <- new.env(parent = baseenv())
    sys.source(file.path(here::here(), "tools/profiling/workflow_perf108.R"), helpers)
    cases <- lapply(c("accepted", "current", "candidate"), function(role)
        lapply(1:2, \(i) list(id = paste("p14", role, i, sep = "-"))))
    cases <- unlist(cases, recursive = FALSE)
    make <- function(seconds, peak = 4 * 1024^3) list(passed = TRUE,
        verifier_completion = list(accepted = TRUE), operation = list(status = "complete",
            measurement_valid = TRUE, lifecycle_passed = TRUE, swap_bytes = 0,
            peak_bytes = peak, elapsed_seconds = seconds))
    receipt <- list(stage = "cache-rollup", cases = list())
    for (i in 1:2) {
        receipt$cases[[paste0("p14-accepted-", i)]] <- make(100)
        receipt$cases[[paste0("p14-current-", i)]] <- make(100)
        receipt$cases[[paste0("p14-candidate-", i)]] <- make(85)
    }
    decision <- helpers$perf108CacheDecision(receipt, cases)
    expect_identical(decision$selected, "candidate")
    expect_true(decision$meaningful_cache_benefit)
    for (i in 1:2) receipt$cases[[paste0("p14-candidate-", i)]] <- make(95)
    expect_identical(helpers$perf108CacheDecision(receipt, cases)$selected, "current")
    receipt$cases[["p14-current-2"]] <- make(110)
    expect_identical(helpers$perf108CacheDecision(receipt, cases)$selected, "candidate")
    receipt$cases[["p14-accepted-2"]]$operation$measurement_valid <- FALSE
    expect_null(helpers$perf108CacheDecision(receipt, cases)$selected)
    receipt$cases[["p14-accepted-2"]] <- make(100)
    receipt$cases[["p14-candidate-2"]] <- make(85, peak = 5 * 1024^3)
    expect_null(helpers$perf108CacheDecision(receipt, cases)$selected)
})
