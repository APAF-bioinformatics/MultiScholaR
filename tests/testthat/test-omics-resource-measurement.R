.omicsResourceRepoPath <- function(...) {
    file.path(normalizePath(testthat::test_path("..", ".."), mustWork = TRUE), ...)
}

source(
    .omicsResourceRepoPath("tools", "profiling", "omics_resource_measurement.R"),
    local = TRUE
)

test_that("resource measurement records a bounded process tree and retention point", {
    skip_on_cran()
    skip_if_not_installed("processx")
    skip_if_not_installed("ps")

    run_dir <- tempfile("omics-resource-run-")
    dir.create(run_dir, recursive = TRUE)
    expression <- sprintf(
        paste(
            "x <- numeric(250000)",
            "writeLines('ready', %s)",
            "Sys.sleep(0.12)",
            sep = "; "
        ),
        dQuote(file.path(run_dir, "retention-ready"))
    )
    measured <- omicsResourceMeasureSubprocess(
        command = file.path(R.home("bin"), "Rscript"),
        command_args = c("--vanilla", "-e", expression),
        run_dir = run_dir,
        execution = list(sampling_interval_ms = 10L),
        categories = list(list(category = "harness", pattern = ".*")),
        repo_root = .omicsResourceRepoPath(),
        timeout_ms = 5000L,
        maximum_retained_samples = 3L,
        require_retention_marker = TRUE
    )

    expect_identical(measured$status, "passed")
    expect_gt(measured$metrics$peak_tree_rss_bytes, 0)
    expect_gt(measured$metrics$retained_tree_rss_bytes, 0)
    expect_gt(measured$metrics$sample_count, 3L)
    expect_lte(length(measured$samples), 3L)
    expect_true(measured$samples_truncated)
    expect_identical(measured$retained_sample_count, length(measured$samples))
    expect_true(measured$retention_marker_observed)
})

test_that("resource measurement fails closed on timeout and process failure", {
    skip_on_cran()
    skip_if_not_installed("processx")
    skip_if_not_installed("ps")

    timeout_dir <- tempfile("omics-resource-timeout-")
    timeout <- omicsResourceMeasureSubprocess(
        command = file.path(R.home("bin"), "Rscript"),
        command_args = c("--vanilla", "-e", "Sys.sleep(2)"),
        run_dir = timeout_dir,
        execution = list(sampling_interval_ms = 10L),
        categories = list(list(category = "harness", pattern = ".*")),
        repo_root = .omicsResourceRepoPath(),
        timeout_ms = 50L,
        maximum_retained_samples = 5L
    )
    expect_identical(timeout$status, "failed")
    expect_identical(timeout$failure, "timeout")

    crash_dir <- tempfile("omics-resource-crash-")
    crash <- omicsResourceMeasureSubprocess(
        command = file.path(R.home("bin"), "Rscript"),
        command_args = c("--vanilla", "-e", "quit(status = 7L)"),
        run_dir = crash_dir,
        execution = list(sampling_interval_ms = 10L),
        categories = list(list(category = "harness", pattern = ".*")),
        repo_root = .omicsResourceRepoPath(),
        timeout_ms = 5000L,
        maximum_retained_samples = 5L
    )
    expect_identical(crash$status, "failed")
    expect_identical(crash$exit_status, 7L)
})

test_that("directory metrics classify files once and retain unassigned totals", {
    run_dir <- tempfile("omics-resource-disk-")
    dir.create(file.path(run_dir, "committed"), recursive = TRUE)
    dir.create(file.path(run_dir, "staging"), recursive = TRUE)
    writeLines("committed", file.path(run_dir, "committed", "a.txt"))
    writeLines("staged", file.path(run_dir, "staging", "b.txt"))

    metrics <- omicsResourceDirMetrics(
        run_dir,
        list(
            list(category = "committed", pattern = "(^|/)committed(/|$)"),
            list(category = "staging", pattern = "(^|/)staging(/|$)"),
            list(category = "other", pattern = ".*")
        )
    )

    expect_identical(metrics$file_count, 2L)
    expect_gt(metrics$bytes$committed, 0)
    expect_gt(metrics$bytes$staging, 0)
    expect_identical(metrics$counts$other, 0L)
    expect_equal(
        metrics$total_bytes,
        metrics$bytes$committed + metrics$bytes$staging
    )
})
