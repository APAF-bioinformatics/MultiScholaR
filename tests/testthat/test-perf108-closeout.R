test_that("retained-runtime closeout uses both populations and preserves failure gates", {
    helpers <- new.env(parent = baseenv())
    sys.source(file.path(here::here(), "tools/profiling/workflow_perf108.R"), helpers)
    bindings <- list(protocol_sha256 = "frozen", full = list(primary = list(), wide = list()))
    for (role in c("accepted", "current")) bindings[[role]] <- list(library = role, build = role)
    cases <- helpers$perf108CacheCases(bindings, "closure-filters")
    expect_length(cases, 8L)
    expect_identical(vapply(cases, `[[`, character(1), "role"),
        rep(c("accepted", "current", "current", "accepted"), 2L))
    expect_identical(vapply(cases, `[[`, character(1), "population"),
        rep(c("primary", "wide"), each = 4L))
    expect_true(all(vapply(cases, function(x) is.null(x$structural_options), logical(1))))
    receipt <- list(stage = "closure-filters", cases = list())
    for (case in cases) receipt$cases[[case$id]] <- list(passed = TRUE,
        verifier_completion = list(accepted = TRUE), operation = list(status = "complete",
            measurement_valid = TRUE, lifecycle_passed = TRUE, swap_bytes = 0,
            peak_bytes = 3 * 1024^3, elapsed_seconds = 100))
    expect_identical(helpers$perf108CacheDecision(receipt, cases)$status, "runtime_non_regression_passed")
    for (field in c("elapsed_seconds", "peak_bytes")) {
        changed <- receipt
        changed$cases[["wide-current-2"]]$operation[[field]] <-
            changed$cases[["wide-current-2"]]$operation[[field]] * 1.051
        expect_identical(helpers$perf108CacheDecision(changed, cases)$status, "held")
    }
    receipt$cases[["primary-current-1"]]$verifier_completion$accepted <- FALSE
    expect_identical(helpers$perf108CacheDecision(receipt, cases)$status, "held")
    expect_false(helpers$perf108CacheDecision(receipt, cases)$science_outputs_and_measurements_valid)
})
