test_that("structural confirmation preserves prospective gates and fixed ordering", {
    helpers <- new.env(parent = baseenv())
    sys.source(file.path(here::here(), "tools/profiling/workflow_perf108_structural.R"), helpers)
    bindings <- list(protocol_sha256 = "frozen", screen = list(p14 = list(rows = 960000L)),
        full = list(p14 = list(rows = 24000000L), primary = list(rows = 6225000L), wide = list(rows = 6225000L)))
    for (role in c("accepted", "current", "candidate")) bindings[[role]] <- list(library = role, build = role)
    cases <- helpers$perf108StructuralCases(bindings, "structural-rollup")
    expect_identical(vapply(cases, `[[`, character(1), "role"),
        c("accepted", "current", "candidate", "candidate", "current", "accepted"))
    expect_length(helpers$perf108StructuralCases(bindings, "structural-screen"), 4L)
    ordered <- bindings
    ordered$worker_options <- list(multischolar.perf108.ordered_groups = TRUE)
    ordered_cases <- helpers$perf108StructuralCases(ordered, "ordered-screen")
    expect_identical(ordered_cases[[2L]]$structural_options, ordered$worker_options)
    expect_null(ordered_cases[[1L]]$structural_options)
    for (stage in c("structural-primary", "structural-wide")) expect_true(all(vapply(
        helpers$perf108StructuralCases(bindings, stage), function(x) is.null(x$structural_options), logical(1))))
    receipt <- list(stage = "structural-rollup", cases = list())
    for (case in cases) receipt$cases[[case$id]] <- list(passed = TRUE,
        verifier_completion = list(accepted = TRUE), operation = list(status = "complete",
            measurement_valid = TRUE, lifecycle_passed = TRUE, swap_bytes = 0,
            peak_bytes = if (case$role == "candidate") 6 * 1024^3 else 8 * 1024^3,
            elapsed_seconds = if (case$role == "candidate") 110 else 100))
    decision <- helpers$perf108StructuralDecision(receipt, cases)
    expect_identical(decision$status, "rollup_candidate_supported")
    expect_false(decision$promotion_authorized)
    changed <- receipt
    changed$cases[["p14-candidate-1"]]$operation$elapsed_seconds <- 115.01
    expect_identical(helpers$perf108StructuralDecision(changed, cases)$status, "no_justified_rollup_gain")
    changed <- receipt
    changed$cases[["p14-accepted-2"]]$operation$elapsed_seconds <- 116
    expect_identical(helpers$perf108StructuralDecision(changed, cases)$status, "inconclusive_reference_variability")
    changed <- receipt
    changed$cases[["p14-candidate-1"]]$verifier_completion$accepted <- FALSE
    expect_identical(helpers$perf108StructuralDecision(changed, cases)$status, "invalid_or_failed_member")
})
