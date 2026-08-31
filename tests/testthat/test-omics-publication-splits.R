test_that("pilot and holdout ownership is source project and seed disjoint", {
    splits <- publicationGovernanceRead("splits-v1.json")
    projects <- publicationGovernanceRead("projects-v1.json")

    expect_silent(publicationGovernanceValidateSplits(splits, projects))
    expect_false(splits$readiness$candidate_access_allowed)
    expect_true(splits$readiness$successor_required_before_candidate)
    expect_true(all(vapply(splits$assignments, \(assignment) {
        !isTRUE(assignment$promotion_eligible)
    }, logical(1))))
})

test_that("split leakage and empty promoted splits fail closed", {
    splits <- publicationGovernanceRead("splits-v1.json")
    projects <- publicationGovernanceRead("projects-v1.json")

    overlap <- publicationGovernanceCopy(splits)
    overlap$assignments[[1L]]$pilot_project_ids <- list("project-1")
    overlap$assignments[[1L]]$holdout_project_ids <- list("project-1")
    expect_error(
        publicationGovernanceValidateSplits(overlap, projects),
        class = "multischolar_publication_governance_error"
    )

    empty_promoted <- publicationGovernanceCopy(splits)
    empty_promoted$assignments[[1L]]$promotion_eligible <- TRUE
    expect_error(
        publicationGovernanceValidateSplits(empty_promoted, projects),
        class = "multischolar_publication_governance_error"
    )
})

test_that("adaptive threshold uses exact bytes and deterministic grid", {
    threshold <- publicationGovernanceRead("threshold-grid-v1.json")

    expect_silent(publicationGovernanceValidateThresholdGrid(threshold))
    expect_identical(
        threshold$size_measure$measure_id,
        "total_uncompressed_input_bytes_v1"
    )
    expect_false(threshold$size_measure$uses_materialized_r_payload)
    expect_false(threshold$size_measure$uses_object_size_multiplier)
    expect_false(threshold$size_measure$uses_fitted_expansion_factor)

    heuristic <- publicationGovernanceCopy(threshold)
    heuristic$size_measure$uses_object_size_multiplier <- TRUE
    expect_error(
        publicationGovernanceValidateThresholdGrid(heuristic),
        class = "multischolar_publication_governance_error"
    )

    non_monotone <- publicationGovernanceCopy(threshold)
    non_monotone$threshold_bytes[[4L]] <- 1L
    expect_error(
        publicationGovernanceValidateThresholdGrid(non_monotone),
        class = "multischolar_publication_governance_error"
    )
})
