test_that("publication protocol freezes comparators design and authorities", {
    protocol <- publicationGovernanceRead("protocol-v1.json")

    expect_silent(publicationGovernanceValidateProtocol(protocol))
    expect_identical(protocol$status, "frozen_governance_no_results")
    expect_null(protocol$comparators[[3L]]$revision)
    expect_false(protocol$amendment_policy$runtime_behavior_changed_by_this_protocol)
    expect_setequal(
        names(protocol$claim_classes),
        c("release_effect", "backend_effect", "repair_effect")
    )
})

test_that("governance manifest byte-binds every version-one authority", {
    manifest <- publicationGovernanceRead("governance-manifest-v1.json")

    expect_silent(publicationGovernanceValidateManifest(manifest))
    expect_length(manifest$records, 12L)
    expect_true(manifest$immutability$additive_successor_required)
    expect_false(manifest$immutability$mutate_version_in_place)

    mutated <- publicationGovernanceCopy(manifest)
    mutated$records[[1L]]$sha256 <- paste(rep("f", 64L), collapse = "")
    expect_error(
        publicationGovernanceValidateManifest(mutated),
        class = "multischolar_publication_governance_error"
    )
})

test_that("protocol rejects weakened independence optional stopping and drift", {
    protocol <- publicationGovernanceRead("protocol-v1.json")

    too_few_projects <- publicationGovernanceCopy(protocol)
    too_few_projects$independence$minimum_projects_per_cross_project_claim <- 2L
    expect_error(
        publicationGovernanceValidateProtocol(too_few_projects),
        class = "multischolar_publication_governance_error"
    )

    optional <- publicationGovernanceCopy(protocol)
    optional$replication$optional_stopping <- TRUE
    expect_error(
        publicationGovernanceValidateProtocol(optional),
        class = "multischolar_publication_governance_error"
    )

    drift <- publicationGovernanceCopy(protocol)
    drift$authorities$normative_plan$sha256 <- paste(rep("f", 64L), collapse = "")
    expect_error(
        publicationGovernanceValidateProtocol(drift),
        class = "multischolar_publication_governance_error"
    )

    leaked_result <- publicationGovernanceCopy(protocol)
    leaked_result$observed_candidate_result <- list(runtime_ratio = 0.5)
    expect_error(
        publicationGovernanceValidateProtocol(leaked_result),
        class = "multischolar_publication_governance_error"
    )
})

test_that("unassigned role principals block implementation and execution", {
    roles <- publicationGovernanceRead("roles-v1.json")

    expect_silent(publicationGovernanceValidateRoles(roles))
    expect_true(roles$readiness$governance_structure_complete)
    expect_false(roles$readiness$principal_assignments_complete)
    expect_false(roles$readiness$runtime_implementation_authorized)
    expect_false(roles$readiness$campaign_execution_authorized)

    unsafe <- publicationGovernanceCopy(roles)
    unsafe$readiness$runtime_implementation_authorized <- TRUE
    expect_error(
        publicationGovernanceValidateRoles(unsafe),
        class = "multischolar_publication_governance_error"
    )
})

test_that("independent reviewer cannot share implementation or analysis ownership", {
    roles <- publicationGovernanceRead("roles-v1.json")
    assigned <- publicationGovernanceCopy(roles)
    for (index in seq_along(assigned$roles)) {
        assigned$roles[[index]]$principal_id <- paste0("principal-", index)
    }
    assigned$roles[[6L]]$principal_id <- assigned$roles[[3L]]$principal_id
    assigned$readiness$principal_assignments_complete <- TRUE

    expect_error(
        publicationGovernanceValidateRoles(assigned),
        class = "multischolar_publication_governance_error"
    )
})
