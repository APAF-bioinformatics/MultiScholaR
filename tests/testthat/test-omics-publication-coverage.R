test_that("publication coverage exactly mirrors current capability inventory", {
    coverage <- publicationGovernanceRead("coverage-v1.json")

    expect_silent(publicationGovernanceValidateCoverage(coverage))
    expect_identical(coverage$counts$full_workflows, 6L)
    expect_identical(coverage$counts$reader_boundaries, 3L)
    expect_identical(coverage$counts$inventory_capabilities, 11L)
    expect_identical(coverage$counts$inventory_formats, 16L)
})

test_that("coverage rejects omitted duplicated and overstated routes", {
    coverage <- publicationGovernanceRead("coverage-v1.json")

    omitted <- publicationGovernanceCopy(coverage)
    omitted$full_workflows <- omitted$full_workflows[-1L]
    expect_error(
        publicationGovernanceValidateCoverage(omitted),
        class = "multischolar_publication_governance_error"
    )

    duplicate <- publicationGovernanceCopy(coverage)
    duplicate$reader_boundaries[[1L]]$capability_id <-
        duplicate$full_workflows[[1L]]$capability_id
    expect_error(
        publicationGovernanceValidateCoverage(duplicate),
        class = "multischolar_publication_governance_error"
    )

    overstated <- publicationGovernanceCopy(coverage)
    overstated$reader_boundaries[[1L]]$full_workflow_promotion_authority <- TRUE
    expect_error(
        publicationGovernanceValidateCoverage(overstated),
        class = "multischolar_publication_governance_error"
    )
})

test_that("missing independent projects remain explicit and non-promotional", {
    coverage <- publicationGovernanceRead("coverage-v1.json")
    projects <- publicationGovernanceRead("projects-v1.json")

    expect_silent(publicationGovernanceValidateProjects(projects, coverage))
    expect_true(all(vapply(projects$full_workflow_claims, \(claim) {
        !isTRUE(claim$promotion_eligible) &&
            identical(
                claim$cross_project_claim_status,
                "insufficient_independent_projects"
            )
    }, logical(1))))

    false_promotion <- publicationGovernanceCopy(projects)
    false_promotion$full_workflow_claims[[2L]]$promotion_eligible <- TRUE
    expect_error(
        publicationGovernanceValidateProjects(false_promotion, coverage),
        class = "multischolar_publication_governance_error"
    )

    synthetic_substitution <- publicationGovernanceCopy(projects)
    synthetic_substitution$full_workflow_claims[[2L]]$current_real_project_authorities <-
        list(list(source_kind = "public_generated"))
    expect_error(
        publicationGovernanceValidateProjects(synthetic_substitution, coverage),
        class = "multischolar_publication_governance_error"
    )
})
