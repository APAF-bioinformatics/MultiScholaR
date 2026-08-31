test_that("estimands freeze every supported subject and measured phase", {
    estimands <- publicationGovernanceRead("estimands-v1.json")
    coverage <- publicationGovernanceRead("coverage-v1.json")

    expect_silent(publicationGovernanceValidateEstimands(estimands, coverage))
    expect_identical(
        estimands$primary_import_work_unit$work_unit_id,
        "validated_input_byte"
    )
    expect_identical(
        estimands$primary_import_work_unit$count_source,
        "total_uncompressed_input_bytes_v1"
    )
    expect_true(estimands$primary_import_work_unit$backend_invariant)
})

test_that("estimands reject unknown phases markers and work units", {
    estimands <- publicationGovernanceRead("estimands-v1.json")
    coverage <- publicationGovernanceRead("coverage-v1.json")

    unknown <- publicationGovernanceCopy(estimands)
    unknown$capability_phase_assignments[[1L]]$phase_ids[[1L]] <- "unknown"
    expect_error(
        publicationGovernanceValidateEstimands(unknown, coverage),
        class = "multischolar_publication_governance_error"
    )

    missing_marker <- publicationGovernanceCopy(estimands)
    missing_marker$phase_definitions[[1L]]$start_marker <- ""
    expect_error(
        publicationGovernanceValidateEstimands(missing_marker, coverage),
        class = "multischolar_publication_governance_error"
    )

    changed_fields <- publicationGovernanceCopy(estimands)
    changed_fields$phase_definitions[[1L]]$manual_denominator <- "rows"
    expect_error(
        publicationGovernanceValidateEstimands(changed_fields, coverage),
        class = "multischolar_publication_governance_error"
    )
})

test_that("historical and matched backend claim classes remain separate", {
    estimands <- publicationGovernanceRead("estimands-v1.json")
    definitions <- stats::setNames(
        estimands$phase_definitions,
        vapply(estimands$phase_definitions, `[[`, character(1), "phase_id")
    )

    expect_true("release_effect" %in%
        unlist(definitions$complete_workflow$allowed_effect_classes))
    expect_true("backend_effect" %in%
        unlist(definitions$complete_workflow$allowed_effect_classes))
    expect_false("release_effect" %in%
        unlist(definitions$scientific_stage$allowed_effect_classes))
})
