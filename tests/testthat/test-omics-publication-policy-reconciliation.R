test_that("current policy reconciliation preserves and supersedes exact authorities", {
    record <- publicationGovernanceRead("current-policy-reconciliation-v1.json")

    expect_silent(publicationGovernanceValidateReconciliation(record))
    expect_identical(record$owner_ticket_id, "OMICS-ART-062")
    expect_identical(
        record$reconciled_current_decision$decision,
        "adaptive_scale_promotion_implemented_engineering_only"
    )
    expect_false(
        record$reconciled_current_decision$publication_decision_supported
    )
    expect_true(record$supersession$preserve_predecessor_bytes)
    expect_setequal(
        vapply(
            record$supersession$later_evidence,
            `[[`,
            character(1),
            "ticket_id"
        ),
        c("OMICS-ART-051", "OMICS-ART-052")
    )
})

test_that("policy reconciliation rejects predecessor and runtime drift", {
    record <- publicationGovernanceRead("current-policy-reconciliation-v1.json")

    bad_hash <- publicationGovernanceCopy(record)
    bad_hash$predecessors[[1L]]$sha256 <- paste(rep("f", 64L), collapse = "")
    expect_error(
        publicationGovernanceValidateReconciliation(bad_hash),
        class = "multischolar_publication_governance_error"
    )

    bad_runtime <- publicationGovernanceCopy(record)
    bad_runtime$current_runtime_source$sha256 <- paste(rep("a", 64L), collapse = "")
    expect_error(
        publicationGovernanceValidateReconciliation(bad_runtime),
        class = "multischolar_publication_governance_error"
    )

    bad_capability <- publicationGovernanceCopy(record)
    bad_capability$implemented_decisions[[1L]]$capability_id <-
        "proteomics.diann.peptide.dia.v1"
    expect_error(
        publicationGovernanceValidateReconciliation(bad_capability),
        class = "multischolar_publication_governance_error"
    )
})

test_that("engineering policy evidence cannot become publication authority", {
    record <- publicationGovernanceRead("current-policy-reconciliation-v1.json")

    promoted_static <- publicationGovernanceCopy(record)
    promoted_static$implemented_decisions[[1L]]$static_certification$auto_eligible <-
        TRUE
    expect_error(
        publicationGovernanceValidateReconciliation(promoted_static),
        class = "multischolar_publication_governance_error"
    )

    publication_claim <- publicationGovernanceCopy(record)
    publication_claim$implemented_decisions[[1L]]$engineering_evidence$publication_authority <-
        TRUE
    expect_error(
        publicationGovernanceValidateReconciliation(publication_claim),
        class = "multischolar_publication_governance_error"
    )

    changed_threshold <- publicationGovernanceCopy(record)
    changed_threshold$implemented_decisions[[1L]]$dynamic_runtime_decision$threshold_bytes <-
        1
    expect_error(
        publicationGovernanceValidateReconciliation(changed_threshold),
        class = "multischolar_publication_governance_error"
    )
})
