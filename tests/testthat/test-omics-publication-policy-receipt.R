publicationReceiptExample <- function() {
    digest_value <- paste(rep("a", 64L), collapse = "")
    list(
        schema = "multischolar.omics_auto_policy_receipt",
        schema_version = "1.0.0",
        receipt_id = "receipt.example.confirmed.v1",
        receipt_kind = "confirmatory_decision",
        owner_ticket_id = "OMICS-ART-078",
        capability_id = "proteomics.diann.peptide.dia.v1",
        capability_version = "1.0.0",
        descriptor_digest = digest_value,
        platforms = list("Linux"),
        size_measure = list(
            measure_id = "total_uncompressed_input_bytes_v1",
            unit = "byte",
            exact = TRUE,
            available_before_full_parse = TRUE
        ),
        threshold_bytes = 67108864L,
        below_threshold = list(backend = "memory", rollout = "none"),
        at_or_above_threshold = list(backend = "artifact", rollout = "evict"),
        decision = "confirmed_promote",
        terminal_reason = NULL,
        claim_scope = list(
            support_tier = "scientifically_supported",
            project_scope = "cross_project",
            host_scope = "validated_linux_hosts",
            publication_authority = TRUE
        ),
        authority_bindings = list(list(
            binding_id = "protocol",
            path_or_uri = "tests/testdata/omics-performance/protocol-v1.json",
            sha256 = digest_value
        )),
        evidence_bindings = list(list(
            binding_id = "raw-evidence",
            path_or_uri = "private://raw-evidence",
            sha256 = digest_value
        )),
        project_ids = list("project-1", "project-2", "project-3"),
        pilot_split_id = "pilot-v1",
        holdout_split_id = "holdout-v1",
        threshold_grid_id = "threshold-grid-v1",
        estimand_ids = list("estimand-1"),
        candidate_revision = paste(rep("b", 40L), collapse = ""),
        role_bindings = list(list(
            role_id = "independent_reviewer",
            principal_id = "reviewer-1",
            handoff_digest = digest_value
        )),
        gate_results = list(list(
            gate_id = "scientific_parity",
            status = "passed",
            evidence_digest = digest_value
        )),
        rollback = list(
            new_work_backend = "memory",
            existing_artifact_mode = "read_only",
            implicit_memory_migration = FALSE,
            delete_artifacts = FALSE
        ),
        digest_contract = list(
            algorithm = "sha256",
            canonicalization = "canonical_json_sorted_keys_utf8_v1",
            excluded_field = "receipt_digest"
        ),
        receipt_digest = digest_value
    )
}

test_that("one receipt schema governs proposed confirmatory and final policy", {
    schema <- publicationGovernanceRead("policy-receipt-schema-v1.json")

    expect_silent(publicationGovernanceValidateReceiptSchema(schema))
    expect_silent(publicationGovernanceValidateReceipt(
        publicationReceiptExample(),
        schema
    ))
    expect_false(schema$additionalProperties)
    expect_identical(
        schema$properties$size_measure$properties$measure_id$const,
        "total_uncompressed_input_bytes_v1"
    )
})

test_that("receipt schema and instances reject heuristic or incomplete promotion", {
    schema <- publicationGovernanceRead("policy-receipt-schema-v1.json")

    weakened <- publicationGovernanceCopy(schema)
    weakened$additionalProperties <- TRUE
    expect_error(
        publicationGovernanceValidateReceiptSchema(weakened),
        class = "multischolar_publication_governance_error"
    )

    heuristic <- publicationReceiptExample()
    heuristic$size_measure$measure_id <- "legacy_r_object_size_x2.v1"
    expect_error(
        publicationGovernanceValidateReceipt(heuristic, schema),
        class = "multischolar_publication_governance_error"
    )

    too_few_projects <- publicationReceiptExample()
    too_few_projects$project_ids <- list("project-1", "project-2")
    expect_error(
        publicationGovernanceValidateReceipt(too_few_projects, schema),
        class = "multischolar_publication_governance_error"
    )

    unknown_field <- publicationReceiptExample()
    unknown_field$manual_threshold_override <- TRUE
    expect_error(
        publicationGovernanceValidateReceipt(unknown_field, schema),
        class = "multischolar_publication_governance_error"
    )
})

test_that("no-promotion receipt must resolve memory with a terminal reason", {
    schema <- publicationGovernanceRead("policy-receipt-schema-v1.json")
    receipt <- publicationReceiptExample()
    receipt$decision <- "confirmed_no_promotion"
    receipt["threshold_bytes"] <- list(NULL)
    receipt$at_or_above_threshold <- list(backend = "memory", rollout = "none")
    receipt$claim_scope$publication_authority <- FALSE
    receipt$terminal_reason <- "insufficient_independent_projects"

    expect_silent(publicationGovernanceValidateReceipt(receipt, schema))

    missing_reason <- publicationGovernanceCopy(receipt)
    missing_reason$terminal_reason <- NULL
    expect_error(
        publicationGovernanceValidateReceipt(missing_reason, schema),
        class = "multischolar_publication_governance_error"
    )
})
