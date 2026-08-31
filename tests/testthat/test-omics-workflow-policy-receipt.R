.workflowPolicyRepoPath <- function(...) {
    file.path(
        normalizePath(testthat::test_path("..", ".."), mustWork = TRUE),
        ...
    )
}

.workflowPolicyEnvelopePath <- function() {
    .workflowPolicyRepoPath(
        "inst",
        "extdata",
        "omics-auto-policy-receipts-v2.json"
    )
}

.workflowPolicyCopy <- function(value) {
    unserialize(serialize(value, NULL, version = 3L))
}

.workflowPolicySealReceipt <- function(receipt) {
    receipt$receipt_digest <- NULL
    receipt$receipt_digest <- workflowPolicyObjectDigest(receipt)
    receipt
}

.workflowPolicySealEnvelope <- function(envelope) {
    envelope$envelope_digest <- NULL
    envelope$envelope_digest <- workflowPolicyObjectDigest(envelope)
    envelope
}

test_that("policy receipt v1 remains immutable and v2 envelope is generated", {
    v1 <- .workflowPolicyRepoPath(
        "tests",
        "testdata",
        "omics-performance",
        "policy-receipt-schema-v1.json"
    )
    expect_identical(
        digest::digest(file = v1, algo = "sha256", serialize = FALSE),
        "20ed982ddfbd4bb0720d4c1ca908bba0754d861028b3329e3c6314fa4fa3c386"
    )

    envelope <- workflowPolicyLoadEnvelope(
        .workflowPolicyEnvelopePath(),
        strict = TRUE
    )
    expect_silent(workflowPolicyValidateEnvelope(envelope))
    expect_length(envelope$receipts, 2L)
    expect_identical(
        envelope$consumer_contract_digest,
        workflowPolicyConsumerContractDigest()
    )
    expect_false(envelope$generation_contract$copied_values_allowed)
    expect_false(envelope$generation_contract$production_defaults_changed)
})

test_that("generated legacy receipts exactly replay routing boundaries", {
    envelope <- workflowPolicyLoadEnvelope(
        .workflowPolicyEnvelopePath(),
        strict = TRUE
    )
    descriptors <- artifactDescriptorCatalogueValues(
        artifactWorkflowDescriptorCatalogue()
    )
    for (receipt in envelope$receipts) {
        descriptor <- descriptors[[receipt$capability_id]]
        expect_identical(receipt$threshold_bytes, 33554432L)
        expect_identical(
            receipt$size_measure$measure_id,
            "legacy_r_object_size_x2.v1"
        )
        below <- workflowPolicyResolve(
            envelope,
            receipt$capability_id,
            receipt$size_measure$measure_id,
            receipt$threshold_bytes - 1L,
            platform = "Linux",
            descriptor_digest = descriptor$descriptor_digest
        )
        equal <- workflowPolicyResolve(
            envelope,
            receipt$capability_id,
            receipt$size_measure$measure_id,
            receipt$threshold_bytes,
            platform = "Linux",
            descriptor_digest = descriptor$descriptor_digest
        )
        expect_identical(below$effective_backend, "memory")
        expect_identical(below$effective_rollout, "none")
        expect_identical(equal$effective_backend, "artifact")
        expect_identical(equal$effective_rollout, "evict")
        expect_false(equal$publication_authority)
    }
    expect_identical(
        artifactPayloadAutoPromotionGate()$minimum_projected_source_bytes,
        33554432
    )
})

test_that("exact proposed receipts are injectable but do not alter defaults", {
    installed <- workflowPolicyLoadEnvelope(
        .workflowPolicyEnvelopePath(),
        strict = TRUE
    )
    proposed <- .workflowPolicyCopy(installed)
    receipt <- proposed$receipts[[1L]]
    receipt$receipt_id <- "test.proposed.exact.v2"
    receipt$receipt_kind <- "proposed_pilot"
    receipt$owner_ticket_id <- "OMICS-ART-077"
    receipt$decision <- "proposed_pilot"
    receipt["terminal_reason"] <- list(NULL)
    receipt$size_measure <- list(
        measure_id = "total_uncompressed_input_bytes_v1",
        unit = "byte",
        exact = TRUE,
        available_before_full_parse = TRUE
    )
    proposed$receipts <- list(.workflowPolicySealReceipt(receipt))
    proposed <- .workflowPolicySealEnvelope(proposed)

    expect_silent(workflowPolicyValidateEnvelope(proposed))
    decision <- workflowPolicyResolve(
        proposed,
        receipt$capability_id,
        "total_uncompressed_input_bytes_v1",
        receipt$threshold_bytes,
        platform = "Linux",
        descriptor_digest = receipt$descriptor_digest
    )
    expect_identical(decision$effective_backend, "artifact")
    expect_false(decision$publication_authority)

    reloaded <- workflowPolicyLoadEnvelope(
        .workflowPolicyEnvelopePath(),
        strict = TRUE
    )
    expect_identical(
        reloaded$envelope_digest,
        installed$envelope_digest
    )
})

test_that("policy drift and inadmissible legacy promotion fail closed", {
    envelope <- workflowPolicyLoadEnvelope(
        .workflowPolicyEnvelopePath(),
        strict = TRUE
    )
    receipt <- envelope$receipts[[1L]]

    stale_descriptor <- workflowPolicyResolve(
        envelope,
        receipt$capability_id,
        receipt$size_measure$measure_id,
        receipt$threshold_bytes,
        platform = "Linux",
        descriptor_digest = strrep("0", 64L)
    )
    expect_identical(stale_descriptor$effective_backend, "memory")
    expect_identical(stale_descriptor$reason_code, "policy_descriptor_stale")

    stale_consumer <- .workflowPolicyCopy(envelope)
    stale_consumer$consumer_contract_digest <- strrep("0", 64L)
    stale_consumer <- .workflowPolicySealEnvelope(stale_consumer)
    unavailable <- workflowPolicyLoadEnvelope(injected = stale_consumer)
    expect_s3_class(unavailable, "WorkflowPolicyEnvelopeUnavailable")
    expect_identical(
        workflowPolicyResolve(
            unavailable,
            receipt$capability_id,
            receipt$size_measure$measure_id,
            receipt$threshold_bytes
        )$effective_backend,
        "memory"
    )

    duplicate <- .workflowPolicyCopy(envelope)
    duplicate$receipts[[2L]] <- duplicate$receipts[[1L]]
    duplicate <- .workflowPolicySealEnvelope(duplicate)
    expect_error(
        workflowPolicyValidateEnvelope(duplicate),
        class = "multischolar_invalid_workflow_policy_envelope"
    )

    promoted_legacy <- .workflowPolicyCopy(receipt)
    promoted_legacy$receipt_kind <- "final_installed"
    promoted_legacy$decision <- "final_installed"
    promoted_legacy$claim_scope$publication_authority <- TRUE
    promoted_legacy <- .workflowPolicySealReceipt(promoted_legacy)
    expect_error(
        workflowPolicyValidateReceipt(promoted_legacy),
        class = "multischolar_invalid_workflow_policy_measure"
    )
})

test_that("receipt loading and resolution are artifact-resource inert", {
    source <- readLines(
        .workflowPolicyRepoPath("R", "utils_workflow_policy_receipt.R"),
        warn = FALSE
    )
    prohibited <- c(
        "newArtifactStore(", "newProjectRegistry(", "DBI::", "duckdb::",
        "newArtifactQuerySession(", "ArtifactWorkflowState$new("
    )
    expect_false(any(vapply(prohibited, function(value) {
        any(grepl(value, source, fixed = TRUE))
    }, logical(1))))

    missing <- workflowPolicyLoadEnvelope(path = tempfile("missing-policy-"))
    expect_s3_class(missing, "WorkflowPolicyEnvelopeUnavailable")
    expect_identical(
        workflowPolicyResolve(
            missing,
            "unknown.capability",
            "total_uncompressed_input_bytes_v1",
            1L
        )$effective_backend,
        "memory"
    )
})

test_that("runtime policy values have one generated authority", {
    sources <- unlist(lapply(c(
        file.path("R", "utils_workflow_policy_receipt.R"),
        file.path("R", "utils_artifact_payload_eviction.R"),
        file.path("tools", "profiling", "generate_omics_policy_receipt.R")
    ), function(path) {
        readLines(.workflowPolicyRepoPath(path), warn = FALSE)
    }), use.names = FALSE)
    forbidden <- c(
        "metabolomics.custom.metabolite.standard.v1",
        "lipidomics.lipidsearch.lipid.standard.v1",
        "33554432",
        "32 * 1024"
    )
    expect_false(any(vapply(forbidden, function(value) {
        any(grepl(value, sources, fixed = TRUE))
    }, logical(1))))
})
