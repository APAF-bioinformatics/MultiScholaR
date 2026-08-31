.preIngressRepoPath <- function(...) {
    file.path(
        normalizePath(testthat::test_path("..", ".."), mustWork = TRUE),
        ...
    )
}

.preIngressCopy <- function(value) {
    unserialize(serialize(value, NULL, version = 3L))
}

.preIngressSealReceipt <- function(receipt) {
    receipt$receipt_digest <- NULL
    receipt$receipt_digest <- workflowPolicyObjectDigest(receipt)
    receipt
}

.preIngressSealEnvelope <- function(envelope) {
    envelope$envelope_digest <- NULL
    envelope$envelope_digest <- workflowPolicyObjectDigest(envelope)
    envelope
}

.preIngressExactEnvelope <- function(capability_id, descriptor_digest) {
    path <- .preIngressRepoPath(
        "inst",
        "extdata",
        "omics-auto-policy-receipts-v2.json"
    )
    envelope <- workflowPolicyLoadEnvelope(path, strict = TRUE)
    receipt <- .preIngressCopy(envelope$receipts[[1L]])
    receipt$receipt_id <- paste0("test.preingress.", capability_id, ".v2")
    receipt$receipt_kind <- "proposed_pilot"
    receipt$owner_ticket_id <- "OMICS-ART-077"
    receipt$capability_id <- capability_id
    receipt$descriptor_digest <- descriptor_digest
    receipt$decision <- "proposed_pilot"
    receipt["terminal_reason"] <- list(NULL)
    receipt$size_measure <- list(
        measure_id = "total_uncompressed_input_bytes_v1",
        unit = "byte",
        exact = TRUE,
        available_before_full_parse = TRUE
    )
    receipt$threshold_bytes <- 1L
    envelope$receipts <- list(.preIngressSealReceipt(receipt))
    .preIngressSealEnvelope(envelope)
}

.preIngressIdentity <- function(capability_id, project_id, omic_label) {
    capability <- workflowCapabilityCatalogue()[[capability_id]]
    c(
        list(project_id = project_id, omic_label = omic_label),
        capability$identity
    )
}

.preIngressContract <- function(
    capability_id,
    descriptor_digest,
    project_id,
    omic_label,
    calls,
    container_type = "plain",
    archive_measure = NULL,
    probe_override = list()
) {
    newWorkflowPreIngressContract(
        contract_id = paste0("test.preingress.", capability_id),
        capability_ids = capability_id,
        member_resolver = function(input) {
            calls$member <- calls$member + 1L
            lapply(seq_along(input$paths), function(index) {
                list(
                    member_id = paste0("member_", index),
                    path = input$paths[[index]],
                    container_type = container_type,
                    semantic_order = as.integer(index)
                )
            })
        },
        identity_probe = function(members) {
            calls$probe <- calls$probe + 1L
            base <- list(
                identity = .preIngressIdentity(
                    capability_id,
                    project_id,
                    omic_label
                ),
                capability_id = capability_id,
                capability_version = "1.0.0",
                descriptor_digest = descriptor_digest,
                bytes_read = 16L,
                lines_read = 1L,
                schema_valid = TRUE,
                ambiguous = FALSE,
                complete_payload_materialized = FALSE
            )
            utils::modifyList(base, probe_override)
        },
        receipt_provider = function(injected) {
            calls$receipt <- calls$receipt + 1L
            workflowPolicyLoadEnvelope(injected = injected)
        },
        archive_measure = archive_measure,
        max_probe_bytes = 64L,
        max_probe_lines = 4L
    )
}

.preIngressCalls <- function() {
    calls <- new.env(parent = emptyenv())
    calls$member <- 0L
    calls$probe <- 0L
    calls$receipt <- 0L
    calls$importer <- 0L
    calls
}

test_that("explicit memory is pre-ingress and artifact-resource inert", {
    root <- withr::local_tempdir()
    path <- file.path(root, "input.tsv")
    writeLines("feature\tvalue\nA\t1", path)
    capability_id <- "metabolomics.custom.metabolite.standard.v1"
    descriptor <- artifactWorkflowDescriptorCatalogue()$descriptors[[capability_id]]
    calls <- .preIngressCalls()
    contract <- .preIngressContract(
        capability_id,
        descriptor$descriptor_digest,
        basename(root),
        "metabolomics-test",
        calls
    )

    outcome <- workflowPreIngressResolve(
        contract,
        list(paths = path),
        root,
        requested_backend = "memory"
    )
    expect_identical(outcome$status, "memory")
    expect_identical(outcome$reason_code, "explicit_memory_without_preingress_probe")
    expect_identical(c(calls$member, calls$probe, calls$receipt), c(0L, 0L, 0L))
    expect_false(outcome$artifact_resources_opened)
})

test_that("exact proposed pre-ingress binds every descriptor identity", {
    capabilities <- c(
        "proteomics.diann.peptide.dia.v1",
        "proteomics.maxquant.protein.lfq.v1",
        "metabolomics.custom.metabolite.standard.v1",
        "lipidomics.lipidsearch.lipid.standard.v1"
    )
    descriptors <- artifactWorkflowDescriptorCatalogue()$descriptors
    for (capability_id in capabilities) {
        root <- withr::local_tempdir()
        path <- file.path(root, "input.tsv")
        writeLines("feature\tvalue\nA\t1", path)
        calls <- .preIngressCalls()
        contract <- .preIngressContract(
            capability_id,
            descriptors[[capability_id]]$descriptor_digest,
            basename(root),
            paste0("label-", capability_id),
            calls
        )
        envelope <- .preIngressExactEnvelope(
            capability_id,
            descriptors[[capability_id]]$descriptor_digest
        )

        outcome <- workflowPreIngressResolve(
            contract,
            list(paths = path),
            root,
            injected_envelope = envelope
        )
        expect_identical(outcome$status, "bound", info = capability_id)
        expect_identical(outcome$effective_backend, "artifact")
        expect_identical(outcome$effective_rollout, "evict")
        expect_identical(outcome$token$identity$omic_type,
            descriptors[[capability_id]]$identity$omic_type)
        expect_true(outcome$token$measure$exact)
        expect_false(outcome$token$probe_evidence$complete_payload_materialized)
        expect_false(outcome$artifact_resources_opened)
        expect_identical(c(calls$member, calls$probe, calls$receipt), c(1L, 1L, 1L))
    }
})

test_that("installed legacy receipts cannot drive exact pre-ingress", {
    root <- withr::local_tempdir()
    path <- file.path(root, "input.tsv")
    writeLines("feature\tvalue\nA\t1", path)
    capability_id <- "metabolomics.custom.metabolite.standard.v1"
    descriptor <- artifactWorkflowDescriptorCatalogue()$descriptors[[capability_id]]
    calls <- .preIngressCalls()
    contract <- .preIngressContract(
        capability_id,
        descriptor$descriptor_digest,
        basename(root),
        "metabolomics-test",
        calls
    )

    outcome <- workflowPreIngressResolve(contract, list(paths = path), root)
    expect_identical(outcome$status, "bound")
    expect_identical(outcome$effective_backend, "memory")
    expect_identical(outcome$reason_code, "policy_exact_measure_unavailable")
    expect_false(outcome$artifact_resources_opened)
})

test_that("plain bundles count unique bytes and reject unsafe members", {
    root <- withr::local_tempdir()
    first <- file.path(root, "first.tsv")
    second <- file.path(root, "second.tsv")
    writeChar("12345", first, eos = NULL)
    writeChar("1234567", second, eos = NULL)
    capability_id <- "lipidomics.lipidsearch.lipid.standard.v1"
    descriptor <- artifactWorkflowDescriptorCatalogue()$descriptors[[capability_id]]
    calls <- .preIngressCalls()
    contract <- .preIngressContract(
        capability_id,
        descriptor$descriptor_digest,
        basename(root),
        "lipidomics-test",
        calls
    )
    envelope <- .preIngressExactEnvelope(
        capability_id,
        descriptor$descriptor_digest
    )

    outcome <- workflowPreIngressResolve(
        contract,
        list(paths = c(first, second)),
        root,
        injected_envelope = envelope
    )
    expect_identical(outcome$token$measure$bytes, 12)

    expect_error(
        workflowPreIngressResolve(
            contract,
            list(paths = c(first, first)),
            root,
            injected_envelope = envelope
        ),
        class = "multischolar_duplicate_workflow_preingress_member"
    )

    link <- file.path(root, "linked.tsv")
    skip_if_not(file.symlink(first, link), "symbolic links unavailable")
    expect_error(
        workflowPreIngressResolve(
            contract,
            list(paths = link),
            root,
            injected_envelope = envelope
        ),
        class = "multischolar_unsafe_workflow_preingress_member"
    )
})

test_that("unsupported archives resolve memory before receipt access", {
    root <- withr::local_tempdir()
    path <- file.path(root, "input.zip")
    writeChar("not-a-real-archive", path, eos = NULL)
    capability_id <- "proteomics.diann.peptide.dia.v1"
    descriptor <- artifactWorkflowDescriptorCatalogue()$descriptors[[capability_id]]
    calls <- .preIngressCalls()
    contract <- .preIngressContract(
        capability_id,
        descriptor$descriptor_digest,
        basename(root),
        "proteomics-test",
        calls,
        container_type = "archive"
    )

    outcome <- workflowPreIngressResolve(contract, list(paths = path), root)
    expect_identical(outcome$status, "memory")
    expect_identical(outcome$reason_code, "exact_archive_size_unavailable")
    expect_identical(calls$receipt, 0L)
    expect_false(outcome$artifact_resources_opened)
})

test_that("unavailable format identities resolve memory before receipts", {
    cases <- list(
        unknown_workflow_format = list(capability_id = "unknown.format.v1"),
        unsupported_workflow_schema = list(schema_valid = FALSE),
        ambiguous_workflow_format = list(ambiguous = TRUE)
    )
    for (reason_code in names(cases)) {
        root <- withr::local_tempdir()
        path <- file.path(root, "input.tsv")
        writeLines("feature\tvalue\nA\t1", path)
        capability_id <- "metabolomics.custom.metabolite.standard.v1"
        descriptor <- artifactWorkflowDescriptorCatalogue()$descriptors[[
            capability_id
        ]]
        calls <- .preIngressCalls()
        contract <- .preIngressContract(
            capability_id,
            descriptor$descriptor_digest,
            basename(root),
            "metabolomics-test",
            calls,
            probe_override = cases[[reason_code]]
        )

        outcome <- workflowPreIngressResolve(
            contract,
            list(paths = path),
            root
        )
        expect_identical(outcome$status, "memory", info = reason_code)
        expect_identical(outcome$reason_code, reason_code)
        expect_identical(calls$receipt, 0L)
        expect_false(outcome$artifact_resources_opened)
    }
})

test_that("bounded probes and importer handoff reject mutation", {
    root <- withr::local_tempdir()
    path <- file.path(root, "input.tsv")
    writeLines("feature\tvalue\nA\t1", path)
    capability_id <- "metabolomics.custom.metabolite.standard.v1"
    descriptor <- artifactWorkflowDescriptorCatalogue()$descriptors[[capability_id]]
    calls <- .preIngressCalls()
    contract <- .preIngressContract(
        capability_id,
        descriptor$descriptor_digest,
        basename(root),
        "metabolomics-test",
        calls
    )
    envelope <- .preIngressExactEnvelope(
        capability_id,
        descriptor$descriptor_digest
    )
    outcome <- workflowPreIngressResolve(
        contract,
        list(paths = path),
        root,
        injected_envelope = envelope
    )
    cat("changed", file = path, append = TRUE)
    expect_error(
        workflowPreIngressConsume(
            contract,
            outcome,
            list(paths = path),
            function(input, token) {
                calls$importer <- calls$importer + 1L
                input
            }
        ),
        class = "multischolar_mutated_workflow_preingress_member"
    )
    expect_identical(calls$importer, 0L)

    unbounded_calls <- .preIngressCalls()
    unbounded <- .preIngressContract(
        capability_id,
        descriptor$descriptor_digest,
        basename(root),
        "metabolomics-test",
        unbounded_calls,
        probe_override = list(bytes_read = 65L)
    )
    expect_error(
        workflowPreIngressResolve(
            unbounded,
            list(paths = path),
            root,
            injected_envelope = envelope
        ),
        class = "multischolar_invalid_workflow_preingress_probe"
    )
})

test_that("valid handoff invokes one importer and binds context first", {
    root <- withr::local_tempdir()
    path <- file.path(root, "input.tsv")
    writeLines("feature\tvalue\nA\t1", path)
    capability_id <- "lipidomics.lipidsearch.lipid.standard.v1"
    descriptor <- artifactWorkflowDescriptorCatalogue()$descriptors[[capability_id]]
    calls <- .preIngressCalls()
    contract <- .preIngressContract(
        capability_id,
        descriptor$descriptor_digest,
        basename(root),
        "lipidomics-test",
        calls
    )
    envelope <- .preIngressExactEnvelope(
        capability_id,
        descriptor$descriptor_digest
    )
    outcome <- workflowPreIngressResolve(
        contract,
        list(paths = path),
        root,
        injected_envelope = envelope
    )
    context <- createWorkflowContext(
        list(base_dir = root, omic_label = "lipidomics-test"),
        "lipidomics",
        storage_policy = list(project_id = basename(root))
    )
    snapshot <- workflowPreIngressBindContext(
        context,
        outcome,
        buildArtifactPaths
    )
    expect_identical(snapshot$resolution$effective_backend, "artifact")
    expect_true(context$isBound())

    consumed <- workflowPreIngressConsume(
        contract,
        outcome,
        list(paths = path),
        function(input, token) {
            calls$importer <- calls$importer + 1L
            expect_true(context$isBound())
            expect_identical(token$state, "bound")
            list(row_count = 1L)
        }
    )
    expect_true(consumed$imported)
    expect_identical(consumed$outcome$status, "consumed")
    expect_identical(consumed$result$row_count, 1L)
    expect_identical(calls$importer, 1L)
})

test_that("pre-ingress tokens reject policy and identity tampering", {
    root <- withr::local_tempdir()
    path <- file.path(root, "input.tsv")
    writeLines("feature\tvalue\nA\t1", path)
    capability_id <- "lipidomics.lipidsearch.lipid.standard.v1"
    descriptor <- artifactWorkflowDescriptorCatalogue()$descriptors[[
        capability_id
    ]]
    calls <- .preIngressCalls()
    contract <- .preIngressContract(
        capability_id,
        descriptor$descriptor_digest,
        basename(root),
        "lipidomics-test",
        calls
    )
    outcome <- workflowPreIngressResolve(
        contract,
        list(paths = path),
        root,
        injected_envelope = .preIngressExactEnvelope(
            capability_id,
            descriptor$descriptor_digest
        )
    )
    mutations <- list(
        backend = function(token) {
            token$decision$effective_backend <- "memory"
            token
        },
        receipt = function(token) {
            token$receipt_binding$receipt_digest <- strrep("0", 64L)
            token
        },
        identity = function(token) {
            token$identity$omic_type <- "proteomics"
            token
        },
        process = function(token) {
            token$creator_pid <- token$creator_pid + 1L
            token
        }
    )
    for (name in names(mutations)) {
        tampered <- outcome
        tampered$token <- mutations[[name]](tampered$token)
        expect_error(
            validateWorkflowPreIngressToken(tampered$token, contract),
            class = "multischolar_invalid_workflow_preingress_token",
            info = name
        )
    }
})
