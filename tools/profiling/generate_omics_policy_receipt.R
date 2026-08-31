#!/usr/bin/env Rscript

policyGeneratorRepoRoot <- function(start = getwd()) {
    path <- normalizePath(start, mustWork = TRUE)
    repeat {
        if (file.exists(file.path(path, "DESCRIPTION")) &&
            dir.exists(file.path(path, "R"))) {
            return(path)
        }
        parent <- dirname(path)
        if (identical(parent, path)) stop("Cannot locate repository root")
        path <- parent
    }
}

.POLICY_GENERATOR_ROOT <- policyGeneratorRepoRoot()

source(file.path(
    .POLICY_GENERATOR_ROOT,
    "R",
    "utils_workflow_policy_receipt.R"
), local = FALSE)

policyGeneratorFileDigest <- function(path) {
    digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

policyGeneratorPath <- function(path) {
    if (startsWith(path, .Platform$file.sep)) return(path)
    file.path(.POLICY_GENERATOR_ROOT, path)
}

policyGeneratorBinding <- function(binding_id, path) {
    list(
        binding_id = binding_id,
        path_or_uri = path,
        sha256 = policyGeneratorFileDigest(policyGeneratorPath(path))
    )
}

policyGeneratorValidateSource <- function(binding) {
    path <- policyGeneratorPath(binding$path)
    valid <- file.exists(path) &&
        identical(policyGeneratorFileDigest(path), binding$sha256)
    if (!valid) stop("Policy generator source binding differs: ", binding$path)
    invisible(binding)
}

policyGeneratorRollback <- function(reconciliation) {
    rollback <- reconciliation$rollback
    list(
        new_work_backend = rollback$new_work_backend,
        existing_artifact_mode = rollback$existing_artifact_mode,
        implicit_memory_migration =
            rollback$existing_memory_projects_migrate_implicitly,
        delete_artifacts = rollback$delete_immutable_artifacts
    )
}

policyGeneratorReceipt <- function(
    decision,
    reconciliation,
    reconciliation_path,
    schema_v1_path,
    runtime_source
) {
    dynamic <- decision$dynamic_runtime_decision
    threshold <- as.integer(dynamic$threshold_bytes)
    receipt <- list(
        schema = "multischolar.omics_auto_policy_receipt",
        schema_version = "2.0.0",
        schema_predecessor = policyGeneratorBinding(
            "policy_receipt_schema_v1",
            schema_v1_path
        ),
        receipt_id = paste0(
            "multischolar.omics_auto_policy_receipt.current.",
            decision$capability_id,
            ".v2"
        ),
        receipt_kind = "reconciled_current",
        owner_ticket_id = "OMICS-ART-069",
        capability_id = decision$capability_id,
        capability_version = NULL,
        descriptor_digest = decision$descriptor_digest,
        platforms = list(dynamic$platform),
        size_measure = list(
            measure_id = dynamic$legacy_size_measure_id,
            unit = "byte",
            exact = FALSE,
            available_before_full_parse = FALSE
        ),
        threshold_bytes = threshold,
        below_threshold = list(
            backend = dynamic$below_threshold_backend,
            rollout = "none"
        ),
        at_or_above_threshold = list(
            backend = dynamic$at_or_above_threshold_backend,
            rollout = dynamic$at_or_above_threshold_rollout
        ),
        decision = "engineering_current_baseline",
        terminal_reason =
            "legacy_transition_preserves_current_routing_until_final_receipt",
        claim_scope = list(
            support_tier = "scientifically_supported",
            project_scope = "project_specific",
            host_scope = "primary_host",
            publication_authority = FALSE
        ),
        authority_bindings = list(
            policyGeneratorBinding(
                "current_policy_reconciliation",
                reconciliation_path
            ),
            list(
                binding_id = "legacy_runtime_source",
                path_or_uri = runtime_source$path,
                sha256 = runtime_source$sha256
            ),
            list(
                binding_id = "descriptor_source",
                path_or_uri = decision$descriptor_source$path,
                sha256 = decision$descriptor_source$sha256
            )
        ),
        evidence_bindings = list(policyGeneratorBinding(
            "scale_workload",
            decision$scale_workload$path
        )),
        role_bindings = list(list(
            role_id = "runtime_policy_implementation",
            principal_id = "OMICS-ART-069",
            handoff_digest = policyGeneratorFileDigest(
                policyGeneratorPath(reconciliation_path)
            )
        )),
        rollback = policyGeneratorRollback(reconciliation),
        digest_contract = list(
            algorithm = "sha256",
            canonicalization = "canonical_json_sorted_keys_utf8_v1",
            excluded_field = "receipt_digest"
        ),
        receipt_digest = NULL
    )
    receipt_basis <- receipt
    receipt_basis$receipt_digest <- NULL
    receipt$receipt_digest <- workflowPolicyObjectDigest(receipt_basis)
    workflowPolicyValidateReceipt(receipt)
    receipt
}

policyGeneratorEnvelope <- function() {
    reconciliation_path <- paste0(
        "tests/testdata/omics-performance/",
        "current-policy-reconciliation-v1.json"
    )
    schema_v1_path <- paste0(
        "tests/testdata/omics-performance/",
        "policy-receipt-schema-v1.json"
    )
    schema_v2_path <- paste0(
        "tests/testdata/omics-performance/",
        "policy-receipt-schema-v2.json"
    )
    generator_path <- "tools/profiling/generate_omics_policy_receipt.R"
    reconciliation <- jsonlite::read_json(
        policyGeneratorPath(reconciliation_path),
        simplifyVector = FALSE
    )
    runtime_source <- reconciliation$current_runtime_source
    if (!workflowPolicyDigestValid(runtime_source$sha256) ||
        !workflowPolicyScalarString(runtime_source$path)) {
        stop("Legacy runtime predecessor binding is malformed")
    }
    lapply(
        reconciliation$implemented_decisions,
        function(decision) {
            policyGeneratorValidateSource(decision$descriptor_source)
            policyGeneratorValidateSource(decision$scale_workload)
        }
    )
    decisions <- reconciliation$implemented_decisions
    receipts <- lapply(decisions, policyGeneratorReceipt,
        reconciliation = reconciliation,
        reconciliation_path = reconciliation_path,
        schema_v1_path = schema_v1_path,
        runtime_source = runtime_source
    )
    envelope <- list(
        schema = "multischolar.omics_auto_policy_receipt_envelope",
        schema_version = "1.0.0",
        envelope_id = paste0(
            "multischolar.omics_auto_policy_receipt_envelope.",
            "2026-08-30.v1"
        ),
        owner_ticket_id = "OMICS-ART-069",
        generated_from = policyGeneratorBinding(
            "current_policy_reconciliation",
            reconciliation_path
        ),
        receipt_schema = policyGeneratorBinding(
            "policy_receipt_schema_v2",
            schema_v2_path
        ),
        consumer_contract_digest = workflowPolicyConsumerContractDigest(),
        receipts = receipts,
        generation_contract = list(
            generator_path = generator_path,
            generator_sha256 = policyGeneratorFileDigest(
                policyGeneratorPath(generator_path)
            ),
            copied_values_allowed = FALSE,
            additive_schema_successor = TRUE,
            production_defaults_changed = FALSE
        ),
        envelope_digest = NULL
    )
    envelope_basis <- envelope
    envelope_basis$envelope_digest <- NULL
    envelope$envelope_digest <- workflowPolicyObjectDigest(envelope_basis)
    workflowPolicyValidateEnvelope(envelope)
    envelope
}

policyGeneratorWriteJson <- function(value, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    temporary <- tempfile("policy-envelope-", tmpdir = dirname(path))
    on.exit(unlink(temporary, force = TRUE), add = TRUE)
    jsonlite::write_json(
        value,
        temporary,
        auto_unbox = TRUE,
        pretty = TRUE,
        null = "null",
        na = "null",
        digits = 17
    )
    if (!file.rename(temporary, path)) stop("Could not publish policy envelope")
    invisible(path)
}

policyGeneratorArgs <- function(argv) {
    values <- list(
        output = file.path(
            .POLICY_GENERATOR_ROOT,
            "inst",
            "extdata",
            "omics-auto-policy-receipts-v2.json"
        ),
        check = "false"
    )
    index <- 1L
    while (index <= length(argv)) {
        token <- argv[[index]]
        key <- gsub("-", "_", sub("^--", "", token), fixed = TRUE)
        if (!startsWith(token, "--") || !key %in% names(values) ||
            index == length(argv)) {
            stop("Policy generator arguments are invalid")
        }
        values[[key]] <- argv[[index + 1L]]
        index <- index + 2L
    }
    values$check <- identical(tolower(values$check), "true")
    values
}

policyGeneratorMain <- function(argv = commandArgs(trailingOnly = TRUE)) {
    args <- policyGeneratorArgs(argv)
    envelope <- policyGeneratorEnvelope()
    if (args$check) {
        if (!file.exists(args$output)) stop("Installed policy envelope is missing")
        installed <- jsonlite::read_json(args$output, simplifyVector = FALSE)
        if (!identical(
            workflowPolicyObjectDigest(installed),
            workflowPolicyObjectDigest(envelope)
        )) {
            stop("Installed policy envelope differs from generated authority")
        }
    } else {
        policyGeneratorWriteJson(envelope, args$output)
    }
    cat(envelope$envelope_digest, "\n")
    invisible(0L)
}

if (identical(environment(), globalenv()) && !interactive()) {
    status <- tryCatch(
        policyGeneratorMain(),
        error = function(error) {
            message(conditionMessage(error))
            1L
        }
    )
    quit(status = as.integer(status), save = "no")
}
