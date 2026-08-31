# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

artifactSettlementAbort <- function(message, class, ...) {
    rlang::abort(
        message,
        class = c(class, "multischolar_artifact_settlement_error"),
        ...
    )
}

artifactSettlementOperationNames <- function() {
    c(
        "current", "encode", "bounded_validate", "publish", "lean_snapshot",
        "release_sources", "close_resources", "rollback", "cleanup_prepared",
        "hydration_active"
    )
}

artifactSettlementValidateParityContract <- function(contract) {
    workflowPolicyRequireNames(contract, c(
        "worker_id", "process_isolation", "required", "input_binding_ids",
        "output_digest_required", "candidate_process_allowed"
    ), "artifact settlement parity contract")
    valid <- workflowPolicyScalarString(contract$worker_id) &&
        identical(contract$process_isolation, "fresh_R_process") &&
        isTRUE(contract$required) && is.character(contract$input_binding_ids) &&
        length(contract$input_binding_ids) > 0L &&
        !anyDuplicated(contract$input_binding_ids) &&
        isTRUE(contract$output_digest_required) &&
        !isTRUE(contract$candidate_process_allowed)
    if (!valid) {
        artifactSettlementAbort(
            "artifact settlement parity contract is invalid",
            "multischolar_invalid_artifact_settlement_contract"
        )
    }
    invisible(contract)
}

newArtifactSettlementContract <- function(
    contract_id,
    descriptor,
    capability_id,
    source_fields,
    stage_roles,
    compatibility_strategy,
    consumer_inventory_digest,
    resolver,
    parity_worker_contract,
    operations,
    max_lean_bytes = 1024^2
) {
    descriptor <- validateArtifactWorkflowDescriptor(descriptor)
    valid <- workflowPolicyScalarString(contract_id) &&
        identical(descriptor$descriptor_id, capability_id) &&
        is.character(source_fields) && length(source_fields) > 0L &&
        !anyNA(source_fields) && all(nzchar(source_fields)) &&
        !anyDuplicated(source_fields) && is.character(stage_roles) &&
        length(stage_roles) > 0L && !anyNA(stage_roles) &&
        all(nzchar(stage_roles)) && !anyDuplicated(stage_roles) &&
        workflowPolicyScalarString(compatibility_strategy) &&
        is.function(consumer_inventory_digest) && is.function(resolver) &&
        is.list(operations) &&
        identical(names(operations), artifactSettlementOperationNames()) &&
        all(vapply(operations, is.function, logical(1))) &&
        is.numeric(max_lean_bytes) && length(max_lean_bytes) == 1L &&
        is.finite(max_lean_bytes) && max_lean_bytes > 0
    if (!valid) {
        artifactSettlementAbort(
            "artifact settlement contract is malformed",
            "multischolar_invalid_artifact_settlement_contract"
        )
    }
    artifactSettlementValidateParityContract(parity_worker_contract)
    digest <- consumer_inventory_digest()
    if (!workflowPolicyDigestValid(digest)) {
        artifactSettlementAbort(
            "artifact settlement consumer inventory digest is invalid",
            "multischolar_invalid_artifact_settlement_contract"
        )
    }
    structure(
        list(
            contract_id = contract_id,
            descriptor = descriptor,
            capability_id = capability_id,
            source_fields = source_fields,
            stage_roles = stage_roles,
            compatibility_strategy = compatibility_strategy,
            consumer_inventory_digest = consumer_inventory_digest,
            resolver = resolver,
            parity_worker_contract = parity_worker_contract,
            operations = operations,
            max_lean_bytes = as.numeric(max_lean_bytes)
        ),
        class = "ArtifactSettlementContract"
    )
}

validateArtifactSettlementContract <- function(contract) {
    expected <- c(
        "contract_id", "descriptor", "capability_id", "source_fields",
        "stage_roles", "compatibility_strategy", "consumer_inventory_digest",
        "resolver", "parity_worker_contract", "operations", "max_lean_bytes"
    )
    valid <- inherits(contract, "ArtifactSettlementContract") &&
        is.list(contract) && identical(names(contract), expected) &&
        workflowPolicyScalarString(contract$contract_id) &&
        identical(
            contract$descriptor$descriptor_id,
            contract$capability_id
        ) && is.character(contract$source_fields) &&
        length(contract$source_fields) > 0L &&
        !anyNA(contract$source_fields) &&
        all(nzchar(contract$source_fields)) &&
        !anyDuplicated(contract$source_fields) &&
        is.character(contract$stage_roles) &&
        length(contract$stage_roles) > 0L &&
        !anyNA(contract$stage_roles) &&
        all(nzchar(contract$stage_roles)) &&
        !anyDuplicated(contract$stage_roles) &&
        workflowPolicyScalarString(contract$compatibility_strategy) &&
        is.function(contract$consumer_inventory_digest) &&
        is.function(contract$resolver) && is.list(contract$operations) &&
        identical(
            names(contract$operations),
            artifactSettlementOperationNames()
        ) && all(vapply(contract$operations, is.function, logical(1))) &&
        is.numeric(contract$max_lean_bytes) &&
        length(contract$max_lean_bytes) == 1L &&
        is.finite(contract$max_lean_bytes) &&
        contract$max_lean_bytes > 0
    if (!valid) {
        artifactSettlementAbort(
            "artifact settlement contract is invalid",
            "multischolar_invalid_artifact_settlement_contract"
        )
    }
    validateArtifactWorkflowDescriptor(contract$descriptor)
    artifactSettlementValidateParityContract(contract$parity_worker_contract)
    contract
}

artifactSettlementForbiddenClass <- function(value) {
    is.data.frame(value) || is.matrix(value) || isS4(value) ||
        is.environment(value) || typeof(value) == "externalptr" ||
        inherits(value, c(
            "DBIConnection", "DBIResult", "ArrowObject", "duckdb_connection",
            "ArtifactResourceScope", "ArtifactQuerySession"
        ))
}

artifactSettlementAssertLeanValue <- function(value, owner = "lean snapshot") {
    if (artifactSettlementForbiddenClass(value)) {
        artifactSettlementAbort(
            paste(owner, "contains a complete or process-bound value"),
            "multischolar_nonlean_artifact_settlement"
        )
    }
    if (is.list(value)) {
        lapply(value, artifactSettlementAssertLeanValue, owner)
    }
    invisible(value)
}

artifactSettlementValidateLeanSnapshot <- function(snapshot, max_bytes) {
    expected <- c(
        "lineage", "current", "bounded_metadata", "audit_state", "config_state",
        "descriptor_pin", "rollback_checkpoint", "locators"
    )
    workflowPolicyRequireNames(snapshot, expected, "artifact lean snapshot")
    forbidden_names <- c(
        "data_tbl", "data_cln", "payload", "hydrated", "result", "writer",
        "lock", "connection", "query", "cache"
    )
    names_recursive <- function(value) {
        if (!is.list(value)) return(character())
        c(names(value), unlist(lapply(value, names_recursive), use.names = FALSE))
    }
    ownership_fields <- snapshot[setdiff(names(snapshot), "locators")]
    nested_names <- names_recursive(ownership_fields)
    valid <- as.numeric(utils::object.size(snapshot)) <= max_bytes &&
        !any(tolower(nested_names) %in% forbidden_names)
    artifactSettlementAssertLeanValue(snapshot)
    if (!valid) {
        artifactSettlementAbort(
            "artifact lean snapshot is unbounded or contains payload ownership",
            "multischolar_nonlean_artifact_settlement"
        )
    }
    invisible(snapshot)
}

artifactSettlementValidatePrepared <- function(prepared, contract) {
    expected <- c(
        "intent_id", "generation_id", "refs", "bounded_metadata",
        "temporary_paths", "complete_payload_hydrated"
    )
    workflowPolicyRequireNames(prepared, expected, "artifact prepared settlement")
    refs <- prepared$refs
    valid <- workflowPolicyScalarString(prepared$intent_id) &&
        workflowPolicyScalarString(prepared$generation_id) && is.list(refs) &&
        setequal(names(refs), contract$stage_roles) &&
        all(vapply(refs, is.list, logical(1))) &&
        is.list(prepared$bounded_metadata) &&
        is.character(prepared$temporary_paths) &&
        !isTRUE(prepared$complete_payload_hydrated)
    if (!valid) {
        artifactSettlementAbort(
            "artifact prepared settlement is invalid",
            "multischolar_invalid_artifact_settlement_preparation"
        )
    }
    artifactSettlementAssertLeanValue(refs, "artifact prepared refs")
    invisible(prepared)
}

artifactSettlementValidateBounded <- function(validation) {
    workflowPolicyRequireNames(validation, c(
        "valid", "refs_validated", "bounded_metadata_validated",
        "complete_payload_hydrated"
    ), "artifact bounded validation")
    valid <- isTRUE(validation$valid) && isTRUE(validation$refs_validated) &&
        isTRUE(validation$bounded_metadata_validated) &&
        !isTRUE(validation$complete_payload_hydrated)
    if (!valid) {
        artifactSettlementAbort(
            "artifact bounded validation failed",
            "multischolar_invalid_artifact_settlement_validation"
        )
    }
    invisible(validation)
}

artifactSettlementValidatePublication <- function(publication, prepared) {
    workflowPolicyRequireNames(publication, c(
        "current", "generation_id", "refs", "published_atomically"
    ), "artifact settlement publication")
    valid <- identical(publication$generation_id, prepared$generation_id) &&
        identical(
            workflowPolicyObjectDigest(publication$refs),
            workflowPolicyObjectDigest(prepared$refs)
        ) && isTRUE(publication$published_atomically)
    if (!valid) {
        artifactSettlementAbort(
            "artifact settlement publication is invalid",
            "multischolar_invalid_artifact_settlement_publication"
        )
    }
    invisible(publication)
}

artifactSettlementSourceValues <- function(workflow_data, source_fields) {
    values <- lapply(source_fields, function(name) workflow_data[[name]])
    names(values) <- source_fields
    if (any(vapply(values, is.null, logical(1)))) {
        artifactSettlementAbort(
            "artifact settlement source ownership is incomplete",
            "multischolar_missing_artifact_settlement_source"
        )
    }
    values
}

artifactSettlementRestoreSources <- function(
    workflow_data,
    source_values,
    source_fields
) {
    for (name in source_fields) workflow_data[[name]] <- source_values[[name]]
    invisible(TRUE)
}

artifactSettlementVerifySourcesReleased <- function(workflow_data, fields) {
    all(vapply(fields, function(name) is.null(workflow_data[[name]]), logical(1)))
}

artifactSettlementVerifyResourcesClosed <- function(info) {
    required <- c(
        "registry_connections", "query_handles", "results", "writers", "locks"
    )
    is.list(info) && all(required %in% names(info)) &&
        all(vapply(info[required], function(value) {
            is.numeric(value) && length(value) == 1L && identical(value, 0L)
        }, logical(1)))
}

artifactSettlementFailure <- function(
    error,
    contract,
    workflow_data,
    source_values,
    previous_current,
    publication,
    prepared,
    transitions
) {
    rollback_error <- tryCatch({
        contract$operations$rollback(
            previous_current,
            publication,
            prepared
        )
        NULL
    }, error = function(error) error)
    cleanup_error <- tryCatch({
        contract$operations$cleanup_prepared(prepared)
        NULL
    }, error = function(error) error)
    artifactSettlementRestoreSources(
        workflow_data,
        source_values,
        contract$source_fields
    )
    transitions[[length(transitions) + 1L]] <- "rolled_back"
    artifactSettlementAbort(
        "artifact settlement transaction failed and was rolled back",
        "multischolar_artifact_settlement_rolled_back",
        parent = error,
        transitions = transitions,
        rollback_error = rollback_error,
        cleanup_error = cleanup_error
    )
}

settleArtifactWorkflow <- function(workflow_data, contract) {
    contract <- validateArtifactSettlementContract(contract)
    consumer_inventory_digest <- contract$consumer_inventory_digest()
    if (!workflowPolicyDigestValid(consumer_inventory_digest)) {
        artifactSettlementAbort(
            "artifact settlement consumer inventory digest is invalid",
            "multischolar_invalid_artifact_settlement_contract"
        )
    }
    if (!is.environment(workflow_data) &&
        !inherits(workflow_data, "reactivevalues")) {
        artifactSettlementAbort(
            "artifact settlement requires mutable workflow data",
            "multischolar_invalid_artifact_settlement_source"
        )
    }
    if (isTRUE(contract$operations$hydration_active(workflow_data))) {
        artifactSettlementAbort(
            "artifact settlement cannot start with active hydration",
            "multischolar_overlapping_artifact_hydration"
        )
    }
    source_values <- artifactSettlementSourceValues(
        workflow_data,
        contract$source_fields
    )
    previous_current <- contract$operations$current()
    transitions <- list("idle")
    prepared <- NULL
    publication <- NULL
    result <- tryCatch({
        prepared <- contract$operations$encode(source_values)
        artifactSettlementValidatePrepared(prepared, contract)
        transitions[[length(transitions) + 1L]] <- "prepared"

        bounded <- contract$operations$bounded_validate(prepared)
        artifactSettlementValidateBounded(bounded)
        transitions[[length(transitions) + 1L]] <- "bounded_validated"

        publication <- contract$operations$publish(
            prepared,
            previous_current
        )
        artifactSettlementValidatePublication(publication, prepared)
        transitions[[length(transitions) + 1L]] <- "published"

        snapshot <- contract$operations$lean_snapshot(
            publication,
            prepared$bounded_metadata
        )
        artifactSettlementValidateLeanSnapshot(
            snapshot,
            contract$max_lean_bytes
        )

        contract$operations$release_sources(
            workflow_data,
            contract$source_fields
        )
        if (!artifactSettlementVerifySourcesReleased(
            workflow_data,
            contract$source_fields
        )) {
            artifactSettlementAbort(
                "artifact settlement source release is incomplete",
                "multischolar_incomplete_artifact_settlement_release"
            )
        }
        transitions[[length(transitions) + 1L]] <- "source_released"

        resources <- contract$operations$close_resources(workflow_data)
        if (!artifactSettlementVerifyResourcesClosed(resources)) {
            artifactSettlementAbort(
                "artifact settlement resources remain open",
                "multischolar_incomplete_artifact_settlement_release"
            )
        }
        transitions[[length(transitions) + 1L]] <- "settled"
        list(
            status = "settled",
            contract_id = contract$contract_id,
            descriptor_id = contract$descriptor$descriptor_id,
            capability_id = contract$capability_id,
            generation_id = prepared$generation_id,
            refs = publication$refs,
            snapshot = snapshot,
            transitions = transitions,
            consumer_inventory_digest = consumer_inventory_digest,
            parity_worker_contract = contract$parity_worker_contract,
            complete_parity_executed_in_commit_process = FALSE,
            source_fields_released = TRUE,
            resources_closed = TRUE,
            promotion_authority = FALSE
        )
    }, error = function(error) error)
    if (inherits(result, "error")) {
        artifactSettlementFailure(
            result,
            contract,
            workflow_data,
            source_values,
            previous_current,
            publication,
            prepared,
            transitions
        )
    }
    result
}
