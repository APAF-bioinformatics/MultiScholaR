# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

workflowPolicyAbort <- function(message, class, ...) {
    rlang::abort(
        message,
        class = c(class, "multischolar_workflow_policy_error"),
        ...
    )
}

workflowPolicyScalarString <- function(value) {
    is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value)
}

workflowPolicyDigestValid <- function(value) {
    workflowPolicyScalarString(value) && grepl("^[0-9a-f]{64}$", value)
}

workflowPolicySortObject <- function(value) {
    if (!is.list(value)) return(value)
    if (!is.null(names(value))) {
        value <- value[order(names(value), method = "radix")]
    }
    lapply(value, workflowPolicySortObject)
}

workflowPolicyObjectDigest <- function(value) {
    encoded <- jsonlite::toJSON(
        workflowPolicySortObject(value),
        auto_unbox = TRUE,
        null = "null",
        na = "null",
        digits = 17,
        pretty = FALSE
    )
    digest::digest(encoded, algo = "sha256", serialize = FALSE)
}

workflowPolicyRequireNames <- function(value, expected, owner) {
    if (!is.list(value) || !setequal(names(value), expected)) {
        workflowPolicyAbort(
            paste(owner, "fields differ"),
            "multischolar_invalid_workflow_policy_receipt"
        )
    }
    invisible(value)
}

workflowPolicyValidateBinding <- function(binding, owner) {
    workflowPolicyRequireNames(
        binding,
        c("binding_id", "path_or_uri", "sha256"),
        owner
    )
    valid <- workflowPolicyScalarString(binding$binding_id) &&
        workflowPolicyScalarString(binding$path_or_uri) &&
        workflowPolicyDigestValid(binding$sha256)
    if (!valid) {
        workflowPolicyAbort(
            paste(owner, "is invalid"),
            "multischolar_invalid_workflow_policy_receipt"
        )
    }
    invisible(binding)
}

workflowPolicyReceiptRequiredFields <- function() {
    c(
        "schema", "schema_version", "schema_predecessor", "receipt_id",
        "receipt_kind", "owner_ticket_id", "capability_id",
        "capability_version", "descriptor_digest", "platforms", "size_measure",
        "threshold_bytes", "below_threshold", "at_or_above_threshold",
        "decision", "terminal_reason", "claim_scope", "authority_bindings",
        "evidence_bindings", "role_bindings", "rollback", "digest_contract",
        "receipt_digest"
    )
}

workflowPolicyReceiptOptionalFields <- function() {
    c(
        "project_ids", "pilot_split_id", "holdout_split_id",
        "threshold_grid_id", "estimand_ids", "candidate_revision",
        "gate_results"
    )
}

workflowPolicyValidateReceiptFields <- function(receipt) {
    required <- workflowPolicyReceiptRequiredFields()
    allowed <- c(required, workflowPolicyReceiptOptionalFields())
    valid <- is.list(receipt) && all(required %in% names(receipt)) &&
        all(names(receipt) %in% allowed) && !anyDuplicated(names(receipt))
    if (!valid) {
        workflowPolicyAbort(
            "workflow policy receipt fields differ",
            "multischolar_invalid_workflow_policy_receipt"
        )
    }
    invisible(receipt)
}

workflowPolicyValidateMeasure <- function(measure, receipt) {
    workflowPolicyRequireNames(
        measure,
        c("measure_id", "unit", "exact", "available_before_full_parse"),
        "workflow policy size measure"
    )
    legacy <- identical(measure$measure_id, "legacy_r_object_size_x2.v1")
    exact <- identical(
        measure$measure_id,
        "total_uncompressed_input_bytes_v1"
    )
    valid <- identical(measure$unit, "byte") && (legacy || exact)
    if (legacy) {
        valid <- valid && !isTRUE(measure$exact) &&
            !isTRUE(measure$available_before_full_parse) &&
            identical(receipt$receipt_kind, "reconciled_current") &&
            identical(receipt$decision, "engineering_current_baseline") &&
            identical(receipt$owner_ticket_id, "OMICS-ART-069") &&
            !isTRUE(receipt$claim_scope$publication_authority)
    } else {
        valid <- valid && isTRUE(measure$exact) &&
            isTRUE(measure$available_before_full_parse)
    }
    if (!valid) {
        workflowPolicyAbort(
            "workflow policy size measure is inadmissible",
            "multischolar_invalid_workflow_policy_measure"
        )
    }
    invisible(measure)
}

workflowPolicyValidateThresholdBranches <- function(receipt) {
    workflowPolicyRequireNames(
        receipt$below_threshold,
        c("backend", "rollout"),
        "below-threshold workflow policy"
    )
    workflowPolicyRequireNames(
        receipt$at_or_above_threshold,
        c("backend", "rollout"),
        "at-or-above workflow policy"
    )
    valid <- identical(receipt$below_threshold$backend, "memory") &&
        identical(receipt$below_threshold$rollout, "none") &&
        receipt$at_or_above_threshold$backend %in% c("memory", "artifact") &&
        receipt$at_or_above_threshold$rollout %in%
            c("none", "dual_write", "read_through", "evict")
    if (!valid) {
        workflowPolicyAbort(
            "workflow policy threshold branches are invalid",
            "multischolar_invalid_workflow_policy_decision"
        )
    }
    invisible(receipt)
}

workflowPolicyValidateDecision <- function(receipt) {
    workflowPolicyValidateThresholdBranches(receipt)
    no_promotion <- identical(receipt$decision, "confirmed_no_promotion")
    threshold_valid <- is.numeric(receipt$threshold_bytes) &&
        length(receipt$threshold_bytes) == 1L &&
        is.finite(receipt$threshold_bytes) && receipt$threshold_bytes > 0
    valid <- if (no_promotion) {
        is.null(receipt$threshold_bytes) &&
            identical(receipt$at_or_above_threshold$backend, "memory") &&
            identical(receipt$at_or_above_threshold$rollout, "none") &&
            workflowPolicyScalarString(receipt$terminal_reason)
    } else {
        threshold_valid
    }
    if (identical(receipt$decision, "engineering_current_baseline")) {
        valid <- valid &&
            identical(receipt$at_or_above_threshold$backend, "artifact") &&
            identical(receipt$at_or_above_threshold$rollout, "evict") &&
            !isTRUE(receipt$claim_scope$publication_authority)
    }
    if (receipt$decision %in% c("confirmed_promote", "final_installed")) {
        promotion_fields <- c(
            "project_ids", "pilot_split_id", "holdout_split_id",
            "threshold_grid_id", "estimand_ids", "candidate_revision",
            "gate_results"
        )
        valid <- valid && all(promotion_fields %in% names(receipt)) &&
            identical(receipt$at_or_above_threshold$backend, "artifact") &&
            identical(receipt$at_or_above_threshold$rollout, "evict") &&
            isTRUE(receipt$claim_scope$publication_authority) &&
            identical(
                receipt$size_measure$measure_id,
                "total_uncompressed_input_bytes_v1"
            ) && length(receipt$project_ids) >= 3L
    }
    if (!valid) {
        workflowPolicyAbort(
            "workflow policy decision is invalid",
            "multischolar_invalid_workflow_policy_decision"
        )
    }
    invisible(receipt)
}

workflowPolicyValidateReceipt <- function(receipt) {
    workflowPolicyValidateReceiptFields(receipt)
    workflowPolicyValidateBinding(
        receipt$schema_predecessor,
        "workflow policy schema predecessor"
    )
    workflowPolicyRequireNames(
        receipt$claim_scope,
        c(
            "support_tier", "project_scope", "host_scope",
            "publication_authority"
        ),
        "workflow policy claim scope"
    )
    workflowPolicyRequireNames(
        receipt$rollback,
        c(
            "new_work_backend", "existing_artifact_mode",
            "implicit_memory_migration", "delete_artifacts"
        ),
        "workflow policy rollback"
    )
    workflowPolicyRequireNames(
        receipt$digest_contract,
        c("algorithm", "canonicalization", "excluded_field"),
        "workflow policy digest contract"
    )
    bindings <- c(receipt$authority_bindings, receipt$evidence_bindings)
    lapply(bindings, workflowPolicyValidateBinding, "workflow policy binding")
    role_valid <- all(vapply(receipt$role_bindings, function(role) {
        is.list(role) && setequal(
            names(role),
            c("role_id", "principal_id", "handoff_digest")
        ) && workflowPolicyScalarString(role$role_id) &&
            workflowPolicyScalarString(role$principal_id) &&
            workflowPolicyDigestValid(role$handoff_digest)
    }, logical(1)))
    platforms <- unlist(receipt$platforms, use.names = FALSE)
    kinds <- c(
        "reconciled_current", "proposed_pilot",
        "confirmatory_decision", "final_installed"
    )
    decisions <- c(
        "engineering_current_baseline", "proposed_pilot",
        "confirmed_promote", "confirmed_no_promotion", "final_installed"
    )
    basis <- receipt
    basis$receipt_digest <- NULL
    valid <- identical(
        receipt$schema,
        "multischolar.omics_auto_policy_receipt"
    ) && identical(receipt$schema_version, "2.0.0") &&
        receipt$receipt_kind %in% kinds && receipt$decision %in% decisions &&
        grepl("^OMICS-ART-[0-9]{3}$", receipt$owner_ticket_id) &&
        workflowPolicyScalarString(receipt$receipt_id) &&
        workflowPolicyScalarString(receipt$capability_id) &&
        workflowPolicyDigestValid(receipt$descriptor_digest) &&
        length(platforms) > 0L && !anyDuplicated(platforms) &&
        all(nzchar(platforms)) && length(bindings) > 0L &&
        length(receipt$role_bindings) > 0L && role_valid &&
        identical(receipt$rollback$new_work_backend, "memory") &&
        identical(receipt$rollback$existing_artifact_mode, "read_only") &&
        !isTRUE(receipt$rollback$implicit_memory_migration) &&
        !isTRUE(receipt$rollback$delete_artifacts) &&
        identical(receipt$digest_contract$algorithm, "sha256") &&
        identical(
            receipt$digest_contract$canonicalization,
            "canonical_json_sorted_keys_utf8_v1"
        ) && identical(
            receipt$digest_contract$excluded_field,
            "receipt_digest"
        ) && identical(
            receipt$receipt_digest,
            workflowPolicyObjectDigest(basis)
        )
    if (!valid) {
        workflowPolicyAbort(
            "workflow policy receipt identity or digest differs",
            "multischolar_invalid_workflow_policy_receipt"
        )
    }
    workflowPolicyValidateMeasure(receipt$size_measure, receipt)
    workflowPolicyValidateDecision(receipt)
    invisible(receipt)
}

workflowPolicyConsumerContractDigest <- function() {
    functions <- c(
        "workflowPolicyValidateMeasure", "workflowPolicyValidateDecision",
        "workflowPolicyValidateReceipt", "workflowPolicyValidateEnvelope",
        "workflowPolicyResolve"
    )
    definitions <- lapply(functions, function(name) {
        fun <- get(name, mode = "function", inherits = TRUE)
        list(
            name = name,
            formals = lapply(formals(fun), function(value) {
                paste(deparse(value, width.cutoff = 500L), collapse = "\n")
            }),
            body = paste(deparse(body(fun), width.cutoff = 500L), collapse = "\n")
        )
    })
    workflowPolicyObjectDigest(definitions)
}

workflowPolicyValidateEnvelope <- function(envelope) {
    workflowPolicyRequireNames(envelope, c(
        "schema", "schema_version", "envelope_id", "owner_ticket_id",
        "generated_from", "receipt_schema", "consumer_contract_digest",
        "receipts", "generation_contract", "envelope_digest"
    ), "workflow policy envelope")
    workflowPolicyValidateBinding(
        envelope$generated_from,
        "workflow policy source reconciliation"
    )
    workflowPolicyValidateBinding(
        envelope$receipt_schema,
        "workflow policy receipt schema"
    )
    workflowPolicyRequireNames(envelope$generation_contract, c(
        "generator_path", "generator_sha256", "copied_values_allowed",
        "additive_schema_successor", "production_defaults_changed"
    ), "workflow policy generation contract")
    lapply(envelope$receipts, workflowPolicyValidateReceipt)
    capability_ids <- vapply(
        envelope$receipts,
        `[[`,
        character(1),
        "capability_id"
    )
    basis <- envelope
    basis$envelope_digest <- NULL
    valid <- identical(
        envelope$schema,
        "multischolar.omics_auto_policy_receipt_envelope"
    ) && identical(envelope$schema_version, "1.0.0") &&
        identical(envelope$owner_ticket_id, "OMICS-ART-069") &&
        workflowPolicyScalarString(envelope$envelope_id) &&
        identical(
            envelope$consumer_contract_digest,
            workflowPolicyConsumerContractDigest()
        ) && length(capability_ids) > 0L && !anyDuplicated(capability_ids) &&
        workflowPolicyScalarString(envelope$generation_contract$generator_path) &&
        workflowPolicyDigestValid(
            envelope$generation_contract$generator_sha256
        ) && !isTRUE(
            envelope$generation_contract$copied_values_allowed
        ) && isTRUE(
            envelope$generation_contract$additive_schema_successor
        ) && !isTRUE(
            envelope$generation_contract$production_defaults_changed
        ) && identical(
            envelope$envelope_digest,
            workflowPolicyObjectDigest(basis)
        )
    if (!valid) {
        workflowPolicyAbort(
            "workflow policy envelope identity or digest differs",
            "multischolar_invalid_workflow_policy_envelope"
        )
    }
    invisible(envelope)
}

workflowPolicyInstalledEnvelopePath <- function() {
    override <- getOption("MultiScholaR.workflow_policy_receipt_path")
    if (workflowPolicyScalarString(override)) return(override)
    installed <- system.file(
        "extdata",
        "omics-auto-policy-receipts-v2.json",
        package = "MultiScholaR"
    )
    if (nzchar(installed)) return(installed)
    source_path <- file.path(
        "inst",
        "extdata",
        "omics-auto-policy-receipts-v2.json"
    )
    if (file.exists(source_path)) return(source_path)
    ""
}

workflowPolicyUnavailableEnvelope <- function(reason, condition = NULL) {
    structure(
        list(
            status = "unavailable",
            reason_code = reason,
            condition_class = if (is.null(condition)) NULL else class(condition)[[1L]]
        ),
        class = "WorkflowPolicyEnvelopeUnavailable"
    )
}

workflowPolicyLoadEnvelope <- function(
    path = workflowPolicyInstalledEnvelopePath(),
    injected = NULL,
    strict = FALSE
) {
    result <- tryCatch({
        envelope <- if (is.null(injected)) {
            if (!workflowPolicyScalarString(path) || !file.exists(path)) {
                workflowPolicyAbort(
                    "workflow policy envelope is missing",
                    "multischolar_missing_workflow_policy_envelope"
                )
            }
            jsonlite::read_json(path, simplifyVector = FALSE)
        } else {
            injected
        }
        workflowPolicyValidateEnvelope(envelope)
        envelope
    }, error = function(error) error)
    if (!inherits(result, "error")) return(result)
    if (isTRUE(strict)) stop(result)
    workflowPolicyUnavailableEnvelope(
        "policy_envelope_invalid_or_missing",
        result
    )
}

workflowPolicyMemoryDecision <- function(reason, capability_id = NULL) {
    list(
        capability_id = capability_id,
        effective_backend = "memory",
        effective_rollout = "none",
        reason_code = reason,
        receipt_id = NULL,
        receipt_digest = NULL,
        publication_authority = FALSE
    )
}

workflowPolicyResolve <- function(
    envelope,
    capability_id,
    measure_id,
    measure_bytes,
    platform = Sys.info()[["sysname"]],
    descriptor_digest = NULL
) {
    if (inherits(envelope, "WorkflowPolicyEnvelopeUnavailable")) {
        return(workflowPolicyMemoryDecision(
            envelope$reason_code,
            capability_id
        ))
    }
    valid <- tryCatch({
        workflowPolicyValidateEnvelope(envelope)
        TRUE
    }, error = function(...) FALSE)
    if (!valid) {
        return(workflowPolicyMemoryDecision(
            "policy_envelope_invalid_or_missing",
            capability_id
        ))
    }
    matches <- Filter(function(receipt) {
        identical(receipt$capability_id, capability_id)
    }, envelope$receipts)
    if (length(matches) != 1L) {
        return(workflowPolicyMemoryDecision(
            "policy_capability_unavailable",
            capability_id
        ))
    }
    receipt <- matches[[1L]]
    platforms <- unlist(receipt$platforms, use.names = FALSE)
    if (!platform %in% platforms) {
        return(workflowPolicyMemoryDecision(
            "policy_platform_unavailable",
            capability_id
        ))
    }
    if (!identical(receipt$size_measure$measure_id, measure_id) ||
        !is.numeric(measure_bytes) || length(measure_bytes) != 1L ||
        !is.finite(measure_bytes) || measure_bytes < 0) {
        return(workflowPolicyMemoryDecision(
            "policy_exact_measure_unavailable",
            capability_id
        ))
    }
    if (workflowPolicyDigestValid(descriptor_digest) && !identical(
        descriptor_digest,
        receipt$descriptor_digest
    )) {
        return(workflowPolicyMemoryDecision(
            "policy_descriptor_stale",
            capability_id
        ))
    }
    branch <- if (measure_bytes < receipt$threshold_bytes) {
        receipt$below_threshold
    } else {
        receipt$at_or_above_threshold
    }
    list(
        capability_id = capability_id,
        effective_backend = branch$backend,
        effective_rollout = branch$rollout,
        reason_code = if (identical(branch$backend, "artifact")) {
            "policy_threshold_artifact"
        } else {
            "policy_threshold_memory"
        },
        receipt_id = receipt$receipt_id,
        receipt_digest = receipt$receipt_digest,
        publication_authority = isTRUE(
            receipt$claim_scope$publication_authority
        )
    )
}

workflowPolicyLegacyReceipts <- function(
    envelope = workflowPolicyLoadEnvelope()
) {
    if (inherits(envelope, "WorkflowPolicyEnvelopeUnavailable")) return(list())
    Filter(function(receipt) {
        identical(
            receipt$size_measure$measure_id,
            "legacy_r_object_size_x2.v1"
        )
    }, envelope$receipts)
}
