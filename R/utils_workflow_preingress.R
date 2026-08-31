# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

workflowPreIngressAbort <- function(message, class, ...) {
    rlang::abort(
        message,
        class = c(class, "multischolar_workflow_preingress_error"),
        ...
    )
}

newWorkflowPreIngressContract <- function(
    contract_id,
    capability_ids,
    member_resolver,
    identity_probe,
    receipt_provider = function(injected) {
        workflowPolicyLoadEnvelope(injected = injected)
    },
    archive_measure = NULL,
    max_probe_bytes = 1024^2,
    max_probe_lines = 1000L
) {
    valid <- workflowPolicyScalarString(contract_id) &&
        is.character(capability_ids) && length(capability_ids) > 0L &&
        !anyNA(capability_ids) && all(nzchar(capability_ids)) &&
        !anyDuplicated(capability_ids) && is.function(member_resolver) &&
        is.function(identity_probe) && is.function(receipt_provider) &&
        (is.null(archive_measure) || is.function(archive_measure)) &&
        is.numeric(max_probe_bytes) && length(max_probe_bytes) == 1L &&
        is.finite(max_probe_bytes) && max_probe_bytes > 0 &&
        is.numeric(max_probe_lines) && length(max_probe_lines) == 1L &&
        is.finite(max_probe_lines) && max_probe_lines > 0
    if (!valid) {
        workflowPreIngressAbort(
            "workflow pre-ingress contract is malformed",
            "multischolar_invalid_workflow_preingress_contract"
        )
    }
    structure(
        list(
            contract_id = contract_id,
            capability_ids = capability_ids,
            member_resolver = member_resolver,
            identity_probe = identity_probe,
            receipt_provider = receipt_provider,
            archive_measure = archive_measure,
            max_probe_bytes = as.numeric(max_probe_bytes),
            max_probe_lines = as.integer(max_probe_lines)
        ),
        class = "WorkflowPreIngressContract"
    )
}

validateWorkflowPreIngressContract <- function(contract) {
    expected <- c(
        "contract_id", "capability_ids", "member_resolver", "identity_probe",
        "receipt_provider", "archive_measure", "max_probe_bytes",
        "max_probe_lines"
    )
    valid <- inherits(contract, "WorkflowPreIngressContract") &&
        is.list(contract) && identical(names(contract), expected) &&
        workflowPolicyScalarString(contract$contract_id) &&
        is.character(contract$capability_ids) &&
        length(contract$capability_ids) > 0L &&
        !anyNA(contract$capability_ids) &&
        all(nzchar(contract$capability_ids)) &&
        !anyDuplicated(contract$capability_ids) &&
        is.function(contract$member_resolver) &&
        is.function(contract$identity_probe) &&
        is.function(contract$receipt_provider) &&
        (is.null(contract$archive_measure) ||
            is.function(contract$archive_measure)) &&
        is.numeric(contract$max_probe_bytes) &&
        length(contract$max_probe_bytes) == 1L &&
        is.finite(contract$max_probe_bytes) &&
        contract$max_probe_bytes > 0 &&
        is.numeric(contract$max_probe_lines) &&
        length(contract$max_probe_lines) == 1L &&
        is.finite(contract$max_probe_lines) &&
        contract$max_probe_lines > 0
    if (!valid) {
        workflowPreIngressAbort(
            "workflow pre-ingress contract is invalid",
            "multischolar_invalid_workflow_preingress_contract"
        )
    }
    contract
}

workflowPreIngressNormalizeRoot <- function(project_root) {
    if (!workflowPolicyScalarString(project_root)) {
        workflowPreIngressAbort(
            "workflow pre-ingress project root is invalid",
            "multischolar_invalid_workflow_preingress_root"
        )
    }
    normalizePath(project_root, winslash = "/", mustWork = TRUE)
}

workflowPreIngressMemberFields <- function() {
    c("member_id", "path", "container_type", "semantic_order")
}

workflowPreIngressMemberSnapshot <- function(member, project_root) {
    workflowPolicyRequireNames(
        member,
        workflowPreIngressMemberFields(),
        "workflow pre-ingress member"
    )
    if (!workflowPolicyScalarString(member$member_id) ||
        !workflowPolicyScalarString(member$path) ||
        !member$container_type %in% c("plain", "archive") ||
        !is.numeric(member$semantic_order) ||
        length(member$semantic_order) != 1L ||
        !is.finite(member$semantic_order)) {
        workflowPreIngressAbort(
            "workflow pre-ingress member is malformed",
            "multischolar_invalid_workflow_preingress_member"
        )
    }
    if (nzchar(Sys.readlink(member$path))) {
        workflowPreIngressAbort(
            "workflow pre-ingress member cannot be a symbolic link",
            "multischolar_unsafe_workflow_preingress_member"
        )
    }
    path <- normalizePath(member$path, winslash = "/", mustWork = TRUE)
    in_root <- identical(path, project_root) ||
        startsWith(path, paste0(project_root, "/"))
    info <- fs::file_info(path)
    if (!in_root || !identical(as.character(info$type[[1L]]), "file")) {
        workflowPreIngressAbort(
            "workflow pre-ingress member is outside its root or not regular",
            "multischolar_unsafe_workflow_preingress_member"
        )
    }
    list(
        member_id = member$member_id,
        canonical_path = path,
        container_type = member$container_type,
        semantic_order = as.integer(member$semantic_order),
        size_bytes = as.numeric(info$size[[1L]]),
        device_id = as.numeric(info$device_id[[1L]]),
        inode = as.numeric(info$inode[[1L]]),
        modification_time = as.numeric(info$modification_time[[1L]]),
        change_time = as.numeric(info$change_time[[1L]])
    )
}

workflowPreIngressCanonicalMembers <- function(members, project_root) {
    if (!is.list(members) || !length(members)) {
        workflowPreIngressAbort(
            "workflow pre-ingress requires at least one member",
            "multischolar_invalid_workflow_preingress_member"
        )
    }
    root <- workflowPreIngressNormalizeRoot(project_root)
    snapshots <- lapply(members, workflowPreIngressMemberSnapshot, root)
    ids <- vapply(snapshots, `[[`, character(1), "member_id")
    paths <- vapply(snapshots, `[[`, character(1), "canonical_path")
    order_values <- vapply(snapshots, `[[`, integer(1), "semantic_order")
    physical <- vapply(snapshots, function(member) {
        paste(member$device_id, member$inode, sep = ":")
    }, character(1))
    valid <- !anyDuplicated(ids) && !anyDuplicated(paths) &&
        !anyDuplicated(physical) && !anyDuplicated(order_values) &&
        identical(sort(order_values), seq_along(order_values))
    if (!valid) {
        workflowPreIngressAbort(
            "workflow pre-ingress members are duplicated or unordered",
            "multischolar_duplicate_workflow_preingress_member"
        )
    }
    snapshots[order(order_values, method = "radix")]
}

workflowPreIngressValidateProbe <- function(probe, contract) {
    workflowPolicyRequireNames(probe, c(
        "identity", "capability_id", "capability_version", "descriptor_digest",
        "bytes_read", "lines_read", "schema_valid", "ambiguous",
        "complete_payload_materialized"
    ), "workflow pre-ingress probe")
    identity_valid <- is.list(probe$identity) &&
        all(.WORKFLOW_IDENTITY_FIELDS %in% names(probe$identity)) &&
        all(vapply(
            probe$identity[.WORKFLOW_IDENTITY_FIELDS],
            workflowPolicyScalarString,
            logical(1)
        ))
    logical_flag <- function(value) {
        is.logical(value) && length(value) == 1L && !is.na(value)
    }
    valid <- identity_valid &&
        workflowPolicyScalarString(probe$capability_id) &&
        workflowPolicyScalarString(probe$capability_version) &&
        workflowPolicyDigestValid(probe$descriptor_digest) &&
        is.numeric(probe$bytes_read) && length(probe$bytes_read) == 1L &&
        is.finite(probe$bytes_read) && probe$bytes_read >= 0 &&
        probe$bytes_read <= contract$max_probe_bytes &&
        is.numeric(probe$lines_read) && length(probe$lines_read) == 1L &&
        is.finite(probe$lines_read) && probe$lines_read >= 0 &&
        probe$lines_read <= contract$max_probe_lines &&
        logical_flag(probe$schema_valid) && logical_flag(probe$ambiguous) &&
        logical_flag(probe$complete_payload_materialized) &&
        !isTRUE(probe$complete_payload_materialized)
    if (!valid) {
        workflowPreIngressAbort(
            "workflow pre-ingress probe is invalid or unbounded",
            "multischolar_invalid_workflow_preingress_probe"
        )
    }
    catalogue <- workflowCapabilityCatalogue()
    capability <- catalogue[[probe$capability_id]]
    if (is.null(capability) ||
        !probe$capability_id %in% contract$capability_ids) {
        return(list(
            available = FALSE,
            reason_code = "unknown_workflow_format"
        ))
    }
    capability_fields <- intersect(
        names(capability$identity),
        names(probe$identity)
    )
    if (!identical(
        probe$identity[capability_fields],
        capability$identity[capability_fields]
    )) {
        return(list(
            available = FALSE,
            reason_code = "workflow_identity_mismatch"
        ))
    }
    if (!isTRUE(probe$schema_valid)) {
        return(list(
            available = FALSE,
            reason_code = "unsupported_workflow_schema"
        ))
    }
    if (isTRUE(probe$ambiguous)) {
        return(list(
            available = FALSE,
            reason_code = "ambiguous_workflow_format"
        ))
    }
    list(available = TRUE, reason_code = "workflow_format_identified")
}

workflowPreIngressMeasure <- function(members, contract) {
    archive <- vapply(members, function(member) {
        identical(member$container_type, "archive")
    }, logical(1))
    if (any(archive) && is.null(contract$archive_measure)) {
        return(list(
            measure_id = "total_uncompressed_input_bytes_v1",
            exact = FALSE,
            available = FALSE,
            bytes = NULL,
            reason_code = "exact_archive_size_unavailable"
        ))
    }
    if (any(archive)) {
        measured <- contract$archive_measure(members)
        valid <- is.list(measured) &&
            identical(names(measured), c("bytes", "exact", "bounded")) &&
            is.numeric(measured$bytes) && length(measured$bytes) == 1L &&
            is.finite(measured$bytes) && measured$bytes >= 0 &&
            isTRUE(measured$exact) && isTRUE(measured$bounded)
        if (!valid) {
            return(list(
                measure_id = "total_uncompressed_input_bytes_v1",
                exact = FALSE,
                available = FALSE,
                bytes = NULL,
                reason_code = "exact_archive_size_unavailable"
            ))
        }
        bytes <- as.numeric(measured$bytes)
    } else {
        bytes <- sum(vapply(members, `[[`, numeric(1), "size_bytes"))
    }
    list(
        measure_id = "total_uncompressed_input_bytes_v1",
        exact = TRUE,
        available = TRUE,
        bytes = bytes,
        reason_code = "exact_input_bytes_available"
    )
}

workflowPreIngressMemoryOutcome <- function(reason, contract_id) {
    structure(
        list(
            status = "memory",
            contract_id = contract_id,
            effective_backend = "memory",
            effective_rollout = "none",
            reason_code = reason,
            token = NULL,
            importer_invoked = FALSE,
            artifact_resources_opened = FALSE
        ),
        class = "WorkflowPreIngressOutcome"
    )
}

workflowPreIngressToken <- function(
    contract,
    project_root,
    members,
    probe,
    measure,
    decision
) {
    member_digest <- workflowPolicyObjectDigest(members)
    receipt_binding <- list(
        receipt_id = decision$receipt_id,
        receipt_digest = decision$receipt_digest
    )
    basis <- list(
        schema = "multischolar.workflow_preingress_token",
        schema_version = "1.0.0",
        contract_id = contract$contract_id,
        state = "bound",
        creator_pid = as.integer(Sys.getpid()),
        project_root = project_root,
        identity = probe$identity,
        capability_id = probe$capability_id,
        capability_version = probe$capability_version,
        descriptor_digest = probe$descriptor_digest,
        members = members,
        member_digest = member_digest,
        measure = measure,
        decision = decision,
        receipt_binding = receipt_binding,
        probe_evidence = list(
            bytes_read = as.numeric(probe$bytes_read),
            lines_read = as.integer(probe$lines_read),
            complete_payload_materialized = FALSE
        ),
        publication_authority = FALSE
    )
    basis$token_id <- workflowPreIngressTokenDigest(basis)
    structure(basis, class = "WorkflowPreIngressToken")
}

workflowPreIngressTokenFields <- function() {
    c(
        "schema", "schema_version", "contract_id", "state", "creator_pid",
        "project_root", "identity", "capability_id", "capability_version",
        "descriptor_digest", "members", "member_digest", "measure",
        "decision", "receipt_binding", "probe_evidence",
        "publication_authority", "token_id"
    )
}

workflowPreIngressTokenDigest <- function(token) {
    workflowPolicyObjectDigest(token[setdiff(names(token), "token_id")])
}

workflowPreIngressValidateTokenIdentity <- function(token, state = "bound") {
    valid <- inherits(token, "WorkflowPreIngressToken") &&
        is.list(token) && setequal(
            names(token),
            workflowPreIngressTokenFields()
        ) && identical(token$state, state) &&
        identical(token$creator_pid, as.integer(Sys.getpid())) &&
        workflowPolicyDigestValid(token$token_id) &&
        identical(token$token_id, workflowPreIngressTokenDigest(token))
    if (!valid) {
        workflowPreIngressAbort(
            "workflow pre-ingress token is invalid or process-stale",
            "multischolar_invalid_workflow_preingress_token"
        )
    }
    invisible(token)
}

workflowPreIngressResolve <- function(
    contract,
    input,
    project_root,
    requested_backend = "auto",
    injected_envelope = NULL
) {
    contract <- validateWorkflowPreIngressContract(contract)
    if (!requested_backend %in% c("auto", "memory")) {
        workflowPreIngressAbort(
            "workflow pre-ingress supports auto or memory routing",
            "multischolar_invalid_workflow_preingress_backend"
        )
    }
    if (identical(requested_backend, "memory")) {
        return(workflowPreIngressMemoryOutcome(
            "explicit_memory_without_preingress_probe",
            contract$contract_id
        ))
    }
    root <- workflowPreIngressNormalizeRoot(project_root)
    members <- workflowPreIngressCanonicalMembers(
        contract$member_resolver(input),
        root
    )
    probe <- contract$identity_probe(members)
    probe_resolution <- workflowPreIngressValidateProbe(probe, contract)
    if (!isTRUE(probe_resolution$available)) {
        return(workflowPreIngressMemoryOutcome(
            probe_resolution$reason_code,
            contract$contract_id
        ))
    }
    measure <- workflowPreIngressMeasure(members, contract)
    if (!isTRUE(measure$available) || !isTRUE(measure$exact)) {
        return(workflowPreIngressMemoryOutcome(
            measure$reason_code,
            contract$contract_id
        ))
    }
    envelope <- contract$receipt_provider(injected_envelope)
    decision <- workflowPolicyResolve(
        envelope,
        capability_id = probe$capability_id,
        measure_id = measure$measure_id,
        measure_bytes = measure$bytes,
        descriptor_digest = probe$descriptor_digest
    )
    token <- workflowPreIngressToken(
        contract,
        root,
        members,
        probe,
        measure,
        decision
    )
    structure(
        list(
            status = "bound",
            contract_id = contract$contract_id,
            effective_backend = decision$effective_backend,
            effective_rollout = decision$effective_rollout,
            reason_code = decision$reason_code,
            token = token,
            importer_invoked = FALSE,
            artifact_resources_opened = FALSE
        ),
        class = "WorkflowPreIngressOutcome"
    )
}

validateWorkflowPreIngressToken <- function(token, contract) {
    contract <- validateWorkflowPreIngressContract(contract)
    workflowPreIngressValidateTokenIdentity(token)
    if (!identical(token$contract_id, contract$contract_id) ||
        !identical(token$member_digest, workflowPolicyObjectDigest(token$members))) {
        workflowPreIngressAbort(
            "workflow pre-ingress token is invalid or process-stale",
            "multischolar_invalid_workflow_preingress_token"
        )
    }
    current <- workflowPreIngressCanonicalMembers(
        lapply(token$members, function(member) {
            list(
                member_id = member$member_id,
                path = member$canonical_path,
                container_type = member$container_type,
                semantic_order = member$semantic_order
            )
        }),
        token$project_root
    )
    if (!identical(
        workflowPolicyObjectDigest(current),
        token$member_digest
    )) {
        workflowPreIngressAbort(
            "workflow pre-ingress members changed before import",
            "multischolar_mutated_workflow_preingress_member"
        )
    }
    invisible(token)
}

workflowPreIngressConsume <- function(contract, outcome, input, importer) {
    contract <- validateWorkflowPreIngressContract(contract)
    if (!is.function(importer)) {
        workflowPreIngressAbort(
            "workflow pre-ingress importer must be a function",
            "multischolar_invalid_workflow_preingress_importer"
        )
    }
    if (!inherits(outcome, "WorkflowPreIngressOutcome") ||
        !identical(outcome$status, "bound") || is.null(outcome$token)) {
        return(list(
            outcome = outcome,
            imported = FALSE,
            result = NULL
        ))
    }
    validateWorkflowPreIngressToken(outcome$token, contract)
    result <- importer(input, outcome$token)
    consumed <- outcome$token
    consumed$state <- "consumed"
    consumed$token_id <- workflowPreIngressTokenDigest(consumed)
    outcome$token <- structure(consumed, class = "WorkflowPreIngressToken")
    outcome$status <- "consumed"
    outcome$importer_invoked <- TRUE
    list(outcome = outcome, imported = TRUE, result = result)
}

workflowPreIngressBindContext <- function(context, outcome, path_builder) {
    if (!inherits(context, "WorkflowContext") ||
        !inherits(outcome, "WorkflowPreIngressOutcome") ||
        !identical(outcome$status, "bound") ||
        !is.function(path_builder)) {
        workflowPreIngressAbort(
            "workflow pre-ingress context binding is invalid",
            "multischolar_invalid_workflow_preingress_context"
        )
    }
    token <- outcome$token
    workflowPreIngressValidateTokenIdentity(token)
    resolution <- list(
        requested_backend = "auto",
        effective_backend = outcome$effective_backend,
        effective_rollout = outcome$effective_rollout,
        capability_id = token$capability_id,
        capability_version = token$capability_version,
        reason_code = outcome$reason_code,
        project_state = "new",
        receipt_id = token$receipt_binding$receipt_id,
        receipt_digest = token$receipt_binding$receipt_digest,
        preingress_token_id = token$token_id
    )
    paths <- if (identical(outcome$effective_backend, "artifact")) {
        path_builder(context$getProjectRoot(), token$identity)
    } else {
        NULL
    }
    context$bind(token$identity, resolution, paths)
    invisible(context$getSnapshot())
}
