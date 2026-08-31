publicationComparatorEnvelopeRoot <- function(root = publicationComparatorRoot()) {
    file.path(root, "comparator-envelopes")
}

publicationComparatorEntry <- function(authority, comparator_id) {
    matches <- vapply(
        authority$comparators,
        \(entry) identical(entry$comparator_id, comparator_id),
        logical(1)
    )
    if (sum(matches) != 1L) {
        publicationComparatorAbort("Comparator envelope id is not unique")
    }
    authority$comparators[[which(matches)]]
}

publicationComparatorEnvelopeDecision <- function(
    comparator_id,
    backend,
    authority
) {
    publicationValidateComparatorAuthority(authority)
    entry <- publicationComparatorEntry(authority, comparator_id)
    if (!backend %in% entry$backend_envelopes) {
        return(list(
            status = "rejected",
            reason = "backend_not_declared_for_comparator",
            runner_invoked = FALSE
        ))
    }
    if (identical(entry$status, "pending_omics_art_077")) {
        return(list(
            status = "rejected",
            reason = "candidate_pending_omics_art_077",
            runner_invoked = FALSE
        ))
    }
    if (!identical(entry$status, "verified")) {
        return(list(
            status = "rejected",
            reason = "comparator_build_not_verified",
            runner_invoked = FALSE
        ))
    }
    list(status = "allowed", reason = NULL, runner_invoked = TRUE)
}

publicationComparatorEnvelopePaths <- function(
    comparator_id,
    backend,
    attempt = 1L,
    root = publicationComparatorRoot()
) {
    envelope_root <- file.path(
        publicationComparatorEnvelopeRoot(root),
        comparator_id,
        backend,
        paste0("attempt-", attempt)
    )
    list(
        root = envelope_root,
        project = file.path(envelope_root, "project"),
        script = file.path(envelope_root, "tiny-workload.R"),
        output = file.path(envelope_root, "result.json"),
        logs = file.path(envelope_root, "logs")
    )
}

publicationComparatorLibraryPaths <- function(build_receipt_path) {
    evidence <- publicationValidateComparatorBuildReceiptEvidence(
        build_receipt_path
    )
    receipt <- evidence$receipt
    package_library <- dirname(receipt$build$installed_inventory$root)
    restore_receipt <- receipt$environment$restore_receipt$path
    dependency_library <- file.path(dirname(restore_receipt), "library")
    if (!dir.exists(package_library) || !dir.exists(dependency_library)) {
        publicationBuildAbort("Comparator envelope library is missing")
    }
    list(
        package_library = package_library,
        dependency_library = dependency_library
    )
}

publicationComparatorMemoryWorkloadLines <- function() {
    c(
        "manager <- MultiScholaR::WorkflowState$new(audit_enabled = FALSE)",
        "payload <- matrix(c(1, 4, 16, 64), nrow = 2L)",
        paste0(
            "manager$saveState(\"tiny\", payload, list(stage = \"tiny\"), ",
            "\"Tiny memory state\")"
        ),
        "stopifnot(identical(manager$getState(), payload))",
        paste0(
            "state_name <- if (is.function(manager$getCurrentStateName)) ",
            "manager$getCurrentStateName() else manager$current_state"
        ),
        "stopifnot(identical(state_name, \"tiny\"))",
        "manager_class <- class(manager)[[1L]]",
        "artifact_file_count <- 0L",
        "settled_payload_free <- FALSE",
        "if (is.function(manager$close)) manager$close()"
    )
}

publicationComparatorArtifactWorkloadLines <- function() {
    c(
        "ns <- asNamespace(\"MultiScholaR\")",
        paste0(
            "capabilities <- Filter(\\(capability) ",
            "identical(capability$identity$input_format, \"diann\") && ",
            "identical(capability$identity$data_level, \"peptide\"), ",
            "ns$workflowCapabilityCatalogue())"
        ),
        "stopifnot(length(capabilities) == 1L)",
        "capability <- capabilities[[1L]]",
        "capability$artifact_eligible <- TRUE",
        "capability$auto_eligible <- TRUE",
        "capability$maximum_artifact_rollout <- \"dual_write\"",
        paste0(
            "context <- ns$createWorkflowContext(",
            "experiment_paths = list(base_dir = project_root, ",
            "omic_label = \"proteomics_study\"), omic_type = \"proteomics\", ",
            "storage_policy = list(requested_backend = \"artifact\", ",
            "requested_rollout = \"dual_write\", migration_requested = TRUE, ",
            "project_id = \"comparator-tiny-artifact\"))"
        ),
        paste0(
            "ns$bindWorkflowContextFromImport(context, workflow_type = \"DIA\", ",
            "input_format = \"diann\", data_level = \"peptide\", ",
            "capabilities = list(capability))"
        ),
        "manager <- ns$newWorkflowState(workflow_context = context)",
        "manager$setWorkflowType(\"DIA\")",
        paste0(
            "peptide_data <- data.frame(Protein.Group = rep(c(\"PG1\", ",
            "\"PG2\"), each = 2L), Protein.Ids = rep(c(\"P1\", \"P2\"), ",
            "each = 2L), Stripped.Sequence = rep(c(\"PEPTIDE_A\", ",
            "\"PEPTIDE_B\"), each = 2L), Run = rep(c(\"S1\", \"S2\"), 2L), ",
            "Q.Value = c(0.001, 0.002, 0.003, 0.004), ",
            "Global.Q.Value = c(0.005, 0.006, 0.007, 0.008), ",
            "Proteotypic = c(1L, 1L, 0L, 0L), ",
            "Precursor.Quantity = c(100, 200, 300, 400), ",
            "Precursor.Normalised = c(10, 20, 30, 40), ",
            "stringsAsFactors = FALSE)"
        ),
        paste0(
            "design <- data.frame(Run = c(\"S1\", \"S2\"), ",
            "group = factor(c(\"control\", \"case\"), ",
            "levels = c(\"case\", \"control\")), replicates = c(\"R1\", ",
            "\"R1\"), stringsAsFactors = FALSE)"
        ),
        paste0(
            "object <- MultiScholaR::PeptideQuantitativeDataDiann(",
            "peptide_data, design, sample_id = \"Run\", group_id = \"group\", ",
            "technical_replicate_id = \"replicates\", ",
            "args = list(branch = \"comparator-tiny\"))"
        ),
        paste0(
            "object@peptide_matrix <- matrix(c(10, 30, 20, 40), nrow = 2L, ",
            "dimnames = list(c(\"PG1%PEPTIDE_A\", \"PG2%PEPTIDE_B\"), ",
            "c(\"S1\", \"S2\")))"
        ),
        "object <- ns$.ensurePeptideFeatureKeyMap(object, record_migration = FALSE)",
        paste0(
            "manager$saveState(\"imported\", object, list(stage = \"import\"), ",
            "\"Tiny artifact state\")"
        ),
        "restored <- manager$getState()",
        "stopifnot(identical(class(restored), class(object)))",
        paste0(
            "stopifnot(all(vapply(methods::slotNames(object), ",
            "\\(slot) identical(methods::slot(restored, slot), ",
            "methods::slot(object, slot)), logical(1))))"
        ),
        "settled_payload_free <- !any(vapply(manager$states, isS4, logical(1)))",
        "stopifnot(settled_payload_free)",
        "manager_class <- class(manager)[[1L]]",
        "state_name <- manager$getCurrentStateName()",
        "manager$close()",
        paste0(
            "artifact_file_count <- length(list.files(project_root, ",
            "recursive = TRUE, all.files = TRUE, no.. = TRUE))"
        ),
        "stopifnot(artifact_file_count > 0L)"
    )
}

publicationComparatorEnvelopeScript <- function(
    comparator_id,
    backend,
    libraries,
    paths
) {
    workload <- if (identical(backend, "artifact")) {
        publicationComparatorArtifactWorkloadLines()
    } else {
        publicationComparatorMemoryWorkloadLines()
    }
    paste(c(
        paste0(
            "package_library <- ",
            publicationBuildRLiteral(libraries$package_library)
        ),
        paste0(
            "dependency_library <- ",
            publicationBuildRLiteral(libraries$dependency_library)
        ),
        paste0("project_root <- ", publicationBuildRLiteral(paths$project)),
        paste0("output <- ", publicationBuildRLiteral(paths$output)),
        ".libPaths(c(package_library, dependency_library, .Library))",
        "dir.create(project_root, recursive = TRUE, showWarnings = FALSE)",
        "invisible(loadNamespace(\"MultiScholaR\"))",
        workload,
        "jsonlite::write_json(list(",
        "    status = \"passed\",",
        paste0("    comparator_id = ", publicationBuildRLiteral(comparator_id), ","),
        paste0("    backend = ", publicationBuildRLiteral(backend), ","),
        "    manager_class = manager_class,",
        "    state_name = state_name,",
        "    artifact_file_count = artifact_file_count,",
        "    settled_payload_free = settled_payload_free",
        "), output, auto_unbox = TRUE, pretty = TRUE, null = \"null\")"
    ), collapse = "\n")
}

publicationRunComparatorTinyEnvelope <- function(
    comparator_id,
    backend,
    build_receipt_path = NULL,
    attempt = 1L,
    root = publicationComparatorRoot()
) {
    authority <- publicationReadJson(
        "tests/testdata/omics-performance/comparators-v1.json"
    )
    decision <- publicationComparatorEnvelopeDecision(
        comparator_id,
        backend,
        authority
    )
    if (!identical(decision$status, "allowed")) {
        return(c(
            list(
                schema = "multischolar.omics_publication_comparator_envelope",
                schema_version = "1.0.0",
                comparator_id = comparator_id,
                backend = backend
            ),
            decision,
            list(publication_authority = FALSE)
        ))
    }
    if (is.null(build_receipt_path)) {
        publicationBuildAbort("Comparator envelope build receipt is required")
    }
    paths <- publicationComparatorEnvelopePaths(
        comparator_id,
        backend,
        attempt,
        root
    )
    if (file.exists(paths$root) || dir.exists(paths$root)) {
        publicationBuildAbort("Comparator envelope root already exists")
    }
    dir.create(paths$root, recursive = TRUE, showWarnings = FALSE)
    libraries <- publicationComparatorLibraryPaths(build_receipt_path)
    script <- publicationBuildWriteScript(
        publicationComparatorEnvelopeScript(
            comparator_id,
            backend,
            libraries,
            paths
        ),
        paths$script
    )
    library_path <- paste(
        c(libraries$package_library, libraries$dependency_library),
        collapse = .Platform$path.sep
    )
    environment <- publicationBuildEnvironment(
        library_path,
        file.path(paths$root, "empty-site-library")
    )
    environment[["HOME"]] <- file.path(paths$root, "home")
    environment[["TMPDIR"]] <- file.path(paths$root, "tmp")
    command <- publicationBuildRun(
        file.path(R.home("bin"), "Rscript"),
        c("--vanilla", paths$script),
        paths$root,
        paths$logs,
        environment,
        timeout_seconds = 600
    )
    publicationBuildRequireSuccess(command, "Comparator tiny envelope")
    evidence <- publicationReadJson(paths$output)
    list(
        schema = "multischolar.omics_publication_comparator_envelope",
        schema_version = "1.0.0",
        comparator_id = comparator_id,
        backend = backend,
        status = "passed",
        runner_invoked = TRUE,
        build_receipt = list(
            path = build_receipt_path,
            sha256 = publicationFileDigest(build_receipt_path)
        ),
        script = script,
        command = command,
        output = list(
            path = paths$output,
            sha256 = publicationFileDigest(paths$output)
        ),
        evidence = evidence,
        publication_authority = FALSE
    )
}

publicationComparatorEnvelopeKey <- function(comparator_id, backend) {
    paste(comparator_id, backend, sep = "::")
}

publicationComparatorEnvelopeExpectedKeys <- function() {
    c(
        "historical_janitor::memory",
        "historical_janitor::artifact",
        "pre_repair_performance::memory",
        "pre_repair_performance::artifact",
        "candidate::memory",
        "candidate::artifact",
        "candidate::proposed_auto"
    )
}

publicationValidateComparatorEnvelopeRecord <- function(record, authority) {
    common <- c(
        "schema", "schema_version", "comparator_id", "backend", "status",
        "runner_invoked", "publication_authority"
    )
    decision <- publicationComparatorEnvelopeDecision(
        record$comparator_id,
        record$backend,
        authority
    )
    if (identical(record$status, "rejected")) {
        publicationRequireNames(record, c(common, "reason"), "Envelope rejection")
        valid <- identical(record$reason, decision$reason) &&
            identical(decision$status, "rejected") &&
            !isTRUE(record$runner_invoked)
    } else {
        publicationRequireNames(record, c(
            common, "build_receipt", "script", "command", "output", "evidence"
        ), "Envelope evidence")
        root <- normalizePath(record$command$workdir, mustWork = TRUE)
        publicationEvidenceFileCurrent(
            record$build_receipt$path,
            record$build_receipt$sha256,
            "envelope build receipt"
        )
        publicationValidateScriptEvidence(record$script, root, "envelope")
        publicationValidateCommandEvidence(record$command, root, "envelope")
        publicationEvidenceFileCurrent(
            record$output$path,
            record$output$sha256,
            "envelope output",
            root
        )
        valid <- identical(decision$status, "allowed") &&
            identical(record$status, "passed") &&
            isTRUE(record$runner_invoked) &&
            identical(publicationReadJson(record$output$path), record$evidence)
    }
    valid <- valid && identical(
        record$schema,
        "multischolar.omics_publication_comparator_envelope"
    ) && identical(record$schema_version, "1.0.0") &&
        !isTRUE(record$publication_authority)
    if (!valid) publicationComparatorAbort("Comparator envelope evidence differs")
    invisible(record)
}

publicationCreateComparatorEnvelopeSummary <- function(entries, authority) {
    keys <- vapply(entries, \(entry) {
        publicationComparatorEnvelopeKey(entry$comparator_id, entry$backend)
    }, character(1))
    if (!identical(keys, publicationComparatorEnvelopeExpectedKeys())) {
        publicationComparatorAbort("Comparator envelope set differs")
    }
    lapply(entries, publicationValidateComparatorEnvelopeRecord, authority)
    list(
        schema = "multischolar.omics_publication_comparator_envelope_summary",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-064",
        comparator_authority = list(
            path = "tests/testdata/omics-performance/comparators-v1.json",
            sha256 = publicationFileDigest(
                "tests/testdata/omics-performance/comparators-v1.json"
            )
        ),
        entries = entries,
        publication_authority = FALSE
    )
}

publicationValidateComparatorEnvelopeSummary <- function(record) {
    publicationRequireNames(record, c(
        "schema", "schema_version", "owner_ticket_id", "comparator_authority",
        "entries", "publication_authority"
    ), "Comparator envelope summary")
    authority <- publicationReadJson(record$comparator_authority$path)
    publicationValidateComparatorAuthority(authority)
    current <- publicationCreateComparatorEnvelopeSummary(
        record$entries,
        authority
    )
    if (!identical(record, current)) {
        publicationComparatorAbort("Comparator envelope summary differs")
    }
    invisible(record)
}
