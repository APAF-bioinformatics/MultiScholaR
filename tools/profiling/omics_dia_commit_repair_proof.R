diaRepairScientificSlotDigests <- function(object, salt) {
    slot_names <- methods::slotNames(object)
    values <- vapply(slot_names, function(slot_name) {
        diaRepairSaltedDigest(
            diaRepairStableDigest(methods::slot(object, slot_name)),
            salt
        )
    }, character(1))
    as.list(values)
}

diaRepairDesignProof <- function(object, salt) {
    design <- object@design_matrix
    column_digests <- vapply(design, function(value) {
        diaRepairSaltedDigest(diaRepairStableDigest(value), salt)
    }, character(1))
    list(
        rows = as.integer(nrow(design)),
        columns = as.integer(ncol(design)),
        column_names = as.list(names(design)),
        column_classes = as.list(vapply(
            design,
            function(value) class(value)[[1L]],
            character(1)
        )),
        salted_column_digests = as.list(column_digests),
        salted_attributes_digest = diaRepairSaltedDigest(
            diaRepairStableDigest(attributes(design)),
            salt
        )
    )
}

diaRepairCompletePairRecords <- function(records) {
    groups <- split(records, vapply(records, `[[`, character(1), "pair_id"))
    complete <- Filter(function(pair) {
        arms <- vapply(pair, `[[`, character(1), "arm")
        proofs <- lapply(pair, function(record) {
            record$measurement$worker$scientific_proof
        })
        length(pair) == 2L && setequal(
            arms,
            c("pre_repair_artifact", "candidate_artifact")
        ) && all(vapply(pair, function(record) {
            measurement <- record$measurement
            identical(measurement$status, "passed") &&
                isTRUE(measurement$publication_certifiable) &&
                isTRUE(measurement$cleanup$valid)
        }, logical(1))) && all(vapply(proofs, is.list, logical(1))) &&
            identical(proofs[[1L]], proofs[[2L]])
    }, groups)
    unlist(complete, recursive = FALSE, use.names = FALSE)
}

diaRepairProofWorker <- function(package_library, run_dir, arm, salt) {
    .libPaths(c(package_library, .libPaths()))
    namespace <- loadNamespace("MultiScholaR", lib.loc = package_library)
    paths <- diaRepairWorkflowPaths(run_dir)
    rollout <- if (identical(arm, "candidate_artifact")) "evict" else "dual_write"
    storage_policy <- list(
        requested_backend = "artifact",
        requested_rollout = rollout,
        migration_requested = TRUE,
        project_id = paths$project_id
    )
    prepared <- diaRepairNamespaceValue(
        namespace,
        "createProtDiaResumeContext"
    )(paths, "dia-repair-canary", storage_policy)
    if (!isTRUE(prepared$enabled)) {
        diaRepairAbort("DIA repair proof context is unavailable")
    }
    manager <- diaRepairNamespaceValue(
        namespace,
        "newProtDiaArtifactStateManager"
    )(prepared$context)
    on.exit(manager$close(), add = TRUE)
    object <- manager$getState()
    list(
        salted_state_digest = diaRepairSaltedDigest(
            diaRepairStableDigest(object),
            salt
        ),
        salted_slot_digests = diaRepairScientificSlotDigests(object, salt),
        design_proof = diaRepairDesignProof(object, salt),
        s4_class = class(object)[[1L]],
        valid_s4 = identical(methods::validObject(object, test = TRUE), TRUE),
        current_state = manager$getCurrentStateName()
    )
}

diaRepairScientificProof <- function(
    workflow,
    commit,
    salt,
    package_library,
    run_dir,
    arm
) {
    verification <- commit$design$hydration_verification
    process_parity <- if (is.null(verification)) {
        TRUE
    } else {
        isTRUE(verification$valid) &&
            identical(
                verification$expected_digest,
                verification$hydrated_digest
            ) && !isTRUE(verification$complete_payload_returned)
    }
    proof_path <- file.path(run_dir, "scientific-proof.json")
    proof_process <- processx::run(
        file.path(R.home("bin"), "Rscript"),
        c(
            "--vanilla", diaRepairRunnerPath(),
            "--mode", "proof",
            "--arm", arm,
            "--run-dir", run_dir,
            "--package-library", package_library
        ),
        error_on_status = FALSE,
        echo = FALSE
    )
    if (proof_process$status != 0L || !file.exists(proof_path)) {
        diaRepairAbort("DIA repair scientific proof worker failed")
    }
    proof <- publicationReadJson(proof_path)
    if (!isTRUE(proof$valid_s4) || !isTRUE(process_parity)) {
        diaRepairAbort("DIA repair benchmark produced an invalid S4 state")
    }
    c(proof, list(
        source_fields_released = is.null(workflow$data_tbl) &&
            is.null(workflow$data_cln)
    ))
}
