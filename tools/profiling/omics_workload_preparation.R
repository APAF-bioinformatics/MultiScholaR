omicsMemoryContainedPath <- function(path, root) {
    resolved_path <- normalizePath(path, mustWork = TRUE)
    resolved_root <- normalizePath(root, mustWork = TRUE)
    identical(resolved_path, resolved_root) || startsWith(
        resolved_path,
        paste0(resolved_root, .Platform$file.sep)
    )
}

omicsMemoryPrepareWorkload <- function(adapter, contract, run_dir, repetition) {
    rng_kind <- unlist(contract$rng$kind, use.names = FALSE)
    do.call(RNGkind, as.list(rng_kind))
    set.seed(as.integer(contract$rng$seed))
    prepared <- adapter$prepare(list(
        contract = contract,
        run_dir = run_dir,
        repetition = repetition
    ))
    omicsWorkloadRequireNames(
        prepared,
        c("payload_path", "truth_path", "metadata"),
        "adapter prepare result"
    )
    for (field in c("payload_path", "truth_path")) {
        omicsWorkloadRequireText(prepared[[field]], paste0("prepared.", field))
        if (!file.exists(prepared[[field]]) || dir.exists(prepared[[field]])) {
            omicsWorkloadAbort(sprintf("prepared %s does not exist", field))
        }
        if (!omicsMemoryContainedPath(prepared[[field]], run_dir)) {
            omicsWorkloadAbort(
                sprintf("prepared %s must be contained by the run directory", field),
                class = "omics_workload_privacy_error"
            )
        }
    }
    if (!is.list(prepared$metadata)) {
        omicsWorkloadAbort("prepared.metadata must be a list")
    }
    prepared
}

omicsMemoryPreparationWorker <- function(spec_path) {
    spec <- jsonlite::read_json(spec_path, simplifyVector = FALSE)
    contract <- omicsWorkloadReadContract(spec$contract_path)
    adapter <- omicsWorkloadLoadAdapter(spec$adapter_path, contract)
    prepared <- omicsMemoryPrepareWorkload(
        adapter,
        contract,
        spec$run_dir,
        as.integer(spec$preparation_index)
    )
    result <- list(
        schema = "multischolar.omics_workload_preparation",
        schema_version = "1.0.0",
        status = "passed",
        payload_path = normalizePath(prepared$payload_path, mustWork = TRUE),
        truth_path = normalizePath(prepared$truth_path, mustWork = TRUE),
        payload_sha256 = omicsWorkloadFileDigest(prepared$payload_path),
        truth_sha256 = omicsWorkloadFileDigest(prepared$truth_path),
        metadata = prepared$metadata
    )
    expected <- contract$expected_digests
    if (!identical(result$payload_sha256, expected$payload_sha256) ||
        !identical(result$truth_sha256, expected$truth_sha256)) {
        omicsWorkloadAbort(
            "prepared workload does not match the contract digests",
            class = "omics_workload_binding_error"
        )
    }
    omicsResourceWriteJson(result, spec$result_path)
    invisible(result)
}

omicsMemoryPreparationWorkerSafely <- function(spec_path) {
    spec <- jsonlite::read_json(spec_path, simplifyVector = FALSE)
    tryCatch(
        omicsMemoryPreparationWorker(spec_path),
        error = \(condition) {
            result <- list(
                schema = "multischolar.omics_workload_preparation",
                schema_version = "1.0.0",
                status = "failed",
                reason = conditionMessage(condition),
                condition_class = class(condition)[[1L]]
            )
            omicsResourceWriteJson(result, spec$result_path)
            invisible(result)
        }
    )
}

omicsMemoryPrepareFresh <- function(
    contract,
    contract_path,
    adapter_path,
    output_dir,
    preparation_index
) {
    run_dir <- file.path(
        output_dir,
        "preparations",
        sprintf("p%02d", preparation_index)
    )
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    result_path <- file.path(run_dir, "final", "preparation-result.json")
    spec_path <- tempfile("omics-memory-preparation-", fileext = ".json")
    on.exit(unlink(spec_path, force = TRUE), add = TRUE)
    spec <- list(
        contract_path = contract_path,
        adapter_path = adapter_path,
        run_dir = normalizePath(run_dir, mustWork = TRUE),
        result_path = result_path,
        preparation_index = preparation_index
    )
    omicsResourceWriteJson(spec, spec_path)
    process <- processx::run(
        command = file.path(R.home("bin"), "Rscript"),
        args = c(
            "--vanilla",
            omicsMemoryScriptPath(),
            "--prepare-spec",
            spec_path
        ),
        wd = .OMICS_MEMORY_REPO_ROOT,
        env = c(
            OMP_NUM_THREADS = as.character(contract$execution$threads),
            OPENBLAS_NUM_THREADS = as.character(contract$execution$threads),
            MKL_NUM_THREADS = as.character(contract$execution$threads),
            ARROW_NUM_THREADS = as.character(contract$execution$threads),
            DUCKDB_THREADS = as.character(contract$execution$threads),
            TZ = contract$execution$timezone
        ),
        timeout = as.numeric(contract$execution$timeout_ms),
        error_on_status = FALSE,
        echo = FALSE
    )
    prepared <- if (file.exists(result_path)) {
        jsonlite::read_json(result_path, simplifyVector = FALSE)
    } else {
        list(status = "missing", reason = "preparation result was not written")
    }
    valid <- identical(process$status, 0L) &&
        identical(prepared$status, "passed") &&
        identical(
            prepared$payload_sha256,
            contract$expected_digests$payload_sha256
        ) && identical(
            prepared$truth_sha256,
            contract$expected_digests$truth_sha256
        )
    if (!valid) {
        reason <- if (is.null(prepared$reason)) {
            "fresh workload preparation failed"
        } else {
            prepared$reason
        }
        omicsWorkloadAbort(reason, class = "omics_workload_binding_error")
    }
    if (!omicsMemoryContainedPath(prepared$payload_path, run_dir) ||
        !omicsMemoryContainedPath(prepared$truth_path, run_dir)) {
        omicsWorkloadAbort(
            "fresh preparation returned a path outside its run directory",
            class = "omics_workload_privacy_error"
        )
    }
    prepared
}
