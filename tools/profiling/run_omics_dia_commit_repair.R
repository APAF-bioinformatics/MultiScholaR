#!/usr/bin/env Rscript

diaRepairRunnerPath <- function() {
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(file_arg)) {
        return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
    }
    normalizePath(
        file.path(
            getwd(),
            "tools",
            "profiling",
            "run_omics_dia_commit_repair.R"
        ),
        mustWork = TRUE
    )
}

.DIA_REPAIR_REPO_ROOT <- normalizePath(
    file.path(dirname(diaRepairRunnerPath()), "..", ".."),
    mustWork = TRUE
)

for (source_file in c(
    "omics_publication_protocol.R",
    "omics_publication_linux_resources.R",
    "omics_publication_retained_resources.R",
    "omics_publication_host_safety.R",
    "omics_dia_commit_repair.R",
    "omics_dia_commit_repair_proof.R"
)) {
    source(
        file.path(
            .DIA_REPAIR_REPO_ROOT,
            "tools",
            "profiling",
            source_file
        ),
        local = FALSE
    )
}

diaRepairRunnerDefaults <- function() {
    list(
        mode = "pilot",
        source = NULL,
        output_root = NULL,
        result = NULL,
        pre_library = NULL,
        candidate_library = NULL,
        pre_revision = NULL,
        candidate_revision = NULL,
        salt_file = NULL,
        prior_result = NULL,
        pairs = 2L,
        arm = NULL,
        run_dir = NULL,
        package_library = NULL,
        dwell_seconds = 5
    )
}

diaRepairRunnerArgs <- function(argv) {
    values <- diaRepairRunnerDefaults()
    index <- 1L
    while (index <= length(argv)) {
        token <- argv[[index]]
        key <- gsub("-", "_", sub("^--", "", token), fixed = TRUE)
        if (!startsWith(token, "--") || !key %in% names(values) ||
            index == length(argv)) {
            diaRepairAbort("DIA repair runner arguments are invalid")
        }
        value <- argv[[index + 1L]]
        values[[key]] <- if (identical(key, "pairs")) {
            as.integer(value)
        } else if (identical(key, "dwell_seconds")) {
            as.numeric(value)
        } else {
            value
        }
        index <- index + 2L
    }
    values
}

diaRepairRunnerRequire <- function(args, fields) {
    missing <- fields[vapply(fields, function(field) {
        is.null(args[[field]]) || !nzchar(args[[field]])
    }, logical(1))]
    if (length(missing)) {
        diaRepairAbort(paste("DIA repair option is missing:", missing[[1L]]))
    }
    invisible(args)
}

diaRepairDirectoryDigest <- function(path) {
    root <- normalizePath(path, mustWork = TRUE)
    files <- list.files(
        root,
        recursive = TRUE,
        full.names = TRUE,
        all.files = TRUE,
        no.. = TRUE
    )
    files <- sort(files[file.exists(files) & !dir.exists(files)], method = "radix")
    relative <- substring(files, nchar(root) + 2L)
    records <- Map(function(file, relative_path) {
        list(path = relative_path, sha256 = publicationFileDigest(file))
    }, files, relative)
    publicationObjectDigest(unname(records))
}

diaRepairHostPreflight <- function(output_root, gates) {
    contract <- publicationReadJson(
        "tests/testdata/omics-performance/host-preflight-contract-v3.json"
    )
    host <- publicationHostEnvelope(output_root)
    preflight <- diaRepairEngineeringPreflight(host, contract, gates)
    list(host = host, preflight = preflight)
}

diaRepairWarmInput <- function(path, chunk_bytes = 1024L * 1024L) {
    connection <- file(path, open = "rb")
    on.exit(close(connection), add = TRUE)
    repeat {
        value <- readBin(connection, what = "raw", n = chunk_bytes)
        if (!length(value)) break
    }
    invisible(TRUE)
}

diaRepairEvidenceSalt <- function(salt_file = NULL) {
    salt <- Sys.getenv("MULTISCHOLAR_DIA_REPAIR_SALT", unset = "")
    if (nzchar(salt)) return(salt)
    if (is.null(salt_file)) {
        salt_file <- Sys.getenv(
            "MULTISCHOLAR_DIA_REPAIR_SALT_FILE",
            unset = ""
        )
    }
    if (!nzchar(salt_file) || !file.exists(salt_file)) {
        diaRepairAbort("DIA repair evidence salt is unavailable")
    }
    salt <- readLines(salt_file, n = 1L, warn = FALSE)
    if (length(salt) != 1L || !nzchar(salt)) {
        diaRepairAbort("DIA repair evidence salt file is invalid")
    }
    salt
}

diaRepairCompactMeasurement <- function(measurement, raw_path) {
    list(
        status = measurement$status,
        publication_certifiable = measurement$publication_certifiable,
        exit_status = measurement$exit_status,
        timed_out = measurement$timed_out,
        safety_aborted = measurement$safety_aborted,
        safety_reason = measurement$safety_reason,
        cgroup_observed = measurement$cgroup_observed,
        cgroup_lost = measurement$cgroup_lost,
        pid_ownership_valid = measurement$pid_ownership_valid,
        retention_acknowledged = measurement$retention_acknowledged,
        worker_resource_ledger_valid = measurement$worker_resource_ledger_valid,
        cleanup = measurement$cleanup,
        metrics = measurement$metrics,
        phase_evidence = measurement$phase_evidence,
        worker = measurement$worker,
        raw_measurement_sha256 = publicationFileDigest(raw_path)
    )
}

diaRepairRunArm <- function(
    slot,
    args,
    host,
    limits,
    salt,
    readiness
) {
    arm <- slot$arm[[1L]]
    package_library <- if (identical(arm, "pre_repair_artifact")) {
        args$pre_library
    } else {
        args$candidate_library
    }
    run_dir <- file.path(
        args$output_root,
        "runs",
        slot$pair_id[[1L]],
        sprintf("%d-%s", slot$sequence_in_pair[[1L]], arm)
    )
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    private_root <- file.path(run_dir, "private-project")
    private_logs <- c(
        file.path(run_dir, "stdout.log"),
        file.path(run_dir, "stderr.log")
    )
    on.exit({
        unlink(private_root, recursive = TRUE, force = TRUE)
        unlink(private_logs, force = TRUE)
    }, add = TRUE)
    home <- file.path(run_dir, "home")
    temp <- file.path(run_dir, "tmp")
    site <- file.path(run_dir, "empty-site-library")
    dir.create(home)
    dir.create(temp)
    dir.create(site)
    library_path <- paste(
        unique(c(normalizePath(package_library), .libPaths())),
        collapse = .Platform$path.sep
    )
    environment <- diaRepairThreadEnvironment(
        home,
        temp,
        library_path,
        site
    )
    environment <- c(
        environment,
        MULTISCHOLAR_DIA_REPAIR_SOURCE = normalizePath(args$source),
        MULTISCHOLAR_DIA_REPAIR_SALT = salt
    )
    worker_args <- c(
        "--vanilla",
        diaRepairRunnerPath(),
        "--mode", "worker",
        "--arm", arm,
        "--run-dir", normalizePath(run_dir, mustWork = FALSE),
        "--package-library", normalizePath(package_library),
        "--dwell-seconds", as.character(args$dwell_seconds)
    )
    measurement <- publicationMeasureCgroupSubprocess(
        command = file.path(R.home("bin"), "Rscript"),
        command_args = worker_args,
        run_dir = run_dir,
        execution = diaRepairExecution(),
        env = environment,
        unit_name = publicationSystemdUnitName("multischolar-dia-repair"),
        require_certified_host = TRUE,
        host_preflight = host$preflight,
        safety_check_fn = publicationRuntimeSafetyMonitor(
            limits,
            args$output_root
        )
    )
    measurement <- tryCatch(
        diaRepairAttachPhase(measurement, run_dir),
        error = function(error) {
            measurement$status <- "failed"
            measurement$publication_certifiable <- FALSE
            measurement$phase_evidence <- list(
                valid = FALSE,
                reason = conditionMessage(error)
            )
            measurement
        }
    )
    unlink(private_root, recursive = TRUE, force = TRUE)
    unlink(private_logs, force = TRUE)
    measurement$stdout_path <- "<private-log-deleted>"
    measurement$stderr_path <- "<private-log-deleted>"
    raw_path <- file.path(run_dir, "measurement.json")
    publicationWriteJson(measurement, raw_path)
    list(
        pair_id = slot$pair_id[[1L]],
        session_id = slot$session_id[[1L]],
        sequence_in_pair = as.integer(slot$sequence_in_pair[[1L]]),
        arm = arm,
        readiness = readiness,
        measurement = diaRepairCompactMeasurement(measurement, raw_path)
    )
}

diaRepairUnavailableResult <- function(args, gates, host) {
    list(
        schema = "multischolar.dia_commit_repair_campaign",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-070",
        status = "host_not_engineering_certified",
        gate_binding = list(
            gate_id = gates$gate_id,
            sha256 = publicationFileDigest(diaRepairGatePath())
        ),
        host_preflight = host$preflight,
        pairs_requested = args$pairs,
        candidate_work_executed = FALSE,
        automatic_policy_authority = FALSE,
        publication_authority = FALSE
    )
}

diaRepairResumeEvidence <- function(args, gates, source_size, salt) {
    if (!identical(args$mode, "resume")) {
        return(list(records = list(), completed_pair_ids = character()))
    }
    diaRepairRunnerRequire(args, "prior_result")
    prior <- publicationReadJson(args$prior_result)
    gate_valid <- identical(prior$gate_binding$gate_id, gates$gate_id) &&
        identical(
            prior$gate_binding$sha256,
            publicationFileDigest(diaRepairGatePath())
        )
    source_valid <- identical(
        as.numeric(prior$source_binding$exact_input_bytes),
        source_size
    ) && identical(
        prior$source_binding$salted_source_fingerprint,
        diaRepairSaltedDigest(publicationFileDigest(args$source), salt)
    )
    pre_binding <- prior$comparator_bindings$pre_repair_artifact
    candidate_binding <- prior$comparator_bindings$candidate_artifact
    comparator_valid <- identical(pre_binding$revision, args$pre_revision) &&
        identical(candidate_binding$revision, args$candidate_revision) &&
        identical(
            pre_binding$installed_package_sha256,
            diaRepairDirectoryDigest(file.path(args$pre_library, "MultiScholaR"))
        ) && identical(
            candidate_binding$installed_package_sha256,
            diaRepairDirectoryDigest(file.path(
                args$candidate_library,
                "MultiScholaR"
            ))
        )
    if (!is.list(prior$records) || !gate_valid || !source_valid ||
        !comparator_valid || isTRUE(prior$automatic_policy_authority) ||
        isTRUE(prior$publication_authority)) {
        diaRepairAbort("DIA repair resume evidence is invalid")
    }
    records <- diaRepairCompletePairRecords(prior$records)
    completed <- unique(vapply(records, `[[`, character(1), "pair_id"))
    expected_prefix <- sprintf("dia-repair-pair-%03d", seq_along(completed))
    if (!identical(sort(completed, method = "radix"), expected_prefix)) {
        diaRepairAbort("DIA repair resume pairs are not one complete prefix")
    }
    list(
        records = records,
        completed_pair_ids = completed,
        binding = list(
            path = normalizePath(args$prior_result),
            sha256 = publicationFileDigest(args$prior_result),
            retained_complete_pairs = length(completed)
        )
    )
}

diaRepairCampaignMain <- function(args) {
    diaRepairRunnerRequire(args, c(
        "source", "output_root", "result", "pre_library",
        "candidate_library", "pre_revision", "candidate_revision",
        "salt_file"
    ))
    gates <- diaRepairReadGates()
    if (args$mode %in% c("campaign", "resume") &&
        !identical(as.integer(args$pairs), 36L)) {
        diaRepairAbort("Confirmatory DIA repair campaign requires 36 pairs")
    }
    if (args$pairs < 2L || args$pairs %% 2L != 0L) {
        diaRepairAbort("DIA repair pilot requires a positive even pair count")
    }
    if (file.exists(args$result) || dir.exists(args$output_root)) {
        diaRepairAbort("DIA repair output already exists")
    }
    source_size <- as.numeric(file.info(args$source)$size)
    if (!is.finite(source_size) || source_size <= 0) {
        diaRepairAbort("DIA repair source is unavailable")
    }
    salt <- diaRepairEvidenceSalt(args$salt_file)
    resumed <- diaRepairResumeEvidence(args, gates, source_size, salt)
    dir.create(args$output_root, recursive = TRUE, showWarnings = FALSE)
    host <- diaRepairHostPreflight(args$output_root, gates)
    initial_readiness <- NULL
    failed_checks <- names(host$preflight$checks)[
        !unlist(host$preflight$checks, use.names = FALSE)
    ]
    if (!isTRUE(host$preflight$certified) &&
        identical(failed_checks, "load")) {
        initial_readiness <- diaRepairWaitForReadiness(host$host, gates)
        if (isTRUE(initial_readiness$ready)) {
            host <- diaRepairHostPreflight(args$output_root, gates)
        }
    }
    if (!isTRUE(host$preflight$certified)) {
        result <- diaRepairUnavailableResult(args, gates, host)
        result$initial_readiness <- initial_readiness
        publicationWriteJson(result, args$result)
        return(invisible(0L))
    }
    session_count <- if (args$pairs %% 3L == 0L) 3L else 1L
    schedule <- diaRepairSchedule(args$pairs, session_count)
    pending <- schedule[
        !schedule$pair_id %in% resumed$completed_pair_ids,
        ,
        drop = FALSE
    ]
    records <- resumed$records
    readiness_failure <- NULL
    for (index in seq_len(nrow(pending))) {
        readiness <- diaRepairWaitForReadiness(host$host, gates)
        if (!isTRUE(readiness$ready)) {
            readiness_failure <- readiness
            break
        }
        limits <- diaRepairRuntimeSafetyLimits(
            host$host,
            host$preflight,
            gates,
            readiness$baseline_one_minute_load
        )
        diaRepairWarmInput(args$source)
        record <- diaRepairRunArm(
            pending[index, , drop = FALSE],
            args,
            host,
            limits,
            salt,
            readiness
        )
        records[[length(records) + 1L]] <- record
        measurement <- record$measurement
        terminal_failure <- isTRUE(measurement$safety_aborted) ||
            isTRUE(measurement$timed_out) ||
            !identical(as.integer(measurement$exit_status), 0L)
        if (terminal_failure) break
    }
    evaluation <- diaRepairEvaluate(records, gates)
    mode_status <- if (args$mode %in% c("campaign", "resume")) {
        evaluation$status
    } else if (!is.null(readiness_failure)) {
        "pilot_readiness_timeout"
    } else if (!all(vapply(records, function(record) {
        identical(record$measurement$status, "passed")
    }, logical(1)))) {
        "pilot_failed_or_safety_aborted"
    } else {
        "pilot_complete_non_authoritative"
    }
    result <- list(
        schema = "multischolar.dia_commit_repair_campaign",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-070",
        status = mode_status,
        mode = args$mode,
        gate_binding = list(
            gate_id = gates$gate_id,
            sha256 = publicationFileDigest(diaRepairGatePath())
        ),
        source_binding = list(
            source_id = "private.dia.cotton_report.v1",
            exact_input_bytes = source_size,
            salted_source_fingerprint = diaRepairSaltedDigest(
                publicationFileDigest(args$source),
                salt
            ),
            private_values_retained = FALSE
        ),
        comparator_bindings = list(
            pre_repair_artifact = list(
                revision = args$pre_revision,
                installed_package_sha256 = diaRepairDirectoryDigest(
                    file.path(args$pre_library, "MultiScholaR")
                )
            ),
            candidate_artifact = list(
                revision = args$candidate_revision,
                installed_package_sha256 = diaRepairDirectoryDigest(
                    file.path(args$candidate_library, "MultiScholaR")
                )
            )
        ),
        host_preflight = list(
            certified = host$preflight$certified,
            sha256 = publicationObjectDigest(host$preflight)
        ),
        initial_readiness = initial_readiness,
        schedule = split(schedule, seq_len(nrow(schedule))),
        resumed_from = resumed$binding,
        records = records,
        readiness_failure = readiness_failure,
        evaluation = evaluation,
        candidate_work_executed = any(vapply(records, function(record) {
            identical(record$arm, "candidate_artifact")
        }, logical(1))),
        automatic_policy_authority = FALSE,
        publication_authority = FALSE
    )
    publicationWriteJson(result, args$result)
    invisible(if (args$mode %in% c("campaign", "resume") &&
        !identical(evaluation$status, "passed")) 1L else 0L)
}

diaRepairWorkerMain <- function(args) {
    diaRepairRunnerRequire(args, c("arm", "run_dir", "package_library"))
    source <- Sys.getenv("MULTISCHOLAR_DIA_REPAIR_SOURCE", unset = "")
    salt <- Sys.getenv("MULTISCHOLAR_DIA_REPAIR_SALT", unset = "")
    diaRepairWorker(
        source,
        args$run_dir,
        args$package_library,
        args$arm,
        salt,
        args$dwell_seconds
    )
}

diaRepairProofMain <- function(args) {
    diaRepairRunnerRequire(args, c("arm", "run_dir", "package_library"))
    salt <- Sys.getenv("MULTISCHOLAR_DIA_REPAIR_SALT", unset = "")
    proof <- diaRepairProofWorker(
        args$package_library,
        args$run_dir,
        args$arm,
        salt
    )
    publicationWriteJson(
        proof,
        file.path(args$run_dir, "scientific-proof.json")
    )
    invisible(0L)
}

diaRepairDiagnosticMain <- function(args) {
    diaRepairRunnerRequire(args, c(
        "source", "result", "arm", "run_dir", "package_library",
        "salt_file"
    ))
    salt <- diaRepairEvidenceSalt(args$salt_file)
    .libPaths(c(args$package_library, .libPaths()))
    namespace <- loadNamespace("MultiScholaR", lib.loc = args$package_library)
    paths <- diaRepairWorkflowPaths(args$run_dir)
    on.exit(unlink(paths$base_dir, recursive = TRUE, force = TRUE), add = TRUE)
    rollout <- if (identical(args$arm, "candidate_artifact")) {
        "evict"
    } else {
        "dual_write"
    }
    workflow <- diaRepairNewWorkflow(namespace, paths, rollout)
    imported <- suppressMessages(diaRepairPrepareImport(
        namespace,
        workflow,
        args$source
    ))
    prepared <- suppressMessages(diaRepairPrepareDesign(
        namespace,
        workflow,
        imported,
        paths
    ))
    commit <- suppressMessages(diaRepairRunCommit(
        namespace,
        workflow,
        paths,
        prepared
    ))
    proof <- diaRepairProofWorker(
        args$package_library,
        args$run_dir,
        args$arm,
        salt
    )
    publicationWriteJson(list(
        status = "diagnostic_non_authoritative",
        arm = args$arm,
        proof = proof,
        complete_payload_returned = isTRUE(
            commit$design$hydration_verification$complete_payload_returned
        ) || isTRUE(commit$settlement$complete_payload_returned),
        promotion_authority = FALSE,
        publication_authority = FALSE
    ), args$result)
    invisible(0L)
}

diaRepairRunnerMain <- function(argv = commandArgs(trailingOnly = TRUE)) {
    args <- diaRepairRunnerArgs(argv)
    if (identical(args$mode, "worker")) {
        diaRepairWorkerMain(args)
    } else if (identical(args$mode, "proof")) {
        diaRepairProofMain(args)
    } else if (identical(args$mode, "diagnostic")) {
        diaRepairDiagnosticMain(args)
    } else if (args$mode %in% c("pilot", "campaign", "resume")) {
        diaRepairCampaignMain(args)
    } else {
        diaRepairAbort("DIA repair runner mode is invalid")
    }
}

if (identical(environment(), globalenv()) && !interactive()) {
    status <- tryCatch(
        diaRepairRunnerMain(),
        error = function(error) {
            message(conditionMessage(error))
            1L
        }
    )
    quit(status = as.integer(status), save = "no")
}
