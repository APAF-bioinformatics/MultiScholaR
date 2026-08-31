#!/usr/bin/env Rscript

lipidPilotRunnerPath <- function() {
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(file_arg)) {
        return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
    }
    normalizePath(
        file.path(
            getwd(),
            "tools",
            "profiling",
            "run_omics_publication_lipidomics_pilot.R"
        ),
        mustWork = TRUE
    )
}

.LIPID_PILOT_REPO_ROOT <- normalizePath(
    file.path(dirname(lipidPilotRunnerPath()), "..", ".."),
    mustWork = TRUE
)

for (source_file in c(
    "omics_publication_protocol.R",
    "omics_publication_comparators.R",
    "omics_publication_lock.R",
    "omics_publication_builds.R",
    "omics_publication_restore_reproducibility.R",
    "omics_publication_comparator_builds.R",
    "omics_publication_comparator_build_reproducibility.R",
    "omics_publication_comparator_evidence.R",
    "omics_publication_measure_spec.R",
    "omics_publication_linux_resources.R",
    "omics_publication_retained_resources.R",
    "omics_publication_host_safety.R",
    "omics_publication_lipidomics_contracts.R",
    "omics_publication_lipidomics_model.R",
    "omics_publication_lipidomics_serializers.R",
    "omics_publication_lipidomics_sources.R",
    "omics_publication_lipidomics_truth.R",
    "omics_publication_lipidomics_governance.R",
    "omics_publication_lipidomics_fixtures.R",
    "omics_publication_lipidomics_freeze.R",
    "omics_publication_workload_lipidomics.R",
    "omics_publication_lipidomics_pilot.R"
)) {
    source(
        file.path(
            .LIPID_PILOT_REPO_ROOT,
            "tools",
            "profiling",
            source_file
        ),
        local = FALSE
    )
}

lipidPilotRunnerDefaults <- function() {
    list(
        mode = "pilot",
        contract = NULL,
        payload_root = NULL,
        truth = NULL,
        build_receipt = NULL,
        output_root = NULL,
        result = NULL,
        run_dir = NULL,
        package_library = NULL,
        dependency_library = NULL,
        historical_revision = NULL,
        dwell_seconds = 5
    )
}

lipidPilotRunnerArgs <- function(argv) {
    values <- lipidPilotRunnerDefaults()
    index <- 1L
    while (index <= length(argv)) {
        token <- argv[[index]]
        key <- gsub("-", "_", sub("^--", "", token), fixed = TRUE)
        if (!startsWith(token, "--") || !key %in% names(values) ||
            index == length(argv)) {
            stop("Pilot runner arguments are invalid", call. = FALSE)
        }
        value <- argv[[index + 1L]]
        values[[key]] <- if (identical(key, "dwell_seconds")) {
            as.numeric(value)
        } else {
            value
        }
        index <- index + 2L
    }
    values
}

lipidPilotRequired <- function(args, fields) {
    missing <- fields[vapply(fields, \(field) {
        is.null(args[[field]]) || !nzchar(args[[field]])
    }, logical(1))]
    if (length(missing)) {
        stop(paste("Pilot option is missing:", missing[[1L]]), call. = FALSE)
    }
    invisible(args)
}

lipidPilotWorkerMain <- function(args) {
    lipidPilotRequired(args, c(
        "contract", "payload_root", "truth", "run_dir", "package_library",
        "dependency_library", "historical_revision"
    ))
    lipidPublicationPilotWorker(
        args$contract,
        args$payload_root,
        args$truth,
        args$run_dir,
        args$package_library,
        args$dependency_library,
        args$historical_revision,
        args$dwell_seconds
    )
    invisible(0L)
}

lipidPilotHostPreflight <- function(output_root) {
    contract <- publicationReadJson(
        "tests/testdata/omics-performance/host-preflight-contract-v3.json"
    )
    host <- publicationHostEnvelope(output_root)
    preflight <- publicationEvaluateHostPreflight(
        host,
        contract,
        expected_peak_bytes = 4 * 1024^3,
        maximum_temporary_bytes = 50 * 1024^3
    )
    list(host = host, preflight = preflight)
}

lipidPilotLibraries <- function(build_receipt_path) {
    evidence <- publicationValidateComparatorBuildReceiptEvidence(
        build_receipt_path
    )
    receipt <- evidence$receipt
    if (!identical(receipt$comparator_id, "historical_janitor")) {
        lipidPublicationAbort("pilot comparator is not historical")
    }
    list(
        receipt = receipt,
        package_library = dirname(receipt$build$installed_inventory$root),
        dependency_library = file.path(
            dirname(receipt$environment$restore_receipt$path),
            "library"
        )
    )
}

lipidPilotUnavailableRecord <- function(
    contract_path,
    build_receipt_path,
    preflight
) {
    list(
        schema = "multischolar.omics_publication_lipidomics_pilot",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-067",
        status = "host_not_publication_certified",
        comparator_role = "historical_janitor",
        comparator_revision = publicationComparatorRevisions()[[
            "historical_janitor"
        ]],
        build_receipt_sha256 = publicationFileDigest(build_receipt_path),
        contract_sha256 = publicationFileDigest(contract_path),
        preflight = preflight,
        candidate_loaded = FALSE,
        promotion_authority = FALSE,
        publication_authority = FALSE
    )
}

lipidPilotParentMain <- function(args) {
    lipidPilotRequired(args, c(
        "contract", "payload_root", "truth", "build_receipt", "output_root",
        "result"
    ))
    if (file.exists(args$output_root) || dir.exists(args$output_root) ||
        file.exists(args$result)) {
        stop("Pilot output already exists", call. = FALSE)
    }
    contract <- publicationReadJson(args$contract)
    lipidPublicationValidatePilotContract(contract)
    members <- file.path(
        args$payload_root,
        unlist(contract$generator$output_members, use.names = FALSE)
    )
    payload <- lipidPublicationPayloadBinding(members)
    if (!identical(
        payload$payload_set_sha256,
        contract$expected_digests$payload_set_sha256
    ) || !identical(
        publicationFileDigest(args$truth),
        contract$expected_digests$truth_sha256
    )) {
        lipidPublicationAbort("pilot input digest differs")
    }
    libraries <- lipidPilotLibraries(args$build_receipt)
    dir.create(args$output_root, recursive = TRUE, showWarnings = FALSE)
    home <- file.path(args$output_root, "home")
    temp <- file.path(args$output_root, "tmp")
    site <- file.path(args$output_root, "empty-site-library")
    dir.create(home)
    dir.create(temp)
    dir.create(site)
    host <- lipidPilotHostPreflight(args$output_root)
    if (!isTRUE(host$preflight$certified)) {
        publicationWriteJson(
            lipidPilotUnavailableRecord(
                args$contract,
                args$build_receipt,
                host$preflight
            ),
            args$result
        )
        return(invisible(0L))
    }
    run_dir <- file.path(args$output_root, "run")
    worker_args <- c(
        "--vanilla",
        lipidPilotRunnerPath(),
        "--mode", "worker",
        "--contract", normalizePath(args$contract),
        "--payload-root", normalizePath(args$payload_root),
        "--truth", normalizePath(args$truth),
        "--run-dir", normalizePath(run_dir, mustWork = FALSE),
        "--package-library", normalizePath(libraries$package_library),
        "--dependency-library", normalizePath(libraries$dependency_library),
        "--historical-revision", publicationComparatorRevisions()[[
            "historical_janitor"
        ]],
        "--dwell-seconds", as.character(args$dwell_seconds)
    )
    limits <- lipidPublicationPilotSafetyLimits(
        host$host,
        host$preflight
    )
    measurement <- publicationMeasureCgroupSubprocess(
        command = file.path(R.home("bin"), "Rscript"),
        command_args = worker_args,
        run_dir = run_dir,
        execution = lipidPublicationPilotExecution(),
        env = lipidPublicationPilotThreadEnvironment(
            home,
            temp,
            libraries$package_library,
            libraries$dependency_library,
            site
        ),
        unit_name = publicationSystemdUnitName("multischolar-lipidomics-pilot"),
        require_certified_host = TRUE,
        host_preflight = host$preflight,
        safety_check_fn = publicationRuntimeSafetyMonitor(
            limits,
            args$output_root
        )
    )
    measurement <- tryCatch(
        lipidPublicationAttachPilotPhase(measurement, run_dir),
        error = \(error) {
            measurement$status <- "failed"
            measurement$publication_certifiable <- FALSE
            measurement$phase_evidence <- list(
                valid = FALSE,
                reason = conditionMessage(error)
            )
            measurement
        }
    )
    raw_path <- file.path(args$output_root, "measurement.json")
    publicationWriteJson(measurement, raw_path)
    record <- lipidPublicationPilotRecord(
        args$contract,
        args$payload_root,
        args$truth,
        args$build_receipt,
        host$host,
        host$preflight,
        measurement,
        raw_path
    )
    publicationWriteJson(record, args$result)
    invisible(0L)
}

lipidPilotRunnerMain <- function(argv = commandArgs(trailingOnly = TRUE)) {
    args <- lipidPilotRunnerArgs(argv)
    if (identical(args$mode, "worker")) {
        lipidPilotWorkerMain(args)
    } else if (identical(args$mode, "pilot")) {
        lipidPilotParentMain(args)
    } else {
        stop("Pilot mode is invalid", call. = FALSE)
    }
}

if (identical(environment(), globalenv()) && !interactive()) {
    status <- tryCatch(
        lipidPilotRunnerMain(),
        error = \(error) {
            message(conditionMessage(error))
            1L
        }
    )
    quit(status = as.integer(status), save = "no")
}
