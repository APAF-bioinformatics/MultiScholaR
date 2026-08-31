#!/usr/bin/env Rscript

metabStressRunnerPath <- function() {
    file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(file_arg)) {
        return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
    }
    normalizePath(
        file.path(
            getwd(),
            "tools",
            "profiling",
            "run_omics_publication_metabolomics_stress.R"
        ),
        mustWork = TRUE
    )
}

.METAB_STRESS_REPO_ROOT <- normalizePath(
    file.path(dirname(metabStressRunnerPath()), "..", ".."),
    mustWork = TRUE
)

source(file.path(
    .METAB_STRESS_REPO_ROOT,
    "tools",
    "profiling",
    "omics_publication_protocol.R"
))

metabStressRunnerArgs <- function(argv) {
    values <- list(
        unit = NULL,
        contract = NULL,
        payload_root = NULL,
        truth = NULL,
        build_receipt = NULL,
        output = NULL,
        receipt = NULL
    )
    index <- 1L
    while (index <= length(argv)) {
        token <- argv[[index]]
        key <- gsub("-", "_", sub("^--", "", token), fixed = TRUE)
        if (!startsWith(token, "--") || !key %in% names(values) ||
            index == length(argv)) {
            stop("Metabolomics stress arguments are invalid", call. = FALSE)
        }
        values[[key]] <- argv[[index + 1L]]
        index <- index + 2L
    }
    if (any(vapply(values, is.null, logical(1)))) {
        stop("Metabolomics stress options are incomplete", call. = FALSE)
    }
    if (!grepl("^[a-z0-9-]+$", values$unit)) {
        stop("Metabolomics stress unit name is invalid", call. = FALSE)
    }
    values
}

metabStressResultField <- function(output, pattern) {
    line <- grep(pattern, output, value = TRUE)
    if (!length(line)) return(NULL)
    trimws(sub(pattern, "", line[[length(line)]]))
}

metabStressRunnerMain <- function(argv = commandArgs(trailingOnly = TRUE)) {
    args <- metabStressRunnerArgs(argv)
    if (file.exists(args$output) || file.exists(args$receipt)) {
        stop("Metabolomics stress output already exists", call. = FALSE)
    }
    run_dir <- dirname(args$receipt)
    dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
    stdout_path <- file.path(run_dir, "stdout.log")
    stderr_path <- file.path(run_dir, "stderr.log")
    import_runner <- file.path(
        .METAB_STRESS_REPO_ROOT,
        "tools",
        "profiling",
        "run_omics_publication_metabolomics_import.R"
    )
    command_args <- c(
        "--vanilla", import_runner,
        "--contract", normalizePath(args$contract, mustWork = TRUE),
        "--payload-root", normalizePath(args$payload_root, mustWork = TRUE),
        "--truth", normalizePath(args$truth, mustWork = TRUE),
        "--build-receipt", normalizePath(args$build_receipt, mustWork = TRUE),
        "--output", normalizePath(args$output, mustWork = FALSE)
    )
    systemd_args <- c(
        "--user", "--wait", "--pipe", "--collect",
        paste0("--unit=", args$unit),
        paste0("--working-directory=", .METAB_STRESS_REPO_ROOT),
        "--property=MemoryAccounting=yes",
        "--property=IOAccounting=yes",
        "--property=TasksMax=512",
        "--property=KillMode=mixed",
        "--property=OOMPolicy=stop",
        "--property=MemoryMax=10737418240",
        "--property=MemorySwapMax=0",
        Sys.which("Rscript"),
        command_args
    )
    result <- processx::run(
        "systemd-run",
        systemd_args,
        stdout = stdout_path,
        stderr = stderr_path,
        error_on_status = FALSE,
        timeout = 1800000,
        echo = FALSE
    )
    combined <- c(
        readLines(stdout_path, warn = FALSE),
        readLines(stderr_path, warn = FALSE)
    )
    success <- identical(result$status, 0L) && file.exists(args$output)
    receipt <- list(
        schema = "multischolar.omics_publication_metabolomics_stress_execution",
        schema_version = "1.0.0",
        owner_ticket_id = "OMICS-ART-066",
        status = if (success) "passed" else "failed",
        unit_name = args$unit,
        contract = list(
            path = args$contract,
            sha256 = publicationFileDigest(args$contract)
        ),
        payload_root = args$payload_root,
        truth_sha256 = publicationFileDigest(args$truth),
        build_receipt_sha256 = publicationFileDigest(args$build_receipt),
        output_sha256 = if (file.exists(args$output)) {
            publicationFileDigest(args$output)
        } else {
            NULL
        },
        systemd_exit_status = as.integer(result$status),
        systemd_result = metabStressResultField(
            combined,
            "^.*Finished with result:[[:space:]]*"
        ),
        service_runtime = metabStressResultField(
            combined,
            "^.*Service runtime:[[:space:]]*"
        ),
        memory_peak = metabStressResultField(
            combined,
            "^.*Memory peak:[[:space:]]*"
        ),
        memory_max_bytes = 10737418240,
        swap_max_bytes = 0,
        stdout_sha256 = publicationFileDigest(stdout_path),
        stderr_sha256 = publicationFileDigest(stderr_path),
        candidate_loaded = FALSE,
        performance_authority = FALSE,
        promotion_authority = FALSE,
        publication_authority = FALSE
    )
    publicationWriteJson(receipt, args$receipt)
    if (!success) {
        stop("Metabolomics stress import failed", call. = FALSE)
    }
    invisible(0L)
}

if (identical(environment(), globalenv()) && !interactive()) {
    status <- tryCatch(
        metabStressRunnerMain(),
        error = function(error) {
            message(conditionMessage(error))
            1L
        }
    )
    quit(status = as.integer(status), save = "no")
}
