# Provenance, one per job (KIT-159 section 5.6, `msrval.provenance.v1`).

PROVENANCE_PHASES <- c("record", "pre", "lifecycle", "counted", "score", "determinism",
                       "correction", "amended", "parity", "check")

#' The runner's identity (G4): image, CPU, cores, memory, kernel; `unknown` where unreadable
#' @noRd
valRunnerInfo <- function() {
    env <- function(x) {
        v <- Sys.getenv(x)
        if (nzchar(v)) v else "unknown"
    }
    cpu <- tryCatch(grep("^model name", readLines("/proc/cpuinfo"), value = TRUE),
                    error = function(e) character())
    mem <- tryCatch(grep("^MemTotal", readLines("/proc/meminfo"), value = TRUE),
                    error = function(e) character())
    list(
        image_os = env("ImageOS"), image_version = env("ImageVersion"),
        cpu_model = if (length(cpu)) trimws(sub("^[^:]*:", "", cpu[[1L]])) else "unknown",
        cores = max(1L, length(cpu)),
        mem_total_mib = if (length(mem)) {
            round(as.numeric(gsub("[^0-9]", "", mem[[1L]])) / 1024)
        } else {
            NA_real_
        },
        kernel = unname(Sys.info()[["release"]])
    )
}

#' The job's run block
#' @noRd
valRunInfo <- function() {
    list(github_run_id = Sys.getenv("GITHUB_RUN_ID"),
         run_attempt = Sys.getenv("GITHUB_RUN_ATTEMPT"),
         job = Sys.getenv("GITHUB_JOB"), workflow = Sys.getenv("GITHUB_WORKFLOW"),
         request_file = Sys.getenv("REQUEST"), request_commit = Sys.getenv("GITHUB_SHA"),
         runner = valRunnerInfo())
}

#' The job's environment block
#' @noRd
valEnvironmentInfo <- function() {
    version_of <- function(p) {
        tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
    }
    pkgs <- c("MultiScholaR", "limma", "edgeR", "fastduckdb", "fastlimpadpc", "duckdb", "arrow",
              "shiny", "testthat", "msrval")
    list(
        r_version = R.version.string, bioc_version = version_of("BiocVersion"),
        lockfile_sha256 = NA_character_, library_sha256 = NA_character_,
        differs_from_lockfile = list(), environment_check = NA_character_,
        packages = as.list(stats::setNames(vapply(pkgs, version_of, character(1)), pkgs)),
        blas = as.character(extSoftVersion()[["BLAS"]]), lapack = La_library(),
        locale = Sys.getlocale(), timezone = Sys.timezone(),
        threads = as.list(Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                                       "MKL_NUM_THREADS"))),
        test_mode = getOption("multischolar.test_mode", NA), python = NULL
    )
}

#' Build a provenance record (section 5.6)
#'
#' @param ticket,stage,label,phase,status Top fields.
#' @param blinded Phase lifecycle scoring blinded.
#' @param started,ended UTC ISO 8601 times.
#' @param refusal `list(class, step)` or `NULL`.
#' @param amended Amendment label or `NULL`.
#' @param code,inputs,truth,cells,chunks,scoring,outputs Blocks of section 5.6.
#' @param environment Environment block (default: read from the job).
#' @param run Run block (default: read from the job).
#' @return The provenance as a list.
#' @export
valProvenance <- function(ticket, stage, label, phase, status, blinded = FALSE,
                          started = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
                          ended = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
                          refusal = NULL, amended = NULL, code = list(), inputs = list(),
                          truth = list(), cells = list(), chunks = list(), scoring = list(),
                          outputs = list(), environment = NULL, run = NULL) {
    list(
        schema = "msrval.provenance.v1", ticket = ticket, stage = stage, label = label,
        phase = phase, amended = amended, blinded = blinded, started = started, ended = ended,
        status = status, refusal = refusal, run = run %||% valRunInfo(), code = code,
        environment = environment %||% valEnvironmentInfo(), inputs = inputs, truth = truth,
        cells = cells, chunks = chunks, scoring = scoring, outputs = outputs
    )
}

#' Write a provenance record as JSON
#'
#' @param prov A provenance list.
#' @param path Output path.
#' @return `path`, invisibly.
#' @export
valWriteProvenance <- function(prov, path) {
    jsonlite::write_json(prov, path, auto_unbox = TRUE, null = "null", na = "null", pretty = TRUE,
                         digits = NA)
    invisible(path)
}

#' Hashes of the per-unit truths the scoring used (5.6; E4 of the re-audit): unit -> truth id ->
#' `sha256`, `rows`, `by_truth_status`, as [valTruthProvenance()] gives for the whole truth
#' @noRd
valUnitTruthsProvenance <- function(truths, work_dir) {
    lapply(truths, function(tu) {
        x <- valTruthProvenance(tu, work_dir, NULL)
        x[setdiff(names(x), "blind")]
    })
}

#' Check a provenance record's required fields and types (an invalid file fails the job)
#'
#' @param x A provenance list or the path of `provenance.json`.
#' @return `TRUE`, or `val_provenance_invalid`.
#' @export
valCheckProvenance <- function(x) {
    if (is.character(x)) x <- jsonlite::fromJSON(x, simplifyVector = FALSE)
    bad <- function(what) valAbort("provenance: {what}", "val_provenance_invalid")
    is_chr <- function(v) is.character(v) && length(v) == 1L
    need <- function(obj, fields, check, where) {
        for (f in fields) {
            if (is.null(obj[[f]]) || !check(obj[[f]])) bad(paste0(where, f))
        }
    }
    need(x, c("schema", "ticket", "stage", "label", "phase", "started", "ended", "status"),
         is_chr, "")
    if (!identical(x$schema, "msrval.provenance.v1")) bad("schema")
    if (!x$phase %in% PROVENANCE_PHASES) bad("phase")
    if (!x$status %in% c("completed", "refused", "failed")) bad("status")
    if (!is.logical(x$blinded) || length(x$blinded) != 1L) bad("blinded")
    need(x, c("run", "code", "environment"), is.list, "")
    need(x$run, c("github_run_id", "run_attempt"), is_chr, "run.")
    need(x$run$runner, c("image_os", "image_version", "cpu_model", "kernel"), is_chr,
         "run.runner.")
    need(x$run$runner, "cores", is.numeric, "run.runner.")
    need(x$code, c("package_commit", "release_commit"), is_chr, "code.")
    need(x$code$harness, c("commit", "tree"), is_chr, "code.harness.")
    need(x$code, "declaration", is.list, "code.")
    need(x$environment, c("r_version", "blas"), is_chr, "environment.")
    for (b in c("inputs", "truth", "cells", "chunks", "scoring", "outputs")) {
        if (!b %in% names(x) || !(is.list(x[[b]]) || is.null(x[[b]]))) bad(b)
    }
    # 5.6's nested fields for the phases that run cells and for phase score (I07).
    if (x$phase %in% c("lifecycle", "counted", "correction", "determinism", "amended", "score")) {
        if (!"product_tree_identical" %in% names(x$code)) bad("code.product_tree_identical")
        need(x$code$declaration, c("stage", "commit", "sha256_md", "sha256_json"), is_chr,
             "code.declaration.")
        need(x$code, c("code_files", "helpers_extracted"), is.list, "code.")
    }
    if (x$phase %in% c("lifecycle", "counted", "correction", "determinism", "amended")) {
        need(x$inputs$universe, "sha256", is_chr, "inputs.universe.")
        need(x$inputs$universe, "rows", is.numeric, "inputs.universe.")
        truths <- x$truth[setdiff(names(x$truth), c("blind", "per_unit"))]
        if (!length(truths)) bad("truth")
        for (t1 in truths) need(t1, "sha256", is_chr, "truth.")
        # E4 (I07): the per-unit truths the scoring used, each hashed (where recorded)
        for (u in x$truth$per_unit) {
            for (t1 in u) need(t1, "sha256", is_chr, "truth.per_unit.")
        }
        need(x$truth, "blind", is.list, "truth.")
        if (!is.list(x$cells) || !length(x$cells)) bad("cells")
        for (cl in x$cells) {
            need(cl, c("unit", "arm", "exit"), is_chr, "cells[].")
            need(cl, "seed", is.numeric, "cells[].")
        }
        need(x$scoring, c("metrics_sha256", "comparisons_sha256"), is_chr, "scoring.")
        if (!is.list(x$outputs) || !length(x$outputs)) bad("outputs")
    }
    if (identical(x$phase, "score")) {
        need(x$scoring, "comparisons_sha256", is_chr, "scoring.")
        need(x$scoring, "combined_runs", is.list, "scoring.")
    }
    TRUE
}
