#!/usr/bin/env Rscript
# DIA-003 measurement infrastructure. Sourcing defines helpers only.

#' Validate one required string in a measurement contract.
#' @param value Candidate scalar.
#' @param label Field name for errors.
#' @return The validated value.
#' @noRd
checkpointString <- function(value, label) {
    if (!is.character(value) || length(value) != 1L ||
        is.na(value) || !nzchar(trimws(value))) {
        stop(paste("missing or invalid", label), call. = FALSE)
    }
    value
}

#' Parse the explicit proof CLI without loading MultiScholaR.
#' @param arguments Character command arguments.
#' @return Named option list; does not execute any workload.
#' @noRd
parseCheckpointProofArgs <- function(arguments) {
    values <- c("mode", "route", "contract", "protocol", "library", "output-root", "diagnostics")
    optional <- "checkpoint-codec"
    flags <- c("help", "validate-only")
    result <- list()
    index <- 1L
    while (index <= length(arguments)) {
        token <- arguments[[index]]
        if (!startsWith(token, "--")) stop("unexpected positional argument", call. = FALSE)
        key <- substring(token, 3L)
        if (!key %in% c(values, optional, flags)) stop(paste("unknown option", token), call. = FALSE)
        if (key %in% names(result)) stop(paste("duplicate option", token), call. = FALSE)
        if (key %in% flags) {
            result[[key]] <- TRUE
        } else {
            index <- index + 1L
            if (index > length(arguments) || startsWith(arguments[[index]], "--")) {
                stop(paste("missing value for", token), call. = FALSE)
            }
            result[[key]] <- checkpointString(arguments[[index]], token)
        }
        index <- index + 1L
    }
    if (isTRUE(result$help)) return(result)
    for (key in values) checkpointString(result[[key]], paste0("--", key))
    legal_modes <- c("baseline", "utility", "consumers", "composed",
        "science-iq", "science-limpa", "gui", "gui-cycles", "gui-failures")
    if (!result$mode %in% legal_modes) stop("unsupported mode", call. = FALSE)
    if (!result$route %in% c("iq", "limpa")) stop("unsupported route", call. = FALSE)
    if (!result$diagnostics %in% c("true", "false")) {
        stop("--diagnostics must be true or false", call. = FALSE)
    }
    if (!is.null(result[["checkpoint-codec"]]) &&
        (!identical(result$mode, "utility") ||
        !result[["checkpoint-codec"]] %in% c("rds", "qs2"))) {
        stop("--checkpoint-codec is rds or qs2, for utility mode only", call. = FALSE)
    }
    if (identical(result$mode, "utility") && is.null(result[["checkpoint-codec"]])) {
        result[["checkpoint-codec"]] <- "rds"
    }
    result$diagnostics <- identical(result$diagnostics, "true")
    result
}

#' Validate the frozen scenario coverage before any workload dispatch.
#' @param contract Parsed DIA-002 recipe contract, without vector simplification.
#' @param route Explicit IQ or Limpa route.
#' @return Frozen recipe identities and capacity-only exclusions, not gate acceptance.
#' @noRd
validateCheckpointScenarioContract <- function(contract, route) {
    if (!identical(contract$schema, "multischolar.workflow_checkpoint_recipes.v1") ||
        !identical(contract$version, "1.0") ||
        !identical(contract$baseline, "c7c12851ea4d9e96f91df0e1f9e7b91c2eb51017")) {
        stop("unsupported scenario contract identity")
    }
    if (!route %in% c("iq", "limpa")) stop("unsupported scenario route")
    checkpointString(contract$routes[[route]], "scenario route declaration")
    if (!is.list(contract$recipes) || !length(contract$recipes)) stop("missing scenario recipes")
    ids <- vapply(contract$recipes, \(recipe) checkpointString(recipe$id, "recipe id"), character(1))
    if (anyDuplicated(ids) || !setequal(ids,
        c("small", "representative", "tall", "wide", "annotation_heavy", "cycles", "capacity"))) {
        stop("incomplete or duplicate frozen recipe inventory")
    }
    for (recipe in contract$recipes) {
        for (key in c("n_proteins", "n_samples", "peptides_per_protein", "forms_per_peptide",
            "history_depth", "annotation_chars", "seed", "precursor_rows")) {
            value <- recipe[[key]]
            if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
                !is.finite(value) || value < 1 || value != floor(value)) stop(paste("invalid recipe", key))
        }
        if (recipe$precursor_rows != recipe$n_proteins * recipe$n_samples *
            recipe$peptides_per_protein * recipe$forms_per_peptide) stop("recipe dimensions disagree")
        capacity <- identical(recipe$id, "capacity")
        if (!identical(recipe$frozen_for_timing, !capacity)) stop("changed timing eligibility")
        if (!is.logical(recipe$control_qualified) || length(recipe$control_qualified) != 1L ||
            is.na(recipe$control_qualified)) stop("missing control outcome")
        # FALSE is expected during baseline preparation, not a prerequisite failure.
        if (capacity && !identical(recipe$preflight, "failed_oom_capacity_only")) {
            stop("missing retained capacity failure")
        }
    }
    scenarios <- contract$scenario_matrix
    if (!is.list(scenarios) || length(scenarios) != 26L) stop("incomplete scenario coverage")
    scenario_ids <- vapply(scenarios, \(scenario) checkpointString(scenario$id, "scenario id"), character(1))
    if (anyDuplicated(scenario_ids) || !setequal(scenario_ids, sprintf("C%02d", 1:26))) {
        stop("invalid scenario identities")
    }
    for (scenario in scenarios) {
        members <- unlist(scenario$recipes, use.names = FALSE)
        if (!is.character(members) || !length(members) || anyNA(members) ||
            anyDuplicated(members) || !all(members %in% ids)) stop("invalid scenario recipe references")
        checkpointString(scenario$expected, "scenario expected outcome")
    }
    list(route = route, recipe_ids = ids, timing_recipe_ids = setdiff(ids, "capacity"),
        capacity_only_recipe_ids = "capacity", scenario_ids = scenario_ids,
        baseline_requires_prior_g0 = FALSE, promotion_authorized = FALSE)
}

#' Validate the frozen decision protocol independently of a baseline outcome.
#' @param protocol Parsed protocol list.
#' @return TRUE invisibly; throws on incomplete or changed fixed contracts.
#' @noRd
validateCheckpointProofProtocol <- function(protocol) {
    required_text <- c("schema", "version", "status", "source_control",
        "independence_unit", "pairing_block", "member_order", "paired_summary",
        "decision_rule", "borderline_rule", "non_regression_rule", "retry_rule",
        "invalid_baseline_rule", "missing_metric_rule")
    for (field in required_text) checkpointString(protocol[[field]], field)
    if (!identical(protocol$schema, "multischolar.dia_measurement_protocol.v1") ||
        !identical(protocol$status, "frozen_for_iq_control_preparation") ||
        !identical(protocol$approval$status, "approved")) {
        stop("protocol is not approved and frozen", call. = FALSE)
    }
    checkpointString(protocol$approval$approved_by, "approval.approved_by")
    for (field in c("uncertainty_method", "combined_peak_definition", "run_budget",
        "correction_budget", "service_limits", "g2", "g4", "boundaries", "routes")) {
        if (!is.list(protocol[[field]]) || !length(protocol[[field]])) {
            stop(paste("missing protocol contract", field), call. = FALSE)
        }
    }
    validateCheckpointProofNumbers(protocol)
    validateCheckpointProofProcedure(protocol)
    invisible(TRUE)
}

#' Validate every required nested procedure field before decisions.
#' @param protocol Frozen protocol.
#' @return TRUE invisibly.
#' @noRd
validateCheckpointProofProcedure <- function(protocol) {
    for (field in c("method", "assumptions")) {
        checkpointString(protocol$uncertainty_method[[field]], field)
    }
    if (!identical(protocol$uncertainty_method$method,
        "paired log-ratio Student-t one-sided upper confidence bound")) {
        stop("unsupported uncertainty procedure", call. = FALSE)
    }
    for (field in c("baseline", "histories", "consumers", "cold_state")) {
        checkpointString(protocol$boundaries[[field]], paste("boundary", field))
    }
    budgets <- c(qualified_blocks_per_build_recipe_route = 7,
        diagnostic_profiles_per_case = 1, fixture_preflights_per_recipe = 1)
    for (field in names(budgets)) {
        value <- protocol$run_budget[[field]]
        if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
            value != budgets[[field]]) stop(paste("invalid run budget", field), call. = FALSE)
    }
    value <- protocol$correction_budget$design_corrections_per_failed_pilot
    if (!is.numeric(value) || length(value) != 1L || is.na(value) || value != 1) {
        stop("invalid correction budget", call. = FALSE)
    }
    checkpointString(protocol$correction_budget$requires, "correction requirements")
    checkpointString(protocol$g2$retained_heap_measurement, "retained heap measurement")
    checkpointString(protocol$g2$scope, "G2 scope")
    checkpointString(protocol$g4$required, "G4 completion boundary")
    for (route in c("iq", "limpa")) checkpointString(protocol$routes[[route]]$status, route)
    if (!identical(protocol$routes$limpa$inherits_iq_evidence, FALSE)) {
        stop("Limpa cannot inherit IQ evidence", call. = FALSE)
    }
    invisible(TRUE)
}

#' Validate fixed numerical limits and complete resource-accounting semantics.
#' @param protocol Parsed protocol list.
#' @return TRUE invisibly.
#' @noRd
validateCheckpointProofNumbers <- function(protocol) {
    limits <- protocol$service_limits
    expected <- c(memory_max_bytes = 8 * 1024^3, memory_swap_max_bytes = 0,
        runtime_max_seconds = 900, optional_cache_max_bytes = 64 * 1024^2)
    for (field in names(expected)) {
        value <- limits[[field]]
        if (!is.numeric(value) || length(value) != 1L ||
            is.na(value) || value != expected[[field]]) {
            stop(paste("changed or missing resource limit", field), call. = FALSE)
        }
    }
    frozen <- list(g2 = c(retained_heap_ratio_max = 0.70,
        sequence_elapsed_ratio_max = 1.10, sequence_cpu_ratio_max = 1.10,
        save_max_seconds = 30, cold_revert_max_seconds = 30),
        g4 = c(elapsed_ratio_max = 0.90, combined_peak_ratio_max = 0.90,
            immediate_comparator_ratio_max = 1.00))
    for (gate in names(frozen)) for (field in names(frozen[[gate]])) {
        value <- protocol[[gate]][[field]]
        if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
            value != frozen[[gate]][[field]]) {
            stop(paste("changed or missing frozen threshold", gate, field), call. = FALSE)
        }
    }
    if (!is.numeric(protocol$repetitions) || length(protocol$repetitions) != 1L ||
        is.na(protocol$repetitions) || protocol$repetitions != 7 ||
        !identical(as.numeric(protocol$uncertainty_method$df), 6) ||
        !identical(as.numeric(protocol$uncertainty_method$confidence), 0.95)) {
        stop("unsupported or missing paired uncertainty contract", call. = FALSE)
    }
    for (field in c("metric", "units", "membership", "interval", "shared_memory",
        "file_cache", "driver", "exclusions")) {
        checkpointString(protocol$combined_peak_definition[[field]], field)
    }
    if (!identical(protocol$combined_peak_definition$metric, "fresh service cgroup v2 memory.peak") ||
        !identical(protocol$combined_peak_definition$units, "bytes")) {
        stop("unsupported combined-peak calculation", call. = FALSE)
    }
    invisible(TRUE)
}

#' Check paired receipt identities before looking at performance ratios.
#' @param receipts Data frame of separate independent control/candidate runs.
#' @param expected Bound route/input/protocol/boundary/dependency and build identities.
#' @param repetitions Required independent pair count.
#' @return NULL if structurally bound, otherwise a held-reason string.
#' @noRd
checkpointPairBindingError <- function(receipts, expected, repetitions) {
    bindings <- c("route", "input_hash", "protocol_hash", "boundary_hash", "dependency_hash")
    columns <- c(bindings, "build_hash", "pair_id", "run_id", "process_instance", "member",
        "status", "science_ok", "outputs_ok", "diagnostics", "elapsed_seconds",
        "combined_peak_bytes", "swap_peak_bytes")
    if (!is.data.frame(receipts) || !all(columns %in% names(receipts))) {
        return("missing receipt fields")
    }
    if (nrow(receipts) != 2 * repetitions || anyNA(receipts[columns])) {
        return("incomplete pairs or missing receipt fields")
    }
    for (field in c(bindings, "control_build", "candidate_build")) {
        if (!is.character(expected[[field]]) || length(expected[[field]]) != 1L ||
            is.na(expected[[field]]) || !nzchar(expected[[field]])) {
            return(paste("missing expected binding", field))
        }
    }
    if (!expected$route %in% c("iq", "limpa")) return("unsupported expected route")
    if (identical(expected$route, "limpa")) {
        return("Limpa evidence qualification not implemented before DIA-017")
    }
    for (field in bindings) {
        if (!is.character(receipts[[field]]) || any(receipts[[field]] != expected[[field]])) {
            return(paste("mismatched binding", field))
        }
    }
    for (field in c("pair_id", "run_id", "process_instance", "build_hash", "member")) {
        if (!is.character(receipts[[field]]) || any(!nzchar(receipts[[field]]))) {
            return(paste("invalid receipt identity", field))
        }
    }
    if (anyDuplicated(receipts$run_id) || anyDuplicated(receipts$process_instance)) {
        return("reused run or process instance: cycles are not independent replicates")
    }
    if (!all(receipts$member %in% c("control", "candidate")) ||
        length(unique(receipts$pair_id)) != repetitions) return("invalid pair membership")
    counts <- table(receipts$pair_id, receipts$member)
    if (ncol(counts) != 2L || any(counts != 1L)) return("unbalanced pair membership")
    builds <- ifelse(receipts$member == "control", expected$control_build, expected$candidate_build)
    if (any(receipts$build_hash != builds)) return("mismatched build")
    for (field in c("science_ok", "outputs_ok", "diagnostics")) {
        if (!is.logical(receipts[[field]])) return(paste("invalid flag", field))
    }
    NULL
}

#' Calculate the frozen paired log-ratio summary.
#' @param candidate,control Positive finite paired measurements.
#' @param confidence One-sided confidence level.
#' @return Point ratio and one-sided upper bound.
#' @noRd
checkpointPairedRatio <- function(candidate, control, confidence) {
    ratios <- candidate / control
    logs <- log(ratios)
    # Preserve exact uniform boundaries without an arbitrary tolerance margin.
    point <- if (length(unique(ratios)) == 1L) ratios[[1L]] else exp(mean(logs))
    upper <- point * exp(stats::qt(confidence, df = length(logs) - 1L) *
        stats::sd(logs) / sqrt(length(logs)))
    list(point = point, upper = upper)
}

#' Evaluate synthetic or verified paired receipts without inferring gate promotion.
#' @param receipts Independent-run receipt data frame.
#' @param protocol Frozen measurement protocol.
#' @param expected Expected identity bindings, including control/candidate builds.
#' @param gate Either g2_history or g4.
#' @param comparison Original control or matched accepted comparator.
#' @return Evidence decision, reasons and ratios; not a signed gate receipt.
#' @noRd
evaluateCheckpointPairs <- function(receipts, protocol, expected,
    gate = c("g2_history", "g4"), comparison = c("original_control", "accepted_comparator")) {
    gate <- match.arg(gate)
    comparison <- match.arg(comparison)
    validateCheckpointProofProtocol(protocol)
    result <- function(outcome, reason, ratios = list()) {
        list(outcome = outcome, reason = reason, ratios = ratios, promotion_authorized = FALSE)
    }
    error <- checkpointPairBindingError(receipts, expected, protocol$repetitions)
    if (!is.null(error)) return(result("held", error))
    if (any(receipts$diagnostics)) return(result("held", "profiled member in timing evidence"))
    controls <- receipts[receipts$member == "control", , drop = FALSE]
    candidates <- receipts[receipts$member == "candidate", , drop = FALSE]
    candidates <- candidates[match(controls$pair_id, candidates$pair_id), , drop = FALSE]
    if (any(controls$status != "complete" | !controls$science_ok | !controls$outputs_ok)) {
        return(result("held", "invalid control; no finite speedup denominator"))
    }
    if (any(candidates$status %in% c("queued", "running", "outputs_pending"))) {
        return(result("held", "candidate is not terminal; no performance decision yet"))
    }
    if (any(candidates$status != "complete" | !candidates$science_ok | !candidates$outputs_ok)) {
        return(result("failed", "candidate execution, science or required outputs failed"))
    }
    metrics <- if (gate == "g4") c("elapsed_seconds", "combined_peak_bytes") else
        c("retained_heap_bytes", "elapsed_seconds", "cpu_seconds")
    required <- unique(c(metrics, "elapsed_seconds", "cpu_seconds", "combined_peak_bytes", "swap_peak_bytes",
        if (gate == "g2_history") c("max_save_seconds", "max_cold_revert_seconds")))
    for (field in required) {
        value <- receipts[[field]]
        if (!is.numeric(value) || length(value) != nrow(receipts) ||
            any(!is.finite(value)) || any(value < 0) ||
            (field %in% metrics && any(value == 0))) {
            return(result("held", paste("missing or invalid metric", field)))
        }
    }
    invalid <- receipts$combined_peak_bytes > protocol$service_limits$memory_max_bytes |
        receipts$elapsed_seconds > protocol$service_limits$runtime_max_seconds |
        receipts$swap_peak_bytes != 0
    if (any(invalid & receipts$member == "control")) return(result("held", "control limit violation"))
    if (any(invalid)) return(result("failed", "candidate limit violation"))
    if (gate == "g2_history" && any(
        candidates$max_save_seconds > protocol$g2$save_max_seconds |
        candidates$max_cold_revert_seconds > protocol$g2$cold_revert_max_seconds)) {
        return(result("failed", "G2 save or cold-revert latency ceiling exceeded"))
    }
    thresholds <- if (gate == "g2_history") c(protocol$g2$retained_heap_ratio_max,
        protocol$g2$sequence_elapsed_ratio_max, protocol$g2$sequence_cpu_ratio_max) else
        if (comparison == "accepted_comparator") rep(protocol$g4$immediate_comparator_ratio_max, 2) else
            c(protocol$g4$elapsed_ratio_max, protocol$g4$combined_peak_ratio_max)
    ratios <- stats::setNames(lapply(seq_along(metrics), \(i) {
        estimate <- checkpointPairedRatio(candidates[[metrics[[i]]]], controls[[metrics[[i]]]],
            protocol$uncertainty_method$confidence)
        c(estimate, list(threshold = thresholds[[i]]))
    }), metrics)
    if (any(vapply(ratios, \(x) x$point > x$threshold, logical(1)))) {
        return(result("failed", "point performance threshold missed", ratios))
    }
    if (any(vapply(ratios, \(x) x$upper > x$threshold, logical(1)))) {
        return(result("held", "uncertainty crosses threshold", ratios))
    }
    result("accepted_evidence", "paired thresholds met; separate gate verification required", ratios)
}


#' Check an explicit filesystem directory without creating it.
#' @param path Requested directory.
#' @param label Option name.
#' @return Canonical existing directory.
#' @noRd
checkpointExistingDirectory <- function(path, label) {
    checkpointString(path, label)
    if (!dir.exists(path)) stop(paste(label, "directory does not exist"), call. = FALSE)
    normalizePath(path, winslash = "/", mustWork = TRUE)
}

#' Validate output storage, not a scientific or performance gate.
#' @param output_root Existing dedicated directory.
#' @param source_root Package source to protect.
#' @return Canonical disk-backed path.
#' @noRd
validateCheckpointOutputRoot <- function(output_root, source_root) {
    output_root <- checkpointExistingDirectory(output_root, "output root")
    source_root <- checkpointExistingDirectory(source_root, "source root")
    protected <- unique(c("/", normalizePath(path.expand("~")), source_root, dirname(source_root)))
    if (output_root %in% protected || startsWith(output_root, paste0(source_root, "/"))) {
        stop("output root must be a dedicated directory outside package source", call. = FALSE)
    }
    if (!identical(Sys.info()[["sysname"]], "Linux")) {
        stop("resource-qualified runner currently requires Linux", call. = FALSE)
    }
    filesystem <- system2("stat", c("-f", "-c", "%T", "--", shQuote(output_root)),
        stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(filesystem, "status")) || length(filesystem) != 1L ||
        !filesystem %in% c("btrfs", "ext2/ext3", "xfs", "zfs")) {
        stop("output root is not verified disk-backed storage", call. = FALSE)
    }
    output_root
}

#' Inspect an explicitly installed package without loading or prewarming it.
#' @param library_path Requested library directory.
#' @return Structural installation metadata; load_verified is always FALSE.
#' @noRd
inspectCheckpointLibrary <- function(library_path) {
    library_path <- checkpointExistingDirectory(library_path, "library")
    package_path <- file.path(library_path, "MultiScholaR")
    required <- c("DESCRIPTION", "NAMESPACE", "Meta/package.rds",
        "R/MultiScholaR.rdb", "R/MultiScholaR.rdx")
    paths <- file.path(package_path, required)
    if (!all(file.exists(paths)) || any(file.info(paths)$isdir) ||
        any(is.na(file.size(paths))) || any(file.size(paths) <= 0)) {
        stop("explicit library does not contain a complete installed MultiScholaR", call. = FALSE)
    }
    description <- read.dcf(file.path(package_path, "DESCRIPTION"))
    if (nrow(description) != 1L || !all(c("Package", "Version", "Built") %in% colnames(description)) ||
        !identical(unname(description[1L, "Package"]), "MultiScholaR")) {
        stop("library metadata is not an installed MultiScholaR build", call. = FALSE)
    }
    native_root <- file.path(package_path, "libs")
    native_paths <- sort(list.files(native_root, pattern = "[.](so|dll|dylib)$",
        recursive = TRUE, full.names = TRUE))
    native_hashes <- stats::setNames(lapply(native_paths, function(path)
        digest::digest(file = path, algo = "sha256")), substring(native_paths, nchar(native_root) + 2L))
    if ("NeedsCompilation" %in% colnames(description) &&
        identical(unname(description[1L, "NeedsCompilation"]), "yes") && !length(native_hashes))
        stop("compiled package native library is missing")
    list(library_path = library_path, native_libraries = native_hashes,
        package_path = normalizePath(package_path, winslash = "/", mustWork = TRUE),
        version = unname(description[1L, "Version"]),
        dia_output_policy = if ("Config/MultiScholaR/DIA-output-policy" %in% colnames(description))
            unname(description[1L, "Config/MultiScholaR/DIA-output-policy"]) else "legacy_s4",
        built = unname(description[1L, "Built"]), load_verified = FALSE)
}

#' Resolve every baseline input from existing preparation receipts, without hydrating data.
#' @param contract Frozen recipe contract.
#' @param output_root Dedicated proof root containing existing preparation receipts.
#' @return Named input bindings for all timing-eligible recipes.
#' @noRd
resolveCheckpointBaselineInputs <- function(contract, output_root) {
    inventory <- validateCheckpointScenarioContract(contract, "iq")
    root <- normalizePath(output_root, mustWork = TRUE)
    prepared_dirs <- list.files(root, "^browser-smoke\\.[A-Za-z0-9]{6}$", full.names = TRUE)
    prepared_receipts <- file.path(prepared_dirs, "inputs.json")
    prepared_receipts <- prepared_receipts[file.exists(prepared_receipts)]
    result <- list()
    for (id in inventory$timing_recipe_ids) {
        direct <- file.path(root, id, "inputs.json")
        candidates <- c(if (file.exists(direct)) direct else character(), prepared_receipts)
        candidates <- candidates[vapply(candidates, \(path) {
            size <- file.size(path)
            if (is.na(size) || size <= 0 || size > 1048576) {
                stop(paste("missing or oversized prepared input receipt:", path))
            }
            identical(jsonlite::read_json(path)$recipe$id, id)
        }, logical(1))]
        if (length(candidates) != 1L) {
            stop(paste("missing or ambiguous prepared input receipt:", id))
        }
        receipt_path <- candidates[[1L]]
        receipt <- jsonlite::read_json(receipt_path)
        recipe <- contract$recipes[[match(id, inventory$recipe_ids)]]
        if (!isTRUE(all.equal(receipt$recipe, recipe, check.attributes = TRUE)) ||
            !receipt$kind %in% c("small_fixture_preparation_not_timing", "frozen_report_companions_not_timing")) {
            stop(paste("prepared inputs do not match frozen recipe:", id))
        }
        files <- receipt$files
        hashes <- receipt$sha256
        if (identical(receipt$kind, "frozen_report_companions_not_timing")) {
            files$report <- receipt$report
            hashes$report <- receipt$report_sha256
        }
        for (key in c("report", "fasta", "design")) {
            checkpointString(files[[key]], paste(id, key, "path"))
            checkpointString(hashes[[key]], paste(id, key, "fingerprint"))
            if (!grepl("^[0-9a-f]{64}$", hashes[[key]])) stop("invalid input fingerprint")
            path <- normalizePath(files[[key]], mustWork = TRUE)
            if (!startsWith(path, paste0(root, "/")) || dir.exists(path) ||
                is.na(file.size(path)) || file.size(path) <= 0) stop("input is not a nonempty owned file")
            files[[key]] <- path
        }
        result[[id]] <- list(recipe = recipe, files = files[c("report", "fasta", "design")],
            sha256 = hashes[c("report", "fasta", "design")],
            preparation_receipt = normalizePath(receipt_path),
            preparation_sha256 = digest::digest(file = receipt_path, algo = "sha256"))
    }
    result
}

#' Dispatch frozen baseline members serially and stop on the first failed member.
#' @param inputs Complete prepared recipe bindings.
#' @param protocol Frozen protocol.
#' @param run_member Existing guarded service adapter, called once per member.
#' @param record_member Persist each member result before another may run.
#' @return Completed member records; never a scientific/performance gate decision.
#' @noRd
dispatchCheckpointBaseline <- function(inputs, protocol, run_member, record_member) {
    validateCheckpointProofProtocol(protocol)
    expected <- c("small", "representative", "tall", "wide", "annotation_heavy", "cycles")
    if (!identical(names(inputs), expected)) stop("baseline dispatch requires the full frozen input inventory")
    results <- list()
    run_ids <- character()
    for (id in names(inputs)) {
        for (block in seq_len(protocol$repetitions)) {
            result <- tryCatch(run_member(inputs[[id]], block), error = function(error) {
                list(status = "failed", error = conditionMessage(error))
            })
            # Record invalid/failing responses too; no second member after a failure.
            if (!is.list(result)) result <- list(status = "failed", error = "invalid member response")
            valid <- identical(result$status, "complete_verified") &&
                is.character(result$run_id) && length(result$run_id) == 1L &&
                !is.na(result$run_id) && nzchar(result$run_id) && !result$run_id %in% run_ids
            if (!valid && identical(result$status, "complete_verified")) {
                result$status <- "failed"
                result$error <- "missing or reused independent run identity"
            }
            record_member(id, block, result)
            if (!valid) stop(paste("baseline member did not complete:", id, block))
            run_ids <- c(run_ids, result$run_id)
            results[[paste(id, block, sep = ":")]] <- result
        }
    }
    results
}

#' Invoke the existing guarded service and resolve its exact owned output directory.
#' @param operation Existing shell service operation.
#' @param library_path Explicit installed library.
#' @param output_root Dedicated proof root.
#' @param arguments Operation-specific paths.
#' @return Service exit status and owned member directory, never completion acceptance.
#' @noRd
runCheckpointService <- function(operation, library_path, output_root, arguments) {
    runner <- file.path(checkpointRunnerDirectory(), "run_workflow_checkpoint_proof.sh")
    lines <- suppressWarnings(system2("bash", shQuote(c(runner, "--service", operation,
        library_path, output_root, arguments)), stdout = TRUE, stderr = TRUE))
    status <- attr(lines, "status")
    if (is.null(status)) status <- 0L
    units <- sub("^unit=", "", grep("^unit=multischolar-dia003-browser-smoke\\.[A-Za-z0-9]{6}$", lines, value = TRUE))
    launches <- grep('^\\{"unit":', lines, value = TRUE)
    for (launch in launches) {
        unit <- jsonlite::fromJSON(launch)$unit
        if (is.character(unit) && length(unit) == 1L) units <- c(units, sub("\\.service$", "", unit))
    }
    units <- unique(units)
    if (length(units) != 1L || !grepl("^multischolar-dia003-browser-smoke\\.[A-Za-z0-9]{6}$", units)) {
        stop(paste("service returned no unique owned run identity; exit", status))
    }
    directory <- normalizePath(file.path(output_root, sub("^multischolar-dia003-", "", units)), mustWork = TRUE)
    if (!startsWith(directory, paste0(normalizePath(output_root), "/"))) stop("service output escaped proof root")
    list(exit_status = as.integer(status), member_dir = directory, unit = paste0(units, ".service"))
}

#' Bind a browser service and its separate output verification to one baseline request.
#' @param input Frozen prepared inputs for one recipe.
#' @param block Independent block index.
#' @param library_path Explicit installed package library.
#' @param output_root Dedicated proof root.
#' @param service Existing service invoker; injectable for orchestration tests only.
#' @return Bound evidence; held until the complete baseline boundary is implemented.
#' @noRd
runCheckpointBaselineMember <- function(input, block, library_path, output_root,
    service = runCheckpointService) {
    result <- list(status = "running", recipe = input$recipe$id, block = block,
        route = "iq", input_sha256 = input$sha256, promotion_authorized = FALSE)
    tryCatch({
        member <- service("browser-baseline", library_path, output_root,
            unlist(input$files[c("report", "fasta", "design")], use.names = FALSE))
        result$member_dir <- member$member_dir
        result$service_exit <- member$exit_status
        if (member$exit_status != 0L) stop("browser service or controller failed; do not verify or retry")
        terminal <- jsonlite::read_json(file.path(member$member_dir, "terminal-interval.json"))
        if (!isTRUE(terminal$service_succeeded) || !identical(terminal$unit, member$unit) ||
            !identical(terminal$status, "service_tree_interval_only") ||
            !is.numeric(terminal$elapsed_seconds) || length(terminal$elapsed_seconds) != 1L ||
            !is.finite(terminal$elapsed_seconds) || terminal$elapsed_seconds <= 0) {
            stop("missing or mismatched successful service-tree completion")
        }
        result$run_id <- checkpointString(terminal$invocation_id, "service invocation")
        result$terminal <- terminal
        browser <- jsonlite::read_json(file.path(member$member_dir, "smoke-receipt.json"))
        if (!identical(browser$status, "passed_summary_ui_only") ||
            !identical(browser$installation$library_path, normalizePath(library_path)) ||
            !identical(browser$input_sha256, input$sha256) ||
            !isTRUE(browser$inputs$qc_defaults) || !identical(browser$diagnostics, FALSE)) {
            stop("browser build/input/settings binding differs from requested baseline")
        }
        verification <- service("verify-da", library_path, output_root, member$member_dir)
        result$verification_dir <- verification$member_dir
        result$verification_exit <- verification$exit_status
        if (verification$exit_status != 0L) stop("independent output verification failed")
        receipt <- jsonlite::read_json(file.path(verification$member_dir, "verification-receipt.json"))
        if (!identical(receipt$status, "passed_small_summary_file_checks") ||
            !identical(receipt$member_dir, member$member_dir)) stop("output verification receipt binding differs")
        if (!identical(receipt$required_iq_qc$status, "passed_iq_qc_presence_only") ||
            !identical(receipt$summary_files$status, "passed_small_summary_files") ||
            !all(c("Volcano_Plots", "NumSigDaMolecules", "Heatmap") %in%
                names(receipt$summary_files$publication_figures))) {
            stop("required baseline output inventory has not been verified")
        }
        result$verified_scope <- "default IQ GUI, full-history session, selected QC/DA/summary outputs and restore"
        result$status <- "complete_verified"
        result
    }, error = function(error) {
        result$status <- "failed"
        result$error <- conditionMessage(error)
        result
    })
}

#' Decode one required nonnegative terminal service counter.
#' @param terminal Parsed retained-service interval receipt.
#' @param name Raw systemd property name.
#' @return Exact numeric value within R's safe integer range.
#' @noRd
checkpointTerminalCounter <- function(terminal, name) {
    value <- terminal$raw_properties[[name]]
    if (!is.character(value) || length(value) != 1L ||
        !grepl("^[0-9]+$", value)) stop(paste("missing terminal counter", name))
    value <- as.numeric(value)
    if (!is.finite(value) || value > 2^53 - 1) stop(paste("unsafe terminal counter", name))
    value
}

#' Bind one utility service member to its separate verifier and accounting.
#' @param recipe Frozen recipe.
#' @param block Independent pair index.
#' @param treatment resident, saved (native RDS) or qs2.
#' @param variant independent, shared_metadata or retained_alias.
#' @param library_path,output_root,contract_path Bound paths.
#' @param service Existing guarded service adapter.
#' @param input_dir Separately prepared immutable input directory.
#' @return Complete scalar member evidence for paired evaluation.
#' @noRd
runCheckpointUtilityServiceMember <- function(recipe, block, treatment, variant,
    library_path, output_root, contract_path, service = runCheckpointService,
    input_dir = NULL) {
    operation <- service("utility-member", library_path, output_root,
        c(contract_path, recipe$id, treatment, variant, as.character(block), input_dir))
    if (operation$exit_status != 0L) stop("utility service or controller failed")
    terminal <- jsonlite::read_json(file.path(operation$member_dir, "terminal-interval.json"))
    member <- jsonlite::read_json(file.path(operation$member_dir, "member-receipt.json"))
    if (isTRUE(member$diagnostics)) stop("diagnostic utility member cannot qualify performance")
    if (!isTRUE(terminal$service_succeeded) ||
        !identical(terminal$status, "service_tree_interval_only") ||
        !identical(member$status, "complete") || !isTRUE(member$science_ok) ||
        !isTRUE(member$outputs_ok) || !identical(member$recipe, recipe$id) ||
        !identical(member$treatment, treatment) || !identical(member$variant, variant) ||
        !identical(member$block, as.integer(block))) stop("utility member binding differs")
    if (!is.null(input_dir) && (isTRUE(member$fixture_generation_in_timed_service) ||
        !identical(member$prepared_inputs_sha256,
            digest::digest(file = file.path(input_dir, "utility-inputs.json"), algo = "sha256")))) {
        stop("utility prepared-input boundary differs")
    }
    verification <- service("utility-verify", library_path, output_root, operation$member_dir)
    if (verification$exit_status != 0L) stop("utility member verification failed")
    verified <- jsonlite::read_json(file.path(verification$member_dir,
        "utility-verification.json"))
    if (!identical(verified$status, "passed_utility_member") ||
        !identical(verified$member_dir, operation$member_dir)) {
        stop("utility verifier binding differs")
    }
    if (!is.null(input_dir) && (is.null(verified$legacy_sha256) ||
        !identical(verified$legacy_sha256, member$legacy_export$sha256))) {
        stop("utility verified legacy file binding differs")
    }
    installation <- inspectCheckpointLibrary(library_path)
    artifacts <- file.path(installation$package_path,
        c("DESCRIPTION", "NAMESPACE", "R/MultiScholaR.rdb", "R/MultiScholaR.rdx"))
    artifact_hashes <- vapply(artifacts, digest::digest, character(1),
        file = TRUE, algo = "sha256")
    build_hash <- digest::digest(list(artifact_hashes, member$source_sha256, treatment),
        algo = "sha256")
    description <- read.dcf(artifacts[[1L]])
    dependency_fields <- intersect(c("Depends", "Imports", "Suggests", "LinkingTo"),
        colnames(description))
    dependency_hash <- digest::digest(description[1L, dependency_fields, drop = FALSE],
        algo = "sha256")
    list(route = "iq", recipe = recipe$id, variant = variant,
        pair_id = paste(recipe$id, variant, sprintf("%02d", block), sep = ":"),
        run_id = terminal$invocation_id,
        process_instance = paste(terminal$invocation_id, member$process_id, sep = ":"),
        treatment = treatment,
        member = if (identical(treatment, "resident")) "control" else "candidate",
        build_hash = build_hash,
        input_hash = digest::digest(list(recipe, variant, member$prepared_inputs_sha256),
            algo = "sha256"),
        dependency_hash = dependency_hash,
        boundary_hash = digest::digest(list(
            "fresh installed utility process through full legacy export and teardown",
            fixture_generation_in_timed_service = member$fixture_generation_in_timed_service),
            algo = "sha256"),
        status = "complete", science_ok = TRUE, outputs_ok = TRUE, diagnostics = FALSE,
        elapsed_seconds = terminal$elapsed_seconds,
        cpu_seconds = checkpointTerminalCounter(terminal, "CPUUsageNSec") / 1e9,
        combined_peak_bytes = checkpointTerminalCounter(terminal, "MemoryPeak"),
        swap_peak_bytes = checkpointTerminalCounter(terminal, "MemorySwapPeak"),
        retained_heap_bytes = member$retained_heap_bytes,
        max_save_seconds = member$max_save_seconds,
        max_cold_revert_seconds = member$max_cold_revert_seconds,
        io_read_bytes = checkpointTerminalCounter(terminal, "IOReadBytes"),
        io_write_bytes = checkpointTerminalCounter(terminal, "IOWriteBytes"),
        payload_read_bytes = member$payload_read_bytes,
        export_heap_bytes = member$export_heap_bytes,
        prepared_inputs_sha256 = member$prepared_inputs_sha256,
        legacy_sha256 = verified$legacy_sha256,
        member_dir = operation$member_dir, verification_dir = verification$member_dir,
        protocol_hash = NA_character_, promotion_authorized = FALSE)
}

#' Convert scalar utility member records to the frozen evaluator frame.
#' @param records Complete member records.
#' @param protocol_hash Frozen protocol fingerprint.
#' @return Data frame with no nested evidence fields.
#' @noRd
checkpointUtilityFrame <- function(records, protocol_hash) {
    fields <- c("route", "input_hash", "protocol_hash", "boundary_hash",
        "dependency_hash", "build_hash", "pair_id", "run_id", "process_instance",
        "member", "status", "science_ok", "outputs_ok", "diagnostics",
        "elapsed_seconds", "cpu_seconds", "combined_peak_bytes", "swap_peak_bytes",
        "retained_heap_bytes", "max_save_seconds", "max_cold_revert_seconds")
    rows <- lapply(records, function(record) {
        record$protocol_hash <- protocol_hash
        as.data.frame(record[fields], stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}

#' Execute the full frozen utility matrix serially.
#' @param recipes Timing-eligible frozen recipes.
#' @param protocol Frozen measurement protocol.
#' @param run_member Member callback.
#' @param record_member Durable per-member callback.
#' @param candidate Explicit saved (native RDS) or qs2 treatment.
#' @return All member records and per-recipe/variant evidence decisions.
#' @noRd
dispatchCheckpointUtility <- function(recipes, protocol, run_member, record_member,
    candidate = "saved") {
    validateCheckpointProofProtocol(protocol)
    if (!candidate %in% c("saved", "qs2")) stop("invalid utility candidate")
    expected <- c("small", "representative", "tall", "wide", "annotation_heavy", "cycles")
    if (!identical(names(recipes), expected)) stop("utility dispatch requires full frozen recipes")
    variants <- c("independent", "shared_metadata", "retained_alias")
    records <- list()
    for (recipe in recipes) for (variant in variants) {
        # DIA-009 revision 1.6: the small path stays resident. This is an explicit
        # frozen-fixture policy, not a production size cutoff or a passed qs2 case.
        if (identical(candidate, "qs2") && identical(recipe$id, "small")) next
        for (block in seq_len(protocol$repetitions)) {
            order <- if (block %% 2L) c("resident", candidate) else c(candidate, "resident")
            for (treatment in order) {
                result <- tryCatch(run_member(recipe, block, treatment, variant),
                    error = function(error) list(status = "failed",
                        error = conditionMessage(error), recipe = recipe$id,
                        block = block, treatment = treatment, variant = variant))
                record_member(result)
                key <- paste(recipe$id, variant, block, treatment, sep = ":")
                records[[key]] <- result
                if (!identical(result$status, "complete")) {
                    stop(paste("utility member did not complete:", key))
                }
            }
            pair <- records[paste(recipe$id, variant, block, order, sep = ":")]
            if (identical(candidate, "qs2") &&
                (is.null(pair[[1L]]$legacy_sha256) ||
                is.null(pair[[1L]]$prepared_inputs_sha256) ||
                !identical(pair[[1L]]$legacy_sha256, pair[[2L]]$legacy_sha256) ||
                !identical(pair[[1L]]$prepared_inputs_sha256, pair[[2L]]$prepared_inputs_sha256))) {
                stop("utility paired legacy outputs or prepared inputs differ; retain both members")
            }
        }
    }
    records
}

#' Display the proof command contract.
#' @return NULL invisibly.
#' @noRd
checkpointProofHelp <- function() {
    cat(paste(
        "DIA installed-package proof harness",
        "Required: --mode --route --contract --protocol --library --output-root --diagnostics",
        "Flags: --help --validate-only",
        "Utility option: --checkpoint-codec rds|qs2 (default rds; explicit qs2 qualification)",
        "qs2 qualification covers the five frozen large recipes; small histories stay resident.",
        "Legal routes: iq, limpa. No default route or library.",
        "Recorded pilot: --service recorded-prepare LIB ROOT SESSION [--validate-only]",
        "--service recorded-member LIB ROOT PREPARED resident|saved [--validate-only]",
        "--service recorded-timing LIB ROOT PREPARED resident|saved|qs2 [--validate-only]",
        "--service recorded-qc LIB ROOT PREPARED [--validate-only]",
        "--service recorded-verify LIB ROOT MEMBER [--validate-only]",
        "Incomplete drivers refuse execution; no fallback to source loading or fake science.",
        sep = "\n"), "\n")
    invisible(NULL)
}


#' Resolve this runner's own directory, independent of the request working directory.
#' @return Canonical tools/profiling directory.
#' @noRd
checkpointRunnerDirectory <- function() {
    files <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
    if (length(files) != 1L) stop("cannot resolve runner location")
    dirname(normalizePath(files, mustWork = TRUE))
}

#' Read bounded launch intent without executing its command or inferring liveness.
#' @param path Request-owned NUL-delimited controller record.
#' @return Decoded identity/limits/argv; this is not a completion or live-job record.
#' @noRd
readCheckpointControllerRequest <- function(path) {
    size <- file.info(path)$size
    if (length(size) != 1L || is.na(size) || size < 1 || size > 1048576) {
        stop("missing or oversized controller request")
    }
    bytes <- readBin(path, "raw", n = size)
    if (!identical(utils::tail(bytes, 1L), as.raw(0))) stop("unterminated controller request")
    connection <- rawConnection(bytes)
    on.exit(close(connection), add = TRUE)
    fields <- readBin(connection, character(), n = 4096L)
    if (length(readBin(connection, "raw", n = 1L)) || length(fields) < 8L) {
        stop("invalid controller request field count")
    }
    if (!identical(fields[[1L]], "dia-controller-request-v1") ||
        !grepl("^multischolar-dia003-browser-smoke\\.[A-Za-z0-9]{6}\\.service$", fields[[2L]]) ||
        !fields[[3L]] %in% c("browser-smoke", "prepare-small", "install-corrected", "browser-import",
            "browser-design", "browser-peptides", "browser-proteins", "browser-normalization",
            "browser-da", "browser-da-accounted", "browser-da-default-qc", "browser-da-profile", "browser-baseline", "verify-normalization", "verify-da", "utility-prepare", "utility-member", "utility-verify",
            "recorded-prepare", "recorded-member", "recorded-timing", "recorded-qc", "recorded-verify",
            "reader-build", "reader-check", "reader-prepare", "reader-inspect", "reader-member", "reader-verify",
            "accounting-probe", "accounting-probe-failure", "accounting-probe-retained", "prepare-companions", "exclusion-benchmark", "exclusion-memory", "exclusion-allocation") ||
        !nzchar(fields[[4L]]) || !nzchar(fields[[8L]]) ||
        !identical(fields[5:7], if (fields[[3L]] == "reader-check" &&
            "--internal-reader-check" %in% fields[-seq_len(7L)] &&
            grepl("^perf107:(B|C|D|E|H|HF|S|L1|L2)$", utils::tail(fields, 1L)))
            c("17179869184", "0", switch(sub("^perf107:", "",
                utils::tail(fields, 1L)), B = "600", S = "600", C = "300", D = "300", "900")) else if (fields[[3L]] == "reader-check" &&
            "--internal-reader-check" %in% fields[-seq_len(7L)] &&
            grepl("^gui-bulk-rollup:(qualification-preflight|qualification-prepare|writer-build|writer-opt-(build|preflight|diagnostic:(primary|contrast))|writer-r-preflight|writer-r-readback:[A-Za-z0-9]{6}|writer-r:(screen|full):(primary|contrast):P1[456]:[12]|io-preflight|io-prepare|io:(duckdb|fastduckdb):[12]|scan-inspection|crossover:(B|P0):(4608|9600|96000|960000):[12]|contrast:(B|P([0-9]|1[0-4]))|diagnostic|diagnostic-preflight|preflight|relational-preflight|materialize-preflight|adaptation-preflight|fast-p3-preflight|writer-preflight|screen:([ABCD]|P([0-9]|1[0-4])):[12]|full:([ABCD]|P([0-9]|1[0-4])))$", utils::tail(fields, 1L)))
            c("17179869184", "0", if (grepl(":(full|contrast):", utils::tail(fields, 1L)) || utils::tail(fields, 1L) %in% c("gui-bulk-rollup:qualification-prepare", "gui-bulk-rollup:writer-build", "gui-bulk-rollup:writer-opt-build")) "1800" else if ((identical(utils::tail(fields, 1L), "gui-bulk-rollup:diagnostic") || startsWith(utils::tail(fields, 1L), "gui-bulk-rollup:writer-opt-diagnostic:"))) "900" else "300") else if (fields[[3L]] == "reader-check" &&
            "--internal-reader-check" %in% fields[-seq_len(7L)] &&
            grepl("^(gui-cohort-size:(256|768|768-threads8)|gui-impute-workers:(preflight|4|8))$", utils::tail(fields, 1L)))
            c("17179869184", "0", "900") else if (fields[[3L]] == "reader-check" &&
            "--internal-reader-check" %in% fields[-seq_len(7L)] &&
            identical(utils::tail(fields, 1L), "gui-raw-replay:300"))
            c("17179869184", "0", "5400") else if (fields[[3L]] == "browser-baseline" &&
            "--internal-browser-baseline" %in% fields[-seq_len(7L)] &&
            identical(utils::tail(fields, 1L), "capacity:qs2_duckdb:ruv-pearson-300:timed:native:reclaim:extended120"))
            c("17179869184", "0", "7200") else if (fields[[3L]] %in% c("browser-baseline", "verify-da") &&
            paste0("--internal-", fields[[3L]]) %in% fields[-seq_len(7L)] &&
            sub(":first-qc$", "", utils::tail(fields, 1L)) %in% c("capacity:qs2_duckdb:ruv-pearson-500:timed:native",
                "diagnostic:ruv-pearson-500:timed:native",
                "capacity:qs2_duckdb:ruv-pearson-300:timed:native",
                "diagnostic:ruv-pearson-300:timed:native",
                "capacity:qs2_duckdb:ruv-pearson-96:timed:native",
                "diagnostic:ruv-pearson-96:timed:native"))
            c("17179869184", "0", "3600") else if (fields[[3L]] == "verify-da" &&
            "--internal-verify-da" %in% fields[-seq_len(7L)] &&
            grepl("^diagnostic:ruv-pearson(:timed)?(:native)?$", sub(":first-qc$", "", utils::tail(fields, 1L))))
            c("19327352832", "0", "3600") else if (fields[[3L]] == "browser-baseline" &&
            "--internal-browser-baseline" %in% fields[-seq_len(7L)] &&
            grepl("^capacity:qs2_duckdb:ruv-pearson(:timed)?(:native)?$", sub(":first-qc$", "", utils::tail(fields, 1L))))
            c("19327352832", "0", "3600") else if (fields[[3L]] == "browser-baseline" &&
            fields[[7L]] == "3600" && "--internal-browser-baseline" %in% fields[-seq_len(7L)] &&
            identical(sub(":first-qc$", "", utils::tail(fields, 1L)), "capacity:qs2_duckdb:diagnostic"))
            c("17179869184", "0", "3600") else if (fields[[3L]] == "exclusion-allocation")
            c("2147483648", "0", "120") else if (fields[[3L]] %in%
                c("browser-baseline", "verify-da") && fields[[5L]] == "17179869184")
            c("17179869184", "0", "900") else c("8589934592", "0", "900"))) {
        stop("invalid controller request identity or limits")
    }
    list(schema = fields[[1L]], unit = fields[[2L]], operation = fields[[3L]],
        library_path = fields[[4L]], command = fields[-seq_len(7L)],
        memory_max_bytes = as.numeric(fields[[5L]]), swap_max_bytes = 0,
        runtime_max_seconds = as.numeric(fields[[7L]]),
        launch_observed = FALSE, completion_observed = FALSE)
}

#' Parse kernel key/value counters without interpreting overlapping components.
#' @param lines Raw cgroup counter lines.
#' @return Named numeric counters; malformed or duplicate keys are errors.
#' @noRd
parseCheckpointCounters <- function(lines) {
    lines <- trimws(lines)
    lines <- lines[nzchar(lines)]
    if (!length(lines)) return(stats::setNames(numeric(), character()))
    parts <- strsplit(lines, "[[:space:]]+")
    if (any(lengths(parts) != 2L)) stop("invalid cgroup counter record")
    keys <- vapply(parts, `[[`, character(1), 1L)
    raw_values <- vapply(parts, `[[`, character(1), 2L)
    if (anyDuplicated(keys) || any(!grepl("^[0-9]+$", raw_values))) {
        stop("invalid or duplicate cgroup counter")
    }
    values <- as.numeric(raw_values)
    if (any(!is.finite(values))) stop("nonfinite cgroup counter")
    stats::setNames(values, keys)
}

#' Read optional accounting evidence, retaining unavailable/error states explicitly.
#' @param path Kernel counter file.
#' @param raw Preserve lines for device-scoped I/O rather than summing devices.
#' @return Availability and values/error, never fabricated zero counters.
#' @noRd
readCheckpointCounters <- function(path, raw = FALSE) {
    tryCatch({
        if (!file.exists(path)) stop("counter file is absent")
        lines <- readLines(path, warn = FALSE)
        list(status = "available", values = if (raw) lines else as.list(parseCheckpointCounters(lines)))
    }, error = function(error) list(status = "unavailable", error = conditionMessage(error)))
}

#' Verify actual service memory constraints and read its combined peak.
#' @param expected_memory_max Exact authorized service cap, in bytes as text.
#' @return Small cgroup record; no RSS substitution.
#' @noRd
checkpointServiceMemory <- function(expected_memory_max = "8589934592") {
    if (!expected_memory_max %in% c("8589934592", "2147483648", "5368709120",
        "6442450944", "11811160064", "17179869184", "19327352832")) {
        stop("unsupported service memory cap")
    }
    member <- grep("^0::", readLines("/proc/self/cgroup", warn = FALSE), value = TRUE)
    if (length(member) != 1L) stop("cgroup v2 is required")
    path <- file.path("/sys/fs/cgroup", sub("^0::/", "", member))
    read_value <- function(name) readLines(file.path(path, name), n = 1L, warn = FALSE)
    if (!identical(read_value("memory.max"), expected_memory_max) ||
        !identical(read_value("memory.swap.max"), "0")) stop("required service limits are absent")
    list(cgroup = path, peak_bytes = as.numeric(read_value("memory.peak")),
        swap_current_bytes = as.numeric(read_value("memory.swap.current")),
        current_bytes = as.numeric(read_value("memory.current")),
        cpu_cumulative = readCheckpointCounters(file.path(path, "cpu.stat")),
        memory_stat = readCheckpointCounters(file.path(path, "memory.stat")),
        memory_events_cumulative = readCheckpointCounters(file.path(path, "memory.events")),
        io_device_cumulative_raw = readCheckpointCounters(file.path(path, "io.stat"), raw = TRUE),
        accounting_note = paste("Raw snapshots are not phase costs or peak breakdowns.",
            "memory.stat mixes byte components and event counters; retain per-field units.",
            "Memory components may overlap; do not add them to memory.peak.",
            "Device I/O counters are retained separately; no cross-device total is inferred.",
            "Exit samples precede final receipt writing and process teardown."))
}

#' Run an installed production landing smoke, explicitly not a full baseline.
#' @param library_path Installed library.
#' @param run_dir Fresh external request directory.
#' @param inputs NULL for landing only, or explicit report/FASTA paths for import.
#' @return Small receipt written even for ordinary startup failure.
#' @noRd
#' Hold a verified, process-owned sleep inhibitor for one guarded browser run.
#' @param run_dir Existing request-owned output directory.
#' @return Process handle; caller must kill its tree on exit.
#' @noRd
checkpointInhibitSleep <- function(run_dir) {
    who <- paste0("MultiScholaR-proof-", Sys.getpid())
    process <- processx::process$new("systemd-inhibit",
        c("--no-ask-password", "--what=sleep:idle", "--mode=block",
            paste0("--who=", who), "--why=Bounded DIA memory measurement", "sleep", "infinity"),
        stdout = file.path(run_dir, "sleep-inhibitor.log"), stderr = "2>&1",
        cleanup = TRUE, cleanup_tree = TRUE)
    ready <- FALSE
    on.exit(if (!ready) process$kill_tree(), add = TRUE)
    for (attempt in seq_len(20L)) {
        if (!process$is_alive()) stop("measurement sleep inhibitor failed")
        listed <- processx::run("systemd-inhibit", c("--list", "--json=short"), timeout = 2000)
        entries <- jsonlite::fromJSON(listed$stdout, simplifyVector = FALSE)
        ready <- any(vapply(entries, \(entry) identical(entry$who, who) &&
            identical(entry$mode, "block") && identical(entry$what, "sleep:idle"), logical(1)))
        if (ready) break
        Sys.sleep(0.05)
    }
    if (!ready) stop("measurement sleep inhibitor was not registered")
    jsonlite::write_json(list(who = who, pid = process$get_pid(), mode = "block",
        what = "sleep:idle", lifetime = "this browser service only"),
        file.path(run_dir, "sleep-inhibitor.json"), auto_unbox = TRUE, pretty = TRUE)
    process
}

runCheckpointBrowserSmoke <- function(library_path, run_dir, inputs = NULL) {
    previous_progress <- getOption("multischolar.dia_progress_directory")
    options(multischolar.dia_progress_directory = run_dir)
    on.exit(options(multischolar.dia_progress_directory = previous_progress), add = TRUE)
    inhibitor <- checkpointInhibitSleep(run_dir)
    on.exit(inhibitor$kill_tree(), add = TRUE)
    requested_cap <- Sys.getenv("MULTISCHOLAR_DIA_SMALL_GUI_CAP")
    diagnostic_cap <- if (is.null(inputs) && identical(requested_cap, "2147483648")) {
        "2147483648"
    } else if (identical(requested_cap, "19327352832") && identical(inputs$diagnostic_profile, "ruv_pearson")) {
        "19327352832"
    } else if (identical(requested_cap, "17179869184")) {
        "17179869184"
    } else if (identical(requested_cap, "11811160064")) {
        "11811160064"
    } else if (identical(requested_cap, "6442450944")) {
        "6442450944"
    } else if (identical(requested_cap, "5368709120")) {
        "5368709120"
    } else {
        "8589934592"
    }
    limits <- checkpointServiceMemory(diagnostic_cap)
    installation <- inspectCheckpointLibrary(library_path)
    receipt_path <- file.path(run_dir, "smoke-receipt.json")
    if (file.exists(receipt_path)) stop("refusing to overwrite a smoke receipt")
    receipt <- list(kind = if (is.null(inputs)) "installed_production_landing_smoke" else
        if (!is.null(inputs$design)) "installed_production_design_smoke" else
            "installed_production_import_smoke", inputs = inputs, status = "running",
        installation = installation, resources_at_start = limits, promotion_authorized = FALSE)
    receipt$output_policy <- if (identical(inputs$session_format, "native_v1"))
        installation$dia_output_policy else "legacy_s4"
    if (diagnostic_cap != "8589934592") {
        receipt$diagnostic_service_cap_bytes <- as.numeric(diagnostic_cap)
    }
    receipt$diagnostics <- isTRUE(inputs$profile_application)
    if (receipt$diagnostics && !is.null(inputs$capacity)) {
        previous_audit <- options(multischolar.dia_memory_audit_directory = normalizePath(run_dir))
        on.exit(options(previous_audit), add = TRUE)
    }
    if (!is.null(inputs$capacity)) receipt$capacity <- inputs$capacity
    if (receipt$diagnostics) {
        receipt$profile <- if (!is.null(inputs$capacity)) list(
            path = file.path(run_dir, "application-phases.jsonl"), mode = "memory_boundaries",
            scope = "named app boundaries and externally sampled process/cgroup memory",
            allocation_profile = FALSE, cpu_profile = FALSE, forced_gc = FALSE,
            gc_scope = "none; passive boundaries do not collect the application heap",
            time_ceiling_seconds = getOption("multischolar.dia_memory_audit_seconds", 900),
            timing_eligible = FALSE) else list(path = file.path(run_dir, "application.Rprof"),
            scope = "application R process only; driver, browser and parallel workers are not profiled",
            event = "cpu", interval_seconds = 0.02, memory = TRUE, gc = TRUE,
            timing_eligible = FALSE)
    }
    jsonlite::write_json(receipt, receipt_path, auto_unbox = TRUE, pretty = TRUE)
    context <- NULL
    on.exit({
        if (!is.null(context)) tryCatch(context$driver$stop(), error = function(error) {
            receipt$cleanup_error <<- conditionMessage(error)
            receipt$status <<- "failed"
        })
        if (!is.null(context)) tryCatch({
            paths <- file.path(run_dir, paste0("installed-", c("load", "exit"), "-identity.json"))
            sizes <- file.size(paths)
            if (any(is.na(sizes)) || any(sizes <= 0) || any(sizes > 65536)) {
                stop("missing or oversized installed runtime identity")
            }
            identity <- lapply(paths, jsonlite::read_json)
            if (!identical(identity[[1L]]$package_path, installation$package_path) ||
                !identical(identity[[2L]]$package_path, installation$package_path) ||
                !identical(identity[[1L]]$namespace_versions$MultiScholaR,
                    installation$version) ||
                !identical(identity[[1L]]$process_id, identity[[2L]]$process_id)) {
                stop("installed runtime identity disagrees with requested library")
            }
            for (record in identity) {
                native <- record$package_native_libraries
                if (is.null(native)) native <- list()
                if ((length(native) || length(installation$native_libraries)) &&
                    !identical(native, installation$native_libraries))
                    stop("loaded package native bytes disagree with requested library")
            }
            if (!identical(inspectCheckpointLibrary(library_path)$native_libraries,
                installation$native_libraries)) stop("package native bytes changed during execution")
            receipt$installation$load_verified <- TRUE
            receipt$installed_runtime <- list(load = list(path = paths[[1L]], bytes = sizes[[1L]],
                sha256 = digest::digest(file = paths[[1L]], algo = "sha256")),
                exit = list(path = paths[[2L]], bytes = sizes[[2L]],
                    sha256 = digest::digest(file = paths[[2L]], algo = "sha256")),
                note = "package native SHA-256 verified at startup and exit; dependency DLL inventory records paths and sizes")
        }, error = function(error) {
            receipt$installed_identity_error <<- conditionMessage(error)
            receipt$status <<- "failed"
        })
        receipt$resources_at_exit <- checkpointServiceMemory(diagnostic_cap)
        jsonlite::write_json(receipt, receipt_path, auto_unbox = TRUE, pretty = TRUE)
    }, add = TRUE)
    environment <- new.env(parent = baseenv())
    sys.source(file.path(checkpointRunnerDirectory(), "workflow_checkpoint_dia_driver.R"),
        envir = environment)
    tryCatch({
        if (!is.null(inputs)) {
            paths <- inputs[c("report", "fasta", "design")]
            paths <- paths[!vapply(paths, is.null, logical(1))]
            receipt$input_sha256 <- lapply(paths, \(path) digest::digest(file = path, algo = "sha256"))
            if (!is.null(inputs$expected_sha256) &&
                !identical(receipt$input_sha256, inputs$expected_sha256)) {
                stop("stress input bytes differ from the frozen preparation")
            }
        }
        browser_args <- list(library_path = library_path, run_dir = run_dir,
            profile_application = isTRUE(inputs$profile_application))
        if (!is.null(inputs$capacity)) browser_args$capacity_treatment <- inputs$capacity$treatment
        if (!is.null(inputs$session_format)) browser_args$session_format <- inputs$session_format
        if (isTRUE(inputs$reclaim_memory)) browser_args$reclaim_memory <- TRUE
        if (isTRUE(inputs$design_memory)) browser_args$design_memory <- TRUE
        if (isTRUE(inputs$profile_application) && !is.null(inputs$capacity)) {
            browser_args$diagnostic_profile <- inputs$diagnostic_profile
            browser_args$allocator_probe <- isTRUE(inputs$allocator_probe)
            browser_args$diagnostic_source <- file.path(checkpointRunnerDirectory(),
                "workflow_checkpoint_gui_diagnostic.R")
            receipt$diagnostic_source_sha256 <- digest::digest(file = browser_args$diagnostic_source,
                algo = "sha256")
            jsonlite::write_json(receipt, receipt_path, auto_unbox = TRUE, pretty = TRUE)
        }
        context <- do.call(environment$newCheckpointBrowser, browser_args)
        context$driver$wait_for_js(
            "document.querySelectorAll('[data-testid=\"tile-proteomics\"]').length === 1",
            timeout = 60000L)
        if (!is.null(inputs)) environment$runCheckpointImportActions(context$driver, run_dir, inputs)
        if (isTRUE(inputs$design_memory)) environment$requestCheckpointDesignMemory(context$driver, run_dir, "import")
        if (!is.null(inputs$design)) {
            receipt$design_selection_map <- if (isTRUE(inputs$design_memory))
                environment$runCheckpointDesignActions(context$driver, inputs$design, stop_before_save = TRUE) else
                environment$runCheckpointDesignActions(context$driver, inputs$design)
            if (isTRUE(inputs$design_memory)) environment$requestCheckpointDesignMemory(context$driver, run_dir, "design_setup")
            if (!isTRUE(inputs$design_memory)) receipt$executed_design <- environment$checkCheckpointSavedDesign(inputs$design,
                file.path(run_dir, "dia_baseline_smoke", "scripts", "proteomics", "design_matrix.tab"))
            jsonlite::write_json(receipt, receipt_path, auto_unbox = TRUE, pretty = TRUE)
        }
        if (isTRUE(inputs$peptide_diagnostic)) {
            receipt$diagnostic_settings <- if (isTRUE(inputs$first_qc))
                "GUI import/design/Q-value defaults; required figures; reference-only Revert" else if (identical(inputs$diagnostic_profile, "ruv_pearson"))
                "GUI QC defaults; cyclicloess; automatic RUV 1-20%; Pearson Apply 0.5; diagnostic only" else if (isTRUE(inputs$qc_defaults))
                "GUI QC defaults; fixed cyclicloess/RUV Skip/correlation Skip diagnostic route" else
                "existing tiny E2E settings; not large-workload acceptance"
            receipt$peptide_states <- environment$runCheckpointPeptideActions(context$driver,
                diagnostic_settings = !isTRUE(inputs$qc_defaults), first_qc = isTRUE(inputs$first_qc))
            if (isTRUE(inputs$first_qc)) receipt$first_qc <- environment$checkCheckpointFirstQc(context$driver, run_dir)
        }
        if (isTRUE(inputs$protein_diagnostic)) {
            receipt$protein_states <- environment$runCheckpointProteinActions(context$driver,
                diagnostic_settings = !isTRUE(inputs$qc_defaults))
        }
        if (isTRUE(inputs$normalization_diagnostic)) {
            receipt$normalization <- environment$runCheckpointNormalizationActions(context$driver,
                profile = inputs$diagnostic_profile)
        }
        if (isTRUE(inputs$da_diagnostic)) {
            receipt$da <- environment$runCheckpointDaActions(context$driver)
        }
        if (isTRUE(inputs$summary_diagnostic)) {
            receipt$summary_actions <- environment$runCheckpointSummaryActions(context$driver)
        }
        if (identical(inputs$capacity$treatment, "qs2_duckdb")) {
            state_path <- file.path(run_dir, "gui-state-latest.json")
            if (!file.exists(state_path)) stop("candidate state adapter was not observed")
            state <- jsonlite::read_json(state_path)
            if (!identical(state$scope, "full-journey-capacity-v1") ||
                state$states != (if (isTRUE(inputs$design_memory)) 1L else if (isTRUE(inputs$first_qc)) 3L else 16L) || state$native_records != 0 || state$working_records != 0 ||
                state$store$bytes != 0 || state$store$open_resources != 0) {
                stop("candidate did not deliver the complete cold-state journey")
            }
            if (!isTRUE(inputs$first_qc) && !isTRUE(inputs$design_memory) && identical(inputs$diagnostic_profile, "ruv_pearson") &&
                !"ruv_corrected" %in% unlist(state$state_names)) stop("real RUV state was not committed")
            if (identical(inputs$session_format, "native_v1") && state$full_reads != 0L)
                stop("native GUI materialized the complete history")
            if (isTRUE(inputs$first_qc) && (state$component_states != 2L ||
                state$full_s4_reads != 0L || !isTRUE(state$reference_states))) {
                stop("first-QC path materialized an assay or missed component checkpoints")
            }
            if (isTRUE(inputs$first_qc) && any(vapply(state$gui_assay_holders,
                function(holder) isTRUE(holder$resident_assay), logical(1)))) {
                stop("first-QC GUI holder retained an assay")
            }
            receipt$capacity$state_ownership <- state
            receipt$capacity$activation_observed <- TRUE
        }
        jsonlite::write_json(context$driver$get_logs(), file.path(run_dir, "browser-logs.json"),
            auto_unbox = TRUE, pretty = TRUE)
        receipt$status <- if (is.null(inputs)) "passed_landing_only" else "passed_import_only"
        if (!is.null(inputs$design)) receipt$status <- "passed_design_only"
        if (isTRUE(inputs$peptide_diagnostic)) receipt$status <- "passed_peptides_only"
        if (isTRUE(inputs$first_qc)) receipt$status <- "passed_first_qc_ui_outputs_revert"
        if (isTRUE(inputs$protein_diagnostic)) receipt$status <- "passed_qc_only"
        if (isTRUE(inputs$normalization_diagnostic)) {
            receipt$status <- "passed_normalization_ui_only"
        }
        if (isTRUE(inputs$da_diagnostic)) receipt$status <- "passed_da_ui_only"
        if (isTRUE(inputs$summary_diagnostic)) receipt$status <- "passed_summary_ui_only"
        if (isTRUE(inputs$design_memory)) receipt$status <- "completed_import_design_memory_diagnostic"
        receipt
    }, error = function(error) {
        receipt$status <<- "failed"
        receipt$error <<- conditionMessage(error)
        if (!is.null(context)) tryCatch({
            jsonlite::write_json(context$driver$get_logs(), file.path(run_dir, "browser-logs.json"),
                auto_unbox = TRUE, pretty = TRUE)
        }, error = function(log_error) {
            receipt$log_error <<- conditionMessage(log_error)
        })
        stop(error)
    })
}


#' Prepare a matching small diagnostic input set outside any timed member.
#' @param library_path Explicit installed baseline library.
#' @param run_dir Fresh externally owned preparation directory.
#' @return Bound file paths and hashes; no scientific acceptance.
#' @noRd
prepareCheckpointSmallInputs <- function(library_path, run_dir) {
    checkpointServiceMemory()
    inspectCheckpointLibrary(library_path)
    loadNamespace("MultiScholaR", lib.loc = library_path)
    source_root <- normalizePath(file.path(checkpointRunnerDirectory(), "..", ".."))
    helpers <- new.env(parent = baseenv())
    helper_path <- file.path(source_root, "tests", "testthat", "helper-workflow-checkpoint-fixtures.R")
    sys.source(helper_path, envir = helpers)
    recipe <- helpers$readWorkflowCheckpointRecipe("small", source_root)
    object <- helpers$makeWorkflowCheckpointPeptide(recipe)
    data <- methods::slot(object, "peptide_data")
    design <- methods::slot(object, "design_matrix")
    design$sample <- design$Run
    paths <- c(report = file.path(run_dir, "report.tsv"), design = file.path(run_dir, "design.tsv"),
        fasta = file.path(run_dir, "proteins.fasta"))
    if (any(file.exists(paths))) stop("refusing to overwrite diagnostic inputs")
    vroom::vroom_write(data, paths[["report"]], progress = FALSE)
    vroom::vroom_write(design, paths[["design"]], progress = FALSE)
    ids <- unique(data$Protein.Ids)
    sequences <- vapply(ids, function(id) {
        paste(unique(data$Stripped.Sequence[data$Protein.Ids == id]), collapse = "")
    }, character(1))
    fasta <- as.vector(rbind(paste0(">", ids, " synthetic diagnostic sequence"), sequences))
    writeLines(fasta, paths[["fasta"]])
    receipt <- list(kind = "small_fixture_preparation_not_timing", recipe = recipe,
        generator_sha256 = digest::digest(file = helper_path, algo = "sha256"),
        files = as.list(paths), sha256 = as.list(vapply(paths, function(path) {
            digest::digest(file = path, algo = "sha256")
        }, character(1))), synthetic_fasta = TRUE, promotion_authorized = FALSE)
    jsonlite::write_json(receipt, file.path(run_dir, "inputs.json"), auto_unbox = TRUE, pretty = TRUE)
    receipt$files
}

#' Inspect input cache residency without reading or warming file contents.
#' @param path Input report path.
#' @return File size/resident bytes and strict full-residency flag.
#' @noRd
checkpointInputCache <- function(path) {
    command <- Sys.which("fincore")
    if (!nzchar(command)) stop("fincore is required for this memory diagnostic")
    output <- system2(command, c("--json", "--bytes", "--output", "SIZE,RES,FILE", "--", shQuote(path)),
        stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(output, "status"))) stop("input cache inspection failed")
    records <- jsonlite::fromJSON(paste(output, collapse = "\n"), simplifyVector = FALSE)$fincore
    if (length(records) != 1L) stop("ambiguous input cache observation")
    record <- records[[1L]]
    if (!identical(record$file, path) || !is.numeric(record$size) || length(record$size) != 1L ||
        !is.numeric(record$res) || length(record$res) != 1L ||
        !is.finite(record$size) || record$size <= 0 || !is.finite(record$res) || record$res < 0) {
        stop("invalid input cache observation")
    }
    list(size_bytes = as.numeric(record$size), resident_bytes = as.numeric(record$res),
        fully_resident = record$res >= record$size, content_read_or_warmed = FALSE)
}

#' Apply the established file-specific cold-cache policy to frozen replay inputs.
#' @param proof_root Existing proof root; only IVL1gM's validated files are eligible.
#' @param run_dir Fresh request-owned receipt directory, before service startup.
#' @param gui_paths Optional frozen 96-sample report/FASTA/design triple.
#' @return Cache observations; no payload reads, file writes or global cache flush.
#' @noRd
conditionCheckpointMechanismCache <- function(proof_root, run_dir, gui_paths = NULL) {
    if (is.null(gui_paths)) {
    member_dir <- normalizePath(file.path(proof_root, "browser-smoke.IVL1gM"), mustWork = TRUE)
    run_dir <- normalizePath(run_dir, mustWork = TRUE)
    stopifnot(startsWith(run_dir, paste0(normalizePath(proof_root), "/browser-smoke.")))
    member <- jsonlite::read_json(file.path(member_dir, "smoke-receipt.json"))
    stopifnot(identical(member$status, "passed_summary_ui_only"),
        identical(member$input_sha256$report, "2414125be92a11fbe4342fdf354e0214915e30592bec66f3e9d6a984f9290841"))
    store <- member$capacity$state_ownership$store_directory
    project <- file.path(member_dir, "dia_baseline_smoke")
    roots <- c(store, file.path(project, "results", "proteomics", c("protein_qc", "peptide_qc")),
        file.path(project, "results", "proteomics", "publication_graphs", "filtering_qc"))
    paths <- unlist(lapply(roots, \(root) list.files(root, recursive = TRUE,
        all.files = TRUE, no.. = TRUE, full.names = TRUE)), use.names = FALSE)
    source <- file.path(project, "scripts", "proteomics")
    paths <- unique(c(paths, list.files(source, "^[^.].*[.]RDS$", full.names = TRUE)))
    paths <- paths[file.info(paths)$size > 0 & !file.info(paths)$isdir]
    } else {
        member_dir <- normalizePath(proof_root, mustWork = TRUE)
        run_dir <- normalizePath(run_dir, mustWork = TRUE)
        paths <- normalizePath(gui_paths, mustWork = TRUE)
        stopifnot(length(paths) == 3L,
            startsWith(run_dir, paste0(member_dir, "/browser-smoke.")),
            identical(basename(paths), c("report.tsv", "proteins.fasta", "design.tsv")))
        inputs <- jsonlite::read_json(file.path(dirname(paths[[1L]]), "inputs.json"))
        stopifnot(identical(inputs$kind, "gui_ruv_pearson_input_preparation"),
            inputs$selection %in% c("ruv_1gb", "ruv_300", "ruv_500"), identical(inputs$status, "prepared"),
            identical(unlist(inputs$files, use.names = FALSE), paths))
        if (inputs$selection == "ruv_1gb") {
            stopifnot(identical(inputs$sha256$report,
                "3d670aa2567263e9ae90e2966b570e4d65200e9d5adac7214b11eb705a01032a"))
        } else {
            samples <- switch(inputs$selection, ruv_300 = 300, ruv_500 = 500)
            stopifnot(inputs$recipe$n_samples == samples, inputs$recipe$n_proteins == 10000,
                inputs$rows_written == samples * 80000, inputs$max_chunk_rows <= 200000,
                all(vapply(inputs$sha256, \(hash) grepl("^[a-f0-9]{64}$", hash), logical(1))))
        }
    }
    paths <- normalizePath(paths, mustWork = TRUE)
    if (!length(paths) || length(paths) > 10000L || any(!startsWith(paths, paste0(member_dir, "/"))))
        stop("mechanism cache conditioning escaped its frozen source")
    observations <- lapply(paths, function(path) {
        result <- system2("dd", shQuote(c(paste0("if=", path), "iflag=nocache", "count=0", "status=none")),
            stdout = TRUE, stderr = TRUE)
        if (!is.null(attr(result, "status"))) stop("file-specific cache conditioning failed")
        observation <- checkpointInputCache(path)
        if (observation$resident_bytes != 0) stop("mechanism source cache is not cold")
        c(list(path = path), observation)
    })
    receipt <- list(policy = "validated source-file data pages cold; dd iflag=nocache count=0; no code-cache claim",
        status = "conditioned_before_service", files = observations,
        setup_excluded_from_timing = TRUE, payload_read_or_warmed = FALSE)
    jsonlite::write_json(receipt, file.path(run_dir, "source-cache.json"), auto_unbox = TRUE,
        pretty = TRUE, digits = NA)
    invisible(receipt)
}

#' Measure the approved exclusion repair on a fixed report prefix, outside the GUI.
#' @param library_path Explicit installed original or candidate library.
#' @param run_dir Fresh request-owned output directory.
#' @param report_path Frozen representative report.
#' @param memory_diagnostic Observe matched resident input and phase memory counters.
#' @param allocation_diagnostic Profile allocations on the approved 38400-row case.
#' @param qvalue_diagnostic Include the production Q-value preprocessing helper.
#' @return Receipt with classifier-only timing and complete output fingerprint.
#' @noRd
runCheckpointExclusionBenchmark <- function(library_path, run_dir, report_path, memory_diagnostic = FALSE,
    allocation_diagnostic = FALSE, qvalue_diagnostic = FALSE) {
    if (qvalue_diagnostic && !allocation_diagnostic) stop("Q-value probe requires guarded allocation diagnostic")
    if (memory_diagnostic && allocation_diagnostic) stop("diagnostic cases cannot be combined")
    rows_requested <- if (allocation_diagnostic) 38400L else 384000L
    installation <- inspectCheckpointLibrary(library_path)
    memory_snapshot <- function() {
        if (allocation_diagnostic) checkpointServiceMemory("2147483648") else checkpointServiceMemory()
    }
    memory_snapshot()
    path <- file.path(run_dir, "exclusion-receipt.json")
    output <- file.path(run_dir, "classification.rds")
    if (any(file.exists(c(path, output)))) stop("refusing to overwrite exclusion benchmark")
    receipt <- list(kind = "isolated_exclusion_repair_not_GUI_performance", status = "running",
        case = "representative_prefix_384000_v1", rows_requested = 384000L,
        installation = installation, diagnostics = FALSE, promotion_authorized = FALSE,
        boundary = "service startup, namespace load, report integrity/read, classifier, RDS/hash, tree exit",
        expected_report_sha256 = "2414125be92a11fbe4342fdf354e0214915e30592bec66f3e9d6a984f9290841")
    if (memory_diagnostic) {
        receipt$case <- "representative_prefix_384000_memory_diagnostic_v1"
        receipt$diagnostics <- TRUE
        receipt$cache_policy <- "input fully resident at entry and before classification; no warming or cache flushing"
        receipt$memory_scope <- "phase snapshots of anonymous process memory and file cache, not exact R heap or a peak decomposition"
    }
    if (allocation_diagnostic) {
        receipt$case <- "representative_prefix_38400_allocation_diagnostic_v1"
        if (qvalue_diagnostic) receipt$case <- "representative_prefix_38400_qvalue_allocation_v1"
        receipt$rows_requested <- rows_requested
        receipt$diagnostics <- TRUE
        receipt$allocation_threshold_bytes <- 100000
        receipt$memory_scope <- "entry/exit service snapshots; allocations are not retained heap or peak attribution"
        if (!isTRUE(capabilities("profmem"))) stop("allocation profiling unavailable")
    }
    if (memory_diagnostic || allocation_diagnostic) receipt$memory <- list()
    jsonlite::write_json(receipt, path, auto_unbox = TRUE, pretty = TRUE)
    on.exit(jsonlite::write_json(receipt, path, auto_unbox = TRUE, pretty = TRUE, digits = NA), add = TRUE)
    tryCatch({
        if (allocation_diagnostic) receipt$memory$entry <- memory_snapshot()
        if (memory_diagnostic) {
            receipt$cache_entry <- checkpointInputCache(report_path)
            if (!isTRUE(receipt$cache_entry$fully_resident)) stop("input cache condition unmatched; diagnostic held")
            receipt$memory$entry <- memory_snapshot()
        }
        namespace <- loadNamespace("MultiScholaR", lib.loc = library_path)
        if (!identical(normalizePath(getNamespaceInfo(namespace, "path")), installation$package_path)) {
            stop("exclusion benchmark loaded a different installed package")
        }
        receipt$installation$load_verified <- TRUE
        if (memory_diagnostic || allocation_diagnostic) receipt$memory$after_namespace <- memory_snapshot()
        if (!identical(digest::digest(file = report_path, algo = "sha256"), receipt$expected_report_sha256)) {
            stop("report differs from frozen representative input")
        }
        if (memory_diagnostic) receipt$memory$after_report_hash <- memory_snapshot()
        input <- vroom::vroom(report_path, delim = "\t", n_max = rows_requested,
            altrep = FALSE, show_col_types = FALSE, progress = FALSE)
        stopifnot(nrow(input) == rows_requested,
            length(unique(input$Protein.Group)) == if (allocation_diagnostic) 100L else 1000L)
        if (memory_diagnostic) {
            receipt$cache_before_classifier <- checkpointInputCache(report_path)
            if (!isTRUE(receipt$cache_before_classifier$fully_resident)) stop("input cache residency lost; diagnostic held")
            receipt$memory$after_input_read <- memory_snapshot()
        }
        if (allocation_diagnostic) {
            receipt$memory$after_input_read <- memory_snapshot()
            allocation_path <- file.path(run_dir, "classifier.Rprofmem")
            utils::Rprofmem(allocation_path, threshold = receipt$allocation_threshold_bytes)
            on.exit(utils::Rprofmem(NULL), add = TRUE)
        }
        if (qvalue_diagnostic) {
            elapsed <- system.time(result <- getExportedValue("MultiScholaR", "srlQvalueProteotypicPeptideCleanHelper")(
                input_table = input, input_matrix_column_ids = c("Run", "Precursor.Id", "Protein.Group",
                    "Protein.Ids", "Stripped.Sequence", "Modified.Sequence", "Precursor.Charge",
                    "Precursor.Quantity", "Precursor.Normalised"),
                protein_id_column = "Protein.Group", peptide_sequence_column = "Stripped.Sequence",
                modified_peptide_sequence_column = "Modified.Sequence", q_value_column = "Q.Value",
                global_q_value_column = "Global.Q.Value", global_pg_q_value_column = "Global.PG.Q.Value",
                proteotypic_peptide_sequence_column = "Proteotypic", return_exclusion_result = TRUE), gcFirst = FALSE)
        } else {
            elapsed <- system.time(result <- getExportedValue("MultiScholaR", "classifyPeptideBiologicalExclusions")(
                input, contaminant_manifest = NULL, exclude_decoys = TRUE, exclude_contaminants = TRUE), gcFirst = FALSE)
        }
        if (allocation_diagnostic) {
            utils::Rprofmem(NULL)
            receipt$allocation_profile <- list(path = allocation_path,
                bytes = file.info(allocation_path)$size,
                sha256 = digest::digest(file = allocation_path, algo = "sha256"))
        }
        receipt$classifier_timing <- as.list(elapsed[c("user.self", "sys.self", "elapsed")])
        if (qvalue_diagnostic) {
            receipt$qvalue_helper_timing <- receipt$classifier_timing
            receipt$classifier_timing <- NULL
            receipt$retained_object_bytes <- list(input = as.numeric(lobstr::obj_size(input)),
                result = as.numeric(lobstr::obj_size(result)),
                shared_union = as.numeric(lobstr::obj_size(input, result)),
                note = "shared R objects counted once in union; not transient allocations or service peak")
        }
        if (memory_diagnostic || allocation_diagnostic) receipt$memory$after_classifier <- memory_snapshot()
        saveRDS(result, output, compress = FALSE, version = 3)
        receipt$output <- list(path = output, bytes = file.info(output)$size,
            sha256 = digest::digest(file = output, algo = "sha256"))
        if (memory_diagnostic || allocation_diagnostic) receipt$memory$after_output <- memory_snapshot()
        receipt$summary <- if (qvalue_diagnostic) result$exclusion_summary else result$summary
        receipt$status <- "complete_unverified"
        receipt
    }, error = function(error) {
        receipt$status <<- "failed"
        receipt$error <<- conditionMessage(error)
        stop(error)
    })
}

#' Prepare one user-requested half-sample/half-form DIA diagnostic by streaming text.
#' @param report_path Frozen full representative report.
#' @param design_path Matching full design TSV.
#' @param fasta_path Matching FASTA, kept unchanged.
#' @param run_dir Fresh request-owned diagnostic directory.
#' @param spec Expected source hashes and exact derived shape.
#' @return Receipt for a modified diagnostic input; never a frozen control.
#' @noRd
prepareCheckpointHalvedInput <- function(report_path, design_path, fasta_path, run_dir, spec) {
    root <- normalizePath(dirname(run_dir), mustWork = TRUE)
    run_dir <- normalizePath(run_dir, mustWork = TRUE)
    for (path in c(report_path, design_path, fasta_path)) {
        resolved <- normalizePath(path, mustWork = TRUE)
        if (!startsWith(resolved, paste0(root, "/")) || dir.exists(resolved) ||
            is.na(file.size(resolved)) || file.size(resolved) <= 0) {
            stop("half-input source is not an owned nonempty file")
        }
    }
    for (key in c("report_sha256", "design_sha256", "fasta_sha256")) {
        if (!is.character(spec[[key]]) || length(spec[[key]]) != 1L ||
            !grepl("^[0-9a-f]{64}$", spec[[key]])) stop("invalid half-input source fingerprint")
    }
    for (key in c("samples_per_group", "proteins", "peptides_per_protein")) {
        value <- spec[[key]]
        if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
            value < 1 || value != floor(value)) stop("invalid half-input shape")
    }
    source_files <- c(report = report_path, design = design_path, fasta = fasta_path)
    expected_hash <- c(report = spec$report_sha256, design = spec$design_sha256,
        fasta = spec$fasta_sha256)
    for (key in names(source_files)) {
        if (!identical(digest::digest(file = source_files[[key]], algo = "sha256"),
            expected_hash[[key]])) stop(paste("frozen source fingerprint differs:", key))
    }
    design_lines <- readLines(design_path, warn = FALSE)
    fields <- strsplit(design_lines, "\t", fixed = TRUE)
    if (!identical(fields[[1L]][1:2], c("Run", "group")) ||
        any(lengths(fields) != length(fields[[1L]]))) stop("unexpected frozen design columns")
    runs <- vapply(fields[-1L], `[[`, character(1), 1L)
    groups <- vapply(fields[-1L], `[[`, character(1), 2L)
    if (anyDuplicated(runs) || !setequal(unique(groups), c("control", "treated"))) {
        stop("frozen design is not a unique two-group DIA design")
    }
    selected <- c(utils::head(runs[groups == "control"], spec$samples_per_group),
        utils::head(runs[groups == "treated"], spec$samples_per_group))
    if (length(selected) != 2L * spec$samples_per_group || anyDuplicated(selected)) {
        stop("cannot select balanced diagnostic samples")
    }
    selected_lines <- design_lines[c(TRUE, runs %in% selected)]
    if ("tech_rep_group" %in% fields[[1L]]) {
        column <- match("tech_rep_group", fields[[1L]])
        pairs <- table(vapply(strsplit(selected_lines[-1L], "\t", fixed = TRUE),
            `[[`, character(1), column))
        if (any(pairs != 2L)) stop("diagnostic selection split a technical-replicate pair")
    }
    outputs <- c(report = file.path(run_dir, "report.tsv"),
        design = file.path(run_dir, "design.tsv"),
        receipt = file.path(run_dir, "halved-diagnostic-inputs.json"))
    stage <- file.path(run_dir, "report.partial.tsv")
    if (any(file.exists(c(outputs, stage)))) stop("refusing to overwrite diagnostic inputs")
    writeLines(selected_lines, outputs[["design"]])
    input <- file(report_path, open = "r")
    output <- file(stage, open = "w")
    on.exit({try(close(input), silent = TRUE); try(close(output), silent = TRUE)}, add = TRUE)
    header <- readLines(input, n = 1L, warn = FALSE)
    if (length(header) != 1L || !identical(strsplit(header, "\t", fixed = TRUE)[[1L]][[1L]],
        "Run")) stop("unexpected frozen report header")
    writeLines(header, output)
    counts <- integer(length(selected))
    proteins <- character()
    peptide_keys <- character()
    repeat {
        lines <- readLines(input, n = 50000L, warn = FALSE)
        if (!length(lines)) break
        run_ids <- sub("\t.*$", "", lines)
        kept <- lines[run_ids %in% selected &
            grepl("^(?:[^\\t]*\\t){5}2\\t", lines, perl = TRUE)]
        if (!length(kept)) next
        writeLines(kept, output)
        kept_runs <- sub("\t.*$", "", kept)
        counts <- counts + tabulate(match(kept_runs, selected), nbins = length(selected))
        first_run <- kept[kept_runs == selected[[1L]]]
        if (length(first_run)) {
            parts <- strsplit(first_run, "\t", fixed = TRUE)
            ids <- vapply(parts, `[[`, character(1), 2L)
            sequences <- vapply(parts, `[[`, character(1), 4L)
            proteins <- unique(c(proteins, ids))
            peptide_keys <- unique(c(peptide_keys, paste(ids, sequences, sep = "\r")))
        }
    }
    close(input)
    close(output)
    expected_per_run <- spec$proteins * spec$peptides_per_protein
    if (any(counts != expected_per_run) || length(proteins) != spec$proteins ||
        length(peptide_keys) != expected_per_run) {
        stop("halved diagnostic does not preserve the expected protein/peptide shape")
    }
    report_hash <- digest::digest(file = stage, algo = "sha256")
    report_bytes <- file.size(stage)
    if (!file.rename(stage, outputs[["report"]]) ||
        !identical(file.size(outputs[["report"]]), report_bytes)) {
        stop("failed to publish request-owned diagnostic report")
    }
    receipt <- list(kind = "half_samples_and_precursor_forms_diagnostic_v1",
        route = "iq", diagnostic_only = TRUE, promotion_authorized = FALSE,
        source_sha256 = as.list(expected_hash), selected_runs = selected,
        sample_groups = as.list(table(groups[runs %in% selected])),
        proteins = length(proteins), peptides = length(peptide_keys),
        precursor_rows = sum(counts), rows_per_run = as.list(stats::setNames(counts, selected)),
        files = list(report = outputs[["report"]], design = outputs[["design"]],
            fasta = normalizePath(fasta_path, mustWork = TRUE)),
        sha256 = list(report = report_hash,
            design = digest::digest(file = outputs[["design"]], algo = "sha256"),
            fasta = expected_hash[["fasta"]]),
        report_bytes = report_bytes,
        note = "modified input; cannot replace the frozen representative G0/G4 control")
    jsonlite::write_json(receipt, outputs[["receipt"]], auto_unbox = TRUE,
        pretty = TRUE, digits = NA)
    receipt
}

#' Derive compact design/FASTA companions from the frozen generator, not the report.
#' @param recipe Frozen, validated workload recipe.
#' @param helpers Environment containing the existing DIA-002 fixture functions.
#' @return Design table and FASTA lines; no precursor or S4 allocation.
#' @noRd
checkpointInputCompanions <- function(recipe, helpers) {
    design <- helpers$makeWorkflowCheckpointDesign(recipe$n_samples)
    design$sample <- design$Run
    sequences <- helpers$workflowCheckpointSequences(
        seq_len(recipe$n_proteins * recipe$peptides_per_protein))
    proteins <- apply(matrix(sequences, nrow = recipe$peptides_per_protein), 2L,
        paste0, collapse = "")
    headers <- paste0(">", sprintf("P%05d", seq_len(recipe$n_proteins)), " synthetic diagnostic sequence")
    list(design = design, fasta = as.vector(rbind(headers, proteins)))
}

#' Bind an existing frozen report to compact GUI companion inputs, outside timing.
#' @param library_path Explicit installed control library.
#' @param run_dir Fresh externally owned preparation directory.
#' @param recipe_id Frozen DIA-002 recipe identifier.
#' @param report_path Existing report, never rewritten here.
#' @param expected_sha256 Previously frozen report fingerprint.
#' @return Preparation receipt; never authorizes performance acceptance.
#' @noRd
prepareCheckpointCompanionInputs <- function(library_path, run_dir, recipe_id, report_path, expected_sha256) {
    installation <- inspectCheckpointLibrary(library_path)
    source_root <- normalizePath(file.path(checkpointRunnerDirectory(), "..", ".."))
    validateCheckpointOutputRoot(run_dir, source_root)
    if (length(expected_sha256) != 1L || is.na(expected_sha256) ||
        !grepl("^[0-9a-f]{64}$", expected_sha256)) stop("invalid frozen report fingerprint")
    report_path <- normalizePath(report_path, mustWork = TRUE)
    if (dir.exists(report_path) || !identical(digest::digest(file = report_path, algo = "sha256"),
        expected_sha256)) stop("report differs from frozen fingerprint")
    helpers <- new.env(parent = baseenv())
    helper_path <- file.path(source_root, "tests", "testthat", "helper-workflow-checkpoint-fixtures.R")
    sys.source(helper_path, envir = helpers)
    recipe <- helpers$readWorkflowCheckpointRecipe(recipe_id, source_root)
    helpers$validateWorkflowCheckpointRecipe(recipe, allow_large = TRUE,
        output_root = run_dir, repo_root = source_root)
    paths <- c(design = file.path(run_dir, "design.tsv"), fasta = file.path(run_dir, "proteins.fasta"))
    receipt_path <- file.path(run_dir, "inputs.json")
    if (any(file.exists(c(paths, receipt_path)))) stop("refusing to overwrite companion inputs")
    companions <- checkpointInputCompanions(recipe, helpers)
    utils::write.table(companions$design, paths[["design"]], sep = "\t", row.names = FALSE, quote = FALSE)
    writeLines(companions$fasta, paths[["fasta"]])
    receipt <- list(kind = "frozen_report_companions_not_timing", recipe = recipe,
        installation = installation, report = report_path, report_sha256 = expected_sha256,
        files = as.list(paths), sha256 = as.list(vapply(paths, \(path) digest::digest(file = path, algo = "sha256"), character(1))),
        generator_sha256 = digest::digest(file = helper_path, algo = "sha256"),
        synthetic_fasta = TRUE, promotion_authorized = FALSE)
    jsonlite::write_json(receipt, receipt_path, auto_unbox = TRUE, pretty = TRUE)
    receipt
}

#' Validate the small diagnostic's self-contained export and table agreement.
#' @param session Decoded legacy session.
#' @param final Decoded standalone final S4 file.
#' @param table Decoded TSV with explicitly assigned column types.
#' @param ruv_diagnostic Explicit automatic-RUV profile, never inferred from science values.
#' @param native_manager Imported native owner; historical values are checked one at a time.
#' @return Small validation summary; not a whole-workflow science gate.
#' @noRd
checkCheckpointSessionContents <- function(session, final, table, ruv_diagnostic = FALSE,
    native_manager = NULL) {
    expected <- c("initial", "raw_data_s4", "qvalue_filtered", "protein_peptide_filtered",
        "precursor_rollup", "intensity_filtered", "sample_filtered", "replicate_filtered",
        "imputed", "protein_s4_created", "protein_accession_cleaned",
        "protein_intensity_filtered", "duplicates_removed", "protein_replicate_filtered",
        "normalized", "correlation_filtered")
    if (ruv_diagnostic) {
        stopifnot(identical(session$ruv_mode, "automatic"), isTRUE(session$ruv_applied),
            is.numeric(session$ruv_k), length(session$ruv_k) == 1L,
            is.finite(session$ruv_k), session$ruv_k >= 1)
        expected[expected == "normalized"] <- "ruv_corrected"
    }
    history <- unlist(session$r6_state_history, use.names = FALSE)
    state_names <- if (is.null(native_manager)) names(session$r6_complete_states) else
        names(native_manager$.context$snapshot$states)
    get_state <- if (is.null(native_manager)) function(name) session$r6_complete_states[[name]]$data else
        native_manager$getState
    stopifnot(identical(history, expected), !anyDuplicated(state_names),
        all(expected %in% state_names),
        identical(session$r6_current_state_name, "correlation_filtered"))
    for (state in setdiff(expected, "initial")) {
        object <- get_state(state)
        stopifnot(isS4(object), methods::validObject(object))
        object <- NULL
    }
    stopifnot(methods::is(final, "ProteinQuantitativeData"), methods::validObject(final),
        identical(serialize(final, NULL), serialize(session$current_s4_object, NULL)),
        identical(serialize(final, NULL), serialize(session$correlation_filtered_s4, NULL)),
        identical(serialize(final, NULL), serialize(get_state("correlation_filtered"), NULL)))
    quant <- final@protein_quant_table
    key <- final@protein_id_column
    samples <- as.character(final@design_matrix[[final@sample_id]])
    if (!identical(as.numeric(session$final_protein_count),
        as.numeric(length(unique(quant[[key]]))))) stop("exported protein count is inconsistent")
    stopifnot(identical(as.numeric(session$final_sample_count), as.numeric(length(samples))),
        identical(names(table), names(quant)), identical(table[[key]], as.character(quant[[key]])))
    for (sample in samples) {
        stopifnot(identical(is.na(table[[sample]]), is.na(quant[[sample]])),
            identical(is.finite(table[[sample]]), is.finite(quant[[sample]])),
            isTRUE(all.equal(table[[sample]], quant[[sample]], tolerance = 1e-12)))
    }
    manager <- native_manager
    if (is.null(manager)) {
        manager <- MultiScholaR::WorkflowState$new()
        manager$states <- session$r6_complete_states
        manager$state_history <- session$r6_state_history
        manager$current_state <- session$r6_current_state_name
    }
    raw <- get_state("raw_data_s4")
    stopifnot(identical(serialize(manager$getState(), NULL), serialize(final, NULL)),
        identical(manager$revertToState("raw_data_s4"), raw))
    if (!is.null(native_manager)) stopifnot(native_manager$.context$full_reads == 0L,
        native_manager$.context$store$open_resources == 0L)
    list(history_entries = length(history), proteins = nrow(quant), samples = length(samples))
}

#' Compare exported DA science columns without changing their row population.
#' @param expected Reference result table from the installed method.
#' @param actual Exported table.
#' @return TRUE invisibly; identifier/order, missingness and nonfinite signs matter.
#' @noRd
checkCheckpointDaTable <- function(expected, actual) {
    stopifnot(nrow(expected) == nrow(actual), all(names(expected) %in% names(actual)))
    for (column in names(expected)) {
        target <- expected[[column]]
        observed <- actual[[column]]
        if (is.numeric(target)) {
            stopifnot(identical(is.na(target), is.na(observed)),
                identical(is.nan(target), is.nan(observed)),
                identical(is.finite(target), is.finite(observed)),
                isTRUE(all.equal(target, observed, tolerance = 1e-10, check.attributes = FALSE)))
        } else {
            stopifnot(identical(as.character(target), as.character(observed)))
        }
    }
    invisible(TRUE)
}

#' Check exported cluster identities and assignments exactly.
#' @param clusters Named integer cluster vector.
#' @param exported Parsed cluster table.
#' @return TRUE invisibly, otherwise errors.
#' @noRd
checkCheckpointClusterTable <- function(clusters, exported) {
    stopifnot(identical(exported$molecule_id, names(clusters)),
        identical(exported$cluster_id, unname(clusters)))
    invisible(TRUE)
}

#' Check the fixed default heatmap export against independently calculated clusters.
#' @param expected Verified DA table.
#' @param object Final quantitative S4 object.
#' @param settings Observed browser heatmap settings.
#' @param project Read-only browser project.
#' @return Small file/cluster summary, not pixel-level graphical agreement.
#' @noRd
verifyCheckpointHeatmapExport <- function(expected, object, settings, project) {
    stopifnot(settings$heatmap_clustering == "both", settings$heatmap_scaling == "row",
        settings$heatmap_distance_method == "euclidean", settings$heatmap_cluster_method == "ward.D2",
        settings$heatmap_tree_cut_method == "k_clusters")
    selected <- expected[which(expected$comparison == settings$heatmap_contrast &
        expected$fdr_qvalue < 0.05), , drop = FALSE]
    selected <- head(selected[order(-abs(selected$log2FC)), , drop = FALSE], as.integer(settings$heatmap_top_n))
    key <- object@protein_id_column
    quant <- object@protein_quant_table
    samples <- intersect(names(quant), as.character(object@design_matrix[[object@sample_id]]))
    quant <- quant[quant[[key]] %in% selected[[key]], , drop = FALSE]
    values <- as.matrix(quant[samples])
    rownames(values) <- quant[[key]]
    values <- t(scale(t(values)))
    values[!is.finite(values)] <- 0
    clusters <- stats::cutree(stats::hclust(stats::dist(values), method = "ward.D2"),
        k = min(as.integer(settings$heatmap_n_clusters), nrow(values)))
    stem <- gsub("[^A-Za-z0-9_]", "_", paste0("prot_", settings$heatmap_contrast))
    files <- file.path(project, "results", "proteomics", "publication_graphs", "Heatmap",
        paste0(stem, c("_heatmap.png", "_heatmap.pdf", "_heatmap_methods.md", "_clusters.csv")))
    stopifnot(all(file.exists(files)), all(file.info(files)$size > 0))
    image <- png::readPNG(files[[1L]], native = TRUE)
    stopifnot(identical(dim(image), c(2400L, 3000L)))
    rm(image)
    stopifnot(identical(rawToChar(readBin(files[[2L]], "raw", n = 5L)), "%PDF-"))
    metadata <- readLines(files[[3L]], warn = FALSE)
    stopifnot(paste0("| contrast | ", settings$heatmap_contrast, " |") %in% metadata,
        "| clustering_method | ward.D2 |" %in% metadata,
        "| scaling | row |" %in% metadata)
    exported <- utils::read.csv(files[[4L]], stringsAsFactors = FALSE)
    checkCheckpointClusterTable(clusters, exported)
    list(rows = nrow(values), samples = ncol(values), clusters = length(unique(clusters)),
        png_dimensions = c(2400L, 3000L),
        sha256 = setNames(vapply(files, \(file) digest::digest(file = file, algo = "sha256"),
            character(1)), basename(files)))
}

#' Replay the same registered DA method outside browser measurement and check exports.
#' @param session Decoded original pre-DA session.
#' @param project Browser-owned project, read only.
#' @param heatmap_settings Optional observed default heatmap settings.
#' @return Per-contrast rows and file hashes. Not an independent estimator oracle.
#' @noRd
verifyCheckpointDaExports <- function(session, project, heatmap_settings = NULL) {
    contrasts <- session$contrasts_tbl
    stopifnot(nrow(contrasts) > 0L, all(c("contrasts", "full_format") %in% names(contrasts)))
    results <- vector("list", nrow(contrasts))
    heatmap_checked <- is.null(heatmap_settings)
    for (index in seq_len(nrow(contrasts))) {
        object <- session$current_s4_object
        contrast <- data.frame(contrasts = contrasts$full_format[[index]])
        parameters <- list(contrasts_tbl = contrast, formula_string = "~ 0 + group",
            da_q_val_thresh = 0.05, treat_lfc_cutoff = 0, eBayes_trend = TRUE,
            eBayes_robust = TRUE, args_group_pattern = "(\\d+)",
            args_row_id = object@protein_id_column, qvalue_column = "fdr_qvalue",
            raw_pvalue_column = "raw_pvalue")
        for (section in c("differentialAbundanceAnalysis", "differentialAbundanceAnalysisHelper")) {
            if (is.null(object@args[[section]])) object@args[[section]] <- parameters
        }
        reference <- MultiScholaR::differentialAbundanceAnalysis(object,
            contrasts_tbl = contrast, formula_string = "~ 0 + group", da_q_val_thresh = 0.05,
            treat_lfc_cutoff = 0, qvalue_column = "fdr_qvalue", raw_pvalue_column = "raw_pvalue")
        expected <- reference$da_proteins_long
        stopifnot(is.data.frame(expected), nrow(expected) > 0L)
        stem <- paste0("da_proteins_", gsub("[^A-Za-z0-9_-]", "_", contrasts$contrasts[[index]]),
            "_long_annot")
        files <- file.path(project, "results", "proteomics", "da_proteins",
            paste0(stem, c(".tsv", ".xlsx")))
        stopifnot(all(file.exists(files)), all(file.info(files)$size > 0))
        types <- do.call(readr::cols, lapply(expected,
            \(column) if (is.numeric(column)) readr::col_double() else readr::col_character()))
        table <- readr::read_tsv(files[[1L]], col_types = types, progress = FALSE, name_repair = "minimal")
        stopifnot(nrow(readr::problems(table)) == 0L)
        checkCheckpointDaTable(expected, table)
        workbook <- readxl::read_xlsx(files[[2L]], .name_repair = "minimal")
        checkCheckpointDaTable(expected, workbook)
        results[[index]] <- list(contrast = contrasts$full_format[[index]], rows = nrow(expected),
            columns_checked = names(expected),
            sha256 = setNames(vapply(files, \(file) digest::digest(file = file, algo = "sha256"),
                character(1)), basename(files)), qvalue_warnings = reference$qvalue_warnings)
        if (!is.null(heatmap_settings) && heatmap_settings$heatmap_contrast %in% expected$comparison) {
            results[[index]]$heatmap <- verifyCheckpointHeatmapExport(expected, object,
                heatmap_settings, project)
            heatmap_checked <- TRUE
        }
    }
    stopifnot(heatmap_checked)
    results
}

#' Check mandatory IQ/QC files for the explicitly selected diagnostic route.
#' @param member_dir Owned browser member directory.
#' @param ruv_diagnostic Require real-RUV filenames and additional QC figures.
#' @return File-presence evidence only; not full baseline/scientific acceptance.
#' @noRd
checkCheckpointIqQcOutputs <- function(member_dir, ruv_diagnostic = FALSE, native_session = FALSE,
    native_first = FALSE) {
    root <- normalizePath(member_dir, mustWork = TRUE)
    project <- file.path(root, "dia_baseline_smoke")
    source_dir <- file.path(project, "scripts", "proteomics")
    results <- file.path(project, "results", "proteomics")
    protein <- file.path(results, "protein_qc")
    final_stub <- paste0(if (ruv_diagnostic) "ruv_" else "",
        "normalised_results_cln_with_replicates")
    fixed <- c(file.path(root, "qc-plots.pdf"),
        file.path(results, "peptide_qc", "peptide_values_imputed.tsv"),
        file.path(protein, c("iq_output_file.txt", "remove_proteins_with_only_one_rep.tsv",
            "normalized_protein_matrix_pre_ruv.tsv", paste0(final_stub, if (native_first) ".tsv" else c(".tsv", ".RDS")),
            paste0("proteomics_composite_QC_figure.", c("rds", "pdf", "png")),
            paste0(rep(c("pre_norm", "post_norm"), each = 4), "_",
                rep(c("pca", "density", "rle", "correlation"), 2), ".png"))),
        file.path(source_dir, c("aa_seq_tbl_final.RDS", "fasta_metadata.RDS", "design_matrix.tab",
            "ruv_optimization_results.RDS",
            if (native_session) "filtered_session_latest.json" else "filtered_session_data_latest.rds",
            "filtered_session_summary.txt")))
    if (ruv_diagnostic) fixed <- c(fixed, file.path(protein,
        paste0("ruv_corrected_", c("pca", "density", "rle", "correlation", "cancor"), ".png")))
    qc_root <- file.path(results, "publication_graphs", "filtering_qc")
    qc_files <- list.files(qc_root, recursive = TRUE, full.names = TRUE)
    steps <- c("2_qval_filtered", "3_protein_peptide_filtered", "4_precursor_rollup",
        "5_intensity_filtered", "6_sample_filtered", "7_replicate_filtered", "8_imputed")
    plot_names <- c("combined_plots", "proteins_total", "proteins_per_run",
        "peptides_total", "peptides_per_protein", "peptides_per_run")
    expected <- unlist(lapply(steps, \(step) unlist(lapply(plot_names,
        \(name) paste0(step, "_", name, ".", c("png", "pdf"))))), use.names = FALSE)
    matches <- lapply(expected, \(name) qc_files[basename(qc_files) == name])
    missing_qc <- expected[lengths(matches) != 1L]
    files <- c(fixed, unlist(matches[lengths(matches) == 1L], use.names = FALSE))
    info <- file.info(files)
    absent <- files[is.na(info$size) | info$isdir | info$size <= 0]
    if (length(absent) || length(missing_qc)) {
        stop(paste("required IQ/QC outputs missing, empty or ambiguous:",
            paste(c(sub(paste0(root, "/"), "", absent, fixed = TRUE), missing_qc), collapse = ", ")))
    }
    canonical <- normalizePath(files, mustWork = TRUE)
    if (any(!startsWith(canonical, paste0(root, "/")))) stop("required output escaped member storage")
    list(status = "passed_iq_qc_presence_only", files = files,
        bytes = unname(info$size), promotion_authorized = FALSE)
}

#' Verify required selected automatic figures were copied into publication storage.
#' @param project Owned project directory for the completed member.
#' @return Source/copy fingerprints for the selected DIA DA views.
#' @noRd
checkCheckpointPublicationFigures <- function(project) {
    original <- file.path(project, "results", "proteomics", "publication_graphs")
    copies <- file.path(project, "results_summary", "proteomics", "Publication_figures")
    sections <- c("Volcano_Plots", "NumSigDaMolecules", "Heatmap")
    result <- list()
    for (section in sections) {
        source <- file.path(original, section)
        target <- file.path(copies, section)
        if (!dir.exists(source) || !dir.exists(target)) {
            stop(paste("missing publication figure directory:", section))
        }
        relative <- list.files(source, recursive = TRUE, all.files = FALSE, no.. = TRUE)
        relative <- relative[!file.info(file.path(source, relative))$isdir]
        if (!length(relative)) stop(paste("no automatic figures produced in:", section))
        originals <- file.path(source, relative)
        copied <- file.path(target, relative)
        info <- file.info(copied)
        if (any(is.na(info$size) | info$isdir | info$size <= 0)) {
            stop(paste("missing or empty publication figure copy:", section))
        }
        source_hash <- vapply(originals, \(path) digest::digest(file = path, algo = "sha256"), character(1))
        copy_hash <- vapply(copied, \(path) digest::digest(file = path, algo = "sha256"), character(1))
        if (!identical(unname(source_hash), unname(copy_hash))) {
            stop(paste("publication figure copy differs:", section))
        }
        result[[section]] <- list(files = relative, sha256 = stats::setNames(unname(source_hash), relative))
    }
    result
}

#' Check the existing Summary tab's file consumers after its three actions.
#' @param member_dir Completed browser member directory.
#' @param label Label entered in the GUI Summary tab.
#' @param ruv_diagnostic Require the applied-RUV publication table filenames.
#' @return File/content binding for this small diagnostic, not report acceptance.
#' @noRd
checkCheckpointSummaryFiles <- function(member_dir, label, ruv_diagnostic = FALSE, native_first = FALSE) {
    checkpointString(label, "summary experiment label")
    root <- normalizePath(member_dir, mustWork = TRUE)
    project <- file.path(root, "dia_baseline_smoke")
    source <- file.path(project, "scripts", "proteomics")
    summary <- file.path(project, "results_summary", "proteomics")
    integration <- file.path(project, "integration")
    session_state <- list.files(source, "^session_state_[0-9]{4}-[0-9]{2}-[0-9]{2}\\.RDS$", full.names = TRUE)
    if (length(session_state) != 1L) stop("missing or ambiguous summary session-state file")
    files <- c(study_source = file.path(source, "study_parameters.txt"),
        integration_s4 = file.path(integration, sprintf(if (native_first) "proteomics_%s_final.msdia.tar" else "proteomics_%s_final_s4.RDS", label)),
        session_state = session_state,
        study_design = file.path(summary, "Study_report", "design_matrix.tab"),
        study_contrasts = file.path(summary, "Study_report", "contrasts_tbl.tab"),
        study_parameters = file.path(summary, "Study_report", "study_parameters.txt"),
        table_tsv = file.path(summary, "Publication_tables",
            if (ruv_diagnostic) "RUV_normalised_results.tsv" else "normalised_results.tsv"),
        table_rds = file.path(summary, "Publication_tables",
            if (native_first) "filtered_session.msdia.tar" else
                if (ruv_diagnostic) "ruv_normalised_results.RDS" else "normalised_results.RDS"),
        da_table = file.path(summary, "Publication_tables", "DA_results_proteomics.xlsx"),
        qc_source_pdf = file.path(project, "results", "proteomics", "protein_qc",
            "proteomics_composite_QC_figure.pdf"),
        qc_source_png = file.path(project, "results", "proteomics", "protein_qc",
            "proteomics_composite_QC_figure.png"),
        qc_publication_pdf = file.path(summary, "QC_figures", "proteomics_composite_QC_figure.pdf"),
        qc_publication_png = file.path(summary, "QC_figures", "proteomics_composite_QC_figure.png"))
    info <- file.info(files)
    missing <- names(files)[is.na(info$size) | info$isdir | info$size <= 0]
    if (length(missing)) stop(paste("required Summary files absent or empty:", paste(missing, collapse = ",")))
    paths <- normalizePath(files, mustWork = TRUE)
    if (any(!startsWith(paths, paste0(root, "/")))) stop("summary file escaped owned project")
    hashes <- vapply(files, \(path) digest::digest(file = path, algo = "sha256"), character(1))
    if (!identical(hashes[["study_source"]], hashes[["study_parameters"]])) {
        stop("publication study-parameter copy differs from its source")
    }
    for (format in c("pdf", "png")) {
        if (!identical(hashes[[paste0("qc_source_", format)]],
            hashes[[paste0("qc_publication_", format)]])) {
            stop(paste("publication QC composite copy differs from source:", format))
        }
    }
    state <- readRDS(files[["session_state"]])
    if (!is.list(state) || !identical(state$experiment_label, label) ||
        !isTRUE(state$workflow_args_saved) || !isTRUE(state$files_copied)) {
        stop("summary session-state metadata disagrees with completed GUI actions")
    }
    if (native_first) {
        archive <- MultiScholaR:::readDIANativeLatestPath(file.path(source, "filtered_session_latest.json"))
        expected <- digest::digest(file = archive, algo = "sha256")
        if (!identical(hashes[["integration_s4"]], expected) ||
            !identical(hashes[["table_rds"]], expected)) stop("Summary native archive copy differs")
    } else {
        object <- readRDS(files[["integration_s4"]])
        if (!methods::is(object, "ProteinQuantitativeData") || !methods::validObject(object)) {
            stop("summary integration output is not a valid protein S4 object")
        }
    }
    figures <- checkCheckpointPublicationFigures(project)
    list(status = "passed_small_summary_files", files = as.list(files), sha256 = as.list(hashes),
        publication_figures = figures,
        promotion_authorized = FALSE)
}

#' Verify a completed small normalization diagnostic in a separate guarded process.
#' @param library_path Explicit installed library used by the browser member.
#' @param run_dir Fresh verification receipt directory.
#' @param member_dir Browser member directory, never modified.
#' @param verify_da Also replay/check the fixed small DA diagnostic.
#' @param comparison_member Optional small native member; never a large resident control.
#' @param ruv_diagnostic Require the separately bound real RUV/Pearson route.
#' @param native_session Require native archive restoration instead of a legacy RDS.
#' @param require_traces Require diagnostic traces; FALSE only for unprofiled timed members.
#' @return Bounded verification receipt; failures remain on disk and are rethrown.
#' @noRd
verifyCheckpointNormalizationExport <- function(library_path, run_dir, member_dir, verify_da = FALSE,
    comparison_member = NULL, ruv_diagnostic = FALSE, native_session = FALSE,
    require_traces = TRUE) {
    requested_verification_cap <- Sys.getenv("MULTISCHOLAR_DIA_SMALL_VERIFY_CAP")
    verification_cap <- if (ruv_diagnostic && identical(requested_verification_cap, "19327352832")) {
        "19327352832"
    } else if (identical(requested_verification_cap, "17179869184")) {
        "17179869184"
    } else if (identical(requested_verification_cap, "2147483648")) {
        "2147483648"
    } else {
        "8589934592"
    }
    receipt <- list(status = "running", promotion_authorized = FALSE,
        member_dir = normalizePath(member_dir, mustWork = TRUE),
        resources_at_start = checkpointServiceMemory(verification_cap))
    path <- file.path(run_dir, "verification-receipt.json")
    if (file.exists(path)) stop("refusing to overwrite verification receipt")
    on.exit({
        receipt$resources_at_exit <- checkpointServiceMemory(verification_cap)
        jsonlite::write_json(receipt, path, auto_unbox = TRUE, pretty = TRUE)
    }, add = TRUE)
    tryCatch({
        member <- jsonlite::read_json(file.path(member_dir, "smoke-receipt.json"))
        expected_status <- if (verify_da) c("passed_da_ui_only", "passed_summary_ui_only") else
            "passed_normalization_ui_only"
        stopifnot(member$status %in% expected_status,
            identical(normalizePath(library_path), member$installation$library_path))
        if (ruv_diagnostic) stopifnot(identical(member$inputs$diagnostic_profile, "ruv_pearson"),
            identical(member$normalization$ruv, "automatic"), identical(member$normalization$correlation, "apply"))
        ns <- loadNamespace("MultiScholaR", lib.loc = library_path)
        stopifnot(identical(normalizePath(getNamespaceInfo(ns, "path")),
            normalizePath(file.path(library_path, "MultiScholaR"))))
        project <- file.path(member_dir, "dia_baseline_smoke")
        native_first <- native_session && identical(member$output_policy, "native_first_v1")
        receipt$required_iq_qc <- checkCheckpointIqQcOutputs(member_dir, ruv_diagnostic, native_session, native_first)
        source_dir <- file.path(project, "scripts", "proteomics")
        quant_dir <- file.path(project, "results", "proteomics", "protein_qc")
        timestamped <- if (native_session) get("readDIANativeLatestPath", ns)(
            file.path(source_dir, "filtered_session_latest.json")) else
            list.files(source_dir, "^filtered_session_data_[0-9]{8}_[0-9]{6}\\.rds$", full.names = TRUE)
        stopifnot(length(timestamped) == 1L)
        final_stub <- paste0(if (ruv_diagnostic) "ruv_" else "",
            "normalised_results_cln_with_replicates")
        files <- c(latest = file.path(source_dir,
                if (native_session) "filtered_session_latest.json" else "filtered_session_data_latest.rds"),
            timestamped = timestamped, summary = file.path(source_dir, "filtered_session_summary.txt"),
            final = file.path(quant_dir, paste0(final_stub, ".RDS")),
            table = file.path(quant_dir, paste0(final_stub, ".tsv")))
        if (native_first) files <- files[names(files) != "final"]
        stopifnot(all(file.exists(files)), all(file.info(files)$size > 0))
        receipt$files <- lapply(files, \(file) list(path = file, bytes = file.info(file)$size,
            sha256 = digest::digest(file = file, algo = "sha256")))
        native_manager <- NULL
        if (native_session) {
            stopifnot(identical(member$inputs$session_format, "native_v1"),
                !file.exists(file.path(source_dir, "filtered_session_data_latest.rds")))
            restore_root <- file.path(run_dir, "native-restore")
            dir.create(restore_root, mode = "0700")
            imported <- get("openDIANativeSession", ns)(timestamped, restore_root)
            native_manager <- imported$manager
            on.exit(get("closeDIAStateAdapter", ns)(native_manager), add = TRUE)
            session <- get("resolveDIANativeSessionView", ns)(native_manager,
                native_manager$.context$native_metadata_reference, imported$extra_references)
            if (get("isSavedWorkflowReference", ns)(session$current_s4_object)) {
                # Independent verification is an explicit compatibility boundary.
                current_reference <- session$current_s4_object
                correlation_reference <- session$correlation_filtered_s4
                session$current_s4_object <- native_manager$getState()
                session$correlation_filtered_s4 <- if (identical(current_reference, correlation_reference))
                    session$current_s4_object else native_manager$getState("correlation_filtered")
            }
            receipt$native_archive <- list(format = "native_v1", archive_schema = imported$manifest$schema,
                files = length(imported$manifest$files),
                dependencies = length(imported$dependencies), external_integrity = "validated")
        } else {
            stopifnot(identical(receipt$files$latest$sha256, receipt$files$timestamped$sha256))
            session <- readRDS(files[["latest"]])
        }
        final <- if (native_first) session$current_s4_object else readRDS(files[["final"]])
        types <- readr::cols(.default = readr::col_double())
        types$cols[[final@protein_id_column]] <- readr::col_character()
        table <- readr::read_tsv(files[["table"]], col_types = types, name_repair = "minimal", progress = FALSE)
        stopifnot(nrow(readr::problems(table)) == 0L)
        receipt$checks <- checkCheckpointSessionContents(session, final, table,
            ruv_diagnostic = ruv_diagnostic, native_manager = native_manager)
        summary <- readLines(files[["summary"]], warn = FALSE)
        stopifnot(sprintf("- Proteins: %d", receipt$checks$proteins) %in% summary,
            sprintf("- Samples: %d", receipt$checks$samples) %in% summary,
            paste("Session File:", basename(timestamped)) %in% summary)
        receipt$status <- "passed_small_export_checks"
        if (verify_da) {
            stopifnot(identical(member$da$formula, "~ 0 + group"),
                identical(member$da$q_value, 0.05), member$da$treat_lfc == 0)
            receipt$da <- verifyCheckpointDaExports(session, project, member$da$heatmap)
            receipt$status <- "passed_small_da_export_checks"
        }
        if (identical(member$status, "passed_summary_ui_only")) {
            receipt$summary_files <- checkCheckpointSummaryFiles(member_dir,
                member$summary_actions$experiment_label, ruv_diagnostic = ruv_diagnostic, native_first = native_first)
            receipt$status <- "passed_small_summary_file_checks"
        }
        if (!is.null(comparison_member)) {
            driver <- new.env(parent = baseenv())
            sys.source(file.path(checkpointRunnerDirectory(), "workflow_checkpoint_dia_driver.R"), driver)
            receipt$small_pair <- driver$compareCheckpointGuiCapacityMembers(
                comparison_member, member_dir, session)
        }
        if (!is.null(member$capacity)) {
            receipt$scope <- "full-journey-capacity-v1"
            receipt$status <- "passed_capacity_output_checks"
            receipt$scientific_coverage <- paste("Original installed DA/output verifier and full legacy restore;",
                "small pair checks all historical scientific states. No large resident journey was repeated.")
        }
        if (ruv_diagnostic && require_traces) {
            helpers <- new.env(parent = baseenv())
            sys.source(file.path(checkpointRunnerDirectory(), "workflow_checkpoint_gui_diagnostic.R"), helpers)
            evidence <- helpers$readCheckpointRuvPearsonEvidence(member_dir)
            stopifnot(isTRUE(evidence$coverage_complete), all(vapply(evidence$calls,
                function(x) x[["enter"]] >= 1L && x[["exit"]] == x[["enter"]], logical(1))))
            receipt$ruv_pearson <- evidence
            receipt$scope <- "ruv-pearson-memory-v1"
            receipt$status <- "passed_ruv_pearson_output_checks"
            receipt$scientific_coverage <- paste("Real automatic RUV/Pearson execution, complete",
                if (native_session) "native history" else "legacy", "restore, TSV/RDS and DA/required-output agreement; diagnostic only.")
        }
        receipt$session_format <- if (native_session) "native_v1" else "legacy"
        receipt$diagnostic_traces_required <- require_traces
        receipt
    }, error = function(error) {
        receipt$status <<- "failed"
        receipt$error <<- conditionMessage(error)
        stop(error)
    })
}

#' Count genuine memory-profile samples without counting multiline Shiny labels.
#' @param path Unmodified diagnostic application.Rprof file.
#' @param frames Exact function names to count on sample-start lines.
#' @return Conservative inclusive frame counts; never combined-memory attribution.
#' @noRd
summarizeCheckpointProfileFrames <- function(path, frames) {
    lines <- readLines(path, warn = FALSE)
    if (!length(lines) || !grepl("^memory profiling:.*sample.interval=[0-9]+$", lines[[1L]])) {
        stop("unsupported profile header")
    }
    interval <- as.numeric(sub(".*sample.interval=", "", lines[[1L]])) / 1e6
    if (!is.finite(interval) || interval <= 0) stop("invalid profile interval")
    if (!is.character(frames) || anyNA(frames) || anyDuplicated(frames) ||
        any(grepl('["\\n\\r]', frames, perl = TRUE))) stop("invalid profile frame names")
    starts <- grepl("^:[0-9]+:[0-9]+:[0-9]+:[0-9]+:", lines)
    samples <- lines[starts]
    if (!length(samples)) stop("profile contains no complete sample prefixes")
    counts <- vapply(frames, \(frame) sum(grepl(paste0('"', frame, '"'), samples, fixed = TRUE)), integer(1))
    list(sample_records = length(samples), continuation_lines = length(lines) - 1L - length(samples),
        interval_seconds = interval, sampled_seconds = length(samples) * interval,
        frames = data.frame(frame = frames, inclusive_samples_lower_bound = unname(counts),
            inclusive_sampled_seconds_lower_bound = unname(counts) * interval,
            sample_share_lower_bound = unname(counts) / length(samples)),
        limitation = paste("Inclusive counts overlap and must not be added. Only sample-start lines",
            "are searched; frames after multiline labels can be undercounted. No heap peak is inferred."))
}

#' Read retained service timestamps and distinguish main exit from cgroup exit.
#' @param path Captured systemctl key=value properties.
#' @param unit Exact request-owned service unit.
#' @return Main and, when established, tree intervals; no acceptance authorization.
#' @noRd
readCheckpointRetainedInterval <- function(path, unit) {
    lines <- readLines(path, warn = FALSE)
    if (!length(lines) || any(!grepl("=", lines, fixed = TRUE))) stop("invalid service properties")
    keys <- sub("=.*$", "", lines)
    if (anyDuplicated(keys)) stop("ambiguous service properties")
    values <- as.list(stats::setNames(sub("^[^=]*=", "", lines), keys))
    if (!identical(values$Id, unit) || !identical(values$MainPID, "0") ||
        !isTRUE(values$SubState %in% c("exited", "failed"))) stop("service main process is not terminal")
    if (is.null(values$InvocationID) || !grepl("^[0-9a-f]{32}$", values$InvocationID)) {
        stop("missing service invocation identity")
    }
    timestamps <- unlist(values[c("ExecMainStartTimestampMonotonic", "ExecMainExitTimestampMonotonic")])
    if (length(timestamps) != 2L || any(!grepl("^[0-9]+$", timestamps))) stop("missing monotonic timestamps")
    timestamps <- as.numeric(timestamps)
    if (any(timestamps <= 0 | timestamps > 2^53 - 1) || timestamps[[2L]] < timestamps[[1L]]) {
        stop("invalid monotonic interval")
    }
    tree_elapsed <- NULL
    if (identical(values$ExitType, "cgroup") && identical(values$SubState, "exited") &&
        identical(values$ControlGroup, "") && !is.null(values$StateChangeTimestampMonotonic) &&
        grepl("^[0-9]+$", values$StateChangeTimestampMonotonic)) {
        tree_end <- as.numeric(values$StateChangeTimestampMonotonic)
        if (tree_end <= 2^53 - 1 && tree_end >= timestamps[[2L]]) {
            tree_elapsed <- (tree_end - timestamps[[1L]]) / 1e6
        }
    }
    list(status = if (is.null(tree_elapsed)) "main_process_interval_only" else "service_tree_interval_only",
        service_result = if (is.null(values$Result)) "unknown" else values$Result,
        service_succeeded = identical(values$Result, "success") &&
            identical(values$ExecMainCode, "1") && identical(values$ExecMainStatus, "0"),
        unit = unit, invocation_id = values$InvocationID, elapsed_seconds = tree_elapsed,
        main_elapsed_seconds = diff(timestamps) / 1e6, promotion_authorized = FALSE,
        raw_properties = values,
        limitation = "No GUI boundary or outer controller-resource qualification is inferred.")
}

#' Decode exact terminal service metrics, without inferring full-boundary acceptance.
#' @param path Captured journal JSON-lines file.
#' @param unit Exact request-owned service unit.
#' @return Partial metrics with unavailable counters left NULL, never zero-filled.
#' @noRd
readCheckpointTerminalAccounting <- function(path, unit) {
    records <- lapply(readLines(path, warn = FALSE), jsonlite::fromJSON, simplifyVector = FALSE)
    matches <- vapply(records, \(record) identical(record$USER_UNIT, unit) &&
        identical(record$MESSAGE_ID, "ae8f7b866b0347b9af31fe1c80b127c0"), logical(1))
    if (sum(matches) != 1L) stop("missing or ambiguous terminal accounting record")
    record <- records[[which(matches)]]
    if (!is.character(record$USER_INVOCATION_ID) || length(record$USER_INVOCATION_ID) != 1L ||
        !grepl("^[0-9a-f]{32}$", record$USER_INVOCATION_ID)) stop("missing service invocation identity")
    counter <- function(name) {
        value <- record[[name]]
        if (!is.character(value) || length(value) != 1L || !grepl("^[0-9]+$", value)) return(NULL)
        number <- as.numeric(value)
        # UINT64_MAX is systemd's unavailable sentinel; values above exact R
        # integer representation are not safe metric observations either.
        if (!is.finite(number) || number > 2^53 - 1) return(NULL)
        number
    }
    list(status = "terminal_metrics_only", unit = unit, invocation_id = record$USER_INVOCATION_ID,
        cpu_usage_nsec = counter("CPU_USAGE_NSEC"), memory_peak_bytes = counter("MEMORY_PEAK"),
        swap_peak_bytes = counter("MEMORY_SWAP_PEAK"), io_read_bytes = counter("IO_METRIC_READ_BYTES"),
        io_write_bytes = counter("IO_METRIC_WRITE_BYTES"), elapsed_seconds = NULL,
        promotion_authorized = FALSE,
        limitation = "No exact elapsed boundary or controller-resource inclusion is established here.")
}

#' Check recorded-pilot input bindings without loading scientific objects.
#' @param operation Recorded service operation.
#' @param input_path Original session, preparation directory or completed member.
#' @param treatment Explicit resident/saved choice for history members.
#' @return TRUE invisibly; no workload is executed.
#' @noRd
validateRecordedCheckpointInput <- function(operation, input_path, treatment = NULL) {
    if (!operation %in% c("recorded-prepare", "recorded-member", "recorded-timing", "recorded-qc", "recorded-verify")) {
        stop("unsupported recorded operation")
    }
    if (identical(operation, "recorded-qc") && !is.null(treatment)) {
        if (!treatment %in% c("sequence-prepare", "sequence:resident:diagnostic",
            "sequence:qs2:diagnostic", "sequence:resident:timing", "sequence:qs2:timing",
            "sequence:qs2_duckdb:timing"))
            stop("unsupported sequential QC treatment")
        if (treatment != "sequence-prepare") {
            path <- file.path(input_path, "qc-sequence-inputs.json")
            if (!file.exists(path) || file.size(path) > 1048576 ||
                !identical(jsonlite::read_json(path)$status, "prepared_qc_sequence"))
                stop("prepared QC sequence manifest is missing or invalid")
            return(invisible(TRUE))
        }
    }
    if (operation == "recorded-prepare") {
        if (!file.exists(input_path) || !identical(as.numeric(file.size(input_path)),
            239745962)) stop("approved recorded session is absent or has changed size")
    } else if (operation %in% c("recorded-member", "recorded-timing", "recorded-qc")) {
        path <- file.path(input_path, "recorded-inputs.json")
        if (!file.exists(path) || file.size(path) > 1048576) stop("prepared manifest is missing or oversized")
        manifest <- jsonlite::read_json(path)
        if (!identical(manifest$status, "prepared") ||
            !identical(manifest$kind, "recorded_dia_history_pilot")) stop("invalid prepared manifest")
        if (operation %in% c("recorded-member", "recorded-timing") &&
            !isTRUE(treatment %in% c("resident", "saved", "qs2"))) {
            stop("recorded history requires an explicit treatment")
        }
    } else {
        paths <- file.path(input_path, c("recorded-member.json", "recorded-qc.json", "qc-sequence-member.json"))
        present <- paths[file.exists(paths)]
        if (length(present) != 1L ||
            !identical(jsonlite::read_json(present)$status, "complete")) {
            stop("recorded verifier requires one completed member")
        }
    }
    invisible(TRUE)
}

#' Execute one recorded pilot inside the existing guarded service.
#' @param operation,input_path,treatment Explicit validated pilot selection.
#' @param library_path,run_dir Installed package and fresh request directory.
#' @return Diagnostic result with durable success/failure envelope.
#' @noRd
runRecordedCheckpointOperation <- function(operation, library_path, run_dir,
    input_path, treatment = NULL) {
    validateRecordedCheckpointInput(operation, input_path, treatment)
    resources <- checkpointServiceMemory()
    installation <- inspectCheckpointLibrary(library_path)
    artifacts <- file.path(installation$package_path,
        c("DESCRIPTION", "NAMESPACE", "R/MultiScholaR.rdb", "R/MultiScholaR.rdx"))
    installation$sha256 <- as.list(stats::setNames(vapply(artifacts,
        \(path) digest::digest(file = path, algo = "sha256"), character(1)), basename(artifacts)))
    path <- file.path(run_dir, "recorded-operation.json")
    receipt <- list(operation = operation, status = "running", installation = installation,
        input_path = normalizePath(input_path), treatment = treatment,
        resources_at_start = resources, promotion_authorized = FALSE)
    on.exit(jsonlite::write_json(receipt, path, auto_unbox = TRUE, pretty = TRUE,
        null = "null", digits = NA), add = TRUE)
    environment <- new.env(parent = baseenv())
    driver <- file.path(checkpointRunnerDirectory(), "workflow_checkpoint_dia_driver.R")
    receipt$driver_sha256 <- digest::digest(file = driver, algo = "sha256")
    sys.source(driver, envir = environment)
    tryCatch({
        result <- switch(operation,
            "recorded-prepare" = environment$prepareRecordedCheckpointSession(
                library_path, run_dir, input_path),
            "recorded-member" = environment$executeRecordedCheckpointReplay(
                library_path, run_dir, input_path, treatment, checkpointServiceMemory),
            "recorded-timing" = environment$executeRecordedCheckpointReplay(
                library_path, run_dir, input_path, treatment, checkpointServiceMemory,
                profile_owners = FALSE),
            "recorded-qc" = if (is.null(treatment)) environment$executeRecordedActiveQc(
                library_path, run_dir, input_path, checkpointServiceMemory) else
                if (treatment == "sequence-prepare") environment$prepareRecordedQcSequence(
                    library_path, run_dir, input_path) else {
                    selection <- strsplit(treatment, ":", fixed = TRUE)[[1L]]
                    environment$executeRecordedQcSequence(library_path, run_dir, input_path,
                        selection[[2L]], selection[[3L]] == "diagnostic", checkpointServiceMemory)
                },
            "recorded-verify" = if (file.exists(file.path(input_path, "qc-sequence-member.json")))
                environment$verifyRecordedQcSequence(library_path, run_dir, input_path) else
                if (file.exists(file.path(input_path, "recorded-member.json")))
                environment$verifyRecordedCheckpointReplay(library_path, run_dir, input_path) else
                environment$verifyRecordedActiveQc(library_path, run_dir, input_path))
        receipt$status <- "complete"
        receipt$installation$load_verified <- identical(
            normalizePath(getNamespaceInfo(asNamespace("MultiScholaR"), "path")),
            installation$package_path)
        receipt$resources_at_exit <- checkpointServiceMemory()
        result
    }, error = function(error) {
        receipt$status <<- "failed"
        receipt$error <<- conditionMessage(error)
        stop(error)
    }, interrupt = function(error) {
        receipt$status <<- "interrupted"
        receipt$error <<- conditionMessage(error)
        stop(error)
    })
}

#' Run the fixed reader screening block using only the existing guarded services.
#' @param library_path,output_root,prepared_dir Explicit installed build, disk root and inputs.
#' @param adapter Run the frozen seven-block installed-adapter comparison.
#' @return Preserved service bindings and descriptive screening decision.
#' @noRd
runCheckpointReaderScreen <- function(library_path, output_root, prepared_dir, adapter = FALSE) {
    inspectCheckpointLibrary(library_path)
    output_root <- checkpointExistingDirectory(output_root, "output-root")
    prepared_dir <- checkpointExistingDirectory(prepared_dir, "reader preparation")
    run_dir <- tempfile("reader-screen.", tmpdir = output_root)
    if (!dir.create(run_dir)) stop("cannot create reader screening directory")
    receipt <- list(status = "running", scope = if (adapter)
        "DIA-013 installed-adapter seven-block decision" else "DIA-013 reader screening only",
        library = library_path, prepared_dir = prepared_dir, members = list(),
        promotion_authorized = FALSE)
    path <- file.path(run_dir, "screen.json")
    flush <- function() jsonlite::write_json(receipt, path,
        auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
    on.exit(flush(), add = TRUE)
    tryCatch({
        for (case in c("small", "wide", "precursor")) {
            for (block in if (adapter) seq_len(7L) else c(1L, 5L)) {
                repetitions <- if (adapter) 1L else block
                engines <- if (adapter) c("adapter_rds", "adapter_duckdb", "adapter_fastduckdb") else
                    c("rds", "arrow", "duckdb", "fastduckdb")
                if (if (adapter) block %% 2L == 0L else repetitions == 5L) engines <- rev(engines)
                for (engine in engines) {
                    key <- paste(case, engine, block, sep = ":")
                    selection <- paste(case, engine, repetitions, sep = ":")
                    member <- runCheckpointService("reader-member", library_path,
                        output_root, c(prepared_dir, selection))
                    receipt$members[[key]] <- list(operation = member, status = "unverified")
                    flush()
                    if (member$exit_status != 0L) stop(paste("reader member failed:", key))
                    verification <- runCheckpointService("reader-verify", library_path,
                        output_root, member$member_dir)
                    receipt$members[[key]]$verification <- verification
                    flush()
                    if (verification$exit_status != 0L) stop(paste("reader verification failed:", key))
                    terminal <- jsonlite::read_json(file.path(member$member_dir, "terminal-interval.json"))
                    if (!isTRUE(terminal$service_succeeded) ||
                        !identical(terminal$status, "service_tree_interval_only")) stop("reader completion is not qualified")
                    counters <- terminal$raw_properties
                    numeric_counter <- function(name) {
                        value <- suppressWarnings(as.numeric(counters[[name]]))
                        if (length(value) != 1L || !is.finite(value) || value < 0) stop("reader metric missing")
                        value
                    }
                    details <- jsonlite::read_json(file.path(member$member_dir, "reader-member.json"))
                    if (!identical(details$case, case) || !identical(details$engine, engine) ||
                        !identical(details$repetitions, as.integer(repetitions)) ||
                        !identical(details$status, "complete") || details$open_resources != 0 ||
                        numeric_counter("MemorySwapPeak") != 0) stop("reader member binding/resources differ")
                    verified <- jsonlite::read_json(file.path(verification$member_dir, "reader-verification.json"))
                    if (!identical(verified$status, "passed_reader_exactness") ||
                        !identical(verified$member, normalizePath(member$member_dir)) ||
                        !identical(verified$input_sha256, details$input_sha256) ||
                        !identical(verified$case, case) || !identical(verified$engine, engine)) {
                        stop("reader exactness verification binding differs")
                    }
                    receipt$members[[key]]$status <- "verified"
                    receipt$members[[key]]$input_sha256 <- details$input_sha256
                    receipt$members[[key]]$run_id <- terminal$invocation_id
                    receipt$members[[key]]$outputs <- details$outputs
                    receipt$members[[key]]$metrics <- list(elapsed_seconds = terminal$elapsed_seconds,
                        cpu_seconds = numeric_counter("CPUUsageNSec") / 1e9,
                        peak_bytes = numeric_counter("MemoryPeak"), swap_bytes = numeric_counter("MemorySwapPeak"),
                        producer_plus_consumer_seconds = terminal$elapsed_seconds + details$producer$elapsed_seconds,
                        producer_plus_consumer_cpu_seconds = numeric_counter("CPUUsageNSec") / 1e9 + details$producer$cpu_seconds,
                        input_file_bytes = details$input_file_bytes,
                        action_seconds = sum(vapply(details$events, \(event) event$elapsed_seconds, numeric(1))),
                        selected_consumer_seconds = sum(vapply(details$events,
                            \(event) sum(vapply(event$value[setdiff(names(event$value), "resources")],
                                \(read) read$elapsed_seconds, numeric(1))), numeric(1))))
                    flush()
                    message("Reader screening verified: ", key)
                }
            }
        }
        receipt$decisions <- list()
        if (adapter) {
            metrics <- c("elapsed_seconds", "cpu_seconds", "peak_bytes")
            for (case in c("small", "wide", "precursor")) {
                comparisons <- list()
                for (control_engine in c("adapter_rds", "adapter_duckdb")) {
                    controls <- lapply(seq_len(7L), \(block)
                        receipt$members[[paste(case, control_engine, block, sep = ":")]])
                    candidates <- lapply(seq_len(7L), \(block)
                        receipt$members[[paste(case, "adapter_fastduckdb", block, sep = ":")]])
                    for (block in seq_len(7L)) {
                        if (!identical(controls[[block]]$input_sha256, candidates[[block]]$input_sha256)) {
                            stop("paired adapter inputs differ")
                        }
                    }
                    ratios <- stats::setNames(lapply(metrics, function(metric) {
                        checkpointPairedRatio(vapply(candidates, \(x) x$metrics[[metric]], numeric(1)),
                            vapply(controls, \(x) x$metrics[[metric]], numeric(1)), 0.95)
                    }), metrics)
                    point_ok <- ratios$elapsed_seconds$point <= 1 && ratios$cpu_seconds$point <= 1 &&
                        ratios$peak_bytes$point < 1
                    upper_ok <- ratios$elapsed_seconds$upper <= 1 && ratios$cpu_seconds$upper <= 1 &&
                        ratios$peak_bytes$upper < 1
                    comparisons[[control_engine]] <- list(ratios = ratios,
                        outcome = if (!point_ok) "not_adopted" else if (!upper_ok) "held" else
                            "eligible_component_evidence")
                }
                outcomes <- vapply(comparisons, \(x) x$outcome, character(1))
                receipt$decisions[[case]] <- list(comparisons = comparisons,
                    outcome = if (any(outcomes == "not_adopted")) "not_adopted" else
                        if (any(outcomes == "held")) "held" else "eligible_component_evidence",
                    gui_adoption = FALSE, resident_replacement = FALSE)
            }
            receipt$status <- "complete"
            message("Installed adapter decision receipt: ", path)
            return(invisible(receipt))
        }
        for (case in c("small", "wide", "precursor")) {
            for (repetitions in c(1L, 5L)) {
                control <- receipt$members[[paste(case, "rds", repetitions, sep = ":")]]$metrics
                candidate <- receipt$members[[paste(case, "fastduckdb", repetitions, sep = ":")]]$metrics
                metrics <- c("elapsed_seconds", "cpu_seconds", "peak_bytes",
                    "producer_plus_consumer_seconds", "producer_plus_consumer_cpu_seconds")
                ratios <- stats::setNames(vapply(metrics,
                    \(metric) candidate[[metric]] / control[[metric]], numeric(1)), metrics)
                eligible <- all(is.finite(ratios)) && all(ratios <= 1) && ratios[["peak_bytes"]] < 1
                receipt$decisions[[paste(case, repetitions, sep = ":")]] <- list(
                    status = if (eligible) "eligible_for_seven_blocks" else "not_adopted_screen",
                    ratios_vs_rds = as.list(ratios), uncertainty_qualified = FALSE)
            }
        }
        receipt$status <- "complete"
        message("Reader screening receipt: ", path)
        invisible(receipt)
    }, error = function(error) {
        receipt$status <<- "failed"
        receipt$error <<- conditionMessage(error)
        stop(error)
    })
}

#' Execute the approved capacity sequence, serially, stopping at the first failure.
#' @param options Existing explicit GUI CLI options.
#' @param contract Separate capacity contract; never supplied to the paired evaluator.
#' @param service Existing guarded service invoker, injectable for tests.
#' @return Capacity receipts without a performance ratio or adoption decision.
#' @noRd
runCheckpointGuiCapacity <- function(options, contract, service = runCheckpointService) {
    source_root <- normalizePath(file.path(checkpointRunnerDirectory(), "..", ".."))
    driver <- new.env(parent = baseenv())
    sys.source(file.path(checkpointRunnerDirectory(), "workflow_checkpoint_dia_driver.R"), driver)
    expected <- driver$readCheckpointGuiCapacitySpec(source_root)
    if (!identical(contract, expected) || !identical(options$route, "iq") ||
        isTRUE(options$diagnostics) || !identical(digest::digest(file = options$protocol,
            algo = "sha256"), expected$protocol_sha256)) stop("capacity scope/boundary differs")
    installation <- inspectCheckpointLibrary(options$library)
    root <- validateCheckpointOutputRoot(options[["output-root"]], source_root)
    if (isTRUE(options[["validate-only"]])) {
        return(invisible(list(status = "eligible_capacity_only_no_workload",
            scope = expected$scope, selections = expected$run_order, promotion_authorized = FALSE)))
    }
    build <- digest::digest(file = file.path(installation$package_path, "R", "MultiScholaR.rdb"), algo = "sha256")
    directory <- file.path(root, paste0("gui-capacity-", substr(build, 1L, 16L)))
    if (dir.exists(directory) || !dir.create(directory)) {
        stop("this build's capacity attempt already exists; inspect it instead of repeating")
    }
    path <- file.path(directory, "capacity.json")
    receipt <- list(scope = expected$scope, status = "running", installation = installation,
        build_sha256 = build, contract_sha256 = digest::digest(file = options$contract, algo = "sha256"),
        protocol_sha256 = expected$protocol_sha256, preparations = list(), members = list(),
        ratios = NULL, promotion_authorized = FALSE)
    flush <- function() jsonlite::write_json(receipt, path, auto_unbox = TRUE, pretty = TRUE,
        null = "null", digits = NA)
    flush()
    on.exit(flush(), add = TRUE)
    invoke <- function(operation, args) service(operation, options$library, root, args)
    tryCatch({
        for (selection in unlist(expected$run_order)) {
            parts <- strsplit(selection, ":", fixed = TRUE)[[1L]]
            recipe <- parts[[1L]]
            treatment <- parts[[2L]]
            if (is.null(receipt$preparations[[recipe]])) {
                prepared <- invoke("prepare-companions", paste0("capacity:", recipe))
                receipt$preparations[[recipe]] <- prepared
                flush()
                if (prepared$exit_status != 0L) stop(paste("capacity preparation failed:", recipe))
            }
            input_path <- file.path(receipt$preparations[[recipe]]$member_dir, "inputs.json")
            input <- jsonlite::read_json(input_path)
            if (!identical(input$status, "prepared") || !identical(input$selection, recipe) ||
                !identical(input$contract_sha256, receipt$contract_sha256)) stop("capacity input binding differs")
            member <- invoke("browser-baseline", c(unlist(input$files[c("report", "fasta", "design")]),
                paste0("--capacity-", treatment)))
            receipt$members[[selection]] <- list(member = member, input_receipt = input_path)
            flush()
            if (member$exit_status != 0L) stop(paste("candidate/control member stopped:", selection))
            reference <- if (identical(selection, "small:qs2_duckdb"))
                receipt$members[["small:resident"]]$member$member_dir else character()
            verified <- invoke("verify-da", c(member$member_dir, reference, "--capacity"))
            receipt$members[[selection]]$verification <- verified
            receipt$members[[selection]]$terminal <- jsonlite::read_json(
                file.path(member$member_dir, "terminal-interval.json"))
            flush()
            if (verified$exit_status != 0L) stop(paste("capacity verification failed:", selection))
        }
        receipt$status <- "complete_capacity_only"
        receipt
    }, error = function(error) {
        receipt$status <<- "stopped"
        receipt$error <<- conditionMessage(error)
        stop(error)
    })
}

#' CLI entry point; unsupported or unqualified routes never fall back silently.
#' @param arguments Command-line options.
#' @return Exits or errors; no workload is run by the current partial implementation.
#' @noRd
checkpointProofMain <- function(arguments) {
    if (length(arguments) == 6L &&
        identical(arguments[[1L]], "--internal-perf107-browser")) {
        environment <- new.env(parent = baseenv())
        sys.source(file.path(checkpointRunnerDirectory(),
            "workflow_checkpoint_dia_driver.R"), envir = environment)
        return(environment$runPerf107ManagedBrowser(arguments[[2L]], arguments[[3L]],
            arguments[[4L]], arguments[[5L]], arguments[[6L]]))
    }
    if (length(arguments) == 5L && identical(arguments[[1L]], "--workflow-lifecycle")) {
        if (!arguments[[5L]] %in% c("protocol", "resources")) stop("invalid lifecycle suite")
        library <- normalizePath(arguments[[2L]], mustWork = TRUE)
        root <- normalizePath(arguments[[3L]], mustWork = TRUE)
        test_root <- normalizePath(arguments[[4L]], mustWork = TRUE)
        installation <- inspectCheckpointLibrary(library)
        withr::local_libpaths(c(library, .libPaths()))
        namespace <- asNamespace("MultiScholaR")
        if (!identical(normalizePath(getNamespaceInfo(namespace, "path")),
            normalizePath(file.path(library, "MultiScholaR")))) stop("installed runtime differs")
        receipt <- list(schema = 1L, status = "starting", suite = arguments[[5L]],
            installation = installation, promotion_authorized = FALSE)
        receipt_path <- file.path(root, "lifecycle-launch.json")
        if (file.exists(receipt_path)) stop("lifecycle receipt already exists")
        on.exit(jsonlite::write_json(receipt, receipt_path, auto_unbox = TRUE, pretty = TRUE), add = TRUE)
        value <- tryCatch({
            config <- get("launchWorkflowApplication", namespace)(root, diagnostic = TRUE,
                mode = arguments[[5L]], test_root = test_root)
            receipt$config <- config
            results <- readRDS(file.path(root, "lifecycle-tests.rds"))
            if (any(results$failed > 0L | results$error)) stop("lifecycle suite failed")
            receipt$status <- "complete"
            receipt$tests <- results[c("test", "passed", "failed", "error", "skipped")]
            get("workflowSystemd", namespace)(c("stop", config$slice), required = FALSE)
            config
        }, error = function(error) {
            receipt$status <<- "failed"
            receipt$error <<- conditionMessage(error)
            stop(error)
        })
        return(invisible(value))
    }

    if (length(arguments) == 6L && identical(arguments[[1L]], "--condition-gui-cache")) {
        return(conditionCheckpointMechanismCache(arguments[[2L]], arguments[[3L]], arguments[4:6]))
    }
    if (length(arguments) == 3L && identical(arguments[[1L]], "--condition-mechanism-cache")) {
        return(conditionCheckpointMechanismCache(arguments[[2L]], arguments[[3L]]))
    }
    if (length(arguments) == 4L && arguments[[1L]] == "--internal-prepare-capacity") {
        checkpointServiceMemory()
        inspectCheckpointLibrary(arguments[[2L]])
        driver <- new.env(parent = baseenv())
        sys.source(file.path(checkpointRunnerDirectory(), "workflow_checkpoint_dia_driver.R"), driver)
        return(driver$prepareCheckpointGuiCapacityInputs(arguments[[2L]], arguments[[3L]],
            arguments[[4L]], normalizePath(file.path(checkpointRunnerDirectory(), "..", ".."))))
    }
    if (length(arguments) == 7L && arguments[[1L]] == "--internal-browser-baseline" &&
        grepl("^capacity:(resident|qs2_duckdb)(:(diagnostic|ruv-pearson|ruv-pearson-small|ruv-pearson-96|ruv-pearson-300|ruv-pearson-500)(:timed)?)?(:native)?(:first-qc)?(:reclaim)?(:design-memory)?(:extended120)?$", arguments[[7L]])) {
        extended120 <- grepl(":extended120$", arguments[[7L]])
        arguments[[7L]] <- sub(":extended120$", "", arguments[[7L]])
        design_memory <- grepl(":design-memory$", arguments[[7L]])
        request_marker <- sub(":design-memory$", "", arguments[[7L]])
        reclaim_memory <- grepl(":reclaim$", request_marker)
        request_marker <- sub(":reclaim$", "", request_marker)
        first_qc <- grepl(":first-qc$", request_marker)
        route_marker <- sub(":first-qc$", "", request_marker)
        native_session <- grepl(":native$", route_marker)
        marker <- sub(":native$", "", route_marker)
        timed <- grepl(":timed$", marker)
        marker <- sub(":timed$", "", marker)
        ruv_small <- identical(marker, "capacity:qs2_duckdb:ruv-pearson-small")
        ruv_stress <- marker %in% paste0("capacity:qs2_duckdb:ruv-pearson-", c(96, 300, 500))
        stress_samples <- if (ruv_stress) sub("^.*-", "", marker) else NULL
        ruv_pearson <- ruv_small || ruv_stress || identical(marker, "capacity:qs2_duckdb:ruv-pearson")
        diagnostic <- design_memory || !first_qc && !timed && (ruv_pearson || identical(marker, "capacity:qs2_duckdb:diagnostic"))
        if (diagnostic || ruv_pearson) {
            ceiling <- Sys.getenv("MULTISCHOLAR_DIA_DIAGNOSTIC_SECONDS", "900")
            if (!ceiling %in% c("900", "3600", if (extended120) "7200") ||
                extended120 && ceiling != "7200") stop("invalid diagnostic time ceiling")
            options(multischolar.dia_extended120 = extended120)
            options(multischolar.dia_memory_audit_seconds = as.numeric(ceiling))
        }
        treatment <- sub(":(diagnostic|ruv-pearson|ruv-pearson-small|ruv-pearson-96|ruv-pearson-300|ruv-pearson-500)$", "", sub("^capacity:", "", marker))
        if (first_qc && (!native_session || !identical(treatment, "qs2_duckdb") ||
            identical(stress_samples, "500"))) stop("invalid first-QC route")
        if (design_memory && (!identical(stress_samples, "300") || !native_session ||
            reclaim_memory || first_qc)) stop("invalid import/design diagnostic route")
        if (native_session && !identical(treatment, "qs2_duckdb")) stop("native session requires the saved-state route")
        if (ruv_pearson) {
            checkpointServiceMemory(if (ruv_small) "8589934592" else if (ruv_stress)
                "17179869184" else "19327352832")
            spec_environment <- new.env(parent = baseenv())
            sys.source(file.path(checkpointRunnerDirectory(), "workflow_checkpoint_dia_driver.R"), spec_environment)
            spec <- spec_environment$readCheckpointRuvPearsonSpec(normalizePath(file.path(checkpointRunnerDirectory(), "..", "..")))
            input <- jsonlite::read_json(file.path(dirname(arguments[[4L]]), "inputs.json"))
            if (!identical(input$kind, "gui_ruv_pearson_input_preparation") ||
                !identical(input$scope, spec$scope) || !identical(input$status, "prepared") ||
                !digest::digest(file = file.path(arguments[[2L]], "MultiScholaR", "R", "MultiScholaR.rdb"),
                    algo = "sha256") %in% if (native_session) c(spec$native_session_rdb_sha256,
                        if (first_qc) unlist(spec$reference_first_qc_rdb_sha256)) else
                        c(spec$runtime_rdb_sha256, spec$export_repair_rdb_sha256))
                stop("RUV/Pearson runtime or input binding differs")
            native <- inspectCheckpointLibrary(arguments[[2L]])$native_libraries
            rdb <- digest::digest(file = file.path(arguments[[2L]], "MultiScholaR", "R", "MultiScholaR.rdb"), algo = "sha256")
            if (extended120 && (!identical(stress_samples, "300") || !native_session ||
                !reclaim_memory || first_qc || design_memory ||
                !rdb %in% unlist(spec$extended_120_rdb_sha256))) stop("120-minute route is not approved for this request")
            if (length(native) && !identical(native, spec$native_library_bindings[[rdb]]))
                stop("native library is not bound to this approved runtime")
            if (ruv_small && (!identical(input$selection, "ruv_small") ||
                input$recipe$precursor_rows != 192000 || input$recipe$n_samples != 24 ||
                file.size(arguments[[4L]]) > 50000000)) stop("small RUV pilot input differs")
            if (ruv_stress) {
                stress <- spec[[paste0("stress_", stress_samples)]]
                allowed_builds <- c(stress$runtime_rdb_sha256,
                    if (stress_samples %in% c("96", "300")) spec$integrated_candidate_rdb_sha256,
                    if (!first_qc && stress_samples %in% c("96", "300"))
                        unlist(spec[[paste0("native_first_", stress_samples, "_rdb_sha256")]]),
                    if (first_qc && stress_samples %in% c("96", "300")) unlist(spec$reference_first_qc_rdb_sha256))
                if (!native_session || !timed || !identical(input$selection, stress$recipe$id) ||
                    !identical(input$recipe, stress$recipe) ||
                    input$rows_written != stress$recipe$precursor_rows || input$max_chunk_rows > 200000 ||
                    !digest::digest(file = file.path(arguments[[2L]], "MultiScholaR",
                        "R", "MultiScholaR.rdb"), algo = "sha256") %in% allowed_builds) {
                    stop("sample stress build, recipe or route differs")
                }
            } else if (input$selection %in% c("ruv_300", "ruv_500")) stop("sample stress requires the 16 GiB route")
        }
        if (identical(treatment, "resident")) {
            input_path <- file.path(dirname(arguments[[4L]]), "inputs.json")
            if (!file.exists(input_path)) stop("capacity resident execution is small-fixture-only")
            input <- jsonlite::read_json(input_path)
            if (!identical(input$kind, "gui_capacity_input_preparation") ||
                !identical(input$selection, "small") || input$recipe$precursor_rows != 24576 ||
                file.size(arguments[[4L]]) > 20000000) stop("large resident capacity runs are forbidden")
        }
        return(runCheckpointBrowserSmoke(arguments[[2L]], arguments[[3L]],
            inputs = list(report = arguments[[4L]], fasta = arguments[[5L]], design = arguments[[6L]],
                design_memory = design_memory,
                peptide_diagnostic = !design_memory, protein_diagnostic = !first_qc && !design_memory,
                normalization_diagnostic = !first_qc && !design_memory,
                da_diagnostic = !first_qc && !design_memory, summary_diagnostic = !first_qc && !design_memory, qc_defaults = TRUE,
                first_qc = first_qc,
                profile_application = diagnostic,
                allocator_probe = ruv_small && diagnostic && native_session &&
                    isTRUE(spec$anonymous_memory_probe) && identical(digest::digest(file =
                        file.path(arguments[[2L]], "MultiScholaR", "R", "MultiScholaR.rdb"),
                        algo = "sha256"), "a91ed1e91c743d8d4734ead86573827e5ef619f84574c5e76246f07dc9e49ef4"),
                diagnostic_profile = if (ruv_pearson) "ruv_pearson" else NULL,
                expected_sha256 = if (ruv_stress) input$sha256 else NULL,
                reclaim_memory = reclaim_memory,
                session_format = if (native_session) "native_v1" else "legacy",
                capacity = list(scope = "full-journey-capacity-v1", treatment = treatment))))
    }
    if (length(arguments) == 5L && identical(arguments[[1L]], "--internal-verify-da") &&
        grepl("^diagnostic:ruv-pearson(-small|-96|-300|-500)?(:timed)?(:native)?$", arguments[[5L]])) {
        return(verifyCheckpointNormalizationExport(arguments[[2L]], arguments[[3L]],
            arguments[[4L]], verify_da = TRUE, ruv_diagnostic = TRUE,
            native_session = grepl(":native$", arguments[[5L]]),
            require_traces = !grepl(":timed", arguments[[5L]], fixed = TRUE)))
    }
    if (length(arguments) == 5L && identical(arguments[c(1L, 5L)],
        c("--internal-verify-da", "native-session"))) {
        return(verifyCheckpointNormalizationExport(arguments[[2L]], arguments[[3L]],
            arguments[[4L]], verify_da = TRUE, native_session = TRUE))
    }
    if (length(arguments) == 5L && arguments[[1L]] == "--internal-verify-da") {
        return(verifyCheckpointNormalizationExport(arguments[[2L]], arguments[[3L]],
            arguments[[4L]], verify_da = TRUE, comparison_member = arguments[[5L]]))
    }
    if (length(arguments) == 5L && identical(arguments[c(1L, 5L)], c("--reader-screen", "adapter"))) {
        return(runCheckpointReaderScreen(arguments[[2L]], arguments[[3L]], arguments[[4L]], adapter = TRUE))
    }
    if (length(arguments) == 4L && arguments[[1L]] == "--reader-screen") {
        return(runCheckpointReaderScreen(arguments[[2L]], arguments[[3L]], arguments[[4L]]))
    }
    if (length(arguments) %in% c(3L, 4L) && arguments[[1L]] == "--validate-reader") {
        environment <- new.env(parent = baseenv())
        sys.source(file.path(checkpointRunnerDirectory(), "workflow_checkpoint_dia_driver.R"),
            envir = environment)
        return(environment$validateCheckpointReaderInput(arguments[[2L]], arguments[[3L]],
            if (length(arguments) == 4L) arguments[[4L]] else NULL))
    }
    if (length(arguments) %in% c(4L, 5L) && arguments[[1L]] %in% paste0(
        "--internal-reader-", c("build", "check", "prepare", "inspect", "member", "verify"))) {
        if (identical(arguments[[1L]], "--internal-reader-check") &&
            startsWith(utils::tail(arguments, 1L), "gui-bulk-rollup:")) {
            inhibitor <- checkpointInhibitSleep(arguments[[3L]])
            on.exit(inhibitor$kill_tree(), add = TRUE)
        }
        checkpointServiceMemory(if (identical(arguments[[1L]], "--internal-reader-check") &&
            (grepl("^(perf108:(B[12]|T[12]|F0|N[12]|I[1-3]|J[12]|K[1-4]|M[12]|Q[CT]|U[BT]|Z[ACOVFSQL]|E[1-8]|R[CS][01]|P[AVO][01]|[YW](S[1-4]|[RPW][1-6])|PC|X(S[1-4]|R[1-6]|F[1-8])|A[SF]1|H[SF]1|C[23]|P[SF][12]|R[SFR][12]|V([1-9]|1[0-6]))|perf101:(D[01]|D0R|V)|perf106:(L[01]|S[1-6]|F[1-8]|R[1-4]|V2-(A|V|I|J|D[12]|Q([1-9]|1[0-9]|2[0-4]))))$", utils::tail(arguments, 1L)) || grepl("^gui-bulk-rollup:(qualification-preflight|qualification-prepare|writer-build|writer-opt-(build|preflight|diagnostic:(primary|contrast))|writer-r-preflight|writer-r-readback:[A-Za-z0-9]{6}|writer-r:(screen|full):(primary|contrast):P1[456]:[12]|io-preflight|io-prepare|io:(duckdb|fastduckdb):[12]|scan-inspection|crossover:(B|P0):(4608|9600|96000|960000):[12]|contrast:(B|P([0-9]|1[0-4]))|diagnostic|diagnostic-preflight|preflight|relational-preflight|materialize-preflight|adaptation-preflight|fast-p3-preflight|writer-preflight|screen:([ABCD]|P([0-9]|1[0-4])):[12]|full:([ABCD]|P([0-9]|1[0-4])))$", utils::tail(arguments, 1L)) || identical(utils::tail(arguments, 1L), "gui-raw-replay:300") ||
                grepl("^(gui-cohort-size:(256|768|768-threads8)|gui-impute-workers:(preflight|4|8))$", utils::tail(arguments, 1L)) ||
                grepl("^perf107:(B|C|D|E|H|HF|S|L1|L2)$", utils::tail(arguments, 1L)) ||
                grepl("^core:((00[6-9]|01[0-4]|017):(control|candidate):[12]|lifecycle:(complete|fail)|011:overflow)$", utils::tail(arguments, 1L))))
            "17179869184" else "8589934592")
        environment <- new.env(parent = baseenv())
        sys.source(file.path(checkpointRunnerDirectory(), "workflow_checkpoint_dia_driver.R"),
            envir = environment)
        return(environment$runCheckpointReaderOperation(sub("^--internal-", "", arguments[[1L]]),
            arguments[[2L]], arguments[[3L]], arguments[[4L]],
            if (length(arguments) == 5L) arguments[[5L]] else NULL))
    }
    if (length(arguments) %in% c(3L, 4L) && arguments[[1L]] == "--validate-recorded") {
        return(validateRecordedCheckpointInput(arguments[[2L]], arguments[[3L]],
            if (length(arguments) == 4L) arguments[[4L]] else NULL))
    }
    if (length(arguments) %in% c(4L, 5L) && arguments[[1L]] %in% paste0(
        "--internal-", c("recorded-prepare", "recorded-member", "recorded-timing", "recorded-qc", "recorded-verify"))) {
        return(runRecordedCheckpointOperation(sub("^--internal-", "", arguments[[1L]]),
            arguments[[2L]], arguments[[3L]], arguments[[4L]],
            if (length(arguments) == 5L) arguments[[5L]] else NULL))
    }
    if (length(arguments) == 5L && arguments[[1L]] == "--internal-utility-prepare") {
        checkpointServiceMemory()
        environment <- new.env(parent = baseenv())
        sys.source(file.path(checkpointRunnerDirectory(), "workflow_checkpoint_dia_driver.R"),
            envir = environment)
        return(environment$prepareCheckpointUtilityInputs(arguments[[2L]], arguments[[3L]],
            arguments[[4L]], arguments[[5L]]))
    }
    utility_diagnostics <- length(arguments) > 0L &&
        identical(arguments[[1L]], "--internal-utility-member") &&
        identical(utils::tail(arguments, 1L), "--diagnostic-memory")
    if (utility_diagnostics) arguments <- utils::head(arguments, -1L)
    if (length(arguments) %in% c(8L, 9L) && arguments[[1L]] == "--internal-utility-member") {
        checkpointServiceMemory()
        environment <- new.env(parent = baseenv())
        sys.source(file.path(checkpointRunnerDirectory(), "workflow_checkpoint_dia_driver.R"),
            envir = environment)
        return(tryCatch(environment$executeCheckpointUtilityMember(arguments[[2L]],
            arguments[[3L]], arguments[[4L]], arguments[[5L]], arguments[[6L]],
            arguments[[7L]], as.integer(arguments[[8L]]),
            if (length(arguments) == 9L) arguments[[9L]] else NULL,
            diagnostics = utility_diagnostics,
            memory_snapshot = checkpointServiceMemory), error = function(error) {
            receipt_path <- file.path(arguments[[3L]], "member-receipt.json")
            receipt <- if (file.exists(receipt_path)) jsonlite::read_json(receipt_path) else list()
            receipt$status <- "failed"
            receipt$error <- conditionMessage(error)
            receipt$diagnostics <- utility_diagnostics
            receipt$promotion_authorized <- FALSE
            jsonlite::write_json(receipt, receipt_path, auto_unbox = TRUE, pretty = TRUE,
                null = "null")
            stop(error)
        }))
    }
    if (length(arguments) == 4L && arguments[[1L]] == "--internal-utility-verify") {
        checkpointServiceMemory()
        environment <- new.env(parent = baseenv())
        sys.source(file.path(checkpointRunnerDirectory(), "workflow_checkpoint_dia_driver.R"),
            envir = environment)
        return(environment$verifyCheckpointUtilityMember(arguments[[2L]],
            arguments[[3L]], arguments[[4L]]))
    }
    if (length(arguments) == 5L && arguments[[1L]] == "--internal-prepare-half") {
        return(prepareCheckpointHalvedInput(arguments[[2L]], arguments[[3L]],
            arguments[[4L]], arguments[[5L]], list(
                report_sha256 = "2414125be92a11fbe4342fdf354e0214915e30592bec66f3e9d6a984f9290841",
                design_sha256 = "4fbd08901905c1339af58ab619d75fade53655484153d627bd7cd39a45bfe6e5",
                fasta_sha256 = "352b56704c330f3e59646101df2c52a469e5ebe17da652fdd0d3ea85e672f6d4",
                samples_per_group = 12L, proteins = 10000L, peptides_per_protein = 4L)))
    }
    if (length(arguments) == 5L && identical(arguments[c(1L, 5L)], c("--internal-exclusion-allocation", "qvalue"))) {
        return(runCheckpointExclusionBenchmark(arguments[[2L]], arguments[[3L]], arguments[[4L]],
            allocation_diagnostic = TRUE, qvalue_diagnostic = TRUE))
    }
    if (length(arguments) == 4L && arguments[[1L]] %in% c("--internal-exclusion-benchmark", "--internal-exclusion-memory", "--internal-exclusion-allocation")) {
        return(runCheckpointExclusionBenchmark(arguments[[2L]], arguments[[3L]], arguments[[4L]],
            memory_diagnostic = identical(arguments[[1L]], "--internal-exclusion-memory"),
            allocation_diagnostic = identical(arguments[[1L]], "--internal-exclusion-allocation")))
    }
    if (length(arguments) == 6L && arguments[[1L]] == "--internal-prepare-companions") {
        return(prepareCheckpointCompanionInputs(arguments[[2L]], arguments[[3L]],
            arguments[[4L]], arguments[[5L]], arguments[[6L]]))
    }
    if (length(arguments) == 3L && arguments[[1L]] == "--internal-retained-interval") {
        if (file.exists(file.path(arguments[[2L]], "terminal-interval.json"))) stop("refusing to overwrite interval receipt")
        result <- readCheckpointRetainedInterval(file.path(arguments[[2L]], "terminal-properties.txt"),
            arguments[[3L]])
        jsonlite::write_json(result, file.path(arguments[[2L]], "terminal-interval.json"),
            auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
        return(invisible(result))
    }
    if (length(arguments) == 3L && arguments[[1L]] == "--internal-accounting-journal") {
        if (file.exists(file.path(arguments[[2L]], "terminal-metrics.json"))) stop("refusing to overwrite metrics receipt")
        result <- readCheckpointTerminalAccounting(file.path(arguments[[2L]], "terminal-journal.jsonl"),
            arguments[[3L]])
        jsonlite::write_json(result, file.path(arguments[[2L]], "terminal-metrics.json"),
            auto_unbox = TRUE, pretty = TRUE, null = "null")
        return(invisible(result))
    }
    if (length(arguments) == 1L && arguments[[1L]] %in%
        c("--internal-accounting-probe", "--internal-accounting-probe-failure", "--internal-accounting-probe-child")) {
        allocation <- raw(64 * 1024^2)
        allocation[] <- as.raw(1)
        cat("bounded accounting probe: ", length(allocation), " bytes touched\n", sep = "")
        if (arguments[[1L]] == "--internal-accounting-probe-child") {
            system2("/usr/bin/sleep", "1", wait = FALSE)
        }
        quit(status = if (arguments[[1L]] == "--internal-accounting-probe-failure") 7L else 0L)
    }
    if (length(arguments) == 4L && arguments[[1L]] == "--internal-verify-da") {
        return(verifyCheckpointNormalizationExport(arguments[[2L]], arguments[[3L]], arguments[[4L]],
            verify_da = TRUE))
    }
    if (length(arguments) == 6L && arguments[[1L]] %in%
        c("--internal-browser-da", "--internal-browser-da-default-qc", "--internal-browser-da-profile")) {
        return(runCheckpointBrowserSmoke(arguments[[2L]], arguments[[3L]],
            inputs = list(report = arguments[[4L]], fasta = arguments[[5L]],
                design = arguments[[6L]], peptide_diagnostic = TRUE,
                protein_diagnostic = TRUE, normalization_diagnostic = TRUE, da_diagnostic = TRUE,
                qc_defaults = arguments[[1L]] %in% c("--internal-browser-da-default-qc", "--internal-browser-da-profile"),
                profile_application = identical(arguments[[1L]], "--internal-browser-da-profile"))))
    }
    if (length(arguments) == 6L && arguments[[1L]] == "--internal-browser-summary") {
        return(runCheckpointBrowserSmoke(arguments[[2L]], arguments[[3L]],
            inputs = list(report = arguments[[4L]], fasta = arguments[[5L]], design = arguments[[6L]],
                peptide_diagnostic = TRUE, protein_diagnostic = TRUE,
                normalization_diagnostic = TRUE, da_diagnostic = TRUE,
                summary_diagnostic = TRUE)))
    }
    if (length(arguments) == 6L && arguments[[1L]] == "--internal-browser-baseline") {
        return(runCheckpointBrowserSmoke(arguments[[2L]], arguments[[3L]],
            inputs = list(report = arguments[[4L]], fasta = arguments[[5L]], design = arguments[[6L]],
                peptide_diagnostic = TRUE, protein_diagnostic = TRUE,
                normalization_diagnostic = TRUE, da_diagnostic = TRUE,
                summary_diagnostic = TRUE, qc_defaults = TRUE, profile_application = FALSE)))
    }
    if (length(arguments) == 4L && arguments[[1L]] == "--internal-verify-normalization") {
        return(verifyCheckpointNormalizationExport(arguments[[2L]], arguments[[3L]], arguments[[4L]]))
    }
    if (length(arguments) == 6L && arguments[[1L]] == "--internal-browser-normalization") {
        return(runCheckpointBrowserSmoke(arguments[[2L]], arguments[[3L]],
            inputs = list(report = arguments[[4L]], fasta = arguments[[5L]],
                design = arguments[[6L]], peptide_diagnostic = TRUE,
                protein_diagnostic = TRUE, normalization_diagnostic = TRUE)))
    }
    if (length(arguments) == 6L && arguments[[1L]] == "--internal-browser-proteins") {
        return(runCheckpointBrowserSmoke(arguments[[2L]], arguments[[3L]],
            inputs = list(report = arguments[[4L]], fasta = arguments[[5L]],
                design = arguments[[6L]], peptide_diagnostic = TRUE, protein_diagnostic = TRUE)))
    }
    if (length(arguments) == 3L && arguments[[1L]] == "--internal-prepare-small") {
        return(prepareCheckpointSmallInputs(arguments[[2L]], arguments[[3L]]))
    }
    if (length(arguments) == 6L && arguments[[1L]] == "--internal-browser-peptides") {
        return(runCheckpointBrowserSmoke(arguments[[2L]], arguments[[3L]],
            inputs = list(report = arguments[[4L]], fasta = arguments[[5L]],
                design = arguments[[6L]], peptide_diagnostic = TRUE)))
    }
    if (length(arguments) == 6L && arguments[[1L]] == "--internal-browser-design") {
        return(runCheckpointBrowserSmoke(arguments[[2L]], arguments[[3L]],
            inputs = list(report = arguments[[4L]], fasta = arguments[[5L]], design = arguments[[6L]])))
    }
    if (length(arguments) == 5L && arguments[[1L]] == "--internal-browser-import") {
        return(runCheckpointBrowserSmoke(arguments[[2L]], arguments[[3L]],
            inputs = list(report = arguments[[4L]], fasta = arguments[[5L]])))
    }
    if (length(arguments) == 3L && arguments[[1L]] == "--check-service-paths") {
        inspectCheckpointLibrary(arguments[[2L]])
        validateCheckpointOutputRoot(arguments[[3L]],
            file.path(checkpointRunnerDirectory(), "..", ".."))
        return(invisible(TRUE))
    }
    if (length(arguments) == 3L && arguments[[1L]] == "--internal-browser-smoke") {
        return(runCheckpointBrowserSmoke(arguments[[2L]], arguments[[3L]]))
    }
    options <- parseCheckpointProofArgs(arguments)
    if (isTRUE(options$help)) return(checkpointProofHelp())
    protocol <- jsonlite::read_json(options$protocol, simplifyVector = TRUE)
    validateCheckpointProofProtocol(protocol)
    contract <- jsonlite::read_json(options$contract, simplifyVector = FALSE)
    if (identical(options$mode, "gui") && identical(contract$schema, "multischolar.dia_gui_capacity.v1")) {
        return(runCheckpointGuiCapacity(options, contract))
    }
    validateCheckpointScenarioContract(contract, options$route)
    if (!identical(digest::digest(file = options$contract, algo = "sha256"),
        "cfdf7da2fc5a6123da269f85011550d8ffd2e8a611de071e0cd1a52a6f35d152")) {
        stop("scenario contract differs from the frozen DIA-002 identity")
    }
    if (!identical(digest::digest(file = options$protocol, algo = "sha256"),
        "517edccb0906165698f3c1f165dfe7d961a17910d93363a60e843fa466c62379")) {
        stop("protocol differs from the frozen DIA-002 identity")
    }
    if (identical(options$route, "limpa")) stop("Limpa execution is held until DIA-017; no IQ fallback")
    inspectCheckpointLibrary(options$library)
    source_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
    if (length(source_file) != 1L) stop("cannot resolve runner source location", call. = FALSE)
    source_root <- file.path(dirname(normalizePath(source_file)), "..", "..")
    validateCheckpointOutputRoot(options[["output-root"]], source_root)
    if (identical(options$mode, "consumers")) {
        if (isTRUE(options[["validate-only"]])) {
            return(invisible(list(status = "eligible_no_workload", route = "iq",
                mode = "consumers", scenarios = c("resident", "complete_rds",
                    "product_lookup", "asset_lookup"), promotion_authorized = FALSE)))
        }
        if (isTRUE(options$diagnostics)) {
            stop("diagnostic profiling is not a consumer correctness member",
                call. = FALSE)
        }
        consumer_dir <- tempfile("consumers.", tmpdir = options[["output-root"]])
        if (!dir.create(consumer_dir)) stop("cannot create owned consumer evidence directory")
        qc_dir <- file.path(consumer_dir, "qc")
        da_dir <- file.path(consumer_dir, "da")
        projection_dir <- file.path(consumer_dir, "projection")
        export_dir <- file.path(consumer_dir, "export")
        if (!dir.create(qc_dir) || !dir.create(da_dir) ||
            !dir.create(projection_dir) || !dir.create(export_dir)) {
            stop("cannot create owned consumer scenario directories")
        }
        environment <- new.env(parent = baseenv())
        sys.source(file.path(checkpointRunnerDirectory(), "workflow_checkpoint_dia_driver.R"),
            envir = environment)
        qc <- environment$runCheckpointQcConsumerScenarios(options$library,
            qc_dir, options$contract)
        da <- environment$runCheckpointDaConsumerScenarios(options$library,
            da_dir, options$contract)
        projection <- environment$runCheckpointProjectionScenarios(options$library,
            projection_dir)
        export <- environment$runCheckpointExportScenarios(options$library,
            export_dir, options$contract)
        if (!identical(qc$status, "passed_qc_consumer_scenarios") ||
            !identical(da$status, "passed_da_consumer_scenarios") ||
            !identical(projection$status, "passed_projection_scenarios") ||
            !identical(export$status, "passed_export_scenarios")) {
            stop("consumer scenarios did not complete")
        }
        result <- list(status = "passed_consumer_scenarios", route = "iq",
            qc = qc, da = da, projection = projection, export = export,
            performance_claim = FALSE,
            promotion_authorized = FALSE)
        jsonlite::write_json(result, file.path(consumer_dir, "completion.json"),
            auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
        return(invisible(result))
    }
    if (identical(options$mode, "utility")) {
        inventory <- validateCheckpointScenarioContract(contract, "iq")
        recipes <- contract$recipes[match(inventory$timing_recipe_ids,
            vapply(contract$recipes, \(recipe) recipe$id, character(1)))]
        names(recipes) <- inventory$timing_recipe_ids
        candidate_treatment <- if (identical(options[["checkpoint-codec"]], "qs2")) "qs2" else "saved"
        retained_recipes <- if (identical(candidate_treatment, "qs2")) "small" else character()
        measured_recipes <- recipes[setdiff(names(recipes), retained_recipes)]
        history_scope <- if (identical(candidate_treatment, "qs2"))
            "standalone-qs2-large-history-v3" else "standalone-rds-history-v1"
        if (isTRUE(options[["validate-only"]])) {
            return(invisible(list(status = "eligible_no_workload", route = "iq",
                mode = "utility", scope_version = history_scope,
                recipes = names(measured_recipes), resident_recipe_ids = retained_recipes,
                promotion_authorized = FALSE)))
        }
        if (isTRUE(options$diagnostics)) {
            stop("diagnostic profiling is not a timed utility member", call. = FALSE)
        }
        utility_dir <- tempfile("utility.", tmpdir = options[["output-root"]])
        if (!dir.create(utility_dir)) stop("cannot create owned utility evidence directory")
        protocol_hash <- digest::digest(file = options$protocol, algo = "sha256")
        prepared <- lapply(measured_recipes, function(recipe) {
            operation <- runCheckpointService("utility-prepare", options$library,
                options[["output-root"]], c(options$contract, recipe$id))
            if (operation$exit_status != 0L ||
                !file.exists(file.path(operation$member_dir, "utility-inputs.json"))) {
                stop("utility input preparation failed; retained service evidence")
            }
            operation$member_dir
        })
        record <- function(result) {
            treatment <- if (is.null(result$treatment)) "unknown" else result$treatment
            path <- file.path(utility_dir, sprintf("%s_%s_%02d_%s.json",
                result$recipe, result$variant, result$block, treatment))
            if (file.exists(path)) stop("refusing to replace utility member receipt")
            jsonlite::write_json(result, path, auto_unbox = TRUE, pretty = TRUE,
                null = "null", digits = NA)
        }
        run <- function(recipe, block, treatment, variant) {
            runCheckpointUtilityServiceMember(recipe, block, treatment, variant,
                options$library, options[["output-root"]], options$contract,
                input_dir = prepared[[recipe$id]])
        }
        records <- dispatchCheckpointUtility(recipes, protocol, run, record, candidate_treatment)
        decisions <- list()
        for (recipe in names(recipes)) for (variant in
            c("independent", "shared_metadata", "retained_alias")) {
            selected <- Filter(\(item) identical(item$recipe, recipe) &&
                identical(item$variant, variant), records)
            key <- paste(recipe, variant, sep = ":")
            if (recipe %in% retained_recipes) {
                decisions[[key]] <- list(outcome = "resident_policy_retained",
                    reason = "small history remains resident; historical qs2 failure is preserved, not accepted",
                    promotion_authorized = FALSE)
                next
            }
            if (identical(variant, "retained_alias")) {
                decisions[[key]] <- list(outcome = "negative_control_recorded",
                    reason = "candidate deliberately retained full aliases; excluded from adoption ratio",
                    promotion_authorized = FALSE)
                next
            }
            frame <- checkpointUtilityFrame(selected, protocol_hash)
            control <- selected[[which(vapply(selected,
                \(item) identical(item$member, "control"), logical(1)))[[1L]]]]
            candidate <- selected[[which(vapply(selected,
                \(item) identical(item$member, "candidate"), logical(1)))[[1L]]]]
            expected <- list(route = "iq", input_hash = control$input_hash,
                protocol_hash = protocol_hash, boundary_hash = control$boundary_hash,
                dependency_hash = control$dependency_hash,
                control_build = control$build_hash, candidate_build = candidate$build_hash)
            decisions[[key]] <- evaluateCheckpointPairs(frame, protocol, expected,
                gate = "g2_history", comparison = "original_control")
        }
        adoption <- Filter(\(decision) !decision$outcome %in%
            c("negative_control_recorded", "resident_policy_retained"), decisions)
        outcomes <- vapply(adoption,
            \(decision) decision$outcome, character(1))
        status <- if (all(outcomes == "accepted_evidence"))
            "all_required_utility_evidence_passed" else if (any(outcomes == "failed"))
                "required_utility_evidence_failed" else "required_utility_evidence_held"
        completion <- list(status = status, route = "iq", mode = "utility",
            scope_version = history_scope, resident_recipe_ids = retained_recipes,
            measured_recipe_ids = names(measured_recipes),
            members = length(records), decisions = decisions,
            performance_gate_authorized = FALSE, promotion_authorized = FALSE)
        jsonlite::write_json(completion, file.path(utility_dir, "completion.json"),
            auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
        return(invisible(completion))
    }
    if (identical(options$mode, "baseline")) {
        input_bindings <- resolveCheckpointBaselineInputs(contract, options[["output-root"]])
        if (isTRUE(options[["validate-only"]])) {
            return(invisible(list(status = "eligible_no_workload", route = "iq",
                recipes = names(input_bindings), promotion_authorized = FALSE)))
        }
        if (isTRUE(options$diagnostics)) {
            stop("diagnostic profiling is not a timed baseline member", call. = FALSE)
        }
        baseline_dir <- tempfile("baseline.", tmpdir = options[["output-root"]])
        if (!dir.create(baseline_dir)) stop("cannot create owned baseline evidence directory")
        installation <- inspectCheckpointLibrary(options$library)
        protocol_hash <- digest::digest(file = options$protocol, algo = "sha256")
        contract_hash <- digest::digest(file = options$contract, algo = "sha256")
        record <- \(id, block, member) {
            path <- file.path(baseline_dir, sprintf("%s_%02d.json", id, block))
            if (file.exists(path)) stop("refusing to replace an existing baseline member receipt")
            jsonlite::write_json(list(recipe = id, block = block, route = "iq",
                installation = installation, protocol_sha256 = protocol_hash,
                contract_sha256 = contract_hash, input_binding = input_bindings[[id]],
                member = member, promotion_authorized = FALSE), path,
                auto_unbox = TRUE, pretty = TRUE, digits = NA)
        }
        run <- \(input, block) runCheckpointBaselineMember(input, block,
            options$library, options[["output-root"]])
        result <- dispatchCheckpointBaseline(input_bindings, protocol, run, record)
        jsonlite::write_json(list(status = "all_members_recorded_unqualified",
            count = length(result), route = "iq", promotion_authorized = FALSE),
            file.path(baseline_dir, "completion.json"), auto_unbox = TRUE, pretty = TRUE)
        return(invisible(result))
    }
    stop(paste("mode implementation not yet available:", options$mode,
        "(DIA-003 in progress; no workload executed)"), call. = FALSE)
}

if (sys.nframe() == 0L) checkpointProofMain(commandArgs(trailingOnly = TRUE))
