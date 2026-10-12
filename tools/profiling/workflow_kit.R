# KIT series harness: bindings, fixed member cases (journey, operation, lifecycle), journey step map,
# phase measures, frozen decisions, cross-run canonical store fingerprints and the independent
# verifier (KIT-001; operation members, the baseline journey pair, the production-default baseline
# and the reproducibility check in KIT-005).
# Pure functions only; execution lives in workflow_kit_coordinator.R and workflow_kit_driver.R.
# The four CORE-pinned harness files are never sourced for logic here, only for the browser helpers.

#' Ordered DIA journey steps in production tab order (audit B33).
#' `state` is the state-manager name the step leaves current (NULL for store-only steps);
#' `on_disk_from` is the first step after which the adopted store holds the state.
#' @noRd
kitJourneySteps <- function() {
    step <- function(id, label, state, kind, module, action) list(id = id, label = label,
        state = state, kind = kind, module = module, action = action)
    list(
        step("import", "Import", "initial", "store", "mod_prot_import", "prot-import-process"),
        step("design", "Design", "raw_data_s4", "store", "mod_prot_design", "builder save_results"),
        step("qvalue", "Q-Value Filter", "qvalue_filtered", "peptide", "mod_prot_qc_peptide_qvalue", "qvalue_filter-apply_qvalue_filter"),
        step("protein_peptide", "Protein Evidence", "protein_peptide_filtered", "peptide", "mod_prot_qc_peptide_protein", "protein_peptide_filter-apply_protein_peptide_filter"),
        step("rollup", "Precursor Rollup", "precursor_rollup", "peptide", "mod_prot_qc_peptide_rollup", "rollup-apply_rollup"),
        step("intensity", "Intensity Filter", "intensity_filtered", "peptide", "mod_prot_qc_peptide_intensity", "intensity_filter-apply_intensity_filter"),
        step("sample", "Sample Quality", "sample_filtered", "peptide", "mod_prot_qc_peptide_sample", "sample_filter-apply_sample_filter"),
        step("replicate", "Replicate Filter", "replicate_filtered", "peptide", "mod_prot_qc_peptide_replicate", "replicate_filter-apply_replicate_filter"),
        step("imputation", "Imputation", "imputed", "peptide", "mod_prot_qc_peptide_impute", "imputation-apply_imputation"),
        step("iq", "IQ Protein Rollup", "protein_s4_created", "protein", "mod_prot_qc_protein_rollup", "rollup-apply_iq_rollup"),
        step("accession", "Accession Cleanup", "protein_accession_cleaned", "protein", "mod_prot_qc_protein_cleanup", "cleanup-apply_accession_cleanup"),
        step("protein_intensity", "Protein Intensity Filter", "protein_intensity_filtered", "protein", "mod_prot_qc_protein_intensity", "intensity_filter-apply_protein_intensity_filter"),
        step("duplicates", "Duplicate Removal", "duplicates_removed", "protein", "mod_prot_qc_protein_dedup", "duplicate_removal-apply_duplicate_removal"),
        step("protein_replicate", "Protein Replicate Filter", "protein_replicate_filtered", "protein", "mod_prot_qc_protein_replicate", "replicate_filter-apply_protein_replicate_filter"),
        step("normalisation_ruv", "Normalisation and RUV", "ruv_corrected", "normalisation", "mod_prot_norm", "run_normalization"),
        step("correlation", "Correlation Filtering", "correlation_filtered", "normalisation", "mod_prot_norm", "apply_correlation_filter"),
        step("session_export", "Filtered session export", NULL, "store", "mod_prot_norm", "export_filtered_session"),
        step("da", "Differential abundance", "correlation_filtered", "da", "mod_prot_da", "run_da_analysis"))
}

#' Journey step identifiers in order.
#' @noRd
kitStepIds <- function() vapply(kitJourneySteps(), `[[`, character(1), "id")

#' Member ceilings in seconds (plan run ledger): lifecycle members 5 min; journey members from the bindings.
#' @noRd
kitMemberSeconds <- function(stage, bindings = NULL) {
    if (grepl("^kit-lifecycle-operation-", stage)) return(as.integer(bindings$operation_lifecycle$ceiling_seconds %||% 300L))
    if (grepl("^kit-lifecycle-", stage)) return(as.integer(bindings$lifecycle$ceiling_seconds %||% 300L))
    if (grepl("^kit-0[0-9]{2}$", stage)) return(600L)
    if (grepl("^kit-gui-0[0-9]{2}$", stage)) return(as.integer(bindings$tickets[[sub("kit-gui-", "", stage)]]$ceiling_seconds %||% 600L))
    if (startsWith(stage, "kit-journey")) {
        seconds <- bindings$journey$ceiling_seconds
        if (!is.numeric(seconds) || length(seconds) != 1L || !is.finite(seconds) || seconds <= 0 || seconds > 14400)
            stop("the journey ceiling is not fixed in the bindings (at most 4 h, protocol v1.2)")
        return(as.integer(seconds))
    }
    600L
}

#' Refuse a counted member unless the host meets the plan's start precondition (as CORE).
#' @noRd
kitSetupRefusal <- function(host) {
    if (!is.numeric(host$available) || length(host$available) != 1L || !is.finite(host$available))
        return("host memory is unobserved")
    if (host$available < 12 * 1024^3) return(sprintf("host MemAvailable %.2f GiB is below 12 GiB",
        host$available / 1024^3))
    if (length(host$units)) return(paste(length(host$units), "other managed unit(s) active"))
    NULL
}

#' Read the tracked KIT bindings.
#' @param path `docs/qa/workflow-execution/kit/bindings.json`.
#' @noRd
kitBindings <- function(path) {
    raw <- jsonlite::read_json(path, simplifyVector = FALSE)
    stopifnot(identical(raw$schema, 1L), is.list(raw$harness), is.list(raw$builds), is.list(raw$corpora),
        is.list(raw$journey), is.character(raw$harness$protocol_sha256))
    list(root = raw$harness$root, protocol_sha256 = raw$harness$protocol_sha256,
        step_map_sha256 = raw$harness$step_map_sha256, files = raw$harness$files,
        builds = raw$builds, outcomes = raw$outcomes, control = raw$control,
        corpora = raw$corpora, journey = raw$journey, lifecycle = raw$lifecycle,
        operation_lifecycle = raw$operation_lifecycle, tickets = raw$tickets,
        populations = raw$populations, oracles = raw$oracles)
}

#' The control build: the latest accepted KIT build, else the bound `control` (core_018 until KIT-006).
#' @noRd
kitControl <- function(bindings) {
    accepted <- Filter(\(record) identical(record$outcome, "accepted"), bindings$outcomes)
    if (length(accepted)) return(accepted[[length(accepted)]]$build)
    if (is.character(bindings$control) && !is.null(bindings$builds[[bindings$control]])) return(bindings$control)
    stop("no control build is bound")
}

#' Expand one KIT stage into its fixed member cases.
#' @param stage `kit-lifecycle-complete|fail` (journey shape, KIT-001), `kit-lifecycle-operation-complete|fail`
#'   (operation shape, KIT-005), `kit-journey` (the baseline pair and the production-default baseline, or a
#'   control/candidate pair) or `kit-0NN` (operation members of ticket NN from `bindings$tickets`).
#' @return Fixed cases in launch order. Every case carries `kind` (`journey` or `operation`), `selector`
#'   and a `binding` digest of everything that fixes it.
#' @noRd
kitCases <- function(bindings, stage) {
    stopifnot(grepl("^kit-(lifecycle-(operation-)?(complete|fail)|journey(-0[0-9]{2})?|gui-0[0-9]{2}|0[0-9]{2})$", stage))
    build_of <- function(build_key) {
        build <- bindings$builds[[build_key]]
        if (is.null(build)) stop("build ", build_key, " is not bound")
        build
    }
    journey_case <- function(id, role, rep, build_key, corpus_key, fault, selector, record = TRUE) {
        build <- build_of(build_key)
        corpus <- bindings$corpora[[corpus_key]]
        if (is.null(corpus)) stop("corpus ", corpus_key, " is not bound")
        list(id = id, kind = "journey", role = role, control = role == "control", rep = rep, fault = fault,
            corpus = corpus_key, record = record, library = build$library, expected_build = build$build,
            exact_sum = build$exact_sum, expect_failure = fault != "none",
            binding = digest::digest(list(corpus = corpus$files, build = build, stage = stage, rep = rep, record = record,
                journey = bindings$journey[setdiff(names(bindings$journey), "pair")], protocol = bindings$protocol_sha256,
                step_map = bindings$step_map_sha256), algo = "sha256"),
            selector = selector)
    }
    operation_case <- function(id, role, rep, build_key, operation, population_key, fault, selector) {
        build <- build_of(build_key)
        population <- bindings$populations[[population_key]]
        if (is.null(population)) stop("population ", population_key, " is not bound")
        parameters <- population$parameters
        if (!length(parameters)) parameters <- list()
        # The bindings are read without simplification, so a JSON array arrives as a list of scalars;
        # the pinned action validator wants atomic vectors (KIT-006 attempt 1: "inclusion entries must be
        # character" refused both sample-filter members before any work).
        parameters <- lapply(parameters, function(value)
            if (is.list(value) && length(value) && all(lengths(value) == 1L) && all(vapply(value, is.atomic, logical(1))))
                unlist(value, use.names = FALSE) else value)
        list(id = id, kind = "operation", role = role, control = role == "control", rep = rep, fault = fault,
            operation = operation, population = population_key, parameters = parameters,
            input = list(store = population$store, head_sha256 = population$head_sha256, rows = population$rows,
                entity_id = population$entity_id, head_entity = population$head_entity, state = population$state),
            library = build$library, expected_build = build$build, exact_sum = build$exact_sum,
            expect_failure = fault != "none",
            binding = digest::digest(list(population = population, build = build, stage = stage, rep = rep,
                operation = operation, protocol = bindings$protocol_sha256), algo = "sha256"),
            selector = selector)
    }
    if (grepl("^kit-lifecycle-operation-", stage)) {
        kind <- sub("kit-lifecycle-operation-", "", stage, fixed = TRUE)
        fault <- bindings$operation_lifecycle$fault
        if (is.null(fault)) stop("the operation lifecycle fault is not bound")
        # The lifecycle pair proves a changed launch path on the build that changed it (KIT-007: `operation_lifecycle$build`);
        # without a named build it runs on the control, as KIT-005 did.
        return(list(operation_case("member", "candidate", 1L, bindings$operation_lifecycle$build %||% kitControl(bindings),
            fault$operation, fault$population,
            fault = if (kind == "fail") fault$name else "none", selector = paste0("kit:lifecycle:operation:", kind))))
    }
    if (startsWith(stage, "kit-lifecycle-")) {
        kind <- sub("kit-lifecycle-", "", stage, fixed = TRUE)
        # `lifecycle$build` names the build whose launch path changed (KIT-011: the production default); else the control.
        return(list(journey_case("member", "candidate", 1L, bindings$lifecycle$build %||% kitControl(bindings), "screen",
            fault = if (kind == "fail") bindings$lifecycle$fault$name else "none",
            selector = paste0("kit:lifecycle:", kind))))
    }
    if (stage == "kit-journey") {
        spec <- bindings$journey$pair
        if (is.null(spec)) stop("the journey pair is defined by KIT-005")
        reps <- as.integer(spec$reps %||% 1L)
        stopifnot(reps %in% 1:2)
        cases <- list()
        for (rep in seq_len(reps)) {
            cases[[length(cases) + 1L]] <- journey_case(paste0("control-", rep), "control", rep, spec$control, "full", "none",
                paste0("kit:journey:control:", rep))
            if (!is.null(spec$candidate)) cases[[length(cases) + 1L]] <- journey_case(paste0("candidate-", rep), "candidate",
                rep, spec$candidate, "full", "none", paste0("kit:journey:candidate:", rep))
        }
        if (isTRUE(spec$production_default)) cases[[length(cases) + 1L]] <- journey_case("default-1", "default", 1L,
            spec$control, "full", "none", "kit:journey:default:1", record = FALSE)
        return(cases)
    }
    if (grepl("^kit-journey-0[0-9]{2}$", stage)) {
        # A ticket's journey pair (KIT-011): control and candidate builds run the whole journey as users run
        # it (`record = FALSE`: no launch record, the build's own production default), one member each.
        ticket <- sub("kit-journey-", "", stage, fixed = TRUE)
        spec <- bindings$tickets[[ticket]]
        if (is.null(spec) || !identical(spec$kind, "journey")) stop("KIT-", ticket, " has no bound journey pair")
        record <- isTRUE(spec$record)
        return(list(
            journey_case("control-1", "control", 1L, spec$control, "full", "none", paste0("kit:journey-", ticket, ":control:1"), record = record),
            journey_case("candidate-1", "candidate", 1L, spec$candidate, "full", "none", paste0("kit:journey-", ticket, ":candidate:1"), record = record)))
    }
    if (grepl("^kit-gui-0[0-9]{2}$", stage)) {
        # A ticket's `gui`-mode managed member (KIT-008): one journey on the bound corpus and build, launched
        # with no record (the build's production default under a managed application), a candidate only.
        ticket <- sub("kit-gui-", "", stage, fixed = TRUE)
        spec <- bindings$tickets[[ticket]]
        if (is.null(spec) || !identical(spec$kind, "gui")) stop("KIT-", ticket, " has no bound gui member")
        return(list(journey_case("member", "candidate", 1L, spec$build, spec$corpus %||% "screen", "none",
            paste0("kit:gui-", ticket, ":candidate:1"), record = FALSE)))
    }
    ticket <- sub("kit-", "", stage, fixed = TRUE)
    spec <- bindings$tickets[[ticket]]
    if (is.null(spec)) stop("KIT-", ticket, " has no bound operation, population, control and candidate builds")
    # A ticket binds one or two pairs (operation, population); pair i is members control-i and
    # candidate-i, control first then candidate, pair after pair (KIT-006: full_p14 then full_primary).
    pairs <- spec$pairs %||% list(list(operation = spec$operation, population = spec$population))
    stopifnot(length(pairs) >= 1L, length(pairs) <= 2L)
    cases <- list()
    for (i in seq_along(pairs)) for (role in c("control", "candidate")) {
        cases[[length(cases) + 1L]] <- operation_case(paste(role, i, sep = "-"), role, i, spec[[role]],
            pairs[[i]]$operation, pairs[[i]]$population, "none", paste("kit", ticket, role, i, sep = ":"))
    }
    cases
}

#' Decide a ticket's operation pairs (protocol section 3): every member eligible and verified exact
#' against its oracle; per pair, the candidate's kernel peak and whole-slice anonymous peak are not
#' worse than 1.02 x the control's. A reduction is reported (the 0.95 rule) and never required here.
#' @return `status` (`accepted`, `rejected` or `invalid`), per-pair ratios and the failing rule.
#' @noRd
kitOperationDecision <- function(receipt, cases, bindings) {
    seconds <- kitMemberSeconds(paste0("kit-", strsplit(cases[[1L]]$selector, ":", fixed = TRUE)[[1L]][[2L]]), bindings)
    members <- list()
    for (case in cases) {
        item <- receipt$cases[[case$id]]
        op <- item$operation
        members[[case$id]] <- list(role = case$role, pair = case$rep, operation = case$operation, population = case$population,
            failing_rule = kitIneligible(item, seconds, "operation"), peak_bytes = op$peak_bytes,
            peak_anon_bytes = op$sampled_anon_bytes, elapsed_seconds = op$elapsed_seconds,
            batches = op$science_execution$admission$batches %||%
                (if (length(op$science_execution$batch_sizes)) length(op$science_execution$batch_sizes) else NULL) %||%
                (if (length(op$science_execution$admission$count_batch_sizes)) length(op$science_execution$admission$count_batch_sizes) else NULL))
    }
    invalid <- Filter(\(m) !is.null(m$failing_rule), members)
    pairs <- list()
    for (i in sort(unique(vapply(members, `[[`, integer(1), "pair")))) {
        control <- members[[paste0("control-", i)]]
        candidate <- members[[paste0("candidate-", i)]]
        if (is.null(control) || is.null(candidate)) next
        ratio <- function(name) candidate[[name]] / control[[name]]
        pairs[[as.character(i)]] <- list(operation = control$operation, population = control$population,
            peak_ratio = ratio("peak_bytes"), peak_anon_ratio = ratio("peak_anon_bytes"), elapsed_ratio = ratio("elapsed_seconds"),
            not_worse = is.finite(ratio("peak_bytes")) && ratio("peak_bytes") <= 1.02 &&
                is.finite(ratio("peak_anon_bytes")) && ratio("peak_anon_bytes") <= 1.02,
            reduction = is.finite(ratio("peak_bytes")) && ratio("peak_bytes") <= 0.95,
            batches = list(control = control$batches, candidate = candidate$batches))
    }
    status <- if (length(invalid)) "invalid" else if (!length(pairs) || !all(vapply(pairs, `[[`, logical(1), "not_worse"))) "rejected" else "accepted"
    list(status = status, members = members, pairs = pairs,
        failing_rule = if (length(invalid)) unlist(lapply(invalid, `[[`, "failing_rule")) else
            if (status == "rejected") "not_worse" else NULL,
        promotion_authorized = FALSE)
}

#' Map a verifier selector back to its stage and case id.
#' @noRd
kitSelectorStage <- function(selection) {
    stopifnot(grepl("^kit:(journey:((control|candidate):[12]|default:1)|journey-0[0-9]{2}:(control|candidate):1|gui-0[0-9]{2}:candidate:1|lifecycle:(operation:)?(complete|fail)|0[0-9]{2}:(control|candidate):[12])$", selection))
    parts <- strsplit(sub("^kit:", "", selection), ":", fixed = TRUE)[[1L]]
    if (parts[[1L]] == "lifecycle") return(list(stage = paste0("kit-lifecycle-", paste(parts[-1L], collapse = "-")), id = "member"))
    # A ticket's gui member is the single case `member` (KIT-008 correction 1: the verifier looked for `candidate-1`).
    if (startsWith(parts[[1L]], "gui-")) return(list(stage = paste0("kit-", parts[[1L]]), id = "member"))
    list(stage = paste0("kit-", parts[[1L]]), id = paste(parts[[2L]], parts[[3L]], sep = "-"))
}

#' Turn the driver's step marks and the coordinator's 50 ms samples into per-step measures.
#' @param phases Data frame `uptime, step, event` (`start`, then `complete` or `failed`; one `teardown` mark).
#' @param trace Data frame `uptime, current, anon` of the coordinator's samples.
#' @return `phase_seconds`, `phase_peak_bytes`, `phase_peak_anon_bytes`, `phase_status`
#'   (`failed` > `open` > `complete`) and `resident_end_bytes` (minimum anonymous memory in the
#'   5 s before the teardown mark; NA without a mark or without samples).
#' @noRd
kitPhaseMeasures <- function(phases, trace) {
    stopifnot(is.data.frame(phases), all(c("uptime", "step", "event") %in% names(phases)),
        is.data.frame(trace), all(c("uptime", "current", "anon") %in% names(trace)))
    starts <- phases[phases$event == "start", , drop = FALSE]
    ends <- phases[phases$event %in% c("complete", "failed"), , drop = FALSE]
    labels <- unique(starts$step)
    seconds <- peaks <- anon <- stats::setNames(rep(NA_real_, length(labels)), labels)
    status <- stats::setNames(rep("complete", length(labels)), labels)
    last <- if (nrow(trace)) max(trace$uptime) else max(c(phases$uptime, -Inf))
    for (i in seq_len(nrow(starts))) {
        label <- starts$step[[i]]
        end <- ends[ends$step == label & ends$uptime >= starts$uptime[[i]], , drop = FALSE]
        closed <- nrow(end) >= 1L
        finish <- if (closed) min(end$uptime) else last
        state <- if (closed) end$event[which.min(end$uptime)] else "open"
        if (state == "failed" || state == "open" && status[[label]] != "failed") status[[label]] <- state
        seconds[[label]] <- sum(seconds[[label]], finish - starts$uptime[[i]], na.rm = TRUE)
        rows <- trace$uptime >= starts$uptime[[i]] & trace$uptime <= finish
        if (any(rows)) {
            peaks[[label]] <- max(peaks[[label]], trace$current[rows], na.rm = TRUE)
            anon[[label]] <- max(anon[[label]], trace$anon[rows], na.rm = TRUE)
        }
    }
    teardown <- phases$uptime[phases$event == "teardown"]
    resident <- NA_real_
    if (length(teardown) == 1L && nrow(trace)) {
        window <- trace$uptime >= teardown - 5 & trace$uptime <= teardown
        if (any(window)) resident <- min(trace$anon[window])
    }
    list(phase_seconds = as.list(seconds), phase_peak_bytes = as.list(peaks),
        phase_peak_anon_bytes = as.list(anon), phase_status = as.list(status),
        resident_end_bytes = resident)
}

#' Measurement validity: coordinator sampling gap at most 250 ms (operation members) or 500 ms
#' (journey members, protocol v1.3: the whole-slice peak is the kernel high-water mark; a gap only
#' blinds the sampled anonymous trace for that window), zero swap, slice stopped.
#' @noRd
kitMeasurementValid <- function(op, kind = "operation") {
    gap <- op$max_sample_interval
    limit <- if (identical(kind, "journey")) 0.5 else 0.25
    is.numeric(gap) && length(gap) == 1L && is.finite(gap) && gap <= limit &&
        is.numeric(op$peak_bytes) && length(op$peak_bytes) == 1L && is.finite(op$peak_bytes) &&
        identical(op$swap_bytes, 0) && isTRUE(op$cleanup$slice_stopped)
}

#' Check one member against the frozen eligibility rules.
#' @return Name of the first failing rule, or NULL when eligible.
#' @noRd
kitIneligible <- function(item, seconds, kind = "operation") {
    op <- item$operation
    rules <- list(
        verifier_accepted = isTRUE(item$verifier_completion$accepted),
        status_complete = identical(op$status, "complete"),
        measurement_valid = kitMeasurementValid(op, kind),
        zero_swap = identical(op$swap_bytes, 0),
        peak_below_early_stop = is.numeric(op$peak_bytes) && length(op$peak_bytes) == 1L &&
            is.finite(op$peak_bytes) && op$peak_bytes > 0 && op$peak_bytes < 15.5 * 1024^3,
        elapsed_within_ceiling = is.numeric(op$elapsed_seconds) && length(op$elapsed_seconds) == 1L &&
            is.finite(op$elapsed_seconds) && op$elapsed_seconds > 0 && op$elapsed_seconds <= seconds)
    failing <- names(rules)[!vapply(rules, isTRUE, logical(1))]
    if (length(failing)) failing[[1L]] else NULL
}

#' Protocol guards the sampler or the application monitor can end a member with.
#' @noRd
kitGuardReasons <- c("memory_guard", "host_guard", "time_guard")

#' The production-default member's outcome: `complete`, the guard that ended it, or NULL.
#' @noRd
kitDefaultOutcome <- function(op) {
    if (identical(op$status, "complete")) return("complete")
    guard <- op$outer_guard %||% op$guard$reason
    if (is.character(guard) && length(guard) == 1L && guard %in% kitGuardReasons) guard else NULL
}

#' Eligibility of the production-default baseline member (KIT-005 acceptance as amended by the user
#' 2026-10-11): the member completes, or ends at a protocol guard, with its per-step anonymous peaks
#' recorded. A completed member follows the journey rules; a guarded member needs a valid measurement
#' (sampling gaps, zero swap, slice stopped) and at least one completed step with a finite peak. Any
#' other end (a driver failure, an application error, a hand stop) is not a recorded baseline.
#' @return Name of the first failing rule, or NULL when eligible.
#' @noRd
kitDefaultIneligible <- function(item, seconds) {
    op <- item$operation
    outcome <- kitDefaultOutcome(op)
    if (identical(outcome, "complete")) return(kitIneligible(item, seconds, "journey"))
    if (is.null(outcome)) return("status_complete_or_guard")
    rules <- list(measurement_valid = kitMeasurementValid(op, "journey"),
        zero_swap = identical(op$swap_bytes, 0),
        phase_peaks_recorded = {
            peaks <- unlist(op$phase_peak_anon_bytes[names(op$phase_status)[unlist(op$phase_status) == "complete"]])
            length(peaks) >= 1L && all(is.finite(peaks))
        })
    failing <- names(rules)[!vapply(rules, isTRUE, logical(1))]
    if (length(failing)) failing[[1L]] else NULL
}

#' Classify a lifecycle member; a failure member is never relabelled as success.
#' @param receipt,cases Stage receipt and its single fixed case.
#' @param bindings KIT bindings (the fault names the step that must fail and the last good step).
#' @return `lifecycle_complete_observed`, `lifecycle_failure_observed` or `invalid` with the failing rule.
#' @noRd
kitLifecycleDecision <- function(receipt, cases, bindings) {
    stopifnot(length(cases) == 1L)
    case <- cases[[1L]]
    item <- receipt$cases[[case$id]]
    op <- item$operation
    result <- list(status = "invalid", failing_rule = NULL, promotion_authorized = FALSE,
        expect_failure = case$expect_failure)
    status <- op$phase_status
    steps <- kitStepIds()
    if (!case$expect_failure) {
        result$failing_rule <- kitIneligible(item, kitMemberSeconds("kit-lifecycle-complete", bindings), "journey")
        if (is.null(result$failing_rule) && !(all(steps %in% names(status)) &&
            all(unlist(status[steps]) == "complete") && all(is.finite(unlist(op$phase_peak_anon_bytes[steps])))))
            result$failing_rule <- "phase_timing"
        if (is.null(result$failing_rule) && !(is.numeric(op$resident_end_bytes) && is.finite(op$resident_end_bytes)))
            result$failing_rule <- "resident_end_unobserved"
        if (is.null(result$failing_rule) && !isTRUE(op$populations_complete)) result$failing_rule <- "populations_incomplete"
        if (is.null(result$failing_rule)) result$status <- "lifecycle_complete_observed"
        return(result)
    }
    fault <- bindings$lifecycle$fault
    rules <- list(verifier_accepted = isTRUE(item$verifier_completion$accepted),
        failure_recorded = identical(op$status, "journey_failed"),
        zero_swap = identical(op$swap_bytes, 0),
        cleanup = isTRUE(op$cleanup$slice_stopped),
        guard_observed = identical(op$guard$reason, "injected_guard_stop"),
        phase_timing = identical(status[[fault$failed_step]], "failed") &&
            identical(status[[fault$after_step]], "complete"),
        head_unchanged = isTRUE(op$head_unchanged))
    failing <- names(rules)[!vapply(rules, isTRUE, logical(1))]
    if (length(failing)) result$failing_rule <- failing[[1L]] else result$status <- "lifecycle_failure_observed"
    result
}

#' P14 phases every completed operation member must record (CORE protocol v3, carried).
#' @noRd
kitRequiredOperationPhases <- c("relational_validation", "working_memory_admission",
    "relational_aggregate_order_summary_write", "ordered_result_stage")

#' Turn worker phase events (`uptime, id, phase, event`) and the coordinator's samples into phase
#' measures for an operation member (the CORE measure, carried unchanged into the kit).
#' @noRd
kitOperationPhaseMeasures <- function(events, trace) {
    stopifnot(is.data.frame(events), all(c("uptime", "id", "phase", "event") %in% names(events)),
        is.data.frame(trace), all(c("uptime", "current") %in% names(trace)))
    starts <- events[events$event == "start", , drop = FALSE]
    ends <- events[events$event %in% c("complete", "failed"), , drop = FALSE]
    labels <- unique(starts$phase)
    seconds <- peaks <- anon <- stats::setNames(rep(NA_real_, length(labels)), labels)
    status <- stats::setNames(rep("complete", length(labels)), labels)
    last <- if (nrow(trace)) max(trace$uptime) else max(c(events$uptime, -Inf))
    for (i in seq_len(nrow(starts))) {
        label <- starts$phase[[i]]
        end <- ends[ends$id == starts$id[[i]], , drop = FALSE]
        closed <- nrow(end) == 1L
        finish <- if (closed) end$uptime[[1L]] else last
        state <- if (closed) end$event[[1L]] else "open"
        if (state == "failed" || state == "open" && status[[label]] != "failed") status[[label]] <- state
        seconds[[label]] <- sum(seconds[[label]], finish - starts$uptime[[i]], na.rm = TRUE)
        rows <- trace$uptime >= starts$uptime[[i]] & trace$uptime <= finish
        if (any(rows)) {
            peaks[[label]] <- max(peaks[[label]], trace$current[rows], na.rm = TRUE)
            if (!is.null(trace$anon)) anon[[label]] <- max(anon[[label]], trace$anon[rows], na.rm = TRUE)
        }
    }
    list(phase_seconds = as.list(seconds), phase_peak_bytes = as.list(peaks),
        phase_peak_anon_bytes = as.list(anon), phase_status = as.list(status))
}

#' Classify an operation lifecycle member (operation-member shape; KIT-005): completion through the
#' installed action dispatch with every P14 phase complete, or the mid-step forced failure
#' (`after_aggregate`: the worker fails inside the ordered result stage) with the head unchanged.
#' @noRd
kitOperationLifecycleDecision <- function(receipt, cases, bindings) {
    stopifnot(length(cases) == 1L)
    case <- cases[[1L]]
    item <- receipt$cases[[case$id]]
    op <- item$operation
    result <- list(status = "invalid", failing_rule = NULL, promotion_authorized = FALSE,
        expect_failure = case$expect_failure)
    status <- op$phase_status
    if (!case$expect_failure) {
        result$failing_rule <- kitIneligible(item, kitMemberSeconds("kit-lifecycle-operation-complete", bindings))
        if (is.null(result$failing_rule) && !(all(kitRequiredOperationPhases %in% names(status)) &&
            all(unlist(status[kitRequiredOperationPhases]) == "complete")))
            result$failing_rule <- "phase_timing"
        if (is.null(result$failing_rule)) result$status <- "lifecycle_complete_observed"
        return(result)
    }
    fault <- bindings$operation_lifecycle$fault
    rules <- list(verifier_accepted = isTRUE(item$verifier_completion$accepted),
        failure_recorded = identical(op$status, "worker_failed"),
        zero_swap = identical(op$swap_bytes, 0),
        cleanup = isTRUE(op$cleanup$slice_stopped),
        phase_timing = identical(status[[fault$failed_phase]], "failed"),
        head_unchanged = isTRUE(op$head_unchanged))
    failing <- names(rules)[!vapply(rules, isTRUE, logical(1))]
    if (length(failing)) result$failing_rule <- failing[[1L]] else result$status <- "lifecycle_failure_observed"
    result
}

#' Compare the per-step fingerprints of two journey members (the oracle reproducibility check).
#' @return Character vector of `step/kind/name` mismatches; empty when every fingerprint reproduces.
#' @noRd
#' The products a declared operation saves under its own scientific identity (KIT-004's runner): their
#' name suffix and content embed the build's software identity, so across builds they differ by
#' construction while the states they reference are identical. A cross-build check reports them
#' (`kitOperationIdentityProducts`) and judges everything else.
#' @noRd
kitOperationIdentityProducts <- "^dia-(normalization-ruv|pearson-pairs|pearson-selection|contrast-references):[0-9a-f]{12}$"

kitReproduceMismatches <- function(reference, candidate, ignore_products = NULL) {
    mismatches <- character()
    keep <- function(names) if (is.null(ignore_products)) names else names[!grepl(ignore_products, names)]
    for (step in names(reference)) {
        expected <- reference[[step]]$fingerprint
        if (is.null(expected)) next
        actual <- candidate[[step]]$fingerprint
        if (is.null(actual)) { mismatches <- c(mismatches, paste0(step, "/missing")); next }
        for (kind in c("components", "support", "products")) {
            names_expected <- names(expected[[kind]])
            names_actual <- names(actual[[kind]])
            if (kind == "products") { names_expected <- keep(names_expected); names_actual <- keep(names_actual) }
            if (!setequal(names_expected, names_actual))
                mismatches <- c(mismatches, paste0(step, "/", kind, "/names"))
            for (name in intersect(names_expected, names_actual)) {
                if (!identical(expected[[kind]][[name]]$sha256, actual[[kind]][[name]]$sha256))
                    mismatches <- c(mismatches, paste0(step, "/", kind, "/", name))
            }
        }
        if (!identical(expected$summary_sha256, actual$summary_sha256)) mismatches <- c(mismatches, paste0(step, "/summary"))
        if (!identical(expected$rows, actual$rows)) mismatches <- c(mismatches, paste0(step, "/rows"))
    }
    mismatches
}

#' Journey stage decision (KIT-005): the baseline pair is recorded when both control members are
#' eligible, verified and the second reproduces the first's fingerprints; the production-default
#' baseline is recorded when it completes with its measures. Nothing is promoted.
#' @noRd
kitJourneyDecision <- function(receipt, cases, bindings) {
    seconds <- kitMemberSeconds("kit-journey", bindings)
    members <- list()
    for (case in cases) {
        item <- receipt$cases[[case$id]]
        op <- item$operation
        failing <- if (case$role == "default") kitDefaultIneligible(item, seconds) else kitIneligible(item, seconds, "journey")
        if (is.null(failing) && case$role != "default" && !isTRUE(op$populations_complete)) failing <- "populations_incomplete"
        if (is.null(failing) && case$rep == 2L && !isTRUE(item$verifier_completion$reproduced)) failing <- "not_reproduced"
        members[[case$id]] <- list(role = case$role, rep = case$rep, eligible = is.null(failing), failing_rule = failing,
            outcome = if (case$role == "default") kitDefaultOutcome(op) else NULL,
            max_sample_interval = op$max_sample_interval, gaps_above_250ms = op$gaps_above_250ms,
            peak_anon_bytes = op$sampled_anon_bytes, peak_bytes = op$peak_bytes,
            resident_end_bytes = op$resident_end_bytes, elapsed_seconds = op$elapsed_seconds,
            phase_peak_anon_bytes = op$phase_peak_anon_bytes, phase_seconds = op$phase_seconds)
    }
    controls <- Filter(\(m) m$role == "control", members)
    noise <- NULL
    if (length(controls) == 2L && all(vapply(controls, `[[`, logical(1), "eligible"))) {
        ratio <- function(name) controls[[2L]][[name]] / controls[[1L]][[name]]
        noise <- list(peak_anon = ratio("peak_anon_bytes"), resident_end = ratio("resident_end_bytes"),
            elapsed = ratio("elapsed_seconds"))
    }
    default <- Filter(\(m) m$role == "default", members)
    list(status = if (length(controls) && all(vapply(controls, `[[`, logical(1), "eligible"))) "baseline_recorded" else "invalid",
        failing_rule = unlist(lapply(members, `[[`, "failing_rule")), members = members, noise = noise,
        production_default = if (length(default)) list(recorded = default[[1L]]$eligible, failing_rule = default[[1L]]$failing_rule,
            outcome = default[[1L]]$outcome) else NULL,
        promotion_authorized = FALSE)
}

#' The KIT-005 per-step oracles as a populations list, for `kitReproduceMismatches` against a member's
#' in-run populations.
#' @noRd
kitOracleJourneyPopulations <- function(bindings) {
    steps <- kitStepIds()
    stats::setNames(lapply(steps, function(step) {
        oracle <- bindings$oracles[[paste0("journey_", step)]]
        if (is.null(oracle) || is.null(oracle$class) || !length(oracle$components)) return(list(fingerprint = NULL))
        list(fingerprint = oracle)
    }), steps)
}

#' Decision for a ticket's journey pair (KIT-011; protocol section 3 and the KIT-005 acceptance as amended
#' by the user): the control is the build users run today and may end at a protocol guard; the candidate
#' must complete the whole journey with its states verified against the KIT-005 oracles. The required
#' reduction (0.95) is judged on the whole-slice anonymous peak: when the control completed, against its
#' peak; when the control ended at a guard, against the anonymous peak it had reached, which is a lower
#' bound of its true peak, so the test is only stricter (`control_peak_is_lower_bound`). Per step, on
#' every step the control completed, the candidate may not be worse than 1.02 x the control (the per-step
#' guard of KIT-017); the ratios are reported. `resident_end_bytes` is reported for both and ratioed when
#' both completed.
#' @return `status` (`accepted`, `rejected` or `invalid`), per-step ratios and the failing rule.
#' @noRd
kitJourneyTicketDecision <- function(receipt, cases, bindings) {
    seconds <- kitMemberSeconds(paste0("kit-", strsplit(cases[[1L]]$selector, ":", fixed = TRUE)[[1L]][[2L]]), bindings)
    members <- list()
    for (case in cases) {
        item <- receipt$cases[[case$id]]
        op <- item$operation
        failing <- if (case$role == "control") kitDefaultIneligible(item, seconds) else kitIneligible(item, seconds, "journey")
        if (is.null(failing) && case$role == "candidate" && !isTRUE(op$populations_complete)) failing <- "populations_incomplete"
        members[[case$id]] <- list(role = case$role, eligible = is.null(failing), failing_rule = failing,
            outcome = if (case$role == "control") kitDefaultOutcome(op) else if (identical(op$status, "complete")) "complete" else NULL,
            peak_anon_bytes = op$sampled_anon_bytes, peak_bytes = op$peak_bytes, resident_end_bytes = op$resident_end_bytes,
            elapsed_seconds = op$elapsed_seconds, max_sample_interval = op$max_sample_interval, gaps_above_250ms = op$gaps_above_250ms,
            phase_peak_anon_bytes = op$phase_peak_anon_bytes, phase_seconds = op$phase_seconds, phase_status = op$phase_status)
    }
    control <- members[["control-1"]]
    candidate <- members[["candidate-1"]]
    invalid <- unlist(lapply(members, `[[`, "failing_rule"))
    if (length(invalid)) return(list(status = "invalid", failing_rule = invalid, members = members, promotion_authorized = FALSE))
    completed <- names(control$phase_status)[unlist(control$phase_status) == "complete"]
    completed <- completed[vapply(completed, function(step) is.finite(control$phase_peak_anon_bytes[[step]] %||% NA_real_) &&
        is.finite(candidate$phase_peak_anon_bytes[[step]] %||% NA_real_), logical(1))]
    step_ratios <- vapply(completed, function(step) candidate$phase_peak_anon_bytes[[step]] / control$phase_peak_anon_bytes[[step]], numeric(1))
    # The per-step guard applies from the first saved state (q-value) onward (user decision 2026-10-11, KIT-011
    # audit PX-01): import and design run the same code on both paths and are reported, not judged.
    reported_only <- c("import", "design")
    guarded <- setdiff(completed, reported_only)
    control_complete <- identical(control$outcome, "complete")
    whole_ratio <- candidate$peak_anon_bytes / control$peak_anon_bytes
    resident_ratio <- if (control_complete && is.finite(control$resident_end_bytes) && is.finite(candidate$resident_end_bytes))
        candidate$resident_end_bytes / control$resident_end_bytes else NA_real_
    not_worse_steps <- length(guarded) >= 1L && all(step_ratios[guarded] <= 1.02)
    reduced <- is.finite(whole_ratio) && whole_ratio <= 0.95
    status <- if (!not_worse_steps) "rejected" else if (!reduced) "rejected" else "accepted"
    list(status = status, failing_rule = if (!not_worse_steps) "step_not_worse" else if (!reduced) "required_reduction" else NULL,
        members = members, control_outcome = control$outcome, control_peak_is_lower_bound = !control_complete,
        steps_compared = completed, steps_guarded = guarded, step_ratios = as.list(step_ratios),
        whole_slice_anon_ratio = whole_ratio, resident_end_ratio = resident_ratio,
        resident_end_bytes = list(control = control$resident_end_bytes, candidate = candidate$resident_end_bytes),
        promotion_authorized = FALSE)
}

#' Decide a ticket's `gui`-mode managed member (KIT-008; protocol: completion and exactness, peak recorded).
#' The member completes every journey step within its ceiling, its populations are complete and verified,
#' every state, support entry and summary reproduces the bound reference screen member (the KIT-011 screen
#' member on the control build; the per-build identity products reported, not judged), and the two
#' production dispatches ran on the managed worker (completed worker journals for both operations in
#' the member's project store).
#' @return `status` (`gui_member_observed`, `rejected` or `invalid`), the failing rule, the dispatch
#'   counts and the reproduction result.
#' @noRd
kitGuiMemberDecision <- function(receipt, cases, bindings, member_dir = NULL, api = NULL) {
    stopifnot(length(cases) == 1L)
    case <- cases[[1L]]
    ticket <- sub("^kit:gui-([0-9]{3}):.*$", "\\1", case$selector)
    spec <- bindings$tickets[[ticket]]
    item <- receipt$cases[[case$id]]
    op <- item$operation
    result <- list(status = "invalid", failing_rule = NULL, promotion_authorized = FALSE, ticket = ticket)
    result$failing_rule <- kitIneligible(item, kitMemberSeconds(paste0("kit-gui-", ticket), bindings), "journey")
    steps <- kitStepIds()
    status <- op$phase_status
    if (is.null(result$failing_rule) && !(all(steps %in% names(status)) && all(unlist(status[steps]) == "complete")))
        result$failing_rule <- "phase_timing"
    if (is.null(result$failing_rule) && !isTRUE(op$populations_complete)) result$failing_rule <- "populations_incomplete"
    result$peak_bytes <- op$peak_bytes
    result$peak_anon_bytes <- op$sampled_anon_bytes
    result$resident_end_bytes <- op$resident_end_bytes
    result$elapsed_seconds <- op$elapsed_seconds
    # Dispatch evidence: completed worker journals under the member's project store, by operation.
    dispatched <- list()
    if (!is.null(member_dir) && dir.exists(member_dir)) {
        requests <- list.files(member_dir, "^request[.]rds$", recursive = TRUE, full.names = TRUE, all.files = TRUE)
        requests <- requests[grepl("/[.]workflow-export/wfx-", requests)]
        for (path in requests) {
            request <- tryCatch(readRDS(path), error = function(error) NULL)
            operation <- request$operation %||% NA_character_
            if (is.na(operation)) next
            journal_path <- file.path(dirname(path), "journal.json")
            status <- if (file.exists(journal_path)) tryCatch(jsonlite::read_json(journal_path)$status, error = function(error) NULL) else NULL
            dispatched[[operation]] <- c(dispatched[[operation]], status %||% NA_character_)
        }
    }
    result$dispatched <- lapply(dispatched, function(statuses) list(requests = length(statuses),
        complete = sum(statuses == "complete", na.rm = TRUE)))
    required <- spec$operations %||% c("precursor_rollup", "minimum_peptides_per_sample")
    worker_ran <- all(vapply(required, function(name) isTRUE(result$dispatched[[name]]$complete >= 1L), logical(1)))
    if (is.null(result$failing_rule) && !worker_ran) result$failing_rule <- "worker_dispatch_missing"
    # Reproduction against the bound reference screen member (populations saved by its driver).
    runtime_root <- bindings$root %||% bindings$harness$root
    reference_path <- if (is.null(runtime_root)) "" else
        file.path(runtime_root, spec$reference %||% "kit-lifecycle-complete", "member", "journey-result.rds")
    candidate_path <- if (is.null(member_dir)) NULL else file.path(member_dir, "journey-result.rds")
    if (is.null(result$failing_rule) && file.exists(reference_path) && !is.null(candidate_path) && file.exists(candidate_path)) {
        reference <- readRDS(reference_path)
        candidate <- readRDS(candidate_path)
        # With an api the state fingerprints of both members are recomputed from their immutable stores under the
        # current canonical form (a correction to that form after the run then applies to both sides, as the
        # journey verifier does); the per-step assets the driver attached stay as recorded.
        if (!is.null(api)) {
            recompute <- function(journey) {
                fresh <- kitFingerprintStore(journey$store, api)
                for (step in names(journey$populations)) {
                    state <- journey$populations[[step]]$state
                    if (is.null(state) || is.null(fresh$states[[state]])) next
                    merged <- journey$populations[[step]]$fingerprint
                    for (field in c("components", "support", "summary_sha256", "rows")) merged[[field]] <- fresh$states[[state]][[field]]
                    for (name in names(fresh$states[[state]]$products)) merged$products[[name]] <- fresh$states[[state]]$products[[name]]
                    journey$populations[[step]]$fingerprint <- merged
                }
                journey
            }
            reference <- recompute(reference)
            candidate <- recompute(candidate)
            result$fingerprints <- "recomputed from both stores under the current canonical form"
        }
        mismatches <- kitReproduceMismatches(reference$populations, candidate$populations, ignore_products = kitOperationIdentityProducts)
        result$reference <- spec$reference %||% "kit-lifecycle-complete"
        result$reproduction <- list(fingerprint_mismatches = length(mismatches), mismatches = mismatches)
        if (length(mismatches)) result$failing_rule <- "reference_mismatch"
    } else if (is.null(result$failing_rule)) result$failing_rule <- "reference_unavailable"
    if (is.null(result$failing_rule)) result$status <- "gui_member_observed" else if (
        result$failing_rule %in% c("worker_dispatch_missing", "reference_mismatch")) result$status <- "rejected"
    result
}

#' Mask run-specific absolute paths so args summaries compare across members.
#' @noRd
#' Only plain lists and data frames are descended, to a bounded depth: plot objects and grob trees
#' (the RUV diagnostic product) nest far deeper and carry no run-specific paths.
#' @noRd
kitMaskPaths <- function(value, depth = 0L) {
    if (is.character(value)) {
        masked <- value
        at <- !is.na(masked) & startsWith(masked, "/")
        masked[at] <- paste0("<path>/", basename(masked[at]))
        return(masked)
    }
    if (is.list(value) && depth < 32L && (!is.object(value) || is.data.frame(value)))
        return(lapply(value, kitMaskPaths, depth = depth + 1L))
    value
}

#' Hash a table component by logical pages (every column, every row, 10,000 rows per page).
#' @noRd
kitFingerprintTable <- function(store, reference, api) {
    scope <- api("newWorkflowTableReadScope")(store, reference, reader = "fastduckdb")
    on.exit(api("releaseWorkflowTableReadScope")(scope), add = TRUE)
    hashes <- kitPageHashes(store, reference, api, scope)
    list(kind = "table", rows = scope$metadata$rows, columns = scope$metadata$columns,
        sha256 = digest::digest(paste(hashes, collapse = "\n"), algo = "sha256", serialize = FALSE))
}

#' Page hashes within an open read scope.
#' @noRd
kitPageHashes <- function(store, reference, api, scope) {
    rows <- scope$metadata$rows
    hashes <- character()
    if (rows > 0L) for (offset in seq.int(0L, rows - 1L, by = 10000L)) {
        # Wide protein tables (one column per sample) exceed the reader's default 64-column cap (KIT-005 correction 6).
        data <- api("readSavedWorkflowTable")(store, reference, scope$metadata$columns,
            offset = offset, limit = min(10000L, rows - offset), max_columns = max(64L, length(scope$metadata$columns)),
            scope = scope)$data
        attr(data, "dia_source_order") <- NULL
        attr(data, "row.names") <- .set_row_names(nrow(data))
        hashes[[length(hashes) + 1L]] <- digest::digest(data, algo = "sha256")
    }
    hashes
}

#' Canonical form of a hydrated value for fingerprints that must reproduce across runs: absolute
#' paths masked, `timestamp`/`mtime` fields dropped, saved references replaced by their payload
#' hash (tables) or their own canonical form (records), plot and grob objects by their class.
#' @noRd
kitCanonical <- function(value, store = NULL, api = NULL, depth = 0L, memo = NULL) {
    if (depth > 48L) return("<depth>")
    if (is.null(value) || is.function(value) || is.environment(value)) return(if (is.null(value)) NULL else "<opaque>")
    if (inherits(value, c("gg", "ggplot", "grob", "gtable", "htmlwidget", "recordedplot"))) return(paste0("<", class(value)[[1L]], ">"))
    if (isS4(value)) return(list(.s4 = class(value)[[1L]]))
    reference <- is.list(value) && !is.data.frame(value) && is.character(value$entity_id) && length(value$entity_id) == 1L &&
        is.list(value$payload) && is.character(value$payload$sha256)
    if (reference) {
        if (identical(value$payload$codec, "parquet_duckdb") || is.null(store) || is.null(api))
            return(list(.ref = value$kind %||% "value", sha256 = value$payload$sha256))
        # Records are hydrated once per fingerprint run (memo by entity id); their reference registries
        # (`dependencies` and the like) accumulate every earlier state and are derived, not science.
        if (!is.null(memo) && !is.null(memo[[value$entity_id]])) return(memo[[value$entity_id]])
        hydrated <- api("hydrateSavedWorkflowEntity")(store, value)
        canonical <- list(.ref = value$kind %||% "value", value = kitCanonical(hydrated, store, api, depth + 1L, memo))
        if (!is.null(memo)) memo[[value$entity_id]] <- canonical
        return(canonical)
    }
    if (is.data.frame(value)) return(list(.df = digest::digest(kitMaskPaths(value), algo = "sha256"), dim = dim(value)))
    if (is.list(value)) {
        names_ <- names(value) %||% rep("", length(value))
        # `declaration` is the storage-format field KIT-002 added to checkpoint records (additive; records without
        # it read as before): a record's science is the same with or without it (KIT-011 correction 1).
        # `core_utilisation` is the recorded worker-cluster specification of a kernel run (NA for a declaration run in
        # a worker, the multidplyr spec for the native loop); `support_table`/`removal_ledger` inside an args section
        # are the native loop's copies of the support entries the checkpoint carries and compares separately
        # (KIT-008 correction 2, user decision 2026-10-12: execution bookkeeping, not science).
        value <- value[!(names_ %in% c("timestamp", "mtime", "dependencies", "product_dependencies", "da_dependencies",
            "owned_components", "component_dependencies", "declaration", "core_utilisation", "support_table", "removal_ledger") |
            grepl("timestamp$", names_))]
        return(lapply(value, kitCanonical, store = store, api = api, depth = depth + 1L, memo = memo))
    }
    kitMaskPaths(value)
}

#' Fingerprint one saved reference: tables by logical pages, other payloads by their canonical value.
#' @noRd
kitFingerprintReference <- function(store, reference, api, memo = NULL) {
    if (identical(reference$payload$codec, "parquet_duckdb")) return(kitFingerprintTable(store, reference, api))
    value <- api("hydrateSavedWorkflowEntity")(store, reference)
    list(kind = reference$kind %||% "value", product_type = if (is.list(value)) value$product_type else NULL,
        sha256 = digest::digest(kitCanonical(value, store, api, memo = memo), algo = "sha256"))
}

#' Fingerprint every component, support entry, product and the args summary of one state record.
#' The fingerprint carries no entity ids (those are run-specific); the population record does.
#' @param store Open entity store; `reference` the state's checkpoint descriptor.
#' @return Named list `components`, `support`, `products`, `summary_sha256`, `class`, `rows`.
#' @noRd
kitFingerprintState <- function(store, reference, api, memo = NULL) {
    record <- api("hydrateSavedWorkflowCheckpoint")(store, reference)
    spec <- record$data
    # The `initial` placeholder carries no data; it is recorded, not fingerprinted.
    if (is.null(spec)) return(list(class = NA_character_, rows = NA_integer_,
        components = list(), support = list(), products = list(), summary_sha256 = NA_character_, empty = TRUE))
    if (!api("isDIAComponentCheckpoint")(spec)) stop("state is not a component checkpoint")
    walk <- function(node, path = character()) {
        if (api("isSavedWorkflowReference")(node)) return(stats::setNames(list(kitFingerprintReference(store, node, api, memo)),
            paste(path, collapse = "/")))
        if (is.list(node) && !is.data.frame(node)) {
            out <- list()
            for (name in names(node)) out <- c(out, walk(node[[name]], c(path, name)))
            return(out)
        }
        list()
    }
    # The step description is display text (the declaration and the native loop word it differently); it is
    # recorded on the state, not part of the exactness summary (KIT-008 correction 2).
    summary <- list(class = spec$class, slots = kitCanonical(spec$slots), dimensions = spec$dimensions,
        config = kitCanonical(record$config))
    list(class = spec$class, rows = spec$dimensions$rows,
        components = walk(spec$components), support = walk(spec$support), products = walk(spec$products),
        summary_sha256 = digest::digest(summary, algo = "sha256"))
}

#' Fingerprint publication files under a project directory (store, logs and scratch excluded).
#' PDF `CreationDate`/`ModDate` lines are stripped before hashing.
#' @return Data frame `path, bytes, sha256, method`.
#' @noRd
kitFingerprintFiles <- function(project, previous = NULL) {
    files <- list.files(project, recursive = TRUE, full.names = TRUE, all.files = TRUE)
    files <- files[!grepl("/[.]workflow-entities/|/[.]dia-(native|gui-adoption)-[^/]*/|[.]log$|/checkpoint-io[.]log$", files)]
    files <- files[!dir.exists(files)]
    relative <- substring(files, nchar(project) + 2L)
    info <- file.info(files, extra_cols = FALSE)
    mtime <- as.numeric(info$mtime)
    bytes <- as.numeric(info$size)
    # Reuse the previous pass's hash when path, size and mtime are unchanged (the verifier recomputes everything).
    reuse <- if (is.null(previous)) rep(FALSE, length(files)) else
        paste(relative, bytes, mtime) %in% paste(previous$path, previous$bytes, previous$mtime)
    hash <- function(path) {
        lower <- tolower(path)
        if (grepl("[.]pdf$", lower)) {
            raw <- readBin(path, "raw", file.size(path))
            text <- rawToChar(raw[raw != as.raw(0L)], multiple = FALSE)
            text <- gsub("/(CreationDate|ModDate) *\\([^)]*\\)", "", text, useBytes = TRUE)
            return(c(digest::digest(text, algo = "sha256", serialize = FALSE), "sha256_pdf_dates_stripped"))
        }
        # Serialized R objects, widgets and workbooks carry run-specific ids or times: hash their content.
        if (grepl("[.]rds$", lower)) {
            value <- tryCatch(readRDS(path), error = function(error) NULL)
            if (!is.null(value)) return(c(digest::digest(kitCanonical(value), algo = "sha256"), "sha256_rds_canonical"))
        }
        if (grepl("[.]html$", lower)) {
            text <- gsub("htmlwidget-[A-Za-z0-9]+", "htmlwidget-X", readLines(path, warn = FALSE))
            return(c(digest::digest(text, algo = "sha256", serialize = FALSE), "sha256_html_widget_ids_masked"))
        }
        if (grepl("[.]xlsx$", lower) && requireNamespace("openxlsx", quietly = TRUE)) {
            sheets <- lapply(openxlsx::getSheetNames(path), function(sheet) openxlsx::read.xlsx(path, sheet = sheet))
            return(c(digest::digest(sheets, algo = "sha256"), "sha256_xlsx_sheets"))
        }
        c(digest::digest(file = path, algo = "sha256"), "sha256")
    }
    hashed <- matrix(NA_character_, nrow = 2L, ncol = length(files))
    if (any(reuse)) {
        at <- match(relative[reuse], previous$path)
        hashed[1L, reuse] <- previous$sha256[at]
        hashed[2L, reuse] <- previous$method[at]
    }
    if (any(!reuse)) hashed[, !reuse] <- vapply(files[!reuse], hash, character(2))
    data.frame(path = relative, bytes = bytes, mtime = mtime, sha256 = hashed[1L, ],
        method = hashed[2L, ], stringsAsFactors = FALSE, row.names = NULL)
}

#' Locate the DIA state store beneath a project root (native session from import, or a GUI adoption store).
#' @noRd
kitFindStore <- function(project) {
    heads <- list.files(project, "^dia-state-head[.]rds$", recursive = TRUE, full.names = TRUE, all.files = TRUE)
    heads <- heads[grepl("/[.]dia-(native|gui-adoption)-[^/]*/[.]workflow-entities/dia-state-head[.]rds$", heads)]
    if (!length(heads)) return(NULL)
    if (length(heads) != 1L) stop("DIA state store is ambiguous")
    dirname(heads)
}

#' Open an existing store read-only through the installed namespace.
#' @noRd
kitOpenStore <- function(store_dir, api) {
    identity <- api("readSavedStoreIdentity")(file.path(store_dir, "store.json"))
    api("newWorkflowEntityStore")(dirname(store_dir), identity$project_id, identity$workflow_id, max_cache_bytes = 0)
}

#' Read the store head and its state references without hydrating any assay.
#' @return `head_sha256`, `current_state`, `history`, `states` (name -> entity id) and `rows` per state.
#' @noRd
kitReadHead <- function(store_dir, api) {
    head_path <- file.path(store_dir, "dia-state-head.rds")
    head <- readRDS(head_path)
    store <- kitOpenStore(store_dir, api)
    on.exit(api("closeWorkflowEntityStore")(store), add = TRUE)
    snapshot <- api("hydrateSavedWorkflowCheckpoint")(store, head)
    list(head_sha256 = digest::digest(file = head_path, algo = "sha256"),
        current_state = snapshot$current_state, history = as.character(unlist(snapshot$state_history)),
        states = lapply(snapshot$states, `[[`, "entity_id"),
        rows = lapply(snapshot$sizes, `[[`, "rows"), references = snapshot$states,
        products = lapply(snapshot$products, `[[`, "entity_id"))
}

#' Find the checkpoint entity whose descriptor is byte-identical to a head file with the given hash
#' (a store's `dia-state-head.rds` is the descriptor of the checkpoint entity it points at, so the head
#' as it stood after any step survives in `entities/<id>/descriptor.rds` of the immutable store).
#' @return The entity directory name, or NULL when no descriptor matches.
#' @noRd
kitFindHeadEntity <- function(store_dir, head_sha256) {
    descriptors <- list.files(file.path(store_dir, "entities"), pattern = "^descriptor\\.rds$", recursive = TRUE, full.names = TRUE)
    for (path in descriptors) if (identical(digest::digest(file = path, algo = "sha256"), head_sha256)) return(basename(dirname(path)))
    NULL
}

#' Validity of a bound population (KIT-005 audit PH-01). A `core` population is a frozen store whose
#' live head is the bound head. A journey population is a step inside a live store whose head has moved
#' on: its bound head is the descriptor of the checkpoint entity `head_entity`, which must still hash to
#' `head_sha256`.
#' @return NULL when valid, else the reason.
#' @noRd
kitPopulationInvalid <- function(input) {
    if (!is.character(input$store) || !dir.exists(input$store)) return("store_missing")
    head <- if (is.null(input$head_entity)) file.path(input$store, "dia-state-head.rds") else
        file.path(input$store, "entities", input$head_entity, "descriptor.rds")
    if (!file.exists(head)) return("head_missing")
    if (!identical(digest::digest(file = head, algo = "sha256"), input$head_sha256)) return("head_hash_mismatch")
    NULL
}

#' The store an operation member works from. A `core` population is used in place (the member copies it).
#' A journey population is materialised under `destination`: the immutable entities and the store identity
#' are reflinked (a plain copy where the filesystem cannot share extents) and the bound checkpoint descriptor
#' becomes the head, so the member sees the store exactly as it stood after that step and never writes into
#' the journey member's evidence store. When `api` is given, the head's current state is checked against the
#' bound state and entity id.
#' @param destination A fresh directory; the snapshot store is its `.workflow-entities` (a store is always
#'   opened from its parent directory, so the store directory keeps the production name).
#' @return The store directory to bind as `source_store`.
#' @noRd
kitPopulationStore <- function(input, destination = NULL, api = NULL) {
    reason <- kitPopulationInvalid(input)
    if (!is.null(reason)) stop("bound population is not usable: ", reason)
    if (is.null(input$head_entity)) return(input$store)
    stopifnot(is.character(destination), length(destination) == 1L, !dir.exists(destination))
    snapshot <- file.path(destination, ".workflow-entities")
    dir.create(snapshot, recursive = TRUE)
    for (name in c("entities", "store.json", "checkpoint-io.log")) {
        source <- file.path(input$store, name)
        if (!file.exists(source)) next
        status <- system2("cp", c("-a", "--reflink=auto", shQuote(source), shQuote(file.path(snapshot, name))))
        if (!identical(status, 0L)) stop("population snapshot copy failed: ", name)
    }
    destination <- snapshot
    head <- file.path(destination, "dia-state-head.rds")
    stopifnot(file.copy(file.path(input$store, "entities", input$head_entity, "descriptor.rds"), head),
        identical(digest::digest(file = head, algo = "sha256"), input$head_sha256))
    if (!is.null(api)) {
        current <- kitReadHead(destination, api)
        stopifnot(identical(current$current_state, input$state), identical(current$states[[input$state]], input$entity_id))
    }
    destination
}

#' The peptide_data component entity id of the store's current state: an operation member whose
#' worker failed must leave it equal to the bound input (the head unchanged).
#' @noRd
kitCurrentInputComponentId <- function(store_dir, api) {
    head <- kitReadHead(store_dir, api)
    store <- kitOpenStore(store_dir, api)
    on.exit(api("closeWorkflowEntityStore")(store), add = TRUE)
    record <- api("hydrateSavedWorkflowCheckpoint")(store, head$references[[head$current_state]])
    record$data$components$peptide_data$entity_id
}

#' Fingerprint every state and product the head references.
#' @return Named list of `kitFingerprintState()` results plus `products`.
#' @noRd
kitFingerprintStore <- function(store_dir, api, states = NULL) {
    head <- kitReadHead(store_dir, api)
    store <- kitOpenStore(store_dir, api)
    on.exit(api("closeWorkflowEntityStore")(store), add = TRUE)
    names <- if (is.null(states)) names(head$references) else intersect(states, names(head$references))
    memo <- new.env(parent = emptyenv())
    result <- lapply(stats::setNames(names, names), \(name) kitFingerprintState(store, head$references[[name]], api, memo))
    snapshot <- api("hydrateSavedWorkflowCheckpoint")(store, readRDS(file.path(store_dir, "dia-state-head.rds")))
    products <- lapply(snapshot$products, \(reference) kitFingerprintReference(store, reference, api, memo))
    list(head_sha256 = head$head_sha256, states = result, products = products)
}

#' Does an in-run state fingerprint match the one recomputed from the immutable store? The in-run
#' record also carries the products that appeared after the step (snapshot products), which a state
#' record does not reference; those are compared only where the recomputed record has them.
#' @noRd
kitStateFingerprintMatches <- function(expected, recomputed) {
    if (is.null(expected) || is.null(recomputed)) return(FALSE)
    for (field in c("class", "rows", "components", "support", "summary_sha256"))
        if (!identical(expected[[field]], recomputed[[field]])) return(FALSE)
    for (name in names(recomputed$products))
        if (!identical(expected$products[[name]], recomputed$products[[name]])) return(FALSE)
    TRUE
}

#' Compare a member's publication files with another member's: content hashes for every file except
#' PDFs, whose compressed metadata (ids, dates) differs between identical runs; PDFs compare by path and size.
#' @return Character vector of differing paths.
#' @noRd
kitFileMismatches <- function(reference, candidate) {
    canonical <- function(files) {
        files <- files[!grepl("/[.]workflow-export/", files$path), , drop = FALSE]
        files$path <- gsub("[0-9]{8}_[0-9]{6}", "<time>", files$path)
        files$path <- gsub("[0-9a-f]{13,}", "<id>", files$path)
        files
    }
    # Compared by path only: PDFs (compressed metadata differs between identical runs) and the session
    # export's bookkeeping files, which carry the export time and the session id (KIT-005 correction 7).
    by_path_only <- function(path) grepl("[.]pdf$|[.]msdia[.]tar$|(^|/)native_export_complete_|(^|/)native-(io-)?preflight-|(^|/)filtered_session_(latest[.]json|summary[.]txt)$",
        path, ignore.case = TRUE)
    key <- function(files) paste(files$path, ifelse(by_path_only(files$path), "", files$sha256))
    reference <- canonical(reference)
    candidate <- canonical(candidate)
    a <- key(reference)
    b <- key(candidate)
    unique(c(reference$path[!a %in% b], candidate$path[!b %in% a]))
}

#' Verify one KIT member outside its measured slice: binding, build, steps, head and fingerprints.
#' Journey members: every in-run population fingerprint is recomputed from the immutable store and the
#' publication files re-hashed; a second repetition must also reproduce the first repetition's
#' fingerprints and files (the oracle reproducibility check); the production-default member has no
#' store and verifies on completion and phase measures. Operation members: the installed action
#' journal, completion and worker result, the P14 page oracle (precursor rollup) or the sample-filter
#' oracle, and the exact-sum identity; a failed member must leave the head on its input.
#' @param run_dir Verifier output directory. `repo` the frozen snapshot holding the bindings.
#' @return Verification receipt saved as `verification.rds` beside the member.
#' @noRd
verifyKitMember <- function(run_dir, repo, selection, api) {
    target <- kitSelectorStage(selection)
    bindings <- kitBindings(file.path(repo, "docs/qa/workflow-execution/kit/bindings.json"))
    cases <- kitCases(bindings, target$stage)
    case <- Filter(\(x) x$id == target$id, cases)[[1L]]
    owned <- file.path(bindings$root, target$stage, case$id)
    path <- file.path(owned, "verification.rds")
    stopifnot(!file.exists(path))
    receipt <- list(status = "starting", run_dir = run_dir, member_id = case$id, selector = selection, kind = case$kind)
    on.exit({
        if (identical(receipt$status, "starting")) receipt$status <- "failed"
        saveRDS(receipt, path)
    }, add = TRUE)
    operation_path <- file.path(owned, "operation-receipt.rds")
    receipt$operation_sha256 <- digest::digest(file = operation_path, algo = "sha256")
    operation <- readRDS(operation_path)
    stopifnot(identical(operation$binding, case$binding), identical(operation$config$build, case$expected_build))
    if (case$kind == "operation") {
        # Assigned, not returned directly: the exit handler saves this variable (KIT-005 correction 2).
        receipt <- verifyKitOperationMember(receipt, owned, case, bindings, operation, api)
        return(receipt)
    }
    journey <- readRDS(file.path(owned, "journey-result.rds"))
    stopifnot(identical(journey$corpus, case$corpus))
    corpus <- bindings$corpora[[case$corpus]]
    for (name in names(corpus$files)) stopifnot(identical(digest::digest(file = corpus$files[[name]]$path,
        algo = "sha256"), corpus$files[[name]]$sha256))
    status <- operation$phase_status
    steps <- kitStepIds()
    if (case$expect_failure) {
        fault <- bindings$lifecycle$fault
        stopifnot(identical(operation$status, "journey_failed"), identical(operation$guard$reason, "injected_guard_stop"),
            identical(status[[fault$failed_step]], "failed"), identical(status[[fault$after_step]], "complete"))
        store_dir <- journey$store
        stopifnot(!is.null(store_dir), dir.exists(store_dir))
        head <- kitReadHead(store_dir, api)
        last_good <- journey$populations[[fault$after_step]]
        stopifnot(identical(head$head_sha256, last_good$head_sha256), identical(head$current_state, last_good$state))
        receipt$failure <- list(failed_step = fault$failed_step, head_sha256 = head$head_sha256,
            current_state = head$current_state, guard = operation$guard)
        receipt$status <- "complete"
        return(receipt)
    }
    if (!isTRUE(case$record) && !identical(operation$status, "complete")) {
        # Production default ended at a protocol guard (KIT-005 acceptance as amended 2026-10-11): the
        # completed steps' peaks are the record; no store, no populations, no file fingerprints.
        outcome <- kitDefaultOutcome(operation)
        stopifnot(!is.null(outcome), is.null(journey$store))
        completed <- names(status)[unlist(status) == "complete"]
        stopifnot(length(completed) >= 1L, all(is.finite(unlist(operation$phase_peak_anon_bytes[completed]))))
        receipt$guard <- list(outcome = outcome, completed_steps = completed,
            failed_step = names(status)[unlist(status) == "failed"])
        receipt$populations <- NULL
        receipt$status <- "complete"
        return(receipt)
    }
    stopifnot(identical(operation$status, "complete"), all(steps %in% names(status)), all(unlist(status[steps]) == "complete"))
    files <- kitFingerprintFiles(journey$project)
    stopifnot(identical(files[c("path", "bytes", "sha256", "method")], journey$files[c("path", "bytes", "sha256", "method")]))
    receipt$files <- nrow(files)
    if (!isTRUE(case$record) && is.null(journey$store)) {
        # Production default on an in-memory build: no store, no populations; completion and measures only.
        stopifnot(is.numeric(operation$resident_end_bytes), is.finite(operation$resident_end_bytes),
            all(is.finite(unlist(operation$phase_peak_anon_bytes[steps]))))
        receipt$populations <- NULL
        receipt$status <- "complete"
        return(receipt)
    }
    if (!isTRUE(case$record) && identical(case$corpus, "full")) {
        # Production default on a build whose default is the on-disk path (KIT-011): the in-run populations are
        # recomputed from the immutable store and every state must equal the KIT-005 per-step oracle. The
        # oracles are the 300-sample corpus's; a screen member without a record (KIT-008's gui member) takes
        # the self-consistency branch below and is reproduced against its reference screen member by its decision.
        store_dir <- journey$store
        stopifnot(dir.exists(store_dir), isTRUE(journey$populations_complete))
        recomputed <- kitFingerprintStore(store_dir, api)
        stopifnot(identical(recomputed$head_sha256, journey$final_head_sha256))
        # Integrity: the states the driver recorded in-run are the ones in the immutable store (class, rows,
        # components, support, summary). Products are judged on the recomputed fingerprints, so a canonical-form
        # correction between the run and this verification (KIT-011 correction 1) is recorded, not fatal.
        recomputed_populations <- journey$populations
        receipt$canonical_form_differences <- character()
        for (step in names(journey$populations)) {
            population <- journey$populations[[step]]
            if (is.null(population$state)) next
            fresh <- recomputed$states[[population$state]]
            for (field in c("class", "rows", "components", "support", "summary_sha256"))
                stopifnot(identical(population$fingerprint[[field]], fresh[[field]]))
            if (!kitStateFingerprintMatches(population$fingerprint, fresh))
                receipt$canonical_form_differences <- c(receipt$canonical_form_differences, step)
            # The in-run fingerprint also carries the per-step assets the driver attached from the head; the
            # recomputed values overlay the entries the state record itself holds.
            merged <- population$fingerprint
            for (name in names(fresh$products)) merged$products[[name]] <- fresh$products[[name]]
            recomputed_populations[[step]]$fingerprint <- merged
        }
        # Cross-build: the operation identity products are reported, every state, support entry, summary and
        # other product is judged (KIT-011 screen check: identical states on kit_001 and kit_011).
        receipt$oracle_mismatches <- kitReproduceMismatches(kitOracleJourneyPopulations(bindings), recomputed_populations,
            ignore_products = kitOperationIdentityProducts)
        receipt$operation_identity_products <- unlist(lapply(journey$populations, function(p)
            grep(kitOperationIdentityProducts, names(p$fingerprint$products), value = TRUE)), use.names = FALSE)
        stopifnot(!length(receipt$oracle_mismatches))
        receipt$populations <- lapply(journey$populations, \(p) list(state = p$state, entity_id = p$entity_id, rows = p$rows))
        receipt$oracle <- "KIT-005 journey_<step> fingerprints"
        receipt$status <- "complete"
        return(receipt)
    }
    store_dir <- journey$store
    stopifnot(!is.null(store_dir), dir.exists(store_dir))
    recomputed <- kitFingerprintStore(store_dir, api)
    stopifnot(identical(recomputed$head_sha256, journey$final_head_sha256))
    mismatches <- character()
    for (step in names(journey$populations)) {
        population <- journey$populations[[step]]
        if (is.null(population$state)) next
        if (!kitStateFingerprintMatches(population$fingerprint, recomputed$states[[population$state]])) mismatches <- c(mismatches, step)
    }
    stopifnot(!length(mismatches), isTRUE(journey$populations_complete))
    receipt$populations <- lapply(journey$populations, \(p) list(state = p$state, entity_id = p$entity_id,
        rows = p$rows, fingerprints = length(p$fingerprint$components) + length(p$fingerprint$support) +
            length(p$fingerprint$products) + 1L))
    receipt$products <- length(recomputed$products)
    if (case$rep == 2L) {
        # Reproducibility: the first repetition of the same role is the oracle (KIT-005 item 5).
        first <- readRDS(file.path(bindings$root, target$stage, paste0(case$role, "-1"), "journey-result.rds"))
        receipt$reproduce_mismatches <- kitReproduceMismatches(first$populations, journey$populations)
        receipt$file_mismatches <- kitFileMismatches(first$files, journey$files)
        receipt$reproduced <- !length(receipt$reproduce_mismatches) && !length(receipt$file_mismatches)
        stopifnot(isTRUE(receipt$reproduced))
    }
    receipt$status <- "complete"
    receipt
}

#' Operation-member verification (CORE's checks carried into the kit; no `population == "p14"` assertion).
#' @noRd
verifyKitOperationMember <- function(receipt, owned, case, bindings, operation, api) {
    member <- readRDS(file.path(owned, "member-result.rds"))
    expected_status <- if (case$expect_failure) "worker_failed" else "complete"
    # The bound population is still valid and the store the member worked from carried the bound head
    # (for a journey population that is the per-step snapshot the coordinator materialised, KIT-005 audit PH-01).
    source_store <- operation$member$source_store %||% case$input$store
    stopifnot(identical(operation$status, expected_status), identical(member$status, expected_status),
        identical(member$member$input_entity_id, case$input$entity_id), is.null(kitPopulationInvalid(case$input)),
        identical(digest::digest(file = file.path(source_store, "dia-state-head.rds"), algo = "sha256"), case$input$head_sha256))
    if (!is.null(case$input$state)) {
        source_head <- kitReadHead(source_store, api)
        stopifnot(identical(source_head$current_state, case$input$state),
            identical(source_head$states[[case$input$state]], case$input$entity_id))
    }
    project <- file.path(owned, "project")
    owner <- api("readSavedStoreIdentity")(file.path(project, ".workflow-entities/store.json"))
    manager <- api("newDIAStateAdapter")(api("newWorkflowEntityStore")(project, owner$project_id, owner$workflow_id), reopen = TRUE)
    on.exit(api("closeDIAStateAdapter")(manager), add = TRUE, after = FALSE)
    checkpoint <- api("readDIAStateMetadata")(manager)$data
    request <- file.path(project, ".workflow-entities/.workflow-export", operation$execution$request_id)
    journal <- jsonlite::read_json(file.path(request, "journal.json"))
    completion <- readRDS(file.path(request, "completion.rds"))
    worker <- readRDS(file.path(request, "worker-result.rds"))
    stopifnot(!identical(journal$status, "conflict"), !identical(completion$status, "guard_payload_invalid"),
        !isTRUE(completion$monitor_lost_after_result), identical(completion$status, expected_status),
        isTRUE(completion$cleanup$worker_stopped), isTRUE(completion$cleanup$monitor_stopped))
    input_component <- checkpoint$components$peptide_data
    if (case$expect_failure) {
        fault <- bindings$operation_lifecycle$fault
        stopifnot(identical(journal$status, "worker_failed"), identical(worker$status, "failed"),
            grepl(fault$error_pattern, worker$error, fixed = TRUE),
            identical(input_component$entity_id, case$input$entity_id))
        receipt$failure <- list(worker_error = worker$error, journal_status = journal$status,
            head_entity = input_component$entity_id, completion_status = completion$status)
        receipt$status <- "complete"
        return(receipt)
    }
    stopifnot(identical(journal$status, "complete"))
    execution <- worker$result$execution
    store <- manager$.context$store
    hashes <- kitFingerprintTable(store, input_component, api)
    if (case$operation == "precursor_rollup") {
        screen <- identical(case$input$rows, bindings$populations$screen_p14$rows)
        oracle <- bindings$oracles$p14
        # `screen_pages`/`full_pages` are the page counts; the hashes are `*_page_hashes` (KIT-005 correction 1).
        expected <- unlist(oracle[[if (screen) "screen_page_hashes" else "full_page_hashes"]], use.names = FALSE)
        stopifnot(length(expected) == oracle[[if (screen) "screen_pages" else "full_pages"]])
        pages <- kitFingerprintTablePages(store, input_component, api)
        stopifnot(identical(unname(pages), unname(expected)), identical(execution$exact_sum_sha256, case$exact_sum),
            identical(execution$query_owner, "fastduckdb_c"), identical(execution$profile$name, "P14"),
            all(kitRequiredOperationPhases %in% names(operation$phase_status)))
        receipt$pages <- pages
    } else if (case$operation == "minimum_peptides_per_sample") {
        # The sample-filter oracle is the native member's verification record bound by CORE
        # (`oracles.filter_primary`): page hashes of the filtered table, design, summary and support.
        oracle <- bindings$oracles[[paste0("filter_", sub("^full_", "", case$population))]]
        stopifnot(!is.null(oracle), file.exists(oracle$path),
            identical(digest::digest(file = oracle$path, algo = "sha256"), oracle$sha256))
        expected <- readRDS(oracle$path)$verification
        actual <- kitFingerprintFilterResult(manager, checkpoint, api)
        stopifnot(identical(actual$page_hashes, expected$page_hashes), identical(actual$rows, expected$rows),
            identical(actual$samples, expected$samples), identical(actual$columns, expected$columns),
            identical(actual$design_sha256, expected$design_sha256), identical(actual$summary_sha256, expected$summary_sha256),
            identical(actual$support_sha256, expected$support_sha256), identical(execution$output_mode, "physical"))
        # The control's fixed-three policy records `count_batch_sizes` only; the adaptive admission records
        # `batches` as well (KIT-006 rerun 1: the verifier refused the control for the missing field).
        batches <- execution$admission$batches %||% length(execution$admission$count_batch_sizes)
        stopifnot(is.numeric(batches), length(batches) == 1L, batches >= 3L)
        receipt$filter <- actual[c("rows", "samples", "design_sha256", "summary_sha256", "support_sha256")]
        receipt$batches <- batches
    } else {
        oracle <- bindings$oracles[[paste0("operation_", case$population)]]
        stopifnot(!is.null(oracle), identical(hashes$sha256, oracle$sha256))
    }
    receipt$output <- hashes
    receipt$execution <- execution
    receipt$status <- "complete"
    receipt
}

#' Fingerprint a sample-filter result exactly as the CORE repair recorded its native oracle
#' (`verifyPerf106Result`): page digests of the raw 10,000-row pages, then the design, the filter
#' summary and the support tables.
#' @noRd
kitFingerprintFilterResult <- function(manager, checkpoint, api) {
    store <- manager$.context$store
    reference <- checkpoint$components$peptide_data
    scope <- api("newWorkflowTableReadScope")(store, reference, reader = "fastduckdb")
    on.exit(api("releaseWorkflowTableReadScope")(scope), add = TRUE)
    hashes <- character()
    for (offset in seq.int(0L, max(0L, scope$metadata$rows - 1L), by = 10000L)) {
        page <- api("readSavedWorkflowTable")(store, reference, scope$metadata$columns, offset = offset,
            limit = min(10000L, scope$metadata$rows - offset), scope = scope)$data
        hashes[[length(hashes) + 1L]] <- digest::digest(page, algo = "sha256")
    }
    design <- api("materializeDIATable")(manager$.context, checkpoint$components$design_matrix)
    support <- lapply(checkpoint$support$filterMinNumPeptidesPerSample,
        function(reference) api("materializeDIATable")(manager$.context, reference))
    list(rows = scope$metadata$rows, samples = nrow(design), columns = scope$metadata$columns, page_hashes = hashes,
        design_sha256 = digest::digest(design, algo = "sha256"),
        summary_sha256 = digest::digest(checkpoint$slots$args$filterMinNumPeptidesPerSample$filter_summary, algo = "sha256"),
        support_sha256 = digest::digest(support, algo = "sha256"))
}

#' The P14 page hashes of a table (every column, 10,000 rows per page), as the bound oracle records them.
#' @noRd
kitFingerprintTablePages <- function(store, reference, api) {
    scope <- api("newWorkflowTableReadScope")(store, reference, reader = "fastduckdb")
    on.exit(api("releaseWorkflowTableReadScope")(scope), add = TRUE)
    kitPageHashes(store, reference, api, scope)
}

`%||%` <- function(x, y) if (is.null(x)) y else x
