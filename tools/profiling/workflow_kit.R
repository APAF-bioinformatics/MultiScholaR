# KIT series harness: bindings, fixed member cases, journey step map, phase measures,
# frozen decisions, store fingerprints and the independent verifier (KIT-001).
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
        step("import", "Import", "raw_data_s4", "store", "mod_prot_import", "prot-import-process"),
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
    if (grepl("^kit-lifecycle-", stage)) return(300L)
    if (stage == "kit-journey") {
        seconds <- bindings$journey$ceiling_seconds
        if (!is.numeric(seconds) || length(seconds) != 1L || !is.finite(seconds) || seconds <= 0 || seconds > 5400)
            stop("the journey ceiling is not fixed in the bindings (at most 90 min)")
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
#' @param stage `kit-lifecycle-complete`, `kit-lifecycle-fail` or `kit-journey`.
#' @return Fixed cases in launch order.
#' @noRd
kitCases <- function(bindings, stage) {
    stopifnot(grepl("^kit-(lifecycle-(complete|fail)|journey|0[0-9]{2})$", stage))
    case <- function(id, role, rep, build_key, corpus_key, fault, selector) {
        build <- bindings$builds[[build_key]]
        if (is.null(build)) stop("build ", build_key, " is not bound")
        corpus <- bindings$corpora[[corpus_key]]
        if (is.null(corpus)) stop("corpus ", corpus_key, " is not bound")
        list(id = id, role = role, control = role == "control", rep = rep, fault = fault,
            corpus = corpus_key, library = build$library, expected_build = build$build,
            exact_sum = build$exact_sum, expect_failure = fault != "none",
            binding = digest::digest(list(corpus = corpus$files, build = build, stage = stage, rep = rep,
                journey = bindings$journey, protocol = bindings$protocol_sha256,
                step_map = bindings$step_map_sha256), algo = "sha256"),
            selector = selector)
    }
    if (startsWith(stage, "kit-lifecycle-")) {
        kind <- sub("kit-lifecycle-", "", stage, fixed = TRUE)
        return(list(case("member", "candidate", 1L, kitControl(bindings), "screen",
            fault = if (kind == "fail") bindings$lifecycle$fault$name else "none",
            selector = paste0("kit:lifecycle:", kind))))
    }
    if (stage == "kit-journey") {
        spec <- bindings$journey$pair
        if (is.null(spec)) stop("the journey pair is defined by KIT-005")
        return(list(
            case("control-1", "control", 1L, spec$control, "full", "none", "kit:journey:control:1"),
            case("candidate-1", "candidate", 1L, spec$candidate, "full", "none", "kit:journey:candidate:1")))
    }
    stop("operation members are defined by KIT-005")
}

#' Map a verifier selector back to its stage and case id.
#' @noRd
kitSelectorStage <- function(selection) {
    stopifnot(grepl("^kit:(journey:(control|candidate):[12]|lifecycle:(complete|fail)|0[0-9]{2}:(control|candidate):[12])$", selection))
    parts <- strsplit(sub("^kit:", "", selection), ":", fixed = TRUE)[[1L]]
    if (parts[[1L]] == "lifecycle") return(list(stage = paste0("kit-lifecycle-", parts[[2L]]), id = "member"))
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

#' Measurement validity: coordinator sampling gap at most 250 ms, zero swap, slice stopped.
#' @noRd
kitMeasurementValid <- function(op) {
    gap <- op$max_sample_interval
    is.numeric(gap) && length(gap) == 1L && is.finite(gap) && gap <= 0.25 &&
        is.numeric(op$peak_bytes) && length(op$peak_bytes) == 1L && is.finite(op$peak_bytes) &&
        identical(op$swap_bytes, 0) && isTRUE(op$cleanup$slice_stopped)
}

#' Check one member against the frozen eligibility rules.
#' @return Name of the first failing rule, or NULL when eligible.
#' @noRd
kitIneligible <- function(item, seconds) {
    op <- item$operation
    rules <- list(
        verifier_accepted = isTRUE(item$verifier_completion$accepted),
        status_complete = identical(op$status, "complete"),
        measurement_valid = kitMeasurementValid(op),
        zero_swap = identical(op$swap_bytes, 0),
        peak_below_early_stop = is.numeric(op$peak_bytes) && length(op$peak_bytes) == 1L &&
            is.finite(op$peak_bytes) && op$peak_bytes > 0 && op$peak_bytes < 15.5 * 1024^3,
        elapsed_within_ceiling = is.numeric(op$elapsed_seconds) && length(op$elapsed_seconds) == 1L &&
            is.finite(op$elapsed_seconds) && op$elapsed_seconds > 0 && op$elapsed_seconds <= seconds)
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
        result$failing_rule <- kitIneligible(item, kitMemberSeconds("kit-lifecycle-complete"))
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

#' Mask run-specific absolute paths so args summaries compare across members.
#' @noRd
kitMaskPaths <- function(value) {
    if (is.character(value)) {
        masked <- value
        at <- !is.na(masked) & startsWith(masked, "/")
        masked[at] <- paste0("<path>/", basename(masked[at]))
        return(masked)
    }
    if (is.list(value)) return(lapply(value, kitMaskPaths))
    value
}

#' Hash a table component by logical pages (every column, every row, 10,000 rows per page).
#' @noRd
kitFingerprintTable <- function(store, reference, api) {
    scope <- api("newWorkflowTableReadScope")(store, reference, reader = "fastduckdb")
    on.exit(api("releaseWorkflowTableReadScope")(scope), add = TRUE)
    rows <- scope$metadata$rows
    hashes <- character()
    if (rows > 0L) for (offset in seq.int(0L, rows - 1L, by = 10000L)) {
        data <- api("readSavedWorkflowTable")(store, reference, scope$metadata$columns,
            offset = offset, limit = min(10000L, rows - offset), scope = scope)$data
        attr(data, "dia_source_order") <- NULL
        attr(data, "row.names") <- .set_row_names(nrow(data))
        hashes[[length(hashes) + 1L]] <- digest::digest(data, algo = "sha256")
    }
    list(kind = "table", rows = rows, columns = scope$metadata$columns,
        sha256 = digest::digest(paste(hashes, collapse = "\n"), algo = "sha256", serialize = FALSE))
}

#' Fingerprint one saved reference: tables by logical pages, other payloads by their hydrated value.
#' @noRd
kitFingerprintReference <- function(store, reference, api) {
    if (identical(reference$payload$codec, "parquet_duckdb")) return(kitFingerprintTable(store, reference, api))
    value <- api("hydrateSavedWorkflowEntity")(store, reference)
    list(kind = reference$kind %||% "value", sha256 = digest::digest(kitMaskPaths(value), algo = "sha256"))
}

#' Fingerprint every component, support entry, product and the args summary of one state record.
#' @param store Open entity store; `reference` the state's checkpoint descriptor.
#' @return Named list `components`, `support`, `products`, `summary_sha256`, `class`, `rows`.
#' @noRd
kitFingerprintState <- function(store, reference, api) {
    record <- api("hydrateSavedWorkflowCheckpoint")(store, reference)
    spec <- record$data
    if (!api("isDIAComponentCheckpoint")(spec)) stop("state is not a component checkpoint")
    walk <- function(node, path = character()) {
        if (api("isSavedWorkflowReference")(node)) return(stats::setNames(list(kitFingerprintReference(store, node, api)),
            paste(path, collapse = "/")))
        if (is.list(node) && !is.data.frame(node)) {
            out <- list()
            for (name in names(node)) out <- c(out, walk(node[[name]], c(path, name)))
            return(out)
        }
        list()
    }
    summary <- list(class = spec$class, slots = kitMaskPaths(spec$slots), dimensions = spec$dimensions,
        config = kitMaskPaths(record$config), description = record$description)
    list(class = spec$class, rows = spec$dimensions$rows, entity_id = reference$entity_id,
        components = walk(spec$components), support = walk(spec$support), products = walk(spec$products),
        summary_sha256 = digest::digest(summary, algo = "sha256"))
}

#' Fingerprint publication files under a project directory (store, logs and scratch excluded).
#' PDF `CreationDate`/`ModDate` lines are stripped before hashing.
#' @return Data frame `path, bytes, sha256, method`.
#' @noRd
kitFingerprintFiles <- function(project, previous = NULL) {
    files <- list.files(project, recursive = TRUE, full.names = TRUE, all.files = TRUE)
    files <- files[!grepl("/[.]workflow-entities/|/[.]dia-gui-adoption-[^/]*/|[.]log$|/checkpoint-io[.]log$", files)]
    files <- files[!dir.exists(files)]
    relative <- substring(files, nchar(project) + 2L)
    info <- file.info(files, extra_cols = FALSE)
    mtime <- as.numeric(info$mtime)
    bytes <- as.numeric(info$size)
    # Reuse the previous pass's hash when path, size and mtime are unchanged (the verifier recomputes everything).
    reuse <- if (is.null(previous)) rep(FALSE, length(files)) else
        paste(relative, bytes, mtime) %in% paste(previous$path, previous$bytes, previous$mtime)
    hash <- function(path) {
        if (grepl("[.]pdf$", path, ignore.case = TRUE)) {
            raw <- readBin(path, "raw", file.size(path))
            text <- rawToChar(raw[raw != as.raw(0L)], multiple = FALSE)
            text <- gsub("/(CreationDate|ModDate) *\\([^)]*\\)", "", text, useBytes = TRUE)
            return(c(digest::digest(text, algo = "sha256", serialize = FALSE), "sha256_pdf_dates_stripped"))
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

#' Locate the adopted DIA store beneath a project root (NULL before adoption).
#' @noRd
kitFindStore <- function(project) {
    heads <- list.files(project, "^dia-state-head[.]rds$", recursive = TRUE, full.names = TRUE, all.files = TRUE)
    heads <- heads[grepl("/[.]dia-gui-adoption-[^/]*/[.]workflow-entities/dia-state-head[.]rds$", heads)]
    if (!length(heads)) return(NULL)
    if (length(heads) != 1L) stop("adopted GUI state store is ambiguous")
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

#' Fingerprint every state and product the head references.
#' @return Named list of `kitFingerprintState()` results plus `products`.
#' @noRd
kitFingerprintStore <- function(store_dir, api, states = NULL) {
    head <- kitReadHead(store_dir, api)
    store <- kitOpenStore(store_dir, api)
    on.exit(api("closeWorkflowEntityStore")(store), add = TRUE)
    names <- if (is.null(states)) names(head$references) else intersect(states, names(head$references))
    result <- lapply(stats::setNames(names, names), \(name) kitFingerprintState(store, head$references[[name]], api))
    snapshot <- api("hydrateSavedWorkflowCheckpoint")(store, readRDS(file.path(store_dir, "dia-state-head.rds")))
    products <- lapply(snapshot$products, \(reference) kitFingerprintReference(store, reference, api))
    list(head_sha256 = head$head_sha256, states = result, products = products)
}

#' Verify one KIT member outside its measured slice: binding, build, steps, head and fingerprints.
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
    receipt <- list(status = "starting", run_dir = run_dir, member_id = case$id, selector = selection)
    on.exit({
        if (identical(receipt$status, "starting")) receipt$status <- "failed"
        saveRDS(receipt, path)
    }, add = TRUE)
    operation_path <- file.path(owned, "operation-receipt.rds")
    receipt$operation_sha256 <- digest::digest(file = operation_path, algo = "sha256")
    operation <- readRDS(operation_path)
    journey <- readRDS(file.path(owned, "journey-result.rds"))
    stopifnot(identical(operation$binding, case$binding), identical(operation$config$build, case$expected_build),
        identical(journey$corpus, case$corpus))
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
    stopifnot(identical(operation$status, "complete"), all(steps %in% names(status)), all(unlist(status[steps]) == "complete"))
    store_dir <- journey$store
    stopifnot(!is.null(store_dir), dir.exists(store_dir))
    recomputed <- kitFingerprintStore(store_dir, api)
    stopifnot(identical(recomputed$head_sha256, journey$final_head_sha256))
    mismatches <- character()
    for (step in names(journey$populations)) {
        population <- journey$populations[[step]]
        if (is.null(population$state)) next
        expected <- population$fingerprint
        actual <- recomputed$states[[population$state]]
        if (!identical(expected, actual)) mismatches <- c(mismatches, step)
    }
    stopifnot(!length(mismatches))
    files <- kitFingerprintFiles(journey$project)
    stopifnot(identical(files[c("path", "bytes", "sha256", "method")], journey$files[c("path", "bytes", "sha256", "method")]))
    receipt$populations <- lapply(journey$populations, \(p) list(state = p$state, entity_id = p$entity_id,
        rows = p$rows, fingerprints = length(p$fingerprint$components) + length(p$fingerprint$support) +
            length(p$fingerprint$products) + 1L))
    receipt$files <- nrow(files)
    receipt$products <- length(recomputed$products)
    receipt$status <- "complete"
    receipt
}

`%||%` <- function(x, y) if (is.null(x)) y else x
