# DIA-003 installed production-browser driver primitives.
# No workload is launched when this file is sourced.

#' Read the separate, explicit capacity-only scope without changing frozen gates.
#' @param source_root Repository root containing the versioned fixture contract.
#' @return Capacity specification with fixed route/order; not adoption eligibility.
#' @noRd
readCheckpointGuiCapacitySpec <- function(source_root) {
    path <- file.path(source_root, "tests", "testdata", "workflow-checkpoints", "gui-capacity.json")
    spec <- jsonlite::read_json(path)
    if (!identical(spec$schema, "multischolar.dia_gui_capacity.v1") ||
        !identical(spec$scope, "full-journey-capacity-v1") || !identical(spec$route, "iq") ||
        !identical(spec$treatment, "qs2_duckdb") || !identical(spec$promotion_authorized, FALSE) ||
        !identical(unlist(spec$run_order), c("small:resident", "small:qs2_duckdb",
            "oom_reference:qs2_duckdb", "report_1gb:qs2_duckdb", "report_2gb:qs2_duckdb"))) {
        stop("unsupported GUI capacity scope")
    }
    dimensions <- lapply(spec$recipes, \(recipe) unlist(recipe[c("n_proteins", "n_samples",
        "peptides_per_protein", "forms_per_peptide", "precursor_rows")], use.names = FALSE))
    if (!identical(vapply(spec$recipes, \(recipe) recipe$id, character(1)),
        c("small", "report_1gb", "report_2gb")) ||
        !isTRUE(all.equal(dimensions, list(c(256, 12, 4, 2, 24576),
            c(10000, 48, 8, 2, 7680000), c(10000, 48, 16, 2, 15360000)), tolerance = 0))) {
        stop("capacity dimensions differ from the approved recipes")
    }
    spec
}

#' Compare one historical record using the established volatile-field policy.
#' @param record Complete historical record.
#' @param root Verified member root or roots for fixture replay.
#' @return Scientific record with only existing volatile metadata normalized.
#' @noRd
comparableCheckpointGuiRecord <- function(record, root) {
    namespace <- asNamespace("MultiScholaR")
    prepare <- get("prepareDIAStateRecord", namespace)
    substantive <- get(".peptideQcSubstantiveAudit", namespace)
    normalize_paths <- function(value, root) {
        if (is.character(value)) {
            original_attributes <- attributes(value)
            for (member_root in root) {
                selected <- !is.na(value) & startsWith(value, member_root)
                value[selected] <- paste0("<MEMBER>", substring(value[selected], nchar(member_root) + 1L))
            }
            value <- gsub("/filtering_qc/[0-9]{8}_[0-9]{6}/", "/filtering_qc/<RUN_TIME>/", value)
            attributes(value) <- original_attributes
            return(value)
        }
        if (is.list(value)) {
            for (index in seq_along(value)) value[index] <- list(normalize_paths(value[[index]], root))
        }
        value
    }
    comparable <- function(record, root) {
        record <- prepare(record)
        record$timestamp <- NULL
        if (isS4(record$data)) {
            args <- methods::slot(record$data, "args")
            if (!is.null(args$peptide_qc_audit)) args$peptide_qc_audit <- substantive(record$data)
            # DIA-025: this exact UI field records Sys.time(), not a scientific
            # parameter. Normalize its value only, preserving presence/type.
            ui_time <- args$daAnalysisUI$timestamp
            if (inherits(ui_time, "POSIXct") && length(ui_time) == 1L && !is.na(ui_time)) {
                ui_time[] <- 0
                args$daAnalysisUI$timestamp <- ui_time
            }
            # Approved scheduling-only differences; scientific settings remain
            # fully compared, including requested core counts in global config.
            for (section in c("rollUpPrecursorToPeptide", "peptideIntensityFiltering",
                "filterMinNumPeptidesPerProtein", "filterMinNumPeptidesPerSample",
                "removePeptidesWithOnlyOneReplicate", "peptideMissingValueImputation",
                "removeProteinsWithOnlyOneReplicate", "removeProteinsWithOnlyOneReplicateFn")) {
                worker <- args[[section]]$core_utilisation
                if (inherits(worker, "DIAWorkerSpec") || identical(worker, NA))
                    args[[section]]$core_utilisation <- "<EXECUTION_BACKEND>"
            }
            methods::slot(record$data, "args") <- normalize_paths(args, root)
        }
        record$config <- normalize_paths(record$config, root)
        if (is.list(record$audit_metadata)) record$audit_metadata$timestamp_utc <- NULL
        record
    }
    comparable(record, root)
}

#' Compare the small actual GUI histories and figures in an independent verifier.
#' @param control_dir,candidate_dir Completed small browser member directories.
#' @param candidate_session Already decoded complete candidate legacy session.
#' @return Bounded scientific/output comparison receipt, not a performance result.
#' @noRd
compareCheckpointGuiCapacityMembers <- function(control_dir, candidate_dir, candidate_session) {
    control <- jsonlite::read_json(file.path(control_dir, "smoke-receipt.json"))
    candidate <- jsonlite::read_json(file.path(candidate_dir, "smoke-receipt.json"))
    if (!identical(control$capacity$treatment, "resident") ||
        !identical(candidate$capacity$treatment, "qs2_duckdb") ||
        !identical(control$input_sha256, candidate$input_sha256) ||
        file.size(control$inputs$report) > 20000000 ||
        !identical(control$installation$package_path, candidate$installation$package_path)) {
        stop("comparison requires the matched small GUI pair")
    }
    comparable <- comparableCheckpointGuiRecord
    original <- readRDS(file.path(control_dir, "dia_baseline_smoke", "scripts", "proteomics",
        "filtered_session_data_latest.rds"))
    if (!identical(names(original$r6_complete_states), names(candidate_session$r6_complete_states)) ||
        !identical(original$r6_state_history, candidate_session$r6_state_history) ||
        !identical(original$r6_current_state_name, candidate_session$r6_current_state_name)) {
        stop("small GUI history topology differs")
    }
    sizes <- list()
    for (name in names(original$r6_complete_states)) {
        expected <- comparable(original$r6_complete_states[[name]], control_dir)
        actual <- comparable(candidate_session$r6_complete_states[[name]], candidate_dir)
        comparison <- all.equal(actual, expected, tolerance = 1e-12, check.attributes = TRUE)
        if (!isTRUE(comparison)) stop(paste("small GUI state differs:", name, paste(comparison, collapse = "; ")))
        if (isS4(actual$data)) {
            field <- if (methods::is(actual$data, "PeptideQuantitativeData")) "peptide_data" else "protein_quant_table"
            a <- methods::slot(actual$data, field)
            b <- methods::slot(expected$data, field)
            if (!identical(vapply(a, typeof, character(1)), vapply(b, typeof, character(1))) ||
                !identical(lapply(a, is.na), lapply(b, is.na)) ||
                !identical(lapply(a, is.nan), lapply(b, is.nan))) stop("small GUI table storage/missingness differs")
            sizes[[name]] <- list(s4_object_size_bytes = as.numeric(utils::object.size(actual$data)),
                table_object_size_bytes = as.numeric(utils::object.size(a)), rows = nrow(a), columns = ncol(a))
        }
    }
    png_inventory <- function(root) {
        project <- file.path(root, "dia_baseline_smoke")
        paths <- list.files(project, "[.]png$", recursive = TRUE, full.names = TRUE)
        names <- gsub("/filtering_qc/[0-9]{8}_[0-9]{6}/", "/filtering_qc/<RUN_TIME>/",
            substring(paths, nchar(project) + 1L))
        if (!length(paths) || anyDuplicated(names)) stop("missing or ambiguous scientific PNG inventory")
        hashes <- stats::setNames(vapply(paths, \(path) digest::digest(file = path, algo = "sha256"),
            character(1)), names)
        hashes[order(names(hashes))]
    }
    if (!identical(png_inventory(control_dir), png_inventory(candidate_dir))) {
        stop("small GUI scientific PNG outputs differ")
    }
    list(status = "passed_small_gui_science", states = names(sizes), sizes = sizes,
        png_contents_identical = TRUE, tolerance = 1e-12,
        allowed_volatile_fields = c("state timestamp", "audit timestamp_utc", "member-root metadata paths",
            "data@args$daAnalysisUI$timestamp wall-clock value; presence/type preserved",
            "filtering_qc timestamp directory", "existing worker/parser portability conversion",
            "approved DIA core_utilisation execution backend; scientific parameters retained"),
        promotion_authorized = FALSE)
}

#' Read the separately frozen sample-scaled RUV/Pearson diagnostic recipes.
#' @param source_root Repository root, never an application runtime path.
#' @return Diagnostic-only specification; historical capacity scenarios are unchanged.
#' @noRd
readCheckpointRuvPearsonSpec <- function(source_root) {
    spec <- jsonlite::read_json(file.path(source_root, "tests", "testdata",
        "workflow-checkpoints", "gui-ruv-pearson-diagnostic.json"))
    stopifnot(identical(spec$scope, "ruv-pearson-memory-v1"), identical(spec$promotion_authorized, FALSE),
        identical(spec$runtime_rdb_sha256, "153db129f2911ff008dbeda2a77b8d229c0816c726d2d9e368e71b7be9e32044"),
        identical(spec$export_repair_rdb_sha256, "8769e614e8b8126c3ed1ccaa6d491fbf857ffad65c988fa595e3cd33f2301a76"),
        identical(spec$integrated_candidate_rdb_sha256, "060455e405f7b7dbd5ee723b2f3e3ba8eff9dc4a4f4833f575ff463696e240d0"),
        identical(spec$final_launch_headroom_bytes, 17179869184),
        identical(unlist(spec$native_first_96_rdb_sha256, use.names = FALSE), c(
            "f529ddc96a9bc22883b2f984ba97edb156bd32284cf1add528cbbadfe0725399",
            "c859941275a7ff2e1986aea235ff0be622a9fdd49d26558979da99c67d5be3de")),
        identical(unlist(spec$native_first_300_rdb_sha256, use.names = FALSE), c(
            "fd69f611d68392b6583187ef0004021f4a589dbb0cac895d7bae165cb42f99db",
            "1dbc963d9eedfc8eba2b168df928cf470b14e04a4e2928ccf781dba7f2814e7b",
            "c55102a392a3384021f979af0536d734a2863350a5321472ffecad9cbfd7f7c8",
            "66aa767cdcc7f7f4628815c758e829d1085b4379864ae009af067f51e2244968",
            "941dc2a7f94175c40139d58e660839ae038ee709d840936871f23ae56810749c")),
        identical(vapply(spec$recipes, `[[`, character(1), "id"), c("ruv_small", "ruv_1gb", "ruv_2gb")),
        identical(vapply(spec$recipes, function(x) as.numeric(x$n_samples), numeric(1)), c(24, 96, 192)),
        spec$limits$memory_max_bytes == 19327352832, spec$limits$runtime_max_seconds == 3600)
    for (samples in c(96, 300, 500)) {
        stress <- spec[[paste0("stress_", samples)]]
        stopifnot(identical(stress$runtime_rdb_sha256,
            "74a7dc9b50821e49e1a6d8c9a4366bdcc879997d61d06a9876acd55387404a08"),
            identical(stress$recipe$id, if (samples == 96) "ruv_1gb" else paste0("ruv_", samples)),
            identical(as.numeric(unlist(stress$recipe[c("n_proteins", "n_samples",
                "peptides_per_protein", "forms_per_peptide", "history_depth", "annotation_chars",
                "seed", "precursor_rows")])), c(10000, samples, 4, 2, 16, 96, 902002, samples * 80000)),
            identical(as.numeric(unlist(stress$limits[c("memory_max_bytes", "early_stop_bytes",
                "swap_max_bytes", "runtime_max_seconds", "launch_headroom_bytes", "host_stop_bytes")])),
                c(17179869184, 16642998272, 0, 3600, 19327352832, 2147483648)))
    }
    stopifnot(identical(spec$stress_96$recipe, spec$recipes[[2L]]),
        spec$integrated_candidate_rdb_sha256 %in% spec$native_session_rdb_sha256)
    spec
}

#' Prepare a bounded synthetic report/companions, or bind the retained OOM inputs.
#' @param library_path,run_dir Installed library and unique guarded preparation directory.
#' @param selection Explicit capacity recipe selector.
#' @param source_root Repository root; generator code is not package runtime.
#' @return Small input receipt; fixture creation remains outside journey timing.
#' @noRd
prepareCheckpointGuiCapacityInputs <- function(library_path, run_dir, selection, source_root) {
    ruv_diagnostic <- selection %in% c("ruv_small", "ruv_1gb", "ruv_2gb", "ruv_300", "ruv_500")
    if (!ruv_diagnostic && !selection %in% c("small", "oom_reference", "report_1gb", "report_2gb")) {
        stop("unsupported capacity input selection")
    }
    spec <- if (ruv_diagnostic) readCheckpointRuvPearsonSpec(source_root) else readCheckpointGuiCapacitySpec(source_root)
    helpers <- new.env(parent = baseenv())
    helper_path <- file.path(source_root, "tests", "testthat", "helper-workflow-checkpoint-fixtures.R")
    sys.source(helper_path, helpers)
    root <- dirname(normalizePath(run_dir, mustWork = TRUE))
    helpers$validateWorkflowFixtureRoot(root, source_root)
    receipt_path <- file.path(run_dir, "inputs.json")
    if (file.exists(receipt_path)) stop("capacity preparation receipt already exists")
    receipt <- list(kind = if (ruv_diagnostic) "gui_ruv_pearson_input_preparation" else "gui_capacity_input_preparation", scope = spec$scope,
        selection = selection, status = "running", synthetic = TRUE, promotion_authorized = FALSE,
        helper_sha256 = digest::digest(file = helper_path, algo = "sha256"),
        contract_sha256 = digest::digest(file = file.path(source_root, "tests", "testdata",
            "workflow-checkpoints", if (ruv_diagnostic) "gui-ruv-pearson-diagnostic.json" else "gui-capacity.json"), algo = "sha256"))
    flush <- function() jsonlite::write_json(receipt, receipt_path, auto_unbox = TRUE,
        pretty = TRUE, null = "null", digits = NA)
    flush()
    on.exit(flush(), add = TRUE)
    tryCatch({
        if (selection == "oom_reference") {
            files <- lapply(spec$reference$files, \(relative) normalizePath(file.path(root, relative), mustWork = TRUE))
            if (!all(vapply(files, \(path) startsWith(path, paste0(root, "/")), logical(1)))) {
                stop("capacity reference input escapes proof root")
            }
            hashes <- lapply(files, \(path) digest::digest(file = path, algo = "sha256"))
            if (!identical(hashes, spec$reference$sha256)) stop("retained OOM inputs differ")
            receipt$recipe <- helpers$readWorkflowCheckpointRecipe("representative", source_root)
            receipt$failed_members <- spec$reference$failed_members
        } else {
            matches <- if (selection %in% c("ruv_300", "ruv_500"))
                list(spec[[sub("^ruv_", "stress_", selection)]]$recipe) else
                Filter(\(recipe) identical(recipe$id, selection), spec$recipes)
            if (length(matches) != 1L) stop("capacity recipe must be unique")
            recipe <- helpers$validateWorkflowCheckpointRecipe(matches[[1L]], TRUE, root, source_root)
            files <- as.list(stats::setNames(file.path(run_dir, c("report.tsv", "proteins.fasta", "design.tsv")),
                c("report", "fasta", "design")))
            if (any(file.exists(unlist(files)))) stop("capacity inputs already exist")
            per_protein <- recipe$n_samples * recipe$peptides_per_protein * recipe$forms_per_peptide
            proteins_per_chunk <- max(1L, floor(200000 / per_protein))
            receipt$rows_written <- 0
            receipt$max_chunk_rows <- 0
            for (start in seq.int(1L, recipe$n_proteins, by = proteins_per_chunk)) {
                ids <- seq.int(start, min(recipe$n_proteins, start + proteins_per_chunk - 1L))
                chunk <- helpers$makeWorkflowPrecursorChunk(recipe, ids)
                vroom::vroom_write(chunk, files$report, delim = "\t", append = start != 1L,
                    col_names = start == 1L, progress = FALSE)
                receipt$rows_written <- receipt$rows_written + nrow(chunk)
                receipt$max_chunk_rows <- max(receipt$max_chunk_rows, nrow(chunk))
                rm(chunk)
                flush()
            }
            stopifnot(receipt$rows_written == recipe$precursor_rows)
            design <- helpers$makeWorkflowCheckpointDesign(recipe$n_samples)
            design$sample <- design$Run
            vroom::vroom_write(design, files$design, progress = FALSE)
            fasta <- file(files$fasta, open = "w")
            tryCatch({
                for (protein in seq_len(recipe$n_proteins)) {
                    peptides <- helpers$workflowCheckpointSequences(
                        (protein - 1) * recipe$peptides_per_protein + seq_len(recipe$peptides_per_protein))
                    writeLines(c(sprintf(">P%05d synthetic diagnostic sequence", protein),
                        paste(peptides, collapse = "")), fasta)
                }
            }, finally = close(fasta))
            receipt$recipe <- recipe
            hashes <- lapply(files, \(path) digest::digest(file = path, algo = "sha256"))
        }
        receipt$files <- files
        receipt$sha256 <- hashes
        receipt$bytes <- lapply(files, file.size)
        receipt$protein_numeric_bytes <- with(receipt$recipe, n_proteins * n_samples * 8)
        receipt$decoded_size_note <- "No full report pre-hydration; actual S4/table sizes belong to execution/verification."
        receipt$status <- "prepared"
        receipt
    }, error = function(error) {
        receipt$status <<- "failed"
        receipt$error <<- conditionMessage(error)
        stop(error)
    })
}

#' Bound waits by the remaining lifetime of the 900-second service process.
#' @param control Use the shorter 30-second missing-control readiness check.
#' @param elapsed_seconds Optional deterministic clock value for tests.
#' @return Positive integer milliseconds; exhaustion is an explicit error.
#' @noRd
checkpointWaitBudget <- function(control = FALSE, elapsed_seconds = NULL) {
    if (is.null(elapsed_seconds)) elapsed_seconds <- unname(proc.time()[["elapsed"]])
    if (!is.logical(control) || length(control) != 1L || is.na(control) ||
        !is.numeric(elapsed_seconds) || length(elapsed_seconds) != 1L ||
        !is.finite(elapsed_seconds) || elapsed_seconds < 0) stop("invalid wait budget inputs")
    ceiling <- getOption("multischolar.dia_memory_audit_seconds", 900)
    if (!is.numeric(ceiling) || length(ceiling) != 1L ||
        !ceiling %in% c(900, 3600, if (isTRUE(getOption("multischolar.dia_extended120"))) 7200))
        stop("invalid diagnostic time ceiling")
    remaining <- ceiling - elapsed_seconds
    if (remaining <= 0) stop("DIA service time budget exhausted")
    milliseconds <- min(remaining * 1000, if (control) 30000 else ceiling * 1000)
    as.integer(max(1, floor(milliseconds)))
}

#' Construct a cold installed-app launch function for shinytest2.
#' @param library_path Explicit installed library.
#' @param identity_root Optional request-owned directory for bounded runtime identity.
#' @param capacity_treatment Explicit experimental treatment, otherwise NULL.
#' @return Zero-argument launch function containing only scalar path bindings.
#' @noRd
checkpointInstalledApp <- function(library_path, identity_root = NULL, capacity_treatment = NULL,
    session_format = "legacy", reclaim_memory = FALSE) {
    library_path <- normalizePath(library_path, mustWork = TRUE)
    package_path <- normalizePath(file.path(library_path, "MultiScholaR"), mustWork = TRUE)
    if (!is.null(identity_root)) identity_root <- normalizePath(identity_root, mustWork = TRUE)
    if (!session_format %in% c("legacy", "native_v1")) stop("invalid session format")
    if (!is.null(capacity_treatment) && (!capacity_treatment %in% c("resident", "qs2_duckdb") ||
        is.null(identity_root))) stop("invalid capacity app activation")
    if (!is.logical(reclaim_memory) || length(reclaim_memory) != 1L || is.na(reclaim_memory) ||
        reclaim_memory && (!identical(session_format, "native_v1") ||
            !identical(capacity_treatment, "qs2_duckdb"))) stop("invalid reclamation activation")
    launch <- function() {
        if ("MultiScholaR" %in% loadedNamespaces()) {
            stop("MultiScholaR was preloaded outside the app startup boundary")
        }
        .libPaths(unique(c(fixed_library, .libPaths())))
        options(multischolar.test_mode = FALSE,
            multischolar.capture_test_checkpoints = FALSE,
            multischolar.dia_capacity_trial = NULL)
        if (identical(fixed_capacity, "qs2_duckdb")) {
            options(multischolar.dia_capacity_trial = list(scope = "full-journey-capacity-v1",
                route = "iq", treatment = fixed_capacity, receipt_dir = fixed_identity_root,
                session_format = fixed_session_format, reclaim_memory = fixed_reclaim))
        }
        Sys.setenv(MULTISCHOLAR_TEST_MODE = "false")
        namespace <- loadNamespace("MultiScholaR", lib.loc = fixed_library)
        if (!identical(normalizePath(getNamespaceInfo(namespace, "path")), fixed_package)) {
            stop("app loaded a different installed package")
        }
        if (!dir.exists(system.file("shiny", "www", package = "MultiScholaR"))) {
            stop("installed static assets are absent; source fallback is forbidden")
        }
        if (getExportedValue("MultiScholaR", "is_test_mode")()) {
            stop("scientific test mode must remain disabled")
        }
        record_runtime <- function(stage) {
            if (is.null(fixed_identity_root)) return(invisible(NULL))
            path <- file.path(fixed_identity_root, paste0("installed-", stage, "-identity.json"))
            if (file.exists(path)) stop("refusing to replace installed runtime identity")
            namespaces <- sort(loadedNamespaces())
            versions <- stats::setNames(vapply(namespaces,
                \(name) as.character(getNamespaceVersion(name)), character(1)), namespaces)
            dlls <- getLoadedDLLs()
            native_paths <- vapply(dlls, \(dll) dll[["path"]], character(1))
            native_bytes <- vapply(native_paths, \(path) if (file.exists(path))
                as.numeric(file.size(path)) else NA_real_, numeric(1))
            own_root <- file.path(fixed_package, "libs")
            own_paths <- sort(list.files(own_root, pattern = "[.](so|dll|dylib)$", recursive = TRUE, full.names = TRUE))
            own_hashes <- stats::setNames(lapply(own_paths, function(path) digest::digest(file = path, algo = "sha256")),
                substring(own_paths, nchar(own_root) + 2L))
            if (length(own_paths) && (!"MultiScholaR" %in% names(native_paths) ||
                !normalizePath(native_paths[["MultiScholaR"]]) %in% normalizePath(own_paths)))
                stop("package native library was not loaded from the requested installation")
            jsonlite::write_json(list(stage = stage, process_id = Sys.getpid(),
                package_native_libraries = own_hashes,
                package_path = fixed_package, r_version = R.version.string,
                namespace_versions = as.list(versions), native_paths = as.list(native_paths),
                native_bytes = as.list(native_bytes),
                native_note = "package libraries have SHA-256 fingerprints; dependency libraries record loaded paths and sizes"),
                path, auto_unbox = TRUE, pretty = TRUE, null = "null")
            if (is.na(file.size(path)) || file.size(path) > 65536) {
                stop("installed runtime identity exceeded 64 KiB")
            }
            invisible(path)
        }
        record_runtime("load")
        app <- getExportedValue("MultiScholaR", "run_app")(test_mode = FALSE)
        if (!inherits(app, "shiny.appobj")) stop("installed launcher did not return a Shiny app")
        shiny::onStop(\() record_runtime("exit"), session = NULL)
        app
    }
    environment(launch) <- list2env(list(fixed_library = library_path,
        fixed_package = package_path, fixed_identity_root = identity_root,
        fixed_capacity = capacity_treatment, fixed_session_format = session_format,
        fixed_reclaim = reclaim_memory), parent = baseenv())
    launch
}

#' Wrap only a diagnostic app launch with CPU/GC/heap sampling.
#' @param library_path Explicit installed control library.
#' @param profile_path Fresh output file in the request-owned directory.
#' @return Launch closure with no scientific objects; never used for timed members.
#' @noRd
checkpointProfiledApp <- function(library_path, profile_path, identity_root = NULL) {
    if (!is.character(profile_path) || length(profile_path) != 1L || is.na(profile_path) ||
        !dir.exists(dirname(profile_path)) || file.exists(profile_path)) stop("invalid or existing profile path")
    launch <- checkpointInstalledApp(library_path, identity_root)
    profiled <- function() {
        if (file.exists(fixed_profile)) stop("refusing to overwrite application profile")
        if (!isTRUE(unname(capabilities("profmem")))) stop("R memory profiling is unavailable")
        utils::Rprof(fixed_profile, interval = 0.02, memory.profiling = TRUE,
            gc.profiling = TRUE, event = "cpu")
        ready <- FALSE
        on.exit(if (!ready) utils::Rprof(NULL), add = TRUE)
        app <- fixed_launch()
        if (!inherits(app, "shiny.appobj")) stop("profiled launcher must return a Shiny app")
        # AppDriver starts runApp after this constructor returns. Keep sampling
        # through that event loop; stop at application shutdown, not construction.
        shiny::onStop(function() utils::Rprof(NULL), session = NULL)
        ready <- TRUE
        app
    }
    environment(profiled) <- list2env(list(fixed_profile = profile_path, fixed_launch = launch),
        parent = baseenv())
    profiled
}

#' Create the production app browser, without source/test-mode fallback.
#' @param library_path Explicit installed package library.
#' @param run_dir Existing request-owned directory for driver artifacts.
#' @param seed Fixed input-independent browser seed.
#' @param profile_application Enable diagnostic-only sampling of the app R process.
#' @param capacity_treatment Explicit experimental treatment; profiling is separate.
#' @param diagnostic_source Harness-only diagnostic helper for the explicit 16 GiB extension.
#' @param diagnostic_profile NULL or the explicit RUV/Pearson stress profile.
#' @param allocator_probe Explicit small diagnostic-only allocator measurements.
#' @param reclaim_memory Opt in to completed-action reclamation for this launch only.
#' @param design_memory Explicit attribution run that stops before Design Save.
#' @return Driver and canonical artifact root; caller owns stop and artifact collection.
#' @noRd
newCheckpointBrowser <- function(library_path, run_dir, seed = 20260428L, profile_application = FALSE,
    capacity_treatment = NULL, diagnostic_source = NULL, diagnostic_profile = NULL, session_format = "legacy", allocator_probe = FALSE, reclaim_memory = FALSE, design_memory = FALSE) {
    if (!is.null(diagnostic_profile) && !identical(diagnostic_profile, "ruv_pearson")) stop("invalid diagnostic profile")
    if (!is.logical(profile_application) || length(profile_application) != 1L || is.na(profile_application)) {
        stop("profile_application must be TRUE or FALSE")
    }
    if (!dir.exists(run_dir)) stop("run directory must already exist")
    if (!requireNamespace("shinytest2", quietly = TRUE) ||
        !requireNamespace("chromote", quietly = TRUE)) {
        stop("shinytest2 and chromote are required; no browser skip is allowed")
    }
    browser <- chromote::find_chrome()
    if (length(browser) != 1L || !nzchar(browser) || !file.exists(browser)) {
        stop("a real browser executable is required")
    }
    capacity_diagnostic <- isTRUE(profile_application) && identical(capacity_treatment, "qs2_duckdb") &&
        !is.null(diagnostic_source) && file.exists(diagnostic_source)
    if (isTRUE(profile_application) && !is.null(capacity_treatment) && !capacity_diagnostic)
        stop("capacity timing cannot enable profiling")
    app <- if (capacity_diagnostic) {
        launch <- if (isTRUE(reclaim_memory)) checkpointInstalledApp(library_path, run_dir,
            capacity_treatment, session_format, reclaim_memory = TRUE) else
            checkpointInstalledApp(library_path, run_dir, capacity_treatment, session_format)
        wrapped <- function() {
            helpers <- new.env(parent = baseenv())
            sys.source(fixed_helper, helpers)
            diagnostic <- helpers$startCheckpointGuiDiagnostic(fixed_directory, collect_heap = FALSE)
            if (isTRUE(fixed_allocator_probe) || isTRUE(fixed_design_memory))
                helpers$enableCheckpointAllocatorProbe(diagnostic, arena_details = isTRUE(fixed_design_memory))
            diagnostic$design_memory <- isTRUE(fixed_design_memory)
            ready <- FALSE
            on.exit(if (!ready) helpers$stopCheckpointGuiDiagnostic(diagnostic), add = TRUE)
            app <- fixed_launch()
            if (isTRUE(fixed_design_memory)) helpers$checkpointAllocatorBoundary(diagnostic, "startup")
            helpers$installCheckpointGuiTraces(diagnostic, asNamespace("MultiScholaR"),
                functions = if (isTRUE(fixed_design_memory)) c("initializeDIACapacitySession",
                    "importProtImportDataByFormat", "applyProtImportResultToWorkflow") else NULL)
            if (!isTRUE(fixed_design_memory) && identical(fixed_profile, "ruv_pearson"))
                helpers$installCheckpointRuvPearsonTraces(diagnostic, asNamespace("MultiScholaR"))
            shiny::onStop(function() helpers$stopCheckpointGuiDiagnostic(diagnostic), session = NULL)
            ready <- TRUE
            app
        }
        environment(wrapped) <- list2env(list(fixed_launch = launch,
            fixed_helper = normalizePath(diagnostic_source), fixed_directory = normalizePath(run_dir),
            fixed_profile = diagnostic_profile, fixed_allocator_probe = isTRUE(allocator_probe),
            fixed_design_memory = isTRUE(design_memory)),
            parent = baseenv())
        wrapped
    } else if (isTRUE(profile_application)) checkpointProfiledApp(library_path,
        file.path(run_dir, "application.Rprof"), run_dir) else
        if (isTRUE(reclaim_memory)) checkpointInstalledApp(library_path, run_dir,
            capacity_treatment, session_format, reclaim_memory = TRUE) else
            checkpointInstalledApp(library_path, run_dir, capacity_treatment, session_format)
    driver <- shinytest2::AppDriver$new(
        app_dir = app,
        name = "dia-installed-baseline", seed = seed,
        timeout = checkpointWaitBudget(), load_timeout = checkpointWaitBudget(),
        width = 1440L, height = 1100L,
        options = list(multischolar.test_mode = FALSE,
            multischolar.capture_test_checkpoints = FALSE)
    )
    list(driver = driver, artifact_root = normalizePath(run_dir, mustWork = TRUE))
}

#' Wait for browser completion and fail on session disconnection.
#' @param driver Live browser; diagnostics use short transport commands.
#' @param script Existing trusted harness completion expression.
#' @param timeout Maximum wait in milliseconds, bounded by the remaining service time.
#' @param stable_ms Required continuous true interval; used only for idle waits.
#' @return NULL invisibly; transport loss is never successful completion.
#' @noRd
checkpointWaitForJs <- function(driver, script, timeout = NULL, stable_ms = 0) {
    if (is.null(getOption("multischolar.dia_memory_audit_directory"))) {
        guarded <- sprintf(paste0("(function(){",
            "if(document.getElementById('shiny-disconnected-overlay') || ",
            "(window.Shiny && Shiny.shinyapp && Shiny.shinyapp.$socket && ",
            "Shiny.shinyapp.$socket.readyState!==1))throw new Error('Shiny session disconnected');",
            "return Boolean((0,eval)(%s));})()"), jsonlite::toJSON(script, auto_unbox = TRUE))
        return(driver$wait_for_js(guarded, timeout = timeout))
    }
    budget <- checkpointWaitBudget()
    if (!is.null(timeout)) budget <- min(budget, timeout)
    deadline <- unname(proc.time()[["elapsed"]]) + budget / 1000
    request <- sprintf(paste0("JSON.stringify({ok:Boolean((0,eval)(%s)),",
        "connected:!!(window.Shiny && Shiny.shinyapp && Shiny.shinyapp.$socket && ",
        "Shiny.shinyapp.$socket.readyState===1)})"), jsonlite::toJSON(script, auto_unbox = TRUE))
    stable_since <- NULL
    repeat {
        now <- unname(proc.time()[["elapsed"]])
        if (now >= deadline) stop("diagnostic browser condition timed out")
        value <- driver$get_js(request, timeout = min(5000, (deadline - now) * 1000))
        if (!is.character(value) || length(value) != 1L || is.na(value))
            stop("diagnostic browser transport returned no response")
        state <- jsonlite::fromJSON(value)
        if (!isTRUE(state$connected)) stop("diagnostic Shiny session disconnected")
        if (isTRUE(state$ok)) {
            if (is.null(stable_since)) stable_since <- now
            if ((now - stable_since) * 1000 >= stable_ms) return(invisible(NULL))
        } else stable_since <- NULL
        Sys.sleep(0.25)
    }
}

#' Preserve normal idle waits; use bounded polling for the memory diagnostic.
#' @param driver Live browser.
#' @param timeout Maximum wait in milliseconds.
#' @return NULL invisibly.
#' @noRd
checkpointWaitForIdle <- function(driver, timeout = NULL) {
    if (is.null(getOption("multischolar.dia_memory_audit_directory"))) {
        return(driver$wait_for_idle(timeout = timeout))
    }
    checkpointWaitForJs(driver, "window.shinytest2 && window.shinytest2.busy !== true",
        timeout = timeout, stable_ms = 200)
}

#' Build a production stepper observation expression.
#' @param label Visible step label, e.g. Import or Normalize.
#' @return JavaScript boolean expression; no test-only digest access.
#' @noRd
checkpointStepCompleteScript <- function(label) {
    if (!is.character(label) || length(label) != 1L || is.na(label) || !nzchar(label)) {
        stop("a single step label is required")
    }
    encoded <- jsonlite::toJSON(label, auto_unbox = TRUE)
    sprintf(paste0(
        "(function(){const root=document.getElementById('proteomics_workflow-workflow_progress');",
        "if(!root)return false;return Array.from(root.querySelectorAll('.stepper-step'))",
        ".some(function(step){const label=step.querySelector('.stepper-label');",
        "return label && label.textContent.trim()===%s && ",
        "step.classList.contains('step-completed');});})()"), encoded)
}

#' Wait for a production completion signal, not output publication.
#' @param driver Live AppDriver.
#' @param label Visible workflow step label.
#' @param timeout Bounded wait in milliseconds.
#' @return NULL invisibly.
#' @noRd
waitCheckpointStepComplete <- function(driver, label, timeout = NULL) {
    timeout <- if (is.null(timeout)) checkpointWaitBudget() else min(timeout, checkpointWaitBudget())
    checkpointWaitForJs(driver, checkpointStepCompleteScript(label), timeout = timeout)
    invisible(NULL)
}

#' List the original peptide actions in current UI order.
#' @return Ordered action/result/state records; no scientific calls.
#' @noRd
checkpointPeptideActions <- function() {
    list(
        list(tab = "Q-Value Filter", action = "qvalue_filter-apply_qvalue_filter",
            result = "qvalue_filter-qvalue_results", state = "qvalue_filtered"),
        list(tab = "Protein Evidence", action = "protein_peptide_filter-apply_protein_peptide_filter",
            result = "protein_peptide_filter-protein_peptida_results", state = "protein_peptide_filtered"),
        list(tab = "Precursor Rollup", action = "rollup-apply_rollup",
            result = "rollup-rollup_results", state = "precursor_rollup"),
        list(tab = "Intensity Filter", action = "intensity_filter-apply_intensity_filter",
            result = "intensity_filter-intensity_results", state = "intensity_filtered"),
        list(tab = "Sample Quality", action = "sample_filter-apply_sample_filter",
            result = "sample_filter-sample_results", state = "sample_filtered"),
        list(tab = "Replicate Filter", action = "replicate_filter-apply_replicate_filter",
            result = "replicate_filter-replicate_results", state = "replicate_filtered"),
        list(tab = "Imputation", action = "imputation-apply_imputation",
            result = "imputation-imputation_results", state = "imputed")
    )
}

#' Record a GUI request boundary without reading any scientific payload.
#' @param stage Fixed action/control label; input values are never recorded.
#' @return NULL invisibly; disabled outside the explicit memory diagnostic.
#' @noRd
checkpointDriverMemoryBoundary <- function(stage) {
    directory <- getOption("multischolar.dia_memory_audit_directory") %||%
        getOption("multischolar.dia_progress_directory")
    if (is.null(directory)) return(invisible(NULL))
    path <- file.path(directory, "driver-phases.jsonl")
    if (file.exists(path) && file.size(path) >= 64 * 1024^2)
        stop("driver memory-boundary log exceeded 64 MiB")
    connection <- file(path, "a")
    on.exit(close(connection), add = TRUE)
    writeLines(jsonlite::toJSON(list(time = as.numeric(Sys.time()), pid = Sys.getpid(),
        stage = stage), auto_unbox = TRUE, digits = NA), connection)
    invisible(NULL)
}

#' Click a real production control by its stable test identifier.
#' @param driver Browser driver.
#' @param test_id Existing data-testid value.
#' @return NULL invisibly.
#' @noRd
checkpointClickControl <- function(driver, test_id) {
    checkpointDriverMemoryBoundary(paste0("enter:control:", test_id))
    on.exit(checkpointDriverMemoryBoundary(paste0("exit:control:", test_id)), add = TRUE)
    id <- jsonlite::toJSON(test_id, auto_unbox = TRUE)
    selector <- sprintf("'[data-testid=\"' + %s + '\"]'", id)
    checkpointWaitForJs(driver, sprintf("document.querySelector(%s) !== null", selector), timeout = checkpointWaitBudget(control = TRUE))
    driver$get_js(sprintf(paste0("(function(){const node=document.querySelector(%s);",
        "const target=node.matches('button,a,input')?node:",
        "(node.querySelector('button,a,input')||node);target.click();return true;})()"), selector))
    checkpointWaitForIdle(driver, timeout = checkpointWaitBudget())
    invisible(NULL)
}

#' Encode a native shinyFiles selection for an existing explicitly supplied file.
#' @param path Existing input file.
#' @param roots Named native file chooser roots.
#' @return The same root/files payload produced by the file chooser.
#' @noRd
checkpointFileSelection <- function(path, roots) {
    path <- normalizePath(path, winslash = "/", mustWork = TRUE)
    if (dir.exists(path)) stop("input must be a file")
    prefixes <- paste0(sub("/+$", "", normalizePath(roots, winslash = "/", mustWork = TRUE)), "/")
    matches <- which(vapply(prefixes, \(prefix) startsWith(path, prefix), logical(1)))
    if (!length(matches)) stop("input is outside the production file chooser roots")
    selected <- matches[which.max(nchar(prefixes[matches]))]
    relative <- substring(path, nchar(prefixes[selected]) + 1L)
    list(root = names(roots)[selected], files = list(strsplit(relative, "/", fixed = TRUE)[[1L]]))
}

#' Drive setup and native file selection through production observers.
#' @param driver Live installed-app browser.
#' @param run_dir Existing unique project parent.
#' @param inputs Explicit report and FASTA file paths.
#' @return NULL invisibly; import completion is not full output publication.
#' @noRd
runCheckpointImportActions <- function(driver, run_dir, inputs) {
    checkpointClickControl(driver, "tile-proteomics")
    checkpointClickControl(driver, "btn-start-analysis")
    driver$set_inputs(experiment_label = "dia_baseline_smoke", project_dir = run_dir)
    checkpointClickControl(driver, "btn-confirm-setup")
    checkpointWaitForJs(driver, "document.querySelector('[data-testid=\"prot-tab-setup\"]') !== null",
        timeout = checkpointWaitBudget())
    roots <- shinyFiles::getVolumes()()
    for (field in c("report", "fasta")) {
        input_name <- if (field == "report") "search_results" else "fasta_file"
        id <- paste0("proteomics_workflow-setup_import-", input_name)
        selection <- checkpointFileSelection(inputs[[field]], roots)
        checkpointWaitForJs(driver, sprintf("document.getElementById(%s) !== null",
            jsonlite::toJSON(id, auto_unbox = TRUE)), timeout = checkpointWaitBudget(control = TRUE))
        driver$get_js(sprintf("Shiny.setInputValue(%s,%s,{priority:'event'});true",
            jsonlite::toJSON(id, auto_unbox = TRUE),
            jsonlite::toJSON(selection, auto_unbox = TRUE)))
        checkpointWaitForIdle(driver, timeout = checkpointWaitBudget())
    }
    args <- list("diann")
    names(args) <- "proteomics_workflow-setup_import-format_override"
    do.call(driver$set_inputs, args)
    checkpointClickControl(driver, "prot-import-process")
    waitCheckpointStepComplete(driver, "Import", timeout = checkpointWaitBudget())
    invisible(NULL)
}

#' Set one real input and wait for its observer.
#' @param driver Browser.
#' @param id Existing input ID.
#' @param value Selected value.
#' @return NULL invisibly.
#' @noRd
checkpointSetInput <- function(driver, id, value) {
    checkpointDriverMemoryBoundary(paste0("enter:input:", id))
    on.exit(checkpointDriverMemoryBoundary(paste0("exit:input:", id)), add = TRUE)
    encoded <- jsonlite::toJSON(id, auto_unbox = TRUE)
    checkpointWaitForJs(driver, sprintf("document.getElementById(%s) !== null", encoded), timeout = checkpointWaitBudget(control = TRUE))
    is_tabs <- driver$get_js(sprintf(
        "document.getElementById(%s).querySelector('a[data-value]') !== null", encoded))
    if (isTRUE(is_tabs)) {
        selected <- driver$get_js(sprintf(paste0("(function(){const el=document.getElementById(%s);",
            "const tab=Array.from(el.querySelectorAll('a[data-value]')).find(x=>x.dataset.value===%s);",
            "if(!tab)return false;tab.click();return true;})()"),
            encoded, jsonlite::toJSON(value, auto_unbox = TRUE)))
        if (!isTRUE(selected)) stop(paste("missing production tab", id, value))
        checkpointWaitForIdle(driver, timeout = checkpointWaitBudget())
        return(invisible(NULL))
    }
    checkpointWaitForJs(driver, sprintf(
        "document.getElementById(%s).classList.contains('shiny-bound-input')", encoded),
        timeout = checkpointWaitBudget(control = TRUE))
    arguments <- stats::setNames(list(value), id)
    # Input-only controls need not render an output until their Apply button.
    # Wait for idle, not an unrelated output change (which can impose 30s delays).
    arguments$wait_ <- FALSE
    do.call(driver$set_inputs, arguments)
    checkpointWaitForIdle(driver, timeout = checkpointWaitBudget())
    invisible(NULL)
}

#' Select only existing options, rejecting ambiguous or invented choices.
#' @param driver Browser.
#' @param id Existing select/selectize input.
#' @param values Requested values.
#' @param multiple Whether this is a multiple selection.
#' @return Actual selected option identities.
#' @noRd
checkpointChooseExisting <- function(driver, id, values, multiple = FALSE) {
    if (!length(values) || anyNA(values) || anyDuplicated(values) ||
        (!isTRUE(multiple) && length(values) != 1L)) {
        stop("invalid or duplicated UI selection")
    }
    encoded_id <- jsonlite::toJSON(id, auto_unbox = TRUE)
    checkpointWaitForJs(driver, sprintf("document.getElementById(%s) !== null", encoded_id), timeout = checkpointWaitBudget(control = TRUE))
    available <- unlist(driver$get_js(sprintf(paste0(
        "(function(){const el=document.getElementById(%s);return el.selectize?",
        "Object.keys(el.selectize.options):Array.from(el.options).map(x=>x.value);})()"), encoded_id)))
    selected <- vapply(values, function(value) {
        exact <- available[available == value]
        if (length(exact) == 1L) return(exact)
        matches <- available[tolower(available) == tolower(value)]
        if (length(matches) != 1L) stop(paste("missing or ambiguous UI option", value))
        matches
    }, character(1), USE.NAMES = FALSE)
    encoded_values <- jsonlite::toJSON(selected, auto_unbox = FALSE)
    driver$get_js(sprintf(paste0("(function(){const el=document.getElementById(%s);",
        "const values=%s;const value=%s?values:values[0];",
        "if(el.selectize){el.selectize.setValue(value,false);}else{",
        "Array.from(el.options).forEach(x=>x.selected=values.includes(x.value));",
        "el.dispatchEvent(new Event('change',{bubbles:true}));}return true;})()"),
        encoded_id, encoded_values, if (multiple) "true" else "false"))
    checkpointWaitForIdle(driver, timeout = checkpointWaitBudget())
    selected
}

#' Drive the pinned two-group diagnostic design through the actual builder.
#' @param driver Browser.
#' @param design_path Explicit pinned sample/group file.
#' @param stop_before_save Diagnostic-only stop after assignments, before assay work.
#' @return Small metadata selection map; not scientific data.
#' @noRd
runCheckpointDesignActions <- function(driver, design_path, stop_before_save = FALSE) {
    design <- utils::read.delim(design_path, stringsAsFactors = FALSE, check.names = FALSE)
    if (!all(c("sample", "group") %in% names(design)) ||
        anyNA(design[c("sample", "group")]) || anyDuplicated(design$sample)) {
        stop("diagnostic design requires unique sample and group metadata")
    }
    groups <- unique(design$group)
    if (length(groups) != 2L) stop("this diagnostic requires an explicit two-group design")
    prefix <- "proteomics_workflow-design_matrix-builder-"
    set <- function(id, value) checkpointSetInput(driver, paste0(prefix, id), value)
    click <- function(id) {
        encoded <- jsonlite::toJSON(paste0(prefix, id), auto_unbox = TRUE)
        checkpointWaitForJs(driver, sprintf("document.getElementById(%s) !== null", encoded), timeout = checkpointWaitBudget(control = TRUE))
        driver$get_js(sprintf("document.getElementById(%s).click();true", encoded))
        # Idle waits include scientific callbacks and their required outputs.
        checkpointWaitForIdle(driver, timeout = checkpointWaitBudget())
    }
    checkpointSetInput(driver, "proteomics_workflow-workflow_tabs", "design")
    set("main_tabs", "Factors")
    for (group in groups) {
        set("new_factor", group)
        click("add_factor")
    }
    set("main_tabs", "Assign Metadata")
    mapping <- list()
    for (group in groups) {
        mapping[[group]] <- checkpointChooseExisting(driver, paste0(prefix, "selected_runs"),
            design$sample[design$group == group], multiple = TRUE)
        checkpointChooseExisting(driver, paste0(prefix, "factor1_select"), group)
        click("assign_metadata")
    }
    if ("replicates" %in% names(design)) {
        if (anyNA(design$replicates) || any(design$replicates < 1)) {
            stop("invalid declared biological replicate numbers")
        }
        set("main_tabs", "Technical Replicates")
        set("tech_rep_assignment_mode", "manual")
        for (group in groups) for (replicate in unique(design$replicates[design$group == group])) {
            samples <- design$sample[design$group == group & design$replicates == replicate]
            if (length(samples) > 1L) {
                driver$get_js(checkpointDesignAssignmentScript(prefix, samples, replicate, prepare = TRUE))
                expected <- sprintf(paste0("Group: %s, Biological Replicate: %s\n",
                    "  Samples: %s\n  Technical Replicates: %s"), group, replicate,
                    paste(samples, collapse = ", "), paste(seq_along(samples), collapse = ", "))
                checkpointWaitForJs(driver, sprintf(paste0("(function(){const el=document.getElementById(%s);",
                    "return !!el && el.textContent.split('\\n\\n').includes(%s) && ",
                    "window.shinytest2 && window.shinytest2.busy!==true;})()"),
                    jsonlite::toJSON(paste0(prefix, "tech_rep_summary"), auto_unbox = TRUE),
                    jsonlite::toJSON(expected, auto_unbox = TRUE)), timeout = checkpointWaitBudget(control = TRUE), stable_ms = 200)
            }
        }
    }
    set("main_tabs", "Contrasts")
    # The pinned diagnostic's KO_vs_WT is second group minus first group.
    checkpointChooseExisting(driver, paste0(prefix, "contrast_group1"), groups[[2L]])
    checkpointChooseExisting(driver, paste0(prefix, "contrast_group2"), groups[[1L]])
    click("add_contrast")
    if (!stop_before_save) {
        click("save_results")
        waitCheckpointStepComplete(driver, "Design", timeout = checkpointWaitBudget())
    }
    mapping
}

#' Flush the selected design inputs before the ordinary production Assign click.
#' @param prefix Existing builder namespace.
#' @param samples,replicate Exact frozen design values.
#' @param prepare Set existing production controls before flushing the assignment.
#' @return Trusted JavaScript; refuses a different or unbound UI selection.
#' @noRd
checkpointDesignAssignmentScript <- function(prefix, samples, replicate, prepare = FALSE) {
    stopifnot(is.character(prefix), length(prefix) == 1L, is.character(samples),
        length(samples) > 1L, !anyNA(samples), !anyDuplicated(samples),
        is.numeric(replicate), length(replicate) == 1L, is.finite(replicate), replicate >= 1)
    sprintf(paste0("(function(){const prefix=%s,samples=%s,replicate=%s;",
        "const selected=document.getElementById(prefix+'tech_rep_samples');",
        "const number=document.getElementById(prefix+'manual_replicate_number');",
        "const mode=document.querySelector('[id=\"'+prefix+'tech_rep_assignment_mode\"] input:checked');",
        "const button=document.getElementById(prefix+'assign_tech_reps');",
        "if(!selected || !selected.selectize || !number || !mode || !button) ",
        "throw new Error('design assignment controls unavailable');",
        "if(%s){if(!samples.every(x=>Object.prototype.hasOwnProperty.call(selected.selectize.options,x))) ",
        "throw new Error('design assignment option absent');",
        "selected.selectize.setValue(samples,true);number.value=String(replicate);}",
        "if(JSON.stringify(selected.selectize.getValue())!==JSON.stringify(samples) || ",
        "Number(number.value)!==replicate || mode.value!=='manual') ",
        "throw new Error('design assignment controls differ from frozen input');",
        "Shiny.setInputValue(prefix+'tech_rep_samples',samples,{priority:'event'});",
        "Shiny.setInputValue(prefix+'tech_rep_assignment_mode','manual',{priority:'event'});",
        "Shiny.setInputValue(prefix+'manual_replicate_number',replicate,{priority:'event'});",
        "button.click();return true;})()"), jsonlite::toJSON(prefix, auto_unbox = TRUE),
        jsonlite::toJSON(samples, auto_unbox = FALSE), jsonlite::toJSON(replicate, auto_unbox = TRUE),
        if (prepare) "true" else "false")
}

#' Reject a mismatched executed design before any QC operation starts.
#' @param design_path Frozen input recipe, including optional replicate assignments.
#' @param saved_path Small design table produced by the real GUI save action.
#' @return Bounded proof of the declared design fields; extra recipe fields are not inferred.
#' @noRd
checkCheckpointSavedDesign <- function(design_path, saved_path) {
    if (!file.exists(saved_path) || file.size(saved_path) > 1048576) stop("saved GUI design is absent or oversized")
    expected <- utils::read.delim(design_path, stringsAsFactors = FALSE, check.names = FALSE)
    actual <- utils::read.delim(saved_path, stringsAsFactors = FALSE, check.names = FALSE)
    stopifnot(all(c("sample", "group") %in% names(expected)))
    fields <- intersect(c("group", "replicates", "tech_rep_group"), names(expected))
    if (!all(c("Run", fields) %in% names(actual)) ||
        !identical(actual$Run, expected$sample)) stop("saved GUI design sample identities/order differ")
    for (field in fields) if (!identical(as.character(actual[[field]]), as.character(expected[[field]]))) {
        stop(paste("saved GUI design differs from frozen input:", field))
    }
    if ("replicates" %in% fields) {
        groups <- interaction(expected$group, expected$replicates, drop = TRUE)
        technical <- stats::ave(seq_len(nrow(expected)), groups, FUN = function(rows) {
            if (length(rows) > 1L) seq_along(rows) else NA_integer_
        })
        if (!"tech_reps" %in% names(actual) ||
            !identical(as.numeric(actual$tech_reps), as.numeric(technical))) {
            stop("saved GUI design differs from frozen input: tech_reps")
        }
        fields <- c(fields, "tech_reps")
    }
    list(status = "passed_before_qc", samples = nrow(actual), checked_fields = c("Run", fields),
        saved_sha256 = digest::digest(file = saved_path, algo = "sha256"))
}

#' Drive the real peptide buttons with the predeclared tiny E2E diagnostic settings.
#' @param driver Browser.
#' @param diagnostic_settings Apply tiny-fixture overrides; FALSE preserves GUI defaults.
#' @return Completed state names; not scientific parity evidence.
#' @noRd
runCheckpointPeptideActions <- function(driver, diagnostic_settings = TRUE, first_qc = FALSE) {
    if (!is.logical(diagnostic_settings) || length(diagnostic_settings) != 1L || is.na(diagnostic_settings)) {
        stop("diagnostic_settings must be TRUE or FALSE")
    }
    checkpointSetInput(driver, "proteomics_workflow-workflow_tabs", "qc")
    prefix <- "proteomics_workflow-quality_control-peptide_qc-"
    checkpointSetInput(driver, "proteomics_workflow-quality_control-qc_tabs_lfq", "Peptide QC")
    settings <- if (diagnostic_settings) list(
        protein_peptide_filtered = list(
            "protein_peptide_filter-min_peptides_per_protein" = 1,
            "protein_peptide_filter-min_peptidoforms_per_protein" = 1),
        intensity_filtered = list("intensity_filter-use_strict_mode" = FALSE,
            "intensity_filter-min_reps_per_group" = 1, "intensity_filter-min_groups" = 1,
            "intensity_filter-intensity_cutoff_percentile" = 0.1),
        sample_filtered = list("sample_filter-min_peptides_per_sample" = 1),
        imputed = list("imputation-proportion_missing_values" = 1)
    ) else list()
    completed <- character()
    actions <- checkpointPeptideActions()
    if (first_qc) actions <- actions[1L]
    for (action in actions) {
        checkpointDriverMemoryBoundary(paste0("peptide:", action$state, ":begin"))
        checkpointSetInput(driver, paste0(prefix, "peptide_filter_tabs"), action$tab)
        for (field in names(settings[[action$state]])) {
            checkpointSetInput(driver, paste0(prefix, field), settings[[action$state]][[field]])
        }
        id <- jsonlite::toJSON(paste0(prefix, action$action), auto_unbox = TRUE)
        checkpointWaitForJs(driver, sprintf("document.getElementById(%s) !== null", id), timeout = checkpointWaitBudget(control = TRUE))
        driver$get_js(sprintf("document.getElementById(%s).click();true", id))
        output <- jsonlite::toJSON(paste0(prefix, action$result), auto_unbox = TRUE)
        state <- jsonlite::toJSON(action$state, auto_unbox = TRUE)
        checkpointWaitForJs(driver, sprintf(paste0("(function(){const el=document.getElementById(%s);",
            "return (el && el.textContent.includes(%s)) || ",
            "document.querySelector('.shiny-notification-error') !== null;})()"), output, state),
            timeout = checkpointWaitBudget())
        error <- driver$get_js(paste0("(function(){const el=document.querySelector(",
            "'.shiny-notification-error');return el?el.textContent:null;})()"))
        if (!is.null(error)) stop(paste("production peptide action failed:", action$state, error))
        completed <- c(completed, action$state)
        checkpointDriverMemoryBoundary(paste0("peptide:", action$state, ":complete"))
    }
    completed
}

#' Drive original protein QC actions using explicit tiny-diagnostic settings.
#' @param driver Installed production browser.
#' @param diagnostic_settings Apply tiny-fixture overrides; FALSE preserves GUI defaults.
#' @return Completed state names, not a scientific gate.
#' @noRd
runCheckpointProteinActions <- function(driver, diagnostic_settings = TRUE) {
    if (!is.logical(diagnostic_settings) || length(diagnostic_settings) != 1L || is.na(diagnostic_settings)) {
        stop("diagnostic_settings must be TRUE or FALSE")
    }
    prefix <- "proteomics_workflow-quality_control-protein_qc-"
    checkpointSetInput(driver, "proteomics_workflow-quality_control-qc_tabs_lfq", "Protein QC")
    actions <- list(
        list(tab = "IQ Protein Rollup", action = "rollup-apply_iq_rollup",
            result = "rollup-iq_rollup_results", state = "protein_s4_created"),
        list(tab = "Accession Cleanup", action = "cleanup-apply_accession_cleanup",
            result = "cleanup-accession_cleanup_results", state = "protein_accession_cleaned"),
        list(tab = "Protein Intensity Filter", action = "intensity_filter-apply_protein_intensity_filter",
            result = "intensity_filter-protein_intensity_filter_results", state = "protein_intensity_filtered"),
        list(tab = "Duplicate Removal", action = "duplicate_removal-apply_duplicate_removal",
            result = "duplicate_removal-duplicate_removal_results", state = "duplicates_removed"),
        list(tab = "Protein Replicate Filter", action = "replicate_filter-apply_protein_replicate_filter",
            result = "replicate_filter-protein_replicate_filter_results", state = "protein_replicate_filtered")
    )
    settings <- if (diagnostic_settings) list(
        protein_s4_created = list("rollup-rollup_method" = "iq"),
        protein_intensity_filtered = list("intensity_filter-use_strict_mode" = FALSE,
            "intensity_filter-min_reps_per_group" = 1, "intensity_filter-min_groups" = 1,
            "intensity_filter-proteins_intensity_cutoff_percentile" = 0.1),
        protein_replicate_filtered = list("replicate_filter-parallel_cores" = 1)
    ) else list(protein_s4_created = list("rollup-rollup_method" = "iq"))
    completed <- character()
    for (action in actions) {
        checkpointSetInput(driver, paste0(prefix, "protein_filter_tabs"), action$tab)
        for (field in names(settings[[action$state]])) {
            checkpointSetInput(driver, paste0(prefix, field), settings[[action$state]][[field]])
        }
        id <- jsonlite::toJSON(paste0(prefix, action$action), auto_unbox = TRUE)
        checkpointWaitForJs(driver, sprintf("document.getElementById(%s) !== null", id), timeout = checkpointWaitBudget(control = TRUE))
        driver$get_js(sprintf("document.getElementById(%s).click();true", id))
        output <- jsonlite::toJSON(paste0(prefix, action$result), auto_unbox = TRUE)
        state <- jsonlite::toJSON(action$state, auto_unbox = TRUE)
        checkpointWaitForJs(driver, sprintf(paste0("(function(){const el=document.getElementById(%s);",
            "return (el && el.textContent.includes(%s)) || ",
            "document.querySelector('.shiny-notification-error') !== null;})()"), output, state),
            timeout = checkpointWaitBudget())
        error <- driver$get_js(paste0("(function(){const el=document.querySelector(",
            "'.shiny-notification-error');return el?el.textContent:null;})()"))
        if (!is.null(error)) stop(paste("production protein action failed:", action$state, error))
        completed <- c(completed, action$state)
    }
    waitCheckpointStepComplete(driver, "QC", timeout = checkpointWaitBudget())
    completed
}

#' Click a production action and wait for its actual UI result or error.
#' @param driver Browser.
#' @param input_id Fully namespaced action id.
#' @param selector CSS selector for completion text.
#' @param text Required completion text.
#' @param timeout Maximum action wait in milliseconds.
#' @return TRUE invisibly; a notification error is never completion.
#' @noRd
checkpointClickAndWaitText <- function(driver, input_id, selector, text, timeout = NULL) {
    checkpointDriverMemoryBoundary(paste0("enter:action:", input_id))
    on.exit(checkpointDriverMemoryBoundary(paste0("exit:action:", input_id)), add = TRUE)
    id <- jsonlite::toJSON(input_id, auto_unbox = TRUE)
    selector <- jsonlite::toJSON(selector, auto_unbox = TRUE)
    text <- jsonlite::toJSON(text, auto_unbox = TRUE)
    checkpointWaitForJs(driver, sprintf("document.getElementById(%s) !== null", id), timeout = checkpointWaitBudget(control = TRUE))
    driver$get_js(sprintf("document.getElementById(%s).click();true", id))
    timeout <- if (is.null(timeout)) checkpointWaitBudget() else min(timeout, checkpointWaitBudget())
    checkpointWaitForJs(driver, sprintf(paste0("(function(){return Array.from(document.querySelectorAll(%s))",
        ".some(el=>el.textContent.includes(%s)) || ",
        "document.querySelector('.shiny-notification-error') !== null;})()"), selector, text),
        timeout = timeout)
    error <- driver$get_js(paste0("(function(){const el=document.querySelector(",
        "'.shiny-notification-error');return el?el.textContent:null;})()"))
    if (!is.null(error)) stop(paste("production action failed:", input_id, error))
    invisible(TRUE)
}

#' Drive frozen normalization settings, with an explicit real RUV/Pearson diagnostic.
#' @param driver Installed production browser.
#' @param profile NULL keeps the historical Skip route; ruv_pearson runs GUI defaults.
#' @return Bounded settings and observed action names, not whole-workflow acceptance.
#' @noRd
runCheckpointNormalizationActions <- function(driver, profile = NULL) {
    stress <- identical(profile, "ruv_pearson")
    if (!is.null(profile) && !stress) stop("invalid normalization diagnostic profile")
    prefix <- "proteomics_workflow-normalization-"
    checkpointSetInput(driver, "proteomics_workflow-workflow_tabs", "normalization")
    checkpointSetInput(driver, paste0(prefix, "norm_method"), "cyclicloess")
    checkpointSetInput(driver, paste0(prefix, "ruv_mode"), if (stress) "automatic" else "skip")
    if (stress) {
        values <- list(auto_percentage_min = 1, auto_percentage_max = 20,
            separation_metric = "max_difference", k_penalty_weight = 0.5,
            adaptive_k_penalty = TRUE, ruv_grouping_variable = "group")
        for (name in names(values)) {
            checkpointSetInput(driver, paste0(prefix, name), values[[name]])
        }
    }
    checkpointClickAndWaitText(driver, paste0(prefix, "run_normalization"),
        ".shiny-notification", if (stress) "Normalization and RUV correction completed!" else "Normalization completed (RUV skipped)!")
    checkpointSetInput(driver, paste0(prefix, "norm_qc_tabs"), "Correlation Filtering")
    if (stress) checkpointSetInput(driver, paste0(prefix, "min_pearson_correlation_threshold"), 0.5)
    checkpointClickAndWaitText(driver, paste0(prefix, if (stress) "apply_correlation_filter" else "skip_correlation_filter"),
        paste0("#", prefix, "correlation_filter_summary"),
        if (stress) "Correlation filtering completed successfully!" else "Correlation filtering SKIPPED.")
    checkpointClickAndWaitText(driver, paste0(prefix, "export_filtered_session"),
        ".shiny-notification", "Filtered session data exported successfully!")
    result <- list(method = "cyclicloess", ruv = if (stress) "automatic" else "skip",
        correlation = if (stress) "apply" else "skip",
        actions = c("run_normalization", if (stress) "apply_correlation_filter" else "skip_correlation_filter", "export_filtered_session"))
    if (stress) result$correlation_threshold <- 0.5
    result
}

#' Verify the shipped Glimma widget renders and table selection updates expression.
#' @param driver Installed browser.
#' @return Bounded point count; no scientific objects.
#' @noRd
runCheckpointInteractiveVolcano <- function(driver) {
    prefix <- "proteomics_workflow-differential_expression-"
    checkpointSetInput(driver, paste0(prefix, "volcano_interactive"), TRUE)
    id <- jsonlite::toJSON(paste0(prefix, "volcano_glimma"), auto_unbox = TRUE)
    checkpointWaitForJs(driver, sprintf(paste0("(function(){const el=document.getElementById(%s);",
        "return el && (el.classList.contains('shiny-output-error') || ",
        "(el.querySelector('.xyContainer svg .mark-symbol.marks path') && ",
        "el.querySelector('.controlContainer table tbody tr')));})()"), id), timeout = checkpointWaitBudget())
    error <- driver$get_js(sprintf(paste0("(function(){const el=document.getElementById(%s);",
        "return el.classList.contains('shiny-output-error')?el.textContent:null;})()"), id))
    if (!is.null(error)) stop(paste("production interactive volcano failed:", error))
    points <- driver$get_js(sprintf(
        "document.getElementById(%s).querySelectorAll('.xyContainer svg .mark-symbol.marks path').length", id))
    driver$get_js(sprintf(
        "document.getElementById(%s).querySelector('.controlContainer table tbody tr').click();true", id))
    checkpointWaitForJs(driver, sprintf(paste0("(function(){const el=document.getElementById(%s);",
        "return el.querySelector('table tbody tr.selected') !== null && ",
        "el.querySelector('.expressionContainer svg .mark-symbol.marks path') !== null;})()"), id),
        timeout = checkpointWaitBudget())
    expression_points <- driver$get_js(sprintf(
        "document.getElementById(%s).querySelectorAll('.expressionContainer svg .mark-symbol.marks path').length", id))
    list(rendered_points = points, expression_points = expression_points,
        table_selection_updates_expression = TRUE)
}

#' Exercise the default heatmap view and its existing save action.
#' @param driver Installed production browser.
#' @return Bounded observed settings; files require separate verification.
#' @noRd
runCheckpointHeatmapActions <- function(driver) {
    prefix <- "proteomics_workflow-differential_expression-"
    checkpointSetInput(driver, paste0(prefix, "da_results_tabs"), "Heatmap")
    id <- jsonlite::toJSON(paste0(prefix, "heatmap_plot"), auto_unbox = TRUE)
    checkpointWaitForJs(driver, sprintf(paste0("(function(){const el=document.getElementById(%s);",
        "if(!el)return false;const img=el.querySelector('img');return ",
        "(img && img.complete && img.naturalWidth>1) || ",
        "el.classList.contains('shiny-output-error');})()"), id), timeout = checkpointWaitBudget())
    error <- driver$get_js(sprintf(paste0("(function(){const el=document.getElementById(%s);",
        "return el.classList.contains('shiny-output-error')?el.textContent:null;})()"), id))
    if (!is.null(error)) stop(paste("production heatmap failed:", error))
    fields <- c("heatmap_contrast", "heatmap_top_n", "heatmap_clustering", "heatmap_scaling",
        "heatmap_cluster_method", "heatmap_distance_method", "heatmap_tree_cut_method",
        "heatmap_n_clusters", "heatmap_color_scheme")
    settings <- driver$get_js(sprintf(paste0("Object.fromEntries(%s.map(name=>",
        "[name,document.getElementById(%s+name).value]))"), jsonlite::toJSON(fields),
        jsonlite::toJSON(prefix, auto_unbox = TRUE)))
    checkpointClickAndWaitText(driver, paste0(prefix, "save_heatmap"),
        ".shiny-notification", "Heatmap and cluster info saved to publication_graphs/Heatmap")
    settings
}

#' Exercise production session restore, DA analysis and the all-results table.
#' @param driver Installed production browser.
#' @return Bounded observed actions/settings, not independent scientific acceptance.
#' @noRd
runCheckpointDaActions <- function(driver) {
    checkpointDriverMemoryBoundary("enter:DA_session_load_and_views")
    on.exit(checkpointDriverMemoryBoundary("exit:DA_session_load_and_views"), add = TRUE)
    prefix <- "proteomics_workflow-differential_expression-"
    checkpointSetInput(driver, "proteomics_workflow-workflow_tabs", "da")
    checkpointClickAndWaitText(driver, paste0(prefix, "load_filtered_session"),
        ".shiny-notification", "Session loaded successfully!")
    checkpointSetInput(driver, paste0(prefix, "formula_string"), "~ 0 + group")
    checkpointSetInput(driver, paste0(prefix, "da_q_val_thresh"), 0.05)
    checkpointSetInput(driver, paste0(prefix, "treat_lfc_cutoff"), 0)
    checkpointClickAndWaitText(driver, paste0(prefix, "run_da_analysis"),
        paste0("#", prefix, "da_status"), "[OK] Analysis Complete")
    checkpointSetInput(driver, paste0(prefix, "da_results_tabs"), "DA Results Table")
    checkpointSetInput(driver, paste0(prefix, "table_significance"), "all")
    selector <- jsonlite::toJSON(paste0("#", prefix, "da_results_table"), auto_unbox = TRUE)
    checkpointWaitForJs(driver, sprintf(paste0("(function(){const el=document.querySelector(%s);",
        "return el && ((el.querySelector('tbody td') !== null && ",
        "!el.querySelector('.dataTables_empty,.dt-empty')) || ",
        "el.classList.contains('shiny-output-error'));})()"), selector), timeout = checkpointWaitBudget())
    error <- driver$get_js(sprintf(paste0("(function(){const el=document.querySelector(%s);",
        "const headers=Array.from(el.querySelectorAll('thead th')).map(x=>x.textContent.trim());",
        "return el.classList.contains('shiny-output-error') || ",
        "!['log2FC','raw_pvalue','fdr_qvalue'].every(x=>headers.includes(x))",
        "?el.textContent.slice(0,2000):null;})()"), selector))
    if (!is.null(error)) stop(paste("production DA table failed:", error))
    checkpointSetInput(driver, paste0(prefix, "da_results_tabs"), "Volcano Plot")
    checkpointSetInput(driver, paste0(prefix, "volcano_interactive"), FALSE)
    plot_id <- jsonlite::toJSON(paste0(prefix, "volcano_plot_static"), auto_unbox = TRUE)
    checkpointWaitForJs(driver, sprintf(paste0("(function(){const el=document.getElementById(%s);",
        "if(!el)return false;const img=el.querySelector('img');return ",
        "(img && img.complete && img.naturalWidth>0) || ",
        "el.classList.contains('shiny-output-error');})()"), plot_id), timeout = checkpointWaitBudget())
    error <- driver$get_js(sprintf(paste0("(function(){const el=document.getElementById(%s);",
        "return el.classList.contains('shiny-output-error')?el.textContent:null;})()"), plot_id))
    if (!is.null(error)) stop(paste("production static volcano failed:", error))
    interactive <- runCheckpointInteractiveVolcano(driver)
    heatmap <- runCheckpointHeatmapActions(driver)
    list(formula = "~ 0 + group", q_value = 0.05, treat_lfc = 0,
        interactive = interactive, heatmap = heatmap,
        actions = c("load_filtered_session", "run_da_analysis", "all_results_table", "static_volcano",
            "interactive_volcano", "heatmap", "save_heatmap"))
}

#' Exercise existing summary publication and session-state actions after DA.
#' @param driver Installed production browser.
#' @return Bounded action names and selected label; files need separate verification.
#' @noRd
runCheckpointSummaryActions <- function(driver) {
    prefix <- "proteomics_workflow-session_summary-"
    checkpointSetInput(driver, "proteomics_workflow-workflow_tabs", "session_summary")
    checkpointSetInput(driver, paste0(prefix, "experiment_label"), "DIA Experiment")
    checkpointClickAndWaitText(driver, paste0(prefix, "save_workflow_args"),
        ".shiny-notification", "Study parameters saved successfully")
    checkpointClickAndWaitText(driver, paste0(prefix, "copy_to_publication"),
        paste0("#", prefix, "copy_status"), "Files copied to publication directory successfully")
    checkpointClickAndWaitText(driver, paste0(prefix, "export_session_state"),
        ".shiny-notification", "Session state exported to:")
    list(experiment_label = "DIA Experiment",
        actions = c("save_workflow_args", "copy_to_publication", "export_session_state"))
}

#' Measure unique live R heap cells after a standardized collection.
#' @return Total bytes represented by used Ncells and Vcells.
#' @noRd
checkpointUtilityHeapBytes <- function() {
    cells <- gc(full = TRUE)
    sizes <- c(Ncells = 56, Vcells = 8)
    used <- cells[names(sizes), "used"]
    as.numeric(sum(used * sizes))
}

#' Measure one utility operation without changing its return value.
#' @param metrics Mutable member-local metrics environment.
#' @param label Unique operation label.
#' @param operation Zero-argument operation.
#' @return Operation result.
#' @noRd
checkpointUtilityTimed <- function(metrics, label, operation) {
    if (is.function(metrics$observe)) metrics$observe(paste0("before:", label))
    start <- proc.time()
    value <- operation()
    finish <- proc.time()
    event <- list(label = label,
        elapsed_seconds = unname((finish - start)[["elapsed"]]),
        cpu_seconds = unname(sum((finish - start)[c("user.self", "sys.self")])))
    metrics$events[[length(metrics$events) + 1L]] <- event
    if (is.function(metrics$observe)) metrics$observe(paste0("after:", label))
    value
}

#' Append a diagnostic-only phase snapshot without collecting or retaining payloads.
#' @param path Request-owned JSONL evidence file.
#' @param label Phase label.
#' @param memory_snapshot Existing service-counter callback.
#' @param heap_bytes Live R heap from an already-required collection, otherwise NULL.
#' @return NULL invisibly. Cgroup components and R heap are not additive.
#' @noRd
observeCheckpointUtilityMemory <- function(path, label, memory_snapshot, heap_bytes = NULL) {
    record <- list(label = label, timestamp = as.numeric(Sys.time()),
        r_live_heap_bytes = heap_bytes, service = memory_snapshot())
    cat(jsonlite::toJSON(record, auto_unbox = TRUE, null = "null", digits = NA),
        "\n", file = path, append = TRUE, sep = "")
    invisible(NULL)
}

#' Build one deterministic legacy state record.
#' @param object Protein S4 state.
#' @param recipe Frozen recipe.
#' @param state_name,previous_state Topology names.
#' @param generation Deterministic generation number.
#' @return Historical record compatible with the legacy states list.
#' @noRd
checkpointUtilityStateRecord <- function(object, recipe, state_name,
                                         previous_state, generation) {
    list(timestamp = as.POSIXct(recipe$seed + generation, origin = "1970-01-01", tz = "UTC"),
        data = object,
        config = list(route = "iq", recipe = recipe$id, generation = generation),
        description = paste("synthetic DIA history", state_name),
        previous_state = previous_state,
        s4_class = class(object)[[1L]],
        audit_metadata = list(provenance_status = "fixture_not_scientific_estimator"))
}

#' Prepare immutable utility inputs outside the measured history process.
#' @param library_path Installed package library used by both treatments.
#' @param run_dir Fresh guarded preparation directory.
#' @param contract_path Frozen recipe contract.
#' @param recipe_id Frozen timing recipe, excluding capacity-only workloads.
#' @return Input manifest; no history is retained in this preparation process.
#' @noRd
prepareCheckpointUtilityInputs <- function(library_path, run_dir, contract_path, recipe_id) {
    recordedCheckpointApi(library_path)
    source_root <- normalizePath(file.path(dirname(contract_path), "..", "..", ".."))
    helpers <- new.env(parent = baseenv())
    sys.source(file.path(source_root, "tests/testthat/helper-workflow-checkpoint-fixtures.R"),
        envir = helpers)
    contract <- jsonlite::read_json(contract_path)
    matches <- Filter(\(recipe) identical(recipe$id, recipe_id), contract$recipes)
    if (length(matches) != 1L || identical(recipe_id, "capacity")) stop("invalid utility recipe")
    recipe <- helpers$validateWorkflowCheckpointRecipe(matches[[1L]], TRUE, run_dir, source_root)
    cycles <- if (is.null(recipe$cycles)) 1L else as.integer(recipe$cycles)
    inputs <- vector("list", (recipe$history_depth + 2L) * cycles + 1L)
    for (generation in seq_along(inputs) - 1L) {
        object <- helpers$makeWorkflowCheckpointProtein(recipe, generation,
            allow_large = TRUE, output_root = run_dir, repo_root = source_root)
        name <- sprintf("generation-%04d.rds", generation)
        path <- file.path(run_dir, name)
        if (file.exists(path)) stop("refusing to replace utility input")
        saveRDS(object, path)
        object <- NULL
        inputs[[generation + 1L]] <- list(path = name, bytes = file.size(path),
            sha256 = digest::digest(file = path, algo = "sha256"))
    }
    package_path <- file.path(normalizePath(library_path), "MultiScholaR")
    artifacts <- c("NAMESPACE", "R/MultiScholaR.rdb", "R/MultiScholaR.rdx")
    result <- list(status = "prepared_utility_inputs", recipe = recipe,
        contract_sha256 = digest::digest(file = contract_path, algo = "sha256"),
        fixture_helper_sha256 = digest::digest(file = file.path(source_root,
            "tests/testthat/helper-workflow-checkpoint-fixtures.R"), algo = "sha256"),
        package_sha256 = as.list(stats::setNames(vapply(file.path(package_path, artifacts),
            digest::digest, character(1), file = TRUE, algo = "sha256"), artifacts)),
        inputs = inputs, fixture_generation_timed = FALSE, promotion_authorized = FALSE)
    jsonlite::write_json(result, file.path(run_dir, "utility-inputs.json"),
        auto_unbox = TRUE, pretty = TRUE)
    result
}

#' Read one integrity-checked prepared input, never a spare full history.
#' @param directory Prepared input directory.
#' @param inputs Bound immutable manifest entries.
#' @param generation Zero-based deterministic fixture generation.
#' @return One complete S4 input; shared variants reuse generation zero in memory.
#' @noRd
readCheckpointUtilityInput <- function(directory, inputs, generation) {
    entry <- inputs[[generation + 1L]]
    expected <- sprintf("generation-%04d.rds", generation)
    if (!identical(entry$path, expected)) stop("utility input generation binding differs")
    path <- file.path(directory, expected)
    if (!file.exists(path) || !identical(as.numeric(file.size(path)), as.numeric(entry$bytes)) ||
        !identical(digest::digest(file = path, algo = "sha256"), entry$sha256)) {
        stop("prepared utility input integrity differs")
    }
    readRDS(path)
}

#' Construct one independently changed or metadata-sharing protein state.
#' @param helpers Frozen DIA fixture helper environment.
#' @param recipe Frozen recipe.
#' @param generation State generation.
#' @param variant independent, shared_metadata or retained_alias.
#' @param template Optional prior state for sharing metadata.
#' @param run_dir External disk-backed member directory.
#' @return Valid ProteinQuantitativeData.
#' @noRd
checkpointUtilityProteinState <- function(helpers, recipe, generation, variant,
                                          template, run_dir) {
    if (identical(variant, "shared_metadata") && !is.null(template)) {
        object <- template
        quant <- methods::slot(object, "protein_quant_table")
        sample <- methods::slot(object, "design_matrix")[[methods::slot(object, "sample_id")]][[1L]]
        quant[[sample]][[1L]] <- quant[[sample]][[1L]] + generation / 1000
        methods::slot(object, "protein_quant_table") <- quant
        args <- methods::slot(object, "args")
        args$fixture$generation <- generation
        methods::slot(object, "args") <- args
    } else if (is.function(helpers$readPreparedState)) {
        object <- helpers$readPreparedState(generation)
    } else {
        object <- helpers$makeWorkflowCheckpointProtein(recipe, generation,
            allow_large = TRUE, output_root = run_dir, repo_root = helpers$source_root)
    }
    methods::validObject(object)
    object
}

#' Save or retain one historical record under the selected treatment.
#' @param context Mutable member context.
#' @param record Complete portable legacy record.
#' @param state_name Logical state name.
#' @param generation Integer generation.
#' @return Accepted descriptor for saved treatment, otherwise the record.
#' @noRd
checkpointUtilityAcceptState <- function(context, record, state_name, generation) {
    if (identical(context$treatment, "resident")) return(record)
    parent <- if (is.null(context$last_reference)) NULL else list(
        parent_id = context$last_reference$entity_id,
        revision = as.integer(generation))
    identity <- list(input = paste0(context$recipe$id, "-", context$variant),
        design = paste0("samples-", context$recipe$n_samples),
        parameters = paste0("generation-", generation),
        selectors = state_name,
        software = context$software)
    metadata <- list(class = class(record), shape = integer(), factors = list())
    reference <- checkpointUtilityTimed(context$metrics, paste0("save:", state_name), \() {
        context$api("saveWorkflowCheckpoint")(context$store, record,
            sprintf("g%04d", generation), identity, metadata, lineage = parent,
            codec = if (identical(context$treatment, "qs2")) "qs2" else "rds")
    })
    context$metrics$max_save_seconds <- max(context$metrics$max_save_seconds,
        context$metrics$events[[length(context$metrics$events)]]$elapsed_seconds)
    context$last_reference <- reference
    reference
}

#' Resolve one logical state from resident memory or the durable store.
#' @param context Mutable member context.
#' @param state_name Logical state name.
#' @param cold Whether to include this read in cold-revert latency.
#' @return Complete legacy state record.
#' @noRd
checkpointUtilityResolveState <- function(context, state_name, cold = FALSE) {
    value <- context$states[[state_name]]
    if (identical(context$treatment, "resident")) return(value)
    record <- checkpointUtilityTimed(context$metrics, paste0("hydrate:", state_name), \() {
        context$api("hydrateSavedWorkflowCheckpoint")(context$store, value)
    })
    event <- context$metrics$events[[length(context$metrics$events)]]
    context$metrics$payload_read_bytes <- context$metrics$payload_read_bytes + value$payload$bytes
    if (isTRUE(cold)) {
        context$metrics$max_cold_revert_seconds <- max(
            context$metrics$max_cold_revert_seconds, event$elapsed_seconds)
    }
    record
}

#' Construct the complete self-contained historical export envelope.
#' @param context Mutable utility member context.
#' @param active_history Ordered active lineage.
#' @param current_state Current logical state name.
#' @return Complete legacy envelope and hydrated records.
#' @noRd
checkpointUtilityExportEnvelope <- function(context, active_history, current_state) {
    records <- lapply(names(context$states), \(name) {
        checkpointUtilityResolveState(context, name, cold = FALSE)
    })
    names(records) <- names(context$states)
    current <- records[[current_state]]$data
    list(envelope = list(r6_complete_states = records,
        r6_state_history = as.list(active_history), current_state = current_state,
        current_s4_object = current, correlation_filtered_s4 = current,
        final_protein_count = nrow(methods::slot(current, "protein_quant_table")),
        final_sample_count = nrow(methods::slot(current, "design_matrix"))),
        records = records)
}

#' Execute one cgroup-owned resident or saved-history member.
#' @param library_path Installed candidate library.
#' @param run_dir Fresh request-owned external directory.
#' @param contract_path Frozen recipe contract.
#' @param recipe_id,treatment,variant,block Explicit member identity.
#' @param input_dir Separately prepared inputs; NULL is diagnostic fixture generation only.
#' @param diagnostics Enable diagnostic-only phase/allocation profiles.
#' @param memory_snapshot Existing service-counter callback, required for diagnostics.
#' @return Completed receipt; full history is exported with the original serializer.
#' @noRd
executeCheckpointUtilityMember <- function(library_path, run_dir, contract_path,
    recipe_id, treatment, variant, block, input_dir = NULL, diagnostics = FALSE,
    memory_snapshot = NULL) {
    stopifnot(treatment %in% c("resident", "saved", "qs2"),
        variant %in% c("independent", "shared_metadata", "retained_alias"),
        is.numeric(block), length(block) == 1L, block >= 1, block == floor(block))
    library_path <- normalizePath(library_path, mustWork = TRUE)
    run_dir <- normalizePath(run_dir, mustWork = TRUE)
    contract_path <- normalizePath(contract_path, mustWork = TRUE)
    stopifnot(is.logical(diagnostics), length(diagnostics) == 1L, !is.na(diagnostics))
    complete <- FALSE
    observe <- function(label, heap_bytes = NULL) {
        if (diagnostics) observeCheckpointUtilityMemory(
            file.path(run_dir, "memory-phases.jsonl"), label, memory_snapshot, heap_bytes)
        invisible(NULL)
    }
    if (diagnostics) {
        stopifnot(is.function(memory_snapshot), isTRUE(capabilities("profmem")))
        observe("before_namespace")
        utils::Rprof(file.path(run_dir, "utility.Rprof"), interval = 0.02,
            memory.profiling = TRUE, gc.profiling = TRUE)
        utils::Rprofmem(file.path(run_dir, "utility.Rprofmem"), threshold = 1000000)
        on.exit({
            utils::Rprofmem(NULL)
            utils::Rprof(NULL)
            observe(if (complete) "complete_before_teardown" else "failure_before_teardown")
        }, add = TRUE)
    }
    .libPaths(unique(c(library_path, .libPaths())))
    namespace <- loadNamespace("MultiScholaR", lib.loc = library_path)
    observe("after_namespace")
    package_path <- normalizePath(getNamespaceInfo(namespace, "path"), mustWork = TRUE)
    if (!identical(package_path, normalizePath(file.path(library_path, "MultiScholaR")))) {
        stop("utility member loaded a different installed package")
    }
    source_root <- normalizePath(file.path(dirname(contract_path), "..", "..", ".."))
    helpers <- new.env(parent = baseenv())
    sys.source(file.path(source_root, "tests", "testthat",
        "helper-workflow-checkpoint-fixtures.R"), envir = helpers)
    helpers$source_root <- source_root
    contract <- jsonlite::read_json(contract_path, simplifyVector = FALSE)
    matches <- Filter(\(recipe) identical(recipe$id, recipe_id), contract$recipes)
    if (length(matches) != 1L || identical(recipe_id, "capacity")) stop("invalid utility recipe")
    recipe <- helpers$validateWorkflowCheckpointRecipe(matches[[1L]], TRUE, run_dir, source_root)
    prepared <- NULL
    if (!is.null(input_dir)) {
        input_dir <- normalizePath(input_dir, mustWork = TRUE)
        prepared <- jsonlite::read_json(file.path(input_dir, "utility-inputs.json"))
        if (!identical(prepared$status, "prepared_utility_inputs") ||
            !identical(prepared$recipe, recipe) || !identical(prepared$contract_sha256,
                digest::digest(file = contract_path, algo = "sha256")) ||
            !identical(prepared$fixture_helper_sha256, digest::digest(file = file.path(source_root,
                "tests/testthat/helper-workflow-checkpoint-fixtures.R"), algo = "sha256")) ||
            !setequal(names(prepared$package_sha256),
                c("NAMESPACE", "R/MultiScholaR.rdb", "R/MultiScholaR.rdx"))) {
            stop("prepared utility recipe binding differs")
        }
        for (name in names(prepared$package_sha256)) {
            if (!identical(digest::digest(file = file.path(package_path, name), algo = "sha256"),
                prepared$package_sha256[[name]])) stop("prepared utility build binding differs")
        }
        helpers$readPreparedState <- function(generation) {
            readCheckpointUtilityInput(input_dir, prepared$inputs, generation)
        }
    }
    api <- \(name) get(name, envir = namespace, inherits = FALSE)
    metrics <- new.env(parent = emptyenv())
    metrics$events <- list()
    metrics$max_save_seconds <- 0
    metrics$max_cold_revert_seconds <- 0
    metrics$payload_read_bytes <- 0
    metrics$observe <- if (diagnostics) observe else NULL
    software <- paste0("MultiScholaR-", as.character(utils::packageVersion("MultiScholaR")))
    if (identical(treatment, "qs2") &&
        !identical(as.character(utils::packageVersion("qs2")), "0.2.2")) {
        stop("utility qs2 qualification requires version 0.2.2")
    }
    store <- if (!identical(treatment, "resident")) api("newWorkflowEntityStore")(
        run_dir, paste(recipe$id, variant, block, sep = "-"), "DIA") else NULL
    on.exit(if (!is.null(store)) api("closeWorkflowEntityStore")(store), add = TRUE)
    context <- new.env(parent = emptyenv())
    list2env(list(treatment = treatment, variant = variant, recipe = recipe,
        software = software, api = api, metrics = metrics, store = store,
        states = list(), last_reference = NULL), envir = context)
    receipt_path <- file.path(run_dir, "member-receipt.json")
    receipt <- list(status = "running", recipe = recipe$id, treatment = treatment,
        variant = variant, block = as.integer(block), package_path = package_path,
        process_id = Sys.getpid(), contract_sha256 = digest::digest(file = contract_path,
            algo = "sha256"), promotion_authorized = FALSE)
    receipt$codec <- if (identical(treatment, "resident")) NULL else
        if (identical(treatment, "qs2")) "qs2" else "rds"
    receipt$codec_settings <- if (identical(treatment, "qs2")) list(
        version = "0.2.2", compress_level = 3L, shuffle = TRUE,
        nthreads = 2L, validate_checksum = TRUE) else NULL
    source_files <- c("tools/profiling/workflow_checkpoint_dia_driver.R",
        "tools/profiling/run_workflow_checkpoint_proof.R",
        "tools/profiling/run_workflow_checkpoint_proof.sh",
        "tests/testthat/helper-workflow-checkpoint-fixtures.R")
    receipt$source_sha256 <- as.list(stats::setNames(vapply(
        file.path(source_root, source_files), digest::digest, character(1),
        file = TRUE, algo = "sha256"), source_files))
    receipt$fixture_generation_in_timed_service <- is.null(input_dir)
    receipt$prepared_inputs_sha256 <- if (is.null(input_dir)) NULL else
        digest::digest(file = file.path(input_dir, "utility-inputs.json"), algo = "sha256")
    receipt$performance_qualified <- FALSE
    receipt$diagnostics <- diagnostics
    jsonlite::write_json(receipt, receipt_path, auto_unbox = TRUE, pretty = TRUE)
    started <- proc.time()
    baseline_heap <- checkpointUtilityHeapBytes()
    observe("baseline_after_gc", baseline_heap)
    cycles <- if (is.null(recipe$cycles)) 1L else as.integer(recipe$cycles)
    cycle_metrics <- vector("list", cycles)
    export_path <- file.path(run_dir, "legacy-session.rds")
    for (cycle in seq_len(cycles)) {
        generation_offset <- (cycle - 1L) * (recipe$history_depth + 2L)
        template <- if (identical(variant, "shared_metadata"))
            checkpointUtilityProteinState(helpers, recipe, 0L, variant, NULL, run_dir) else NULL
        aliases <- list()
        previous <- NA_character_
        for (state_index in seq_len(recipe$history_depth)) {
            generation <- generation_offset + state_index
            state_name <- sprintf("state_%03d", state_index)
            object <- checkpointUtilityProteinState(helpers, recipe, generation,
                variant, template, run_dir)
            record <- checkpointUtilityStateRecord(object, recipe, state_name,
                previous, generation)
            record <- api("prepareDIAWorkflowHistory")(stats::setNames(list(record), state_name))[[1L]]
            context$states[[state_name]] <- checkpointUtilityAcceptState(context,
                record, state_name, generation)
            if (!identical(treatment, "resident") && identical(variant, "retained_alias")) {
                aliases[[state_name]] <- record
            }
            previous <- state_name
        }
        current_state <- previous
        current_record <- record
        object <- record <- template <- NULL
        observe("forward_before_gc")
        retained_total <- checkpointUtilityHeapBytes()
        observe("forward_after_gc", retained_total)
        retained_heap <- max(1, retained_total - baseline_heap)
        active_history <- sprintf("state_%03d", seq_len(recipe$history_depth))
        immediate <- active_history[[length(active_history) - 1L]]
        reverted <- checkpointUtilityResolveState(context, immediate, cold = TRUE)
        active_history <- active_history[seq_len(length(active_history) - 1L)]
        deep <- active_history[[1L]]
        reverted <- checkpointUtilityResolveState(context, deep, cold = TRUE)
        active_history <- deep
        observe("after_reverts")
        # Reset the active state/lineage to the initial checkpoint, preserving dormant
        # historical records just as WorkflowState$revertToState does.
        current_state <- deep
        current_record <- reverted
        if (!identical(treatment, "resident")) context$last_reference <- context$states[[deep]]
        generation <- generation_offset + recipe$history_depth + 1L
        reapply_name <- "state_002"
        object <- checkpointUtilityProteinState(helpers, recipe, generation,
            variant, reverted$data, run_dir)
        record <- checkpointUtilityStateRecord(object, recipe, reapply_name, deep, generation)
        record <- api("prepareDIAWorkflowHistory")(stats::setNames(list(record), reapply_name))[[1L]]
        context$states[[reapply_name]] <- checkpointUtilityAcceptState(context,
            record, reapply_name, generation)
        active_history <- c(active_history, reapply_name)
        generation <- generation + 1L
        branch_name <- "branch_001"
        object <- checkpointUtilityProteinState(helpers, recipe, generation,
            variant, record$data, run_dir)
        record <- checkpointUtilityStateRecord(object, recipe, branch_name,
            reapply_name, generation)
        record <- api("prepareDIAWorkflowHistory")(stats::setNames(list(record), branch_name))[[1L]]
        context$states[[branch_name]] <- checkpointUtilityAcceptState(context,
            record, branch_name, generation)
        active_history <- c(active_history, branch_name)
        current_state <- branch_name
        current_record <- record
        object <- record <- reverted <- NULL
        observe("after_reapply_before_export")
        export <- checkpointUtilityTimed(metrics, "complete_legacy_export_hydration", \() {
            checkpointUtilityExportEnvelope(context, active_history, current_state)
        })
        # Receipt-only object hashes allocate serialized buffers. Complete the
        # same checks before the existing collection, not after the RDS write.
        state_hashes <- checkpointUtilityTimed(metrics, "complete_legacy_export_state_hashes", \() {
            vapply(export$records, digest::digest, character(1), algo = "sha256")
        })
        current_data_sha256 <- digest::digest(current_record$data, algo = "sha256")
        observe("after_export_state_hashes")
        export_heap <- checkpointUtilityHeapBytes()
        observe("export_after_gc", export_heap)
        checkpointUtilityTimed(metrics, "complete_legacy_export_serialization", \() {
            saveRDS(export$envelope, export_path)
        })
        cycle_metrics[[cycle]] <- list(cycle = cycle, retained_heap_bytes = retained_heap,
            retained_total_heap_bytes = retained_total, export_heap_bytes = export_heap,
            reset_state = deep, current_state = current_state,
            active_history = active_history, states_exported = length(export$records))
        retained_alias_count <- length(aliases)
        # Do not keep an export or obsolete active object alive into the next cycle.
        export <- current_record <- aliases <- NULL
        observe("after_export_release")
    }
    references_path <- NULL
    if (!identical(treatment, "resident")) {
        references_path <- file.path(run_dir, "references.rds")
        saveRDS(context$states, references_path)
        api("closeWorkflowEntityStore")(store)
        store <- NULL
    }
    finished <- proc.time()
    events <- metrics$events
    receipt$status <- "complete"
    receipt$science_ok <- TRUE
    receipt$outputs_ok <- TRUE
    receipt$history_depth <- recipe$history_depth
    receipt$state_names <- names(context$states)
    receipt$active_history <- active_history
    receipt$current_state <- current_state
    receipt$state_hashes <- as.list(state_hashes)
    receipt$current_data_sha256 <- current_data_sha256
    receipt$cycles_completed <- cycles
    receipt$cycle_metrics <- cycle_metrics
    receipt$retained_heap_bytes <- max(vapply(cycle_metrics,
        \(item) item$retained_heap_bytes, numeric(1)))
    receipt$retained_total_heap_bytes <- retained_total
    receipt$pre_history_heap_bytes <- baseline_heap
    receipt$export_heap_bytes <- max(vapply(cycle_metrics,
        \(item) item$export_heap_bytes, numeric(1)))
    receipt$operation_elapsed_seconds <- unname((finished - started)[["elapsed"]])
    receipt$operation_cpu_seconds <- unname(sum((finished - started)[c("user.self", "sys.self")]))
    receipt$max_save_seconds <- metrics$max_save_seconds
    receipt$max_cold_revert_seconds <- metrics$max_cold_revert_seconds
    receipt$payload_read_bytes <- metrics$payload_read_bytes
    receipt$events <- events
    receipt$legacy_export <- list(path = export_path, bytes = file.size(export_path),
        sha256 = digest::digest(file = export_path, algo = "sha256"))
    receipt$references <- if (is.null(references_path)) NULL else list(path = references_path,
        bytes = file.size(references_path), sha256 = digest::digest(file = references_path,
            algo = "sha256"))
    receipt$retained_alias_count <- retained_alias_count
    receipt$heap_accounting <- paste("retained_heap_bytes is the positive incremental live R heap",
        "after full GC; Ncells use 56 bytes and Vcells 8 bytes; it is not RSS or cgroup peak")
    jsonlite::write_json(receipt, receipt_path, auto_unbox = TRUE, pretty = TRUE,
        null = "null", digits = NA)
    complete <- TRUE
    receipt
}

#' Verify one completed utility member in a separate installed-package process.
#' @param library_path Installed candidate library.
#' @param verify_dir Fresh verifier directory.
#' @param member_dir Completed utility member directory.
#' @return Verification receipt.
#' @noRd
verifyCheckpointUtilityMember <- function(library_path, verify_dir, member_dir) {
    library_path <- normalizePath(library_path, mustWork = TRUE)
    verify_dir <- normalizePath(verify_dir, mustWork = TRUE)
    member_dir <- normalizePath(member_dir, mustWork = TRUE)
    .libPaths(unique(c(library_path, .libPaths())))
    namespace <- loadNamespace("MultiScholaR", lib.loc = library_path)
    api <- \(name) get(name, envir = namespace, inherits = FALSE)
    member <- jsonlite::read_json(file.path(member_dir, "member-receipt.json"))
    if (!identical(member$status, "complete") || !isTRUE(member$science_ok) ||
        !isTRUE(member$outputs_ok)) stop("utility member is not complete")
    legacy <- readRDS(file.path(member_dir, "legacy-session.rds"))
    if (!identical(digest::digest(file = file.path(member_dir, "legacy-session.rds"),
        algo = "sha256"), member$legacy_export$sha256)) stop("legacy export fingerprint differs")
    states <- legacy$r6_complete_states
    hashes <- vapply(states, digest::digest, character(1), algo = "sha256")
    if (!identical(as.list(hashes), member$state_hashes) ||
        !identical(unlist(legacy$r6_state_history, use.names = FALSE),
            unlist(member$active_history, use.names = FALSE)) ||
        !identical(legacy$current_state, member$current_state) ||
        !identical(digest::digest(legacy$current_s4_object, algo = "sha256"),
            member$current_data_sha256)) stop("legacy history verification failed")
    for (state in states) methods::validObject(state$data)
    persistent_checked <- 0L
    if (!identical(member$treatment, "resident")) {
        if (!identical(digest::digest(file = file.path(member_dir, "references.rds"),
            algo = "sha256"), member$references$sha256)) stop("checkpoint references fingerprint differs")
        references <- readRDS(file.path(member_dir, "references.rds"))
        store <- api("newWorkflowEntityStore")(member_dir,
            paste(member$recipe, member$variant, member$block, sep = "-"), "DIA")
        on.exit(api("closeWorkflowEntityStore")(store), add = TRUE)
        for (name in names(references)) {
            if (!identical(references[[name]]$payload$codec, member$codec)) {
                stop("checkpoint codec differs from the selected treatment")
            }
            restored <- api("hydrateSavedWorkflowCheckpoint")(store, references[[name]])
            if (!identical(digest::digest(restored, algo = "sha256"),
                member$state_hashes[[name]])) stop("persistent checkpoint differs")
            persistent_checked <- persistent_checked + 1L
        }
    }
    receipt <- list(status = "passed_utility_member", member_dir = member_dir,
        recipe = member$recipe, treatment = member$treatment, variant = member$variant,
        block = member$block, states_checked = length(states),
        persistent_states_checked = persistent_checked,
        legacy_sha256 = digest::digest(file = file.path(member_dir, "legacy-session.rds"),
            algo = "sha256"), promotion_authorized = FALSE)
    jsonlite::write_json(receipt, file.path(verify_dir, "utility-verification.json"),
        auto_unbox = TRUE, pretty = TRUE)
    receipt
}

#' Measure one standalone consumer operation with CPU and elapsed time.
#' @param label Operation label.
#' @param operation Zero-argument operation.
#' @return Value and descriptive cost; never a performance-gate estimate.
#' @noRd
checkpointConsumerOperation <- function(label, operation) {
    start <- proc.time()
    value <- operation()
    finish <- proc.time()
    list(label = label, value = value,
        elapsed_seconds = unname((finish - start)[["elapsed"]]),
        cpu_seconds = unname(sum((finish - start)[c("user.self", "sys.self")])))
}

#' Load exactly one installed package for the recorded-state pilot.
#' @param library_path Explicit installed library.
#' @return Internal function lookup closure.
#' @noRd
recordedCheckpointApi <- function(library_path) {
    library_path <- normalizePath(library_path, mustWork = TRUE)
    .libPaths(unique(c(library_path, .libPaths())))
    namespace <- loadNamespace("MultiScholaR", lib.loc = library_path)
    if (!identical(normalizePath(getNamespaceInfo(namespace, "path")),
        normalizePath(file.path(library_path, "MultiScholaR")))) {
        stop("recorded pilot loaded a different installed package")
    }
    \(name) get(name, envir = namespace, inherits = FALSE)
}

#' Preserve each recorded section's worker settings without merging runtime handles.
#' @param record Complete legacy checkpoint.
#' @param api Installed preparation helpers.
#' @return Portable record; no worker is created or stopped.
#' @noRd
prepareRecordedCheckpointRecord <- function(record, api) {
    record <- api("mapDIAStateConfig")(record, function(config) {
        for (section in names(config)) {
            value <- config[[section]]
            if (!is.list(value) || !inherits(value$core_utilisation, "multidplyr_cluster")) next
            partial <- stats::setNames(list(value), section)
            # A historical handle records that section's actual worker count.
            # Preserve the separate global setting, even when it differs; this
            # storage replay does not choose or launch an execution backend.
            config[[section]] <- api("splitDIAWorkerConfig")(partial)$config[[section]]
        }
        config
    })
    api("prepareDIAWorkflowHistory")(list(recorded = record))[[1L]]
}

#' Extract complete portable records without changing the source session.
#' @param session Decoded verified legacy session, owned by preparation only.
#' @param run_dir Fresh preparation directory.
#' @param api Installed internal lookup closure.
#' @param source_hash Verified original session SHA-256.
#' @return Small manifest; no checkpoint values are returned.
#' @noRd
splitRecordedCheckpointSession <- function(session, run_dir, api, source_hash) {
    states <- session$r6_complete_states
    names_expected <- unlist(session$r6_state_history, use.names = FALSE)
    if (!is.list(states) || anyDuplicated(names(states)) ||
        !all(c("initial", "raw_data_s4", "qvalue_filtered", "correlation_filtered") %in%
            names(states)) || !identical(names(states), names_expected)) {
        stop("recorded session history differs from the declared replay")
    }
    source_bytes <- as.numeric(lobstr::obj_size(states))
    records <- vector("list", length(states))
    names(records) <- names(states)
    for (index in seq_along(states)) {
        name <- names(states)[[index]]
        record <- prepareRecordedCheckpointRecord(states[[index]], api)
        api("checkSavedWorkflowPayload")(record)
        file <- sprintf("state-%02d.rds", index)
        path <- file.path(run_dir, file)
        if (file.exists(path)) stop("refusing to overwrite a prepared checkpoint")
        saveRDS(record, path, compress = FALSE, version = 3)
        shape <- if (methods::is(record$data, "PeptideQuantitativeData"))
            dim(methods::slot(record$data, "peptide_data")) else
            if (methods::is(record$data, "ProteinQuantitativeData"))
                dim(methods::slot(record$data, "protein_quant_table")) else integer()
        records[[index]] <- list(file = file, bytes = file.size(path),
            file_sha256 = digest::digest(file = path, algo = "sha256"),
            record_sha256 = digest::digest(record, algo = "sha256"),
            object_bytes = as.numeric(lobstr::obj_size(record)),
            class = class(record$data), shape = as.list(shape))
    }
    session$r6_complete_states <- NULL
    session$current_s4_object <- NULL
    session$correlation_filtered_s4 <- NULL
    template_path <- file.path(run_dir, "session-template.rds")
    saveRDS(session, template_path, compress = FALSE, version = 3)
    qc_config_path <- file.path(run_dir, "qvalue-config.rds")
    saveRDS(list(config = states$qvalue_filtered$config,
        timestamp = states$qvalue_filtered$timestamp), qc_config_path)
    manifest <- list(kind = "recorded_dia_history_pilot", status = "prepared",
        source_sha256 = source_hash, state_names = names_expected, records = records,
        current_state = session$r6_current_state_name, source_history_bytes = source_bytes,
        template_sha256 = digest::digest(file = template_path, algo = "sha256"),
        qvalue_config_sha256 = digest::digest(file = qc_config_path, algo = "sha256"),
        preparation_note = "Per-section actual worker counts and separate global CPU settings preserved independently; existing parser portability applied identically to both arms; no live worker rebind",
        scope = "separately restored recorded states; no forward-science or GUI performance claim",
        promotion_authorized = FALSE)
    jsonlite::write_json(manifest, file.path(run_dir, "recorded-inputs.json"),
        auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
    manifest
}

#' Prepare the exact approved 16-state DIA session in its own process.
#' @param library_path,run_dir Installed library and fresh preparation directory.
#' @param session_path Approved source session.
#' @return Prepared manifest.
#' @noRd
prepareRecordedCheckpointSession <- function(library_path, run_dir, session_path) {
    expected <- "4671bdfab031cd880ab2bd2b414c4625cc30902a247e7d60f1c8dc11411a50c3"
    if (!identical(digest::digest(file = session_path, algo = "sha256"), expected)) {
        stop("source session differs from approved rzEEWj session")
    }
    api <- recordedCheckpointApi(library_path)
    session <- readRDS(session_path)
    if (length(session$r6_complete_states) != 16L) stop("expected sixteen recorded states")
    splitRecordedCheckpointSession(session, run_dir, api, expected)
}

#' Read the small prepared-input manifest without hydrating a checkpoint.
#' @param prepared_dir Verified preparation directory.
#' @return Manifest with validated file names.
#' @noRd
readRecordedCheckpointManifest <- function(prepared_dir) {
    manifest <- jsonlite::read_json(file.path(prepared_dir, "recorded-inputs.json"))
    names_expected <- unlist(manifest$state_names, use.names = FALSE)
    if (!identical(manifest$kind, "recorded_dia_history_pilot") ||
        !identical(manifest$status, "prepared") || !length(names_expected) ||
        anyDuplicated(names_expected) || !identical(names(manifest$records), names_expected)) {
        stop("invalid prepared history manifest")
    }
    for (record in manifest$records) {
        if (!is.character(record$file) || length(record$file) != 1L ||
            !grepl("^state-[0-9]{2}\\.rds$", record$file)) stop("invalid prepared record path")
    }
    manifest
}

#' Load exactly one checksum-bound input record.
#' @param prepared_dir,manifest Prepared input directory and manifest.
#' @param name Exact checkpoint name.
#' @return Complete record, owned only by the caller.
#' @noRd
readRecordedCheckpoint <- function(prepared_dir, manifest, name) {
    spec <- manifest$records[[name]]
    if (is.null(spec)) stop("unknown recorded checkpoint")
    path <- file.path(prepared_dir, spec$file)
    if (!identical(as.numeric(file.size(path)), as.numeric(spec$bytes)) ||
        !identical(digest::digest(file = path, algo = "sha256"), spec$file_sha256)) {
        stop("prepared checkpoint integrity differs")
    }
    readRDS(path)
}

#' Observe retained live R allocations after obsolete action owners were released.
#' @param context Replay owner environment.
#' @param label Completed action name.
#' @param memory_fn Existing service snapshot function.
#' @param inspect_graph Traverse the complete object graph only in ownership diagnostics.
#' @return Small scalar-only observation.
#' @noRd
sampleRecordedCheckpointOwners <- function(context, label, memory_fn, inspect_graph = TRUE) {
    started <- proc.time()
    heap <- checkpointUtilityHeapBytes()
    owned <- if (inspect_graph)
        as.numeric(lobstr::obj_size(context$states, context$current)) else NULL
    elapsed <- proc.time() - started
    list(label = label, heap_bytes = heap, incremental_heap_bytes = heap - context$baseline_heap,
        owned_state_bytes = owned, history_entries = length(context$states),
        observation_elapsed_seconds = unname(elapsed[["elapsed"]]),
        observation_cpu_seconds = unname(sum(elapsed[c("user.self", "sys.self")])),
        service = memory_fn())
}

#' Replay the same recorded states under resident or persisted ownership.
#' @param library_path,run_dir Installed library and fresh member directory.
#' @param prepared_dir Accepted individual checkpoint fixtures.
#' @param treatment resident, saved (native RDS), or qs2 (explicit checkpoint experiment).
#' @param memory_fn Existing guarded-service observation function.
#' @param profile_owners Collect per-state GC/unique-owner observations only in diagnostics.
#' @return Diagnostic receipt; formal G2 is never inferred from this pair.
#' @noRd
executeRecordedCheckpointReplay <- function(library_path, run_dir, prepared_dir,
    treatment, memory_fn, profile_owners = TRUE) {
    if (!treatment %in% c("resident", "saved", "qs2")) stop("invalid recorded replay treatment")
    if (!is.logical(profile_owners) || length(profile_owners) != 1L || is.na(profile_owners)) {
        stop("invalid ownership-observation policy")
    }
    api <- recordedCheckpointApi(library_path)
    manifest <- readRecordedCheckpointManifest(prepared_dir)
    state_names <- unlist(manifest$state_names, use.names = FALSE)
    metrics <- new.env(parent = emptyenv())
    metrics$events <- list()
    context <- new.env(parent = emptyenv())
    context$baseline_heap <- checkpointUtilityHeapBytes()
    context$current <- NULL
    context$states <- list()
    store <- if (treatment != "resident") api("newWorkflowEntityStore")(
        run_dir, "recorded-replay", "DIA", max_cache_bytes = 0) else NULL
    if (!is.null(store)) on.exit(api("closeWorkflowEntityStore")(store), add = TRUE)
    accept <- function(name, record, index) {
        if (treatment == "resident") return(record)
        identity <- list(input = manifest$source_sha256, design = "recorded-design",
            parameters = paste0("recorded-state-", index), selectors = name,
            software = as.character(utils::packageVersion("MultiScholaR")))
        checkpointUtilityTimed(metrics, paste0("save:", name), \() {
            arguments <- list(store, record, sprintf("record-%02d", index),
                identity, list(class = "list", shape = integer(), factors = list()))
            if (treatment == "qs2") arguments$codec <- "qs2"
            do.call(api("saveWorkflowCheckpoint"), arguments)
        })
    }
    resolve <- function(name) {
        if (treatment == "resident") return(context$states[[name]])
        checkpointUtilityTimed(metrics, paste0("hydrate:", name), \() {
            api("hydrateSavedWorkflowCheckpoint")(store, context$states[[name]])
        })
    }
    observations <- list()
    for (index in seq_along(state_names)) {
        name <- state_names[[index]]
        record <- checkpointUtilityTimed(metrics, paste0("input:", name), \() {
            readRecordedCheckpoint(prepared_dir, manifest, name)
        })
        context$states[[name]] <- accept(name, record, index)
        context$current <- record
        record <- NULL
        observations[[length(observations) + 1L]] <- if (profile_owners)
            sampleRecordedCheckpointOwners(context, name, memory_fn) else
            list(label = name, history_entries = length(context$states), service = memory_fn())
    }
    retained <- if (profile_owners) observations[[length(observations)]] else
        sampleRecordedCheckpointOwners(context, "forward_complete", memory_fn, inspect_graph = FALSE)
    for (name in c(state_names[[length(state_names) - 1L]], "raw_data_s4",
        "correlation_filtered")) {
        context$current <- checkpointUtilityTimed(metrics, paste0("revert:", name),
            \() resolve(name))
    }
    record <- context$current
    record$description <- "Recorded-state persistence metadata branch"
    record$previous_state <- "correlation_filtered"
    record$config$recorded_history_branch <- TRUE
    branch <- "replay_metadata_branch"
    context$states[[branch]] <- accept(branch, record, length(state_names) + 1L)
    context$current <- record
    record <- NULL
    template_path <- file.path(prepared_dir, "session-template.rds")
    if (!identical(digest::digest(file = template_path, algo = "sha256"),
        manifest$template_sha256)) stop("session template integrity differs")
    session <- readRDS(template_path)
    session$r6_complete_states <- checkpointUtilityTimed(metrics, "export_hydration", \() {
        stats::setNames(lapply(names(context$states), resolve), names(context$states))
    })
    session$r6_current_state_name <- branch
    session$r6_state_history <- as.list(names(context$states))
    session$current_s4_object <- context$current$data
    session$correlation_filtered_s4 <- session$r6_complete_states$correlation_filtered$data
    export_heap <- checkpointUtilityHeapBytes()
    export_memory <- memory_fn()
    export_dir <- file.path(run_dir, "export")
    if (!dir.create(export_dir)) stop("cannot create fresh recorded export directory")
    export_time <- session$export_timestamp
    files <- checkpointUtilityTimed(metrics, "legacy_export", \() {
        api("saveProtNormExportArtifacts")(session, export_dir,
            timeFn = \() export_time, messageFn = \(...) invisible(NULL))
    })
    if (treatment != "resident") saveRDS(context$states, file.path(run_dir, "references.rds"))
    session <- NULL
    context$current <- NULL
    context$states <- list()
    if (!is.null(store)) api("closeWorkflowEntityStore")(store)
    released_heap <- checkpointUtilityHeapBytes()
    receipt <- list(kind = "recorded_dia_history_pilot", status = "complete",
        treatment = treatment, prepared_dir = normalizePath(prepared_dir),
        checkpoint_codec = if (treatment == "qs2") "qs2" else "rds",
        qs2_version = if (treatment == "qs2") as.character(utils::packageVersion("qs2")) else NULL,
        input_manifest_sha256 = digest::digest(file = file.path(prepared_dir,
            "recorded-inputs.json"), algo = "sha256"), source_sha256 = manifest$source_sha256,
        observations = observations, retained_heap_bytes = retained$incremental_heap_bytes,
        retained_observation = retained, profile_owners = profile_owners,
        retained_owned_bytes = retained$owned_state_bytes, export_heap_bytes = export_heap,
        export_memory = export_memory, released_heap_bytes = released_heap,
        events = metrics$events, files = files, final_state = branch,
        source_states = state_names, history_entries = length(state_names) + 1L,
        optional_cache_bytes = 0, open_resources = if (is.null(store)) 0L else store$open_resources,
        scientific_execution = FALSE, independently_verified = FALSE,
        promotion_authorized = FALSE)
    jsonlite::write_json(receipt, file.path(run_dir, "recorded-member.json"),
        auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
    receipt
}

#' Independently verify recorded history, native exports and durable checkpoints.
#' @param library_path,run_dir Installed library and fresh verifier directory.
#' @param member_dir Completed replay directory.
#' @return Exactness receipt without a performance claim.
#' @noRd
verifyRecordedCheckpointReplay <- function(library_path, run_dir, member_dir) {
    api <- recordedCheckpointApi(library_path)
    member <- jsonlite::read_json(file.path(member_dir, "recorded-member.json"))
    if (!identical(member$status, "complete")) stop("recorded member is incomplete")
    manifest <- readRecordedCheckpointManifest(member$prepared_dir)
    if (!identical(digest::digest(file = file.path(member$prepared_dir,
        "recorded-inputs.json"), algo = "sha256"), member$input_manifest_sha256)) {
        stop("recorded member input binding differs")
    }
    latest <- member$files$latestFilepath
    timestamped <- member$files$sessionFilepath
    if (!identical(digest::digest(file = latest, algo = "sha256"),
        digest::digest(file = timestamped, algo = "sha256"))) stop("legacy exports differ")
    session <- readRDS(latest)
    template_path <- file.path(member$prepared_dir, "session-template.rds")
    if (!identical(digest::digest(file = template_path, algo = "sha256"),
        manifest$template_sha256)) stop("session template integrity differs")
    template <- readRDS(template_path)
    metadata_fields <- setdiff(names(template), c("r6_current_state_name", "r6_state_history"))
    if (!identical(serialize(session[metadata_fields], NULL),
        serialize(template[metadata_fields], NULL))) stop("legacy session metadata differs")
    metadata_files <- c(accession_cleanup_results = "accession_cleanup_results.RDS",
        qc_params = "qc_params.RDS", fasta_metadata = "fasta_metadata.RDS",
        ruv_optimization_result = "ruv_optimization_results.RDS",
        mixed_species_analysis = "mixed_species_analysis.RDS",
        limpa_parameters = "limpa_parameters.RDS", limpa_dpc_quant_results = "limpa_dpc_quant_results.RDS")
    for (field in names(metadata_files)) if (!is.null(template[[field]])) {
        path <- file.path(dirname(latest), metadata_files[[field]])
        if (!file.exists(path) || !identical(serialize(readRDS(path), NULL),
            serialize(template[[field]], NULL))) stop(paste("export metadata differs:", field))
    }
    expected_summary <- api("buildProtNormExportSummaryContent")(session,
        basename(timestamped), timeFn = \() template$export_timestamp)
    if (!identical(paste(readLines(member$files$summaryFilepath), collapse = "\n"),
        expected_summary)) stop("legacy summary content differs")
    expected_names <- c(names(manifest$records), "replay_metadata_branch")
    if (!identical(names(session$r6_complete_states), expected_names) ||
        !identical(unlist(session$r6_state_history, use.names = FALSE), expected_names) ||
        !identical(session$r6_current_state_name, "replay_metadata_branch")) {
        stop("exported recorded lineage differs")
    }
    references <- if (member$treatment != "resident")
        readRDS(file.path(member_dir, "references.rds")) else NULL
    store <- if (is.null(references)) NULL else api("newWorkflowEntityStore")(
        member_dir, "recorded-replay", "DIA", max_cache_bytes = 0)
    if (!is.null(store)) on.exit(api("closeWorkflowEntityStore")(store), add = TRUE)
    for (name in expected_names) {
        expected <- readRecordedCheckpoint(member$prepared_dir, manifest,
            if (name == "replay_metadata_branch") "correlation_filtered" else name)
        if (name == "replay_metadata_branch") {
            expected$description <- "Recorded-state persistence metadata branch"
            expected$previous_state <- "correlation_filtered"
            expected$config$recorded_history_branch <- TRUE
        }
        actual <- session$r6_complete_states[[name]]
        if (!identical(actual, expected)) stop(paste("legacy record differs:", name))
        if (isS4(actual$data)) methods::validObject(actual$data)
        if (!is.null(store)) {
            restored <- api("hydrateSavedWorkflowCheckpoint")(store, references[[name]])
            if (!identical(restored, expected)) stop(paste("durable record differs:", name))
            restored <- NULL
        }
        actual <- expected <- NULL
        gc()
    }
    if (!identical(session$current_s4_object,
        session$r6_complete_states$replay_metadata_branch$data) ||
        !identical(session$correlation_filtered_s4,
            session$r6_complete_states$correlation_filtered$data)) stop("export aliases differ")
    manager <- api("WorkflowState")$new()
    manager$states <- session$r6_complete_states
    manager$state_history <- session$r6_state_history
    manager$current_state <- session$r6_current_state_name
    if (!identical(manager$revertToState("raw_data_s4"),
        session$r6_complete_states$raw_data_s4$data)) stop("legacy revert differs")
    result <- list(status = "passed_recorded_history", member_dir = normalizePath(member_dir),
        treatment = member$treatment, checked_states = length(expected_names),
        full_metadata_verified = TRUE,
        legacy_sha256 = digest::digest(file = latest, algo = "sha256"),
        promotion_authorized = FALSE)
    jsonlite::write_json(result, file.path(run_dir, "recorded-verification.json"),
        auto_unbox = TRUE, pretty = TRUE)
    result
}

#' Extract only settings for the approved forward QC sequence outside timing.
#' @param library_path,run_dir,prepared_dir Installed package, output and recorded inputs.
#' @return Small immutable settings manifest; downstream objects are not retained.
#' @noRd
prepareRecordedQcSequence <- function(library_path, run_dir, prepared_dir) {
    api <- recordedCheckpointApi(library_path)
    manifest <- readRecordedCheckpointManifest(prepared_dir)
    states <- c("qvalue_filtered", "protein_peptide_filtered", "precursor_rollup")
    settings <- stats::setNames(lapply(states, function(name) {
        record <- readRecordedCheckpoint(prepared_dir, manifest, name)
        list(config = record$config, timestamp = record$timestamp)
    }), states)
    file <- file.path(run_dir, "qc-sequence-settings.rds")
    saveRDS(settings, file, version = 3)
    if (file.size(file) > 1048576) stop("sequence settings must remain bounded")
    receipt <- list(status = "prepared_qc_sequence", prepared_dir = normalizePath(prepared_dir),
        source_sha256 = manifest$source_sha256,
        recorded_manifest_sha256 = digest::digest(file = file.path(prepared_dir,
            "recorded-inputs.json"), algo = "sha256"),
        settings_sha256 = digest::digest(file = file, algo = "sha256"),
        states = states, raw_shape = manifest$records$raw_data_s4$shape,
        promotion_authorized = FALSE)
    jsonlite::write_json(receipt, file.path(run_dir, "qc-sequence-inputs.json"),
        auto_unbox = TRUE, pretty = TRUE, digits = NA)
    receipt
}

#' Create the proof-only state facade; public WorkflowState remains unchanged.
#' @param context Member-owned environment with installed API, store and settings.
#' @return Minimal interface consumed by the three production Apply helpers.
#' @noRd
recordedQcSequenceManager <- function(context) {
    force(context)
    facade <- new.env(parent = emptyenv())
    facade$current_state <- "raw_data_s4"
    facade$audit_enabled <- TRUE
    facade$getState <- function(state_name = NULL) {
        name <- if (is.null(state_name)) facade$current_state else state_name
        record <- recordedQcSequenceResolve(context, name)
        context$current <- record
        record$data
    }
    facade$saveState <- function(state_name, s4_data_object, config_object, description,
        audit_metadata = NULL) {
        native <- context$native
        prior <- list(current = native$current_state, history = native$state_history,
            states = native$states, audit = native$audit_records)
        committed <- FALSE
        on.exit(if (!committed) {
            native$current_state <- prior$current
            native$state_history <- prior$history
            native$states <- prior$states
            native$audit_records <- prior$audit
        }, add = TRUE)
        native$current_state <- facade$current_state
        native$saveState(state_name, s4_data_object, config_object, description, audit_metadata)
        native$states[[state_name]]$timestamp <- context$settings[[state_name]]$timestamp
        record <- prepareRecordedCheckpointRecord(native$states[[state_name]], context$api)
        native$states[[state_name]] <- record
        recordedQcSequenceAccept(context, state_name, record)
        facade$current_state <- state_name
        committed <- TRUE
        invisible(NULL)
    }
    facade$revertToState <- function(state_name) {
        if (!state_name %in% context$names) stop("sequence state not found")
        record <- recordedQcSequenceResolve(context, state_name)
        facade$current_state <- state_name
        context$native$current_state <- state_name
        index <- match(state_name, unlist(context$native$state_history))
        if (!is.na(index)) context$native$state_history <- context$native$state_history[seq_len(index)]
        context$current <- record
        record$data
    }
    facade
}

#' Persist one complete checkpoint or retain its genuine resident allocation.
#' @param context Member owner.
#' @param name,record Logical state and complete portable legacy record.
#' @return NULL invisibly; publish references only after durable acceptance.
#' @noRd
recordedQcSequenceAccept <- function(context, name, record) {
    if (identical(context$treatment, "qs2")) {
        reference <- checkpointUtilityTimed(context$metrics, paste0("save:", name), function() {
            context$api("saveWorkflowCheckpoint")(context$store, record,
                paste0("sequence-", name),
                list(input = context$source_sha256, design = "recorded-design",
                    parameters = name, selectors = "qc-three-step",
                    software = as.character(utils::packageVersion("MultiScholaR"))),
                list(class = "list", shape = integer(), factors = list()), codec = "qs2")
        })
        context$references[[name]] <- reference
        saveRDS(context$references, file.path(context$run_dir, "qc-sequence-references.rds"))
    } else context$native$states[[name]] <- record
    context$names <- unique(c(context$names, name))
    invisible(NULL)
}

#' Resolve one existing state without ever loading the recorded downstream oracle.
#' @param context,name Member owner and state name.
#' @return Complete current record.
#' @noRd
recordedQcSequenceResolve <- function(context, name) {
    if (!name %in% context$names) stop("sequence state not found")
    if (identical(context$treatment, "resident")) return(context$native$states[[name]])
    checkpointUtilityTimed(context$metrics, paste0("hydrate:", name), function() {
        context$api("hydrateSavedWorkflowCheckpoint")(context$store, context$references[[name]])
    })
}

#' Release every action-scoped payload, including configuration and audit aliases.
#' @param context Member owner.
#' @return NULL invisibly; resident histories are deliberately untouched.
#' @noRd
releaseRecordedQcSequence <- function(context) {
    context$current <- NULL
    context$rollup_result <- NULL
    context$config_env$config_list <- NULL
    if (identical(context$treatment, "qs2")) {
        context$native$states <- list()
        context$native$audit_records <- list()
    }
    invisible(NULL)
}

#' Bind only the executing section's unchanged worker specification, owning startup.
#' @param context,object Member resources and current S4 object.
#' @param section Scientific function being executed.
#' @return Executable S4 with the original worker count, or unchanged serial input.
#' @noRd
bindRecordedQcSequence <- function(context, object, section) {
    args <- methods::slot(object, "args")
    spec <- args[[section]]$core_utilisation
    if (inherits(spec, "DIAWorkerSpec")) {
        key <- as.character(spec$workers)
        if (is.null(context$clusters[[key]])) {
            context$clusters[[key]] <- context$api("startDIAWorkerCluster")(spec)
        }
        args[[section]] <- context$api("bindDIAWorkerConfig")(
            stats::setNames(list(args[[section]]), section), context$clusters[[key]])[[section]]
        methods::slot(object, "args") <- args
    }
    object
}

#' Run one real forward scientific action, audit, state commit and existing outputs.
#' @param context,workflow Member resources and facade-backed workflow.
#' @param name Exact stage name.
#' @return NULL invisibly; only display handles survive the call.
#' @noRd
runRecordedQcSequenceStep <- function(context, workflow, name) {
    api <- context$api
    setting <- context$settings[[name]]
    config <- setting$config
    update_config <- function(...) api("updateConfigParameter")(..., env = context$config_env)
    now <- function() setting$timestamp
    section <- switch(name, qvalue_filtered = "srlQvalueProteotypicPeptideClean",
        protein_peptide_filtered = "filterMinNumPeptidesPerProtein",
        precursor_rollup = "rollUpPrecursorToPeptide")
    if (isTRUE(context$sql_rollup) && name == "precursor_rollup") {
        context$rollup_result <- checkpointUtilityTimed(context$metrics,
            "query:precursor_rollup", function() {
                api("computeSavedWorkflowPrecursorRollup")(context$store, context$rollup_input)
            })
        if (context$store$open_resources != 0L) stop("rollup query retained a database lease")
    }
    # Getter owns one hydration for this action; the Apply helper consumes that object.
    current <- workflow$state_manager$getState()
    current <- bindRecordedQcSequence(context, current, section)
    context$config_env$config_list <- methods::slot(current, "args")
    original_get <- workflow$state_manager$getState
    workflow$state_manager$getState <- function(state_name = NULL) {
        if (is.null(state_name)) current else original_get(state_name)
    }
    on.exit(workflow$state_manager$getState <- original_get, add = TRUE)
    qvalueFilterFn <- api("srlQvalueProteotypicPeptideClean")
    filterMinNumPeptidesPerProteinFn <- api("filterMinNumPeptidesPerProtein")
    rollupFn <- api(if (isTRUE(context$sql_rollup)) "applySavedWorkflowPrecursorRollup" else
        "rollUpPrecursorToPeptide")
    scientific <- function(theObject, ...) {
        checkpointUtilityTimed(context$metrics, paste0("method:", name), function() {
            # Existing parameter provenance reads the call symbol, not the generic.
            # Preserve each production Apply helper's exact invocation name.
            switch(name,
                qvalue_filtered = qvalueFilterFn(theObject, ...),
                protein_peptide_filtered = filterMinNumPeptidesPerProteinFn(theObject, ...),
                precursor_rollup = if (isTRUE(context$sql_rollup))
                    rollupFn(theObject, rollup_result = context$rollup_result, ...) else
                    rollupFn(theObject, ...))
        })
    }
    result <- checkpointUtilityTimed(context$metrics, paste0("apply:", name), function() {
        switch(name,
            qvalue_filtered = api("runPeptideQvalueApplyStep")(workflow,
                config$qvalue_threshold, config$global_qvalue_threshold, config$proteotypic_only,
                updateConfigParameterFn = update_config, qvalueFilterFn = scientific,
                globalPGQvalueThreshold = config$global_pg_qvalue_threshold,
                confidenceProvenanceMode = config$confidence_provenance_mode, nowFn = now),
            protein_peptide_filtered = api("runProteinPeptideApplyStep")(workflow,
                config$min_peptides_per_protein, config$min_peptidoforms_per_protein,
                updateConfigParameterFn = update_config,
                filterMinNumPeptidesPerProteinFn = scientific, nowFn = now),
            precursor_rollup = api("runPeptideRollupApplyStep")(workflow,
                rollupFn = scientific, nowFn = now))
    })
    if (isTRUE(context$sql_rollup) && name == "protein_peptide_filtered") {
        context$rollup_input <- checkpointUtilityTimed(context$metrics,
            "persist:rollup_input", function() {
                api("saveWorkflowPrecursorRollupInput")(context$store, result$filteredS4,
                    "sequence-rollup-input", list(
                        input = context$references[[name]]$payload$sha256,
                        design = "recorded-design", parameters = "native-rollup-contract",
                        selectors = "complete-population-rollup-columns",
                        software = as.character(utils::packageVersion("MultiScholaR"))))
            })
    }
    checkpointUtilityTimed(context$metrics, paste0("outputs:", name), function() {
        writer <- switch(name, qvalue_filtered = "updatePeptideQvalueOutputs",
            protein_peptide_filtered = "updateProteinPeptideOutputs",
            precursor_rollup = "updatePeptideRollupOutputs")
        api(writer)(context$output, function(value) context$panels[[name]] <- value,
            result, "DIA", "sequence", projectDirs = list(DIA_sequence = context$paths))
    })
    if (is.null(context$panels[[name]])) stop("required QC plot construction failed")
    invisible(NULL)
}

#' Detect assay aliases in retained display handles without scanning namespace roots.
#' @param value Display owners or an intentional alias negative control.
#' @return TRUE when a quantitative S4 or abundance table remains reachable.
#' @noRd
recordedQcSequenceHasPayload <- function(value) {
    seen <- character()
    root <- environment(recordedQcSequenceHasPayload)
    visit <- function(x) {
        if (methods::is(x, "PeptideQuantitativeData") || methods::is(x, "ProteinQuantitativeData"))
            return(TRUE)
        if (is.data.frame(x) && any(c("Precursor.Quantity", "Peptide.RawQuantity", "Q.Value") %in% names(x)))
            return(TRUE)
        if (is.function(x)) return(!is.primitive(x) && visit(environment(x)))
        if (is.environment(x)) {
            if (isNamespace(x) || isBaseNamespace(x) || identical(x, globalenv()) ||
                identical(x, baseenv()) || identical(x, emptyenv()) || identical(x, root)) return(FALSE)
            address <- lobstr::obj_addr(x)
            if (address %in% seen) return(FALSE)
            seen <<- c(seen, address)
            if (any(vapply(as.list(x, all.names = TRUE), visit, logical(1)))) return(TRUE)
            return(visit(parent.env(x)))
        }
        if (isS4(x)) return(any(vapply(methods::slotNames(x),
            function(name) visit(methods::slot(x, name)), logical(1))))
        # Classed grid units subset to another unit, including at length one.
        # Inspect their storage instead of recursively dispatching [[.unit.
        if (is.list(x) || is.pairlist(x)) return(any(vapply(unclass(x), visit, logical(1))))
        FALSE
    }
    visit(value)
}

#' Observe owner graphs only in the separately declared diagnostic pair.
#' @param context,label Member owner and completed boundary.
#' @return Small scalar observation.
#' @noRd
observeRecordedQcSequence <- function(context, label) {
    display_payload <- if (context$diagnostic) recordedQcSequenceHasPayload(
        list(context$output, context$panels)) else NULL
    if (isTRUE(display_payload)) stop("sequential QC display retains an assay payload")
    heap <- if (context$diagnostic) checkpointUtilityHeapBytes() else NULL
    owners <- if (context$diagnostic) as.numeric(lobstr::obj_size(
        context$native, context$current, context$config_env, context$output,
        context$panels, context$references, context$store)) else NULL
    list(label = label, timestamp = as.numeric(Sys.time()), heap_bytes = heap,
        unique_owner_bytes = owners, service = context$memory_fn(),
        display_payload_detected = display_payload,
        resident_history_entries = length(context$native$states),
        resident_audit_entries = length(context$native$audit_records),
        current_released = is.null(context$current),
        config_released = is.null(context$config_env$config_list))
}

#' Execute the approved three-stage history pilot, not a full GUI qualification.
#' @param library_path,run_dir,input_dir Installed library, output and settings preparation.
#' @param treatment resident or qs2; diagnostic enables ownership/profiling only.
#' @param memory_fn Existing guarded-service accounting callback.
#' @return Durable receipt with references and separately verified outputs.
#' @noRd
executeRecordedQcSequence <- function(library_path, run_dir, input_dir, treatment,
    diagnostic, memory_fn) {
    stopifnot(treatment %in% c("resident", "qs2", "qs2_duckdb"), is.logical(diagnostic), length(diagnostic) == 1L)
    api <- recordedCheckpointApi(library_path)
    if (treatment != "resident" && as.character(utils::packageVersion("qs2")) != "0.2.2")
        stop("sequence requires frozen qs2 0.2.2")
    input <- jsonlite::read_json(file.path(input_dir, "qc-sequence-inputs.json"))
    prepared <- input$prepared_dir
    if (!identical(input$status, "prepared_qc_sequence") ||
        !identical(digest::digest(file = file.path(prepared, "recorded-inputs.json"),
            algo = "sha256"), input$recorded_manifest_sha256) ||
        !identical(digest::digest(file = file.path(input_dir, "qc-sequence-settings.rds"),
            algo = "sha256"), input$settings_sha256)) stop("sequence input integrity differs")
    manifest <- readRecordedCheckpointManifest(prepared)
    context <- new.env(parent = emptyenv())
    context$api <- api
    context$run_dir <- run_dir
    context$source_sha256 <- manifest$source_sha256
    context$treatment <- if (treatment == "qs2_duckdb") "qs2" else treatment
    context$sql_rollup <- treatment == "qs2_duckdb"
    context$rollup_input <- context$rollup_result <- NULL
    context$diagnostic <- diagnostic
    context$memory_fn <- memory_fn
    context$settings <- readRDS(file.path(input_dir, "qc-sequence-settings.rds"))
    context$native <- api("WorkflowState")$new()
    context$references <- context$panels <- context$clusters <- list()
    context$names <- character()
    context$config_env <- new.env(parent = emptyenv())
    context$output <- new.env(parent = emptyenv())
    context$store <- if (context$treatment == "qs2") api("newWorkflowEntityStore")(
        run_dir, "qc-sequence", "DIA", max_cache_bytes = 0) else NULL
    context$metrics <- new.env(parent = emptyenv())
    context$metrics$events <- list()
    context$metrics$observe <- function(label) observeCheckpointUtilityMemory(
        file.path(run_dir, "qc-sequence-phases.jsonl"), label, memory_fn)
    receipt <- list(kind = "real_qc_three_step", status = "running", treatment = treatment,
        diagnostic = diagnostic, input_dir = normalizePath(input_dir),
        input_sha256 = digest::digest(file = file.path(input_dir, "qc-sequence-inputs.json"),
            algo = "sha256"), prepared_dir = prepared, completed_steps = character(),
        observations = list(), promotion_authorized = FALSE, independently_verified = FALSE)
    receipt_path <- file.path(run_dir, "qc-sequence-member.json")
    flush_receipt <- function() jsonlite::write_json(receipt, receipt_path,
        auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
    # Private process: preserve any pre-existing progress owner used by small tests.
    had_progress <- exists("filtering_progress", envir = .GlobalEnv, inherits = FALSE)
    prior_progress <- if (had_progress) get("filtering_progress", envir = .GlobalEnv) else NULL
    assign("filtering_progress", methods::new("FilteringProgress"), envir = .GlobalEnv)
    grDevices::pdf(NULL)
    device <- grDevices::dev.cur()
    on.exit({
        if (diagnostic) utils::Rprof(NULL)
        if (device %in% grDevices::dev.list()) grDevices::dev.off(device)
        for (cluster in context$clusters) api("stopDIAWorkerCluster")(cluster)
        if (!is.null(context$store)) api("closeWorkflowEntityStore")(context$store)
        releaseRecordedQcSequence(context)
        if (had_progress) assign("filtering_progress", prior_progress, envir = .GlobalEnv) else
            rm("filtering_progress", envir = .GlobalEnv)
        receipt$events <- context$metrics$events
        receipt$open_resources <- if (is.null(context$store)) 0L else context$store$open_resources
        flush_receipt()
    }, add = TRUE)
    context$paths <- list(time_dir = file.path(run_dir, "plots"),
        publication_graphs_dir = file.path(run_dir, "plots"))
    dir.create(context$paths$time_dir)
    flush_receipt()
    tryCatch({
        receipt$observations$baseline <- observeRecordedQcSequence(context, "baseline")
        for (name in c("initial", "raw_data_s4")) {
            record <- readRecordedCheckpoint(prepared, manifest, name)
            context$native$states[[name]] <- record
            recordedQcSequenceAccept(context, name, record)
        }
        record <- NULL
        context$native$current_state <- "raw_data_s4"
        context$native$state_history <- as.list(context$names)
        releaseRecordedQcSequence(context)
        workflow <- new.env(parent = emptyenv())
        workflow$state_manager <- recordedQcSequenceManager(context)
        if (diagnostic) utils::Rprof(file.path(run_dir, "qc-sequence.Rprof"), interval = 0.02,
            memory.profiling = TRUE, gc.profiling = TRUE)
        for (name in names(context$settings)) {
            runRecordedQcSequenceStep(context, workflow, name)
            releaseRecordedQcSequence(context)
            receipt$completed_steps <- c(receipt$completed_steps, name)
            receipt$observations[[name]] <- observeRecordedQcSequence(context, name)
            flush_receipt()
        }
        if (diagnostic) utils::Rprof(NULL)
        context$metrics$observe("before:revert")
        invisible(workflow$state_manager$revertToState("raw_data_s4"))
        releaseRecordedQcSequence(context)
        invisible(workflow$state_manager$getState("precursor_rollup"))
        releaseRecordedQcSequence(context)
        context$native$state_history <- as.list(context$names)
        workflow$state_manager$current_state <- "precursor_rollup"
        context$native$current_state <- "precursor_rollup"
        context$metrics$observe("after:revert")
        context$metrics$observe("before:export")
        records <- stats::setNames(lapply(context$names, function(name) {
            recordedQcSequenceResolve(context, name)
        }), context$names)
        final <- records$precursor_rollup$data
        session <- list(r6_complete_states = records, r6_state_history = as.list(context$names),
            r6_current_state_name = "precursor_rollup", current_s4_object = final,
            workflow_type = "DIA", qc_params = workflow$qc_params,
            normalization_method = "not_run", ruv_mode = "not_run",
            final_protein_count = api(".countPeptideProteinGroups")(final),
            final_sample_count = nrow(methods::slot(final, "design_matrix")),
            scope = "three-step peptide QC; not a DA-ready normalized session")
        export_dir <- file.path(run_dir, "export")
        dir.create(export_dir)
        export_time <- context$settings$precursor_rollup$timestamp
        files <- checkpointUtilityTimed(context$metrics, "legacy_export", function() {
            api("saveProtNormExportArtifacts")(session, export_dir, timeFn = function() export_time)
        })
        records <- final <- session <- NULL
        context$metrics$observe("after:export")
        receipt$files <- recordedQcSequenceInventory(run_dir)
        receipt$status <- "complete"
        receipt$state_names <- context$names
        receipt$export <- files
        receipt$worker_counts <- as.integer(names(context$clusters))
        receipt$phase_sampling_note <- "100 ms controller target, actual sample gaps retained; phase maxima sampled, service peak exact"
        receipt
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

#' Require the three production plot inventories and all native export artifacts.
#' @param run_dir Member directory.
#' @return Content fingerprints; missing required writes are fatal.
#' @noRd
recordedQcSequenceInventory <- function(run_dir) {
    prefixes <- c("2_qval_filtered", "3_protein_peptide_filtered", "4_precursor_rollup")
    names <- c("proteins_total", "proteins_per_run", "peptides_total", "peptides_per_protein",
        "peptides_per_run", "combined_plots")
    plots <- unlist(lapply(prefixes, function(prefix) outer(paste0(prefix, "_", names),
        c("png", "pdf"), paste, sep = ".")), use.names = FALSE)
    paths <- c(file.path("plots", plots), file.path("export", "filtered_session_data_latest.rds"),
        file.path("export", "qc_params.RDS"), file.path("export", "filtered_session_summary.txt"),
        file.path("export", list.files(file.path(run_dir, "export"),
            pattern = "^filtered_session_data_[0-9]+_[0-9]+[.]rds$")))
    if (length(paths) != 40L || any(!file.exists(file.path(run_dir, paths))) ||
        any(file.size(file.path(run_dir, paths)) <= 0)) stop("required sequential QC output is missing")
    stats::setNames(lapply(paths, function(path) list(bytes = file.size(file.path(run_dir, path)),
        sha256 = digest::digest(file = file.path(run_dir, path), algo = "sha256"))), paths)
}

#' Verify newly calculated states against the independent recorded oracle.
#' @param library_path,run_dir,member_dir Installed library, verification output and member.
#' @return Exact-state/file verification; no throughput or GUI acceptance.
#' @noRd
verifyRecordedQcSequence <- function(library_path, run_dir, member_dir) {
    api <- recordedCheckpointApi(library_path)
    receipt <- jsonlite::read_json(file.path(member_dir, "qc-sequence-member.json"))
    if (!identical(receipt$status, "complete") || receipt$open_resources != 0L)
        stop("sequence did not complete with released resources")
    actual_files <- recordedQcSequenceInventory(member_dir)
    if (!identical(names(actual_files), names(receipt$files)) ||
        any(!vapply(names(actual_files), function(name) {
            identical(actual_files[[name]]$sha256, receipt$files[[name]]$sha256) &&
                identical(as.numeric(actual_files[[name]]$bytes), as.numeric(receipt$files[[name]]$bytes))
        }, logical(1)))) stop("sequence output fingerprints differ")
    manifest <- readRecordedCheckpointManifest(receipt$prepared_dir)
    session <- readRDS(file.path(member_dir, "export", "filtered_session_data_latest.rds"))
    names <- c("initial", "raw_data_s4", "qvalue_filtered", "protein_peptide_filtered", "precursor_rollup")
    if (!identical(names(session$r6_complete_states), names) ||
        !identical(unlist(session$r6_state_history), names)) stop("sequence history topology differs")
    comparable <- function(record) {
        if (isS4(record$data)) {
            args <- methods::slot(record$data, "args")
            args$peptide_qc_audit <- api(".peptideQcSubstantiveAudit")(record$data)
            methods::slot(record$data, "args") <- args
        }
        if (is.list(record$audit_metadata)) record$audit_metadata$timestamp_utc <- NULL
        record
    }
    for (name in names) {
        actual <- comparable(session$r6_complete_states[[name]])
        expected <- comparable(readRecordedCheckpoint(receipt$prepared_dir, manifest, name))
        if (!identical(actual, expected)) stop(paste("sequential state differs:", name,
            paste(utils::capture.output(all.equal(actual, expected)), collapse = "; ")))
    }
    if (!identical(session$current_s4_object, session$r6_complete_states$precursor_rollup$data))
        stop("sequence current state differs")
    timestamps <- list.files(file.path(member_dir, "export"),
        pattern = "^filtered_session_data_[0-9]+_[0-9]+[.]rds$", full.names = TRUE)
    if (!identical(readRDS(timestamps), session) ||
        !identical(readRDS(file.path(member_dir, "export", "qc_params.RDS")), session$qc_params))
        stop("sequence legacy export content differs")
    result <- list(status = "passed_real_qc_sequence", states = names, required_files = 40L,
        full_attributes = TRUE, substantive_audit = TRUE, promotion_authorized = FALSE)
    jsonlite::write_json(result, file.path(run_dir, "qc-sequence-verification.json"),
        auto_unbox = TRUE, pretty = TRUE)
    result
}

#' Observe the complete installed Q-value Apply operation with only active inputs.
#' @param library_path,run_dir Installed library and fresh diagnostic directory.
#' @param prepared_dir Accepted individual state fixtures.
#' @param memory_fn Existing guarded-service observation function.
#' @return Scalar observations and output reference, never a performance gate.
#' @noRd
executeRecordedActiveQc <- function(library_path, run_dir, prepared_dir, memory_fn) {
    entry <- memory_fn()
    api <- recordedCheckpointApi(library_path)
    manifest <- readRecordedCheckpointManifest(prepared_dir)
    raw_record <- readRecordedCheckpoint(prepared_dir, manifest, "raw_data_s4")
    config_path <- file.path(prepared_dir, "qvalue-config.rds")
    if (!identical(digest::digest(file = config_path, algo = "sha256"),
        manifest$qvalue_config_sha256)) stop("Q-value config integrity differs")
    settings <- readRDS(config_path)
    config <- settings$config
    required <- c("qvalue_threshold", "global_qvalue_threshold", "proteotypic_only",
        "global_pg_qvalue_threshold", "confidence_provenance_mode")
    if (!all(required %in% names(config))) stop("recorded Q-value settings are incomplete")
    manager <- api("WorkflowState")$new()
    manager$states <- list(initial = readRecordedCheckpoint(prepared_dir, manifest, "initial"),
        raw_data_s4 = raw_record)
    manager$state_history <- as.list(names(manager$states))
    manager$current_state <- "raw_data_s4"
    config_env <- new.env(parent = emptyenv())
    config_env$config_list <- methods::slot(raw_record$data, "args")
    raw_record <- NULL
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    observations <- list(entry = entry)
    observe <- function(label, extra = NULL) {
        observations[[label]] <<- list(service = memory_fn(),
            unique_owner_bytes = as.numeric(lobstr::obj_size(manager$states,
                config_env$config_list, extra)))
        invisible(NULL)
    }
    observe("active_input")
    allocation_path <- file.path(run_dir, "active-qc.Rprofmem")
    profile_path <- file.path(run_dir, "active-qc.Rprof")
    if (!isTRUE(capabilities("profmem"))) stop("allocation profiling is unavailable")
    apply_step <- api("runPeptideQvalueApplyStep")
    utils::Rprofmem(allocation_path, threshold = 100000)
    utils::Rprof(profile_path, interval = 0.02, memory.profiling = TRUE, gc.profiling = TRUE)
    on.exit({ utils::Rprofmem(NULL); utils::Rprof(NULL) }, add = TRUE)
    elapsed <- system.time(result <- apply_step(
        workflowData = workflow, qvalueThreshold = config$qvalue_threshold,
        globalQvalueThreshold = config$global_qvalue_threshold,
        proteotypicOnly = config$proteotypic_only,
        globalPGQvalueThreshold = config$global_pg_qvalue_threshold,
        confidenceProvenanceMode = config$confidence_provenance_mode,
        updateConfigParameterFn = function(...) {
            api("updateConfigParameter")(..., env = config_env)
        },
        qvalueFilterFn = function(theObject) {
            observe("method_entry", theObject)
            value <- api("srlQvalueProteotypicPeptideClean")(theObject)
            observe("method_return_before_audit_commit", list(theObject, value))
            value
        }, nowFn = \() settings$timestamp,
        logInfoFn = \(...) invisible(NULL), logWarnFn = \(...) invisible(NULL)),
        gcFirst = FALSE)
    utils::Rprofmem(NULL)
    utils::Rprof(NULL)
    observe("after_audit_commit", result$filteredS4)
    output <- file.path(run_dir, "qvalue-output.rds")
    saveRDS(result$filteredS4, output, compress = FALSE, version = 3)
    result <- NULL
    manager$states <- list()
    workflow$state_manager <- NULL
    config_env$config_list <- NULL
    manager <- NULL
    released_heap <- checkpointUtilityHeapBytes()
    observations$after_release <- list(service = memory_fn(), live_heap_bytes = released_heap)
    receipt <- list(kind = "recorded_active_qvalue_diagnostic", status = "complete",
        prepared_dir = normalizePath(prepared_dir), source_sha256 = manifest$source_sha256,
        input_manifest_sha256 = digest::digest(file = file.path(prepared_dir,
            "recorded-inputs.json"), algo = "sha256"), observations = observations,
        operation_elapsed_seconds = unname(elapsed[["elapsed"]]),
        operation_cpu_seconds = unname(sum(elapsed[c("user.self", "sys.self")])),
        output = output, output_sha256 = digest::digest(file = output, algo = "sha256"),
        allocation_profile = allocation_path, cpu_profile = profile_path,
        profiling = TRUE, independently_verified = FALSE, promotion_authorized = FALSE,
        limitation = "allocation totals are not peak residency; GUI, plots and later states excluded")
    jsonlite::write_json(receipt, file.path(run_dir, "recorded-qc.json"),
        auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
    receipt
}

#' Compare the complete QC S4 while excluding only documented audit timestamps.
#' @param library_path,run_dir Installed library and fresh verifier directory.
#' @param member_dir Completed active-QC diagnostic.
#' @return Scientific-content verification for this one recorded operation.
#' @noRd
verifyRecordedActiveQc <- function(library_path, run_dir, member_dir) {
    api <- recordedCheckpointApi(library_path)
    member <- jsonlite::read_json(file.path(member_dir, "recorded-qc.json"))
    if (!identical(member$status, "complete") || !identical(
        digest::digest(file = member$output, algo = "sha256"), member$output_sha256)) {
        stop("active QC output is incomplete or changed")
    }
    manifest <- readRecordedCheckpointManifest(member$prepared_dir)
    expected <- readRecordedCheckpoint(member$prepared_dir, manifest, "qvalue_filtered")$data
    actual <- readRDS(member$output)
    comparable <- function(object) {
        args <- methods::slot(object, "args")
        args$peptide_qc_audit <- api(".peptideQcSubstantiveAudit")(object)
        methods::slot(object, "args") <- args
        object
    }
    methods::validObject(actual)
    if (!identical(comparable(actual), comparable(expected))) {
        differences <- utils::capture.output(all.equal(comparable(actual), comparable(expected)))
        stop(paste("recorded Q-value result differs:", paste(differences, collapse = " ")))
    }
    result <- list(status = "passed_recorded_active_qc", member_dir = normalizePath(member_dir),
        exact_scientific_content = TRUE, excluded = "audit timestamp_utc only",
        promotion_authorized = FALSE)
    jsonlite::write_json(result, file.path(run_dir, "recorded-verification.json"),
        auto_unbox = TRUE, pretty = TRUE)
    result
}

#' Execute small installed DIA QC resident/RDS/product/asset scenarios.
#' @param library_path Installed candidate library.
#' @param run_dir Fresh request-owned output directory.
#' @param contract_path Frozen recipe contract.
#' @return Bound descriptive receipt with actual producer/renderer/decoder counters.
#' @noRd
runCheckpointQcConsumerScenarios <- function(library_path, run_dir, contract_path) {
    library_path <- normalizePath(library_path, mustWork = TRUE)
    run_dir <- normalizePath(run_dir, mustWork = TRUE)
    contract_path <- normalizePath(contract_path, mustWork = TRUE)
    .libPaths(unique(c(library_path, .libPaths())))
    namespace <- loadNamespace("MultiScholaR", lib.loc = library_path)
    if (!identical(normalizePath(getNamespaceInfo(namespace, "path")),
        normalizePath(file.path(library_path, "MultiScholaR")))) {
        stop("consumer scenarios loaded a different installed package")
    }
    api <- \(name) get(name, envir = namespace, inherits = FALSE)
    source_root <- normalizePath(file.path(dirname(contract_path), "..", "..", ".."))
    helpers <- new.env(parent = baseenv())
    sys.source(file.path(source_root, "tests", "testthat",
        "helper-workflow-checkpoint-fixtures.R"), envir = helpers)
    contract <- jsonlite::read_json(contract_path, simplifyVector = FALSE)
    recipe <- Filter(\(item) identical(item$id, "small"), contract$recipes)[[1L]]
    object <- helpers$makeWorkflowCheckpointProtein(recipe)
    aesthetics <- list(color_var = "group", shape_var = "batch")
    software <- paste0("MultiScholaR-", as.character(utils::packageVersion("MultiScholaR")))
    producer_calls <- 0L
    renderer_calls <- 0L
    produce <- function(type) {
        producer_calls <<- producer_calls + 1L
        api("buildProtNormQcPlotProduct")(object, aesthetics, type)
    }
    save_plot <- function(...) {
        renderer_calls <<- renderer_calls + 1L
        ggplot2::ggsave(...)
    }
    style <- list(title = "DIA QC consumer proof", font_size = 8)
    dimensions <- c(800L, 600L)
    store <- api("newWorkflowEntityStore")(run_dir, "qc-consumers", "DIA",
        max_cache_bytes = 67108864)
    events <- list()
    results <- list()
    for (type in c("pca", "rle", "density", "correlation")) {
        resident <- checkpointConsumerOperation(paste0(type, ":resident_produce_render"), \() {
            product <- produce(type)
            api("renderWorkflowDiaQcProduct")(product, style, "png", dimensions,
                savePlotFn = save_plot)
        })
        events[[resident$label]] <- resident[c("elapsed_seconds", "cpu_seconds")]
        rds_path <- file.path(run_dir, paste0(type, "-source.rds"))
        complete_rds <- checkpointConsumerOperation(paste0(type, ":complete_rds"), \() {
            saveRDS(object, rds_path)
            restored <- readRDS(rds_path)
            product <- api("buildProtNormQcPlotProduct")(restored, aesthetics, type)
            producer_calls <<- producer_calls + 1L
            api("renderWorkflowDiaQcProduct")(product, style, "png", dimensions,
                savePlotFn = save_plot)
        })
        events[[complete_rds$label]] <- complete_rds[c("elapsed_seconds", "cpu_seconds")]
        identity <- api("newWorkflowDiaQcProductIdentity")(object,
            "post_normalization", type, aesthetics, software)
        creation <- checkpointConsumerOperation(paste0(type, ":product_create_write"), \() {
            product <- produce(type)
            reference <- api("saveWorkflowDiaQcProduct")(store, product,
                paste0(type, "-product"), identity, "post_normalization")
            list(product = product, reference = reference)
        })
        events[[creation$label]] <- creation[c("elapsed_seconds", "cpu_seconds")]
        api("releaseWorkflowEntityCache")(store)
        product_hit <- checkpointConsumerOperation(paste0(type, ":product_lookup"), \() {
            api("resolveWorkflowDiaQcProduct")(store, identity,
                "post_normalization", type)
        })
        events[[product_hit$label]] <- product_hit[c("elapsed_seconds", "cpu_seconds")]
        render_identity <- api("newSavedWorkflowRenderIdentity")(
            creation$value$reference, style, "png", dimensions, "ggplot2-consumer-v1")
        asset_creation <- checkpointConsumerOperation(paste0(type, ":asset_create_write"), \() {
            asset <- api("renderWorkflowDiaQcProduct")(product_hit$value$value$value,
                style, "png", dimensions, savePlotFn = save_plot)
            reference <- api("saveWorkflowDiaQcAsset")(store, asset,
                paste0(type, "-asset"), render_identity)
            list(asset = asset, reference = reference)
        })
        events[[asset_creation$label]] <- asset_creation[c("elapsed_seconds", "cpu_seconds")]
        api("releaseWorkflowEntityCache")(store)
        before_asset <- api("workflowEntityCacheInfo")(store)
        asset_hit <- checkpointConsumerOperation(paste0(type, ":asset_lookup"), \() {
            api("resolveSavedWorkflowAsset")(store, creation$value$reference,
                render_identity)
        })
        events[[asset_hit$label]] <- asset_hit[c("elapsed_seconds", "cpu_seconds")]
        after_asset <- api("workflowEntityCacheInfo")(store)
        stale <- identity
        stale$software <- paste0(stale$software, "-changed")
        miss <- api("resolveWorkflowDiaQcProduct")(store, stale,
            "post_normalization", type)
        if (!identical(product_hit$value$status, "hit") ||
            !identical(asset_hit$value$status, "hit") ||
            !identical(miss$status, "miss") ||
            after_asset$payload_reads != before_asset$payload_reads + 1L) {
            stop("QC consumer scenario hit/miss counters disagree")
        }
        results[[type]] <- list(product_entity = creation$value$reference$entity_id,
            asset_entity = asset_creation$value$reference$entity_id,
            product_bytes = creation$value$reference$payload$bytes,
            asset_bytes = asset_creation$value$reference$payload$bytes,
            source_rds_bytes = file.size(rds_path),
            asset_content_bytes = length(asset_hit$value$value$bytes))
    }
    counters <- api("workflowEntityCacheInfo")(store)
    api("closeWorkflowEntityStore")(store)
    receipt <- list(status = "passed_qc_consumer_scenarios", route = "iq",
        recipe = "small", package_path = normalizePath(getNamespaceInfo(namespace, "path")),
        contract_sha256 = digest::digest(file = contract_path, algo = "sha256"),
        producer_calls = producer_calls, renderer_calls = renderer_calls,
        payload_reads = counters$payload_reads, payload_decodes = counters$payload_decodes,
        cache_hits = counters$cache_hits, events = events, results = results,
        includes_creation_write_and_miss = TRUE, performance_claim = FALSE,
        promotion_authorized = FALSE)
    jsonlite::write_json(receipt, file.path(run_dir, "qc-consumers.json"),
        auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
    receipt
}

#' Execute small installed DIA fit/result/view reuse scenarios.
#' @param library_path Installed candidate library.
#' @param run_dir Fresh request-owned output directory.
#' @param contract_path Frozen recipe contract.
#' @return Descriptive receipt with separate statistical-operation counters.
#' @noRd
runCheckpointDaConsumerScenarios <- function(library_path, run_dir, contract_path) {
    library_path <- normalizePath(library_path, mustWork = TRUE)
    run_dir <- normalizePath(run_dir, mustWork = TRUE)
    contract_path <- normalizePath(contract_path, mustWork = TRUE)
    .libPaths(unique(c(library_path, .libPaths())))
    namespace <- loadNamespace("MultiScholaR", lib.loc = library_path)
    if (!identical(normalizePath(getNamespaceInfo(namespace, "path")),
        normalizePath(file.path(library_path, "MultiScholaR")))) {
        stop("DA consumer scenarios loaded a different installed package")
    }
    api <- \(name) get(name, envir = namespace, inherits = FALSE)
    source_root <- normalizePath(file.path(dirname(contract_path), "..", "..", ".."))
    helpers <- new.env(parent = baseenv())
    sys.source(file.path(source_root, "tests", "testthat",
        "helper-workflow-checkpoint-fixtures.R"), envir = helpers)
    contract <- jsonlite::read_json(contract_path, simplifyVector = FALSE)
    recipe <- Filter(\(item) identical(item$id, "small"), contract$recipes)[[1L]]
    object <- helpers$makeWorkflowCheckpointProtein(recipe)
    quant <- methods::slot(object, "protein_quant_table")
    key <- methods::slot(object, "protein_id_column")
    design <- methods::slot(object, "design_matrix")
    sample_id <- methods::slot(object, "sample_id")
    data <- as.matrix(quant[setdiff(names(quant), key)])
    rownames(data) <- quant[[key]]
    rownames(design) <- as.character(design[[sample_id]])
    formula <- stats::as.formula("~0+group")
    model <- stats::model.matrix(formula, stats::model.frame(formula, design))
    data <- data[, rownames(model), drop = FALSE]
    contrast <- "treated_vs_control=grouptreated-groupcontrol"
    counts <- new.env(parent = emptyenv())
    count <- function(operation) {
        assign(operation, get0(operation, envir = counts, ifnotfound = 0L) + 1L,
            envir = counts)
    }
    contrast_matrix <- limma::makeContrasts(contrasts = contrast,
        levels = colnames(model))
    events <- list()
    resident <- checkpointConsumerOperation("resident_fit_contrast_results", \() {
        base <- api("fitProtDaBaseModel")(data, model, design, rownames(model),
            counterFn = count, messageFn = \(...) invisible(NULL))
        contrasted <- api("applyProtDaContrastModel")(base, contrast_matrix,
            counterFn = count, messageFn = \(...) invisible(NULL))
        results <- api("buildProtDaResultTables")(contrasted, contrast,
            counterFn = count, messageFn = \(...) invisible(NULL))
        list(base = base, contrasted = contrasted, results = results)
    })
    events[[resident$label]] <- resident[c("elapsed_seconds", "cpu_seconds")]
    base <- resident$value$base
    results <- resident$value$results
    software <- paste0("MultiScholaR-", as.character(utils::packageVersion("MultiScholaR")))
    base_identity <- api("newWorkflowDiaBaseFitIdentity")(data, design,
        "~0+group", base$block, NA, "fixture-none", "ruv-skip", software)
    store <- api("newWorkflowEntityStore")(run_dir, "da-consumers", "DIA",
        max_cache_bytes = 67108864)
    base_write <- checkpointConsumerOperation("base_fit_write", \() {
        api("saveWorkflowDiaBaseFit")(store, base, "base-fit", base_identity)
    })
    events[[base_write$label]] <- base_write[c("elapsed_seconds", "cpu_seconds")]
    base_ref <- base_write$value
    result_identity <- api("newWorkflowDiaResultIdentity")(base_ref, contrast,
        NA_real_, FALSE, FALSE, "qvalue-then-BH-v1", software)
    result_write <- checkpointConsumerOperation("full_result_write", \() {
        api("saveWorkflowDiaResult")(store, results, "full-result", result_identity)
    })
    events[[result_write$label]] <- result_write[c("elapsed_seconds", "cpu_seconds")]
    result_ref <- result_write$value
    api("releaseWorkflowEntityCache")(store)
    counts_before_hit <- as.list(counts)
    result_lookup <- checkpointConsumerOperation("full_result_lookup", \() {
        api("resolveWorkflowDiaResult")(store, result_identity)
    })
    events[[result_lookup$label]] <- result_lookup[c("elapsed_seconds", "cpu_seconds")]
    result_hit <- result_lookup$value
    if (!identical(as.list(counts), counts_before_hit) ||
        !identical(result_hit$value$value$results, results$results)) {
        stop("full DA result hit recomputed or changed statistics")
    }
    reverse <- "control_vs_treated=groupcontrol-grouptreated"
    reverse_matrix <- limma::makeContrasts(contrasts = reverse,
        levels = colnames(model))
    base_reuse <- checkpointConsumerOperation("base_fit_lookup_changed_contrast", \() {
        base_hit <- api("resolveWorkflowDiaBaseFit")(store, base_identity,
            base$uses_duplicate_correlation, base$consensus_correlation)
        restored_base <- base
        restored_base$fit <- base_hit$value$fit
        reversed <- api("applyProtDaContrastModel")(restored_base, reverse_matrix,
            counterFn = count, messageFn = \(...) invisible(NULL))
        api("buildProtDaResultTables")(reversed, reverse,
            counterFn = count, messageFn = \(...) invisible(NULL))
    })
    events[[base_reuse$label]] <- base_reuse[c("elapsed_seconds", "cpu_seconds")]
    reverse_results <- base_reuse$value
    if (get0("lmFit", envir = counts, ifnotfound = 0L) != 1L) {
        stop("base-fit hit unexpectedly refit the model")
    }
    envelope <- api("newSavedWorkflowResultTable")(results$results[[1L]],
        ".rownames", contrast, "fdr_qvalue")
    views <- checkpointConsumerOperation("table_volcano_heatmap_views", \() {
        list(page = api("pageWorkflowDiaResult")(envelope, 1L, 25L),
            volcano = api("filterWorkflowDiaVolcano")(envelope, 0.05, 1),
            heatmap = api("selectWorkflowDiaHeatmap")(envelope, data, 50L,
                scaling = "row", clustering = "both"))
    })
    events[[views$label]] <- views[c("elapsed_seconds", "cpu_seconds")]
    page <- views$value$page
    volcano <- views$value$volcano
    heatmap <- views$value$heatmap
    changed_design <- design
    changed_design$group <- rev(changed_design$group)
    changed_identity <- api("newWorkflowDiaBaseFitIdentity")(data, changed_design,
        "~0+group", base$block, NA, "fixture-none", "ruv-skip", software)
    if (!identical(api("resolveWorkflowDiaBaseFit")(store, changed_identity,
        base$uses_duplicate_correlation, base$consensus_correlation)$status, "miss")) {
        stop("changed DA design reused an incompatible fit")
    }
    s4_path <- file.path(run_dir, "source-s4.rds")
    s4_rds <- checkpointConsumerOperation("complete_s4_rds_write_read", \() {
        saveRDS(object, s4_path)
        readRDS(s4_path)
    })
    events[[s4_rds$label]] <- s4_rds[c("elapsed_seconds", "cpu_seconds")]
    restored_object <- s4_rds$value
    methods::validObject(restored_object)
    projection <- tryCatch({
        api("hydrateSavedWorkflowProjection")(store, result_ref, c("logFC"))
        "unexpected"
    }, saved_entity_unsupported = \(error) "unsupported")
    cache <- api("workflowEntityCacheInfo")(store)
    api("closeWorkflowEntityStore")(store)
    receipt <- list(status = "passed_da_consumer_scenarios", route = "iq",
        recipe = "small", package_path = normalizePath(getNamespaceInfo(namespace, "path")),
        contract_sha256 = digest::digest(file = contract_path, algo = "sha256"),
        counters = as.list(counts), base_fit_bytes = base_ref$payload$bytes,
        result_bytes = result_ref$payload$bytes, source_s4_bytes = file.size(s4_path),
        events = events,
        result_hit_without_statistics = TRUE, changed_design_miss = TRUE,
        reverse_contrast_rows = nrow(reverse_results$results[[1L]]),
        page_rows = nrow(page$data), full_population_rows = page$full_population_rows,
        volcano_rows = nrow(volcano$data), heatmap_rows = nrow(heatmap$matrix),
        heatmap_samples = ncol(heatmap$matrix), projection = projection,
        payload_reads = cache$payload_reads, payload_decodes = cache$payload_decodes,
        performance_claim = FALSE, promotion_authorized = FALSE)
    jsonlite::write_json(receipt, file.path(run_dir, "da-consumers.json"),
        auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
    receipt
}

#' Build a deterministic DIA-shaped DA table for access-path comparison.
#' @param rows,samples Declared protein and sample dimensions.
#' @return Table with matching normalized/raw sample columns.
#' @noRd
checkpointProjectionTable <- function(rows = 8020L, samples = 24L) {
    set.seed(13013L)
    sample_ids <- sprintf("run_%04d", seq_len(samples))
    table <- data.frame(Protein.Group = sprintf("PG_%06d", seq_len(rows)),
        logFC = stats::rnorm(rows), raw_pvalue = stats::runif(rows),
        stringsAsFactors = FALSE)
    table$fdr_qvalue <- stats::p.adjust(table$raw_pvalue, method = "BH")
    table$comparison <- factor("treated_vs_control",
        levels = c("treated_vs_control", "unused"))
    table$numerator <- "treated"
    table$denominator <- "control"
    table$uniprot_acc <- sprintf("P%05d", seq_len(rows))
    table$annotation <- paste0("synthetic_", seq_len(rows))
    for (sample in sample_ids) {
        table[[paste0("log2norm.", sample)]] <- stats::rnorm(rows, 20, 2)
    }
    for (sample in sample_ids) {
        table[[paste0("raw.", sample)]] <- 2^table[[paste0("log2norm.", sample)]]
    }
    stopifnot(ncol(table) == 9L + 2L * samples)
    rownames(table) <- table$Protein.Group
    attr(table, "provenance") <- "DIA reduced-control shape proxy; access only"
    table
}

#' Validate the fixed, profiling-only reader operations without loading tables.
#' @param operation,input_path,selection Runner operation, owned input and fixed case selector.
#' @return TRUE invisibly.
#' @noRd
validateCheckpointReaderInput <- function(operation, input_path, selection = NULL) {
    if (!operation %in% paste0("reader-", c("build", "check", "prepare", "inspect", "member", "verify")) ||
        !dir.exists(input_path)) stop("invalid reader operation or input directory")
    if (operation == "reader-check" && !is.null(selection) &&
        !selection %in% c("perf101:D0", "perf101:D1", "perf101:D0R", "perf101:V", "gui-sharing-repair", "gui-export-history-member", "gui-iq-history-member", "gui-ruv-pearson", "gui-native-session", "gui-native-outputs", "gui-reference-first-qc", "gui-reference-core", "gui-reference-routing", "gui-mechanism-preflight", "gui-rollup-profile", "gui-rollup-profile-preflight", "gui-raw-replay:small", "gui-raw-replay:outputs", "gui-raw-replay:verify", "gui-cohort-profile", "gui-impute-workers:preflight", "gui-impute-workers:4", "gui-impute-workers:8", "gui-cohort-size:256", "gui-cohort-size:768", "gui-cohort-size:768-threads8", "gui-cohort-threads:1", "gui-cohort-threads:2", "gui-cohort-threads:4", "gui-cohort-threads:8", "gui-raw-replay:300") &&
        !grepl("^gui-(rollup|session)-performance:(control|candidate)$", selection) &&
        !grepl("^gui-bulk-rollup:(qualification-preflight|qualification-prepare|writer-build|writer-opt-(build|preflight|diagnostic:(primary|contrast))|writer-r-preflight|writer-r-readback:[A-Za-z0-9]{6}|writer-r:(screen|full):(primary|contrast):P1[456]:[12]|io-preflight|io-prepare|io:(duckdb|fastduckdb):[12]|scan-inspection|crossover:(B|P0):(4608|9600|96000|960000):[12]|contrast:(B|P([0-9]|1[0-4]))|diagnostic|diagnostic-preflight|preflight|relational-preflight|materialize-preflight|adaptation-preflight|fast-p3-preflight|writer-preflight|screen:([ABCD]|P([0-9]|1[0-4])):[12]|full:([ABCD]|P([0-9]|1[0-4])))$", selection) &&
        !grepl("^gui-rollup-compare:(control|candidate)$", selection) &&
        !grepl("^gui-first-qc-verify:[A-Za-z0-9]{6}:[A-Za-z0-9]{6}$", selection) &&
        !grepl("^gui-history-verify:[A-Za-z0-9]{6}(:[A-Za-z0-9]{6})?$", selection) &&
        !selection %in% c("sequence", "gui-capacity", "gui-capacity-verifier", "gui-diagnostic", "gui-passive", "gui-lifetime-repair", "gui-memory-repair", "gui-memory-verification", "gui-audit-repair", "gui-audit-diagnostic", "gui-audit-member", "gui-audit-toggle", "gui-sleep-guard", "gui-integration-repair", "gui-worker-member", "gui-export-member"))
        stop("unsupported installed check selection")
    if (operation == "reader-build") {
        lock <- file.path(input_path, "Cargo.lock")
        if (!file.exists(lock) || !identical(digest::digest(file = lock, algo = "sha256"),
            "0dfffbc58e67f7c3eff793b6b7fdbe17ebe8733838e52297aebeef02dc02e073") ||
            !file.exists(file.path(input_path, "fastduckdb", "DESCRIPTION"))) {
            stop("reader source does not match the pinned upstream lockfile")
        }
    } else if (operation %in% c("reader-prepare", "reader-inspect")) {
        if (!file.exists(file.path(input_path, "recorded-inputs.json"))) stop("missing recorded input")
        if (!is.null(selection) && (operation != "reader-prepare" ||
            !identical(selection, "adapter"))) stop("unsupported reader preparation")
    } else if (operation == "reader-member") {
        if (!file.exists(file.path(input_path, "reader-inputs.json")) ||
            !is.character(selection) || length(selection) != 1L ||
            !grepl(paste0("^(small|wide|precursor):((rds|arrow|duckdb|fastduckdb):(1|5)|",
                "adapter_(rds|duckdb|fastduckdb):1)$"), selection)) {
            stop("invalid frozen reader selection")
        }
    } else if (operation == "reader-verify" &&
        !file.exists(file.path(input_path, "reader-member.json"))) stop("missing reader member")
    invisible(TRUE)
}

#' Install the pinned upstream reader into the explicit isolated proof library.
#' @param library_path,run_dir,source Installed library, owned run directory and source checkout.
#' @return Bound native build receipt.
#' @noRd
buildCheckpointReader <- function(library_path, run_dir, source) {
    if (dir.exists(file.path(library_path, "fastduckdb"))) stop("refusing to replace installed reader")
    lock <- file.path(source, "Cargo.lock")
    lock_hash <- digest::digest(file = lock, algo = "sha256")
    Sys.setenv(CARGO_BUILD_JOBS = "2", MAKEFLAGS = "-j2")
    status <- system2(file.path(R.home("bin"), "R"),
        c("CMD", "INSTALL", shQuote(paste0("--library=", library_path)),
            shQuote(file.path(source, "fastduckdb"))))
    if (status != 0L) stop("pinned fastduckdb installation failed; retain build log")
    if (!identical(digest::digest(file = lock, algo = "sha256"), lock_hash)) {
        stop("upstream dependency lockfile changed during build")
    }
    artifacts <- list.files(file.path(library_path, "fastduckdb"),
        pattern = "(\\.so$|\\.rdb$|\\.rdx$|^DESCRIPTION$|^NAMESPACE$)",
        recursive = TRUE, full.names = TRUE)
    receipt <- list(status = "built", commit = "d78511504fbdb89fbf34bc26b22b3ff1524b0555",
        source = normalizePath(source), lock_sha256 = lock_hash,
        rustc = system2("rustc", "--version", stdout = TRUE),
        cargo = system2("cargo", "--version", stdout = TRUE),
        build_jobs = 2L, R = R.version.string, library = normalizePath(library_path),
        artifacts = as.list(stats::setNames(vapply(artifacts,
            \(path) digest::digest(file = path, algo = "sha256"), character(1)),
            substring(artifacts, nchar(library_path) + 2L))))
    jsonlite::write_json(receipt, file.path(run_dir, "reader-build.json"),
        auto_unbox = TRUE, pretty = TRUE)
    receipt
}

#' Persist one reader input without changing the application's storage contract.
#' @param table,case,run_dir,api Original table, fixed case name, owned directory and installed helpers.
#' @param adapter Use the installed store's explicit metadata-v2 writer.
#' @return Small input binding and measured producer costs.
#' @noRd
prepareCheckpointReaderTable <- function(table, case, run_dir, api, adapter = FALSE) {
    directory <- file.path(run_dir, case)
    if (!dir.create(directory)) stop("reader input directory already exists")
    large <- identical(case, "precursor")
    plain <- as.data.frame(table)
    metadata <- api("workflowSavedTableMetadata")(if (large) plain[0, , drop = FALSE] else plain,
        ".rownames")
    metadata$rows <- nrow(table)
    # Native I/O-only input: keep automatic row names compact, not a new store codec.
    row_names <- .row_names_info(table, type = 0L)
    metadata$row_names <- NULL
    named_columns <- names(table)[vapply(table, \(column) !is.null(names(column)), logical(1))]
    name_fields <- stats::setNames(paste0(".multischolar_names_", seq_along(named_columns)),
        named_columns)
    if (!length(named_columns)) name_fields <- stats::setNames(character(), character())
    restoration <- list(metadata = metadata, row_names = row_names,
        table_class = class(table), name_fields = name_fields)
    if (as.numeric(utils::object.size(restoration)) > 8 * 1024^2) stop("reader metadata exceeds 8 MiB")
    metadata_path <- file.path(directory, "metadata.rds")
    events <- list()
    rds_path <- file.path(directory, "table.rds")
    events$rds <- checkpointConsumerOperation("rds_write_hash", function() {
        saveRDS(table, rds_path)
        digest::digest(file = rds_path, algo = "sha256")
    })
    parquet_path <- if (adapter) file.path(directory,
        ".workflow-entities/entities/adapter-table/payload.parquet") else file.path(directory, "table.parquet")
    reference <- NULL
    events$parquet <- checkpointConsumerOperation("parquet_write_hash", function() {
        if (adapter) {
            store <- api("newWorkflowEntityStore")(directory, paste0("adapter-", case), "DIA")
            on.exit(api("closeWorkflowEntityStore")(store), add = TRUE)
            identity <- api("newSavedWorkflowScientificIdentity")(
                events$rds$value, "frozen-DIA-reader-shape", "metadata-v2",
                paste0("adapter-", case), "MultiScholaR-installed-adapter-v1")
            reference <<- api("saveWorkflowTable")(store, table, "adapter-table", identity,
                ".rownames", metadata_version = 2L)
            return(reference$payload$sha256)
        } else if (large || length(name_fields) || !identical(class(table), "data.frame")) {
            storage <- plain
            storage$.multischolar_order <- seq_len(nrow(table))
            for (column in names(name_fields)) storage[[name_fields[[column]]]] <- names(table[[column]])
            arrow::write_parquet(storage, parquet_path, compression = "snappy")
        } else api("writeSavedWorkflowParquet")(table, ".rownames", parquet_path)
        digest::digest(file = parquet_path, algo = "sha256")
    })
    if (adapter) restoration$reference <- reference
    saveRDS(restoration, metadata_path, compress = FALSE)
    files <- c(rds = rds_path, parquet = parquet_path, metadata = metadata_path)
    list(case = case, rows = nrow(table), columns = names(table),
        directory = normalizePath(directory), native_io_only = large && !adapter,
        store_adapter = adapter,
        files = Map(function(path, name) list(name = if (adapter && name == "parquet")
            ".workflow-entities/entities/adapter-table/payload.parquet" else basename(path), bytes = file.size(path),
            sha256 = if (name == "metadata") digest::digest(file = path, algo = "sha256") else
                events[[name]]$value), files, names(files)),
        producer = lapply(events, \(event) event[c("elapsed_seconds", "cpu_seconds")]))
}

#' Prepare the three approved shapes once in an independent guarded process.
#' @param library_path,run_dir,prepared_dir Installed build, owned output and recorded source.
#' @param adapter Explicit installed-store comparison, not the historical native-I/O trial.
#' @return Immutable trial manifest.
#' @noRd
prepareCheckpointReader <- function(library_path, run_dir, prepared_dir, adapter = FALSE) {
    api <- recordedCheckpointApi(library_path)
    manifest <- readRecordedCheckpointManifest(prepared_dir)
    if (!identical(manifest$source_sha256,
        "4671bdfab031cd880ab2bd2b414c4625cc30902a247e7d60f1c8dc11411a50c3")) {
        stop("recorded source differs from approved trial")
    }
    cases <- list()
    for (case in c("small", "wide", "precursor")) {
        table <- switch(case, small = checkpointProjectionTable(),
            wide = checkpointProjectionTable(10000L, 48L), precursor = {
                record <- readRecordedCheckpoint(prepared_dir, manifest, "raw_data_s4")
                value <- methods::slot(record$data, "peptide_data")
                if (!identical(dim(value), c(960000L, 14L))) stop("precursor shape changed")
                rm(record)
                value
            })
        cases[[case]] <- prepareCheckpointReaderTable(table, case, run_dir, api, adapter)
        rm(table)
        gc()
    }
    result <- list(status = "prepared", schema = if (adapter) "dia013_store_adapter_v1" else
        "dia013_reader_trial_v1", cases = cases,
        source_sha256 = manifest$source_sha256,
        settings = list(threads = 2L, lazy_strings = FALSE, compression = "snappy",
            page_rows = 25L, range_start = 480001L, range_end = 489600L),
        promotion_authorized = FALSE)
    jsonlite::write_json(result, file.path(run_dir, "reader-inputs.json"),
        auto_unbox = TRUE, pretty = TRUE, digits = NA)
    result
}

#' Validate a prepared case's paths and payloads once before opening it.
#' @param spec Case descriptor.
#' @param root Explicit owning preparation directory.
#' @param engine Which payload must be verified.
#' @return Named canonical paths.
#' @noRd
checkpointReaderPaths <- function(spec, root, engine) {
    adapter <- startsWith(engine, "adapter_")
    if (adapter && !isTRUE(spec$store_adapter)) stop("adapter requires store-backed inputs")
    engine <- sub("^adapter_", "", engine)
    if (!engine %in% c("rds", "arrow", "duckdb", "fastduckdb") ||
        !spec$case %in% c("small", "wide", "precursor", "test")) stop("unsupported reader case")
    directory <- normalizePath(file.path(root, spec$case), mustWork = TRUE)
    if (!identical(directory, spec$directory) ||
        !startsWith(directory, paste0(normalizePath(root), "/"))) stop("foreign reader input")
    payload <- if (engine == "rds") "rds" else "parquet"
    components <- if (adapter && engine == "rds") payload else c(payload, "metadata")
    paths <- vapply(components, function(name) {
        file <- spec$files[[name]]
        expected <- c(rds = "table.rds", parquet = "table.parquet", metadata = "metadata.rds")[[name]]
        if (isTRUE(spec$store_adapter) && name == "parquet") expected <-
            ".workflow-entities/entities/adapter-table/payload.parquet"
        if (!identical(file$name, expected)) stop("invalid reader payload name")
        path <- file.path(directory, expected)
        if (nzchar(Sys.readlink(path)) || !file.exists(path) ||
            !identical(as.numeric(file.size(path)), as.numeric(file$bytes)) ||
            !(adapter && name == "parquet") &&
                !identical(digest::digest(file = path, algo = "sha256"), file$sha256)) {
            stop("reader input integrity failure")
        }
        normalizePath(path)
    }, character(1))
    stats::setNames(paths, if (length(paths) == 1L) "payload" else c("payload", "metadata"))
}

#' Open one action with no resident result cache; never starts DuckDB in another arm.
#' @param spec,root,engine,api Frozen case, owner, reader and installed helpers.
#' @return Explicit action owner released by closeCheckpointReader.
#' @noRd
openCheckpointReader <- function(spec, root, engine, api) {
    paths <- checkpointReaderPaths(spec, root, engine)
    scope <- new.env(parent = emptyenv())
    scope$paths <- paths
    scope$adapter <- startsWith(engine, "adapter_")
    engine <- sub("^adapter_", "", engine)
    scope$engine <- engine
    scope$api <- api
    scope$restoration <- if (scope$adapter && engine == "rds")
        list(metadata = list(columns = unlist(spec$columns), rows = spec$rows)) else
            readRDS(paths[["metadata"]])
    scope$identity <- api("workflowTableFileIdentity")(paths[["payload"]])
    scope$connection <- NULL
    scope$table <- NULL
    scope$table_scope <- NULL
    scope$store <- NULL
    scope$read_metrics <- NULL
    scope$closed <- FALSE
    opened <- FALSE
    on.exit(if (!opened) closeCheckpointReader(scope), add = TRUE)
    if (scope$adapter && engine != "rds") {
        scope$store <- api("newWorkflowEntityStore")(spec$directory, paste0("adapter-", spec$case), "DIA")
        reference <- scope$restoration$reference
        if (!identical(reference$payload$sha256, spec$files$parquet$sha256)) {
            stop("adapter reference differs from prepared payload")
        }
        scope$table_scope <- api("newWorkflowTableReadScope")(scope$store, reference, reader = engine)
    } else if (engine == "duckdb") {
        scope$connection <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE), dbdir = ":memory:")
        DBI::dbExecute(scope$connection, "SET threads = 2")
    } else if (engine == "rds") scope$table <- readRDS(paths[["payload"]])
    else if (engine == "arrow") arrow::set_cpu_count(2L)
    opened <- TRUE
    scope
}

#' Release the action even when a source mutation or consumer fails.
#' @param scope Explicit action owner.
#' @return NULL invisibly.
#' @noRd
closeCheckpointReader <- function(scope) {
    if (scope$closed) return(invisible(NULL))
    if (!is.null(scope$connection)) DBI::dbDisconnect(scope$connection, shutdown = TRUE)
    problem <- NULL
    if (!is.null(scope$store)) {
        if (!is.null(scope$table_scope)) problem <- tryCatch({
            scope$api("releaseWorkflowTableReadScope")(scope$table_scope)
            NULL
        }, error = identity)
        scope$read_metrics <- scope$api("workflowEntityCacheInfo")(scope$store)
        scope$api("closeWorkflowEntityStore")(scope$store)
    }
    if (scope$engine == "fastduckdb" && !isTRUE(scope$adapter)) fastduckdb::clear_parquet_cache()
    scope$connection <- NULL
    scope$table <- NULL
    scope$restoration <- NULL
    scope$table_scope <- NULL
    scope$store <- NULL
    scope$closed <- TRUE
    if (inherits(problem, "error")) stop(problem)
    invisible(NULL)
}

#' Select from an already-resident native table without reconstructing its types.
#' @param table Original typed R table.
#' @param columns Requested unique columns.
#' @param rows NULL for the full population, otherwise one inclusive range.
#' @return Native selected table preserving scientific attributes and element names.
#' @noRd
subsetCheckpointAdapterTable <- function(table, columns, rows = NULL) {
    if (is.null(rows)) {
        result <- table[columns]
        for (extra in setdiff(names(attributes(table)), c("names", "row.names", "class"))) {
            attr(result, extra) <- attr(table, extra)
        }
        return(result)
    }
    selected <- if (rows[[1L]] > nrow(table)) integer() else
        seq.int(rows[[1L]], min(rows[[2L]], nrow(table)))
    result <- unclass(table[selected, columns, drop = FALSE])
    for (column in columns) {
        attrs <- attributes(table[[column]])
        if (!is.null(attrs$names)) attrs$names <- attrs$names[selected]
        attributes(result[[column]]) <- attrs
    }
    for (extra in setdiff(names(attributes(table)), c("names", "row.names", "class"))) {
        attr(result, extra) <- attr(table, extra)
    }
    class(result) <- class(table)
    result
}

#' Read only the fixed projection/range capability, restoring exact R metadata.
#' @param scope Action owner.
#' @param columns Requested unique column names.
#' @param rows Inclusive contiguous range, or NULL for all rows.
#' @return Fully materialized data frame; no arbitrary key or rank queries.
#' @noRd
readCheckpointReader <- function(scope, columns, rows = NULL) {
    if (scope$closed) stop("reader scope is closed")
    metadata <- scope$restoration$metadata
    if (!is.character(columns) || !length(columns) || anyNA(columns) ||
        anyDuplicated(columns) || !all(columns %in% metadata$columns) ||
        !is.null(rows) && (!is.numeric(rows) || length(rows) != 2L || anyNA(rows) ||
            any(!is.finite(rows)) || any(rows != floor(rows)) || rows[[1L]] < 1 ||
            rows[[2L]] < rows[[1L]])) stop("unsupported reader selection")
    identity <- scope$api("workflowTableFileIdentity")
    if (!identical(identity(scope$paths[["payload"]]), scope$identity)) stop("reader input changed")
    if (isTRUE(scope$adapter) && scope$engine == "rds") {
        return(subsetCheckpointAdapterTable(scope$table, columns, rows))
    }
    if (isTRUE(scope$adapter) && scope$engine != "rds") {
        offset <- if (is.null(rows)) 0L else rows[[1L]] - 1L
        limit <- if (is.null(rows)) metadata$rows else rows[[2L]] - rows[[1L]] + 1L
        return(scope$api("readSavedWorkflowTable")(scope$store,
            scope$restoration$reference, columns, offset = offset, limit = limit,
            scope = scope$table_scope)$data)
    }
    path <- scope$paths[["payload"]]
    name_fields <- scope$restoration$name_fields
    name_fields <- name_fields[names(name_fields) %in% columns]
    selected <- c(".multischolar_order", columns, unname(name_fields))
    data <- switch(scope$engine,
        rds = {
            value <- scope$table
            indices <- if (is.null(rows)) seq_len(nrow(value)) else
                which(seq_len(nrow(value)) >= rows[[1L]] & seq_len(nrow(value)) <= rows[[2L]])
            value <- value[indices, columns, drop = FALSE]
            value$.multischolar_order <- indices
            for (column in names(name_fields)) value[[name_fields[[column]]]] <-
                names(scope$table[[column]])[indices]
            value
        },
        fastduckdb = fastduckdb::read_parquet_fast(path, columns = selected,
            n_threads = 2L, lazy_strings = FALSE,
            filter = if (is.null(rows)) NULL else list(.multischolar_order = rows)),
        arrow = {
            dataset <- arrow::open_dataset(path, format = "parquet")
            query <- dplyr::select(dataset, dplyr::all_of(selected))
            if (!is.null(rows)) query <- dplyr::filter(query,
                .multischolar_order >= rows[[1L]], .multischolar_order <= rows[[2L]])
            as.data.frame(dplyr::collect(query))
        },
        duckdb = {
            quoted <- as.character(DBI::dbQuoteIdentifier(scope$connection, selected))
            sql <- paste0("SELECT ", paste(quoted, collapse = ", "), " FROM read_parquet(?)",
                if (is.null(rows)) "" else " WHERE \".multischolar_order\" BETWEEN ? AND ?",
                " ORDER BY \".multischolar_order\"")
            DBI::dbGetQuery(scope$connection, sql,
                params = if (is.null(rows)) list(path) else list(path, rows[[1L]], rows[[2L]]))
        })
    if (is.unsorted(data$.multischolar_order)) {
        data <- data[order(data$.multischolar_order), , drop = FALSE]
    }
    source_order <- as.integer(data$.multischolar_order)
    # Reuse the package's name restoration. This adapts the historical trial
    # sidecar only; it does not relabel those files as accepted store payloads.
    restoration <- metadata
    restoration$schema <- 2L
    restoration$name_columns <- name_fields
    restoration$class <- scope$restoration$table_class
    result <- unclass(scope$api("restoreSavedWorkflowTable")(data, restoration, columns))
    original_names <- scope$restoration$row_names
    automatic <- is.integer(original_names) && length(original_names) == 2L && is.na(original_names[[1L]])
    attr(result, "row.names") <- if ("tbl_df" %in% scope$restoration$table_class)
        .set_row_names(length(source_order)) else if (automatic)
            as.character(source_order) else original_names[source_order]
    class(result) <- scope$restoration$table_class
    if (!identical(identity(path), scope$identity)) stop("reader input changed")
    result
}

#' Fixed access requests; the large case is not an application store capability.
#' @param spec Prepared table descriptor.
#' @return Three explicit selections.
#' @noRd
checkpointReaderRequests <- function(spec) {
    projection <- if (spec$case == "precursor") c("Run", "Precursor.Id", "Precursor.Quantity") else
        c("Protein.Group", "logFC", "raw_pvalue", "fdr_qvalue", "comparison", "uniprot_acc")
    if (isTRUE(spec$store_adapter)) return(list(projected = list(columns = projection,
        rows = if (spec$case == "precursor") c(480001L, 489600L) else NULL)))
    list(full = list(columns = unlist(spec$columns), rows = NULL),
        projected = list(columns = projection, rows = NULL),
        page = list(columns = projection, rows = if (spec$case == "precursor")
            c(480001L, 489600L) else c(1L, 25L)))
}

#' Consume results inside the boundary without retaining a second full assay.
#' @param data Materialized reader result.
#' @param api Optional installed helpers for full-population DIA display preparation.
#' @param volcano Whether this is the complete projected DA population.
#' @return Small summary; exact files are independently checked afterward.
#' @noRd
consumeCheckpointReader <- function(data, api = NULL, volcano = FALSE) {
    numeric <- vapply(data, is.numeric, logical(1))
    result <- list(rows = nrow(data), columns = ncol(data),
        numeric_sums = lapply(data[numeric], \(column) sum(column, na.rm = TRUE)),
        text_bytes = sum(vapply(data[!numeric],
            \(column) sum(nchar(as.character(column), type = "bytes"), na.rm = TRUE), numeric(1))))
    if (volcano) {
        envelope <- api("newSavedWorkflowResultTable")(data, "Protein.Group",
            "treated_vs_control", "fdr_qvalue")
        view <- api("filterWorkflowDiaVolcano")(envelope, 0.05, 1)
        page <- api("pageWorkflowDiaResult")(envelope, 1L, 25L)
        result$volcano <- list(selected_ids = view$data$Protein.Group,
            full_population_rows = view$full_population_rows)
        result$page <- page$data$Protein.Group
    }
    result
}

#' Run one fresh-process reader member with one-off or five-action reuse.
#' @param library_path,run_dir,prepared_dir,selection Installed build and frozen member arguments.
#' @return Timings and exact-output bindings; no performance promotion.
#' @noRd
executeCheckpointReader <- function(library_path, run_dir, prepared_dir, selection) {
    api <- recordedCheckpointApi(library_path)
    parts <- strsplit(selection, ":", fixed = TRUE)[[1L]]
    manifest_path <- file.path(prepared_dir, "reader-inputs.json")
    manifest <- jsonlite::read_json(manifest_path)
    spec <- manifest$cases[[parts[[1L]]]]
    engine <- parts[[2L]]
    if (startsWith(engine, "adapter_") != isTRUE(spec$store_adapter)) {
        stop("reader engine and preparation scope differ")
    }
    repetitions <- as.integer(parts[[3L]])
    requests <- checkpointReaderRequests(spec)
    events <- list()
    summaries <- list()
    outputs <- list()
    for (index in seq_len(repetitions)) {
        action <- checkpointConsumerOperation(paste0("action_", index), function() {
            scope <- openCheckpointReader(spec, prepared_dir, engine, api)
            on.exit(closeCheckpointReader(scope), add = TRUE)
            values <- list()
            for (name in names(requests)) {
                request <- requests[[name]]
                read <- checkpointConsumerOperation(name, function() {
                    value <- readCheckpointReader(scope, request$columns, request$rows)
                    summary <- consumeCheckpointReader(value, api,
                        volcano = name == "projected" && spec$case != "precursor")
                    list(data = value, summary = summary)
                })
                values[[name]] <- read[c("elapsed_seconds", "cpu_seconds")]
                summaries[[paste(index, name)]] <<- read$value$summary
                if (isTRUE(spec$store_adapter)) closeCheckpointReader(scope)
                if (index == 1L) {
                    path <- file.path(run_dir, paste0(name, ".rds"))
                    saveRDS(read$value$data, path, compress = FALSE)
                    outputs[[name]] <<- list(file = basename(path),
                        sha256 = digest::digest(file = path, algo = "sha256"),
                        result_bytes = as.numeric(utils::object.size(read$value$data)))
                }
                rm(read)
            }
            closeCheckpointReader(scope)
            values$resources <- list(closed = scope$closed,
                open_connections = as.integer(!is.null(scope$connection)),
                retained_tables = as.integer(!is.null(scope$table)),
                store_metrics = scope$read_metrics)
            values
        })
        events[[action$label]] <- action[c("elapsed_seconds", "cpu_seconds", "value")]
        rm(action)
        gc()
    }
    receipt <- list(status = "complete", case = spec$case, engine = engine,
        repetitions = repetitions, events = events, summaries = summaries, outputs = outputs,
        prepared_dir = normalizePath(prepared_dir),
        input_sha256 = digest::digest(file = manifest_path, algo = "sha256"),
        producer = spec$producer[[if (engine %in% c("rds", "adapter_rds")) "rds" else "parquet"]],
        input_file_bytes = spec$files[[if (engine %in% c("rds", "adapter_rds")) "rds" else "parquet"]]$bytes,
        native_io_only = spec$native_io_only,
        store_adapter = isTRUE(spec$store_adapter),
        open_resources = sum(vapply(events, \(event) event$value$resources$open_connections +
            event$value$resources$retained_tables, integer(1))),
        independently_verified = FALSE, promotion_authorized = FALSE,
        threads = 2L, lazy_strings = FALSE, R = R.version.string,
        versions = lapply(c("MultiScholaR", "arrow", "duckdb", "fastduckdb"),
            \(package) list(package = package, version = as.character(utils::packageVersion(package)))))
    jsonlite::write_json(receipt, file.path(run_dir, "reader-member.json"),
        auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
    receipt
}

#' Compare materialized outputs independently with the frozen original R table.
#' @param library_path,run_dir,member_dir Explicit installed build, output and member.
#' @return Exactness receipt, never a timing result.
#' @noRd
verifyCheckpointReader <- function(library_path, run_dir, member_dir) {
    api <- recordedCheckpointApi(library_path)
    member <- jsonlite::read_json(file.path(member_dir, "reader-member.json"))
    manifest_path <- file.path(member$prepared_dir, "reader-inputs.json")
    if (!identical(member$status, "complete") ||
        !identical(digest::digest(file = manifest_path, algo = "sha256"), member$input_sha256)) {
        stop("reader verification input changed")
    }
    spec <- jsonlite::read_json(manifest_path)$cases[[member$case]]
    paths <- checkpointReaderPaths(spec, member$prepared_dir, "rds")
    original <- readRDS(paths[["payload"]])
    requests <- checkpointReaderRequests(spec)
    for (name in names(requests)) {
        request <- requests[[name]]
        selected <- if (is.null(request$rows)) seq_len(nrow(original)) else
            seq.int(request$rows[[1L]], min(request$rows[[2L]], nrow(original)))
        expected <- original[selected, request$columns, drop = FALSE]
        for (column in request$columns) {
            attributes <- attributes(original[[column]])
            attributes$names <- NULL
            attributes(expected[[column]]) <- attributes
            if (!is.null(names(original[[column]]))) names(expected[[column]]) <- names(original[[column]])[selected]
        }
        extras <- setdiff(names(attributes(original)), c("names", "row.names", "class"))
        for (extra in extras) attr(expected, extra) <- attr(original, extra)
        if (isTRUE(spec$store_adapter)) expected <- subsetCheckpointAdapterTable(
            original, request$columns, request$rows)
        path <- file.path(member_dir, paste0(name, ".rds"))
        if (!identical(digest::digest(file = path, algo = "sha256"), member$outputs[[name]]$sha256) ||
            !identical(readRDS(path), expected)) stop(paste("reader output differs:", name))
        consumed <- consumeCheckpointReader(expected, api,
            volcano = name == "projected" && spec$case != "precursor")
        # JSON changes integer storage and represents nonfinite sums as null;
        # compare using the same public receipt encoding, not relaxed table equality.
        encoded <- jsonlite::fromJSON(jsonlite::toJSON(consumed, auto_unbox = TRUE,
            null = "null", digits = NA), simplifyVector = FALSE)
        for (index in seq_len(member$repetitions)) {
            if (!identical(encoded, member$summaries[[paste(index, name)]])) {
                stop(paste("reader consumption differs:", name))
            }
        }
        rm(expected)
    }
    result <- list(status = "passed_reader_exactness", member = normalizePath(member_dir),
        input_sha256 = member$input_sha256, case = member$case, engine = member$engine,
        full_attributes = TRUE, promotion_authorized = FALSE)
    jsonlite::write_json(result, file.path(run_dir, "reader-verification.json"),
        auto_unbox = TRUE, pretty = TRUE)
    result
}

#' Prove the installed readers on a tiny exactness/lifecycle fixture before large work.
#' @param library_path,run_dir Installed proof library and owned output.
#' @return Small correctness receipt.
#' @noRd
checkCheckpointReader <- function(library_path, run_dir) {
    api <- recordedCheckpointApi(library_path)
    table <- checkpointProjectionTable(300L, 4L)
    table$edge <- rep(c(NA_real_, NaN, Inf, -Inf, 0), 60L)
    table$direction <- ordered(rep(c("up", "down"), 150L), levels = c("down", "up", "unused"))
    attr(table$logFC, "units") <- "log2"
    spec <- prepareCheckpointReaderTable(table, "test", run_dir, api)
    verified <- list()
    for (engine in c("rds", "arrow", "duckdb", "fastduckdb")) {
        scope <- openCheckpointReader(spec, run_dir, engine, api)
        tryCatch({
            for (rows in list(NULL, c(1L, 25L), c(299L, 302L), c(301L, 305L))) {
                columns <- c("Protein.Group", "edge", "direction", "comparison", "logFC")
                indices <- if (is.null(rows)) seq_len(nrow(table)) else
                    which(seq_len(nrow(table)) >= rows[[1L]] & seq_len(nrow(table)) <= rows[[2L]])
                expected <- table[indices, columns, drop = FALSE]
                attr(expected$logFC, "units") <- "log2"
                attr(expected, "provenance") <- attr(table, "provenance")
                if (!identical(readCheckpointReader(scope, columns, rows), expected)) {
                    stop(paste("tiny reader exactness failed:", engine))
                }
            }
            invalid <- tryCatch(readCheckpointReader(scope, "missing"), error = identity)
            if (!inherits(invalid, "error")) stop("invalid columns were accepted")
        }, finally = closeCheckpointReader(scope))
        if (!scope$closed || !is.null(scope$connection) || !is.null(scope$table)) stop("reader leaked")
        verified[[engine]] <- list(exact = TRUE, closed = scope$closed)
    }
    # The recorded precursor table is a tibble with element names on Run.
    # Those names are row data, not bounded metadata, and must survive projection.
    tibble_root <- file.path(run_dir, "tibble")
    dir.create(tibble_root)
    table <- tibble::as_tibble(table)
    names(table$Protein.Group) <- paste0("element_", seq_len(nrow(table)))
    spec <- prepareCheckpointReaderTable(table, "test", tibble_root, api)
    for (engine in c("rds", "arrow", "duckdb", "fastduckdb")) {
        scope <- openCheckpointReader(spec, tibble_root, engine, api)
        tryCatch({
            for (range in list(NULL, c(7L, 12L))) {
                indices <- if (is.null(range)) seq_len(nrow(table)) else seq.int(range[[1L]], range[[2L]])
                expected <- table[indices, "Protein.Group", drop = FALSE]
                names(expected$Protein.Group) <- names(table$Protein.Group)[indices]
                if (!identical(readCheckpointReader(scope, "Protein.Group", range), expected)) {
                    stop(paste("tibble reader exactness failed:", engine))
                }
            }
        }, finally = closeCheckpointReader(scope))
    }
    receipt <- list(status = "passed_reader_checks", engines = verified,
        rows = nrow(table), promotion_authorized = FALSE)
    jsonlite::write_json(receipt, file.path(run_dir, "reader-check.json"),
        auto_unbox = TRUE, pretty = TRUE)
    receipt
}

#' Dispatch one reader experiment under the existing guarded service.
#' @param operation,library_path,run_dir,input_path,selection Explicit runner arguments.
#' @return Operation result; failure never publishes successful completion.
#' @noRd
runCheckpointReaderOperation <- function(operation, library_path, run_dir,
    input_path, selection = NULL) {
    validateCheckpointReaderInput(operation, input_path, selection)
    .libPaths(unique(c(normalizePath(library_path), .libPaths())))
    receipt <- list(status = "running", operation = operation, selection = selection,
        library = normalizePath(library_path), input = normalizePath(input_path),
        promotion_authorized = FALSE)
    package_files <- file.path(library_path, "MultiScholaR",
        c("DESCRIPTION", "NAMESPACE", "R/MultiScholaR.rdb", "R/MultiScholaR.rdx"))
    receipt$installed_build <- as.list(stats::setNames(vapply(package_files,
        \(path) digest::digest(file = path, algo = "sha256"), character(1)), basename(package_files)))
    driver <- file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))),
        "workflow_checkpoint_dia_driver.R")
    if (length(driver) == 1L && file.exists(driver)) receipt$driver_sha256 <-
        digest::digest(file = driver, algo = "sha256")
    native <- file.path(library_path, "fastduckdb", "libs", "fastduckdb.so")
    if (file.exists(native)) receipt$reader_binary_sha256 <- digest::digest(file = native, algo = "sha256")
    on.exit(jsonlite::write_json(receipt, file.path(run_dir, "reader-operation.json"),
        auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA), add = TRUE)
    tryCatch({
        if (operation == "reader-check" && isTRUE(selection %in% c("perf101:D0", "perf101:D1", "perf101:D0R", "perf101:V"))) {
            snapshot <- file.path(run_dir, "source")
            receipt$execution_sources <- bindPerf101Sources(snapshot)
            saveRDS(receipt$execution_sources, file.path(run_dir, "source-manifest.rds"))
            diagnostic <- new.env(parent = asNamespace("MultiScholaR"))
            sys.source(file.path(snapshot, "tools/profiling/workflow_perf101.R"), diagnostic)
            value <- diagnostic$runPerf101Diagnostic(library_path, run_dir, snapshot, selection)
            receipt$status <- "complete"
            return(value)
        }
        value <- switch(operation,
            "reader-build" = buildCheckpointReader(library_path, run_dir, input_path),
            "reader-check" = {
                verifier_check <- identical(selection, "gui-capacity-verifier")
                diagnostic_check <- identical(selection, "gui-diagnostic")
                passive_check <- identical(selection, "gui-passive")
                lifetime_check <- identical(selection, "gui-lifetime-repair")
                sharing_check <- identical(selection, "gui-sharing-repair")
                ruv_check <- identical(selection, "gui-ruv-pearson")
                native_check <- identical(selection, "gui-native-session")
                native_outputs <- identical(selection, "gui-native-outputs")
                reference_first_qc <- isTRUE(selection %in% c("gui-reference-first-qc", "gui-reference-core", "gui-reference-routing"))
                repair_check <- identical(selection, "gui-memory-repair")
                memory_verification <- identical(selection, "gui-memory-verification")
                audit_repair <- isTRUE(selection %in% c("gui-audit-repair", "gui-audit-diagnostic"))
                audit_member <- identical(selection, "gui-audit-member")
                audit_toggle <- identical(selection, "gui-audit-toggle")
                sleep_guard <- identical(selection, "gui-sleep-guard")
                integration_check <- identical(selection, "gui-integration-repair")
                mechanism_member <- isTRUE(selection %in% c("gui-worker-member", "gui-export-member",
                    "gui-export-history-member", "gui-iq-history-member", "gui-mechanism-preflight", "gui-rollup-profile", "gui-rollup-profile-preflight", "gui-raw-replay:small", "gui-raw-replay:outputs", "gui-raw-replay:verify", "gui-cohort-profile", "gui-impute-workers:preflight", "gui-impute-workers:4", "gui-impute-workers:8", "gui-cohort-size:256", "gui-cohort-size:768", "gui-cohort-size:768-threads8", "gui-cohort-threads:1", "gui-cohort-threads:2", "gui-cohort-threads:4", "gui-cohort-threads:8", "gui-raw-replay:300")) ||
                    (!is.null(selection) && grepl("^gui-(rollup|session)-performance:(control|candidate)$", selection)) ||
                    (!is.null(selection) && grepl("^gui-bulk-rollup:(qualification-preflight|qualification-prepare|writer-build|writer-opt-(build|preflight|diagnostic:(primary|contrast))|writer-r-preflight|writer-r-readback:[A-Za-z0-9]{6}|writer-r:(screen|full):(primary|contrast):P1[456]:[12]|io-preflight|io-prepare|io:(duckdb|fastduckdb):[12]|scan-inspection|crossover:(B|P0):(4608|9600|96000|960000):[12]|contrast:(B|P([0-9]|1[0-4]))|diagnostic|diagnostic-preflight|preflight|relational-preflight|materialize-preflight|adaptation-preflight|fast-p3-preflight|writer-preflight|screen:([ABCD]|P([0-9]|1[0-4])):[12]|full:([ABCD]|P([0-9]|1[0-4])))$", selection)) ||
                    (!is.null(selection) && grepl("^gui-rollup-compare:(control|candidate)$", selection)) ||
                    (!is.null(selection) && grepl("^gui-first-qc-verify:[A-Za-z0-9]{6}:[A-Za-z0-9]{6}$", selection)) ||
                    (!is.null(selection) && grepl("^gui-history-verify:[A-Za-z0-9]{6}(:[A-Za-z0-9]{6})?$", selection))
                capacity_check <- identical(selection, "gui-capacity")
                sequence_check <- identical(selection, "sequence") || capacity_check || verifier_check || diagnostic_check || passive_check || lifetime_check || sharing_check || ruv_check || native_check || native_outputs || reference_first_qc || repair_check || memory_verification || audit_repair || audit_member || audit_toggle || sleep_guard || integration_check || mechanism_member
                checked <- if (sequence_check) {
                    recordedCheckpointApi(library_path)
                    list(status = "sequence_checks", promotion_authorized = FALSE)
                } else checkCheckpointReader(library_path, run_dir)
                test_dir <- file.path(input_path, "tests", "testthat")
                if (!dir.exists(test_dir)) stop("reader checks require the repository tests")
                origin <- jsonlite::read_json(file.path(dirname(library_path), "library-origin.json"))
                withr::with_dir(input_path, {
                    here::i_am("DESCRIPTION")
                    withr::with_envvar(c(MULTISCHOLAR_DIA008_LIBRARY = library_path,
                        MULTISCHOLAR_DIA025_AUDIT_MEASURE = if (audit_repair) "true" else "false",
                        MULTISCHOLAR_DIA025_AUDIT_DIAGNOSTIC = if (identical(selection, "gui-audit-diagnostic")) "true" else "false",
                        MULTISCHOLAR_DIA025_AUDIT_MEMBER = if (audit_member) "true" else "false",
                        MULTISCHOLAR_DIA025_WORKER_PROBE = if (diagnostic_check) "true" else "false",
                        MULTISCHOLAR_DIA025_MECHANISM_MEMBER = if (mechanism_member) selection else "",
                        MULTISCHOLAR_DIA008_B0_LIBRARY = origin$baseline_library,
                        MULTISCHOLAR_DIA008_RECEIPT = file.path(run_dir, "installed-compatibility.json")), {
                        results <- testthat::test_dir(test_dir,
                            filter = if (identical(selection, "gui-reference-core")) "workflow-saved-projection|workflow-saved-qc-rollup|workflow-checkpoint-dia-state|prot-02b-qc-peptide-module-contracts" else if (identical(selection, "gui-reference-routing")) "workflow-ruv-pearson-diagnostic" else if (reference_first_qc) "workflow-saved-projection|workflow-checkpoint-dia-state|workflow-checkpoint-harness|workflow-ruv-pearson-diagnostic|prot-02b-qc-peptide-module-contracts|prot-04b-design-server-shared" else if (native_outputs) "workflow-native-session|prot-05[bd]-norm-module|prot-12.*summary|module-ci-prot-summary|general-filemgmt-report-helpers-shared|workflow-gui-capacity-verifier|workflow-checkpoint-harness" else if (native_check) "^workflow-native-session$|^workflow-saved-qc$|^workflow-saved-da$" else if (ruv_check) "workflow-ruv-pearson-diagnostic|workflow-gui-passive-diagnostic|workflow-sleep-guard|workflow-checkpoints$|workflow-saved-recovery|workflow-saved-exports" else if (sharing_check)
                                "workflow-gui-memory-lifetimes|workflow-saved-exports|workflow-checkpoint-dia-state|prot-02c-qc-protein-module-contracts|prot-02h-peptide-audit-shared|prot-05b-norm-module-contracts|workflow-gui-capacity-verifier|workflow-memory-mechanism-member"
                            else if (passive_check) "workflow-gui-passive-diagnostic" else if (lifetime_check)
                                "workflow-gui-memory-lifetimes|workflow-checkpoint-dia-state|workflow-saved-exports|prot-01c-import-module-contracts|prot-fasta-processing-direct-shared|prot-04-design$|prot-04a-design-module-characterization|prot-04b-design-server-shared|prot-02b-qc-peptide-module-contracts|prot-02c-qc-protein-module-contracts|prot-02h-peptide-audit-shared|prot-05b-norm-module-contracts|prot-05d-norm-module-characterization|prot-07b-da-handlers-compat|module-ci-prot-da|protein-da-golden-master|workflow-gui-capacity-verifier|workflow-gui-passive-diagnostic"
                            else if (identical(selection, "gui-bulk-rollup:writer-r-preflight")) "workflow-bulk-relational" else if (!is.null(selection) && grepl("^gui-bulk-rollup:(qualification-preflight|qualification-prepare|writer-build|writer-opt-(build|preflight|diagnostic:(primary|contrast))|writer-r-preflight|writer-r-readback:[A-Za-z0-9]{6}|writer-r:(screen|full):(primary|contrast):P1[456]:[12]|io-preflight|io-prepare|io:(duckdb|fastduckdb):[12]|scan-inspection|crossover:|contrast:)", selection)) "workflow-rollup-qualification" else if (!is.null(selection) && selection %in% c("gui-bulk-rollup:relational-preflight", "gui-bulk-rollup:materialize-preflight", "gui-bulk-rollup:adaptation-preflight", "gui-bulk-rollup:fast-p3-preflight", "gui-bulk-rollup:writer-preflight")) "workflow-bulk-relational" else if (mechanism_member) "workflow-memory-mechanism-member" else if (integration_check)
                                "workflow-saved-resources|workflow-gui-memory-lifetimes|workflow-saved-da$|workflow-saved-qc$|workflow-saved-exports|workflow-checkpoint-dia-state|workflow-gui-capacity-verifier|general-filemgmt-config-contracts|prot-01c-import-module-contracts|prot-04b-design-server-shared|prot-02c-qc-protein-module-contracts|prot-02h-qc-protein-replicate-shared|prot-05b-norm-module-contracts|prot-07b-da-handlers-compat|module-ci-prot-da|protein-da-golden-master|prot-summary-publication"
                            else if (sleep_guard) "workflow-sleep-guard" else if (audit_toggle) "workflow-server-shared|prot-02h-peptide-audit-shared|workflow-checkpoint-dia-state" else if (audit_member) "workflow-audit-memory-member" else if (audit_repair) "prot-02h-peptide-audit-shared" else if (memory_verification) "workflow-gui-memory-(browser|repairs)|protein-da-golden-master" else if (repair_check) "workflow-gui-memory-(repairs|browser)|workflow-checkpoint-dia-state|prot-02h-peptide-audit-shared|prot-02c-qc-protein-module-contracts|prot-05b-norm-module-contracts|prot-05d-norm-module-characterization|prot-08-volcano|module-ci-prot-protein-qc|module-ci-prot-da|protein-da-golden-master" else if (diagnostic_check) "workflow-gui-diagnostic" else if (verifier_check) "workflow-gui-capacity-verifier|workflow-memory-mechanism-member" else if (capacity_check)
                                "workflow-checkpoint-dia-state|workflow-checkpoint-harness|workflow-saved-installed|prot-01c-import-module-contracts|prot-04b-design-server-shared|prot-02b-qc-peptide-module-contracts|prot-02c-qc-protein-module-contracts|prot-05b-norm-module-contracts|prot-07b-da-handlers-compat" else if (sequence_check)
                                "workflow-checkpoint-harness|prot-02b-qc-peptide-module-contracts|workflow-saved-qc-rollup|workflow-saved-installed" else
                                "workflow-saved-projection|workflow-saved-resources|workflow-checkpoint-harness|workflow-saved-installed",
                            reporter = testthat::SummaryReporter$new(max_reports = Inf), package = "MultiScholaR",
                            load_package = "installed", stop_on_failure = FALSE)
                        results <- as.data.frame(results)
                        saveRDS(results, file.path(run_dir, "reader-tests.rds"))
                        if (any(results$failed > 0L | results$error)) stop("installed checks failed")
                    })
                })
                checked
            },
            "reader-prepare" = prepareCheckpointReader(library_path, run_dir, input_path,
                adapter = identical(selection, "adapter")),
            "reader-inspect" = {
                recordedCheckpointApi(library_path)
                manifest <- readRecordedCheckpointManifest(input_path)
                record <- readRecordedCheckpoint(input_path, manifest, "raw_data_s4")
                table <- methods::slot(record$data, "peptide_data")
                result <- list(class = class(table), shape = dim(table),
                    attributes = lapply(attributes(table), function(value) list(
                        type = typeof(value), class = class(value), length = length(value))),
                    structure = utils::capture.output(utils::str(table[0, ], max.level = 3)))
                jsonlite::write_json(result, file.path(run_dir, "reader-schema.json"),
                    auto_unbox = TRUE, pretty = TRUE)
                result
            },
            "reader-member" = executeCheckpointReader(library_path, run_dir, input_path, selection),
            "reader-verify" = verifyCheckpointReader(library_path, run_dir, input_path))
        receipt$status <- "complete"
        value
    }, error = function(error) {
        receipt$status <<- "failed"
        receipt$error <<- conditionMessage(error)
        stop(error)
    }, interrupt = function(error) {
        receipt$status <<- "interrupted"
        stop(error)
    })
}

#' Execute full-RDS versus Parquet/DuckDB access comparisons.
#' @param library_path Installed candidate library.
#' @param run_dir Fresh request-owned output directory.
#' @param rows,samples Frozen access-table dimensions; defaults preserve DIA-013.
#' @return Descriptive construction, persistence and repeated-access receipt.
#' @noRd
runCheckpointProjectionScenarios <- function(library_path, run_dir,
    rows = 8020L, samples = 24L) {
    if (!is.numeric(rows) || length(rows) != 1L || !is.finite(rows) ||
        rows < 50L || rows > 10000L || rows != floor(rows) ||
        !is.numeric(samples) || length(samples) != 1L || !is.finite(samples) ||
        samples < 2L || samples > 63L || samples != floor(samples)) {
        stop("unsupported access-table dimensions")
    }
    library_path <- normalizePath(library_path, mustWork = TRUE)
    run_dir <- normalizePath(run_dir, mustWork = TRUE)
    .libPaths(unique(c(library_path, .libPaths())))
    namespace <- loadNamespace("MultiScholaR", lib.loc = library_path)
    api <- \(name) get(name, envir = namespace, inherits = FALSE)
    events <- list()
    constructed <- checkpointConsumerOperation("table_construction",
        \() checkpointProjectionTable(rows, samples))
    events[[constructed$label]] <- constructed[c("elapsed_seconds", "cpu_seconds")]
    table <- constructed$value
    identity <- api("newSavedWorkflowScientificIdentity")(
        digest::digest(table, algo = "sha256"), "projection-design-v1",
        "adjusted-full-population-v1", "page-volcano-heatmap-v1",
        paste0("MultiScholaR-", as.character(utils::packageVersion("MultiScholaR"))))
    rds_path <- file.path(run_dir, "da-table.rds")
    rds_write <- checkpointConsumerOperation("rds_persist_fingerprint", \() {
        saveRDS(table, rds_path)
        digest::digest(file = rds_path, algo = "sha256")
    })
    events[[rds_write$label]] <- rds_write[c("elapsed_seconds", "cpu_seconds")]
    store <- api("newWorkflowEntityStore")(run_dir, "projection-consumers", "DIA")
    parquet_write <- checkpointConsumerOperation("parquet_persist_fingerprint", \() {
        api("saveWorkflowTable")(store, table, "da-table", identity, ".rownames",
            factors = list(route = "iq", layer = "full_results"))
    })
    events[[parquet_write$label]] <- parquet_write[c("elapsed_seconds", "cpu_seconds")]
    reference <- parquet_write$value
    page_columns <- c("Protein.Group", "logFC", "raw_pvalue", "fdr_qvalue",
        "comparison", "uniprot_acc")
    volcano_columns <- page_columns
    norm_columns <- grep("^log2norm\\.", names(table), value = TRUE)
    heatmap_columns <- c("Protein.Group", norm_columns)
    repetitions <- 5L
    rds_elapsed <- numeric()
    parquet_elapsed <- numeric()
    parquet_collected_bytes <- 0
    first <- list()
    for (index in seq_len(repetitions)) {
        rds_page <- checkpointConsumerOperation(paste0("rds_page_", index), \() {
            value <- readRDS(rds_path)
            value[1:25, page_columns, drop = FALSE]
        })
        events[[rds_page$label]] <- rds_page[c("elapsed_seconds", "cpu_seconds")]
        rds_elapsed <- c(rds_elapsed, rds_page$elapsed_seconds)
        parquet_page <- checkpointConsumerOperation(paste0("parquet_page_", index), \() {
            api("readSavedWorkflowTable")(store, reference, page_columns,
                offset = 0L, limit = 25L)
        })
        events[[parquet_page$label]] <- parquet_page[c("elapsed_seconds", "cpu_seconds")]
        parquet_elapsed <- c(parquet_elapsed, parquet_page$elapsed_seconds)
        parquet_collected_bytes <- parquet_collected_bytes +
            parquet_page$value$collected_bytes
        rds_volcano <- checkpointConsumerOperation(paste0("rds_volcano_", index), \() {
            value <- readRDS(rds_path)
            value[, volcano_columns, drop = FALSE]
        })
        events[[rds_volcano$label]] <- rds_volcano[c("elapsed_seconds", "cpu_seconds")]
        rds_elapsed <- c(rds_elapsed, rds_volcano$elapsed_seconds)
        parquet_volcano <- checkpointConsumerOperation(
            paste0("parquet_volcano_", index), \() {
                api("readSavedWorkflowTable")(store, reference, volcano_columns,
                    offset = 0L, limit = nrow(table), max_rows = 10000L)
            })
        events[[parquet_volcano$label]] <-
            parquet_volcano[c("elapsed_seconds", "cpu_seconds")]
        parquet_elapsed <- c(parquet_elapsed, parquet_volcano$elapsed_seconds)
        parquet_collected_bytes <- parquet_collected_bytes +
            parquet_volcano$value$collected_bytes
        rds_heatmap <- checkpointConsumerOperation(paste0("rds_heatmap_", index), \() {
            value <- readRDS(rds_path)
            order <- order(value$fdr_qvalue, -abs(value$logFC), seq_len(nrow(value)))
            value[order[1:50], heatmap_columns, drop = FALSE]
        })
        events[[rds_heatmap$label]] <- rds_heatmap[c("elapsed_seconds", "cpu_seconds")]
        rds_elapsed <- c(rds_elapsed, rds_heatmap$elapsed_seconds)
        parquet_heatmap <- checkpointConsumerOperation(
            paste0("parquet_heatmap_", index), \() {
                stats <- api("readSavedWorkflowTable")(store, reference,
                    c("Protein.Group", "logFC", "fdr_qvalue"), offset = 0L,
                    limit = nrow(table), max_rows = 10000L)
                order <- order(stats$data$fdr_qvalue, -abs(stats$data$logFC),
                    seq_len(nrow(stats$data)))
                keys <- stats$data$Protein.Group[order[1:50]]
                matrix <- api("readSavedWorkflowTable")(store, reference,
                    heatmap_columns, keys = keys, limit = 50L, max_rows = 50L)
                data <- matrix$data
                rownames(data) <- matrix$source_row_names
                list(data = data, keys = keys,
                    collected_bytes = stats$collected_bytes + matrix$collected_bytes)
            })
        events[[parquet_heatmap$label]] <-
            parquet_heatmap[c("elapsed_seconds", "cpu_seconds")]
        parquet_elapsed <- c(parquet_elapsed, parquet_heatmap$elapsed_seconds)
        parquet_collected_bytes <- parquet_collected_bytes +
            parquet_heatmap$value$collected_bytes
        if (index == 1L) {
            first <- list(rds_page = rds_page$value, parquet_page = parquet_page$value$data,
                rds_volcano = rds_volcano$value,
                parquet_volcano = parquet_volcano$value$data,
                rds_heatmap = rds_heatmap$value,
                parquet_heatmap = parquet_heatmap$value$data)
        }
    }
    if (!isTRUE(all.equal(unclass(first$parquet_page), unclass(first$rds_page),
            check.attributes = FALSE)) ||
        !isTRUE(all.equal(unclass(first$parquet_volcano), unclass(first$rds_volcano),
            check.attributes = FALSE)) ||
        !identical(rownames(first$parquet_heatmap), rownames(first$rds_heatmap)) ||
        !isTRUE(all.equal(unclass(first$parquet_heatmap),
            unclass(first$rds_heatmap), check.attributes = FALSE))) {
        stop("projection consumer results differ from full RDS")
    }
    full_decoded_bytes <- as.numeric(utils::object.size(table))
    rm(table, first, constructed, rds_page, parquet_page, rds_volcano,
        parquet_volcano, rds_heatmap, parquet_heatmap)
    gc()
    scoped_elapsed <- numeric(repetitions)
    compound_rds_elapsed <- numeric(repetitions)
    scoped_collected_bytes <- 0
    for (index in seq_len(repetitions)) {
        rds_action <- checkpointConsumerOperation(
            paste0("rds_compound_", index), \() {
                value <- readRDS(rds_path)
                ordered <- order(value$fdr_qvalue, -abs(value$logFC),
                    seq_len(nrow(value)), na.last = TRUE)
                list(page = value[1:25, page_columns, drop = FALSE],
                    volcano = value[, volcano_columns, drop = FALSE],
                    heatmap = value[ordered[1:50], heatmap_columns, drop = FALSE])
            })
        events[[rds_action$label]] <- rds_action[c("elapsed_seconds", "cpu_seconds")]
        compound_rds_elapsed[[index]] <- rds_action$elapsed_seconds
        scoped_action <- checkpointConsumerOperation(
            paste0("parquet_compound_", index), \() {
                scope <- api("newWorkflowTableReadScope")(store, reference)
                on.exit(api("releaseWorkflowTableReadScope")(scope), add = TRUE)
                page <- api("readSavedWorkflowTable")(store, reference, page_columns,
                    offset = 0L, limit = 25L, scope = scope)
                volcano <- api("readSavedWorkflowTable")(store, reference,
                    volcano_columns, offset = 0L, limit = rows, scope = scope)
                heatmap <- api("readSavedWorkflowDiaHeatmapTable")(store, reference,
                    heatmap_columns, 50L, scope = scope)
                list(page = page$data, volcano = volcano$data,
                    heatmap = heatmap$data, collected_bytes =
                        page$collected_bytes + volcano$collected_bytes +
                        heatmap$collected_bytes)
            })
        events[[scoped_action$label]] <-
            scoped_action[c("elapsed_seconds", "cpu_seconds")]
        scoped_elapsed[[index]] <- scoped_action$elapsed_seconds
        scoped_collected_bytes <- scoped_collected_bytes +
            scoped_action$value$collected_bytes
        if (!isTRUE(all.equal(unclass(scoped_action$value$page),
                unclass(rds_action$value$page), check.attributes = FALSE)) ||
            !isTRUE(all.equal(unclass(scoped_action$value$volcano),
                unclass(rds_action$value$volcano), check.attributes = FALSE)) ||
            !identical(rownames(scoped_action$value$heatmap),
                rownames(rds_action$value$heatmap)) ||
            !isTRUE(all.equal(unclass(scoped_action$value$heatmap),
                unclass(rds_action$value$heatmap), check.attributes = FALSE))) {
            stop("scoped projection differs from full RDS")
        }
        rm(rds_action, scoped_action)
        gc()
    }
    cache <- api("workflowEntityCacheInfo")(store)
    api("closeWorkflowEntityStore")(store)
    receipt <- list(status = "passed_projection_scenarios", rows = rows,
        columns = 9L + 2L * samples, repetitions = repetitions, events = events,
        rds_file_bytes = file.size(rds_path), parquet_file_bytes = reference$payload$bytes,
        full_decoded_bytes = full_decoded_bytes,
        repeated_full_rds_decoded_bytes = full_decoded_bytes * repetitions * 3,
        repeated_parquet_collected_bytes = parquet_collected_bytes,
        median_rds_elapsed_seconds = stats::median(rds_elapsed),
        median_parquet_elapsed_seconds = stats::median(parquet_elapsed),
        parquet_to_rds_elapsed_ratio = stats::median(parquet_elapsed) /
            stats::median(rds_elapsed),
        median_compound_rds_elapsed_seconds = stats::median(compound_rds_elapsed),
        median_compound_scoped_elapsed_seconds = stats::median(scoped_elapsed),
        scoped_to_rds_compound_ratio = stats::median(scoped_elapsed) /
            stats::median(compound_rds_elapsed),
        repeated_scoped_collected_bytes = scoped_collected_bytes,
        table_reads = cache$table_reads,
        table_rows = cache$table_rows, table_columns = cache$table_columns,
        table_bytes = cache$table_bytes, open_resources = cache$open_resources,
        exact_results = TRUE, performance_claim = FALSE,
        promotion_authorized = FALSE)
    jsonlite::write_json(receipt, file.path(run_dir, "projection-consumers.json"),
        auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
    receipt
}

#' Execute original double serialization versus export-scoped single serialization.
#' @param library_path Installed candidate library.
#' @param run_dir Fresh request-owned output directory.
#' @param contract_path Frozen recipe contract.
#' @return Descriptive complete-export and restore receipt.
#' @noRd
runCheckpointExportScenarios <- function(library_path, run_dir, contract_path) {
    library_path <- normalizePath(library_path, mustWork = TRUE)
    run_dir <- normalizePath(run_dir, mustWork = TRUE)
    contract_path <- normalizePath(contract_path, mustWork = TRUE)
    .libPaths(unique(c(library_path, .libPaths())))
    namespace <- loadNamespace("MultiScholaR", lib.loc = library_path)
    api <- \(name) get(name, envir = namespace, inherits = FALSE)
    source_root <- normalizePath(file.path(dirname(contract_path), "..", "..", ".."))
    helpers <- new.env(parent = baseenv())
    sys.source(file.path(source_root, "tests", "testthat",
        "helper-workflow-checkpoint-fixtures.R"), envir = helpers)
    contract <- jsonlite::read_json(contract_path, simplifyVector = FALSE)
    recipe <- Filter(\(item) identical(item$id, "small"), contract$recipes)[[1L]]
    store_root <- file.path(run_dir, "store")
    original_root <- file.path(run_dir, "original")
    optimized_root <- file.path(run_dir, "optimized")
    dir.create(store_root)
    dir.create(original_root)
    dir.create(optimized_root)
    store <- api("newWorkflowEntityStore")(store_root, "export-consumers", "DIA")
    references <- list()
    for (generation in seq_len(recipe$history_depth)) {
        object <- helpers$makeWorkflowCheckpointProtein(recipe, generation)
        references[[generation]] <- api("saveWorkflowCheckpoint")(store, object,
            sprintf("state-%02d", generation),
            identity = list(input = digest::digest(object, algo = "sha256"),
                design = "small-design", parameters = paste0("generation-", generation),
                selectors = sprintf("state-%02d", generation),
                software = "export-consumer-v1"),
            metadata = list(class = class(object),
                shape = as.integer(dim(methods::slot(object, "protein_quant_table"))),
                factors = list(route = "iq")))
    }
    state_names <- sprintf("state_%02d", seq_len(recipe$history_depth))
    states <- lapply(seq_along(state_names), function(index) {
        list(timestamp = as.POSIXct(index, origin = "1970-01-01"),
            data = references[[index]], config = list(route = "iq"),
            description = state_names[[index]],
            previous_state = if (index == 1L) NA_character_ else state_names[[index - 1L]],
            s4_class = "ProteinQuantitativeData",
            audit_metadata = list(provenance_status = "fixture"))
    })
    names(states) <- state_names
    states$state_alias <- states[[1L]]
    state_names <- names(states)
    template <- list(r6_current_state_name = state_names[[recipe$history_depth]],
        r6_complete_states = states, r6_state_history = as.list(state_names),
        current_s4_object = references[[recipe$history_depth]],
        correlation_filtered_s4 = references[[recipe$history_depth]],
        config_list = list(globalParameters = list(workflow_type = "DIA")),
        export_timestamp = as.POSIXct("2026-09-16 16:30:00", tz = "UTC"),
        final_protein_count = recipe$n_proteins,
        final_sample_count = recipe$n_samples)
    scope <- api("newWorkflowExportHydrationScope")(store)
    session <- api("hydrateWorkflowExportSession")(scope, template)
    hydration <- api("workflowExportHydrationInfo")(scope)
    api("releaseWorkflowExportHydrationScope")(scope)
    rm(object, template, references)
    gc(full = TRUE)
    heap_before_serialization <- checkpointUtilityHeapBytes()
    original <- checkpointConsumerOperation("original_double_serialization", \() {
        first <- file.path(original_root, "filtered_session_data_20260916_163000.rds")
        latest <- file.path(original_root, "filtered_session_data_latest.rds")
        saveRDS(session, first)
        saveRDS(session, latest)
        list(first = first, latest = latest)
    })
    tsv <- file.path(run_dir, "normalised_results.tsv")
    figure <- file.path(run_dir, "qc.png")
    writeLines(c("Protein.Group\tS1", "PG_000001\t1"), tsv)
    writeBin(as.raw(c(137L, 80L, 78L, 71L, 1L, 2L, 3L)), figure)
    serializations <- 0L
    optimized <- checkpointConsumerOperation("single_serialization_verified_copy", \() {
        api("runWorkflowLegacyExport")(session, optimized_root, "export-consumer",
            "20260916_163000",
            additional_files = list(table = list(source = tsv,
                destination = "Publication_tables/normalised_results.tsv"),
                figure = list(source = figure,
                    destination = "Publication_figures/QC/qc.png")),
            serializeFn = function(value, path) {
                serializations <<- serializations + 1L
                saveRDS(value, path)
            })
    })
    restored_timestamped <- api("validateExternalWorkflowLegacySession")(
        optimized$value$timestamped_path, state_names)
    restored_latest <- api("validateExternalWorkflowLegacySession")(
        optimized$value$latest_path, state_names)
    original_first <- readRDS(original$value$first)
    original_latest <- readRDS(original$value$latest)
    if (!identical(original_first, session) || !identical(original_latest, session) ||
        !identical(restored_timestamped, session) || !identical(restored_latest, session)) {
        stop("legacy export comparison changed session contents")
    }
    cache <- api("workflowEntityCacheInfo")(store)
    api("closeWorkflowEntityStore")(store)
    receipt <- list(status = "passed_export_scenarios", route = "iq",
        recipe = "small", states = length(state_names), hydration = hydration,
        heap_before_serialization = heap_before_serialization,
        original_elapsed_seconds = original$elapsed_seconds,
        optimized_elapsed_seconds = optimized$elapsed_seconds,
        elapsed_ratio = optimized$elapsed_seconds / original$elapsed_seconds,
        serialization_count = serializations,
        scientific_checkpoint_count = optimized$value$scientific_checkpoint_count,
        timestamped_sha256 = digest::digest(file = optimized$value$timestamped_path,
            algo = "sha256"),
        latest_sha256 = digest::digest(file = optimized$value$latest_path,
            algo = "sha256"),
        timestamped_bytes = file.size(optimized$value$timestamped_path),
        latest_bytes = file.size(optimized$value$latest_path),
        additional_files = names(optimized$value$additional),
        store_payload_reads = cache$payload_reads,
        store_payload_decodes = cache$payload_decodes,
        open_resources = cache$open_resources, exact_restore = TRUE,
        performance_claim = FALSE, promotion_authorized = FALSE)
    jsonlite::write_json(receipt, file.path(run_dir, "export-consumers.json"),
        auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
    receipt
}

#' Verify first-QC automatic figures and exercise the actual GUI Revert button.
#' @param driver Installed browser.
#' @param run_dir Existing guarded member directory.
#' @return Bounded output and Revert receipt, without source hydration.
#' @noRd
checkCheckpointFirstQc <- function(driver, run_dir) {
    plots <- file.path(run_dir, "dia_baseline_smoke", "results", "proteomics", "publication_graphs", "filtering_qc")
    files <- list.files(plots, pattern = "^2_qval_filtered_.*[.](png|pdf)$", recursive = TRUE, full.names = TRUE)
    files <- files[!grepl("/.workflow-export/", files, fixed = TRUE)]
    expected <- as.vector(outer(paste0("2_qval_filtered_", c("proteins_total", "proteins_per_run",
        "peptides_total", "peptides_per_protein", "peptides_per_run", "combined_plots")),
        c("png", "pdf"), paste, sep = "."))
    if (!identical(sort(basename(files)), sort(expected)) || any(file.size(files) <= 0)) {
        stop("first-QC required figure inventory is incomplete")
    }
    for (file in files[endsWith(files, ".png")]) {
        picture <- png::readPNG(file)
        if (length(dim(picture)) < 2L || any(dim(picture)[1:2] < 100L)) stop("invalid first-QC PNG")
        picture <- NULL
    }
    id <- "proteomics_workflow-quality_control-peptide_qc-qvalue_filter-revert_qvalue"
    driver$get_js(sprintf("document.getElementById(%s).click();true", jsonlite::toJSON(id, auto_unbox = TRUE)))
    output_id <- "proteomics_workflow-quality_control-peptide_qc-qvalue_filter-qvalue_results"
    checkpointWaitForJs(driver, sprintf("document.getElementById(%s).textContent.includes('Reverted to raw data state')",
        jsonlite::toJSON(output_id, auto_unbox = TRUE)), timeout = checkpointWaitBudget())
    driver$wait_for_idle()
    state <- jsonlite::read_json(file.path(run_dir, "gui-state-latest.json"))
    if (!identical(state$current_state, "raw_data_s4") || state$full_s4_reads != 0L ||
        state$store$open_resources != 0L) stop("first-QC Revert is not reference-only")
    list(status = "passed_first_qc_outputs_revert", files = basename(files),
        bytes = file.size(files), current_state = state$current_state, full_s4_reads = state$full_s4_reads)
}

#' Request bounded import/design attribution at an idle diagnostic boundary.
#' @param driver Installed diagnostic browser.
#' @param run_dir Existing evidence directory.
#' @param label Import or completed Design setup, before scientific Save.
#' @return NULL after the application records completion; never inferred from startup.
#' @noRd
requestCheckpointDesignMemory <- function(driver, run_dir, label) {
    if (!label %in% c("import", "design_setup")) stop("invalid design memory boundary")
    driver$get_js(paste0("Shiny.setInputValue('.dia_memory_probe', ",
        jsonlite::toJSON(label, auto_unbox = TRUE), ", {priority:'event'});true"))
    path <- file.path(run_dir, paste0("memory-boundary-", label, ".json"))
    deadline <- unname(proc.time()[["elapsed"]]) + 60
    repeat {
        if (file.exists(path)) {
            receipt <- jsonlite::read_json(path)
            if (!identical(receipt$status, "complete")) stop("memory diagnostic did not complete")
            return(invisible(NULL))
        }
        if (unname(proc.time()[["elapsed"]]) >= deadline) stop("memory diagnostic acknowledgement timed out")
        Sys.sleep(0.1)
    }
}

#' Bind the frozen source set used by a PERF-101 guarded process.
#' @param snapshot Request-owned source directory created before service launch.
#' @return Named SHA256 identities of the executed source and protocol bytes.
#' @noRd
bindPerf101Sources <- function(snapshot) {
    required <- c("tools/profiling/run_workflow_checkpoint_proof.R",
        "tools/profiling/run_workflow_checkpoint_proof.sh",
        "tools/profiling/workflow_checkpoint_dia_driver.R", "tools/profiling/workflow_perf101.R",
        "tests/testthat/helper-workflow-bulk-rollup.R",
        "tests/testdata/workflow-execution/qualification.json",
        "tickets/openai/verification/pilot-source-bindings.json",
        "tickets/openai/verification/pilot-protocol.json", "R/func_prot_qc_peptide_methods.R")
    paths <- file.path(snapshot, required)
    if (!all(file.exists(paths))) stop("PERF-101 source snapshot is incomplete")
    stats::setNames(vapply(paths, function(path) digest::digest(file = path,
        algo = "sha256"), character(1)), required)
}
