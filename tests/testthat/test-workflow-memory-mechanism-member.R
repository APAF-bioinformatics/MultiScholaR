#' Hash fixed complete-group blocks independently of the cohort's storage boundaries.
#' @param value Restored complete-group cohort.
#' @param automatic_names Whether row names are implicit rather than identifiers.
#' @return Ordered fingerprints; no scientific values or attributes are discarded.
#' @noRd
fingerprintDIACohortBlocks <- function(value, automatic_names) {
    positions <- attr(value, "workflow_source_order")
    stopifnot(nrow(value) %% 9600L == 0L)
    vapply(seq.int(1L, nrow(value), by = 9600L), function(start) {
        at <- seq.int(start, length.out = 9600L)
        part <- value[at, , drop = FALSE]
        attr(part, "workflow_source_order") <- positions[at]
        if (automatic_names) attr(part, "row.names") <- .set_row_names(length(at))
        digest::digest(part, algo = "sha256")
    }, character(1))
}

#' Execute one installed roll-up or session-export operation on an owned replay.
#' @param manager Freshly copied adapter, never the original evidence owner.
#' @param mechanism Either rollup or session.
#' @param native_session Use the approved native export for the candidate only.
#' @param source_dir Original frozen member, read only.
#' @param output Fresh request-owned output directory.
#' @return Small output receipt; no assay escapes the operation.
#' @noRd
runDIAIsolatedMechanism <- function(manager, mechanism, native_session, source_dir, output) {
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    if (mechanism == "rollup") {
        manager$current_state <- "protein_peptide_filtered"
        history <- manager$state_history
        manager$state_history <- history[seq_len(match(manager$current_state, unlist(history)))]
        manager$.context$sql_rollup <- TRUE
        withDIAStateScope(manager, function() prepareDIACapacityRollupInput(manager, manager$getState()))
        withDIAStateScope(manager, function() {
            runPeptideRollupApplyStep(workflow)
            invisible(NULL)
        })
        return(list(state = "precursor_rollup", boundary = "projection/envelope preparation, scientific Apply and durable checkpoint; GUI rendering excluded"))
    }
    stopifnot(mechanism == "session")
    descriptors <- list.files(manager$.context$store$directory, "^descriptor[.]rds$",
        recursive = TRUE, full.names = TRUE)
    stopifnot(all(file.size(descriptors) <= 65536))
    descriptors <- lapply(descriptors, readRDS)
    selected <- Filter(function(value) identical(value$identity$parameters, "export-metadata"), descriptors)
    stopifnot(length(selected) == 1L)
    metadata <- resolveDIAStateValue(manager$.context, selected[[1L]])
    for (field in c("contrasts_tbl", "design_matrix", "config_list", "fasta_metadata",
        "accession_cleanup_results", "ruv_optimization_result", "qc_params", "protein_counts",
        "mixed_species_analysis")) workflow[[field]] <- metadata[[field]]
    project <- file.path(source_dir, "dia_baseline_smoke")
    for (field in c("aa_seq_tbl_final", "uniprot_dat_cln")) {
        path <- file.path(project, "scripts", "proteomics", paste0(field, ".RDS"))
        if (file.exists(path)) workflow[[field]] <- readRDS(path)
    }
    if (native_session) {
        manager$.context$session_format <- "native_v1"
        manager$.context$experiment_paths <- list(source_dir = output,
            protein_qc_dir = file.path(project, "results", "proteomics", "protein_qc"),
            peptide_qc_dir = file.path(project, "results", "proteomics", "peptide_qc"))
        config <- getDIAStateRecord(manager, "protein_s4_created")$config
        snapshot <- manager$.context$snapshot
        for (field in c("iq_output_file", "peptide_input_file")) {
            path <- config[[field]]
            stopifnot(file.exists(path))
            snapshot$file_dependencies[[field]] <- c(list(source = path,
                destination = paste("execution", field, basename(path), sep = "/")), describeDIANativeFile(path))
        }
        publishDIAStateHead(manager$.context, snapshot)
    }
    result <- withDIAStateScope(manager, function() {
        norm <- list(correlation_filtered_obj = manager$getState(), best_k = metadata$ruv_k,
            correlation_threshold = metadata$correlation_threshold)
        runProtNormExportSessionWorkflow(workflow, norm,
            list(norm_method = metadata$normalization_method, ruv_mode = metadata$ruv_mode), output,
            withProgressFn = function(message, value, expr) force(expr), incProgressFn = function(...) NULL)
    })
    stopifnot(is.null(result$sessionData), identical(result$exportArtifacts$publication$status, "complete"))
    list(outputs = result$exportArtifacts, native_session = native_session,
        boundary = "metadata/dependency preparation, complete session export, packaging/copies/fingerprints/publication")
}

#' Reopen an owned copy, leaving all original scientific evidence untouched.
#' @param source_store Existing complete saved store directory.
#' @param destination Fresh empty fixture root under the guarded member.
#' @return Reopened manager over the copied store.
#' @noRd
copyDIAHistoryFixture <- function(source_store, destination) {
    source_store <- normalizePath(source_store, mustWork = TRUE)
    destination <- normalizePath(destination, mustWork = TRUE)
    if (length(list.files(destination, all.files = TRUE, no.. = TRUE)) ||
        identical(source_store, destination) || startsWith(destination, paste0(source_store, "/"))) {
        stop("history fixture requires a fresh independent destination")
    }
    identity <- readSavedStoreIdentity(file.path(source_store, "store.json"))
    if (!file.copy(source_store, destination, recursive = TRUE, copy.date = TRUE)) {
        stop("history fixture copy failed")
    }
    newDIAStateAdapter(newWorkflowEntityStore(destination, identity$project_id,
        identity$workflow_id), reopen = TRUE)
}

#' Compare one state without retaining a second complete history.
#' @param expected,actual Historical state records.
#' @param expected_roots,actual_roots Verified volatile member-path roots.
#' @param comparable Existing harness scientific comparison policy.
#' @return NULL invisibly; any scientific or metadata difference raises an error.
#' @noRd
checkDIAHistoryFixtureRecord <- function(expected, actual, expected_roots, actual_roots, comparable) {
    expected <- comparable(expected, expected_roots)
    actual <- comparable(actual, actual_roots)
    comparison <- all.equal(actual, expected, tolerance = 1e-12, check.attributes = TRUE)
    if (!isTRUE(comparison)) stop(paste("historical state differs:", paste(comparison, collapse = "; ")))
    if (isS4(actual$data)) {
        for (field in intersect(c("peptide_data", "protein_quant_table", "design_matrix", "protein_id_table"),
            methods::slotNames(actual$data))) {
            left <- methods::slot(actual$data, field)
            right <- methods::slot(expected$data, field)
            if (!identical(attributes(left), attributes(right))) stop(paste(
                "historical table attributes differ:", field,
                "; actual attribute order:", paste(names(attributes(left)), collapse = ","),
                "; expected attribute order:", paste(names(attributes(right)), collapse = ",")))
            for (index in seq_along(left)) {
                if (!identical(typeof(left[[index]]), typeof(right[[index]])) ||
                    !identical(attributes(left[[index]]), attributes(right[[index]])) ||
                    !identical(is.na(left[[index]]), is.na(right[[index]])) ||
                    !identical(is.nan(left[[index]]), is.nan(right[[index]]))) {
                    stop("historical column type, attributes or missingness differ")
                }
            }
        }
    }
    invisible(NULL)
}

#' Map only the two already fingerprinted IQ files from the replay layout.
#' @param record Protein rollup record whose other fields remain untouched.
#' @param files Exact owned paths checked against the original files by the caller.
#' @return Record with checked IQ paths represented by the same artifact names.
#' @noRd
normalizeDIAIqReplayFiles <- function(record, files) {
    stopifnot(identical(names(files), c("peptide_input_file", "iq_output_file")))
    for (field in names(files)) {
        path <- record$config[[field]]
        if (!is.character(path) || length(path) != 1L ||
            !identical(normalizePath(path, mustWork = TRUE), normalizePath(files[[field]], mustWork = TRUE))) {
            stop("IQ replay file location differs from the verified artifact")
        }
        record$config[[field]] <- paste0("<IQ_FILE>/", basename(files[[field]]))
    }
    record
}

#' Bind this continuation to its previously verified complete large DIA history.
#' @param proof_root Existing guarded proof root, derived from the member receipt.
#' @return Frozen member metadata and its store directory.
#' @noRd
readDIAHistoryFixture <- function(proof_root) {
    member_dir <- file.path(proof_root, "browser-smoke.IVL1gM")
    member <- jsonlite::read_json(file.path(member_dir, "smoke-receipt.json"))
    stopifnot(identical(member$status, "passed_summary_ui_only"),
        identical(member$capacity$treatment, "qs2_duckdb"),
        identical(member$input_sha256$report,
            "2414125be92a11fbe4342fdf354e0214915e30592bec66f3e9d6a984f9290841"),
        identical(digest::digest(file = file.path(member$installation$package_path,
            "R", "MultiScholaR.rdb"), algo = "sha256"),
            "153db129f2911ff008dbeda2a77b8d229c0816c726d2d9e368e71b7be9e32044"))
    store <- normalizePath(member$capacity$state_ownership$store_directory, mustWork = TRUE)
    stopifnot(startsWith(store, paste0(normalizePath(member_dir), "/")))
    list(member_dir = member_dir, member = member, store = store,
        head_sha256 = digest::digest(file = file.path(store, "dia-state-head.rds"), algo = "sha256"))
}

test_that("isolated rollup and session performance members retain complete operation boundaries", {
    mode <- Sys.getenv("MULTISCHOLAR_DIA025_MECHANISM_MEMBER")
    pilot <- identical(mode, "gui-mechanism-preflight")
    selected <- grepl("^gui-(rollup|session)-performance:(control|candidate)$", mode)
    skip_if(!pilot && !selected, "explicit isolated performance or preflight selector only")
    root <- dirname(Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT"))
    library <- normalizePath(Sys.getenv("MULTISCHOLAR_DIA008_LIBRARY"))
    rdb <- digest::digest(file = file.path(library, "MultiScholaR", "R", "MultiScholaR.rdb"), algo = "sha256")
    if (pilot) {
        source_dir <- file.path(dirname(root), "browser-smoke.vOSqsa")
        member <- jsonlite::read_json(file.path(source_dir, "smoke-receipt.json"))
        source <- list(member_dir = source_dir, member = member,
            store = member$capacity$state_ownership$store_directory)
        source$head_sha256 <- digest::digest(file = file.path(source$store, "dia-state-head.rds"), algo = "sha256")
        selections <- c("rollup", "session")
    } else {
        source <- readDIAHistoryFixture(dirname(root))
        selections <- sub("^gui-(rollup|session)-performance:.*$", "\\1", mode)
        candidate <- grepl(":candidate$", mode)
        expected <- if (!candidate) "8769e614e8b8126c3ed1ccaa6d491fbf857ffad65c988fa595e3cd33f2301a76" else
            if (selections == "rollup") "2aead2cd7912c9303331e2a5a0bb67f599583b40871745643e0b88fdaf0c8f7c" else
                "b6c7ef9aff9546d2dbb758931be6d7095277dc613919ca4f50e7dd9190f6c187"
        if (candidate && selections == "session") expected <- c(expected,
            "74a7dc9b50821e49e1a6d8c9a4366bdcc879997d61d06a9876acd55387404a08")
        stopifnot(rdb %in% expected)
    }
    receipt <- list(status = "running", mode = mode, source_member = source$member_dir,
        source_head_sha256 = source$head_sha256, input_sha256 = source$member$input_sha256,
        library = library, rdb_sha256 = rdb, profiling = FALSE, performance_acceptance = FALSE,
        r_version = R.version.string,
        namespace_versions = as.list(stats::setNames(vapply(sort(loadedNamespaces()),
            \(name) as.character(getNamespaceVersion(name)), character(1)), sort(loadedNamespaces()))),
        boundary = "fresh installed service includes startup, owned store copy/reopen, operation and teardown; oracle separate")
    withr::defer(jsonlite::write_json(receipt, file.path(root, "history-mechanism.json"),
        auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA))
    driver <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "workflow_checkpoint_dia_driver.R"), driver)
    for (mechanism in selections) {
        destination <- file.path(root, paste0(mechanism, "-store"))
        output <- file.path(root, paste0(mechanism, "-outputs"))
        stopifnot(dir.create(destination), dir.create(output))
        started <- proc.time()[["elapsed"]]
        manager <- copyDIAHistoryFixture(source$store, destination)
        withr::defer(closeDIAStateAdapter(manager))
        receipt$copy_reopen_seconds <- proc.time()[["elapsed"]] - started
        receipt$store_directory <- manager$.context$store$directory
        started <- proc.time()[["elapsed"]]
        result <- runDIAIsolatedMechanism(manager, mechanism,
            native_session = mechanism == "session" && (pilot || candidate),
            source_dir = source$member_dir, output = output)
        receipt$operation_seconds <- proc.time()[["elapsed"]] - started
        receipt$result <- result
        receipt$ownership <- diaStateOwnership(manager)
        expect_identical(receipt$ownership$working_records, 0L)
        expect_identical(receipt$ownership$native_records, 0L)
        expect_identical(receipt$ownership$store$open_resources, 0L)
        expect_equal(receipt$ownership$store$bytes, 0)
        if (pilot) {
            if (mechanism == "rollup") {
                identity <- readSavedStoreIdentity(file.path(source$store, "store.json"))
                original <- newDIAStateAdapter(newWorkflowEntityStore(dirname(source$store),
                    identity$project_id, identity$workflow_id), reopen = TRUE)
                withr::defer(closeDIAStateAdapter(original))
                checkDIAHistoryFixtureRecord(getDIAStateRecord(original, "precursor_rollup"),
                    getDIAStateRecord(manager, "precursor_rollup"), source$member_dir,
                    c(root, source$member_dir), driver$comparableCheckpointGuiRecord)
            } else {
                restored <- openDIANativeSession(result$outputs$sessionFilepath, output)
                expect_identical(restored$manager$getState(), manager$getState())
                closeDIAStateAdapter(restored$manager)
            }
        }
        closeDIAStateAdapter(manager)
    }
    expect_identical(digest::digest(file = file.path(source$store, "dia-state-head.rds"), algo = "sha256"), source$head_sha256)
    receipt$status <- if (pilot) "passed_mechanism_preflight" else "complete_pending_independent_verification"
})

test_that("history replay copies are independent through completion and failed ownership checks", {
    skip_if(nzchar(Sys.getenv("MULTISCHOLAR_DIA025_MECHANISM_MEMBER")),
        "small lifecycle is checked separately from measured members")
    skip_if_not(exists("writeProteinIqInput", asNamespace("MultiScholaR"), inherits = FALSE),
        "withdrawn DIA-025 combined candidate; retained for its installed trial")
    source <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "replay", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(source))
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    source$saveState("first", object, list(), "first")
    head <- source$.context$head
    destination <- withr::local_tempdir()
    copy <- copyDIAHistoryFixture(source$.context$store$directory, destination)
    withr::defer(closeDIAStateAdapter(copy))
    expect_identical(copy$getState(), source$getState())
    copy$saveState("second", object, list(), "second")
    expect_identical(source$current_state, "first")
    expect_identical(source$.context$head, head)
    expect_identical(readRDS(source$.context$head_path), head)
    expect_error(copyDIAHistoryFixture(source$.context$store$directory, destination), "fresh independent")
    expect_identical(diaStateOwnership(copy)$working_records, 0L)
    expect_identical(diaStateOwnership(copy)$store$open_resources, 0L)
    driver <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "workflow_checkpoint_dia_driver.R"), driver)
    record <- getDIAStateRecord(source, "first")
    expect_null(checkDIAHistoryFixtureRecord(record, record, "source", "copy",
        driver$comparableCheckpointGuiRecord))
    changed <- record
    changed$data@protein_quant_table[[2L]][[1L]] <- 999999
    expect_error(checkDIAHistoryFixtureRecord(record, changed, "source", "copy",
        driver$comparableCheckpointGuiRecord), "historical state differs")
    files <- stats::setNames(file.path(destination, c("peptide_values_imputed.tsv", "iq_output_file.txt")),
        c("peptide_input_file", "iq_output_file"))
    for (file in files) writeLines("small verifier fixture", file)
    record$config <- c(as.list(files), list(scientific_parameter = 0.01))
    mapped <- normalizeDIAIqReplayFiles(record, files)
    expect_identical(mapped$data, record$data)
    expect_identical(mapped$config$scientific_parameter, 0.01)
    expect_identical(mapped$config$iq_output_file, "<IQ_FILE>/iq_output_file.txt")
    record$config$iq_output_file <- files[[1L]]
    expect_error(normalizeDIAIqReplayFiles(record, files), "file location differs")
})

test_that("recorded-history verification compares science outside measured members", {
    mode <- Sys.getenv("MULTISCHOLAR_DIA025_MECHANISM_MEMBER")
    first_qc <- grepl("^gui-first-qc-verify:[A-Za-z0-9]{6}:[A-Za-z0-9]{6}$", mode)
    skip_if(!first_qc && !grepl("^gui-history-verify:[A-Za-z0-9]{6}(:[A-Za-z0-9]{6})?$", mode),
        "explicit independent history verifier only")
    root <- dirname(Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT"))
    identifiers <- strsplit(mode, ":", fixed = TRUE)[[1L]][-1L]
    member_dir <- file.path(dirname(root), paste0("browser-smoke.", identifiers[[1L]]))
    driver <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "workflow_checkpoint_dia_driver.R"), driver)
    mechanism_path <- file.path(member_dir, "history-mechanism.json")
    mechanism <- file.exists(mechanism_path)
    session <- NULL
    performance_mechanism <- FALSE
    if (mechanism) {
        member <- jsonlite::read_json(mechanism_path)
        stopifnot(identical(member$status, "complete_pending_independent_verification"),
            identical(member$library, normalizePath(Sys.getenv("MULTISCHOLAR_DIA008_LIBRARY"))))
        source <- readDIAHistoryFixture(dirname(root))
        stopifnot(identical(member$source_head_sha256, source$head_sha256),
            identical(member$input_sha256, source$member$input_sha256))
        source_dir <- source$member_dir
        expected_store <- source$store
        actual_store <- member$store_directory
        performance_mechanism <- grepl("^gui-(rollup|session)-performance:(control|candidate)$", member$mode)
    } else {
        stopifnot(length(identifiers) == 2L)
        source_dir <- file.path(dirname(root), paste0("browser-smoke.", identifiers[[2L]]))
        member <- jsonlite::read_json(file.path(member_dir, "smoke-receipt.json"))
        control <- jsonlite::read_json(file.path(source_dir, "smoke-receipt.json"))
        stopifnot(identical(member$status, if (first_qc) "passed_first_qc_ui_outputs_revert" else "passed_summary_ui_only"),
            identical(control$status, "passed_summary_ui_only"),
            identical(member$capacity$treatment, "qs2_duckdb"),
            identical(control$capacity$treatment, "qs2_duckdb"),
            identical(member$input_sha256, control$input_sha256),
            identical(member$installation$library_path, normalizePath(Sys.getenv("MULTISCHOLAR_DIA008_LIBRARY"))))
        expected_store <- control$capacity$state_ownership$store_directory
        actual_store <- member$capacity$state_ownership$store_directory
    }
    stopifnot(startsWith(normalizePath(expected_store), paste0(normalizePath(source_dir), "/")),
        startsWith(normalizePath(actual_store), paste0(normalizePath(member_dir), "/")))
    read_only_comparison <- performance_mechanism || !mechanism
    head_paths <- c(file.path(expected_store, "dia-state-head.rds"), file.path(actual_store, "dia-state-head.rds"))
    head_hashes <- vapply(head_paths, \(path) digest::digest(file = path, algo = "sha256"), character(1))
    if (read_only_comparison) {
        # Oracle getters are read-only; avoid copying two entire histories just
        # to compare their immutable records one at a time.
        identity <- readSavedStoreIdentity(file.path(expected_store, "store.json"))
        expected <- newDIAStateAdapter(newWorkflowEntityStore(dirname(expected_store),
            identity$project_id, identity$workflow_id), reopen = TRUE)
        identity <- readSavedStoreIdentity(file.path(actual_store, "store.json"))
        actual <- newDIAStateAdapter(newWorkflowEntityStore(dirname(actual_store),
            identity$project_id, identity$workflow_id), reopen = TRUE)
    } else {
        expected_root <- file.path(root, "expected-store")
        actual_root <- file.path(root, "actual-store")
        stopifnot(dir.create(expected_root), dir.create(actual_root))
        expected <- copyDIAHistoryFixture(expected_store, expected_root)
        actual <- copyDIAHistoryFixture(actual_store, actual_root)
    }
    withr::defer(closeDIAStateAdapter(expected))
    withr::defer(closeDIAStateAdapter(actual))
    states <- if (first_qc) c("initial", "raw_data_s4", "qvalue_filtered") else
        names(expected$.context$snapshot$states)
    iq_expected_files <- iq_actual_files <- NULL
    if (performance_mechanism && startsWith(member$mode, "gui-rollup-performance:")) {
        states <- "precursor_rollup"
    } else if (performance_mechanism) {
        if (isTRUE(member$result$native_session)) {
            archive <- openDIANativeSession(member$result$outputs$sessionFilepath, root)
            closeDIAStateAdapter(actual)
            actual <- archive$manager
            expect_identical(actual$getHistory(), expected$getHistory())
        } else {
            session <- validateExternalWorkflowLegacySession(member$result$outputs$latestFilepath, states)
            expect_identical(digest::digest(file = member$result$outputs$latestFilepath, algo = "sha256"),
                digest::digest(file = member$result$outputs$sessionFilepath, algo = "sha256"))
        }
    } else if (mechanism && identical(member$mode, "gui-export-history-member")) {
        session <- validateExternalWorkflowLegacySession(member$outputs$latestFilepath, states)
        stopifnot(identical(digest::digest(file = member$outputs$latestFilepath, algo = "sha256"),
            digest::digest(file = member$outputs$sessionFilepath, algo = "sha256")))
        stopifnot(identical(session$r6_state_history, expected$state_history))
    } else if (mechanism) {
        stopifnot(identical(member$mode, "gui-iq-history-member"))
        states <- "protein_s4_created"
        iq_expected_files <- iq_actual_files <- list()
        for (file in c("peptide_values_imputed.tsv", "iq_output_file.txt")) {
            section <- if (startsWith(file, "peptide")) "peptide_qc" else "protein_qc"
            original <- file.path(source_dir, "dia_baseline_smoke", "results", "proteomics", section, file)
            copied <- file.path(member_dir, "outputs", file)
            stopifnot(identical(digest::digest(file = original, algo = "sha256"),
                digest::digest(file = copied, algo = "sha256")))
            field <- if (startsWith(file, "peptide")) "peptide_input_file" else "iq_output_file"
            iq_expected_files[[field]] <- original
            iq_actual_files[[field]] <- copied
        }
    } else if (first_qc) {
        stopifnot(identical(names(actual$.context$snapshot$states), states),
            identical(actual$getHistory(), c("initial", "raw_data_s4")),
            identical(actual$current_state, "raw_data_s4"))
    } else {
        stopifnot(identical(names(actual$.context$snapshot$states), states),
            identical(actual$state_history, expected$state_history))
    }
    if (!mechanism && !first_qc && usesDIAReferenceStates(actual)) {
        expected_config <- getDIAStateRecord(expected, "protein_s4_created")$config
        actual_config <- readDIAStateMetadata(actual, getDIAStateReference(actual, "protein_s4_created"))$config
        iq_expected_files <- iq_actual_files <- list()
        for (field in c("peptide_input_file", "iq_output_file")) {
            original <- normalizePath(expected_config[[field]], mustWork = TRUE)
            staged <- normalizePath(actual_config[[field]], mustWork = TRUE)
            stopifnot(startsWith(original, paste0(normalizePath(source_dir), "/")),
                startsWith(staged, paste0(normalizePath(member_dir), "/")),
                identical(basename(original), basename(staged)))
            section <- if (field == "peptide_input_file") "peptide_qc" else "protein_qc"
            published <- file.path(member_dir, "dia_baseline_smoke", "results", "proteomics", section, basename(staged))
            hashes <- vapply(c(original, staged, published), function(path) digest::digest(file = path, algo = "sha256"), character(1))
            stopifnot(length(unique(hashes)) == 1L)
            iq_expected_files[[field]] <- original
            iq_actual_files[[field]] <- staged
        }
    }
    differences <- character()
    for (name in states) {
        observed <- if (is.null(session)) getDIAStateRecord(actual, name) else session$r6_complete_states[[name]]
        expected_record <- getDIAStateRecord(expected, name)
        if (!is.null(iq_expected_files) && name == "protein_s4_created") {
            expected_record <- normalizeDIAIqReplayFiles(expected_record, iq_expected_files)
            observed <- normalizeDIAIqReplayFiles(observed, iq_actual_files)
        }
        tryCatch(checkDIAHistoryFixtureRecord(expected_record, observed,
            source_dir, c(member_dir, source_dir), driver$comparableCheckpointGuiRecord),
            error = function(error) differences <<- c(differences, paste("state", name, ":", conditionMessage(error))))
        observed <- expected_record <- NULL
        gc(full = TRUE) # Independent verification only, never application timing.
    }
    if (length(differences)) stop(paste(differences, collapse = "\n"), call. = FALSE)
    expect_gt(length(states), 0L)
    if (first_qc) {
        find_pngs <- function(directory) {
            paths <- list.files(file.path(directory, "dia_baseline_smoke", "results", "proteomics",
                "publication_graphs", "filtering_qc"), pattern = "^2_qval_filtered_.*[.]png$",
                recursive = TRUE, full.names = TRUE)
            paths <- paths[!grepl("/.workflow-export/", paths, fixed = TRUE)]
            stopifnot(length(paths) == 6L, !anyDuplicated(basename(paths)))
            stats::setNames(vapply(paths, function(path) digest::digest(file = path, algo = "sha256"),
                character(1)), basename(paths))
        }
        expect_identical(find_pngs(member_dir), find_pngs(source_dir))
    }
    if (read_only_comparison) expect_identical(vapply(head_paths,
        \(path) digest::digest(file = path, algo = "sha256"), character(1)), head_hashes)
    if (performance_mechanism) expect_identical(digest::digest(
        file = file.path(expected_store, "dia-state-head.rds"), algo = "sha256"), member$source_head_sha256)
    jsonlite::write_json(list(status = "passed_history_science", member = member_dir,
        reference = source_dir, states = states, tolerance = 1e-12,
        normalized_da_ui_clock = "data@args$daAnalysisUI$timestamp value only; presence/type preserved",
        read_only_comparison = read_only_comparison,
        complete_legacy_restored = !is.null(session), performance_acceptance = FALSE),
        file.path(root, "history-verification.json"), auto_unbox = TRUE, pretty = TRUE)
})

test_that("recorded large export or IQ member measures the actual saved-history mechanism", {
    mode <- Sys.getenv("MULTISCHOLAR_DIA025_MECHANISM_MEMBER")
    skip_if(!mode %in% c("gui-export-history-member", "gui-iq-history-member"),
        "explicit guarded recorded-history member only")
    receipt_path <- Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT")
    expect_true(nzchar(receipt_path))
    root <- dirname(receipt_path)
    source <- readDIAHistoryFixture(dirname(root))
    diagnostic <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "workflow_checkpoint_gui_diagnostic.R"), diagnostic)
    observation <- diagnostic$startCheckpointGuiDiagnostic(root)
    withr::defer(diagnostic$stopCheckpointGuiDiagnostic(observation))
    event <- function(phase) diagnostic$checkpointGuiDiagnosticEvent(observation, phase)
    receipt <- list(status = "running", mode = mode, source_member = source$member_dir,
        source_head_sha256 = source$head_sha256, input_sha256 = source$member$input_sha256,
        library = normalizePath(Sys.getenv("MULTISCHOLAR_DIA008_LIBRARY")),
        performance_acceptance = FALSE,
        boundary = "fresh installed service includes copy/reopen/preparation; operation includes outputs; independent readback")
    withr::defer(jsonlite::write_json(receipt, file.path(root, "history-mechanism.json"),
        auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA))
    destination <- file.path(root, "fixture-store")
    expect_true(dir.create(destination))
    event("fixture:before_copy")
    started <- proc.time()[["elapsed"]]
    manager <- copyDIAHistoryFixture(source$store, destination)
    withr::defer(closeDIAStateAdapter(manager))
    receipt$copy_reopen_seconds <- proc.time()[["elapsed"]] - started
    receipt$store_directory <- manager$.context$store$directory
    receipt$copied_bytes <- sum(file.info(list.files(source$store, recursive = TRUE,
        all.files = TRUE, no.. = TRUE, full.names = TRUE))$size)
    event("fixture:after_copy")
    directory <- file.path(root, "outputs")
    expect_true(dir.create(directory))
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    export_member <- identical(mode, "gui-export-history-member")
    if (export_member) {
        descriptors <- list.files(manager$.context$store$directory, "^descriptor[.]rds$",
            recursive = TRUE, full.names = TRUE)
        expect_true(all(file.size(descriptors) <= 65536))
        descriptors <- lapply(descriptors, readRDS)
        selected <- Filter(function(value) identical(value$identity$parameters, "export-metadata"), descriptors)
        expect_length(selected, 1L)
        metadata <- resolveDIAStateValue(manager$.context, selected[[1L]])
        for (field in c("contrasts_tbl", "design_matrix", "config_list", "fasta_metadata",
            "accession_cleanup_results", "ruv_optimization_result", "qc_params", "protein_counts",
            "mixed_species_analysis")) workflow[[field]] <- metadata[[field]]
        event("operation:before_export")
        started <- proc.time()[["elapsed"]]
        result <- withDIAStateScope(manager, function() {
            norm <- list(correlation_filtered_obj = manager$getState(), best_k = metadata$ruv_k,
                correlation_threshold = metadata$correlation_threshold)
            runProtNormExportSessionWorkflow(workflow, norm,
                list(norm_method = metadata$normalization_method, ruv_mode = metadata$ruv_mode), directory,
                withProgressFn = function(message, value, expr) force(expr), incProgressFn = function(...) NULL)
        })
        receipt$operation_seconds <- proc.time()[["elapsed"]] - started
        receipt$outputs <- result$exportArtifacts
        receipt$shared_tables <- manager$.context$export_shared_tables
        receipt$shared_columns <- manager$.context$export_shared_columns
        receipt$shared_support_tables <- manager$.context$export_shared_support_tables
        expect_null(result$sessionData)
        expect_identical(receipt$outputs$publication$status, "complete")
        expect_identical(receipt$outputs$publication$serialization_count, 1L)
        event("operation:after_export")
    } else {
        saveRDS(manager$.context$snapshot$states$protein_s4_created,
            file.path(root, "expected-iq-reference.rds"))
        # Metadata selection only; input decoding remains inside the measured Apply.
        manager$current_state <- "imputed"
        history <- manager$state_history
        manager$state_history <- history[seq_len(match("imputed", unlist(history)))]
        event("operation:before_iq")
        started <- proc.time()[["elapsed"]]
        withr::with_dir(directory, withDIAStateScope(manager, function() {
            result <- runProteinIqRollupApplyStep(workflow,
                list(peptide_qc_dir = directory, protein_qc_dir = directory),
                showNotificationFn = function(...) NULL)
            expect_s4_class(result$proteinObj, "ProteinQuantitativeData")
            invisible(NULL)
        }))
        receipt$operation_seconds <- proc.time()[["elapsed"]] - started
        receipt$outputs <- as.list(file.path(directory,
            c("peptide_values_imputed.tsv", "iq_output_file.txt", "qc-plots.pdf")))
        expect_true(all(file.exists(unlist(receipt$outputs))))
        event("operation:after_iq")
    }
    receipt$ownership <- diaStateOwnership(manager)
    expect_identical(receipt$ownership$working_records, 0L)
    expect_identical(receipt$ownership$native_records, 0L)
    expect_identical(receipt$ownership$store$open_resources, 0L)
    expect_equal(receipt$ownership$store$bytes, 0)
    expect_identical(digest::digest(file = file.path(source$store, "dia-state-head.rds"), algo = "sha256"),
        source$head_sha256)
    receipt$status <- "complete_pending_independent_verification"
})
test_that("fresh worker or export member records exact science and complete operation cost", {
    mode <- Sys.getenv("MULTISCHOLAR_DIA025_MECHANISM_MEMBER")
    skip_if(!mode %in% c("gui-worker-member", "gui-export-member"), "explicit guarded member only")
    receipt <- Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT")
    expect_true(nzchar(receipt))
    root <- dirname(receipt)
    export_member <- identical(mode, "gui-export-member")
    recipe <- readWorkflowCheckpointRecipe(if (export_member) "representative" else "small")
    object <- makeWorkflowCheckpointProtein(recipe, allow_large = export_member,
        output_root = root)
    input_hash <- digest::digest(object, algo = "sha256")
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "member", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    workflow$design_matrix <- object@design_matrix
    workflow$contrasts_tbl <- data.frame(contrasts = "grouptreated-groupcontrol")
    paths <- list(source_dir = file.path(root, "outputs"), protein_qc_dir = file.path(root, "outputs"))
    dir.create(paths$source_dir)
    serializations <- 0L
    if (export_member) {
        for (index in seq_len(recipe$history_depth)) {
            object@args$export_revision <- index
            if (index %% 5L == 0L) object@protein_quant_table[[2L]] <- object@protein_quant_table[[2L]] + 0.01
            manager$saveState(if (index == recipe$history_depth) "correlation_filtered" else
                paste0("step_", index), object, list(), "component export fixture")
        }
        norm <- list(correlation_filtered_obj = object, best_k = NA_integer_)
    } else manager$saveState("raw", object, list(), "source")
    started <- proc.time()[["elapsed"]]
    if (export_member) {
        writer <- function(value, path) {
            if (is.list(value) && "r6_complete_states" %in% names(value))
                serializations <<- serializations + 1L
            saveRDS(value, path)
        }
        export_args <- list(workflowData = workflow, normData = norm,
            input = list(norm_method = "none", ruv_mode = "skip"), sourceDir = paths$source_dir,
            withProgressFn = function(message, value, expr) force(expr), incProgressFn = function(...) NULL,
            saveExportArtifactsFn = function(...) saveProtNormExportArtifacts(..., saveRdsFn = writer))
        if ("serializeSessionFn" %in% names(formals(runProtNormExportSessionWorkflow))) {
            export_args$serializeSessionFn <- writer
        }
        result <- do.call(runProtNormExportSessionWorkflow, export_args)
        session <- readRDS(result$exportArtifacts$latestFilepath)
        expect_identical(session$current_s4_object, object)
        hashes <- lapply(session$r6_complete_states, function(record) {
            if (is.null(record$data)) return(NULL)
            digest::digest(list(record$data@protein_quant_table, record$data@design_matrix,
                record$data@protein_id_table, record$data@args), algo = "sha256")
        })
    } else {
        result <- runProteinReplicateFilterApplyStep(workflow, paths, "group", 4L,
            isTestModeFn = function() FALSE)
        expected <- removeProteinsWithOnlyOneReplicate(object, NA, "group")
        expect_identical(result$filteredS4@protein_quant_table, expected@protein_quant_table)
        expect_identical(result$filteredS4@design_matrix, expected@design_matrix)
        hashes <- digest::digest(list(result$filteredS4@protein_quant_table,
            result$filteredS4@design_matrix), algo = "sha256")
    }
    elapsed <- proc.time()[["elapsed"]] - started
    heap <- gc()
    info <- diaStateOwnership(manager)
    expect_identical(info$working_records, 0L)
    expect_identical(info$store$open_resources, 0L)
    jsonlite::write_json(list(mode = mode, input_sha256 = input_hash,
        science_sha256 = hashes, serialization_count = serializations,
        operation_and_readback_seconds = elapsed, vector_heap_high_water_mib = unname(heap["Vcells", 6L]),
        shared_export_tables = manager$.context$export_shared_tables,
        ownership = info, boundary = "fresh installed service; operation includes outputs and readback; service also includes preparation/startup",
        performance_acceptance = FALSE), file.path(root, "mechanism-member.json"),
        auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
})

test_that("retained rollup attribution separates bounded read and kernel work", {
    mode <- Sys.getenv("MULTISCHOLAR_DIA025_MECHANISM_MEMBER")
    preflight <- identical(mode, "gui-rollup-profile-preflight")
    comparison <- grepl("^gui-rollup-compare:(control|candidate)$", mode)
    skip_if(!comparison && !mode %in% c("gui-rollup-profile", "gui-rollup-profile-preflight"),
        "explicit retained rollup diagnostic only")
    output <- dirname(Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT"))
    library_path <- Sys.getenv("MULTISCHOLAR_DIA008_LIBRARY")
    if (!identical(mode, "gui-rollup-compare:candidate")) expect_identical(digest::digest(file = file.path(library_path, "MultiScholaR", "R", "MultiScholaR.rdb"),
        algo = "sha256"), "c55102a392a3384021f979af0536d734a2863350a5321472ffecad9cbfd7f7c8")
    if (preflight) {
        fixture <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
        fixture_root <- withr::local_tempdir()
        fixture_store <- newWorkflowEntityStore(fixture_root, "profile-preflight", "DIA")
        reference <- saveWorkflowTable(fixture_store, fixture@peptide_data, "source",
            newSavedWorkflowScientificIdentity("profile", "small", "preflight", "all", "frozen"),
            ".rownames", metadata_version = 3L)
        fixture_scope <- newWorkflowTableReadScope(fixture_store, reference, reader = "arrow")
        for (key in c("Protein.Group", "Run")) saveWorkflowGroupIndex(fixture_scope, key)
        releaseWorkflowTableReadScope(fixture_scope)
        closeWorkflowEntityStore(fixture_store)
        directory <- file.path(fixture_root, ".workflow-entities")
        request_id <- "source"
    } else {
        member <- file.path(dirname(output), "browser-smoke.cuurUw")
        state <- jsonlite::read_json(file.path(member, "gui-state-latest.json"))
        directory <- state$store_directory
        request_id <- "dia-data-326b617be1c6a8"
    }
    owner <- jsonlite::read_json(file.path(directory, "store.json"))
    source_paths <- list.files(directory, recursive = TRUE, full.names = TRUE, all.files = TRUE)
    source_paths <- source_paths[!dir.exists(source_paths)]
    before <- file.info(source_paths)[c("size", "mtime")]
    store <- newWorkflowEntityStore(dirname(directory), owner$project_id, owner$workflow_id)
    withr::defer(closeWorkflowEntityStore(store))
    # All integrity checks remain enabled; their log belongs to this diagnostic.
    testthat::local_mocked_bindings(hashWorkflowOwnedFile = function(store, path)
        hashWorkflowFile(path, file.path(output, "source-validation.log")),
        releaseWorkflowStoreFilePages = function(store, paths) {
            releaseWorkflowFileCache(paths, file.path(output, "source-release.log"))
            invisible(NULL)
        },
        saveWorkflowConsumerEntity = function(...) stop("diagnostic attempted to mutate frozen evidence"),
        saveWorkflowTable = function(...) stop("diagnostic attempted to mutate frozen evidence"))
    for (request in list.dirs(file.path(directory, "entities"), recursive = FALSE, full.names = FALSE))
        readCommittedWorkflowRequest(store, request, payload = FALSE)
    reference <- readCommittedWorkflowRequest(store, request_id, payload = FALSE)
    stopifnot(!is.null(reference), preflight || reference$metadata$shape[[1L]] == 24000000)
    receipt <- list(status = "running", source = reference$entity_id, source_sha256 = reference$payload$sha256,
        cases = list(), phases = list(), mode = mode,
        runtime_rdb_sha256 = digest::digest(file = file.path(library_path, "MultiScholaR", "R", "MultiScholaR.rdb"), algo = "sha256"),
        profile_scope = if (comparison) "disabled" else "first complete-sample kernel only",
        comparison_scope = "bounded selected groups; not full-workflow timing")
    flush <- function() jsonlite::write_json(receipt, file.path(output, "rollup-attribution.json"),
        auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
    flush()
    gc.time(on = TRUE)
    measure <- function(label, expr) {
        started <- proc.time()
        before_gc <- gc.time()
        receipt$active_phase <<- label
        flush()
        value <- tryCatch(force(expr), error = function(error) {
            receipt$status <<- "failed"
            receipt$error <<- conditionMessage(error)
            flush()
            stop(error)
        })
        elapsed <- proc.time() - started
        receipt$phases[[length(receipt$phases) + 1L]] <<- list(label = label,
            elapsed = unname(elapsed[["elapsed"]]), cpu = unname(sum(elapsed[c("user.self", "sys.self")])),
            gc_cpu = sum(gc.time() - before_gc))
        flush()
        value
    }
    scope <- measure("validate_source", newWorkflowTableReadScope(store, reference, reader = "arrow"))
    withr::defer(releaseWorkflowTableReadScope(scope))
    original_stream <- streamSavedWorkflowTable
    testthat::local_mocked_bindings(streamSavedWorkflowTable = function(scope, columns, consume, ...) {
        label <- paste0(receipt$active_phase, if (identical(scope$reference, reference)) ":assay" else ":index")
        callbacks <- 0
        rows <- 0
        started <- proc.time()[["elapsed"]]
        value <- original_stream(scope, columns, function(data, offset, row_names) {
            at <- proc.time()[["elapsed"]]
            consume(data, offset, row_names)
            callbacks <<- callbacks + proc.time()[["elapsed"]] - at
            rows <<- rows + nrow(data)
        }, ...)
        total <- proc.time()[["elapsed"]] - started
        receipt$phases[[length(receipt$phases) + 1L]] <<- list(label = label,
            stream_elapsed = total, callback_assembly_elapsed = callbacks,
            decoding_restoration_elapsed = total - callbacks, rows = rows)
        flush()
        value
    })
    for (iteration in seq_len(if (comparison) 3L else 1L)) for (key in c("Protein.Group", "Run")) {
        index <- measure(paste0(key, ":index"), {
            if (!comparison) saveWorkflowGroupIndex(scope, key) else {
                # Compare readers/kernels against the same frozen index artifacts.
                parameters <- digest::digest(list(key, scope$metadata$column_metadata[[key]]), algo = "sha256")
                records <- Filter(function(record) identical(record$kind, "scientific_product") &&
                    identical(record$identity$input, reference$entity_id) &&
                    identical(record$identity$selectors, "stable-group-ordinals-v1") &&
                    identical(record$identity$parameters, parameters), as.list(store$records, all.names = TRUE))
                stopifnot(length(records) == 1L)
                value <- hydrateSavedWorkflowEntity(store, records[[1L]])$value
                stopifnot(identical(value$source_id, reference$entity_id), identical(value$group, key))
                list(reference = records[[1L]], value = value)
            }
        })
        groups <- readWorkflowIndexDirectory(store, index$value$groups)
        selected <- unique(c(1L, ceiling(nrow(groups) / 2), nrow(groups)))
        columns <- c("Run", "Protein.Group", "Protein.Ids", "Stripped.Sequence", "Modified.Sequence",
            "Precursor.Id", "Precursor.Quantity", "Precursor.Normalised",
            "identification_peptide_count", "identification_peptidoform_count",
            "peptides_for_protein_count", "peptidoforms_for_protein_count")
        columns <- intersect(columns, scope$metadata$columns)
        for (i in selected) {
            label <- paste(key, i, iteration, sep = ":")
            part <- measure(paste0(label, ":gather"), readWorkflowGroupCohort(scope, columns,
                index$value$positions, groups[i, , drop = FALSE]))
            if (key == "Run") {
                annotations <- unique(part[c("Protein.Group", "Protein.Ids")])
                stopifnot(!anyDuplicated(annotations$Protein.Group))
                part$Protein.Ids <- NULL
                attr(part, "workflow_source_order") <- NULL
                profiled <- !comparison && i == selected[[1L]]
                profile <- file.path(output, "rollup-kernel.Rprof")
                if (profiled) utils::Rprof(profile, interval = 0.02, gc.profiling = TRUE)
                result <- tryCatch(measure(paste0(label, ":kernel"), rollUpPrecursorToPeptideHelper(part,
                    sample_id_column = "Run", protein_id_column = "Protein.Group",
                    peptide_sequence_column = "Stripped.Sequence", modified_peptide_sequence_column = "Modified.Sequence",
                    precursor_quantity_column = "Precursor.Quantity", precursor_normalised_column = "Precursor.Normalised",
                    precursor_id_column = "Precursor.Id", precursor_charge_column = "Precursor.Charge",
                    core_utilisation = NA, is_logged_data = FALSE, return_rollup_result = TRUE)),
                    finally = if (profiled) utils::Rprof(NULL))
                result$data$Protein.Ids <- annotations$Protein.Ids[match(result$data$Protein.Group, annotations$Protein.Group)]
                if (comparison) {
                    oracle <- readRDS(file.path(dirname(output), "browser-smoke.y0oLNN",
                        paste0("sample-", i, "-oracle.rds")))
                    expect_equal(result, oracle, tolerance = 1e-12)
                    expect_identical(result$summary, oracle$summary)
                }
                accumulator <- newDIAFilteringAccumulator(result$data[0, , drop = FALSE], "Protein.Group")
                measure(paste0(label, ":summary"), {
                    for (start in seq.int(1L, nrow(result$data), by = 10000L))
                        accumulator$push(result$data[seq.int(start, min(nrow(result$data), start + 9999L)), , drop = FALSE])
                    NULL
                })
                measure(paste0(label, ":write"), writeSavedWorkflowParquet(result$data, ".rownames",
                    file.path(output, paste0("sample-", i, if (comparison) paste0("-rep-", iteration), ".parquet")), metadata_version = 3L))
                saveRDS(result, file.path(output, paste0("sample-", i,
                    if (comparison) paste0("-rep-", iteration), "-oracle.rds")))
                counted <- accumulator$finish()
                counted$summary <- unclass(counted$summary)
                receipt$cases[[label]] <- list(input_rows = nrow(part), output_rows = nrow(result$data),
                    source_group = as.character(groups$value[[i]]), summary = counted)
                result <- accumulator <- NULL
            }
            part <- NULL
            flush()
        }
    }
    releaseWorkflowTableReadScope(scope)
    expect_identical(file.info(source_paths)[c("size", "mtime")], before)
    expect_identical(list.files(directory, recursive = TRUE, full.names = TRUE, all.files = TRUE)[
        !dir.exists(list.files(directory, recursive = TRUE, full.names = TRUE, all.files = TRUE))], source_paths)
    if (preflight) {
        expect_error(measure("injected_failure", stop("bounded diagnostic failure")), "bounded diagnostic failure")
        expect_identical(jsonlite::read_json(file.path(output, "rollup-attribution.json"))$status, "failed")
        receipt$failure_observed <- TRUE
        receipt$error <- NULL
    }
    receipt$status <- "complete"
    receipt$active_phase <- NULL
    receipt$source_unchanged <- TRUE
    if (!comparison) {
        stopifnot(file.size(file.path(output, "rollup-kernel.Rprof")) < 67108864)
        profile <- summaryRprof(file.path(output, "rollup-kernel.Rprof"))
        receipt$profile <- list(self = utils::head(profile$by.self, 15L), total = utils::head(profile$by.total, 15L))
    }
    flush()
})

#' Replay production scientific/publication actions from a committed raw checkpoint.
#' @param source_member Frozen small or 300-sample member.
#' @param root Fresh guarded-member output directory.
#' @return Bounded phase receipt, with no GUI-success claim.
#' @noRd
runDIARawCheckpointReplay <- function(source_member, root) {
    source <- jsonlite::read_json(file.path(source_member, "gui-state-latest.json"))
    source_store <- source$store_directory
    source_head <- digest::digest(file = file.path(source_store, "dia-state-head.rds"), algo = "sha256")
    project <- file.path(root, "dia_baseline_smoke")
    dir.create(project)
    directories <- setupDirectories(project, "proteomics")
    paths <- directories$proteomics
    directories$proteomics_dia_baseline_smoke <- paths
    state_root <- file.path(paths$source_dir, ".checkpoint-replay")
    dir.create(state_root)
    manager <- copyDIAReplayCheckpoint(source_store, state_root)
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    raw <- getDIAStateReference(manager, "raw_data_s4")
    manager$.context$reference_states <- TRUE
    manager$.context$session_format <- "native_v1"
    manager$.context$sql_rollup <- TRUE
    if (!identical(raw$payload$codec, "rds")) {
        # Explicit compatibility conversion of the one required small legacy record.
        old <- withDIAStateScope(manager, function() getDIAStateRecord(manager, "raw_data_s4"))
        manager$saveState("raw_data_s4", old$data, old$config, old$description)
        old <- NULL
    }
    selectDIAStateReference(manager, "raw_data_s4")
    record <- readDIAStateMetadata(manager)
    stopifnot(isDIAComponentCheckpoint(record$data))
    config <- record$config
    source_scripts <- file.path(source_member, "dia_baseline_smoke", "scripts", "proteomics")
    sidecars <- c("aa_seq_tbl_final", "uniprot_dat_cln", "fasta_metadata")
    annotation <- stats::setNames(lapply(sidecars, function(name) {
        path <- file.path(source_scripts, paste0(name, ".RDS"))
        stopifnot(file.exists(path), file.size(path) < 67108864)
        readRDS(path)
    }), sidecars)
    contrasts <- readProtDesignImportedContrasts(file.path(source_scripts, "contrast_strings.tab"))
    for (name in c("config.ini", "contrast_strings.tab", "design_matrix.tab", paste0(sidecars, ".RDS")))
        stopifnot(file.copy(file.path(source_scripts, name), paths$source_dir, overwrite = TRUE))
    if (file.exists(file.path(source_scripts, "data_cln.tab")))
        stopifnot(file.link(file.path(source_scripts, "data_cln.tab"), file.path(paths$source_dir, "data_cln.tab")))
    receipt <- list(status = "running", source_member = basename(source_member), source_head_sha256 = source_head,
        input_rows = record$data$dimensions$rows, phases = list(), kind = "raw_checkpoint_replay")
    flush <- function() jsonlite::write_json(receipt, file.path(root, "checkpoint-replay.json"),
        auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
    flush()
    withr::local_options(multischolar.test_mode = FALSE, multischolar.dia_native_sessions = TRUE,
        multischolar.capture_test_checkpoints = FALSE)
    withr::local_envvar(MULTISCHOLAR_TEST_MODE = "false")
    rlang::local_bindings(config_list = config, project_dirs = directories,
        filtering_progress = methods::new("FilteringProgress"), .env = .GlobalEnv)
    grDevices::pdf(NULL)
    device <- grDevices::dev.cur()
    on.exit(if (device %in% grDevices::dev.list()) grDevices::dev.off(device), add = TRUE)
    tryCatch(shiny::testServer(function(input, output, session) {
        design <- readDIAExecutionComponents(manager, fields = "design_matrix")@design_matrix
        workflow <- shiny::reactiveValues(state_manager = manager, design_matrix = design,
            config_list = config, contrasts_tbl = contrasts, data_type = "peptide", data_format = "diann",
            aa_seq_tbl_final = annotation$aa_seq_tbl_final, uniprot_dat_cln = annotation$uniprot_dat_cln,
            fasta_metadata = annotation$fasta_metadata)
        manager$.context$experiment_paths <- paths
        manager$.context$reclaim_enabled <- TRUE
        attachDIAStateSession(manager, session, workflow, receipt_dir = root)
        norm <- shiny::reactiveValues()
        da <- shiny::reactiveValues()
        values <- shiny::reactiveValues()
        phase <- function(name, fn) {
            receipt$active_phase <<- name
            flush()
            started <- proc.time()
            shiny::isolate(withDIAStateScope(manager, fn))
            session$flushReact()
            ownership <- diaStateOwnership(manager)
            stopifnot(ownership$store$open_resources == 0L, ownership$working_records == 0L)
            receipt$phases[[length(receipt$phases) + 1L]] <<- list(name = name,
                elapsed = unname((proc.time() - started)[["elapsed"]]), state = manager$current_state,
                state_reference = getDIAStateReference(manager)$entity_id)
            flush()
            invisible(NULL)
        }
        publish <- function(result, label) {
            publishDIAProteinProgress(result, label, "proteomics", "dia_baseline_smoke", updateProteinFiltering)
            invisible(NULL)
        }
        session$userData$execute_replay <- function() {
            phase("qvalue", function() publish(runPeptideQvalueApplyStep(workflow, 0.01, 0.01, TRUE), "2_qval_filtered"))
            phase("evidence", function() publish(runProteinPeptideApplyStep(workflow, 1, 2), "3_protein_peptide_filtered"))
            phase("precursor_rollup", function() publish(runPeptideRollupApplyStep(workflow), "4_precursor_rollup"))
            phase("peptide_intensity", function() publish(runPeptideIntensityApplyStep(workflow, FALSE, 2, 2, 1), "5_intensity_filtered"))
            phase("sample_filter", function() publish(runPeptideSampleApplyStep(workflow, 500), "6_sample_filtered"))
            phase("peptide_replicates", function() publish(runPeptideReplicateApplyStep(workflow, "replicates"), "7_replicate_filtered"))
            phase("imputation", function() publish(runPeptideImputationStep(workflow, 0.5), "8_imputed"))
            phase("IQ_MaxLFQ", function() publish(runProteinIqRollupApplyStep(workflow, paths,
                processLongFormatFn = function(...) iq::process_long_format(...,
                    pdf_out = file.path(root, "qc-plots.pdf")),
                showNotificationFn = function(...) NULL), "9_protein_s4_created"))
            phase("protein_cleanup", function() publish(runProteinAccessionCleanupStep(workflow, ";", "mean"), "10_protein_accession_cleaned"))
            phase("protein_intensity", function() publish(runProteinIntensityFilterApplyStep(workflow, FALSE, 2, 2, 1), "11_protein_intensity_filtered"))
            phase("protein_duplicates", function() publish(runProteinDuplicateRemovalStep(workflow, "mean"), "12_duplicates_removed"))
            phase("protein_replicates", function() publish(runProteinReplicateFilterApplyStep(workflow, paths,
                "group", 8), "13_protein_replicate_filtered"))
            settings <- list(norm_method = "cyclicloess", ruv_mode = "automatic",
                auto_percentage_min = 1, auto_percentage_max = 20, separation_metric = "max_difference",
                k_penalty_weight = 0.5, adaptive_k_penalty = TRUE, color_variable = "group", shape_variable = NULL)
            aesthetics <- function() list(color_var = "group", shape_var = NULL)
            phase("normalization_RUV", function() {
                generateProtNormPreNormalizationQcArtifacts(manager, paths$protein_qc_dir, aesthetics())
                post <- function(reference) {
                    norm$qc_plot_paths <- generateProtNormPostNormalizationQcArtifacts(reference,
                        norm$publication_paths$protein_qc_dir %||% paths$protein_qc_dir, aesthetics(), norm$qc_plot_paths)
                }
                ruv <- function(reference) {
                    norm$qc_plot_paths <- generateProtNormRuvCorrectedQcArtifacts(reference,
                        norm$publication_paths$protein_qc_dir %||% paths$protein_qc_dir, aesthetics(), norm$qc_plot_paths)
                }
                stopifnot(runDIAReferenceNormalization(settings, workflow, norm, paths,
                    "proteomics", "dia_baseline_smoke", post, ruv, function() "group"))
            })
            phase("Pearson", function() {
                runDIAReferenceCorrelation(workflow, norm, paths, 0.5, FALSE,
                    "proteomics", "dia_baseline_smoke", function() "group", aesthetics)
            })
            phase("native_session", function() {
                runProtNormExportSessionWorkflow(workflow, norm, settings, paths$source_dir,
                    withProgressFn = function(message, value, expr) force(expr), incProgressFn = function(...) NULL)
            })
            phase("DA", function() {
                da$current_s4_object <- getDIAStateReference(manager)
                runDIAReferenceDaAnalysis(list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05,
                    treat_lfc_cutoff = 0), da, workflow, paths)
            })
            phase("DA_views", function() {
                manifest <- da$da_results_list
                contrast <- manifest$contrasts$friendly_names[[1L]]
                readDIAResultPage(manifest, contrast, "all", 0.05, 0, 0, 50)
                widget <- generateProtDAVolcanoPlotGlimma(manifest, contrast,
                    args_row_id = manifest$protein_id_column, output_dir = paths$da_output_dir)
                stopifnot(!is.null(widget))
                htmltools::renderTags(widget)
                widget <- NULL
                heatmap <- generateProtDAHeatmap(manifest, contrast, tree_cut_method = "k_clusters")
                stopifnot(!is.null(heatmap))
                save_heatmap_products(heatmap$plot, heatmap$row_clusters,
                    list(contrast = contrast, top_n = 50, clustering_method = "ward.D2",
                        distance_method = "euclidean", cluster_rows = TRUE, cluster_cols = TRUE,
                        scaling = "row", color_scheme = "RdBu", tree_cut_method = "k_clusters", n_clusters = 4),
                    paths$publication_graphs_dir, gsub("[^A-Za-z0-9_]", "_", paste0("prot_", contrast)))
            })
            phase("Summary", function() {
                completeProtSummaryWorkflowArgsSave(output, values, directories, workflow,
                    experimentLabel = "dia_baseline_smoke", description = "DIA IQ checkpoint replay",
                    showNotificationFn = function(...) NULL)
                stopifnot(isTRUE(values$workflow_args_saved))
                runProtSummaryPublicationCopy(output, values, directories,
                    experimentLabel = "dia_baseline_smoke", description = "DIA IQ checkpoint replay",
                    contrastsTbl = workflow$contrasts_tbl, designMatrix = workflow$design_matrix,
                    showNotificationFn = function(...) NULL)
                stopifnot(isTRUE(values$files_copied))
                stopifnot(completeProtSummarySessionStateExport(directories,
                    experimentLabel = "dia_baseline_smoke", description = "DIA IQ checkpoint replay",
                    workflowArgsSaved = values$workflow_args_saved, filesCopied = values$files_copied,
                    reportGenerated = FALSE, reportPath = NULL, workflowData = workflow,
                    showNotificationFn = function(...) NULL))
            })
            receipt$ownership <<- diaStateOwnership(manager)
        }
    }, { session$userData$execute_replay() }), error = function(error) {
        receipt$status <<- "failed"
        receipt$error <<- conditionMessage(error)
        flush()
        stop(error)
    })
    finishDIARawReplayReceipt(receipt, root, source_store, source_head)
}

#' Finalize a replay using its pre-close ownership and session-end receipt.
#' @param receipt Bounded phase and ownership record captured before session close.
#' @param root Replay output directory.
#' @param source_store,source_head Original store and its unchanged head digest.
#' @return Completed replay receipt; no access to a closed store.
#' @noRd
finishDIARawReplayReceipt <- function(receipt, root, source_store, source_head) {
    path <- file.path(root, "checkpoint-replay.json")
    tryCatch({
        stopifnot(identical(digest::digest(file = file.path(source_store, "dia-state-head.rds"),
            algo = "sha256"), source_head))
        ended <- jsonlite::read_json(file.path(root, "gui-state-latest.json"))
        stopifnot(identical(ended$stage, "session_end"), ended$store$open_resources == 0L,
            ended$working_records == 0L, !is.null(receipt$ownership),
            receipt$ownership$full_s4_reads == 0L)
        receipt$status <- "complete_pending_output_verification"
        receipt$active_phase <- NULL
    }, error = function(error) {
        receipt$status <<- "failed"
        receipt$error <<- conditionMessage(error)
        jsonlite::write_json(receipt, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
        stop(error)
    })
    jsonlite::write_json(receipt, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
    invisible(receipt)
}


test_that("raw-checkpoint replay uses production scientific and publication actions", {
    mode <- Sys.getenv("MULTISCHOLAR_DIA025_MECHANISM_MEMBER")
    skip_if(!mode %in% c("gui-raw-replay:small", "gui-raw-replay:300"), "explicit replay only")
    root <- dirname(Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT"))
    source <- file.path(dirname(root), "browser-smoke.FKZX69")
    if (mode == "gui-raw-replay:300") {
        spec <- jsonlite::read_json(test_path("..", "testdata", "workflow-checkpoints",
            "gui-ruv-pearson-diagnostic.json"))$raw_replay_300
        stopifnot(identical(digest::digest(file = system.file("R", "MultiScholaR.rdb",
            package = "MultiScholaR"), algo = "sha256"), spec$runtime_rdb_sha256))
        source <- file.path(dirname(root), paste0("browser-smoke.", spec$source_member))
        state <- jsonlite::read_json(file.path(source, "gui-state-latest.json"))
        stopifnot(identical(digest::digest(file = file.path(state$store_directory,
            "dia-state-head.rds"), algo = "sha256"), spec$source_head_sha256),
            state$state_sizes$raw_data_s4$rows == spec$rows)
    }
    result <- runDIARawCheckpointReplay(source, root)
    expect_identical(result$status, "complete_pending_output_verification")
    expect_identical(result$ownership$full_s4_reads, 0L)
})

test_that("retained replay outputs and teardown are checked without repeating science", {
    skip_if(Sys.getenv("MULTISCHOLAR_DIA025_MECHANISM_MEMBER") != "gui-raw-replay:outputs",
        "explicit retained-output check only")
    root <- dirname(Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT"))
    source <- file.path(dirname(root), "browser-smoke.Vyr2MC")
    ended <- jsonlite::read_json(file.path(source, "gui-state-latest.json"))
    original_head <- file.path(ended$store_directory, "dia-state-head.rds")
    head_hash <- digest::digest(file = original_head, algo = "sha256")
    owned <- file.path(root, "replay-store")
    dir.create(owned)
    manager <- copyDIAReplayCheckpoint(ended$store_directory, owned)
    paths <- list(publication_graphs_dir = file.path(root, "graphs"))
    dir.create(paths$publication_graphs_dir)
    manager$.context$reference_states <- TRUE
    receipt <- list(status = "running", phases = list(), source_member = basename(source))
    ownership <- new.env(parent = emptyenv())
    shiny::testServer(function(input, output, session) {
        workflow <- shiny::reactiveValues(state_manager = manager)
        attachDIAStateSession(manager, session, workflow, receipt_dir = root)
        session$userData$check_outputs <- function() withDIAStateScope(manager, function() {
            manifest <- resolveDIAStateValue(manager$.context, manager$.context$snapshot$da_results)
            before <- getDIAStateReference(manager)
            files <- testthat::with_mocked_bindings(
                writeDIAInteractiveVolcanoOutputs(manifest, paths$publication_graphs_dir, NULL, 0.05),
                differentialAbundanceAnalysis = function(...) stop("output retry refitted DA"),
                materializeDIAStateRecord = function(...) stop("output retry hydrated S4"))
            expect_identical(getDIAStateReference(manager), before)
            expect_true(length(files) > 0L && all(file.size(files) > 0))
            html <- files[endsWith(files, ".html")]
            expect_length(html, length(manifest$results))
            expect_identical(files, html)
            expect_false(any(vapply(html, function(path)
                any(grepl('(src|href)="[^"]*_files/', readLines(path, warn = FALSE))), logical(1))))
            expect_true(all(vapply(html, function(path)
                any(grepl("glimma", readLines(path, warn = FALSE), ignore.case = TRUE)), logical(1))))
        })
    }, {
        session$userData$check_outputs()
        ownership$value <- diaStateOwnership(manager)
    })
    receipt$ownership <- ownership$value
    expect_error(workflowEntityCacheInfo(manager$.context$store), class = "saved_entity_scope")
    completed <- finishDIARawReplayReceipt(receipt, root, ended$store_directory, head_hash)
    expect_identical(completed$status, "complete_pending_output_verification")
    expect_error(finishDIARawReplayReceipt(receipt, root, ended$store_directory, "wrong-head"))
    expect_identical(jsonlite::read_json(file.path(root, "checkpoint-replay.json"))$status, "failed")
    finishDIARawReplayReceipt(receipt, root, ended$store_directory, head_hash)
    expect_identical(digest::digest(file = original_head, algo = "sha256"), head_hash)
})

#' Recover only presentation outputs, then independently verify retained science.
#' @param source Original failed small replay, never modified.
#' @param root Fresh guarded output directory.
#' @return Verification receipt; not a completed GUI/timing claim.
#' @noRd
verifyDIARetainedReplayOutputs <- function(source, root) {
    runner <- new.env(parent = globalenv())
    sys.source(test_path("..", "..", "tools", "profiling", "run_workflow_checkpoint_proof.R"), runner)
    project_source <- file.path(source, "dia_baseline_smoke")
    project <- file.path(root, "dia_baseline_smoke")
    dir.create(project)
    directories <- setupDirectories(project, "proteomics")
    paths <- directories$proteomics
    directories$proteomics_dia_baseline_smoke <- paths
    relative <- list.files(project_source, recursive = TRUE, all.files = FALSE)
    original <- file.path(project_source, relative)
    original_hash <- vapply(original, function(path) digest::digest(file = path, algo = "sha256"), character(1))
    for (name in relative) {
        destination <- file.path(project, name)
        dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
        stopifnot(file.copy(file.path(project_source, name), destination, overwrite = TRUE))
    }
    time_dirs <- list.dirs(file.path(project_source, "results", "proteomics", "publication_graphs",
        "filtering_qc"), recursive = FALSE, full.names = FALSE)
    stopifnot(length(time_dirs) == 1L)
    paths$time_dir <- file.path(paths$publication_graphs_dir, "filtering_qc", time_dirs)
    directories$proteomics <- directories$proteomics_dia_baseline_smoke <- paths
    archive <- readDIANativeLatestPath(file.path(project_source, "scripts", "proteomics", "filtered_session_latest.json"))
    restore <- file.path(root, "native-restore")
    dir.create(restore)
    imported <- openDIANativeSession(archive, restore)
    manager <- imported$manager
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    manager$.context$experiment_paths <- paths
    view <- withDIAStateScope(manager, function() resolveDIANativeSessionView(manager,
        manager$.context$native_metadata_reference, imported$extra_references))
    withr::local_options(multischolar.dia_native_sessions = TRUE, multischolar.test_mode = FALSE)
    rlang::local_bindings(config_list = view$config_list, project_dirs = directories,
        filtering_progress = methods::new("FilteringProgress"), .env = .GlobalEnv)
    result <- new.env(parent = emptyenv())
    grDevices::pdf(NULL)
    device <- grDevices::dev.cur()
    on.exit(if (device %in% grDevices::dev.list()) grDevices::dev.off(device), add = TRUE)
    shiny::testServer(function(input, output, session) {
        workflow <- do.call(shiny::reactiveValues, c(view, list(state_manager = manager)))
        workflow$ruv_normalised_for_da_analysis_obj <- getDIAStateReference(manager)
        attachDIAStateSession(manager, session, workflow, receipt_dir = root)
        session$userData$recover_outputs <- function() withDIAStateScope(manager, function() {
            before <- getDIAStateReference(manager)
            manifest <- resolveDIAStateValue(manager$.context, manager$.context$snapshot$da_results)
            testthat::with_mocked_bindings({
                writeDIAInteractiveVolcanoOutputs(manifest, paths$publication_graphs_dir,
                    view$uniprot_dat_cln, 0.05)
                states <- c("qvalue_filtered", "protein_peptide_filtered", "precursor_rollup",
                    "intensity_filtered", "sample_filtered", "replicate_filtered", "imputed",
                    "protein_s4_created", "protein_accession_cleaned", "protein_intensity_filtered",
                    "duplicates_removed", "protein_replicate_filtered", "ruv_corrected", "correlation_filtered")
                labels <- c(paste0(2:13, "_", states[1:12]), "11_RUV_filtered", "12_correlation_filtered")
                labels[1L] <- "2_qval_filtered"
                for (index in seq_along(states)) {
                    reference <- getDIAStateReference(manager, states[[index]])
                    spec <- readDIAStateMetadata(manager, reference)$data
                    if (!is.null(spec$components$peptide_data)) {
                        scope <- newWorkflowTableReadScope(manager$.context$store, spec$components$peptide_data)
                        counts <- tryCatch(summarizeSavedDIAFiltering(scope,
                            readSavedWorkflowTablePrototype(scope), spec$slots$protein_id_column)$summary,
                            finally = releaseWorkflowTableReadScope(scope))
                    } else counts <- readDIAProteinProgress(manager, reference)
                    updateProteinFiltering(counts, labels[[index]], "proteomics", "dia_baseline_smoke",
                        overwrite = TRUE, return_grid = TRUE, project_dirs = directories)
                }
                workflow$da_analysis_results_list <- manifest
                norm <- shiny::reactiveValues(correlation_filtered_obj = getDIAStateReference(manager),
                    best_k = view$ruv_k, correlation_threshold = view$correlation_threshold)
                runProtNormExportSessionWorkflow(workflow, norm,
                    list(norm_method = view$normalization_method, ruv_mode = view$ruv_mode), paths$source_dir,
                    withProgressFn = function(message, value, expr) force(expr), incProgressFn = function(...) NULL)
                values <- shiny::reactiveValues()
                completeProtSummaryWorkflowArgsSave(output, values, directories, workflow,
                    experimentLabel = "dia_baseline_smoke", description = "DIA IQ retained-output recovery",
                    showNotificationFn = function(...) NULL)
                stopifnot(isTRUE(values$workflow_args_saved))
                runProtSummaryPublicationCopy(output, values, directories, experimentLabel = "dia_baseline_smoke",
                    description = "DIA IQ retained-output recovery", contrastsTbl = workflow$contrasts_tbl,
                    designMatrix = workflow$design_matrix, showNotificationFn = function(...) NULL)
                stopifnot(isTRUE(values$files_copied))
                stopifnot(completeProtSummarySessionStateExport(directories,
                    experimentLabel = "dia_baseline_smoke", description = "DIA IQ retained-output recovery",
                    workflowArgsSaved = TRUE, filesCopied = TRUE, reportGenerated = FALSE,
                    reportPath = NULL, workflowData = workflow, showNotificationFn = function(...) NULL))
            }, differentialAbundanceAnalysis = function(...) stop("output recovery refitted DA"),
                materializeDIAStateRecord = function(...) stop("output recovery hydrated a full S4"))
            stopifnot(identical(getDIAStateReference(manager), before))
        })
    }, {
        session$userData$recover_outputs()
        result$ownership <- diaStateOwnership(manager)
    })
    # IQ's original QC plot uses its preprocessing input; no MaxLFQ fit is repeated.
    input <- file.path(paths$peptide_qc_dir, "peptide_values_imputed.tsv")
    iq_data <- iq::fast_read(input, sample_id = "Run", primary_id = "Protein.Group",
        secondary_id = "Stripped.Sequence", intensity_col = "Peptide.Imputed", annotation_col = NULL,
        filter_string_equal = NULL, filter_double_less = c("Q.Value" = "0.01", "PG.Q.Value" = "0.01"),
        intensity_col_sep = ";", na_string = "0")
    iq::fast_preprocess(iq_data$quant_table, median_normalization = FALSE,
        pdf_out = file.path(root, "qc-plots.pdf"))
    iq_data <- NULL
    presence <- runner$checkCheckpointIqQcOutputs(root, TRUE, TRUE, TRUE)
    extra_steps <- c("9_protein_s4_created", "10_protein_accession_cleaned",
        "11_protein_intensity_filtered", "12_duplicates_removed", "13_protein_replicate_filtered",
        "11_RUV_filtered", "12_correlation_filtered")
    extra_names <- as.vector(outer(paste0(extra_steps, "_combined_plots"), c("png", "pdf"), paste, sep = "."))
    stopifnot(all(file.exists(file.path(paths$time_dir, extra_names))),
        all(file.size(file.path(paths$time_dir, extra_names)) > 0))
    summary <- runner$checkCheckpointSummaryFiles(root, "dia_baseline_smoke", TRUE, TRUE)
    verify_root <- file.path(root, "verify-restore")
    dir.create(verify_root)
    verified <- openDIANativeSession(readDIANativeLatestPath(file.path(paths$source_dir,
        "filtered_session_latest.json")), verify_root)
    on.exit(closeDIAStateAdapter(verified$manager), add = TRUE)
    session_view <- withDIAStateScope(verified$manager, function() resolveDIANativeSessionView(verified$manager,
        verified$manager$.context$native_metadata_reference, verified$extra_references))
    # Explicit compatibility verification, one required checkpoint at a time.
    final <- verified$manager$getState()
    session_view$current_s4_object <- session_view$correlation_filtered_s4 <- final
    types <- readr::cols(.default = readr::col_double())
    types$cols[[final@protein_id_column]] <- readr::col_character()
    table <- readr::read_tsv(file.path(paths$protein_qc_dir, "ruv_normalised_results_cln_with_replicates.tsv"),
        col_types = types, progress = FALSE)
    science <- runner$checkCheckpointSessionContents(session_view, final, table, TRUE, verified$manager)
    da <- runner$verifyCheckpointDaExports(session_view, project)
    stopifnot(identical(vapply(original, function(path) digest::digest(file = path, algo = "sha256"),
        character(1)), original_hash))
    list(status = "passed_retained_output_recovery_and_verification", source_member = basename(source),
        source_unchanged = TRUE, ownership = result$ownership, iq_qc = presence, summary = summary,
        science = science, da = da, full_journey_timing_claim = FALSE)
}

test_that("retained replay verifies required outputs and science without rerunning the journey", {
    skip_if(Sys.getenv("MULTISCHOLAR_DIA025_MECHANISM_MEMBER") != "gui-raw-replay:verify",
        "explicit retained recovery and verification only")
    root <- dirname(Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT"))
    receipt <- list(status = "running")
    path <- file.path(root, "retained-output-verification.json")
    on.exit(jsonlite::write_json(receipt, path, auto_unbox = TRUE, pretty = TRUE, null = "null"), add = TRUE)
    receipt <- tryCatch(verifyDIARetainedReplayOutputs(file.path(dirname(root), "browser-smoke.Vyr2MC"), root),
        error = function(error) {
            receipt$status <<- "failed"
            receipt$error <<- conditionMessage(error)
            stop(error)
        })
    expect_identical(receipt$status, "passed_retained_output_recovery_and_verification")
})

test_that("real cohort attribution separates decoding from destination-buffer copies", {
    mode <- Sys.getenv("MULTISCHOLAR_DIA025_MECHANISM_MEMBER")
    threaded <- grepl("^gui-cohort-threads:(1|2|4|8)$", mode)
    skip_if(mode != "gui-cohort-profile" && !threaded, "explicit bounded cohort attribution only")
    if (threaded) {
        threads <- as.integer(sub(".*:", "", mode))
        previous_cpu <- arrow::cpu_count()
        previous_io <- arrow::io_thread_count()
        withr::defer({ arrow::set_cpu_count(previous_cpu); arrow::set_io_thread_count(previous_io) })
        arrow::set_cpu_count(threads)
        arrow::set_io_thread_count(1L)
        cursor <- openWorkflowArrowCursor
        changed <- 0L
        for (i in seq_along(body(cursor))) {
            expression <- body(cursor)[[i]]
            if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
                identical(expression[[2L]], as.name("scanner")) &&
                "use_threads" %in% names(expression[[3L]])) {
                stopifnot(identical(expression[[3L]][["use_threads"]], FALSE))
                expression[[3L]][["use_threads"]] <- threads > 1L
                body(cursor)[[i]] <- expression
                changed <- changed + 1L
            }
        }
        stopifnot(changed == 1L)
        testthat::local_mocked_bindings(openWorkflowArrowCursor = cursor)
    }
    root <- dirname(Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT"))
    source <- file.path(dirname(root), "browser-smoke.BpSyxB")
    state <- jsonlite::read_json(file.path(source, "gui-state-latest.json"))
    head <- file.path(state$store_directory, "dia-state-head.rds")
    head_hash <- digest::digest(file = head, algo = "sha256")
    expect_identical(head_hash, "5d64ca562669d94eea1677a1a03d5b410a8a81d8ea59add8cabd6065deb91353")
    owned <- file.path(root, "cohort-store")
    dir.create(owned)
    manager <- copyDIAReplayCheckpoint(state$store_directory, owned)
    withr::defer(closeDIAStateAdapter(manager))
    spec <- readDIAStateMetadata(manager, getDIAStateReference(manager, "replicate_filtered"))$data
    scope <- newWorkflowTableReadScope(manager$.context$store, spec$components$peptide_data, reader = "arrow")
    withr::defer(releaseWorkflowTableReadScope(scope))
    columns <- intersect(scope$metadata$columns, unique(c(unlist(spec$slots[c("sample_id",
        "protein_id_column", "peptide_sequence_column", "raw_quantity_column", "norm_quantity_column")]),
        "Modified.Sequence", "identification_peptide_count", "identification_peptidoform_count")))
    index <- saveWorkflowGroupIndex(scope, spec$slots$protein_id_column)
    groups <- readWorkflowIndexDirectory(scope$store, index$value$groups)
    take <- sum(cumsum(as.double(groups$size)) * workflowProjectedRowBytes(scope$metadata, columns) <= 268435456)
    groups <- groups[seq_len(take), , drop = FALSE]
    stopifnot(sum(groups$size) == 4473600L, nrow(groups) == 3728L)
    receipt <- list(status = "running", source_member = "BpSyxB", source_head_sha256 = head_hash,
        source = scope$reference$entity_id, source_rows = scope$metadata$rows,
        rows = sum(groups$size), complete_groups = nrow(groups), columns = columns, streams = list(),
        mode = mode, cpu_threads = arrow::cpu_count(), io_threads = arrow::io_thread_count(),
        tracing = !threaded, repetitions = list())
    path <- file.path(root, "cohort-attribution.json")
    flush <- function() jsonlite::write_json(receipt, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
    on.exit(flush(), add = TRUE)
    flush()
    stream <- streamSavedWorkflowTable
    allocate <- newWorkflowTableBuffer
    tracked <- character()
    testthat::local_mocked_bindings(newWorkflowTableBuffer = function(metadata, columns, rows) {
        buffer <- allocate(metadata, columns, rows)
        quantity <- spec$slots$norm_quantity_column
        if (!threaded) tracked <<- c(quantity = tracemem(buffer$columns[[quantity]]), positions = tracemem(buffer$positions))
        buffer
    }, streamSavedWorkflowTable = function(scope, columns, consume, ...) {
        callback <- 0
        batches <- 0L
        started <- proc.time()[["elapsed"]]
        value <- stream(scope, columns, function(data, offset, row_names) {
            at <- proc.time()[["elapsed"]]
            consume(data, offset, row_names)
            callback <<- callback + proc.time()[["elapsed"]] - at
            batches <<- batches + 1L
        }, ...)
        elapsed <- proc.time()[["elapsed"]] - started
        receipt$streams[[length(receipt$streams) + 1L]] <<- list(columns = columns,
            total_seconds = elapsed, callback_seconds = callback, decode_restore_seconds = elapsed - callback,
            batches = batches)
        flush()
        value
    })
    for (iteration in seq_len(if (threaded) 3L else 1L)) {
        started <- proc.time()
        traces <- tryCatch(if (threaded) {
            value <- readWorkflowGroupCohort(scope, columns, index$value$positions, groups)
            character()
        } else utils::capture.output(invisible(value <- readWorkflowGroupCohort(scope, columns,
            index$value$positions, groups))), error = function(error) {
                receipt$status <<- "failed"
                receipt$error <<- conditionMessage(error)
                stop(error)
            })
        stopifnot(length(traces) <= 2000L, sum(nchar(traces, type = "bytes")) <= 1048576)
        current <- gsub("[<>]", "", tracked)
        copies <- stats::setNames(integer(length(current)), names(current))
        for (line in traces[startsWith(traces, "tracemem[")]) {
            addresses <- regmatches(line, gregexpr("0x[0-9a-f]+", line))[[1L]]
            if (length(addresses) < 2L) next
            matched <- which(current == addresses[[1L]])
            if (length(matched)) {
                copies[matched] <- copies[matched] + 1L
                current[matched] <- addresses[[2L]]
            }
        }
        receipt$elapsed_seconds <- unname((proc.time() - started)[["elapsed"]])
        receipt$traced_copies <- as.list(copies)
        receipt$traced_copy_bytes_lower_bound <- sum(copies * c(8, 4) * nrow(value))
        receipt$result_sha256 <- digest::digest(value, algo = "sha256")
        receipt$result_bytes <- as.numeric(object.size(value))
        if (!identical(receipt$result_sha256,
            "f4e9951c1ab59eb9e541aca2c33c3302ea2c2aea7730de9a51c6e1b74868d92a")) {
            receipt$status <- "failed"
            receipt$error <- "cohort fingerprint differs from frozen result"
            stop(receipt$error)
        }
        expect_equal(nrow(value), sum(groups$size))
        expect_identical(names(value), columns)
        receipt$repetitions[[iteration]] <- list(elapsed = receipt$elapsed_seconds,
            sha256 = receipt$result_sha256)
        value <- NULL
        flush()
    }
    releaseWorkflowTableReadScope(scope)
    expect_identical(workflowEntityCacheInfo(manager$.context$store)$open_resources, 0L)
    expect_identical(digest::digest(file = head, algo = "sha256"), head_hash)
    receipt$status <- "complete"
    receipt$median_cohort_seconds <- stats::median(vapply(receipt$repetitions, `[[`, numeric(1), "elapsed"))
    receipt$performance_scope <- "one fixed real cohort; no full-workflow or statistical qualification claim"
    flush()
})

test_that("a 768 MiB cohort reads all retained peptide groups with bounded batches", {
    mode <- Sys.getenv("MULTISCHOLAR_DIA025_MECHANISM_MEMBER")
    skip_if(!mode %in% c("gui-cohort-size:256", "gui-cohort-size:768", "gui-cohort-size:768-threads8"), "explicit cohort-size experiment only")
    budget <- as.numeric(sub("-threads8$", "", sub(".*:", "", mode))) * 1024^2
    threads <- if (endsWith(mode, "-threads8")) 8L else 2L
    root <- dirname(Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT"))
    source <- jsonlite::read_json(file.path(dirname(root), "browser-smoke.BpSyxB", "gui-state-latest.json"))
    head <- file.path(source$store_directory, "dia-state-head.rds")
    head_hash <- digest::digest(file = head, algo = "sha256")
    stopifnot(identical(head_hash, "5d64ca562669d94eea1677a1a03d5b410a8a81d8ea59add8cabd6065deb91353"))
    previous_cpu <- arrow::cpu_count()
    previous_io <- arrow::io_thread_count()
    withr::defer({ arrow::set_cpu_count(previous_cpu); arrow::set_io_thread_count(previous_io) })
    arrow::set_cpu_count(threads)
    arrow::set_io_thread_count(1L)
    cursor <- openWorkflowArrowCursor
    changed <- 0L
    for (i in seq_along(body(cursor))) {
        expression <- body(cursor)[[i]]
        if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
            identical(expression[[2L]], as.name("scanner")) && "use_threads" %in% names(expression[[3L]])) {
            stopifnot(identical(expression[[3L]][["use_threads"]], FALSE))
            expression[[3L]][["use_threads"]] <- TRUE
            body(cursor)[[i]] <- expression
            changed <- changed + 1L
        }
    }
    stopifnot(changed == 1L)
    reader <- readWorkflowGroupCohort
    code <- deparse(body(reader), width.cutoff = 500L)
    stopifnot(sum(lengths(regmatches(code, gregexpr("268435456", code, fixed = TRUE)))) == 3L)
    body(reader) <- parse(text = gsub("268435456", format(budget, scientific = FALSE), code, fixed = TRUE))[[1L]]
    testthat::local_mocked_bindings(openWorkflowArrowCursor = cursor, readWorkflowGroupCohort = reader)
    owned <- file.path(root, "cohort-store")
    dir.create(owned)
    manager <- copyDIAReplayCheckpoint(source$store_directory, owned)
    withr::defer(closeDIAStateAdapter(manager))
    spec <- readDIAStateMetadata(manager, getDIAStateReference(manager, "replicate_filtered"))$data
    scope <- newWorkflowTableReadScope(manager$.context$store, spec$components$peptide_data, reader = "arrow")
    withr::defer(releaseWorkflowTableReadScope(scope))
    columns <- intersect(scope$metadata$columns, unique(c(unlist(spec$slots[c("sample_id",
        "protein_id_column", "peptide_sequence_column", "raw_quantity_column", "norm_quantity_column")]),
        "Modified.Sequence", "identification_peptide_count", "identification_peptidoform_count")))
    index <- saveWorkflowGroupIndex(scope, spec$slots$protein_id_column)
    groups <- readWorkflowIndexDirectory(scope$store, index$value$groups)
    row_bytes <- workflowProjectedRowBytes(scope$metadata, columns)
    stopifnot(nrow(groups) == 10000L, all(groups$size == 1200L), sum(groups$size) == 12000000L,
        max(groups$size) * row_bytes <= 67108864)
    receipt <- list(status = "running", mode = mode, budget_bytes = budget, cpu_threads = threads,
        io_threads = 1L, source_head_sha256 = head_hash, source_entity = scope$reference$entity_id,
        runtime_rdb_sha256 = digest::digest(file = system.file("R", "MultiScholaR.rdb", package = "MultiScholaR"), algo = "sha256"),
        cohorts = list(), read_seconds = 0, block_hashes = character())
    path <- file.path(root, "cohort-size.json")
    flush <- function() jsonlite::write_json(receipt, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
    on.exit(flush(), add = TRUE)
    flush()
    tryCatch({
        next_group <- 1L
        while (next_group <= nrow(groups)) {
            remaining <- seq.int(next_group, nrow(groups))
            take <- sum(cumsum(as.double(groups$size[remaining])) * row_bytes <= budget)
            stopifnot(take > 0L)
            chosen <- groups[remaining[seq_len(take)], , drop = FALSE]
            started <- proc.time()[["elapsed"]]
            value <- readWorkflowGroupCohort(scope, columns, index$value$positions, chosen)
            elapsed <- proc.time()[["elapsed"]] - started
            receipt$read_seconds <- receipt$read_seconds + elapsed
            receipt$cohorts[[length(receipt$cohorts) + 1L]] <- list(groups = nrow(chosen),
                rows = nrow(value), result_bytes = as.numeric(object.size(value)), read_seconds = elapsed)
            receipt$block_hashes <- c(receipt$block_hashes,
                fingerprintDIACohortBlocks(value, isTRUE(scope$metadata$automatic_row_names)))
            value <- NULL
            next_group <- next_group + take
            flush()
        }
        releaseWorkflowTableReadScope(scope)
        stopifnot(length(receipt$block_hashes) == 1250L,
            workflowEntityCacheInfo(manager$.context$store)$open_resources == 0L,
            identical(digest::digest(file = head, algo = "sha256"), head_hash))
        if (budget == 805306368) {
            files <- Sys.glob(file.path(dirname(root), "browser-smoke.*", "cohort-size.json"))
            controls <- Filter(function(x) identical(x$status, "complete") && identical(x$mode, "gui-cohort-size:256") &&
                identical(x$runtime_rdb_sha256, receipt$runtime_rdb_sha256) && identical(x$source_head_sha256, head_hash),
                lapply(files, jsonlite::read_json, simplifyVector = TRUE))
            stopifnot(length(controls) == 1L, length(receipt$cohorts) == 1L,
                identical(unname(controls[[1L]]$block_hashes), unname(receipt$block_hashes)))
            receipt$control_read_seconds <- controls[[1L]]$read_seconds
            receipt$exact_block_parity <- TRUE
        } else stopifnot(length(receipt$cohorts) == 3L)
        receipt$status <- "complete"
        receipt$scope <- "full 12-million-row projected peptide table; assembly only, not scientific kernel capacity"
    }, error = function(error) {
        receipt$status <<- "failed"
        receipt$error <<- conditionMessage(error)
        stop(error)
    })
    flush()
})

#' Execute one complete-group imputation task using the unchanged public method.
#' @param part Bounded input with stable source ordinals.
#' @param prototype Zero-assay scientific object with the full design.
#' @param fail Explicit lifecycle failure injection.
#' @return Bounded scientific tables and worker CPU accounting.
#' @noRd
runDIAImputationWorkerBatch <- function(part, prototype, fail = FALSE) {
    if (fail) stop("injected imputation worker failure")
    started <- proc.time()
    methods::slot(prototype, "peptide_data") <- part
    result <- get("peptideMissingValueImputation", asNamespace("MultiScholaR"))(
        prototype, proportion_missing_values = 0.5, core_utilisation = NA)
    args <- methods::slot(result, "args")$peptideMissingValueImputation
    scientific <- list(data = methods::slot(result, "peptide_data"),
        support_table = args$support_table, summary = args$imputation_summary)
    bytes <- as.numeric(utils::object.size(scientific))
    stopifnot(bytes <= 67108864)
    list(scientific = scientific, bytes = bytes, pid = Sys.getpid(),
        cpu = unname(sum((proc.time() - started)[c("user.self", "sys.self")])))
}

#' Run bounded waves on the existing DIA workers and stream ordered outputs.
#' @param cohort,prototype Real selected data and full scientific design.
#' @param workers Existing backend worker count.
#' @param root Fresh output directory.
#' @param observe Bounded progress callback.
#' @param inject_failure Whether the first worker should fail.
#' @return Exact job fingerprints and aggregate worker accounting.
#' @noRd
runDIAImputationWorkerWaves <- function(cohort, prototype, workers, root, observe,
    inject_failure = FALSE) {
    withr::local_envvar(c(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1"))
    cluster <- startDIAWorkerCluster(newDIAWorkerSpec(workers))
    writers <- new.env(parent = emptyenv())
    on.exit({
        stopDIAWorkerCluster(cluster)
        for (process in cluster) process$wait(timeout = 5000)
        stopifnot(!any(vapply(cluster, function(x) x$is_alive(), logical(1))))
        for (writer in as.list(writers)) abortWorkflowTableWriter(writer)
    }, add = TRUE)
    multidplyr::cluster_send(cluster, {
        .libPaths(!!.libPaths())
        loadNamespace("MultiScholaR")
        arrow::set_cpu_count(1L)
        arrow::set_io_thread_count(1L)
        NULL
    })
    worker <- runDIAImputationWorkerBatch
    environment(worker) <- baseenv()
    positions <- attr(cohort, "workflow_source_order")
    starts <- seq.int(1L, nrow(cohort), by = 9600L)
    fingerprints <- character(length(starts))
    pids <- integer()
    cpu <- numeric()
    max_input <- max_output <- 0
    for (wave in seq.int(1L, length(starts), by = workers)) {
        jobs <- seq.int(wave, min(length(starts), wave + workers - 1L))
        for (slot in seq_along(jobs)) {
            job <- jobs[[slot]]
            at <- seq.int(starts[[job]], length.out = min(9600L, nrow(cohort) - starts[[job]] + 1L))
            part <- cohort[at, , drop = FALSE]
            attr(part, "workflow_source_order") <- NULL
            attr(part, "row.names") <- .set_row_names(nrow(part))
            part$.dia_source_order <- positions[at]
            stopifnot(nrow(part) <= 10000L, as.numeric(object.size(part)) <= 4194304)
            max_input <- max(max_input, as.numeric(object.size(part)))
            cluster[[slot]]$call(worker, list(part = part, prototype = prototype,
                fail = inject_failure && job == 1L))
            part <- NULL
        }
        for (slot in seq_along(jobs)) {
            cluster[[slot]]$poll_process(-1)
            response <- cluster[[slot]]$read()
            if (!is.null(response$error)) stop(response$error)
            result <- response$result
            key <- as.character(result$pid)
            pids <- union(pids, result$pid)
            cpu[key] <- (if (key %in% names(cpu)) cpu[[key]] else 0) + result$cpu
            max_output <- max(max_output, result$bytes)
            fingerprints[[jobs[[slot]]]] <- digest::digest(result$scientific, algo = "sha256")
            for (field in c("data", "support_table")) {
                table <- result$scientific[[field]]
                if (!is.data.frame(table) || !ncol(table)) next
                if (!exists(field, writers, inherits = FALSE)) {
                    metadata <- workflowSavedTableMetadata(table, ".rownames", 3L)
                    metadata$rows <- NULL
                    writers[[field]] <- openWorkflowTableWriter(file.path(root, paste0(field, ".parquet")), metadata)
                }
                for (first in seq.int(1L, max(1L, nrow(table)), by = 10000L)) {
                    at <- if (nrow(table)) seq.int(first, min(nrow(table), first + 9999L)) else integer()
                    appendWorkflowTableWriter(writers[[field]], table, at)
                }
                table <- NULL
            }
            result <- response <- NULL
        }
        observe(max(jobs), length(starts))
    }
    outputs <- lapply(as.list(writers), finishWorkflowTableWriter)
    stopDIAWorkerCluster(cluster)
    for (process in cluster) process$wait(timeout = 5000)
    stopifnot(!any(vapply(cluster, function(x) x$is_alive(), logical(1))))
    list(job_hashes = fingerprints, worker_pids = pids, worker_cpu_seconds = as.list(cpu),
        max_input_bytes = max_input, max_output_bytes = max_output,
        output_rows = lapply(outputs, function(x) x$rows), workers_stopped = TRUE)
}

test_that("real imputation executes on four then eight bounded scientific workers", {
    mode <- Sys.getenv("MULTISCHOLAR_DIA025_MECHANISM_MEMBER")
    skip_if(!mode %in% c("gui-impute-workers:preflight", "gui-impute-workers:4", "gui-impute-workers:8"),
        "explicit scientific-worker experiment only")
    preflight <- endsWith(mode, "preflight")
    workers <- if (preflight) 4L else as.integer(sub(".*:", "", mode))
    root <- dirname(Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT"))
    source <- jsonlite::read_json(file.path(dirname(root), "browser-smoke.BpSyxB", "gui-state-latest.json"))
    head <- file.path(source$store_directory, "dia-state-head.rds")
    head_hash <- digest::digest(file = head, algo = "sha256")
    stopifnot(identical(head_hash, "5d64ca562669d94eea1677a1a03d5b410a8a81d8ea59add8cabd6065deb91353"))
    owned <- file.path(root, "worker-store")
    dir.create(owned)
    manager <- copyDIAReplayCheckpoint(source$store_directory, owned)
    withr::defer(closeDIAStateAdapter(manager))
    reference <- getDIAStateReference(manager, "replicate_filtered")
    spec <- readDIAStateMetadata(manager, reference)$data
    prototype <- withDIAStateScope(manager, function() readDIAExecutionComponents(manager, reference, fields = "design_matrix"))
    stopifnot(nrow(prototype@peptide_data) == 0L, nrow(prototype@design_matrix) == 300L,
        as.numeric(object.size(prototype)) <= 8388608)
    scope <- newWorkflowTableReadScope(manager$.context$store, spec$components$peptide_data, reader = "arrow")
    withr::defer(releaseWorkflowTableReadScope(scope))
    columns <- intersect(scope$metadata$columns, unique(c(unlist(spec$slots[c("sample_id",
        "protein_id_column", "peptide_sequence_column", "raw_quantity_column", "norm_quantity_column")]),
        "Modified.Sequence", "identification_peptide_count", "identification_peptidoform_count")))
    index <- saveWorkflowGroupIndex(scope, spec$slots$protein_id_column)
    groups <- readWorkflowIndexDirectory(scope$store, index$value$groups)
    stopifnot(nrow(groups) == 10000L, all(groups$size == 1200L))
    if (preflight) groups <- groups[seq_len(32L), , drop = FALSE]
    reader <- readWorkflowGroupCohort
    code <- deparse(body(reader), width.cutoff = 500L)
    stopifnot(sum(lengths(regmatches(code, gregexpr("268435456", code, fixed = TRUE)))) == 3L)
    body(reader) <- parse(text = gsub("268435456", "805306368", code, fixed = TRUE))[[1L]]
    receipt <- list(status = "running", mode = mode, workers = workers, source_head_sha256 = head_hash,
        input_rows = sum(groups$size), full_design_rows = 300L, completed_jobs = 0L)
    path <- file.path(root, "scientific-workers.json")
    flush <- function() jsonlite::write_json(receipt, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
    on.exit(flush(), add = TRUE)
    flush()
    tryCatch({
        started <- proc.time()[["elapsed"]]
        cohort <- reader(scope, columns, index$value$positions, groups)
        receipt$read_seconds <- proc.time()[["elapsed"]] - started
        releaseWorkflowTableReadScope(scope)
        started <- proc.time()[["elapsed"]]
        result <- runDIAImputationWorkerWaves(cohort, prototype, workers, root, function(done, total) {
            receipt$completed_jobs <<- done
            receipt$total_jobs <<- total
            flush()
        })
        receipt$worker_start_compute_write_stop_seconds <- proc.time()[["elapsed"]] - started
        stopifnot(length(result$worker_pids) == workers, result$workers_stopped)
        receipt$result <- result
        # Independent serial checks occur after the measured worker phase.
        selected <- if (preflight) seq_along(result$job_hashes) else unique(c(1L,
            ceiling(length(result$job_hashes) / 2), length(result$job_hashes)))
        positions <- attr(cohort, "workflow_source_order")
        for (job in selected) {
            at <- seq.int((job - 1L) * 9600L + 1L, min(job * 9600L, nrow(cohort)))
            part <- cohort[at, , drop = FALSE]
            attr(part, "workflow_source_order") <- NULL
            attr(part, "row.names") <- .set_row_names(nrow(part))
            part$.dia_source_order <- positions[at]
            expected <- runDIAImputationWorkerBatch(part, prototype)$scientific
            stopifnot(identical(digest::digest(expected, algo = "sha256"), result$job_hashes[[job]]))
            part <- expected <- NULL
        }
        receipt$serial_oracle_jobs <- selected
        if (preflight) {
            failure_root <- file.path(root, "failure")
            dir.create(failure_root)
            expect_error(runDIAImputationWorkerWaves(cohort, prototype, workers, failure_root,
                function(...) NULL, inject_failure = TRUE), "injected imputation worker failure")
            receipt$failure_observed <- TRUE
        }
        if (workers == 8L) {
            controls <- Filter(function(x) identical(x$status, "complete") && identical(x$mode, "gui-impute-workers:4") &&
                identical(x$source_head_sha256, head_hash), lapply(Sys.glob(file.path(dirname(root),
                    "browser-smoke.*", "scientific-workers.json")), jsonlite::read_json, simplifyVector = TRUE))
            stopifnot(length(controls) == 1L,
                identical(unname(controls[[1L]]$result$job_hashes), unname(result$job_hashes)))
            receipt$all_jobs_match_four_workers <- TRUE
        }
        stopifnot(identical(digest::digest(file = head, algo = "sha256"), head_hash),
            workflowEntityCacheInfo(manager$.context$store)$open_resources == 0L)
        receipt$status <- "complete"
        receipt$scope <- "real imputation kernel, complete protein groups/full design, ordered streamed data and support outputs"
    }, error = function(error) {
        receipt$status <<- "failed"
        receipt$error <<- conditionMessage(error)
        stop(error)
    })
    flush()
})

test_that("bulk rollup experiment preserves complete operations and rejects relational rounding", {
    mode <- Sys.getenv("MULTISCHOLAR_DIA025_MECHANISM_MEMBER")
    skip_if(!grepl("^gui-bulk-rollup:(diagnostic|diagnostic-preflight|preflight|screen:([ABCD]|P([0-9]|1[0-4])):[12]|full:([ABCD]|P([0-9]|1[0-4])))$", mode),
        "explicit bounded bulk experiment only")
    root <- dirname(Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT"))
    proof <- dirname(root)
    native_preflight <- identical(mode, "gui-bulk-rollup:preflight")
    diagnostic <- mode %in% c("gui-bulk-rollup:diagnostic", "gui-bulk-rollup:diagnostic-preflight")
    preflight <- native_preflight || endsWith(mode, "diagnostic-preflight")
    withr::local_envvar(c(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
        DUCKPLYR_FALLBACK_COLLECT = "0", DUCKPLYR_FALLBACK_AUTOUPLOAD = "0"))
    receipt <- list(status = "running", mode = mode, phases = list(), failures = list(),
        profiling = diagnostic, promotion_authorized = FALSE)
    path <- file.path(root, "bulk-rollup.json")
    flush <- function() jsonlite::write_json(receipt, path, auto_unbox = TRUE, pretty = TRUE,
        null = "null", digits = NA)
    on.exit(flush(), add = TRUE)
    phase <- function(label, fn) {
        receipt$active_phase <<- label
        flush()
        started <- proc.time()
        value <- tryCatch(fn(), error = function(error) {
            receipt$failures[[length(receipt$failures) + 1L]] <<- list(phase = label,
                error = conditionMessage(error), elapsed = unname((proc.time() - started)[["elapsed"]]))
            flush()
            stop(error)
        })
        receipt$phases[[length(receipt$phases) + 1L]] <<- list(phase = label,
            elapsed = unname((proc.time() - started)[["elapsed"]]),
            cpu = unname(sum((proc.time() - started)[c("user.self", "sys.self")])))
        flush()
        value
    }
    if (diagnostic) {
        receipt$memory_boundaries <- list()
        bindings <- makeDIARollupDiagnostic(root, function(stage) {
            receipt$memory_boundaries[[length(receipt$memory_boundaries) + 1L]] <<-
                list(stage = stage, observation = snapshotWorkflowMemory())
            flush()
        })
        testthat::local_mocked_bindings(newDIARelationalScope = bindings$scope,
            streamDIARelationalQuery = bindings$stream)
    }
    tryCatch({
        frozen <- bulkRollupFrozenKernel(proof)
        receipt$frozen_source_sha256 <- digest::digest(file = file.path(proof,
            "browser-smoke.x3Dkg2", "source", "R", "func_prot_rollup_methods.R"), algo = "sha256")
        receipt$runtime_rdb_sha256 <- digest::digest(file = file.path(Sys.getenv("MULTISCHOLAR_DIA008_LIBRARY"),
            "MultiScholaR", "R", "MultiScholaR.rdb"), algo = "sha256")
        receipt$source_hashes <- vapply(c(
            here::here("tests/testthat/test-workflow-memory-mechanism-member.R"),
            here::here("tests/testthat/helper-workflow-bulk-rollup.R"),
            here::here("tests/testthat/helper-workflow-relational-rollup.R"),
            here::here("src/dia_rsum.c"),
            here::here("tools/profiling/run_workflow_checkpoint_proof.R"),
            here::here("tools/profiling/run_workflow_checkpoint_proof.sh"),
            here::here("tools/profiling/workflow_checkpoint_dia_driver.R")),
            \(path) digest::digest(file = path, algo = "sha256"), character(1))
        withr::local_libpaths(c(file.path(proof, "bulk-rollup-library"), .libPaths()))
        receipt$packages <- lapply(c("duckplyr", "duckdb", "dplyr", "arrow", "multidplyr", "nanoarrow", "fastduckdb"), function(package) {
            directory <- find.package(package)
            binaries <- list.files(directory, "[.](so|rdb|rdx)$", recursive = TRUE, full.names = TRUE)
            list(package = package, version = as.character(utils::packageVersion(package)),
                files = stats::setNames(vapply(binaries, \(path) digest::digest(file = path, algo = "sha256"),
                    character(1)), substring(binaries, nchar(directory) + 2L)))
        })
        receipt$duckplyr_source_sha256 <- digest::digest(file = file.path(proof,
            "bulk-rollup-library", "duckplyr_1.2.1.tar.gz"), algo = "sha256")
        receipt$exact_sum_extension_sha256 <- diaRollupNativeInfo()$manifest$binary_sha256
        if (native_preflight) {
            withr::local_libpaths(c(file.path(proof, "bulk-rollup-library"), .libPaths()))
            stopifnot(as.character(utils::packageVersion("duckdb")) == "1.5.5",
                as.character(utils::packageVersion("duckplyr")) == "1.2.1")
            dir.create(file.path(root, "numerical"))
            store <- newWorkflowEntityStore(file.path(root, "numerical"), "bulk-numeric", "DIA")
            on.exit(closeWorkflowEntityStore(store), add = TRUE)
            input <- data.frame(group = c(1L, 1L, 1L, 2L, 2L, 2L, 2L),
                ordinal = 1:7, quantity = c(1e16, 1, 1, 1e16, 1, 1e16, 2))
            ref <- saveWorkflowTable(store, input, "numeric", list(input = "exact", design = "two",
                parameters = "ordered-sum", selectors = "all", software = "frozen"), ".rownames")
            scope <- newWorkflowTableReadScope(store, ref)
            on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
            con <- scope$connection
            DBI::dbExecute(con, "SET threads = 8")
            DBI::dbExecute(con, "SET memory_limit = '12GiB'")
            query <- paste0("SELECT * FROM read_parquet(", DBI::dbQuoteString(con, scope$path), ")")
            direct <- DBI::dbGetQuery(con, paste0('SELECT "group", fsum(quantity ORDER BY ordinal) AS total FROM (',
                query, ') GROUP BY "group" ORDER BY "group"'))
            lazy <- duckplyr::read_sql_duckdb(query, con = scope$connection, prudence = "stingy")
            translated <- lazy |> dplyr::arrange(ordinal) |>
                dplyr::summarise(total = sum(quantity, na.rm = TRUE), .by = group) |> dplyr::arrange(group)
            translated <- as.data.frame(translated)
            native <- vapply(split(input$quantity, input$group), sum, numeric(1))
            receipt$relational <- list(native = unname(native), C = direct$total, D = translated$total,
                C_eligible = identical(direct$total, unname(native)),
                D_eligible = identical(translated$total, unname(native)),
                conclusion = "numerical mismatch disqualifies; no timed relational members")
            expect_false(receipt$relational$C_eligible)
            expect_false(receipt$relational$D_eligible)
            releaseWorkflowTableReadScope(scope)
            expect_identical(store$open_resources, 0L)
            expect_error(local({
                cancelled <- newWorkflowTableReadScope(store, ref, reader = "arrow")
                on.exit(releaseWorkflowTableReadScope(cancelled))
                streamSavedWorkflowTable(cancelled, "quantity", function(...) stop("injected cancellation"))
            }), "injected cancellation")
            expect_identical(store$open_resources, 0L)
            closeWorkflowEntityStore(store)
            table <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))@peptide_data
            invoke <- function(fn, data, logged = FALSE) fn(data, sample_id_column = "Run",
                protein_id_column = "Protein.Group", peptide_sequence_column = "Stripped.Sequence",
                modified_peptide_sequence_column = "Modified.Sequence", precursor_quantity_column = "Precursor.Quantity",
                precursor_normalised_column = "Precursor.Normalised", precursor_id_column = "Precursor.Id",
                core_utilisation = NA, is_logged_data = logged, return_rollup_result = TRUE)
            variants <- list(ordinary = table, empty = table[0, , drop = FALSE])
            missing <- table
            missing$Precursor.Quantity[1:2] <- NA_real_
            missing$Precursor.Normalised[3] <- NA_real_
            missing$Precursor.Quantity[4] <- 0
            missing$Run <- factor(missing$Run, levels = c(rev(unique(missing$Run)), "unused"))
            missing$Protein.Ids[1:3] <- c("P1; P2", " P2;P3;", "P1; P2")
            variants$missing <- missing
            threshold <- table[rep(1L, 4L), , drop = FALSE]
            threshold$Precursor.Id <- paste0("threshold-", 1:4)
            threshold$Precursor.Quantity <- threshold$Precursor.Normalised <- c(1e16, 1, 1e16, 2)
            variants$compensated_threshold <- threshold
            threshold <- threshold[1:3, , drop = FALSE]
            threshold$Precursor.Quantity <- threshold$Precursor.Normalised <- c(1e16, 1, 1)
            variants$threshold <- threshold
            for (data in variants) expect_identical(invoke(rollUpPrecursorToPeptideHelper, data), invoke(frozen, data))
            for (fn in list(frozen, rollUpPrecursorToPeptideHelper)) {
                for (bad in c(-1, Inf, -Inf, NaN)) {
                    invalid <- table
                    invalid$Precursor.Quantity[[1L]] <- bad
                    expect_error(invoke(fn, invalid), "negative or non-finite")
                }
                expect_error(invoke(fn, table, TRUE), "logged abundance")
                expect_error(invoke(fn, rbind(table, table[1L, ])), "duplicate precursor")
                invalid <- table
                invalid$identification_peptide_count <- 3L
                invalid$identification_peptide_count[[1L]] <- 4L
                expect_error(invoke(fn, invalid), "conflicting frozen")
            }
            cluster <- startDIAWorkerCluster(newDIAWorkerSpec(8L))
            expect_error(multidplyr::cluster_call(cluster, stop("injected bulk worker failure")),
                "injected bulk worker failure")
            stopDIAWorkerCluster(cluster)
            expect_true(all(vapply(cluster, \(worker) !worker$is_alive(), logical(1))))
        }
        arms <- if (diagnostic) "P0" else if (preflight) c("A", "B") else strsplit(mode, ":", fixed = TRUE)[[1L]][[3L]]
        receipt$results <- list()
        for (arm in arms) {
            owned <- file.path(root, paste0("arm-", arm))
            dir.create(owned)
            samples <- NULL
            if (preflight) {
                object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
                manager <- newDIAStateAdapter(newWorkflowEntityStore(owned, paste0("bulk-", arm), "DIA"),
                    WorkflowState$new(audit_enabled = FALSE))
                manager$.context$reference_states <- TRUE
                manager$saveState("protein_peptide_filtered", object, list(), "bulk fixture")
            } else {
                state <- jsonlite::read_json(file.path(proof, "browser-smoke.BpSyxB", "gui-state-latest.json"))
                head <- file.path(state$store_directory, "dia-state-head.rds")
                receipt$source_head_sha256 <- digest::digest(file = head, algo = "sha256")
                stopifnot(identical(receipt$source_head_sha256,
                    "5d64ca562669d94eea1677a1a03d5b410a8a81d8ea59add8cabd6065deb91353"))
                manager <- copyDIAReplayCheckpoint(state$store_directory, owned)
                selectDIAStateReference(manager, "protein_peptide_filtered")
                spec <- readDIAStateMetadata(manager)$data
                stopifnot(spec$dimensions$rows == 24000000L)
                receipt$source_reference <- spec$components$peptide_data
                if (grepl(":screen:", mode)) {
                    design <- withDIAStateScope(manager, function()
                        readDIAExecutionComponents(manager, fields = "design_matrix"))@design_matrix
                    stopifnot(nrow(design) == 300L)
                    samples <- as.character(design[[spec$slots$sample_id]][c(1:4, 149:152, 297:300)])
                }
            }
            on.exit(closeDIAStateAdapter(manager), add = TRUE)
            started <- proc.time()[["elapsed"]]
            value <- runDIABulkOperation(manager, owned, proof, arm, samples, phase, failure = preflight)
            phase("teardown", function() closeDIAStateAdapter(manager))
            elapsed <- proc.time()[["elapsed"]] - started
            receipt$results[[arm]] <- c(value, list(complete_operation_seconds = elapsed, samples = samples))
            # Independent bounded readback follows the complete-operation timer.
            identity <- readSavedStoreIdentity(file.path(owned, ".workflow-entities", "store.json"))
            check_store <- newWorkflowEntityStore(owned, identity$project_id, identity$workflow_id)
            check_scope <- newWorkflowTableReadScope(check_store, value$reference, reader = "arrow")
            hashes <- character()
            streamSavedWorkflowTable(check_scope, check_scope$metadata$columns, function(data, offset, row_names) {
                attr(data, "dia_source_order") <- NULL
                attr(data, "row.names") <- .set_row_names(nrow(data))
                hashes[[length(hashes) + 1L]] <<- digest::digest(data, algo = "sha256")
            })
            releaseWorkflowTableReadScope(check_scope)
            closeWorkflowEntityStore(check_store)
            receipt$results[[arm]]$output_hashes <- hashes
            if (!preflight) stopifnot(identical(digest::digest(file = head, algo = "sha256"), receipt$source_head_sha256))
            flush()
        }
        if (native_preflight) {
            expect_identical(receipt$results$A$output_hashes, receipt$results$B$output_hashes)
            expect_identical(receipt$results$A$summary, receipt$results$B$summary)
            expect_length(receipt$failures, 2L)
        }
        receipt$status <- "complete"
    }, error = function(error) {
        receipt$status <<- "failed"
        receipt$error <<- conditionMessage(error)
        flush()
        stop(error)
    })
    flush()
})
