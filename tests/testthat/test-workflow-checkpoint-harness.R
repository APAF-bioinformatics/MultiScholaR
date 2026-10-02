checkpointHarnessTestEnvironment <- function() {
    environment <- new.env(parent = baseenv())
    sys.source(testthat::test_path("..", "..", "tools", "profiling",
        "run_workflow_checkpoint_proof.R"), envir = environment)
    environment
}

test_that("installed identities require compiled bytes and retain pure R compatibility", {
    runner <- checkpointHarnessTestEnvironment()
    root <- withr::local_tempdir()
    package <- file.path(root, "MultiScholaR")
    for (folder in c("Meta", "R", "libs")) dir.create(file.path(package, folder), recursive = TRUE)
    for (path in c("NAMESPACE", "Meta/package.rds", "R/MultiScholaR.rdb", "R/MultiScholaR.rdx"))
        writeLines("identity fixture", file.path(package, path))
    description <- c("Package: MultiScholaR", "Version: 0.5.0", "Built: identity fixture")
    writeLines(description, file.path(package, "DESCRIPTION"))
    expect_length(runner$inspectCheckpointLibrary(root)$native_libraries, 0L)
    writeLines(c(description, "NeedsCompilation: yes"), file.path(package, "DESCRIPTION"))
    expect_error(runner$inspectCheckpointLibrary(root), "native library is missing")
    binary <- file.path(package, "libs", "MultiScholaR.so")
    writeBin(as.raw(1:3), binary)
    first <- runner$inspectCheckpointLibrary(root)$native_libraries
    expect_identical(first, workflowNativeLibraryIdentity(package))
    expect_identical(first, jsonlite::fromJSON(jsonlite::toJSON(first, auto_unbox = TRUE),
        simplifyVector = FALSE))
    writeBin(as.raw(3:1), binary)
    expect_false(identical(first, runner$inspectCheckpointLibrary(root)$native_libraries))
})

test_that("capacity scope keeps frozen controls and explicit installed activation separate", {
    root <- normalizePath(test_path("..", ".."))
    driver <- new.env(parent = baseenv())
    sys.source(file.path(root, "tools", "profiling", "workflow_checkpoint_dia_driver.R"), driver)
    spec <- driver$readCheckpointGuiCapacitySpec(root)
    expect_identical(unlist(spec$run_order), c("small:resident", "small:qs2_duckdb",
        "oom_reference:qs2_duckdb", "report_1gb:qs2_duckdb", "report_2gb:qs2_duckdb"))
    expect_false(spec$promotion_authorized)
    expect_identical(vapply(spec$recipes, \(x) x$precursor_rows, numeric(1)),
        c(24576, 7680000, 15360000))
    small <- makeWorkflowPrecursorChunk(spec$recipes[[1L]], seq_len(256L))
    expect_equal(nrow(small), 24576)
    expect_equal(length(unique(small$Run)), 12L)
    expect_gte(length(unique(small$Stripped.Sequence)), 500L)
    library_path <- withr::local_tempdir()
    dir.create(file.path(library_path, "MultiScholaR"))
    run <- withr::local_tempdir()
    native <- driver$checkpointInstalledApp(library_path, run, "resident")
    saved <- driver$checkpointInstalledApp(library_path, run, "qs2_duckdb")
    expect_identical(environment(native)$fixed_capacity, "resident")
    expect_identical(environment(saved)$fixed_capacity, "qs2_duckdb")
    expect_error(driver$checkpointInstalledApp(library_path, run, "limpa"), "invalid capacity")
    expect_identical(digest::digest(file = file.path(root, "tests", "testdata",
        "workflow-checkpoints", "scenarios.json"), algo = "sha256"),
        "cfdf7da2fc5a6123da269f85011550d8ffd2e8a611de071e0cd1a52a6f35d152")
})

test_that("capacity dispatcher stops immediately instead of repeating or scaling a failed member", {
    environment <- checkpointHarnessTestEnvironment()
    source_root <- normalizePath(test_path("..", ".."))
    environment$checkpointRunnerDirectory <- function() file.path(source_root, "tools", "profiling")
    root <- withr::local_tempdir()
    library_path <- file.path(root, "library")
    dir.create(library_path)
    dir.create(file.path(library_path, "R"))
    writeLines("installed identity fixture", file.path(library_path, "R", "MultiScholaR.rdb"))
    environment$inspectCheckpointLibrary <- function(...) list(package_path = library_path)
    environment$validateCheckpointOutputRoot <- function(...) root
    contract_path <- file.path(source_root, "tests", "testdata", "workflow-checkpoints", "gui-capacity.json")
    contract <- jsonlite::read_json(contract_path)
    options <- list(mode = "gui", route = "iq", diagnostics = FALSE,
        library = library_path, contract = contract_path,
        protocol = file.path(source_root, "tests", "testdata", "workflow-checkpoints", "protocol.json"),
        "output-root" = root)
    operations <- character()
    service <- function(operation, library_path, output_root, arguments) {
        operations <<- c(operations, operation)
        directory <- tempfile("member-", tmpdir = root)
        dir.create(directory)
        if (operation == "prepare-companions") {
            jsonlite::write_json(list(status = "prepared", selection = "small",
                contract_sha256 = digest::digest(file = contract_path, algo = "sha256"),
                files = list(report = "r", fasta = "f", design = "d")),
                file.path(directory, "inputs.json"), auto_unbox = TRUE)
        }
        list(exit_status = if (operation == "browser-baseline") 76L else 0L,
            member_dir = directory, unit = "owned-fixture.service")
    }
    expect_error(environment$runCheckpointGuiCapacity(options, contract, service), "member stopped")
    expect_identical(operations, c("prepare-companions", "browser-baseline"))
    expect_error(environment$runCheckpointGuiCapacity(options, contract, service), "already exists")
    expect_identical(operations, c("prepare-companions", "browser-baseline"))
})

test_that("reader runner validates fixed selectors and retains late failures", {
    driver <- new.env(parent = baseenv())
    sys.source(testthat::test_path("..", "..", "tools", "profiling",
        "workflow_checkpoint_dia_driver.R"), envir = driver)
    root <- withr::local_tempdir()
    jsonlite::write_json(list(status = "prepared"), file.path(root, "reader-inputs.json"))
    expect_true(driver$validateCheckpointReaderInput("reader-member", root, "small:fastduckdb:1"))
    expect_true(driver$validateCheckpointReaderInput("reader-member", root, "small:adapter_fastduckdb:1"))
    expect_error(driver$validateCheckpointReaderInput("reader-member", root,
        "small:adapter_fastduckdb:5"), "frozen reader selection")
    expect_error(driver$validateCheckpointReaderInput("reader-member", root, "huge:fastduckdb:50"),
        "frozen reader selection")
    expect_error(driver$validateCheckpointReaderInput("reader-build", root), "pinned")
    library_path <- dirname(find.package("MultiScholaR"))
    run_dir <- withr::local_tempdir()
    driver$checkCheckpointReader <- function(...) stop("late reader consumer failure")
    expect_error(driver$runCheckpointReaderOperation("reader-check", library_path, run_dir, root),
        "late reader consumer failure")
    receipt <- jsonlite::read_json(file.path(run_dir, "reader-operation.json"))
    expect_identical(receipt$status, "failed")
    expect_false(receipt$promotion_authorized)
    expect_false(file.exists(file.path(run_dir, "reader-member.json")))
})

recordedPilotTestEnvironment <- function() {
    driver <- new.env(parent = baseenv())
    sys.source(testthat::test_path("..", "..", "tools", "profiling",
        "workflow_checkpoint_dia_driver.R"), envir = driver)
    driver$recordedCheckpointApi <- function(library_path) {
        namespace <- asNamespace("MultiScholaR")
        \(name) get(name, envir = namespace, inherits = FALSE)
    }
    driver
}

makeRecordedPilotTestSession <- function(sequence = FALSE) {
    recipe <- readWorkflowCheckpointRecipe("small")
    peptide <- makeWorkflowCheckpointPeptide(recipe)
    if (sequence) {
        peptide@args$filterMinNumPeptidesPerProtein <- list(
            peptides_per_protein_cutoff = 2, peptidoforms_per_protein_cutoff = 2)
        peptide@args$rollUpPrecursorToPeptide <- list(core_utilisation = NA)
    }
    peptide@args$srlQvalueProteotypicPeptideClean <- list(
        input_matrix_column_ids = c("Run", "Precursor.Id", "Protein.Group", "Protein.Ids",
            "Stripped.Sequence", "Modified.Sequence", "Precursor.Charge",
            "Precursor.Quantity", "Precursor.Normalised"),
        exclude_contaminants = FALSE, qvalue_threshold = 0.01,
        global_qvalue_threshold = 0.01, global_pg_qvalue_threshold = 0.01,
        choose_only_proteotypic_peptide = 1, confidence_provenance_mode = "diann_main_report")
    manager <- WorkflowState$new()
    manager$saveState("raw_data_s4", peptide, list(), "raw recorded test")
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    config_env <- new.env(parent = emptyenv())
    config_env$config_list <- peptide@args
    timestamp <- as.POSIXct("2026-09-16 00:00:00", tz = "UTC")
    invisible(runPeptideQvalueApplyStep(workflow, 0.01, 0.01, TRUE,
        updateConfigParameterFn = function(...) updateConfigParameter(..., env = config_env),
        nowFn = \() timestamp, logInfoFn = \(...) invisible(NULL),
        logWarnFn = \(...) invisible(NULL)))
    manager$states$qvalue_filtered$timestamp <- timestamp
    if (sequence) {
        invisible(runProteinPeptideApplyStep(workflow, 2, 2,
            updateConfigParameterFn = function(...) updateConfigParameter(..., env = config_env),
            nowFn = \() timestamp, logInfoFn = \(...) invisible(NULL)))
        manager$states$protein_peptide_filtered$timestamp <- timestamp
        invisible(runPeptideRollupApplyStep(workflow, nowFn = \() timestamp,
            logInfoFn = \(...) invisible(NULL)))
        manager$states$precursor_rollup$timestamp <- timestamp
    }
    protein <- makeWorkflowCheckpointProtein(recipe)
    manager$saveState("correlation_filtered", protein, list(), "final recorded test")
    list(r6_complete_states = manager$states, r6_state_history = manager$state_history,
        r6_current_state_name = manager$current_state, current_s4_object = protein,
        correlation_filtered_s4 = protein, export_timestamp = timestamp,
        final_protein_count = recipe$n_proteins, final_sample_count = recipe$n_samples,
        normalization_method = "none", ruv_mode = "skip", workflow_type = "DIA",
        qc_params = list(fixed = "recorded test settings"),
        contrasts_tbl = data.frame(friendly_names = "treated_vs_control"))
}

test_that("three QC text handoffs do not retain their S4 payload", {
    object <- methods::new("PeptideQuantitativeData")
    object@args <- list(retention_probe = double(1000000))
    for (writer in list(updatePeptideQvalueOutputs, updateProteinPeptideOutputs,
        updatePeptideRollupOutputs)) {
        output <- new.env(parent = emptyenv())
        result <- list(resultText = "complete", filteredS4 = object, rolledUpS4 = object)
        writer(output, function(value) invisible(NULL), result, "DIA", "test",
            updateProteinFilteringFn = function(...) list())
        renderer <- as.list(output)[[1L]]
        expect_lt(as.numeric(lobstr::obj_size(renderer)), 1000000)
        expect_identical(shiny::isolate(renderer()), "complete")
        driver <- recordedPilotTestEnvironment()
        expect_false(driver$recordedQcSequenceHasPayload(renderer))
        expect_true(driver$recordedQcSequenceHasPayload(list(leaked_result = result)))
    }
})

test_that("QC ownership inspection descends through grid-unit storage", {
    driver <- recordedPilotTestEnvironment()
    layout <- grid::unit.c(grid::unit(1, "cm"), grid::unit(1, "null"))
    expect_true(is.list(layout[[1L]]))
    expect_false(driver$recordedQcSequenceHasPayload(list(layout = layout)))
    expect_false(driver$recordedQcSequenceHasPayload(
        grid::grobTree(grid::rectGrob(width = layout[1L]))))
    object <- methods::new("PeptideQuantitativeData")
    expect_true(driver$recordedQcSequenceHasPayload(list(layout = layout,
        hidden = list(assay = object))))
    expect_true(driver$recordedQcSequenceHasPayload(structure(list(object),
        class = "qc_classed_list_probe")))
})

test_that("real sequential QC preserves all stages and strict checkpoint release", {
    skip_if_not_installed("qs2")
    driver <- recordedPilotTestEnvironment()
    prepared <- withr::local_tempdir()
    driver$splitRecordedCheckpointSession(makeRecordedPilotTestSession(sequence = TRUE), prepared,
        driver$recordedCheckpointApi("unused"), strrep("c", 64L))
    settings <- withr::local_tempdir()
    driver$prepareRecordedQcSequence("unused", settings, prepared)
    sessions <- list()
    for (treatment in c("resident", "qs2", "qs2_duckdb")) {
        member <- withr::local_tempdir()
        driver$executeRecordedQcSequence("unused", member, settings, treatment, treatment == "qs2",
            function() list(peak_bytes = 0, current_bytes = 0))
        receipt <- jsonlite::read_json(file.path(member, "qc-sequence-member.json"))
        expect_identical(receipt$status, "complete")
        expect_identical(unlist(receipt$completed_steps),
            c("qvalue_filtered", "protein_peptide_filtered", "precursor_rollup"))
        expect_identical(receipt$open_resources, 0L)
        expect_length(receipt$files, 40L)
        for (name in unlist(receipt$completed_steps)) {
            observation <- receipt$observations[[name]]
            expect_true(observation$current_released)
            expect_true(observation$config_released)
            if (treatment != "resident") {
                expect_identical(observation$resident_history_entries, 0L)
                expect_identical(observation$resident_audit_entries, 0L)
            }
        }
        verified <- driver$verifyRecordedQcSequence("unused", withr::local_tempdir(), member)
        expect_identical(verified$status, "passed_real_qc_sequence")
        sessions[[treatment]] <- readRDS(file.path(member, "export", "filtered_session_data_latest.rds"))
        if (treatment != "resident") {
            refs <- readRDS(file.path(member, "qc-sequence-references.rds"))
            expect_named(refs, c("initial", "raw_data_s4", "qvalue_filtered",
                "protein_peptide_filtered", "precursor_rollup"))
        }
        unlink(file.path(member, "plots", "4_precursor_rollup_combined_plots.png"))
        expect_error(driver$recordedQcSequenceInventory(member), "required sequential QC output")
    }
    expect_identical(sessions$resident, sessions$qs2)
    expect_identical(sessions$resident, sessions$qs2_duckdb)
})

test_that("sequential QC selectors do not replace the old recorded operation", {
    harness <- checkpointHarnessTestEnvironment()
    input <- withr::local_tempdir()
    jsonlite::write_json(list(status = "prepared_qc_sequence"),
        file.path(input, "qc-sequence-inputs.json"), auto_unbox = TRUE)
    expect_true(harness$validateRecordedCheckpointInput("recorded-qc", input,
        "sequence:resident:diagnostic"))
    expect_true(harness$validateRecordedCheckpointInput("recorded-qc", input,
        "sequence:qs2:timing"))
    expect_true(harness$validateRecordedCheckpointInput("recorded-qc", input,
        "sequence:qs2_duckdb:timing"))
    expect_error(harness$validateRecordedCheckpointInput("recorded-qc", input,
        "sequence:qs2:seven-blocks"), "unsupported sequential")
    expect_error(harness$validateRecordedCheckpointInput("recorded-qc", input), "manifest")
})

test_that("sequential QC failures retain committed checkpoints and close their owners", {
    skip_if_not_installed("qs2")
    driver <- recordedPilotTestEnvironment()
    api <- driver$recordedCheckpointApi("unused")
    prepared <- withr::local_tempdir()
    driver$splitRecordedCheckpointSession(makeRecordedPilotTestSession(sequence = TRUE), prepared,
        api, strrep("d", 64L))
    settings <- withr::local_tempdir()
    driver$prepareRecordedQcSequence("unused", settings, prepared)
    for (failure in c("save", "outputs", "interrupt")) {
        store <- NULL
        driver$recordedCheckpointApi <- function(library_path) function(name) {
            if (name == "newWorkflowEntityStore") return(function(...) {
                store <<- api(name)(...)
                store
            })
            if (failure == "save" && name == "saveWorkflowCheckpoint") return(function(...) {
                if (store$count == 2L) stop("injected checkpoint save failure")
                api(name)(...)
            })
            if (failure != "save" && name == "updatePeptideQvalueOutputs") return(function(...) {
                if (failure == "interrupt") stop(structure(list(message = "injected cancellation"),
                    class = c("interrupt", "condition")))
                stop("injected required output failure")
            })
            api(name)
        }
        member <- withr::local_tempdir()
        error <- tryCatch(driver$executeRecordedQcSequence("unused", member, settings,
            "qs2", FALSE, function() list(peak_bytes = 0)), error = identity, interrupt = identity)
        expect_s3_class(error, if (failure == "interrupt") "interrupt" else "error")
        expect_true(store$closed)
        expect_identical(store$open_resources, 0L)
        receipt <- jsonlite::read_json(file.path(member, "qc-sequence-member.json"))
        expect_identical(receipt$status, if (failure == "interrupt") "interrupted" else "failed")
        expect_length(receipt$completed_steps, 0L)
        refs <- readRDS(file.path(member, "qc-sequence-references.rds"))
        expect_length(refs, if (failure == "save") 2L else 3L)
    }
})

test_that("recorded history replay releases owners and preserves native exports", {
    driver <- recordedPilotTestEnvironment()
    session <- makeRecordedPilotTestSession()
    prepared <- withr::local_tempdir()
    manifest <- driver$splitRecordedCheckpointSession(session, prepared,
        driver$recordedCheckpointApi("unused"), strrep("a", 64L))
    expect_identical(manifest$state_names, names(session$r6_complete_states))
    expect_false(manifest$promotion_authorized)
    outputs <- list()
    treatments <- c("resident", "saved")
    if (requireNamespace("qs2", quietly = TRUE)) treatments <- c(treatments, "qs2")
    for (treatment in treatments) {
        member <- withr::local_tempdir()
        result <- driver$executeRecordedCheckpointReplay("unused", member, prepared,
            treatment, \() list(peak_bytes = 0), profile_owners = treatment == "saved")
        expect_identical(result$status, "complete")
        expect_identical(result$history_entries, 5L)
        expect_identical(result$open_resources, 0L)
        expect_identical(result$profile_owners, treatment == "saved")
        expect_gte(result$retained_observation$observation_elapsed_seconds, 0)
        if (treatment == "resident") expect_null(result$retained_owned_bytes)
        expect_false(result$independently_verified)
        verify <- driver$verifyRecordedCheckpointReplay("unused", withr::local_tempdir(), member)
        expect_identical(verify$status, "passed_recorded_history")
        expect_identical(verify$checked_states, 5L)
        expect_true(verify$full_metadata_verified)
        outputs[[treatment]] <- verify$legacy_sha256
        saveRDS(list(changed = TRUE), file.path(member, "export", "qc_params.RDS"))
        expect_error(driver$verifyRecordedCheckpointReplay("unused", withr::local_tempdir(), member),
            "export metadata differs")
    }
    expect_identical(outputs$resident, outputs$saved)
    if (!is.null(outputs$qs2)) expect_identical(outputs$resident, outputs$qs2)
    base_api <- driver$recordedCheckpointApi("unused")
    failed_store <- NULL
    driver$recordedCheckpointApi <- function(library_path) function(name) {
        if (name == "saveProtNormExportArtifacts") return(function(...) stop("late export failure"))
        if (name == "newWorkflowEntityStore") return(function(...) {
            failed_store <<- base_api(name)(...)
            failed_store
        })
        base_api(name)
    }
    failed_member <- withr::local_tempdir()
    expect_error(driver$executeRecordedCheckpointReplay("unused", failed_member, prepared,
        "saved", \() list(peak_bytes = 0)), "late export failure")
    expect_true(failed_store$closed)
    expect_identical(failed_store$open_resources, 0L)
    expect_false(file.exists(file.path(failed_member, "recorded-member.json")))
    expect_length(list.dirs(file.path(failed_member, ".workflow-entities", "entities"),
        recursive = FALSE), 5L)
    changed <- manifest
    changed$records$raw_data_s4$bytes <- 1
    expect_error(driver$readRecordedCheckpoint(prepared, changed, "raw_data_s4"), "integrity")
})

test_that("recorded preparation preserves distinct per-section worker specifications", {
    driver <- recordedPilotTestEnvironment()
    api <- driver$recordedCheckpointApi("unused")
    first <- structure(list(new.env(parent = emptyenv()), new.env(parent = emptyenv())),
        class = "multidplyr_cluster")
    second <- structure(list(new.env(parent = emptyenv()), new.env(parent = emptyenv()),
        new.env(parent = emptyenv())),
        class = "multidplyr_cluster")
    config <- list(globalParameters = list(number_of_cpus = 2),
        one = list(core_utilisation = first, threshold = 0.1),
        two = list(core_utilisation = second, method = "unchanged"))
    expect_error(api("splitDIAWorkerConfig")(config), "distinct worker clusters")
    record <- list(data = makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small")),
        config = config)
    prepared <- driver$prepareRecordedCheckpointRecord(record, api)
    expect_identical(prepared$data, record$data)
    expect_identical(prepared$config$one$core_utilisation, api("newDIAWorkerSpec")(2))
    expect_identical(prepared$config$two$core_utilisation, api("newDIAWorkerSpec")(3))
    expect_identical(prepared$config$globalParameters, config$globalParameters)
    expect_identical(prepared$config$one$threshold, config$one$threshold)
    expect_identical(prepared$config$two$method, config$two$method)
    expect_identical(record$config$one$core_utilisation, first)
})

test_that("recorded active QC profiles the complete operation and verifies its result", {
    skip_if_not(isTRUE(capabilities("profmem")))
    driver <- recordedPilotTestEnvironment()
    prepared <- withr::local_tempdir()
    driver$splitRecordedCheckpointSession(makeRecordedPilotTestSession(), prepared,
        driver$recordedCheckpointApi("unused"), strrep("b", 64L))
    member <- withr::local_tempdir()
    result <- driver$executeRecordedActiveQc("unused", member, prepared, \() list(peak_bytes = 0))
    expect_identical(result$status, "complete")
    expect_named(result$observations, c("entry", "active_input", "method_entry",
        "method_return_before_audit_commit", "after_audit_commit", "after_release"))
    expect_gt(file.size(result$allocation_profile), 0)
    expect_true(any(grepl("apply_step",
        readLines(result$allocation_profile), fixed = TRUE)))
    expect_identical(driver$verifyRecordedActiveQc("unused", withr::local_tempdir(),
        member)$status, "passed_recorded_active_qc")
})

test_that("recorded pilot guard checks peak and host headroom without allocations", {
    runner <- testthat::test_path("..", "..", "tools", "profiling", "run_workflow_checkpoint_proof.sh")
    lines <- readLines(runner)
    guard <- paste(lines[seq_len(grep("^proof_script_dir=", lines) - 1L)], collapse = "\n")
    check <- function(arguments) processx::run("bash", c("-c", paste(guard,
        paste("checkpointRecordedGuard", arguments), sep = "\n")), error_on_status = FALSE)
    expect_identical(check("10737418240 0 0 start")$status, 0L)
    expect_identical(check("10737418239 0 0 start")$status, 1L)
    expect_identical(check("2147483648 8053063679 8053063679 running")$status, 0L)
    expect_identical(check("2147483647 0 0 running")$status, 1L)
    expect_identical(check("10737418240 0 8053063680 running")$status, 1L)
    expect_identical(check("10737418240 8053063680 0 running")$status, 1L)
    expect_match(check("10737418240 unknown 0 running")$stdout, "missing_memory")
})

test_that("recorded dispatch and validation refuse incomplete pilot inputs", {
    harness <- checkpointHarnessTestEnvironment()
    harness$runRecordedCheckpointOperation <- function(...) list(...)
    result <- harness$checkpointProofMain(c("--internal-recorded-member", "library",
        "run", "prepared", "saved"))
    expect_identical(result, list("recorded-member", "library", "run", "prepared", "saved"))
    timed <- harness$checkpointProofMain(c("--internal-recorded-timing", "library",
        "run", "prepared", "resident"))
    expect_identical(timed, list("recorded-timing", "library", "run", "prepared", "resident"))
    expect_error(harness$validateRecordedCheckpointInput("recorded-member",
        withr::local_tempdir(), "saved"), "manifest")
    expect_error(harness$validateRecordedCheckpointInput("recorded-verify",
        withr::local_tempdir()), "completed member")
})

test_that("halved DIA preparation preserves exact rows and balanced paired samples", {
    environment <- checkpointHarnessTestEnvironment()
    root <- withr::local_tempdir()
    source <- file.path(root, "source")
    output <- file.path(root, "browser-smoke.ABC123")
    dir.create(source)
    dir.create(output)
    runs <- c(paste0("a", 1:4), paste0("b", 1:4))
    design <- data.frame(Run = runs, group = rep(c("control", "treated"), each = 4L),
        replicates = rep(1:2, each = 2L, times = 2L),
        tech_rep_group = rep(c("control_1", "control_2", "treated_1", "treated_2"),
            each = 2L))
    report <- expand.grid(Run = runs, Protein.Group = c("PG1", "PG2"),
        Stripped.Sequence = c("AA", "BB"), Precursor.Charge = c(2L, 3L),
        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    report <- report[c("Run", "Protein.Group", "Stripped.Sequence", "Precursor.Charge")]
    report$Protein.Ids <- sub("PG", "P", report$Protein.Group)
    report$Modified.Sequence <- ifelse(report$Precursor.Charge == 2L,
        report$Stripped.Sequence, paste0("M", report$Stripped.Sequence))
    report$Precursor.Id <- paste(report$Modified.Sequence, report$Precursor.Charge, sep = "_")
    report <- report[c("Run", "Protein.Group", "Protein.Ids", "Stripped.Sequence",
        "Modified.Sequence", "Precursor.Charge", "Precursor.Id")]
    paths <- c(report = file.path(source, "report.tsv"),
        design = file.path(source, "design.tsv"), fasta = file.path(source, "proteins.fasta"))
    utils::write.table(report, paths[["report"]], sep = "\t", quote = FALSE, row.names = FALSE)
    utils::write.table(design, paths[["design"]], sep = "\t", quote = FALSE, row.names = FALSE)
    writeLines(c(">P1", "AAAA", ">P2", "BBBB"), paths[["fasta"]])
    hashes <- vapply(paths, \(path) digest::digest(file = path, algo = "sha256"), character(1))
    result <- environment$prepareCheckpointHalvedInput(paths[["report"]], paths[["design"]],
        paths[["fasta"]], output, list(report_sha256 = hashes[["report"]],
            design_sha256 = hashes[["design"]], fasta_sha256 = hashes[["fasta"]],
            samples_per_group = 2L, proteins = 2L, peptides_per_protein = 2L))
    expect_true(result$diagnostic_only)
    expect_false(result$promotion_authorized)
    expect_identical(result$precursor_rows, 16L)
    expect_identical(result$proteins, 2L)
    expect_identical(result$peptides, 4L)
    expect_setequal(result$selected_runs, c("a1", "a2", "b1", "b2"))
    source_lines <- readLines(paths[["report"]])
    source_rows <- source_lines[-1L]
    expected <- c(source_lines[[1L]], source_rows[report$Run %in% result$selected_runs &
        report$Precursor.Charge == 2L])
    expect_identical(readLines(result$files$report), expected)
    expect_identical(digest::digest(file = result$files$report, algo = "sha256"),
        result$sha256$report)
    expect_identical(readLines(result$files$design),
        readLines(paths[["design"]])[c(TRUE, design$Run %in% result$selected_runs)])
    expect_error(environment$prepareCheckpointHalvedInput(paths[["report"]], paths[["design"]],
        paths[["fasta"]], output, list(report_sha256 = hashes[["report"]],
            design_sha256 = hashes[["design"]], fasta_sha256 = hashes[["fasta"]],
            samples_per_group = 2L, proteins = 2L, peptides_per_protein = 2L)), "overwrite")
})

test_that("summary browser actions use existing controls and require all three acknowledgements", {
    environment <- new.env(parent = baseenv())
    sys.source(testthat::test_path("..", "..", "tools", "profiling",
        "workflow_checkpoint_dia_driver.R"), envir = environment)
    selected <- list()
    clicked <- list()
    environment$checkpointSetInput <- function(driver, id, value) {
        selected[[id]] <<- value
    }
    environment$checkpointClickAndWaitText <- function(driver, id, selector, text, ...) {
        clicked[[id]] <<- list(selector = selector, text = text)
    }
    result <- environment$runCheckpointSummaryActions(list())
    prefix <- "proteomics_workflow-session_summary-"
    expect_identical(selected[["proteomics_workflow-workflow_tabs"]], "session_summary")
    expect_identical(selected[[paste0(prefix, "experiment_label")]], "DIA Experiment")
    expect_identical(names(clicked), paste0(prefix, result$actions))
    expect_match(clicked[[paste0(prefix, "save_workflow_args")]]$text, "successfully")
    expect_match(clicked[[paste0(prefix, "copy_to_publication")]]$text, "Files copied")
    expect_match(clicked[[paste0(prefix, "export_session_state")]]$text, "Session state exported")
})

test_that("summary file check rejects an incomplete action rather than calling it complete", {
    environment <- checkpointHarnessTestEnvironment()
    root <- withr::local_tempdir()
    expect_error(environment$checkCheckpointSummaryFiles(root, "DIA Experiment"), "session-state")
    project <- file.path(root, "dia_baseline_smoke")
    source <- file.path(project, "scripts", "proteomics")
    dir.create(source, recursive = TRUE)
    saveRDS(list(experiment_label = "DIA Experiment", workflow_args_saved = FALSE,
        files_copied = FALSE), file.path(source, "session_state_2026-09-16.RDS"))
    expect_error(environment$checkCheckpointSummaryFiles(root, "DIA Experiment"), "absent or empty")
})

test_that("selected publication figure copies match their source files", {
    environment <- checkpointHarnessTestEnvironment()
    root <- withr::local_tempdir()
    sections <- c("Volcano_Plots", "NumSigDaMolecules", "Heatmap")
    for (section in sections) {
        source <- file.path(root, "results", "proteomics", "publication_graphs", section)
        target <- file.path(root, "results_summary", "proteomics", "Publication_figures", section)
        dir.create(source, recursive = TRUE)
        dir.create(target, recursive = TRUE)
        writeLines("source fixture", file.path(source, "plot.png"))
        writeLines("source fixture", file.path(target, "plot.png"))
    }
    result <- environment$checkCheckpointPublicationFigures(root)
    expect_setequal(names(result), sections)
    target <- file.path(root, "results_summary", "proteomics", "Publication_figures",
        "Heatmap", "plot.png")
    writeLines("changed copy", target)
    expect_error(environment$checkCheckpointPublicationFigures(root), "differs: Heatmap")
})

test_that("required IQ QC inventory rejects absent empty and duplicate automatic plots", {
    environment <- checkpointHarnessTestEnvironment()
    root <- withr::local_tempdir()
    expect_error(environment$checkCheckpointIqQcOutputs(root), "2_qval_filtered_combined_plots")
    project <- file.path(root, "dia_baseline_smoke")
    source_dir <- file.path(project, "scripts", "proteomics")
    results <- file.path(project, "results", "proteomics")
    protein <- file.path(results, "protein_qc")
    qc <- file.path(results, "publication_graphs", "filtering_qc", "fixture_time")
    paths <- c(file.path(root, "qc-plots.pdf"),
        file.path(results, "peptide_qc", "peptide_values_imputed.tsv"),
        file.path(source_dir, c("aa_seq_tbl_final.RDS", "fasta_metadata.RDS", "design_matrix.tab",
            "ruv_optimization_results.RDS", "filtered_session_data_latest.rds", "filtered_session_summary.txt")),
        file.path(protein, c("iq_output_file.txt", "remove_proteins_with_only_one_rep.tsv",
            "normalized_protein_matrix_pre_ruv.tsv", "normalised_results_cln_with_replicates.tsv",
            "normalised_results_cln_with_replicates.RDS", "proteomics_composite_QC_figure.rds",
            "proteomics_composite_QC_figure.pdf", "proteomics_composite_QC_figure.png")))
    for (stage in c("pre_norm", "post_norm")) for (plot in c("pca", "density", "rle", "correlation")) {
        paths <- c(paths, file.path(protein, paste0(stage, "_", plot, ".png")))
    }
    for (step in c("2_qval_filtered", "3_protein_peptide_filtered", "4_precursor_rollup",
        "5_intensity_filtered", "6_sample_filtered", "7_replicate_filtered", "8_imputed")) {
        for (plot in c("combined_plots", "proteins_total", "proteins_per_run",
            "peptides_total", "peptides_per_protein", "peptides_per_run")) {
            paths <- c(paths, file.path(qc, paste0(step, "_", plot, ".", c("png", "pdf"))))
        }
    }
    for (path in paths) {
        dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
        writeLines("presence fixture only; not valid scientific content", path)
    }
    result <- environment$checkCheckpointIqQcOutputs(root)
    expect_identical(result$status, "passed_iq_qc_presence_only")
    expect_false(result$promotion_authorized)
    expect_setequal(result$files, paths)
    correlation_rds <- file.path(protein, "normalised_results_cln_with_replicates.RDS")
    expect_true(file.rename(correlation_rds, paste0(correlation_rds, ".compatibility")))
    writeLines("native latest presence fixture", file.path(source_dir, "filtered_session_latest.json"))
    native <- environment$checkCheckpointIqQcOutputs(root, native_session = TRUE, native_first = TRUE)
    expect_false(correlation_rds %in% native$files)
    expect_true(file.path(source_dir, "filtered_session_latest.json") %in% native$files)
    expect_error(environment$checkCheckpointIqQcOutputs(root), "normalised_results_cln_with_replicates.RDS")
    expect_true(file.rename(paste0(correlation_rds, ".compatibility"), correlation_rds))
    writeBin(raw(), file.path(root, "qc-plots.pdf"))
    expect_error(environment$checkCheckpointIqQcOutputs(root), "qc-plots.pdf")
    writeLines("presence fixture", file.path(root, "qc-plots.pdf"))
    duplicate <- file.path(dirname(qc), "second_time", "2_qval_filtered_combined_plots.png")
    dir.create(dirname(duplicate))
    writeLines("duplicate", duplicate)
    expect_error(environment$checkCheckpointIqQcOutputs(root), "ambiguous")
})

test_that("scenario validation preserves all frozen members without requiring its own G0", {
    environment <- checkpointHarnessTestEnvironment()
    contract <- jsonlite::read_json(testthat::test_path("..", "testdata", "workflow-checkpoints", "scenarios.json"))
    result <- environment$validateCheckpointScenarioContract(contract, "iq")
    expect_false(result$baseline_requires_prior_g0)
    expect_false(result$promotion_authorized)
    expect_length(result$scenario_ids, 26L)
    expect_length(result$timing_recipe_ids, 6L)
    expect_identical(result$capacity_only_recipe_ids, "capacity")
    changed <- contract
    changed$recipes <- changed$recipes[-2L]
    expect_error(environment$validateCheckpointScenarioContract(changed, "iq"), "inventory")
    changed <- contract
    changed$recipes[[1L]]$precursor_rows <- 1
    expect_error(environment$validateCheckpointScenarioContract(changed, "iq"), "dimensions")
    changed <- contract
    changed$recipes[[7L]]$frozen_for_timing <- TRUE
    expect_error(environment$validateCheckpointScenarioContract(changed, "iq"), "eligibility")
    changed <- contract
    changed$scenario_matrix[[1L]]$recipes <- list("unknown")
    expect_error(environment$validateCheckpointScenarioContract(changed, "iq"), "references")
    changed <- contract
    changed$scenario_matrix <- changed$scenario_matrix[-1L]
    expect_error(environment$validateCheckpointScenarioContract(changed, "iq"), "coverage")
    expect_error(environment$validateCheckpointScenarioContract(contract, "unknown"), "route")
})

test_that("baseline input resolution binds every recipe and refuses partial inventories", {
    environment <- checkpointHarnessTestEnvironment()
    contract <- jsonlite::read_json(testthat::test_path("..", "testdata", "workflow-checkpoints", "scenarios.json"))
    root <- withr::local_tempdir()
    expect_error(environment$resolveCheckpointBaselineInputs(contract, root), "receipt: small")
    for (recipe in contract$recipes) {
        if (!recipe$frozen_for_timing) next
        directory <- file.path(root, recipe$id)
        dir.create(directory)
        paths <- stats::setNames(file.path(directory, c("report.tsv", "proteins.fasta", "design.tsv")),
            c("report", "fasta", "design"))
        for (path in paths) writeLines("metadata-resolution fixture, not scientific data", path)
        receipt <- list(kind = "small_fixture_preparation_not_timing", recipe = recipe,
            files = as.list(paths), sha256 = as.list(vapply(paths,
                \(path) digest::digest(file = path, algo = "sha256"), character(1))))
        jsonlite::write_json(receipt, file.path(directory, "inputs.json"), auto_unbox = TRUE)
    }
    bound <- environment$resolveCheckpointBaselineInputs(contract, root)
    expect_identical(names(bound), c("small", "representative", "tall", "wide", "annotation_heavy", "cycles"))
    expect_false("capacity" %in% names(bound))
    archived <- file.path(root, "browser-smoke.ABC123")
    dir.create(archived)
    expect_true(file.rename(file.path(root, "small", "inputs.json"),
        file.path(archived, "inputs.json")))
    discovered <- environment$resolveCheckpointBaselineInputs(contract, root)
    expect_identical(discovered$small$preparation_receipt,
        normalizePath(file.path(archived, "inputs.json")))
    expect_true(file.copy(file.path(archived, "inputs.json"),
        file.path(root, "small", "inputs.json")))
    expect_error(environment$resolveCheckpointBaselineInputs(contract, root), "ambiguous.*small")
    expect_true(file.rename(file.path(archived, "inputs.json"),
        file.path(archived, "retained_original.json")))
    path <- file.path(root, "small", "inputs.json")
    receipt <- jsonlite::read_json(path)
    receipt$recipe$seed <- 1
    jsonlite::write_json(receipt, path, auto_unbox = TRUE)
    expect_error(environment$resolveCheckpointBaselineInputs(contract, root), "match frozen recipe")
})

test_that("baseline dispatcher is serial, records failures and never retries or skips members", {
    environment <- checkpointHarnessTestEnvironment()
    protocol <- jsonlite::read_json(testthat::test_path("..", "testdata", "workflow-checkpoints", "protocol.json"),
        simplifyVector = TRUE)
    ids <- c("small", "representative", "tall", "wide", "annotation_heavy", "cycles")
    inputs <- stats::setNames(lapply(ids, \(id) list(id = id)), ids)
    calls <- character()
    recorded <- list()
    run <- function(input, block) {
        id <- paste(input$id, block, sep = ":")
        # Previous receipt must be durable before the next member starts.
        expect_length(recorded, length(calls))
        calls <<- c(calls, id)
        list(status = "complete_verified", run_id = id)
    }
    record <- function(id, block, result) recorded[[paste(id, block, sep = ":")]] <<- result
    result <- environment$dispatchCheckpointBaseline(inputs, protocol, run, record)
    expect_length(result, 42L)
    expect_identical(names(result), calls)
    expect_identical(calls[1:8], c(paste0("small:", 1:7), "representative:1"))
    calls <- character()
    recorded <- list()
    fail <- function(input, block) {
        calls <<- c(calls, paste(input$id, block, sep = ":"))
        stop("fixture service failure")
    }
    expect_error(environment$dispatchCheckpointBaseline(inputs, protocol, fail, record), "did not complete")
    expect_length(calls, 1L)
    expect_match(recorded[[1L]]$error, "fixture service failure")
    expect_error(environment$dispatchCheckpointBaseline(inputs[-2L], protocol, run, record), "full frozen")
    calls <- character()
    recorded <- list()
    duplicate <- function(input, block) list(status = "complete_verified", run_id = "reused")
    expect_error(environment$dispatchCheckpointBaseline(inputs, protocol, duplicate, record), "did not complete")
    expect_length(recorded, 2L)
    expect_identical(recorded[[2L]]$status, "failed")
})

test_that("baseline service handoff verifies only matching completed members without gate promotion", {
    environment <- checkpointHarnessTestEnvironment()
    root <- withr::local_tempdir()
    member <- file.path(root, "browser-smoke.ABC123")
    verification <- file.path(root, "browser-smoke.DEF456")
    dir.create(member)
    dir.create(verification)
    hashes <- list(report = strrep("a", 64), fasta = strrep("b", 64), design = strrep("c", 64))
    input <- list(recipe = list(id = "small"), files = list(report = "r", fasta = "f", design = "d"), sha256 = hashes)
    terminal <- list(service_succeeded = TRUE, status = "service_tree_interval_only",
        unit = "multischolar-dia003-browser-smoke.ABC123.service", elapsed_seconds = 1,
        invocation_id = strrep("d", 32))
    browser <- list(status = "passed_summary_ui_only", installation = list(library_path = normalizePath(root)),
        input_sha256 = hashes, inputs = list(qc_defaults = TRUE), diagnostics = FALSE)
    oracle <- list(status = "passed_small_summary_file_checks", member_dir = member,
        required_iq_qc = list(status = "passed_iq_qc_presence_only"),
        summary_files = list(status = "passed_small_summary_files",
            publication_figures = list(Volcano_Plots = list(), NumSigDaMolecules = list(),
                Heatmap = list())))
    write_receipts <- function() {
        jsonlite::write_json(terminal, file.path(member, "terminal-interval.json"), auto_unbox = TRUE)
        jsonlite::write_json(browser, file.path(member, "smoke-receipt.json"), auto_unbox = TRUE)
        jsonlite::write_json(oracle, file.path(verification, "verification-receipt.json"), auto_unbox = TRUE)
    }
    operations <- character()
    service_exit <- 0L
    verify_exit <- 0L
    service <- function(operation, library_path, output_root, arguments) {
        operations <<- c(operations, operation)
        if (operation == "browser-baseline") {
            expect_identical(arguments, c("r", "f", "d"))
            return(list(exit_status = service_exit, member_dir = member, unit = terminal$unit))
        }
        expect_identical(arguments, member)
        list(exit_status = verify_exit, member_dir = verification)
    }
    run <- function() environment$runCheckpointBaselineMember(input, 1L, root, root, service)
    write_receipts()
    result <- run()
    expect_identical(operations, c("browser-baseline", "verify-da"))
    expect_identical(result$status, "complete_verified")
    expect_identical(result$run_id, terminal$invocation_id)
    expect_false(result$promotion_authorized)
    service_exit <- 75L
    operations <- character()
    result <- run()
    expect_identical(result$status, "failed")
    expect_identical(operations, "browser-baseline")
    service_exit <- 0L
    browser$input_sha256$report <- strrep("e", 64)
    write_receipts()
    operations <- character()
    result <- run()
    expect_match(result$error, "binding differs")
    expect_identical(operations, "browser-baseline")
    browser$input_sha256 <- hashes
    oracle$member_dir <- root
    write_receipts()
    result <- run()
    expect_match(result$error, "receipt binding differs")
    oracle$member_dir <- member
    write_receipts()
    verify_exit <- 1L
    result <- run()
    expect_match(result$error, "verification failed")
    verify_exit <- 0L
    oracle$summary_files$publication_figures$Heatmap <- NULL
    write_receipts()
    result <- run()
    expect_match(result$error, "required baseline output inventory")
})

test_that("utility dispatch is serial, paired and stops after a recorded failure", {
    environment <- checkpointHarnessTestEnvironment()
    protocol <- jsonlite::read_json(testthat::test_path("..", "testdata",
        "workflow-checkpoints", "protocol.json"), simplifyVector = TRUE)
    ids <- c("small", "representative", "tall", "wide", "annotation_heavy", "cycles")
    recipes <- stats::setNames(lapply(ids, \(id) list(id = id)), ids)
    calls <- list()
    recorded <- list()
    run <- function(recipe, block, treatment, variant) {
        calls[[length(calls) + 1L]] <<- c(recipe$id, as.character(block), treatment, variant)
        list(status = "complete", recipe = recipe$id, block = block,
            treatment = treatment, variant = variant, legacy_sha256 = strrep("a", 64L),
            prepared_inputs_sha256 = strrep("b", 64L))
    }
    record <- function(result) recorded[[length(recorded) + 1L]] <<- result
    result <- environment$dispatchCheckpointUtility(recipes, protocol, run, record)
    expect_length(result, 6L * 3L * 7L * 2L)
    expect_length(recorded, length(result))
    expect_identical(calls[[1L]][3:4], c("resident", "independent"))
    expect_identical(calls[[3L]][3:4], c("saved", "independent"))
    expect_identical(calls[[5L]][3:4], c("resident", "independent"))
    calls <- list()
    recorded <- list()
    qs2 <- environment$dispatchCheckpointUtility(recipes, protocol, run, record,
        candidate = "qs2")
    expect_length(qs2, 5L * 3L * 7L * 2L)
    expect_false(any(vapply(calls, \(item) item[[1L]] == "small", logical(1))))
    expect_setequal(vapply(calls, \(item) item[[1L]], character(1)), ids[-1L])
    expect_identical(calls[[2L]][3L], "qs2")
    expect_identical(calls[[3L]][3L], "qs2")
    expect_false(any(vapply(calls, \(item) item[[3L]] == "saved", logical(1))))
    expect_error(environment$dispatchCheckpointUtility(recipes, protocol, run, record,
        candidate = "unknown"), "invalid utility candidate")
    mismatched <- function(recipe, block, treatment, variant) {
        result <- run(recipe, block, treatment, variant)
        if (identical(treatment, "qs2")) result$legacy_sha256 <- strrep("c", 64L)
        result
    }
    recorded <- list()
    expect_error(environment$dispatchCheckpointUtility(recipes, protocol, mismatched,
        record, candidate = "qs2"), "paired legacy outputs or prepared inputs differ")
    expect_length(recorded, 2L)

    calls <- list()
    recorded <- list()
    fail <- function(recipe, block, treatment, variant) {
        list(status = "failed", error = "injected member failure", recipe = recipe$id,
            block = block, treatment = treatment, variant = variant)
    }
    expect_error(environment$dispatchCheckpointUtility(recipes, protocol, fail, record),
        "did not complete")
    expect_length(recorded, 1L)
    expect_identical(recorded[[1L]]$error, "injected member failure")
    expect_error(environment$dispatchCheckpointUtility(recipes[-1L], protocol, run, record),
        "full frozen")
})

test_that("utility histories preserve sharing, codecs, ten cycles and failure cleanup", {
    skip_if_not_installed("qs2", minimum_version = "0.2.2")
    root <- withr::local_tempdir()
    library_path <- file.path(root, "library")
    dir.create(file.path(library_path, "MultiScholaR"), recursive = TRUE)
    source_root <- normalizePath(testthat::test_path("..", ".."))
    files <- c("tools/profiling/workflow_checkpoint_dia_driver.R",
        "tools/profiling/run_workflow_checkpoint_proof.R",
        "tools/profiling/run_workflow_checkpoint_proof.sh",
        "tests/testthat/helper-workflow-checkpoint-fixtures.R")
    for (path in files) {
        dir.create(dirname(file.path(root, path)), recursive = TRUE, showWarnings = FALSE)
        expect_true(file.copy(file.path(source_root, path), file.path(root, path)))
    }
    contract <- file.path(root, "tests/testdata/workflow-checkpoints/scenarios.json")
    dir.create(dirname(contract), recursive = TRUE)
    recipe <- readWorkflowCheckpointRecipe("small")
    recipe$cycles <- 10L
    jsonlite::write_json(list(recipes = list(recipe)), contract, auto_unbox = TRUE)
    driver <- new.env(parent = baseenv())
    sys.source(file.path(root, files[[1L]]), envir = driver)
    # Namespace/path doubles are unit-test-only; installed service proof is separate.
    driver$loadNamespace <- function(...) asNamespace("MultiScholaR")
    driver$getNamespaceInfo <- function(...) file.path(library_path, "MultiScholaR")
    driver$recordedCheckpointApi <- function(...) NULL
    withr::local_libpaths(.libPaths())
    helpers <- new.env(parent = baseenv())
    sys.source(file.path(root, files[[4L]]), envir = helpers)
    helpers$source_root <- root
    first <- makeWorkflowCheckpointProtein(recipe)
    second <- driver$checkpointUtilityProteinState(helpers, recipe, 1L,
        "shared_metadata", first, root)
    prepared <- prepareDIAWorkflowHistory(list(second = list(data = second)))$second$data
    expect_identical(rlang::obj_address(first@protein_id_table),
        rlang::obj_address(prepared@protein_id_table))
    expect_identical(rlang::obj_address(first@design_matrix),
        rlang::obj_address(prepared@design_matrix))
    expect_identical(rlang::obj_address(first@protein_quant_table[[3L]]),
        rlang::obj_address(prepared@protein_quant_table[[3L]]))
    expect_false(identical(rlang::obj_address(first@protein_quant_table[[2L]]),
        rlang::obj_address(prepared@protein_quant_table[[2L]])))
    dir.create(file.path(library_path, "MultiScholaR", "R"))
    for (path in c("NAMESPACE", "R/MultiScholaR.rdb", "R/MultiScholaR.rdx")) {
        writeLines("bound unit fixture", file.path(library_path, "MultiScholaR", path))
    }
    input_dir <- file.path(root, "prepared")
    dir.create(input_dir)
    manifest <- driver$prepareCheckpointUtilityInputs(library_path, input_dir, contract, "small")
    expect_identical(manifest$status, "prepared_utility_inputs")
    expect_length(manifest$inputs, (recipe$history_depth + 2L) * 10L + 1L)
    expect_identical(driver$readCheckpointUtilityInput(input_dir, manifest$inputs, 0L), first)
    driver$sys.source <- function(file, envir, ...) {
        base::sys.source(file, envir, ...)
        envir$makeWorkflowCheckpointProtein <- function(...) {
            stop("fixture construction reached timed member")
        }
    }
    hashes <- list()
    for (treatment in c("resident", "saved", "qs2")) {
        member_dir <- file.path(root, treatment)
        verify_dir <- file.path(root, paste0(treatment, "-verify"))
        dir.create(member_dir)
        dir.create(verify_dir)
        result <- driver$executeCheckpointUtilityMember(library_path, member_dir,
            contract, "small", treatment, "shared_metadata", 1L, input_dir,
            diagnostics = identical(treatment, "qs2"),
            memory_snapshot = \() list(current_bytes = 100, peak_bytes = 120,
                memory_stat = list(status = "available", values = list(anon = 60, file = 30))))
        expect_identical(result$diagnostics, identical(treatment, "qs2"))
        if (isTRUE(result$diagnostics)) {
            phases <- lapply(readLines(file.path(member_dir, "memory-phases.jsonl")),
                jsonlite::fromJSON)
            labels <- vapply(phases, \(item) item$label, character(1))
            expect_true(all(c("before_namespace", "forward_after_gc", "after_reverts",
                "before:complete_legacy_export_serialization", "after_export_release",
                "complete_before_teardown") %in% labels))
            expect_false(any(grepl("failure", labels)))
            expect_lt(match("after_export_state_hashes", labels),
                match("export_after_gc", labels))
            expect_lt(match("export_after_gc", labels),
                match("before:complete_legacy_export_serialization", labels))
            expect_gt(file.size(file.path(member_dir, "utility.Rprof")), 0)
            expect_true(file.exists(file.path(member_dir, "utility.Rprofmem")))
        } else {
            expect_false(file.exists(file.path(member_dir, "memory-phases.jsonl")))
        }
        expect_identical(result$cycles_completed, 10L)
        expect_length(result$cycle_metrics, 10L)
        expect_true(all(vapply(result$cycle_metrics,
            \(item) identical(item$reset_state, "state_001"), logical(1))))
        expect_identical(result$current_state, "branch_001")
        expect_identical(result$active_history, c("state_001", "state_002", "branch_001"))
        expect_length(result$state_names, recipe$history_depth + 1L)
        expect_false(result$fixture_generation_in_timed_service)
        expect_false(result$performance_qualified)
        verified <- driver$verifyCheckpointUtilityMember(library_path, verify_dir, member_dir)
        expect_identical(verified$states_checked, recipe$history_depth + 1L)
        expect_identical(verified$persistent_states_checked,
            if (treatment == "resident") 0L else recipe$history_depth + 1L)
        hashes[[treatment]] <- result$state_hashes
        if (treatment != "resident") {
            refs <- readRDS(file.path(member_dir, "references.rds"))
            expect_true(all(vapply(refs, \(ref) identical(ref$payload$codec,
                if (treatment == "qs2") "qs2" else "rds"), logical(1))))
        }
    }
    expect_identical(hashes$resident, hashes$saved)
    expect_identical(hashes$resident, hashes$qs2)
    bad <- manifest$inputs
    bad[[1L]]$bytes <- 1
    expect_error(driver$readCheckpointUtilityInput(input_dir, bad, 0L), "integrity differs")
    bad[[1L]]$path <- "../different.rds"
    expect_error(driver$readCheckpointUtilityInput(input_dir, bad, 0L), "generation binding differs")
    driver$sys.source <- base::sys.source
    failed_store <- NULL
    driver$checkpointUtilityAcceptState <- function(context, ...) {
        failed_store <<- context$store
        stop("injected serialization failure")
    }
    failed_dir <- file.path(root, "failed")
    dir.create(failed_dir)
    expect_error(driver$executeCheckpointUtilityMember(library_path, failed_dir,
        contract, "small", "qs2", "independent", 1L, diagnostics = TRUE,
        memory_snapshot = \() list(current_bytes = 100)), "injected serialization failure")
    expect_true(failed_store$closed)
    expect_identical(failed_store$cache_bytes, 0)
    expect_false(file.exists(file.path(failed_dir, "legacy-session.rds")))
    failed_phases <- lapply(readLines(file.path(failed_dir, "memory-phases.jsonl")),
        jsonlite::fromJSON)
    expect_identical(tail(failed_phases, 1L)[[1L]]$label, "failure_before_teardown")
})

test_that("utility service handoff binds accounting and separate verification", {
    environment <- checkpointHarnessTestEnvironment()
    root <- withr::local_tempdir()
    library_path <- file.path(root, "library")
    package <- file.path(library_path, "MultiScholaR")
    dir.create(file.path(package, "Meta"), recursive = TRUE)
    dir.create(file.path(package, "R"))
    writeLines(c("Package: MultiScholaR", "Version: 0.5.0",
        "Built: R 4.6.1; ; 2026-09-16; unix", "Imports: digest"),
        file.path(package, "DESCRIPTION"))
    for (path in c("NAMESPACE", "Meta/package.rds", "R/MultiScholaR.rdb",
        "R/MultiScholaR.rdx")) writeLines("fixture", file.path(package, path))
    member <- file.path(root, "browser-smoke.ABC123")
    verifier <- file.path(root, "browser-smoke.DEF456")
    dir.create(member)
    dir.create(verifier)
    terminal <- list(status = "service_tree_interval_only", service_succeeded = TRUE,
        invocation_id = strrep("a", 32), elapsed_seconds = 2,
        raw_properties = list(CPUUsageNSec = "1000000000", MemoryPeak = "1048576",
            MemorySwapPeak = "0", IOReadBytes = "1024", IOWriteBytes = "2048"))
    receipt <- list(status = "complete", science_ok = TRUE, outputs_ok = TRUE,
        recipe = "small", treatment = "saved", variant = "independent", block = 1L,
        process_id = 42L, retained_heap_bytes = 1000, max_save_seconds = 0.1,
        max_cold_revert_seconds = 0.2, payload_read_bytes = 512,
        export_heap_bytes = 2000)
    verified <- list(status = "passed_utility_member", member_dir = member)
    jsonlite::write_json(terminal, file.path(member, "terminal-interval.json"),
        auto_unbox = TRUE)
    jsonlite::write_json(receipt, file.path(member, "member-receipt.json"),
        auto_unbox = TRUE)
    jsonlite::write_json(verified, file.path(verifier, "utility-verification.json"),
        auto_unbox = TRUE)
    operations <- character()
    service <- function(operation, library_path, output_root, arguments) {
        operations <<- c(operations, operation)
        if (operation == "utility-member") {
            expect_identical(arguments[-1L], c("small", "saved", "independent", "1"))
            return(list(exit_status = 0L, member_dir = member,
                unit = "multischolar-dia003-browser-smoke.ABC123.service"))
        }
        expect_identical(arguments, member)
        list(exit_status = 0L, member_dir = verifier,
            unit = "multischolar-dia003-browser-smoke.DEF456.service")
    }
    contract <- testthat::test_path("..", "testdata", "workflow-checkpoints",
        "scenarios.json")
    result <- environment$runCheckpointUtilityServiceMember(list(id = "small"),
        1L, "saved", "independent", library_path, root, contract, service)
    expect_identical(operations, c("utility-member", "utility-verify"))
    expect_identical(result$status, "complete")
    expect_identical(result$member, "candidate")
    expect_identical(result$combined_peak_bytes, 1048576)
    expect_identical(result$swap_peak_bytes, 0)
    expect_identical(result$retained_heap_bytes, 1000L)
    expect_false(result$promotion_authorized)
    receipt$diagnostics <- TRUE
    jsonlite::write_json(receipt, file.path(member, "member-receipt.json"), auto_unbox = TRUE)
    expect_error(environment$runCheckpointUtilityServiceMember(list(id = "small"),
        1L, "saved", "independent", library_path, root, contract, service),
        "diagnostic utility member cannot qualify performance")
})

test_that("service invoker binds terminal markers without guessing the newest directory", {
    environment <- checkpointHarnessTestEnvironment()
    root <- withr::local_tempdir()
    member <- file.path(root, "browser-smoke.ABC123")
    dir.create(member)
    environment$checkpointRunnerDirectory <- function() root
    lines <- c('{"unit":"multischolar-dia003-browser-smoke.ABC123.service"}',
        "unit=multischolar-dia003-browser-smoke.ABC123")
    environment$system2 <- function(...) lines
    result <- environment$runCheckpointService("browser-da-default-qc", "library", root, c("r", "f", "d"))
    expect_identical(result$member_dir, normalizePath(member))
    expect_identical(result$exit_status, 0L)
    attr(lines, "status") <- 75L
    result <- environment$runCheckpointService("browser-da-default-qc", "library", root, c("r", "f", "d"))
    expect_identical(result$exit_status, 75L)
    lines <- c(lines, "unit=multischolar-dia003-browser-smoke.DEF456")
    expect_error(environment$runCheckpointService("browser-da-default-qc", "library", root, c("r", "f", "d")),
        "no unique owned")
})

test_that("allocation diagnostic guards refuse low headroom and stop only the owned service", {
    skip_if_not(identical(Sys.info()[["sysname"]], "Linux"))
    root <- withr::local_tempdir()
    bin <- file.path(root, "bin")
    library_path <- file.path(root, "library")
    dir.create(bin)
    dir.create(library_path)
    report <- file.path(root, "report.tsv")
    writeLines("tiny orchestration fixture, never scientific data", report)
    stop_log <- file.path(root, "stopped")
    launch_log <- file.path(root, "launched")
    writeLines(c("#!/usr/bin/env bash", "exit 0"), file.path(bin, "Rscript"))
    writeLines(c("#!/usr/bin/env bash", 'printf "%s\\n" "$@" >> "$DIA_GUARD_LAUNCH"'),
        file.path(bin, "systemd-run"))
    writeLines(c("#!/usr/bin/env bash", 'echo "$DIA_GUARD_AVAILABLE"'), file.path(bin, "awk"))
    writeLines(c("#!/usr/bin/env bash",
        'if [[ "$2" == show ]]; then',
        'printf "SubState=%s\\nMemoryCurrent=%s\\nMemoryPeak=%s\\nResult=success\\nExecMainStatus=0\\n" "$DIA_GUARD_STATE" "$DIA_GUARD_CURRENT" "${DIA_GUARD_PEAK:-0}"',
        'elif [[ "$2" == stop ]]; then printf "%s\\n" "$3" >> "$DIA_GUARD_STOP"; else exit 3; fi'),
        file.path(bin, "systemctl"))
    Sys.chmod(file.path(bin, c("Rscript", "systemd-run", "systemctl", "awk")), "0755")
    withr::local_envvar(c(PATH = paste(bin, Sys.getenv("PATH"), sep = .Platform$path.sep),
        CHROMOTE_CHROME = Sys.which("true"), DIA_GUARD_STOP = stop_log,
        DIA_GUARD_LAUNCH = launch_log, DIA_GUARD_AVAILABLE = "8589934592",
        DIA_GUARD_CURRENT = "1610612736", DIA_GUARD_STATE = "running"))
    runner <- normalizePath(testthat::test_path("..", "..", "tools", "profiling",
        "run_workflow_checkpoint_proof.sh"))
    run <- function() processx::run("bash", c(runner, "--service", "exclusion-allocation",
        library_path, root, report), error_on_status = FALSE, timeout = 5000)
    for (available in c("8589934591", "", "unknown")) {
        result <- withr::with_envvar(c(DIA_GUARD_AVAILABLE = available), run())
        expect_identical(result$status, 1L)
        expect_false(file.exists(launch_log))
    }
    for (current in c("1610612736", "2147483648", "[not set]")) {
        before <- list.dirs(root, recursive = FALSE)
        result <- withr::with_envvar(c(DIA_GUARD_CURRENT = current), run())
        expect_identical(result$status, 76L)
        directory <- setdiff(list.dirs(root, recursive = FALSE), before)
        expect_length(directory, 1L)
        expect_true(file.exists(file.path(directory, "guard-stop.txt")))
        request <- checkpointHarnessTestEnvironment()$readCheckpointControllerRequest(
            file.path(directory, "controller-request.nul"))
        expect_identical(request$memory_max_bytes, 2147483648)
        expect_identical(request$runtime_max_seconds, 120)
        expect_identical(tail(readLines(stop_log), 1L), request$unit)
    }
    expect_true(all(c("--property=MemoryMax=2147483648", "--property=MemorySwapMax=0",
        "--property=RuntimeMaxSec=120") %in% readLines(launch_log)))
    stops_before <- length(readLines(stop_log))
    result <- withr::with_envvar(c(DIA_GUARD_STATE = "exited", DIA_GUARD_CURRENT = "0"), run())
    expect_identical(result$status, 0L)
    expect_identical(length(readLines(stop_log)), stops_before + 1L)
    recorded <- function() processx::run("bash", c(runner, "--service", "recorded-prepare",
        library_path, root, report), error_on_status = FALSE, timeout = 5000)
    for (state in c("running", "exited")) {
        before <- list.dirs(root, recursive = FALSE)
        result <- withr::with_envvar(c(DIA_GUARD_AVAILABLE = "10737418240",
            DIA_GUARD_CURRENT = "0", DIA_GUARD_PEAK = "8053063680",
            DIA_GUARD_STATE = state), recorded())
        expect_identical(result$status, 76L)
        directory <- setdiff(list.dirs(root, recursive = FALSE), before)
        expect_length(directory, 1L)
        expect_true(file.exists(file.path(directory, "guard-stop.txt")))
        request <- checkpointHarnessTestEnvironment()$readCheckpointControllerRequest(
            file.path(directory, "controller-request.nul"))
        expect_identical(request$memory_max_bytes, 8589934592)
        expect_identical(request$runtime_max_seconds, 900)
        expect_identical(tail(readLines(stop_log), 1L), request$unit)
    }
    # Exercise exact safe/unsafe boundaries without allocating memory.
    lines <- readLines(runner)
    guard <- paste(lines[seq_len(grep("^proof_script_dir=", lines) - 1L)], collapse = "\n")
    for (current in c("0", "1610612735")) {
        expect_identical(processx::run("bash", c("-c", paste(guard,
            paste("checkpointDiagnosticGuard 8589934592", current), sep = "\n")))$status, 0L)
    }
})

test_that("allocation profiling dispatch is separate from matched-cache and acceptance cases", {
    environment <- checkpointHarnessTestEnvironment()
    environment$runCheckpointExclusionBenchmark <- function(...) list(...)
    result <- environment$checkpointProofMain(c("--internal-exclusion-allocation", "library", "run", "report"))
    expect_true(result$allocation_diagnostic)
    expect_false(result$memory_diagnostic)
    qvalue <- environment$checkpointProofMain(c("--internal-exclusion-allocation", "library", "run", "report", "qvalue"))
    expect_true(qvalue$allocation_diagnostic)
    expect_true(qvalue$qvalue_diagnostic)
})

test_that("allocation profiler captures a tiny allocation and stops before later work", {
    skip_if_not(isTRUE(capabilities("profmem")))
    path <- withr::local_tempfile()
    utils::Rprofmem(path, threshold = 100000)
    withr::defer(utils::Rprofmem(NULL))
    allocateFixture <- function() numeric(20000L)
    retained <- allocateFixture()
    utils::Rprofmem(NULL)
    records <- readLines(path)
    expect_true(any(grepl("allocateFixture", records, fixed = TRUE)))
    expect_length(retained, 20000L)
    after <- numeric(30000L)
    expect_identical(readLines(path), records)
})

test_that("service memory reader validates the diagnostic cap without relaxing baseline limits", {
    environment <- checkpointHarnessTestEnvironment()
    actual_cap <- "2147483648"
    environment$readLines <- function(path, ...) {
        if (identical(path, "/proc/self/cgroup")) return("0::/owned-fixture")
        if (basename(path) == "memory.max") return(actual_cap)
        "0"
    }
    environment$readCheckpointCounters <- function(...) list(status = "fixture")
    expect_error(environment$checkpointServiceMemory(), "required service limits")
    expect_identical(environment$checkpointServiceMemory("2147483648")$peak_bytes, 0)
    actual_cap <- "8589934592"
    expect_error(environment$checkpointServiceMemory("2147483648"), "required service limits")
    expect_identical(environment$checkpointServiceMemory()$peak_bytes, 0)
    actual_cap <- "5368709120"
    expect_identical(environment$checkpointServiceMemory("5368709120")$peak_bytes, 0)
    expect_error(environment$checkpointServiceMemory(), "required service limits")
    actual_cap <- "6442450944"
    expect_identical(environment$checkpointServiceMemory("6442450944")$peak_bytes, 0)
    expect_error(environment$checkpointServiceMemory(), "required service limits")
    actual_cap <- "11811160064"
    expect_identical(environment$checkpointServiceMemory("11811160064")$peak_bytes, 0)
    expect_error(environment$checkpointServiceMemory(), "required service limits")
    actual_cap <- "17179869184"
    expect_identical(environment$checkpointServiceMemory("17179869184")$peak_bytes, 0)
    expect_error(environment$checkpointServiceMemory(), "required service limits")
    expect_error(environment$checkpointServiceMemory("max"), "unsupported")
    root <- withr::local_tempdir()
    environment$inspectCheckpointLibrary <- function(...) list(package_path = root)
    observed_caps <- character()
    environment$checkpointServiceMemory <- function(expected_memory_max = "8589934592") {
        observed_caps <<- c(observed_caps, expected_memory_max)
        list()
    }
    environment$loadNamespace <- function(...) stop("stop before loading scientific code")
    expect_error(environment$runCheckpointExclusionBenchmark("unused", root, "unused",
        allocation_diagnostic = TRUE), "stop before loading")
    expect_identical(observed_caps, rep("2147483648", 2L))
})

test_that("memory diagnostic observes cache residency and refuses an unmatched entry", {
    skip_if_not(identical(Sys.info()[["sysname"]], "Linux"))
    skip_if_not(nzchar(Sys.which("fincore")))
    environment <- checkpointHarnessTestEnvironment()
    report <- withr::local_tempfile()
    writeLines("cache observation fixture", report)
    observation <- environment$checkpointInputCache(report)
    expect_true(observation$fully_resident)
    expect_identical(observation$size_bytes, as.numeric(file.info(report)$size))
    expect_false(observation$content_read_or_warmed)
    root <- withr::local_tempdir()
    environment$inspectCheckpointLibrary <- function(...) list(package_path = root)
    environment$checkpointServiceMemory <- function() list()
    environment$checkpointInputCache <- function(...) list(fully_resident = FALSE)
    environment$loadNamespace <- function(...) stop("scientific load must not occur")
    expect_error(environment$runCheckpointExclusionBenchmark("unused", root, report,
        memory_diagnostic = TRUE), "cache condition unmatched")
    receipt <- jsonlite::read_json(file.path(root, "exclusion-receipt.json"))
    expect_identical(receipt$status, "failed")
    expect_true(receipt$diagnostics)
    expect_false(file.exists(file.path(root, "classification.rds")))
})

test_that("profile accounting excludes multiline label continuations from sample counts", {
    environment <- checkpointHarnessTestEnvironment()
    path <- withr::local_tempfile()
    writeLines(c("memory profiling: GC profiling: sample.interval=20000",
        ':1:2:3:4:"sub" "classify" "bindEvent(observe({',
        '  label <- "classify"', '  }))" "serviceApp"',
        ':1:2:3:4:"other" "serviceApp"'), path)
    result <- environment$summarizeCheckpointProfileFrames(path, c("classify", "missing"))
    expect_identical(result$sample_records, 2L)
    expect_identical(result$continuation_lines, 2L)
    expect_identical(result$sampled_seconds, 0.04)
    expect_identical(result$frames$inclusive_samples_lower_bound, c(1L, 0L))
    expect_identical(result$frames$sample_share_lower_bound, c(0.5, 0))
})

test_that("scientific waits share the remaining service budget and controls stay bounded", {
    environment <- new.env(parent = baseenv())
    sys.source(testthat::test_path("..", "..", "tools", "profiling",
        "workflow_checkpoint_dia_driver.R"), envir = environment)
    environment$proc.time <- function() c(elapsed = 0)
    path <- withr::local_tempfile()
    utils::write.table(data.frame(sample = c("s1", "s2"), group = c("control", "treated")),
        path, sep = "\t", row.names = FALSE, quote = FALSE)
    environment$checkpointSetInput <- function(...) NULL
    environment$checkpointChooseExisting <- function(driver, id, values, ...) values
    environment$waitCheckpointStepComplete <- function(...) TRUE
    waits <- integer()
    controls <- integer()
    environment$proc.time <- function() c(elapsed = 120)
    driver <- list(wait_for_js = function(script, timeout) controls <<- c(controls, timeout), get_js = function(...) TRUE,
        wait_for_idle = function(timeout) waits <<- c(waits, timeout))
    environment$runCheckpointDesignActions(driver, path)
    expect_true(all(waits == 780000L))
    expect_true(all(controls == 30000L))
    expect_identical(environment$checkpointWaitBudget(elapsed_seconds = 899), 1000L)
    expect_identical(environment$checkpointWaitBudget(control = TRUE, elapsed_seconds = 899), 1000L)
    expect_error(environment$checkpointWaitBudget(elapsed_seconds = 900), "exhausted")
    expect_error(environment$checkpointWaitBudget(elapsed_seconds = Inf), "invalid")
    expect_error(environment$checkpointWaitBudget(control = NA), "invalid")
})

test_that("default-QC dispatch cannot silently reuse relaxed tiny-fixture settings", {
    environment <- checkpointHarnessTestEnvironment()
    environment$runCheckpointBrowserSmoke <- function(library_path, run_dir, inputs) inputs
    paths <- c("library", "run", "report", "fasta", "design")
    defaults <- environment$checkpointProofMain(c("--internal-browser-da-default-qc", paths))
    tiny <- environment$checkpointProofMain(c("--internal-browser-da", paths))
    expect_true(defaults$qc_defaults)
    expect_false(tiny$qc_defaults)
    expect_false(tiny$profile_application)
    expect_false(defaults$profile_application)
    profile <- environment$checkpointProofMain(c("--internal-browser-da-profile", paths))
    expect_true(profile$qc_defaults)
    expect_true(profile$profile_application)
    baseline <- environment$checkpointProofMain(c("--internal-browser-baseline", paths))
    expect_true(baseline$qc_defaults)
    expect_true(baseline$summary_diagnostic)
    expect_false(baseline$profile_application)
    defaults$qc_defaults <- tiny$qc_defaults
    expect_identical(defaults, tiny)
})

test_that("compact companions match the frozen report generator without precursor materialization", {
    environment <- checkpointHarnessTestEnvironment()
    helpers <- new.env(parent = baseenv())
    source_root <- normalizePath(testthat::test_path("..", ".."))
    sys.source(file.path(source_root, "tests", "testthat", "helper-workflow-checkpoint-fixtures.R"),
        envir = helpers)
    recipe <- helpers$readWorkflowCheckpointRecipe("small", source_root)
    companions <- environment$checkpointInputCompanions(recipe, helpers)
    report <- helpers$makeWorkflowPrecursorChunk(recipe, seq_len(recipe$n_proteins))
    expected_sequences <- vapply(unique(report$Protein.Ids), \(id) {
        paste(unique(report$Stripped.Sequence[report$Protein.Ids == id]), collapse = "")
    }, character(1))
    expect_identical(companions$fasta[seq(2L, length(companions$fasta), by = 2L)], unname(expected_sequences))
    expect_identical(companions$design$Run, unique(report$Run))
    expect_identical(companions$design$sample, companions$design$Run)
    expect_identical(levels(companions$design$batch), c("batch_1", "batch_2", "unused_batch"))
    helpers$makeWorkflowPrecursorChunk <- function(...) stop("must not read/build precursor rows")
    helpers$makeWorkflowCheckpointPeptide <- function(...) stop("must not build S4")
    large <- helpers$readWorkflowCheckpointRecipe("representative", source_root)
    large_companions <- environment$checkpointInputCompanions(large, helpers)
    expect_length(large_companions$fasta, 2L * large$n_proteins)
    expect_equal(nrow(large_companions$design), large$n_samples)
})

test_that("retained interval binds monotonic timestamps to a terminal main process", {
    environment <- checkpointHarnessTestEnvironment()
    path <- withr::local_tempfile()
    fields <- c("Id=owned.service", paste0("InvocationID=", strrep("a", 32)),
        "SubState=exited", "MainPID=0", "ExecMainStartTimestampMonotonic=1976097228",
        "ExecMainExitTimestampMonotonic=1976666656", "TasksCurrent=[not set]")
    writeLines(fields, path)
    result <- environment$readCheckpointRetainedInterval(path, "owned.service")
    expect_identical(result$main_elapsed_seconds, 0.569428)
    expect_null(result$elapsed_seconds)
    expect_false(result$promotion_authorized)
    expect_identical(result$raw_properties$TasksCurrent, "[not set]")
    expect_identical(result$service_result, "unknown")
    expect_false(result$service_succeeded)
    expect_error(environment$readCheckpointRetainedInterval(path, "another.service"), "terminal")
    writeLines(c(fields, "Id=owned.service"), path)
    expect_error(environment$readCheckpointRetainedInterval(path, "owned.service"), "ambiguous")
    writeLines(sub("MainPID=0", "MainPID=123", fields, fixed = TRUE), path)
    expect_error(environment$readCheckpointRetainedInterval(path, "owned.service"), "terminal")
    writeLines(sub("1976666656", "0", fields, fixed = TRUE), path)
    expect_error(environment$readCheckpointRetainedInterval(path, "owned.service"), "interval")
    tree_fields <- c(fields, "ExitType=cgroup", "ControlGroup=", "StateChangeTimestampMonotonic=1977666656")
    writeLines(tree_fields, path)
    tree <- environment$readCheckpointRetainedInterval(path, "owned.service")
    expect_identical(tree$elapsed_seconds, 1.569428)
    expect_gt(tree$elapsed_seconds, tree$main_elapsed_seconds)
    writeLines(sub("ExitType=cgroup", "ExitType=main", tree_fields, fixed = TRUE), path)
    expect_null(environment$readCheckpointRetainedInterval(path, "owned.service")$elapsed_seconds)
    writeLines(sub("ControlGroup=", "ControlGroup=/still-populated", tree_fields, fixed = TRUE), path)
    expect_null(environment$readCheckpointRetainedInterval(path, "owned.service")$elapsed_seconds)
    root <- withr::local_tempdir()
    writeLines(tree_fields, file.path(root, "terminal-properties.txt"))
    arguments <- c("--internal-retained-interval", root, "owned.service")
    environment$checkpointProofMain(arguments)
    persisted <- jsonlite::read_json(file.path(root, "terminal-interval.json"))
    expect_identical(persisted$elapsed_seconds, 1.569428)
    expect_error(environment$checkpointProofMain(arguments), "overwrite")
    oom_fields <- c(sub("SubState=exited", "SubState=failed", tree_fields, fixed = TRUE),
        "Result=oom-kill", "ExecMainCode=2", "ExecMainStatus=15")
    writeLines(oom_fields, path)
    oom <- environment$readCheckpointRetainedInterval(path, "owned.service")
    expect_identical(oom$service_result, "oom-kill")
    expect_false(oom$service_succeeded)
    expect_null(oom$elapsed_seconds)
})

test_that("terminal accounting uses exact unit-bound counters and preserves missing metrics", {
    environment <- checkpointHarnessTestEnvironment()
    path <- withr::local_tempfile()
    record <- list(USER_UNIT = "owned.service", USER_INVOCATION_ID = strrep("a", 32),
        MESSAGE_ID = "ae8f7b866b0347b9af31fe1c80b127c0", CPU_USAGE_NSEC = "476526000",
        MEMORY_PEAK = "397221888", MEMORY_SWAP_PEAK = "0",
        IO_METRIC_READ_BYTES = "18446744073709551615")
    jsonlite::write_json(record, path, auto_unbox = TRUE)
    result <- environment$readCheckpointTerminalAccounting(path, "owned.service")
    expect_identical(result$memory_peak_bytes, 397221888)
    expect_identical(result$cpu_usage_nsec, 476526000)
    expect_identical(result$swap_peak_bytes, 0)
    expect_null(result$io_read_bytes)
    expect_null(result$io_write_bytes)
    expect_null(result$elapsed_seconds)
    expect_false(result$promotion_authorized)
    expect_error(environment$readCheckpointTerminalAccounting(path, "different.service"), "missing")
    line <- readLines(path)
    writeLines(rep(line, 2), path)
    expect_error(environment$readCheckpointTerminalAccounting(path, "owned.service"), "ambiguous")
    record$USER_INVOCATION_ID <- NULL
    jsonlite::write_json(record, path, auto_unbox = TRUE)
    expect_error(environment$readCheckpointTerminalAccounting(path, "owned.service"), "invocation")
})

test_that("accounted browser dispatch retains terminal evidence and stops its owned unit", {
    skip_if_not(identical(Sys.info()[["sysname"]], "Linux"))
    root <- withr::local_tempdir()
    bin <- file.path(root, "bin")
    library_path <- file.path(root, "library")
    dir.create(bin)
    dir.create(library_path)
    inputs <- file.path(root, c("report.tsv", "proteins.fasta", "design.tsv"))
    for (path in inputs) writeLines("orchestration fixture only", path)
    stop_log <- file.path(root, "stopped-unit")
    writeLines(c("#!/usr/bin/env bash", "exit 0"), file.path(bin, "Rscript"))
    writeLines(c("#!/usr/bin/env bash", "exit 0"), file.path(bin, "systemd-run"))
    writeLines(c("#!/usr/bin/env bash", "echo 21474836480"), file.path(bin, "awk"))
    writeLines(c("#!/usr/bin/env bash",
        'if [[ "$2" == show ]]; then',
        'if [[ "${DIA_STUB_GUARD_TRIP:-0}" == 1 ]]; then',
        'printf "Id=%s\\nInvocationID=aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\\nSubState=running\\nMemoryCurrent=12000000000\\nMemoryPeak=16642998272\\nResult=success\\nExecMainStatus=0\\n" "$3"',
        'exit 0; fi',
        'printf "Id=%s\\nInvocationID=aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\\nSubState=exited\\nMainPID=0\\nMemoryCurrent=0\\nMemoryPeak=1048576\\nResult=%s\\nExecMainStatus=%s\\n" "$3" "${DIA_STUB_RESULT:-success}" "${DIA_STUB_STATUS:-0}"',
        'elif [[ "$2" == stop ]]; then printf "%s\\n" "$3" >> "$DIA_STUB_STOP_LOG"; else exit 3; fi'),
        file.path(bin, "systemctl"))
    Sys.chmod(file.path(bin, c("Rscript", "systemd-run", "systemctl", "awk")), "0755")
    withr::local_envvar(c(PATH = paste(bin, Sys.getenv("PATH"), sep = .Platform$path.sep),
        CHROMOTE_CHROME = Sys.which("true"), DIA_STUB_STOP_LOG = stop_log))
    runner <- normalizePath(testthat::test_path("..", "..", "tools", "profiling",
        "run_workflow_checkpoint_proof.sh"))
    for (failed in c(FALSE, TRUE)) {
        before <- list.dirs(root, recursive = FALSE)
        result <- withr::with_envvar(c(DIA_STUB_RESULT = if (failed) "exit-code" else "success",
            DIA_STUB_STATUS = if (failed) "7" else "0"), processx::run("bash",
            c(runner, "--service", "browser-da-accounted", library_path, root, inputs), error_on_status = FALSE))
        directory <- setdiff(list.dirs(root, recursive = FALSE), before)
        expect_length(directory, 1L)
        expect_identical(result$status, if (failed) 1L else 0L)
        request <- checkpointHarnessTestEnvironment()$readCheckpointControllerRequest(
            file.path(directory, "controller-request.nul"))
        expect_identical(request$operation, "browser-da-accounted")
        expect_true("--internal-browser-da" %in% request$command)
        expect_identical(tail(request$command, 3L), inputs)
        expect_identical(tail(readLines(stop_log), 1L), request$unit)
        expect_true(file.exists(file.path(directory, "terminal-properties.txt")))
    }
    before <- list.dirs(root, recursive = FALSE)
    result <- processx::run("bash", c(runner, "--service", "browser-baseline",
        library_path, root, inputs), error_on_status = FALSE)
    directory <- setdiff(list.dirs(root, recursive = FALSE), before)
    expect_length(directory, 1L)
    expect_identical(result$status, 0L)
    request <- checkpointHarnessTestEnvironment()$readCheckpointControllerRequest(
        file.path(directory, "controller-request.nul"))
    expect_identical(request$operation, "browser-baseline")
    expect_true("--internal-browser-baseline" %in% request$command)
    expect_identical(request$memory_max_bytes, 8589934592)
    before <- list.dirs(root, recursive = FALSE)
    contract <- normalizePath(testthat::test_path("..", "testdata",
        "workflow-checkpoints", "scenarios.json"))
    result <- processx::run("bash", c(runner, "--service", "utility-member",
        library_path, root, contract, "small", "resident", "independent", "1"),
        error_on_status = FALSE)
    directory <- setdiff(list.dirs(root, recursive = FALSE), before)
    expect_length(directory, 1L)
    expect_identical(result$status, 0L)
    request <- checkpointHarnessTestEnvironment()$readCheckpointControllerRequest(
        file.path(directory, "controller-request.nul"))
    expect_identical(request$operation, "utility-member")
    expect_true("--internal-utility-member" %in% request$command)
    expect_identical(tail(request$command, 4L), c("small", "resident", "independent", "1"))
    expect_identical(request$memory_max_bytes, 8589934592)
    before <- list.dirs(root, recursive = FALSE)
    result <- withr::with_envvar(c(DIA_STUB_GUARD_TRIP = "1"), processx::run("bash",
        c(runner, "--service", "utility-member", library_path, root, contract,
            "small", "qs2", "independent", "1"), error_on_status = FALSE))
    directory <- setdiff(list.dirs(root, recursive = FALSE), before)
    expect_identical(result$status, 76L)
    request <- checkpointHarnessTestEnvironment()$readCheckpointControllerRequest(
        file.path(directory, "controller-request.nul"))
    expect_identical(tail(request$command, 4L), c("small", "qs2", "independent", "1"))
    expect_identical(tail(readLines(stop_log), 1L), request$unit)
    expect_true(file.exists(file.path(directory, "guard-stop.txt")))
    before <- list.dirs(root, recursive = FALSE)
    result <- processx::run("bash", c(runner, "--service", "browser-baseline",
        library_path, root, inputs, "--capped-sixteen-gib"), error_on_status = FALSE)
    directory <- setdiff(list.dirs(root, recursive = FALSE), before)
    expect_length(directory, 1L)
    expect_identical(result$status, 0L)
    request <- checkpointHarnessTestEnvironment()$readCheckpointControllerRequest(
        file.path(directory, "controller-request.nul"))
    expect_identical(request$operation, "browser-baseline")
    expect_identical(request$memory_max_bytes, 17179869184)
    expect_identical(request$swap_max_bytes, 0)
    expect_true("--internal-browser-baseline" %in% request$command)
    before <- list.dirs(root, recursive = FALSE)
    result <- withr::with_envvar(c(DIA_STUB_GUARD_TRIP = "1"), processx::run("bash",
        c(runner, "--service", "browser-baseline", library_path, root, inputs,
            "--capped-sixteen-gib"), error_on_status = FALSE))
    directory <- setdiff(list.dirs(root, recursive = FALSE), before)
    expect_length(directory, 1L)
    expect_identical(result$status, 76L)
    expect_match(paste(readLines(file.path(directory, "guard-stop.txt")), collapse = "\n"),
        "reason=early_memory_stop")
    request <- checkpointHarnessTestEnvironment()$readCheckpointControllerRequest(
        file.path(directory, "controller-request.nul"))
    expect_identical(tail(readLines(stop_log), 1L), request$unit)
})

test_that("controller preserves service and log failures in durable terminal records", {
    skip_if_not(identical(Sys.info()[["sysname"]], "Linux"))
    root <- withr::local_tempdir()
    bin <- file.path(root, "bin")
    library_path <- file.path(root, "library with spaces\nand 'quotes'")
    dir.create(bin)
    dir.create(library_path)
    # Shell orchestration doubles only: no scientific/service qualification.
    writeLines(c("#!/usr/bin/env bash", "exit 0"), file.path(bin, "Rscript"))
    writeLines(c("#!/usr/bin/env bash", "echo retained-service-output",
        'if [[ -n "${DIA_STUB_READY:-}" ]]; then touch -- "$DIA_STUB_READY"; sleep 2; fi',
        'exit "${DIA_STUB_SERVICE_EXIT:-0}"'), file.path(bin, "systemd-run"))
    real_tee <- Sys.which("tee")
    writeLines(c("#!/usr/bin/env bash", paste(shQuote(real_tee), '"$@"'),
        'if [[ "$1" == */service.log && "${DIA_STUB_LOG_FAIL:-0}" == 1 ]]; then exit 1; fi'),
        file.path(bin, "tee"))
    Sys.chmod(file.path(bin, c("Rscript", "systemd-run", "tee")), "0755")
    withr::local_envvar(c(PATH = paste(bin, Sys.getenv("PATH"), sep = .Platform$path.sep),
        CHROMOTE_CHROME = Sys.which("true"), DIA_STUB_SERVICE_EXIT = "0",
        DIA_STUB_LOG_FAIL = "0", DIA_STUB_READY = NA_character_))
    runner <- normalizePath(testthat::test_path("..", "..", "tools", "profiling",
        "run_workflow_checkpoint_proof.sh"))
    for (case in list(c(service = 0, log = 0, expected = 0),
        c(service = 7, log = 0, expected = 7), c(service = 0, log = 1, expected = 74))) {
        before <- list.dirs(root, recursive = FALSE)
        result <- withr::with_envvar(c(DIA_STUB_SERVICE_EXIT = as.character(case[["service"]]),
            DIA_STUB_LOG_FAIL = as.character(case[["log"]])), processx::run("bash",
            c(runner, "--service", "browser-smoke", library_path, root), error_on_status = FALSE))
        directory <- setdiff(list.dirs(root, recursive = FALSE), before)
        expect_length(directory, 1L)
        expect_equal(result$status, case[["expected"]])
        expect_identical(readLines(file.path(directory, "service.log")), "retained-service-output")
        status <- readLines(file.path(directory, "controller-status.txt"))
        expect_true(paste0("service_exit=", case[["service"]]) %in% status)
        expect_true(paste0("log_exit=", case[["log"]]) %in% status)
        request <- readBin(file.path(directory, "controller-request.nul"), character(), n = 1000)
        expect_identical(request[[1L]], "dia-controller-request-v1")
        expect_identical(request[[3L]], "browser-smoke")
        expect_identical(request[[4L]], normalizePath(library_path))
        expect_identical(request[5:7], c("8589934592", "0", "900"))
        expect_identical(tail(request, 2), c(normalizePath(library_path), directory))
        decoded <- checkpointHarnessTestEnvironment()$readCheckpointControllerRequest(
            file.path(directory, "controller-request.nul"))
        expect_identical(decoded$command, request[-seq_len(7L)])
        expect_false(decoded$launch_observed)
    }
    # A readiness marker avoids probing flock itself and racing the lock owner.
    ready <- file.path(root, "service-ready")
    first <- withr::with_envvar(c(DIA_STUB_READY = ready), processx::process$new("bash",
        c(runner, "--service", "browser-smoke", library_path, root), stdout = "|", stderr = "|"))
    withr::defer(if (first$is_alive()) first$kill_tree())
    deadline <- Sys.time() + 5
    while (!file.exists(ready) && first$is_alive() && Sys.time() < deadline) Sys.sleep(0.01)
    expect_true(file.exists(ready))
    expect_true(first$is_alive())
    second <- processx::run("bash", c(runner, "--service", "browser-smoke", library_path, root),
        error_on_status = FALSE)
    expect_equal(second$status, 2)
    first$wait(timeout = 5000)
    expect_equal(first$get_exit_status(), 0)
    expect_equal(processx::run("bash", c(runner, "--service", "browser-smoke", library_path, root),
        error_on_status = FALSE)$status, 0)
})

test_that("controller request reader refuses invalid, oversized and unterminated records", {
    environment <- checkpointHarnessTestEnvironment()
    path <- withr::local_tempfile()
    fields <- c("dia-controller-request-v1", "multischolar-dia003-browser-smoke.ABC123.service",
        "browser-smoke", "/fixture library", "8589934592", "0", "900",
        "Rscript", "argument with\nquotes' and spaces")
    writeBin(fields, path)
    request <- environment$readCheckpointControllerRequest(path)
    expect_identical(request$command, fields[8:9])
    expect_false(request$completion_observed)
    for (index in c(1, 2, 3, 5, 6, 7)) {
        changed <- fields
        changed[[index]] <- "invalid"
        writeBin(changed, path)
        expect_error(environment$readCheckpointControllerRequest(path), "invalid controller request")
    }
    writeBin(charToRaw("unterminated"), path)
    expect_error(environment$readCheckpointControllerRequest(path), "unterminated")
    writeBin(fields[1:3], path)
    expect_error(environment$readCheckpointControllerRequest(path), "field count")
    writeBin(raw(1048577L), path)
    expect_error(environment$readCheckpointControllerRequest(path), "oversized")
    expect_error(environment$readCheckpointControllerRequest(tempfile()), "missing")
})

test_that("resource counters retain missing evidence and reject malformed values", {
    environment <- checkpointHarnessTestEnvironment()
    expect_identical(environment$parseCheckpointCounters(c("usage_usec 1200000", "user_usec 900000")),
        c(usage_usec = 1200000, user_usec = 900000))
    expect_length(environment$parseCheckpointCounters(character()), 0L)
    for (lines in list("x -1", "x NaN", "x 1.5", "x 1 extra", c("x 1", "x 2"))) {
        expect_error(environment$parseCheckpointCounters(lines))
    }
    missing <- environment$readCheckpointCounters(tempfile())
    expect_identical(missing$status, "unavailable")
    expect_null(missing$values)
    path <- withr::local_tempfile()
    writeLines(c("8:0 rbytes=10 wbytes=20", "8:1 rbytes=5 wbytes=0"), path)
    raw <- environment$readCheckpointCounters(path, raw = TRUE)
    expect_identical(raw$status, "available")
    expect_identical(raw$values, readLines(path))
    expect_identical(environment$readCheckpointCounters(path)$status, "unavailable")
    writeLines("oom_kill 1", path)
    expect_identical(environment$readCheckpointCounters(path)$values$oom_kill, 1)
})

test_that("DA export comparisons preserve order, populations and nonfinite values", {
    environment <- checkpointHarnessTestEnvironment()
    reference <- data.frame(Protein.Group = c("PG_B", "PG_A", "PG_C", "PG_D"),
        log2FC = c(1, NA_real_, Inf, NaN), comparison = "treated_vs_control")
    expect_true(environment$checkCheckpointDaTable(reference, reference))
    changed <- reference[4:1, ]
    expect_error(environment$checkCheckpointDaTable(reference, changed))
    changed <- reference
    changed$log2FC[[1L]] <- 2
    expect_error(environment$checkCheckpointDaTable(reference, changed))
    changed <- reference
    changed$log2FC[[3L]] <- -Inf
    expect_error(environment$checkCheckpointDaTable(reference, changed))
    changed <- reference
    changed$log2FC[[4L]] <- NA_real_
    expect_error(environment$checkCheckpointDaTable(reference, changed))
    expect_error(environment$checkCheckpointDaTable(reference, reference[-1L, ]))
    expect_error(environment$checkCheckpointDaTable(reference, reference["Protein.Group"]))
})

test_that("heatmap cluster checks reject changed identities, order and assignments", {
    environment <- checkpointHarnessTestEnvironment()
    clusters <- c(G2 = 1L, G1 = 2L, G3 = 1L)
    table <- data.frame(molecule_id = names(clusters), cluster_id = unname(clusters))
    expect_true(environment$checkCheckpointClusterTable(clusters, table))
    expect_error(environment$checkCheckpointClusterTable(clusters, table[3:1, ]))
    changed <- table
    changed$cluster_id[[1L]] <- 2L
    expect_error(environment$checkCheckpointClusterTable(clusters, changed))
    changed <- table
    changed$molecule_id[[1L]] <- "wrong"
    expect_error(environment$checkCheckpointClusterTable(clusters, changed))
})

checkpointHarnessProtocol <- function() {
    jsonlite::read_json(testthat::test_path("..", "testdata", "workflow-checkpoints",
        "protocol.json"), simplifyVector = TRUE)
}

checkpointHarnessReceiptFixture <- function(ratios = rep(0.8, 7)) {
    expected <- list(route = "iq", input_hash = "input", protocol_hash = "protocol",
        boundary_hash = "boundary", dependency_hash = "dependencies",
        control_build = "control-build", candidate_build = "candidate-build")
    control <- data.frame(pair_id = paste0("pair-", seq_len(7)), member = "control",
        run_id = paste0("control-", seq_len(7)), process_instance = paste0("process-c-", seq_len(7)),
        build_hash = expected$control_build, route = expected$route,
        input_hash = expected$input_hash, protocol_hash = expected$protocol_hash,
        boundary_hash = expected$boundary_hash, dependency_hash = expected$dependency_hash,
        status = "complete", science_ok = TRUE, outputs_ok = TRUE, diagnostics = FALSE,
        elapsed_seconds = 100, combined_peak_bytes = 1000000000, swap_peak_bytes = 0,
        cpu_seconds = 100, retained_heap_bytes = 100000000,
        max_save_seconds = 1, max_cold_revert_seconds = 1, stringsAsFactors = FALSE)
    candidate <- control
    candidate$member <- "candidate"
    candidate$run_id <- paste0("candidate-", seq_len(7))
    candidate$process_instance <- paste0("process-t-", seq_len(7))
    candidate$build_hash <- expected$candidate_build
    candidate$elapsed_seconds <- control$elapsed_seconds * ratios
    candidate$combined_peak_bytes <- control$combined_peak_bytes * ratios
    candidate$retained_heap_bytes <- control$retained_heap_bytes * 0.6
    list(receipts = rbind(control, candidate), expected = expected)
}

test_that("proof CLI requires explicit options and never loads the runtime for help", {
    environment <- checkpointHarnessTestEnvironment()
    expect_true(environment$parseCheckpointProofArgs("--help")$help)
    expect_output(environment$checkpointProofHelp(), "installed-package")
    expect_error(environment$parseCheckpointProofArgs(character()), "--mode")
    expect_error(environment$parseCheckpointProofArgs(c("--route", "iq", "--route", "iq")), "duplicate")
    expect_error(environment$parseCheckpointProofArgs("--mode"), "missing value")
    expect_error(environment$parseCheckpointProofArgs("--unknown"), "unknown option")
    args <- c("--mode", "baseline", "--route", "iq", "--contract", "contract.json",
        "--protocol", "protocol.json", "--library", "/not/a/library",
        "--output-root", "/not/a/root", "--diagnostics", "false", "--validate-only")
    parsed <- environment$parseCheckpointProofArgs(args)
    expect_false(parsed$diagnostics)
    expect_true(parsed[["validate-only"]])
    expect_error(environment$parseCheckpointProofArgs(c(args, "--checkpoint-codec", "qs2")),
        "utility mode only")
    args[2] <- "utility"
    expect_identical(environment$parseCheckpointProofArgs(args)[["checkpoint-codec"]], "rds")
    expect_identical(environment$parseCheckpointProofArgs(
        c(args, "--checkpoint-codec", "qs2"))[["checkpoint-codec"]], "qs2")
    expect_error(environment$parseCheckpointProofArgs(c(args, "--checkpoint-codec", "blob")),
        "rds or qs2")
    args[4] <- "lipidomics"
    expect_error(environment$parseCheckpointProofArgs(args), "unsupported route")
})

test_that("frozen protocol rejects missing fields and weakened limits", {
    environment <- checkpointHarnessTestEnvironment()
    protocol <- checkpointHarnessProtocol()
    expect_true(environment$validateCheckpointProofProtocol(protocol))
    for (field in c("independence_unit", "borderline_rule", "combined_peak_definition",
        "correction_budget", "run_budget")) {
        bad <- protocol
        bad[[field]] <- NULL
        expect_error(environment$validateCheckpointProofProtocol(bad), "missing|invalid")
    }
    bad <- protocol
    bad$approval$status <- "requested"
    expect_error(environment$validateCheckpointProofProtocol(bad), "approved and frozen")
    bad <- protocol
    bad$g4$elapsed_ratio_max <- 0.95
    expect_error(environment$validateCheckpointProofProtocol(bad), "threshold")
    bad <- protocol
    bad$service_limits$memory_max_bytes <- 16 * 1024^3
    expect_error(environment$validateCheckpointProofProtocol(bad), "resource limit")
    bad <- protocol
    bad$combined_peak_definition$metric <- "sum of process RSS plus cgroup peak"
    expect_error(environment$validateCheckpointProofProtocol(bad), "combined-peak")
    bad <- protocol
    bad$uncertainty_method$method <- "outcome-selected alternative"
    expect_error(environment$validateCheckpointProofProtocol(bad), "uncertainty procedure")
    bad <- protocol
    bad$run_budget$qualified_blocks_per_build_recipe_route <- NULL
    expect_error(environment$validateCheckpointProofProtocol(bad), "run budget")
    bad <- protocol
    bad$boundaries$cold_state <- NULL
    expect_error(environment$validateCheckpointProofProtocol(bad), "boundary")
})

test_that("paired evidence distinguishes pass fail borderline and exact boundaries", {
    environment <- checkpointHarnessTestEnvironment()
    protocol <- checkpointHarnessProtocol()
    fixture <- checkpointHarnessReceiptFixture()
    evaluate <- function(fixture) environment$evaluateCheckpointPairs(
        fixture$receipts, protocol, fixture$expected, gate = "g4")
    pass <- evaluate(fixture)
    expect_identical(pass$outcome, "accepted_evidence")
    expect_false(pass$promotion_authorized)
    expect_equal(pass$ratios$elapsed_seconds$point, 0.8)
    expect_identical(evaluate(checkpointHarnessReceiptFixture(rep(0.91, 7)))$outcome, "failed")
    expect_identical(evaluate(checkpointHarnessReceiptFixture(rep(0.90, 7)))$outcome,
        "accepted_evidence")
    borderline <- evaluate(checkpointHarnessReceiptFixture(c(0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1)))
    expect_identical(borderline$outcome, "held")
    expect_match(borderline$reason, "uncertainty")
    shuffled <- fixture
    shuffled$receipts <- shuffled$receipts[rev(seq_len(14)), ]
    expect_identical(evaluate(shuffled)$outcome, pass$outcome)
})

test_that("invalid control failures missing metrics and stale identities cannot pass", {
    environment <- checkpointHarnessTestEnvironment()
    protocol <- checkpointHarnessProtocol()
    fixture <- checkpointHarnessReceiptFixture()
    evaluate <- function(x) environment$evaluateCheckpointPairs(
        x, protocol, fixture$expected, gate = "g4")
    for (status in c("oom-kill", "timeout", "failed", "running")) {
        bad <- fixture$receipts
        bad$status[1] <- status
        expect_identical(evaluate(bad)$outcome, "held")
        bad <- fixture$receipts
        bad$status[8] <- status
        expect_identical(evaluate(bad)$outcome, if (status == "running") "held" else "failed")
    }
    for (field in c("route", "input_hash", "protocol_hash", "boundary_hash",
        "dependency_hash", "build_hash")) {
        bad <- fixture$receipts
        bad[[field]][8] <- "stale"
        expect_identical(evaluate(bad)$outcome, "held")
    }
    bad <- fixture$receipts
    bad$elapsed_seconds[8] <- NA_real_
    expect_identical(evaluate(bad)$outcome, "held")
    bad$elapsed_seconds[8] <- Inf
    expect_identical(evaluate(bad)$outcome, "held")
    bad <- fixture$receipts
    bad$outputs_ok[8] <- FALSE
    expect_identical(evaluate(bad)$outcome, "failed")
    bad <- fixture$receipts
    bad$diagnostics[8] <- TRUE
    expect_identical(evaluate(bad)$outcome, "held")
    bad <- fixture$receipts
    bad$process_instance[8:14] <- "one-session-seven-cycles"
    expect_identical(evaluate(bad)$outcome, "held")
    expect_identical(evaluate(fixture$receipts[-1, ])$outcome, "held")
    limpa <- fixture
    limpa$expected$route <- "limpa"
    limpa$receipts$route <- "limpa"
    expect_identical(environment$evaluateCheckpointPairs(limpa$receipts, protocol,
        limpa$expected, gate = "g4")$outcome, "held")
})

test_that("G2 allowance never leaks into G4 or ignores latency and service limits", {
    environment <- checkpointHarnessTestEnvironment()
    protocol <- checkpointHarnessProtocol()
    fixture <- checkpointHarnessReceiptFixture(rep(1.05, 7))
    evaluate <- function(x, gate = "g2_history", comparison = "original_control") {
        environment$evaluateCheckpointPairs(x, protocol, fixture$expected,
            gate = gate, comparison = comparison)
    }
    expect_identical(evaluate(fixture$receipts)$outcome, "accepted_evidence")
    expect_identical(evaluate(fixture$receipts, "g4")$outcome, "failed")
    expect_identical(evaluate(fixture$receipts, "g4", "accepted_comparator")$outcome, "failed")
    bad <- fixture$receipts
    bad$max_cold_revert_seconds[8] <- 31
    expect_identical(evaluate(bad)$outcome, "failed")
    bad <- fixture$receipts
    bad$swap_peak_bytes[8] <- 1
    expect_identical(evaluate(bad)$outcome, "failed")
    bad <- fixture$receipts
    bad$combined_peak_bytes[8] <- 9 * 1024^3
    expect_identical(evaluate(bad)$outcome, "failed")
    bad <- fixture$receipts
    bad$elapsed_seconds[1] <- 901
    expect_identical(evaluate(bad)$outcome, "held")
})

test_that("path checks refuse source-tree output and incomplete installations", {
    environment <- checkpointHarnessTestEnvironment()
    source_root <- testthat::test_path("..", "..")
    expect_error(environment$checkpointExistingDirectory(
        file.path(tempdir(), "missing-proof-directory"), "library"), "does not exist")
    expect_error(environment$validateCheckpointOutputRoot(source_root, source_root), "outside package")
    expect_error(environment$validateCheckpointOutputRoot(dirname(normalizePath(source_root)),
        source_root), "dedicated directory")
    empty_library <- withr::local_tempdir()
    expect_error(environment$inspectCheckpointLibrary(empty_library), "complete installed")
    if (identical(Sys.info()[["sysname"]], "Linux")) {
        filesystem <- system2("stat", c("-f", "-c", "%T", shQuote(tempdir())), stdout = TRUE)
        if (filesystem %in% c("tmpfs", "ramfs")) {
            expect_error(environment$validateCheckpointOutputRoot(tempdir(), source_root), "disk-backed")
        }
    }
})

test_that("diagnostic profiling wraps only the app launch and refuses profile overwrite", {
    skip_if_not(identical(Sys.info()[["sysname"]], "Linux"))
    skip_if_not(capabilities("profmem"))
    environment <- new.env(parent = baseenv())
    sys.source(testthat::test_path("..", "..", "tools", "profiling",
        "workflow_checkpoint_dia_driver.R"), envir = environment)
    root <- withr::local_tempdir()
    dir.create(file.path(root, "MultiScholaR"))
    path <- file.path(root, "application.Rprof")
    launch <- environment$checkpointProfiledApp(root, path)
    bindings <- base::environment(launch)
    expect_identical(sort(ls(bindings)), c("fixed_launch", "fixed_profile"))
    expect_identical(parent.env(bindings), baseenv())
    expect_identical(sort(ls(base::environment(bindings$fixed_launch))),
        c("fixed_capacity", "fixed_identity_root", "fixed_library", "fixed_package", "fixed_session_format"))
    expect_identical(base::environment(bindings$fixed_launch)$fixed_session_format, "legacy")
    writeLines("existing profile fixture", path)
    expect_error(launch(), "overwrite")
    expect_error(environment$checkpointProfiledApp(root, path), "existing")
    expect_error(environment$newCheckpointBrowser(root, root, profile_application = NA), "TRUE or FALSE")
})

test_that("application profiling remains active for work after the Shiny constructor returns", {
    skip_if_not(identical(Sys.info()[["sysname"]], "Linux"))
    skip_if_not(capabilities("profmem"))
    source_file <- normalizePath(testthat::test_path("..", "..", "tools", "profiling",
        "workflow_checkpoint_dia_driver.R"))
    path <- file.path(withr::local_tempdir(), "application.Rprof")
    result <- callr::r(function(source_file, path) {
        api <- new.env(parent = baseenv())
        sys.source(source_file, envir = api)
        workAfterConstruction <- function() {
            start <- proc.time()[["user.self"]]
            while (proc.time()[["user.self"]] - start < 0.3) sum(sqrt(seq_len(10000L)))
        }
        api$checkpointInstalledApp <- function(...) function() {
            shiny::shinyApp(ui = shiny::div(), server = function(input, output, session) {},
                onStart = function() later::later(function() {
                    workAfterConstruction()
                    shiny::stopApp(42L)
                }, 0.01))
        }
        runOnce <- function(profile_path, stop_early) {
            app <- api$checkpointProfiledApp("unused-fixture-library", profile_path)()
            if (stop_early) utils::Rprof(NULL)
            value <- shiny::runApp(app, host = "127.0.0.1", port = httpuv::randomPort(),
                launch.browser = FALSE, quiet = TRUE)
            profile <- utils::summaryRprof(profile_path)
            list(value = value, captured = any(grepl("workAfterConstruction", rownames(profile$by.total))))
        }
        list(construction_only = runOnce(paste0(path, "-construction-only"), TRUE),
            corrected = runOnce(path, FALSE))
    }, args = list(source_file, path))
    expect_identical(result$construction_only$value, 42L)
    expect_false(result$construction_only$captured)
    expect_identical(result$corrected$value, 42L)
    expect_true(result$corrected$captured)
})

test_that("installed launch closure cannot capture source scientific objects", {
    environment <- new.env(parent = baseenv())
    sys.source(testthat::test_path("..", "..", "tools", "profiling",
        "workflow_checkpoint_dia_driver.R"), envir = environment)
    environment$proc.time <- function() c(elapsed = 0)
    library_path <- withr::local_tempdir()
    dir.create(file.path(library_path, "MultiScholaR"))
    launch <- environment$checkpointInstalledApp(library_path)
    expect_identical(sort(ls(base::environment(launch))),
        c("fixed_capacity", "fixed_identity_root", "fixed_library", "fixed_package", "fixed_session_format"))
    expect_identical(base::environment(launch)$fixed_session_format, "legacy")
    native_launch <- environment$checkpointInstalledApp(library_path,
        identity_root = library_path, capacity_treatment = "qs2_duckdb", session_format = "native_v1")
    expect_identical(base::environment(native_launch)$fixed_session_format, "native_v1")
    expect_identical(sort(ls(base::environment(native_launch))), sort(ls(base::environment(launch))))
    expect_identical(parent.env(base::environment(launch)), baseenv())
    calls <- all.names(body(launch))
    expect_false("load_all" %in% calls)
    expect_true("loadNamespace" %in% calls)
    expect_match(paste(deparse(body(launch)), collapse = "\n"),
        "multischolar.capture_test_checkpoints = FALSE", fixed = TRUE)
    script <- environment$checkpointStepCompleteScript('Normalize"quoted')
    expect_match(script, "step-completed", fixed = TRUE)
    expect_false(grepl("test_state_digest", script, fixed = TRUE))
    expect_match(script, 'Normalize\\\"quoted', fixed = TRUE)
    actions <- environment$checkpointPeptideActions()
    expect_identical(vapply(actions, \(x) x$state, character(1)),
        c("qvalue_filtered", "protein_peptide_filtered", "precursor_rollup",
            "intensity_filtered", "sample_filtered", "replicate_filtered", "imputed"))
    expect_identical(actions[[2]]$tab, "Protein Evidence")
})

test_that("baseline QC leaves scientific defaults intact while tiny diagnostics remain explicit", {
    environment <- new.env(parent = baseenv())
    sys.source(testthat::test_path("..", "..", "tools", "profiling",
        "workflow_checkpoint_dia_driver.R"), envir = environment)
    environment$proc.time <- function() c(elapsed = 0)
    inputs <- list()
    environment$checkpointSetInput <- function(driver, id, value) inputs[[id]] <<- value
    completion_timeout <- NULL
    environment$waitCheckpointStepComplete <- function(driver, label, timeout) completion_timeout <<- timeout
    driver <- list(wait_for_js = function(...) NULL, get_js = function(...) NULL)
    peptide <- environment$runCheckpointPeptideActions(driver, diagnostic_settings = FALSE)
    protein <- environment$runCheckpointProteinActions(driver, diagnostic_settings = FALSE)
    expect_length(peptide, 7L)
    expect_length(protein, 5L)
    expect_identical(completion_timeout, 900000L)
    expect_false(any(grepl("min_peptides|min_peptidoforms|strict_mode|min_reps|min_groups|percentile|missing_values|parallel_cores",
        names(inputs))))
    expect_identical(inputs[["proteomics_workflow-quality_control-protein_qc-rollup-rollup_method"]], "iq")
    inputs <- list()
    expect_identical(environment$runCheckpointPeptideActions(driver), peptide)
    expect_identical(environment$runCheckpointProteinActions(driver), protein)
    expect_identical(inputs[["proteomics_workflow-quality_control-peptide_qc-imputation-proportion_missing_values"]], 1)
    expect_identical(inputs[["proteomics_workflow-quality_control-protein_qc-replicate_filter-parallel_cores"]], 1)
    expect_error(environment$runCheckpointPeptideActions(driver, NA), "TRUE or FALSE")
    expect_error(environment$runCheckpointProteinActions(driver, "false"), "TRUE or FALSE")
})

test_that("native chooser payload preserves exact filesystem paths", {
    environment <- new.env(parent = baseenv())
    sys.source(testthat::test_path("..", "..", "tools", "profiling",
        "workflow_checkpoint_dia_driver.R"), envir = environment)
    environment$proc.time <- function() c(elapsed = 0)
    root <- withr::local_tempdir()
    path <- file.path(root, "input's file.tsv")
    writeLines("column\nvalue", path)
    roots <- c(Test = root)
    selection <- environment$checkpointFileSelection(path, roots)
    expect_identical(selection$root, "Test")
    parsed <- shinyFiles::parseFilePaths(roots, selection)
    expect_identical(normalizePath(parsed$datapath), normalizePath(path))
    expect_error(environment$checkpointFileSelection(root, roots), "must be a file")
})

test_that("production option selection rejects missing ambiguous and duplicate values", {
    environment <- new.env(parent = baseenv())
    sys.source(testthat::test_path("..", "..", "tools", "profiling",
        "workflow_checkpoint_dia_driver.R"), envir = environment)
    environment$proc.time <- function() c(elapsed = 0)
    driver <- list(wait_for_js = function(...) NULL, wait_for_idle = function(...) NULL,
        get_js = function(script) {
            if (grepl("Object.keys", script, fixed = TRUE)) c("wt_1", "ko_1") else TRUE
        })
    expect_identical(environment$checkpointChooseExisting(driver, "input", "WT_1"), "wt_1")
    expect_error(environment$checkpointChooseExisting(driver, "input", "absent"), "missing")
    expect_error(environment$checkpointChooseExisting(driver, "input", c("wt_1", "wt_1"),
        multiple = TRUE), "duplicated")
    ambiguous <- driver
    ambiguous$get_js <- function(script) c("WT", "wt")
    expect_error(environment$checkpointChooseExisting(ambiguous, "input", "Wt"), "ambiguous")
})

test_that("input-only edits wait for idle rather than an unrelated rendered output", {
    environment <- new.env(parent = baseenv())
    sys.source(testthat::test_path("..", "..", "tools", "profiling",
        "workflow_checkpoint_dia_driver.R"), envir = environment)
    environment$proc.time <- function() c(elapsed = 0)
    captured <- new.env(parent = emptyenv())
    captured$idle <- 0L
    driver <- list(
        wait_for_js = function(...) NULL,
        get_js = function(...) FALSE,
        set_inputs = function(...) captured$arguments <- list(...),
        wait_for_idle = function(...) captured$idle <- captured$idle + 1L
    )
    environment$checkpointSetInput(driver, "manual_replicate_number", 1L)
    expect_false(captured$arguments$wait_)
    expect_identical(captured$arguments$manual_replicate_number, 1L)
    expect_identical(captured$idle, 1L)
})

test_that("normalization diagnostic uses explicit production controls and completion text", {
    environment <- new.env(parent = baseenv())
    sys.source(testthat::test_path("..", "..", "tools", "profiling",
        "workflow_checkpoint_dia_driver.R"), envir = environment)
    environment$proc.time <- function() c(elapsed = 0)
    actions <- list()
    inputs <- list()
    environment$checkpointSetInput <- function(driver, id, value) inputs[[id]] <<- value
    environment$checkpointClickAndWaitText <- function(driver, input_id, selector, text) {
        actions[[length(actions) + 1L]] <<- list(id = input_id, selector = selector, text = text)
    }
    result <- environment$runCheckpointNormalizationActions(NULL)
    expect_identical(result$method, "cyclicloess")
    expect_identical(result$ruv, "skip")
    expect_identical(inputs[["proteomics_workflow-normalization-norm_method"]], "cyclicloess")
    expect_identical(inputs[["proteomics_workflow-normalization-norm_qc_tabs"]],
        "Correlation Filtering")
    expect_identical(vapply(actions, `[[`, character(1), "id"),
        paste0("proteomics_workflow-normalization-", result$actions))
    expect_identical(actions[[3L]]$text, "Filtered session data exported successfully!")
})

test_that("production notification errors cannot satisfy action completion", {
    environment <- new.env(parent = baseenv())
    sys.source(testthat::test_path("..", "..", "tools", "profiling",
        "workflow_checkpoint_dia_driver.R"), envir = environment)
    environment$proc.time <- function() c(elapsed = 0)
    driver <- list(wait_for_js = function(...) NULL, get_js = function(script) {
        if (grepl("return el?el.textContent:null", script, fixed = TRUE)) return("failed output")
        TRUE
    })
    expect_error(environment$checkpointClickAndWaitText(driver, "action", ".result", "done"),
        "production action failed: action failed output")
    driver$get_js <- function(...) NULL
    expect_true(environment$checkpointClickAndWaitText(driver, "action", ".result", "done"))
})

test_that("DA diagnostic loads the real session before analysis and table display", {
    environment <- new.env(parent = baseenv())
    sys.source(testthat::test_path("..", "..", "tools", "profiling",
        "workflow_checkpoint_dia_driver.R"), envir = environment)
    environment$proc.time <- function() c(elapsed = 0)
    inputs <- list()
    actions <- character()
    environment$checkpointSetInput <- function(driver, id, value) inputs[[id]] <<- value
    environment$checkpointClickAndWaitText <- function(driver, input_id, selector, text) {
        actions <<- c(actions, input_id)
    }
    driver <- list(wait_for_js = function(...) NULL, get_js = function(...) NULL)
    result <- environment$runCheckpointDaActions(driver)
    expect_identical(actions, paste0("proteomics_workflow-differential_expression-",
        c("load_filtered_session", "run_da_analysis", "save_heatmap")))
    expect_identical(inputs[["proteomics_workflow-workflow_tabs"]], "da")
    expect_identical(inputs[["proteomics_workflow-differential_expression-formula_string"]], "~ 0 + group")
    expect_identical(inputs[["proteomics_workflow-differential_expression-table_significance"]], "all")
    expect_identical(result$treat_lfc, 0)
    expect_identical(inputs[["proteomics_workflow-differential_expression-volcano_interactive"]], TRUE)
    expect_true("static_volcano" %in% result$actions)
    expect_true("interactive_volcano" %in% result$actions)
    expect_identical(tail(result$actions, 2), c("heatmap", "save_heatmap"))
    driver$get_js <- function(...) "table error"
    expect_error(environment$runCheckpointDaActions(driver), "production DA table failed")
})

test_that("small export verifier rejects wrong counts, missing history and changed quantities", {
    environment <- checkpointHarnessTestEnvironment()
    raw <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    final <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    history <- c("initial", "raw_data_s4", "qvalue_filtered", "protein_peptide_filtered",
        "precursor_rollup", "intensity_filtered", "sample_filtered", "replicate_filtered",
        "imputed", "protein_s4_created", "protein_accession_cleaned",
        "protein_intensity_filtered", "duplicates_removed", "protein_replicate_filtered",
        "normalized", "correlation_filtered")
    # Validator unit inputs only: repeated fixture objects are not a scientific replay.
    states <- setNames(lapply(seq_along(history), \(i) list(data =
        if (i == 1L) NULL else if (i < 10L) raw else final)), history)
    session <- list(r6_state_history = as.list(history), r6_complete_states = states,
        r6_current_state_name = "correlation_filtered", current_s4_object = final,
        correlation_filtered_s4 = final, final_protein_count = 64, final_sample_count = 12)
    table <- final@protein_quant_table
    checked <- environment$checkCheckpointSessionContents(session, final, table)
    expect_identical(checked$history_entries, 16L)
    expect_equal(checked$proteins, 64)
    # Session provenance may contain environments. Independent hydration changes
    # their addresses, not their serialized contents.
    final@args$fixture_environment <- list2env(list(label = "immutable"), parent = baseenv())
    hydrated <- unserialize(serialize(final, NULL))
    expect_false(identical(final, hydrated))
    with_environment <- session
    with_environment$current_s4_object <- final
    with_environment$correlation_filtered_s4 <- final
    with_environment$r6_complete_states$correlation_filtered$data <- final
    expect_equal(environment$checkCheckpointSessionContents(with_environment, hydrated, table)$proteins, 64)
    final <- session$current_s4_object
    broken <- session
    broken$final_protein_count <- 0
    expect_error(environment$checkCheckpointSessionContents(broken, final, table),
        "exported protein count is inconsistent")
    broken <- session
    broken$r6_complete_states$raw_data_s4 <- NULL
    expect_error(environment$checkCheckpointSessionContents(broken, final, table))
    broken <- session
    broken$r6_complete_states$raw_data_s4$data@peptide_data$Q.Value <- NULL
    expect_error(environment$checkCheckpointSessionContents(broken, final, table), "Q value column")
    table[[2L]][1L] <- 999
    expect_error(environment$checkCheckpointSessionContents(session, final, table))
})
