#' Build one small native-session owner with an inactive retained branch.
#' @param root Fresh fixture project directory.
#' @return Manager and original scientific values; no GUI is started.
#' @noRd
makeDIANativeSessionFixture <- function(root) {
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "native-test", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    manager$setWorkflowType("DIA")
    raw <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    manager$saveState("raw_data_s4", raw, list(threshold = 0.01), "raw")
    final <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    manager$saveState("inactive", final, list(), "retained branch")
    manager$revertToState("raw_data_s4")
    manager$saveState("correlation_filtered", final, list(), "final")
    list(manager = manager, raw = raw, final = final,
        metadata = list(design_matrix = final@design_matrix,
            contrasts_tbl = data.frame(contrasts = "B-A", friendly_names = "test"),
            config_list = final@args, final_protein_count = nrow(final@protein_quant_table),
            final_sample_count = nrow(final@design_matrix), workflow_type = "DIA",
            normalization_method = "cyclicloess", ruv_mode = "skip", ruv_applied = FALSE))
}

test_that("component archives preserve complete closure and fail on a missing component", {
    fixture <- makeDIANativeSessionFixture(withr::local_tempdir())
    manager <- fixture$manager
    withr::defer(closeDIAStateAdapter(manager))
    convertDIANativeLegacyCheckpoints(manager)
    source <- getDIAStateReference(manager)
    acknowledgeDIAOutputs(manager, source, "correlation")
    before <- manager$.context$full_s4_reads
    result <- testthat::with_mocked_bindings(
        runDIANativeSessionExport(manager, fixture$metadata, withr::local_tempdir()),
        materializeDIAStateRecord = function(...) stop("native export hydrated a state"))
    expect_identical(manager$.context$full_s4_reads, before)
    imported <- openDIANativeSession(result$publication$archive_path, withr::local_tempdir())
    withr::defer(closeDIAStateAdapter(imported$manager))
    expect_true(usesDIAReferenceStates(imported$manager))
    expect_identical(imported$manifest$schema, "multischolar.dia_native_session.v3")
    expect_identical(imported$manager$getState(), fixture$final)
    view <- resolveDIANativeSessionView(imported$manager,
        imported$manager$.context$native_metadata_reference, imported$extra_references)
    expect_true(isSavedWorkflowReference(view$current_s4_object))
    expect_true(validateProtDaFilteredSession(view, manager = imported$manager))
    component <- readDIAStateMetadata(imported$manager)$data$components$protein_quant_table
    path <- file.path(imported$manager$.context$store$directory, component$payload$relative_path)
    unlink(path)
    expect_error(readDIAExecutionComponents(imported$manager), class = "saved_entity_corrupt")
})

test_that("logical table bindings export their physical dependency and restore attributes", {
    fixture <- makeDIANativeSessionFixture(withr::local_tempdir())
    manager <- fixture$manager
    withr::defer(closeDIAStateAdapter(manager))
    convertDIANativeLegacyCheckpoints(manager)
    expected <- fixture$final
    sample <- as.character(expected@design_matrix[[expected@sample_id]][[1L]])
    at <- which(is.finite(expected@protein_quant_table[[sample]]))[[1L]]
    expected@protein_quant_table[[sample]][[at]] <- expected@protein_quant_table[[sample]][[at]] + 1
    physical <- saveDIATableComponent(manager$.context, expected@protein_quant_table)
    scope <- newWorkflowTableReadScope(manager$.context$store, physical, reader = "arrow", metadata_only = TRUE)
    metadata <- scope$metadata
    metadata$attributes$logical_note <- "metadata-only restoration"
    metadata$attribute_order <- c(metadata$attribute_order, "logical_note")
    binding <- saveWorkflowTableComponent(scope, metadata, "logical-view",
        newSavedWorkflowScientificIdentity(physical$entity_id, "view", "note", "metadata", "test"))
    releaseWorkflowTableReadScope(scope)
    attr(expected@protein_quant_table, "logical_note") <- "metadata-only restoration"
    checkpoint <- readDIAStateMetadata(manager)$data
    checkpoint$components$protein_quant_table <- binding
    state <- commitDIAReferenceState(manager, "correlation_filtered", checkpoint, list(), "logical view")
    acknowledgeDIAOutputs(manager, state, "correlation")
    result <- testthat::with_mocked_bindings(
        runDIANativeSessionExport(manager, fixture$metadata, withr::local_tempdir()),
        materializeDIAStateRecord = function(...) stop("native export hydrated an assay"))
    imported <- openDIANativeSession(result$publication$archive_path, withr::local_tempdir())
    withr::defer(closeDIAStateAdapter(imported$manager))
    expect_identical(imported$manager$getState(), expected)
    expect_true(exists(physical$entity_id, imported$manager$.context$store$records, inherits = FALSE))
    unlink(file.path(imported$manager$.context$store$directory, physical$payload$relative_path))
    expect_error(imported$manager$getState(), class = "saved_entity_corrupt")
})

test_that("native export streams keep named file caches bounded and propagate I/O failures", {
    skip_if_not(identical(Sys.info()[["sysname"]], "Linux"))
    receipt <- Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT")
    skip_if(!nzchar(receipt), "disk-backed guarded proof root required for cache residency")
    disk_root <- dirname(receipt)
    filesystem <- system2("stat", c("-f", "-c", "%T", "--", shQuote(disk_root)), stdout = TRUE)
    expect_true(filesystem %in% c("btrfs", "ext2/ext3", "xfs", "zfs"))
    root <- withr::local_tempdir(tmpdir = disk_root)
    output <- withr::local_tempdir(tmpdir = disk_root)
    fixture <- makeDIANativeSessionFixture(root)
    withr::defer(closeDIAStateAdapter(fixture$manager))
    # Non-aligned length and shell metacharacters exercise exact bytes/path quoting.
    source <- file.path(root, "input ' $ sample.bin")
    connection <- file(source, "wb")
    block <- rep(as.raw(0:255), 16384L)
    for (index in seq_len(8L)) writeBin(block, connection)
    writeBin(as.raw(0:122), connection)
    close(connection)
    expected <- digest::digest(file = source, algo = "sha256")
    cache_bytes <- function(path) {
        command <- unname(Sys.which("fincore"))
        expect_true(nzchar(command))
        result <- system2(command, c("--json", "--bytes", "--output", "RES", "--", shQuote(path)),
            stdout = TRUE, stderr = TRUE)
        expect_null(attr(result, "status"))
        jsonlite::fromJSON(paste(result, collapse = "\n"))$fincore$res[[1L]]
    }
    result <- runDIANativeSessionExport(fixture$manager, fixture$metadata, output,
        dependencies = list(input = list(source = source, destination = "inputs/input.bin")))
    archive <- result$publication$archive_path
    staged <- file.path(output, ".workflow-export", result$publication$request_id, "files", expected)
    # Inspect residency before correctness reads deliberately warm these files.
    residency <- vapply(c(source, staged, archive), cache_bytes, numeric(1))
    expect_true(all(residency <= 4194304))
    expect_identical(digest::digest(file = staged, algo = "sha256"), expected)
    expect_identical(file.size(staged), file.size(source))
    imported <- openDIANativeSession(archive, withr::local_tempdir())
    withr::defer(closeDIAStateAdapter(imported$manager))
    expect_identical(digest::digest(file = imported$dependencies$input$path, algo = "sha256"), expected)
    expect_identical(imported$manager$getState(), fixture$final)
    log_path <- file.path(output, "io-failure.log")
    expect_error(runDIANativeIO(list(c("dd", paste0("if=", file.path(root, "missing")),
        "status=none"), "sha256sum"), log_path), class = "saved_entity_write")
    expect_error(runDIANativeTar(c("--create", "--file", "-", "--directory", root,
        basename(source)), log_path, output_path = output), class = "saved_entity_write")
    if (nzchar(receipt)) jsonlite::write_json(list(source_bytes = file.size(source),
        archive_bytes = file.size(archive), resident_bytes = residency,
        full_reads = diaStateOwnership(fixture$manager)$full_reads,
        content_sha256 = expected), file.path(dirname(receipt), "native-io-residency.json"),
        auto_unbox = TRUE, pretty = TRUE)
})

test_that("native archives restore complete references without hydrating history", {
    skip_if_not(identical(Sys.info()[["sysname"]], "Linux"))
    root <- withr::local_tempdir()
    output <- withr::local_tempdir()
    restored_root <- withr::local_tempdir()
    fixture <- makeDIANativeSessionFixture(root)
    manager <- fixture$manager
    withr::defer(closeDIAStateAdapter(manager))
    expected_names <- names(manager$.context$snapshot$states)
    expected_history <- manager$getHistory()
    source_file <- file.path(root, "input.tsv")
    writeLines(c("id\tvalue", "P1\t3"), source_file)
    before <- diaStateOwnership(manager)
    result <- testthat::with_mocked_bindings(
        runDIANativeSessionExport(manager, fixture$metadata, output,
            dependencies = list(input = list(source = source_file, destination = "inputs/input.tsv"))),
        materializeDIAStates = function(...) stop("history hydration is forbidden"),
        hydrateSavedWorkflowCheckpoint = function(...) stop("scientific hydration is forbidden"))
    expect_identical(result$publication$status, "complete")
    expect_identical(result$publication$scientific_checkpoint_count, 0L)
    expect_true(file.exists(result$publication$archive_path))
    expect_false(file.exists(file.path(output, "filtered_session_data_latest.rds")))
    expect_null(manager$.context$export_pending)
    expect_identical(manager$getHistory(), expected_history)
    expect_identical(diaStateOwnership(manager)$full_reads, before$full_reads)
    expect_identical(diaStateOwnership(manager)$hydrates, before$hydrates)
    pointer <- jsonlite::read_json(result$publication$latest_path)
    expect_identical(pointer$archive, basename(result$publication$archive_path))
    expect_lt(file.size(result$publication$latest_path), 8192)
    closeDIAStateAdapter(manager)
    # The source store is deliberately unavailable to the importer.
    moved <- paste0(root, "-unavailable")
    expect_true(file.rename(root, moved))
    withr::defer(file.rename(moved, root))
    imported <- openDIANativeSession(result$publication$archive_path, restored_root)
    withr::defer(closeDIAStateAdapter(imported$manager))
    expect_identical(imported$metadata, fixture$metadata)
    expect_identical(names(imported$manager$.context$snapshot$states), expected_names)
    expect_identical(imported$manager$getHistory(), expected_history)
    expect_identical(imported$manager$getState(), fixture$final)
    expect_identical(imported$manager$getState("inactive"), fixture$final)
    expect_identical(imported$manager$revertToState("raw_data_s4"), fixture$raw)
    expect_identical(readLines(imported$dependencies$input$path), c("id\tvalue", "P1\t3"))
    expect_identical(diaStateOwnership(imported$manager)$full_reads, 0L)
    expect_identical(diaStateOwnership(imported$manager)$working_records, 0L)
    expect_identical(imported$manager$.context$store$open_resources, 0L)
    expect_null(imported$manager$.context$owned_export)
    expect_length(imported$manager$.context$owned_workers, 0L)
})

test_that("native output retries preserve the previous latest and do not repeat science", {
    root <- withr::local_tempdir()
    output <- withr::local_tempdir()
    fixture <- makeDIANativeSessionFixture(root)
    manager <- fixture$manager
    withr::defer(closeDIAStateAdapter(manager))
    accepted <- runDIANativeSessionExport(manager, fixture$metadata, output)
    old_pointer <- readBin(accepted$publication$latest_path, "raw",
        n = file.size(accepted$publication$latest_path))
    manager$saveState("next_state", fixture$final, list(), "new head")
    history <- manager$getHistory()
    for (boundary in c("before_native_archive", "after_native_archive_publication",
        "before_native_latest_publication")) {
        failed <- FALSE
        expect_error(testthat::with_mocked_bindings(
            runDIANativeSessionExport(manager, fixture$metadata, output),
            workflowExportBoundary = function(phase) {
                if (identical(phase, boundary)) { failed <<- TRUE; stop("injected native failure") }
            }), "injected native failure")
        expect_true(failed)
        expect_identical(manager$getHistory(), history)
        expect_identical(readBin(accepted$publication$latest_path, "raw",
            n = file.size(accepted$publication$latest_path)), old_pointer)
        expect_error(assertDIAExportReady(manager), class = "saved_entity_outputs_pending")
        # Simulate loss of the in-memory request; recovery reads the durable journal.
        manager$.context$export_pending <- NULL
        recovered <- testthat::with_mocked_bindings(
            runDIANativeSessionExport(manager, NULL, output),
            materializeDIAStates = function(...) stop("unexpected history read"),
            saveDIAStateValue = function(...) stop("unexpected scientific/metadata save"))
        expect_identical(recovered$publication$status, "complete")
        expect_identical(recovered$publication$scientific_checkpoint_count, 0L)
        expect_identical(manager$getHistory(), history)
        expect_null(manager$.context$export_pending)
        accepted <- recovered
        old_pointer <- readBin(accepted$publication$latest_path, "raw",
            n = file.size(accepted$publication$latest_path))
    }
})

test_that("native imports reject corrupt files and unsafe archive members", {
    root <- withr::local_tempdir()
    output <- withr::local_tempdir()
    fixture <- makeDIANativeSessionFixture(root)
    withr::defer(closeDIAStateAdapter(fixture$manager))
    result <- runDIANativeSessionExport(fixture$manager, fixture$metadata, output)
    staging <- withr::local_tempdir()
    runDIANativeTar(c("--extract", "--file", result$publication$archive_path,
        "--directory", staging), file.path(staging, "unpack.log"))
    manifest <- jsonlite::read_json(file.path(staging, "manifest.json"))
    payload <- grep("payload[.]qs2$", names(manifest$files), value = TRUE)[[1L]]
    path <- file.path(staging, payload)
    bytes <- readBin(path, "raw", n = file.size(path))
    bytes[[length(bytes)]] <- as.raw(bitwXor(as.integer(bytes[[length(bytes)]]), 1L))
    writeBin(bytes, path)
    members <- c("manifest.json", names(manifest$files))
    member_list <- file.path(staging, "members.nul")
    writeBin(members, member_list)
    corrupted <- file.path(output, "corrupted.msdia.tar")
    runDIANativeTar(c("--create", "--format=pax", "--no-recursion", "--file", corrupted,
        "--directory", staging, "--null", "--files-from", member_list), file.path(staging, "pack.log"))
    expect_error(openDIANativeSession(corrupted, withr::local_tempdir()), class = "saved_entity_corrupt")
    for (bad in list("../escape", "/absolute", "files/../escape", c("manifest.json", "manifest.json"))) {
        expect_error(checkDIANativeMemberNames(bad), class = "saved_entity_path")
    }
    link <- file.path(staging, "link")
    expect_true(file.symlink(path, link))
    unsafe <- file.path(output, "unsafe.msdia.tar")
    runDIANativeTar(c("--create", "--format=pax", "--file", unsafe, "--directory", staging,
        "link"), file.path(staging, "link-pack.log"))
    expect_error(openDIANativeSession(unsafe, withr::local_tempdir()), class = "saved_entity_path")
    expect_error(runDIANativeSessionExport(fixture$manager,
        list(r6_complete_states = list()), output), class = "saved_entity_invalid")
    expect_error(collectDIANativeReferences(list(value = fixture$raw)), class = "saved_entity_invalid")
})

test_that("the independent native verifier checks history without collecting it", {
    runner <- new.env(parent = baseenv())
    sys.source(test_path("..", "..", "tools", "profiling", "run_workflow_checkpoint_proof.R"), runner)
    fixture <- makeDIANativeSessionFixture(withr::local_tempdir())
    manager <- fixture$manager
    withr::defer(closeDIAStateAdapter(manager))
    states <- c("raw_data_s4", "qvalue_filtered", "protein_peptide_filtered", "precursor_rollup",
        "intensity_filtered", "sample_filtered", "replicate_filtered", "imputed", "protein_s4_created",
        "protein_accession_cleaned", "protein_intensity_filtered", "duplicates_removed",
        "protein_replicate_filtered", "ruv_corrected", "correlation_filtered")
    manager$revertToState("raw_data_s4")
    for (name in states[-1L]) manager$saveState(name, fixture$final, list(), name)
    session <- c(fixture$metadata, list(r6_state_history = manager$state_history,
        r6_current_state_name = manager$current_state, current_s4_object = fixture$final,
        correlation_filtered_s4 = fixture$final, ruv_k = 1L))
    session$ruv_applied <- TRUE
    session$ruv_mode <- "automatic"
    table <- fixture$final@protein_quant_table
    table[[fixture$final@protein_id_column]] <- as.character(table[[fixture$final@protein_id_column]])
    result <- testthat::with_mocked_bindings(
        runner$checkCheckpointSessionContents(session, fixture$final, table,
            ruv_diagnostic = TRUE, native_manager = manager),
        materializeDIAStates = function(...) stop("unexpected full history collection"))
    expect_equal(result$history_entries, 16)
    expect_identical(manager$current_state, "raw_data_s4")
    expect_identical(diaStateOwnership(manager)$full_reads, 0L)
    expect_identical(diaStateOwnership(manager)$working_records, 0L)
    session$r6_state_history <- session$r6_state_history[-2L]
    expect_error(runner$checkCheckpointSessionContents(session, fixture$final, table,
        ruv_diagnostic = TRUE, native_manager = manager))
})

#' Create only the existing output directories used by the native GUI fixtures.
#' @param root Temporary fixture project.
#' @return Experiment-path list.
#' @noRd
makeDIANativeGuiPaths <- function(root) {
    paths <- list(source_dir = file.path(root, "scripts", "proteomics"),
        protein_qc_dir = file.path(root, "results", "proteomics", "protein_qc"),
        peptide_qc_dir = file.path(root, "results", "proteomics", "peptide_qc"))
    for (path in c(unlist(paths), file.path(dirname(paths$protein_qc_dir), "publication_graphs"))) {
        dir.create(path, recursive = TRUE, showWarnings = FALSE)
    }
    paths
}

#' Provide session-lifetime callbacks without starting a browser for pure fixtures.
#' @return Minimal session-like owner.
#' @noRd
makeDIANativeFixtureSession <- function() {
    session <- new.env(parent = emptyenv())
    session$userData <- new.env(parent = emptyenv())
    session$onFlush <- session$onFlushed <- session$onSessionEnded <- function(...) function() invisible(NULL)
    session
}

test_that("native GUI export and owned DA load do not materialize historical states", {
    paths <- makeDIANativeGuiPaths(withr::local_tempdir())
    fixture <- makeDIANativeSessionFixture(withr::local_tempdir())
    manager <- fixture$manager
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$session_format <- "native_v1"
    manager$.context$experiment_paths <- paths
    workflow <- list2env(fixture$metadata, parent = emptyenv())
    workflow$state_manager <- manager
    workflow$ruv_optimization_result <- list(ruv_skipped = TRUE)
    workflow$data_format <- "diann"
    workflow$data_type <- "peptide"
    norm <- new.env(parent = emptyenv())
    norm$correlation_filtered_obj <- saveDIAOwnedValue(fixture$final, "correlation_filtered_obj", manager)
    norm$correlation_threshold <- 0.5
    testthat::with_mocked_bindings(
        saveProtNormCorrelationResults(fixture$final, paths, workflow$ruv_optimization_result,
            saveRdsFn = function(...) stop("automatic S4 export is forbidden")),
        getDIAActiveManager = function() manager)
    expect_false(file.exists(file.path(paths$protein_qc_dir,
        "normalised_results_cln_with_replicates.RDS")))
    inputs <- list(norm_method = "cyclicloess", ruv_mode = "skip")
    reader <- hydrateSavedWorkflowCheckpoint
    result <- testthat::with_mocked_bindings(
        runProtNormExportSessionWorkflow(workflow, norm, inputs, paths$source_dir,
            withProgressFn = function(message, value, expr) force(expr), incProgressFn = function(...) NULL),
        materializeDIAStates = function(...) stop("unexpected history hydration"),
        hydrateSavedWorkflowCheckpoint = function(store, reference) {
            if (!identical(reference$identity$parameters, "native-session-metadata")) stop("unexpected assay read")
            reader(store, reference)
        })
    expect_null(result$sessionData)
    expect_identical(result$exportArtifacts$publication$status, "complete")
    expect_match(result$exportArtifacts$sessionFilepath, "[.]msdia[.]tar$")
    expect_true(file.exists(result$exportArtifacts$summaryFilepath))
    required_files <- collectDIANativeGuiFiles(manager, collectDIANativeSessionMetadata(workflow, norm, inputs)$metadata)
    expect_true("source/ruv_optimization_results.RDS" %in%
        vapply(required_files, `[[`, character(1), "destination"))
    expect_true(all(file.exists(vapply(required_files, `[[`, character(1), "source"))))
    expect_identical(manager$.context$full_reads, 0L)
    expect_identical(diaStateOwnership(manager)$working_records, 0L)
    view <- testthat::with_mocked_bindings(
        withDIAStateScope(manager, function() resolveDIANativeOwnedExport(manager,
            result$exportArtifacts$latestFilepath)),
        materializeDIAStates = function(...) stop("unexpected history hydration"),
        describeDIANativeFile = function(...) stop("owned DA must not rehash the archive"))
    expect_identical(view$current_s4_object, fixture$final)
    expect_identical(view$correlation_filtered_s4, fixture$final)
    expect_null(view$r6_complete_states)
    fresh_paths <- makeDIANativeGuiPaths(withr::local_tempdir())
    imported <- openDIANativeSession(result$exportArtifacts$sessionFilepath, fresh_paths$source_dir)
    withr::defer(closeDIAStateAdapter(imported$manager))
    restored <- restoreDIANativeGuiFiles(imported, fresh_paths)
    expect_type(restored, "character")
    expect_true(length(restored) > 0L && all(file.exists(restored)))
    expect_false(file.exists(file.path(fresh_paths$protein_qc_dir,
        "normalised_results_cln_with_replicates.RDS")))
    expect_identical(imported$manager$getState(), fixture$final)
    portable <- file.path(fresh_paths$source_dir, "integration.msdia.tar")
    copyDIANativeSessionArtifact(paths$source_dir, portable)
    expect_identical(digest::digest(file = portable, algo = "sha256"),
        digest::digest(file = result$exportArtifacts$sessionFilepath, algo = "sha256"))
    testthat::with_mocked_bindings(
        expect_identical(copyDIANativeSessionArtifact(paths$source_dir, portable), portable),
        runDIANativeIO = function(...) stop("unchanged archive was copied again"))
    expect_identical(digest::digest(file = portable, algo = "sha256"),
        digest::digest(file = result$exportArtifacts$sessionFilepath, algo = "sha256"))
    expect_error(copyDIANativeSessionArtifact(paths$source_dir, sub("msdia.tar$", "RDS", portable)),
        class = "saved_entity_path")
    integration <- file.path(fresh_paths$source_dir, "integration")
    values <- new.env(parent = emptyenv())
    output <- new.env(parent = emptyenv())
    expect_true(completeProtSummaryWorkflowArgsSave(output, values,
        projectDirs = list(proteomics = c(paths, list(integration_dir = integration))),
        workflowData = workflow, experimentLabel = "native_fixture", description = "native output",
        assignFn = function(...) NULL,
        createWorkflowArgsFn = function(...) file.path(paths$source_dir, "study_parameters.txt"),
        saveRDSFn = function(...) stop("automatic Summary S4 export is forbidden"),
        renderTextFn = function(expr) force(expr), showNotificationFn = function(...) NULL,
        catFn = function(...) NULL))
    expect_identical(digest::digest(file = file.path(integration,
        "proteomics_native_fixture_final.msdia.tar"), algo = "sha256"),
        digest::digest(file = portable, algo = "sha256"))
    expect_false(file.exists(file.path(integration, "proteomics_native_fixture_final_s4.RDS")))
    expect_identical(restoreDIANativeGuiFiles(imported, fresh_paths), character())
    writeLines("conflict", file.path(fresh_paths$protein_qc_dir, "normalised_results_cln_with_replicates.tsv"))
    expect_error(restoreDIANativeGuiFiles(imported, fresh_paths), class = "saved_entity_conflict")
})

test_that("native legacy import validates and persists before replacing the GUI owner", {
    withr::local_options(multischolar.dia_native_sessions = TRUE)
    paths <- makeDIANativeGuiPaths(withr::local_tempdir())
    fixture <- makeDIANativeSessionFixture(withr::local_tempdir())
    withr::defer(closeDIAStateAdapter(fixture$manager))
    legacy <- c(fixture$metadata, list(r6_complete_states = fixture$manager$states,
        r6_state_history = fixture$manager$getHistory(), r6_current_state_name = fixture$manager$current_state,
        current_s4_object = fixture$final, correlation_filtered_s4 = fixture$final))
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- WorkflowState$new(audit_enabled = FALSE)
    previous <- workflow$state_manager
    bad <- legacy
    bad$r6_complete_states$raw_data_s4$data <- fixture$manager$.context$snapshot$states$raw_data_s4
    session <- makeDIANativeFixtureSession()
    expect_error(ingestDIANativeLegacySession(bad, workflow, session, paths$source_dir, paths),
        class = "saved_entity_corrupt")
    expect_identical(workflow$state_manager, previous)
    expect_true(ingestDIANativeLegacySession(legacy, workflow, session, paths$source_dir, paths))
    withr::defer(closeDIAStateAdapter(workflow$state_manager))
    expect_true(isDIANativeSession(workflow$state_manager))
    # KIT-018: the adopted history is a reference history (every state a saved reference, no resident S4).
    expect_true(usesDIAReferenceStates(workflow$state_manager))
    for (state in workflow$state_manager$getHistory())
        expect_true(isSavedWorkflowReference(getDIAStateReference(workflow$state_manager, state)))
    expect_identical(workflow$state_manager$getState(), fixture$final)
    expect_identical(workflow$state_manager$revertToState("raw_data_s4"), fixture$raw)
    expect_identical(workflow$state_manager$.context$full_reads, 0L)
    expect_false(file.exists(file.path(paths$protein_qc_dir, "normalised_results_cln_with_replicates.RDS")))
    expect_true(file.exists(file.path(paths$protein_qc_dir, "normalised_results_cln_with_replicates.tsv")))
    saveProtNormCorrelationResults(fixture$final, paths,
        normRuvOptimizationResult = list(ruv_skipped = TRUE), compatibilityExport = TRUE)
    expect_identical(readRDS(file.path(paths$protein_qc_dir,
        "normalised_results_cln_with_replicates.RDS")), fixture$final)
    expect_null(workflow$state_manager$.context$export_pending)
})

test_that("the design checkpoint on the DIA default commits a reference state through newDIAPeptideCheckpoint (KIT-018)", {
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "design-default", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$reference_states <- TRUE
    design <- object@design_matrix
    design$tech_rep_group <- NULL
    workflow <- list2env(list(state_manager = manager, design_matrix = design, config_list = object@args,
        data_cln = saveDIATableComponent(manager$.context, object@peptide_data)), parent = emptyenv())
    expect_true(usesDIAReferenceStates(manager))
    expect_identical(buildProtDesignStateCheckpoint(workflow, "DIA", "Import"), "raw_data_s4")
    expect_identical(manager$current_state, "raw_data_s4")
    expect_true(isSavedWorkflowReference(getDIAStateReference(manager, "raw_data_s4")))
    expect_identical(workflow$design_matrix$tech_rep_group, paste(design$group, design$replicates, sep = "_"))
    restored <- manager$getState("raw_data_s4")
    expect_s4_class(restored, "PeptideQuantitativeData")
    expect_identical(nrow(restored@peptide_data), nrow(object@peptide_data))
    expect_identical(restored@design_matrix$tech_rep_group, workflow$design_matrix$tech_rep_group)
})

test_that("native import is the DIA-NN/IQ default and does not migrate Limpa or other formats", {
    paths <- makeDIANativeGuiPaths(withr::local_tempdir())
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- WorkflowState$new(audit_enabled = FALSE)
    session <- makeDIANativeFixtureSession()
    # KIT-011: no option gate; the former `multischolar.dia_native_sessions = FALSE` arm no longer exists.
    withr::local_options(multischolar.dia_native_sessions = NULL)
    expect_true(isDIANativeDefaultEnabled())
    original <- workflow$state_manager
    initializeDIANativeImport(workflow, paths, session, "fragpipe", "iq")
    expect_identical(workflow$state_manager, original)
    initializeDIANativeImport(workflow, paths, session, "diann", "limpa")
    expect_identical(workflow$state_manager, original)
    initializeDIANativeImport(workflow, paths, session, "diann", "iq")
    withr::defer(closeDIAStateAdapter(workflow$state_manager))
    expect_true(isDIANativeSession(workflow$state_manager))
    expect_identical(workflow$state_manager$.context$experiment_paths, paths)
    expect_error(initializeDIANativeImport(workflow, paths, session, "diann", "limpa"),
        class = "saved_entity_scope")
    expect_identical(workflow$state_manager$.context$store$max_cache_bytes, 0)
})

test_that("ten native normalization revert style export cycles release GUI-owned payloads", {
    withr::local_options(multischolar.dia_native_sessions = TRUE,
        multischolar.capture_test_checkpoints = FALSE)
    paths <- makeDIANativeGuiPaths(withr::local_tempdir())
    paths$publication_graphs_dir <- file.path(dirname(paths$protein_qc_dir), "publication_graphs")
    paths$time_dir <- file.path(paths$publication_graphs_dir, "filtering_qc")
    dir.create(paths$time_dir, recursive = TRUE)

    original <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    receipts <- list()
    receipt_path <- Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT")
    inspect_display_owners <- function(plot) {
        if (!inherits(plot, c("gg", "ggplot"))) return(list())
        frames <- c(list(plot$plot_env), lapply(plot$mapping, function(value)
            if (rlang::is_quosure(value)) rlang::quo_get_env(value) else NULL))
        lapply(frames, function(frame) {
            if (!is.environment(frame) || identical(frame, baseenv()) ||
                identical(frame, globalenv()) || isNamespace(frame)) return(list())
            fields <- intersect(c("theObject", "s4Object", "normalized", "object", "Y",
                "data", "data_filtered", "pca.res", "da_results_list"), ls(frame, all.names = TRUE))
            stats::setNames(lapply(fields, function(field) {
                value <- get(field, frame, inherits = FALSE)
                list(class = class(value), dimensions = dim(value), bytes = as.numeric(object.size(value)))
            }), fields)
        })
    }
    shiny::testServer(function(input, output, session) {
        workflow <- shiny::reactiveValues(state_manager = WorkflowState$new(audit_enabled = FALSE),
            design_matrix = original@design_matrix, config_list = original@args,
            contrasts_tbl = data.frame(contrasts = "grouptreated-groupcontrol"))
        initializeDIANativeImport(workflow, paths, session, "diann", "iq")
        manager <- workflow$state_manager
        manager$.context$reclaim_enabled <- TRUE
        manager$saveState("protein_replicate_filtered", original, list(), "source")
        norm <- shiny::reactiveValues()
        rlang::local_bindings(config_list = original@args, .env = .GlobalEnv)
        shiny::observeEvent(input$cycle, {
          withDIAStateScope(manager, function() {
            testthat::with_mocked_bindings({
              selectDIAStateReference(manager, "protein_replicate_filtered")
              resetProtNormReactiveState(norm, workflow)
              aesthetics <- function() list(color_var = "group", shape_var = NULL)
              generateProtNormPreNormalizationQcArtifacts(manager, paths$protein_qc_dir, aesthetics())
              ruv_mode <- if (input$cycle == 10L) "manual" else "skip"
              settings <- list(norm_method = input$method, ruv_mode = ruv_mode, ruv_k = 1,
                ruv_percentage = 50, color_variable = "group", shape_variable = NULL)
              post <- function(reference) {
                norm$qc_plot_paths <- generateProtNormPostNormalizationQcArtifacts(reference,
                  norm$publication_paths$protein_qc_dir %||% paths$protein_qc_dir,
                  aesthetics(), norm$qc_plot_paths)
              }
              ruv <- function(reference) {
                norm$qc_plot_paths <- generateProtNormRuvCorrectedQcArtifacts(reference,
                  norm$publication_paths$protein_qc_dir %||% paths$protein_qc_dir,
                  aesthetics(), norm$qc_plot_paths)
              }
              expect_true(runDIAReferenceNormalization(settings, workflow, norm, paths,
                "proteomics", "cycles", post, ruv, function() "group"))
            }, materializeDIAStateRecord = function(...) stop("GUI lifecycle hydrated an S4 state"))
          })
        })
        shiny::observeEvent(input$finish, {
          withDIAStateScope(manager, function() {
            testthat::with_mocked_bindings({
              aesthetics <- function() list(color_var = "group", shape_var = NULL)
              settings <- list(norm_method = input$method,
                ruv_mode = if (input$cycle == 10L) "manual" else "skip", ruv_k = 1,
                ruv_percentage = 50, color_variable = "group", shape_variable = NULL)
              skipped <- input$cycle %% 2L != 0L
              final <- runDIAReferenceCorrelation(workflow, norm, paths,
                0.5, skipped, "proteomics", "cycles", function() "group", aesthetics)
              expect_true(isSavedWorkflowReference(final))
              # The module edge completes the GUI state after the declaration (KIT-004).
              completeProtNormCorrelationWorkflow(final, workflow, output, norm,
                correlationThreshold = if (skipped) 0 else 0.5, skipped = skipped,
                successNotification = "done", completionMessage = "done", messagePrefix = "Correlation",
                showNotificationFn = function(...) NULL, messageFn = function(...) NULL)
              first <- getDIAStateReference(manager)
              if (input$cycle %% 2L == 0L) {
                testthat::with_mocked_bindings(runDIAReferenceCorrelation(workflow, norm, paths,
                  0.5, FALSE, "proteomics", "cycles", function() "group", aesthetics),
                  pearsonCorForSamplePairs = function(...) stop("Pearson refitted on reapply"))
                expect_identical(getDIAStateReference(manager), first)
              }
              plot <- resolveDIAReferenceQcPlot(final, aesthetics(), "correlation", buildProtNormQcPlotProduct,
                plotPca, plotRle, buildProtNormDensityPlot, buildProtNormCorrelationPlot)
              norm$render_warnings <- character()
              withCallingHandlers(saveProtNormQcPlotArtifact(paths$protein_qc_dir,
                "post_norm_correlation.png", plot, 4, 3, dpi = 100), warning = function(condition) {
                  norm$render_warnings <- c(norm$render_warnings, conditionMessage(condition))
                  invokeRestart("muffleWarning")
                })
              norm$publication <- runProtNormExportSessionWorkflow(workflow, norm, settings, paths$source_dir,
                withProgressFn = function(message, value, expr) force(expr), incProgressFn = function(...) NULL)
            }, materializeDIAStateRecord = function(...) stop("GUI lifecycle hydrated an S4 state"))
          })
        })
    }, {
        for (cycle in seq_len(10L)) {
            method <- c("cyclicloess", "quantile", "scale", "none")[(cycle - 1L) %% 4L + 1L]
            session$setInputs(method = method, cycle = cycle)
            session$flushReact()
            release <- manager$.context$reclaim_result
            expect_true(release$status %in% c("released", "no_pages_released", "unsupported"))
            expect_identical(release$label, "normalization")
            expect_identical(release$state, getDIAStateReference(manager)$entity_id)
            expect_null(manager$.context$reclaim_pending)
            normalization_owners <- inspect_display_owners(ggplot2::get_last_plot())
            session$flushReact()
            expect_identical(manager$.context$reclaim_result, release)
            session$setInputs(finish = cycle)
            session$flushReact()
            expect_true(isSavedWorkflowReference(norm$normalized_protein_obj))
            expect_true(isSavedWorkflowReference(norm$ruv_normalized_obj))
            expect_true(isSavedWorkflowReference(norm$correlation_filtered_obj))
            expect_identical(norm$publication$exportArtifacts$publication$status, "complete")
            imported <- openDIANativeSession(norm$publication$exportArtifacts$sessionFilepath,
                withr::local_tempdir())
            final <- manager$getState()
            expected <- normaliseBetweenSamples(original, normalisation_method = method)
            if (cycle == 10L) expected <- ruvIII_C_Varying(expected, ruv_grouping_variable = "group",
                ruv_number_k = 1L, ctrl = norm$control_genes_index)
            if (cycle == 10L) expected <- removeRowsWithMissingValuesPercent(expected)
            if (cycle %% 2L == 0L) expected <- filterSamplesByProteinCorrelationThreshold(expected,
                pearson_correlation_per_pair = pearsonCorForSamplePairs(expected,
                    tech_rep_remove_regex = "pool", correlation_group = "group"),
                min_pearson_correlation_threshold = 0.5)
            expect_equal(final@protein_quant_table, expected@protein_quant_table, tolerance = 1e-12)
            expect_identical(final@design_matrix, expected@design_matrix)
            expect_identical(imported$manager$getState(), final)
            expect_false(isS4(workflow$ruv_normalised_for_da_analysis_obj))
            expect_true(isSavedWorkflowReference(workflow$ruv_normalised_for_da_analysis_obj))
            expect_identical(imported$manager$getHistory(), manager$getHistory())
            expect_identical(imported$metadata$normalization_method, method)
            expect_identical(imported$metadata$ruv_mode, if (cycle == 10L) "manual" else "skip")
            png_copy <- imported$dependencies[vapply(imported$dependencies,
                \(entry) identical(entry$destination, "protein_qc/post_norm_correlation.png"), logical(1))]
            expect_length(png_copy, 1L)
            expect_identical(digest::digest(file = png_copy[[1L]]$path, algo = "sha256"),
                digest::digest(file = file.path(paths$protein_qc_dir, "post_norm_correlation.png"), algo = "sha256"))
            closeDIAStateAdapter(imported$manager)
            session$flushReact()
            info <- diaStateOwnership(manager)
            expect_identical(info$full_reads, 0L)
            expect_identical(info$working_records, 0L)
            expect_identical(info$native_records, 0L)
            expect_identical(info$store$open_resources, 0L)
            expect_equal(info$store$bytes, 0)
            expect_lte(info$store$metadata_bytes, 64 * 1024^2)
            expect_length(manager$.context$owned_workers, 0L)
            expect_null(manager$.context$export_pending)
            rss <- grep("^VmRSS:", readLines("/proc/self/status"), value = TRUE)
            receipts[[cycle]] <- list(cycle = cycle, normalization = method,
                reclamation = release,
                normalization_global_owners = normalization_owners,
                final_global_owners = inspect_display_owners(ggplot2::get_last_plot()),
                final_qc_owners = inspect_display_owners(norm$final_qc_plot),
                render_warnings = norm$render_warnings,
                rss_bytes = as.numeric(sub("^VmRSS:[[:space:]]*([0-9]+).*", "\\1", rss)) * 1024,
                ownership = info)
            if (nzchar(receipt_path)) jsonlite::write_json(list(cycles_completed = cycle,
                observations = receipts, performance_acceptance = FALSE),
                file.path(dirname(receipt_path), "native-cycle-checks.json"), auto_unbox = TRUE,
                pretty = TRUE, null = "null", digits = NA)
        }
        expect_length(receipts, 10L)
    })
})

test_that("a native archive supports real IQ reapply after its source project is unavailable", {
    root <- withr::local_tempdir()
    paths <- makeDIANativeGuiPaths(root)
    output <- withr::local_tempdir()
    fixture <- makeDIANativeSessionFixture(paths$source_dir)
    manager <- fixture$manager
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$session_format <- "native_v1"
    peptide <- fixture$raw
    keys <- c(peptide@protein_id_column, peptide@sample_id, peptide@peptide_sequence_column)
    peptide@peptide_data <- peptide@peptide_data[!duplicated(peptide@peptide_data[keys]), ]
    peptide@peptide_data$Peptide.Imputed <- 1000 + seq_len(nrow(peptide@peptide_data))
    manager$saveState("imputed", peptide, list(), "synthetic post-imputation IQ input")
    workflow <- list2env(list(state_manager = manager), parent = emptyenv())
    apply_iq <- function(workflow, paths) withDIAStateScope(workflow$state_manager,
        function() runProteinIqRollupApplyStep(workflow, paths,
            captureCheckpointFn = function(...) NULL, showNotificationFn = function(...) NULL))
    original <- apply_iq(workflow, paths)$proteinObj
    metadata <- fixture$metadata
    metadata$design_matrix <- original@design_matrix
    metadata$final_protein_count <- nrow(original@protein_quant_table)
    metadata$final_sample_count <- nrow(original@design_matrix)
    result <- runDIANativeSessionExport(manager, metadata, output,
        dependencies = resolveDIANativeDependencies(manager))
    closeDIAStateAdapter(manager)
    moved <- paste0(root, "-unavailable")
    expect_true(file.rename(root, moved))
    withr::defer(file.rename(moved, root))
    fresh_paths <- makeDIANativeGuiPaths(withr::local_tempdir())
    imported <- openDIANativeSession(result$publication$archive_path, fresh_paths$source_dir)
    workflow$state_manager <- WorkflowState$new(audit_enabled = FALSE)
    view <- adoptDIANativeSession(imported, workflow, makeDIANativeFixtureSession(),
        fresh_paths$source_dir, fresh_paths)
    withr::defer(closeDIAStateAdapter(workflow$state_manager))
    expect_true(isSavedWorkflowReference(view$current_s4_object))
    expect_identical(workflow$state_manager$getState(), original)
    expect_identical(workflow$state_manager$revertToState("imputed"), peptide)
    expect_true(all(file.exists(vapply(imported$dependencies, `[[`, character(1), "path"))))
    result <- apply_iq(workflow, fresh_paths)
    expect_true(isSavedWorkflowReference(result$stateReference))
    actual <- workflow$state_manager$getState()
    expect_identical(actual@protein_quant_table, original@protein_quant_table)
    expect_identical(actual@protein_id_table, original@protein_id_table)
    expect_identical(actual@design_matrix, original@design_matrix)
    expect_identical(actual@args, original@args)
    expect_identical(diaStateOwnership(workflow$state_manager)$full_reads, 0L)
    expect_identical(diaStateOwnership(workflow$state_manager)$working_records, 0L)
    expect_identical(workflow$state_manager$.context$store$open_resources, 0L)
})

test_that("native publication recovers after process loss without another scientific checkpoint", {
    skipIfNoInstalledRuntime()
    library <- dirname(getNamespaceInfo(asNamespace("MultiScholaR"), "path"))
    writer <- function(library, root, output, metadata, phase, ready) {
        ns <- loadNamespace("MultiScholaR", lib.loc = library)
        store <- get("newWorkflowEntityStore", ns)(root, "native-test", "DIA")
        manager <- get("newDIAStateAdapter", ns)(store, reopen = TRUE)
        testthat::with_mocked_bindings(
            get("runDIANativeSessionExport", ns)(manager, metadata, output),
            workflowExportBoundary = function(current) {
                if (identical(current, phase)) {
                    writeLines("ready", ready)
                    repeat Sys.sleep(0.05)
                }
            }, .package = "MultiScholaR")
    }
    environment(writer) <- baseenv()
    for (phase in c("before_native_archive", "after_native_archive_publication",
        "before_native_completion_acknowledgement")) {
        root <- withr::local_tempdir()
        output <- withr::local_tempdir()
        fixture <- makeDIANativeSessionFixture(root)
        manager <- fixture$manager
        accepted <- runDIANativeSessionExport(manager, fixture$metadata, output)
        old_latest <- digest::digest(file = accepted$publication$latest_path, algo = "sha256")
        manager$saveState("next", fixture$final, list(), "committed science before failed export")
        head <- manager$.context$head
        history <- manager$getHistory()
        closeDIAStateAdapter(manager)
        ready <- file.path(root, "writer-ready")
        process <- callr::r_bg(writer, list(library, root, output, fixture$metadata, phase, ready),
            stdout = file.path(root, "writer-out.log"), stderr = file.path(root, "writer-error.log"))
        withr::defer(if (process$is_alive()) process$kill())
        deadline <- Sys.time() + 30
        while (!file.exists(ready) && process$is_alive() && Sys.time() < deadline) Sys.sleep(0.05)
        expect_true(file.exists(ready), info = phase)
        if (!file.exists(ready)) stop(paste("native writer did not reach", phase,
            paste(readLines(file.path(root, "writer-error.log"), warn = FALSE), collapse = "\n")))
        process$signal(9L)
        process$wait(timeout = 5000)
        expect_false(process$is_alive())
        expect_equal(process$get_exit_status(), -9L)
        if (phase != "before_native_completion_acknowledgement") {
            expect_identical(digest::digest(file = accepted$publication$latest_path, algo = "sha256"), old_latest)
        }
        recovered <- newDIAStateAdapter(newWorkflowEntityStore(root, "native-test", "DIA"), reopen = TRUE)
        withr::defer(closeDIAStateAdapter(recovered))
        request <- findDIANativePendingExport(recovered, output)
        expect_identical(request$status, "outputs_pending")
        expect_identical(recovered$.context$head, head)
        result <- testthat::with_mocked_bindings(runDIANativeSessionExport(recovered, NULL, output),
            saveDIAStateValue = function(...) stop("recovery must not save science or metadata"),
            materializeDIAStates = function(...) stop("recovery must not hydrate history"))
        expect_identical(result$publication$status, "complete")
        expect_identical(result$publication$scientific_checkpoint_count, 0L)
        expect_identical(recovered$.context$head, head)
        expect_identical(recovered$getHistory(), history)
        expect_identical(recovered$.context$writes, 0L)
        expect_null(recovered$.context$export_pending)
        expect_identical(diaStateOwnership(recovered)$working_records, 0L)
        expect_identical(recovered$.context$store$open_resources, 0L)
        closeDIAStateAdapter(recovered)
    }
})

# KIT-013: the hand-off's checkpoint defects (items 2, 3 and 5).

makeKit013CheckpointWorkflow <- function(native, object) {
    shiny::reactiveValues(state_manager = native, data_tbl = object@peptide_data,
        data_cln = object@peptide_data, data_type = "peptide", data_format = "diann",
        config_list = list(globalParameters = list(workflow_type = "DIA")),
        design_matrix = object@design_matrix,
        column_mapping = list(run_col = "Run", protein_col = "Protein.Group",
            peptide_col = "Stripped.Sequence", raw_quantity_col = "Precursor.Quantity",
            norm_quantity_col = "Precursor.Normalised"))
}

makeKit013CheckpointPaths <- function(root, name) {
    paths <- list(source_dir = file.path(root, paste0(name, "-source")),
        peptide_qc_dir = file.path(root, paste0(name, "-peptide")))
    lapply(paths, dir.create)
    paths
}

test_that("KIT-013 item 3: a peptide checkpoint taken at every accepted state reopens", {
    # The hand-off claimed raw_data_s4 and qvalue_filtered could not be reopened because the
    # validator requires a design. The design module commits raw_data_s4 at design time, so a
    # checkpoint at any accepted peptide state carries one; this reproduces the claim and shows
    # every state reopening. Direction recorded in kit-013/outcome.md: no validator change.
    recipe <- readWorkflowCheckpointRecipe("small")
    recipe$n_proteins <- 8L; recipe$n_samples <- 8L
    recipe$peptides_per_protein <- 4L; recipe$forms_per_peptide <- 2L
    object <- makeWorkflowCheckpointPeptide(recipe)
    root <- withr::local_tempdir()
    states <- c("raw_data_s4", "qvalue_filtered", "protein_peptide_filtered", "precursor_rollup",
        "intensity_filtered", "sample_filtered", "replicate_filtered", "imputed")
    writeLines("accepted output", file.path(makeKit013CheckpointPaths(root, "save")$peptide_qc_dir, "peptide.tsv"))
    save_paths <- list(source_dir = file.path(root, "save-source"), peptide_qc_dir = file.path(root, "save-peptide"))
    shiny::testServer(function(input, output, session) {}, {
        session$userData$peptide_checkpoint_design <- function() list(
            design_matrix = object@design_matrix, groups = unique(object@design_matrix$group),
            factors = character(), contrasts = data.frame(), removed_samples = character(),
            formula = "~ 0 + group")
        for (state in states) {
            native <- WorkflowState$new(audit_enabled = FALSE)
            native$setWorkflowType("DIA")
            for (earlier in states[seq_len(match(state, states))]) native$saveState(earlier, object, list(), earlier)
            saving <- makeKit013CheckpointWorkflow(native, object)
            publication <- saveDIAPeptideCheckpoint(saving, session, save_paths)
            opening <- makeKit013CheckpointWorkflow(WorkflowState$new(audit_enabled = FALSE), object)
            opening$state_manager$setWorkflowType("DIA")
            open_paths <- makeKit013CheckpointPaths(root, paste0("open-", state))
            expect_identical(openDIAPeptideCheckpoint(publication$archive_path, opening, session, open_paths), state)
            manager <- opening$state_manager
            expect_true(usesDIAReferenceStates(manager))
            expect_identical(manager$getHistory()[[length(manager$getHistory())]], state)
            expect_identical(manager$getState(state)@peptide_data, object@peptide_data)
            closeDIAStateAdapter(manager)
        }
    })
})

test_that("KIT-013 items 2 and 5: a failed Open leaves no restored file and a successful Open removes the replaced folder", {
    recipe <- readWorkflowCheckpointRecipe("small")
    recipe$n_proteins <- 8L; recipe$n_samples <- 8L
    recipe$peptides_per_protein <- 4L; recipe$forms_per_peptide <- 2L
    object <- makeWorkflowCheckpointPeptide(recipe)
    root <- withr::local_tempdir()
    save_paths <- makeKit013CheckpointPaths(root, "save")
    writeLines("accepted output", file.path(save_paths$peptide_qc_dir, "peptide.tsv"))
    native <- WorkflowState$new(audit_enabled = FALSE)
    native$setWorkflowType("DIA")
    native$saveState("intensity_filtered", object, list(), "before")
    native$saveState("sample_filtered", object, list(), "after")
    shiny::testServer(function(input, output, session) {}, {
        session$userData$peptide_checkpoint_design <- function() list(
            design_matrix = object@design_matrix, groups = unique(object@design_matrix$group),
            factors = character(), contrasts = data.frame(), removed_samples = character(),
            formula = "~ 0 + group")
        publication <- saveDIAPeptideCheckpoint(makeKit013CheckpointWorkflow(native, object), session, save_paths)
        workflow <- makeKit013CheckpointWorkflow(native, object)
        paths <- makeKit013CheckpointPaths(root, "project")
        restored <- file.path(paths$peptide_qc_dir, "peptide.tsv")
        private_folders <- function() list.files(paths$source_dir, "^[.]dia-", all.files = TRUE, include.dirs = TRUE)
        # Item 2: the restore runs after the rollback is registered; a failure after the owner
        # swap puts the bus back and removes the file this Open restored.
        with_mocked_bindings({
            expect_error(openDIAPeptideCheckpoint(publication$archive_path, workflow, session, paths),
                "injected restore failure")
        }, workflowExecutionBoundary = function(stage) {
            if (stage == "peptide_checkpoint_after_owner_swap") stop("injected restore failure")
        }, .package = "MultiScholaR")
        expect_identical(workflow$state_manager, native)
        expect_false(file.exists(restored))
        expect_identical(private_folders(), character())
        # A stale project (the bus changed under the Open) is refused before any file is restored.
        stale <- makeKit013CheckpointWorkflow(native, object)
        with_mocked_bindings({
            expect_error(openDIAPeptideCheckpoint(publication$archive_path, stale, session, paths),
                class = "workflow_execution_stale")
        }, validateSavedWorkflowReference = function(store, reference) {
            stale$state_manager <- WorkflowState$new(audit_enabled = FALSE)
            invisible(TRUE)
        }, .package = "MultiScholaR")
        expect_false(file.exists(restored))
        expect_identical(private_folders(), character())
        # Item 5: the first Open installs its private folder; the second removes the first.
        expect_identical(openDIAPeptideCheckpoint(publication$archive_path, workflow, session, paths), "sample_filtered")
        expect_identical(readLines(restored), "accepted output")
        first_root <- dirname(workflow$state_manager$.context$store$directory)
        expect_true(dir.exists(first_root))
        expect_true(startsWith(basename(first_root), ".dia-native-import-"))
        expect_identical(private_folders(), basename(first_root))
        expect_identical(openDIAPeptideCheckpoint(publication$archive_path, workflow, session, paths), "sample_filtered")
        second_root <- dirname(workflow$state_manager$.context$store$directory)
        expect_false(identical(first_root, second_root))
        expect_false(dir.exists(first_root))
        expect_true(dir.exists(second_root))
        expect_identical(private_folders(), basename(second_root))
        # The rule never touches a user's folder, a folder elsewhere, or the adopted owner.
        user_folder <- file.path(paths$source_dir, "my-results")
        dir.create(user_folder)
        expect_false(removeDIAReplacedSessionFolder(user_folder, second_root, paths))
        expect_true(dir.exists(user_folder))
        expect_false(removeDIAReplacedSessionFolder(second_root, second_root, paths))
        expect_true(dir.exists(second_root))
        elsewhere <- file.path(root, ".dia-native-import-elsewhere")
        dir.create(elsewhere)
        expect_false(removeDIAReplacedSessionFolder(elsewhere, second_root, paths))
        expect_true(dir.exists(elsewhere))
        closeDIAStateAdapter(workflow$state_manager)
    })
})

test_that("KIT-013 audit PX-01: a QC-time configuration update does not fire the design sync observer", {
    recipe <- readWorkflowCheckpointRecipe("small")
    recipe$n_proteins <- 8L; recipe$n_samples <- 8L
    recipe$peptides_per_protein <- 4L; recipe$forms_per_peptide <- 2L
    object <- makeWorkflowCheckpointPeptide(recipe)
    root <- withr::local_tempdir()
    native <- WorkflowState$new(audit_enabled = FALSE)
    native$saveState("intensity_filtered", object, list(), "input")
    workflow <- shiny::reactiveValues(state_manager = native, data_tbl = object@peptide_data,
        data_cln = object@peptide_data, config_list = list(globalParameters = list(workflow_type = "DIA"),
            srlQvalueProteotypicPeptideClean = list(qvalue_threshold = 0.05)),
        data_type = "peptide", column_mapping = list(run_col = "Run", protein_col = "Protein.Group",
            peptide_col = "Stripped.Sequence", raw_quantity_col = "Precursor.Quantity",
            norm_quantity_col = "Precursor.Normalised"))
    shiny::testServer(function(input, output, session) {}, {
        session$userData$workflow_data <- workflow
        manager <- stageDIAGUIReferenceOwner(native, list(source_dir = root))
        adoptDIAGUIReferenceOwner(manager, workflow, session)
        builder <- createProtDesignBuilderState()
        initial_runs <- 0L
        sync_runs <- 0L
        initial_state <- shiny::reactive({
            initial_runs <<- initial_runs + 1L
            buildDIAReferenceDesignState(manager, workflow$data_tbl, workflow$config_list, workflow$column_mapping)
        })
        count_design <- function(value) { sync_runs <<- sync_runs + 1L; builder$designMatrix(value) }
        observer <- registerProtDesignInitialStateSyncObserver(input, session,
            initialState = initial_state, dataTbl = function() workflow$data_tbl,
            designMatrix = count_design, dataClnReactive = builder$dataClnReactive,
            groups = builder$groups, factors = builder$factors, contrasts = builder$contrasts,
            removedSamples = builder$removedSamples, resultRv = builder$resultRv)
        session$flushReact()
        # On adoption the observer consumes the pending Design result and returns before writing selections.
        expect_identical(sync_runs, 0L)
        expect_identical(initial_runs, 1L)
        # A QC Apply writes the configuration through the bus environment (correction 1): neither the
        # initial-state reactive nor the sync observer runs again.
        state <- manager$getState("intensity_filtered")
        updated <- suppressMessages(updateConfigParameter(state, "srlQvalueProteotypicPeptideClean", "qvalue_threshold", 0.01))
        expect_equal(shiny::isolate(workflow$config_list)$srlQvalueProteotypicPeptideClean$qvalue_threshold, 0.01)
        session$flushReact()
        expect_identical(sync_runs, 0L)
        expect_identical(initial_runs, 1L)
        # The design initial state is rebuilt only on its own events (a new data reference or Reset).
        session$setInputs(reset_changes = 1)
        session$flushReact()
        expect_identical(sync_runs, 1L)
        expect_identical(initial_runs, 2L)
        observer$destroy()
        closeDIAStateAdapter(manager)
    })
})
