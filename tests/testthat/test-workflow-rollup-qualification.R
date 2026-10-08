test_that("completed writer measurements support independent readback recovery", {
    mode <- Sys.getenv("MULTISCHOLAR_DIA025_MECHANISM_MEMBER")
    skip_if(!grepl("^gui-bulk-rollup:writer-r-readback:[A-Za-z0-9]{6}$", mode),
        "explicit verification of a retained completed operation only")
    root <- dirname(Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT"))
    proof <- dirname(root)
    member <- sub("^gui-bulk-rollup:writer-r-readback:", "", mode)
    original <- file.path(proof, paste0("browser-smoke.", member))
    receipt_path <- file.path(original, "writer-comparison.json")
    before <- digest::digest(file = receipt_path, algo = "sha256")
    receipt <- jsonlite::read_json(receipt_path, simplifyVector = TRUE)
    stopifnot(identical(receipt$case, "contrast"), !is.null(receipt$result$complete_operation_seconds),
        "teardown" %in% receipt$phases$phase)
    actual <- fingerprintDIAWriterOutput(file.path(original, paste0("arm-", receipt$arm)),
        receipt$result$reference)
    oracle_root <- file.path(proof, "browser-smoke.EYBmKx")
    oracle <- jsonlite::read_json(file.path(oracle_root, "rollup-qualification.json"), simplifyVector = TRUE)
    expected <- fingerprintDIAWriterOutput(file.path(oracle_root, "arm-B"), oracle$result$reference,
        receipt$samples)
    expect_identical(actual, expected)
    expect_identical(actual, receipt$result$output_hashes)
    expect_length(receipt$result$figures, 12L)
    expect_identical(digest::digest(file = receipt_path, algo = "sha256"), before)
    jsonlite::write_json(list(status = "complete", exact = identical(actual, expected),
        source_member = member, source_receipt_sha256 = before, computation_repeated = FALSE,
        output_hashes = actual), file.path(root, "writer-readback.json"), auto_unbox = TRUE, pretty = TRUE)
})

test_that("writer candidates use identical complete sample populations", {
    mode <- Sys.getenv("MULTISCHOLAR_DIA025_MECHANISM_MEMBER")
    skip_if(!grepl("^gui-bulk-rollup:writer-r:(screen|full):(primary|contrast):P1[456]:[12]$", mode),
        "explicit unprofiled writer comparison only")
    fields <- strsplit(mode, ":", fixed = TRUE)[[1L]]
    screen <- fields[[3L]] == "screen"
    case <- fields[[4L]]
    arm <- fields[[5L]]
    root <- dirname(Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT"))
    proof <- dirname(root)
    withr::local_envvar(c(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
        DUCKPLYR_FALLBACK_AUTOUPLOAD = "0", DUCKPLYR_FALLBACK_COLLECT = "0"))
    if (case == "primary") source <- jsonlite::read_json(file.path(proof,
        "browser-smoke.BpSyxB", "gui-state-latest.json")) else {
        manifests <- lapply(Sys.glob(file.path(proof, "browser-smoke.*",
            "rollup-qualification-inputs.json")), jsonlite::read_json)
        valid <- which(vapply(manifests, \(x) identical(x$status, "complete") &&
            identical(x$revision, "p3-adaptations-v1"), logical(1)))
        stopifnot(length(valid) == 1L)
        source <- manifests[[valid]]$cases$contrast
    }
    head <- file.path(source$store_directory, "dia-state-head.rds")
    head_hash <- digest::digest(file = head, algo = "sha256")
    stopifnot(identical(head_hash, if (case == "primary")
        "5d64ca562669d94eea1677a1a03d5b410a8a81d8ea59add8cabd6065deb91353" else source$head_sha256))
    owned <- file.path(root, paste0("arm-", arm))
    dir.create(owned)
    manager <- copyDIAReplayCheckpoint(source$store_directory, owned)
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    selectDIAStateReference(manager, "protein_peptide_filtered")
    checkpoint <- readDIAStateMetadata(manager)$data
    stopifnot(checkpoint$dimensions$rows == if (case == "primary") 24000000L else 3840000L)
    design <- withDIAStateScope(manager, function()
        readDIAExecutionComponents(manager, fields = "design_matrix"))@design_matrix
    samples <- if (!screen) NULL else as.character(design[[checkpoint$slots$sample_id]][
        if (case == "primary") c(1:4, 149:152, 297:300) else c(1:2, 24:25, 47:48)])
    receipt <- list(status = "running", mode = mode, arm = arm, case = case, samples = samples,
        source_head_sha256 = head_hash, profiling = FALSE, phases = list(), failures = list())
    receipt_path <- file.path(root, "writer-comparison.json")
    flush <- function() jsonlite::write_json(receipt, receipt_path, auto_unbox = TRUE,
        pretty = TRUE, null = "null", digits = NA)
    on.exit(flush(), add = TRUE)
    phase <- function(label, fn) {
        receipt$active_phase <<- label
        flush()
        started <- proc.time()
        value <- fn()
        receipt$phases[[length(receipt$phases) + 1L]] <<- list(phase = label,
            seconds = unname((proc.time() - started)[["elapsed"]]),
            cpu = unname(sum((proc.time() - started)[c("user.self", "sys.self")])))
        flush()
        value
    }
    tryCatch({
        started <- proc.time()[["elapsed"]]
        value <- runDIABulkOperation(manager, owned, proof, arm, samples, phase)
        phase("teardown", function() closeDIAStateAdapter(manager))
        receipt$result <- c(value, list(complete_operation_seconds = proc.time()[["elapsed"]] - started))
        hashes <- fingerprintDIAWriterOutput(owned, value$reference)
        receipt$result$output_hashes <- hashes
        oracle_id <- if (case == "primary") if (screen) "QGWkjJ" else "60v50g" else "EYBmKx"
        oracle_root <- file.path(proof, paste0("browser-smoke.", oracle_id))
        oracle <- jsonlite::read_json(file.path(oracle_root, if (case == "primary")
            "bulk-rollup.json" else "rollup-qualification.json"), simplifyVector = TRUE)
        expected <- if (case == "primary") oracle$results[[1L]] else oracle$result
        if (screen && case == "contrast") {
            # Reuse the independently checked complete-sample native fingerprints.
            # Re-reading the entire wide native reference in every service would
            # measure oracle construction in its memory peak, not this operation.
            verified <- jsonlite::read_json(file.path(proof, "browser-smoke.1wXXDD",
                "writer-readback.json"), simplifyVector = TRUE)
            original_path <- file.path(proof, "browser-smoke.CHnC1c", "writer-comparison.json")
            original <- jsonlite::read_json(original_path, simplifyVector = TRUE)
            stopifnot(identical(verified$status, "complete"), isTRUE(verified$exact),
                identical(verified$source_member, "CHnC1c"),
                identical(verified$source_receipt_sha256, digest::digest(file = original_path, algo = "sha256")),
                identical(original$samples, samples))
            expected$output_hashes <- verified$output_hashes
        }
        expect_identical(hashes, expected$output_hashes)
        if (!screen || case == "primary") {
            expect_identical(jsonlite::toJSON(value$summary, auto_unbox = TRUE),
                jsonlite::toJSON(expected$summary, auto_unbox = TRUE))
        }
        expect_length(value$figures, 12L)
        expect_identical(value$execution$engine_threads, 8L)
        expect_identical(value$execution$writer_threads, 8L)
        expect_equal(value$execution$intermediate_bytes_written, 0)
        expect_identical(digest::digest(file = head, algo = "sha256"), head_hash)
        receipt$oracle <- oracle_id
        receipt$exact <- identical(hashes, expected$output_hashes)
        receipt$status <- "complete"
    }, error = function(error) {
        receipt$status <<- "failed"
        receipt$error <<- conditionMessage(error)
        flush()
        stop(error)
    })
    flush()
})

test_that("P14 writer attribution observes completion and failure", {
    mode <- Sys.getenv("MULTISCHOLAR_DIA025_MECHANISM_MEMBER")
    skip_if(!mode %in% c("gui-bulk-rollup:writer-opt-preflight",
        "gui-bulk-rollup:writer-opt-diagnostic:primary",
        "gui-bulk-rollup:writer-opt-diagnostic:contrast"), "explicit writer attribution only")
    root <- dirname(Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT"))
    proof <- dirname(root)
    withr::local_envvar(c(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
        DUCKPLYR_FALLBACK_AUTOUPLOAD = "0", DUCKPLYR_FALLBACK_COLLECT = "0"))
    diagnostic <- makeDIAWriterDiagnostic()
    testthat::local_mocked_bindings(appendWorkflowTableWriter = diagnostic$append,
        updateWorkflowFastWriter = diagnostic$update)
    receipt <- list(status = "running", mode = mode, profiling = TRUE, phases = list(), failures = list())
    receipt$source_hashes <- as.list(stats::setNames(vapply(c(
        "tests/testthat/helper-workflow-bulk-rollup.R",
        "tests/testthat/test-workflow-rollup-qualification.R",
        "tools/profiling/fastduckdb-writer-diagnostic.patch"),
        \(path) digest::digest(file = here::here(path), algo = "sha256"), character(1)),
        c("bulk_helper", "qualification_test", "native_diagnostics")))
    receipt_path <- file.path(root, "writer-attribution.json")
    flush <- function() {
        receipt$writer <- as.list(diagnostic$observations)
        jsonlite::write_json(receipt, receipt_path, auto_unbox = TRUE, pretty = TRUE,
            null = "null", digits = NA)
    }
    on.exit(flush(), add = TRUE)
    phase <- function(label, fn) {
        receipt$active_phase <<- label
        flush()
        started <- proc.time()[["elapsed"]]
        value <- tryCatch(fn(), error = function(error) {
            receipt$failures[[length(receipt$failures) + 1L]] <<- list(phase = label,
                message = conditionMessage(error))
            stop(error)
        })
        receipt$phases[[length(receipt$phases) + 1L]] <<- list(phase = label,
            seconds = proc.time()[["elapsed"]] - started)
        flush()
        value
    }
    tryCatch({
        if (mode == "gui-bulk-rollup:writer-opt-preflight") {
            table <- data.frame(key = seq_len(10000L), quantity = as.double(seq_len(10000L)),
                annotation = rep(c("ordinary", NA_character_), 5000L))
            metadata <- workflowSavedTableMetadata(table[0, , drop = FALSE], ".rownames", 3L)
            writer <- openWorkflowTableWriter(file.path(root, "completed.parquet"), metadata,
                backend = "fast_r")
            on.exit(abortWorkflowTableWriter(writer), add = TRUE)
            phase("append_and_finish", function() {
                appendWorkflowTableWriter(writer, table)
                invisible(gc())
                appendWorkflowTableWriter(writer, table)
                expect_identical(finishWorkflowTableWriter(writer)$rows, 20000L)
            })
            native <- diagnostic$observations$native[[writer$path]]
            expect_true(native$closed)
            expect_gt(native$conversion_seconds, 0)
            expect_gt(native$concat_calls, 0)
            expect_gt(native$concat_bytes, 0)
            expect_gt(native$flush_seconds, 0)
            failed <- openWorkflowTableWriter(file.path(root, "failed.parquet"), metadata,
                backend = "fast_r")
            on.exit(abortWorkflowTableWriter(failed), add = TRUE)
            appendWorkflowTableWriter(failed, table)
            temporary <- failed$write_path
            table$quantity <- as.character(table$quantity)
            expect_error(phase("injected_schema_failure", function()
                appendWorkflowTableWriter(failed, table)), "schema differs")
            abortWorkflowTableWriter(failed)
            expect_false(file.exists(temporary))
            expect_false(file.exists(failed$path))
            expect_identical(diagnostic$observations$completed, 3L)
            expect_identical(diagnostic$observations$failed, 1L)
            expect_length(receipt$failures, 1L)
        } else {
            case <- if (endsWith(mode, ":primary")) "primary" else "contrast"
            if (case == "primary") {
                source <- jsonlite::read_json(file.path(proof, "browser-smoke.BpSyxB", "gui-state-latest.json"))
            } else {
                manifests <- lapply(Sys.glob(file.path(proof, "browser-smoke.*",
                    "rollup-qualification-inputs.json")), jsonlite::read_json)
                valid <- which(vapply(manifests, \(x) identical(x$status, "complete") &&
                    identical(x$revision, "p3-adaptations-v1"), logical(1)))
                stopifnot(length(valid) == 1L)
                source <- manifests[[valid]]$cases$contrast
            }
            head <- file.path(source$store_directory, "dia-state-head.rds")
            receipt$source_head_sha256 <- digest::digest(file = head, algo = "sha256")
            stopifnot(identical(receipt$source_head_sha256, if (case == "primary")
                "5d64ca562669d94eea1677a1a03d5b410a8a81d8ea59add8cabd6065deb91353" else source$head_sha256))
            owned <- file.path(root, "arm-P14")
            dir.create(owned)
            manager <- copyDIAReplayCheckpoint(source$store_directory, owned)
            on.exit(closeDIAStateAdapter(manager), add = TRUE)
            selectDIAStateReference(manager, "protein_peptide_filtered")
            checkpoint <- readDIAStateMetadata(manager)$data
            stopifnot(checkpoint$dimensions$rows == if (case == "primary") 24000000L else 3840000L)
            design <- withDIAStateScope(manager, function()
                readDIAExecutionComponents(manager, fields = "design_matrix"))@design_matrix
            samples <- as.character(design[[checkpoint$slots$sample_id]][if (case == "primary")
                c(1:4, 149:152, 297:300) else c(1:2, 24:25, 47:48)])
            receipt$samples <- samples
            result <- new.env(parent = emptyenv())
            profile <- profvis::profvis({
                result$value <- runDIABulkOperation(manager, owned, proof, "P14", samples, phase)
                phase("teardown", function() closeDIAStateAdapter(manager))
            }, interval = 0.01, prof_output = file.path(root, "writer.Rprof"), rerun = FALSE)
            value <- result$value
            if (is.null(value)) stop(if (length(receipt$failures))
                receipt$failures[[length(receipt$failures)]]$message else "profiled operation did not complete")
            saveRDS(profile, file.path(root, "writer-profvis.rds"), compress = FALSE)
            receipt$result <- value
            expect_identical(value$summary$input_precursor_rows,
                if (case == "primary") 960000L else 480000L)
            expect_gt(diagnostic$observations$completed, 0L)
            expect_identical(diagnostic$observations$failed, 0L)
            expect_identical(digest::digest(file = head, algo = "sha256"), receipt$source_head_sha256)
        }
        receipt$status <- "complete"
    }, error = function(error) {
        receipt$status <<- "failed"
        receipt$error <<- conditionMessage(error)
        stop(error)
    })
    flush()
})

test_that("production rollup qualification uses frozen prepared cases", {
    mode <- Sys.getenv("MULTISCHOLAR_DIA025_MECHANISM_MEMBER")
    skip_if(!grepl("^gui-bulk-rollup:(qualification-preflight|qualification-prepare|crossover:(B|P0):(4608|9600|96000|960000):[12]|contrast:(B|P([0-9]|1[0-4])))$", mode),
        "explicit qualification member only")
    root <- dirname(Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT"))
    proof <- dirname(root)
    withr::local_envvar(c(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
        DUCKPLYR_FALLBACK_AUTOUPLOAD = "0", DUCKPLYR_FALLBACK_COLLECT = "0"))
    if (mode == "gui-bulk-rollup:qualification-preflight") {
        value <- repairDIARollupQualification(proof, root, preflight = TRUE)
        expect_identical(value$status, "preflight_complete")
    } else if (mode == "gui-bulk-rollup:qualification-prepare") {
        value <- repairDIARollupQualification(proof, root)
        expect_identical(value$status, "complete")
    } else {
        fields <- strsplit(mode, ":", fixed = TRUE)[[1L]]
        kind <- fields[[2L]]
        arm <- fields[[3L]]
        case <- if (kind == "contrast") "contrast" else fields[[4L]]
        paths <- Sys.glob(file.path(proof, "browser-smoke.*", "rollup-qualification-inputs.json"))
        candidates <- lapply(paths, jsonlite::read_json)
        valid <- which(vapply(candidates, function(x) identical(x$status, "complete") &&
            identical(x$revision, "p3-adaptations-v1"), logical(1)))
        stopifnot(length(valid) == 1L)
        preparation <- candidates[[valid]]
        selected <- preparation$cases[[case]]
        source_head <- file.path(selected$store_directory, "dia-state-head.rds")
        stopifnot(identical(digest::digest(file = source_head, algo = "sha256"), selected$head_sha256))
        owned <- file.path(root, paste0("arm-", arm))
        dir.create(owned)
        manager <- copyDIAReplayCheckpoint(selected$store_directory, owned)
        on.exit(closeDIAStateAdapter(manager), add = TRUE)
        receipt <- list(status = "running", mode = mode, case = case, arm = arm,
            source = selected, phases = list(), failures = list(), profiling = FALSE,
            runtime_rdb_sha256 = digest::digest(file = file.path(Sys.getenv("MULTISCHOLAR_DIA008_LIBRARY"),
                "MultiScholaR", "R", "MultiScholaR.rdb"), algo = "sha256"),
            extension_sha256 = diaRollupNativeInfo()$manifest$binary_sha256)
        path <- file.path(root, "rollup-qualification.json")
        flush <- function() jsonlite::write_json(receipt, path, pretty = TRUE, auto_unbox = TRUE,
            null = "null", digits = NA)
        on.exit(flush(), add = TRUE)
        phase <- function(label, fn) {
            started <- proc.time()
            receipt$active_phase <<- label
            flush()
            value <- fn()
            receipt$phases[[length(receipt$phases) + 1L]] <<- list(phase = label,
                seconds = unname((proc.time() - started)[["elapsed"]]),
                cpu = unname(sum((proc.time() - started)[c("user.self", "sys.self")])))
            flush()
            value
        }
        tryCatch({
            started <- proc.time()[["elapsed"]]
            value <- runDIABulkOperation(manager, owned, proof, arm, NULL, phase)
            phase("teardown", function() closeDIAStateAdapter(manager))
            receipt$result <- c(value, list(complete_operation_seconds = proc.time()[["elapsed"]] - started))
            identity <- readSavedStoreIdentity(file.path(owned, ".workflow-entities", "store.json"))
            check_store <- newWorkflowEntityStore(owned, identity$project_id, identity$workflow_id)
            on.exit(closeWorkflowEntityStore(check_store), add = TRUE)
            scope <- newWorkflowTableReadScope(check_store, value$reference, reader = "arrow")
            on.exit(releaseWorkflowTableReadScope(scope), add = TRUE, after = FALSE)
            hashes <- character()
            streamSavedWorkflowTable(scope, scope$metadata$columns, function(data, offset, row_names) {
                attr(data, "dia_source_order") <- NULL
                attr(data, "row.names") <- .set_row_names(nrow(data))
                hashes[[length(hashes) + 1L]] <<- digest::digest(data, algo = "sha256")
            })
            releaseWorkflowTableReadScope(scope)
            closeWorkflowEntityStore(check_store)
            receipt$result$output_hashes <- hashes
            stopifnot(identical(digest::digest(file = source_head, algo = "sha256"), selected$head_sha256))
            receipt$status <- "complete"
        }, error = function(error) {
            receipt$status <<- "failed"
            receipt$error <<- conditionMessage(error)
            flush()
            stop(error)
        })
        flush()
    }
})

test_that("bounded eight-thread reader qualification preserves the existing codec", {
    mode <- Sys.getenv("MULTISCHOLAR_DIA025_MECHANISM_MEMBER")
    skip_if(!grepl("^gui-bulk-rollup:(io-preflight|io-prepare|io:(duckdb|fastduckdb):[12])$", mode),
        "explicit bounded reader qualification only")
    root <- dirname(Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT"))
    proof <- dirname(root)
    receipt <- list(status = "running", mode = mode, results = list(), failures = list())
    path <- file.path(root, "rollup-io.json")
    flush <- function() jsonlite::write_json(receipt, path, pretty = TRUE, auto_unbox = TRUE,
        null = "null", digits = NA)
    on.exit(flush(), add = TRUE)
    tryCatch({
        if (mode == "gui-bulk-rollup:io-prepare") {
            input <- prepareDIARollupIO(proof, root)
            receipt$input <- list(index = input$index$entity_id, wide = input$wide$entity_id)
        } else {
            if (mode == "gui-bulk-rollup:io-preflight") {
                expect_false(identical(fingerprintDIARollupIO(data.frame(value = c(NA_real_, NaN))),
                    fingerprintDIARollupIO(data.frame(value = c(NA_real_, NA_real_)))))
                owned <- file.path(root, "fixture")
                dir.create(owned)
                store <- newWorkflowEntityStore(owned, "io-preflight", "DIA")
                on.exit(closeWorkflowEntityStore(store), add = TRUE)
                table <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))@peptide_data[1:30, ]
                table$Run <- factor(table$Run, levels = c(unique(table$Run), "unused"))
                rownames(table) <- NULL
                index <- data.frame(source_order = seq_len(30L))
                identity <- list(input = "io-preflight", design = "fixed", parameters = "eight-threads",
                    selectors = "bounded", software = "bvT5nP")
                input <- list(root = owned, identity = store$identity,
                    index = saveWorkflowTable(store, index, "index", identity, ".rownames", metadata_version = 3L),
                    wide = saveWorkflowTable(store, table, "wide", identity, ".rownames", metadata_version = 3L),
                    requests = list(index = list(offset = c(0L, 10L, 20L), limit = c(5L, 5L, 5L)),
                        wide = list(offset = c(0L, 10L, 20L), limit = c(5L, 5L, 5L))))
                input$expected <- lapply(list(index = index, wide = table), function(data)
                    vapply(c(0L, 10L, 20L), function(offset) {
                        value <- data[offset + seq_len(5L), , drop = FALSE]
                        fingerprintDIARollupIO(restoreWorkflowAttributeOrder(value, names(attributes(data))))
                    }, character(1)))
                closeWorkflowEntityStore(store)
                readers <- c("duckdb", "fastduckdb")
            } else {
                inputs <- Sys.glob(file.path(proof, "browser-smoke.*", "rollup-io-inputs-v2.rds"))
                stopifnot(length(inputs) == 1L)
                input <- readRDS(inputs[[1L]])
                stopifnot(identical(input$revision, "codec-attributes-v2"))
                readers <- strsplit(mode, ":", fixed = TRUE)[[1L]][[3L]]
                receipt$input_sha256 <- digest::digest(file = inputs[[1L]], algo = "sha256")
            }
            started <- proc.time()[["elapsed"]]
            for (reader in readers) for (workload in c("index", "wide")) {
                key <- paste(reader, workload, sep = ":")
                receipt$active <- key
                flush()
                receipt$results[[key]] <- runDIARollupIO(input, workload, reader, root)
                if (mode == "gui-bulk-rollup:io-preflight") {
                    expect_error(runDIARollupIO(input, workload, reader, root, failure = TRUE),
                        "injected bounded reader failure")
                    cleanup <- jsonlite::read_json(file.path(root,
                        paste0("io-cleanup-", reader, "-", workload, "-TRUE.json")))
                    expect_identical(cleanup$open_resources, 0L)
                    receipt$failures[[key]] <- list(injected = TRUE, cleanup = cleanup)
                }
                flush()
            }
            receipt$complete_operation_seconds <- proc.time()[["elapsed"]] - started
            if (mode == "gui-bulk-rollup:io-preflight") {
                real <- repairDIARollupIOInputs(proof, root)
                # Check every frozen request before admitting any timed member.
                for (reader in readers) for (workload in c("index", "wide")) {
                    key <- paste("prepared", reader, workload, sep = ":")
                    receipt$active <- key
                    flush()
                    receipt$results[[key]] <- runDIARollupIO(real, workload, reader, root)
                }
                saveRDS(real, file.path(root, "rollup-io-inputs-v2.rds"), compress = FALSE)
                receipt$verified_prepared_requests <- TRUE
            }
        }
        receipt$status <- "complete"
    }, error = function(error) {
        receipt$status <<- "failed"
        receipt$error <<- conditionMessage(error)
        stop(error)
    })
    flush()
})


test_that("rollup scan inspection preserves immutable source inputs", {
    skip_if(Sys.getenv("MULTISCHOLAR_DIA025_MECHANISM_MEMBER") != "gui-bulk-rollup:scan-inspection",
        "explicit untimed physical-plan inspection only")
    root <- dirname(Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT"))
    value <- inspectDIARollupScans(dirname(root), root)
    expect_identical(value$status, "complete")
})

test_that("the direct native writer build is isolated and source-bound", {
    skip_if(!Sys.getenv("MULTISCHOLAR_DIA025_MECHANISM_MEMBER") %in% c("gui-bulk-rollup:writer-build", "gui-bulk-rollup:writer-opt-build"),
        "explicit isolated native build only")
    root <- dirname(Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT"))
    proof <- dirname(root)
    request_name <- if (Sys.getenv("MULTISCHOLAR_DIA025_MECHANISM_MEMBER") == "gui-bulk-rollup:writer-opt-build") "fast-writer-opt-source.json" else "fast-writer-source.json"
    request <- jsonlite::read_json(file.path(proof, request_name))
    sources <- names(request$files)
    stopifnot(all(vapply(sources, function(path) identical(
        digest::digest(file = file.path(request$source, path), algo = "sha256"), request$files[[path]]), logical(1))))
    destination <- file.path(root, "source")
    dir.create(destination)
    for (path in list.files(request$source, all.files = TRUE, no.. = TRUE, full.names = TRUE))
        stopifnot(file.copy(path, destination, recursive = TRUE))
    package <- file.path(destination, "fastduckdb")
    cache <- file.path(package, "src", ".cargo")
    dir.create(cache)
    for (path in list.files(request$cargo_cache, all.files = TRUE, no.. = TRUE, full.names = TRUE))
        stopifnot(file.copy(path, cache, recursive = TRUE))
    library <- file.path(root, "library")
    dir.create(library)
    withr::local_envvar(c(CARGO_BUILD_JOBS = "8", CARGO_NET_OFFLINE = "true"))
    status <- system2(file.path(R.home("bin"), "R"), c("CMD", "INSTALL",
        shQuote(paste0("--library=", library)), shQuote(package)))
    stopifnot(status == 0L)
    jsonlite::write_json(list(status = "complete", source = request,
        library = library, source_copy = destination, threads = 8L,
        version = read.dcf(file.path(library, "fastduckdb", "DESCRIPTION"))[1L, "Version"],
        binaries = stats::setNames(vapply(list.files(file.path(library, "fastduckdb"),
            "[.](so|rdb|rdx)$", recursive = TRUE, full.names = TRUE),
            function(path) digest::digest(file = path, algo = "sha256"), character(1)),
            list.files(file.path(library, "fastduckdb"), "[.](so|rdb|rdx)$", recursive = TRUE))),
        file.path(root, "fast-writer-build.json"), auto_unbox = TRUE, pretty = TRUE)
})
