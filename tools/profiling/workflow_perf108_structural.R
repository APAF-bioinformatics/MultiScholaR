# PERF-108 fixed component probes. Diagnostic evidence only, never a full-operation denominator.

#' Expand the fixed screen or full confirmation schedule.
#' @param bindings Frozen installed libraries, inputs and original scientific oracles.
#' @param stage Structural screen, rollup, primary or wide stage.
#' @return Fixed ordered cases; O is enabled only for the rollup candidate.
#' @noRd
perf108StructuralCases <- function(bindings, stage) {
    kind <- sub("^ordered-", "structural-", stage)
    stopifnot(kind %in% c("structural-screen", "structural-rollup", "structural-primary", "structural-wide"))
    screen <- kind == "structural-screen"
    population <- switch(kind, `structural-primary` = "primary", `structural-wide` = "wide", "p14")
    order <- if (screen) c("current", "candidate", "candidate", "current") else
        c("accepted", "current", "candidate", "candidate", "current", "accepted")
    repetitions <- c(accepted = 0L, current = 0L, candidate = 0L)
    input <- bindings[[if (screen) "screen" else "full"]][[population]]
    lapply(seq_along(order), function(i) {
        role <- order[[i]]
        repetitions[[role]] <<- repetitions[[role]] + 1L
        rep <- repetitions[[role]]
        options <- if (role == "candidate" && population == "p14") {
            if (is.null(bindings$worker_options)) list(multischolar.perf108.native_output = TRUE) else bindings$worker_options
        } else NULL
        list(id = paste(population, role, rep, sep = "-"), role = role, rep = rep,
            variant = if (population == "p14") "p14" else "physical", population = population,
            control = role == "accepted", fault = "none", profile_worker = FALSE,
            structural_options = options, library = bindings[[role]]$library,
            expected_build = bindings[[role]]$build, input = input, verification_slot = i,
            binding = digest::digest(list(input, bindings[[role]], stage, rep, options,
                protocol = bindings$protocol_sha256), algo = "sha256"))
    })
}

#' Apply prospective full-operation gates without overriding historical holds.
#' @param receipt,cases Completed stage receipt and fixed cases.
#' @return Separate science, variability, benefit and runtime-eligibility decisions.
#' @noRd
perf108StructuralDecision <- function(receipt, cases) {
    screen <- receipt$stage %in% c("structural-screen", "ordered-screen")
    valid <- vapply(cases, function(case) {
        item <- receipt$cases[[case$id]]
        op <- item$operation
        isTRUE(item$passed) && isTRUE(item$verifier_completion$accepted) &&
            identical(op$status, "complete") && isTRUE(op$measurement_valid) &&
            isTRUE(op$lifecycle_passed) && identical(op$swap_bytes, 0) &&
            is.finite(op$peak_bytes) && op$peak_bytes > 0 && op$peak_bytes < 15.5 * 1024^3 &&
            is.finite(op$elapsed_seconds) && op$elapsed_seconds > 0 &&
            op$elapsed_seconds <= if (screen) 300 else 1800
    }, logical(1))
    result <- list(status = "invalid_or_failed_member", promotion_authorized = FALSE,
        science_outputs_and_measurements_valid = all(valid), valid_members = valid,
        historical_primary_hold_unchanged = TRUE)
    if (!all(valid)) return(result)
    if (screen) { result$status <- "screen_valid"; return(result) }
    operation <- function(role, rep) {
        case <- Filter(function(x) x$role == role && x$rep == rep, cases)[[1L]]
        receipt$cases[[case$id]]$operation
    }
    result$reference_variability <- lapply(c("accepted", "current"), function(role) {
        time <- vapply(1:2, function(rep) operation(role, rep)$elapsed_seconds, numeric(1))
        peak <- vapply(1:2, function(rep) operation(role, rep)$peak_bytes, numeric(1))
        list(role = role, elapsed = max(time) / min(time) - 1, peak = max(peak) / min(peak) - 1)
    })
    comparisons <- lapply(c("accepted", "current"), function(role) {
        pairs <- lapply(1:2, function(rep) {
            a <- operation("candidate", rep); b <- operation(role, rep)
            c(elapsed = a$elapsed_seconds / b$elapsed_seconds, peak = a$peak_bytes / b$peak_bytes)
        })
        time <- vapply(pairs, `[[`, numeric(1), "elapsed")
        peak <- vapply(pairs, `[[`, numeric(1), "peak")
        list(reference = role, pairs = pairs,
            memory_winner = all(peak <= 0.8 & time <= 1.15),
            time_winner = all(time <= 0.9 & peak <= 1.05),
            non_regression = all(time <= 1.05 & peak <= 1.05))
    })
    names(comparisons) <- c("accepted", "current")
    result$comparisons <- comparisons
    if (any(vapply(result$reference_variability, function(x)
        x$elapsed > 0.15 || x$peak > 0.10, logical(1)))) {
        result$status <- "inconclusive_reference_variability"
        return(result)
    }
    accepted <- comparisons$accepted; current <- comparisons$current
    result$runtime_eligible <- accepted$non_regression || accepted$memory_winner
    result$meaningful_benefit <- current$memory_winner || current$time_winner
    result$status <- if (receipt$stage %in% c("structural-rollup", "ordered-rollup")) {
        if (result$runtime_eligible && result$meaningful_benefit) "rollup_candidate_supported" else "no_justified_rollup_gain"
    } else if (result$runtime_eligible) "consumer_non_regression_passed" else "consumer_regression"
    result
}

#' Prepare one immutable representative cohort from accepted artifacts.
#' @param run_dir,repo Existing guarded member and its source snapshot.
#' @return Frozen fixture manifest; original payload bytes and heads remain unchanged.
#' @noRd
preparePerf108StructuralFixture <- function(run_dir, repo) {
    bindings <- readRDS(file.path(repo, "docs/qa/workflow-execution/perf-108/cache-qualification/bindings.rds"))
    manifest <- file.path(dirname(run_dir), "structural-fixture.rds")
    stopifnot(!file.exists(manifest))
    helpers <- new.env(parent = asNamespace("MultiScholaR"))
    sys.source(file.path(repo, "tests/testthat/helper-workflow-bulk-rollup.R"), helpers)
    project <- file.path(run_dir, "project")
    stopifnot(dir.create(project))
    manager <- helpers$copyDIAReplayCheckpoint(bindings$full$p14$store, project)
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    record <- readDIAStateMetadata(manager, bindings$full$p14$state)$data
    object <- readDIAExecutionComponents(manager, bindings$full$p14$state, fields = "design_matrix")
    roles <- list(sample = object@sample_id, protein = object@protein_id_column,
        peptide = object@peptide_sequence_column, modified = "Modified.Sequence", precursor = "Precursor.Id",
        raw = object@raw_quantity_column, normalised = object@norm_quantity_column)
    columns <- names(object@peptide_data)
    input <- list(reference = record$components$peptide_data, roles = roles,
        evidence = intersect(c("identification_peptide_count", "identification_peptidoform_count",
            "peptides_for_protein_count", "peptidoforms_for_protein_count"), columns),
        annotation = if (roles$protein != "Protein.Ids" && "Protein.Ids" %in% columns) "Protein.Ids" else character())
    stopifnot(identical(input$reference, bindings$full$p14$reference))
    accepted <- readRDS(file.path(bindings$root, "cache-rollup/p14-accepted-1/member-result.rds"))
    original <- file.path(accepted$project_root, ".workflow-entities", accepted$output$reference$payload$relative_path)
    stopifnot(identical(digest::digest(file = original, algo = "sha256"), accepted$output$reference$payload$sha256))
    metadata <- readSavedWorkflowTableMetadata(original)
    metadata["rows"] <- list(NULL)
    encoded <- jsonlite::base64_enc(serialize(metadata, NULL, version = 3))
    path <- file.path(run_dir, "accepted-100-samples.parquet")
    con <- fastduckdb::open_query_fast()
    class(con) <- "workflow_fast_query"
    on.exit(fastduckdb::close_query_fast(con), add = TRUE)
    q <- function(x) as.character(workflowQueryQuoteString(con, x))
    workflowQueryExecute(con, paste0("COPY (SELECT * FROM read_parquet(", q(original),
        ') WHERE ".multischolar_order" <= 4000000 ORDER BY ".multischolar_order") TO ', q(path),
        " (FORMAT PARQUET, COMPRESSION SNAPPY, ROW_GROUP_SIZE 10240, KV_METADATA {multischolar_table: ",
        q(encoded), "})"))
    samples <- workflowQueryCollect(con, paste0("SELECT DISTINCT ",
        workflowQueryQuoteIdentifier(con, roles$sample), " AS sample FROM read_parquet(", q(path), ") ORDER BY sample"))$sample
    stopifnot(identical(samples, bindings$p14_oracle$samples[seq_len(100L)]),
        readSavedWorkflowTableMetadata(path)$rows == 4000000L,
        identical(digest::digest(file = file.path(bindings$full$p14$store, "dia-state-head.rds"),
            algo = "sha256"), bindings$full$p14$head_sha256))
    fixture <- list(project = project, input = input, samples = samples, output_path = path,
        output_sha256 = digest::digest(file = path, algo = "sha256"), original_output_sha256 = accepted$output$reference$payload$sha256,
        input_sha256 = input$reference$payload$sha256, input_rows = 24000000L,
        cohort_input_rows = 8000000L, cohort_output_rows = 4000000L,
        input_bindings_sha256 = digest::digest(bindings, algo = "sha256"))
    fixture$output_fingerprints <- fingerprintPerf108ProbeFile(path, ".multischolar_order", 4000000L,
        metadata = readSavedWorkflowTableMetadata(path))
    saveRDS(fixture, manifest)
    fixture
}

#' Run a fixed diagnostic component comparison member.
#' @param stage One of PA0/PA1, PV0/PV1, PO0/PO1; zero is its reference.
#' @param run_dir,repo Existing guarded member and source snapshot.
#' @return Diagnostic receipt with phase boundaries and bounded output fingerprints.
#' @noRd
runPerf108StructuralProbe <- function(stage, run_dir, repo) {
    stopifnot(stage %in% c("PA0", "PA1", "PV0", "PV1", "PO0", "PO1"))
    fixture <- readRDS(file.path(dirname(run_dir), "structural-fixture.rds"))
    candidate <- endsWith(stage, "1")
    opportunity <- substr(stage, 2L, 2L)
    withr::local_options(multischolar.perf108.paired_state = opportunity == "A" && candidate,
        multischolar.perf108.partition_validation = opportunity == "V" && candidate,
        multischolar.perf108.native_output = opportunity == "O" && candidate,
        multischolar.perf108.persistent_reader = FALSE)
    observer <- new.env(parent = asNamespace("MultiScholaR"))
    sys.source(file.path(repo, "tools/profiling/workflow_perf108.R"), observer)
    finish <- observer$startPerf108Attribution(run_dir)
    on.exit(finish(), add = TRUE)
    phase <- getOption("multischolar.dia_rollup_phase")
    saveRDS(list(epoch = as.numeric(Sys.time()),
        uptime = as.numeric(strsplit(readLines("/proc/uptime", n = 1L), " ")[[1L]][[1L]])),
        file.path(run_dir, "clock-anchor.rds"))
    receipt <- list(status = "failed", stage = stage, library = find.package("MultiScholaR"),
        fixture_sha256 = digest::digest(fixture, algo = "sha256"), threads = 8L,
        measurement_kind = "diagnostic component; not complete operation")
    on.exit(saveRDS(receipt, file.path(run_dir, "component-result.rds")), add = TRUE)
    if (opportunity %in% c("A", "V")) {
        owner <- readSavedStoreIdentity(file.path(fixture$project, ".workflow-entities/store.json"))
        store <- newWorkflowEntityStore(fixture$project, owner$project_id, owner$workflow_id)
        on.exit(closeWorkflowEntityStore(store), add = TRUE)
        scope <- newDIARelationalScope(store, fixture$input, profile = "P14", native = TRUE)
        on.exit(releaseWorkflowTableReadScope(scope), add = TRUE, after = FALSE)
        con <- scope$connection
        rows <- phase("relational_validation", function() validateDIARelationalInput(scope,
            fixture$input, if (opportunity == "A") fixture$samples else NULL))
        stopifnot(rows == if (opportunity == "A") fixture$cohort_input_rows else fixture$input_rows)
        if (opportunity == "V") {
            receipt$admission <- phase("working_memory_admission", function()
                estimateDIARollupWorkingMemory(scope, fixture$input, 16 * 1024^3))
        } else {
            workflowQueryExecute(con, "CREATE TEMP VIEW qc_aggregate_input AS SELECT * FROM qc_input")
            prepareDIARelationalAggregate(scope, fixture$input, "F")
            keys <- unlist(fixture$input$roles[c("sample", "protein", "peptide")])
            order <- paste(workflowQueryQuoteIdentifier(con, keys), collapse = ",")
            query <- paste0('SELECT *, row_number() OVER (ORDER BY ', order,
                ') AS ".probe_order" FROM qc_aggregated ORDER BY ', order)
            path <- file.path(run_dir, "aggregate.parquet")
            phase("ordered_result_stage", function() workflowQueryExecute(con, paste0("COPY (", query,
                ") TO ", workflowQueryQuoteString(con, path), " (FORMAT PARQUET, ROW_GROUP_SIZE 10240)")))
            fastduckdb::close_query_fast(con)
            scope$connection <- scope$action$connection <- NULL
            receipt$fingerprints <- fingerprintPerf108ProbeFile(path, ".probe_order", fixture$cohort_output_rows)
        }
    } else {
        path <- fixture$output_path
        stopifnot(identical(digest::digest(file = path, algo = "sha256"), fixture$output_sha256))
        metadata <- readSavedWorkflowTableMetadata(path)
        schema <- lapply(metadata$column_metadata, function(spec) vector(spec$typeof, 0L))
        schema <- as.data.frame(schema, check.names = FALSE)
        prototype <- restoreSavedWorkflowTable(schema, metadata, metadata$columns)
        metadata["rows"] <- list(NULL)
        destination <- file.path(run_dir, "output.parquet")
        counts <- newDIAFilteringAccumulator(prototype, fixture$input$roles$protein)
        writer <- openWorkflowTableWriter(destination, metadata, backend = "fast_r")
        writer$native_query <- candidate
        on.exit(abortWorkflowTableWriter(writer), add = TRUE)
        appendWorkflowTableWriter(writer, prototype)
        phase("output_read_summary_write", function() {
            if (candidate) {
                con <- fastduckdb::open_query_fast()
                class(con) <- "workflow_fast_query"
                on.exit(fastduckdb::close_query_fast(con))
                fields <- unique(c("Run", "Protein.Ids", "Stripped.Sequence", fixture$input$roles$protein))
                query <- paste0("SELECT * FROM read_parquet(", workflowQueryQuoteString(con, path),
                    ') ORDER BY ".multischolar_order"')
                receipt$calls <<- streamWorkflowQueryToWriter(con, query, writer, fields,
                    function(part) counts$push(restoreSavedWorkflowTable(part, metadata, fields)))
            } else {
                # Match the complete producer boundary: the reference must
                # create its private ordered stage before reading it back.
                con <- fastduckdb::open_query_fast()
                class(con) <- "workflow_fast_query"
                on.exit(fastduckdb::close_query_fast(con), add = TRUE)
                staged <- file.path(run_dir, "ordered-result.parquet")
                on.exit(unlink(staged), add = TRUE)
                workflowQueryExecute(con, paste0("COPY (SELECT *, row_number() OVER (ORDER BY ",
                    '".multischolar_order") AS ".dia_output_order" FROM read_parquet(',
                    workflowQueryQuoteString(con, path), ') ORDER BY ".multischolar_order") TO ',
                    workflowQueryQuoteString(con, staged),
                    " (FORMAT PARQUET, COMPRESSION SNAPPY, ROW_GROUP_SIZE 10240)"))
                receipt$intermediate_bytes_written <<- file.size(staged)
                fastduckdb::close_query_fast(con)
                fastduckdb::clear_parquet_cache()
                on.exit(fastduckdb::clear_parquet_cache(), add = TRUE)
                for (start in seq.int(1L, fixture$cohort_output_rows, by = 10000L)) {
                    part <- fastduckdb::read_parquet_fast(staged, n_threads = 8L, lazy_strings = FALSE,
                        filter = list(.dia_output_order = c(start, start + 9999L)))
                    stopifnot(identical(as.integer(part$.multischolar_order), start:(start + 9999L)))
                    data <- restoreSavedWorkflowTable(part, metadata, metadata$columns)
                    counts$push(data)
                    appendWorkflowTableWriter(writer, data)
                }
                receipt$calls <<- 400L
                unlink(staged)
            }
            finishWorkflowTableWriter(writer)
        })
        receipt$counts <- counts$finish()
        receipt$fingerprints <- fingerprintPerf108ProbeFile(destination, ".multischolar_order", fixture$cohort_output_rows,
            metadata = readSavedWorkflowTableMetadata(destination))
        stopifnot(identical(receipt$fingerprints, fixture$output_fingerprints))
    }
    receipt$status <- "complete"
    receipt
}

#' Fingerprint fixed logical pages outside the observed component phase.
#' @param path,ordinal Output file and dense ordering field.
#' @param rows Expected complete population.
#' @param metadata Optional existing restoration contract for logical fingerprints.
#' @return One SHA-256 per 10,000 rows, with bounded R reads.
#' @noRd
fingerprintPerf108ProbeFile <- function(path, ordinal, rows, metadata = NULL) {
    fastduckdb::clear_parquet_cache()
    on.exit(fastduckdb::clear_parquet_cache())
    vapply(seq.int(1L, rows, by = 10000L), function(start) {
        part <- fastduckdb::read_parquet_fast(path, n_threads = 8L, lazy_strings = FALSE,
            filter = stats::setNames(list(c(start, min(rows, start + 9999L))), ordinal))
        stopifnot(identical(as.integer(part[[ordinal]]), seq.int(start, min(rows, start + 9999L))))
        if (!is.null(metadata)) part <- restoreSavedWorkflowTable(part, metadata, metadata$columns)
        rownames(part) <- NULL
        digest::digest(part, algo = "sha256")
    }, character(1))
}
