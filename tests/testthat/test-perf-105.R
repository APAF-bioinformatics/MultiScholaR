testthat::skip_if_not(identical(Sys.info()[["sysname"]], "Linux"))

#' Execute the unchanged native oracle and qualified P14 boundary on one fixture.
#' @param object Linear peptide object.
#' @param worker_cap Explicit admitted process bytes.
#' @return Native and persisted results plus their owned store.
#' @noRd
runPerf105Rollup <- function(object, worker_cap = 12 * 1024^3) {
    artifact <- tryCatch(diaRollupExtensionInfo(), error = function(error) NULL)
    testthat::skip_if(is.null(artifact), "qualified exact-sum runtime unavailable")
    roles <- list(sample = object@sample_id, protein = object@protein_id_column,
        peptide = object@peptide_sequence_column, modified = "Modified.Sequence",
        precursor = "Precursor.Id", raw = object@raw_quantity_column,
        normalised = object@norm_quantity_column)
    native <- rollUpPrecursorToPeptideHelper(object@peptide_data,
        sample_id_column = roles$sample, protein_id_column = roles$protein,
        peptide_sequence_column = roles$peptide,
        modified_peptide_sequence_column = roles$modified,
        precursor_quantity_column = roles$raw,
        precursor_normalised_column = roles$normalised,
        precursor_id_column = roles$precursor, core_utilisation = NA,
        is_logged_data = FALSE, return_rollup_result = TRUE)
    root <- withr::local_tempdir()
    store <- newWorkflowEntityStore(root, basename(root), "DIA")
    identity <- newSavedWorkflowScientificIdentity(
        digest::digest(object@peptide_data, algo = "sha256"), "design",
        "linear-sum", "complete-samples", "perf105")
    input <- saveWorkflowPrecursorRollupInput(store, object, "input", identity)
    result <- computeDIARelationalRollup(store, input,
        list(prototype = native$data, identity = identity), profile = "P14",
        worker_cap = worker_cap)
    actual <- materializeDIATable(list(store = store), result$reference)
    list(store = store, input = input, native = native, result = result, actual = actual)
}

test_that("bulk query helpers retain fixed engine and metadata limits", {
    store <- newWorkflowEntityStore(withr::local_tempdir(), "perf105-query", "DIA")
    on.exit(closeWorkflowEntityStore(store), add = TRUE)
    table <- data.frame(id = letters[1:4], quantity = as.double(1:4))
    reference <- saveWorkflowTable(store, table, "input",
        newSavedWorkflowScientificIdentity("input", "design", "query", "all", "perf105"),
        "id", metadata_version = 3L)
    scope <- newWorkflowTableReadScope(store, reference, reader = "duckdb")
    on.exit(releaseWorkflowTableReadScope(scope), add = TRUE)
    expect_identical(configureWorkflowBulkQueryScope(scope),
        list(threads = 8L, engine_limit = "12GiB", spill = FALSE))
    settings <- collectWorkflowBulkMetadata(scope,
        "SELECT current_setting('threads') AS threads, current_setting('max_temp_directory_size') AS spill",
        max_rows = 1L, max_bytes = 65536)
    expect_identical(as.integer(settings$threads), 8L)
    expect_identical(settings$spill, "0 bytes")
    expect_identical(collectWorkflowBulkMetadata(scope,
        "SELECT CAST(i AS INTEGER) AS i FROM range(4) input(i) ORDER BY i", max_rows = 4L)$i,
        0:3)
    expect_error(collectWorkflowBulkMetadata(scope,
        "SELECT CAST(i AS INTEGER) AS i FROM range(5) input(i)", max_rows = 4L),
        class = "saved_entity_limit")
    expect_identical(DBI::dbGetQuery(scope$connection, "SELECT 1 AS ok")$ok, 1L)
    releaseWorkflowTableReadScope(scope)
    expect_identical(store$open_resources, 0L)
})

test_that("complete scientific groups split stably without fragmentation", {
    groups <- factor(sprintf("S%02d", seq_len(10L)),
        levels = c(sprintf("S%02d", seq_len(10L)), "unused"), ordered = TRUE)
    batches <- partitionWorkflowCompleteGroups(groups, 3L)
    expect_identical(as.integer(lengths(batches)), c(3L, 3L, 4L))
    expect_identical(do.call(c, batches), groups)
    expect_identical(lapply(batches, levels), rep(list(levels(groups)), 3L))
    expect_identical(partitionWorkflowCompleteGroups(character(), 3L), list(character()))
    expect_error(partitionWorkflowCompleteGroups(c("A", "A"), 3L),
        class = "saved_entity_limit")
    expect_error(partitionWorkflowCompleteGroups(c("A", NA_character_), 3L),
        class = "saved_entity_limit")
    expect_error(partitionWorkflowCompleteGroups(sprintf("S%05d", seq_len(10001L)), 3L),
        class = "saved_entity_limit")
})

test_that("P14 preserves native order sums missingness annotations and three batches", {
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    table <- object@peptide_data
    table$identification_peptide_count <- 3L
    table$identification_peptidoform_count <- 6L
    table$peptides_for_protein_count <- 3L
    table$peptidoforms_for_protein_count <- 6L
    key <- paste(table$Run, table$Protein.Group, table$Stripped.Sequence)
    groups <- unique(key)
    table$Precursor.Quantity[key == groups[[1L]]] <- NA_real_
    table$Precursor.Normalised[which(key == groups[[2L]])[[1L]]] <- NA_real_
    table$Precursor.Quantity[which(key == groups[[3L]])[[1L]]] <- 0
    table$Protein.Ids[1:3] <- c("P1; P2", " P2;P3;", "P1; P2")
    table$Run <- factor(table$Run, levels = c(rev(unique(table$Run)), "unused"), ordered = TRUE)
    object@peptide_data <- table
    value <- runPerf105Rollup(object)
    on.exit(closeWorkflowEntityStore(value$store), add = TRUE)
    expect_identical(value$actual, value$native$data)
    expect_identical(value$result$summary, value$native$summary)
    expect_identical(value$result$execution$independence_unit, "complete_sample")
    expect_identical(value$result$execution$global_validation,
        "completed_before_batching")
    expect_identical(length(value$result$execution$batch_sizes), 3L)
    expect_identical(sum(value$result$execution$batch_sizes),
        length(unique(table$Run)))
    expect_identical(value$result$execution$writer_backend, "fast_r")
    expect_identical(value$result$execution$writer_threads, 8L)
    expect_identical(value$result$execution$engine_threads, 8L)
    expect_gt(value$result$summary$raw_all_missing_outputs, 0L)
    expect_gt(value$result$summary$normalised_partial_missing_outputs, 0L)
    expect_identical(value$store$open_resources, 0L)
})

test_that("P14 retains the ordered native threshold and groups crossing I/O batches", {
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    table <- object@peptide_data[rep(1L, 10001L), , drop = FALSE]
    table$Precursor.Id <- sprintf("boundary-%05d", seq_len(nrow(table)))
    table$Precursor.Quantity <- table$Precursor.Normalised <- 1
    table$Precursor.Quantity[1:3] <- table$Precursor.Normalised[1:3] <- c(1e16, 1, 1)
    object@peptide_data <- table
    value <- runPerf105Rollup(object)
    on.exit(closeWorkflowEntityStore(value$store), add = TRUE)
    expect_identical(value$actual, value$native$data)
    expect_identical(value$result$summary, value$native$summary)
    expect_identical(value$result$execution$max_group_rows, 10001L)
    expect_identical(value$actual$Peptide.RawQuantity,
        value$native$data$Peptide.RawQuantity)
    expect_identical(value$store$open_resources, 0L)
})

test_that("P14 globally rejects invalid science before batch aggregation", {
    base <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    cases <- list(
        integer_quantity = function(x) { x$Precursor.Quantity <- as.integer(x$Precursor.Quantity); x },
        missing_key = function(x) { x$Precursor.Id[[1L]] <- NA_character_; x },
        blank_key = function(x) { x$Stripped.Sequence[[1L]] <- " "; x },
        duplicate = function(x) rbind(x, x[1L, ]),
        negative = function(x) { x$Precursor.Quantity[[1L]] <- -1; x },
        nonfinite = function(x) { x$Precursor.Normalised[[1L]] <- Inf; x },
        evidence = function(x) {
            x$identification_peptide_count <- 3L
            x$identification_peptide_count[[1L]] <- 4L
            x
        })
    patterns <- c(integer_quantity = "double-precision", missing_key = "missing or blank",
        blank_key = "missing or blank", duplicate = "duplicate precursor",
        negative = "negative or non-finite", nonfinite = "negative or non-finite",
        evidence = "conflicting frozen")
    for (name in names(cases)) {
        object <- base
        object@peptide_data <- cases[[name]](base@peptide_data)
        store <- newWorkflowEntityStore(withr::local_tempdir(), paste0("invalid-", name), "DIA")
        identity <- newSavedWorkflowScientificIdentity(name, "design", "invalid", "all", "perf105")
        input <- saveWorkflowPrecursorRollupInput(store, object, "input", identity)
        phases <- character()
        expect_error(computeDIARelationalRollup(store, input,
            list(prototype = base@peptide_data[0, ], identity = identity), profile = "P14",
            phase = function(label, fn) { phases <<- c(phases, label); fn() }), patterns[[name]])
        expect_identical(phases, "relational_validation")
        expect_identical(store$open_resources, 0L)
        closeWorkflowEntityStore(store)
    }
    native <- base
    native@peptide_data$Precursor.Quantity <- as.integer(native@peptide_data$Precursor.Quantity)
    expect_no_error(rollUpPrecursorToPeptide(native, core_utilisation = NA))
})

test_that("P14 refuses one indivisible oversized group before aggregation", {
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    table <- object@peptide_data[rep(1L, 2000L), , drop = FALSE]
    table$Precursor.Id <- paste0("large-group-", seq_len(nrow(table)))
    table$Modified.Sequence <- rep(strrep("A", 40000L), nrow(table))
    object@peptide_data <- table
    store <- newWorkflowEntityStore(withr::local_tempdir(), "oversized", "DIA")
    identity <- newSavedWorkflowScientificIdentity("oversized", "design", "group", "all", "perf105")
    input <- saveWorkflowPrecursorRollupInput(store, object, "input", identity)
    expect_error(computeDIARelationalRollup(store, input,
        list(prototype = object@peptide_data[0, ], identity = identity), profile = "P14",
        worker_cap = 12 * 1024^3), "indivisible rollup group")
    expect_identical(store$open_resources, 0L)
    closeWorkflowEntityStore(store)
})

test_that("sample filtering uses one global bulk population and matches native", {
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    object@peptide_data$Q.Value <- 0
    object@q_value_column <- "Q.Value"
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(),
        "perf105-filter", "DIA"), WorkflowState$new(audit_enabled = FALSE))
    on.exit(closeDIAStateAdapter(manager), add = TRUE)
    manager$.context$reference_states <- TRUE
    manager$saveState("protein_peptide_filtered", object, list(), "source")
    parameters <- list(peptides_per_sample_cutoff = 2L,
        inclusion_list = as.character(object@design_matrix[[object@sample_id]][[1L]]))
    native <- filterMinNumPeptidesPerSample(object,
        peptides_per_sample_cutoff = parameters$peptides_per_sample_cutoff,
        inclusion_list = parameters$inclusion_list, core_utilisation = NA)
    before <- getDIAStateReference(manager)
    value <- computeDIAReferenceSampleFilter(manager, before, parameters)
    actual <- materializeDIATable(manager$.context, value$reference)
    expect_identical(actual, native@peptide_data)
    expect_identical(value$summary,
        native@args$filterMinNumPeptidesPerSample$filter_summary)
    expect_identical(value$execution$bulk_query$settings,
        list(threads = 8L, engine_limit = "12GiB", spill = FALSE))
    expect_identical(value$execution$bulk_query$global_statistics,
        "single_complete_population_pass")
    expect_identical(value$execution$bulk_query$independence_unit,
        "sample_decision_no_batching")
    expect_identical(getDIAStateReference(manager), before)
    expect_identical(manager$.context$store$open_resources, 0L)
})
