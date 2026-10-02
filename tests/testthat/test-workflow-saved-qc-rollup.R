test_that("saved rollup matches native populations, types, missingness and annotations", {
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    table <- object@peptide_data
    key <- paste(table$Run, table$Protein.Group, table$Stripped.Sequence)
    groups <- unique(key)
    table$Precursor.Quantity[key == groups[[1L]]] <- NA_real_
    table$Precursor.Normalised[which(key == groups[[2L]])[[1L]]] <- NA_real_
    table$identification_peptide_count <- 3L
    table$identification_peptidoform_count <- 6L
    table$peptides_for_protein_count <- 3L
    table$peptidoforms_for_protein_count <- 6L
    attr(table$Protein.Group, "identity_role") <- "quantitative protein group"
    object@peptide_data <- table
    for (factors in c(FALSE, TRUE)) {
        current <- object
        if (factors) current@peptide_data$Run <- factor(table$Run,
            levels = c(rev(unique(table$Run)), "unused_run"), ordered = TRUE)
        expected <- rollUpPrecursorToPeptideHelper(current@peptide_data,
            sample_id_column = "Run", protein_id_column = "Protein.Group",
            peptide_sequence_column = "Stripped.Sequence",
            modified_peptide_sequence_column = "Modified.Sequence",
            precursor_quantity_column = "Precursor.Quantity",
            precursor_normalised_column = "Precursor.Normalised",
            core_utilisation = NA, precursor_id_column = "Precursor.Id",
            return_rollup_result = TRUE)
        store <- newWorkflowEntityStore(withr::local_tempdir(), "rollup-test", "DIA", max_cache_bytes = 0)
        input <- saveWorkflowPrecursorRollupInput(store, current, "input",
            list(input = "fixture", design = "paired", parameters = "native-rollup",
                selectors = "all", software = "test"))
        actual <- computeSavedWorkflowPrecursorRollup(store, input)
        expect_identical(actual, expected)
        expect_identical(store$open_resources, 0L)
        expect_gt(actual$summary$raw_all_missing_outputs, 0)
        expect_gt(actual$summary$normalised_partial_missing_outputs, 0)
        expect_identical(levels(actual$data$Run), levels(expected$data$Run))
        closeWorkflowEntityStore(store)
    }
})

test_that("saved rollup refuses invalid science and releases every query lease", {
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    cases <- list(
        missing_key = function(table) { table$Precursor.Id[[1L]] <- NA_character_; table },
        duplicate = function(table) rbind(table, table[1L, ]),
        negative = function(table) { table$Precursor.Quantity[[1L]] <- -1; table },
        infinite = function(table) { table$Precursor.Normalised[[1L]] <- Inf; table },
        nan = function(table) { table$Precursor.Quantity[[1L]] <- NaN; table },
        evidence = function(table) {
            table$identification_peptide_count <- 3L
            table$identification_peptide_count[[1L]] <- 4L
            table
        })
    patterns <- c(missing_key = "missing or blank", duplicate = "duplicate precursor",
        negative = "negative or non-finite", infinite = "negative or non-finite",
        nan = "negative or non-finite", evidence = "conflicting frozen")
    for (name in names(cases)) {
        invalid <- object
        invalid@peptide_data <- cases[[name]](object@peptide_data)
        store <- newWorkflowEntityStore(withr::local_tempdir(), "invalid-rollup", "DIA", max_cache_bytes = 0)
        input <- saveWorkflowPrecursorRollupInput(store, invalid, "input",
            list(input = name, design = "paired", parameters = "native-rollup",
                selectors = "all", software = "test"))
        expect_error(computeSavedWorkflowPrecursorRollup(store, input), patterns[[name]])
        expect_identical(store$open_resources, 0L)
        expect_true(file.exists(file.path(store$directory, input$reference$payload$relative_path)))
        closeWorkflowEntityStore(store)
    }
    store <- newWorkflowEntityStore(withr::local_tempdir(), "logged-rollup", "DIA")
    object@is_logged_data <- TRUE
    expect_error(saveWorkflowPrecursorRollupInput(store, object, "input",
        list()), "linear peptide")
    closeWorkflowEntityStore(store)
})

test_that("rollup collection fills one typed result and removes its temporary query", {
    connection <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE), dbdir = ":memory:")
    withr::defer(DBI::dbDisconnect(connection, shutdown = TRUE))
    sql <- paste("SELECT CAST(i AS INTEGER) AS id, CAST(i / 7 AS DOUBLE) AS quantity,",
        "CASE WHEN i % 3 = 0 THEN NULL ELSE 'peptide' END AS label",
        "FROM range(25005) AS input(i) ORDER BY i")
    expected <- DBI::dbGetQuery(connection, sql)
    expect_identical(collectWorkflowRollupQuery(connection, sql, 25005), expected)
    empty <- "SELECT CAST(NULL AS INTEGER) AS id, CAST(NULL AS DOUBLE) AS quantity WHERE FALSE"
    expect_identical(collectWorkflowRollupQuery(connection, empty, 0), DBI::dbGetQuery(connection, empty))
    expect_error(collectWorkflowRollupQuery(connection, sql, 25004), class = "saved_entity_limit")
    expect_false("qc_collection" %in% DBI::dbListTables(connection))
    expect_identical(DBI::dbGetQuery(connection, "SELECT 1 AS ok")$ok, 1L)
})

test_that("rollup collection does not allocate a complete numeric column per batch", {
    skip_if_not(isTRUE(capabilities("profmem")))
    connection <- DBI::dbConnect(duckdb::duckdb(shared_home = FALSE), dbdir = ":memory:")
    withr::defer(DBI::dbDisconnect(connection, shutdown = TRUE))
    sql <- "SELECT CAST(i AS DOUBLE) AS value FROM range(200000) AS input(i)"
    # Warm dispatch on tiny data; only the following collection is observed.
    collectWorkflowRollupQuery(connection, "SELECT 1.0::DOUBLE AS value", 1)
    profile <- withr::local_tempfile()
    withr::defer(utils::Rprofmem(NULL))
    utils::Rprofmem(profile, threshold = 1000000)
    value <- collectWorkflowRollupQuery(connection, sql, 200000)
    utils::Rprofmem(NULL)
    expect_identical(value$value, as.numeric(0:199999))
    allocations <- readLines(profile, warn = FALSE)
    sizes <- suppressWarnings(as.numeric(sub(" .*", "", allocations)))
    # One destination plus bounded driver overhead, not twenty destination copies.
    expect_lte(sum(sizes >= 200000 * 8, na.rm = TRUE), 6L)
    receipt <- Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT")
    if (nzchar(receipt)) expect_true(file.copy(profile,
        file.path(dirname(receipt), "rollup-collection.Rprofmem"), overwrite = FALSE))
})

test_that("saved rollup uses the matching envelope after reopen without decoding the predecessor", {
    root <- withr::local_tempdir()
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    object@args$rollUpPrecursorToPeptide$core_utilisation <- newDIAWorkerSpec(2L)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "envelope", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    manager$.context$sql_rollup <- TRUE
    manager$saveState("protein_peptide_filtered", object, list(), "source")
    prepareDIACapacityRollupInput(manager, object)
    source <- manager$.context$snapshot$states$protein_peptide_filtered
    preparation <- manager$.context$snapshot$rollup_input
    expect_identical(manager$getState(), object)
    expect_identical(manager$getHistory(), c("initial", "protein_peptide_filtered"))
    closeDIAStateAdapter(manager)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, "envelope", "DIA"), reopen = TRUE)
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$sql_rollup <- TRUE
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- manager
    original_reader <- hydrateSavedWorkflowCheckpoint
    decoded <- character()
    now <- function() as.POSIXct("2026-09-23 00:00:00", tz = "UTC")
    result <- testthat::with_mocked_bindings(
        withDIAStateScope(manager, function() runPeptideRollupApplyStep(workflow, nowFn = now)),
        hydrateSavedWorkflowCheckpoint = function(store, reference) {
            if (identical(reference$entity_id, source$entity_id)) stop("full predecessor decoded")
            decoded <<- c(decoded, reference$entity_id)
            original_reader(store, reference)
        })
    expect_true(preparation$envelope$entity_id %in% decoded)
    expect_false(source$entity_id %in% decoded)
    expect_identical(manager$getHistory(), c("initial", "protein_peptide_filtered", "precursor_rollup"))
    expect_identical(diaStateOwnership(manager)$working_records, 0L)
    expect_identical(manager$.context$store$open_resources, 0L)
    expect_identical(manager$revertToState("protein_peptide_filtered"), object)
    expected <- testthat::with_mocked_bindings(
        withDIAStateScope(manager, function() runPeptideRollupApplyStep(workflow, nowFn = now)),
        getDIARollupEnvelope = function(...) NULL)
    expect_identical(result, expected)
    expect_identical(manager$revertToState("protein_peptide_filtered"), object)
    manager$.context$rollup_input$source <- "different-generation"
    expect_error(getDIARollupEnvelope(manager, source$entity_id), class = "saved_entity_scope")
    manager$.context$rollup_input <- manager$.context$snapshot$rollup_input
    manager$audit_enabled <- TRUE
    expect_null(getDIARollupEnvelope(manager, source$entity_id))
    manager$audit_enabled <- FALSE
    manager$.context$rollup_input$envelope <- NULL
    expect_null(getDIARollupEnvelope(manager, source$entity_id))
    manager$.context$rollup_input <- preparation
    path <- file.path(manager$.context$store$directory, preparation$envelope$payload$relative_path)
    writeBin(as.raw(1:8), path)
    before <- manager$.context$head
    expect_error(withDIAStateScope(manager, function() runPeptideRollupApplyStep(workflow)),
        class = "saved_entity_corrupt")
    expect_identical(manager$.context$head, before)
    expect_identical(manager$current_state, "protein_peptide_filtered")
    expect_identical(manager$.context$store$open_resources, 0L)
})

test_that("optimized R rollup matches the frozen pre-repair kernel on compatibility cases", {
    receipt <- Sys.getenv("MULTISCHOLAR_DIA008_RECEIPT")
    skip_if(!nzchar(receipt), "frozen installed-kernel oracle is exercised under the guarded runner")
    frozen_path <- file.path(dirname(dirname(receipt)), "browser-smoke.SW6Naw", "source", "R",
        "func_prot_rollup_methods.R")
    expect_true(file.exists(frozen_path))
    frozen <- new.env(parent = asNamespace("MultiScholaR"))
    helpers <- c(".rollupSumObservedOrNA", ".rollupAssertPopulatedKey", ".rollupAssertLinearQuantity",
        ".rollupDuplicateExamples", ".rollupSummarisePrepared", "rollUpPrecursorToPeptideHelper")
    for (expr in parse(frozen_path)) {
        if (is.call(expr) && identical(expr[[1L]], as.name("<-")) && is.symbol(expr[[2L]]) &&
            as.character(expr[[2L]]) %in% helpers) eval(expr, frozen)
    }
    expect_true(all(vapply(helpers, exists, logical(1), envir = frozen, inherits = FALSE)))
    table <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))@peptide_data
    table$identification_peptide_count <- 3L
    table$identification_peptidoform_count <- 6L
    invoke <- function(fn, data) fn(data, sample_id_column = "Run", protein_id_column = "Protein.Group",
        peptide_sequence_column = "Stripped.Sequence", modified_peptide_sequence_column = "Modified.Sequence",
        precursor_quantity_column = "Precursor.Quantity", precursor_normalised_column = "Precursor.Normalised",
        precursor_id_column = "Precursor.Id", precursor_charge_column = "Precursor.Charge",
        core_utilisation = NA, return_rollup_result = TRUE)
    variants <- list(original = table, empty = table[0, , drop = FALSE])
    missing <- table
    missing$Precursor.Quantity[1:2] <- NA_real_
    missing$Precursor.Normalised[3L] <- NA_real_
    missing$Precursor.Quantity[4L] <- 0
    missing$Run <- factor(missing$Run, levels = c(rev(unique(missing$Run)), "unused"), ordered = TRUE)
    missing$Protein.Ids[1:3] <- c("P1; P2", " P2;P3;", "P1; P2")
    variants$missing_factor_annotation <- missing
    integer_data <- table
    integer_data$Precursor.Quantity <- as.integer(seq_len(nrow(table)) %% 37L)
    integer_data$Precursor.Normalised <- as.integer(seq_len(nrow(table)) %% 19L)
    variants$integer <- integer_data
    integer_data$Precursor.Quantity[1:2] <- NA_integer_
    variants$integer_missing <- integer_data
    large <- table
    large$Precursor.Quantity[1:2] <- c(1e16, 1)
    extra <- large[1L, , drop = FALSE]
    extra$Precursor.Id <- paste0(extra$Precursor.Id, "_extra")
    extra$Precursor.Quantity <- 1
    variants$rounding <- rbind(large, extra)
    for (name in names(variants)) {
        expected <- invoke(frozen$rollUpPrecursorToPeptideHelper, variants[[name]])
        actual <- invoke(rollUpPrecursorToPeptideHelper, variants[[name]])
        expect_equal(actual$data, expected$data, tolerance = 1e-12, info = name)
        expect_identical(actual$summary, expected$summary, info = name)
        expect_identical(actual$data$Peptide.RawQuantity, expected$data$Peptide.RawQuantity, info = name)
        expect_identical(actual$data$Peptide.Normalised, expected$data$Peptide.Normalised, info = name)
        expect_identical(lapply(actual$data, typeof), lapply(expected$data, typeof), info = name)
        keys <- c("Run", "Protein.Group", "Stripped.Sequence")
        expect_identical(actual$data[keys], expected$data[keys], info = name)
    }
})


test_that("rollup retains native summation at a floating threshold boundary", {
    prepared <- data.frame(.rollup_sample = rep("S1", 3L), .rollup_protein = rep("P1", 3L),
        .rollup_peptide = rep("PEPTIDE", 3L), .rollup_modified = c("A", "B", "C"),
        .rollup_raw = c(1e16, 1, 1), .rollup_normalised = c(1e16, 1, 1))
    value <- .rollupSummarisePrepared(prepared, NA)$data
    expect_identical(value$Peptide.RawQuantity, sum(prepared$.rollup_raw))
    expect_identical(value$Peptide.Normalised, sum(prepared$.rollup_normalised))
    expect_true(value$Peptide.RawQuantity > 1e16)
})
