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
        actual <- computeSavedWorkflowPrecursorRollup(store, input, non_exact = TRUE)
        expect_identical(actual$execution$aggregation_frontend, "duckdb_sum_non_exact")
        actual$execution <- NULL
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
        expect_error(computeSavedWorkflowPrecursorRollup(store, input, non_exact = TRUE), patterns[[name]])
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

coreRollupIdentity <- function(name) list(input = name, design = "paired",
    parameters = "native-rollup", selectors = "all", software = "test")

coreRollupKernel <- function(table) rollUpPrecursorToPeptideHelper(table,
    sample_id_column = "Run", protein_id_column = "Protein.Group",
    peptide_sequence_column = "Stripped.Sequence",
    modified_peptide_sequence_column = "Modified.Sequence",
    precursor_quantity_column = "Precursor.Quantity",
    precursor_normalised_column = "Precursor.Normalised",
    core_utilisation = NA, precursor_id_column = "Precursor.Id",
    return_rollup_result = TRUE)

test_that("the DuckDB sum() rollup arm refuses unless explicitly marked non-exact", {
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    store <- newWorkflowEntityStore(withr::local_tempdir(), "non-exact", "DIA", max_cache_bytes = 0)
    withr::defer(closeWorkflowEntityStore(store))
    input <- saveWorkflowPrecursorRollupInput(store, object, "input", coreRollupIdentity("non-exact"))
    expect_error(computeSavedWorkflowPrecursorRollup(store, input), class = "dia_rollup_non_exact")
    expect_identical(store$open_resources, 0L)
    result <- computeSavedWorkflowPrecursorRollup(store, input, non_exact = TRUE)
    expect_identical(result$execution$aggregation_frontend, "duckdb_sum_non_exact")
    expect_identical(store$open_resources, 0L)
    # Production routes never select the non-exact arm.
    expect_false("non_exact" %in% all.names(body(runDIAReferencePrecursorRollup)))
    expect_false("non_exact" %in% all.names(body(computeDIAReferencePrecursorRollup)))
    expect_false("non_exact" %in% all.names(body(computeDIARelationalRollup)))
})

test_that("indexed rollup keeps kernel order when sample levels are not alphabetical", {
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    table <- object@peptide_data
    runs <- sort(unique(as.character(table$Run)))
    for (ordered in c(FALSE, TRUE)) {
        current <- object
        current@peptide_data$Run <- factor(as.character(table$Run), levels = rev(runs), ordered = ordered)
        expected <- coreRollupKernel(current@peptide_data)
        store <- newWorkflowEntityStore(withr::local_tempdir(), "indexed-order", "DIA", max_cache_bytes = 0)
        identity <- coreRollupIdentity(paste0("indexed-order-", ordered))
        input <- saveWorkflowPrecursorRollupInput(store, current, "input", identity)
        result <- computeSavedWorkflowPrecursorRollup(store, input,
            list(prototype = expected$data[0, , drop = FALSE], identity = identity))
        actual <- readSavedWorkflowTable(store, result$reference, names(expected$data),
            limit = max(1L, nrow(expected$data)))$data
        expect_identical(as.character(actual$Run), as.character(expected$data$Run), info = paste("ordered", ordered))
        expect_identical(actual, expected$data, info = paste("ordered", ordered))
        closeWorkflowEntityStore(store)
    }
})

test_that("blank keys follow the native kernel: factor labels pass, character blanks are refused", {
    # The kernel flags blanks only in character columns (.rollupAssertPopulatedKey);
    # both SQL validators must reproduce that rule, not tighten it.
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    labels <- as.character(object@peptide_data$Protein.Group)
    labels[[1L]] <- ""
    cases <- list(
        factor = factor(labels, levels = c("", setdiff(unique(labels), ""))),
        character = labels)
    for (kind in names(cases)) {
        current <- object
        current@peptide_data$Protein.Group <- cases[[kind]]
        store <- newWorkflowEntityStore(withr::local_tempdir(), paste0("blank-", kind), "DIA",
            max_cache_bytes = 0)
        input <- saveWorkflowPrecursorRollupInput(store, current, "input",
            coreRollupIdentity(paste0("blank-", kind)))
        validate <- function() {
            scope <- newWorkflowTableReadScope(store, input$reference)
            on.exit(releaseWorkflowTableReadScope(scope))
            validateDIARelationalInput(scope, input, NULL)
        }
        if (kind == "factor") {
            expected <- coreRollupKernel(current@peptide_data)
            actual <- computeSavedWorkflowPrecursorRollup(store, input, non_exact = TRUE)
            actual$execution <- NULL
            expect_identical(actual, expected)
            expect_identical(validate(), nrow(current@peptide_data))
        } else {
            expect_error(coreRollupKernel(current@peptide_data), "missing/blank rollup keys")
            expect_error(computeSavedWorkflowPrecursorRollup(store, input, non_exact = TRUE),
                "missing or blank")
            expect_error(validate(), "missing or blank identity keys")
        }
        expect_identical(store$open_resources, 0L, info = kind)
        closeWorkflowEntityStore(store)
    }
})
