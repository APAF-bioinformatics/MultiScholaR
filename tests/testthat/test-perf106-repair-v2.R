# Reuse the frozen native-parity assertion helper without replaying its tests.
perf106_parity_helper <- parse(testthat::test_path("test-perf-106.R"))[[1L]]
stopifnot(identical(perf106_parity_helper[[2L]], as.name("checkPerf106Candidate")))
eval(perf106_parity_helper, envir = environment())

test_that("streamed column identities preserve digest bytes without a serialized buffer", {
    store <- newWorkflowEntityStore(withr::local_tempdir(), "perf106-hash", "DIA")
    withr::defer(closeWorkflowEntityStore(store))
    values <- list(integer(), c(NA_real_, NaN, -0, Inf, 1e16),
        c("", NA_character_, "\u03b2", "\u4f60\u597d"), factor(c("a", NA), levels = c("unused", "a")),
        stats::setNames(1:4, c("one", "two", "three", "four")))
    for (value in values) expect_identical(hashWorkflowSerializedValue(value, store),
        digest::digest(value, algo = "sha256", serializeVersion = 2L))
    value <- list(value = rep(strrep("x", 4096L), 10000L), names = NULL)
    expected <- digest::digest(value, algo = "sha256", serializeVersion = 2L)
    if (isTRUE(capabilities("profmem"))) {
        profile <- withr::local_tempfile()
        utils::Rprofmem(profile, threshold = 1048576)
        actual <- tryCatch(hashWorkflowSerializedValue(value, store),
            finally = utils::Rprofmem(NULL))
        bytes <- suppressWarnings(as.numeric(sub(" .*", "", readLines(profile))))
        expect_lt(max(c(0, bytes), na.rm = TRUE), 4 * 1024^2)
        expect_identical(actual, expected)
    } else expect_identical(hashWorkflowSerializedValue(value, store), expected)
    expect_length(list.files(file.path(store$directory, ".pending")), 0L)
    expect_identical(store$open_resources, 0L)
})

test_that("repaired physical selection and native modes agree on wide shuffled inputs", {
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    original <- object@peptide_data
    samples <- unique(original$Run)
    keys <- paste(original$Protein.Group, original$Stripped.Sequence)
    # Remove identities from one sample; zero/missing quantity never changes membership.
    table <- original[original$Run != samples[[1L]] | keys == keys[[1L]], , drop = FALSE]
    table <- table[rev(seq_len(nrow(table))), , drop = FALSE]
    rownames(table) <- NULL
    table$Protein.Ids <- paste0(table$Protein.Ids, strrep("x", 4096L))
    table$Precursor.Quantity[] <- NA_real_
    table$Precursor.Normalised[] <- 0
    object@peptide_data <- table
    for (mode in c("native", "physical", "selection")) {
        checkPerf106Candidate(object,
            list(peptides_per_sample_cutoff = 3L, inclusion_list = c(samples[[1L]], "absent")), mode)
    }
})

test_that("wide admission bounds transfers and validates globally before three count batches", {
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    object@peptide_data$Protein.Ids <- paste0(object@peptide_data$Protein.Ids, strrep("x", 4096L))
    store <- newWorkflowEntityStore(withr::local_tempdir(), "perf106-transfer", "DIA")
    manager <- newDIAStateAdapter(store, WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$reference_states <- TRUE
    manager$saveState("protein_peptide_filtered", object, list(), "fixture")
    reference <- getDIAStateReference(manager)
    value <- computeDIAReferenceSampleFilter(manager, reference,
        list(peptides_per_sample_cutoff = 3L), output_mode = "physical")
    expect_lt(value$execution$admission$batch_rows, 10000L)
    expect_lte(value$execution$admission$batch_rows * value$execution$admission$max_row_bytes,
        33554432)
    expect_length(value$execution$admission$count_batch_sizes, 3L)
    expect_equal(sum(value$execution$admission$count_batch_sizes), length(unique(object@peptide_data$Run)))
    expect_identical(value$execution$bulk_query$independence_unit,
        "complete_sample_after_global_validation")
    expect_identical(store$open_resources, 0L)
    expect_identical(getDIAStateReference(manager), reference)
    expect_error(computeDIAReferenceSampleFilter(manager, reference,
        list(peptides_per_sample_cutoff = 3L), output_mode = "physical", worker_cap = 1),
        class = "workflow_execution_unsupported")
    expect_identical(store$open_resources, 0L)
})

test_that("worker count references register in an already open controller store", {
    root <- withr::local_tempdir()
    controller <- newWorkflowEntityStore(root, "perf106-count-handoff", "DIA")
    withr::defer(closeWorkflowEntityStore(controller))
    worker <- newWorkflowEntityStore(root, "perf106-count-handoff", "DIA")
    withr::defer(closeWorkflowEntityStore(worker))
    data <- data.frame(Run = c("r1", "r2"), Protein.Ids = c("p1", "p2"),
        Stripped.Sequence = c("s1", "s2"), Protein.Group = c("g1", "g2"))
    identity <- newSavedWorkflowScientificIdentity("input", "design", "params", "rows", "software")
    table <- saveWorkflowTable(worker, data, "handoff-table", identity, "Run")
    accumulator <- newDIAFilteringAccumulator(data[0, ], "Protein.Group")
    accumulator$push(data)
    counts <- accumulator$finish()
    counts$padding <- raw(5 * 1024^2)
    reference <- saveDIALargeFilteringCountsReference(worker,
        list(reference = table, counts = counts), "Protein.Group")
    expect_null(describeSavedWorkflowEntity(controller, reference$entity_id))
    scope <- newWorkflowTableReadScope(controller, table, reader = "arrow")
    withr::defer(releaseWorkflowTableReadScope(scope))
    actual <- summarizeSavedDIAFiltering(scope, NULL, "Protein.Group", reference)
    expect_identical(actual$summary, counts$summary)
    expect_identical(actual$dependencies[[1L]], reference)
    expect_identical(describeSavedWorkflowEntity(controller, reference$entity_id), reference)
})

test_that("qualification variants are explicit immutable server choices", {
    # KIT-008: an unqualified launch runs the registry's production variant (KIT-034, quarantine Q-09).
    expect_identical(workflowQualificationVariant(list(mode = "gui"),
        "minimum_peptides_per_sample"), workflowActionSpec("minimum_peptides_per_sample")$default_variant)
    config <- list(mode = "adapters", qualification = list(programme = "perf106-repair-v2",
        variant = "physical", member_sha256 = strrep("a", 64), deadline = 300L))
    for (variant in c("native", "physical", "selection")) {
        config$qualification$variant <- variant
        expect_identical(workflowQualificationVariant(config, "minimum_peptides_per_sample"), variant)
        expect_error(workflowQualificationVariant(config, "precursor_rollup"),
            class = "workflow_execution_invalid")
    }
    config$qualification$variant <- "p14"
    expect_identical(workflowQualificationVariant(config, "precursor_rollup"), "p14")
    config$mode <- "gui"
    expect_error(workflowQualificationVariant(config, "precursor_rollup"),
        class = "workflow_execution_invalid")
})

test_that("repaired ranking uses both widths and no superseded relative speed gate", {
    bindings <- populations <- list(primary = list(), wide = list())
    for (population in names(populations)) for (arm in c("native", "physical", "selection")) {
        binding <- digest::digest(list(population, arm), algo = "sha256")
        bindings[[population]][[arm]] <- binding
        populations[[population]][[arm]] <- lapply(1:2, function(rep) list(
            status = "complete", bindings = binding, rep = rep,
            request_id = paste(population, arm, rep),
            science = TRUE, outputs = TRUE, measurement_valid = TRUE,
            integrity = TRUE, guard = TRUE, admission = TRUE,
            peak_bytes = if (arm == "native") 3e9 else if (arm == "physical") 2e9 else 1.5e9,
            elapsed_seconds = if (arm == "selection") 100 else 20,
            swap_bytes = 0, spill_bytes = 0))
    }
    expect_identical(rankPerfRepairScreen(populations, bindings)$finalist, "selection")
    for (i in 1:2) populations$wide$selection[[i]]$peak_bytes <- 2.5e9
    expect_identical(rankPerfRepairScreen(populations, bindings)$finalist, "physical")
    populations$wide$physical[[1L]]$status <- "memory_guard"
    expect_identical(rankPerfRepairScreen(populations, bindings)$finalist, "selection")
    populations$primary$selection[[1L]]$status <- "time_guard"
    expect_identical(rankPerfRepairScreen(populations, bindings)$status, "no_finalist")
    populations$wide$native[[1L]]$status <- "memory_guard"
    expect_identical(rankPerfRepairScreen(populations, bindings)$status, "held")
})
