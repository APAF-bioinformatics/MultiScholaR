#' Compare one explicit candidate with the unchanged native sample-filter method.
#' @param object Native fixture with declared design and column roles.
#' @param parameters Native cutoff and inclusion arguments.
#' @param mode Physical or selection candidate.
#' @return Invisible NULL after parity and ownership assertions.
#' @noRd
checkPerf106Candidate <- function(object, parameters, mode) {
    native <- filterMinNumPeptidesPerSample(object,
        peptides_per_sample_cutoff = parameters$peptides_per_sample_cutoff,
        inclusion_list = parameters$inclusion_list, core_utilisation = NA)
    root <- withr::local_tempdir()
    store <- newWorkflowEntityStore(root, basename(root), "DIA")
    manager <- newDIAStateAdapter(store, WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$reference_states <- TRUE
    manager$saveState("protein_peptide_filtered", object, list(), "input")
    before <- getDIAStateReference(manager)
    input <- readDIAStateMetadata(manager, before)$data$components$peptide_data
    input_path <- file.path(store$directory, input$payload$relative_path)
    original <- digest::digest(file = input_path, algo = "sha256")
    value <- computeDIAReferenceSampleFilter(manager, before, parameters, output_mode = mode)
    expect_identical(materializeDIATable(manager$.context, value$reference), native@peptide_data)
    expected <- native@args$filterMinNumPeptidesPerSample
    expect_identical(value$summary, expected$filter_summary)
    for (field in c("support_table", "removal_ledger")) expect_identical(
        materializeDIATable(manager$.context,
            value$checkpoint$support$filterMinNumPeptidesPerSample[[field]]), expected[[field]])
    expect_identical(value$counts$rows, as.double(nrow(native@peptide_data)))
    expect_equal(value$counts$summary$total_peptides,
        nrow(unique(native@peptide_data[c("Protein.Ids", "Stripped.Sequence")])))
    expect_equal(value$counts$proteins,
        length(unique(native@peptide_data[[object@protein_id_column]])))
    expect_identical(value$execution$writer_backend, "fast_r")
    expect_identical(value$execution$bulk_query$settings$threads, 8L)
    expect_identical(getDIAStateReference(manager), before)
    expect_identical(digest::digest(file = input_path, algo = "sha256"), original)
    if (mode == "selection") {
        component <- readWorkflowTableComponent(store, value$reference)
        expect_identical(component$source, input)
        expect_false(is.null(component$view))
        scope <- newWorkflowTableReadScope(store, value$reference, reader = "arrow")
        exported <- saveWorkflowTablePhysicalExport(scope, basename(tempfile("export-")),
            newSavedWorkflowScientificIdentity(value$reference$entity_id, "design",
                "export", "all", "perf106"))
        releaseWorkflowTableReadScope(scope)
        expect_identical(materializeDIATable(manager$.context, exported), native@peptide_data)
    }
    expect_identical(store$open_resources, 0L)
    invisible(NULL)
}

test_that("both sample-filter candidates preserve native decisions and scientific outputs", {
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    table <- object@peptide_data
    samples <- unique(table$Run)
    pieces <- lapply(seq_along(samples), function(i) {
        part <- table[table$Run == samples[[i]], , drop = FALSE]
        keys <- paste(part$Protein.Group, part$Stripped.Sequence)
        part[keys %in% head(unique(keys), c(2L, 3L, 4L)[(i - 1L) %% 3L + 1L]), , drop = FALSE]
    })
    object@peptide_data <- do.call(rbind, pieces)
    rownames(object@peptide_data) <- NULL
    object@peptide_data$Precursor.Quantity[] <- NA_real_
    object@peptide_data$Precursor.Normalised[] <- 0
    variants <- list(mixed = object, factor = object, empty = object, boundary = object)
    variants$factor@peptide_data$Run <- factor(object@peptide_data$Run,
        levels = c(rev(samples), "unused"), ordered = TRUE)
    variants$empty@peptide_data <- object@peptide_data[0, , drop = FALSE]
    repeated <- object@peptide_data[rep(seq_len(nrow(object@peptide_data)),
        length.out = 12001L), , drop = FALSE]
    rownames(repeated) <- NULL
    variants$boundary@peptide_data <- repeated
    parameters <- list(
        mixed = list(peptides_per_sample_cutoff = 3L,
            inclusion_list = c(samples[[1L]], samples[[1L]], "absent", " ", NA_character_)),
        all = list(peptides_per_sample_cutoff = 1L, inclusion_list = NULL),
        none = list(peptides_per_sample_cutoff = 100L, inclusion_list = character()))
    for (mode in c("physical", "selection")) {
        for (setting in parameters) checkPerf106Candidate(object, setting, mode)
        for (variant in variants[-1L]) checkPerf106Candidate(variant, parameters$mixed, mode)
    }
    native <- filterMinNumPeptidesPerSample(object, peptides_per_sample_cutoff = 3L,
        inclusion_list = samples[[1L]], core_utilisation = NA)
    without_override <- filterMinNumPeptidesPerSample(object,
        peptides_per_sample_cutoff = 3L, inclusion_list = character(), core_utilisation = NA)
    expect_false(identical(native@peptide_data, without_override@peptide_data))
    # Counting observed quantities would wrongly reject all these real identities.
    expect_gt(nrow(native@peptide_data), 0L)
    expect_identical(sum(!is.na(object@peptide_data$Precursor.Quantity)), 0L)
})

test_that("candidate failures preserve the input head and release owned resources", {
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    root <- withr::local_tempdir()
    store <- newWorkflowEntityStore(root, "perf106-failure", "DIA")
    manager <- newDIAStateAdapter(store, WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$reference_states <- TRUE
    manager$saveState("protein_peptide_filtered", object, list(), "input")
    before <- getDIAStateReference(manager)
    parameters <- list(peptides_per_sample_cutoff = 1L, inclusion_list = NULL)
    for (mode in c("physical", "selection")) {
        original <- appendWorkflowTableWriter
        expect_error(testthat::with_mocked_bindings(
            computeDIAReferenceSampleFilter(manager, before, parameters, output_mode = mode),
            appendWorkflowTableWriter = function(writer, data, ...) {
                if (nrow(data)) stop("injected candidate writer failure")
                original(writer, data, ...)
            }), "injected candidate writer failure")
        expect_identical(getDIAStateReference(manager), before)
        expect_identical(store$open_resources, 0L)
    }
    invalid <- object
    invalid@design_matrix <- rbind(invalid@design_matrix, invalid@design_matrix[1L, ])
    expect_error(filterMinNumPeptidesPerSample(invalid, peptides_per_sample_cutoff = 1L,
        core_utilisation = NA))
    manager$saveState("invalid-design", invalid, list(), "invalid fixture")
    for (mode in c("physical", "selection")) expect_error(computeDIAReferenceSampleFilter(
        manager, getDIAStateReference(manager), parameters, output_mode = mode))
    expect_identical(store$open_resources, 0L)
})

test_that("candidate completion reuses producer counts and preserves native design cleanup", {
    object <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    sample <- unique(object@peptide_data$Run)[[1L]]
    object@peptide_data <- object@peptide_data[object@peptide_data$Run == sample, , drop = FALSE]
    rownames(object@peptide_data) <- NULL
    parameters <- list(peptides_per_sample_cutoff = 1L, inclusion_list = "absent")
    native <- filterMinNumPeptidesPerSample(object, peptides_per_sample_cutoff = 1L,
        inclusion_list = "absent", core_utilisation = NA)
    for (mode in c("physical", "selection")) {
        root <- withr::local_tempdir()
        manager <- newDIAStateAdapter(newWorkflowEntityStore(root, basename(root), "DIA"),
            WorkflowState$new(audit_enabled = FALSE))
        withr::defer(closeDIAStateAdapter(manager))
        manager$.context$reference_states <- TRUE
        manager$saveState("protein_peptide_filtered", object, list(), "input")
        value <- computeDIAReferenceSampleFilter(manager, getDIAStateReference(manager),
            parameters, output_mode = mode)
        workflow <- new.env(parent = emptyenv())
        workflow$state_manager <- manager
        workflow$qc_params <- list(peptide_qc = list())
        product <- testthat::with_mocked_bindings(
            withDIAStateScope(manager, function() finishDIAReferencePeptideStep(workflow,
                value$checkpoint, value$state_name, value$description, value$config, value$section,
                value$counts, observed_inclusion_absence = value$summary$inclusion_entries_absent_from_data)),
            newDIAFilteringAccumulator = function(...) stop("producer counts were discarded"))
        checkpoint <- readDIAStateMetadata(manager, product$stateReference)$data
        expect_identical(materializeDIATable(manager$.context, checkpoint$components$design_matrix),
            native@design_matrix)
        expect_identical(materializeDIATable(manager$.context, checkpoint$components$peptide_data),
            native@peptide_data)
        expect_identical(checkpoint$slots$args$filterMinNumPeptidesPerSample$filter_summary,
            native@args$filterMinNumPeptidesPerSample$filter_summary)
        expect_equal(product$filteringSummary$proteins, countUniqueProteins(native))
        expect_equal(product$filteringSummary$total_peptides, calcTotalPeptides(native))
        # The existing saved-summary contract uses data.frames; native counters
        # return tibbles. Compare every value/column after that container conversion.
        expect_equal(product$filteringSummary$peptides_per_protein,
            as.data.frame(calcPeptidesPerProtein(native)), tolerance = 0)
        expect_equal(product$filteringSummary$peptides_per_run,
            as.data.frame(countPeptidesPerRun(native)), tolerance = 0)
        expect_equal(product$filteringSummary$proteins_per_run,
            as.data.frame(countProteinsPerRun(native)), tolerance = 0)
        expect_identical(manager$.context$store$open_resources, 0L)
        closeDIAStateAdapter(manager)
    }
})
