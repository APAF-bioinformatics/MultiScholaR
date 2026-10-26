test_that("protein QC counts declared samples without treating the key as a run", {
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    before <- serialize(object, NULL)
    samples <- as.character(object@design_matrix[[object@sample_id]])
    expected <- vapply(samples, \(sample) sum(!is.na(object@protein_quant_table[[sample]])),
        integer(1))
    counts <- countProteinsPerRun(object)
    expect_identical(counts$Run, sort(samples))
    expect_equal(counts$n_proteins, unname(expected[sort(samples)]))
    expect_equal(countUniqueProteins(object), nrow(object@protein_quant_table))
    expect_identical(serialize(object, NULL), before)

    legacy <- object
    names(legacy@protein_quant_table)[names(legacy@protein_quant_table) == "Protein.Group"] <- "Protein.Ids"
    names(legacy@protein_id_table)[names(legacy@protein_id_table) == "Protein.Group"] <- "Protein.Ids"
    legacy@protein_id_column <- "Protein.Ids"
    expect_equal(countProteinsPerRun(legacy), counts)
    expect_equal(countProteinsPerRun(legacy@protein_quant_table), counts)
    missing_sample <- object
    missing_sample@protein_quant_table[[samples[[1L]]]] <- NULL
    expect_error(countProteinsPerRun(missing_sample), "Declared sample columns")
})

test_that("rolled peptide validity preserves raw confidence requirements and full history", {
    raw <- makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small"))
    raw_before <- serialize(raw, NULL)
    broken_raw <- raw
    broken_raw@peptide_data$Q.Value <- NULL
    expect_error(methods::validObject(broken_raw), "Q value column must be")
    rolled <- rollUpPrecursorToPeptide(raw, core_utilisation = 1)
    expect_false("Q.Value" %in% names(rolled@peptide_data))
    expect_true(methods::validObject(rolled))
    expect_identical(serialize(raw, NULL), raw_before)
    broken <- rolled
    broken@args$rollUpPrecursorToPeptide <- NULL
    expect_error(methods::validObject(broken), "Q value column must be")
    broken <- rolled
    broken@args$rollUpPrecursorToPeptide$rollup_summary$active_protein_key <- "wrong"
    expect_error(methods::validObject(broken), "Q value column must be")
    broken <- rolled
    broken@peptide_data$Peptide.RawQuantity <- NULL
    expect_error(methods::validObject(broken), "Q value column must be")

    manager <- WorkflowState$new()
    manager$saveState("raw", raw, list(), "raw confidence retained")
    manager$saveState("rolled", rolled, list(), "actual precursor rollup")
    path <- withr::local_tempfile(fileext = ".rds")
    saveRDS(manager, path)
    restored <- readRDS(path)
    expect_identical(restored$getState("raw"), raw)
    expect_identical(restored$getState("rolled"), rolled)
    expect_true(methods::validObject(restored$getState("rolled")))
    expect_identical(restored$revertToState("raw"), raw)
})

test_that("protein tracking accepts the real S4 and preserves it", {
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    before <- serialize(object, NULL)
    rlang::local_bindings(filtering_progress = methods::new("FilteringProgress"),
        .env = .GlobalEnv)
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    grid <- updateProteinFiltering(object, "protein", return_grid = TRUE)
    expect_s3_class(grid, "gtable")
    progress <- get("filtering_progress", envir = .GlobalEnv)
    expect_equal(progress@proteins, nrow(object@protein_quant_table))
    expect_identical(progress@proteins_per_run[[1L]]$Run,
        sort(as.character(object@design_matrix$Run)))
    expect_identical(serialize(object, NULL), before)
})

test_that("exported session counts use the same declared key as the saved object", {
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- WorkflowState$new()
    workflow$state_manager$saveState("correlation_filtered", object, list(), "final")
    norm <- list(correlation_filtered_obj = object, best_k = NA_integer_)
    session <- MultiScholaR:::collectProtNormExportSessionData(workflow, norm,
        list(norm_method = "cyclicloess", ruv_mode = "skip"),
        messageFn = \(...) invisible(NULL))
    expect_equal(session$final_protein_count, nrow(object@protein_quant_table))
    expect_identical(session$current_s4_object, object)
    expect_identical(session$correlation_filtered_s4, object)
    expect_identical(session$r6_complete_states, workflow$state_manager$states)
})

test_that("duplicate removal uses the declared protein group without changing the estimator", {
    skip("quarantined: Q-22 (KIT-048) the step on an in-memory manager is retired; the owner is the one path")
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    object@protein_quant_table <- rbind(object@protein_quant_table,
        object@protein_quant_table[1L, , drop = FALSE])
    original <- serialize(object, NULL)
    manager <- WorkflowState$new()
    manager$saveState("protein", object, list(), "duplicate test")
    result <- runProteinDuplicateRemovalStep(list(state_manager = manager), "median",
        logInfoFn = \(message) invisible(NULL))
    expected <- object@protein_quant_table |>
        dplyr::group_by(.data$Protein.Group) |>
        dplyr::summarise(dplyr::across(dplyr::all_of(object@design_matrix$Run),
            \(values) stats::median(values, na.rm = TRUE)), .groups = "drop")
    expect_identical(result$deduplicatedS4@protein_quant_table, expected)
    expect_identical(result$duplicates, object@protein_quant_table$Protein.Group[[1L]])
    expect_identical(result$deduplicatedS4@design_matrix, object@design_matrix)
    expect_identical(result$deduplicatedS4@protein_id_table, object@protein_id_table)
    expect_true(methods::validObject(result$deduplicatedS4))
    expect_identical(serialize(object, NULL), original)
})

test_that("RUV and correlation tracking retain declared protein identities", {
    skip("quarantined: Q-22 (KIT-048) the step on an in-memory manager is retired; the owner is the one path")
    object <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    before <- serialize(object, NULL)
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- WorkflowState$new()
    workflow$tab_status <- list(normalization = "pending")
    workflow$state_update_trigger <- 0
    norm <- new.env(parent = emptyenv())
    norm$best_k <- 1L
    tracking <- function(data, ...) {
        expect_identical(data, object)
        "tracking_plot"
    }
    cleaned <- MultiScholaR:::finalizeProtNormRuvCleanupStep(
        object, list(ruv_percentage = 10, norm_method = "median", ruv_mode = "manual"),
        norm, workflow, "proteomics", "test",
        removeRowsWithMissingValuesPercentFn = \(theObject) theObject,
        updateProteinFilteringFn = tracking, messageFn = \(...) invisible(NULL))
    expect_identical(cleaned, object)
    expect_equal(workflow$protein_counts$after_ruv_filtering, nrow(object@protein_quant_table))
    expect_identical(norm$post_norm_filtering_plot, "tracking_plot")
    MultiScholaR:::updateProtNormFinalFilteringPlot(object, norm, "proteomics", "test",
        updateProteinFilteringFn = tracking, messageFn = \(...) invisible(NULL))
    for (skipped in c(FALSE, TRUE)) {
        result <- MultiScholaR:::finalizeProtNormCorrelationWorkflowState(
            object, workflow, correlationThreshold = 0.7, skipped = skipped,
            messageFn = \(...) invisible(NULL), catFn = \(...) invisible(NULL))
        expect_equal(result$finalProteinCount, nrow(object@protein_quant_table))
        expect_equal(result$finalSampleCount, nrow(object@design_matrix))
    }
    expect_identical(serialize(object, NULL), before)
})
