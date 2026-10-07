test_that("synthetic cohorts are reproducible, bounded and preserve native sample counts", {
    code <- new.env(parent = asNamespace("MultiScholaR"))
    root <- normalizePath(test_path("..", ".."))
    sys.source(file.path(root, "tools/profiling/workflow_perf107_synthetic.R"), code)
    sys.source(file.path(root, "tests/testthat/helper-workflow-checkpoint-fixtures.R"), code)
    recipe <- code$perf107SyntheticRecipe("calibration", small = TRUE)
    catalogue <- code$perf107SyntheticCatalogue(recipe)
    parts <- lapply(seq_len(recipe$samples), function(i) code$perf107SyntheticSample(recipe, catalogue, i))
    expect_identical(parts[[1L]], code$perf107SyntheticSample(recipe, catalogue, 1L))
    report <- do.call(rbind, parts)
    expect_false(anyDuplicated(report[c("Run", "Precursor.Id")]) > 0L)
    expect_true(anyNA(report$Precursor.Quantity))
    expect_true(any(report$Precursor.Quantity == 0, na.rm = TRUE))
    design <- data.frame(Run = unique(report$Run), group = "control", tech_rep_group = "A")
    object <- MultiScholaR::PeptideQuantitativeDataDiann(peptide_data = report,
        design_matrix = design, technical_replicate_id = "tech_rep_group",
        args = list(globalParameters = list(workflow_type = "DIA")))
    filtered <- MultiScholaR::filterMinNumPeptidesPerSample(object,
        peptides_per_sample_cutoff = recipe$cutoff, core_utilisation = NA)
    expect_identical(unique(as.character(filtered@peptide_data$Run)),
        design$Run[recipe$counts >= recipe$cutoff])
    directory <- withr::local_tempdir()
    result <- code$writePerf107SyntheticCohort(file.path(directory, "complete"), recipe)
    expect_true(result$verified_bounded_scan)
    expect_identical(result$rows, sum(vapply(parts, nrow, integer(1))))
    expect_error(code$writePerf107SyntheticCohort(file.path(directory, "complete"), recipe), "already exists")
    expect_error(code$writePerf107SyntheticCohort(file.path(directory, "failed"), recipe,
        fail_after_sample = 1L), "injected generation failure")
    status <- jsonlite::read_json(file.path(directory, "failed", "generation-status.json"))
    expect_identical(status$status, "failed")
    expect_false(file.exists(file.path(directory, "failed", "manifest.json")))
    other <- code$perf107SyntheticCatalogue(code$perf107SyntheticRecipe("heldout", small = TRUE))
    expect_length(intersect(catalogue$Stripped.Sequence, other$Stripped.Sequence), 0L)
})
