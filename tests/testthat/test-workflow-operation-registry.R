# KIT-007: the operation registry behind the execution layer.

test_that("the registry carries the three built-in operations with every fact the execution layer reads", {
    names <- workflowOperationNames()
    expect_identical(names, c("minimum_peptides_per_sample", "precursor_rollup", "wfx_fixture"))
    fixture <- workflowActionSpec("wfx_fixture")
    expect_true(fixture$diagnostic_only)
    expect_identical(fixture$compute, "fixture")
    expect_identical(fixture$handler, "computeWorkflowFixture")
    rollup <- workflowActionSpec("precursor_rollup")
    expect_identical(rollup$step_label, "4_precursor_rollup")
    expect_identical(rollup$variants, "p14")
    expect_identical(rollup$declaration$unit, "complete_sample")
    expect_identical(rollup$parameter_limit, 65536)
    filter <- workflowActionSpec("minimum_peptides_per_sample")
    expect_identical(filter$step_label, "6_sample_filtered")
    expect_identical(filter$parameter_limit, 262144)
    expect_identical(filter$gui_variants, c("native", "physical"))
    expect_identical(filter$resolve, "resolveWorkflowSampleFilterParameters")
    expect_true(filter$inclusion_absence)
    expect_identical(workflowGuiQualifiedOperation(), "minimum_peptides_per_sample")
    for (name in names) expect_true(exists(workflowActionSpec(name)$handler, mode = "function"))
    for (name in names) expect_true(exists(workflowActionSpec(name)$validate, mode = "function"))
})

test_that("unregistered or malformed operations refuse with the execution layer's class; entries are validated", {
    expect_error(workflowActionSpec("arbitrary_R_function"), class = "workflow_execution_unsupported")
    expect_error(workflowActionSpec(c("a", "b")), class = "workflow_execution_unsupported")
    expect_error(workflowActionSpec(NA_character_), class = "workflow_execution_unsupported")
    expect_error(validateWorkflowActionParameters("unregistered", list(), list()), class = "workflow_execution_unsupported")
    expect_error(workflowQualificationVariant(list(mode = "adapters", qualification = list(programme = "perf106-repair-v2",
        variant = "physical", member_sha256 = strrep("a", 64), deadline = 600)), "precursor_rollup"),
        class = "workflow_execution_invalid")
    expect_identical(workflowQualificationVariant(list(mode = "adapters", qualification = list(programme = "perf106-repair-v2",
        variant = "p14", member_sha256 = strrep("a", 64), deadline = 600)), "precursor_rollup"), "p14")
    entry <- workflowActionSpec("precursor_rollup")
    expect_error(registerWorkflowOperation(entry), class = "workflow_operation_invalid")
    broken <- entry
    broken$completion <- "protein"
    expect_error(registerWorkflowOperation(broken, replace = TRUE), class = "workflow_operation_invalid")
    mismatched <- entry
    mismatched$name <- "renamed"
    expect_error(registerWorkflowOperation(mismatched), class = "workflow_operation_invalid")
    # A registered input class other than the record's refuses before any worker is dispatched.
    record <- list(class = "ProteinQuantitativeData", components = list(peptide_data = list(), design_matrix = list()),
        slots = list(sample_id = "Run", protein_id_column = "Protein", peptide_sequence_column = "Peptide"),
        dimensions = list(columns = c("Run", "Protein", "Peptide")))
    testthat::with_mocked_bindings(
        expect_error(validateWorkflowActionInput("precursor_rollup", NULL, list(state = "owned")), class = "workflow_execution_unsupported"),
        isSavedWorkflowReference = function(x) identical(x, "owned"), isDIAComponentCheckpoint = function(x) TRUE,
        readDIAStateMetadata = function(...) list(data = record))
    expect_identical(sort(workflowOperationNames()), sort(c("minimum_peptides_per_sample", "precursor_rollup", "wfx_fixture")))
})

test_that("a fresh R process resolves the same registry from the operation name alone", {
    skip_if_not_installed("callr")
    library_path <- .libPaths()
    root <- normalizePath(testthat::test_path("..", ".."))
    remote <- callr::r(function(root, library_path) {
        .libPaths(library_path)
        if (file.exists(file.path(root, "DESCRIPTION")) && requireNamespace("pkgload", quietly = TRUE))
            pkgload::load_all(root, quiet = TRUE) else library(MultiScholaR)
        ns <- asNamespace("MultiScholaR")
        lapply(get("workflowOperationNames", ns)(), function(name) {
            entry <- get("workflowActionSpec", ns)(name)
            entry[c("input", "input_class", "compute", "handler", "validate", "resolve", "parameter_limit",
                "completion", "publication", "step_label", "variants", "gui_variants", "diagnostic_only")]
        })
    }, args = list(root = root, library_path = library_path))
    local <- lapply(workflowOperationNames(), function(name) workflowActionSpec(name)[c("input", "input_class", "compute",
        "handler", "validate", "resolve", "parameter_limit", "completion", "publication", "step_label", "variants",
        "gui_variants", "diagnostic_only")])
    expect_identical(remote, local)
})
