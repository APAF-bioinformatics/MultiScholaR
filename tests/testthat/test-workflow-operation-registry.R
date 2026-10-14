# KIT-007: the operation registry behind the execution layer.

test_that("the registry carries the built-in operations with every fact the execution layer reads", {
    names <- workflowOperationNames()
    expect_identical(names, c("DA-contrast-output", "da_publication", "enrichment", "minimum_peptides_per_sample",
        "normalization-ruv", "pearson-pairs", "pearson-selection", "peptide_imputation", "peptide_intensity_filter",
        "peptide_protein_evidence_filter", "peptide_qvalue_filter", "peptide_replicate_filter", "precursor_rollup",
        "protein_accession_cleanup", "protein_duplicate_removal", "protein_intensity_filter", "protein_iq_rollup",
        "protein_replicate_filter", "wfx_fixture"))
    # KIT-010: the protein replicate filter as an operation on the row unit with the protein completion and publication.
    replicate <- workflowActionSpec("protein_replicate_filter")
    expect_identical(replicate$input, "checkpoint")
    expect_identical(replicate$input_class, "ProteinQuantitativeData")
    expect_identical(replicate$compute, "protein")
    expect_identical(replicate$completion, "protein")
    expect_identical(replicate$publication, "protein")
    expect_identical(replicate$variants, c("whole", "rows"))
    expect_identical(replicate$default_variant, "rows")
    expect_identical(replicate$declaration$unit, "row")
    expect_error(validateWorkflowProteinReplicateParameters(list(parallel_cores = 2L), list()), "required")
    expect_error(validateWorkflowProteinReplicateParameters(list(grouping_variable = "group", parallel_cores = 0), list()), "cores")
    expect_null(validateWorkflowProteinReplicateParameters(list(grouping_variable = "group"), list()))
    # KIT-009: the imputation step as an operation on the feature-group unit.
    imputation <- workflowActionSpec("peptide_imputation")
    expect_identical(imputation$step_label, "8_imputed")
    expect_identical(imputation$variants, c("grouped", "partitioned"))
    expect_identical(imputation$default_variant, "partitioned")
    expect_identical(imputation$declaration$unit, "feature_group")
    expect_identical(imputation$declaration$parameters, c("proportion_missing_values", "exclusion_column"))
    expect_identical(workflowQualificationVariant(list(mode = "gui"), "peptide_imputation"), "partitioned")
    expect_error(validateWorkflowPeptideImputationParameters(list(proportion_missing_values = 1.5), list()), "proportion")
    expect_error(validateWorkflowPeptideImputationParameters(list(exclusion_column = "x"), list()), "required")
    expect_null(validateWorkflowPeptideImputationParameters(list(proportion_missing_values = 0.5, exclusion_column = NULL), list()))
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
    # KIT-008: the production variant an unqualified launch runs, in the worker and in process.
    expect_identical(filter$default_variant, "physical")
    expect_identical(rollup$default_variant, "p14")
    expect_null(fixture$default_variant)
    expect_identical(workflowQualificationVariant(list(mode = "gui"), "minimum_peptides_per_sample"), "physical")
    expect_identical(workflowQualificationVariant(list(mode = "gui"), "precursor_rollup"), "p14")
    expect_identical(filter$resolve, "resolveWorkflowSampleFilterParameters")
    expect_true(filter$inclusion_absence)
    for (name in names) if (!identical(workflowActionSpec(name)$compute, "in_process")) expect_true(exists(workflowActionSpec(name)$handler, mode = "function"), info = name)
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
    broken$completion <- "lipid"
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
    expect_identical(sort(workflowOperationNames()), sort(c("DA-contrast-output", "da_publication", "enrichment",
        "minimum_peptides_per_sample", "normalization-ruv", "pearson-pairs", "pearson-selection", "peptide_imputation",
        "peptide_intensity_filter", "peptide_protein_evidence_filter", "peptide_qvalue_filter", "peptide_replicate_filter",
        "precursor_rollup", "protein_accession_cleanup", "protein_duplicate_removal", "protein_intensity_filter",
        "protein_iq_rollup", "protein_replicate_filter", "wfx_fixture")))
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
                "completion", "publication", "step_label", "variants", "default_variant", "diagnostic_only")]
        })
    }, args = list(root = root, library_path = library_path))
    local <- lapply(workflowOperationNames(), function(name) workflowActionSpec(name)[c("input", "input_class", "compute",
        "handler", "validate", "resolve", "parameter_limit", "completion", "publication", "step_label", "variants",
        "default_variant", "diagnostic_only")])
    expect_identical(remote, local)
})

test_that("a worker request's publication context supplies the experiment paths the protein compute records (KIT-010 audit PE-03)", {
    request_dir <- withr::local_tempdir()
    expect_null(workflowRequestExperimentPaths(list(directory = request_dir)))
    expect_null(workflowRequestExperimentPaths(list()))
    saveRDS(list(experiment_label = "x", omic_type = "proteomics", experiment_paths = list(protein_qc_dir = "/p/protein_qc")),
        file.path(request_dir, "publication-context.rds"))
    expect_identical(workflowRequestExperimentPaths(list(directory = request_dir))$protein_qc_dir, "/p/protein_qc")
    saveRDS(list(proteomics = list(protein_qc_dir = "/q/protein_qc")), file.path(request_dir, "publication-context.rds"))
    expect_identical(workflowRequestExperimentPaths(list(directory = request_dir))$protein_qc_dir, "/q/protein_qc")
})

test_that("the DA publication is a registered worker operation with its own completion and publication roles (KIT-019)", {
    entry <- workflowActionSpec("da_publication")
    expect_identical(entry$input, "checkpoint")
    expect_identical(entry$input_class, "ProteinQuantitativeData")
    expect_identical(entry$compute, "protein")
    expect_identical(entry$completion, "da")
    expect_identical(entry$publication, "da")
    expect_identical(entry$parameters, "annotations")
    expect_null(entry$declaration)
    for (name in c(entry$handler, entry$validate)) expect_true(is.function(get(name, mode = "function")))
    expect_error(validateWorkflowDaPublicationParameters(list(other = 1), NULL), "annotations")
    expect_error(validateWorkflowDaPublicationParameters(list(annotations = "a-path"), NULL), "saved table reference")
    expect_null(validateWorkflowDaPublicationParameters(list(), NULL))
})

test_that("the whole-matrix declarations are registered as in-process operations that refuse worker dispatch (KIT-024)", {
    for (name in c("normalization-ruv", "pearson-pairs", "pearson-selection", "DA-contrast-output", "enrichment")) {
        entry <- workflowActionSpec(name)
        expect_identical(entry$compute, "in_process")
        expect_null(entry$handler)
        expect_identical(entry$completion, "in_process")
        expect_identical(entry$declaration$name, name)
        expect_error(assertWorkflowOperationDispatchable(name), class = "workflow_execution_unsupported")
    }
    expect_identical(assertWorkflowOperationDispatchable("protein_intensity_filter")$compute, "protein")
    # An in-process entry must carry its declaration and no handler; a worker entry must name a handler.
    expect_error(validateWorkflowOperationEntry(list(name = "x", input = "checkpoint", input_class = "ProteinQuantitativeData",
        compute = "in_process", handler = "h", validate = "validateWorkflowInProcessParameters", resolve = NULL,
        parameter_limit = 1, parameters = character(), completion = "in_process", publication = "in_process",
        omic_type = "proteomics", variants = "whole", diagnostic_only = FALSE, inclusion_absence = FALSE)),
        class = "workflow_operation_invalid")
})
