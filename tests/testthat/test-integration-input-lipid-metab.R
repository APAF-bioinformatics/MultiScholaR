# KIT-135: lipidomics and metabolomics write the declared integration input (KIT-061's contract,
# utils_workflow_integration_input.R) beside today's S4 RDS: in process at the Summary hook (`saveWorkflowSummaryParameters`
# under the omic's group) and on the worker (`<prefix>_integration_export`). Fixtures: the declarations' owners and
# journeys (helper-workflow-lipid-declarations.R, helper-workflow-metab-declarations.R).

integrationOmicsQuiet <- function(expr) {
    value <- NULL
    utils::capture.output(value <- suppressWarnings(suppressMessages(expr)))
    value
}

integrationOmicsOwner <- function(omic, root, env = parent.frame()) {
    fixture <- if (identical(omic, "lipidomics")) lipidDeclarationsOwner(root) else metabDeclarationsOwner(root)
    withr::defer(closeDIAStateAdapter(fixture$owner), envir = env)
    fixture
}

integrationOmicsSummarySave <- function(omic, fixture, root, label = "kit135") {
    paths <- list(source_dir = file.path(root, "scripts", omic), integration_dir = file.path(root, "integration"), base_dir = root)
    dir.create(paths$source_dir, recursive = TRUE, showWarnings = FALSE)
    shown <- character()
    saved <- NULL
    failure <- NULL
    integrationOmicsQuiet(saveWorkflowSummaryParameters(omicGuiConfig(omic)$summary, fixture$workflow,
        projectDirs = stats::setNames(list(paths), omic), omicType = omic, label = label, description = "d",
        notify = function(message, ...) shown <<- c(shown, message), done = function(value) saved <<- value,
        failed = function(error) failure <<- error, createWorkflowArgsFn = function(...) "study_parameters.txt",
        writeAnalysisRecordFn = function(...) NULL, catFn = function(...) NULL))
    if (!is.null(failure)) stop(failure)
    list(shown = shown, saved = saved, dir = file.path(paths$integration_dir, omic), integration_dir = paths$integration_dir,
        rds = file.path(paths$integration_dir, sprintf("%s_%s_final_s4.RDS", omic, label)))
}

#' The state's assays as the input's matrices: `feature_id`, then the samples of the input in its order.
integrationOmicsExpected <- function(omic, owner, state, samples) {
    binding <- omicsDeclarationBindings()[[omic]]
    object <- readWorkflowOmicsState(binding, owner, state)
    id <- methods::slot(object, binding$id_slot)
    lapply(methods::slot(object, binding$data_slot), function(table) {
        table <- as.data.frame(table, check.names = FALSE)
        data.frame(feature_id = as.character(table[[id]]), table[samples], check.names = FALSE, stringsAsFactors = FALSE)
    })
}

#' The journey up to the between-sample normalisation (the declarations' journeys' first steps).
integrationOmicsNormalised <- function(fixture, itsd_ids = NULL) {
    owner <- fixture$owner
    config <- fixture$workflow$config_list
    integrationOmicsQuiet({
        runLipidReferenceIntensityFilter(owner, 1, 0.5, config, hydrate = FALSE)
        runLipidReferenceDuplicateResolution(owner, config, hydrate = FALSE)
        commitLipidReferenceStoreStep(owner, "lipid_qc_finalize", "lipid_qc_complete", config, "QC processing complete")
        runLipidReferenceItsdNormalisation(owner, "median", itsd_ids, config, hydrate = FALSE)
        runLipidReferenceLog2(owner, 1, config, hydrate = FALSE)
        runLipidReferenceBetweenSampleNormalisation(owner, "cyclicloess", config, hydrate = FALSE)
    })
    invisible(fixture)
}

for (omic in c("lipidomics", "metabolomics")) local({
    omic <- omic
    prefix <- if (identical(omic, "lipidomics")) "lipid" else "metab"
    test_that(sprintf("KIT-135: %s with today's RUV writes the corrected matrix per assay as omic_assay views, beside the RDS", omic), {
        root <- withr::local_tempdir()
        fixture <- integrationOmicsOwner(omic, file.path(root, "owner"))
        journey <- integrationOmicsQuiet(if (identical(omic, "lipidomics")) lipidDeclarationsJourney(fixture) else
            metabDeclarationsJourney(fixture))
        owner <- fixture$owner
        final <- paste0(prefix, "_correlation_filtered")
        expect_identical(owner$current_state, final)
        run <- integrationOmicsSummarySave(omic, fixture, root)
        # Today's file and notice are unchanged; the input is the omic's folder beside it.
        expect_true("Saved Integration S4 Object" %in% run$shown)
        expect_identical(run$saved$integrationFile, sprintf("%s_kit135_final_s4.RDS", omic))
        expect_s4_class(readRDS(run$rds), omicsDeclarationBindings()[[omic]]$class)
        expect_setequal(list.files(run$integration_dir), c(omic, basename(run$rds)))
        input <- readIntegrationInputWithMatrices(run$dir)
        manifest <- input$manifest
        assays <- names(methods::slot(fixture$object, omicsDeclarationBindings()[[omic]]$data_slot))
        views <- paste(omic, assays, sep = "_")
        expect_identical(names(input$views), views)
        expect_identical(unlist(lapply(manifest$source_state$assays, `[[`, "assay")), stats::setNames(assays, views))
        expect_identical(manifest$omic_type, omic)
        # KIT-153: the role of each data set's matrix file (the normalised one); the corrected matrix beside it.
        expect_identical(manifest$matrix_role, "normalised")
        expect_identical(manifest$scale, "log2")
        expect_identical(manifest$normalisation$source$name, paste0(prefix, "_normalized"))
        expect_identical(manifest$normalisation$source$entity_id, getWorkflowStateReference(owner, paste0(prefix, "_normalized"))$entity_id)
        expect_identical(manifest$source_state$name, final)
        expect_identical(manifest$source_state$entity_id, getWorkflowStateReference(owner, final)$entity_id)
        correction <- manifest$correction
        expect_identical(correction[c("state", "id", "method", "ruv_mode", "replicate_structure", "grouping_variable", "matrix",
            "da_tests")], list(state = paste0(prefix, "_ruv_corrected"), id = "C4", method = "RUVIII_C_Varying", ruv_mode = "manual",
            replicate_structure = "group", grouping_variable = "group", matrix = "ruv_corrected_values", da_tests = "this_matrix"))
        expect_identical(unlist(correction$grouping_levels), c("WT", "KO"))
        expect_identical(names(correction$views), views)
        for (view in views) {
            expect_true(correction$views[[view]]$applied, info = view)
            expect_identical(correction$views[[view]]$k, 1L, info = view)
        }
        # A legacy design names no new-design field.
        expect_false("ruv_design" %in% names(correction))
        # KIT-184 K4 (gate G1): each data set's written factors are what its correction removed.
        for (view in views) expectIntegrationRemovedFactors(run$dir, view, correction$views[[view]]$k, input = input)
        # The design: the samples kept by the correlation filter, the biological sample column named (INT-2).
        expect_identical(manifest$sample_id_column, "Run")
        expect_identical(manifest$group_column, "group")
        expect_identical(manifest$technical_replicate_column, "tech_rep_group")
        expect_identical(manifest$biological_sample_column, "tech_rep_group")
        expected <- integrationOmicsExpected(omic, owner, final, input$samples$Run)
        # KIT-153: the data set's own matrix is the between-sample normalisation the RUV step read, on the corrected
        # matrix's features and runs; and RUV-III-C with the step's k and controls, run again on that normalised state,
        # gives the corrected matrix.
        binding <- omicsDeclarationBindings()[[omic]]
        normalised_state <- getWorkflowStateReference(owner, paste0(prefix, "_normalized"))
        normalised <- integrationOmicsExpected(omic, owner, normalised_state, input$samples$Run)
        again <- integrationOmicsQuiet(ruvIII_C_Varying(theObject = readWorkflowOmicsState(binding, owner, normalised_state),
            ruv_grouping_variable = "group", ruv_number_k = binding$extract_best_k(journey$ruv$optimisation),
            ctrl = binding$extract_controls(journey$ruv$optimisation)))
        id_column <- methods::slot(again, binding$id_slot)
        for (i in seq_along(views)) {
            view <- input$views[[views[[i]]]]
            expect_true(view$corrected, info = views[[i]])
            expect_equal(view$corrected_matrix, expected[[assays[[i]]]], tolerance = 1e-12, ignore_attr = TRUE, info = views[[i]])
            at <- match(view$matrix$feature_id, normalised[[assays[[i]]]]$feature_id)
            expect_false(anyNA(at), info = views[[i]])
            expect_identical(view$matrix$feature_id, view$corrected_matrix$feature_id, info = views[[i]])
            expect_equal(unname(as.matrix(view$matrix[-1L])), unname(as.matrix(normalised[[assays[[i]]]][at, -1L])), tolerance = 0,
                info = views[[i]])
            recomputed <- as.data.frame(methods::slot(again, binding$data_slot)[[assays[[i]]]], check.names = FALSE)
            at <- match(view$corrected_matrix$feature_id, as.character(recomputed[[id_column]]))
            expect_false(anyNA(at), info = views[[i]])
            expect_equal(unname(as.matrix(view$corrected_matrix[-1L])), unname(as.matrix(recomputed[at, input$samples$Run])),
                tolerance = 1e-10, info = views[[i]])
            features <- input$views[[views[[i]]]]$features
            expect_identical(features$assay, rep(assays[[i]], nrow(features)), info = views[[i]])
            # Internal standards are kept and marked (INT-3), here by the object's regex `^IS_` on the id.
            standard <- startsWith(features$feature_id, "IS_")
            expect_true(any(standard), info = views[[i]])
            expect_identical(features$is_internal_standard, ifelse(standard, "TRUE", "FALSE"), info = views[[i]])
        }
        expect_identical(manifest$internal_standards$regex, "^IS_")
    })
})

test_that("KIT-184 K4: lipidomics with today's RUV at two factors: the written factors are what was removed", {
    root <- withr::local_tempdir()
    fixture <- integrationOmicsOwner("lipidomics", file.path(root, "owner"))
    integrationOmicsQuiet(lipidDeclarationsJourney(fixture, ruv = utils::modifyList(lipidDeclarationsRuvParameters(),
        list(manual_k = 2L))))
    run <- integrationOmicsSummarySave("lipidomics", fixture, root)
    input <- readIntegrationInput(run$dir)
    features <- 0L
    rows <- 0L
    for (view in names(input$views)) {
        expect_identical(input$manifest$correction$views[[view]]$k, 2L, info = view)
        expectIntegrationRemovedFactors(run$dir, view, 2L, input = input)
        entry <- Filter(function(item) identical(item$name, view), input$manifest$views)[[1L]]
        features <- features + as.integer(entry$removed_factors$features)
        rows <- rows + input$views[[view]]$rows
    }
    # Features with a missing run were corrected one by one and are left out of the factors (gate G2).
    expect_lt(features, rows)
})

test_that("KIT-135: lipidomics and metabolomics views do not collide: both LCMS_Pos assays are qualified by the omic (INT-4)", {
    expect_identical(integrationInputViewName("lipidomics", "LCMS_Pos"), "lipidomics_LCMS_Pos")
    expect_identical(integrationInputViewName("metabolomics", "LCMS_Pos"), "metabolomics_LCMS_Pos")
    expect_identical(integrationInputViewName("metabolomics", "LC-MS Pos (C18)"), "metabolomics_LC-MS_Pos__C18_")
})

test_that("KIT-135: a new RUV design writes the RUV state's corrected table and names the design (M10); design none is normalised", {
    withr::local_options(MultiScholaR.ruv_design_panel = TRUE)
    root <- withr::local_tempdir()
    fixture <- integrationOmicsNormalised(integrationOmicsOwner("lipidomics", file.path(root, "owner")),
        itsd_ids = list(LCMS_Pos = "pos_L05", LCMS_Neg = "neg_L05"))
    owner <- fixture$owner
    parameters <- list(ruv_design = "groups_all_features", ruv_grouping_variable = "group", ruv_k = 1L, ruv_k_rule = "user",
        ruv_k_seed = NULL, ruv_control_ids = NULL)
    ruv <- integrationOmicsQuiet(runWorkflowOmicsReferenceRuv(lipidomicsDeclarationBinding, owner, parameters, list()))
    expect_identical(owner$current_state, "lipid_ruv_corrected")
    run <- integrationOmicsSummarySave("lipidomics", fixture, root)
    input <- readIntegrationInputWithMatrices(run$dir)
    manifest <- input$manifest
    expect_identical(manifest$matrix_role, "normalised")
    expect_identical(manifest$source_state$name, "lipid_ruv_corrected")
    correction <- manifest$correction
    expect_identical(correction[c("id", "method", "matrix", "da_tests", "ruv_design", "ruv_replicate_structure", "ruv_k_rule",
        "ruv_da_method")], list(id = "C4", method = "RUVIII", matrix = "ruv_corrected_values",
        da_tests = "uncorrected_values_with_factors_in_model", ruv_design = "groups_all_features", ruv_replicate_structure = "group",
        ruv_k_rule = "user", ruv_da_method = "limma_ruv_factors_exact_se"))
    expect_true(all(RUV_DESIGN_RECORD_FIELDS %in% names(correction)))
    expect_true(all(vapply(correction$views, `[[`, logical(1), "applied")))
    # The corrected matrix is the state's own (corrected) table, not the uncorrected assays the DA tests.
    state <- withWorkflowStateScope(owner, function() readWorkflowExecutionComponents(owner, ruv$reference,
        fields = lipidomicsDeclarationBinding$fields))
    expected <- integrationOmicsExpected("lipidomics", owner, "lipid_ruv_corrected", input$samples$Run)
    expect_equal(input$views$lipidomics_LCMS_Pos$corrected_matrix, expected$LCMS_Pos, tolerance = 1e-12, ignore_attr = TRUE)
    uncorrected <- as.data.frame(state@args$ruv_uncorrected_assays$LCMS_Pos, check.names = FALSE)
    expect_false(isTRUE(all.equal(unname(as.matrix(input$views$lipidomics_LCMS_Pos$corrected_matrix[input$samples$Run])),
        unname(as.matrix(uncorrected[input$samples$Run])))))
    # KIT-153: the data set's own matrix is the matrix that entered RUV: the uncorrected assays the RUV kernel kept for the
    # DA (its input), and the normalised state it ran on, in every assay.
    for (assay in c("LCMS_Pos", "LCMS_Neg")) {
        written <- input$views[[paste0("lipidomics_", assay)]]$matrix
        kept <- as.data.frame(state@args$ruv_uncorrected_assays[[assay]], check.names = FALSE)
        at <- match(written$feature_id, as.character(kept[[state@lipid_id_column]]))
        expect_false(anyNA(at), info = assay)
        expect_equal(unname(as.matrix(written[-1L])), unname(as.matrix(kept[at, input$samples$Run])), tolerance = 0, info = assay)
        normalised <- integrationOmicsExpected("lipidomics", owner, "lipid_normalized", input$samples$Run)[[assay]]
        expect_equal(unname(as.matrix(written[-1L])), unname(as.matrix(normalised[match(written$feature_id, normalised$feature_id), -1L])),
            tolerance = 0, info = assay)
    }
    # KIT-184 K4 (gate G1): each assay's written factors are what was removed; the largest principal angle to the
    # factors the DA keeps (computed on the features observed in every sample) is printed, not asserted.
    for (assay in c("LCMS_Pos", "LCMS_Neg")) {
        factors <- expectIntegrationRemovedFactors(run$dir, paste0("lipidomics_", assay), 1L, input = input)
        stored <- as.matrix(state@args$ruv_factors[[assay]]$W)
        if (identical(nrow(stored), nrow(factors))) {
            if (!is.null(rownames(stored)) && all(rownames(factors) %in% rownames(stored))) {
                stored <- stored[rownames(factors), , drop = FALSE]
            }
            basis <- qr.Q(qr(scale(stored, scale = FALSE)))
            cat(sprintf("KIT-184 G1 lipidomics_%s: largest sine of the principal angle to the DA's stored factors %.3g\n", assay,
                max(abs(basis - factors %*% crossprod(factors, basis)))))
        }
    }
    # The internal-standard step removes the standards it used (`remove_itsd_after_norm`), so the hand-selected one is
    # gone; the selection is recorded, and the regex's standards it did not use stay, marked.
    features <- input$views$lipidomics_LCMS_Pos$features
    expect_false("pos_L05" %in% features$feature_id)
    expect_true(manifest$internal_standards$selected_by_hand)
    expect_identical(features$is_internal_standard, ifelse(startsWith(features$feature_id, "IS_"), "TRUE", "FALSE"))

    # Design none: the same state name, no factor; the input is normalised, with no correction.
    none <- integrationOmicsNormalised(integrationOmicsOwner("lipidomics", file.path(root, "none")))
    integrationOmicsQuiet(runWorkflowOmicsReferenceRuv(lipidomicsDeclarationBinding, none$owner, list(ruv_design = "none",
        ruv_grouping_variable = "group", ruv_k = NULL, ruv_k_rule = NULL, ruv_k_seed = NULL, ruv_control_ids = NULL), list()))
    expect_identical(none$owner$current_state, "lipid_ruv_corrected")
    after <- readIntegrationInput(integrationOmicsSummarySave("lipidomics", none, file.path(root, "none-project"))$dir)$manifest
    expect_identical(after$matrix_role, "normalised")
    expect_null(after$correction)
    # KIT-184 K4: design none removes nothing: no factors.
    expect_true(all(vapply(after$views, function(view) is.null(view$removed_factors), logical(1))))
    expect_identical(after$source_state$name, "lipid_ruv_corrected")
})

test_that("KIT-135: without RUV the normalised matrix is written; with no final state nothing is, and an earlier input is removed", {
    root <- withr::local_tempdir()
    fixture <- integrationOmicsNormalised(integrationOmicsOwner("lipidomics", file.path(root, "owner")))
    owner <- fixture$owner
    expect_identical(owner$current_state, "lipid_normalized")
    input <- readIntegrationInput(integrationOmicsSummarySave("lipidomics", fixture, root)$dir)
    expect_identical(input$manifest$matrix_role, "normalised")
    expect_null(input$manifest$correction)
    expect_identical(input$manifest$source_state$name, "lipid_normalized")
    expect_false(input$manifest$internal_standards$selected_by_hand)

    empty <- integrationOmicsOwner("metabolomics", file.path(root, "empty"))
    project <- file.path(root, "empty-project")
    # An earlier save's input (K061-3): it must not stay to be read as this session's.
    stale <- file.path(project, "integration", "metabolomics")
    dir.create(stale, recursive = TRUE)
    writeLines("{}", file.path(stale, INTEGRATION_INPUT_MANIFEST))
    run <- integrationOmicsSummarySave("metabolomics", empty, project)
    # Today's RDS is written as before (here of no state); the input is not, and the earlier one is gone.
    expect_true(file.exists(run$rds))
    expect_false(dir.exists(run$dir))
})

test_that("KIT-135: the worker's integration export writes the same input beside its RDS", {
    skip_if_not(workflowOwnerCapable())
    # The in-test worker under load_all (test-workflow-omics-application-reads.R's `omicsReadMocks`).
    testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
        finishWorkflowExecutionJob = function(owner, status) { owner$job <- NULL; list(status = status) }, .package = "MultiScholaR")
    root <- withr::local_tempdir()
    fixture <- integrationOmicsOwner("metabolomics", file.path(root, "owner"))
    integrationOmicsQuiet(metabDeclarationsJourney(fixture))
    owner <- fixture$owner
    in_process <- readIntegrationInputWithMatrices(integrationOmicsSummarySave("metabolomics", fixture, file.path(root, "in-process"))$dir)
    integration <- file.path(root, "worker", "integration")
    dir.create(integration, recursive = TRUE)
    path <- file.path(integration, "metabolomics_kit135_final_s4.RDS")
    done <- integrationOmicsQuiet(omicsWorkerAction(fixture, omicsWorkerController(fixture), "metab_integration_export",
        list(path = path, state_name = owner$current_state)))
    expect_identical(done$acceptance, "complete")
    expect_s4_class(readRDS(path), "MetaboliteAssayData")
    worker <- readIntegrationInputWithMatrices(file.path(integration, "metabolomics"))
    expect_identical(names(worker$views), c("metabolomics_LCMS_Pos", "metabolomics_GCMS"))
    drop <- c("written_at", "software")
    expect_identical(worker$manifest[setdiff(names(worker$manifest), drop)], in_process$manifest[setdiff(names(in_process$manifest), drop)])
    expect_identical(worker$views, in_process$views)
    expect_identical(worker$samples, in_process$samples)
})
