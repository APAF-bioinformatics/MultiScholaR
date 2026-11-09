# KIT-105 W3: proteomics' normalisation back-references #10 to #22, their six module-file callees and the DIA_LIMPA
# pass-through moved, unchanged and with their names kept, into R/utils_workflow_dia_norm.R (design 2.2, 4.7). These
# blocks show where they live, that their definitions are the ones the module files held, and that their core
# callers (the normalisation and correlation actions, the pass-through, the legacy restore, the summary archive and
# the lipid/metab RUV declaration) give the output recorded before the move through them.

protNormCoreMoved <- c(
    "buildProtNormSkippedRuvResult", "resolveProtRuvParameters", "validateProtNormRuvK", "stripProtNormQcPlotOwners",
    "buildProtNormQcPlotProduct", "buildProtNormDensityPlot", "buildProtNormCorrelationPlot",
    "resolveProtNormCorrelationResultFilenames", "saveProtNormCorrelationResults", "saveProtNormExportMetadataFiles",
    "buildProtNormExportSummaryContent", "generateProtNormCompositeQcFigure", "finalizeProtNormWorkflowState",
    "buildProtNormManualRuvResult", "resolveProtNormCompositeFigureInputs", "generateProtNormCompositeFromFiles",
    "buildProtNormTitlePlot", "buildProtNormLabelPlot", "loadProtNormImageAsPlot",
    "runDIAReferenceLimpaNormalisationPassThrough"
)

# Recorded at aa79fe3f (before W3): the sha256 of each definition deparsed without source references.
protNormCoreBodyDigests <- c(
    buildProtNormSkippedRuvResult = "ffcb7b36bb44464f316a128153cdfa907751c53382560291f82778d2aa92790c",
    resolveProtRuvParameters = "25f9253c1c4b9013cb6e2e3b3ae3659c78a61034c7b1179ddec5cf3204c0dfac",
    validateProtNormRuvK = "85dcecdc2ee46797f765ba734d4a62c5e185be2bb4ad5dd469b240eeae1d9002",
    stripProtNormQcPlotOwners = "9e678b9accd43a77d40d7d25995f72736ff8597590c3451a3131162b0c4a3a7c",
    buildProtNormQcPlotProduct = "d55fad6feaaa048a5e87da556b282d6b8b3a87d87d15f186bc2b793be6b4dce1",
    buildProtNormDensityPlot = "4c7d554497d04192405497f345627ea5d5a4a88d03740f67a97d3163f39fb704",
    buildProtNormCorrelationPlot = "ab28205cedb12f6c58867c055a6a309f7abe023d49fef596cdf347836491f843",
    resolveProtNormCorrelationResultFilenames = "55944dc746f87220718ccf901be206cb4db0c5186df5081736a71f766e10a7e6",
    saveProtNormCorrelationResults = "4b3a2391cea3ffa07f61d865eb60a872402a4862ec1ceee9970862c49348a58c",
    saveProtNormExportMetadataFiles = "51b26cb0b588270032a2f22f98350dbfebf324eed5b5ab76745daae45bab28f3",
    buildProtNormExportSummaryContent = "6d9f6dfae03215a14af70506e13000499f519a2338c57414cb109e0461d9f1c9",
    generateProtNormCompositeQcFigure = "e8a3c81e06f2ee660b1b48e365ce5025e23503c1b746ef070093375ec0547df7",
    finalizeProtNormWorkflowState = "6348516fc7f8f3adfef0c9250aa6e92a616a0a9d6a5710992b1ed49bb8f93f89",
    buildProtNormManualRuvResult = "f83b5129eee4bae937f4c2462d229c187fa87d93326f25ec12058108ede24595",
    resolveProtNormCompositeFigureInputs = "9771a4bee6c40b1fb647d38cfe944b4cbd5d2c664f200b6cb17c25d56e736944",
    generateProtNormCompositeFromFiles = "1ad17f33d5ed25d7928653b7c6636e80a99594dd02bba89987ecdced9b4ce72e",
    buildProtNormTitlePlot = "5993ca10f9e4080ef7ec0b14dfb4de7e56628c532ca0a2eda8b726b83441c538",
    buildProtNormLabelPlot = "3f45f763cb05b83ed873e523c33d019800e5719c06c1cbf6b771b98d8fa9bf69",
    loadProtNormImageAsPlot = "6b9882f8594382fdabaf741c05fad43452517a5a37fdd07cad65b2a87e02e7ff",
    runDIAReferenceLimpaNormalisationPassThrough = "63aea4d328e48d7d42e5230e89d05cea930220b8f389942017c8f3c9bb607729"
)

# Recorded at aa79fe3f (before W3) through the core callers on the fixtures below.
protNormCoreOutputDigests <- c(
    normalization_skip = "d21daf95829cd776ee2495206c2de1f40eec75da29c1396ec92db7142e2fd6cc",
    normalization_manual = "4f46dfedb8dbdd3a05f5e95c65f94b7f723f11b5f101296a0c2b13e27f2a4b8c",
    correlation = "14c84122c28b8bb98b1919009d93fefbd6cf70a0d25daaa76b11e43b7581dd08",
    pass_through = "415313cc2e24e23303f2e4283028a5a51af954102859e1861dd1ffcaf59923d9",
    legacy_restore = "453dfbdd95c8eb8ab856f72c9ca131e091c43d026848cf2115bdbbfeb46856da",
    summary_archive = "cb4e0787fc2fb873b79f7c8858e6ac8671948894b85310e72583f1bed5289e41",
    lipid_ruv = "0921c3d1d2dae12528e3b3690680e84dc1d9749fea2226a419aeade1b11390f6",
    metab_ruv = "dea4b04738e8b80a58bbd8639b8eef9eb90ae27d6c6da30da0d434592f9d10df"
)

# MULTISCHOLAR_KIT105_W3_DUMP=1 prints each digest (how the values above were recorded).
protNormCoreExpectDigest <- function(key, value, pinned = protNormCoreOutputDigests) {
    found <- digest::digest(value, algo = "sha256")
    if (nzchar(Sys.getenv("MULTISCHOLAR_KIT105_W3_DUMP"))) {
        cat(sprintf("\nKIT105W3 %s = \"%s\"\n", key, found))
        parts <- rapply(value, function(x) digest::digest(x, algo = "sha256"), how = "unlist")
        cat(sprintf("PART %s %s %s\n", key, names(parts), substr(parts, 1, 12)), sep = "")
    }
    expect_identical(found, pinned[[key]], info = key)
}

protNormCoreSourceFile <- function(name) {
    ns <- asNamespace("MultiScholaR")
    if (!exists(name, envir = ns, inherits = FALSE)) return(NA_character_)
    file <- utils::getSrcFilename(get(name, envir = ns, inherits = FALSE))
    if (length(file)) file else NA_character_
}

protNormCoreDefinitionDigest <- function(name) {
    digest::digest(deparse(utils::removeSource(get(name, envir = asNamespace("MultiScholaR")))), algo = "sha256")
}

# Counts the calls made by name to namespace functions while `code` runs.
protNormCoreCountCalls <- function(names, code) {
    calls <- stats::setNames(as.list(integer(length(names))), names)
    counters <- lapply(names, function(name) {
        original <- get(name, envir = asNamespace("MultiScholaR"))
        function(...) {
            calls[[name]] <<- calls[[name]] + 1L
            original(...)
        }
    })
    names(counters) <- names
    value <- rlang::inject(testthat::with_mocked_bindings(code, !!!counters, .package = "MultiScholaR"))
    list(value = value, calls = unlist(calls))
}

protNormCoreQuiet <- function(code) suppressWarnings(suppressMessages(code))

# Files under `root` (relative names), with the md5 of each text table; images, PDFs and RDS files by name only.
protNormCoreFiles <- function(root) {
    files <- sort(list.files(root, recursive = TRUE))
    files <- files[!grepl("preflight", files, fixed = TRUE)]
    tables <- files[grepl("[.](tsv|txt)$", files)]
    list(files = files, tables = stats::setNames(unname(tools::md5sum(file.path(root, tables))), tables))
}

protNormCoreWithoutPlots <- function(value) {
    value[!vapply(value, function(x) inherits(x, c("gg", "ggplot", "grob", "gtable")), logical(1))]
}

# A native DIA owner (the iq route) on the "small" protein recipe of test-workflow-checkpoint-dia-state.R, with
# the paths and the plain-environment workflow the GUI's actions read.
protNormCoreOwner <- function(envir = parent.frame()) {
    withr::local_options(multischolar.dia_native_sessions = TRUE, multischolar.capture_test_checkpoints = FALSE,
        .local_envir = envir)
    root <- withr::local_tempdir(.local_envir = envir)
    paths <- list(source_dir = file.path(root, "scripts", "proteomics"),
        protein_qc_dir = file.path(root, "results", "proteomics", "protein_qc"),
        peptide_qc_dir = file.path(root, "results", "proteomics", "peptide_qc"))
    paths$publication_graphs_dir <- file.path(root, "results", "proteomics", "publication_graphs")
    paths$time_dir <- file.path(paths$publication_graphs_dir, "filtering_qc")
    for (path in unlist(paths)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    session <- new.env(parent = emptyenv())
    session$userData <- new.env(parent = emptyenv())
    session$onFlush <- session$onFlushed <- session$onSessionEnded <- function(...) function() invisible(NULL)
    original <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    workflow <- new.env(parent = emptyenv())
    workflow$state_manager <- WorkflowState$new(audit_enabled = FALSE)
    workflow$design_matrix <- original@design_matrix
    workflow$config_list <- original@args
    initializeDIANativeImport(workflow, paths, session, "diann", "iq")
    manager <- workflow$state_manager
    withr::defer(closeDIAStateAdapter(manager), envir = envir)
    manager$saveState("protein_replicate_filtered", original, list(), "source")
    list(root = root, paths = paths, workflow = workflow, manager = manager, original = original)
}

protNormCoreSettings <- function(ruv_mode) {
    list(norm_method = "quantile", ruv_mode = ruv_mode, auto_percentage_min = NULL, auto_percentage_max = NULL,
        separation_metric = NULL, k_penalty_weight = NULL, adaptive_k_penalty = NULL, ruv_percentage = 50, ruv_k = 1,
        color_variable = "group", shape_variable = NULL)
}

# A saved state's science (its slots without `args`) in place of its entity id, which is local to the store.
protNormCoreState <- function(manager, reference) {
    if (is.null(reference)) return(NULL)
    withWorkflowStateScope(manager, function() {
        object <- resolveDIAStateValue(manager$.context, reference)$data
        if (!isS4(object)) object <- rebuildDIACheckpointObject(manager$.context, object, TRUE, list())
        slots <- setdiff(methods::slotNames(object), "args")
        lapply(stats::setNames(nm = slots), function(slot) methods::slot(object, slot))
    })
}

protNormCoreNormState <- function(norm, manager) {
    table <- if (!is.null(norm$ruv_optimization_table)) materializeDIATable(manager$.context, norm$ruv_optimization_table)
    list(
        normalized = protNormCoreState(manager, norm$normalized_protein_obj),
        final = protNormCoreState(manager, norm$ruv_normalized_obj),
        input = protNormCoreState(manager, norm$normalization_input_reference), best_k = norm$best_k,
        controls = norm$control_genes_index, ruv = protNormCoreWithoutPlots(norm$ruv_optimization_result),
        table = table, complete = norm$normalization_complete,
        ruv_complete = norm$ruv_complete, text = norm$filtering_summary_text,
        filtering_plot = class(norm$post_norm_filtering_plot), qc_paths = norm$qc_plot_paths,
        state = manager$current_state, history = unlist(manager$state_history)
    )
}

test_that("back-references #10 to #22, their callees and the pass-through are defined in the proteomics norm core", {
    if (is.null(utils::getSrcref(runDIAReferenceNormalization))) {
        testthat::skip("source references are kept under load_all only")
    }
    found <- vapply(protNormCoreMoved, protNormCoreSourceFile, character(1))
    expect_identical(found, stats::setNames(rep("utils_workflow_dia_norm.R", length(protNormCoreMoved)), protNormCoreMoved))
    # The module keeps its adapters and readers in its own files.
    expect_identical(protNormCoreSourceFile("resolveProtNormRuvParameters"), "mod_prot_norm_ruv_helpers.R")
    expect_identical(protNormCoreSourceFile("runProtNormNormalizationWorkflow"), "mod_prot_norm_workflow_helpers.R")
})

test_that("the moved definitions are the module files' definitions, in the namespace and not exported", {
    ns <- asNamespace("MultiScholaR")
    for (name in protNormCoreMoved) {
        expect_identical(environment(get(name, envir = ns, inherits = FALSE)), ns, info = name)
    }
    found <- vapply(protNormCoreMoved, protNormCoreDefinitionDigest, character(1))
    if (nzchar(Sys.getenv("MULTISCHOLAR_KIT105_W3_DUMP"))) {
        cat("\n", sprintf("    %s = \"%s\",\n", names(found), found), sep = "")
    }
    expect_identical(found, protNormCoreBodyDigests[protNormCoreMoved])
    expect_false(any(protNormCoreMoved %in% getNamespaceExports("MultiScholaR")))
})

test_that("the normalisation action and the correlation action give the recorded output through the moved functions", {
    fixture <- protNormCoreOwner()
    workflow <- fixture$workflow
    manager <- fixture$manager
    norm <- new.env(parent = emptyenv())
    aesthetics <- function() list(color_var = "group", shape_var = NULL)
    post <- function(reference) invisible(NULL)
    run <- function(mode) {
        protNormCoreCountCalls(c("buildProtNormSkippedRuvResult", "resolveProtRuvParameters", "stripProtNormQcPlotOwners",
            "generateProtNormCompositeQcFigure", "resolveProtNormCompositeFigureInputs",
            "generateProtNormCompositeFromFiles", "finalizeProtNormWorkflowState"),
            protNormCoreQuiet(runDIAReferenceNormalization(protNormCoreSettings(mode), workflow, norm, fixture$paths,
                "proteomics", "core", post, post, function() "group")))
    }

    skipped <- run("skip")
    expect_true(skipped$value)
    expect_identical(skipped$calls, c(buildProtNormSkippedRuvResult = 1L, resolveProtRuvParameters = 0L,
        stripProtNormQcPlotOwners = 0L, generateProtNormCompositeQcFigure = 1L, resolveProtNormCompositeFigureInputs = 1L,
        generateProtNormCompositeFromFiles = 1L, finalizeProtNormWorkflowState = 1L))
    protNormCoreExpectDigest("normalization_skip", list(protNormCoreNormState(norm, manager),
        protNormCoreFiles(fixture$root)))

    manual <- run("manual")
    expect_true(manual$value)
    expect_identical(manual$calls[c("buildProtNormSkippedRuvResult", "resolveProtRuvParameters",
        "generateProtNormCompositeQcFigure", "finalizeProtNormWorkflowState")],
        c(buildProtNormSkippedRuvResult = 0L, resolveProtRuvParameters = 1L, generateProtNormCompositeQcFigure = 1L,
            finalizeProtNormWorkflowState = 1L))
    expect_true(inherits(norm$ruv_optimization_result$best_cancor_plot, "grob"))
    protNormCoreExpectDigest("normalization_manual", list(protNormCoreNormState(norm, manager),
        protNormCoreFiles(fixture$root)))

    correlation <- protNormCoreCountCalls(c("buildProtNormQcPlotProduct", "buildProtNormDensityPlot",
        "buildProtNormCorrelationPlot", "stripProtNormQcPlotOwners", "resolveProtNormCorrelationResultFilenames"),
        protNormCoreQuiet(runDIAReferenceCorrelation(workflow, norm, fixture$paths, 0.5, FALSE, "proteomics", "core",
            function() "group", aesthetics)))
    expect_true(isSavedWorkflowReference(correlation$value))
    expect_identical(correlation$calls[c("buildProtNormQcPlotProduct", "stripProtNormQcPlotOwners",
        "resolveProtNormCorrelationResultFilenames")],
        # The owner stripper is counted with its recursive call (it recurses by name), as recorded before the move.
        c(buildProtNormQcPlotProduct = 1L, stripProtNormQcPlotOwners = 2L, resolveProtNormCorrelationResultFilenames = 1L))
    expect_identical(norm$final_qc_plot$plot_env, rlang::base_env())
    protNormCoreExpectDigest("correlation", list(
        final = protNormCoreState(manager, correlation$value),
        vector = hydrateSavedWorkflowEntity(manager$.context$store, norm$correlation_vector)$value,
        threshold = norm$correlation_threshold, state = manager$current_state,
        plot_data = as.data.frame(norm$final_qc_plot$data), title = norm$final_qc_plot$labels$title,
        filtering_plot = class(norm$final_filtering_plot), files = protNormCoreFiles(fixture$root)))
})

test_that("the DIA_LIMPA pass-through designates the current state as the correlation input, as recorded", {
    fixture <- protNormCoreOwner()
    norm <- new.env(parent = emptyenv())
    seen <- NULL
    value <- protNormCoreQuiet(runDIAReferenceLimpaNormalisationPassThrough(protNormCoreSettings("covariates"),
        fixture$workflow, norm, function(reference) seen <<- reference))
    expect_true(value)
    current <- getDIAStateReference(fixture$manager)
    expect_identical(seen$entity_id, current$entity_id)
    for (field in c("normalization_input_reference", "normalized_protein_obj", "ruv_normalized_obj")) {
        expect_identical(norm[[field]]$entity_id, current$entity_id, info = field)
    }
    protNormCoreExpectDigest("pass_through", list(protNormCoreState(fixture$manager, current),
        mget(c("ruv_optimization_result", "best_k", "control_genes_index", "normalization_complete",
            "filtering_summary_text"), envir = norm),
        fixture$manager$current_state))
})

test_that("the legacy restore writes the recorded correlation table and metadata files through the moved functions", {
    root <- withr::local_tempdir()
    original <- makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small"))
    value <- list(current_s4_object = original, ruv_applied = FALSE, qc_params = list(min = 1),
        accession_cleanup_results = list(cleaned = 2L), limpa_parameters = list(quantification = "dpcquant"),
        mixed_species_analysis = list(enabled = FALSE))
    captured <- NULL
    restored <- protNormCoreCountCalls(c("saveProtNormCorrelationResults", "resolveProtNormCorrelationResultFilenames",
        "saveProtNormExportMetadataFiles"),
        testthat::with_mocked_bindings(
            protNormCoreQuiet(restoreDIANativeLegacyArtifacts(value, NULL, NULL, root)),
            restoreDIANativeGuiFiles = function(request, paths) {
                captured <<- vapply(request$dependencies, function(entry) entry$destination, character(1))
                invisible(TRUE)
            }))
    expect_identical(restored$calls, c(saveProtNormCorrelationResults = 1L,
        resolveProtNormCorrelationResultFilenames = 1L, saveProtNormExportMetadataFiles = 1L))
    contents <- lapply(stats::setNames(nm = list.files(file.path(root, "legacy-source"))),
        function(file) readRDS(file.path(root, "legacy-source", file)))
    protNormCoreExpectDigest("legacy_restore", list(destinations = captured, metadata = contents,
        table = unname(tools::md5sum(list.files(file.path(root, "legacy-protein-qc"), full.names = TRUE)))))
})

test_that("the summary archive writes the recorded metadata files and summary text through the moved functions", {
    source_dir <- withr::local_tempdir()
    metadata <- list(final_protein_count = 120L, final_sample_count = 12L,
        contrasts_tbl = data.frame(friendly_names = c("B_vs_A", "C_vs_A")), normalization_method = "cyclicloess",
        ruv_mode = "automatic", ruv_k = 2L, correlation_threshold = 0.5, workflow_type = "DIA",
        report_template = "default", qc_params = list(min = 1), ruv_optimization_result = list(best_k = 2L))
    manager <- new.env(parent = emptyenv())
    manager$.context <- list(owned_export = list(metadata_reference = list(entity_id = "metadata"),
        generation = "earlier", extra_references = list()), head = list(entity_id = "current"))
    manager$current_state <- "correlation_filtered"
    manager$state_history <- list("protein_replicate_filtered", "ruv_corrected", "correlation_filtered")
    workflow <- list(state_manager = manager, da_ui_params = list(q = 0.05))
    stamp <- as.POSIXct("2026-11-08 12:00:00", tz = "UTC")
    archived <- protNormCoreCountCalls(c("saveProtNormExportMetadataFiles", "buildProtNormExportSummaryContent"),
        testthat::with_mocked_bindings(
            protNormCoreQuiet(refreshDIANativeSummaryArchive(workflow, source_dir)),
            assertDIAExportReady = function(...) invisible(TRUE),
            resolveDIAStateValue = function(context, reference) metadata,
            collectDIANativeGuiFiles = function(...) list(),
            runDIANativeSessionExport = function(manager, metadata, source_dir, extra_references, dependencies,
                prepare_outputs) {
                metadata$export_timestamp <- stamp
                prepare_outputs(metadata, list(archive_name = "filtered_session_data_core.zip"))
                list(status = "complete", state = metadata$r6_current_state_name)
            }))
    expect_identical(archived$calls, c(saveProtNormExportMetadataFiles = 1L, buildProtNormExportSummaryContent = 1L))
    expect_identical(archived$value$state, "correlation_filtered")
    files <- sort(list.files(source_dir))
    protNormCoreExpectDigest("summary_archive", list(files = files,
        rds = lapply(stats::setNames(nm = files[grepl("[.]RDS$", files)]), function(file) readRDS(file.path(source_dir, file))),
        summary = readLines(file.path(source_dir, "filtered_session_summary.txt"))))
})

# The lipid/metab RUV declaration strips the RUV plots' owners with the moved function (omicsRuvOperation,
# utils_workflow_omics_declarations.R:295); its corrected state and the saved optimisation are as recorded.
protNormCoreOmicsRuv <- function(omic) {
    root <- withr::local_tempdir(.local_envir = parent.frame())
    fixture <- if (omic == "lipid") lipidDeclarationsOwner(root) else metabDeclarationsOwner(root)
    owner <- fixture$owner
    withr::defer(closeDIAStateAdapter(owner), envir = parent.frame())
    config <- fixture$workflow$config_list
    counted <- protNormCoreQuiet(if (omic == "lipid") {
        runLipidReferenceIntensityFilter(owner, 1, 0.5, config, hydrate = FALSE)
        runLipidReferenceDuplicateResolution(owner, config, hydrate = FALSE)
        commitLipidReferenceStoreStep(owner, "lipid_qc_finalize", "lipid_qc_complete", config, "QC processing complete - ready for normalization")
        runLipidReferenceItsdNormalisation(owner, "median", NULL, config, hydrate = FALSE)
        runLipidReferenceLog2(owner, 1, config, hydrate = FALSE)
        runLipidReferenceBetweenSampleNormalisation(owner, "cyclicloess", config, hydrate = FALSE)
        protNormCoreCountCalls("stripProtNormQcPlotOwners",
            runLipidReferenceRuv(owner, lipidDeclarationsRuvParameters(), config, hydrate = FALSE))
    } else {
        runMetabReferenceIntensityFilter(owner, 1, 0.5, config, hydrate = FALSE)
        runMetabReferenceDuplicateResolution(owner, config, hydrate = FALSE)
        commitWorkflowOmicsStoreStep(owner, "metab_qc_finalize", "metab_qc_complete", config, "QC processing complete - ready for normalization")
        runMetabReferenceItsdNormalisation(owner, "median", NULL, config, hydrate = FALSE)
        runMetabReferenceLog2(owner, 1, config, hydrate = FALSE)
        runMetabReferenceBetweenSampleNormalisation(owner, "cyclicloess", config, hydrate = FALSE)
        protNormCoreCountCalls("stripProtNormQcPlotOwners",
            runMetabReferenceRuv(owner, metabDeclarationsRuvParameters(), config, hydrate = FALSE))
    })
    ruv <- counted$value
    plots <- 0L
    optimisation <- lapply(ruv$optimisation, function(result) lapply(result, function(field) {
        if (!inherits(field, c("gg", "ggplot"))) return(field)
        plots <<- plots + 1L
        expect_identical(field$plot_env, rlang::base_env())
        list(data = as.data.frame(field$data), labels = field$labels)
    }))
    expect_true(plots > 0L)
    expect_identical(counted$calls[["stripProtNormQcPlotOwners"]], plots)
    list(state = owner$current_state, final = protNormCoreState(owner, ruv$value$final),
        values = ruv$value[c("best_k", "control_count", "percentage")], optimisation = optimisation)
}

test_that("the lipid and metab RUV declaration give the recorded result through the moved owner stripper", {
    skip_if_not(workflowOwnerCapable())
    protNormCoreExpectDigest("lipid_ruv", protNormCoreOmicsRuv("lipid"))
    protNormCoreExpectDigest("metab_ruv", protNormCoreOmicsRuv("metab"))
})

# ------------------------------------------------------------------------------------------------
# KIT-111 W1 (a): the three normalisation-tab runners (runDIAReferenceNormalization, runDIAReferenceCorrelation,
# runDIAReferenceRuvCovariates) and the DIA_LIMPA pass-through, in process at core level, captured on the start
# commit (7c7c599e) before the split (W2) and the dispatch (W3), and frozen from two byte-identical dumps made in
# separate processes (docs/qa/workflow-execution/kit/kit-111/w1/test-first.txt). They pass before and after; W2 and
# W3 compare the worker path with these values. Per case, after each click: the committed states (their science,
# their args' limpa and RUV covariates entries, config, description, previous state and product names), the head's
# history and current state, the head's products fingerprinted by type, the step provenance and the analysis
# record's steps (store-local ids, times and builds masked), the completed outputs, the published files,
# protNormCoreNormState, every norm_data field, the workflow's RUV result without plots and its other fields, the
# runner's value or refusal, the callbacks' inputs and the call counts of the halves' callees.
# MULTISCHOLAR_KIT111_W1_DUMP=<file> appends each case's capture as R code (how the frozen values were made).
# ------------------------------------------------------------------------------------------------

protNormW1Callees <- c("assertDIAExportReady", "runSavedOperation", "commitDIAReferenceState", "saveDIADisplayTable",
    "publishDIAActionOutputs", "readDIAProteinProgress", "updateProteinFiltering", "writeDIADesignTable",
    "resolveDIAReferenceQcPlot", "resolveProtNormCorrelationResultFilenames", "generateProtNormCompositeQcFigure",
    "finalizeProtNormWorkflowState")

protNormW1Sha <- function(x) digest::digest(x, algo = "sha256")

#' Store-local ids and the fixture's temporary root masked in every character value of plain lists.
protNormW1Mask <- function(x, root) {
    if (is.character(x)) {
        x[!is.na(x) & grepl("^[a-f0-9]{64}$", x)] <- "<id>"
        return(gsub(root, "<root>", x, fixed = TRUE))
    }
    if (is.list(x) && !is.object(x)) {
        x[] <- lapply(x, protNormW1Mask, root = root)
        if (!is.null(names(x))) names(x) <- protNormW1Mask(names(x), root)
    }
    x
}

#' A saved checkpoint's science: protNormCoreState's slots, with the limpa and RUV covariates entries of its args.
protNormW1State <- function(manager, reference) {
    withWorkflowStateScope(manager, function() {
        object <- resolveDIAStateValue(manager$.context, reference)$data
        if (!isS4(object)) object <- rebuildDIACheckpointObject(manager$.context, object, TRUE, list())
        slots <- setdiff(methods::slotNames(object), "args")
        list(slots = lapply(stats::setNames(nm = slots), function(slot) methods::slot(object, slot)),
            args = object@args[intersect(c("limpa_dpc_quant_results", "ruvCovariates"), names(object@args))])
    })
}

#' A value without store-local ids: a checkpoint by its science, a table by its rows, any other saved entity by its
#' hydrated envelope; plots by class (a ggplot also by its data and title); bookkeeping fields dropped.
protNormW1Canonical <- function(value, manager, root, depth = 0L) {
    if (depth > 32L) return("<depth>")
    if (is.null(value)) return(NULL)
    if (is.function(value) || is.environment(value)) return("<opaque>")
    if (inherits(value, c("gg", "ggplot"))) return(list(.gg = class(value),
        data = protNormW1Sha(as.data.frame(value$data)), title = value$labels$title))
    if (inherits(value, c("grob", "gtable", "htmlwidget", "recordedplot"))) return(paste0("<", class(value)[[1L]], ">"))
    if (isSavedWorkflowReference(value)) {
        if (identical(value$kind, "checkpoint")) return(list(.state = protNormW1Sha(protNormW1State(manager, value))))
        if (identical(value$payload$codec, "parquet_duckdb")) {
            table <- as.data.frame(materializeWorkflowTable(manager$.context, value))
            attr(table, "dia_source_order") <- NULL
            rownames(table) <- NULL
            return(list(.table = protNormW1Sha(table)))
        }
        if (identical(value$kind, "rendered_asset")) return(list(.rendered = value$kind))
        hydrated <- hydrateSavedWorkflowEntity(manager$.context$store, value)
        return(list(.entity = value$kind, value = protNormW1Canonical(hydrated, manager, root, depth + 1L)))
    }
    if (isS4(value)) return(list(.s4 = class(value)[[1L]]))
    if (is.data.frame(value)) {
        value <- as.data.frame(value)
        rownames(value) <- NULL
        return(list(.df = protNormW1Sha(protNormW1Mask(as.list(value), root)), dim = dim(value)))
    }
    if (is.list(value)) {
        keep <- !(names(value) %||% rep("", length(value))) %in% c("timestamp", "mtime", "dependencies",
            "product_dependencies", "owned_components", "component_dependencies", "committed_at")
        return(lapply(value[keep], protNormW1Canonical, manager = manager, root = root, depth = depth + 1L))
    }
    protNormW1Mask(value, root)
}

#' The head's products (and its output products), each by its type and a fingerprint of its canonical value; a
#' rendered asset by its kind only.
protNormW1Products <- function(manager, root) {
    snapshot <- manager$.context$snapshot
    one <- function(reference, prefix) {
        if (identical(reference$kind, "rendered_asset")) return(paste0(prefix, "rendered_asset"))
        canonical <- protNormW1Canonical(reference, manager, root)
        type <- canonical$value$product_type
        if (!is.character(type)) type <- reference$kind
        paste0(prefix, type, ":", substr(protNormW1Sha(canonical), 1L, 16L))
    }
    sort(c(vapply(unname(snapshot$products), one, "", prefix = ""),
        vapply(unname(snapshot$output_products), one, "", prefix = "output/")))
}

#' One committed state's record: its science, config, description, previous state and product names.
protNormW1Record <- function(manager, reference, root) {
    record <- readDIAStateMetadata(manager, reference)
    list(state = protNormW1State(manager, reference), config = protNormW1Canonical(record$config, manager, root),
        description = record$description, previous = record$previous_state,
        products = sort(names(record$data$products)))
}

#' The step provenance without its times, builds and commits.
protNormW1Provenance <- function(snapshot, root) {
    protNormW1Mask(lapply(snapshot$step_provenance, function(entry) entry[setdiff(names(entry),
        c("committed_at", "build", "source_commit"))]), root)
}

#' The analysis record's steps without ids, times and software, each product by its kind and type.
protNormW1RecordSteps <- function(manager, root) {
    steps <- tryCatch(buildWorkflowAnalysisRecord(manager)$steps, error = function(e) paste("error:", conditionMessage(e)))
    if (is.character(steps)) return(steps)
    protNormW1Mask(lapply(steps, function(step) {
        step$products <- lapply(step$products, function(product) product[c("kind", "product_type")])
        step[setdiff(names(step), c("committed_at", "result_entity_id", "input_entity_ids", "software"))]
    }), root)
}

#' Run one click with the halves' callees counted; a refusal is returned as its message.
protNormW1Run <- function(code) {
    calls <- stats::setNames(integer(length(protNormW1Callees)), protNormW1Callees)
    counters <- lapply(protNormW1Callees, function(name) {
        original <- get(name, envir = asNamespace("MultiScholaR"))
        function(...) {
            calls[[name]] <<- calls[[name]] + 1L
            original(...)
        }
    })
    names(counters) <- protNormW1Callees
    outcome <- rlang::inject(testthat::with_mocked_bindings(
        tryCatch(list(value = protNormCoreQuiet(code)), error = function(e) list(error = conditionMessage(e))),
        !!!counters, .package = "MultiScholaR"))
    c(outcome, list(calls = calls))
}

#' The capture after one click (see the section's header).
protNormW1Capture <- function(fixture, norm, run, seen) {
    force(run)
    manager <- fixture$manager
    root <- fixture$root
    # The head as an action scope reads it (the context's copy outside a scope can lag the last commit).
    snapshot <- withWorkflowStateScope(manager, function() manager$.context$snapshot)
    history <- as.character(unlist(snapshot$state_history))
    committed <- setdiff(history, "initial")
    # The workflow's fields but its owners and its state_update_trigger (a time).
    workflow <- mget(setdiff(sort(ls(fixture$workflow)), c("state_manager", "design_matrix", "config_list",
        "state_update_trigger")), envir = fixture$workflow)
    value <- if (isSavedWorkflowReference(run$value)) protNormW1Canonical(run$value, manager, root) else run$value
    list(
        current = manager$current_state,
        history = as.character(unlist(manager$state_history)),
        head_current = snapshot$current_state,
        head_history = history,
        states = vapply(stats::setNames(nm = committed), function(state)
            protNormW1Sha(protNormW1Record(manager, snapshot$states[[state]], root)), ""),
        outputs = lapply(stats::setNames(nm = committed), function(state)
            unlist(snapshot$completed_outputs[[snapshot$states[[state]]$entity_id]])),
        products = protNormW1Products(manager, root),
        provenance = protNormW1Sha(protNormW1Provenance(snapshot, root)),
        record = protNormW1Sha(protNormW1RecordSteps(manager, root)),
        files = protNormW1Sha(protNormCoreFiles(root)),
        norm_state = protNormW1Sha(protNormW1Mask(protNormCoreNormState(norm, manager), root)),
        norm_fields = vapply(mget(sort(ls(norm)), envir = norm), function(field)
            substr(protNormW1Sha(protNormW1Canonical(field, manager, root)), 1L, 16L), ""),
        ruv_result = protNormW1Sha(protNormW1Mask(protNormCoreWithoutPlots(fixture$workflow$ruv_optimization_result), root)),
        workflow = protNormW1Sha(protNormW1Canonical(workflow, manager, root)),
        value = protNormW1Sha(value),
        error = if (is.null(run$error)) NA_character_ else protNormW1Mask(run$error, root),
        calls = run$calls,
        callbacks = seen)
}

protNormW1DumpPath <- function() {
    path <- Sys.getenv("MULTISCHOLAR_KIT111_W1_DUMP")
    if (nzchar(path)) path else NULL
}

#' Compare one case with its frozen value; under the dump switch, append it as R code first.
protNormW1Expect <- function(case, capture) {
    path <- protNormW1DumpPath()
    if (!is.null(path)) {
        fields <- vapply(names(capture), function(field) paste0("        ", field, " = ", paste(deparse(capture[[field]],
            width.cutoff = 100L, control = c("keepNA", "keepInteger", "niceNames", "showAttributes")),
            collapse = "\n            ")), "")
        cat(sprintf("    %s = list(\n%s),\n", case, paste(fields, collapse = ",\n")), file = path, append = TRUE)
    }
    frozen <- protNormW1Frozen()[[case]]
    expect_false(is.null(frozen), info = case)
    for (field in union(names(frozen), names(capture))) {
        expect_identical(capture[[field]], frozen[[field]], info = paste(case, field))
    }
}

#' A DIA_LIMPA owner (the KIT-053 record fixture of test-workflow-analysis-record.R): the "small" recipe rolled up and
#' quantified by the limpa DPC engine, then the protein QC stage's between-sample normalisation (none) and peptide
#' count filter, committed as the route commits them, with protNormCoreOwner's paths.
protNormW1LimpaOwner <- function(envir = parent.frame()) {
    withr::local_options(multischolar.dia_native_sessions = TRUE, multischolar.capture_test_checkpoints = FALSE,
        .local_envir = envir)
    root <- withr::local_tempdir(.local_envir = envir)
    paths <- list(source_dir = file.path(root, "scripts", "proteomics"),
        protein_qc_dir = file.path(root, "results", "proteomics", "protein_qc"),
        peptide_qc_dir = file.path(root, "results", "proteomics", "peptide_qc"))
    paths$publication_graphs_dir <- file.path(root, "results", "proteomics", "publication_graphs")
    paths$time_dir <- file.path(paths$publication_graphs_dir, "filtering_qc")
    for (path in unlist(paths)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(paths$source_dir, "kit111-limpa", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager), envir = envir)
    manager$.context$reference_states <- TRUE
    manager$setWorkflowType("DIA_LIMPA")
    manager$.context$experiment_paths <- paths
    rolled <- rollUpPrecursorToPeptide(makeWorkflowCheckpointPeptide(readWorkflowCheckpointRecipe("small")), core_utilisation = NA)
    quantified <- protNormCoreQuiet(proteinMissingValueImputationLimpa(log2TransformPeptideMatrix(calcPeptideMatrix(rolled)),
        dpc_slope = 0.8, verbose = FALSE, chunk = 1000, engine = "limpa"))
    setDIAStepProvenance(manager, "design_builder")
    manager$saveState("protein_dpc_quant", quantified, list(), "fixture")
    workflow <- list2env(list(state_manager = manager, qc_params = list(protein_qc = list()), config_list = list()))
    protNormCoreQuiet(withDIAStateScope(manager, function() {
        v1 <- computeDIAReferenceProteinNormalisation(manager, getDIAStateReference(manager), list(norm_method = "none"))
        setDIAStepProvenance(manager, "protein_between_sample_normalisation", "whole", "norm_method", NULL,
            values = list(norm_method = "none"))
        finishDIAReferenceProteinStep(workflow, v1, publish = FALSE)
        v2 <- computeDIAReferenceProteinPeptideCountFilter(manager, getDIAStateReference(manager),
            list(num_peptides_per_protein_thresh = 1, num_peptidoforms_per_protein_thresh = 2))
        setDIAStepProvenance(manager, "protein_peptide_count_filter", "whole",
            c("num_peptides_per_protein_thresh", "num_peptidoforms_per_protein_thresh"), NULL, values = v2$config[1:2])
        finishDIAReferenceProteinStep(workflow, v2, publish = FALSE)
    }))
    list(root = root, paths = paths, workflow = workflow, manager = manager)
}

#' The tab's clicks on one owner, each captured as protNormW1Capture describes. The callbacks record the science
#' (12 characters of its digest) of the reference each receives.
protNormW1Clicks <- function(fixture) {
    norm <- new.env(parent = emptyenv())
    seen <- character()
    manager <- fixture$manager
    callback <- function(kind) function(reference) {
        seen <<- c(seen, paste0(kind, ":", substr(protNormW1Sha(protNormW1State(manager, reference)), 1L, 12L)))
        invisible(NULL)
    }
    capture <- function(run) {
        out <- protNormW1Capture(fixture, norm, run, seen)
        seen <<- character()
        out
    }
    aesthetics <- function() list(color_var = "group", shape_var = NULL)
    list(
        norm = norm,
        normalization = function(mode, ...) {
            settings <- utils::modifyList(protNormCoreSettings(mode), list(...))
            capture(protNormW1Run(runDIAReferenceNormalization(settings, fixture$workflow, norm, fixture$paths,
                "proteomics", "core", callback("post"), callback("ruv"), function() "group")))
        },
        correlation = function(threshold, skipped = FALSE) {
            capture(protNormW1Run(runDIAReferenceCorrelation(fixture$workflow, norm, fixture$paths, threshold, skipped,
                "proteomics", "core", if (skipped) function() NULL else function() "group", aesthetics)))
        },
        pass_through = function() {
            capture(protNormW1Run(runDIAReferenceLimpaNormalisationPassThrough(protNormCoreSettings("covariates"),
                fixture$workflow, norm, callback("post"))))
        },
        ruv_covariates = function(mode) {
            capture(protNormW1Run(runDIAReferenceRuvCovariates(list(ruv_mode = mode, color_variable = "group",
                shape_variable = NULL), fixture$workflow, norm, fixture$paths, "proteomics", "core", callback("ruv"),
                function() "group")))
        })
}

# The frozen values: two dumps on the start commit in separate processes, byte-identical (cmp), pasted from the first.
protNormW1Frozen <- function() list(
    dia_normalization_skip = list(
        current = "normalized",
        history = c("initial", "protein_replicate_filtered", "normalized"),
        head_current = "normalized",
        head_history = c("initial", "protein_replicate_filtered", "normalized"),
        states = c(protein_replicate_filtered = "932a9ccbf7a5c866708b9960912487faa3fa398571c005986822e936be65ca78", normalized = "87f9eba78ad3cf1d3961f3c9dbaa34f6fb1880e450c5e7018abe2c6fdfd68692"
            ),
        outputs = list(protein_replicate_filtered = NULL, normalized = c(normalization = TRUE)),
        products = "dia-normalization-ruv:fc952b8790c07658",
        provenance = "241946e9b6426a565a2d1574529eeece1ab4882f18afa37b820bc8562f25b6cb",
        record = "dad690818a6ce5c56ab230e01c35799f994d0c8e63948665c94c38593a10acea",
        files = "10f5f602bae9c8e28239b40f613dfa3179733b587f2c1cb5ce7299a8f2c9e868",
        norm_state = "0809b8018bafde2e009b93e4b479b9716af885734a5e949f687648f1663d46db",
        norm_fields = c(QC_composite_figure = "71557d1c8bac9bbe", best_k = "3fe0774c75a5f7ec", control_genes_index = "3fe0774c75a5f7ec", 
            filtering_summary_text = "71557d1c8bac9bbe", normalization_complete = "21dd4c7bf548536a", normalization_input_reference = "0c78c0cd259e526a", 
            normalized_protein_obj = "3ea964fe2023dddb", post_norm_filtering_plot = "71557d1c8bac9bbe", publication_paths = "71557d1c8bac9bbe", 
            qc_plot_paths = "28df520b992b4a4b", ruv_complete = "21dd4c7bf548536a", ruv_normalized_obj = "3ea964fe2023dddb", 
            ruv_optimization_result = "21d2cfada7829073", ruv_optimization_table = "71557d1c8bac9bbe"),
        ruv_result = "21d2cfada7829073d467cb158e5062c42c6507f5639f35da2b03362d6efb4581",
        workflow = "fbea021eb73497238620c200d02c93de5ed4ab1caa7a29be69dd1f7604357566",
        value = "21dd4c7bf548536a0a6a9811519c53a3223f1a4f0357aa5382dd8a49aad258b6",
        error = NA_character_,
        calls = c(assertDIAExportReady = 1L, runSavedOperation = 1L, commitDIAReferenceState = 1L, saveDIADisplayTable = 1L, 
            publishDIAActionOutputs = 1L, readDIAProteinProgress = 0L, updateProteinFiltering = 0L, writeDIADesignTable = 1L, 
            resolveDIAReferenceQcPlot = 0L, resolveProtNormCorrelationResultFilenames = 0L, generateProtNormCompositeQcFigure = 1L, 
            finalizeProtNormWorkflowState = 1L),
        callbacks = "post:32d9093cb75d"),
    dia_normalization_skip_repeat = list(
        current = "normalized",
        history = c("initial", "protein_replicate_filtered", "normalized"),
        head_current = "normalized",
        head_history = c("initial", "protein_replicate_filtered", "normalized"),
        states = c(protein_replicate_filtered = "932a9ccbf7a5c866708b9960912487faa3fa398571c005986822e936be65ca78", normalized = "87f9eba78ad3cf1d3961f3c9dbaa34f6fb1880e450c5e7018abe2c6fdfd68692"
            ),
        outputs = list(protein_replicate_filtered = NULL, normalized = c(normalization = TRUE)),
        products = "dia-normalization-ruv:fc952b8790c07658",
        provenance = "241946e9b6426a565a2d1574529eeece1ab4882f18afa37b820bc8562f25b6cb",
        record = "dad690818a6ce5c56ab230e01c35799f994d0c8e63948665c94c38593a10acea",
        files = "10f5f602bae9c8e28239b40f613dfa3179733b587f2c1cb5ce7299a8f2c9e868",
        norm_state = "0809b8018bafde2e009b93e4b479b9716af885734a5e949f687648f1663d46db",
        norm_fields = c(QC_composite_figure = "71557d1c8bac9bbe", best_k = "3fe0774c75a5f7ec", control_genes_index = "3fe0774c75a5f7ec", 
            filtering_summary_text = "71557d1c8bac9bbe", normalization_complete = "21dd4c7bf548536a", normalization_input_reference = "0c78c0cd259e526a", 
            normalized_protein_obj = "3ea964fe2023dddb", post_norm_filtering_plot = "71557d1c8bac9bbe", publication_paths = "71557d1c8bac9bbe", 
            qc_plot_paths = "28df520b992b4a4b", ruv_complete = "21dd4c7bf548536a", ruv_normalized_obj = "3ea964fe2023dddb", 
            ruv_optimization_result = "21d2cfada7829073", ruv_optimization_table = "71557d1c8bac9bbe"),
        ruv_result = "21d2cfada7829073d467cb158e5062c42c6507f5639f35da2b03362d6efb4581",
        workflow = "fbea021eb73497238620c200d02c93de5ed4ab1caa7a29be69dd1f7604357566",
        value = "21dd4c7bf548536a0a6a9811519c53a3223f1a4f0357aa5382dd8a49aad258b6",
        error = NA_character_,
        calls = c(assertDIAExportReady = 1L, runSavedOperation = 1L, commitDIAReferenceState = 0L, saveDIADisplayTable = 1L, 
            publishDIAActionOutputs = 1L, readDIAProteinProgress = 0L, updateProteinFiltering = 0L, writeDIADesignTable = 0L, 
            resolveDIAReferenceQcPlot = 0L, resolveProtNormCorrelationResultFilenames = 0L, generateProtNormCompositeQcFigure = 0L, 
            finalizeProtNormWorkflowState = 1L),
        callbacks = character(0)),
    dia_normalization_manual = list(
        current = "ruv_corrected",
        history = c("initial", "protein_replicate_filtered", "normalized", "ruv_corrected"),
        head_current = "ruv_corrected",
        head_history = c("initial", "protein_replicate_filtered", "normalized", "ruv_corrected"),
        states = c(protein_replicate_filtered = "932a9ccbf7a5c866708b9960912487faa3fa398571c005986822e936be65ca78", normalized = "87f9eba78ad3cf1d3961f3c9dbaa34f6fb1880e450c5e7018abe2c6fdfd68692", 
            ruv_corrected = "89eaade153d179c743701f10610a68323d35b21ce15c8a62b2391b89cc55dbe8"),
        outputs = list(protein_replicate_filtered = NULL, normalized = c(normalization = TRUE), ruv_corrected = c(normalization = TRUE)),
        products = c("dia-normalization-ruv:a2ff7f83981cff63", "dia-normalization-ruv:fc952b8790c07658", "output/dia-protein-filtering-counts:2b3a70b7fdc2233a"
            ),
        provenance = "4916dc4d938cd3913c9c10d9e01cac46959630236d768aa52eca9ff151d35826",
        record = "d19acc85da5ba7999f4af330499f1c93d518fa9cd94ce02091286207b8c9867e",
        files = "25c82c19abbd083fc8dfea532ce38b6f0a80eaadc861d6c905ae4846b2695fbc",
        norm_state = "8892d37afe006d0afc191cdd8bc3a4168ff6e471d1a688aad331501ca4bae334",
        norm_fields = c(QC_composite_figure = "71557d1c8bac9bbe", best_k = "60235a10c1ff424d", control_genes_index = "f1eb628a13bb1aeb", 
            filtering_summary_text = "f5162d65b6daff8e", normalization_complete = "21dd4c7bf548536a", normalization_input_reference = "0c78c0cd259e526a", 
            normalized_protein_obj = "3ea964fe2023dddb", post_norm_filtering_plot = "1cc1786a46386211", publication_paths = "71557d1c8bac9bbe", 
            qc_plot_paths = "28df520b992b4a4b", ruv_complete = "21dd4c7bf548536a", ruv_normalized_obj = "8b5b1c46fa8314c7", 
            ruv_optimization_result = "ffc12bf635e5c112", ruv_optimization_table = "c7841b609b82546d"),
        ruv_result = "330ecfe26fe03dbb4a4374188c4a208ecb7e13b46ed8772d2884d547a88b28ad",
        workflow = "458611adfe717bdfaabe100a28e62d4a729d7974a0ebfca2ed4643e561ede894",
        value = "21dd4c7bf548536a0a6a9811519c53a3223f1a4f0357aa5382dd8a49aad258b6",
        error = NA_character_,
        calls = c(assertDIAExportReady = 1L, runSavedOperation = 1L, commitDIAReferenceState = 1L, saveDIADisplayTable = 1L, 
            publishDIAActionOutputs = 1L, readDIAProteinProgress = 1L, updateProteinFiltering = 1L, writeDIADesignTable = 1L, 
            resolveDIAReferenceQcPlot = 0L, resolveProtNormCorrelationResultFilenames = 0L, generateProtNormCompositeQcFigure = 1L, 
            finalizeProtNormWorkflowState = 1L),
        callbacks = c("post:32d9093cb75d", "ruv:a50c91a55ed5")),
    dia_correlation_apply = list(
        current = "correlation_filtered",
        history = c("initial", "protein_replicate_filtered", "normalized", "ruv_corrected", "correlation_filtered"),
        head_current = "correlation_filtered",
        head_history = c("initial", "protein_replicate_filtered", "normalized", "ruv_corrected", "correlation_filtered"),
        states = c(protein_replicate_filtered = "932a9ccbf7a5c866708b9960912487faa3fa398571c005986822e936be65ca78", normalized = "87f9eba78ad3cf1d3961f3c9dbaa34f6fb1880e450c5e7018abe2c6fdfd68692", 
            ruv_corrected = "89eaade153d179c743701f10610a68323d35b21ce15c8a62b2391b89cc55dbe8", correlation_filtered = "91a967cfb1269399f1f49928f51211e97f373a82715700d7c0018c8f604823a8"
            ),
        outputs = list(protein_replicate_filtered = NULL, normalized = c(normalization = TRUE), ruv_corrected = c(normalization = TRUE), 
                correlation_filtered = c(correlation = TRUE)),
        products = c("dia-normalization-ruv:a2ff7f83981cff63", "dia-normalization-ruv:fc952b8790c07658", "dia-pearson-pairs:e4d56e08c2383db2", 
            "dia-pearson-selection:ac622a900a2062ab", "dia_normalization_pca:a15ef8c115aaba13", "output/dia-protein-filtering-counts:2b3a70b7fdc2233a"
            ),
        provenance = "68f3678ffee19ddb999f986041ed409dca3a439f0ddb7c8bb2ea8e3bb415e7f0",
        record = "df134a21406aa67bb9db7a625937f8248f3f9b8e52f9a30cec91f6f95363a7d8",
        files = "5b238a84991d085e2867dcc8a02dbde4e206504d440a971e42fa22662e741d23",
        norm_state = "fb69ce836f260cd45d6b3c5c219cabaca658005228a0634cab95720c76acfec2",
        norm_fields = c(QC_composite_figure = "71557d1c8bac9bbe", best_k = "60235a10c1ff424d", control_genes_index = "f1eb628a13bb1aeb", 
            correlation_filtered_obj = "8b5b1c46fa8314c7", correlation_threshold = "acf06dd9ce0cc864", correlation_vector = "e4d56e08c2383db2", 
            filtering_summary_text = "f5162d65b6daff8e", final_filtering_plot = "1cc1786a46386211", final_qc_plot = "0833312183929473", 
            normalization_complete = "21dd4c7bf548536a", normalization_input_reference = "0c78c0cd259e526a", normalized_protein_obj = "3ea964fe2023dddb", 
            post_norm_filtering_plot = "1cc1786a46386211", publication_paths = "71557d1c8bac9bbe", qc_plot_paths = "28df520b992b4a4b", 
            ruv_complete = "21dd4c7bf548536a", ruv_normalized_obj = "8b5b1c46fa8314c7", ruv_optimization_result = "ffc12bf635e5c112", 
            ruv_optimization_table = "c7841b609b82546d"),
        ruv_result = "330ecfe26fe03dbb4a4374188c4a208ecb7e13b46ed8772d2884d547a88b28ad",
        workflow = "458611adfe717bdfaabe100a28e62d4a729d7974a0ebfca2ed4643e561ede894",
        value = "8b5b1c46fa8314c7646cf463787c44ff6ee60255f50cadbc5f5f26f5b9290034",
        error = NA_character_,
        calls = c(assertDIAExportReady = 1L, runSavedOperation = 2L, commitDIAReferenceState = 1L, saveDIADisplayTable = 0L, 
            publishDIAActionOutputs = 1L, readDIAProteinProgress = 1L, updateProteinFiltering = 1L, writeDIADesignTable = 1L, 
            resolveDIAReferenceQcPlot = 1L, resolveProtNormCorrelationResultFilenames = 1L, generateProtNormCompositeQcFigure = 0L, 
            finalizeProtNormWorkflowState = 0L),
        callbacks = character(0)),
    dia_correlation_reapply = list(
        current = "correlation_filtered",
        history = c("initial", "protein_replicate_filtered", "normalized", "ruv_corrected", "correlation_filtered"),
        head_current = "correlation_filtered",
        head_history = c("initial", "protein_replicate_filtered", "normalized", "ruv_corrected", "correlation_filtered"),
        states = c(protein_replicate_filtered = "932a9ccbf7a5c866708b9960912487faa3fa398571c005986822e936be65ca78", normalized = "87f9eba78ad3cf1d3961f3c9dbaa34f6fb1880e450c5e7018abe2c6fdfd68692", 
            ruv_corrected = "89eaade153d179c743701f10610a68323d35b21ce15c8a62b2391b89cc55dbe8", correlation_filtered = "416d38a77b1e8bb581ad6f46b7101c18d653bc2691a5d139976e143dab4ebc9e"
            ),
        outputs = list(protein_replicate_filtered = NULL, normalized = c(normalization = TRUE), ruv_corrected = c(normalization = TRUE), 
                correlation_filtered = c(correlation = TRUE)),
        products = c("dia-normalization-ruv:a2ff7f83981cff63", "dia-normalization-ruv:fc952b8790c07658", "dia-pearson-pairs:e4d56e08c2383db2", 
            "dia-pearson-selection:ac622a900a2062ab", "dia-pearson-selection:ac622a900a2062ab", "dia_normalization_pca:a15ef8c115aaba13", 
            "dia_normalization_pca:a15ef8c115aaba13", "output/dia-protein-filtering-counts:2b3a70b7fdc2233a"),
        provenance = "a181cb9aa9e95f4f52221586dfc19c10d6e3197d780ee3a1d7fc4f559fa2e579",
        record = "4c68dd10bf802b0822f4a5fc48adcbbb2010e98162a4e9f018005174d820e834",
        files = "5b238a84991d085e2867dcc8a02dbde4e206504d440a971e42fa22662e741d23",
        norm_state = "fb69ce836f260cd45d6b3c5c219cabaca658005228a0634cab95720c76acfec2",
        norm_fields = c(QC_composite_figure = "71557d1c8bac9bbe", best_k = "60235a10c1ff424d", control_genes_index = "f1eb628a13bb1aeb", 
            correlation_filtered_obj = "8b5b1c46fa8314c7", correlation_threshold = "29a39c60c37c4f96", correlation_vector = "e4d56e08c2383db2", 
            filtering_summary_text = "f5162d65b6daff8e", final_filtering_plot = "1cc1786a46386211", final_qc_plot = "0833312183929473", 
            normalization_complete = "21dd4c7bf548536a", normalization_input_reference = "0c78c0cd259e526a", normalized_protein_obj = "3ea964fe2023dddb", 
            post_norm_filtering_plot = "1cc1786a46386211", publication_paths = "71557d1c8bac9bbe", qc_plot_paths = "28df520b992b4a4b", 
            ruv_complete = "21dd4c7bf548536a", ruv_normalized_obj = "8b5b1c46fa8314c7", ruv_optimization_result = "ffc12bf635e5c112", 
            ruv_optimization_table = "c7841b609b82546d"),
        ruv_result = "330ecfe26fe03dbb4a4374188c4a208ecb7e13b46ed8772d2884d547a88b28ad",
        workflow = "458611adfe717bdfaabe100a28e62d4a729d7974a0ebfca2ed4643e561ede894",
        value = "8b5b1c46fa8314c7646cf463787c44ff6ee60255f50cadbc5f5f26f5b9290034",
        error = NA_character_,
        calls = c(assertDIAExportReady = 1L, runSavedOperation = 2L, commitDIAReferenceState = 1L, saveDIADisplayTable = 0L, 
            publishDIAActionOutputs = 1L, readDIAProteinProgress = 1L, updateProteinFiltering = 1L, writeDIADesignTable = 1L, 
            resolveDIAReferenceQcPlot = 1L, resolveProtNormCorrelationResultFilenames = 1L, generateProtNormCompositeQcFigure = 0L, 
            finalizeProtNormWorkflowState = 0L),
        callbacks = character(0)),
    dia_correlation_skip = list(
        current = "correlation_filtered",
        history = c("initial", "protein_replicate_filtered", "normalized", "ruv_corrected", "correlation_filtered"),
        head_current = "correlation_filtered",
        head_history = c("initial", "protein_replicate_filtered", "normalized", "ruv_corrected", "correlation_filtered"),
        states = c(protein_replicate_filtered = "932a9ccbf7a5c866708b9960912487faa3fa398571c005986822e936be65ca78", normalized = "87f9eba78ad3cf1d3961f3c9dbaa34f6fb1880e450c5e7018abe2c6fdfd68692", 
            ruv_corrected = "89eaade153d179c743701f10610a68323d35b21ce15c8a62b2391b89cc55dbe8", correlation_filtered = "d36fa4aa2a6fabdae4f4cd9690cb1b8c74f363d72e834ba7c434614b1fea0c4d"
            ),
        outputs = list(protein_replicate_filtered = NULL, normalized = c(normalization = TRUE), ruv_corrected = c(normalization = TRUE), 
                correlation_filtered = c(correlation = TRUE)),
        products = c("dia-normalization-ruv:a2ff7f83981cff63", "dia-normalization-ruv:fc952b8790c07658", "dia-pearson-pairs:e4d56e08c2383db2", 
            "dia-pearson-selection:ac622a900a2062ab", "dia-pearson-selection:ac622a900a2062ab", "dia_normalization_pca:a15ef8c115aaba13", 
            "dia_normalization_pca:a15ef8c115aaba13", "dia_normalization_pca:a15ef8c115aaba13", "output/dia-protein-filtering-counts:2b3a70b7fdc2233a"
            ),
        provenance = "a181cb9aa9e95f4f52221586dfc19c10d6e3197d780ee3a1d7fc4f559fa2e579",
        record = "aad6d183e1fb0e4cded898b0b032d6725b09fab5d64ca72911faef202f0b63d3",
        files = "5b238a84991d085e2867dcc8a02dbde4e206504d440a971e42fa22662e741d23",
        norm_state = "fb69ce836f260cd45d6b3c5c219cabaca658005228a0634cab95720c76acfec2",
        norm_fields = c(QC_composite_figure = "71557d1c8bac9bbe", best_k = "60235a10c1ff424d", control_genes_index = "f1eb628a13bb1aeb", 
            correlation_filtered_obj = "8b5b1c46fa8314c7", correlation_threshold = "a991209291aaafd6", correlation_vector = "e4d56e08c2383db2", 
            filtering_summary_text = "f5162d65b6daff8e", final_filtering_plot = "1cc1786a46386211", final_qc_plot = "b9f17e080d57b5bc", 
            normalization_complete = "21dd4c7bf548536a", normalization_input_reference = "0c78c0cd259e526a", normalized_protein_obj = "3ea964fe2023dddb", 
            post_norm_filtering_plot = "1cc1786a46386211", publication_paths = "71557d1c8bac9bbe", qc_plot_paths = "28df520b992b4a4b", 
            ruv_complete = "21dd4c7bf548536a", ruv_normalized_obj = "8b5b1c46fa8314c7", ruv_optimization_result = "ffc12bf635e5c112", 
            ruv_optimization_table = "c7841b609b82546d"),
        ruv_result = "330ecfe26fe03dbb4a4374188c4a208ecb7e13b46ed8772d2884d547a88b28ad",
        workflow = "458611adfe717bdfaabe100a28e62d4a729d7974a0ebfca2ed4643e561ede894",
        value = "8b5b1c46fa8314c7646cf463787c44ff6ee60255f50cadbc5f5f26f5b9290034",
        error = NA_character_,
        calls = c(assertDIAExportReady = 1L, runSavedOperation = 0L, commitDIAReferenceState = 1L, saveDIADisplayTable = 0L, 
            publishDIAActionOutputs = 1L, readDIAProteinProgress = 1L, updateProteinFiltering = 1L, writeDIADesignTable = 1L, 
            resolveDIAReferenceQcPlot = 1L, resolveProtNormCorrelationResultFilenames = 1L, generateProtNormCompositeQcFigure = 0L, 
            finalizeProtNormWorkflowState = 0L),
        callbacks = character(0)),
    dia_normalization_refused = list(
        current = "correlation_filtered",
        history = c("initial", "protein_replicate_filtered", "normalized", "ruv_corrected", "correlation_filtered"),
        head_current = "correlation_filtered",
        head_history = c("initial", "protein_replicate_filtered", "normalized", "ruv_corrected", "correlation_filtered"),
        states = c(protein_replicate_filtered = "932a9ccbf7a5c866708b9960912487faa3fa398571c005986822e936be65ca78", normalized = "87f9eba78ad3cf1d3961f3c9dbaa34f6fb1880e450c5e7018abe2c6fdfd68692", 
            ruv_corrected = "89eaade153d179c743701f10610a68323d35b21ce15c8a62b2391b89cc55dbe8", correlation_filtered = "d36fa4aa2a6fabdae4f4cd9690cb1b8c74f363d72e834ba7c434614b1fea0c4d"
            ),
        outputs = list(protein_replicate_filtered = NULL, normalized = c(normalization = TRUE), ruv_corrected = c(normalization = TRUE), 
                correlation_filtered = c(correlation = TRUE)),
        products = c("dia-normalization-ruv:a2ff7f83981cff63", "dia-normalization-ruv:fc952b8790c07658", "dia-pearson-pairs:e4d56e08c2383db2", 
            "dia-pearson-selection:ac622a900a2062ab", "dia-pearson-selection:ac622a900a2062ab", "dia_normalization_pca:a15ef8c115aaba13", 
            "dia_normalization_pca:a15ef8c115aaba13", "dia_normalization_pca:a15ef8c115aaba13", "output/dia-protein-filtering-counts:2b3a70b7fdc2233a"
            ),
        provenance = "a181cb9aa9e95f4f52221586dfc19c10d6e3197d780ee3a1d7fc4f559fa2e579",
        record = "aad6d183e1fb0e4cded898b0b032d6725b09fab5d64ca72911faef202f0b63d3",
        files = "5b238a84991d085e2867dcc8a02dbde4e206504d440a971e42fa22662e741d23",
        norm_state = "fb69ce836f260cd45d6b3c5c219cabaca658005228a0634cab95720c76acfec2",
        norm_fields = c(QC_composite_figure = "71557d1c8bac9bbe", best_k = "60235a10c1ff424d", control_genes_index = "f1eb628a13bb1aeb", 
            correlation_filtered_obj = "8b5b1c46fa8314c7", correlation_threshold = "a991209291aaafd6", correlation_vector = "e4d56e08c2383db2", 
            filtering_summary_text = "f5162d65b6daff8e", final_filtering_plot = "1cc1786a46386211", final_qc_plot = "b9f17e080d57b5bc", 
            normalization_complete = "21dd4c7bf548536a", normalization_input_reference = "0c78c0cd259e526a", normalized_protein_obj = "3ea964fe2023dddb", 
            post_norm_filtering_plot = "1cc1786a46386211", publication_paths = "71557d1c8bac9bbe", qc_plot_paths = "28df520b992b4a4b", 
            ruv_complete = "21dd4c7bf548536a", ruv_normalized_obj = "8b5b1c46fa8314c7", ruv_optimization_result = "ffc12bf635e5c112", 
            ruv_optimization_table = "c7841b609b82546d"),
        ruv_result = "330ecfe26fe03dbb4a4374188c4a208ecb7e13b46ed8772d2884d547a88b28ad",
        workflow = "458611adfe717bdfaabe100a28e62d4a729d7974a0ebfca2ed4643e561ede894",
        value = "71557d1c8bac9bbe3cbec8d00bb223a2f372279827064095447e569fbf5a760a",
        error = "RUV-III cannot estimate 50 factors from 32 control proteins observed in every sample\ni choose a smaller k or a larger control set, filter for completeness first, or run normalisation without RUV",
        calls = c(assertDIAExportReady = 1L, runSavedOperation = 1L, commitDIAReferenceState = 0L, saveDIADisplayTable = 0L, 
            publishDIAActionOutputs = 0L, readDIAProteinProgress = 0L, updateProteinFiltering = 0L, writeDIADesignTable = 0L, 
            resolveDIAReferenceQcPlot = 0L, resolveProtNormCorrelationResultFilenames = 0L, generateProtNormCompositeQcFigure = 0L, 
            finalizeProtNormWorkflowState = 0L),
        callbacks = character(0)),
    dia_normalization_retry = list(
        current = "ruv_corrected",
        history = c("initial", "protein_replicate_filtered", "normalized", "ruv_corrected", "correlation_filtered"),
        head_current = "ruv_corrected",
        head_history = c("initial", "protein_replicate_filtered", "normalized", "ruv_corrected", "correlation_filtered"),
        states = c(protein_replicate_filtered = "932a9ccbf7a5c866708b9960912487faa3fa398571c005986822e936be65ca78", normalized = "87f9eba78ad3cf1d3961f3c9dbaa34f6fb1880e450c5e7018abe2c6fdfd68692", 
            ruv_corrected = "2c7fbb38bdb82f2d9f78a3296aea6df1afe0e35fa5fa1127041ec805df74cad7", correlation_filtered = "d36fa4aa2a6fabdae4f4cd9690cb1b8c74f363d72e834ba7c434614b1fea0c4d"
            ),
        outputs = list(protein_replicate_filtered = NULL, normalized = c(normalization = TRUE), ruv_corrected = c(normalization = TRUE), 
                correlation_filtered = c(correlation = TRUE)),
        products = c("dia-normalization-ruv:a2ff7f83981cff63", "dia-normalization-ruv:fc952b8790c07658", "dia-pearson-pairs:e4d56e08c2383db2", 
            "dia-pearson-selection:ac622a900a2062ab", "dia-pearson-selection:ac622a900a2062ab", "dia_normalization_pca:a15ef8c115aaba13", 
            "dia_normalization_pca:a15ef8c115aaba13", "dia_normalization_pca:a15ef8c115aaba13", "output/dia-protein-filtering-counts:2b3a70b7fdc2233a"
            ),
        provenance = "a181cb9aa9e95f4f52221586dfc19c10d6e3197d780ee3a1d7fc4f559fa2e579",
        record = "aad6d183e1fb0e4cded898b0b032d6725b09fab5d64ca72911faef202f0b63d3",
        files = "5b238a84991d085e2867dcc8a02dbde4e206504d440a971e42fa22662e741d23",
        norm_state = "0a885333fe431e018f77a3073f63740973a79c88b15b3807dfa0ee58bfef2fe6",
        norm_fields = c(QC_composite_figure = "71557d1c8bac9bbe", best_k = "60235a10c1ff424d", control_genes_index = "f1eb628a13bb1aeb", 
            correlation_filtered_obj = "8b5b1c46fa8314c7", correlation_threshold = "a991209291aaafd6", correlation_vector = "e4d56e08c2383db2", 
            filtering_summary_text = "f5162d65b6daff8e", final_filtering_plot = "1cc1786a46386211", final_qc_plot = "b9f17e080d57b5bc", 
            normalization_complete = "21dd4c7bf548536a", normalization_input_reference = "0c78c0cd259e526a", normalized_protein_obj = "3ea964fe2023dddb", 
            post_norm_filtering_plot = "1cc1786a46386211", publication_paths = "71557d1c8bac9bbe", qc_plot_paths = "28df520b992b4a4b", 
            ruv_complete = "21dd4c7bf548536a", ruv_normalized_obj = "8b5b1c46fa8314c7", ruv_optimization_result = "ffc12bf635e5c112", 
            ruv_optimization_table = "c7841b609b82546d"),
        ruv_result = "330ecfe26fe03dbb4a4374188c4a208ecb7e13b46ed8772d2884d547a88b28ad",
        workflow = "458611adfe717bdfaabe100a28e62d4a729d7974a0ebfca2ed4643e561ede894",
        value = "21dd4c7bf548536a0a6a9811519c53a3223f1a4f0357aa5382dd8a49aad258b6",
        error = NA_character_,
        calls = c(assertDIAExportReady = 1L, runSavedOperation = 1L, commitDIAReferenceState = 1L, saveDIADisplayTable = 1L, 
            publishDIAActionOutputs = 1L, readDIAProteinProgress = 1L, updateProteinFiltering = 1L, writeDIADesignTable = 1L, 
            resolveDIAReferenceQcPlot = 0L, resolveProtNormCorrelationResultFilenames = 0L, generateProtNormCompositeQcFigure = 1L, 
            finalizeProtNormWorkflowState = 1L),
        callbacks = c("post:32d9093cb75d", "ruv:a50c91a55ed5")),
    limpa_covariates_pass_through = list(
        current = "peptide_count_filtered",
        history = c("initial", "protein_dpc_quant", "protein_normalised", "peptide_count_filtered"),
        head_current = "peptide_count_filtered",
        head_history = c("initial", "protein_dpc_quant", "protein_normalised", "peptide_count_filtered"),
        states = c(protein_dpc_quant = "41dfd92a6007ea8f051ed9a447d3688f61902363a2b9aa84aedf596e4c0e9549", protein_normalised = "ab8d0aac28c16422de25afd8b64efc63a80bdf5670f34d7bfe526bc0d88edd64", 
            peptide_count_filtered = "75c25ed13486d80bb56c8e3691369dde2d61de37bcacecf3224cdc3573723fd0"),
        outputs = list(protein_dpc_quant = NULL, protein_normalised = NULL, peptide_count_filtered = NULL),
        products = c("output/dia-protein-filtering-counts:67049a631f7b624a", "output/dia-protein-qc-output:8b8101851eb54168", 
            "output/dia-protein-qc-output:deef48a7ebdbcb3a"),
        provenance = "b735d98157cc15fc39d183d25e744262d612b0f19d994ce5b073da963f5b60be",
        record = "6d3d46d9c054af2bb3b2f07469d89981e294f90d4f84a5d713938372b89e86d2",
        files = "8ab794fecc83e91d92f2fee89cd8fdc808bd44a3e8c693472e90be9cbf830156",
        norm_state = "d61f7e02b4c6a6d12fca78b9106ca07e9c39c3d7f406541c74a2ddc2947c3167",
        norm_fields = c(best_k = "3fe0774c75a5f7ec", control_genes_index = "3fe0774c75a5f7ec", filtering_summary_text = "1ba7030e9325c39d", 
            normalization_complete = "21dd4c7bf548536a", normalization_input_reference = "7e1a356bda20ce87", normalized_protein_obj = "7e1a356bda20ce87", 
            ruv_normalized_obj = "7e1a356bda20ce87", ruv_optimization_result = "71557d1c8bac9bbe"),
        ruv_result = "71557d1c8bac9bbe3cbec8d00bb223a2f372279827064095447e569fbf5a760a",
        workflow = "31f37cefd94f75f4f1e164b03125e1825178d30df82fc99a89aa16b7ba947e3a",
        value = "21dd4c7bf548536a0a6a9811519c53a3223f1a4f0357aa5382dd8a49aad258b6",
        error = NA_character_,
        calls = c(assertDIAExportReady = 1L, runSavedOperation = 0L, commitDIAReferenceState = 0L, saveDIADisplayTable = 0L, 
            publishDIAActionOutputs = 0L, readDIAProteinProgress = 0L, updateProteinFiltering = 0L, writeDIADesignTable = 0L, 
            resolveDIAReferenceQcPlot = 0L, resolveProtNormCorrelationResultFilenames = 0L, generateProtNormCompositeQcFigure = 0L, 
            finalizeProtNormWorkflowState = 0L),
        callbacks = "post:f08bb05015f2"),
    limpa_covariates_correlation = list(
        current = "correlation_filtered",
        history = c("initial", "protein_dpc_quant", "protein_normalised", "peptide_count_filtered", "correlation_filtered"
            ),
        head_current = "correlation_filtered",
        head_history = c("initial", "protein_dpc_quant", "protein_normalised", "peptide_count_filtered", "correlation_filtered"
            ),
        states = c(protein_dpc_quant = "41dfd92a6007ea8f051ed9a447d3688f61902363a2b9aa84aedf596e4c0e9549", protein_normalised = "ab8d0aac28c16422de25afd8b64efc63a80bdf5670f34d7bfe526bc0d88edd64", 
            peptide_count_filtered = "75c25ed13486d80bb56c8e3691369dde2d61de37bcacecf3224cdc3573723fd0", correlation_filtered = "312ee5d40ac0d3f3c809c3864a1d440e1dbc000c8cc75a98cfdd9c95b96f5935"
            ),
        outputs = list(protein_dpc_quant = NULL, protein_normalised = NULL, peptide_count_filtered = NULL, correlation_filtered = c(correlation = TRUE)),
        products = c("dia-pearson-pairs:9eaece3602f93d46", "dia-pearson-selection:36ea8f5f2e373570", "dia_normalization_pca:d9b5720af8d8a298", 
            "output/dia-protein-filtering-counts:67049a631f7b624a", "output/dia-protein-qc-output:8b8101851eb54168", 
            "output/dia-protein-qc-output:deef48a7ebdbcb3a"),
        provenance = "9201164d5d72b3fe7f7ee72607a2e03edc64d9eb2a7d6c94b87505b9bd14cc62",
        record = "181955afdf5b3fc2465446465437e447f422ace7d0f85293e0d5a1e0fb446ff3",
        files = "34e64dc1e1e246298b07d49ed56b42a6743ca9e662a09925d3be94dc69c4735c",
        norm_state = "3b0fbc63724b4e8bf08fb33060738e1459d38205eb8bd692f03e358d03181fa1",
        norm_fields = c(best_k = "3fe0774c75a5f7ec", control_genes_index = "3fe0774c75a5f7ec", correlation_filtered_obj = "7e1a356bda20ce87", 
            correlation_threshold = "acf06dd9ce0cc864", correlation_vector = "9eaece3602f93d46", filtering_summary_text = "1ba7030e9325c39d", 
            final_filtering_plot = "1cc1786a46386211", final_qc_plot = "1252e1714ac6c9a2", normalization_complete = "21dd4c7bf548536a", 
            normalization_input_reference = "7e1a356bda20ce87", normalized_protein_obj = "7e1a356bda20ce87", ruv_normalized_obj = "7e1a356bda20ce87", 
            ruv_optimization_result = "71557d1c8bac9bbe"),
        ruv_result = "71557d1c8bac9bbe3cbec8d00bb223a2f372279827064095447e569fbf5a760a",
        workflow = "31f37cefd94f75f4f1e164b03125e1825178d30df82fc99a89aa16b7ba947e3a",
        value = "7e1a356bda20ce879fd8f423db352ae28d5e6b171ffb73d4cabec81acfa36136",
        error = NA_character_,
        calls = c(assertDIAExportReady = 1L, runSavedOperation = 2L, commitDIAReferenceState = 1L, saveDIADisplayTable = 0L, 
            publishDIAActionOutputs = 1L, readDIAProteinProgress = 1L, updateProteinFiltering = 1L, writeDIADesignTable = 1L, 
            resolveDIAReferenceQcPlot = 1L, resolveProtNormCorrelationResultFilenames = 1L, generateProtNormCompositeQcFigure = 0L, 
            finalizeProtNormWorkflowState = 0L),
        callbacks = character(0)),
    limpa_covariates_ruv = list(
        current = "ruv_covariates",
        history = c("initial", "protein_dpc_quant", "protein_normalised", "peptide_count_filtered", "correlation_filtered", 
            "ruv_covariates"),
        head_current = "ruv_covariates",
        head_history = c("initial", "protein_dpc_quant", "protein_normalised", "peptide_count_filtered", "correlation_filtered", 
            "ruv_covariates"),
        states = c(protein_dpc_quant = "41dfd92a6007ea8f051ed9a447d3688f61902363a2b9aa84aedf596e4c0e9549", protein_normalised = "ab8d0aac28c16422de25afd8b64efc63a80bdf5670f34d7bfe526bc0d88edd64", 
            peptide_count_filtered = "75c25ed13486d80bb56c8e3691369dde2d61de37bcacecf3224cdc3573723fd0", correlation_filtered = "312ee5d40ac0d3f3c809c3864a1d440e1dbc000c8cc75a98cfdd9c95b96f5935", 
            ruv_covariates = "b2442669d93f3fbceeca3c61238b3fef1991c236c7ab4391e9f2dbe650ec072e"),
        outputs = list(protein_dpc_quant = NULL, protein_normalised = NULL, peptide_count_filtered = NULL, correlation_filtered = c(correlation = TRUE), 
                ruv_covariates = c("ruv-covariates" = TRUE)),
        products = c("dia-pearson-pairs:9eaece3602f93d46", "dia-pearson-selection:36ea8f5f2e373570", "dia-ruv-covariates:52b9cf17c41aa8c0", 
            "dia_normalization_pca:d9b5720af8d8a298", "output/dia-protein-filtering-counts:67049a631f7b624a", "output/dia-protein-qc-output:8b8101851eb54168", 
            "output/dia-protein-qc-output:deef48a7ebdbcb3a"),
        provenance = "b277eb8197ce4789839d312cf291fed25eae7f1c6b5450795b7b6706b9c33067",
        record = "4a271d4fbb4ee0107c841aa8738df37ceaa193138ffcbc27f137f8c354350b38",
        files = "b8e5c65e86cc0c32a8345221bfe52f4ecf30b8bdf26b4a213f5cd207ba3e4803",
        norm_state = "abc3041ac1452514dd38bea7dac56eaeb710139787ec11909029bb4ca7387977",
        norm_fields = c(best_k = "a7b8541c0f557ee0", control_genes_index = "1c497b478576b300", correlation_filtered_obj = "8f53dc8138ce299d", 
            correlation_threshold = "acf06dd9ce0cc864", correlation_vector = "9eaece3602f93d46", filtering_summary_text = "1ba7030e9325c39d", 
            final_filtering_plot = "1cc1786a46386211", final_qc_plot = "1252e1714ac6c9a2", normalization_complete = "21dd4c7bf548536a", 
            normalization_input_reference = "7e1a356bda20ce87", normalized_protein_obj = "7e1a356bda20ce87", publication_paths = "71557d1c8bac9bbe", 
            ruv_covariates_complete = "21dd4c7bf548536a", ruv_normalized_obj = "8f53dc8138ce299d", ruv_optimization_result = "eeae0389f229a0bd", 
            ruv_optimization_table = "bf419098d53c4e33"),
        ruv_result = "dd5a614c8a2494b1c66e75e08ce2a65e9259d557dcfa5a4058b9f720b5e8f4e6",
        workflow = "1aea9a4f60d498c1d539b6ced76fa6d4522eb4875db793ea53cf2a73fe2e049d",
        value = "8f53dc8138ce299d67fb2dbeb73332a20afa4f5e15c31dd7f3aefd09892e4f32",
        error = NA_character_,
        calls = c(assertDIAExportReady = 1L, runSavedOperation = 1L, commitDIAReferenceState = 1L, saveDIADisplayTable = 1L, 
            publishDIAActionOutputs = 1L, readDIAProteinProgress = 0L, updateProteinFiltering = 0L, writeDIADesignTable = 2L, 
            resolveDIAReferenceQcPlot = 0L, resolveProtNormCorrelationResultFilenames = 0L, generateProtNormCompositeQcFigure = 0L, 
            finalizeProtNormWorkflowState = 0L),
        callbacks = "ruv:e1c1cd4c99f5"),
    limpa_covariates_ruv_repeat = list(
        current = "ruv_covariates",
        history = c("initial", "protein_dpc_quant", "protein_normalised", "peptide_count_filtered", "correlation_filtered", 
            "ruv_covariates"),
        head_current = "ruv_covariates",
        head_history = c("initial", "protein_dpc_quant", "protein_normalised", "peptide_count_filtered", "correlation_filtered", 
            "ruv_covariates"),
        states = c(protein_dpc_quant = "41dfd92a6007ea8f051ed9a447d3688f61902363a2b9aa84aedf596e4c0e9549", protein_normalised = "ab8d0aac28c16422de25afd8b64efc63a80bdf5670f34d7bfe526bc0d88edd64", 
            peptide_count_filtered = "75c25ed13486d80bb56c8e3691369dde2d61de37bcacecf3224cdc3573723fd0", correlation_filtered = "312ee5d40ac0d3f3c809c3864a1d440e1dbc000c8cc75a98cfdd9c95b96f5935", 
            ruv_covariates = "49a8e114895c4afeb722b4f41bb56e96c18da6abc6351d3014c5f8c571b92722"),
        outputs = list(protein_dpc_quant = NULL, protein_normalised = NULL, peptide_count_filtered = NULL, correlation_filtered = c(correlation = TRUE), 
                ruv_covariates = c("ruv-covariates" = TRUE)),
        products = c("dia-pearson-pairs:9eaece3602f93d46", "dia-pearson-selection:36ea8f5f2e373570", "dia-ruv-covariates:52b9cf17c41aa8c0", 
            "dia-ruv-covariates:63c261200c1ced1e", "dia_normalization_pca:d9b5720af8d8a298", "output/dia-protein-filtering-counts:67049a631f7b624a", 
            "output/dia-protein-qc-output:8b8101851eb54168", "output/dia-protein-qc-output:deef48a7ebdbcb3a"),
        provenance = "b277eb8197ce4789839d312cf291fed25eae7f1c6b5450795b7b6706b9c33067",
        record = "4a271d4fbb4ee0107c841aa8738df37ceaa193138ffcbc27f137f8c354350b38",
        files = "b8e5c65e86cc0c32a8345221bfe52f4ecf30b8bdf26b4a213f5cd207ba3e4803",
        norm_state = "abc3041ac1452514dd38bea7dac56eaeb710139787ec11909029bb4ca7387977",
        norm_fields = c(best_k = "a7b8541c0f557ee0", control_genes_index = "1c497b478576b300", correlation_filtered_obj = "8f53dc8138ce299d", 
            correlation_threshold = "acf06dd9ce0cc864", correlation_vector = "9eaece3602f93d46", filtering_summary_text = "1ba7030e9325c39d", 
            final_filtering_plot = "1cc1786a46386211", final_qc_plot = "1252e1714ac6c9a2", normalization_complete = "21dd4c7bf548536a", 
            normalization_input_reference = "7e1a356bda20ce87", normalized_protein_obj = "7e1a356bda20ce87", publication_paths = "71557d1c8bac9bbe", 
            ruv_covariates_complete = "21dd4c7bf548536a", ruv_normalized_obj = "8f53dc8138ce299d", ruv_optimization_result = "eeae0389f229a0bd", 
            ruv_optimization_table = "bf419098d53c4e33"),
        ruv_result = "dd5a614c8a2494b1c66e75e08ce2a65e9259d557dcfa5a4058b9f720b5e8f4e6",
        workflow = "1aea9a4f60d498c1d539b6ced76fa6d4522eb4875db793ea53cf2a73fe2e049d",
        value = "8f53dc8138ce299d67fb2dbeb73332a20afa4f5e15c31dd7f3aefd09892e4f32",
        error = NA_character_,
        calls = c(assertDIAExportReady = 1L, runSavedOperation = 1L, commitDIAReferenceState = 1L, saveDIADisplayTable = 1L, 
            publishDIAActionOutputs = 1L, readDIAProteinProgress = 0L, updateProteinFiltering = 0L, writeDIADesignTable = 2L, 
            resolveDIAReferenceQcPlot = 0L, resolveProtNormCorrelationResultFilenames = 0L, generateProtNormCompositeQcFigure = 0L, 
            finalizeProtNormWorkflowState = 0L),
        callbacks = "ruv:60a06748b674"),
    limpa_correction_pass_through = list(
        current = "peptide_count_filtered",
        history = c("initial", "protein_dpc_quant", "protein_normalised", "peptide_count_filtered"),
        head_current = "peptide_count_filtered",
        head_history = c("initial", "protein_dpc_quant", "protein_normalised", "peptide_count_filtered"),
        states = c(protein_dpc_quant = "41dfd92a6007ea8f051ed9a447d3688f61902363a2b9aa84aedf596e4c0e9549", protein_normalised = "ab8d0aac28c16422de25afd8b64efc63a80bdf5670f34d7bfe526bc0d88edd64", 
            peptide_count_filtered = "75c25ed13486d80bb56c8e3691369dde2d61de37bcacecf3224cdc3573723fd0"),
        outputs = list(protein_dpc_quant = NULL, protein_normalised = NULL, peptide_count_filtered = NULL),
        products = c("output/dia-protein-filtering-counts:67049a631f7b624a", "output/dia-protein-qc-output:8b8101851eb54168", 
            "output/dia-protein-qc-output:deef48a7ebdbcb3a"),
        provenance = "b735d98157cc15fc39d183d25e744262d612b0f19d994ce5b073da963f5b60be",
        record = "6d3d46d9c054af2bb3b2f07469d89981e294f90d4f84a5d713938372b89e86d2",
        files = "8ab794fecc83e91d92f2fee89cd8fdc808bd44a3e8c693472e90be9cbf830156",
        norm_state = "d61f7e02b4c6a6d12fca78b9106ca07e9c39c3d7f406541c74a2ddc2947c3167",
        norm_fields = c(best_k = "3fe0774c75a5f7ec", control_genes_index = "3fe0774c75a5f7ec", filtering_summary_text = "1ba7030e9325c39d", 
            normalization_complete = "21dd4c7bf548536a", normalization_input_reference = "7e1a356bda20ce87", normalized_protein_obj = "7e1a356bda20ce87", 
            ruv_normalized_obj = "7e1a356bda20ce87", ruv_optimization_result = "71557d1c8bac9bbe"),
        ruv_result = "71557d1c8bac9bbe3cbec8d00bb223a2f372279827064095447e569fbf5a760a",
        workflow = "31f37cefd94f75f4f1e164b03125e1825178d30df82fc99a89aa16b7ba947e3a",
        value = "21dd4c7bf548536a0a6a9811519c53a3223f1a4f0357aa5382dd8a49aad258b6",
        error = NA_character_,
        calls = c(assertDIAExportReady = 1L, runSavedOperation = 0L, commitDIAReferenceState = 0L, saveDIADisplayTable = 0L, 
            publishDIAActionOutputs = 0L, readDIAProteinProgress = 0L, updateProteinFiltering = 0L, writeDIADesignTable = 0L, 
            resolveDIAReferenceQcPlot = 0L, resolveProtNormCorrelationResultFilenames = 0L, generateProtNormCompositeQcFigure = 0L, 
            finalizeProtNormWorkflowState = 0L),
        callbacks = "post:f08bb05015f2"),
    limpa_correction_correlation = list(
        current = "correlation_filtered",
        history = c("initial", "protein_dpc_quant", "protein_normalised", "peptide_count_filtered", "correlation_filtered"
            ),
        head_current = "correlation_filtered",
        head_history = c("initial", "protein_dpc_quant", "protein_normalised", "peptide_count_filtered", "correlation_filtered"
            ),
        states = c(protein_dpc_quant = "41dfd92a6007ea8f051ed9a447d3688f61902363a2b9aa84aedf596e4c0e9549", protein_normalised = "ab8d0aac28c16422de25afd8b64efc63a80bdf5670f34d7bfe526bc0d88edd64", 
            peptide_count_filtered = "75c25ed13486d80bb56c8e3691369dde2d61de37bcacecf3224cdc3573723fd0", correlation_filtered = "312ee5d40ac0d3f3c809c3864a1d440e1dbc000c8cc75a98cfdd9c95b96f5935"
            ),
        outputs = list(protein_dpc_quant = NULL, protein_normalised = NULL, peptide_count_filtered = NULL, correlation_filtered = c(correlation = TRUE)),
        products = c("dia-pearson-pairs:9eaece3602f93d46", "dia-pearson-selection:36ea8f5f2e373570", "dia_normalization_pca:d9b5720af8d8a298", 
            "output/dia-protein-filtering-counts:67049a631f7b624a", "output/dia-protein-qc-output:8b8101851eb54168", 
            "output/dia-protein-qc-output:deef48a7ebdbcb3a"),
        provenance = "9201164d5d72b3fe7f7ee72607a2e03edc64d9eb2a7d6c94b87505b9bd14cc62",
        record = "181955afdf5b3fc2465446465437e447f422ace7d0f85293e0d5a1e0fb446ff3",
        files = "34e64dc1e1e246298b07d49ed56b42a6743ca9e662a09925d3be94dc69c4735c",
        norm_state = "3b0fbc63724b4e8bf08fb33060738e1459d38205eb8bd692f03e358d03181fa1",
        norm_fields = c(best_k = "3fe0774c75a5f7ec", control_genes_index = "3fe0774c75a5f7ec", correlation_filtered_obj = "7e1a356bda20ce87", 
            correlation_threshold = "acf06dd9ce0cc864", correlation_vector = "9eaece3602f93d46", filtering_summary_text = "1ba7030e9325c39d", 
            final_filtering_plot = "1cc1786a46386211", final_qc_plot = "1252e1714ac6c9a2", normalization_complete = "21dd4c7bf548536a", 
            normalization_input_reference = "7e1a356bda20ce87", normalized_protein_obj = "7e1a356bda20ce87", ruv_normalized_obj = "7e1a356bda20ce87", 
            ruv_optimization_result = "71557d1c8bac9bbe"),
        ruv_result = "71557d1c8bac9bbe3cbec8d00bb223a2f372279827064095447e569fbf5a760a",
        workflow = "31f37cefd94f75f4f1e164b03125e1825178d30df82fc99a89aa16b7ba947e3a",
        value = "7e1a356bda20ce879fd8f423db352ae28d5e6b171ffb73d4cabec81acfa36136",
        error = NA_character_,
        calls = c(assertDIAExportReady = 1L, runSavedOperation = 2L, commitDIAReferenceState = 1L, saveDIADisplayTable = 0L, 
            publishDIAActionOutputs = 1L, readDIAProteinProgress = 1L, updateProteinFiltering = 1L, writeDIADesignTable = 1L, 
            resolveDIAReferenceQcPlot = 1L, resolveProtNormCorrelationResultFilenames = 1L, generateProtNormCompositeQcFigure = 0L, 
            finalizeProtNormWorkflowState = 0L),
        callbacks = character(0)),
    limpa_correction_ruv = list(
        current = "ruv_corrected",
        history = c("initial", "protein_dpc_quant", "protein_normalised", "peptide_count_filtered", "correlation_filtered", 
            "ruv_corrected"),
        head_current = "ruv_corrected",
        head_history = c("initial", "protein_dpc_quant", "protein_normalised", "peptide_count_filtered", "correlation_filtered", 
            "ruv_corrected"),
        states = c(protein_dpc_quant = "41dfd92a6007ea8f051ed9a447d3688f61902363a2b9aa84aedf596e4c0e9549", protein_normalised = "ab8d0aac28c16422de25afd8b64efc63a80bdf5670f34d7bfe526bc0d88edd64", 
            peptide_count_filtered = "75c25ed13486d80bb56c8e3691369dde2d61de37bcacecf3224cdc3573723fd0", correlation_filtered = "312ee5d40ac0d3f3c809c3864a1d440e1dbc000c8cc75a98cfdd9c95b96f5935", 
            ruv_corrected = "e1206824d8302e3581c57eeb39d5e7535da216dec6881f1b0575d03077085598"),
        outputs = list(protein_dpc_quant = NULL, protein_normalised = NULL, peptide_count_filtered = NULL, correlation_filtered = c(correlation = TRUE), 
                ruv_corrected = c("ruv-covariates" = TRUE)),
        products = c("dia-pearson-pairs:9eaece3602f93d46", "dia-pearson-selection:36ea8f5f2e373570", "dia-ruv-covariates:8cc4809c5a4a073b", 
            "dia_normalization_pca:d9b5720af8d8a298", "output/dia-protein-filtering-counts:67049a631f7b624a", "output/dia-protein-qc-output:8b8101851eb54168", 
            "output/dia-protein-qc-output:deef48a7ebdbcb3a"),
        provenance = "22d6232c9a55eab50e0d9d1c186d266c5e8bd33368999fa3af0e8ac982064470",
        record = "8c98ca3f561320405afd44b86f07d53fa0ceda8f03b8fc7e635869dfb2bcfeb1",
        files = "b8e5c65e86cc0c32a8345221bfe52f4ecf30b8bdf26b4a213f5cd207ba3e4803",
        norm_state = "25f3c43b9d26b351078901bd372606692d82831c945aa7393aa91c10e6d3b672",
        norm_fields = c(best_k = "a7b8541c0f557ee0", control_genes_index = "1c497b478576b300", correlation_filtered_obj = "8422092ca23581c1", 
            correlation_threshold = "acf06dd9ce0cc864", correlation_vector = "9eaece3602f93d46", filtering_summary_text = "1ba7030e9325c39d", 
            final_filtering_plot = "1cc1786a46386211", final_qc_plot = "1252e1714ac6c9a2", normalization_complete = "21dd4c7bf548536a", 
            normalization_input_reference = "7e1a356bda20ce87", normalized_protein_obj = "7e1a356bda20ce87", publication_paths = "71557d1c8bac9bbe", 
            ruv_covariates_complete = "21dd4c7bf548536a", ruv_normalized_obj = "8422092ca23581c1", ruv_optimization_result = "d91a80603a5612cd", 
            ruv_optimization_table = "bf419098d53c4e33"),
        ruv_result = "1a23a1e795d23d37e619847429154a4ed2c177b269d98f0e463d4e920589b81f",
        workflow = "3707fc6f2eed26843082dd0112a28f8f6ac8ebb623036e6eb1d80a9782cf17c0",
        value = "8422092ca23581c14e7503c24ccac8c235a1f244336f28163c3fe6d6bd7533f7",
        error = NA_character_,
        calls = c(assertDIAExportReady = 1L, runSavedOperation = 1L, commitDIAReferenceState = 1L, saveDIADisplayTable = 1L, 
            publishDIAActionOutputs = 1L, readDIAProteinProgress = 0L, updateProteinFiltering = 0L, writeDIADesignTable = 2L, 
            resolveDIAReferenceQcPlot = 0L, resolveProtNormCorrelationResultFilenames = 0L, generateProtNormCompositeQcFigure = 0L, 
            finalizeProtNormWorkflowState = 0L),
        callbacks = "ruv:e1c1cd4c99f5"),
    limpa_correction_ruv_repeat = list(
        current = "ruv_corrected",
        history = c("initial", "protein_dpc_quant", "protein_normalised", "peptide_count_filtered", "correlation_filtered", 
            "ruv_corrected"),
        head_current = "ruv_corrected",
        head_history = c("initial", "protein_dpc_quant", "protein_normalised", "peptide_count_filtered", "correlation_filtered", 
            "ruv_corrected"),
        states = c(protein_dpc_quant = "41dfd92a6007ea8f051ed9a447d3688f61902363a2b9aa84aedf596e4c0e9549", protein_normalised = "ab8d0aac28c16422de25afd8b64efc63a80bdf5670f34d7bfe526bc0d88edd64", 
            peptide_count_filtered = "75c25ed13486d80bb56c8e3691369dde2d61de37bcacecf3224cdc3573723fd0", correlation_filtered = "312ee5d40ac0d3f3c809c3864a1d440e1dbc000c8cc75a98cfdd9c95b96f5935", 
            ruv_corrected = "0c680e295e60bf593b7864d35762f69e81dd99629cd7662359be16c81285cb8d"),
        outputs = list(protein_dpc_quant = NULL, protein_normalised = NULL, peptide_count_filtered = NULL, correlation_filtered = c(correlation = TRUE), 
                ruv_corrected = c("ruv-covariates" = TRUE)),
        products = c("dia-pearson-pairs:9eaece3602f93d46", "dia-pearson-selection:36ea8f5f2e373570", "dia-ruv-covariates:8cc4809c5a4a073b", 
            "dia-ruv-covariates:8f30ef662db854de", "dia_normalization_pca:d9b5720af8d8a298", "output/dia-protein-filtering-counts:67049a631f7b624a", 
            "output/dia-protein-qc-output:8b8101851eb54168", "output/dia-protein-qc-output:deef48a7ebdbcb3a"),
        provenance = "22d6232c9a55eab50e0d9d1c186d266c5e8bd33368999fa3af0e8ac982064470",
        record = "548a545e77a0fb5d09953a780b89dc4ffbf5d9dd5cb217eee5a4847fab89fe01",
        files = "09d8cbcc30ea9902537d19ee36c7e73843b1f1638ce516c1853f6b39a0545f64",
        norm_state = "1d6de9a50ebc818aa87294c77c8f5847eed7f02e38bcbeb691544c0c5510aa75",
        norm_fields = c(best_k = "a7b8541c0f557ee0", control_genes_index = "19d490e85cf85a3f", correlation_filtered_obj = "aebfcf21f9e83b32", 
            correlation_threshold = "acf06dd9ce0cc864", correlation_vector = "9eaece3602f93d46", filtering_summary_text = "1ba7030e9325c39d", 
            final_filtering_plot = "1cc1786a46386211", final_qc_plot = "1252e1714ac6c9a2", normalization_complete = "21dd4c7bf548536a", 
            normalization_input_reference = "7e1a356bda20ce87", normalized_protein_obj = "7e1a356bda20ce87", publication_paths = "71557d1c8bac9bbe", 
            ruv_covariates_complete = "21dd4c7bf548536a", ruv_normalized_obj = "aebfcf21f9e83b32", ruv_optimization_result = "8faef6229524df22", 
            ruv_optimization_table = "4dd1aab74e24f986"),
        ruv_result = "3b26a46b119a2ad250fc602170716a5883a82760a3b572c920ff3d5f5766fe5f",
        workflow = "87e84fcf3901962127a4b85055f3d423299d2eadabe3987a1a98c90eb4f91d4e",
        value = "aebfcf21f9e83b325a4c92d9d1d81c8b1cda3d82cf6a607e1cc5e84a71a2e59b",
        error = NA_character_,
        calls = c(assertDIAExportReady = 1L, runSavedOperation = 1L, commitDIAReferenceState = 1L, saveDIADisplayTable = 1L, 
            publishDIAActionOutputs = 1L, readDIAProteinProgress = 0L, updateProteinFiltering = 0L, writeDIADesignTable = 2L, 
            resolveDIAReferenceQcPlot = 0L, resolveProtNormCorrelationResultFilenames = 0L, generateProtNormCompositeQcFigure = 0L, 
            finalizeProtNormWorkflowState = 0L),
        callbacks = "ruv:5e9abad6fbec")
)

test_that("KIT-111 W1 (a): the DIA normalisation and correlation clicks give the recorded states, products, provenance, files and display", {
    fixture <- protNormCoreOwner()
    clicks <- protNormW1Clicks(fixture)
    protNormW1Expect("dia_normalization_skip", clicks$normalization("skip"))
    # The same click again: the saved product is reused and nothing is committed.
    protNormW1Expect("dia_normalization_skip_repeat", clicks$normalization("skip"))
    protNormW1Expect("dia_normalization_manual", clicks$normalization("manual"))
    protNormW1Expect("dia_correlation_apply", clicks$correlation(0.5))
    # Re-applied at another threshold: the source is the normalisation's final, not the head.
    protNormW1Expect("dia_correlation_reapply", clicks$correlation(0.7))
    protNormW1Expect("dia_correlation_skip", clicks$correlation(0, skipped = TRUE))
    # A refusal (manual k above what the data allow), then a valid run on the same owner.
    protNormW1Expect("dia_normalization_refused", clicks$normalization("manual", ruv_k = 50))
    protNormW1Expect("dia_normalization_retry", clicks$normalization("manual"))
})

test_that("KIT-111 W1 (a): the DIA_LIMPA pass-through, correlation and RUV covariates in both modes give the recorded values", {
    skip_if_not_installed("limpa")
    skip_if_not_installed("mixOmics")
    for (mode in c("covariates", "correction")) {
        fixture <- protNormW1LimpaOwner()
        clicks <- protNormW1Clicks(fixture)
        protNormW1Expect(paste0("limpa_", mode, "_pass_through"), clicks$pass_through())
        protNormW1Expect(paste0("limpa_", mode, "_correlation"), clicks$correlation(0.5))
        protNormW1Expect(paste0("limpa_", mode, "_ruv"), clicks$ruv_covariates(mode))
        # The same click again, as today: its source is the committed RUV state (the runner sets
        # correlation_filtered_obj to it), so the declaration runs again and commits.
        protNormW1Expect(paste0("limpa_", mode, "_ruv_repeat"), clicks$ruv_covariates(mode))
    }
})

# ------------------------------------------------------------------------------------------------
# KIT-111 W2: the split runners on the managed worker, at core level. Each W1 (a) click runs here as the managed tab
# runs it: the gate and the request in the controller, the declaration through the worker's entry point
# (`workflowWorkerMain`) on the request's named input, the controller's acceptance and completion, then the
# publication half on the delivered product. The DIA_LIMPA pass-through and the skipped correlation stay in process on
# both paths. The captures must equal W1 (a)'s frozen values; the halves' call counts differ by design (the worker runs
# the compute half and the completion commits) and are not compared.
# ------------------------------------------------------------------------------------------------

#' One managed action, in this process: the frozen request (with its named input) through the worker's entry point,
#' the acceptance and what the tab reads. A worker failure is raised with the worker's message.
protNormW2Dispatch <- function(fixture, controller, operation, parameters, input, context = NULL) {
    store <- fixture$manager$.context$store
    id <- paste0("wfx-", substr(digest::digest(list(tempfile(), operation), algo = "sha256"), 1L, 32L))
    paths <- workflowExportPaths(store$directory, id)
    token <- list(request_id = id, incarnation = store$identity$incarnation, expected_head = workflowExecutionHead(store)$entity_id,
        generation = controller$control$generation, controller_id = controller$control$controller_id)
    writeWorkflowExecutionMessage(file.path(paths$request_dir, "publication-context.rds"),
        list(experiment_label = "core", experiment_paths = fixture$paths, omic_type = "proteomics"))
    writeWorkflowExportJournal(paths, list(schema = "multischolar.workflow_action.v1", request_id = id, status = "admitting",
        expected_head = token$expected_head, token_hash = digest::digest(token, algo = "sha256")))
    request <- list(schema = 1L, build = workflowExecutionBuildIdentity(), operation = operation,
        input_references = list(state = input), parameters = parameters, token = token, owner = store$identity,
        project_root = dirname(store$directory), directory = paths$request_dir, config = list(mode = "adapters"),
        worker_cap = 12 * 1024^3)
    path <- file.path(paths$request_dir, "request.rds")
    writeWorkflowExecutionMessage(path, request)
    withr::with_envvar(c(MULTISCHOLAR_WORKFLOW_WORKER = "1"), suppressWarnings(suppressMessages(utils::capture.output(try(
        workflowWorkerMain(c(path, digest::digest(file = path, algo = "sha256"))), silent = TRUE)))))
    message <- readRDS(file.path(paths$request_dir, "worker-result.rds"))
    if (!identical(message$status, "complete")) stop(message$error, call. = FALSE)
    controller$job <- list(request = list(token = token, operation = operation, config = list(mode = "adapters"),
        parameters = parameters), paths = paths, completion_context = validateWorkflowOmicsCompletionContext(operation, context))
    receipt <- acceptWorkflowActionResult(controller, message[c("schema", "status", "token", "result")])
    accepted <- resumeWorkflowAction(controller, id)
    c(readManagedOmicsResult(fixture$manager, accepted, NULL), list(acceptance = receipt$status, result = message$result))
}

#' The managed normalisation click: the gate and the request in the controller, the worker, the publication half.
protNormW2Normalization <- function(fixture, dispatch, norm, settings, post_qc, ruv_qc) {
    manager <- fixture$manager
    request <- withDIAStateScope(manager, function() {
        assertDIAExportReady(manager, allow_action = "normalization")
        diaNormalizationRequest(manager, settings, norm, "group")
    })
    state <- if (identical(settings$ruv_mode, "skip")) "normalized" else "ruv_corrected"
    result <- dispatch("normalization-ruv", request$parameters, request$source, state)
    withDIAStateScope(manager, function() publishDIAReferenceNormalization(manager, settings, result[c("value", "reference")],
        fixture$workflow, norm, fixture$paths, "proteomics", "core", post_qc, ruv_qc))
}

#' The managed correlation click: the pairs, then the selection on them, each on the normalisation's final.
protNormW2Correlation <- function(fixture, dispatch, norm, threshold, aesthetics) {
    manager <- fixture$manager
    source <- withDIAStateScope(manager, function() {
        assertDIAExportReady(manager, allow_action = "correlation")
        source <- norm$ruv_normalized_obj
        if (!isSavedWorkflowReference(source)) cli::cli_abort("normalization must complete before correlation filtering")
        readDIAStateMetadata(manager, source)
        source
    })
    pairs <- dispatch("pearson-pairs", diaPearsonPairsParameters("group"), source)
    selection <- dispatch("pearson-selection", list(pairs = pairs$reference$entity_id, threshold = threshold), source,
        "correlation_filtered")
    withDIAStateScope(manager, function() publishDIAReferenceCorrelation(manager, list(pairs = pairs$reference,
        final = selection$value), fixture$workflow, norm, fixture$paths, threshold, FALSE, "proteomics", "core", aesthetics))
}

#' The managed RUV covariates click (DIA_LIMPA).
protNormW2RuvCovariates <- function(fixture, dispatch, norm, mode, ruv_qc) {
    manager <- fixture$manager
    request <- withDIAStateScope(manager, function() {
        assertDIAExportReady(manager, allow_action = "ruv-covariates")
        source <- norm$correlation_filtered_obj %||% getDIAStateReference(manager)
        if (!isSavedWorkflowReference(source)) cli::cli_abort("the correlation step must complete before the RUV covariates")
        list(source = source, parameters = resolveDIARuvCovariatesParameters(fixture$workflow$config_list, "group", mode))
    })
    result <- dispatch("protein_ruv_covariates", request$parameters, request$source,
        if (identical(mode, "correction")) "ruv_corrected" else "ruv_covariates")
    withDIAStateScope(manager, function() publishDIAReferenceRuvCovariates(manager, list(ruv_mode = mode,
        color_variable = "group", shape_variable = NULL), result[c("value", "reference")], fixture$workflow, norm,
        fixture$paths, ruv_qc))
}

#' The tab's clicks on one owner through the worker, each captured as protNormW1Capture describes. The dispatch
#' context names the state with a placeholder config and description, which the handler's commit description replaces.
#' `results()` lists each worker value (its operation, size, names, commit description and input).
protNormW2Clicks <- function(fixture) {
    norm <- new.env(parent = emptyenv())
    seen <- character()
    results <- list()
    manager <- fixture$manager
    controller <- omicsWorkerController(list(owner = manager))
    callback <- function(kind) function(reference) {
        seen <<- c(seen, paste0(kind, ":", substr(protNormW1Sha(protNormW1State(manager, reference)), 1L, 12L)))
        invisible(NULL)
    }
    capture <- function(run) {
        out <- protNormW1Capture(fixture, norm, run, seen)
        seen <<- character()
        out
    }
    dispatch <- function(operation, parameters, input, state_name = NULL) {
        context <- if (!is.null(state_name)) list(state_name = state_name, config = list(), description = "dispatch")
        out <- protNormW2Dispatch(fixture, controller, operation, parameters, input, context)
        results[[length(results) + 1L]] <<- list(operation = operation, bytes = as.numeric(object.size(out$result)),
            names = names(out$result), commit = out$result$commit, input = input$entity_id)
        out
    }
    aesthetics <- function() list(color_var = "group", shape_var = NULL)
    list(norm = norm, results = function() results,
        normalization = function(mode, ...) {
            settings <- utils::modifyList(protNormCoreSettings(mode), list(...))
            capture(protNormW1Run(protNormW2Normalization(fixture, dispatch, norm, settings, callback("post"), callback("ruv"))))
        },
        correlation = function(threshold, skipped = FALSE) {
            capture(protNormW1Run(if (skipped) runDIAReferenceCorrelation(fixture$workflow, norm, fixture$paths, threshold,
                TRUE, "proteomics", "core", function() NULL, aesthetics) else
                protNormW2Correlation(fixture, dispatch, norm, threshold, aesthetics)))
        },
        pass_through = function() {
            capture(protNormW1Run(runDIAReferenceLimpaNormalisationPassThrough(protNormCoreSettings("covariates"),
                fixture$workflow, norm, callback("post"))))
        },
        ruv_covariates = function(mode) capture(protNormW1Run(protNormW2RuvCovariates(fixture, dispatch, norm, mode, callback("ruv")))))
}

#' Compare one managed case with W1 (a)'s frozen in-process value, every field but the halves' call counts.
protNormW2Expect <- function(case, capture) {
    frozen <- protNormW1Frozen()[[case]]
    expect_false(is.null(frozen), info = case)
    for (field in setdiff(union(names(frozen), names(capture)), "calls")) {
        expect_identical(capture[[field]], frozen[[field]], info = paste("worker", case, field))
    }
}

protNormW2Mocks <- function(env = parent.frame()) {
    testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
        finishWorkflowExecutionJob = function(owner, status) { owner$job <- NULL; list(status = status) },
        .package = "MultiScholaR", .env = env)
}

test_that("KIT-111 W2 (W1 (b)): the DIA normalisation and correlation clicks through the worker give W1 (a)'s values", {
    skip_if_not(workflowOwnerCapable())
    protNormW2Mocks()
    fixture <- protNormCoreOwner()
    clicks <- protNormW2Clicks(fixture)
    protNormW2Expect("dia_normalization_skip", clicks$normalization("skip"))
    # A repeated click: the worker's hit carries a NULL commit description and the completion commits nothing.
    protNormW2Expect("dia_normalization_skip_repeat", clicks$normalization("skip"))
    protNormW2Expect("dia_normalization_manual", clicks$normalization("manual"))
    protNormW2Expect("dia_correlation_apply", clicks$correlation(0.5))
    # R6: the pairs' provenance stays pending until the selection commits, and is recorded before it.
    provenance <- withDIAStateScope(fixture$manager, function() fixture$manager$.context$snapshot$step_provenance$correlation_filtered)
    expect_identical(provenance$operation, "pearson-selection")
    expect_identical(lapply(provenance$preceded_by, `[[`, "operation"), list("pearson-pairs"))
    expect_identical(provenance$preceded_by[[1L]]$values, list(group = "group", tech_rep_remove_regex = "pool"))
    # The re-run names its source (the normalisation's final), which is not the head.
    protNormW2Expect("dia_correlation_reapply", clicks$correlation(0.7))
    protNormW2Expect("dia_correlation_skip", clicks$correlation(0, skipped = TRUE))
    protNormW2Expect("dia_normalization_refused", clicks$normalization("manual", ruv_k = 50))
    protNormW2Expect("dia_normalization_retry", clicks$normalization("manual"))
    results <- clicks$results()
    expect_identical(vapply(results, `[[`, "", "operation"), c("normalization-ruv", "normalization-ruv", "normalization-ruv",
        "pearson-pairs", "pearson-selection", "pearson-pairs", "pearson-selection", "normalization-ruv"))
    expect_null(results[[2L]]$commit)
    expect_true("commit" %in% results[[2L]]$names)
    expect_identical(results[[1L]]$commit$description, "Post-normalization complete: RUV-III skipped by user")
    # The reapply's input is the normalisation's final, as the first apply's was; neither is the head.
    expect_identical(results[[6L]]$input, results[[4L]]$input)
    # R4: the values cross the boundary as references and the commit description.
    expect_true(all(vapply(results, `[[`, 0, "bytes") < 8388608))
})

test_that("KIT-111 W2 (W1 (b)): the DIA_LIMPA correlation and RUV covariates through the worker give W1 (a)'s values", {
    skip_if_not(workflowOwnerCapable())
    skip_if_not_installed("limpa")
    skip_if_not_installed("mixOmics")
    protNormW2Mocks()
    for (mode in c("covariates", "correction")) {
        fixture <- protNormW1LimpaOwner()
        clicks <- protNormW2Clicks(fixture)
        protNormW2Expect(paste0("limpa_", mode, "_pass_through"), clicks$pass_through())
        protNormW2Expect(paste0("limpa_", mode, "_correlation"), clicks$correlation(0.5))
        protNormW2Expect(paste0("limpa_", mode, "_ruv"), clicks$ruv_covariates(mode))
        # F111-02, today's behaviour on both paths: the repeat's source is the committed RUV state, and it commits again.
        protNormW2Expect(paste0("limpa_", mode, "_ruv_repeat"), clicks$ruv_covariates(mode))
        covariates <- Filter(function(result) identical(result$operation, "protein_ruv_covariates"), clicks$results())
        expect_length(covariates, 2L)
        # R4, on the covariates case: under the worker value's 8 MiB.
        expect_true(all(vapply(covariates, `[[`, 0, "bytes") < 8388608), info = mode)
        expect_identical(covariates[[1L]]$commit$state_name, if (identical(mode, "correction")) "ruv_corrected" else "ruv_covariates")
    }
})

test_that("KIT-111 W2: a DIA worker value whose commit description is NULL commits nothing; one that carries it replaces the context's (extension 1)", {
    skip_if_not(workflowOwnerCapable())
    protNormW2Mocks()
    fixture <- protNormCoreOwner()
    manager <- fixture$manager
    clicks <- protNormW2Clicks(fixture)
    clicks$normalization("skip")
    record <- readDIAStateMetadata(manager)
    expect_identical(record$description, "Post-normalization complete: RUV-III skipped by user")
    expect_identical(record$config$ruv_applied, FALSE)
    head <- withDIAStateScope(manager, function() manager$.context$snapshot$states$normalized)
    counted <- protNormCoreCountCalls("commitWorkflowReferenceState", clicks$normalization("skip"))
    expect_identical(counted$calls[["commitWorkflowReferenceState"]], 0L)
    expect_identical(withDIAStateScope(manager, function() manager$.context$snapshot$states$normalized), head)
    expect_identical(as.character(unlist(manager$state_history)), c("initial", "protein_replicate_filtered", "normalized"))
    # The completion's own rule: a NULL commit takes the product-only branch whatever the context names.
    entry <- workflowActionSpec("normalization-ruv")
    value <- list(checkpoint = record$data, reference = record$data$products$normalization_science, commit = NULL)
    counted <- protNormCoreCountCalls("commitWorkflowReferenceState", withDIAStateScope(manager, function()
        completeWorkflowOmicsAction(manager, entry, value, list(state_name = "normalized", config = list(), description = "x"))))
    expect_null(counted$value)
    expect_identical(counted$calls[["commitWorkflowReferenceState"]], 0L)
})

test_that("KIT-111 W2 (R7): a managed publication that fails is completed by the repeated click, a hit that commits nothing", {
    skip_if_not(workflowOwnerCapable())
    protNormW2Mocks()
    fixture <- protNormCoreOwner()
    clicks <- protNormW2Clicks(fixture)
    failed <- testthat::with_mocked_bindings(clicks$normalization("skip"),
        generateProtNormCompositeQcFigure = function(...) stop("forced publication failure"), .package = "MultiScholaR")
    expect_identical(failed$error, "forced publication failure")
    expect_identical(failed$head_history, c("initial", "protein_replicate_filtered", "normalized"))
    protNormW2Expect("dia_normalization_skip", clicks$normalization("skip"))
    results <- clicks$results()
    expect_length(results, 2L)
    expect_null(results[[2L]]$commit)
})

test_that("KIT-111 W2: a managed action names its input; without one its input is the head (extension 2)", {
    fixture <- protNormCoreOwner()
    manager <- fixture$manager
    protNormCoreQuiet(runDIAReferenceNormalization(protNormCoreSettings("skip"), fixture$workflow, new.env(), fixture$paths,
        "proteomics", "core", function(reference) NULL, function(reference) NULL, function() "group"))
    source <- getDIAStateReference(manager, "protein_replicate_filtered")
    head <- getDIAStateReference(manager)
    expect_false(identical(source$entity_id, head$entity_id))
    seen <- list()
    testthat::local_mocked_bindings(
        workflowSessionApplication = function(session) list(),
        beginWorkflowGUIAction = function(application, manager, operation, publication) list2env(list(id = operation, closed = FALSE)),
        finishWorkflowGUIAction = function(reservation, status, ...) { reservation$closed <- TRUE; list() },
        newWorkflowExecutionController = function(manager, application) list(manager = manager),
        runWorkflowAction = function(operation, input_references, parameters, owner, reservation, completion_context = NULL) {
            seen[[length(seen) + 1L]] <<- list(operation = operation, input = input_references$state$entity_id)
            paste0("worker-", length(seen))
        },
        pollWorkflowAction = function(owner) list(status = "running"),
        .package = "MultiScholaR")
    shiny::testServer(function(input, output, session) {
        session$userData$runner <- installWorkflowOmicsActions(session, fixture$workflow, "proteomics", "core", fixture$paths)
    }, {
        runner <- session$userData$runner
        pairs <- list(group = "group", tech_rep_remove_regex = "pool")
        runner$run(list(list(operation = "pearson-pairs", parameters = pairs, input = source, deliver = function(result) NULL)))
        runner$step$cancel("test")
        runner$run(list(list(operation = "pearson-pairs", parameters = pairs, deliver = function(result) NULL)))
        runner$step$cancel("test")
    })
    expect_identical(seen, list(list(operation = "pearson-pairs", input = source$entity_id),
        list(operation = "pearson-pairs", input = head$entity_id)))
})

test_that("KIT-111 W2: the DIA adaptor's norm entries name registered omics operations and existing functions (D4)", {
    norm <- proteomicsDiaBinding$norm
    expect_false(is.null(norm))
    operations <- c(norm$run$operation, norm$correlation$pairs, norm$correlation$selection, norm$ruv_covariates$operation)
    expect_identical(operations, c("normalization-ruv", "pearson-pairs", "pearson-selection", "protein_ruv_covariates"))
    for (operation in operations) {
        entry <- workflowActionSpec(operation)
        expect_identical(entry$compute, "omics", info = operation)
        expect_true(exists(entry$handler %||% "", envir = asNamespace("MultiScholaR"), mode = "function"), info = operation)
    }
    calls <- c(run = "runDIAReferenceNormalization", run_publish = "publishDIAReferenceNormalization",
        pass_through = "runDIAReferenceLimpaNormalisationPassThrough", correlation = "runDIAReferenceCorrelation",
        correlation_publish = "publishDIAReferenceCorrelation", ruv_covariates = "runDIAReferenceRuvCovariates",
        ruv_covariates_publish = "publishDIAReferenceRuvCovariates")
    entries <- list(run = norm$run$in_process, run_publish = norm$run$publish, pass_through = norm$pass_through,
        correlation = norm$correlation$in_process, correlation_publish = norm$correlation$publish,
        ruv_covariates = norm$ruv_covariates$in_process, ruv_covariates_publish = norm$ruv_covariates$publish)
    for (name in names(calls)) {
        expect_true(is.function(entries[[name]]), info = name)
        if (!is.function(entries[[name]])) next
        expect_identical(as.character(body(entries[[name]])[[1L]]), calls[[name]], info = name)
        expect_true(exists(calls[[name]], envir = asNamespace("MultiScholaR"), mode = "function"), info = name)
    }
    # The configuration stays valid and proteomics does not join the shared normalisation frame (D5).
    expect_null(omicGuiConfig("proteomics")$norm)
})

# ------------------------------------------------------------------------------------------------
# KIT-111 W3: the tab through its public entry point (`mod_prot_norm_server`) in testServer. Unmanaged, the session has
# no runner and the tab runs in process (today's path). Managed, the shell installs the omics runner and the tab
# dispatches through it under `omicsManagedChainMocks`, with `runWorkflowAction` running each request on its named input
# through the worker's entry point (`protNormW2Dispatch`), as the real one does since F111-03. MCI-008's module fixture is
# an in-memory manager, which the managed step refuses ("the managed worker needs the on-disk state owner"), so the tab
# runs on W1's owners with W1 (a)'s settings. Captured after each click: the history, the committed states as W1 (a)
# records them, the head's products, the provenance, the record's steps, the files, the tab's fields, the bus, the
# summary texts and the notices; managed, each action's operation and input.
# ------------------------------------------------------------------------------------------------

protNormW3Inputs <- function(mode, ...) utils::modifyList(list(norm_method = "quantile", ruv_mode = mode,
    ruv_percentage = 50, ruv_k = 1, color_variable = "group", ruv_grouping_variable = "group",
    min_pearson_correlation_threshold = 0.5), list(...))

#' The notices shown (`[type duration] text`) and the ids removed, recorded in `log`.
protNormW3Notices <- function(log, root, env = parent.frame()) {
    notify <- function(ui, action = NULL, duration = 5, closeButton = TRUE, id = NULL,
        type = c("default", "message", "warning", "error"), session = NULL) {
        log$notes <- c(log$notes, sprintf("[%s %s] %s", type[[1L]], if (is.null(duration)) "NULL" else format(duration),
            protNormW1Mask(paste(as.character(ui), collapse = " "), root)))
        invisible(id %||% "notification")
    }
    remove <- function(id, session = NULL) {
        log$removed <- c(log$removed, id)
        invisible(id)
    }
    testthat::local_mocked_bindings(showNotification = notify, removeNotification = remove, .package = "shiny", .env = env)
    testthat::local_mocked_bindings(showNotification = notify, removeNotification = remove, .package = "MultiScholaR", .env = env)
}

#' The managed session: the chain mocks, with each action's request run on its named input and recorded in `probe`
#' (`actions`: operation and input entity id; `values`: the worker's values). `probe$hold` keeps the action running.
protNormW3Managed <- function(fixture, probe, env = parent.frame()) {
    controller <- omicsManagedChainMocks(list(owner = fixture$manager, paths = fixture$paths), probe, env = env)
    testthat::local_mocked_bindings(
        runWorkflowAction = function(operation, input_references, parameters, owner, reservation, completion_context = NULL) {
            probe$actions <- c(probe$actions, list(list(operation = operation, input = input_references$state$entity_id)))
            out <- protNormW2Dispatch(fixture, controller, operation, parameters, input_references$state, completion_context)
            probe$values <- c(probe$values, list(out$result))
            probe$handoff <- list(accepted = list(operation = operation, reference = out$reference,
                state_reference = out$stateReference, display = out$display))
            paste0("worker-", operation)
        },
        pollWorkflowAction = function(owner) list(status = if (isTRUE(probe$hold)) "running" else "complete"),
        .package = "MultiScholaR", .env = env)
    controller
}

#' The capture after one click (see the section's header). `head_id` is local to the store.
protNormW3Capture <- function(fixture, norm, output, log, actions) {
    manager <- fixture$manager
    root <- fixture$root
    snapshot <- withWorkflowStateScope(manager, function() manager$.context$snapshot)
    history <- as.character(unlist(snapshot$state_history))
    committed <- setdiff(history, "initial")
    fields <- shiny::isolate(shiny::reactiveValuesToList(norm))
    # The RUV covariates' publication leaves its QC paths in the staging directory, named by store-local ids (draft
    # F111-05), on both paths; masked to the staged role.
    if (!is.null(fields$qc_plot_paths)) fields$qc_plot_paths <- rapply(list(fields$qc_plot_paths), function(path)
        sub("^.*/[.]workflow-export/[^/]+/files/", "<staged>/", path), classes = "character", how = "replace")[[1L]]
    workflow <- mget(setdiff(sort(ls(fixture$workflow)), c("state_manager", "design_matrix", "config_list",
        "state_update_trigger")), envir = fixture$workflow)
    text <- function(name) tryCatch(protNormW1Mask(output[[name]], root), error = function(e) NA_character_)
    list(
        current = manager$current_state,
        history = as.character(unlist(manager$state_history)),
        head_current = snapshot$current_state,
        head_history = history,
        states = vapply(stats::setNames(nm = committed), function(state)
            protNormW1Sha(protNormW1Record(manager, snapshot$states[[state]], root)), ""),
        products = protNormW1Products(manager, root),
        provenance = protNormW1Sha(protNormW1Provenance(snapshot, root)),
        record = protNormW1Sha(protNormW1RecordSteps(manager, root)),
        files = protNormCoreFiles(root),
        norm = vapply(fields[sort(names(fields))], function(field)
            substr(protNormW1Sha(protNormW1Canonical(field, manager, root)), 1L, 16L), ""),
        workflow = protNormW1Sha(protNormW1Canonical(workflow, manager, root)),
        outputs = c(correlation = text("correlation_filter_summary"), ruv_covariates = text("ruv_covariates_summary")),
        notes = log$notes,
        head_id = snapshot$states[[snapshot$current_state]]$entity_id,
        actions = actions)
}

#' Drive the tab through `steps` (named `function(session)`), capturing after each; managed under the runner.
protNormW3Drive <- function(fixture, steps, managed, probe = new.env()) {
    log <- new.env(parent = emptyenv())
    protNormW3Notices(log, fixture$root)
    if (managed) protNormW3Managed(fixture, probe)
    acc <- new.env(parent = emptyenv())
    acc$captures <- list()
    protNormCoreQuiet(utils::capture.output(shiny::testServer(mod_prot_norm_server, args = list(
        workflow_data = fixture$workflow, experiment_paths = fixture$paths, omic_type = "proteomics",
        experiment_label = "core"), {
        norm <- session$getReturned()
        for (case in names(steps)) {
            seen <- length(probe$actions)
            log$notes <- character()
            steps[[case]](session)
            session$flushReact()
            acc$captures[[case]] <- protNormW3Capture(fixture, norm, output, log,
                probe$actions[seq_along(probe$actions) > seen])
        }
    })))
    acc$captures
}

#' A click on `button` (after setting `inputs`), counted per drive.
protNormW3Clicker <- function() {
    counts <- new.env(parent = emptyenv())
    function(button, ...) {
        inputs <- list(...)
        function(session) {
            if (length(inputs)) do.call(session$setInputs, inputs)
            counts[[button]] <- (counts[[button]] %||% 0) + 1
            do.call(session$setInputs, stats::setNames(list(counts[[button]]), button))
        }
    }
}

#' W1 (a)'s DIA clicks as the tab's buttons.
protNormW3DiaSteps <- function() {
    click <- protNormW3Clicker()
    list(
        dia_normalization_skip = do.call(click, c(list("run_normalization"), protNormW3Inputs("skip"))),
        dia_normalization_skip_repeat = click("run_normalization"),
        dia_normalization_manual = click("run_normalization", ruv_mode = "manual"),
        dia_correlation_apply = click("apply_correlation_filter"),
        dia_correlation_reapply = click("apply_correlation_filter", min_pearson_correlation_threshold = 0.7),
        dia_correlation_skip = click("skip_correlation_filter"),
        dia_normalization_refused = click("run_normalization", ruv_k = 50),
        dia_normalization_retry = click("run_normalization", ruv_k = 1))
}

#' W1 (a)'s DIA_LIMPA clicks in one RUV mode as the tab's buttons.
protNormW3LimpaSteps <- function(mode) {
    click <- protNormW3Clicker()
    stats::setNames(list(do.call(click, c(list("run_normalization"), protNormW3Inputs(mode))),
        click("apply_correlation_filter"), click("apply_ruv_covariates"), click("apply_ruv_covariates")),
        paste0("limpa_", mode, c("_pass_through", "_correlation", "_ruv", "_ruv_repeat")))
}

#' The fields of a module capture that W1 (a)'s core capture also records and the tab's own work does not change.
protNormW3W1Fields <- c("current", "history", "head_current", "head_history", "states")

protNormW3ExpectW1 <- function(captures, label) {
    for (case in names(captures)) {
        frozen <- protNormW1Frozen()[[case]]
        expect_false(is.null(frozen), info = case)
        for (field in protNormW3W1Fields) {
            expect_identical(captures[[case]][[field]], frozen[[field]], info = paste(label, case, field))
        }
    }
}

#' The notices without the runner's own (its working notice, and its error before the tab's), each error by its text up
#' to the first colon: a worker's error reaches the tab as the controller words it ("managed action ended: <status>"
#' through the service), not as the kernel's condition.
protNormW3TabNotes <- function(notes) {
    notes <- notes[!grepl("^\\[(default NULL|error 15)\\] (Applying|Error applying) proteomics step", notes)]
    sub("(?s)^(\\[error [0-9]+\\] [^:]+:).*$", "\\1", notes, perl = TRUE)
}

#' Managed against unmanaged: every field but the store-local head id and the actions; the notices without the runner's.
protNormW3ExpectSame <- function(managed, unmanaged) {
    expect_identical(names(managed), names(unmanaged))
    for (case in names(unmanaged)) {
        for (field in setdiff(names(unmanaged[[case]]), c("notes", "head_id", "actions", "norm"))) {
            expect_identical(managed[[case]][[field]], unmanaged[[case]][[field]], info = paste("managed", case, field))
        }
        expect_identical(names(managed[[case]]$norm), names(unmanaged[[case]]$norm), info = paste("managed", case, "norm"))
        for (field in names(unmanaged[[case]]$norm)) {
            expect_identical(managed[[case]]$norm[[field]], unmanaged[[case]]$norm[[field]], info = paste("managed", case, "norm", field))
        }
        expect_identical(protNormW3TabNotes(managed[[case]]$notes), protNormW3TabNotes(unmanaged[[case]]$notes),
            info = paste("managed", case, "notes"))
    }
}

#' The unmanaged drive on a fresh owner, once per key in the file (the managed blocks compare with it).
protNormW3Unmanaged <- local({
    cache <- list()
    function(key, owner, steps) {
        if (is.null(cache[[key]])) {
            drive <- function() protNormW3Drive(owner(environment()), steps, managed = FALSE)
            cache[[key]] <<- drive()
        }
        cache[[key]]
    }
})

test_that("KIT-111 W3: unmanaged, the tab's run, apply, skip and refusal give W1 (a)'s states and today's notices", {
    captures <- protNormW3Unmanaged("dia", protNormCoreOwner, protNormW3DiaSteps())
    protNormW3ExpectW1(captures, "unmanaged")
    expect_identical(captures$dia_normalization_skip$notes, paste("[message 10]", buildProtNormCompletionNotification("skip")))
    expect_identical(captures$dia_normalization_manual$notes, paste("[message 10]", buildProtNormCompletionNotification("manual")))
    expect_identical(captures$dia_correlation_apply$notes,
        "[message 5] Correlation filtering completed! Ready for differential expression analysis.")
    expect_identical(captures$dia_correlation_skip$notes,
        "[message 5] Correlation filtering skipped! Ready for differential expression analysis.")
    expect_length(captures$dia_normalization_refused$notes, 1L)
    expect_match(captures$dia_normalization_refused$notes, "^\\[error 10\\] Error in normalization: ")
    for (case in names(captures)) expect_length(captures[[case]]$actions, 0L)
})

test_that("KIT-111 W3: managed, the tab's run, apply and refusal run on the worker on their named inputs and equal the unmanaged tab", {
    skip_if_not(workflowOwnerCapable())
    unmanaged <- protNormW3Unmanaged("dia", protNormCoreOwner, protNormW3DiaSteps())
    fixture <- protNormCoreOwner()
    probe <- new.env()
    managed <- protNormW3Drive(fixture, protNormW3DiaSteps(), managed = TRUE, probe = probe)
    protNormW3ExpectW1(managed, "managed")
    protNormW3ExpectSame(managed, unmanaged)
    operations <- lapply(managed, function(capture) vapply(capture$actions, `[[`, "", "operation"))
    expect_identical(operations, list(dia_normalization_skip = "normalization-ruv", dia_normalization_skip_repeat = "normalization-ruv",
        dia_normalization_manual = "normalization-ruv", dia_correlation_apply = c("pearson-pairs", "pearson-selection"),
        dia_correlation_reapply = c("pearson-pairs", "pearson-selection"), dia_correlation_skip = character(),
        dia_normalization_refused = "normalization-ruv", dia_normalization_retry = "normalization-ruv"))
    # The completion texts are today's, after the runner's working notice.
    working <- "[default NULL] Applying proteomics step ..."
    expect_identical(managed$dia_normalization_skip$notes, c(working, paste("[message 10]", buildProtNormCompletionNotification("skip"))))
    expect_identical(managed$dia_normalization_manual$notes, c(working, paste("[message 10]", buildProtNormCompletionNotification("manual"))))
    expect_identical(managed$dia_correlation_apply$notes, c(working, working,
        "[message 5] Correlation filtering completed! Ready for differential expression analysis."))
    # The worker's failure: the runner's error, then the tab's; nothing committed (W1 (a)'s states, above); the retry completes.
    refused <- managed$dia_normalization_refused$notes
    expect_length(refused, 3L)
    expect_match(refused[[2L]], "^\\[error 15\\] Error applying proteomics step: ")
    expect_match(refused[[3L]], "^\\[error 10\\] Error in normalization: ")
    # Each action's input: the normalisation's recorded source; the normalisation's final for the pairs and the selection,
    # on the reapply too, when the head is correlation_filtered.
    inputs <- lapply(managed, function(capture) vapply(capture$actions, `[[`, "", "input"))
    states <- withDIAStateScope(fixture$manager, function() fixture$manager$.context$snapshot$states)
    expect_identical(unique(unlist(inputs[grep("^dia_normalization", names(inputs))])), states$protein_replicate_filtered$entity_id)
    expect_identical(unique(inputs$dia_correlation_apply), inputs$dia_correlation_reapply[[1L]])
    expect_identical(unique(inputs$dia_correlation_reapply), inputs$dia_correlation_apply[[1L]])
    expect_false(identical(inputs$dia_correlation_reapply[[1L]], managed$dia_correlation_apply$head_id))
    # The repeated click is a hit whose commit description is NULL.
    expect_identical(vapply(probe$values[1:2], function(value) is.null(value$commit), NA), c(FALSE, TRUE))
})

test_that("KIT-111 W3: the DIA_LIMPA RUV covariates, unmanaged and managed, give W1 (a)'s states, the summary and the inputs", {
    skip_if_not_installed("limpa")
    skip_if_not_installed("mixOmics")
    for (mode in c("covariates", "correction")) {
        unmanaged <- protNormW3Unmanaged(mode, protNormW1LimpaOwner, protNormW3LimpaSteps(mode))
        protNormW3ExpectW1(unmanaged, "unmanaged")
        state <- if (identical(mode, "correction")) "ruv_corrected" else "ruv_covariates"
        summary <- unmanaged[[paste0("limpa_", mode, "_ruv")]]$outputs[["ruv_covariates"]]
        expect_true(endsWith(summary, paste0("State saved as: '", state, "'")), info = mode)
        if (!workflowOwnerCapable()) next
        fixture <- protNormW1LimpaOwner()
        managed <- protNormW3Drive(fixture, protNormW3LimpaSteps(mode), managed = TRUE)
        protNormW3ExpectW1(managed, "managed")
        protNormW3ExpectSame(managed, unmanaged)
        actions <- lapply(managed, `[[`, "actions")
        expect_identical(unname(lapply(actions, function(list) vapply(list, `[[`, "", "operation"))), list(character(),
            c("pearson-pairs", "pearson-selection"), "protein_ruv_covariates", "protein_ruv_covariates"), info = mode)
        # The covariates' input is the tab's correlation result (the selection's checkpoint, as in process); F111-02, as
        # today: the repeat's is the committed RUV state.
        expect_false(identical(actions[[3L]][[1L]]$input, actions[[4L]][[1L]]$input), info = mode)
        expect_identical(actions[[4L]][[1L]]$input, managed[[3L]]$head_id, info = mode)
    }
})

test_that("KIT-111 W3 (R5): a click while the tab's action runs waits for it; a click the application refuses commits nothing", {
    skip_if_not(workflowOwnerCapable())
    fixture <- protNormCoreOwner()
    probe <- new.env()
    probe$hold <- TRUE
    log <- new.env(parent = emptyenv())
    protNormW3Notices(log, fixture$root)
    protNormW3Managed(fixture, probe)
    done <- paste("[message 10]", buildProtNormCompletionNotification("skip"))
    protNormCoreQuiet(utils::capture.output(shiny::testServer(mod_prot_norm_server, args = list(
        workflow_data = fixture$workflow, experiment_paths = fixture$paths, omic_type = "proteomics",
        experiment_label = "core"), {
        do.call(session$setInputs, protNormW3Inputs("skip"))
        session$setInputs(run_normalization = 1)
        session$setInputs(run_normalization = 2)
        # The second click waits: one action reserved, nothing delivered.
        expect_identical(probe$operations, "normalization-ruv")
        expect_identical(log$notes, "[default NULL] Applying proteomics step ...")
        probe$hold <- FALSE
        session$elapse(100)
        session$flushReact()
        expect_identical(probe$operations, c("normalization-ruv", "normalization-ruv"))
        expect_identical(sum(log$notes == done), 2L)
        expect_identical(as.character(unlist(fixture$manager$state_history)), c("initial", "protein_replicate_filtered", "normalized"))
        expect_null(probe$values[[2L]]$commit)
        # Another heavy action holds the application: the click is refused with that notice and the tab's error.
        log$notes <- character()
        testthat::with_mocked_bindings({
            session$setInputs(ruv_mode = "manual")
            session$setInputs(run_normalization = 3)
        }, beginWorkflowGUIAction = function(application, manager, operation, publication)
            cli::cli_abort("another heavy action is active", class = "workflow_execution_busy"), .package = "MultiScholaR")
        expect_identical(log$notes, c("[error 15] another heavy action is active",
            "[error 10] Error in normalization: another heavy action is active"))
        expect_identical(fixture$manager$current_state, "normalized")
        expect_length(probe$actions, 2L)
    })))
})
