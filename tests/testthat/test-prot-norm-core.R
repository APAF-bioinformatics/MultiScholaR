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
