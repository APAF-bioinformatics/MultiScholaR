# KIT-105 W2: the lipid/metab normalisation computation in the core, written before it existed. Each core function is
# compared with today's lipid and metab functions on the MCI-021/MCI-015 fixtures (design D7: lipid's and metab's
# differences are kept per omic): the composite figure and the composite from files (utils_workflow_omics_figures.R),
# the final-QC PCA and the in-process normalisation sequence (utils_workflow_omics_actions.R), the session-export
# collector and saver (utils_workflow_saved_exports.R), and the bindings' `composite_figure` and `final_pca`. The RUV
# kernels are the MCI lanes' stubs; every other kernel is real. Today's module copies stay until W4/W5 delete them;
# then the "today" sides of these blocks are quarantined with them.
# KIT-105 W4 deleted the lipid copies: the lipid halves are quarantined (Q-86 to Q-89); the metab halves stay live until W5.
# KIT-105 W5 deleted the metab copies: the metab halves are quarantined (Q-108 to Q-111).

normCoreSpec <- function(omic) {
    switch(omic,
        lipid = list(omic_type = "lipidomics", binding = lipidomicsDeclarationBinding, object = module_ci_lipid_norm_object,
            workflow = module_ci_lipid_norm_workflow, data = module_ci_lipid_norm_data, ruv_results = module_ci_lipid_norm_ruv_results,
            qc_key = "lipid_qc_dir", ruv_kernel = "runLipidPerAssayRuvOptimization", qc_kernel = "generateLipidQcPlots",
            complete = "lipid_norm_complete"),
        metab = list(omic_type = "metabolomics", binding = metabolomicsDeclarationBinding, object = module_ci_metab_norm_object,
            workflow = module_ci_metab_norm_workflow, data = module_ci_metab_norm_data, ruv_results = module_ci_metab_norm_ruv_results,
            qc_key = "metabolite_qc_dir", ruv_kernel = "runPerAssayRuvOptimization", qc_kernel = "generateMetabQcPlots",
            complete = "metab_norm_complete"))
}

normCoreAssays <- c("LCMS_Pos", "GCMS")
normCoreTime <- as.POSIXct("2026-11-08 12:00:00", tz = "UTC")

#' Small distinct PNGs at the composite's file names: `<assay>_<stage>_<kind>.png`.
normCoreWritePngs <- function(dir, assays = normCoreAssays) {
    i <- 0L
    for (assay in assays) for (stage in c("pre_norm", "post_norm", "ruv_corrected")) for (kind in c("pca", "density", "rle", "correlation")) {
        i <- i + 1L
        grDevices::png(file.path(dir, sprintf("%s_%s_%s.png", gsub("[^A-Za-z0-9]", "_", tolower(assay)), stage, kind)), width = 60, height = 40)
        graphics::par(mar = c(0, 0, 0, 0))
        graphics::plot.new()
        graphics::rect(0, 0, i / 30, 1, col = "grey40")
        grDevices::dev.off()
    }
    invisible(dir)
}

#' Run one composite writer into a fresh output directory: its value, log, the annotation labels it drew (in order),
#' the files it wrote and the bytes of its PNG.
normCoreComposite <- function(write) {
    out <- tempfile("kit105-w2-composite-")
    dir.create(out)
    on.exit(unlink(out, recursive = TRUE), add = TRUE)
    acc <- new.env(parent = emptyenv())
    acc$log <- character()
    acc$labels <- character()
    annotate <- ggplot2::annotate
    testthat::local_mocked_bindings(annotate = function(geom, x = NULL, y = NULL, ..., label = NULL, size = NULL) {
        if (!is.null(label)) acc$labels <- c(acc$labels, paste0(paste(label, collapse = " "), " [", size %||% "", "]"))
        annotate(geom, x = x, y = y, ..., label = label, size = size)
    }, .package = "ggplot2")
    value <- suppressMessages(write(out, function(entry) acc$log <- c(acc$log, entry)))
    files <- sort(list.files(out))
    png <- file.path(out, files[endsWith(files, ".png")])
    list(value = value, log = acc$log, labels = acc$labels, files = files,
        png = if (length(png)) unname(tools::md5sum(png)))
}

#' A plot's PNG bytes (its digest), at a fixed size.
normCorePlotPng <- function(plot) {
    path <- tempfile(fileext = ".png")
    on.exit(unlink(path), add = TRUE)
    suppressMessages(ggplot2::ggsave(path, plot, width = 6, height = 8, dpi = 72))
    unname(tools::md5sum(path))
}

#' The step provenance entries without their build identity and time (W1's reading; a stored product's entity id
#' differs between owners, so it is read as `<ENTITY>`).
normCoreProvenance <- function(manager) {
    lapply(manager$.context$snapshot$step_provenance, function(entry) {
        entry <- entry[intersect(c("operation", "variant", "parameter_names", "fixed", "values", "preceded_by"), names(entry))]
        rapply(entry, function(value) sub("^[0-9a-f]{64}$", "<ENTITY>", value), classes = "character", how = "replace")
    })
}

#' The committed state (or states) a `normData` slot's reference names.
normCoreSlot <- function(value, manager) {
    if (is.null(value)) return("<NULL>")
    if (!isSavedWorkflowReference(value)) return(paste0("<", class(value)[[1L]], ">"))
    states <- manager$.context$snapshot$states
    named <- names(Filter(function(ref) identical(ref$entity_id, value$entity_id), states))
    if (length(named)) paste(named, collapse = "|") else "<unlisted reference>"
}

#' The normalisation run's inputs for a scenario, by today's input names.
normCoreInput <- function(scenario) {
    list(apply_itsd = !identical(scenario, "skip"), itsd_aggregation = "median", log_offset = 1,
        norm_method = if (identical(scenario, "manual")) "quantile" else "cyclicloess",
        ruv_mode = switch(scenario, skip = "skip", manual = "manual", "automatic"),
        auto_percentage_min = 5, auto_percentage_max = 10, ruv_grouping_variable = "tech_rep_group",
        separation_metric = "max_difference", k_penalty_weight = 0.5, adaptive_k_penalty = TRUE, ruv_k = 2L, ruv_percentage = 5,
        min_pearson_correlation_threshold = 0.5)
}

#' Run the in-process normalisation once, today's (`handleLipidNormRunNormalization`, metab's observer wrapper) or the
#' core's (the frame's part as today's: the start line, then the core under the progress wrapper), and record what it
#' left: history, state digests, provenance, `normData`, the log, notifications, progress and the QC files.
#' Scenarios: skip (ITSD off, cyclicloess, RUV skip), manual (a manual ITSD selection, quantile, RUV manual k 2 / 5%),
#' automatic (ITSD by pattern, cyclicloess, RUV automatic, one post-normalisation figure removed before the composite),
#' error (as automatic, the RUV optimisation refused).
normCoreRun <- function(omic, scenario, core) {
    spec <- normCoreSpec(omic)
    object <- spec$object(positive_only = TRUE)
    root <- tempfile(paste0("kit105-w2-", omic, "-"))
    dir.create(root)
    on.exit(unlink(root, recursive = TRUE), add = TRUE)
    paths <- stats::setNames(list(file.path(root, "qc"), file.path(root, "source"), file.path(root, "export")),
        c(spec$qc_key, "source_dir", "export_dir"))
    workflow <- spec$workflow(object)
    norm_data <- spec$data(object)
    norm_data$post_filter_obj <- NULL
    norm_data$itsd_selections <- if (identical(scenario, "manual")) list(LCMS_Pos = 2L, GCMS = 1L) else list()
    rec <- new.env(parent = emptyenv())
    for (field in c("log", "notes", "progress")) rec[[field]] <- character()
    ruv_results <- spec$ruv_results
    stubs <- list(function(theObject, ruv_mode, params, experiment_paths) {
        if (identical(scenario, "error")) stop("kit105 w2: the RUV optimisation is refused")
        ruv_results(theObject, mode = ruv_mode)
    }, function(theObject, ruv_grouping_variable, ruv_number_k, ctrl) {
        theObject@args$ruvIII_C_Varying <- list(ruv_grouping_variable = ruv_grouping_variable, ruv_number_k = ruv_number_k, ctrl = ctrl)
        theObject
    })
    names(stubs) <- c(spec$ruv_kernel, "ruvIII_C_Varying")
    qc_kernel <- get(spec$qc_kernel, envir = asNamespace("MultiScholaR"))
    stubs[[spec$qc_kernel]] <- function(...) {
        result <- qc_kernel(...)
        if (scenario %in% c("automatic", "error") && identical(list(...)$stage, "post_norm")) {
            unlink(file.path(paths[[1L]], "lcms_pos_post_norm_density.png"))
        }
        result
    }
    do.call(testthat::local_mocked_bindings, c(stubs, list(.package = "MultiScholaR")))
    testthat::local_mocked_bindings(withProgress = function(expr, ...) expr,
        incProgress = function(amount = NULL, message = NULL, detail = NULL, ...) rec$progress <- c(rec$progress, paste(amount, detail)),
        showNotification = function(ui, ..., type = "default") {
            rec$notes <- c(rec$notes, paste0("[", type, "] ", ui))
            invisible("note")
        }, .package = "shiny")
    log <- function(entry) rec$log <- c(rec$log, entry)
    input <- normCoreInput(scenario)
    aesthetics <- list(color_var = "group", shape_var = "batch")
    withr::local_seed(20261108L)
    suppressWarnings(suppressMessages(if (!core) {
        if (identical(omic, "lipid")) {
            # Read only by the quarantined lipid block (Q-88): handleLipidNormRunNormalization was deleted by KIT-105 W4.
            handleLipidNormRunNormalization(input, workflow, paths, spec$omic_type, norm_data, log,
                getPlotAestheticsFn = function() aesthetics)
        } else {
            # Read only by the quarantined metab block (Q-110): runMetabNormNormalizationObserverWrapper was deleted by
            # KIT-105 W5.
            runMetabNormNormalizationObserverWrapper(workflow, input, paths, spec$omic_type, norm_data, addLogFn = log,
                getPlotAestheticsFn = function(colorVariable, shapeVariable) aesthetics)
        }
    } else {
        log("Starting normalization pipeline...")
        shiny::withProgress(message = "Running normalization pipeline...", value = 0, {
            runWorkflowOmicsInProcessNormalisation(spec$binding, workflow$state_manager, list(apply_itsd = input$apply_itsd,
                itsd_aggregation = input$itsd_aggregation, itsd_feature_ids = NULL, log_offset = input$log_offset,
                norm_method = input$norm_method, ruv_mode = input$ruv_mode, ruv_parameters = list(ruv_mode = input$ruv_mode,
                    percentage_min = input$auto_percentage_min, percentage_max = input$auto_percentage_max,
                    ruv_grouping_variable = input$ruv_grouping_variable, separation_metric = input$separation_metric,
                    k_penalty_weight = input$k_penalty_weight, adaptive_k_penalty = input$adaptive_k_penalty, manual_k = input$ruv_k,
                    manual_percentage = input$ruv_percentage),
                grouping_variable = aesthetics$color_var, shape_variable = aesthetics$shape_var, assay_names = norm_data$assay_names,
                qc_dir = paths[[1L]], experiment_paths = paths), norm_data, workflow$config_list, log,
                advance = function(detail) shiny::incProgress(1 / 6, detail = detail), notify = shiny::showNotification,
                finish = function() {
                    norm_data$plot_refresh_trigger <- norm_data$plot_refresh_trigger + 1
                    log("Normalization pipeline complete!")
                    shiny::showNotification("Normalization pipeline complete!", type = "message")
                }, itsd_table = omicGuiConfig(spec$omic_type)$norm$itsd_table)
        })
    }))
    manager <- workflow$state_manager
    nd <- as.list(norm_data)
    qc_files <- sort(list.files(paths[[1L]]))
    pngs <- file.path(paths[[1L]], qc_files[endsWith(qc_files, ".png")])
    list(history = as.character(manager$getHistory()), states = guiEntryStateDigests(manager),
        provenance = normCoreProvenance(manager),
        flags = vapply(c("normalization_complete", "ruv_complete", "correlation_filtering_complete"), function(flag) isTRUE(nd[[flag]]), TRUE),
        slots = vapply(c("post_filter_obj", "post_itsd_obj", "post_log2_obj", "post_norm_obj", "ruv_corrected_obj", "correlation_filtered_obj"),
            function(slot) normCoreSlot(nd[[slot]], manager), ""),
        ruv = nd$ruv_optimization_results, refresh = nd$plot_refresh_trigger,
        log = gsub(root, "<ROOT>", rec$log, fixed = TRUE), notes = gsub(root, "<ROOT>", rec$notes, fixed = TRUE),
        progress = rec$progress, qc_files = qc_files, png = stats::setNames(unname(tools::md5sum(pngs)), basename(pngs)))
}

test_that("the W2 core functions are defined in their core files, in the namespace", {
    files <- c(writeWorkflowOmicsNormComposite = "utils_workflow_omics_figures.R",
        omicsNormCompositeFromFiles = "utils_workflow_omics_figures.R",
        omicsNormFinalQcPcaPlot = "utils_workflow_omics_actions.R",
        runWorkflowOmicsInProcessNormalisation = "utils_workflow_omics_actions.R",
        collectWorkflowOmicsNormSession = "utils_workflow_saved_exports.R",
        saveWorkflowOmicsNormSessionArtifacts = "utils_workflow_saved_exports.R")
    namespace <- asNamespace("MultiScholaR")
    for (name in names(files)) {
        fn <- get0(name, envir = namespace, inherits = FALSE)
        expect_true(is.function(fn), info = name)
        if (!is.function(fn)) next
        expect_identical(environment(fn), namespace, info = name)
        # Source references are kept under load_all only; the check layout's installed package has none.
        if (!is.null(utils::getSrcref(fn))) expect_identical(basename(utils::getSrcFilename(fn, full.names = TRUE)), files[[name]], info = name)
    }
})

# KIT-105 W4: each comparison block that covered lipid and metab is split into a lipid block and a metab block, with the
# expectations unchanged; each expectation that reads only the core stays in the live (metab) block. The lipid blocks'
# "today" side called the lipid module functions W4 deleted, so they are quarantined (Q-86 to Q-89,
# docs/qa/workflow-execution/kit/quarantined-tests.md); the metab blocks stay live until W5. KIT-105 W5 deleted the metab
# module functions, so the metab blocks are quarantined too (Q-108 to Q-111).

test_that("writeWorkflowOmicsNormComposite is lipid's composite step: files, labels, log and returned fields (lipid)", {
    testthat::skip("quarantined: Q-86 (KIT-105 W4) it compares the core with writeLipidNormCompositeFigure and buildLipidNormCompositeFromFilesGenerator, deleted with mod_lipid_norm_{workflow,support}_helpers.R; the lipid tab draws the core's writeWorkflowOmicsNormComposite and omicsNormCompositeFromFiles, covered by test-omic-norm-module.R (b) and (c) (the composite with a missing file) and the lipid_canonical norm step")
    qc <- withr::local_tempdir(pattern = "kit105-w2-qc-")
    normCoreWritePngs(qc)
    # A missing cell: the composite is still written, the cell blank (W1 (c)).
    unlink(file.path(qc, "lcms_pos_post_norm_density.png"))
    for (mode in c("skip", "automatic")) {
        lipid <- normCoreComposite(function(out, log) writeLipidNormCompositeFigure(qc, normCoreAssays, mode, "lipidomics",
            buildLipidNormCompositeFromFilesGenerator(), output_dir = out, addLog = log))
        core_lipid <- normCoreComposite(function(out, log) writeWorkflowOmicsNormComposite(qc, normCoreAssays, mode, "lipidomics",
            output_dir = out, addLog = log))
        expect_identical(core_lipid[c("log", "labels", "files", "png")], lipid[c("log", "labels", "files", "png")], info = mode)
        expect_true(isTRUE(lipid$value$compositeSaved))
        expect_identical(core_lipid$value[names(lipid$value)], lipid$value, info = mode)
        expect_identical(lipid$files, c("lipidomics_composite_QC_figure.pdf", "lipidomics_composite_QC_figure.png",
            "lipidomics_composite_QC_figure.rds"))
    }
    # No QC directory: nothing is drawn and nothing is written.
    lipid <- normCoreComposite(function(out, log) writeLipidNormCompositeFigure(NULL, normCoreAssays, "skip", "lipidomics",
        buildLipidNormCompositeFromFilesGenerator(), output_dir = out, addLog = log))
    core <- normCoreComposite(function(out, log) writeWorkflowOmicsNormComposite(NULL, normCoreAssays, "skip", "lipidomics",
        output_dir = out, addLog = log))
    expect_identical(core$value[names(lipid$value)], lipid$value)
    expect_identical(core[c("log", "files")], lipid[c("log", "files")])
    # The figure's write failing: the warning in the log and the message returned.
    refuse <- function(...) stop("kit105 w2: the figure is refused")
    lipid <- normCoreComposite(function(out, log) writeLipidNormCompositeFigure(qc, normCoreAssays, "skip", "lipidomics",
        buildLipidNormCompositeFromFilesGenerator(), output_dir = out, addLog = log, savePlotFn = refuse, logWarnFn = function(...) NULL))
    core <- normCoreComposite(function(out, log) writeWorkflowOmicsNormComposite(qc, normCoreAssays, "skip", "metabolomics",
        output_dir = out, addLog = log, savePlotFn = refuse, logWarnFn = function(...) NULL))
    expect_identical(core$value[names(lipid$value)], lipid$value)
    expect_identical(core$log, lipid$log)
})

test_that("writeWorkflowOmicsNormComposite is metab's composite step: files, labels, log and returned fields (metab)", {
    testthat::skip("quarantined: Q-108 (KIT-105 W5) it compares the core with runMetabNormCompositeQcFigureStep and generateMetabNormCompositeFromFiles, deleted with mod_metab_norm_{pipeline,render}_helpers.R; the metab tab and its binding draw the core's writeWorkflowOmicsNormComposite and omicsNormCompositeFromFiles, covered by test-omic-norm-module.R (b) and (c) (the composite with a missing file), MCI-015 and the metab_lc, metab_gc and metab_combined norm steps")
    qc <- withr::local_tempdir(pattern = "kit105-w2-qc-")
    normCoreWritePngs(qc)
    # A missing cell: the composite is still written, the cell blank (W1 (c)).
    unlink(file.path(qc, "lcms_pos_post_norm_density.png"))
    for (mode in c("skip", "automatic")) {
        metab <- normCoreComposite(function(out, log) runMetabNormCompositeQcFigureStep(list(metabolite_qc_dir = qc), normCoreAssays, mode,
            "metabolomics", addLogFn = log, generateCompositeFromFilesFn = generateMetabNormCompositeFromFiles, outputDir = out))
        core_lipid <- normCoreComposite(function(out, log) writeWorkflowOmicsNormComposite(qc, normCoreAssays, mode, "lipidomics",
            output_dir = out, addLog = log))
        core_metab <- normCoreComposite(function(out, log) writeWorkflowOmicsNormComposite(qc, normCoreAssays, mode, "metabolomics",
            output_dir = out, addLog = log))
        expect_identical(core_metab[c("log", "labels", "files", "png")], metab[c("log", "labels", "files", "png")], info = mode)
        # Metab's extra fields: the core returns the superset (D7).
        expect_identical(core_metab$value, metab$value, info = mode)
        expect_identical(core_lipid$value[setdiff(names(core_lipid$value), "qcDir")], core_metab$value[setdiff(names(core_metab$value), "qcDir")])
    }
    # No QC directory: nothing is drawn and nothing is written.
    metab <- normCoreComposite(function(out, log) runMetabNormCompositeQcFigureStep(list(metabolite_qc_dir = NULL), normCoreAssays, "skip",
        "metabolomics", addLogFn = log, generateCompositeFromFilesFn = generateMetabNormCompositeFromFiles, outputDir = out))
    core <- normCoreComposite(function(out, log) writeWorkflowOmicsNormComposite(NULL, normCoreAssays, "skip", "lipidomics",
        output_dir = out, addLog = log))
    expect_identical(core$value, metab$value)
    # The figure's write failing: the warning in the log and the message returned.
    refuse <- function(...) stop("kit105 w2: the figure is refused")
    metab <- normCoreComposite(function(out, log) runMetabNormCompositeQcFigureStep(list(metabolite_qc_dir = qc), normCoreAssays, "skip",
        "metabolomics", addLogFn = log, generateCompositeFromFilesFn = generateMetabNormCompositeFromFiles, savePlotFn = refuse,
        logWarnFn = function(...) NULL, outputDir = out))
    core <- normCoreComposite(function(out, log) writeWorkflowOmicsNormComposite(qc, normCoreAssays, "skip", "metabolomics",
        output_dir = out, addLog = log, savePlotFn = refuse, logWarnFn = function(...) NULL))
    expect_identical(core$value, metab$value)
    expect_identical(core$log[[2L]], "Warning: Could not generate composite QC figure: kit105 w2: the figure is refused")
})

#' The composite-from-files comparison's inputs: one assay's files (one absent), its row labels and column titles, and a
#' drawing that records the annotation labels.
normCoreCompositeFromFilesCase <- function(qc) {
    normCoreWritePngs(qc, "LCMS_Pos")
    files <- file.path(qc, sprintf("lcms_pos_%s_%s.png", c("pre_norm", "post_norm", "ruv_corrected"), rep(c("pca", "density"), each = 3)))
    files[[2L]] <- file.path(qc, "absent.png")
    rows <- list(lcms_pos_pca = c("a)", "b)", "c)"), lcms_pos_density = c("d)", "e)", "f)"))
    columns <- c("Pre-Normalisation", "Post-Normalisation", "RUV-Corrected")
    function(generate) {
        labels <- character()
        annotate <- ggplot2::annotate
        testthat::local_mocked_bindings(annotate = function(geom, x = NULL, y = NULL, ..., label = NULL, size = NULL) {
            labels <<- c(labels, paste0(label, " [", size, "]"))
            annotate(geom, x = x, y = y, ..., label = label, size = size)
        }, .package = "ggplot2")
        result <- suppressMessages(generate(files, 3, rows, columns))
        list(size = c(result$width, result$height), labels = labels, png = normCorePlotPng(result$plot))
    }
}

test_that("omicsNormCompositeFromFiles is lipid's composite from files: the plot, its size and its labels (lipid)", {
    testthat::skip("quarantined: Q-86 (KIT-105 W4) it compares the core with buildLipidNormCompositeFromFilesGenerator, deleted with mod_lipid_norm_support_helpers.R; the lipid tab and its binding draw the core's omicsNormCompositeFromFiles, covered by the metab block below (the same core), test-workflow-omics-figures.R (KIT-079) and the lipid_canonical norm step")
    draw <- normCoreCompositeFromFilesCase(withr::local_tempdir(pattern = "kit105-w2-qc-"))
    core <- draw(omicsNormCompositeFromFiles)
    expect_identical(core, draw(buildLipidNormCompositeFromFilesGenerator()))
})

test_that("omicsNormCompositeFromFiles is metab's composite from files: the plot, its size and its labels (metab)", {
    testthat::skip("quarantined: Q-108 (KIT-105 W5) it compares the core with generateMetabNormCompositeFromFiles, deleted with mod_metab_norm_render_helpers.R; the metab tab and its binding draw the core's omicsNormCompositeFromFiles, covered by test-omic-norm-module.R (b) and (c) (the composite with a missing file), test-workflow-omics-figures.R (KIT-079), test-omics-audit-carry-overs.R (BUG-16) and the metab_lc, metab_gc and metab_combined norm steps")
    qc <- withr::local_tempdir(pattern = "kit105-w2-qc-")
    draw <- normCoreCompositeFromFilesCase(qc)
    core <- draw(omicsNormCompositeFromFiles)
    expect_identical(core, draw(generateMetabNormCompositeFromFiles))
    expect_identical(core$size, c(13, 14))
    # No file and no column titles: no composite, with the same warning.
    none <- file.path(qc, c("x.png", "y.png", "z.png"))
    expect_warning(expect_null(suppressMessages(omicsNormCompositeFromFiles(none, 3, list(r = c("a)", "b)", "c)"))))),
        "No valid plot sections to combine")
    expect_warning(expect_null(suppressMessages(generateMetabNormCompositeFromFiles(none, 3, list(r = c("a)", "b)", "c)"))))),
        "No valid plot sections to combine")
})

test_that("omicsNormFinalQcPcaPlot is lipid's final-QC PCA builder (lipid)", {
    testthat::skip("quarantined: Q-87 (KIT-105 W4) it compares the core with buildLipidNormFinalQcPcaPlot, deleted with mod_lipid_norm_support_helpers.R; the lipid tab and its binding draw the core's omicsNormFinalQcPcaPlot, covered by the metab block below (the core through both bindings), test-omic-norm-module.R (b) and (c) (the final QC plot and the final PCA failing)")
    for (omic in c("lipid", "metab")) {
        spec <- normCoreSpec(omic)
        for (layout in c("combined", "lc")) {
            # The PCA of a normalised state (the final QC plot's); the raw fixture has no PCA (W1's observation).
            object <- suppressMessages(suppressWarnings(logTransformAssays(spec$object(layout = layout, positive_only = TRUE), offset = 1)))
            core <-suppressMessages(suppressWarnings(omicsNormFinalQcPcaPlot(object, "group", "batch")))
            lipid <- suppressMessages(suppressWarnings(buildLipidNormFinalQcPcaPlot(object, "group", "batch")))
            expect_identical(class(core), class(lipid), info = paste(omic, layout))
            if (inherits(core, "patchwork")) {
                expect_identical(normCorePlotPng(core), normCorePlotPng(lipid), info = paste(omic, layout))
            } else {
                expect_identical(core$data, lipid$data, info = paste(omic, layout))
                expect_identical(core$labels, lipid$labels, info = paste(omic, layout))
            }
        }
    }
    # The PCA failing (D7): lipid's builder raises.
    testthat::local_mocked_bindings(plotPca = function(...) stop("kit105 w2: the PCA is refused"), .package = "MultiScholaR")
    object <- module_ci_metab_norm_object(positive_only = TRUE)
    expect_error(buildLipidNormFinalQcPcaPlot(object, "group", "batch"), "kit105 w2: the PCA is refused")
})

test_that("omicsNormFinalQcPcaPlot is metab's final-QC PCA builder; metab's binding keeps its error placeholder (metab)", {
    testthat::skip("quarantined: Q-109 (KIT-105 W5) it compares the core with buildMetabNormFinalQcPcaPlot, deleted with mod_metab_norm_session_helpers.R; the metab tab and its binding draw the core's omicsNormFinalQcPcaPlot (the binding keeping metab's placeholder), covered by test-omic-norm-module.R (b) and (c) (the final QC plot and the final PCA failing) and the metab_lc, metab_gc and metab_combined norm steps")
    for (omic in c("lipid", "metab")) {
        spec <- normCoreSpec(omic)
        for (layout in c("combined", "lc")) {
            # The PCA of a normalised state (the final QC plot's); the raw fixture has no PCA (W1's observation).
            object <- suppressMessages(suppressWarnings(logTransformAssays(spec$object(layout = layout, positive_only = TRUE), offset = 1)))
            core <-suppressMessages(suppressWarnings(omicsNormFinalQcPcaPlot(object, "group", "batch")))
            metab <- suppressMessages(suppressWarnings(buildMetabNormFinalQcPcaPlot(object, "group", "batch")))
            if (inherits(core, "patchwork")) {
                expect_identical(normCorePlotPng(core), normCorePlotPng(metab), info = paste(omic, layout))
            } else {
                expect_identical(core$data, metab$data, info = paste(omic, layout))
            }
            # Through the bindings, the managed path's builder (the worker's `final_pca`).
            through <- suppressMessages(suppressWarnings(spec$binding$final_pca(object, "group", "batch")))
            expect_identical(normCorePlotPng(through), normCorePlotPng(core), info = paste(omic, layout))
        }
    }
    expect_true("omicsNormFinalQcPcaPlot" %in% all.names(body(lipidomicsDeclarationBinding$final_pca)))
    expect_true("omicsNormFinalQcPcaPlot" %in% all.names(body(metabolomicsDeclarationBinding$final_pca)))
    # The PCA failing (D7): the core raises; metab's binding returns metab's placeholder.
    testthat::local_mocked_bindings(plotPca = function(...) stop("kit105 w2: the PCA is refused"), .package = "MultiScholaR")
    object <- module_ci_metab_norm_object(positive_only = TRUE)
    expect_error(omicsNormFinalQcPcaPlot(object, "group", "batch"), "kit105 w2: the PCA is refused")
    expect_error(lipidomicsDeclarationBinding$final_pca(module_ci_lipid_norm_object(positive_only = TRUE), "group", "batch"),
        "kit105 w2: the PCA is refused")
    placeholder <- metabolomicsDeclarationBinding$final_pca(object, "group", "batch")
    today <- buildMetabNormFinalQcPcaPlot(object, "group", "batch")
    expect_identical(placeholder$layers[[1L]]$aes_params, today$layers[[1L]]$aes_params)
    expect_identical(placeholder$layers[[1L]]$aes_params$label, "Error: kit105 w2: the PCA is refused")
    expect_identical(placeholder$layers[[1L]]$aes_params$size, 4)
    expect_identical(normCorePlotPng(placeholder), normCorePlotPng(today))
})

test_that("the lipid binding's composite_figure gives today's product through the core (lipid)", {
    testthat::skip("quarantined: Q-86 (KIT-105 W4) it compares the binding's composite_figure with writeLipidNormCompositeFigure, deleted with mod_lipid_norm_workflow_helpers.R; the binding calls the core (asserted by the metab block below), covered by the lipid journeys' managed composite (test-workflow-omics-worker.R)")
    qc <- withr::local_tempdir(pattern = "kit105-w2-qc-")
    normCoreWritePngs(qc, "LCMS_Pos")
    for (mode in c("skip", "manual")) {
        lipid <- normCoreComposite(function(out, log) writeLipidNormCompositeFigure(qc, "LCMS_Pos", mode, "lipidomics",
            buildLipidNormCompositeFromFilesGenerator(), output_dir = out))
        through_lipid <- normCoreComposite(function(out, log) lipidomicsDeclarationBinding$composite_figure(qc, out, "LCMS_Pos", mode))
        # Lipid's step returned `compositeSaved` (and `errorMessage`), what the worker reads; the core returns the superset.
        expect_identical(through_lipid$value[names(lipid$value)], lipid$value, info = mode)
        expect_identical(through_lipid[c("files", "png", "labels")], lipid[c("files", "png", "labels")], info = mode)
    }
})

test_that("the metab binding's composite_figure gives today's product through the core (metab)", {
    testthat::skip("quarantined: Q-108 (KIT-105 W5) it compares the binding's composite_figure with runMetabNormCompositeQcFigureStep, deleted with mod_metab_norm_pipeline_helpers.R; the binding calls the core's writeWorkflowOmicsNormComposite, covered by the metab journeys' managed composite (test-workflow-omics-worker.R), covered by test-omic-norm-module.R (b) and (c) and the metab_lc, metab_gc and metab_combined norm steps")
    qc <- withr::local_tempdir(pattern = "kit105-w2-qc-")
    normCoreWritePngs(qc, "LCMS_Pos")
    for (mode in c("skip", "manual")) {
        metab <- normCoreComposite(function(out, log) runMetabNormCompositeQcFigureStep(list(metabolite_qc_dir = qc), "LCMS_Pos", mode,
            "metabolomics", generateCompositeFromFilesFn = generateMetabNormCompositeFromFiles, outputDir = out))
        through_metab <- normCoreComposite(function(out, log) metabolomicsDeclarationBinding$composite_figure(qc, out, "LCMS_Pos", mode))
        expect_identical(through_metab[c("value", "files", "png", "labels")], metab[c("value", "files", "png", "labels")], info = mode)
    }
    expect_true("writeWorkflowOmicsNormComposite" %in% all.names(body(lipidomicsDeclarationBinding$composite_figure)))
    expect_true("writeWorkflowOmicsNormComposite" %in% all.names(body(metabolomicsDeclarationBinding$composite_figure)))
})

test_that("runWorkflowOmicsInProcessNormalisation is lipid's in-process sequence: states, provenance, normData, log (lipid)", {
    testthat::skip("quarantined: Q-88 (KIT-105 W4) its \"today\" side is handleLipidNormRunNormalization, deleted with mod_lipid_norm_workflow_helpers.R; the lipid tab runs the core's runWorkflowOmicsInProcessNormalisation, covered by test-omic-norm-module.R (b) (the three runs), MCI-021.1 to .3 and the lipid_canonical norm step")
    runs <- list()
    for (scenario in c("skip", "manual", "automatic", "error")) {
        today <- normCoreRun("lipid", scenario, core = FALSE)
        core <- normCoreRun("lipid", scenario, core = TRUE)
        for (field in names(today)) expect_identical(core[[field]], today[[field]], info = paste("lipid", scenario, field))
        runs[[paste("lipid", scenario)]] <- core
    }
    # What the scenarios exercised (so that an empty run cannot pass).
    expect_true(any(grepl("ITSD features selected", runs[["lipid manual"]]$log)), info = "lipid")
    automatic <- runs[["lipid automatic"]]
    expect_identical(automatic$history, c("initial", "lipid_qc_complete", "lipid_itsd_norm", "lipid_log2", "lipid_normalized",
        "lipid_ruv_corrected"))
    expect_identical(automatic$notes, "[message] Normalization pipeline complete!")
    expect_false("lcms_pos_post_norm_density.png" %in% automatic$qc_files)
    expect_true("lipidomics_composite_QC_figure.png" %in% automatic$qc_files)
})

test_that("runWorkflowOmicsInProcessNormalisation is metab's in-process sequence: states, provenance, normData, log (metab)", {
    testthat::skip("quarantined: Q-110 (KIT-105 W5) its \"today\" side is runMetabNormNormalizationObserverWrapper, deleted with mod_metab_norm_binding_helpers.R (and its pipeline shells with mod_metab_norm_pipeline_helpers.R); the metab tab runs the core's runWorkflowOmicsInProcessNormalisation, covered by test-omic-norm-module.R (b) and (c) (the three runs), MCI-015.1 to .3 and the metab_lc, metab_gc and metab_combined norm steps")
    runs <- list()
    for (scenario in c("skip", "manual", "automatic", "error")) {
        today <- normCoreRun("metab", scenario, core = FALSE)
        core <- normCoreRun("metab", scenario, core = TRUE)
        for (field in names(today)) expect_identical(core[[field]], today[[field]], info = paste("metab", scenario, field))
        runs[[paste("metab", scenario)]] <- core
    }
    # What the scenarios exercised (so that an empty run cannot pass).
    expect_true(any(grepl("ITSD features selected", runs[["metab manual"]]$log)), info = "metab")
    expect_identical(runs[["metab skip"]]$slots[["ruv_corrected_obj"]], "metab_normalized")
    error <- runs[["metab error"]]
    expect_identical(error$notes[[1L]], "[error] Error: kit105 w2: the RUV optimisation is refused")
    expect_false(error$flags[["ruv_complete"]])
})

#' The normalisation export's fixture: the run's object committed under the omic's completion state, `normData` as a
#' finished run leaves it (RUV automatic results), and today's inputs.
normCoreExportFixture <- function(omic) {
    spec <- normCoreSpec(omic)
    object <- spec$object(positive_only = TRUE)
    workflow <- spec$workflow(object, current_state_name = spec$complete)
    norm_data <- spec$data(object)
    norm_data$normalization_complete <- TRUE
    norm_data$ruv_complete <- TRUE
    norm_data$correlation_filtering_complete <- TRUE
    norm_data$ruv_optimization_results <- spec$ruv_results(object, mode = "automatic")
    norm_data$correlation_results <- module_ci_lipid_norm_correlation_pairs("pass_all", normCoreAssays)
    list(spec = spec, workflow = workflow, norm_data = norm_data, input = normCoreInput("automatic"),
        export = omicGuiConfig(spec$omic_type)$norm$export)
}

#' One omic's collector comparison (`today`: the omic's builder): the in-process list, the counts field's place, the
#' samples count rule (P-18) and the managed list.
normCoreCollectCase <- function(omic, today) {
    fixture <- normCoreExportFixture(omic)
    core <- function() collectWorkflowOmicsNormSession(fixture$spec$binding, fixture$workflow, fixture$norm_data, fixture$input,
        "kit105", fixture$export, normCoreTime)
    in_process <- core()
    expect_identical(in_process, today(fixture), info = omic)
    expect_identical(names(in_process)[[21L]], fixture$export$counts_field)
    expect_true(isS4(in_process$current_s4_object))
    # P-18, per omic: lipid counts the design's samples among the columns, metab the columns less the id and
    # annotation columns.
    expect_identical(vapply(in_process$feature_counts, function(counts) as.integer(counts$samples), 1L),
        stats::setNames(if (identical(omic, "lipid")) c(4L, 4L) else c(5L, 5L), normCoreAssays), info = omic)
    # A managed session: the head's shell, and no state in the list.
    local({
        testthat::local_mocked_bindings(workflowOmicsManagedActions = function(omic_type) list(), .package = "MultiScholaR")
        managed <- core()
        expect_identical(managed, today(fixture), info = omic)
        expect_null(managed$current_s4_object)
        expect_true("current_s4_object" %in% names(managed))
    })
}

test_that("collectWorkflowOmicsNormSession is lipid's session list: fields, order, counts (P-18) (lipid)", {
    testthat::skip("quarantined: Q-89 (KIT-105 W4) it compares the core with buildLipidNormExportSessionData, deleted with mod_lipid_norm_workflow_helpers.R; the lipid tab collects through the core's collectWorkflowOmicsNormSession, covered by test-omic-norm-module.R (b) (the export's session-RDS fields and counts), MCI-021.5 and .6, test-workflow-omics-sessions.R and the lipid_canonical session_export step")
    normCoreCollectCase("lipid", function(fixture) {
        buildLipidNormExportSessionData(fixture$workflow, fixture$norm_data, fixture$input, "kit105", exportTimestamp = normCoreTime)
    })
})

test_that("collectWorkflowOmicsNormSession is metab's session list: fields, order, counts (P-18) (metab)", {
    testthat::skip("quarantined: Q-111 (KIT-105 W5) it compares the core with buildMetabNormExportSessionData, deleted with mod_metab_norm_session_helpers.R; the metab tab collects through the core's collectWorkflowOmicsNormSession, covered by test-omic-norm-module.R (b) (the export's session-RDS fields and counts), MCI-015.5 and .6, test-workflow-omics-sessions.R and the metab_lc, metab_gc and metab_combined norm steps")
    normCoreCollectCase("metab", function(fixture) {
        buildMetabNormExportSessionData(fixture$workflow, fixture$norm_data, fixture$input, "kit105", exportTimestamp = normCoreTime)
    })
})

#' One omic's saver comparison (`today(session, dir, log, progress)`: the omic's savers): the files, their bytes, the
#' summary, the log and the progress.
normCoreSaveCase <- function(omic, today) {
    request <- list(request_id = "dia-export-kit105w2", timestamp = "20261108_120000")
    fixture <- normCoreExportFixture(omic)
    session <- collectWorkflowOmicsNormSession(fixture$spec$binding, fixture$workflow, fixture$norm_data, fixture$input, "kit105",
        fixture$export, normCoreTime)
    run <- function(save) {
        dir <- tempfile(paste0("kit105-w2-export-", omic, "-"))
        dir.create(dir)
        on.exit(unlink(dir, recursive = TRUE), add = TRUE)
        acc <- new.env(parent = emptyenv())
        acc$log <- character()
        acc$progress <- character()
        log <- function(message) acc$log <- c(acc$log, gsub(dir, "<DIR>", message, fixed = TRUE))
        progress <- function(amount, detail) acc$progress <- c(acc$progress, paste(amount, detail))
        value <- save(dir, log, progress)
        # F105-03 (KIT-105 W5): the returned publication's completion time is the clock's, as the journal's is below.
        if (!is.null(value$publication$completed_at)) value$publication$completed_at <- "<TIME>"
        files <- sort(list.files(dir, recursive = TRUE, all.files = TRUE))
        rds <- files[grepl("\\.(rds|RDS)$", files)]
        text <- files[grepl("\\.(txt|json)$", files)]
        list(value = rapply(value, function(x) gsub(dir, "<DIR>", x, fixed = TRUE), classes = "character", how = "replace"),
            files = files, rds = lapply(stats::setNames(file.path(dir, rds), rds), readRDS),
            # The journal's completion time is the clock's.
            text = lapply(stats::setNames(file.path(dir, text), text), function(path) {
                gsub("\"completed_at\": *\"[^\"]*\"", "\"completed_at\":\"<TIME>\"", gsub(dir, "<DIR>", readLines(path, warn = FALSE), fixed = TRUE))
            }),
            log = acc$log, progress = acc$progress)
    }
    today_run <- run(function(dir, log, progress) today(session, dir, request, log, progress))
    core <- run(function(dir, log, progress) {
        saveWorkflowOmicsNormSessionArtifacts(fixture$spec$binding, session, dir, request, fixture$export,
            advance = if (identical(omic, "metab")) progress, summaryTime = normCoreTime, logInfoFn = log)
    })
    for (field in c("files", "rds", "text", "log", "progress")) expect_identical(core[[field]], today_run[[field]], info = paste(omic, field))
    expect_identical(core$value[names(today_run$value)], today_run$value, info = omic)
    summary <- core$text[[paste0(fixture$export$prefix, "filtered_session_summary.txt")]]
    expect_identical(summary[[1L]], paste(if (identical(omic, "lipid")) "Lipidomics" else "Metabolomics",
        "Normalized Session Data Export Summary"))
    expect_true(paste0(fixture$export$prefix, "ruv_optimization_results.RDS") %in% core$files)
}

test_that("saveWorkflowOmicsNormSessionArtifacts is lipid's saver: files, bytes, summary, log and progress (lipid)", {
    testthat::skip("quarantined: Q-89 (KIT-105 W4) it compares the core with saveLipidNormExportArtifacts, deleted with mod_lipid_norm_workflow_helpers.R; the lipid tab saves through the core's saveWorkflowOmicsNormSessionArtifacts, covered by test-omic-norm-module.R (b) (the export's files and summary text), MCI-021.5 and .6, test-workflow-omics-sessions.R and the lipid_canonical session_export step")
    normCoreSaveCase("lipid", function(session, dir, request, log, progress) saveLipidNormExportArtifacts(session, dir, request, logInfoFn = log))
})

test_that("saveWorkflowOmicsNormSessionArtifacts is metab's savers: files, bytes, summary, log and progress (metab)", {
    testthat::skip("quarantined: Q-111 (KIT-105 W5) it compares the core with saveMetabNormExportSessionRdsFiles, saveMetabNormExportMetadataFiles and saveMetabNormExportSummaryFile, deleted with mod_metab_norm_session_helpers.R; the metab tab saves through the core's saveWorkflowOmicsNormSessionArtifacts, covered by test-omic-norm-module.R (b) (the export's files and summary text), test-omic-norm-module.R's W5 (D7) export block, MCI-015.5, test-workflow-omics-sessions.R and the metab_lc, metab_gc and metab_combined norm steps")
    normCoreSaveCase("metab", function(session, dir, request, log, progress) {
        # runMetabNormExportSessionWorkflow's saveFn (mod_metab_norm_session_helpers.R:741-756), with the summary's
        # time fixed.
        files <- saveMetabNormExportSessionRdsFiles(sessionData = session, sourceDir = dir, logInfoFn = log,
            incProgressFn = function(amount, detail) progress(amount, detail), publicationRequest = request)
        progress(0.1, "Saving metadata files...")
        saveMetabNormExportMetadataFiles(sessionData = session, sourceDir = dir, logInfoFn = log)
        progress(0.2, "Creating summary...")
        saveMetabNormExportSummaryFile(sessionData = session, sourceDir = dir, sessionFilename = files$sessionFilename,
            timeFn = function() normCoreTime, logInfoFn = log)
        files
    })
    # Without a publication request the saver refuses, as both savers do.
    fixture <- normCoreExportFixture("lipid")
    expect_error(saveWorkflowOmicsNormSessionArtifacts(fixture$spec$binding, list(), tempdir(), NULL, fixture$export),
        class = "saved_entity_scope")
})
