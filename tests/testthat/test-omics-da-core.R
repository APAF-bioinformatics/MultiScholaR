# KIT-106 W2: the DA tab's core moves (D2, D3, D4), written before them.
#   - each moved function is in the package namespace and its source is its new file;
#   - `workflowOmicsInProcessVolcanoPageFactory` gives the managed factory's iframe from the binding's
#     `<prefix>_volcano_page`, run in process (lipid, metab and trx on their MCI fixtures), and NULL when
#     the widget reports nothing;
#   - `loadWorkflowOmicsDaSession` restores W1 (b)'s bus, tab values and choices, and gives its errors
#     (test-omic-da-module.R: lipid_main, lipid_mismatch, lipid_reordered, lipid_refusals and the metab twins);
#   - the D4 equivalence: `runWorkflowOmicsPublicationInProcess(..., "<prefix>_da_publication")` makes the export's
#     per-assay volcano and heatmap calls and writes the file set and tables of W1's lipid_main and metab_main run
#     step, with the "DA results saved to disk" notice. Its comparison side (the in-process writers of the start
#     commit, `outputLipidDaResultsAllContrasts` on `lipidDaResultsWithSource` and `writeMetabDaResultsArtifacts`
#     on `metabDaResultsWithSource`) is one block per omic, quarantined when W3 (lipid) and W4 (metab) delete those
#     mod_ helpers; the core-side block stays live.

omicsDaCoreSrcFile <- function(name) {
    fn <- get(name, envir = asNamespace("MultiScholaR"), inherits = FALSE)
    list(env = environmentName(environment(fn)), file = utils::getSrcFilename(fn))
}

test_that("(W2) the moved functions are in the namespace, each in its new file", {
    moved <- c(
        protDaScalarString = "utils_workflow_dia_da_validation.R",
        normaliseProtDaContrastsTable = "utils_workflow_dia_da_validation.R",
        protDaAnalysisContrastStrings = "utils_workflow_dia_da_validation.R",
        protDaContrastExpression = "utils_workflow_dia_da_validation.R",
        protDaContrastTerms = "utils_workflow_dia_da_validation.R",
        validateProtDaFilteredSession = "utils_workflow_dia_da_validation.R",
        validateProtDaModelAndContrasts = "utils_workflow_dia_da_validation.R",
        validateLipidDaAnalysisPreflight = "func_omics_design_shared.R",
        trxDaModelText = "utils_workflow_trx_da.R",
        trxDaReadMetadata = "utils_workflow_trx_da.R",
        trxDaContrastEntries = "utils_workflow_trx_da.R",
        trxDaContrastEntry = "utils_workflow_trx_da.R",
        trxDaStaticVolcano = "utils_workflow_trx_da.R",
        trxDaExport = "utils_workflow_trx_da.R",
        workflowOmicsInProcessVolcanoPageFactory = "utils_workflow_omics_da.R",
        loadWorkflowOmicsDaSession = "utils_workflow_omics_da.R")
    for (name in names(moved)) {
        expect_true(exists(name, envir = asNamespace("MultiScholaR"), inherits = FALSE), info = name)
        if (!exists(name, envir = asNamespace("MultiScholaR"), inherits = FALSE)) next
        where <- omicsDaCoreSrcFile(name)
        expect_identical(where$env, "MultiScholaR", info = name)
        expect_identical(where$file, moved[[name]], info = name)
    }
    expect_false(file.exists(testthat::test_path("..", "..", "R", "mod_prot_da_validation.R")))
})

# ------------------------------------------------------------------------------------------------
# Fixtures: W1 (b)'s lipid and metab fixtures, without the module.
# ------------------------------------------------------------------------------------------------

#' Lipid's or metab's fixture as W1 (b) makes it (test-omic-da-module.R `omicDaW1FeatureSetup`): the MCI object
#' on an owner of the omic (`state = FALSE`: no owner), the multi-group contrasts and the MCI legacy session
#' file (`session_file = FALSE`: none). The tab's values are an environment, as the MCI buses are.
omicsDaCoreFeature <- function(omic, state = TRUE, session_file = TRUE, session_assays = NULL, env = parent.frame()) {
    spec <- switch(omic,
        lipid = list(omic_type = "lipidomics", object = module_ci_lipid_da_object, contrasts = module_ci_lipid_da_contrasts,
            session = module_ci_lipid_da_session_data, state_name = "lipid_correlation_filtered",
            binding = lipidomicsDeclarationBinding),
        metab = list(omic_type = "metabolomics", object = module_ci_metab_da_object, contrasts = module_ci_metab_da_contrasts,
            session = module_ci_metab_da_session_data, state_name = "metab_correlation_filtered",
            binding = metabolomicsDeclarationBinding))
    root <- normalizePath(withr::local_tempdir(pattern = paste0("kit106w2-", omic, "-"), .local_envir = env))
    paths <- list(source_dir = file.path(root, "source"), da_output_dir = file.path(root, "da_output"),
        publication_graphs_dir = file.path(root, "publication_graphs"))
    for (path in paths) dir.create(path, recursive = TRUE)
    object <- spec$object()
    bus <- list2env(list(contrasts_tbl = spec$contrasts("multi_group"), design_matrix = object@design_matrix,
        config_list = list(da = list(mode = "module-ci"))), parent = emptyenv())
    if (state) {
        bus <- omicsTestOwner(object, spec$state_name, workflow = bus, omic_type = spec$omic_type, experiment_paths = paths,
            env = env)
    }
    session_load <- omicGuiConfig(spec$omic_type)$da$session_load
    # The assay check's order (D7): metab's before the adoption, lipid's after the contrasts. Reported to the lead
    # as the `da$session_load$assay_check_first` value to wire.
    session_load$assay_check_first <- identical(omic, "metab")
    file <- file.path(paths$source_dir, session_load$file)
    if (session_file) saveRDS(spec$session(assay_names = session_assays), file)
    tab <- new.env(parent = emptyenv())
    list(omic = omic, binding = spec$binding, bus = bus, tab = tab, paths = paths, root = root, file = file,
        session_load = session_load, state_name = spec$state_name)
}

#' The session load as the frame calls it, in the module's session (the adopted owner attaches to it), with its
#' "Loading" notification recorded.
omicsDaCoreLoad <- function(fixture) {
    notes <- character()
    out <- suppressMessages(shiny::withReactiveDomain(omicsOwnerSession(), loadWorkflowOmicsDaSession(fixture$binding,
        fixture$bus, fixture$tab, fixture$paths, fixture$session_load, reading = function() notes <<- c(notes, "reading"))))
    out$notes <- notes
    out
}

#' The tab's values as W1 (b) describes them (`omicDaW1TabValues`): a vector by its values, a data frame by its shape.
omicsDaCoreTab <- function(tab) {
    describe <- function(value) {
        if (is.null(value)) return("<NULL>")
        if (is.data.frame(value)) return(sprintf("<data.frame %dx%d>", nrow(value), ncol(value)))
        paste(as.character(value), collapse = "|")
    }
    vapply(c("assays_available", "contrasts_available", "current_state_name", "formula_from_s4"),
        function(name) describe(tab[[name]]), "")
}

test_that("(W2, D4) the session load restores W1 (b)'s bus, tab values and choices for lipid and metab", {
    for (omic in c("lipid", "metab")) {
        fixture <- omicsDaCoreFeature(omic)
        out <- omicsDaCoreLoad(fixture)
        expect_identical(out$status, "success", info = omic)
        expect_null(out$message, info = omic)
        expect_identical(out$notes, "reading", info = omic)
        expect_identical(out$contrast_choices, c("KO_vs_WT", "RES_vs_WT"), info = omic)
        expect_identical(out$assay_choices, c("Combined", "LCMS_Pos", "GCMS"), info = omic)
        expect_identical(out$table_assay_choices, c("All", "LCMS_Pos", "GCMS"), info = omic)
        expect_identical(out$formula, "~ 0 + group", info = omic)
        expect_identical(omicsDaCoreTab(fixture$tab), c(assays_available = "LCMS_Pos|GCMS",
            contrasts_available = "<data.frame 2x3>", current_state_name = fixture$state_name,
            formula_from_s4 = "~ 0 + group"), info = omic)
        expect_identical(fixture$bus$contrasts_tbl, fixture$tab$contrasts_available, info = omic)
        # The legacy session is adopted into a fresh owner on the bus (W1: the `.dia-native-legacy-` store).
        expect_identical(fixture$bus$state_manager$current_state, fixture$state_name, info = omic)
        expect_true(any(grepl(".dia-native-legacy-", list.files(fixture$paths$source_dir, all.files = TRUE), fixed = TRUE)),
            info = omic)
    }
})

test_that("(W2, D4, D7) the session load's refusals: a missing file, and the assay check's rule, text and order", {
    for (omic in c("lipid", "metab")) {
        # W1 lipid_refusals / metab_refusals load_missing: no owner and no file; nothing is read.
        missing <- omicsDaCoreFeature(omic, state = FALSE, session_file = FALSE)
        out <- omicsDaCoreLoad(missing)
        expect_identical(out$status, "missing", info = omic)
        expect_identical(out$message, sprintf("Session file not found: %s", missing$file), info = omic)
        expect_identical(out$notes, character(), info = omic)
        expect_identical(unname(omicsDaCoreTab(missing$tab)), rep("<NULL>", 4L), info = omic)
    }
    # A missing source and export directory.
    nowhere <- omicsDaCoreFeature("lipid", state = FALSE, session_file = FALSE)
    nowhere$paths$source_dir <- file.path(nowhere$root, "absent")
    expect_identical(omicsDaCoreLoad(nowhere)[c("status", "message")],
        list(status = "missing", message = "Could not find source directory for session data."))

    # W1 lipid_mismatch: lipid adopts the session and restores its contrasts before the check refuses.
    lipid <- omicsDaCoreFeature("lipid", session_assays = c("LCMS_Pos", "LCMS_Neg"))
    out <- omicsDaCoreLoad(lipid)
    expect_identical(out$status, "error")
    expect_identical(out$message, paste("Loaded lipid session assay_names do not match S4 lipid_data assays:",
        "session [LCMS_Pos, LCMS_Neg], S4 [LCMS_Pos, GCMS]"))
    expect_identical(out$contrast_choices, c("KO_vs_WT", "RES_vs_WT"))
    expect_null(out$assay_choices)
    expect_null(out$formula)
    expect_identical(omicsDaCoreTab(lipid$tab), c(assays_available = "<NULL>", contrasts_available = "<data.frame 2x3>",
        current_state_name = "lipid_correlation_filtered", formula_from_s4 = "<NULL>"))
    # W1 lipid_reordered: a set comparison, so the same names in another order load.
    reordered <- omicsDaCoreFeature("lipid", session_assays = c("GCMS", "LCMS_Pos"))
    out <- omicsDaCoreLoad(reordered)
    expect_identical(out$status, "success")
    expect_identical(out$assay_choices, c("Combined", "GCMS", "LCMS_Pos"))
    expect_identical(out$table_assay_choices, c("All", "GCMS", "LCMS_Pos"))

    # W1 metab_mismatch and metab_reordered: metab compares in order and refuses before the adoption.
    for (assays in list(c("LCMS_Pos", "LCMS_Neg"), c("GCMS", "LCMS_Pos"))) {
        metab <- omicsDaCoreFeature("metab", session_assays = assays)
        head <- metab$bus$state_manager
        out <- omicsDaCoreLoad(metab)
        expect_identical(out$status, "error")
        expect_identical(out$message, sprintf("Session assay_names do not match current_s4_object assays: session=%s; object=%s",
            paste(assays, collapse = ", "), "LCMS_Pos, GCMS"))
        expect_null(out$contrast_choices)
        expect_identical(unname(omicsDaCoreTab(metab$tab)), rep("<NULL>", 4L))
        expect_identical(metab$bus$state_manager, head)
    }
})

# ------------------------------------------------------------------------------------------------
# The in-process volcano page (D3).
# ------------------------------------------------------------------------------------------------

#' The managed factory's iframe (utils_workflow_omics_actions.R `workflowOmicsVolcanoPageFactory`) of a saved page.
omicsDaCoreIframe <- function(manager, page, contrast) {
    shiny::tags$iframe(srcdoc = rawToChar(hydrateSavedWorkflowEntity(manager$.context$store, page)$bytes),
        style = "width:100%;height:820px;border:0;", title = paste("Interactive volcano", contrast))
}

#' The page the factory read: the `<prefix>_volcano_page` display, recorded as the runner returns it.
omicsDaCoreRecordPages <- function(env = parent.frame()) {
    seen <- new.env(parent = emptyenv())
    seen$displays <- list()
    original <- runWorkflowOmicsPublicationInProcess
    testthat::local_mocked_bindings(runWorkflowOmicsPublicationInProcess = function(binding, manager, operation, parameters,
        experiment_paths) {
        display <- original(binding, manager, operation, parameters, experiment_paths)
        seen$displays[[length(seen$displays) + 1L]] <- list(operation = operation, parameters = parameters, display = display)
        display
    }, .env = env)
    seen
}

#' Lipid's or metab's DA in process on the loaded session, as the tab runs it (W1 main: the load, then the run).
omicsDaCoreRunFeature <- function(fixture) {
    expect_identical(omicsDaCoreLoad(fixture)$status, "success")
    manager <- fixture$bus$state_manager
    run <- if (identical(fixture$omic, "lipid")) runLipidReferenceDaAnalysis else runMetabReferenceDaAnalysis
    source <- getWorkflowStateReference(manager, fixture$tab$current_state_name)
    suppressWarnings(suppressMessages(run(manager, source, list(formula_string = "~ 0 + group", da_q_val_thresh = 0.05,
        treat_lfc_cutoff = 0), fixture$bus$contrasts_tbl)$results))
}

test_that("(W2, D3) the in-process page factory gives the managed factory's iframe for lipid and metab, NULL for no widget", {
    testthat::skip_if_not_installed("Glimma")
    for (omic in c("lipid", "metab")) {
        fixture <- omicsDaCoreFeature(omic)
        results <- omicsDaCoreRunFeature(fixture)
        manager <- fixture$bus$state_manager
        seen <- omicsDaCoreRecordPages()
        factory <- workflowOmicsInProcessVolcanoPageFactory(fixture$binding, manager, fixture$paths)
        page <- suppressWarnings(suppressMessages(factory(results, "KO_vs_WT", "LCMS_Pos", 0.05)))
        expect_length(seen$displays, 1L)
        expect_identical(seen$displays[[1L]]$operation, paste0(fixture$binding$prefix, "_volcano_page"), info = omic)
        expect_identical(seen$displays[[1L]]$parameters, list(contrast = "KO_vs_WT", assay = "LCMS_Pos", qvalue = 0.05),
            info = omic)
        reference <- seen$displays[[1L]]$display$page
        expect_false(is.null(reference), info = omic)
        expect_identical(page, omicsDaCoreIframe(manager, reference, "KO_vs_WT"), info = omic)
        expect_match(page$attribs$srcdoc, "<html", fixed = TRUE, info = omic)
    }
    # The widget reports nothing: no page, so the tab says it could not draw one.
    fixture <- omicsDaCoreFeature("lipid")
    results <- omicsDaCoreRunFeature(fixture)
    testthat::local_mocked_bindings(generateLipidDAVolcanoPlotGlimma = function(...) NULL)
    factory <- workflowOmicsInProcessVolcanoPageFactory(fixture$binding, fixture$bus$state_manager, fixture$paths)
    expect_null(suppressMessages(factory(results, "KO_vs_WT", "LCMS_Pos", 0.05)))
    # The worker's message is raised, as the managed factory raises it.
    testthat::local_mocked_bindings(generateLipidDAVolcanoPlotGlimma = function(...) stop("no volcano here"))
    expect_error(suppressMessages(factory(results, "RES_vs_WT", "LCMS_Pos", 0.05)), "no volcano here")
})

test_that("(W2, D3) the in-process page factory gives trx's page for the contrast's entry (MCI-030's project)", {
    testthat::skip_if_not_installed("Glimma")
    testthat::skip_if_not_installed("edgeR")
    fixture <- module_ci_trx_project(module_ci_trx_corpus())
    owner <- fixture$owner
    suppressWarnings(suppressMessages(runTrxReferenceRuvs(owner, trxRuvGridParameters(), list())))
    results <- suppressWarnings(suppressMessages(runTrxReferenceDaAnalysis(owner, getWorkflowStateReference(owner),
        list(da_q_val_thresh = 0.05, treat_lfc_cutoff = 0), trxDaContrastEntries(fixture$workflow$contrasts_tbl))$results))
    seen <- omicsDaCoreRecordPages()
    factory <- workflowOmicsInProcessVolcanoPageFactory(transcriptomicsDeclarationBinding, owner, fixture$paths)
    # The tab maps the selected name to its entry before the factory (D3).
    entry <- trxDaContrastEntry(fixture$workflow, "Sera_vs_RPMI")
    page <- suppressWarnings(suppressMessages(factory(results, entry, "counts", 0.05)))
    expect_identical(seen$displays[[1L]]$operation, "trx_volcano_page")
    expect_identical(seen$displays[[1L]]$parameters, list(contrast = entry, assay = "counts", qvalue = 0.05))
    expect_identical(page, omicsDaCoreIframe(owner, seen$displays[[1L]]$display$page, entry))
})

# ------------------------------------------------------------------------------------------------
# The in-process export (D4).
# ------------------------------------------------------------------------------------------------

#' Pass-through recorders on the export's volcano and heatmap kernels, in W1 (b)'s text (`omicDaW1CallText`): the
#' kernel's name, its scalar arguments and a digest of what it returned. Masked (`guiEntryMask`), the digests read
#' `<HEX>`, as W1's frozen values do.
omicsDaCoreRecordKernels <- function(omic, env = parent.frame()) {
    acc <- new.env(parent = emptyenv())
    acc$calls <- character()
    short <- function(x) substr(guiEntryDigest(x), 1L, 16L)
    digest_of <- function(value) {
        if (inherits(value, "ggplot")) return(paste("ggplot", short(value$data)))
        if (is.list(value) && methods::is(value$plot, "Heatmap")) {
            return(sprintf("list(plot = Heatmap %s, row_clusters = %s)", short(value$plot@matrix), short(value$row_clusters)))
        }
        paste0("<", class(value)[[1L]], ">")
    }
    kernels <- if (identical(omic, "lipid")) c("generateLipidDAVolcanoStatic", "generateLipidDAHeatmap") else
        c("generateMetabDAVolcanoStatic", "generateMetabDAHeatmap")
    recorder <- function(name) {
        force(name)
        original <- get(name, envir = asNamespace("MultiScholaR"))
        function(...) {
            value <- original(...)
            args <- list(...)
            keep <- vapply(args, function(arg) is.atomic(arg) && !is.null(arg) && length(arg) <= 3L && is.null(dim(arg)), logical(1))
            named <- names(args) %||% rep("", length(args))
            shown <- vapply(which(keep), function(i) paste0(if (nzchar(named[[i]])) paste0(named[[i]], "=") else "",
                paste(format(args[[i]]), collapse = "|")), "")
            acc$calls <- c(acc$calls, sprintf("%s(%s) -> %s", name, paste(shown, collapse = ", "), digest_of(value)))
            value
        }
    }
    do.call(testthat::local_mocked_bindings, c(lapply(stats::setNames(kernels, kernels), recorder),
        list(.package = "MultiScholaR", .env = env)))
    acc
}

#' W1 (b)'s export calls of the run step (lipid_main, metab_main): per assay and contrast, the volcano and the heatmap.
omicsDaCoreExportCalls <- function(omic) {
    kernel <- if (identical(omic, "lipid")) "Lipid" else "Metab"
    labels <- if (identical(omic, "lipid")) "show_lipid_names" else "show_metabolite_names"
    calls <- character()
    for (assay in c("LCMS_Pos", "GCMS")) for (contrast in c("groupKO-groupWT", "groupRES-groupWT")) {
        calls <- c(calls,
            sprintf(paste0("generate%sDAVolcanoStatic(selected_contrast=%s, selected_assay=%s, da_q_val_thresh=0.05, ",
                "lfc_threshold=0, show_labels=TRUE, n_labels=15) -> ggplot <HEX>"), kernel, contrast, assay),
            sprintf(paste0("generate%sDAHeatmap(selected_contrast=%s, selected_assay=%s, top_n=50, clustering_method=ward.D2, ",
                "distance_method=euclidean, cluster_rows=TRUE, cluster_cols=TRUE, scale_data=row, color_scheme=RdBu, ",
                "%s=FALSE, da_q_val_thresh=0.05) -> list(plot = Heatmap <HEX>, row_clusters = <HEX>)"), kernel, contrast,
                assay, labels))
    }
    calls
}

#' W1 (b)'s files of the run step (lipid_main, metab_main) under `da_output/` and `publication_graphs/`: the tables
#' by digest, every other file present.
omicsDaCoreExportFiles <- function(omic) {
    noun <- if (identical(omic, "lipid")) "lipids" else "metabolites"
    tables <- if (identical(omic, "lipid")) {
        c(gcms_KO = "074f019210bb05e4db534de339e88fba7e321a2f2e9b0d50bb503d7d5bb8eb98",
            gcms_RES = "8597e926b23cddc4514ddec594433c871289662093c37c5f845074adafcdc8e3",
            posmode_KO = "3ef61cc44f484ab6f9f7c6a1fef40960ffe082360f5bed58ddbff740cdbcdd23",
            posmode_RES = "c10aa8e923a83d86c58b8e06db8c79d9917eeb185799fdb802622f5286c9e7e8")
    } else {
        c(gcms_KO = "8182dd547c4ecf1bfcd64eede3d2ca7fed28796be30db143270b4e1938b81f4a",
            gcms_RES = "751ac64ceb5f2d1bb8602c40f16806d39341ff1e4c25cdf5b9e31fb66604cfce",
            posmode_KO = "cc3b4d6c4323e14b01e7caa7f564239f0dccf2816bee6c3c4d191ad2d90b273b",
            posmode_RES = "b557e8c011c4f6809ca4d940a91686b44fdc3ed0b2af331174826e78ff3dcb7b")
    }
    files <- character()
    for (assay in c("gcms", "posmode")) for (contrast in c("KO", "RES")) {
        stem <- sprintf("de_%s_%s_group%s-groupWT_long_annot", assay, noun, contrast)
        files[[sprintf("da_output/%s.tsv", stem)]] <- tables[[paste0(assay, "_", contrast)]]
        files[[sprintf("da_output/%s.xlsx", stem)]] <- "<present>"
    }
    files[["publication_graphs/Heatmaps/all_heatmaps_combined.pdf"]] <- "<present>"
    for (assay in c("gcms", "posmode")) for (contrast in c("KO", "RES")) for (ext in c("pdf", "png")) {
        files[[sprintf("publication_graphs/Heatmaps/%s_group%s-groupWT_heatmap.%s", assay, contrast, ext)]] <- "<present>"
    }
    files[[sprintf("publication_graphs/NumSigDeMolecules/%s_num_sig_barplot.png", noun)]] <- "<present>"
    files[[sprintf("publication_graphs/NumSigDeMolecules/%s_num_sig_de_molecules.tab", noun)]] <-
        "a36fe768533ac33524d9bd3b2697cf5a138ab567180f4533e9cfd91d13a15713"
    files[[sprintf("publication_graphs/NumSigDeMolecules/%s_num_sig_de_molecules.xlsx", noun)]] <- "<present>"
    files[["publication_graphs/Volcano_Plots/all_volcano_plots_combined.pdf"]] <- "<present>"
    for (assay in c("gcms", "posmode")) for (contrast in c("KO", "RES")) for (ext in c("pdf", "png")) {
        files[[sprintf("publication_graphs/Volcano_Plots/%s_group%s-groupWT.%s", assay, contrast, ext)]] <- "<present>"
    }
    files[order(names(files))]
}

#' The export's files under `base` (W1 (b)'s `guiEntryArtifacts`), without the store and the session file.
omicsDaCoreFiles <- function(base, roots) {
    files <- unlist(guiEntryArtifacts(base, roots))
    files <- files[grepl("^(da_output|publication_graphs)/", names(files))]
    stats::setNames(as.character(files), names(files))
}

#' Fresh export directories under `root/<name>`.
omicsDaCoreExportPaths <- function(root, name) {
    paths <- list(da_output_dir = file.path(root, name, "da_output"),
        publication_graphs_dir = file.path(root, name, "publication_graphs"))
    for (path in paths) dir.create(path, recursive = TRUE)
    paths
}

#' The notice the tab shows for an export's display (`reportWorkflowOmicsDaExport`), as W1 (b) records it.
omicsDaCoreNotice <- function(display) {
    notes <- character()
    notify <- function(ui, duration = 5, type = "default", ...) notes <<- c(notes, sprintf("[%s %s -] %s", type, duration, ui))
    reportWorkflowOmicsDaExport(display, notify, logInfo = function(...) NULL, logWarn = function(...) NULL)
    notes
}

test_that("(W2, D4) the in-process DA publication makes W1 (b)'s export calls, files and notice for lipid and metab", {
    for (omic in c("lipid", "metab")) {
        fixture <- omicsDaCoreFeature(omic)
        omicsDaCoreRunFeature(fixture)
        acc <- omicsDaCoreRecordKernels(omic)
        paths <- omicsDaCoreExportPaths(fixture$root, "core")
        display <- suppressWarnings(suppressMessages(runWorkflowOmicsPublicationInProcess(fixture$binding,
            fixture$bus$state_manager, paste0(fixture$binding$prefix, "_da_publication"), list(), paths)))
        expect_identical(display, list(success = TRUE), info = omic)
        expect_identical(vapply(acc$calls, guiEntryMask, "", USE.NAMES = FALSE), omicsDaCoreExportCalls(omic), info = omic)
        expect_identical(omicsDaCoreFiles(file.path(fixture$root, "core"), fixture$root), omicsDaCoreExportFiles(omic),
            info = omic)
        expect_identical(omicsDaCoreNotice(display), "[message 5 -] DA results saved to disk (tables, volcano plots, heatmaps)",
            info = omic)
    }
})

#' The D4 comparison: the start commit's in-process write (`write(results, fixture, paths, notify)`) and the generic
#' runner on the same results make the same calls with the same arguments and data digests (unmasked), write the same
#' files and table bytes, and show the same notice.
omicsDaCoreCompareExport <- function(omic, write) {
    fixture <- omicsDaCoreFeature(omic)
    results <- omicsDaCoreRunFeature(fixture)
    acc <- omicsDaCoreRecordKernels(omic)
    before <- omicsDaCoreExportPaths(fixture$root, "before")
    notes <- character()
    notify <- function(ui, duration = 5, type = "default", ...) notes <<- c(notes, sprintf("[%s %s -] %s", type, duration, ui))
    suppressWarnings(suppressMessages(write(results, fixture, before, notify)))
    calls_before <- acc$calls
    acc$calls <- character()
    after <- omicsDaCoreExportPaths(fixture$root, "after")
    display <- suppressWarnings(suppressMessages(runWorkflowOmicsPublicationInProcess(fixture$binding,
        fixture$bus$state_manager, paste0(fixture$binding$prefix, "_da_publication"), list(), after)))
    expect_length(calls_before, 8L)
    expect_identical(acc$calls, calls_before)
    expect_identical(omicsDaCoreFiles(file.path(fixture$root, "after"), fixture$root),
        omicsDaCoreFiles(file.path(fixture$root, "before"), fixture$root))
    expect_identical(omicsDaCoreNotice(display), notes)
}

# Comparison side (lipid): quarantined when W3 deletes mod_lipid_da_session_helpers.R (`lipidDaResultsWithSource`).
test_that("(W2, D4; comparison side) lipid: the in-process publication equals the start commit's write", {
    testthat::skip("quarantined: Q-130 (KIT-106 W3) it compared the in-process DA publication with the start commit's lipid write, `outputLipidDaResultsAllContrasts` on `lipidDaResultsWithSource` (mod_lipid_da_session_helpers.R), which W3 deleted; the tab's export is `runWorkflowOmicsPublicationInProcess(..., \"lipid_da_publication\")` in the shared frame (mod_omic_da.R), whose core side stays live in the block above and in test-omic-da-module.R (b) lipid_main's run step")
    omicsDaCoreCompareExport("lipid", function(results, fixture, paths, notify) {
        # mod_lipid_da_analysis_helpers.R `finalizeLipidDaRunAnalysisSuccess`'s in-process write and its notice.
        written <- outputLipidDaResultsAllContrasts(da_results_list = lipidDaResultsWithSource(results, fixture$tab, fixture$bus),
            da_output_dir = paths$da_output_dir, publication_graphs_dir = paths$publication_graphs_dir,
            da_q_val_thresh = 0.05, lfc_threshold = 0, heatmap_top_n = 50, heatmap_clustering = "both",
            heatmap_color_scheme = "RdBu")
        if (isTRUE(written)) notify("DA results saved to disk (tables, volcano plots, heatmaps)", type = "message", duration = 5)
    })
})

# Comparison side (metab): quarantined when W4 deletes mod_metab_da_observer_helpers.R (`writeMetabDaResultsArtifacts`,
# `metabDaResultsWithSource`).
test_that("(W2, D4; comparison side) metab: the in-process publication equals the start commit's write", {
    omicsDaCoreCompareExport("metab", function(results, fixture, paths, notify) {
        writeMetabDaResultsArtifacts(results = metabDaResultsWithSource(results, fixture$tab, fixture$bus),
            experimentPaths = paths, daQValThresh = 0.05, treatLfcCutoff = 0, showNotification = notify)
    })
})
