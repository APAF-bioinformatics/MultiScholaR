# KIT-259, second correction (VAL-287; AUDIT-KIT-259-impl-20261126 I02, I04, I05, I11; the owner's item (ii)).
#
# I02: the lipid/metab histogram's x breaks for a lowest correlation in [-0.3, -0.2) hold 5.55e-17 for 0; the labels
#      are expected to turn scientific (by reading). A3-E5 establishes it on the old code by the printed label.
# I04: a placeholder written without a notice: the DIA-with-limpa RUV covariates step, a managed pre-figure redraw.
# I05: the figure catch also took the head-ownership conditions and R's allocation error.
# I11: the partial PNG of a failed ggsave is replaced by the placeholder's bytes (a guard of the merged fix).
# (ii): the proteomics report names the Figure 1 panels replaced by a placeholder; with none its text is unchanged.
#
# E rows fail on the code before this correction; G rows pass on both.

k259bQuiet <- function(code) suppressWarnings(suppressMessages(code))

k259bRun <- function(code) {
    warnings <- list()
    value <- tryCatch(withCallingHandlers(code, warning = function(w) {
        warnings[[length(warnings) + 1L]] <<- w
        invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage")), error = function(e) e)
    list(value = value, warnings = warnings)
}

# ------------------------------------------------------------------------------------------------
# I02: A3-E5.
# ------------------------------------------------------------------------------------------------

k259bOmicsPlot <- function(omic) {
    # Group A: correlation -0.2085 (lowest; the axis's lower end is -0.3); group B: 0.997.
    assay <- data.frame(feature_id = paste0("F", 1:5), S1 = c(1, 2, 3, 4, 5), S2 = c(3, 5, 2.5, 1, 4),
        S3 = c(1, 2, 3, 4, 5), S4 = c(1.1, 2.2, 2.9, 4.2, 5.1), stringsAsFactors = FALSE)
    design <- data.frame(Run = c("S1", "S2", "S3", "S4"), group = c("A", "A", "B", "B"),
        tech_rep = c("t1", "t2", "t3", "t4"), stringsAsFactors = FALSE)
    object <- if (identical(omic, "lipid")) {
        createLipidomicsAssayData(lipid_data = list(Assay1 = assay), design_matrix = design, lipid_id_column = "feature_id",
            sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep")
    } else {
        createMetaboliteAssayData(metabolite_data = list(Assay1 = assay), design_matrix = design,
            metabolite_id_column = "feature_id", sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep")
    }
    k259bQuiet(plotPearson(object, correlation_group = "group"))$Assay1
}

for (omic in c("lipid", "metab")) {
    test_that(sprintf("KIT-259 A3-E5 (%s): a lowest correlation in [-0.3, -0.2) labels the x axis in plain numbers", omic), {
        plot <- k259bOmicsPlot(omic)
        expect_s3_class(plot, "ggplot")
        expect_equal(plot$scales$get_scales("x")$limits[[1L]], -0.3)
        labels <- as.character(ggplot2::get_guide_data(plot, "x")$.label)
        expect_false(any(grepl("e[-+]", labels)), label = paste(labels, collapse = " "))
        expect_true(any(suppressWarnings(as.numeric(labels)) == 0, na.rm = TRUE), label = paste(labels, collapse = " "))
    })
}

# ------------------------------------------------------------------------------------------------
# I05: the catch re-raises the head-ownership conditions and R's allocation error.
# ------------------------------------------------------------------------------------------------

test_that("KIT-259 I05-E1: workflow_execution_stale, workflow_execution_busy and an allocation error are not replaced", {
    for (case in list(
        list(class = "workflow_execution_stale", fn = function() rlang::abort("workflow head changed before mutation",
            class = "workflow_execution_stale")),
        list(class = "workflow_execution_busy", fn = function() rlang::abort("another controller holds the project",
            class = "workflow_execution_busy")),
        list(class = "error", fn = function() stop("cannot allocate vector of size 8.0 Gb")))) {
        dir <- withr::local_tempdir()
        got <- k259bRun(drawProtNormQcFigure(dir, "post_norm_pca.png", 8, 6, case$fn))
        expect_true(inherits(got$value, case$class) && inherits(got$value, "error"), info = case$class)
        expect_false(file.exists(file.path(dir, "post_norm_pca.png")), info = case$class)
    }
})

test_that("KIT-259 I05-G1: a drawing error is still replaced; with no folder it is logged only, nothing written", {
    dir <- withr::local_tempdir()
    got <- k259bRun(drawProtNormQcFigure(dir, "post_norm_pca.png", 8, 6, function() stop("drawing failed")))
    expect_identical(got$value, file.path(dir, "post_norm_pca.png"))
    expect_true(file.exists(file.path(dir, "post_norm_pca.png")))
    none <- k259bRun(drawProtNormQcFigure(NULL, "post_norm_pca.png", 8, 6, function() stop("drawing failed")))
    expect_null(none$value)
    expect_length(Filter(function(w) inherits(w, "prot_norm_figure_unavailable"), none$warnings), 0L)
})

# ------------------------------------------------------------------------------------------------
# I11: the partial file of a failed ggsave is replaced by the placeholder's bytes.
# ------------------------------------------------------------------------------------------------

test_that("KIT-259 I11-G1: 16 groups with the default 'Shape by': the PCA file is the placeholder's bytes", {
    recipe <- readWorkflowCheckpointRecipe("small")
    recipe$n_samples <- 32L
    object <- makeWorkflowCheckpointProtein(recipe)
    groups <- sprintf("g%02d", rep(seq_len(16L), each = 2L))
    design <- object@design_matrix
    design$group <- groups
    design$tech_rep_group <- paste(groups, seq_along(groups), sep = "_")
    object@design_matrix <- design
    dir <- withr::local_tempdir()
    got <- k259bRun(generateProtNormPostNormalizationQcArtifacts(object, dir, list(color_var = "group", shape_var = "group"),
        messageFn = function(...) invisible(NULL), gcFn = function() invisible(NULL)))
    expect_false(inherits(got$value, "error"))
    unavailable <- Filter(function(w) inherits(w, "prot_norm_figure_unavailable"), got$warnings)
    expect_length(unavailable, 1L)
    if (length(unavailable) == 1L) {
        expect_identical(unavailable[[1L]]$figure, "post_norm_pca.png")
        expected <- file.path(withr::local_tempdir(), "post_norm_pca.png")
        k259bQuiet(writeProtNormFigurePlaceholder(expected, "post_norm_pca.png", unavailable[[1L]]$reason, 8, 6))
        expect_identical(unname(tools::md5sum(file.path(dir, "post_norm_pca.png"))), unname(tools::md5sum(expected)))
    }
})

# ------------------------------------------------------------------------------------------------
# I04: a notice where a placeholder was written without one.
# ------------------------------------------------------------------------------------------------

k259bFailCorrelation <- function(state, env = parent.frame()) {
    original <- buildProtNormCorrelationPlot
    testthat::local_mocked_bindings(buildProtNormCorrelationPlot = function(...) {
        if (isTRUE(state$fail)) stop("injected drawing failure")
        original(...)
    }, .package = "MultiScholaR", .env = env)
}

#' The tab through `steps` (each a function(session)); the notices after each.
k259bDrive <- function(fixture, steps, managed = FALSE) {
    log <- new.env(parent = emptyenv())
    protNormW3Notices(log, fixture$root)
    probe <- new.env()
    if (managed) protNormW3Managed(fixture, probe)
    acc <- new.env(parent = emptyenv())
    acc$out <- list()
    protNormCoreQuiet(utils::capture.output(shiny::testServer(mod_prot_norm_server, args = list(
        workflow_data = fixture$workflow, experiment_paths = fixture$paths, omic_type = "proteomics",
        experiment_label = "core"), {
        for (case in names(steps)) {
            log$notes <- character()
            steps[[case]](session)
            session$flushReact()
            acc$out[[case]] <- log$notes
        }
    })))
    acc$out
}

k259bCovariatesSteps <- function(state) {
    click <- protNormW3Clicker()
    covariates <- click("apply_ruv_covariates")
    list(run = do.call(click, c(list("run_normalization"), protNormW3Inputs("covariates"))),
        correlation = click("apply_correlation_filter"),
        covariates = function(session) {
            state$fail <- TRUE
            covariates(session)
        })
}

# Third pass (VAL-307, W11, declared edit): at the covariates step the figures are named as the RUV-III-adjusted matrix
# drawn for quality control (they were "RUV-corrected").
k259bCovariatesWarning <- "^\\[warning [^]]+\\] .*'RUV-III-adjusted \\(for quality control\\) Pearson correlation' could not be drawn: injected drawing failure"

for (managed in c(FALSE, TRUE)) {
    test_that(sprintf("KIT-259 I04-E1 (%s): the DIA-with-limpa RUV covariates step names a figure it replaced",
        if (managed) "managed" else "in process"), {
        skip_if_not_installed("limpa")
        skip_if_not_installed("mixOmics")
        if (managed) skip_if_not(workflowOwnerCapable())
        state <- new.env()
        k259bFailCorrelation(state)
        notes <- k259bDrive(protNormW1LimpaOwner(), k259bCovariatesSteps(state), managed = managed)
        expect_false(any(grepl("Error", notes$covariates, fixed = TRUE)), label = paste(notes$covariates, collapse = " | "))
        expect_identical(sum(grepl(k259bCovariatesWarning, notes$covariates)), 1L, label = paste(notes$covariates, collapse = " | "))
    })
}

test_that("KIT-259 I04-E2: managed, a pre-normalisation figure redrawn because the input changed is named in a notice", {
    captured <- new.env()
    shown <- character()
    testthat::local_mocked_bindings(showNotification = function(ui, ..., type = "default") {
        shown <<- c(shown, paste0("[", type, "] ", paste(as.character(ui), collapse = " ")))
        invisible("id")
    }, .package = "shiny")
    testthat::local_mocked_bindings(
        withDIAStateScope = function(manager, fn) fn(),
        assertDIAExportReady = function(...) invisible(NULL),
        diaNormalizationRequest = function(...) list(source = "source", parameters = list(), redraw_pre = TRUE),
        protNormFiguresAction = function(source, aesthetics, figures = "stage", deliver = function(result) invisible(NULL),
            composite = NULL) list(source = source, deliver = deliver),
        .package = "MultiScholaR")
    # KIT-257: the runner may receive a function giving the actions when they start.
    runner <- list(run = function(actions, done = NULL, failed = NULL)
        captured$actions <- if (is.function(actions)) actions() else actions)
    preQc <- function() {
        rlang::warn("The QC figure 'pre_norm_pca.png' could not be drawn: injected", class = "prot_norm_figure_unavailable",
            figure = "pre_norm_pca.png", reason = "injected")
        invisible(NULL)
    }
    dispatchProtNormNormalization(runner, list(ruv_mode = "skip", color_variable = "group", shape_variable = "group"),
        list(state_manager = NULL), NULL, NULL, "proteomics", "core", function(...) NULL, function(...) NULL,
        function() "group", function(...) NULL, function(...) NULL, preQc = preQc)
    expect_identical(captured$actions[[1L]]$source, "source")
    got <- k259bRun(captured$actions[[1L]]$deliver(NULL))
    expect_length(Filter(function(w) inherits(w, "prot_norm_figure_unavailable"), got$warnings), 0L)
    expect_identical(sum(grepl("^\\[warning\\] The figure 'Pre-normalisation PCA' could not be drawn: injected", shown)), 1L)
})

# ------------------------------------------------------------------------------------------------
# (ii): the placeholder's marker and the report sentence.
# ------------------------------------------------------------------------------------------------

test_that("KIT-259 (ii)-E1: a placeholder leaves a marker beside it, which a later drawn figure removes", {
    dir <- withr::local_tempdir()
    marker <- file.path(dir, "post_norm_pca.png.placeholder")
    k259bRun(drawProtNormQcFigure(dir, "post_norm_pca.png", 8, 6, function() stop("drawing failed")))
    expect_true(file.exists(marker))
    k259bRun(drawProtNormQcFigure(dir, "post_norm_pca.png", 8, 6, function() {
        k259bQuiet(ggplot2::ggsave(file.path(dir, "post_norm_pca.png"), ggplot2::ggplot(), width = 8, height = 6, dpi = 150))
        file.path(dir, "post_norm_pca.png")
    }))
    expect_false(file.exists(marker))
})

test_that("KIT-259 (ii)-E2: a normalisation run removes the earlier run's post-normalisation and RUV markers, not the pre ones", {
    dir <- withr::local_tempdir()
    names <- c("post_norm_correlation.png.placeholder", "ruv_corrected_pca.png.placeholder", "pre_norm_pca.png.placeholder")
    file.create(file.path(dir, names))
    removeProtNormOwnFiles(list(protein_qc_dir = dir), "normalization")
    expect_identical(file.exists(file.path(dir, names)), c(FALSE, FALSE, TRUE))
})

# ------------------------------------------------------------------------------------------------
# Third pass (VAL-307; closure check C01): the replaced panels are recorded with Figure 1. The publication writes
# `proteomics_composite_QC_figure.placeholders` (one line per replaced panel: its letter and its figure) beside the
# combined figure when it builds it; KIT-257's start-of-run removal takes it; Copy and the session archive carry it; the
# four proteomics templates read it beside the combined figure they show.
# ------------------------------------------------------------------------------------------------

k259bList <- "proteomics_composite_QC_figure.placeholders"

#' The drive's capture of the list in the project's QC folder after each step (its lines, or NULL).
k259bListAfter <- function(fixture, steps, managed = FALSE) {
    path <- file.path(fixture$paths$protein_qc_dir, k259bList)
    wrapped <- lapply(steps, function(step) function(session) {
        step(session)
        session$flushReact()
        fixture$lists[[length(fixture$lists) + 1L]] <- if (file.exists(path)) readLines(path) else NULL
    })
    fixture$lists <- list()
    k259bDrive(fixture, wrapped, managed = managed)
    fixture$lists
}

k259bFailingRun <- function(state, click, fail = TRUE) {
    run <- do.call(click, c(list("run_normalization"), protNormW3Inputs("skip")))
    function(session) {
        state$fail <- fail
        run(session)
    }
}

for (managed in c(FALSE, TRUE)) {
    test_that(sprintf("KIT-259 C01-E1 (%s): the run writes the replaced panels beside the combined figure; a redraw leaves them",
        if (managed) "managed" else "in process"), {
        if (managed) skip_if_not(workflowOwnerCapable())
        state <- new.env()
        k259bFailCorrelation(state)
        click <- protNormW3Clicker()
        fixture <- protNormCoreOwner()
        fixture <- list2env(fixture)
        lists <- k259bListAfter(fixture, list(
            run = k259bFailingRun(state, click),
            redraw = function(session) {
                state$fail <- FALSE
                session$setInputs(color_variable = "batch")
            }), managed = managed)
        expect_identical(lists[[1L]], "h\tpost_norm_correlation")
        expect_identical(lists[[2L]], "h\tpost_norm_correlation")
        expect_true(file.exists(file.path(fixture$paths$protein_qc_dir, "proteomics_composite_QC_figure.png")))
    })
}

test_that("KIT-259 C01-E2: a redraw that fails after a run whose figures drew writes no list; pre markers enter the list", {
    state <- new.env()
    k259bFailCorrelation(state)
    click <- protNormW3Clicker()
    fixture <- list2env(protNormCoreOwner())
    lists <- k259bListAfter(fixture, list(
        run = k259bFailingRun(state, click, fail = FALSE),
        redraw = function(session) {
            state$fail <- TRUE
            session$setInputs(color_variable = "batch")
        }))
    expect_null(lists[[1L]])
    expect_null(lists[[2L]])
    expect_true(file.exists(file.path(fixture$paths$protein_qc_dir, "post_norm_correlation.png.placeholder")))
    # A pre-normalisation figure replaced earlier (its marker in the project folder) is one of Figure 1's panels.
    state$fail <- FALSE
    pre <- list2env(protNormCoreOwner())
    k259bQuiet(ggplot2::ggsave(file.path(pre$paths$protein_qc_dir, "pre_norm_pca.png"), ggplot2::ggplot(),
        width = 8, height = 6, dpi = 150))
    writeLines("injected", file.path(pre$paths$protein_qc_dir, "pre_norm_pca.png.placeholder"))
    lists <- k259bListAfter(pre, list(run = k259bFailingRun(new.env(), protNormW3Clicker(), fail = FALSE)))
    expect_identical(lists[[1L]], "a\tpre_norm_pca")
})

test_that("KIT-259 C01-E3: a normalisation run removes the earlier list; the pre markers stay", {
    dir <- withr::local_tempdir()
    file.create(file.path(dir, c(k259bList, "pre_norm_pca.png.placeholder")))
    removeProtNormOwnFiles(list(protein_qc_dir = dir), "normalization")
    expect_identical(file.exists(file.path(dir, c(k259bList, "pre_norm_pca.png.placeholder"))), c(FALSE, TRUE))
})

test_that("KIT-259 C01-E4: Copy carries the list into QC_figures beside the combined figure", {
    root <- withr::local_tempdir()
    paths <- list(results_dir = file.path(root, "results"),
        results_summary_dir = file.path(root, "results_summary", "proteomics_kit259"),
        publication_graphs_dir = file.path(root, "publication_graphs"), time_dir = file.path(root, "time"),
        qc_dir = file.path(root, "qc"), da_output_dir = file.path(root, "da"), pathway_dir = file.path(root, "pathway"),
        source_dir = file.path(root, "source"), feature_qc_dir = file.path(root, "feature_qc"),
        subfeature_qc_dir = file.path(root, "peptide_qc"))
    for (path in unlist(paths)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    for (ext in c("png", "pdf")) k259bQuiet(ggplot2::ggsave(file.path(paths$feature_qc_dir,
        paste0("proteomics_composite_QC_figure.", ext)), ggplot2::ggplot(), width = 4, height = 3, dpi = 50))
    writeLines("h\tpost_norm_correlation", file.path(paths$feature_qc_dir, k259bList))
    testthat::local_mocked_bindings(getProjectPaths = function(...) paths, .package = "MultiScholaR")
    k259bQuiet(utils::capture.output(copyToResultsSummary(omic_type = "proteomics", experiment_label = "kit259",
        contrasts_tbl = data.frame(contrasts = "groupB-groupA"), design_matrix = data.frame(Run = "S1", group = "A"),
        force = TRUE)))
    copied <- file.path(paths$results_summary_dir, "QC_figures", k259bList)
    expect_true(file.exists(copied))
    if (file.exists(copied)) expect_identical(readLines(copied), "h\tpost_norm_correlation")
})

test_that("KIT-259 C01-E5: a session exported after the run and imported again keeps the list", {
    state <- new.env()
    k259bFailCorrelation(state)
    click <- protNormW3Clicker()
    fixture <- list2env(protNormCoreOwner())
    k259bDrive(fixture, list(run = k259bFailingRun(state, click),
        correlation = function(session) {
            state$fail <- FALSE
            click("apply_correlation_filter")(session)
        },
        export = click("export_filtered_session")))
    archives <- list.files(fixture$root, pattern = "[.]msdia[.]tar$", recursive = TRUE, full.names = TRUE)
    expect_true(length(archives) >= 1L)
    skip_if(!length(archives))
    archive <- archives[order(file.mtime(archives), decreasing = TRUE)][[1L]]
    fresh <- file.path(withr::local_tempdir(), "fresh")
    paths <- list(source_dir = file.path(fresh, "scripts", "proteomics"),
        protein_qc_dir = file.path(fresh, "results", "proteomics", "protein_qc"),
        peptide_qc_dir = file.path(fresh, "results", "proteomics", "peptide_qc"))
    for (path in c(unlist(paths), file.path(dirname(paths$protein_qc_dir), "publication_graphs"))) {
        dir.create(path, recursive = TRUE, showWarnings = FALSE)
    }
    imported <- k259bQuiet(openDIANativeSession(archive, paths$source_dir))
    withr::defer(closeDIAStateAdapter(imported$manager))
    k259bQuiet(restoreDIANativeGuiFiles(imported, paths))
    restored <- file.path(paths$protein_qc_dir, k259bList)
    expect_true(file.exists(restored))
    if (file.exists(restored)) expect_identical(readLines(restored), "h\tpost_norm_correlation")
})

test_that("KIT-259 C01-G1: DIA with limpa builds no combined figure, so no list is written even when a figure is replaced", {
    skip_if_not_installed("limpa")
    skip_if_not_installed("mixOmics")
    state <- new.env()
    k259bFailCorrelation(state)
    fixture <- protNormW1LimpaOwner()
    k259bDrive(fixture, k259bCovariatesSteps(state))
    expect_false(file.exists(file.path(fixture$paths$protein_qc_dir, k259bList)))
    expect_false(file.exists(file.path(fixture$paths$protein_qc_dir, "proteomics_composite_QC_figure.png")))
})

k259bTemplates <- c(LFQ = "LFQ_report.rmd", TMT = "TMT_report.rmd", DIANN = "DIANN_report.rmd",
    DIANN_limpa = "DIANN_limpa_report.rmd")

#' The lines of one chunk of a proteomics template, or NULL when the template has none of that label.
k259bChunk <- function(file, label) {
    lines <- readLines(system.file("reports", "proteomics", file, package = "MultiScholaR"), warn = FALSE, encoding = "UTF-8")
    start <- grep(paste0("^```\\{r ", label, "[,}]"), lines)
    if (length(start) != 1L) return(NULL)
    end <- start + which(lines[-seq_len(start)] == "```")[[1L]]
    lines[start:end]
}

#' What the template's placeholder chunk prints: `list` written beside a combined figure (when `composite`) in the
#' report's QC_figures folder, and per-figure `markers` in the project's QC folder (the second correction's source).
k259bSentence <- function(file, list = NULL, composite = TRUE, markers = character()) {
    chunk <- k259bChunk(file, "QC_metrics_placeholders")
    if (is.null(chunk)) return(NULL)
    figures <- withr::local_tempdir()
    feature <- withr::local_tempdir()
    if (composite) file.create(file.path(figures, "proteomics_composite_QC_figure.png"))
    if (length(list)) writeLines(list, file.path(figures, k259bList))
    if (length(markers)) file.create(file.path(feature, paste0(markers, ".png.placeholder")))
    env <- new.env(parent = globalenv())
    env$params <- list(omic_type = "proteomics")
    env$report_qc_figures_dir <- figures
    env$current_project_paths <- list(feature_qc_dir = feature)
    env$results_data <- list(numbers = list(ruv_applied = FALSE))
    rmd <- tempfile(fileext = ".Rmd")
    writeLines(chunk, rmd, useBytes = TRUE)
    code <- tempfile(fileext = ".R")
    suppressMessages(knitr::purl(rmd, output = code, quiet = TRUE, documentation = 0))
    paste(utils::capture.output(for (e in parse(code, keep.source = FALSE, encoding = "UTF-8")) eval(e, env)),
        collapse = "\n")
}

for (name in names(k259bTemplates)) {
    test_that(sprintf("KIT-259 (ii)-E3 (%s): the report names Figure 1's replaced panels from the list beside it", name), {
        file <- k259bTemplates[[name]]
        expect_false(is.null(k259bChunk(file, "QC_metrics_placeholders")))
        skip_if(is.null(k259bChunk(file, "QC_metrics_placeholders")))
        expect_match(k259bSentence(file, c("a\tpre_norm_pca", "h\tpost_norm_correlation")), paste0("In Figure 1, panels a ",
            "\\(principal component analysis after log2 transformation\\) and h \\(Pearson correlation after between-sample ",
            "normalisation\\) could not be drawn by the application; each shows a note giving the reason instead of the plot."))
        expect_match(k259bSentence(file, c("a\tpre_norm_pca", "k\tpost_norm_correlation", "m\truv_corrected_cancor")),
            "panels a \\(.*\\), k \\(Pearson correlation after between-sample normalisation\\) and m \\(RUV-III canonical correlation\\)")
        expect_match(k259bSentence(file, "f\tpost_norm_rle"), paste0("In Figure 1, panel f \\(relative log expression ",
            "after between-sample normalisation\\) could not be drawn by the application; it shows a note"))
    })

    test_that(sprintf("KIT-259 C01-E6 (%s): the tab's current markers are not read; a list without a combined figure prints nothing",
        name), {
        file <- k259bTemplates[[name]]
        skip_if(is.null(k259bChunk(file, "QC_metrics_placeholders")))
        expect_identical(k259bSentence(file, markers = c("pre_norm_pca", "post_norm_correlation")), "")
        expect_identical(k259bSentence(file, "h\tpost_norm_correlation", composite = FALSE), "")
        expect_identical(k259bSentence(file, NULL), "")
    })

    test_that(sprintf("KIT-259 (ii)-G1 (%s): with no placeholder the Figure 1 section knits to the same text", name), {
        file <- k259bTemplates[[name]]
        lines <- readLines(system.file("reports", "proteomics", file, package = "MultiScholaR"), warn = FALSE,
            encoding = "UTF-8")
        start <- grep("^```\\{r QC_metrics,", lines)
        end <- grep("^# Differential Abundance Analysis", lines)
        expect_length(start, 1L)
        expect_true(length(end) >= 1L)
        region <- lines[start:(end[end > start][[1L]] - 1L)]
        chunk <- k259bChunk(file, "QC_metrics_placeholders")
        without <- if (is.null(chunk)) region else region[!seq_along(region) %in%
            (match(chunk[[1L]], region) + seq_along(chunk) - 1L)]
        # One folder for both knits: the chunk's message names it.
        qc_figures <- withr::local_tempdir()
        feature_qc <- withr::local_tempdir()
        knit <- function(text) {
            env <- new.env(parent = globalenv())
            env$params <- list(omic_type = "proteomics")
            env$report_qc_figures_dir <- qc_figures
            env$current_project_paths <- list(feature_qc_dir = feature_qc)
            env$results_data <- list(numbers = list(ruv_applied = TRUE))
            out <- suppressMessages(knitr::knit(text = paste(text, collapse = "\n"), envir = env, quiet = TRUE))
            gsub("\n{2,}", "\n\n", trimws(out))
        }
        expect_identical(knit(region), knit(without))
    })
}
