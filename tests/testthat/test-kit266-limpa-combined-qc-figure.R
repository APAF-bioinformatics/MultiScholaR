# KIT-266 (VAL-308, VAL-316; AUDIT-KIT-266-plan-20261126): the DIA-with-limpa route builds the combined quality-control
# figure (Figure 1), so the limpa report has the Figure 1 its Results text describes.
#
# Rows E1 to E21 fail on the code before this change, at the values stated beside them, by assertion; rows G1 to G6
# pass on both sides. The new functions are looked up by name (`k266Fn`), so a row on the old code fails by its
# assertion and does not error. The fixtures and the drive are copied file-locally from test-kit259b-second-correction.R,
# test-prot-norm-worker-figures.R and test-kit227-report-ruv-sentences.R, and are never sourced from them.

k266Quiet <- function(code) suppressWarnings(suppressMessages(code))

#' A function of the package by name, or NULL when it does not exist (the old code).
k266Fn <- function(name) get0(name, envir = asNamespace("MultiScholaR"), inherits = FALSE)

k266Needs <- function() {
    skip_if_not_installed("limpa")
    skip_if_not_installed("mixOmics")
}

# The image titles of the ticket's section 5, as drafted (the owner sets the final wording, VAL-41).
k266Left <- "DPC-Quant estimates\n(log2)"
k266CentreNone <- "After the peptide count filter\n(no between-sample\nnormalisation applied)"
k266CentreMethod <- function(name) paste0("After ", name, " normalisation\nand the peptide count filter")
k266RightCovariates <- "RUV-III-adjusted\n(for quality control only;\nnot the data analysed)"
k266RightCorrection <- "RUV-III-corrected\n(the data analysed)"

k266Kinds <- c("pca", "density", "rle", "correlation")
k266TwoColumnFiles <- c("pre_norm_pca.png", "post_norm_pca.png", "pre_norm_density.png", "post_norm_density.png",
    "pre_norm_rle.png", "post_norm_rle.png", "pre_norm_correlation.png", "post_norm_correlation.png")
k266ThreeColumnFiles <- c("pre_norm_pca.png", "post_norm_pca.png", "ruv_corrected_pca.png",
    "pre_norm_density.png", "post_norm_density.png", "ruv_corrected_density.png",
    "pre_norm_rle.png", "post_norm_rle.png", "ruv_corrected_rle.png",
    "pre_norm_correlation.png", "post_norm_correlation.png", "ruv_corrected_correlation.png")
k266List <- "proteomics_composite_QC_figure.placeholders"

k266Composite <- function(fixture) file.path(fixture$paths$protein_qc_dir,
    paste0("proteomics_composite_QC_figure.", c("png", "pdf", "rds")))

#' The width and height in pixels of a PNG, from its header (NA when there is no file).
k266PngSize <- function(path) {
    if (!file.exists(path)) return(c(NA_integer_, NA_integer_))
    readBin(readBin(path, "raw", n = 24L)[17:24], "integer", n = 2L, size = 4L, endian = "big")
}

k266Dpc <- function(fixture) withDIAStateScope(fixture$manager, function()
    getDIAStateReference(fixture$manager, "protein_dpc_quant")$entity_id)

#' The persisted source of the pre-normalisation figures (its entity id), or NA.
k266PreSource <- function(fixture) withDIAStateScope(fixture$manager, function()
    fixture$manager$.context$snapshot$qc_sources$pre_normalization$entity_id %||% NA_character_)

#' The rendered assets the head roots.
k266Assets <- function(fixture) withDIAStateScope(fixture$manager, function() sum(vapply(
    unname(fixture$manager$.context$snapshot$products), function(reference) identical(reference$kind, "rendered_asset"), NA)))

#' A step of the drive from a clicker's step: the tab's session and its selected-tab value.
k266Step <- function(step) function(session, selected) step(session)

#' The tab driven through `steps` (named `function(session, selected)`), on `fixture`, in process or managed (each action
#' through the worker's entry point in this process, as protNormK122Drive does; the chained runs are given time to
#' finish, S14). Recorded: the composite builder's calls by process (`MULTISCHOLAR_WORKFLOW_WORKER`), the sources of
#' the QC plot products, and per step the composite's files, size and bytes, the pre-normalisation files, the list, the
#' managed actions and the rendered assets the head roots.
k266Drive <- function(fixture, steps, managed = FALSE, mocks = list()) {
    probe <- new.env()
    log <- new.env(parent = emptyenv())
    protNormW3Notices(log, fixture$root)
    if (managed) protNormW3Managed(fixture, probe)
    out <- new.env(parent = emptyenv())
    out$builds <- list()
    out$sources <- character()
    out$notes <- list()
    original <- generateProtNormCompositeFromFiles
    plot_of <- resolveDIAReferenceQcPlot
    testthat::local_mocked_bindings(
        generateProtNormCompositeFromFiles = function(plotFiles, ncol = 3, rowLabels = NULL, columnLabels = NULL, ...) {
            process <- if (identical(Sys.getenv("MULTISCHOLAR_WORKFLOW_WORKER"), "1")) "worker" else "application"
            out$builds[[length(out$builds) + 1L]] <- list(process = process, ncol = ncol, columns = columnLabels,
                letters = unlist(rowLabels, use.names = FALSE),
                files = unname(vapply(plotFiles, function(f) if (is.na(f)) NA_character_ else basename(f), "")))
            original(plotFiles = plotFiles, ncol = ncol, rowLabels = rowLabels, columnLabels = columnLabels, ...)
        },
        resolveDIAReferenceQcPlot = function(source, ...) {
            out$sources <- c(out$sources, source$entity_id)
            plot_of(source, ...)
        }, .package = "MultiScholaR")
    if (length(mocks)) do.call(testthat::local_mocked_bindings, c(mocks, list(.package = "MultiScholaR")))
    composite <- k266Composite(fixture)
    qc <- fixture$paths$protein_qc_dir
    selected <- shiny::reactiveVal("import")
    k266Quiet(utils::capture.output(shiny::testServer(mod_prot_norm_server, args = list(
        workflow_data = fixture$workflow, experiment_paths = fixture$paths, omic_type = "proteomics",
        experiment_label = "core", selected_tab = selected), {
        session$userData$dia_capacity_manager <- fixture$manager
        norm <- session$getReturned()
        for (case in names(steps)) {
            seen <- length(probe$actions)
            built <- length(out$builds)
            steps[[case]](session, selected)
            for (i in 1:3) {
                session$elapse(100)
                session$flushReact()
            }
            list_path <- file.path(qc, k266List)
            out$step[[case]] <- list(
                exists = file.exists(composite), md5 = unname(tools::md5sum(composite[[1L]])),
                bytes = file.size(composite), size = k266PngSize(composite[[1L]]),
                builds = out$builds[seq_along(out$builds) > built],
                actions = vapply(probe$actions[seq_along(probe$actions) > seen], `[[`, "", "operation"),
                list = if (file.exists(list_path)) readLines(list_path) else NULL,
                pre = file.exists(file.path(qc, paste0("pre_norm_", k266Kinds, ".png"))),
                pre_md5 = unname(tools::md5sum(file.path(qc, "pre_norm_pca.png"))),
                markers = file.exists(file.path(qc, paste0("pre_norm_", k266Kinds, ".png.placeholder"))),
                ruv_files = file.exists(file.path(qc, paste0("ruv_corrected_", k266Kinds, ".png"))),
                cancor = file.exists(file.path(qc, "ruv_covariates_cancor.png")),
                source = k266PreSource(fixture), generated = shiny::isolate(norm$pre_norm_qc_generated),
                assets = k266Assets(fixture))
        }
    })))
    list(steps = out$step, builds = out$builds, sources = out$sources, notes = log$notes)
}

#' The limpa lane in one RUV mode: the pass-through, the correlation filter, the RUV step (and `extra` steps after it).
k266LaneSteps <- function(mode, extra = list()) {
    click <- protNormW3Clicker()
    c(list(pass_through = k266Step(do.call(click, c(list("run_normalization"), protNormW3Inputs(mode)))),
        correlation = k266Step(click("apply_correlation_filter")), ruv = k266Step(click("apply_ruv_covariates"))),
        lapply(extra, function(make) make(click)))
}

#' The cached lane: a fresh owner per key, driven once per file (the rows compare with it).
k266Lane <- local({
    cache <- list()
    function(mode, managed = FALSE, again = NULL) {
        key <- paste(mode, managed, again)
        if (is.null(cache[[key]])) {
            drive <- function() {
                fixture <- protNormW1LimpaOwner(environment())
                extra <- list(ruv_repeat = function(click) k266Step(click("apply_ruv_covariates")))
                if (identical(again, "apply")) extra <- list(correlation_again = function(click)
                    k266Step(click("apply_correlation_filter")))
                if (identical(again, "skip")) extra <- list(correlation_skip = function(click)
                    k266Step(click("skip_correlation_filter")))
                k266Drive(fixture, k266LaneSteps(mode, extra), managed)
            }
            cache[[key]] <<- drive()
        }
        cache[[key]]
    }
})

k266First <- function(builds) if (length(builds)) builds[[1L]] else list()
k266Last <- function(builds) if (length(builds)) builds[[length(builds)]] else list()

# ------------------------------------------------------------------------------------------------
# E1 to E4: the build after the pass-through and after the RUV step, and its layout.
# ------------------------------------------------------------------------------------------------

for (managed in c(FALSE, TRUE)) {
    test_that(sprintf("KIT-266 E1 (%s): the pass-through press builds the combined figure, once", if (managed) "managed" else "in process"), {
        k266Needs()
        if (managed) skip_if_not(workflowOwnerCapable())
        step <- k266Lane("covariates", managed)$steps$pass_through
        expect_identical(step$exists, c(TRUE, TRUE, TRUE))
        expect_length(step$builds, 1L)
        expect_identical(k266First(step$builds)$process, if (managed) "worker" else "application")
    })
}

test_that("KIT-266 E2: the pass-through build is the two-column layout with the limpa image titles and letters a to h", {
    k266Needs()
    build <- k266First(k266Lane("covariates")$steps$pass_through$builds)
    expect_identical(build$ncol, 2)
    expect_identical(build$columns, c(k266Left, k266CentreNone))
    expect_identical(build$letters, paste0(letters[1:8], ")"))
    expect_identical(build$files, k266TwoColumnFiles)
    resolver <- k266Fn("resolveDiaLimpaCompositeFigureInputs")
    expect_true(is.function(resolver))
    skip_if(!is.function(resolver))
    qc <- withr::local_tempdir()
    for (method in c("none", "cyclicloess", "quantile", "scale")) {
        got <- resolver("skip", qc, method)
        centre <- if (identical(method, "none")) k266CentreNone else
            k266CentreMethod(c(cyclicloess = "cyclic loess", quantile = "quantile", scale = "scale")[[method]])
        expect_identical(got$columnLabels, c(k266Left, centre), info = method)
        expect_identical(got$ncol, 2, info = method)
    }
})

test_that("KIT-266 E3: the covariates step rebuilds the figure with three columns; the cancor panel m) only when its image exists", {
    k266Needs()
    lane <- k266Lane("covariates")
    step <- lane$steps$ruv
    build <- k266Last(step$builds)
    expect_identical(build$ncol, 3)
    expect_identical(build$columns, c(k266Left, k266CentreNone, k266RightCovariates))
    expect_identical(build$files[1:12], k266ThreeColumnFiles)
    expect_identical(any(build$files == "ruv_covariates_cancor.png", na.rm = TRUE), step$cancor)
    expect_identical(build$letters[1:12], paste0(letters[1:12], ")"))
    if (isTRUE(step$cancor)) expect_identical(build$letters[13:15], c("", "m)", ""))
    resolver <- k266Fn("resolveDiaLimpaCompositeFigureInputs")
    expect_true(is.function(resolver))
    skip_if(!is.function(resolver))
    qc <- withr::local_tempdir()
    without <- resolver("covariates", qc, "none")
    expect_length(without$plotFiles, 12L)
    expect_false("cancor" %in% names(without$rowLabels))
    file.create(file.path(qc, "ruv_covariates_cancor.png"))
    with <- resolver("covariates", qc, "none")
    expect_length(with$plotFiles, 15L)
    expect_identical(basename(with$plotFiles[[14L]]), "ruv_covariates_cancor.png")
    expect_identical(with$rowLabels$cancor, c("", "m)", ""))
    expect_identical(resolver("covariates", qc, "quantile")$columnLabels[[2L]], k266CentreMethod("quantile"))
})

test_that("KIT-266 E4: the correction step's third image title is 'RUV-III-corrected (the data analysed)'", {
    k266Needs()
    build <- k266Last(k266Lane("correction")$steps$ruv$builds)
    expect_identical(build$ncol, 3)
    expect_identical(build$columns, c(k266Left, k266CentreNone, k266RightCorrection))
})

# ------------------------------------------------------------------------------------------------
# E5 to E7: RUV skipped, no factor, Reset.
# ------------------------------------------------------------------------------------------------

test_that("KIT-266 E5: RUV skipped after the pass-through leaves the two-column figure, bytes unchanged", {
    k266Needs()
    fixture <- protNormW1LimpaOwner()
    click <- protNormW3Clicker()
    got <- k266Drive(fixture, list(pass_through = k266Step(do.call(click, c(list("run_normalization"), protNormW3Inputs("skip")))),
        correlation = k266Step(click("apply_correlation_filter")),
        ruv_skip = k266Step(click("apply_ruv_covariates", ruv_mode = "skip"))))
    expect_identical(got$steps$pass_through$exists, c(TRUE, TRUE, TRUE))
    expect_identical(got$steps$ruv_skip$md5, got$steps$pass_through$md5)
    expect_identical(got$steps$ruv_skip$ruv_files, rep(FALSE, 4L))
    expect_true(all(vapply(got$builds, function(build) identical(build$ncol, 2), NA)) && length(got$builds) > 0L)
})

test_that("KIT-266 E6: a covariates step that found no factor (k = 0) leaves the two-column figure", {
    k266Needs()
    fixture <- protNormW1LimpaOwner()
    click <- protNormW3Clicker()
    none <- function(settings, workflow, norm_data, ...) {
        norm_data$best_k <- 0L
        norm_data$control_genes_index <- NA
        norm_data$ruv_optimization_result <- list(best_percentage = NA, replicate_structure = "group")
        TRUE
    }
    got <- k266Drive(fixture, list(pass_through = k266Step(do.call(click, c(list("run_normalization"), protNormW3Inputs("covariates")))),
        correlation = k266Step(click("apply_correlation_filter")),
        ruv_none = k266Step(click("apply_ruv_covariates"))), mocks = list(runDIAReferenceRuvCovariates = none))
    expect_identical(got$steps$ruv_none$exists, c(TRUE, TRUE, TRUE))
    expect_identical(got$steps$ruv_none$md5, got$steps$pass_through$md5)
    expect_true(all(vapply(got$builds, function(build) identical(build$ncol, 2), NA)) && length(got$builds) > 0L)
})

test_that("KIT-266 E7: Reset, then a head press after a covariates run: the RUV figures go (KIT-257) and the figure is two-column again", {
    k266Needs()
    fixture <- protNormW1LimpaOwner()
    click <- protNormW3Clicker()
    got <- k266Drive(fixture, list(pass_through = k266Step(do.call(click, c(list("run_normalization"), protNormW3Inputs("covariates")))),
        correlation = k266Step(click("apply_correlation_filter")), ruv = k266Step(click("apply_ruv_covariates")),
        reset = k266Step(click("reset_normalization")), press = k266Step(click("run_normalization"))))
    expect_identical(got$steps$press$exists, c(TRUE, TRUE, TRUE))
    expect_identical(got$steps$press$ruv_files, rep(FALSE, 4L))
    expect_identical(k266Last(got$steps$press$builds)$ncol, 2)
    expect_false(identical(got$steps$press$md5, got$steps$ruv$md5))
})

# ------------------------------------------------------------------------------------------------
# E8, E8b: the left column.
# ------------------------------------------------------------------------------------------------

for (managed in c(FALSE, TRUE)) {
    test_that(sprintf("KIT-266 E8 (%s): the head press draws the pre-normalisation figures from the DPC-Quant state", if (managed) "managed" else "in process"), {
        k266Needs()
        if (managed) skip_if_not(workflowOwnerCapable())
        fixture <- protNormW1LimpaOwner()
        click <- protNormW3Clicker()
        got <- k266Drive(fixture, list(pass_through = k266Step(do.call(click, c(list("run_normalization"), protNormW3Inputs("covariates"))))), managed)
        step <- got$steps$pass_through
        expect_identical(step$pre, rep(TRUE, 4L))
        expect_identical(step$source, k266Dpc(fixture))
        expect_true(isTRUE(step$generated))
        expect_identical(got$sources[[1L]] %||% NA_character_, k266Dpc(fixture))
    })
}

test_that("KIT-266 E8b: after a head press a tab re-entry at correlation_filtered redraws nothing; the source stays the DPC-Quant state", {
    k266Needs()
    fixture <- protNormW1LimpaOwner()
    click <- protNormW3Clicker()
    got <- k266Drive(fixture, list(pass_through = k266Step(do.call(click, c(list("run_normalization"), protNormW3Inputs("covariates")))),
        correlation = k266Step(click("apply_correlation_filter")),
        tab_entry = function(session, selected) selected("normalization")))
    expect_identical(got$steps$tab_entry$source, k266Dpc(fixture))
    expect_identical(got$steps$tab_entry$pre_md5, got$steps$correlation$pre_md5)
    expect_false(is.na(got$steps$tab_entry$pre_md5))
})

# ------------------------------------------------------------------------------------------------
# E9, E9b: left cells from another table, or no DPC-Quant state: labelled notes.
# ------------------------------------------------------------------------------------------------

#' The lane to the covariates step on a fresh owner, then `tamper(fixture)` in the owner's scope, then the wrapper
#' called directly for the three-column layout. Returns the warnings the call raised, its error and the project's files.
k266Tampered <- function(tamper) {
    fixture <- protNormW1LimpaOwner(parent.frame())
    k266Drive(fixture, k266LaneSteps("covariates"))
    withDIAStateScope(fixture$manager, function() {
        snapshot <- fixture$manager$.context$snapshot
        snapshot <- tamper(fixture$manager, snapshot)
        publishDIAStateHead(fixture$manager$.context, snapshot)
    })
    build <- k266Fn("buildDiaLimpaCombinedQcFigure")
    warnings <- list()
    error <- if (is.function(build)) tryCatch({
        withCallingHandlers(k266Quiet(build(fixture$paths, "proteomics", "covariates", fixture$manager)),
            warning = function(w) {
                warnings[[length(warnings) + 1L]] <<- w
                invokeRestart("muffleWarning")
            })
        NULL
    }, error = function(e) e) else simpleError("buildDiaLimpaCombinedQcFigure does not exist")
    qc <- fixture$paths$protein_qc_dir
    list(error = error, warnings = warnings, fixture = fixture,
        markers = file.exists(file.path(qc, paste0("pre_norm_", k266Kinds, ".png.placeholder"))),
        list = if (file.exists(file.path(qc, k266List))) readLines(file.path(qc, k266List)) else NULL,
        composite = file.exists(k266Composite(fixture)))
}

k266LeftLines <- c("a\tpre_norm_pca", "d\tpre_norm_density", "g\tpre_norm_rle", "j\tpre_norm_correlation")

test_that("KIT-266 E9: with no DPC-Quant state in the history the four left cells are labelled notes and the list names a, d, g, j", {
    k266Needs()
    got <- k266Tampered(function(manager, snapshot) {
        snapshot$state_history <- setdiff(snapshot$state_history, "protein_dpc_quant")
        snapshot$qc_sources$pre_normalization <- NULL
        snapshot
    })
    expect_null(got$error, label = conditionMessage(got$error %||% simpleError("")))
    expect_identical(got$markers, rep(TRUE, 4L))
    expect_true(all(got$composite))
    expect_identical(got$list, k266LeftLines)
    expect_length(Filter(function(w) inherits(w, "prot_norm_figure_unavailable"), got$warnings), 0L)
})

test_that("KIT-266 E9b: a session whose pre-normalisation source is the correlation-filtered table gets notes, not the stale images", {
    k266Needs()
    got <- k266Tampered(function(manager, snapshot) {
        snapshot$qc_sources$pre_normalization <- getDIAStateReference(manager, "correlation_filtered")
        snapshot
    })
    expect_null(got$error, label = conditionMessage(got$error %||% simpleError("")))
    expect_identical(got$markers, rep(TRUE, 4L))
    expect_true(all(got$composite))
    expect_identical(got$list, k266LeftLines)
    expect_length(Filter(function(w) inherits(w, "prot_norm_figure_unavailable"), got$warnings), 0L)
})

# ------------------------------------------------------------------------------------------------
# E10, E11: Copy to Results Summary and the session archive carry the figure.
# ------------------------------------------------------------------------------------------------

test_that("KIT-266 E10: Copy to Results Summary after the covariates step delivers the combined figure", {
    k266Needs()
    fixture <- protNormW1LimpaOwner()
    k266Drive(fixture, k266LaneSteps("covariates"))
    root <- withr::local_tempdir()
    paths <- list(results_dir = file.path(root, "results"),
        results_summary_dir = file.path(root, "results_summary", "proteomics_kit266"),
        publication_graphs_dir = file.path(root, "publication_graphs"), time_dir = file.path(root, "time"),
        qc_dir = file.path(root, "qc"), da_output_dir = file.path(root, "da"), pathway_dir = file.path(root, "pathway"),
        source_dir = file.path(root, "source"), feature_qc_dir = fixture$paths$protein_qc_dir,
        subfeature_qc_dir = file.path(root, "peptide_qc"))
    for (path in setdiff(unlist(paths), fixture$paths$protein_qc_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    testthat::local_mocked_bindings(getProjectPaths = function(...) paths, .package = "MultiScholaR")
    printed <- tryCatch(k266Quiet(utils::capture.output(copyToResultsSummary(omic_type = "proteomics", experiment_label = "kit266",
        contrasts_tbl = data.frame(contrasts = "groupB-groupA"), design_matrix = data.frame(Run = "S1", group = "A"),
        force = TRUE))), error = function(e) conditionMessage(e))
    figures <- file.path(paths$results_summary_dir, "QC_figures")
    expect_true(file.exists(file.path(figures, "proteomics_composite_QC_figure.png")))
    expect_true(file.exists(file.path(figures, "proteomics_composite_QC_figure.pdf")))
    expect_false(any(grepl("Source file not found.*[Cc]omposite|[Cc]omposite.*Source file not found", printed)))
})

test_that("KIT-266 E11: a session exported after the covariates step and imported again carries the combined figure, bytes equal", {
    k266Needs()
    fixture <- protNormW1LimpaOwner()
    k266Drive(fixture, k266LaneSteps("covariates"))
    manager <- fixture$manager
    manager$.context$session_format <- "native_v1"
    restored <- tryCatch({
        # The native export as KIT-271's test runs it (`collectDIANativeSessionMetadata`, `collectDIANativeGuiFiles`,
        # `runDIANativeSessionExport`), then the import and the restore of the GUI files.
        norm <- new.env(parent = emptyenv())
        archive <- NULL
        k266Quiet(withDIAStateScope(manager, function() {
            norm$correlation_filtered_obj <- getDIAStateReference(manager)
            collected <- collectDIANativeSessionMetadata(fixture$workflow, norm, list(norm_method = "none", ruv_mode = "covariates"))
            files <- collectDIANativeGuiFiles(manager, collected$metadata)
            source_dir <- fixture$paths$source_dir
            result <- runDIANativeSessionExport(manager, collected$metadata, source_dir, extra_references = collected$references,
                dependencies = files, prepare_outputs = function(metadata, request) {
                    saveProtNormExportMetadataFiles(metadata, source_dir, messageFn = function(...) NULL, strict = TRUE)
                    writeLines("kit266", file.path(source_dir, "filtered_session_summary.txt"))
                })
            archive <<- result$publication$archive_path
        }))
        fresh <- file.path(withr::local_tempdir(), "fresh")
        paths <- list(source_dir = file.path(fresh, "scripts", "proteomics"),
            protein_qc_dir = file.path(fresh, "results", "proteomics", "protein_qc"),
            peptide_qc_dir = file.path(fresh, "results", "proteomics", "peptide_qc"))
        for (path in c(unlist(paths), file.path(dirname(paths$protein_qc_dir), "publication_graphs")))
            dir.create(path, recursive = TRUE, showWarnings = FALSE)
        imported <- k266Quiet(openDIANativeSession(archive, paths$source_dir))
        withr::defer(closeDIAStateAdapter(imported$manager))
        k266Quiet(restoreDIANativeGuiFiles(imported, paths))
        file.path(paths$protein_qc_dir, paste0("proteomics_composite_QC_figure.", c("png", "pdf", "rds")))
    }, error = function(e) e)
    expect_false(inherits(restored, "error"), label = if (inherits(restored, "error")) conditionMessage(restored) else "")
    skip_if(inherits(restored, "error"))
    expect_identical(file.exists(restored), c(TRUE, TRUE, TRUE))
    expect_identical(unname(tools::md5sum(restored)), unname(tools::md5sum(k266Composite(fixture))))
})

# ------------------------------------------------------------------------------------------------
# E12: KIT-259's list of replaced panels on limpa.
# ------------------------------------------------------------------------------------------------

k266FailCorrelation <- function(state, env = parent.frame()) {
    original <- buildProtNormCorrelationPlot
    testthat::local_mocked_bindings(buildProtNormCorrelationPlot = function(...) {
        if (isTRUE(state$fail)) stop("injected drawing failure")
        original(...)
    }, .package = "MultiScholaR", .env = env)
}

test_that("KIT-266 E12 (a), (c): a replaced RUV panel writes the list; a second run that draws it leaves no stale marker or list", {
    k266Needs()
    state <- new.env()
    k266FailCorrelation(state)
    fixture <- protNormW1LimpaOwner()
    click <- protNormW3Clicker()
    got <- k266Drive(fixture, list(pass_through = k266Step(do.call(click, c(list("run_normalization"), protNormW3Inputs("covariates")))),
        correlation = k266Step(click("apply_correlation_filter")),
        ruv_failing = function(session, selected) {
            state$fail <- TRUE
            click("apply_ruv_covariates")(session)
        },
        ruv_drawn = function(session, selected) {
            state$fail <- FALSE
            click("apply_ruv_covariates")(session)
        }))
    expect_identical(got$steps$ruv_failing$list, "l\truv_corrected_correlation")
    expect_null(got$steps$ruv_drawn$list)
    expect_false(file.exists(file.path(fixture$paths$protein_qc_dir, "ruv_corrected_correlation.png.placeholder")))
})

test_that("KIT-266 E12 (b): a replaced canonical-correlation panel is named through the limpa resolver", {
    k266Needs()
    resolver <- k266Fn("resolveDiaLimpaCompositeFigureInputs")
    expect_true(is.function(resolver))
    skip_if(!is.function(resolver))
    qc <- withr::local_tempdir()
    file.create(file.path(qc, c(k266ThreeColumnFiles, "ruv_covariates_cancor.png", "proteomics_composite_QC_figure.png")))
    file.create(file.path(qc, "ruv_covariates_cancor.png.placeholder"))
    written <- tryCatch(writeProtNormCompositePlaceholders("covariates", qc, "proteomics",
        resolveInputsFn = function(ruvMode, qcDir) resolver(ruvMode, qcDir, "none")), error = function(e) e)
    expect_false(inherits(written, "error"), label = if (inherits(written, "error")) conditionMessage(written) else "")
    path <- file.path(qc, k266List)
    expect_true(file.exists(path))
    if (file.exists(path)) expect_identical(readLines(path), "m\truv_covariates_cancor")
})

# ------------------------------------------------------------------------------------------------
# E13, E14, G3: the limpa template's Results text, and the covariates summary text.
# ------------------------------------------------------------------------------------------------

#' A chunk of a template, extracted with knitr::purl and parsed.
k266Exprs <- function(path, label) {
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    start <- grep(paste0("^```\\{r ", label, "[,}]"), lines)
    stopifnot(length(start) == 1L)
    end <- start + which(lines[-seq_len(start)] == "```")[[1L]]
    chunk <- tempfile(fileext = ".Rmd")
    writeLines(lines[start:end], chunk, useBytes = TRUE)
    code <- tempfile(fileext = ".R")
    suppressMessages(knitr::purl(chunk, output = code, quiet = TRUE, documentation = 0))
    parse(code, keep.source = FALSE, encoding = "UTF-8")
}

#' The functions the limpa template defines in its Result Section and Main chunks, in one environment.
k266LimpaEnv <- function() {
    path <- system.file("reports", "proteomics", "DIANN_limpa_report.rmd", package = "MultiScholaR")
    env <- new.env(parent = globalenv())
    env$params <- list(workflow_name = "maxquant", omic_type = "proteomics", experiment_label = "kit266")
    env$d66_log <- function(...) invisible(NULL)
    env$tibble <- tibble::tibble
    for (label in c("Result Section", "Main Execution Function")) {
        for (e in k266Exprs(path, label)) {
            if (is.call(e) && identical(e[[1L]], as.name("<-")) && is.call(e[[3L]]) &&
                identical(e[[3L]][[1L]], as.name("function"))) eval(e, env)
        }
    }
    env
}

#' The `results` element of Main's report content (the Results call as the template makes it).
k266ResultsCall <- function(env) {
    find <- function(expr) {
        if (!is.call(expr)) return(NULL)
        if (!is.null(names(expr)) && all(c("results", "method_details") %in% names(expr))) return(expr[["results"]])
        for (i in seq_along(expr)) if (is.call(expr[[i]])) {
            hit <- find(expr[[i]])
            if (!is.null(hit)) return(hit)
        }
        NULL
    }
    call <- find(body(env$Main))
    stopifnot(!is.null(call))
    call
}

#' The Results text for `numbers`, as the template's Main prints it.
k266Results <- function(env, numbers) {
    results_env <- new.env(parent = env)
    results_env$analysis_results <- list(numbers = numbers, table = data.frame())
    results_env$correlation_data <- list()
    paste(as.character(eval(k266ResultsCall(env), results_env)), collapse = "")
}

k266LeftWords <- "DPC-Quant's protein estimates (a complete log2 matrix; model-based where peptides were not detected)"

#' The panel letters the text names for each column (left, centre or right), in the order they are printed.
k266TextLetters <- function(text, column) {
    pattern <- paste0("\\(", column, " column: (panels )?([a-z, ]+)\\)")
    hit <- regmatches(text, regexec(pattern, text))[[1L]]
    if (length(hit) < 3L) return(NA_character_)
    strsplit(hit[[3L]], ", ", fixed = TRUE)[[1L]]
}

#' The letters of a layout's columns, from its resolver: left, centre (and right), without the cancor panel.
k266ResolverLetters <- function(resolver, mode) {
    qc <- withr::local_tempdir()
    got <- resolver(mode, qc, "none")
    rows <- got$rowLabels[setdiff(names(got$rowLabels), "cancor")]
    columns <- lapply(seq_len(got$ncol), function(j) sub("[)]$", "", vapply(rows, function(row) row[[j]], "")))
    stats::setNames(lapply(columns, unname), c("left", "centre", "right")[seq_len(got$ncol)])
}

k266Recorded <- list(dpc_engine = "fastlimpadpc", limpa_norm_method = "none")

test_that("KIT-266 E13: the limpa template's Results text names the stages, the method and the panel letters the figure has", {
    k266Needs()
    resolver <- k266Fn("resolveDiaLimpaCompositeFigureInputs")
    expect_true(is.function(resolver))
    skip_if(!is.function(resolver))
    env <- k266LimpaEnv()
    centre_none <- "after the peptide count filter, with no between-sample normalisation applied"
    covariates <- k266Results(env, c(k266Recorded, list(ruvcov_k = 2, ruvcov_mode = "covariates")))
    expect_true(grepl(paste0("Figure 1 shows four quality-control views at three stages of processing: ", k266LeftWords,
        " (left column: panels a, d, g, j); ", centre_none, ", before the sample-correlation filter (centre column: b, e, h, k); ",
        "and for the RUV-III-adjusted matrix after the sample-correlation filter (drawn for quality control only; the ",
        "differential abundance model carried the RUV-III factors as covariates) (right column: c, f, i, l)."),
        covariates, fixed = TRUE), label = covariates)
    letters3 <- k266ResolverLetters(resolver, "covariates")
    for (column in c("left", "centre", "right")) expect_identical(k266TextLetters(covariates, column), letters3[[column]], info = column)
    correction <- k266Results(env, c(k266Recorded, list(ruvcov_k = 2, ruvcov_mode = "correction")))
    expect_true(grepl(paste0("and after RUV-III correction (the corrected matrix entered the differential abundance model), ",
        "applied after the sample-correlation filter (right column: c, f, i, l)."), correction, fixed = TRUE), label = correction)
    expect_false(grepl("carried the RUV-III factors as covariates", correction, fixed = TRUE))
    two <- paste0("Figure 1 shows four quality-control views at two stages of processing: ", k266LeftWords,
        " (left column: panels a, c, e, g) and ", centre_none, ", before the sample-correlation filter (right column: b, d, f, h). ",
        "RUV-III correction was not applied to this dataset.")
    letters2 <- k266ResolverLetters(resolver, "skip")
    # Every two-column state: no covariates step, a record without the factor count, k = 0 (mode none), a skipped RUV.
    for (numbers in list(k266Recorded, c(k266Recorded, list(ruvcov_k = "not recorded")),
        c(k266Recorded, list(ruvcov_k = 0, ruvcov_mode = "none")))) {
        text <- k266Results(env, numbers)
        expect_true(grepl(two, text, fixed = TRUE), label = text)
        expect_identical(k266TextLetters(text, "left"), letters2$left)
        expect_identical(k266TextLetters(text, "right"), letters2$centre)
    }
    # The method branch: the centre phrase follows the recorded between-sample method.
    quantile <- k266Results(env, c(k266Recorded[1L], list(limpa_norm_method = "quantile", ruvcov_k = 2, ruvcov_mode = "covariates")))
    expect_true(grepl("after the quantile normalisation and the peptide count filter, before the sample-correlation filter (centre column: b, e, h, k)",
        quantile, fixed = TRUE), label = quantile)
    loess <- k266Results(env, c(k266Recorded[1L], list(limpa_norm_method = "cyclicloess")))
    expect_true(grepl("after the cyclic loess normalisation and the peptide count filter, before the sample-correlation filter (right column: b, d, f, h)",
        loess, fixed = TRUE), label = loess)
    # No DPC engine recorded: the text is as it was.
    expect_true(grepl("at each stage of processing drawn for this analysis", k266Results(env, list(ruv_applied = NA)), fixed = TRUE))
})

test_that("KIT-266 E13 (caption): the limpa template's Figure 1 caption names the sample-correlation filter; the other three keep theirs", {
    caption <- function(file) {
        lines <- readLines(system.file("reports", "proteomics", file, package = "MultiScholaR"), warn = FALSE, encoding = "UTF-8")
        grep("^```\\{r QC_metrics,", lines, value = TRUE)
    }
    expect_match(caption("DIANN_limpa_report.rmd"), paste0('fig.cap="Figure 1: QC metrics across normalisation stages. On DIA with limpa ',
        'the left and centre columns are before the sample-correlation filter and the right column is after it."'), fixed = TRUE)
    for (file in c("LFQ_report.rmd", "TMT_report.rmd", "DIANN_report.rmd")) {
        expect_match(caption(file), 'fig.cap="Figure 1: QC metrics across normalisation stages"', fixed = TRUE, info = file)
    }
})

test_that("KIT-266 E14: the covariates summary text names the RUV-III-adjusted QC figures", {
    texts <- list()
    normData <- list2env(list(best_k = 2L, control_genes_index = c(TRUE, FALSE, TRUE),
        ruv_optimization_result = list(best_percentage = 10, replicate_structure = "group")))
    output <- list2env(list())
    k266Quiet(finishProtNormRuvCovariates("accepted", "covariates", output, NULL, normData,
        renderTextFn = function(text) {
            texts$summary <<- text
            text
        }, showNotificationFn = function(...) NULL, messageFn = function(...) NULL))
    expect_match(texts$summary, "RUV-III-adjusted QC figures", fixed = TRUE)
    expect_false(grepl("RUV-corrected QC plots", texts$summary, fixed = TRUE))
})

# Declared edit after KIT-272 (VAL-332; AUDIT-KIT-272-impl I07): KIT-272 replaced the PCA bullet and the density bullet in
# every state of the three classic templates, so the base copies (tests/testdata/kit236/base) are compared after a base-side
# substitution of those two bullets by the head's, each with an exactly-one-match assertion. Everything else must be equal.
k266Bullets <- function(text, last) {
    first <- gregexpr("• Principal component analysis", text, fixed = TRUE)[[1L]]
    end <- gregexpr(paste0(last, "[^\n]*"), text)[[1L]]
    list(first = first, end = end, count = c(sum(first > 0L), sum(end > 0L)),
        span = if (sum(first > 0L) == 1L && sum(end > 0L) == 1L)
            substr(text, first[[1L]], end[[1L]] + attr(end, "match.length")[[1L]] - 1L) else NA_character_)
}

test_that("KIT-266 G3: the three classic templates' Results text is the base's for the three RUV states, but for KIT-272's two bullets", {
    states <- list(TRUE, FALSE, NA)
    for (file in c("LFQ_report.rmd", "TMT_report.rmd", "DIANN_report.rmd")) {
        texts <- lapply(c(head = system.file("reports", "proteomics", file, package = "MultiScholaR"),
            base = testthat::test_path("..", "testdata", "kit236", "base", file)), function(path) {
            env <- new.env(parent = globalenv())
            for (e in k266Exprs(path, "Result Section")) eval(e, env)
            lapply(states, function(state) paste(as.character(env$FormatResults(list(), ruv_applied = state)), collapse = ""))
        })
        for (i in seq_along(states)) {
            info <- paste(file, "state", i)
            new <- k266Bullets(texts$head[[i]], "• Principal component scores by group")
            old <- k266Bullets(texts$base[[i]], "• Density plots")
            expect_identical(new$count, c(1L, 1L), info = info)
            expect_identical(old$count, c(1L, 1L), info = info)
            substituted <- if (is.na(new$span) || is.na(old$span)) NA_character_ else
                paste0(substr(texts$base[[i]], 1L, old$first[[1L]] - 1L), new$span,
                    substring(texts$base[[i]], old$end[[1L]] + attr(old$end, "match.length")[[1L]]))
            expect_identical(substituted, texts$head[[i]], info = info)
        }
        expect_true(all(nzchar(unlist(texts$head))), info = file)
    }
})

# ------------------------------------------------------------------------------------------------
# E15: the correlation filter re-applied after a covariates run (the figure only; the record is F8, outside KIT-266).
# ------------------------------------------------------------------------------------------------

for (case in list(list(again = "apply", managed = FALSE), list(again = "apply", managed = TRUE),
    list(again = "skip", managed = FALSE), list(again = "skip", managed = TRUE))) {
    test_that(sprintf("KIT-266 E15 (%s, %s): after a covariates run the re-applied correlation filter leaves the two-column figure",
        case$again, if (case$managed) "managed" else "in process"), {
        k266Needs()
        if (case$managed) skip_if_not(workflowOwnerCapable())
        lane <- k266Lane("covariates", case$managed, case$again)
        step <- lane$steps[[if (identical(case$again, "apply")) "correlation_again" else "correlation_skip"]]
        expect_identical(step$exists, c(TRUE, TRUE, TRUE))
        expect_length(step$builds, 1L)
        expect_identical(k266First(step$builds)$ncol, 2)
        expect_identical(k266First(step$builds)$process, if (case$managed) "worker" else "application")
        expect_false(identical(step$md5, lane$steps$ruv$md5))
        expect_identical(step$size[[1L]], 1500L)
    })
}

# ------------------------------------------------------------------------------------------------
# E16 to E21: the managed route builds the figure on the worker.
# ------------------------------------------------------------------------------------------------

test_that("KIT-266 E16: managed pass-through: the pre figures, the source's figures, then one limpa_composite action built by the worker", {
    k266Needs()
    skip_if_not(workflowOwnerCapable())
    managed <- k266Lane("covariates", managed = TRUE)
    unmanaged <- k266Lane("covariates", managed = FALSE)
    step <- managed$steps$pass_through
    expect_identical(step$actions, rep("normalization-qc-figures", 3L))
    workers <- Filter(function(build) identical(build$process, "worker"), step$builds)
    expect_length(workers, 1L)
    expect_length(Filter(function(build) identical(build$process, "application"), step$builds), 0L)
    expect_identical(step$exists, c(TRUE, TRUE, TRUE))
    expect_identical(step$md5, unmanaged$steps$pass_through$md5)
    expect_identical(step$assets - unmanaged$steps$pass_through$assets, 3L)
})

test_that("KIT-266 E17: managed covariates step: a limpa_composite action after the stage action, built by the worker; a repeat finds it saved", {
    k266Needs()
    skip_if_not(workflowOwnerCapable())
    managed <- k266Lane("covariates", managed = TRUE)
    unmanaged <- k266Lane("covariates", managed = FALSE)
    step <- managed$steps$ruv
    expect_identical(step$actions, c("protein_ruv_covariates", "normalization-qc-figures", "normalization-qc-figures"))
    expect_identical(vapply(step$builds, `[[`, "", "process"), "worker")
    expect_identical(k266Last(step$builds)$ncol, 3)
    expect_identical(step$md5, unmanaged$steps$ruv$md5)
    repeated <- managed$steps$ruv_repeat
    expect_identical(repeated$actions, c("protein_ruv_covariates", "normalization-qc-figures", "normalization-qc-figures"))
    expect_length(repeated$builds, 0L)
    expect_identical(repeated$assets, step$assets)
    expect_identical(step$assets - unmanaged$steps$ruv$assets, 6L)
})

test_that("KIT-266 E18: managed correlation re-apply after a covariates run dispatches one limpa_composite (two columns); a first correlation none", {
    k266Needs()
    skip_if_not(workflowOwnerCapable())
    lane <- k266Lane("covariates", managed = TRUE, again = "apply")
    figures <- "normalization-qc-figures"
    expect_identical(lane$steps$correlation$actions, c("pearson-pairs", "pearson-selection", figures))
    expect_length(lane$steps$correlation$builds, 0L)
    again <- lane$steps$correlation_again
    expect_identical(again$actions, c("pearson-pairs", "pearson-selection", figures, figures))
    expect_identical(vapply(again$builds, `[[`, "", "process"), "worker")
    expect_identical(k266First(again$builds)$ncol, 2)
    skipped <- k266Lane("covariates", managed = TRUE, again = "skip")$steps$correlation_skip
    expect_identical(skipped$actions, figures)
    expect_identical(k266First(skipped$builds)$ncol, 2)
})

test_that("KIT-266 E19: when the worker saved nothing the delivery builds in the application and logs why", {
    k266Needs()
    skip_if_not(workflowOwnerCapable())
    saving <- k266Fn("saveWorkflowDiaLimpaCompositeFigure")
    expect_true(is.function(saving))
    skip_if(!is.function(saving))
    logged <- character()
    old <- logger::log_appender(function(lines) logged <<- c(logged, lines))
    withr::defer(logger::log_appender(old))
    fixture <- protNormW1LimpaOwner()
    click <- protNormW3Clicker()
    got <- k266Drive(fixture, list(pass_through = k266Step(do.call(click, c(list("run_normalization"), protNormW3Inputs("covariates"))))),
        managed = TRUE, mocks = list(saveWorkflowDiaLimpaCompositeFigure = function(manager, request) invisible(NULL)))
    step <- got$steps$pass_through
    expect_identical(vapply(step$builds, `[[`, "", "process"), "application")
    expect_identical(step$exists, c(TRUE, TRUE, TRUE))
    expect_true(any(grepl("not found saved", logged, fixed = TRUE)), label = paste(logged, collapse = " | "))
})

test_that("KIT-266 E20: GUI audit A111-01 on limpa: the final figure has the layout of the final head", {
    k266Needs()
    skip_if_not(workflowOwnerCapable())
    fixture <- protNormW1LimpaOwner()
    probe <- new.env()
    log <- new.env(parent = emptyenv())
    protNormW3Notices(log, fixture$root)
    protNormW3Managed(fixture, probe)
    final <- new.env(parent = emptyenv())
    k266Quiet(utils::capture.output(shiny::testServer(mod_prot_norm_server, args = list(
        workflow_data = fixture$workflow, experiment_paths = fixture$paths, omic_type = "proteomics",
        experiment_label = "core"), {
        session$userData$dia_capacity_manager <- fixture$manager
        flush <- function() for (i in 1:3) {
            session$elapse(100)
            session$flushReact()
        }
        do.call(session$setInputs, protNormW3Inputs("covariates"))
        session$setInputs(run_normalization = 1)
        flush()
        session$setInputs(apply_correlation_filter = 1)
        flush()
        session$setInputs(apply_ruv_covariates = 1)
        flush()
        # The correlation re-apply is in flight when the covariates button is pressed.
        probe$hold <- TRUE
        session$setInputs(apply_correlation_filter = 2)
        session$setInputs(apply_ruv_covariates = 2)
        probe$hold <- FALSE
        flush()
        final$size <- k266PngSize(k266Composite(fixture)[[1L]])
        final$operations <- vapply(probe$actions, `[[`, "", "operation")
    })))
    expect_identical(final$size[[1L]], 1950L)
    expect_length(final$operations, 16L)
})

test_that("KIT-266 E21: the labelled notes are in the project before the managed composite action starts; the figure carries them", {
    k266Needs()
    skip_if_not(workflowOwnerCapable())
    fixture <- protNormW1LimpaOwner()
    click <- protNormW3Clicker()
    got <- k266Drive(fixture, list(pass_through = k266Step(do.call(click, c(list("run_normalization"), protNormW3Inputs("covariates")))),
        correlation = k266Step(click("apply_correlation_filter")),
        tampered = function(session, selected) withDIAStateScope(fixture$manager, function() {
            snapshot <- fixture$manager$.context$snapshot
            snapshot$qc_sources$pre_normalization <- getDIAStateReference(fixture$manager, "correlation_filtered")
            publishDIAStateHead(fixture$manager$.context, snapshot)
        }),
        ruv = k266Step(click("apply_ruv_covariates"))), managed = TRUE)
    step <- got$steps$ruv
    expect_identical(step$markers, rep(TRUE, 4L))
    expect_identical(step$list, k266LeftLines)
    expect_false(identical(step$md5, got$steps$correlation$md5))
    expect_identical(vapply(step$builds, `[[`, "", "process"), "worker")
})

# ------------------------------------------------------------------------------------------------
# Guards G1 to G6 (they pass on the old code too).
# ------------------------------------------------------------------------------------------------

test_that("KIT-266 G1: the classic DIA route calls the builder with the classic layout, captions and letters", {
    fixture <- protNormCoreOwner()
    click <- protNormW3Clicker()
    got <- k266Drive(fixture, list(skip = k266Step(do.call(click, c(list("run_normalization"), protNormW3Inputs("skip")))),
        manual = k266Step(click("run_normalization", ruv_mode = "manual"))))
    skipped <- k266First(got$steps$skip$builds)
    expect_length(got$steps$skip$builds, 1L)
    expect_identical(skipped$columns, c("Pre-Normalisation", "Post-Normalisation"))
    expect_identical(skipped$files, k266TwoColumnFiles)
    expect_identical(skipped$letters, paste0(letters[1:8], ")"))
    manual <- k266First(got$steps$manual$builds)
    expect_length(got$steps$manual$builds, 1L)
    expect_identical(manual$columns, c("Pre-Normalisation", "Post-Normalisation", "RUV-Corrected"))
    expect_identical(manual$files, c(k266ThreeColumnFiles, NA, "ruv_corrected_cancor.png", NA))
    expect_identical(manual$letters, c(paste0(letters[1:12], ")"), "", "m)", ""))
})

test_that("KIT-266 G2: Reset after a limpa head press returns to the table the press designated", {
    k266Needs()
    fixture <- protNormW1LimpaOwner()
    click <- protNormW3Clicker()
    before <- fixture$manager$current_state
    k266Drive(fixture, list(pass_through = k266Step(do.call(click, c(list("run_normalization"), protNormW3Inputs("covariates")))),
        reset = k266Step(click("reset_normalization"))))
    expect_identical(fixture$manager$current_state, before)
})

test_that("KIT-266 G2b: a restored session with no recorded input and a persisted source at correlation_filtered resets to correlation_filtered", {
    k266Needs()
    fixture <- protNormW1LimpaOwner()
    click <- protNormW3Clicker()
    got <- k266Drive(fixture, list(pass_through = k266Step(do.call(click, c(list("run_normalization"), protNormW3Inputs("covariates")))),
        correlation = k266Step(click("apply_correlation_filter")),
        restored = function(session, selected) {
            norm <- session$getReturned()
            norm$normalization_input_reference <- NULL
            withDIAStateScope(fixture$manager, function() {
                snapshot <- fixture$manager$.context$snapshot
                snapshot$qc_sources$pre_normalization <- getDIAStateReference(fixture$manager, "correlation_filtered")
                publishDIAStateHead(fixture$manager$.context, snapshot)
            })
        },
        reset = k266Step(click("reset_normalization"))))
    expect_identical(fixture$manager$current_state, "correlation_filtered")
})

test_that("KIT-266 G4: a press on the tab's own RUV state rebuilds nothing; the figure's bytes are unchanged", {
    k266Needs()
    fixture <- protNormW1LimpaOwner()
    got <- k266Drive(fixture, k266LaneSteps("covariates", list(press = function(click) k266Step(click("run_normalization")))))
    expect_identical(got$steps$press$md5, got$steps$ruv$md5)
    expect_length(got$steps$press$builds, 0L)
})

test_that("KIT-266 G5: the covariates observer with no project paths raises no error notice and builds nothing", {
    captured <- new.env(parent = emptyenv())
    fake_runner <- function(settings, workflowData, normData, experimentPaths, omicType, experimentLabel, ruv_qc, grouping) {
        normData$best_k <- 1L
        normData$control_genes_index <- c(TRUE, FALSE, TRUE)
        normData$ruv_optimization_result <- list(best_percentage = 5, replicate_structure = "group")
        "accepted"
    }
    testthat::local_mocked_bindings(runDIAReferenceRuvCovariates = fake_runner, updateProtNormRuvAuditTrail = function(...) invisible(NULL))
    workflow <- list2env(list(state_manager = list(workflow_type = "DIA_LIMPA")))
    for (mode in c("correction", "covariates")) {
        norm_data <- list2env(list(correlation_filtering_complete = TRUE))
        notes <- character()
        runProtNormRuvCovariatesObserver(list(ruv_mode = mode, color_variable = NULL, shape_variable = NULL), list2env(list()),
            workflow, norm_data, list(), "proteomics", "x", function(reference) NULL, function() "group",
            renderTextFn = function(text) {
                captured$text <- text
                text
            }, showNotificationFn = function(message, ...) notes <<- c(notes, message), messageFn = function(...) NULL)
        expect_false(any(grepl("^Error", notes)), info = mode)
        expect_match(captured$text, "State saved as:", fixed = TRUE)
    }
})

test_that("KIT-266 G6: a first correlation filter after the pass-through, and again after a Reset, rebuilds nothing", {
    k266Needs()
    fixture <- protNormW1LimpaOwner()
    click <- protNormW3Clicker()
    got <- k266Drive(fixture, list(pass_through = k266Step(do.call(click, c(list("run_normalization"), protNormW3Inputs("covariates")))),
        correlation = k266Step(click("apply_correlation_filter")),
        reset = k266Step(click("reset_normalization")),
        press = k266Step(click("run_normalization")),
        correlation_after_reset = k266Step(click("apply_correlation_filter"))))
    expect_length(got$steps$correlation$builds, 0L)
    expect_identical(got$steps$correlation$md5, got$steps$pass_through$md5)
    expect_length(got$steps$correlation_after_reset$builds, 0L)
    expect_identical(got$steps$correlation_after_reset$md5, got$steps$press$md5)
})
