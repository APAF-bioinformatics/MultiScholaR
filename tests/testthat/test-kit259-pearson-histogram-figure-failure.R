# KIT-259: the normalisation QC figures that silently dropped data, and a proteomics normalisation figure that could not
# be drawn after the state was committed.
#
# A1: the proteomics within-group Pearson histogram fixed its count axis at c(0, 4), so a 0.001 bin holding more than 4
#     pairs was removed at drawing (rebuildWorkflowDiaQcPlot, the figure every stage draws; plotPearson, the producer).
# A2: no within-group pair with a finite correlation gave an error at seq(Inf, ...) in both.
# A3: the lipidomics and metabolomics histogram drew [0, 1] only (a negative correlation was not drawn) and dropped the
#     pairs with no correlation without a word.
# B2: a QC figure of the normalisation that cannot be drawn is replaced, under its own file name and size, by a
#     placeholder image naming it and the reason; the step completes with one warning notice.
#
# E rows fail on the code before the fix at the value the ticket states (section 10); G rows pass on both.

# ------------------------------------------------------------------------------------------------
# Fixtures.
# ------------------------------------------------------------------------------------------------

#' A protein object: each sample a common profile (SD 1 across features) plus its own noise.
k259Protein <- function(groups, noise_sd, n_features = 400L, seed = 259L, constant = FALSE) {
    samples <- sprintf("S%02d", seq_along(groups))
    noise_sd <- rep_len(noise_sd, length(groups))
    values <- withr::with_seed(seed, {
        profile <- stats::rnorm(n_features, 20, 1)
        vapply(seq_along(groups), function(i) profile + stats::rnorm(n_features, sd = noise_sd[[i]]), numeric(n_features))
    })
    if (constant) values[] <- 5
    colnames(values) <- samples
    table <- data.frame(Protein.Group = sprintf("PG%04d", seq_len(n_features)), values, check.names = FALSE)
    design <- data.frame(Run = samples, group = groups, batch = NA_character_,
        tech_rep_group = paste(groups, seq_along(groups), sep = "_"), row.names = samples, stringsAsFactors = FALSE)
    methods::new("ProteinQuantitativeData", protein_quant_table = table, protein_id_column = "Protein.Group",
        design_matrix = design, sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep_group")
}

k259Quiet <- function(code) suppressWarnings(suppressMessages(code))

#' The value of `code`, or the error it raised; the warnings it signalled, muffled.
k259Run <- function(code) {
    warnings <- list()
    value <- tryCatch(withCallingHandlers(code, warning = function(w) {
        warnings[[length(warnings) + 1L]] <<- w
        invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage")), error = function(e) e)
    list(value = value, warnings = warnings)
}

#' The labels of every text layer of a plot, white space collapsed (a label may be wrapped).
k259Labels <- function(plot) {
    labels <- unlist(lapply(seq_along(plot$layers), function(i) {
        data <- tryCatch(ggplot2::layer_data(plot, i), error = function(e) NULL)
        if (!is.null(data) && "label" %in% names(data)) as.character(data$label)
    }))
    gsub("[[:space:]]+", " ", labels)
}

#' The width and height in pixels from a PNG file's header.
k259PngSize <- function(path) {
    readBin(readBin(path, "raw", n = 24L)[17:24], "integer", n = 2L, size = 4L, endian = "big")
}

#' The legacy histogram (the code before KIT-259), built from the same data.
k259LegacyHistogram <- function(data) {
    breaks <- seq(min(round(data$pearson_correlation - 0.5, 2), na.rm = TRUE), 1, 0.001)
    ggplot2::ggplot(data, ggplot2::aes(.data[["pearson_correlation"]])) +
        ggplot2::geom_histogram(breaks = breaks) +
        ggplot2::scale_y_continuous(breaks = seq(0, 4, 1), limits = c(0, 4)) +
        ggplot2::xlab("Pearson Correlation") + ggplot2::ylab("Counts") +
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            panel.background = ggplot2::element_blank())
}

k259Png <- function(plot, dir, name) {
    path <- file.path(dir, name)
    k259Quiet(ggplot2::ggsave(path, plot, width = 10, height = 8, dpi = 150))
    path
}

# Two groups of five, noise SD 0.01: the 20 within-group correlations fall in one 0.001 bin.
k259Tight <- function() k259Protein(rep(c("A", "B"), each = 5L), 0.01)

k259ExpectAllDrawn <- function(plot, pairs) {
    ld <- ggplot2::layer_data(plot, 1L)
    expect_identical(as.integer(sum(ld$count[!is.na(ld$ymax)])), as.integer(pairs))
    grDevices::pdf(NULL)
    withr::defer(grDevices::dev.off())
    expect_no_warning(ggplot2::ggplotGrob(plot))
}

# ------------------------------------------------------------------------------------------------
# A1: every within-group pair is drawn.
# ------------------------------------------------------------------------------------------------

test_that("KIT-259 A-E1: the rebuilt proteomics histogram draws all 20 pairs of a bin holding more than 4", {
    plot <- k259Quiet(buildProtNormQcPlotProduct(k259Tight(), list(color_var = "group", shape_var = "group"),
        "correlation")$plot)
    expect_identical(max(ggplot2::layer_data(plot, 1L)$count), 20)
    k259ExpectAllDrawn(plot, 20L)
})

test_that("KIT-259 A-E2: the native plotPearson histogram draws all 20 pairs of a bin holding more than 4", {
    plot <- k259Quiet(plotPearson(k259Tight(), correlation_group = "group"))
    k259ExpectAllDrawn(plot, 20L)
})

test_that("KIT-259 A-G1: a histogram whose bins hold at most 4 draws ticks 0 to 4 and the legacy figure's bytes", {
    object <- k259Protein(rep(c("A", "B"), each = 3L), c(0.05, 0.2, 0.6, 0.1, 0.35, 0.9))
    dir <- withr::local_tempdir()
    rebuilt <- k259Quiet(buildProtNormQcPlotProduct(object, list(color_var = "group", shape_var = "group"),
        "correlation")$plot)
    native <- k259Quiet(plotPearson(object, correlation_group = "group"))
    for (case in list(list(name = "rebuilt", plot = rebuilt), list(name = "native", plot = native))) {
        legacy <- k259LegacyHistogram(case$plot$data)
        expect_lte(max(ggplot2::layer_data(case$plot, 1L)$count), 4, label = case$name)
        expect_equal(ggplot2::get_guide_data(case$plot, "y")$.value, 0:4, info = case$name)
        expect_identical(ggplot2::layer_data(case$plot, 1L), ggplot2::layer_data(legacy, 1L), info = case$name)
        expect_identical(unname(tools::md5sum(k259Png(case$plot, dir, paste0(case$name, ".png")))),
            unname(tools::md5sum(k259Png(legacy, dir, paste0(case$name, "-legacy.png")))), info = case$name)
    }
})

# ------------------------------------------------------------------------------------------------
# A2: no within-group pair with a finite correlation gives a labelled figure, not an error.
# ------------------------------------------------------------------------------------------------

test_that("KIT-259 A-E3: plotPearson with no within-group pair (an unset batch) returns a labelled empty figure", {
    got <- k259Run(plotPearson(k259Tight(), correlation_group = "batch"))
    expect_s3_class(got$value, "ggplot")
    if (inherits(got$value, "ggplot")) {
        expect_match(paste(k259Labels(got$value), collapse = " "), "No two samples share a value of 'batch'", fixed = TRUE)
        grDevices::pdf(NULL)
        withr::defer(grDevices::dev.off())
        expect_no_error(ggplot2::ggplotGrob(got$value))
    }
})

test_that("KIT-259 A-E4: the figure drawn from a saved reference with no within-group pair is labelled and its product saved", {
    object <- k259Tight()
    manager <- newDIAStateAdapter(newWorkflowEntityStore(withr::local_tempdir(), "kit259-a-e4", "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager))
    manager$.context$reference_states <- TRUE
    manager$saveState("protein_replicate_filtered", object, list(), "input")
    reference <- getDIAStateReference(manager)
    got <- testthat::with_mocked_bindings(k259Run(resolveDIAReferenceQcPlot(reference,
        list(color_var = "batch", shape_var = "group"), "correlation", buildProtNormQcPlotProduct, plotPca, plotRle,
        buildProtNormDensityPlot, buildProtNormCorrelationPlot)), getDIAActiveManager = function() manager)
    expect_s3_class(got$value, "ggplot")
    if (inherits(got$value, "ggplot")) {
        expect_match(paste(k259Labels(got$value), collapse = " "), "No two samples share a value of 'batch'", fixed = TRUE)
        expect_true(isSavedWorkflowReference(attr(got$value, "dia_product_reference")))
    }
})

test_that("KIT-259 A-E5: every within-group correlation NA gives a labelled empty figure, native and rebuilt", {
    object <- k259Protein(rep(c("A", "B"), each = 5L), 0.01, constant = TRUE)
    native <- k259Run(plotPearson(object, correlation_group = "group"))
    rebuilt <- k259Run(buildProtNormQcPlotProduct(object, list(color_var = "group", shape_var = "group"), "correlation"))
    expect_s3_class(native$value, "ggplot")
    expect_false(inherits(rebuilt$value, "error"))
    for (plot in list(native$value, if (is.list(rebuilt$value) && !inherits(rebuilt$value, "error")) rebuilt$value$plot)) {
        if (!inherits(plot, "ggplot")) next
        labels <- paste(k259Labels(plot), collapse = " ")
        expect_match(labels, "20 within-group sample pairs")
        expect_no_match(labels, "No two samples share")
    }
})

# ------------------------------------------------------------------------------------------------
# A3: lipidomics and metabolomics: negative correlations are drawn; the pairs with no correlation are counted.
# ------------------------------------------------------------------------------------------------

k259Omics <- function(omic, assay) {
    design <- data.frame(Run = c("S1", "S2", "S3", "S4"), group = c("A", "A", "B", "B"),
        tech_rep = c("t1", "t2", "t3", "t4"), stringsAsFactors = FALSE)
    if (identical(omic, "lipid")) {
        createLipidomicsAssayData(lipid_data = list(Assay1 = assay), design_matrix = design, lipid_id_column = "feature_id",
            sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep")
    } else {
        createMetaboliteAssayData(metabolite_data = list(Assay1 = assay), design_matrix = design,
            metabolite_id_column = "feature_id", sample_id = "Run", group_id = "group", technical_replicate_id = "tech_rep")
    }
}

k259Assay <- function(S1, S2, S3, S4) {
    data.frame(feature_id = c("F1", "F2", "F3", "F4"), S1 = S1, S2 = S2, S3 = S3, S4 = S4, stringsAsFactors = FALSE)
}

k259OmicsPlot <- function(omic, assay) {
    k259Quiet(plotPearson(k259Omics(omic, assay), correlation_group = "group"))$Assay1
}

#' The legacy lipid/metab histogram (the code before KIT-259), built from the same data.
k259LegacyOmicsHistogram <- function(data) {
    ggplot2::ggplot(data, ggplot2::aes(pearson_correlation)) +
        ggplot2::geom_histogram(breaks = seq(0, 1, 0.01), na.rm = TRUE) +
        ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
        ggplot2::xlab("Pearson Correlation") + ggplot2::ylab("Counts") + ggplot2::theme_bw() +
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())
}

for (omic in c("lipid", "metab")) {
    test_that(sprintf("KIT-259 A3-E1/E2 (%s): a negative within-group correlation is drawn", omic), {
        plot <- k259OmicsPlot(omic, k259Assay(c(1, 2, 3, 4), c(4, 3, 2, 1), c(1, 2, 4, 3), c(2, 4, 8, 6.5)))
        expect_s3_class(plot, "ggplot")
        expect_lt(min(plot$data$pearson_correlation), 0)
        expect_identical(sum(ggplot2::layer_data(plot, 1L)$count), 2)
        expect_lte(plot$scales$get_scales("x")$limits[[1L]], min(plot$data$pearson_correlation))
    })

    test_that(sprintf("KIT-259 A3-E4 (%s): a pair with no correlation is counted in a subtitle", omic), {
        plot <- k259OmicsPlot(omic, k259Assay(c(1, 2, 3, 4), c(1.1, 2.2, 2.9, 4.2), c(5, 5, 5, 5), c(2, 4, 8, 6.5)))
        expect_s3_class(plot, "ggplot")
        expect_identical(sum(is.na(plot$data$pearson_correlation)), 1L)
        expect_identical(plot$labels$subtitle,
            "1 sample pairs shown; 1 have no correlation (too few shared values, or a sample with no variation) and are not shown.")
    })

    test_that(sprintf("KIT-259 A3-G1 (%s): correlations at or above 0 and none missing draw as before, with no subtitle", omic), {
        plot <- k259OmicsPlot(omic, k259Assay(c(1, 2, 3, 4), c(1.1, 2.2, 2.9, 4.2), c(1, 2, 4, 3), c(2, 4, 8, 6.5)))
        expect_s3_class(plot, "ggplot")
        expect_identical(plot$scales$get_scales("x")$limits, c(0, 1))
        expect_equal(plot$scales$get_scales("x")$breaks, seq(0, 1, 0.1))
        expect_null(plot$labels$subtitle)
        expect_identical(ggplot2::layer_data(plot, 1L), ggplot2::layer_data(k259LegacyOmicsHistogram(plot$data), 1L))
    })
}

# ------------------------------------------------------------------------------------------------
# B2: a figure that cannot be drawn is replaced by a placeholder; the step completes with one warning notice.
# ------------------------------------------------------------------------------------------------

k259FailCorrelation <- function(state, env = parent.frame()) {
    original <- buildProtNormCorrelationPlot
    testthat::local_mocked_bindings(buildProtNormCorrelationPlot = function(...) {
        if (isTRUE(state$fail)) stop("injected drawing failure")
        original(...)
    }, .package = "MultiScholaR", .env = env)
}

test_that("KIT-259 B-E1: the post-normalisation generator replaces a figure it cannot draw by a placeholder", {
    dir <- withr::local_tempdir()
    got <- k259Run(generateProtNormPostNormalizationQcArtifacts(k259Tight(), dir, list(color_var = "group", shape_var = "group"),
        messageFn = function(...) invisible(NULL), gcFn = function() invisible(NULL),
        buildCorrelationFn = function(...) stop("injected drawing failure")))
    expect_false(inherits(got$value, "error"))
    expect_setequal(list.files(dir), c("post_norm_pca.png", "post_norm_rle.png", "post_norm_density.png",
        "post_norm_correlation.png"))
    unavailable <- Filter(function(w) inherits(w, "prot_norm_figure_unavailable"), got$warnings)
    expect_length(unavailable, 1L)
    if (length(unavailable) == 1L && file.exists(file.path(dir, "post_norm_correlation.png"))) {
        expect_identical(unavailable[[1L]]$figure, "post_norm_correlation.png")
        expect_match(unavailable[[1L]]$reason, "injected drawing failure", fixed = TRUE)
        expect_identical(k259PngSize(file.path(dir, "post_norm_correlation.png")), c(1500L, 1200L))
        writer <- get0("writeProtNormFigurePlaceholder", envir = asNamespace("MultiScholaR"), inherits = FALSE)
        expect_true(is.function(writer))
        if (is.function(writer)) {
            expected <- file.path(withr::local_tempdir(), "post_norm_correlation.png")
            k259Quiet(writer(expected, "post_norm_correlation.png", unavailable[[1L]]$reason, 10, 8))
            expect_identical(unname(tools::md5sum(file.path(dir, "post_norm_correlation.png"))),
                unname(tools::md5sum(expected)))
        }
    }
})

test_that("KIT-259 B-G2: a figure that could not be written and a store condition still stop the generator", {
    for (class in c("workflow_figure_unwritten", "saved_entity_corrupt")) {
        dir <- withr::local_tempdir()
        expect_error(k259Quiet(generateProtNormPostNormalizationQcArtifacts(k259Tight(), dir,
            list(color_var = "group", shape_var = "group"), messageFn = function(...) invisible(NULL),
            gcFn = function() invisible(NULL), buildCorrelationFn = function(...) rlang::abort("injected", class = class))),
            class = class)
        expect_false(file.exists(file.path(dir, "post_norm_correlation.png")))
    }
})

#' The tab on a native DIA owner (helper-prot-norm-core.R), through `steps` (each step and the flush after it); the
#' notices, the head, the completion flag and the post-normalisation correlation image's md5 after each.
k259Drive <- function(fixture, steps, managed = FALSE) {
    log <- new.env(parent = emptyenv())
    protNormW3Notices(log, fixture$root)
    probe <- new.env()
    if (managed) protNormW3Managed(fixture, probe)
    acc <- new.env(parent = emptyenv())
    acc$out <- list()
    protNormCoreQuiet(utils::capture.output(shiny::testServer(mod_prot_norm_server, args = list(
        workflow_data = fixture$workflow, experiment_paths = fixture$paths, omic_type = "proteomics",
        experiment_label = "core"), {
        norm <- session$getReturned()
        for (case in names(steps)) {
            log$notes <- character()
            steps[[case]](session)
            session$flushReact()
            correlation <- file.path(fixture$paths$protein_qc_dir, "post_norm_correlation.png")
            acc$out[[case]] <- list(notes = log$notes, head = fixture$manager$current_state,
                complete = isTRUE(shiny::isolate(norm$normalization_complete)),
                md5 = if (file.exists(correlation)) unname(tools::md5sum(correlation)) else NA_character_)
        }
    })))
    acc$out
}

k259Warning <- "^\\[warning [^]]+\\] The normalisation finished, but the figure 'Post-normalisation Pearson correlation' could not be drawn: injected drawing failure"
k259Completion <- paste("[message 10]", buildProtNormCompletionNotification("skip"))

k259FailingSteps <- function(state) {
    click <- protNormW3Clicker()
    run <- do.call(click, c(list("run_normalization"), protNormW3Inputs("skip")))
    apply <- click("apply_correlation_filter")
    list(
        fail = function(session) {
            state$fail <- TRUE
            run(session)
        },
        apply = function(session) {
            state$fail <- FALSE
            apply(session)
        })
}

test_that("KIT-259 B-E2, B-E6: in process, a figure that cannot be drawn completes the step with one warning; Apply then runs", {
    state <- new.env()
    k259FailCorrelation(state)
    out <- k259Drive(protNormCoreOwner(), k259FailingSteps(state))
    expect_true(out$fail$complete)
    expect_identical(out$fail$head, "normalized")
    expect_true(k259Completion %in% out$fail$notes)
    expect_identical(sum(grepl(k259Warning, out$fail$notes)), 1L)
    expect_false(any(grepl("Error in normalization", out$fail$notes, fixed = TRUE)))
    expect_false(any(grepl("outputs are pending", out$apply$notes, fixed = TRUE)))
    expect_true("[message 5] Correlation filtering completed! Ready for differential expression analysis." %in% out$apply$notes)
})

test_that("KIT-259 B-E3: managed, a stage figure failing on the worker completes the chain with one warning", {
    skip_if_not(workflowOwnerCapable())
    state <- new.env()
    k259FailCorrelation(state)
    out <- k259Drive(protNormCoreOwner(), k259FailingSteps(state)["fail"], managed = TRUE)
    expect_true(out$fail$complete)
    expect_identical(out$fail$head, "normalized")
    expect_true(k259Completion %in% out$fail$notes)
    expect_identical(sum(grepl(k259Warning, out$fail$notes)), 1L)
    expect_false(any(grepl("Error in normalization", out$fail$notes, fixed = TRUE)))
})

test_that("KIT-259 B-E4: a failing run replaces the earlier run's figure by its placeholder", {
    state <- new.env()
    k259FailCorrelation(state)
    click <- protNormW3Clicker()
    out <- k259Drive(protNormCoreOwner(), list(
        first = do.call(click, c(list("run_normalization"), protNormW3Inputs("skip"))),
        second = function(session) {
            state$fail <- TRUE
            click("run_normalization", norm_method = "cyclicloess")(session)
        }))
    expect_false(is.na(out$first$md5))
    expect_true(out$second$complete)
    expect_false(identical(out$second$md5, out$first$md5))
})

#' protNormCoreOwner (helper-prot-norm-core.R) with 32 samples in 16 groups of 2.
k259Owner16 <- function(envir = parent.frame()) {
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
    recipe <- readWorkflowCheckpointRecipe("small")
    recipe$n_samples <- 32L
    original <- makeWorkflowCheckpointProtein(recipe)
    groups <- sprintf("g%02d", rep(seq_len(16L), each = 2L))
    design <- original@design_matrix
    design$group <- groups
    design$tech_rep_group <- paste(groups, seq_along(groups), sep = "_")
    original@design_matrix <- design
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

test_that("KIT-259 B-E5: 16 groups with the default 'Shape by' give a PCA placeholder naming the shape limit; the step completes", {
    click <- protNormW3Clicker()
    fixture <- k259Owner16()
    out <- k259Drive(fixture, list(run = do.call(click, c(list("run_normalization"), protNormW3Inputs("skip")))))
    expect_true(out$run$complete)
    expect_identical(out$run$head, "normalized")
    expect_true(k259Completion %in% out$run$notes)
    expect_identical(sum(grepl(paste0("^\\[warning [^]]+\\] The normalisation finished, but the figure ",
        "'Post-normalisation PCA' could not be drawn: Insufficient values in manual scale"), out$run$notes)), 1L)
    pca <- file.path(fixture$paths$protein_qc_dir, "post_norm_pca.png")
    expect_true(file.exists(pca))
    if (file.exists(pca)) expect_identical(k259PngSize(pca), c(1200L, 900L))
})

test_that("KIT-259 B-G1: a run whose figures draw shows no warning and writes every figure", {
    click <- protNormW3Clicker()
    fixture <- protNormCoreOwner()
    out <- k259Drive(fixture, list(run = do.call(click, c(list("run_normalization"), protNormW3Inputs("skip")))))
    expect_true(out$run$complete)
    expect_identical(out$run$notes, k259Completion)
    expect_true(all(file.exists(file.path(fixture$paths$protein_qc_dir,
        paste0("post_norm_", c("pca", "rle", "density", "correlation"), ".png")))))
})
