# KIT-122: the proteomics normalisation tab's QC figures as a worker product, driven through the tab unmanaged and
# managed. The first two blocks moved here unchanged from test-prot-norm-core.R (which had reached the suite's
# per-file ceiling); the fixtures they share with that file are in helper-prot-norm-core.R.

# ------------------------------------------------------------------------------------------------
# KIT-122 (audit E1, F3): the tab's QC figures as a worker product. On a managed session the worker saves the QC plot
# products and the stage images (`normalization-qc-figures`) before the tab's generators and publication halves run, at
# every place the tab draws them: tab entry, the normalisation's publication, an aesthetics change, the correlation's
# publication (applied and skipped) and, on DIA_LIMPA, the pass-through and the RUV covariates. The tab is driven
# unmanaged (it computes and draws in the application, as before) and managed, on the same fixture. Compared after each
# step: the md5 of every stage image in the QC directory (the files are byte-identical) and the correlation tab's PCA.
# Counted, in the application and in the worker (`MULTISCHOLAR_WORKFLOW_WORKER`, set around the worker's entry point):
# the saves of a QC plot product (one per product computed from a hydrated state) and of a rendered image (one per
# image drawn).
# ------------------------------------------------------------------------------------------------

#' The md5 of every stage image in the fixture's protein QC directory.
protNormK122Images <- function(fixture) {
    directory <- fixture$paths$protein_qc_dir
    files <- sort(list.files(directory, pattern = "^(pre_norm|post_norm|ruv_corrected)_(pca|rle|density|correlation)[.]png$"))
    stats::setNames(unname(tools::md5sum(file.path(directory, files))), files)
}

#' Drive the tab through `steps` (named `function(session, selected)`; `selected` is the workflow's selected tab).
#' @return `images` and `final` (the correlation tab's PCA: its data and title) after each step, the `counts`, and the
#'   managed actions' operations per step.
protNormK122Drive <- function(fixture, steps, managed) {
    probe <- new.env()
    log <- new.env(parent = emptyenv())
    protNormW3Notices(log, fixture$root)
    if (managed) protNormW3Managed(fixture, probe)
    counts <- new.env(parent = emptyenv())
    counted <- function(name) {
        original <- get(name, envir = asNamespace("MultiScholaR"))
        function(...) {
            key <- paste(if (identical(Sys.getenv("MULTISCHOLAR_WORKFLOW_WORKER"), "1")) "worker" else "application", name)
            counts[[key]] <- (counts[[key]] %||% 0L) + 1L
            original(...)
        }
    }
    testthat::local_mocked_bindings(saveWorkflowDiaQcProduct = counted("saveWorkflowDiaQcProduct"),
        saveWorkflowDiaQcAsset = counted("saveWorkflowDiaQcAsset"), .package = "MultiScholaR")
    selected <- shiny::reactiveVal("import")
    out <- new.env(parent = emptyenv())
    protNormCoreQuiet(utils::capture.output(shiny::testServer(mod_prot_norm_server, args = list(
        workflow_data = fixture$workflow, experiment_paths = fixture$paths, omic_type = "proteomics",
        experiment_label = "core", selected_tab = selected), {
        # The application's session names its owner (the generators resolve it outside a publication's scope).
        session$userData$dia_capacity_manager <- fixture$manager
        norm <- session$getReturned()
        for (case in names(steps)) {
            seen <- length(probe$actions)
            steps[[case]](session, selected)
            for (i in 1:3) {
                session$elapse(100)
                session$flushReact()
            }
            out$images[[case]] <- protNormK122Images(fixture)
            plot <- shiny::isolate(norm$final_qc_plot)
            out$final[[case]] <- if (!is.null(plot)) list(data = protNormW1Sha(as.data.frame(plot$data)), title = plot$labels$title)
            out$actions[[case]] <- vapply(probe$actions[seq_along(probe$actions) > seen], `[[`, "", "operation")
            out$generated[[case]] <- shiny::isolate(norm$pre_norm_qc_generated)
        }
    })))
    list(images = out$images, final = out$final, actions = out$actions, generated = out$generated,
        counts = unlist(mget(sort(ls(counts)), envir = counts)), notes = log$notes)
}

#' The DIA route's steps: tab entry, two normalisations, an aesthetics change, the correlation applied and skipped.
protNormK122DiaSteps <- function() {
    click <- protNormW3Clicker()
    press <- function(...) {
        step <- click(...)
        function(session, selected) step(session)
    }
    list(
        tab_entry = function(session, selected) {
            do.call(session$setInputs, protNormW3Inputs("skip"))
            selected("normalization")
        },
        normalization_skip = press("run_normalization"),
        normalization_manual = press("run_normalization", ruv_mode = "manual"),
        aesthetics = function(session, selected) session$setInputs(color_variable = "batch"),
        correlation_apply = press("apply_correlation_filter"),
        correlation_skip = press("skip_correlation_filter"))
}

test_that("KIT-122 (E1): managed, the worker saves the tab's QC figures and the application draws none; the images are byte-identical to the unmanaged tab's (fails before)", {
    skip_if_not(workflowOwnerCapable())
    unmanaged <- protNormK122Drive(protNormCoreOwner(), protNormK122DiaSteps(), managed = FALSE)
    managed <- protNormK122Drive(protNormCoreOwner(), protNormK122DiaSteps(), managed = TRUE)
    # Unmanaged: every figure is computed and drawn in the application, at each of the tab's places (today's path).
    expect_identical(names(unmanaged$images$tab_entry), paste0("pre_norm_", c("correlation", "density", "pca", "rle"), ".png"))
    expect_length(unmanaged$images$normalization_skip, 8L)
    expect_length(unmanaged$images$normalization_manual, 12L)
    expect_false(identical(unmanaged$images$aesthetics, unmanaged$images$normalization_manual))
    expect_true(all(unmanaged$counts[c("application saveWorkflowDiaQcAsset", "application saveWorkflowDiaQcProduct")] > 0L))
    expect_false(any(startsWith(names(unmanaged$counts), "worker")))
    for (case in names(unmanaged$actions)) expect_length(unmanaged$actions[[case]], 0L)
    # Managed: the same files, byte for byte, after every step; the same PCA on the correlation tab.
    expect_identical(managed$images, unmanaged$images)
    expect_identical(managed$final, unmanaged$final)
    expect_identical(managed$generated, unmanaged$generated)
    # The worker computed every product and drew every image; the application saved none (it found each one saved).
    expect_identical(names(managed$counts), c("worker saveWorkflowDiaQcAsset", "worker saveWorkflowDiaQcProduct"))
    expect_identical(managed$counts[["worker saveWorkflowDiaQcAsset"]], unmanaged$counts[["application saveWorkflowDiaQcAsset"]])
    expect_identical(managed$counts[["worker saveWorkflowDiaQcProduct"]], unmanaged$counts[["application saveWorkflowDiaQcProduct"]])
    # The actions: each place's figures after its science and before its publication; the skipped correlation's PCA is
    # the normalisation's final's, which the worker saved with that step, so it dispatches nothing.
    figures <- "normalization-qc-figures"
    expect_identical(managed$actions, list(tab_entry = figures, normalization_skip = c("normalization-ruv", figures),
        normalization_manual = c("normalization-ruv", figures, figures), aesthetics = c(figures, figures, figures),
        correlation_apply = c("pearson-pairs", "pearson-selection", figures), correlation_skip = character()))
})

test_that("KIT-122 (E1): on DIA_LIMPA the pass-through's and the RUV covariates' figures are the worker's, byte-identical (fails before)", {
    skip_if_not(workflowOwnerCapable())
    skip_if_not_installed("limpa")
    skip_if_not_installed("mixOmics")
    steps <- function() {
        click <- protNormW3Clicker()
        press <- function(...) {
            step <- click(...)
            function(session, selected) step(session)
        }
        list(pass_through = do.call(press, c(list("run_normalization"), protNormW3Inputs("covariates"))),
            correlation = press("apply_correlation_filter"), ruv = press("apply_ruv_covariates"))
    }
    unmanaged <- protNormK122Drive(protNormW1LimpaOwner(), steps(), managed = FALSE)
    managed <- protNormK122Drive(protNormW1LimpaOwner(), steps(), managed = TRUE)
    expect_length(unmanaged$images$pass_through, 4L)
    expect_length(unmanaged$images$ruv, 8L)
    expect_identical(managed$images, unmanaged$images)
    expect_identical(managed$final, unmanaged$final)
    expect_false(any(startsWith(names(unmanaged$counts), "worker")))
    expect_identical(names(managed$counts), c("worker saveWorkflowDiaQcAsset", "worker saveWorkflowDiaQcProduct"))
    expect_identical(unname(managed$counts), unname(unmanaged$counts))
    figures <- "normalization-qc-figures"
    expect_identical(managed$actions, list(pass_through = figures, correlation = c("pearson-pairs", "pearson-selection", figures),
        ruv = c("protein_ruv_covariates", figures)))
})
