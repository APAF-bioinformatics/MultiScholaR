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
            # KIT-182: each action's input entity, beside its operation.
            out$inputs[[case]] <- vapply(probe$actions[seq_along(probe$actions) > seen], function(action)
                action$input %||% NA_character_, "")
            out$generated[[case]] <- shiny::isolate(norm$pre_norm_qc_generated)
        }
    })))
    list(images = out$images, final = out$final, actions = out$actions, inputs = out$inputs, generated = out$generated,
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
    # KIT-122 (composite; re-frozen: the operations only): a normalisation's last figures action builds its composite
    # QC figure, after the stage figures.
    figures <- "normalization-qc-figures"
    expect_identical(managed$actions, list(tab_entry = figures, normalization_skip = c("normalization-ruv", figures, figures),
        normalization_manual = c("normalization-ruv", figures, figures, figures), aesthetics = c(figures, figures, figures),
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
    expect_length(unmanaged$images$pass_through, 8L)
    expect_length(unmanaged$images$ruv, 12L)
    expect_identical(managed$images, unmanaged$images)
    expect_identical(managed$final, unmanaged$final)
    expect_false(any(startsWith(names(unmanaged$counts), "worker")))
    expect_identical(names(managed$counts), c("worker saveWorkflowDiaQcAsset", "worker saveWorkflowDiaQcProduct"))
    expect_identical(unname(managed$counts), unname(unmanaged$counts))
    figures <- "normalization-qc-figures"
    # KIT-266 (declared edit): the pass-through also saves the pre figures' state and builds the combined figure on the worker
    # (three figures actions); the RUV step builds it after its stage figures.
    expect_identical(managed$actions, list(pass_through = c(figures, figures, figures),
        correlation = c("pearson-pairs", "pearson-selection", figures), ruv = c("protein_ruv_covariates", figures, figures)))
})

# ------------------------------------------------------------------------------------------------
# KIT-122 (composite): the normalisation's composite QC figure as a worker product. Assembled in the application, its
# rasters were the normalisation step's largest memory event there. On a managed session the worker builds it after the
# stage figures and before the publication half (the `composite` figures of `normalization-qc-figures`) and saves its
# three files; the publication half writes their bytes inside its journaled `produce`, where the builder wrote them.
# The tab is driven unmanaged and managed on the DIA fixture through the tab entry (the pre-normalisation images), a
# normalisation with RUV skipped (two columns), one with RUV applied (three columns and the canonical-correlation
# image) and that one pressed again. Counted, in the application and in the worker: the builder's calls. Recorded: the
# staged file set of each normalisation publication, and the published figure's files after each step (the image by
# its md5; the PDF and the RDS, which are not byte-stable between two builds, by their size).
# ------------------------------------------------------------------------------------------------

#' Drive the tab through `steps` (as `protNormK122Drive`), counting the composite's builder by process.
#' @return `builds` (the builder's calls: `application`, `worker`), `staged` (per normalisation publication, the files
#'   of each staged role), `figure` (after each step: the image's md5 and the three files' sizes) and the managed
#'   actions' operations per step.
protNormK122CompositeDrive <- function(fixture, steps, managed) {
    probe <- new.env()
    log <- new.env(parent = emptyenv())
    protNormW3Notices(log, fixture$root)
    if (managed) protNormW3Managed(fixture, probe)
    out <- new.env(parent = emptyenv())
    build <- generateProtNormCompositeQcFigure
    publish <- publishDIAActionOutputs
    testthat::local_mocked_bindings(
        generateProtNormCompositeQcFigure = function(...) {
            process <- if (identical(Sys.getenv("MULTISCHOLAR_WORKFLOW_WORKER"), "1")) "worker" else "application"
            out$builds[[process]] <- (out$builds[[process]] %||% 0L) + 1L
            build(...)
        },
        publishDIAActionOutputs = function(manager, reference, label, destinations, produce, identity = NULL) {
            publish(manager, reference, label, destinations, function(staged) {
                display <- produce(staged)
                if (identical(label, "normalization")) out$staged <- c(out$staged,
                    list(lapply(staged, function(directory) sort(list.files(directory, recursive = TRUE)))))
                display
            }, identity)
        }, .package = "MultiScholaR")
    selected <- shiny::reactiveVal("import")
    files <- file.path(fixture$paths$protein_qc_dir, paste0("proteomics_composite_QC_figure.", c("png", "pdf", "rds")))
    protNormCoreQuiet(utils::capture.output(shiny::testServer(mod_prot_norm_server, args = list(
        workflow_data = fixture$workflow, experiment_paths = fixture$paths, omic_type = "proteomics",
        experiment_label = "core", selected_tab = selected), {
        session$userData$dia_capacity_manager <- fixture$manager
        for (case in names(steps)) {
            seen <- length(probe$actions)
            steps[[case]](session, selected)
            for (i in 1:3) {
                session$elapse(100)
                session$flushReact()
            }
            out$figure[[case]] <- list(image = unname(tools::md5sum(files[[1L]])), bytes = file.size(files))
            out$actions[[case]] <- vapply(probe$actions[seq_along(probe$actions) > seen], `[[`, "", "operation")
        }
    })))
    list(builds = out$builds, staged = out$staged, figure = out$figure, actions = out$actions)
}

test_that("KIT-122 (composite): managed, the worker builds the composite QC figure and the application writes its saved files; the image is byte-identical to the unmanaged tab's (fails before)", {
    skip_if_not(workflowOwnerCapable())
    # The last step presses the button again with nothing changed.
    steps <- function() {
        dia <- protNormK122DiaSteps()
        c(dia[c("tab_entry", "normalization_skip", "normalization_manual")], list(normalization_repeat = dia$normalization_skip))
    }
    runs <- c("normalization_skip", "normalization_manual", "normalization_repeat")
    unmanaged <- protNormK122CompositeDrive(protNormCoreOwner(), steps(), managed = FALSE)
    managed <- protNormK122CompositeDrive(protNormCoreOwner(), steps(), managed = TRUE)
    # Unmanaged: the application builds the figure in each publication (today's path); no figure before the first, and
    # none for the repeated click, whose publication is already journaled.
    expect_identical(unmanaged$builds, list(application = 2L))
    expect_identical(unmanaged$figure$tab_entry, list(image = NA_character_, bytes = rep(NA_real_, 3L)))
    for (case in runs) {
        expect_false(is.na(unmanaged$figure[[case]]$image), info = case)
        expect_true(all(unmanaged$figure[[case]]$bytes > 0), info = case)
    }
    # RUV skipped and RUV applied are two layouts: two images.
    expect_false(identical(unmanaged$figure$normalization_skip$image, unmanaged$figure$normalization_manual$image))
    expect_identical(unmanaged$figure$normalization_repeat, unmanaged$figure$normalization_manual)
    # Managed: the worker builds each figure and the application builds none; for the repeated click the worker finds
    # the figure saved and builds nothing.
    expect_identical(managed$builds, list(worker = 2L))
    # The image is the unmanaged tab's, byte for byte, after every step; the PDF and the RDS are there and not empty.
    expect_identical(lapply(managed$figure, `[[`, "image"), lapply(unmanaged$figure, `[[`, "image"))
    for (case in runs) expect_true(all(managed$figure[[case]]$bytes > 0), info = case)
    # The staged file set of each normalisation publication is the same: the figure's three files are in the journal's
    # inventory on both paths, with the stage images and, when RUV is applied, the canonical-correlation image.
    expect_length(unmanaged$staged, 2L)
    expect_identical(managed$staged, unmanaged$staged)
    composite <- paste0("proteomics_composite_QC_figure.", c("pdf", "png", "rds"))
    for (publication in managed$staged) expect_true(all(composite %in% publication$protein_qc_dir))
    expect_false("ruv_corrected_cancor.png" %in% managed$staged[[1L]]$protein_qc_dir)
    expect_true("ruv_corrected_cancor.png" %in% managed$staged[[2L]]$protein_qc_dir)
    # The actions: the composite's after the stage figures of the states the publication draws.
    figures <- "normalization-qc-figures"
    expect_identical(managed$actions, list(tab_entry = figures, normalization_skip = c("normalization-ruv", figures, figures),
        normalization_manual = c("normalization-ruv", figures, figures, figures),
        normalization_repeat = c("normalization-ruv", figures, figures, figures)))
    for (case in names(unmanaged$actions)) expect_length(unmanaged$actions[[case]], 0L)
})

test_that("KIT-122 (composite): a figure is found saved only for the stage images it was built from", {
    manager <- list(.context = list(software = "test-build"))
    product <- list(entity_id = strrep("a", 64L))
    directory <- withr::local_tempdir()
    images <- c("pre_norm_pca.png", "post_norm_pca.png", "ruv_corrected_pca.png", "ruv_corrected_cancor.png")
    for (image in images) writeLines(image, file.path(directory, image))
    identities <- diaCompositeFigureIdentities(manager, product, "manual", directory)
    expect_identical(names(identities), c("png", "pdf", "rds"))
    expect_identical(anyDuplicated(identities), 0L)
    expect_identical(diaCompositeFigureIdentities(manager, product, "manual", directory), identities)
    # Another directory with the same images is the same figure; the canonical-correlation image is not read.
    copy <- withr::local_tempdir()
    file.copy(file.path(directory, images[1:3]), copy)
    writeLines("another rendering", file.path(copy, images[[4L]]))
    expect_identical(diaCompositeFigureIdentities(manager, product, "manual", copy), identities)
    # Another image (other aesthetics), a missing image, another layout or another run is another figure.
    writeLines("redrawn", file.path(copy, images[[2L]]))
    expect_false(identical(diaCompositeFigureIdentities(manager, product, "manual", copy)$png, identities$png))
    file.copy(file.path(directory, images[[2L]]), copy, overwrite = TRUE)
    expect_identical(diaCompositeFigureIdentities(manager, product, "manual", copy), identities)
    unlink(file.path(copy, images[[1L]]))
    expect_false(identical(diaCompositeFigureIdentities(manager, product, "manual", copy)$png, identities$png))
    expect_false(identical(diaCompositeFigureIdentities(manager, product, "skip", directory)$png, identities$png))
    expect_false(identical(diaCompositeFigureIdentities(manager, list(entity_id = strrep("b", 64L)), "manual", directory)$png,
        identities$png))
    # With no directory nothing is looked up or written: the caller builds, as before.
    expect_false(writeDIACompositeFigureFiles(manager, product, "manual", NULL, "proteomics"))
})

# ------------------------------------------------------------------------------------------------
# KIT-182 T8 (managed; F-e): the pre-normalisation figures follow the normalisation's source. The tab is driven
# unmanaged and managed through the tab entry (figures of PRF1), a normalisation, a Quality Control Revert and a QC
# change (the Protein Intensity Filter through its kernel and the step's display publication, as KIT-182's
# test-prot-norm-input-current.R makes it), and the normalisation again. On the fix the run reads the QC result, the
# managed actions start with that table's figures, and the pre-normalisation PCA is redrawn from it, byte-identical to
# the unmanaged tab's.
# ------------------------------------------------------------------------------------------------

protNormK182Steps <- function(seen) {
    # R1 instrument repair (VAL-54): the managed chain's read counter passes `fields = NULL` on, which reads no assay; the
    # scripted QC change in the application reads through the package's own function.
    read <- readWorkflowExecutionComponents
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
        qc_change = function(session, selected) {
            manager <- getDIAActiveManager()
            rlang::local_bindings(filtering_progress = methods::new("FilteringProgress"), .env = .GlobalEnv)
            grDevices::pdf(NULL)
            device <- grDevices::dev.cur()
            on.exit(if (device %in% grDevices::dev.list()) grDevices::dev.off(device), add = TRUE)
            testthat::with_mocked_bindings(withDIAStateScope(manager, function() {
                selectDIAStateReference(manager, "protein_replicate_filtered")
                value <- computeDIAReferenceProteinIntensityFilter(manager, getDIAStateReference(manager),
                    list(strict_mode = FALSE, min_reps_per_group = 2, min_groups = 2,
                        proteins_intensity_cutoff_percentile = 10))
                result <- finishDIAReferenceProteinStep(list2env(list(state_manager = manager, qc_params = list())),
                    value, publish = FALSE)
                publishDIAProteinProgress(result, "11_protein_intensity_filtered", "proteomics", "core", updateProteinFiltering)
            }), readWorkflowExecutionComponents = read, .package = "MultiScholaR")
            seen$q <- getDIAStateReference(manager)$entity_id
            seen$prf1 <- getDIAStateReference(manager, "protein_replicate_filtered")$entity_id
        },
        normalization_again = press("run_normalization"))
}

test_that("KIT-182 T8 (managed): after a QC change the run's actions start with the QC result's figures, and the pre-normalisation PCA is redrawn, byte-identical to the unmanaged tab's (fails before)", {
    skip_if_not(workflowOwnerCapable())
    unmanaged_seen <- new.env(parent = emptyenv())
    managed_seen <- new.env(parent = emptyenv())
    # The steps are built before the drive installs its mocks (they keep the package's own component read).
    unmanaged_steps <- protNormK182Steps(unmanaged_seen)
    managed_steps <- protNormK182Steps(managed_seen)
    unmanaged <- protNormK122Drive(protNormCoreOwner(), unmanaged_steps, managed = FALSE)
    managed <- protNormK122Drive(protNormCoreOwner(), managed_steps, managed = TRUE)
    figures <- "normalization-qc-figures"
    # The QC change removed proteins and is the head before the second run.
    expect_false(is.null(managed_seen$q))
    expect_false(identical(managed_seen$q, managed_seen$prf1))
    # Managed: the source's figures first, then the normalisation, both on the QC result.
    expect_identical(managed$actions$normalization_again, c(figures, "normalization-ruv", figures, figures))
    expect_identical(unname(managed$inputs$normalization_again[1:2]), rep(managed_seen$q, 2L))
    # The pre-normalisation PCA was redrawn after the QC change, on both paths, and is the same file on both.
    pca <- "pre_norm_pca.png"
    expect_false(identical(unmanaged$images$normalization_again[[pca]], unmanaged$images$tab_entry[[pca]]))
    expect_false(identical(managed$images$normalization_again[[pca]], managed$images$tab_entry[[pca]]))
    expect_identical(managed$images$normalization_again[[pca]], unmanaged$images$normalization_again[[pca]])
    # The steps before the QC change are as KIT-122 records them.
    expect_identical(managed$actions[c("tab_entry", "normalization_skip")],
        list(tab_entry = figures, normalization_skip = c("normalization-ruv", figures, figures)))
})
