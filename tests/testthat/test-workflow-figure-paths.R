# KIT-052: figures and over-long paths. The PNG device on Windows cannot write a path longer than 259
# characters and fails with a warning only; the doubles below reproduce that on any platform.

testthat::skip_if_not(workflowOwnerCapable(), "the saved-entity storage packages are required")

figurePathManager <- function(env = parent.frame()) {
    root <- withr::local_tempdir(.local_envir = env)
    manager <- newDIAStateAdapter(newWorkflowEntityStore(root, basename(root), "DIA"),
        WorkflowState$new(audit_enabled = FALSE))
    withr::defer(closeDIAStateAdapter(manager), envir = env)
    manager$.context$reference_states <- TRUE
    manager$saveState("protein_s4_created", makeWorkflowCheckpointProtein(readWorkflowCheckpointRecipe("small")),
        list(), "fixture")
    manager
}

figurePathNames <- function(step_name) {
    as.vector(outer(paste0(step_name, "_", c("proteins_total", "proteins_per_run", "peptides_total",
        "peptides_per_protein", "peptides_per_run", "combined_plots")), c("png", "pdf"), paste, sep = "."))
}

# A renderer that, like the Windows PNG device, leaves a file unwritten when `writes(path)` is FALSE.
figurePathRenderer <- function(writes) {
    function(step_name, omic_type, experiment_label, project_dirs, ...) {
        directory <- project_dirs[[paste(omic_type, experiment_label, sep = "_")]]$time_dir
        for (path in file.path(directory, figurePathNames(step_name))) {
            if (writes(path)) writeBin(charToRaw(basename(path)), path)
        }
        "grid"
    }
}

figurePathArguments <- function(time_dir) {
    list(data = NULL, step_name = "9_protein_s4_created", omic_type = "proteomics", experiment_label = "fixture",
        return_grid = TRUE, overwrite = TRUE, progress_owner = new.env(parent = emptyenv()),
        project_dirs = list(proteomics_fixture = list(time_dir = time_dir, publication_graphs_dir = time_dir)))
}

test_that("Q-value figures publish from a project whose staged paths exceed 259 characters", {
    manager <- figurePathManager()
    time_dir <- file.path(withr::local_tempdir(), strrep("a", 100L), strrep("b", 100L))
    dir.create(time_dir, recursive = TRUE)
    time_dir <- normalizePath(time_dir, winslash = "/")
    reference <- getDIAStateReference(manager)
    names <- figurePathNames("9_protein_s4_created")
    request <- file.path(time_dir, ".workflow-export", paste0("q-", substr(reference$entity_id, 1L, 16L)))
    expect_true(all(nchar(file.path(request, "f", names)) > 259L))
    scratch_before <- list.files(tempdir(), pattern = "^qvalue-")
    grid <- publishDIAQvalueOutput(manager, list(stateReference = reference), figurePathArguments(time_dir),
        figurePathRenderer(function(path) nchar(path) <= 259L))
    expect_identical(grid, "grid")
    expect_identical(vapply(file.path(time_dir, names), function(path) rawToChar(readBin(path, "raw", 200L)),
        character(1), USE.NAMES = FALSE), names)
    journal <- jsonlite::read_json(file.path(request, "journal.json"))
    expect_identical(journal$status, "complete")
    expect_identical(sort(names(journal$files)), sort(names))
    expect_identical(vapply(names, function(name) journal$files[[name]]$sha256, character(1), USE.NAMES = FALSE),
        vapply(file.path(time_dir, names), digest::digest, character(1), algo = "sha256", file = TRUE,
            USE.NAMES = FALSE))
    expect_true(manager$.context$snapshot$completed_outputs[[reference$entity_id]]$qc)
    expect_identical(list.files(tempdir(), pattern = "^qvalue-"), scratch_before)
})

test_that("a figure left unwritten by the renderer aborts the publication and is named", {
    manager <- figurePathManager()
    time_dir <- withr::local_tempdir()
    reference <- getDIAStateReference(manager)
    error <- expect_error(publishDIAQvalueOutput(manager, list(stateReference = reference),
        figurePathArguments(time_dir), figurePathRenderer(function(path) !grepl("[.]png$", path))),
        class = "saved_entity_outputs_pending")
    expect_s3_class(error, "workflow_figure_unwritten")
    text <- cli::ansi_strip(paste(conditionMessage(error), collapse = "\n"))
    expect_match(text, "required Q-value figures are incomplete", fixed = TRUE)
    for (name in grep("[.]png$", figurePathNames("9_protein_s4_created"), value = TRUE)) {
        expect_match(text, name, fixed = TRUE)
    }
    expect_no_match(text, "proteins_total.pdf", fixed = TRUE)
    expect_false(isTRUE(manager$.context$snapshot$completed_outputs[[reference$entity_id]]$qc))
})

test_that("the figure check names the missing files and, on Windows, the over-long path", {
    directory <- withr::local_tempdir()
    written <- file.path(directory, "written.pdf")
    writeLines("figure", written)
    expect_identical(checkWorkflowFigureFiles(written), written)
    long <- file.path(directory, strrep("c", 150L), paste0(strrep("d", 150L), ".png"))
    empty <- file.path(directory, "empty.png")
    file.create(empty)
    message_of <- function(...) {
        error <- expect_error(checkWorkflowFigureFiles(...), class = "workflow_figure_unwritten")
        gsub("\\s+", " ", cli::ansi_strip(paste(conditionMessage(error), collapse = " ")))
    }
    windows <- message_of(c(written, long, empty), windows = TRUE)
    expect_match(windows, basename(long), fixed = TRUE)
    expect_match(windows, "empty.png", fixed = TRUE)
    expect_no_match(windows, "written.pdf", fixed = TRUE)
    expect_match(windows, paste(nchar(long), "characters"), fixed = TRUE)
    expect_match(windows, "longer than 259 characters", fixed = TRUE)
    expect_no_match(message_of(c(written, long), windows = FALSE), "259", fixed = TRUE)
    expect_no_match(message_of(empty, windows = TRUE), "259", fixed = TRUE)
    expect_error(checkWorkflowFigureFiles(empty, class = "saved_entity_outputs_pending"),
        class = "saved_entity_outputs_pending")
})

test_that("the filtering renderer fails when the device warns and leaves a figure unwritten", {
    manager <- figurePathManager()
    summary <- withDIAStateScope(manager, function() readDIAProteinProgress(manager, getDIAStateReference(manager)))
    time_dir <- withr::local_tempdir()
    grDevices::pdf(NULL)
    device <- grDevices::dev.cur()
    withr::defer(if (device %in% grDevices::dev.list()) grDevices::dev.off(device))
    render <- function(unwritten) {
        testthat::local_mocked_bindings(ggsave = function(filename, ...) {
            if (grepl(unwritten, filename)) warning("agg could not write to the given file") else
                writeLines("figure", filename)
        }, .package = "MultiScholaR")
        suppressWarnings(suppressMessages(updateProteinFiltering(summary, "9_protein_s4_created", "proteomics",
            "fixture", overwrite = TRUE, return_grid = TRUE, progress_owner = new.env(parent = emptyenv()),
            project_dirs = list(proteomics_fixture = list(time_dir = time_dir, publication_graphs_dir = time_dir)))))
    }
    error <- expect_error(render("proteins_per_run[.]png$"), class = "workflow_figure_unwritten")
    expect_match(cli::ansi_strip(paste(conditionMessage(error), collapse = "\n")),
        "9_protein_s4_created_proteins_per_run.png", fixed = TRUE)
    expect_error(render("combined_plots[.]png$"), class = "workflow_figure_unwritten")
    expect_s3_class(render("^$"), "gtable")
})

# Decision M18: the staging folders have short names; one journaled under the earlier name is still found.

test_that("the staging folders have the short names, built by one helper", {
    root <- withr::local_tempdir()
    id <- strrep("0123456789abcdef", 4L)
    expect_identical(workflowStagedRequest(root, "qvalue", id), list(request_id = "q-0123456789abcdef", staged = "f"))
    expect_identical(workflowStagedRequest(root, "action", id, "da"), list(request_id = "a-0123456789abcdef", staged = "f"))
    earlier <- file.path(root, ".workflow-export", c(paste0("qvalue-", substr(id, 1L, 40L)), paste0("dia-da-", substr(id, 1L, 32L))))
    for (directory in earlier) dir.create(directory, recursive = TRUE)
    # An earlier folder without a journal belongs to no publication.
    expect_identical(workflowStagedRequest(root, "qvalue", id)$request_id, "q-0123456789abcdef")
    for (directory in earlier) writeLines("{}", file.path(directory, "journal.json"))
    expect_identical(workflowStagedRequest(root, "qvalue", id), list(request_id = basename(earlier[[1L]]), staged = "figures"))
    expect_identical(workflowStagedRequest(root, "action", id, "da"), list(request_id = basename(earlier[[2L]]), staged = "files"))
    expect_identical(workflowStagedRequest(root, "action", id, "correlation")$request_id, "a-0123456789abcdef")
})

test_that("Q-value figures stage under the short name and resume reads the same journal", {
    manager <- figurePathManager()
    time_dir <- normalizePath(withr::local_tempdir(), winslash = "/")
    reference <- getDIAStateReference(manager)
    manager$.context$experiment_paths <- list(time_dir = time_dir)
    manager$.context$snapshot$output_products[[reference$entity_id]] <- "product"
    record <- list(data = list(products = list()))
    expect_true(resumeDIAQvalueOutput(manager, record, hydrate = FALSE))
    seen <- NULL
    renderer <- figurePathRenderer(function(path) TRUE)
    publishDIAQvalueOutput(manager, list(stateReference = reference), figurePathArguments(time_dir),
        function(...) {
            seen <<- list.files(file.path(time_dir, ".workflow-export"))
            renderer(...)
        })
    expect_identical(seen, paste0("q-", substr(reference$entity_id, 1L, 16L)))
    expect_identical(list.files(file.path(time_dir, ".workflow-export")), seen)
    expect_true(dir.exists(file.path(time_dir, ".workflow-export", seen, "f")))
    expect_identical(jsonlite::read_json(file.path(time_dir, ".workflow-export", seen, "journal.json"))$status, "complete")
    # The acknowledgement aside, the journal alone says the publication is complete.
    manager$.context$snapshot$completed_outputs[[reference$entity_id]]$qc <- NULL
    expect_null(resumeDIAQvalueOutput(manager, record, hydrate = FALSE))
})

test_that("a Q-value publication left pending under the earlier name is found and finished there", {
    manager <- figurePathManager()
    time_dir <- normalizePath(withr::local_tempdir(), winslash = "/")
    reference <- getDIAStateReference(manager)
    manager$.context$experiment_paths <- list(time_dir = time_dir)
    manager$.context$snapshot$output_products[[reference$entity_id]] <- "product"
    names <- figurePathNames("9_protein_s4_created")
    earlier <- file.path(time_dir, ".workflow-export", paste0("qvalue-", substr(reference$entity_id, 1L, 40L)))
    dir.create(file.path(earlier, "figures"), recursive = TRUE)
    for (name in names) writeBin(charToRaw(name), file.path(earlier, "figures", name))
    saveRDS(list(grid = "earlier grid", progress = methods::new("FilteringProgress")), file.path(earlier, "display.rds"))
    jsonlite::write_json(list(schema = "multischolar.dia_qc_publication.v1", status = "outputs_pending",
        request_id = basename(earlier), state_id = reference$entity_id,
        incarnation = manager$.context$store$identity$incarnation,
        files = stats::setNames(lapply(names, function(name) list(
            sha256 = digest::digest(file = file.path(earlier, "figures", name), algo = "sha256"), bytes = nchar(name))), names),
        display_sha256 = digest::digest(file = file.path(earlier, "display.rds"), algo = "sha256"),
        scientific_checkpoint_count = 0L), file.path(earlier, "journal.json"), auto_unbox = TRUE, null = "null")
    expect_true(resumeDIAQvalueOutput(manager, list(data = list(products = list())), hydrate = FALSE))
    grid <- publishDIAQvalueOutput(manager, list(stateReference = reference), figurePathArguments(time_dir),
        function(...) stop("render repeated"))
    expect_identical(grid, "earlier grid")
    expect_identical(list.files(file.path(time_dir, ".workflow-export")), basename(earlier))
    expect_identical(jsonlite::read_json(file.path(earlier, "journal.json"))$status, "complete")
    expect_identical(vapply(file.path(time_dir, names), function(path) rawToChar(readBin(path, "raw", 200L)),
        character(1), USE.NAMES = FALSE), names)
    expect_true(manager$.context$snapshot$completed_outputs[[reference$entity_id]]$qc)
})

test_that("action outputs stage under the short name; one pending under the earlier name is finished there", {
    manager <- figurePathManager()
    reference <- getDIAStateReference(manager)
    destination <- normalizePath(withr::local_tempdir(), winslash = "/")
    key <- digest::digest(list(state = reference$entity_id, action = "protein-qc", identity = NULL), algo = "sha256")
    staged_at <- NULL
    publishWorkflowActionOutputs(manager, reference, "protein-qc", list(protein_qc_dir = destination), function(staged) {
        staged_at <<- staged$protein_qc_dir
        writeLines("new", file.path(staged$protein_qc_dir, "a.tsv"))
        list(completed = TRUE)
    })
    short <- paste0("a-", substr(key, 1L, 16L))
    expect_identical(staged_at, file.path(destination, ".workflow-export", short, "f", "protein_qc_dir"))
    expect_identical(manager$.context$snapshot$action_outputs[[reference$entity_id]]$request_id, short)
    expect_identical(readLines(file.path(destination, "a.tsv")), "new")
    expect_null(assertWorkflowExportReady(manager))

    earlier_root <- normalizePath(withr::local_tempdir(), winslash = "/")
    earlier <- file.path(earlier_root, ".workflow-export", paste0("dia-protein-qc-", substr(key, 1L, 32L)))
    dir.create(file.path(earlier, "files", "protein_qc_dir"), recursive = TRUE)
    writeLines("earlier", file.path(earlier, "files", "protein_qc_dir", "a.tsv"))
    saveRDS(list(completed = "earlier"), file.path(earlier, "display.rds"), version = 3)
    jsonlite::write_json(list(schema = "multischolar.dia_qc_publication.v1", status = "staged",
        request_id = basename(earlier), state_id = reference$entity_id,
        incarnation = manager$.context$store$identity$incarnation,
        files = list(`protein_qc_dir/a.tsv` = list(sha256 = digest::digest(
            file = file.path(earlier, "files", "protein_qc_dir", "a.tsv"), algo = "sha256"),
            bytes = file.size(file.path(earlier, "files", "protein_qc_dir", "a.tsv")))),
        display_sha256 = digest::digest(file = file.path(earlier, "display.rds"), algo = "sha256")),
        file.path(earlier, "journal.json"), auto_unbox = TRUE, null = "null")
    display <- publishWorkflowActionOutputs(manager, reference, "protein-qc", list(protein_qc_dir = earlier_root),
        function(staged) stop("produce repeated"))
    expect_identical(display, list(completed = "earlier"))
    expect_identical(readLines(file.path(earlier_root, "a.tsv")), "earlier")
    expect_identical(list.files(file.path(earlier_root, ".workflow-export")), basename(earlier))
    expect_identical(manager$.context$snapshot$action_outputs[[reference$entity_id]]$request_id, basename(earlier))
    expect_identical(jsonlite::read_json(file.path(earlier, "journal.json"))$status, "complete")
    expect_null(assertWorkflowExportReady(manager))
})

# The worker's request folder (KIT-134 run 44241609555, defect #5): `<root>/scripts/<omic>/.dia-native-<14>/.workflow-entities/
# .workflow-export/wfx-<32>/{da,omics}-publication/<role>/<file>`. Its figures are drawn in a short folder and copied in; a
# copy that fails (Windows with long paths disabled) is named. The device double below refuses, silently, as ggsave's PNG
# device (ragg) on Windows does, any path of 260 characters or more.

figurePathLongDevice <- function(env = parent.frame()) {
    drawn <- new.env(parent = emptyenv())
    drawn$paths <- character()
    original <- ggplot2::ggsave
    testthat::local_mocked_bindings(ggsave = function(filename, ...) {
        drawn$paths <- c(drawn$paths, filename)
        if (nchar(filename) < 260L) original(filename, ...)
        invisible(filename)
    }, .package = "ggplot2", .env = env)
    drawn
}

# A folder no file can be copied into, as a path over 259 characters is on Windows with long paths disabled.
figurePathReadOnly <- function(path, env = parent.frame()) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    Sys.chmod(path, "0555")
    withr::defer(Sys.chmod(path, "0755"), envir = env)
    testthat::skip_if(file.access(path, 2L) == 0L, "the folder stays writable (root?)")
    path
}

figurePathText <- function(error) gsub("\\s+", " ", cli::ansi_strip(paste(conditionMessage(error), collapse = " ")))

test_that("the omics worker draws its figures in a short folder: from a long project root none is lost", {
    testthat::local_mocked_bindings(workflowExecutionBuildIdentity = function() strrep("f", 64),
        finishWorkflowExecutionJob = function(owner, status) {
            owner$job <- NULL
            list(status = status)
        }, .package = "MultiScholaR")
    fixture <- lipidDeclarationsOwner(file.path(withr::local_tempdir(), strrep("p", 60L), strrep("q", 60L)))
    withr::defer(closeDIAStateAdapter(fixture$owner))
    suppressWarnings(suppressMessages(runWorkflowOmicsReferenceDuplicateResolution(lipidomicsDeclarationBinding, fixture$owner, list())))
    controller <- omicsWorkerController(fixture)
    drawn <- figurePathLongDevice()
    scratch_before <- list.files(tempdir(), pattern = "^wfx-")
    done <- omicsWorkerAction(fixture, controller, "lipid_qc_figures",
        list(stage = "post_filter", grouping_variable = "group", shape_variable = NULL))
    expect_identical(done$acceptance, "complete")
    names <- basename(drawn$paths)
    expect_gt(length(names), 0L)
    # Drawn in the request folder, every one of these paths would be 260 characters or more.
    request <- file.path(done$paths$request_dir, "omics-publication", "lipid_qc_dir", names)
    expect_true(all(nchar(request) >= 260L))
    expect_true(all(nchar(drawn$paths) < 260L))
    expect_true(all(file.exists(request)))
    expect_setequal(list.files(fixture$paths$lipid_qc_dir, "[.]png$"), names)
    expect_identical(list.files(tempdir(), pattern = "^wfx-"), scratch_before)
})

figurePathDaWorker <- function(request_dir, draw, env = parent.frame()) {
    head <- list(entity_id = "head")
    testthat::local_mocked_bindings(
        getWorkflowStateReference = function(manager) head,
        resolveWorkflowStateValue = function(context, reference) list(source = head, parameters = list(qvalue = 0.05)),
        readWorkflowStateMetadata = function(manager, reference) list(data = "checkpoint"),
        isDIAResultManifest = function(x) TRUE,
        withWorkflowStateScope = function(manager, fn) fn(),
        produceDIADaPublication = function(manager, manifest, annotations, qvalue, staged, save_assets = TRUE) {
            for (name in c("Volcano_Plots/groupB-groupA.png", "Volcano_Plots/groupB-groupA.pdf",
                "NumSigDaMolecules/da_proteins_num_sig_da_molecules.png")) {
                path <- file.path(staged$publication_graphs_dir, name)
                dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
                draw(path)
            }
            writeLines("table", file.path(staged$da_output_dir, "da_proteins_groupB-groupA_long_annot.tsv"))
            list(completed = TRUE)
        }, .package = "MultiScholaR", .env = env)
    manager <- new.env()
    manager$.context <- new.env()
    manager$.context$snapshot <- list(da_results = list(entity_id = "manifest"))
    manager$current_state <- "da"
    request <- list(operation = "da_publication", input_references = list(state = head), parameters = list(),
        directory = request_dir, worker_cap = 12 * 1024^3)
    function() computeWorkflowDaPublicationAction(manager, request, "whole")
}

test_that("the proteomics worker's DA publication draws in a short folder and copies into the request folder", {
    request_dir <- file.path(withr::local_tempdir(), strrep("p", 100L), strrep("q", 100L))
    dir.create(request_dir, recursive = TRUE)
    drawn <- character()
    run <- figurePathDaWorker(request_dir, function(path) {
        drawn <<- c(drawn, path)
        if (nchar(path) < 260L) writeBin(charToRaw(basename(path)), path)
    })
    value <- run()
    expected <- c("da_output_dir/da_proteins_groupB-groupA_long_annot.tsv", "publication_graphs_dir/NumSigDaMolecules/da_proteins_num_sig_da_molecules.png",
        "publication_graphs_dir/Volcano_Plots/groupB-groupA.pdf", "publication_graphs_dir/Volcano_Plots/groupB-groupA.png")
    expect_setequal(value$files, expected)
    published <- file.path(request_dir, "da-publication", expected[-1L])
    expect_true(all(nchar(published) >= 260L))
    expect_true(all(nchar(drawn) < 260L))
    expect_identical(vapply(published, function(path) rawToChar(readBin(path, "raw", 200L)), character(1), USE.NAMES = FALSE),
        basename(published))
})

test_that("a copy into the worker's request folder that fails stops the step and names the files", {
    request_dir <- withr::local_tempdir()
    figurePathReadOnly(file.path(request_dir, "da-publication", "publication_graphs_dir", "Volcano_Plots"))
    run <- figurePathDaWorker(request_dir, function(path) writeBin(charToRaw(basename(path)), path))
    error <- expect_error(suppressWarnings(run()), class = "workflow_figure_unwritten")
    expect_s3_class(error, "saved_entity_outputs_pending")
    text <- figurePathText(error)
    expect_match(text, "could not be copied into its request folder", fixed = TRUE)
    expect_match(text, "groupB-groupA.png", fixed = TRUE)
    expect_match(text, "groupB-groupA.pdf", fixed = TRUE)
    expect_no_match(text, "num_sig", fixed = TRUE)
})

test_that("a copy that fails is named in process too: the omics publication and the action outputs", {
    root <- withr::local_tempdir()
    testthat::local_mocked_bindings(
        computeWorkflowTrxQcFiguresAction = function(manager, request, context) {
            target <- file.path(request$directory, "omics-publication", "gene_qc_dir")
            dir.create(target, recursive = TRUE)
            writeLines("figure", file.path(target, "trx_raw_pca.png"))
            list(display = list(saved = TRUE))
        },
        getWorkflowStateReference = function(manager) list(entity_id = "head"),
        withWorkflowStateScope = function(manager, fn) fn(),
        .package = "MultiScholaR")
    target <- figurePathReadOnly(file.path(root, "gene_qc"))
    error <- expect_error(suppressWarnings(runWorkflowOmicsPublicationInProcess(transcriptomicsDeclarationBinding, "owner",
        "trx_qc_figures", list(stage = "raw"), list(gene_qc_dir = target))), class = "workflow_figure_unwritten")
    expect_match(figurePathText(error), "trx_raw_pca.png", fixed = TRUE)

    manager <- figurePathManager()
    reference <- getDIAStateReference(manager)
    project <- normalizePath(withr::local_tempdir(), winslash = "/")
    manager$.context$experiment_paths <- list(source_dir = project)
    destination <- figurePathReadOnly(file.path(project, "qc"))
    error <- expect_error(suppressWarnings(publishWorkflowActionOutputs(manager, reference, "protein-qc",
        list(protein_qc_dir = destination), function(staged) {
            writeLines("new", file.path(staged$protein_qc_dir, "a.png"))
            list(completed = TRUE)
        }, identity = "read-only")), class = "workflow_figure_unwritten")
    expect_match(figurePathText(error), "cannot stage action output publication of protein_qc_dir/a.png", fixed = TRUE)
    expect_false(isTRUE(manager$.context$snapshot$completed_outputs[[reference$entity_id]][["protein-qc"]]))
})

test_that("a failed copy says, on Windows, that the project's path is too long with long paths disabled", {
    directory <- withr::local_tempdir()
    long <- file.path(directory, strrep("c", 150L), paste0(strrep("d", 150L), ".png"))
    copied <- file.path(directory, "copied.png")
    file.create(copied)
    message_of <- function(...) figurePathText(expect_error(checkWorkflowFigureFiles(...), class = "workflow_figure_unwritten"))
    windows <- message_of(c(copied, long), "files could not be copied", windows = TRUE, copied = c(TRUE, FALSE))
    expect_match(windows, basename(long), fixed = TRUE)
    # A copy may carry an empty file: it is not named.
    expect_no_match(windows, "copied.png", fixed = TRUE)
    expect_match(windows, paste(nchar(long), "characters"), fixed = TRUE)
    expect_match(windows, "too long for Windows with long paths disabled", fixed = TRUE)
    expect_match(windows, "enable Windows long paths", fixed = TRUE)
    expect_no_match(message_of(long, windows = FALSE, copied = FALSE), "259", fixed = TRUE)
    # A copy that reported success but left no file is named as well.
    expect_match(message_of(file.path(directory, "gone.png"), copied = TRUE), "gone.png", fixed = TRUE)
})

# KIT-134 run 44275615462: the in-process DA publication drew in place in the action staging folder; at a project root of
# 125 characters `da_proteins_num_sig_da_molecules.png` (260) was missing and nothing said so.

test_that("the in-process DA publication draws in a short folder: from a long staging path none is lost", {
    manager <- figurePathManager()
    root <- normalizePath(withr::local_tempdir(), winslash = "/")
    project <- file.path(root, strrep("p", 150L - nchar(root)))
    paths <- list(da_output_dir = file.path(project, "da"), publication_graphs_dir = file.path(project, "graphs"))
    for (path in paths) dir.create(path, recursive = TRUE)
    names <- c("Volcano_Plots/groupB-groupA.png", "Volcano_Plots/groupB-groupA.pdf",
        "NumSigDaMolecules/da_proteins_num_sig_da_molecules.png")
    drawn <- character()
    manifest <- list(contrasts = data.frame(friendly_names = "groupB-groupA"))
    testthat::local_mocked_bindings(
        diaDaRequest = function(...) list(source = NULL, steps = list(), manifest = manifest, annotations = data.frame()),
        acceptDIADaManifest = function(manager, manifest, results) manifest,
        produceDIADaPublication = function(manager, manifest, annotations, qvalue, staged, save_assets = TRUE) {
            for (name in names) {
                path <- file.path(staged$publication_graphs_dir, name)
                dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
                drawn <<- c(drawn, path)
                # As ggsave's PNG device (ragg) on Windows: nothing written, nothing said, from 260 characters.
                if (nchar(path) < 260L) writeBin(charToRaw(basename(path)), path)
            }
            writeLines("table", file.path(staged$da_output_dir, "da_proteins_groupB-groupA_long_annot.tsv"))
            list(completed = TRUE)
        }, .package = "MultiScholaR")
    da_data <- new.env()
    workflow <- new.env()
    workflow$state_manager <- manager
    workflow$tab_status <- list()
    runDIAReferenceDaAnalysis(list(da_q_val_thresh = 0.05), da_data, workflow, paths)
    reference <- getDIAStateReference(manager)
    key <- digest::digest(list(state = reference$entity_id, action = "da",
        identity = diaDaPublicationIdentity(manifest, NULL)), algo = "sha256")
    # Drawn in the action staging folder, the number-significant figure's path would be 260 characters or more.
    expect_gte(nchar(file.path(paths$da_output_dir, ".workflow-export", paste0("a-", substr(key, 1L, 16L)), "f",
        "publication_graphs_dir", names[[3L]])), 260L)
    expect_true(all(nchar(drawn) < 260L))
    published <- file.path(paths$publication_graphs_dir, names)
    expect_identical(vapply(published, function(path) if (file.exists(path)) rawToChar(readBin(path, "raw", 200L)) else NA_character_,
        character(1), USE.NAMES = FALSE), basename(published))
    expect_true(file.exists(file.path(paths$da_output_dir, "da_proteins_groupB-groupA_long_annot.tsv")))
    expect_true(manager$.context$snapshot$completed_outputs[[reference$entity_id]]$da)
    expect_identical(da_data$da_results_list, manifest)
})
