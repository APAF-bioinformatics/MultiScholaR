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
    request <- file.path(time_dir, ".workflow-export", paste0("qvalue-", substr(reference$entity_id, 1L, 40L)))
    expect_true(all(nchar(file.path(request, "figures", names)) > 259L))
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
