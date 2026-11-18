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
