# KIT-234 (VAL-209): no R Markdown workbook ships, and the application copies none into a project. Declared before the
# run (ticket section 6, revision 1.2; each case's failing value on the base is stated beside it there):
#   W1  system.file("workbooks") of the source tree under load_all is "" (base: the folder of the eight workbooks);
#   W2  a new proteomics and metabolomics project holds no .rmd under scripts/<omic>, and both folders exist (base: the
#       five DIA and the two metabolomics workbooks);
#   W3  "Use Existing" leaves a user's edited copy and own file byte for byte as they were (base: the edited copy is
#       overwritten by the shipped workbook, KIT-231 G0);
#   W4  with the integration tab on, scripts/integration exists and holds no .rmd (base: the integration workbook);
#   W5  (guard) the Integration tab's session store opens on W4's project (passes on the base too);
#   W7  (guard, E-5) "Use Existing" on a project whose scripts/<omic> folders were deleted recreates them, empty, before
#       each omic's module server starts (base: recreated holding the workbooks);
#   W8  "Push to GitHub" with no R Markdown file says so; with one report template it uses that file (base: several
#       workbooks, so utils::menu stops in a non-interactive R);
#   W9  ?integrationTrainModel holds a Limitations section with INT-28's content terms ("synthetic" or "simulated"
#       accepted alike), and ?integrationFindInputs names integrationTrainModel (base: absent).
# W6 (the installed library and the built tarball) is the script docs/qa/workflow-execution/kit/kit-234/w6.R in
# run-tests' scripts job. The application cases drive app_server with testServer as test-app-shell-shared.R does, with
# the real setupDirectories under a temporary root and the module servers and UI stood in.

# Copied from tests/testthat/test-app-shell-shared.R:5-42 (localNamespaceBinding, localNamespaceBindings), renamed.
k234LocalBinding <- function(env, name, value, .local_envir = parent.frame()) {
    had_binding <- exists(name, envir = env, inherits = FALSE)
    old_value <- if (had_binding) get(name, envir = env, inherits = FALSE) else NULL
    was_locked <- had_binding && bindingIsLocked(name, env)
    if (was_locked) unlockBinding(name, env)
    assign(name, value, envir = env)
    if (was_locked) lockBinding(name, env)
    withr::defer({
        if (exists(name, envir = env, inherits = FALSE) && bindingIsLocked(name, env)) unlockBinding(name, env)
        if (had_binding) {
            assign(name, old_value, envir = env)
        } else if (exists(name, envir = env, inherits = FALSE)) {
            rm(list = name, envir = env)
        }
        if (was_locked && exists(name, envir = env, inherits = FALSE)) lockBinding(name, env)
    }, envir = .local_envir)
}

k234LocalBindings <- function(env, bindings, .local_envir = parent.frame()) {
    for (name in names(bindings)) k234LocalBinding(env, name, bindings[[name]], .local_envir = .local_envir)
}

k234Rmd <- function(dir) sort(list.files(dir, pattern = "\\.[Rr][Mm][Dd]$", recursive = TRUE))

# The shell's stand-ins (test-app-shell-shared.R:172-269): logger, modals, notifications, the folder browser, the tab
# switch, and per omic a UI and a module server that records, at its call, whether scripts/<omic> exists and which
# R Markdown files it holds.
k234LocalShell <- function(project, omics, captured, .local_envir = parent.frame()) {
    stand_ins <- list(setup_shiny_logger = function() invisible(NULL), log_messages = function() "")
    for (omic in omics) {
        stand_ins[[paste0("mod_", omic, "_ui")]] <- local({
            name <- omic
            function(id) shiny::div(`data-id` = id, name)
        })
        stand_ins[[paste0("mod_", omic, "_server")]] <- local({
            name <- omic
            function(id, project_dirs, omic_type, experiment_label, volumes = NULL) {
                dir <- file.path(project, "scripts", name)
                captured$at_call[[name]] <- list(exists = dir.exists(dir), rmd = k234Rmd(dir))
                captured$project_dirs[[name]] <- project_dirs[[name]]
                list(kind = name)
            }
        })
    }
    k234LocalBindings(asNamespace("MultiScholaR"), stand_ins, .local_envir = .local_envir)
    testthat::local_mocked_bindings(showModal = function(...) invisible(NULL), removeModal = function(...) invisible(NULL),
        showNotification = function(...) invisible(NULL), .package = "shiny", .env = .local_envir)
    testthat::local_mocked_bindings(getVolumes = function() function() c(Home = dirname(project)),
        shinyDirChoose = function(...) invisible(NULL), parseDirPath = function(...) dirname(project),
        .package = "shinyFiles", .env = .local_envir)
    testthat::local_mocked_bindings(updateTabItems = function(...) invisible(NULL), .package = "shinydashboard",
        .env = .local_envir)
    testthat::local_mocked_bindings(log_info = function(...) invisible(NULL), log_warn = function(...) invisible(NULL),
        log_error = function(...) invisible(NULL), .package = "logger", .env = .local_envir)
}

# One session of the application on <root>/<label>: the omics chosen, the project set up (a new one, or "Use Existing"
# when the folder exists), the watchdog run. Returns what the stand-in module servers recorded.
k234Session <- function(root, label, omics, reuse = FALSE) {
    project <- file.path(root, label)
    captured <- new.env(parent = emptyenv())
    captured$at_call <- list()
    captured$project_dirs <- list()
    k234LocalShell(project, omics, captured)
    utils::capture.output(shiny::testServer(app_server, {
        session$setInputs(selected_omics = omics)
        session$setInputs(start_analysis = 1)
        session$flushReact()
        session$setInputs(experiment_label = label, project_dir = root)
        session$setInputs(confirm_setup = 1)
        session$flushReact()
        if (reuse) {
            session$setInputs(reuse_project = 1)
            session$flushReact()
        }
        session$flushReact()
        captured$initialized <- values$initialized_omics
    }))
    captured
}

# A project made directly by setupDirectories(force = TRUE), as an earlier session would have made it.
k234ExistingProject <- function(root, label, omics) {
    project <- file.path(root, label)
    dir.create(project, recursive = TRUE)
    testthat::local_mocked_bindings(log_info = function(...) invisible(NULL), log_warn = function(...) invisible(NULL),
        .package = "logger")
    utils::capture.output(setupDirectories(base_dir = project, omic_types = omics, force = TRUE))
    project
}

test_that("KIT-234 W1: the source tree holds no workbooks folder that system.file finds", {
    expect_identical(system.file("workbooks", package = "MultiScholaR"), "")
})

test_that("KIT-234 W2: a new proteomics and metabolomics project gets its scripts folders and no workbook", {
    root <- withr::local_tempdir("kit234-w2-")
    captured <- k234Session(root, "kit234_project", c("proteomics", "metabolomics"))
    expect_setequal(captured$initialized, c("proteomics", "metabolomics"))
    for (omic in c("proteomics", "metabolomics")) {
        dir <- file.path(root, "kit234_project", "scripts", omic)
        expect_true(dir.exists(dir), info = omic)
        expect_identical(k234Rmd(dir), character(), info = omic)
    }
})

test_that("KIT-234 W3: Use Existing leaves a user's edited workbook copy and own file byte for byte", {
    root <- withr::local_tempdir("kit234-w3-")
    project <- k234ExistingProject(root, "kit234_project", "proteomics")
    scripts <- file.path(project, "scripts", "proteomics")
    writeLines("user edit", file.path(scripts, "DIA_workflow_starter.rmd"))
    writeLines("mine <- TRUE", file.path(scripts, "mine.R"))
    before <- tools::md5sum(file.path(scripts, c("DIA_workflow_starter.rmd", "mine.R")))
    captured <- k234Session(root, "kit234_project", "proteomics", reuse = TRUE)
    expect_identical(captured$initialized, "proteomics")
    after <- tools::md5sum(file.path(scripts, c("DIA_workflow_starter.rmd", "mine.R")))
    expect_identical(unname(after), unname(before))
})

test_that("KIT-234 W4: with the integration tab on, scripts/integration is made and holds no workbook", {
    withr::local_options(MultiScholaR.integration_gui = TRUE)
    root <- withr::local_tempdir("kit234-w4-")
    captured <- k234Session(root, "kit234_project", c("proteomics", "integration"))
    expect_true(isTRUE(captured$at_call$integration$exists))
    dir <- file.path(root, "kit234_project", "scripts", "integration")
    expect_true(dir.exists(dir))
    expect_identical(k234Rmd(dir), character())
})

test_that("KIT-234 W5 (guard): the Integration tab's session store opens on a new project's integration folder", {
    testthat::skip_if_not(workflowOwnerCapable(), "the integration store needs qs2, DuckDB and Arrow")
    withr::local_options(MultiScholaR.integration_gui = TRUE)
    root <- withr::local_tempdir("kit234-w5-")
    captured <- k234Session(root, "kit234_project", c("proteomics", "integration"))
    owner <- NULL
    expect_no_error(owner <- integrationSessionOwner(captured$project_dirs$integration))
    if (!is.null(owner)) closeWorkflowStateOwner(owner)
})

test_that("KIT-234 W7 (guard, E-5): Use Existing recreates deleted scripts folders, empty, before the module servers", {
    root <- withr::local_tempdir("kit234-w7-")
    project <- k234ExistingProject(root, "kit234_project", c("proteomics", "metabolomics"))
    for (omic in c("proteomics", "metabolomics")) {
        unlink(file.path(project, "scripts", omic), recursive = TRUE)
        stopifnot(!dir.exists(file.path(project, "scripts", omic)), dir.exists(file.path(project, "results", omic)))
    }
    captured <- k234Session(root, "kit234_project", c("proteomics", "metabolomics"), reuse = TRUE)
    for (omic in c("proteomics", "metabolomics")) {
        expect_true(isTRUE(captured$at_call[[omic]]$exists), info = omic)
        expect_identical(captured$at_call[[omic]]$rmd, character(), info = omic)
    }
})

# Copied from tests/testthat/test-github-management-direct-shared.R:10-37 (localNamespaceBinding), renamed.
k234PackageBinding <- function(name, value, package, env = parent.frame()) {
    ns <- asNamespace(package)
    had_value <- exists(name, envir = ns, inherits = FALSE)
    old_value <- if (had_value) get(name, envir = ns, inherits = FALSE) else NULL
    was_locked <- had_value && bindingIsLocked(name, ns)
    if (was_locked) unlockBinding(name, ns)
    assign(name, value, envir = ns)
    if (was_locked) lockBinding(name, ns)
    withr::defer({
        if (exists(name, envir = ns, inherits = FALSE) && bindingIsLocked(name, ns)) unlockBinding(name, ns)
        if (had_value) {
            assign(name, old_value, envir = ns)
        } else if (exists(name, envir = ns, inherits = FALSE)) {
            rm(list = name, envir = ns)
        }
        if (was_locked) lockBinding(name, ns)
    }, envir = env)
}

test_that("KIT-234 W8: Push to GitHub refuses a project with no R Markdown file and uses its one report template", {
    testthat::skip_if_not_installed("git2r")
    testthat::skip_if_not_installed("gh")
    testthat::skip_if(interactive())
    root <- withr::local_tempdir("kit234-w8-")
    k234Session(root, "kit234_project", c("proteomics", "metabolomics"))
    project <- file.path(root, "kit234_project")
    scripts <- file.path(project, "scripts", "proteomics")
    # The token stood in as test-github-management-direct-shared.R:124 does; nothing reaches GitHub before the
    # .Rproj check that stops (b).
    k234PackageBinding("gh_token", function(...) "kit234-not-a-token", package = "gh")
    withr::local_options(github_org = "kit234-org")
    withr::defer(unlink(file.path(tempdir(), "kit234"), recursive = TRUE, force = TRUE))
    push <- function() {
        messages <- character()
        error <- tryCatch(withCallingHandlers(utils::capture.output(pushProjectToGithub(base_dir = project,
            source_dir = scripts, project_id = "kit234")), message = function(m) {
            messages <<- c(messages, conditionMessage(m))
            invokeRestart("muffleMessage")
        }), error = function(e) e)
        list(messages = messages, error = error)
    }
    # (a) No R Markdown file in the project.
    a <- push()
    expect_s3_class(a$error, "error")
    expect_match(conditionMessage(a$error), "No R Markdown files (.Rmd, .rmd, .RMD) found recursively in", fixed = TRUE)
    # (b) One report template, as Generate Report copies it.
    template <- system.file("reports", "proteomics", "DIANN_report.rmd", package = "MultiScholaR")
    stopifnot(nzchar(template), file.copy(template, scripts))
    b <- push()
    expect_true(any(grepl("Using the only R Markdown file found: scripts/proteomics/DIANN_report.rmd", b$messages,
        fixed = TRUE)), info = paste(b$messages, collapse = " | "))
    expect_s3_class(b$error, "error")
    expect_match(conditionMessage(b$error), "R project file does not exist", fixed = TRUE)
    # KIT-245 Part S (declared edit, ticket 6.1): the staged copy is removed on every exit, so the md5 of the staged
    # template is no longer there to read; the message above names the file chosen.
    expect_false(dir.exists(file.path(tempdir(), "kit234")))
})

test_that("KIT-234 W9: ?integrationTrainModel states INT-28's limitations and ?integrationFindInputs points to it", {
    root <- normalizePath(testthat::test_path("..", ".."), mustWork = FALSE)
    rd <- function(name) file.path(root, "man", paste0(name, ".Rd"))
    testthat::skip_if_not(file.exists(rd("integrationTrainModel")), "the package source tree is not present")
    flat <- function(path) gsub("[[:space:]]+", " ", paste(readLines(path, warn = FALSE), collapse = " "))
    train <- flat(rd("integrationTrainModel"))
    expect_true(grepl("\\section{Limitations}", train, fixed = TRUE))
    section <- if (grepl("\\section{Limitations}", train, fixed = TRUE)) sub("^.*\\\\section\\{Limitations\\}", "", train) else ""
    for (term in c("20 samples", "missing at random", "Gaussian", "STRING")) {
        expect_true(grepl(term, section, fixed = TRUE), info = term)
    }
    # The check against the reference script on generated data, in either word (the owner sets the wording).
    expect_true(grepl("synthetic|simulated", section), info = "synthetic or simulated")
    expect_true(grepl("integrationTrainModel", flat(rd("integrationFindInputs")), fixed = TRUE))
})
