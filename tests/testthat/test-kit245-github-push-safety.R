# KIT-245 Part S (ticket revision 1.2, section 10): pushProjectToGithub makes every local check before it creates a
# repository on GitHub, never writes the token to a file, removes its staged folder on every exit, validates the
# Project ID before any path is built, and names a repository it created but could not fill. Declared before the run,
# each case's base value beside it (ticket section 10):
#   S1  a push GitHub refuses (system answers 128): classed multischolar_github_push_failed naming <org>/<id>; at the
#       push the clone's .git/config does not hold the token and the credential comes from environment-borne git
#       configuration; afterwards no file under tempdir() holds the token and the environment is as before
#       (base: "Failed to push to repository ..."; the clone's .git/config holds the token);
#   S2  no git program: classed multischolar_github_git_missing, nothing created, nothing staged (base: POST, then the
#       push error, the staged folder left);
#   S3  a 101 MiB file in the source folder: classed multischolar_github_file_too_large naming it by its path in the
#       project (base: POST, then success);
#   S4  a Project ID that climbs out of tempdir() ("../<sentinel>"): classed multischolar_github_project_id_invalid,
#       the sentinel folder beside tempdir() untouched (base: the sentinel's file deleted by the clean-up of an
#       existing staged folder);
#   S5  git older than 2.31: classed multischolar_github_git_too_old, nothing created (base: POST, then success);
#   G1  no .Rproj: "R project file does not exist", nothing created (passes on the base too);
#   G2  two R Markdown files in a non-interactive R: the menu() error, nothing created (passes on the base too).
# Correction pass (AUDIT-KIT-245-impl-20261125, I01 to I03; base = 41171e21, the audited head):
#   S1  (revised, I01) at the push the credential helper is on git's command line, `-c credential.helper=` first and
#       the helper single-quoted (Unix) with the variable's name, never the token; MSR_GITHUB_PUSH_TOKEN holds the
#       token; no GIT_CONFIG_COUNT is set by the function (base: no -c options; GIT_CONFIG_COUNT "2");
#   S6  (I02, N01) a real `git credential fill` with no network, through the function's own git options, under a
#       global configuration holding two decoy helpers: the token answers for https://github.com, the decoy never
#       does; `git credential approve` leaves the decoy store's file unwritten; for https://example.com the token is
#       not given (base: the token given for example.com too; the github.com and store checks pass on Linux);
#   S7  (I03) refused Project IDs, with unlink and dir.create stood in so that no base can delete or create anything:
#       each gets classed multischolar_github_project_id_invalid with no unlink, dir.create or GitHub call; "-x" and
#       "a.b" are accepted (base: '...', 'a.', 'nul', 'NUL.txt' accepted and reach unlink or dir.create).
# No row reaches GitHub or runs a real git push: gh::gh is a recorder that never calls the original, gh::gh_token is
# stood in (gh 1.6.1 does not return a token-variable value that is not shaped like a token: KIT-245 G0), and
# system() is replaced in every row (S6's stand-in runs `git credential`, never `push`). The token is the plainly fake
# string below.

K245_DUMMY <- "kit245-dummy-not-a-token"
K245_ORG <- "kit245-test-org"

# As test-github-management-direct-shared.R:4-8 (makeFunctionWithOverrides), renamed.
k245WithOverrides <- function(fun, replacements) {
    environment(fun) <- list2env(replacements, parent = environment(fun))
    fun
}

# As test-github-management-direct-shared.R:10-37 (localNamespaceBinding), renamed.
k245Binding <- function(name, value, package, env = parent.frame()) {
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

# The GitHub API recorder (never calls the original) and the token stand-in.
k245Recorder <- function(env = parent.frame()) {
    rec <- new.env(parent = emptyenv())
    rec$calls <- character()
    k245Binding("gh", function(endpoint, ..., .token = NULL) {
        if (!identical(as.character(unclass(.token)), K245_DUMMY)) stop("kit245 recorder: a token other than the dummy")
        rec$calls <- c(rec$calls, endpoint)
        list(name = list(...)$name, private = TRUE)
    }, package = "gh", env = env)
    k245Binding("gh_token", function(...) K245_DUMMY, package = "gh", env = env)
    rec
}

# A project with one report template in scripts/<omic>/ and, unless asked otherwise, its <id>.Rproj.
k245Project <- function(id, omic = "proteomics", rproj = TRUE, second_template = FALSE, env = parent.frame()) {
    root <- withr::local_tempdir("kit245-project-", .local_envir = env)
    scripts <- file.path(root, "scripts", omic)
    dir.create(scripts, recursive = TRUE)
    template <- system.file("reports", "proteomics", "DIANN_report.rmd", package = "MultiScholaR")
    stopifnot(nzchar(template), file.copy(template, scripts))
    if (second_template) {
        stopifnot(file.copy(system.file("reports", "proteomics", "DIANN_limpa_report.rmd", package = "MultiScholaR"),
            scripts))
    }
    writeLines(c("Run\tgroup", "s1\tA", "s2\tB"), file.path(scripts, "design_matrix.tab"))
    if (rproj) writeLines("Version: 1.0", file.path(root, paste0(id, ".Rproj")))
    withr::defer(unlink(file.path(tempdir(), id), recursive = TRUE, force = TRUE, expand = FALSE), envir = env)
    list(base_dir = root, source_dir = scripts)
}

k245Push <- function(fun, fx, id) {
    tryCatch({
        utils::capture.output(result <- suppressMessages(fun(base_dir = fx$base_dir, source_dir = fx$source_dir,
            project_id = id, github_org = K245_ORG, github_user_email = "kit245@example.invalid",
            github_user_name = "kit245-test")))
        result
    }, error = function(e) e)
}

k245Text <- function(e) {
    if (!inherits(e, "condition")) return(paste("no condition; the call returned", format(e)))
    gsub("\\s+", " ", conditionMessage(e))
}

k245HoldsDummy <- function(path) {
    size <- file.size(path)
    !is.na(size) && size > 0 && length(grepRaw(K245_DUMMY, readBin(path, "raw", size), fixed = TRUE)) > 0L
}

# Every file under root (recursive, hidden files included) whose raw bytes hold the literal dummy.
k245FilesHoldingDummy <- function(root) {
    files <- list.files(root, recursive = TRUE, all.files = TRUE, full.names = TRUE, no.. = TRUE)
    files[vapply(files, k245HoldsDummy, logical(1))]
}

k245GitPath <- function(...) c(git = "/usr/bin/git")
k245PushOk <- function(command, intern = TRUE, ignore.stderr = FALSE) character()

test_that("KIT-245 S1: a push GitHub refuses leaves no token on disk and names the repository", {
    testthat::skip_if_not_installed("git2r")
    rec <- k245Recorder()
    id <- "k245-s1"
    fx <- k245Project(id)
    seen <- new.env(parent = emptyenv())
    push <- k245WithOverrides(pushProjectToGithub, list(
        Sys.which = k245GitPath,
        system = function(command, intern = TRUE, ignore.stderr = FALSE) {
            config <- file.path(tempdir(), id, ".git", "config")
            seen$config_holds_dummy <- file.exists(config) && k245HoldsDummy(config)
            seen$command <- command
            seen$count <- Sys.getenv("GIT_CONFIG_COUNT", unset = NA)
            seen$token_variable_is_dummy <- identical(Sys.getenv("MSR_GITHUB_PUSH_TOKEN", unset = NA), K245_DUMMY)
            structure(character(), status = 128L)
        }))
    variables <- c("GIT_CONFIG_COUNT", "GIT_CONFIG_KEY_0", "GIT_CONFIG_VALUE_0", "GIT_CONFIG_VALUE_1",
        "MSR_GITHUB_PUSH_TOKEN")
    before <- Sys.getenv(variables, unset = NA)

    err <- k245Push(push, fx, id)

    expect_s3_class(err, "multischolar_github_push_failed")
    expect_match(k245Text(err), paste0(K245_ORG, "/", id, " was created on GitHub"), fixed = TRUE)
    expect_false(dir.exists(file.path(tempdir(), id)))
    expect_false(isTRUE(seen$config_holds_dummy))
    # I01: the helper settings are on git's command line, the empty reset first; the token never is.
    expect_match(seen$command, "^git -c credential\\.helper= -c ")
    expect_match(seen$command, " push origin refs/heads/main$")
    if (.Platform$OS.type == "unix") {
        expect_match(seen$command, "-c 'credential.helper=!f() {", fixed = TRUE)
    }
    expect_match(seen$command, "$MSR_GITHUB_PUSH_TOKEN", fixed = TRUE)
    expect_false(grepl(K245_DUMMY, seen$command, fixed = TRUE))
    expect_true(isTRUE(seen$token_variable_is_dummy))
    expect_identical(seen$count, before[["GIT_CONFIG_COUNT"]])
    expect_identical(Sys.getenv(variables, unset = NA), before)
    expect_identical(k245FilesHoldingDummy(tempdir()), character())
    expect_identical(rec$calls, "POST /orgs/{org}/repos")
})

test_that("KIT-245 S2: without the git program nothing is created on GitHub and nothing is staged", {
    testthat::skip_if_not_installed("git2r")
    rec <- k245Recorder()
    id <- "k245-s2"
    fx <- k245Project(id)
    push <- k245WithOverrides(pushProjectToGithub, list(
        Sys.which = function(...) c(git = ""),
        system = function(command, intern = TRUE, ignore.stderr = FALSE) structure(character(), status = 127L)))

    err <- k245Push(push, fx, id)

    expect_s3_class(err, "multischolar_github_git_missing")
    expect_length(rec$calls, 0L)
    expect_false(dir.exists(file.path(tempdir(), id)))
})

test_that("KIT-245 S3: a file over GitHub's 100 MiB limit is refused before anything is created, by its project path", {
    testthat::skip_if_not_installed("git2r")
    rec <- k245Recorder()
    id <- "k245-s3"
    fx <- k245Project(id, omic = "lipidomics")
    big <- file.path(fx$source_dir, "kit245-big.bin")
    con <- file(big, "wb")
    seek(con, 101 * 1024^2 - 1, rw = "write")
    writeBin(as.raw(0L), con)
    close(con)
    stopifnot(file.size(big) == 101 * 1024^2)
    push <- k245WithOverrides(pushProjectToGithub, list(Sys.which = k245GitPath, system = k245PushOk))

    err <- k245Push(push, fx, id)

    expect_s3_class(err, "multischolar_github_file_too_large")
    expect_match(k245Text(err), "scripts/lipidomics/kit245-big.bin", fixed = TRUE)
    expect_false(grepl("scripts/proteomics/kit245-big.bin", k245Text(err), fixed = TRUE))
    expect_length(rec$calls, 0L)
    expect_false(dir.exists(file.path(tempdir(), id)))
})

test_that("KIT-245 S4: a Project ID that leaves tempdir() is refused and deletes nothing", {
    testthat::skip_if_not_installed("git2r")
    rec <- k245Recorder()
    sentinel_name <- paste0("kit245-sentinel-", basename(tempfile("")))
    sentinel <- file.path(dirname(tempdir()), sentinel_name)
    dir.create(sentinel)
    writeLines("kit245 sentinel", file.path(sentinel, "keep.txt"))
    withr::defer(unlink(sentinel, recursive = TRUE, force = TRUE, expand = FALSE))
    fx <- k245Project("k245-s4")
    push <- k245WithOverrides(pushProjectToGithub, list(Sys.which = k245GitPath, system = k245PushOk))

    err <- k245Push(push, fx, paste0("../", sentinel_name))

    expect_s3_class(err, "multischolar_github_project_id_invalid")
    expect_true(file.exists(file.path(sentinel, "keep.txt")))
    expect_length(rec$calls, 0L)
})

test_that("KIT-245 S5: git older than 2.31 is refused before anything is created", {
    testthat::skip_if_not_installed("git2r")
    rec <- k245Recorder()
    id <- "k245-s5"
    fx <- k245Project(id)
    push <- k245WithOverrides(pushProjectToGithub, list(
        Sys.which = k245GitPath,
        system2 = function(command, args = character(), ...) "git version 2.30.0",
        system = k245PushOk))

    err <- k245Push(push, fx, id)

    expect_s3_class(err, "multischolar_github_git_too_old")
    expect_match(k245Text(err), "2.30.0", fixed = TRUE)
    expect_length(rec$calls, 0L)
    expect_false(dir.exists(file.path(tempdir(), id)))
})

test_that("KIT-245 G1: without the project's .Rproj file nothing is created on GitHub", {
    testthat::skip_if_not_installed("git2r")
    rec <- k245Recorder()
    id <- "k245-g1"
    fx <- k245Project(id, rproj = FALSE)
    push <- k245WithOverrides(pushProjectToGithub, list(Sys.which = k245GitPath, system = k245PushOk))

    err <- k245Push(push, fx, id)

    expect_s3_class(err, "error")
    expect_match(conditionMessage(err), "R project file does not exist", fixed = TRUE)
    expect_length(rec$calls, 0L)
})

test_that("KIT-245 G2: two R Markdown files in a non-interactive R stop before anything is created on GitHub", {
    testthat::skip_if_not_installed("git2r")
    testthat::skip_if(interactive())
    rec <- k245Recorder()
    id <- "k245-g2"
    fx <- k245Project(id, second_template = TRUE)
    push <- k245WithOverrides(pushProjectToGithub, list(Sys.which = k245GitPath, system = k245PushOk))

    err <- k245Push(push, fx, id)

    expect_s3_class(err, "error")
    expect_match(conditionMessage(err), "menu() cannot be used non-interactively", fixed = TRUE)
    expect_length(rec$calls, 0L)
})

test_that("KIT-245 S6: the function's git options give the token to github.com only, past the user's own helpers", {
    testthat::skip_if_not_installed("git2r")
    testthat::skip_if(!nzchar(Sys.which("git")), "no git program")
    rec <- k245Recorder()
    id <- "k245-s6"
    fx <- k245Project(id)
    home <- withr::local_tempdir("kit245-home-")
    store <- file.path(home, "decoy-store")
    global <- file.path(home, "gitconfig")
    # Two decoy helpers in the user's global configuration: one that answers every request, one that stores.
    writeLines(c("[credential]",
        "\thelper = \"!f() { echo username=decoy; echo password=decoy; }; f\"",
        sprintf("\thelper = \"store --file %s\"", store)), global)
    withr::local_envvar(HOME = home, GIT_CONFIG_GLOBAL = global, GIT_CONFIG_NOSYSTEM = "1", GIT_TERMINAL_PROMPT = "0",
        GIT_ASKPASS = NA, SSH_ASKPASS = NA, GIT_CONFIG_COUNT = NA, GIT_CONFIG_PARAMETERS = NA)
    seen <- new.env(parent = emptyenv())
    # The function's own command, with `push origin refs/heads/main` replaced by a credential operation: no network.
    credential <- function(command, operation, input) {
        stopifnot(grepl(" push origin refs/heads/main$", command))
        suppressWarnings(base::system(sub("push origin refs/heads/main$", paste("credential", operation), command),
            intern = TRUE, input = input, ignore.stderr = TRUE))
    }
    push <- k245WithOverrides(pushProjectToGithub, list(
        system = function(command, intern = TRUE, ignore.stderr = FALSE) {
            seen$github <- credential(command, "fill", c("protocol=https", "host=github.com", ""))
            seen$other <- credential(command, "fill", c("protocol=https", "host=example.com", ""))
            credential(command, "approve", c("protocol=https", "host=github.com", "username=x-access-token",
                paste0("password=", K245_DUMMY), ""))
            structure(character(), status = 128L)
        }))

    k245Push(push, fx, id)

    expect_true(paste0("password=", K245_DUMMY) %in% seen$github)
    expect_false("password=decoy" %in% seen$github)
    expect_false(any(grepl(K245_DUMMY, seen$other, fixed = TRUE)))
    expect_false(file.exists(store))
    expect_identical(rec$calls, "POST /orgs/{org}/repos")
})

test_that("KIT-245 S7: Project IDs that name no plain folder are refused before anything is deleted or created", {
    testthat::skip_if_not_installed("git2r")
    rec <- k245Recorder()
    fx <- k245Project("k245-s7", rproj = FALSE)
    calls <- new.env(parent = emptyenv())
    calls$n <- 0L
    # unlink and dir.create stood in: on a base without the check nothing can be deleted or created.
    push <- k245WithOverrides(pushProjectToGithub, list(Sys.which = k245GitPath, system = k245PushOk,
        unlink = function(x, recursive = FALSE, force = FALSE, expand = TRUE) {
            calls$n <- calls$n + 1L
            0L
        },
        dir.create = function(path, showWarnings = TRUE, recursive = FALSE, mode = "0777") {
            calls$n <- calls$n + 1L
            TRUE
        }))
    refused <- list("..", ".", "...", "a.", "nul", "NUL.txt", "con", "COM1", "lpt9.log", "a/b", "a\\b", "/x", "C:x", "",
        NA_character_, c("a", "b"), strrep("a", 101), "é", "а", "a b", "a*")
    for (project_id in refused) {
        calls$n <- 0L
        label <- paste(deparse(project_id), collapse = "")
        err <- k245Push(push, fx, project_id)
        expect_true(inherits(err, "multischolar_github_project_id_invalid"), info = label)
        expect_identical(calls$n, 0L, info = label)
    }
    expect_length(rec$calls, 0L)
    # Accepted IDs stay inside tempdir(): the real unlink and dir.create, so the call reaches the .Rproj check
    # (run 45617558634: with the stand-ins it stopped earlier, writing .gitignore into a folder never created).
    accept <- k245WithOverrides(pushProjectToGithub, list(Sys.which = k245GitPath, system = k245PushOk))
    for (project_id in c("-x", "a.b")) {
        withr::defer(unlink(file.path(tempdir(), project_id), recursive = TRUE, force = TRUE, expand = FALSE))
        err <- k245Push(accept, fx, project_id)
        expect_false(inherits(err, "multischolar_github_project_id_invalid"), info = project_id)
        expect_match(k245Text(err), "R project file does not exist", fixed = TRUE, info = project_id)
    }
})
