# fidelity-coverage-compare: shared
library(testthat)

makeFunctionWithOverrides <- function(fun, replacements) {
  fun_override <- fun
  environment(fun_override) <- list2env(replacements, parent = environment(fun))
  fun_override
}

localNamespaceBinding <- function(name, value, package, env = parent.frame()) {
  ns <- asNamespace(package)
  had_value <- exists(name, envir = ns, inherits = FALSE)
  old_value <- if (had_value) get(name, envir = ns, inherits = FALSE) else NULL
  was_locked <- had_value && bindingIsLocked(name, ns)

  if (was_locked) {
    unlockBinding(name, ns)
  }
  assign(name, value, envir = ns)
  if (was_locked) {
    lockBinding(name, ns)
  }

  withr::defer({
    if (exists(name, envir = ns, inherits = FALSE) && bindingIsLocked(name, ns)) {
      unlockBinding(name, ns)
    }
    if (had_value) {
      assign(name, old_value, envir = ns)
    } else if (exists(name, envir = ns, inherits = FALSE)) {
      rm(list = name, envir = ns)
    }
    if (was_locked) {
      lockBinding(name, ns)
    }
  }, envir = env)
}

pushProjectToGithubFromDirs <- get(
  "pushProjectToGithubFromDirs",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
pushProjectToGithub <- get(
  "pushProjectToGithub",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)

test_that("github management helpers preserve directory extraction, validation, and local repo staging flow", {
  skip_if_not_installed("gh")
  skip_if_not_installed("git2r")
  skip_if_not_installed("fs")

  forwarded <- new.env(parent = emptyenv())
  push_from_dirs <- makeFunctionWithOverrides(
    pushProjectToGithubFromDirs,
    list(
      pushProjectToGithub = function(base_dir,
                                     source_dir,
                                     project_id,
                                     github_org,
                                     github_user_email,
                                     github_user_name,
                                     commit_message,
                                     private) {
        forwarded$args <- list(
          base_dir = base_dir,
          source_dir = source_dir,
          project_id = project_id,
          github_org = github_org,
          github_user_email = github_user_email,
          github_user_name = github_user_name,
          commit_message = commit_message,
          private = private
        )
        invisible(TRUE)
      }
    )
  )

  expect_error(
    push_from_dirs(
      project_dirs = list(),
      omic_type = "proteomics",
      experiment_label = "pilot",
      project_id = "proj"
    ),
    "Cannot find directory structure",
    fixed = TRUE
  )
  expect_error(
    push_from_dirs(
      project_dirs = list(proteomics_pilot = list()),
      omic_type = "proteomics",
      experiment_label = "pilot",
      project_id = "proj"
    ),
    "Missing required paths",
    fixed = TRUE
  )

  push_from_dirs(
    project_dirs = list(
      proteomics_pilot = list(
        base_dir = "/tmp/base",
        source_dir = "/tmp/source"
      )
    ),
    omic_type = "proteomics",
    experiment_label = "pilot",
    project_id = "proj",
    github_org = "test-org",
    github_user_email = "user@example.com",
    github_user_name = "tester",
    commit_message = "Initial sync",
    private = FALSE
  )
  expect_identical(forwarded$args$base_dir, "/tmp/base")
  expect_identical(forwarded$args$source_dir, "/tmp/source")
  expect_identical(forwarded$args$project_id, "proj")
  expect_false(forwarded$args$private)

  localNamespaceBinding("gh_token", function() "TOKEN123", package = "gh")
  github_calls <- new.env(parent = emptyenv())
  localNamespaceBinding(
    "gh",
    function(endpoint, ...) {
      github_calls$endpoint <- endpoint
      github_calls$args <- list(...)
      list(status = "ok")
    },
    package = "gh"
  )

  project_root <- tempfile("github-management-")
  base_dir <- file.path(project_root, "base")
  source_dir <- file.path(project_root, "source")
  dir.create(base_dir, recursive = TRUE)
  dir.create(source_dir, recursive = TRUE)
  withr::defer(unlink(project_root, recursive = TRUE, force = TRUE))

  project_id <- "unitproj"
  writeLines(c("---", "title: Workflow"), file.path(base_dir, "workflow.Rmd"))
  writeLines("Version: 1.0", file.path(base_dir, paste0(project_id, ".Rproj")))
  writeLines("helper <- TRUE", file.path(source_dir, "helper.R"))
  writeLines("ignore me", file.path(source_dir, "data_cln.tab"))

  captured <- new.env(parent = emptyenv())
  push_main <- makeFunctionWithOverrides(
    pushProjectToGithub,
    list(
      system = function(command, intern = TRUE, ignore.stderr = FALSE) {
        captured$git_cmd <- command
        character()
      },
      unlink = function(x, recursive = FALSE, force = FALSE) {
        captured$unlinked <- c(captured$unlinked, x)
        if (grepl(project_id, x, fixed = TRUE)) {
          captured$temp_dir <- x
        }
        0L
      }
    )
  )

  expect_error(
    push_main(
      base_dir = base_dir,
      source_dir = source_dir,
      project_id = project_id,
      github_org = "your_organization_here"
    ),
    "Please specify your GitHub organization or username",
    fixed = TRUE
  )

  empty_base_dir <- file.path(project_root, "empty-base")
  dir.create(empty_base_dir)
  expect_error(
    push_main(
      base_dir = empty_base_dir,
      source_dir = source_dir,
      project_id = project_id,
      github_org = "test-org",
      github_user_email = "user@example.com",
      github_user_name = "tester"
    ),
    "No R Markdown files",
    fixed = TRUE
  )

  no_rproj_base <- file.path(project_root, "no-rproj")
  dir.create(no_rproj_base)
  writeLines(c("---", "title: Workflow"), file.path(no_rproj_base, "workflow.Rmd"))
  expect_error(
    push_main(
      base_dir = no_rproj_base,
      source_dir = source_dir,
      project_id = project_id,
      github_org = "test-org",
      github_user_email = "user@example.com",
      github_user_name = "tester"
    ),
    "R project file does not exist",
    fixed = TRUE
  )

  result <- push_main(
    base_dir = base_dir,
    source_dir = source_dir,
    project_id = project_id,
    github_org = "test-org",
    github_user_email = "user@example.com",
    github_user_name = "tester",
    commit_message = "Initial sync",
    private = TRUE
  )

  expect_true(isTRUE(result))
  expect_identical(github_calls$endpoint, "POST /orgs/{org}/repos")
  expect_identical(github_calls$args$org, "test-org")
  expect_identical(github_calls$args$.token, "TOKEN123")
  expect_identical(github_calls$args$name, project_id)
  expect_match(captured$git_cmd, "refs/heads/main")

  staged_dir <- captured$temp_dir
  expect_true(dir.exists(staged_dir))
  expect_true(file.exists(file.path(staged_dir, ".gitignore")))
  expect_true(file.exists(file.path(staged_dir, paste0(project_id, ".Rproj"))))
  expect_true(file.exists(file.path(staged_dir, "scripts", "proteomics", paste0(project_id, ".Rmd"))))
  expect_true(file.exists(file.path(staged_dir, "scripts", "proteomics", "helper.R")))
  expect_false(file.exists(file.path(staged_dir, "scripts", "proteomics", "data_cln.tab")))
})
