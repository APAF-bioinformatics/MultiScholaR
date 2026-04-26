# fidelity-coverage-compare: shared
library(testthat)

loadDependencies <- get(
  "loadDependencies",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)

makeFunctionWithOverrides <- function(fun, replacements) {
  fun_override <- fun
  environment(fun_override) <- list2env(replacements, parent = environment(fun))
  fun_override
}

localNamespaceBinding <- function(env, name, value, .local_envir = parent.frame()) {
  had_binding <- exists(name, envir = env, inherits = FALSE)
  old_value <- if (had_binding) get(name, envir = env, inherits = FALSE) else NULL
  was_locked <- had_binding && bindingIsLocked(name, env)

  if (was_locked) {
    unlockBinding(name, env)
  }
  assign(name, value, envir = env)
  if (was_locked) {
    lockBinding(name, env)
  }

  withr::defer({
    if (exists(name, envir = env, inherits = FALSE) && bindingIsLocked(name, env)) {
      unlockBinding(name, env)
    }
    if (had_binding) {
      assign(name, old_value, envir = env)
    } else if (exists(name, envir = env, inherits = FALSE)) {
      rm(list = name, envir = env)
    }
    if (was_locked && exists(name, envir = env, inherits = FALSE)) {
      lockBinding(name, env)
    }
  }, envir = .local_envir)
}

test_that("loadDependencies preserves already-installed loading flow", {
  calls <- new.env(parent = emptyenv())
  calls$library <- character()
  calls$p_load <- character()
  calls$message <- character()

  localNamespaceBinding(
    asNamespace("pacman"),
    "p_load",
    function(char, character.only = TRUE) {
      calls$p_load <- c(calls$p_load, char)
      invisible(TRUE)
    }
  )

  helper <- makeFunctionWithOverrides(
    loadDependencies,
    list(
      requireNamespace = function(pkg, quietly = TRUE) TRUE,
      library = function(package, ...) {
        calls$library <- c(calls$library, deparse(substitute(package)))
        invisible(TRUE)
      },
      message = function(...) {
        calls$message <- c(calls$message, paste(..., collapse = ""))
        invisible(NULL)
      }
    )
  )

  expect_invisible(helper(verbose = TRUE))
  expect_identical(calls$library, "pacman")
  expect_true("tidyverse" %in% calls$p_load)
  expect_true("limma" %in% calls$p_load)
  expect_true("RUVIIIC" %in% calls$p_load)
  expect_true(any(grepl("already installed, loading", calls$message, fixed = TRUE)))
  expect_true(any(grepl("All dependencies processed successfully", calls$message, fixed = TRUE)))
})

test_that("loadDependencies preserves installation branches across CRAN, Bioconductor, and GitHub", {
  calls <- new.env(parent = emptyenv())
  calls$library <- character()
  calls$cran <- character()
  calls$bioc <- character()
  calls$github <- character()
  calls$p_load <- character()
  calls$message <- character()

  localNamespaceBinding(
    asNamespace("utils"),
    "install.packages",
    function(pkgs, ...) {
      calls$cran <- c(calls$cran, pkgs)
      invisible(TRUE)
    }
  )
  localNamespaceBinding(
    asNamespace("BiocManager"),
    "install",
    function(pkgs, ...) {
      calls$bioc <- c(calls$bioc, pkgs)
      invisible(TRUE)
    }
  )
  localNamespaceBinding(
    asNamespace("devtools"),
    "install_github",
    function(repo, force = FALSE, ...) {
      calls$github <- c(calls$github, repo)
      invisible(TRUE)
    }
  )
  localNamespaceBinding(
    asNamespace("pacman"),
    "p_load",
    function(char, character.only = TRUE) {
      calls$p_load <- c(calls$p_load, char)
      invisible(TRUE)
    }
  )

  helper <- makeFunctionWithOverrides(
    loadDependencies,
    list(
      requireNamespace = function(pkg, quietly = TRUE) FALSE,
      library = function(package, ...) {
        calls$library <- c(calls$library, deparse(substitute(package)))
        invisible(TRUE)
      },
      message = function(...) {
        calls$message <- c(calls$message, paste(..., collapse = ""))
        invisible(NULL)
      }
    )
  )

  expect_invisible(helper(verbose = TRUE))
  expect_identical(calls$library, "pacman")
  expect_true("pacman" %in% calls$cran)
  expect_true("BiocManager" %in% calls$cran)
  expect_true("devtools" %in% calls$cran)
  expect_true("tidyverse" %in% calls$cran)
  expect_true("limma" %in% calls$bioc)
  expect_true("ComplexHeatmap" %in% calls$bioc)
  expect_true("cran/RUVIIIC" %in% calls$github)
  expect_true("APAF-bioinformatics/GlimmaV2" %in% calls$github)
  expect_true("tidyverse" %in% calls$p_load)
  expect_true("limma" %in% calls$p_load)
  expect_true("Glimma" %in% calls$p_load)
  expect_true(any(grepl("Installing pacman", calls$message, fixed = TRUE)))
  expect_true(any(grepl("Installing limma from Bioconductor", calls$message, fixed = TRUE)))
  expect_true(any(grepl("Installing Glimma from GitHub", calls$message, fixed = TRUE)))
})
