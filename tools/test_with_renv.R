#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

repo_root <- normalizePath(".", mustWork = TRUE)
Sys.setenv(
  RENV_PATHS_ROOT = file.path(repo_root, ".renv-paths"),
  RENV_PATHS_CACHE = file.path(repo_root, ".renv-paths", "cache"),
  RENV_PATHS_LIBRARY = file.path(repo_root, "renv", "library"),
  RENV_CONFIG_SANDBOX_ENABLED = "FALSE"
)

source(file.path(repo_root, "renv", "activate.R"))

run_restore <- "--restore" %in% args
args <- setdiff(args, "--restore")
run_load_all <- !"--no-load-all" %in% args
args <- setdiff(args, "--no-load-all")

suppressPackageStartupMessages({
  library(renv)
})

if (run_restore) {
  renv::restore(prompt = FALSE)
}

suppressPackageStartupMessages({
  if (run_load_all) {
    library(pkgload)
  }
  library(testthat)
})

if (run_load_all) {
  pkgload::load_all(repo_root, quiet = TRUE, export_all = TRUE)
}

test_paths <- args

if (!length(test_paths)) {
  testthat::test_dir(file.path(repo_root, "tests", "testthat"), stop_on_failure = TRUE)
} else {
  for (path in test_paths) {
    full_path <- if (grepl("^/", path)) path else file.path(repo_root, path)
    testthat::test_file(full_path, stop_on_failure = TRUE)
  }
}
