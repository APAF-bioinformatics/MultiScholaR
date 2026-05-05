module_ci_artifact_root <- function() {
  root <- Sys.getenv("MULTISCHOLAR_MODULE_CI_ARTIFACT_DIR", unset = "")
  if (!nzchar(root)) {
    root <- testthat::test_path("_module_ci_artifacts")
  }
  root
}

module_ci_case_artifact_dir <- function(case_id, create = TRUE) {
  safe_id <- gsub("[^A-Za-z0-9_.-]+", "_", as.character(case_id))
  path <- file.path(module_ci_artifact_root(), safe_id)
  if (isTRUE(create) && !dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  path
}

module_ci_browser_dependencies_available <- function() {
  if (exists("e2e_browser_dependencies_available", mode = "function")) {
    return(e2e_browser_dependencies_available())
  }
  requireNamespace("shinytest2", quietly = TRUE) &&
    requireNamespace("chromote", quietly = TRUE) &&
    nzchar(tryCatch(chromote::find_chrome(), error = function(e) ""))
}

module_ci_skip_if_no_browser <- function() {
  testthat::skip_if_not(
    module_ci_browser_dependencies_available(),
    "module-CI browser dependencies are not available"
  )
}

module_ci_new_app_driver <- function(case_id, ...) {
  if (!exists("e2e_new_app_driver", mode = "function")) {
    stop("e2e_new_app_driver is not loaded; source helper-e2e-browser.R first", call. = FALSE)
  }
  e2e_new_app_driver(
    case_id = case_id,
    artifact_dir = module_ci_case_artifact_dir(case_id),
    ...
  )
}
