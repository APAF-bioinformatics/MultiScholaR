library(testthat)

test_that("module CI artifact helpers create deterministic case directories", {
  old_root <- Sys.getenv("MULTISCHOLAR_MODULE_CI_ARTIFACT_DIR", unset = NA_character_)
  root <- tempfile("module-ci-artifacts-")
  Sys.setenv(MULTISCHOLAR_MODULE_CI_ARTIFACT_DIR = root)
  withr::defer({
    if (is.na(old_root)) {
      Sys.unsetenv("MULTISCHOLAR_MODULE_CI_ARTIFACT_DIR")
    } else {
      Sys.setenv(MULTISCHOLAR_MODULE_CI_ARTIFACT_DIR = old_root)
    }
  })

  path <- module_ci_case_artifact_dir("MCI 001/foundation")
  expect_true(dir.exists(path))
  expect_match(path, "MCI_001_foundation", fixed = TRUE)
})

test_that("module CI browser availability returns a scalar logical", {
  available <- module_ci_browser_dependencies_available()
  expect_type(available, "logical")
  expect_length(available, 1L)
})
