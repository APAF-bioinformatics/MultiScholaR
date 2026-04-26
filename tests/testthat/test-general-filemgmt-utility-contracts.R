# fidelity-coverage-compare: shared
library(testthat)

makeFunctionWithOverrides <- function(fun, replacements) {
  fun_override <- fun
  environment(fun_override) <- list2env(replacements, parent = environment(fun))
  fun_override
}

test_that("testRequiredFiles logs and quits once per missing file", {
  fixture_dir <- tempfile("general-filemgmt-required-files-")
  dir.create(fixture_dir, recursive = TRUE)
  on.exit(unlink(fixture_dir, recursive = TRUE, force = TRUE), add = TRUE)

  existing_file <- file.path(fixture_dir, "present.txt")
  writeLines("ok", existing_file)
  missing_files <- file.path(fixture_dir, c("missing-a.txt", "missing-b.txt"))

  logged_messages <- character()
  quit_calls <- 0L
  helper_under_test <- makeFunctionWithOverrides(
    testRequiredFiles,
    list(
      logerror = function(fmt, file) {
        logged_messages <<- c(logged_messages, sprintf(fmt, file))
        invisible(NULL)
      },
      q = function(...) {
        quit_calls <<- quit_calls + 1L
        invisible(NULL)
      }
    )
  )

  result <- withVisible(helper_under_test(c(existing_file, missing_files)))

  expect_false(result$visible)
  expect_identical(
    logged_messages,
    sprintf("Missing required file: %s", missing_files)
  )
  expect_equal(quit_calls, 2)
})

test_that("testRequiredFilesWarning logs once per missing file without quitting", {
  fixture_dir <- tempfile("general-filemgmt-required-files-warning-")
  dir.create(fixture_dir, recursive = TRUE)
  on.exit(unlink(fixture_dir, recursive = TRUE, force = TRUE), add = TRUE)

  existing_file <- file.path(fixture_dir, "present.txt")
  writeLines("ok", existing_file)
  missing_files <- file.path(fixture_dir, c("missing-a.txt", "missing-b.txt"))

  logged_messages <- character()
  helper_under_test <- makeFunctionWithOverrides(
    testRequiredFilesWarning,
    list(
      logwarn = function(fmt, file) {
        logged_messages <<- c(logged_messages, sprintf(fmt, file))
        invisible(NULL)
      }
    )
  )

  result <- withVisible(helper_under_test(c(existing_file, missing_files)))

  expect_false(result$visible)
  expect_identical(
    logged_messages,
    sprintf("Missing required file: %s", missing_files)
  )
})

test_that("testRequiredArguments logs and quits once per missing argument", {
  logged_messages <- character()
  quit_calls <- 0L
  helper_under_test <- makeFunctionWithOverrides(
    testRequiredArguments,
    list(
      logerror = function(fmt, parameter) {
        logged_messages <<- c(logged_messages, sprintf(fmt, parameter))
        invisible(NULL)
      },
      q = function(...) {
        quit_calls <<- quit_calls + 1L
        invisible(NULL)
      }
    )
  )

  result <- withVisible(helper_under_test(
    arg_list = list(alpha = 1, gamma = 3),
    parameters = c("alpha", "beta", "gamma", "delta")
  ))

  expect_false(result$visible)
  expect_identical(
    logged_messages,
    c(
      "Missing required argument: beta",
      "Missing required argument: delta"
    )
  )
  expect_equal(quit_calls, 2)
})

test_that("parseType currently returns the original list without caller-visible coercion", {
  original <- list(alpha = "1.5", beta = "2", gamma = "keep")
  parsed <- parseType(
    arg_list = original,
    parameters = c("alpha", "beta"),
    functType = as.numeric
  )

  expect_identical(parsed, original)
})

test_that("parseString currently returns the original list without caller-visible quote stripping", {
  original <- list(
      double_quoted = "\"alpha\"",
      single_quoted = "'beta'",
      untouched = "gamma"
    )
  parsed <- parseString(
    arg_list = original,
    parameters = c("double_quoted", "single_quoted", "missing")
  )

  expect_identical(parsed, original)
})

test_that("parseList currently returns the original list without caller-visible splitting", {
  original <- list(formats = "pdf, png, svg", untouched = "keep")
  parsed <- parseList(
    arg_list = original,
    parameters = "formats"
  )

  expect_identical(parsed, original)
})

test_that("isArgumentDefined treats present null entries as defined but rejects missing and empty values", {
  arg_list <- list(
    present = "value",
    empty = "",
    null_value = NULL
  )

  expect_true(isArgumentDefined(arg_list, "present"))
  expect_false(isArgumentDefined(arg_list, "empty"))
  expect_true(isArgumentDefined(arg_list, "null_value"))
  expect_false(isArgumentDefined(arg_list, "missing"))
})
