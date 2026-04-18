#!/usr/bin/env Rscript

resolve_refactor_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)))
  }

  frames <- sys.frames()
  for (idx in rev(seq_along(frames))) {
    ofile <- frames[[idx]]$ofile
    if (!is.null(ofile)) {
      return(dirname(normalizePath(ofile, mustWork = TRUE)))
    }
  }

  normalizePath("tools/refactor", mustWork = TRUE)
}

refactor_dir <- resolve_refactor_dir()
source(file.path(refactor_dir, "fidelity_inventory.R"))

parse_args <- function(args) {
  opts <- list(
    repo_root = ".",
    test_file = NULL,
    load_package = FALSE,
    pretty = FALSE
  )

  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]

    if (identical(arg, "--repo-root")) {
      i <- i + 1L
      opts$repo_root <- args[[i]]
    } else if (identical(arg, "--test-file")) {
      i <- i + 1L
      opts$test_file <- args[[i]]
    } else if (identical(arg, "--load-package")) {
      opts$load_package <- TRUE
    } else if (identical(arg, "--pretty")) {
      opts$pretty <- TRUE
    } else {
      abort("Unknown argument: ", arg)
    }

    i <- i + 1L
  }

  if (is.null(opts$test_file)) {
    abort("--test-file is required")
  }

  opts
}

expectation_counts <- function(results) {
  counts <- c(success = 0L, failure = 0L, skip = 0L, warning = 0L)
  for (item in results) {
    item_class <- class(item)
    if ("expectation_failure" %in% item_class) {
      counts[["failure"]] <- counts[["failure"]] + 1L
    } else if ("expectation_skip" %in% item_class) {
      counts[["skip"]] <- counts[["skip"]] + 1L
    } else if ("expectation_warning" %in% item_class) {
      counts[["warning"]] <- counts[["warning"]] + 1L
    } else if ("expectation_success" %in% item_class) {
      counts[["success"]] <- counts[["success"]] + 1L
    }
  }
  counts
}

test_status <- function(counts) {
  if (counts[["failure"]] > 0L) {
    return("failed")
  }
  if (counts[["skip"]] > 0L && counts[["success"]] == 0L) {
    return("skipped")
  }
  "passed"
}

serialize_test_results <- function(results) {
  lapply(results, function(item) {
    counts <- expectation_counts(item$results)
    list(
      test_name = item$test,
      status = test_status(counts),
      expectation_counts = as.list(counts),
      messages = unname(vapply(
        item$results,
        function(result) {
          result$message %||% ""
        },
        character(1)
      ))
    )
  })
}

run_contract_file <- function(repo_root, test_file, load_package) {
  if (!requireNamespace("testthat", quietly = TRUE)) {
    abort("testthat package is required for contract execution")
  }

  warnings <- character()
  messages <- character()
  output <- character()
  top_error <- NULL
  reporter_results <- NULL

  with_handlers <- function(expr) {
    withCallingHandlers(
      expr,
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      },
      message = function(m) {
        messages <<- c(messages, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    )
  }

  output <- tryCatch(
    with_handlers(capture.output({
      old_wd <- getwd()
      on.exit(setwd(old_wd), add = TRUE)
      setwd(repo_root)

      if (isTRUE(load_package)) {
        if (!requireNamespace("devtools", quietly = TRUE)) {
          abort("devtools package is required when --load-package is set")
        }
        devtools::load_all(repo_root, quiet = TRUE)
      }

      reporter_results <- testthat::test_file(test_file, reporter = testthat::ListReporter$new())
    }, type = "output")),
    error = function(e) {
      top_error <<- list(
        message = conditionMessage(e),
        class = unname(class(e))
      )
      character()
    }
  )

  if (!is.null(top_error)) {
    return(list(
      file_path = test_file,
      status = "runner_error",
      load_package = load_package,
      test_count = 0L,
      failure_count = 0L,
      skip_count = 0L,
      tests = list(),
      warnings = warnings,
      messages = messages,
      output = output,
      error = top_error
    ))
  }

  tests <- serialize_test_results(reporter_results)
  statuses <- unname(vapply(tests, function(item) item$status, character(1)))
  list(
    file_path = test_file,
    status = if (any(statuses == "failed")) "failed" else "passed",
    load_package = load_package,
    test_count = length(tests),
    failure_count = sum(statuses == "failed"),
    skip_count = sum(statuses == "skipped"),
    tests = tests,
    warnings = warnings,
    messages = messages,
    output = output,
    error = NULL
  )
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  repo_root <- normalizePath(opts$repo_root, mustWork = TRUE)
  test_file <- opts$test_file
  if (!grepl("^/", test_file)) {
    test_file <- file.path(repo_root, test_file)
  }
  test_file <- normalizePath(test_file, mustWork = TRUE)
  result <- run_contract_file(repo_root, test_file, opts$load_package)

  cat(
    jsonlite::toJSON(
      result,
      auto_unbox = TRUE,
      pretty = isTRUE(opts$pretty),
      null = "null",
      na = "null"
    )
  )
  cat("\n")
}

tryCatch(
  main(),
  error = function(e) {
    writeLines(conditionMessage(e), con = stderr())
    quit(status = 1)
  }
)
